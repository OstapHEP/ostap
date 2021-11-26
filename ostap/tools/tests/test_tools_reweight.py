#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/tests/test_tools_reweight.py
#  Test for reweighting machinery
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-05-11
# =============================================================================
"""Test for reweighting machinery in  Ostap
"""
# =============================================================================
import ROOT, random, math, os, time  
from   builtins               import range
from   ostap.core.pyrouts     import *
import ostap.io.zipshelve     as     DBASE
import ostap.io.root_file
import ostap.trees.trees
import ostap.parallel.kisa 
from   ostap.utils.timing     import timing
from   ostap.logger.colorized import attention, allright  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_tools_reweight' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================    
logger.info ( 'Test for Reweighting machinery')
# =============================================================================
from ostap.utils.cleanup import CleanUp
testdata   = CleanUp.tempfile ( suffix = '.root' , prefix ='ostap-test-tools-reweight-' )
tag_data   = 'DATA_histogram'
tag_mc     = 'MC_tree'
dbname     = CleanUp.tempfile ( suffix = '.db'   , prefix ='ostap-test-tools-reweight-' )

if os.path.exists ( testdata ) : os.remove ( testdata ) 
if os.path.exists ( dbname   ) : os.remove ( dbname   ) 

## prepare data for tests 
seed = 1234567890
logger.info ( 'Test RANDOM data will be generated/seed=%s' % seed  )
random.seed ( seed )

## prepare"data" histogram: lets use simple exponential and non-equidistance bins 
hdata = h1_axis( [    i   for i in range(0,20) ] +
                 [ 20+i*2 for i in range(0,20) ] +
                 [ 60+i*4 for i in range(0,10) ] + [ 100 ] )

for i in range( 0, 500000 ) :
    v = 100 + random.expovariate ( -1.0/60 ) 
    while v <   0 : v +=100
    while v > 100 : v -=100        
    hdata.Fill(v)
for i in range( 0, 500000 ) :
    v = random.gauss(50,10) 
    while v <   0 : v +=100
    while v > 100 : v -=100        
    hdata.Fill(v)
    
with ROOT.TFile.Open( testdata ,'recreate') as mc_file:
    mc_file.cd() 

    mctree  = ROOT.TTree ( tag_mc , 'mc-tree' )
    mctree .SetDirectory ( mc_file ) 
    
    from array import array 
    xvar = array  ( 'f', [0])
    mctree.Branch ( 'x' , xvar , 'x/F' )
    
    for i in range ( 400000 ) : 
        xvar[0] = random.expovariate(1.0/80)            
        mctree.Fill()
        
    mctree.Write()
    
    mc_file[tag_data] = hdata
    
## Read data from 
dbroot = ROOT.TFile.open ( testdata , 'r' ) 
logger.info ( 'Test data is fetched from "%s"' % testdata )   
dbroot.ls()
hdata  = dbroot[ tag_data ]
mctree = dbroot[ tag_mc   ]

## prepare template histogram for MC 
hmc = ROOT.TH1D('hMC','histo-template for MC', 77 ,0,100 ) ; hmc.Sumw2() 

## prepare re-weighting machinery 
maxIter = 10  

## check database 
import os
if not os.path.exists( dbname ) :
    logger.info('Create new weights DBASE') 
    db = DBASE.open ( dbname , 'c' ) ##  create new empty db 
    db.close()
else :
    logger.info('Existing weights DBASE will be used') 

#
## make reweighting iterations
# 
from   ostap.tools.reweight         import Weight, makeWeights, WeightingPlot, W2Data
from   ostap.fitting.pyselectors    import SelectorWithVars, Variable 
import ostap.parallel.parallel_fill

# =============================================================================
## weighting configuration: ## variable     address in DB    
weighting = ( Weight.Var     ( 'x'       ,  address = 'x-reweight'  ) , )
# ============================================================================
## variables to be used in MC-dataset 
variables  = [ Variable ( 'x'  , 'x-variable' , 0  , 100 ) ]

# =============================================================================
## start iterations:
for iter in range ( 1 , maxIter + 1  ) :    
    
    tag = 'Reweighting iteration #%d' % iter
    logger.info ( allright ( tag ) ) 

    with timing ( tag + ': prepare MC-dataset:' , logger = logger ) : 

        # =============================================================================
        ## 0) The weighter object
        weighter = Weight ( dbname , weighting )
        
        # ===============================================================================
        ## 1a) create mcdataset
        selector = SelectorWithVars ( variables , '0<x && x<100 ' , silence = True )
        mctree.process ( selector , silent = True )
        mcds     = selector.data ## dataset
        
    with timing ( tag + ': add weight to MC-dataset' , logger = logger ) :
      
        ## 1b) add "weight" variable to the dataset
        mcds.add_reweighting ( weighter , name = 'weight' )
        if 1 == iter % 10  : logger.info ( ( tag + ' MCDATA:\n%s' ) %  mcds )
        
    with timing ( tag + ': make actual reweighting:' , logger = logger ) :
        # ==============================================================================
        ## 2) update weights
        plots   = [ WeightingPlot( 'x'   , 'weight' , 'x-reweight'  , hdata , hmc ) ]    
        more    = makeWeights ( mcds               ,
                                plots              ,
                                dbname             ,
                                delta      = 0.005 ,
                                minmax     = 0.01  ,
                                power      = 1.05  , ## tiny ``overreweighting''
                                make_plots = False , 
                                tag        = tag   )
        
    with timing ( tag + ': project weighted MC-dataset:' , logger = logger ) : 
       # ==============================================================================
       ## 3) make MC-histogram 
       mcds .project  ( hmc , 'x' , 'weight'  )

    with timing ( tag + ': compare DATA and MC distributions:' , logger = logger ) :  
        # ==============================================================================
        ## 4) compare "Data" and "MC"  after the reweighting on the given iteration    
        logger.info  ( tag + ': compare DATA and MC for iteration #%d' % iter )

        hh = 'Iteration#%d: ' % iter 
        
        ## 4a) compare the basic properties: mean, rms, skewness and kurtosis 
        title = tag + ': DATA vs MC comparison'
        logger.info ( '%s:\n%s' % ( title , hdata.cmp_prnt
                                    ( hmc , density = True , title = title , prefix = '# ' ) ) )
        
        ## 4b) compare them        
        title = tag + ': DATA vs MC difference'
        logger.info ( '%s:\n%s' % ( title , hdata.cmp_diff_prnt
                                    ( hmc , density = True , title = title , prefix = '# ' ) ) )
        
    # =========================================================================
    ## prepare the plot of weighted MC for the given iteration
    
    ## final density on data 
    data_density = hdata.density()
    
    ## final density on mc 
    mc_density   = hmc.density()
    
    data_density.red  ()
    mc_density  .blue ()
    data_density.draw ('e1')
    mc_density  .draw ('e1 same')
    time.sleep ( 5 ) 
    
    if not more and iter > 3 : 
        logger.info    ( allright ( 'No more iterations, converged after #%d' % iter ) )
        break

    cvs_file = CleanUp.tempfile ( suffix = '.csv' , prefix ='ostap-test-tools-reweight-' )
    mcds.to_csv ( cvs_file , dialect = 'excel-tab' )

    mcds.clear () 
    del mcds , selector

    
else :

    logger.error ( "No convergency!" )



data_density.draw ('e1')
mc_density  .draw ('e1 same')
time.sleep(10)

# =============================================================================
##                                                                      The END 
# =============================================================================
