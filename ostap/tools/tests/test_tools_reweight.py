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
from   ostap.utils.timing     import timing
from   ostap.core.core        import Ostap 
from   ostap.logger.colorized import attention, allright
from   ostap.plotting.canvas  import use_canvas
from   ostap.utils.root_utils import batch_env 
from   ostap.utils.cleanup    import CleanUp
from   ostap.histos.histos    import h1_axis
from   ostap.logger.symbols   import iteration 
import ostap.io.zipshelve     as     DBASE
import ostap.io.root_file
import ostap.trees.trees
import ostap.parallel.kisa 
import ostap.core.pyrouts 
import ROOT, random, math, os, time  
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
## set batch from environment 
batch_env ( logger )
# =============================================================================

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


N1 = 500000
N2 = 500000
N3 = 400000

for i in range( 0, N1 ) :
    v = 100 + random.expovariate ( -1.0/60 ) 
    while v <   0 : v +=100
    while v > 100 : v -=100        
    hdata.Fill(v)
for i in range( 0, N2 ) :
    v = random.gauss(50,10) 
    while v <   0 : v +=100
    while v > 100 : v -=100        
    hdata.Fill(v)

# ========================================================================================
with ROOT.TFile.Open( testdata ,'recreate') as mc_file:
    mc_file.cd() 

    mctree  = ROOT.TTree ( tag_mc , 'mc-tree' )
    mctree .SetDirectory ( mc_file ) 
    
    from array import array 
    xvar = array  ( 'f', [0])
    mctree.Branch ( 'x' , xvar , 'x/F' )
    
    for i in range ( N3 ) : 
        xvar[0] = random.expovariate(1.0/80)            
        mctree.Fill()
        
    mctree.Write()
    
    mc_file[tag_data] = hdata
    
## Read data from 
dbroot = ROOT.TFile.open ( testdata , 'r' ) 
logger.info ( 'Test data is fetched from "%s"' % testdata )   
dbroot.ls()
hdata  = dbroot[ tag_data ]

mctree = ROOT.TChain ( tag_mc ) ; mctree.Add ( testdata ) 

## prepare template histogram for MC 
hmc = ROOT.TH1D('hMC','histo-template for MC', 77 ,0,100 ) ; hmc.Sumw2() 

## prepare re-weighting machinery 
maxIter = 10  

## check database 
if not os.path.exists ( dbname ) :
    with DBASE.open   ( dbname , 'c' ) :
        logger.info("Create new weights DBASE '%s'" % dbname ) 
else :
    logger.info("Existing weights DBASE '%s' will be used" % dbname ) 

# =============================================================================
## make reweighting iterations
# =============================================================================
from   ostap.tools.reweight         import Weight, makeWeights, WeightingPlot, W2Data
from   ostap.fitting.pyselectors    import SelectorWithVars, Variable 
import ostap.parallel.parallel_fill

# =============================================================================
## weighting configuration: ## variable     address in DB    
weighting = ( Weight.Var     ( 'x'       ,  address = 'x-reweight'  ) , )
# ============================================================================
## variables to be used in MC-dataset 
variables = [ Variable ( 'x'  , 'x-variable' , 0  , 100 ) ]

plots      = [ WeightingPlot( 'x'   , 'weight' , 'x-reweight'  , hdata , hmc ) ]    

converged = False
active    = len ( plots )
# =============================================================================
## start iterations:
for iter in range ( 1 , maxIter + 1  ) :    
    
    tag = 'Reweighting iteration #%d%s' %  ( iter , iteration ) 
    logger.info ( allright ( tag ) ) 

    with timing ( tag + ': prepare MC-dataset:' , logger = logger ) : 

        # =============================================================================
        ## 0) The weighter object
        weighter = Weight ( dbname , weighting )
        
        # ===============================================================================
        ## 1a) create MC-dataset
        
        ## fill dataset from input CONTROL tree
        mcds , _ = mctree.make_dataset ( variables = variables      ,
                                         selection = '0<x && x<100' , 
                                         silent    = True           ) 
        
    with timing ( tag + ': add weight to MC-dataset' , logger = logger ) :
      
        ## 1b) add "weight" variable to the dataset
        mcds.add_reweighting ( weighter , name = 'weight' , progress = True , report = False )
        if 1 == iter % 10  : logger.info ( ( tag + ' MCDATA:\n%s' ) % mcds )
        
    with timing ( tag + ': make the actual reweighting:' , logger = logger ) :
        # ==============================================================================
        ## 2) update weights
        active , _ = makeWeights ( mcds               ,
                                   plots              ,
                                   dbname             ,
                                   delta      = 0.001 ,
                                   minmax     = 0.002 ,
                                   power      = 1.05  , ## tiny `overreweighting'
                                   make_plots = True  ,
                                   wtruncate  = ()    , 
                                   tag        = tag   )
        
    with timing ( tag + ': project weighted MC-dataset:' , logger = logger ) : 
        # ==============================================================================
        ## 3) make MC-histogram 
        hmc = mcds .project  ( hmc , 'x' , 'weight'  )

    with timing ( tag + ': compare DATA and MC distributions:' , logger = logger ) :  
        # ==============================================================================
        ## 4) compare "Data" and "MC"  after the reweighting on the given iteration    
        logger.info  ( tag + ': compare DATA and MC for iteration #%d%s' %  ( iter , iteration ) ) 

        hh = 'Iteration#%d: ' % iter 

        ## 4a) compare the basic properties: mean, rms, skewness and kurtosis 
        title = tag + ': DATA vs MC comparison'
        logger.info ( '%s:\n%s' % ( title , hdata.cmp_prnt
                                    ( hmc , density = True , title = title , prefix = '# ' ) ) )

        ## 4b) compare them        
        title = tag + ': DATA vs MC difference'
        logger.info ( '%s:\n%s' % ( title , hdata.cmp_diff_prnt
                                    ( hmc , density = True , title = title , prefix = '# ' ) ) )

    if not active and 5 < iter :
        
        logger.info    ( allright ( 'No more iterations, converged after #%d%s' % ( iter , iteration ) ) )
        title = 'Reweighted dataset after #%d iterations' % iter 
        logger.info ( '%s:\n%s' % ( title , mcds.table2 ( variables = [ 'x' ]  ,
                                                          title     = title    ,
                                                          cuts      = 'weight' , 
                                                          prefix    = '# '     ) ) )
        ## dump data as CVS file 
        cvs_file = CleanUp.tempfile ( suffix = '.csv' , prefix ='ostap-test-tools-reweight-' )
        mcds.to_csv ( cvs_file , dialect = 'excel-tab' )

        converged = True 
        break

    ## delete the dataset 
    mcds.clear() 

else :

    converged = False 
    logger.error ( "No convergency!" )

# ===========================================================================
title = 'Weighter object'
logger.info ( '%s:\n%s' % ( title , weighter.table ( prefix = '# ' ) ) )
# ============================================================================

# ============================================================================
## draw the convergency graphs 
graphs = weighter.graphs ()
for key in graphs : 
    with use_canvas ( "Convergency graph for '%s'" % key , wait = 4 ) :
        graph = graphs [ key ]
        graph.draw ( 'a' )
        
# =============================================================================
if converged :
    with timing ( "Add reweighting results to original MC-tree" , logger = logger ) :
        mctree.add_reweighting ( weighter , name = 'weight' )
    title = 'MC-tree with weights'
    logger.info ( '%s:\n%s' % ( title , mctree.table ( title = title , prefix = '# ') ) ) 
        
# ===========================================================================
with DBASE.open   ( dbname , 'r' ) as db :
    logger.info("(Reweighting database %s " % dbname ) 
    db.ls ()

# ===========================================================================
## convert to ROOT and back
from   ostap.tools.reweight import backup_to_ROOT, restore_from_ROOT
root_file = backup_to_ROOT    ( dbname     )
new_db    = restore_from_ROOT ( root_file  )

# ============================================================================

# =============================================================================
##                                                                      The END 
# =============================================================================
