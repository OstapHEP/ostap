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
try:
    from builtins import range
except ImportError:
    from __builtin__ import range 
from   ostap.core.pyrouts import *
import ostap.io.zipshelve as     DBASE
import ostap.io.root_file
import ostap.trees.trees
import ostap.parallel.kisa 
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
from ostap.utils.utils import CleanUp
testdata   = CleanUp.tempfile ( suffix = '.root' , prefix ='test_tools_reweight_' )
tag_data   = 'DATA_histogram'
tag_mc     = 'MC_tree'
dbname     = CleanUp.tempfile ( suffix = '.db'   , prefix ='test_tools_reweight_' )

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
    
    for i in range ( 200000 ) : 
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
from ostap.tools.reweight     import Weight, makeWeights, WeightingPlot
from ostap.fitting.selectors  import SelectorWithVars, Variable 

## start iterations:
for iter in range ( 0 , maxIter ) :

    weighting = (
        ## variable          address in DB    
        Weight.Var( 'x' , address = 'x-reweight'  ) , 
        )
    
    weighter   = Weight( dbname , weighting )
    ## variables to be used in MC-dataset 
    variables  = [
        ## Variable ( 'x'  , 'x-variable' , 0  , 100 , lambda s : s.x ) , 
        Variable ( 'x'  , 'x-variable' , 0  , 100 ) , 
        Variable ( 'weight' , 'weight' , accessor =  weighter      )  
        ]
    
    #
    ## create new "weighted" mcdataset
    # 
    selector = SelectorWithVars (
        variables ,
        '0<x && x<100 '
        )
    
    mctree.pprocess ( selector , chunk_size = len ( mctree ) // 20 )
    ##mctree.process ( selector )
    mcds = selector.data             ## new reweighted dataset

    #
    ## update weights
    #    
    plots    = [
        WeightingPlot( 'x'   , 'weight' , 'x-reweight'  , hdata , hmc )  
        ]
    
    more = makeWeights ( mcds , plots , dbname , delta = 0.001 )

    ## make MC-histogram 
    mcds .project  ( hmc , 'x' , 'weight'  )
    
    logger.info    ( 'Compare DATA and MC for iteration #%d' % iter )
    #
    ## compare the basic properties: mean, rms, skewness and kurtosis
    # 
    hdata.cmp_prnt ( hmc , 'DATA' , 'MC' , 'DATA vs MC' )
    #
    ## calculate the distances
    #
    dist = hdata.cmp_dist ( hmc , density = True )
    logger.info ('DATA-MC "distance"      %s' % dist )
    #
    ## calculate the 'orthogonality'
    #  
    cost = hdata.cmp_cos  ( hmc , density = True )
    logger.info ('DATA-MC "orthogonality" %s' % cost )
    #
    ## try to fit it DATA with MC and vice versa 
    #
    fit1 = hdata.cmp_fit ( hmc   , density = True )
    if fit1 and 0 == fit1.Status() :
        logger.info ( 'Fit DATA with MC   Prob=%.3g[%%] ' % ( fit1.Prob() * 100 ) )
    fit2 = hmc  .cmp_fit ( hdata , density  = True )
    if fit2 and 0 == fit2.Status() :
        logger.info ( 'Fit MC   with DATA Prob=%.3g[%%] ' % ( fit2.Prob() * 100 ) )
    #
    ## make chi2-comparison between data and MC
    #
    c2ndf,prob = hdata.cmp_chi2 ( hmc   , density = True )
    logger.info ( 'DATA/MC: chi2/ndf (%.4g) and Prob %.5g%% ' % ( c2ndf , prob*100 ) )
    c2ndf,prob = hmc  .cmp_chi2 ( hdata , density = True )
    logger.info ( 'MC/DATA: chi2/ndf (%.4g) and Prob %.5g%% ' % ( c2ndf , prob*100 ) )
    
    mn,mx = hdata.cmp_minmax ( hmc   , diff = lambda a,b : a/b , density = True )
    logger.info ("DATA*/MC   ``min/max-distance''[%%] (%s)/(%s) at x=%.1f/%.1f" % (
        (100*mn[1]-100).toString ( '%+.1f+-%.1f' ) ,
        (100*mx[1]-100).toString ( '%+.1f+-%.1f' ) , mn[0]  , mx[0] ) )
    mn,mx = hmc  .cmp_minmax ( hdata , diff = lambda a,b : b/a , density = True )
    logger.info ("DATA/MC*   ``min/max-distance''[%%] (%s)/(%s) at x=%.1f/%.1f" % (
        (100*mn[1]-100).toString ( '%+.1f+-%.1f' ) ,
        (100*mx[1]-100).toString ( '%+.1f+-%.1f' ) , mn[0]  , mx[0] ) )                 
     
    ## final density on data 
    data_density = hdata.density()
    ## final density on mc 
    mc_density   = hmc.density()
    data_density.red  ()
    mc_density  .blue ()
    data_density.draw ('e1')
    mc_density  .draw ('e1 same')
    time.sleep ( 5 ) 
    
    if not more : 
        logger.info    ( 'No more iterations are needed #%d' % iter )
        break
    
    if iter + 1 != maxIter :
        mcds.clear() 
        del mcds , selector
    else :
        del selector 


data_density.draw ('e1')
mc_density  .draw ('e1 same')
time.sleep(60)
dbroot.Close()


# =============================================================================
# The END 
# =============================================================================
