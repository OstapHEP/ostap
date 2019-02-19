#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/tests/test_reweight2.py
#
#  Test for 2D-reweighting machinery
# 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-05-11
# =============================================================================
"""Test for 2D-reweighting machinery
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-05-10"
__all__     = ()  ## nothing to be imported 
# =============================================================================
import ROOT, random, math, os, time 
from   ostap.core.pyrouts import *
import ostap.io.zipshelve as     DBASE
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_reweight2' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================    
logger.info ( 'Test for 2D-Reweighting machinery')
# =============================================================================
import  tempfile
testdata   = os.path.join ( tempfile.gettempdir() , 'ostap_test_reweight2.root' )
tag_data   = 'DATA2_histogram'
tag_datax  = 'DATAX_histogram'
tag_datay  = 'DATAY_histogram'
tag_mc     = 'MC2_tree'
dbname     = os.path.join ( tempfile.gettempdir() , 'ostap_test_reweight2.db' )
 
if os.path.exists ( testdata ) : os.remove ( testdata ) 
if os.path.exists ( dbname   ) : os.remove ( dbname   )

import ostap.parallel.kisa

if not os.path.exists( testdata ) :
    #
    seed =  1234567890L 
    random.seed (  1234567890L ) 
    logger.info ( 'Test *RANDOM* data will be generated/seed=%s' % seed  )   
    ## prepare "data" histograms:
    # 1) 2D hstograms
    ix,iy  = 20 , 15
    hdata  = h2_axes ( [ 20.0/ix*i for i in range(ix+1)  ] ,
                       [ 15.0/iy*i for i in range(iy+1)  ] )
    # 2) non-equal binning 1D histogramm for x-component    
    hxdata = h1_axis ( [    i     for i in  range (5 ) ] +
                       [  5+i*0.2 for i in  range (50) ] +
                       [ 15+i     for i in  range (6 ) ] )
    # 2) equal binning 1D histogramm for y-component    
    hydata = h1_axis ( [ i*0.5    for i in  range (31) ]   )


    with ROOT.TFile.Open( testdata ,'recreate') as mc_file:
        mc_file.cd() 
        
        datatree  = ROOT.TTree ( 'DATA_tree', 'data-tree' )
        datatree .SetDirectory ( mc_file ) 
        
        from array import array 
        xvar = array  ( 'f', [0])
        yvar = array  ( 'f', [0])
        datatree.Branch ( 'x' , xvar , 'x/F' )
        datatree.Branch ( 'y' , yvar , 'y/F' )
    
        
        for i in range( 0, 3000000 ) :
            v1 = random.gauss ( 0 , 3 )
            v2 = random.gauss ( 0 , 2 )
            
            x  =  4 + v1 + v2   
            y  = 12      - v2
            
            while 20 < x : x-=20 
            while 15 < y : y-=15

            while 0 > x  : x+=20
            while 0 > y  : y+=15
            
            hdata .Fill(x,y)
            hxdata.Fill(x)
            hydata.Fill(y)
 
            xvar[0] = x 
            yvar[0] = y
            
            datatree.Fill()
            
           
        for i in range( 0, 3000000 ) :
            v1 = random.gauss ( 0 , 3 )
            v2 = random.gauss ( 0 , 2 )
            
            x  = 14 + v1 + v2   
            y  =  3      + v2
            
            while 20 < x : x-=20 
            while 15 < y : y-=15
            
            while 0 > x  : x+=20
            while 0 > y  : y+=15
            
            hdata .Fill(x,y)
            hxdata.Fill(x)
            hydata.Fill(y)
            
            xvar[0] = x 
            yvar[0] = y
            
            datatree.Fill()
            
        for i in range( 0, 4000000 ) :
            x = random.uniform ( 0 , 20 ) 
            y = random.uniform ( 0 , 15 )
            
            hdata .Fill(x,y)
            hxdata.Fill(x)
            hydata.Fill(y)
            
            xvar[0] = x 
            yvar[0] = y
            
            datatree.Fill()
            
        datatree.Write()
        
        ## write the histogram 
        mc_file[tag_data]  = hdata
        mc_file[tag_datax] = hxdata
        mc_file[tag_datay] = hydata
        
        mctree  = ROOT.TTree ( tag_mc , 'mc-tree' )
        mctree .SetDirectory ( mc_file ) 
        
        from array import array 
        xvar = array  ( 'f', [0])
        yvar = array  ( 'f', [0])
        mctree.Branch ( 'x' , xvar , 'x/F' )
        mctree.Branch ( 'y' , yvar , 'y/F' )
        
        for i in  xrange ( 500000 ) :

            xv = random.uniform ( 0 , 20 ) 
            yv = random.uniform ( 0 , 15 ) 

            while 20 < xv : xv-=20 
            while  0 > xv : xv+=20 

            xvar[0] = xv 
            yvar[0] = yv
            
            mctree.Fill()
            
        mctree.Write()
        mc_file.ls()
        
## Read data from DB
dbroot = ROOT.TFile.open ( testdata , 'r' )
logger.info ( 'Test data is fetched from DBASE "%s"' % testdata )   
dbroot.ls()
hdata    = dbroot[ tag_data    ]
hxdata   = dbroot[ tag_datax   ]
hydata   = dbroot[ tag_datay   ]
mctree   = dbroot[ tag_mc      ]
datatree = dbroot[ 'DATA_tree' ]
datastat = datatree.statCov('x','y')
#
## prebook random MC histograms
#
## ix = random.randint ( 35 , 40 ) 
## iy = random.randint ( 25 , 45 )
## ix,iy  = 16 , 16 
ix,iy  = 40 , 30
hmc  = h2_axes ( [ 20.0/ix*i for i in range(ix+1)  ] ,
                 [ 15.0/iy*i for i in range(iy+1)  ] )
## ix = random.randint ( 35 , 50 ) 
## iy = random.randint ( 30 , 50 )     
ix,iy  = 50 , 45 
hmcx = h1_axis ( [ 20.0/ix*i for i in  range(ix+1) ] )
hmcy = h1_axis ( [ 15.0/iy*i for i in  range(iy+1) ] )

## prepare re-weighting machinery 
maxIter = 50

## check database 
import os
if not os.path.exists( dbname ) :
    logger.info('Create new weights DBASE') 
    db = DBASE.open ( dbname , 'c' ) ##  create new empty db 
    db.close()
else :
    logger.info('Existing weights DBASE will be used') 
    
#
## make reweigthing iterations
# 
from ostap.tools.reweight      import Weight, makeWeights,  WeightingPlot  
from ostap.fitting.selectors   import SelectorWithVars, Variable 

## start iterations:
for iter in range ( 0 , maxIter ) :

    weightings = (
        ## variable          address in DB    
        #Weight.Var ( lambda s : s.x       , 'x-reweight'  ) , 
        #Weight.Var ( lambda s : s.y       , 'y-reweight'  ) , 
        #Weight.Var ( lambda s : (s.x,s.y) , '2D-reweight' ) , 
        Weight.Var (  'x'      , 'x-reweight'  ) , 
        Weight.Var (  'y'      , 'y-reweight'  ) , 
        Weight.Var ( ('x','y') , '2D-reweight' ) , 
        )
    
    weighter   = Weight( dbname , weightings )
    ## variables to be used in MC-dataset 
    variables  = [
        Variable ( 'x'      , 'x-var'  , 0  , 20 ) , 
        Variable ( 'y'      , 'y-var'  , 0  , 15 ) ,
        Variable ( 'weight' , 'weight' , accessor = weighter )  
        ]
    
    #
    ## create new "weighted" mcdataset
    # 
    selector = SelectorWithVars (
        variables ,
        '0<x && x<20 && 0<y && y<20'
        )

    mctree.pprocess ( selector , chunk_size = len(mctree) // 20 )
    
    mcds = selector.data             ## new reweighted dataset

    logger.info ('MCDATA: %s' %  mcds )
    
    #
    ## update weights
    #

    plots = [ WeightingPlot ( 'y:x' , 'weight' , '2D-reweight' , hdata  , hmc  ) ]
    if 3 < iter: 
        plots  = [
            WeightingPlot ( 'x'   , 'weight' , 'x-reweight'  , hxdata , hmcx       ) ,  
            WeightingPlot ( 'y'   , 'weight' , 'y-reweight'  , hydata , hmcy       ) , 
            WeightingPlot ( 'y:x' , 'weight' , '2D-reweight' , hdata  , hmc  , 0.5 ) , 
            ]
    
    ## more iteration?  number of ``active'' reweightings    
    more = makeWeights ( mcds , plots , dbname , delta = 0.015 , power = 2 if 1 != len(plots) else 1 ) 
    
    ## make MC-histogram 
    mcds .project  ( hmcx , 'x'   , 'weight'  )
    mcds .project  ( hmcy , 'y'   , 'weight'  )
    mcds .project  ( hmc  , 'y:x' , 'weight'  )
    
    logger.info    ( 'Compare DATA and MC for iteration #%d' % iter )
    #
    ## compare the basic properties: mean, rms, skewness and kurtosis
    # 
    hxdata.cmp_prnt ( hmcx , 'DATA' , 'MC' , 'DATA(x) vs MC(x)' )
    hydata.cmp_prnt ( hmcy , 'DATA' , 'MC' , 'DATA(y) vs MC(y)' )
    
    #
    ## calculate the distances
    #
    dist = hxdata.cmp_dist ( hmcx , density = True )
    logger.info ("DATA(x)-MC(x)  ``distance''        %s" % dist )
    dist = hydata.cmp_dist ( hmcy , density = True )
    logger.info ("DATA(y)-MC(y)  ``distance''        %s" % dist )

    #
    ## calculate the 'orthogonality'
    #  
    cost = hxdata.cmp_cos  ( hmcx , density = True )
    logger.info ("DATA(x)-MC(x)  ``orthogonality'' %s" % cost )
    cost = hydata.cmp_cos  ( hmcy , density = True )
    logger.info ("DATA(y)-MC(y)  ``orthogonality'' %s" % cost )
    #
    
    mn,mx = hxdata.cmp_minmax ( hmcx   , diff = lambda a,b : a/b , density = True )
    logger.info ("DATA*(x)/MC(x)   ``min/max-distance''[%%] (%s)/(%s) at x=%.1f/%.1f" % (
        (100*mn[1]-100).toString ( '%+.1f+-%.1f' ) ,
        (100*mx[1]-100).toString ( '%+.1f+-%.1f' ) , mn[0]  , mx[0] ) )
    mn,mx = hmcx  .cmp_minmax ( hxdata , diff = lambda a,b : b/a , density = True )
    logger.info ("DATA(x)/MC*(x)   ``min/max-distance''[%%] (%s)/(%s) at x=%.1f/%.1f" % (
        (100*mn[1]-100).toString ( '%+.1f+-%.1f' ) ,
        (100*mx[1]-100).toString ( '%+.1f+-%.1f' ) , mn[0]  , mx[0] ) )                 
    mn,mx = hydata.cmp_minmax ( hmcy   , diff = lambda a,b : a/b , density = True )
    logger.info ("DATA*(y)/MC(y)   ``min/max-distance''[%%] (%s)/(%s) at y=%.1f/%.1f" % (
        (100*mn[1]-100).toString ( '%+.1f+-%.1f' ) ,
        (100*mx[1]-100).toString ( '%+.1f+-%.1f' ) , mn[0]  , mx[0] ) ) 
    mn,mx = hmcy  .cmp_minmax ( hydata , diff = lambda a,b : b/a , density = True )
    logger.info ("DATA(y)/MC*(y)   ``min/max-distance''[%%] (%s)/(%s) at y=%.1f/%.1f" % (
        (100*mn[1]-100).toString ( '%+.1f+-%.1f' ) ,
        (100*mx[1]-100).toString ( '%+.1f+-%.1f' ) , mn[0]  , mx[0] ) )            
    mn,mx = hdata .cmp_minmax ( hmc   , diff = lambda a,b : a/b  , density = True )
    logger.info ("DATA*(xy)/MC(xy) ``min/max-distance''[%%] (%s)/(%s) at (x,y)=(%.1f,%.1f)/(%.1f,%.1f)" % (
        (100*mn[2]-100).toString ( '%+.1f+-%.1f' ) ,
        (100*mx[2]-100).toString ( '%+.1f+-%.1f' ) , mn[0]  , mn[1] , mx[0]  , mx[1] ) )            
    mn,mx = hmc   .cmp_minmax ( hdata , diff = lambda a,b : b/a  , density = True )
    logger.info ("DATA(xy)/MC*(xy) ``min/max-distance''[%%] (%s)/(%s) at (x,y)=(%.1f,%.1f)/(%.1f,%.1f)" % (
        (100*mn[2]-100).toString ( '%+.1f+-%.1f' ) ,
        (100*mx[2]-100).toString ( '%+.1f+-%.1f' ) , mn[0]  , mn[1] , mx[0]  , mx[1] ) )            


    mcstat = mcds.statCov('x','y','weight')
    logger.info ('MCSTAT:\nx=%s\ny=%s\ncov2:\n%s'   % mcstat  [:3] ) 
    logger.info ('DATASTAT:\nx=%s\ny=%s\ncov2:\n%s' %datastat[:3] ) 
    
    ## final density on data 
    datax_density = hxdata.density()
    ## final density on mc 
    mcx_density   = hmcx.density()

    ## final density on data 
    datay_density = hydata.density()
    ## final density on mc 
    mcy_density   = hmcy.density()
    
    datax_density.red   ()
    mcx_density  .blue  ()
    
    datay_density.green ()
    mcy_density  .yellow()

    datax_density.SetMinimum(0.00)
    datax_density.SetMaximum(0.12)
    datax_density.draw ('e1')
    mcx_density  .draw ('e1 same')
    datay_density.draw ('e1 same')
    mcy_density  .draw ('e1 same')
    time.sleep(5)

    if not more and iter > 6 : 
        logger.info    ( 'No more iterations are needed #%d' % iter )
        break
    
    if iter + 1 != maxIter :
        mcds.clear() 
        del mcds , selector
    else :
        del selector 


logger.info ('MCSTAT:\nx=%s\ny=%s\ncov2:\n%s'   %mcstat  [:3] ) 
logger.info ('DATASTAT:\nx=%s\ny=%s\ncov2:\n%s' %datastat[:3] ) 

datax_density.draw ('e1')
mcx_density  .draw ('e1 same')
datay_density.draw ('e1 same')
mcy_density  .draw ('e1 same')
time.sleep(60)
dbroot.Close() 

# =============================================================================
# The END 
# =============================================================================
