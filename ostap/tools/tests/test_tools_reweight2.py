#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/tests/test_tools_reweight2.py
#  Test for 2D-reweighting machinery
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
from   builtins               import range
from   ostap.core.pyrouts     import *
import ostap.io.zipshelve     as     DBASE
from   ostap.utils.timing     import timing
from   ostap.logger.colorized import attention, allright  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_tools_reweight2' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================    
logger.info ( 'Test for 2D-Reweighting machinery')
# ============================================================================
from ostap.utils.cleanup import CleanUp
testdata   = CleanUp.tempfile ( suffix = '.root' , prefix ='test_tools_reweight2_' )
tag_data   = 'DATA2_histogram'
tag_datax  = 'DATAX_histogram'
tag_datay  = 'DATAY_histogram'
tag_mc     = 'MC2_tree'
dbname     = CleanUp.tempfile ( suffix = '.db' , prefix ='test_tools_reweight2_'   )
 
if os.path.exists ( testdata ) : os.remove ( testdata ) 
if os.path.exists ( dbname   ) : os.remove ( dbname   )

import ostap.parallel.kisa


if not os.path.exists( testdata ) :
    #
    seed =  1234567890 
    random.seed ( seed ) 
    logger.info ( 'Test *RANDOM* data will be generated/seed=%s' % seed  )   
    ## prepare "data" histograms:
    # 1) 2D histograms
    ix , iy  = 45 , 25
    hdata  = h2_axes ( [ 20.0/ix*i for i in range ( ix + 1 ) ] ,
                       [ 15.0/iy*i for i in range ( iy + 1 ) ] )
    # 2) non-equal binning 1D histograms for x-component    
    hxdata = h1_axis ( [    i      for i in range ( 5  ) ] +
                       [  5+i*0.2  for i in range ( 50 ) ] +
                       [ 15+i      for i in range ( 6  ) ] )
    # 2) equal binning 1D histograms for y-component    
    hydata = h1_axis ( [ i*0.5     for i in range ( 31 ) ] )


    with ROOT.TFile.Open ( testdata ,'recreate') as mc_file:
        mc_file.cd() 
        
        datatree  = ROOT.TTree ( 'DATA_tree', 'data-tree' )
        datatree .SetDirectory ( mc_file ) 
        
        from array import array 
        xvar = array  ( 'f', [0])
        yvar = array  ( 'f', [0])
        datatree.Branch ( 'x' , xvar , 'x/F' )
        datatree.Branch ( 'y' , yvar , 'y/F' )
    
        N1 = 1000000
        for i in range ( 0, N1 ) :
            
            v1 = random.gauss ( 0 , 3 )
            v2 = random.gauss ( 0 , 2 )
            
            x  =  4 + v1 + v2   
            y  = 12      - v2
            
            while 20 < x : x -= 20 
            while 15 < y : y -= 15

            while  0 > x : x += 20
            while  0 > y : y += 15
            
            hdata .Fill ( x , y )
            hxdata.Fill ( x )
            hydata.Fill ( y )
 
            xvar [ 0 ] = x 
            yvar [ 0 ] = y
            
            datatree.Fill()
            
           
        for i in range ( 0 , N1 ) :
            
            v1 = random.gauss ( 0 , 3 )
            v2 = random.gauss ( 0 , 2 )
            
            x  = 14 + v1 + v2   
            y  =  3      + v2
            
            while 20 < x : x -= 20 
            while 15 < y : y -= 15
            
            while  0 > x : x += 20
            while  0 > y : y += 15
            
            hdata .Fill ( x , y )
            hxdata.Fill ( x )
            hydata.Fill ( y )
            
            xvar [ 0 ] = x 
            yvar [ 0 ] = y
            
            datatree.Fill()
            
        for i in range ( 0 , 2 * N1 ) :
            
            x = random.uniform ( 0 , 20 ) 
            y = random.uniform ( 0 , 15 )
            
            hdata .Fill ( x , y )
            hxdata.Fill ( x )
            hydata.Fill ( y )
            
            xvar [ 0 ] = x 
            yvar [ 0 ] = y
            
            datatree.Fill()
            
        datatree.Write()
        
        ## write the histogram 
        mc_file [ tag_data  ] = hdata
        mc_file [ tag_datax ] = hxdata
        mc_file [ tag_datay ] = hydata
        
        mctree  = ROOT.TTree ( tag_mc , 'mc-tree' )
        mctree .SetDirectory ( mc_file ) 
        
        from array import array 
        xvar = array  ( 'f', [ 0.0 ] )
        yvar = array  ( 'f', [ 0.0 ] )
        mctree.Branch ( 'x' , xvar , 'x/F' )
        mctree.Branch ( 'y' , yvar , 'y/F' )

        N2 = 200000
        for i in  range ( N2 ) :

            xv = random.uniform ( 0 , 20 ) 
            yv = random.uniform ( 0 , 15 ) 

            xvar [ 0 ] = xv 
            yvar [ 0 ] = yv
            
            mctree.Fill()
            
        fx = lambda x : ( x/10.0-1 ) ** 2
        fy = lambda y : ( y/ 7.5-1 ) ** 2
        
        for i in  range ( N2 ) :
            
            while True : 
                xv = random.uniform ( 0 , 20 )
                if random.uniform   ( 0 ,  1 ) < fx ( xv ) : break
                
            while True : 
                yv = random.uniform ( 0 , 15 )
                if random.uniform   ( 0 ,  1 ) < fy ( yv ) : break 
                
            xvar [ 0 ] = xv 
            yvar [ 0 ] = yv
            
            mctree.Fill()
            
        mctree.Write()
        mc_file.ls()

# =============================================================================
## Read data from DB
# =============================================================================
dbroot = ROOT.TFile.open ( testdata , 'r' )
logger.info ( 'Test data is fetched from DBASE "%s"' % testdata )   
dbroot.ls()
hdata    = dbroot [ tag_data    ]
hxdata   = dbroot [ tag_datax   ]
hydata   = dbroot [ tag_datay   ]
mctree   = dbroot [ tag_mc      ]
datatree = dbroot [ 'DATA_tree' ]
datastat = datatree.statCov('x','y')
# =============================================================================
## prebook MC histograms
# =============================================================================
ix , iy  = 45 , 25  # #DATA
ix , iy  = 60 , 50 
## ix , iy  = 35 , 22  
## ix , iy  = 60 , 40
hmc  = h2_axes ( [ 20.0/ix*i for i in range ( ix + 1 ) ] ,
                 [ 15.0/iy*i for i in range ( iy + 1 ) ] )

ix , iy  = 60 , 36 
hmcx = h1_axis ( [ 20.0/ix*i for i in range ( ix + 1 ) ] )
hmcy = h1_axis ( [ 15.0/iy*i for i in range ( iy + 1 ) ] )

## prepare re-weighting machinery 
maxIter = 25

## check database 
import os
if not os.path.exists( dbname ) :
    logger.info('Create new weights DBASE') 
    db = DBASE.open ( dbname , 'c' ) ##  create new empty db 
    db.close()
else :
    logger.info('Existing weights DBASE will be used') 
    
# =============================================================================
## make reweighting iterations

from   ostap.tools.reweight         import Weight, makeWeights,  WeightingPlot, W2Data  
from   ostap.fitting.selectors      import SelectorWithVars, Variable 
import ostap.parallel.parallel_fill

# =============================================================================
## configuration of reweighting 
weightings = (
    ## variable          address in DB    
    Weight.Var (  'x'      , 'x-reweight'  ) , 
    Weight.Var (  'y'      , 'y-reweight'  ) , 
    Weight.Var ( ('x','y') , '2D-reweight' ) , 
    )

# =============================================================================
## variables to be used in MC-dataset 
variables  = [
    Variable ( 'x'      , 'x-var'  , 0  , 20 ) , 
    Variable ( 'y'      , 'y-var'  , 0  , 15 ) ,
    ]
selector = SelectorWithVars ( variables , '0<x && x<20 && 0<y && y<20' , silence = True )
mctree.process ( selector , silent = True )
mcds_ = selector.data             ## dataset 
# =============================================================================
## start reweighting iterations:
for iter in range ( 1 , maxIter + 1 ) :

    logger.info ( allright ( 'Reweighting iteration %d ' % iter ) ) 

    with timing ( 'Prepare MC-dataset:' , logger = logger ) : 
        # =========================================================================
        ## 0) The weighter object
        weighter = Weight ( dbname , weightings )
        
        # =========================================================================
        ## 1a) create new "weighted" mcdataset
       mcds=mcds_.clone()

    with timing ( 'Add weight to MC-dataset' , logger = logger ) :
        ## 1b) add  "weight" variable to dataset 
        mcds.add_reweighting ( weighter ,  name = 'weight' ) 
        
        logger.info ('MCDATA:\n%s' %  mcds )
    
    # =========================================================================
    ## 2) update weights
    plots = [ WeightingPlot ( 'y:x' , 'weight' , '2D-reweight' , hdata  , hmc  ) ]
    if 2 < iter: 
        plots  = [
            WeightingPlot ( 'x'     , 'weight' , 'x-reweight'  , hxdata , hmcx        ) ,  
            WeightingPlot ( 'y'     , 'weight' , 'y-reweight'  , hydata , hmcy        ) , 
            WeightingPlot ( 'y:x'   , 'weight' , '2D-reweight' , hdata  , hmc  , 0.99 ) , 
            ]

    with timing ( 'Make one reweighting iteration:' , logger = logger ) : 
        # =========================================================================
        ## 2a) the most important line: perform single iteration step  
        more = makeWeights (
            mcds                                   , ## what to be reweighted
            plots                                  , ## reweighting plots/setup
            dbname                                 , ## DBASE with reweigting constant 
            delta  = 0.04                          , ## stopping criteria
            minmax = 0.08                          , ## stopping criteria  
            power = 2 if 1 != len ( plots ) else 1 , ## tune: effective power
            tag = "Reweight/%d" % iter             ) ## tag for printout
        
    with timing ( 'Project weighted MC-dataset:' , logger = logger ) : 
        # =========================================================================
        ## 3) make MC-histograms  
        mcds .project  ( hmcx , 'x'   , 'weight'  )
        mcds .project  ( hmcy , 'y'   , 'weight'  )
        mcds .project  ( hmc  , 'y:x' , 'weight'  )

    with timing ( 'Compare DATA and MC distributions:' , logger = logger ) :  
        # ==============================================================================
        ## 4) compare "Data" and "MC"  after the reweighting on the given iteration    
        logger.info    ( 'Compare DATA and MC for iteration #%d' % iter )

        hh = 'Iteration#%d: ' % iter 
        ## 4a) compare the basic properties: mean, rms, skewness and kurtosis&moments
        logger.info ( hh + 'DATA(x)  %% MC(x)  comparison:\n%s' % hxdata.cmp_prnt ( hmcx , 'DATA' , 'MC' , 'DATA(x)  vs MC(x)'  , prefix = '# ') ) 
        logger.info ( hh + 'DATA(y)  %% MC(y)  comparison:\n%s' % hydata.cmp_prnt ( hmcy , 'DATA' , 'MC' , 'DATA(y)  vs MC(y)'  , prefix = '# ') ) 
        logger.info ( hh + 'DATA(xy) %% MC(xy) comparison:\n%s' % hdata .cmp_prnt ( hmc  , 'DATA' , 'MC' , 'DATA(xy) vs MC(xy)' , prefix = '# ') ) 
        
        ## 4b) calculate the ``distances''
        logger.info ( hh + "DATA(x)  - MC(x)  ``distance''         %s" % hxdata.cmp_dist ( hmcx , density = True ) )
        logger.info ( hh + "DATA(y)  - MC(y)  ``distance''         %s" % hydata.cmp_dist ( hmcy , density = True ) )
        
        ## 4c) calculate the ``orthogonality''
        logger.info ( hh + "DATA(x)  - MC(x)  ``orthogonality''    %s" % hxdata.cmp_cos  ( hmcx , density = True ) )
        logger.info ( hh + "DATA(y)  - MC(y)  ``orthogonality''    %s" % hydata.cmp_cos  ( hmcy , density = True ) )
        
        ## 4d) get min/max difference between data and MC 
        mn , mx = hxdata.cmp_minmax ( hmcx   , diff = lambda a,b : a/b , density = True )
        logger.info ( hh + "DATA(x)  / MC(x)  ``min/max-distance'' (%s)/(%s)[%%] at xmin/xmax=%.1f/%.1f" % (
            (100*mn[1]-100).toString ( '%+.1f+-%.1f' ) ,
            (100*mx[1]-100).toString ( '%+.1f+-%.1f' ) , mn[0]  , mx[0] ) )
        mn , mx = hydata.cmp_minmax ( hmcy   , diff = lambda a,b : a/b , density = True )
        logger.info ( hh + "DATA(y)  / MC(y)  ``min/max-distance'' (%s)/(%s)[%%] at ymin/ymax=%.1f/%.1f" % (
            (100*mn[1]-100).toString ( '%+.1f+-%.1f' ) ,
            (100*mx[1]-100).toString ( '%+.1f+-%.1f' ) , mn[0]  , mx[0] ) ) 
        mn , mx = hdata .cmp_minmax ( hmc   , diff = lambda a,b : a/b  , density = True )
        logger.info ( hh + "DATA(xy) / MC(xy) ``min/max-distance'' (%s)/(%s)[%%] at (x,y)min/max=(%.1f,%.1f)/(%.1f,%.1f)" % (
            (100*mn[2]-100).toString ( '%+.1f+-%.1f' ) ,
            (100*mx[2]-100).toString ( '%+.1f+-%.1f' ) , mn[0]  , mn[1] , mx[0]  , mx[1] ) )

        h1 = hdata.density()
        h2 = hmc  .density()
        
        logger.info  ('MIN DATA/MC : %s,%s' %( h1 ( mn [0] , mn [1] ), h2 ( mn [0] , mn[1] )))
        logger.info  ('MAX DATA/MC : %s,%s' % (h1 ( mx [0] , mx [1] ), h2 ( mx [0] , mx[1] )))
        
        ## 4e) 2D-statistics 
        mcstat = mcds.statCov('x','y','weight')
        logger.info  ( hh + 'x/y covariance DATA (unbinned):\n# %s' % ( str( datastat [2] ).replace ( '\n' , '\n# ' ) ) )
        logger.info  ( hh + 'x/y covariance MC   (unbinned):\n# %s' % ( str(   mcstat [2] ).replace ( '\n' , '\n# ' ) ) )
        
    # =========================================================================
    ## prepare the plot of weighted MC for the given iteration
    
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
    time.sleep ( 5 )

    if not more and iter > 6 : 
        logger.info    ( allright ( 'No more iterations, converged after #%d' % iter ) )
        break
    
    mcds.clear()
    del mcds
    

else :

    logger.error ( "No convergency!" )

    
del selector   
logger.info ('MCSTAT:\nx=%s\ny=%s\ncov2:\n%s'   %mcstat  [:3] ) 
logger.info ('DATASTAT:\nx=%s\ny=%s\ncov2:\n%s' %datastat[:3] ) 

datax_density.draw ('e1')
mcx_density  .draw ('e1 same')
datay_density.draw ('e1 same')
mcy_density  .draw ('e1 same')
time.sleep(10)

# =============================================================================
# The END 
# =============================================================================
