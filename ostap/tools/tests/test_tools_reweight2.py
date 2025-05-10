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
from   ostap.core.core        import Ostap 
from   ostap.utils.timing     import timing
from   ostap.logger.colorized import attention, allright  
from   ostap.plotting.canvas  import use_canvas
from   ostap.utils.cleanup    import CleanUp
from   ostap.utils.root_utils import batch_env
from   ostap.histos.histos    import h1_axis, h2_axes 
import ostap.io.zipshelve     as     DBASE
import ostap.logger.table     as     T
import  ostap.core.pyrouts    
import ROOT, random, math, os, time 
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
## set batch from environment 
batch_env ( logger )
# =============================================================================


testdata   = CleanUp.tempfile ( suffix = '.root' , prefix ='ostap-test-tools-reweight2-' )
tag_data   = 'DATA2_histogram'
tag_datax  = 'DATAX_histogram'
tag_datay  = 'DATAY_histogram'
tag_mc     = 'MC2_tree'
dbname     = CleanUp.tempfile ( suffix = '.db' , prefix ='ostap-test-tools-reweight2-'   )
 
if os.path.exists ( testdata ) : os.remove ( testdata ) 
if os.path.exists ( dbname   ) : os.remove ( dbname   )

import ostap.parallel.kisa


N1 = 1000000
N2 = 50000

xmax     = 20.0
ymax     = 15.0 

def prepare_data ( ) : 
    #
        
    seed =  1234567890 
    random.seed ( seed ) 
    logger.info ( 'Test *RANDOM* data will be generated/seed=%s' % seed  )   
    ## prepare "data" histograms:
    # 1) 2D histograms
    ix , iy  = 30 , 30
    hdata    = h2_axes ( [ xmax/ix*i for i in range ( ix + 1 ) ] ,
                         [ ymax/iy*i for i in range ( iy + 1 ) ] )
    # 2) non-equal binning 1D histograms for x-component    
    hxdata   = h1_axis ( [    i      for i in range ( 5  ) ] +
                         [  5+i*0.2  for i in range ( 50 ) ] +
                         [ 15+i      for i in range ( 6  ) ] )
    # 2) equal binning 1D histograms for y-component    
    hydata   = h1_axis ( [ i*0.5     for i in range ( 31 ) ] )

    assert hdata .xmax() == xmax , 'XMAX is invalid!'
    assert hdata .ymax() == ymax , 'YMAX is invalid!'
    assert hxdata.xmax() == xmax , 'XMAX is invalid!'
    assert hydata.xmax() == ymax , 'XMAX is invalid!'

    with ROOT.TFile.Open ( testdata ,'recreate') as mc_file:
        
        mc_file.cd() 
        
        datatree  = ROOT.TTree ( 'DATA_tree', 'data-tree' )
        datatree .SetDirectory ( mc_file ) 
        
        from array import array 
        xvar = array  ( 'f', [0])
        yvar = array  ( 'f', [0])
        datatree.Branch ( 'x' , xvar , 'x/F' )
        datatree.Branch ( 'y' , yvar , 'y/F' )
    
        for i in  range ( N1 ) :

            x , y = -1, -1

            while not 0 <= x < xmax or not 0 <= y < ymax :
                
                v1 = random.gauss ( 0 , 3 )
                v2 = random.gauss ( 0 , 2 )
                
                x  =  5 + v1 + v2   
                y  = 12 + v1 - v2
                
            hdata .Fill ( x , y )
            hxdata.Fill ( x )
            hydata.Fill ( y )
 
            xvar [ 0 ] = x 
            yvar [ 0 ] = y
            
            datatree.Fill()

            
        for i in  range ( N1 ) :

            x , y = -1 , -1

            while not 0 <= x < xmax or not 0 <= y < ymax :
                
                v1 = random.gauss ( 0 , 3 )
                v2 = random.gauss ( 0 , 2 )

                x  = 15 + v1 
                y  =  4 + v2    
                
            hdata .Fill ( x , y )
            hxdata.Fill ( x )
            hydata.Fill ( y )
 
            xvar [ 0 ] = x 
            yvar [ 0 ] = y
            
            datatree.Fill()

            
        for i in  range ( N1 ) :

            x , y = -1 , -1

            while not 0 <= x < xmax or not 0 <= y < ymax :
                
                v1 = random.gauss ( 0 , 3 )
                v2 = random.gauss ( 0 , 2 )

                x  = 15 + v1 - v2  
                y  = 12 + v1 + v2     
                
            hdata .Fill ( x , y )
            hxdata.Fill ( x )
            hydata.Fill ( y )
 
            xvar [ 0 ] = x 
            yvar [ 0 ] = y
            
            datatree.Fill()


        for i in range ( 0 ,  4 * N1 ) :
            
            x = random.uniform ( 0 , xmax ) 
            y = random.uniform ( 0 , ymax )
            
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

        for i in  range ( 2 * N2 ) :

            xv = random.uniform ( 0 , xmax ) 
            yv = random.uniform ( 0 , ymax ) 

            xvar [ 0 ] = xv 
            yvar [ 0 ] = yv
            
            mctree.Fill()
            
        fx = lambda x : ( x/10.0-1 ) ** 2
        fy = lambda y : ( y/ 7.5-1 ) ** 2
        

        for i in  range ( N2 ) :
            
            while True : 
                xv = random.uniform ( 0 , xmax )
                if random.uniform   ( 0 ,  1 ) < fx ( xv ) : break
                
            while True : 
                yv = random.uniform ( 0 , ymax )
                if random.uniform   ( 0 ,  1 ) < fy ( yv ) : break 
                
            xvar [ 0 ] = xv 
            yvar [ 0 ] = yv
            
            mctree.Fill()

            
        mctree.Write()
        mc_file.ls()


if not os.path.exists( testdata ) :
    with timing ( "Prepare input data" , logger = logger ) :
        prepare_data ()
        
# =============================================================================
## Read data from DB
# =============================================================================
with ROOT.TFile.open ( testdata , 'r' ) as dbroot : 
    logger.info ( 'Test data is fetched from DBASE "%s"' % testdata )   
    dbroot.ls()
    hdata    = dbroot [ tag_data  ].clone() 
    hxdata   = dbroot [ tag_datax ].clone()
    hydata   = dbroot [ tag_datay ].clone()
    
datatree = ROOT.TChain ( 'DATA_tree' ) ; datatree.Add ( testdata ) 
mctree   = ROOT.TChain ( tag_mc      ) ; mctree.Add   ( testdata ) 
datastat = datatree.statCov ( 'x' , 'y' )        
# =============================================================================
## prebook MC histograms
# =============================================================================
ix , iy  = 45 , 25  # #DATA
ix , iy  = 25 , 25 
## ix , iy  = 35 , 22
## ix , iy  = 60 , 40
hmc  = h2_axes ( [ 20.0/ix*i for i in range ( ix + 1 ) ] ,
                 [ 15.0/iy*i for i in range ( iy + 1 ) ] )
## hmc = hdata.clone()

ix , iy  = 65 , 60
ix , iy  = 50 , 50
hmcx = h1_axis ( [ 20.0/ix*i for i in range ( ix + 1 ) ] )
hmcy = h1_axis ( [ 15.0/iy*i for i in range ( iy + 1 ) ] )

assert hmc .xmax() == xmax , 'XMAX is invalid!'
assert hmc .ymax() == ymax , 'YMAX is invalid!'
assert hmcx.xmax() == xmax , 'XMAX is invalid!'
assert hmcy.xmax() == ymax , 'XMAX is invalid!'

## prepare re-weighting machinery 
maxIter = 15

## check database 
import os
if not os.path.exists( dbname ) :
    with DBASE.open ( dbname , 'c' ) : ##  create new empty db 
        logger.info("Create new weights DBASE '%s'" % dbname ) 
else :
    logger.info("Existing weights DBASE '%s' will be used" % dbname ) 
    
# =============================================================================
## make reweighting iterations

from   ostap.tools.reweight           import Weight, makeWeights,  WeightingPlot, W2Data  
from   ostap.fitting.pyselectors      import Variable 
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

with timing ( 'Prepare initial MC-dataset:' , logger = logger ) :

    mctree    = ROOT.TChain ( tag_mc      ) ; mctree.Add   ( testdata ) 
    mcds_ , _ = mctree.make_dataset ( variables = variables      ,
                                      selection = '0<x && x<20 && 0<y && y<20' ,
                                      silent    = True           ) 

# =============================================================================
## Configurtaion of reweighting plots 
# =
plots  = [
    WeightingPlot ( 'x'     , 'weight' , 'x-reweight'  , hxdata , hmcx ) ,  
    WeightingPlot ( 'y'     , 'weight' , 'y-reweight'  , hydata , hmcy ) , 
    WeightingPlot ( 'x:y'   , 'weight' , '2D-reweight' , hdata  , hmc  ) , 
    ]

# ============================================================================
## table of global statistics 
glob_stat   = [ ( '#' , 'Mahalanobis' , 'Hotelling' , 'KL/S-C' , 'KL/C-S' , 'KL-sym' ) ]
n_data      = len ( datatree ) 
n_mc        = len ( mctree   )
##
datatree = ROOT.TChain ( 'DATA_tree' ) ; datatree.Add ( testdata ) 
mctree   = ROOT.TChain ( tag_mc      ) ; mctree.Add   ( testdata ) 
##
vct_data    = datatree.statVct ( 'x,y' )
vct_init    = mctree  .statVct ( 'x,y' )
trow = ( '%d'    % 0 ,
         '%.4g'  % vct_init.mahalanobis                 ( vct_data ) ,
         '%.4g'  % Ostap.Math.hotelling ( vct_init, n_mc , vct_data , n_data ) ,         
         '%+.4g' % vct_init.asymmetric_kullback_leibler ( vct_data ) , 
         '%+.4g' % vct_data.asymmetric_kullback_leibler ( vct_init ) , 
         '%+.4g' % vct_data.           kullback_leibler ( vct_init ) )
glob_stat.append ( trow )
# =============================================================================

# =============================================================================
## start reweighting iterations:
for iter in range ( 1 , maxIter + 1 ) :

    tag = 'Reweighting iteration #%d' % iter
    logger.info ( allright ( tag ) ) 
    
    with timing ( tag + ': prepare MC-dataset:' , logger = logger ) : 
        # =========================================================================
        ## 0) The weighter object
        weighter = Weight ( dbname , weightings )
        
        # =========================================================================
        ## 1a) create new "weighted" mcdataset
        mcds = mcds_.Clone()

    with timing ( tag + ': add weight to MC-dataset' , logger = logger ) :
        ## 1b) add  "weight" variable to dataset 
        mcds.add_reweighting ( weighter ,  name = 'weight' ) 
        if 1 == iter % 10  : logger.info ( ( tag + ' MCDATA:\n%s' ) %  mcds )
    
    # =========================================================================
    ## 2) update weights
    
    if    iter <=  3 : power = 0.75
    elif  iter <= 10 : power = lambda nactive : 1.5 / nactive if 1 < nactive else 1.05
    elif  iter <= 15 : power = lambda nactive : 1.3 / nactive if 1 < nactive else 1.05 
    else             : power = lambda nactive : 1.1 / nactive if 1 < nactive else 1.05 
    
    with timing ( tag + ': make actual reweighting:' , logger = logger ) :
        
        # =========================================================================
        ## 2a) the most important line: perform single iteration step  
        active , cmp_plots  = makeWeights (
            mcds               , ## what to be reweighted
            plots              , ## reweighting plots/setup
            dbname             , ## DBASE with reweigting constant 
            delta      = 0.02  , ## stopping criteria
            minmax     = 0.10  , ## stopping criteria  
            power      = power , ## tune: effective power
            make_plots = True  , 
            tag        = tag   ) ## tag for printout
        
    with timing ( tag + ': project weighted MC-dataset:' , logger = logger ) : 
        # =========================================================================
        ## 3) make MC-histograms  
        mcds .project  ( hmcx , 'x'   , 'weight'  )
        mcds .project  ( hmcy , 'y'   , 'weight'  )
        mcds .project  ( hmc  , 'x:y' , 'weight'  )
        
        ## 3.1) compare control and signal samples  
        vct_i = mcds.statVct ( 'x,y' , 'weight' )
        n_mc  = int ( mcds.nEff('weight') ) 
        trow = ( '%d'    % iter ,
                 '%.4g'  % vct_i   .mahalanobis                 ( vct_data ) ,
                 '%.4g'  % Ostap.Math.hotelling ( vct_i , n_mc  , vct_data , n_data ) ,         
                 '%+.4g' % vct_i   .asymmetric_kullback_leibler ( vct_data ) , 
                 '%+.4g' % vct_data.asymmetric_kullback_leibler ( vct_i    ) , 
                 '%+.4g' % vct_data.           kullback_leibler ( vct_i    ) )
        glob_stat.append ( trow )
  
    rows = [] 
    with timing ( tag + ': compare DATA and MC distributions:' , logger = logger ) :
        
        # ==============================================================================
        ## 4) compare "Data" and "MC"  after the reweighting on the given iteration    
        logger.info    ( tag + ': compare DATA and MC for iteration #%d' % iter )
        
        ## 4a) compare the basic properties: mean, rms, skewness and kurtosis&moments
        logger.info ( tag + ': DATA(x)   vs MC(x)  comparison:\n%s'  % hxdata.cmp_prnt ( hmcx , 'DATA' , 'MC' , 'DATA(x)  vs MC(x)'  , prefix = '# ' , density = True ) )
        logger.info ( tag + ': DATA(y)   vs MC(y)  comparison:\n%s'  % hydata.cmp_prnt ( hmcy , 'DATA' , 'MC' , 'DATA(y)  vs MC(y)'  , prefix = '# ' , density = True ) )        
        logger.info ( tag + ': DATA(x,y) vs MC(x,y) comparison:\n%s' % hdata .cmp_prnt ( hmc  , 'DATA' , 'MC' , 'DATA(xy) vs MC(xy)' , prefix = '# ' , density = True ) )
                
        title = tag + ': DATA(x)   vs MC(x) difference'
        logger.info ( '%s:\n%s' % ( title , hxdata.cmp_diff_prnt ( hmcx , density = True , title = title , prefix = '# ' ) ) )
        
        title = tag + ': DATA(y)   vs MC(y) difference'
        logger.info ( '%s:\n%s' % ( title , hydata.cmp_diff_prnt ( hmcy , density = True , title = title , prefix = '# ' ) ) )
        
        title = tag + ': DATA(x,y) vs MC(x,y) difference'
        logger.info ( '%s:\n%s' % ( title , hdata .cmp_diff_prnt ( hmc  , density = True , title = title , prefix = '# ' ) ) )
        
        ## 4e) 2D-statistics
        mcstat = mcds    .statCov ( 'x' , 'y' , 'weight' )
        logger.info  ( tag + ': x/y correlation DATA (unbinned): %+.2f' % datastat.correlation () ) 
        logger.info  ( tag + ': x/y correlation MC   (unbinned): %+.2f' %   mcstat.correlation () )

    if not active and 3 < iter : 
        logger.info    ( allright ( 'No more iterations, converged after #%d' % iter ) )
        title = 'Reweighted dataset after #%d iterations' % iter 
        logger.info ( '%s:\n%s' % ( title , mcds.table2 ( variables = [ 'x' , 'y' ] ,
                                                          title     = title    ,
                                                          cuts      = 'weight' , 
                                                          prefix    = '# '     ) ) )
        break
    
    mcds.clear()
    del mcds
    
else :

    logger.error ( "No convergency!" )


# ===========================================================================
title = 'Weighter object'
logger.info ( '%s:\n%s' % ( title , weighter.table ( prefix = '# ' ) ) )
# ============================================================================
## draw the convergency graphs 
graphs = weighter.graphs ()
for key in graphs : 
    with use_canvas ( "Convergency graph for '%s'" % key ) :
        graph = graphs [ key ]
        graph.draw ( 'a' )
# =============================================================================
    
with ROOT.TFile.open ( testdata , 'r' ) as dbroot : 
    logger.info ( 'Test data is picked from DBASE "%s" for reweighting' % testdata )   
    dbroot.ls()
    mctree   = dbroot [ tag_mc ]    
    ## 0) The weighter object
    with timing ( "Add weight column to MC-tree" , logger = logger ) : 
        weighter = Weight ( dbname , weightings )
        mctree.add_reweighting ( weighter ,  name = 'weight' )

# =============================================================================
data_tree = ROOT.TChain  ( 'DATA_tree' ) ; data_tree.Add ( testdata ) 
mc_tree   = ROOT.TChain  ( tag_mc      ) ;   mc_tree.Add ( testdata ) 

# =============================================================================
vct_final = mc_tree.statVct ( 'x,y' , 'weight' )
n_mc = int ( mc_tree.nEff('weight')  )
trow = ( '*' ,
         '%.4g'  % vct_final .mahalanobis                 ( vct_data  ) , 
         '%.4g'  % Ostap.Math.hotelling ( vct_final , n_mc ,  vct_data  , n_data ) ,         
         '%+.4g' % vct_final .asymmetric_kullback_leibler ( vct_data  ) , 
         '%+.4g' % vct_data  .asymmetric_kullback_leibler ( vct_final ) , 
         '%+.4g' % vct_data  .           kullback_leibler ( vct_final ) )
glob_stat.append ( trow )

title = 'Global DATA/MC similarity '
table = T.table ( glob_stat , title = title , prefix = '# ' ) 
logger.info ( '%s\n%s' % ( title , table ) ) 

# =============================================================================
title = 'Data/target dataset'
logger.info ( '%s:\n%s' % ( title , data_tree.table2 ( variables = [ 'x' , 'y' ] ,
                                                       title     = title    ,
                                                       prefix    = '# '     ) ) )
# =============================================================================
title = 'MC dataset before reweighting' 
logger.info ( '%s:\n%s' % ( title , mc_tree.table2   ( variables = [ 'x' , 'y' ] ,
                                                       title     = title    ,
                                                       prefix    = '# '     ) ) )
# =============================================================================
title = 'MC dataset after reweighting' 
logger.info ( '%s:\n%s' % ( title , mc_tree.table2   ( variables = [ 'x' , 'y' ] ,
                                                       title     = title    ,
                                                       cuts      = 'weight' , 
                                                       prefix    = '# '     ) ) )

# =============================================================================
##                                                                      The END 
# =============================================================================
