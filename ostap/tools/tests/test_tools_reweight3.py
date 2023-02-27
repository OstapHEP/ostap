1#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/tests/test_tools_reweight3.py
#  Test for 3D-reweighting machinery
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2023-01-20
# =============================================================================
"""Test for 3D-reweighting machinery
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2023-01-20"
__all__     = ()  ## nothing to be imported 
# =============================================================================
from   builtins               import range
from   ostap.core.pyrouts     import Ostap 
import ostap.io.root_file 
import ostap.parallel.kisa
import ostap.io.zipshelve     as     DBASE
from   ostap.histos.histos    import h1_axis, h2_axes, h3_axes 
from   ostap.utils.timing     import timing
from   ostap.logger.colorized import attention, allright  
import ROOT, random, math, os, time 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_tools_reweight3' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================    
logger.info ( 'Test for 3D-Reweighting machinery')
# ============================================================================
from ostap.utils.cleanup import CleanUp
testdata   = CleanUp.tempfile ( suffix = '.root' , prefix ='ostap-test-tools-reweight3-' )

tag_data_3d = 'DATA_3D_histogram'

tag_data_xy = 'DATA_XY_histogram'
tag_data_yz = 'DATA_YZ_histogram'
tag_data_xz = 'DATA_XZ_histogram'

tag_data_x  = 'DATA_X_histogram'
tag_data_y  = 'DATA_Y_histogram'
tag_data_z  = 'DATA_Z_histogram'

tag_data    = 'DATA_tree'
tag_mc      = 'MC_tree'

dbname      = CleanUp.tempfile ( suffix = '.db' , prefix ='ostap-test-tools-reweight3-'   )

NDATA1      = 500000
NDATA2      = 400000
NDATA3      = 400000
NMC         = 500000

xmax        = 15.0
ymax        = 12.0 
zmax        = 10.0 

def prepare_data ( ) : 
    #
        
    seed =  1234567890 
    random.seed ( seed ) 
    logger.info ( 'Test *RANDOM* data will be generated/seed=%s' % seed  )   
    ## prepare "data" histograms:

    # 1) 3D histogram
    ix , iy , iz  = 5 , 5 , 5 
    h3_histo  = h3_axes ( [ xmax/ix*i for i in range ( ix + 1 ) ] ,
                          [ ymax/iy*i for i in range ( iy + 1 ) ] ,
                          [ zmax/iz*i for i in range ( iz + 1 ) ] )
    
    # 2) XY,YX &XZ histograms

    ix , iy , iz  =  7 ,  7 ,  7
    hxy_histo = h2_axes ( [ xmax/ix*i for i in range ( ix + 1 ) ] ,
                          [ ymax/iy*i for i in range ( iy + 1 ) ] )

    ix , iy , iz  =  8 ,  8 ,  8 
    hyz_histo = h2_axes ( [ ymax/iy*i for i in range ( iy + 1 ) ] ,
                          [ zmax/iz*i for i in range ( iz + 1 ) ] )

    ix , iy , iz  =  9 ,  9 ,  9 
    hxz_histo = h2_axes ( [ xmax/ix*i for i in range ( ix + 1 ) ] ,
                          [ zmax/iz*i for i in range ( iz + 1 ) ] )
    
    ## 3) X,Y,Z histogramss
    
    ix , iy , iz  = 100 , 100 , 100 
    hx_histo  = h1_axis ( [ xmax/ix*i for i in range ( ix + 1 ) ] )
    hy_histo  = h1_axis ( [ ymax/iy*i for i in range ( iy + 1 ) ] )
    hz_histo  = h1_axis ( [ zmax/iz*i for i in range ( iz + 1 ) ] )

    with ROOT.TFile.Open ( testdata ,'recreate') as mc_file:
        
        mc_file.cd() 
        
        datatree  = ROOT.TTree ( tag_data , 'data-tree' )
        datatree .SetDirectory ( mc_file ) 
        
        from array import array 
        xvar = array    ( 'f', [0])
        yvar = array    ( 'f', [0])
        zvar = array    ( 'f', [0])
        datatree.Branch ( 'x' , xvar , 'x/F' )
        datatree.Branch ( 'y' , yvar , 'y/F' )
        datatree.Branch ( 'z' , zvar , 'z/F' )

        ## Gaussian 3D-component with correlations 
        
        g3d = Ostap.Math.Gauss3D ( 0.1 * xmax  , 0.5 * ymax  , 0.8 * zmax  ,
                                   8           , 5           , 4           ,
                                   math.pi / 5 , math.pi / 5 , math.pi / 5 )
        
        for x, y , z in g3d.random ( NDATA1 , 0 , xmax , 0 , ymax , 0 , zmax ) : 

            h3_histo .Fill ( x , y , z )
            
            hxy_histo.Fill ( x , y )
            hyz_histo.Fill ( y , z )
            hxz_histo.Fill ( x , z )

            hx_histo .Fill ( x )
            hy_histo .Fill ( y )
            hz_histo .Fill ( z )
            
            xvar [ 0 ] = x 
            yvar [ 0 ] = y
            zvar [ 0 ] = z 
            
            datatree.Fill()

        for i in range ( NDATA2 ) :
            
            x , y, z  = -1, -1, -1 

            while not 0 < x < xmax : x = xmax - random.expovariate ( 1.0 / xmax ) 
            while not 0 < y < ymax : y =        random.expovariate ( 1.0 / ymax ) 
            while not 0 < z < zmax : z =        random.expovariate ( 1.0 / zmax ) 
                
            h3_histo .Fill ( x , y , z )
            
            hxy_histo.Fill ( x , y )
            hyz_histo.Fill ( y , z )
            hxz_histo.Fill ( x , z )

            hx_histo .Fill ( x )
            hy_histo .Fill ( y )
            hz_histo .Fill ( z )
            
            xvar [ 0 ] = x 
            yvar [ 0 ] = y
            zvar [ 0 ] = z 
            
            datatree.Fill()

        for i in range ( NDATA3 ) :
            
            x , y, z  = -1, -1, -1 

            x = random.uniform ( 0 , xmax ) 
            y = random.uniform ( 0 , ymax ) 
            z = random.uniform ( 0 , zmax ) 
                
            h3_histo .Fill ( x , y , z )
            
            hxy_histo.Fill ( x , y )
            hyz_histo.Fill ( y , z )
            hxz_histo.Fill ( x , z )

            hx_histo .Fill ( x )
            hy_histo .Fill ( y )
            hz_histo .Fill ( z )
            
            xvar [ 0 ] = x 
            yvar [ 0 ] = y
            zvar [ 0 ] = z 
            
            datatree.Fill()
            

        datatree.Write()
        
        ## write the histogram 
        mc_file [ tag_data_3d  ] = h3_histo
        
        mc_file [ tag_data_xy  ] = hxy_histo
        mc_file [ tag_data_yz  ] = hyz_histo
        mc_file [ tag_data_xz  ] = hxz_histo

        mc_file [ tag_data_x   ] = hx_histo
        mc_file [ tag_data_y   ] = hy_histo
        mc_file [ tag_data_z   ] = hz_histo

        ## 
        mctree  = ROOT.TTree ( tag_mc , 'mc-tree' )
        mctree .SetDirectory ( mc_file ) 
        
        from array import array 
        xvar = array  ( 'f', [ 0.0 ] )
        yvar = array  ( 'f', [ 0.0 ] )
        zvar = array  ( 'f', [ 0.0 ] )

        mctree.Branch ( 'x' , xvar , 'x/F' )
        mctree.Branch ( 'y' , yvar , 'y/F' )
        mctree.Branch ( 'z' , zvar , 'z/F' )
        
        for i in  range ( NMC ) :

            xv = random.uniform  ( 0 , xmax )
            yv = random.uniform  ( 0 , ymax )
            zv = random.uniform  ( 0 , zmax )
            
            xvar [ 0 ] = xv 
            yvar [ 0 ] = yv
            zvar [ 0 ] = zv
            
            mctree.Fill()
            
        mctree.Write()

if not os.path.exists( testdata ) :
    with timing ( "Prepare input data" , logger = logger ) :
        prepare_data ()
        
# =============================================================================
## Read data from DB
# =============================================================================
with ROOT.TFile.open ( testdata , 'r' ) as dbroot :
    
    logger.info ( 'Test data is fetched from DBASE "%s"' % testdata )   
    dbroot.ls_table ()
    
    h3d_data = dbroot [ tag_data_3d ] . clone()
    hxy_data = dbroot [ tag_data_xy ] . clone()
    hyz_data = dbroot [ tag_data_yz ] . clone() 
    hxz_data = dbroot [ tag_data_xz ] . clone() 
    hx_data  = dbroot [ tag_data_x  ] . clone()
    hy_data  = dbroot [ tag_data_y  ] . clone()
    hz_data  = dbroot [ tag_data_z  ] . clone()


# =============================================================================
## prebook MC histograms
# =============================================================================
ix , iy , iz  = 5 , 5 , 5 
mc_3d = h3_axes ( [ xmax/ix*i for i in range ( ix + 1 ) ] ,
                  [ ymax/iy*i for i in range ( iy + 1 ) ] ,
                  [ zmax/iz*i for i in range ( iz + 1 ) ] )
ix , iy , iz  =  8 ,  8 ,  8 
mc_xy = h2_axes ( [ xmax/ix*i for i in range ( ix + 1 ) ] ,
                  [ ymax/iy*i for i in range ( iy + 1 ) ] )
ix , iy , iz  =  9 ,  9 ,  9 
mc_yz = h2_axes ( [ ymax/iy*i for i in range ( iy + 1 ) ] ,
                  [ zmax/iz*i for i in range ( iz + 1 ) ] )
ix , iy , iz  =  7 ,  7 ,  7 
mc_xz = h2_axes ( [ xmax/ix*i for i in range ( ix + 1 ) ] ,
                  [ zmax/iz*i for i in range ( iz + 1 ) ] )
ix , iy , iz  = 95 , 95 , 95 
mc_x  = h1_axis ( [ xmax/ix*i for i in range ( ix + 1 ) ] )
mc_y  = h1_axis ( [ ymax/iy*i for i in range ( iy + 1 ) ] )
mc_z  = h1_axis ( [ zmax/iz*i for i in range ( iz + 1 ) ] )


## prepare re-weighting machinery 
maxIter = 15

## check database 
if not os.path.exists( dbname ) :
    logger.info('Create new weights DBASE') 
    with DBASE.open ( dbname , 'c' ) : ##  create new empty db 
        pass 
else :
    logger.info('Existing weights DBASE will be used') 
    with DBASE.open ( dbname , 'r' ) as db : 
        db.ls() 
        
# =============================================================================
## make reweighting iterations

from   ostap.tools.reweight           import Weight, makeWeights,  WeightingPlot, W2Data  
from   ostap.fitting.pyselectors      import Variable 
import ostap.parallel.parallel_fill

# =============================================================================
## configuration of reweighting 
weightings = (
    ## variable          address in DB    
    Weight.Var (  'x'          , 'x-reweight'  ) , 
    Weight.Var (  'y'          , 'y-reweight'  ) , 
    Weight.Var (  'z'          , 'z-reweight'  ) , 
    Weight.Var ( ('x','y')     , 'XY-reweight' ) , 
    Weight.Var ( ('y','z')     , 'YZ-reweight' ) , 
    Weight.Var ( ('x','z')     , 'XZ-reweight' ) , 
    Weight.Var ( ('x','y','z') , '3D-reweight' ) , 
    )

# =============================================================================
## variables to be used in MC-dataset 
variables  = [
    Variable ( 'x' , 'x-var'  , 0  , xmax ) , 
    Variable ( 'y' , 'y-var'  , 0  , ymax ) ,
    Variable ( 'z' , 'z-var'  , 0  , zmax ) ,
    ]


# =============================================================================
datatree   = ROOT.TChain ( tag_data ) ; datatree.Add ( testdata )  
title = 'Data/target dataset'
logger.info ( '%s:\n%s' % ( title , datatree.table2 ( variables = [ 'x' , 'y' , 'z' ] ,
                                                      title     = title    ,
                                                      prefix    = '# '     ) ) )

# =============================================================================
with timing ( 'Prepare initial MC-dataset:' , logger = logger ) :
    mctree   = ROOT.TChain ( tag_mc   ) ; mctree   .Add ( testdata )
    ## fill dataset from input MC tree
    mcds_ , _ = mctree.make_dataset ( variables = variables ,
                                      selection = '0<x && x<%s && 0<y && y<%s && 0<z && z<%s' % ( xmax , ymax, zmax ) )
                                  
##    selector = SelectorWithVars (
##        variables ,
##        '0<x && x<%s && 0<y && y<%s && 0<z && z<%s' % ( xmax , ymax, zmax ) , silence = True )
##    mctree.process ( selector , silent = True )
##    mcds_ = selector.data             ## dataset


# =============================================================================
## Configuration of reweighting plots 
# =============================================================================
plots  = [
    WeightingPlot ( 'z:y:x' , 'weight' , '3D-reweight' , h3d_data , mc_3d ) , 
    WeightingPlot ( 'z:x'   , 'weight' , 'XZ-reweight' , hxz_data , mc_xz ) , 
    WeightingPlot ( 'z:y'   , 'weight' , 'YZ-reweight' , hyz_data , mc_yz ) , 
    WeightingPlot ( 'y:x'   , 'weight' , 'XY-reweight' , hxy_data , mc_xy ) , 
    WeightingPlot ( 'x'     , 'weight' , 'x-reweight'  , hx_data  , mc_x  ) ,  
    WeightingPlot ( 'z'     , 'weight' , 'z-reweight'  , hz_data  , mc_z  ) ,  
    WeightingPlot ( 'y'     , 'weight' , 'y-reweight'  , hy_data  , mc_y  ) ,  
    ]

# =============================================================================
## start reweighting iterations:
for iter in range ( 1 , maxIter + 1 ) :

    tag = 'Reweighting iteration #%d' % iter
    logger.info ( allright ( tag ) ) 
    
    # =========================================================================
    ## 0) The weighter object
    weighter = Weight ( dbname , weightings )
    ## 1a) create new "weighted" mcdataset
    mcds = mcds_.Clone()
    
    with timing ( tag + ': add weight to MC-dataset' , logger = logger ) :
        # =========================================================================
        ## 1b) add  "weight" variable to dataset 
        mcds.add_reweighting ( weighter ,  name = 'weight' , progress = True )
        title = 'Reweighted dataset at iteration #%d' % iter 
        logger.info ( '%s:\n%s' % ( title , mcds.table2 ( variables = [ 'x' , 'y' , 'z' ] ,
                                                          title     = title    ,
                                                          cuts      = 'weight' , 
                                                          prefix    = '# '     ) ) )        
    
    # =========================================================================
    ## 2) update weights

    pmax  = 1.0 if iter <= 4 else 0.66 + 0.34 * math.exp ( - ( iter - 4.0 ) / 4.0 )
    power = lambda n : pmax if n <= 1 else pmax * ( 0.66666 / n + 0.33333 )

    if   iter < 3 : the_plots = plots [ : 1 ]
    elif iter < 4 : the_plots = plots [ : 4 ]
    else          : the_plots = plots

    ## weigth truncation: avoid very large change of weights for  single iteration 
    wtruncate = ( 0.1 , 10 ) if iter < 4 else ( 0.5 , 2.0 )  
    
    with timing ( tag + ': make actual reweighting:' , logger = logger ) :
        
        # =========================================================================
        ## 2a) the most important line: perform single iteration step  
        active , _ = makeWeights (
            mcds                   , ## what to be reweighted
            the_plots              , ## reweighting plots/setup
            dbname                 , ## DBASE with reweighting constant 
            delta      = 0.02      , ## stopping criteria
            minmax     = 0.08      , ## stopping criteria  
            power      = power     , ## tune: effective power
            wtruncate  = wtruncate , ## truncate weights 
            make_plots = True      , ## make control plots 
            tag        = tag       ) ## tag for printout

    if not active and 5 <= iter : 
        logger.info    ( allright ( 'No more iterations, converged after #%d' % iter ) )
        break
    
    mcds.clear()
    del mcds
    
else :

    logger.error ( "No convergency!" )

# =============================================================================
with timing ( "Add weight column to initial MC-tree" , logger = logger ) : 
    mctree   = ROOT.TChain ( tag_mc   ) ; mctree   .Add ( testdata )  
    weighter = Weight ( dbname , weightings )
    mctree   = mctree.add_reweighting ( weighter ,  name = 'weight' )
    mctree   = ROOT.TChain ( tag_mc   ) ; mctree   .Add ( testdata )  

# =============================================================================
## compare DATA  and MC before and after reweighting
# =============================================================================

datatree   = ROOT.TChain ( tag_data ) ; datatree.Add ( testdata )  
mctree     = ROOT.TChain ( tag_mc   ) ; mctree  .Add ( testdata )  
# =============================================================================
title = 'Data/target dataset'
logger.info ( '%s:\n%s' % ( title , datatree.table2 ( variables = [ 'x' , 'y' , 'z' ] ,
                                                      title     = title    ,
                                                      prefix    = '# '     ) ) )
# =============================================================================
title = 'MC-tree before reweighting' 
logger.info ( '%s:\n%s' % ( title , mctree.table2   ( variables = [ 'x' , 'y' , 'z' ] ,
                                                      title     = title    ,
                                                      prefix    = '# '     ) ) )
# =============================================================================
title = 'MC-tree after reweighting' 
logger.info ( '%s:\n%s' % ( title , mctree.table2   ( variables = [ 'x' , 'y' , 'z' ] ,
                                                      title     = title    ,
                                                      cuts      = 'weight' , 
                                                      prefix    = '# '     ) ) )

# =============================================================================
##                                                                      The END 
# =============================================================================
