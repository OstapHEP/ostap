#!/usr/bin/env python
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
from   ostap.core.pyrouts       import Ostap, VE  
from   ostap.histos.histos      import h1_axis, h2_axes, h3_axes 
from   ostap.utils.timing       import timing
from   ostap.logger.colorized   import attention, allright
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.root_utils   import batch_env 
from   ostap.utils.cleanup      import CleanUp
from   ostap.logger.symbols     import iteration, plus_minus 
from   ostap.utils.memory       import memory_usage, delta_ram
from   ostap.utils.basic        import numcpu 
from   ostap.utils.progress_bar import progress_bar 
import ostap.io.zipshelve       as     DBASE
import ostap.logger.table       as     T 
import ostap.logger.table       as     T 
import ostap.io.root_file 
import ostap.parallel.kisa
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
## set batch from environment 
batch_env ( logger )
# =============================================================================
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
dbname      = 'reweight3.db'

if numcpu () <= 8 : 
    NDATA1      = 2000 ## 00
    NDATA2      = 2000 ## 00
    NDATA3      = 2000 ## 00
    NMC         = 2000 ## 00
else : 
    NDATA1      =  5000 ## 00
    NDATA2      =  5000 ## 00
    NDATA3      =  5000 ## 00
    NMC         = 10000 ## 00
    
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
        xvar = array    ( 'f',  [ 0 ] )
        yvar = array    ( 'f',  [ 0 ] )
        zvar = array    ( 'f',  [ 0 ] )
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

        for i in progress_bar ( range ( NDATA2 ) ) :
            
            x , y, z  = -1, -1, -1 

            while not 0 <= x < xmax : x = xmax - random.expovariate ( 1.0 / xmax ) 
            while not 0 <= y < ymax : y =        random.expovariate ( 1.0 / ymax ) 
            while not 0 <= z < zmax : z =        random.expovariate ( 1.0 / zmax ) 
                
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

        for i in progress_bar ( range ( NDATA3 ) ) :
            
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
        
        for i in progress_bar (  range ( NMC ) ) :

            xv = random.uniform  ( 0 , xmax )
            yv = random.uniform  ( 0 , ymax )
            zv = random.uniform  ( 0 , zmax )
            
            xvar [ 0 ] = xv 
            yvar [ 0 ] = yv
            zvar [ 0 ] = zv
            
            mctree.Fill()
            
        mctree.Write()

if not os.path.exists ( testdata ) :
    with timing ( "Prepare input data" , logger = logger ) :
        prepare_data ()

# ==============================================================================
## Compare datasets using several methods 
# ==============================================================================
from ostap.stats.adval        import ADVAL_LGBM  as COMPARATOR5
from ostap.stats.gof_np       import  ( MIXnp           as COMPARATOR4 , 
                                        KullbackLeibler as COMPARATOR3 , 
                                        Mahalanobis     as COMPARATOR2 , 
                                        Hotelling       as COMPARATOR1 )

from ostap.stats.data_compare import data_compare     
comparators = ( COMPARATOR1 ( parallel = True , nToys = 100 ) ,
                COMPARATOR2 ( parallel = True , nToys = 100 ) ,
                COMPARATOR3 ( parallel = True , nToys = 100 ) ,
                COMPARATOR4 ( parallel = True , nToys = 100 ) , 
                COMPARATOR5 ( parallel = True , nToys = 100 ) )

if numcpu () <= 8 : comparators = comparators[ :3 ]
    

# ============================================================================
## The table of global comparison statistics 
header    = ( '#' , '#eff' ) + tuple ( c.method for c in comparators ) 
glob_stat = [ header ]
alignment = 'lc' + 'c' * len ( comparators )           

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
    logger.attention ( "Create new weights DBASE `%s'" % dbname ) 
    with DBASE.open ( dbname , 'c' ) : ##  create new empty db 
        pass 
else :
    logger.attention ( "Existing weights DBASE `%s' will be (re)used" % dbname ) 
    with DBASE.open ( dbname , 'r' ) as db : 
        db.ls() 

# =============================================================================
## name of weight variable 
weight_name = 'weight'
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
title      = 'Data/target dataset'
logger.info ( '%s:\n%s' % ( title , datatree.table2 ( variables = [ 'x' , 'y' , 'z' ] ,
                                                      title     = title    ,
                                                      prefix    = '# '     ) ) )

# =============================================================================
with timing ( 'Prepare initial MC-dataset:' , logger = logger ) :
    mctree   = ROOT.TChain ( tag_mc   ) ; mctree   .Add ( testdata )
    ## fill dataset from input MC tree
    mcds_init , _ = mctree.make_dataset ( variables = variables ,
                                          selection = '0<=x && x<%s && 0<=y && y<%s && 0<=z && z<%s' % ( xmax , ymax, zmax ) )

# =============================================================================
## Configuration of reweighting plots 
# =============================================================================
plots  = [
    WeightingPlot ( 'x:y:z' , address = '3D-reweight' , data = h3d_data , mc = mc_3d ) , 
    WeightingPlot ( 'x:z'   , address = 'XZ-reweight' , data = hxz_data , mc = mc_xz ) , 
    WeightingPlot ( 'y:z'   , address = 'YZ-reweight' , data = hyz_data , mc = mc_yz ) , 
    WeightingPlot ( 'x:y'   , address = 'XY-reweight' , data = hxy_data , mc = mc_xy ) , 
    WeightingPlot ( 'x'     , address = 'x-reweight'  , data = hx_data  , mc = mc_x  ) ,  
    WeightingPlot ( 'z'     , address = 'z-reweight'  , data = hz_data  , mc = mc_z  ) ,  
    WeightingPlot ( 'y'     , address = 'y-reweight'  , data = hy_data  , mc = mc_y  ) ,  
    ]

memory_init = memory_usage() 
converged   = False

maxIter = 2
# =============================================================================
## start reweighting iterations:
for iter in range ( 1 , maxIter + 1 ) :
    
    tag = 'Reweighting iteration #%d%s' %  ( iter , iteration )
    mem = ''
    if 1 < iter : mem = ' Memory:%s=%+.2f[MB]' % ( delta_ram , memory_usage () - memory_init )
    logger.info ( allright ( tag + mem ) )
    
    # =========================================================================
    ## 0) The weighter object
    weighter = Weight ( dbname , weightings )
    ## 1a) create new "weighted" mcdataset
    mcds = mcds_init.copy()
    
    with timing ( tag + ': add weight to MC-dataset' , logger = logger ) :
        # =========================================================================
        ## 1b) add  "weight" variable to dataset 
        mcds  = mcds.add_reweighting ( weighter ,  name = 'weight' , progress = True )
        ## 1c) make MC dataset "weighted" 
        wmcds = mcds.makeWeighted ( weight_name )
        mcds.clear () ; del mcds
        mcds  = wmcds 
        ##
        if 1 == iter % 10  :
            title = '%s MC dataset' % tag 
            logger.info ( '%s:\n%s' % ( title , mcds.table ( title = title , prefix = '# ' ) ) ) 

    # =========================================================================
    ## 2) update weights

    pmax  = 1.0 if iter <= 4 else 0.66 + 0.34 * math.exp ( - ( iter - 4.0 ) / 4.0 )
    power = lambda n : pmax if n <= 1 else pmax * ( 0.66666 / n + 0.33333 )

    if   iter < 3 : the_plots = plots [ : 1 ]
    elif iter < 4 : the_plots = plots [ : 4 ]
    else          : the_plots = plots

    ## the_plots = plots
    
    ## weight truncation: avoid very large change of weights for  single iteration 
    wtruncate = ( 0.05 , 3 ) if iter < 4 else ( 0.5 , 1.3 )  
    
    with timing ( tag + ': make actual reweighting:' , logger = logger ) :
        
        # =========================================================================
        ## 2a) the most important line: perform single iteration step  
        active , _ = makeWeights (
            mcds                   , ## what to be reweighted
            the_plots              , ## reweighting plots/setup
            dbname                 , ## DBASE with reweighting constant 
            delta      = 0.05      , ## stopping criteria
            minmax     = 0.10      , ## stopping criteria
            maxchi2    = 0.50      , ## stopping criteria             
            power      = power     , ## tune: effective power
            wtruncate  = wtruncate , ## truncate weights 
            make_plots = True      , ## make control plots 
            tag        = tag       ) ## tag for printout

    # =============================================================================
    if 1 == iter % 3  or True : 
        
        with timing ( tag + ': compare DATA & weighted MC-dataset:' , logger = logger ) :

            ## 3) compare control and signal samples        
            datatree    = ROOT.TChain ( tag_data , files = testdata ) 
            comparisons = data_compare ( comparators  ,
                                         datatree     ,
                                         mcds         ,
                                         expressions  =  ( 'x' , 'y' , 'z' ) ) 
            ## effective MC statistics at this step 
            n_eff   = mcds.nEff () 
            row     = ( '%d'  % iter , '%.1f' % n_eff )
            pvalues = tuple ( VE ( r.pvalue ) * 100 for r in comparisons  )
            row     += tuple ( '%.2f%s%.2f' % ( pv.value () , plus_minus , pv.error () ) for pv in pvalues)
            glob_stat.append ( row )
            
            title = 'Global DATA/MC similarity p-values [%]'
            table = T.table ( glob_stat , title = title , prefix = '# ' , alignment = alignment )
            logger.info ( '%s\n%s' % ( title , table ) )
                   
    ## explicitely clear and delete dataset
    mcds.clear()
    del mcds
    
    if not active and 5 <= iter : 
        logger.info    ( allright ( 'No more iterations, converged after #%d' % iter ) )
        converged = True 
        break
    
else :
    
    del mcds_init 
    converged = False 
    logger.error ( "No convergency!" )

# ===========================================================================
logger.attention ( 'Memory:%+.2f[MB]' % ( memory_usage () - memory_init ) )                            
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

# =============================================================================
## Add obtained weight into MC-tree
# =============================================================================
with timing ( "Add `%s' column to initial MC-tree" % weight_name , logger = logger ) : 
    mctree   = ROOT.TChain ( tag_mc  , files = testdata )  
    weighter = Weight ( dbname , weightings )
    mctree   = mctree.add_reweighting ( weighter ,  name = weight_name )
    mctree   = ROOT.TChain ( tag_mc  , files = testdata )  
    
# =============================================================================
## For comparison try to use "other" reweighters 
# =============================================================================

datatree   = ROOT.TChain ( tag_data , files = testdata )  
mctree     = ROOT.TChain ( tag_mc   , files = testdata )

from ostap.tools.data_reweighter import DataReweighter

weights = [ '' , weight_name ]

# =============================================================================
## (1) GBReweighter by Alex Rogozhnikov from hep_ml 
# =============================================================================
from ostap.tools.reweighter      import Reweighter     as GBRW 
rw1 = DataReweighter ( GBRW                        , ## reweighter type 
                       original         = mctree   ,
                       target           = datatree ,
                       target_variables = 'x,y,z'  ) 

weight_GBRW = 'weight_GBRW'
with timing ( "GBRW reweight" , logger = logger ) :    
    rw1.reweight ( mctree , name = weight_GBRW ) 
    weights.append ( weight_GBRW )
    
# =============================================================================
## (2) home-made reweighter based on LightGBM
# =============================================================================
from ostap.tools.reweighters     import Reweighter_LGBM     as LGBM
rw2 = DataReweighter ( LGBM                        , ## reweighter type 
                       original         = mctree   ,
                       target           = datatree ,
                       target_variables = 'x,y,z'  ) 

weight_LGBM = 'weight_LGBM'
with timing ( "LGBM reweight" , logger = logger ) :    
    rw2.reweight ( mctree , name = weight_LGBM ) 
    weights.append ( weight_LGBM )

# =============================================================================
## (3) home-made reweighter based on XGBoost  
# =============================================================================
from ostap.tools.reweighters     import Reweighter_XGB     as XGB
rw3 = DataReweighter ( XGB                         , ## reweighter type 
                       original         = mctree   ,
                       target           = datatree ,
                       target_variables = 'x,y,z'  ) 

weight_XGB = 'weight_XGB'
with timing ( "XGB reweight" , logger = logger ) :    
    rw3.reweight ( mctree , name = weight_XGB  ) 
    weights.append ( weight_XGB )

# =============================================================================
## (4) home-made reweighter based on CatBoost  
# =============================================================================
from ostap.core.cpu_info import HAS_AVX2 
if HAS_AVX2 :
    from ostap.tools.reweighters     import Reweighter_CAT     as CAT
    rw3 = DataReweighter ( CAT                         , ## reweighter type 
                           original         = mctree   ,
                           target           = datatree ,
                           target_variables = 'x,y,z'  ) 
    
    weight_XGB = 'weight_CAT'
    with timing ( "XGB reweight" , logger = logger ) :    
        rw3.reweight ( mctree , name = weight_CAT ) 
        weights.append ( weight_CAT )

# ============================================================================
## Compare the quality of all reweighters 
# ============================================================================
glob_stat = [ header ]
for weight in weights :
    
    mctree      = ROOT.TChain ( tag_mc   , files = testdata )
    title       = "MC-tree after %s reweighting" %  ( weight if weight else "NO" ) 
    logger.info ( '%s:\n%s' % ( title , mctree.table2   ( variables = [ 'x' , 'y' , 'z' ] ,
                                                          title     = title    ,
                                                          cuts      = weight   , 
                                                          prefix    = '# '     ) ) )

    datatree    = ROOT.TChain  ( tag_data    , files = testdata ) 
    mctree      = ROOT.TChain  ( tag_mc      , files = testdata )
    comparisons = data_compare ( comparators ,
                                 datatree    ,
                                 mctree      ,
                                 expressions = ( 'x' , 'y' , 'z' ) ,
                                 cuts2       = weight ) 
    ## effective MC statistics at this step 
    n_eff   = mctree.nEff ( weight ) 
    row     = ( weight , '%.1f'  % n_eff )
    pvalues = tuple ( VE ( r.pvalue ) * 100 for r in comparisons  )
    row    += tuple ( '%.2f%s%.2f' % ( pv.value () , plus_minus , pv.error () ) for pv in pvalues)
    glob_stat.append ( row )
    
    title = 'Global DATA/MC similarity: p-values [%]'
    table = T.table ( glob_stat , title = title , prefix = '# ' , alignment = alignment )
    logger.info ( '%s\n%s' % ( title , table ) )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
