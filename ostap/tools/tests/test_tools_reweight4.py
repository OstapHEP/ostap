#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/tests/test_tools_reweight4.py
#  Test for 3D-reweighting machinery
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2023-01-20
# =============================================================================
"""Test for 3D-reweighting machinery
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2025-09-15"
__all__     = ()  ## nothing to be imported 
# =============================================================================
from   ostap.core.pyrouts     import Ostap, VE 
from   ostap.histos.histos    import h1_axis, h2_axes, h3_axes 
from   ostap.utils.timing     import timing
from   ostap.utils.basic      import numcpu 
from   ostap.logger.colorized import attention, allright
from   ostap.plotting.canvas  import use_canvas
from   ostap.utils.root_utils import batch_env 
from   ostap.utils.cleanup    import CleanUp
from   ostap.utils.memory     import memory_usage, delta_ram
from   ostap.logger.symbols   import ( iteration      , sup_eff  ,
                                       plus_minus     , script_p ,
                                       triangular_flag as start_symbol , 
                                       trophy          as final_symbol ) 
from   ostap.stats.tools      import ( hasLightGBM , hasXGBoost ,
                                       hasCatBoost , hasSkLearn ,
                                       hasHepML    ) 
import ostap.io.zipshelve     as     DBASE
import ostap.logger.table     as     T 
import ostap.logger.table     as     T
import ostap.histos.histos  
import ostap.io.root_file 
import ostap.parallel.kisa
import ROOT, random, math, os, time 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_tools_reweight4' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================    
logger.info ( 'Test for ND-Reweighting machinery')
# ============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

if 8 <= numcpu () : 
    
    NDATA  =  500000
    NMC    = 1000000

    NDATA  =  5000
    NMC    = 10000
    
else :
     
    NDATA  =   1200
    NMC    =   1200

     
NDATA  =   2000
NMC    =   2000
    
tag_data_r1  = 'DATA_R1_histogram'
tag_data_r2  = 'DATA_R2_histogram'
tag_data_r3  = 'DATA_R3_histogram'
tag_data_r   = 'DATA_R_histogram'
tag_data_3d  = 'DATA_3D_histogram'
tag_data_2d  = 'DATA_2D_histogram'
tag_data_12  = 'DATA_12_histogram'
tag_data_23  = 'DATA_23_histogram'
tag_data_31  = 'DATA_31_histogram'

tag_data     = 'DATA_tree'
tag_mc       = 'MC_tree'

dbname       = CleanUp.tempfile ( suffix = '.db' , prefix ='ostap-test-tools-reweight4-'   )

rmax = 20 

ir              = 50
ix2 , iy2 , iz2 = 8 , 8 , 8 
ix3 , iy3 , iz3 = 3 , 3 , 3 
def prepare_data ( nF = 5 ) :
    ## 
    assert 1 <= nF , "Invalid number of files is specified!"
        
    seed =  1234567890 
    random.seed ( seed ) 
    logger.info ( 'Test *RANDOM* data will be generated/seed=%s' % seed  )   
    ## prepare "data" histograms:

    ## X,Y,Z histograms

    hr1_histo  = h1_axis ( [ rmax/ir*i  for i in range ( ir + 1 ) ] )
    hr2_histo  = h1_axis ( [ rmax/ir*i  for i in range ( ir + 1 ) ] )
    hr3_histo  = h1_axis ( [ rmax/ir*i  for i in range ( ir + 1 ) ] )
    hr_histo   = h1_axis ( [ rmax/ir*i  for i in range ( ir + 1 ) ] )

    h12_histo  = h2_axes ( [ rmax/ix2*i for i in range ( ix2 + 1 ) ] ,
                           [ rmax/iy2*i for i in range ( iy2 + 1 ) ] )
    
    h23_histo  = h2_axes ( [ rmax/iy2*i for i in range ( iz2 + 1 ) ] ,
                           [ rmax/iz2*i for i in range ( iz2 + 1 ) ] )
    
    h31_histo  = h2_axes ( [ rmax/iz2*i for i in range ( iz2 + 1 ) ] ,
                           [ rmax/ix2*i for i in range ( ix2 + 1 ) ] )
    
    h3d_histo  = h3_axes ( [ rmax/ix3*i for i in range ( ix3 + 1 ) ] ,
                           [ rmax/iy3*i for i in range ( iy3 + 1 ) ] ,
                           [ rmax/iz3*i for i in range ( iz3 + 1 ) ] )

    

    ND1 = 1 + NDATA // nF 
    NM1 = 1 + NMC   // nF

    files = []
    
    for f in range ( nF ) :
        
        testfile   = CleanUp.tempfile ( suffix = '.root' , prefix ='ostap-test-tools-reweight4-' )        
        with ROOT.TFile.Open ( testfile ,'recreate') as mc_file :
            
            mc_file.cd() 
            
            datatree  = ROOT.TTree ( tag_data , 'data-tree' )
            datatree .SetDirectory ( mc_file ) 
            
            from array import array 
            r1var = array    ( 'f', [0])
            r2var = array    ( 'f', [0])
            r3var = array    ( 'f', [0])
            
            datatree.Branch ( 'r1' , r1var , 'r1/F' )
            datatree.Branch ( 'r2' , r2var , 'r2/F' )
            datatree.Branch ( 'r3' , r3var , 'r3/F' )
            
            ## Gaussian 3D-component with correlations 
            
            nentries = 0  
            while nentries <  ND1 :
                
                r1 = random.uniform ( 0 , rmax ) 
                r2 = random.uniform ( 0 , rmax ) 
                r3 = random.uniform ( 0 , rmax ) 

                rr = r1 , r2 , r3 
                if 0.8 <= random.uniform ( 0 , 1 ) :
                    r1 , r2 , r3 = sorted ( rr ) 
                                    
                hr1_histo.Fill ( r1 ) 
                hr2_histo.Fill ( r2 ) 
                hr3_histo.Fill ( r3 ) 

                hr_histo.Fill  ( r1 ) 
                hr_histo.Fill  ( r2 ) 
                hr_histo.Fill  ( r3 ) 

                h12_histo.Fill ( r1 , r2 ) 
                h23_histo.Fill ( r2 , r3 ) 
                h31_histo.Fill ( r3 , r1 ) 

                h3d_histo.Fill ( r1 , r2 , r3 )
                
                r1var [ 0 ] = r1 
                r2var [ 0 ] = r2 
                r3var [ 0 ] = r3 
                
                datatree.Fill()
                nentries += 1 
            
            datatree.Write()

            ## write the histogram 
            
            mc_file [ tag_data_r1 ] = hr1_histo
            mc_file [ tag_data_r2 ] = hr2_histo
            mc_file [ tag_data_r3 ] = hr3_histo
            mc_file [ tag_data_3d ] = h3d_histo
            
            mc_file [ tag_data_12 ] = h12_histo
            mc_file [ tag_data_23 ] = h23_histo
            mc_file [ tag_data_31 ] = h31_histo
            
            mc_file [ tag_data_r  ] = hr_histo
            
            ## 
            mctree  = ROOT.TTree ( tag_mc , 'mc-tree' )
            mctree .SetDirectory ( mc_file ) 
            
            from array import array 

            r1var = array  ( 'f', [ 0.0 ] )
            r2var = array  ( 'f', [ 0.0 ] )
            r3var = array  ( 'f', [ 0.0 ] )
            
            mctree.Branch ( 'r1' , r1var , 'r1/F' )
            mctree.Branch ( 'r2' , r2var , 'r2/F' )
            mctree.Branch ( 'r3' , r3var , 'r3/F' )
            
            nentries = 0 
            while nentries <  NM1 :
                
                r1 = random.gauss ( 0.5  * rmax , 0.6 * rmax )
                if not 0 < r1 < rmax : continue
                
                r2 = random.gauss ( 0.5  * rmax , 0.6 * rmax )
                if not 0 < r2 < rmax : continue
                
                r3 = random.gauss ( 0.5  * rmax , 0.6 * rmax )
                if not 0 < r3 < rmax : continue
                
                rrr = r1 , r2  , r3 
                
                ## if 0.50 <= random.uniform ( 0 , 1 ) : rrr = sorted ( rrr )
                
                r1var [ 0 ] = rrr [ 0 ] 
                r2var [ 0 ] = rrr [ 1 ]
                r3var [ 0 ] = rrr [ 2 ]
                
                mctree.Fill()                
                nentries += 1 
                
            mctree .Write()
            mc_file.Write() 

        files.append ( testfile )

    return files


## list of test-files 
testdata = prepare_data ( nF = 5 )
    
# =============================================================================
## Read data from DB
# =============================================================================
hr1_data  = None 
hr2_data  = None 
hr3_data  = None 
h3d_data  = None
h12_data  = None
h23_data  = None
h31_data  = None
hr_data   = None

for test_file in testdata :
    
    with ROOT.TFile.open ( test_file , 'r' ) as dbroot :

        hr1 = dbroot [ tag_data_r1 ]
        hr2 = dbroot [ tag_data_r2 ]
        hr3 = dbroot [ tag_data_r3 ]
        h3d = dbroot [ tag_data_3d ]
        h12 = dbroot [ tag_data_12 ]
        h23 = dbroot [ tag_data_23 ]
        h31 = dbroot [ tag_data_31 ]
        hr  = dbroot [ tag_data_r ]
        
        if not hr1_data : hr1_data  = hr1.clone()
        else            : hr1_data += hr1

        if not hr2_data : hr2_data  = hr2.clone()
        else            : hr2_data += hr2
        
        if not hr3_data : hr3_data  = hr3.clone()
        else            : hr3_data += hr3

        if not hr_data  : hr_data   = hr.clone()
        else            : hr_data  += hr
        
        h3d_data  = h3d.clone()
        h12_data  = h12.clone()
        h23_data  = h23.clone()
        h31_data  = h31.clone()
        
# =============================================================================
## prebook MC histograms
# =============================================================================

ir     = 32
mc_r1  = h1_axis ( [ rmax/ir*i  for i in range ( ir  + 1 ) ] )
mc_r2  = h1_axis ( [ rmax/ir*i  for i in range ( ir  + 1 ) ] )
mc_r3  = h1_axis ( [ rmax/ir*i  for i in range ( ir  + 1 ) ] )
mc_r   = h1_axis ( [ rmax/ir*i  for i in range ( ir  + 1 ) ] )
mc_3d  = h3_axes ( [ rmax/ix3*i for i in range ( ix3 + 1 ) ] ,
                   [ rmax/iy3*i for i in range ( iy3 + 1 ) ] ,
                   [ rmax/iz3*i for i in range ( iz3 + 1 ) ] )

mc_12 = h12_data.clone()
mc_23 = h23_data.clone()
mc_31 = h31_data.clone()

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
    Weight.Var ( 'r1'  , 'r1-reweight' ) , 
    Weight.Var ( 'r2'  , 'r2-reweight' ) , 
    Weight.Var ( 'r3'  , 'r3-reweight' ) , 
    )

# =============================================================================
## variables to be used in MC-dataset 
variables  = [
    Variable ( 'r1' , 'r1-var'  , 0  , rmax ) , 
    Variable ( 'r2' , 'r2-var'  , 0  , rmax ) ,
    Variable ( 'r3' , 'r3-var'  , 0  , rmax ) ,
    ]


vars = 'r1' , 'r2' , 'r3'

# =============================================================================
datatree   = ROOT.TChain ( tag_data , files = testdata )  
title      = 'Data/target dataset %d' % len ( datatree ) 
logger.info ( '%s:\n%s' % ( title , datatree.table ( variables = vars  ,
                                                     title     = title ,
                                                     prefix    = '# '  ) ) )

# ==============================================================================
## Compare datasets using several methods 
# ==============================================================================
import ostap.stats.gof_np as GnP
cconf       = { 'parallel' : True , 'nToys' : 20 , 'silent' : True , 'progress' : True } 
comparators = (
    GnP.Chi2            ( **cconf ) ,
    GnP.KullbackLeibler ( **cconf ) ,
    GnP.Jeffrey         ( **cconf ) ,
    GnP.JensenShannon   ( **cconf ) ,
    GnP.Mahalanobis     ( **cconf ) ,
    GnP.Hotelling       ( **cconf ) ,
    GnP.Bhattacharyya   ( **cconf ) ,
    GnP.Wasserstein     ( **cconf ) ,
    GnP.Hellinger       ( **cconf ) )

# =========================================================================
has_lightgbm  = hasLightGBM  ()
if has_lightgbm :  logger.attention ( 'USE LightGBM!'              )
else            :  logger.warning   ( 'LightGBM is not available!' )
            
has_xgboost   = hasXGBoost  ()
if has_xgboost  :  logger.attention ( 'USE XGBoost!'               )
else            :  logger.warning   ( 'XGBoost is not available!'  )

has_catboost  = hasCatBoost  ()
if has_catboost :  logger.attention ( 'USE CatBoost!'              )
else            :  logger.warning   ( 'CatBoost is not available!' )

has_sklearn   = hasSkLearn  ()
if has_sklearn  :  logger.attention ( 'USE SkLearn!'              )
else            :  logger.warning   ( 'SkLearn  is not available!' )

if has_lightgbm :  
    from ostap.stats.adval        import ADVAL_LGBM  as CMP 
    comparators += ( CMP ( **cconf ) , ) 

if has_xgboost:  
    from ostap.stats.adval        import ADVAL_XGB  as CMP
    comparators += ( CMP ( **cconf ) , ) 

if has_catboost:  
    from ostap.stats.adval        import ADVAL_CATB  as CMP 
    comparators += ( CMP ( **cconf ) , ) 

if False and has_sklearn:
    
    from ostap.stats.adval        import ADVAL_HGBC  as CMP1 
    comparators += ( CMP1 ( **cconf ) , )
    
    from ostap.stats.adval        import ADVAL_GBC   as CMP2
    comparators += ( CMP2 ( **cconf ) , ) 

comparators = comparators
# ============================================================================
## The table of global comparison statistics 
header    = ( '#%s' % iteration , '#%s' % sup_eff ) + tuple ( c.method for c in comparators ) 
glob_stat = [ header ]
alignment = 'lc' + 'c' * len ( comparators )           


# =============================================================================
with timing ( 'Prepare initial MC-dataset:' , logger = logger ) :
    mctree    = ROOT.TChain ( tag_mc , files = testdata )
    ## fill dataset from input MC tree
    mcds_init , _ = mctree.make_dataset ( variables = variables ,
                                          selection = '0<=r1 && r1<%s && 0<=r2 && r2<%s && 0<=r3 && r3<%s' % ( rmax , rmax, rmax ) )
    
# =============================================================================
## Configuration of reweighting plots 
# =============================================================================
plots  = [
    WeightingPlot ( 'r1:r2:r3' , address = '3D-reweight'  , data = h3d_data , mc = mc_3d  ) ,     
    WeightingPlot ( 'r1:r2'    , address = 'r12-reweight' , data = h12_data , mc = mc_12  ) ,
    WeightingPlot ( 'r2:r3'    , address = 'r23-reweight' , data = h23_data , mc = mc_23  ) ,
    WeightingPlot ( 'r3:r1'    , address = 'r31-reweight' , data = h31_data , mc = mc_31  ) ,    
    WeightingPlot ( 'r1'       , address = 'r1-reweight'  , data = hr1_data , mc = mc_r1  ) , ## ignore = True ) ,  
    WeightingPlot ( 'r2'       , address = 'r2-reweight'  , data = hr2_data , mc = mc_r2  ) , ## ignore = True ) ,  
    WeightingPlot ( 'r3'       , address = 'r3-reweight'  , data = hr3_data , mc = mc_r3  ) , ## ignore = True ) ,  
    ## WeightingPlot ( 'r1,r2,r3' , address = 'r-reweight'  , data = hr_data  , mc = mc_r   ) ,
]


from ostap.stats.data_compare import data_compare     
# =============================================================================
## name of weight variable 
weight_name = 'weight'
# =============================================================================
with timing ( 'Compare DATA & initial MC-dataset:' , logger = logger ) :
    
    datatree = ROOT.TChain  ( tag_data , files = testdata )            
    mctree   = ROOT.TChain  ( tag_mc   , files = testdata )            
    results  = data_compare ( comparators ,
                              datatree    ,
                              mctree      ,
                              expressions = vars  , 
                              importance  = False , 
                              silent      = True  )
    
    n_eff = mctree.nEff () 
    trow  = [ start_symbol , '%.1f' % n_eff ] 
    for r in results :
        pv100 = VE ( r.pvalue ) * 100            
        trow .append ( '%6.2f%s%.2f' % ( pv100.value () , plus_minus , pv100.error () ) )
    glob_stat.append ( tuple ( trow ) )
    
    title = 'Global DATA/MC similarity %s-values [%%] before reweighting' % script_p
    table = T.table ( glob_stat , title = title , prefix = '# ' , alignment = alignment )
    logger.info ( '%s\n%s' % ( title , table ) )
    

converged   = False
maxIter     = 6  
memory_init = memory_usage() 

# =============================================================================
## start reweighting iterations:
for iter in range ( 1 , maxIter + 1 ) :

    tag = 'Reweighting iteration #%d%s' %  ( iter , iteration )
    mem = ''
    if 1 < iter : mem = ' Memory: %s=%+.2f[MB]' % ( delta_ram , memory_usage () - memory_init )
    logger.info ( allright ( tag + mem ) )
    
    # =========================================================================
    ## 0) The weighter object
    weighter = Weight ( dbname , weightings )    
    ## 1a) create new "weighted" mcdataset
    mcds = mcds_init.copy() 
    
    with timing ( tag + ': add weight to MC-dataset' , logger = logger ) :
        # =========================================================================
        ## 1b) add  "weight" variable to dataset 
        mcds  = mcds.add_reweighting ( weighter ,  name = weight_name , progress = True )        
        ## 1c) make MC dataset "weighted" 
        wmcds = mcds.makeWeighted ( weight_name  )
        mcds.clear () ; del mcds
        mcds  = wmcds 
        title = 'Reweighted dataset at iteration #%d' % iter 
        logger.info ( '%s:\n%s' % ( title , mcds.table ( variables = vars  ,
                                                         title     = title ,
                                                         prefix    = '# '  ) ) )        
        
    # =========================================================================
    with timing ( tag + ': compare DATA & weighted MC-dataset:' , logger = logger ) :
        
        datatree = ROOT.TChain  ( tag_data , files = testdata )            
        results  = data_compare ( comparators ,
                                  datatree    ,
                                  mcds        ,
                                  expressions = vars , 
                                  importance  = True , 
                                  silent      = True )
        
        n_eff = mcds.nEff () 
        trow  = [ '%d' % iter , '%.1f' % n_eff ] 
        for r in results :
            pv100 = VE ( r.pvalue ) * 100            
            trow .append ( '%6.2f%s%.2f' % ( pv100.value () , plus_minus , pv100.error () ) )
        glob_stat.append ( tuple ( trow ) )
        
        title = 'Global DATA/MC similarity %s-values [%%]' % script_p
        table = T.table ( glob_stat , title = title , prefix = '# ' , alignment = alignment )
        logger.info ( '%s\n%s' % ( title , table ) )

        
    if   iter < 2 : the_plots = plots [ : 1 ]
    elif iter < 4 : the_plots = plots [ : 4 ]
    else          : the_plots = plots 
        
    with timing ( tag + ': make actual reweighting:' , logger = logger ) :
        
        # =========================================================================
        ## 2a) the most important line: perform single iteration step  
        active , _ = makeWeights (
            mcds                      , ## what to be reweighted
            the_plots                 , ## reweighting plots/setup
            dbname                    , ## DBASE with reweighting constant
            delta      = 0.10         , ## stopping criteria
            minmax     = 0.20         , ## stopping criteria  
            maxchi2    = 1.50         , ## stopping criteria
            wtruncate  = ( 0.05 , 7 ) , ## truncate small/large weights? 
            make_plots = True         , ## make control plots 
            tag        = tag          ) ## tag for printout

    if not active and 5 <= iter : 
        logger.info    ( allright ( 'No more iterations, converged after #%d' % iter ) )
        converged = True 
        break
    
    del mcds
    
else :
    
    converged = False 
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

# =============================================================================
with timing ( "Add weight column to initial MC-tree/chain" , logger = logger ) : 
    mctree   = ROOT.TChain ( tag_mc   , files = testdata )  
    weighter = Weight ( dbname , weightings )
    mctree   = mctree.add_reweighting ( weighter               ,
                                        name     = weight_name ,
                                        report   = True        ,
                                        progress = True        , 
                                        parallel = 5000 < len ( mctree ) )
    mctree   = ROOT.TChain ( tag_mc   , files = testdata )  
    
#  (1) compare
with timing ( 'Compare DATA & final weighted MC-tree:' , logger = logger ) :
    
    datatree = ROOT.TChain  ( tag_data , files = testdata )            
    mctree   = ROOT.TChain  ( tag_mc   , files = testdata )  
    results  = data_compare ( comparators ,
                              datatree    ,
                              mctree      ,
                              expressions = vars  ,
                              cuts2        = weight_name , 
                              importance  = False , 
                              silent      = True  )
    
    n_eff = mctree.nEff () 
    trow  = [ final_symbol , '%.1f' % n_eff ] 
    for r in results :
        pv100 = VE ( r.pvalue ) * 100            
        trow .append ( '%6.2f%s%.2f' % ( pv100.value () , plus_minus , pv100.error () ) )
    glob_stat.append ( tuple ( trow ) )
    
    title = 'Global DATA/MC similarity %s-values [%%]' % script_p
    table = T.table ( glob_stat , title = title , prefix = '# ' , alignment = alignment )
    logger.info ( '%s\n%s' % ( title , table ) )

# ===========================================================================
## from   ostap.tools.reweight import backup_to_ROOT, restore_from_ROOT
##  root_file = backup_to_ROOT    ( dbname     )
## new_db    = restore_from_ROOT ( root_file  )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
