#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/tests/test_tools_reweight5.py
#  Test for nD-reweighting machinery
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2023-01-20
# =============================================================================
""" Test for nD-reweighting machinery
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2023-01-20"
__all__     = ()  ## nothing to be imported 
# =============================================================================
## ATTENTION! 
import os 
os.environ [ "MKL_NUM_THREADS"      ]  = "2"
os.environ [ "OMP_NUM_THREADS"      ]  = "2"
os.environ [ "OPENBLAS_NUM_THREADS" ]  = "2"
# =============================================================================
os.environ [ 'KMP_DUPLICATE_LIB_OK' ] = 'TRUE'
os.environ [ 'KMP_INIT_AT_FORK'     ] = 'FALSE' 
os.environ [ 'KMP_BLOCKTIME'        ] = '0'
# =============================================================================
os.environ [ 'OMP_PROC_BIND'        ] = 'false'
os.environ [ 'OMP_WAIT_POLICY'      ] = 'PASSIVE' 
# =============================================================================
os.environ [ 'MKL_THREADING_LAYER'  ] = 'GNU' 
# =============================================================================
from   ostap.utils.root_utils   import batch_env 
from   ostap.utils.core         import typename 
from   ostap.utils.basic        import numcpu, is_ci 
from   ostap.logger.pretty      import nice_print, pretty_floats, format_pow10 
from   ostap.stats.utils        import nEff, num_samples, w_header    
from   ostap.stats.tools        import ( hasLightGBM , hasXGBoost ,
                                         hasCatBoost , hasSkLearn ,
                                         hasPyTorch  , hasHepML   )
from   ostap.logger.symbols     import p_value_percent_expr 
from   ostap.tools.reweighter   import CascadeReweighter as  CASCADE
from   ostap.utils.utils4tests  import light_tests,   active_test
import ostap.io.zipshelve       as     DBASE
import ostap.logger.table       as     T 
import ostap.io.root_file 
import ostap.parallel.kisa
import os, numpy, math , random  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
logger = getLogger ( 'ostap.test_tools_reweight5' )
# =============================================================================    
logger.info ( 'Test for nD-Reweighting machinery')
# ============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
light, weak = light_tests ( 8 )

random_freq = 0.10
 
# =============================================================================
if   light : logger.warning ( "CONTINUOUS_INTEGRATION: run *LIGHT* version of tests" ) 
elif weak  : logger.warning ( "Machine is not powerful enough: reduce some tests" )

# =============================================================================
import numpy
import numpy

def make_datasets ( n_target    = 5000 ,
                    n_original1 = None ,
                    n_original2 = None ,
                    n_dim       = 3    ,
                    seed        = None ,
                    low         = -2.0 , 
                    high        = +2.0 ,
                    scale_targ  =  0.5 ,
                    scale_orig  =  1.0 ) :
    """
    Generates target, original1, and original2 datasets strictly bounded within 
    a hypercube [low, high]^d using smoothly decaying, truncated Gaussians.
    
    - Originals: Shared wide truncated Gaussian pool split into original1 & original2.
    - Target:    Narrower truncated Gaussian with local event shifts and 
                 non-linear quadratic correlations.
    """
    rng = numpy.random.default_rng ( seed )
    
    if scale_orig <= scale_targ :
        raise ValueError ( "scale_orig must be greater than scale_targ to keep weights stable." )

    if isinstance ( n_target, int ) and 1000 <= n_target : pass
    else : raise ValueError ( "Invalid n_target %s" % n_target )
        
    if n_original1 is None or n_original1 < n_target : n_original1 = n_target 
    if n_original2 is None or n_original2 < n_target : n_original2 = n_original1 

    def generate_filtered ( n_needed , scale , is_target = False , rr = None ) :
        
        collected = []
        factor    = 2.0
        
        while sum ( len ( chunk ) for chunk in collected ) < n_needed :
            remaining = n_needed - sum ( len ( chunk ) for chunk in collected )
            n_gen     = int ( remaining * factor )
            
            if not is_target : data = rng.normal ( loc = 0.0, scale = scale, size = ( n_gen, n_dim ) )
            else :
                event_locs = rng.normal ( loc = 0.0        , scale = 0.10  , size = ( n_gen, n_dim ) )
                data       = rng.normal ( loc = event_locs , scale = scale , size = ( n_gen, n_dim ) )
                
                for n in range ( 1, n_dim ) :
                    data [ : , n ] = data [ : , n ] + rr [ n ] * ( data [ : , n - 1 ] ** 2 )
            
            mask       = numpy.all ( ( data >= low ) & ( data <= high ) , axis = 1 )
            valid_data = data [ mask ]
            
            if 0 < len ( valid_data ) : collected.append ( valid_data )
                
            factor *= 1.5
            
        result = numpy.vstack ( collected )
        return result [ : n_needed ]

    # 1. Generate a single combined pool for both originals to ensure they share 
    # the exact same distribution pool, then split them into original1 and original2.
    total_originals_needed = n_original1 + n_original2
    combined_originals     = generate_filtered ( total_originals_needed , scale = scale_orig , is_target = False )
    
    original1 = combined_originals [ : n_original1 ]
    original2 = combined_originals [ n_original1 : ]

    # 2. Target dataset
    rr     = rng.uniform ( low = -0.5 , high = +0.5 , size = n_dim )
    target = generate_filtered ( n_target , scale = scale_targ, is_target = True, rr = rr )

    return target, original1, original2

# ==============================================================================
## Compare datasets using several methods 
# ==============================================================================
from ostap.stats.gof_np       import  ( KullbackLeibler as COMPARATORA , 
                                        Jeffrey         as COMPARATORB , 
                                        JensenShannon   as COMPARATORC ,
                                        Mahalanobis     as COMPARATORD ,
                                        Hotelling       as COMPARATORE , 
                                        Bhattacharyya   as COMPARATORF , 
                                        Wasserstein     as COMPARATORG , 
                                        Hellinger       as COMPARATORH ) 

from ostap.stats.data_compare import data_compare
comparators = []

if active_test ( random_freq ) :  comparators .append ( COMPARATORA ( parallel = True , nToys = 100 ) )
if active_test ( random_freq ) :  comparators .append ( COMPARATORB ( parallel = True , nToys = 100 ) ) 
if active_test ( random_freq ) :  comparators .append ( COMPARATORC ( parallel = True , nToys = 100 ) )
if active_test ( random_freq ) :  comparators .append ( COMPARATORD ( parallel = True , nToys = 100 ) )
if active_test ( random_freq ) :  comparators .append ( COMPARATORE ( parallel = True , nToys = 100 ) )
if active_test ( random_freq ) :  comparators .append ( COMPARATORF ( parallel = True , nToys = 100 ) )
if active_test ( random_freq ) :  comparators .append ( COMPARATORG ( parallel = True , nToys = 100 ) )
if active_test ( random_freq ) :  comparators .append ( COMPARATORH ( parallel = True , nToys = 100 ) )

if not comparators : comparators .append ( COMPARATORA ( parallel = True , nToys = 100 ) )

if active_test ( random_freq ) and hasCatBoost  () :
    
    from ostap.stats.adval        import ADVAL_CATB  as COMPARATOR7
    comparators .append ( COMPARATOR7  ( parallel = True , nToys = 30 ) ) 

if active_test ( random_freq ) and hasLightGBM () :
    
    from ostap.stats.adval        import ADVAL_LGBM  as COMPARATOR5
    comparators .append ( COMPARATOR5  ( parallel = True , nToys = 30 ) ) 

if active_test ( random_freq ) and hasXGBoost  () :
    
    from ostap.stats.adval        import ADVAL_XGB  as COMPARATOR6
    comparators .append ( COMPARATOR6  ( parallel = False , nToys = 30 , n_jobs = 1 ) ) 

if active_test  ( random_freq ) and hasSkLearn ():
    
    from ostap.stats.adval        import ADVAL_HGBC  as COMPARATOR8
    comparators .append ( COMPARATOR8  ( parallel = True , nToys = 30 ) )
    
if False and active_test ( random_freq ) and hasSkLearn ():

    from ostap.stats.adval        import ADVAL_GBC   as COMPARATOR9
    comparators .appemd  ( COMPARATOR9  ( parallel = True , nToys = 25 ) ) 

if False and active_test ( random_freq ) and hasSkLearn ():

    from ostap.stats.adval        import ADVAL_RF   as COMPARATOR10
    comparators .append ( COMPARATOR10 ( parallel = True , nToys = 25 ) ) 

comparators = tuple ( comparators )

# ============================
def run_reweight ( n_dim       = 3     ,
                   n_target    = 10000 , 
                   n_original1 = None  , 
                   n_original2 = None  , 
                   seed        = None  ) :
    
    ## generate samples
    target, original1 , original2 = None , None, None 
    dbname = os.path.basename ( __file__ ) .replace ( '.py' , '.db' )
    if os.path.exists ( dbname ) : 
        with DBASE.open ( dbname , 'READ' ) as db :
            db.ls () 
            target    = db.get ( 'target'    , None )
            original1 = db.get ( 'original1' , None )
            original2 = db.get ( 'original2' , None )
            seed      = db.get ( 'seed'      , seed )
            if target is not None and original1 is not None and original2 is not None :
                logger.attention ( 'Read datasets from %s' % dbname )
                
    if target is None or original1 is None or original2 is None : 
        target , original1 , original2  = make_datasets ( n_target    = n_target    ,
                                                          n_original1 = n_original1 ,
                                                          n_original2 = n_original2 ,                                                          
                                                          n_dim       = n_dim       ,
                                                          seed        = seed        ) 
        logger.attention ( 'New datasets are regenerated and saved into %s' % dbname )
        with DBASE.open ( dbname ) as db :
            db [ 'target'    ] = target
            db [ 'original1' ] = original1
            db [ 'original2' ] = original2
            db [ 'seed'      ] = seed
            db.ls()

    ## (LARGE) combined original 
    original3  = numpy.vstack ( ( original1 , original2 ) )

    
    reweighters = [ None ]
    reweighters = [] 
    
    if active_test ( random_freq ) and hasCatBoost () :
        
        from ostap.tools.reweighters  import CatBoostDensityReweighter as CATB
        rw = CATB ( target   = target      , 
                    original = original1   )                           
        reweighters.append ( rw )
        
        
    if active_test ( random_freq ) and hasCatBoost () :
        
        from ostap.tools.reweighters  import CatBoostIterativeReweighter as CATB
        rw = CATB ( target   = target    , 
                    original = original1 )        
        reweighters.append ( rw )
        
        
    if active_test ( random_freq ) and hasLightGBM () :
        
        from ostap.tools.reweighters  import LightGBMDensityReweighter as  LGBM
        rw = LGBM ( target                = target    ,
                    original              = original1 )
        
        reweighters.append ( rw )

    if active_test ( random_freq ) and hasLightGBM () :
        
        from ostap.tools.reweighters  import LightGBMIterativeReweighter as  LGBM
        rw = LGBM ( target                = target    ,
                    original              = original1 , 
                    progress              = False    )
        
        reweighters.append ( rw )
        

    if active_test ( random_freq ) and hasXGBoost  () :
        
        from ostap.tools.reweighters  import XGBoostDensityReweighter as  XGB 
        rw = XGB  ( target                = target    ,
                    original              = original1 ,
                   )
        
        reweighters.append ( rw )
        
    if active_test ( random_freq ) and hasXGBoost  () :

        from ostap.tools.reweighters  import XGBoostIterativeReweighter as  XGB 
        rw = XGB  ( target                = target    ,
                     original              = original1 ,
                     progress              = False     ,
                    )
        
        reweighters.append ( rw )
        
    if active_test ( random_freq ) and hasPyTorch () :

        from ostap.tools.reweighters  import PyTorchDensityReweighter as TORCH
        rw = TORCH ( target        = target        ,
                     original      = original1     ,
                     hidden_layers = ( 256 , 256 ) )                       
        reweighters.append ( rw )

        
    if active_test ( random_freq ) and hasPyTorch () :

        from ostap.tools.reweighters  import PyTorchIterativeReweighter as TORCH
        rw = TORCH ( target        = target        ,
                     original      = original1     ,
                     ## use very simple construction:
                     hidden_layers           = ( 16 , )  ,
                     progress                = False     ,
                     )
        
        reweighters.append ( rw )
                      

    if active_test ( random_freq ) and hasHepML    () :
        from ostap.tools.reweighters  import GBReweighter   as GBRW
        rw = GBRW ( target              = target         ,
                    original            = original1      ,
                   ) 
        
        reweighters.append ( rw )
        
    if active_test ( random_freq ) and hasSkLearn ():
        
        from ostap.tools.reweighters  import LogRegressionDensityReweighter as LOGREG 
        
        rw = LOGREG ( target       = target    ,                       
                      original     = original1 , 
                      polynomials2 = 2         )
        
        reweighters.append ( rw )
        
        rw = CASCADE ( stages   = ( ( LOGREG , { 'polynomials' : 4 } ) ,
                                    ( LOGREG , { 'polynomials' : 4 } ) ) , 
                       target   = target ,
                       original = original1 )
        reweighters.append ( rw )

        rw = CASCADE ( stages   = ( ( LOGREG , { 'polynomials' : 4 } ) ,
                                    ( LOGREG , { 'polynomials' : 4 } ) , 
                                    ( LOGREG , { 'polynomials' : 4 } ) ) , 
                       target   = target ,
                       original = original1 )
        reweighters.append ( rw )

    # =========================================================================
    ## Summary table 
    # =========================================================================
    
    header    = [] 
    rows      = []  

    orginals = ( original2 , ) if ( light or weak ) else ( original2 , ) 
    
    for rw in reweighters :
            
        for original in originals : 
        
            if rw : logger.info ( 'Using REWEIGHTER %s %s ' % ( typename ( rw ) , rw.method ) )
            
            if rw is None : original_weight = None
            else          : original_weight = rw.weights ( original )
            
            wmin, wmax , wunit = '' , '' , '' 
            if not original_weight is None :
                wmin = float ( numpy.min ( original_weight ) ) 
                wmax = float ( numpy.max ( original_weight ) )
                ww   , expo = pretty_floats ( wmin , wmax , precision = 4 , with_sign = False , strip = True )
                wmin , wmax = ww
                wunit       =  ( '[%s]' % format_pow10 ( expo ) ) if expo else '' 
                               
            for c in comparators :
                
                tv , pv = c.pvalue ( data1   = target          ,
                                     data2   = original        , 
                                     weight2 = original_weight ) 
                
                header , row = c.the_row ()
                
                neff    = nEff ( original , original_weight ) 
                rw_type = '' if rw is None else rw.method
                rw_symb = '' if rw is None else rw.symbol                 
                neff    = nice_print ( neff , precision = 2 , with_sign = False , strip = True ) 
                row     = ( rw_symb , rw_type , c.symbol, c.method , neff ) + tuple ( row ) + ( wmin , wmax , wunit )             
                rows.append ( row )

            header =  ( '' , 'Reweighter' , '' , 'Comparator'  , '#eff' ) + header + w_header 
            if rows [ 0 ] != header :  rows         = [ header ] + rows
            title = 'Test results : %s' % p_value_percent_expr
            table = T.remove_empty_columns ( rows ) 
            table = T.table ( table , title = title , prefix = '# ' , alignment = 'clclccccccccccccccccccccccccc' )
            logger.info ( '%s\n%s' % ( title , table ) )

    # =================================================================
    # can we save reweighters into database ?
    with DBASE.tmpdb () as db :
        
        for rw in reweighters : db [  typename ( rw ) ] = rw                
        db [ 'reweighters' ] = reweighters
        
        db.ls () 
        
# =============================================================================
if '__main__' == __name__ :

    
    run_reweight ( 3 , 3000 , 7000 , 8000 ) 
    
# =============================================================================
        
# =============================================================================
##                                                                      The END 
# =============================================================================
