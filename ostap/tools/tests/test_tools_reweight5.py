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
from   ostap.utils.basic        import numcpu
from   ostap.logger.pretty      import nice_print
from   ostap.stats.utils        import nEff, num_samples, p_value_percent, w_header    
from   ostap.stats.tools        import ( hasLightGBM , hasXGBoost ,
                                         hasCatBoost , hasSkLearn ,
                                         hasPyTorch  , hasHepML   )
from   ostap.tools.reweighter   import CascadeReweighter as  CASCADE
import ostap.io.zipshelve       as     DBASE
import ostap.logger.table       as     T 
import ostap.io.root_file 
import ostap.parallel.kisa
import os, numpy, random  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_tools_reweight5' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================    
logger.info ( 'Test for nD-Reweighting machinery')
# ============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
def make_datasets ( n_target    = 5000 ,
                    n_original1 = None ,
                    n_original2 = None ,
                    n_dim       =    3 ,
                    seed        = None ) :
    """ Temporarily simplified datasets for baseline validation:
        smooth distributions with moderate shifts and no sharp boundaries.
    """
    numpy.random.seed ( seed )
    
    if isinstance ( n_target , int ) and 1000 <= n_target : pass
    else : raise ValueError ( "Invalis n_target %s" % n_target )
        
    if n_original1 is None or n_original1 < n_target : n_original1 = n_target 
    if n_original2 is None or n_original2 < n_target : n_original2 = n_original1 

    original1 = numpy.random.normal ( loc = 0.0   , scale = 1.3 , size = ( n_original1 , n_dim ) )    
    original2 = numpy.random.normal ( loc = 0.0   , scale = 1.3 , size = ( n_original2 , n_dim ) )
    
    # 2. Target dataset: slightly shifted normal distribution
    target = numpy.random.normal    ( loc = 0.10  , scale = 0.9 , size = ( n_target , n_dim ) )
    
    # 3.Add some mild & smooth correlations     
    rr = numpy.random.uniform ( low = -0.40 , high = 0.30 , size= n_dim )    
    for n in range ( 1 , n_dim ) :
        target [ : , n ] = target [ : , n ] + rr [ n ]  * ( target [ : , n - 1 ] ** 2 )        

    return target , original1, original2 

# ==============================================================================
## Compare datasets using several methods 
# ==============================================================================
from ostap.stats.gof_np       import  ( KullbackLeibler as COMPARATORA , 
                                        Jeffrey         as COMPARATORB , 
                                        JensenShannon   as COMPARATORC ,
                                        Mahalanobis     as COMPARATORD ,
                                        Hotelling       as COMPARATORE , 
                                        Bhattacharyya   as COMPARATORF , 
                                        Hellinger       as COMPARATORG ) 

from ostap.stats.data_compare import data_compare     
comparators = ( COMPARATORA ( parallel = True , nToys = 100 ) ,
                COMPARATORB ( parallel = True , nToys = 100 ) ,
                COMPARATORC ( parallel = True , nToys = 100 ) , 
                COMPARATORD ( parallel = True , nToys = 100 ) , 
                COMPARATORE ( parallel = True , nToys = 100 ) , 
                COMPARATORF ( parallel = True , nToys = 100 ) , 
                COMPARATORG ( parallel = True , nToys = 100 ) )

if hasCatBoost () :  
    from ostap.stats.adval        import ADVAL_CATB  as COMPARATOR7
    comparators += ( COMPARATOR7  ( parallel = True , nToys = 25 ) , ) 
    
if hasLightGBM () :
    from ostap.stats.adval        import ADVAL_LGBM  as COMPARATOR5
    comparators += ( COMPARATOR5  ( parallel = True , nToys = 25 ) , ) 

if hasXGBoost  () :  
    from ostap.stats.adval        import ADVAL_XGB  as COMPARATOR6
    comparators += ( COMPARATOR6  ( parallel = False , nToys = 25 , n_jobs = 1) , ) 

if hasSkLearn (): 
    from ostap.stats.adval        import ADVAL_HGBC  as COMPARATOR8
    comparators += ( COMPARATOR8  ( parallel = True , nToys = 25 ) , )
    
    from ostap.stats.adval        import ADVAL_GBC   as COMPARATOR9
    comparators += ( COMPARATOR9  ( parallel = True , nToys = 25 ) , ) 

    from ostap.stats.adval        import ADVAL_RF   as COMPARATOR10
    comparators += ( COMPARATOR10 ( parallel = True , nToys = 25 ) , ) 
    
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

    ## combined original 
    original3 = numpy.vstack ( ( original1 , original2 ) )


    reweighters = [ None ]
    reweighters = [] 
    
    if False and hasCatBoost () :

        """
        from ostap.tools.reweighters  import CatBoostDensityReweighter as CATB
        rw4 = CATB ( target               = target      , 
                    original              = original    , 
                    random_state          = seed        , 
                    learning_rate         = 0.02        ,   
                    l2_leaf_reg           = 0.1         ,  
                    early_stopping_rounds = None        , 
                    ## bootstrap_type        = 'Bernoulli' ,
                    ## subsample             = 0.8         , 
                    bootstrap_type        = 'Bayesian'     , 
                    bagging_temperature   = 1.0            , 
                    ##   
                    boosting_type         = 'Ordered'
                    )                           
        reweighters.append ( rw4 )
        """

        from ostap.tools.reweighters  import CatBoostIterativeReweighter as CATB
        rw4 = CATB ( target        = target      , 
                     original      = original    ,
                     ## 
                     leaf_estimation_iterations =    1 , 
                     min_data_in_leaf           =    2 ,
                     learning_rate              = 0.05 ,
                     max_depth                  =    8 , 
                     n_estimators               =  300 ,
                     max_delsta_step            =    1 ,
                     l2_leaf_reg                =  1.0 , 
                     ##
                    )       

        reweighters.append ( rw4 )
        
    if False and hasLightGBM () :
        
        ## from ostap.tools.reweighters  import LightGBMDensityReweighter as  LGBM
        ## rw2 = LGBM ( target                = target   ,
        ##             original              = original ,
        ##             early_stopping_rounds = None     ,
        ##             min_child_samples     = 10       ,
        ##             subsample             = 0.8      ,
        ##             path_smooth           = 0.5      )
        ## 
        ## reweighters.append ( rw2 )

        from ostap.tools.reweighters  import LightGBMIterativeReweighter as  LGBM
        rw2 = LGBM ( target                = target    ,
                     original              = original1 )

        reweighters.append ( rw2 )
        

    if  hasXGBoost  () :
        
        ## from ostap.tools.reweighters  import XGBoostDensityReweighter as  XGB 
        ## rw3 = XGB  ( target                = target   ,
        ##             original              = original ,
        ##             subsample             = 0.8      ,
        ##             early_stopping_rounds = None     , 
        ##             tree_method           = 'exact'  )
        
        ## reweighters.append ( rw3 )

        from ostap.tools.reweighters  import XGBoostIterativeReweighter as  XGB 
        rw3 = XGB  ( target                = target    ,
                     original              = original1 ,
                     ##
                     n_splits              = 5       , 
                     n_estimators          = 300     ,
                     learning_rate         = 0.02    ,
                     ## max_depth             = 7       ,
                     ## 
                     ## tree_method           = 'exact' ,
                     ## max_bin               = 1023    ,
                     ##
                     min_child_weight      = 1e-3    ,
                     max_delta_step        = 9.0     , 
                     subsample             = 1.0     ,
                     ## 
                     gamma                 = 0.0     , 
                     reg_alpha             = 0.0     ,
                     reg_lambda            = 0.0000  , ## 0.0     , ## 1e-3    , ## 0.01 
                     ##
                     grow_policy           = 'depthwise' ,
                     max_depth             =  7          ,
                     ## 
                     ## max_leaves            =  0      ,
                     ##
                     ## grow_policy           = 'lossguide' , 
                     max_leaves            =  63         ,
                     ## tree_method           = 'hist'   ,
                     tree_method           = 'exact'      ,
                     max_bin               = 1023        ,
                     ## 
                     voting                = 'trim:1' ,
                     oof_mode              = False    , 
                     progress              = False    ,
                    )
        
        reweighters.append ( rw3 )
        
    if False and hasPyTorch () :

        """
        from ostap.tools.reweighters  import PyTorchDensityReweighter as TORCH
        rw5 = TORCH ( target        = target        ,
                      original      = original      ,
                      n_splits      = 5             ,
                      epochs        = 1000          ,
                      bce_epochs    = 30            , 
                      min_epochs    = 100           ,
                      eval_freq     = 5             ,  
                      patience      = 20            ,  
                      batch_size    = 2 ** 6        ,
                      dropout       = 0.0           , 
                      learning_rate = 1e-4          , 
                      logit_reg     = 0             ,
                      max_logit     = 10            , 
                      ma_window     = 20            ,
                      weight_decay  = 0             ,
                      min_delta     = 1e-6          ,
                      hidden_dims   = ( 256 , 256 , 256 ) )                       
        reweighters.append ( rw5 )
        """

        from ostap.tools.reweighters  import PyTorchIterativeReweighter as TORCH
        rw5 = TORCH ( target        = target        ,
                      original      = original      )
        
        reweighters.append ( rw5 )
                      

    if False and hasHepML    () :
        from ostap.tools.reweighters  import GBReweighter   as GBRW
        rw1 = GBRW ( target              = target         ,
                     original            = original1      ,
                     n_estimators        = 300            ,
                     learning_rate       = 0.01           , 
                     max_depth           = 7              ,
                     min_samples_leaf    = 5              ,
                     loss_regularization = 1              ) 
        
        reweighters.append ( rw1 )
        
    if False and hasSkLearn ():
        
        from ostap.tools.reweighters  import LogRegressionDensityReweighter as LOGREG 

        ## rw6 = LOGREG ( target = target , original = original )
        ## reweighters.append ( rw6 )
        
        
        ## rw7 = CASCADE ( classes = [ LOGREG , LOGREG ] , target = target , original = original )
        ## reweighters.append ( rw7 )
        
        rw8 = CASCADE ( classes     = [ LOGREG , LOGREG , LOGREG ] ,
                        target      = target   ,
                        original    = original ,
                        ##
                        polynomials = 2 , 
                        ##
                       )
        
        reweighters.append ( rw8 )


    # =========================================================================
    ## Summary table 
    # =========================================================================
    
    header    = [] 
    rows      = []  
    
    for rw in reweighters :
            
        for original in ( original1 , 
                          original2 ,
                          original3 ) :
        
            if rw : logger.info ( 'Using REWEIGHTER %s %s ' % ( typename ( rw ) , rw.method ) )
            
            if rw is None : original_weight = None
            else          : original_weight = rw.weights ( original )
            
            wmin, wmax = '' , '' 
            if not original_weight is None :
                wmin = float ( numpy.min ( original_weight ) ) 
                wmax = float ( numpy.max ( original_weight ) )
                wmin , wmax = '%.5g' % wmin , '%.5g' % wmax 
                               
            for c in comparators :
                
                tv , pv = c.pvalue ( data1   = target          ,
                                     data2   = original        , 
                                     weight2 = original_weight ) 
                
                header , row = c.the_row ()
                
                neff    = nEff ( original , original_weight ) 
                rw_type = '' if rw is None else rw.method            
                neff    = nice_print ( neff , precision = 2 , with_sign = False , strip = True ) 
                row     = ( rw_type , c.method , neff ) + tuple ( row ) + ( wmin , wmax )             
                rows.append ( row )

            header =  ( 'Reweighter' , 'Comparator' , '#eff' ) + header + w_header 
            if rows [ 0 ] != header :  rows         = [ header ] + rows
            title = 'Test results : %s' % p_value_percent 
            table = T.table ( rows , title = title , prefix = '# ' )
            logger.info ( '%s\n%s' % ( title , table ) )


# =============================================================================
if '__main__' == __name__ :

    ## if numcpu() <= 8 : run_reweight ( 3 ,  3000 , 123456 )
    ## else             : run_reweight ( 3 ,  3000 )
    
    run_reweight ( 3 , 3000 , 7000 , 8000 ) 
    
# =============================================================================
        
# =============================================================================
##                                                                      The END 
# =============================================================================
