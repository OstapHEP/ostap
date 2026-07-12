3#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/adversarial_vaildation.py
#  Implement adversarial vaildation to probe
#  the difference between two weighted datasets
#
#  @see ADVAL_LGBM , class based on lightgbm, the most CPU efficient 
#  @see ADVAL_HGBC , class based on HistGradientBoosterClassifier
#  @see ADVAL_GBC  , class based on GradientBoosterClassifier
#  @see ADVAL_RF   , class based on RandomForestClassifier 
#  @see ADVAL_XGB  , class based on XGBoost 
#  @see ADVAL_CATB , class based on CatBoost 
#
#  As t-value \f$ 100 \times \left( 1 - 2 \times AUC \right)^2\f$ is used
#  To estimate the~p-value permutations are used.
# 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-04
# =============================================================================
""" Implement adversarial vaildation to probe
    the difference between two weighted datasets

 - see ADVAL_LGBM , class based on lightgbm, the most CPU efficient 
 - see ADVAL_HGBC , class based on HistGradientBoosterClassifier
 - see ADVAL_GBC  , class based on GradientBoosterClassifier
 - see ADVAL_RF   , class based on RandomForestClassifier 
 - see ADVAL_XGB  , class based on XGBoost 
 - see ADVAL_CATB , class based on CatBoost 

As t-value 100 * ( 1 - 2 * AUC AUC )**2  is used
To estimate the~p-value permutations are used. 

"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2026-07-04"
__all__     = (
    'ADVAL_LGBM' , ## LightGBM-based class for Adversarial Validation for the differecne between two (weighted) dataset
    'ADVAL_XGB'  , ## XGBoost-based class for Adversarial Validation for the difference between two (weighted) dataset
    'ADVAL_HGBC' , ## HGBC-based class for Adversarial Validation for the difference between two (weighted) dataset
    'ADVAL_GBC'  , ## GBC-based class for Adversarial Validation for the difference between two (weighted) dataset
    'ADVAL_RF'   , ## RandomForst-based class for Adversarial Validation for the difference between two (weighted) dataset
)
# =============================================================================
from   ostap.core.ostap_types   import string_types
from   ostap.core.cpu_info      import HAS_AVX2
from   ostap.utils.basic        import typename , numcpu 
from   ostap.stats.gof_np       import GoFnp 
from   ostap.stats.gof_utils    import num_jobs, weight_trivial
import ROOT, numpy, abc, os   
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.adversarial_validation' )
else                       : logger = getLogger ( __name__ )
# =============================================================================
logger.debug ( 'Implement adversarial validation for Goodness-of-fit & Two-Samples test' )
# =============================================================================
## t-value from A
#  t-value is defined as \f$  100 \times \left( 1 - 2 \times AUC \right)^2 \f$ 
def tvaleu_from_AUC ( auc ) :
    """ t-value is defined as 100 * abs(1-2*AUC)**2
    """
    return 100 * ( 1.0 - 2.0 * auc ) ** 2 
# =============================================================================
## allow parallel run
#  - check "OMP_NUM_THREADS"
#  - check "MKL_NUM_THREADS"
#  - check "OPENBLAS_NUM_THREADS"
def run_parallel ( parallel ) : 
    """ Allow parallel run
    - check "OMP_NUM_THREADS"
    - check "MKL_NUM_THREADS"
    - check "OPENBLAS_NUM_THREADS"
    """
    ##
    if not parallel : return False
    ## 
    if   '1' != os.environ.get ( "OMP_NUM_THREADS"      , "" ) : return False 
    elif '1' != os.environ.get ( "MKL_NUM_THREADS"      , "" ) : return False 
    elif '1' != os.environ.get ( "OPENBLAS_NUM_THREADS" , "" ) : return False
    ## 
    return True 

# ==============================================================================
## @class ADVAL_base 
#  Base class for adversarial validation for the difference between two (weighted) dataset
class ADVAL_base (GoFnp):
    """ Base class for Adversarial validation for the differece between two (weighted) dataset
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  ,
                   method   = "Adversarial Validation" , **params ) :
        
        n_splits = params.pop ( 'n_splits' , 5 ) 
        assert isinstance ( n_splits , int ) and 2 <= n_splits , "Invalid n_splits:%s" % n_splits
        
        GoFnp.__init__ ( self                 ,
                         nToys     = nToys    ,
                         parallel  = parallel , 
                         silent    = silent   ,
                         progress  = progress ,
                         normalize = False    , 
                         method    = method   , **params )

        self.__n_splits = n_splits 

    # ============================================================================
    @property
    def n_splits ( self ) :
        """`n_splits`: # of splits for cross-validation"""
        return self.__n_splits 

    # ==================================================================================
    @property
    def config ( self ) :
        """`config` : get all configuratino parameters"""
        conf = super().config
        conf [ 'n_splits' ] = self.n_splits 
        return conf
    
    # =========================================================================
    ## Are weigths supported by this estimator?
    @property 
    def weights_supported ( self ) :
        """`weights_supported` : Are weigths supported by this estimator?"""
        return True 

    # =========================================================================
    ## Good for two-samples comparison?
    #  Can this estimator be used for comparison of two samples?
    @property 
    def two_samples ( self ) :
        """`two_samples`: Can this estimator be used for comparison of two samples?
        """
        return True 
        
    # ==========================================================================
    ## Calculate t-value for two non-structured (weighted) datasets 
    #   @param data1   the first dataset
    #   @param data2   the second dataset
    #   @param weight1 the first array of weights 
    #   @param weight3 the second array of weights
    #   tvalue is defined as \f$  100 \times \left( 1 - 2 \times AUC \right)^2 \f$ 
    def tvalue ( self              ,
                 data1             ,
                 data2             ,  * , 
                 weight1   = None  ,
                 weight2   = None  ,
                 normalize = False ) :
        """ Calculate t-value for two non-structured (weighted) arrays 
        - data1   : the first dataset
        - data2   : the second dataset
        - weight1 : the first array of weights 
        - weight3 : the second array of weights 
         t-value is defined as 100 * abs(1-2*AUC)**2
        """
        ## 
        sh1 = data1.shape 
        sh2 = data2.shape
        ## 
        assert 2 == len ( sh1 ) and 2 == len ( sh2 ) and sh1 [ 1 ] == sh2 [ 1 ] and sh1 [ 0 ] and sh2 [ 0 ] , \
            "Invalid arrays: %s , %s" % ( sh1 , sh2 )
        
        ## convert numpy arrays into pandas dataframes
        import pandas as pd 

        df_1 = pd.DataFrame ( data1 , dtype = float ) 
        df_2 = pd.DataFrame ( data2 , dtype = float )

        column_target = 'column_target'        
        column_weight = 'column_weight'
        
        df_1 [ column_target ] = 1
        df_2 [ column_target ] = 0

        w1_trivial = weight_trivial ( weight1 )
        w2_trivial = weight_trivial ( weight2 )
        
        if w1_trivial and w2_trivial : weights = False 
        else                         :            
            weights = True
            df_1 [ column_weight ] = 1 if w1_trivial else weight1 
            df_2 [ column_weight ] = 1 if w2_trivial else weight2 
            
        ## merge datasets together
        dataset = pd.concat ( [ df_1 , df_2 ] , axis = 0 ).reset_index ( drop = True )

        X = dataset . drop ( columns = [ column_target , column_weight ] if weights else [ column_target ] )
        Y = dataset [ column_target ]
        W = dataset [ column_weight ] if weights else None 
        N = len ( dataset )
        
        ## cross-validation
        from sklearn.model_selection import StratifiedKFold
        from sklearn.metrics         import roc_auc_score
        
        random_state = self.params.get ( 'random_state' )
        cv           = StratifiedKFold ( n_splits = self.n_splits , shuffle = True , random_state = random_state )
        oof_preds    = numpy.zeros( N )

        ## 
        for fold , ( train_idx , val_idx ) in enumerate ( cv.split ( X, Y ) ):

            X_train , Y_train , W_train = X.iloc [ train_idx ] , Y.iloc [ train_idx ] , W.iloc [ train_idx ] if weights else None 
            X_val   , Y_val   , W_val   = X.iloc [ val_idx   ] , Y.iloc [ val_idx   ] , W.iloc [ val_idx   ] if weights else None 
                            
            ## train model and make predictions & predict 
            oof_preds [ val_idx ] = self.work ( X_train , Y_train , W_train , X_val  , Y_val , W_val   )
            
        ## (weighted) ROC-AUC
        auc_score = roc_auc_score ( Y , oof_preds , sample_weight = W )

        ## 
        return tvaleu_from_AUC ( auc_score )
        
# =======================================================================================
## @class ADVAL_LGBM
#  LightGBM-based lass for Adversarial Validation for the difference between two (weighted) dataset
#  @see lightgbm
class ADVAL_LGBM (ADVAL_base) : 
    """ LightGBM-based class for Adversarial Validation for the differece between two (weighted) dataset
    - see lightgbm 
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  , **params ) :
        
        config = {
            'objective'        : 'binary' ,
            'metric'           : 'auc'    ,
            'learning_rate'    :  0.05    , 
            'num_leaves'       : 24       ,
            'max_depth'        :  4       ,
            'min_data_in_leaf' : 10       ,
            'verbose'          : -1       ,
        }
        ##
        
        # =================================================================================
        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "Parallel processing is switched OFF! (OMP/MKL/OPENBLAS)_NUM_THREADS" ) 
            parallel = False 
        
        ## Attention! 
        params [ 'num_threads' ] = 1 if parallel else num_jobs ( params , numcpu () - 1 )
        
        config.update ( params ) 
        ## 
        import lightgbm as LightGBM
        ADVAL_base.__init__ ( self, 
                              nToys    = nToys    ,
                              parallel = parallel ,
                              silent   = silent   , 
                              progress = progress , 
                              method   = "Adversarial Validation/LightGBM" , **config   ) 

    # ==================================================================================
    ## Train the model and make predictions
    def work ( self    ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   ) :
        """" Train the model and make predictions
        """
        
        ## 
        import lightgbm as LightGBM
        ## 
        train_data = LightGBM.Dataset ( X_train , label = Y_train , weight = W_train )
        val_data   = LightGBM.Dataset ( X_val   , label = Y_val   , weight = W_val   , reference = train_data )
        ## 
        ## train the model
        model = LightGBM.train (
            self.params ,
            train_data  ,
            num_boost_round = 500 ,
            valid_sets = [ val_data ],
            callbacks  = [ LightGBM.early_stopping ( stopping_rounds = 20 , verbose = False ) ]
        )
        
        ## predict the results 
        return model.predict ( X_val , num_iteration = model.best_iteration )
    
# =======================================================================================
## @class ADVAL_XGB
#  XGBoost-based lass for Adversarial Validation for the difference between two (weighted) dataset
#  @see xgboost 
class ADVAL_XGB (ADVAL_base) : 
    """ XGBoost-based class for Adversarial Validation for the differece between two (weighted) dataset
    - see xgboost
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  , **params ) :
                
        config = {
            'objective'     : 'binary:logistic' ,
            'eval_metric'   : 'auc'             ,
            'tree_method'   : 'hist'            , # Histogram method, similar to LigthGBM 
            'learning_rate' : 0.05              ,
            'max_depth'     : 5                 ,
            'seed'          : 42
        }
        ##

        ## Attention! 
        params [ 'n_jobs' ] = 1 if parallel else num_jobs ( params , numcpu() - 1 )
        
        config.update ( params ) 

        import xgboost as XGBoost             
        ADVAL_base.__init__ ( self, 
                              nToys    = nToys    ,
                              parallel = parallel ,
                              silent   = silent   , 
                              progress = progress , 
                              method   = "Adversarial Validation/XGBoost" , **config ) 


    # ==================================================================================
    ## Train the model and make predictions
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   ) :
        """" Train the model and make predictions
        """

        ## 
        import xgboost as XGBoost 

        dtrain = XGBoost.DMatrix ( X_train , label = Y_train , weight = W_train )
        dval   = XGBoost.DMatrix ( X_val   , label = Y_val   , weight = W_val   )
        
        # Train the model 
        model = XGBoost.train (
            self.params                  ,
            dtrain                       ,
            num_boost_round = 500        ,
            evals = [ ( dval , 'val' ) ] ,
            early_stopping_rounds = 20   ,
            verbose_eval = False 
        )
        
        # predict 
        return model.predict ( dval ) ## , iteration_range= ( 0 , model.best_iteration + 1 ) )

# ======================================================================================
if HAS_AVX2 : 
    # ==================================================================================
    ## @class ADVAL_CATB
    #  CatBoost-based lass for Adversarial Validation for the difference between two (weighted) dataset
    #  @see catboost 
    class ADVAL_CATB (ADVAL_base) : 
        """ CatBoost-based class for Adversarial Validation for the differece between two (weighted) dataset
        - see Catboost
        """
        def __init__ ( self             ,
                       nToys    = 400   ,
                       parallel = False ,
                       silent   = False ,
                       progress = True  , **params ) :
            
            config = { 'iterations'    : 500    ,
                       'learning_rate' : 0.05   ,
                       'depth'         : 5      ,
                       'eval_metric'   : 'AUC'  ,
                       'verbose'       : False  }
            
            
            # =================================================================================
            if parallel and not run_parallel ( parallel ) :
                logger.warning ( "Parallel processing is switched OFF! (OMP/MKL/OPENBLAS)_NUM_THREADS" ) 
                parallel = False 

            ## Attention! 
            params [ 'thread_count' ] = 1 if parallel else num_jobs ( params , numcpu () - 1 )  

            
            if   'random_seed'  in params :                            params.pop ( 'random_state'      )
            elif 'random_state' in params : params [ 'random_seed' ] = params.pop ( 'random_state' , 42 )
            
            ## 
            import catboost as CatBoost 
            config.update ( params ) 
            ADVAL_base.__init__ ( self, 
                                  nToys    = nToys    ,
                                  parallel = parallel ,
                                  silent   = silent   , 
                                  progress = progress , 
                                  method   = "Adversarial Validation/CatBoost" , **config   ) 
            
                
        # ==================================================================================
        ## Train the model and make predictions
        def work ( self ,
                   X_train , Y_train , W_train ,
                   X_val   , Y_val   , W_val   ) :
            """ Train the model and make predictions
            """
            
            import catboost as CatBoost 
            
            train_pool = CatBoost.Pool ( data = X_train , label = Y_train , weight = W_train )
            val_pool   = CatBoost.Pool ( data = X_val   , label = Y_val   , weight = W_val   )
            
            ## create the model 
            model = CatBoost.CatBoostClassifier( **self.params )
            
            # train it 
            model.fit ( train_pool , eval_set=val_pool, early_stopping_rounds = 20 )
            
            return model.predict_proba(X_val)[:, 1]
        
    __all__ += ( 'ADVAL_CATB' , )
        
# =============================================================================
## @class ADVAL_HGBC
#  Class for adversarial validation for the difference between two (weighted) dataset
#  based HistoGradientBoosterClassifier 
#  @see HistoGradientBoosterClassifier 
class ADVAL_HGBC (ADVAL_base) : 
    """ HGBC-based class for Adversarial validation for the differece between two (weighted) dataset
    @see HistoGradientBoosterClassifier 
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  , **params ) :
        
        config =  { 'learning_rate'    : 0.05  ,
                    'max_depth'        : 5     ,
                    'max_iter'         : 200   ,  ## number of trees (num_boost_round)
                    'early_stopping'   : False ,  ## embedded early stopping
                    'n_iter_no_change' : 15    } 
        
        # =================================================================================
        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "Parallel processing is switched OFF! (OMP/MKL/OPENBLAS)_NUM_THREADS" ) 
            parallel = False 
                
        config.update ( params )
        
        import sklearn.ensemble 
        ADVAL_base.__init__ ( self, 
                              nToys    = nToys    ,
                              parallel = parallel ,
                              silent   = silent   , 
                              progress = progress , 
                              method   = "Adversarial Validation/HGBC" , **config   ) 

    # ==================================================================================
    ## Train the model and make predictions
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   ) :
        """" Train the model and make predictions
        """
        
        from sklearn.ensemble import HistGradientBoostingClassifier as CLASSIFIER
        ## 
        model = CLASSIFIER ( **self.params )
        model.fit ( X_train , Y_train , sample_weight = W_train )
        
        return model.predict_proba ( X_val ) [ : , 1 ]
            
# =============================================================================
## @class ADVAL_GBC
#  Class for adversarial validation for the difference between two (weighted) dataset
#  based GradientBoosterClassifier 
#  @see GradientBoosterClassifier 
class ADVAL_GBC (ADVAL_base) : 
    """ GBC-based class for Adversarial validation for the differece between two (weighted) dataset
    @see GradientBoosterClassifier 
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  , **params   ) :
        
        config =  { 'learning_rate'       : 0.05 ,
                    'max_depth'           : 5    ,
                    'n_estimators'        : 100  , ## number of trees 
                    'validation_fraction' : 0.10 , ## internal validation for early stopping
                    'n_iter_no_change'    : 15   }

        # =================================================================================
        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "Parallel processing is switched OFF! (OMP/MKL/OPENBLAS)_NUM_THREADS" ) 
            parallel = False 
        
        config.update ( params ) 

        import sklearn.ensemble 
        ADVAL_base.__init__ ( self, 
                              nToys    = nToys    ,
                              parallel = parallel ,
                              silent   = silent   , 
                              progress = progress , 
                              method   = "Adversarial Validation/GBC" , **config   ) 

    # ==================================================================================
    ## Train the model and make predictions
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   ) :
        """" Train the model and make predictions
        """
        ## 
        from sklearn.ensemble import GradientBoostingClassifier as CLASSIFIER 
        ##
        model = CLASSIFIER ( **self.params )
        ##
        model.fit ( X_train , Y_train , sample_weight = W_train )
        ## 
        return model.predict_proba ( X_val ) [ : , 1 ]

# =============================================================================
## @class ADVAL_RF
#  Class for adversarial validation for the difference between two (weighted) dataset
#  based on RandomForestClassifier 
#  @see RandomForestClassifier 
class ADVAL_RF (ADVAL_base) : 
    """ GBC-based class for Adversarial validation for the differece between two (weighted) dataset
    @see GradientBoosterClassifier 
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  , **params   ) :

        config =  {
            'n_estimators' : 100 ,
            'max_depth'    :   6 ,  
            }
        
        # =================================================================================
        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "Parallel processing is switched OFF! (OMP/MKL/OPENBLAS)_NUM_THREADS" ) 
            parallel = False 

        ## Attention! 
        params [ 'n_jobs' ] = 1 if parallel else num_jobs ( params , numcpu() - 1  ) 
        
        config.update ( params ) 

        import sklearn.ensemble 
        ADVAL_base.__init__ ( self, 
                              nToys    = nToys    ,
                              parallel = parallel ,
                              silent   = silent   , 
                              progress = progress , 
                              method   = "Adversarial Validation/RandomForest" , **config  ) 

    # ==================================================================================
    ## Train the model and make predictions
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   ) :
        """" Train the model and make predictions
        """
        
        from sklearn.ensemble import RandomForestClassifier as CLASSIFIER 
        ## 
        model = CLASSIFIER ( **self.params )
        model.fit ( X_train , Y_train , sample_weight = W_train )
        ## 
        return model.predict_proba ( X_val ) [ : , 1 ]


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import pandas # =======================================================
        # =====================================================================
    except ImportError : # ====================================================
        # =====================================================================
        logger.error ( 'pandas cannot be imported!' ) # =======================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================                       
        from sklearn.ensemble        import HistGradientBoostingClassifier
        from sklearn.ensemble        import GradientBoostingClassifier
        from sklearn.model_selection import StratifiedKFold
        from sklearn.metrics         import roc_auc_score
        # =====================================================================
    except ImportError : # ====================================================
        # =====================================================================
        logger.error ( "sklearn cannot be imported!" ) # ======================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import lightgbm
        # =====================================================================
    except ImportError : # ====================================================
        # =====================================================================
        logger.error ( "lightgbm cannot be imported!" ) # =====================
        # =====================================================================
    try : # ===================================================================
        # =====================================================================
        import xgboost 
        # =====================================================================
    except ImportError : # ====================================================
        # =====================================================================
        logger.error ( "xgboost cannot be imported!" ) # ======================
    # =========================================================================
    if HAS_AVX2 : # ===========================================================
        # =====================================================================
        try : # ===============================================================
            # =================================================================
            import catboost 
            # =================================================================
        except ImportError : # ================================================
            # =================================================================
            logger.error ( "catboost cannot be imported!" ) # =================
        
# =============================================================================
##                                                                      The END 
# =============================================================================
