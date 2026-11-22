#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/adversarial_validation.py
#  Implement adversarial validation (regression mode) to probe
#  the difference between two weighted datasets (supporting negative weights via sPlot)
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-04
# =============================================================================
""" Implement adversarial validation (regression mode) to probe
    the difference between two weighted datasets (supporting negative weights via sPlot)

 - see ADVAL_LGBM        , class based on LightGBM
 - see ADVAL_XGB         , class based on XGBoost 
 - see ADVAL_CATB        , class based on CatBoost 
 - see ADVAL_HGBC        , class based on HistGradientBoostingRegressor
 - see ADVAL_GBC         , class based on GradientBoostingRegressor
 - see ADVAL_RF          , class based on RandomForestRegressor 
 - see ADVAL_TORCH_ADAM  , class based on PyTorch (Regression)
 - see ADVAL_TORCH_LBFGS , class based on PyTorch (Regression)

As t-value: based on Weighted MSE / RMSE is used.
To estimate the p-value, permutations are used. 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2026-07-04"
__all__     = (
    'ADVAL_LGBM'        , 
    'ADVAL_XGB'         , 
    'ADVAL_CATB'        , 
    'ADVAL_HGBC'        , 
    'ADVAL_GBC'         , 
    'ADVAL_RF'          , 
    'ADVAL_TORCH_ADAM'  , 
    'ADVAL_TORCH_LBFGS' , 
)
# =============================================================================
from   collections            import defaultdict 
from   ostap.core.ostap_types import string_types, sequence_types, integer_types 
from   ostap.utils.core       import typename
from   ostap.utils.basic      import classprop 
from   ostap.stats.utils      import ( weight_trivial , nEff         , 
                                       num_samples    , num_features ,
                                       check_all      ,
                                       optimal_batch_size            )
from   ostap.utils.basic      import num_jobs, run_parallel 
from   ostap.stats.gof_np     import GoFnp
from   ostap.logger.pretty    import nice_print 
from   ostap.stats.tools      import hasSkLearn 
## 
import ostap.stats.utils      as     U
import numpy, abc
# =============================================================================
# Logging setup 
# =============================================================================
from ostap.logger.logger import getLogger 
logger = getLogger ( 'ostap.stats.adval' )
# =============================================================================
logger.debug ( 'Implement Adversarial Validation (Regression mode)' )
# =============================================================================
has_sklearn = hasSkLearn()
from   sklearn.metrics        import mean_squared_error
# =============================================================================
DEFAULT_ESTIMATORS         = 500
DEFAULT_EARLY_STOPPING     =  50 
MAX_DEPTH                  =   5
# =============================================================================
## convert MSE to t-value 
def tvalue_from_MSE ( mse ) :
    """ Convert MST to tvalue"
    - MSE = 0.25  means no difference 
    """
    return 1 - 4 * float ( mse ) 

# =============================================================================
## Helper to handle negative weights for regression-based models via target inversion
def transform_weights_and_targets ( Y , W ) :
    """ Transform negative weights into positive ones by flipping the regression target 
        Y -> (1 - Y) when W < 0, preserving the physical meaning for MSE optimization.
    """
    if W is None: return Y, None
    
    w_arr = numpy.asarray ( W , dtype = numpy.float32 )
    signs = numpy.sign ( w_arr )
    signs [ signs == 0 ] = 1.0
    
    Y_mod = numpy.where ( signs < 0 , 1.0 - Y , Y )
    W_mod = numpy.abs   ( w_arr )
    
    return Y_mod, W_mod

# =============================================================================
## @class ADVAL_base 
#  Base class for adversarial validation (regression mode)
class ADVAL_base (GoFnp):
    """ Base class for Adversarial Validation (regression mode) for weighted datasets
    """
    def __init__ ( self              ,
                   nToys      = 400  ,
                   parallel   = True ,
                   normalize  = True ,
                   n_splits   = 5    , **params ) :

        # =====================================================================
        if   not isinstance ( n_splits , int ) : raise TypeError ( "Invalid `n_splits` type  : %s" % typename ( n_splits ) )
        elif n_splits < 1                      : n_splits = 1 
        ## store it! 
        self.__n_splits = n_splits
        
        ## switch off parallel processing
        parallel = True if parallel else False 
        if parallel and not run_parallel ( parallel ) :
            logger.info ( '%s: (internal) parallel processing is OFF' % typename ( self ) )
            parallel = False
            
        ## (re)define n_jobs     
        params [ 'n_jobs'    ] = 1 if parallel else num_jobs ( params )        
        params [ 'parallel'  ] = parallel
        params [ 'normalize' ] = True if normalize else False 

        self.__n_splits            = n_splits 
        self.__importance_features = {}

        GoFnp.__init__ ( self            ,
                         nToys  = nToys  , **params )                
    # =========================================================================
    @property
    def n_splits ( self ) :
        """`n_splits`: Number of splits for cross-validation
        """
        return self.__n_splits 

    # =========================================================================
    @property
    def config ( self ) :
        """`config`: Get all configuration parameters"""
        conf = {}
        conf.update ( super().config ) 
        conf [ 'n_splits'    ] = self.n_splits
        return conf
    
    @property 
    def weights_supported ( self ) :
        return True 

    @property 
    def two_samples ( self ) :
        return True

    # =========================================================================
    ## max-dimensionality
    @property 
    def max_dimensionality(self):
        """`max_dimensionality` : maximal dimensionality of dataset
        """
        return -1

    ## importance features (if suppported by backend) 
    @property 
    def importance_features ( self ) :
        """`importance_features`: Importance features (if suppported by backend)
        """
        return self.__importance_features

    ## prepare the table of importance features 
    def importance_table  ( self ,
                            title  = '' ,
                            prefix = '' ,
                            style  = '' ) :
        """ Prepare the table of importance features
        """
        
        rows  = [ ( 'Feature/#' , 'Importance [%]' ) ]
        rows += [ ( str ( feature ) , '%.2f' % float ( gain ) ) for feature, gain in self.importance_features.items () ] 
        title = title if title else "%s importance" % self.method
        import ostap.logger.table as T
        return T.table ( rows               ,
                         title     = title  ,
                         prefix    = prefix ,
                         alignment = 'cc'   ,  
                         style     = style  )
    
    @abc.abstractmethod 
    def work ( self    ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        return NotImplemented

    
    def tvalue ( self               ,
                 data1              ,
                 data2              ,  * , 
                 weight1    = None  ,
                 weight2    = None  ,
                 normalize  = True  ,
                 importance = False ) :
        
        data1, data2 = self.unpack ( data1 , data2 )

        ## check input data 
        check_all ( data1 , data2 , weight1 , weight2 , typename ( self ) )
        
        if self.normalize and normalize :
            uds1 , uds2 = self.normalize_pooled ( data1   = data1   ,
                                                  data2   = data2   ,
                                                  weight1 = weight1 ,
                                                  weight2 = weight2 )
            return self.tvalue ( uds1       , uds2       ,
                                 weight1    = weight1    ,
                                 weight2    = weight2    ,
                                 normalize  = False      ,
                                 importance = importance )
            
        N1 = num_samples ( data1 )
        N2 = num_samples ( data2 ) 
        N  = N1 + N2

        X = numpy.vstack ( [ data1 , data2 ] )
        Y = numpy.zeros  ( N , dtype = numpy.float32 )
        Y [ : N1 ] = 1.0  # Test gets 1, Train gets 0

        w1_trivial = weight_trivial ( weight1 )
        w2_trivial = weight_trivial ( weight2 )

        if w1_trivial and w2_trivial:
            weights = False
            W       = None
        else:
            weights = True
            w1      = numpy.ones ( N1 , dtype = numpy.float32 ) if w1_trivial else numpy.asarray ( weight1 , dtype = numpy.float32 )
            w2      = numpy.ones ( N2 , dtype = numpy.float32 ) if w2_trivial else numpy.asarray ( weight2 , dtype = numpy.float32 )

            sumw1   = numpy.sum ( w1 , dtype = numpy.float64 )
            sumw2   = numpy.sum ( w2 , dtype = numpy.float64 )            
            w2      = w2 * numpy.float32 ( sumw1 / sumw2 ) 
                        
            W       = numpy.concatenate ( [ w1 , w2 ] )

        from sklearn.model_selection import KFold
        
        random_state = self.params.get ( 'random_state' , None )
        
        oof_preds    = numpy.zeros ( N , dtype = numpy.float32 )
        importances  = numpy.zeros ( X.shape [ 1 ] , dtype = numpy.float32 ) if importance else None

        if 1 == self.n_splits :
            full_idx = numpy.arange ( N , dtype = numpy.intp )
            splits   = [ ( full_idx , full_idx ) ]
        else :
            from sklearn.model_selection import KFold
            cv     = KFold ( n_splits = self.n_splits , shuffle = True , random_state = random_state )
            splits = cv.split ( X , Y )
            
        for fold, ( train_idx , val_idx ) in enumerate ( splits ) :
            
            X_train, Y_train_f = X[train_idx], Y[train_idx]
            X_val,   Y_val_f   = X[val_idx],   Y[val_idx]
            
            W_train = W[train_idx] if weights else None
            W_val   = W[val_idx]   if weights else None

            fold_predictions, fold_importances = self.work ( X_train , Y_train_f , W_train ,
                                                             X_val   , Y_val_f   , W_val   ,
                                                             importance = importance     )

            oof_preds [ val_idx ] = fold_predictions

            if importance and importances is not None and fold_importances is not None:
                importances += fold_importances

        if importance and importances is not None and 0 < numpy.sum ( importances ) :
            imp_sum = numpy.sum ( importances )
            importances *= 100.0 / imp_sum
            
            feature_names = [f'f{i}' for i in range(X.shape[1])]
            sorted_pairs = sorted(zip(feature_names, importances), key=lambda x: x[1], reverse=True)
            self.__importance_features = {feat: val for feat, val in sorted_pairs}


        Y_eval, eval_weights = transform_weights_and_targets ( Y , W )
        ## mse_score = mean_squared_error ( Y , oof_preds , sample_weight = eval_weights )
        mse_score = mean_squared_error ( Y_eval , oof_preds , sample_weight = eval_weights )

        return tvalue_from_MSE ( mse_score )

# =======================================================================================
## @class ADVAL_LGBM (Regression)
class ADVAL_LGBM (ADVAL_base) :
    
    def __init__ ( self             ,
                   nToys    = 400   , **params ) :
        
        config = {
            'objective'             : 'regression'       ,
            'metric'                : 'rmse'             ,
            ## 
            'n_estimators'          : DEFAULT_ESTIMATORS     ,
            'max_depth'             : MAX_DEPTH              ,
            'early_stopping_rounds' : DEFAULT_EARLY_STOPPING , 
            'learning_rate'         :  0.03                  ,
            ## 
            'num_leaves'           : 24                 ,
            'min_child_samples'    : 50                 ,
            'reg_alpha'            : 0.1                ,
            'reg_lambda'           : 1.0                ,
            'force_row_wise'       : True               ,
            'force_col_wise'       : False              ,
            'max_bin'              : 1023               , 
            'n_jobs'               : -1                 ,  
            'verbosity'            : -1                 ,
            'verbose'              : -1                 ,
        }
        
        config.update ( params ) 

        row_wise = bool ( config.get ( 'force_row_wise' , True  ) ) 
        col_wise = bool ( config.get ( 'force_col_wise' , False ) ) 

        if row_wise ^ col_wise : pass 
        else :
            row_wise, col_wise = True, False 

        config [ 'force_col_wise' ] = col_wise 
        config [ 'force_row_wise' ] = row_wise

        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' , DEFAULT_ESTIMATORS ) 
        if 'iterations'      in config : config [ 'n_estimators' ] = config.pop ( 'iterations'      , DEFAULT_ESTIMATORS ) 
        if 'max_iter'        in config : config [ 'n_estimators' ] = config.pop ( 'max_iter'        , DEFAULT_ESTIMATORS ) 
            
        config [ 'normalize'      ] = False         
        ADVAL_base.__init__ ( self  , nToys = nToys , **config   )
        
        if 'n_jobs' in self.params : self.params [ 'num_threads' ] = self.params.pop ( 'n_jobs' , 1 )
        if self.silent :
            self.params [ 'verbosity'   ] = -1 
            self.params [ 'verbose'     ] = -1 

    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_LGBM  
    
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_LGBM

    def work ( self    ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        
        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        Y_val_mod,   W_val_mod   = transform_weights_and_targets ( Y_val  , W_val   )

        params = {}
        params.update ( self.params )

        nf =  num_features ( X_train )
        ns =  num_samples  ( X_train )

        # ====================================================================
        ## num-boost-rounds & early stopping 
        # ====================================================================
        num_boost_round  = ( params.pop ( 'iterations'      , None ) or
                             params.pop ( 'max_iter'        , None ) or
                             params.pop ( 'num_boost_round' , None ) or
                             params.pop ( 'n_estimators'    , None ) or DEFAULT_ESTIMATORS )         
        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , DEFAULT_EARLY_STOPPING )
        if   early_stopping_rounds is None             : early_stopping_rounds = 0 
        elif early_stopping_rounds <= 0                : early_stopping_rounds = 0
        elif num_boost_round <= early_stopping_rounds  : early_stopping_rounds = 0
        

        import lightgbm as LightGBM
        callbacks = []
        if 0 < early_stopping_rounds :
            callbacks.append ( LightGBM.early_stopping ( stopping_rounds = early_stopping_rounds , verbose = False ) )

        train_data = LightGBM.Dataset ( X_train , label = Y_train_mod , weight = W_train_mod , free_raw_data = False , params    = params     )
        val_data   = LightGBM.Dataset ( X_val   , label = Y_val_mod   , weight = W_val_mod   , free_raw_data = False , reference = train_data )
        
        model      = LightGBM.train ( params          = params          ,
                                      train_set       = train_data      ,
                                      num_boost_round = num_boost_round ,
                                      valid_sets      = [ val_data ]    ,
                                      callbacks       = callbacks       )
        
        predictions = model.predict ( X_val , num_iteration = model.best_iteration )
        imps = model.feature_importance ( importance_type = 'gain') if importance else None
        
        return predictions, imps

# =======================================================================================
## @class ADVAL_XGB (Regression)
class ADVAL_XGB (ADVAL_base) : 
    def __init__ ( self             ,
                   nToys    = 400   , **params ) :
        
        config = {  'objective'             : 'reg:squarederror',
                    'eval_metric'           : 'rmse',
                    'tree_method'           : 'hist',
                    ##
                    'n_estimators'          : DEFAULT_ESTIMATORS     ,
                    'max_depth'             : MAX_DEPTH              ,
                    'early_stopping_rounds' : DEFAULT_EARLY_STOPPING , 
                    'learning_rate'         : 0.03                   ,
                    ## 
                    'max_leaves'       : 15    ,
                    'min_child_weight' :  1.0  ,
                    'subsample'        :  0.8  ,
                    'colsample_bytree' :  0.8  ,
                    'alpha'            :  0.1  ,
                    'lambda'           :  1.0  ,
                    'n_jobs'           : -1    ,
                    'verbosity'        :  0
                  }
        config.update ( params ) 

        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' , DEFAULT_ESTIMATORS ) 
        if 'iterations'      in config : config [ 'n_estimators' ] = config.pop ( 'iterations'      , DEFAULT_ESTIMATORS ) 
        if 'max_iter'        in config : config [ 'n_estimators' ] = config.pop ( 'max_iter'        , DEFAULT_ESTIMATORS ) 
            
        config [ 'normalize'      ] = False         
        ADVAL_base.__init__ ( self , nToys = nToys , **config   ) 
        
        if 'n_jobs' in self.params : self.params [ 'nthread' ] = self.params.pop ( 'n_jobs' , 1 )            
        if self.silent : self.params [ 'verbosity' ] = 0

        
    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_XGB 
    
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_XGB 

    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        
        import xgboost as XGBoost             

        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        Y_val_mod,   W_val_mod   = transform_weights_and_targets ( Y_val  , W_val   )

        dtrain = XGBoost.DMatrix ( X_train , label = Y_train_mod , weight = W_train_mod )
        dval   = XGBoost.DMatrix ( X_val   , label = Y_val_mod   , weight = W_val_mod   )
        evals  = [ ( dtrain , 'train' ) , ( dval , 'val' ) ]

        params = {}
        
        params.update ( self.params )

        # ====================================================================
        ## num-boost-rounds & early stopping 
        # ====================================================================
        num_boost_round  = ( params.pop ( 'iterations'      , None ) or
                             params.pop ( 'max_iter'        , None ) or
                             params.pop ( 'num_boost_round' , None ) or
                             params.pop ( 'n_estimators'    , None ) or DEFAULT_ESTIMATORS )         
        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , DEFAULT_EARLY_STOPPING )
        if   early_stopping_rounds is None             : early_stopping_rounds = 0 
        elif early_stopping_rounds <= 0                : early_stopping_rounds = 0
        elif num_boost_round <= early_stopping_rounds  : early_stopping_rounds = 0

        
        model = XGBoost.train ( params                = params,
                                dtrain                = dtrain,
                                num_boost_round       = num_boost_round,
                                evals                 = evals,
                                early_stopping_rounds = early_stopping_rounds,
                                verbose_eval          = False )    

        predictions = model.predict ( dval )
        
        if importance :
            score_dict = model.get_score ( importance_type = 'gain' )
            if hasattr ( X_train , 'columns' ) :
                feature_names = X_train.columns.tolist()
            else :
                feature_names = [ f'f{i}' for i in range ( num_features ( X_train ) ) ]
            imps = numpy.array ( [ score_dict.get ( col , 0.0 ) for col in feature_names ] , dtype = numpy.float32 )
        else :
            imps = None

        return predictions, imps

# ======================================================================================
## @class ADVAL_CATB (Regression)
class ADVAL_CATB (ADVAL_base) : 
    def __init__ ( self             ,
                   nToys    = 400   , **params ) :
        
        config = {  'loss_function'         : 'RMSE',
                    'eval_metric'           : 'RMSE',
                    ## 
                    'n_estimators'          : DEFAULT_ESTIMATORS     , 
                    'depth'                 : MAX_DEPTH              ,
                    'early_stopping_rounds' : DEFAULT_EARLY_STOPPING ,
                    'learning_rate'         : 0.03,
                    ## 
                    'min_data_in_leaf'      :  20   ,
                    'l2_leaf_reg'           :   3.0 ,
                    'subsample'             :   0.8 ,
                    'bootstrap_type'        : 'Bernoulli',
                    'verbose'               : False
            }

        
        if   'random_seed'  in params : pass
        elif 'random_state' in params : params [ 'random_seed' ] = params.pop ( 'random_state' , None  )

        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' , DEFAULT_ESTIMATORS ) 
        if 'iterations'      in config : config [ 'n_estimators' ] = config.pop ( 'iterations'      , DEFAULT_ESTIMATORS ) 
        if 'max_iter'        in config : config [ 'n_estimators' ] = config.pop ( 'max_iter'        , DEFAULT_ESTIMATORS ) 
                    
        config.update ( params )

        config [ 'normalize'      ] = False                 
        ADVAL_base.__init__ ( self,  nToys = nToys , **config   ) 

        ## 
        if 'n_jobs' in self.params :
            self.params [ 'thread_count' ] = self.params.pop ( 'n_jobs' , 1 )
            
    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_CATB 
    
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_CATB

    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        
        import catboost as CatBoost
        
        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        Y_val_mod,   W_val_mod   = transform_weights_and_targets ( Y_val  , W_val   )

        train_pool = CatBoost.Pool ( X_train , label = Y_train_mod , weight = W_train_mod )
        val_pool   = CatBoost.Pool ( X_val   , label = Y_val_mod   , weight = W_val_mod   )

        params = {}
        params.update ( self.params )

        # ====================================================================
        ## num-boost-rounds & early stopping 
        # ====================================================================
        num_boost_round  = ( params.pop ( 'iterations'      , None ) or
                             params.pop ( 'max_iter'        , None ) or
                             params.pop ( 'num_boost_round' , None ) or
                             params.pop ( 'n_estimators'    , None ) or DEFAULT_ESTIMATORS )         
        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , DEFAULT_EARLY_STOPPING )
        if   early_stopping_rounds is None             : early_stopping_rounds = 0 
        elif early_stopping_rounds <= 0                : early_stopping_rounds = 0
        elif num_boost_round <= early_stopping_rounds  : early_stopping_rounds = 0
        
        iterations = num_boost_round

        params [ 'iterations'            ] = iterations
        params [ 'early_stopping_rounds' ] = early_stopping_rounds  
        params [ 'use_best_model'        ] = True if early_stopping_rounds else False
            
        model = CatBoost.CatBoostRegressor ( **params ) 
        model.fit ( train_pool , eval_set = val_pool )

        predictions = model.predict                ( val_pool )
        imps        = model.get_feature_importance ( data = train_pool ) if importance else None

        return predictions, imps

# =============================================================================
## @class ADVAL_HGBC (Regression)
class ADVAL_HGBC (ADVAL_base) : 
    def __init__ ( self           ,
                   nToys    = 400 , **params ) :
        
        config = {  'loss'              : 'squared_error',
                    'learning_rate'     : 0.03   ,
                    'max_iter'          : DEFAULT_ESTIMATORS ,
                    'max_depth'         : MAX_DEPTH ,
                    'max_leaf_nodes'    : 31    ,
                    'min_samples_leaf'  : 20    ,
                    'l2_regularization' :  0.1  ,
                    'early_stopping'    :  True ,
                    'n_iter_no_change'  : 20    ,
                    'scoring'           : 'neg_root_mean_squared_error',
                }

        config.update ( params )

        
        import sklearn
        from   packaging.version        import Version
        if Version ( sklearn.__version__ ) < Version ( '1.0' ) :
            if 'squared_error'               == config.get ( 'loss'   ) : config [ 'loss'    ] = 'least_squares' 
            if 'neg_root_mean_squared_error' == config.get ( 'scoring') : config [ 'scoring' ] = 'loss'  
            
        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' , DEFAULT_ESTIMATORS ) 
        if 'iterations'      in config : config [ 'n_estimators' ] = config.pop ( 'iterations'      , DEFAULT_ESTIMATORS ) 
        if 'max_iter'        in config : config [ 'n_estimators' ] = config.pop ( 'max_iter'        , DEFAULT_ESTIMATORS ) 
        
        config [ 'normalize'  ] = False                 
        ADVAL_base.__init__ ( self , nToys = nToys , **config   ) 

        if 'n_jobs' in self.params : self.params.pop ( 'n_jobs' , None )
        

    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :

        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_HGBC
    
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_HGBC


    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        
        from sklearn.ensemble import HistGradientBoostingRegressor
        from threadpoolctl    import threadpool_limits

        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        w_tr_arr = W_train_mod.values if hasattr ( W_train_mod , 'values' ) else W_train_mod

        nf =  num_features ( X_train )
        ns =  num_samples  ( X_train )
        
        params = {}
        params.update ( self.params )

        # ====================================================================
        ## num-boost-rounds & early stopping 
        # ====================================================================
        num_boost_round  = ( params.pop ( 'iterations'      , None ) or
                             params.pop ( 'max_iter'        , None ) or
                             params.pop ( 'num_boost_round' , None ) or
                             params.pop ( 'n_estimators'    , None ) or DEFAULT_ESTIMATORS )         
        
        max_iter = num_boost_round
        
        
        params [ 'max_iter' ] = max_iter

        X_tr_c  = numpy.ascontiguousarray ( X_train     , dtype = numpy.float32 )
        X_val_c = numpy.ascontiguousarray ( X_val       , dtype = numpy.float32 )
        Y_tr_c  = numpy.ascontiguousarray ( Y_train_mod , dtype = numpy.float32 )
        w_tr_c  = numpy.ascontiguousarray ( w_tr_arr    , dtype = numpy.float32 ) if w_tr_arr is not None else None

        with threadpool_limits ( limits = 1 ):
    
            model = HistGradientBoostingRegressor ( **params )
            model.fit ( X_tr_c , Y_tr_c , sample_weight = w_tr_c )
            predictions = model.predict ( X_val_c )
            
        return predictions, None 
            
# =============================================================================
## @class ADVAL_GBC (Regression)
class ADVAL_GBC (ADVAL_base) : 
    def __init__ ( self             ,
                   nToys    = 400   , **params   ) :
        
        config = {  'loss'              : 'squared_error',
                    'learning_rate'     : 0.05  ,
                    'n_estimators'      : DEFAULT_ESTIMATORS ,
                    'max_depth'         : MAX_DEPTH ,
                    'min_samples_split' :  10   ,
                    'min_samples_leaf'  :   5   ,
                    'subsample'         :   1.0 ,
                    'max_features'      :   1.0 ,
                }
        
        config.update ( params ) 
        
        import sklearn
        from   packaging.version        import Version
        if Version ( sklearn.__version__ ) < Version ( '1.0' ) :
            if 'squared_error' == config.get ( 'loss' ) : config [ 'loss'    ] = 'ls' 

        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' , DEFAULT_ESTIMATORS ) 
        if 'iterations'      in config : config [ 'n_estimators' ] = config.pop ( 'iterations'      , DEFAULT_ESTIMATORS ) 
        if 'max_iter'        in config : config [ 'n_estimators' ] = config.pop ( 'max_iter'        , DEFAULT_ESTIMATORS ) 
            
        config [ 'normalize'  ] = False                 
        ADVAL_base.__init__ ( self  , nToys = nToys , **config ) 
        
        if 'n_jobs' in self.params : self.params.pop ( 'n_jobs' , None )
        
    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_GBC 
    
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_GBC

    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        
        from sklearn.ensemble import GradientBoostingRegressor
        from threadpoolctl    import threadpool_limits

        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        w_tr_arr = W_train_mod.values if hasattr(W_train_mod, 'values') else W_train_mod

        nf =  num_features ( X_train )
        ns =  num_samples  ( X_train )
        
        params = {}
        params.update ( self.params )

        # ====================================================================
        ## num-boost-rounds & early stopping 
        # ====================================================================
        num_boost_round  = ( params.pop ( 'iterations'      , None ) or
                             params.pop ( 'max_iter'        , None ) or
                             params.pop ( 'num_boost_round' , None ) or
                             params.pop ( 'n_estimators'    , None ) or DEFAULT_ESTIMATORS )         
        
        n_estimators = num_boost_round

                                                     
        params [ 'n_estimators' ] = n_estimators
        
        with threadpool_limits ( limits = 1 ):
            
            model = GradientBoostingRegressor ( **params )
            model.fit ( X_train , Y_train_mod , sample_weight = w_tr_arr )            
            predictions = model.predict ( X_val )
            imps = model.feature_importances_ if importance else None

        return predictions, imps

# =============================================================================
## @class ADVAL_RF (Regression)
class ADVAL_RF (ADVAL_base) : 
    def __init__ ( self           ,
                   nToys    = 400 , **params   ) :
  
        config = {  'n_estimators'      : DEFAULT_ESTIMATORS ,
                    'n_jobs'            : -1    ,
                    'criterion'         : 'squared_error',
                    'max_depth'         : MAX_DEPTH ,
                    'min_samples_split' : 10    ,
                    'min_samples_leaf'  :  5    ,
                    'max_features'      :  1.0  ,
                    'bootstrap'         :  True ,
                    'max_samples'       :  0.8  ,
                }  
        config.update ( params )
        
        import sklearn
        from   packaging.version        import Version
        if Version ( sklearn.__version__ ) < Version ( '1.0' ) :
            if 'squared_error' == config.get ( 'criterion' ) : config [ 'criterion' ] = 'mse' 
        
        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' , DEFAULT_ESTIMATORS ) 
        if 'iterations'      in config : config [ 'n_estimators' ] = config.pop ( 'iterations'      , DEFAULT_ESTIMATORS ) 
        if 'max_iter'        in config : config [ 'n_estimators' ] = config.pop ( 'max_iter'        , DEFAULT_ESTIMATORS ) 

        config [ 'normalize'  ] = False                 
        ADVAL_base.__init__ ( self, nToys = nToys , **config ) 
        
    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_RF 
    
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_RF

         
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        
        from sklearn.ensemble import RandomForestRegressor
        from threadpoolctl    import threadpool_limits

        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        w_tr_arr = W_train_mod.values if hasattr(W_train_mod, 'values') else W_train_mod

        params = {}
        params.update ( self.params )

        # ====================================================================
        ## num-boost-rounds & early stopping 
        # ====================================================================
        num_boost_round  = ( params.pop ( 'iterations'      , None ) or
                             params.pop ( 'max_iter'        , None ) or
                             params.pop ( 'num_boost_round' , None ) or
                             params.pop ( 'n_estimators'    , None ) or DEFAULT_ESTIMATORS )         
        
        n_estimators = num_boost_round

        params [ 'n_estimators' ] = n_estimators
        
        with threadpool_limits ( limits = 1 ):
            
            model = RandomForestRegressor ( **params )
            model.fit ( X_train , Y_train_mod , sample_weight = w_tr_arr )            
            predictions = model.predict ( X_val )
            imps = model.feature_importances_ if importance else None
            
        return predictions, imps

# =============================================================================
TORCH_ADAM_DEFAULT_EPOCHS = 200 
TORCH_HIDDEN_DIMS         = 16, 
TORCH_DROPOUT             = 0.05
TORCH_MAX_BATCH_SIZE      = 2 ** 16

# =============================================================================
## @class ADVAL_TORCH_base
#  (Abstract) base class for PyTorch-based Adversarial Validation in regression mode.
# 
#           Provides shared infrastructure for evaluating distribution shifts
#           between two weighted datasets using neural regression on binary labels
#           (with sign-flipping target transformation for negative weights via
#           @c transform_weights_and_targets) and weighted MSE scoring
#           Handles device resolution, thread configuration, model construction,
#           tensor preparation, inference, and explicit memory cleanup, while
#           delegating optimizer-specific training to subclasses via @_fit_fold_model.
class ADVAL_TORCH_base ( ADVAL_base ) :
    """ Base class for PyTorch-based Adversarial Validation (regression mode).
    """

    # =========================================================================
    ## Initialize the PyTorch adversarial validation base estimator.
    #
    #  @param nToys Number of toy permutations used to estimate p-values
    #  @param hidden_dim Hidden layer width (or None for heuristic sizing based on feature count).
    #  @param dropout Dropout probability applied after hidden activation (0 disables dropout).
    #  @param params Additional keyword arguments forwarded to @c ADVAL_base and backends.
    def __init__ ( self        ,
                   nToys = 100 , **params ) :
        """ Initialize the PyTorch adversarial validation base estimator.
        """
        config = {
            'learning_rate' : 1.e-2      ,
            'weight_decay'  : 1.e-4      ,
            'hidden_dims'   : TORCH_HIDDEN_DIMS ,
            'dropout'       : TORCH_DROPOUT     ,
            'n_jobs'        : 1          ,
            'device'        : None       , ## autodetect CUDA
        }
        config.update ( params )    
        ## 
        if 'lr' in config : config  [ 'learning_rate' ] = config.pop ( 'lr' )
        
        import torch as Torch
        
        device = config .get ( 'device' , None )
        if   device is None : device = 'cuda' if Torch.cuda.is_available () else 'cpu'
        elif isinstance ( device, string_types ) :
            device = device.strip().lower()
            device = 'cuda' if 'cuda' == device and Torch.cuda.is_available () else 'cpu'
        else : device = 'cpu'
        
        ## store device 
        self.__device = device 
        config [ 'device'    ] = self.__device 
        
        hidden_dims = config.pop ( 'hidden_dims' , TORCH_HIDDEN_DIMS )
        if     isinstance   ( hidden_dims , integer_types  ) : hidden_dims = hidden_dims ,
        if not isinstance   ( hidden_dims , sequence_types ) :
            raise TypeError ( "Invalid `hidden dims' type: %s" % typename (  hidden_dims ) ) 
        if not hidden_dims : 
            raise TypeError ( "Invalid `hidden dims': %s" % str ( hidden_dims ) )
        if not all ( isinstance ( d , integer_types ) and 1 < d for d in hidden_dims ) :            
            raise TypeError ( "Invalid `hidden dims' type: %s" % typename (  hidden_dims ) ) 
                
        config [ 'hidden_dims' ] = tuple ( int ( d ) for d in hidden_dims )
        self.__hidden_dims = hidden_dims

        dropout = config.get ( 'dropout' , 0 )
        if not 0 <= dropout < 1 : raise TypeError ( "Invalid `dropout': %s" % dropout )
            
        self.__dropout     = dropout   
                
        # Input normalization is enabled by default to stabilize neural training
        # across observables with different scales and units
        config [ 'normalize' ] = True

        super ().__init__ ( nToys = nToys , **config )
        
        n_jobs = self.params.get ( 'n_jobs' , 1 )
        Torch.set_num_threads ( n_jobs if 0 < n_jobs else 1 )
        
    # =========================================================================
    ## Short symbol/tag representing the PyTorch adversarial validator.
    #  @return String identifier symbol used for logging and display.
    @classprop
    def symbol ( self ) :
        """ `symbol`: Short symbol/tag representing the PyTorch adversarial validator.
        """
        return U.symbol_TORCH

    # =========================================================================
    ## Full method identifier name for reporting and tables.
    #  @return Method string tag.
    @property
    def method ( self ) :
        """ `method`: Full method identifier name for reporting and tables.
        """
        return U.method_TORCH

    # =========================================================================
    ## Resolve execution device and configure PyTorch CPU thread count.
    #  @return Normalized @c torch.device instance (CUDA if available, otherwise CPU).
    def get_device ( self ) :
        """ Resolve execution device and configure PyTorch CPU thread count.
        """
        import torch as Torch
        return Torch.device ( self.__device )

    # =========================================================================
    ## Construct a feed-forward regression MLP on target device.
    #
    #           Builds a sequential stack of @c Linear -> @c ReLU -> (optional @c Dropout)
    #           blocks according to @c hidden_dims, followed by a scalar output @c Linear layer.
    #
    #  @param n_features Number of input features.
    #  @param device Target @c torch.device for placing model parameters.
    #  @return Instantiated @c torch.nn.Sequential regression model on @c device.
    def build_model ( self , n_features , device ) :
        """ Construct a feed-forward regression MLP on target device.
        """
        import torch.nn as NN

        hidden_dims  = self.__hidden_dims 
        dropout_rate = self.__dropout     

        layers   = []
        in_dim   = int ( n_features )
        
        for h_dim in hidden_dims :
            
            layers.append ( NN.Linear ( in_dim , h_dim ) )
            layers.append ( NN.ReLU () )
            if 0.0 < dropout_rate : layers.append ( NN.Dropout ( dropout_rate ) )
            in_dim = h_dim

        layers.append ( NN.Linear ( in_dim , 1 ) )
        return NN.Sequential ( *layers ).to ( device )
    
    # =========================================================================
    ## Compute sample-weighted Mean Squared Error loss.
    #  @param preds Predicted continuous values tensor of shape (N, 1).
    #  @param targets Target labels tensor of shape (N, 1).
    #  @param weights Non-negative sample weights tensor of shape (N, 1).
    #  @return Scalar weighted MSE loss tensor.
    @staticmethod
    def weighted_mse ( preds , targets , weights ) :
        """ Compute sample-weighted Mean Squared Error loss.
        """
        import torch.nn.functional as F

        loss = F.mse_loss ( preds , targets , reduction = 'none' )
        return ( loss * weights ).mean ()

    # =========================================================================
    ## Convert feature matrix, targets, and weights into float32 device tensors.
    #  @param X Feature matrix (array-like or DataFrame).
    #  @param Y Target vector/array.
    #  @param W Optional sample weights (array-like, Series, or None).
    #  @param device Target @c torch.device.
    #  @param is_features If True, sanitize NaN/Inf values in feature matrix X.
    #  @return Tuple of contiguous float32 tensors @c (x_t, y_t, w_t) on @c device.
    @staticmethod
    def prepare_xyw ( X , Y , W , device , is_features = False ) :
        """ Convert feature matrix, targets, and weights into float32 device tensors.
        """
        import torch as Torch

        arr_x = numpy.asarray ( X , dtype = numpy.float32 )
        if is_features :
            # Guard against non-finite values before transferring to PyTorch tensors
            arr_x = numpy.nan_to_num ( arr_x , nan = 0.0 , posinf = 0.0 , neginf = 0.0 )

        y_arr = numpy.asarray ( Y , dtype = numpy.float32 ).reshape ( -1 , 1 )
        if W is None :
            w_arr = numpy.ones_like ( y_arr )
        else :
            w_src = W.values if hasattr ( W , 'values' ) else W
            w_arr = numpy.asarray ( w_src , dtype = numpy.float32 ).reshape ( -1 , 1 )

        x_t = Torch.from_numpy ( numpy.ascontiguousarray ( arr_x ) ).to ( device , non_blocking = True )
        y_t = Torch.from_numpy ( numpy.ascontiguousarray ( y_arr ) ).to ( device , non_blocking = True )
        w_t = Torch.from_numpy ( numpy.ascontiguousarray ( w_arr ) ).to ( device , non_blocking = True )

        return x_t , y_t , w_t

    # =========================================================================
    ## Fit the fold model using a concrete optimization strategy.
    #
    #  @param model Neural network model placed on @c device.
    #  @param X_tr_t Training feature tensor on @c device.
    #  @param Y_tr_t Training target tensor on @c device.
    #  @param W_tr_t Training sample weight tensor on @c device.
    #  @param X_va_t Validation feature tensor on @c device.
    #  @param Y_va_t Validation target tensor on @c device.
    #  @param W_va_t Validation sample weight tensor on @c device.
    #  @param device Target @c torch.device.
    @abc.abstractmethod
    def _fit_fold_model ( self   ,
                          model  ,
                          X_tr_t , Y_tr_t , W_tr_t ,
                          X_va_t , Y_va_t , W_va_t , device ) :
        """ Fit the fold model using a concrete optimization strategy.
        """
        raise NotImplementedError

    # =========================================================================
    ## Train on one CV fold and return validation predictions.
    #
    #           Transforms signed weights into non-negative weights with flipped
    #           targets for negative entries, trains @c model via @_fit_fold_model,
    #           evaluates predictions on @c X_val, and explicitly frees tensors
    #           and temporary buffers to keep memory usage low across folds and toys
    #
    #  @param X_train Training features for the current fold.
    #  @param Y_train Training binary labels for the current fold.
    #  @param W_train Training weights for the current fold (or None).
    #  @param X_val Validation features for the current fold.
    #  @param Y_val Validation binary labels for the current fold.
    #  @param W_val Validation weights for the current fold (or None).
    #  @param importance Unused flag kept for uniform @c ADVAL_base interface.
    #  @return Tuple @c (predictions, None) with 1D validation prediction array.
    def work ( self       ,
               X_train    , Y_train    , W_train ,
               X_val      , Y_val      , W_val   , importance = False ) :
        """ Train on one CV fold and return validation predictions.
        """
        import torch as Torch

        Y_tr_mod , W_tr_mod = transform_weights_and_targets ( Y_train , W_train )
        Y_va_mod , W_va_mod = transform_weights_and_targets ( Y_val   , W_val   )

        device     = self.get_device ()
        n_features = num_features ( X_train )

        X_tr_t , Y_tr_t , W_tr_t = self.prepare_xyw ( X_train , Y_tr_mod , W_tr_mod , device , is_features = True )
        X_va_t , Y_va_t , W_va_t = self.prepare_xyw ( X_val   , Y_va_mod , W_va_mod , device , is_features = True )

        model = self.build_model ( n_features , device )
        self._fit_fold_model ( model  ,
                               X_tr_t , Y_tr_t , W_tr_t ,
                               X_va_t , Y_va_t , W_va_t , device )

        model.eval ()
        with Torch.no_grad () :
            predictions = model ( X_va_t ).detach ().cpu ().numpy ().ravel ()

        # Explicitly release fold tensors and model before next fold / toy iteration
        del model , X_tr_t , Y_tr_t , W_tr_t , X_va_t , Y_va_t , W_va_t
        if device.type == 'cuda' : Torch.cuda.empty_cache ()

        return predictions , None

# =============================================================================
## @class ADVAL_TORCH_ADAM
#  Adversarial validation estimator trained with mini-batch AdamW and early stopping.
#
#           Optimizes the regression MLP over epochs using mini-batches from the
#           training fold and tracks validation weighted MSE to restore the best
#           checkpoint while avoiding unnecessary GPU memory retention.
class ADVAL_TORCH_ADAM ( ADVAL_TORCH_base ) :
    """ Adversarial validation estimator trained with mini-batch AdamW and early stopping.
    """

    # =========================================================================
    ## Initialize the AdamW-based adversarial validation estimator.
    #
    #  @param nToys Number of toy permutations used to estimate p-values
    #  @param params Additional keyword arguments (e.g. @c epochs, @c batch_size,
    #         @c early_stopping_rounds, @c learning_rate, @c weight_decay).
    def __init__ ( self            ,
                   nToys    = 100  , **params ) :
        """ Initialize the AdamW-based adversarial validation estimator.
        """
        config = {
            'epochs'                : TORCH_ADAM_DEFAULT_EPOCHS ,
            'batch_size'            : 'auto' ,
            'early_stopping_rounds' : 20     ,
        }
        config.update ( params )
        super ().__init__ ( nToys = nToys , **config )

    # =========================================================================
    ## Short symbol/tag representing the PyTorch adversarial validator.
    #  @return String identifier symbol used for logging and display.
    @classprop
    def symbol ( self ) :
        """ `symbol`: Short symbol/tag representing the PyTorch adversarial validator.
        """
        return U.symbol_TORCH_ADAM

    # =========================================================================
    ## Full method identifier name for reporting and tables.
    #  @return Method string tag.
    @property
    def method ( self ) :
        """ `method`: Full method identifier name for reporting and tables.
        """
        return U.method_TORCH_ADAM

    # =========================================================================
    ## Fit fold model using mini-batch AdamW with validation early stopping.
    # 
    #           Stores best weights as detached CPU tensor clones instead of
    #           deep-copying dictionaries on GPU, then reloads them onto @c device
    #           at the end of training.
    #
    #  @param model Neural network model placed on @c device.
    #  @param X_tr_t Training feature tensor on @c device.
    #  @param Y_tr_t Training target tensor on @c device.
    #  @param W_tr_t Training sample weight tensor on @c device.
    #  @param X_va_t Validation feature tensor on @c device.
    #  @param Y_va_t Validation target tensor on @c device.
    #  @param W_va_t Validation sample weight tensor on @c device.
    #  @param device Target @c torch.device.
    def _fit_fold_model ( self   ,
                          model  ,
                          X_tr_t , Y_tr_t , W_tr_t ,
                          X_va_t , Y_va_t , W_va_t , device ) :
        """ Fit fold model using mini-batch AdamW with validation early stopping.
        """
        import torch as Torch
        from torch.utils.data import DataLoader , TensorDataset

        epochs       = self.params.get ( 'epochs'                , TORCH_ADAM_DEFAULT_EPOCHS ) 
        lr           = self.params.get ( 'learning_rate'         , None   ) or self.params.get ( 'lr' , 1.e-3 ) 
        weight_decay = self.params.get ( 'weight_decay'          , 1.e-3  ) 
        patience     = self.params.get ( 'early_stopping_rounds' , 20     )

        n_rows       = num_samples  ( X_tr_t )
        n_features   = num_features ( X_tr_t )
        
        batch_size   = self.params.get ( 'batch_size'            , 'auto' )
        
        if isinstance ( batch_size, integer_types ) and 2 <= batch_size <= min ( n_rows, TORCH_MAX_BATCH_SIZE ) :
            if batch_size & ( batch_size - 1 ) : batch_size = 2 ** math.floor ( math.log2 ( batch_size ) )
        elif isinstance ( batch_size, string_types ) and batch_size.strip ().lower () == 'full' : batch_size = n_rows
        else : batch_size = optimal_batch_size ( n_rows, n_features, TORCH_MAX_BATCH_SIZE )


        loader    = DataLoader ( TensorDataset ( X_tr_t , Y_tr_t , W_tr_t ) ,
                                 batch_size = batch_size , shuffle = True )
        optimizer = Torch.optim.AdamW ( model.parameters () ,
                                        lr           = lr     ,
                                        weight_decay = weight_decay )
        
        best_val_loss    = float ( 'inf' )
        best_cpu_state   = None
        patience_counter = 0

        for _ in range ( epochs ) :
            
            model.train ()
            for bx , by , bw in loader :
                # set_to_none=True avoids zero-filling gradient buffers every step
                optimizer.zero_grad ( set_to_none = True )
                loss = self.weighted_mse ( model ( bx ) , by , bw )
                loss.backward ()
                optimizer.step ()

            model.eval ()
            with Torch.no_grad () :
                val_loss = float ( self.weighted_mse ( model ( X_va_t ) , Y_va_t , W_va_t ).item () )

            if val_loss < best_val_loss :
                best_val_loss    = val_loss
                # Keep best checkpoint on CPU to save device memory and avoid deepcopy overhead
                best_cpu_state   = { k : v.detach ().cpu ().clone () for k , v in model.state_dict ().items () }
                patience_counter = 0
            else :
                patience_counter += 1
                if 0 < patience and patience <= patience_counter :
                    break

        if best_cpu_state is not None :
            model.load_state_dict ( { k : v.to ( device ) for k , v in best_cpu_state.items () } )
            del best_cpu_state

# =============================================================================
## @class ADVAL_TORCH_LBFGS
#  Adversarial validation estimator trained with quasi-Newton L-BFGS.
#
#           Optimizes the regression MLP using @c torch.optim.LBFGS with strong
#           Wolfe line search over the full training fold by default to preserve
#           sensitivity to subtle or localized differences. Optional random chunk
#           subsampling (@c lbfgs_chunk) can be enabled explicitly when a fold
#           is very large and speed/memory outweigh maximum tail sensitivity.
class ADVAL_TORCH_LBFGS ( ADVAL_TORCH_base ) :
    """ Adversarial validation estimator trained with quasi-Newton L-BFGS.
    """

    # =========================================================================
    ## Initialize the L-BFGS-based adversarial validation estimator.
    # 
    #  @param nToys Number of toy permutations used to estimate p-values
    #  @param params Additional keyword arguments (e.g. @c max_iter, @c lbfgs_chunk,
    #         @c learning_rate, @c weight_decay).
    def __init__ ( self            ,
                   nToys    = 100  , **params ) :
        """ Initialize the L-BFGS-based adversarial validation estimator.
        """
        config = {
            'max_iter'    : 50   ,
            'lbfgs_chunk' : None ,
            'dropout'     : 0    ,  
        }
        
        config.update ( params )
        
        if 'epochs' in config : config  [ 'max_iter' ] = config.pop ( 'epochs' , 25 )
        
        super ().__init__ ( nToys = nToys , **config )

    # =========================================================================
    ## Short symbol/tag representing the PyTorch adversarial validator.
    #  @return String identifier symbol used for logging and display.
    @classprop
    def symbol ( self ) :
        """ `symbol`: Short symbol/tag representing the PyTorch adversarial validator.
        """
        return U.symbol_TORCH_LBFGS

    # =========================================================================
    ## Full method identifier name for reporting and tables.
    #  @return Method string tag.
    @property
    def method ( self ) :
        """ `method`: Full method identifier name for reporting and tables.
        """
        return U.method_TORCH_LBFGS

    # =========================================================================
    ## Fit fold model using L-BFGS on the full fold or an explicit chunk.
    #
    #           Uses the entire training fold by default so subtle distribution
    #           shifts are not diluted by subsampling; if @c lbfgs_chunk is set
    #           to a positive integer smaller than the fold size, draws a random
    #           subset without replacement.
    #
    #  @param model Neural network model placed on @c device.
    #  @param X_tr_t Training feature tensor on @c device.
    #  @param Y_tr_t Training target tensor on @c device.
    #  @param W_tr_t Training sample weight tensor on @c device.
    #  @param X_va_t Validation feature tensor (unused directly in L-BFGS step).
    #  @param Y_va_t Validation target tensor (unused directly in L-BFGS step).
    #  @param W_va_t Validation sample weight tensor (unused directly in L-BFGS step).
    #  @param device Target @c torch.device.
    def _fit_fold_model ( self   ,
                          model  ,
                          X_tr_t , Y_tr_t , W_tr_t ,
                          X_va_t , Y_va_t , W_va_t , device ) :
        """ Fit fold model using L-BFGS on the full fold or an explicit chunk.
        """
        import torch as Torch
        ## 
        max_iter     = self.params.get ( 'max_iter'      , None  ) or self.params.get ( 'epochs' , 25   ) 
        lr           = self.params.get ( 'learning_rate' , None  ) or self.params.get ( 'lr'     , 0.01 ) 
        weight_decay = self.params.get ( 'weight_decay'  , 1.e-4 ) 
        chunk_raw    = self.params.get ( 'lbfgs_chunk'   , None  )
        chunk_size   = int ( chunk_raw ) if chunk_raw is not None else None

        n_tr = X_tr_t.shape [ 0 ]
        if chunk_size is not None and 0 < chunk_size < n_tr :
            # Explicit subsampling only when requested for very large folds
            idx   = Torch.randperm ( n_tr , device = device ) [ : chunk_size ]
            x_opt = X_tr_t.index_select ( 0 , idx )
            y_opt = Y_tr_t.index_select ( 0 , idx )
            w_opt = W_tr_t.index_select ( 0 , idx )
            del idx
        else :
            x_opt , y_opt , w_opt = X_tr_t , Y_tr_t , W_tr_t
            
        optimizer = Torch.optim.LBFGS ( model.parameters () ,
                                        lr             = lr ,
                                        max_iter       = max_iter ,
                                        line_search_fn = 'strong_wolfe' )
        
        def closure () :
            optimizer.zero_grad ( set_to_none = True )
            loss = self.weighted_mse ( model ( x_opt ) , y_opt , w_opt )
            if 0.0 < weight_decay :
                l2_reg = sum ( Torch.sum ( p * p ) for p in model.parameters () )
                loss   = loss + 0.5 * weight_decay * l2_reg
            loss.backward ()
            return loss

        model.train ()
        optimizer.step ( closure )

        if x_opt is not X_tr_t :
            del x_opt , y_opt , w_opt

            
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    from ostap.stats.tools import ( hasLightGBM , hasXGBoost , hasCatBoost ,
                                    hasSkLearn  , hasPyTorch )
    
    if not hasLightGBM ( False ) : logger.warning  ( "No LightGBM available!" ) 
    if not hasXGBoost  ( False ) : logger.warning  ( "No XGBoost  available!" ) 
    if not hasCatBoost ( False ) : logger.warning  ( "No CatBoost available!" ) 
    if not hasSkLearn  ( False ) : logger.warning  ( "No sklearn  available!" ) 
    if not hasPyTorch  ( False ) : logger.warning  ( "No PyTorch  available!" ) 
  
# =============================================================================
##                                                                      The END 
# =============================================================================
    
