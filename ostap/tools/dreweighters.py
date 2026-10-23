#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/dreweghters.py
#  Density ratio reweighters based on LightGBM, XGBoost, CatBoost, SkLearn and HepML
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @Date  2026-07-18
# =============================================================================
""" Density ratio reweighters based on
 - LightGBM
 - XGBoost
 - CatBoost
 - SkLearn 
 - HepML.
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    'LightGBMDensityReweighter'      , # LightGBM-based density reweighter
    'XGBoostDensityReweighter'       , # XGBoost-based density reweighter
    'CatBoostDensityReweighter'      , # CatBoost-based density reweighter
    'PyTorchDensityReweighter'       , # PyTorch-based density reweighter
    'LogRegressionDensityReweighter' , # Density reweighter based on Logistic Regression
    'GBReweighter'                   , # Reweighter based on GBReweighter from hep_ml    
) 
# =============================================================================
from   ostap.core.ostap_types   import num_types, string_types, integer_types  
from   ostap.utils.core         import typename
from   ostap.utils.basic        import numcpu, NoContext
from   ostap.logger.utils       import map2table_ex
from   ostap.logger.pretty      import nice_print 
from   ostap.utils.progress_bar import progress_bar, text_bar
from   ostap.tools.reweighter   import Reweighter
from   ostap.stats.utils        import ( weight_trivial     ,
                                         valid_weight       ,
                                         valid_data_shape   ,
                                         compatible_weights , 
                                         num_features       ,
                                         num_samples        ,
                                         check_all          ,
                                         nEff               )
from   ostap.stats.tools            import hasSkLearn 
from   ostap.tools.base_reweighters import DensityReweighter  
import ostap.logger.symbols     as     S
import numpy, abc, math, gc  
# =============================================================================
# Logging setup
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.dreweighters' )
else                       : logger = getLogger( __name__ )
# =============================================================================
has_sklearn = hasSkLearn() 
# =============================================================================
## Decorations 
# =============================================================================
method_LGBM  = 'DRW/%s'   % ( S.light_bulb if S.light_bulb  else 'LightGBM' ) 
method_XGB   = 'DRW/%s'   % ( S.rocket     if S.rocket      else 'XGBoost'  ) 
method_CATB  = 'DRW/%s'   % ( S.cat_face   if S.cat_face    else 'CatBoost' ) 
method_TORCH = 'DRW/%s'   % ( S.flashlight if S.flashlight  else 'TORCH'    ) 
method_LR    = 'DRW/%s'   % ( S.ruler      if S.ruler       else 'LOGREG'   ) 
method_GBRW  = 'HepML/%s' % ( S.wood       if S.wood        else 'GBRW'     )
#
epoch_symbol = ( '%s :' % S.repeat ) if S.show else 'Epoch:'
# =============================================================================
# Global Configuration Constants
# =============================================================================
DEFAULT_ESTIMATORS             = 1000  # High capacity: deep ensemble with low learning rate
REGULARIZED_ESTIMATORS         =  500  # Strong regularization: constrained number of trees
MAX_DEPTH                      =   10  # Deep tree capacity for high dimensions / large statistics
REGULARIZED_DEPTH              =    5  # Shallow depth (8 leaves max) to enforce heavy smoothing
LEARNING_RATE                  =  0.02 # Small step for stable density ratio convergence
REGULARIZED_LEARNING_RATE      =  0.04 # Slightly larger step for shallow regularized trees

LGBM_DEFAULT_ESTIMATORS        = 1600 # High capacity: deep ensemble with low learning rate
LGBM_MAX_DEPTH                 =   10 # Deep tree capacity for high dimensions / large statistics
LGBM_LEARNING_RATE             = 0.01 # Small step for stable density ratio convergence

XGB_DEFAULT_ESTIMATORS         = 1600 # High capacity: deep ensemble with low learning rate
XGB_MAX_DEPTH                  =   10 # Deep tree capacity for high dimensions / large statistics
XGB_LEARNING_RATE              = 0.01 # Small step for stable density ratio convergence

CATB_DEFAULT_ESTIMATORS        =  1600 # High capacity: deep ensemble with low learning rate
CATB_MAX_DEPTH                 =    10 # Deep tree capacity for high dimensions / large statistics
CATB_LEARNING_RATE             = 0.005 # Small step for stable density ratio convergence

# =============================================================================
## @class  LightGBMDensityReweighter
#  Density ratio reweighter using LightGBM as the underlying classifier.
class LightGBMDensityReweighter ( DensityReweighter ) :
    """ Density ratio reweighter using LightGBM as the underlying classifier.
    """
    # =========================================================================
    ## Initialize LightGBM density reweighter.
    #  @param original Features array for original sample.
    #  @param target Features array for target sample.
    #  @param original_weight Initial weights for original sample (optional).
    #  @param target_weight Initial weights for target sample (optional).
    #  @param kwargs Additional LightGBM parameters.
    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        """ Initialize LightGBM density reweighter (High Capacity Mode Defaults).
        """
        config = {
            'objective'             : 'binary'                ,
            'metric'                : 'binary_logloss'        ,
            'n_estimators'          : LGBM_DEFAULT_ESTIMATORS ,
            'learning_rate'         : LGBM_LEARNING_RATE      ,
            'max_depth'             : LGBM_MAX_DEPTH          ,
            'num_leaves'            : min ( 1023 , 2 ** LGBM_MAX_DEPTH - 1 ) , 
            'max_bin'               : 2047                    , # 
            'boost_from_average'    : True                    ,
            'early_stopping_rounds' : LGBM_DEFAULT_ESTIMATORS // 2 ,
            'verbosity'             : -1                      ,
            'n_jobs'                : -1                      ,
        }
        
        config.update ( kwargs )

        if 'silent' in config and config.get ( 'silent' ) : config [ 'verbose' ] = -1 

        if 'num_boost_round' in config :
            config [ 'n_estimators' ] = config.pop ( 'num_boost_round' , LGBM_DEFAULT_ESTIMATORS ) 
            
        super ().__init__ (
            original        = original        ,
            target          = target          ,
            original_weight = original_weight ,
            target_weight   = target_weight   ,
            **config
        )

    # =========================================================================
    ## Return the method identifier name.
    #  @return Method string identifier.
    @property
    def method ( self ) :
        """ Return the method identifier name.
        """
        return method_LGBM 

    # =========================================================================
    ## Train single LightGBM model on fold data.
    #  @param X_train Training features array.
    #  @param y_train Training binary labels.
    #  @param w_train Training event weights or None.
    #  @param X_val Validation features array.
    #  @param y_val Validation binary labels.
    #  @param w_val Validation event weights or None.
    #  @return Tuple of (trained_lightgbm_booster, val_predictions).
    def _train_single_model ( self    ,
                              X_train , y_train , w_train ,
                              X_val   , y_val   , w_val   ) : 
        """ Train single LightGBM model on fold data.
        """
        import lightgbm as LightGBM
        import gc
        
        X_tr     = numpy.ascontiguousarray ( X_train , dtype = numpy.float32 )
        X_va     = numpy.ascontiguousarray ( X_val   , dtype = numpy.float32 )
        w_tr     = numpy.ascontiguousarray ( w_train , dtype = numpy.float32 ) if w_train is not None else None
        w_va     = numpy.ascontiguousarray ( w_val   , dtype = numpy.float32 ) if w_val   is not None else None

        trn_data = LightGBM.Dataset    ( X_tr    , label = y_train , weight = w_tr ,                        free_raw_data = True )
        val_data = LightGBM.Dataset    ( X_va    , label = y_val   , weight = w_va , reference = trn_data , free_raw_data = True )

        params = self.params.copy ()
        num_boost_round          = params.pop ( 'num_boost_round'       , None ) or params.pop ( 'n_estimators' , LGBM_DEFAULT_ESTIMATORS )
        early_stopping_rounds    = params.pop ( 'early_stopping_rounds' , None )
        params [ 'num_threads' ] = params.get ( 'num_threads'           , params.pop ( 'n_jobs' , max ( 1 , numcpu () // 2 ) ) ) 

        
        if early_stopping_rounds :
            callbacks  = [ LightGBM.early_stopping ( stopping_rounds = early_stopping_rounds , verbose = False ) ]
            valid_sets = [ val_data ]
        else : 
            callbacks  = []
            valid_sets = None 

        model = LightGBM.train ( params          = params          ,
                                 train_set       = trn_data        ,
                                 num_boost_round = num_boost_round ,
                                 valid_sets      = valid_sets      ,
                                 callbacks       = callbacks       )
            
        # =====================================================================
        ## PHOENIX
        # =====================================================================
        if True : # ===========================================================
            # =================================================================
            model_bytes = model.model_to_string()
            best_iter   = getattr ( model , 'best_iteration_' , None ) or getattr ( model , 'best_iteration' , None )                    
            model.free_dataset() 
            del model
            del trn_data
            del val_data
            del X_tr, X_va, w_tr, w_va
            
            model = LightGBM.Booster  ( model_str = model_bytes  )
            if not best_iter is None : model.best_iteration = best_iter 
            del model_bytes
            
            gc.collect ()                

            
        val_preds = self._predict_single_model ( model , X_val )

        return model , val_preds

    # =========================================================================
    ## Predict probabilities using a LightGBM model.
    #  @param model Trained LightGBM booster instance.
    #  @param X Features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model( self, model, X ):
        """ Predict probabilities using a LightGBM model.
        """
        best_iter  = getattr ( model , 'best_iteration_' , None ) or getattr ( model , 'best_iteration' , None )        
        kwargs     = {}
        if best_iter is not None and 0 < best_iter : kwargs[ 'num_iteration' ] = best_iter            
        num_threads = self.params.get ( 'num_threads' , self.params.get ( 'n_jobs' , max ( 1 , numcpu () // 2 ) ) )                                    
        p = model.predict ( X , num_threads = num_threads , **kwargs )
        return p.astype( numpy.float32, copy = False )

# =============================================================================
## @class  XGBoostDensityReweighter
#  Density ratio reweighter using XGBoost as the underlying classifier.
class XGBoostDensityReweighter ( DensityReweighter ) :
    """ Density ratio reweighter using XGBoost as the underlying classifier.
    """

    # =========================================================================
    ## Initialize XGBoost density reweighter.
    #  @param original Features array for original sample.
    #  @param target Features array for target sample.
    #  @param original_weight Initial weights for original sample (optional).
    #  @param target_weight Initial weights for target sample (optional).
    #  @param kwargs Additional XGBoost parameters.
    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        """ Initialize XGBoost density reweighter (High Capacity Mode Defaults).
        """
        config = {
            'objective'             : 'binary:logistic'   ,
            'eval_metric'           : 'logloss'           ,
            'n_estimators'          : XGB_DEFAULT_ESTIMATORS  ,
            'learning_rate'         : XGB_LEARNING_RATE       ,
            'max_depth'             : XGB_MAX_DEPTH           ,
            'tree_method'           : 'hist'                  , # Fast histogram-based algorithm (similar to LGBM)
            'max_bin'               : 1023                    , # Limit binning to prevent noise fitting
            'min_child_weight'      : 1.e-3                   , # Minimum sum of instance weight (hessian) in a child
            'early_stopping_rounds' : XGB_DEFAULT_ESTIMATORS // 2 ,            
            'verbosity'             :  0                      ,      
            'n_jobs'                : -1                      ,
        }
        config.update ( kwargs )
        
        if 'silent' in config and config.get ( 'silent' ) : config [ 'verbosity' ] = 0 
        
        super ().__init__ (
            original        = original        ,
            target          = target          ,
            original_weight = original_weight ,
            target_weight   = target_weight   ,
            **config
        )

    # =========================================================================
    ## Return the method identifier name.
    #  @return Method string identifier.
    @property
    def method ( self ) :
        """ Return the method identifier name.
        """
        return method_XGB 


    # =========================================================================
    ## Train single XGBoost booster on fold data.
    #  @param X_train Training features array.
    #  @param y_train Training binary labels.
    #  @param w_train Training event weights or None.
    #  @param X_val Validation features array.
    #  @param y_val Validation binary labels.
    #  @param w_val Validation event weights or None.
    #  @return Tuple of (fitted_xgb_booster, val_predictions).
    def _train_single_model ( self,
                              X_train , y_train , w_train ,
                              X_val   , y_val   , w_val   ) :
        """ Train single XGBoost booster on fold data.
        """
        import xgboost as XGBoost
        
        # Ensure contiguous C-style memory layout and explicit float32 data types
        X_tr = numpy.ascontiguousarray ( X_train , dtype = numpy.float32 )
        y_tr = numpy.ascontiguousarray ( y_train , dtype = numpy.float32 )
        w_tr = numpy.ascontiguousarray ( w_train , dtype = numpy.float32 ) if w_train is not None else None

        X_v  = numpy.ascontiguousarray ( X_val   , dtype = numpy.float32 )
        y_v  = numpy.ascontiguousarray ( y_val   , dtype = numpy.float32 )
        w_v  = numpy.ascontiguousarray ( w_val   , dtype = numpy.float32 ) if w_val is not None else None

        # Construct C++ DMatrix objects for training and validation
        dtrain = XGBoost.DMatrix ( X_tr , label = y_tr , weight = w_tr )
        dval   = XGBoost.DMatrix ( X_v  , label = y_v  , weight = w_v  )
        
        params = self.params.copy()
        params [ 'nthread' ] = params.get ( 'nthread' , params.pop ( 'n_jobs' , max ( 1 , numcpu () // 2 ) ) ) 
        
        n_estimators          = params.pop ( 'n_estimators'           , 1000 )
        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , None )

        evals = [ ( dval , 'val' ) ] if early_stopping_rounds is not None else []
        
        # Train model via low-level C++ API
        model = XGBoost.train (
            params,
            dtrain,
            num_boost_round      = n_estimators,
            evals                = evals,
            early_stopping_rounds= early_stopping_rounds,
            verbose_eval         = False
        )

        # =====================================================================
        ## PHOENIX
        # =====================================================================
        if True : # ===========================================================
            # =================================================================
            
            # Export tree structure into a lightweight binary/JSON byte array
            model_bytes = model.save_raw ( raw_format = "json" )
            best_iter   = getattr ( model , 'best_iteration_' , None ) or getattr ( model , 'best_iteration' , None )                    
            
            # Destroy heavy training DMatrix objects and training C++ booster instance
            del dtrain, dval, model,
            del X_tr, y_tr , w_tr, X_v , y_v, w_v 
            
            # Re-instantiate clean C++ booster dedicated strictly to inference
            model = XGBoost.Booster()
            model.load_model ( model_bytes )
            if not best_iter is None : model.best_iteration = best_iter
            
            del model_bytes            
            gc.collect()

        # Compute validation predictions safely using pre-allocated contiguous array
        val_preds = self._predict_single_model ( model , X_val )

        return model, val_preds
    
    # =========================================================================
    ## Predict probabilities using an XGBoost model.
    #  @param model Trained XGBoost Booster.
    #  @param X Input features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model( self , model, X ):
        """ Predict probabilities using an XGBoost model.
        """
        import xgboost as XGBoost
        ## 
        dmat      = XGBoost.DMatrix ( X )
        best_iter = getattr ( model , 'best_iteration_' ,  None ) or getattr ( model , 'best_iteration' , None )
        kwargs    = {}
        if best_iter is not None and 0 < best_iter : kwargs [ 'iteration_range' ] = ( 0, best_iter + 1 )
        p         = model.predict ( dmat , **kwargs )
        return p.astype ( numpy.float32, copy = False )

# =============================================================================
## @class CatBoostDensityReweighter
#  Density ratio reweighter using CatBoost as the underlying classifier.
#
# Tested/tuned configurations for "small" test"
# - 
class CatBoostDensityReweighter(DensityReweighter):
    """ Density ratio reweighter using CatBoost as the underlying classifier.
    """
    
    # =========================================================================
    def __init__ ( self                   , * , 
                   original               ,
                   target                 ,
                   original_weight        = None ,
                   target_weight          = None ,
                   store_original_weights = True , **params ) :
        """ Initialize CatBoost density reweighter.
        """
        
        from catboost.utils import get_gpu_device_count
        gpu = get_gpu_device_count() 
        self.__gpu = gpu 
        
        config = { 'loss_function'         : 'Logloss'                        ,
                   'eval_metric'           : 'Logloss'                        ,
                   'iterations'            : CATB_DEFAULT_ESTIMATORS          ,
                   'learning_rate'         : CATB_LEARNING_RATE               ,
                   'depth'                 : CATB_MAX_DEPTH                   ,
                   'l2_leaf_reg'           : 0.5                              ,
                   'min_data_in_leaf'      : 5                                ,
                   'bootstrap_type'        : 'MVS'                            ,
                   'mvs_reg'               : 2.0                              , 
                   'random_strength'       : 0.01                             ,
                   'early_stopping_rounds' : CATB_DEFAULT_ESTIMATORS // 2     ,
                   'boosting_type'         : 'Plain'                          ,
                   'device'                : 'GPU' if 0 < self.gpu else 'CPU' , 
                   'verbose'               : False                            }
        
        config.update ( params )
        
        # Alias mappings to native CatBoost parameters
        if 'min_child_samples' in config : config [ 'min_data_in_leaf' ] = config.pop ( 'min_child_samples' )
        if 'n_jobs'            in config : config [ 'thread_count'     ] = config.pop ( 'n_jobs'            )
        if 'n_estimators'      in config : config [ 'iterations'       ] = config.pop ( 'n_estimators'      )
        
        ## if 'silent' in config :
        ##    silent = config.get ( 'silent' , True ) 
        ## config [ 'verbose' ] = not silent 

        device = config.get ( 'device' , 'GPU' if 0 < self.gpu else 'CPU' )
        if self.gpu <= 0 : config [ 'device' ] = 'CPU'
        
        super().__init__( original               = original               ,
                          target                 = target                 ,
                          original_weight        = original_weight        ,
                          target_weight          = target_weight          ,
                          store_original_weights = store_original_weights , **config ) 
        
    # =========================================================================
    @property
    def gpu ( self ) :
        """`gpu` : result of `catboost.utils.get_gpu+_device_count"""
        return self.__gpu
    
    # =========================================================================
    @property
    def method(self):
        return method_CATB
    
    # =========================================================================
    def _train_single_model ( self    ,
                              X_train , y_train , w_train ,
                              X_val   , y_val   , w_val   ) :
        """ Train single CatBoost model on memory fold data without forced contiguous copies.
        """

        import gc
        import catboost as CatBoost
        
        # Ensure contiguous C-style memory layout and explicit data types
        X_tr = numpy.ascontiguousarray ( X_train , dtype = numpy.float32  )
        y_tr = numpy.ascontiguousarray ( y_train , dtype = numpy.int32    )
        w_tr = numpy.ascontiguousarray ( w_train , dtype = numpy.float32  ) if w_train is not None else None
    
        X_v  = numpy.ascontiguousarray ( X_val   ,   dtype=numpy.float32  )
        y_v  = numpy.ascontiguousarray ( y_val   ,   dtype=numpy.int32    )
        w_v  = numpy.ascontiguousarray ( w_val   ,   dtype=numpy.float32  ) if w_val is not None else None
        
        # Instantiate heavy C++ Data Pools
        trn_pool = CatBoost.Pool ( X_tr , label = y_tr , weight = w_tr )
        val_pool = CatBoost.Pool ( X_v  , label = y_v  , weight = w_v  )
        
        params = self.params.copy()
        
        params [ 'thread_count'  ] = params.get ( 'thread_count' , params.pop ( 'n_jobs' , max ( 1 , numcpu () // 2 ) ) ) 
        params [ 'boosting_type' ] = params.get ( 'boosting_type' , 'Plain')
        
        iterations = params.pop ('n_estimators' , params.pop ( 'iterations' , CATB_DEFAULT_ESTIMATORS ) ) 
        if not isinstance(iterations, int) or not 10 <= iterations <= 10000:
            iterations = CATB_DEFAULT_ESTIMATORS
        
        early_stopping_rounds = params.pop('early_stopping_rounds', None)
        if isinstance(early_stopping_rounds, int) and 1 < early_stopping_rounds < iterations : pass
        else    : early_stopping_rounds = None
        
        params [ 'iterations' ] = iterations
        
        fit_kwargs = {}
        if early_stopping_rounds is not None:
            params     [ 'early_stopping_rounds' ] = early_stopping_rounds
            params     [ 'use_best_model'        ] = True
            fit_kwargs [ 'eval_set'              ] = val_pool
            
        params [ 'task_type' ] = params.pop ( 'device', 'GPU' if 0 < self.gpu else 'CPU' )
            
        model = CatBoost.CatBoostClassifier(**params)
        
        verbose_level = params.get ('verbose' , False)
        model.fit ( trn_pool , verbose = verbose_level, **fit_kwargs)

        # =====================================================================
        ## PHOENIX
        # =====================================================================
        if True : # ===========================================================
            # =================================================================
            import pickle 
                        
            # Extract tree geometry into a lightweight binary blob
            model_bytes = pickle.dumps ( model ) 
            best_iter   = getattr ( model , 'best_iteration_' , None ) or getattr ( model , 'best_iteration' , None )                    
            
            # Destroy training pools and heavy training model instance
            del trn_pool, val_pool, model, X_tr, y_tr, w_tr, X_v, y_v, w_v
            
            # Re-instantiate clean, standalone inference engine
            model = pickle.loads ( model_bytes )
            if not best_iter is None : model.best_iteration = best_iter
            
            del model_bytes 

            gc.collect()
            
        # Compute validation predictions using pre-allocated contiguous array
        val_preds = self._predict_single_model ( model,  X_val )
        
        return model, val_preds

    # =========================================================================
    def _predict_single_model ( self , model, X ) :
        """ Predict probabilities using a CatBoost model.
        """
        X_clean   = numpy.ascontiguousarray ( X , dtype = numpy.float32  )
        best_iter = getattr ( model , 'best_iteration_' , None ) or getattr ( model , 'best_iteration' , None )        
        kwargs    = {}
        if best_iter is not None and 0 <= best_iter : kwargs [ 'ntree_end' ] = best_iter + 1
        thread_count = self.params.get ( 'thread_count' , self.params.get ( 'n_jobs' , max ( 1 , numcpu () // 2 ) ) )                                    
        p = model.predict_proba ( X_clean , thread_count = thread_count ,  **kwargs)[:, 1]
        
        return p.astype (numpy.float32 , copy = False )


# =============================================================================
## @class PyTorchDensityReweighter
#  Optimized PyTorch MLP density ratio reweighter.
#  
#  This implementation utilizes index-based zero-copy batching, Automatic 
#  Mixed Precision (AMP), Exponential Moving Average (EMA) for weights, 
#  strict gradient clipping, and bounded direct density ratio loss for 
#  maximum training stability and hardware utilization.
# =============================================================================
class PyTorchDensityReweighter ( DensityReweighter ) :
    """ Highly optimized and stabilized PyTorch MLP density ratio reweighter.
    """

    # =========================================================================
    ## Constructor.
    #
    #  @param original        Original dataset.
    #  @param target          Target dataset.
    #  @param original_weight Optional weights for the original dataset.
    #  @param target_weight   Optional weights for the target dataset.
    #  @param progress        Boolean flag to show epoch progress bar.
    #  @param kwargs          Additional model hyperparameters.
    # =========================================================================
    def __init__ ( self                   , * ,
                   original               ,
                   target                 ,
                   original_weight = None ,
                   target_weight   = None ,
                   progress        = True , **kwargs ) :
                   
        import torch
        cuda = torch.cuda.is_available () 

        # ---------------------------------------------------------------------
        # Default configuration aligned for readability
        # ---------------------------------------------------------------------
        config = { 'hidden_dims'   : ( 256 , 256 , 256 )       ,
                   'dropout'       : 0.35                      ,               
                   'learning_rate' : 2e-4                      ,              
                   'weight_decay'  : 1e-5                      ,              
                   'batch_size'    : 256                       ,
                   'epochs'        : 2000                      ,         
                   'patience'      :   50                      ,
                   'min_epochs'    :   50                      ,
                   'min_delta'     : 1e-5                      , 
                   'eval_freq'     : 5                         , 
                   'max_logit'     : 8.0                       , # Cap w = exp(z) <= ~400
                   'logit_reg'     : 1e-4                      , # L2 penalty on logits
                   'grad_clip'     : 0.2                       ,
                   'ema_decay'     : 0.999                     , 
                   'ma_window'     : 5                         , 
                   'device'        : 'cuda' if cuda else 'cpu' ,
                   'compile'       : False                     ,
                   'jit'           : True                      ,
                   'n_jobs'        : -1                        }

        config.update ( kwargs )
        self.__progress_epochs = config.pop ( 'progress' , True )
        
        device = config.get ( 'device' , 'cuda' if cuda else 'cpu' )        
        if not cuda : config [ 'device' ] = 'cpu'

        batch_size = config.get ( 'batch_size' , None )

        n = num_samples ( target ) + num_samples ( original )
        
        if   isinstance ( batch_size , string_types  ) and batch_size.lower() == 'auto' : batch_size = None
        elif isinstance ( batch_size , integer_types ) and not 1 <= batch_size  <= n    : batch_size = None
        
        if batch_size is None : 
            batch_size = math.floor ( math.log2 ( batch_size ) )
            batch_size = 2 ** batch_size  
            config [ 'batch_size' ] = batch_size 
            
        super().__init__ ( original        = original        ,
                           target          = target          ,
                           original_weight = original_weight ,
                           target_weight   = target_weight   ,
                           progress        = False           , **config )

        torch.set_num_threads ( 1  )

    # =========================================================================
    ## Returns the backend method identifier.
    # =========================================================================
    @property
    def method ( self ) :
        return method_TORCH 

    # =========================================================================
    ## Trains a single neural network model using Tabular ResNet architecture.
    #
    #  @param X_train Features for training.
    #  @param y_train Labels for training.
    #  @param w_train Sample weights for training.
    #  @param X_val   Features for validation.
    #  @param y_val   Labels for validation.
    #  @param w_val   Sample weights for validation.
    #  @return Tuple of ( trained_model , validation_predictions ).
    # =========================================================================
    def _train_single_model ( self , X_train , y_train , w_train , X_val , y_val , w_val ) :

        import torch
        
        n_threads = self.params.get ( 'n_jobs' , max ( 1 , numcpu () // 2 ) )
        n_threads = min ( 4  , n_threads ,       max ( 1 , numcpu () // 2 ) )
        n_threads = max ( 1  , n_threads )
        torch.set_num_threads ( 4 )

        import torch.nn                        as nn
        import torch.nn.functional             as F
        from   torch.optim.swa_utils           import AveragedModel , get_ema_multi_avg_fn
        from   sklearn.preprocessing           import QuantileTransformer, StandardScaler 
        from   collections                     import deque
        import time

        device  = torch.device ( self.params.get ( 'device' , 'cuda' if torch.cuda.is_available () else 'cpu' ) )
        is_cuda = device.type == 'cuda'
        
        n_samples  = num_samples  ( X_train )
        n_features = num_features ( X_train )
        
        # =====================================================================
        # Extract and validate hyperparameters
        # =====================================================================
        batch_size = self.params.get ( 'batch_size' , 256 )
        if not isinstance ( batch_size , int ) or batch_size <= 1 : batch_size = 256
            
        batch_size = min ( max ( 2 , ( n_samples + 1 ) // 2 ) , batch_size ) 
        batch_size = max ( 2 , 2 ** math.floor ( math.log2 ( batch_size ) ) )

        nepochs       = self.params.get ( 'epochs'        , 3000          )
        patience      = self.params.get ( 'patience'      , 300           )
        min_epochs    = self.params.get ( 'min_epochs'    , 100           )
        min_delta     = self.params.get ( 'min_delta'     , 1.e-4         ) 
        eval_freq     = self.params.get ( 'eval_freq'     , 5             )
        max_logit     = self.params.get ( 'max_logit'     , 8.0           )
        logit_reg     = self.params.get ( 'logit_reg'     , 1e-3          )
        grad_clip     = self.params.get ( 'grad_clip'     , 0.2           )
        ema_decay     = self.params.get ( 'ema_decay'     , 0.9           )
        ma_window     = self.params.get ( 'ma_window'     , 5             )
        learning_rate = self.params.get ( 'learning_rate' , 3e-4          )
        dropout       = self.params.get ( 'dropout'       , 0.1           )
        hidden_dims   = self.params.get ( 'hidden_dims'   , ( 256 , 256 ) )
        weight_decay  = self.params.get ( 'weight_decay'  , 1.e-5         )
        bce_epochs    = self.params.get ( 'bce_epochs'    , 0             )
        use_compile   = self.params.get ( 'compile'       , False         ) and hasattr ( torch , 'compile' )
        use_jit       = self.params.get ( 'jit'           , False         ) 
        
        # =====================================================================
        # Data preprocessing
        # =====================================================================

        ## n_quantiles = min ( 10000 , n_samples ) 
        ## scaler      = QuantileTransformer  ( output_distribution = 'normal' , n_quantiles = n_quantiles , random_state = self.random_state )
        scaler      = StandardScaler () 

        X_tr_scaled = scaler.fit_transform ( numpy.ascontiguousarray ( X_train , dtype = numpy.float32 ) ) .astype ( numpy.float32 , copy = False)
        X_va_scaled = scaler.transform     ( numpy.ascontiguousarray ( X_val   , dtype = numpy.float32 ) ) .astype ( numpy.float32 , copy = False)

        w_tr = numpy.ascontiguousarray ( w_train , dtype = numpy.float32 ) if w_train is not None else numpy.ones ( n_samples   , dtype = numpy.float32 )
        w_va = numpy.ascontiguousarray ( w_val   , dtype = numpy.float32 ) if w_val   is not None else numpy.ones ( len ( y_val ) , dtype = numpy.float32 )
        y_tr = y_train.astype ( numpy.float32 )

        # =====================================================================
        # Move data directly to device tensors
        # =====================================================================
        X_tr_t = torch.from_numpy ( X_tr_scaled ).to ( device )
        y_tr_t = torch.from_numpy ( y_tr        ).to ( device )
        w_tr_t = torch.from_numpy ( w_tr        ).to ( device )
        
        X_va_t = torch.from_numpy ( X_va_scaled ).to ( device )
        y_va_t = torch.from_numpy ( y_val.astype ( numpy.float32 ) ).to ( device )
        w_va_t = torch.from_numpy ( w_va        ).to ( device )

        # =====================================================================
        ## @class ResBlock
        
        class DensityRatioLoss(nn.Module):
            def __init__(self, max_logit=12.0, logit_reg=0.0):
                super().__init__()
                self.max_logit = max_logit
                self.logit_reg = logit_reg

            def forward(self, logits, targets, weights):
                logits  = logits.float()
                targets = targets.float()
                weights = weights.float()
                
                mask_orig = (targets == 0)
                mask_targ = (targets == 1)
                
                z_orig = logits[mask_orig]
                w_orig = weights[mask_orig]
                z_targ = logits[mask_targ]
                w_targ = weights[mask_targ]
                
                w_orig_sum = w_orig.abs().sum() + 1e-8
                w_targ_sum = w_targ.abs().sum() + 1e-8
                
                # Ограничиваем ТОЛЬКО вход в exp(), чтобы избегать inf во float32
                exp_z = torch.exp(torch.clamp(z_orig, max=self.max_logit))
                
                loss_orig = torch.sum(w_orig * exp_z) / w_orig_sum
                loss_targ = torch.sum(w_targ * z_targ) / w_targ_sum
                
                loss = loss_orig - loss_targ
                
                if self.logit_reg > 0:
                    loss += self.logit_reg * torch.mean(logits ** 2)
                    
                return loss
                
        
        # =====================================================================
        ## @class ResBlock
        class ResBlock ( nn.Module ) :
            def __init__ ( self , dim , dropout = 0.1 ) :
                super().__init__()

                ## Build layers 
                layers = [ nn.Linear      ( dim , dim ) ,
                           nn.LayerNorm   ( dim )       ,
                           ## nn.BatchNorm1d ( dim )       ,
                           nn.SiLU        ()            ] 
                
                if 0 < dropout < 1 : layers.append ( nn.Dropout ( p = dropout ) ) 
                                    
                layers.append ( nn.Linear    ( dim , dim ) )
                layers.append ( nn.LayerNorm ( dim )       ) 
                
                self.block = nn.Sequential ( *layers ) 
                self.act   = nn.SiLU ()

            def forward ( self , x ) :
                return self.act ( x + self.block ( x ) )

        # =====================================================================
        ## @class TabularResNet
        class TabularResNet ( nn.Module ) :
            def __init__ ( self , n_features , hidden_dims = ( 256 , 256 ) , dropout = 0.1 ) :
                super().__init__()
                layers = []
                
                first_dim = hidden_dims [ 0 ]
                layers.append ( nn.Linear ( n_features , first_dim ) )
                
                in_dim = first_dim
                for h_dim in hidden_dims :
                    if in_dim != h_dim :
                        layers.append ( nn.Linear ( in_dim , h_dim ) )
                        in_dim = h_dim
                    layers.append ( ResBlock ( h_dim , dropout = dropout ) )

                if 0.0 < dropout < 1.0 : layers.append ( nn.Dropout ( p = dropout ) )
            
                layers.append ( nn.Linear ( in_dim , 1 ) )
                self.net = nn.Sequential ( *layers )
                
            def forward ( self , x ) : 
                return self.net ( x ).squeeze ( -1 )

        model_config = { 'n_features'  : n_features  ,
                         'hidden_dims' : hidden_dims ,
                         'dropout'     : dropout     }
                         
        model     = TabularResNet ( **model_config ).to ( device )
        ema_model = AveragedModel ( model , multi_avg_fn = get_ema_multi_avg_fn ( ema_decay ) )

        # ======================================================================================
        if False and use_compile : # ===========================================================
            # ==================================================================================
            try : # ============================================================================
                # ==============================================================================
                model     = torch.compile ( model     )
                ema_model = torch.compile ( ema_model )
                # ==============================================================================
            except Exception as e : # ==========================================================
                # ==============================================================================
                logger.warning ( f"torch.compile failed, falling back to eager mode: {e}" )
            # ==================================================================================
        elif use_jit : # =======================================================================
            # ==================================================================================
            try : # ============================================================================
                # ==============================================================================
                model     = torch.jit.script ( model     )
                # ==============================================================================
            except Exception as e : # ==========================================================
                # ==============================================================================
                logger.warning ( f"torch.jit.script failed, falling back to eager mode: {e}" )
            # ==================================================================================
            
        optimizer     = torch.optim.AdamW                          ( model.parameters () , lr      = learning_rate , weight_decay =       weight_decay  )
        scheduler     = torch.optim.lr_scheduler.CosineAnnealingLR ( optimizer           , T_max   = nepochs       , eta_min      = 0.1 * learning_rate )
        scaler_amp    = torch.amp.GradScaler                       ( 'cuda'              , enabled = is_cuda )
        criterion_drw = DensityRatioLoss                           ( max_logit = max_logit , logit_reg = logit_reg )

        w_va_sum         = w_va_t.abs().sum() + 1e-8
        best_loss        = float ( 'inf' )
        best_epoch       = 0 
        best_state       = None
        patience_counter = 0

        # Moving average window for smooth validation monitoring
        val_loss_history = deque ( maxlen = ma_window )

        if not self.silent :
            data_size = n_features * n_samples * 4
            if   1024 ** 3 < data_size : data_size = '%.0fGB' % ( data_size / 1024 ** 3 )
            elif 1024 ** 2 < data_size : data_size = '%.0fMB' % ( data_size / 1024 ** 2 )
            elif 1024      < data_size : data_size = '%.0fkB' % ( data_size / 1024      )
            else                       : data_size = '%dB'    % ( data_size             )                                    
            steps_epoch = math.ceil  ( n_samples  / batch_size )
            logger.info ( "Train single model #events=%d total_data=%s batch_size=%d steps/epoch=%.0f" % ( n_samples   ,
                                                                                                           data_size   ,
                                                                                                           batch_size  ,
                                                                                                           steps_epoch ) )
            
        def get_autocast_config ( device ) :
            if   'cuda' == device.type : return True , torch.bfloat16
            elif 'cpu'  == device.type :
                bf16_supported = hasattr ( torch.cpu , 'is_bf16_supported' ) and torch.cpu.is_bf16_supported ()
                return bf16_supported , torch.bfloat16
            return False , torch.float32

        ## Autocast ? 
        use_autocast , autocast_dtype = get_autocast_config ( device )
        if not is_cuda :
            use_autocast, autocast_dtype  = False, torch.float32
            
        # =====================================================================
        # Main training loop
        # =====================================================================
        show_bar = self.__progress_epochs and not ( self.silent or self.progress )
        ## for epoch in progress_bar ( nepochs , silent = not show_bar , description = epoch_symbol ) :
        for _ , epoch , bar in text_bar ( nepochs + bce_epochs , silent = not show_bar , description = epoch_symbol ) :
            
            if bce_epochs == epoch : 
                val_loss_history.clear()
                best_loss        = float ( 'inf' )
                best_epoch       = 0 
                best_state       = None
                patience_counter = 0
                
                
            t_start = time.time() 

            model.train ()
            
            perm = torch.randperm ( n_samples , device=device)
            X_sh = X_tr_t [ perm ]
            y_sh = y_tr_t [ perm ]
            w_sh = w_tr_t [ perm ]
            
            t0 = time.time ()
            for i in range ( 0 , n_samples , batch_size ) :
                
                b_x = X_sh [ i : i + batch_size ]  # Zero-copy slicing!
                b_y = y_sh [ i : i + batch_size ]
                b_w = w_sh [ i : i + batch_size ]
            
                optimizer.zero_grad ( set_to_none = True )
                
                with torch.autocast ( device_type = device.type , dtype = autocast_dtype , enabled = use_autocast ) :
                    logits = model ( b_x )
                    
                    if bce_epochs <= epoch : loss = criterion_drw ( logits , b_y , b_w )
                    else : 
                        bce_vals = nn.functional.binary_cross_entropy_with_logits ( logits , b_y , reduction = 'none' )
                        loss     = torch.sum ( bce_vals * b_w ) / ( b_w.abs().sum() + 1e-8 )
                        if 0 < logit_reg : loss += logit_reg * torch.mean ( logits ** 2 ) 
                        
                if is_cuda :
                    scaler_amp.scale ( loss ).backward()
                    scaler_amp.unscale_( optimizer )
                    if 0 < grad_clip : torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=grad_clip)
                    scaler_amp.step   ( optimizer)
                    scaler_amp.update ( )
                else :
                    loss.backward()
                    if 0 < grad_clip : torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=grad_clip)
                    optimizer.step()
                        
                
            t_train = time.time() - t0

            ## EMA update
            ema_model.update_parameters ( model )
            
            scheduler.step ()

            # =================================================================
            # Validation using EMA model
            # =================================================================
            if  ( epoch      <= bce_epochs and epoch <= ma_window ) or \
                ( bce_epochs <= epoch   <= bce_epochs + ma_window ) or \
                ( epoch % eval_freq == 0 or nepochs <= epoch + 10 ) : 
                
                ema_model.eval ()
                val_loss_raw = 0.0
                
                with torch.inference_mode () :
                    
                    for i in range ( 0 , len ( y_val ) , batch_size ) :
                        v_x = X_va_t [ i : i + batch_size ]
                        v_y = y_va_t [ i : i + batch_size ]
                        v_w = w_va_t [ i : i + batch_size ]
                        
                        with torch.autocast ( device_type = 'cuda' if is_cuda else 'cpu' , enabled = is_cuda ) :
                            v_logits     = ema_model ( v_x )
                            
                            if epoch < bce_epochs :
                                v_bce = nn.functional.binary_cross_entropy_with_logits ( v_logits , v_y , reduction = 'none' )
                                v_loss_batch = torch.sum ( v_bce * v_w ) / ( v_w.abs().sum() + 1e-8 )
                                if 0 < logit_reg : v_loss_batch += logit_reg * torch.mean ( v_logits ** 2 )
                            else :
                                v_loss_batch = criterion_drw ( v_logits , v_y , v_w )
                            
                            val_loss_raw += v_loss_batch.item () * v_w.abs().sum().item()
                            
                    val_loss_raw /= w_va_sum.item ()

                val_loss_history.append ( val_loss_raw )
                while ma_window < len ( val_loss_history ) : val_loss_history.popleft() 
                
                val_loss_smoothed = sum ( val_loss_history ) / len ( val_loss_history )

                if show_bar :
                    rl = nice_print ( val_loss_raw               , precision = 3 )
                    sl = nice_print ( val_loss_smoothed          , precision = 3 )
                    lr = nice_print ( scheduler.get_last_lr()[0] , precision = 2 , with_sign = False ) 
                    bar.text = f"Epoch: {epoch:3d}/{nepochs} | Raw/Smoothed Losses:{rl}/{sl} | LR:{lr} | Best epoch: #{best_epoch} | Time:{time.time()-t_start:.2f}s/epoch" 

                if   epoch <= ma_window or  bce_epochs <= epoch <= bce_epochs + ma_window : continue
                elif val_loss_smoothed + min_delta < best_loss or best_state is None      :
                    best_loss        = val_loss_smoothed
                    unwrapped_ema    = ema_model.module if hasattr ( ema_model , 'module' ) else ema_model
                    best_state       = { k : v.clone().detach() for k , v in unwrapped_ema.state_dict().items() }
                    best_epoch       = epoch 
                    patience_counter = 0
                else :
                    patience_counter += eval_freq
                    if min_epochs <= epoch and patience <= patience_counter :
                        if show_bar :
                            bl       = nice_print ( best_loss , precision = 4 )
                            bar.text = f"Early stopping at epoch #{epoch}/{nepochs} | Best Smoothed Loss: {bl} | Best epoch: #{best_epoch}" 
                        break
                    
        # =====================================================================
        ## PHOENIX 
        # =====================================================================
        if True : # ===========================================================
            # =================================================================
            del model , ema_model , optimizer , scaler_amp , scheduler
            del X_tr_t , y_tr_t , w_tr_t , X_va_t , y_va_t , w_va_t
            
            gc.collect ()
            if is_cuda : torch.cuda.empty_cache ()
    
            final_model = TabularResNet ( **model_config ).to ( device )
            
            if best_state is not None :
                final_model.load_state_dict ( best_state )
                del best_state
                gc.collect ()
    
            final_model.eval ()
            final_model.scaler = scaler

        val_preds = self._predict_single_model ( final_model , X_val )
        return final_model , val_preds

    # =========================================================================
    ## Predicts probabilities for the given dataset.
    #
    #  @param model Trained PyTorch model.
    #  @param X     Features to predict on.
    #  @return Numpy array of predicted probabilities.
    # =========================================================================
    def _predict_single_model ( self , model , X ) :
        
        import torch
        
        model.eval ()
        device  = next ( model.parameters () ).device
        is_cuda = device.type == 'cuda'

        n_samples = num_samples ( X )
        
        batch_size = 2**17 if is_cuda else 2**16
        batch_size = max ( 1024 , 2 ** math.floor ( math.log2 ( max ( 1 , min ( batch_size , n_samples ) ) ) ) )

        X_scaled = model.scaler.transform ( numpy.ascontiguousarray ( X , dtype = numpy.float32 ) ) .astype ( numpy.float32 , copy = False)        
        preds    = numpy.empty ( n_samples , dtype = numpy.float32 )

        # =====================================================================
        # Batched inference loop
        # =====================================================================
        with torch.inference_mode () :
            for i in range ( 0 , n_samples , batch_size ) :
                chunk = X_scaled [ i : i + batch_size ]
                b_x   = torch.from_numpy ( chunk ).to ( device )
                
                with torch.autocast ( device_type = 'cuda' if is_cuda else 'cpu' , enabled = is_cuda ) :
                    logits  = model ( b_x )
                    b_preds = torch.sigmoid ( logits )
                
                preds [ i : i + batch_size ] = b_preds.cpu().numpy ()

        return preds

# =============================================================================
# Filter parameters to keep only those accepted by scikit-learn LogisticRegression
valid_LR_params = (
    'penalty'           ,
    'dual'              ,
    'tol'               ,
    'C'                 ,
    'fit_intercept'     ,
    'intercept_scaling' ,
    'class_weight'      ,
    'random_state'      ,
    'solver'            ,
    'max_iter'          ,
    'multi_class'       ,
    'verbose'           ,
    'warm_start'        ,
    'l1_ratio'          ,
    'n_jobs' 
)
# =============================================================================
## @class LogRegressionDensityReweighter
#  Density ratio reweighter using Logistic Regression as the underlying classifier.
class LogRegressionDensityReweighter ( DensityReweighter ) :
    """ Density ratio reweighter using Logistic Regression as the underlying classifier.
    """
    # =========================================================================
    ## Initialize Logistic Regression density reweighter.
    #  @param original Features array for original sample.
    #  @param target Features array for target sample.
    #  @param original_weight Initial weights for original sample (optional).
    #  @param target_weight Initial weights for target sample (optional).
    #  @param store_original_weights If True, store computed results for original sample.
    #  @param params Additional Logistic Regression parameters.
    def __init__( self , * , 
                  original               ,
                  target                 ,
                  original_weight        = None ,
                  target_weight          = None ,
                  store_original_weights = True , **params ) :        
        """ Initialize Logistic Regression density reweighter.
        """
        config = { 'C'             : 1.0     ,
                   'solver'        : 'lbfgs' ,
                   'max_iter'      : 2000    ,
                   'polynomials'   :    2    }
        
        config.update ( params )
        
        poly       = config.get   ( 'polynomials', 2 )
        n_features = num_features ( original         )
        
        if 15 <= n_features and 1 < poly :
            config [ 'polynomials' ] = 1
        
        super().__init__ ( original               = original               ,
                           target                 = target                 ,
                           original_weight        = original_weight        ,
                           target_weight          = target_weight          ,
                           store_original_weights = store_original_weights , **config )
        
    # =========================================================================
    @property
    def method ( self ) :
        return method_LR 

    # =========================================================================
    ## Train single Logistic Regression model on fold data.
    def _train_single_model ( self    ,
                              X_train , y_train , w_train ,
                              X_val   , y_val   , w_val   ) :
        """ Train single Logistic Regression model on fold data.
        """
        from sklearn.linear_model  import LogisticRegression
        from sklearn.pipeline      import Pipeline        
        from sklearn.preprocessing import StandardScaler, PolynomialFeatures
        
        params = {}
        params.update ( self.params    )
        params.pop    ( 'n_jobs', None ) 

        poly = params.pop ( 'polynomials', 2 )
        n_features = num_features ( X_train )
        if 20 <= n_features and 1 < poly : poly = 1
        
        lr_kwargs = { k: v for k, v in params.items() if k in valid_LR_params }

        # Build the pipeline 
        steps = [ ( 'scaler1' , StandardScaler () ) ] 

        if poly and isinstance ( poly , int ) and 1 < poly <= 5 :
            steps += [ ( 'poly'    , PolynomialFeatures ( degree = poly  , include_bias = False ) ) ]
            steps += [ ( 'scaler2' , StandardScaler () ) ]
            
        lr_model = LogisticRegression ( **lr_kwargs )
        steps += [ ( 'logistic', lr_model ) ]
        
        model = Pipeline ( steps  )

        w_tr = numpy.ascontiguousarray ( w_train , dtype = numpy.float32 ) if w_train is not None else None

        if w_tr is not None : model.fit ( X_train , y_train , logistic__sample_weight = w_tr )
        else                : model.fit ( X_train , y_train )

        # =====================================================================
        ## PHOENIX
        # =====================================================================
        if True : # ===========================================================
            # =================================================================
            import pickle 
        
            model_bytes = pickle.dumps ( model ) 
            del model
            del w_tr
            
            model = pickle.loads ( model_bytes )
            del model_bytes 

            gc.collect()
        
        val_preds = self._predict_single_model ( model , X_val )
        return model , val_preds.astype ( numpy.float32 , copy = False )

    # =========================================================================
    ## Predict probabilities using a Logistic Regression model.
    def _predict_single_model ( self , model , X ) :
        """ Predict probabilities using a Logistic Regression model.
        """
        X_clean = numpy.ascontiguousarray ( X , dtype = numpy.float32 )
        p       = model.predict_proba ( X_clean ) [ : , 1 ]
        return p.astype ( numpy.float32 , copy = False )

# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    from ostap.stats.tools import ( hasLightGBM ,
                                    hasXGBoost  ,
                                    hasCatBoost ,
                                    hasPyTorch  ,
                                    hasSkLearn  )

    if not hasLightGBM ( False ) : logger.warning  ( "No LightGBM available!" ) 
    if not hasXGBoost  ( False ) : logger.warning  ( "No XGBoost  available!" ) 
    if not hasCatBoost ( False ) : logger.warning  ( "No CatBoost available!" ) 
    if not hasPyTorch  ( False ) : logger.warning  ( "No PyTorch  available!" ) 
    if not hasSkLearn  ( False ) : logger.warning  ( "No SkLearn  available!" ) 

# =============================================================================
#                                                                       The END 
# =============================================================================
