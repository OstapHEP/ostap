#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/dreweighters.py
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
    'PyTorchAdamDensityReweighter'   , # PyTorch-based density reweighter
    'PyTorchLBFGSDensityReweighter'  , # PyTorch-based density reweighter
    'LogRegressionDensityReweighter' , # Density reweighter based on Logistic Regression
) 
# =============================================================================
from   ostap.core.ostap_types       import string_types, integer_types  
from   ostap.utils.basic            import numcpu, classprop, safe_loop
from   ostap.stats.utils            import num_features, num_samples        
from   ostap.stats.tools            import hasSkLearn 
from   ostap.tools.base_reweighters import DensityReweighter  
from   ostap.stats.utils            import optimal_batch_size
from   ostap.utils.progress_bar     import progress_bar 
import ostap.stats.utils            as     U
import ostap.logger.symbols         as     S
import numpy, math, abc 
# =============================================================================
# Logging setup
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
logger = getLogger ( 'ostap.tools.dreweighters' )
# =============================================================================
has_sklearn = hasSkLearn () 
# =============================================================================
## Decorations 
# =============================================================================
method_LGBM        =  U.method_name ( 'LightGBM/Density' )
method_XGB         =  U.method_name ( 'XGBoost/Density'  )
method_CATB        =  U.method_name ( 'CatBoost/Density' )
method_LOGREG      =  U.method_name ( 'Logistic/Density' )
method_TORCH       =  U.method_name ( 'PyTorch/Density'  )
method_TORCH_ADAM  =  U.method_name ( 'PyTorch/Adam/Density'   )
method_TORCH_LBFGS =  U.method_name ( 'PyTorch/L-BFGS/Density' )
train_symbol       =  S.hammer_and_wrench if S.show and S.hammer_and_wrench else U.symbol_name ( 'Train' ) 
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

CATB_DEFAULT_ESTIMATORS        = 1600 # High capacity: deep ensemble with low learning rate
CATB_MAX_DEPTH                 =   10 # Deep tree capacity for high dimensions / large statistics
CATB_LEARNING_RATE             = 0.005 # Small step for stable density ratio convergence

# =============================================================================
THRESHOLD_NSAMPLES    = 50000
THRESHOLD_NFEATURES   = 10

# =============================================================================
## @class LightGBMDensityReweighter
#  Density ratio reweighter using LightGBM as the underlying classifier.
class LightGBMDensityReweighter ( DensityReweighter ) :
    """ Density ratio reweighter using LightGBM as the underlying classifier
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
        """ Initialize LightGBM density reweighter (High Capacity Mode Defaults)
        """
        config = {
            'objective'          : 'binary'                ,
            'metric'             : 'binary_logloss'        ,
            'n_estimators'       : LGBM_DEFAULT_ESTIMATORS ,
            'learning_rate'      : LGBM_LEARNING_RATE      ,
            'max_depth'          : LGBM_MAX_DEPTH          ,
            'num_leaves'         : min ( 1023 , 2 ** LGBM_MAX_DEPTH - 1 ) , 
            'max_bin'            : 2047                    ,
            'min_child_samples'  : 10                      , 
            'boost_from_average' : True                    ,
            'verbosity'          : -1                      ,
            'n_jobs'             : -1                      ,
        }
        
        config.update ( kwargs )

        if 'silent' in config and config.get ( 'silent' ) : config [ 'verbose' ] = -1 
        
        if 'num_boost_round' in config  : config [ 'n_estimators'     ] = config.pop ( 'num_boost_round' , LGBM_DEFAULT_ESTIMATORS ) 
        if 'bagging_fraction' in config : config [ 'subsample'        ] = config.pop ( 'bagging_fraction'  , 1.0 )         
        if 'colsample_bytree' in config : config [ 'feature_fraction' ] = config.pop ( 'colsample_by_tree' , 1.0 )         
        
        n_rows = num_samples ( original ) + num_samples ( target )
        if n_rows < THRESHOLD_NSAMPLES or 1 < config.get ( 'n_splits' , 5 )     : config [ 'subsample'        ] = 1.0    
        if 1 != config.get ( 'subsample' , 1 ) and not 'bagging_freq' in config : config [ 'bagging_freq'     ] = 1            
        if num_features    ( original ) < THRESHOLD_NFEATURES                   : config [ 'feature_fraction' ] = 1.0
            
        super ().__init__ (
            original        = original        ,
            target          = target          ,
            original_weight = original_weight ,
            target_weight   = target_weight   ,
            **config
        )

    # =========================================================================
    ## Short symbol or tag representing the LightGBM reweighter method.
    #  @return String identifier/symbol used for logging and visual presentation.
    @classprop
    def symbol ( self ) :
        """ Short symbol/tag for reweighter. """
        return U.symbol_LGBM
        
    # =========================================================================
    ## Return the method identifier name. 
    @property
    def method ( self ) :
        """ Return the method identifier name. """
        return method_LGBM 

    # =========================================================================
    ## Train single LightGBM model on fold data.
    #  @param X_train Training features array (N-1)/N.
    #  @param y_train Training binary labels.
    #  @param w_train Training event weights or None.
    #  @return Trained LightGBM booster instance.
    def _train_single_model ( self, X_train, y_train, w_train ) : 
        """ Train single LightGBM model on fold data. """
        import lightgbm as LightGBM
    
        X_tr = numpy.ascontiguousarray ( X_train , dtype = numpy.float32 )
        w_tr = numpy.ascontiguousarray ( w_train , dtype = numpy.float32 ) if w_train is not None else None

        trn_data = LightGBM.Dataset ( X_tr, label = y_train, weight = w_tr, free_raw_data = True )

        params = self.params.copy ()
        num_boost_round          = params.pop ( 'num_boost_round', None ) or params.pop ( 'n_estimators', LGBM_DEFAULT_ESTIMATORS )
        params [ 'num_threads' ] = params.get ( 'num_threads', params.pop ( 'n_jobs', max ( 1, numcpu () // 2 ) ) ) 

        model = LightGBM.train ( params          = params          ,
                                 train_set       = trn_data        ,
                                 num_boost_round = num_boost_round )

        # PHOENIX memory optimization: export raw C++ model string and reload clean booster
        if self.phoenix :
            model_bytes = model.model_to_string ()
            best_iter   = getattr ( model, 'best_iteration_', None ) or getattr ( model, 'best_iteration', None )                    
            model.free_dataset ()
            del model 
            model       = LightGBM.Booster ( model_str = model_bytes )
            if best_iter is not None : model.best_iteration = best_iter 
            del model_bytes
            
        del trn_data, X_tr, w_tr

        return model

    # =========================================================================
    ## Predict probabilities using a LightGBM model.
    #  @param model Trained LightGBM booster instance.
    #  @param X Features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model ( self, model, X ) :
        """ Predict probabilities using a LightGBM model.
        """
        best_iter   = getattr ( model, 'best_iteration_', None ) or getattr ( model, 'best_iteration', None )        
        kwargs      = {}
        if best_iter is not None and 0 < best_iter : kwargs [ 'num_iteration' ] = best_iter            
        num_threads = self.params.get ( 'num_threads', self.params.get ( 'n_jobs', max ( 1, numcpu () // 2 ) ) )                                    
        p           = model.predict ( X, num_threads = num_threads, **kwargs )
        return p.astype ( numpy.float32, copy = False )

# =============================================================================
## @class XGBoostDensityReweighter
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
            'objective'        : 'binary:logistic'     ,
            'eval_metric'      : 'logloss'             ,
            'n_estimators'     : XGB_DEFAULT_ESTIMATORS,
            'learning_rate'    : XGB_LEARNING_RATE     ,
            'max_depth'        : XGB_MAX_DEPTH         ,
            'tree_method'      : 'hist'                ,
            'max_bin'          : 1023                  ,
            'min_child_weight' : 1.e-3                 ,
            'verbosity'        : 0                     ,      
            'n_jobs'           : -1                    ,
        }
        config.update ( kwargs )
        
        if 'silent' in config and config.get ( 'silent' ) : config [ 'verbosity' ] = 0 

        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' )         
        n_rows = num_samples ( original ) + num_samples ( target )
        if n_rows < THRESHOLD_NSAMPLES or 1 < config.get ( 'n_splits' , 5 ) : config [ 'subsample'         ] = 1.0
        if num_features      ( original ) < THRESHOLD_NFEATURES             : config [ 'colsamples_bytree' ] = 1.0
        
        super ().__init__ (
            original        = original        ,
            target          = target          ,
            original_weight = original_weight ,
            target_weight   = target_weight   ,
            **config
        )

    # =========================================================================
    ## Short symbol or tag representing the LightGBM reweighter method.
    #  @return String identifier/symbol used for logging and visual presentation.
    @classprop
    def symbol ( self ) :
        """ Short symbol/tag for reweighter. """
        return U.symbol_XGB
        
    # =========================================================================
    ## Return the method identifier name.
    @property
    def method ( self ) :
        """ Return the method identifier name. """
        return method_XGB 

    # =========================================================================
    ## Train single XGBoost booster on fold data.
    #  @param X_train Training features array (N-1)/N.
    #  @param y_train Training binary labels.
    #  @param w_train Training event weights or None.
    #  @return Trained XGBoost Booster instance.
    def _train_single_model ( self, X_train, y_train, w_train ) :
        """ Train single XGBoost booster on fold data.
        """
        import xgboost as XGBoost
        
        X_tr   = numpy.ascontiguousarray ( X_train, dtype = numpy.float32 )
        y_tr   = numpy.ascontiguousarray ( y_train, dtype = numpy.float32 )
        w_tr   = numpy.ascontiguousarray ( w_train, dtype = numpy.float32 ) if w_train is not None else None

        dtrain = XGBoost.DMatrix ( X_tr, label = y_tr, weight = w_tr )
        
        params = self.params.copy ()
        params [ 'nthread' ] = params.get ( 'nthread', params.pop ( 'n_jobs', max ( 1, numcpu () // 2 ) ) ) 
        n_estimators         = params.pop ( 'n_estimators', 1000 )

        model = XGBoost.train (
            params,
            dtrain,
            num_boost_round = n_estimators,
            verbose_eval    = False
        )

        # PHOENIX memory optimization: export raw JSON model bytes and reload clean inference engine
        if self.phoenix : 
            model_bytes = model.save_raw ( raw_format = "json" )
            best_iter   = getattr ( model, 'best_iteration_', None ) or getattr ( model, 'best_iteration', None )                    
            del model
            model = XGBoost.Booster ()
            model.load_model ( model_bytes )
            del model_bytes            
            if best_iter is not None : model.best_iteration = best_iter
            
            
        del dtrain, X_tr, y_tr, w_tr

        return model
    
    # =========================================================================
    ## Predict probabilities using an XGBoost model.
    #  @param model Trained XGBoost Booster.
    #  @param X Input features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model ( self, model, X ) :
        """ Predict probabilities using an XGBoost model.
        """
        import xgboost as XGBoost
        dmat      = XGBoost.DMatrix ( X )
        best_iter = getattr ( model, 'best_iteration_', None ) or getattr ( model, 'best_iteration', None )
        kwargs    = {}
        if best_iter is not None and 0 < best_iter : kwargs [ 'iteration_range' ] = ( 0, best_iter + 1 )
        p         = model.predict ( dmat, **kwargs )
        return p.astype ( numpy.float32, copy = False )

# =============================================================================
## @class CatBoostDensityReweighter
#  Density ratio reweighter using CatBoost as the underlying classifier.
class CatBoostDensityReweighter ( DensityReweighter ) :
    """ Density ratio reweighter using CatBoost as the underlying classifier.
    """
    
    def __init__ ( self                   , * , 
                   original               ,
                   target                 ,
                   original_weight        = None ,
                   target_weight          = None ,
                   store_original_weights = True , **params ) :
        """ Initialize CatBoost density reweighter.
        """
        
        from catboost.utils import get_gpu_device_count
        gpu        = get_gpu_device_count () 
        self.__gpu = gpu 
        
        config = { 'loss_function'       : 'Logloss'                        ,
                   'eval_metric'         : 'Logloss'                        ,
                   'iterations'          : CATB_DEFAULT_ESTIMATORS          ,
                   'learning_rate'       : CATB_LEARNING_RATE               ,
                   'max_depth'           : CATB_MAX_DEPTH                   ,
                   'l2_leaf_reg'         : 0.1                              ,
                   'min_data_in_leaf'    : 5                                ,
                   'mvs_reg'             : 2.0                              , 
                   'random_strength'     : 0.01                             ,
                   'boosting_type'       : 'Plain'                          ,
                   ## 'bootstrap_type'      : 'MVS'                         ,
                   'bootstrap_type'      : 'Bayesian'                       , ## MVS ? 
                   'bagging_temperature' : 1.0                              , 
                   'device'              : 'GPU' if 0 < self.gpu else 'CPU' , 
                   'verbose'             : False
                  }
        
        config.update ( params )
        
        if 'min_child_samples' in config : config [ 'min_data_in_leaf' ] = config.pop ( 'min_child_samples' )
        if 'n_jobs'            in config : config [ 'thread_count'     ] = config.pop ( 'n_jobs'            )
        if 'n_estimators'      in config : config [ 'iterations'       ] = config.pop ( 'n_estimators'      )
        if 'depth'             in config : config [ 'max_depth'        ] = config.pop ( 'depth'             ) 
        
        if 'num_boost_round'  in config : config [ 'n_estimators'      ] = config.pop ( 'num_boost_round'      ) 
        if 'colsample_bytree' in config : config [ 'rsm'               ] = config.pop ( 'colsample_by_tree' , 1.0 )         

        n_rows = num_samples ( original ) + num_samples ( target )
        if n_rows < THRESHOLD_NSAMPLES or 1 < config.get ( 'n_splits' , 5 ) : config.pop ( 'subsample' , None )
        if num_features    ( original ) < THRESHOLD_NFEATURES               : config.pop ( 'rsm'       , None ) 

        if self.gpu <= 0 : config [ 'device' ] = 'CPU'

        if 'bayesian' == config.get( 'bootstrap_type' , '' ).strip().lower()  :
            config.pop ( 'subsample' , None )
            config.pop ( 'rsm'       , None )
        
            
        super ().__init__ ( original               = original               ,
                            target                 = target                 ,
                            original_weight        = original_weight        ,
                            target_weight          = target_weight          ,
                            store_original_weights = store_original_weights , **config ) 
        
    # =========================================================================
    @property
    def gpu ( self ) :
        """ `gpu`: result of catboost.utils.get_gpu_device_count """
        return self.__gpu
    
    # =========================================================================
    ## Short symbol or tag representing the LightGBM reweighter method.
    #  @return String identifier/symbol used for logging and visual presentation.
    @classprop
    def symbol ( self ) :
        """ Short symbol/tag for reweighter
        """
        return U.symbol_CATB

    # =========================================================================
    ## Return method name identifier
    @property
    def method ( self ) :
        """ Return method name identifier
        """
        return method_CATB
    
    # =========================================================================
    ## Train single CatBoost model on fold data.
    #  @param X_train Training features array (N-1)/N.
    #  @param y_train Training binary labels.
    #  @param w_train Training event weights or None.
    #  @return Trained CatBoostClassifier instance.
    def _train_single_model ( self, X_train, y_train, w_train ) :
        """ Train single CatBoost model on fold data. """
        import gc
        import pickle
        import catboost as CatBoost
        
        X_tr = numpy.ascontiguousarray ( X_train, dtype = numpy.float32 )
        y_tr = numpy.ascontiguousarray ( y_train, dtype = numpy.int32   )
        w_tr = numpy.ascontiguousarray ( w_train, dtype = numpy.float32 ) if w_train is not None else None
    
        trn_pool = CatBoost.Pool ( X_tr, label = y_tr, weight = w_tr )
        params   = self.params.copy ()
        
        params [ 'depth'         ] = params.pop ( 'max_depth'    , CATB_MAX_DEPTH ) 
        params [ 'thread_count'  ] = params.get ( 'thread_count' , params.pop ( 'n_jobs', max ( 1, numcpu () // 2 ) ) ) 
        params [ 'boosting_type' ] = params.get ( 'boosting_type', 'Plain' )
        
        iterations = params.pop ( 'n_estimators', params.pop ( 'iterations', CATB_DEFAULT_ESTIMATORS ) ) 
        if not isinstance ( iterations, int ) or not 10 <= iterations <= 10000 :
            iterations = CATB_DEFAULT_ESTIMATORS
        
        params [ 'iterations' ] = iterations
        params [ 'task_type'  ] = params.pop ( 'device', 'GPU' if 0 < self.gpu else 'CPU' )
            
        model         = CatBoost.CatBoostClassifier ( **params )
        verbose_level = params.get ( 'verbose', False )
        model.fit ( trn_pool, verbose = verbose_level )

        # PHOENIX memory optimization: serialize model blob and purge heavy training pools
        if self.phoenix : 
            model_bytes = pickle.dumps ( model ) 
            best_iter   = getattr ( model, 'best_iteration_', None ) or getattr ( model, 'best_iteration', None )                    
            del model
            model = pickle.loads ( model_bytes )
            if best_iter is not None : model.best_iteration = best_iter
            del model_bytes  
            
        del trn_pool, X_tr, y_tr, w_tr
        
        return model

    # =========================================================================
    ## Predict probabilities using a CatBoost model.
    #  @param model Trained CatBoost model.
    #  @param X Features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model ( self, model, X ) :
        """ Predict probabilities using a CatBoost model.
        """
        X_clean   = numpy.ascontiguousarray ( X, dtype = numpy.float32 )
        best_iter = getattr ( model, 'best_iteration_', None ) or getattr ( model, 'best_iteration', None )        
        kwargs    = {}
        if best_iter is not None and 0 <= best_iter : kwargs [ 'ntree_end' ] = best_iter + 1
        thread_count = self.params.get ( 'thread_count', self.params.get ( 'n_jobs', max ( 1, numcpu () // 2 ) ) ) 
        p            = model.predict_proba ( X_clean, thread_count = thread_count, **kwargs ) [ :, 1 ]
        return p.astype ( numpy.float32, copy = False )


    
    
# ==========================================================================
TORCH_OPTIMIZER_ADAM  = 'adam'
TORCH_OPTIMIZER_LBFGS = 'lbfgs'
TORCH_OPTIMIZERS      = TORCH_OPTIMIZER_ADAM , TORCH_OPTIMIZER_LBFGS
TORCH_OPTIMIZER       = TORCH_OPTIMIZER_LBFGS 
TORCH_HIDDEN_LAYERS   = ( 256 , 256 ) ## single  hiddle layer with 32 nodes 
TORCH_MAX_BATCH_SIZE  = 2 ** 16       ## maximal batch size
TORCH_LEARNING_RATE   = 1e-4     
TORCH_WEIGHT_DECAY    = 1e-5     
TORCH_LBFGS_CHUNK     = 200_000
TORCH_EPOCHS          = 1000    

# =============================================================================
## @class PyTorchBaseDensityReweighter
#  (Abstract) base class for PyTorch MLP-based density ratio reweighters.
#
#           Implements a single-stage neural density-ratio estimator on top
#           of @c DensityReweighter by training a Multi-Layer Perceptron
#           (@c WeakMLP) per cross-validation fold against a direct density-ratio
#           objective (@c DensityRatioLoss) rather than standard binary cross-
#           entropy
#
#           Shared responsibilities—device resolution, batch-size calculation,
#           network and loss construction, fold model orchestration, and
#           batched probability inference with @c torch.inference_mode and
#           CUDA @c torch.autocast—live in this base class, while optimizer-
#           specific training loops are delegated to subclasses via @c _fit_model.
# =============================================================================
class PyTorchDensityReweighter (DensityReweighter ) :
    """ Abstract base class for PyTorch MLP-based density ratio reweighters

    Coordinates device normalization, batch size selection, ``WeakMLP`` and
    ``DensityRatioLoss`` instantiation, fold-level training via ``_fit_model``,
    and batched probability prediction via ``_predict_single_model``
    """
    # =========================================================================
    ## Initialize the PyTorch density ratio reweighter base class.
    #
    #           Merges caller keyword arguments with default PyTorch density
    #           reweighting parameters, normalizes the target device to a
    #           @c torch.device instance (CUDA GPU or CPU), computes the
    #           effective training batch size, and invokes @c DensityReweighter.__init__
    #           with @c progress=False so fold-level progress bars do not collide
    #           with per-epoch bars when present
    #
    #  @param original Source feature matrix (N-dimensional array or dataframe)
    #  @param target Target feature matrix (N-dimensional array or dataframe)
    #  @param original_weight Optional source observation weights (supports signed weights)
    #  @param target_weight Optional target observation weights (supports signed weights)
    #  @param progress Optional progress flag retained for subclass/epoch progress control.
    #  @param kwargs Additional keyword arguments controlling network architecture,
    #         learning rate, regularization, logit bounds, and execution device.
    def __init__ ( self                      , * ,
                   original                  ,
                   target                    ,
                   original_weight = None    ,
                   target_weight   = None    ,
                   progress        = True    , **kwargs ) :
        """ Initialize ``PyTorchDensityReweighter`` base configuration

        Resolves target ``torch.device``, determines batch size, and delegates
        cross-validation setup and fitting to ``DensityReweighter``
        """
        import torch

        config = { 'hidden_layers' : TORCH_HIDDEN_LAYERS ,
                   'learning_rate' : TORCH_LEARNING_RATE ,
                   'weight_decay'  : TORCH_WEIGHT_DECAY  ,
                   'epochs'        : TORCH_EPOCHS        ,
                   'max_logit'     : 10                  ,
                   'logit_reg'     : 0                   ,
                   'device'        : None                ,
                   'n_jobs'        : -1                  ,
                  }
        
        config.update ( kwargs )
        self.__progress_epochs = bool ( config.get ( 'progress', progress ) )

        # Normalize device specification into an explicit torch.device object
        device = config.get ( 'device', None )
        if device is None : config [ 'device' ] = torch.device ( 'cuda' if torch.cuda.is_available () else 'cpu' )
        elif isinstance ( device, string_types ) and device.strip ().lower () in ( 'gpu', 'cuda' ) :
            config [ 'device' ] = torch.device ( 'cuda' if torch.cuda.is_available () else 'cpu' )
        else :
            config [ 'device' ] = torch.device ( device )

        n_rows     = num_samples  ( target ) + num_samples ( original )
        n_features = num_features ( target )

        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   ,
                            progress        = False           , **config )

    # =========================================================================
    ## Return whether epoch-level progress reporting is enabled.
    #  @return Boolean flag controlling epoch progress display in subclasses
    @property
    def progress_epochs ( self ) :
        """ Return epoch-level progress bar flag
        """
        return self.__progress_epochs

    # =========================================================================
    ## Get the short symbol/tag representing the PyTorch reweighter method.
    #  @return String identifier symbol used for logging and visual presentation
    @classprop
    def symbol ( self ) :
        """ Return short symbol/tag for PyTorch density reweighter
        """
        return U.symbol_TORCH

    # =========================================================================
    ## @brief Get the full method identifier name.
    #  @return Method string tag used in summaries and diagnostic logs
    @property
    def method ( self ) :
        """ Return method identifier name for PyTorch density reweighter
        """
        return method_TORCH

    # =========================================================================
    ## Build the direct density-ratio loss module (@c DensityRatioLoss).
    #
    #           Estimates log density ratios directly by minimizing
    #           \f$\frac{\sum_{i \in \text{orig}} w_i \exp(\text{clip}(z_i, z_{\max}))}
    #           {\sum_{i \in \text{orig}} |w_i|} - \frac{\sum_{j \in \text{targ}} w_j z_j}{\sum_{j \in \text{targ}} |w_j|}\f$
    #           plus an optional quadratic logit penalty @c logit_reg * mean(z^2)
    #           to prevent logit blow-up in low-density regions
    #
    #  @param max_logit Upper clamp bound applied to source logits before @c exp
    #  @param logit_reg Weight of the quadratic logit regularization term
    #  @return Instantiated @c torch.nn.Module computing the density-ratio loss
    def _build_loss_module ( self             ,
                             max_logit = 8.0  ,
                             logit_reg = 1e-4 ) :
        """ Construct and return the ``DensityRatioLoss`` ``torch.nn.Module``
        
        Minimizes expected exponential source logits minus expected target logits
        with optional quadratic regularization on raw logits
        """
        import torch
        import torch.nn as nn

        ## Direct density-ratio estimation loss with logit clamping and L2 penalty
        class DensityRatioLoss ( nn.Module ) :
            """ Direct density-ratio estimation loss with logit clamping and L2 penalty
            """

            ## Initialize clamping bound and quadratic logit regularization strength
            def __init__ ( self, max_logit = 8.0, logit_reg = 1e-4 ) :
                """ Initialize clamping bound and quadratic logit regularization strength
                """
                super ().__init__ ()
                self.max_logit = max_logit
                self.logit_reg = logit_reg

            ## Compute weighted density-ratio loss over combined batch samples
            def forward ( self, logits, targets, weights ) :
                """ Compute weighted density-ratio loss over combined batch samples
                """
                logits  = logits.float ()
                targets = targets.float ()
                weights = weights.float ()

                mask_orig = ( targets == 0 )
                mask_targ = ( targets == 1 )

                z_orig = logits  [ mask_orig ]
                w_orig = weights [ mask_orig ]
                z_targ = logits  [ mask_targ ]
                w_targ = weights [ mask_targ ]

                # Guard denominators against zero or degenerate weight subsets
                w_orig_sum = w_orig.abs ().sum () + 1e-8
                w_targ_sum = w_targ.abs ().sum () + 1e-8

                exp_z     = torch.exp ( torch.clamp ( z_orig, max = self.max_logit ) )
                loss_orig = torch.sum ( w_orig * exp_z ) / w_orig_sum
                loss_targ = torch.sum ( w_targ * z_targ ) / w_targ_sum

                loss = loss_orig - loss_targ
                if 0 < self.logit_reg :
                    loss = loss + self.logit_reg * torch.mean ( logits ** 2 )
                return loss

        return DensityRatioLoss ( max_logit = max_logit, logit_reg = logit_reg )

    # =========================================================================
    ## Build the @c WeakMLP neural network module class.
    #
    #           Constructs a @c torch.nn.Module comprising a @c BatchNorm1d
    #           input normalizer, hidden @c Linear + @c LeakyReLU(0.1) blocks
    #           initialized via Kaiming uniform initialization, and a final
    #           @c Linear output head with @c bias=False so outputs remain centered
    #           around zero without an unconstrained global bias term
    #
    #  @return Subclass of @c torch.nn.Module implementing @c WeakMLP
    def _build_mlp_class ( self ) :
        """ Construct and return the ``WeakMLP`` ``torch.nn.Module``

        Builds an MLP with ``BatchNorm1d`` input normalization, ``LeakyReLU``
        hidden layers, and a bias-free linear output layer
        """
        import torch
        import torch.nn as nn

        ## Multi-Layer Perceptron module predicting raw log density ratios
        class WeakMLP ( nn.Module ) :
            """ Multi-Layer Perceptron module predicting raw log density ratios
            """

            ## Initialize BatchNorm input stage, hidden layers, and output head
            def __init__ ( self, input_dim, hidden_layers = TORCH_HIDDEN_LAYERS ) :
                """ Initialize BatchNorm input stage, hidden layers, and output head
                """
                super ().__init__ ()

                # Allow passing either a single layer width or an iterable/tuple of widths
                if isinstance ( hidden_layers, integer_types ) : hidden_layers = ( hidden_layers, )

                self.bn = nn.BatchNorm1d ( input_dim )

                layers   = []
                curr_dim = input_dim
                for h_dim in hidden_layers :
                    linear = nn.Linear ( curr_dim, h_dim )
                    nn.init.kaiming_uniform_ ( linear.weight, a = 0.1 )
                    nn.init.zeros_ ( linear.bias )

                    layers.append ( linear )
                    layers.append ( nn.LeakyReLU ( 0.1 ) )
                    curr_dim = h_dim

                self.hidden_net = nn.Sequential ( *layers )

                # bias=False keeps the linear readout anchored without a free offset
                self.fc_out = nn.Linear ( curr_dim, 1, bias = False )
                nn.init.xavier_uniform_ ( self.fc_out.weight )

            ## Execute forward pass and return 1D raw logit predictions
            def forward ( self, x ) :
                """ Execute forward pass and return 1D raw logit predictions
                """
                h = self.bn ( x )
                h = self.hidden_net ( h )
                return self.fc_out ( h ).squeeze ( -1 )

        return WeakMLP

    # =========================================================================
    ## Train ``model`` on fold data using the subclass's optimizer strategy.
    #
    #           Hook method invoked by @c _train_single_model after contiguous
    #           NumPy arrays and device tensors have been prepared.
    #           Subclasses implement optimizer-specific parameter updates here
    #
    #  @param model Instantiated @c WeakMLP network placed on @c device
    #  @param criterion Loss module returned by @c _build_loss_module
    #  @param X_tr Contiguous float32 NumPy feature array of shape @c (n_samples, n_features)
    #  @param y_tr Float32 NumPy binary label array of shape @c (n_samples,)
    #  @param w_tr Contiguous float32 NumPy weight array of shape @c (n_samples,)
    #  @param X_tr_t Feature tensor on @c device
    #  @param y_tr_t Label tensor on @c device
    #  @param w_tr_t Sample weight tensor on @c device
    #  @param device Target @c torch.device for training
    @abc.abstractmethod
    def _fit_model ( self      ,
                     model     ,
                     criterion ,
                     X_tr      ,
                     y_tr      ,
                     w_tr      ,
                     X_tr_t    ,
                     y_tr_t    ,
                     w_tr_t    ,
                     device    ) :
        """ Optimize ``model`` on the training fold using a concrete optimizer
        """
        raise NotImplementedError

    # =========================================================================
    ## Train a single PyTorch density-ratio model on one CV fold.
    #
    #           Converts fold arrays to contiguous @c float32 tensors on the
    #           configured device, instantiates @c WeakMLP and @c DensityRatioLoss,
    #           delegates parameter optimization to @c _fit_model, and switches
    #           the trained model into @c eval() mode before returning it
    #
    #  @param X_train Training features matrix for the current fold @c (N-1)/N
    #  @param y_train Training binary labels (0 for original, 1 for target)
    #  @param w_train Training sample weights or @c None
    #  @return Trained @c WeakMLP model in evaluation mode
    def _train_single_model ( self, X_train, y_train, w_train ) :
        """ Train a single fold ``WeakMLP`` model using ``DensityRatioLoss``
        """
        import torch

        device = torch.device ( self.params.get ( 'device', 'cuda' if torch.cuda.is_available () else 'cpu' ) )

        n_samples  = num_samples  ( X_train )
        n_features = num_features ( X_train )

        max_logit     = self.params.get ( 'max_logit'    , 8.0                 )
        logit_reg     = self.params.get ( 'logit_reg'    , 1e-4                )
        hidden_layers = self.params.get ( 'hidden_layers', TORCH_HIDDEN_LAYERS )

        X_tr = numpy.ascontiguousarray ( X_train, dtype = numpy.float32 )
        w_tr = (
            numpy.ascontiguousarray ( w_train, dtype = numpy.float32 )
            if w_train is not None
            else numpy.ones ( n_samples, dtype = numpy.float32 )
        )
        y_tr = y_train.astype ( numpy.float32 )

        X_tr_t = torch.from_numpy ( X_tr ).to ( device )
        y_tr_t = torch.from_numpy ( y_tr ).to ( device )
        w_tr_t = torch.from_numpy ( w_tr ).to ( device )

        WeakMLPClass = self._build_mlp_class ()
        model        = WeakMLPClass ( n_features, hidden_layers ).to ( device )
        criterion    = self._build_loss_module ( max_logit = max_logit, logit_reg = logit_reg )

        self._fit_model ( model     = model     ,
                          criterion = criterion ,
                          X_tr      = X_tr      ,
                          y_tr      = y_tr      ,
                          w_tr      = w_tr      ,
                          X_tr_t    = X_tr_t    ,
                          y_tr_t    = y_tr_t    ,
                          w_tr_t    = w_tr_t    ,
                          device    = device    )

        model.eval ()
        return model

    # =========================================================================
    ## Predict probabilities @c p(x) using a trained fold @c WeakMLP.
    #
    #           Evaluates @c model in chunks under @c torch.inference_mode and
    #           @c torch.autocast (when running on CUDA), converting raw logits
    #           via @c torch.sigmoid into a @c float32 NumPy probability array
    #
    #  @param model Trained @c WeakMLP model instance
    #  @param X Feature matrix of samples to score
    #  @return @c float32 NumPy array of predicted probabilities @c p(x)
    def _predict_single_model ( self, model, X ) :
        """ Predict probabilities ``p(x)`` in batches using a trained model
        """
        import torch

        model.eval ()
        device     = next ( model.parameters () ).device
        is_cuda    = device.type == 'cuda'
        n_samples  = num_samples ( X )

        max_bs     = 2 ** 17 if is_cuda else 2 ** 16
        batch_size = max ( 1024, 2 ** math.floor ( math.log2 ( max ( 1, min ( max_bs, n_samples ) ) ) ) )

        X_clean    = numpy.ascontiguousarray ( X, dtype = numpy.float32 )
        preds      = numpy.empty ( n_samples, dtype = numpy.float32 )

        with torch.inference_mode () :
            for i in range ( 0, n_samples, batch_size ) :
                chunk = X_clean [ i : i + batch_size ]
                b_x   = torch.from_numpy ( chunk ).to ( device )
                with torch.autocast ( device_type = device.type, enabled = is_cuda ) :
                    logits  = model ( b_x )
                    b_preds = torch.sigmoid ( logits )
                preds [ i : i + batch_size ] = b_preds.cpu ().numpy ()

        del X_clean 
        return preds.astype ( numpy.float32, copy = False )


# =============================================================================
## @class PyTorchAdamDensityReweighter
#  Density ratio reweighter training @c WeakMLP fold models with @c Adam.
#           Specializes @c PyTorchDensityReweighter to optimize each fold's
#           @c WeakMLP network over @c epochs using mini-batch stochastic gradient
#           descent via @c torch.optim.Adam and @c DataLoader with shuffling
class PyTorchAdamDensityReweighter ( PyTorchDensityReweighter ) :
    """ Density ratio reweighter training ``WeakMLP`` models with ``Adam``

    Iterates over mini-batches of ``(X_tr_t, y_tr_t, w_tr_t)`` and updates
    network parameters using ``torch.optim.Adam``
    """
    # =========================================================================
    ## Initialize the PyTorch Adam density ratio reweighter.
    #
    #           Sets default @c Adam learning rate, weight decay, epoch count,
    #           and batch size before delegating to @c PyTorchDensityReweighter
    #
    #  @param original Source feature matrix
    #  @param target Target feature matrix
    #  @param original_weight Optional source observation weights
    #  @param target_weight Optional target observation weights
    #  @param progress Optional progress flag for epoch loop display
    #  @param kwargs Additional keyword arguments passed to the base class
    def __init__ ( self                     , * ,
                   original                  ,
                   target                    ,
                   original_weight = None    ,
                   target_weight   = None    ,
                   progress        = True    , **kwargs ) :
        """ Initialize ``PyTorchAdamDensityReweighter`` with Adam defaults
        """
        config = {
            'batch_size'    : 'auto' ,
        }
        config.update ( kwargs )
        
        n_features = num_features ( target   ) 
        n_target   = num_samples  ( target   ) 
        n_original = num_samples  ( original ) 
        n_rows     = n_target + n_original
        
        # Honor explicit integer batch sizes (rounded down to a power of two),
        # 'full' batch mode, or fall back to Ostap's heuristic batch sizing
        batch_size = config.get ( 'batch_size', 'auto' )
        if isinstance ( batch_size, integer_types ) and 2 <= batch_size <= min ( n_rows, TORCH_MAX_BATCH_SIZE ) :
            if batch_size & ( batch_size - 1 ) : batch_size = 2 ** math.floor ( math.log2 ( batch_size ) )
        elif isinstance ( batch_size, string_types ) and batch_size.strip ().lower () == 'full' : batch_size = n_rows
        else : batch_size = optimal_batch_size ( n_rows, n_features, TORCH_MAX_BATCH_SIZE )
        
        config [ 'batch_size' ] = batch_size

        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   ,
                            progress        = progress        , **config )

    # =========================================================================
    ## Get the short symbol/tag representing the PyTorch reweighter method.
    #  @return String identifier symbol used for logging and console visual display
    @classprop
    def symbol ( self ) :
        """ Return short symbol/tag for PyTorch reweighter console display
        """
        return U.symbol_TORCH_ADAM

    # =========================================================================
    ## Get the full method identifier name.
    #  @return Method string tag used in summaries and warnings
    @property
    def method ( self ) :
        """ Return method identifier name for PyTorch iterative reweighter
        """
        return method_TORCH_ADAM

    # =========================================================================
    ## Train a fold @c WeakMLP model using mini-batch @c Adam.
    # 
    #           Creates a @c TensorDataset and shuffled @c DataLoader, runs
    #           @c epochs passes over the fold data, and applies @c optim.Adam
    #           updates using unnormalized sample weights (which @c DensityRatioLoss
    #           normalizes per class subset inside @c forward)
    #
    #  @param model @c WeakMLP instance to train in-place
    #  @param criterion @c DensityRatioLoss instance
    #  @param X_tr Contiguous float32 NumPy feature array (unused; tensors used)
    #  @param y_tr Float32 NumPy binary label array (unused; tensors used)
    #  @param w_tr Contiguous float32 NumPy weight array (unused; tensors used)
    #  @param X_tr_t Feature tensor on @c device
    #  @param y_tr_t Label tensor on @c device
    #  @param w_tr_t Weight tensor on @c device
    #  @param device Target @c torch.device
    def _fit_model ( self      ,
                     model     ,
                     criterion ,
                     X_tr      ,
                     y_tr      ,
                     w_tr      ,
                     X_tr_t    ,
                     y_tr_t    ,
                     w_tr_t    ,
                     device    ) :
        """ Train a single fold model using mini-batch ``Adam``
        """
        import torch.optim as optim
        from torch.utils.data import TensorDataset, DataLoader

        batch_size    = self.params.get ( 'batch_size'   , 256                 )
        nepochs       = self.params.get ( 'epochs'       , TORCH_EPOCHS        )
        learning_rate = self.params.get ( 'learning_rate', TORCH_LEARNING_RATE )
        weight_decay  = self.params.get ( 'weight_decay' , TORCH_WEIGHT_DECAY  )

        optimizer = optim.Adam ( model.parameters () ,
                                 lr           = learning_rate ,
                                 weight_decay = weight_decay  )
        
        dataset = TensorDataset ( X_tr_t, y_tr_t, w_tr_t )
        loader  = DataLoader    ( dataset, batch_size = batch_size, shuffle = True )

        model.train ()

        description = S.fixed_width ( train_symbol, 2 ) + S.fixed_width ( self.symbol, 2 ) + ':'
        the_bar     = progress_bar ( nepochs, description = description, silent = not self.progress_epochs )

        for epoch in safe_loop ( the_bar ) :
            for bx, by, bw in loader :
                optimizer.zero_grad ()
                preds = model ( bx )
                loss  = criterion ( preds, by, bw )
                loss.backward ()
                optimizer.step ()

# =============================================================================
## @class PyTorchLBFGSDensityReweighter
#  Density ratio reweighter training @c WeakMLP fold models with @c L-BFGS.
#
#           Specializes @c PyTorchDensityReweighter to optimize each fold's
#           @c WeakMLP network using @c torch.optim.LBFGS with strong Wolfe line
#           search.
#           When @c n_samples exceeds @c lbfgs_chunk, draws a
#           random subsample of size @c lbfgs_chunk without replacement so that
#           quasi-Newton memory and step cost stay bounded
class PyTorchLBFGSDensityReweighter ( PyTorchDensityReweighter ) :
    """ Density ratio reweighter training ``WeakMLP`` models with ``L-BFGS``

    Optimizes each fold's ``WeakMLP`` using ``torch.optim.LBFGS`` and
    ``strong_wolfe`` line search over the full fold or a chunk subsample
    """

    # =========================================================================
    ## Initialize the PyTorch L-BFGS density ratio reweighter.
    #
    #           Sets default @c lbfgs_chunk, @c learning_rate, @c weight_decay,
    #           and @c epochs (@c max_iter) before delegating to
    #           @c PyTorchDensityReweighter
    #
    #  @param original Source feature matrix
    #  @param target Target feature matrix
    #  @param original_weight Optional source observation weights
    #  @param target_weight Optional target observation weights
    #  @param progress Optional progress flag
    #  @param kwargs Additional keyword arguments passed to the base class
    def __init__ ( self                     , * ,
                   original                  ,
                   target                    ,
                   original_weight = None    ,
                   target_weight   = None    ,
                   progress        = True    , **kwargs ) :
        """ Initialize ``PyTorchLBFGSDensityReweighter`` with L-BFGS defaults
        """
        config = {
            'lbfgs_chunk'   : TORCH_LBFGS_CHUNK   ,
        }
        ## 
        config.update ( kwargs )

        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   ,
                            progress        = progress        , **config )

    # =========================================================================
    ## Get the short symbol/tag representing the PyTorch reweighter method.
    #  @return String identifier symbol used for logging and console visual display
    @classprop
    def symbol ( self ) :
        """ Return short symbol/tag for PyTorch reweighter console display
        """
        return U.symbol_TORCH_LBFGS

    # =========================================================================
    ## Get the full method identifier name.
    #  @return Method string tag used in summaries and warnings
    @property
    def method ( self ) :
        """ Return method identifier name for PyTorch iterative reweighter
        """
        return method_TORCH_LBFGS
    
    # =========================================================================
    ## Train a fold @c WeakMLP model using @c L-BFGS with strong Wolfe search.
    #           Prepares either the full fold tensors or a random chunk of size
    #           @c lbfgs_chunk, normalizes chunk weights, defines the @c closure
    #           evaluating @c criterion plus optional explicit L2 penalty, and
    #           runs @c optimizer.step(closure)
    #
    #  @param model @c WeakMLP instance to train in-place
    #  @param criterion @c DensityRatioLoss instance
    #  @param X_tr Contiguous float32 NumPy feature array of shape @c (n_samples, n_features)
    #  @param y_tr Float32 NumPy binary label array of shape @c (n_samples,)
    #  @param w_tr Contiguous float32 NumPy weight array of shape @c (n_samples,)
    #  @param X_tr_t Full-fold feature tensor on @c device
    #  @param y_tr_t Full-fold label tensor on @c device
    #  @param w_tr_t Full-fold weight tensor on @c device
    #  @param device Target @c torch.device
    def _fit_model ( self      ,
                     model     ,
                     criterion ,
                     X_tr      ,
                     y_tr      ,
                     w_tr      ,
                     X_tr_t    ,
                     y_tr_t    ,
                     w_tr_t    ,
                     device    ) :
        """ Train a single fold model using ``L-BFGS`` and ``strong_wolfe``
        """
        import torch
        import torch.optim as optim

        nepochs       = self.params.get ( 'epochs'       , TORCH_EPOCHS        )
        learning_rate = self.params.get ( 'learning_rate', TORCH_LEARNING_RATE )
        weight_decay  = self.params.get ( 'weight_decay' , TORCH_WEIGHT_DECAY  )
        lbfgs_chunk   = self.params.get ( 'lbfgs_chunk'  , TORCH_LBFGS_CHUNK   )

        n_samples = len ( X_tr )

        # Subsample without replacement when the fold exceeds lbfgs_chunk
        if lbfgs_chunk < n_samples :
            chunk_idx = numpy.random.choice ( n_samples, size = lbfgs_chunk, replace = False )
            X_ch      = X_tr [ chunk_idx ]
            y_ch      = y_tr [ chunk_idx ]
            w_ch      = w_tr [ chunk_idx ]
            w_ch_norm = w_ch / ( numpy.sum ( w_ch ) + 1e-12 )

            X_ch_t = torch.from_numpy ( X_ch      ).to ( device )
            y_ch_t = torch.from_numpy ( y_ch      ).to ( device )
            w_ch_t = torch.from_numpy ( w_ch_norm ).to ( device )
        else :
            X_ch_t = X_tr_t
            y_ch_t = y_tr_t
            w_ch_t = w_tr_t

        optimizer = optim.LBFGS ( model.parameters ()      ,
                                  lr             = learning_rate ,
                                  max_iter       = nepochs       ,
                                  line_search_fn = 'strong_wolfe' )

        ## Zero gradients, evaluate density-ratio + L2 loss, and backpropagate
        def closure () :
            """ Zero gradients, evaluate density-ratio + L2 loss, and backpropagate
            """
            optimizer.zero_grad ()
            preds = model ( X_ch_t )
            loss  = criterion ( preds, y_ch_t, w_ch_t )
            if 0 < weight_decay :
                l2_reg = sum ( torch.sum ( p ** 2 ) for p in model.parameters () )
                loss   = loss + 0.5 * weight_decay * l2_reg
            loss.backward ()
            return loss

        model.train ()
        optimizer.step ( closure )

        if lbfgs_chunk < n_samples :
            del X_ch_t, y_ch_t, w_ch_t

# =============================================================================
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
#  Density ratio reweighter using Logistic Regression as the underlying classifier
#
#  With \f$d\f$ polynomial features  and \$n$\ features, number of "new"
# features id  \f$ C^{n+d}_{d} = \left( \begin{array}{c} n+d \\ d \end{array}\right) \f$
# 
class LogRegressionDensityReweighter ( DensityReweighter ) :
    """ Density ratio reweighter using Logistic Regression as the underlying classifier.
    """

    def __init__ ( self , * , 
                  original               ,
                  target                 ,
                  original_weight        = None ,
                  target_weight          = None ,
                  store_original_weights = True , **params ) :        
        """ Initialize Logistic Regression density reweighter.
        """
        config = { 'C'           : 1.0     ,
                   'solver'      : 'lbfgs' ,
                   'max_iter'    : 2000    ,
                   'polynomials' : 2       }
        
        config.update ( params )

        if 'polynomial' in config : config [ 'polynomials' ] = config.pop ( 'polynomial' )
        
        poly       = config.get   ( 'polynomials', 2 )
        n_features = num_features ( original )
        
        if 15 <= n_features and 1 < poly : config [ 'polynomials' ] = 1
        
        super ().__init__ ( original               = original               ,
                            target                 = target                 ,
                            original_weight        = original_weight        ,
                            target_weight          = target_weight          ,
                            store_original_weights = store_original_weights , **config )
        
    # =========================================================================
    ## Short symbol or tag representing the LightGBM reweighter method.
    #  @return String identifier/symbol used for logging and visual presentation.
    @classprop
    def symbol ( self ) :
        """ Short symbol/tag for reweighter.
        """
        return U.symbol_LOGREG

    # =========================================================================
    ## Return method identifier name. 
    @property
    def method ( self ) :
        """ Return method identifier name.
        """
        return method_LOGREG

    # =========================================================================
    ## Train single Logistic Regression model on fold data.
    #  @param X_train Training features array (N-1)/N.
    #  @param y_train Training binary labels.
    #  @param w_train Training event weights or None.
    #  @return Trained Scikit-Learn Pipeline instance.
    def _train_single_model ( self, X_train, y_train, w_train ) :
        """ Train single Logistic Regression model on fold data
        """
        from sklearn.linear_model  import LogisticRegression
        from sklearn.pipeline      import Pipeline        
        from sklearn.preprocessing import StandardScaler, PolynomialFeatures
        import pickle
        
        params = {}
        params.update ( self.params )
        params.pop ( 'n_jobs', None ) 

        poly       = params.pop ( 'polynomials', 2 )
        n_features = num_features ( X_train )
        if 20 <= n_features and 1 < poly : poly = 1
        
        lr_kwargs = { k: v for k, v in params.items () if k in valid_LR_params }
        steps     = [ ( 'scaler1', StandardScaler () ) ] 

        if poly and isinstance ( poly, int ) and 1 < poly <= 5 :
            steps += [ ( 'poly'   , PolynomialFeatures ( degree = poly, include_bias = False ) ) ]
            steps += [ ( 'scaler2', StandardScaler () ) ]
            
        lr_model = LogisticRegression ( **lr_kwargs )
        steps   += [ ( 'logistic', lr_model ) ]
        model    = Pipeline ( steps )

        w_tr = numpy.ascontiguousarray ( w_train, dtype = numpy.float32 ) if w_train is not None else None

        if w_tr is not None : model.fit ( X_train, y_train, logistic__sample_weight = w_tr )
        else                : model.fit ( X_train, y_train )

        if self.phoenix :            
            model_bytes = pickle.dumps ( model )
            del model
            model       = pickle.loads ( model_bytes )            
            del model_bytes 
        
        return model

    # =========================================================================
    ## Predict probabilities using a Logistic Regression model.
    #  @param model Trained Logistic Regression pipeline.
    #  @param X Input features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model ( self, model, X ) :
        """ Predict probabilities using a Logistic Regression model.
        """
        X_clean = numpy.ascontiguousarray ( X, dtype = numpy.float32 )
        p       = model.predict_proba ( X_clean ) [ :, 1 ]
        return p.astype ( numpy.float32, copy = False )

# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__, logger = logger )

    from ostap.stats.tools import ( hasLightGBM ,
                                    hasXGBoost  ,
                                    hasCatBoost ,
                                    hasPyTorch  ,
                                    hasSkLearn  )

    if not hasLightGBM ( False ) : logger.warning ( "No LightGBM available!" ) 
    if not hasXGBoost  ( False ) : logger.warning ( "No XGBoost  available!" ) 
    if not hasCatBoost ( False ) : logger.warning ( "No CatBoost available!" ) 
    if not hasPyTorch  ( False ) : logger.warning ( "No PyTorch  available!" ) 
    if not hasSkLearn  ( False ) : logger.warning ( "No SkLearn  available!" ) 

# =============================================================================
#                                                                       The END 
# =============================================================================
