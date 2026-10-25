
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
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.tools.dreweighters' )
else                       : logger = getLogger ( __name__ )
# =============================================================================
has_sklearn = hasSkLearn () 
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

CATB_DEFAULT_ESTIMATORS        = 1600 # High capacity: deep ensemble with low learning rate
CATB_MAX_DEPTH                 =   10 # Deep tree capacity for high dimensions / large statistics
CATB_LEARNING_RATE             = 0.005 # Small step for stable density ratio convergence

# =============================================================================
## @class LightGBMDensityReweighter
#  Density ratio reweighter using LightGBM as the underlying classifier.
class LightGBMDensityReweighter ( DensityReweighter ) :
    """ Density ratio reweighter using LightGBM as the underlying classifier. """

    ## @brief Initialize LightGBM density reweighter.
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
        """ Initialize LightGBM density reweighter (High Capacity Mode Defaults). """
        config = {
            'objective'          : 'binary'                ,
            'metric'             : 'binary_logloss'        ,
            'n_estimators'       : LGBM_DEFAULT_ESTIMATORS ,
            'learning_rate'      : LGBM_LEARNING_RATE      ,
            'max_depth'          : LGBM_MAX_DEPTH          ,
            'num_leaves'         : min ( 1023 , 2 ** LGBM_MAX_DEPTH - 1 ) , 
            'max_bin'            : 2047                    , 
            'boost_from_average' : True                    ,
            'verbosity'          : -1                      ,
            'n_jobs'             : -1                      ,
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

    ## @property method
    #  @brief Return the method identifier name.
    @property
    def method ( self ) :
        """ Return the method identifier name. """
        return method_LGBM 

    ## @brief Train single LightGBM model on fold data.
    #  @param X_train Training features array (N-1)/N.
    #  @param y_train Training binary labels.
    #  @param w_train Training event weights or None.
    #  @return Trained LightGBM booster instance.
    def _train_single_model ( self, X_train, y_train, w_train ) : 
        """ Train single LightGBM model on fold data. """
        import lightgbm as LightGBM
        import gc
        
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
        model_bytes = model.model_to_string ()
        best_iter   = getattr ( model, 'best_iteration_', None ) or getattr ( model, 'best_iteration', None )                    
        model.free_dataset () 
        del model, trn_data, X_tr, w_tr
        
        model = LightGBM.Booster ( model_str = model_bytes )
        if best_iter is not None : model.best_iteration = best_iter 
        del model_bytes
        
        gc.collect ()                

        return model

    ## @brief Predict probabilities using a LightGBM model.
    #  @param model Trained LightGBM booster instance.
    #  @param X Features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model ( self, model, X ) :
        """ Predict probabilities using a LightGBM model. """
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
    """ Density ratio reweighter using XGBoost as the underlying classifier. """

    ## @brief Initialize XGBoost density reweighter.
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
        """ Initialize XGBoost density reweighter (High Capacity Mode Defaults). """
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
        
        super ().__init__ (
            original        = original        ,
            target          = target          ,
            original_weight = original_weight ,
            target_weight   = target_weight   ,
            **config
        )

    ## @property method
    #  @brief Return the method identifier name.
    @property
    def method ( self ) :
        """ Return the method identifier name. """
        return method_XGB 

    ## @brief Train single XGBoost booster on fold data.
    #  @param X_train Training features array (N-1)/N.
    #  @param y_train Training binary labels.
    #  @param w_train Training event weights or None.
    #  @return Trained XGBoost Booster instance.
    def _train_single_model ( self, X_train, y_train, w_train ) :
        """ Train single XGBoost booster on fold data. """
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
        model_bytes = model.save_raw ( raw_format = "json" )
        best_iter   = getattr ( model, 'best_iteration_', None ) or getattr ( model, 'best_iteration', None )                    
        
        del dtrain, model, X_tr, y_tr, w_tr
        
        model = XGBoost.Booster ()
        model.load_model ( model_bytes )
        if best_iter is not None : model.best_iteration = best_iter
        
        del model_bytes            
        gc.collect ()

        return model
    
    ## @brief Predict probabilities using an XGBoost model.
    #  @param model Trained XGBoost Booster.
    #  @param X Input features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model ( self, model, X ) :
        """ Predict probabilities using an XGBoost model. """
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
    """ Density ratio reweighter using CatBoost as the underlying classifier. """
    
    def __init__ ( self                   , * , 
                   original               ,
                   target                 ,
                   original_weight        = None ,
                   target_weight          = None ,
                   store_original_weights = True , **params ) :
        """ Initialize CatBoost density reweighter. """
        
        from catboost.utils import get_gpu_device_count
        gpu        = get_gpu_device_count () 
        self.__gpu = gpu 
        
        config = { 'loss_function'    : 'Logloss'                        ,
                   'eval_metric'      : 'Logloss'                        ,
                   'iterations'       : CATB_DEFAULT_ESTIMATORS          ,
                   'learning_rate'    : CATB_LEARNING_RATE               ,
                   'depth'            : CATB_MAX_DEPTH                   ,
                   'l2_leaf_reg'      : 0.5                              ,
                   'min_data_in_leaf' : 5                                ,
                   'bootstrap_type'   : 'MVS'                            ,
                   'mvs_reg'          : 2.0                              , 
                   'random_strength'  : 0.01                             ,
                   'boosting_type'    : 'Plain'                          ,
                   'device'           : 'GPU' if 0 < self.gpu else 'CPU' , 
                   'verbose'          : False                            }
        
        config.update ( params )
        
        if 'min_child_samples' in config : config [ 'min_data_in_leaf' ] = config.pop ( 'min_child_samples' )
        if 'n_jobs'            in config : config [ 'thread_count'     ] = config.pop ( 'n_jobs'            )
        if 'n_estimators'      in config : config [ 'iterations'       ] = config.pop ( 'n_estimators'      )
        
        if self.gpu <= 0 : config [ 'device' ] = 'CPU'
        
        super ().__init__ ( original               = original               ,
                            target                 = target                 ,
                            original_weight        = original_weight        ,
                            target_weight          = target_weight          ,
                            store_original_weights = store_original_weights , **config ) 
        
    @property
    def gpu ( self ) :
        """ `gpu`: result of catboost.utils.get_gpu_device_count """
        return self.__gpu
    
    @property
    def method ( self ) :
        """ Return method name identifier. """
        return method_CATB
    
    ## @brief Train single CatBoost model on fold data.
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
        
        params [ 'thread_count'  ] = params.get ( 'thread_count', params.pop ( 'n_jobs', max ( 1, numcpu () // 2 ) ) ) 
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
        model_bytes = pickle.dumps ( model ) 
        best_iter   = getattr ( model, 'best_iteration_', None ) or getattr ( model, 'best_iteration', None )                    
        
        del trn_pool, model, X_tr, y_tr, w_tr
        
        model = pickle.loads ( model_bytes )
        if best_iter is not None : model.best_iteration = best_iter
        
        del model_bytes 
        gc.collect ()
        
        return model

    ## @brief Predict probabilities using a CatBoost model.
    #  @param model Trained CatBoost model.
    #  @param X Features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model ( self, model, X ) :
        """ Predict probabilities using a CatBoost model. """
        X_clean   = numpy.ascontiguousarray ( X, dtype = numpy.float32 )
        best_iter = getattr ( model, 'best_iteration_', None ) or getattr ( model, 'best_iteration', None )        
        kwargs    = {}
        if best_iter is not None and 0 <= best_iter : kwargs [ 'ntree_end' ] = best_iter + 1
        thread_count = self.params.get ( 'thread_count', self.params.get ( 'n_jobs', max ( 1, numcpu () // 2 ) ) )                                    
        p            = model.predict_proba ( X_clean, thread_count = thread_count, **kwargs ) [ :, 1 ]
        return p.astype ( numpy.float32, copy = False )

# =============================================================================
## @class PyTorchDensityReweighter
#  Optimized PyTorch MLP density ratio reweighter.
class PyTorchDensityReweighter ( DensityReweighter ) :
    """ Highly optimized and stabilized PyTorch MLP density ratio reweighter. """

    def __init__ ( self                   , * ,
                   original               ,
                   target                 ,
                   original_weight = None ,
                   target_weight   = None ,
                   progress        = True , **kwargs ) :
        """ Constructor. """
        import torch
        cuda = torch.cuda.is_available () 

        config = { 'hidden_dims'   : ( 256, 256, 256 )         ,
                   'dropout'       : 0.35                      ,               
                   'learning_rate' : 2e-4                      ,              
                   'weight_decay'  : 1e-5                      ,              
                   'batch_size'    : 256                       ,
                   'epochs'        : 2000                      ,         
                   'patience'      :   50                      ,
                   'min_epochs'    :   50                      ,
                   'min_delta'     : 1e-5                      , 
                   'eval_freq'     : 5                         , 
                   'max_logit'     : 8.0                       , 
                   'logit_reg'     : 1e-4                      , 
                   'grad_clip'     : 0.2                       ,
                   'ema_decay'     : 0.999                     , 
                   'device'        : 'cuda' if cuda else 'cpu' ,
                   'compile'       : False                     ,
                   'jit'           : True                      ,
                   'n_jobs'        : -1                        }

        config.update ( kwargs )
        self.__progress_epochs = config.pop ( 'progress', True )
        
        device = config.get ( 'device', 'cuda' if cuda else 'cpu' )        
        if not cuda : config [ 'device' ] = 'cpu'

        batch_size = config.get ( 'batch_size', None )
        n          = num_samples ( target ) + num_samples ( original )
        
        if   isinstance ( batch_size, string_types  ) and batch_size.lower () == 'auto' : batch_size = None
        elif isinstance ( batch_size, integer_types ) and not 1 <= batch_size <= n      : batch_size = None
        
        if batch_size is None : 
            batch_size = math.floor ( math.log2 ( batch_size ) )
            batch_size = 2 ** batch_size  
            config [ 'batch_size' ] = batch_size 
            
        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   ,
                            progress        = False           , **config )

        torch.set_num_threads ( 1 )

    @property
    def method ( self ) :
        """ Return method identifier name. """
        return method_TORCH 

    ## @brief Trains a single neural network model on fold data.
    #  @param X_train Features matrix for training (N-1)/N.
    #  @param y_train Binary labels for training fold.
    #  @param w_train Sample weights for training fold.
    #  @return Trained PyTorch model instance.
    def _train_single_model ( self, X_train, y_train, w_train ) :
        """ Trains a single neural network model on fold data. """
        import torch
        import torch.nn            as nn
        import torch.nn.functional as F
        from   torch.optim.swa_utils import AveragedModel, get_ema_multi_avg_fn
        from   sklearn.preprocessing import StandardScaler 
        import time

        n_threads = self.params.get ( 'n_jobs', max ( 1, numcpu () // 2 ) )
        n_threads = min ( 4, n_threads, max ( 1, numcpu () // 2 ) )
        torch.set_num_threads ( max ( 1, n_threads ) )

        device  = torch.device ( self.params.get ( 'device', 'cuda' if torch.cuda.is_available () else 'cpu' ) )
        is_cuda = device.type == 'cuda'
        
        n_samples  = num_samples  ( X_train )
        n_features = num_features ( X_train )
        
        batch_size = self.params.get ( 'batch_size', 256 )
        if not isinstance ( batch_size, int ) or batch_size <= 1 : batch_size = 256
            
        batch_size = min ( max ( 2, ( n_samples + 1 ) // 2 ), batch_size ) 
        batch_size = max ( 2, 2 ** math.floor ( math.log2 ( batch_size ) ) )

        nepochs       = self.params.get ( 'epochs'       , 2000          )
        max_logit     = self.params.get ( 'max_logit'    , 8.0           )
        logit_reg     = self.params.get ( 'logit_reg'    , 1e-4          )
        grad_clip     = self.params.get ( 'grad_clip'    , 0.2           )
        ema_decay     = self.params.get ( 'ema_decay'    , 0.999         )
        learning_rate = self.params.get ( 'learning_rate', 2e-4          )
        dropout       = self.params.get ( 'dropout'      , 0.35          )
        hidden_dims   = self.params.get ( 'hidden_dims'  , ( 256, 256, 256 ) )
        weight_decay  = self.params.get ( 'weight_decay' , 1e-5          )

        scaler      = StandardScaler () 
        X_tr_scaled = scaler.fit_transform ( numpy.ascontiguousarray ( X_train, dtype = numpy.float32 ) ).astype ( numpy.float32, copy = False )
        w_tr        = numpy.ascontiguousarray ( w_train, dtype = numpy.float32 ) if w_train is not None else numpy.ones ( n_samples, dtype = numpy.float32 )
        y_tr        = y_train.astype ( numpy.float32 )

        X_tr_t = torch.from_numpy ( X_tr_scaled ).to ( device )
        y_tr_t = torch.from_numpy ( y_tr        ).to ( device )
        w_tr_t = torch.from_numpy ( w_tr        ).to ( device )

        class DensityRatioLoss ( nn.Module ) :
            def __init__ ( self, max_logit = 8.0, logit_reg = 1e-4 ) :
                super ().__init__ ()
                self.max_logit = max_logit
                self.logit_reg = logit_reg

            def forward ( self, logits, targets, weights ) :
                logits     = logits.float ()
                targets    = targets.float ()
                weights    = weights.float ()
                
                mask_orig  = ( targets == 0 )
                mask_targ  = ( targets == 1 )
                
                z_orig     = logits [ mask_orig ]
                w_orig     = weights[ mask_orig ]
                z_targ     = logits [ mask_targ ]
                w_targ     = weights[ mask_targ ]
                
                w_orig_sum = w_orig.abs ().sum () + 1e-8
                w_targ_sum = w_targ.abs ().sum () + 1e-8
                
                exp_z      = torch.exp ( torch.clamp ( z_orig, max = self.max_logit ) )
                loss_orig  = torch.sum ( w_orig * exp_z ) / w_orig_sum
                loss_targ  = torch.sum ( w_targ * z_targ ) / w_targ_sum
                
                loss       = loss_orig - loss_targ
                if 0 < self.logit_reg : loss += self.logit_reg * torch.mean ( logits ** 2 )
                return loss

        class ResBlock ( nn.Module ) :
            def __init__ ( self, dim, dropout = 0.1 ) :
                super ().__init__ ()
                layers = [ nn.Linear    ( dim, dim ) ,
                           nn.LayerNorm ( dim )       ,
                           nn.SiLU      ()            ] 
                if 0 < dropout < 1 : layers.append ( nn.Dropout ( p = dropout ) ) 
                layers.append ( nn.Linear    ( dim, dim ) )
                layers.append ( nn.LayerNorm ( dim )       ) 
                self.block = nn.Sequential ( *layers ) 
                self.act   = nn.SiLU ()

            def forward ( self, x ) :
                return self.act ( x + self.block ( x ) )

        class TabularResNet ( nn.Module ) :
            def __init__ ( self, n_features, hidden_dims = ( 256, 256 ), dropout = 0.1 ) :
                super ().__init__ ()
                layers    = []
                first_dim = hidden_dims [ 0 ]
                layers.append ( nn.Linear ( n_features, first_dim ) )
                in_dim    = first_dim
                for h_dim in hidden_dims :
                    if in_dim != h_dim :
                        layers.append ( nn.Linear ( in_dim, h_dim ) )
                        in_dim = h_dim
                    layers.append ( ResBlock ( h_dim, dropout = dropout ) )
                if 0.0 < dropout < 1.0 : layers.append ( nn.Dropout ( p = dropout ) )
                layers.append ( nn.Linear ( in_dim, 1 ) )
                self.net = nn.Sequential ( *layers )
                
            def forward ( self, x ) : 
                return self.net ( x ).squeeze ( -1 )

        model_config  = { 'n_features': n_features, 'hidden_dims': hidden_dims, 'dropout': dropout }
        model         = TabularResNet ( **model_config ).to ( device )
        ema_model     = AveragedModel ( model, multi_avg_fn = get_ema_multi_avg_fn ( ema_decay ) )

        optimizer     = torch.optim.AdamW                          ( model.parameters (), lr = learning_rate, weight_decay = weight_decay )
        scheduler     = torch.optim.lr_scheduler.CosineAnnealingLR ( optimizer, T_max = nepochs, eta_min = 0.1 * learning_rate )
        scaler_amp    = torch.amp.GradScaler                       ( 'cuda', enabled = is_cuda )
        criterion_drw = DensityRatioLoss                           ( max_logit = max_logit, logit_reg = logit_reg )

        use_autocast  = is_cuda

        for epoch in range ( nepochs ) :
            model.train ()
            perm = torch.randperm ( n_samples, device = device )
            X_sh = X_tr_t [ perm ]
            y_sh = y_tr_t [ perm ]
            w_sh = w_tr_t [ perm ]
            
            for i in range ( 0, n_samples, batch_size ) :
                b_x = X_sh [ i : i + batch_size ]
                b_y = y_sh [ i : i + batch_size ]
                b_w = w_sh [ i : i + batch_size ]
            
                optimizer.zero_grad ( set_to_none = True )
                with torch.autocast ( device_type = device.type, enabled = use_autocast ) :
                    logits = model ( b_x )
                    loss   = criterion_drw ( logits, b_y, b_w )
                        
                if is_cuda :
                    scaler_amp.scale ( loss ).backward ()
                    scaler_amp.unscale_ ( optimizer )
                    if 0 < grad_clip : torch.nn.utils.clip_grad_norm_ ( model.parameters (), max_norm = grad_clip )
                    scaler_amp.step   ( optimizer )
                    scaler_amp.update ()
                else :
                    loss.backward ()
                    if 0 < grad_clip : torch.nn.utils.clip_grad_norm_ ( model.parameters (), max_norm = grad_clip )
                    optimizer.step ()
                        
            ema_model.update_parameters ( model )
            scheduler.step ()

        unwrapped_ema = ema_model.module if hasattr ( ema_model, 'module' ) else ema_model
        best_state    = { k : v.clone ().detach () for k, v in unwrapped_ema.state_dict ().items () }

        del model, ema_model, optimizer, scaler_amp, scheduler, X_tr_t, y_tr_t, w_tr_t
        gc.collect ()
        if is_cuda : torch.cuda.empty_cache ()

        final_model = TabularResNet ( **model_config ).to ( device )
        final_model.load_state_dict ( best_state )
        final_model.eval ()
        final_model.scaler = scaler

        return final_model

    ## @brief Predicts probabilities for the given dataset.
    #  @param model Trained PyTorch model.
    #  @param X Features array to predict on.
    #  @return Numpy array of predicted probabilities.
    def _predict_single_model ( self, model, X ) :
        """ Predicts probabilities for the given dataset. """
        import torch
        model.eval ()
        device    = next ( model.parameters () ).device
        is_cuda   = device.type == 'cuda'
        n_samples = num_samples ( X )
        
        batch_size = 2**17 if is_cuda else 2**16
        batch_size = max ( 1024, 2 ** math.floor ( math.log2 ( max ( 1, min ( batch_size, n_samples ) ) ) ) )

        X_scaled = model.scaler.transform ( numpy.ascontiguousarray ( X, dtype = numpy.float32 ) ).astype ( numpy.float32, copy = False )        
        preds    = numpy.empty ( n_samples, dtype = numpy.float32 )

        with torch.inference_mode () :
            for i in range ( 0, n_samples, batch_size ) :
                chunk   = X_scaled [ i : i + batch_size ]
                b_x     = torch.from_numpy ( chunk ).to ( device )
                with torch.autocast ( device_type = device.type, enabled = is_cuda ) :
                    logits  = model ( b_x )
                    b_preds = torch.sigmoid ( logits )
                preds [ i : i + batch_size ] = b_preds.cpu ().numpy ()

        return preds

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
#  Density ratio reweighter using Logistic Regression as the underlying classifier.
class LogRegressionDensityReweighter ( DensityReweighter ) :
    """ Density ratio reweighter using Logistic Regression as the underlying classifier. """

    def __init__ ( self , * , 
                  original               ,
                  target                 ,
                  original_weight        = None ,
                  target_weight          = None ,
                  store_original_weights = True , **params ) :        
        """ Initialize Logistic Regression density reweighter. """
        config = { 'C'           : 1.0     ,
                   'solver'      : 'lbfgs' ,
                   'max_iter'    : 2000    ,
                   'polynomials' : 2       }
        
        config.update ( params )
        
        poly       = config.get   ( 'polynomials', 2 )
        n_features = num_features ( original )
        
        if 15 <= n_features and 1 < poly :
            config [ 'polynomials' ] = 1
        
        super ().__init__ ( original               = original               ,
                            target                 = target                 ,
                            original_weight        = original_weight        ,
                            target_weight          = target_weight          ,
                            store_original_weights = store_original_weights , **config )
        
    @property
    def method ( self ) :
        """ Return method identifier name. """
        return method_LR 

    ## @brief Train single Logistic Regression model on fold data.
    #  @param X_train Training features array (N-1)/N.
    #  @param y_train Training binary labels.
    #  @param w_train Training event weights or None.
    #  @return Trained Scikit-Learn Pipeline instance.
    def _train_single_model ( self, X_train, y_train, w_train ) :
        """ Train single Logistic Regression model on fold data. """
        from sklearn.linear_model  import LogisticRegression
        from sklearn.pipeline      import Pipeline        
        from sklearn.preprocessing import StandardScaler, PolynomialFeatures
        import pickle
        import gc
        
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

        model_bytes = pickle.dumps ( model ) 
        del model, w_tr
        
        model = pickle.loads ( model_bytes )
        del model_bytes 
        gc.collect ()
        
        return model

    ## @brief Predict probabilities using a Logistic Regression model.
    #  @param model Trained Logistic Regression pipeline.
    #  @param X Input features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model ( self, model, X ) :
        """ Predict probabilities using a Logistic Regression model. """
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