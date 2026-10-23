#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/reweghters.py
#  Iterative Density rreweighters based on LightGBM, XGBoost, CatBoost & PyTorch 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @Date  2026-07-18
# =============================================================================
""" Iterative Density ratio reweighters based on
 - LightGBM
 - XGBoost
 - CatBoost
 - PyTorch 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    'LightGBMIterativeReweighter'  , ## LightGBM-based iterative density reweighter
    'XGBoostIterativeReweighter'   , ## XGBoost-based  iterative density reweighter
    'CatBoostIterativeReweighter'  , ## CatBoost-based iterative density reweighter
    'PyTorchIterativeReweighter'   , ## PyTorch-based iterative density reweighter
) 
# =============================================================================
from   ostap.core.ostap_types       import num_types, string_types, integer_types  
from   ostap.utils.core             import typename
from   ostap.utils.basic            import numcpu, num_jobs, NoContext
from   ostap.logger.pretty          import nice_print 
from   ostap.utils.progress_bar     import progress_bar, text_bar
from   ostap.stats.utils            import num_samples
from   ostap.tools.base_reweighters import IterativeReweighter  
import ostap.logger.symbols         as     S
import numpy, abc, math, gc  
# =============================================================================
# Logging setup
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighters' )
else                       : logger = getLogger( __name__ )
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
LGBM_ESTIMATORS       =  300 # High capacity: deep ensemble with low learning rate
LGBM_MAX_DEPTH        =    6 # Deep tree capacity for high dimensions / large statistics
LGBM_LEARNING_RATE    = 0.05 # Small step for stable density ratio convergence

XGB_ESTIMATORS        =  300 # High capacity: deep ensemble with low learning rate
XGB_MAX_DEPTH         =    6 # Deep tree capacity for high dimensions / large statistics
XGB_LEARNING_RATE     = 0.05 # Small step for stable density ratio convergence

CATB_ESTIMATORS       =  300 # High capacity: deep ensemble with low learning rate
CATB_MAX_DEPTH        =    6 # Deep tree capacity for high dimensions / large statistics
CATB_LEARNING_RATE    = 0.10 # Small step for stable density ratio convergence

TORCH_ESTIMATORS      =  500 # Number of boosting steps (epochs/models)
TORCH_LEARNING_RATE   = 0.02 # Step size for weight boosting
TORCH_HIDDEN_DIM      =  128 # Hidden layer size for small MLP weak learners
TORCH_EPOCHS_PER_STEP =  100 # Training epochs per boosting step

# =============================================================================
## LigthGBM
# =============================================================================

# =============================================================================
## @class LightGBMIterativeReweighter
#  Iterative weight boosting density reweighter using LightGBM.
class LightGBMIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using LightGBM.
    """
    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        
        config = {
            ## 
            'n_estimators'            :  LGBM_ESTIMATORS    ,
            'learning_rate'           :  LGBM_LEARNING_RATE ,
            'max_depth'               :  LGBM_MAX_DEPTH     ,
            ## 
            'max_bins'                :  1023 ,
            'min_data_in_leaf'        :  30   ,
            'min_sum_hessian_in_leaf' :  1e-2 , 
            ## 
            'bagging_fraction'        :  0.9  , 
            'bagging_freq'            :  1    ,
            'lambda_l1'               :  0    ,
            'lambda_l2'               :  0    ,                        
            ## 
            'verbosity'               : -1    ,
            'n_jobs'                  :  1    , ## we are working with a single tree 
        }
        
        config.update ( kwargs )
        if 'silent' in config and config.get ( 'silent' ) : config [ 'verbosity' ] = -1        
        if 'num_boost_round'  in config : config [ 'n_estimators'            ] = config.pop ( 'num_boost_round'  , LGBM_ESTIMATORS ) 
        if 'min_child_weight' in config : config [ 'min_sum_hessian_in_leaf' ] = config.pop ( 'min_child_weight' , 1e-2 )         
        
        super().__init__ ( original        = original        ,
                           target          = target          ,
                           original_weight = original_weight ,
                           target_weight   = target_weight   , **config )

    # =========================================================================
    @property
    def method ( self ) :
        return 'DRW/LGBM-Iterative'

    # =========================================================================
    def _fit_iterative_stream ( self               ,
                                X_orig             ,
                                w_orig             ,
                                X_targ             ,
                                w_targ             ) :
        
        import lightgbm as LightGBM

        n_orig = len ( X_orig )
        n_targ = len ( X_targ )

        X_comb = numpy.vstack ( [ X_orig, X_targ ] )
        
        abs_w_orig = numpy.abs ( w_orig )
        abs_w_targ = numpy.abs ( w_targ )
        
        sum_w_orig = numpy.sum ( abs_w_orig )
        sum_w_targ = numpy.sum ( abs_w_targ )
        
        norm_targ = w_targ * ( sum_w_orig / ( sum_w_targ + 1e-12 ) )

        y_comb = numpy.hstack ( [ numpy.zeros ( n_orig , dtype = numpy.float32 ) ,
                                  numpy.ones  ( n_targ , dtype = numpy.float32 ) ] )

        current_w_orig = abs_w_orig.copy ()
        
        n_estimators  = self.params.get ( 'n_estimators'  , LGBM_ESTIMATORS    )
        learning_rate = self.params.get ( 'learning_rate' , LGBM_LEARNING_RATE )
        max_depth     = self.params.get ( 'max_depth'     , LGBM_MAX_DEPTH     )
        num_leaves    = self.params.get ( 'num_leaves'    , 2 ** max_depth - 1 )
    
        models = []

        params = {
            'objective'               : self.params.get ( 'objective'                , 'binary'         ) ,     
            'metric'                  : self.params.get ( 'metric'                   , 'binary_logloss' ) , 
            ## 
            'max_bins'                : self.params.get ( 'max_bins'                 , 1023 ) , 
            'min_data_in_leaf'        : self.params.get ( 'min_data_in_leaf'         , 30   ) , 
            'min_sum_hessian_in_leaf' : self.params.get ( 'min_sum_hessian_in_leaf'  , 1e-2 ) , 
            ##
            'feature_fraction'        : self.params.get ( 'feature_fraction'         , 1.0  ) , 
            'bagging_fraction'        : self.params.get ( 'bagging_fraction'         , 0.9  ) , 
            'bagging_freq'            : self.params.get ( 'bagging_freq'             , 1    ) , 
            ## 
            'lambda_l1'               : self.params.get ( 'lambda_l1'                , 0    ) , 
            'lambda_l2'               : self.params.get ( 'lambda_l2'                , 0    ) , 
            ##
            'max_depth'               : max_depth  ,
            'num_leaves'              : num_leaves ,
            ## 
            'learning_rate'           : 1.0               , ## ATTENTION!!! 
            'verbosity'               : self.params.get ( 'verbosity' , -1 ) , 
            'n_jobs'                  : self.params.get ( 'n_jobs'    ,  1 ) , 
        }

        eps = 1e-4 
        for it in progress_bar ( n_estimators , description = 'Boosting:' , silent = not self.progress ) :
            
            comb_weights  = numpy.hstack ( [ current_w_orig , numpy.abs ( norm_targ ) ] )
            comb_weights /= ( numpy.sum ( comb_weights ) + 1e-12 )

            trn_data = LightGBM.Dataset ( X_comb , label = y_comb , weight = comb_weights, free_raw_data = True )            
            model    = LightGBM.train ( params          = params   ,
                                        train_set       = trn_data ,
                                        num_boost_round = 1        )

            p_orig          = model.predict ( X_orig )            
            p_orig          = numpy.clip    ( p_orig , eps , 1.0 - eps )
            r_step          = p_orig / ( 1.0 - p_orig )
            current_w_orig *= numpy.power ( numpy.clip ( r_step , 1e-2 , 1e+2 ) , learning_rate )

            # =================================================================
            if True : # =======================================================
                # =============================================================
                model_bytes = model.model_to_string ()
                model.free_dataset ()
                del model
                del trn_data
                ## recreate model 
                model = LightGBM.Booster ( model_str = model_bytes )
                del model_bytes
            
            models.append ( model )
            
        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_iter_clip ( final_ratios , w_orig , norm_targ )

        return models, final_ratios

    # =========================================================================
    def _predict_iterative_stream ( self, models, X ) :
        
        learning_rate = self.params.get ( 'learning_rate' , 0.05 )
        log_ratios    = numpy.zeros ( len ( X ) , dtype = numpy.float32 )

        for model in models :
            p           = model.predict ( X )
            p           = numpy.clip ( p , 1e-4 , 1.0 - 1e-4 )
            r_step      = p / ( 1.0 - p )
            log_ratios += numpy.log ( numpy.clip ( r_step , 1e-2 , 1e+2 ) ) * learning_rate

        return numpy.exp ( log_ratios )

    
# =============================================================================
## XGBoost 
# =============================================================================

# =============================================================================
# True Density Ratio Objective for XGBoost (Exact Exponential DRE)
# =============================================================================
def xgboost_density_ratio_obj ( preds , dtrain ) :
    """ Exact Density Ratio Estimation objective for XGBoost. """
    labels  = dtrain.get_label  ()
    weights = dtrain.get_weight ()
    
    if weights is None or len ( weights ) == 0:
        weights = numpy.ones_like ( labels )
        
    # Защита от зануления градиентов (baseline weight)
    weights = numpy.maximum ( weights , 1e-4 )
        
    margin = numpy.clip ( preds  , -10.0, 10.0 )
    exp_m  = numpy.exp  ( margin )
    
    # Для Target (1) истинный гессиан равен 0.0 (стабильность обеспечит reg_lambda в параметрах)
    grad = numpy.where ( labels == 1 , -1.0 , exp_m )
    hess = numpy.where ( labels == 1 ,  0.0 , exp_m )
    
    grad *= weights
    hess *= weights
    
    return grad, hess

# =============================================================================
## @class XGBoostIterativeReweighter
#  Iterative weight boosting density reweighter using XGBoost.
class XGBoostIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using XGBoost.
    """
    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        
        config = {
            ## 
            'n_estimators'            :  XGB_ESTIMATORS     ,
            'learning_rate'           :  XGB_LEARNING_RATE  ,
            'max_depth'               :  XGB_MAX_DEPTH      ,
            ## 
            'tree_method'             :  'hist'     , 
            'max_bin'                 :  1023       ,
            ##
            'grow_policy'             : 'depthwise' ,
            'max_leaves'              :  31         ,                     
            ## 
            'min_child_weight'        :  1e-6   ,
            'max_delta_step'          :  0.0    , ## UNSET: Allow unrestricted logits for high density ratios
            ##
            'gamma'                   :  0.0   , 
            'subsample'               :  1.0    , ## UNSET: Retain full sample to resolve sharp density peaks
            'reg_alpha'               :  0.0   ,
            'reg_lambda'              :  0.0   ,                        
            ## 
            'verbosity'               :  0     ,
            'n_jobs'                  :  1     ,
        }

        config.update ( kwargs )
        if 'silent' in config and config.get ( 'silent' ) : config [ 'verbosity' ] = 0        
        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' , XGB_ESTIMATORS ) 
        
        super().__init__ ( original        = original        ,
                           target          = target          ,
                           original_weight = original_weight ,
                           target_weight   = target_weight   , **config )

    # =========================================================================
    @property
    def method ( self ) :
        return method_XGB

    # =========================================================================
    def _fit_iterative_stream ( self               ,
                                X_orig             ,
                                w_orig             ,
                                X_targ             ,
                                w_targ             ) :
        
        import xgboost as XGBoost

        n_orig = len ( X_orig )
        n_targ = len ( X_targ )

        X_comb = numpy.vstack ( [ X_orig, X_targ ] )
        
        abs_w_orig = numpy.abs ( w_orig )
        abs_w_targ = numpy.abs ( w_targ )
        
        norm_targ = w_targ * ( numpy.sum ( abs_w_orig ) / ( numpy.sum ( abs_w_targ ) + 1e-12 ) )

        y_comb = numpy.hstack ( [ numpy.zeros ( n_orig , dtype = numpy.float32 ) ,
                                  numpy.ones  ( n_targ , dtype = numpy.float32 ) ] )

        current_w_orig = abs_w_orig.copy ()

        
        n_estimators  = self.params.get ( 'n_estimators'  , XGB_ESTIMATORS    )
        learning_rate = self.params.get ( 'learning_rate' , XGB_LEARNING_RATE )
        max_depth     = self.params.get ( 'max_depth'     , XGB_MAX_DEPTH     )
    
        models = []

        params = {
            'eval_metric'             : 'logloss'         , 
            'tree_method'             : self.params.get ( 'tree_method'        , 'exact'     ) , 
            'max_bin'                 : self.params.get ( 'max_bin'            , 1023        ) ,
            'min_child_weight'        : self.params.get ( 'min_child_weight'   , 0.001       ) ,
            'max_delta_step'          : self.params.get ( 'max_delta_step'     , 0.0         ) ,
            'subsample'               : self.params.get ( 'subsample'          , 1.0         ) , 
            'gamma'                   : self.params.get ( 'gamma'              , 0.0         ) , 
            'reg_alpha'               : self.params.get ( 'reg_alpha'          , 0.0         ) , 
            'reg_lambda'              : self.params.get ( 'reg_lambda'         , 0.0         ) , 
            'max_depth'               : max_depth ,
            'learning_rate'           : 1.0       , ## Margin is multiplied by external learning_rate
            ## 
            'grow_policy'             : self.params.get ( 'grow_policy' , 'depthwise' ) , 
            'max_leaves'              : self.params.get ( 'max_leaves'  , 0           ) , 
            'disable_default_eval_metric' : 1 ,
            ##             
            'verbosity'               : self.params.get ( 'verbosity' , 0 ) , 
            'n_jobs'                  : self.params.get ( 'n_jobs'    , 1 ) , 
        }
        
        dmat_orig = XGBoost.DMatrix ( X_orig )
        
        prev_model = None
        for it in progress_bar ( n_estimators , description = 'Boosting:' , silent = not self.progress ) :
            
            comb_weights  = numpy.hstack ( [ current_w_orig , numpy.abs ( norm_targ ) ] )
            comb_weights /= ( numpy.sum ( comb_weights ) + 1e-12 )
            
            trn_data = XGBoost.DMatrix ( X_comb , label = y_comb , weight = comb_weights )
            
            # Train next tree using custom_density_ratio_obj instead of default objective
            model    = XGBoost.train   ( params          = params,
                                         dtrain          = trn_data,
                                         num_boost_round = 1,
                                         obj             = xgboost_density_ratio_obj )
            
            prev_model = model

            # Obtain raw margin step exp(lr * margin)
            margin = model.predict ( dmat_orig, output_margin = True )

            # Apply external learning_rate update step
            step_ratio      = numpy.exp   ( numpy.clip ( margin , -5.0 , 5.0 ) )
            current_w_orig *= numpy.power ( step_ratio, learning_rate)
            
            # 2. НОРМИРОВКА (убирает дрейф среднего веса):
            current_w_orig *= (numpy.sum ( abs_w_orig ) / ( numpy.sum ( current_w_orig ) + 1e-12))

            # Log progress
            if ( ( it + 1 ) % 10 == 0 or ( it + 1 ) == n_estimators ) : 
                ratios_curr = current_w_orig / ( abs_w_orig + 1e-12 )
                logger.info ( f" Iteration {it+1}/{n_estimators} -> weights min/max: {ratios_curr.min():.4f} / {ratios_curr.max():.4f}" )
            
            models.append ( model )
            
        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_iter_clip ( final_ratios , w_orig , norm_targ )
        
        return models, final_ratios

    
    # =========================================================================
    def _predict_iterative_stream(self, models, X):
        import xgboost as XGBoost
        dmat = XGBoost.DMatrix(X)
    
        learning_rate = self.params.get('learning_rate', XGB_LEARNING_RATE)
        margin        = numpy.zeros(len(X), dtype=numpy.float32)
    
        for model in models:
            step_margin = model.predict(dmat, output_margin=True)
            # Зажимаем шаг дерева аналогично fit_stream
            step_margin = numpy.clip(step_margin, -5.0 , 5.0)
            margin     += step_margin * learning_rate
        
        weights = numpy.exp(margin)
        # Нормируем предсказанные веса, чтобы их среднее равнялось 1.0 (как при fit)
        weights /= numpy.mean(weights)
    
        return weights

# =============================================================================
## CatBoost 
# =============================================================================

# =============================================================================
## @class CatBoostIterativeReweighter
#  Iterative weight boosting density reweighter using CatBoost.
class CatBoostIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using CatBoost.
    """

    # =========================================================================
    # Custom Exact Density Ratio Objective for CatBoost (Numba JIT Compatible)
    # =========================================================================
    class DensityRatioObj ( object ) :
        
        def calc_ders_range ( self , approxes , targets , weights ) :
            margin = numpy.clip ( approxes , -10.0 , 10.0 )
            exp_m  = numpy.exp  ( margin )
            
            grad = numpy.where ( targets == 1.0 , -1.0  , exp_m )
            hess = numpy.where ( targets == 1.0 ,  1e-5 , exp_m )
            
            if weights is not None and 0 < len ( weights ) :
                # Защита от зануления градиентов (baseline weight)
                weights = numpy.maximum ( weights , 1e-4 )
                grad    = grad * weights
                hess    = hess * weights
                
            der1 = -grad
            der2 =  hess
            
            return list ( zip ( der1 , der2 ) )
        

    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        
        config = {
            ## Core boosting hyperparameters
            'n_estimators'               :  CATB_ESTIMATORS    ,
            'learning_rate'              :  CATB_LEARNING_RATE ,
            'max_depth'                  :  CATB_MAX_DEPTH     ,
            ## 
            'grow_policy'                :  'Lossguide' , # Asymmetric tree growth to capture high ratio peaks
            'max_leaves'                 :  31    ,
            'l2_leaf_reg'                :  0.0   , # Unregularized for exact density ratio convergence
            'subsample'                  :  1.0   , # Full sampling to retain sample statistics
            'random_strength'            :  0.0   ,
            'leaf_estimation_iterations' :  1     ,
            'min_data_in_leaf'           :  1     ,
            ## Execution settings
            'verbose'                    :  False ,
            'thread_count'               :  1     , # Single-threaded weak learner inside Ostap loops
        }
        
        config.update ( kwargs )
        if 'silent' in config and config.get ( 'silent' ) : config [ 'verbose' ] = False        
        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' , CATB_ESTIMATORS ) 
        
        super().__init__ ( original        = original        ,
                           target          = target          ,
                           original_weight = original_weight ,
                           target_weight   = target_weight   , **config )

    # =========================================================================
    @property
    def method ( self ) :
        return method_CATB

    # =========================================================================
    def _fit_iterative_stream ( self               ,
                                X_orig             ,
                                w_orig             ,
                                X_targ             ,
                                w_targ             ) :
        
        import catboost as CatBoost

        n_orig = len ( X_orig )
        n_targ = len ( X_targ )

        X_comb = numpy.vstack ( [ X_orig, X_targ ] )
        
        abs_w_orig = numpy.abs ( w_orig )
        abs_w_targ = numpy.abs ( w_targ )
        
        sum_w_orig = numpy.sum ( abs_w_orig )
        sum_w_targ = numpy.sum ( abs_w_targ )
        
        norm_targ = w_targ * ( sum_w_orig / ( sum_w_targ + 1e-12 ) )

        y_comb = numpy.hstack ( [ numpy.zeros ( n_orig , dtype = numpy.float32 ) ,
                                  numpy.ones  ( n_targ , dtype = numpy.float32 ) ] )

        current_w_orig = abs_w_orig.copy ()
        
        n_estimators  = self.params.get ( 'n_estimators'  , CATB_ESTIMATORS    )
        learning_rate = self.params.get ( 'learning_rate' , CATB_LEARNING_RATE )
        max_depth     = self.params.get ( 'max_depth'     , CATB_MAX_DEPTH     )

        models = []
        pool_orig = CatBoost.Pool ( data = X_orig )

        # Configure CatBoost parameters for single-tree iterative DRE boosting
        params = {
            'loss_function'             : self.DensityRatioObj()                          , # Custom exact DRE objective
            'eval_metric'               : self.params.get ( 'eval_metric' , 'Logloss' )   , # Required by CatBoost for user-defined loss
            'iterations'                : 1                                               , # Train a single tree per iteration
            'learning_rate'             : 1.0                                             , # Step size scaled externally by learning_rate
            'depth'                     : max_depth                                       ,
            'grow_policy'               : self.params.get ( 'grow_policy' , 'Lossguide' ) , # Asymmetric tree growth
            'max_leaves'                : self.params.get ( 'max_leaves'  , 31          ) ,
            'l2_leaf_reg'               : self.params.get ( 'l2_leaf_reg' , 0.0         ) , 
            'subsample'                 : self.params.get ( 'subsample'   , 1.0         ) , 
            'random_strength'           : 0.0                                             , # Disable perturbation noise
            'leaf_estimation_iterations': 1                                               ,
            'min_data_in_leaf'          : self.params.get ( 'min_data_in_leaf' , 1  )     , 
            'verbose'                   : False                                           , 
            'thread_count'              : self.params.get ( 'n_jobs' , 1 )                ,
            'leaf_estimation_method'    : 'Newton'                                        ,
        }

        for it in progress_bar ( n_estimators , description = 'Boosting:' , silent = not self.progress ) :
            
            # Combine current original weights with initial target weights
            comb_weights  = numpy.hstack ( [ current_w_orig , numpy.abs ( norm_targ ) ] )
            comb_weights /= ( numpy.sum ( comb_weights ) + 1e-12 )

            trn_pool = CatBoost.Pool ( data = X_comb , label = y_comb , weight = comb_weights )
            
            # Train weak learner with custom objective
            model = CatBoost.train ( pool    = trn_pool ,
                                     params  = params   ,
                                     verbose = False    )

            # Retrieve raw formula values (margins/logits) directly
            margin = model.predict ( pool_orig , prediction_type = 'RawFormulaVal' )
            if 1 < margin.ndim : margin = margin[:, 0]
            
            # Apply external learning rate scaling to margin update
            step_ratio      = numpy.exp   ( numpy.clip ( margin , -5.0 , 5.0 ) )
            current_w_orig *= numpy.power ( step_ratio , learning_rate )
            
            # Per-iteration weight normalization to prevent constant logit drift
            current_w_orig *= ( sum_w_orig / ( numpy.sum ( current_w_orig ) + 1e-12 ) )
            
            if ( ( it + 1 ) % 10 == 0 or ( it + 1 ) == n_estimators ) : 
                ratios_curr = current_w_orig / ( abs_w_orig + 1e-12 )
                logger.info ( f" Iteration {it+1}/{n_estimators} -> weights min/max: {ratios_curr.min():.4f} / {ratios_curr.max():.4f}" )
            
            models.append ( model )

        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_iter_clip ( final_ratios , w_orig , norm_targ )

        return models, final_ratios

    # =========================================================================
    def _predict_iterative_stream ( self, models, X ) :
        """ Accumulate raw logit predictions across all trained CatBoost weak learners
        and evaluate normalized reweighting factors for unseen dataset X.
        """
        import catboost as CatBoost
        
        learning_rate = self.params.get ( 'learning_rate' , CATB_LEARNING_RATE )
        margin        = numpy.zeros ( len ( X ) , dtype = numpy.float32 )
        pool          = CatBoost.Pool ( data = X )

        for model in models :
            # Extract raw logit prediction
            step_margin = model.predict ( pool , prediction_type = 'RawFormulaVal' )
            if 1 < step_margin.ndim : step_margin = step_margin[:, 0]
            
            # Clip step margin to match fit_stream bounding box
            step_margin = numpy.clip ( step_margin , -5.0 , 5.0 )
            margin     += step_margin * learning_rate
            
        # Exponentiate logits to ratios and normalize mean to 1.0
        weights  = numpy.exp ( margin )
        weights /= ( numpy.mean ( weights ) + 1e-12 )
        
        return weights

"""                
    # =================================================================
            if True : # =======================================================
                # =============================================================
                # Трюк «Феникса» для CatBoost: сохранение в байты/диск и очистка
                import tempfile, os
                with tempfile.NamedTemporaryFile(suffix='.cbm', delete=False) as tmp:
                    tmp_name = tmp.name
                
                current_model.save_model(tmp_name)
                del current_model
                del trn_pool
                
                current_model = CatBoost.CatBoostClassifier()
                current_model.load_model(tmp_name)
                try:
                    os.unlink(tmp_name)
                except OSError:
                    pass
"""

# =============================================================================
## PyTorch  
# =============================================================================

# =============================================================================
## @class PyTorchIterativeReweighter
#  Iterative weight boosting density reweighter using PyTorch (MLP weak learners).
class PyTorchIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using PyTorch.
    """
    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        
        config = {
            ## 
            'n_estimators'            :  TORCH_ESTIMATORS       ,
            'learning_rate'           :  TORCH_LEARNING_RATE    ,
            'hidden_dim'              :  TORCH_HIDDEN_DIM       ,
            'epochs_per_step'         :  TORCH_EPOCHS_PER_STEP  ,
            'batch_size'              :  'auto'                 ,
            'weight_decay'            :  1e-4                   ,
            ## 
            'verbosity'               :  0    ,
            'n_jobs'                  :  1    , 
        }
        
        config.update ( kwargs )
        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' , TORCH_ESTIMATORS ) 

        n_rows     = num_samples ( target ) + num_samples ( original )        
        batch_size = config.get ( 'batch_size' , None ) 
        if   isinstance ( batch_size , string_types  ) and batch_size.lower()  == 'auto' : batch_size = None
        elif isinstance ( batch_size , integer_types ) and not 2 <= batch_size <  n_rows : batch_size = None        
        if batch_size is None : 
            batch_size = 2 ** math.floor ( max ( 2 , math.log2 ( n_rows / 25 ) ) ) 
            config [ 'batch_size' ] = batch_size 

        super().__init__ ( original        = original        ,
                           target          = target          ,
                           original_weight = original_weight ,
                           target_weight   = target_weight   , **config )

    # =========================================================================
    @property
    def method ( self ) :
        return method_TORCH

    # =========================================================================
    def _fit_iterative_stream ( self               ,
                                X_orig             ,
                                w_orig             ,
                                X_targ             ,
                                w_targ             ) :
        
        import torch
        import torch.nn as nn
        import torch.optim as optim

        # Set device (CPU is usually optimal for small reweighting MLPs due to overhead)
        device = torch.device ( 'cpu' )

        n_orig = len ( X_orig )
        n_targ = len ( X_targ )
        n_features = X_orig.shape[ 1 ]

        X_comb = numpy.vstack ( [ X_orig, X_targ ] ).astype ( numpy.float32 )
        
        abs_w_orig = numpy.abs ( w_orig )
        abs_w_targ = numpy.abs ( w_targ )
        
        sum_w_orig = numpy.sum ( abs_w_orig )
        sum_w_targ = numpy.sum ( abs_w_targ )
        
        norm_targ = w_targ * ( sum_w_orig / ( sum_w_targ + 1e-12 ) )

        y_comb = numpy.hstack ( [ numpy.zeros ( n_orig , dtype = numpy.float32 ) ,
                                  numpy.ones  ( n_targ , dtype = numpy.float32  ) ] )

        current_w_orig = abs_w_orig.copy ()
        
        n_estimators    = self.params.get ( 'n_estimators'    , TORCH_ESTIMATORS      )
        learning_rate   = self.params.get ( 'learning_rate'   , TORCH_LEARNING_RATE   )
        hidden_dim      = self.params.get ( 'hidden_dim'      , TORCH_HIDDEN_DIM      )
        epochs_per_step = self.params.get ( 'epochs_per_step' , TORCH_EPOCHS_PER_STEP )
        batch_size      = self.params.get ( 'batch_size'      , 256                   )
        weight_decay    = self.params.get ( 'weight_decay'    , 1e-4                  )

        # Define a simple weak learner MLP architecture
        class WeakMLP ( nn.Module ) :
            def __init__ ( self, input_dim, hidden_dim ) :
                super().__init__ ()
                self.net = nn.Sequential (
                    nn.Linear ( input_dim, hidden_dim ),
                    nn.Tanh   (),
                    nn.Linear ( hidden_dim, hidden_dim ),
                    nn.Tanh   (),
                    nn.Linear ( hidden_dim, 1 ) # No Sigmoid!
                )
            def forward ( self, x ) :
                return self.net ( x )
            
            
        models = []
        X_orig_tensor = torch.tensor ( X_orig.astype ( numpy.float32 ), dtype = torch.float32 ).to ( device )
        X_comb_tensor = torch.tensor ( X_comb, dtype = torch.float32 ).to ( device )
        y_comb_tensor = torch.tensor ( y_comb, dtype = torch.float32 ).unsqueeze ( 1 ).to ( device )

        eps = 1e-4 

        for it in progress_bar ( n_estimators , description = 'Boosting:' , silent = not self.progress ) :
            
            comb_weights  = numpy.hstack ( [ current_w_orig , numpy.abs ( norm_targ ) ] )
            comb_weights /= ( numpy.sum ( comb_weights ) + 1e-12 )
            weights_tensor = torch.tensor ( comb_weights, dtype = torch.float32 ).unsqueeze ( 1 ).to ( device )

            # Instantiate a fresh weak learner for this boosting step
            model = WeakMLP ( n_features, hidden_dim ).to ( device )
            optimizer = optim.Adam ( model.parameters (), lr = 0.01, weight_decay = weight_decay )
            
            ## criterion = nn.BCELoss ( reduction = 'none' )
            criterion = nn.BCEWithLogitsLoss ( reduction = 'none' )

            # Train the weak learner for a few epochs on weighted combined data
            model.train ()
            dataset = torch.utils.data.TensorDataset ( X_comb_tensor, y_comb_tensor, weights_tensor )
            loader  = torch.utils.data.DataLoader ( dataset, batch_size = batch_size, shuffle = True )

            for epoch in range ( epochs_per_step ) :
                for bx, by, bw in loader :
                    optimizer.zero_grad ()
                    preds = model ( bx )
                    loss  = torch.mean ( criterion ( preds, by ) * bw )
                    loss.backward ()
                    optimizer.step ()

            # Evaluate predictions on original sample
            model.eval ()
            with torch.no_grad () :
                p_orig = model ( X_orig_tensor ).cpu ().numpy ().ravel ()

            p_orig = numpy.clip ( p_orig , eps , 1.0 - eps )
            r_step = p_orig / ( 1.0 - p_orig )

            current_w_orig *= numpy.power ( numpy.clip ( r_step , 1e-2 , 1e+2 ) , learning_rate )

            # Store model state dict (lightweight storage instead of full object)
            models.append ( { 'state_dict' : model.state_dict (), 'input_dim' : n_features, 'hidden_dim' : hidden_dim } )
            
        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_iter_clip ( final_ratios , w_orig , norm_targ )

        return models, final_ratios

    
    # =========================================================================
    def _predict_iterative_stream ( self, models, X ) :
        import torch
        import torch.nn as nn

        device = torch.device ( 'cpu' )
        learning_rate = self.params.get ( 'learning_rate' , TORCH_LEARNING_RATE )
        log_ratios    = numpy.zeros ( len ( X ) , dtype = numpy.float32 )
        X_tensor      = torch.tensor ( X.astype ( numpy.float32 ), dtype = torch.float32 ).to ( device )

        # No Sigmoid
        class WeakMLP ( nn.Module ) :
            def __init__ ( self, input_dim, hidden_dim ) :
                super().__init__ ()
                self.net = nn.Sequential (
                    nn.Linear ( input_dim, hidden_dim ),
                    nn.Tanh   (),
                    nn.Linear ( hidden_dim, hidden_dim ),
                    nn.Tanh   (),
                    nn.Linear ( hidden_dim, 1 )
                )
            def forward ( self, x ) :
                return self.net ( x )

        for m_dict in models :
            model = WeakMLP ( m_dict[ 'input_dim' ], m_dict[ 'hidden_dim' ] ).to ( device )
            model.load_state_dict ( m_dict[ 'state_dict' ] )
            model.eval ()

            with torch.no_grad () :
                logits = model ( X_tensor )
                p      = torch.sigmoid ( logits ).cpu ().numpy ().ravel ()

            p           = numpy.clip ( p , 1e-4 , 1.0 - 1e-4 )
            r_step      = p / ( 1.0 - p )
            log_ratios += numpy.log ( numpy.clip ( r_step , 1e-2 , 1e+2 ) ) * learning_rate

        return numpy.exp ( log_ratios )

    
# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    from ostap.stats.tools import ( hasLightGBM ,
                                    hasXGBoost  ,
                                    hasCatBoost ,
                                    hasPyTorch  ) 

    if not hasLightGBM ( False ) : logger.warning  ( "No LightGBM available!" ) 
    if not hasXGBoost  ( False ) : logger.warning  ( "No XGBoost  available!" ) 
    if not hasCatBoost ( False ) : logger.warning  ( "No CatBoost available!" ) 
    if not hasPyTorch  ( False ) : logger.warning  ( "No PyTorch  available!" ) 

# =============================================================================
#                                                                       The END 
# =============================================================================
