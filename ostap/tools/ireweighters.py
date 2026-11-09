#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/ireweighters.py
#  Iterative Density ratio reweighters based on LightGBM, XGBoost, CatBoost & PyTorch 
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
    'LightGBMIterativeReweighter'  , # LightGBM-based iterative density reweighter
    'XGBoostIterativeReweighter'   , # XGBoost-based  iterative density reweighter
    'CatBoostIterativeReweighter'  , # CatBoost-based iterative density reweighter
    'PyTorchIterativeReweighter'   , # PyTorch-based iterative density reweighter
) 
# =============================================================================
from   ostap.core.ostap_types       import num_types, string_types, integer_types
from   ostap.utils.core             import typename
from   ostap.utils.basic            import safe_loop
from   ostap.utils.utils            import classprop 
from   ostap.logger.pretty          import nice_print
from   ostap.math.math_base         import pos_infinity  
from   ostap.utils.progress_bar     import progress_bar, text_bar 
from   ostap.stats.utils            import num_samples, num_features, nEff
from   ostap.tools.base_reweighters import IterativeReweighter  
import ostap.stats.utils            as     U 
import ostap.logger.symbols         as     S
import numpy, math, time   
# =============================================================================
# Logging setup
# =============================================================================
from ostap.logger.logger import getLogger, logAttention
logger = getLogger ( 'ostap.tools.ireweighters' )
# =============================================================================
## Decorations 
# =============================================================================
method_LGBM     =  U.method_name ( 'LightGBM/Iterative' )
method_XGB      =  U.method_name ( 'XGBoost/Iterative'  )
method_CATB     =  U.method_name ( 'CatBoost/Iterative' )
method_TORCH    =  U.method_name ( 'PyTorch/Iterative'  )
##
epoch_symbol    =  U.symbol_EPOCH 
## 
boosting_symbol = S.ladder      if S.show and S.ladder     else 'Boosting'
# =============================================================================
# Global Configuration Constants
# =============================================================================
LGBM_ESTIMATORS       =  300 # High capacity: deep ensemble with low learning rate
LGBM_MAX_DEPTH        =    6 # Deep tree capacity for high dimensions / large statistics
LGBM_LEARNING_RATE    = 0.02 # Small step for stable density ratio convergence

XGB_ESTIMATORS        =  300 # High capacity: deep ensemble with low learning rate
XGB_MAX_DEPTH         =    7 # Deep tree capacity for high dimensions / large statistic
XGB_LEARNING_RATE     = 0.02 # Small step for stable density ratio convergence
XGB_MAX_LEAVES        =   63 # 

CATB_ESTIMATORS       =  300 # High capacity: deep ensemble with low learning rate
CATB_MAX_DEPTH        =    6 # Deep tree capacity for high dimensions / large statistics
CATB_LEARNING_RATE    = 0.10 # Small step for stable density ratio convergence

TORCH_ESTIMATORS      =  300 # Number of boosting steps (epochs/models)
TORCH_EPOCHS_PER_STEP =  100 # Training epochs per boosting step
TORCH_LEARNING_RATE   = 0.02 # Step size for weight boosting
TORCH_HIDDEN_DIM      =  128 # Hidden layer size for small MLP weak learners
TORCH_WEIGHT_DECAY    = 1e-5 
TORCH_BATCH_SIZE      = 'auto'
TORCH_MAX_BATCH_SIZE  = 2 ** 16 
TORCH_HIDDEN_LAYERS   = 32 ,  # ## single hiddle layer with 32 nodes 

TORCH_OPTIMIZER_ADAM  = 'adam'
TORCH_OPTIMIZER_LBFGS = 'lbfgs'
TORCH_OPTIMIZERS      = TORCH_OPTIMIZER_ADAM , TORCH_OPTIMIZER_LBFGS
TORCH_OPTIMIZER       = TORCH_OPTIMIZER_LBFGS 
TORCH_LBFGS_CHUNK     = 200_000

THRESHOLD_NSAMPLES    = 50000
THRESHOLD_NFEATURES   = 10

PRINT_FREQ            = 10 
# =============================================================================
LOG_STEP_CLIPS        = -7.0 , 7.0 
# =============================================================================
## @class LightGBMIterativeReweighter
#  @brief Iterative density ratio reweighter based on LightGBM gradient boosting trees.
#  @details Implements sequential gradient boosting for density ratio estimation 
#           \f$r(x) = \frac{p(x)}{1-p(x)}\f$ using raw margin predictions (log-odds)
#           with dynamically updated sample weights on each boosting iteration.
# =============================================================================
class LightGBMIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using LightGBM backend.
    """

    # =========================================================================
    ## @brief Initialize the LightGBM iterative density reweighter instance.
    #  @param self Instance of the class.
    #  @param original Source feature matrix (N-dimensional array or dataframe).
    #  @param target Target feature matrix (N-dimensional array or dataframe).
    #  @param original_weight Optional observation weights for the source dataset.
    #  @param target_weight Optional observation weights for the target dataset.
    #  @param kwargs Additional keyword arguments for LightGBM configuration.
    # =========================================================================
    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        """ Initialize LightGBM iterative density reweighter with custom defaults. """
        config = {
            'n_estimators'            :  LGBM_ESTIMATORS    ,
            'learning_rate'           :  LGBM_LEARNING_RATE ,
            'max_depth'               :  LGBM_MAX_DEPTH     ,
            'max_bins'                :  1023               ,
            'min_data_in_leaf'        :  10                 ,
            'min_sum_hessian_in_leaf' :  1e-4               , 
            'subsample'               :  1.0                ,
            'bagging_freq'            :  0                  , 
            'lambda_l1'               :  0.0                ,
            'feature_fraction'        :  1.0                , 
            'lambda_l2'               :  0.0                ,                        
            'verbosity'               : -1                  ,
            'print_freq'              : PRINT_FREQ          , 
            'n_jobs'                  :  1                  ,
            'n_splits'                :  5                  ,
            
        }
        
        config.update ( kwargs )
        if 'silent' in config and config.get ( 'silent' ) : config [ 'verbosity' ] = -1        
        if 'num_boost_round'  in config : config [ 'n_estimators'            ] = config.pop ( 'num_boost_round'   , LGBM_ESTIMATORS ) 
        if 'min_child_weight' in config : config [ 'min_sum_hessian_in_leaf' ] = config.pop ( 'min_child_weight'  , 1e-2 )         
        if 'bagging_fraction' in config : config [ 'subsample'               ] = config.pop ( 'bagging_fraction'  , 1.0 )         
        if 'colsample_bytree' in config : config [ 'feature_fraction'        ] = config.pop ( 'colsample_by_tree' , 1.0 )         
        
        n_rows = num_samples ( original ) + num_samples ( target )
        if n_rows < THRESHOLD_NSAMPLES or 1 < config.get ( 'n_splits' , 5 )     : config [ 'subsample'        ] = 1.0    
        if 1 != config.get ( 'subsample' , 1 ) and not 'bagging_freq' in config : config [ 'bagging_freq'     ] = 1            
        if num_features    ( original ) < THRESHOLD_NFEATURES                   : config [ 'feature_fraction' ] = 1.0

        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **config )

    # =========================================================================
    ## @brief Get the short symbol/tag representing the LightGBM reweighter method.
    #  @return String identifier symbol used for logging and console visual display.
    # =========================================================================
    @classprop
    def symbol ( self ) :
        """ Short symbol/tag for reweighter. """
        return U.symbol_LGBM
    
    # =========================================================================
    ## @brief Get the full method identifier name.
    #  @return Method string tag.
    # =========================================================================
    @property
    def method ( self ) :
        """ Return method identifier name. """
        return method_LGBM 
    
    # =========================================================================
    ## @brief Execute the iterative LightGBM boosting stream on a specific cross-validation fold.
    #  @details Sequentially fits LightGBM weak learners (1 tree per step), evaluates raw margins 
    #           (log-odds), applies log-space clipping to [-7.0, 7.0], and updates sample weights 
    #           via exponential scaling.
    #  @param self Instance of the class.
    #  @param X_orig Original feature subset matrix for current fold.
    #  @param w_orig Original sample weights subset for current fold.
    #  @param X_targ Target feature subset matrix for current fold.
    #  @param w_targ Target sample weights subset for current fold.
    #  @return Tuple `(models, final_ratios)` containing trained boosters and final density ratio array.
    # =========================================================================
    def _fit_iterative_stream ( self   ,
                                X_orig ,
                                w_orig ,
                                X_targ ,
                                w_targ ) :
        """ Execute iterative LightGBM boosting stream on fold data. """
        import lightgbm as LightGBM

        n_orig = len ( X_orig )
        n_targ = len ( X_targ )

        X_comb     = numpy.vstack ( [ X_orig, X_targ ] )
        abs_w_orig = numpy.abs ( w_orig )
        abs_w_targ = numpy.abs ( w_targ )
        
        sum_w_orig = numpy.sum ( abs_w_orig )
        sum_w_targ = numpy.sum ( abs_w_targ )
        norm_targ  = w_targ * ( sum_w_orig / ( sum_w_targ + 1e-12 ) )

        y_comb         = numpy.hstack ( [ numpy.zeros ( n_orig, dtype = numpy.float32 ),
                                          numpy.ones  ( n_targ, dtype = numpy.float32 ) ] )
        current_w_orig = abs_w_orig.copy ()
        
        n_estimators  = self.params.get ( 'n_estimators' , LGBM_ESTIMATORS    )
        learning_rate = self.params.get ( 'learning_rate', LGBM_LEARNING_RATE )
        
        max_depth     = self.params.get ( 'max_depth'    , LGBM_MAX_DEPTH     )
        num_leaves    = self.params.get ( 'num_leaves'   , 2 ** max_depth - 1 )
    
        models = []

        params = {
            'objective'               : self.params.get ( 'objective'               , 'binary'         ) ,     
            'metric'                  : self.params.get ( 'metric'                  , 'binary_logloss' ) , 
            'max_bins'                : self.params.get ( 'max_bins'                , 1023             ) , 
            'min_data_in_leaf'        : self.params.get ( 'min_data_in_leaf'        , 30               ) , 
            'min_sum_hessian_in_leaf' : self.params.get ( 'min_sum_hessian_in_leaf' , 1e-3             ) , 
            'feature_fraction'        : self.params.get ( 'feature_fraction'        , 1.0              ) , 
            'subsample'               : self.params.get ( 'subsample'               , 1.0              ) , 
            'bagging_freq'            : self.params.get ( 'bagging_freq'            , 0                ) , 
            'lambda_l1'               : self.params.get ( 'lambda_l1'               , 0.0              ) , 
            'lambda_l2'               : self.params.get ( 'lambda_l2'               , 0.001            ) , 
            'max_depth'               : max_depth                                                      ,
            'num_leaves'              : num_leaves                                                     ,
            'learning_rate'           : 1.0                                                            , 
            'verbosity'               : self.params.get ( 'verbosity'               , -1               ) , 
            'n_jobs'                  : self.params.get ( 'n_jobs'                  ,  1               ) , 
        }
        
        print_freq = self.params.get ( 'print_freq' , PRINT_FREQ )

        neff0 = nEff ( X_orig, current_w_orig )

        description = boosting_symbol + ' ' + S.fixed_width ( self.symbol , 3 ) + ':'        
        if   self.progress : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = not self.progress ) 
        elif self.verbose  : the_bar = text_bar     ( n_estimators , full = True , description = description , silent =     self.silent   )
        else               : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = True  ) 

        t_start = time.time () 

        for it , _ ,  bar in safe_loop ( the_bar ) :
            
            comb_weights  = numpy.hstack ( [ current_w_orig, numpy.abs ( norm_targ ) ] )
            comb_weights /= ( numpy.sum ( comb_weights ) + 1e-12 )

            trn_data = LightGBM.Dataset ( X_comb, label = y_comb, weight = comb_weights, free_raw_data = True )            
            model    = LightGBM.train   ( params          = params   ,
                                          train_set       = trn_data ,
                                          num_boost_round = 1        )

            raw_vals = model.predict ( X_orig, raw_score = True )
            if 1 < raw_vals.ndim : raw_vals = raw_vals [ :, 0 ]
            
            log_r_step      = numpy.clip ( raw_vals   , *LOG_STEP_CLIPS )
            current_w_orig *= numpy.exp  ( log_r_step * learning_rate   )

            if bar and not self.silent and ( it < 5  or ( it + 1 ) % print_freq == 0 or ( it + 5 ) == n_estimators ) :
                ratios_curr = current_w_orig / ( abs_w_orig + 1e-12 )
                neff        = nEff ( X_orig, current_w_orig )
                wmin        = nice_print ( float ( ratios_curr.min () ) , precision = 2, with_sign = False, strip = True )
                wmax        = nice_print ( float ( ratios_curr.max () ) , precision = 2, with_sign = False, strip = True )
                ne0         = nice_print ( float ( neff0              ) , precision = 2, with_sign = False, strip = True )
                ne          = nice_print ( float ( neff               ) , precision = 2, with_sign = False, strip = True )

                dt       = ( time.time () - t_start ) / ( it + 1 )
                if dt < 1 : dt = '%dms'  % ( math.ceil ( 1000 * dt ) )
                else      : dt = '%.1fs' %                      dt                 
                message  = f"Step {it+1:4d}/{n_estimators:<4d} | min/max: {wmin}/{wmax} | #eff: {ne0} {S.arrow_right} {ne} | Time/step:{dt}" 
                bar     += message

            model_bytes = model.model_to_string ()
            model.free_dataset ()
            del model, trn_data
            
            model = LightGBM.Booster ( model_str = model_bytes )
            del model_bytes
            
            models.append ( model )
            
        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_clip ( final_ratios, w_orig )

        return models, final_ratios

    # =========================================================================
    ## @brief Accumulate raw margin predictions and compute log-ratios across sequential LightGBM models.
    #  @details Evaluates each trained LightGBM booster on feature matrix X in raw log-odds space,
    #           clips log-ratio steps to [-7.0, 7.0], accumulates them scaled by learning rate,
    #           and returns final exponential density ratio predictions.
    #  @param self Instance of the class.
    #  @param models List of trained LightGBM booster models comprising sequential chain.
    #  @param X Feature matrix of unseen samples to evaluate.
    #  @return Array of predicted density ratios r(x).
    # =========================================================================
    def _predict_iterative_stream ( self, models, X ) :
        """ Accumulate predictions across sequential LightGBM models. """
        
        learning_rate = self.params.get ( 'learning_rate', LGBM_LEARNING_RATE )
        log_ratios    = numpy.zeros ( len ( X ), dtype = numpy.float32 )

        for model in models :
            raw_vals    = model.predict ( X, raw_score = True )
            if 1 < raw_vals.ndim : raw_vals = raw_vals [ :, 0 ]
            
            log_r_step  = numpy.clip ( raw_vals, *LOG_STEP_CLIPS )
            log_ratios += log_r_step * learning_rate

        return numpy.exp ( log_ratios )

# =============================================================================
## @class XGBoostIterativeReweighter
#  @brief Iterative weight boosting density reweighter using XGBoost backend.
#  @details Implements sequential gradient boosting for density ratio estimation 
#           \f$r(x) = \frac{p(x)}{1-p(x)}\f$ using raw margin predictions (output_margin)
#           with dynamically updated sample weights on each boosting iteration.
# =============================================================================
class XGBoostIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using XGBoost backend.
    """

    # =========================================================================
    ## @brief Initialize the XGBoost iterative density reweighter instance.
    #  @param self Instance of the class.
    #  @param original Source feature matrix (N-dimensional array or dataframe).
    #  @param target Target feature matrix (N-dimensional array or dataframe).
    #  @param original_weight Optional observation weights for the source dataset.
    #  @param target_weight Optional observation weights for the target dataset.
    #  @param kwargs Additional keyword arguments for XGBoost configuration.
    # =========================================================================
    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        """ Initialize XGBoost iterative density reweighter with custom defaults. """
        config = {
            'n_estimators'      :  XGB_ESTIMATORS    ,
            'learning_rate'     :  XGB_LEARNING_RATE ,
            'max_depth'         :  XGB_MAX_DEPTH     ,
            'tree_method'       :  'hist'            , 
            'max_bin'           :  1023              , 
            'grow_policy'       :  'depthwise'       , 
            'max_leaves'        :  XGB_MAX_LEAVES    , 
            'min_child_weight'  :  1e-3              , 
            'max_delta_step'    :  9.0               , 
            'gamma'             :  0.0               , 
            'subsample'         :  1.0               , 
            'colsamples_bytree' :  1.0               , 
            'reg_alpha'         :  0.0               , 
            'reg_lambda'        :  0.0               , 
            'verbosity'         :  0                 ,
            'print_freq'        : PRINT_FREQ         , 
            'n_jobs'            :  1                 ,
            'n_splits'          :  5                 ,
            
        }

        config.update ( kwargs )
        if 'silent' in config and config.get ( 'silent' ) : config [ 'verbosity' ] = 0        
        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round', XGB_ESTIMATORS ) 
        
        n_rows = num_samples ( original ) + num_samples ( target )
        if n_rows < THRESHOLD_NSAMPLES or 1 < config.get ( 'n_splits' , 5 ) : config [ 'subsample'         ] = 1.0
        if num_features      ( original ) < THRESHOLD_NFEATURES             : config [ 'colsamples_bytree' ] = 1.0

        if n_rows < 10000 : config [ 'tree_method' ]  = 'exact'

        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **config )

    # =========================================================================
    ## @brief Get the short symbol/tag representing the XGBoost reweighter method.
    #  @return String identifier symbol used for logging and visual presentation.
    # =========================================================================
    @classprop
    def symbol ( self ) :
        """ Short symbol/tag for reweighter. """
        return U.symbol_XGB

    # =========================================================================
    ## @brief Get the full method identifier name.
    #  @return Method string tag.
    # =========================================================================
    @property
    def method ( self ) :
        """ Return method identifier name. """
        return method_XGB

    # =========================================================================
    ## @brief Execute the iterative XGBoost boosting stream on a specific cross-validation fold.
    #  @details Sequentially fits XGBoost weak learners (1 tree per step), evaluates raw margins 
    #           (output_margin), applies log-space clipping to [-7.0, 7.0], and updates sample weights 
    #           via exponential scaling.
    #  @param self Instance of the class.
    #  @param X_orig Original feature subset matrix for current fold.
    #  @param w_orig Original sample weights subset for current fold.
    #  @param X_targ Target feature subset matrix for current fold.
    #  @param w_targ Target sample weights subset for current fold.
    #  @return Tuple `(models, final_ratios)` containing trained boosters and final density ratio array.
    # =========================================================================
    def _fit_iterative_stream ( self   ,
                                X_orig ,
                                w_orig ,
                                X_targ ,
                                w_targ ) :
        """ Execute iterative XGBoost boosting stream on fold data. """
        import xgboost as XGBoost

        n_orig     = len ( X_orig )
        n_targ     = len ( X_targ )
        X_comb     = numpy.vstack ( [ X_orig, X_targ ] )
        abs_w_orig = numpy.abs ( w_orig )
        abs_w_targ = numpy.abs ( w_targ )
        
        sum_w_orig = numpy.sum ( abs_w_orig )
        sum_w_targ = numpy.sum ( abs_w_targ )
        norm_targ  = w_targ * ( sum_w_orig / ( sum_w_targ + 1e-12 ) )

        y_comb         = numpy.hstack ( [ numpy.zeros ( n_orig, dtype = numpy.float32 ),
                                          numpy.ones  ( n_targ, dtype = numpy.float32 ) ] )
        current_w_orig = abs_w_orig.copy ()

        n_estimators  = self.params.get ( 'n_estimators' , XGB_ESTIMATORS    )
        learning_rate = self.params.get ( 'learning_rate', XGB_LEARNING_RATE )
        
        models = []

        params = {
            'objective'        : 'binary:logistic' , 
            'eval_metric'      : 'logloss'         ,
            'learning_rate'    : 1.0               ,
            ##
            'tree_method'      : self.params.get ( 'tree_method'      , 'hist'      ) , 
            'max_bin'          : self.params.get ( 'max_bin'          , 1023        ) ,
            'min_child_weight' : self.params.get ( 'min_child_weight' , 1e-3        ) ,
            'max_delta_step'   : self.params.get ( 'max_delta_step'   , 9.0         ) , 
            'subsample'        : self.params.get ( 'subsample'        , 1.0         ) ,
            'colsample_bytree' : self.params.get ( 'colsample_bytree' , 1.0         ) ,
            'gamma'            : self.params.get ( 'gamma'            , 0.0         ) ,             
            'reg_alpha'        : self.params.get ( 'reg_alpha'        , 0.0         ) , 
            'reg_lambda'       : self.params.get ( 'reg_lambda'       , 0.0         ) ,
            ##
            'max_depth'        : self.params.get ( 'max_depth'        , XGB_MAX_DEPTH   ) , 
            'grow_policy'      : self.params.get ( 'grow_policy'      , 'depthwise'     ) , 
            'max_leaves'       : self.params.get ( 'max_leaves'       , XGB_MAX_LEAVES  ) ,
            ## 
            'n_jobs'           : self.params.get ( 'n_jobs'           , 1               ) , 
            'verbosity'        : 0                 , 
        }
        print_freq        = self.params.get ( 'print_freq' , PRINT_FREQ )
        
        neff0 = nEff ( X_orig, current_w_orig )
        
        description = boosting_symbol + ' ' + S.fixed_width ( self.symbol , 3 ) + ':'        
        if   self.progress : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = not self.progress ) 
        elif self.verbose  : the_bar = text_bar     ( n_estimators , full = True , description = description , silent =     self.silent   )
        else               : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = True  ) 

        t_start = time.time () 
        for it , _ ,  bar in safe_loop ( the_bar ) :
                        
            comb_weights  = numpy.hstack ( [ current_w_orig, numpy.abs ( norm_targ ) ] )
            comb_weights /= ( numpy.sum ( comb_weights ) + 1e-12 )
            
            dmat_comb = XGBoost.DMatrix ( X_comb, label = y_comb, weight = comb_weights )
            
            model = XGBoost.train ( params          = params,
                                    dtrain          = dmat_comb,
                                    num_boost_round = 1 )
            
            dmat_orig = XGBoost.DMatrix ( X_orig )
            # Retrieve raw margins (output_margin=True) for exact mathematical consistency
            raw_vals  = model.predict ( dmat_orig, output_margin = True )
            if 1 < raw_vals.ndim : raw_vals = raw_vals [ :, 0 ]
            
            log_r_step      = numpy.clip ( raw_vals, *LOG_STEP_CLIPS )
            current_w_orig *= numpy.exp ( log_r_step * learning_rate )
        
            if bar and not self.silent and ( it < 5  or ( it + 1 ) % print_freq == 0 or ( it + 5 ) == n_estimators ) :
                ratios_curr = current_w_orig / ( abs_w_orig + 1e-12 )
                neff        = nEff ( X_orig, current_w_orig )
                wmin        = nice_print ( float ( ratios_curr.min () ) , precision = 2, with_sign = False, strip = True )
                wmax        = nice_print ( float ( ratios_curr.max () ) , precision = 2, with_sign = False, strip = True )
                ne0         = nice_print ( float ( neff0              ) , precision = 2, with_sign = False, strip = True )
                ne          = nice_print ( float ( neff               ) , precision = 2, with_sign = False, strip = True )

                dt       = ( time.time () - t_start ) / ( it + 1 )
                if dt < 1 : dt = '%dms'  % ( math.ceil ( 1000 * dt ) )
                else      : dt = '%.1fs' %                      dt                 
                message  = f"Step {it+1:4d}/{n_estimators:<4d} | min/max: {wmin}/{wmax} | #eff: {ne0} {S.arrow_right} {ne} | Time/step:{dt}" 
                bar     += message
                
            models.append ( model )
            
        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_clip ( final_ratios, w_orig )

        return models, final_ratios

    # =========================================================================
    ## @brief Accumulate raw margin predictions across sequential XGBoost models.
    #  @details Evaluates each trained XGBoost booster on feature matrix X in raw log-odds space,
    #           clips log-ratio steps to [-7.0, 7.0], accumulates them scaled by learning rate,
    #           and returns final exponential density ratio predictions.
    #  @param self Instance of the class.
    #  @param models List of trained XGBoost booster models comprising sequential chain.
    #  @param X Feature matrix of unseen samples to evaluate.
    #  @return Array of predicted density ratios r(x).
    # =========================================================================
    def _predict_iterative_stream ( self, models, X ) :
        """ Accumulate predictions across sequential XGBoost models. """
        import xgboost as XGBoost
        dmat          = XGBoost.DMatrix ( X )
        learning_rate = self.params.get ( 'learning_rate', XGB_LEARNING_RATE )
        log_ratios    = numpy.zeros ( len ( X ), dtype = numpy.float32 )

        for model in models :
            raw_vals    = model.predict ( dmat, output_margin = True )
            if 1 < raw_vals.ndim : raw_vals = raw_vals [ :, 0 ]
            
            log_r_step  = numpy.clip ( raw_vals, *LOG_STEP_CLIPS )
            log_ratios += log_r_step * learning_rate

        return numpy.exp ( log_ratios )


# =============================================================================
## @class CatBoostIterativeReweighter
#  @brief Iterative density ratio reweighter based on CatBoost gradient boosting trees.
#  @details Implements sequential gradient boosting for density ratio estimation 
#           \f$r(x) = \frac{p(x)}{1-p(x)}\f$ using Symmetric (Oblivious) Trees and 
#           raw margin outputs (RawFormulaVal) with dynamically updated sample weights.
# =============================================================================
class CatBoostIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using CatBoost with Oblivious Trees.
    """

    # =========================================================================
    ## @brief Initialize the CatBoost iterative density reweighter instance.
    #  @param self Instance of the class.
    #  @param original Source feature matrix (N-dimensional array or dataframe).
    #  @param target Target feature matrix (N-dimensional array or dataframe).
    #  @param original_weight Optional observation weights for the source dataset.
    #  @param target_weight Optional observation weights for the target dataset.
    #  @param kwargs Additional keyword arguments for CatBoost configuration.
    # =========================================================================
    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        """ Initialize CatBoost iterative density reweighter with custom defaults. """
        config = {
            'n_estimators'               :  CATB_ESTIMATORS    ,
            'learning_rate'              :  CATB_LEARNING_RATE ,
            'max_depth'                  :  CATB_MAX_DEPTH     ,
            'grow_policy'                :  'SymmetricTree'    ,
            'l2_leaf_reg'                :  0.001              ,
            'subsample'                  :  1.0                ,
            'rsm'                        :  1.0                , 
            'random_strength'            :  0.0                ,
            'leaf_estimation_iterations' :  1                  ,
            'min_data_in_leaf'           :  1                  ,
            'verbose'                    :  False              ,
            'print_freq'                 : PRINT_FREQ          , 
            'thread_count'               :  1                  ,
            'n_splits'                   :  5                  , 
        }
        
        config.update ( kwargs )
        if 'silent'           in config and config.get ( 'silent' ) : config [ 'verbose' ] = 0        
        if 'num_boost_round'  in config : config [ 'n_estimators'   ] = config.pop ( 'num_boost_round'   , CATB_ESTIMATORS ) 
        if 'colsample_bytree' in config : config [ 'rsm'            ] = config.pop ( 'colsample_by_tree' , 1.0 )         
        if 'depth'            in config : config [ 'max_depth'      ] = config.pop ( 'depth'             , CATB_MAX_DEPTH  )            

        n_rows = num_samples ( original ) + num_samples ( target )
        if n_rows < THRESHOLD_NSAMPLES or 1 < config.get ( 'n_splits' , 5 ) : config [ 'subsample' ] = 1.0
        if num_features    ( original ) < THRESHOLD_NFEATURES               : config [ 'rsm'       ] = 1.0

        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **config )

    # =========================================================================
    ## @brief Get the short symbol/tag representing the CatBoost reweighter method.
    #  @return String identifier symbol used for logging and visual presentation.
    # =========================================================================
    @classprop
    def symbol ( self ) :
        """ Short symbol/tag for reweighter. """
        return U.symbol_CATB

    # =========================================================================
    ## @brief Get the full method identifier name.
    #  @return Method string tag.
    # =========================================================================
    @property
    def method ( self ) :
        """ Return method identifier name. """
        return method_CATB

    # =========================================================================
    ## @brief Execute the iterative CatBoost boosting stream on a specific cross-validation fold.
    #  @details Sequentially fits CatBoost weak learners (1 tree per step), evaluates raw margins 
    #           (RawFormulaVal), applies log-space clipping to [-7.0, 7.0], and updates sample weights 
    #           via exponential scaling.
    #  @param self Instance of the class.
    #  @param X_orig Original feature subset matrix for current fold.
    #  @param w_orig Original sample weights subset for current fold.
    #  @param X_targ Target feature subset matrix for current fold.
    #  @param w_targ Target sample weights subset for current fold.
    #  @return Tuple `(models, final_ratios)` containing trained boosters and final density ratio array.
    # =========================================================================
    def _fit_iterative_stream ( self   ,
                                X_orig ,
                                w_orig ,
                                X_targ ,
                                w_targ ) :
        """ Execute iterative CatBoost boosting stream on fold data. """
        import catboost as CatBoost

        n_orig     = len ( X_orig )
        n_targ     = len ( X_targ )
        X_comb     = numpy.vstack ( [ X_orig, X_targ ] )
        abs_w_orig = numpy.abs ( w_orig )
        abs_w_targ = numpy.abs ( w_targ )
        
        sum_w_orig = numpy.sum ( abs_w_orig )
        sum_w_targ = numpy.sum ( abs_w_targ )
        norm_targ  = w_targ * ( sum_w_orig / ( sum_w_targ + 1e-12 ) )

        y_comb         = numpy.hstack ( [ numpy.zeros ( n_orig, dtype = numpy.float32 ),
                                          numpy.ones  ( n_targ, dtype = numpy.float32 ) ] )
        current_w_orig = abs_w_orig.copy ()
        
        n_estimators  = self.params.get ( 'n_estimators' , CATB_ESTIMATORS    )
        learning_rate = self.params.get ( 'learning_rate', CATB_LEARNING_RATE )

        models    = []
        pool_orig = CatBoost.Pool ( data = X_orig )

        params = {
            'loss_function'             : 'Logloss'                                       ,
            'eval_metric'               : 'Logloss'                                       ,
            'iterations'                : 1                                               ,
            'learning_rate'             : 1.0                                             ,
            'depth'                     : self.params.get ( 'max_depth'                   , CATB_MAX_DEPTH  ) , 
            'grow_policy'               : self.params.get ( 'grow_policy'                 , 'SymmetricTree' ) , 
            'l2_leaf_reg'               : self.params.get ( 'l2_leaf_reg'                 , 0.001           ) , 
            'subsample'                 : self.params.get ( 'subsample'                   , 1.0             ) ,
            'rsm'                       : self.params.get ( 'rsm'                         , 1.0             ) ,            
            'random_strength'           : self.params.get ( 'random_strength'             , 0.0             ) , 
            'leaf_estimation_iterations': self.params.get ( 'leaf_estimation_iterations'  , 1               ) , 
            'min_data_in_leaf'          : self.params.get ( 'min_data_in_leaf'            , 1               ) , 
            'verbose'                   : False                                           , 
            'thread_count'              : self.params.get ( 'n_jobs', 1 )                 ,
            'leaf_estimation_method'    : 'Newton'                                        ,
        }
        print_freq        = self.params.get ( 'print_freq' , PRINT_FREQ )
        
        neff0 = nEff ( X_orig, current_w_orig )

        description = boosting_symbol + ' ' + S.fixed_width ( self.symbol , 3 ) + ':'        
        if   self.progress : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = not self.progress ) 
        elif self.verbose  : the_bar = text_bar     ( n_estimators , full = True , description = description , silent =     self.silent   )
        else               : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = True  ) 

        t_start = time.time () 

        for it , _ ,  bar in safe_loop ( the_bar ) :
            
            comb_weights  = numpy.hstack ( [ current_w_orig, numpy.abs ( norm_targ ) ] )
            comb_weights /= ( numpy.sum ( comb_weights ) + 1e-12 )

            trn_pool = CatBoost.Pool ( data = X_comb, label = y_comb, weight = comb_weights )
            
            model    = CatBoost.train ( pool    = trn_pool ,
                                        params  = params   ,
                                        verbose = False    )

            # Retrieve raw margins (RawFormulaVal) for exact mathematical consistency
            raw_vals = model.predict ( pool_orig, prediction_type = 'RawFormulaVal' )
            if 1 < raw_vals.ndim : raw_vals = raw_vals [ :, 0 ]
            
            log_r_step      = numpy.clip ( raw_vals   , *LOG_STEP_CLIPS )
            current_w_orig *= numpy.exp  ( log_r_step * learning_rate   )

            if bar and not self.silent and ( it < 5  or ( it + 1 ) % print_freq == 0 or ( it + 5 ) == n_estimators ) :
                ratios_curr = current_w_orig / ( abs_w_orig + 1e-12 )
                neff        = nEff ( X_orig, current_w_orig )
                wmin        = nice_print ( float ( ratios_curr.min () ) , precision = 2, with_sign = False, strip = True )
                wmax        = nice_print ( float ( ratios_curr.max () ) , precision = 2, with_sign = False, strip = True ) 
                ne0         = nice_print ( float ( neff0              ) , precision = 2, with_sign = False, strip = True )
                ne          = nice_print ( float ( neff               ) , precision = 2, with_sign = False, strip = True )

                dt       = ( time.time () - t_start ) / ( it + 1 )
                if dt < 1 : dt = '%dms'  % ( math.ceil ( 1000 * dt ) )
                else      : dt = '%.1fs' %                      dt                 
                message  = f"Step {it+1:4d}/{n_estimators:<4d} | min/max: {wmin}/{wmax} | #eff: {ne0} {S.arrow_right} {ne} | Time/step:{dt}" 
                bar     += message

            models.append ( model )

        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_clip ( final_ratios, w_orig )

        return models, final_ratios

    # =========================================================================
    ## @brief Accumulate raw margin predictions across sequential CatBoost models.
    #  @details Evaluates each trained CatBoost booster on feature matrix X in raw log-odds space
    #           (RawFormulaVal), clips log-ratio steps to [-7.0, 7.0], accumulates them scaled 
    #           by learning rate, and returns final exponential density ratio predictions.
    #  @param self Instance of the class.
    #  @param models List of trained CatBoost booster models comprising sequential chain.
    #  @param X Feature matrix of unseen samples to evaluate.
    #  @return Array of predicted density ratios r(x).
    # =========================================================================
    def _predict_iterative_stream ( self, models, X ) :
        """ Accumulate predictions across sequential CatBoost models. """
        import catboost as CatBoost
        
        learning_rate = self.params.get ( 'learning_rate', CATB_LEARNING_RATE )
        pool          = CatBoost.Pool ( data = X )
        log_ratios    = numpy.zeros ( len ( X ), dtype = numpy.float32 )

        for model in models :
            raw_vals    = model.predict ( pool, prediction_type = 'RawFormulaVal' )
            if 1 < raw_vals.ndim : raw_vals = raw_vals [ :, 0 ]
            
            log_r_step  = numpy.clip ( raw_vals , *LOG_STEP_CLIPS )
            log_ratios += log_r_step * learning_rate

        return numpy.exp ( log_ratios )


# =============================================================================
## maximal batch size for PyTorch 
# =============================================================================
## @brief Calculate optimal batch size for PyTorch iterative boosting.
#  @details Dynamically balances Full-Batch stability for smaller datasets with 
#           SIMD efficiency, RAM limits, and BatchNorm1d stability for larger datasets.
#  @param n_samples  Total number of samples (N_orig + N_targ)
#  @param n_features Number of input features (D)
#  @param max_val    Maximum allowed batch size cap (default: 2^16 = 65536)
#  @return Optimal integer batch size (power of 2)
# =============================================================================
def optimal_batch_size ( n_samples      ,
                         n_features     ,
                         max_batch_size = TORCH_MAX_BATCH_SIZE ) :
    """ Calculate optimal batch size dynamically.
    """
    if n_samples <= max_batch_size  :return n_samples

    # Ensure at least 64 samples per feature axis to stabilize BatchNorm statistics
    target_bs = max ( n_samples // 2 , 64 * n_features )
    
    # Smoothly round down to the nearest power of 2
    batch_size = 2 ** math.floor ( math.log2 ( target_bs ) )
    
    return min ( batch_size , max_batch_size )


# =============================================================================
## @class PyTorchIterativeReweighter
#  @brief Iterative weight boosting density reweighter using PyTorch (MLP weak learners).
#  @details Implements sequential gradient boosting for density ratio estimation 
#           \f$r(x) = \frac{p(x)}{1-p(x)}\f$ using feature-dependent Multi-Layer 
#           Perceptron (MLP) weak learners. Supports flexible optimization backends 
#           (Adam and quasi-Newton L-BFGS) with automatic static chunking for large datasets.
class PyTorchIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using PyTorch backend.
    
    This class sequentially trains PyTorch MLP weak learners to estimate 
    density ratios between target and original feature distributions. It 
    supports dynamic batch size optimization, spatial variance control via 
    zero-mean centering, and early stopping based on sliding loss improvement.
    """

    # =========================================================================
    ## @brief Initialize the PyTorch iterative density reweighter instance.
    #  @details Sets up hyperparameter configuration defaults, validates batch size, 
    #           configures the computation device (CPU or GPU), and initializes 
    #           the underlying IterativeReweighter framework.
    #  @param self Instance of the class.
    #  @param original Source feature matrix (N-dimensional array or dataframe).
    #  @param target Target feature matrix (N-dimensional array or dataframe).
    #  @param original_weight Optional observation weights for the source dataset.
    #  @param target_weight Optional observation weights for the target dataset.
    #  @param kwargs Additional keyword arguments for PyTorch configuration:
    #         - n_estimators (int): Number of boosting steps (default: 300).
    #         - hidden_layers (tuple): Architecture of hidden layers (default: (128, 64)).
    #         - learning_rate (float): Learning rate for density ratio boosting step (default: 0.02).
    #         - device (str or torch.device): Computation device ('cpu', 'cuda', or None for auto-detection).
    # =========================================================================
    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        """ Initialize PyTorch iterative density reweighter with customized hyperparameter defaults.

        Merges user arguments with default configurations, resolves batch size constraints,
        configures execution device, and initializes parent cross-validation infrastructure.
        """
        import torch

        config = {
            'n_estimators'    :  TORCH_ESTIMATORS       ,
            'hidden_layers'   :  TORCH_HIDDEN_LAYERS    ,
            'learning_rate'   :  TORCH_LEARNING_RATE    ,
            'weight_decay'    :  TORCH_WEIGHT_DECAY     , 
            'batch_size'      :  TORCH_BATCH_SIZE       ,
            'lbfgs_chunk'     :  TORCH_LBFGS_CHUNK      ,
            'device'          :  None                   , # Auto-detect GPU/CPU if not specified
            ##
            'epochs_per_step' :  TORCH_EPOCHS_PER_STEP  ,
            'min_epoches'     :  5                      ,
            'min_delta'       :  1e-4                   ,
            'patience'        :  5                      ,
            ## 
            'optimizer'       :  TORCH_OPTIMIZER        , # Default optimizer type ('adam' or 'lbfgs')
            # 
            'verbosity'       :  0                      ,
            'print_freq'      : PRINT_FREQ              , 
            'n_jobs'          :  1                      , 
        }
        
        config.update ( kwargs )
        # Alias resolution for backward compatibility
        if 'num_boost_round' in config : config [ 'n_estimators'  ] = config.pop ( 'num_boost_round', TORCH_ESTIMATORS ) 
        if 'hidden_dim'      in config : config [ 'hidden_layers' ] = config.pop ( 'hidden_dim'     , ( 128, 64 )      )

        # Automatic GPU / CPU device selection in constructor
        user_device = config.get ( 'device', None )
        if user_device is None                        : config [ 'device' ] = torch.device ( 'cuda' if torch.cuda.is_available () else 'cpu' )
        elif isinstance ( user_device, string_types ) : config [ 'device' ] = torch.device ( user_device )

        n_rows     = num_samples  ( target ) + num_samples ( original )
        n_features = num_features ( target )

        batch_size = config.get ( 'batch_size' , 'auto' ) 
        
        # Validate user-provided integer batch size or dynamically compute optimal size
        if isinstance ( batch_size, integer_types ) and 2 <= batch_size <= min ( n_rows , TORCH_MAX_BATCH_SIZE ) :
            if batch_size & ( batch_size - 1 ) : batch_size = 2 ** math.floor ( math.log2 ( batch_size ) )
        else : 
            batch_size = optimal_batch_size ( n_rows, n_features, TORCH_MAX_BATCH_SIZE )
            
        config [ 'batch_size' ] = batch_size
        
        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **config )

    # =========================================================================
    ## @brief Get the short symbol/tag representing the PyTorch reweighter method.
    #  @return String identifier symbol used for logging and console visual display.
    # =========================================================================
    @classprop
    def symbol ( self ) :
        """ Return short symbol/tag for reweighter console visual representation. """
        return U.symbol_TORCH

    # =========================================================================
    ## @brief Get the full method identifier name.
    #  @return Method string tag.
    # =========================================================================
    @property
    def method ( self ) :
        """ Return method identifier name. """
        return method_TORCH

    # =========================================================================
    ## @brief Factory method generating a dynamic PyTorch Module class for weak learners.
    #  @details Constructs a dynamic `WeakMLP` subclass supporting arbitrary hidden layer 
    #           architectures. Includes BatchNorm1d for input normalization and forces 
    #           `bias=False` on the output linear layer to eliminate global offset learning.
    #  @return WeakMLP class inheriting from `torch.nn.Module`.
    # =========================================================================
    def get_MLP_class ( self ) :
        """ Dynamically import torch and construct the WeakMLP PyTorch Module class.
        
        Returns:
            class: A PyTorch nn.Module class tailored for spatial density ratio step learning.
        """
        import torch

        class WeakMLP ( torch.nn.Module ) :
            """ Weak Multi-Layer Perceptron neural network module. """
            def __init__ ( self, input_dim, hidden_layers = ( 128, ) ) :
                super ().__init__ ()
                
                if isinstance ( hidden_layers, int ) :
                    hidden_layers = ( hidden_layers, )

                # Input feature batch normalization to prevent gradient scaling anomalies
                self.bn = torch.nn.BatchNorm1d ( input_dim )
                
                layers = []
                curr_dim = input_dim
                for h_dim in hidden_layers :
                    linear = torch.nn.Linear ( curr_dim, h_dim )
                    torch.nn.init.kaiming_uniform_ ( linear.weight, a = 0.1 )
                    torch.nn.init.zeros_ ( linear.bias )
                    
                    layers.append ( linear )
                    layers.append ( torch.nn.LeakyReLU ( 0.1 ) )
                    curr_dim = h_dim
                    
                self.hidden_net = torch.nn.Sequential ( *layers )
                
                # CRITICAL: Setting bias=False forces the model to rely purely on spatial 
                # feature correlations rather than fitting a constant background offset.
                self.fc_out = torch.nn.Linear ( curr_dim, 1, bias = False ) 
                torch.nn.init.xavier_uniform_ ( self.fc_out.weight )

            def forward ( self, x ) :
                """ Execute forward pass through BatchNorm, hidden layers, and linear head. """
                h = self.bn ( x )
                h = self.hidden_net ( h )
                return self.fc_out ( h )
                
        return WeakMLP

    # =========================================================================
    ## @brief Execute the iterative PyTorch boosting stream on a specific cross-validation fold.
    #  @details Sequentially fits PyTorch MLP weak learners, computes logit predictions,
    #           applies zero-mean centering, clips log-ratio steps to [-7.0, 7.0], 
    #           and updates sample weights via exponential scaling. Supports static chunking for L-BFGS.
    #  @param self Instance of the class.
    #  @param X_orig Original feature subset matrix for current fold.
    #  @param w_orig Original sample weights subset for current fold.
    #  @param X_targ Target feature subset matrix for current fold.
    #  @param w_targ Target sample weights subset for current fold.
    #  @return Tuple `(models, final_ratios)` containing list of model dicts and final ratio array.
    # =========================================================================
    def _fit_iterative_stream ( self   ,
                                X_orig ,
                                w_orig ,
                                X_targ ,
                                w_targ ) :
        """ Execute iterative PyTorch boosting stream on fold data.
        
        Performs sequential weak learner optimization with early stopping, dynamic loss 
        tracking, weight normalization, and zero-mean spatial calibration.
        """
        import torch
        import torch.nn as nn
        import torch.optim as optim

        device     = self.params.get ( 'device' , 'cuda' if torch.cuda.is_available () else 'cpu' )        
        device     = torch.device ( device )
        
        n_orig     = len ( X_orig )
        n_targ     = len ( X_targ )
        n_features = X_orig.shape [ 1 ]

        X_comb     = numpy.vstack ( [ X_orig, X_targ ] ).astype ( numpy.float32 )
        abs_w_orig = numpy.abs ( w_orig )
        abs_w_targ = numpy.abs ( w_targ )
        
        sum_w_orig = numpy.sum ( abs_w_orig )
        sum_w_targ = numpy.sum ( abs_w_targ )
        norm_targ  = w_targ * ( sum_w_orig / ( sum_w_targ + 1e-12 ) )

        y_comb         = numpy.hstack ( [ numpy.zeros ( n_orig, dtype = numpy.float32 ),
                                          numpy.ones  ( n_targ, dtype = numpy.float32 ) ] )
        current_w_orig = abs_w_orig.copy ()
        
        n_estimators      = self.params.get ( 'n_estimators'       , TORCH_ESTIMATORS      )
        learning_rate     = self.params.get ( 'learning_rate'      , TORCH_LEARNING_RATE   )
        hidden_layers     = self.params.get ( 'hidden_layers'      , TORCH_HIDDEN_LAYERS   )
        batch_size        = self.params.get ( 'batch_size'         , TORCH_BATCH_SIZE      )
        lbfgs_chunk       = self.params.get ( 'lbfgs_chunk'        , TORCH_LBFGS_CHUNK     )
        weight_decay      = self.params.get ( 'weight_decay'       , TORCH_WEIGHT_DECAY    )
        
        optimizer_type    = self.params.get ( 'optimizer'          , TORCH_OPTIMIZER       )
        optimizer_type    = optimizer_type.strip().lower().replace('-','')        
        if not optimizer_type in TORCH_OPTIMIZERS :
            raise TypeError ( "Unknown `optimizer': %s" % optimizer_type ) 

        optimizer_learning_rate = self.params.get ( 'optimizer_learning_rate' , 0.01 )
        print_freq        = self.params.get ( 'print_freq'         , PRINT_FREQ )

        epochs_per_step   = self.params.get ( 'epochs_per_step'    , TORCH_EPOCHS_PER_STEP )
        min_epoches       = self.params.get ( 'min_epoches'        , 5                     )
        min_delta         = self.params.get ( 'min_delta'          , 1e-4                  )
        patience          = self.params.get ( 'patience'           , 5                     )
        
        models        = []
        X_orig_tensor = torch.tensor ( X_orig.astype ( numpy.float32 ), dtype = torch.float32 ).to ( device )
        X_comb_tensor = torch.tensor ( X_comb, dtype = torch.float32 ).to ( device )
        y_comb_tensor = torch.tensor ( y_comb, dtype = torch.float32 ).unsqueeze ( 1 ).to ( device )

        WeakMLPClass = self.get_MLP_class ()

        neff0       = nEff ( X_orig, current_w_orig )        
        description = boosting_symbol + ' ' + S.fixed_width ( self.symbol , 3 ) + ':'
        
        if   self.progress : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = not self.progress ) 
        elif self.verbose  : the_bar = text_bar     ( n_estimators , full = True , description = description , silent =     self.silent   )
        else               : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = True  )
        
        t_start   = time.time ()
        prev_loss = 0  
        for step , index , bar in safe_loop ( the_bar ) : 
            
            comb_weights   = numpy.hstack ( [ current_w_orig, numpy.abs ( norm_targ ) ] )
            comb_weights  /= ( numpy.sum ( comb_weights ) + 1e-12 )
            weights_tensor = torch.tensor ( comb_weights, dtype = torch.float32 ).unsqueeze ( 1 ).to ( device )

            model     = WeakMLPClass ( n_features, hidden_layers ).to ( device )
            criterion = nn.BCEWithLogitsLoss ( reduction = 'none' )
            model.train ()

            if optimizer_type == 'lbfgs' :
                
                optimizer = optim.LBFGS ( model.parameters ()                      ,
                                          lr             = optimizer_learning_rate ,
                                          max_iter       = epochs_per_step         ,
                                          line_search_fn = 'strong_wolfe'          )
                
                if lbfgs_chunk < n_orig :
                    orig_idx = numpy.random.choice ( n_orig, size = lbfgs_chunk , replace = False )
                    t_size   = min ( n_targ, int ( lbfgs_chunk * n_targ / n_orig ) )
                    targ_idx = numpy.random.choice ( n_targ, size = t_size, replace = False )
                    
                    X_chunk = numpy.vstack ( [ X_orig[orig_idx], X_targ[targ_idx] ] ).astype ( numpy.float32 )
                    y_chunk = numpy.hstack ( [ numpy.zeros ( len(orig_idx), dtype=numpy.float32 ), numpy.ones ( len(targ_idx), dtype=numpy.float32 ) ] )
                    w_chunk = numpy.hstack ( [ current_w_orig[orig_idx], numpy.abs ( norm_targ[targ_idx] ) ] )
                    w_chunk /= ( numpy.sum ( w_chunk ) + 1e-12 )
                else :
                    X_chunk = X_comb
                    y_chunk = y_comb
                    w_chunk = comb_weights

                X_chunk_tensor = torch.tensor ( X_chunk, dtype = torch.float32 ).to ( device )
                y_chunk_tensor = torch.tensor ( y_chunk, dtype = torch.float32 ).unsqueeze ( 1 ).to ( device )
                w_chunk_tensor = torch.tensor ( w_chunk, dtype = torch.float32 ).unsqueeze ( 1 ).to ( device )

                # 1. Compute initial loss on the chunk BEFORE L-BFGS optimization
                with torch.no_grad () :
                    preds_init = model ( X_chunk_tensor )
                    loss_init  = torch.sum ( criterion ( preds_init, y_chunk_tensor ) * w_chunk_tensor )
                    if 0 < weight_decay :
                        l2_reg_init = sum ( torch.sum ( p ** 2 ) for p in model.parameters() )
                        loss_init = loss_init + 0.5 * weight_decay * l2_reg_init
                    initial_loss = float ( loss_init.item () )

                def closure () :
                    optimizer.zero_grad ()
                    preds = model ( X_chunk_tensor )
                    loss  = torch.sum ( criterion ( preds, y_chunk_tensor ) * w_chunk_tensor )
                    if 0 < weight_decay :
                        l2_reg = sum ( torch.sum ( p ** 2 ) for p in model.parameters() )
                        loss = loss + 0.5 * weight_decay * l2_reg
                    loss.backward ()
                    return loss

                optimizer.step ( closure )
                
                # 2. Compute final loss on the chunk AFTER L-BFGS optimization
                with torch.no_grad () :
                    preds_final = model ( X_chunk_tensor )
                    loss_final  = torch.sum ( criterion ( preds_final, y_chunk_tensor ) * w_chunk_tensor )
                    if 0 < weight_decay :
                        l2_reg_final = sum ( torch.sum ( p ** 2 ) for p in model.parameters() )
                        loss_final = loss_final + 0.5 * weight_decay * l2_reg_final
                    running_loss = float ( loss_final.item () )
                    best_loss    = running_loss
                    nepoches     = 1
                    
            else :
                
                optimizer = optim.Adam  ( model.parameters ()                      ,
                                          lr             = optimizer_learning_rate ,
                                          weight_decay   = weight_decay            )
                
                dataset   = torch.utils.data.TensorDataset ( X_comb_tensor, y_comb_tensor, weights_tensor )
                loader    = torch.utils.data.DataLoader    ( dataset, batch_size = batch_size, shuffle = True )

                initial_loss     = None 
                running_loss     = 0.0
                best_loss        = pos_infinity 
                patience_counter = 0
                nepoches         = 0
                
                for epoch in safe_loop ( epochs_per_step ) :
                    nepoches   += 1 
                    epoch_loss  = 0.0
                    batches     = 0  
                    for bx, by, bw in loader :
                        optimizer.zero_grad ()
                        preds   = model ( bx )
                        bw_norm = bw / ( torch.sum ( bw ) + 1e-12 )
                        loss    = torch.sum ( criterion ( preds, by ) * bw_norm )
                        
                        current_loss = loss.item ()
                        loss.backward  ()
                        optimizer.step ()
                        
                        if initial_loss is None : 
                            initial_loss = float ( current_loss )
                                            
                        epoch_loss += current_loss 
                        batches    += 1
                        
                    running_loss = float ( epoch_loss / max ( batches, 1 ) ) 
                    
                    if min_epoches <= epoch :
                        if running_loss + min_delta < best_loss :
                            best_loss        =  running_loss
                            patience_counter =  0  
                        else:
                            patience_counter += 1 
                            
                        if patience <= patience_counter : break

            delta_loss = running_loss - prev_loss
            prev_loss  = running_loss

            model.eval ()
            with torch.no_grad () :
                raw_vals = model ( X_orig_tensor ).cpu ().numpy ().ravel ()

            # -----------------------------------------------------------------
            # ZERO-MEAN CENTERING
            # Extract the global shift to ensure the boosting step only applies 
            # spatial (feature-dependent) variance. Save this shift for inference.
            # -----------------------------------------------------------------
            train_shift = float ( numpy.mean ( raw_vals ) )
            raw_vals    = raw_vals - train_shift

            log_r_step      = numpy.clip ( raw_vals, *LOG_STEP_CLIPS )
            current_w_orig *= numpy.exp  ( log_r_step * learning_rate )

            if bar and not self.silent and ( step < 5  or ( step + 1 ) % print_freq == 0 or ( step + 5 ) == n_estimators ) :
                ratios_curr = current_w_orig / ( abs_w_orig + 1e-12 )
                neff        = nEff ( X_orig, current_w_orig )
                wmin        = nice_print ( float ( ratios_curr.min () ) , precision = 2, with_sign = False, strip = True )
                wmax        = nice_print ( float ( ratios_curr.max () ) , precision = 2, with_sign = False, strip = True )
                ne0         = nice_print ( float ( neff0              ) , precision = 2, with_sign = False, strip = True )
                ne          = nice_print ( float ( neff               ) , precision = 2, with_sign = False, strip = True )

                dt       = ( time.time () - t_start ) / ( step + 1 )
                if dt < 0.1 : dt = '%dms'  % ( math.ceil ( 1000 * dt ) )
                else        : dt = '%.1fs' %                      dt

                init_str  = nice_print ( initial_loss , precision = 3 , with_sign = False , strip = True )
                loss_str  = nice_print ( running_loss , precision = 3 , with_sign = False , strip = True )
                best_str  = nice_print ( best_loss    , precision = 3 , with_sign = False , strip = True ) if best_loss != running_loss else "="
                delta_str = nice_print ( delta_loss   , precision = 3 , with_sign = True  , strip = True )
                
                message   = f"Step {step+1:4d}/{n_estimators:<4d}"
                if best_str != loss_str : message  += f" | loss: {init_str} {S.arrow_right} {loss_str}/{best_str}"
                else                    : message  += f" | loss: {init_str} {S.arrow_right} {loss_str}"                
                message  += f" | {S.delta_symbol}: {delta_str}"
                if nepoches != epochs_per_step : message += f" | #epochs/step:{nepoches}"
                message  += f" | wmin/wmax: {wmin}/{wmax}"
                message  += f" | #eff: {ne0} {S.arrow_right} {ne} | Time/step:{dt}"
                
                bar      += message
    
            models.append ( { 
                'state_dict'    : model.state_dict () , 
                'input_dim'     : n_features          , 
                'hidden_layers' : hidden_layers       ,
                'shift'         : train_shift         
            } )
            
        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_clip ( final_ratios, w_orig )

        return models, final_ratios

    # =========================================================================
    ## @brief Accumulate raw logit predictions across sequential PyTorch models.
    #  @details Evaluates each trained PyTorch weak learner on unseen feature matrix X, 
    #           subtracts its stored spatial zero-centering shift, clips log-ratio steps 
    #           to [-7.0, 7.0], accumulates them scaled by learning rate, and returns 
    #           final exponential density ratio predictions.
    #  @param self Instance of the class.
    #  @param models List of trained PyTorch model state dicts comprising sequential chain.
    #  @param X Feature matrix of unseen samples to evaluate.
    #  @return Array of predicted density ratios r(x).
    # =========================================================================
    def _predict_iterative_stream ( self, models, X ) :
        """ Accumulate raw logit predictions across sequential PyTorch models.
        
        Evaluates saved weak learners, applies zero-mean spatial corrections, 
        and computes cumulative density ratios in log-exponential space.
        """
        import torch

        # Retrieve device dynamically from model parameters
        device        = self.params.get ( 'device', 'cuda' if torch.cuda.is_available () else 'cpu' )
        device        = torch.device ( device )
        
        learning_rate = self.params.get ( 'learning_rate', TORCH_LEARNING_RATE )
        log_ratios    = numpy.zeros ( len ( X ), dtype = numpy.float32 )
        X_tensor      = torch.tensor ( X.astype ( numpy.float32 ), dtype = torch.float32 ).to ( device )

        WeakMLPClass = self.get_MLP_class ()

        for m_dict in models :
            h_layers = m_dict.get ( 'hidden_layers', TORCH_HIDDEN_LAYERS  )
            
            model = WeakMLPClass ( m_dict [ 'input_dim' ], h_layers ).to ( device )
            model.load_state_dict ( m_dict [ 'state_dict' ] )
            model.eval ()

            with torch.no_grad () :
                raw_vals = model ( X_tensor ).cpu ().numpy ().ravel ()

            train_shift = m_dict.get ( 'shift', 0.0 )
            raw_vals    = raw_vals - train_shift
            
            log_r_step  = numpy.clip ( raw_vals, *LOG_STEP_CLIPS )
            log_ratios += log_r_step * learning_rate

        return numpy.exp ( log_ratios )


# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__, logger = logger )

    from ostap.stats.tools import ( hasLightGBM ,
                                    hasXGBoost  ,
                                    hasCatBoost ,
                                    hasPyTorch  ) 

    if not hasLightGBM ( False ) : logger.warning ( "No LightGBM available!" ) 
    if not hasXGBoost  ( False ) : logger.warning ( "No XGBoost  available!" ) 
    if not hasCatBoost ( False ) : logger.warning ( "No CatBoost available!" ) 
    if not hasPyTorch  ( False ) : logger.warning ( "No PyTorch  available!" ) 

# =============================================================================
#                                                                       The END 
# =============================================================================
