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
    'LightGBMIterativeReweighter'     , # LightGBM-based iterative density reweighter
    'XGBoostIterativeReweighter'      , # XGBoost-based  iterative density reweighter
    'CatBoostIterativeReweighter'     , # CatBoost-based iterative density reweighter
    'PyTorchAdamIterativeReweighter'  , # PyTorch-based Adam iterative density reweighter
    'PyTorchLBFGSIterativeReweighter' , # PyTorch-based L-BFGS iterative density reweighter
) 
# =============================================================================
from   ostap.core.ostap_types       import num_types, string_types, integer_types
from   ostap.utils.core             import typename
from   ostap.utils.basic            import safe_loop, LoopGC 
from   ostap.utils.utils            import classprop 
from   ostap.utils.memory           import use_gc 
from   ostap.logger.pretty          import nice_print, pretty_floats, format_pow10
from   ostap.math.math_base         import pos_infinity  
from   ostap.utils.progress_bar     import progress_bar, text_bar 
from   ostap.stats.utils            import ( num_samples        ,
                                             num_features       ,
                                             nEff               ,
                                             optimal_batch_size ) 
from   ostap.tools.base_reweighters import IterativeReweighter  
import ostap.stats.utils            as     U 
import ostap.logger.symbols         as     S
import numpy, math, time, pickle, gc 
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
boosting_symbol = S.ladder if S.show and S.ladder  else 'Boosting'
# =============================================================================
# Global Configuration Constants
# =============================================================================
LGBM_ESTIMATORS       =  300 # High capacity: deep ensemble with low learning rate
LGBM_MAX_DEPTH        =    6 # Deep tree capacity for high dimensions / large statistics
LGBM_LEARNING_RATE    = 0.02 # Small step for stable density ratio convergence

XGB_ESTIMATORS        =  300 # High capacity: deep ensemble with low learning rat
XGB_MAX_DEPTH         =    6 # Deep tree capacity for high dimensions / large statistic
XGB_LEARNING_RATE     = 0.02 # Small step for stable density ratio convergence
XGB_MAX_LEAVES        =   63 # 

CATB_ESTIMATORS       =  300 # High capacity: deep ensemble with low learning rate
CATB_MAX_DEPTH        =    6 # Deep tree capacity for high dimensions / large statistics
CATB_LEARNING_RATE    = 0.02 # Small step for stable density ratio convergence

TORCH_ESTIMATORS      =  300 # Number of boosting steps (epochs/models)
TORCH_LEARNING_RATE   = 0.02 # Step size for weight boosting
TORCH_WEIGHT_DECAY    = 1e-5 # regularization parameter  
TORCH_BATCH_SIZE      = 'auto'
TORCH_HIDDEN_LAYERS   = 32 ,  # ## single hiddle layer with 32 nodes 
TORCH_MAX_BATCH_SIZE  = 2 ** 16

TORCH_ADAM_EPOCHS_PER_STEP  = 100 # 
TORCH_LBFGS_EPOCHS_PER_STEP =  25 # 

TORCH_OPTIMIZER_ADAM  = 'adam'
TORCH_OPTIMIZER_LBFGS = 'lbfgs'
TORCH_OPTIMIZERS      = TORCH_OPTIMIZER_ADAM , TORCH_OPTIMIZER_LBFGS
TORCH_OPTIMIZER       = TORCH_OPTIMIZER_LBFGS 
TORCH_LBFGS_CHUNK     = 200_000
TORCH_LBFGS_SYMBOL    = S.gear     
TORCH_ADAM_SYMBOL     = S.lightning

# =============================================================================
THRESHOLD_NSAMPLES    = 50000
THRESHOLD_NFEATURES   = 10

PRINT_FREQ            = 10 
# =============================================================================
LOG_STEP_CLIPS        = -7.0 , 7.0 
# =============================================================================
## @class LightGBMIterativeReweighter
#  Iterative density ratio reweighter based on LightGBM gradient boosting trees.
#
#  Implements sequential gradient boosting for density ratio estimation 
#           \f$r(x) = \frac{p(x)}{1-p(x)}\f$ using raw margin predictions (log-odds)
#           with dynamically updated sample weights on each boosting iteration.
# =============================================================================
class LightGBMIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using LightGBM backend.
    """

    # =========================================================================
    ## Initialize the LightGBM iterative density reweighter instance.
    #  @param original Source feature matrix (N-dimensional array or dataframe).
    #  @param target Target feature matrix (N-dimensional array or dataframe).
    #  @param original_weight Optional observation weights for the source dataset.
    #  @param target_weight Optional observation weights for the target dataset.
    #  @param kwargs Additional keyword arguments for LightGBM configuration.
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
            ## 
            'lambda_l1'               :  0.0                ,
            'lambda_l2'               :  0.0                ,                        
            ##
            'bagging_freq'            :  0                  ,
            'subsample'               :  1.0                ,
            'feature_fraction'        :  1.0                ,
            ##
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
    ## Get the short symbol/tag representing the LightGBM reweighter method.
    #  @return String identifier symbol used for logging and console visual display.
    @classprop
    def symbol ( self ) :
        """ Short symbol/tag for reweighter.
        """
        return U.symbol_LGBM
    
    # =========================================================================
    ## Get the full method identifier name.
    #  @return Method string tag.
    @property
    def method ( self ) :
        """ Return method identifier name.
        """
        return method_LGBM 
    
    # =========================================================================
    ## Execute the iterative LightGBM boosting stream on a specific cross-validation fold.
    #  @details Sequentially fits LightGBM weak learners (1 tree per step), evaluates raw margins 
    #           (log-odds), applies log-space clipping to [-7.0, 7.0], and updates sample weights 
    #           via exponential scaling.
    #  @param X_orig Original feature subset matrix for current fold.
    #  @param w_orig Original sample weights subset for current fold.
    #  @param X_targ Target feature subset matrix for current fold.
    #  @param w_targ Target sample weights subset for current fold.
    #  @return Tuple `(models, final_ratios)` containing trained boosters and final density ratio array.
    def _fit_iterative_stream ( self   ,
                                X_orig ,
                                w_orig ,
                                X_targ ,
                                w_targ ) :
        """ Execute iterative LightGBM boosting stream on fold data.
        """
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
        
        abs_norm_targ = numpy.abs   ( norm_targ ).astype ( numpy.float32 , copy=False)
        sum_abs_targ  = numpy.sum   ( abs_norm_targ   , dtype = numpy.float64 )
        comb_weights  = numpy.empty ( n_orig + n_targ , dtype = numpy.float32 )

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

        description = S.fixed_width ( boosting_symbol , 2 ) + S.fixed_width ( self.symbol , 2 ) + ' :'        
        if   self.progress : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = not self.progress ) 
        elif self.verbose  : the_bar = text_bar     ( n_estimators , full = True , description = description , silent =     self.silent   )
        else               : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = True  ) 

        t_start = time.time () 

        for it , _ ,  bar in LoopGC ( safe_loop ( the_bar ) ) :

            w_sum    = numpy.sum ( current_w_orig , dtype=numpy.float64 ) + sum_abs_targ + 1e-12
            inv_sum  = numpy.float32 ( 1.0 / w_sum )
            comb_weights [ : n_orig ] = current_w_orig * inv_sum
            comb_weights [ n_orig : ] = abs_norm_targ  * inv_sum
            
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

                wmin    , wmax  = float ( ratios_curr.min () ) , float ( ratios_curr.max () ) 
                wminmax , wexpo = pretty_floats ( wmin , wmax , precision = 2 , with_sign = False , strip = True )
                wminmax         = '/'.join ( wminmax )
                if wexpo : wminmax = '(%s)%s%s' % ( wminmax , S.times , format_pow ( wexpo ) )
                
                ne0         = nice_print ( float ( neff0              ) , precision = 2, with_sign = False, strip = True )
                ne          = nice_print ( float ( neff               ) , precision = 2, with_sign = False, strip = True )

                dt       = ( time.time () - t_start ) / ( it + 1 )
                if dt < 1 : dt = '%dms'  % ( math.ceil ( 1000 * dt ) )
                else      : dt = '%.1fs' %                      dt                 
                message  = f"Step {it+1:4d}/{n_estimators:<4d} | wmin/wmax: {wminmax} | #eff: {ne0} {S.arrow_right} {ne} | Time/step:{dt}" 
                bar     += message

                ## clean-up 
                del ratios_curr
                
            if self.phoenix :
                
                model_bytes = model.model_to_string ()
                model.free_dataset ()                
                del model
                model = LightGBM.Booster ( model_str = model_bytes )
                del model_bytes
                
            models.append ( model )
            
            ## clean-up 
            del trn_data
            del raw_vals
            del log_r_step
            
        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_clip ( final_ratios, w_orig )

        ## more cleanup 
        del X_comb, y_comb 
        del current_w_orig
        del abs_norm_targ
        
        return models, final_ratios

    # =========================================================================
    ## Accumulate raw margin predictions and compute log-ratios across sequential LightGBM models.
    #  @details Evaluates each trained LightGBM booster on feature matrix X in raw log-odds space,
    #           clips log-ratio steps to [-7.0, 7.0], accumulates them scaled by learning rate,
    #           and returns final exponential density ratio predictions.
    #  @param models List of trained LightGBM booster models comprising sequential chain.
    #  @param X Feature matrix of unseen samples to evaluate.
    #  @return Array of predicted density ratios r(x).
    # =========================================================================
    def _predict_iterative_stream ( self, models, X ) :
        """ Accumulate predictions across sequential LightGBM models.
        """
        
        learning_rate = self.params.get ( 'learning_rate', LGBM_LEARNING_RATE )
        log_ratios    = numpy.zeros ( len ( X ), dtype = numpy.float32 )

        for model in models :
            raw_vals    = model.predict ( X, raw_score = True )
            if 1 < raw_vals.ndim : raw_vals = raw_vals [ :, 0 ]
            
            log_r_step  = numpy.clip ( raw_vals, *LOG_STEP_CLIPS )
            
            log_r_step *= learning_rate
            log_ratios += log_r_step
            
            del log_r_step

        rexp = numpy.exp ( log_ratios )
            
        del log_ratios
        
        return rexp

# =============================================================================
## @class XGBoostIterativeReweighter
#  Iterative weight boosting density reweighter using XGBoost backend.
#
#  Implements sequential gradient boosting for density ratio estimation 
#           \f$r(x) = \frac{p(x)}{1-p(x)}\f$ using raw margin predictions (output_margin)
#           with dynamically updated sample weights on each boosting iteration.
# =============================================================================
class XGBoostIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using XGBoost backend.
    """

    # =========================================================================
    ## Initialize the XGBoost iterative density reweighter instance.
    #  @param original Source feature matrix (N-dimensional array or dataframe).
    #  @param target Target feature matrix (N-dimensional array or dataframe).
    #  @param original_weight Optional observation weights for the source dataset.
    #  @param target_weight Optional observation weights for the target dataset.
    #  @param kwargs Additional keyword arguments for XGBoost configuration.
    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        """ Initialize XGBoost iterative density reweighter with custom defaults.
        """
        config = {
            'n_estimators'      :  XGB_ESTIMATORS    ,
            'learning_rate'     :  XGB_LEARNING_RATE ,
            'max_depth'         :  XGB_MAX_DEPTH     ,
            'tree_method'       :  'hist'            , 
            'max_bin'           :  1023              , 
            'grow_policy'       :  'depthwise'       , 
            'max_leaves'        :  XGB_MAX_LEAVES    , 
            'min_child_weight'  :  1e-4              , 
            'max_delta_step'    :  0.0               ,
            ## 
            'gamma'             :  0.0               , 
            'reg_alpha'         :  0.0               , 
            'reg_lambda'        :  1e-4              ,
            ##
            'subsample'         :  1.0               , 
            'colsamples_bytree' :  1.0               , 
            'verbosity'         :  0                 ,
            'print_freq'        : PRINT_FREQ         , 
            'n_jobs'            :  1                 ,
            'n_splits'          :  5                 ,
            
        }

        config.update ( kwargs )
        if 'silent' in config and config.get ( 'silent' ) : config [ 'verbosity' ] = 0        
        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round', XGB_ESTIMATORS ) 
        
        n_rows = num_samples ( original ) + num_samples ( target )
        if n_rows < THRESHOLD_NSAMPLES or 1 < config.get ( 'n_splits' , 5 ) : config.pop ( 'subsample'         , None ) 
        if num_features      ( original ) < THRESHOLD_NFEATURES             : config.pop ( 'colsamples_bytree' , None )

        if n_rows < 10000 : config [ 'tree_method' ]  = 'exact'

        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **config )

    # =========================================================================
    ## Get the short symbol/tag representing the XGBoost reweighter method.
    #  @return String identifier symbol used for logging and visual presentation.
    @classprop
    def symbol ( self ) :
        """ Short symbol/tag for reweighter.
        """
        return U.symbol_XGB

    # =========================================================================
    ## Get the full method identifier name.
    #  @return Method string tag.
    @property
    def method ( self ) :
        """ Return method identifier name.
        """
        return method_XGB

    # =========================================================================
    ## Execute the iterative XGBoost boosting stream on a specific cross-validation fold.
    #
    #  Sequentially fits XGBoost weak learners (1 tree per step), evaluates raw margins 
    #           (output_margin), applies log-space clipping to [-7.0, 7.0], and updates sample weights 
    #           via exponential scaling.
    #  @param X_orig Original feature subset matrix for current fold.
    #  @param w_orig Original sample weights subset for current fold.
    #  @param X_targ Target feature subset matrix for current fold.
    #  @param w_targ Target sample weights subset for current fold.
    #  @return Tuple `(models, final_ratios)` containing trained boosters and final density ratio array.
    def _fit_iterative_stream ( self   ,
                                X_orig ,
                                w_orig ,
                                X_targ ,
                                w_targ ) :
        """ Execute iterative XGBoost boosting stream on fold data.
        """
        import xgboost as XGBoost

        n_orig     = len ( X_orig )
        n_targ     = len ( X_targ )
        X_comb     = numpy.vstack ( [ X_orig, X_targ ] )

        dmat_orig = XGBoost.DMatrix ( X_orig )

        abs_w_orig = numpy.abs ( w_orig )
        abs_w_targ = numpy.abs ( w_targ )
        
        sum_w_orig = numpy.sum ( abs_w_orig )
        sum_w_targ = numpy.sum ( abs_w_targ )
        norm_targ  = w_targ * ( sum_w_orig / ( sum_w_targ + 1e-12 ) )

        y_comb         = numpy.hstack ( [ numpy.zeros ( n_orig, dtype = numpy.float32 ),
                                          numpy.ones  ( n_targ, dtype = numpy.float32 ) ] )
        current_w_orig = abs_w_orig.copy ()
        
        abs_norm_targ = numpy.abs   ( norm_targ ).astype ( numpy.float32 , copy=False)
        sum_abs_targ  = numpy.sum   ( abs_norm_targ   , dtype = numpy.float64 )
        comb_weights  = numpy.empty ( n_orig + n_targ , dtype = numpy.float32 )

        n_estimators  = self.params.get ( 'n_estimators' , XGB_ESTIMATORS    )
        learning_rate = self.params.get ( 'learning_rate', XGB_LEARNING_RATE )
        
        dmat_comb = XGBoost.DMatrix ( X_comb , label = y_comb )
        models    = []

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
        
        description = S.fixed_width ( boosting_symbol , 2 ) + S.fixed_width ( self.symbol , 2 ) + ' :'        
        if   self.progress : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = not self.progress ) 
        elif self.verbose  : the_bar = text_bar     ( n_estimators , full = True , description = description , silent =     self.silent   )
        else               : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = True  ) 

        t_start = time.time () 
        for it , _ ,  bar in LoopGC ( safe_loop ( the_bar ) ) :

            w_sum    = numpy.sum ( current_w_orig , dtype=numpy.float64 ) + sum_abs_targ + 1e-12
            inv_sum  = numpy.float32 ( 1.0 / w_sum )
            comb_weights [ : n_orig ] = current_w_orig * inv_sum
            comb_weights [ n_orig : ] = abs_norm_targ  * inv_sum

            ## just change the weights  
            dmat_comb.set_weight ( comb_weights )
            
            model     = XGBoost.train ( params          = params,
                                        dtrain          = dmat_comb,
                                        num_boost_round = 1 )
            
            # Retrieve raw margins (output_margin=True) for exact mathematical consistency
            raw_vals  = model.predict ( dmat_orig, output_margin = True )
            if 1 < raw_vals.ndim : raw_vals = raw_vals [ :, 0 ]
            
            log_r_step      = numpy.clip ( raw_vals, *LOG_STEP_CLIPS )
            current_w_orig *= numpy.exp ( log_r_step * learning_rate )
        
            if bar and not self.silent and ( it < 5  or ( it + 1 ) % print_freq == 0 or ( it + 5 ) == n_estimators ) :
                ratios_curr = current_w_orig / ( abs_w_orig + 1e-12 )
                neff        = nEff ( X_orig, current_w_orig )

                wmin    , wmax  = float ( ratios_curr.min () ) , float ( ratios_curr.max () ) 
                wminmax , wexpo = pretty_floats ( wmin , wmax , precision = 2 , with_sign = False , strip = True )
                wminmax         = '/'.join ( wminmax )
                if wexpo : wminmax = '(%s)%s%s' % ( wminmax , S.times , format_pow ( wexpo ) )

                ne0         = nice_print ( float ( neff0              ) , precision = 2, with_sign = False, strip = True )
                ne          = nice_print ( float ( neff               ) , precision = 2, with_sign = False, strip = True )

                dt       = ( time.time () - t_start ) / ( it + 1 )
                if dt < 1 : dt = '%dms'  % ( math.ceil ( 1000 * dt ) )
                else      : dt = '%.1fs' %                      dt                 
                message  = f"Step {it+1:4d}/{n_estimators:<4d} | wmin/wmax: {wminmax} | #eff: {ne0} {S.arrow_right} {ne} | Time/step:{dt}" 
                bar     += message

                ## cleanup
                del ratios_curr
                
            # PHOENIX memory optimization: export raw JSON model bytes and reload clean inference engine
            if self.phoenix : 
                model_bytes = model.save_raw ( raw_format = "json" )
                best_iter   = getattr ( model, 'best_iteration_', None ) or getattr ( model, 'best_iteration', None )                    
                del model
                model = XGBoost.Booster ()
                model.load_model ( model_bytes )
                del model_bytes            
                if best_iter is not None : model.best_iteration = best_iter
                
            models.append ( model )

            ## cleanup 
            del raw_vals
            del log_r_step 

        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_clip ( final_ratios, w_orig )
        
        ## more cleanup
        del dmat_orig 
        del dmat_comb
        del X_comb, y_comb
        del current_w_orig
        del abs_norm_targ
        
        return models, final_ratios

    # =========================================================================
    ## @brief Accumulate raw margin predictions across sequential XGBoost models.
    #  @details Evaluates each trained XGBoost booster on feature matrix X in raw log-odds space,
    #           clips log-ratio steps to [-7.0, 7.0], accumulates them scaled by learning rate,
    #           and returns final exponential density ratio predictions.
    #  @param models List of trained XGBoost booster models comprising sequential chain.
    #  @param X Feature matrix of unseen samples to evaluate.
    #  @return Array of predicted density ratios r(x).
    # =========================================================================
    def _predict_iterative_stream ( self, models, X ) :
        """ Accumulate predictions across sequential XGBoost models.
        """
        import xgboost as XGBoost

        with use_gc () :
            
            dmat          = XGBoost.DMatrix ( X )
            learning_rate = self.params.get ( 'learning_rate', XGB_LEARNING_RATE )
            log_ratios    = numpy.zeros ( len ( X ), dtype = numpy.float32 )
            
            for model in models :
                
                raw_vals    = model.predict ( dmat, output_margin = True )
                if 1 < raw_vals.ndim : raw_vals = raw_vals [ :, 0 ]
                
                log_r_step  = numpy.clip ( raw_vals, *LOG_STEP_CLIPS )
                log_r_step *= learning_rate
                log_ratios += log_r_step
                
                del log_r_step
                del raw_vals
                
            rexp = numpy.exp ( log_ratios )
            
            del dmat 
            del log_ratios
            
            return rexp


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
    #  @param original Source feature matrix (N-dimensional array or dataframe).
    #  @param target Target feature matrix (N-dimensional array or dataframe).
    #  @param original_weight Optional observation weights for the source dataset.
    #  @param target_weight Optional observation weights for the target dataset.
    #  @param kwargs Additional keyword arguments for CatBoost configuration.
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
            'l2_leaf_reg'                :  1e-3               ,
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

        if 'bayesian' == config.get ( 'bootstrap_type' , '' ).strip().lower()  :
            config.pop ( 'subsample'  , None )
            config.pop ( 'rsm'        , None )
                    
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
    ## Execute the iterative CatBoost boosting stream on a specific cross-validation fold.
    #  @details Sequentially fits CatBoost weak learners (1 tree per step), evaluates raw margins 
    #           (RawFormulaVal), applies log-space clipping to [-7.0, 7.0], and updates sample weights 
    #           via exponential scaling.
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
        """ Execute iterative CatBoost boosting stream on fold data.
        """
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

        abs_norm_targ = numpy.abs   ( norm_targ ).astype ( numpy.float32 , copy=False)
        sum_abs_targ  = numpy.sum   ( abs_norm_targ   , dtype = numpy.float64 )
        comb_weights  = numpy.empty ( n_orig + n_targ , dtype = numpy.float32 )
 
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

        description = S.fixed_width ( boosting_symbol , 2 ) + S.fixed_width ( self.symbol , 2 ) + ' :'        
        if   self.progress : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = not self.progress ) 
        elif self.verbose  : the_bar = text_bar     ( n_estimators , full = True , description = description , silent =     self.silent   )
        else               : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = True  ) 

        t_start = time.time () 

        for it , _ ,  bar in LoopGC ( safe_loop ( the_bar ) ) :
            
            w_sum    = numpy.sum ( current_w_orig , dtype=numpy.float64 ) + sum_abs_targ + 1e-12
            inv_sum  = numpy.float32 ( 1.0 / w_sum )
            comb_weights [ : n_orig ] = current_w_orig * inv_sum
            comb_weights [ n_orig : ] = abs_norm_targ  * inv_sum
            
            trn_pool = CatBoost.Pool ( data = X_comb, label = y_comb, weight = comb_weights )
            
            model    = CatBoost.train ( pool    = trn_pool ,
                                        params  = params   ,
                                        verbose = False    )

            # Retrieve raw margins (RawFormulaVal) for exact mathematical consistency
            raw_vals = model.predict ( pool_orig , prediction_type = 'RawFormulaVal' )
            if 1 < raw_vals.ndim : raw_vals = raw_vals [ :, 0 ]
            
            log_r_step      = numpy.clip ( raw_vals   , *LOG_STEP_CLIPS )
            current_w_orig *= numpy.exp  ( log_r_step * learning_rate   )

            if bar and not self.silent and ( it < 5  or ( it + 1 ) % print_freq == 0 or ( it + 5 ) == n_estimators ) :
                
                ratios_curr = current_w_orig / ( abs_w_orig + 1e-12 )
                neff        = nEff ( X_orig, current_w_orig )

                wmin    , wmax  = float ( ratios_curr.min () ) , float ( ratios_curr.max () ) 
                wminmax , wexpo = pretty_floats ( wmin , wmax , precision = 2 , with_sign = False , strip = True )
                wminmax         = '/'.join ( wminmax )
                if wexpo : wminmax = '(%s)%s%s' % ( wminmax , S.times , format_pow ( wexpo ) )

                ne0         = nice_print ( float ( neff0              ) , precision = 2, with_sign = False, strip = True )
                ne          = nice_print ( float ( neff               ) , precision = 2, with_sign = False, strip = True )

                dt       = ( time.time () - t_start ) / ( it + 1 )
                if dt < 1 : dt = '%dms'  % ( math.ceil ( 1000 * dt ) )
                else      : dt = '%.1fs' %                      dt                 
                message  = f"Step {it+1:4d}/{n_estimators:<4d} | wmin/wmax: {wminmax} | #eff: {ne0} {S.arrow_right} {ne} | Time/step:{dt}" 
                bar     += message

                del ratios_curr

            # PHOENIX memory optimization: serialize model blob and purge heavy training pools
            if self.phoenix :  
                model_bytes = pickle.dumps ( model       )                        
                model       = pickle.loads ( model_bytes )
                del model_bytes
            
            models.append ( model )
            
            ## some clean up
            del trn_pool
            del raw_vals
            del log_r_step
                
            
        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_clip ( final_ratios, w_orig )

        ## more cleanup 
        del pool_orig
        del X_comb, y_comb    
        del current_w_orig
        del abs_norm_targ
        
        return models, final_ratios

    # =========================================================================
    ## @brief Accumulate raw margin predictions across sequential CatBoost models.
    #  @details Evaluates each trained CatBoost booster on feature matrix X in raw log-odds space
    #           (RawFormulaVal), clips log-ratio steps to [-7.0, 7.0], accumulates them scaled 
    #           by learning rate, and returns final exponential density ratio predictions.
    #  @param models List of trained CatBoost booster models comprising sequential chain.
    #  @param X Feature matrix of unseen samples to evaluate.
    #  @return Array of predicted density ratios r(x).
    # =========================================================================
    def _predict_iterative_stream ( self, models, X ) :
        """ Accumulate predictions across sequential CatBoost models.
        """
        import catboost as CatBoost
        
        learning_rate = self.params.get ( 'learning_rate', CATB_LEARNING_RATE )

        with use_gc () :
            
            pool          = CatBoost.Pool ( data = X )
            log_ratios    = numpy.zeros ( len ( X ), dtype = numpy.float32 )
            
            for model in models :
                raw_vals    = model.predict ( pool, prediction_type = 'RawFormulaVal' )
                if 1 < raw_vals.ndim : raw_vals = raw_vals [ :, 0 ]
                
                log_r_step  = numpy.clip ( raw_vals , *LOG_STEP_CLIPS )
                log_r_step *= learning_rate
                log_ratios += log_r_step
                
                del log_r_step
                del raw_vals

            rexp = numpy.exp ( log_ratios )
            
            del pool
            del log_ratios
            
            return rexp

# =============================================================================
## @class PyTorchBaseIterativeReweighter
#  (Abstract) base class for PyTorch-based iterative weight boosting
#            density ratio reweighters.
#
#  -         Implements the shared sequential gradient-boosting pipeline for
#           estimating density ratios \f$r(x) = \frac{p(x)}{1 - p(x)}\f$ using
#           feature-dependent Multi-Layer Perceptron (MLP) weak learners
#           At each boosting step, a fresh weak learner is trained on the
#           combined original and target datasets with dynamically updated
#           source weights, evaluated on the source subset, zero-mean centered
#           so that the boosting update acts purely on spatial variation,
#           clipped in log space via @c LOG_STEP_CLIPS, and accumulated
#           multiplicatively into @c current_w_orig
#
#           Concrete optimizer backends (such as first-order stochastic
#           @c Adam or quasi-Newton @c L-BFGS) are encapsulated in subclasses
#           by overriding @c _fit_step and @c optimizer_symbol
class PyTorchBaseIterativeReweighter ( IterativeReweighter  ) :
    """ (Abstract) base class for PyTorch-based iterative weight boosting reweighters.

    Encapsulates device selection, MLP weak-learner construction, sequential
    boosting iterations, zero-mean spatial centering, log-step clipping,
    progress reporting, model state serialization, and out-of-sample
    inference across saved weak learners. Subclasses implement
    optimizer-specific single-step training via ``_fit_step``.
    """

    # =========================================================================
    ## Initialize the PyTorch iterative density ratio reweighter base.
    #  - Merges caller-supplied keyword arguments with default PyTorch
    #           boosting parameters, resolves legacy option aliases
    #           (@c num_boost_round and @c hidden_dim), normalizes the execution
    #           device to a @c torch.device instance (CPU or CUDA), and invokes
    #           the parent @c IterativeReweighter initializer
    #
    #  @param original Source feature matrix (N-dimensional array or dataframe)
    #  @param target Target feature matrix (N-dimensional array or dataframe)
    #  @param original_weight Optional observation weights for the source dataset
    #  @param target_weight Optional observation weights for the target dataset
    #  @param kwargs Additional keyword arguments controlling boosting iterations,
    #         learning rate, network architecture, regularization, and device
    def __init__ ( self            ,
                   original        ,
                   target          ,
                   original_weight = None ,
                   target_weight   = None ,
                   **kwargs        ) :
        """ Initialize PyTorch iterative density ratio reweighter base configuration.

        Resolves aliases, selects and normalizes the target ``torch.device``,
        and delegates cross-validation and stream fitting to ``IterativeReweighter``.
        """
        import torch

        config = {
            'n_estimators'  : TORCH_ESTIMATORS    ,
            'learning_rate' : TORCH_LEARNING_RATE ,
            'hidden_layers' : TORCH_HIDDEN_LAYERS ,
            'weight_decay'  : TORCH_WEIGHT_DECAY  ,
            'device'        : None                , # Auto-detect CUDA GPU or fallback to CPU
            'verbosity'     : 0                   ,
            'print_freq'    : PRINT_FREQ          ,
            'n_jobs'        : 1                   ,
        }

        config.update ( kwargs )

        # Preserve backward compatibility with legacy parameter names used elsewhere in Ostap
        if 'num_boost_round' in config : config [ 'n_estimators'  ] = config.pop ( 'num_boost_round' , TORCH_ESTIMATORS    )
        if 'hidden_dim'      in config : config [ 'hidden_layers' ] = config.pop ( 'hidden_dim'      , TORCH_HIDDEN_LAYERS )

        # Normalize device specification into an explicit torch.device object
        device = config.get ( 'device', None )
        if device is None : config [ 'device' ] = torch.device ( 'cuda' if torch.cuda.is_available () else 'cpu' )
        elif isinstance ( device, string_types ) and device.strip ().lower () in ( 'gpu', 'cuda' ) :
            config [ 'device' ] = torch.device ( 'cuda' if torch.cuda.is_available () else 'cpu' )
        else :
            config [ 'device' ] = torch.device ( device )

        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **config )

    # =========================================================================
    ## Get the short symbol/tag representing the PyTorch reweighter method.
    #  @return String identifier symbol used for logging and console visual display
    @classprop
    def symbol ( self ) :
        """ Return short symbol/tag for PyTorch reweighter console display
        """
        return U.symbol_TORCH

    # =========================================================================
    ## Get the full method identifier name.
    #  @return Method string tag used in summaries and warnings
    @property
    def method ( self ) :
        """ Return method identifier name for PyTorch iterative reweighter
        """
        return method_TORCH

    # =========================================================================
    ## Get the symbol/tag representing the concrete optimizer backend.
    #
    #  Subclasses override this property to return the visual tag for
    #  their optimization algorithm (e.g. L-BFGS or Adam)
    #  @return Optimizer symbol string used in progress bar descriptions
    @property
    @abc.abstractmethod
    def optimizer_symbol ( self ) :
        """ Return optimizer symbol/tag for progress bar description
        """
        raise NotImplementedError

    # =========================================================================
    ## Factory method generating a PyTorch Module class for weak learners.
    #  - Builds a @c WeakMLP subclass of @c torch.nn.Module that first
    #           normalizes input features via @c BatchNorm1d, passes them through
    #           a sequence of @c Linear + @c LeakyReLU layers, and terminates in
    #           a final @c Linear layer with @c bias=False so that predictions
    #           capture relative spatial structure rather than a constant shift.
    #  @return Subclass of @c torch.nn.Module implementing @c WeakMLP.
    def get_MLP_class ( self ) :
        """ Construct and return the ``WeakMLP`` ``torch.nn.Module``.

        Creates a lightweight neural network with input ``BatchNorm1d``,
        ``LeakyReLU`` hidden activations, and a bias-free linear output layer
        """
        import torch

        class WeakMLP ( torch.nn.Module ) :
            """ Weak Multi-Layer Perceptron neural network module for density steps
            """

            def __init__ ( self, input_dim, hidden_layers = TORCH_HIDDEN_LAYERS ) :
                """ Initialize BatchNorm input stage, hidden blocks, and bias-free head
                """
                super ().__init__ ()

                # Accept either a single integer layer width or a sequence/tuple of widths.
                if isinstance ( hidden_layers, integer_types ) :
                    hidden_layers = ( hidden_layers, )

                # BatchNorm1d stabilizes feature scales across heterogeneous observables
                self.bn = torch.nn.BatchNorm1d ( input_dim )

                layers   = []
                curr_dim = input_dim
                for h_dim in hidden_layers :
                    linear = torch.nn.Linear ( curr_dim, h_dim )
                    torch.nn.init.kaiming_uniform_ ( linear.weight, a = 0.1 )
                    torch.nn.init.zeros_ ( linear.bias )

                    layers.append ( linear )
                    layers.append ( torch.nn.LeakyReLU ( 0.1 ) )
                    curr_dim = h_dim

                self.hidden_net = torch.nn.Sequential ( *layers )

                # bias=False prevents learning a global constant offset inside a single weak step;
                # any residual mean shift is explicitly removed via zero-mean centering
                self.fc_out = torch.nn.Linear ( curr_dim, 1, bias = False )
                torch.nn.init.xavier_uniform_ ( self.fc_out.weight )

            def forward ( self, x ) :
                """ Execute forward pass through BatchNorm, hidden layers, and output head
                """
                h = self.bn ( x )
                h = self.hidden_net ( h )
                return self.fc_out ( h )

        return WeakMLP

    # =========================================================================
    ## Format optional optimizer-specific diagnostics for the progress bar.
    #  - Subclasses can override this hook to append extra status fields
    #           (for example, actual epochs executed when early stopping triggers)
    #  @param step_info Dictionary of metrics returned by @c _fit_step for the
    #         current boosting iteration.
    #  @return Formatted diagnostic string (empty string if no extra info applies).
    def _format_extra_step_info ( self, step_info ) :
        """ Format optional optimizer-specific progress diagnostics for the step
        """
        return ""

    # =========================================================================
    ## Execute a single weak-learner training step for the active optimizer.
    #  - Called once per boosting step inside @c _fit_iterative_stream
    #    after the combined sample weights @c comb_weights have been
    #           renormalized for the current iteration.
    #           Implementations must
    #           train @c model in-place and return a dictionary containing at
    #           least @c 'initial_loss' and @c 'running_loss'
    #
    #  @param model Freshly initialized @c WeakMLP instance on @c device
    #  @param criterion Element-wise binary cross-entropy loss module
    #         (@c reduction='none')
    #  @param X_comb Combined NumPy feature matrix of shape @c (n_orig + n_targ, n_features)
    #  @param y_comb Combined NumPy binary labels array (0 for original, 1 for target)
    #  @param comb_weights Normalized NumPy sample weights array for @c X_comb
    #  @param X_comb_tensor Combined feature tensor already placed on @c device
    #  @param y_comb_tensor Combined label tensor of shape @c (n_comb, 1) on @c device
    #  @param weights_tensor Normalized sample weight tensor of shape @c (n_comb, 1)
    #         on @c device
    #  @param device Target @c torch.device on which tensors and model reside
    #  @return Dictionary with step loss statistics (@c initial_loss, @c running_loss,
    #          and optional extra diagnostics).
    @abc.abstractmethod
    def _fit_step ( self           ,
                    model          ,
                    criterion      ,
                    X_comb         ,
                    y_comb         ,
                    comb_weights   ,
                    X_comb_tensor  ,
                    y_comb_tensor  ,
                    weights_tensor ,
                    device         ) :
        """ Execute one optimizer step on ``model`` and return loss/step metrics
        """
        raise NotImplementedError

    # =========================================================================
    ## Execute the iterative PyTorch boosting stream on a single CV fold.
    #
    #          Prepares combined training arrays and tensors, iterates for
    #           @c n_estimators boosting steps, delegates parameter optimization
    #           of each @c WeakMLP weak learner to @c _fit_step, applies zero-mean
    #           centering on @c X_orig_tensor, clips log-ratio updates to
    #           @c LOG_STEP_CLIPS, updates @c current_w_orig multiplicatively,
    #           logs periodic progress metrics, and records model weights and
    #           centering shifts
    #
    #  @param X_orig Source feature subset matrix for the current fold
    #  @param w_orig Source sample weights subset for the current fold
    #  @param X_targ Target feature subset matrix for the current fold
    #  @param w_targ Target sample weights subset for the current fold
    #  @return Tuple @c (models, final_ratios) where @c models is a list of
    #          serialized weak-learner dictionaries and @c final_ratios is the
    #          normalized density ratio array for @c X_orig
    def _fit_iterative_stream ( self   ,
                                X_orig ,
                                w_orig ,
                                X_targ ,
                                w_targ ) :
        """ Execute iterative PyTorch boosting stream on fold data

        Runs sequential weak-learner optimization via ``_fit_step``, applies
        zero-mean spatial centering and bounded exponential weight updates,
        and returns saved model dictionaries together with final ratios
        """
        import torch
        import torch.nn as nn

        device = torch.device ( self.params.get ( 'device', 'cuda' if torch.cuda.is_available () else 'cpu' ) )

        n_orig     = len ( X_orig )
        n_targ     = len ( X_targ )
        n_features = X_orig.shape [ 1 ]
        n_comb     = n_orig + n_targ

        X_comb     = numpy.vstack ( [ X_orig, X_targ ] ).astype ( numpy.float32 )
        abs_w_orig = numpy.abs ( w_orig )
        abs_w_targ = numpy.abs ( w_targ )

        # Match total target weight mass to total source weight mass before each step
        sum_w_orig = numpy.sum ( abs_w_orig )
        sum_w_targ = numpy.sum ( abs_w_targ )
        norm_targ  = w_targ * ( sum_w_orig / ( sum_w_targ + 1e-12 ) )

        y_comb = numpy.hstack ( [ numpy.zeros ( n_orig, dtype = numpy.float32 ),
                                  numpy.ones  ( n_targ, dtype = numpy.float32 ) ] )

        abs_norm_targ  = numpy.abs   ( norm_targ ).astype ( numpy.float32, copy = False )
        sum_abs_targ   = numpy.sum   ( abs_norm_targ, dtype = numpy.float64 )
        comb_weights   = numpy.empty ( n_comb, dtype = numpy.float32 )
        current_w_orig = abs_w_orig.copy ()

        n_estimators  = self.params.get ( 'n_estimators' , TORCH_ESTIMATORS    )
        learning_rate = self.params.get ( 'learning_rate', TORCH_LEARNING_RATE )
        hidden_layers = self.params.get ( 'hidden_layers', TORCH_HIDDEN_LAYERS )
        print_freq    = self.params.get ( 'print_freq'   , PRINT_FREQ          )

        models        = []
        X_orig_tensor = torch.tensor ( X_orig.astype ( numpy.float32 ), dtype = torch.float32 ).to ( device )
        X_comb_tensor = torch.tensor ( X_comb, dtype = torch.float32 ).to ( device )
        y_comb_tensor = torch.tensor ( y_comb, dtype = torch.float32 ).unsqueeze ( 1 ).to ( device )

        WeakMLPClass = self.get_MLP_class ()
        criterion    = nn.BCEWithLogitsLoss ( reduction = 'none' )

        neff0 = nEff ( X_orig, current_w_orig )
        description = (
            S.fixed_width ( boosting_symbol      , 2 ) +
            S.fixed_width ( self.symbol          , 2 ) +
            S.fixed_width ( self.optimizer_symbol, 2 ) + ' :'
        )

        if   self.progress : the_bar = progress_bar ( n_estimators, full = True, description = description, silent = not self.progress )
        elif self.verbose  : the_bar = text_bar     ( n_estimators, full = True, description = description, silent =     self.silent   )
        else               : the_bar = progress_bar ( n_estimators, full = True, description = description, silent = True              )

        t_start   = time.time ()
        prev_loss = 0.0

        for step , index , bar in LoopGC ( safe_loop ( the_bar ) ) :

            # Recompute normalized sample weights reflecting updated source weights
            w_sum    = numpy.sum ( current_w_orig, dtype = numpy.float64 ) + sum_abs_targ + 1e-12
            inv_sum  = numpy.float32 ( 1.0 / w_sum )
            comb_weights [ : n_orig ] = current_w_orig * inv_sum
            comb_weights [ n_orig : ] = abs_norm_targ  * inv_sum

            weights_tensor = torch.tensor ( comb_weights, dtype = torch.float32 ).unsqueeze ( 1 ).to ( device )

            model = WeakMLPClass ( n_features, hidden_layers ).to ( device )
            model.train ()

            step_info    = self._fit_step ( model          = model          ,
                                            criterion      = criterion      ,
                                            X_comb         = X_comb         ,
                                            y_comb         = y_comb         ,
                                            comb_weights   = comb_weights   ,
                                            X_comb_tensor  = X_comb_tensor  ,
                                            y_comb_tensor  = y_comb_tensor  ,
                                            weights_tensor = weights_tensor ,
                                            device         = device         )

            initial_loss = float ( step_info [ 'initial_loss' ] )
            running_loss = float ( step_info [ 'running_loss' ] )

            delta_loss = running_loss - prev_loss
            prev_loss  = running_loss

            model.eval ()
            with torch.no_grad () :
                raw_vals = model ( X_orig_tensor ).cpu ().numpy ().ravel ()

            # Zero-mean centering removes global constant shift so the step only
            # shapes relative feature-dependent density variations
            train_shift = float ( numpy.mean ( raw_vals ) )
            raw_vals    = raw_vals - train_shift

            log_r_step      = numpy.clip ( raw_vals, *LOG_STEP_CLIPS )
            current_w_orig *= numpy.exp  ( log_r_step * learning_rate )

            if bar and not self.silent and ( step < 5 or ( step + 1 ) % print_freq == 0 or ( step + 5 ) == n_estimators ) :
                ratios_curr = current_w_orig / ( abs_w_orig + 1e-12 )
                neff        = nEff ( X_orig, current_w_orig )

                wmin    , wmax  = float ( ratios_curr.min () ), float ( ratios_curr.max () )
                wminmax , wexpo = pretty_floats ( wmin, wmax, precision = 2, with_sign = False, strip = True )
                wminmax         = '/'.join ( wminmax )
                if wexpo :
                    wminmax = '(%s)%s%s' % ( wminmax, S.times, format_pow ( wexpo ) )

                ne0 = nice_print ( float ( neff0 ), precision = 2, with_sign = False, strip = True )
                ne  = nice_print ( float ( neff  ), precision = 2, with_sign = False, strip = True )

                dt = ( time.time () - t_start ) / ( step + 1 )
                if dt < 0.1 : dt = '%dms'  % ( math.ceil ( 1000 * dt ) )
                else        : dt = '%.1fs' %                      dt

                init_str  = nice_print ( initial_loss, precision = 4, with_sign = False, strip = True )
                loss_str  = nice_print ( running_loss, precision = 4, with_sign = False, strip = True )
                delta_str = nice_print ( delta_loss  , precision = 4, with_sign = True , strip = True )
                extra_str = self._format_extra_step_info ( step_info )

                message  = f"Step {step+1:4d}/{n_estimators:<4d}"
                message += f" | loss: {init_str} {S.arrow_right} {loss_str}"
                message += f" | {S.delta_symbol}: {delta_str}"
                message += extra_str
                message += f" | wmin/wmax: {wminmax}"
                message += f" | #eff: {ne0} {S.arrow_right} {ne} | Time/step:{dt}"

                bar += message
                del ratios_curr

            models.append ( {
                'state_dict'    : model.state_dict () ,
                'input_dim'     : n_features          ,
                'hidden_layers' : hidden_layers       ,
                'shift'         : train_shift         ,
            } )

            del weights_tensor
            del model
            del raw_vals
            del log_r_step

        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_clip ( final_ratios, w_orig )

        del X_orig_tensor, X_comb_tensor, y_comb_tensor
        del X_comb, y_comb, comb_weights
        del current_w_orig, abs_norm_targ

        return models, final_ratios

    # =========================================================================
    ## Accumulate raw logit predictions across sequential PyTorch models.
    #
    #           Reconstructs each saved @c WeakMLP weak learner, loads its
    #           learned @c state_dict, evaluates unseen features @c X, subtracts
    #           the stored per-step zero-mean @c shift, clips step log-ratios
    #           to @c LOG_STEP_CLIPS, scales by @c learning_rate, accumulates
    #           in log-ratio space, and returns @c exp(log_ratios)
    #
    #  @param models List of model dictionaries produced by @c _fit_iterative_stream
    #  @param X Feature matrix of unseen samples to evaluate
    #  @return Array of predicted density ratios @c r(x)
    def _predict_iterative_stream ( self, models, X ) :
        """ Accumulate predictions across sequential PyTorch weak learners

        Evaluates every saved weak learner on ``X``, applies its stored
        zero-mean shift, clips log-ratio increments, and returns ``exp``
        of the cumulative log-ratio sum
        """
        import torch

        device        = torch.device ( self.params.get ( 'device', 'cuda' if torch.cuda.is_available () else 'cpu' ) )
        learning_rate = self.params.get ( 'learning_rate', TORCH_LEARNING_RATE )
        log_ratios    = numpy.zeros ( len ( X ), dtype = numpy.float32 )
        X_tensor      = torch.tensor ( X.astype ( numpy.float32 ), dtype = torch.float32 ).to ( device )

        WeakMLPClass = self.get_MLP_class ()

        for m_dict in models :
            h_layers = m_dict.get ( 'hidden_layers', TORCH_HIDDEN_LAYERS )

            model = WeakMLPClass ( m_dict [ 'input_dim' ], h_layers ).to ( device )
            model.load_state_dict ( m_dict [ 'state_dict' ] )
            model.eval ()

            with torch.no_grad () :
                raw_vals = model ( X_tensor ).cpu ().numpy ().ravel ()

            train_shift = m_dict.get ( 'shift', 0.0 )
            raw_vals    = raw_vals - train_shift

            log_r_step  = numpy.clip ( raw_vals, *LOG_STEP_CLIPS )
            log_r_step *= learning_rate
            log_ratios += log_r_step

            del model
            del log_r_step
            del raw_vals

        rexp = numpy.exp ( log_ratios )

        del log_ratios
        del X_tensor
        return rexp

# =============================================================================
## @class PyTorchAdamIterativeReweighter
#  Iterative weight boosting density ratio reweighter using the @c Adam
#  stochastic optimizer.
#
#  Specializes @c PyTorchBaseIterativeReweighter to train each @c WeakMLP
#           weak learner via mini-batch stochastic gradient descent with
#           @c torch.optim.Adam. Supports power-of-two / automatic /
#           full-batch sizing and early stopping governed by @c min_epochs,
#           @c min_delta, and @c patience over @c epochs_per_step epochs
class PyTorchAdamIterativeReweighter ( PyTorchBaseIterativeReweighter ) :
    """ Iterative weight boosting density reweighter backed by PyTorch ``Adam``

    Trains each weak MLP learner with mini-batch ``Adam`` optimization and
    early stopping per boosting iteration
    """
    # =========================================================================
    ## Initialize the PyTorch Adam iterative density ratio reweighter.
    #
    #           Sets default mini-batch, learning-rate, weight-decay, and
    #           early-stopping options for @c Adam, normalizes integer batch
    #           sizes down to the largest power of two not exceeding the input
    #           limit, and delegates initialization to @c PyTorchIterativeReweighter
    #
    #  @param original Source feature matrix (N-dimensional array or dataframe)
    #  @param target Target feature matrix (N-dimensional array or dataframe)
    #  @param original_weight Optional observation weights for the source dataset
    #  @param target_weight Optional observation weights for the target dataset
    #  @param kwargs Additional keyword arguments including @c batch_size,
    #         @c optimizer_learning_rate, @c epochs_per_step, @c min_epochs,
    #         @c min_delta, and @c patience
    def __init__ ( self            ,
                   original        ,
                   target          ,
                   original_weight = None ,
                   target_weight   = None ,
                   **kwargs        ) :
        """ Initialize ``PyTorchAdamIterativeReweighter`` with Adam and early-stopping defaults
        """
        config = {
            'batch_size'              : TORCH_BATCH_SIZE           ,
            'optimizer_learning_rate' : 0.01                       ,
            'epochs_per_step'         : TORCH_ADAM_EPOCHS_PER_STEP ,
            'min_epochs'              : 5                          ,
            'min_delta'               : 1e-4                       ,
            'patience'                : 5                          ,
        }

        config.update ( kwargs )

        n_rows     = num_samples  ( target ) + num_samples ( original )
        n_features = num_features ( target )

        batch_size = config.get ( 'batch_size', 'auto' )
        
        if isinstance ( batch_size, integer_types ) and 2 <= batch_size <= min ( n_rows, TORCH_MAX_BATCH_SIZE ) :
            if batch_size & ( batch_size - 1 ) : batch_size = 2 ** math.floor ( math.log2 ( batch_size ) )
        elif isinstance ( batch_size, string_types ) and batch_size.strip ().lower () == 'full' : batch_size = n_rows
        else : batch_size = optimal_batch_size ( n_rows, n_features, TORCH_MAX_BATCH_SIZE )

        config [ 'batch_size' ] = batch_size

        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **config )

    # =========================================================================
    ## Get the symbol/tag representing the Adam optimizer backend.
    #  @return Adam visual symbol used in progress bar descriptions
    @property
    def optimizer_symbol ( self ) :
        """ Return visual symbol representing the ``Adam`` optimizer
        """
        return TORCH_ADAM_SYMBOL

    # =========================================================================
    ## Format epoch count when early stopping terminates before @c epochs_per_step.
    #  @param step_info Dictionary of metrics returned by @c _fit_step.
    #  @return Formatted string with actual epoch count, or empty string.
    def _format_extra_step_info ( self, step_info ) :
        """ Append actual epoch count when early stopping stops early
        """
        epochs_per_step = self.params.get ( 'epochs_per_step', TORCH_ADAM_EPOCHS_PER_STEP )
        nepochs         = step_info.get ( 'nepochs', epochs_per_step )
        if nepochs != epochs_per_step : return f" | #epochs/step:{nepochs}"
        return ""

    # =========================================================================
    ## Train one @c WeakMLP weak learner using mini-batch @c Adam.
    #
    #          Builds a @c TensorDataset and shuffled @c DataLoader over the
    #           combined dataset, normalizes weights within each mini-batch,
    #           updates parameters with @c torch.optim.Adam, and stops early once
    #           mean epoch loss fails to improve by at least @c min_delta for
    #           @c patience consecutive epochs after @c min_epochs
    #
    #  @param model @c WeakMLP instance to train in-place
    #  @param criterion Element-wise @c BCEWithLogitsLoss (@c reduction='none')
    #  @param X_comb Combined NumPy feature matrix (unused directly; tensors are used)
    #  @param y_comb Combined NumPy label array (unused directly; tensors are used)
    #  @param comb_weights Normalized NumPy weight array (unused directly; tensors are used)
    #  @param X_comb_tensor Combined feature tensor on @c device
    #  @param y_comb_tensor Combined label tensor on @c device
    #  @param weights_tensor Normalized weight tensor on @c device
    #  @param device Target @c torch.device
    #  @return Dictionary with @c initial_loss, @c running_loss, and @c nepochs
    def _fit_step ( self           ,
                    model          ,
                    criterion      ,
                    X_comb         ,
                    y_comb         ,
                    comb_weights   ,
                    X_comb_tensor  ,
                    y_comb_tensor  ,
                    weights_tensor ,
                    device         ) :
        """ Train a single weak learner step with mini-batch ``Adam`` and early stopping
        """
        import torch.optim as optim
        from torch.utils.data import TensorDataset, DataLoader

        batch_size      = self.params.get ( 'batch_size'             , TORCH_BATCH_SIZE      )
        lr              = self.params.get ( 'optimizer_learning_rate', 0.01                  )
        weight_decay    = self.params.get ( 'weight_decay'           , TORCH_WEIGHT_DECAY    )
        epochs_per_step = self.params.get ( 'epochs_per_step'        , TORCH_EPOCHS_PER_STEP )
        min_epochs      = self.params.get ( 'min_epochs'             , 5                     )
        min_delta       = self.params.get ( 'min_delta'              , 1e-4                  )
        patience        = self.params.get ( 'patience'               , 5                     )

        optimizer = optim.Adam ( model.parameters (),
                               lr           = lr           ,
                               weight_decay = weight_decay )

        dataset = TensorDataset ( X_comb_tensor, y_comb_tensor, weights_tensor )
        loader  = DataLoader    ( dataset, batch_size = batch_size, shuffle = True )

        initial_loss     = None
        running_loss     = 0.0
        best_loss        = pos_infinity
        patience_counter = 0
        nepochs          = 0

        for epoch in safe_loop ( epochs_per_step ) :
            
            nepochs   += 1
            epoch_loss = 0.0
            batches    = 0

            for bx, by, bw in loader :
                optimizer.zero_grad ()
                preds   = model ( bx )
                bw_norm = bw / ( bw.sum () + 1e-12 )
                loss    = ( criterion ( preds, by ) * bw_norm ).sum ()

                current_loss = float ( loss.item () )
                loss.backward ()
                optimizer.step ()

                if initial_loss is None :
                    initial_loss = current_loss

                epoch_loss += current_loss
                batches    += 1

            running_loss = float ( epoch_loss / max ( batches, 1 ) )

            if min_epochs <= epoch :
                if running_loss + min_delta < best_loss :
                    best_loss        = running_loss
                    patience_counter = 0
                else :
                    patience_counter += 1

                if patience <= patience_counter :
                    break

        return { 'initial_loss' : 0.0 if initial_loss is None else initial_loss ,
                 'running_loss' : running_loss                                  ,
                 'nepochs'      : nepochs                                       }


# =============================================================================
## @class PyTorchLBFGSIterativeReweighter
#  Iterative weight boosting density ratio reweighter using the
#  quasi-Newton @c L-BFGS optimizer.
#
#           Specializes @c PyTorchIterativeReweighter to optimize each @c WeakMLP
#           weak learner using @c torch.optim.LBFGS with strong Wolfe line
#           search. When the combined dataset exceeds @c lbfgs_chunk,
#           draws a random subset without replacement for the step so memory
#           footprint and quasi-Newton step cost remain bounded
class PyTorchLBFGSIterativeReweighter ( PyTorchBaseIterativeReweighter ) :
    """ Iterative weight boosting density reweighter backed by PyTorch ``L-BFGS`

    Trains each weak MLP learner using ``L-BFGS`` optimization with optional
    random chunk subsampling for large datasets
    """
    # =========================================================================
    ## Initialize the PyTorch L-BFGS iterative density ratio reweighter.
    #
    #           Sets default @c L-BFGS step budget (@c epochs_per_step),
    #           step size (@c optimizer_learning_rate), L2 regularization
    #           (@c weight_decay), and chunk size (@c lbfgs_chunk), then
    #           delegates initialization to @c PyTorchIterativeReweighter.
    #
    #  @param original Source feature matrix (N-dimensional array or dataframe)
    #  @param target Target feature matrix (N-dimensional array or dataframe)
    #  @param original_weight Optional observation weights for the source dataset
    #  @param target_weight Optional observation weights for the target dataset
    #  @param kwargs Additional keyword arguments including @c lbfgs_chunk,
    #         @c epochs_per_step, @c optimizer_learning_rate, and @c weight_decay
    def __init__ ( self            ,
                   original        ,
                   target          ,
                   original_weight = None ,
                   target_weight   = None ,
                   **kwargs        ) :
        """ Initialize ``PyTorchLBFGSIterativeReweighter`` with L-BFGS defaults
        """
        config = {
            'lbfgs_chunk'             : TORCH_LBFGS_CHUNK           ,
            'optimizer_learning_rate' : 0.01                        ,
            'epochs_per_step'         : TORCH_LBFGS_EPOCHS_PER_STEP ,
        }

        config.update ( kwargs )

        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **config )

    # =========================================================================
    ## Get the symbol/tag representing the L-BFGS optimizer backend.
    #  @return L-BFGS visual symbol used in progress bar descriptions
    @property
    def optimizer_symbol ( self ) :
        """ Return visual symbol representing the ``L-BFGS`` optimizer
        """
        return TORCH_LBFGS_SYMBOL

    # =========================================================================
    ## Train one @c WeakMLP weak learner using @c L-BFGS on a chunk.
    #
    #           Selects either the full dataset or a random sample of size
    #           @c lbfgs_chunk, evaluates initial weighted loss plus optional
    #           L2 penalty, runs @c optimizer.step(closure) with
    #           @c line_search_fn='strong_wolfe', and evaluates final loss
    #
    #  @param model @c WeakMLP instance to optimize in-place
    #  @param criterion Element-wise @c BCEWithLogitsLoss (@c reduction='none')
    #  @param X_comb Combined NumPy feature matrix of shape @c (n_comb, n_features)
    #  @param y_comb Combined NumPy binary label array of shape @c (n_comb,)
    #  @param comb_weights Normalized NumPy sample weights of shape @c (n_comb,)
    #  @param X_comb_tensor Combined feature tensor on @c device
    #  @param y_comb_tensor Combined label tensor on @c device
    #  @param weights_tensor Normalized weight tensor on @c device
    #  @param device Target @c torch.device
    #  @return Dictionary containing @c initial_loss and @c running_loss
    def _fit_step ( self           ,
                    model          ,
                    criterion      ,
                    X_comb         ,
                    y_comb         ,
                    comb_weights   ,
                    X_comb_tensor  ,
                    y_comb_tensor  ,
                    weights_tensor ,
                    device         ) :
        """ Train a single weak learner step with ``L-BFGS`` and strong Wolfe line search
        """
        import torch
        import torch.optim as optim

        lbfgs_chunk  = self.params.get ( 'lbfgs_chunk'            , TORCH_LBFGS_CHUNK           )
        max_iter     = self.params.get ( 'epochs_per_step'        , TORCH_LBFGS_EPOCHS_PER_STEP )
        lr           = self.params.get ( 'optimizer_learning_rate', 0.01                        )
        weight_decay = self.params.get ( 'weight_decay'           , TORCH_WEIGHT_DECAY          )

        n_comb = len ( X_comb )

        if lbfgs_chunk < n_comb :
            chunk_idx = numpy.random.choice ( n_comb, size = lbfgs_chunk, replace = False )
            X_ch      = X_comb [ chunk_idx ]
            y_ch      = y_comb [ chunk_idx ]
            w_ch      = comb_weights [ chunk_idx ]
            w_ch      = w_ch / ( numpy.sum ( w_ch ) + 1e-12 )

            X_ch_tensor = torch.tensor ( X_ch, dtype = torch.float32 ).to ( device )
            y_ch_tensor = torch.tensor ( y_ch, dtype = torch.float32 ).unsqueeze ( 1 ).to ( device )
            w_ch_tensor = torch.tensor ( w_ch, dtype = torch.float32 ).unsqueeze ( 1 ).to ( device )
        else :
            X_ch_tensor = X_comb_tensor
            y_ch_tensor = y_comb_tensor
            w_ch_tensor = weights_tensor

        ## Evaluate weighted BCE loss plus optional explicit L2 penalty
        def _compute_loss () :
            """ Evaluate weighted BCE loss plus optional explicit L2 penalty
            """
            preds = model ( X_ch_tensor )
            loss  = ( criterion ( preds, y_ch_tensor ) * w_ch_tensor ).sum ()
            if 0 < weight_decay :
                l2_reg = sum ( torch.sum ( p ** 2 ) for p in model.parameters () )
                loss   = loss + 0.5 * weight_decay * l2_reg
            return loss

        with torch.no_grad () :
            initial_loss = float ( _compute_loss ().item () )

        optimizer = optim.LBFGS ( model.parameters ()             ,
                                  lr             = lr             ,
                                  max_iter       = max_iter       ,
                                  line_search_fn = 'strong_wolfe' )

        ## Zero gradients, compute loss, and backpropagate for L-BFGS
        def closure () :
            """ Zero gradients, compute loss, and backpropagate for L-BFGS
            """
            optimizer.zero_grad ()
            loss = _compute_loss ()
            loss.backward ()
            return loss

        optimizer.step ( closure )

        with torch.no_grad () :
            running_loss = float ( _compute_loss ().item () )

        if lbfgs_chunk < n_comb :
            del X_ch_tensor, y_ch_tensor, w_ch_tensor

        return { 'initial_loss' : initial_loss ,
                 'running_loss' : running_loss }

    
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
