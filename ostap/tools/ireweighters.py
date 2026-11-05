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
XGB_MAX_DEPTH         =    7 # Deep tree capacity for high dimensions / large statistics
XGB_LEARNING_RATE     = 0.02 # Small step for stable density ratio convergence

CATB_ESTIMATORS       =  300 # High capacity: deep ensemble with low learning rate
CATB_MAX_DEPTH        =    6 # Deep tree capacity for high dimensions / large statistics
CATB_LEARNING_RATE    = 0.10 # Small step for stable density ratio convergence

TORCH_ESTIMATORS      =  500 # Number of boosting steps (epochs/models)
TORCH_LEARNING_RATE   = 0.02 # Step size for weight boosting
TORCH_HIDDEN_DIM      =  128 # Hidden layer size for small MLP weak learners
TORCH_EPOCHS_PER_STEP =  100 # Training epochs per boosting step


THRESHOLD_NSAMPLES    = 50000
THRESHOLD_NFEATURES   = 10

PRINT_FREQ            = 10 
# =============================================================================
## @class LightGBMIterativeReweighter
#  @brief Iterative density ratio reweighter based on LightGBM gradient boosting trees.
#  @details Implements sequential gradient boosting for density ratio estimation $r(x) = p(x)/(1-p(x))$
#           using LightGBM datasets with dynamically updated sample weights on each boosting iteration.
class LightGBMIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using LightGBM backend
    """
    # =========================================================================
    ## @brief Initialize the LightGBM iterative density reweighter instance.
    #  @param self Instance of the class.
    #  @param original Source feature matrix (N-dimensional array or dataframe).
    #  @param target Target feature matrix (N-dimensional array or dataframe).
    #  @param original_weight Optional observation weights for the source dataset (supports signed weights).
    #  @param target_weight Optional observation weights for the target dataset (supports signed weights).
    #  @param kwargs Additional keyword arguments for LightGBM configuration and parent initializers.
    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        """ Initialize LightGBM iterative density reweighter with custom defaults
        """
        config = {
            'n_estimators'            :  LGBM_ESTIMATORS    ,
            'learning_rate'           :  LGBM_LEARNING_RATE ,
            'max_depth'               :  LGBM_MAX_DEPTH     ,
            ## 
            'max_bins'                :  1023               ,
            'min_data_in_leaf'        :  10                 ,
            'min_sum_hessian_in_leaf' :  1e-4               , 
            'subsample'               :  1.0                ,
            'bagging_freq'            :  0                  , 
            'lambda_l1'               :  0.0                ,
            'feature_fraction'        :  1.0                , 
            'lambda_l2'               :  0.0                ,                        
            'verbosity'               : -1                  ,
            'n_jobs'                  :  1                  ,
            'n_splits'                :  5                  , 
        }
        
        config.update ( kwargs )
        if 'silent' in config and config.get ( 'silent' ) : config [ 'verbosity' ] = -1        
        if 'num_boost_round'  in config : config [ 'n_estimators'            ] = config.pop ( 'num_boost_round'   , LGBM_ESTIMATORS ) 
        if 'min_child_weight' in config : config [ 'min_sum_hessian_in_leaf' ] = config.pop ( 'min_child_weight'  , 1e-2 )         
        if 'bagging_fraction' in config : config [ 'subsample'               ] = config.pop ( 'bagging_fraction'  , 1.0 )         
        if 'colsample_bytree' in config : config [ 'feature_fraction'        ] = config.pop ( 'colsample_by_tree' , 1.0 )         
        
        ## Ignore 'subsample' argument for small smapels or splitting/folsig regime 
        n_rows = num_samples ( original ) + num_samples ( target )
        if n_rows < THRESHOLD_NSAMPLES or 1 < config.get ( 'n_splits' , 5 )     : config [ 'subsample'        ] = 1.0    
        if 1 != config.get ( 'subsample' , 1 ) and not 'bagging_freq' in config : config [ 'bagging_freq'     ] = 1            
        if num_features    ( original ) < THRESHOLD_NFEATURES                   : config [ 'feature_fraction' ] = 1.0

        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **config )

    # =========================================================================
    ## @property symbol
    #  @brief Short symbol or tag representing the LightGBM reweighter method.
    #  @return String identifier/symbol used for logging and visual presentation.
    @classprop
    def symbol ( self ) :
        """ Short symbol/tag for reweighter. """
        return U.symbol_LGBM
    
    @property
    def method ( self ) :
        """ Return method identifier name. """
        return method_LGBM 
    
    # =========================================================================
    ## @brief Execute the iterative LightGBM boosting stream on a specific cross-validation fold.
    #  @details Sequentially fits LightGBM weak learners, updates the weights of original events 
    #           via exponential learning-rate scaling steps, monitors effective sample size convergence, 
    #           and renders live status updates via progress bars and safe loops.
    #  @param self Instance of the class.
    #  @param X_orig Original feature subset matrix for the current fold.
    #  @param w_orig Original sample weights subset for the current fold.
    #  @param X_targ Target feature subset matrix for the current fold.
    #  @param w_targ Target sample weights subset for the current fold.
    #  @return A tuple `(models, final_ratios)` containing the trained boosting chain and computed ratios.
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

        neff0 = nEff ( X_orig, current_w_orig )

        description = boosting_symbol + ' ' + S.fixed_width ( self.symbol , 3 ) + ':'        
        if   self.progress : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = not self.progress ) 
        elif self.verbose  : the_bar = text_bar     ( n_estimators , full = True , description = description , silent =     self.silent   )
        else               : the_bar = progress_bar ( n_estimators , full = True , description = descritpion , silent = True  ) 

        t_start = time.time () 
        eps     = 1e-4 

        for it , _ ,  bar in safe_loop ( the_bar ) :
            
            comb_weights  = numpy.hstack ( [ current_w_orig, numpy.abs ( norm_targ ) ] )
            comb_weights /= ( numpy.sum ( comb_weights ) + 1e-12 )

            trn_data = LightGBM.Dataset ( X_comb, label = y_comb, weight = comb_weights, free_raw_data = True )            
            model    = LightGBM.train   ( params          = params   ,
                                          train_set       = trn_data ,
                                          num_boost_round = 1        )

            p_orig          = model.predict ( X_orig )            
            p_orig          = numpy.clip ( p_orig, eps, 1.0 - eps )
            r_step          = p_orig / ( 1.0 - p_orig )
            current_w_orig *= numpy.power ( numpy.clip ( r_step, 1e-2, 1e+2 ), learning_rate )

            if bar and not self.silent and ( it < 5  or ( it + 1 ) % PRINT_FREQ == 0 or ( it + 5 ) == n_estimators ) :
                ratios_curr = current_w_orig / ( abs_w_orig + 1e-12 )
                neff        = nEff ( X_orig, current_w_orig )
                wmin        = nice_print ( float ( ratios_curr.min () ) , precision = 2, with_sign = False, strip = True )
                wmax        = nice_print ( float ( ratios_curr.max () ) , precision = 2, with_sign = False, strip = True )
                ne0         = nice_print ( float ( neff0              ) , precision = 2, with_sign = False, strip = True )
                ne          = nice_print ( float ( neff               ) , precision = 2, with_sign = False, strip = True )

                dt       = ( time.time () - t_start ) / ( it + 1 )
                if dt < 1 : dt = '%dms'  % ( math.ceil ( 1000 * dt ) )
                else      : dt = '%.1fs' %                      dt                 
                message  = f"Iteration {it+1:4d}/{n_estimators:-4d} | min/max: {wmin}/{wmax} | #eff: {ne0} {S.arrow_right} {ne} | Time/tree:{dt}" 
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
    ## @brief Accumulate predictions and compute log-ratios across sequential LightGBM models.
    #  @details Evaluates each trained LightGBM booster on feature matrix X, clips predicted probabilities 
    #           to prevent numerical divergence, converts them to log-ratios, and accumulates them 
    #           scaled by the boosting learning rate.
    #  @param self Instance of the class.
    #  @param models List of trained LightGBM booster models comprising the sequential chain.
    #  @param X Feature matrix of unseen samples to evaluate.
    #  @return An array of predicted density ratios.
    def _predict_iterative_stream ( self, models, X ) :
        """ Accumulate predictions across sequential LightGBM models. """
        
        learning_rate = self.params.get ( 'learning_rate', LGBM_LEARNING_RATE )
        log_ratios    = numpy.zeros ( len ( X ), dtype = numpy.float32 )

        for model in models :
            p           = model.predict ( X )
            p           = numpy.clip ( p, 1e-4, 1.0 - 1e-4 )
            r_step      = p / ( 1.0 - p )
            log_ratios += numpy.log ( numpy.clip ( r_step, 1e-2, 1e+2 ) ) * learning_rate

        return numpy.exp ( log_ratios )


# =============================================================================
## @class XGBoostIterativeReweighter
#  Iterative weight boosting density reweighter using XGBoost.
class XGBoostIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using XGBoost backend
    """

    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        """ Initialize XGBoost iterative density reweighter
        """
        config = {
            'n_estimators'      :  XGB_ESTIMATORS    ,
            'learning_rate'     :  XGB_LEARNING_RATE ,
            'max_depth'         :  XGB_MAX_DEPTH     ,
            ## 
            'tree_method'       :  'hist'            , 
            'max_bin'           :  1023              , 
            'grow_policy'       : 'depthwise'        , 
            ##
            'max_leaves'        :  63                , 
            'min_child_weight'  :  1e-3              , 
            'max_delta_step'    :  9.0               , 
            'gamma'             :  0.0               , 
            'subsample'         :  1.0               , 
            'colsamples_bytree' :  1.0               , 
            'reg_alpha'         :  0.0               , 
            'reg_lambda'        :  0.0               , 
            'verbosity'         :  0                 , 
            'n_jobs'            :  1                 ,
            'n_splits'          :  5                 ,
        }

        config.update ( kwargs )
        if 'silent' in config and config.get ( 'silent' ) : config [ 'verbosity' ] = 0        
        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round', XGB_ESTIMATORS ) 
        
        ## Ignore 'subsample' argument for small samples or splitting/folds regime 
        n_rows = num_samples ( original ) + num_samples ( target )
        if n_rows < THRESHOLD_NSAMPLES or 1 < config.get ( 'n_splits' , 5 ) : config [ 'subsample'         ] = 1.0

        if num_features      ( original ) < THRESHOLD_NFEATURES             : config [ 'colsamples_bytree' ] = 1.0

        ## Use exact method if statistics is not large enought  
        if n_rows < 10000 : config [ 'tree_method' ]  = 'exact'

        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **config )

    # =========================================================================
    ## Short symbol/tag for reweighter 
    @classprop
    def symbol ( self ) :
        """ `symbol` : short symbol/tag for reweighter """
        return U.symbol_XGB

    @property
    def method ( self ) :
        """ Return method identifier name. """
        return method_XGB

    ## @brief Execute iterative XGBoost boosting stream on fold data using margins.
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
        max_depth     = self.params.get ( 'max_depth'    , XGB_MAX_DEPTH     )
        max_leaves    = self.params.get ( 'max_leaves'   , 63                )
    
        models = []

        params = {
            ## 
            'objective'        : 'binary:logistic' , 
            'eval_metric'      : 'logloss'         ,
            ## 
            'tree_method'      : self.params.get ( 'tree_method'      , 'hist'      ) , 
            'max_bin'          : self.params.get ( 'max_bin'          , 1023        ) ,
            ## 
            'min_child_weight' : self.params.get ( 'min_child_weight' , 1e-3        ) ,
            'max_delta_step'   : self.params.get ( 'max_delta_step'   , 9.0         ) , 
            'subsample'        : self.params.get ( 'subsample'        , 1.0         ) ,
            'colsample_bytree' : self.params.get ( 'colsample_bytree' , 1.0         ) ,
            ##
            'gamma'            : self.params.get ( 'gamma'            , 0.0         ) ,             
            'reg_alpha'        : self.params.get ( 'reg_alpha'        , 0.0         ) , 
            'reg_lambda'       : self.params.get ( 'reg_lambda'       , 0.0         ) ,
            ## 
            'max_depth'        : max_depth ,
            'learning_rate'    : 1.0       ,
            ##
            'grow_policy'      : self.params.get ( 'grow_policy', 'depthwise' ) , 
            'max_leaves'       : max_leaves , 
            'verbosity'        : 0                 , 
            'n_jobs'           : self.params.get ( 'n_jobs'   , 1 ) , 
        }
        
        neff0 = nEff ( X_orig, current_w_orig )
        
        description = boosting_symbol + ' ' + S.fixed_width ( self.symbol , 3 ) + ':'        
        if   self.progress : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = not self.progress ) 
        elif self.verbose  : the_bar = text_bar     ( n_estimators , full = True , description = description , silent =     self.silent   )
        else               : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = True  ) 

        t_start = time.time () 
        eps     = 1e-4

        for it , _ ,  bar in safe_loop ( the_bar ) :
                        
            comb_weights  = numpy.hstack ( [ current_w_orig, numpy.abs ( norm_targ ) ] )
            comb_weights /= ( numpy.sum ( comb_weights ) + 1e-12 )
            
            dmat_comb = XGBoost.DMatrix ( X_comb, label = y_comb, weight = comb_weights )
            
            model = XGBoost.train ( params          = params,
                                    dtrain          = dmat_comb,
                                    num_boost_round = 1 )
            
            dmat_orig = XGBoost.DMatrix ( X_orig )
            p_orig    = model.predict ( dmat_orig )
            p_orig    = numpy.clip ( p_orig, eps, 1.0 - eps )
            r_step    = p_orig / ( 1.0 - p_orig )
            
            current_w_orig *= numpy.power ( numpy.clip ( r_step, 1e-2, 1e+2 ), learning_rate )
        
            if bar and not self.silent and ( it < 5  or ( it + 1 ) % PRINT_FREQ == 0 or ( it + 5 ) == n_estimators ) :
                ratios_curr = current_w_orig / ( abs_w_orig + 1e-12 )
                neff        = nEff ( X_orig, current_w_orig )
                wmin        = nice_print ( float ( ratios_curr.min () ) , precision = 2, with_sign = False, strip = True )
                wmax        = nice_print ( float ( ratios_curr.max () ) , precision = 2, with_sign = False, strip = True )
                ne0         = nice_print ( float ( neff0              ) , precision = 2, with_sign = False, strip = True )
                ne          = nice_print ( float ( neff               ) , precision = 2, with_sign = False, strip = True )

                dt       = ( time.time () - t_start ) / ( it + 1 )
                if dt < 1 : dt = '%dms'  % ( math.ceil ( 1000 * dt ) )
                else      : dt = '%.1fs' %                      dt                 
                message  = f"Iteration {it+1:4d}/{n_estimators:-4d} | min/max: {wmin}/{wmax} | #eff: {ne0} {S.arrow_right} {ne} | Time/tree:{dt}" 
                bar     += message
                
            models.append ( model )
            
        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_clip ( final_ratios, w_orig )

        return models, final_ratios

    # =========================================================================
    ## @brief Accumulate predictions across sequential XGBoost models.
    def _predict_iterative_stream ( self, models, X ) :
        """ Accumulate predictions across sequential XGBoost models. """
        import xgboost as XGBoost
        dmat          = XGBoost.DMatrix ( X )
        learning_rate = self.params.get ( 'learning_rate', XGB_LEARNING_RATE )
        log_ratios    = numpy.zeros ( len ( X ), dtype = numpy.float32 )

        for model in models :
            p           = model.predict ( dmat )
            p           = numpy.clip ( p, 1e-4, 1.0 - 1e-4 )
            r_step      = p / ( 1.0 - p )
            log_ratios += numpy.log ( numpy.clip ( r_step, 1e-2, 1e+2 ) ) * learning_rate

        return numpy.exp ( log_ratios )

    
# =============================================================================
## @class CatBoostIterativeReweighter
#  @brief Iterative density ratio reweighter based on CatBoost gradient boosting trees.
#  @details Implements sequential gradient boosting for density ratio estimation $r(x) = p(x)/(1-p(x))$
#           using custom objective functions adapted for weighted classification tasks.
#           Leverages symmetric trees (oblivious trees) with custom Newton leaf estimation.
class CatBoostIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using CatBoost.
    """

    # =========================================================================
    ## @class DensityRatioObj
    #  @brief Custom density ratio estimation (DRE) objective function for CatBoost.
    #  @details Computes the first-order (gradient) and second-order (Hessian) derivatives
    #           of the custom loss function tailored for direct density ratio estimation 
    #           with weighted original and target samples.
    class DensityRatioObj ( object ) :
        """ Custom density ratio estimation objective function for CatBoost. """

        ## @brief Calculate first and second derivatives for a range of approxes, targets, and weights.
        #  @param self Instance of the objective class.
        #  @param approxes Raw formula values (margins) predicted by the boosting model.
        #  @param targets True binary classification labels (0 for original stream, 1 for target stream).
        #  @param weights Observation weights associated with the samples.
        #  @return A list of tuples containing `(der1, der2)` for each sample.
        def calc_ders_range ( self, approxes, targets, weights ) :
            
            margin = numpy.clip ( approxes, -15.0, 15.0 )
            exp_m  = numpy.exp  ( margin )
            p      = exp_m / ( 1.0 + exp_m )
            
            # Correct gradients: p - y (where y=1 for target, y=0 for original)
            grad   = numpy.where ( targets == 1.0, p - 1.0, p )
            
            # Correct Hessian: p * (1 - p) for both classes (ensures target points contribute to curvature)
            hess   = p * ( 1.0 - p ) + 1e-5
            
            if weights is not None and 0 < len ( weights ) :
                weights = numpy.maximum ( weights, 1e-5 )
                grad    = grad * weights
                hess    = hess * weights
                
            return list ( zip ( grad, hess ) )
        
    # =========================================================================
    ## @brief Initialize the CatBoost iterative density reweighter instance.
    #  @param self Instance of the class.
    #  @param original Source feature matrix (N-dimensional array or dataframe).
    #  @param target Target feature matrix (N-dimensional array or dataframe).
    #  @param original_weight Optional observation weights for the source dataset (supports signed weights).
    #  @param target_weight Optional observation weights for the target dataset (supports signed weights).
    #  @param kwargs Additional keyword arguments for CatBoost configuration and parent initializers.
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
            'grow_policy'                :  'Lossguide'        ,
            'max_leaves'                 :  31                 ,
            'l2_leaf_reg'                :  0.001              ,
            'subsample'                  :  1.0                ,
            'rsm'                        :  1.0                , 
            'random_strength'            :  0.0                ,
            'leaf_estimation_iterations' :  1                  ,
            'min_data_in_leaf'           :  1                  ,
            'verbose'                    :  False              ,
            'thread_count'               :  1                  ,
            'n_splits'                   :  5                  , 
        }
        
        config.update ( kwargs )
        if 'silent'           in config and config.get ( 'silent' ) : config [ 'verbose' ] = 0        
        if 'num_boost_round'  in config : config [ 'n_estimators'   ] = config.pop ( 'num_boost_round', CATB_ESTIMATORS ) 
        if 'colsample_bytree' in config : config [ 'rsm'            ] = config.pop ( 'colsample_by_tree' , 1.0 )         
        
        ## Ignore 'subsample' argument for small smapels or splitting/folsig regime 
        n_rows = num_samples ( original ) + num_samples ( target )
        if n_rows < THRESHOLD_NSAMPLES or 1 < config.get ( 'n_splits' , 5 ) : config [ 'subsample' ] = 1.0
        if num_features    ( original ) < THRESHOLD_NFEATURES               : config [ 'rsm'       ] = 1.0

        
        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **config )

    # =========================================================================
    ## @property symbol
    #  @brief Short symbol or tag representing the CatBoost reweighter method.
    #  @return String identifier/symbol used for logging and visual presentation.
    @classprop
    def symbol ( self ) :
        """ Short symbol/tag for reweighter. """
        return U.symbol_CATB

    @property
    def method ( self ) :
        """ Return method identifier name. """
        return method_CATB

    # =========================================================================
    ## @brief Execute the iterative CatBoost boosting stream on a specific cross-validation fold.
    #  @details Sequentially fits CatBoost weak learners using the custom DRE objective,
    #           updates the weights of original events via exponential learning-rate scaling steps,
    #           monitors effective sample size convergence, and renders live status updates via progress bars.
    #  @param self Instance of the class.
    #  @param X_orig Original feature subset matrix for the current fold.
    #  @param w_orig Original sample weights subset for the current fold.
    #  @param X_targ Target feature subset matrix for the current fold.
    #  @param w_targ Target sample weights subset for the current fold.
    #  @return A tuple `(models, final_ratios)` containing the trained boosting chain and computed ratios.
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
        max_depth     = self.params.get ( 'max_depth'    , CATB_MAX_DEPTH     )
        max_leaves    = self.params.get ( 'max_leaves'   , 31                 )

        models    = []
        pool_orig = CatBoost.Pool ( data = X_orig )

        params = {
            'loss_function'             : self.DensityRatioObj ()                         ,
            'eval_metric'               : self.params.get ( 'eval_metric', 'Logloss' )    ,
            'iterations'                : 1                                               ,
            'learning_rate'             : 1.0                                             ,
            'depth'                     : max_depth                                       ,
            'grow_policy'               : self.params.get ( 'grow_policy', 'Lossguide' )  ,
            'max_leaves'                : max_leaves                                      ,
            'l2_leaf_reg'               : self.params.get ( 'l2_leaf_reg', 0.001       )  , 
            'subsample'                 : self.params.get ( 'subsample'  , 1.0         )  ,
            'rsm'                       : self.params.get ( 'rsm'        , 1.0         )  ,            
            'random_strength'           : 0.0                                             ,
            'leaf_estimation_iterations': 1                                               ,
            'min_data_in_leaf'          : self.params.get ( 'min_data_in_leaf', 1 )       , 
            'verbose'                   : False                                           , 
            'thread_count'              : self.params.get ( 'n_jobs', 1 )                 ,
            'leaf_estimation_method'    : 'Newton'                                        ,
        }

        neff0 = nEff ( X_orig, current_w_orig )

        description = boosting_symbol + ' ' + S.fixed_width ( self.symbol , 3 ) + ':'        
        if   self.progress : the_bar = progress_bar ( n_estimators , full = True , description = description , silent = not self.progress ) 
        elif self.verbose  : the_bar = text_bar     ( n_estimators , full = True , description = description , silent =     self.silent   )
        else               : the_bar = progress_bar ( n_estimators , full = True , description = descritpion , silent = True  ) 

        t_start = time.time () 

        for it , _ ,  bar in safe_loop ( the_bar ) :
            
            comb_weights  = numpy.hstack ( [ current_w_orig, numpy.abs ( norm_targ ) ] )
            comb_weights /= ( numpy.sum ( comb_weights ) + 1e-12 )

            trn_pool = CatBoost.Pool ( data = X_comb, label = y_comb, weight = comb_weights )
            
            model    = CatBoost.train ( pool    = trn_pool ,
                                        params  = params   ,
                                        verbose = False    )

            margin = model.predict ( pool_orig, prediction_type = 'RawFormulaVal' )
            if 1 < margin.ndim : margin = margin [ :, 0 ]
            
            margin          = numpy.clip ( margin, -15.0, 15.0 )
            step_ratios     = numpy.exp  ( margin )
            current_w_orig *= numpy.power ( numpy.clip ( step_ratios, 1e-6, 1e+6 ), learning_rate )

            if bar and not self.silent and ( it < 5  or ( it + 1 ) % PRINT_FREQ == 0 or ( it + 5 ) == n_estimators ) :
                ratios_curr = current_w_orig / ( abs_w_orig + 1e-12 )
                neff        = nEff ( X_orig, current_w_orig )
                wmin        = nice_print ( float ( ratios_curr.min () ) , precision = 2, with_sign = False, strip = True )
                wmax        = nice_print ( float ( ratios_curr.max () ) , precision = 2, with_sign = False, strip = True )
                ne0         = nice_print ( float ( neff0              ) , precision = 2, with_sign = False, strip = True )
                ne          = nice_print ( float ( neff               ) , precision = 2, with_sign = False, strip = True )

                dt       = ( time.time () - t_start ) / ( it + 1 )
                if dt < 1 : dt = '%dms'  % ( math.ceil ( 1000 * dt ) )
                else      : dt = '%.1fs' %                      dt                 
                message  = f"Iteration {it+1:4d}/{n_estimators:-4d} | min/max: {wmin}/{wmax} | #eff: {ne0} {S.arrow_right} {ne} | Time/tree:{dt}" 
                bar     += message

            models.append ( model )

        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_clip ( final_ratios, w_orig )

        return models, final_ratios

    # =========================================================================
    ## @brief Predict density ratios directly using the sequential CatBoost model chain.
    #  @details Accumulates raw formula values (margins) from all trained trees in log-space,
    #           scales them by the learning rate, applies clipping boundaries for numerical stability,
    #           and converts them back via exponential transformation.
    #  @param self Instance of the class.
    #  @param models List of trained CatBoost models comprising the boosting chain.
    #  @param X Feature matrix of unseen samples to evaluate.
    #  @return An array of predicted density ratios normalized to a local mean of 1.0.
    def _predict_iterative_stream ( self, models, X ) :
        """ Predict ratios directly using iterative CatBoost stream. """
        import catboost as CatBoost
        
        learning_rate = self.params.get ( 'learning_rate', CATB_LEARNING_RATE )
        pool          = CatBoost.Pool ( data = X )
        log_ratios    = numpy.zeros ( len ( X ), dtype = numpy.float32 )

        for model in models :
            margin = model.predict ( pool, prediction_type = 'RawFormulaVal' )
            if 1 < margin.ndim : margin = margin [ :, 0 ]
            
            margin      = numpy.clip ( margin, -15.0, 15.0 )
            log_ratios += margin * learning_rate
            
        ratios = numpy.exp ( numpy.clip ( log_ratios, -15.0, 15.0 ) )
        
        mean_r = numpy.mean ( ratios )
        if 0 < mean_r : ratios /= mean_r
            
        return ratios
    
# =============================================================================
## @class PyTorchIterativeReweighter
#  Iterative weight boosting density reweighter using PyTorch (MLP weak learners).
class PyTorchIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using PyTorch. """

    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        """ Initialize PyTorch iterative density reweighter. """
        config = {
            'n_estimators'    :  TORCH_ESTIMATORS       ,
            'learning_rate'   :  TORCH_LEARNING_RATE    ,
            'hidden_dim'      :  TORCH_HIDDEN_DIM       ,
            'epochs_per_step' :  TORCH_EPOCHS_PER_STEP  ,
            'batch_size'      :  'auto'                 ,
            'weight_decay'    :  1e-4                   ,
            'verbosity'       :  0                      ,
            'n_jobs'          :  1                      , 
        }
        
        config.update ( kwargs )
        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round', TORCH_ESTIMATORS ) 

        n_rows     = num_samples ( target ) + num_samples ( original )        
        batch_size = config.get ( 'batch_size', None ) 
        if   isinstance ( batch_size, string_types  ) and batch_size.lower () == 'auto' : batch_size = None
        elif isinstance ( batch_size, integer_types ) and not 2 <= batch_size < n_rows   : batch_size = None        
        if batch_size is None : 
            batch_size = 2 ** math.floor ( max ( 2, math.log2 ( n_rows / 25 ) ) ) 
            config [ 'batch_size' ] = batch_size 

        super ().__init__ ( original        = original        ,
                            target          = target          ,
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **config )

    # =========================================================================
    ## short symbol/tag for reweighter 
    @classprop
    def symbol ( self ) :
        """`symbol` : short symbol/tag for reweighter 
        """
        return U.symbol_TORCH

    @property
    def method ( self ) :
        """ Return method identifier name. """
        return method_TORCH

    ## @brief Execute iterative PyTorch boosting stream on fold data.
    def _fit_iterative_stream ( self   ,
                                X_orig ,
                                w_orig ,
                                X_targ ,
                                w_targ ) :
        """ Execute iterative PyTorch boosting stream on fold data. """
        import torch
        import torch.nn as nn
        import torch.optim as optim

        device     = torch.device ( 'cpu' )
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
        
        n_estimators    = self.params.get ( 'n_estimators'   , TORCH_ESTIMATORS      )
        learning_rate   = self.params.get ( 'learning_rate'  , TORCH_LEARNING_RATE   )
        hidden_dim      = self.params.get ( 'hidden_dim'     , TORCH_HIDDEN_DIM      )
        epochs_per_step = self.params.get ( 'epochs_per_step', TORCH_EPOCHS_PER_STEP )
        batch_size      = self.params.get ( 'batch_size'     , 256                   )
        weight_decay    = self.params.get ( 'weight_decay'   , 1e-4                  )

        class WeakMLP ( nn.Module ) :
            def __init__ ( self, input_dim, hidden_dim ) :
                super ().__init__ ()
                self.net = nn.Sequential (
                    nn.Linear ( input_dim, hidden_dim ),
                    nn.Tanh   (),
                    nn.Linear ( hidden_dim, hidden_dim ),
                    nn.Tanh   (),
                    nn.Linear ( hidden_dim, 1 )
                )
            def forward ( self, x ) :
                return self.net ( x )
            
        models        = []
        X_orig_tensor = torch.tensor ( X_orig.astype ( numpy.float32 ), dtype = torch.float32 ).to ( device )
        X_comb_tensor = torch.tensor ( X_comb, dtype = torch.float32 ).to ( device )
        y_comb_tensor = torch.tensor ( y_comb, dtype = torch.float32 ).unsqueeze ( 1 ).to ( device )

        eps = 1e-4 

        for it in safe_progress_bar ( n_estimators, description = 'Boosting:', silent = not self.progress ) :
            
            comb_weights   = numpy.hstack ( [ current_w_orig, numpy.abs ( norm_targ ) ] )
            comb_weights  /= ( numpy.sum ( comb_weights ) + 1e-12 )
            weights_tensor = torch.tensor ( comb_weights, dtype = torch.float32 ).unsqueeze ( 1 ).to ( device )

            model     = WeakMLP ( n_features, hidden_dim ).to ( device )
            optimizer = optim.Adam ( model.parameters (), lr = 0.01, weight_decay = weight_decay )
            criterion = nn.BCEWithLogitsLoss ( reduction = 'none' )

            model.train ()
            dataset   = torch.utils.data.TensorDataset ( X_comb_tensor, y_comb_tensor, weights_tensor )
            loader    = torch.utils.data.DataLoader ( dataset, batch_size = batch_size, shuffle = True )

            for epoch in range ( epochs_per_step ) :
                for bx, by, bw in loader :
                    optimizer.zero_grad ()
                    preds = model ( bx )
                    loss  = torch.mean ( criterion ( preds, by ) * bw )
                    loss.backward ()
                    optimizer.step ()

            model.eval ()
            with torch.no_grad () :
                p_orig = model ( X_orig_tensor ).cpu ().numpy ().ravel ()

            p_orig          = numpy.clip ( p_orig, eps, 1.0 - eps )
            r_step          = p_orig / ( 1.0 - p_orig )
            current_w_orig *= numpy.power ( numpy.clip ( r_step, 1e-2, 1e+2 ), learning_rate )

            models.append ( { 'state_dict' : model.state_dict (), 'input_dim' : n_features, 'hidden_dim' : hidden_dim } )
            
        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self._normalize_and_clip ( final_ratios, w_orig )

        return models, final_ratios

    ## @brief Accumulate predictions across sequential PyTorch models.
    def _predict_iterative_stream ( self, models, X ) :
        """ Accumulate predictions across sequential PyTorch models. """
        import torch
        import torch.nn as nn

        device        = torch.device ( 'cpu' )
        learning_rate = self.params.get ( 'learning_rate', TORCH_LEARNING_RATE )
        log_ratios    = numpy.zeros ( len ( X ), dtype = numpy.float32 )
        X_tensor      = torch.tensor ( X.astype ( numpy.float32 ), dtype = torch.float32 ).to ( device )

        class WeakMLP ( nn.Module ) :
            def __init__ ( self, input_dim, hidden_dim ) :
                super ().__init__ ()
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
            model = WeakMLP ( m_dict [ 'input_dim' ], m_dict [ 'hidden_dim' ] ).to ( device )
            model.load_state_dict ( m_dict [ 'state_dict' ] )
            model.eval ()

            with torch.no_grad () :
                logits = model ( X_tensor )
                p      = torch.sigmoid ( logits ).cpu ().numpy ().ravel ()

            p           = numpy.clip ( p, 1e-4, 1.0 - 1e-4 )
            r_step      = p / ( 1.0 - p )
            log_ratios += numpy.log ( numpy.clip ( r_step, 1e-2, 1e+2 ) ) * learning_rate

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
