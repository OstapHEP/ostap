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

XGB_ESTIMATORS        = 1000 # High capacity: deep ensemble with low learning rate
XGB_MAX_DEPTH         =    8 # Deep tree capacity for high dimensions / large statistics
XGB_LEARNING_RATE     = 0.05 # Small step for stable density ratio convergence

CATB_ESTIMATORS       = 1000 # High capacity: deep ensemble with low learning rate
CATB_MAX_DEPTH        =    8 # Deep tree capacity for high dimensions / large statistics
CATB_LEARNING_RATE    = 0.05 # Small step for stable density ratio convergence

TORCH_ESTIMATORS      =  500 # Number of boosting steps (epochs/models)
TORCH_LEARNING_RATE   = 0.02 # Step size for weight boosting
TORCH_HIDDEN_DIM      =  128 # Hidden layer size for small MLP weak learners
TORCH_EPOCHS_PER_STEP =  100 # Training epochs per boosting step

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
            'feature_fraction'        :  0.9  , 
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
            'feature_fraction'        : self.params.get ( 'feature_fraction'         , 0.9  ) , 
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
        final_ratios = self.__normalize_and_iter_clip ( final_ratios , w_orig , norm_targ )

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

    # =========================================================================
    def __normalize_and_iter_clip ( self, ratios, w_orig, w_targ ) :
        clipped  = numpy.clip ( ratios, 0.0 , self.config [ 'clip_threshold' ] )
        orig_sum = numpy.sum ( numpy.abs ( w_orig ) * clipped , dtype = numpy.float64 )
        targ_sum = numpy.sum ( numpy.abs ( w_targ )          , dtype = numpy.float64 )
        norm     = targ_sum / ( orig_sum + 1e-12 )
        return clipped * norm


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
            'tree_method'             :  'hist' , 
            'max_bin'                 :  1023   ,
            'min_child_weight'        :  1e-6   , 
            ##
            'gamma'                   :  0     , 
            'subsample'               :  0.95  , 
            'colsample_bytree'        :  1.0   , 
            'reg_alpha'               :  0     ,
            'reg_lambda'              :  0     ,                        
            ## 
            'verbosity'               :  0     ,
            'n_jobs'                  :  1     , ## we are working with a single tree 
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
        
        sum_w_orig = numpy.sum ( abs_w_orig )
        sum_w_targ = numpy.sum ( abs_w_targ )
        
        norm_targ = w_targ * ( sum_w_orig / ( sum_w_targ + 1e-12 ) )

        y_comb = numpy.hstack ( [ numpy.zeros ( n_orig , dtype = numpy.float32 ) ,
                                  numpy.ones  ( n_targ , dtype = numpy.float32 ) ] )

        current_w_orig = abs_w_orig.copy ()
        
        n_estimators  = self.params.get ( 'n_estimators'  , XGB_ESTIMATORS    )
        learning_rate = self.params.get ( 'learning_rate' , XGB_LEARNING_RATE )
        max_depth     = self.params.get ( 'max_depth'     , XGB_MAX_DEPTH     )
    
        models = []

        params = {
            'objective'               : self.params.get ( 'objective'          , 'binary:logistic' ) ,     
            'eval_metric'             : self.params.get ( 'eval_metric'        , 'logloss'         ) , 
            'tree_method'             : self.params.get ( 'tree_method'        , 'hist' ) , 
            'max_bin'                 : self.params.get ( 'max_bin'            , 1023   ) ,
            'min_child_weight'        : self.params.get ( 'min_child_weight'   , 1e-4   ) , 
            'subsample'               : self.params.get ( 'subsample'          , 0.95   ) , 
            'colsample_bytree'        : self.params.get ( 'colsample_bytree'   , 0.9    ) ,        
            'gamma'                   : self.params.get ( 'gamma'              , 0      ) , 
            'reg_alpha'               : self.params.get ( 'reg_alpha'          , 0      ) , 
            'reg_lambda'              : self.params.get ( 'reg_lambda'         , 0      ) , 
            'max_depth'               : max_depth ,
            'learning_rate'           : 1.0       , 
            'verbosity'               : self.params.get ( 'verbosity' , 0 ) , 
            'n_jobs'                  : self.params.get ( 'n_jobs'    , 1 ) , 
        }

        dmat_orig = XGBoost.DMatrix ( X_orig )

        for it in progress_bar ( n_estimators , description = 'Boosting:' , silent = not self.progress ) :
            
            comb_weights  = numpy.hstack ( [ current_w_orig , numpy.abs ( norm_targ ) ] )
            comb_weights /= ( numpy.sum ( comb_weights ) + 1e-12 )

            trn_data = XGBoost.DMatrix ( X_comb , label = y_comb , weight = comb_weights )
            model    = XGBoost.train   ( params          = params   ,
                                         dtrain          = trn_data ,
                                         num_boost_round = 1        )

            ## # БЕРЕМ СЫРЫЕ МАРГИНЫ ВМЕСТО ВЕРОЯТНОСТЕЙ:
            ## margins        = model.predict ( dmat_orig , output_margin = True )
            ## log_step       = numpy.clip    ( margins , -10.0 , 10.0 )
            ## current_w_orig *= numpy.exp    ( log_step * learning_rate )


            # Для binary:logitraw метод predict сразу возвращает сырые логиты!
            log_step       = model.predict ( dmat_orig )
            log_step       = numpy.clip    ( log_step , -10.0 , 10.0 )
            current_w_orig *= numpy.exp    ( log_step * learning_rate )

            
            # Очистка памяти (механизм Феникса)
            model_bytes = model.save_raw ()
            del model
            del trn_data
            model = XGBoost.Booster ()
            model.load_model ( model_bytes )
            del model_bytes
            
            models.append ( model )
            
        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self.__normalize_and_iter_clip ( final_ratios , w_orig , norm_targ )
        
        return models, final_ratios

    # =========================================================================
    def _predict_iterative_stream ( self, models, X ) :
        import xgboost as XGBoost
        learning_rate = self.params.get ( 'learning_rate' , XGB_LEARNING_RATE )
        log_ratios    = numpy.zeros ( len ( X ) , dtype = numpy.float32 )
        dmat          = XGBoost.DMatrix ( X )
        
        for model in models :
            margins     = model.predict ( dmat , output_margin = True )
            log_ratios += numpy.clip ( margins , -10.0 , 10.0 ) * learning_rate
            
        return numpy.exp ( log_ratios )
    
    # =========================================================================
    def __normalize_and_iter_clip ( self, ratios, w_orig, w_targ ) :
        clipped  = numpy.clip ( ratios, 0.0 , self.config [ 'clip_threshold' ] )
        orig_sum = numpy.sum ( numpy.abs ( w_orig ) * clipped , dtype = numpy.float64 )
        targ_sum = numpy.sum ( numpy.abs ( w_targ )          , dtype = numpy.float64 )
        norm     = targ_sum / ( orig_sum + 1e-12 )
        return clipped * norm

# =============================================================================
## @class CatBoostIterativeReweighter
#  Iterative weight boosting density reweighter using CatBoost.
class CatBoostIterativeReweighter ( IterativeReweighter ) :
    """ Iterative weight boosting density reweighter using CatBoost.
    """
    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        
        config = {
            ## 
            'n_estimators'            :  CATB_ESTIMATORS    ,
            'learning_rate'           :  CATB_LEARNING_RATE ,
            'max_depth'               :  CATB_MAX_DEPTH     ,
            ## 
            'l2_leaf_reg'             :  0.0   , 
            'subsample'               :  0.9   , 
            'random_strength'         :  0.0   ,
            ## 
            'verbose'                 :  False ,
            'thread_count'            :  1     , ## we are working with a single tree 
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
        l2_leaf_reg   = self.params.get ( 'l2_leaf_reg'   , 0.0                )
        subsample     = self.params.get ( 'subsample'     , 0.9                )

        models = []

        params = {
            'loss_function'   : 'Logloss'   ,
            'iterations'      : 1           ,
            'learning_rate'   : 1.0         , # ATTENTION!!! Managed externally via weights
            'depth'           : max_depth   ,
            'l2_leaf_reg'     : l2_leaf_reg ,
            'subsample'       : subsample   ,
            'random_strength' : self.params.get ( 'random_strength' , 0.0   ) ,
            'verbose'         : self.params.get ( 'verbose'         , False ) , 
            'thread_count'    : self.params.get ( 'thread_count'    , 1     ) ,
        }

        eps = 1e-4 
        pool_orig = CatBoost.Pool ( data = X_orig )

        for it in progress_bar ( n_estimators , description = 'Boosting:' , silent = not self.progress ) :
            
            comb_weights  = numpy.hstack ( [ current_w_orig , numpy.abs ( norm_targ ) ] )
            comb_weights /= ( numpy.sum ( comb_weights ) + 1e-12 )

            trn_pool = CatBoost.Pool ( data = X_comb , label = y_comb , weight = comb_weights )
            
            # В CatBoost для дообучения по одному шагу используется метод train с parameters
            model = CatBoost.train ( pool             = trn_pool ,
                                     params           = params   ,
                                     verbose          = False    )


            ## p_orig = model.predict ( pool_orig , prediction_type = 'Probability' )
            ## if 1 < p_orig.ndim : p_orig = p_orig [ :, 1 ]
            ## p_orig          = numpy.clip    ( p_orig , eps , 1.0 - eps )
            ## r_step          = p_orig / ( 1.0 - p_orig )
            ## current_w_orig *= numpy.power ( numpy.clip ( r_step , 1e-3 , 1e+3 ) , learning_rate )

            margins         = model.predict ( pool_orig , prediction_type = 'RawFormulaVal' )
            log_step        = numpy.clip    ( margins , -10.0 , 10.0 )
            current_w_orig *= numpy.exp     ( log_step * learning_rate )
            

            # =================================================================
            if True : # =======================================================
                # =============================================================
                # Трюк «Феникса» для CatBoost: сохранение в байты/диск и очистка
                import tempfile, os
                with tempfile.NamedTemporaryFile(suffix='.cbm', delete=False) as tmp:
                    tmp_name = tmp.name
                
                model.save_model(tmp_name)
                del model
                del trn_pool
                
                model = CatBoost.CatBoostClassifier()
                model.load_model(tmp_name)
                try:
                    os.unlink(tmp_name)
                except OSError:
                    pass
            
            models.append ( model )
            
        final_ratios = current_w_orig / ( abs_w_orig + 1e-12 )
        final_ratios = self.__normalize_and_iter_clip ( final_ratios , w_orig , norm_targ )

        return models, final_ratios

    # =========================================================================
    def _predict_iterative_stream ( self, models, X ) :
        
        import catboost as CatBoost
        learning_rate = self.params.get ( 'learning_rate' , CATB_LEARNING_RATE )
        log_ratios    = numpy.zeros ( len ( X ) , dtype = numpy.float32 )
        pool          = CatBoost.Pool ( data = X )

        ## for model in models :
        ##    p           = model.predict ( pool , prediction_type = 'Probability' )
        ##    if 1 < p.ndim : p = p [ :, 1 ]        
        ##    p           = numpy.clip ( p , 1e-4 , 1.0 - 1e-4 )
        ##    r_step      = p / ( 1.0 - p )
        ##    log_ratios += numpy.log ( numpy.clip ( r_step , 1e-2 , 1e+2 ) ) * learning_rate
        ## 
        ## return numpy.exp ( log_ratios )
    
        for model in models :
            margins     = model.predict ( pool , prediction_type = 'RawFormulaVal' )
            log_ratios += numpy.clip    ( margins , -10.0 , 10.0 ) * learning_rate

        return numpy.exp ( log_ratios )
    
    # =========================================================================
    def __normalize_and_iter_clip ( self, ratios, w_orig, w_targ ) :
        clipped  = numpy.clip ( ratios, 0.0 , self.config [ 'clip_threshold' ] )
        orig_sum = numpy.sum ( numpy.abs ( w_orig ) * clipped , dtype = numpy.float64 )
        targ_sum = numpy.sum ( numpy.abs ( w_targ )           , dtype = numpy.float64 )
        norm     = targ_sum / ( orig_sum + 1e-12 )
        return clipped * norm


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
        final_ratios = self.__normalize_and_iter_clip ( final_ratios , w_orig , norm_targ )

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

    # =========================================================================
    def __normalize_and_iter_clip ( self, ratios, w_orig, w_targ ) :
        clipped  = numpy.clip ( ratios, 0.0 , self.config [ 'clip_threshold' ] )
        orig_sum = numpy.sum ( numpy.abs ( w_orig ) * clipped , dtype = numpy.float64 )
        targ_sum = numpy.sum ( numpy.abs ( w_targ )          , dtype = numpy.float64 )
        norm     = targ_sum / ( orig_sum + 1e-12 )
        return clipped * norm

    
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
