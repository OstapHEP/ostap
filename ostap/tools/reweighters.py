#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with various "reweighters"
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-18
# =============================================================================
""" Module with various `reweighters'
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    'Reweighter_LGBM' , ## LightGBM-based reweighter
    'Reweighter_XGB'  , ## XGBoost-based reweighter
    'Reweighter_CATB' , ## CatBoost-based reweighter
) 
# =============================================================================
from   ostap.math.math_base  import weight_trivial
from   ostap.utils.core      import typename
from   ostap.utils.basic     import numcpu, num_jobs
from   typing                import Any, Dict, Tuple, Optional
import numpy, abc, warnings  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighters' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## @class BsseDensityReweighter
#  Abstract base class for immutable, adaptive density-ratio reweighting.
#  Implements 1-, 2-, and 4-stream signed-measure density estimations.
# 
#  The instance executes training immediately upon instantiation (__init__)
#  and seals the resulting weights for the original sample as read-only attributes
#  (if store_original_weights=True).
class BaseDensityReweighter(abc.ABC):
    """ Abstract base class for immutable, adaptive density-ratio reweighting.
    Implements 1-, 2-, and 4-stream signed-measure density estimations.
    
    The instance executes training immediately upon instantiation (__init__)
    and seals the resulting weights for the original sample as read-only attributes
    (if store_original_weights=True).
    """

    def __init__( self, * , 
                  original               : numpy.ndarray,
                  target                 : numpy.ndarray,
                  original_weight        : Optional [ numpy.ndarray ] = None ,
                  target_weight          : Optional [ numpy.ndarray ] = None ,
                  clip_threshold         : float = 10.0 ,
                  n_splits               : int   = 5    ,
                  store_original_weights : bool = True ) :
        
        self.__clip_threshold = float(clip_threshold)
        self.__n_splits = int(n_splits)
        self.__store_original_weights = bool(store_original_weights)

        # Internal storage for fitted models and stream parameters
        self.__fitted_models = {}
        self.__priors = {}
        self.__target_weights_info = {}
        self.__norm_factor = numpy.float32(1.0)
        self.__mode = None
        
        # Execute training immediately and compute original sample weights (float32)
        original_ratios, original_reweighted_weights = self.__fit_and_compute (
            original.astype        ( numpy.float32 , copy = False ) ,
            original_weight.astype ( numpy.float32 , copy = False ) if original_weight is not None else None,
            target.astype          ( numpy.float32 , copy = False ) ,
            target_weight.astype   ( numpy.float32 , copy = False ) if target_weight is not None else None,
        )
        
        # Optionally store original sample reweighting results
        if self.__store_original_weights:
            self.__original_ratios             = original_ratios
            self.__original_reweighted_weights = original_reweighted_weights
        else:
            self.__original_ratios             = None
            self.__original_reweighted_weights = None

    # ==================================================================
    # Public Read-Only Properties
    # ==================================================================
    @property
    def original_ratios(self) -> Optional [ numpy.ndarray ]:
        """ Computed density ratio factors r(x) for the original sample.
        Returns None if store_original_weights=False during initialization.
        """
        return self.__original_ratios

    @property
    def original_reweighted_weights(self) -> Optional [ numpy.ndarray ] : 
        """ Final reweighted event weights (w * r) for the original sample.
        Returns None if store_original_weights=False during initialization.
        """
        return self.__original_reweighted_weights

    @property
    def store_original_weights(self) -> bool:
        """ Flag indicating whether original sample weights are stored in memory."""
        return self.__store_original_weights

    @property
    def mode(self) -> str:
        """ Active stream decomposition scheme ('1-stream', '2-stream_orig', '2-stream_targ', or '4-stream')."""
        return self.__mode

    # ==================================================================
    # Abstract Methods (Must be implemented by child framework wrappers)
    # ==================================================================
    @abc.abstractmethod
    def _train_single_model( self ,
                             X_train : numpy.ndarray,
                             y_train : numpy.ndarray,
                             w_train : Optional [ numpy.ndarray ] ,
                             X_val   : numpy.ndarray,
                             y_val   : numpy.ndarray,
                             w_val   : Optional [ numpy.ndarray ] ,
                            ) -> Tuple [ Any , numpy.ndarray ] :
        """ Train a single base classifier fold and return (model, val_predictions)."""
        raise NotImplementedError

    ## Predict probability p(y=1|x) using a single trained model.
    @abc.abstractmethod
    def _predict_single_model ( self        , 
                                model : Any  , 
                                X     : numpy.ndarray ) -> numpy.ndarray:
        """ Predict probability p(y=1|x) using a single trained model.

        Parameters
        ----------
        model : Any
            Trained booster or classifier instance (e.g., lgb.Booster, xgb.Booster, CatBoost).
        X : numpy.ndarray
            Feature matrix to make predictions on.

        Returns
        -------
        numpy.ndarray
            Predicted probabilities p(y=1|x) with shape (n_samples,), dtype float32/float64.
        """
        raise NotImplementedError
    
    # ==================================================================
    # Internal Fitting Pipeline
    # ==================================================================
    def __fit_and_compute ( self, 
                            original        : numpy.ndarray, 
                            original_weight : Optional[ numpy.ndarray ] , 
                            target          : numpy.ndarray, 
                            target_weight   : Optional[ numpy.ndarray ] ) :

        # ==============================================================
        # FAST PATH: Standard positive unit weights for both samples
        # ==============================================================
        if original_weight is None and target_weight is None:
            
            self.__mode  = '1-stream'
            ratios_valid = self.__fit_eval_stream ( original , None , target , None, 'base' )
            
            # if all weights are equalto 1: sum ( w ) / sum ( w * r )
            clipped = numpy.clip ( ratios_valid ,
                                   numpy.float32 ( 0.0 ) ,
                                   numpy.float32 ( self.__clip_threshold ) )
            self.__norm_factor = numpy.float32 ( len ( original ) / ( numpy.sum ( clipped, dtype =numpy.float64) + 1e-12 ) )
            
            original_ratios    = clipped * self.__norm_factor
            return original_ratios, original_ratios.copy()

        # ==============================================================
        # SLOW PATH: Explicit handling of provided weight arrays
        # ==============================================================
        w_orig = original_weight if original_weight is not None else numpy.ones ( len ( original ) , dtype = numpy.float32 )
        w_targ = target_weight   if target_weight   is not None else numpy.ones ( len ( target   ) , dtype = numpy.float32 )
        
        # ==============================================================
        # 1. Filter out zero-weight events
        nz_orig_mask = w_orig != 0
        nz_targ_mask = w_targ != 0

        X_orig_valid , w_orig_valid = original [ nz_orig_mask ] , w_orig [ nz_orig_mask ]
        X_targ_valid , w_targ_valid = target   [ nz_targ_mask ] , w_targ [ nz_targ_mask ]

        # ==============================================================
        # 2. Strict sign masks
        mask_orig_pos , mask_orig_neg = w_orig_valid > 0 , w_orig_valid < 0
        mask_targ_pos , mask_targ_neg = w_targ_valid > 0 , w_targ_valid < 0

        has_neg_orig = numpy.any ( mask_orig_neg )
        has_neg_targ = numpy.any ( mask_targ_neg )
        
        # ==============================================================
        # SCENARIO 1: Standard positive weights only
        # ==============================================================
        if not has_neg_orig and not has_neg_targ:
            self.__mode  = '1-stream'
            ratios_valid = self.__fit_eval_stream ( X_orig_valid ,
                                                    w_orig_valid ,
                                                    X_targ_valid ,
                                                    w_targ_valid , 'base' )

        # ==============================================================
        # SCENARIO 2: Negative weights ONLY in Original sample
        # ==============================================================
        elif has_neg_orig and not has_neg_targ:
            self.__mode  = '2-stream_orig'
            ratios_valid = numpy.zeros ( len ( X_orig_valid ), dtype = numpy.float32 )
            ratios_valid [ mask_orig_pos ] = self.__fit_eval_stream ( X_orig_valid [ mask_orig_pos ] ,
                                                                      w_orig_valid [ mask_orig_pos ] ,
                                                                      X_targ_valid ,
                                                                      w_targ_valid , 'orig_pos' ) 
            ratios_valid [ mask_orig_neg ] = self.__fit_eval_stream ( X_orig_valid [ mask_orig_neg ] ,
                                                                      w_orig_valid [ mask_orig_neg ] ,
                                                                      X_targ_valid ,
                                                                      w_targ_valid , 'orig_neg' )
        # ==============================================================
        # SCENARIO 3: Negative weights ONLY in Target sample
        # ==============================================================
        elif not has_neg_orig and has_neg_targ:
            self.__mode = '2-stream_targ'
            X_targ_pos , w_targ_pos = X_targ_valid [ mask_targ_pos ] , w_targ_valid [ mask_targ_pos ]
            X_targ_neg , w_targ_neg = X_targ_valid [ mask_targ_neg ] , w_targ_valid [ mask_targ_neg ]

            w_pos_sum = numpy.float32 ( numpy.sum ( w_targ_pos , dtype = numpy.float64 ) )
            w_neg_sum = numpy.float32 ( numpy.abs ( numpy.sum ( w_targ_neg , dtype = numpy.float64 ) ) )
            self.__target_weights_info = { 'W_pos': w_pos_sum, 'W_neg': w_neg_sum }
            w_total   = w_pos_sum + w_neg_sum

            r_pos_pos = self.__fit_eval_stream ( X_orig_valid , w_orig_valid , X_targ_pos , w_targ_pos , 'pos_pos' )
            r_pos_neg = self.__fit_eval_stream ( X_orig_valid , w_orig_valid , X_targ_neg , w_targ_neg , 'pos_neg' )

            ratios_valid = (r_pos_pos * w_pos_sum - r_pos_neg * w_neg_sum) / (w_total + numpy.float32(1e-12))
    
        # ==============================================================
        # SCENARIO 4: Fully symmetric 4-stream (Negative weights in both)
        # ==============================================================
        else:
            
            self.__mode = '4-stream'
            X_orig_pos , w_orig_pos = X_orig_valid [ mask_orig_pos ] , w_orig_valid [ mask_orig_pos ]
            X_orig_neg , w_orig_neg = X_orig_valid [ mask_orig_neg ] , w_orig_valid [ mask_orig_neg ]

            X_targ_pos , w_targ_pos = X_targ_valid [ mask_targ_pos ] , w_targ_valid [ mask_targ_pos ]
            X_targ_neg , w_targ_neg = X_targ_valid [ mask_targ_neg ] , w_targ_valid [ mask_targ_neg ]

            w_pos_sum = numpy.float32 ( numpy.sum ( w_targ_pos , dtype = numpy.float64 ) )
            w_neg_sum = numpy.float32 ( numpy.abs ( numpy.sum ( w_targ_neg, dtype = numpy.float64 ) ) )
            self.__target_weights_info = {'W_pos': w_pos_sum, 'W_neg': w_neg_sum}
            w_total   = w_pos_sum + w_neg_sum

            r_pos_pos = self.__fit_eval_stream ( X_orig_pos , w_orig_pos , X_targ_pos , w_targ_pos , 'pos_pos')
            r_pos_neg = self.__fit_eval_stream ( X_orig_pos , w_orig_pos , X_targ_neg , w_targ_neg , 'pos_neg')
            r_neg_pos = self.__fit_eval_stream ( X_orig_neg , w_orig_neg , X_targ_pos , w_targ_pos , 'neg_pos')
            r_neg_neg = self.__fit_eval_stream ( X_orig_neg , w_orig_neg , X_targ_neg , w_targ_neg , 'neg_neg') 

            ratios_valid = numpy.zeros ( len ( X_orig_valid ) , dtype = numpy.float32 )
            denom       = w_total + numpy.float32 ( 1e-12 )
            ratios_valid [ mask_orig_pos ] = ( r_pos_pos * w_pos_sum - r_pos_neg * w_neg_sum ) / denom
            ratios_valid [ mask_orig_neg ] = ( r_neg_pos * w_pos_sum - r_neg_neg * w_neg_sum ) / denom


        # ==============================================================
        # 3. Clip, normalize, and construct final arrays
        # ==============================================================
        ratios_valid_norm = self.__normalize_and_clip(ratios_valid, w_orig_valid)

        original_ratios   = numpy.zeros ( len ( w_orig ) , dtype = numpy.float32 )
        original_ratios [ nz_orig_mask ] = ratios_valid_norm
        original_reweighted_weights      = w_orig * original_ratios

        return original_ratios, original_reweighted_weights

    
    def __fit_eval_stream ( self, 
                            X_orig_sub: numpy.ndarray        , 
                            w_orig_sub: Optional[ numpy.ndarray ] , 
                            X_targ_sub: Optional[ numpy.ndarray ] , 
                            w_targ_sub: Optional[ numpy.ndarray ] , 
                            stream_key: str                  ) -> numpy.ndarray:
        
        X_comb = numpy.vstack ( [ X_orig_sub , X_targ_sub ] )
        y_comb = numpy.hstack ( [ numpy.zeros ( len ( X_orig_sub ), dtype = numpy.float32 ), 
                                  numpy.ones  ( len ( X_targ_sub ), dtype = numpy.float32 ) ] )

        # add weights to the model only t aat least one if not trivial
        if w_orig_sub is None and w_targ_sub is None:
            w_comb   = None
            sum_orig = len ( X_orig_sub )
            sum_targ = len ( X_targ_sub )
        else:
            w_orig_abs = numpy.abs ( w_orig_sub ) if w_orig_sub is not None else numpy.ones ( len ( X_orig_sub ) , dtype = numpy.float32 )
            w_targ_abs = numpy.abs ( w_targ_sub ) if w_targ_sub is not None else numpy.ones ( len ( X_targ_sub ) , dtype = numpy.float32 )
            w_comb     = numpy.hstack ( [w_orig_abs , w_targ_abs ] )
            sum_orig   = numpy.sum ( w_orig_abs, dtype = numpy.float64 )
            sum_targ   = numpy.sum ( w_targ_abs, dtype = numpy.float64 )

        oof_preds = numpy.zeros(len(X_comb), dtype=numpy.float32)

        from sklearn.model_selection import StratifiedKFold        
        skf       = StratifiedKFold ( n_splits     = self.__n_splits     ,
                                      shuffle      = True                ,
                                      random_state = self.__random_state )
        
        stream_models = []
        for train_idx, val_idx in skf.split ( X_comb , y_comb ):
            X_tr , y_tr = X_comb [ train_idx ] , y_comb [ train_idx ]
            X_va , y_va = X_comb [ val_idx   ] , y_comb [ val_idx   ]
            
            w_tr = w_comb [ train_idx ] if w_comb is not None else None
            w_va = w_comb [ val_idx   ] if w_comb is not None else None

            model , val_preds = self._train_single_model ( X_tr , y_tr , w_tr, X_va, y_va, w_va )
            oof_preds [ val_idx ] = val_preds.astype ( numpy.float32 , copy = False )
            stream_models.append ( model )

        self.__fitted_models [ stream_key ] = stream_models

        prior  = numpy.float32 ( sum_orig / ( sum_targ + 1e-12 ) )
        self.__priors [ stream_key ] = prior

        p_orig         = oof_preds[:len(X_orig_sub)]
        p_orig_clipped = numpy.clip ( p_orig , numpy.float32 ( 1e-7 ) , numpy.float32 ( 1.0 - 1e-7 ) )
        return ( p_orig_clipped / ( numpy.float32(1.0) - p_orig_clipped ) ) * prior

    
    def __normalize_and_clip ( self ,
                               ratios: numpy.ndarray ,
                               w_orig: numpy.ndarray ) -> numpy.ndarray:
        
        clipped_ratios = numpy.clip ( ratios , numpy.float32(0.0), numpy.float32 ( self.__clip_threshold ) )        
        orig_sum       = numpy.sum  ( w_orig ,                  dtype = numpy.float64 )
        reweighted_sum = numpy.sum  ( w_orig * clipped_ratios , dtype = numpy.float64 )
        
        self.__norm_factor = numpy.float32 ( orig_sum / ( reweighted_sum + 1e-12 ) )
        return clipped_ratios * self.__norm_factor

    def __predict_stream_ratios ( self             ,
                                  stream_key: str  ,
                                  X: numpy.ndarray ) -> numpy.ndarray:
        
        models = self.__fitted_models [ stream_key ]
        prior  = self.__priors        [ stream_key ]

        preds_list = [ self._predict_single_model ( m , X ) for m in models ]
        p_avg      = numpy.mean ( preds_list , axis = 0 , dtype = numpy.float32 )        
        p_avg      = numpy.clip ( p_avg , numpy.float32 ( 1e-7 ) , numpy.float32 (1.0 - 1e-7 ) )
        return ( p_avg / ( numpy.float32 ( 1.0 ) - p_avg ) ) * prior

    # ==================================================================
    # Public Inference Method
    # ==================================================================
    def predict_weights ( self,
                          X_new : numpy.ndarray               ,
                          w_new : Optional[ numpy.ndarray ] = None ) -> numpy.ndarray:
        """ Apply the fitted immutable model ensemble to evaluate ratio factors r(x) on new data.
        """
        
        X_new_f32 = X_new.astype ( numpy.float32 , copy = False )
        
        # ==============================================================
        # FAST PATH for Inference without sample weights
        # ==============================================================
        if w_new is None and self.__mode == '1-stream':
            ratios  = self.__predict_stream_ratios('base', X_new_f32)
            clipped = numpy.clip ( ratios , numpy.float32 ( 0.0 ), numpy.float32 ( self.__clip_threshold ) )
            return clipped * self.__norm_factor
        
        # ==============================================================
        # SLOW PATH for explicit weights or multi-stream modes
        # ==============================================================
        w_new_f32 = w_new.astype ( numpy.float32 , copy = False ) if w_new is not None else numpy.ones ( len ( X_new ), dtype = numpy.float32 )

        nz_mask  = w_new_f32 != 0
        X_valid  , w_valid = X_new_f32 [ nz_mask ], w_new_f32 [ nz_mask ]

        mask_pos     = w_valid > 0
        mask_neg     = w_valid < 0
        ratios_valid = numpy.zeros ( len ( X_valid ), dtype = numpy.float32 )

        if self.__mode == '1-stream':
            ratios_valid = self.__predict_stream_ratios ( 'base', X_valid )

        elif self.__mode == '2-stream_orig':
            if numpy.any ( mask_pos ):
                ratios_valid [ mask_pos ] = self.__predict_stream_ratios ( 'orig_pos' , X_valid [ mask_pos ] )
            if numpy.any(mask_neg):
                ratios_valid [ mask_neg ] = self.__predict_stream_ratios ( 'orig_neg' , X_valid [ mask_neg ] )

        elif self.__mode == '2-stream_targ':
            w_pos   = self.__target_weights_info [ 'W_pos' ]
            w_neg   = self.__target_weights_info [ 'W_neg' ]
            w_total = w_pos + w_neg + numpy.float32(1e-12)

            r_pos_pos    = self.__predict_stream_ratios ( 'pos_pos' , X_valid )
            r_pos_neg    = self.__predict_stream_ratios ( 'pos_neg' , X_valid )
            ratios_valid = ( r_pos_pos * w_pos - r_pos_neg * w_neg ) / w_total

        elif self.__mode == '4-stream':
            w_pos = self.__target_weights_info [ 'W_pos' ]
            w_neg = self.__target_weights_info [ 'W_neg' ]
            w_total = w_pos + w_neg + numpy.float32(1e-12)

            if numpy.any ( mask_pos ):
                X_p      = X_valid [ mask_pos ]
                r_pos_pos = self.__predict_stream_ratios ( 'pos_pos' , X_p )
                r_pos_neg = self.__predict_stream_ratios ( 'pos_neg' , X_p )
                ratios_valid [ mask_pos ] = ( r_pos_pos * w_pos - r_pos_neg * w_neg ) / w_total

            if numpy.any ( mask_neg ):
                X_n       = X_valid [ mask_neg ]
                r_neg_pos = self.__predict_stream_ratios ( 'neg_pos' , X_n )
                r_neg_neg = self.__predict_stream_ratios ( 'neg_neg' , X_n )
                ratios_valid [ mask_neg ] = ( r_neg_pos * w_pos - r_neg_neg * w_neg ) / w_total

        clipped      = numpy.clip ( ratios_valid , numpy.floay32 ( 0.0 ), numpy.float32 ( self.__clip_threshold ) )
        final_ratios = numpy.zeros(len(w_new_f32), dtype=numpy.float32)
        final_ratios [ nz_mask ] = clipped * self.__norm_factor
        return final_ratios


# =================================================================
## @class LightGBMDensityReweighter
#   Density-ratio reweighter implementation using LightGBM as the underlying classifier.
#   Inherits from BaseDensityReweighter. Automatically handles 1-, 2-, and 4-stream
#   decompositions for positive and negative sample weights.
class LightGBMDensityReweighter ( BaseDensityReweighter ): 
    """ Density-ratio reweighter implementation using LightGBM as the underlying classifier.
    Inherits from BaseDensityReweighter. Automatically handles 1-, 2-, and 4-stream
    decompositions for positive and negative sample weights.
    """    
    def __init__(  self , * , 
                   original               : numpy.ndarray,
                   target                 : numpy.ndarray,
                   original_weight        : Optional[numpy.ndarray] = None,
                   target_weight          : Optional[numpy.ndarray] = None,
                   lgb_params             : Optional[Dict[str, Any]] = None,
                   num_boost_round        : int = 100,
                   early_stopping_rounds  : Optional[int] = 10,
                   clip_threshold         : float = 10.0,
                   n_splits               : int = 5,
                   store_original_weights : bool = True,
                 ):        
        """
        Parameters
        ----------
        original : numpy.ndarray
            Features of the original (source) dataset.
        target : numpy.ndarray
            Features of the target dataset.
        original_weight : numpy.ndarray or None, optional
            Event weights for the original dataset. If None, unweighted (all 1.0).
        target_weight : numpy.ndarray or None, optional
            Event weights for the target dataset. If None, unweighted (all 1.0).
        lgb_params : dict or None, optional
            Hyperparameters for LightGBM Booster (e.g., {'learning_rate': 0.05, 'max_depth': 6}).
        num_boost_round : int, default=100
            Number of boosting iterations.
        early_stopping_rounds : int or None, default=10
            Activates early stopping on validation fold. Set to None to disable.
        clip_threshold : float, default=10.0
            Upper bound limit for raw density ratio values r(x).
        n_splits : int, default=5
            Number of folds for StratifiedKFold cross-validation.
        random_state : int, default=42
            Random seed for reproducibility.
        store_original_weights : bool, default=True
            Whether to keep reweighted values of original sample in memory.
        """
        
        self.__num_boost_round       = int(num_boost_round)
        self.__early_stopping_rounds = early_stopping_rounds

        # Base LightGBM parameters setup
        default_lgb_params = { "objective" : "binary"          ,
                               "metric"    : "binary_logloss"  ,
                               "verbosity" : -1 ,
                               "n_jobs"    : -1 }

        if lgb_params is not None:
            default_lgb_params.update(lgb_params)
            
        self.__lgb_params = default_lgb_params
        
        # Delegate execution to BaseDensityReweighter __init__
        super().__init__ ( original               = original               ,
                           target                 = target                 ,
                           original_weight        = original_weight        ,
                           target_weight          = target_weight          ,
                           clip_threshold         = clip_threshold         ,
                           n_splits               = n_splits               ,
                           store_original_weights = store_original_weights )
        
    # =============================================================
    # Abstract Method Implementations
    # =============================================================
    def _train_single_model ( self,
                              X_train : numpy.ndarray,
                              y_train : numpy.ndarray,
                              w_train : Optional [ numpy.ndarray ] ,
                              X_val   : numpy.ndarray,
                              y_val   : numpy.ndarray,
                              w_val   : Optional [ numpy.ndarray ],
                             ) -> Tuple[ Any , numpy.ndarray ]:
        """ Trains a single LightGBM booster on train fold and evaluates on validation fold.

        Returns
        -------
        model : lgb.Booster
            Fitted LightGBM Booster instance.
        val_preds : numpy.ndarray
            Predicted probabilities p(y=1|x) on the validation fold (float32).
        """
            
        import lightgbm as LightGBM 
        
        # Construct LightGBM Datasets (pass weight only if not None)
        trn_data = LightGBM . Dataset ( X_train , label = y_train , weight = w_train , free_raw_data = False )
        val_data = LightGBM . Dataset ( X_val   , label = y_val   , weight = w_val   , free_raw_data = False , reference = trn_data)

        callbacks = []
        if self.__early_stopping_rounds is not None and self.__early_stopping_rounds > 0:
            callbacks.append ( LightGBM.early_stopping ( stopping_rounds = self.__early_stopping_rounds, verbose = False ) )

        booster = LightGBM.train ( params          = self.__lgb_params,
                                   train_set       = trn_data,
                                   num_boost_round = self.__num_boost_round,
                                   valid_sets      = [ val_data ] ,
                                   callbacks       = callbacks    )
        
        # Predict raw probabilities p(y=1|x) for validation fold
        # Best iteration is used automatically if early stopping triggered
        val_preds = booster.predict ( X_val , num_iteration = booster.best_iteration )
        ## 
        return booster, val_preds.astype ( numpy.float32 , copy = False )

    ## Generates probability predictions p(y=1|x) for input features X.
    def _predict_single_model ( self  ,
                                model ,
                                X : numpy.ndarray ) -> numpy.ndarray :
        """ Generates probability predictions p(y=1|x) for input features X.
        """
        preds = model.predict ( X , num_iteration = model.best_iteration )
        return preds.astype ( numpy.float32 , copy = False )

    
## @class Reweighter_base
#  Base class for various reweighters 
class Reweighter_base(abc.ABC) :
    """ Base class for various reweighters    
    """    
    def __init__ ( self                   , * ,
                   original               , ## original/MC sample 
                   target                 , ## target/DATA sample
                   original_weight = None , ## original/MC weights
                   target_weight   = None , ## target/DATA weights
                   ## 
                   n_splits        = 0             , ## Use Cross-vaildation for train? 
                   silent          = True          ,
                   method          = "UNSPECIFIED" , **params ) :
        
        
        config = { 'max_depth' : 3 }
        config.update ( params )
        
        assert isinstance ( n_splits , int ) and 0 <= n_splits , \
            "Invalid `n_splits':%s" % n_splits 
        
        self.__n_splits    = n_splits
        self.__method      = method
        self.__params      = config 
        self.__silent      = True if silent else False 
        
        ## MC weights are used in training?
        self.__weight_used = not weight_trivial ( original_weight )
        ## 
        self.__models      = self.__train ( original        = original        ,
                                            target          = target          ,
                                            original_weight = original_weight ,
                                            target_weight   = target_weight   )
                
    @property
    def n_splits ( self ) :
        """`n_splits` : number of splits for cross-validation (XV):  no XV if n_splits == 0 """
        return self.__n_splits
    
    @property
    def method ( self ) :
        """`method` : underlying method/engine"""
        return self.__method

    @property
    def params ( self ) :
        """`params` : configuration parameters for underlying engine"""
        return self.__params

    @property
    def models ( self ) :
        """`models`: list/tuple of trained models"""
        return self.__models 

    @property
    def weight_used ( self ) :
        """`weight_used` : was orignal weight used for training?"""
        return self.__weight_used 

    @property
    def silent ( self ) :
        """`silent` : silent processing?"""
        return self.__silent
    
    @property 
    def config ( self ) :
        """`config` : Reweighter configuraton"""
        conf = {}
        conf.update ( self.__params )
        conf [ 'method'      ] = self.__method
        conf [ 'n_splits'    ] = self.__n_splits 
        conf [ 'weight_used' ] = self.__weight_used
        conf [ 'silent'      ] = self.__silent 

    # =========================================================================
    ## self-print get the configuration 
    def table (  self , prefix = '# ') : 
        """ print configuration """
        from ostap.logger.utils import map2table_ex
        title = title if title else "%s configuration " % typename ( self ) 
        return map2table_ex ( self.config , 
                              header      = ( 'Parameter' , 'type' , 'value' ) ,
                              ailgnment   = 'rcw'  , 
                              prefix      = prefix ,
                              title       = title  )
    
    def __str__  ( self ) : return self.table ( prefix = '' ) 
    def __repr__ ( self ) : return self.__str__ ()

    # =======================================================================
    ## create & train/fit the actual model
    #  @return the model
    @abc.abstractmethod
    def model ( self   ,
                X      ,
                Y      ,
                weight = None ) :
        """ Create & train/fit the actual classifier model
        - return the model 
        """
        return NotImplemented 

    # =======================================================================
    ## call model.predict_proba 
    def predict_proba ( self  ,
                        model ,
                        sample ) :
        """ Call `model.predict_proba` """
        return model.predict_proba ( sample ) [ : , 1 ]
    
    # ======================================================================
    ## train the model
    def __train ( self                   , * , 
                  original               ,   ## original/MC sample 
                  target                 ,   ## target/DATA sample
                  original_weight = None ,   ## original/MC weights
                  target_weight   = None ) : ## target/DATA weights
        
        n_original  = len ( original )
        n_target    = len ( target   )
        
        X = numpy.vstack      ( [ original , target ] )
        Y = numpy.concatenate ( [ numpy.zeros ( n_original ), numpy.ones  ( n_target   ) ] )
        
        ow_trivial = weight_trivial ( original_weight ) 
        tw_trivial = weight_trivial ( target_weight   ) 
        
        if ow_trivial and tw_trivial : weights = None        
        else :
            
            o_weight  = numpy.ones ( n_original ) if original_weight is None  else original_weight / numpy.sum ( original_weight ) 
            t_weight  = numpy.ones ( n_target   ) if target_weight   is None  else target_weight   / numpy.sum ( targer_weight   ) 
            o_weight /= numpy.sum  ( o_weight   )
            t_weight /= numpy.sum  ( t_weight   )
            weights   = numpy.concatenate( [ o_weight , t_weight ] )

        ## clear the list of trained models 
        models = []

        # ==================================================================
        if not self.n_splits : # ===========================================
            # ==============================================================
            ## create & train/fit the model 
            model       = self.model ( X , Y , weight = weights  )        
            predictions = self.predict_proba ( model , X )             
            ## 
            models.append ( model ) 
            # ==============================================================
        else : ## Use cross-vaildation # ===================================
            # ==============================================================
            
            from sklearn.model_selection import KFold
            
            predictions = numpy.zeros ( len ( X ) )
            splits      = KFold ( n_splits     = self.n_splits ,
                                  shuffle      = True          ,
                                  random_state = self.params.get ( 'random_state' , None ) )
            
            for train_idx, val_idx in splits.split ( X , Y ) :
                X_train , Y_train , W_train = X [ train_idx ] , Y [ train_idx ], weights [ train_idx ] if not weights is None else None 
                X_val = X [ val_idx ]
                ## 
                ## create & train/fit the model 
                model = self.model ( X_train , Y_train , weight = W_train , **self.params )
                predictions [ val_idx ] = self.predict_proba ( model , X_val )                
                ## 
                models.append ( model ) 
                                
        # ============================================================================
        ## common part
        # ===========================================================================

        ## from sklearn.metrics import roc_auc_score
        ## auc          = roc_auc_score ( Y , predictions , sample_weight = weights )
        
        from ostap.stats.adval import safe_roc_auc_score
        auc          = safe_roc_auc_score ( Y , predictions , sample_weight = weights ) 
        
        print(f"[{self.method.upper()}] OOF Weighted ROC AUC: {auc:.4f}")
        
        return tuple ( models )    

    # ========================================================================
    ## get/predict new weights for (new) originals
    def weight ( self                   ,
                 original               ,
                 original_weight = None ) :
        """ Get/predict new weights for (new) originals
        """        
        if not self.silent and self.weight_used and weight_trivial ( original_weight ) :
            logger.warning ( "Reweighter: `original-weight' was used for training but not provided for evaluation" ) 
            
        n_original   = len ( original ) 
        predictions  = numpy.zeros ( n_original )        
        for model in self.models : predictions += self.predict_proba ( model , original ) 
        predictions /= len ( self.models )

        ## clip: zero is allowed! 
        predictions  = numpy.clip ( predictions , 0.0 , 1.0 - 1e-6 )
        predictions  = predictions / ( 1.0 - predictions ) 
        
        if not weight_trivial ( original_weight ) :
            predictions  = predictions * original_weight 
            predictions *= numpy.sum ( original_weight ) / np.sum( predictions )
        else :             
            predictions *= 1.0 * n_original / numpy.sum ( predictions )
            
        return predictions

    weights     = weight
    new_weight  = weight
    new_weights = weight

    # ==========================================================================
    ## Get/predict new weights for (new) original
    def __call__ ( self                   ,
                   original               ,
                   original_weight = None ) :
        """ Get/predict  new weights for (new) original
        """
        return self.weight ( original , original_weight = original_weight )
    
# ============================================================================
## @class Reweighter_LGBM
#  Actual Reweighter class based on LightGBM 
class Reweighter_LGBM(Reweighter_base):
    """ Actual Reweighted class based on LigthGBM
    """
    def __init__ ( self , * , 
                   original               , ## original/MC sample 
                   target                 , ## target/DATA sample
                   original_weight = None , ## original/MC weights
                   target_weight   = None , ## target/DATA weights
                   ## 
                   n_splits       = 0     , ## Use Cross-vaildation for train? 
                   silent         = True  , **params ) :
        
        config = { 'objective'        : 'binary',
                   'metric'           : 'binary_logloss',   # LogLoss optimizes predicted probability calibration
                    # Slow training for smooth probability transition
                    'learning_rate'    : 0.01,               # Small step size prevents aggressive jumps
                    'n_estimators'     : 300,
                    # Drastically simplified trees (forces smooth density estimation)
                    'max_depth'        : 3,                  # Shallow depth prevents isolated outlier leaves
                    'num_leaves'       : 8,                  # Small leaf count keeps probabilities close to baseline
                    'min_data_in_leaf' : 100,                # Large leaf capacity averages out extreme predictions
                    # Heavy regularization (dampens logit outputs)
                    'subsample'        : 0.8,
                    'subsample_freq'   : 1,
                    'colsample_bytree' : 0.8,
                    'reg_alpha'        : 1.0,                # Higher L1 penalty
                    'reg_lambda'       : 10.0,               # Strong L2 penalty to pull leaf outputs toward center
                    'n_jobs'           : -1,
                    'verbosity'        : -1
                }
        
        ## Attention! 
        params [ 'n_jobs' ] = num_jobs ( params , numcpu () - 1 )
        config.update ( params )
        
        ## initiailze the baze 
        Reweighter_base.__init__ ( self ,
                                   original         = original       ,
                                   target          = target          ,
                                   original_weight = original_weight ,
                                   target_weight   = target_weight   ,
                                   ##
                                   n_splits        = n_splits ,
                                   silent          = silent   ,
                                   method          = "Reweighter/LigthGBM" , **config )

    # =======================================================================
    ## call model.predict_proba 
    #  - suppress warnings from LightGBM
    def predict_proba ( self  ,
                        model ,
                        sample ) :
        """ Call `model.predict_proba`
        - suppress warnings from LightGBM
        """
        with warnings.catch_warnings():
            warnings.simplefilter ( "ignore", category = UserWarning )            
            return super().predict_proba ( model , sample )
        
    # =======================================================================
    ## create & train/fit the actual model
    #  @return the model
    def model ( self   ,
                X      ,
                Y      ,
                weight = None , **config ) :
        """ Create & train/fit the actual classifier model
        - return the model 
        """
        import lightgbm as LightGBM
        ## 
        model = LightGBM.LGBMClassifier ( **self.params )
        model.fit ( X , Y , sample_weight = weight )
        ## 
        return model

# ============================================================================
## @class Reweighter_XGB
#  Actual Reweighter class based on XGBoost 
class Reweighter_XGB(Reweighter_base):
    """ Actual Reweighted class based on XGBoost
    """
    def __init__ ( self                   , * , 
                   original               , ## original/MC sample 
                   target                 , ## target/DATA sample
                   original_weight = None , ## original/MC weights
                   target_weight   = None , ## target/DATA weights
                   ## 
                   n_splits       = 0     , ## Use Cross-vaildation for train? 
                   silent         = True  , **params ) :
           
        config = {  'objective'        : 'binary:logistic',
                    'eval_metric'      : 'logloss',          # Optimizes probability calibration
                    'tree_method'      : 'hist',
                    # Slow training for smooth probability transition
                    'learning_rate'    : 0.01,               # Small learning rate avoids aggressive probability spikes
                    'n_estimators'     : 300,
                    # Drastically simplified trees (forces broad density estimation)
                    'max_depth'        : 3,                  # Shallow depth prevents isolated outlier leaves
                    'max_leaves'       : 8,                  # Equivalent to num_leaves in LightGBM
                    'min_child_weight' : 10.0,               # High leaf threshold (~min_data_in_leaf=100 in LightGBM)
                    # Heavy regularization (dampens logit outputs)
                    'subsample'        : 0.8,
                    'colsample_bytree' : 0.8,
                    'alpha'            : 1.0,                # Stronger L1 penalty
                    'lambda'           : 10.0,               # Strong L2 penalty to pull leaf scores toward center
                    ##
                    'n_jobs'           : -1,
                    'verbosity'        : 0 }
        
        ## Attention! 
        params [ 'n_jobs' ] = num_jobs ( params , numcpu () - 1 )
        config.update ( params )
        
        ## initiailze the baze 
        Reweighter_base.__init__ ( self ,
                                   original         = original       ,
                                   target          = target          ,
                                   original_weight = original_weight ,
                                   target_weight   = target_weight   ,
                                   ##
                                   n_splits        = n_splits ,
                                   silent          = silent   ,
                                   method          = "Reweighter/XGBoost" , **config )

    # =======================================================================
    ## create & train/fit the actual model
    #  @return the model
    def model ( self   ,
                X      ,
                Y      ,
                weight = None ) :
        """ Create & train/fit the actual classifier model
        - return the model 
        """
        import xgboost as XGBoost
        ## 
        model = XGBoost.XGBClassifier ( **self.params )
        model.fit ( X , Y , sample_weight = weight )
        ## 
        return model
    
# ============================================================================
## @class Reweighter_CATB
#  Actual Reweighter class based on CatBoost
class Reweighter_CATB(Reweighter_base):
    """ Actual Reweighted class based on CatBoost
    """
    def __init__ ( self                   , * , 
                   original               , ## original/MC sample 
                   target                 , ## target/DATA sample
                   original_weight = None , ## original/MC weights
                   target_weight   = None , ## target/DATA weights
                   ## 
                   n_splits       = 0     , ## Use Cross-vaildation for train? 
                   silent         = True  , **params ) :

        
        config = {  'loss_function'    : 'Logloss',
                    'eval_metric'      : 'Logloss',          # Ensures probability calibration
                    # Slow training for smooth probability transitions
                    'learning_rate'    : 0.01,               # Small learning rate prevents sudden probability jumps
                    'iterations'       : 300,
                    # Drastically simplified trees (forces smooth density estimation)
                    'max_depth'        : 3,                  # Shallow depth prevents isolated outlier leaves
                    'min_data_in_leaf' : 100,                # High leaf threshold averages out extreme probabilities
                    # Heavy L2 regularization (dampens logit outputs)
                    'l2_leaf_reg'      : 15.0,               # Strong L2 penalty to pull leaf scores toward the center
                    'subsample'        : 0.8,
                    'bootstrap_type'   : 'Bernoulli',
                    ##
                    'thread_count'     : -1,
                    'verbose'          : False
                }
        
        ## Attention! 
        params [ 'thread_count' ] = num_jobs ( params , numcpu () - 1 )
        config.update ( params )
        
        ## initialize the baze 
        Reweighter_base.__init__ ( self ,
                                   original         = original       ,
                                   target          = target          ,
                                   original_weight = original_weight ,
                                   target_weight   = target_weight   ,
                                   ##
                                   n_splits        = n_splits ,
                                   silent          = silent   ,
                                   method          = "Reweighter/CatBoost" , **config )
        
    # =======================================================================
    ## create & train/fit the actual model
    #  @return the model
    def model ( self   ,
                X      ,
                Y      ,
                weight = None ) :
        """ Create & train/fit the actual classifier model
        - return the model 
        """
        import catboost as CatBoost
        ##   
        model = CatBoost.CatBoostClassifier ( **self.params )    
        model.fit ( CatBoost.Pool ( X , Y , weight = weight ) ) 
        ## 
        return model

# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        
# =============================================================================
##                                                                      The END 
# =============================================================================
