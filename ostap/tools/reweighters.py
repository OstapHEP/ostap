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
from   ostap.math.math_base   import weight_trivial
from   ostap.core.ostap_types import num_types 
from   ostap.utils.core       import typename
from   ostap.utils.basic      import numcpu, num_jobs
import numpy, abc, warnings  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighters' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## Number of features for training data
#  Safe extraction of n_features (supports DataFrame, 2D array, 1D vector)
def num_features ( X ) :
    """ Number of features for training data
    - Safe extraction of n_features (supports DataFrame, 2D array, 1D vector)
    """
    return X.shape[1] if hasattr ( X , 'shape' ) and 1 < len ( X.shape ) else 1
# ============================================================================
## Returns the number of samples/events (rows) in the dataset.
#  Supports NumPy arrays, Pandas DataFrames/Series, SciPy sparse matrices,
#  PyTorch/TensorFlow tensors, and standard Python collections.
def num_samples ( X ) :
    """ Returns the number of samples/events (rows) in the dataset.

    Supports NumPy arrays, Pandas DataFrames/Series, SciPy sparse matrices,
    PyTorch/TensorFlow tensors, and standard Python collections.
    """
    if X is None : return 0

    # Handle objects with a .shape attribute (NumPy, Pandas, SciPy, PyTorch, TensorFlow)
    if hasattr ( X ,  "shape" ) :
        shape = X.shape
        if 0 == len ( shape ) :  return 1
        return int ( shape [ 0 ] )

    # Handle lists, tuples, and other standard Python sequences
    if hasattr ( X , "__len__" ) : return len ( X )

    raise TypeError ( "Unsupported data structure for determining sample count: %" % typename ( X )  )
# =============================================================================
## @class DensityReweighter
#  Abstract base class for immutable, adaptive density-ratio reweighting.
#  Implements 1-, 2-, and 4-stream signed-measure density estimations.
# 
#  The instance executes training immediately upon instantiation (__init__)
#  and seals the resulting weights for the original sample as read-only attributes
#  (if store_original_weights=True).
class DensityReweighter(abc.ABC):
    """ Abstract base class for immutable, adaptive density-ratio reweighting.
    Implements 1-, 2-, and 4-stream signed-measure density estimations.
    
    The instance executes training immediately upon instantiation (__init__)
    and seals the resulting weights for the original sample as read-only attributes
    (if store_original_weights=True).
    """

    def __init__( self, * , 
                  original               ,
                  target                 ,
                  original_weight        = None ,
                  target_weight          = None ,
                  clip_threshold         = 10.0 ,
                  n_splits               = 5    ,
                  random_state           = 42   , 
                  store_original_weights = True , **params ) :
        
        if not isinstance ( n_splits       , int )       : raise TypeError  ( "Invalid `n_splits' type  %s"       % typename ( n_splits ) )
        if not 0 <= n_splits <= 1000                     : raise ValueError ( "Invalid `n_splits' value %s"       % n_splits )     
        if not isinstance ( random_state   , int )       : raise TypeError  ( "Invalid `random_state' type %s"    % typename ( random_state ) )
        if not isinstance ( clip_threshold , num_types ) : raise TypeError  ( "Invalid `clip_threshold' type %s"  % typename ( clip_threshold ) )
        if not 0 < clip_threshold                        : raise ValueError ( "Invalid `clip_threshold' value %s" % clip_threshold )
        
        self.__clip_threshold  = float ( clip_threshold )  
        self.__n_splits        = n_splits  
        self.__random_state    = random_state
        self.__params          = params 
        
        # Internal storage for fitted models and stream parameters
        self.__fitted_models       = {}
        self.__priors              = {}
        self.__target_weights_info = {}
        self.__norm_factor         = numpy.float32(1.0)
        self.__mode                = None
        
        # Execute training immediately and compute original sample weights (float32)
        original_ratios, original_reweighted_weights = self.__fit_and_compute (
            original        .astype ( numpy.float32 , copy = False ) ,
            original_weight .astype ( numpy.float32 , copy = False ) if original_weight is not None else None,
            target          .astype ( numpy.float32 , copy = False ) ,
            target_weight   .astype ( numpy.float32 , copy = False ) if target_weight is not None else None,
        )
        
        self.__original_ratios             = None
        self.__original_reweighted_weights = None
        
        # Optionally store original sample reweighting results
        if store_original_weights :
            self.__original_ratios             = original_ratios
            self.__original_reweighted_weights = original_reweighted_weights

    # ==================================================================
    # Public Read-Only Properties
    # ==================================================================
    @property
    def original_ratios ( self ) :
        """ Computed density ratio factors r(x) for the original sample.
        Returns None if store_original_weights=False during initialization.
        """
        return self.__original_ratios

    @property
    def original_reweighted_weights ( self ) : 
        """ Final reweighted event weights (w * r) for the original sample.
        Returns None if store_original_weights=False during initialization.
        """
        return self.__original_reweighted_weights

    @property
    def mode ( self ) :
        """ Active stream decomposition scheme ('1-stream', '2-stream_orig', '2-stream_targ', or '4-stream')."""
        return self.__mode
    
    @property 
    def params ( self ) :       
        """`params` : dictionary of hyperparameters for the underlying classifier."""
        return self.__params
    
    @property 
    def config ( self ) :
        """`config` : Reweighter configuraton"""
        conf = {}
        conf.update ( self.__params )
        conf [ 'mode'            ] = self.__mode
        conf [ 'clip_threshold'  ] = self.__clip_threshold
        conf [ 'n_splits'        ] = self.__n_splits 
        conf [ 'original_ratios' ] = True if self.__original_ratios is not None else False
        conf [ 'original_reweighted_weights' ] = True if self.__original_reweighted_weights is not None else False 
        conf [ 'random_state'    ] = self.__random_state 
    
    # =========================================================================
    ## self-print get the configuration 
    def table (  self , title = '' , prefix = '# ') : 
            """ print configuration """
            from ostap.logger.utils import map2table_ex
            title = title if title else "%s configuration " % typename ( self ) 
            return map2table_ex ( self.config , 
                                  header      = ( 'Parameter' , 'type' , 'value' ) ,
                                  alignment   = 'rcw'  , 
                                  prefix      = prefix ,
                                  title       = title  )
        
    def __str__  ( self ) : return self.table ( prefix = '' ) 
    def __repr__ ( self ) : return self.__str__ ()
    
    # ==================================================================
    # Abstract Methods (Must be implemented by child framework wrappers)
    # ==================================================================
    @abc.abstractmethod
    def _train_single_model( self ,
                             X_train , y_train , w_train ,
                             X_val   , y_val   , w_val   ) : 
        """ Train a single base classifier fold and return (model, val_predictions)."""
        raise NotImplementedError

    ## Predict probability p(y=1|x) using a single trained model.
    @abc.abstractmethod
    def _predict_single_model ( self  , model , X ) : 
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

    def use_strong_regularization ( self , X ) :
        """ Checks if dataset size is small relative to feature count,
        indicating a need for stronger regularization.
        """
        ns = num_samples   ( X )
        nf = num_features  ( X )

        return ns < max ( 300 if nf <= 3  else 1000 , nf * 50 )
    
    # ==================================================================
    # Internal Fitting Pipeline
    # ==================================================================
    def __fit_and_compute ( self, 
                            original , original_weight ,
                            target   , target_weight   ) : 

        # ==============================================================
        # FAST PATH: Standard positive unit weights for both samples
        # ==============================================================
        if original_weight is None and target_weight is None:
            
            self.__mode  = '1-stream'
            ratios_valid = self.__fit_eval_stream ( 'base' , original , None , target , None  )
            
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
            ratios_valid = self.__fit_eval_stream ( 'base' , X_orig_valid , w_orig_valid , X_targ_valid , w_targ_valid )

        # ==============================================================
        # SCENARIO 2: Negative weights ONLY in Original sample
        # ==============================================================
        elif has_neg_orig and not has_neg_targ:
            self.__mode  = '2-stream_orig'
            ratios_valid = numpy.zeros ( len ( X_orig_valid ), dtype = numpy.float32 )
            ratios_valid [ mask_orig_pos ] = self.__fit_eval_stream ( 'orig_pos' , X_orig_valid [ mask_orig_pos ] ,
                                                                      w_orig_valid [ mask_orig_pos ] ,
                                                                      X_targ_valid ,
                                                                      w_targ_valid ) 
            ratios_valid [ mask_orig_neg ] = self.__fit_eval_stream ( 'orig_neg' , X_orig_valid [ mask_orig_neg ] ,
                                                                      w_orig_valid [ mask_orig_neg ] ,
                                                                      X_targ_valid ,
                                                                      w_targ_valid )
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

            r_pos_pos = self.__fit_eval_stream ( 'pos_pos' , X_orig_valid , w_orig_valid , X_targ_pos , w_targ_pos )
            r_pos_neg = self.__fit_eval_stream ( 'pos_neg' , X_orig_valid , w_orig_valid , X_targ_neg , w_targ_neg )

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

            r_pos_pos = self.__fit_eval_stream ( 'pos_pos' , X_orig_pos , w_orig_pos , X_targ_pos , w_targ_pos )
            r_pos_neg = self.__fit_eval_stream ( 'pos_neg' , X_orig_pos , w_orig_pos , X_targ_neg , w_targ_neg )
            r_neg_pos = self.__fit_eval_stream ( 'neg_pos' , X_orig_neg , w_orig_neg , X_targ_pos , w_targ_pos )
            r_neg_neg = self.__fit_eval_stream ( 'neg_neg' , X_orig_neg , w_orig_neg , X_targ_neg , w_targ_neg )

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

    def __fit_eval_stream ( self       , stream_key ,
                            X_orig_sub , w_orig_sub ,
                            X_targ_sub , w_targ_sub ) :  
        
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
    def predict_weights ( self, X_new , w_new = None ) :
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

        clipped      = numpy.clip ( ratios_valid , numpy.float32 ( 0.0 ), numpy.float32 ( self.__clip_threshold ) )
        final_ratios = numpy.zeros(len(w_new_f32), dtype=numpy.float32)
        final_ratios [ nz_mask ] = clipped * self.__norm_factor
        return final_ratios


    weight      = predict_weights
    weights     = predict_weights
    new_weight  = predict_weights
    new_weights = predict_weights

    # ==========================================================================
    ## Get/predict new weights for (new) original
    def __call__ ( self                   ,
                   original               ,
                   original_weight = None ) :
        """ Get/predict  new weights for (new) original
        """
        return self.weight ( original , original_weight = original_weight )
    
# =================================================================
## @class LightGBMDensityReweighter
#   Density-ratio reweighter implementation using LightGBM as the underlying classifier.
#   Inherits from BaseDensityReweighter. Automatically handles 1-, 2-, and 4-stream
#   decompositions for positive and negative sample weights.
class LightGBMDensityReweighter ( DensityReweighter ): 
    """ Density-ratio reweighter implementation using LightGBM as the underlying classifier.
    Inherits from BaseDensityReweighter. Automatically handles 1-, 2-, and 4-stream
    decompositions for positive and negative sample weights.
    """    
    def __init__(  self , * , 
                   original               ,
                   target                 ,
                   original_weight        = None  ,
                   target_weight          = None  ,
                   store_original_weights = True  , **params ,
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
        num_boost_round : int, default=100
            Number of boosting iterations.  
    
        """    
        ## Recommended LightGBM parameters for density ratio reweighting (covariate shift)
        config = { 
            # Core objective and evaluation metric
            "objective"             : "binary" ,
            "metric"                : "binary_logloss" ,
            "n_estimators"          : 500 ,
            "early_stopping_rounds" :  20 ,
            # 1. Tree structure regularization (prevents overconfident density ratios)
            "max_depth"             : 5  , # Shallow trees to avoid overfitting local sample variations
            "num_leaves"            : 31 , # Upper bound on leaf nodes (<= 2^max_depth)
            "min_child_samples"     : 50 , # Minimum samples per leaf (smooths probability estimates)
            # 2. Learning rate control
            "learning_rate"         : 0.03 , # Small step size for smoother probability boundaries
            # 3. Regularization on leaf weights
            "reg_alpha"             : 0.1 , # L1 regularization to prune noisy feature splits
            "reg_lambda"            : 1.0 , # L2 regularization to stabilize raw ratios p / (1 - p)
            # 4. Stochastic subsampling
            "subsample"             : 0.8 , # Row subsampling fraction (bagging)
            "subsample_freq"        : 1   , # Perform bagging at every iteration
            "colsample_bytree"      : 0.8 , # Feature subsampling fraction per tree
            # 6. Runtime parameters
            "verbosity"             : -1 ,
            "n_jobs"                : -1 }
        
        config.update ( params )
        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' )
        
        # Delegate execution to BaseDensityReweighter __init__
        super().__init__ ( original               = original               ,
                           target                 = target                 ,
                           original_weight        = original_weight        ,
                           target_weight          = target_weight          ,
                           store_original_weights = store_original_weights , **config )
        
    # =============================================================
    # Abstract Method Implementations
    # =============================================================
    def _train_single_model ( self,
                              X_train , y_train , w_train ,
                              X_val   , y_val   , w_val   ) : 
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

        params = {}
        params.update ( self.params ) 
        
        num_boost_round       = params.pop ( 'num_boost_round' , None ) or params.pop ( 'n_estimators' ,  None  ) or 500 
        if isinstance ( num_boost_round , int ) and 10 < num_boost_round < 10000 : pass 
        else                                                                     : num_boost_round = 500 
        
        early_stopping_rounds = params.pop  ( 'early_stopping_rounds' , 10 )
        if isinstance ( early_stopping_rounds , int ) and 1 < early_stopping_rounds < num_boost_round : pass
        else                                                                      :  early_stopping_rounds = 10  
        
        if self.use_strong_regularization ( X_train ) :  
            max_depth  = min ( 3 , params.pop ( 'max_depth' ,  5 ) ) 
            num_leaves = 2 ** max_depth - 1
            params [ 'max_depth'         ] = max_depth 
            params [ 'num_leaves'        ] = min ( num_leaves , params.pop ( 'num_leaves'  , 31 ) )
            params [ 'min_child_samples' ] = max ( 75 , params.pop ( 'min_child_samples' , 50 ) )
            params [ 'colsample_bytree'  ] = 1.0 
            params [ 'reg_alpha'         ] = max ( 0.5  , params.pop ( 'reg_alpha'     , 0.1  ) )
            params [ 'reg_lambda'        ] = max ( 2.0  , params.pop ( 'reg_lambda'    , 1.0  ) )
            params [ "learning_rate"     ] = min ( 0.02 , params.get ( "learning_rate" , 0.03 ) )
            num_boost_round       = min ( 150 , num_boost_round       )
            early_stopping_rounds = min (  10 , early_stopping_rounds )
            
        callbacks = []
        if 0 < early_stopping_rounds : 
            callbacks.append ( LightGBM.early_stopping ( stopping_rounds = early_stopping_rounds, verbose = False ) )

        model = LightGBM.train ( params          = params          ,
                                 train_set       = trn_data        ,
                                 num_boost_round = num_boost_round ,
                                 valid_sets      = [ val_data ]    ,
                                 callbacks       = callbacks       )
        
        # Predict raw probabilities p(y=1|x) for validation fold
        # Best iteration is used automatically if early stopping triggered
        best_iteration = model.best_iteration if 0 < getattr ( model , 'best_iteration' , 0 ) else None
        val_preds = model.predict ( X_val , num_iteration = best_iteration )
        ## 
        return model, val_preds.astype ( numpy.float32 , copy = False )

    ## Generates probability predictions p(y=1|x) for input features X.
    def _predict_single_model ( self  , model , X ) : 
        """ Generates probability predictions p(y=1|x) for input features X.
        """
        best_iteration = model.best_iteration if 0 < getattr ( model, 'best_iteration', 0 ) else None
        preds          = model.predict ( X , num_iteration = best_iteration )
        return preds.astype ( numpy.float32 , copy = False )

# =================================================================
## @class XGBoostDensityReweighter
#   Density-ratio reweighter implementation using XGBoost as the underlying classifier.
#   Inherits from BaseDensityReweighter. Automatically handles 1-, 2-, and 4-stream
#   decompositions for positive and negative sample weights.
class XGBoostDensityReweighter ( DensityReweighter ): 
    """ Density-ratio reweighter implementation using XGBoost as the underlying classifier.
    Inherits from BaseDensityReweighter. Automatically handles 1-, 2-, and 4-stream
    decompositions for positive and negative sample weights.
    """    
    def __init__(  self , * , 
                   original               ,
                   target                 ,
                   original_weight        = None ,
                   target_weight          = None ,
                   store_original_weights = True , **params ) :
        
        # Recommended XGBoost parameters for density ratio reweighting (covariate shift)
        config = { # Core objective and metric
            "objective"             : "binary:logistic" ,
            "eval_metric"           : "logloss" ,
            "n_estimators"          : 500 ,
            "early_stopping_rounds" :  20 ,
            # 1. Tree structure regularization (prevents overconfident density ratios)
            "max_depth"             : 5   , # Shallow trees prevent overfitting to local sample noise
            "min_child_weight"      : 5.0 , # Minimum sum of instance weight in a child (smooths density estimation)    
            # 2. Learning rate control
            "learning_rate"         : 0.03 , # Small step size yields smoother probability boundaries
            # 3. Regularization on leaf weights (L1 / L2)
            "alpha"                 : 0.1 , # L1 regularization (reg_alpha equivalent)
            "lambda"                : 1.0 , # L2 regularization (reg_lambda equivalent) to stabilize raw ratios
            # 4. Stochastic subsampling
            "subsample"             : 0.8 , # Row subsampling fraction
            "colsample_bytree"      : 0.8 , # Feature subsampling fraction per tree
            # 5. Runtime parameters
            "verbosity"             :  0  ,
            "n_jobs"                : -1  }
        
        config.update ( params )

        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' )

        # Delegate execution to BaseDensityReweighter __init__
        super().__init__ ( original               = original               ,
                           target                 = target                 ,
                           original_weight        = original_weight        ,
                           target_weight          = target_weight          ,
                           store_original_weights = store_original_weights , **config )
        
    # =============================================================
    # Abstract Method Implementations
    # =============================================================
    def _train_single_model ( self,
                              X_train , y_train , w_train ,
                              X_val  , y_val    , w_val   ) :
        """ Trains a single XGBoost booster on train fold and evaluates on validation fold.

        Returns
        -------
        model : xgb.Booster
            Fitted XGBoost Booster instance.
        val_preds : numpy.ndarray
            Predicted probabilities p(y=1|x) on the validation fold (float32).
        """
            
        import xgboost as XGBoost
        
        # Construct XGBoost DMatrix instances
        dtrain = XGBoost.DMatrix ( X_train , label = y_train , weight = w_train )
        dval   = XGBoost.DMatrix ( X_val   , label = y_val   , weight = w_val   )

        evals = [ ( dval , "val" ) ]

        params = {}
        params.update ( self.params )
        
        num_boost_round       = params.pop ( 'num_boost_round' , None ) or params.pop ( 'n_estimators' ,  None  ) or 500 
        if isinstance ( num_boost_round , int ) and 10 < num_boost_round < 10000 : pass 
        else                                                                     : num_boost_round = 500 
                
        early_stopping_rounds = params.pop  ( 'early_stopping_rounds' , 10 )
        if isinstance ( early_stopping_rounds , int ) and 1 < early_stopping_rounds < num_boost_round : pass
        else : early_stopping_rounds = 10 
        
        if self.use_strong_regularization ( X_train ) :  
            params [ 'max_depth'         ] = min (  3 , params.pop ( 'max_depth'   ,  5 ) ) 
            params [ 'min_child_weight'  ] = max ( 20 , params.pop ( 'min_child_weight' , 5 ) )
            params [ 'colsample_bytree'  ] = 1.0 
            params [ 'alpha'             ] = max ( 0.5 , params.pop ( 'alpha'  , 0.1 ) )
            params [ 'lambda'            ] = max ( 2.0 , params.pop ( 'lambda' , 1.0 ) )
            params [ "learning_rate"     ] = min ( 0.02 , params.get ( "learning_rate" , 0.03 ) )
            num_boost_round       = min ( 150 , num_boost_round       )
            early_stopping_rounds = min (  10 , early_stopping_rounds )

        model = XGBoost.train ( params                = params,
                                dtrain                = dtrain,
                                num_boost_round       = num_boost_round,
                                evals                 = evals,
                                early_stopping_rounds = early_stopping_rounds,
                                verbose_eval          = False )
    
        
        best_iter  = getattr ( model , "best_iteration" , None)
        iter_range = ( 0 , best_iter + 1 ) if best_iter is not None and 0 < best_iter  else (0, 0)
        val_preds  = model.predict ( dval, iteration_range = iter_range )
        return model, val_preds.astype ( numpy.float32 , copy = False )
    
    def _predict_single_model ( self  ,
                                model ,
                                X     ):
        """ Generates probability predictions p(y=1|x) for input features X.
        """
        
        import xgboost as XGBoost

        dmatrix    = XGBoost.DMatrix(X)
        best_iter  = getattr ( model , "best_iteration" , None )
        iter_range = ( 0 , best_iter + 1 ) if best_iter is not None and 0 < best_iter  else (0, 0)
        preds      = model.predict(dmatrix, iteration_range=iter_range)
        return preds.astype(numpy.float32, copy=False)
        
# =================================================================
## @class CatBoostDensityReweighter
#   Density-ratio reweighter implementation using CatBoost as the underlying classifier.
#   Inherits from BaseDensityReweighter. Automatically handles 1-, 2-, and 4-stream
#   decompositions for positive and negative sample weights.
class CatBoostDensityReweighter ( DensityReweighter ): 
    """ Density-ratio reweighter implementation using CatBoost as the underlying classifier.
    Inherits from BaseDensityReweighter. Automatically handles 1-, 2-, and 4-stream
    decompositions for positive and negative sample weights.
    """    
    def __init__(  self , * , 
                   original               ,
                   target                 ,
                   original_weight        = None ,
                   target_weight          = None ,
                   store_original_weights = True  , **params ) :        
        
        # Recommended CatBoost parameters for density ratio reweighting (covariate shift)
        config = { # Core objective and metric
            "loss_function"         : "Logloss" ,
            "eval_metric"           : "Logloss" ,
            "n_estimators"          : 500 ,
            "early_stopping_rounds" :  20 , 
            "min_child_samples"     :   5 , 
            # 1. Tree structure regularization (prevents overconfident density ratios)
            "max_depth"             : 6    , # Standard depth for oblivious trees
            "l2_leaf_reg"           : 3.0  , # L2 regularization on leaf values (smooths probability estimates)
            "random_strength"       : 1.0 , # Randomness scoring split candidates to combat overfitting
            # 2. Learning rate control
            "learning_rate"         : 0.03 , # Small step size for smoother probability boundaries
            # 3. Stochastic subsampling
            "subsample"             : 0.8  , # Row subsampling fraction
            # 5. Runtime settings
            "verbose"               : False ,
            "thread_count"          : -1    }

        config.update ( params )
        
         # CatBoost uses 'thread_count' instead of 'n_jobs'
        if 'n_jobs'     in config : config [ 'thread_count' ] = config.pop ( 'n_jobs' )  
        if 'iterations' in config : config [ 'n_estimators' ] = config.pop ( 'iterations' ) 
                   
        # Delegate execution to BaseDensityReweighter __init__
        super().__init__ ( original               = original               ,
                           target                 = target                 ,
                           original_weight        = original_weight        ,
                           target_weight          = target_weight          ,
                           store_original_weights = store_original_weights , **config )
        
    # =============================================================
    # Abstract Method Implementations
    # =============================================================
    def _train_single_model ( self    ,
                              X_train , y_train , w_train ,
                              X_val  , y_val    ,  w_val  ) :
        """ Trains a single CatBoostClassifier on train fold and evaluates on validation fold.

        Returns
        -------
        model : cb.CatBoostClassifier
            Fitted CatBoostClassifier instance.
        val_preds : numpy.ndarray
            Predicted probabilities p(y=1|x) on the validation fold (float32).
        """
            
        import catboost as CatBoost
        
        trn_pool = CatBoost.Pool ( X_train , label = y_train , weight = w_train )
        val_pool = CatBoost.Pool ( X_val   , label = y_val   , weight = w_val   )

        params = {}
        params.update ( self.params )
        
        iterations       = params.pop ( 'iterations' , None ) or params.pop ( 'n_estimators' ,  None  ) or 500 
        if isinstance ( iterations , int ) and 10 < iterations < 10000 : pass 
        else   : iterations = 500 
                
        early_stopping_rounds = params.pop  ( 'early_stopping_rounds' , 10 )
        if isinstance ( early_stopping_rounds , int ) and 1 < early_stopping_rounds < iterations : pass
        else  : early_stopping_rounds = 10
        
        if self.use_strong_regularization ( X_train ) : 
            params [ 'max_depth'         ] = min (  3   , params.pop ( 'max_depth'         , 6 ) ) 
            params [ 'min_child_samples' ] = max ( 20   , params.pop ( 'min_child_samples' , 5 ) )
            params [ 'rsm'               ] = 1.0
            params [ 'l2_leaf_reg'       ] = max ( 2.0  , params.pop ( 'l2_leaf_reg'       , 3.0 ) )
            params [ "learning_rate"     ] = min ( 0.02 , params.get ( "learning_rate"     , 0.03 ) )
            iterations       = min ( 150 , iterations       )
            early_stopping_rounds = min (  10 , early_stopping_rounds )
        
        params [ 'iterations'            ] = iterations
        params [ 'early_stopping_rounds' ] = early_stopping_rounds  
        params [ 'use_best_model'        ] = True
                
        model = CatBoost.CatBoostClassifier ( **params ) 

        model.fit ( trn_pool ,
                    eval_set = val_pool ,
                    verbose  = False    )

        # Predict probability for positive class p(y=1|x)
        val_preds = model.predict_proba ( val_pool ) [ : , 1 ]
        return model, val_preds.astype ( numpy.float32 , copy = False )

    def _predict_single_model ( self  ,
                                model ,
                                X     ):
        """ Generates probability predictions p(y=1|x) for input features X.
        """
        preds = model.predict_proba ( X ) [ : , 1 ]
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
    def table ( self , title = '' , prefix = '# ') : 
        """ print configuration """
        from ostap.logger.utils import map2table_ex
        title = title if title else "%s configuration " % typename ( self ) 
        return map2table_ex ( self.config , 
                              header    = ( 'Parameter' , 'type' , 'value' ) ,
                              alignment = 'rcw'  , 
                              prefix    = prefix ,
                              title     = title  )
    
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
        else : ## Use cross-validaton # ====================================
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
                   n_splits       = 0     , ## Use Cross-validation for train? 
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
                   n_splits       = 0     , ## Use Cross-validation for train? 
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
