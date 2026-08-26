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
    'LightGBMDensityReweighter' , ## LightGBM-based density reweighter
    'XGBoostDensityReweighter'  , ## XGboost-based  density reweighter
    'CatBoostDensityReweighter' , ## CatBoost       density reweighter
    'GBReweighter'              , ## Reweighter based on GBReweighter from hep_ml 
) 
# =============================================================================
from   ostap.core.ostap_types import num_types 
from   ostap.utils.core       import typename
from   ostap.utils.basic      import numcpu, num_jobs, NoContext 
from   ostap.tools.reweighter import Reweighter
from   ostap.stats.utils      import ( weight_trivial     ,
                                       valid_weight       ,
                                       valid_data_shape   ,
                                       compatible_weights , 
                                       num_features       ,
                                       num_samples        ,
                                       check_all          ,
                                       nEff               ) 
import numpy, abc, warnings
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighters' )
else                       : logger = getLogger( __name__ )
# =============================================================================
DEFAULT_ESTIMATORS     = 500
REGULARIZED_ESTIMATORS = 400
EPS1                   = 1.e-3
EPS2                   = 1.e-3
# =============================================================================
## Check if strong regularization is needed considering BOTH original and target samples.
#  Evaluates sample dimensionality, Kish's effective sample size (nEff),
#  and the fraction of negative sPlot weights.
def RW_needs_regularization ( original                       ,
                              target                         ,
                              original_weight        = None  , 
                              target_weight          = None  ) :
    """ Check if strong regularization is needed for BDT-based reweighting tests.
        Evaluates the limiting effective statistics between original and target samples,
        phase space dimensionality, and background subtraction weight fluctuations.
    """
    check_all ( data1   = original                  ,
                data2   = target                    ,
                weight1 = original_weight           ,
                weight2 = target_weight             ,
                where   = "RW_needs_regularization" )   


    nf = num_features ( original )
    
    # 1. Low dimensionality (<= 4 features) ALWAYS needs strong regularization
    if nf <= 4 : return True 
    
    # 2. Raw sample sizes (nRaw)
    nraw_orig = num_samples ( original )
    nraw_targ = num_samples ( target   )
    
    # 3. Calculate effective statistics (nEff)
    neff_orig = nEff ( original , original_weight )
    neff_targ = nEff ( target   , target_weight   )
    neff      = min  ( neff_orig , neff_targ       )
    
    # 4. Non-linear density threshold (phase space growth)
    required_stats = 1500.0 * ( nf ** 1.8 )
    if neff < required_stats : return True
    
    # 5. Check weight efficiency ratios (nEff / nRaw)
    #    An efficiency below 65% corresponds to >10% negative sPlot weight noise or high weight variance
    eff_orig = neff_orig / float ( nraw_orig ) if 0 < nraw_orig else 0.0
    eff_targ = neff_targ / float ( nraw_targ ) if 0 < nraw_targ else 0.0
    
    if eff_orig < 0.65 or eff_targ < 0.65 : return True

    return False

# =============================================================================
## @class DensityReweighter
#  Abstract base class for immutable, adaptive density-ratio reweighting.
#  Implements 1-, 2-, and 4-stream signed-measure density estimations.
# 
#  The instance executes training immediately upon instantiation (__init__)
#  and seals the resulting weights for the original sample as read-only attributes
#  (if store_original_weights=True).
# =============================================================================
class DensityReweighter ( Reweighter, abc.ABC ) :
    """ Abstract base class for immutable, adaptive density-ratio reweighting.
    Implements 1-, 2-, and 4-stream signed-measure density estimations.

    The instance executes training immediately upon instantiation (__init__)
    and seals the resulting weights for the original sample as read-only attributes
    (if store_original_weights=True).
    """

    def __init__( self                          , * ,
                 original                       ,
                 target                         ,
                 original_weight        = None  , 
                 target_weight          = None  ,
                 clip_threshold         = 1.e+4 ,
                 n_splits               = 5     , ## n-fold!
                 random_state           = 42    ,
                 store_original_weights = True  , **params ) :

        if not isinstance ( n_splits, int ) : raise TypeError ( "Invalid `n_splits' type %s" % typename( n_splits ) )
        if not 0 <= n_splits <= 1000        : raise ValueError ( "Invalid `n_splits' value %s" % n_splits )
        if not isinstance( clip_threshold, num_types ) :
            raise TypeError ( "Invalid `clip_threshold' type %s" % typename( clip_threshold ) )
        if not 0 < clip_threshold :
            raise ValueError ( "Invalid `clip_threshold' value %s" % clip_threshold )
        
        self.__clip_threshold = float ( clip_threshold )
        self.__n_splits       = n_splits
        self.__random_state   = random_state

        # Internal storage for fitted models and stream parameters
        self.__fitted_models       = {}
        self.__priors              = {}
        self.__target_weights_info = {}
        self.__norm_factor         = numpy.float32( 1.0 )
        self.__mode                = None
        self.__scale_factors       = {}

        ## perform regularization!! 
        if self.needs_regularization ( original        = original        ,
                                       target          = target          ,
                                       original_weight = original_weight ,
                                       target_weight   = target_weight   ) :
            
            n_features = num_features ( original )
            neff_orig  = nEff ( original  , original_weight )
            neff_targ  = nEff ( target    , target_weight   ) 
            neff       = min  ( neff_orig , neff_targ )
            
            params.update ( self.regularization ( params , n_features , neff ) )

            ESR = params.get ( 'early_stopping_rounds' ,  15 )
            if  ESR is None : params [ 'early_stopping_rounds' ] = None 
            else            : params [ 'early_stopping_rounds' ] = min ( 15 , ESR )
            
            logger.attention ( "%s: strong regularization is applied" % typename ( self ) ) 

        self.__original_ratios             = None
        self.__original_reweighted_weights = None
            
        # =====================================================================
        ## Initialize the base: check input data & print config
        # =====================================================================
        Reweighter.__init__ ( self            ,
                              original        = original        ,
                              target          = target          , 
                              original_weight = original_weight ,
                              target_weight   = target_weight   , **params )

        # =====================================================================
        ## Execute training immediately and compute original sample weights (float32)
        # =====================================================================
        original_ratios, original_reweighted_weights = (
            self.__fit_and_compute(
                original.astype( numpy.float32, copy = False ),
                (
                    original_weight.astype( numpy.float32, copy = False )
                    if original_weight is not None
                    else None
                ),
                target.astype( numpy.float32, copy = False ),
                (
                    target_weight.astype( numpy.float32, copy = False )
                    if target_weight is not None
                    else None
                ),
            )
        )

        # =====================================================================
        ## Optionally store original sample reweighting results
        # =====================================================================
        if store_original_weights :
            self.__original_ratios             = original_ratios
            self.__original_reweighted_weights = original_reweighted_weights

    @property
    def n_splits ( self ) :
        """`n_splits` : number of splits/folds, same as `n_folds`"""
        return self.__n_splits

    # =============================================================================
    ## Check if strong regularization is needed considering BOTH original and target samples.
    #  Evaluates sample dimensionality, Kish's effective sample size (nEff),
    #  and the fraction of negative sPlot weights.
    def needs_regularization ( self ,
                               original                       ,
                               target                         ,
                               original_weight        = None  , 
                               target_weight          = None  ) :
        """ Check if strong regularization is needed for BDT-based reweighting tests.
        Evaluates the limiting effective statistics between original and target samples,
        phase space dimensionality, and background subtraction weight fluctuations.
        """
        check_all ( data1   = original            ,
                    data2   = target              ,
                    weight1 = original_weight     ,
                    weight2 = target_weight       ,
                    where   = "%s:needs_regularizatino" % typename ( self ) ) 
        return RW_needs_regularization ( original        = original        ,
                                         target          = target          ,
                                         original_weight = original_weight ,
                                         target_weight   = target_weight   )
    
    # =========================================================================
    # Public Read-Only Properties
    # =========================================================================

    @property
    def original_ratios( self ) :
        """Computed density ratio factors r(x) for the original sample.
        Returns None if store_original_weights=False during initialization.
        """
        return self.__original_ratios

    @property
    def original_reweighted_weights( self ) :
        """Final reweighted event weights (w * r) for the original sample.
        Returns None if store_original_weights=False during initialization.
        """
        return self.__original_reweighted_weights

    @property
    def mode ( self ) :
        """Active stream decomposition scheme ('1-stream', '2-stream_orig', '2-stream_targ', or '4-stream')."""
        return self.__mode

    @property
    def config ( self ) :
        """`config` : Reweighter configuration"""
        conf = {}
        conf.update( super().config if hasattr( super(), "config" ) else {} )
        conf [ "mode" ]                        = self.__mode
        conf [ "clip_threshold" ]              = self.__clip_threshold
        conf [ "n_splits" ]                    = self.__n_splits
        conf [ "original_ratios" ]             = self.__original_ratios             is not None
        conf [ "original_reweighted_weights" ] = self.__original_reweighted_weights is not None
        return conf

    # =========================================================================
    # Abstract Methods (Must be implemented by child framework wrappers)
    # =========================================================================
    @abc.abstractmethod
    def _train_single_model( self    ,
                             X_train ,
                             y_train ,
                             w_train ,
                             X_val   ,
                             y_val   ,
                             w_val   ) :
        """ Train a single base classifier fold and return (model, val_predictions)."""
        raise NotImplementedError

    @abc.abstractmethod
    def _predict_single_model( self  ,
                               model ,
                               X     ) :
        """ Predict probability p(y=1|x) using a single trained model."""
        raise NotImplementedError
    
    # ========================================================================
    ## Checks if dataset size is small relative to feature count indicating a need for stronger regularization.
    def use_strong_regularization ( self , X ) :
        """ Checks if dataset size is small relative to feature count,
        indicating a need for stronger regularization.
        """
        ns = num_samples  ( X )
        nf = num_features ( X )
        return ns < max ( 300 if nf <= 3 else 1000 , nf * 50 )

    # =========================================================================
    ## get the regularized params 
    def regularized_config ( self , config , n_features , n_samples ) :
        return None 

    # =========================================================================
    # Internal Fitting Pipeline
    # =========================================================================
    def __fit_and_compute ( self            ,
                            original        ,
                            original_weight ,
                            target          ,
                            target_weight   ) :
        
        # =====================================================================
        ## FAST PATH: Standard positive unit weights for both samples
        # =====================================================================
        if original_weight is None and target_weight is None :
            original_weight = numpy.ones ( len ( original ), dtype = numpy.float32 )
            target_weight   = numpy.ones ( len ( target   ), dtype = numpy.float32 )
            
        if original_weight is None and target_weight is None :
            self.__mode  = "1-stream"
            ratios_valid = self.__fit_eval_stream ( "base", original, None, target, None )

            clipped            = numpy.clip    ( ratios_valid,
                                                 numpy.float32 ( 0.0 ) ,
                                                 numpy.float32 ( self.__clip_threshold ) , )
            self.__norm_factor = numpy.float32 ( len( original )
                                                 / ( numpy.sum ( clipped, dtype = numpy.float64 ) + 1e-12 ) )

            original_ratios = clipped * self.__norm_factor
            return original_ratios, original_ratios.copy()

        # =====================================================================
        ## SLOW PATH: Explicit handling of provided weight arrays
        # =====================================================================
        w_orig = original_weight if original_weight is not None else numpy.ones ( len ( original ) , dtype = numpy.float32 )
        w_targ = target_weight   if target_weight   is not None else numpy.ones ( len ( target   ) , dtype = numpy.float32 )

        # 1. Filter out zero-weight events
        nz_orig_mask = w_orig != 0
        nz_targ_mask = w_targ != 0

        X_orig_valid, w_orig_valid = original [ nz_orig_mask ] , w_orig [ nz_orig_mask ]
        X_targ_valid, w_targ_valid = target   [ nz_targ_mask ] , w_targ [ nz_targ_mask ]

        # 2. Strict sign masks
        mask_orig_pos, mask_orig_neg = w_orig_valid > 0, w_orig_valid < 0
        mask_targ_pos, mask_targ_neg = w_targ_valid > 0, w_targ_valid < 0

        has_neg_orig = numpy.any ( mask_orig_neg )
        has_neg_targ = numpy.any ( mask_targ_neg )

        # SCENARIO 1: Standard positive weights only
        if not has_neg_orig and not has_neg_targ :
            self.__mode  = "1-stream"
            ratios_valid = self.__fit_eval_stream(
                "base",
                X_orig_valid,
                w_orig_valid,
                X_targ_valid,
                w_targ_valid,
            )

        # SCENARIO 2: Negative weights ONLY in Original sample
        elif has_neg_orig and not has_neg_targ :
            self.__mode                   = "2-stream_orig"
            ratios_valid                  = numpy.zeros ( len( X_orig_valid ) , dtype = numpy.float32 )
            ratios_valid[ mask_orig_pos ] = self.__fit_eval_stream(
                "orig_pos",
                X_orig_valid[ mask_orig_pos ],
                w_orig_valid[ mask_orig_pos ],
                X_targ_valid,
                w_targ_valid,
            )
            ratios_valid[ mask_orig_neg ] = self.__fit_eval_stream(
                "orig_neg",
                X_orig_valid[ mask_orig_neg ],
                w_orig_valid[ mask_orig_neg ],
                X_targ_valid,
                w_targ_valid,
            )

        # SCENARIO 3: Negative weights ONLY in Target sample
        elif not has_neg_orig and has_neg_targ :
            self.__mode            = "2-stream_targ"
            X_targ_pos, w_targ_pos = X_targ_valid [ mask_targ_pos ] , w_targ_valid [ mask_targ_pos ] 
            X_targ_neg, w_targ_neg = X_targ_valid [ mask_targ_neg ] , w_targ_valid [ mask_targ_neg ]

            w_pos_sum                  = numpy.float32 ( numpy.sum ( w_targ_pos, dtype = numpy.float64 ) )
            w_neg_sum                  = numpy.float32 ( numpy.abs ( numpy.sum( w_targ_neg, dtype = numpy.float64 ) ) )
            
            self.__target_weights_info = { "W_pos": w_pos_sum, "W_neg": w_neg_sum }
            w_total                    = w_pos_sum + w_neg_sum

            r_pos_pos = self.__fit_eval_stream(
                "pos_pos",
                X_orig_valid,
                w_orig_valid,
                X_targ_pos,
                w_targ_pos,
            )
            r_pos_neg = self.__fit_eval_stream(
                "pos_neg",
                X_orig_valid,
                w_orig_valid,
                X_targ_neg,
                w_targ_neg,
            )

            ratios_valid = (
                r_pos_pos * w_pos_sum - r_pos_neg * w_neg_sum
            ) / ( w_total + numpy.float32( 1e-12 ) )

        # SCENARIO 4: Fully symmetric 4-stream (Negative weights in both)
        else :
            self.__mode            = "4-stream"
            X_orig_pos , w_orig_pos = X_orig_valid [ mask_orig_pos ] , w_orig_valid [ mask_orig_pos ]
            X_orig_neg , w_orig_neg = X_orig_valid [ mask_orig_neg ] , w_orig_valid [ mask_orig_neg ]
            X_targ_pos , w_targ_pos = X_targ_valid [ mask_targ_pos ] , w_targ_valid [ mask_targ_pos ]
            X_targ_neg , w_targ_neg = X_targ_valid [ mask_targ_neg ] , w_targ_valid [ mask_targ_neg ]

            w_pos_sum = numpy.float32 ( numpy.sum ( w_targ_pos, dtype = numpy.float64 ) )
            w_neg_sum = numpy.float32 ( numpy.abs ( numpy.sum ( w_targ_neg, dtype = numpy.float64 ) ) ) 
            self.__target_weights_info = { "W_pos": w_pos_sum, "W_neg": w_neg_sum }
            w_total                    = w_pos_sum + w_neg_sum

            r_pos_pos = self.__fit_eval_stream(
                "pos_pos", X_orig_pos, w_orig_pos, X_targ_pos, w_targ_pos
            )
            r_pos_neg = self.__fit_eval_stream(
                "pos_neg", X_orig_pos, w_orig_pos, X_targ_neg, w_targ_neg
            )
            r_neg_pos = self.__fit_eval_stream(
                "neg_pos", X_orig_neg, w_orig_neg, X_targ_pos, w_targ_pos
            )
            r_neg_neg = self.__fit_eval_stream(
                "neg_neg", X_orig_neg, w_orig_neg, X_targ_neg, w_targ_neg
            )

            ratios_valid = numpy.zeros ( len ( X_orig_valid ) , dtype = numpy.float32 )
            denom                         = w_total + numpy.float32( 1e-12 )
            ratios_valid[ mask_orig_pos ] = (
                r_pos_pos * w_pos_sum - r_pos_neg * w_neg_sum
            ) / denom
            ratios_valid[ mask_orig_neg ] = (
                r_neg_pos * w_pos_sum - r_neg_neg * w_neg_sum
            ) / denom

        # 3. Clip, normalize, and construct final arrays
        ratios_valid_norm = self.__normalize_and_clip ( ratios_valid, w_orig_valid )

        original_ratios                 = numpy.zeros( len( w_orig ), dtype = numpy.float32 )
        original_ratios[ nz_orig_mask ] = ratios_valid_norm
        original_reweighted_weights     = w_orig * original_ratios

        return original_ratios, original_reweighted_weights




    def __fit_eval_stream( self       ,
                           stream_key ,
                           X_orig_sub ,
                           w_orig_sub ,
                           X_targ_sub ,
                           w_targ_sub ) :

        X_comb = numpy.vstack( [ X_orig_sub, X_targ_sub ] )
        y_comb = numpy.hstack( [ numpy.zeros ( len( X_orig_sub ) , dtype = numpy.float32 ) ,
                                 numpy.ones  ( len( X_targ_sub ) , dtype = numpy.float32 ) ] )

        if w_orig_sub is None and w_targ_sub is None :
            w_comb   = None
            scale_factor = numpy.float32 ( len ( X_targ_sub ) / len ( X_orig_sub ) )
        else :
            w_orig_abs = numpy.abs ( w_orig_sub ) if w_orig_sub is not None else numpy.ones ( len ( X_orig_sub ), dtype = numpy.float32 )
            w_targ_abs = numpy.abs ( w_targ_sub ) if w_targ_sub is not None else numpy.ones ( len ( X_targ_sub ), dtype = numpy.float32 )
            
            sum_w_orig = numpy.sum ( w_orig_abs )
            sum_w_targ = numpy.sum ( w_targ_abs )
            
            scale_factor = ( sum_w_targ / sum_w_orig ) if sum_w_orig > 0 else numpy.float32 ( 1.0 )
            
            w_orig_scaled = w_orig_abs * scale_factor
            w_comb        = numpy.hstack ( [ w_orig_scaled , w_targ_abs ] )

        oof_raw = numpy.zeros( len( X_comb ), dtype = numpy.float32 )

        from sklearn.model_selection import StratifiedKFold
        from sklearn.linear_model import LogisticRegression

        skf = StratifiedKFold ( n_splits     = self.n_splits     ,
                                shuffle      = True              ,
                                random_state = self.random_state )

        stream_models = []
        for train_idx, val_idx in skf.split( X_comb, y_comb ) :
            X_tr, y_tr = X_comb[ train_idx ], y_comb[ train_idx ]
            X_va, y_va = X_comb[ val_idx ], y_comb[ val_idx ]

            w_tr = w_comb[ train_idx ] if w_comb is not None else None
            w_va = w_comb[ val_idx ] if w_comb is not None else None

            # 1. Train base booster model
            model, _ = self._train_single_model ( X_tr, y_tr, w_tr, X_va, y_va, w_va )
            stream_models.append( model )
            
            # 2. Get raw probabilities via subclass polymorphism
            raw_val_p = self._predict_single_model ( model, X_va )
            if raw_val_p.ndim > 1 :
                raw_val_p = raw_val_p[:, 1]
                
            oof_raw[ val_idx ] = raw_val_p.astype ( numpy.float32, copy = False )

        # 3. Fit global Platt scaling on OOF predictions
        eps = 1e-6
        oof_raw_clipped = numpy.clip ( oof_raw, eps, 1.0 - eps )
        oof_logits = numpy.log ( oof_raw_clipped / ( 1.0 - oof_raw_clipped ) )

        platt = LogisticRegression()
        platt.fit ( oof_logits.reshape(-1, 1), y_comb, sample_weight = w_comb )

        # Store models and calibrator tuple
        self.__fitted_models[ stream_key ] = ( stream_models, platt )
        
        self.__scale_factors[ stream_key ] = scale_factor

        self.__priors[ stream_key ] = numpy.float32( 1.0 )

        # Compute OOF calibrated predictions for original sample
        oof_calib = platt.predict_proba ( oof_logits.reshape(-1, 1) )[:, 1]
        p_orig = oof_calib[ : len ( X_orig_sub ) ]
        p_orig_clipped = numpy.clip ( p_orig, numpy.float32 ( 1e-4 ) , numpy.float32( 1.0 - 1e-4 ) )
        
        ratios = p_orig_clipped / ( numpy.float32( 1.0 ) - p_orig_clipped )
        return ratios / scale_factor

    def __normalize_and_clip( self   ,
                              ratios ,
                              w_orig ) :
        
        clipped_ratios = numpy.clip ( ratios, numpy.float32 ( 0.0 ) , numpy.float32 ( self.__clip_threshold ) )
        orig_sum       = numpy.sum  ( w_orig, dtype = numpy.float64 )
        reweighted_sum = numpy.sum  ( w_orig * clipped_ratios, dtype = numpy.float64 )
        self.__norm_factor = numpy.float32 ( orig_sum / ( reweighted_sum + 1e-12 ) )
        return clipped_ratios * self.__norm_factor
    

    def __predict_stream_ratios( self, stream_key, X ):
        models, platt = self.__fitted_models[ stream_key ]
        scale_factor  = self.__scale_factors[ stream_key ]

        ratios_list = []
        for model in models :
            p = self._predict_single_model ( model, X )
            if p.ndim > 1 :
                p = p[:, 1]
                
            eps = 1e-4
            p_clipped = numpy.clip ( p, eps, 1.0 - eps )
            logits = numpy.log ( p_clipped / ( 1.0 - p_clipped ) )
            
            p_calib = platt.predict_proba ( logits.reshape(-1, 1) )[:, 1]
            p_calib_clipped = numpy.clip ( p_calib, numpy.float32 ( eps ) , numpy.float32 ( 1.0 - eps ) )
            
            stream_ratios = ( p_calib_clipped / ( numpy.float32( 1.0 ) - p_calib_clipped ) ) / scale_factor
            ratios_list.append( stream_ratios )
            
        return numpy.mean ( ratios_list, axis = 0, dtype = numpy.float32 )

    
    # =========================================================================
    # Public Inference Method
    # =========================================================================

    def weights ( self                   ,
                  original               ,
                  original_weight = None ) :
        """ Apply the fitted immutable model ensemble to evaluate reweighted events for a new original sample.

        Returns
        -------
        numpy.ndarray
            Calculated event weights (factors * original_weight).
        """

        X_new_f32 = original.astype ( numpy.float32, copy = False )

        if not valid_data_shape   ( X_new_f32                   ) : raise TypeError ( "Invalid `original` type/shape: %s" % typename ( original ) )
        if not valid_weight       ( original_weight             ) : raise TypeError ( "Invalid `original_weight`!" )        
        if not compatible_weights ( X_new_f32 , original_weight ) : raise TypeError ( "Incompatible `original` data/weight!" )
        if self.n_features != num_features ( X_new_f32 )          : raise TypeError ( "Invalid #features!!")
        

        # FAST PATH: Single-stream without input weights
        if original_weight is None and self.__mode == "1-stream" :
            ratios  = self.__predict_stream_ratios( "base", X_new_f32 )
            clipped = numpy.clip ( ratios,
                                   numpy.float32 ( 0.0 ) ,
                                   numpy.float32 ( self.__clip_threshold ) )
            return clipped * self.__norm_factor

        # SLOW PATH: Multi-stream or explicit input weights
        w_new_f32 = ( original_weight.astype( numpy.float32, copy = False )
                      if original_weight is not None
                      else numpy.ones( len( original ), dtype = numpy.float32 ) )

        nz_mask          = w_new_f32 != 0
        X_valid, w_valid = X_new_f32[ nz_mask ], w_new_f32[ nz_mask ]

        mask_pos     = w_valid > 0
        mask_neg     = w_valid < 0
        ratios_valid = numpy.zeros( len( X_valid ), dtype = numpy.float32 )

        if self.__mode == "1-stream" :
            ratios_valid = self.__predict_stream_ratios( "base", X_valid )

        elif self.__mode == "2-stream_orig" :
            if numpy.any( mask_pos ) :
                ratios_valid[ mask_pos ] = self.__predict_stream_ratios(
                    "orig_pos", X_valid[ mask_pos ]
                )
            if numpy.any( mask_neg ) :
                ratios_valid[ mask_neg ] = self.__predict_stream_ratios(
                    "orig_neg", X_valid[ mask_neg ]
                )

        elif self.__mode == "2-stream_targ" :
            w_pos   = self.__target_weights_info[ "W_pos" ]
            w_neg   = self.__target_weights_info[ "W_neg" ]
            w_total = w_pos + w_neg + numpy.float32( 1e-12 )

            r_pos_pos    = self.__predict_stream_ratios( "pos_pos", X_valid )
            r_pos_neg    = self.__predict_stream_ratios( "pos_neg", X_valid )
            ratios_valid = (
                r_pos_pos * w_pos - r_pos_neg * w_neg
            ) / w_total

        elif self.__mode == "4-stream" :
            w_pos   = self.__target_weights_info[ "W_pos" ]
            w_neg   = self.__target_weights_info[ "W_neg" ]
            w_total = w_pos + w_neg + numpy.float32( 1e-12 )

            if numpy.any( mask_pos ) :
                X_p                      = X_valid[ mask_pos ]
                r_pos_pos                = self.__predict_stream_ratios( "pos_pos", X_p )
                r_pos_neg                = self.__predict_stream_ratios( "pos_neg", X_p )
                ratios_valid[ mask_pos ] = (
                    r_pos_pos * w_pos - r_pos_neg * w_neg
                ) / w_total

            if numpy.any( mask_neg ) :
                X_n                      = X_valid[ mask_neg ]
                r_neg_pos                = self.__predict_stream_ratios( "neg_pos", X_n )
                r_neg_neg                = self.__predict_stream_ratios( "neg_neg", X_n )
                ratios_valid[ mask_neg ] = (
                    r_neg_pos * w_pos - r_neg_neg * w_neg
                ) / w_total

        clipped                = numpy.clip ( ratios_valid,
                                              numpy.float32 ( 0.0 ),
                                              numpy.float32 ( self.__clip_threshold ) )
        final_ratios           = numpy.zeros( len( w_new_f32 ), dtype = numpy.float32 )
        final_ratios[ nz_mask ] = clipped * self.__norm_factor

        return final_ratios if original_weight is None else final_ratios * w_new_f32

# =================================================================
## @class LightGBMDensityReweighter
#   Density-ratio reweighter implementation using LightGBM as the underlying classifier.
#   Inherits from BaseDensityReweighter. Automatically handles 1-, 2-, and 4-stream
#   decompositions for positive and negative sample weights.
class LightGBMDensityReweighter ( DensityReweighter ) :
    """
    LightGBM-based density ratio reweighter for Ostap
    """
    def __init__ ( self            , 
                   original        , 
                   target          , 
                   original_weight = None , 
                   target_weight   = None , 
                   **kwargs        ) :
        
        # --- Default High-Capacity configuration for large statistics
        config = {
            'objective'             : 'binary'            ,
            'metric'                : 'binary_logloss'    ,
            'n_estimators'          : DEFAULT_ESTIMATORS  , ## Default 400 trees budget
            'learning_rate'         : 0.03                , ## Smooth updates for KDE-like density ratio
            'max_depth'             : 4                   , ## Default depth for rich phase space
            'num_leaves'            : 15                  ,
            'min_child_samples'     : 30                  ,
            'min_child_weight'      : 1e-3                ,
            'reg_alpha'             : 0.1                 ,
            'reg_lambda'            : 2.0                 , ## Moderate L2 penalty on leaf values
            'subsample'             : 0.8                 ,
            'subsample_freq'        : 1                   ,
            'colsample_bytree'      : 0.8                 ,
            'path_smooth'           : 1.0                 , ## Node smoothing suppresses leaf-boundary spikiness
            'boost_from_average'    : True                ,
            'early_stopping_rounds' : None                ,
            'verbosity'             : -1                  ,
            'n_jobs'                : -1                  ,
        }
        config.update ( kwargs )
        
        super ( LightGBMDensityReweighter , self ).__init__ (
            original        = original        ,
            target          = target          ,
            original_weight = original_weight ,
            target_weight   = target_weight   ,
            **config
        )

    @property
    def method ( self ) :
        """ Return the reweighter method name
        """
        return 'Density Reweighter/LightGBM'

    def regularization ( self                ,
                         params     , 
                         n_features ,
                         n_samples  ) :
        """ Match hep_ml.GBReweighter's flexibility by avoiding aggressive leaf constraints 
        and allowing the model to reach extreme probabilities (high weights).
        """


        params [ 'n_estimators'      ] = min ( REGULARIZED_ESTIMATORS , params.get ( 'n_estimators' , REGULARIZED_ESTIMATORS ) )
        params [ 'learning_rate'     ] = 0.1     
        params [ 'max_depth'         ] = 5       
        params [ 'num_leaves'        ] = 31      
        
        # Soft dynamic limits for leaves to allow deep splits in rare tails
        leaf_samples = max ( 2, min ( 30, int ( n_samples  * 0.0001 ) ) )
        params [ 'min_child_samples' ] = leaf_samples      
        params [ 'min_child_weight'  ] = max ( 1e-4, float ( n_samples * 0.00005 ) )  
        
        params [ 'reg_alpha'         ] = 0.0    
        params [ 'reg_lambda'        ] = 0.0 
        params [ 'subsample'         ] = 1.0    
        params [ 'colsample_bytree'  ] = 1.0  
        params [ 'early_stopping_rounds' ] = None
        
        if 'path_smooth' in params : 
            params.pop ( 'path_smooth' )

        return params
    
    def _train_single_model ( self    ,
                              X_train , y_train , w_train ,
                              X_val   , y_val   , w_val   ) : 
        """
        Train single LightGBM estimator and return validation density ratio weights
        """
        import lightgbm as LightGBM 

        w_tr = numpy.ascontiguousarray ( w_train , dtype = numpy.float32 ) if w_train is not None else None
        w_va = numpy.ascontiguousarray ( w_val   , dtype = numpy.float32 ) if w_val   is not None else None

        trn_data = LightGBM.Dataset ( X_train , label = y_train , weight = w_tr , free_raw_data = False )
        val_data = LightGBM.Dataset ( X_val   , label = y_val   , weight = w_va , reference = trn_data , free_raw_data = False )

        params = self.params.copy ()
        
        num_boost_round       = params.pop ( 'num_boost_round' , None ) or params.pop ( 'n_estimators' , 400 )
        
        # --- Extract early_stopping_rounds, defaulting to None for unbiased OOF estimation
        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , None )

        callbacks = []
        if early_stopping_rounds is not None and eval_set is not None :
            callbacks.append ( lightgbm.early_stopping ( stopping_rounds = early_stopping_rounds , verbose = False ) )

        model = LightGBM.train ( params          = params          ,
                                 train_set       = trn_data        ,
                                 num_boost_round = num_boost_round ,
                                 valid_sets      = [ val_data ]    ,
                                 callbacks       = callbacks       )
        
        val_preds = self._predict_single_model ( model , X_val )
        
        return model , val_preds


    def _predict_single_model( self, model, X ):
        best_iter = getattr( model, 'best_iteration', 0 )
        kwargs = {}
        if best_iter is not None and best_iter > 0 :
            kwargs[ 'num_iteration' ] = best_iter            
        p = model.predict( X, **kwargs )
        return p.astype( numpy.float32, copy = False )
    
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
        
        # --- Default High-Capacity configuration for large statistics
        config = {
            'objective'             : 'binary:logistic'   ,
            'eval_metric'           : 'logloss'           ,
            'n_estimators'          : DEFAULT_ESTIMATORS  , ## Default 400 trees budget
            'learning_rate'         : 0.03                , ## Smooth updates for KDE-like density ratio
            'max_depth'             : 4                   , ## Default depth for rich phase space
            'min_child_weight'      : 5.0                 , ## Minimum sum of hessians per leaf
            'gamma'                 : 0.05                , ## Minimum loss reduction to force split
            'reg_alpha'             : 0.1                 ,
            'reg_lambda'            : 2.0                 , ## Moderate L2 penalty on leaf values
            'subsample'             : 0.8                 ,
            'colsample_bytree'      : 0.8                 ,
            'tree_method'           : 'hist'              , ## Fast histogram-based binning
            'early_stopping_rounds' : None                ,
            'verbosity'             : 0                   ,
            'n_jobs'                : -1                  ,
        }
        
        config.update ( params )
        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' )

        # Delegate execution to BaseDensityReweighter __init__
        super().__init__ ( original               = original               ,
                           target                 = target                 ,
                           original_weight        = original_weight        ,
                           target_weight          = target_weight          ,
                           store_original_weights = store_original_weights , **config )
        
    @property
    def method ( self ) :
        """`method` : underlying method/engine"""
        return "Density Reweighter/XGBoost"
    
    def regularization ( self , params , n_features , n_samples ) :
        """ Dynamic regularization tuned for XGBoost to match LightGBM performance."""
        
        params [ 'n_estimators'          ] = min ( REGULARIZED_ESTIMATORS , params.get ( 'n_estimators' , REGULARIZED_ESTIMATORS ) )
        params [ 'learning_rate'         ] = 0.1
        params [ 'max_depth'             ] = 5
        
        # Soft dynamic min_child_weight scaling with neff, keeping a low lower bound
        params [ 'min_child_weight'      ] = max ( 1e-4, float ( n_samples * 0.0001 ) )
        
        params [ 'alpha'                 ] = 0.0
        params [ 'lambda'                ] = 0.0
        params [ 'subsample'             ] = 1.0
        params [ 'colsample_bytree'      ] = 1.0
        params [ 'early_stopping_rounds' ] = None
    
        return params

    def _predict_single_model( self, model, X ):
        import xgboost as XGBoost
        return model.predict( XGBoost.DMatrix ( X ) )
    
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

        # --- Construct XGBoost DMatrix instances
        dtrain = XGBoost.DMatrix ( X_train , label = y_train , weight = w_train )
        dval   = XGBoost.DMatrix ( X_val   , label = y_val   , weight = w_val   )

        params = {}
        params.update ( self.params )

        # --- Determine number of boosting rounds
        num_boost_round = params.pop ( 'num_boost_round' , None ) or params.pop ( 'n_estimators' , None ) or 500
        if isinstance ( num_boost_round , int ) and 10 < num_boost_round < 10000 : pass
        else                                                                    : num_boost_round = 500

        # --- Extract early_stopping_rounds, defaulting to None for unbiased OOF estimation
        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , None )
        if isinstance ( early_stopping_rounds , int ) and 1 < early_stopping_rounds < num_boost_round : pass
        else                                                                                          : early_stopping_rounds = None

        train_kwargs = {}
        if early_stopping_rounds is not None :
            train_kwargs [ 'evals' ]                 = [ ( dval , "val" ) ]
            train_kwargs [ 'early_stopping_rounds' ] = early_stopping_rounds

        # --- Train XGBoost booster
        model = XGBoost.train ( params          = params,
                                dtrain          = dtrain,
                                num_boost_round = num_boost_round,
                                verbose_eval    = False,
                                **train_kwargs )

        # --- Predict probabilities on validation fold
        predict_kwargs = {}
        best_iter = getattr ( model , "best_iteration" , None )
        if early_stopping_rounds is not None and best_iter is not None and 0 < best_iter :
            predict_kwargs [ 'iteration_range' ] = ( 0 , best_iter + 1 )

        val_preds = model.predict ( dval , **predict_kwargs )
        return model , val_preds.astype ( numpy.float32 , copy = False )


    def _predict_single_model( self, model, X ):
        
        best_iter = getattr ( model, 'best_iteration', 0 )
        kwargs = {}
        if best_iter is not None and 0 < best_iter :
            kwargs [ 'ntree_limit' ] = best_iter # or num_iteration
            
        import xgboost as XGBoost
        dmat = XGBoost.DMatrix ( X    )
        p    = model.predict   ( dmat , **kwargs )
        return p.astype ( numpy.float32, copy = False )
        
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


        # --- Default High-Capacity configuration for large statistics
        config = {
            'loss_function'         : 'Logloss'           ,
            'eval_metric'           : 'Logloss'           ,
            'n_estimators'          : DEFAULT_ESTIMATORS  , ## Default 400 trees budget
            'learning_rate'         : 0.03                , ## Smooth updates for KDE-like density ratio
            'depth'                 : 4                   , ## Default depth for rich phase space
            'l2_leaf_reg'           : 2.0                 , ## Moderate L2 penalty on leaf values
            'min_child_samples'     : 30                  ,
            'subsample'             : 0.8                 ,
            'random_strength'       : 1.0                 , ## Adds randomness to split scoring
            'early_stopping_rounds' : None                ,
            'verbose'               : False               ,
            'thread_count'          : -1                  ,
        }
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
        
    @property
    def method ( self ) :
        """`method` : underlying method/engine"""
        return "Density Reweighter/CatBoost"


    def regularization ( self , params , n_features , n_samples ) :
        """ Dynamic regularization tuned for CatBoost to allow high density ratios."""
        
        params [ 'n_estimators'         ] = min ( REGULARIZED_ESTIMATORS , params.get ( 'n_estimators' , REGULARIZED_ESTIMATORS ) )
        params [ 'learning_rate'        ] = 0.1
        params [ 'depth'                ] = 6
        
        if 'min_data_in_leaf' in params :  params.pop ( 'min_data_in_leaf' , None )
        
        # Soft dynamic min_child_samples scaling with neff
        params [ 'min_child_samples'    ] = max ( 2, min ( 30, int ( n_samples * 0.0001 ) ) )
        
        params [ 'l2_leaf_reg'           ] = 1.0
        params [ 'subsample'             ] = 1.0
        params [ 'early_stopping_rounds' ] = None
        
        return params

    # =============================================================
    # Abstract Method Implementations
    # =============================================================
    def _train_single_model ( self    ,
                              X_train , y_train , w_train ,
                              X_val   , y_val   , w_val   ) :
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

        # --- Determine number of iterations
        iterations = params.pop ( 'iterations' , None ) or params.pop ( 'n_estimators' , None ) or 500
        if isinstance ( iterations , int ) and 10 < iterations < 10000 : pass
        else                                                           : iterations = 500

        # --- Extract early_stopping_rounds, defaulting to None for unbiased OOF estimation
        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , None )
        if isinstance ( early_stopping_rounds , int ) and 1 < early_stopping_rounds < iterations : pass
        else                                                                                      : early_stopping_rounds = None

        params [ 'iterations' ] = iterations

        fit_kwargs = {}
        if early_stopping_rounds is not None :
            params [ 'early_stopping_rounds' ] = early_stopping_rounds
            params [ 'use_best_model'        ] = True
            fit_kwargs [ 'eval_set'          ] = val_pool

        model = CatBoost.CatBoostClassifier ( **params )

        model.fit ( trn_pool ,
                    verbose = False ,
                    **fit_kwargs )

        # --- Predict probability for positive class p(y=1|x)
        val_preds = model.predict_proba ( val_pool ) [ : , 1 ]
        return model , val_preds.astype ( numpy.float32 , copy = False )


    def _predict_single_model( self, model, X ):
        best_iter = getattr( model, 'best_iteration', 0 )
        kwargs = {}
        if best_iter is not None and best_iter > 0 :
            kwargs[ 'num_iteration' ] = best_iter            
        p = model.predict( X, **kwargs )
        return p.astype( numpy.float32, copy = False )

# ==============================================================================
## @class GBReweighter
#  Helper class for reweighting using <code>GBReweighter</code> from hep_ml by Alex Rogozhnikov 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-09-22 
class GBReweighter(Reweighter) :
    """ Helper class for reweighting using `GBReweighter` from hep_ml by Alex Rogozhnikov
    - see hep_ml.reweight.GBReweighter
        """
    def __init__ ( self                   , * , 
                   original               ,
                   target                 ,
                   original_weight        = None  ,
                   target_weight          = None  , 
                   silent                 = False ,
                   n_splits               = 5     ,
                   store_original_weights = True  , **params ) :
        """ Helper class for reweighting using `GBReweighter`
        - see hep_ml.reweight.GBReweighter
        """
        
        if not isinstance ( n_splits  , int ) : raise TypeError  ( "Invalid `n_splits' type  %s"       % typename ( n_splits ) )
        if not 0 <= n_splits <= 1000          : raise ValueError ( "Invalid `n_splits' value %s"       % n_splits )     

        self.__n_splits = n_splits
        
        # Baseline GBReweighter configuration (already achieving good p-value)
        config = {            
            "n_estimators"      : 150    , 
            "learning_rate"     :   0.03 ,
            "max_depth"         :   5    ,
            "min_samples_leaf"  :  30    ,            
            # Advanced Scikit-Learn GradientBoosting/GBReweighter arguments
            "gb_args"           : {
                "subsample"     : 0.8    ,
                "max_features"  : 0.65   }
        }
        
        config.update ( params ) 
        
        if self.needs_regularization ( original        = original        ,
                                       target          = target          ,
                                       original_weight = original_weight ,
                                       target_weight   = target_weight   ) :
            
            n_features = num_features ( original )
            
            neff_orig  = nEff ( original  , original_weight )
            neff_targ  = nEff ( target    , target_weight   ) 
            neff       = min  ( neff_orig , neff_targ       )
            
            config.update ( self.regularization ( config , n_features , neff ) )
            logger.attention ( "%s: strong regularization is applied" % typename ( self ) ) 
            
        # =====================================================================
        ## Initialize the base: check input data & print config 
        # =====================================================================
        Reweighter.__init__ ( self            ,
                              original        = original        ,
                              target          = target          , 
                              original_weight = original_weight ,
                              target_weight   = target_weight   , **config )

        # =====================================================================
        ## some patch: needed? 
        try : # ===============================================================
            # =================================================================
            import sklearn.tree._classes as _cl
            if hasattr ( _cl , 'CRITERIA_REG' ) :
                if 'mse' in _cl.CRITERIA_REG and 'squared_error' not in _cl.CRITERIA_REG:
                    _cl.CRITERIA_REG [ 'squared_error' ] = _cl.CRITERIA_REG [ 'mse']
                elif 'squared_error' in _cl.CRITERIA_REG and 'mse' not in _cl.CRITERIA_REG:
                    _cl.CRITERIA_REG [ 'mse'] = _cl.CRITERIA_REG [ 'squared_error' ]
            # =================================================================
        except ( ImportError , AttributeError ) : # ============================
            # =================================================================
            pass

        if not hasattr ( numpy , 'float' ) :
            logger.warning ( 'No `numpy.float` ... add it!')
            numpy.float = numpy.float64
            
        # =====================================================================
        ## Execute training immediately and compute original sample weights (float32) 
        # =====================================================================

        from hep_ml.reweight import GBReweighter as GBRW
        self.__reweighter = GBRW ( **self.params  )
        
        if 0 < self.n_splits :
            from hep_ml.reweight import FoldingReweighter as FRW
            self.__reweighter = FRW ( self.reweighter ,
                                      n_folds         = self.n_splits     , 
                                      random_state    = self.random_state )
               
        with logAttention() if self.silent else NoContext() : 
            self.__reweighter.fit ( original        ,
                                    target          ,
                                    original_weight = original_weight , 
                                    target_weight   = target_weight   )


        self.__original_ratios             = None
        self.__original_reweighted_weights = None
        if store_original_weights :
            factors = self.__reweighter.predict_weights (
                original, 
                original_weight = original_weight )
            self.__original_ratios             = factors 
            self.__original_reweighted_weights = factors if weight_trivial ( original_weight ) else factors * original_weight

    # =============================================================================
    ## Check if strong regularization is needed considering BOTH original and target samples.
    #  Evaluates sample dimensionality, Kish's effective sample size (nEff),
    #  and the fraction of negative sPlot weights.
    def needs_regularization ( self ,
                               original                       ,
                               target                         ,
                               original_weight        = None  , 
                               target_weight          = None  ) :
        """ Check if strong regularization is needed for BDT-based reweighting tests.
        Evaluates the limiting effective statistics between original and target samples,
        phase space dimensionality, and background subtraction weight fluctuations.
        """
        check_all ( data1   = original            ,
                    data2   = target              ,
                    weight1 = original_weight     ,
                    weight2 = target_weight       ,
                    where   = "%s:needs_regularizatino" % typename ( self ) ) 
        return RW_needs_regularization ( original        = original        ,
                                         target          = target          ,
                                         original_weight = original_weight ,
                                         target_weight   = target_weight   )
    
    def regularization ( self       ,
                         params     , 
                         n_features ,
                         n_samples  ) : 
        N = n_samples
        params [ 'n_estimators'     ] = 100
        params [ 'learning_rate'    ] = 0.05
        params [ 'max_depth'        ] = 3
        params [ 'min_samples_leaf' ] = max ( 100, int ( N * 0.005 ) )
        
        gb_args = params.get ( 'gb_args', {} ).copy ()
        gb_args.update ( { 'subsample': 0.8, 'max_features': 0.8 } )
        params [ 'gb_args'          ] = gb_args
        return params
    
    @property
    def method ( self ) :
        """`method` : underlying method/engine"""
        return "HepML/GBReweighter"
                        
    @property
    def n_splits ( self ) :
        """`n_splits` : number of splits/folds, same as `n_folds`"""
        return self.__n_splits
    
    @property
    def n_folds ( self ) :
        """`n_splits` : number of splits/folds, same as `n_splits`"""
        return self.n_splits 
    
    @property 
    def config ( self ) :
        """`config` : Reweighter configuraton"""
        conf = {}
        conf.update ( super().config )
        conf [ 'n_splits'  ] = self.n_splits
        return conf
    
    @property
    def original_ratios ( self ) :
        """ Computed ratio factors r(x) for the original sample.
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
    def reweighter ( self ) :
        """`reweighter' : get the underlying reweighter object"""
        return self.__reweighter
    
    # =========================================================================
    ## Get/predict new weights for (new) original
    def weights ( self                   ,
                  original               ,
                  original_weight = None ) :
        """ Get/predict  new weights for (new) original
        """
        ## 
        if not valid_data_shape   ( original                   ) : raise TypeError ( "Invalid `original` type/shape: %s" % typename ( original ) )
        if not valid_weight       ( original_weight            ) : raise TypeError ( "Invalid `original_weight`!" )        
        if not compatible_weights ( original , original_weight ) : raise TypeError ( "Incompatible `original` data/weight!" )
        if self.n_features != num_features ( original )          : raise TypeError ( "Invalid #features!!")
        ##

        
        factors = self.reweighter.predict_weights ( original        = original        ,
                                                    original_weight = original_weight )
        ## 
        return factors if weight_trivial ( original_weight ) else factors * original_weight
    
        
# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    from ostap.stats.tools import ( hasLightGBM ,
                                    hasXGBoost  ,
                                    hasCatBoost ,
                                    hasHepML    )

    if not hasLightGBM ( False ) : logger.warning  ( "No LightGBM available!" ) 
    if not hasXGBoost  ( False ) : logger.warning  ( "No XGBoost  available!" ) 
    if not hasCatBoost ( False ) : logger.warning  ( "No CatBoost available!" ) 
    if not hasHepML    ( False ) : logger.warning  ( "No HepMC    available!" ) 

    
# =============================================================================
##                                                                      The END 
# =============================================================================
