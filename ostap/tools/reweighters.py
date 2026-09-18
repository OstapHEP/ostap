#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/reweghters.py
#  Density ratio reweighters based on LightGBM, XGBoost, CatBoost
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @Date  2026-07-18
# =============================================================================
""" Density ratio reweighters based on LightGBM, XGBoost, CatBoost, and HepML.
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    'LightGBMDensityReweighter' , # LightGBM-based density reweighter
    'XGBoostDensityReweighter'  , # XGBoost-based density reweighter
    'CatBoostDensityReweighter' , # CatBoost-based density reweighter
    'PyTorchDensityReweighter'  , # PyTorch-based density reweighter 
    'GBReweighter'              , # Reweighter based on GBReweighter from hep_ml 
) 
# =============================================================================
from   ostap.core.ostap_types   import num_types 
from   ostap.utils.core         import typename
from   ostap.utils.basic        import numcpu, num_jobs, NoContext
from   ostap.logger.utils       import map2table_ex
from   ostap.logger.pretty      import nice_print 
from   ostap.logger.symbols     import arrow_right  
from   ostap.utils.progress_bar import progress_bar, ProgressBar  
from   ostap.tools.reweighter   import Reweighter
from   ostap.stats.counters     import SE, table_counters  
from   ostap.stats.utils        import ( weight_trivial     ,
                                         valid_weight       ,
                                         valid_data_shape   ,
                                         compatible_weights , 
                                         num_features       ,
                                         num_samples        ,
                                         check_all          ,
                                         nEff               )
import ostap.logger.symbols   as     S
import numpy, abc, copy, warnings
# =============================================================================
# Logging setup
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighters' )
else                       : logger = getLogger( __name__ )
# =============================================================================
# Global Configuration Constants
# =============================================================================
DEFAULT_ESTIMATORS     = 500
REGULARIZED_ESTIMATORS = 250
MAX_DEPTH              =   5 
REG_DEPTH              =   3 
# =============================================================================
method_LGBM  = 'DRW/%s'   % ( S.light_bulb           if S.light_bulb  else 'LightGBM' ) 
method_XGB   = 'DRW/%s'   % ( S.rocket               if S.rocket      else 'XGBoost'  ) 
method_CATB  = 'DRW/%s'   % ( S.cat_face             if S.cat_face    else 'CatBoost' ) 
method_TORCH = 'DRW/%s'   % ( S.flashlight           if S.flashlight  else 'TORCH'    ) 
method_GBRW  = 'HepML/%s' % ( S.wood                 if S.wood        else 'GBRW'     ) 
# =============================================================================
## @brief Check if strong regularization is needed for BDT-based reweighting.
#  @param original Features array for the original sample.
#  @param target Features array for the target sample.
#  @param original_weight Event weights for the original sample (optional).
#  @param target_weight Event weights for the target sample (optional).
#  @return Explanation string if needed, else empty string.
def RW_needs_regularization ( original                       ,
                              target                         ,
                              original_weight        = None  , 
                              target_weight          = None  ) :
    """Check if strong regularization is needed for BDT-based reweighting.

    Evaluates phase space dimensionality, Kish's effective sample statistics (nEff),
    and weight efficiency ratios to decide whether tight constraints should be enforced.

    :param original: Features array for original sample.
    :param target: Features array for target sample.
    :param original_weight: Event weights for original sample.
    :param target_weight: Event weights for target sample.
    :return: String with reason for regularization if needed; empty string otherwise.
    """
    nf = num_features ( original )
    
    # 1. Raw sample sizes
    nraw_orig = num_samples ( original )
    nraw_targ = num_samples ( target   )
    
    # 2. Calculate effective statistics (nEff)
    neff_orig = nEff ( original  , original_weight )
    neff_targ = nEff ( target    , target_weight   ) 
    
    # Two-sample effective size: Harmonic pooled nEff
    if neff_orig <= 0 or neff_targ <= 0 : neff = 0.0
    else : neff = 4.0 * ( neff_orig * neff_targ ) / ( neff_orig + neff_targ )

    # 3. Check weight efficiency ratios (nEff / nRaw)
    eff_orig = neff_orig / float ( nraw_orig ) if 0 < nraw_orig else 0.0
    eff_targ = neff_targ / float ( nraw_targ ) if 0 < nraw_targ else 0.0
    
    if eff_orig < 0.65 : return 'eff_orig<65%'
    if eff_targ < 0.65 : return 'eff_targ<65%'

    # 4. Low dimensionality (<= 4 features) needs regularization under limited statistics
    threshold = 50000.0 
    if nf <= 4 and neff < threshold :
        th = nice_print ( threshold ) 
        return 'nf<=4&neff_pooled<%s' % th 
    
    # 5. Non-linear density threshold for multidimensional phase space growth
    required_stats = 1500.0 * ( nf ** 1.8 )
    if neff < required_stats :
        rs = nice_print ( required_stats ) 
        return 'neff_pooled<%s' % rs       

    return ''

# =============================================================================
## @class DensityReweighter
#  Abstract base class for immutable, adaptive density-ratio reweighting.
#  Implements 1-, 2-, and 4-stream signed-measure density estimations.
#  Executes training immediately upon instantiation and seals resulting weights.
class DensityReweighter ( Reweighter, abc.ABC ) :
    """ Abstract base class for immutable, adaptive density-ratio reweighting.
    
    Implements 1-, 2-, and 4-stream signed-measure density estimations.
    Executes training immediately upon instantiation and seals resulting weights.
    """

    # ==========================================================================
    ## @brief Initialize and fit the density ratio reweighter ensemble.
    #  @param original Features array for original dataset.
    #  @param target Features array for target dataset.
    #  @param original_weight Event weights for original sample.
    #  @param target_weight Event weights for target sample.
    #  @param clip_threshold Maximum allowable weight ratio before clipping.
    #  @param n_splits Number of cross-validation folds.
    #  @param random_state Seed for reproducible fold splitting.
    #  @param store_original_weights Store computed ratios/weights for original sample.
    #  @param progress Display progress bar during CV training loops.
    #  @param params Keyword parameters forwarded to base classifier.
    def __init__( self                           , * ,
                  original                       ,
                  target                         ,
                  original_weight        = None  , 
                  target_weight          = None  ,
                  clip_threshold         = 1.e+4 ,
                  n_splits               = 5     ,
                  random_state           = 42    ,
                  store_original_weights = True  ,
                  progress               = True  , **params ) :
        """ Initialize and fit the density ratio reweighter ensemble.

        :param original: Features array for original dataset.
        :param target: Features array for target dataset.
        :param original_weight: Initial weights for original sample.
        :param target_weight: Initial weights for target sample.
        :param clip_threshold: Upper boundary for weight ratios.
        :param n_splits: Number of CV folds.
        :param random_state: Random state for KFold splits.
        :param store_original_weights: Cache initial reweighted weights.
        :param progress: Enable/disable progress bar.
        :param params: Model parameters.
        """
        if not isinstance ( n_splits, int ) : raise TypeError  ( "Invalid `n_splits' type %s"  % typename( n_splits ) )
        if not 0 <= n_splits <= 1000        : raise ValueError ( "Invalid `n_splits' value %s" % n_splits )
        if not isinstance( clip_threshold, num_types ) :
            raise TypeError ( "Invalid `clip_threshold' type %s" % typename( clip_threshold ) )
        if not 0 < clip_threshold :
            raise ValueError ( "Invalid `clip_threshold' value %s" % clip_threshold )
        
        check_all ( data1   = original          ,
                    data2   = target            ,
                    weight1 = original_weight   ,
                    weight2 = target_weight     ,
                    where   = typename ( self ) )   
        
        self.__progress       = True if progress else False 
        self.__clip_threshold = float ( clip_threshold )
        self.__n_splits       = n_splits
        self.__random_state   = random_state

        self.__fitted_models       = {}
        self.__priors              = {}
        self.__target_weights_info = {}
        self.__norm_factor         = numpy.float32( 1.0 )
        self.__mode                = None
        self.__scale_factors       = {}

        params [ 'n_jobs' ] = num_jobs ( params )

        reg_case = self.needs_regularization ( original        = original        ,
                                               target          = target          ,
                                               original_weight = original_weight ,
                                               target_weight   = target_weight   ) 
        
        if reg_case :
            n_features = num_features ( original )
            neff_orig  = nEff ( original  , original_weight )
            neff_targ  = nEff ( target    , target_weight   ) 
            neff       = min  ( neff_orig , neff_targ )
            
            params.update ( self.regularization ( params , n_features , neff ) )

            ESR = params.get ( 'early_stopping_rounds' , 15 )
            if ESR is None : params [ 'early_stopping_rounds' ] = None 
            else           : params [ 'early_stopping_rounds' ] = min ( 15 , ESR )

            title = '%s strong regularization' % typename ( self ) 
            table = map2table_ex ( params    , 
                                   header    = ( 'Parameter' , 'type' , 'value' ) ,
                                   alignment = 'rcw'  , 
                                   prefix    = '# '   ,
                                   title     = title  )
            
            logger.info ( "%s is applied, case %s:\n%s" % ( title , reg_case , table ) )

        self.__original_ratios             = None
        self.__original_reweighted_weights = None

        Reweighter.__init__ ( self            ,
                              original        = original        ,
                              target          = target          , 
                              original_weight = original_weight ,
                              target_weight   = target_weight   , **params )

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

        if original_ratios is not None and 0 < len ( original_ratios ) : 
            if not numpy.all ( numpy.isfinite ( original_ratios ) ) :
                logger.error ( "%s: NaN/Inf found in original_ratios" % typename ( self ) )
            if numpy.any ( original_ratios < 0 ):
                logger.error ( "%s: Negative values found for ratio r(x)" % typename ( self ) )

            r_min = float ( numpy.min ( original_ratios ) )
            r_max = float ( numpy.max ( original_ratios ) )
            r_std = float ( numpy.std ( original_ratios ) )

            if numpy.isclose ( r_min, r_max, atol = 1e-5 ) or r_std < 1e-6 :
                r1 = nice_print ( r_min )
                r2 = nice_print ( r_std ) 
                logger.warning( "%s: All ratios are constant r(x) = %s (std=%s)" % ( self.method , r1 , r2 ) ) 
                
            neff_before = nEff ( original , original_weight             )
            neff_after  = nEff ( original , original_reweighted_weights )
            
            if neff_after < 0.10 * neff_before :
                n1 = nice_print ( neff_before )
                n2 = nice_print ( neff_after  )                
                logger.warning ( "%s: Large degradation of nEff %s %s %s" % ( self.method , n1 , arrow_right , n2 ) ) 

            clipped_count = numpy.sum ( original_ratios >= ( self.__clip_threshold * 0.99 ) )
            if 0 < clipped_count : 
                frac = 100.0 * clipped_count / len ( original_ratios ) 
                if 1 < frac : logger.warning ( "%s: Too many clipped events %.1f[%%]" % ( self.method , frac ) )
                                 
        cnt1 = SE()
        cnt2 = SE()
        for r in original_ratios             : cnt1.add ( r )
        for w in original_reweighted_weights : cnt2.add ( w )
        
        if not self.silent :
            counters = { 'Ratios' : cnt1 , 'Weights' : cnt2 }
            title    = '%s weight info' % self.method 
            table    = table_counters ( counters , prefix = '# ' , title = title )
            logger.info ( '%s:\n%s' % ( title , table ) )  
                                                
        if store_original_weights :
            self.__original_ratios             = original_ratios
            self.__original_reweighted_weights = original_reweighted_weights

    # =================================================================================
    ## Get the number of splits/folds used for cross-validation.
    #  @return Integer count of CV folds.
    @property
    def n_splits ( self ) :
        """ Get the number of splits/folds used for cross-validation.
        """
        return self.__n_splits

    # =================================================================================
    ## Check if strong regularization is required considering both samples.
    #  @param original Features array for original sample.
    #  @param target Features array for target sample.
    #  @param original_weight Event weights for original sample.
    #  @param target_weight Event weights for target sample.
    #  @return Reason string if regularization is needed, else empty string.
    def needs_regularization ( self ,
                               original                       ,
                               target                         ,
                               original_weight        = None  , 
                               target_weight          = None  ) :
        """ Check if strong regularization is required considering both samples.
        """
        check_all ( data1   = original            ,
                    data2   = target              ,
                    weight1 = original_weight     ,
                    weight2 = target_weight       ,
                    where   = "%s:needs_regularization" % typename ( self ) )
        return RW_needs_regularization ( original        = original        ,
                                         target          = target          ,
                                         original_weight = original_weight ,
                                         target_weight   = target_weight   )

    # =================================================================================
    ## Get computed density ratio factors r(x) for original sample.
    #  @return Array of density ratios r(x).
    @property
    def original_ratios( self ) :
        """ Get computed density ratio factors r(x) for original sample.
        """
        return self.__original_ratios

    # =================================================================================
    ## Get final reweighted event weights (w * r) for original sample.
    #  @return Array of reweighted weights.
    @property
    def original_reweighted_weights( self ) :
        """ Get final reweighted event weights (w * r) for original sample.
        """
        return self.__original_reweighted_weights

    # =================================================================================
    ## Get active stream decomposition scheme string.
    #  @return Stream mode identifier string.
    @property
    def mode ( self ) :
        """ Get active stream decomposition scheme string.\
        """
        return self.__mode

    # =================================================================================
    ## Get progress bar display setting.
    #  @return True if progress bar is enabled.
    @property
    def progress ( self ) :
        """ Get progress bar display setting.
        """
        return self.__progress
    
    # =================================================================================
    ## Get full reweighter configuration dictionary.
    #  @return Configuration parameters dict.
    @property
    def config ( self ) :
        """ Get full reweighter configuration dictionary.
        """
        conf = {}
        conf.update( super().config if hasattr( super(), "config" ) else {} )
        conf [ "progress" ]                    = self.progress
        conf [ "mode" ]                        = self.mode
        conf [ "clip_threshold" ]              = self.__clip_threshold
        conf [ "n_splits" ]                    = self.__n_splits
        conf [ "original_ratios" ]             = self.__original_ratios             is not None
        conf [ "original_reweighted_weights" ] = self.__original_reweighted_weights is not None
        return conf

    # =================================================================================
    ## Abstract method to train a single base classifier on one fold.
    #  @param X_train Training features.
    #  @param y_train Training labels.
    #  @param w_train Training weights or None.
    #  @param X_val Validation features.
    #  @param y_val Validation labels.
    #  @param w_val Validation weights or None.
    #  @return Tuple of (trained_model, val_predictions).
    @abc.abstractmethod
    def _train_single_model ( self    ,
                              X_train ,
                              y_train ,
                              w_train ,
                              X_val   ,
                              y_val   ,
                              w_val   ) :
        """ Train a single base classifier model on a fold dataset.
        """
        raise NotImplementedError

    # =================================================================================
    ## Abstract method to predict target probability p(y=1|x) using a single model.
    #  @param model Trained model instance.
    #  @param X Features matrix to evaluate.
    #  @return Array of predicted target probabilities.
    @abc.abstractmethod
    def _predict_single_model( self  ,
                               model ,
                               X     ) :
        """ Predict target probabilities using a trained single model.
        """
        raise NotImplementedError

    # =================================================================================
    ## Execute fitting pipeline and calculate ratios for original sample.
    #  @param original Features array for original sample.
    #  @param original_weight Weights for original sample or None.
    #  @param target Features array for target sample.
    #  @param target_weight Weights for target sample or None.
    #  @return Tuple of (original_ratios, original_reweighted_weights).
    def __fit_and_compute ( self            ,
                            original        ,
                            original_weight ,
                            target          ,
                            target_weight   ) :
        """ Fit stream models and compute initial reweighting factors.
        """
        if original_weight is None and target_weight is None :
            self.__mode  = "1-stream"
            ratios_valid = self.__fit_eval_stream ( "base", original, None, target, None )

            clipped            = numpy.clip    ( ratios_valid,
                                                 numpy.float32 ( 0.0 ) ,
                                                 numpy.float32 ( self.__clip_threshold ) )
            self.__norm_factor = numpy.float32 ( len( original )
                                                 / ( numpy.sum ( clipped, dtype = numpy.float64 ) + 1e-12 ) )

            original_ratios = clipped * self.__norm_factor
            return original_ratios, original_ratios.copy()

        w_orig = original_weight if original_weight is not None else numpy.ones ( len ( original ) , dtype = numpy.float32 )
        w_targ = target_weight   if target_weight   is not None else numpy.ones ( len ( target   ) , dtype = numpy.float32 )

        nz_orig_mask = w_orig != 0
        nz_targ_mask = w_targ != 0

        X_orig_valid, w_orig_valid = original [ nz_orig_mask ] , w_orig [ nz_orig_mask ]
        X_targ_valid, w_targ_valid = target   [ nz_targ_mask ] , w_targ [ nz_targ_mask ]

        mask_orig_pos, mask_orig_neg = w_orig_valid > 0, w_orig_valid < 0
        mask_targ_pos, mask_targ_neg = w_targ_valid > 0, w_targ_valid < 0

        has_neg_orig = numpy.any ( mask_orig_neg )
        has_neg_targ = numpy.any ( mask_targ_neg )

        if not has_neg_orig and not has_neg_targ :
            self.__mode  = "1-stream"
            ratios_valid = self.__fit_eval_stream ( "base"       ,
                                                    X_orig_valid ,
                                                    w_orig_valid ,
                                                    X_targ_valid ,
                                                    w_targ_valid )

        elif has_neg_orig and not has_neg_targ :
            self.__mode                   = "2-stream_orig"
            ratios_valid                  = numpy.zeros ( len( X_orig_valid ) , dtype = numpy.float32 )
            ratios_valid[ mask_orig_pos ] = self.__fit_eval_stream ( "orig_pos"                     ,
                                                                     X_orig_valid [ mask_orig_pos ] ,
                                                                     w_orig_valid [ mask_orig_pos ] ,
                                                                     X_targ_valid                   ,
                                                                     w_targ_valid                   )
            ratios_valid[ mask_orig_neg ] = self.__fit_eval_stream ( "orig_neg"                     ,
                                                                     X_orig_valid [ mask_orig_neg ] ,
                                                                     w_orig_valid [ mask_orig_neg ] ,
                                                                     X_targ_valid                   ,
                                                                     w_targ_valid                   )

        elif not has_neg_orig and has_neg_targ :
            self.__mode            = "2-stream_targ"
            X_targ_pos, w_targ_pos = X_targ_valid [ mask_targ_pos ] , w_targ_valid [ mask_targ_pos ] 
            X_targ_neg, w_targ_neg = X_targ_valid [ mask_targ_neg ] , w_targ_valid [ mask_targ_neg ]

            w_pos_sum                  = numpy.float32 ( numpy.sum ( w_targ_pos, dtype = numpy.float64 ) )
            w_neg_sum                  = numpy.float32 ( numpy.abs ( numpy.sum ( w_targ_neg, dtype = numpy.float64 ) ) )
            
            self.__target_weights_info = { "W_pos": w_pos_sum, "W_neg": w_neg_sum }
            w_total                    = w_pos_sum + w_neg_sum

            r_pos_pos = self.__fit_eval_stream ( "pos_pos"    ,
                                                 X_orig_valid ,
                                                 w_orig_valid ,
                                                 X_targ_pos   ,
                                                 w_targ_pos   )
            r_pos_neg = self.__fit_eval_stream ( "pos_neg"    ,
                                                 X_orig_valid ,
                                                 w_orig_valid ,
                                                 X_targ_neg   ,
                                                 w_targ_neg   )

            ratios_valid = ( r_pos_pos * w_pos_sum - r_pos_neg * w_neg_sum ) / ( w_total + numpy.float32( 1e-12 ) )

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

            r_pos_pos = self.__fit_eval_stream ( "pos_pos", X_orig_pos, w_orig_pos, X_targ_pos, w_targ_pos )
            r_pos_neg = self.__fit_eval_stream ( "pos_neg", X_orig_pos, w_orig_pos, X_targ_neg, w_targ_neg )
            r_neg_pos = self.__fit_eval_stream ( "neg_pos", X_orig_neg, w_orig_neg, X_targ_pos, w_targ_pos )
            r_neg_neg = self.__fit_eval_stream ( "neg_neg", X_orig_neg, w_orig_neg, X_targ_neg, w_targ_neg )

            ratios_valid = numpy.zeros ( len ( X_orig_valid ) , dtype = numpy.float32 )
            denom                         = w_total + numpy.float32( 1e-12 )
            ratios_valid[ mask_orig_pos ] = ( r_pos_pos * w_pos_sum - r_pos_neg * w_neg_sum ) / denom
            ratios_valid[ mask_orig_neg ] = ( r_neg_pos * w_pos_sum - r_neg_neg * w_neg_sum ) / denom

        ratios_valid_norm = self.__normalize_and_clip ( ratios_valid, w_orig_valid )

        original_ratios                 = numpy.zeros( len( w_orig ), dtype = numpy.float32 )
        original_ratios[ nz_orig_mask ] = ratios_valid_norm
        original_reweighted_weights     = w_orig * original_ratios

        return original_ratios, original_reweighted_weights

    # =================================================================================
    ## Fit fold models for a specific stream decomposition component.
    #  @param stream_key Identifier key for the stream component.
    #  @param X_orig_sub Sub-dataset features for original sample.
    #  @param w_orig_sub Sub-dataset weights for original sample.
    #  @param X_targ_sub Sub-dataset features for target sample.
    #  @param w_targ_sub Sub-dataset weights for target sample.
    #  @return Array of out-of-fold calculated density ratios.
    def __fit_eval_stream ( self       ,
                            stream_key ,
                            X_orig_sub ,
                            w_orig_sub ,
                            X_targ_sub ,
                            w_targ_sub ) :
        """ Fit stream models using KFold cross validation.
        """
        X_comb = numpy.vstack( [ X_orig_sub, X_targ_sub ] )
        y_comb = numpy.hstack( [ numpy.zeros ( len( X_orig_sub ) , dtype = numpy.float32 ) ,
                                 numpy.ones  ( len( X_targ_sub ) , dtype = numpy.float32 ) ] )

        if w_orig_sub is None and w_targ_sub is None :
            scale_factor = numpy.float32 ( len ( X_targ_sub ) / len ( X_orig_sub ) )
            w_orig_scaled_for_tr = numpy.full(len(X_orig_sub), scale_factor, dtype=numpy.float32)
            w_comb = numpy.hstack([w_orig_scaled_for_tr, numpy.ones(len(X_targ_sub), dtype=numpy.float32)])
        else :
            w_orig_abs = numpy.abs ( w_orig_sub ) if w_orig_sub is not None else numpy.ones ( len ( X_orig_sub ), dtype = numpy.float32 )
            w_targ_abs = numpy.abs ( w_targ_sub ) if w_targ_sub is not None else numpy.ones ( len ( X_targ_sub ), dtype = numpy.float32 )
            
            sum_w_orig = numpy.sum ( w_orig_abs )
            sum_w_targ = numpy.sum ( w_targ_abs )
            
            scale_factor = ( sum_w_targ / sum_w_orig ) if sum_w_orig > 0 else numpy.float32 ( 1.0 )
            
            w_orig_scaled = w_orig_abs * scale_factor
            w_comb        = numpy.hstack ( [ w_orig_scaled , w_targ_abs ] )

        oof_raw = numpy.zeros( len( X_comb ), dtype = numpy.float32 )

        if 1 < self.n_splits :
            from sklearn.model_selection import StratifiedKFold
            skf     = StratifiedKFold ( n_splits     = self.n_splits     , 
                                        shuffle      = True              , 
                                        random_state = self.random_state )
            splits  = skf.split(X_comb, y_comb)
            n_folds = self.n_splits
        else:
            indices = numpy.arange ( len ( X_comb ) )
            splits  = [ ( indices , indices ) ]
            n_folds = 1

        stream_models = []
        for train_idx, val_idx in progress_bar ( splits      ,
                                                 max_value   = n_folds        ,
                                                 description = 'Folds:'       , 
                                                 silent      = not self.progress or self.silent or 1 >= self.n_splits ) :
            
            X_tr, y_tr = X_comb [ train_idx ] , y_comb [ train_idx ]
            X_va, y_va = X_comb [ val_idx   ] , y_comb [ val_idx ]

            w_tr = w_comb [ train_idx ] if w_comb is not None else None
            w_va = w_comb [ val_idx   ] if w_comb is not None else None

            model, _ = self._train_single_model ( X_tr, y_tr, w_tr, X_va, y_va, w_va )
            stream_models.append( model )
            
            raw_val_p = self._predict_single_model ( model, X_va )
            if raw_val_p.ndim > 1 :
                raw_val_p = raw_val_p[:, 1]
                
            oof_raw[ val_idx ] = raw_val_p.astype ( numpy.float32, copy = False )

        self.__fitted_models[ stream_key ] = stream_models 
        self.__scale_factors[ stream_key ] = scale_factor
        self.__priors[ stream_key ]        = numpy.float32( 1.0 )

        p_orig = oof_raw[ : len ( X_orig_sub ) ]
        p_orig_clipped = numpy.clip ( p_orig, numpy.float32 ( 1e-4 ) , numpy.float32( 1.0 - 1e-4 ) )
        
        ratios = p_orig_clipped / ( numpy.float32( 1.0 ) - p_orig_clipped )
        return ratios

    # =================================================================================
    ## Calculate average stream density ratios r(x) for unseen sample features.
    #  @param stream_key Stream identifier key.
    #  @param X Feature matrix to evaluate.
    #  @return Array of calculated stream density ratios.
    def __predict_stream_ratios( self, stream_key, X ):
        """ Calculate average stream density ratios for unseen features.
        """
        models = self.__fitted_models[ stream_key ]

        ratios_list = []
        for model in models :
            p = self._predict_single_model ( model, X )
            if p.ndim > 1 :
                p = p[:, 1]
                
            eps = 1e-4
            p_clipped = numpy.clip ( p, numpy.float32 ( eps ), numpy.float32 ( 1.0 - eps ) )
            
            stream_ratios = p_clipped / ( numpy.float32( 1.0 ) - p_clipped )
            ratios_list.append( stream_ratios )
            
        return numpy.mean ( ratios_list, axis = 0, dtype = numpy.float32 )

    # =================================================================================
    ## Default implementation to predict class probabilities from a single model.
    #  @param model Trained model instance.
    #  @param X Input feature array.
    #  @return Array of positive class probabilities p(y=1|x).
    def _predict_single_model ( self, model, X ):
        """ Predict class probabilities from a single model instance.
        """
        best_iter = getattr ( model, 'best_iteration_', None ) or getattr( model, 'best_iteration', None )
        kwargs = {}
        if best_iter is not None and best_iter > 0 :
            kwargs[ 'ntree_end' ] = best_iter            
        p = model.predict_proba( X, **kwargs )[:, 1]
        return p.astype( numpy.float32, copy = False )

    # =================================================================================
    ## Clip extreme ratios to upper limit and normalize weight integral.
    #  @param ratios Computed raw density ratios.
    #  @param w_orig Original sample event weights.
    #  @return Normalized and clipped density ratios array.
    def __normalize_and_clip ( self   ,
                               ratios ,
                               w_orig ) :
        """ Clip extreme ratios and normalize sum of weights.
        """
        clipped_ratios = numpy.clip ( ratios, numpy.float32 ( 0.0 ) , numpy.float32 ( self.__clip_threshold ) )
        orig_sum       = numpy.sum  ( w_orig, dtype = numpy.float64 )
        reweighted_sum = numpy.sum  ( w_orig * clipped_ratios, dtype = numpy.float64 )
        self.__norm_factor = numpy.float32 ( orig_sum / ( reweighted_sum + 1e-12 ) )
        return clipped_ratios * self.__norm_factor

    # =================================================================================
    ## Apply fitted model ensemble to compute reweighted event weights for new data.
    #  @param original Features matrix for new sample.
    #  @param original_weight Initial event weights for new sample.
    #  @return Array of calculated reweighted event weights.
    def weights ( self                   ,
                  original               ,
                  original_weight = None ) :
        """ Apply fitted model ensemble to compute reweighted event weights for new data.

        :param original: Features array for target/original sample.
        :param original_weight: Optional initial weights array.
        :return: Array of reweighted event weights.
        """
        X_new_f32 = original.astype ( numpy.float32, copy = False )

        if not valid_data_shape   ( X_new_f32                   ) : raise TypeError ( "Invalid `original` type/shape: %s" % typename ( original ) )
        if not valid_weight       ( original_weight             ) : raise TypeError ( "Invalid `original_weight`!" )        
        if not compatible_weights ( X_new_f32 , original_weight ) : raise TypeError ( "Incompatible `original` data/weight!" )
        if self.n_features != num_features ( X_new_f32 )          : raise TypeError ( "Invalid #features!!")

        if original_weight is None and self.__mode == "1-stream" :
            ratios  = self.__predict_stream_ratios( "base", X_new_f32 )
            clipped = numpy.clip ( ratios,
                                   numpy.float32 ( 0.0 ) ,
                                   numpy.float32 ( self.__clip_threshold ) )
            return clipped * self.__norm_factor

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
                ratios_valid[ mask_pos ] = self.__predict_stream_ratios( "orig_pos", X_valid[ mask_pos ] )
            if numpy.any( mask_neg ) :
                ratios_valid[ mask_neg ] = self.__predict_stream_ratios( "orig_neg", X_valid[ mask_neg ] )

        elif self.__mode == "2-stream_targ" :
            w_pos   = self.__target_weights_info[ "W_pos" ]
            w_neg   = self.__target_weights_info[ "W_neg" ]
            w_total = w_pos + w_neg + numpy.float32( 1e-12 )

            r_pos_pos    = self.__predict_stream_ratios( "pos_pos", X_valid )
            r_pos_neg    = self.__predict_stream_ratios( "pos_neg", X_valid )
            ratios_valid = ( r_pos_pos * w_pos - r_pos_neg * w_neg ) / w_total

        elif self.__mode == "4-stream" :
            w_pos   = self.__target_weights_info[ "W_pos" ]
            w_neg   = self.__target_weights_info[ "W_neg" ]
            w_total = w_pos + w_neg + numpy.float32( 1e-12 )

            if numpy.any( mask_pos ) :
                X_p                      = X_valid[ mask_pos ]
                r_pos_pos                = self.__predict_stream_ratios( "pos_pos", X_p )
                r_pos_neg                = self.__predict_stream_ratios( "pos_neg", X_p )
                ratios_valid[ mask_pos ] = ( r_pos_pos * w_pos - r_pos_neg * w_neg ) / w_total

            if numpy.any( mask_neg ) :
                X_n                      = X_valid[ mask_neg ]
                r_neg_pos                = self.__predict_stream_ratios( "neg_pos", X_n )
                r_neg_neg                = self.__predict_stream_ratios( "neg_neg", X_n )
                ratios_valid[ mask_neg ] = ( r_neg_pos * w_pos - r_neg_neg * w_neg ) / w_total

        clipped                = numpy.clip ( ratios_valid,
                                              numpy.float32 ( 0.0 ),
                                              numpy.float32 ( self.__clip_threshold ) )
        final_ratios           = numpy.zeros( len( w_new_f32 ), dtype = numpy.float32 )
        final_ratios[ nz_mask ] = clipped * self.__norm_factor

        return final_ratios if original_weight is None else final_ratios * w_new_f32

# =============================================================================
## @class  LightGBMDensityReweighter
#  Density ratio reweighter using LightGBM as the underlying classifier.
class LightGBMDensityReweighter ( DensityReweighter ) :
    """ Density ratio reweighter using LightGBM as the underlying classifier.
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
        """ Initialize LightGBM density reweighter.
        """
        config = {
            'objective'             : 'binary'            ,
            'metric'                : 'binary_logloss'    ,
            'n_estimators'          : DEFAULT_ESTIMATORS  ,
            'learning_rate'         : 0.03                ,
            'max_depth'             : MAX_DEPTH           ,
            'max_bin'               : 2048                ,
            'num_leaves'            : 31                  ,
            'min_child_samples'     : 30                  ,
            'min_child_weight'      : 1e-3                ,
            'reg_alpha'             : 0.1                 ,
            'reg_lambda'            : 2.0                 ,
            'subsample'             : 0.8                 ,
            'subsample_freq'        : 1                   ,
            'colsample_bytree'      : 0.8                 ,
            'path_smooth'           : 1.0                 ,
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

    # =========================================================================
    ## Return the method identifier name.
    #  @return Method string identifier.
    @property
    def method ( self ) :
        """ Return the method identifier name.
        """
        return method_LGBM 

    # =========================================================================
    ## Apply soft regularization rules for low statistics or low dimensions.
    #  @param params Current parameters dictionary.
    #  @param n_features Number of phase space features.
    #  @param n_samples Effective sample size.
    #  @return Updated parameters dictionary.
    def regularization ( self       ,
                         params     , 
                         n_features ,
                         n_samples  ) :
        """ Apply soft regularization rules for low statistics or low dimensions.
        """
        params [ 'n_estimators'      ] = min ( REGULARIZED_ESTIMATORS , params.get ( 'n_estimators' , REGULARIZED_ESTIMATORS ) )
        params [ 'learning_rate'     ] = 0.05
        params [ 'max_depth'         ] = REG_DEPTH       
        params [ 'num_leaves'        ] = 2**REG_DEPTH - 1       
        params [ 'min_child_samples' ] = max ( 50 , int ( n_samples  * 0.01 ) )
        params [ 'min_child_weight'  ] = 1.e-7 
        params [ 'reg_alpha'         ] = 1.0    
        params [ 'reg_lambda'        ] = 5.0 
        params [ 'subsample'         ] = 1.0    
        params [ 'colsample_bytree'  ] = 1.0  
        params [ 'early_stopping_rounds' ] = None
        params [ 'min_data_in_bin'   ] = 1      
        
        if 'path_smooth' in params : 
            params.pop ( 'path_smooth' )

        return params
    
    # =========================================================================
    ## Train single LightGBM model on fold data.
    #  @param X_train Training features array.
    #  @param y_train Training binary labels.
    #  @param w_train Training event weights or None.
    #  @param X_val Validation features array.
    #  @param y_val Validation binary labels.
    #  @param w_val Validation event weights or None.
    #  @return Tuple of (trained_lightgbm_booster, val_predictions).
    def _train_single_model ( self    ,
                              X_train , y_train , w_train ,
                              X_val   , y_val   , w_val   ) : 
        """ Train single LightGBM model on fold data.
        """
        import lightgbm as LightGBM 

        w_tr = numpy.ascontiguousarray ( w_train , dtype = numpy.float32 ) if w_train is not None else None
        w_va = numpy.ascontiguousarray ( w_val   , dtype = numpy.float32 ) if w_val   is not None else None

        trn_data = LightGBM.Dataset ( X_train , label = y_train , weight = w_tr , free_raw_data = False )
        val_data = LightGBM.Dataset ( X_val   , label = y_val   , weight = w_va , reference = trn_data , free_raw_data = False )

        params = self.params.copy ()
        num_boost_round       = params.pop ( 'num_boost_round' , None ) or params.pop ( 'n_estimators' , 400 )
        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , None )

        callbacks = []
        if early_stopping_rounds is not None : 
            callbacks.append ( LightGBM.early_stopping ( stopping_rounds = early_stopping_rounds , verbose = False ) )
            
        model = LightGBM.train ( params          = params          ,
                                 train_set       = trn_data        ,
                                 num_boost_round = num_boost_round ,
                                 valid_sets      = [ val_data ]    ,
                                 callbacks       = callbacks       )
        
        val_preds = self._predict_single_model ( model , X_val )
        return model , val_preds

    # =========================================================================
    ## Predict probabilities using a LightGBM model.
    #  @param model Trained LightGBM booster instance.
    #  @param X Features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model( self, model, X ):
        """ Predict probabilities using a LightGBM model.
        """
        best_iter = getattr( model, 'best_iteration', 0 )
        kwargs = {}
        if best_iter is not None and best_iter > 0 :
            kwargs[ 'num_iteration' ] = best_iter            
        p = model.predict( X, **kwargs )
        return p.astype( numpy.float32, copy = False )

# =============================================================================
## @class XGBoostDensityReweighter 
#  Density ratio reweighter using XGBoost as the underlying classifier.
class XGBoostDensityReweighter ( DensityReweighter ): 
    """ Density ratio reweighter using XGBoost as the underlying classifier.
    """

    ## @brief Initialize XGBoost density reweighter.
    #  @param original Features array for original dataset.
    #  @param target Features array for target dataset.
    #  @param original_weight Initial weights for original sample (optional).
    #  @param target_weight Initial weights for target sample (optional).
    #  @param store_original_weights If True, store computed results for original sample.
    #  @param params Additional XGBoost parameters.
    def __init__(  self , * , 
                   original               ,
                   target                 ,
                   original_weight        = None ,
                   target_weight          = None ,
                   store_original_weights = True , **params ) :
        """Initialize XGBoost density reweighter."""
        config = {
            'objective'             : 'binary:logistic'   ,
            'eval_metric'           : 'logloss'           ,
            'n_estimators'          : DEFAULT_ESTIMATORS  ,
            'learning_rate'         : 0.03                ,
            'max_depth'             : MAX_DEPTH           ,
            'min_child_weight'      : 0.1                 ,
            'gamma'                 : 0.001               ,
            'reg_alpha'             : 0.1                 ,
            'reg_lambda'            : 2.0                 ,
            'subsample'             : 0.8                 ,
            'colsample_bytree'      : 0.8                 ,
            'tree_method'           : 'hist'              ,
            'early_stopping_rounds' : None                ,
            'verbosity'             : 0                   ,
            'n_jobs'                : -1                  ,
        }
        
        config.update ( params )
        if 'num_boost_round' in config : config [ 'n_estimators' ] = config.pop ( 'num_boost_round' )

        super().__init__ ( original               = original               ,
                           target                 = target                 ,
                           original_weight        = original_weight        ,
                           target_weight          = target_weight          ,
                           store_original_weights = store_original_weights , **config )
        
    # =========================================================================
    ## Return the method identifier name.
    #  @return Method string identifier.
    @property
    def method ( self ) :
        """ Return the method identifier name.
        """
        return method_XGB 
    
    # =========================================================================
    ## Dynamic regularization rules for XGBoost.
    #  @param params Current parameters dictionary.
    #  @param n_features Number of features.
    #  @param n_samples Effective sample statistics.
    #  @return Updated parameters dictionary.
    def regularization ( self , params , n_features , n_samples ) :
        """ Dynamic regularization rules for XGBoost.
        """
        params [ 'n_estimators'          ] = min ( REGULARIZED_ESTIMATORS , params.get ( 'n_estimators' , REGULARIZED_ESTIMATORS ) )
        params [ 'learning_rate'         ] = 0.05
        params [ 'max_depth'             ] = REG_DEPTH
        params [ 'min_child_weight'      ] = 1.e-7 
        params [ 'gamma'                 ] = 0.0
        params [ 'reg_alpha'             ] = 0.0
        params [ 'reg_lambda'            ] = 0.0
        params [ 'subsample'             ] = 1.0
        params [ 'colsample_bytree'      ] = 1.0
        params [ 'early_stopping_rounds' ] = None
        params [ 'tree_method'           ] = 'exact'
        
        return params

    # =========================================================================
    ## Train single XGBoost booster on fold data.
    #  @param X_train Training features array.
    #  @param y_train Training binary labels.
    #  @param w_train Training event weights or None.
    #  @param X_val Validation features array.
    #  @param y_val Validation binary labels.
    #  @param w_val Validation event weights or None.
    #  @return Tuple of (fitted_xgb_booster, val_predictions).
    def _train_single_model ( self,
                              X_train , y_train , w_train ,
                              X_val   , y_val   , w_val   ) :
        """ Train single XGBoost booster on fold data.
        """
        import xgboost as XGBoost

        dtrain = XGBoost.DMatrix ( X_train , label = y_train , weight = w_train )
        dval   = XGBoost.DMatrix ( X_val   , label = y_val   , weight = w_val   )

        params = {}
        params.update ( self.params )

        num_boost_round = params.pop ( 'num_boost_round' , None ) or params.pop ( 'n_estimators' , None ) or 500
        if not isinstance ( num_boost_round , int ) or num_boost_round <= 10 or num_boost_round >= 10000 :
            num_boost_round = 500

        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , None )
        if not isinstance ( early_stopping_rounds , int ) or early_stopping_rounds <= 1 or early_stopping_rounds >= num_boost_round :
            early_stopping_rounds = None

        train_kwargs = {}
        if early_stopping_rounds is not None :
            train_kwargs [ 'evals' ]                 = [ ( dval , "val" ) ]
            train_kwargs [ 'early_stopping_rounds' ] = early_stopping_rounds

        model = XGBoost.train ( params          = params,
                                dtrain          = dtrain,
                                num_boost_round = num_boost_round,
                                verbose_eval    = False,
                                **train_kwargs )

        predict_kwargs = {}
        best_iter = getattr ( model , "best_iteration" , None )
        if early_stopping_rounds is not None and best_iter is not None and 0 < best_iter :
            predict_kwargs [ 'iteration_range' ] = ( 0 , best_iter + 1 )
            
        val_preds = model.predict ( dval , **predict_kwargs )
        return model , val_preds.astype ( numpy.float32 , copy = False )
    
    # =========================================================================
    ## Predict probabilities using an XGBoost model.
    #  @param model Trained XGBoost Booster.
    #  @param X Input features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model( self, model, X ):
        """ Predict probabilities using an XGBoost model.
        """
        import xgboost as XGBoost
        dmat = XGBoost.DMatrix ( X )
        best_iter = getattr ( model, 'best_iteration', None )
        kwargs = {}
        if best_iter is not None and 0 < best_iter :
            kwargs [ 'iteration_range' ] = ( 0, best_iter + 1 )
        p = model.predict ( dmat , **kwargs )
        return p.astype ( numpy.float32, copy = False )

# =============================================================================
## @class CatBoostDensityReweighter
#  Density ratio reweighter using CatBoost as the underlying classifier.
class CatBoostDensityReweighter ( DensityReweighter ) : 
    """ Density ratio reweighter using CatBoost as the underlying classifier.
    """
    # =========================================================================
    ## Initialize CatBoost density reweighter.
    #  @param original Features array for original dataset.
    #  @param target Features array for target dataset.
    #  @param original_weight Initial weights for original sample (optional).
    #  @param target_weight Initial weights for target sample (optional).
    #  @param store_original_weights If True, store computed results for original sample.
    #  @param params Additional CatBoost parameters.
    def __init__( self , * , 
                   original               ,
                   target                 ,
                   original_weight        = None ,
                   target_weight          = None ,
                   store_original_weights = True  , **params ) :        
        """ Initialize CatBoost density reweighter.
        """
        config = {
            'loss_function'         : 'Logloss'           ,
            'eval_metric'           : 'Logloss'           ,
            'n_estimators'          : DEFAULT_ESTIMATORS  ,
            'learning_rate'         : 0.03                ,
            'depth'                 : MAX_DEPTH           ,
            'l2_leaf_reg'           : 2.0                 ,
            'min_child_samples'     : 30                  ,
            'subsample'             : 0.8                 ,
            'random_strength'       : 1.0                 ,
            'early_stopping_rounds' : None                ,
            'verbose'               : False               ,
            'thread_count'          : 1                   ,
            'boosting_type'         : 'Plain'             ,
        }
        config.update ( params )
        
        if 'n_jobs'     in config : config [ 'thread_count' ] = config.pop ( 'n_jobs'     )  
        if 'iterations' in config : config [ 'n_estimators' ] = config.pop ( 'iterations' ) 

        config [ 'boosting_type' ] = config.get ( 'boosting_type' , 'Plain' )
        
        super().__init__ ( original               = original               ,
                           target                 = target                 ,
                           original_weight        = original_weight        ,
                           target_weight          = target_weight          ,
                           store_original_weights = store_original_weights , **config )
        
    # =========================================================================
    ## Return the method identifier name.
    #  @return Method string identifier.
    @property
    def method ( self ) :
        """ Return the method identifier name.
        """
        return method_CATB 

    # =========================================================================
    ## Dynamic regularization rules for CatBoost.
    #  @param params Current parameters dictionary.
    #  @param n_features Number of features.
    #  @param n_samples Effective sample size.
    #  @return Updated parameters dictionary.
    def regularization ( self , params , n_features , n_samples ) :
        """ Dynamic regularization rules for CatBoost.
        """
        params [ 'n_estimators'          ] = min ( REGULARIZED_ESTIMATORS , params.get ( 'n_estimators' , REGULARIZED_ESTIMATORS ) )
        params [ 'learning_rate'         ] = 0.1
        params [ 'depth'                 ] = REG_DEPTH
        
        if 'min_data_in_leaf' in params : params.pop ( 'min_data_in_leaf' , None )
        
        params [ 'min_child_samples'     ] = max ( 2, min ( 30, int ( n_samples * 0.0001 ) ) )
        params [ 'l2_leaf_reg'           ] = 1.0
        params [ 'subsample'             ] = 1.0
        params [ 'early_stopping_rounds' ] = None
        params [ 'thread_count'          ] = 1
        params [ 'boosting_type'         ] = 'Plain'
        
        return params

    # =========================================================================
    ## Train single CatBoost model on contiguous memory fold data.
    #  @param X_train Training features array.
    #  @param y_train Training binary labels.
    #  @param w_train Training weights array or None.
    #  @param X_val Validation features array.
    #  @param y_val Validation binary labels.
    #  @param w_val Validation weights array or None.
    #  @return Tuple of (fitted_catboost_model, val_predictions).
    def _train_single_model ( self    ,
                              X_train , y_train , w_train ,
                              X_val   , y_val   , w_val   ) :
        """ Train single CatBoost model on contiguous memory fold data.
        """
        import catboost as CatBoost

        X_tr = numpy.ascontiguousarray ( X_train , dtype = numpy.float32 )
        y_tr = numpy.ascontiguousarray ( y_train )
        w_tr = numpy.ascontiguousarray ( w_train , dtype = numpy.float32 ) if w_train is not None else None

        X_v  = numpy.ascontiguousarray ( X_val   , dtype = numpy.float32 )
        y_v  = numpy.ascontiguousarray ( y_val   )
        w_v  = numpy.ascontiguousarray ( w_val   , dtype = numpy.float32 ) if w_val is not None else None

        trn_pool = CatBoost.Pool ( X_tr , label = y_tr , weight = w_tr )
        val_pool = CatBoost.Pool ( X_v  , label = y_v  , weight = w_v  )

        params = {}
        params.update ( self.params )

        params [ 'thread_count'  ] = params.pop ( 'n_jobs'        , 1       ) 
        params [ 'boosting_type' ] = params.get ( 'boosting_type' , 'Plain' )

        iterations = params.pop ( 'iterations' , None ) or params.pop ( 'n_estimators' , None ) or 500
        if not isinstance ( iterations , int ) or iterations <= 10 or iterations >= 10000 :
            iterations = 500

        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , None )
        if not isinstance ( early_stopping_rounds , int ) or early_stopping_rounds <= 1 or early_stopping_rounds >= iterations :
            early_stopping_rounds = None

        params [ 'iterations' ] = iterations

        fit_kwargs = {}
        if early_stopping_rounds is not None :
            params     [ 'early_stopping_rounds' ] = early_stopping_rounds
            params     [ 'use_best_model'        ] = True
            fit_kwargs [ 'eval_set'              ] = val_pool

        model = CatBoost.CatBoostClassifier ( **params )
        model.fit ( trn_pool , verbose = False , **fit_kwargs )

        val_preds = model.predict_proba ( val_pool ) [ : , 1 ]
        return model , val_preds.astype ( numpy.float32 , copy = False )

    # =========================================================================
    ## Predict probabilities using a CatBoost model.
    #  @param model Fitted CatBoostClassifier instance.
    #  @param X Input features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model ( self , model , X ) :
        """ Predict probabilities using a CatBoost model.
        """
        X_clean   = numpy.ascontiguousarray ( X , dtype = numpy.float32 )
        best_iter = getattr ( model , 'best_iteration_' , None ) or getattr ( model , 'best_iteration' , None )
        kwargs    = {}
        if best_iter is not None and best_iter > 0 :
            kwargs [ 'ntree_end' ] = best_iter            
        p = model.predict_proba ( X_clean , **kwargs ) [ : , 1 ]
        return p.astype ( numpy.float32 , copy = False )

# ==============================================================================

# =============================================================================
## @class PyTorchDensityReweighter
#  Density ratio reweighter using PyTorch MLP as the underlying classifier
class PyTorchDensityReweighter(DensityReweighter):
    """ Density ratio reweighter using PyTorch MLP as the underlying classifier.
    """

    # =========================================================================
    ## Initialize PyTorch density reweighter.
    #  @param original Features array for original sample.
    #  @param target Features array for target sample.
    #  @param original_weight Initial weights for original sample (optional).
    #  @param target_weight Initial weights for target sample (optional).
    #  @param kwargs Additional training and architecture parameters.
    def __init__( self            , * ,
                  original        ,
                  target          ,
                  original_weight = None ,
                  target_weight   = None ,
                  progress        = True , **kwargs):
        """ Initialize PyTorch density reweighter.
        """
        import torch

        config = { 'hidden_dims'   : (128, 64, 32),
                   'dropout'       : 0.1,
                   'learning_rate' : 1e-3,
                   'weight_decay'  : 1e-4,
                   'batch_size'    : 1024,
                   'epochs'        : 150,
                   'patience'      : 15,
                   'device'        : 'cuda' if torch.cuda.is_available() else 'cpu',
                  }
        config.update(kwargs)

        self.__progress_epochs = True if progress else False 
        
        super().__init__ ( original        = original        ,
                           target          = target          ,
                           original_weight = original_weight ,
                           target_weight   = target_weight   ,
                           progress        = False           , **config )

    # =========================================================================
    ## Return method identifier name.
    #  @return Method string identifier.
    @property
    def method(self):
        """ Return method identifier name.
        """
        return method_TORCH 

    # =========================================================================
    ## Apply regularization rules for low statistics or low dimensions.
    #  @param params Current parameters dictionary.
    #  @param n_features Number of phase space features.
    #  @param n_samples Effective sample size.
    #  @return Updated parameters dictionary.
    def regularization ( self       ,
                         params     ,
                         n_features ,
                         n_samples  ) :
        """ Apply regularization rules for low statistics or low dimensions.
        """
        params [ 'hidden_dims'  ] = ( 64 , 32 )
        params [ 'dropout'      ] = 0.2
        params [ 'weight_decay' ] = 1e-2
        params [ 'patience'     ] = 10
        return params

    # =========================================================================
    ## Factory method to create a DensityMLP instance with lazy torch.nn import.
    #  @param in_features Number of input features.
    #  @return Instantiated PyTorch DensityMLP module.
    def _create_model ( self , in_features ):
        """ Factory method to create a DensityMLP instance with lazy torch.nn import.
        """
        import torch.nn as nn

        # =====================================================================
        ## @class DentiyMLP
        #  Multilayer Perceptron for binary classification and density ratio estimation.
        class DensityMLP(nn.Module):
            """ Multilayer Perceptron for binary classification and density ratio estimation.
            """
            # =================================================================
            def __init__(self, in_features, hidden_dims=(128, 64, 32), dropout=0.1):
                super().__init__()
                layers = []
                curr_dim = in_features
                for h_dim in hidden_dims:
                    layers.extend([
                        nn.Linear(curr_dim, h_dim),
                        nn.BatchNorm1d(h_dim),
                        nn.SiLU(),
                        nn.Dropout(dropout)
                    ])
                    curr_dim = h_dim
                layers.append(nn.Linear(curr_dim, 1))
                self.net = nn.Sequential(*layers)
                
            # =================================================================
            def forward(self, x):
                return self.net(x).squeeze(-1)  # Return raw logits

        return DensityMLP ( in_features = in_features                 ,
                            hidden_dims = self.params ['hidden_dims'] ,
                            dropout     = self.params [ 'dropout'   ] )

    # =========================================================================
    ## Train single PyTorch model on fold data with Early Stopping.
    #  @param X_train Training features array.
    #  @param y_train Training binary labels.
    #  @param w_train Training event weights or None.
    #  @param X_val Validation features array.
    #  @param y_val Validation binary labels.
    #  @param w_val Validation event weights or None.
    #  @return Tuple of (trained_pytorch_model, val_predictions).
    def _train_single_model(self, X_train, y_train, w_train, X_val, y_val, w_val):
        """ Train single PyTorch model on fold data with Early Stopping.
        """
        import torch
        import torch.nn as nn
        from   torch.utils.data      import DataLoader, TensorDataset
        from   sklearn.preprocessing import StandardScaler

        device     = torch.device(self.params.get('device', 'cpu'))
        batch_size = self.params['batch_size']

        # 1. Mandatory feature scaling (fit on train set only to prevent data leakage)
        scaler = StandardScaler()
        X_tr_scaled = scaler.fit_transform(X_train)
        X_va_scaled = scaler.transform(X_val)

        # 2. Prepare weights and ensure 1D shape using numpy.ravel()
        w_tr = w_train if w_train is not None else numpy.ones(len(y_train), dtype=numpy.float32)
        w_va = w_val   if w_val   is not None else numpy.ones(len(y_val),   dtype=numpy.float32)

        y_tr_flat = numpy.ravel ( y_train )
        w_tr_flat = numpy.ravel ( w_tr    )
        y_va_flat = numpy.ravel ( y_val   )
        w_va_flat = numpy.ravel ( w_va    )

        # 3. Create datasets and loaders
        ds_tr = TensorDataset  ( torch.tensor ( X_tr_scaled , dtype = torch.float32 ) ,
                                 torch.tensor ( y_tr_flat   , dtype = torch.float32 ) ,
                                 torch.tensor ( w_tr_flat   , dtype = torch.float32 )  )
        loader_tr = DataLoader ( ds_tr , batch_size = batch_size , shuffle = True, drop_last = True )
        
        ds_va = TensorDataset  ( torch.tensor ( X_va_scaled , dtype = torch.float32 ) ,
                                 torch.tensor ( y_va_flat   , dtype = torch.float32 ) ,
                                 torch.tensor ( w_va_flat   , dtype = torch.float32 ) )
        loader_va = DataLoader ( ds_va , batch_size = batch_size , shuffle = False  )

        # 4. Instantiate network via lazy factory method
        model = self._create_model ( in_features = X_train.shape [ 1 ] ) .to( device )

        optimizer = torch.optim.AdamW ( model.parameters(),
                                        lr          = self.params [ 'learning_rate' ] ,
                                        weight_decay= self.params [ 'weight_decay'  ] )
        criterion = nn.BCEWithLogitsLoss ( reduction = 'none' )

        # 5. Training loop with early stopping
        best_loss        = float('inf')
        best_state       = None
        patience_counter = 0

        patience         = self.params [ 'patience' ]
        
        nepochs = self.params [ 'epochs' ]
        no_pbar = self.silent or self.progress or not self.__progress_epochs

        with ProgressBar ( max_value = nepochs , silent = no_pbar , description = 'Epochs:' ) as pbar : 
        
            for epoch in range ( nepochs )  :
                
                model.train()
                
                for bx, by, bw in loader_tr:
                    bx, by, bw = bx.to(device), by.to(device), bw.to(device)
                    optimizer.zero_grad()
                    logits = model(bx)
                    loss   = (criterion(logits, by) * bw).mean()
                    loss.backward()
                    optimizer.step()

                # Validation step
                model.eval()
                val_loss_sum   = 0.0
                val_weight_sum = 0.0
                
                with torch.no_grad():
                    for bx, by, bw in loader_va:
                        bx, by, bw      = bx.to(device), by.to(device), bw.to(device)
                        val_logits      = model(bx)
                        batch_loss      = (criterion(val_logits, by) * bw).sum()
                        val_loss_sum   += batch_loss.item()
                        val_weight_sum += bw.sum().item()
                        
                val_loss = val_loss_sum / val_weight_sum if val_weight_sum > 0 else float('inf')

                pbar += 1                                
                
                if val_loss < best_loss:
                    best_loss        = val_loss
                    best_state       = copy.deepcopy ( model.state_dict() )
                    patience_counter = 0
                else:
                    patience_counter += 1
                    if patience <= patience_counter : break

        if best_state is not None:
            model.load_state_dict(best_state)

        model.eval()
        model.scaler = scaler

        val_preds = self._predict_single_model(model, X_val)
        return model, val_preds

    # =========================================================================
    ## Predict target probabilities using a trained PyTorch model.
    #  @param model Trained PyTorch DensityMLP instance.
    #  @param X Input features matrix.
    #  @return Array of predicted target probabilities p(y=1|x).
    def _predict_single_model(self, model, X):
        """ Predict target probabilities using a trained PyTorch model.
        """
        import torch
        from torch.utils.data import DataLoader, TensorDataset

        device = next(model.parameters()).device
        model.eval()

        X_scaled = model.scaler.transform(X)
        X_t = torch.tensor(X_scaled, dtype=torch.float32)

        dataset = TensorDataset(X_t)
        loader = DataLoader(dataset, batch_size=self.params['batch_size'], shuffle=False)

        probs = []
        with torch.no_grad():
            for (bx,) in loader:
                bx = bx.to(device)
                logits = model(bx)
                batch_probs = torch.sigmoid(logits).cpu().numpy()
                probs.append(batch_probs)

        probs = numpy.concatenate(probs, axis=0)
        return probs.astype(numpy.float32, copy=False)


# ==============================================================================
## @class GBReweighter
#  Helper wrapper class for reweighting using <code>hep_ml.reweight.GBReweighter</code>
#  by Alex Rogozhnikov 
class GBReweighter(Reweighter) :
    """ Helper wrapper class for reweighting using
    `hep_ml.reweight.GBReweighter` by Alex  Rogozhnikov 
    """

    # =========================================================================
    ## Initialize hep_ml GBReweighter wrapper.
    #  @param original Features array for original sample.
    #  @param target Features array for target sample.
    #  @param original_weight Initial weights for original sample (optional).
    #  @param target_weight Initial weights for target sample (optional).
    #  @param n_splits Number of cross-validation folds.
    #  @param store_original_weights If True, store computed original results.
    #  @param params Additional parameters for hep_ml GBReweighter.
    def __init__ ( self                   , * , 
                   original               ,
                   target                 ,
                   original_weight        = None  ,
                   target_weight          = None  , 
                   n_splits               = 5     ,
                   store_original_weights = True  , **params ) :
        """ Initialize hep_ml GBReweighter wrapper.
        """
        if not isinstance ( n_splits  , int ) : raise TypeError  ( "Invalid `n_splits' type %s" % typename ( n_splits ) )
        if not 0 <= n_splits <= 1000          : raise ValueError ( "Invalid `n_splits' value %s" % n_splits )     

        self.__n_splits = n_splits
        
        config = {            
            "n_estimators"      : 150       , 
            "learning_rate"     : 0.03      ,
            "max_depth"         : MAX_DEPTH ,
            "min_samples_leaf"  : 30        ,            
            "gb_args"           : {
                "subsample"     : 0.8    ,
                "max_features"  : 0.65   }
        }
        config.update ( params ) 

        check_all ( data1   = original          ,
                    data2   = target            ,
                    weight1 = original_weight   ,
                    weight2 = target_weight     ,
                    where   = typename ( self ) )   
        
        reg_case = self.needs_regularization ( original        = original        ,
                                               target          = target          ,
                                               original_weight = original_weight ,
                                               target_weight   = target_weight   )

        if reg_case :
            n_features = num_features ( original )
            neff_orig  = nEff ( original  , original_weight )
            neff_targ  = nEff ( target    , target_weight   ) 
            neff       = min  ( neff_orig , neff_targ       )
            
            config.update ( self.regularization ( config , n_features , neff ) )
            
            title = '%s strong regularization' % typename ( self ) 
            table = map2table_ex ( config    , 
                                   header    = ( 'Parameter' , 'type' , 'value' ) ,
                                   alignment = 'rcw'  , 
                                   prefix    = '# '   ,
                                   title     = title  )
            
            logger.info ( "%s is applied, case %s:\n%s" % ( title , reg_case , table ) )

        Reweighter.__init__ ( self            ,
                              original        = original        ,
                              target          = target          , 
                              original_weight = original_weight ,
                              target_weight   = target_weight   , **config )

        # =====================================================================
        try : # ===============================================================
            # =================================================================
            import sklearn.tree._classes as _cl
            if hasattr ( _cl , 'CRITERIA_REG' ) :
                if 'mse' in _cl.CRITERIA_REG and 'squared_error' not in _cl.CRITERIA_REG:
                    _cl.CRITERIA_REG [ 'squared_error' ] = _cl.CRITERIA_REG [ 'mse']
                elif 'squared_error' in _cl.CRITERIA_REG and 'mse' not in _cl.CRITERIA_REG:
                    _cl.CRITERIA_REG [ 'mse'] = _cl.CRITERIA_REG [ 'squared_error' ]
            # =================================================================
        except ( ImportError , AttributeError ) : # ===========================
            # =================================================================
            pass

        if not hasattr ( numpy , 'float' ) :
            logger.info ( 'No `numpy.float` found, aliasing `numpy.float64` as `numpy.float`')
            numpy.float = numpy.float64

        from hep_ml.reweight import GBReweighter as GBRW
        self.__reweighter = GBRW ( **self.params )
        
        if 1 < self.n_splits :
            from hep_ml.reweight import FoldingReweighter as FRW
            self.__reweighter = FRW ( self.reweighter ,
                                      n_folds      = self.n_splits     , 
                                      random_state = self.random_state ,
                                      verbose      = not self.silent  )
               
        with logAttention() if self.silent else NoContext() : 
            self.__reweighter.fit ( original        ,
                                    target          ,
                                    original_weight = original_weight , 
                                    target_weight   = target_weight   )

        self.__original_ratios             = None
        self.__original_reweighted_weights = None
        if store_original_weights :
            factors = self.__reweighter.predict_weights (
                original , 
                original_weight = original_weight )
            self.__original_ratios             = factors 
            self.__original_reweighted_weights = factors if weight_trivial ( original_weight ) else factors * original_weight

    # =========================================================================
    ## Check if strong regularization is needed for GBReweighter.
    #  @param original Features array for original sample.
    #  @param target Features array for target sample.
    #  @param original_weight Optional weights for original sample.
    #  @param target_weight Optional weights for target sample.
    #  @return String reason for regularization if needed, else empty string.
    def needs_regularization ( self ,
                               original                       ,
                               target                         ,
                               original_weight        = None  , 
                               target_weight          = None  ) :
        """ Check if strong regularization is needed for GBReweighter.
        """
        check_all ( data1   = original            ,
                    data2   = target              ,
                    weight1 = original_weight     ,
                    weight2 = target_weight       ,
                    where   = "%s:needs_regularization" % typename ( self ) ) 
        return RW_needs_regularization ( original        = original        ,
                                         target          = target          ,
                                         original_weight = original_weight ,
                                         target_weight   = target_weight   )
    
    # =========================================================================
    ## Apply regularization rules for GBReweighter.
    #  @param params Parameter dictionary.
    #  @param n_features Number of features.
    #  @param n_samples Effective sample size.
    #  @return Updated parameter dictionary.
    def regularization ( self       ,
                         params     , 
                         n_features ,
                         n_samples  ) : 
        """ Apply regularization rules for GBReweighter.
        """
        N = n_samples
        params [ 'n_estimators'     ] = 100
        params [ 'learning_rate'    ] = 0.05
        params [ 'max_depth'        ] = REG_DEPTH 
        params [ 'min_samples_leaf' ] = max ( 100, int ( N * 0.005 ) )
        
        gb_args = params.get ( 'gb_args', {} ).copy ()
        gb_args.update ( { 'subsample': 0.8, 'max_features': 0.8 } )
        params [ 'gb_args'          ] = gb_args
        return params
    
    # =========================================================================
    ## Get method identifier string.
    #  @return Method string identifier.
    @property
    def method ( self ) :
        """ Get method identifier string.
        """
        return method_GBRW
                        
    # =========================================================================
    ## Get number of cross-validation splits.
    #  @return Integer count of CV folds.
    @property
    def n_splits ( self ) :
        """`n_splits` : Get number of cross-validation splits.
        """
        return self.__n_splits
    
    # =========================================================================
    ## Alias property for n_splits.
    #  @return Integer count of CV folds.
    @property
    def n_folds ( self ) :
        """ `n_folds` : Alias property for n_splits.
        """
        return self.n_splits 
    
    # =========================================================================
    ## Get configuration dictionary.
    #  @return Configuration parameters dict.
    @property 
    def config ( self ) :
        """ Get configuration dictionary.
        """
        conf = {}
        conf.update ( super().config )
        conf [ 'n_splits' ] = self.n_splits
        return conf
    
    # =========================================================================
    ## Get density ratio factors r(x) for original sample.
    #  @return Array of density ratios.
    @property
    def original_ratios ( self ) :
        """ Get density ratio factors r(x) for original sample.
        """
        return self.__original_ratios
    
    # =========================================================================
    ## Get final reweighted weights for original sample.
    #  @return Array of reweighted weights.
    @property
    def original_reweighted_weights ( self ) : 
        """ Get final reweighted weights for original sample.
        """
        return self.__original_reweighted_weights
                
    # =========================================================================
    ## Get underlying hep_ml reweighter object.
    #  @return Internal GBReweighter or FoldingReweighter instance.
    @property
    def reweighter ( self ) :
        """ Get underlying hep_ml reweighter object.
        """
        return self.__reweighter
    
    # =========================================================================
    ## Compute reweighted event weights for new original data.
    #  @param original Features array for new dataset.
    #  @param original_weight Initial event weights (optional).
    #  @return Array of calculated reweighted weights.
    def weights ( self                   ,
                  original               ,
                  original_weight = None ) :
        """ Compute reweighted event weights for new original data.
        """
        if not valid_data_shape   ( original                   ) : raise TypeError ( "Invalid `original` type/shape: %s" % typename ( original ) )
        if not valid_weight       ( original_weight            ) : raise TypeError ( "Invalid `original_weight`!" )        
        if not compatible_weights ( original , original_weight ) : raise TypeError ( "Incompatible `original` data/weight!" )
        if self.n_features != num_features ( original )          : raise TypeError ( "Invalid #features!!")

        factors = self.reweighter.predict_weights ( original        = original        ,
                                                    original_weight = original_weight )
        return factors if weight_trivial ( original_weight ) else factors * original_weight

# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    from ostap.stats.tools import ( hasLightGBM ,
                                    hasXGBoost  ,
                                    hasCatBoost ,
                                    hasPyTorch  , 
                                    hasHepML    )

    if not hasLightGBM ( False ) : logger.warning  ( "No LightGBM available!" ) 
    if not hasXGBoost  ( False ) : logger.warning  ( "No XGBoost  available!" ) 
    if not hasCatBoost ( False ) : logger.warning  ( "No CatBoost available!" ) 
    if not hasPyTorch  ( False ) : logger.warning  ( "No PyTorch  available!" ) 
    if not hasHepML    ( False ) : logger.warning  ( "No HepMC    available!" ) 

# =============================================================================
#                                                                       The END 
# =============================================================================
