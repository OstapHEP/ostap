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
    'LightGBMDensityReweighter'      , # LightGBM-based density reweighter
    'XGBoostDensityReweighter'       , # XGBoost-based density reweighter
    'CatBoostDensityReweighter'      , # CatBoost-based density reweighter
    'PyTorchDensityReweighter'       , # PyTorch-based density reweighter
    'LogRegressionDensityReweighter' , # Density reweighter based on Logistic Regression
    'GBReweighter'                   , # Reweighter based on GBReweighter from hep_ml    
) 
# =============================================================================
from   ostap.core.ostap_types   import num_types 
from   ostap.utils.core         import typename
from   ostap.utils.basic        import numcpu, num_jobs, NoContext
from   ostap.logger.utils       import map2table_ex
from   ostap.logger.pretty      import nice_print 
from   ostap.utils.progress_bar import progress_bar
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
from   ostap.stats.tools      import hasSkLearn 
from   ostap.utils.memory     import memory
from   ostap.utils.timing     import timing 
import ostap.logger.symbols   as     S
import numpy, abc, math, gc  
# =============================================================================
# Logging setup
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighters' )
else                       : logger = getLogger( __name__ )
# =============================================================================
has_sklearn = hasSkLearn() 
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
DEFAULT_ESTIMATORS        = 1000  # High capacity: deep ensemble with low learning rate
REGULARIZED_ESTIMATORS    =  500  # Strong regularization: constrained number of trees
MAX_DEPTH                 =   10  # Deep tree capacity for high dimensions / large statistics
REGULARIZED_DEPTH         =    5  # Shallow depth (8 leaves max) to enforce heavy smoothing
LEARNING_RATE             =  0.02 # Small step for stable density ratio convergence
REGULARIZED_LEARNING_RATE =  0.04 # Slightly larger step for shallow regularized trees

LGBM_DEFAULT_ESTIMATORS        = 1600 # High capacity: deep ensemble with low learning rate
LGBM_REGULARIZED_ESTIMATORS    = 1000 # Strong regularization: constrained number of trees
LGBM_MAX_DEPTH                 =   10 # Deep tree capacity for high dimensions / large statistics
LGBM_REGULARIZED_DEPTH         =   10 # Shallow depth (8 leaves max) to enforce heavy smoothing
LGBM_LEARNING_RATE             = 0.01 # Small step for stable density ratio convergence
LGBM_REGULARIZED_LEARNING_RATE = 0.02 # Slightly larger step for shallow regularized trees

XGB_DEFAULT_ESTIMATORS         = 1600 # High capacity: deep ensemble with low learning rate
XGB_REGULARIZED_ESTIMATORS     = 1000 # Strong regularization: constrained number of trees
XGB_MAX_DEPTH                  =   10 # Deep tree capacity for high dimensions / large statistics
XGB_REGULARIZED_DEPTH          =   10 # Shallow depth (8 leaves max) to enforce heavy smoothing
XGB_LEARNING_RATE              = 0.01 # Small step for stable density ratio convergence
XGB_REGULARIZED_LEARNING_RATE  = 0.02 # Slightly larger step for shallow regularized trees

CATB_DEFAULT_ESTIMATORS        =  1600 # High capacity: deep ensemble with low learning rate
CATB_REGULARIZED_ESTIMATORS    =  1000 # Strong regularization: constrained number of trees
CATB_MAX_DEPTH                 =    10 # Deep tree capacity for high dimensions / large statistics
CATB_REGULARIZED_DEPTH         =    10 # Shallow depth (8 leaves max) to enforce heavy smoothing
CATB_LEARNING_RATE             = 0.005 # Small step for stable density ratio convergence
CATB_REGULARIZED_LEARNING_RATE = 0.010 # Slightly larger step for shallow regularized trees

## DEFAULT_ESTIMATORS        = 1000  # High capacity: deep ensemble with low learning rate
## REGULARIZED_ESTIMATORS    = 1000  # Strong regularization: constrained number of trees
## MAX_DEPTH                 =   10  # Deep tree capacity for high dimensions / large statistics
## REGULARIZED_DEPTH         =    5  # Shallow depth (8 leaves max) to enforce heavy smoothing

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

    neff_orig = float ( neff_orig )
    neff_targ = float ( neff_targ )
    
    # Two-sample effective size: Harmonic pooled nEff
    if neff_orig <= 0 or neff_targ <= 0 : neff = 0.0
    else : neff = 4.0 * ( neff_orig * neff_targ ) / ( neff_orig + neff_targ )

    # 3. Check weight efficiency ratios (nEff / nRaw)
    eff_orig = neff_orig / float ( nraw_orig ) if 0 < nraw_orig else 0.0
    eff_targ = neff_targ / float ( nraw_targ ) if 0 < nraw_targ else 0.0

    eff_orig = float ( eff_orig )
    eff_targ = float ( eff_targ )

    threshold0 = 5000.0  
    if eff_orig < 0.65 and neff_orig < threshold0 :
        effp   = eff_orig * 100.0
        v1     = nice_print ( neff_orig  , precision = 1 , width = 2 , with_sign = False )
        v2     = nice_print ( threshold0 , precision = 1 , width = 2 , with_sign = False )        
        result = 'eff_orig[%.0f%%]<65%%&neff_orig[%s]<%s' % ( effp , v1 , v2 )
        return result.replace ( ' ' , '' )
    
    if eff_targ < 0.65 and neff_targ < threshold0 :
        effp   = eff_targ * 100.0
        v1     = nice_print ( neff_targ  , precision = 1 , width = 2 , with_sign = False )
        v2     = nice_print ( threshold0 , precision = 1 , width = 2 , with_sign = False )        
        result = 'eff_targ[%.0f%%]<65%%&neff_targ[%s]<%s' % ( effp , v1 , v2 )
        return result.replace ( ' ' , '' )

    # 4. Low dimensionality (<= 4 features) needs regularization under limited statistics
    threshold = 50000.0 
    if nf <= 4 and neff < threshold :
        th     = nice_print ( threshold , precision = 1 , width = 2 , with_sign = False )        
        vv     = nice_print ( neff      , precision = 1 , width = 2 , with_sign = False )        
        result = 'nf[%d]<=4&neff[%s]<%s' % ( nf , vv , th ) 
        return result.replace ( ' ' , '' )
    
    # 5. Non-linear density threshold for multidimensional phase space growth
    required_stats = 1500.0 * ( nf ** 1.8 )
    if neff < required_stats :
        rs     = nice_print ( required_stats , precision = 1 , width = 2 , with_sign = False )
        vv     = nice_print ( neff           , precision = 1 , width = 2 , with_sign = False )        
        result = 'neff[%s]<%s' % ( vv , rs )
        return result.replace ( ' ' , '' )
        
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
        # ======================================================================
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

        self.__fitted_models       = {}
        self.__priors              = {}
        self.__target_weights_info = {}
        self.__norm_factor         = numpy.float32( 1.0 )
        self.__mode                = None
        self.__scale_factors       = {}
        

        ## n_features = num_features ( target ) 
        ## n_rows     = num_samples  ( target ) + num_samples  ( original )
        ## data_size  = n_rows * n_features 
        ## K          = 100_000         
        ## n_jobs     = num_jobs     ( params )
        ## n_jobs     = min ( n_jobs , max ( 1 , math.floor ( data_size / K ) ) )
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
            
            logger.info ( "%s is applied, case '%s':\n%s" % ( title , reg_case , table ) )

        self.__original_ratios             = None
        self.__original_reweighted_weights = None

        ## initialize the base 
        super ().__init__ ( original        = original        ,
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
                logger.error ( "%s: NaN/Inf found in original_ratios"     % typename ( self ) )
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
                logger.warning ( "%s: Large degradation of nEff %s %s %s" % ( self.method , n1 , S.arrow_right , n2 ) ) 

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
    

    @property 
    def scale_factors ( self ) :
        """`scale_factors` : scale factors for traing data(per stream)"""
        return self.__scale_factors 
        
        

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
            w_comb        = w_comb * ( len ( X_comb ) / ( numpy.sum ( w_comb ) + 1e-12 ) )
                        
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
            if 1 < raw_val_p.ndim : raw_val_p = raw_val_p[:, 1]
                
            oof_raw [ val_idx ] = raw_val_p.astype ( numpy.float32, copy = False )

        self.__fitted_models [ stream_key ] = stream_models 
        self.__scale_factors [ stream_key ] = scale_factor
        self.__priors        [ stream_key ] = numpy.float32 ( 1.0 )

        p_orig = oof_raw[ : len ( X_orig_sub ) ]
        eps    = 1.e-4 
        p_orig_clipped = numpy.clip ( p_orig, numpy.float32 ( eps ) , numpy.float32( 1.0 - eps  ) )
        
        ratios = p_orig_clipped / ( numpy.float32( 1.0 ) - p_orig_clipped )
        return ratios

    # =================================================================================
    ## Calculate average stream density ratios r(x) for unseen sample features.
    #  @param stream_key Stream identifier key.
    #  @param X Feature matrix to evaluate.
    #  @return Array of calculated stream density ratios.
    def __predict_stream_ratios( self, stream_key, X ):
        """ Calculate average stream density ratios r(x) for unseen sample features.
        """
        
        models = self.__fitted_models[ stream_key ]
        ratios_sum = numpy.zeros(len(X), dtype=numpy.float32)
        
        eps = 1.e-4 
        for model in models :
            p = self._predict_single_model ( model, X )
            if p.ndim > 1 :p = p[:, 1]
            
            p_clipped   = numpy.clip ( p, numpy.float32 ( eps ), numpy.float32 ( 1.0 - eps ) )
            ratios_sum += p_clipped / ( numpy.float32 ( 1.0 ) - p_clipped )

        return ratios_sum / len(models)

    # =================================================================================
    ## Default implementation to predict class probabilities from a single model.
    #  @param model Trained model instance.
    #  @param X Input feature array.
    #  @return Array of positive class probabilities p(y=1|x).
    def _predict_single_model ( self, model, X ):
        """ Predict class probabilities from a single model instance.
        """
        best_iter = getattr ( model, 'best_iteration_' , None ) or getattr( model, 'best_iteration', None )
        kwargs    = {}
        if best_iter is not None and 0 < best_iter : kwargs [ 'ntree_end' ] = best_iter            
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
            ratios  = self.__predict_stream_ratios ( "base", X_new_f32 )
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
                ratios_valid[ mask_pos ] = self.__predict_stream_ratios ( "orig_pos", X_valid[ mask_pos ] )
            if numpy.any( mask_neg ) :
                ratios_valid[ mask_neg ] = self.__predict_stream_ratios ( "orig_neg", X_valid[ mask_neg ] )

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
                r_pos_pos                = self.__predict_stream_ratios ( "pos_pos", X_p )
                r_pos_neg                = self.__predict_stream_ratios ( "pos_neg", X_p )
                ratios_valid[ mask_pos ] = ( r_pos_pos * w_pos - r_pos_neg * w_neg ) / w_total

            if numpy.any( mask_neg ) :
                X_n                      = X_valid[ mask_neg ]
                r_neg_pos                = self.__predict_stream_ratios ( "neg_pos", X_n )
                r_neg_neg                = self.__predict_stream_ratios ( "neg_neg", X_n )
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
        """ Initialize LightGBM density reweighter (High Capacity Mode Defaults).
        """
        config = {
            'objective'             : 'binary'                ,
            'metric'                : 'binary_logloss'        ,
            'n_estimators'          : LGBM_DEFAULT_ESTIMATORS ,
            'learning_rate'         : LGBM_LEARNING_RATE      ,
            'max_depth'             : LGBM_MAX_DEPTH          ,
            'num_leaves'            : min ( 1023 , 2 ** LGBM_MAX_DEPTH - 1 ) , 
            'max_bin'               : 2047                    , # 
            'min_child_samples'     : 10                      , # Slightly increased from 10 to stabilize leaf estimation
            'min_child_weight'      : 1e-3                    , # Increased from 1e-4 for numerical stability in ratios
            'min_split_gain'        : 0.0                     , # 
            'reg_alpha'             : 0.0                     , # L1 penalty to prune non-informative splits
            'reg_lambda'            : 0.0                     , # Increased L2 penalty (was 0.01) for smoother probabilities
            'subsample'             : 1.0                     , # Bagging fraction to reduce variance
            'subsample_freq'        : 1                       ,
            'colsample_bytree'      : 1.0                     , # Feature subsampling to increase ensemble diversity (was 1.0)
            'path_smooth'           : 1.0                     , # LightGBM tree smoothing for density ratio stability
            'boost_from_average'    : True                    ,
            'early_stopping_rounds' : LGBM_DEFAULT_ESTIMATORS // 2 ,
            'verbosity'             : -1                      ,
            'n_jobs'                : -1                      ,
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
        """ Apply strong regularization rules for limited statistics or low dimensions.
        """

                
        n_estimators = params.get ( 'n_estimators' , LGBM_REGULARIZED_ESTIMATORS )
        n_estimators = min        (  n_estimators  , LGBM_REGULARIZED_ESTIMATORS , 1 + math.floor ( n_samples / 2 ) ) 

        learning_rate                      = LGBM_REGULARIZED_LEARNING_RATE
               
        params [ 'n_estimators'          ] = n_estimators
        params [ 'learning_rate'         ] = learning_rate 
        params [ 'max_depth'             ] = LGBM_REGULARIZED_DEPTH

        num_leaves = min ( 1023 , 2 ** LGBM_REGULARIZED_DEPTH - 1    ) # maximl number of leaves for the given depth
        num_leaves = min ( num_leaves , 1 + math.floor ( n_samples ) )
        
        params [ 'num_leaves'            ] = num_leaves 

        params [ 'min_child_samples'     ] = max ( 10 , int ( n_samples * 0.001 ) ) # Scale with sample size (~0.5% min per leaf)
        params [ 'min_child_weight'      ] = 1.e-3 # Substantially increased from 1e-7 to prevent unstable leaves
        params [ 'min_split_gain'        ] = 0.0   # 
        params [ 'reg_alpha'             ] = 0.0   # Active L1 regularization
        params [ 'reg_lambda'            ] = 0.0   # Strong L2 regularization 
        params [ 'subsample'             ] = 1.0   # Active row subsampling (was incorrectly set to 1.0)
        params [ 'subsample_freq'        ] = 1
        params [ 'colsample_bytree'      ] = 1.0   # Active feature subsampling
        params [ 'min_data_in_bin'       ] = 1     # Increased from 1 to smooth out histogram binning
        
        raw_patience          = int ( 3   / learning_rate     )
        max_limit             = min ( 200 , n_estimators // 5 )
        early_stopping_rounds = max ( 15  , min ( raw_patience , max_limit ) )

        params [ 'early_stopping_rounds' ] = None ## 250 ## None ##  100 ## early_stopping_rounds 
        
        ## current_max_bin      = params.get ( 'max_bin', 4095 )
        ## params [ 'max_bin' ] = min ( current_max_bin , max ( 31 , 1 + math.floor ( n_samples / 5 ) ) ) 
        
        # Enforce tree-path smoothing parameter if supported by LightGBM
        ## params [ 'path_smooth' ] = 2.0

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
        import gc
        
        X_tr     = numpy.ascontiguousarray ( X_train , dtype = numpy.float32 )
        X_va     = numpy.ascontiguousarray ( X_val   , dtype = numpy.float32 )
        w_tr     = numpy.ascontiguousarray ( w_train , dtype = numpy.float32 ) if w_train is not None else None
        w_va     = numpy.ascontiguousarray ( w_val   , dtype = numpy.float32 ) if w_val   is not None else None

        trn_data = LightGBM.Dataset    ( X_tr    , label = y_train , weight = w_tr ,                        free_raw_data = True )
        val_data = LightGBM.Dataset    ( X_va    , label = y_val   , weight = w_va , reference = trn_data , free_raw_data = True )

        params = self.params.copy ()
        num_boost_round          = params.pop ( 'num_boost_round'       , None ) or params.pop ( 'n_estimators' , LGBM_DEFAULT_ESTIMATORS )
        early_stopping_rounds    = params.pop ( 'early_stopping_rounds' , None )
        params [ 'num_threads' ] = params.get ( 'num_threads'           , params.pop ( 'n_jobs' , max ( 1 , numcpu () // 2 ) ) ) 

        
        if early_stopping_rounds :
            callbacks  = [ LightGBM.early_stopping ( stopping_rounds = early_stopping_rounds , verbose = False ) ]
            valid_sets = [ val_data ]
        else : 
            callbacks  = []
            valid_sets = None 

        model = LightGBM.train ( params          = params          ,
                                 train_set       = trn_data        ,
                                 num_boost_round = num_boost_round ,
                                 valid_sets      = valid_sets      ,
                                 callbacks       = callbacks       )
            
        # =====================================================================
        ## PHOENIX
        # =====================================================================
        if True : # ===========================================================
            # =================================================================
            model_bytes = model.model_to_string()
            best_iter   = getattr ( model , 'best_iteration_' , None ) or getattr ( model , 'best_iteration' , None )                    
            model.free_dataset() 
            del model
            del trn_data
            del val_data
            del X_tr, X_va, w_tr, w_va
            
            model = LightGBM.Booster  ( model_str = model_bytes  )
            if not best_iter is None : model.best_iteration = best_iter 
            del model_bytes
            
            gc.collect ()                

            
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
        best_iter  = getattr ( model , 'best_iteration_' , None ) or getattr ( model , 'best_iteration' , None )        
        kwargs     = {}
        if best_iter is not None and 0 < best_iter : kwargs[ 'num_iteration' ] = best_iter            
        num_threads = self.params.get ( 'num_threads' , self.params.get ( 'n_jobs' , max ( 1 , numcpu () // 2 ) ) )                                    
        p = model.predict ( X , num_threads = num_threads , **kwargs )
        return p.astype( numpy.float32, copy = False )

# =============================================================================
## @class  XGBoostDensityReweighter
#  Density ratio reweighter using XGBoost as the underlying classifier.
class XGBoostDensityReweighter ( DensityReweighter ) :
    """ Density ratio reweighter using XGBoost as the underlying classifier.
    """

    # =========================================================================
    ## Initialize XGBoost density reweighter.
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
        """ Initialize XGBoost density reweighter (High Capacity Mode Defaults).
        """
        config = {
            'objective'             : 'binary:logistic'   ,
            'eval_metric'           : 'logloss'           ,
            'n_estimators'          : XGB_DEFAULT_ESTIMATORS  ,
            'learning_rate'         : XGB_LEARNING_RATE       ,
            'max_depth'             : XGB_MAX_DEPTH           ,
            'tree_method'           : 'hist'                  , # Fast histogram-based algorithm (similar to LGBM)
            'max_bin'               : 1023                    , # Limit binning to prevent noise fitting
            'min_child_weight'      : 1.e-3                   , # Minimum sum of instance weight (hessian) in a child
            'gamma'                 : 0.0                     , # Minimum loss reduction required for partition
            'alpha'                 : 0.0                     , # L1 penalty on leaf weights
            'lambda'                : 0.0                     , # L2 penalty on leaf weights
            'subsample'             : 0.8                     , # Row subsampling to reduce variance
            'colsample_bytree'      : 1.0                     , # Feature subsampling to increase ensemble diversity
            'early_stopping_rounds' : XGB_DEFAULT_ESTIMATORS // 2 ,            
            'verbosity'             :  0                      ,      
            'n_jobs'                : -1                      ,
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

    # =========================================================================
    ## Return the method identifier name.
    #  @return Method string identifier.
    @property
    def method ( self ) :
        """ Return the method identifier name.
        """
        return method_XGB 

    # =========================================================================
    ## Apply strong regularization rules for low statistics or low dimensions.
    #  @param params Current parameters dictionary.
    #  @param n_features Number of phase space features.
    #  @param n_samples Effective sample size.
    #  @return Updated parameters dictionary.
    def regularization ( self       ,
                         params     , 
                         n_features ,
                         n_samples  ) :
        """ Apply strong regularization rules for limited statistics or low dimensions.
        """
        
        n_estimators = params.get ( 'n_estimators' , XGB_REGULARIZED_ESTIMATORS )
        n_estimators = min        (  n_estimators  , XGB_REGULARIZED_ESTIMATORS , 1 + math.floor ( n_samples / 2 ) ) 

        learning_rate                      = XGB_REGULARIZED_LEARNING_RATE

        params [ 'max_depth'             ] = XGB_REGULARIZED_DEPTH       

        params [ 'n_estimators'          ] = n_estimators         
        params [ 'learning_rate'         ] = learning_rate
        
        # Scale min_child_weight (sum of hessian) with sample size to prevent isolated leaves
        ## params [ 'min_child_weight'      ] = max ( 10 , int ( n_samples * 0.001 ) ) 
        params [ 'min_child_weight'      ] = 1.e-3 
        
        params [ 'gamma'                 ] = 0.0    # More aggressive complexity control/pruning
        params [ 'alpha'                 ] = 0.0    # Active L1 regularization to encourage sparsity
        params [ 'lambda'                ] = 0.0    # Strong L2 regularization for smoother weights
        params [ 'subsample'             ] = 1.0    # Active row subsampling 
        params [ 'colsample_bytree'      ] = 1.0    # Active feature subsampling 
        params [ 'early_stopping_rounds' ] = None   # Tighter early stopping

        if n_samples <= 5000                                 :  params [ 'tree_method' ] = 'exact'
        elif 'hist' == params.get ( 'tree_method' , 'hist' ) :
            params [ 'tree_method' ] = 'hist'
            current_max_bin          = params.get ( 'max_bin' , 1023 )
            params [ 'max_bin'     ] = min ( current_max_bin , max ( 31 , int ( n_samples / 10 ) ) )
                    
        # Dynamic binning adjustment, only applicable if using histogram-based tree method
        ## if params.get ( 'tree_method' , 'hist' ) == 'hist' :
        ## current_max_bin      = params.get ( 'max_bin' , 255 )
        ## params [ 'max_bin' ] = min ( current_max_bin , max ( 31 , int ( n_samples / 20 ) ) )

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
        
        # Ensure contiguous C-style memory layout and explicit float32 data types
        X_tr = numpy.ascontiguousarray ( X_train , dtype = numpy.float32 )
        y_tr = numpy.ascontiguousarray ( y_train , dtype = numpy.float32 )
        w_tr = numpy.ascontiguousarray ( w_train , dtype = numpy.float32 ) if w_train is not None else None

        X_v  = numpy.ascontiguousarray ( X_val   , dtype = numpy.float32 )
        y_v  = numpy.ascontiguousarray ( y_val   , dtype = numpy.float32 )
        w_v  = numpy.ascontiguousarray ( w_val   , dtype = numpy.float32 ) if w_val is not None else None

        # Construct C++ DMatrix objects for training and validation
        dtrain = XGBoost.DMatrix ( X_tr , label = y_tr , weight = w_tr )
        dval   = XGBoost.DMatrix ( X_v  , label = y_v  , weight = w_v  )
        
        params = self.params.copy()
        params [ 'nthread' ] = params.get ( 'nthread' , params.pop ( 'n_jobs' , max ( 1 , numcpu () // 2 ) ) ) 
        
        n_estimators          = params.pop ( 'n_estimators'           , 1000 )
        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , None )

        evals = [ ( dval , 'val' ) ] if early_stopping_rounds is not None else []
        
        # Train model via low-level C++ API
        model = XGBoost.train (
            params,
            dtrain,
            num_boost_round      = n_estimators,
            evals                = evals,
            early_stopping_rounds= early_stopping_rounds,
            verbose_eval         = False
        )

        # =====================================================================
        ## PHOENIX
        # =====================================================================
        if True : # ===========================================================
            # =================================================================
            
            # Export tree structure into a lightweight binary/JSON byte array
            model_bytes = model.save_raw ( raw_format = "json" )
            best_iter   = getattr ( model , 'best_iteration_' , None ) or getattr ( model , 'best_iteration' , None )                    
            
            # Destroy heavy training DMatrix objects and training C++ booster instance
            del dtrain, dval, model,
            del X_tr, y_tr , w_tr, X_v , y_v, w_v 
            
            # Re-instantiate clean C++ booster dedicated strictly to inference
            model = XGBoost.Booster()
            model.load_model ( model_bytes )
            if not best_iter is None : model.best_iteration = best_iter
            
            del model_bytes            
            gc.collect()

        # Compute validation predictions safely using pre-allocated contiguous array
        val_preds = self._predict_single_model ( model , X_val )

        return model, val_preds
    
    # =========================================================================
    ## Predict probabilities using an XGBoost model.
    #  @param model Trained XGBoost Booster.
    #  @param X Input features array.
    #  @return Array of predicted probabilities.
    def _predict_single_model( self , model, X ):
        """ Predict probabilities using an XGBoost model.
        """
        import xgboost as XGBoost
        ## 
        dmat      = XGBoost.DMatrix ( X )
        best_iter = getattr ( model , 'best_iteration_' ,  None ) or getattr ( model , 'best_iteration' , None )
        kwargs    = {}
        if best_iter is not None and 0 < best_iter : kwargs [ 'iteration_range' ] = ( 0, best_iter + 1 )
        p         = model.predict ( dmat , **kwargs )
        return p.astype ( numpy.float32, copy = False )

# =============================================================================
## @class CatBoostDensityReweighter
#  Density ratio reweighter using CatBoost as the underlying classifier.
class CatBoostDensityReweighter(DensityReweighter):
    """ Density ratio reweighter using CatBoost as the underlying classifier.
    """
    
    # =========================================================================
    def __init__ ( self                   , * , 
                   original               ,
                   target                 ,
                   original_weight        = None ,
                   target_weight          = None ,
                   store_original_weights = True , **params ) :
        """ Initialize CatBoost density reweighter.
        """
        
        from catboost.utils import get_gpu_device_count
        gpu = get_gpu_device_count() 
        self.__gpu = gpu 
        
        config = { 'loss_function'         : 'Logloss'                        ,
                   'eval_metric'           : 'Logloss'                        ,
                   'iterations'            : CATB_DEFAULT_ESTIMATORS          ,
                   'learning_rate'         : CATB_LEARNING_RATE               ,
                   'depth'                 : CATB_MAX_DEPTH                   ,
                   'l2_leaf_reg'           : 0.5                              ,
                   'min_data_in_leaf'      : 5                                ,
                   'bootstrap_type'        : 'MVS'                            ,
                   'mvs_reg'               : 2.0                              , 
                   'subsample'             : 1.0                              ,
                   'random_strength'       : 0.01                             ,
                   'early_stopping_rounds' : CATB_DEFAULT_ESTIMATORS // 2     ,
                   'boosting_type'         : 'Plain'                          ,
                   'device'                : 'GPU' if 0 < self.gpu else 'CPU' , 
                   'verbose'               : False                            }
        
        config.update ( params )
        
        # Alias mappings to native CatBoost parameters
        if 'min_child_samples' in config : config [ 'min_data_in_leaf' ] = config.pop ( 'min_child_samples' )
        if 'n_jobs'            in config : config [ 'thread_count'     ] = config.pop ( 'n_jobs'            )
        if 'n_estimators'      in config : config [ 'iterations'       ] = config.pop ( 'n_estimators'      )
        
        if 'silent' in config :
            silent = config.get ( 'silent' , True ) 
            config [ 'verbose' ] = not silent 

        device = config.get ( 'device' , 'GPU' if 0 < self.gpu else 'CPU' )
        if self.gpu <= 0 : config [ 'device' ] = 'CPU'
        
        super().__init__( original               = original               ,
                          target                 = target                 ,
                          original_weight        = original_weight        ,
                          target_weight          = target_weight          ,
                          store_original_weights = store_original_weights , **config ) 
        
    # =========================================================================
    @property
    def gpu ( self ) :
        """`gpu` : result of `catboost.utils.get_gpu+_device_count"""
        return self.__gpu
    
    # =========================================================================
    @property
    def method(self):
        return method_CATB
    
    # =========================================================================
    def regularization ( self , params , n_features , n_samples ):
        """ Dynamic regularization rules for CatBoost.
        """
        iterations = params.pop ( 'n_estimators' , params.pop ( 'iterations' , CATB_REGULARIZED_ESTIMATORS ) )         
        ## iterations = min ( iterations , 1 + math.floor(n_samples / 2 ) ) 
        
        params [ 'iterations'            ] = iterations
        params [ 'learning_rate'         ] = CATB_REGULARIZED_LEARNING_RATE
        params [ 'depth'                 ] = CATB_REGULARIZED_DEPTH

        # the penalty smoothly increases to prevent overfitting.
        dynamic_l2 = max ( 3.0 , 10000.0 / ( n_samples + 1 ) ) 
        min_data   = max ( 1, int ( n_samples * 0.001 ) )
        
        params [ 'l2_leaf_reg'           ] = 0 ## dynamic_l2 
        params [ 'min_data_in_leaf'      ] = 1 ## min_data 
        
        ## params [ 'bootstrap_type'        ] = 'Bernoulli'
        ## params [ 'subsample'             ] = 1.0
        ## params [ 'random_strength'       ] = 0.001
        
        params [ 'early_stopping_rounds' ] = None
        params [ 'boosting_type'         ] = 'Plain'

        ## MVS - ОЧЕНЬ ХОРОЩО!!! 
        params [ 'bootstrap_type' ] = 'MVS'
        params [ 'mvs_reg'        ] =  1.0 

        ## НЕПЛОХО!
        ## params [ 'bootstrap_type' ] = 'Bernoulli'
        ## params [ 'subsample'      ] =  1.0

        return params


    
    # =========================================================================
    def _train_single_model ( self    ,
                              X_train , y_train , w_train ,
                              X_val   , y_val   , w_val   ) :
        """ Train single CatBoost model on memory fold data without forced contiguous copies.
        """

        import gc
        import catboost as CatBoost
        
        # Ensure contiguous C-style memory layout and explicit data types
        X_tr = numpy.ascontiguousarray ( X_train , dtype = numpy.float32  )
        y_tr = numpy.ascontiguousarray ( y_train , dtype = numpy.int32    )
        w_tr = numpy.ascontiguousarray ( w_train , dtype = numpy.float32  ) if w_train is not None else None
    
        X_v  = numpy.ascontiguousarray ( X_val   ,   dtype=numpy.float32  )
        y_v  = numpy.ascontiguousarray ( y_val   ,   dtype=numpy.int32    )
        w_v  = numpy.ascontiguousarray ( w_val   ,   dtype=numpy.float32  ) if w_val is not None else None
        
        # Instantiate heavy C++ Data Pools
        trn_pool = CatBoost.Pool ( X_tr , label = y_tr , weight = w_tr )
        val_pool = CatBoost.Pool ( X_v  , label = y_v  , weight = w_v  )
        
        params = self.params.copy()
        
        params [ 'thread_count'  ] = params.get ( 'thread_count' , params.pop ( 'n_jobs' , max ( 1 , numcpu () // 2 ) ) ) 
        params [ 'boosting_type' ] = params.get ( 'boosting_type' , 'Plain')
        
        iterations = params.pop ('n_estimators' , params.pop ( 'iterations' , CATB_DEFAULT_ESTIMATORS ) ) 
        if not isinstance(iterations, int) or not 10 <= iterations <= 10000:
            iterations = CATB_DEFAULT_ESTIMATORS
        
        early_stopping_rounds = params.pop('early_stopping_rounds', None)
        if isinstance(early_stopping_rounds, int) and 1 < early_stopping_rounds < iterations : pass
        else    : early_stopping_rounds = None
        
        params [ 'iterations' ] = iterations
        
        fit_kwargs = {}
        if early_stopping_rounds is not None:
            params     [ 'early_stopping_rounds' ] = early_stopping_rounds
            params     [ 'use_best_model'        ] = True
            fit_kwargs [ 'eval_set'              ] = val_pool
            
        params [ 'task_type' ] = params.pop ( 'device', 'GPU' if 0 < self.gpu else 'CPU' )
            
        model = CatBoost.CatBoostClassifier(**params)
        
        verbose_level = params.get ('verbose' , False)
        model.fit ( trn_pool , verbose = verbose_level, **fit_kwargs)

        # =====================================================================
        ## PHOENIX
        # =====================================================================
        if True : # ===========================================================
            # =================================================================
            import pickle 
                        
            # Extract tree geometry into a lightweight binary blob
            model_bytes = pickle.dumps ( model ) 
            best_iter   = getattr ( model , 'best_iteration_' , None ) or getattr ( model , 'best_iteration' , None )                    
            
            # Destroy training pools and heavy training model instance
            del trn_pool, val_pool, model, X_tr, y_tr, w_tr, X_v, y_v, w_v
            
            # Re-instantiate clean, standalone inference engine
            model = pickle.loads ( model_bytes )
            if not best_iter is None : model.best_iteration = best_iter
            
            del model_bytes 

            gc.collect()
            
        # Compute validation predictions using pre-allocated contiguous array
        val_preds = self._predict_single_model ( model,  X_val )
        
        return model, val_preds

    # =========================================================================
    def _predict_single_model ( self , model, X ) :
        """ Predict probabilities using a CatBoost model.
        """
        X_clean   = numpy.ascontiguousarray ( X , dtype = numpy.float32  )
        best_iter = getattr ( model , 'best_iteration_' , None ) or getattr ( model , 'best_iteration' , None )        
        kwargs    = {}
        if best_iter is not None and 0 <= best_iter : kwargs [ 'ntree_end' ] = best_iter + 1
        thread_count = self.params.get ( 'thread_count' , self.params.get ( 'n_jobs' , max ( 1 , numcpu () // 2 ) ) )                                    
        p = model.predict_proba ( X_clean , thread_count = thread_count ,  **kwargs)[:, 1]
        
        return p.astype (numpy.float32 , copy = False )


# =============================================================================
## @class PyTorchDensityReweighter
#  Optimized PyTorch MLP density ratio reweighter.
#  
#  This implementation utilizes index-based zero-copy batching, Automatic 
#  Mixed Precision (AMP), and optional torch.compile capabilities for 
#  maximum hardware utilization and memory safety.
# =============================================================================
class PyTorchDensityReweighter ( DensityReweighter ) :
    """ 
    Highly optimized PyTorch MLP density ratio reweighter.
    """

    # =========================================================================
    ## Constructor.
    #
    #  @param original        Original dataset.
    #  @param target          Target dataset.
    #  @param original_weight Optional weights for the original dataset.
    #  @param target_weight   Optional weights for the target dataset.
    #  @param progress        Boolean flag to show epoch progress bar.
    #  @param kwargs          Additional model hyperparameters.
    # =========================================================================
    def __init__ ( self            , * ,
                   original        ,
                   target          ,
                   original_weight = None ,
                   target_weight   = None ,
                   progress        = True , **kwargs ) :
                   
        import torch
        cuda = torch.cuda.is_available () 

        # ---------------------------------------------------------------------
        # Default configuration aligned for readability
        # ---------------------------------------------------------------------
        config = { 'hidden_dims'   : ( 512 , 512 , 512 , 512 ) ,
                   'dropout'       : 0.0                       ,               
                   'learning_rate' : 5e-4                      ,              
                   'weight_decay'  : 1e-5                      ,              
                   'batch_size'    : 2**16                     ,
                   'epochs'        : 2000                      ,              
                   'patience'      : 400                       ,
                   'min_epochs'    : 25                        ,
                   'min_delta'     : 2.e-4                     , 
                   'eval_freq'     : 5                         , 
                   'device'        : 'cuda' if cuda else 'cpu' ,
                   'compile'       : False                     , 
                   'n_jobs'        : -1                        ,
                 }

        config.update ( kwargs )
        self.__progress_epochs = config.pop ( 'progress' , True )
        
        device = config.get ( 'device' , 'cuda' if cuda else 'cpu' )        
        if not cuda : config [ 'device' ] = 'cpu'
        
        super().__init__ ( original        = original        ,
                           target          = target          ,
                           original_weight = original_weight ,
                           target_weight   = target_weight   ,
                           progress        = False           , **config )

        torch.set_num_threads ( 1  )

    # =========================================================================
    ## Returns the backend method identifier.
    # =========================================================================
    @property
    def method ( self ) :
        return method_TORCH 

    # =========================================================================
    ## Adjusts hyperparameters based on dataset size for regularization.
    # =========================================================================
    def regularization ( self , params , n_features , n_samples ) :
        """ Adjusts hyperparameters based on dataset size for regularization
        """
        params [ 'hidden_dims'   ] = ( 512 , 512 , 512 )
        params [ 'hidden_dims'   ] = ( 64  , 64  , 64  )
        params [ 'hidden_dims'   ] = ( 256 , 256 , 256 )
        params [ 'hidden_dims'   ] = ( 256 , 256 )
        params [ 'hidden_dims'   ] = ( 128 , 128 , 128 )
        params [ 'learning_rate' ] = 2e-4
        
        epochs                     = max ( 5000 , params.get ( 'epochs'     , 5000 ) )
        params [ 'epochs'        ] = epochs
        min_epochs                 = max ( 200  , params.get ( 'min_epochs' ,  200 ) )
        params [ 'min_epochs'    ] = min ( epochs   , min_epochs ) 
        patience                   = max ( 300 , params.get ( 'patience'    ,  300 ) )
        params [ 'patience'      ] = min ( patience , epochs ) 
        return params

    # =========================================================================
    ## Trains a single neural network model using Tabular ResNet architecture.
    #
    #  @param X_train Features for training.
    #  @param y_train Labels for training.
    #  @param w_train Sample weights for training.
    #  @param X_val   Features for validation.
    #  @param y_val   Labels for validation.
    #  @param w_val   Sample weights for validation.
    #  @return Tuple of ( trained_model , validation_predictions ).
    # =========================================================================
    def _train_single_model ( self , X_train , y_train , w_train , X_val , y_val , w_val ) :

        ## import os
        # ---------------------------------------------------------------------
        # Unlock CPU Multithreading before PyTorch backend initialization
        # Set to 4-8 threads depending on your batch system allocation
        # ---------------------------------------------------------------------
        ## n_threads = "8"
        ## os.environ [ 'OMP_NUM_THREADS'      ] = n_threads
        ## os.environ [ 'MKL_NUM_THREADS'      ] = n_threads
        ## os.environ [ 'OPENBLAS_NUM_THREADS' ] = n_threads

        import torch
        
        n_threads = self.params.get ( 'n_jobs' , max ( 1 , numcpu () // 2 ) )
        n_threads = min ( 4  , n_threads ,       max ( 1 , numcpu () // 2 ) )
        n_threads = max ( 1  , n_threads )

        torch.set_num_threads ( 8 )
        ### torch.set_num_interop_threads ( int ( n_threads ) )

        import torch.nn            as nn
        import torch.nn.functional as F
        from   sklearn.preprocessing import QuantileTransformer
        from   sklearn.preprocessing import MinMaxScaler

        import time

        device  = torch.device ( self.params.get ( 'device' , 'cuda' if torch.cuda.is_available () else 'cpu' ) )
        is_cuda = device.type == 'cuda'
        
        n_samples  = num_samples  ( X_train )
        n_features = num_features ( X_train )
        
        # =====================================================================
        # Extract and validate hyperparameters
        # =====================================================================
        batch_size = self.params.get ( 'batch_size' , 2**17 if is_cuda else 2**15 )
        if not isinstance ( batch_size , int ) or batch_size <= 0 :
            batch_size = 2**17 if is_cuda else 2**15
            
        batch_size = min ( max ( 2 , ( n_samples + 1 ) // 2 ) , batch_size ) 
        batch_size = max ( 2 , 2 ** math.floor ( math.log2 ( batch_size ) ) )

        nepochs       = self.params.get ( 'epochs'        , 2500                  )
        patience      = self.params.get ( 'patience'      , 400                   )
        min_epochs    = self.params.get ( 'min_epochs'    , 20                    )
        min_delta     = self.params.get ( 'min_delta'     , 2.e-4                 ) 
        eval_freq     = self.params.get ( 'eval_freq'     , 10                    )
        max_logit     = self.params.get ( 'max_logit'     , 10                    )
        learning_rate = self.params.get ( 'learning_rate' , 3e-3                  )
        hidden_dims   = self.params.get ( 'hidden_dims'   , ( 512 , 256 , 256 )   )
        weight_decay  = self.params.get ( 'weight_decay'  , 1.e-5                 )
        use_compile   = self.params.get ( 'compile'       , False                 ) and hasattr ( torch , 'compile' )
        
        # =====================================================================
        # Data preprocessing: QuantileTransformer maps complex kinematics
        # and resonance peaks into a smooth Gaussian space
        # =====================================================================
        n_quantiles = min ( 10000 , n_samples ) 
        scaler      = QuantileTransformer  ( output_distribution = 'normal' , n_quantiles = n_quantiles , random_state = self.random_state )
        ## scaler      = MinMaxScaler         ( feature_range = ( -1 , 1 ) )

        X_tr_scaled = scaler.fit_transform ( numpy.ascontiguousarray ( X_train , dtype = numpy.float32 ) )
        X_va_scaled = scaler.transform     ( numpy.ascontiguousarray ( X_val   , dtype = numpy.float32 ) )

        w_tr = numpy.ascontiguousarray ( w_train , dtype = numpy.float32 ) if w_train is not None else numpy.ones ( n_samples   , dtype = numpy.float32 )
        w_va = numpy.ascontiguousarray ( w_val   , dtype = numpy.float32 ) if w_val   is not None else numpy.ones ( len ( y_val ) , dtype = numpy.float32 )
        y_tr = y_train.astype ( numpy.float32 )

        # =====================================================================
        # Move data directly to device tensors
        # =====================================================================
        X_tr_t = torch.from_numpy ( X_tr_scaled ).to ( device )
        y_tr_t = torch.from_numpy ( y_tr        ).to ( device )
        w_tr_t = torch.from_numpy ( w_tr        ).to ( device )
        
        X_va_t = torch.from_numpy ( X_va_scaled ).to ( device )
        y_va_t = torch.from_numpy ( y_val.astype ( numpy.float32 ) ).to ( device )
        w_va_t = torch.from_numpy ( w_va        ).to ( device )

        class DensityRatioLoss(nn.Module):
            def __init__(self, max_logit=8.0):
                super().__init__()
                self.max_logit = max_logit

            def forward(self, logits, targets, weights):
                mask_orig = (targets == 0)
                mask_targ = (targets == 1)

                z_orig = logits[mask_orig]
                w_orig = weights[mask_orig]
                z_targ = logits[mask_targ]
                w_targ = weights[mask_targ]

                # Жесткий зажим логитов от взрыва экспоненты
                z_orig_safe = torch.clamp(z_orig, max=self.max_logit)

                # Вычисление с зажимом
                w_orig_sum = w_orig.abs().sum() + 1e-8
                loss_orig  = torch.sum(w_orig * torch.exp(z_orig_safe)) / w_orig_sum

                w_targ_sum = w_targ.abs().sum() + 1e-8
                loss_targ  = torch.sum(w_targ * z_targ) / w_targ_sum

                return loss_orig - loss_targ
    
        # =====================================================================
        ## @class ResBlock
        #  Inner class defining a Residual Block to prevent gradient vanishing
        #  and improve fitting of sharp kinematic edges.
        class ResBlock ( nn.Module ) :
            def __init__ ( self , dim ) :
                super().__init__()
                self.block = nn.Sequential (
                    nn.Linear    ( dim , dim ) ,
                    nn.LayerNorm ( dim )       ,
                    nn.SiLU      ()            ,
                    nn.Linear    ( dim , dim ) ,
                    nn.LayerNorm ( dim )
                )
                self.act = nn.SiLU ()

            def forward ( self , x ) :
                return self.act ( x + self.block ( x ) )

        # =====================================================================
        ## @class TabularResNet
        #  Inner class defining the network architecture with skip-connections.
        class TabularResNet ( nn.Module ) :
            def __init__ ( self , n_features , hidden_dims = ( 512 , 256 , 256 ) ) :
                super().__init__()
                layers = []
                
                first_dim = hidden_dims [ 0 ]
                layers.append ( nn.Linear ( n_features , first_dim ) )
                
                in_dim = first_dim
                for h_dim in hidden_dims :
                    if in_dim != h_dim :
                        layers.append ( nn.Linear ( in_dim , h_dim ) )
                        in_dim = h_dim
                    layers.append ( ResBlock ( h_dim ) )
                    
                layers.append ( nn.Linear ( in_dim , 1 ) )
                self.net = nn.Sequential ( *layers )
                
            def forward ( self , x ) : 
                return self.net ( x ).squeeze ( -1 )

        model_config = { 'n_features'  : n_features , 
                         'hidden_dims' : hidden_dims }
                         
        model = TabularResNet ( **model_config ).to ( device )
        
        if False and use_compile :
            try :
                model = torch.compile ( model )
            except Exception as e :
                logger.warning ( f"torch.compile failed, falling back to eager mode: {e}" )


        ## JIT?
        dummy_input = torch.randn ( batch_size , n_features , device = device )
        model       = torch.jit.trace ( model , dummy_input )

                
        optimizer  = torch.optim.AdamW                          ( model.parameters () , lr      = learning_rate , weight_decay =       weight_decay  )
        scheduler  = torch.optim.lr_scheduler.CosineAnnealingLR ( optimizer           , T_max   = nepochs       , eta_min      = 0.1 * learning_rate )
        scaler_amp = torch.amp.GradScaler                       ( 'cuda'              , enabled = is_cuda )

        w_va_sum         = w_va_t.abs().sum() + 1e-8
        best_loss        = float ( 'inf' )
        best_state       = None
        patience_counter = 0

        if not self.silent :
            data_size = n_features * n_samples * 4 ## in bytes 
            if   1024 ** 3 < data_size : data_size = '%.0fGB' % ( data_size / 1024 ** 3 )
            elif 1024 ** 2 < data_size : data_size = '%.0fMB' % ( data_size / 1024 ** 2 )
            elif 1024      < data_size : data_size = '%.0fkB' % ( data_size / 1024      )
            else                       : data_size = '%dB'    % ( data_size             )                                    
            steps_epoch = math.ceil  ( n_samples  / batch_size )
            logger.info ( "Train single model #events=%d total_data=%s batch_size=%d steps/epoch=%.0f" % ( n_samples   ,
                                                                                                           data_size   ,
                                                                                                           batch_size  ,
                                                                                                           steps_epoch ) )
            
        def get_autocast_config(device):
            if   'cuda' == device.type : return True, torch.bfloat16
            elif 'cpu'  == device.type :
                bf16_supported = hasattr ( torch.cpu, 'is_bf16_supported') and torch.cpu.is_bf16_supported ()
                return bf16_supported , torch.bfloat16
            return False, torch.float32

        use_autocast, autocast_dtype = get_autocast_config(device)

        criterion = DensityRatioLoss(max_logit=10.0)
        
        best_epoch = None 
        # =====================================================================
        # Main training loop
        # =====================================================================
        show_bar = self.__progress_epochs and not ( self.silent or self.progress )
        for epoch in progress_bar ( nepochs , silent = not show_bar , description = epoch_symbol ) :
            
            t0 = time.time ()
            
            model.train    ()

            perm = torch.randperm ( n_samples , device = device )

            for i in range ( 0 , n_samples , batch_size ) :
                idx = perm [ i : i + batch_size ]
                
                b_x = X_tr_t [ idx ]
                b_y = y_tr_t [ idx ]
                b_w = w_tr_t [ idx ]

                optimizer.zero_grad ( set_to_none = True )
                b_w_sum = b_w.abs().sum() + 1e-8
                
                with torch.autocast ( device_type = device.type, dtype = autocast_dtype , enabled = use_autocast ) :
                    
                    logits = model ( b_x )
                    ## loss   = F.binary_cross_entropy_with_logits ( logits , b_y , weight = b_w , reduction = 'sum' ) / b_w_sum
                    
                    loss    = criterion ( logits , b_y , b_w ) 

                    
                scaler_amp.scale  ( loss      ).backward ()
                scaler_amp.step   ( optimizer )
                scaler_amp.update ()

            # Step the learning rate scheduler
            scheduler.step ()

            # =================================================================
            # Validation
            # =================================================================
            if epoch % eval_freq == 0 or nepochs <= epoch + 10 : 
                
                model.eval ()
                val_loss = 0.0
                
                with torch.inference_mode () :
                    
                    for i in range ( 0 , len ( y_val ) , batch_size ) :
                        v_x = X_va_t [ i : i + batch_size ]
                        v_y = y_va_t [ i : i + batch_size ]
                        v_w = w_va_t [ i : i + batch_size ]
                        
                        ##  with torch.autocast ( device_type = 'cuda' if is_cuda else 'cpu' , enabled = is_cuda ) :
                        ##    v_logits     = model ( v_x )
                        ##    v_loss_batch = F.binary_cross_entropy_with_logits ( v_logits , v_y , weight = v_w , reduction = 'sum' )
                        ##    val_loss    += v_loss_batch.item ()
                        with torch.autocast ( device_type = 'cuda' if is_cuda else 'cpu' , enabled = is_cuda ) :
                            v_logits     = model ( v_x )
                            # Compute validation loss using custom density ratio loss
                            v_loss_batch = criterion ( v_logits , v_y , v_w )
                            val_loss    += v_loss_batch.item () * v_w.abs().sum().item()
                            
                    val_loss /= w_va_sum.item ()

                if not self.silent :
                    logger.info ( f"Epoch {epoch:03d}/{nepochs} | Val Loss: {val_loss:.5f} | LR: {scheduler.get_last_lr()[0]:.2e} | Time: {time.time()-t0:.3f} s" )

                if val_loss + min_delta < best_loss or best_state is None :
                    best_loss        = val_loss
                    best_epoch       = epoch 
                    unwrapped_model  = model._orig_mod if hasattr ( model , '_orig_mod' ) else model
                    best_state       = { k : v.clone().detach() for k , v in unwrapped_model.state_dict().items() }
                    patience_counter = 0
                else :
                    patience_counter += eval_freq
                    if min_epochs <= epoch and patience <= patience_counter :
                        if not self.silent : 
                            logger.info ( f"[Info] Early stopping at epoch {epoch}. Best Val Loss: {best_loss:.5f}" )
                            
                        if best_state is None :
                            unwrapped_model  = model._orig_mod if hasattr ( model , '_orig_mod' ) else model
                            best_state       = { k : v.clone().detach() for k , v in unwrapped_model.state_dict().items() }
                            
                        break

        # =====================================================================
        # Cleanup memory and restore best weights
        # =====================================================================
        del model , optimizer , scaler_amp , scheduler
        del X_tr_t , y_tr_t , w_tr_t , X_va_t , y_va_t , w_va_t
        
        gc.collect ()
        if is_cuda : torch.cuda.empty_cache ()

        model = TabularResNet ( **model_config ).to ( device )
        
        if best_state is not None :
            model.load_state_dict ( best_state )
            del best_state
            gc.collect ()

        model.eval ()
        model.scaler = scaler

        val_preds = self._predict_single_model ( model , X_val )
        return model , val_preds

    # =========================================================================
    ## Predicts probabilities for the given dataset.
    #
    #  @param model Trained PyTorch model.
    #  @param X     Features to predict on.
    #  @return Numpy array of predicted probabilities.
    # =========================================================================
    def _predict_single_model ( self , model , X ) :
        
        import torch
        
        model.eval ()
        device  = next ( model.parameters () ).device
        is_cuda = device.type == 'cuda'

        n_samples = num_samples ( X )
        
        batch_size = 2**17 if is_cuda else 2**16
        batch_size = max ( 1024 , 2 ** math.floor ( math.log2 ( max ( 1 , min ( batch_size , n_samples ) ) ) ) )

        X_scaled = model.scaler.transform ( numpy.ascontiguousarray ( X , dtype = numpy.float32 ) )
        preds    = numpy.empty ( n_samples , dtype = numpy.float32 )

        # =====================================================================
        # Batched inference loop
        # =====================================================================
        with torch.inference_mode () :
            for i in range ( 0 , n_samples , batch_size ) :
                chunk = X_scaled [ i : i + batch_size ]
                b_x   = torch.from_numpy ( chunk ).to ( device )
                
                with torch.autocast ( device_type = 'cuda' if is_cuda else 'cpu' , enabled = is_cuda ) :
                    logits  = model ( b_x )
                    b_preds = torch.sigmoid ( logits )
                
                preds [ i : i + batch_size ] = b_preds.cpu().numpy ()

        return preds
    
# =============================================================================
# Filter parameters to keep only those accepted by scikit-learn LogisticRegression
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
    """ Density ratio reweighter using Logistic Regression as the underlying classifier.
    """
    # =========================================================================
    ## Initialize Logistic Regression density reweighter.
    #  @param original Features array for original sample.
    #  @param target Features array for target sample.
    #  @param original_weight Initial weights for original sample (optional).
    #  @param target_weight Initial weights for target sample (optional).
    #  @param store_original_weights If True, store computed results for original sample.
    #  @param params Additional Logistic Regression parameters.
    def __init__( self , * , 
                  original               ,
                  target                 ,
                  original_weight        = None ,
                  target_weight          = None ,
                  store_original_weights = True , **params ) :        
        """ Initialize Logistic Regression density reweighter.
        """
        config = { 'C'             : 1.0     ,
                   'solver'        : 'lbfgs' ,
                   'max_iter'      : 2000    ,
                   'polynomials'   :    2    }
        
        config.update ( params )
        
        poly       = config.get   ( 'polynomials', 2 )
        n_features = num_features ( original         )
        
        if 15 <= n_features and 1 < poly :
            config [ 'polynomials' ] = 1
        
        super().__init__ ( original               = original               ,
                           target                 = target                 ,
                           original_weight        = original_weight        ,
                           target_weight          = target_weight          ,
                           store_original_weights = store_original_weights , **config )
        
    # =========================================================================
    @property
    def method ( self ) :
        return method_LR 

    # =========================================================================
    ## Dynamic regularization rules for Logistic Regression.
    def regularization ( self , params , n_features , n_samples ) :
        """ Dynamic regularization rules for Logistic Regression.
        """
        penalty = params.get ( 'penalty' , 'l2' )
        if penalty not in  ( None, 'none' ) :
            current_c = params.get ( 'C', 1.0 )
            params [ 'C' ] = max ( 0.1, min ( current_c, 0.5 ) ) 
        return params

    # =========================================================================
    ## Train single Logistic Regression model on fold data.
    def _train_single_model ( self    ,
                              X_train , y_train , w_train ,
                              X_val   , y_val   , w_val   ) :
        """ Train single Logistic Regression model on fold data.
        """
        from sklearn.linear_model  import LogisticRegression
        from sklearn.pipeline      import Pipeline        
        from sklearn.preprocessing import StandardScaler, PolynomialFeatures
        
        params = {}
        params.update ( self.params    )
        params.pop    ( 'n_jobs', None ) 

        poly = params.pop ( 'polynomials', 2 )
        n_features = num_features ( X_train )
        if 20 <= n_features and 1 < poly : poly = 1
        
        lr_kwargs = { k: v for k, v in params.items() if k in valid_LR_params }

        # Build the pipeline 
        steps = [ ( 'scaler1' , StandardScaler () ) ] 

        if poly and isinstance ( poly , int ) and 1 < poly <= 5 :
            steps += [ ( 'poly'    , PolynomialFeatures ( degree = poly  , include_bias = False ) ) ]
            steps += [ ( 'scaler2' , StandardScaler () ) ]
            
        lr_model = LogisticRegression ( **lr_kwargs )
        steps += [ ( 'logistic', lr_model ) ]
        
        model = Pipeline ( steps  )

        w_tr = numpy.ascontiguousarray ( w_train , dtype = numpy.float32 ) if w_train is not None else None

        if w_tr is not None : model.fit ( X_train , y_train , logistic__sample_weight = w_tr )
        else                : model.fit ( X_train , y_train )

        # =====================================================================
        ## PHOENIX
        # =====================================================================
        if True : # ===========================================================
            # =================================================================
            import pickle 
        
            model_bytes = pickle.dumps ( model ) 
            del model
            del w_tr
            
            model = pickle.loads ( model_bytes )
            del model_bytes 

            gc.collect()
        
        val_preds = self._predict_single_model ( model , X_val )
        return model , val_preds.astype ( numpy.float32 , copy = False )

    # =========================================================================
    ## Predict probabilities using a Logistic Regression model.
    def _predict_single_model ( self , model , X ) :
        """ Predict probabilities using a Logistic Regression model.
        """
        X_clean = numpy.ascontiguousarray ( X , dtype = numpy.float32 )
        p       = model.predict_proba ( X_clean ) [ : , 1 ]
        return p.astype ( numpy.float32 , copy = False )

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
            "n_estimators"      : 150           , 
            "learning_rate"     : LEARNING_RATE ,
            "max_depth"         : MAX_DEPTH     ,
            "min_samples_leaf"  : 30        ,            
            "gb_args"           : {
                "subsample"     : 1.0    ,
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
            
            logger.info ( "%s is applied, case '%s':\n%s" % ( title , reg_case , table ) )

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
                                      verbose      = not self.silent   )
               
        with logAttention() if self.silent else NoContext() : 
            self.__reweighter.fit ( original        ,
                                    target          ,
                                    original_weight = original_weight , 
                                    target_weight   = target_weight   )

        # =====================================================================
        ## PHOENIX
        # =====================================================================
        if True : # ===========================================================
            # =================================================================
            import pickle
            
            model_bytes       = pickle.dumps ( self.__reweighter ) 
            del self.__reweighter
            self.__reweighter = pickle.loads ( model_bytes )
            del model_bytes 
            
            gc.collect()

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
        params [ 'learning_rate'    ] = REGULARIZED_LEARNING_RATE
        params [ 'max_depth'        ] = REGULARIZED_DEPTH 
        params [ 'min_samples_leaf' ] = max ( 100, int ( N * 0.005 ) )
        
        gb_args = params.get ( 'gb_args', {} ).copy ()
        gb_args.update ( { 'subsample': 1.0, 'max_features': 0.8 } )
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
                                    hasSkLearn  , 
                                    hasHepML    )

    if not hasLightGBM ( False ) : logger.warning  ( "No LightGBM available!" ) 
    if not hasXGBoost  ( False ) : logger.warning  ( "No XGBoost  available!" ) 
    if not hasCatBoost ( False ) : logger.warning  ( "No CatBoost available!" ) 
    if not hasPyTorch  ( False ) : logger.warning  ( "No PyTorch  available!" ) 
    if not hasSkLearn  ( False ) : logger.warning  ( "No SkLearn  available!" ) 
    if not hasHepML    ( False ) : logger.warning  ( "No HepML    available!" ) 

# =============================================================================
#                                                                       The END 
# =============================================================================
