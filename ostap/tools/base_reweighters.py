#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/base_reweighters.py
#  Base density ratio reweighters and iterative weight boosting implementations
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @Date  2026-07-19
# =============================================================================
""" Base density ratio reweighters and iterative weight boosting classes
    (compatible with signed weights and Ostap infrastructure).
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__all__     = (
    'DensityReweighter'        , # Abstract base class for density reweighters
    'IterativeReweighter'      , # Base class for iterative reweighters
) 
# =============================================================================
from   ostap.core.ostap_types   import num_types, string_types  
from   ostap.utils.core         import typename
from   ostap.utils.basic        import num_jobs
from   ostap.logger.pretty      import nice_print 
from   ostap.utils.progress_bar import progress_bar
from   ostap.tools.reweighter   import Reweighter
from   ostap.stats.counters     import SE, table_counters  
from   ostap.stats.utils        import ( valid_weight       ,
                                         valid_data_shape   ,
                                         compatible_weights , 
                                         num_features       ,
                                         check_all          ,
                                         nEff               ,
                                         truncated_mean     )
import ostap.logger.symbols   as     S
import numpy, abc, re  
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.tools.base_reweighters' )
else                       : logger = getLogger ( __name__ )
# =============================================================================
## @class Voting
#  @brief Encapsulates fold prediction aggregation (voting) strategies.
#  @details Acts as a callable strategy pattern object that reduces a stacked fold 
#           predictions matrix of shape (N_folds, N_samples) into an aggregated array 
#           of shape (N_samples,) using predefined statistical reduction rules or custom callables,
#           including arithmetic, geometric, median, and truncated mean aggregation.
class Voting ( object ) :
    """ Encapsulates fold prediction aggregation (voting) strategies. """
    def __init__ ( self, strategy = 'geometric' ) :
        """ Initialize voting strategy.
            @param strategy Aggregation strategy name (e.g. 'arithmetic', 'mean', 'median', 'geometric', 
                            'trim:1', 'truncated = 2') or a custom callable expecting `(array, axis)`.
        """
        # Default trim count if not specified explicitly
        self.__trim_low  = 1
        self.__trim_high = None
        
        if callable ( strategy ) :
            
            self.__func = strategy
            self.__name = getattr ( strategy, '__name__', 'custom' )
            
            return

        if not isinstance ( strategy , string_types ) :
            raise TypeError ( "Invalid `voting` type: %s" % type ( strategy ) )
    
        
        v_lower = strategy.lower ().strip ()
        
        if   v_lower in ( 'arithmetic', 'mean' ) : 
            self.__func = numpy.mean
            self.__name = v_lower
        elif v_lower == 'median' : 
            self.__func = numpy.median
            self.__name = v_lower
        elif v_lower == 'geometric' :
            # Compute geometric mean across fold predictions in log-space to prevent numerical underflow
            self.__func = lambda x, axis=0 : numpy.exp ( 
                numpy.mean ( numpy.log ( numpy.clip ( x, 1e-12, 1e+12 ) ) , axis = axis ) 
            )
            self.__name = v_lower
                    
        else :
            
            # Match pattern: optional spaces, trim or truncate(d), optional spaces, optional (: or =), digit
            match = re.match ( r'^\s*(trim|truncate|truncated)\s*(?:([:=])\s*(\d+))?\s*$', strategy, re.IGNORECASE )
            if not match : raise ValueError ( f"Invalid `voting` strategy string: {strategy}" )
            
            digit_str = match.group ( 3 )
            if digit_str is not None :
                try :
                    self.__trim_low = int ( digit_str )
                except ValueError :
                    raise ValueError ( f"Invalid trim value in voting strategy string: {strategy}" )
                
                # Bind the standalone truncated_mean function
                self.__func = lambda x, axis = 0 : truncated_mean ( 
                    x, axis = axis, trim_low = self.__trim_low, trim_high = self.__trim_high 
                )
                self.__name = v_lower

    ## @property name
    #  @brief String identifier of the configured voting strategy.
    @property
    def name ( self ) :
        """ `name`: String identifier of the configured voting strategy. """
        return self.__name

    ## @brief Apply the voting aggregation function to stacked fold predictions.
    #  @param x Matrix of fold predictions with shape (N_folds, N_samples).
    #  @param axis Dimension along which to aggregate predictions across folds.
    #  @return Aggregated prediction array of shape (N_samples,).
    def __call__ ( self, x, axis = 0 ) :
        """ Apply the voting aggregation function to stacked fold predictions. """
        return self.__func ( x, axis = axis )

# ===============================================================================
## @class DensityReweighter
#  @brief Abstract base class for immutable, adaptive density-ratio reweighting.
#  @details Implements 1-, 2-, and 4-stream signed-measure density estimations.
#           Executes training immediately upon instantiation and seals resulting weights.
#           During initialization, K-Fold cross-validation splits the data to train N independent 
#           models on (N-1)/N fractions. Initial weights for the training original events are 
#           computed strictly using Out-Of-Fold (OOF) predictions from their respective fold models.
#           Subsequent calls to `weights()` on unseen data leverage all N trained models 
#           aggregated via the configured `Voting` strategy.
class DensityReweighter ( Reweighter, abc.ABC ) :
    """ Abstract base class for immutable, adaptive density-ratio reweighting. """

    ## @brief Initialize and fit the density ratio reweighter ensemble.
    #  @param original Source feature matrix.
    #  @param target Target feature matrix.
    #  @param original_weight Source observation weights (supports signed weights).
    #  @param target_weight Target observation weights (supports signed weights).
    #  @param clip_threshold Upper threshold for clipping extreme density ratios.
    #  @param n_splits Number of cross-validation splits (K-Fold).
    #  @param voting Aggregation strategy for combining fold models during prediction.
    #  @param oof_mode Flag controlling initial evaluation strategy: 
    #                  True  - strict Out-Of-Fold evaluation (each sample evaluated by its validation fold model);
    #                  False - training-ensemble smoothing (each sample averaged over all N-1 models where it was in train set).
    #  @param store_original_weights Flag to retain fitted original ratio arrays.
    #  @param progress Enable progress bar display during fold fitting.
    #  @param params Additional parameters passed to the parent class.
    def __init__ ( self                                  , * ,
                   original                              ,
                   target                                ,
                   original_weight        = None         , 
                   target_weight          = None         ,
                   clip_threshold         = 1.e+4        ,
                   n_splits               = 5            ,
                   voting                 = 'arithmetic' , 
                   oof_mode               = True         , 
                   store_original_weights = False        ,
                   progress               = True         , **params ) :
        """ Initialize and fit the density ratio reweighter ensemble. """
        
        if not isinstance ( n_splits, int ) : raise TypeError  ( "Invalid `n_splits' type %s"  % typename ( n_splits ) )
        if not 0 <= n_splits <= 1000        : raise ValueError ( "Invalid `n_splits' value %s" % n_splits )
        if not isinstance ( clip_threshold, num_types ) :
            raise TypeError ( "Invalid `clip_threshold' type %s" % typename ( clip_threshold ) )
        if not 0 < clip_threshold :
            raise ValueError ( "Invalid `clip_threshold' value %s" % clip_threshold )
        
        check_all ( data1   = original          ,
                    data2   = target            ,
                    weight1 = original_weight   ,
                    weight2 = target_weight     ,
                    where   = typename ( self ) )   
        
        self.__progress            = True if progress else False 
        self.__clip_threshold      = float ( clip_threshold )
        self.__n_splits            = n_splits
        self.__voting              = Voting ( voting ) 
        self.__oof_mode            = True if oof_mode else False
        
        # Public attributes for clean and direct access in subclasses
        self.fitted_models         = {}
        self.priors                = {}
        self.__target_weights_info = {}
        self._norm_factor          = numpy.float32 ( 1.0 )
        self.__mode                = None
        self.scale_factors         = {}
        
        params [ 'n_jobs' ]        = num_jobs ( params )
                
        self.__original_reweighted_weights = None

        super ().__init__ ( original        = original        ,
                            target          = target          , 
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **params )

        # Execute immediate fitting across positive streams and compute OOF initial weights
        original_ratios, final_weight =  self.__fit_and_compute (
            original.astype        ( numpy.float32 , copy = False ) ,
            original_weight.astype ( numpy.float32 , copy = False ) if original_weight is not None else None , 
            target.astype          ( numpy.float32 , copy = False ) ,
            target_weight.astype   ( numpy.float32 , copy = False ) if target_weight   is not None  else None )
                
        ## final weights 
        final_weight  = self.preserve_mass ( original , original_weight , final_weight ) 

        if not numpy.all ( numpy.isfinite ( original_ratios ) ) : logger.error ( "%s: NaN/Inf found in original_ratios" % typename ( self ) )
        if not numpy.all ( numpy.isfinite ( final_weight    ) ) : logger.error ( "%s: NaN/Inf found in final_weight"    % typename ( self ) )
        
        r_min = float ( numpy.min ( original_ratios ) )
        r_max = float ( numpy.max ( original_ratios ) )
        r_std = float ( numpy.std ( original_ratios ) )
        
        if numpy.isclose ( r_min, r_max, atol = 1e-5 ) or r_std < 1e-6 :
            r1 = nice_print ( r_min , precision = 4 , with_sign = False , strip = True )
            r2 = nice_print ( r_std , precision = 4 , with_sign = False , strip = True )
            logger.warning ( "%s: All ratios are constant r(x) = %s (std=%s)" % ( self.method , r1 , r2 ) ) 
            
        neff_before = nEff ( original , original_weight )
        neff_after  = nEff ( original , final_weight    )
        if neff_after < 0.10 * neff_before :
            nb = nice_print ( neff_before , precision = 2 , with_sign = False , strip = True )
            na = nice_print ( neff_after  , precision = 2 , with_sign = False , strip = True )
            logger.warning ( "%s: Huge degradation of #eff %s %s %s" % ( self.method , nb , S.arrow_right , na ) ) 
            
        clipped_count = numpy.sum ( original_ratios >= ( self.__clip_threshold * 0.999 ) )
        if 0 < clipped_count : 
            frac = 100.0 * clipped_count / len ( original_ratios ) 
            if 1 < frac : logger.warning ( "%s: Too many clipped events %.2f[%%]" % ( self.method , frac ) )

        if store_original_weights : 
            self.__original_ratios             = original_ratios 
            self.__original_reweighted_weights = final_weight 
        
        ## As a last step: produce summary of weights 
        if self.verbose : self.weight_summary ( original        ,
                                                original_weight ,
                                                final_weight    ) 

    ## @property n_splits
    #  @brief Number of cross-validation folds.
    @property
    def n_splits ( self ) :
        """ `n_splits`: Number of cross-validation folds. """
        return self.__n_splits

    ## @property voting
    #  @brief Configured voting strategy helper instance.
    @property 
    def voting ( self ) :
        """ `voting`: Configured voting strategy helper instance. """
        return self.__voting
    
    ## @property oof_mode 
    #  @brief OOF-mode
    @property 
    def oof_mode( self ) :
        """ `oof_mode`: OOF-mode"""
        return self.__oof_mode 
    
    ## @property original_ratios
    #  @brief Out-of-fold ratios for training original dataset.
    @property
    def original_ratios ( self ) :
        """ `original_ratios`: Out-of-fold ratios for training original dataset. """
        return self.__original_ratios

    ## @property original_reweighted_weights
    #  @brief Final reweighted weight values for training original dataset.
    @property
    def original_reweighted_weights ( self ) :
        """ `original_reweighted_weights`: Final reweighted weight values for training original dataset. """
        return self.__original_reweighted_weights
    
    ## @property norm_factor
    #  @brief Global normalization factor conserving total reweighted mass.
    @property
    def norm_factor ( self ) :
        """ `norm_factor`: Global normalization factor conserving total reweighted mass. """
        return self._norm_factor
        
    ## @property mode
    #  @brief Partition mode ('1-stream', '2-stream_orig', '2-stream_targ', '4-stream').
    @property
    def mode ( self ) :
        """ `mode`: Partition mode ('1-stream', '2-stream_orig', '2-stream_targ', '4-stream'). """
        return self.__mode

    ## @property progress
    #  @brief Progress bar status flag.
    @property
    def progress ( self ) :
        """ `progress`: Progress bar status flag. """
        return self.__progress

    @property
    def method ( self ) :
        """ Return method identifier name. """
        return '%s DRW' % self.symbol 
    
    ## @property config
    #  @brief Serialization configuration parameters dictionary.
    @property
    def config ( self ) :
        """ `config`: Serialization configuration parameters dictionary. """
        conf = {}
        conf.update ( super ().config )
        conf [ "progress" ]                    = self.progress
        conf [ "mode"     ]                    = self.mode
        conf [ "clip_threshold"  ]             = self.__clip_threshold
        conf [ "voting"          ]             = self.__voting.name 
        conf [ "n_splits"        ]             = self.__n_splits
        conf [ "oof_mode"        ]             = self.__oof_mode
        conf [ "original_reweighted_weights" ] = self.__original_reweighted_weights is not None
        
        return conf

    ## @brief Abstract method to train a single backend estimator on fold data.
    #  @param X_train Training features matrix for current fold (N-1)/N.
    #  @param y_train Binary labels for training fold (0=original, 1=target).
    #  @param w_train Training sample weights.
    #  @return Trained model instance.
    @abc.abstractmethod
    def _train_single_model ( self    ,
                              X_train ,
                              y_train ,
                              w_train ) :
        """ Abstract method to train a single backend estimator on fold data. """
        raise NotImplementedError

    ## @brief Abstract method to return probability predictions p(x) for a single model.
    #  @param model Backend estimator instance.
    #  @param X Features matrix.
    #  @return Predicted probabilities array p(x).
    @abc.abstractmethod
    def _predict_single_model ( self  ,
                                model ,
                                X     ) :
        """ Abstract method to return probability predictions p(x) for a single model. """
        raise NotImplementedError

    ## @brief Fit streams based on sign decomposition of original and target weights.
    #  @details Decomposes signed distributions into non-negative streams:
    #           - 1-stream: both original and target weights are non-negative.
    #           - 2-stream_orig: original has signed weights; target is non-negative.
    #           - 2-stream_targ: original is non-negative; target has signed weights.
    #           - 4-stream: both original and target contain signed weights.
    #  @param original Source feature matrix.
    #  @param original_weight Source observation weights.
    #  @param target Target feature matrix.
    #  @param target_weight Target observation weights.
    #  @return Tuple of (original_ratios, original_reweighted_weights).
    def __fit_and_compute ( self            ,
                            original        ,
                            original_weight ,
                            target          ,
                            target_weight   ) :
        """ Fit streams based on sign decomposition of original and target weights. """
        
        if original_weight is None and target_weight is None :
            
            self.__mode  = "1-stream"
            ratios_valid = self._fit_eval_stream ( "base", original, None, target, None )

            clipped           = numpy.clip ( ratios_valid,
                                             numpy.float32 ( 0.0 ) ,
                                             numpy.float32 ( self.__clip_threshold ) )
            self._norm_factor = numpy.float32 ( len ( original ) / ( numpy.sum ( clipped, dtype = numpy.float64 ) + 1e-12 ) )

            original_ratios = self._apply_norm_factor ( clipped ) 
            return original_ratios, original_ratios.copy ()

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
            ratios_valid = self._fit_eval_stream ( "base"       ,
                                                   X_orig_valid ,
                                                   w_orig_valid ,
                                                   X_targ_valid ,
                                                   w_targ_valid )

        elif has_neg_orig and not has_neg_targ :
            
            self.__mode                   = "2-stream_orig"
            ratios_valid                  = numpy.zeros ( len ( X_orig_valid ) , dtype = numpy.float32 )
            ratios_valid[ mask_orig_pos ] = self._fit_eval_stream ( "orig_pos"                     ,
                                                                    X_orig_valid [ mask_orig_pos ] ,
                                                                    w_orig_valid [ mask_orig_pos ] ,
                                                                    X_targ_valid                   ,
                                                                    w_targ_valid                   )
            ratios_valid[ mask_orig_neg ] = self._fit_eval_stream ( "orig_neg"                     ,
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

            r_pos_pos = self._fit_eval_stream ( "pos_pos"    ,
                                                X_orig_valid ,
                                                w_orig_valid ,
                                                X_targ_pos   ,
                                                w_targ_pos   )
            r_pos_neg = self._fit_eval_stream ( "pos_neg"    ,
                                                X_orig_valid ,
                                                w_orig_valid ,
                                                X_targ_neg   ,
                                                w_targ_neg   )

            ratios_valid = ( r_pos_pos * w_pos_sum - r_pos_neg * w_neg_sum ) / ( w_total + numpy.float32 ( 1e-12 ) )

        else :
            
            self.__mode            = "4-stream"
            X_orig_pos , w_orig_pos = X_orig_valid [ mask_orig_pos ] , w_orig_valid [ mask_orig_pos ]
            X_orig_neg , w_orig_neg = X_orig_valid [ mask_orig_neg ] , w_orig_valid [ mask_orig_neg ]
            X_targ_pos , w_targ_pos = X_targ_valid [ mask_targ_pos ] , w_targ_valid [ mask_targ_pos ]
            X_targ_neg , w_targ_neg = X_targ_valid [ mask_targ_neg ] , w_targ_valid [ mask_targ_neg ]

            w_pos_sum                  = numpy.float32 ( numpy.sum ( w_targ_pos, dtype = numpy.float64 ) )
            w_neg_sum                  = numpy.float32 ( numpy.abs ( numpy.sum ( w_targ_neg, dtype = numpy.float64 ) ) ) 
            self.__target_weights_info = { "W_pos": w_pos_sum, "W_neg": w_neg_sum }
            w_total                    = w_pos_sum + w_neg_sum

            r_pos_pos = self._fit_eval_stream ( "pos_pos", X_orig_pos, w_orig_pos, X_targ_pos, w_targ_pos )
            r_pos_neg = self._fit_eval_stream ( "pos_neg", X_orig_pos, w_orig_pos, X_targ_neg, w_targ_neg )
            r_neg_pos = self._fit_eval_stream ( "neg_pos", X_orig_neg, w_orig_neg, X_targ_pos, w_targ_pos )
            r_neg_neg = self._fit_eval_stream ( "neg_neg", X_orig_neg, w_orig_neg, X_targ_neg, w_targ_neg )

            ratios_valid                  = numpy.zeros ( len ( X_orig_valid ) , dtype = numpy.float32 )
            denom                         = w_total + numpy.float32 ( 1e-12 )
            ratios_valid[ mask_orig_pos ] = ( r_pos_pos * w_pos_sum - r_pos_neg * w_neg_sum ) / denom
            ratios_valid[ mask_orig_neg ] = ( r_neg_pos * w_pos_sum - r_neg_neg * w_neg_sum ) / denom

        ratios_valid_norm = self._normalize_and_clip ( ratios_valid, w_orig_valid )

        original_ratios                 = numpy.zeros ( len ( w_orig ), dtype = numpy.float32 )
        original_ratios[ nz_orig_mask ] = ratios_valid_norm
        original_reweighted_weights     = w_orig * original_ratios


        # --- DEBUG PRINTS ---
        print("DEBUG __fit_and_compute:")
        print("  len(w_orig):", len(w_orig))
        print("  len(nz_orig_mask):", len(nz_orig_mask), "True count:", numpy.sum(nz_orig_mask))
        print("  len(ratios_valid):", len(ratios_valid))
        print("  len(ratios_valid_norm):", len(ratios_valid_norm))
        print("  ratios_valid_norm min/max:", numpy.min(ratios_valid_norm), numpy.max(ratios_valid_norm))
        print("  original_ratios min/max:", numpy.min(original_ratios), numpy.max(original_ratios))
        print("----------------------")
        
        return original_ratios, original_reweighted_weights

    ## @brief Execute K-Fold cross-validation for a specific non-negative stream.
    #  @details Supports two evaluation modes depending on `oof_mode`:
    #           - Strict OOF (`oof_mode=True`): Each sample gets a prediction strictly from its single 
    #             validation fold model. Pure methodology, but prone to variance spikes on sparse distribution tails.
    #           - Training Ensemble (`oof_mode=False`): Each sample is evaluated by averaging across all 
    #             N-1 models for which it served as a training instance, bypassing single-model tail extrapolation spikes.
    #  @param stream_key Identifier for the stream (e.g. 'orig_pos').
    #  @param X_orig_sub Original feature subset.
    #  @param w_orig_sub Original weight subset.
    #  @param X_targ_sub Target feature subset.
    #  @param w_targ_sub Target weight subset.
    #  @return Array of predicted density ratios.
    def _fit_eval_stream ( self       ,
                           stream_key ,
                           X_orig_sub ,
                           w_orig_sub ,
                           X_targ_sub ,
                           w_targ_sub ) :
        """ Execute K-Fold cross-validation with switchable OOF / training-ensemble mode. """

        # Combine original and target feature matrices into a unified cross-validation dataset
        X_comb = numpy.vstack ( [ X_orig_sub, X_targ_sub ] )
        y_comb = numpy.hstack ( [ numpy.zeros ( len ( X_orig_sub ) , dtype = numpy.float32 ) ,
                                  numpy.ones  ( len ( X_targ_sub ) , dtype = numpy.float32 ) ] )

        w_orig_abs = numpy.abs ( w_orig_sub ) if w_orig_sub is not None else numpy.ones ( len ( X_orig_sub ), dtype = numpy.float32 )
        w_targ_abs = numpy.abs ( w_targ_sub ) if w_targ_sub is not None else numpy.ones ( len ( X_targ_sub ), dtype = numpy.float32 )
        
        sum_w_orig = numpy.sum ( w_orig_abs, dtype = numpy.float64 )
        sum_w_targ = numpy.sum ( w_targ_abs, dtype = numpy.float64 )
        
        # Scale factor aligning total target mass with total original mass
        scale_factor = ( sum_w_targ / ( sum_w_orig + 1e-12 ) ) if sum_w_orig > 0 else numpy.float32 ( 1.0 )
        w_orig_scaled = w_orig_abs * scale_factor
        
        w_comb_base     = numpy.hstack ( [ w_orig_scaled, w_targ_abs ] )
        total_len       = len ( X_orig_sub ) + len ( X_targ_sub )
        global_comb_sum = numpy.sum ( w_comb_base, dtype = numpy.float64 )
        w_comb_global   = w_comb_base * ( total_len / ( global_comb_sum + 1e-12 ) )
                        
        if 1 < self.n_splits :
            from sklearn.model_selection import StratifiedKFold
            skf     = StratifiedKFold ( n_splits     = self.n_splits     , 
                                        shuffle      = True              , 
                                        random_state = self.random_state )
            splits  = list ( skf.split ( X_comb, y_comb ) )
            n_folds = self.n_splits
        else:
            indices = numpy.arange ( len ( X_comb ) )
            splits  = [ ( indices , indices ) ]
            n_folds = 1

        stream_models   = []
        fold_preds_list = []
        n_samples       = len ( X_orig_sub )
        ratios          = numpy.ones ( n_samples, dtype = numpy.float32 )

        if self.__oof_mode :
            # --- MODE 1: Strict Out-Of-Fold (OOF) Evaluation ---
            # Each original sample is evaluated exclusively by the fold model where it was part of the validation set.
            for train_idx, valid_idx in progress_bar ( splits      ,
                                                       max_value   = n_folds           ,
                                                       description = 'Folds (OOF):'    , 
                                                       silent      = not self.progress ) :
                X_tr, y_tr = X_comb [ train_idx ] , y_comb [ train_idx ]
                w_tr       = w_comb_global [ train_idx ] if w_comb_global is not None else None
                
                model = self._train_single_model ( X_tr, y_tr, w_tr )
                stream_models.append ( model )

                orig_mask = valid_idx < n_samples
                if not numpy.any ( orig_mask ) : continue
                    
                orig_indices_global = valid_idx [ orig_mask ]
                orig_indices_local  = orig_indices_global

                X_val_fold = X_comb [ orig_indices_global ]
                raw_p      = self._predict_single_model ( model, X_val_fold )
                if 1 < raw_p.ndim : raw_p = raw_p [ :, 1 ]
                
                eps         = numpy.finfo ( numpy.float32 ).eps 
                p_clipped   = numpy.clip ( raw_p, eps, numpy.float32 ( 1.0 ) - eps )
                fold_ratios = p_clipped / ( numpy.float32 ( 1.0 ) - p_clipped )
                
                ratios [ orig_indices_local ] = fold_ratios.astype ( numpy.float32, copy = False )
        else :
            # --- MODE 2: Training-Ensemble Averaging (Cheating / Stabilized Mode) ---
            # Each sample's ratio is averaged across all N-1 models for which it served as a training point.
            for train_idx, valid_idx in progress_bar ( splits      ,
                                                       max_value   = n_folds           ,
                                                       description = 'Folds (Train Ens):', 
                                                       silent      = not self.progress ) :
                X_tr, y_tr = X_comb [ train_idx ] , y_comb [ train_idx ]
                w_tr       = w_comb_global [ train_idx ] if w_comb_global is not None else None
                
                model = self._train_single_model ( X_tr, y_tr, w_tr )
                stream_models.append ( model )

                # Predict on the entire original subset using the current fold model
                raw_p = self._predict_single_model ( model, X_orig_sub )
                if 1 < raw_p.ndim : raw_p = raw_p [ :, 1 ]
                
                eps         = numpy.finfo ( numpy.float32 ).eps 
                p_clipped   = numpy.clip ( raw_p, eps, numpy.float32 ( 1.0 ) - eps )
                fold_ratios = p_clipped / ( numpy.float32 ( 1.0 ) - p_clipped )
                fold_preds_list.append ( fold_ratios )

            if 1 < n_folds :
                all_fold_ratios = numpy.vstack ( fold_preds_list )
                
                # Build mask identifying where each model was validation (to exclude it)
                valid_mask_per_fold = numpy.zeros ( ( n_folds, n_samples ), dtype = bool )
                for fold_idx, (train_idx, valid_idx) in enumerate ( splits ) :
                    orig_mask = valid_idx < n_samples
                    if numpy.any ( orig_mask ) :
                        valid_mask_per_fold [ fold_idx, valid_idx [ orig_mask ] ] = True
                
                # Invert mask to keep only models where the sample was part of the training set
                train_mask_per_fold = ~valid_mask_per_fold
                sum_ratios   = numpy.zeros ( n_samples, dtype = numpy.float32 )
                count_models = numpy.zeros ( n_samples, dtype = numpy.int32 )
                
                for fold_idx in range ( n_folds ) :
                    mask = train_mask_per_fold [ fold_idx ]
                    sum_ratios [ mask ]   += all_fold_ratios [ fold_idx, mask ]
                    count_models [ mask ] += 1
                    
                count_models = numpy.maximum ( count_models, 1 )
                ratios       = sum_ratios / count_models
            else :
                ratios = fold_preds_list [ 0 ]

        # Save fitted state attributes for subsequent inference/evaluation
        self.fitted_models [ stream_key ] = stream_models 
        self.scale_factors [ stream_key ] = scale_factor
        self.priors        [ stream_key ] = numpy.float32 ( 1.0 )
        
        # Global post-processing & normalization of the ratios array
        mean_ratio = numpy.mean ( ratios, dtype = numpy.float64 )
        if 0 < mean_ratio > 0 : ratios = ratios / numpy.float32 ( mean_ratio )
        
        return ratios


    ## @brief Extract density ratios r(x) from a single fold's model.
    #  @details Acts as a polymorphic bridge for converting raw probability outputs 
    #           into ratio estimates r(x) = p(x) / (1 - p(x)).
    #  @param model Trained fold estimator.
    #  @param X Features matrix.
    #  @return Array of predicted density ratios.
    def _predict_fold_ratios ( self, model, X ) :
        """ Extract density ratios r(x) from a single fold's model. """
        p = self._predict_single_model ( model, X )
        if 1 < p.ndim : p = p [ :, 1 ]
        eps       = numpy.finfo ( numpy.float32 ).eps
        p_clipped = numpy.clip ( p, eps, numpy.float32 ( 1.0 ) - eps )
        return p_clipped / ( numpy.float32 ( 1.0 ) - p_clipped )

    ## @brief Hook for post-processing aggregated ensemble ratios.
    #  @details Can be overridden by subclasses (e.g., IterativeReweighter) to apply 
    #           custom normalization or scaling after fold prediction voting.
    #  @param ratios Aggregated ratios array.
    #  @return Post-processed ratios array.
    def _postprocess_ensemble_ratios ( self, ratios ) :
        """ Hook for post-processing aggregated ensemble ratios. """
        return ratios

    ## @brief Predict density ratios for an unseen feature matrix X using all trained N fold models.
    #  @details Evaluates each fold model independently and aggregates their predictions 
    #           matrix of shape (N_splits, N_samples) using the configured `Voting` strategy.
    #  @param stream_key Identifier key of the stream being evaluated.
    #  @param X Features matrix of unseen events.
    #  @return Aggregated density ratio predictions for the stream.
    def _predict_stream_ratios ( self, stream_key, X ) :
        """ Predict density ratios for an unseen feature matrix X using all trained N fold models. """
        models_list = self.fitted_models [ stream_key ]
        if not models_list :
            return numpy.ones ( len ( X ) , dtype = numpy.float32 )

        fold_preds = []
        for model in models_list :
            ratios_fold = self._predict_fold_ratios ( model, X )
            fold_preds.append ( ratios_fold )

        # Stack fold predictions into matrix of shape (N_splits, N_samples)
        stacked_ratios = numpy.vstack ( fold_preds )
        
        # Apply Voting strategy helper to aggregate predictions across the fold axis
        ensemble_ratios = self.voting ( stacked_ratios, axis = 0 )
        ensemble_ratios = ensemble_ratios.astype ( numpy.float32, copy = False )
        
        return self._postprocess_ensemble_ratios ( ensemble_ratios )

    ## @brief Clip ratios at threshold and compute global normalization scale factor.
    #  @param ratios Density ratios array.
    #  @param w_orig Original source weights array.
    #  @return Normalized and clipped density ratios array.
    def _normalize_and_clip ( self   ,
                              ratios ,
                              w_orig ) :
        """ Clip ratios at threshold and compute global normalization scale factor. """
        clipped_ratios    = numpy.clip ( ratios, numpy.float32 ( 0.0 ) , numpy.float32 ( self.__clip_threshold ) )
        orig_sum          = numpy.sum  ( w_orig, dtype = numpy.float64 )
        reweighted_sum    = numpy.sum  ( w_orig * clipped_ratios, dtype = numpy.float64 )
        self._norm_factor = numpy.float32 ( orig_sum / ( reweighted_sum + 1e-12 ) )
        return clipped_ratios * self._norm_factor

    ## @brief Evaluate reweighting factors for unseen source feature samples.
    #  @details Executes the stream prediction pipeline, combines streams according to sign mode, 
    #           applies clipping threshold, and scales the result by global normalization factor.
    #  @param original Source feature matrix of unseen samples.
    #  @param original_weight Source observation weights.
    #  @return Final reweighted factors for input samples.
    def weights ( self                   ,
                  original               ,
                  original_weight = None ) :
        """ Evaluate reweighting factors for unseen source feature samples. """
        
        X_new_f32 = original.astype ( numpy.float32, copy = False )

        if not valid_data_shape   ( X_new_f32                   ) : raise TypeError ( "Invalid `original` type/shape: %s" % typename ( original ) )
        if not valid_weight       ( original_weight             ) : raise TypeError ( "Invalid `original_weight`!" )        
        if not compatible_weights ( X_new_f32 , original_weight ) : raise TypeError ( "Incompatible `original` data/weight!" )
        if self.n_features != num_features ( X_new_f32 )          : raise TypeError ( "Invalid #features!!")

        if original_weight is None and self.__mode == "1-stream" :
            
            ratios  = self._predict_stream_ratios ( "base", X_new_f32 )
            clipped = numpy.clip ( ratios,
                                   numpy.float32 ( 0.0 ) ,
                                   numpy.float32 ( self.__clip_threshold ) )
            
            result  = clipped * self._norm_factor
            
            ## final weights 
            final_weight  = result 
            final_weight  = self.preserve_mass ( original , original_weight , final_weight )  
            
            ## summary table 
            if self.verbose : self.weight_summary ( original        ,
                                                    original_weight ,
                                                    final_weight    ) 
            return final_weight
                    
        w_new_f32 = ( original_weight.astype ( numpy.float32, copy = False )
                      if original_weight is not None
                      else numpy.ones ( len ( original ), dtype = numpy.float32 ) )

        nz_mask          = w_new_f32 != 0
        X_valid, w_valid = X_new_f32[ nz_mask ], w_new_f32[ nz_mask ]

        mask_pos     = w_valid > 0
        mask_neg     = w_valid < 0
        ratios_valid = numpy.zeros ( len ( X_valid ), dtype = numpy.float32 )

        if self.__mode == "1-stream" :
            ratios_valid = self._predict_stream_ratios ( "base", X_valid )

        elif self.__mode == "2-stream_orig" :
            if numpy.any ( mask_pos ) :
                ratios_valid[ mask_pos ] = self._predict_stream_ratios ( "orig_pos", X_valid[ mask_pos ] )
            if numpy.any ( mask_neg ) :
                ratios_valid[ mask_neg ] = self._predict_stream_ratios ( "orig_neg", X_valid[ mask_neg ] )

        elif self.__mode == "2-stream_targ" :
            w_pos   = self.__target_weights_info[ "W_pos" ]
            w_neg   = self.__target_weights_info[ "W_neg" ]
            w_total = w_pos + w_neg + numpy.float32 ( 1e-12 )

            r_pos_pos    = self._predict_stream_ratios ( "pos_pos", X_valid )
            r_pos_neg    = self._predict_stream_ratios ( "pos_neg", X_valid )
            ratios_valid = ( r_pos_pos * w_pos - r_pos_neg * w_neg ) / w_total

        elif self.__mode == "4-stream" :
            w_pos   = self.__target_weights_info[ "W_pos" ]
            w_neg   = self.__target_weights_info[ "W_neg" ]
            w_total = w_pos + w_neg + numpy.float32 ( 1e-12 )

            if numpy.any ( mask_pos ) :
                X_p                      = X_valid[ mask_pos ]
                r_pos_pos                = self._predict_stream_ratios ( "pos_pos", X_p )
                r_pos_neg                = self._predict_stream_ratios ( "pos_neg", X_p )
                ratios_valid[ mask_pos ] = ( r_pos_pos * w_pos - r_pos_neg * w_neg ) / w_total

            if numpy.any ( mask_neg ) :
                X_n                      = X_valid[ mask_neg ]
                r_neg_pos                = self._predict_stream_ratios ( "neg_pos", X_n )
                r_neg_neg                = self._predict_stream_ratios ( "neg_neg", X_n )
                ratios_valid[ mask_neg ] = ( r_neg_pos * w_pos - r_neg_neg * w_neg ) / w_total

        clipped                 = numpy.clip ( ratios_valid,
                                               numpy.float32 ( 0.0 ),
                                               numpy.float32 ( self.__clip_threshold ) )
        final_ratios            = numpy.zeros ( len ( w_new_f32 ), dtype = numpy.float32 )
        final_ratios[ nz_mask ] = self._apply_norm_factor ( clipped ) 

        ## summary table 
        if self.verbose : self.weight_summary ( original       ,
                                                orginal_weight ,
                                                final_ratios   ) 
        
        return final_ratios if original_weight is None else final_ratios * w_new_f32

    ## @brief Apply global normalization factor to clipped ratios.
    #  @param clipped_ratios Density ratios array post-clipping.
    #  @return Normalized ratios array.
    def _apply_norm_factor ( self, clipped_ratios ) :
        """ Apply global normalization factor to clipped ratios. """
        return clipped_ratios * self._norm_factor

# =============================================================================
## @class IterativeReweighter
#  @brief Base class for iterative weight boosting density reweighters.
#  @details Adapts sequential iterative density ratio estimation algorithms 
#           (LightGBM, XGBoost, CatBoost, and PyTorch iterative wrappers) 
#           to fit seamlessly into the cross-validation and evaluation pipeline 
#           provided by the DensityReweighter base class.
class IterativeReweighter ( DensityReweighter ) :
    """ Base class for iterative weight boosting density reweighters (Gradient Weight Boosting). """
    
    ## @brief Initialize the iterative weight boosting ensemble.
    #  @param original Source feature matrix.
    #  @param target Target feature matrix.
    #  @param original_weight Source observation weights.
    #  @param target_weight Target observation weights.
    #  @param n_estimators Number of boosting iterations per fold.
    #  @param learning_rate Boosting learning rate.
    #  @param max_depth Maximum tree depth for tree-based backends.
    #  @param clip_threshold Ratio clipping threshold.
    #  @param store_original_weights Retain original reweighted arrays.
    #  @param progress Enable progress bar display.
    #  @param params Additional backend boosting parameters.
    def __init__ ( self                   , * ,
                   original               ,
                   target                 ,
                   original_weight        = None  ,
                   target_weight          = None  ,
                   n_estimators           = 150   ,
                   learning_rate          = 0.05  ,
                   max_depth              = 4     ,
                   clip_threshold         = 1.e+4 ,
                   store_original_weights = True  ,
                   progress               = True  , **params ) :
        """ Initialize the iterative weight boosting ensemble. """
        
        # Inject core boosting hyperparameters into parameter dictionary
        params [ 'n_estimators'  ] = n_estimators
        params [ 'learning_rate' ] = learning_rate
        params [ 'max_depth'     ] = max_depth
        
        # Delegate initialization and immediate fitting to parent DensityReweighter
        super ().__init__ ( original               = original               ,
                            target                 = target                 ,
                            original_weight        = original_weight        ,
                            target_weight          = target_weight          ,
                            clip_threshold         = clip_threshold         ,
                            store_original_weights = store_original_weights ,
                            progress               = progress               , **params )

    @property
    def method ( self ) :
        """ Return method identifier name. """
        return '%s IDRW' % self.symbol 

    ## @property norm_factor
    #  @brief Normalization scale handled internally per boosting step.
    @property
    def norm_factor ( self ) :
        """ `norm_factor`: Normalization scale handled internally per boosting step. """
        return numpy.float32 ( 1.0 )
    
    ## @brief Abstract method for concrete backends to execute sequential boosting steps.
    #  @param X_orig Original feature subset for current fold.
    #  @param w_orig Original weight subset for current fold.
    #  @param X_targ Target feature subset for current fold.
    #  @param w_targ Target weight subset for current fold.
    #  @return Tuple of (models_list, final_ratios_array).
    @abc.abstractmethod
    def _fit_iterative_stream ( self, X_orig, w_orig, X_targ, w_targ ) :
        """ Abstract method for concrete backends to execute sequential boosting steps. """
        raise NotImplementedError

    ## @brief Abstract method to accumulate log-ratios across sequential models for matrix X.
    #  @param models List of trained sequential boosting models for a fold.
    #  @param X Unseen feature matrix.
    #  @return Density ratios array predicted by boosting chain.
    @abc.abstractmethod
    def _predict_iterative_stream ( self, models, X ) :
        """ Abstract method to accumulate log-ratios across sequential models for matrix X. """
        raise NotImplementedError

    ## @brief Predict ratios directly using the iterative sequential chain.
    #  @param models Trained model chain.
    #  @param X Features matrix.
    #  @return Array of density ratios.
    def _predict_fold_ratios ( self, models, X ) :
        """ Predict ratios directly using the iterative sequential chain. """
        return self._predict_iterative_stream ( models, X )

    ## @brief Normalize aggregated ensemble response mean value to 1.0.
    #  @param ratios Aggregated ratios array.
    #  @return Normalized ratios array.
    def _postprocess_ensemble_ratios ( self, ratios ) :
        """ Normalize aggregated ensemble response mean value to 1.0. """
        mean_ens = numpy.mean ( ratios )
        if 0 < mean_ens : ratios /= mean_ens
        return ratios

    ## @brief Train an iterative boosting chain on the training partition of the fold.
    #  @param X_train Combined training feature matrix (N-1)/N.
    #  @param y_train Binary labels for training fold.
    #  @param w_train Training sample weights.
    #  @return List of trained sequential models for the fold.
    def _train_single_model ( self, X_train, y_train, w_train ) :
        """ Train an iterative boosting chain on the training partition of the fold. """
        orig_mask = ( y_train == 0 )
        targ_mask = ( y_train == 1 )

        X_orig_tr = X_train [ orig_mask ]
        w_orig_tr = w_train [ orig_mask ] if w_train is not None else None
        X_targ_tr = X_train [ targ_mask ]
        w_targ_tr = w_train [ targ_mask ] if w_train is not None else None

        models, _ = self._fit_iterative_stream ( X_orig_tr , w_orig_tr ,
                                                 X_targ_tr , w_targ_tr )
        return models

    ## @brief Map predicted ratios r(x) to pseudo-probabilities p = r / (1 + r).
    #  @param models List of sequential boosting models.
    #  @param X Features matrix.
    #  @return Array of pseudo-probabilities.
    def _predict_single_model ( self, models, X ) :
        """ Map predicted ratios r(x) to pseudo-probabilities p = r / (1 + r). """
        ratios = self._predict_iterative_stream ( models, X )
        p      = ratios / ( 1.0 + ratios )
        return p.astype ( numpy.float32, copy = False )

    ## @brief Soft normalization for iterative fold ensembles keeping mean weight stable.
    #  @param clipped_ratios Clipped ratios array.
    #  @return Normalized ratios array.
    def _apply_norm_factor ( self, clipped_ratios ) :
        """ Soft normalization for iterative fold ensembles keeping mean weight stable. """
        mean_val = numpy.mean ( clipped_ratios )
        if 0 < mean_val :  return clipped_ratios / mean_val
        return clipped_ratios

    ## @brief Clip extreme ratios and compute mass conservation norm factor for logging.
    #  @param ratios Unclipped ratios array.
    #  @param w_orig Original source weights.
    #  @return Clipped and locally normalized ratios array.
    def _normalize_and_clip ( self, ratios, w_orig = None ) :
        """ Clip extreme ratios and compute mass conservation norm factor for logging. """
        clipped = numpy.clip ( ratios, 0.0, self.params.get ( 'clip_threshold', 1.e4 ) )
        
        if w_orig is None :
            mass_before = len ( ratios )
            mass_after  = numpy.sum ( clipped )
        else :
            mass_before = numpy.sum ( numpy.abs ( w_orig ) )
            mass_after  = numpy.sum ( numpy.abs ( w_orig * clipped ) )
            
        norm_factor = mass_before / ( mass_after + 1e-12 )
        
        logger.debug (
            "[IterativeReweighter] mass_before=%.2f, mass_after=%.2f, norm_factor=%.4f" % 
            ( mass_before, mass_after, norm_factor )
        )
        
        return clipped * norm_factor

# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
#                                                                       The END 
# =============================================================================
