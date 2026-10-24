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
from   ostap.core.ostap_types   import num_types 
from   ostap.utils.core         import typename
from   ostap.utils.basic        import numcpu, num_jobs, NoContext
from   ostap.logger.utils       import map2table_ex
from   ostap.logger.pretty      import nice_print 
from   ostap.utils.progress_bar import progress_bar, text_bar
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
from   ostap.utils.memory       import memory
from   ostap.utils.timing       import timing 
import ostap.logger.symbols   as     S
import numpy, abc, math, gc  
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.base_reweighters' )
else                       : logger = getLogger( __name__ )
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
    def __init__ ( self                           , * ,
                   original                       ,
                   target                         ,
                   original_weight        = None  , 
                   target_weight          = None  ,
                   clip_threshold         = 1.e+4 ,
                   n_splits               = 5     ,
                   store_original_weights = False ,
                   progress               = True  , **params ) :
        """ Initialize and fit the density ratio reweighter ensemble.
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
        
        self.__progress            = True if progress else False 
        self.__clip_threshold      = float ( clip_threshold )
        self.__n_splits            = n_splits

        # Public attributes for clean and direct access in subclasses
        self.fitted_models         = {}
        self.priors                = {}
        self.__target_weights_info = {}
        self.__norm_factor         = numpy.float32( 1.0 )
        self.__mode                = None
        self.scale_factors         = {}
        
        params [ 'n_jobs' ] = num_jobs ( params )
                
        self.__original_ratios             = None
        self.__original_reweighted_weights = None

        super ().__init__ ( original        = original        ,
                            target          = target          , 
                            original_weight = original_weight ,
                            target_weight   = target_weight   , **params )
        
        original_ratios, original_reweighted_weights = (
            self.__fit_and_compute (
                original.astype ( numpy.float32 , copy = False ),
                (
                    original_weight.astype ( numpy.float32, copy = False )
                    if original_weight is not None
                    else None
                ),
                target.astype   ( numpy.float32 , copy = False ),
                (
                    target_weight.astype  ( numpy.float32, copy = False )
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

    @property
    def n_splits ( self ) :
        return self.__n_splits

    @property
    def original_ratios( self ) :
        return self.__original_ratios

    @property
    def original_reweighted_weights( self ) :
        return self.__original_reweighted_weights
    
    @property
    def norm_factor ( self ) :
        return self._norm_factor
     
     
    ## @property 
    ## def scale_factors ( self ) :
    ##    return self.scale_factors 
        
    @property
    def mode ( self ) :
        return self.__mode

    @property
    def progress ( self ) :
        return self.__progress
    
    @property
    def config ( self ) :
        conf = {}
        conf.update( super().config )
        conf [ "progress" ]                    = self.progress
        conf [ "mode"     ]                    = self.mode
        conf [ "clip_threshold" ]              = self.__clip_threshold
        conf [ "n_splits" ]                    = self.__n_splits
        conf [ "original_ratios" ]             = self.__original_ratios             is not None
        conf [ "original_reweighted_weights" ] = self.__original_reweighted_weights is not None
        return conf

    @abc.abstractmethod
    def _train_single_model ( self    ,
                              X_train ,
                              y_train ,
                              w_train ,
                              X_val   ,
                              y_val   ,
                              w_val   ) :
        raise NotImplementedError

    @abc.abstractmethod
    def _predict_single_model( self  ,
                               model ,
                               X     ) :
        raise NotImplementedError

    def __fit_and_compute ( self            ,
                            original        ,
                            original_weight ,
                            target          ,
                            target_weight   ) :
        if original_weight is None and target_weight is None :
            self.__mode  = "1-stream"
            ratios_valid = self._fit_eval_stream ( "base", original, None, target, None )

            clipped            = numpy.clip    ( ratios_valid,
                                                 numpy.float32 ( 0.0 ) ,
                                                 numpy.float32 ( self.__clip_threshold ) )
            self.__norm_factor = numpy.float32 ( len( original )
                                                 / ( numpy.sum ( clipped, dtype = numpy.float64 ) + 1e-12 ) )

            original_ratios = self._apply_norm_factor ( clipped ) 
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
            ratios_valid = self._fit_eval_stream ( "base"       ,
                                                   X_orig_valid ,
                                                   w_orig_valid ,
                                                   X_targ_valid ,
                                                   w_targ_valid )

        elif has_neg_orig and not has_neg_targ :
            self.__mode                   = "2-stream_orig"
            ratios_valid                  = numpy.zeros ( len( X_orig_valid ) , dtype = numpy.float32 )
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

            r_pos_pos = self._fit_eval_stream ( "pos_pos", X_orig_pos, w_orig_pos, X_targ_pos, w_targ_pos )
            r_pos_neg = self._fit_eval_stream ( "pos_neg", X_orig_pos, w_orig_pos, X_targ_neg, w_targ_neg )
            r_neg_pos = self._fit_eval_stream ( "neg_pos", X_orig_neg, w_orig_neg, X_targ_pos, w_targ_pos )
            r_neg_neg = self._fit_eval_stream ( "neg_neg", X_orig_neg, w_orig_neg, X_targ_neg, w_targ_neg )

            ratios_valid = numpy.zeros ( len ( X_orig_valid ) , dtype = numpy.float32 )
            denom                         = w_total + numpy.float32( 1e-12 )
            ratios_valid[ mask_orig_pos ] = ( r_pos_pos * w_pos_sum - r_pos_neg * w_neg_sum ) / denom
            ratios_valid[ mask_orig_neg ] = ( r_neg_pos * w_pos_sum - r_neg_neg * w_neg_sum ) / denom

        ratios_valid_norm = self._normalize_and_clip ( ratios_valid, w_orig_valid )

        original_ratios                 = numpy.zeros( len( w_orig ), dtype = numpy.float32 )
        original_ratios[ nz_orig_mask ] = ratios_valid_norm
        original_reweighted_weights     = w_orig * original_ratios

        return original_ratios, original_reweighted_weights

    def _fit_eval_stream ( self       ,
                           stream_key ,
                           X_orig_sub ,
                           w_orig_sub ,
                           X_targ_sub ,
                           w_targ_sub ) :
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
                                                 max_value   = n_folds           ,
                                                 description = 'Folds:'          , 
                                                 silent      = not self.progress ) :
            
            X_tr, y_tr = X_comb [ train_idx ] , y_comb [ train_idx ]
            X_va, y_va = X_comb [ val_idx   ] , y_comb [ val_idx ]

            w_tr = w_comb [ train_idx ] if w_comb is not None else None
            w_va = w_comb [ val_idx   ] if w_comb is not None else None
            
            model, _ = self._train_single_model ( X_tr, y_tr, w_tr, X_va, y_va, w_va )
            stream_models.append( model )

            raw_val_p = self._predict_single_model ( model, X_va )
            if 1 < raw_val_p.ndim : raw_val_p = raw_val_p[:, 1]
                
            oof_raw [ val_idx ] = raw_val_p.astype ( numpy.float32, copy = False )

        self.fitted_models [ stream_key ] = stream_models 
        self.scale_factors [ stream_key ] = scale_factor
        self.priors        [ stream_key ] = numpy.float32 ( 1.0 )

        p_orig = oof_raw[ : len ( X_orig_sub ) ]
        eps    = 1.e-4 
        p_orig_clipped = numpy.clip ( p_orig, numpy.float32 ( eps ) , numpy.float32( 1.0 - eps  ) )
        
        ratios = p_orig_clipped / ( numpy.float32( 1.0 ) - p_orig_clipped )
        return ratios

    def _predict_stream_ratios( self, stream_key, X ):
        models = self.fitted_models[ stream_key ]
        ratios_sum = numpy.zeros(len(X), dtype=numpy.float32)
        
        eps = 1.e-4 
        for model in models :
            p = self._predict_single_model ( model, X )
            if p.ndim > 1 : p = p[:, 1]
            
            p_clipped   = numpy.clip ( p, numpy.float32 ( eps ), numpy.float32 ( 1.0 - eps ) )
            ratios_sum += p_clipped / ( numpy.float32 ( 1.0 ) - p_clipped )

        return ratios_sum / len(models)

    def _normalize_and_clip ( self   ,
                              ratios ,
                              w_orig ) :
        clipped_ratios = numpy.clip ( ratios, numpy.float32 ( 0.0 ) , numpy.float32 ( self.__clip_threshold ) )
        orig_sum       = numpy.sum  ( w_orig, dtype = numpy.float64 )
        reweighted_sum = numpy.sum  ( w_orig * clipped_ratios, dtype = numpy.float64 )
        self.__norm_factor = numpy.float32 ( orig_sum / ( reweighted_sum + 1e-12 ) )
        return clipped_ratios * self.__norm_factor

    def weights ( self                   ,
                  original               ,
                  original_weight = None ) :
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
            ratios_valid = self._predict_stream_ratios( "base", X_valid )

        elif self.__mode == "2-stream_orig" :
            if numpy.any( mask_pos ) :
                ratios_valid[ mask_pos ] = self._predict_stream_ratios ( "orig_pos", X_valid[ mask_pos ] )
            if numpy.any( mask_neg ) :
                ratios_valid[ mask_neg ] = self._predict_stream_ratios ( "orig_neg", X_valid[ mask_neg ] )

        elif self.__mode == "2-stream_targ" :
            w_pos   = self.__target_weights_info[ "W_pos" ]
            w_neg   = self.__target_weights_info[ "W_neg" ]
            w_total = w_pos + w_neg + numpy.float32( 1e-12 )

            r_pos_pos    = self._predict_stream_ratios( "pos_pos", X_valid )
            r_pos_neg    = self._predict_stream_ratios( "pos_neg", X_valid )
            ratios_valid = ( r_pos_pos * w_pos - r_pos_neg * w_neg ) / w_total

        elif self.__mode == "4-stream" :
            w_pos   = self.__target_weights_info[ "W_pos" ]
            w_neg   = self.__target_weights_info[ "W_neg" ]
            w_total = w_pos + w_neg + numpy.float32( 1e-12 )

            if numpy.any( mask_pos ) :
                X_p                      = X_valid[ mask_pos ]
                r_pos_pos                = self._predict_stream_ratios ( "pos_pos", X_p )
                r_pos_neg                = self._predict_stream_ratios ( "pos_neg", X_p )
                ratios_valid[ mask_pos ] = ( r_pos_pos * w_pos - r_pos_neg * w_neg ) / w_total

            if numpy.any( mask_neg ) :
                X_n                      = X_valid[ mask_neg ]
                r_neg_pos                = self._predict_stream_ratios ( "neg_pos", X_n )
                r_neg_neg                = self._predict_stream_ratios ( "neg_neg", X_n )
                ratios_valid[ mask_neg ] = ( r_neg_pos * w_pos - r_neg_neg * w_neg ) / w_total

        clipped                 = numpy.clip ( ratios_valid,
                                               numpy.float32 ( 0.0 ),
                                               numpy.float32 ( self.__clip_threshold ) )
        final_ratios            = numpy.zeros( len( w_new_f32 ), dtype = numpy.float32 )
        final_ratios[ nz_mask ] = self._apply_norm_factor ( clipped ) 

        return final_ratios if original_weight is None else final_ratios * w_new_f32

    def _apply_norm_factor ( self, clipped_ratios ) :
        """ Apply global normalization factor to clipped ratios """
        return clipped_ratios * self.__norm_factor
    
# =============================================================================
## @class IterativeReweighter
#  Base class for iterative weight boosting density reweighters.
#  Implements Gradient Weight Boosting (GWB) combined with cross-validation 
#  via the Template Method pattern inherited from DensityReweighter.
class IterativeReweighter ( DensityReweighter ) :
    """Base class for iterative weight boosting density reweighters (Gradient Weight Boosting).
    
    This class adapts sequential iterative density ratio estimation algorithms 
    (such as LightGBM, XGBoost, CatBoost, and PyTorch iterative wrappers) 
    to fit seamlessly into the cross-validation and evaluation pipeline 
    provided by the DensityReweighter base class[cite: 2].
    """
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
        
        # Inject core boosting hyperparameters into the parameter dictionary
        params [ 'n_estimators'  ] = n_estimators
        params [ 'learning_rate' ] = learning_rate
        params [ 'max_depth'     ] = max_depth
        
        # Delegate initialization and immediate fitting to the parent DensityReweighter[cite: 2]
        super().__init__ ( original               = original               ,
                           target                 = target                 ,
                           original_weight        = original_weight        ,
                           target_weight          = target_weight          ,
                           clip_threshold         = clip_threshold         ,
                           store_original_weights = store_original_weights ,
                           progress               = progress               , **params )
    @property
    def norm_factor ( self ) :
        """ For iterative reweighters, normalization is handled internally per-step. """
        return numpy.float32 ( 1.0 )
    

    # =========================================================================
    ## @brief Abstract method to execute sequential boosting steps.
    #  @details To be implemented by concrete backend reweighters (e.g., LightGBM, XGBoost).
    #  @param self Object instance.
    #  @param X_orig Original features subset.
    #  @param w_orig Original weights subset.
    #  @param X_targ Target features subset.
    #  @param w_targ Target weights subset.
    #  @return Tuple of (models list, final ratios array).
    @abc.abstractmethod
    def _fit_iterative_stream      ( self      ,
                                     X_orig    ,
                                     w_orig    ,
                                     X_targ    ,
                                     w_targ    ) :
        """Abstract method to be implemented by concrete backend reweighters 
        to execute sequential boosting steps.
        """
        raise NotImplementedError

    # =========================================================================
    ## @brief Abstract method to accumulate log-ratios or margins.
    #  @details Evaluates predictions across all trained sequential models for an unseen feature matrix X.
    #  @param self Object instance.
    #  @param models List of trained sequential models.
    #  @param X Unseen feature matrix.
    #  @return Array of density ratios predicted by the iterative chain.
    @abc.abstractmethod
    def _predict_iterative_stream ( self       ,
                                    models     ,
                                    X          ) :
        """Abstract method to accumulate log-ratios or margins across all trained 
        sequential models for an unseen feature matrix X.
        """
        raise NotImplementedError


    # =========================================================================
    ## @brief Template Method Hook: Predicting stream ratios across CV folds.
    #  @details Computes the geometric mean of predicted density ratios across all 
    #           cross-validation folds by accumulating log-ratios before exponentiation.
    #           This matches the additive log-margin mechanics of gradient boosting,
    #           maintains total mass consistency, and prevents upper-bound divergence.
    #  @param self Object instance.
    #  @param stream_key Stream identifier string.
    #  @param X Unseen feature matrix to predict density ratios for.
    #  @return Averaged numpy array of density ratios.
    def _predict_stream_ratios ( self, stream_key, X ) :
        """ Predict density ratios across all CV folds and aggregate them using the median 
        to guarantee that the ensemble maximum never exceeds individual fold limits.
        """
        models_list = self.fitted_models [ stream_key ]
        if not models_list :
            return numpy.ones ( len ( X ) , dtype = numpy.float32 )

        # Collect predictions from every single fold model chain
        fold_preds = []
        for models in models_list :
            ratios_fold = self._predict_iterative_stream ( models, X )
            fold_preds.append ( ratios_fold )

        # Stack into shape (n_folds, len(X))
        stacked = numpy.vstack ( fold_preds )

        # Take the median across folds. 
        # Median is mathematically bounded by the min and max of the individual folds,
        # preventing any out-of-fold tail accumulation or explosion.
        ensemble_ratios = numpy.median ( stacked, axis = 0 )

        return ensemble_ratios.astype ( numpy.float32, copy = False )

    # =========================================================================
    ## @brief Template Method Hook: Training a single fold.
    #  @details Trains an iterative boosting chain on the training partition of the current fold.
    #           Acts as a bridge, unpacking binary-labeled combined fold data (y=0 for original, y=1 for target)
    #           and passing them separately into the backend's iterative training routine.
    #  @param self Object instance.
    #  @param X_train Combined training feature matrix for the fold.
    #  @param y_train Combined binary labels (0 for original, 1 for target).
    #  @param w_train Observation weights for the training partition.
    #  @param X_val Combined validation feature matrix.
    #  @param y_val Combined validation labels.
    #  @param w_val Validation observation weights.
    #  @return Tuple of (list of trained sequential models, dummy validation predictions array).
    def _train_single_model ( self, X_train, y_train, w_train, X_val, y_val, w_val ) :
        """Train an iterative boosting chain on the training partition of the current fold.
        
        This method acts as a bridge, unpacking the binary-labeled combined fold data 
        (where y=0 represents original samples and y=1 represents target samples) 
        and passing them separately into the backend's iterative training routine.
        """
        # Separate original (y=0) and target (y=1) samples within the current fold
        orig_mask = ( y_train == 0 )
        targ_mask = ( y_train == 1 )

        X_orig_tr = X_train [ orig_mask ]
        w_orig_tr = w_train [ orig_mask ] if w_train is not None else None
        X_targ_tr = X_train [ targ_mask ]
        w_targ_tr = w_train [ targ_mask ] if w_train is not None else None

        # Execute sequential gradient weight boosting on the fold's training data
        models, _ = self._fit_iterative_stream ( X_orig_tr , w_orig_tr ,
                                                 X_targ_tr , w_targ_tr )
        
        # Generate dummy placeholder predictions to satisfy DensityReweighter's unpacking requirement 
        # (`model, val_preds = self._train_single_model(...)`). 
        # Actual fold ensembling is handled natively via _predict_stream_ratios.
        dummy_val_preds = numpy.zeros ( len ( X_val ) , dtype = numpy.float32 )
        
        return models, dummy_val_preds

    # =========================================================================
    ## @brief Template Method Hook: Predicting fold validation probabilities.
    #  @details Iterative reweighters naturally compute density ratios r(x). 
    #           To fit into the standard classification-based pipeline of DensityReweighter, 
    #           these ratios are mapped to equivalent pseudo-probabilities via p = r / (1 + r).
    #  @param self Object instance.
    #  @param models List of trained sequential weak learners for the fold.
    #  @param X Feature matrix to predict on.
    #  @return Numpy array of float32 pseudo-probabilities.
    def _predict_single_model ( self, models, X ) :
        """Predict pseudo-probabilities for validation/evaluation data.
        
        Iterative reweighters naturally compute density ratios r(x). 
        To fit into the standard classification-based pipeline of DensityReweighter, 
        these ratios are mapped to equivalent pseudo-probabilities via p = r / (1 + r).
        """
        # Obtain direct density ratio estimations from the iterative model chain
        ratios = self._predict_iterative_stream ( models, X )
        
        # Clip ratios to a physically reasonable range to prevent overflow in pseudo-probability conversion
        ratios_clipped = numpy.clip ( ratios, 1.e-5, self.config [ 'clip_threshold' ] )
        
        # Convert density ratios r(x) back to equivalent classification probabilities p(y=1|x)
        # Mathematical identity: if r = p / (1 - p), then p = r / (1 + r)
        p = ratios_clipped / ( 1.0 + ratios_clipped )
        
        return p.astype ( numpy.float32, copy = False )


    def _normalize_and_iter_clip ( self, ratios, w_orig, w_targ ) :
        """Clip extreme density ratio values and normalize total statistical mass 
        between original and target distributions.
        """
        clipped  = numpy.clip ( ratios , 0.0 , self.params.get ( 'clip_threshold', 1e4 ) ) 
        orig_sum = numpy.sum  ( numpy.abs ( w_orig ), dtype = numpy.float64 )
        targ_sum = numpy.sum  ( numpy.abs ( w_targ ), dtype = numpy.float64 )
        norm = targ_sum / ( orig_sum + 1e-12 )
        return clipped * norm

    def _normalize_and_clip ( self, ratios, w_orig ) :
        """Clip extreme density ratio values for iterative streams."""
        return numpy.clip ( ratios, numpy.float32 ( 0.0 ) , numpy.float32 ( self.params.get ( 'clip_threshold', 1e4 ) ) )
       
    def _apply_norm_factor ( self, clipped_ratios ) :
        """Мягкая нормализация для итеративного ансамбля фолдов, 
           сохраняющая среднее значение весов стабильным.
        """
        ## mean_val = numpy.mean ( clipped_ratios )
        ## if 0 < mean_val : return clipped_ratios / mean_val
        return clipped_ratios

# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
#                                                                       The END 
# =============================================================================
