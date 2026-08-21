#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/adversarial_validation.py
#  Implement adversarial validation to probe
#  the difference between two weighted datasets
#
#  @see ADVAL_LGBM , class based on LightGBM, the most CPU efficient 
#  @see ADVAL_HGBC , class based on HistGradientBoostingClassifier
#  @see ADVAL_GBC  , class based on GradientBoostingClassifier
#  @see ADVAL_RF   , class based on RandomForestClassifier 
#  @see ADVAL_XGB  , class based on XGBoost 
#  @see ADVAL_CATB , class based on CatBoost 
#  @see ADVAL_TORCH, class based on PyTorch
#  @see ADVAL_KERAS, class based on Keras
#
#  As t-value \f$ 100 \times \left( 1 - 2 \times AUC \right)^2\f$ is used.
#  To estimate the p-value, permutations are used.
#
#  This adversarial validation pipeline and LightGBM multiprocessing
#  architecture were co-developed and optimized with Gemini (Google AI). 
# 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-04
# =============================================================================
""" Implement adversarial validation to probe
    the difference between two weighted datasets

 - see ADVAL_LGBM  , class based on LightGBM, the most CPU efficient 
 - see ADVAL_HGBC  , class based on HistGradientBoostingClassifier
 - see ADVAL_GBC   , class based on GradientBoostingClassifier
 - see ADVAL_RF    , class based on RandomForestClassifier 
 - see ADVAL_XGB   , class based on XGBoost 
 - see ADVAL_CATB  , class based on CatBoost 
 - see ADVAL_TORCH , class based on PyTorch
 - see ADVAL_KERAS , class based on Keras

As t-value: 100 * ( 1 - 2 * AUC )**2 is used.
To estimate the p-value, permutations are used. 

  This adversarial validation pipeline and LightGBM multiprocessing
  architecture were co-developed and optimized with Gemini (Google AI).
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2026-07-04"
__all__     = (
    'ADVAL_LGBM'  , ## LightGBM-based class for Adversarial Validation
    'ADVAL_XGB'   , ## XGBoost-based class for Adversarial Validation
    'ADVAL_CATB'  , ## CatBoost-based class for Adversarial Validation
    'ADVAL_HGBC'  , ## HGBC-based class for Adversarial Validation
    'ADVAL_GBC'   , ## GBC-based class for Adversarial Validation
    'ADVAL_RF'    , ## RandomForest-based class for Adversarial Validation
    'ADVAL_TORCH' , ## PyTorch-based class for Adversarial Validation
    'ADVAL_KERAS' , ## Keras-based class for Adversarial Validation    
)
# =============================================================================
from   ostap.core.ostap_types import string_types
from   ostap.utils.core       import typename
from   ostap.utils.basic      import numcpu, num_jobs, run_parallel
from   ostap.stats.utils      import ( weight_trivial , compatible_shapes ,
                                       num_samples    , num_features      ) 
from   ostap.stats.gof_np     import GoFnp 
from   sklearn.metrics        import roc_auc_score
import ROOT, numpy, abc, os   
# =============================================================================
# Logging setup 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.adval' )
else                       : logger = getLogger ( __name__ )
# =============================================================================
logger.debug ( 'Implement Adversarial Validation for Goodness-of-fit & Two-Samples test' )
# =============================================================================
DEFAULT_ESTIMATORS         = 500
MAX_REGULARIZED_ESTIMATORS = 100

# =============================================================================
## Calculate t-value from AUC metric
#  t-value is defined as \f$ 100 \times \left( 1 - 2 \times AUC \right)^2 \f$ 
def tvalue_from_AUC ( auc ) :
    """ Calculate t-value defined as 100 * (1 - 2 * AUC)**2
    """
    return 100.0 * ( 1.0 - 2.0 * auc ) ** 2


# =============================================================================
## Invert binary labels or probabilities (1 - x) where sample_weight < 0.
def invert_if_negative_weight ( values , sample_weight ) :
    """ Invert binary labels or probabilities (1 - x) where sample_weight < 0."""
    if weight_trivial ( sample_weight ) : return values
    
    weights = sample_weight.values if hasattr(sample_weight, 'values') else sample_weight
    
    if ( weights < 0 ).any():
        return numpy.where(weights < 0, 1.0 - values, values)
    
    return values

# =============================================================================
## Transform targets and sample weights to resolve negative weight constraints.
# 
#   If sample_weight contains negative values (w < 0), flips the binary target 
#   (y = 1 - y) and takes the absolute value of the weight (|w|).
#  
#   Returns:
#   y_adj, w_adj : tuple of (np.ndarray or pd.Series, np.ndarray or pd.Series or None)
def transform_weighted_target ( y , sample_weight):
    """ Transform targets and sample weights to resolve negative weight constraints.
    
    If sample_weight contains negative values (w < 0), flips the binary target 
    (y = 1 - y) and takes the absolute value of the weight (|w|).
    
    Returns:
    --------
    y_adj, w_adj : tuple of (np.ndarray or pd.Series, np.ndarray or pd.Series or None)
    """
    if sample_weight is None : return y , None
    
    # Handle pandas Series or numpy arrays consistently
    weights = sample_weight.values if hasattr ( sample_weight , 'values') else sample_weight
    
    if ( weights < 0 ).any():
        y_adj = invert_if_negative_weight ( y , weights )
        w_adj = numpy.abs ( weights )
        return y_adj, w_adj

    return y , sample_weight

# =============================================================================
## Compute ROC-AUC score safely handling negative, missing, or trivial sample weights.
#
#  Standard scikit-learn `roc_auc_score` raises a ValueError if sample_weight contains
#  negative values. This wrapper mathematically resolves negative weights by transforming 
#  instances with w < 0 into their equivalent dual form: y_new = 1 - y and w_new = |w|.
#
#  It also optimizes execution by bypassing sample weighting entirely when weights 
#  are trivial (None, empty, or constant).
#
#  @param y_true : array-like of shape (n_samples,)
#                  True binary target labels (0 or 1).
#  @param y_pred : array-like of shape (n_samples,)
#                  Target scores or predicted probabilities P(Class=1).
#  @param sample_weight : array-like of shape (n_samples,) or None, optional
#                         Sample weights. Can contain negative or zero values.
#  @param kwargs : dict
#                 Additional keyword arguments passed directly to sklearn.metrics.roc_auc_score
#               (e.g., max_fpr, multi_class).
#  @return score : float Calculated ROC-AUC score.
def safe_roc_auc_score ( y_true , y_pred , sample_weight = None , **kwargs):
    """ Compute ROC-AUC score safely handling negative, missing, or trivial sample weights.
    
    Standard scikit-learn `roc_auc_score` raises a ValueError if sample_weight contains
    negative values. This wrapper mathematically resolves negative weights by transforming 
    instances with w < 0 into their equivalent dual form: y_new = 1 - y and w_new = |w|.
    
    It also optimizes execution by bypassing sample weighting entirely when weights 
    are trivial (None, empty, or constant).

    Parameters:
    -----------
    y_true : array-like of shape (n_samples,)
        True binary target labels (0 or 1).
    y_pred : array-like of shape (n_samples,)
        Target scores or predicted probabilities P(Class=1).
    sample_weight : array-like of shape (n_samples,) or None, optional
        Sample weights. Can contain negative or zero values.
    **kwargs : dict
        Additional keyword arguments passed directly to sklearn.metrics.roc_auc_score
        (e.g., max_fpr, multi_class).

    Returns:
    --------
    score : float
        Calculated ROC-AUC score.
    """
    # 1. Fast path: If weights are missing or constant (e.g. all 1s),
    # compute standard unweighted ROC-AUC to avoid array allocations & overhead.
    if weight_trivial ( sample_weight ) :
        return roc_auc_score ( y_true , y_pred , **kwargs)

    # 2. Handle negative weights if present
    if ( sample_weight < 0 ).any ():
        y_auc, w_auc = transform_weighted_target ( y_true , sample_weight ) 
        return roc_auc_score ( y_auc, y_pred , sample_weight = w_auc , **kwargs )

    # 3. Standard weighted ROC-AUC when all weights are non-negative (w >= 0)
    return roc_auc_score ( y_true , y_pred , sample_weight = sample_weight , **kwargs)

# =============================================================================
## @class ADVAL_base 
#  Base class for adversarial validation for the difference between two (weighted) datasets
class ADVAL_base (GoFnp):
    """ Base class for Adversarial Validation for the difference between two (weighted) datasets
    """
    def __init__ ( self               ,
                   nToys      = 400   ,
                   parallel   = False ,
                   silent     = False ,
                   progress   = True  ,
                   normalize  = True  , 
                   method     = "Adversarial Validation" , **params ) :
        
        n_splits = params.pop ( 'n_splits' , 5 ) 
        assert isinstance ( n_splits , int ) and 2 <= n_splits , "Invalid n_splits:%s" % n_splits
        
        self.__n_splits            = n_splits 
        self.__importance_features = {}
        self.__regularized         = None  
        
        GoFnp.__init__ ( self                  ,
                         nToys     = nToys     ,
                         parallel  = parallel  , 
                         silent    = silent    ,
                         progress  = progress  ,
                         normalize = normalize , 
                         method    = method    , **params )
                
    # ============================================================================
    @property
    def n_splits ( self ) :
        """`n_splits`: Number of splits for cross-validation"""
        return self.__n_splits 

    # ==================================================================================
    @property
    def config ( self ) :
        """`config`: Get all configuration parameters"""
        conf = {}
        conf.update ( super().config ) 
        conf [ 'n_splits' ] = self.n_splits
        for key, value in self.regularized.items() :  conf [ '%s [reg]' % key ] = value
        return conf
    
    # ==========================================================================
    @property 
    def regularized ( self ) :
        """`regularized`: regularized configuration
        - attention: copy is returned!
        """
        if self.__regularized is None :
            self.__regularized = self.regularization ()
        ## ATTENTION HERE: COPY IS RETURNED 
        return self.__regularized.copy() ## ATTENTION HERE: COPY IS RETURNED 

    # =========================================================================
    ##  Regularize 
    #   @code
    #   gof = ...
    #   gof.regularization () 
    #   @endcode
    @abc.abstractmethod 
    def regularization ( self ) :
        """ Create the regularized configuration
        >>> gof = ...
        >>> gof.regularization () 
        """
        params = {} 
        return params 
        
    # =========================================================================
    ## Are weights supported by this estimator?
    @property 
    def weights_supported ( self ) :
        """`weights_supported`: Are weights supported by this estimator?"""
        return True 

    # =========================================================================
    ## Are negative weights supported by this GoF estimator?
    @property 
    def negative_weights_supported ( self ) :
        """`negative_weghts_supported`: Are weights supported by this estimator?
        """
        return True 
        
    # =========================================================================
    ## Good for two-samples comparison?
    #  Can this estimator be used for comparison of two samples?
    @property 
    def two_samples ( self ) :
        """`two_samples`: Can this estimator be used for comparison of two samples?
        """
        return True 

    # ==========================================================================
    ## Feature importances
    @property 
    def importance_features ( self ) :
        """`importance_features`: Dictionary of feature importances"""
        return self.__importance_features
    
    # ==========================================================================
    ## Get formatted table of feature importances
    def importance_table  ( self ,
                            title  = '' ,
                            prefix = '' ,
                            style  = '' ) : 
        """ Get formatted table of feature importances
        """
        rows  = [ ( 'Feature/#' , 'Importance [%]' ) ]
        rows += [ ( str ( feature ) , '%.1f' % gain ) for feature, gain in self.importance_features.items () ] 
        title = title if title else "%s importance" % typename ( self )
        import ostap.logger.table as T
        return T.table ( rows               ,
                         title     = title  ,
                         prefix    = prefix ,
                         alignment = 'cc'   ,  
                         style     = style  )
    
    # ==================================================================================
    ## Train the model and make predictions
    #  @code
    #  gof = ...
    #  pred , imps = gof.work ( .... ) 
    #  @endcode 
    #  @return predictions and feature importances
    @abc.abstractmethod 
    def work ( self    ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        """ Train the model and make predictions
        - return predictions and feature importances
        >>> gof  = ...
        >>> pred , imps = gof.work ( .... ) 
        """
        return NotImplemented

    # ========================================================================
    ## Checks if dataset requires strong regularization for GoF / ADVAL
    def use_strong_regularization ( self , X ) :
        """ Checks if dataset size is small relative to feature count,
        or if low dimensionality (1D-3D) requires strict tree depth
        smoothing to prevent noise overfitting in Goodness-of-Fit tests.
        """
        ns = num_samples  ( X )
        nf = num_features ( X )

        # Always enforce strong regularization for low-dimensional spaces (1D-3D)
        if nf <= 3 : return True

        # Standard check for small sample sizes relative to feature space dimension
        return ns < max ( 300 if nf <= 3 else 1000 , nf * 50 )
    
    # ==========================================================================
    ## Calculate t-value for two non-structured (weighted) datasets 
    #   @param data1   the first dataset
    #   @param data2   the second dataset
    #   @param weight1 the first array of weights 
    #   @param weight2 the second array of weights
    #   tvalue is defined as \f$ 100 \times \left( 1 - 2 \times AUC \right)^2 \f$ 
    def tvalue ( self               ,
                 data1              ,
                 data2              ,  * , 
                 weight1    = None  ,
                 weight2    = None  ,
                 normalize  = True  ,
                 importance = False ) :
        """ Calculate t-value metric for two datasets under cross-validation.
        
        Parameters:
        -----------
        data1, data2 : array-like
            First and second dataset arrays (e.g., Train vs Test).
        weight1, weight2 : array-like, optional
            Sample weights corresponding to data1 and data2 (supports signed weights).
        normalize : bool, optional
            If True, applies pooled normalization to datasets.
        importance : bool, optional
            If True, aggregates feature importances across cross-validation folds.

        Returns:
        --------
        tv : float
            Calculated t-value metric: 100 * (1 - 2 * AUC)^2
        """
        data1, data2 = self.unpack ( data1 , data2 )
        sh1, sh2     = data1.shape , data2.shape
        assert len ( sh1 ) == 2 and len ( sh2 ) == 2 and sh1 [ 1 ] == sh2 [ 1 ] and sh1 [ 0 ] and sh2 [ 0 ], \
            f"Invalid dataset shapes: {sh1}, {sh2}"
        
        ## Apply normalization if requested
        if self.normalize and normalize :
            uds1 , uds2 = normalize_pooled ( data1 , data2 )
            return self.tvalue ( uds1       , uds2       ,
                                 weight1    = weight1    ,
                                 weight2    = weight2    ,
                                 normalize  = False      ,
                                 importance = importance )
            
        N1, N2 = sh1 [ 0 ] , sh2 [ 0 ]
        N      = N1 + N2

        # Combine datasets and labels using fast NumPy operations
        X = numpy.vstack ( [ data1 , data2 ] )
        Y = numpy.zeros  ( N , dtype = numpy.int8 )
        Y [ : N1 ] = 1

        w1_trivial = weight_trivial ( weight1 )
        w2_trivial = weight_trivial ( weight2 )

        if w1_trivial and w2_trivial:
            weights = False
            W       = None
        else:
            weights = True
            w1      = numpy.ones ( N1 , dtype = numpy.float32 ) if w1_trivial else numpy.asarray ( weight1 , dtype = numpy.float32 )
            w2      = numpy.ones ( N2 , dtype = numpy.float32 ) if w2_trivial else numpy.asarray ( weight2 , dtype = numpy.float32 )

            # ==============================================================================
            ## ATTENTION! 
            # ==============================================================================
            sumw1   = numpy.sum ( w1 , dtype = numpy.float64 )
            sumw2   = numpy.sum ( w2 , dtype = numpy.float64 )            
            w2      = w2 * numpy.float32 ( sumw1 / sumw2 ) 
                        
            W       = numpy.concatenate ( [ w1 , w2 ] )

        ## Cross-validation loop
        from sklearn.model_selection import StratifiedKFold
        
        random_state = self.params.get ( 'random_state' )
        cv           = StratifiedKFold ( n_splits = self.n_splits , shuffle = True , random_state = random_state )
        oof_preds    = numpy.zeros ( N , dtype = numpy.float32 )
        importances  = numpy.zeros ( X.shape [ 1 ] , dtype = numpy.float32 ) if importance else None
        
        for fold, ( train_idx , val_idx ) in enumerate ( cv.split ( X , Y ) ) :
            X_train, Y_train = X[train_idx], Y[train_idx]
            X_val,   Y_val   = X[val_idx],   Y[val_idx]
            
            W_train = W[train_idx] if weights else None
            W_val   = W[val_idx]   if weights else None

            fold_predictions, fold_importances = self.work ( X_train , Y_train , W_train ,
                                                             X_val   , Y_val   , W_val   ,
                                                             importance = importance     )

            oof_preds [ val_idx ] = fold_predictions

            if importance and importances is not None and fold_importances is not None:
                importances += fold_importances

        # Aggregate Feature Importances
        if importance and importances is not None and 0 < numpy.sum ( importances ) :
            imp_sum = numpy.sum ( importances )
            importances *= 100.0 / imp_sum
            
            feature_names = [f'f{i}' for i in range(X.shape[1])]
            sorted_pairs = sorted(zip(feature_names, importances), key=lambda x: x[1], reverse=True)
            self.__importance_features = {feat: val for feat, val in sorted_pairs}

        # Calculate Global ROC-AUC safely using standalone utility
        auc_score = safe_roc_auc_score ( Y , oof_preds , sample_weight = W if weights else None )

        return tvalue_from_AUC ( auc_score )
            
# =======================================================================================
## @class ADVAL_LGBM
#  LightGBM-based class for Adversarial Validation
#  @see lightgbm
class ADVAL_LGBM (ADVAL_base) : 
    """ LightGBM-based class for Adversarial Validation for the difference between two (weighted) datasets
    - see lightgbm 
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  , **params ) :
            
        config = {  'objective'         : 'binary',
                    'metric'            : 'auc',               # Explicitly measure ROC-AUC score
                    # Speed and convergence
                    'learning_rate'     : 0.03,                # Standard learning rate for stable training
                    'n_estimators'      : DEFAULT_ESTIMATORS , # High upper limit (MUST be used with early stopping!)
                    # Tree complexity (allows model to catch subtle feature shifts)
                    'max_depth'         : 5,                   # Moderate tree depth
                    'num_leaves'        : 24,                  # Sufficient leaf capacity
                    'min_child_samples' : 50,                  # Small threshold to detect fine-grained patterns
                    # Mild regularization to prevent overfitting on random noise
                    'subsample'         : 0.8,                # Row subsampling (80% per tree)
                    'subsample_freq'    : 1,
                    'colsample_bytree'  : 0.8,                # Feature subsampling (80% per tree)
                    'reg_alpha'         : 0.1,                # Slight L1 regularization
                    'reg_lambda'        : 1.0,                # Slight L2 regularization
                    ##  
                    'n_jobs'            : -1,                 # Utilize all CPU cores
                    'verbosity'         : -1                  # Suppress warning output
                }
        
        # Check parallel processing setup
        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "Parallel processing is switched OFF! (OMP/MKL/OPENBLAS)_NUM_THREADS" ) 
            parallel = False 
        
        ## Configure CPU jobs allocation
        params [ 'n_jobs' ] = 1 if parallel else num_jobs ( params , numcpu () - 1 )
        
        config.update ( params ) 
        ## 
        import lightgbm as LightGBM
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              parallel  = parallel ,
                              silent    = silent   , 
                              progress  = progress ,
                              normalize = False    ,
                              method    = "Adversarial Validation/LightGBM" , **config   ) 

    # =========================================================================
    ### get regularized configuraton
    #   @code
    #   gof = ...
    #   params = gof.regularization () 
    #   @endcode
    def regularization ( self ) : 
        """ Get regularized configuration
        >>> gof = ...
        >>> params = gof.regularized_parameters() 
        """
        params = {}
        
        # Enforce max_depth = 2 for 1D-3D datasets to prevent step-like discretization noise
        max_depth  = min ( 3 , self.params.get ( 'max_depth' , 5 ) )        
        num_leaves = 2 ** max_depth - 1
        
        params [ 'max_depth'         ] = max_depth
        params [ 'num_leaves'        ] = min ( num_leaves , self.params.get ( 'num_leaves'        , 31 ) )
        params [ 'min_child_samples' ] = max ( 100        , self.params.get ( 'min_child_samples' , 50 ) )
        params [ 'colsample_bytree'  ] = 1.0
        params [ 'reg_alpha'         ] = max ( 1.0        , self.params.get ( 'reg_alpha'         , 0.1 ) )
        params [ 'reg_lambda'        ] = max ( 5.0        , self.params.get ( 'reg_lambda'        , 1.0 ) )
        params [ 'learning_rate'     ] = min ( 0.02       , self.params.get ( 'learning_rate'     , 0.03 ) )
        
        return params 
            
    # ==================================================================================
    ## Train the LightGBM model and make predictions
    def work ( self    ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        """ Train LightGBM model on fold data and return validation predictions.
        """
        
        # =====================================================================
        ## 0. Parameters Preparation & Regularization Enforcement
        # =====================================================================
        params = {}
        params.update ( self.params )

        num_boost_round = params.pop ( 'num_boost_round' , None ) or params.pop ( 'n_estimators' , None ) or DEFAULT_ESTIMATORS 
        if isinstance ( num_boost_round , int ) and 10 < num_boost_round < 10000 : pass
        else : num_boost_round = DEFAULT_ESTIMATORS 

        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , 10 )
        if isinstance ( early_stopping_rounds , int ) and 1 < early_stopping_rounds < num_boost_round : pass
        else : early_stopping_rounds = 10

        # Enforce strong regularization for GoF (especially strict for 1D-3D features)
        if self.use_strong_regularization ( X_train ) :

            ## get the regularized parameters 
            regpars = self.regularized
            
            # Enforce max_depth = 2 for 1D-3D datasets to prevent step-like discretization noise
            nf         = num_features ( X_train )
            max_depth  = min ( 2 if nf <= 3 else 3       , regpars.get ( 'max_depth' , 5 ) )
            num_leaves = 2 ** max_depth - 1
            
            regpars [ 'max_depth'  ] = max_depth
            regpars [ 'num_leaves' ] = min ( num_leaves , regpars.get ( 'num_leaves', 31 ) )
            
            params.update ( regpars )
            
            num_boost_round        = min ( MAX_REGULARIZED_ESTIMATORS , num_boost_round       )
            early_stopping_rounds  = min (                         10 , early_stopping_rounds )
            
            
        import lightgbm as LightGBM
        
        callbacks = []
        if 0 < early_stopping_rounds :
            callbacks.append ( LightGBM.early_stopping ( stopping_rounds = early_stopping_rounds , verbose = False ) )

        # =====================================================================
        # 1. Transform targets and weights for train & val sets
        # =====================================================================
        Y_tr, W_tr = transform_weighted_target ( Y_train , W_train )
        Y_va, W_va = transform_weighted_target ( Y_val   , W_val   )
        
        # =====================================================================
        ## 2. Create native LightGBM Datasets
        # =====================================================================
        train_data = LightGBM.Dataset ( X_train , label = Y_tr , weight = W_tr , free_raw_data = False )
        val_data   = LightGBM.Dataset ( X_val   , label = Y_va , weight = W_va , free_raw_data = False , reference = train_data )

        # =====================================================================
        ## 3. Model Training
        # =====================================================================
        model = LightGBM.train ( params          = params          ,
                                 train_set       = train_data      ,
                                 num_boost_round = num_boost_round ,
                                 valid_sets      = [ val_data ]    ,
                                 callbacks       = callbacks       )

        # =====================================================================
        ## 4. Predict and restore predictions to original target probability space
        # =====================================================================
        raw_predictions = model.predict ( X_val , num_iteration = model.best_iteration )
        predictions     = invert_if_negative_weight ( raw_predictions , W_val )
        
        # 5. Extract Feature Importances (Gain)
        imps = model.feature_importance ( importance_type = 'gain') if importance else None
        
        return predictions, imps

# =======================================================================================
## @class ADVAL_XGB
#  XGBoost-based class for Adversarial Validation
#  @see xgboost 
class ADVAL_XGB (ADVAL_base) : 
    """ XGBoost-based class for Adversarial Validation for the difference between two (weighted) datasets
    - see xgboost
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  , **params ) :

        config = {  'objective'        : 'binary:logistic',   # Standard binary classification
                    'eval_metric'      : 'logloss',           # LogLoss objective and evaluation
                    'tree_method'      : 'hist',              # Fast histogram-based tree building algorithm
                    # Speed and step size
                    'learning_rate'    : 0.03,                # Standard step size
                    'n_estimators'     : DEFAULT_ESTIMATORS , # High upper bound
                    # Tree complexity
                    'max_depth'        : 5,                   # Equivalent to max_depth=5 in LightGBM
                    'max_leaves'       : 24,                  # Limit total leaves per tree
                    'min_child_weight' : 1.0,                 # Minimum sum of instance weight needed in a child
                    # Regularization
                    'subsample'        : 0.8,                 # Row subsampling (80% per tree)
                    'colsample_bytree' : 0.8,                 # Feature subsampling (80% per tree)
                    'alpha'            : 0.1,                 # L1 regularization
                    'lambda'           : 1.0,                 # L2 regularization
                    ##
                    'n_jobs'           : -1,                  # Utilize all CPU threads
                    'verbosity'        : 0                    # Suppress warning output
                }

        ## Configure CPU jobs allocation
        params [ 'n_jobs' ] = 1 if parallel else num_jobs ( params , numcpu() - 1 )
        
        config.update ( params ) 

        import xgboost as XGBoost             
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              parallel  = parallel ,
                              silent    = silent   , 
                              progress  = progress ,
                              normalize = False    ,
                              method    = "Adversarial Validation/XGBoost" , **config ) 

    # =========================================================================
    ##  Regularize 
    #   @code
    #   gof = ...
    #   gof.regularization () 
    #   @endcode
    def regularization ( self ) :
        """ Create the regularized configuration
        >>> gof = ...
        >>> gof.regularization () 
        """
        params = {} 
        
        # Enforce max_depth = 2 for 1D-3D datasets to prevent step-like discretization noise
        max_depth  = min ( 3 , self.params.get ( 'max_depth' , 5 ) )        
   
        params [ 'max_depth'         ] = max_depth 
        params [ 'min_child_weight'  ] = max ( 20   , self.params.get ( 'min_child_weight' , 5 ) )
        params [ 'colsample_bytree'  ] = 1.0
        params [ 'alpha'             ] = max ( 0.5  , self.params.get ( 'alpha'  , 0.1 ) )
        params [ 'lambda'            ] = max ( 2.0  , self.params.get ( 'lambda' , 1.0 ) )
        params [ "learning_rate"     ] = min ( 0.02 , self.params.get ( "learning_rate" , 0.03 ) )

        return params


    ## Train the XGBoost model and make predictions
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        """ Train XGBoost model on fold data and return validation predictions.
        
        Handles negative sample weights using transform_weighted_target and 
        restores predictions using invert_if_negative_weight.
        """
        import xgboost as XGBoost             

        # ===============================================================================
        # 1. Transform targets and weights for train & val sets (handling w < 0)
        # ===============================================================================
        Y_tr, W_tr = transform_weighted_target ( Y_train , W_train )
        Y_va, W_va = transform_weighted_target ( Y_val   , W_val   )

        # ===============================================================================
        # 2. Create native XGBoost DMatrix objects
        # ===============================================================================
        dtrain = XGBoost.DMatrix ( X_train , label = Y_tr , weight = W_tr )
        dval   = XGBoost.DMatrix ( X_val   , label = Y_va , weight = W_va )

        # Define evals list for early stopping monitoring
        evals  = [ ( dtrain , 'train' ) , ( dval , 'val' ) ]

        params = {}
        params.update ( self.params )
        
        num_boost_round = params.pop ( 'num_boost_round' , None ) or params.pop ( 'n_estimators' ,  None  ) or DEFAULT_ESTIMATORS 
        if isinstance ( num_boost_round , int ) and 10 < num_boost_round < 10000 : pass 
        else : num_boost_round = DEFAULT_ESTIMATORS 
                
        early_stopping_rounds = params.pop  ( 'early_stopping_rounds' , 10 )
        if isinstance ( early_stopping_rounds , int ) and 1 < early_stopping_rounds < num_boost_round : pass
        else : early_stopping_rounds = 10 
        
        if self.use_strong_regularization ( X_train ) :
            
            ## get the regularized parameters 
            regpars = self.regularized 
            
            nf         = num_features ( X_train )
            max_depth  = min ( 2 if nf <= 3 else 3 , regpars.get ( 'max_depth' , 5 ) )
            
            regpars [ 'max_depth'  ] = max_depth

            params.update ( regpars ) 

            num_boost_round        = min ( MAX_REGULARIZED_ESTIMATORS , num_boost_round       )
            early_stopping_rounds  = min (                         10 , early_stopping_rounds )

        # ===============================================================================
        # 3. Model Training
        # ===============================================================================
        model = XGBoost.train ( params                = params,
                                dtrain                = dtrain,
                                num_boost_round       = num_boost_round,
                                evals                 = evals,
                                early_stopping_rounds = early_stopping_rounds,
                                verbose_eval          = False )    

        # ===============================================================================
        # 4. Predict and restore predictions to original target probability space
        # ===============================================================================
        raw_predictions = model.predict ( dval )
        predictions     = invert_if_negative_weight ( raw_predictions , W_val )
        
        # ===============================================================================
        # 5. Extract Feature Importances (Gain)
        # ===============================================================================
        if importance :
            score_dict    = model.get_score ( importance_type = 'gain' )
            
            # Map feature names correctly using num_features helper
            if hasattr ( X_train , 'columns' ) :
                feature_names = X_train.columns.tolist()
            else :
                n_feats       = num_features ( X_train )
                feature_names = [ f'f{i}' for i in range ( n_feats ) ]

            imps = numpy.array ( [ score_dict.get ( col , 0.0 ) for col in feature_names ] , dtype = numpy.float32 )
        else :
            imps = None

        return predictions, imps

# ======================================================================================
## @class ADVAL_CATB
#  CatBoost-based class for Adversarial Validation
#  @see catboost 
class ADVAL_CATB (ADVAL_base) : 
    """ CatBoost-based class for Adversarial Validation for the difference between two (weighted) datasets
    - see CatBoost
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  , **params ) :
        
        config = {  'loss_function'         : 'Logloss',           # Standard binary classification loss
                    'eval_metric'           : 'Logloss',           # Metric monitored for early stopping
                    # Speed and convergence
                    'learning_rate'         : 0.03,                # Standard step size
                    'n_estimators'          : DEFAULT_ESTIMATORS , # Maximum boosting rounds
                    'early_stopping_rounds' : 20 ,
                    # Tree complexity
                    'depth'                 : 5,                   # Moderate depth (CatBoost trees are symmetric)
                    'min_data_in_leaf'      : 20,                  # Minimum samples per leaf
                    # Regularization & Subsampling
                    'l2_leaf_reg'           : 3.0,                 # L2 regularization for leaf values
                    'subsample'             : 0.8,                 # Row subsampling ratio
                    'bootstrap_type'        : 'Bernoulli',         # Enables subsampling
                    ##   
                    'thread_count'          : -1,                  # Utilize all CPU cores
                    'verbose'               : False
            }

        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "Parallel processing is switched OFF! (OMP/MKL/OPENBLAS)_NUM_THREADS" ) 
            parallel = False 
            
        ## Configure CPU threads allocation
        params [ 'thread_count' ] = 1 if parallel else num_jobs ( params , numcpu () - 1 )          
        
        if   'random_seed'  in params :                            params.pop ( 'random_state'      )
        elif 'random_state' in params : params [ 'random_seed' ] = params.pop ( 'random_state' , 42 )
        
        config.update ( params )
        ## 
        import catboost as CatBoost 
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              parallel  = parallel ,
                              silent    = silent   , 
                              progress  = progress ,
                              normalize = False    ,
                              method    = "Adversarial Validation/CatBoost" , **config   ) 


    # =========================================================================
    ##  Regularize 
    #   @code
    #   gof = ...
    #   gof.regularization () 
    #   @endcode
    def regularization ( self ) :
        """ Create the regularized configuration
        >>> gof = ...
        >>> gof.regularization () 
        """
        params = {} 
        
        # Enforce max_depth = 2 for 1D-3D datasets to prevent step-like discretization noise
        max_depth  = min ( 3 , params.pop ( 'depth' , 5 ) )
                
        params [ 'depth'             ] = max_depth  
        params [ 'min_data_in_leaf'  ] = max ( 20   , params.pop ( 'min_data_in_leaf' , 5 ) )
        params [ 'rsm'               ] = 1.0
        params [ 'l2_leaf_reg'       ] = max ( 5.0  , params.pop ( 'l2_leaf_reg'       , 3.0 ) )
        params [ "learning_rate"     ] = min ( 0.02 , params.get ( "learning_rate"     , 0.03 ) )

        return params 
        
    # ==================================================================================
    ## Train the CatBoost model and make predictions
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        """ Train CatBoost model on fold data and return validation predictions.
        
        Handles negative sample weights using transform_weighted_target and 
        restores predictions using invert_if_negative_weight.
        """
        import catboost as CatBoost
        
        # 1. Transform targets and weights for train & val sets (handling w < 0)
        Y_tr, W_tr = transform_weighted_target ( Y_train , W_train )
        Y_va, W_va = transform_weighted_target ( Y_val   , W_val   )

        # 2. Create native CatBoost Pool objects
        train_pool = CatBoost.Pool ( X_train , label = Y_tr , weight = W_tr )
        val_pool   = CatBoost.Pool ( X_val   , label = Y_va , weight = W_va )

        params = {}
        params.update ( self.params )
        
        iterations = params.pop ( 'iterations' , None ) or params.pop ( 'n_estimators' ,  None  ) or DEFAULT_ESTIMATORS 
        if isinstance ( iterations , int ) and 10 < iterations < 10000 : pass 
        else   : iterations = DEFAULT_ESTIMATORS  
                
        early_stopping_rounds = params.pop  ( 'early_stopping_rounds' , 10 )
        if isinstance ( early_stopping_rounds , int ) and 1 < early_stopping_rounds < iterations : pass
        else  : early_stopping_rounds = 10
        
        if self.use_strong_regularization ( X_train ) :

            regpars = self.regularized 
            
            nf = num_features ( X_train )            
            # Enforce max_depth = 2 for 1D-3D datasets to prevent step-like discretization noise
            max_depth  = min ( 2 if nf <= 3 else 3 , params.pop ( 'depth' , 5 ) )
            regpars [ 'depth'             ] = max_depth

            params.update ( regpars )
            
            iterations            = min ( MAX_REGULARIZED_ESTIMATORS , iterations            )
            early_stopping_rounds = min (                         10 , early_stopping_rounds )
        
        params [ 'iterations'            ] = iterations
        params [ 'early_stopping_rounds' ] = early_stopping_rounds  
        params [ 'use_best_model'        ] = True
                
        # 3. Model Training
        model = CatBoost.CatBoostClassifier ( **params ) 
        model.fit ( train_pool , eval_set = val_pool )

        # 4. Predict probabilities P(Class=1) and restore original probability space
        raw_predictions = model.predict_proba ( val_pool ) [ : , 1 ]
        predictions     = invert_if_negative_weight ( raw_predictions , W_val )

        # 5. Extract Feature Importances (PredictionValuesChange)
        imps = model.get_feature_importance(data=train_pool) if importance else None

        return predictions, imps

# =============================================================================
## @class ADVAL_HGBC
#  Class for adversarial validation based on HistGradientBoostingClassifier 
#  @see HistGradientBoostingClassifier 
class ADVAL_HGBC (ADVAL_base) : 
    """ HGBC-based class for Adversarial Validation for the difference between two (weighted) datasets
    @see HistGradientBoostingClassifier 
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  , **params ) :
        
        config = {  'loss'              : 'log_loss',          # Binary cross-entropy
                    'learning_rate'     : 0.03,                # Moderate step size
                    'max_iter'          : DEFAULT_ESTIMATORS , # Upper bound for boosting iterations
                    # Tree complexity
                    'max_depth'         : 5,                   # Restrict max tree depth
                    'max_leaf_nodes'    : 24,                  # Controls tree capacity
                    'min_samples_leaf'  : 20,                  # Minimum samples required in a leaf
                    # Regularization
                    'l2_regularization' : 0.1,                 # Mild L2 regularization
                    # Early stopping (built-in)
                    'early_stopping'    : True,                # Enables built-in early stopping
                    'n_iter_no_change'  : 20,                  # Number of iterations to wait for improvement
                    'scoring'           : 'roc_auc',           # Optimize directly for ROC-AUC
                }
        
        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "Parallel processing is switched OFF! (OMP/MKL/OPENBLAS)_NUM_THREADS" ) 
            parallel = False 
                
        config.update ( params )

        import sklearn.ensemble 
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              parallel  = parallel ,
                              silent    = silent   , 
                              progress  = progress , 
                              normalize = False    ,
                              method    = "Adversarial Validation/HGBC" , **config   ) 

    # =========================================================================
    ##  Regularize 
    #   @code
    #   gof = ...
    #   gof.regularization () 
    #   @endcode
    def regularization ( self ) :
        """ Create the regularized configuration
        >>> gof = ...
        >>> gof.regularization () 
        """
        params = {} 
        
        # Depth & leaf count limits
        max_depth = 3 
        params [ 'max_depth'         ] = min ( max_depth , self.params.get ( 'max_depth' , 5 ) or 5 )
        params [ 'max_leaf_nodes'    ] = 2 ** params [ 'max_depth' ] - 1  # 3 leaves for depth=2
        
        # Leaf size & L2 regularization
        orig_min_leaf = params.pop ( 'min_samples_leaf' , self.params.get( 'min_child_samples' , 20 ) )
        params [ 'min_samples_leaf'  ] = max ( 50  , orig_min_leaf )
        params [ 'l2_regularization' ] = max ( 5.0 , self.params.get ( 'l2_regularization' , 0.0 ) )
        
        # Learning rate & iterations
        params [ 'learning_rate'     ] = min ( 0.02 , self.params.get ( 'learning_rate' , 0.03 ) )

        return params 
            
    # ==================================================================================
    ## Train the HistGradientBoostingClassifier model and make predictions
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        """ Train sklearn's HistGradientBoostingClassifier on fold data.
        
        Handles negative sample weights using transform_weighted_target and 
        restores predictions using invert_if_negative_weight.
        """
        from sklearn.ensemble import HistGradientBoostingClassifier

        # 1. Transform targets and weights for train & val sets (handling w < 0)
        Y_tr, W_tr = transform_weighted_target ( Y_train , W_train )
        Y_va, W_va = transform_weighted_target ( Y_val   , W_val   )

        w_tr_arr = W_tr.values if hasattr ( W_tr , 'values' ) else W_tr


        params = {}
        params.update ( self.params )
        
        max_iter = params.pop ( 'max_iter' , None ) or params.pop ( 'n_estimators' , None ) or DEFAULT_ESTIMATORS 
        if isinstance ( max_iter , int ) and 10 < max_iter < 10000 : pass
        else  : max_iter = DEFAULT_ESTIMATORS 
        
        if self.use_strong_regularization ( X_train ) :

            regpars = self.regularized 
            
            nf = num_features ( X_train )
            # Depth & leaf count limits
            max_depth = 2 if nf <= 3 else 3
            regpars [ 'max_depth'         ] = min ( max_depth , regpars.get ( 'max_depth' , 5 ) or 5 )
            regpars [ 'max_leaf_nodes'    ] = 2 ** regpars[ 'max_depth' ] - 1  # 3 leaves for depth=2
            
            params.update ( regpars )
            
            max_iter                       = min ( max_iter , MAX_REGULARIZED_ESTIMATORS )
            
        params [ 'max_iter' ] = max_iter

        # 2. Model Training
        model = HistGradientBoostingClassifier ( **params )
        model.fit ( X_train , Y_tr , sample_weight = w_tr_arr )

        # 3. Predict probabilities P(Class=1) and restore original probability space
        raw_predictions = model.predict_proba(X_val)[:, 1]
        predictions     = invert_if_negative_weight ( raw_predictions , W_val )
        
        return predictions, None 
            
# =============================================================================
## @class ADVAL_GBC
#  Class for adversarial validation based on GradientBoostingClassifier 
#  @see GradientBoostingClassifier 
class ADVAL_GBC (ADVAL_base) : 
    """ GBC-based class for Adversarial Validation for the difference between two (weighted) datasets
    @see GradientBoostingClassifier 
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  , **params   ) :
        
        config = {  'loss'              : 'log_loss',          # Binary cross-entropy
                    'learning_rate'     : 0.05,                # Slightly higher rate for faster convergence
                    'n_estimators'      : DEFAULT_ESTIMATORS , # Sequential trees
                    # Tree complexity
                    'max_depth'         : 4,               # Restrict tree depth
                    'min_samples_split' : 20,              # Minimum samples to split an internal node
                    'min_samples_leaf'  : 10,              # Minimum samples required in a leaf node
                    # Regularization / Subsampling
                    'subsample'         : 0.8,             # Stochastic Gradient Boosting
                    'max_features'      : 'sqrt',          # Feature subsampling per split
                }
        
        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "Parallel processing is switched OFF! (OMP/MKL/OPENBLAS)_NUM_THREADS" ) 
            parallel = False 
        
        config.update ( params ) 

        import sklearn.ensemble 
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              parallel  = parallel ,
                              silent    = silent   , 
                              progress  = progress ,
                              normalize = False    ,
                              method    = "Adversarial Validation/GBC" , **config   ) 
        
    # =========================================================================
    ##  Regularize 
    #   @code
    #   gof = ...
    #   gof.regularization () 
    #   @endcode
    def regularization ( self ) :
        """ Create the regularized configuration
        >>> gof = ...
        >>> gof.regularization () 
        """
        params = {} 

        # Enforce max_depth = 2 for 1D-3D datasets to prevent step-like discretization noise
        max_depth  = min ( 3 , self.params.get ( 'max_depth' , 5 ) )        
        
        # 1. Depth & node constraints
        params [ 'max_depth'         ] = min ( max_depth , self.params.get ( 'max_depth' , 3 ) )
        
        orig_min_leaf  = params.pop ( 'min_samples_leaf'  , self.params.get ( 'min_child_samples' , 1 ) )
        orig_min_split = params.pop ( 'min_samples_split' , 2 )
        
        params [ 'min_samples_leaf'  ] = max ( 50  , orig_min_leaf  )
        params [ 'min_samples_split' ] = max ( 100 , orig_min_split )
        
        # 2. Pruning & Subsampling
        params [ 'ccp_alpha'         ] = max ( 0.01 , self.params.get ( 'ccp_alpha' , 0.0 ) )
        params [ 'subsample'         ] = min ( 0.8  , self.params.get ( 'subsample' , 1.0 ) )
        
        # 3. Feature sampling & learning rate
        params [ 'max_features'      ] = 1.0
        params [ 'learning_rate'     ] = min ( 0.02 , self.params.get ( 'learning_rate' , 0.03 ) )
        
        return params 

    # ==================================================================================
    ## Train the GradientBoostingClassifier model and make predictions
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        """ Train sklearn's standard GradientBoostingClassifier on fold data.
        
        Handles negative sample weights using transform_weighted_target and 
        restores predictions using invert_if_negative_weight.
        """
        from sklearn.ensemble import GradientBoostingClassifier

        # 1. Transform targets and weights for train & val sets (handling w < 0)
        Y_tr, W_tr = transform_weighted_target(Y_train, W_train)
        Y_va, W_va = transform_weighted_target(Y_val, W_val)

        w_tr_arr = W_tr.values if hasattr(W_tr, 'values') else W_tr

        params = {}
        params.update ( self.params )
        
        n_estimators = params.pop ( 'num_boost_round' , None ) or params.pop ( 'n_estimators' , None ) or DEFAULT_ESTIMATORS 
        if isinstance ( n_estimators , int ) and 10 < n_estimators < 10000 : pass
        else : n_estimators = DEFAULT_ESTIMATORS 
        
        if self.use_strong_regularization ( X_train ) :

            regpars = self.regularized 

            nf = num_features ( X_train )
            
            # 1. Depth & node constraints
            max_depth = 2 if nf <= 3 else 3
            regpars [ 'max_depth'         ] = min ( max_depth , regpars.pop ( 'max_depth' , 3 ) )

            params.update ( regpars )
            
            n_estimators                    = min ( MAX_REGULARIZED_ESTIMATORS  , n_estimators )
            
        params [ 'n_estimators' ] = n_estimators
        
        # 2. Model Training
        model = GradientBoostingClassifier ( **params )
        model.fit ( X_train , Y_tr , sample_weight = w_tr_arr)

        # 3. Predict probabilities P(Class=1) and restore original probability space
        raw_predictions = model.predict_proba(X_val)[:, 1]
        predictions     = invert_if_negative_weight ( raw_predictions , W_val )

        # 4. Extract Impurity-based Feature Importances
        imps = model.feature_importances_ if importance else None

        return predictions, imps

# =============================================================================
## @class ADVAL_RF
#  Class for adversarial validation based on RandomForestClassifier 
#  @see RandomForestClassifier 
class ADVAL_RF (ADVAL_base) : 
    """ RandomForest-based class for Adversarial Validation for the difference between two (weighted) datasets
    @see RandomForestClassifier 
    """
    def __init__ ( self             ,
                   nToys    = 500   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  , **params   ) :
  
        config = {  'n_estimators'      : DEFAULT_ESTIMATORS , # High tree count reduces variance
                    'n_jobs'            : -1,                  # Utilize all CPU cores
                    'criterion'         : 'log_loss',          # Optimizes cross-entropy
                    'max_depth'         : 10,                  # Restricted depth prevents overfitting
                    'min_samples_split' : 10,                  # Minimum samples required to split an internal node
                    'min_samples_leaf'  : 5,                   # Minimum samples required at a leaf node
                    'max_features'      : 'sqrt',              # Feature subsampling ratio per split
                    'bootstrap'         : True,                # Enables bootstrap sampling
                    'max_samples'       : 0.8,                 # Subsample 80% of rows per tree
                }  

        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "Parallel processing is switched OFF! (OMP/MKL/OPENBLAS)_NUM_THREADS" ) 
            parallel = False 

        ## Configure CPU jobs allocation
        params [ 'n_jobs' ] = 1 if parallel else num_jobs ( params , numcpu() - 1  ) 
        
        config.update ( params ) 

        import sklearn.ensemble 
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              parallel  = parallel ,
                              silent    = silent   , 
                              progress  = progress ,
                              normalize = False    ,                              
                              method    = "Adversarial Validation/RandomForest" , **config  )
        
    # =========================================================================
    ##  Regularize 
    #   @code
    #   gof = ...
    #   gof.regularization () 
    #   @endcode
    def regularization ( self ) :
        """ Create the regularized configuration
        >>> gof = ...
        >>> gof.regularization () 
        """
        params = {} 
        
        params [ 'max_depth'         ] = min ( 3 , self.params.get ( 'max_depth' , 3 ) or 3 )
        params [ 'max_leaf_nodes'    ] = 2 ** params [ 'max_depth' ]  # max 4 leaves for depth=2
        
        orig_min_leaf  = self.params.get ( 'min_samples_leaf'  , params.get ( 'min_child_samples' , 1 ) )
        orig_min_split = self.params.get ( 'min_samples_split' , 2 )
        
        params [ 'min_samples_leaf'  ] = max ( 50  , orig_min_leaf  )
        params [ 'min_samples_split' ] = max ( 100 , orig_min_split )
        
        # 2. Pruning, Sampling & Bootstrap
        params [ 'ccp_alpha'         ] = max ( 0.01 , self.params.get ( 'ccp_alpha' , 0.0 ) )
        params [ 'bootstrap'         ] = True
        params [ 'max_samples'       ] = min ( 0.8  , self.params.get ( 'max_samples' , 0.8 ) )
        params [ 'max_features'      ] = 1.0

        return params 
        
    # ==================================================================================
    ## Train the RandomForestClassifier model and make predictions
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        """ Train sklearn's RandomForestClassifier on fold data.
        
        Handles negative sample weights using transform_weighted_target and 
        restores predictions using invert_if_negative_weight.
        """
        from sklearn.ensemble import RandomForestClassifier

        # 1. Transform targets and weights for train & val sets (handling w < 0)
        Y_tr, W_tr = transform_weighted_target ( Y_train , W_train )
        Y_va, W_va = transform_weighted_target ( Y_val   , W_val   )

        w_tr_arr = W_tr.values if hasattr(W_tr, 'values') else W_tr


        params = {}
        params.update ( self.params )
        
        n_estimators = params.pop ( 'n_estimators' , None ) or params.pop ( 'num_boost_round' , None ) or DEFAULT_ESTIMATORS 
        if isinstance ( n_estimators , int ) and 10 < n_estimators < 10000 : pass
        else : n_estimators = DEFAULT_ESTIMATORS 
        
        if self.use_strong_regularization ( X_train ) :

            regpars = self.regularized
            
            nf = num_features ( X_train )
            
            # 1. Depth & Leaf Constraints
            max_depth = 2 if nf <= 3 else 3
            regpars [ 'max_depth'         ] = min ( max_depth , regpars.get ( 'max_depth' , 3 ) or 3 )
            regpars [ 'max_leaf_nodes'    ] = 2 ** regpars [ 'max_depth' ]  # max 4 leaves for depth=2

            params.update ( regpars )
            
            n_estimators                   = min ( MAX_REGULARIZED_ESTIMATORS , n_estimators )
            
        params [ 'n_estimators' ] = n_estimators
        
        # 2. Model Training
        model = RandomForestClassifier ( **params )
        model.fit(X_train, Y_tr, sample_weight=w_tr_arr)

        # 3. Predict probabilities P(Class=1) and restore original probability space
        raw_predictions = model.predict_proba(X_val)[:, 1]
        predictions     = invert_if_negative_weight(raw_predictions, W_val)

        # 4. Extract Impurity-based Feature Importances
        imps = model.feature_importances_ if importance else None

        return predictions, imps

# =============================================================================
## @class ADVAL_TORCH
#  Class for adversarial validation based on PyTorch
#  @see PyTorch
class ADVAL_TORCH (ADVAL_base) : 
    """ PyTorch-based class for Adversarial Validation for the difference between two (weighted) datasets
    @see PyTorch
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  , **params   ) :

        config =  { 'epochs'                : 100   ,
                    'batch_size'            : 256   ,
                    'lr'                    : 1.e-3 ,
                    'weight_decay'          : 1.e-4 ,
                    'early_stopping_rounds' :  15   , 
                    'dropout'               : 0.1   }

        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "Parallel processing is switched OFF! (OMP/MKL/OPENBLAS)_NUM_THREADS" ) 
            parallel = False 

        params [ 'n_jobs' ] = 1 if parallel else num_jobs ( params , numcpu() - 1  ) 
        
        config.update ( params ) 

        import torch as Torch 
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              parallel  = parallel ,
                              silent    = silent   , 
                              progress  = progress ,
                              normalize = True     , ## ATTENTION: Normalize MUST be True for PyTorch!
                              method    = "Adversarial Validation/PyTorch" , **config  ) 

    # =========================================================================
    ##  Regularize 
    #   @code
    #   gof = ...
    #   gof.regularization () 
    #   @endcode
    def regularization ( self ) :
        """ Create the regularized configuration
        >>> gof = ...
        >>> gof.regularization () 
        """
        params = {} 
        
        # 1. Network Architecture (1D-3D -> minimal representation)
        ##  params [ 'hidden_dim'    ] = 8 if nf <= 3 else 16
        params [ 'num_layers'    ] = 1   # Single hidden layer
        params [ 'dropout'       ] = max ( 0.3  , self.params.get ( 'dropout' , 0.1 ) )
        
        # 2. Optimizer Regularization (L2 Penalty)
        l2_reg = params.pop ( 'weight_decay'    , self.params.get ( 'reg_lambda' , 1e-4 ) )
        params [ 'weight_decay'  ] = max ( 1e-2 , l2_reg )
        
        # 3. Learning dynamics
        params [ 'learning_rate' ] = min ( 0.005 , self.params.get ( 'learning_rate' , 0.01 ) )
        params [ 'batch_size'    ] = max ( 128   , self.params.get ( 'batch_size'    , 64   ) )

        return params
    
    # ==================================================================================
    ## Train the PyTorch MLP model and make predictions
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        """ Train a PyTorch Multi-Layer Perceptron (MLP) on fold data.
        
        Handles negative sample weights using transform_weighted_target and 
        restores predictions using invert_if_negative_weight.
        """
        import copy
        import torch    as Torch 
        import torch.nn as NN 
        from   torch.utils.data import DataLoader, TensorDataset


        params = {}
        params.update ( self.params )
        
        epochs = params.pop ( 'epochs' , None ) or params.pop ( 'n_estimators' , None ) or 200
        if isinstance ( epochs , int ) and 10 < epochs < 10000 : pass
        else                                                   : epochs = 200

        if self.use_strong_regularization ( X_train ) :

            regpars = self.regularized 
            
            nf = num_features ( X_train )

            # 1. Network Architecture (1D-3D -> minimal representation)
            regpars [ 'hidden_dim'    ] = 8 if nf <= 3 else 16

            params.update ( regpars ) 
            epochs                     = min ( 50    , epochs )

        params [ 'epochs' ] = epochs
    
        
        # Training hyperparameters
        epochs       = params.get ( 'epochs'                , 100    )
        batch_size   = params.get ( 'batch_size'            , 256    )
        lr           = params.get ( 'lr'                    ,   1e-3 )
        weight_decay = params.get ( 'weight_decay'          ,   1e-4 )
        patience     = params.get ( 'early_stopping_rounds' ,  15    )        
        dropout_rate = params.get ( 'dropout'               ,   0.1  )

        n_features = num_features ( X_train )

        hidden_dim = 256 
        if   n_features <=   4 : hidden_dim =   4
        elif n_features <=   8 : hidden_dim =   8
        elif n_features <=  16 : hidden_dim =  16
        elif n_features <=  24 : hidden_dim =  24
        elif n_features <=  32 : hidden_dim =  32
        elif n_features <=  64 : hidden_dim =  64
        elif n_features <= 128 : hidden_dim = 128

        n_jobs = params.get ( 'n_jobs', -1 )
        if 0 < n_jobs : Torch.set_num_threads ( n_jobs )
        else          : Torch.set_num_threads ( 1      )
            
        device = Torch.device ('cuda' if Torch.cuda.is_available() else 'cpu' )

        # 1. Transform targets and weights for train & val sets (handling w < 0)
        Y_tr, W_tr = transform_weighted_target ( Y_train , W_train)
        Y_va, W_va = transform_weighted_target ( Y_val   , W_val)

        # Ensure float32 numpy arrays for PyTorch tensor compatibility
        X_tr_arr = numpy.asarray ( X_train , dtype = numpy.float32 )
        Y_tr_arr = numpy.asarray ( Y_tr    , dtype = numpy.float32 ) . reshape ( -1 , 1 )        
        W_tr_arr = ( numpy.ones_like ( Y_tr_arr ) if W_tr is None 
                     else numpy.asarray ( W_tr , dtype = numpy.float32 ) . reshape  (-1 , 1 ) )
        
        X_va_arr = numpy.asarray ( X_val , dtype = numpy.float32 )
        Y_va_arr = numpy.asarray ( Y_va  , dtype = numpy.float32 ) . reshape ( -1 , 1 )
        W_va_arr = ( numpy.ones_like ( Y_va_arr) if W_va is None 
                     else numpy.asarray ( W_va , dtype = numpy.float32 ) . reshape (-1 , 1 ) )

        # Fill NaNs with 0.0
        X_tr_arr = numpy.nan_to_num ( X_tr_arr , nan = 0.0 )
        X_va_arr = numpy.nan_to_num ( X_va_arr , nan = 0.0 )

        # 2. Data Loaders & Model Construction
        train_dataset = TensorDataset ( Torch.from_numpy ( X_tr_arr ) , 
                                        Torch.from_numpy ( Y_tr_arr ) , 
                                        Torch.from_numpy ( W_tr_arr ) )
        
        train_loader = DataLoader ( train_dataset , batch_size = batch_size , shuffle = True )

        val_x_tensor = Torch.from_numpy ( X_va_arr ).to ( device )
        val_y_tensor = Torch.from_numpy ( Y_va_arr ).to ( device )
        val_w_tensor = Torch.from_numpy ( W_va_arr ).to ( device )

        # Adaptive MLP Architecture
        if n_features == 1:
            model = NN.Sequential ( NN.Linear ( 1 , 2 ) ,
                                    NN.ReLU   ()        ,
                                    NN.Linear ( 2 , 1 ) ).to ( device )
        elif n_features <= 5:
            model = NN.Sequential ( NN.Linear ( n_features , 4 ),
                                    NN.ReLU   ()  ,
                                    NN.Linear ( 4 , 2 ) ,
                                    NN.ReLU   ()  ,
                                    NN.Linear ( 2 , 1 ) ).to ( device )
        else :         
            model = NN.Sequential ( NN.Linear      ( n_features   , hidden_dim ) ,
                                    NN.BatchNorm1d ( hidden_dim             ) ,
                                    NN.ReLU        (),
                                    NN.Dropout     ( dropout_rate           ) ,
                                    NN.Linear      ( hidden_dim   , hidden_dim // 2 ),
                                    NN.BatchNorm1d ( hidden_dim // 2        ) ,
                                    NN.ReLU        () ,
                                    NN.Dropout     ( dropout_rate           ),
                                    NN.Linear      ( hidden_dim // 2 , 1    ) ).to ( device )
            
            
        optimizer = Torch.optim.AdamW ( model.parameters() , lr = lr , weight_decay = weight_decay )
        
        # Helper function for weighted BCE loss
        def compute_weighted_bce ( logits , targets , weights ):
            loss = NN.functional.binary_cross_entropy_with_logits ( logits , targets , reduction = 'none' )
            return (loss * weights).mean()

        # 3. Training Loop with Early Stopping
        best_val_loss      = float('inf')
        best_model_weights = None
        patience_counter   = 0

        for epoch in range ( epochs ) :
            model.train()
            for bx, by, bw in train_loader:
                bx, by, bw = bx.to ( device ) , by.to ( device ) , bw.to ( device )
                optimizer.zero_grad()
                logits = model(bx)
                loss   = compute_weighted_bce(logits, by, bw)
                loss.backward()
                optimizer.step()

            # Validation evaluation
            model.eval()
            with Torch.no_grad():
                val_logits = model ( val_x_tensor )
                val_loss   = compute_weighted_bce ( val_logits , val_y_tensor, val_w_tensor).item()

            if val_loss < best_val_loss:
                best_val_loss      = val_loss
                best_model_weights = copy.deepcopy ( model.state_dict() )
                patience_counter   = 0
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    break

        # Restore best weights
        if best_model_weights is not None:
            model.load_state_dict ( best_model_weights )

        # 4. Predict probabilities P(Class=1) and restore target probability space
        model.eval()
        with Torch.no_grad():
            raw_predictions = Torch.sigmoid ( model ( val_x_tensor) ) .cpu().numpy().ravel()
        
        predictions = invert_if_negative_weight(raw_predictions, W_val)

        return predictions, None 

# =============================================================================
## @class ADVAL_KERAS
#  Class for adversarial validation based on Keras 
#  @see Keras
class ADVAL_KERAS (ADVAL_base) : 
    """ Keras-based class for Adversarial Validation for the difference between two (weighted) datasets
    @see Keras
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  , **params   ) :

        config =  { 'epochs'                : 100   ,
                    'batch_size'            : 256   ,
                    'lr'                    : 2.e-3 ,
                    'early_stopping_rounds' :  15   , 
                    'dropout'               : 0.1   }

        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "Parallel processing is switched OFF! (OMP/MKL/OPENBLAS)_NUM_THREADS" ) 
            parallel = False 

        params [ 'n_jobs' ] = 1 if parallel else num_jobs ( params , numcpu() - 1  ) 
        
        config.update ( params ) 

        import keras
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              parallel  = parallel ,
                              silent    = silent   , 
                              progress  = progress , 
                              normalize = True     , ## ATTENTION: Normalize MUST be True for Keras!
                              method    = "Adversarial Validation/Keras" , **config  ) 

    # =========================================================================
    ##  Regularize 
    #   @code
    #   gof = ...
    #   gof.regularization () 
    #   @endcode
    def regularization ( self ) :
        """ Create the regularized configuration
        >>> gof = ...
        >>> gof.regularization () 
        """
        params = {} 
        
        # 1. Network Architecture (minimal representation for 1D-3D)
        ## params [ 'hidden_dim'    ] = 8 if nf <= 3 else 16
        params [ 'num_layers'    ] = 1   # Single hidden layer
        params [ 'dropout'       ] = max ( 0.3  , self.params.get( 'dropout' , 0.1 ) )
        
        # 2. L2 Regularization
        l2_reg = params.pop ( 'l2_reg' , self.params.get ( 'reg_lambda' , 1e-4 ) )
        params [ 'l2_reg'        ] = max ( 1e-2 , l2_reg )
        
        # 3. Training Dynamics
        params [ 'learning_rate' ] = min ( 0.005 , self.params.get ( 'learning_rate' , 0.01 ) )
        params [ 'batch_size'    ] = max ( 128   , self.params.get ( 'batch_size'    , 64   ) )

        return params
    
    # ==================================================================================
    ## Train the Keras model and make predictions
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        """ Train an adaptive Keras 3 model on fold data.
        
        Handles negative sample weights using transform_weighted_target and 
        restores predictions using invert_if_negative_weight.
        """
        import keras as     Keras 
        from   keras import layers as Layers 


        params = {}
        params.update ( self.params )
        
        epochs = params.pop ( 'epochs' , None ) or params.pop ( 'n_estimators' , None ) or 200
        if isinstance ( epochs , int ) and 10 < epochs < 10000 : pass
        else                                                   : epochs = 200

        if self.use_strong_regularization ( X_train ) :

            regpars = self.regularized
            
            nf = num_features ( X_train )

            # 1. Network Architecture (minimal representation for 1D-3D)
            regpars [ 'hidden_dim'    ] = 8 if nf <= 3 else 16
            
            params.update ( regpars ) 
            epochs                     = min ( 50    , epochs )

        params [ 'epochs' ] = epochs
                
        n_features = num_features ( X_train )
        
        hidden_dim = 256 
        if   n_features <=   4 : hidden_dim =   4
        elif n_features <=   8 : hidden_dim =   8
        elif n_features <=  16 : hidden_dim =  16
        elif n_features <=  24 : hidden_dim =  24
        elif n_features <=  32 : hidden_dim =  32
        elif n_features <=  64 : hidden_dim =  64
        elif n_features <= 128 : hidden_dim = 128

        # Extract hyperparameters
        epochs       = params.get ( 'epochs'               , 100 )
        batch_size   = params.get ( 'batch_size'           , 256 )
        patience     = params.get ( 'early_stopping_rounds',  15 )
        dropout_rate = params.get ( 'dropout'              , 0.1 if n_features > 5 else 0.0 )
        lr           = params.get ( 'lr'                   , 1e-3 if n_features > 5 else 5e-3 )
        
        # 1. Transform targets and weights for train & val sets (handling w < 0)
        Y_tr, W_tr = transform_weighted_target ( Y_train , W_train )
        Y_va, W_va = transform_weighted_target ( Y_val   , W_val   )

        # Ensure float32 numpy arrays
        X_tr_arr = numpy.nan_to_num ( numpy.asarray ( X_train , dtype=numpy.float32 ).reshape(-1, n_features ), nan = 0.0 )
        Y_tr_arr = numpy.asarray ( Y_tr, dtype = numpy.float32 ) . reshape ( -1 , 1 )
        W_tr_arr = ( numpy.ones_like ( Y_tr_arr ) if W_tr is None 
                     else numpy.asarray( W_tr , dtype = numpy.float32 ) . reshape ( -1 , 1 ) )

        X_va_arr = numpy.nan_to_num ( numpy.asarray ( X_val , dtype = numpy.float32 ) . reshape ( -1 , n_features ), nan = 0.0 )
        Y_va_arr = numpy.asarray ( Y_va, dtype = numpy.float32 ) . reshape ( -1 , 1 )
        W_va_arr = ( numpy.ones_like  (Y_va_arr ) if W_va is None 
                     else numpy.asarray ( W_va, dtype = numpy.float32 ) . reshape (-1, 1 ) )
        
        # 2. Construct Keras Model Adaptive to N Features
        inputs = Keras.Input ( shape = ( n_features , ) )

        if n_features == 1:
            x       = Layers.Dense ( 2 , activation = 'relu'    )( inputs )
            outputs = Layers.Dense ( 1 , activation = 'sigmoid' )( x )
        elif n_features <= 5:
            x       = Layers.Dense ( hidden_dim , activation='relu' ) ( inputs )
            x       = Layers.Dense ( max ( 2, hidden_dim // 2) , activation = 'relu'    ) ( x )
            outputs = Layers.Dense ( 1                         , activation = 'sigmoid' ) ( x )
        else:
            x       = Layers.Dense ( hidden_dim ) ( inputs )
            x       = Layers.BatchNormalization() ( x )
            x       = Layers.Activation ('relu' ) ( x )
            if 0 < dropout_rate : x = Layers.Dropout ( dropout_rate ) ( x )
            x       = Layers.Dense ( hidden_dim // 2 ) ( x )
            x       = Layers.BatchNormalization()      ( x )
            x       = Layers.Activation('relu')        ( x )
            if 0 < dropout_rate : x = Layers.Dropout ( dropout_rate ) ( x )            
            outputs = Layers.Dense ( 1 , activation = 'sigmoid' ) ( x )

        model = Keras.Model(inputs=inputs, outputs=outputs)

        # 3. Compile Model
        optimizer = Keras.optimizers.AdamW ( learning_rate = lr )
        model.compile ( optimizer        = optimizer,
                        loss             = Keras.losses.BinaryCrossentropy(),
                        weighted_metrics = [] )
        
        # 4. Early Stopping Callback
        early_stopping = Keras.callbacks.EarlyStopping ( monitor              = 'val_loss' ,
                                                         patience             = patience   ,
                                                         restore_best_weights = True       ,
                                                         verbose              = 0          )
        
        # 5. Fit Model
        model.fit( X_tr_arr        ,
                   Y_tr_arr        ,
                   sample_weight   = W_tr_arr.ravel()   ,
                   validation_data = ( X_va_arr, Y_va_arr , W_va_arr.ravel() ),
                   epochs          = epochs             ,
                   batch_size      = batch_size         ,
                   callbacks       = [ early_stopping ] ,
                   verbose         = 0 )

        # 6. Predict probabilities P(Class=1) and restore target probability space
        raw_predictions = model.predict ( X_va_arr , batch_size = batch_size , verbose = 0 ).ravel()
        predictions     = invert_if_negative_weight(raw_predictions, W_val)

        return predictions, None 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    from ostap.stats.tools import ( hasLightGBM , hasXGBoost , hasCatBoost ,
                                    hasSkLearn  , hasPyTorch , hasKeras    )

    if not hasLightGBM ( False ) : logger.warning  ( "No LightGBM available!" ) 
    if not hasXGBoost  ( False ) : logger.warning  ( "No XGBoost  available!" ) 
    if not hasCatBoost ( False ) : logger.warning  ( "No CatBoost available!" ) 
    if not hasSkLearn  ( False ) : logger.warning  ( "No scikit-learn available!" ) 
    if not hasPyTorch  ( False ) : logger.warning  ( "No PyTorch available!" ) 
    if not hasKeras    ( False ) : logger.warning  ( "No Keras   available!" ) 
            
# =============================================================================
##                                                                      The END 
# =============================================================================

