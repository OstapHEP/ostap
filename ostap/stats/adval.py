#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/adversarial_validation.py
#  Implement adversarial validation (regression mode) to probe
#  the difference between two weighted datasets (supporting negative weights via sPlot)
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-04
# =============================================================================
""" Implement adversarial validation (regression mode) to probe
    the difference between two weighted datasets (supporting negative weights via sPlot)

 - see ADVAL_LGBM  , class based on LightGBM
 - see ADVAL_XGB   , class based on XGBoost 
 - see ADVAL_CATB  , class based on CatBoost 
 - see ADVAL_HGBC  , class based on HistGradientBoostingRegressor
 - see ADVAL_GBC   , class based on GradientBoostingRegressor
 - see ADVAL_RF    , class based on RandomForestRegressor 
 - see ADVAL_TORCH , class based on PyTorch (Regression)
 - see ADVAL_KERAS , class based on Keras (Regression)

As t-value: based on Weighted MSE / RMSE is used.
To estimate the p-value, permutations are used. 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2026-07-04"
__all__     = (
    'ADVAL_LGBM'  , 
    'ADVAL_XGB'   , 
    'ADVAL_CATB'  , 
    'ADVAL_HGBC'  , 
    'ADVAL_GBC'   , 
    'ADVAL_RF'    , 
    'ADVAL_TORCH' , 
    'ADVAL_KERAS' , 
)
# =============================================================================
from   collections            import defaultdict 
from   ostap.core.ostap_types import string_types
from   ostap.utils.core       import typename
from   ostap.stats.utils      import ( weight_trivial , nEff         , 
                                       num_samples    , num_features , 
                                       check_all      )
from   ostap.utils.basic      import numcpu, num_jobs, run_parallel 
from   ostap.stats.gof_np     import GoFnp
from   ostap.stats.counters   import EffCounter
from   ostap.logger.pretty    import nice_print 
from   sklearn.metrics        import mean_squared_error, r2_score 
import ostap.logger.symbols   as     S
import ROOT, numpy, abc, os   
# =============================================================================
# Logging setup 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.adval' )
else                       : logger = getLogger ( __name__ )
# =============================================================================
logger.debug ( 'Implement Adversarial Validation (Regression mode)' )
# =============================================================================
DEFAULT_ESTIMATORS         = 500
MAX_REGULARIZED_ESTIMATORS = 100
MAX_DEPTH                  =   5
# =============================================================================
method_LGBM  = 'AdVal/%s' % ( S.light_bulb           if S.light_bulb  else 'LightGBM' ) 
method_XGB   = 'AdVal/%s' % ( S.rocket               if S.rocket      else 'XGBoost'  ) 
method_CATB  = 'AdVal/%s' % ( S.cat_face             if S.cat_face    else 'CatBoost' ) 
method_HGBC  = 'AdVal/%s' % ( S.proportion + S.nabla if S.show        else 'HGBC'     ) 
method_GBC   = 'AdVal/%s' % ( S.nabla                if S.show        else 'GBC'      )
method_RF    = 'AdVal/%s' % ( S.tree                 if S.show        else 'RF'       ) 
method_TORCH = 'AdVal/%s' % ( S.flashlight           if S.flashlight  else 'TORCH'    ) 
method_KERAS = 'AdVal/%s' % ( S.postal_horn          if S.postal_horn else 'KERAS'    ) 
# =============================================================================
## Evaluate whether strong regularization is required for BDT-based models.
#  Evaluates feature dimension, Kish's effective sample size, and sPlot weight noise.
#
#  @param X Array-like feature matrix (n_samples, n_features).
#  @param W Optional sample weights (supports negative sPlot weights).
#  @return bool True if strong regularization should be applied, False otherwise.
# =============================================================================
def BDT_needs_regularization ( X , W = None ) :
    """ Determines whether strong regularization is necessary for tree-based 
        Adversarial Validation models.
    """
    nf   = num_features ( X )
    nraw = num_samples  ( X )
    neff = nEff         ( X , W )
    
    # Since W is class-balanced (sumw1 == sumw2), this global nEff 
    # mathematically equals the two-sample harmonic effective size:
    # 4 * (N_eff1 * N_eff2) / (N_eff1 + N_eff2).
    # It correctly bottlenecks capacity by the weaker sample.
    # 
    # 1. Check weight efficiency (nEff / nRaw)
    #    Efficiency < 65% indicates significant sPlot negative weight fluctuations
    eff = neff / float ( nraw ) if 0 < nraw else 0.0
    if eff < 0.65 :
        effp = 100 * eff 
        return 'eff[%.0f%%]<65%%' % effp 

    # 2. Low dimensionality (<= 3) requires regularization only under limited statistics
    threshold = 50000.0 
    if nf <= 3 and neff < threshold :
        th     = nice_print ( threshold , precision = 1 , width = 2 , with_sign = False )
        vv     = nice_print ( neff      , precision = 1 , width = 2 , with_sign = False ) 
        result = "nf[%d]<=3&neff[%s]<%s" % ( nf , vv , th )
        return result.replace ( ' ' , '' ) 

    # 3. Dimensionality-dependent minimum statistical threshold
    required_stats = 1000.0 * ( nf ** 1.5 )
    if neff < required_stats :
        rs     = nice_print ( required_stats , precision = 1 , width = 2 , with_sign = False )
        vv     = nice_print ( neff           , precision = 1 , width = 2 , with_sign = False ) 
        result = "neff[%s]<%s" % ( vv ,rs )
        return result.replace ( ' ' , '' ) 
        
    return ''

# =============================================================================
## Determine whether strong regularization is required specifically for Neural Networks
#  (PyTorch / Keras) when handling weighted or sPlot datasets.
#
#  Neural networks require significantly higher effective statistics and stricter 
#  bounds on weight dispersion to prevent overfitting on tabular sPlot data.
#
#  @param X Array-like feature matrix (n_samples, n_features).
#  @param W Optional sample weights (supports negative sPlot weights).
#  @return bool True if strong regularization should be applied, False otherwise.
# =============================================================================
def NN_needs_regularization ( X , W = None ):
    """ Determines whether strong regularization (e.g., dropout, weight decay) 
        is required specifically for Neural Networks handling weighted datasets.
    """
    nf   = num_features ( X )
    nraw = num_samples  ( X )
    neff = nEff         ( X , W )

    # 1. Strict absolute statistic threshold for Neural Networks
    threshold1 = 10000.0 
    if neff < threshold1 :
        th1    = nice_print ( threshold1 , precision = 1 , width = 2 , with_sign = False )
        vv     = nice_print ( neff       , precision = 1 , width = 2 , with_sign = False )         
        result ='neff[%s]<%s' % ( vv , th1 )
        return result.replace ( ' ' , '' ) 
        
    
    # 2. Higher statistical coverage required for low-dimensional spaces to ensure smooth boundaries
    threshold2 = 100000.0     
    if nf <= 3 and neff < threshold2 :
        th2    = nice_print ( threshold2 )
        vv     = nice_print ( neff           , precision = 1 , width = 2 , with_sign = False )                 
        result = 'nf[%d]<=3&neffp[%s]<%s' % ( nf , vv , th2 )
        return result.replace ( ' ' , '' ) 
        
    # 3. Non-linear event density threshold per feature
    required_stats = 3000.0 * ( nf ** 1.5 )
    if neff < required_stats :
        rs     = nice_print ( required_stats , precision = 1 , width = 2 , with_sign = False )
        vv     = nice_print ( neff           , precision = 1 , width = 2 , with_sign = False )         
        result = 'neff[%s]<%s'% ( vv , rs )
        return result.replace ( ' ' , '' ) 
        
    # 4. Weight efficiency and dispersion checks for non-trivial weights
    if not weight_trivial ( W ) :
        
        eff = neff / float ( nraw ) if 0 < nraw else 0.0
        if eff < 0.65 :
            effp = 100 * eff 
            return 'eff[%.1f%%]<65%%' % effp 

        w_arr        = numpy.asarray ( W     , dtype = numpy.float32 )
        mean_w       = numpy.mean    ( w_arr , dtype = numpy.float64 )
        w_dispersion = ( numpy.std   ( w_arr , dtype = numpy.float64 ) / mean_w ) if 0.0 < mean_w else 0.0
        
        # Enforce strict weight dispersion limit (> 1.5) for network convergence stability
        if w_dispersion > 1.5 : return 'w_dispersion[%.1f]>1.5'

    return '' 

# =============================================================================
## convert MSE to t-value 
def tvalue_from_MSE ( mse ) :
    """ Convert MST to tvalue"
    - MSE = 0.25  means no difference 
    """
    return 1 - 4 * float ( mse ) 

# =============================================================================
## Helper to handle negative weights for regression-based models via target inversion
def transform_weights_and_targets ( Y , W ) :
    """ Transform negative weights into positive ones by flipping the regression target 
        Y -> (1 - Y) when W < 0, preserving the physical meaning for MSE optimization.
    """
    if W is None: return Y, None
    
    w_arr = numpy.asarray ( W , dtype = numpy.float32 )
    signs = numpy.sign ( w_arr )
    signs [ signs == 0 ] = 1.0
    
    Y_mod = numpy.where ( signs < 0 , 1.0 - Y , Y )
    W_mod = numpy.abs   ( w_arr )
    
    return Y_mod, W_mod

# =============================================================================
## @class ADVAL_base 
#  Base class for adversarial validation (regression mode)
class ADVAL_base (GoFnp):
    """ Base class for Adversarial Validation (regression mode) for weighted datasets
    """
    def __init__ ( self               ,
                   nToys      = 400   ,
                   parallel   = True  ,
                   normalize  = True  ,
                   n_splits   = 5     , 
                   method     = "Adversarial Validation (Regression)" , **params ) :

        # =====================================================================
        if   not isinstance ( n_splits , int ) : raise TypeError ( "Invalid `n_splits` type  : %s" % typename ( n_splits ) )
        elif n_splits < 2                      : raise TypeError ( "Invalid `n_splits` value : %d" %            n_splits   )
        ## store it! 
        self.__n_splits = n_splits
        
        ## switch off parallel processing
        parallel = True if parallel else False 
        if parallel and not run_parallel ( parallel ) :
            logger.info ( '%s: (internal) parallel processing is OFF' % typename ( self ) )
            parallel = False
            
        ## (re)define n_jobs     
        params [ 'n_jobs'    ] = 1 if parallel else num_jobs ( params )        
        params [ 'parallel'  ] = parallel
        params [ 'normalize' ] = True if normalize else False 

        self.__n_splits            = n_splits 
        self.__importance_features = {}

        ## how often the regularization applied?
        self.__regularized         = defaultdict(int)

        GoFnp.__init__ ( self            ,
                         nToys  = nToys  ,
                         method = method , **params )
                
    # =========================================================================
    @property
    def n_splits ( self ) :
        """`n_splits`: Number of splits for cross-validation"""
        return self.__n_splits 

    # =========================================================================
    @property
    def config ( self ) :
        """`config`: Get all configuration parameters"""
        conf = {}
        conf.update ( super().config ) 
        conf [ 'n_splits'    ] = self.n_splits
        if self.regularized  :
            regs = {}
            regs.update ( self.regularized )
            conf [ 'regularized' ] = regs 
        return conf
    
    @property 
    def regularized ( self ) :
        return self.__regularized 

    @property 
    def weights_supported ( self ) :
        return True 

    @property 
    def two_samples ( self ) :
        return True

    @property 
    def importance_features ( self ) :
        return self.__importance_features

    ## prepare the table of importance features 
    def importance_table  ( self ,
                            title  = '' ,
                            prefix = '' ,
                            style  = '' ) :
        """ Prepare the table of importance features
        """
        
        rows  = [ ( 'Feature/#' , 'Importance [%]' ) ]
        rows += [ ( str ( feature ) , '%.2f' % float ( gain ) ) for feature, gain in self.importance_features.items () ] 
        title = title if title else "%s importance" % typename ( self )
        import ostap.logger.table as T
        return T.table ( rows               ,
                         title     = title  ,
                         prefix    = prefix ,
                         alignment = 'cc'   ,  
                         style     = style  )
    

    @abc.abstractmethod 
    def work ( self    ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        return NotImplemented

    # =========================================================================
    ## Parameters for strong regularization
    #  @code
    #  params gof.regularization ( params , n_features , n_samples ) 
    #  @endcode
    @abc.abstractmethod 
    def regularization ( self       ,
                         params     ,
                         n_features , 
                         n_samples  ) :
        
        """ Parameters for strong regularization
        >>> params = gof.regularization ( params , n_features , n_samples ) 
        """
        return NotImplemented
        
    def tvalue ( self               ,
                 data1              ,
                 data2              ,  * , 
                 weight1    = None  ,
                 weight2    = None  ,
                 normalize  = True  ,
                 importance = False ) :
        
        data1, data2 = self.unpack ( data1 , data2 )

        ## check input data 
        check_all ( data1 , data2 , weight1 , weight2 , typename ( self ) )
        
        if self.normalize and normalize :
            uds1 , uds2 = self.normalize_pooled ( data1 , data2 )
            return self.tvalue ( uds1       , uds2       ,
                                 weight1    = weight1    ,
                                 weight2    = weight2    ,
                                 normalize  = False      ,
                                 importance = importance )
            
        N1 = num_samples ( data1 )
        N2 = num_samples ( data2 ) 
        N  = N1 + N2

        X = numpy.vstack ( [ data1 , data2 ] )
        Y = numpy.zeros  ( N , dtype = numpy.float32 )
        Y [ : N1 ] = 1.0  # Test gets 1, Train gets 0

        w1_trivial = weight_trivial ( weight1 )
        w2_trivial = weight_trivial ( weight2 )

        if w1_trivial and w2_trivial:
            weights = False
            W       = None
        else:
            weights = True
            w1      = numpy.ones ( N1 , dtype = numpy.float32 ) if w1_trivial else numpy.asarray ( weight1 , dtype = numpy.float32 )
            w2      = numpy.ones ( N2 , dtype = numpy.float32 ) if w2_trivial else numpy.asarray ( weight2 , dtype = numpy.float32 )

            sumw1   = numpy.sum ( w1 , dtype = numpy.float64 )
            sumw2   = numpy.sum ( w2 , dtype = numpy.float64 )            
            w2      = w2 * numpy.float32 ( sumw1 / sumw2 ) 
                        
            W       = numpy.concatenate ( [ w1 , w2 ] )

        from sklearn.model_selection import KFold
        
        random_state = self.params.get ( 'random_state' )
        cv           = KFold ( n_splits = self.n_splits , shuffle = True , random_state = random_state )
        oof_preds    = numpy.zeros ( N , dtype = numpy.float32 )
        importances  = numpy.zeros ( X.shape [ 1 ] , dtype = numpy.float32 ) if importance else None
        
        for fold, ( train_idx , val_idx ) in enumerate ( cv.split ( X , Y ) ) :
            X_train, Y_train_f = X[train_idx], Y[train_idx]
            X_val,   Y_val_f   = X[val_idx],   Y[val_idx]
            
            W_train = W[train_idx] if weights else None
            W_val   = W[val_idx]   if weights else None

            fold_predictions, fold_importances = self.work ( X_train , Y_train_f , W_train ,
                                                             X_val   , Y_val_f   , W_val   ,
                                                             importance = importance     )

            oof_preds [ val_idx ] = fold_predictions

            if importance and importances is not None and fold_importances is not None:
                importances += fold_importances

        if importance and importances is not None and 0 < numpy.sum ( importances ) :
            imp_sum = numpy.sum ( importances )
            importances *= 100.0 / imp_sum
            
            feature_names = [f'f{i}' for i in range(X.shape[1])]
            sorted_pairs = sorted(zip(feature_names, importances), key=lambda x: x[1], reverse=True)
            self.__importance_features = {feat: val for feat, val in sorted_pairs}


        Y_eval, eval_weights = transform_weights_and_targets ( Y , W )
        ## mse_score = mean_squared_error ( Y , oof_preds , sample_weight = eval_weights )
        mse_score = mean_squared_error ( Y_eval , oof_preds , sample_weight = eval_weights )

        return tvalue_from_MSE ( mse_score )

    # ===================================================================================
    ## print regularized parameters in "no-silent" regime
    def report_regularization ( self , params = {} , **kwargs ) :
        """ print regularized paramters in "no-silent" regime
        """
        
        if self.silent               : return
        if not params and not kwargs : return
        
        rpars = {}
        rpars.update ( params , **kwargs )
        title = "%s regularized" % typename ( self )
        from ostap.logger.utils import map2table_ex
        logger.info ( '%s:\n%s' % ( title , map2table_ex ( rpars       , 
                                                           header      = ( 'Parameter' , 'type' , 'value' ) ,
                                                           alignment   = 'rcw'  , 
                                                           prefix      = "# "   ,
                                                           title       = title  ) ) ) 

# =======================================================================================
## @class ADVAL_LGBM (Regression)
class ADVAL_LGBM (ADVAL_base) : 
    def __init__ ( self             ,
                   nToys    = 400   , **params ) :
        
        config = {
            'objective'         : 'regression',
            'metric'            : 'rmse',
            'learning_rate'     :  0.03,
            'n_estimators'      : DEFAULT_ESTIMATORS ,
            'max_depth'         : MAX_DEPTH ,
            'num_leaves'        : 24   ,
            'min_child_samples' : 50   ,
            'subsample'         :  0.8 ,
            'subsample_freq'    :  1   ,
            'colsample_bytree'  :  0.8 ,
            'reg_alpha'         :  0.1 ,
            'reg_lambda'        :  1.0 ,
            'force_col_wise'    : True ,
            'max_bin'           :  127 , 
            'n_jobs'            : -11  ,  
            'verbosity'         : -1   ,
            'verbose'           : -1   ,
        }
            
        config.update ( params ) 
        
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              normalize = False    ,
                              method    = method_LGBM , **config   ) 

        if 'n_jobs' in self.params :
            self.params [ 'num_thread' ] = self.params.pop ( 'n_jobs' , 1 )
        if self.silent :
            self.params [ 'verbosity' ] = -1 
            self.params [ 'verbose'   ] = -1 

    # =========================================================================
    ## Parameters for strong regularization
    #  @code
    #  params = gof.regularization ( params , n_features , n_samples ) 
    #  @endcode
    def regularization ( self       ,
                         params     ,
                         n_features , 
                         n_samples  ) :
        
        """ Parameters for strong regularization
        >>> params = gof.regularization ( params , n_features , n_samples ) 
        """
        
        # ADAPTIVE DEPTH    
        max_depth  = min ( 6 , max ( 1 , n_features ) )
        
        # SAFE NUMBER OF LEAVES 
        num_leaves = max ( 2 , int ( 2 ** ( max_depth * 0.8 ) ) ) 
        
        params [ 'max_depth'         ] = max_depth
        params [ 'num_leaves'        ] = num_leaves
        
        # DYNAMIC LEAF MASS 
        params [ 'min_child_samples' ] = max ( 20 , int ( n_samples * 0.001 ) )
        params [ 'min_child_weight'  ] = 1e-3
        
        params [ 'colsample_bytree'  ] = 1.0
        params [ 'subsample'         ] = 1.0
        params [ 'reg_alpha'         ] = 0.0
        params [ 'reg_lambda'        ] = 0.0
        
        params [ 'learning_rate'     ] = min (  0.05 , params.get ( 'learning_rate' ,   0.05 ) )
        params [ 'max_bins'          ] = min ( 63    , params.get ( 'mas_bin'       , 127    ) )
        
        return params
    
    def work ( self    ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        
        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        Y_val_mod,   W_val_mod   = transform_weights_and_targets ( Y_val  , W_val   )

        params = {}
        params.update ( self.params )

        nf =  num_features ( X_train )
        ns =  num_samples  ( X_train )

        num_boost_round       = params.pop ( 'num_boost_round' , None ) or params.pop ( 'n_estimators' , None ) or DEFAULT_ESTIMATORS 
        early_stopping_rounds = params.pop ( 'early_stopping_rounds' ,  10 )

        reg_case = BDT_needs_regularization ( X_train , W_train )        
        if reg_case :
            self.regularized [ reg_case ] += 1
            
            ## update parameters 
            params.update ( self.regularization ( params , nf , ns ) )
            
            num_boost_round        = min ( 20 if 1 == nf else 50 , num_boost_round )
            early_stopping_rounds  = max ( 10 , early_stopping_rounds )
            
            ## print regularized parameters in "no-silent" regime
            self.report_regularization ( params                ,
                                         reg_case              = reg_case              , 
                                         num_features          = nf                    , 
                                         num_samples           = ns                    , 
                                         num_boost_round       = num_boost_round       ,
                                         early_stopping_rounds = early_stopping_rounds )
            
        import lightgbm as LightGBM
        callbacks = []
        if 0 < early_stopping_rounds :
            callbacks.append ( LightGBM.early_stopping ( stopping_rounds = early_stopping_rounds , verbose = False ) )

        train_data = LightGBM.Dataset ( X_train , label = Y_train_mod , weight = W_train_mod , free_raw_data = False , params    = params     )
        val_data   = LightGBM.Dataset ( X_val   , label = Y_val_mod   , weight = W_val_mod   , free_raw_data = False , reference = train_data )

        model = LightGBM.train ( params          = params          ,
                                 train_set       = train_data      ,
                                 num_boost_round = num_boost_round ,
                                 valid_sets      = [ val_data ]    ,
                                 callbacks       = callbacks       )

        predictions = model.predict ( X_val , num_iteration = model.best_iteration )
        imps = model.feature_importance ( importance_type = 'gain') if importance else None
        
        return predictions, imps

# =======================================================================================
## @class ADVAL_XGB (Regression)
class ADVAL_XGB (ADVAL_base) : 
    def __init__ ( self             ,
                   nToys    = 400   , **params ) :

        config = {  'objective'        : 'reg:squarederror',
                    'eval_metric'      : 'rmse',
                    'tree_method'      : 'hist',
                    'learning_rate'    : 0.03  ,
                    'n_estimators'     : DEFAULT_ESTIMATORS ,
                    'max_depth'        : MAX_DEPTH ,
                    'max_leaves'       : 15    ,
                    'min_child_weight' :  1.0  ,
                    'subsample'        :  0.8  ,
                    'colsample_bytree' :  0.8  ,
                    'alpha'            :  0.1  ,
                    'lambda'           :  1.0  ,
                    'n_jobs'           : -1    ,
                    'verbosity'        :  0
                }
        config.update ( params ) 
        
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              normalize = False    ,
                              method    = method_XGB , **config   ) 

        if 'n_jobs' in self.params :
            self.params [ 'nthread' ] = self.params.pop ( 'n_jobs' , 1 )
            
        if self.silent : self.params [ 'verbosity' ] = 0
            
            
    # =========================================================================
    ## Parameters for strong regularization
    #  @code
    #  params = gof.regularization ( params , n_features , n_samples ) 
    #  @endcode
    def regularization ( self       ,
                         params     ,
                         n_features , 
                         n_samples  ) :
        
        """ Parameters for strong regularization
        >>> params = gof.regularization ( params , n_features , n_samples ) 
        """

        max_depth = min ( 6 , max ( 1 , n_features ) )
        params [ 'max_depth'         ] = max_depth
        
        params [ 'min_child_weight'  ] = 1e-3
        
        params [ 'colsample_bytree'  ] = 1.0
        params [ 'subsample'         ] = 1.0
        params [ 'alpha'             ] = 0.0
        params [ 'lambda'            ] = 0.0
        params [ 'learning_rate'     ] = min ( 0.05 , params.get ( 'learning_rate' , 0.05 ) )
        
        return params
    
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        
        import xgboost as XGBoost             

        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        Y_val_mod,   W_val_mod   = transform_weights_and_targets ( Y_val  , W_val   )

        dtrain = XGBoost.DMatrix ( X_train , label = Y_train_mod , weight = W_train_mod )
        dval   = XGBoost.DMatrix ( X_val   , label = Y_val_mod   , weight = W_val_mod   )
        evals  = [ ( dtrain , 'train' ) , ( dval , 'val' ) ]

        params = {}
        
        params.update ( self.params )

        num_boost_round       = params.pop ( 'num_boost_round' , None ) or params.pop ( 'n_estimators' , None ) or DEFAULT_ESTIMATORS 
        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , 10 )
        
        nf = num_features ( X_train )
        ns = num_samples  ( X_train )
        
        reg_case = BDT_needs_regularization ( X_train , W_train )        
        if reg_case :
            self.regularized [ reg_case ] += 1
            
            ## update parameters 
            params.update ( self.regularization ( params , nf , ns ) ) 

            num_boost_round        = min ( 20 if 1 == nf else 50 , num_boost_round )
            early_stopping_rounds  = max ( 10 , early_stopping_rounds )

            ## print regularized parameters in "no-silent" regime
            self.report_regularization ( params                ,
                                         reg_case              = reg_case              ,                                          
                                         num_features          = nf                    , 
                                         num_samples           = ns                    , 
                                         num_boost_round       = num_boost_round       ,
                                         early_stopping_rounds = early_stopping_rounds )
            
            
        model = XGBoost.train ( params                = params,
                                dtrain                = dtrain,
                                num_boost_round       = num_boost_round,
                                evals                 = evals,
                                early_stopping_rounds = early_stopping_rounds,
                                verbose_eval          = False )    

        predictions = model.predict ( dval )
        
        if importance :
            score_dict = model.get_score ( importance_type = 'gain' )
            if hasattr ( X_train , 'columns' ) :
                feature_names = X_train.columns.tolist()
            else :
                feature_names = [ f'f{i}' for i in range ( num_features ( X_train ) ) ]
            imps = numpy.array ( [ score_dict.get ( col , 0.0 ) for col in feature_names ] , dtype = numpy.float32 )
        else :
            imps = None

        return predictions, imps

# ======================================================================================
## @class ADVAL_CATB (Regression)
class ADVAL_CATB (ADVAL_base) : 
    def __init__ ( self             ,
                   nToys    = 400   , **params ) :
        
        config = {  'loss_function'         : 'RMSE',
                    'eval_metric'           : 'RMSE',
                    'learning_rate'         : 0.03,
                    'n_estimators'          : DEFAULT_ESTIMATORS ,
                    'early_stopping_rounds' :  20   ,
                    'depth'                 :  MAX_DEPTH ,
                    'min_data_in_leaf'      :  20   ,
                    'l2_leaf_reg'           :   3.0 ,
                    'subsample'             :   0.8 ,
                    'bootstrap_type'        : 'Bernoulli',
                    'verbose'               : False
            }

        
        if   'random_seed'  in params : pass
        elif 'random_state' in params : params [ 'random_seed' ] = params.pop ( 'random_state' , None  )

        config.update ( params )
        import catboost as CatBoost 
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              normalize = False    ,
                              method    = method_CATB , **config   ) 

        ## 
        if 'n_jobs' in self.params :
            self.params [ 'thread_count' ] = self.params.pop ( 'n_jobs' , 1 )

    # =========================================================================
    ## Parameters for strong regularization
    #  @code
    #  params = gof.regularization ( params , n_features , n_samples ) 
    #  @endcode
    def regularization ( self       ,
                         params     ,
                         n_features , 
                         n_samples  ) :
        
        """ Parameters for strong regularization
        >>> params = gof.regularization ( params , n_features , n_samples ) 
        """
    
        depth = min ( 6 , max ( 1 , n_features ) )
        params [ 'depth'             ] = depth
        
        # 
        params [ 'min_data_in_leaf'  ] = max ( 20 , int ( n_samples * 0.001 ) )
        
        params [ 'rsm'               ] = 1.0
        params [ 'subsample'         ] = 1.0
        params [ 'l2_leaf_reg'       ] = 1.0
        params [ 'learning_rate'     ] = min ( 0.05 , params.get ( 'learning_rate' , 0.05 ) )
        
        return params
    
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        
        import catboost as CatBoost
        
        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        Y_val_mod,   W_val_mod   = transform_weights_and_targets ( Y_val  , W_val   )

        train_pool = CatBoost.Pool ( X_train , label = Y_train_mod , weight = W_train_mod )
        val_pool   = CatBoost.Pool ( X_val   , label = Y_val_mod   , weight = W_val_mod   )

        params = {}
        params.update ( self.params )
        
        iterations            = params.pop ( 'iterations' , None ) or params.pop ( 'n_estimators' , None ) or DEFAULT_ESTIMATORS 
        early_stopping_rounds = params.pop ( 'early_stopping_rounds' , 10 )
        
        nf =  num_features ( X_train )
        ns =  num_samples  ( X_train )

        reg_case = BDT_needs_regularization ( X_train , W_train )
        if reg_case : 

            self.regularized [ reg_case ] += 1
            
            ## update parameters 
            params.update ( self.regularization ( params , nf , ns ) )
            
            iterations             = min ( 20 if 1 == nf else 50 , iterations )
            
            # --- CRITICAL FIX FOR CATBOOST: MUST BE NONE (OR > 0), NEVER 0 ---
            early_stopping_rounds  = None

            ## print regularized parameters in "no-silent" regime
            self.report_regularization ( params                ,
                                         reg_case              = reg_case              ,
                                         num_features          = nf                    , 
                                         num_samples           = ns                    , 
                                         iterations            = iterations            ,
                                         early_stopping_rounds = early_stopping_rounds )

        params [ 'iterations'            ] = iterations
        params [ 'early_stopping_rounds' ] = early_stopping_rounds  
        params [ 'use_best_model'        ] = True if early_stopping_rounds else False
            
        model = CatBoost.CatBoostRegressor ( **params ) 
        model.fit ( train_pool , eval_set = val_pool )

        predictions = model.predict ( val_pool )
        imps = model.get_feature_importance ( data = train_pool ) if importance else None

        return predictions, imps

# =============================================================================
## @class ADVAL_HGBC (Regression)
class ADVAL_HGBC (ADVAL_base) : 
    def __init__ ( self             ,
                   nToys    = 400   , **params ) :
        
        config = {  'loss'              : 'squared_error',
                    'learning_rate'     : 0.03   ,
                    'max_iter'          : DEFAULT_ESTIMATORS ,
                    'max_depth'         : MAX_DEPTH ,
                    'max_leaf_nodes'    : 31    ,
                    'min_samples_leaf'  : 20    ,
                    'l2_regularization' :  0.1  ,
                    'early_stopping'    :  True ,
                    'n_iter_no_change'  : 20    ,
                    'scoring'           : 'neg_root_mean_squared_error',
                }
        
        config.update ( params )
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              normalize = False    ,
                              method    = method_HGBC , **config   ) 

        if 'n_jobs' in self.params : self.params.pop ( 'n_jobs' , None )
        

    # =========================================================================
    ## Parameters for strong regularization
    #  @code
    #  params = gof.regularization ( params , n_features , n_samples ) 
    #  @endcode
    def regularization ( self       ,
                         params     ,
                         n_features , 
                         n_samples  ) :
        
        """ Parameters for strong regularization
        >>> params = gof.regularization ( params , n_features , n_samples ) 
        """
        nf = n_features
        
        max_depth      = 1 if 1 == nf else min ( 2 if nf <= 3 else 3 , params.get ( 'max_depth' , 5 ) )
        max_leaf_nodes = 2 if max_depth == 1 else min ( 2 ** max_depth , params.get ( 'max_leaf_nodes' , 31 ) )
        
        params [ 'max_depth'            ] = max_depth
        params [ 'max_leaf_nodes'       ] = max_leaf_nodes
        
        params [ 'min_samples_leaf'     ] = max ( 50 , params.get ( 'min_samples_leaf' , 20 ) )
        
        # --- Moderate L2 regularization to suppress sPlot noise in binning ---
        params [ 'l2_regularization'    ] = 1.0
        
        params [ 'learning_rate'        ] = min ( 0.05 , params.get ( 'learning_rate' , 0.05 ) )

        return params


    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        
        from sklearn.ensemble import HistGradientBoostingRegressor

        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        w_tr_arr = W_train_mod.values if hasattr ( W_train_mod , 'values' ) else W_train_mod

        nf =  num_features ( X_train )
        ns =  num_samples  ( X_train )
        
        params = {}
        params.update ( self.params )
        max_iter = params.pop ( 'max_iter' , None ) or params.pop ( 'n_estimators' , None ) or DEFAULT_ESTIMATORS 
        
        reg_case = BDT_needs_regularization ( X_train , W_train )
        if reg_case : 

            self.regularized [ reg_case ] += 1
            
            ## update parameters 
            params.update ( self.regularization ( params , nf , ns ) )

            max_iter = min ( 20 if 1 == nf else 50 , max_iter )
            
            ## print regularized parameters in "no-silent" regime
            self.report_regularization ( params                  ,
                                         reg_case     = reg_case , 
                                         num_features = nf       , 
                                         num_samples  = ns       , 
                                         max_iter     = max_iter )
            
        params [ 'max_iter' ] = max_iter

        model = HistGradientBoostingRegressor ( **params )
        model.fit ( X_train , Y_train_mod , sample_weight = w_tr_arr )

        predictions = model.predict ( X_val )
        return predictions, None 
            
# =============================================================================
## @class ADVAL_GBC (Regression)
class ADVAL_GBC (ADVAL_base) : 
    def __init__ ( self             ,
                   nToys    = 400   , **params   ) :
        
        config = {  'loss'              : 'squared_error',
                    'learning_rate'     : 0.05  ,
                    'n_estimators'      : DEFAULT_ESTIMATORS ,
                    'max_depth'         : MAX_DEPTH ,
                    'min_samples_split' :  10   ,
                    'min_samples_leaf'  :   5   ,
                    'subsample'         :   1.0 ,
                    'max_features'      :   1.0 ,
                }
        
        config.update ( params ) 
        
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              normalize = False    ,
                              method    = method_GBC , **config   ) 
        
        if 'n_jobs' in self.params : self.params.pop ( 'n_jobs' , None )
        
    # =========================================================================
    ## Parameters for strong regularization
    #  @code
    #  params = gof.regularization ( params , n_features , n_samples ) 
    #  @endcode
    def regularization ( self       ,
                         params     ,
                         n_features , 
                         n_samples  ) :
        
        """ Parameters for strong regularization
        >>> params = gof.regularization ( params , n_features , n_samples ) 
        """
        nf = n_features
        
        max_depth = 1 if 1 == nf else min ( 2 , params.get ( 'max_depth' , 5 ) )
        
        params [ 'max_depth'                ] = max_depth
        
        params [ 'min_samples_leaf'         ] = 5
        params [ 'min_samples_split'        ] = 10
        
        params [ 'subsample'                ] = 1.0
        params [ 'ccp_alpha'                ] = 0.0
        
        params [ 'learning_rate'            ] = min ( 0.05 , params.get ( 'learning_rate' , 0.05 ) )

        return params 

    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        from sklearn.ensemble import GradientBoostingRegressor

        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        w_tr_arr = W_train_mod.values if hasattr(W_train_mod, 'values') else W_train_mod

        nf =  num_features ( X_train )
        ns =  num_samples  ( X_train )
        
        params = {}
        params.update ( self.params )
        n_estimators = params.pop ( 'num_boost_round' , None ) or params.pop ( 'n_estimators' , None ) or DEFAULT_ESTIMATORS 
        
        reg_case = BDT_needs_regularization ( X_train , W_train )
        if reg_case : 

            self.regularized [ reg_case ] += 1

            ## update parameters 
            params.update ( self.regularization ( params , nf , ns ) )

            n_estimators = min ( 20 if 1 == nf else 50 , n_estimators , DEFAULT_ESTIMATORS )
            
            ## print regularized parameters in "no-silent" regime
            self.report_regularization ( params                      ,
                                         reg_case     = reg_case     ,  
                                         num_features = nf           , 
                                         num_samples  = ns           ,
                                         n_estimators = n_estimators )
                                                     
        params [ 'n_estimators' ] = n_estimators
        
        model = GradientBoostingRegressor ( **params )
        model.fit ( X_train , Y_train_mod , sample_weight = w_tr_arr )

        predictions = model.predict ( X_val )
        imps = model.feature_importances_ if importance else None

        return predictions, imps

# =============================================================================
## @class ADVAL_RF (Regression)
class ADVAL_RF (ADVAL_base) : 
    def __init__ ( self           ,
                   nToys    = 400 , **params   ) :
  
        config = {  'n_estimators'      : DEFAULT_ESTIMATORS ,
                    'n_jobs'            : -1    ,
                    'criterion'         : 'squared_error',
                    'max_depth'         : MAX_DEPTH ,
                    'min_samples_split' : 10    ,
                    'min_samples_leaf'  :  5    ,
                    'max_features'      :  1.0  ,
                    'bootstrap'         :  True ,
                    'max_samples'       :  0.8  ,
                }  
        config.update ( params ) 
        
        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              normalize = False    ,                              
                              method    = method_RF, **config   ) 
        
    # =========================================================================
    ## Parameters for strong regularization
    #  @code
    #  params = gof.regularization ( params , n_features , n_samples ) 
    #  @endcode
    def regularization ( self       ,
                         params     ,
                         n_features , 
                         n_samples  ) :
        
        """ Parameters for strong regularization
        >>> params = gof.regularization ( params , n_features , n_samples ) 
        """
        nf = n_features 
        ns = n_samples 
        
        max_depth = 1 if 1 == nf else min ( 2 if nf <= 3 else 3 , params.get ( 'max_depth' , 5 ) )
        
        params [ 'max_depth'          ] = max_depth
        
        # --- Dynamic minimum samples leaf scaling with increased Asimov dataset size ---
        params [ 'min_samples_leaf'   ] = max ( int ( 0.01 * ns ) , 100 )
        params [ 'min_samples_split'  ] = max ( int ( 0.02 * ns ) , 200 )
        
        # --- Use full bootstrap samples to eliminate mean prediction bias on large Asimov datasets ---
        params [ 'max_samples'        ] = None
        params [ 'bootstrap'          ] = True
        
        # --- Disable CCP pruning; sPlot negative weights distort complexity-cost pruning ---
        params [ 'ccp_alpha'          ] = 0.0

        return params 
         
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        from sklearn.ensemble import RandomForestRegressor

        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        w_tr_arr = W_train_mod.values if hasattr(W_train_mod, 'values') else W_train_mod

        params = {}
        params.update ( self.params )
        n_estimators = params.pop ( 'n_estimators' , None ) or params.pop ( 'num_boost_round' , None ) or DEFAULT_ESTIMATORS 

        nf =  num_features ( X_train )
        ns =  num_samples  ( X_train )
        
        reg_case = BDT_needs_regularization ( X_train , W_train )
        if reg_case : 

            self.regularized [ reg_case ] += 1

            
            ## update parameters 
            params.update ( self.regularization ( params , nf , ns ) )
            
            n_estimators = min ( 20 if 1 == nf else 50 , n_estimators , MAX_REGULARIZED_ESTIMATORS ) 
                                 
            ## print regularized parameters in "no-silent" regime
            self.report_regularization ( params                      ,
                                         reg_case     = reg_case     ,  
                                         num_features = nf           , 
                                         num_samples  = ns           ,
                                         n_estimators = n_estimators )

            
        params [ 'n_estimators' ] = n_estimators
        
        model = RandomForestRegressor ( **params )
        model.fit ( X_train , Y_train_mod , sample_weight = w_tr_arr )

        predictions = model.predict ( X_val )
        imps = model.feature_importances_ if importance else None

        return predictions, imps

# =============================================================================
## @class ADVAL_TORCH (Regression)
class ADVAL_TORCH (ADVAL_base) : 
    def __init__ ( self             ,
                   nToys    = 100   , **params   ) :

        config =  { 'epochs'                : 100   ,
                    'batch_size'            : 256   ,
                    'learning_rate'         : 1.e-3 ,
                    'weight_decay'          : 1.e-4 ,
                    'early_stopping_rounds' :  15   , 
                    'dropout'               : 0.1   }

        config.update ( params ) 

        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              normalize = True     , 
                              method    = method_TORCH, **config   ) 

    # =========================================================================
    ## Parameters for strong regularization
    def regularization ( self       ,
                         params     ,
                         n_features , 
                         n_samples  ) :
        """ Parameters for strong regularization
        """
        nf = n_features
        ns = n_samples
        
        params [ 'batch_size'            ] = max ( int ( 0.02 * ns ) , 256 )
        params [ 'weight_decay'          ] = 1e-3
        params [ 'hidden_dim'            ] = 8 if 1 == nf else 16
        params [ 'learning_rate'         ] = min ( 0.005 , params.get ( 'learning_rate' , params.get ( 'lr' , 1e-3 ) ) )
        params [ 'early_stopping_rounds' ] = 0
        
        return params 
    
    def work ( self ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        import copy
        import torch    as Torch 
        import torch.nn as NN 
        from   torch.utils.data import DataLoader, TensorDataset

        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        Y_val_mod,   W_val_mod   = transform_weights_and_targets ( Y_val  , W_val   )

        params = {}
        params.update ( self.params )
        epochs = params.pop ( 'epochs' , None ) or 100

        nf = num_features ( X_train )
        ns = num_samples  ( X_train )
        
        reg_case = NN_needs_regularization ( X_train , W_train )
        if reg_case : 
            self.regularized [ reg_case ] += 1
            
            params.update ( self.regularization ( params , nf , ns ) )
            epochs = min ( 30 if 1 == nf else 50 , epochs ) 
            self.report_regularization ( params       ,
                                         reg_case     = reg_case ,         
                                         num_features = nf       ,
                                         num_samples  = ns       ,
                                         epochs       = epochs   )
            
        batch_size   = params.get ( 'batch_size'            , 256  )
        lr           = params.get ( 'learning_rate'         , params.get ( 'lr' , 1e-3 ) )
        weight_decay = params.get ( 'weight_decay'          , 1e-4 )
        patience     = params.get ( 'early_stopping_rounds' ,  15  )        
        dropout_rate = params.get ( 'dropout'               , 0.1  )
        hidden_dim   = params.get ( 'hidden_dim'            , 64 if 15 < nf else 16 )

        n_jobs = params.get ( 'n_jobs', -1 )
        if 0 < n_jobs : Torch.set_num_threads ( n_jobs )
        else          : Torch.set_num_threads ( 1      )
            
        device = Torch.device ( 'cuda' if Torch.cuda.is_available() else 'cpu' )

        X_tr_arr = numpy.nan_to_num ( numpy.asarray ( X_train , dtype = numpy.float32 ), nan = 0.0 )
        Y_tr_arr = numpy.asarray ( Y_train_mod , dtype = numpy.float32 ) . reshape ( -1 , 1 )        
        W_tr_arr = ( numpy.ones_like ( Y_tr_arr ) if W_train_mod is None 
                     else numpy.asarray ( W_train_mod , dtype = numpy.float32 ) . reshape  (-1 , 1 ) )
        
        X_va_arr = numpy.nan_to_num ( numpy.asarray ( X_val , dtype = numpy.float32 ), nan = 0.0 )
        Y_va_arr = numpy.asarray ( Y_val_mod  , dtype = numpy.float32 ) . reshape ( -1 , 1 )
        W_va_arr = ( numpy.ones_like ( Y_va_arr ) if W_val_mod is None 
                     else numpy.asarray ( W_val_mod , dtype = numpy.float32 ) . reshape (-1 , 1 ) )

        train_dataset = TensorDataset ( Torch.from_numpy ( X_tr_arr ) , 
                                        Torch.from_numpy ( Y_tr_arr ) , 
                                        Torch.from_numpy ( W_tr_arr ) )
        train_loader = DataLoader ( train_dataset , batch_size = batch_size , shuffle = True )

        val_x_tensor = Torch.from_numpy ( X_va_arr ).to ( device )
        val_y_tensor = Torch.from_numpy ( Y_va_arr ).to ( device )
        val_w_tensor = Torch.from_numpy ( W_va_arr ).to ( device )

        model = NN.Sequential ( NN.Linear ( nf , hidden_dim ) ,
                                NN.ReLU   ()                   ,
                                NN.Dropout( dropout_rate )     ,
                                NN.Linear ( hidden_dim , 1 )   ).to ( device )
            
        optimizer = Torch.optim.AdamW ( model.parameters() , lr = lr , weight_decay = weight_decay )
        
        def compute_weighted_mse ( preds , targets , weights ):
            loss = NN.functional.mse_loss ( preds , targets , reduction = 'none' )
            return ( loss * weights ).mean()

        best_val_loss      = float('inf')
        best_model_weights = None
        patience_counter   = 0

        for epoch in range ( epochs ) :
            model.train()
            for bx, by, bw in train_loader:
                bx, by, bw = bx.to ( device ) , by.to ( device ) , bw.to ( device )
                optimizer.zero_grad()
                preds = model(bx)
                loss   = compute_weighted_mse(preds, by, bw)
                loss.backward()
                optimizer.step()

            model.eval()
            with Torch.no_grad():
                val_preds = model ( val_x_tensor )
                val_loss   = compute_weighted_mse ( val_preds , val_y_tensor, val_w_tensor).item()

            if val_loss < best_val_loss:
                best_val_loss      = val_loss
                best_model_weights = copy.deepcopy ( model.state_dict() )
                patience_counter   = 0
            else:
                patience_counter += 1
                if patience and patience_counter >= patience:
                    break

        if best_model_weights is not None:
            model.load_state_dict ( best_model_weights )

        model.eval()
        with Torch.no_grad():
            predictions = model ( val_x_tensor ) .cpu().numpy().ravel()

        return predictions, None


# =============================================================================
## @class ADVAL_KERAS (Regression)
class ADVAL_KERAS (ADVAL_base) : 
    def __init__ ( self             ,
                   nToys    = 100   , **params   ) :

        config =  { 'epochs'                :  50    ,
                    'batch_size'            : 256    ,
                    'learning_rate'         :   0.01 ,
                    'weight_decay'          :   1e-4 , 
                    'early_stopping_rounds' :  15    , 
                    'dropout'               :   0.1  }

        ADVAL_base.__init__ ( self, 
                              nToys     = nToys    ,
                              normalize = True     , 
                              method    = method_KERAS, **config   ) 
        
    # =========================================================================
    ## Parameters for strong regularization
    def regularization ( self       ,
                         params     ,
                         n_features , 
                         n_samples  ) :
        """ Parameters for strong regularization
        """
        nf = n_features
        ns = n_samples

        params [ 'batch_size'            ] = max ( int ( 0.02 * ns ) , 256 )                
        params [ 'weight_decay'          ] = 1e-3
        params [ 'hidden_dim'            ] = 8 if 1 == nf else 16
        params [ 'learning_rate'         ] = min ( 0.005 , params.get ( 'learning_rate' , 0.01 ) ) 
        params [ 'early_stopping_rounds' ] = 0

        return params 
        
    def work ( self    ,
               X_train , Y_train , W_train ,
               X_val   , Y_val   , W_val   , importance = False ) :
        
        import keras as     Keras 
        from   keras import layers as Layers 

        Y_train_mod, W_train_mod = transform_weights_and_targets ( Y_train, W_train )
        Y_val_mod,   W_val_mod   = transform_weights_and_targets ( Y_val  , W_val   )

        params = {}
        params.update ( self.params )
        
        nf     = num_features ( X_train )
        ns     = num_samples  ( X_train )
        epochs = params.pop ( 'epochs' , None ) or 100
        
        reg_case = NN_needs_regularization ( X_train , W_train )
        if reg_case : 
            self.regularized [ reg_case ] += 1
            params.update ( self.regularization ( params , nf , ns ) )
            epochs = min ( 30 if 1 == nf else 50 , epochs )
            self.report_regularization ( params       ,
                                         reg_case     = reg_case ,                                                  
                                         num_features = nf       ,
                                         num_samples  = ns       ,
                                         epochs       = epochs   )

        n_features   = nf 
        hidden_dim   = params.get ( 'hidden_dim'           , 64 if 15 < n_features else 16 )
        batch_size   = params.get ( 'batch_size'           , max ( int ( 0.02 * ns ) , 256 ) ) 

        patience     = params.get ( 'early_stopping_rounds',  15   )
        dropout_rate = params.get ( 'dropout'              ,   0.1 )
        lr           = params.get ( 'learning_rate'        , params.get ( 'lr' , 0.01 ) )
        weight_decay = params.get ( 'weight_decay'         , 1e-4  )
        
        X_tr_arr = numpy.nan_to_num ( numpy.asarray ( X_train , dtype = numpy.float32 ).reshape ( -1 , n_features ) , nan = 0.0 )
        Y_tr_arr = numpy.asarray ( Y_train_mod , dtype = numpy.float32 ).reshape ( -1 , 1 )
        W_tr_arr = ( numpy.ones_like ( Y_tr_arr ) if W_train_mod is None 
                     else numpy.asarray ( W_train_mod , dtype = numpy.float32 ).reshape ( -1 , 1 ) )

        X_va_arr = numpy.nan_to_num ( numpy.asarray ( X_val , dtype = numpy.float32 ).reshape ( -1 , n_features ) , nan = 0.0 )
        Y_va_arr = numpy.asarray ( Y_val_mod , dtype = numpy.float32 ).reshape ( -1 , 1 )
        W_va_arr = ( numpy.ones_like ( Y_va_arr ) if W_val_mod is None 
                     else numpy.asarray ( W_val_mod , dtype = numpy.float32 ).reshape ( -1 , 1 ) )
        
        inputs = Keras.Input ( shape = ( n_features , ) )
        x      = Layers.Dense ( hidden_dim , activation = 'relu' )( inputs )
        if 0 < dropout_rate:
            x = Layers.Dropout ( dropout_rate )( x )
        outputs = Layers.Dense ( 1 )( x )

        model = Keras.Model ( inputs = inputs , outputs = outputs )

        optimizer = Keras.optimizers.AdamW ( learning_rate = lr , weight_decay = weight_decay )
        model.compile ( optimizer        = optimizer ,
                        loss             = Keras.losses.MeanSquaredError() ,
                        weighted_metrics = [] )

        verbose = 0 if self.silent else 1
        
        callbacks = []
        if patience and 0 < patience:
            callbacks.append (
                Keras.callbacks.EarlyStopping ( monitor              = 'val_loss' ,
                                                 patience             = patience  ,
                                                 restore_best_weights = True      ,
                                                 verbose              = verbose   )
            )
        
        model.fit ( X_tr_arr        ,
                    Y_tr_arr        ,
                    sample_weight   = W_tr_arr.ravel() ,
                    validation_data = ( X_va_arr , Y_va_arr , W_va_arr.ravel() ) ,
                    epochs          = epochs             ,
                    batch_size      = batch_size         ,
                    callbacks       = callbacks          ,
                    verbose         = verbose            )

        predictions = model.predict ( X_va_arr , batch_size = batch_size , verbose = 0 ).ravel()

        Keras.backend.clear_session()
        
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
    
