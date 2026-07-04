#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/adversarial_vaildation.py
#  Implement adversarial vaildation to probe
#  the difference between two weighted datasets
#
#  @see ADVAL_LGBM, class based on lightgbm, the most CPU efficient 
#  @see ADVAL_HGBC, class based on HistGradientBoosterClassifier, slower 
#  @see ADVAL_GBC
#
#  As t-value \f$ 400 \times \left( \frac{1}{2} - AUC \right)^2\f$ is used
#  To estimate the~p-value permutations are used.
# 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-04
# =============================================================================
""" Implement adversarial vaildation to probe
    the difference between two weighted datasets

 - see ADVAL_LGBM, class based on lightgbm, the most CPU efficient 
 - see ADVAL_HGBC, class based on HistGradientBoosterClassifier, slower 
 - see ADVAL_GBC

As t-value 400*( AUC - 0.5)**2  is used
To estimate the~p-value permutations are used. 

"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2026-07-04"
__all__     = (
    'ADVAL_LGBM' , ## LightGBM-based class for Adversarial Validation for the differecne between two (weighted) dataset
    'ADVAL_HGBC' , ## HGBC-based class for Adversarial Validation for the difference between two (weighted) dataset
    'ADVAL_GBC'  , ## GBC-based class for Adversarial Validation for the difference between two (weighted) dataset
)
# =============================================================================
from   ostap.core.ostap_types   import string_types
from   ostap.stats.gof_np       import GoFnp , s2u 
from   ostap.stats.gof_utils    import PERMUTATOR, normalize as ds_normalize
import ROOT, numpy 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.adversarial_validation' )
else                       : logger = getLogger ( __name__ )
# =============================================================================
logger.debug ( 'Implement adversarial vaildation for Goodness-of-fit and Two-Sample discrimination ' )
# =============================================================================
## @class ADVAL_base 
#  Base class for adversarial validation for the difference between two (weighted) dataset
class ADVAL_base (GoFnp):
    """ Base class for Adversarial validation for the differece between two (weighted) dataset
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  ,
                   method   = "Adversarial Validation" ,
                   n_splits = 8     , ##  number fo splits for cross-validation 
                   **params         ) :
        
        GoFnp.__init__ ( self                ,
                         nToys    = nToys    ,
                         parallel = parallel , 
                         silent   = silent   ,
                         progress = progress , 
                         method   = method   )
        
        ## update parameters from argument 
        self.__params = {} 
        self.__params.update ( params ) 
        ## Empirical CDF for t-value distribution from permutations"""
        self.__ecdf    = None

        assert isinstance ( n_splits , int ) and 2 <= n_splits , "Invalid n_splits:%s" % n_splits
        self.__n_splits = n_splits 
        
    # =========================================================================    
    ## Calculate T-value for two (structured) datasets 
    #  @code
    #  adval  = ...
    #  data1 = ... ## the first  data set 
    #  data2 = ... ## the second data set
    #  t = adval ( data1 , data1 , normalize = False ) 
    #  t = adval ( data1 , data1 , normalize = True  ) 
    #  @endcode
    def __call__ ( self              ,
                   data1             ,
                   data2             , * ,
                   weight1   = None  ,
                   weight2   = None  ,
                   normalize = False ) :
        
        """ Calculate T-value for two (STRUCTURED) data sets 
        >>> adval = ...
        >>> data1 = ... ## the first  data set 
        >>> data2 = ... ## the second data set
        >>> t = adval ( data1 , data1 , normalize = False ) 
        >>> t = adval ( data1 , data1 , normalize = True  ) 
        """

        ## transform/normalize ? 
        if normalize : ds1 , ds2 = self.normalize ( data1 , data2 )
        else         : ds1 , ds2 = data1 , data2 

        ## convert to unstructured datasets 
        uds1 = s2u ( ds1 , copy = False ) if ds1.dtype.fields else ds1
        uds2 = s2u ( ds2 , copy = False ) if ds2.dtype.fields else ds2

        return self.t_value ( uds1 , uds2 , weight1 = weight1 , weight2 = weight2 )
    
    # =========================================================================
    ## Calculate the t & p-values
    #  @code
    #  adval = ...
    #  ds1 , ds2 = ...
    #  t , p = adval.pvalue ( ds1 , ds2 , normalize = True ) 
    #  @endcode 
    def pvalue ( self              ,
                 data1             ,
                 data2             , * ,
                 weight1   = None  ,
                 weight2   = None  ,
                 normalize = False ) :
        """ Calculate the t & p-values
        >>> adval = ...
        >>> ds1 , ds2 = ...
        >>> t , p = adval.pvalue ( ds1 , ds2 , normalize = True ) 
        """
        
        ## transform/normalize ? 
        if normalize : ds1 , ds2 = self.normalize ( data1 , data2 )
        else         : ds1 , ds2 = data1 , data2 

        ## convert to unstructured datasets 
        uds1 = s2u ( ds1 , copy = False ) if ds1.dtype.fields else ds1
        uds2 = s2u ( ds2 , copy = False ) if ds2.dtype.fields else ds2

        t_value    = self.t_value ( uds1 , uds2 , weight1 = weight1 , weight2 = weight2 )

        permutator = PERMUTATOR ( self , t_value , uds1 , uds2 , weight1 = weight1 , weight2 = weight2 )

        if self.parallel and permutator.run :
            counter = permutator.run ( self.nToys , progress = self.progress )            
        else :
            counter = permutator     ( self.nToys , progress = self.progress )

        self.__ecdf = permutator.ecdf
        
        p_value = counter.eff

        ## if self.__increasing : p_value = 1 - p_value

        return t_value , p_value
    
    @property
    def ecdf ( self ) :
        """`ecdf` : empirical CDF for t-value distribution from permutations"""
        return self.__ecdf
    
    @property
    def params ( self ) :
        """`params` : configuration of underlying classifier"""
        return self.__params
    
    @property
    def n_splits ( self ) :
        """`n_splits`: # of splits for cross-validation"""
        return self.__n_splits 

# =======================================================================================
## @class ADVAL_LGBM
#  ligрtgbm-based lass for adversarial validation for the difference between two (weighted) dataset
#  @see lightgbm
class ADVAL_LGBM (ADVAL_base) : 
    """ LightGBM-based class for Adversarial validation for the differece between two (weighted) dataset
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  ,
                   n_splits = 8     , **params ) :
        
        config = {
            'objective'     : 'binary' ,
            'metric'        : 'auc'    ,
            'learning_rate' :  0.05    , 
            'max_depth'     :  5       ,
            'verbose'       : -1       ,
            'random_state'  : 42
        }
        config.update ( params ) 
        
        ADVAL_base.__init__ ( self, 
                              nToys    = nToys    ,
                              parallel = parallel ,
                              silent   = silent   , 
                              progress = progress , 
                              method   = "Adversarial Validation/LightGBM" ,
                              n_splits = n_splits , **config   ) 
        
    # ==========================================================================
    ## Calculate t-value for two non-structured (weighted) datasets 
    #   @param data1   the first dataset
    #   @param data2   the second dataset
    #   @param weight1 the first array of weights 
    #   @param weight3 the second array of weights
    #   tvalue is defined as \f$  400 \times \left( \frac{1}{2} - AUC \right)^2 \f$ 
    def t_value ( self  ,
                  data1          ,
                  data2          , * , 
                  weight1 = None ,
                  weight2 = None ) :
        """ Calculate t-value for two non-structured (weighted) arrays 
        - data1   : the first dataset
        - data2   : the second dataset
        - weight1 : the first array of weights 
        - weight3 : the second array of weights 
         t-value is defined as 400 * abs(0.5-AUC)**2
        """

        sh1 = data1.shape 
        sh2 = data2.shape
        
        assert 2 == len ( sh1 ) and 2 == len ( sh2 ) and sh1 [ 1 ] == sh2 [ 1 ] and sh1 [ 0 ] and sh2 [ 0 ] , \
            "Invalid arrays: %s , %s" % ( sh1 , sh2 )
        
        if weight1 is None : weigth1 = numpy.ones ( sh1 [ 0 ] )
        if weight2 is None : weigth2 = numpy.ones ( sh2 [ 0 ] )
        
        ## convert numpy arrays into pandas dataframes
        import pandas as pd 

        df_1 = pd.DataFrame ( data1 ) 
        df_2 = pd.DataFrame ( data2 )

        column_target = 'column_target'
        column_weight = 'column_weight'
        
        df_1 [ column_target ] = 1
        df_2 [ column_target ] = 0
        
        df_1 [ column_weight ] = weight1 
        df_2 [ column_weight ] = weight2 

        ## merge datasets together
        dataset = pd.concat ( [ df_1 , df_2 ] , axis = 0 ).reset_index ( drop = True )

        X       = dataset . drop ( columns = [ column_target , column_weight ]  )
        Y       = dataset [ column_target ]
        W       = dataset [ column_weight ]
        N       = len ( dataset )
        
        ## cross-validation
        from sklearn.model_selection import StratifiedKFold
        from sklearn.metrics         import roc_auc_score 
        cv        = StratifiedKFold ( n_splits = self.n_splits , shuffle=True , random_state = 42 )
        oof_preds = numpy.zeros( N )

        import lightgbm as LGB 

        ## model parameters
        params = {
            'objective'     : 'binary' ,
            'metric'        : 'auc'    ,
            'learning_rate' :  0.05    , 
            'max_depth'     :  5       ,
            'verbose'       : -1       ,
            'random_state'  : 42
        }
        
        for fold, ( train_idx , val_idx ) in enumerate ( cv.split ( X, Y ) ):
            
            X_train , Y_train , W_train = X.iloc [ train_idx ] , Y.iloc [ train_idx ] , W.iloc [ train_idx ]
            X_val   , Y_val   , W_val   = X.iloc [   val_idx ] , Y.iloc [   val_idx ] , W.iloc [   val_idx ]

            train_data = LGB.Dataset ( X_train , label = Y_train , weight = W_train)
            val_data   = LGB.Dataset (   X_val , label =   Y_val , weight = W_val , reference = train_data )
        
            # train the model
            model = LGB.train (
                self.params ,
                train_data  ,
                num_boost_round = 500 ,
                valid_sets = [ val_data ],
                callbacks  = [ LGB.early_stopping ( stopping_rounds=20 , verbose=False ) ]
            )
            
            # predictions 
            oof_preds [ val_idx ] = model.predict ( X_val , num_iteration = model.best_iteration )
            
        #  weighted ROC-AUC
        auc_score = roc_auc_score ( Y , oof_preds , sample_weight = W )

        ## 
        return  100 * ( 1.0 - 2.0 * auc_score ) ** 2 

# =============================================================================
## @class ADVAL_HGBC
#  Class for adversarial validation for the difference between two (weighted) dataset
#  based HistoGradientBoosterClassifier 
#  @see HistoGradientBoosterClassifier 
class ADVAL_HGBC (ADVAL_base) : 
    """ HGBC-based class for Adversarial validation for the differece between two (weighted) dataset
    @see HistoGradientBoosterClassifier 
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  ,
                   n_splits = 8     , **params ) :
        
        config =  { 'learning_rate'    : 0.05 ,
                    'max_depth'        : 5    ,
                    'max_iter'         : 200  , # Количество деревьев (аналог num_boost_round)
                    'early_stopping'   : True ,  # Включаем встроенный early stopping
                    'n_iter_no_change' : 15   ,
                    'random_state'     : 42   }
        config.update ( params ) 
        
        ADVAL_base.__init__ ( self, 
                              nToys    = nToys    ,
                              parallel = parallel ,
                              silent   = silent   , 
                              progress = progress , 
                              method   = "Adversarial Validation/HGBC" ,
                              n_splits = n_splits  , **config   ) 
        
    # ==========================================================================
    ## Calculate t-value for two non-structured (weighted) datasets 
    #   @param data1   the first dataset
    #   @param data2   the second dataset
    #   @param weight1 the first array of weights 
    #   @param weight3 the second array of weights
    #   tvalue is defined as \f$  400 \times \left( \frac{1}{2} - AUC \right)^2 \f$ 
    def t_value ( self  ,
                  data1          ,
                  data2          , * , 
                  weight1 = None ,
                  weight2 = None ) :
        """ Calculate t-value for two non-structured (weighted) arrays 
        - data1   : the first dataset
        - data2   : the second dataset
        - weight1 : the first array of weights 
        - weight3 : the second array of weights 
         t-value is defined as 400 * abs(0.5-AUC)**2
        """

        sh1 = data1.shape 
        sh2 = data2.shape
        
        assert 2 == len ( sh1 ) and 2 == len ( sh2 ) and sh1 [ 1 ] == sh2 [ 1 ] and sh1 [ 0 ] and sh2 [ 0 ] , \
            "Invalid arrays: %s , %s" % ( sh1 , sh2 )
        
        if weight1 is None : weigth1 = numpy.ones ( sh1 [ 0 ] )
        if weight2 is None : weigth2 = numpy.ones ( sh2 [ 0 ] )
        
        ## convert numpy arrays into pandas dataframes
        import pandas as pd 

        df_1 = pd.DataFrame ( data1 ) 
        df_2 = pd.DataFrame ( data2 )

        column_target = 'column_target'
        column_weight = 'column_weight'
        
        df_1 [ column_target ] = 1
        df_2 [ column_target ] = 0
        
        df_1 [ column_weight ] = weight1 
        df_2 [ column_weight ] = weight2 

        ## merge datasets together
        dataset = pd.concat ( [ df_1 , df_2 ] , axis = 0 ).reset_index ( drop = True )

        n_splits = 5
                
        X       = dataset . drop ( columns = [ column_target , column_weight ]  )
        Y       = dataset [ column_target ]
        W       = dataset [ column_weight ]
        N       = len ( dataset )
        
        ## cross-validation
        from sklearn.model_selection import StratifiedKFold
        from sklearn.metrics         import roc_auc_score 
        cv        = StratifiedKFold ( n_splits = n_splits , shuffle=True , random_state = 42 )
        oof_preds = numpy.zeros( N )
        
        from sklearn.ensemble import HistGradientBoostingClassifier

        ## model parameters
        params = {
            'objective'     : 'binary' ,
            'metric'        : 'auc'    ,
            'learning_rate' :  0.05    , 
            'max_depth'     :  5       ,
            'verbose'       : -1       ,
            'random_state'  : 42
        }
        
        for fold, ( train_idx , val_idx ) in enumerate ( cv.split ( X, Y ) ):
            
            X_train , Y_train , W_train = X.iloc [ train_idx ] , Y.iloc [ train_idx ] , W.iloc [ train_idx ]
            X_val   , Y_val   , W_val   = X.iloc [   val_idx ] , Y.iloc [   val_idx ] , W.iloc [   val_idx ]

            # train the model

            model = HistGradientBoostingClassifier ( **self.params )
            model.fit ( X_train , Y_train , sample_weight = W_train )
            
            # Предсказываем вероятность принадлежности к классу 1 (выборке P)
            oof_preds [ val_idx ] = model.predict_proba ( X_val ) [ : , 1 ]
            
        #  weighted ROC-AUC
        auc_score = roc_auc_score ( Y , oof_preds , sample_weight = W )

        ## 
        return  100 * ( 1.0 - 2.0 * auc_score ) ** 2 


# =============================================================================
## @class ADVAL_GBC
#  Class for adversarial validation for the difference between two (weighted) dataset
#  based GradientBoosterClassifier 
#  @see GradientBoosterClassifier 
class ADVAL_GBC (ADVAL_base) : 
    """ GBC-based class for Adversarial validation for the differece between two (weighted) dataset
    @see GradientBoosterClassifier 
    """
    def __init__ ( self             ,
                   nToys    = 400   ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  ,
                   n_splits = 8     , **params   ) :
        
        config =  { 'learning_rate'       : 0.05 ,
                    'max_depth'           : 5    ,
                    'n_estimators'        : 200  , ## number of trees 
                    'validation_fraction' : 0.10 , ## internal validation for early stopping
                    'n_iter_no_change'    : 15   ,
                    'random_state'        : 42   }
        config.update ( params ) 
        
        ADVAL_base.__init__ ( self, 
                              nToys    = nToys    ,
                              parallel = parallel ,
                              silent   = silent   , 
                              progress = progress , 
                              method   = "Adversarial Validation/GBC" ,
                              n_splits = n_splits , **config   ) 
            
    # ==========================================================================
    ## Calculate t-value for two non-structured (weighted) datasets 
    #   @param data1   the first dataset
    #   @param data2   the second dataset
    #   @param weight1 the first array of weights 
    #   @param weight3 the second array of weights
    #   tvalue is defined as \f$  400 \times \left( \frac{1}{2} - AUC \right)^2 \f$ 
    def t_value ( self  ,
                  data1          ,
                  data2          , * , 
                  weight1 = None ,
                  weight2 = None ) :
        """ Calculate t-value for two non-structured (weighted) arrays 
        - data1   : the first dataset
        - data2   : the second dataset
        - weight1 : the first array of weights 
        - weight3 : the second array of weights 
         t-value is defined as 400 * abs(0.5-AUC)**2
        """

        sh1 = data1.shape 
        sh2 = data2.shape
        
        assert 2 == len ( sh1 ) and 2 == len ( sh2 ) and sh1 [ 1 ] == sh2 [ 1 ] and sh1 [ 0 ] and sh2 [ 0 ] , \
            "Invalid arrays: %s , %s" % ( sh1 , sh2 )
        
        if weight1 is None : weigth1 = numpy.ones ( sh1 [ 0 ] )
        if weight2 is None : weigth2 = numpy.ones ( sh2 [ 0 ] )
        
        ## convert numpy arrays into pandas dataframes
        import pandas as pd 

        df_1 = pd.DataFrame ( data1 ) 
        df_2 = pd.DataFrame ( data2 )

        column_target = 'column_target'
        column_weight = 'column_weight'
        
        df_1 [ column_target ] = 1
        df_2 [ column_target ] = 0
        
        df_1 [ column_weight ] = weight1 
        df_2 [ column_weight ] = weight2 

        ## merge datasets together
        dataset = pd.concat ( [ df_1 , df_2 ] , axis = 0 ).reset_index ( drop = True )

        X       = dataset . drop ( columns = [ column_target , column_weight ]  )
        Y       = dataset [ column_target ]
        W       = dataset [ column_weight ]
        N       = len ( dataset )
        
        ## cross-validation
        from sklearn.model_selection import StratifiedKFold
        from sklearn.metrics         import roc_auc_score 
        cv        = StratifiedKFold ( n_splits = self.n_splits , shuffle=True , random_state = 42 )
        oof_preds = numpy.zeros( N )
        
        from sklearn.ensemble import GradientBoostingClassifier

        for fold, ( train_idx , val_idx ) in enumerate ( cv.split ( X, Y ) ):
            
            X_train , Y_train , W_train = X.iloc [ train_idx ] , Y.iloc [ train_idx ] , W.iloc [ train_idx ]
            X_val   , Y_val   , W_val   = X.iloc [   val_idx ] , Y.iloc [   val_idx ] , W.iloc [   val_idx ]

            # train the model

            model = GradientBoostingClassifier ( **self.params )
            model.fit ( X_train , Y_train , sample_weight = W_train )
            
            # Предсказываем вероятность принадлежности к классу 1 (выборке P)
            oof_preds [ val_idx ] = model.predict_proba ( X_val ) [ : , 1 ]
            
        #  weighted ROC-AUC
        auc_score = roc_auc_score ( Y , oof_preds , sample_weight = W )

        ## 
        return  100 * ( 1.0 - 2.0 * auc_score ) ** 2 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import pandas # =======================================================
        # =====================================================================
    except ImportError : # ====================================================
        # =====================================================================
        logger.error ( 'pandas cannot be imported!' ) # =======================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================                       
        from sklearn.ensemble        import HistGradientBoostingClassifier
        from sklearn.model_selection import StratifiedKFold
        from sklearn.metrics         import roc_auc_score
        # =====================================================================
    except ImportError : # ====================================================
        # =====================================================================
        logger.error ( "sklearn cannot be imported!" ) # ======================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import lightgbm
        # =====================================================================
    except ImportError : # ====================================================
        # =====================================================================
        logger.error ( "lightgbm cannot be imported!" ) # =====================
        

# =============================================================================
##                                                                      The END 
# =============================================================================
