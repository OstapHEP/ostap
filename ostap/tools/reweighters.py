#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with various "reweighters"
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-18
# =============================================================================
""" Module with variosu `reweighters'
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    ) 
# =============================================================================
from   ostap.math.math_base  import weight_trivial
from   ostap.utils.basic     import typename, numcpu, num_jobs 
import numpy, abc, warnings  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighters' )
else                       : logger = getLogger( __name__ )
# =============================================================================
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
        
        assert isinstance ( n_splits , int ) and 0 <= n_splits , \
            "Invalid `n_splits':%s" % n_splits 
        
        params [ 'max_depth' ] = params.pop ( 'max_depth' , 5 if n_splits else 3 ) 
        self.__n_splits    = n_splits
        self.__method      = method
        self.__params      = params
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

        from sklearn.metrics import roc_auc_score
        auc          = roc_auc_score ( Y , predictions , sample_weight = weights )
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
        
        conf = {
            'n_estimators'  : 200 , 
            'learning_rate' : 0.1 , 
            'verbose'       : -1 if silent else 0       
        }
        ## Attention! 
        params [ 'n_jobs' ] = num_jobs ( params , numcpu () - 1 )
        conf.update ( params )
        
        ## initiailze the baze 
        Reweighter_base.__init__ ( self ,
                                   original         = original       ,
                                   target          = target          ,
                                   original_weight = original_weight ,
                                   target_weight   = target_weight   ,
                                   ##
                                   n_splits        = n_splits ,
                                   silent          = silent   ,
                                   method          = "Reweighter/LigthGBM" , **conf )

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
        
        conf = { 'n_estimators'  : 200       , 
                 'learning_rate' : 0.10      , 
                 'eval_metric'   : 'logloss' , 
                 'verbosity'     : 0 if silent else 1       
                }
        ## Attention! 
        params [ 'n_jobs' ] = num_jobs ( params , numcpu () - 1 )
        conf.update ( params )
        
        ## initiailze the baze 
        Reweighter_base.__init__ ( self ,
                                   original         = original       ,
                                   target          = target          ,
                                   original_weight = original_weight ,
                                   target_weight   = target_weight   ,
                                   ##
                                   n_splits        = n_splits ,
                                   silent          = silent   ,
                                   method          = "Reweighter/XGBoost" , **conf )

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
## @class Reweighter_CAT
#  Actual Reweighter class based on CatBoost
class Reweighter_CAT(Reweighter_base):
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

        conf = { 'iterations'    : 200 , 
                 'learning_rate' : 0.1 , 
                 'verbose'       : 0 if silent else 1 }
        ## Attention! 
        params [ 'thread_count' ] = num_jobs ( params , numcpu () - 1 )
        conf.update ( params )
        
        ## initialize the baze 
        Reweighter_base.__init__ ( self ,
                                   original         = original       ,
                                   target          = target          ,
                                   original_weight = original_weight ,
                                   target_weight   = target_weight   ,
                                   ##
                                   n_splits        = n_splits ,
                                   silent          = silent   ,
                                   method          = "Reweighter/CatBoost" , **conf )
        
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
