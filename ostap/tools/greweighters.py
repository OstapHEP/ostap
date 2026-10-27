#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/greweghters.py
#  Density ratio reweighters based on HepML
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @Date  2026-07-18
# =============================================================================
""" Density ratio reweighters based on
 - HepML.
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    'GBReweighter'                   , # Reweighter based on GBReweighter from hep_ml    
) 
# =============================================================================
from   ostap.core.ostap_types   import num_types, string_types, integer_types  
from   ostap.utils.core         import typename
from   ostap.utils.basic        import numcpu, NoContext
from   ostap.logger.utils       import map2table_ex
from   ostap.logger.pretty      import nice_print 
from   ostap.utils.progress_bar import progress_bar, text_bar
from   ostap.tools.reweighter   import Reweighter
from   ostap.stats.utils        import ( weight_trivial     ,
                                         valid_weight       ,
                                         valid_data_shape   ,
                                         compatible_weights , 
                                         num_features       ,
                                         num_samples        ,
                                         check_all          ,
                                         nEff               )
from   ostap.stats.tools            import hasSkLearn 
from   ostap.tools.base_reweighters import DensityReweighter  
import ostap.logger.symbols     as     S
import numpy, abc, math, gc  
# =============================================================================
# Logging setup
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.greweighters' )
else                       : logger = getLogger( __name__ )
# =============================================================================
has_sklearn = hasSkLearn() 
# =============================================================================
method_GBRW  = 'HepML/%s' % ( S.wood       if S.wood        else 'GBRW'     )
# 
DEFAULT_ESTIMATORS   =   300 # High capacity: deep ensemble with low learning rate
MAX_DEPTH            =     6 # Deep tree capacity for high dimensions / large statistics
LEARNING_RATE        =  0.01 # Small step for stable density ratio convergence

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
            "n_estimators"         : DEFAULT_ESTIMATORS , 
            "learning_rate"        : LEARNING_RATE      ,
            "max_depth"            : MAX_DEPTH          ,
            "min_samples_leaf"     : 30                 ,
            "loss_regularization"  : 5.0                ,
            ##
            "gb_args"              : { "subsample"      : 0.95 ,
                                    "max_features"   : 0.95 }
        }
        
        config.update ( params ) 

        check_all ( data1   = original          ,
                    data2   = target            ,
                    weight1 = original_weight   ,
                    weight2 = target_weight     ,
                    where   = typename ( self ) )   

        ## for small number of features ignore 'max_features' argument 
        if num_features ( target ) <= 6 :
            gb_args = config.get ( 'gb_args' , {} )
            gb_args.pop ( 'max_features' , None ) 
            config  [ 'gb_args'      ] = gb_args 

        ## for small number of samples ignore 'subsample' argument 
        n_rows = num_samples ( original ) + num_samples ( target ) 
        if 0 < self.n_splits and n_rows <= 10000 :
            gb_args = config.get ( 'gb_args' , {} )
            gb_args.pop ( 'subsample' , None ) 
            config  [ 'gb_args'   ] = gb_args 
                       
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
               
        with logAttention() if self.silent else NoContext () :
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

    
        self.__original_reweighted_weights = None
        
        if store_original_weights :
            
            factors = self.__reweighter.predict_weights (
                original , 
                original_weight = original_weight )
            
            ## final weights 
            final_weight  = factors if weight_trivial ( original_weight ) else factors * original_weight
            
            ## scale them to the same sum: "preserve the mass"
            final_weight  = self.preserve_mass ( original , original_weight , final_weight ) 
            
            self.__original_reweighted_weights = final_weight 

            ## produce summary of weights 
            if self.verbose : self.weight_summary ( original        ,
                                                    original_weight ,
                                                    final_weight    ) 
            
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

        ## final weights 
        final_weight = factors if weight_trivial ( original_weight ) else factors * original_weight
        final_weight = self.preserve_mass ( original , original_weight , final_weight ) 
        
        ## produce summary table of weights 
        if self.verbose : self.weight_summary ( original        ,
                                                original_weight ,
                                                final_weight    ) 
    
        return final_weight 
    

# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    from ostap.stats.tools import ( hasSkLearn  , hasHepML    )

    if not hasSkLearn  ( False ) : logger.warning  ( "No SkLearn  available!" ) 
    if not hasHepML    ( False ) : logger.warning  ( "No HepML    available!" ) 

# =============================================================================
#                                                                       The END 
# =============================================================================
