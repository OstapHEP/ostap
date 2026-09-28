#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Utilities for advanced reweighting
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-09-22
# =============================================================================
""" Utilities for advanced reweighting
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    # =========================================================================
    'Reweighter'        , ## the abstract base clss for advanced reweighters
    'CascadeReweighter' , ## Cascade/Iterative reweighter 
    # =========================================================================
) # ===========================================================================
# =============================================================================
from   ostap.utils.core       import typename
from   ostap.core.ostap_types import dictlike_types 
from   ostap.stats.utils      import weight_trivial, num_samples, num_features, check_all  
from   ostap.utils.config     import Config
from   itertools              import zip_longest 
import numpy, abc 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighter' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## @class Reweighter
#  Abstract base class class for Advanced Reweighters: 
#  - custom family of DensityReweighters 
#  - GBReweighter from hep_ml
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-09-22 
class Reweighter(abc.ABC,Config) :
    """ Abstract base class for Advanced Reweighters: 
    - Custom family of DensityReweighters 
    - GBReweighter from hep_ml
     """
    def __init__ ( self                    , * , 
                   original                ,
                   target                  ,
                   original_weight = None  ,
                   target_weight   = None  , 
                   silent          = False ,
                   random_state    = None  , **params  ) :
        """ Abstract base class for Advanced Reweigters
        """

        if not random_state is None and not isinstance ( random_state , int ) :
            raise TypeError  ( "Invalid `random_state' type %s" % typename ( random_state ) )
        
        self.__random_state = random_state

        # =====================================================================
        ## Check the input data
        # ====================================================================
        check_all ( original        ,
                    target          ,
                    original_weight ,
                    target_weight   , typename ( self ) )

        self.__n_features = num_features ( target ) 

        # =====================================================================
        ## initialize the base 
        # =====================================================================
        Config.__init__ ( self , silent = silent , **params ) 

    @property
    @abc.abstractmethod
    def method ( self ) :
        """`method` : underlying method/engine"""
        pass 

    @property
    def random_state ( self ) :
        """`random_state` : random number seed for shuffling/splitting/k-fold/..."""
        return self.__random_state

    @property
    def n_features ( self ) :
        """`n_features` : number of features for this reweighter"""
        return self.__n_features

    @property 
    def config ( self ) :
        """`config` : Reweighter configuraton"""
        conf = {}
        conf.update ( super().config  )
        conf [ 'method'       ] = self.method
        conf [ 'random_state' ] = self.random_state 
        conf [ 'n_features'   ] = self.n_features 
        return conf 
                    
    # =========================================================================
    ## Get/predict new weights for (new) original
    @abc.abstractmethod 
    def weights ( self                   ,
                  original               ,
                  original_weight = None ) :
        """ Get/predict  new weights for (new) original
        """
        pass

    # =========================================================================
    ## alias 
    def new_weights ( self                   ,
                      original               ,
                      original_weight = None ) :
        """ Get/predict  new weights for (new) original
        """
        return self.weights ( original , original_weight )

    # =========================================================================
    ## alias 
    def get_weights ( self                   ,
                      original               ,
                      original_weight = None ) :
        """ Get/predict  new weights for (new) original
        """
        return self.weights ( original , original_weight )
    
    # =========================================================================
    ## Get/predict new weights for (new) original
    def __call__ ( self                   ,
                   original               ,
                   original_weight = None ) :
        """ Get/predict  new weights for (new) original
            """
        return self.weights ( original , original_weight )
     

# =============================================================================
## @class CascadeReweighter
#  Cascade density reweighter supporting heterogeneous estimators with smart padding.
#  Sequentially trains multiple estimators, passing accumulated weights 
#  from previous stages to eliminate severe underfitting.
#  @author Personal AI Collaborator
#  @date 2026-09-28
class CascadeReweighter(Reweighter):
    """ Cascade density reweighter supporting heterogeneous estimators with smart padding.
    Sequentially trains multiple estimators, passing accumulated weights 
    from previous stages to eliminate severe underfitting.
    """    
    # =========================================================================
    ## Constructor for CascadeDensityReweighter
    #  @param stages List of tuples/dictionaries defining (ReweighterClass, config_dict) for each stage
    #  @param classes List of reweighter classes for separate configuration mode
    #  @param configs List of configuration dictionaries corresponding to classes
    #  @param original Original dataset to be reweighted
    #  @param target Target dataset to match against
    #  @param original_weight Initial weights for the original dataset
    #  @param target_weight Initial weights for the target dataset
    #  @param silent Suppress output/logging if True
    #  @param random_state Random seed for reproducibility
    #  @param params Additional keyword arguments for base classes
    def __init__ ( self, *, 
                   stages          = ()    ,
                   classes         = ()    ,
                   configs         = ()    ,
                   original                , 
                   target                  , 
                   original_weight = None  , 
                   target_weight   = None  , 
                   silent          = False , 
                   random_state    = None  , 
                  **params ) :

        ## list of reweigters 
        self.__reweighters = []
        
        # =====================================================================
        # Initialize the base class (performs check_all validations and sets n_features)[cite: 1]
        super().__init__( original        = original        ,
                          target          = target          , 
                          original_weight = original_weight , 
                          target_weight   = target_weight   , 
                          silent          = silent          , 
                          random_state    = random_state    , **params)
        
        if   isinstance ( stages  , type ) and issubclass ( stages  , Reweighter ) : stages  = [ ( stages, {} ) ]
        elif isinstance ( stages  , Reweighter )                                   : stages  = [   stages       ]
        if   isinstance ( classes , type ) and issubclass ( classes , Reweighter ) : classes = [   classes      ] 
        elif isinstance ( classes , Reweighter )                                   : classes = [   classes      ] 

        if not stages and not classes : 
            raise ValueError ( "Neither 'stages' nor 'classes' are provided!")
        
        the_stages = []

        ## copy stages 
        for stage in stages :
            if   isinstance ( stage , type ) and issubclass ( stage , Reweighter ) : stage = stage, {}
            elif isinstance ( stage , Reweighter ) :
                pars = { 'silent' : stage.silent , 'random_state' : stage.random_state }                
                pars.update ( stage.params )                
                stage = type ( stage ) , pars
                
            cl , cnf = stage 
            if not cnf : cnf = {} 
            if isinstance  ( cl , type ) and not issubclass ( cl  , Reweighter ) : raise TypeError ( "Invalid `classes` type!" )
            if not isinstance ( cnf , dictlike_types )                           : raise TypeError ( "Invalid `stages` type!"  )
            the_stages.append ( ( cl , cnf ) ) 
            
        last = None 
        for cl, cnf in zip_longest ( classes , configs , fillvalue = None ) :
            if cl is None : cl  = last
            if not cnf    : cnf = {}
            if  not isinstance ( cnf , dictlike_types ) : raise TypeError ( "Invalid `configs` type!" )

            ## instance ? 
            if cl and isinstance ( cl , Reweighter ) :
                pars = { 'silent' : cl.silent , 'random_state' : cl.random_state }                
                pars.update ( cl.params )
                pars.update ( cnf  ) 
                cl  = type  ( cl   )
                cnf = pars
            elif not issubclass ( cl , Reweighter      ) : raise TypeError ( "Invalid `classes` type!" )
            ## 
            the_stages.append ( ( cl , cnf ) ) 
            last = cl


        n_orig = num_samples ( original )
        n_targ = num_samples ( target   )
        
        w_orig = numpy.ones  ( n_orig , dtype = numpy.float32 ) if original_weight is None else numpy.array ( original_weight , dtype = numpy.float32 )
        w_targ = numpy.ones  ( n_targ , dtype = numpy.float32 ) if target_weight   is None else numpy.array ( target_weight   , dtype = numpy.float32 )

        
        # Train heterogeneous cascade stages sequentially with safety checks
        for i , stage_info in enumerate ( the_stages ) :
            
            REWEIGHTER , rwcnf = stage_info 
            
            conf = { 'silent' : self.silent , 'random_state' : self.random_state }
            conf.update ( params )
            conf.update ( rwcnf  ) 
            
            reweighter = REWEIGHTER ( original        = original ,
                                      target          = target   ,
                                      original_weight = w_orig   ,
                                      target_weight   = w_targ   , **conf )

            # 2. Get weights from this stage, explicitly passing current w_orig 
            # to the weights function 
            w_orig = reweighter . weights ( original , original_weight = w_orig )
            
            # Store the trained stage instance
            self.__reweighters.append ( reweighter )

        ## convert to tuple 
        self.__reweighters = tuple ( self.__reweighters )

    # =========================================================================
    ## @property
    #  @brief Underlying methods/engines of all cascade stages
    #  @return List of method names/identifiers for each stage
    @property
    def method(self):
        """Implementation of the abstract 'method' property"""
        from ostap.logger.symbols import bold_arrow_right 
        s = ' %s ' % bold_arrow_right  
        return s.join ( [ r.method for r in self.__reweighters ] )

    # =========================================================================
    ## @brief Get/predict new weights for (new) original dataset
    #  @param original Dataset to compute weights for
    #  @param original_weight Base weights to scale
    #  @return Final compounded weights as a numpy array
    def weights ( self     ,
                  original ,
                  original_weight = None ):
        """ Implementation of the abstract 'weights' method.
        Applies each stage sequentially, passing the running weight vector.
        """
        w_final = original_weight
        
        for stage in self.__reweighters:
            w_final = stage.weights ( original ,
                                      original_weight = w_final )
            
        return w_final

# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        
# =============================================================================
##                                                                      The END 
# =============================================================================
