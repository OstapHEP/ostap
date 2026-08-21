#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Utilities for advanced reweighting
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-09-22
# =============================================================================
""" Utilities for advanced reweighting
- GBReweighter from hep_ml
- custom family of DensityReweighters 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    'num_features' , ## number of features for trainig data
    'num_samples'  , ## number of samples/events (rows) in the dataset.
    'Reweighter'   , ## the abstract base clss for advanced reweighters 
) 
# =============================================================================
from   ostap.utils.core     import typename
from   ostap.stats.utils    import weight_trivial, num_samples, num_features, check_all  
from   ostap.utils.config   import Config
import numpy, os, abc 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighter' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## @class Reweighter
#  Abstract base class class for Advanced Reweighters
#  - GBReweighter from hep_ml
#  - custom family of DensityReweighters 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-09-22 
class Reweighter(abc.ABC,Config) :
    """ Abstract base class for Advanced Reweighters: 
    - GBReweighter from hep_ml
    - Custom family of DensityReweighters 
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

        # ==========================================================================
        ## Check input data
        # ==========================================================================
        check_all ( original        ,
                    target          ,
                    original_weight ,
                    target_weight   , typename ( self ) )

        self.__n_features = num_features ( target ) 

        # ==========================================================================
        ## initialize the base 
        # ==========================================================================
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
                    
    # ==============================================================================
    ## Get/predict new weights for (new) original
    @abc.abstractmethod 
    def weights ( self                   ,
                  original               ,
                  original_weight = None ) :
        """ Get/predict  new weights for (new) original
        """
        pass

    # ==============================================================================
    ## alias 
    def new_weights ( self                   ,
                      original               ,
                      original_weight = None ) :
        """ Get/predict  new weights for (new) original
        """
        return self.weights ( original , original_weight )

    # ==============================================================================
    ## alias 
    def get_weights ( self                   ,
                      original               ,
                      original_weight = None ) :
        """ Get/predict  new weights for (new) original
        """
        return self.weights ( original , original_weight )
    
    # ==============================================================================
    ## Get/predict new weights for (new) original
    def __call__ ( self                   ,
                   original               ,
                   original_weight = None ) :
        """ Get/predict  new weights for (new) original
            """
        return self.weights ( original , original_weight )
     
# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        
# =============================================================================
##                                                                      The END 
# =============================================================================
