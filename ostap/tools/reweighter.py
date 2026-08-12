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
from   ostap.utils.basic    import typename
from   ostap.math.math_base import weight_trivial
from   ostap.utils.config   import Config
import numpy, os, abc 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighter' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## Number of features for training data
#  Safe extraction of n_features (supports DataFrame, 2D array, 1D vector)
def num_features ( X ) :
    """ Number of features for training data
    - Safe extraction of n_features (supports DataFrame, 2D array, 1D vector)
    """
    return X.shape [ 1 ] if hasattr ( X , 'shape' ) and 1 < len ( X.shape ) else 1
# ============================================================================
## Returns the number of samples/events (rows) in the dataset.
#  Supports NumPy arrays, Pandas DataFrames/Series, SciPy sparse matrices,
#  PyTorch/TensorFlow tensors, and standard Python collections.
def num_samples ( X ) :
    """ Returns the number of samples/events (rows) in the dataset.

    Supports NumPy arrays, Pandas DataFrames/Series, SciPy sparse matrices,
    PyTorch/TensorFlow tensors, and standard Python collections.
    """
    if X is None : return 0

    # Handle objects with a .shape attribute (NumPy, Pandas, SciPy, PyTorch, TensorFlow)
    if hasattr ( X ,  "shape" ) :
        shape = X.shape
        if 0 == len ( shape ) :  return 1
        return int ( shape [ 0 ] )

    # Handle lists, tuples, and other standard Python sequences
    if hasattr ( X , "__len__" ) : return len ( X )

    raise TypeError ( "Unsupported data structure for determining sample count: %" % typename ( X )  )
# =============================================================================
## @class Reweighter
#  Abstract base class class for Advanced Reweighters
#  - GBReweighter from hep_ml
#  - custom family of DensityReweighters 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-09-22 
class Reweighter(abc.ABC,Config) :
    """ Abstract base class for Advanced Reweighters 
    - GBReweighter from hep_ml
    - custom family of DensityReweighters 
    """
    def __init__ ( self                    , * , 
                   original                ,
                   target                  ,
                   original_weight = None  ,
                   target_weight   = None  , 
                   silent          = False ,
                   random_state    = 42    , **params  ) :
        """ Abstract base class for Advanced Reweigters
        """

        if not isinstance ( random_state , int ) :
            raise TypeError  ( "Invalid `random_state' type %s"    % typename ( random_state ) )
        
        self.__random_state = random_state

        # ==========================================================================
        ## Check input data
        # ==========================================================================
        
        shape1  = original.shape
        shape2  = target.shape 

        nf_orig = num_features ( original )
        nf_targ = num_features ( target   )
        
        if nf_orig != nf_targ : raise TypeError ( "Inconsistent original/targer shapes: %s vs %s " % ( original.shape , target.shape ) )

        ns_orig = num_samples  ( original )
        ns_targ = num_features ( target   )
        
        if original_weight is None or len ( original_weight ) == ns_orig : pass
        else : raise TypeError ( "Inconsistent original/weight : %s vs %s " % ( ns_orig , len ( original_weight ) ) )

        if target_weight   is None or len ( target_weight   ) == ns_targ : pass
        else : raise TypeError ( "Inconsistent target/weight   : %s vs %s " % ( ns_targ , len ( target_weight ) ) )

        if not weight_trivial ( target_weight ) :
            sw = numpy.sum ( target_weight )
            if sw <= 0 : logger.error ( "Sum of target   weigths is non-positive: %s" % float ( sw ) )
            
        if not weight_trivial ( original_weight ) :
            sw = numpy.sum ( original_weight )
            if sw <= 0 : logger.error ( "Sum of original weigths is non-positive: %s" % float ( sw ) )
            
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
    def config ( self ) :
        """`config` : Reweighter configuraton"""
        conf = {}
        conf.update ( super().config  )
        conf [ 'method'       ] = self.method
        conf [ 'random_state' ] = self.random_state 
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
