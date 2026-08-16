#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Small utilities for statistical calcualtios 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-18
# =============================================================================
""" Small utilities for statistical calcualtios 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    'weight_trivial'    , ## Is weight(1D numpy array) "trivial" ?
    'num_features'      , ## Number of features for training data
    'num_samples'       , ## Number of samples/events (rows) in the dataset.
    'compatible_shapes' , ## Check if input datasets have compatible shapes 
) 
# =============================================================================
from   ostap.core.ostap_types import num_types, numpy_buffer_types 
import numpy
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.utils' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## Check if input dataset shapes have compatible feature dimensions
def compatible_shapes ( ds1 , ds2 ) :
    """ Check if ds1 and ds2 have compatible feature dimensions """
    s1 = numpy.shape ( ds1 )
    s2 = numpy.shape ( ds2 )

    if not ( len ( s1 ) > 0 and len ( s2 ) > 0 and s1 [ 0 ] > 0 and s2 [ 0 ] > 0 ) :
        return False

    # Get feature dimension (1 for 1D arrays, shape[1:] for N-D arrays)
    dim1 = s1 [ 1 : ] if len ( s1 ) > 1 else ( 1 , )
    dim2 = s2 [ 1 : ] if len ( s2 ) > 1 else ( 1 , )

    return dim1 == dim2

# =============================================================================
## Trvial weight ? 
#  - None
#  - one 
#  - all ones
def weight_trivial ( weight ) :
    """ Trivial weight ? 
    - None
    - one 
    - all ones
    """
    ## 
    if    weight is None                             : return True
    elif  isinstance ( weight , num_types          ) : return 1 == weight 
    elif  isinstance ( weight , numpy_buffer_types ) : return 1 == weight 
    elif  isinstance ( weight , numpy.ndarray      ) : return numpy.all ( weight == 1 ) 
    return False
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


# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        
# =============================================================================
##                                                                      The END 
# =============================================================================
