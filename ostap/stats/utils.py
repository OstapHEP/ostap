#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Small utilities for statistical calculations 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-18
# =============================================================================
""" Small utilities for statistical calculations 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    'weight_trivial'      , ## Is weight(1D numpy array) "trivial" ?
    'valid_weights_shape' , ## Valid weights shape ? 
    'valid_data_shape'    , ## Valid data shape ?
    'valid_weight'        , ## Valid weights (valid shape and positibe sum)
    'num_features'        , ## Number of features for training data
    'num_samples'         , ## Number of samples/events (rows) in the dataset.
    'compatible_shapes'   , ## Check if input datasets have compatible shapes 
    'compatible_weights'  , ## Check for data and weights compatibility
    'check_all'           , ## check ALL
    'nEff'                , ## Compute effective sample size (Kish's design effect formula)
) 
# =============================================================================
from   ostap.core.ostap_types import num_types, numpy_buffer_types, sized_types 
from   ostap.utils.core       import typename 
import numpy
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.utils' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## Trivial weight ? 
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
    if   weight is None                             : return True
    elif isinstance ( weight , num_types          ) : return 1 == weight 
    elif isinstance ( weight , numpy_buffer_types ) : return 1 == weight 
    elif isinstance ( weight , numpy.ndarray      ) : return numpy.all ( weight == 1 ) 
    elif isinstance ( weight , sized_types        ) : return all ( 1 == w for w in weight )
    return False

# ==============================================================================
## Compute effective sample size (Kish's design effect formula)
#  taking into account sample weights (including negative sPlot weights).
def nEff ( X , W = None ) :
    """ Compute effective sample size (Kish's design effect formula)
        taking into account sample weights (including negative sPlot weights).
    """
    ns = num_samples  ( X )
    if weight_trivial ( W  ) : return float ( ns )
    
    w_arr  = numpy.asarray ( W     , dtype = numpy.float32 )
    sum_w  = numpy.sum     ( w_arr , dtype = numpy.float64 )
    sum_w2 = numpy.dot     ( w_arr.astype (  numpy.float64 ) , w_arr ) 
    
    return ( sum_w** 2 ) / sum_w2 if 0 < sum_w2 else float ( ns )

# =============================================================================
## Check if weights array has a valid 1D/column vector shape
# - Supports NumPy arrays, lists, tuples, and custom containers.
def valid_weights_shape ( weights ) :
    """ Check if weights have a valid 1D array or (N, 1) column shape.
    - Supports NumPy arrays, lists, tuples, and custom containers.
    """
    if weight_trivial ( weights ) : return True

    # Extract shape attribute if present (NumPy arrays, Tensors, DataFrames)
    shape = getattr ( weights , 'shape' , None )

    if shape is not None:
        ndim = len ( shape )
        # Accept 1D vectors (N,) or 0D scalars
        if ndim <= 1 : return True
        # Accept (N, 1) column vectors, reject any multi-column (N, >1) or >2D arrays
        if ndim == 2 and shape[1] == 1: return True
        return False
    
    # For standard Python sequences (lists, tuples): check nested dimensions
    if hasattr ( weights , '__len__' ) and 0 < len  ( weights ) :
        first_element = weights [ 0 ]
        # Reject lists of lists/tuples with length > 1
        if isinstance ( first_element, sized_types ) :
            return 1 == len ( first_element )

    return True

# =============================================================================
## Check if dataset has a valid structure/shape
# - Supports NumPy arrays, DataFrames, PyROOT C++ containers, and sequences.
def valid_data_shape(data):
    """ Check if data has a valid 1D, 2D, or Structured Array shape.
    - Supports NumPy arrays, DataFrames, PyROOT C++ containers, and sequences.
    """
    if data is None: return False

    # Check shape attribute if present (NumPy arrays, Pandas, Tensors)
    shape = getattr ( data , 'shape' , None )
    if shape is not None:
        ndim = len ( shape )
        # Allow 1D vectors (N,) and 2D tabular features (N, F)
        return True if 1 <= ndim <= 2 else False 

    # For standard Python sequences / PyROOT containers
    return True if hasattr ( data , '__len__' ) else False 

# =============================================================================
## Number of features for training data
#  - Safe extraction of n_features (supports DataFrame, 2D array, 1D vector, Structured Array)
def num_features ( X ) :
    """ Number of features for training data
    - Safe extraction of n_features (supports DataFrame, 2D array, 1D vector, Structured Array)
    """
    # 1. Structured Arrays (check field names in dtype)
    dtype = getattr ( X , 'dtype' , None )
    if dtype is not None and dtype.names : return len ( dtype.names )

    # 2. Standard NumPy arrays, DataFrames, etc.
    shape = getattr ( X , 'shape' , None )
    if shape is not None and 1 < len ( shape ) : return shape [ 1 ]

    # 3. 1D vectors and fallback
    return 1

# =============================================================================
## Returns the number of samples/events (rows) in the dataset.
#  Supports NumPy arrays, Pandas DataFrames/Series, SciPy sparse matrices,
#  PyTorch/TensorFlow tensors, C++ containers (via PyROOT), and standard Python sequences.
def num_samples(X):
    """ Number of samples/events (rows) in the dataset.

    Supports NumPy arrays, Pandas DataFrames/Series, SciPy sparse matrices,
    PyTorch/TensorFlow tensors, C++ containers (via PyROOT), and standard Python sequences.
    """
    if X is None : return 0

    # Handle standard Python sequences, DataFrames, ROOT objects, PyROOT C++ vectors
    if hasattr ( X, "__len__" ) : return len ( X )

    # Handle objects with .shape attribute (NumPy, PyTorch, TensorFlow, SciPy )
    shape = getattr ( X , "shape" , None )
    if shape is not None:
        if 0 == len ( shape ) : return 0 # 0D scalar has 0 samples 
        return int( shape  [ 0 ] )

    raise TypeError ( "Unsupported data structure for determining sample count: %s" % typename ( X ) )

# =============================================================================
## Check if input dataset shapes have compatible feature dimensions
def compatible_shapes ( ds1 , ds2 ) :
    """ Check if ds1 and ds2 have compatible feature dimensions"""
    
    # Extract dtype for structured arrays (if applicable)
    dt1 = getattr ( ds1 , 'dtype' , None )
    dt2 = getattr ( ds2 , 'dtype' , None )

    # If both are structured arrays, compare their field names/types
    if dt1 is not None and dt1.names and dt2 is not None and dt2.names:
        return dt1.names == dt2.names

    s1 = numpy.shape ( ds1 )
    s2 = numpy.shape ( ds2 )

    # Validate non-empty datasets
    if not ( s1 and s2 and s1 [ 0 ] > 0 and s2 [ 0 ] > 0 ):
        return False

    # Get feature dimensions (ignoring event count at index 0)
    dim1 = s1 [ 1 : ] if 1 < len ( s1 ) else ()
    dim2 = s2 [ 1 : ] if 1 < len ( s2 ) else ()

    return dim1 == dim2

# =============================================================================
## Check if dataset and weights are compatible
#  - Ensures weights are 1D (or 2D column vector) matching the data length.
def compatible_weights ( data , weights = None ) :
    """Check for data and weights compatibility.
    - Ensures weights are 1D (or 2D column vector) matching the data length.
    """
    if weight_trivial ( weights ) : return True

    # 1. Reject data shapes first
    if not valid_data_shape    ( data    ) : return False
    
    # 2. Reject invalid weight matrix shapes second 
    if not valid_weights_shape ( weights ) : return False

    # 3. Check event count compatibility using num_samples
    return num_samples ( data ) == num_samples ( weights )

# ===========================================================================
## Valid weight?
#  - trivial or valid shape 
#  - trivial or sum is positive
def valid_weight ( weight ) :
    """ Valid weight?
    - trivial or valid shape 
    - trivial or sum is positive
    """
    if     weight_trivial      ( weight ) : return True
    if not valid_weights_shape ( weight ) : return False
    
    # =========================================================================
    # Safely evaluate sum for NumPy arrays, lists, ROOT containers, etc.
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        return 0 < float ( numpy.sum ( weight ) )
        # ===================================================================== 
    except Exception : # ======================================================
        # =====================================================================
        return False

# ============================================================================
## check data and weights
#  - <code>data1</code>   is valid 
#  - <code>data2</code>   is valid 
#  - <code>data1</code>   is "compatible" with <code>data2</code>
#  - <code>weight1</code> is valid
#  - <code>weight2</code> is valid
#  - <code>data1</code>   and <code>weight1</code> are "compatible"
#  - <code>data2</code>   and <code>weight2</code> are "compatible"
def check_all ( data1   ,
                data2   ,
                weight1 = None        ,
                weight2 = None        ,
                where   = "check_all" ) :
    
    if not valid_weight       ( weight1         ) : raise TypeError ( "%s: invalid `weight1`" % where ) 
    if not valid_weight       ( weight2         ) : raise TypeError ( "%s: invalid `weight2`" % where ) 
    if not valid_data_shape   ( data1           ) : raise TypeError ( "%s: invalid `data1`"   % where ) 
    if not valid_data_shape   ( data2           ) : raise TypeError ( "%s: invalid `data2`"   % where ) 
    if not compatible_shapes  ( data1 , data2   ) : raise TypeError ( "%s: incompatible data shape`"  % where ) 
    if not compatible_weights ( data1 , weight1 ) : raise TypeError ( "%s: incompatible `data1/weight1`"  % where ) 
    if not compatible_weights ( data2 , weight2 ) : raise TypeError ( "%s: incompatible `data2/weight2`"  % where ) 
    return True 
    
# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        
# =============================================================================
##                                                                      The END 
# =============================================================================
