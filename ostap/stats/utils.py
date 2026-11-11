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
    'np2vct'              , ## numpy arrays into SVectorWithError
    'eff_size'            , ## effective size of two samples
    ##
    'full_header'         , ## the header for statistical table (full) 
    'short_header'        , ## the header for statistical table (short) 
    't_header'            , ## the header for statistical table (only t-values)
    'w_header'            , ## "header" for min/max weights
    ##
    'symbol_LGBM'         , ## symbold for LigthGBM 
    'symbol_XGB'          , ## symbol for XGBoost 
    'symbol_CATB'         , ## symbol for CatBoost 
    'symbol_HGBC'         , ## symbol for HistoGradinetBoost 
    'symbol_GBC'          , ## symbol for GradientBoost 
    'symbol_RF'           , ## symbol for RandomForest 
    'symbol_TORCH'        , ## symbol for PyTorch 
    'symbol_KERAS'        , ## symbol for Keras
    'symbol_LOGREG'       , ## symbol fo rLogistic Regression 
    'symbol_USTAT'        , ## symbol for uStat
    'symbol_CHI2'         , ## symbol for chi2
    ##
    'symbol_MIX'          , ## symbol for NMixed-Sample
    'symbol_PPD'          , ## symbol for Point-to-Point Dissimilarity 
    'symbol_DNN'          , ## symbol for Distance-to-Nearest Neighbor 
    ## 
    'symbol_KL'           , ## symbol for Kullback-Leibler 
    'symbol_J'            , ## symbol for Jeffrey 
    'symbol_JS'           , ## symbol for Jensen-Shannon 
    'symbol_M'            , ## symbol for Mahalanobis 
    'symbol_T2'           , ## symbol for Hotelling T2 
    'symbol_B'            , ## symbol for Bhattacharyya'               
    'symbol_H'            , ## symbol for Hellinger 
    'symbol_W'            , ## symbol for Wasserstein 
    'symbol_W1'           , ## symbol for Wasserstein-1 
    'symbol_W2'           , ## symbol for Wasserstein-2
    'symbol_R'            , ## symbol for Rosenblatt 
    'symbol_U2'           , ## Symbol for Watson U2 
    'symbol_MW'           , ## Symbol for Mann-Whintey 
    'symbol_S'            , ## Symbol for Savage 
    ##
    'symbol_NLL'          , ## symbol for NLL 
    'symbol_AIC'          , ## symbol for Akainke IC 
    'symbol_BIC'          , ## symbol for Bayesian IS     
    ##
    'symbol_KS'           , ## symbol for Kolmogorov-Smirnov 
    'symbol_K'            , ## symbol for Kuiper 
    'symbol_AD'           , ## symbol for Anderson-Darling 
    'symbol_CVM'          , ## symbol for Cramers-von Mises
    'symbol_BJ'           , ## symbol for Berk-Jones 
    'symbol_ZA'           , ## symbol for Zhang ZA 
    'symbol_ZC'           , ## symbol for Zhang ZC 
    'symbol_ZK'           , ## symbol fo Zhang ZK 
    ##
    'method_LGBM'         , ##            LightGBM method 
    'method_XGB'          , ##             XGBoost method 
    'method_CATB'         , ##            CatBoost method 
    'method_HGBC'         , ##  HistoGradientBoost method 
    'method_GBC'          , ##       GradientBoost method 
    'method_RF'           , ##        RandomForest method 
    'method_TORCH'        , ##             PyTorch method 
    'method_KERAS'        , ##               Keras method
    'method_LOGREG'       , ## Logistic Regression method 
    'method_USTAT'        , ##              uStat method 
    'method_CHI2'         , ##               Chi2 method 
    ##
    'method_MIX'          , ##                 Mixed-Sample method 
    'method_PPD'          , ## Point-to-Point Dissimilarity method 
    'method_DNN'          , ## Distance-to-Nearest Neighbor method
    ## 
    'method_KL'           , ##             Kullback-Leibler method 
    'method_J'            , ##                      Jeffrey method  
    'method_JS'           , ##               Jensen-Shannon method 
    'method_M'            , ##                  Mahalanobis method 
    'method_B'            , ##                Bhattacharyya method 
    'method_H'            , ##                    Hellinger method 
    'method_T2'           , ##                    Hotteling method                       
    'method_W'            , ##                  Wasserstein method
    'method_W1'           , ##                 Wassertein/1 method 
    'method_W2'           , ##                 Wassertein/2 method 
    'method_R'            , ##                   Rosenblatt method 
    'method_U2'           , ##                    Watson U2 method 
    'method_MW'           , ##                Mann-Whitney  method 
    ##
    'method_NLL'          , ## Negative log-likelihood 
    'method_AIC'          , ## Akaike IC
    'method_BIC'          , ## Bayesian IC'
    ##
    'method_KS'           , ## Kolmogorov-Smirnov     
    'method_K'            , ## Kuiper                
    'method_AD'           , ## Anderson-Darling
    'method_CVM'          , ## Cramer-von Mises'        
    'method_BJ'           , ## Berk-Kones'             
    'method_ZA'           , ## zhang ZA 
    'method_ZC'           , ## zhang ZC 
    'method_ZK'           , ## zhang ZK 
) 
# =============================================================================
from   ostap.core.ostap_types import num_types, numpy_buffer_types, sized_types 
from   ostap.utils.core       import typename
import ostap.logger.symbols   as     S 
import numpy
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
logger = getLogger ( 'ostap.stats.utils' )
# ============================================================================
symbol_name  = S.bolditalic 
method_name  = S.italic

symbol_EPOCH = f'{S.repeat} :' if S.show and S.repeat else 'Epoch:'

symbol_LGBM   = symbol_name ( S.light_bulb   if S.show and S.light_bulb     else 'LightGBM' ) 
symbol_XGB    = symbol_name ( S.rocket       if S.show and S.rocket         else 'XGBoost'  ) 
symbol_CATB   = symbol_name ( S.cat_face     if S.show and S.cat_face       else 'CatBoost' ) 
symbol_HGBC   = symbol_name ( S.bar_chart    if S.show and S.bar_chart      else 'HGBC'     ) 
symbol_GBC    = symbol_name ( S.nabla        if S.show and S.nabla          else 'GBC'      ) 
symbol_RF     = symbol_name ( S.palm         if S.show and S.palm           else 'RF'       )
symbol_TORCH  = symbol_name ( S.flashlight   if S.show and S.flashlight     else 'TORCH'    ) 
symbol_KERAS  = symbol_name ( S.postal_horn  if S.show and S.postal_horn    else 'KERAS'    ) 
symbol_LOGREG = symbol_name ( S.ruler        if S.show and S.ruler          else 'LOGREG'   ) 
symbol_USTAT  = symbol_name ( S.pushpin      if S.show and S.pushpin        else 'USTAT'    ) 
symbol_CHI2   = symbol_name ( S.chi2         if S.show and S.chi2           else 'Chi2'     ) 

symbol_MIX   = symbol_name ( S.shuffle      if S.show and S.shuffle        else 'MIX'      ) 
symbol_PPD   = symbol_name ( S.ruler        if S.show and S.ruler          else 'PPD'      ) 
symbol_DNN   = symbol_name ( S.compass      if S.show and S.compass        else 'DNN'      ) 

symbol_KL    = symbol_name ( 'K-L' )
symbol_J     = symbol_name ( 'J'   )
symbol_JS    = symbol_name ( 'J-S' )
symbol_M     = symbol_name ( 'M'   )
symbol_T2    = symbol_name ( f'T{S.sup_two}' if S.show else 'T2'  )
symbol_B     = symbol_name ( 'B'   )
symbol_H     = symbol_name ( 'H'   )
symbol_W     = symbol_name ( f'W{S.sup_two}' if S.show else 'W2'  )
symbol_W1    = symbol_name ( f'W{S.sub_one}' if S.show else 'W_1' )
symbol_W2    = symbol_name ( f'W{S.sub_two}' if S.show else 'W_2' )
symbol_R     = symbol_name ( 'R')
symbol_U2    = symbol_name ( f'U{S.sup_two}' if S.show else 'U2'  )
symbol_MW    = symbol_name ( 'M-W' )
symbol_S     = symbol_name ( 'S'   )

symbol_NLL   = symbol_name ( f'-{S.delta_symbol}log{S.script_L}' if S.show else 'NLL' ) 
symbol_AIC   = symbol_name ( 'AIC' ) 
symbol_BIC   = symbol_name ( 'BIC' ) 

symbol_KS    = symbol_name ( 'K-S' )
symbol_K     = symbol_name ( 'K'   )
symbol_AD    = symbol_name ( 'A-D' )
symbol_CVM   = symbol_name ( 'CvM' )
symbol_BJ    = symbol_name ( 'B-J' )
symbol_ZA    = symbol_name ( f'Z{S.subscript_A}' )
symbol_ZC    = symbol_name ( f'Z{S.subscript_C}' )
symbol_ZK    = symbol_name ( f'Z{S.subscript_K}' )


method_LGBM   = method_name ( 'LightGBM'            )
method_XGB    = method_name ( 'XGBoost'             ) 
method_CATB   = method_name ( 'CatBoost'            ) 
method_HGBC   = method_name ( 'HistoGradientBoost'  ) 
method_GBC    = method_name ( 'GradientBoost'       ) 
method_RF     = method_name ( 'RandomForest'        )
method_TORCH  = method_name ( 'PyTorch'             )
method_KERAS  = method_name ( 'Keras'               ) 
method_LOGREG = method_name ( 'Logistic Regression' ) 
method_USTAT  = method_name ( 'uStat'               )
method_CHI2   = method_name ( S.chi2  if S.show and S.chi2  else 'Chi2' ) 

method_MIX   = method_name ( 'Mixed-Sample'                  )
method_PPD   = method_name ( 'Point-to-Point Dissimilarity'  )
method_DNN   = method_name ( 'Distance-to-Nearest Neighbour' )

method_KL    = method_name ( 'Kullback-Leibler'              )
method_J     = method_name ( 'Jeffrey'                       )
method_JS    = method_name ( 'Jensen-Shannon'                )
method_M     = method_name ( 'Mahalanobis'                   )
method_T2    = method_name ( f'Hotteling T{S.sup_two}'       )
method_B     = method_name ( 'Bhattacharyya'                 )
method_H     = method_name ( 'Hellinger'                     )
method_W     = method_name ( f'Wasserstein W{S.sup_two}'     )
method_W1    = method_name ( f'Wasserstein W{S.sub_one}'     )
method_W2    = method_name ( f'Wasserstein W{S.sub_two}'     )
method_R     = method_name ( 'Rosenblatt'                    )
method_U2    = method_name ( f'Watson U{S.sup_two}'          )
method_MW    = method_name ( 'Mann-Whitney'                  )
method_S     = method_name ( 'Savage'                        )

method_NLL   = method_name ( f'-{S.delta_symbol}log{S.script_L}' if S.show else 'Negative Delta Likelihood' ) 
method_AIC   = method_name ( 'Akaike IC'   )
method_BIC   = method_name ( 'Bayesian IC' )

method_KS    = method_name ( 'Kolmogorov-Smirnov'      )
method_K     = method_name ( 'Kuiper'                  )
method_AD    = method_name ( 'Anderson-Darling'        )
method_CVM   = method_name ( 'Cramer-von Mises'        )
method_BJ    = method_name ( 'Berk-Jones'              )
method_ZA    = method_name ( f"Zhang Z{S.subscript_A}" )
method_ZC    = method_name ( f"Zhang Z{S.subscript_C}" )
method_ZK    = method_name ( f"Zhang Z{S.subscript_K}" )

# ============================================================================
## the headers to statistics table
full_header     = ( S.t_value_expr         , 
                    S.t_mean_expr          , 
                    S.t_rms_expr           ,
                    S.t_min_expr           ,
                    S.t_max_expr           ,
                    S.t_unit_expr          ,
                    S.p_value_percent_expr ,
                    S.n_sigma_expr         )
short_header    = ( S.t_value_expr         , 
                    S.t_unit_expr          ,
                    S.p_value_percent_expr ,
                    S.n_sigma_expr         )
t_header        = S.t_value_expr  , S.t_unit_expr           
w_header        = S.w_min_expr    , S.w_max_expr , S.w_unit_expr
# =============================================================================

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
def nEff ( X , W = None ):
    """Compute the effective sample size using Kish's design effect formula.
    
    This function takes sample weights into account, properly handling negative 
    weights which commonly occur in sPlot background-subtraction techniques.

    Parameters:
    -----------
    X : array-like or dataset container
        The input dataset (used solely to determine the total number of samples).
    W : array-like, optional
        Sample weights corresponding to the dataset. Can contain negative values 
        (e.g., sPlot weights). If None or trivial, unweighted sample size is returned.

    Returns:
    --------
    float
        The effective sample size (N_eff). Falls back to the raw number of samples 
        if weights are trivial or if the sum of squared weights is non-positive.
    """
    # Step 1: Retrieve the raw number of samples from the dataset
    ns = float ( num_samples ( X ) )
    
    # Step 2: Check if weights are trivial (e.g., None, all ones, or empty)
    if weight_trivial ( W ) : return ns

    # Step 3: Convert weights directly to float64 to prevent precision loss.
    # High-precision accumulation is critical here because sPlot weights can have 
    # large dynamic ranges and negative values, making them sensitive to rounding errors.
    w_arr = numpy.asarray ( W , dtype = numpy.float64 )
    
    # Step 4: Compute the sum of weights (Sum of W)
    sum_w = float ( numpy.sum ( w_arr ) ) 
    
    # Step 5: Compute the sum of squared weights (Sum of W^2) 
    # using an optimized dot product for performance.
    sum_w2 = float ( numpy.dot ( w_arr, w_arr ) ) 

    # Step 6: Apply Kish's formula: N_eff = (sum_w)^2 / sum_w2.
    # Guard against division by zero or negative sum of squares (pathological weight cases).
    return float ( sum_w ** 2 ) / sum_w2 if 0  < sum_w2 else ns 

# =============================================================================
## Compute the effective size of two samples
#  \f$ N = \frac{N^{eff}_1N^{eff}_2}{N^{eff}_1+N^{eff}_2}\f$
def eff_size ( data1   , 
               data2   , 
               weight1 = None ,
               weight2 = None ) :
    """ Compute the effective size of two samples
    N = Neff_1 * Neff_2 / ( Neff_1 + Neff_2 )
    """
    n1 = nEff ( data1 , weight1 )
    n2 = nEff ( data2 , weight2 )
    return float ( n1 * n2 / ( n1 + n2 ) ) 
    
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
        return int ( shape  [ 0 ] )

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

# ==============================================================================
## Convert numpy-array statistics into Ostap::SVectorWithError
#  @see Ostap::Math::SVectorWithError
def np2vct ( data , weight = None ) :
    """ Convert numpy-array statistics into `Ostap.Math.SVectorWithError`
    - see `Ostap.Math.SVectorWithError`
    """

    shape     = data.shape
    n , N     = shape
    
    w         = None if weight_trivial ( weight ) else weight

    ## 
    from statsmodels.stats.weightstats import DescrStatsW as DSW 
    dsw       = DSW  ( data , weights = w )
    mean      = dsw.mean
    covmtrx   = dsw.cov

    
    ## load corresponding linear-algebra tricks:
    import ostap.math.linalg
    from   ostap.core.core   import Ostap
    
    MN        = Ostap.Math.SymMatrix ( N )
    
    ## prepare output 
    RT        = Ostap.Math.SVectorWithError [ N ]
    values    = RT.Value      () 
    covs      = RT.Covariance () 
    
    for i in range ( N  ) :
        values [ i     ] = float ( mean    [ i ]       )
        covs   [ i , i ] = float ( covmtrx [ i ] [ i ] )
        for j in range (  0 , i  ) :
            cij = float ( covmtrx [ i ] [ j ]  )
            cji = float ( covmtrx [ j ] [ i ]  )                
            cc  = 0.5 * ( cij + cji ) 
            covs [ i , j ] = cc
            covs [ j , i ] = cc  
        
    return RT ( values ,  covs ) 

# =============================================================================
## @brief Calculate optimal batch size for PyTorch iterative boosting.
#  @details Dynamically balances Full-Batch stability for smaller datasets with 
#           SIMD efficiency, RAM limits, and BatchNorm1d stability for larger datasets.
#  @param n_samples  Total number of samples (N_orig + N_targ)
#  @param n_features Number of input features (D)
#  @param max_val    Maximum allowed batch size cap (default: 2^16 = 65536)
#  @return Optimal integer batch size (power of 2)
# =============================================================================
def optimal_batch_size ( n_samples      ,
                         n_features     ,
                         max_batch_size ) :
    """ Calculate optimal batch size dynamically.
    """
    if n_samples <= max_batch_size  : return n_samples

    # Ensure at least 64 samples per feature axis to stabilize BatchNorm statistics
    target_bs = max ( n_samples // 2 , 64 * n_features )
    
    # Smoothly round down to the nearest power of 2
    batch_size = 2 ** math.floor ( math.log2 ( target_bs ) )
    
    return min ( batch_size , max_batch_size )

# =============================================================================
## @brief Compute truncated mean with asymmetric low/high trimming along the specified axis.
#  @param x Input array (e.g. fold predictions matrix of shape (N_folds, N_samples)).
#  @param axis Dimension along which to sort and trim extremes.
#  @param trim_low Number or fraction of elements to trim from the bottom (low end).
#  @param trim_high Number or fraction of elements to trim from the top (high end). 
#                 If None, defaults to `trim_low` (symmetric trimming).
#  @return Aggregated array with trimmed mean.
def truncated_mean ( x         , *    ,
                     axis      = 0    ,
                     trim_low  = 1    ,
                     trim_high = None ) :
    """ Compute asymmetric truncated mean along the specified axis. """
    
    x = numpy.asarray ( x, dtype = numpy.float32 )
    n_folds = x.shape [ axis ]
    
    # Fallback to standard mean if there are too few elements to trim safely
    if n_folds < 3 : return numpy.mean ( x, axis = axis )
    
    # If trim_high is not specified, make it symmetric to trim_low
    if trim_high is None : trim_high = trim_low
        
    # Sort predictions along the fold axis
    sorted_x = numpy.sort ( x, axis = axis )
    
    # Determine number of elements to trim from the bottom (k_low)
    k_low  = int ( trim_low  )
    k_high = int ( trim_high )
        
    # Safety bounds: ensure we don't trim everything or overlap extremes
    k_low  = max ( 0 , k_low   )
    k_high = max ( 0 , k_high  )
    if k_low + k_high >= n_folds : k_low = k_high = 0
    
    if k_low > 0 or k_high > 0 :
        # Construct a dynamic slice tuple to support arbitrary dimensions safely
        slices   = [ slice ( None ) ] * x.ndim
        end_idx  = n_folds - k_high if k_high > 0 else None
        slices [ axis ] = slice ( k_low, end_idx )
        trimmed = sorted_x [ tuple ( slices ) ]
    else :
        trimmed = sorted_x
        
    return numpy.mean ( trimmed, axis = axis )
    
# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        
# =============================================================================
##                                                                      The END 
# =============================================================================
