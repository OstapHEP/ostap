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
    'structured2unstructured' , ## convert structured numpy-array to unstructured one
    ##
    'TwoDataSets'           , ## helper keeper of two dataset with optional weights
    'ComparisonResult'      , ## helper keeper of Comparison results for two datasets 
    'HistoComparisonResult' , ## helper keeper of Comparison results for two histograms                                       
    'prepare_data'          , ## Prepare two datasets for comparison    
    ##
    'clip_pvalue'         , ## clip-value
    'format_pvalue'       , ## format pvalue as nice value in [%]
    'format_nsigma'       , ## format #sigma as nice value in [#sigma]
    ##
    'full_header'         , ## the header for statistical table (full) 
    'short_header'        , ## the header for statistical table (short) 
    't_header'            , ## the header for statistical table (only t-values)
    'p_header'            , ## the header for statistical table (only p-values&nsigma
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
    'symbol_AB'           , ## Symbol for Ansari-Bradley
    'symbol_ST'           , ## Symbol for Siegel-Tukey
    'symbol_MD'           , ## Symbol for Mood
    'symbol_VDW'          , ## Symbol for Van der Waerden 
    'symbol_BM'           , ## Symbol for Brunner-Munzel 
    'symbol_ES'           , ## Symbol for Epps-Singleton
    'symbol_RA'           , ## Symbol for Random Axes/Cramer-Wold 
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
    'method_S'            , ##                       Savage method 
    'method_AB'           , ##               Ansari-Bradley method 
    'method_ST'           , ##                 Siegel-Tukey method 
    'method_MD'           , ##                         Mood method 
    'method_VDW'          , ##              Van der Waerden method
    'method_BM'           , ##               Brunner-Munzel method
    'method_ES'           , ##               Epps-Singleton method
    'method_RA'           , ##      Random Axes/Cramer-Wold method 
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
    ## 
) 
# =============================================================================
from   collections            import namedtuple 
from   ostap.core.ostap_types import num_types, numpy_buffer_types, sized_types 
from   ostap.utils.core       import typename
from   ostap.stats.counters   import ECDF, WECDF 
from   ostap.math.math_base   import data2vct, isfinite , FIRST_ENTRY , LAST_ENTRY
from   ostap.math.ve          import VE
import ostap.logger.symbols   as     S 
import numpy, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
logger = getLogger ( 'ostap.stats.utils' )
# ============================================================================
symbol_name  = S.bolditalic 
method_name  = S.italic

symbol_EPOCH = f'{S.repeat} :' if S.show and S.repeat else 'Epoch:'

symbol_LGBM   = S.light_bulb   if S.show and S.light_bulb     else symbol_name ( 'LightGBM' ) 
symbol_XGB    = S.rocket       if S.show and S.rocket         else symbol_name ( 'XGBoost'  ) 
symbol_CATB   = S.cat_face     if S.show and S.cat_face       else symbol_name ( 'CatBoost' ) 
symbol_HGBC   = S.bar_chart    if S.show and S.bar_chart      else symbol_name ( 'HGBC'     ) 
symbol_GBC    = S.nabla        if S.show and S.nabla          else symbol_name ( 'GBC'      ) 
symbol_RF     = S.palm         if S.show and S.palm           else symbol_name ( 'RF'       )
symbol_TORCH  = S.flashlight   if S.show and S.flashlight     else symbol_name ( 'TORCH'    ) 
symbol_KERAS  = S.postal_horn  if S.show and S.postal_horn    else symbol_name ( 'KERAS'    ) 
symbol_LOGREG = S.ruler        if S.show and S.ruler          else symbol_name ( 'LOGREG'   ) 
symbol_USTAT  = S.pushpin      if S.show and S.pushpin        else symbol_name ( 'USTAT'    ) 
symbol_CHI2   = S.chi2         if S.show and S.chi2           else symbol_name ( 'Chi2'     ) 

symbol_MIX   = S.shuffle       if S.show and S.shuffle        else symbol_name ( 'MIX'      ) 
symbol_PPD   = S.ruler         if S.show and S.ruler          else symbol_name ( 'PPD'      ) 
symbol_DNN   = S.compass       if S.show and S.compass        else symbol_name ( 'DNN'      ) 

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
symbol_AB    = symbol_name ( 'A-B' )
symbol_ST    = symbol_name ( 'S-T' )
symbol_MD    = symbol_name ( 'Md'  )
symbol_VDW   = symbol_name ( 'VdW' )
symbol_BM    = symbol_name ( 'B-M' )
symbol_ES    = symbol_name ( 'E-S' )
symbol_RA    = S.globe if S.show and S.globe else symbol_name ( 'C-W' )   


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
method_AB    = method_name ( 'Ansari-Bradley'                )
method_ST    = method_name ( 'Siegel-Tukey'                  )
method_MD    = method_name ( 'Mood'                          )
method_VDW   = method_name ( 'Van der Waerden'               )
method_BM    = method_name ( 'Brunner-Munzel'                )
method_ES    = method_name ( 'Epps-Singleton'                )
method_RA    = method_name ( 'Cramer-Wold'                   )

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
p_header        = ( S.p_value_percent_expr ,
                    S.n_sigma_expr         )
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
    w_arr = numpy.asarray ( W , dtype = numpy.float64 )
    
    # Step 4: Compute the sum of weights (Sum of W)
    sum_w = float ( numpy.sum ( w_arr ) ) 
    
    # Step 5: Compute the sum of squared weights (Sum of W^2) 
    # using an optimized dot product for performance.
    sum_w2 = float ( numpy.dot ( w_arr, w_arr ) ) 

    # Step 6: Apply Kish's formula: N_eff = (sum_w)^2 / sum_w2.
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
def valid_weights_shape ( weights ) :
    """ Check if weights have a valid 1D array or (N, 1) column shape.
    - Supports NumPy arrays, lists, tuples, and custom containers.
    """
    if weight_trivial ( weights ) : return True

    shape = getattr ( weights , 'shape' , None )

    if shape is not None:
        ndim = len ( shape )
        if ndim <= 1 : return True
        if ndim == 2 and shape[1] == 1: return True
        return False
    
    if hasattr ( weights , '__len__' ) and 0 < len  ( weights ) :
        first_element = weights [ 0 ]
        if isinstance ( first_element, sized_types ) :
            return 1 == len ( first_element )

    return True

# =============================================================================
## Check if dataset has a valid structure/shape
def valid_data_shape(data):
    """ Check if data has a valid 1D, 2D, or Structured Array shape.
    """
    if data is None: return False

    shape = getattr ( data , 'shape' , None )
    if shape is not None:
        ndim = len ( shape )
        return True if 1 <= ndim <= 2 else False 

    return True if hasattr ( data , '__len__' ) else False 

# =============================================================================
## Number of features for training data
def num_features ( X ) :
    """ Number of features for training data
    """
    dtype = getattr ( X , 'dtype' , None )
    if dtype is not None and dtype.names : return len ( dtype.names )

    shape = getattr ( X , 'shape' , None )
    if shape is not None and 1 < len ( shape ) : return shape [ 1 ]

    return 1

# =============================================================================
## Returns the number of samples/events (rows) in the dataset.
def num_samples(X):
    """ Number of samples/events (rows) in the dataset.
    """
    if X is None : return 0

    if hasattr ( X, "__len__" ) : return len ( X )

    shape = getattr ( X , "shape" , None )
    if shape is not None:
        if 0 == len ( shape ) : return 0 
        return int ( shape  [ 0 ] )

    raise TypeError ( "Unsupported data structure for determining sample count: %s" % typename ( X ) )

# =============================================================================
## Check if input dataset shapes have compatible feature dimensions
def compatible_shapes ( ds1 , ds2 ) :
    """ Check if ds1 and ds2 have compatible feature dimensions"""
    dt1 = getattr ( ds1 , 'dtype' , None )
    dt2 = getattr ( ds2 , 'dtype' , None )

    if dt1 is not None and dt1.names and dt2 is not None and dt2.names:
        return dt1.names == dt2.names

    s1 = numpy.shape ( ds1 )
    s2 = numpy.shape ( ds2 )

    if not ( s1 and s2 and s1 [ 0 ] > 0 and s2 [ 0 ] > 0 ):
        return False

    dim1 = s1 [ 1 : ] if 1 < len ( s1 ) else ()
    dim2 = s2 [ 1 : ] if 1 < len ( s2 ) else ()

    return dim1 == dim2

# =============================================================================
## Check if dataset and weights are compatible
def compatible_weights ( data , weights = None ) :
    """Check for data and weights compatibility.
    """
    if weight_trivial ( weights ) : return True

    if not valid_data_shape    ( data    ) : return False
    if not valid_weights_shape ( weights ) : return False

    return num_samples ( data ) == num_samples ( weights )

# ===========================================================================
## Valid weight?
def valid_weight ( weight ) :
    """ Valid weight?
    """
    if     weight_trivial      ( weight ) : return True
    if not valid_weights_shape ( weight ) : return False
    
    try : 
        return 0 < float ( numpy.sum ( weight ) )
    except Exception : 
        return False
    
# =============================================================================
## transform structured array to unstructured one
from numpy.lib.recfunctions import structured_to_unstructured as s2u
# =============================================================================
## Convert structured array to unstructured one
def structured2unstructured ( data, copy = False , dtype = numpy.float64  ) :
    """ Convert structured array to unstructured one
    """
    if not data.dtype.names : return data
    return s2u ( data , copy = copy , dtype = dtype ) 
    
# =============================================================================
from ostap.math.math_ve import gauss_cdf 
one_sigma_minus = gauss_cdf ( -1.0 ) 
one_sigma_plus  = gauss_cdf ( +1.0 )  

# =============================================================================
## Internal implementation for single or pooled normalization 
def _normalize_impl ( data1          ,
                      data2   = None ,
                      weight1 = None ,
                      weight2 = None ) :
    """ Internal implementation for single or pooled normalization """
    
    w1_trivial = weight_trivial ( weight1 )
    w2_trivial = weight_trivial ( weight2 ) if data2 is not None else True

    n1   = num_samples ( data1 )
    
    # ========================================================================
    ## Case with a single dataset
    # ========================================================================
    if data2 is None :

                
        if w1_trivial :
            
            mu  = numpy.mean ( data1 , axis = 0 , keepdims = True )
            var = numpy.var  ( data1 , axis = 0 , keepdims = True )  
            
            var   = numpy.where ( var > 0 , var , 1.0 )
            scale = numpy.sqrt  ( var )
            
            return mu, scale
        
        else :
            
            w1    = numpy.asarray ( weight1 , dtype = numpy.float64 ).ravel() if weight1 is not None else numpy.ones ( n1 , dtype = numpy.float64 )
            w1    = numpy.ascontiguousarray ( w1 )
            sum_w = numpy.sum ( w1 )
            
            # 1. Centering: weighted mean (1D weight @ 2D data)
            mean = ( w1 @ data1 ) / sum_w
            if mean.ndim == 1:
                mean = mean.reshape(1, -1)
            
            x_centered = data1 - mean
            
            # 2. Weighted variance O(N)
            var = ( w1 @ ( x_centered * x_centered ) ) / sum_w
            if var.ndim == 1:
                var = var.reshape(1, -1)

            if numpy.all ( 0 <= w1 ) or numpy.all ( 0 < var ) :
                var   = numpy.where ( var > 0 , var , 1.0 )
                scale = numpy.sqrt  ( var )
                return mean, scale
            else :
                scale = numpy.empty_like ( mean )
                pos_mask = var > 0
                scale [ pos_mask ] = numpy.sqrt ( var [ pos_mask ] )
                
                if not numpy.all ( pos_mask ) :
                    n_features = num_features ( data1 )
                    for j in range ( n_features ) :
                        if not pos_mask [ 0 , j ] :
                            col_data    = data1 [ : , j ] if 1 < data1.ndim else data1
                            col_weights = w1
                            
                            ecdf = data2ecdf ( col_data , col_weights )
                            q1   = ecdf.quantile ( one_sigma_minus )
                            q2   = ecdf.quantile ( one_sigma_plus  )
                            
                            robust_sigma = ( q2 - q1 ) * 0.5
                            scale [ 0 , j ] = robust_sigma
                            
                scale = numpy.where ( scale > 0 , scale , 1.0 )
                return mean, scale

    # ========================================================================
    ## Case with two datasets (pooling)
    # ========================================================================
    n2   = num_samples ( data2 )
    ntot = n1 + n2

    data1 = structured2unstructured ( data1 , copy = False ) 
    data2 = structured2unstructured ( data2 , copy = False ) 
    
    if w1_trivial and w2_trivial :
        mu1  = numpy.mean ( data1 , axis = 0 , keepdims = True )
        mu2  = numpy.mean ( data2 , axis = 0 , keepdims = True )
        
        var1 = numpy.var  ( data1 , axis = 0 , keepdims = True )  
        var2 = numpy.var  ( data2 , axis = 0 , keepdims = True ) 

        pooled_mean = ( n1 * mu1  + n2 * mu2  ) / ntot 
        
        diff1 = mu1 - pooled_mean
        diff2 = mu2 - pooled_mean
        
        pooled_var  = ( n1 * ( var1 + diff1 * diff1 ) + 
                        n2 * ( var2 + diff2 * diff2 ) ) / ntot 
        
        pooled_var   = numpy.where ( pooled_var > 0 , pooled_var , 1.0 )
        pooled_scale = numpy.sqrt  ( pooled_var )

        return pooled_mean, pooled_scale

    # ========================================================================
    ## Case with weights (two datasets) - 1D arrays  + оператор @
    # ========================================================================
    w1 = numpy.asarray ( weight1 , dtype = numpy.float64 ).ravel() if weight1 is not None else numpy.ones ( n1 , dtype = numpy.float64 )
    w2 = numpy.asarray ( weight2 , dtype = numpy.float64 ).ravel() if weight2 is not None else numpy.ones ( n2 , dtype = numpy.float64 )

    w1 = numpy.ascontiguousarray ( w1 )
    w2 = numpy.ascontiguousarray ( w2 )

    x_pool = numpy.concatenate ( ( data1 , data2 ) , axis = 0 )
    w_pool = numpy.concatenate ( ( w1    , w2    ) , axis = 0 )
    
    sum_w  = numpy.sum ( w_pool )
    
    # 1. Centering: weighted mean of the pool (1D w_pool @ 2D x_pool)
    pooled_mean = ( w_pool @ x_pool ) / sum_w
    if pooled_mean.ndim == 1:
        pooled_mean = pooled_mean.reshape(1, -1)
        
    x_centered  = x_pool - pooled_mean
    
    # 2. Weighted variance O(N)
    pooled_var = ( w_pool @ ( x_centered * x_centered ) ) / sum_w
    if pooled_var.ndim == 1:
        pooled_var = pooled_var.reshape(1, -1)

    # ========================================================================
    ## Case with all weights non-negative or all pooled var non-negative 
    # ========================================================================     
    if numpy.all ( 0 <= w_pool ) or numpy.all ( 0 < pooled_var ) :
        pooled_var   = numpy.where ( pooled_var > 0 , pooled_var , 1.0 )
        pooled_scale = numpy.sqrt  ( pooled_var )
        return pooled_mean, pooled_scale
    
    # ========================================================================
    ## General case: negative weights and some pooled variances are <= 0
    # ========================================================================
    pooled_scale = numpy.empty_like ( pooled_mean )
    
    pos_mask     = pooled_var > 0
    pooled_scale [ pos_mask ] = numpy.sqrt ( pooled_var [ pos_mask ] )
    
    if not numpy.all ( pos_mask ) :
        n_features = num_features ( x_pool )
        
        for j in range ( n_features ) :
            if not pos_mask [ 0 , j ] :
                col_data   = x_pool [ : , j ] if 1 < x_pool.ndim else x_pool
                col_weight = w_pool
                
                ecdf = data2ecdf ( col_data , col_weight )
                q1   = ecdf.quantile ( one_sigma_minus )
                q2   = ecdf.quantile ( one_sigma_plus  )
                
                robust_sigma = ( q2 - q1 ) * 0.5
                pooled_scale [ 0 , j ] = robust_sigma 
                
    pooled_scale = numpy.where ( pooled_scale > 0 , pooled_scale , 1.0 )

    return pooled_mean, pooled_scale


# =============================================================================
## Normalize ("standardize") two input datasets such that pooled dataset has
#  zero mean and unit variance
# =============================================================================
def normalize_pooled ( data1         ,
                       data2         ,
                       weight1 = None,
                       weight2 = None ) :
    """ Normalize ("Standardize") two input datasets such that pooled
    dataset has zero means and unit variances for each variable
    
    >>> data1, data2, weight1, weight2 = ...
    >>> ds1 , ds2 = normalize_pooled ( data1 , data2 , weight1 , weight2 )
    """
    ##
    
    data1 = structured2unstructured ( data1 , copy = False ) 
    data2 = structured2unstructured ( data2 , copy = False ) 
    
    check_all ( data1   = data1              ,
                data2   = data2              ,
                weight1 = weight1            ,
                weight2 = weight2            ,
                where   = "normalize_pooled" )

    mean, scale = _normalize_impl ( data1   = data1   ,
                                    data2   = data2   ,
                                    weight1 = weight1 ,
                                    weight2 = weight2 )

    ds1 = ( data1 - mean ) / scale
    ds2 = ( data2 - mean ) / scale
    
    return ds1 , ds2


# =============================================================================
## Normalize ("standardize") single input dataset
# =============================================================================
def normalize_single ( data       ,
                       weight = None ) :
    """ Normalize ("Standardize") single input dataset such that
    dataset has zero mean and unit variance for each variable
    
    >>> data, weight = ...
    >>> ds = normalize ( data, weight )
    """     
    ##
    data = structured2unstructured ( data , copy = False ) 
    ##    
    mean, scale = _normalize_impl ( data1   = data   ,
                                    weight1 = weight )
    return ( data - mean ) / scale

# ============================================================================
## check data and weights
def check_all ( data1   ,
                data2   ,
                weight1 = None        ,
                weight2 = None        ,
                where   = "check_all" ) :
    
    if not valid_data_shape   ( data1           ) : raise TypeError ( "%s: invalid `data1`"   % where ) 
    if not valid_data_shape   ( data2           ) : raise TypeError ( "%s: invalid `data2`"   % where ) 
    if not compatible_shapes  ( data1 , data2   ) : raise TypeError ( "%s: incompatible data shape`"  % where ) 
    if not valid_weight       ( weight1         ) : raise TypeError ( "%s: invalid `weight1`" % where ) 
    if not valid_weight       ( weight2         ) : raise TypeError ( "%s: invalid `weight2`" % where ) 
    if not compatible_weights ( data1 , weight1 ) : raise TypeError ( "%s: incompatible `data1/weight1`"  % where ) 
    if not compatible_weights ( data2 , weight2 ) : raise TypeError ( "%s: incompatible `data2/weight2`"  % where ) 
    return True 

# ==============================================================================
## Convert numpy-array statistics into Ostap::SVectorWithError
def np2vct ( data , weight = None ) :
    """ Convert numpy-array statistics into `Ostap.Math.SVectorWithError`
    """

    shape     = data.shape
    n , N     = shape
    
    w         = None if weight_trivial ( weight ) else weight

    from statsmodels.stats.weightstats import DescrStatsW as DSW 
    dsw       = DSW  ( data , weights = w )
    mean      = dsw.mean
    covmtrx   = dsw.cov
    
    import ostap.math.linalg
    from   ostap.core.core   import Ostap
    
    MN        = Ostap.Math.SymMatrix ( N )
    
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
## clip p-value 
def clip_pvalue ( pvalue , clip = 0.5 ) :
    """ Clip p-value
    """
    if 0 < float ( pvalue ) < 1 : return pvalue 
    pv   = VE ( pvalue )
    cov2 = max (  0 , pv.cov2() ) 
    clip = min ( 0.5 , abs ( clip ) * pv.error () ) if 0 < cov2 else min ( clip , 1e-6 )  
    value = pv.value () 
    if   1 <= value : pv = VE ( 1 - clip , cov2 )
    elif 0 >= value : pv = VE (     clip , cov2 )
    return pv 

# ==============================================================================
## format p-value [in %]  as nice string  
def format_pvalue ( pvalue , clip = 0.5 ) :
    """ Format p-value [in %]  as nice string
    """
    pv = clip_pvalue ( pvalue , clip = clip ) 
    pv = VE ( pv )
    if pv.cov2 () <= 0 : 
        pv = 100 *  float ( pv )
        return '%.2f' % pv 
     
    pv = str ( ( 100 * pv ) .toString ( '%% 5.2f %s %%-.2f' % S.plus_minus ) )
    return pv 

# =============================================================================
## format signifcance [#sigma] as nice string 
def format_nsigma ( nsigma ) : 
    """ Format significance [#sigma] as nce string 
    """
    ns = VE  ( nsigma )
    if ns.cov2 () <= 0 : 
        ns = float ( ns )
        return ( '%.2f' % ns ) if isfinite ( ns ) and ns < 1000 else  S.infinity_pos 
    
    return str ( ns.toString ( '%%.2f %s %%-.2f' % S.plus_minus ) if float ( ns ) < 1000 else S.infinity_pos  ) 

# =============================================================================
## Calculate optimal batch size for PyTorch iterative boosting.
def optimal_batch_size ( n_samples      ,
                         n_features     ,
                         max_batch_size ) :
    """ Calculate optimal batch size dynamically.
    """
    if n_samples <= max_batch_size  : return n_samples

    target_bs = max ( n_samples // 2 , 64 * n_features )
    batch_size = 2 ** math.floor ( math.log2 ( target_bs ) )
    
    return min ( batch_size , max_batch_size )

# =============================================================================
## Compute truncated mean with asymmetric low/high trimming along the specified axis.
def truncated_mean ( x         , *    ,
                     axis      = 0    ,
                     trim_low  = 1    ,
                     trim_high = None ) :
    """ Compute asymmetric truncated mean along the specified axis. """
    
    x = numpy.asarray ( x, dtype = numpy.float32 )
    n_folds = x.shape [ axis ]
    
    if n_folds < 3 : return numpy.mean ( x, axis = axis )
    
    if trim_high is None : trim_high = trim_low
        
    sorted_x = numpy.sort ( x, axis = axis )
    
    k_low  = int ( trim_low  )
    k_high = int ( trim_high )
        
    k_low  = max ( 0 , k_low   )
    k_high = max ( 0 , k_high  )
    if k_low + k_high >= n_folds : k_low = k_high = 0
    
    if k_low > 0 or k_high > 0 :
        slices   = [ slice ( None ) ] * x.ndim
        end_idx  = n_folds - k_high if k_high > 0 else None
        slices [ axis ] = slice ( k_low, end_idx )
        trimmed = sorted_x [ tuple ( slices ) ]
    else :
        trimmed = sorted_x
        
    return numpy.mean ( trimmed, axis = axis )
    

# ============================================================================
ECDF_types = ECDF , WECDF
# ============================================================================
## Convert data/weight pair into (W)ECDF
def data2ecdf ( data , weight = None ) :
    """ Convert data/weight pair into (W)ECDF
    """
    if isinstance ( data , ECDF_types ) : return data
    
    data  = numpy.asarray   ( data   , dtype = numpy.float64 ).ravel()
    data  = numpy.ascontiguousarray ( data )
    
    if 1 != num_features        ( data   )          : raise TypeError  ( "Invalid `data` type!"                   )
    
    if       weight_trivial     ( weight )          : return ECDF ( data2vct ( data ) )
    elif not valid_weight       ( weight )          : raise ValueError ( "Invalid `weight` type/value!"           )
    elif not compatible_weights ( data   , weight ) : raise TypeError  ( "Incompatible `data/weight` structures!" )
    
    weight = numpy.asarray ( weight , dtype = numpy.float64 ).ravel()
    weight = numpy.ascontiguousarray ( weight )
    
    return WECDF ( data2vct ( data ) , data2vct ( weight ) ) 

# =============================================================================
## comparison result  
ComparisonResult     = namedtuple ( 'ComparisonResult'  ,
                                    ( 'symbol'     ,
                                      'method'     ,
                                      'tvalue'     ,
                                      'pvalue'     ,
                                      'nsigma'     ,
                                      'importance' ,
                                      'table_row'  ) ,
                                    defaults = ( '' , '<UNKNOWN>' , -1e+9 , -1 , -1 , {} , '' ) )
# =============================================================================
## 1D-histogram comparison result  
HistoComparisonResult = namedtuple ( 'HistoComparisonResult'  ,
                                     ( 'histo1'    ,
                                       'histo2'    ,
                                       'chi2'      ,
                                       'nDoF'      , 
                                       'pvalue'    ,
                                       'nsigma'    ) )

# =============================================================================
## helper keeper of two datasets
#  - convert to dataclass(frozen=True) ??
TwoDataSets = namedtuple ( 'TwoDataSets' ,
                           ( 'data1'   ,
                             'data2'   ,
                             'weight1' ,
                             'weight2' ) ,
                           defaults = ( None , None ) )

# =============================================================================
## Prepare two datasets for comparison:
#  Each dataset is converted to *STRUCTURED* numpy array 
#  @see numpy_compare 
#  @code
#  datasets = prepare_data ( data1       = data1    ,
#                            data2       = data2    ,
#                            expressions = 'X,y,z'  ,  
#                            cuts        =   'pt>1' )
#  @endcode 
def prepare_data ( data         ,
                   data2        ,
                   expressions  ,                   ## variables in data1 
                   cuts         = ''          , * , ## cuts for data1
                   expressions2 = None        ,
                   cuts2        = ''          ,                            
                   first        = FIRST_ENTRY ,
                   last         =  LAST_ENTRY ,                                                                      
                   first2       = FIRST_ENTRY ,
                   last2        =  LAST_ENTRY ,                                                                      
                   cut_range    = ''          ,
                   cut_range2   = None        , 
                   prescale     = None        , 
                   prescale2    = None        , 
                   use_frame    = False       , 
                   use_frame2   = None        , 
                   parallel     = False       , 
                   parallel2    = None        , 
                   progress     = False       ,
                   progress2    = None        ) :
    """ Prepare two datasets for comparison:
    
    - Each dataset is converted to *STRUCTURED* numpy array 
    
    >>> data1, data2 , weight1 , weight2 = data_prepare ( data1       = data1    ,
    ...                                                   data2       = data2    ,
    ...                                                   expressinos = 'X,y,z'  ,  
    ...                                                   cuts        =   'pt>1' )
    """
    if expressions2 is None : expressions2 = expressions
    if cuts2        is None : cuts2        = cuts
    
    if first2       is None : first2       = first
    if last2        is None : last2        = last
    
    if cut_range2   is None : cut_range2   = cut_range
    if use_frame2   is None : use_frame2   = use_frame
    if parallel2    is None : parallel2    = parallel
    if progress2    is None : progress2    = progress
    
    from   ostap.trees.cuts     import vars_and_cuts
    from   ostap.stats.statvars import data_slice
    
    varlst1 , cuts  , _  = vars_and_cuts ( expressions  , cuts )
    varlst1 = tuple ( varlst1 ) 
    
    varlst2 , cuts2 , _  = vars_and_cuts ( expressions2 , cuts2 )
    varlst2 = tuple ( varlst2 ) 
    
    assert varlst1 and len ( varlst1 ) == len ( varlst2 ) , \
        "Different lenghts for variable lists!"

    ## (1) create numpy datasets 
    npdata1 , weight1 = data_slice ( data                   ,
                                     varlst1                , 
                                     cuts       = cuts      ,
                                     first      = first     ,
                                     last       = last      ,
                                     cut_range  = cut_range ,
                                     prescale   = prescale  , 
                                     structured = True      , ## ATTENTION!  
                                     transpose  = True      ,
                                     progress   = progress  ,
                                     use_frame  = use_frame , 
                                     parallel   = parallel  )
    
    npdata2 , weight2 = data_slice ( data2                   ,
                                     varlst2                 , 
                                     cuts       = cuts2      ,
                                     first      = first2     ,
                                     last       = last2      ,
                                     cut_range  = cut_range2 ,
                                     prescale   = prescale2  , 
                                     structured = True       , ## ATTENTION!  
                                     transpose  = True       ,
                                     progress   = progress2  ,
                                     use_frame  = use_frame2 , 
                                     parallel   = parallel2  )

    ## 
    return TwoDataSets ( data1   = npdata1 ,
                         data2   = npdata2 ,
                         weight1 = weight1 , 
                         weight2 = weight2 ) 

# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        
# =============================================================================
##                                                                      The END 
# =============================================================================



