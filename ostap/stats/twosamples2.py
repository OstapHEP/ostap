#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/twosamples2.py
#  Two-Sample Weighted Tests
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-10-16
# =============================================================================
""" Two-Sample Weighted Tests
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-16"
__all__     = (
    'kolmogorov_smirnov' , ## Kolmogorov-Smirnov 2-samples test 
    'kuiper'             , ## Kuiper             2-samples test 
    'cramer_von_mises'   , ## Cramer-von Mises   2-sample  test 
    'anderson_darling'   , ## Anderson-Darling   2-sample  test 
    'berk_jones'         , ## Berk-Jones         2-sample  test 
    'ZA'                 , ## Zhang's ZA         2-sample  test 
    'ZC'                 , ## Zhang's ZC         2-sample  test 
    'ZK'                 , ## Zhang's ZK         2-sample  test
) 
# =============================================================================
from   ostap.stats.utils import weight_trivial 
import numpy 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.twosamples2' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Two-sample Weighted test' )
# =============================================================================
## Calculates the two-sample Kolmogorov-Smirnov D-statistic for 1D datasets
#  with weights.
# 
#  Supports positive weights and sPlot weights. Returns raw D-statistic
#  (maximum absolute difference between empirical CDFs)
# 
#   To obtain a valid p-value for sPlot data (with negative weights), use the
#   returned 'd_stat' with an sPlot-safe Bootstrap.
def kolmogorov_smirnov ( data1    ,
                         data2    ,
                         weight1 = None ,
                         weight2 = None ):
    """ Calculates the two-sample Kolmogorov-Smirnov D-statistic for 1D datasets
    with weights.

    Supports positive weights and sPlot weights. Returns raw D-statistic
    (maximum absolute difference between empirical CDFs)

    To obtain a valid p-value for sPlot data (with negative weights), use the
    returned 'd_stat' with an sPlot-safe Bootstrap.
    """
    x1 = numpy.asarray ( data1 , dtype = float ).ravel() 
    x2 = numpy.asarray ( data2 , dtype = float ).ravel() 

    w1_trivial = weight_trivial ( weight1 )
    w2_trivial = weight_trivial ( weight2 )
    
    w1 = numpy.asarray ( weight1 , dtype = float ) if not w1_trivial else numpy.ones_like ( x1 )
    w2 = numpy.asarray ( weight2 , dtype = float ) if not w2_trivial else numpy.ones_like ( x2 )

    N1 = len ( x1 )
    N2 = len ( x2 )
    if N1 < 2 or N2 < 2: return  0

    # 1. Sort data and corresponding weights
    idx1 = numpy.argsort ( x1 )
    idx2 = numpy.argsort ( x2 )

    x1, w1 = x1 [ idx1 ] , w1 [ idx1 ]
    x2, w2 = x2 [ idx2 ] , w2 [ idx2 ]

    # 2. Compute normalized cumulative sum of weights (weighted ECDFs)
    sum_w1 = numpy.sum ( w1 )
    sum_w2 = numpy.sum ( w2 )

    cdf1 = numpy.cumsum ( w1 ) / sum_w1
    cdf2 = numpy.cumsum ( w2 ) / sum_w2

    # 3. Evaluate both ECDFs on the union of all data points
    all_x = numpy.union1d ( x1, x2 )

    cdf1_all = numpy.searchsorted ( x1 , all_x , side = "right" )
    cdf2_all = numpy.searchsorted ( x2 , all_x , side = "right" )

    cdf1_eval = numpy.pad ( cdf1 , ( 1 , 0 ) ) [ cdf1_all ]
    cdf2_eval = numpy.pad ( cdf2 , ( 1 , 0 ) ) [ cdf2_all ]

    # 4. Compute directional deviations and maximum absolute difference (D-statistic)
    diff      = cdf1_eval - cdf2_eval
    d_plus    = float ( numpy.max ( diff ) )
    d_minus   = float ( numpy.max (-diff ) )
    d_stat    = float ( numpy.max ( numpy.abs ( diff ) ) )

    return d_stat

# =============================================================================
## Calculates the two-sample Kuiper statistic (V-statistic) for 1D datasets
#  with weights.
# 
#  Kuiper's statistic V = D+ + D- measures the sum of maximum deviations
#  above and below the ECDFs. It is invariant under cyclic transformations
#  and equally sensitive at the tails as in the median area.
# 
#  Supports positive weights and sPlot weights.
def kuiper ( data1          ,
             data2          ,
             weight1 = None ,
             weight2 = None ) :
    """ Calculates the two-sample Kuiper statistic for 1D datasets with weights.

    Kuiper's statistic V = D+ + D- measures the sum of maximum deviations
    above and below the ECDFs. It is invariant under cyclic transformations.

    Supports positive weights and sPlot weights.
    """
    x1 = numpy.asarray ( data1 , dtype = float ).ravel()
    x2 = numpy.asarray ( data2 , dtype = float ).ravel() 

    w1_trivial = weight_trivial ( weight1 )
    w2_trivial = weight_trivial ( weight2 )
    
    w1 = numpy.asarray ( weight1 , dtype = float ) if not w1_trivial else numpy.ones_like ( x1 )
    w2 = numpy.asarray ( weight2 , dtype = float ) if not w2_trivial else numpy.ones_like ( x2 )

    N1 = len ( x1 )
    N2 = len ( x2 )
    if N1 < 2 or N2 < 2 : return 0.0

    # 1. Sort data and corresponding weights
    idx1 = numpy.argsort ( x1 )
    idx2 = numpy.argsort ( x2 )

    x1, w1 = x1 [ idx1 ] , w1 [ idx1 ]
    x2, w2 = x2 [ idx2 ] , w2 [ idx2 ]

    # 2. Compute normalized cumulative sum of weights (weighted ECDFs)
    sum_w1 = numpy.sum ( w1 )
    sum_w2 = numpy.sum ( w2 )

    cdf1 = numpy.cumsum ( w1 ) / sum_w1
    cdf2 = numpy.cumsum ( w2 ) / sum_w2

    # 3. Evaluate both ECDFs on the union of all data points
    all_x = numpy.union1d ( x1 , x2 )

    cdf1_all = numpy.searchsorted ( x1 , all_x , side = "right" )
    cdf2_all = numpy.searchsorted ( x2 , all_x , side = "right" )

    cdf1_eval = numpy.pad ( cdf1 , ( 1 , 0 ) ) [ cdf1_all ]
    cdf2_eval = numpy.pad ( cdf2 , ( 1 , 0 ) ) [ cdf2_all ]

    # 4. Compute directional deviations and Kuiper statistic V = D+ + D-
    diff   = cdf1_eval - cdf2_eval
    d_plus = float ( numpy.max ( diff ) )
    d_minus= float ( numpy.max ( -diff ) )
    
    # Kuiper statistic is non-negative
    v_stat = max ( 0.0 , d_plus ) + max ( 0.0 , d_minus )

    return float ( v_stat )

# =============================================================================
## Calculates the two-sample Cramér-von Mises statistic (T / W^2) for 1D datasets
#  with weights.
# 
#  Measures the integrated squared difference between empirical CDFs:
#  T = integral ( F1(x) - F2(x) )^2 d H(x), where H(x) is the pooled CDF.
# 
#  Supports positive weights and sPlot weights.
def cramer_von_mises ( data1          ,
                       data2          ,
                       weight1 = None ,
                       weight2 = None ) :
    """ Calculates the two-sample Cramér-von Mises statistic for 1D datasets
    with weights.

    Measures the integrated squared difference between empirical CDFs weighted
    by the pooled empirical distribution.

    Supports positive weights and sPlot weights.
    """
    x1 = numpy.asarray ( data1 , dtype = float ).ravel()
    x2 = numpy.asarray ( data2 , dtype = float ).ravel() 

    w1_trivial = weight_trivial ( weight1 )
    w2_trivial = weight_trivial ( weight2 )
    
    w1 = numpy.asarray ( weight1 , dtype = float ) if not w1_trivial else numpy.ones_like ( x1 )
    w2 = numpy.asarray ( weight2 , dtype = float ) if not w2_trivial else numpy.ones_like ( x2 )

    N1 = len ( x1 )
    N2 = len ( x2 )
    if N1 < 2 or N2 < 2 : return 0.0

    # 1. Sort data and corresponding weights
    idx1 = numpy.argsort ( x1 )
    idx2 = numpy.argsort ( x2 )

    x1, w1 = x1 [ idx1 ] , w1 [ idx1 ]
    x2, w2 = x2 [ idx2 ] , w2 [ idx2 ]

    # 2. Compute normalized cumulative sum of weights
    sum_w1 = numpy.sum ( w1 )
    sum_w2 = numpy.sum ( w2 )

    cdf1 = numpy.cumsum ( w1 ) / sum_w1
    cdf2 = numpy.cumsum ( w2 ) / sum_w2

    # 3. Evaluate ECDFs on pooled data points
    all_x = numpy.union1d ( x1 , x2 )

    cdf1_all = numpy.searchsorted ( x1 , all_x , side = "right" )
    cdf2_all = numpy.searchsorted ( x2 , all_x , side = "right" )

    cdf1_eval = numpy.pad ( cdf1 , ( 1 , 0 ) ) [ cdf1_all ]
    cdf2_eval = numpy.pad ( cdf2 , ( 1 , 0 ) ) [ cdf2_all ]

    # 4. Compute pooled weights dH(x) for integration
    # Find weights associated with pooled unique points
    w1_pooled = numpy.bincount ( numpy.searchsorted ( all_x , x1 ) , weights = w1 , minlength = len ( all_x ) )
    w2_pooled = numpy.bincount ( numpy.searchsorted ( all_x , x2 ) , weights = w2 , minlength = len ( all_x ) )
    
    dH = ( w1_pooled + w2_pooled ) / ( sum_w1 + sum_w2 )

    # 5. Compute integrated squared deviation integral (F1 - F2)^2 dH
    diff = cdf1_eval - cdf2_eval
    cvm_stat = numpy.sum ( diff * diff * dH )

    return float ( cvm_stat )

# =============================================================================
## Calculates the two-sample Anderson-Darling statistic (A^2) for 1D datasets
#  with weights.
# 
#  Measures the integrated squared difference weighted by H(x)*(1 - H(x)),
#  giving much higher weight to the tails of the distribution.
# 
#  Supports positive weights and sPlot weights.
def anderson_darling ( data1          ,
                       data2          ,
                       weight1 = None ,
                       weight2 = None ,
                       eps      = 1e-10 ) :
    """ Calculates the two-sample Anderson-Darling statistic for 1D datasets
    with weights.

    Weighted integrated squared difference with higher sensitivity to tail discrepancies.

    Supports positive weights and sPlot weights.
    """
    x1 = numpy.asarray ( data1 , dtype = float ).ravel() 
    x2 = numpy.asarray ( data2 , dtype = float ).ravel() 

    w1_trivial = weight_trivial ( weight1 )
    w2_trivial = weight_trivial ( weight2 )
    
    w1 = numpy.asarray ( weight1 , dtype = float ) if not w1_trivial else numpy.ones_like ( x1 )
    w2 = numpy.asarray ( weight2 , dtype = float ) if not w2_trivial else numpy.ones_like ( x2 )

    N1 = len ( x1 )
    N2 = len ( x2 )
    if N1 < 2 or N2 < 2 : return 0.0

    # 1. Sort data and corresponding weights
    idx1 = numpy.argsort ( x1 )
    idx2 = numpy.argsort ( x2 )

    x1, w1 = x1 [ idx1 ] , w1 [ idx1 ]
    x2, w2 = x2 [ idx2 ] , w2 [ idx2 ]

    # 2. Compute normalized cumulative sum of weights
    sum_w1 = numpy.sum ( w1 )
    sum_w2 = numpy.sum ( w2 )

    cdf1 = numpy.cumsum ( w1 ) / sum_w1
    cdf2 = numpy.cumsum ( w2 ) / sum_w2

    # 3. Evaluate ECDFs on pooled data points
    all_x = numpy.union1d ( x1 , x2 )

    cdf1_all = numpy.searchsorted ( x1 , all_x , side = "right" )
    cdf2_all = numpy.searchsorted ( x2 , all_x , side = "right" )

    cdf1_eval = numpy.pad ( cdf1 , ( 1 , 0 ) ) [ cdf1_all ]
    cdf2_eval = numpy.pad ( cdf2 , ( 1 , 0 ) ) [ cdf2_all ]

    # 4. Pooled CDF H(x) and step sizes dH(x)
    w1_pooled = numpy.bincount ( numpy.searchsorted ( all_x , x1 ) , weights = w1 , minlength = len ( all_x ) )
    w2_pooled = numpy.bincount ( numpy.searchsorted ( all_x , x2 ) , weights = w2 , minlength = len ( all_x ) )
    
    dH = ( w1_pooled + w2_pooled ) / ( sum_w1 + sum_w2 )
    H  = numpy.cumsum ( dH )

    # 5. Tail weighting factor 1 / (H*(1-H)) with epsilon clamping
    variance_weight = H * ( 1.0 - H )
    variance_weight = numpy.clip ( variance_weight , eps , 0.25 )

    # 6. Compute integrated weighted squared deviation
    diff = cdf1_eval - cdf2_eval
    ad_stat = numpy.sum ( ( diff * diff / variance_weight ) * dH )

    return float ( ad_stat )

# =============================================================================
## Calculates the two-sample Berk-Jones statistic (R_BJ) for 1D datasets
#  with weights.
# 
#  Measures the supremum of point-wise Kullback-Leibler divergences between
#  empirical CDFs. Exceptionally powerful for sparse tail signals.
# 
#  Supports positive weights and sPlot weights.
def berk_jones ( data1          ,
                 data2          ,
                 weight1 = None ,
                 weight2 = None ,
                 eps      = 1e-12 ) :
    """ Calculates the two-sample Berk-Jones statistic for 1D datasets
    with weights.

    Based on supremum of point-wise Kullback-Leibler divergence.
    Optimal power for rare and tail discrepancies.

    Supports positive weights and sPlot weights.
    """
    x1 = numpy.asarray ( data1 , dtype = float ).ravel()
    x2 = numpy.asarray ( data2 , dtype = float ).ravel()

    w1_trivial = weight_trivial ( weight1 )
    w2_trivial = weight_trivial ( weight2 )
    
    w1 = numpy.asarray ( weight1 , dtype = float ) if not w1_trivial else numpy.ones_like ( x1 )
    w2 = numpy.asarray ( weight2 , dtype = float ) if not w2_trivial else numpy.ones_like ( x2 )

    N1 = len ( x1 )
    N2 = len ( x2 )
    if N1 < 2 or N2 < 2 : return 0.0

    # 1. Sort data and corresponding weights
    idx1 = numpy.argsort ( x1 )
    idx2 = numpy.argsort ( x2 )

    x1, w1 = x1 [ idx1 ] , w1 [ idx1 ]
    x2, w2 = x2 [ idx2 ] , w2 [ idx2 ]

    # 2. Compute normalized cumulative sum of weights
    sum_w1 = numpy.sum ( w1 )
    sum_w2 = numpy.sum ( w2 )

    cdf1 = numpy.cumsum ( w1 ) / sum_w1
    cdf2 = numpy.cumsum ( w2 ) / sum_w2

    # 3. Evaluate ECDFs on pooled points
    all_x = numpy.union1d ( x1 , x2 )

    cdf1_all = numpy.searchsorted ( x1 , all_x , side = "right" )
    cdf2_all = numpy.searchsorted ( x2 , all_x , side = "right" )

    cdf1_eval = numpy.pad ( cdf1 , ( 1 , 0 ) ) [ cdf1_all ]
    cdf2_eval = numpy.pad ( cdf2 , ( 1 , 0 ) ) [ cdf2_all ]

    # Clamp ECDFs to interval [eps, 1 - eps] to prevent log(0)
    F1 = numpy.clip ( cdf1_eval , eps , 1.0 - eps )
    F2 = numpy.clip ( cdf2_eval , eps , 1.0 - eps )

    # 4. Point-wise Binary Kullback-Leibler divergence K(F1, F2)
    # K(p, q) = p * log(p / q) + (1 - p) * log((1 - p) / (1 - q))
    kl_div = F1 * numpy.log ( F1 / F2 ) + ( 1.0 - F1 ) * numpy.log ( ( 1.0 - F1 ) / ( 1.0 - F2 ) )

    # 5. Berk-Jones statistic is the maximum KL divergence
    bj_stat = float ( numpy.max ( kl_div ) )

    return bj_stat

# =============================================================================
## Calculates the two-sample Zhang's Z_A statistic for 1D datasets with weights.
# 
#  Z_A is an Anderson-Darling type log-likelihood ratio statistic.
#  Gives maximum weight to discrepancies in the tails of distributions.
# 
#  Supports positive weights and sPlot weights.
def ZA ( data1            ,
         data2            ,
         weight1 = None  ,
         weight2 = None  ,
         eps      = 1e-12 ) :
    """ Calculates the two-sample Zhang's Z_A statistic for 1D datasets
    with weights.

    Log-likelihood ratio statistic sensitive to tail differences.

    Supports positive weights and sPlot weights.
    """
    x1 = numpy.asarray ( data1 , dtype = float ).ravel()
    x2 = numpy.asarray ( data2 , dtype = float ).ravel() 

    w1_trivial = weight_trivial ( weight1 )
    w2_trivial = weight_trivial ( weight2 )
    
    w1 = numpy.asarray ( weight1 , dtype = float ) if not w1_trivial else numpy.ones_like ( x1 )
    w2 = numpy.asarray ( weight2 , dtype = float ) if not w2_trivial else numpy.ones_like ( x2 )

    N1 = len ( x1 )
    N2 = len ( x2 )
    if N1 < 2 or N2 < 2 : return 0.0

    # 1. Sort data and corresponding weights
    idx1 = numpy.argsort ( x1 )
    idx2 = numpy.argsort ( x2 )

    x1, w1 = x1 [ idx1 ] , w1 [ idx1 ]
    x2, w2 = x2 [ idx2 ] , w2 [ idx2 ]

    # 2. Compute normalized cumulative sum of weights
    sum_w1 = numpy.sum ( w1 )
    sum_w2 = numpy.sum ( w2 )

    cdf1 = numpy.cumsum ( w1 ) / sum_w1
    cdf2 = numpy.cumsum ( w2 ) / sum_w2

    # 3. Evaluate ECDFs on pooled points
    all_x = numpy.union1d ( x1 , x2 )

    cdf1_all = numpy.searchsorted ( x1 , all_x , side = "right" )
    cdf2_all = numpy.searchsorted ( x2 , all_x , side = "right" )

    cdf1_eval = numpy.pad ( cdf1 , ( 1 , 0 ) ) [ cdf1_all ]
    cdf2_eval = numpy.pad ( cdf2 , ( 1 , 0 ) ) [ cdf2_all ]

    # 4. Pooled weights dH(x) and Pooled CDF H(x)
    w1_pooled = numpy.bincount ( numpy.searchsorted ( all_x , x1 ) , weights = w1 , minlength = len ( all_x ) )
    w2_pooled = numpy.bincount ( numpy.searchsorted ( all_x , x2 ) , weights = w2 , minlength = len ( all_x ) )
    
    dH = ( w1_pooled + w2_pooled ) / ( sum_w1 + sum_w2 )

    # Clamp ECDFs to interval [eps, 1 - eps]
    F1 = numpy.clip ( cdf1_eval , eps , 1.0 - eps )
    F2 = numpy.clip ( cdf2_eval , eps , 1.0 - eps )

    # 5. Point-wise Zhang log-likelihood terms for Z_A
    # Sum of log(F1) + log(1 - F1) weighted by pooled step measure dH
    ll_terms = numpy.log ( F1 ) + numpy.log ( 1.0 - F1 ) + numpy.log ( F2 ) + numpy.log ( 1.0 - F2 )
    za_stat  = - numpy.sum ( ll_terms * dH )

    return float ( za_stat )

# =============================================================================
## Calculates the two-sample Zhang's Z_C statistic for 1D datasets with weights.
# 
#  Z_C is a Cramér-von Mises type log-likelihood ratio statistic.
#  Measures integrated divergence evenly across the distribution core.
# 
#  Supports positive weights and sPlot weights.
def ZC ( data1            ,
         data2            ,
         weight1 = None  ,
         weight2 = None  ,
         eps      = 1e-12 ) :
    """ Calculates the two-sample Zhang's Z_C statistic for 1D datasets
    with weights.

    Log-likelihood ratio statistic with uniform weight distribution across support.

    Supports positive weights and sPlot weights.
    """
    x1 = numpy.asarray ( data1 , dtype = float ).ravel()
    x2 = numpy.asarray ( data2 , dtype = float ).ravel() 

    w1_trivial = weight_trivial ( weight1 )
    w2_trivial = weight_trivial ( weight2 )
    
    w1 = numpy.asarray ( weight1 , dtype = float ) if not w1_trivial else numpy.ones_like ( x1 )
    w2 = numpy.asarray ( weight2 , dtype = float ) if not w2_trivial else numpy.ones_like ( x2 )

    N1 = len ( x1 )
    N2 = len ( x2 )
    if N1 < 2 or N2 < 2 : return 0.0

    # 1. Sort data and corresponding weights
    idx1 = numpy.argsort ( x1 )
    idx2 = numpy.argsort ( x2 )

    x1, w1 = x1 [ idx1 ] , w1 [ idx1 ]
    x2, w2 = x2 [ idx2 ] , w2 [ idx2 ]

    # 2. Compute normalized cumulative sum of weights
    sum_w1 = numpy.sum ( w1 )
    sum_w2 = numpy.sum ( w2 )

    cdf1 = numpy.cumsum ( w1 ) / sum_w1
    cdf2 = numpy.cumsum ( w2 ) / sum_w2

    # 3. Evaluate ECDFs on pooled points
    all_x = numpy.union1d ( x1 , x2 )

    cdf1_all = numpy.searchsorted ( x1 , all_x , side = "right" )
    cdf2_all = numpy.searchsorted ( x2 , all_x , side = "right" )

    cdf1_eval = numpy.pad ( cdf1 , ( 1 , 0 ) ) [ cdf1_all ]
    cdf2_eval = numpy.pad ( cdf2 , ( 1 , 0 ) ) [ cdf2_all ]

    # Clamp ECDFs to interval [eps, 1 - eps]
    F1 = numpy.clip ( cdf1_eval , eps , 1.0 - eps )
    F2 = numpy.clip ( cdf2_eval , eps , 1.0 - eps )

    # 4. Discrete step-weighted log-likelihood divergence
    # Sum over pooled points scaled by 1 / N_pooled
    rank_weights = 1.0 / ( numpy.arange ( 1 , len ( all_x ) + 1 , dtype = float ) - 0.5 )
    ll_terms     = ( F1 - F2 ) ** 2 / ( F1 * ( 1.0 - F1 ) + F2 * ( 1.0 - F2 ) )
    zc_stat      = numpy.sum ( ll_terms * rank_weights )

    return float ( zc_stat )

# =============================================================================
## Calculates the two-sample Zhang's Z_K statistic for 1D datasets with weights.
# 
#  Z_K is a Kolmogorov-Smirnov type log-likelihood ratio statistic.
#  Computes the supremum of point-wise log-likelihood ratio distances.
# 
#  Supports positive weights and sPlot weights.
def ZK ( data1            , 
         data2            ,
         weight1 = None  ,
         weight2 = None  ,
         eps      = 1e-12 ) :
    """ Calculates the two-sample Zhang's Z_K statistic for 1D datasets
    with weights.

    Supremum-based log-likelihood ratio statistic.

    Supports positive weights and sPlot weights.
    """
    x1 = numpy.asarray ( data1 , dtype = float ).ravel()
    x2 = numpy.asarray ( data2 , dtype = float ).ravel()

    w1_trivial = weight_trivial ( weight1 )
    w2_trivial = weight_trivial ( weight2 )
    
    w1 = numpy.asarray ( weight1 , dtype = float ) if not w1_trivial else numpy.ones_like ( x1 )
    w2 = numpy.asarray ( weight2 , dtype = float ) if not w2_trivial else numpy.ones_like ( x2 )

    N1 = len ( x1 )
    N2 = len ( x2 )
    if N1 < 2 or N2 < 2 : return 0.0

    # 1. Sort data and corresponding weights
    idx1 = numpy.argsort ( x1 )
    idx2 = numpy.argsort ( x2 )

    x1, w1 = x1 [ idx1 ] , w1 [ idx1 ]
    x2, w2 = x2 [ idx2 ] , w2 [ idx2 ]

    # 2. Compute normalized cumulative sum of weights
    sum_w1 = numpy.sum ( w1 )
    sum_w2 = numpy.sum ( w2 )

    cdf1 = numpy.cumsum ( w1 ) / sum_w1
    cdf2 = numpy.cumsum ( w2 ) / sum_w2

    # 3. Evaluate ECDFs on pooled points
    all_x = numpy.union1d ( x1 , x2 )

    cdf1_all = numpy.searchsorted ( x1 , all_x , side = "right" )
    cdf2_all = numpy.searchsorted ( x2 , all_x , side = "right" )

    cdf1_eval = numpy.pad ( cdf1 , ( 1 , 0 ) ) [ cdf1_all ]
    cdf2_eval = numpy.pad ( cdf2 , ( 1 , 0 ) ) [ cdf2_all ]

    # Clamp ECDFs to interval [eps, 1 - eps]
    F1 = numpy.clip ( cdf1_eval , eps , 1.0 - eps )
    F2 = numpy.clip ( cdf2_eval , eps , 1.0 - eps )

    # 4. Point-wise log-likelihood ratio supremum (max over all points)
    ll_terms = numpy.abs ( numpy.log ( F1 / F2 ) + numpy.log ( ( 1.0 - F1 ) / ( 1.0 - F2 ) ) )
    zk_stat  = float ( numpy.max ( ll_terms ) )

    return zk_stat

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
