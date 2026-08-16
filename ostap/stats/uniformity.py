#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/uniformity.py
#  Simple test to check uniform distribution
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-09-16
# =============================================================================
""" Simple test to check uniform distributiona
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2026-07-04"
__all__     = (
    'greenwood_t_value'          , ## Greenwood          uniformity test 
    'cramer_von_mises_t_value'   , ## Cramer-von Mises   unifromity test 
    'logarithmic_tail_t_value'   , ## Logarithmis tail   uniformity test 
    'kolmogorov_smirnov_t_value' , ## Kolmogorov-Smirnov uniformity test 
)
# =============================================================================
import numpy, math 
# =============================================================================
# Logging setup 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.uniformity' )
else                       : logger = getLogger ( __name__ )
# =============================================================================
logger.debug ( "Test to check uniform distributions" )
# =============================================================================
## Calculates the standardized t-value (z-score) using Greenwood's Spacing Statistic 
#  to test if the u-values are Uniformly distributed on [0, 1].
# 
#  Ideal for Distance-to-Nearest-Neighbor (DNN) Goodness-of-Fit tests, 
# as it has maximum power against local clustering and density flaws.
#
# @param uvalues Vector of u-values in range [0, 1].
# @return Standardized t-value (~ N(0,1) under H0). Values near 0 indicate good fit,
#         large positive values (|t| >> 2) indicate rejection of uniformity.
def greenwood_t_value ( uvalues ) :
    """ Calculates the standardized t-value (z-score) using Greenwood's Spacing
    Statistic to test if the u-values are Uniformly distributed on [0, 1].

    Ideal for Distance-to-Nearest-Neighbor (DNN) Goodness-of-Fit tests.
    """
    u = numpy.asarray ( uvalues , dtype = numpy.float64 )
    N = len ( u )
    if N < 2 : return 0.0

    # 1. Sort u-values in ascending order
    u_sorted = numpy.sort ( u )

    # 2. Compute the sum of squared spacings: S = sum(D_i^2)
    # Account for boundaries: u_0 = 0.0 and u_{N+1} = 1.0
    
    u_ext  = numpy.empty(N + 2, dtype=float)
    u_ext [0   ] = 0.0
    u_ext [1:-1] = u_sorted
    u_ext [  -1] = 1.0

    spacings        = numpy.diff ( u_ext )
    sum_sq_spacings = numpy.sum  ( spacings * spacings )

    # 3. Theoretical Mean and Variance under H0 ~ U(0,1)
    n_dbl  = float ( N )
    mean_S = 2.0 / (n_dbl + 2.0)
    var_S  = ( 4.0 * n_dbl ) / ( ( ( n_dbl + 2.0 ) ** 2 ) * ( n_dbl + 3.0 ) )

    # 4. Return standardized t-value (z-score)
    return float ( ( sum_sq_spacings - mean_S ) / math.sqrt ( var_S ) )

# ============================================================================
## Calculate standardized t-value using the Cramér–von Mises statistic (W^2).
#  Excellent all-round Goodness-of-Fit test for general shape deviations.
# 
#  @param uvalues Vector of u-values in range [0, 1].
#  @return Standardized t-value (~ N(0,1) under H0).
def cramer_von_mises_t_value(uvalues):
    """ Calculates standardized t-value using the Cramér–von Mises statistic
    (W^2).

    Excellent all-round Goodness-of-Fit test for general shape deviations.
    """
    u = numpy.asarray ( uvalues , dtype = numpy.float64 )
    N = len ( u )
    if N < 2 : return 0.0

    # 1. Sort u-values in ascending order
    u_sorted = numpy.sort ( u )

    # 2. Compute Cramér–von Mises W^2 statistic
    n_dbl = float ( N )
    i     = numpy.arange ( N , dtype = numpy.float64 )
    e_i   = ( 2.0 * i + 1.0 ) / ( 2.0 * n_dbl )

    diff  = u_sorted - e_i
    sum_sq_diff = numpy.sum ( diff * diff )

    W2 = sum_sq_diff + ( 1.0 / ( 12.0 * n_dbl ) )

    # 3. Theoretical Mean and Variance under H0 ~ U(0,1)
    mean_W2 = 1.0 / 6.0
    var_W2  = ( 4.0 * n_dbl - 3.0 ) / ( 180.0 * n_dbl * n_dbl )

    # 4. Return standardized t-value
    return float ( ( W2 - mean_W2 ) / math.sqrt ( var_W2 ) )

# ============================================================================
## Calculate standardized t-value using Logarithmic Tail Mean (-E[ln(u) + ln(1-u)]).
#  Extremely sensitive to tail behavior (u -> 0 or u -> 1).
#  Fast O(N) performance — does NOT require sorting.
#
#  @param uvalues Vector of u-values in range [0, 1].
#  @param delta Safety margin to avoid log(0).
#  @return Standardized t-value (~ N(0,1) under H0).
def logarithmic_tail_t_value (uvalues, delta = 1e-12 ) :
    """ Calculate standardized t-value using Logarithmic Tail Mean
    (-E[ln(u) + ln(1-u)]).

    Extremely sensitive to tail behavior (u -> 0 or u -> 1).
    Fast O(N) performance — does NOT require sorting.
    """
    u = numpy.asarray ( uvalues , dtype = numpy.float64 )
    N = len  ( u )
    if N < 2 : return 0.0

    n_dbl = float( N )

    # 1. Accumulate -ln(u) - ln(1 - u) with clamp to safe interval [delta, 1 - delta]
    u_clamped = numpy.clip ( u, delta, 1.0 - delta )
    sum_log   = numpy.sum  ( numpy.log( u_clamped ) + numpy.log ( 1.0 - u_clamped ) ) 

    mean_log_stat = -sum_log / n_dbl

    # 2. Theoretical Mean = 2.0, Variance = (2*pi^2/3 - 4) / N approx 2.5797 / N
    mean_H0 = 2.0
    var_H0  = ( 2.0 * ( math.pi**2 ) / 3.0 - 4.0 ) / n_dbl

    # 3. Return standardized t-value
    return float ( ( mean_log_stat - mean_H0 ) / math.sqrt ( var_H0 ) )


# ============================================================================
## Calculate standardized t-value (Stephens-adjusted) for Kolmogorov-Smirnov test.
#   Measures maximum absolute deviation between empirical and theoretical CDFs.
#   @param uvalues Vector of u-values in range [0, 1].
#   @return Standardized KS statistic (~ 0 for uniform, > 1.358 for p < 0.05).
def kolmogorov_smirnov_t_value ( uvalues ) :
    """ Calculates standardized t-value (Stephens-adjusted) for Kolmogorov-
    Smirnov test.

    Measures maximum absolute deviation between empirical and theoretical CDFs.
    """
    u = numpy.asarray ( uvalues , dtype = numpy.float64 )
    N = len ( u )
    if N < 2 : return 0.0

    # 1. Sort u-values
    u_sorted = numpy.sort ( u )

    n_dbl = float ( N )
    i_curr = numpy.arange ( 1, N + 1, dtype = numpy.float64 )

    # 2. Compute max positive and negative deviations: D+ and D-
    d_plus  = ( i_curr / n_dbl ) - u_sorted
    d_minus = u_sorted -  ( ( i_curr - 1.0 ) / n_dbl )

    max_d   = max ( numpy.max ( d_plus ), numpy.max ( d_minus ) )

    # 3. Apply Stephens (1970) modification factor for standardized T-value
    sqrt_n  = math.sqrt ( n_dbl )
    t_value = max_d * ( sqrt_n + 0.12 + 0.11 / sqrt_n )

    return float ( t_value )
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# ============================================================================
#                                                                      The END
# ============================================================================
