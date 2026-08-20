#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/twosamples2.py
#  Two-Samples 1D-Weighted Tests
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-10-16
# =============================================================================
""" Two-Samples 1D Weighted Tests
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
from   ostap.stats.utils import weight_trivial , check_all , num_features

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
## Helper function to validate, sort, compute weighted ECDFs and evaluate them
#  on pooled unique data points.
def _prepare_ecdfs ( data1   ,
                     data2   ,
                     weight1 = None ,
                     weight2 = None ) :
    """ Helper function to validate, sort, compute weighted ECDFs and evaluate them
    on pooled unique data points.
    """
    
    x1 = numpy.asarray ( data1 , dtype = float )
    x2 = numpy.asarray ( data2 , dtype = float )

    check_all ( x1 , x2 , weight1 , weight2 )    
    if 1 != num_features ( x1 ) : raise ValueError ( "Invalid `num_features`" )
    
    x1 = x1.ravel () 
    x2 = x2.ravel () 

    w1_trivial = weight_trivial ( weight1 )
    w2_trivial = weight_trivial ( weight2 )
    
    w1 = numpy.asarray ( weight1 , dtype = float ) if not w1_trivial else numpy.ones_like ( x1 )
    w2 = numpy.asarray ( weight2 , dtype = float ) if not w2_trivial else numpy.ones_like ( x2 )
    
    N1 = len ( x1 )
    N2 = len ( x2 )
    if 2 > N1  or 2 > N2 : return None

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

    # 3. Evaluate ECDFs on pooled data points
    all_x = numpy.union1d ( x1 , x2 )

    cdf1_all = numpy.searchsorted ( x1 , all_x , side = "right" )
    cdf2_all = numpy.searchsorted ( x2 , all_x , side = "right" )

    cdf1_eval = numpy.pad ( cdf1 , ( 1 , 0 ) ) [ cdf1_all ]
    cdf2_eval = numpy.pad ( cdf2 , ( 1 , 0 ) ) [ cdf2_all ]

    # 4. Pooled weight step sizes dH(x)
    w1_pooled = numpy.bincount ( numpy.searchsorted ( all_x , x1 ) , weights = w1 , minlength = len ( all_x ) )
    w2_pooled = numpy.bincount ( numpy.searchsorted ( all_x , x2 ) , weights = w2 , minlength = len ( all_x ) )
    
    dH = ( w1_pooled + w2_pooled ) / ( sum_w1 + sum_w2 )

    return cdf1_eval , cdf2_eval , dH , all_x

# =============================================================================
## Calculates the two-sample Kolmogorov-Smirnov D-statistic for 1D datasets
#  with weights.
def kolmogorov_smirnov ( data1    ,
                         data2    ,
                         weight1 = None ,
                         weight2 = None ):
    """ Calculates the two-sample Kolmogorov-Smirnov D-statistic for 1D datasets
    with weights.
    """
    prepared = _prepare_ecdfs ( data1 , data2 , weight1 , weight2 )
    if prepared is None : return 0.0

    cdf1_eval , cdf2_eval , _ , _ = prepared

    diff   = cdf1_eval - cdf2_eval
    d_stat = float ( numpy.max ( numpy.abs ( diff ) ) )

    return d_stat


# =============================================================================
## Calculates the two-sample Kuiper statistic (V-statistic) for 1D datasets
#  with weights.
def kuiper ( data1          ,
             data2          ,
             weight1 = None ,
             weight2 = None ) :
    """ Calculates the two-sample Kuiper statistic for 1D datasets with weights.
    """
    prepared = _prepare_ecdfs ( data1 , data2 , weight1 , weight2 )
    if prepared is None : return 0.0

    cdf1_eval , cdf2_eval , _ , _ = prepared

    diff    = cdf1_eval - cdf2_eval
    d_plus  = float ( numpy.max ( diff ) )
    d_minus = float ( numpy.max ( -diff ) )
    
    v_stat  = max ( 0.0 , d_plus ) + max ( 0.0 , d_minus )

    return float ( v_stat )


# =============================================================================
## Calculates the two-sample Cramér-von Mises statistic (T / W^2) for 1D datasets
#  with weights.
def cramer_von_mises ( data1          ,
                       data2          ,
                       weight1 = None ,
                       weight2 = None ) :
    """ Calculates the two-sample Cramér-von Mises statistic for 1D datasets
    with weights.
    """
    prepared = _prepare_ecdfs ( data1 , data2 , weight1 , weight2 )
    if prepared is None : return 0.0

    cdf1_eval , cdf2_eval , dH , _ = prepared

    diff     = cdf1_eval - cdf2_eval
    cvm_stat = numpy.sum ( diff * diff * numpy.abs ( dH ) )

    return float ( cvm_stat )


# =============================================================================
## Calculates the two-sample Anderson-Darling statistic (A^2) for 1D datasets
#  with weights.
def anderson_darling ( data1           ,
                       data2           ,
                       weight1 = None  ,
                       weight2 = None  ,
                       eps     = 1e-10 ) :
    """ Calculates the two-sample Anderson-Darling statistic for 1D datasets
    with weights.
    """
    prepared = _prepare_ecdfs ( data1 , data2 , weight1 , weight2 )
    if prepared is None : return 0.0

    cdf1_eval , cdf2_eval , dH , _ = prepared

    # Monotonic accumulation for H(x) to prevent negative variance on sPlot weights
    H = numpy.cumsum ( dH )
    H = numpy.maximum.accumulate ( H )
    if H [ -1 ] != 0 : H = H / H [ -1 ]

    H_clipped       = numpy.clip ( H , eps , 1.0 - eps )
    variance_weight = H_clipped * ( 1.0 - H_clipped )

    diff    = cdf1_eval - cdf2_eval
    ad_stat = numpy.sum ( ( diff * diff / variance_weight ) * numpy.abs ( dH ) )

    return float ( ad_stat )


# =============================================================================
## Calculates the two-sample Berk-Jones statistic (R_BJ) for 1D datasets
#  with weights.
def berk_jones ( data1          ,
                 data2          ,
                 weight1 = None ,
                 weight2 = None ,
                 eps     = 1e-12 ) :
    """ Calculates the two-sample Berk-Jones statistic for 1D datasets
    with weights.
    """
    prepared = _prepare_ecdfs ( data1 , data2 , weight1 , weight2 )
    if prepared is None : return 0.0

    cdf1_eval , cdf2_eval , _ , _ = prepared

    F1 = numpy.clip ( cdf1_eval , eps , 1.0 - eps )
    F2 = numpy.clip ( cdf2_eval , eps , 1.0 - eps )

    kl_div  = F1 * numpy.log ( F1 / F2 ) + ( 1.0 - F1 ) * numpy.log ( ( 1.0 - F1 ) / ( 1.0 - F2 ) )
    bj_stat = float ( numpy.max ( kl_div ) )

    return bj_stat


# =============================================================================
## Calculates the two-sample Zhang's Z_A statistic for 1D datasets with weights.
def ZA ( data1           ,
         data2           ,
         weight1 = None  ,
         weight2 = None  ,
         eps     = 1e-12 ) :
    """ Calculates the two-sample Zhang's Z_A statistic for 1D datasets
    with weights.
    """
    prepared = _prepare_ecdfs ( data1 , data2 , weight1 , weight2 )
    if prepared is None : return 0.0

    cdf1_eval , cdf2_eval , dH , _ = prepared

    F1 = numpy.clip ( cdf1_eval , eps , 1.0 - eps )
    F2 = numpy.clip ( cdf2_eval , eps , 1.0 - eps )

    ll_terms = numpy.log ( F1 ) + numpy.log ( 1.0 - F1 ) + numpy.log ( F2 ) + numpy.log ( 1.0 - F2 )
    za_stat  = - numpy.sum ( ll_terms * numpy.abs ( dH ) )

    return float ( za_stat )


# =============================================================================
## Calculates the two-sample Zhang's Z_C statistic for 1D datasets with weights.
def ZC ( data1           ,
         data2           ,
         weight1 = None  ,
         weight2 = None  ,
         eps     = 1e-12 ) :
    """ Calculates the two-sample Zhang's Z_C statistic for 1D datasets
    with weights.
    """
    prepared = _prepare_ecdfs ( data1 , data2 , weight1 , weight2 )
    if prepared is None : return 0.0

    cdf1_eval , cdf2_eval , dH , _ = prepared

    F1 = numpy.clip ( cdf1_eval , eps , 1.0 - eps )
    F2 = numpy.clip ( cdf2_eval , eps , 1.0 - eps )

    ll_terms = ( F1 - F2 ) ** 2 / ( F1 * ( 1.0 - F1 ) + F2 * ( 1.0 - F2 ) )
    zc_stat  = numpy.sum ( ll_terms * numpy.abs ( dH ) )

    return float ( zc_stat )


# =============================================================================
## Calculates the two-sample Zhang's Z_K statistic for 1D datasets with weights.
def ZK ( data1           , 
         data2           ,
         weight1 = None  ,
         weight2 = None  ,
         eps     = 1e-12 ) :
    """ Calculates the two-sample Zhang's Z_K statistic for 1D datasets
    with weights.
    """
    prepared = _prepare_ecdfs ( data1 , data2 , weight1 , weight2 )
    if prepared is None : return 0.0

    cdf1_eval , cdf2_eval , _ , _ = prepared

    F1 = numpy.clip ( cdf1_eval , eps , 1.0 - eps )
    F2 = numpy.clip ( cdf2_eval , eps , 1.0 - eps )

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

