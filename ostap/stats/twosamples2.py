#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/twosamples2.py
#  Two-Samples Test for 1D *weighted* data 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-10-16
# =============================================================================
""" Two-Samples Test for 1D *weighted* data 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-16"
__all__     = (
    'kolmogorov_smirnov' , ## Kolmogorov-Smirnov 2-samples test 
    'kuiper'             , ## Kuiper             2-samples test 
    'cramer_von_mises'   , ## Cramer-von Mises   2-samples test 
    'anderson_darling'   , ## Anderson-Darling   2-samples test 
    'berk_jones'         , ## Berk-Jones         2-samples test 
    'ZA'                 , ## Zhang's ZA         2-samples test 
    'ZC'                 , ## Zhang's ZC         2-samples test 
    'ZK'                 , ## Zhang's ZK         2-samples test
    'wasserstein1'       , ## The first Wasserstein 2-sample test
    'wasserstein2'       , ## The second Wasserstein 2-sample test
    'rosenblatt'         , ## Rosenblatt 2-sample test 
    'watson_u2'          , ## Watson U2 2-sample test
    ## Rank-based tests:  *NO WEIGHTS* 
    'mann_whitney'       , ## Mann-Whitney       2-samples test (NO WEIGHTS!)
    'savage'             , ## symmetrized Savage 2-samples test (NO WEIGHTS!)
    ## Rank-based tests:  *WITH WEIGHTS* 
    'ansari_bradley'     , ## Ansari-Bradley     2-samples test (WITH WEIGHTS!)
    'siegel_tukey'       , ## Siegel-Tukey       2-samples test (WITH WEIGHTS!)
    'mood'               , ##  Mood              2-samples test (WITH WEIGHTS!)
    'van_der_waerden'    , ## Van der Waerden    2-samples test (WITH WEIGHTS!)
    'brunner_munzel'     , ## Brunner-Munzel     2-samples test (WITH WEIGHTS!)
    'epps_singleton'     , ## Epps-Singleton     2-samples test (WITH WEIGHTS!)
    ## 
) 
# =============================================================================
from   ostap.stats.utils    import weight_trivial, data2ecdf
from   ostap.math.math_base import Ostap
from   ostap.utils.ranges   import lrange 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
logger = getLogger ( 'ostap.stats.twosamples2' )
# =============================================================================
logger.debug ( 'Two-Samples 1D (weighted) Test' )
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
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.kolmogorov_smirnov ( ecdf1 , ecdf2 )

# =============================================================================
## Calculates the two-sample Kuiper statistic (V-statistic) for 1D datasets
#  with weights.
def kuiper ( data1          ,
             data2          ,
             weight1 = None ,
             weight2 = None ) :
    """ Calculates the two-sample Kuiper statistic for 1D datasets with weights.
    """
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.kuiper ( ecdf1 , ecdf2 )

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
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.cramer_von_mises  ( ecdf1 , ecdf2 )

# =============================================================================
## Calculates the two-sample Anderson-Darling statistic (A^2) for 1D datasets
#  with weights.
def anderson_darling ( data1           ,
                       data2           ,
                       weight1 = None  ,
                       weight2 = None  ,
                       eps     = 1e-12 ) :
    """ Calculates the two-sample Anderson-Darling statistic for 1D datasets
    with weights.
    """
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.anderson_darling ( ecdf1 , ecdf2 , eps )

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
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.berk_jones ( ecdf1 , ecdf2 , eps )

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
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.ZA ( ecdf1 , ecdf2 , eps )

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
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.ZC ( ecdf1 , ecdf2 , eps )

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
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.ZK ( ecdf1 , ecdf2 , eps )

# =============================================================================
## Calculates First Wasserstein distance \f$W_1\f$ / Earth Mover's Distance.
#
#  Calculates the 1st Wasserstein distance between two empirical distributions:
#  \f[ W_1(F_1, F_2) = \int_{-\infty}^{\infty} |F_1(x) - F_2(x)| \, dx \f]
#
#  @see https://en.wikipedia.org/wiki/Wasserstein_metric
def wasserstein1 ( data1           , 
                   data2           ,
                   weight1 = None  ,
                   weight2 = None  ) : 
    """ Calculates the First Wasserstein distance / Earth Mover's Distance).
    @see https://en.wikipedia.org/wiki/Wasserstein_metric
    """
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.wasserstein1 ( ecdf1 , ecdf2 )

# =============================================================================
## Calculates First Wasserstein distance (\f$W_1\f$ / Earth Mover's Distance).
#
#  Calculates the square root of the integrated squared differences between empirical CDFs:
#  \f[ W_2(F_1, F_2) = \left( \int_{-\infty}^{\infty} (F_1(x) - F_2(x))^2 \, dx \right)^{1/2} \f]
#
#
#  @see https://en.wikipedia.org/wiki/Wasserstein_metric
def wasserstein2 ( data1           , 
                   data2           ,
                   weight1 = None  ,
                   weight2 = None  ) : 
    """ Calculates the second Wasserstein distance .
    @see https://en.wikipedia.org/wiki/Wasserstein_metric
    """
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.wasserstein2 ( ecdf1 , ecdf2 )

# =============================================================================
## @brief Rosenblatt-type two-sample test statistic (\f$R\f$).
#  Calculates cumulative squared discrepancies between two empirical CDFs 
#  evaluated on the pooled sample points:
#   \f[ R = \sum_i (F_1(x_i) - F_2(x_i))^2 \f]
# 
#   @see Rosenblatt, M. (1952)
def rosenblatt ( data1           , 
                 data2           ,
                 weight1 = None  ,
                 weight2 = None  ) : 
    """ Calculates Rosenblatt distance .
    """
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.rosenblatt ( ecdf1 , ecdf2 )

# =============================================================================
## @brief Watson's two-sample test statistic (\f$U^2\f$).
#
#  A Cramér-von Mises type variant invariant under cyclic transformations:
#  \f[ U^2 = \int_{-\infty}^{\infty} \left( (F_1(x) - F_2(x)) - \bar{D} \right)^2 dH(x) \f]
#  @see https://en.wikipedia.org/wiki/Watson_test
#   @note Watson, G.S. (1961, 1962).
def watson_u2 ( data1           , 
                data2           ,
                weight1 = None  ,
                weight2 = None  ) : 
    """ Calculates Watson U2 distance .
    """
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.watson_u2 ( ecdf1 , ecdf2 )

# =============================================================================
## Matt-Whitney statistics 
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  mw    = mann_whitney ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Matt-Whitney statistics 
def mann_whitney ( data1   ,
                   data2   ,
                   weight1 = None , 
                   weight2 = None ) : 
    """ Get Mann-Whitney statistics 
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ##
    if not weight_trivial ( weight1 ) : raise ValueError ( "Mann-Whitney does not support weights!" )
    if not weight_trivial ( weight2 ) : raise ValueError ( "Mann-Whitney does not support weights!" )    
    ## 
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.mann_whitney ( ecdf1 , ecdf2 )

# =============================================================================
## Symmetrized Savage statistics 
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  ss    = savage ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Savage statistics 
def savage ( data1   ,
             data2   ,
             weight1 = None , 
             weight2 = None ) :              
    """ Get Savage statistics 
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ##
    if not weight_trivial ( weight1 ) : raise ValueError ( "Savage does not support weights!" )
    if not weight_trivial ( weight2 ) : raise ValueError ( "Savage does not support weights!" )    
    ## 
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.savage ( ecdf1 , ecdf2 )


# =============================================================================
## Ansari-Bradley two-sample scale test statistic.
# 
# Rank-based nonparametric test for equality of scale (dispersion) parameters
# assuming equal medians. Assigns smaller ranks to extreme observations and
# larger ranks to central observations. Returns the absolute standardized
# deviation from the null expectation to provide a two-sided test statistic.
# 
# @see https://en.wikipedia.org/wiki/Ansari%E2%80%93Bradley_test
# @note Ansari, A.R., Bradley, R.A. (1960). "Rank-sum tests for dispersions".
#       The Annals of Mathematical Statistics, 31(4): 1174–1189.
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Andari-Bradley statistics 
def ansari_bradley ( data1   ,
                     data2   ,
                     weight1 = None , 
                     weight2 = None ) :              
    """ Get Ansari-Bradley statistics 
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.ansari_bradley ( ecdf1 , ecdf2 )

# =============================================================================
##  Siegel-Tukey two-sample scale test statistic.
#
# Rank-based nonparametric test for differences in scale/variability.
# Assigns ranks in a serpentine pattern from the extremes toward the center
#  (1, 4, 5, 8, 9... from left; 2, 3, 6, 7, 10... from right). Returns the
# absolute standardized deviation from the null expectation.
#
# @see https://en.wikipedia.org/wiki/Siegel%E2%80%93Tukey_test
# @note Siegel, S., Tukey, J.W. (1960). "A nonparametric sum of ranks procedure
#         for relative spread in unpaired samples". JASA, 55(291): 429–445.
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Siegel-Tukey statistics 
def siegel_tukey ( data1   ,
                   data2   ,
                   weight1 = None , 
                   weight2 = None ) :              
    """ Get Siegel-Tukey statistics 
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.siegel_tukey ( ecdf1 , ecdf2 )

# =============================================================================
## Mood's two-sample scale test statistic.
#
# Rank-based nonparametric test for equality of scale parameters based on
# squared deviations of ordinary ranks from the median rank (N+1)/2.
# Returns the absolute standardized deviation from the null expectation.
#  
# @note Mood, A.M. (1954). "On the asymptotic efficiency of certain nonparametric
#   two-sample tests". The Annals of Mathematical Statistics, 25(3): 514–522.
#
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Mood statistics 
def mood  ( data1   ,
            data2   ,
            weight1 = None , 
            weight2 = None ) :              
    """ Get Mood statistics 
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.mood  ( ecdf1 , ecdf2 )

# =============================================================================
## Van der Waerden (normal scores) two-sample location test statistic.
#  
# Nonparametric rank-based test that converts midranks \f$ r_i \f$ into 
# standard normal quantiles (probit transformation) using Ostap::Math::probit:
# \f[ a(r_i) = \Phi^{-1}\left( \frac{r_i}{N + 1} \right) \f]
#  where \f$ \Phi^{-1} \f$ is the inverse standard normal CDF.
#
# Achieves the same asymptotic relative efficiency as the parametric Student's
# t-test under normality while remaining strictly distribution-free. Returns 
# the standardized absolute deviation from the null expectation (asymptotic 
# \f$ |Z| \f$-score) with exact finite-sample variance accounting for ties.
#
# Standardized two-sided Van der Waerden statistic \f$ |Z| \f$
# @see https://en.wikipedia.org/wiki/Van_der_Waerden_test
# @note Van der Waerden, B.L. (1952). "Order tests for the two-sample problem 
#       and their power". Indagationes Mathematicae, 14: 453–458.
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Van der Waerden statistics 
def van_der_waerden ( data1   ,
                      data2   ,
                      weight1 = None , 
                      weight2 = None ) :              
    """ Get Van der Waerden statistics 
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.van_der_waerden ( ecdf1 , ecdf2 )


# =============================================================================
## Brunner-Munzel two-sample rank test statistic.
#
# Nonparametric test for stochastic equality (generalized Mann-Whitney test)
# that remains robust under unequal variances and heterogeneous distributions 
# (Behrens-Fisher problem for ranks). 
#  Supports both unweighted (\f$ ECDF \f$) and weighted (\f$ WECDF \f$) data containers.
#
# @see https://en.wikipedia.org/wiki/Brunner%E2%80%93Munzel_test
# @note Brunner, E., Munzel, U. (2000). "The nonparametric Behrens-Fisher problem: 
#      asymptotic theory and a small-sample approximation". 
#      Biometrical Journal, 42(1): 17–25.
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Standardized two-sided Brunner-Munzel test statistic (|Z|-score) 
def brunner_munzel ( data1   ,
                     data2   ,
                     weight1 = None , 
                     weight2 = None ) :              
    """ Get Brunner-Munzel statistics 
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.brunner_munzel ( ecdf1 , ecdf2 )


# =============================================================================
## Get an array of Epps-Singleton frequencins
#  N+1 points from low to high uniformly in log-scale  
def epps_singleton_frequencies ( * , 
                                 low  = 0.2 , 
                                 high = 10  , 
                                 N    = 12  ) :
    """ Get an array  of Epps-Singleton frequencies
    - N+1 points from low to high uniformly in log-scale 
    """
    if not   0 < low < high                      : raise ValueError ( "Epps-Sigleton frequencies: invalid `low/high' %s,%s" % ( low , high ) )
    if not ( isinstance ( N , int ) and 1 <= N ) : raise ValueError ( "Epps-Sigleton frequencies: invalid `N' %s" % N )
    return tuple (  i for i in lrange ( low , high , N ) )

# =============================================================================
## Epps-Singleton two-sample test statistic.
# 
# Omnibus nonparametric test based on the empirical characteristic function (ECF) 
# evaluated at specific frequency points. Highly sensitive to differences in 
# distribution shape, skewness, and tails.
# Supports both unweighted (\f$ ECDF \f$) and weighted (\f$ WECDF \f$) data containers.
#  
# @attention Frequencies are defined in terms of standardized data! 
# @see https://en.wikipedia.org/wiki/Epps%E2%80%93Singleton_test
# @note Epps, T.W., Singleton, K.J. (1986). "An omnibus test for the two-sample 
#       problem using the empirical characteristic function". 
#       Journal of Statistical Computation and Simulation, 26(3-4): 177–203.
#
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Epps-Singleton test statistic
def epps_singleton ( data1       ,
                     data2       ,
                     weight1     = None , 
                     weight2     = None , 
                     frequencies = epps_singleton_frequencies() ) :              
    """ Get Epps-Singleton statistics 
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    return Ostap.Math.epps_singleton ( ecdf1 , ecdf2 , frequencies )


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
