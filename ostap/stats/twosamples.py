#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/twosamples.py
#  Two-Sample Tests for 1D datasets 
#  @see https://www.jstor.org/stable/25471118
#  For dataset with weights, see 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-10-16
# =============================================================================
""" Two-Sample Tests for 1D datasets 
 - see https://www.jstor.org/stable/25471118
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-16"
__all__     = (
    ##  Classical Two-Samples Tests: 
    'kolmogorov_smirnov' , ## Kolmogorov-Smirnov     2-samples test 
    'kuiper'             , ## Kuiper                 2-samples test 
    'cramer_von_mises'   , ## Cramer-von Mises       2-samples test 
    'anderson_darling'   , ## Anderson-Darling       2-samples test 
    'berk_jones'         , ## Berk-Jones             2-samples test 
    'ZA'                 , ## Zhang's ZA             2-samples test 
    'ZC'                 , ## Zhang's ZC             2-samples test 
    'ZK'                 , ## Zhang's ZK             2-samples test
    'wasserstein1'       , ## The first  Wasserstein 2-samples test
    'wasserstein2'       , ## The second Wasserstein 2-samples test
    'rosenblatt'         , ## Rosenblatt             2-samples test 
    'watson_u2'          , ## Watson U2              2-samples test
    ## Rank-based Two-Samples Tests:  *WITH WEIGHTS* 
    'ansari_bradley'     , ## Ansari-Bradley         2-samples test (WITH WEIGHTS!)
    'siegel_tukey'       , ## Siegel-Tukey           2-samples test (WITH WEIGHTS!)
    'mood'               , ## Mood                   2-samples test (WITH WEIGHTS!)
    'van_der_waerden'    , ## Van der Waerden        2-samples test (WITH WEIGHTS!)
    'brunner_munzel'     , ## Brunner-Munzel         2-samples test (WITH WEIGHTS!)
    'epps_singleton'     , ## Epps-Singleton         2-samples test (WITH WEIGHTS!) 
    ## Rank-based Two-Samples Tests:  *NO WEIGHTS* 
    'mann_whitney'       , ## Mann-Whitney           2-samples test (NO WEIGHTS!)
    'savage'             , ## (symmetrized) Savage   2-samples test (NO WEIGHTS!)
    ## 
    "TwoSamplesTest"     , ## Two-Samples Test
    "TwoSamplesToys"     , ## Toys/permutations for Two-Samples Test 
    ## 
)
# =============================================================================
from   ostap.core.core         import VE, Ostap
from   ostap.math.math_base    import axis_range, data2vct, pos_infinity
from   ostap.utils.cidict      import cidict_fun
from   ostap.logger.pretty     import pretty_float, format_pow10 
from   ostap.stats.gof_utils   import Labels, Keys, clip_pvalue 
from   ostap.stats.utils       import full_header, t_header  
from   ostap.math.math_ve      import significance
from   ostap.math.ve           import fmt_pretty_ve
from   ostap.stats.counters    import SE, EffCounter, ECDF
from   ostap.stats.utils       import ( data2ecdf     , weight_trivial , 
                                        clip_pvalue   , 
                                        format_pvalue , format_nsigma  )  
from   ostap.stats.pvalue      import make_permutations
from   ostap.utils.ranges      import lrange 
from   collections             import defaultdict, namedtuple
import ostap.logger.table      as     T
import ostap.math.models 
import ostap.fitting.ds2numpy 
import ostap.fitting.roofit
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
logger = getLogger ( 'ostap.stats.twosamples' )
# =============================================================================
logger.debug ( 'Two-Samples Test for 1D datasets' )
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
## @class TwoSamplesTest 
#  Two-Samples Test
class TwoSamplesTest ( object ) :
    """ Two-sample test 
    """
    def __init__ ( self        ,
                   data1       ,
                   data2       , 
                   weight1     = None ,
                   weight2     = None , *  , 
                   frequencies = epps_singleton_frequencies () ) :
        
        ## transform input data 
        ecdf1 = data2ecdf ( data1 , weight1 )
        ecdf2 = data2ecdf ( data2 , weight2 )
        
        self.__frequencies = tuple  ( frequencies )
            
        self.__ecdf1   = ecdf1
        self.__ecdf2   = ecdf2
        self.__weight1 = weight1 
        self.__weight2 = weight2 
        self.__ecdf    = ecdf1.merge ( ecdf2 ) 
        
        self.__estimators = {            
            'KS'  : kolmogorov_smirnov ( self.ecdf1 , self.ecdf2 ) , 
            'K'   : kuiper             ( self.ecdf1 , self.ecdf2 ) , 
            'AD'  : anderson_darling   ( self.ecdf1 , self.ecdf2 ) , 
            'CM'  : cramer_von_mises   ( self.ecdf1 , self.ecdf2 ) , 
            'BJ'  : berk_jones         ( self.ecdf1 , self.ecdf2 ) , 
            'ZK'  : ZK                 ( self.ecdf1 , self.ecdf2 ) , 
            'ZA'  : ZA                 ( self.ecdf1 , self.ecdf2 ) , 
            'ZC'  : ZC                 ( self.ecdf1 , self.ecdf2 ) , 
            'W1'  : wasserstein1       ( self.ecdf1 , self.ecdf2 ) , 
            'W2'  : wasserstein2       ( self.ecdf1 , self.ecdf2 ) , 
            'R'   : rosenblatt         ( self.ecdf1 , self.ecdf2 ) , 
            'U2'  : watson_u2          ( self.ecdf1 , self.ecdf2 ) ,
            'AB'  : ansari_bradley     ( self.ecdf1 , self.ecdf2 ) , 
            'ST'  : siegel_tukey       ( self.ecdf1 , self.ecdf2 ) , 
            'MD'  : mood               ( self.ecdf1 , self.ecdf2 ) , 
            'VDW' : van_der_waerden    ( self.ecdf1 , self.ecdf2 ) , 
            'BM'  : brunner_munzel     ( self.ecdf1 , self.ecdf2 ) ,
            'ES'  : epps_singleton     ( self.ecdf1 , self.ecdf2 , frequencies = self.frequencies ) ,   
        }
        if  weight_trivial ( weight1 ) and weight_trivial ( weight2 ) :
            self.__estimators [ 'MW' ] = mann_whitney ( self.ecdf1 , self.ecdf2 )
            self.__estimators [ 'S'  ] = savage       ( self.ecdf1 , self.ecdf2 )
        
    @property
    def ecdf1 ( self ) :
        """`ecdf1` : empirical CDF for the 1st dataset
        """
        return self.__ecdf1 
    @property
    def ecdf2 ( self ) :
        """`ecdf2` : empirical CDF for the 2nd dataset
        """
        return self.__ecdf2
    
    @property
    def ecdf ( self ) :
        """`ecdf` : total/combined CDF for both datasets
        """
        return self.__ecdf
    
    @property
    def weight1  ( self  ) :
        """`weight1`: weight for the 1st dataset"""
        return self.__weight1 
    
    @property
    def weight2  ( self  ) :
        """`weight2`: weight for the 2nd dataset"""
        return self.__weight2 
        
    @property
    def frequencies ( self ) :
        """`frequencies`: frequences for Epps-Singleton method"""
        return self.__frequencies 
    
    # =========================================================================
    ## Print the summary as Table
    def table ( self             , * , 
                title     = ''   ,
                prefix    = ''   ,
                precision = 4    , 
                width     = 6    ,
                style     = None ) :
        """ Print the summary as Table
        """
        ##
        header = ( 'Statistics' , )  + t_header 
        rows   = [] 
        
        for label , value  in self.estimators.items () :
            
            the_label = Labels.get ( label , label )
            
            result , expo = pretty_float ( value , width = width , precision = precision )            
            row = the_label , result , ( '[%s]' % format_pow10 ( expo ) )  if expo else ''
            
            rows.append ( row )
            
        rows  = [ header ] + sorted ( rows ) 
        rows  = T.remove_empty_columns ( rows )
        title = title if title else 'Two-Samples Test'    
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lcl' , style = style  )

    ## print estimator as table 
    __repr__ = table
    __str__  = table
    
    # =========================================================================
    @property
    def estimators ( self ) :
        """`estimators` : get a dictionary of all estimators"""
        return self.__estimators

    # =========================================================================
    ## draw all involved empirical CDFS
    def draw  ( self , option = '' , *options , **kwargs ) :
        """ Draw fit CDF & empirical CDF
        """

        ecdf1 = self.ecdf1
        ecdf2 = self.ecdf2
        ecdf  = self.ecdf

        xmin  , xmax  = ecdf .xmin (), ecdf .xmax()
        xmin1 , xmax1 = ecdf1.xmin (), ecdf1.xmax()
        xmin2 , xmax2 = ecdf2.xmin (), ecdf2.xmax()

        xmin = min ( xmin , xmin1 , xmin2 )
        xmax = max ( xmax , xmax1 , xmax2 )

        xmin , xmax = axis_range ( xmin , xmax , delta = 0.10 )

        xmin = kwargs.pop ( 'xmin' , xmin )
        xmax = kwargs.pop ( 'xmax' , xmax )
        
        option  = option.strip() 
        optsame = 'same' if not option else '%s %s' % ( 'same' , option ) 
        r  = ecdf .draw ( option  , *options , color = 8 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )
        r1 = ecdf1.draw ( optsame , *options , color = 2 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )
        r2 = ecdf2.draw ( optsame , *options , color = 4 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )
        
        return r, r1, r2

# =============================================================================
## @class TwoSamplesToys
#  Toys for Two-Samples Test
class TwoSamplesToys(TwoSamplesTest):
    """ Toys for Two-sample test 
    """
    ## result of 2-samples test
    Result = namedtuple ( 'Result' , 'statistics counter pvalue nsigma' )
    # =========================================================================
    ## Initialize TSToys object 
    def __init__ ( self            ,
                   data1           ,
                   data2           ,
                   weight1 = None  ,
                   weight2 = None  , * , 
                   nToys   = 100   , 
                   silent  = False , **kwargs ) :
        
        ## initialize the base
        super().__init__ (  data1   = data1   , 
                            data2   = data2   , 
                            weight1 = weight1 , 
                            weight2 = weight2 , **kwargs ) 

        self.__counters = defaultdict(SE) 
        self.__ecdfs    = {} 
        self.__silent   = True if silent else False
        
        self.__counter  = EffCounter() 
        ## self.__counter += True 
        
        self.__nToys  = 0 
        if 0 < nToys : self.run ( nToys , silent = silent ) 
    
    # =========================================================================
    ## number of toys/permutations 
    @property
    def nToys ( self ) :
        """`nToys` : number of toys/permutations"""
        return self.__nToys
    
    # =========================================================================
    ## efficienty counter for p-value
    @property
    def counter ( self ) : 
        """`counter`: efficinect counter for p-value
        """
        return self.__counter
     
    # =========================================================================
    ## ECDFs
    @property 
    def ecdfs ( self ) :
        """`ecdfs` : toys/permutations results as empirical cumulative distribution functions"""
        return self.__ecdfs
    
    # =========================================================================
    ## Counters  
    @property 
    def counters ( self ) :
        """`counters` : toy results as counters"""
        return self.__counters
    
    # ===============================================================================
    ## run toys 
    def run ( self , nToys = 1000 , silent = False ) :
        """ Run toys 
        """ 
        assert isinstance ( nToys , int ) and 0 < nToys , "Invalid `nToys` argument!"

        n1   = len ( self.ecdf1 )
        n2   = len ( self.ecdf1 )
        
        data1 , weight1 = self.ecdf1.raw_data ()
        data2 , weight2 = self.ecdf2.raw_data ()
         
        results    = defaultdict(list)
        counters   = self.counters
        estimators = self.estimators 
        
        no_weight = weight_trivial ( weight1 ) and weight_trivial ( weight2 )
        
        ref_ks  = estimators.get ( 'KS'  , pos_infinity )
        ref_k   = estimators.get ( 'K'   , pos_infinity )
        ref_ad  = estimators.get ( 'AD'  , pos_infinity )
        ref_cm  = estimators.get ( 'CM'  , pos_infinity )
        ref_bj  = estimators.get ( 'BJ'  , pos_infinity )
        ref_zk  = estimators.get ( 'ZK'  , pos_infinity )
        ref_za  = estimators.get ( 'ZA'  , pos_infinity )
        ref_zc  = estimators.get ( 'ZC'  , pos_infinity )
        ref_w1  = estimators.get ( 'W1'  , pos_infinity )
        ref_w2  = estimators.get ( 'W2'  , pos_infinity )
        ref_r   = estimators.get ( 'R'   , pos_infinity )
        ref_u2  = estimators.get ( 'U2'  , pos_infinity )
        ref_ab  = estimators.get ( 'AB'  , pos_infinity )
        ref_st  = estimators.get ( 'ST'  , pos_infinity )
        ref_md  = estimators.get ( 'MD'  , pos_infinity ) 
        ref_vdw = estimators.get ( 'VDW' , pos_infinity )
        ref_bm  = estimators.get ( 'BM'  , pos_infinity )
        ref_es  = estimators.get ( 'ES'  , pos_infinity )
        ref_mw  = estimators.get ( 'MW'  , pos_infinity )
        ref_s   = estimators.get ( 'S'   , pos_infinity ) 
                             
                                
        for ds1, ds2, w1, w2 in make_permutations ( nToys , 
                                                      data1   = data1   , 
                                                      data2   = data1   , 
                                                      weight1 = weight1 , 
                                                      weight2 = weight2 ) :
        
            ## transform input data 
            ecdf1 = data2ecdf ( ds1 , w1  )
            ecdf2 = data2ecdf ( ds2 , w2  )
                
            ks  = kolmogorov_smirnov ( ecdf1 , ecdf2 )
            k   = kuiper             ( ecdf1 , ecdf2 )
            ad  = anderson_darling   ( ecdf1 , ecdf2 )            
            cm  = cramer_von_mises   ( ecdf1 , ecdf2 )             
            bj  = berk_jones         ( ecdf1 , ecdf2 )             
            zk  = ZK                 ( ecdf1 , ecdf2 )
            za  = ZA                 ( ecdf1 , ecdf2 )
            zc  = ZC                 ( ecdf1 , ecdf2 )
            w1  = wasserstein1       ( ecdf1 , ecdf2 )             
            w2  = wasserstein2       ( ecdf1 , ecdf2 )             
            r   = rosenblatt         ( ecdf1 , ecdf2 )             
            u2  = watson_u2          ( ecdf1 , ecdf2 )             
            ab  = ansari_bradley     ( ecdf1 , ecdf2 )
            st  = siegel_tukey       ( ecdf1 , ecdf2 )
            md  = mood               ( ecdf1 , ecdf2 )
            vdw = van_der_waerden    ( ecdf1 , ecdf2 )
            bm  = brunner_munzel     ( ecdf1 , ecdf2 )
            es  = epps_singleton     ( ecdf1 , ecdf2 , frequencies = self.frequencies )
            
            
            counters [ 'KS'  ] += ks 
            counters [ 'K'   ] += k 
            counters [ 'AD'  ] += ad
            counters [ 'CM'  ] += cm
            counters [ 'BJ'  ] += bj
            counters [ 'ZK'  ] += zk
            counters [ 'ZA'  ] += za
            counters [ 'ZC'  ] += zc
            counters [ 'W1'  ] += w1 
            counters [ 'W2'  ] += w2 
            counters [ 'R'   ] += r
            counters [ 'U2'  ] += u2
            counters [ 'AB'  ] += ab 
            counters [ 'ST'  ] += st  
            counters [ 'MD'  ] += md 
            counters [ 'VDW' ] += vdw 
            counters [ 'BM'  ] += bm 
            counters [ 'ES'  ] += es 
            
            results  [ 'KS'  ] .append ( ks  )    
            results  [ 'K'   ] .append ( k   ) 
            results  [ 'AD'  ] .append ( ad  ) 
            results  [ 'CM'  ] .append ( cm  ) 
            results  [ 'BJ'  ] .append ( bj  ) 
            results  [ 'ZK'  ] .append ( zk  ) 
            results  [ 'ZA'  ] .append ( za  ) 
            results  [ 'ZC'  ] .append ( zc  ) 
            results  [ 'W1'  ] .append ( w1  ) 
            results  [ 'W2'  ] .append ( w2  ) 
            results  [ 'R'   ] .append ( r   ) 
            results  [ 'U2'  ] .append ( u2  ) 
            results  [ 'AB'  ] .append ( ab  ) 
            results  [ 'ST'  ] .append ( st  ) 
            results  [ 'MD'  ] .append ( md  ) 
            results  [ 'VDW' ] .append ( vdw ) 
            results  [ 'BM'  ] .append ( bm  ) 
            results  [ 'ES'  ] .append ( es  ) 
            
            in_excess =   ( ref_ks  <= ks  or 
                            ref_k   <= k   or 
                            ref_ad  <= ad  or 
                            ref_cm  <= cm  or 
                            ref_bj  <= bj  or 
                            ref_zk  <= zk  or 
                            ref_za  <= za  or 
                            ref_zc  <= zc  or 
                            ref_w1  <= w1  or 
                            ref_w2  <= w2  or  
                            ref_r   <= r   or
                            ref_u2  <= u2  or  
                            ref_ab  <= ab  or 
                            ref_st  <= st  or 
                            ref_md  <= md  or 
                            ref_vdw <= vdw or 
                            ref_bm  <= bm  or 
                            ref_es  <= es  ) 
                                                                                                            
            if no_weight : 
                mw = mann_whitney  ( ecdf1 , ecdf2 )             
                s  = savage        ( ecdf1 , ecdf2 )     
                counters [ 'MW'  ] += mw
                counters [ 'S'   ] += s 
                results  [ 'MW'  ] .append ( mw ) 
                results  [ 'S'   ] .append ( s  ) 
                in_excess = in_excess or ref_mw <= mw or ref_s <= s
                
            self.__counter += in_excess 
            self.__nToys   += 1 
            
        ## accumulate number of to

        ##  ECDF = Ostap.Math.ECDF 
        for key in results :
            data = results [ key ]
            if not data : continue
            if not key in self.__ecdfs : self.__ecdfs [ key ] = ECDF ( data , True ) ## complementary ECDF!
            else                       : self.__ecdfs [ key ]  .add  ( data2vct ( data ) ) 
            
        del results

    # ========================================================================
    ## Helper method to get the result
    def result ( self , label ) :
        """ Helper method to get the result 
        """
        if not label in self.estimators : return None
        if not label in self.ecdfs      : return None
        if not label in self.counters   : return None
        ##
        value   = self.estimators   [ label ]
        ecdf    = self.ecdfs        [ label ] 
        counter = self.counters     [ label ] 
        ##
        pvalue  = ecdf. estimate ( value  ) ## estimate the p-value
        #
        pv      = clip_pvalue    ( pvalue , 0.5 ) 
        nsigma  = significance   ( pv           ) ## convert  it to significace
         
        return self.Result ( value   ,
                             counter ,
                             pvalue  ,
                             nsigma  )
    
    # =========================================================================
    ## format a row in the summary table
    def row  ( self , what , result , width = 6 , precision = 4 ) :
        """ Format a row in the sumamry table
        """
        value      = result.statistics
        counter    = result.counter
        pvalue     = result.pvalue
        nsigma     = result.nsigma
        
        mean       = counter.mean   ()
        rms        = counter.rms    () 
        vmin, vmax = counter.minmax () 
        
        mxv = max ( abs ( value        ) ,
                    abs ( mean.value() ) ,
                    mean.error()         , rms ,
                    abs ( vmin )  , abs ( vmax ) ) 
        
        fmt, fmtv , fmte , expo = fmt_pretty_ve ( VE ( mxv ,  mean.cov2() ) ,
                                                  width       = width       ,
                                                  precision   = precision   , 
                                                  parentheses = False       )
        
        if expo : scale = 10**expo
        else    : scale = 1
        
        fmt2 = '%s/%s' % ( fmtv , fmtv ) 

        vs  = value / scale
        vm  = mean  / scale
        vr  = rms   / scale
        vmn = vmin  / scale
        vmx = vmax  / scale
        
        pvalue = format_pvalue ( pvalue )
        nsigma = format_nsigma ( nsigma )
        
        return ( what  ,
                 fmtv  % vs ,
                 fmt   % ( vm.value() , vm.error() ) ,
                 fmtv  % vr                          ,
                 fmtv  % vmn                         ,
                 fmtv  % vmx                         ,
                 ( '[%s]' % format_pow10 ( expo ) )  if expo else '' , pvalue , nsigma )
    
    # =========================================================================
    ## Make a summary table
    def table ( self , title = '' , prefix = '' , width = 6 , precision = 4 , style = None ) :
        """ Make a summary table
        """
        import ostap.logger.table  as     T                 
        header = ( 'Statistics'  , ) + full_header
        rows   = [] 
        
        for label in self.ecdfs :
            
            result  = self.result ( label )
            if not result : continue

            the_label = Labels.get ( label , label )
            row = self.row ( the_label , result , width = width , precision = precision )
            rows.append ( row ) 
                    
        if   not title and self.nToys : title = 'Goodness of 1D-fit with #%d toys' % self.nToys  
        elif not title                : title = 'Goodness of 1D-fit'

        rows = [ header ] + sorted ( rows ) 
        
        ## ## the last row 
        ##  pvalue = self.counter.eff 
        ##  pv     = clip_pvalue   ( pvalue ) 
        ##  nsigma = significance  ( pv     ) 
        ##  pvalue = format_pvalue ( pvalue )
        ##  nsigma = format_nsigma ( nsigma ) 
        ##  last   = 'COMBINED' , '' , '' , '' ,'' , '' , '' , pvalue  , nsigma  
        ##  rows.append ( last )
        
        rows = T.remove_empty_columns ( rows ) 
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lcccccccccc' , style = style )

    __repr__ = table
    __str__  = table

    # =========================================================================
    ## Draw ECDF for toys & statistical estimator 
    def draw  ( self , what , option = '' , *options , **kwargs ) :
        """ Draw ECDF for toys & statistical estgimator 
        """
        key = cidict_fun ( what ) 
        if   key in Keys [ 'KS' ] and 'KS' in self.ecdfs :             
            result = self.result ( 'KS' )
            ecdf   = self.ecdfs  [ 'KS' ]
        elif key in Keys [ 'K'  ] and 'K'  in self.ecdfs : 
            result = self.result ( 'K' )
            ecdf   = self.ecdfs  [ 'K' ]
        elif key in Keys [ 'AD' ] and 'AD' in self.ecdfs :             
            result = self.result ( 'AD' )
            ecdf   = self.ecdfs  [ 'AD' ]
        elif key in Keys [ 'CM' ] and 'CM' in self.ecdfs : 
            result = self.result  ( 'CM' )
            ecdf   = self.ecdfs   [ 'CM' ]
        elif key in Keys [ 'BJ' ] and 'BJ' in self.ecdfs : 
            result = self.result  ( 'BJ' )
            ecdf   = self.ecdfs   [ 'BJ' ]            
        elif key in Keys [ 'ZK' ]  and 'ZK' in self.ecdfs : 
            result = self.result  ( 'ZK' )
            ecdf   = self.ecdfs   [ 'ZK' ]
        elif key in Keys [ 'ZA' ]   and 'ZA' in self.ecdfs :  
            result = self.result  ( 'ZA' )
            ecdf   = self.ecdfs   [ 'ZA' ]
        elif key in Keys [ 'ZC' ] and 'ZC' in self.ecdfs : 
            result = self.result  ( 'ZC' )
            ecdf   = self.ecdfs   [ 'ZC' ]
        elif key in Keys [ 'W1' ] and 'W1' in self.ecdfs : 
            result = self.result  ( 'W1' )
            ecdf   = self.ecdfs   [ 'W1' ]
        elif key in Keys [ 'W2' ] and 'W2' in self.ecdfs : 
            result = self.result  ( 'W2' )
            ecdf   = self.ecdfs   [ 'W2' ]
        elif key in Keys [ 'R' ] and 'R' in self.ecdfs : 
            result = self.result  ( 'R' )
            ecdf   = self.ecdfs   [ 'R' ]
        elif key in Keys [ 'U2' ] and 'U2' in self.ecdfs : 
            result = self.result  ( 'U2' )
            ecdf   = self.ecdfs   [ 'U2' ]
        elif key in Keys [ 'AB' ] and 'AB' in self.ecdfs : 
            result = self.result  ( 'AB' )
            ecdf   = self.ecdfs   [ 'AB' ]
        elif key in Keys [ 'ST' ] and 'ST' in self.ecdfs : 
            result = self.result  ( 'ST' )
            ecdf   = self.ecdfs   [ 'ST' ]
        elif key in Keys [ 'MD' ] and 'MD' in self.ecdfs : 
            result = self.result  ( 'MD' )
            ecdf   = self.ecdfs   [ 'MD' ] 
        elif key in Keys [ 'VDW' ] and 'VDW' in self.ecdfs : 
            result = self.result  ( 'VDW' )
            ecdf   = self.ecdfs   [ 'VDW' ]    
        elif key in Keys [ 'BM' ] and 'BM' in self.ecdfs : 
             result = self.result  ( 'BM' )
             ecdf   = self.ecdfs   [ 'BM' ] 
        elif key in Keys [ 'ES' ] and 'ES' in self.ecdfs : 
             result = self.result  ( 'ES' )
             ecdf   = self.ecdfs   [ 'ES' ]    
        elif key in Keys [ 'MW' ] and 'MW' in self.ecdfs : 
            result = self.result  ( 'MW' )
            ecdf   = self.ecdfs   [ 'MW' ]
        elif key in Keys [ 'S' ] and 'S' in self.ecdfs : 
            result = self.result  ( 'S' )
            ecdf   = self.ecdfs   [ 'S' ]
        else :
            raise KeyError (  "draw: Invalid `what`:%s" % what )

        xmin , xmax = ecdf.low_edge () , ecdf.high_edge ()
        value       = result.statistics
        xmin        = min ( xmin , value )
        xmax        = max ( xmax , value )
        xmin , xmax = axis_range ( xmin , xmax , delta = 0.20 )

        kwargs [ 'xmin' ] = kwargs.get ( 'xmin' , xmin ) 
        kwargs [ 'xmax' ] = kwargs.get ( 'xmax' , xmax )

        result    = ecdf.draw  ( option , *options , **kwargs ) 
        line      = ROOT.TLine ( value  , 1e-3 , value , 1 - 1e-3 )
        ## 
        line.SetLineWidth ( 4 ) 
        line.SetLineColor ( 8 ) 
        line.draw ( 'same' )
        ##
        self._line = line
        return result, line
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================




