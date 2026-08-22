#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/twosamples2.py
#  Two-Samples 1D-Weighted Tests
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-10-16
# =============================================================================
""" Two-Samples 1D (weighted) Test
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
from   ostap.stats.utils    import ( weight_trivial     ,
                                     valid_weight       ,
                                     compatible_weights , 
                                     check_all          ,
                                     num_features       )
from   ostap.math.math_base import Ostap, data2vct 
from   ostap.stats.counters import ECDF, WECDF 
import numpy 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.twosamples2' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Two-Samples 1D (weighted) Test' )
# =============================================================================
ecdf_types = ECDF , WECDF
# =============================================================================
## Convert data/weight pair into (W)ECDF
#  @see Ostap::Math::ECDF
#  @see Ostap::Math::WECDF
def data2ecdf ( data , weight = None ) :
    """ Convert data/weight pair into (W)ECDF
    - see Ostap::Math::ECDF
    - see Ostap::Math::WECDF
    """
    ##
    if isinstance ( data , ecdf_types ) : return data
    ## 
    data      = numpy.asarray ( data , dtype = float ).ravel() 
    if 1 != num_features      ( data   )          : raise TypeError  ( "Invalid `data` type!"                   )
    if not valid_weight       ( weight )          : raise ValueError ( "Invalid `weight` type/value!"           )
    if not compatible_weights ( data   , weight ) : raise TypeError  ( "Incompatible `data/weight` structures!" )
    ##
    if weight_trivial ( weight ) : return ECDF ( data2vct ( data ) )
    ## 
    weight = numpy.asarray ( weight , dtype = float ).ravel() 
    return WECDF ( data2vct ( data ) , data2vct ( weight ) ) 

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
    ## 
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
    ## 
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
    ## 
    return Ostap.Math.cramer_von_mises  ( ecdf1 , ecdf2 )

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
    ##
    ecdf1 = data2ecdf ( data1 , weight1 )
    ecdf2 = data2ecdf ( data2 , weight2 )
    ## 
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
    ## 
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
    ## 
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
    ## 
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
    ## 
    return Ostap.Math.ZC ( ecdf1 , ecdf2 , eps )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
