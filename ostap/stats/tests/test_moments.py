#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
""" Test module for ostap/stat/moments.py.
"""
# =============================================================================
import math
from   ostap.stats.moments import  ( Mean          ,
                                     Variance      ,
                                     Median        ,
                                     Mode          ,
                                     Width         ,
                                     RMS           ,
                                     Moment        ,
                                     CentralMoment ,
                                     Skewness      ,
                                     Kurtosis      ,
                                     Quantile      ,
                                     cl_symm       ,
                                     cl_asymm      ,
                                     skewness      ,
                                     quantile      ,
                                     kurtosis      )
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_moments' )
else                       : logger = getLogger ( __name__        )
# ============================================================================= 

# ============================================================================= 
def test_moments1() :
    #
    import math
    mean_ = Mean     (0, math.pi)
    logger.info ( 'sin@[0,pi]                mean: %s ' % mean_ (math.sin) ) 

    var2 = Variance (0, math.pi)
    logger.info ( 'sin@[0,pi]            variance: %s ' % var2 (math.sin) ) 

    rms_ = RMS     (0, math.pi)
    logger.info ( 'sin@[0,pi]                 rms: %s ' % rms_  (math.sin) ) 

    mom5 = Moment  ( 5, 0, math.pi)
    logger.info ( 'sin@[0,pi]          5th moment: %s ' % mom5 (math.sin) ) 

    mom5 = CentralMoment  ( 5, 0, math.pi)
    logger.info ( 'sin@[0,pi] 5th central moment : %s ' % mom5 (math.sin) ) 

    mom1 = Moment  ( 1, 0, math.pi)
    logger.info ( 'sin@[0,pi]          1st moment: %s ' % mom1 (math.sin) ) 

    mom1 = CentralMoment  ( 1, 0, math.pi)
    logger.info ( 'sin@[0,pi] 1st central moment : %s ' % mom1 (math.sin) ) 

    s    = Skewness ( 0 , math.pi )
    logger.info ( 'sin@[0,pi]            skewness: %s ' % s    (math.sin) ) 
    
    k    = Kurtosis ( 0 , math.pi )
    logger.info ( 'sin@[0,pi]            kurtosis: %s ' % k    (math.sin) ) 

    from math import exp 
    gau = lambda x : exp(-0.5*x*x)

    logger.info ( 'Skewness(gauss,-10,10) %s ' % skewness ( gau  , -10 , 10 ) ) 
    logger.info ( 'Kurtosis(gauss,-10,10) %s ' % kurtosis ( gau  , -10 , 10 ) )

    gau1 = lambda x : exp(-0.5*x*x) if x > 0 else 0

    logger.info ( 'Skewness(agau ,-10,10) %s ' % skewness ( gau1 , -10 , 10 ) ) 
    logger.info ( 'Kurtosis(agau ,-10,10) %s ' % kurtosis ( gau1 , -10 , 10 ) ) 

    logger.info ( 80*'*' ) 


# =============================================================================
def test_moments2 () :
    #
    import math

    try :
        from scipy.optimize import brentq
    except ImportError :
        logger.warning('Scipy.optimize.brentq is not availabe, skip test')
        return
        
    med  = Median   (0, math.pi)
    logger.info ( 'sin@[0,pi]              median: %s ' % med  (math.sin) ) 

    mode_ = Mode     (0, math.pi)
    logger.info ( 'sin@[0,pi]                mode: %s ' % mode_ (math.sin)  ) 

    def fwhm ( fun ) :
        _w    =   Width   (0, math.pi)
        x1,x2 = _w ( fun )
        return x2-x1
    logger.info ( 'sin@[0,pi]                fwhm: %s ' % fwhm (math.sin)  ) 

    quan = Quantile ( 0.980 , 0, 10)
    logger.info ( 'sin@[0,pi]      0.980-quantile: %s ' % quan (math.sin) ) 

    quan = Quantile ( 0.252 , 0, 10)
    logger.info ( '1@[0,10]        0.252-quantile: %s ' % quan ( lambda x : 1 ) ) 
    logger.info ( '1@[0,1]         0.501-quantile: %s ' % quantile ( lambda x : 1 , 0.501 , 0 , 1 ) )
    logger.info ( '1@[0,1]         0.201-quantile: %s ' % quantile ( lambda x : 1 , 0.201 , 0 , 1 ) ) 
    
    from math import exp 
    gau = lambda x : exp(-0.5*x*x)

    logger.info ( 'CL(gauss,0.68,-10,10)  %s ' % cl_symm ( gau , 0.68 , -10 , 10 ) ) 
    logger.info ( 'CL(gauss,0.68,0,10,0)  %s ' % cl_symm ( gau , 0.68 ,   0 , 10 , x0 = 0 ) ) 

    gau1 = lambda x : exp(-0.5*x*x) if x > 0 else 0

    logger.info ( 'CLa(gauss,0.68,-10,10) (%.3f,%.3f) ' % cl_asymm ( gau  , 0.68 , -10 , 10 ) ) 
    logger.info ( 'CLa(aga  ,0.68,-10,10) (%.3f,%.3f) ' % cl_asymm ( gau1 , 0.68 , -10 , 10 ) )
    
    logger.info ( 80*'*' ) 

# =============================================================================
if '__main__' == __name__ :

    test_moments1()
    test_moments2()
        
# =============================================================================
# The END 
# ============================================================================= 

