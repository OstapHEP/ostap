#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_covtransform.py
#  Test module for the file ostap/math/covtransfrom.py
# ============================================================================= 
""" Test module for ostap/math/covtransform.py
"""
# ============================================================================= 
from   ostap.core.core         import Ostap, VE
from   ostap.math.covtransform import transform
from   ostap.utils.utils       import batch_env 
import ostap.math.linalg 
import ROOT, math 
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'tests_math_covtranform'  )
else                       : logger = getLogger ( __name__                  )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
#


def test_covtransform1 ()  :

    logger.info ( 'Descartes to polar tranformation' )

    ## Descartes coordinates
    X = 1 , 4

    ##  their covarinac ematrix
    C = Ostap.SymMatrix(2)()
    C [ 0, 0 ] = 0.20
    C [ 0, 1 ] = 0.05
    C [ 1, 1 ] = 0.30
    
    R  = C.correlations()
    
    ##  Polar coordinated
    r   = lambda x,y : math.sqrt(x*x+y*y)
    phi = lambda x,y : math.atan2(y,x) 

    
    C_polar = transform ( C , X , r , phi  )

    R_polar = C_polar.correlations()

    logger.info ( 'Descartes Covariance  :\n%s' % C       )
    logger.info ( 'Descartes Correlation :\n%s' % R       )
    logger.info ( 'Polar     Covariance  :\n%s' % C_polar )
    logger.info ( 'Polar     Correlation :\n%s' % R_polar )


# =============================================================================
def test_covtransform2 ()  :

    logger.info ( 'Absolute yields to the ratios' )

    ##  three independent variables:
    X = VE ( 81143.24 , 286.58 ** 2  ) , \
        VE (  4231.37 ,  72.57 ** 2  ) , \
        VE (   137.11 ,  25.79 ** 2  )

    ## their ratios 
    
    alpha = lambda x0 , x1 , x2 : x2 * 1.0 / x1
    beta  = lambda x0 , x1 , x2 : x2 * 1.0 / x0
    gamma = lambda x0 , x1 , x2 : x1 * 1.0 / x0
    
    C = transform ( None , X , alpha , beta , gamma )

    R = C.correlations()

    logger.info ( 'Covariance  :\n%s' % C )
    logger.info ( 'Correlation :\n%s' % R )


# =============================================================================
def test_covtransform3 ()  :

    logger.info ( 'Absolute values  to the differences' )

    ##  three independent variables:
    X = VE ( 3686.05 , 0.01 ** 2  ) , \
        VE ( 3871.55 , 0.06 ** 2  ) , \
        VE ( 3824.04 , 0.53 ** 2  )

    ## their differences 
    
    alpha = lambda x0 , x1 , x2 : x1 - x2 
    beta  = lambda x0 , x1 , x2 : x2 - x0 
    gamma = lambda x0 , x1 , x2 : x1 - x0 
    
    C = transform ( None , X , alpha , beta , gamma )

    R = C.correlations()

    logger.info ( 'Covariance  :\n%s' % C )
    logger.info ( 'Correlation :\n%s' % R )

    
# =============================================================================
if '__main__' == __name__ :

    test_covtransform1 () 
    test_covtransform2 () 
    test_covtransform3 () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
