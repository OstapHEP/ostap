#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_linalg2.py
#  Test module for the file ostap/math/linalg2.py
# ============================================================================= 
""" Test module for ostap/math/linalg2.py
"""
# =============================================================================
from   ostap.math.base    import Ostap
from   ostap.math.linalgg import Matrix
from   ostap.utils.utils  import batch_env
import math, random  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'tests_math_linalgg'  )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

def test_linalg_PLU ( M = 4 , N = 4 ) :
    
    logger = getLogger ( 'test_linalg_(P)LU(%s,%s)' % ( M , N )  )
    
    A = Matrix ( M , N )
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A.set ( i , j , i + j + random.gauss ( 1 , 1 ) )
            
    logger.info ( '(P)LU The matrix is:\n%s' % A )

    ## PLU decomposiiton 
    P, L, U = A.PLU () 
    
    logger.info ( '(P)LU decomposition: P :\n%s' % P )
    logger.info ( '(P)LU decomposition: L :\n%s' % L )
    logger.info ( '(P)LU decomposition: U :\n%s' % U )

    D     = L * U - P * A
    delta = Ostap.Math.maxabs_element ( D ) 

    logger.info ( '(P)LU max-difference %.3g :\n%s' % ( delta , D ) ) 

    
# =============================================================================
if '__main__' == __name__ :
    
    test_linalg_PLU ( 3, 6 )
    test_linalg_PLU ( 3, 3 )
    test_linalg_PLU ( 6, 3 )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
