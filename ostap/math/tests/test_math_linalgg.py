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
tolerance1 = 1.e-10
tolerance2 = 1.e-12 

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

    D      = L * U - P * A
    delta1 = Ostap.Math.maxabs_element ( D )
    
    logger.info ( '(P)LU max-difference %.3g :\n%s' % ( delta1 , D ) ) 

    assert delta1 < tolerance1 , '(P)LU: result are inconsistent delta1=%.3g' % delta1
    if not delta1 < tolerance2 : logger.error ( "delta1 is too large: %.3g"   % delta1 )

# ==========================================================================
def test_linalg_PQR ( M = 4 , N = 4 ) :
    
    logger = getLogger ( 'test_linalg_(P)QR(%s,%s)' % ( M , N )  )
    
    A = Matrix ( M , N )
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A.set ( i , j , i + j + random.gauss ( 1 , 1 ) )
            
    logger.info ( '(P)QR The matrix is:\n%s' % A )

    ## PLU decomposiiton 
    P, Q, R = A.PQR () 
    
    logger.info ( '(P)QR decomposition: P :\n%s' % P )
    logger.info ( '(P)QR decomposition: Q :\n%s' % Q )
    logger.info ( '(P)QR decomposition: R :\n%s' % R )

    D      = Q * R - A * P
    delta1 = Ostap.Math.maxabs_element ( D ) 
    logger.info ( '(P)QR max-difference %.3g :\n%s' % ( delta1 , D ) ) 

    QQ     = Q*Q.t()
    QQ    -= 1 
    delta2 = Ostap.Math.maxabs_element ( QQ )    
    logger.info ( '(P)QR non-orthogonality of Q %.3g \n%s' % ( delta2 , QQ ) ) 

    assert delta1 < tolerance1 , '(P)QR: result are inconsistent delta1=%.3g' % delta1
    assert delta2 < tolerance1 , '(P)QR: result are inconsistent delta2=%.3g' % delta2
    if not delta1 < tolerance2 : logger.error ( "delta1 is too large: %.3g"   % delta1 )
    if not delta2 < tolerance2 : logger.error ( "delta2 is too large: %.3g"   % delta2 )
    
# ==========================================================================
def test_linalg_LQ ( M = 4 , N = 4 ) :
    
    logger = getLogger ( 'test_linalg_LQ(%s,%s)' % ( M , N )  )
    
    A = Matrix ( M , N )
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A.set ( i , j , i + j + random.gauss ( 1 , 1 ) )
            
    logger.info ( 'LQ The matrix is:\n%s' % A )

    ## PLU decomposiiton 
    L, Q  = A.LQ () 
    
    logger.info ( 'LQ decomposition: L :\n%s' % L )
    logger.info ( 'LQ decomposition: Q :\n%s' % Q )
 
    D      = L * Q - A 
    delta1 = Ostap.Math.maxabs_element ( D ) 
    
    logger.info ( 'LQ max-difference %.3g :\n%s' % ( delta1 , D ) ) 

    QQ     = Q*Q.t()
    QQ    -= 1 
    delta2 = Ostap.Math.maxabs_element ( QQ )    
    logger.info ( 'LQ non-orthogonality of Q %.3g \n%s' % ( delta2 , QQ ) ) 
    
    assert delta1 < tolerance1 , 'LQ: result are inconsistent delta1=%.3g'  % delta1
    assert delta2 < tolerance1 , 'LQ: result are inconsistent delta2=%.3g'  % delta2
    if not delta1 < tolerance2 : logger.error ( "delta1 is too large: %.3g" % delta1 )
    if not delta2 < tolerance2 : logger.error ( "delta2 is too large: %.3g" % delta2 )
    
# ==========================================================================
def test_linalg_QL ( M = 4 , N = 4 ) :
    
    logger = getLogger ( 'test_linalg_QL(%s,%s)' % ( M , N )  )
    
    A = Matrix ( M , N )
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A.set ( i , j , i + j + random.gauss ( 1 , 1 ) )
            
    logger.info ( 'QL The matrix is:\n%s' % A )

    ## QL-decomposiiton 
    Q , L = A.QL () 
    
    logger.info ( 'QL decomposition: Q :\n%s' % Q )
    logger.info ( 'QL decomposition: L :\n%s' % L )

    D      = Q * L - A 
    delta1 = Ostap.Math.maxabs_element ( D ) 
    
    logger.info ( 'QL max-difference %.3g :\n%s' % ( delta1 , D ) ) 

    QQ     = Q*Q.t()
    QQ    -= 1 
    delta2 = Ostap.Math.maxabs_element ( QQ )    
    logger.info ( 'QL non-orthogonality of Q %.3g \n%s' % ( delta2 , QQ ) ) 

    assert delta1 < tolerance1 , 'QL: result are inconsistent delta1=%.3g'   % delta1
    assert delta2 < tolerance1 , 'QL: result are inconsistent delta2=%.3g'   % delta2
    if not delta1 < tolerance2 : logger.error ( "delta1 is too large: %.3g"  % delta1 )
    if not delta2 < tolerance2 : logger.error ( "delta2 is too large: %.3g"  % delta2 )
           
# ==========================================================================
def test_linalg_COD ( M = 4 , N = 4 ) :
    
    logger = getLogger ( 'test_linalg_COD(%s,%s)' % ( M , N )  )
    
    A = Matrix ( M , N )
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A.set ( i , j , i + j + random.gauss ( 1 , 1 ) )
            
    logger.info ( 'COD The matrix is:\n%s' % A )

    ## COD-decomposition 
    P, Q , R , Z  = A.COD () 
    
    logger.info ( 'COD decomposition: P :\n%s' % P )
    logger.info ( 'COD decomposition: Q :\n%s' % Q )
    logger.info ( 'COD decomposition: R :\n%s' % R )
    logger.info ( 'COD decomposition: Z :\n%s' % Z )
    
    D      = Q * R * Z.t() - A * P  
    delta1 = Ostap.Math.maxabs_element ( D ) 
    logger.info ( 'COD max-difference %.3g :\n%s' % ( delta1 , D ) ) 
    
    QQ      = Q*Q.t()
    QQ     -= 1 
    delta2  = Ostap.Math.maxabs_element ( QQ )    
    logger.info ( 'COD non-orthogonality of Q %.3g \n%s' % ( delta2 , QQ ) ) 
    
    ZZ      = Z*Z.t()
    ZZ     -= 1 
    delta3  = Ostap.Math.maxabs_element ( ZZ )    
    logger.info ( 'COD non-orthogonality of Z %.3g \n%s' % ( delta3 , ZZ ) ) 
    
    assert delta1 < tolerance1 , 'COD: result are inconsistent delta1=%.3g' % delta1
    assert delta2 < tolerance1 , 'COD: result are inconsistent delta2=%.3g' % delta2
    assert delta3 < tolerance1 , 'COD: result are inconsistent delta3=%.3g' % delta3 
    
    if not delta1 < tolerance2 : logger.error ( "delta1 is too large: %.3g" % delta1 )
    if not delta2 < tolerance2 : logger.error ( "delta2 is too large: %.3g" % delta2 )
    if not delta3 < tolerance2 : logger.error ( "delta3 is too large: %.3g" % delta3 )

# ==========================================================================
def test_linalg_SVD ( M = 4 , N = 4 ) :
    
    logger = getLogger ( 'test_linalg_SVD(%s,%s)' % ( M , N )  )
    
    A = Matrix ( M , N )
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A.set ( i , j , i + j + random.gauss ( 1 , 1 ) )
            
    logger.info ( 'SVD The matrix is:\n%s' % A )

    ## COD-decomposition 
    S , U , V = A.SVD () 
    
    logger.info ( 'SVD decomposition: S :\n%s' % S )
    logger.info ( 'SVD decomposition: U :\n%s' % U )
    logger.info ( 'SVD decomposition: V :\n%s' % V )
    
    D      = U * Matrix ( S ) * V.t() - A 
    delta1 = Ostap.Math.maxabs_element ( D ) 
    logger.info ( 'SVD max-difference %.3g :\n%s' % ( delta1 , D ) ) 

    
    UU      = U.t() * U if M >= N else U * U.t() 
    UU     -= 1 
    delta2  = Ostap.Math.maxabs_element ( UU )    
    logger.info ( 'SVD non-orthogonality of U %.3g \n%s' % ( delta2 , UU ) ) 

    VV      = V*V.t() if M >= N else V.t() * V 
    VV     -= 1 
    delta3  = Ostap.Math.maxabs_element ( VV )    
    logger.info ( 'SVD non-orthogonality of V %.3g \n%s' % ( delta3 , VV ) ) 
    
    assert delta1 < tolerance1 , 'SVD: result are inconsistent delta1=%.3g' % delta1
    assert delta2 < tolerance1 , 'SVD: result are inconsistent delta2=%.3g' % delta2
    assert delta3 < tolerance1 , 'SVD: result are inconsistent delta3=%.3g' % delta3 
    
    if not delta1 < tolerance2 : logger.error ( "delta1 is too large: %.3g" % delta1 )
    if not delta2 < tolerance2 : logger.error ( "delta2 is too large: %.3g" % delta2 )
    if not delta3 < tolerance2 : logger.error ( "delta3 is too large: %.3g" % delta3 )

# ==========================================================================
def test_linalg_POLAR( M = 4 ) :
    
    logger = getLogger ( 'test_linalg_POLAR(%s)' % ( M )  )
    
    A = Matrix ( M , M )
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A.set ( i , j , i + j + random.gauss ( 1 , 1 ) )
            
    logger.info ( 'POLAR The matrix is:\n%s' % A )

    ## COD-decomposition 
    U , P = A.POLAR() 
    
    logger.info ( 'POLAR decomposition: U :\n%s' % U )
    logger.info ( 'POLAR decomposition: P :\n%s' % P )
    
    D      = U * P - A 
    delta1 = Ostap.Math.maxabs_element ( D ) 
    logger.info ( 'POLAR max-difference %.3g :\n%s' % ( delta1 , D ) ) 
    
    UU      = U * U.t() 
    UU     -= 1 
    delta2  = Ostap.Math.maxabs_element ( UU )    
    logger.info ( 'POLAR non-orthogonality of U %.3g \n%s' % ( delta2 , UU ) ) 
    
    assert delta1 < tolerance1 , 'SVD: result are inconsistent delta1=%.3g' % delta1
    assert delta2 < tolerance1 , 'SVD: result are inconsistent delta2=%.3g' % delta2
    
    if not delta1 < tolerance2 : logger.error ( "delta1 is too large: %.3g" % delta1 )
    if not delta2 < tolerance2 : logger.error ( "delta2 is too large: %.3g" % delta2 )
    
# =============================================================================
if '__main__' == __name__ :

    """
    test_linalg_PLU   ( 3, 6 )
    test_linalg_PLU   ( 3, 3 )
    test_linalg_PLU   ( 6, 3 )

    test_linalg_PQR   ( 3, 6 )
    test_linalg_PQR   ( 3, 3 )
    test_linalg_PQR   ( 6, 3 )

    test_linalg_LQ    ( 3, 6 )
    test_linalg_LQ    ( 3, 3 )
    test_linalg_LQ    ( 6, 3 )

    test_linalg_QL    ( 3, 6 )
    test_linalg_QL    ( 3, 3 )
    test_linalg_QL    ( 6, 3 )
    
    test_linalg_COD   ( 3, 6 )
    test_linalg_COD   ( 3, 3 )
    test_linalg_COD   ( 6, 3 )
    """
    
    test_linalg_SVD   ( 3, 6 )
    test_linalg_SVD   ( 3, 3 )
    test_linalg_SVD   ( 6, 3 )

    test_linalg_POLAR ( 3 )
    test_linalg_POLAR ( 6 )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
