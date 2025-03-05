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
from   ostap.math.linalg    import checkops, gsl_info 
from   ostap.core.core      import Ostap
from   ostap.math.base      import numpy 
from   ostap.utils.utils    import batch_env 
import math, random  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'tests_math_linalg2'  )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
tolerance1 = 1.e-10
tolerance2 = 1.e-12 

# =============================================================================
def test_linalg2_vct () :

    logger = getLogger ( "test_linalg2_vct" )
    
    logger.info ( 'Basic operations with vectors' )

    LA3 = Ostap.Vector(3)
    l1  = LA3 ( 0 , 1 , 2 )
    l2  = LA3 ( 3 , 4 , 5 )
    
    logger.info ( 'l1 ,  l2    : \n%s  \n%s      '  % ( l1 , l2  ) )
    logger.info ( 'l1 +  l2    : \n%s '             % ( l1 + l2  ) )
    logger.info ( 'l1 -  l2    : \n%s '             % ( l1 - l2  ) )
    logger.info ( 'l1 *   2    : \n%s '             % ( l1 *  2  ) )
    logger.info ( 'l1 /   3    : \n%s '             % ( l1 /  3  ) )
    logger.info ( ' 2 *  l1    : \n%s '             % ( 2  *  l1 ) )

    logger.info ( 'l1 *  l2    : %s '               % ( l1 * l2  ) )

    logger.info ( 'l1 x  l2    : \n%s '             % ( l1.cross ( l2 ) ) )
    logger.info ( 'l2 x  l1    : \n%s '             % ( l2.cross ( l1 ) ) )

    logger.info ( 'l1 == l2    : %s '               % ( l1 == l2  ) )
    logger.info ( 'l1 != l2    : %s '               % ( l1 != l2  ) )
    logger.info ( 'l1 == l1    : %s '               % ( l1 == l1  ) )
    logger.info ( 'l1 != l1    : %s '               % ( l1 != l1  ) )

    v1 = ( 0 , 1 , 2 ) 
    logger.info ( 'l1 == tuple : %s '               % (  l1 == v1 ) )
    logger.info ( 'tuple == l1 : %s '               % (  v1 == l1 ) )
    v1 = [ 0 , 1 , 2 ] 
    logger.info ( 'l1 == list  : %s '               % (  l1 == v1 ) )
    logger.info ( 'list == l1  : %s '               % (  v1 == l1 ) )

    l1 *= 2 
    logger.info ( 'l1 *= 2     : \n%s '  % l1 )

    l1 /= 2 
    logger.info ( 'l1 /= 2     : \n%s '  % l1 )
    

    logger.info ( 'l1 @ l2 : %s    '  % ( l1 @ l2  ) )
    logger.info ( 'l1 @  2 : %s    '  % ( l1 @  2  ) )
    logger.info ( ' 2 @ l2 : %s    '  % ( 2  @ l2  ) )

    ##summary table for allowed  binary operations
    checkops ( l1 , l2 , logger = logger ) 

# =============================================================================         
def test_linalg2_mtrx () :
        
    logger = getLogger ( "test_linalg2_mtrx" )

    logger.info ( 'Basic operations with matrices' )
    
    logger.info('TEST matrices')
    
    m22 = Ostap.Math.Matrix(2,2)  ()
    m23 = Ostap.Math.Matrix(2,3)  ()
    s22 = Ostap.Math.SymMatrix(2) ()
    s33 = Ostap.Math.SymMatrix(3) ()
    
    l2  = Ostap.Math.Vector(2)( 1 , 2 )
    l3  = Ostap.Math.Vector(3)( 1 , 2 , 3 )
        
    logger.info ( 'initial    l2 , l3 : \n%s \n%s '  % ( l2 , l3  ) )
    
    m22 [ 0 , 0 ] = 1
    m22 [ 0 , 1 ] = 1
    m22 [ 1 , 1 ] = 1
    
    m23 [ 0 , 0 ] = 1
    m23 [ 1 , 1 ] = 1
    m23 [ 0 , 2 ] = 1
    
    s22 [ 0 , 0 ] = 2
    s22 [ 1 , 0 ] = 1
    s22 [ 1 , 1 ] = 3
    
    logger.info ( 'initial    m22\n%s'    % m22     ) 
    logger.info ( 'initial    s22\n%s'    % s22     ) 
    logger.info ( 'initial    m23\n%s'    % m23     )

    
    logger.info ( 'scale      m22/3\n%s'  % (m22/3) )    
    logger.info ( 'scale      m23*3\n%s'  % (m23*3) ) 
    logger.info ( 'scale      3*m23\n%s'  % (3*m23) )
    
    logger.info ( 'plus       m22+s22\n%s'  % (m22+s22) )
    logger.info ( 'plus       s22+m22\n%s'  % (s22+s22) )
    logger.info ( 'minus      s22-m22\n%s'  % (s22-m22) )
    logger.info ( 'minus      m22-s22\n%s'  % (m22-s22) )

    logger.info ( 'add        m22+2\n%s'  % (m22+2) )
    logger.info ( 'add        2+m22\n%s'  % (2+m22) )
    logger.info ( 'add        2-m22\n%s'  % (2-m22) )
    logger.info ( 'add        s22+2\n%s'  % (s22+2) )
    logger.info ( 'add        2+s22\n%s'  % (2+s22) )
    logger.info ( 'add        2-s22\n%s'  % (2-s22) )

    s22 += 2 
    logger.info ( 'iadd       s22+=2\n%s' % (s22) )
    s22 -= 2 
    logger.info ( 'isub       s22+=2\n%s' % (s22) )
    
    m22 += 2 
    logger.info ( 'iadd       m22+=2\n%s' % (m22) )
    m22 -= 2 
    logger.info ( 'isub       m22+=2\n%s' % (m22) )

    m22 += s22
    logger.info ( 'iadd       m22+=s22\n%s' % (m22) )
    m22 -= s22 
    logger.info ( 'isub       m22-=s22\n%s' % (m22) )


    logger.info ( 'power      m22**3\n%s'  % m22 ** 3 ) 
    logger.info ( 'power      s22**4\n%s'  % s22 ** 4 )
    logger.info ( 'power      m22**-1\n%s' % m22 **-1 ) 
    logger.info ( 'power      m22**-3\n%s' % m22 **-3 ) 
    logger.info ( 'power      m22** 1\n%s' % m22 ** 1 ) 
    logger.info ( 'power      m22** 0\n%s' % m22 ** 0 ) 
    logger.info ( 'power      (m22**-3)*(m22**3)\n%s' % ( (m22**-3)*(m22**3) ) ) 
    
    logger.info ( 'multiply   m22 * m23 :\n%s' % ( m22 * m23 ) ) 
    logger.info ( 'multiply   m22 *  l2 :\n%s ' % ( m22 * l2  ) ) 
    logger.info ( 'multiply   s22 *  l2 :\n%s ' % ( s22 * l2  ) ) 
    logger.info ( 'multiply   l2  * m22 :\n%s ' % ( l2  * m22 ) ) 
    logger.info ( 'multiply   m23 *  l3 :\n%s ' % ( m23 * l3  ) ) 
    logger.info ( 'multiply   l2  * m23 :\n%s ' % ( l2  * m23 ) )
    logger.info ( 'multiply   l2  * s22 :\n%s ' % ( l2  * s22 ) )
    
    logger.info ( 'expression m22 * s22 + 2 * m22 :\n%s ' % ( m22*s22 + 2*m22  ) )
    logger.info ( 'equality   m22 == m22*1.0 : %s '       % ( m22 == m22 * 1.0 ) )
    
    logger.info ( 'equality   m22 != m22*1.1 : %s ' % (  m22 != m22 * 1.1 ) )
    logger.info ( 'equality   m23 == m23*1.0 : %s ' % (  m23 == m23 * 1.0 ) )
    logger.info ( 'equality   m23 != m23*1.1 : %s ' % (  m23 != m23 * 1.1 ) )
    logger.info ( 'equality   s22 == s22*1.0 : %s ' % (  s22 == s22 * 1.0 ) )
    logger.info ( 'equality   s22 != s22*1.1 : %s ' % (  s22 != s22 * 1.1 ) )

        
    logger.info ( 'm23 @ 3   :\n%s' % ( m23 @ 3   ) ) 
    logger.info ( 'm22 @ m23 :\n%s' % ( m22 @ m23 ) ) 
    logger.info ( 'm22 @  l2 : %s ' % ( m22 @ l2  ) ) 
    logger.info ( 'm23 @  l3 : %s ' % ( m23 @ l3  ) ) 

    logger.info ( 'sim        s22.sim(l2)   : %s  ' % (  s22.sim ( l2  ) ) )
    logger.info ( 'sim        s22.sim(s22)  :\n%s ' % (  s22.sim ( s22 ) ) )
    logger.info ( 'sim        s22.sim(m22)  :\n%s ' % (  s22.sim ( m22 ) ) )
    logger.info ( 'simT       s22.simT(m23) :\n%s ' % (  s22.simT( m23 ) ) )

    ## eigen values and eigen vectors 
    values, vectors = s22.eigenVectors( sorted = True , ascending = True )
    logger.info ( 'Eigen values  \n%s'           % values  )
    logger.info ( 'Eigen vectors (columns) \n%s' % vectors )
    t22 = s22.simT ( vectors )
    logger.info ( ' s22.simT (V) (must be diagonal) \n%s' % t22 )
    d22 = t22 - t22.diagonal () 
    ln  = d22.lnorm ()
    mn  = d22.mnorm ()    
    logger.info ( ' s22.simT (V) (off-diagonal) %+.5g/%+.5g\n%s'  % ( ln , mn , d22 ) ) 

        
    ## eigenvalues and eigenvectors
    for e , v in s22.eigenitems ( sorted = True , ascending = False ) :
        logger.info ( 'Eigenvalue:%+.5g '     % e  )
        logger.info ( ' - Eigenvector    |v|= %+.5g |s22*v|=%+.5g\n%s ' % ( abs ( v ) , abs ( s22*v) , v ) )
        logger.info ( ' - check: |s22*v-e*v|= %+.5g' % abs ( s22 * v - e * v ) ) 


    s33 [ 0 , 0 ] = 1
    s33 [ 1 , 1 ] = 10
    s33 [ 2 , 2 ] = 100

    s33 [ 0 , 1 ] = 2 
    s33 [ 0 , 2 ] = 20 

    s33 [ 1 , 2 ] = 50 

    logger.info ( 'Matrix 3x3\n%s' % s33 )
    
    ## eigen values and eigen vectors 
    values, vectors = s33.eigenVectors( sorted = True , ascending = True )
    logger.info ( 'Eigen values  \n%s'                    % values  )
    logger.info ( 'Eigen vectors (columns) \n%s'          % vectors )
    t33 = s33.simT ( vectors )
    logger.info ( ' s33.simT (V) (must be diagonal) \n%s' % t33 )
    d33 = t33 - t33.diagonal ()
    ln  = d33.lnorm ()
    mn  = d33.mnorm ()    
    logger.info ( ' s33.simT (V) (off-diagonal) %+.5g/%+.5g\n%s'  % ( ln , mn , d33 ) ) 
    
    ## eigenvalues and eigenvectors
    for e , v in s33.eigenitems ( sorted = True , ascending = False ) :
        logger.info ( 'Eigenvalue:%+.5g '     % e  )
        logger.info ( ' - Eigenvector    |v|= %+.5g |s33*v|=%+.5g\n%s ' % ( abs ( v ) , abs ( s33 * v ) , v ) )
        logger.info ( ' - check: |s33*v-e*v|= %+.5g' % abs ( s33 * v - e * v ) ) 

    return 


## powers of the nilponent matrix
    N   = 4 
    nNN = Ostap.Math.Matrix(N,N) ()
    for i in range ( N ) :
        for j in range ( i + 1 , N ) :
            nNN [ i,j ] = 1
    logger.info ( 'Powers for nilpotent matrix:\ns' % nNN ) 
    for i in range ( N + 1 ) : 
        logger.info ( '   %2d power of matrix is:\n%s' % ( i , nNN ** i ) ) 
    

    logger.info ( 'simT       s22.simT(m23) :\n%s ' % (  s22.simT( m23 ) ) )


    logger.info ( 'Equality with a number' ) 
    nn  = Ostap.Math.Matrix(3,3) ()
    logger.info ( 'nn      :\n%s ' % (  nn  ) )
    nn += 2
    logger.info ( 'nn += 2 :\n%s ' % (  nn       ) )
    logger.info ( 'equality   nn == 2       : %s '  % (  nn == 2  ) )
    logger.info ( 'equality   nn != 2       : %s '  % (  nn != 2  ) )
    logger.info ( 'equality   2  == nn      : %s '  % (  2  == nn ) )
    logger.info ( 'equality   2 != nn       : %s '  % (  2  != nn ) )
    
    
    
# =============================================================================
def test_linalg2_np () :

    logger = getLogger("test_linalg2_np")
    
    logger.info ( 'Test combined operations with numpy objects')
    
    if not numpy :
        logger.warning  ( 'No Numpy, test is disabled ')
        return
    
    LA2 = Ostap.Vector(2)
    vla = LA2 ( 1 , 2 )
    vnp = numpy.array ( [1.0 , 2.0] )

    logger.info    ( 'Scalar product of two vectors (dot)' ) 
    logger.info    ( 'la  * np  : \n%s' % ( vla  * vnp ) )
    logger.info    ( 'np  + la  : \n%s' % ( vnp  + vla ) )
    logger.info    ( 'la  + np  : \n%s' % ( vla  + vnp ) )
    logger.info    ( 'np  - la  : \n%s' % ( vnp  - vla ) )
    logger.info    ( 'la  - np  : \n%s' % ( vla  - vnp ) )
    
    logger.info    ( 'la == np  : %s' % ( vla == vnp ) )
    
    vnp += vla
    logger.info    ( 'np += la  : \n%s' % ( vnp ) ) 
    vnp -= vla
    logger.info    ( 'np -= la  :\n%s' % ( vnp ) ) 

    vla += vnp
    logger.info    ( 'la += np  :\n%s' % ( vla ) ) 
    vla -= vnp
    logger.info    ( 'la -= np  :\n%s' % ( vla ) ) 


    logger.warning ('Wrong results with np as a left argument:' ) 
    logger.warning ( 'np  * la  : %s' % ( vnp  * vla ) )
    logger.warning ( 'np == la  : %s' % ( vnp == vla ) )


    m22 = Ostap.Math.Matrix(2,2)  ()
    m23 = Ostap.Math.Matrix(2,3)  ()
    s22 = Ostap.Math.SymMatrix(2) ()

    m22 [ 0 , 0 ] = 1
    m22 [ 0 , 1 ] = 1
    m22 [ 1 , 1 ] = 1
    
    m23 [ 0 , 0 ] = 1
    m23 [ 1 , 1 ] = 1
    m23 [ 0 , 2 ] = 1
    
    s22 [ 0 , 0 ] = 2
    s22 [ 1 , 0 ] = 1
    s22 [ 1 , 1 ] = 3
    
    logger.info    ( 's22 * la  : \n%s' % ( s22 * vla  ) )
    logger.info    ( 's22 * np  : \n%s' % ( s22 * vnp  ) )
    logger.info    ( 'la * s22  : \n%s' % ( vla * s22  ) )
    
    logger.info    ( 'm22 * la  : \n%s' % ( m22 * vla  ) )    
    logger.info    ( 'm22 * np  : \n%s' % ( m22 * vnp  ) )    
    logger.info    ( 'la * m22  : \n%s' % ( vla * m22  ) )

    logger.warning ('Wrong results with np as left argument:' ) 
    logger.warning ( 'np * s22  : \n%s' % ( vnp * s22  ) )
    logger.warning ( 'np * m22  : \n%s' % ( vnp * m22  ) )

    logger.info    ( 's22 as np : \n%s' % ( s22.to_numpy() ) )
    logger.info    ( 'm22 as np : \n%s' % ( m22.to_numpy() ) )
    logger.info    ( 'm23 as np : \n%s' % ( m23.to_numpy() ) )

    
    logger.info    ( 'm22  + m22(np)    :\n%s' % ( m22 + m22.to_numpy () ) )
    logger.info    ( 'm22  + s22(np)    :\n%s' % ( m22 + s22.to_numpy () ) )
    logger.info    ( 's22  + s22(np)    :\n%s' % ( s22 + s22.to_numpy () ) )
    logger.info    ( 's22  + m22(np)    :\n%s' % ( s22 + s22.to_numpy () ) )
    
    logger.info    ( 'm22  * m22(np)    :\n%s' % ( m22 * m22.to_numpy () ) )
    logger.info    ( 's22  * s22(np)    :\n%s' % ( s22 * s22.to_numpy () ) )
    logger.info    ( 's22  * m23(np)    :\n%s' % ( s22 * m23.to_numpy () ) )
    
    logger.info    ( 'la   * m22(np)    :\n%s'  % ( vla * m22.to_numpy () ) )
    logger.info    ( 'la   * s22(np)    :\n%s'  % ( vla * s22.to_numpy () ) )
    logger.info    ( 'la   * m23(np)    :\n%s'  % ( vla * m23.to_numpy () ) )

    m22 += s22 
    logger.info    ( 'm22 += s22        :\n%s' % ( m22 ) )
    m22 -= s22 
    logger.info    ( 'm22 -= s22        :\n%s' % ( m22 ) )
    m22 += s22.to_numpy()  
    logger.info    ( 'm22 += s22 (np)   :\n%s' % ( m22 ) )
    m22 -= s22.to_numpy()  
    logger.info    ( 'm22 -= s22 (np)   :\n%s' % ( m22 ) )

    logger.info    ( 's22.sim(la)       : %s  ' % (  s22.sim ( vla             ) ) )
    logger.info    ( 's22.sim(la)  (np) : %s  ' % (  s22.sim ( vla.to_numpy () ) ) )
    
    logger.info    ( 's22.sim(s22)      :\n%s ' % (  s22.sim ( s22             ) ) )
    logger.info    ( 's22.sim(s22) (np) :\n%s ' % (  s22.sim ( s22.to_numpy () ) ) )
    logger.info    ( 's22.sim(m22)      :\n%s ' % (  s22.sim ( m22             ) ) )
    logger.info    ( 's22.sim(m22) (np) :\n%s ' % (  s22.sim ( m22.to_numpy () ) ) )
    logger.info    ( 's22.simT(m23)     :\n%s ' % (  s22.simT( m23             ) ) )
    logger.info    ( 's22.simT(m23)(np) :\n%s ' % (  s22.simT( m23.to_numpy () ) ) )

# =============================================================================
def test_linalg2_ve () :

    logger = getLogger ( "test_linalg2_ve" ) 
    logger.info ( 'Test SVectorWithErrors')

    v2  = Ostap.Math.VectorE (2)()

    v2 [ 0 ] = 3
    v2 [ 1 ] = 4
    
    v2 . cov2 () [ 0 , 0 ] = 0.10
    v2 . cov2 () [ 0 , 1 ] = 0.05
    v2 . cov2 () [ 1 , 1 ] = 0.20

    rho = lambda x,y : ( x * x + y * y ) **  0.5
    phi = lambda x,y : math.atan2 ( y , x ) 
    

    r1 = v2.transform ( rho , phi )
    logger.info ( " -> rho, phi \n%s " % r1 )

    r2 = v2.transform ( rho  )
    logger.info ( " -> rho      \n%s " % r2 )
    

# ==================================================================================    
def test_linalg2_PLU ( M = 4 , N = 4 ) :
    
    logger = getLogger ( 'test_linalg_(P)LU(%s,%s)' % ( M , N )  )
    
    A = Ostap.Math.Matrix ( M , N )()
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A[ i , j ] = i + j + random.gauss ( 1 , 1 ) 
            
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
    
    
# =========================================================================================
def test_linalg2_PQR ( M = 4 , N = 4 ) :
    
    logger = getLogger ( 'test_linalg2_(P)QR(%s,%s)' % ( M , N )  )
    
    A = Ostap.Math.Matrix ( M , N )()
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A [ i , j ] = i + j + random.gauss ( 1 , 1 ) 
            
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
def test_linalg2_LQ ( M = 4 , N = 4 ) :
    
    logger = getLogger ( 'test_linalg2_LQ(%s,%s)' % ( M , N )  )
    
    A = Ostap.Math.Matrix ( M , N )()
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A [ i , j ] = i + j + random.gauss ( 1 , 1 ) 
            
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
def test_linalg2_QL ( M = 4 , N = 4 ) :
    
    logger = getLogger ( 'test_linalg2_QL(%s,%s)' % ( M , N )  )
    
    if gsl_info < ( 2 , 7 ) :
        logger.info ( 'Test is disbaled for GSL<2.7')
        return 
    
    A = Ostap.Math.Matrix ( M , N )()
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A [ i , j ] = i + j + random.gauss ( 1 , 1 ) 

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


## ==========================================================================
def test_linalg2_COD ( M = 4 , N = 4 ) :
    
    logger = getLogger ( 'test_linalg2_COD(%s,%s)' % ( M , N )  )
    
    A = Ostap.Math.Matrix ( M , N )()
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A [ i , j ] = i + j + random.gauss ( 1 , 1 ) 
            
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
def test_linalg2_SVD ( M = 4 , N = 4 ) :
    
    logger = getLogger ( 'test_linalg2_SVD(%s,%s)' % ( M , N )  )
    
    A = Ostap.Math.Matrix ( M , N )()
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A [ i , j ] = i + j + random.gauss ( 1 , 1 ) 
            
    logger.info ( 'SVD The matrix is:\n%s' % A )

    ## COD-decomposition 
    S , U , V = A.SVD () 
    
    logger.info ( 'SVD decomposition: S :\n%s' % S )
    logger.info ( 'SVD decomposition: U :\n%s' % U )
    logger.info ( 'SVD decomposition: V :\n%s' % V )
    
    D      = U * S * V.t() - A 
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
def test_linalg2_POLAR( M = 4 ) :
    
    logger = getLogger ( 'test_linalg2_POLAR(%s)' % ( M )  )
    
    A = Ostap.Math.Matrix ( M , M )()
    for i in range ( A.kRows ) :
        for j in range ( A.kCols ) :
            A [ i , j ] = i + j + random.gauss ( 1 , 1 ) 
            
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
    test_linalg2_vct  ()
    test_linalg2_mtrx ()
    test_linalg2_np   ()
    test_linalg2_ve   ()
    """
    
    test_linalg2_PLU   ( 3 , 6 )
    test_linalg2_PLU   ( 3 , 3 )
    test_linalg2_PLU   ( 6 , 3 )

    test_linalg2_PQR   ( 3 , 6 )
    test_linalg2_PQR   ( 3 , 3 )
    test_linalg2_PQR   ( 6 , 3 )

    test_linalg2_LQ    ( 3 , 6 )
    test_linalg2_LQ    ( 3 , 3 )
    test_linalg2_LQ    ( 6 , 3 )

    test_linalg2_QL    ( 3 , 6 )
    test_linalg2_QL    ( 3 , 3 )
    test_linalg2_QL    ( 6 , 3 )

    test_linalg2_COD   ( 3 , 6 )
    test_linalg2_COD   ( 3 , 3 )
    test_linalg2_COD   ( 6 , 3 )

    test_linalg2_SVD   ( 3 , 6 )
    test_linalg2_SVD   ( 3 , 3 )
    test_linalg2_SVD   ( 6 , 3 )

    test_linalg2_POLAR ( 3 )
    test_linalg2_POLAR ( 6 )

# =============================================================================
##                                                                      The END 
# =============================================================================
