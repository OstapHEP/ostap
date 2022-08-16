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
from __future__             import print_function
from   sys                  import version_info as python_version
import math 
from   ostap.math.linalg    import checkops 
from   ostap.core.core      import Ostap 
from   ostap.core.meta_info import root_version_int 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'tests_math_linalg2'  )
else                       : logger = getLogger ( __name__              )
# ============================================================================= 
try :
    import numpy as np
except ImportError :
    np = None


# =============================================================================
def test_linalg2_vct () :

    logger = getLogger ( "test_linalg2_vct" )
    
    logger.info ( 'Basic operations with vectors' )

    LA3 = Ostap.Vector(3)
    l1  = LA3 ( 0 , 1 , 2 )
    l2  = LA3 ( 3 , 4 , 5 )

    
    logger.info ( 'l1 ,  l2    : %s ,  %s      '    % ( l1 , l2  ) )
    logger.info ( 'l1 +  l2    : %s +  %s = %s '    % ( l1 , l2 , l1 + l2  ) )
    logger.info ( 'l1 -  l2    : %s -  %s = %s '    % ( l1 , l2 , l1 - l2  ) )
    logger.info ( 'l1 *   2    : %s *  %s = %s '    % ( l1 ,  2 , l1 *  2  ) )
    logger.info ( 'l1 /   3    : %s /  %s = %s '    % ( l1 ,  3 , l1 *  3  ) )
    logger.info ( ' 2 *  l1    : %s *  %s = %s '    % (  2 , l1 ,  2 *  l1 ) )

    logger.info ( 'l1 *  l2    : %s *  %s = %s '    % ( l1 , l2 , l1 * l2  ) )

    logger.info ( 'l1 x  l2    : %s x  %s = \n%s '  % ( l1 , l2 , l1.cross ( l2 ) ) )
    logger.info ( 'l2 x  l1    : %s x  %s = \n%s '  % ( l2 , l1 , l2.cross ( l1 ) ) )

    logger.info ( 'l1 == l2    : %s == %s = %s '    % ( l1 , l2 , l1 == l2  ) )
    logger.info ( 'l1 != l2    : %s != %s = %s '    % ( l1 , l2 , l1 != l2  ) )
    logger.info ( 'l1 == l1    : %s == %s = %s '    % ( l1 , l1 , l1 == l1  ) )
    logger.info ( 'l1 != l1    : %s != %s = %s '    % ( l1 , l1 , l1 != l1  ) )

    v1 = ( 0 , 1 , 2 ) 
    logger.info ( 'l1 == tuple : %s == %s = %s '    % (  l1 , v1 , l1 == v1 ) )
    logger.info ( 'tuple == l1 : %s == %s = %s '    % (  v1 , l1 , v1 == l1 ) )
    v1 = [ 0 , 1 , 2 ] 
    logger.info ( 'l1 == list  : %s == %s = %s '    % (  l1 , v1 , l1 == v1 ) )
    logger.info ( 'list == l1  : %s == %s = %s '    % (  v1 , l1 , v1 == l1 ) )

    l1 *= 2 
    logger.info ( 'l1 *= 2     : %s '  % l1 )

    l1 /= 2 
    logger.info ( 'l1 /= 2     : %s '  % l1 )
    

    ## if ( 3 , 5 ) <= python_version :
    
    ##     logger.info ( 'l1 @ l2 : %s    '  % ( l1 @ l2  ) )
    ##     logger.info ( 'l1 @  2 : %s    '  % ( l1 @  2  ) )
    ##     logger.info ( ' 2 @ l2 : %s    '  % ( 2  @ l2  ) )

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
    
    l2  = Ostap.Math.Vector(2)( 1 , 2 )
    l3  = Ostap.Math.Vector(3)( 1 , 2 , 3 )
        
    logger.info ( 'initial    l2 , l3 : %s , %s '  % ( l2 , l3  ) )
    
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


    logger.info ( 'power      m22**3\n%s'  % m22**3  ) 
    logger.info ( 'power      s22**4\n%s'  % s22**4  )
    logger.info ( 'power      m22**-1\n%s' % m22**-1 ) 
    logger.info ( 'power      m22**-3\n%s' % m22**-3 ) 
    logger.info ( 'power      (m22**-3)*(m22**3)\n%s' % ( (m22**-3)*(m22**3) ) ) 
    
    logger.info ( 'multiply   m22 * m23 :\n%s' % ( m22 * m23 ) ) 
    logger.info ( 'multiply   m22 *  l2 : %s ' % ( m22 * l2  ) ) 
    logger.info ( 'multiply   s22 *  l2 : %s ' % ( s22 * l2  ) ) 
    logger.info ( 'multiply   l2  * m22 : %s ' % ( l2  * m22 ) ) 
    logger.info ( 'multiply   m23 *  l3 : %s ' % ( m23 * l3  ) ) 
    logger.info ( 'multiply   l2  * m23 : %s ' % ( l2  * m23 ) )
    logger.info ( 'multiply   l2  * s22 : %s ' % ( l2  * s22 ) )
    
    logger.info ( 'expression m22 * s22 + 2 * m22 :\n%s ' % ( m22*s22 + 2*m22  ) )
    logger.info ( 'equality   m22 == m22*1.0 : %s '       % ( m22 == m22 * 1.0 ) )
    
    logger.info ( 'equality   m22 != m22*1.1 : %s ' % (  m22 != m22 * 1.1 ) )
    logger.info ( 'equality   m23 == m23*1.0 : %s ' % (  m23 == m23 * 1.0 ) )
    logger.info ( 'equality   m23 != m23*1.1 : %s ' % (  m23 != m23 * 1.1 ) )
    logger.info ( 'equality   s22 == s22*1.0 : %s ' % (  s22 == s22 * 1.0 ) )
    logger.info ( 'equality   s22 != s22*1.1 : %s ' % (  s22 != s22 * 1.1 ) )
    
    ## if ( 3 , 5 ) <= python_version :
        
    ##     logger.info ( 'm23 @ 3   :\n%s' % ( m23 @ 3   ) ) 
    ##     logger.info ( 'm22 @ m23 :\n%s' % ( m22 @ m23 ) ) 
    ##     logger.info ( 'm22 @  l2 : %s ' % ( m22 @ l2  ) ) 
    ##     logger.info ( 'm23 @  l3 : %s ' % ( m23 @ l3  ) ) 

    logger.info ( 'sim        s22.sim(l2)   : %s  ' % (  s22.sim ( l2  ) ) )
    logger.info ( 'sim        s22.sim(s22)  :\n%s ' % (  s22.sim ( s22 ) ) )
    logger.info ( 'sim        s22.sim(m22)  :\n%s ' % (  s22.sim ( m22 ) ) )
    logger.info ( 'simT       s22.simT(m23) :\n%s ' % (  s22.simT( m23 ) ) )



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
    
    if not np :
        logger.warning  ( 'No Numpy, test is disabled ')
        return
    
    LA2 = Ostap.Vector(2)
    vla = LA2 ( 1 , 2 )
    vnp = np.array ( [1.0 , 2.0] )

    logger.info    ( 'Scalar product of two vectors (dot)' ) 
    logger.info    ( 'la  * np  : %s *  %s = %s' % ( vla , vnp , vla  * vnp ) )
    logger.info    ( 'np  + la  : %s +  %s = %s' % ( vnp , vla , vnp  + vla ) )
    logger.info    ( 'la  + np  : %s +  %s = %s' % ( vla , vnp , vla  + vnp ) )
    logger.info    ( 'np  - la  : %s -  %s = %s' % ( vnp , vla , vnp  - vla ) )
    logger.info    ( 'la  - np  : %s -  %s = %s' % ( vla , vnp , vla  - vnp ) )
    
    logger.info    ( 'la == np  : %s == %s = %s' % ( vla , vnp , vla == vnp ) )
    
    vnp += vla
    logger.info    ( 'np += la  : %s' % ( vnp ) ) 
    vnp -= vla
    logger.info    ( 'np -= la  : %s' % ( vnp ) ) 

    vla += vnp
    logger.info    ( 'la += np  : %s' % ( vla ) ) 
    vla -= vnp
    logger.info    ( 'la -= np  : %s' % ( vla ) ) 


    logger.warning ('Wrong results with np as a left argument:' ) 
    logger.warning ( 'np  * la  : %s *  %s = %s' % ( vnp , vla , vnp  * vla ) )
    logger.warning ( 'np == la  : %s == %s = %s' % ( vnp , vla , vnp == vla ) )


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
    
    logger.info    ( 's22 * la  : %s' % ( s22 * vla  ) )
    logger.info    ( 's22 * np  : %s' % ( s22 * vnp  ) )
    logger.info    ( 'la * s22  : %s' % ( vla * s22  ) )
    
    logger.info    ( 'm22 * la  : %s' % ( m22 * vla  ) )    
    logger.info    ( 'm22 * np  : %s' % ( m22 * vnp  ) )    
    logger.info    ( 'la * m22  : %s' % ( vla * m22  ) )

    logger.warning ('Wrong results with np as left argument:' ) 
    logger.warning ( 'np * s22  :\n%s' % ( vnp * s22  ) )
    logger.warning ( 'np * m22  :\n%s' % ( vnp * m22  ) )

    logger.info    ( 's22 as np :\n%s' % ( s22.to_numpy() ) )
    logger.info    ( 'm22 as np :\n%s' % ( m22.to_numpy() ) )
    logger.info    ( 'm23 as np :\n%s' % ( m23.to_numpy() ) )

    
    logger.info    ( 'm22  + m22(np)    :\n%s' % ( m22 + m22.to_numpy () ) )
    logger.info    ( 'm22  + s22(np)    :\n%s' % ( m22 + s22.to_numpy () ) )
    logger.info    ( 's22  + s22(np)    :\n%s' % ( s22 + s22.to_numpy () ) )
    logger.info    ( 's22  + m22(np)    :\n%s' % ( s22 + s22.to_numpy () ) )
    
    logger.info    ( 'm22  * m22(np)    :\n%s' % ( m22 * m22.to_numpy () ) )
    logger.info    ( 's22  * s22(np)    :\n%s' % ( s22 * s22.to_numpy () ) )
    logger.info    ( 's22  * m23(np)    :\n%s' % ( s22 * m23.to_numpy () ) )
    
    logger.info    ( 'la   * m22(np)    : %s'  % ( vla * m22.to_numpy () ) )
    logger.info    ( 'la   * s22(np)    : %s'  % ( vla * s22.to_numpy () ) )
    logger.info    ( 'la   * m23(np)    : %s'  % ( vla * m23.to_numpy () ) )

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
    logger.info ( " -> rho, phi %s " % r1 )

    r2 = v2.transform ( rho  )
    logger.info ( " -> rho      %s " % r2 )
    
    
# =============================================================================
if '__main__' == __name__ :

    test_linalg2_vct  ()
    test_linalg2_mtrx ()
    test_linalg2_np   ()
    test_linalg2_ve   ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
