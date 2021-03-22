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
from __future__ import print_function
from   sys      import version_info as python_version
import math 
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'tests_math_linalg2'  )
else                       : logger = getLogger ( __name__              )
# ============================================================================= 
import ostap.math.linalg
from   ostap.core.core      import Ostap 
from   ostap.core.meta_info import root_version_int 
try :
    import numpy as np
except ImportError :
    np = None


# =============================================================================
## The main function to test linear algebra 
def test_linalg2() :
    """The main function to test linear algebra
    """
    
    logger.info('Test Linaear Algebra: ')
    
    LA3 = Ostap.Vector(3)
    l1  = LA3(0,1,2)
    l2  = LA3(3,4,5)
    
    logger.info ( 'l1 , l2 : %s %s '  % ( l1 , l2  ) )
    logger.info ( 'l1 + l2 : %s    '  % ( l1 + l2  ) )
    
    logger.info ( 'l1 - l2 : %s    '  % ( l1 - l2  ) )
    logger.info ( 'l1 * l2 : %s    '  % ( l1 * l2  ) )
    logger.info ( 'l1 *  2 : %s    '  % ( l1 *  2  ) )
    logger.info ( ' 2 * l2 : %s    '  % ( 2  * l2  ) )
    logger.info ( 'l1 /  2 : %s    '  % ( l1 /  2  ) )
    
    l1 /= 2 
    logger.info ( 'l1 /= 2 : %s    '  % l1 )
    
    l1 *= 2 
    logger.info ( 'l1 *= 2 : %s    '  % l1 )


    ## if ( 3 , 5 ) <= python_version :
        
    ##     logger.info ( 'l1 @ l2 : %s    '  % ( l1 @ l2  ) )
    ##     logger.info ( 'l1 @  2 : %s    '  % ( l1 @  2  ) )
    ##     logger.info ( ' 2 @ l2 : %s    '  % ( 2  @ l2  ) )
        
    logger.info('TEST matrices: ')
    
    m22 = Ostap.Math.Matrix(2,2) ()
    m23 = Ostap.Math.Matrix(2,3) ()
    s22 = Ostap.Math.SymMatrix(2)()
    
    l2  = Ostap.Math.Vector(2)()
    l3  = Ostap.Math.Vector(3)()
    
    l2[0]    = 1
    l2[1]    = 2
    
    l3[0]    = 1
    l3[1]    = 2
    l3[1]    = 3
    
    logger.info ( 'l2 , l3 : %s %s '  % ( l2 , l3  ) )
    
    m22[0,0] = 1
    m22[0,1] = 1
    m22[1,1] = 1
    
    m23[0,0] = 1
    m23[1,1] = 1
    m23[0,2] = 1
    
    s22[0,0] = 2
    s22[1,0] = 1
    s22[1,1] = 3
    
    logger.info ( 'm22\n%s'    % m22     ) 
    logger.info ( 's22\n%s'    % s22     ) 
    logger.info ( 'm23\n%s'    % m23     ) 
    logger.info ( 'm22/3\n%s'  % (m22/3) )
    
    logger.info ( 'm23*3\n%s'  % (m23*3) ) 

    logger.info ( 'm22**3\n%s' % m22**3  ) 
    logger.info ( 's22**4\n%s' % s22**4  ) 

    logger.info ( 'm22 * m23 :\n%s' % ( m22 * m23 ) ) 
    logger.info ( 'm22 *  l2 : %s ' % ( m22 * l2  ) ) 
    logger.info ( 'l2  * m22 : %s ' % ( l2  * m22 ) ) 
    logger.info ( 'm23 *  l3 : %s ' % ( m23 * l3  ) ) 
    logger.info ( 'l2  * m23 : %s ' % ( l2  * m23 ) )
    
    logger.info ( 'm22 * s22 + 2 * m22 :\n%s ' %  ( m22*s22 + 2*m22  ) )
    logger.info ( 'm22 == m22*1.0 : %s ' % (  m22 == m22 * 1.0 ) )
    logger.info ( 'm22 != m22*1.1 : %s ' % (  m22 != m22 * 1.1 ) )
    logger.info ( 'm23 == m23*1.0 : %s ' % (  m23 == m23 * 1.0 ) )
    logger.info ( 'm23 != m23*1.1 : %s ' % (  m23 != m23 * 1.1 ) )
    logger.info ( 'l1  == l1 *1.0 : %s ' % (  l1  == l1  * 1.0 ) )
    logger.info ( 'l1  != l1 *1.1 : %s ' % (  l1  != l1  * 1.1 ) )
    logger.info ( 's22 == s22*1.0 : %s ' % (  s22 == s22 * 1.0 ) )
    logger.info ( 's22 != s22*1.1 : %s ' % (  s22 != s22 * 1.1 ) )
    
    logger.info ( ' l1 == (0,1,2) : %s ' % (  l1 == ( 0 , 1 , 2 ) ) )
    logger.info ( ' l1 == [0,1,2] : %s ' % (  l1 == [ 0 , 1 , 2 ] ) )

    ## if ( 3 , 5 ) <= python_version :
        
    ##     logger.info ( 'm23 @ 3   :\n%s' % ( m23 @ 3   ) ) 
    ##     logger.info ( 'm22 @ m23 :\n%s' % ( m22 @ m23 ) ) 
    ##     logger.info ( 'm22 @  l2 : %s ' % ( m22 @ l2  ) ) 
    ##     logger.info ( 'm23 @  l3 : %s ' % ( m23 @ l3  ) ) 
         
    m22[0,0] = 1
    m22[0,1] = 2
    m22[1,0] = 2
    m22[1,1] = 3
    
    s22[0,0] = 1
    s22[0,1] = 2
    s22[1,1] = 3
    
    logger.info ( ' m22 == s22     : %s ' % ( m22 == s22       ) )
    logger.info ( ' m22 == s22*1.0 : %s ' % ( m22 == s22 * 1.0 ) )
    logger.info ( ' m22 != s22*1.1 : %s ' % ( m22 != s22 * 1.1 ) )

    ## ok 
    m22 + m22

    ## crash 
    m22 += m22

    ## crash
    m22 += Ostap.Math.Matrix(2,2) ()

    logger.info ( ' m22 += m22  :\n%s ' % m22 )

    m22 -= m22*2

    logger.info ( ' m22 += m22*2 :\n%s ' % m22 )


    m22 += s22*0
    m22 += s22
    m22 = m22 + s22


    logger.info ( ' m22 += s22*0 :\n%s ' % m22 )

    m22 -= s22*2
    logger.info ( ' m22 -= s22*2 :\n%s ' % m22 )

    s22 += s22*2
    logger.info ( ' s22 += s22*2 :\n%s ' % s22 )

    s22 -= s22*2
    logger.info ( ' s22 -= s22*2 :\n%s ' % s22 )
    
    if np :
        logger.info ( 'Operations with numpy objects')
        
        v2 = np.array ( [1.0,2.0]      )
        v3 = np.array ( [1.0,2.0,3.0 ] )
        
        logger.info ( 'v2  * l2  : %s' % ( v2  * l2  ) )
        logger.info ( 'l3  * v3  : %s' % ( l3  * v3  ) )
        logger.info ( 's22 * v2  : %s' % ( s22 * v2  ) )
        logger.info ( 'm22 * v2  : %s' % ( m22 * v2  ) )
        logger.info ( 'm23 * v3  : %s' % ( m23 * v3  ) )

        logger.info ( 'm22 as np : %s' % ( m22.to_numpy() ) )
        logger.info ( 's22 as np : %s' % ( s22.to_numpy() ) )
        logger.info ( 'm23 as np : %s' % ( m23.to_numpy() ) )
        
        if 62006 <= root_version_int :
            logger.warning ("Tests with numpy are broken for ROOT %s" %  root_version_int ) 
        else :
            
            logger.info ( 'm22  + m22(np) :\n%s' % ( m22 + m22.to_numpy () ) )
            logger.info ( 'm22  + s22(np) :\n%s' % ( m22 + s22.to_numpy () ) )
            logger.info ( 's22  + s22(np) :\n%s' % ( s22 + s22.to_numpy () ) )
            logger.info ( 's22  + m22(np) :\n%s' % ( s22 + s22.to_numpy () ) )
            
            logger.info ( 'm22  * m22(np) :\n%s' % ( m22 * m22.to_numpy () ) )
            logger.info ( 's22  * s22(np) :\n%s' % ( s22 * s22.to_numpy () ) )
            logger.info ( 's22  * m23(np) :\n%s' % ( s22 * m23.to_numpy () ) )        
            logger.info ( 'l2   * m22(np) :\n%s' % ( l2  * m22.to_numpy () ) )

        
    logger.info ( 'SVector with errors')

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

    test_linalg2 ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
