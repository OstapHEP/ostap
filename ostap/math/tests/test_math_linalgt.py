#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_linalgt.py
#  Test module for the file ostap/math/linalgt.py
# ============================================================================= 
""" Test module for ostap/math/linalgt.py
"""
# ============================================================================= 
from __future__ import print_function
from   sys      import version_info as python_version
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'tests_math_linalgt'  )
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
def test_linalgt() :
    """The main function to test linear algebra
    """
    
    logger.info('Test Linear Algebra: ')

    logger.info('TEST vectors: ')
    
    l1 = Ostap.TVector(3)
    l2 = Ostap.TVector(3)

    l1[0],l1[1],l1[2] = 0,1,2
    l2[0],l2[1],l2[2] = 3,4,5
    
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
    
    m22 = Ostap.Math.TMatrix(2,2)
    m23 = Ostap.Math.TMatrix(2,3) 
    s22 = Ostap.Math.TMatrixSym(2)
    
    l2  = Ostap.TVector(2)
    l3  = Ostap.TVector(3)
    
    l2[0]    = 1
    l2[1]    = 2
    
    l3[0]    = 1
    l3[1]    = 2
    l3[1]    = 3
    
    logger.info ( 'l2 , l3 : %s %s '  % ( l2 , l3  ) )

    
    ## if ( 3 , 5 ) <= python_version :
        
    ##     logger.info ( 'm23 @ 3   :\n%s' % ( m23 @ 3   ) ) 
    ##     logger.info ( 'm22 @ m23 :\n%s' % ( m22 @ m23 ) ) 
    ##     logger.info ( 'm22 @  l2 : %s ' % ( m22 @ l2  ) ) 
    ##     logger.info ( 'm23 @  l3 : %s ' % ( m23 @ l3  ) ) 
         

    m22[0,0] = 1
    m22[0,1] = 1
    m22[1,1] = 1
    
    m23[0,0] = 1
    m23[1,1] = 1
    m23[0,2] = 1
    
    s22[0,0] = 2
    s22[1,0] = 1
    s22[1,1] = 3
    
    logger.info ( 'm22\n%s'   % m22     ) 
    logger.info ( 's22\n%s'   % s22     ) 
    logger.info ( 'm23\n%s'   % m23     ) 
    logger.info ( 'm22/3\n%s' % (m22/3) ) 
    logger.info ( 'm23*3\n%s' % (m23*3) ) 

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

    m22 += m22*2
    m22 -= m22*1

    m22 += s22*2
    m22 -= s22*1

    s22 += s22*2
    s22 -= s22*1
    
    
    if np :
        logger.info ( 'Operations with numpy objects')
        
        v2 = np.array ( [1.0,2.0]      )
        v3 = np.array ( [1.0,2.0,3.0 ] )

        logger.info ( 'v2  * l2  : %s' % ( v2  * l2  ) )
        logger.info ( 'l3  * v3  : %s' % ( l3  * v3  ) )
        logger.info ( 's22 * v2  : %s' % ( s22 * v2  ) )
        logger.info ( 'm22 * v2  : %s' % ( m22 * v2  ) )
        logger.info ( 'm23 * v3  : %s' % ( m23 * v3  ) )
        

        n22_m = m22.to_numpy ()
        n22_s = s22.to_numpy ()
        n23   = m23.to_numpy ()
        
        if 62006 <= root_version_int :
            logger.warning ("Tests with numpy are broken for ROOT %s" %  root_version_int ) 
        else : 
            logger.info ( 'm22  * m22(np) :\n%s' % ( m22 * m22.to_numpy() ) )
            logger.info ( 's22  * s22(np) :\n%s' % ( s22 * s22.to_numpy() ) )
            logger.info ( 's22  * m23(np) :\n%s' % ( s22 * m23.to_numpy() ) )        
            logger.info ( 'l2   * m22(np) :\n%s' % ( l2  * m22.to_numpy() ) )


# =============================================================================
if '__main__' == __name__ :

    test_linalgt ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
