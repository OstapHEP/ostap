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
from   sys                  import version_info as python_version
from   ostap.core.meta_info import root_info
from   ostap.math.linalg    import checkops 
from   ostap.core.core      import Ostap
from   ostap.utils.basic    import typename
from   ostap.utils.utils    import batch_env 
import ostap.logger.table   as     T 
import ROOT, array, random   
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'tests_math_linalgt'  )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
try : # =======================================================================
    import numpy as np
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    np = None

# =============================================================================
def test_linalgt_vct () :

    logger = getLogger ( "test_linalgt_vct" )
    
    logger.info ( 'Basic operations with vectors' )

    LA3 = Ostap.Vector(3)
    
    l1  = LA3 ( 0 , 1 , 2 )
    l2  = LA3 ( 3 , 4 , 5 )
    
    t1  = ROOT.TVectorD(3)
    t2  = ROOT.TVectorD(3)
    
    t1 [ 0 ] , t1 [ 1 ] , t1 [ 2]  = 0 , 1 , 2
    t2 [ 0 ] , t2 [ 1 ] , t2 [ 2 ] = 3 , 4 , 5

    
    logger.info ( 'l1 ,  l2    : \n%s \n%s      '    % ( l1 , l2  ) )
    logger.info ( 't1 ,  t2    : \n%s \n%s      '    % ( t1 , t2  ) )
    
    logger.info ( 'l1 +  l2    : \n%s '    % ( l1 + l2  ) )
    logger.info ( 'l1 +  t2    : \n%s '    % ( l1 + t2  ) )
    logger.info ( 't1 +  l2    : \n%s '    % ( t1 + l2  ) )
    logger.info ( 't1 +  t2    : \n%s '    % ( t1 + t2  ) )

    logger.info ( 'l1 -  l2    : \n%s '    % ( l1 - l2  ) )
    logger.info ( 'l1 -  t2    : \n%s '    % ( l1 - t2  ) )
    logger.info ( 't1 -  l2    : \n%s '    % ( t1 - l2  ) )
    logger.info ( 't1 -  t2    : \n%s '    % ( t1 - t2  ) )

    
    logger.info ( 'l1 *   2    : \n%s '    % ( l1 *  2  ) )
    logger.info ( 't1 *   2    : \n%s '    % ( t1 *  2  ) )
    
    logger.info ( 'l1 /   3    : \n%s '    % ( l1 *  3  ) )
    logger.info ( 't1 /   3    : \n%s '    % ( t1 *  3  ) )

    
    logger.info ( ' 2 *  l1    : \n%s '    % (  2 *  l1 ) )
    logger.info ( ' 2 *  t1    : \n%s '    % (  2 *  t1 ) )

    
    logger.info ( 'l1 *  l2    : \n%s '    % ( l1 * l2  ) )
    logger.info ( 'l1 *  t2    : \n%s '    % ( l1 * t2  ) )
    logger.info ( 't1 *  l2    : \n%s '    % ( t1 * l2  ) )
    logger.info ( 't1 *  t2    : \n%s '    % ( t1 * t2  ) )

    logger.info ( 'l1 x  l2    : \n%s '  % ( l1.cross ( l2 ) ) )
    logger.info ( 'l1 x  t2    : \n%s '  % ( l1.cross ( t2 ) ) )
    logger.info ( 't1 x  l2    : \n%s '  % ( t1.cross ( l2 ) ) )
    logger.info ( 't1 x  t2    : \n%s '  % ( t1.cross ( t2 ) ) )


    logger.info ( 'l1 == l2    : %s '    % ( l1 == l2  ) )
    logger.info ( 'l1 == t2    : %s '    % ( l1 == t2  ) )
    logger.info ( 't1 == l2    : %s '    % ( t1 == l2  ) )
    logger.info ( 't1 == t2    : %s '    % ( t1 == t2  ) )

    logger.info ( 'l1 != l2    : %s '    % ( l1 != l2  ) )
    logger.info ( 'l1 != t2    : %s '    % ( l1 != t2  ) )
    logger.info ( 't1 != l2    : %s '    % ( t1 != l2  ) )
    logger.info ( 't1 != t2    : %s '    % ( t1 != t2  ) )

    
    v1 = ( 0 , 1 , 2 ) 
    logger.info ( 'l1 == tuple : %s '    % (  l1 == v1 ) )
    logger.info ( 't1 == tuple : %s '    % (  t1 == v1 ) )
    logger.info ( 'tuple == l1 : %s '    % (  v1 == l1 ) )
    logger.info ( 'tuple == t1 : %s '    % (  v1 == t1 ) )

    v1 = [ 0 , 1 , 2 ] 
    logger.info ( 'l1 == list  : %s '    % (  l1 == v1 ) )
    logger.info ( 't1 == list  : %s '    % (  t1 == v1 ) )    
    logger.info ( 'list == l1  : %s '    % (  v1 == l1 ) )
    logger.info ( 'list == t1  : %s '    % (  v1 == t1 ) )

    l1 *= 2 
    t1 *= 2 
    logger.info ( 'l1 *= 2     : \n%s '  % l1 )
    logger.info ( 't1 *= 2     : \n%s '  % t1 )

    l1 /= 2 
    t1 /= 2 
    logger.info ( 'l1 /= 2     : \n%s '  % l1 )
    logger.info ( 't1 /= 2     : \n%s '  % t1 )
    

    ## if ( 3 , 5 ) <= python_version :
    
    ##     logger.info ( 'l1 @ l2 : %s    '  % ( l1 @ l2  ) )
    ##     logger.info ( 'l1 @  2 : %s    '  % ( l1 @  2  ) )
    ##     logger.info ( ' 2 @ l2 : %s    '  % ( 2  @ l2  ) )

    ##summary table for allowed  binary operations
    checkops ( l1 , l2 , logger = logger ) 
    checkops ( l1 , t2 , logger = logger ) 
    checkops ( t1 , l2 , logger = logger ) 
    checkops ( t1 , t2 , logger = logger ) 



# =============================================================================         
def test_linalgt_mtrx () :
        
    logger = getLogger ( "test_linalgt_mtrx" )

    logger.info ( 'Basic operations with matrices' )
    
    rows = [ ( 'Type1' , 'Op' , 'Type2' , 'Result' ) ]


    m22 = Ostap.Math.Matrix(2,2)  ()
    m23 = Ostap.Math.Matrix(2,3)  ()
    s22 = Ostap.Math.SymMatrix(2) ()

    t22 = ROOT.TMatrixD    ( 2 , 2 )
    t23 = ROOT.TMatrixD    ( 2 , 3 )
    r22 = ROOT.TMatrixDSym ( 2     ) 
    
    l2  = Ostap.Math.Vector(2)( 1 , 2 )
    l3  = Ostap.Math.Vector(3)( 1 , 2 , 3 )
    t2  = ROOT.TVectorD ( 2 )
    t3  = ROOT.TVectorD ( 3 )

    vectors2    = ( l2 , t2 )
    vectors3    = ( l3 , t3 )
    vectorsa    = vectors2  + vectors3 
    matrices22  = s22 , m22 , r22 , t22
    matrices23  = m23 , t23 
    symmatrices = s22 , r22
    
    t2 [ 0 ] , t2 [ 1 ]   = 1 , 2
    t3 [ 0 ] , t3 [ 1 ] , t3 [ 2 ] = 1 , 2 , 3 
    
    logger.info ( 'initial    l2 , l3 : \n%s \n%s '  % ( l2 , l3  ) )
    logger.info ( 'initial    t2 , t3 : \n%s \n%s '  % ( t2 , t3  ) )
    
    m22 [ 0 , 0 ] = 1
    m22 [ 0 , 1 ] = 1
    m22 [ 1 , 1 ] = 1
    
    m23 [ 0 , 0 ] = 1
    m23 [ 1 , 1 ] = 1
    m23 [ 0 , 2 ] = 1
    
    s22 [ 0 , 0 ] = 2
    s22 [ 1 , 0 ] = 1
    s22 [ 1 , 1 ] = 3

    t22 [ 0 , 0 ] = 1
    t22 [ 0 , 1 ] = 1
    t22 [ 1 , 1 ] = 1
    
    t23 [ 0 , 0 ] = 1
    t23 [ 1 , 1 ] = 1
    t23 [ 0 , 2 ] = 1
    
    r22 [ 0 , 0 ] = 2
    r22 [ 1 , 0 ] = 1
    r22 [ 1 , 1 ] = 3
    
    logger.info ( 'initial    m22\n%s'    % m22     ) 
    logger.info ( 'initial    s22\n%s'    % s22     ) 
    logger.info ( 'initial    m23\n%s'    % m23     )

    logger.info ( 'initial    t22\n%s'    % t22     ) 
    logger.info ( 'initial    r22\n%s'    % r22     ) 
    logger.info ( 'initial    t23\n%s'    % t23     )


    logger.info ( 'scale      m22/3\n%s'  % (m22/3) )    
    logger.info ( 'scale      m23*3\n%s'  % (m23*3) ) 
    logger.info ( 'scale      3*m23\n%s'  % (3*m23) )

    logger.info ( 'scale      t22/3\n%s'  % (t22/3) )    
    logger.info ( 'scale      t23*3\n%s'  % (t23*3) ) 
    logger.info ( 'scale      3*t23\n%s'  % (3*t23) )


    s22 += 2 
    logger.info ( 'iadd       s22+=2\n%s'   % (s22) )
    s22 -= 2 
    logger.info ( 'isub       s22+=2\n%s'   % (s22) )
    
    m22 += 2 
    logger.info ( 'iadd       m22+=2\n%s'   % (m22) )
    m22 -= 2 
    logger.info ( 'isub       m22+=2\n%s'   % (m22) )

    t22 += 2 
    logger.info ( 'iadd       t22+=2\n%s'   % (t22) )
    t22 -= 2 
    logger.info ( 'isub       t22+=2\n%s'   % (t22) )

    r22 += 2 
    logger.info ( 'iadd       r22+=2\n%s'   % (r22) )
    r22 -= 2 
    logger.info ( 'isub       r22+=2\n%s'   % (r22) )

    m22 += s22
    logger.info ( 'iadd       m22+=s22\n%s' % (m22) )
    m22 -= s22 
    logger.info ( 'isub       m22-=s22\n%s' % (m22) )

    m22 += r22
    logger.info ( 'iadd       m22+=r22\n%s' % (m22) )
    m22 -= r22 
    logger.info ( 'isub       m22-=r22\n%s' % (m22) )

    m22 += t22
    logger.info ( 'iadd       m22+=t22\n%s' % (m22) )
    m22 -= t22 
    logger.info ( 'isub       m22-=t22\n%s' % (m22) )

    t22 += r22
    logger.info ( 'iadd       t22+=r22\n%s' % (t22) )
    t22 -= r22 
    logger.info ( 'isub       r22-=r22\n%s' % (t22) )

    t22 += s22
    logger.info ( 'iadd       t22+=s22\n%s' % (t22) )
    t22 -= s22 
    logger.info ( 'isub       r22-=s22\n%s' % (t22) )

    t22 += m22
    logger.info ( 'iadd       t22+=m22\n%s' % (t22) )
    t22 -= m22 
    logger.info ( 'isub       r22-=m22\n%s' % (t22) )

    ## the basic operations 
    for v1 in vectors2 :
        n1 = typename ( v1 ) 
        for v2 in vectors2 :
            n2 = typename ( v1 ) 

            r = v1 + v2
            row     = n1 , '+' , n2 , typename ( r )
            rows.append ( row )
            
            r = v1 - v2
            row     = n1 , '-' , n2 , typename ( r )
            rows.append ( row )
            
            v1 += v2
            row     = n1 , '+=' , n2 , typename ( v1 )
            rows.append ( row )
            
            v1 -= v2
            row     = n1 , '-=' , n2 , typename ( v1 )
            rows.append ( row )

            r = v1 * v2
            row     = n1 , '*' , n2 , typename ( r )
            rows.append ( row )

    ## cross-product (tensor product) 
    for v1 in vectorsa :
        n1 = typename ( v1 ) 
        for v2 in vectorsa :
            n2 = typename ( v2 ) 
            r = v1.cross  ( v1 )
            row     = n1 , 'x' , n2 , typename ( r )
            rows.append ( row )
            logger.info ( ' %-25s x   %25s \n%s' % ( n1 , n2 , r ) )

    ## matrix-vector multiplications
    for m in matrices22 :
        n1 = typename ( m )
        for v in vectors2 :
            n2 = typename ( v  ) 
            r   = m * v
            logger.info ( ' %-25s +   %25s = \n%s' % ( n1 , n2 , r ) )
            row = n1 , '*' , n2 , typename ( r )
            rows.append ( row )

    ## vector-matrix multiplications
    for v in vectors2 :
        n1 = typename ( v  ) 
        for m in matrices22 :
            n2 = typename ( m )
            r   = v * m 
            logger.info ( ' %-25s +   %25s = \n%s' % ( n1 , n2 , r ) )
            row = n1 ,  '*' , n2 ,  typename ( r )
            rows.append ( row )

    ## matrix addition
    for m1 in matrices22 + ( 1.0 , ) :
        n1 = typename ( m1 ) 
        for m2 in matrices22 + ( 1.0 , ) :
            n2 = typename ( m2 )        
            if isinstance ( m1 , float ) and isinstance ( m2 , float  ) : continue 
            
            r   = m1 + m2
            logger.info ( ' %-25s +   %25s \n%s' % ( n1 , n2 , r ) )
            row = n1 , '+' , n2 , typename ( r )
            rows.append ( row )

    ## matrix subtraction 
    for m1 in matrices22 + ( 1.0 , )  :
        n1 = typename ( m1 ) 
        for m2 in matrices22 + ( 1.0 , ) :
            n2 = typename ( m2 ) 
            if isinstance ( m1 , float ) and isinstance ( m2 , float  ) : continue 

            r   = m1 - m2
            logger.info ( ' %-25s -  %25s \n%s' % ( n1 , n2 , r ) )
            row = n1 , '-' , n2 , typename ( r )
            rows.append ( row )

    ## matrix addtion/subtraction in place 
    for m1 in matrices22 :
        n1 = typename ( m1 ) 
        for m2 in matrices22 + ( 1.0, ):
            n2 = typename ( m2 ) 

            m1 += m2 
            logger.info ( ' %-25s +=  %25s \n%s' % ( n1 , n2 , m1 ) )
            row = n1 , '+=' , n2 , typename ( m1 )
            rows.append ( row )

            m1 -= m2 
            logger.info ( ' %-25s -=  %25s \n%s' % ( n1 , n2 , m1 ) )
            row = n1 , '-=' , n2 , typename ( m1 )
            rows.append ( row )

    ## matrix-matrix multiplication 
    for m1 in matrices22 :
        n1 = typename ( m1 ) 
        for m2 in matrices22 + matrices23 :
            n2 = typename ( m2 )
            
            r_mult   = m1 * m2
            logger.info ( ' %-25s *   %25s \n%s' % ( n1 , n2 , r_mult  ) )
            
            row3 = n1 , '*' , n2 , typename ( r_mult   )
            rows.append ( row3 )
            
    ## pow for matrices 
    for m1 in matrices22 :

        n1 = typename ( m1 )
        n2 = typename ( 3  )

        print ( 'm1:' , m1 )
        
        r1 = m1 **  3
        row = n1 , '**' , n2,  typename ( r1 )
        rows.append ( row )

        mm = m1 + 1 
        for d in range  ( 3 ) :
            
            logger.info ( 'power  M**%d\n%s'          % (  d , mm** d ) )
            logger.info ( 'power  M**%d\n%s'          % ( -d , mm**-d ) )      
            logger.info ( 'power  (M**%d)*(M**%d)\n%s' % ( d , -d , (mm**d)*(mm**-d)) ) 
                                        
    ## sim for matrices 
    for m1 in ( s22 , r22 ) :
        n1 = typename ( m1 ) 
        for m2 in matrices22 + ( l2 , t2 ) :
            n2 = typename ( m2 )             
            r_sim   = m1 .sim (  m2 ) 
            logger.info ( ' %-25s sim %25s = \n%s' % ( n1 , n2 , r_sim  ) )
            row     = n1 , 'sim' , n2 , typename ( r_sim )
            rows.append ( row )

    ## simT for matrices 
    for m1 in ( s22 , r22 ) :
        n1 = typename ( m1 ) 
        for m2 in matrices22 + matrices23 : 
            n2 = typename ( m2 )             
            r_sim   = m1 .simT (  m2 ) 
            logger.info ( ' %-25s simT %25s = \n%s' % ( n1 , n2 , r_sim  ) )
            row     = n1 , 'simT' , n2 , typename ( r_sim )
            rows.append ( row )
            
        
    title = 'Binary operations'
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lcll' )
    logger.info ( '%s:\n%s' % ( title , table ) ) 
        
            

# =============================================================================
## The main function to test linear algebra 
def test_linalgt_old () :
    """ The main function to test linear algebra
    """

    logger = getLogger( 'test_linangt_old')
    
    logger.info('Test Linear Algebra: ')

    logger.info('TEST vectors: ')
    
    l1 = Ostap.TVector(3)
    l2 = Ostap.TVector(3)

    l1[0],l1[1],l1[2] = 0,1,2
    l2[0],l2[1],l2[2] = 3,4,5
    
    logger.info ( 'l1 , l2 : \n%s \n%s '  % ( l1 , l2  ) )
    logger.info ( 'l1 + l2 : \n%s    '    % ( l1 + l2  ) )
    
    logger.info ( 'l1 - l2 : \n%s    '  % ( l1 - l2  ) )
    logger.info ( 'l1 * l2 : \n%s    '  % ( l1 * l2  ) )
    logger.info ( 'l1 *  2 : \n%s    '  % ( l1 *  2  ) )
    logger.info ( ' 2 * l2 : \n%s    '  % ( 2  * l2  ) )
    logger.info ( 'l1 /  2 : \n%s    '  % ( l1 /  2  ) )
    
    l1 /= 2 
    logger.info ( 'l1 /= 2 : \n%s    '  % l1 )
    l1 *= 2 
    logger.info ( 'l1 *= 2 : \n%s    '  % l1 )

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
    
    logger.info ( 'l2 , l3 : \n%s \n%s '  % ( l2 , l3  ) )

    
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

    logger.info ( 'm22 * m23 : \n%s' % ( m22 * m23 ) ) 
    logger.info ( 'm22 *  l2 : \n%s ' % ( m22 * l2  ) ) 
    logger.info ( 'l2  * m22 : \n%s ' % ( l2  * m22 ) ) 
    logger.info ( 'm23 *  l3 : \n%s ' % ( m23 * l3  ) ) 
    logger.info ( 'l2  * m23 : \n%s ' % ( l2  * m23 ) )
    
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
    
    ## DISABLE!!!
    if np : ## and False :

        logger.info ( 'Operations with numpy objects')
        
        v2 = np.array ( [1.0,2.0]      )
        v3 = np.array ( [1.0,2.0,3.0 ] )

        logger.info ( 'v2  * l2  : \n%s' % ( v2  * l2  ) )
        logger.info ( 'l3  * v3  : \n%s' % ( l3  * v3  ) )
        logger.info ( 's22 * v2  : \n%s' % ( s22 * v2  ) )
        logger.info ( 'm22 * v2  : \n%s' % ( m22 * v2  ) )
        logger.info ( 'm23 * v3  : \n%s' % ( m23 * v3  ) )
        

        n22_m = m22.to_numpy ()
        n22_s = s22.to_numpy ()
        n23   = m23.to_numpy ()
        
        if ( 6 , 20 ) <= root_info < ( 6 ,33 ) : 
            logger.warning ("Tests with numpy are broken for ROOT %s" %  str ( root_info ) ) 
        else :
            
            logger.info ( 'm22  * m22(np) :\n%s' % ( m22 * m22.to_numpy() ) )
            logger.info ( 's22  * s22(np) :\n%s' % ( s22 * s22.to_numpy() ) )
            logger.info ( 's22  * m23(np) :\n%s' % ( s22 * m23.to_numpy() ) )        
            logger.info ( 'l2   * m22(np) :\n%s' % ( l2  * m22.to_numpy() ) )



# ==============================================================================
## test to check TMatrix <-> array interplay 
def test_linalgt_arr () :
    """ Test to check TMatrix <-> array interplay"""
    

    logger = getLogger( 'test_linangt_2')

    M = ROOT.TMatrixD
    
    for k in range ( 1 , 11  ) :
        for n in range ( 1 , 11  ) :

            m = M ( k , n )

            for r in range ( m.GetNrows() ) :
                for c in range ( m.GetNcols () ) :
                    m [ r , c ] = random.uniform ( -10 , 10 )

            ## create the array from matrix 
            a = array.array ( 'd' , m )
            assert len(a) == m.GetNoElements() , \
                   'Invalid array is created!'

            ## recteare it from the array 
            m2 = M ( k , n , a )
            
            assert m == m2 , 'Matrix is nor recreated!'

            logger.info ('Test with TMatrixD(%2d,%2d) is %s' % ( k , n , m == m2 ) ) 
            
    MS = ROOT.TMatrixDSym

    for k in range ( 1 , 11  ) :

        m = MS ( k)

        for r in range ( m.GetNrows() ) :
            for c in range ( m.GetNcols () ) :
                m [ r , c ] = random.uniform ( -10 , 10 )
                
        ## create the array from matrix 
        a = array.array ( 'd' , m )
        assert len ( a ) == m.GetNoElements () , 'Invalid array is created!'
    
        ## recreate it from the array 
        m2 = MS ( k , a )

        assert m == m2 , 'Matrix is nor recreated!'

        logger.info ('Test with TMatrixDSym(%2d) is %s' % ( k , m == m2 ) ) 

# =============================================================================
if '__main__' == __name__ :
    

    test_linalgt_vct  () 
    test_linalgt_mtrx () 
    test_linalgt_old  ()
    test_linalgt_arr  ()


# =============================================================================
##                                                                      The END 
# =============================================================================
