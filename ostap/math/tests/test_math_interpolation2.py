#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# =============================================================================
## @file ostap/math/tests/test_math_interpolation3.py
#  Test module for the file ostap/math/interpolation.py
# ============================================================================= 
""" Test module for ostap/math/interpolation.py
"""
# ============================================================================= 
import ROOT, random, math  
import ostap.math.linalg

import ostap.math.models
from   ostap.math.interpolation import ( Berrut1st      ,
                                         Berrut2nd      ,
                                         Barycentric    ,
                                         FloaterHormann )

from   ostap.core.core          import Ostap
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_interpolation2' ) 
else                       : logger = getLogger ( __name__                         )
# =============================================================================

# =============================================================================
def test_matrices0 () :

    logger = getLogger ( 'test_matrices0' )
    logger.info ( 'Test matrices' )
    
    MT = Ostap.Math.Matrix(4,4) 

    fun = lambda x : x 

    x0  = 0
    m0  = MT()
    for i in range( m0.kRows ) :
        for j in range ( i ,  m0.kCols ) :
            m0 [i,j] = (i+j)**2 

    x1 = 1
    m1  = MT()
    for i in range( m0.kRows ) :
        for j in range ( min ( i + 1 , m1.kCols ) ) :
            m1 [i,j] = (i+j)**2

    data = {} 
    data [ x0 ] = m0
    data [ x1 ] = m1
    
    interpolants = [
        ( 'Berrut1st'   , Berrut1st   ( data ) ) ,
        ( 'Berrut2nd'   , Berrut2nd   ( data ) ) ,
        ( 'Barycentric' , Barycentric ( data ) ) ,
        ]
    for d in range ( 2 ) :
         item = 'FloaterHormann/%s' % d , FloaterHormann ( data , degree = d ) 
         interpolants.append ( item )
         
    xs = [ random.uniform ( 0 , 1  ) for i in range ( 3 ) ]
    xs.sort()
     
    for x in xs : 
        for n, i in interpolants :            
            logger.info ( 'Interpolant %15s x=%.4g \n%s' % ( n , x , i ( x ) ) ) 
    
# =============================================================================
def test_matrices1 () :

    logger = getLogger ( 'test_matrices1' )
    logger.info ( 'Test matrices1' )
    
    MT = Ostap.Math.Matrix(3,3) 

    fun = lambda x : x 

    N , low , high = 5 , 0 , 10 

    xs   = [ random.uniform ( low , high ) for i in range ( N ) ]
    xs.append ( 0       )
    xs.append ( math.pi )

    data = {}
    for k, x in enumerate ( xs ) : 
        m = MT ()
        for i in range( m.kRows ) :
            for j in range( m.kCols ) :
                m [i,j] = fun ( x ) 
        data [ x ] = m

    interpolants = [
        ( 'Berrut1st'   , Berrut1st   ( data ) ) ,
        ( 'Berrut2nd'   , Berrut2nd   ( data ) ) ,
        ( 'Barycentric' , Barycentric ( data ) ) ,
        ]
    for d in range ( 5 ) :
         item = 'FloaterHormann/%s' % d , FloaterHormann ( data , degree = d ) 
         interpolants.append ( item )

         
    xs = [ random.uniform ( low , high  ) for i in range ( 3 ) ]
    xs.sort()
    
    for x in xs : 
        for n, i in interpolants :            
            logger.info ( 'Interpolant %15s x=%.4g \n%s' % ( n , x , i ( x ) ) ) 
        
# =============================================================================
if '__main__' == __name__ :
    
    test_matrices0    ()
    test_matrices1    ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
