#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_derivative.py
#  Test module for the file ostap/math/derivative.py
# ============================================================================= 
""" Test module for ostap/math/derivative.py

It tests local implementation of numerical derivatives 
"""
# ============================================================================= 
from __future__ import print_function
# ============================================================================= 
import random
from   math                  import sin, cos, pi 
from   ostap.math.derivative import derivative, iszero 
from   ostap.stats.counters  import SE
from   ostap.utils.timing    import timing
import ostap.logger.table    as     T 
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_derivative' )
else                       : logger = getLogger ( __name__                    )
# ============================================================================= 


def test_derivative_1 ():

    logger = getLogger ( 'test_derivative_1' )
    
    cnt  = {} 
    cntE = {} 
    for I in range(1,9) : 
        cnt  [I] = SE()

    
    func = lambda x :           sin(x)/x 
    deri = lambda x : (cos(x) - sin(x)/x)/x 

    with timing() :
        
        for i in range(10000) :
            
            x = random.uniform ( 0 , pi )
            
            d_true = deri ( x )
            
            for I in range ( 1 , 9 ) :
                delta    = derivative ( func , x , I = I ) - d_true
                cnt [I] += delta

    rows = [ ('Order' , '#' , 'min' , 'max' ) ]
    for i in cnt :
        c = cnt[I] 
        row = '%d' % i , '%d' % c.nEntries() , '%+.5g' % c.min() , '%+.5g' % c.max() 
        rows.append ( row ) 
        
    table = T.table ( rows ,
                      title  = 'Test numerical derivatives' ,
                      prefix = '# ' , alignment = 'crll'     )
    
    logger.info ( 'Test numerical derivatives\n%s' % table )


def test_derivative_2 ():

    logger = getLogger ( 'test_derivative_2' )
    
    cnt  = {} 
    cntE = {} 
    for I in range(1,9) : 
        cnt  [I] = SE()

    func = lambda x :     sin(10*x)+x
    deri = lambda x :  10*cos(10*x)+1
    

    with timing() :
        
        for i in range(10000) :
            
            x = random.uniform ( 0 , pi )
            
            d_true = deri ( x )
            
            for I in range ( 1 , 9 ) : 
                delta    = derivative ( func , x , I = I ) - d_true
                cnt [I] += delta

    rows = [ ('Order' , '#' , 'min' , 'max' ) ]
    for i in cnt :
        c = cnt[I] 
        row = '%d' % i , '%d' % c.nEntries() , '%+.5g' % c.min() , '%+.5g' % c.max() 
        rows.append ( row ) 

    table = T.table ( rows ,
                      title  = 'Test numerical derivatives' ,
                      prefix = '# ' , alignment = 'crll'     )
    
    logger.info ( 'Test numerical derivatives\n%s' % table )

        
# =============================================================================
if '__main__' == __name__ :

    test_derivative_1 ()
    test_derivative_2 ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
