#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_derivative.py
#  Test module for the file ostap/math/derivative.py
# ============================================================================= 
""" Test module for ostap/math/derivative.py

It tests local implementation of numerical derivatives 
"""
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_deriv' )
else                       : logger = getLogger ( __name__        )
# ============================================================================= 


def test_derivative ():

    from math import sin, cos, pi 
    from ostap.math.derivative import derivative, iszero 
    import random
    from ostap.stats.counters import SE
    from ostap.utils.timing   import timing
    
    cnt  = {} 
    cntE = {} 
    for I in range(1,9) : 
        cnt  [I] = SE()

    func = lambda x :     sin(10*x)+x
    deri = lambda x :  10*cos(10*x)+1
    
    func = lambda x :           sin(x)/x 
    deri = lambda x : (cos(x) - sin(x)/x)/x 

    ## func = lambda x :   x**9
    ## deri = lambda x : 9*x**8  

    ## from math import sinh, cosh
    ## func = lambda x :    sinh(2*x)
    ## deri = lambda x :  2*cosh(2*x)

    ## from math import tanh,cosh
    ## func = lambda x :      tanh(3*x)
    ## deri = lambda x :  3./(cosh(3*x)**2)

    with timing() :
        
        for i in range(50000) :
            
            x = random.uniform ( 0 , pi )
            
            d_true = deri ( x )
            
            for I in range(1,9) : 
                delta = derivative ( func , x , I = I , err = True ) - d_true
                ## delta = derivative ( sin , x , I = I ) - d_true
                ## cnt [I] += delta*1.e+10  
                ## cnt [I] += abs(delta/d_true)*1.e+10
                if not iszero ( delta.cov2() ) :   
                    cnt [I] += abs(delta.value())/delta.error()
                
    for I in range(1,9) : 
        print I , cnt  [I]
        
        
# =============================================================================
if '__main__' == __name__ :

    test_derivative ()
    
# =============================================================================
# The END 
# =============================================================================
