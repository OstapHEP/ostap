#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
""" Test module for ostap/math/bspline.py.
"""
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_bspline'  ) 
else                       : logger = getLogger ( __name__        )
# ============================================================================= 
import random  
import ostap.math.models 
import ostap.math.bspline
from   ostap.core.core  import Ostap

# ============================================================================
##  test solution of equation  B(x) = c
def test_solve ():
    """Test solution of equation  B(x) = c 
    """

    # 1) construct spline with known roots
    
    ## roots in [0,1]
    troots  =          [  random.uniform(0.01,0.99) for i in  range(4) ]
    ## roots in [1,10]
    roots   = troots + [  random.uniform(1.01,9.99) for i in  range(4) ]
    ## complex roots  
    croots  = [  complex ( random.uniform(-1,3),random.gauss(0,1) )  for i in  range(4) ]

    ## Bernstein polynomial with known roots
    bs = Ostap.Math.Bernstein (  0 , 1 , roots , croots )

    ## convert Bernstein polynomial into B-spline 
    bs = Ostap.Math.BSpline   ( bs )

    ## add several internal knots ...
    for i in range(6) :
        x = random.uniform ( 0.001 , 0.999 )
        bs.insert ( x )
        
    ##  find root of B-spline:
        
    rr = bs.solve()
    logger.info ('Roots found : %s' % list( rr ) )
    
    troots.sort() 
    logger.info ('Roots true  : %s' % troots     )


# ============================================================================
##  test spline interpolation 
def test_interpolation ():
    """Test spline interpolation
    """
    from math import sin,pi, sqrt

    fun = lambda x  : sin(2*pi*x) 
    bs  =  ostap.math.bspline.interpolate ( fun , None, [0] + [  random.uniform(0.01,0.99) for i in range(50) ] + [1] , 2 )

    from ostap.stats.counters import SE
    s = SE()
    for i in range(10000) :
        x = random.uniform ( 0 , 1 )
        vf = fun(x)
        vb = bs (x) 
        s += vf-vb
            
    logger.info ('Interpolation quality %s' % s )


# ============================================================================
##  test spline approxmation
def test_approximation  ():
    """Test spline approximation
    """
    from math import sin,pi, sqrt

    fun = lambda x  : sin(2*pi*x) 
    bs  =  ostap.math.bspline.approximate ( fun , [0] + [  random.uniform(0.01,0.99) for i in range(50) ] + [1] , 2 )

    from ostap.stats.counters import SE
    s = SE()
    for i in range(10000) :
        x = random.uniform ( 0 , 1 )
        vf = fun(x)
        vb = bs (x) 
        s += vf-vb
            
    logger.info ('Approximation quality %s' % s )

    
    
    
    
    
# =============================================================================
if '__main__' == __name__ :

    test_solve         ()
    test_interpolation ()
    test_approximation ()
    
# =============================================================================
# The END 
# =============================================================================
