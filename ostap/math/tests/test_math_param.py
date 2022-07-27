#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file ostap/math/tests/test_math_parampy
# Test module for ostap/math/param.py
# - It tests parameterisation of functions
# ============================================================================= 
""" Test module for ostap/math/param.py
- It tests parameterisations of functions
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import 
# ============================================================================= 
from   builtins               import range
import ostap.histos.param
import ostap.histos.histos
import ostap.fitting.funcs
from   ostap.plotting.canvas  import use_canvas
from   ostap.utils.utils      import wait 
from   ostap.utils.timing     import timing
from   ostap.math.models      import f1_draw 
import ROOT, random, math,  time
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'ostap.test_math_param' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( 'Test for function parameterisation')
# =============================================================================
use_scipy = False 
try :
    import scipy
    use_scipy = True 
except ImportError :
    use_scipy = False 

# =============================================================================
## make a single test 
def make_test ( func , xmin , xmax , ptype , orders , logger ) :
    

    results = {}
    for order in orders :
        r = ptype ( func = func  ,
                    xmin = xmin  ,
                    xmax = xmax  ,
                    N    = order ) 
        results [ order ] = r
        

    f1_draw ( func , xmin = xmin , xmax = xmax , linecolor=2 , linewidth = 3 )
    logger.info ( 'Thick red line - original function' ) 
    for o in results : 
        r = results [o]
        r.draw ( 'same' , xmin = xmin , xmax = xmax , linecolor = o ) 
        logger.info ( 'Lone color %2d : approximation with order %s' % ( o , o ) )  
    

# =============================================================================
def test_legendre_sum () :
    """Test Legendre sum
    """

    logger = getLogger('test_legendre_sum')
    logger.info ( 'Represent the function as Legendre sum')

    
    func = lambda x : 1.0 if abs(x)<1.e-6  else math.sin(x)/x 

    xmin = 0
    xmax = 2 * math.pi 

    from ostap.math.param import legendre_sum

    with wait ( 3 ) , use_canvas ( 'test_legendre_sum' ) , timing ( logger = logger ) :
        make_test ( func , xmin , xmax , legendre_sum , range ( 1 , 11 ) , logger )
        

# =============================================================================
def test_chebyshev_sum () :
    """Test Chebyshev sum
    """

    logger = getLogger('test_chebyshev_sum')
    logger.info ( 'Represent the function as Chebyshev sum')

    
    func = lambda x : 1.0 if abs(x)<1.e-6  else math.sin(x)/x 

    xmin = 0
    xmax = 2 * math.pi 

    from ostap.math.param import chebyshev_sum

    
    with wait ( 3 ) , use_canvas ( 'test_chebyshev_sum' ) , timing ( logger = logger ) :
        make_test ( func , xmin , xmax , chebyshev_sum , range ( 1 , 11 ) , logger )

        
# =============================================================================
def test_fourier_sum () :
    """Test Fourier sum
    """

    logger = getLogger('test_fourier_sum')
    logger.info ( 'Represent the function as Fourier sum')


    func = lambda x : 1.0 if abs(x)<1.e-6  else math.sin(x)/x 
    
    xmin = 0 
    xmax = 2 * math.pi 

    try :
        from ostap.math.param import fourier_sum
    except ImportError :
        logger.error ("Can't import fourier_sum!")
        return 

    with wait ( 3 ) , use_canvas ( 'test_fourier_sum' ) , timing ( logger = logger ) :
        make_test ( func , xmin , xmax , fourier_sum , range ( 1 , 11 ) , logger )

# =============================================================================
def test_cosine_sum () :
    """Test Cosine sum
    """

    logger = getLogger('test_cosine_sum')
    logger.info ( 'Represent the function as a sum of cosine functions')
    
    func = lambda x : 1.0 if abs(x)<1.e-6  else math.sin(x)/x 

    xmin = 0 
    xmax = 2 * math.pi 

    try :        
        from ostap.math.param import cosine_sum
    except ImportError :
        logger.error ("Can't import cosine_sum!")
        return 

    
    with wait ( 3 ) , use_canvas ( 'test_cosine_sum' ) , timing ( logger = logger ) :
        make_test ( func , xmin , xmax , cosine_sum , range ( 1 , 11 ) , logger )

# =============================================================================
def test_bernstein_sum () :
    """Test Bernstein sum
    """

    logger = getLogger('test_bernstein_sum')
    logger.info ( 'Represent the function as a sum of Bernstein polynomials')

    
    func = lambda x : 1.0 if abs(x)<1.e-6  else math.sin(x)/x 

    xmin = 0 
    xmax = 2 * math.pi 

    from ostap.math.param import bernstein_sum

    
    with wait ( 3 ) , use_canvas ( 'test_bernstein_sum' ) , timing ( logger = logger ) :
        make_test ( func , xmin , xmax , bernstein_sum , range ( 1 , 11 ) , logger )


# =============================================================================
def test_bernsteineven_sum1 () :
    """Test Bernstein Even sum/1
    """

    logger = getLogger('test_bernsteineven_sum1')
    logger.info ( 'Represent the function as a sum of even Bernstein polynomials')
    
    func = lambda x : 1.0 if abs(x)<1.e-6  else math.sin(x)/x 
    xmin = -2 * math.pi 
    xmax =  2 * math.pi 

    from ostap.math.param import bernsteineven_sum

    
    with wait ( 3 ) , use_canvas ( 'test_bernsteineven_sum' ) , timing ( logger = logger ) :
        make_test ( func , xmin , xmax , bernsteineven_sum , range ( 1 , 11 ) , logger )

# =============================================================================
def test_bernsteineven_sum2 () :
    """Test Bernstein Even sum/2
    """

    logger = getLogger('test_bernsteineven_sum2')
    logger.info ( 'Represent the function as a sum of even Bernstein polynomials')
    
    func = lambda x : 1.0 if abs(x-2*math.pi)<1.e-6  else math.sin(x-2*math.pi)/(x-2*math.pi)
    xmin =  0 * math.pi 
    xmax =  4 * math.pi 

    from ostap.math.param import bernsteineven_sum

    
    with wait ( 3 ) , use_canvas ( 'test_bernsteineven_sum2' ) , timing ( logger = logger ) :
        make_test ( func , xmin , xmax , bernsteineven_sum , range ( 1 , 21 ) , logger )

        
# =============================================================================
if '__main__' == __name__ :

    test_legendre_sum       ()
    test_chebyshev_sum      ()
    test_fourier_sum        ()
    test_cosine_sum         ()
    test_bernstein_sum      ()
    test_bernsteineven_sum1 ()
    test_bernsteineven_sum2 ()

    
# =============================================================================
##                                                                      The END
# =============================================================================
