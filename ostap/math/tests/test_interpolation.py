#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# =============================================================================
## @file ostap/math/tests/test_interpolation.py
#  Test module for the file ostap/math/interpolation.py
# ============================================================================= 
""" Test module for ostap/math/interpolation.py
"""
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_interpolation' ) 
else                       : logger = getLogger ( __name__             )
# =============================================================================
import random,math 
import ostap.math.models 
from   ostap.math.interpolation import interpolate, interpolate_bernstein  
from   ostap.core.core          import Ostap,  SE 
# =============================================================================

# =============================================================================
## interpolate the function 
def run_func_interpolation (  fun , N , low , high , scale = 1.e-8 ) :
    
    Abscissas =  Ostap.Math.Interpolation.Abscissas

    ##  uniform abscissas 
    i0 = interpolate ( fun , Abscissas ( N , low , high , 0 ) )

    ## Chebyshev abscissas 
    i1 = interpolate ( fun , Abscissas ( N , low , high , 1 ) )

    ## lobatto abscissas 
    i2 = interpolate ( fun , Abscissas ( N , low , high , 2 ) )

    ## random abscissas 
    i3 = interpolate ( fun , [ random.uniform ( low , high ) for i in   range (N ) ] )

    ## bernstein interpolant 
    i4 = interpolate_bernstein ( fun ,  Abscissas ( N , low , high , 2 ) , low , high )
    
    xx = []
    for i in range ( 100000 ) : xx.append ( random.uniform ( low , high ) ) 
    xx.sort()

    c0 = SE()
    c1 = SE()
    c2 = SE()
    c3 = SE()
    c4 = SE()
    
    for x in xx :
        
        f  = fun      ( x ) 
        f0 = i0       ( x ) 
        f1 = i1       ( x ) 
        f2 = i2       ( x ) 
        f3 = i3       ( x ) 
        f4 = i4       ( x ) 

        d0 = f0 - f 
        d1 = f1 - f
        d2 = f2 - f 
        d3 = f3 - f 
        d4 = f4 - f 
        
        c0 += abs ( d0 ) / scale  
        c1 += abs ( d1 ) / scale 
        c2 += abs ( d2 ) / scale
        c3 += abs ( d3 ) / scale
        c4 += abs ( d4 ) / scale
        
    logger.info ( 'Uniform   precision: mean/max[%s] = %9.2f +- %-09.1f/%-9.1f' % ( scale , c0.mean () . value () , c0.rms () , c0.max () ) )
    logger.info ( 'Chebyshev precision: mean/max[%s] = %9.2f +- %-09.1f/%-9.1f' % ( scale , c1.mean () . value () , c1.rms () , c1.max () ) )
    logger.info ( 'Lobatto   precision: mean/max[%s] = %9.2f +- %-09.1f/%-9.1f' % ( scale , c2.mean () . value () , c2.rms () , c2.max () ) )
    logger.info ( 'Random    precision: mean/max[%s] = %9.2f +- %-09.1f/%-9.1f' % ( scale , c3.mean () . value () , c3.rms () , c3.max () ) )
    logger.info ( 'Bernstein precision: mean/max[%s] = %9.2f +- %-09.1f/%-9.1f' % ( scale , c4.mean () . value () , c4.rms () , c4.max () ) )


# ==========================================================================================================
## interpolate the function 
def run_grid_interpolation ( tfunc , dct , N , low , high , scale = 1.e-8 ) :
    
    Abscissas =  Ostap.Math.Interpolation.Abscissas

    ##  uniform abscissas 
    i0 = interpolate           ( dct )
    
    ## bernstein interpolant 
    i1 = interpolate_bernstein ( dct , None , low , high )

    xx = []
    for i in range ( 100000 ) : xx.append ( random.uniform ( low , high ) ) 
    xx.sort()

    c0 = SE()
    c1 = SE()
    
    for x in xx :

        f  = tfunc    ( x )
        
        f0 = i0       ( x ) 
        f1 = i1       ( x ) 

        d0 = f0 - f 
        d1 = f1 - f
        
        c0 += abs ( d0 ) / scale  
        c1 += abs ( d1 ) / scale 
        
    logger.info ( 'Grid      precision: mean/max[%s] = %9.2f +- %-09.1f/%-9.1f' % ( scale , c0.mean () . value () , c0.rms () , c0.max () ) )
    logger.info ( 'Bernstein precision: mean/max[%s] = %9.2f +- %-09.1f/%-9.1f' % ( scale , c1.mean () . value () , c1.rms () , c1.max () ) )


# =============================================================================
## interpolate cos function 
def test_cos () :
    
    fun , N , low , high = math.cos , 12 , 0 , 2 * math.pi
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( 'cos' , N , low , high )  ) 
    
    random.seed ( 9876543210 )
    return run_func_interpolation ( fun , N , low , high ) 

# =============================================================================
## interpolate cos function 
def test_abssin () :
    
    fun = lambda x : abs( math.sin ( x ) )

    N , low , high = 12 , -0.1 , 0.1
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( '|sin|' , N , low , high )  ) 
    
    random.seed ( 50948524584 )
    return run_func_interpolation ( fun , N , low , high , scale = 1.e-5 ) 

# =============================================================================
## interpolate the table of values 
def test_dict () :

    random.seed ( 50948524584 )
    
    tfun = math.sin 
    
    N , low , high = 12 , -0.5 , 3.0

    dct = {}
    while 12 > len ( dct ) :
        x = random.uniform ( low , high )
        dct[x] = tfun ( x) 
        
    dct[low]  = tfun( low  )
    dct[high] = tfun( high )
    
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( 'sin' , N , low , high )  ) 
    return run_grid_interpolation ( tfun , dct , N , low , high , scale = 1.e-10 ) 


    
# =============================================================================
if '__main__' == __name__ :
    
    test_cos    () 
    test_abssin ()
    test_dict   () 
    
    
# =============================================================================
# The END 
# =============================================================================
