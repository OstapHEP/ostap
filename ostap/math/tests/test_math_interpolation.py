#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# =============================================================================
## @file ostap/math/tests/test_math_interpolation.py
#  Test module for the file ostap/math/interpolation.py
# ============================================================================= 
""" Test module for ostap/math/interpolation.py
"""
# ============================================================================= 
import random,math 
import ostap.math.models
from   ostap.math.interpolation import ( interpolate , points  ,
                                         interpolate_bernstein ,
                                         interpolate_bspline   ) 
from   ostap.core.core          import Ostap,  SE
from   ostap.math.models        import f1_draw 
from   ostap.utils.utils        import wait
from   ostap.plotting.canvas    import use_canvas
import ostap.logger.table       as     T 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_interpolation' ) 
else                       : logger = getLogger ( __name__                        )
# =============================================================================


# =============================================================================
## interpolate the function 
def run_func_interpolation (  fun , N , low , high , scale = 1.e-5 , logger = logger , name = 'Interpolation') :

    Abscissas =  Ostap.Math.Interpolation.Abscissas

    ##  uniform abscissas 
    i0 = interpolate ( fun , Abscissas ( N , low , high , 0 ) )

    ## Chebyshev abscissas 
    i1 = interpolate ( fun , Abscissas ( N , low , high , 1 ) )

    ## lobatto abscissas 
    i2 = interpolate ( fun , Abscissas ( N , low , high , 2 ) )

    ## random abscissas 
    i3 = interpolate ( fun , [ random.uniform ( low , high ) for i in range ( N ) ] )

    ## bernstein interpolant 
    i4 = interpolate_bernstein ( fun ,  Abscissas ( N , low , high , 2 ) , low , high )

    ## bspline interpolation
    degree = 4 
    bs = Ostap.Math.BSpline   ( low , high ,  N - 1 - degree , degree  )
    i5 = interpolate_bspline  ( fun , Abscissas ( N , low , high , 1 ) , bs ) 
    
    config =  [ ( i0 , 'Uniform   abscissas'   ) ,
                ( i1 , 'Chebyshev abscissas'   ) ,
                ( i2 , 'Lobatto   abscissas'   ) ,
                ( i3 , 'Random    abscissas'   ) ,
                ( i4 , 'Bernstein interpolant' ) ,
                ( i5 , 'bSpline   interpolant' ) ]
    
    with wait ( 1 ) , use_canvas ( name ) :
        ff = lambda x : fun  ( x )
        f1_draw ( ff , xmin = low , xmax = high , linecolor = 2 , linewidth = 2 )
        for i , c in enumerate ( config ) :
            f , n = c
            color = i + 3 
            f.draw ( 'same' , linecolor = color )
            if   1 == color : color = 'Black'
            elif 2 == color : color = 'Red'
            elif 3 == color : color = 'Green'
            elif 4 == color : color = 'Blue'
            elif 5 == color : color = 'Yellow'
            elif 6 == color : color = 'Magenta'
            elif 7 == color : color = 'Cyan'
            elif 8 == color : color = 'DarkGreen'
            
            logger.info ( 'Color %10s for %s' % ( color , n ) ) 

    xx = []
    NP = 100000
    for i in range ( NP ) : xx.append ( random.uniform ( low , high ) ) 
    xx.sort()

    c0 = SE()
    c1 = SE()
    c2 = SE()
    c3 = SE()
    c4 = SE()
    c5 = SE()
    
    for x in xx :
        
        f  = fun      ( x ) 
        f0 = i0       ( x ) 
        f1 = i1       ( x ) 
        f2 = i2       ( x ) 
        f3 = i3       ( x ) 
        f4 = i4       ( x ) 
        f5 = i5       ( x ) 

        d0 = f0 - f 
        d1 = f1 - f
        d2 = f2 - f 
        d3 = f3 - f 
        d4 = f4 - f 
        d5 = f5 - f 
        
        c0 += abs ( d0 ) / scale  
        c1 += abs ( d1 ) / scale 
        c2 += abs ( d2 ) / scale
        c3 += abs ( d3 ) / scale
        c4 += abs ( d4 ) / scale
        c5 += abs ( d5 ) / scale

    config =  [ ( c0 , 'Uniform   abscissas'   ) ,
                ( c1 , 'Chebyshev abscissas'   ) ,
                ( c2 , 'Lobatto   abscissas'   ) ,
                ( c3 , 'Random    abscissas'   ) ,
                ( c4 , 'Bernstein interpolant' ) ,
                ( c5 , 'bSpline   interpolant' ) ]
    
    rows = [ ( 'Configuration' , 'mean+/-rms'  , 'max' ) ]
    for c , n in config :
        row = n , '%9.2f +/- %-09.1f' % ( c.mean().value() , c.rms() ) , '%-9.1f' % c.max()
        rows.append ( row )

    title = 'Interpolation precision (%d random points)[x%s]' % ( NP , scale  ) 
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lll' )
    
    logger.info ( '%s:\n%s' % ( title , table ) ) 
        

# ==========================================================================================================
## interpolate the grid 
def run_grid_interpolation ( tfunc , dct , N , low , high , scale = 1.e-8 , logger = logger , name = 'interpolation' ) :
    
    Abscissas =  Ostap.Math.Interpolation.Abscissas
        
    data = points ( dct )

    ##  uniform abscissas 
    i0 = interpolate            ( data )
    
    ## bernstein interpolant 
    i1 = interpolate_bernstein  ( data , None , low , high )

    ## neville interpolant
    i2 = Ostap.Math.Neville     ( data )

    ## largange interpolant 
    i3 = Ostap.Math.Lagrange    ( data )

    ## Barycentric interpolant 
    i4 = Ostap.Math.Barycentric ( data )

    ## newton interpolant 
    i5 = Ostap.Math.Newton      ( data )

    ## bspline interpolant
    degree = 4
    
    ## bs = Ostap.Math.BSpline   ( low  , high ,  len ( data ) - 1 - degree , degree  )
    i6 = interpolate_bspline  ( data , None , degree ) 

    config = [ ( i0 , 'Uniform     abscissas'   ) ,
               ( i1 , 'Bernstein   interpolant' ) ,
               ( i2 , 'Neville     interpolant' ) ,
               ( i3 , 'Lagrange    interpolant' ) ,
               ( i4 , 'Barycentric interpolant' ) ,
               ( i5 , 'Newton      interpolant' ) ,
               ( i6 , 'bSpline     interpolant' ) ]
               
    with wait ( 1 ) , use_canvas ( name ) :
        ff = lambda x : tfunc  ( x )
        f1_draw ( ff , xmin = low , xmax = high , linecolor = 2 , linewidth = 2 )
        for i , c in enumerate ( config ) :
            f , n = c

            print (i,c) 

                
            color = i + 3 
            f.draw ( 'same' , linecolor = color )
            if   1 == color : color = 'Black'
            elif 2 == color : color = 'Red'
            elif 3 == color : color = 'Green'
            elif 4 == color : color = 'Blue'
            elif 5 == color : color = 'Yellow'
            elif 6 == color : color = 'Magenta'
            elif 7 == color : color = 'Cyan'
            elif 8 == color : color = 'DarkGreen'
            
            logger.info ( 'Color %10s for %s' % ( color , n ) ) 
    
        
    xx = []
    NP = 100000
    for i in range ( NP ) : xx.append ( random.uniform ( low , high ) ) 
    xx.sort()

    c0 = SE ()
    c1 = SE ()
    c2 = SE ()
    c3 = SE ()
    c4 = SE ()
    c5 = SE ()
    c6 = SE ()
    
    for x in xx :

        f  = tfunc    ( x )
        
        f0 = i0       ( x ) 
        f1 = i1       ( x ) 
        f2 = i2       ( x )
        f3 = i3       ( x )
        f4 = i4       ( x )
        f5 = i5       ( x )
        f6 = i6       ( x )
        
        d0 = f0 - f 
        d1 = f1 - f
        d2 = f2 - f
        d3 = f3 - f
        d4 = f4 - f
        d5 = f5 - f
        d6 = f6 - f
        
        c0 += abs ( d0 ) / scale  
        c1 += abs ( d1 ) / scale 
        c2 += abs ( d2 ) / scale 
        c3 += abs ( d3 ) / scale 
        c4 += abs ( d4 ) / scale 
        c5 += abs ( d5 ) / scale 
        c6 += abs ( d6 ) / scale 

    config = [ ( c0 , 'Uniform     abscissas'   ) ,
               ( c1 , 'Bernstein   interpolant' ) ,
               ( c2 , 'Neville     interpolant' ) ,
               ( c3 , 'Lagrange    interpolant' ) ,
               ( c4 , 'Barycentric interpolant' ) ,
               ( c5 , 'Newton      interpolant' ) ,
               ( c6 , 'bSpline     interpolant' ) ]
    
    rows = [ ( 'Configuration' , 'mean+/-rms'  , 'max' ) ]
    for c , n in config :
        row = n , '%9.2f +/- %-09.1f' % ( c.mean().value() , c.rms() ) , '%-9.1f' % c.max()
        rows.append ( row )

    title = 'Interpolation precision (%d random points)[x%s]' % ( NP , scale  ) 
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lll' )
    
    logger.info ( '%s:\n%s' % ( title , table ) ) 
        

# =============================================================================
## interpolate cos function 
def test_cos () :

    logger = getLogger ( 'test_cos' ) 
    fun , N , low , high = math.cos , 8 , 0 , 2 * math.pi
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( 'cos' , N , low , high )  ) 
    
    return run_func_interpolation ( fun , N , low , high , logger = logger , name = 'cos(x)') 


# =============================================================================
## interpolate |sin| function 
def test_abssin () :
    
    logger = getLogger ( 'test_abssin' ) 

    fun = lambda x : abs( math.sin ( x ) )

    N , low , high = 12 , -0.1 , 0.1
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( '|sin|' , N , low , high )  ) 
    
    return run_func_interpolation ( fun , N , low , high , scale = 1.e-5 , logger = logger , name = '|sin(x)|' ) 

# =============================================================================
## interpolate |sin(2x)| function 
def test_abs2sin () :
    
    logger = getLogger ( 'test_abs2sin' ) 

    fun = lambda x : abs( math.sin ( 2*x ) )

    N , low , high = 18 , 0 , math.pi  
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( '|sin|' , N , low , high )  ) 
    
    random.seed ( 50948524584 )
    return run_func_interpolation ( fun , N , low , high , scale = 1.e-3 , logger = logger , name = '|sin(2x)|') 



# =============================================================================
## interpolate the table of values 
def test_random_sin () :

    
    logger = getLogger ( 'test_random_sin' ) 

    tfun = math.sin 
    
    N , low , high = 10 , 0 , 4 * math.pi  

    dct = {} 
    while N > len ( dct ) :
        xi = random.uniform ( low , high )
        dct [ xi ] = tfun ( xi ) 

    dct [ low  ] = tfun ( low  )
    dct [ high ] = tfun ( high )
    
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( 'sin' , N , low , high )  )
    
    return run_grid_interpolation ( tfun , dct , N , low , high , scale = 1.e-4 , logger = logger , name = 'sin(x)') 


# =============================================================================
## interpolate the table of values 
def test_random_gauss () :

    
    logger = getLogger ( 'test_random_sin' ) 

    tfun = lambda x : Ostap.Math.gauss_pdf ( x ) 
    
    N , low , high = 12 , -5 , 5 

    dct = {} 
    while N > len ( dct ) :
        xi = random.uniform ( low , high )
        dct [ xi ] = tfun ( xi ) 

    dct [ low  ] = tfun ( low  )
    dct [ high ] = tfun ( high )
    
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( 'sin' , N , low , high )  )
    
    return run_grid_interpolation ( tfun , dct , N , low , high , scale = 1.e-3 , logger = logger , name = 'gauss') 
    
# =============================================================================
if '__main__' == __name__ :
    
    ## test_cos       () 
    ## test_abssin    ()
    ## test_abs2sin   ()
    test_random_sin   ()
    test_random_gauss ()
    
    
# =============================================================================
##                                                                      The END 
# =============================================================================
