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
from   ostap.math.base          import doubles 
from   ostap.core.core          import Ostap,  SE
from   ostap.math.models        import f1_draw 
from   ostap.utils.utils        import wait
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
import ostap.logger.table       as     T 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_interpolation' ) 
else                       : logger = getLogger ( __name__                        )
# =============================================================================

## calcualet "distance" between two functions 
def distance ( fun1 , fun2 , low , high ) :
    """calculate ``distance'' between two functions"""

    df = lambda x : abs(fun1(x)-fun2(x))

    from ostap.math.integral import integral 
    di = integral ( df , low , high ) 

    return di / ( high - low ) 
    
# =============================================================================
## interpolate the function 
def run_func_interpolation ( fun , N , low , high , scale = 1.e-5 , logger = logger , name = 'Interpolation') :
    """Interpolate the function"""

    Abscissas =  Ostap.Math.Interpolation.Abscissas
    
    abscissas = ( ( 'Uniform'   , Abscissas ( N , low , high, 0 ) ) ,
                  ( 'Chebyshev' , Abscissas ( N , low , high, 1 ) ) ,
                  ( 'Lobatto'   , Abscissas ( N , low , high, 2 ) ) ,
                  ## ( 'Random'    , Abscissas ( doubles (  random.uniform ( low , high ) for i in range ( N ) ) ) )
                  ) 
    
    tables       = [ ( a[0] , points ( fun , a[1] ) )  for a in abscissas ]
    
    interpolants = []

    for i , t in enumerate ( tables ) :
        
        item = ( 'Bernstein'   , t[0] ) , interpolate_bernstein  ( t[1] , None , low , high )
        interpolants.append ( item )
        
        item = ( 'Neville'     , t[0] ) ,  Ostap.Math.Neville     ( t[1] ) 
        interpolants.append ( item )
        
        item = ( 'Lagrange'    , t[0] ) , Ostap.Math.Lagrange     ( t[1] ) 
        interpolants.append ( item )
        
        item = ( 'Newton'      , t[0] )  , Ostap.Math.Newton      ( t[1] ) 
        interpolants.append ( item )
        
        item = ( 'Berrut1st'   , t[0] )  , Ostap.Math.Berrut1st   ( t[1] ) 
        interpolants.append ( item )

        item = ( 'Berrut2nd'   , t[0] )  , Ostap.Math.Berrut2nd   ( t[1] ) 
        interpolants.append ( item )
        
        item = ( 'Barycentric' , t[0] )  , Ostap.Math.Barycentric ( t[1] ) 
        interpolants.append ( item )

        for d in range ( 0 , 9 ) :
            item = ( 'FloaterHormann%d' % d , t[0] ) , Ostap.Math.FloaterHormann ( t[1] , d )
            interpolants.append ( item )
        
        for d in range ( 1 , 5 ) :
            item = ( 'BSpline%d' % d , t[0] ) , interpolate_bspline  ( t[1] , None , d )
            interpolants.append ( item )
            
 
    with wait ( 3 ) , use_canvas ( name ) :
        ff = lambda x : fun  ( x )
        f1_draw ( ff , xmin = low , xmax = high , linecolor = 2 , linewidth = 2 )

        for i , item in enumerate ( interpolants ) :

            p  , f  = item
            n1 , n2 = p 
            
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
            
            logger.info ( 'Color %10s for %s:%s' % ( color , n1 , n2  ) ) 

    xx = []
    NP = 50000
    for i in range ( NP ) : xx.append ( random.uniform ( low , high ) ) 
    xx.sort()

    
    from collections import defaultdict
    counters = defaultdict(SE) 
    
    cpu = {}
    ## loop over all interpolants 
    for n , fi in interpolants :

        n1 , n2 = n
        
        cnt = counters [  ( n1 , n2 , fi ) ]

        with timing ( '' , logger = None )  as t :
            for x in xx :            
                v  = fun   ( x )
                vi = fi    ( x )
                cnt += abs ( vi - v ) / scale
                
        cpu [ (n1,n2) ] = t.delta

    rows = [ ( 'Interpolant' , 'Grid' , 'mean+/-rms'  , 'max' , 'distance') ]
    for item in counters :

        n1 , n2 , ff = item
        
        c  = counters[item]
        
        d  = distance ( ff , fun , low , high ) 
        row = n1 , n2  , '%9.2f +/- %-09.1f' % ( c.mean().value() , c.rms() ) , '%-9.1f' % c.max() , '%.3g' % d 
        rows.append ( row )

    title = 'Interpolation precision (%d random points)[x%s]' % ( NP , scale  ) 
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lllll' )
    
    logger.info ( '%s:\n%s' % ( title , table ) ) 

    lst  = [] 
    for k in cpu :
        item = cpu[k] , k[0], k[1] 
        lst.append ( item )
    lst.sort()
    
    rows = [ ( 'Interpolant' , 'Grid' , 'CPU [s]' ) ]
    for t,k0,k1 in lst :
        row = k0, k1 , '%.4g' % cpu[ (k0,k1) ]
        rows.append ( row )
        
    title = 'CPU: %d points' % NP 
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'll' )
    
    logger.info ( '%s:\n%s' % ( title , table ) ) 


# ==========================================================================================================
## interpolate the grid 
def run_grid_interpolation ( tfunc , dct , N , low , high , scale = 1.e-8 , logger = logger , name = 'interpolation' ) :
    """Interpolate the grid"""
    
    Abscissas =  Ostap.Math.Interpolation.Abscissas
        
    data = points ( dct )

    ## list of interpolants 
    interpolants = []
    
    ## bernstein interpolant 
    interpolants.append ( ( 'Bernstein'   , interpolate_bernstein  ( data , None , low , high ) ) ) 
    
    ## neville interpolant
    interpolants.append ( ( 'Neville'     , Ostap.Math.Neville     ( data ) ) )
    
    ## largange interpolant 
    interpolants.append ( ( 'Lagrange'    , Ostap.Math.Lagrange    ( data ) ) ) 
    
    ## (true) Barycentric interpolant 
    interpolants.append ( ( 'Barycentric' , Ostap.Math.Barycentric ( data ) ) ) 
    
    ## Newton interpolant 
    interpolants.append ( ( 'Newton'      , Ostap.Math.Newton      ( data ) ) ) 
    
    ## 1st Berrut interpolant 
    interpolants.append ( ( 'Berrut 1st'  , Ostap.Math.Berrut1st    ( data ) ) ) 
    
    ## 2nd Berrut interpolant 
    interpolants.append ( ( 'Berrut 2nd'  , Ostap.Math.Berrut2nd    ( data ) ) ) 
    
    for d in range ( 10 ) :
        interpolants.append ( ( 'FloaterHormann/%d' % d  , Ostap.Math.FloaterHormann ( data , d ) ) ) 
        
    ## bspline interpolant
    ## bs = Ostap.Math.BSpline   ( low  , high ,  len ( data ) - 1 - degree , degree  )
    for d in range ( 1 , 5 ) : 
        interpolants.append ( ( 'BSpline/%s' % d , interpolate_bspline  ( data , None , d ) ) ) 
    
    with wait ( 1 ) , use_canvas ( name ) :
        
        ff = lambda x : tfunc  ( x )
        f1_draw ( ff , xmin = low , xmax = high , linecolor = 1 , linewidth = 3 )
        
        for i , c in enumerate ( interpolants ) :
            
            n , f  = c

            color = i + 2 
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
    NP = 50000
    for i in range ( NP ) : xx.append ( random.uniform ( low , high ) ) 
    xx.sort()

    from collections import defaultdict
    counters = defaultdict(SE) 
    
    cpu = {}
    ## loop over all interpolants 
    for n , fi in interpolants :

        cnt = counters [  ( n , fi ) ]

        with timing ( '' , logger = None )  as t :
            for x in xx :            
                v  = tfunc ( x )
                vi = fi    ( x )
                cnt += abs ( vi - v ) / scale
                
        cpu [ n ] = t.delta

    rows = [ ( 'Configuration' , 'mean+/-rms'  , 'max' , 'distance') ]
    for item in counters :

        n , ff = item
        c     = counters[item] 

        d  = distance ( ff , tfunc , low , high ) 
        row = n , '%9.2f +/- %-09.1f' % ( c.mean().value() , c.rms() ) , '%-9.1f' % c.max() , '%.3g' % d 
        rows.append ( row )

    title = 'Interpolation precision (%d random points)[x%s]' % ( NP , scale  ) 
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'llll' )
    
    logger.info ( '%s:\n%s' % ( title , table ) ) 

    lst  = [] 
    for k in cpu :
        item = cpu[k] , k
        lst.append ( item )
    lst.sort()
    
    rows = [ ( 'Interpolant' , 'CPU [s]' ) ]
    for t,k in lst :
        row = k , '%.4g' % cpu[k]
        rows.append ( row )
        
    title = 'CPU: %d points' % NP 
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'll' )
    
    logger.info ( '%s:\n%s' % ( title , table ) ) 

    
# =============================================================================
## interpolate cos function 
def test_cos () :

    logger = getLogger ( 'test_cos' ) 
    fun , N , low , high = math.cos , 10 , 0 , 2 * math.pi
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( 'cos' , N , low , high )  ) 
    
    return run_func_interpolation ( fun , N , low , high , logger = logger , name = 'cos(x)') 


# =============================================================================
## interpolate |sin| function 
def test_abssin () :
    
    logger = getLogger ( 'test_abssin' ) 

    fun = lambda x : abs( math.sin ( x ) )

    N , low , high = 20 , -0.1 , 0.1
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( '|sin|' , N , low , high )  ) 
    
    return run_func_interpolation ( fun , N , low , high , scale = 1.e-5 , logger = logger , name = '|sin(x)|' ) 

# =============================================================================
## interpolate |sin(2x)| function 
def test_abs2sin () :
    
    logger = getLogger ( 'test_abs2sin' ) 

    fun = lambda x : abs( math.sin ( 2*x ) )

    N , low , high = 24 , 0 , math.pi  
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( '|sin(2x)|' , N , low , high )  ) 
    
    random.seed ( 50948524584 )
    return run_func_interpolation ( fun , N , low , high , scale = 1.e-3 , logger = logger , name = '|sin(2x)|') 



# =============================================================================
## interpolate the table of values 
def test_random_grid_sin () :

    
    logger = getLogger ( 'test_random_grid_sin' ) 

    tfun = math.sin 
    
    N , low , high = 14 , 0 , 3 * math.pi  

    dct = {} 
    mid = 0.5  * ( low + high ) 
    while N > len ( dct ) :
        xi = random.uniform ( low , mid  )
        dct [ xi ] = tfun ( xi ) 
        xi = random.uniform ( mid , high )
        dct [ xi ] = tfun ( xi ) 

    dct [ low  ] = tfun ( low  )
    dct [ high ] = tfun ( high )
    
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( 'sin' , N , low , high )  )
    
    return run_grid_interpolation ( tfun , dct , N , low , high , scale = 1.e-4 , logger = logger , name = 'sin(x)') 


# =============================================================================
## interpolate the table of values 
def test_random_grid_abssin () :

    
    logger = getLogger ( 'test_random_grid_abssin' ) 

    tfun = lambda x : abs ( math.sin ( x ) ) 
    
    N , low , high = 14 , 0 , math.pi  

    dct = {} 
    mid = 0.5  * ( low + high ) 
    while N > len ( dct ) :
        xi = random.uniform ( low , mid  )
        dct [ xi ] = tfun ( xi ) 
        xi = random.uniform ( mid , high )
        dct [ xi ] = tfun ( xi ) 

    dct [ low  ] = tfun ( low  )
    dct [ high ] = tfun ( high )
    
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( '|sin|' , N , low , high )  )
    
    return run_grid_interpolation ( tfun , dct , N , low , high , scale = 1.e-4 , logger = logger , name = '|sin(x)|') 

# =============================================================================
## interpolate the table of values 
def test_random_grid_sin2 () :

    
    logger = getLogger ( 'test_random_grid_sin2' ) 

    tfun = lambda x : math.sin ( x ) **2  
    
    N , low , high = 14 , 0 , 2 * math.pi  

    dct = {} 
    mid = 0.5  * ( low + high ) 
    while N > len ( dct ) :
        xi = random.uniform ( low , mid  )
        dct [ xi ] = tfun ( xi ) 
        xi = random.uniform ( mid , high )
        dct [ xi ] = tfun ( xi ) 

    dct [ low  ] = tfun ( low  )
    dct [ high ] = tfun ( high )
    
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( '|sin|' , N , low , high )  )
    
    return run_grid_interpolation ( tfun , dct , N , low , high , scale = 1.e-4 , logger = logger , name = 'sin^2(x)') 


# =============================================================================
## interpolate the table of values 
def test_random_grid_gauss () :

    logger = getLogger ( 'test_random_grid_gauss' ) 

    tfun = lambda x : Ostap.Math.gauss_pdf ( x ) 
    
    N , low , high = 14 , -4 , 4

    dct = {} 
    mid = 0.5  * ( low + high ) 
    while N > len ( dct ) :
        xi = random.uniform ( low , mid  )
        dct [ xi ] = tfun ( xi ) 
        xi = random.uniform ( mid , high )
        dct [ xi ] = tfun ( xi ) 

    dct [ low  ] = tfun ( low  )
    dct [ high ] = tfun ( high )
    
    logger.info ( 'Interpolate %12s, %3d points, (%s,%s) interval' %  ( 'gauss' , N , low , high )  )
    
    return run_grid_interpolation ( tfun , dct , N , low , high , scale = 1.e-3 , logger = logger , name = 'gauss') 
    
# =============================================================================
if '__main__' == __name__ :
    
    ## test_cos                () 
    ## test_abssin             ()
    ## test_abs2sin            ()
    test_random_grid_sin    ()
    test_random_grid_abssin ()
    test_random_grid_sin2   ()
    test_random_grid_gauss  ()
    
    
# =============================================================================
##                                                                      The END 
# =============================================================================
