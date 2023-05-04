#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_integral.py
#  Test module for the file ostap/math/integral.py
# 
# It tests local implementation of numerical integratuon using Romberg's method
# - see https://en.wikipedia.org/wiki/Romberg's_method
# - see https://en.wikipedia.org/wiki/Richardson_extrapolation
# - see https://en.wikipedia.org/wiki/Numerical_integration#Adaptive_algorithms
# CPU performance is not superb, but it is numerically stable.
#
# Also it tests Ostap::Math::Integrator
# - 1D,2D&3D integrations
#  @see Ostap::Math::Integrator   
# ============================================================================= 
""" Test module for ostap/math/integral.py

It tests local implementation of numerical integratuon using Romberg's method
- see https://en.wikipedia.org/wiki/Romberg's_method
- see https://en.wikipedia.org/wiki/Richardson_extrapolation
- see https://en.wikipedia.org/wiki/Numerical_integration#Adaptive_algorithms
CPU performance is not superb, but it is numerically stable.

Also it tests Ostap::Math::Integrator
 - 1D,2D&3D integrations
 
"""
# =============================================================================
from   ostap.core.meta_info     import root_info 
from   ostap.core.pyrouts       import Ostap, SE  
from   ostap.utils.timing       import timing
from   ostap.utils.progress_bar import progress_bar
from   ostap.math.integral      import ( integral  , romberg     ,
                                         clenshaw_curtis         , 
                                         integral2 , genzmalik2  ,  
                                         integral3 , genzmalik3  ,
                                       complex_circle_integral )
from   ostap.math.make_fun      import make_fun1, make_fun2 , make_fun3
import ostap.math.integrator 
import ostap.logger.table       as     T
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_math_integral' )
else                       : logger = getLogger ( __name__             )
# ============================================================================= 
    
# =============================================================================
def test_integral ():

    logger = getLogger('test_integral')
    logger.info ( 'Simple test for 1D-integrations' )
    
    from math import sin, cos , exp, log, pi, e  

    funcs = [
        ( sin , 0      , pi , 2   ) ,
        ( cos , 0      , pi , 0   ) ,
        ( exp , 0      , 1  , e-1 ) ,
        ( log , 1.e-8  , 1  , -1  )
        ]

    scale = 1.e+12    
    rows  = [  ( 'Function' ,  'I' , 'R' , 'CC'    , \
                 'd(I) [%.0e]'  % ( 1.0 / scale )  , \
                 'd(R) [%.0e]'  % ( 1.0 / scale )  , \
                 'd(CC) [%.0e]' % ( 1.0 / scale )  , \
                 'r(I)' , 'r(R)' , 'r(CC)' ) ]
    
    for entry in funcs :

        args = entry[:3]
        vi   = integral        ( *entry[:3] , err = True )
        vr   = romberg         ( *entry[:3] , err = True , maxdepth = 200 , epsrel = 1.e-8 , epsabs = 1.e-8 )
        vc   = clenshaw_curtis ( *entry[:3] , err = True , maxdepth = 200 , epsrel = 1.e-8 , epsabs = 1.e-8 )
        
        value   = entry[ 3]
        
        func    = 'int(%s,%+.5f,%+.5f)' % ( entry[0].__name__ , entry[1] , entry[2] )

        row = func                                       , \
              '%+.5f'  % vi                              , \
              '%+.5f'  % vr                              , \
              '%+.5f'  % vc                              , \
              '%+.5f'  % ( ( vi - value ) * scale )      , \
              '%+.5f'  % ( ( vr - value ) * scale )      , \
              '%+.5f'  % ( ( vc - value ) * scale )      , \
              '%+5.3f' % ( ( vi - value ) / vi.error() ) , \
              '%+3.3f' % ( ( vr - value ) / vr.error() ) , \
              '%+3.3f' % ( ( vc - value ) / vr.error() )
        
        rows.append ( row )
        
    title = '1D integrations '
    table = T.table ( rows , title = title ,  prefix = '# ' , alignment = 'lllc' )
    logger.info ( '%s\n%s' % ( title , table ) ) 
    


# =============================================================================
def test_inf_integrals ():

    logger = getLogger('test_inf_integrals')

    logger.info ( 'Test 1D-integrals for (semi)infinite intervals' )
    ## if root_info < ( 6 , 18 ) :
    ##     logger.warning ( "Test is disabled for %s" % str ( root_info ) )
    ##    return 

    from math               import pi, exp 
    from ostap.math.math_ve import sech 

    scale = 1.e+12    
    rows  = [  ( 'Function' , 'I' , 'd(I) [%.0e]' % ( 1.0 / scale )  ) ]

    I    = Ostap.Math.Integrator()
    
    aprec = 1.e-12
    rprec = 1.e-12

    def sech ( x ) :
        return  0 if abs (  x ) > 200 else 2.0 / ( exp ( x ) + exp ( -x ) )
    
    # (1) sech 
    
    fun = lambda x : sech ( x ) 
    ff  = make_fun1 ( fun )
    
    vi = I.integrate_infinity      ( ff ,       0 , aprec , rprec )
    vl = I.integrate_to_infinity   ( ff , 0.0 , 0 , aprec , rprec )
    vr = I.integrate_from_infinity ( ff , 0.0 , 0 , aprec , rprec )
    
    row = 'sech  (-inf,+inf)' , '%+.5f' % vi , '%+.5f' % ( ( vi - pi   ) * scale )  
    rows.append ( row )
    
    row = 'sech  (   0,+inf)' , '%+.5f' % vl , '%+.5f' % ( ( vl - pi/2 ) * scale )  
    rows.append ( row )
    
    row = 'sech  (-inf,   0)' , '%+.5f' % vr , '%+.5f' % ( ( vr - pi/2 ) * scale ) 
    rows.append ( row )

    # (2) Cauchy
    
    fun = lambda x : 1.0 / ( x * x + 1.0 ) 
    ff  = make_fun1 ( fun )
    
    vi = I.integrate_infinity      ( ff ,       0 , aprec , rprec )
    vl = I.integrate_to_infinity   ( ff , 0.0 , 0 , aprec , rprec )
    vr = I.integrate_from_infinity ( ff , 0.0 , 0 , aprec , rprec )
    
    row = 'Cauchy(-inf,+inf)' , '%+.5f' % vi , '%+.5f' % ( ( vi - pi   ) * scale )  
    rows.append ( row )
    
    row = 'Cauchy(   0,+inf)' , '%+.5f' % vl , '%+.5f' % ( ( vl - pi/2 ) * scale )  
    rows.append ( row )
    
    row = 'Cauchy(-inf,   0)' , '%+.5f' % vr , '%+.5f' % ( ( vr - pi/2 ) * scale ) 
    rows.append ( row )
    

    title = '1D integrations'
    table = T.table ( rows , title = title ,  prefix = '# ' , alignment = 'lllc' )
    logger.info ( '%s\n%s' % ( title , table ) ) 


# =============================================================================
def test_cauchy_integrals ():

    logger = getLogger('test_cauchy_integrals')

    logger.info ( 'Test Cauchy principal value intervals' )
    ## if root_info < ( 6 , 18 ) :
    ##    logger.warning ( "Test is disabled for %s" % str ( root_info ) )
    ##    return 

    from math               import pi, exp 
    from ostap.math.math_ve import sech 

    scale = 1.e+12    
    rows  = [  ( 'Function' , 'I' ) ]

    I    = Ostap.Math.Integrator()
    
    aprec = 1.e-12
    rprec = 1.e-12

    def sech ( x ) :
        return  0 if abs (  x ) > 200 else 2.0 / ( exp ( x ) + exp ( -x ) )
    
    # (1) sech 
    
    fun = lambda x : sech ( x ) 
    ff  = make_fun1 ( fun )
    
    v1 = I.cauchy_pv                ( ff , 0.0 , -10 , 10 , 0 , 0 , aprec , rprec )
    v2 = I.cauchy_pv                ( ff , 1.0 , -10 , 10 , 0 , 0 , aprec , rprec )
    v3 = I.cauchy_pv                ( ff , 2.0 , -10 , 10 , 0 , 0 , aprec , rprec )
    
    row = 'sech/pv   (0,-10,+10)' , '%+.5f' % v1
    rows.append ( row )
    
    row = 'sech/pv   (1,-10,+10)' , '%+.5f' % v2 
    rows.append ( row )
    
    row = 'sech/pv   (2,-10,+10)' , '%+.5f' % v3 
    rows.append ( row )

    # (2) Cauchy
    
    fun = lambda x : 1.0 / ( x * x + 1.0 ) 
    ff  = make_fun1 ( fun )
    
    v1 = I.cauchy_pv                ( ff , 0.0 , -10 , 10 , 0 , 0 , aprec , rprec )
    v2 = I.cauchy_pv                ( ff , 1.0 , -10 , 10 , 0 , 0 , aprec , rprec )
    v3 = I.cauchy_pv                ( ff , 2.0 , -10 , 10 , 0 , 0 , aprec , rprec )
    
    row = 'cauchy/pv (0,-10,+10)' , '%+.5f' % v1
    rows.append ( row )
    
    row = 'cauchy/pv (1,-10,+10)' , '%+.5f' % v2 
    rows.append ( row )
    
    row = 'cauchy/pv (2,-10,+10)' , '%+.5f' % v3 
    rows.append ( row )

    title = '1D integrations'
    table = T.table ( rows , title = title ,  prefix = '# ' , alignment = 'lllc' )
    logger.info ( '%s\n%s' % ( title , table ) ) 
    

        
# =============================================================================
def test_integrators ():

    logger = getLogger('test_integrators')

    from math import sin, pi 

    ## function to be integrated 
    ff    = lambda x : sin ( x )  
    low   = 0.0
    high  = 1.0 * pi 
    exact = 2.0
    
    ## create the function object 
    f1 = make_fun1 ( ff ) 


    N    = 100000
    
    scale = 1.e+12

    def my_romberg         ( *args ) :
        return romberg ( *args , epsrel = 1.e-11 , epsabs = 1.e-11 )
    
    def my_clenshaw_curtis ( *args ) :
        return clenshaw_curtis ( *args , epsrel = 1.e-11 , epsabs = 1.e-11 ) 
    
    results = []

    ## if root_info < ( 6 , 18 ) :
    ##     for_test = ( ( 'Native/1'  , f1 ,   integral          ) ,
    ##                  ( 'Native/2'  , ff ,   integral          ) ,
    ##                  ( 'Romberg/1' , f1 ,  my_romberg         ) ,
    ##                  ( 'Romberg/2' , ff ,  my_romberg         ) )
    ## else :
    I = Ostap.Math.Integrator() 
    for_test = ( ( 'QAG'              , f1 , I.integrate         ) ,
                 ( 'CQUAD'            , f1 , I.integrate_cquad   ) ,
                 ( 'Romberg'          , f1 , I.integrate_romberg ) ,
                 ( 'Native/1'         , f1 ,   integral          ) ,
                 ( 'Native/2'         , ff ,   integral          ) ,
                 ( 'Romberg/1'        , f1 ,  my_romberg         ) ,
                 ( 'Romberg/2'        , ff ,  my_romberg         ) , 
                 ( 'ClenshawCurtis/1' , f1 ,  my_clenshaw_curtis ) ,
                 ( 'ClenshawCurtis/2' , ff ,  my_clenshaw_curtis ) )
                 
            
    for name , fun , func in for_test : 
        cnt = SE()
        with timing ( '%9s integrator' % name , logger = logger ) as t :  
            for i in progress_bar ( range ( N ) ) : cnt += abs( func ( fun , low , high ) - exact ) * scale 
        results.append ( ( name , cnt , t.delta ) )

    rows = [ ( 'Integrator' , 'CPU [s]' , 'delta [%.0e]' % ( 1.0/scale ) , 'max [%.0e]' % ( 1.0/scale ) ) ]

    for name , cnt, td in results :        
        row = name                         , \
              '%.2f'  % td                 , \
              '%+.4f' % cnt.mean().value() , \
              '%+.4f' % ( cnt.max() )
        rows.append ( row )

    title = 'Compare different integrators'
    table = T.table ( rows , title = title ,  prefix = '# ' , alignment = 'lllc' )
    logger.info ( '%s\n%s' % ( title , table ) ) 
    

# =============================================================================
def test_integral_2D ():

    logger = getLogger('test_integral_2D')

    logger.info ( 'Test for 2D-integration' )

    from math import sin, cos , exp, log, pi, e  

    funcs = [
        ( 'x*x+y*y'                     , lambda x,y: x*x+y*y                     , -1 , 1 , -1  , 1 , 8/3. )  ,
        ( 'x*x+y*sin(y)'                , lambda x,y: x*x+y*sin(y)                , -1 , 1 , -2  , 1 , 6.085519557719447 ) , 
        ( 'x*exp(x)+y*sin(u)'           , lambda x,y: x*exp(x)+y*sin(y)           , -2 , 1 , -2  , 2 , 12.07356999835915445374 ) ,
        ( 'sin(x)*cos(x)+exp(y)*sin(y)' , lambda x,y: sin(x)*cos(x)+exp(y)*sin(y) ,  0 , 3 ,  0  , 3 , 35.60837524996321690196 ) 
        ]

    scale = 1.e+12    
    rows  = [  ( 'Function' , 'GM' , 'I2' , \
                 'd(GM) [%.0e]' % ( 1.0 / scale )  , \
                 'd(I2) [%.0e]' % ( 1.0 / scale )  , 'r(GM)' , 'r(I2)' ) ]
    
    for entry in funcs :

        v1    = genzmalik2 ( *entry[1:6] , err = True )
        v2    = integral2  ( *entry[1:6] , err = True )
        vv    = entry[-1]

        row = '%s(%+.1f,%+.1f,%+.1f,%+.1f)' % ( entry[0] , entry[2] , entry[3] , entry[4] , entry [5] ) , \
              '%+.5f'  % v1                             , '%+.5f'  % v2                             , \
              '%+.5f'  % ( ( v1 - vv   ) * scale )      , '%+.5f'  % ( ( v2 - vv ) * scale        ) , \
              '%+5.3f' % ( ( v1 - vv   ) / v1.error() ) , '%+3.3f' % ( ( v2 - vv ) / v2.error() )

        rows.append ( row )


    title = '2D integrations '
    table = T.table ( rows , title = title ,  prefix = '# ' , alignment = 'lllc' )
    logger.info ( '%s\n%s' % ( title , table ) ) 


# =============================================================================
def test_integrators_2D ():
    
    logger = getLogger('test_integrators_2D')

    logger.info ( 'Test for 2D-integrators' )

    from math import sin, cos , exp, log, pi, e  

    funcs = [
        ( 'x*x+y*y'                     , lambda x,y: x*x+y*y                     , -1 , 1 , -1  , 1 , 8/3. )  ,
        ( 'x*x+y*sin(y)'                , lambda x,y: x*x+y*sin(y)                , -1 , 1 , -2  , 1 , 6.085519557719447 ) , 
        ( 'x*exp(x)+y*sin(u)'           , lambda x,y: x*exp(x)+y*sin(y)           , -2 , 1 , -2  , 2 , 12.07356999835915445374 ) ,
        ( 'sin(x)*cos(x)+exp(y)*sin(y)' , lambda x,y: sin(x)*cos(x)+exp(y)*sin(y) ,  0 , 3 ,  0  , 3 , 35.60837524996321690196 ) 
        ]

    scale = 1.e+12    
    rows  = [  ( 'Function' , 'GM' , 'I' , \
                 'd(GM) [%.0e]' % ( 1.0 / scale )  , \
                 'd(I) [%.0e]'  % ( 1.0 / scale )  , 'r(GM)' , 'r(I)' ) ]

    N = 1000

    results = []
    
    I = Ostap.Math.Integrator()

    for entry in funcs :

        vv    = entry[-1]

        name  = '%s(%+.1f,%+.1f.%+.1f.%+.1f)' % ( entry[0] , entry[2] , entry[3] , entry[4] , entry [5] )
        
        with timing ( 'GenzMalik2' ) as td :
            cnt  = SE() 
            for i in progress_bar ( range  ( N ) ) :
                cnt +=  abs ( genzmalik2 ( *entry[1:6] , err = True ) - vv ) * scale                
        results.append (  ( name , 'GenzMalik2' , cnt , td.delta ) ) 

        with timing ( 'Integral' ) as td :
            cnt  = SE() 
            for i in progress_bar ( range  ( N ) ) :
                cnt +=  abs ( integral2 ( *entry[1:6] , err = True ) - vv ) * scale                
        results.append (  ( name , 'Integral2' , cnt , td.delta ) ) 

        ## if (6,18) <= root_info : 
        with timing ( 'Cubature2' ) as td :
            cnt  = SE()
            fn   = make_fun2 ( entry[1] )
            for i in progress_bar ( range  ( N ) ) :
                args = entry[2:6] + ( 0,0,0) 
                cnt +=  abs ( I.integrate2 (  fn , *args ) - vv ) * scale                
        results.append (  ( name , 'Cubature2' , cnt , td.delta ) ) 
                    
    rows = [ ( 'Function' , 'Integrator' , 'CPU [s]' , 'delta [%.0e]' % ( 1.0/scale ) , 'max [%.0e]' % ( 1.0/scale ) ) ]

    for name , method, cnt , td in results :
        row = name , method , '%.2f' % td , \
              '%+.4f' % cnt.mean().value() , \
              '%+.3f' % ( cnt.max() )
        rows.append ( row )

    title = 'Compare different 2D integrators'
    table = T.table ( rows , title = title ,  prefix = '# ' , alignment = 'lllc' )
    logger.info ( '%s\n%s' % ( title , table ) ) 
    
# =============================================================================
def test_integral_3D ():

    logger = getLogger('test_integral_3D')

    logger.info ( 'Test for 3D-integration' )

    from math import sin, cos , exp, log, pi, e  

    from ostap.math.integral import  genzmalik3, integral3  
    funcs = [
        ( 'x*x+y*|y|'                      , lambda x,y,z: x*x+y*abs(y)                       , -1. , 2. , -1.  , 2. , -1. , 2. , 48 ) ,  
        ( 'x*x+|y|+z*z*cos(z)'             , lambda x,y,z: x*x+y*abs(y)+z*z*cos(z)            , -1. , 2. , -1.  , 2. , -1. , 2. , 51.53827020952059001502 ) ,       
        ( 'x*exp(x)+sin(y)*|y|+z*z*exp(z)' , lambda x,y,z: x*exp(x)+sin(y)*abs(y)+z*z*exp(z)  , -1. , 2. , -1.  , 2. , -1. , 2. , 202.53557154832049036486) ,       
        ]
    
    scale = 1.e+12    
    rows  = [  ( 'Function' , 'GM' , 'I3' , \
                 'd(GM) [%.0e]' % ( 1.0 / scale )  , \
                 'd(I3) [%.0e]' % ( 1.0 / scale )  , 'r(GM)' , 'r(I3)' ) ]
    
    for entry in funcs :
        
        v1    = genzmalik3 ( *entry[1:8] , err = True )
        v2    = integral3  ( *entry[1:8] , err = True )
        vv    = entry[-1]

        row = '%s(%+.1f,%+.1f.%+.1f.%+.1f,%+.1f,%+.1f)' % ( entry[0] , entry[2] , entry[3] , entry[4] , entry [5] , entry[6] , entry[7] ) , \
              '%+.5f'  % v1                             , '%+.5f'  % v2                             , \
              '%+.5f'  % ( ( v1 - vv   ) * scale )      , '%+.5f'  % ( ( v2 - vv ) * scale        ) , \
              '%+5.3f' % ( ( v1 - vv   ) / v1.error() ) , '%+3.3f' % ( ( v2 - vv ) / v2.error() )

        rows.append ( row )

    title = '3D integrations '
    table = T.table ( rows , title = title ,  prefix = '# ' , alignment = 'lllc' )
    logger.info ( '%s\n%s' % ( title , table ) ) 



# =============================================================================
def test_integrators_3D ():
    
    logger = getLogger('test_integrators_3D')

    logger.info ( 'Test for 3D-integrators' )

    from math import sin, cos , exp, log, pi, e  

    from ostap.math.integral import  genzmalik3, integral3  
    funcs = [
        ( 'x*x+y*|y|'                      , lambda x,y,z: x*x+y*abs(y)                       , -1. , 2. , -1.  , 2. , -1. , 2. , 48 ) ,  
        ( 'x*x+|y|+z*z*cos(z)'             , lambda x,y,z: x*x+y*abs(y)+z*z*cos(z)            , -1. , 2. , -1.  , 2. , -1. , 2. , 51.53827020952059001502 ) ,       
        ( 'x*exp(x)+sin(y)*|y|+z*z*exp(z)' , lambda x,y,z: x*exp(x)+sin(y)*abs(y)+z*z*exp(z)  , -1. , 2. , -1.  , 2. , -1. , 2. , 202.53557154832049036486) ,       
        ]
    

    scale = 1.e+10 
    rows  = [  ( 'Function' , 'GM' , 'I' , \
                 'd(GM) [%.0e]' % ( 1.0 / scale )  , \
                 'd(I) [%.0e]'  % ( 1.0 / scale )  , 'r(GM)' , 'r(I)' ) ]

    N = 20

    results = []
    
    I = Ostap.Math.Integrator()
    
    for entry in funcs :

        vv    = entry[-1]

        name  = '%s(%+.1f,%+.1f,%+.1f.%+.1f,%+.1f,%.1f)' % ( entry[0] , entry[2] , entry[3] , entry[4] , entry [5] , entry[6] , entry[7])
        
        with timing ( 'GenzMalik3' ) as td :
            cnt  = SE() 
            for i in progress_bar ( range  ( N ) ) :
                cnt +=  abs ( genzmalik3 ( *entry[1:8] , err = True ) - vv ) * scale                
        results.append (  ( name , 'GenzMalik3' , cnt , td.delta ) ) 

        with timing ( 'Integral' ) as td :
            cnt  = SE() 
            for i in progress_bar ( range  ( N ) ) :
                cnt +=  abs ( integral3 ( *entry[1:8] , err = True ) - vv ) * scale                
        results.append (  ( name , 'Integral3' , cnt , td.delta ) ) 

        ## if ( 6,18) <= root_info : 
        with timing ( 'Cubature3' ) as td :
            cnt  = SE()
            fn   = make_fun3 ( entry[1] )
            for i in progress_bar ( range  ( N ) ) :
                args = entry[2:8] + ( 0, 0.0 , 0.0 ) 
                cnt +=  abs ( I.integrate3 (  fn , *args ) - vv ) * scale                
        results.append (  ( name , 'Cubature3' , cnt , td.delta ) ) 
                    
    rows = [ ( 'Function' , 'Integrator' , 'CPU [s]' , 'delta [%.0e]' % ( 1.0/scale ) , 'max [%.0e]' % ( 1.0/scale ) ) ]

    for name , method, cnt , td in results :
        row = name , method , '%.2f' % td , \
              '%+.4f' % cnt.mean().value() , \
              '%+.3f' % ( cnt.max() )
        rows.append ( row )

    title = 'Compare different 3D integrators'
    table = T.table ( rows , title = title ,  prefix = '# ' , alignment = 'lllc' )
    logger.info ( '%s\n%s' % ( title , table ) ) 
 

# =============================================================================
def test_integral_contour ():

    logger = getLogger('test_integral_contour')
    logger.info ( 'Test for complex contour integration' )

    import cmath
    
    from   math import pi

    ## function with simple poles

    poles = ( -1j , +1j )     
    cc    = 1.0 / ( 2 * pi ) 
    func  = lambda x : cc * sum ( 1.0 / ( x - p ) for p in poles )

    rows = [  ( 'center' , 'radius' , 'result' , 'delta [10^-12]' ) ] 
              
    ## radius 
    for radius in ( 0.5 , 5.5 ) :
        ## center
        for real in ( -1 , 0 , 1 ) :
            for imag in ( -1 , 0 , 1 ) :
                center = real+imag*1j

                result = complex_circle_integral ( func , center = center , radius = radius )

                exact = 0
                for p in poles :
                    if abs ( center - p ) < radius : exact += 1.0j

                row = '%-+.1f%+.1fj' % ( center.real , center.imag )  , \
                      '%.2f' % radius , \
                      '%-+.6f%+.6fj' % ( result.real , result.imag )  , \
                      '%.3e'   %  ( abs ( exact - result ) * 1e+12 ) 
                       
                rows.append ( row )


    title = 'Complex Contour integrals'
    table = T.table ( rows , title = title ,  prefix = '# ' )
    logger.info ( '%s\n%s' % ( title , table ) ) 


# ==============================================================================
if '__main__' == __name__ :

    test_integral         ()
    test_integral_2D      ()
    test_integral_3D      ()
    test_integral_contour ()
    
    test_integrators      ()
    test_integrators_2D   ()
    test_integrators_3D   ()

    test_inf_integrals    ()
    test_cauchy_integrals ()


# =============================================================================
##                                                                      The END 
# =============================================================================
