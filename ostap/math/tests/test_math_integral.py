#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_integral.py
#  Test module for the file ostap/math/integral.py
# ============================================================================= 
""" Test module for ostap/math/integral.py

It tests local implementation of numerical integrtauon using Romberg's method
- see https://en.wikipedia.org/wiki/Romberg's_method
- see https://en.wikipedia.org/wiki/Richardson_extrapolation
- see https://en.wikipedia.org/wiki/Numerical_integration#Adaptive_algorithms
CPU performance is not superb, but it is numerically stable.
"""
# =============================================================================
from ostap.core.pyrouts  import Ostap, SE  
from ostap.utils.timing  import timing 
from ostap.math.integral import ( integral  , romberg     , 
                                  integral2 , genzmalik2  ,  
                                  integral3 , genzmalik3  ,
                                  complex_circle_integral ) 
import ostap.logger.table  as     T
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
    logger.info ( 'Simple test for 1D-inttegrtaions' )
    
    from math import sin, cos , exp, log, pi, e  

    funcs = [
        ( sin , 0      , pi , 2   ) ,
        ( cos , 0      , pi , 0   ) ,
        ( exp , 0      , 1  , e-1 ) ,
        ( log , 1.e-50 , 1  , -1  )
        ]
    
    for entry in funcs :

        args = entry[:3]
        vi   = integral ( *entry[:3] , err = True )
        vr   = romberg  ( *entry[:3] , err = True , maxdepth = 100 )
        
        value   = entry[ 3]
        
        func    = 'int(%s,%g,%g)' % ( entry[0].__name__ , entry[1] , entry[2] )
        logger.info ( '%20s: I        %-20s %-20s' % ( func , vi       , vr   ) ) 
        logger.info ( '%20s: Delta    %-20s %-20s' % ( func , vi-value , vr - value  ) )
        logger.info ( '%20s: Delta/I  %-20s %-20s' % ( func , (vi-value)/vi.error() ,
                                                       (vr - value)/vr.error() ) )
        
# =============================================================================
def test_integrator ():

    logger = getLogger('test_integrator')

    logger.info ( 'Simple test for Ostap.Math.Integrator' )

    from math import sin, pi 

    ## function to be integrated 
    ff    = lambda x : sin ( x )  
    low   = 0.0
    high  = 1.0 * pi 
    exact = 2.0
    
    ## create the function object 
    f1 = Ostap.Math.Apply ( ff ) 

    I = Ostap.Math.Integrator() 

    cnt1 = SE ()
    cnt2 = SE ()
    cnt3 = SE ()
    cnt4 = SE ()
    cnt5 = SE ()
    N    = 100000
    
    scale = 1.e+12

    def my_romberg ( *args ) :
        return romberg ( *args , epsrel = 1.e-11 , epsabs = 1.e-11 ) 
    
    results = []
    for name , fun , func in (  ( 'GAC'       , f1 , I.integrate         ) ,
                                ( 'CQUAD'     , f1 , I.integrate_cquad   ) ,
                                ( 'Romberg'   , f1 , I.integrate_romberg ) ,
                                ( 'Native/1'  , f1 ,   integral          ) ,
                                ( 'Native/2'  , ff ,   integral          ) ,
                                ( 'Romberg/1' , f1 ,  my_romberg         ) ,
                                ( 'Romberg/2' , ff ,  my_romberg         ) ,
                                ) :
        cnt = SE()
        with timing ( '%9s integrator' % name , logger = logger ) as t :  
            for i in range ( N ) : cnt += abs( func ( fun , low , high ) - exact ) * scale 
        results.append ( ( name , cnt , t.delta ) )

    rows = [ ( 'Integrator' , 'CPU [s]' , 'delta [%.0e]' % ( 1.0/scale ) , 'max [%.0e]' % ( 1.0/scale ) ) ]

    for name , cnt, td in results :
        
        row = name , '%.2f' % td , \
              '%+.4f' % cnt.mean().value() , \
              '%+.3f' % ( cnt.max() )
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
        ( lambda x,y: x*x+y*y                     , -1 , 1 , -1  , 1 , 8/3. )  ,
        ( lambda x,y: x*x+y*sin(y)                , -1 , 1 , -2  , 1 , 6.085519557719447 ) , 
        ( lambda x,y: x*exp(x)+y*sin(y)           , -2 , 1 , -2  , 2 , 12.07356999835915445374 ) ,
        ( lambda x,y: sin(x)*cos(x)+exp(y)*sin(y) ,  0 , 3 ,  0  , 3 , 35.60837524996321690196 ) 
        ]
    
    for entry in funcs :
        
        v1    = genzmalik2 ( *entry[:5] , err = True )
        v2    = integral2  ( *entry[:5] , err = True )
        vv    = entry[-1]

        logger.info ( '%20s: I        %-20s %-20s  %-20s ' % ( entry[0] , v1 , v2  , vv ) ) 
        logger.info ( '%20s: Delta    %-20s %-20s'         % ( entry[0] , v1-vv , v2 - vv ) )
        logger.info ( '%20s: Delta/I  %-20s %-20s'         % ( entry[0] , (v1-vv)/vv         , (v2 - vv)/vv         ) ) 
        logger.info ( '%20s: Delta/E  %-20s %-20s'         % ( entry[0] , (v1-vv)/v1.error() , (v2 - vv)/v2.error() ) ) 
                
        
# =============================================================================
def test_integral_3D ():

    logger = getLogger('test_integral_3D')

    logger.info ( 'Test for 3D-integration' )

    from math import sin, cos , exp, log, pi, e  

    from ostap.math.integral import  genzmalik3, integral3  
    funcs = [
        ( lambda x,y,z: x*x+y*abs(y)                       , -1 , 2 , -1  , 2 , -1 , 2 , 48 ) ,  
        ( lambda x,y,z: x*x+y*abs(y)+z*z*cos(z)            , -1 , 2 , -1  , 2 , -1 , 2 , 51.53827020952059001502 ) ,       
        ( lambda x,y,z: x*exp(x)+sin(y)*abs(y)+z*z*exp(z)  , -1 , 2 , -1  , 2 , -1 , 2 , 202.53557154832049036486) ,       
        ]
    
    for entry in funcs :
        
        v1    = genzmalik3 ( *entry[:7] , err = True )
        v2    = integral3  ( *entry[:7] , err = True )
        vv    = entry[-1]
        logger.info ( '%20s: I        %-20s %-20s  %-20s ' % ( entry[0] , v1 , v2  , vv ) ) 
        logger.info ( '%20s: Delta    %-20s %-20s'         % ( entry[0] , v1-vv , v2 - vv ) )
        logger.info ( '%20s: Delta/I  %-20s %-20s'         % ( entry[0] , (v1-vv)/vv         , (v2 - vv)/vv         ) ) 
        logger.info ( '%20s: Delta/E  %-20s %-20s'         % ( entry[0] , (v1-vv)/v1.error() , (v2 - vv)/v2.error() ) ) 
                

# =============================================================================
def test_integral_contour ():

    logger = getLogger('test_integral_contour')
    logger.info ( 'Test for coplex contour integration' )

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


# =============================================================================
if '__main__' == __name__ :

    test_integral         ()
    test_integrator       ()
    test_integral_2D      ()
    test_integral_3D      ()
    test_integral_contour ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
