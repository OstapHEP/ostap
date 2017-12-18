#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
""" Test module for ostap/math/integral.py.
It tests local implementation of numerical integrtauon using Romberg's method
- see https://en.wikipedia.org/wiki/Romberg's_method
- see https://en.wikipedia.org/wiki/Richardson_extrapolation
- see https://en.wikipedia.org/wiki/Numerical_integration#Adaptive_algorithms
CPU performance is not superb, but it is numerically stable.
"""
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_integral' )
else                       : logger = getLogger ( __name__        )
# ============================================================================= 


def test_integral ():

    from math import sin, cos , exp, log, pi, e  

    from ostap.math.integral import integral, romberg 
    funcs = [
        ( sin , 0      , pi , 2   ) ,
        ( cos , 0      , pi , 0   ) ,
        ( exp , 0      , 1  , e-1 ) ,
        ( log , 1.e-50 , 1  , -1  )
        ]
    
    for entry in funcs :

        args = entry[:3]
        vi = integral ( *entry[:3] , err = True )
        vr = romberg  ( *entry[:3] , err = True , maxdepth = 100 )
        
        value   = entry[ 3]
        func    = 'int(%s,%g,%g)' % ( entry[0].__name__ , entry[1] , entry[2] )
        logger.info ( '%20s: I        %-20s %-20s' % ( func , vi       , vr   ) ) 
        logger.info ( '%20s: Delta    %-20s %-20s' % ( func , vi-value , vr - value  ) )
        logger.info ( '%20s: Delta/I  %-20s %-20s' % ( func , (vi-value)/vi.error() ,
                                                       (vr - value)/vr.error() ) )
        
def test_integral_2D ():

    from math import sin, cos , exp, log, pi, e  

    from ostap.math.integral import  genzmalik2, integral2  
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
        print vv, '%.20f'    % v2.value()  
        logger.info ( '%20s: I        %-20s %-20s  %-20s ' % ( entry[0] , v1 , v2  , vv ) ) 
        logger.info ( '%20s: Delta    %-20s %-20s'         % ( entry[0] , v1-vv , v2 - vv ) )
        logger.info ( '%20s: Delta/I  %-20s %-20s'         % ( entry[0] , (v1-vv)/vv         , (v2 - vv)/vv         ) ) 
        logger.info ( '%20s: Delta/E  %-20s %-20s'         % ( entry[0] , (v1-vv)/v1.error() , (v2 - vv)/v2.error() ) ) 
                

        
def test_integral_3D ():

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
        print vv, '%.20f'    % v2.value()  
        logger.info ( '%20s: I        %-20s %-20s  %-20s ' % ( entry[0] , v1 , v2  , vv ) ) 
        logger.info ( '%20s: Delta    %-20s %-20s'         % ( entry[0] , v1-vv , v2 - vv ) )
        logger.info ( '%20s: Delta/I  %-20s %-20s'         % ( entry[0] , (v1-vv)/vv         , (v2 - vv)/vv         ) ) 
        logger.info ( '%20s: Delta/E  %-20s %-20s'         % ( entry[0] , (v1-vv)/v1.error() , (v2 - vv)/v2.error() ) ) 
                


        
# =============================================================================
if '__main__' == __name__ :

    test_integral    ()
    test_integral_2D ()
    test_integral_3D ()
    
# =============================================================================
# The END 
# =============================================================================
