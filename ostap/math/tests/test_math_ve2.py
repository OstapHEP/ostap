#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# =============================================================================
## @file ostap/math/tests/test_math_ve2.py
#  Test module for the file ostap/math/math_ve.py
# ============================================================================= 
""" Test module for ostap/math/math_ve.py
"""
# ============================================================================= 
import math, random 
from   ostap.math.ve       import VE
from   ostap.math.base     import cpp, iszero, isequal
import ostap.logger.table  as     T
from   ostap.math.math_ve  import *
# ============================================================================= 
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_math_ve2' )
else                       : logger = getLogger ( __name__             )
# ============================================================================= 


vars  = [ VE ( 0.001 , 0.0001**2 ) , VE(1,0) , VE(1,0.1**2) , VE(10,0.01**2) ]

funcs = [ exp    , expm1   ,
          log    , log10   , log1p    ,
          sqrt   , cbrt    ,
          sin    , cos     , tan      ,
          sinh   , cosh    , tanh     , sech   ,
          asin   , acos    , atan     ,
          asinh  , acosh   , atanh    ,
          erf    , erfc    , erfi     , erfcx  ,
          sinc   , 
          probit , 
          gamma  , tgamma  , lgamma   , igamma ,
          psi    , digamma , trigamma ,
          ]

## use helper object:
from ostap.math.derivative import EvalVE
funcs += [ EvalVE ( math.sin , math.cos ) ,
           EvalVE ( math.sin )            ]

# =============================================================================
def test_math_ve1 ():

    logger = getLogger("test_math_ve1")

    rows = [ ( 'Function' , 'argument' , 'value' ) ]
    
    for v in vars :
        for f in funcs :
            row = f.__name__ , str ( v ) , str ( f ( v ) )
            rows.append ( row )

    title = 'Functions'
    table = T.table ( rows, title = title , prefix = '# ' , alignment = 'lcc' ) 
    logger.info ( '%s\n%s' % ( title , table ) ) 



# =============================================================================
def test_math_ve2 ():

    logger = getLogger("test_math_ve2")


    funcs = [ ( exp      , -1 , 1 ) ,
              ( expm1    , -1 , 1 ) ,
              ( log      ,  1 , 2 ) ,
              ( log10    ,  1 , 2 ) ,
              ( log1p    ,  1 , 2 ) ,
              ( sqrt     ,  1 , 2 ) ,
              ( cbrt     ,  1 , 2 ) ,              
              ( sin      , -1 , 1 ) ,
              ( cos      , -1 , 1 ) ,
              ( tan      , -1 , 1 ) ,
              ( sinh     , -1 , 1 ) ,
              ( cosh     , -1 , 1 ) ,
              ( tanh     , -1 , 1 ) ,
              ( sech     , -1 , 1 ) ,
              ( asin     , -1 , 1 ) ,
              ( acos     , -1 , 1 ) ,
              ( atan     , -1 , 1 ) ,
              ( asinh    , -1 , 1 ) ,
              ( acosh    ,  1 , 2 ) ,
              ( atanh    , -1 , 1 ) ,
              ( erf      , -1 , 1 ) ,
              ( erfc     , -1 , 1 ) ,
              ( erfi     , -1 , 1 ) ,
              ( erfcx    , -1 , 1 ) ,
              ( sinc     , -1 , 1 ) ,
              ( probit   , -1 , 1 ) ,
              ( gamma    ,  1 , 2 ) ,
              ( tgamma   ,  1 , 2 ) ,
              ( lgamma   ,  1 , 2 ) ,
              ( igamma   ,  1 , 2 ) ,
              ( psi      ,  1 , 2 ) ,
              ( digamma  ,  1 , 2 ) ,
              ( trigamma ,  1 , 2 ) ]

    rows = [ ( 'Function' , 'argument' , 'value' ) ]
    
    for f,a,b in funcs :

        x = random.uniform ( a , b ) 
        row = f.__name__ , '%s' % x , '%s' % f ( x ) 
        ## row = f.__name__ , '%+5.3f'  % x , '+.3g' % f ( x ) )
        rows.append ( row )

    title = 'Functions'
    table = T.table ( rows, title = title , prefix = '# ' , alignment = 'lcc' ) 
    logger.info ( '%s\n%s' % ( title , table ) ) 


# =============================================================================
def test_math_ve3 ():

    logger = getLogger("test_math_ve3")


    funcs = ( bessel_J , bessel_Y ,bessel_I , bessel_K )

    rows = [ ( 'Function' , 'order' ,  'argument' , 'value' ) ]
    
    from   ostap.utils.gsl       import gslCount

    with gslCount() :
        
        for f in funcs :
            
            for n in range ( 0 , 4 ) : 
                x   = random.uniform ( 0 , 1 )
                row = f.__name__ , '%s' % n , '%s' % x , '%s' % f ( n , x )        
                rows.append ( row )
                
                x   = VE ( random.uniform ( 1 , 2 ) , 0.1 **2 ) 
                row = f.__name__ , '%s' % n , '%s' % x , '%s' % f ( n , x )        
                rows.append ( row )
                
            n   = random.uniform ( 1 , 3 )
            row = f.__name__ , '%+.3f' % n , '%s' % x , '%s' % f ( n , x )        
            rows.append ( row )
            
            x   = VE ( random.uniform ( 1 , 2 ) , 0.1 **2 ) 
            row = f.__name__ , '%+.3f' % n , '%s' % x , '%s' % f ( n , x )        
            rows.append ( row )

    title = 'Functions'
    table = T.table ( rows, title = title , prefix = '# ' , alignment = 'lccc' ) 
    logger.info ( '%s\n%s' % ( title , table ) ) 


# =============================================================================
def test_math_ve4 ():

    logger = getLogger("test_math_ve4")

    import ostap.math.derivative as     D
    import ostap.math.integral   as     I
    from   ostap.core.core       import Ostap, SE  
    from   ostap.math.models     import f1_draw
    from   ostap.plotting.canvas import use_canvas
    from   ostap.utils.utils     import wait 
    from   ostap.utils.gsl       import gslCount
    
    functions = [
        ##
        ( bessel_J , 0   , Ostap.Math.der_bessel_Jn  ) ,
        ( bessel_J , 0.0 , Ostap.Math.der_bessel_Jnu ) ,
        ( bessel_J , 2   , Ostap.Math.der_bessel_Jn  ) ,
        ( bessel_J , 2.0 , Ostap.Math.der_bessel_Jnu ) ,
        ( bessel_J , 1.5 , Ostap.Math.der_bessel_Jnu ) ,
        ( bessel_J , 2.5 , Ostap.Math.der_bessel_Jnu ) ,
        ##
        ( bessel_Y , 0   , Ostap.Math.der_bessel_Yn  ) ,
        ( bessel_Y , 0.0 , Ostap.Math.der_bessel_Ynu ) ,
        ( bessel_Y , 2   , Ostap.Math.der_bessel_Yn  ) ,
        ( bessel_Y , 2.0 , Ostap.Math.der_bessel_Ynu ) ,
        ( bessel_Y , 1.5 , Ostap.Math.der_bessel_Ynu ) ,
        ( bessel_Y , 2.5 , Ostap.Math.der_bessel_Ynu ) ,
        ##
        ( bessel_I , 0   , Ostap.Math.der_bessel_In  ) ,
        ( bessel_I , 0.0 , Ostap.Math.der_bessel_Inu ) ,
        ( bessel_I , 2   , Ostap.Math.der_bessel_In  ) ,
        ( bessel_I , 2.0 , Ostap.Math.der_bessel_Inu ) ,
        ( bessel_I , 1.5 , Ostap.Math.der_bessel_Inu ) ,
        ( bessel_I , 2.5 , Ostap.Math.der_bessel_Inu ) ,
        ##
        ##
        ( bessel_K , 0   , Ostap.Math.der_bessel_Kn  ) ,
        ( bessel_K , 0.0 , Ostap.Math.der_bessel_Knu ) ,
        ( bessel_K , 2   , Ostap.Math.der_bessel_Kn  ) ,
        ( bessel_K , 2.0 , Ostap.Math.der_bessel_Knu ) ,
        ( bessel_K , 1.5 , Ostap.Math.der_bessel_Knu ) ,
        ( bessel_K , 2.5 , Ostap.Math.der_bessel_Knu ) ,
        ##
        ]

    keep = []
    
    scale1 = 1.e+15 
    scale2 = 1.e+12 
    rows = [ ( 'Derivative' , 'nu' ,
               'delta [%.0g]' % ( 1/scale1 )  ,
               'amean [%.0g]' % ( 1/scale1 )  , 'amax [%.0g]' % ( 1/scale1 )  , 
               'rmean [%.0g]' % ( 1/scale2 )  , 'rmax [%.0g]' % ( 1/scale2 )  ) ]
    
    xmin, xmax = 1 , 7
    
    with gslCount () :
        
        for fun , nu, der in functions :
            
            ff = lambda x : fun ( nu , x )
            fd = lambda x : der ( nu , x )
            dn = D.Derivative  ( ff )
            
            title = '%s_%s' % ( fun.__name__ , nu )
            
            with use_canvas ( title ) , wait ( 1 ) :
                f1_draw ( fd , xmin = xmin , xmax = xmax , linecolor = 2 )
                f1_draw ( dn , xmin = xmin , xmax = xmax , linecolor = 4 )
                
                dd = lambda x : ( fd ( x ) - dn ( x ) ) ** 2
                d2 = I.integral ( dd , xmin = xmin , xmax = xmax ) / ( xmax - xmin ) 
                d  = math.sqrt ( d2 ) / scale1 
                
                cnt1 = SE()  
                cnt2 = SE()  
                for i in range  ( 10000 ) :
                    xx   = random.uniform ( xmin , xmax )
                    
                    v1  = fd ( xx )
                    v2  = dn ( xx ) 
                    cnt1 += abs ( v1 - v2 )                                * scale1 
                    cnt2 += abs ( v1 - v2 ) / ( abs  ( v1 ) + abs ( v2 ) ) * scale2 
                    
                dmax1  = cnt1.max()          
                dmean1 = cnt1.mean().value()  
                dmax2  = cnt2.max()          
                dmean2 = cnt2.mean().value()  
                row  = fun.__name__ , '%4s' % nu , '%.4g' % d  , \
                       '%.4g' % dmax1 , '%.4g' % dmean1 , \
                       '%.4g' % dmax2 , '%.4g' % dmean2 
                rows.append ( row )
                
            
            keep.append (   ( ff, fd , dn , dd ) ) 
        
            
    title = 'Derivatives'
    table = T.table ( rows, title = title , prefix = '# ' , alignment = 'lcccc' ) 
    logger.info ( '%s\n%s' % ( title , table ) ) 
        
        
    
# =============================================================================
if '__main__' == __name__ :

    test_math_ve1 ()
    test_math_ve2 ()
    test_math_ve3 ()
    test_math_ve4 ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
