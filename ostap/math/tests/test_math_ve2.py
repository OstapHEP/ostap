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

    logger = getLogger("test_math_ve2")


    funcs = ( bessel_J , bessel_Y ,bessel_I , bessel_K )

    rows = [ ( 'Function' , 'order' ,  'argument' , 'value' ) ]
    
    for f in funcs :

        for n in range ( 0 , 4 ) : 
            x   = random.uniform ( 0 , 1 )
            row = f.__name__ , '%s' % n , '%s' % x , '%s' % f ( n , x )        
            rows.append ( row )

            x   = VE ( random.uniform ( 1 , 2 ) , 0.1 **2 ) 
            row = f.__name__ , '%s' % n , '%s' % x , '%s' % f ( n , x )        
            rows.append ( row )

        n   = random.uniform ( -1 , 3 )
        row = f.__name__ , '%+.3f' % n , '%s' % x , '%s' % f ( n , x )        
        rows.append ( row )
        
        x   = VE ( random.uniform ( 1 , 2 ) , 0.1 **2 ) 
        row = f.__name__ , '%+.3f' % n , '%s' % x , '%s' % f ( n , x )        
        rows.append ( row )

    title = 'Functions'
    table = T.table ( rows, title = title , prefix = '# ' , alignment = 'lccc' ) 
    logger.info ( '%s\n%s' % ( title , table ) ) 

    
# =============================================================================
if '__main__' == __name__ :

    test_math_ve1 ()
    test_math_ve2 ()
    test_math_ve3 ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
