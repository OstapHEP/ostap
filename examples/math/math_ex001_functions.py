#!/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright Ostap developers
# =============================================================================
## @file examples/math/math_ex001_functions.py
#  Example of basic operations with VE objects and functions
# =============================================================================
"""Example of basic operations with VE objects and functions
"""
# =============================================================================
from   __future__ import print_function
import math 
from   ostap.math.ve         import VE
from   ostap.math.math_ve    import *  ## get functions 

a =  VE ( 2 , 0.2**2 )
b =  VE ( 3 , 0.3**2 )

print ( 'a,b  : %s   %s' % (  a ,  b ) )
print ( 'a+b  : %s '     % (  a +  b ) )
print ( 'a-b  : %s '     % (  a -  b ) ) 
print ( 'a*b  : %s '     % (  a *  b ) )
print ( 'a/b  : %s '     % (  a /  b ) )
print ( 'a**b : %s '     % (  a ** b ) )

print ( 'a+2  : %s '     % (  a +  2 ) )
print ( 'a-2  : %s '     % (  a -  2 ) )
print ( 'a*2  : %s '     % (  a *  2 ) )
print ( 'a/2  : %s '     % (  a /  2 ) )
print ( 'a**2 : %s '     % (  a ** 2 ) )

print ( '2+a  : %s '     % (  2 +  a ) )
print ( '2-a  : %s '     % (  2 -  a ) )
print ( '2*a  : %s '     % (  2 *  a ) )
print ( '2/a  : %s '     % (  2 /  a ) )
print ( '2**a : %s '     % (  2 ** a ) )

print ( 'a+a  : %s '     % (  a +  a ) ) 
print ( 'a-a  : %s '     % (  a -  a ) ) 
print ( 'a*a  : %s '     % (  a *  a ) )
print ( 'a/a  : %s '     % (  a /  a ) )
print ( 'a**a : %s '     % (  a ** a ) )

funcs = [ exp    , expm1  ,
          log    , log10  , log1p  ,
          sqrt   , cbrt   ,
          sin    , cos    , tan    ,
          sinh   , cosh   , tanh   , sech   ,
          asin   , acos   , atan   ,
          asinh  , acosh  , atanh  ,
          erf    , erfc   , erfi   , erfcx  ,
          probit ,
          gamma  , tgamma , lgamma , igamma ]

from ostap.math.derivative import EvalVE
funcs += [ EvalVE ( math.sin , math.cos ) ,  ## use function and known derivative 
           EvalVE ( math.sin )            ]  #  use numerical derivative

v  = 0.25 * a 
print ( 'x=%s'  % v )
for f in funcs :
    print ( "\t%s\t%s = %s " % ( f.__name__ , v ,  f ( v ) ) )
    
# =============================================================================
# The END 
# =============================================================================
