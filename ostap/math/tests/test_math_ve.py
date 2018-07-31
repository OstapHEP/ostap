#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
""" Test module for ostap/math/math_ve.py.
"""
# ============================================================================= 
import math

from ostap.math.ve      import VE
from ostap.math.base    import cpp, iszero, isequal
from ostap.math.math_ve import *

from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_math_ve' )
else                       : logger = getLogger ( __name__             )

vars  = [ VE ( 0.001 , 0.0001**2 ) , VE(1,0) , VE(1,0.1**2) , VE(10,0.01**2) ]

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

## use helper object:
from ostap.math.derivative import EvalVE
funcs += [ EvalVE ( math.sin , math.cos ) ,
           EvalVE ( math.sin )            ]

def test_math_ve():
    for v in vars :
        logger.info ( 'Var = %s ' % v )
        for f in funcs :
            logger.info ( "\t%s\t%s = %s " % ( f.__name__ , v ,  f(v) ) )

            
# =============================================================================
if '__main__' == __name__ :

    test_math_ve()
    
# =============================================================================
# The END 
# =============================================================================
