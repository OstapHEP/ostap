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
# =============================================================================
if '__main__' == __name__ :

    test_integral ()
    
# =============================================================================
# The END 
# =============================================================================
