#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_rootfinder.py
#  Test module for the file ostap/math/rotfinder.py
# ============================================================================= 
""" Test module for ostap/math/rootfinder.py.
"""
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_rootfinder' ) 
else                       : logger = getLogger ( __name__          )
# =============================================================================
import random,math 
import ostap.math.models 
from   ostap.math.rootfinder    import find_root, findroot 
# =============================================================================


def test_root_sin () :

    fun  = lambda  x : -1.00 * math.sin ( 0.5 * ( x - 1.0 ) ) 
    der1 = lambda  x : -0.50 * math.cos ( 0.5 * ( x - 1.0 ) ) 
    der2 = lambda  x : +0.25 * math.sin ( 0.5 * ( x - 1.0 ) ) 
    
    logger.info ( 'sin/Halley:  %.15g\n%s' % find_root ( fun , -0.5 , 2.5   ,
                                                         deriv1      = der1 ,
                                                         deriv2      = der2 ,
                                                         full_output = True ) )
    logger.info ( 'sin/Newton:  %.15g\n%s' % find_root ( fun , -0.5 , 2.5   ,
                                                         deriv2      = der1 ,
                                                         full_output = True ) )
    logger.info ( 'sin/Plain :  %.15g\n%s' % find_root ( fun , -0.5 , 2.5   ,
                                                         full_output = True ) )
    logger.info ( 'sin/Brent :  %.15g\n%s' % findroot  ( fun , -0.5 , 2.5   ,
                                                         full_output = True ) )
    
def test_root_mult () :

    K = 1
    N = 2*K + 1 
    
    fun  = lambda  x : 1.0                * ( x - 1.0 ) **   N 
    der1 = lambda  x : 1.0 * N            * ( x - 1.0 ) ** ( N - 1 )
    der2 = lambda  x : 1.0 * N * ( N -1 ) * ( x - 1.0 ) ** ( N - 2 )
    
    logger.info ( 'mult/Halley:  %.15g\n%s' % find_root ( fun , -1 , 10       ,
                                                          deriv1      = der1 ,
                                                          deriv2      = der2 ,
                                                          full_output = True  ,
                                                          disp        = False ) )
    logger.info ( 'mult/Newton:  %.15g\n%s' % find_root ( fun , -1 , 10       ,
                                                          deriv1      = der1 ,
                                                          full_output = True  ,
                                                          disp        = False ) )
    logger.info ( 'mult/Plain :  %.15g\n%s' % find_root ( fun , -1 , 10       ,
                                                          full_output = True  ,
                                                          disp        = False ) )
    logger.info ( 'mult/Brent :  %.15g\n%s' % findroot  ( fun , -1 , 10       ,
                                                          full_output = True  ,
                                                          disp        = False ) )
    

# =============================================================================
if '__main__' == __name__ :
    
    test_root_sin  () 
    test_root_mult () 
    
    
# =============================================================================
# The END 
# =============================================================================
