#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file ostap/math/tests/test_math_comvolution.py
# Test module for ostap/math/param.py
# - It tests parameterisation of functions
# ============================================================================= 
""" Test module for ostap/math/param.py
- It tests parameterisations of functions
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import 
# ============================================================================= 
from   ostap.plotting.canvas  import use_canvas
from   ostap.math.models      import f1_draw 
import ROOT, math
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'ostap.test_math_convolution' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( 'Test for function convolution')
# =============================================================================

# =============================================================================
def test_convolution() :

    logger = getLogger ( 'test_convolution' )
    logger.info ( 'Test convolution' ) 

    try :
        # =====================================================================
        from ostap.math.sp_convolution import GaussConvolution
        # =====================================================================
    except ImportError :
        # =====================================================================
        logger.error ( "Gauss convolution is not available!" )
        return

    fun = lambda x : x - math.floor ( x )

    g1 = GaussConvolution ( fun , xmin = 1 , xmax = 5 , sigma = 0.05 , N = 2**14 )
    with use_canvas ( 'Gauss convolution with sigma=0.05' , wait = 3 ) :
        f1_draw ( g1 )

    g2 = GaussConvolution ( fun , xmin = 1 , xmax = 5 , sigma = 0.10 , N = 2**14 )
    with use_canvas ( 'Gauss convolution with sigma=0.10' , wait = 3 ) :
        f1_draw ( g2 )

    g3 = GaussConvolution ( fun , xmin = 1 , xmax = 5 , sigma = 0.15 , N = 2**14 )
    with use_canvas ( 'Gauss convolution with sigma=0.15' , wait = 3 ) :
        f1_draw ( g3 )

    g4 = GaussConvolution ( fun , xmin = 1 , xmax = 5 , sigma = 0.20 , N = 2**14 )
    with use_canvas ( 'Gauss convolution with sigma=0.20' , wait = 3 ) :
        f1_draw ( g4 )
        
    
# ============================================================================
if '__main__' == __name__ :

    test_convolution  () 
                

                  
# =============================================================================
##                                                                      The END 
# =============================================================================
                  
