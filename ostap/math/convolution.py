#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/convolution.py
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-02-27
# =============================================================================
""" Module with some useful utilities for dealing with BSplines
- control_polygon      : get a control polygon for BSpline
- upper_convex_hull    : upper convex hull for BSpline
- lower_convex_hull    : lower convex hull for BSpline
- convex_hull          :       convex hull for BSpline
- crossing_points      : get crossing points of control polygon with x-axis
- solve                : solve equation B(x) = C
- interpolate          : construct interpolating B-spline 
- approximate          : construct approximating B-spline
- generate&shoot       : generate random numbers 
#"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = (
    'control_polygon'      , ## get a control polygon for BSpline
    'upper_convex_hull'    , ## upper convex hull for BSpline
    'lower_convex_hull'    , ## lower convex hull for BSpline
    'convex_hull'          , ##       convex hull for BSpline
    'crossing_points'      , ## get crossing points of control polygon with x-axis
    'solve'                , ## solve equation B(x) = C
    'interpolate'          , ## spline interpolation
    'approximate'          , ## variation diminishing approximation 
    'generate'             , ## generate random numbers 
    'shoot'                , ## generate random numbers 
    )
# =============================================================================
import  ROOT, math
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.convolution' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================

def convolution  ( fun  , xmin  , xmax  , npoints ,
                   fun2 , xmin2 , xmax2           ) :

    import numpy as np
    
    xmn     = float ( xmin )
    xmx     = float ( xmin )    
    dx      = ( xmx - xmn ) / npoints

    vfun1   = np.vectorize ( fun )
    vfun2   = np.vectorize ( fun )

    ps1     = np.arange ( xmn  , xmx  + dx , dx )
    ps2     = np.arange ( xmn2 , xmx2 + dx , dx )

    fval1   = vfun1 ( ps1 )
    fval2   = vfun1 ( ps1 )

    cnv     = np.convolve ( fval1 , fval2 , 'same' )

    
    
    
    


# =============================================================================
##                                                                      The END 
# =============================================================================
