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
## __all__     = (
##    )
# =============================================================================
import  ROOT, math
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.convolution' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================

try :
    from ostap.math.sp_convolution import ( ArrayConvolution ,
                                            GaussConvolution , Convolution )
    __all__  = ( 'ArrayConvolution' ,
                 'GaussConvolution' ,
                 'Convolution'      )
    
except ImportError :
    pass
    
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
