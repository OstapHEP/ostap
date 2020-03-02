#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/sp_interpolation.py
#  Module with some useful scipy-based utilities for dealing with interpolation.
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2020-02-28
# =============================================================================
"""Useful scipy-based utilities for dealing with spline interpolation.
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2020-02-28"
__all__     = (
    )
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.sp_interpolation' )
else                       : logger = getLogger ( __name__                      )
# =============================================================================
import scipy.interpolate as SI

# ==============================================================================
## simple class for scipy-based interpolation  
class SplineInterpolator(object) :
    """Simple class for scipy-based interpolation
    >>> x = ... ##
    >>> y = ... ##
    >>> sp = SplineInterpolator ( (x,y) , 3 )
    """
    # ===========================================================================
    ## Create the interpolator
    #  @code
    #  >>> x = ... ##
    #  >>> y = ... ##
    #  >>> sp = SplineInterpolator ( (x,y) , 3 )
    #  @endcode
    def __init__ ( self , data , order = 3 , ext = 0 ) :
        """ Create the interpolator
        >>> x = ... ##
        >>> y = ... ##
        >>> sp = SplineInterpolator ( ( x , y ) , 3 )
        """

        x , y = data

        ## create the interpolation spline 
        self.__spline = SI.InterpolatedUnivariateSpline ( x   = x     ,
                                                          y   = y     ,
                                                          k   = order ,
                                                          ext = ext   )
        self.__antiderivative = None
        self.__derivative     = None

    # =========================================================================
    ## the main method:
    #  @code
    #  >>> spline = ...
    #  >>> result = spline ( x ) 
    #  @endcode 
    def __call__ ( self , x ) :
        return self.__spline ( x )

    @property
    def derivative ( self ) :
        """``derivative'' : get the spline object for defivative"""
        if not self.__derivative :
            self.__derivative = self.__spline.derivative()
        return self.__derivative

    @property
    def antiderivative ( self ) :
        """``antiderivative'' : get the spline object for indefinite itegral"""
        if not self.__antiderivative :
            self.__antiderivative = self.__spline.antiderivative()
        return self.__antiderivative
    
    # =========================================================================
    ## get the integral
    #  @code 
    #  >>> spline = ...
    #  >>> o = spline.integral()        ## indefinite integral
    #  >>> r = spline.integral( a , b ) ##   definite integral form a to b
    #  @endcode 
    def integral ( self , a = None , b = None ) :
        """Get the integral
        >>> spline = ...
        >>> o = spline.integral()        ## indefinite integral
        >>> r = spline.integral( a , b ) ##   definite integral form a to b         
        """
        if   a is None and b is None : return self.antiderivative

        assert ( not a is None ) and ( not b is None ) ,\
               'lower and upper edged must be splicified!'

        return self.__spline.integral ( a , b )  
        
    @property
    def spline ( self ) :
        """``spline'' : get the underlying bspline/scipy object"""
        return self.__spline

    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
