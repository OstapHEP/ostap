#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file funstats.py 
#  Utilities to get moments for various functions/distributions/pdfs
"""
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-06"
__all__     = (
    ##
    'FunARGS'     , ## Helper base class to keep addtional arguments for function
    'FunBASE1D'   , ## Helper base class to keep addtional arguments & range for 1D-function
    'FunBASE2D'   , ## Helper base class to keep addtional arguments & range for 2D-function
    'FunBASE3D'   , ## Helper base class to keep addtional arguments & range for 3D-function
    ##
    'FunMean'     , ## mean value for function values over the 1D-interval
    'FunVariance' , ## variance  for function values over the 1D-interval
    'FunRms'      , ## RMS for function values over the 1D-interval
) 
# =============================================================================
from   ostap.core.ostap_types import integer_types, num_types
from   ostap.math.base        import pos_infinity, neg_infinity
from   ostap.utils.basic      import typename
import math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.funstats' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
## @calss FuARGS
#  Helper base class to keep addtional arguments for function
class FunARGS(object) :
    """ Helper base class to keep addtional arguments for function
    """
    def __init__ ( self , *args ) :
        self.__args   = args
    @property
    def args ( self ) :
        """`args` : additional positional arguments for function call
        """
        return self.__args

# =============================================================================
## @class FunBASE1D
#  Helper base class to keep xmin/xmax/args/kwargs for 1D function 
class FunBASE1D(FunARGS) :
    """ Helper base class to keep xmin/xmax/args/kwargs for 1D function 
    """
    def __init__ ( self , xmin , xmax , *args ) :
        FunARGS.__init__ ( self , *args )
        assert isinstance ( xmin , num_types ) , "Invalid `xmin' type" % typename ( xmin )
        assert isinstance ( xmax , num_types ) , "Invalid `xmax' type" % typename ( xmax )
        xmin = float ( xmin )
        xmax = float ( xmax )        
        self.__xmin = min ( xmin , xmax ) 
        self.__xmax = max ( xmin , xmax )

    @property
    def  xmin ( self ) :
        "`xmin'- low edge of the X-interval"
        return self.__xmin
    @property
    def  xmax ( self ) :
        "`xmax'- high edge of the X-interval"
        return self.__xmax

# =============================================================================
## @class FunBASE2D
#  Helper base class to keep xmin/xmax/ymin/ymax/args/kwargs for 2D function 
class FunBASE2D(FunBASE1D) :
    """ Helper base class to keep xmin/xmax/args/kwargs for 1D function 
    """
    def __init__ ( self , xmin , xmax , ymin, ymax , *args , **kwargs ) :
        FunBASED1D.__init__ ( self , xmin , xmax , *args , **kwargs )
        assert isinstance ( ymin , num_types ) , "Invalid `ymin' type" % typename ( ymin )
        assert isinstance ( ymax , num_types ) , "Invalid `ymax' type" % typename ( ymax )
        ymin = float ( ymin )
        ymax = float ( ymax )        
        self.__ymin = min ( ymin , ymax ) 
        self.__ymax = max ( ymin , ymax )
        
    @property
    def  ymin ( self ) :
        "`ymin'- low edge of the Y-interval"
        return self.__ymin
    @property
    def  ymax ( self ) :
        "`ymax'- high edge of the Y-interval"
        return self.__ymax

# =============================================================================
## @class FunBASE3D
#  Helper base class to keep xmin/xmax/ymin/ymax/args/kwargs for 2D function 
class FunBASE3D(FunBASE2D) :
    """ Helper base class to keep xmin/xmax/args/kwargs for 1D function 
    """
    def __init__ ( self , xmin , xmax , ymin, ymax , zmin , zmax , *args , **kwargs ) :
        FunBASED2D.__init__ ( self , xmin , xmax , ymin , ymax , *args , **kwargs )
        assert isinstance ( zmin , num_types ) , "Invalid `zmin' type" % typename ( zmin )
        assert isinstance ( zmax , num_types ) , "Invalid `zmax' type" % typename ( zmax )
        zmin = float ( zmin )
        zmax = float ( zmax )        
        self.__zmin = min ( zmin , zmax ) 
        self.__zmax = max ( zmin , zmax )
        
    @property
    def  zmin ( self ) :
        "`zmin'- low edge of the Z-interval"
        return self.__zmin
    @property
    def  zmax ( self ) :
        "`zmax'- high edge of the Z-interval"
        return self.__zmax

# ============================================================================
## @class FunMean
#  Mean-value of function over 1D-interval
#  @code
#  mean  = FunMean( fun ,xmin,xmax) 
#  value = mean  ( math.sin )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class FunMean(FunBASE1D) :
    """ Calculate the mean-value of function over the 1D-interval
    >>> mean  = FunMean( fun ,xmin,xmax) 
    >>> value = mean  ( math.sin )
    """
    # =========================================================================
    ## integrate the function between xmin and xmax 
    def integral ( self , func , xmin , xmax , *args ) :
        """ Integrate the function between xmin and xmax"""
        args  = args if args else self.args
        from ostap.math.integral import Integral as II  
        integrator = II ( func , xmin , err = False , args = args )
        return integrator ( xmax , *args )

    ## mean value of function over the interval 
    def fun_mean ( self , func , xmin , xmax , *args ) :
        args  = args if args else self.args
        value = self.integral ( func , xmin, xmax , *args )
        return value / ( xmax - xmin ) 
                               
    def __call__ ( self , func , *arg ) :
        args  = args if args else self.args        
        return self.fun_mean ( func , self.xmin , self.xmax , *args )

# ==============================================================================
## @class FunVariance
#  Calculate variance of the function over 1D-interval
class FunVariance(FunMean) :
    """ Calculate variance of the function over 1D-interval
    """
    ## Calculate variance of the function over interval
    def fun_variance ( self , func , xmin , xmax , *args ) :
        """ Calculate variance of the function over interval
        """        
        args  = args if args else self.args        

        ## get mean value
        fn_mean = self.fun_mean ( func , *args )

        ## delta squared 
        def f2 ( x , *args ) :
            fv = func ( x , *args ) - fn_mean 
            return fv * fv
        
        value = self.integral ( f2 , xmin , xmax , *args )
        return value / ( xmax - xmin )
    
    def __call__ ( self , func , *arg ) :
        args  = args if args else self.args                
        return self.fun_variance ( func , self.xmin , self.xmax , *args )
    
# ==============================================================================
## @class FunRms
#  Calculate Rms of the function over 1D-interval
class FunRms(FunVariance) :
    """ Calculate Rms of the function over the 1D-interval
    """
    def __call__ ( self , func , *arg ) :
        args  = args if args else self.args                
        value = self.fun_variance ( func , self.xmin , self.xmax , *args )
        return math.sqrt ( value )
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
