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
    'FunMEAN'     , ## mean value for function values over the 1D-interval
    'FunVariance' , ## variance  for function values over the 1D-interval
    'FunRMS'      , ## RMS for function values over the 1D-interval
    'FunMNMX'     , ## MINMAX for function using the brute force approach
)
# =============================================================================
from   ostap.core.ostap_types import integer_types, num_types
from   ostap.math.base        import pos_infinity, neg_infinity, isequal
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
    def __init__ ( self , *args , **kwargs ) :
        self.__args   = args
        self.__kwargs = kwargs
    @property
    def args ( self ) :
        """`args` : additional positional arguments for function call
        """
        return self.__args
    def kwargs ( self ) :
        """`kwargs` : additional keyword arguments for function call
        """
        return self.__kwargs

# =============================================================================
## @class FunBASE1D
#  Helper base class to keep xmin/xmax/args/kwargs for 1D function 
class FunBASE1D(FunARGS) :
    """ Helper base class to keep xmin/xmax/args/kwargs for 1D function 
    """
    def __init__ ( self , xmin , xmax , *args , **kwargs ) :
        FunARGS.__init__ ( self , *args , **kwargs )
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
## @class FunMEAN
#  Mean-value of function over 1D-interval
#  @code
#  mean  = FunMEAN( fun ,xmin,xmax) 
#  value = mean  ( math.sin )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class FunMEAN(FunBASE1D) :
    """ Calculate the mean-value of function over the 1D-interval
    >>> mean  = FunMean( fun ,xmin,xmax) 
    >>> value = mean  ( math.sin )
    """
    # =========================================================================
    ## integrate the function between xmin and xmax 
    def integral ( self , func , xmin , xmax , *args , **kwargs ) :
        """ Integrate the function between xmin and xmax"""
        args   =   args if   args else self.args
        kwargs = kwargs if akwrgs else self.kwargs         
        from ostap.math.integral import Integral as II  
        integrator = II ( func , err = False , args = args , kwargs = kwargs )
        return integrator.integrate ( xmin , xmax )

    ## mean value of function over the interval 
    def fun_mean ( self , func , xmin , xmax , *args , **kwargs ) :
        value = self.integral ( func , xmin, xmax , *args , **kwargs )
        return value / ( xmax - xmin ) 
                               
    def __call__ ( self , func , *arg , **kwrgs) :
        return self.fun_mean ( func , self.xmin , self.xmax , *args , **kwargs )

# ==============================================================================
## @class FunVariance
#  Calculate variance of the function over 1D-interval
class FunVariance(FunMEAN) :
    """ Calculate variance of the function over 1D-interval
    """
    ## Calculate variance of the function over interval
    def fun_variance ( self , func , xmin , xmax , *args , **kwargs ) :
        """ Calculate variance of the function over interval
        """        
        ## get mean value
        fn_mean = self.fun_mean ( func , *args , **kwargs )

        ## delta squared 
        def f2 ( x , *a , **kw ) :
            fv = func ( x , *a , **kw ) - fn_mean 
            return fv * fv
        
        value = self.integral ( f2 , xmin , xmax , *args , **kwargs )
        return value / ( xmax - xmin )
    
    def __call__ ( self , func , *arg , **kwargs ) :
        return self.fun_variance ( func , self.xmin , self.xmax , *args , **kwargs )
    
# ==============================================================================
## @class FunRMS
#  Calculate Rms of the function over 1D-interval
class FunRMS(FunVariance) :
    """ Calculate RMS of the function over the 1D-interval
    """
    def __call__ ( self , func , *args , **kwargs ) :
        value = self.fun_variance ( func , self.xmin , self.xmax , *args , **kwargs )
        return math.sqrt ( value )
    
# =============================================================================
## Estimate the min/max for the function using brute-force approach
#  @attentoon the approximaiton coudl be rather rude 
class FunMNMX(FunRMS) :
    """ Estimate the min/max for the function using brute-force approach
    - attentoon the approximation could be rather rude 
    """
    
    def __init__ ( self , xmin , xmax , N = 100 , NMAX = 100000 , *args , **kwargs ) :
        
        super(FunMNMX,self).__init__  ( xmin , xmax , *args , **kwargs )
        
        self.__N    = N 
        self.__NMAX = NMAX 
        
    def __call__ ( self , func , *args ) :
        
        args   = args   if args   else self.args 
        kwargs = kwargs if kwargs else self.kwargs                
        
        from ostap.utils.ranges import vrange 
        from ostap.utils.utils  import memoize 
        
        memfun = memoize ( func )
                
        IMAX = 15  
        NMAX = 100000
        
        N = self.__N 
        
        xmin, xmax = self.xmin , self.xmax 
        vmn0 = min ( memfun ( x , *args , **kwargs ) for x in vrange ( xmin , xmax , N ) ) 
        vmx0 = max ( memfun ( x , *args , **kwargs ) for x in vrange ( xmin , xmax , N ) ) 
        
        vmin = vmn0 
        vmax = vmx0 
        
        A, B = 2**5 , 1.0 / ( 2**5 - 1 )
        
        for i in range ( IMAX ):
            
            N *= 2
            if self.__NMAX < N : break 
            
            vmn  = min ( memfun ( x , *args , **kwargs ) for x in vrange ( xmin , xmax , N ) ) 
            vmx  = max ( memfun ( x , *args , **kwargs ) for x in vrange ( xmin , xmax , N ) )
        
            ## a kind of Richardson's extrapolation    
            vmin = ( A * vmn - vmn0 ) * B  
            vmax = ( A * vmx - vmx0 ) * B 
             
            if 2 < i :
                mnok = isequal ( vmin , vmn ) or isequal ( vmn , vmn0 )
                if mnok : 
                    mxok = isequal ( vmax , vmx ) or isequal ( vmx, vmx0 )
                    if mxok : return vmin , vmax 
             
            vmn0 = vmn 
            vmx0 = vmx 
            
        return vmin, vmax 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
