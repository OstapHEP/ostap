#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file  integral.py
#  Simple wrapper over scipy integration  in a spirit of derivative.py 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
# =============================================================================
"""Simple wrapper over scipy integration  in a spirit of derivative.py 

>>> func = lambda x : x*x
>>> print integral ( func , 0 , 1 )

there are also object form:

>>> func = ...
>>> integ = Integral   ( func , 0 )
>>> print integ ( 0.1 )  
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-06"
__all__     = (
    ##
    "integral"      , ## numerical integration     (as function,  using scipy)
    "Integral"      , ## numerical integration     (as as object, using scipy)
    ##
    "IntegralCache" , ## numerical integration     (as object, using scipy and cache)
    ) 
# =============================================================================
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.integral' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
try :
    from scipy import integrate
except ImportError :
    logger.warning ('scipy.integrate is not availabe')

# =============================================================================
## Calculate the integral (from x0 to x) for the 1D-function 
#  @code 
#  func = ...
#  v    = integral ( func , 0 , 1 )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
def integral ( fun , x0 , x , err = False , *args ) :
    """Calculate the integral for the 1D-function using scipy
    
    >>> func = lambda x : x * x 
    >>> v = integral(func,0,1)
    """
    from scipy import integrate
    func   = lambda x : float ( fun ( x ) ) 
    result = integrate.quad ( func , x0 , x , args = args )
    return VE( result[0] , result[1] * result[1] ) if err else result[0] 

# =============================================================================
## @class Integral
#  Calculate the integral (from x0 to x) for the 1D-function 
#  @code 
#  func  = lambda x : x * x 
#  iint  = Integral ( func , 0 ) ## specify x_low 
#  value = iiint (  10  )        ## specify x_high 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class Integral(object) :
    """Calculate the integral for the 1D-function using scipy
    
    >>> func  = lambda x : x * x      ## define function 
    >>> iint  = Integral ( func , 0 ) ## specify x_low 
    >>> value = iint (  10  )         ## specify x_high 
    """
    ## Calculate the integral for the 1D-function using scipy
    def __init__ ( self , func , xlow = 0 , err = False , *args ) :
        """Calculate the integral for the 1D-function using scipy
        
        >>> func   = ...
        >>> func_0 = Integral(func,0)
        >>> value  = func_
        """
        self._func   = func 
        self._xmin   = float ( xlow ) if isinstance ( xlow , (int,long) ) else xlow 
        self._err    = err
        self._args   = args
        
    ## Calculate the integral for the 1D-function using scipy
    def _integrate_ ( self , xmn , xmx , *args ) :
        args   = args if args else self._args
        ## try : 
        return integral ( self._func , xmn , xmx , self._err , *args )
        ##except :
        ##    print 'EXCEPT' , xmn, xmx , type(xmn), type(xmx)
            
    ## Calculate the integral for the 1D-function using scipy
    def __call__ ( self , x , *args ) :
        """Calculate the integral for the 1D-function using scipy
        
        >>> func = ...
        >>> func_0 = Integral(func,0)
        >>> func_0 ( 10 ) 
        """
        return self._integrate_ ( self._xmin , x , *args )

# =============================================================================
## @class IntegralCache
#  Calculate the integral (from x0 to x) for the 1D-function 
#  @code 
#  >>> func  = lambda x : x * x 
#  >>> iint  = IntegralCache(func,0) ## specify x_low 
#  >>> value = iint  ( 10 )          ## specify x_hgh 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class IntegralCache(Integral) :
    """Calculate the integral for the 1D-function using scipy
    >>> func   = lambda x : x*x 
    >>> iint   = IntegralCache ( func , 0 ) ## specify x_low 
    >>> value  = iint ( 10 )                ## specify x_high
    """
    ## Calculate the integral for the 1D-function using scipy
    def __init__ ( self , func , xlow = 0 , err = False , *args ) :
        """Calculate the integral for the 1D-function using scipy
        
        >>> func = ...
        >>> func_0 = Integral(func,0)
        """
        Integral.__init__ ( self , func , xlow , err , *args )
        self._prev   = None 
        
    ## Calculate the integral for the 1D-function using scipy
    def __call__ ( self , x , *args ) :
        """Calculate the integral for the 1D-function using scipy
        
        >>> func = ...
        >>> func_0 = Integral(func,0)
        >>> func_0 ( 10 ) 
        """
        x      = float ( x ) if isinstance ( x , ( int , long ) ) else x  
        args   = args if args else self._args

        xmn = self._xmin
        dlt = 0
        
        ## check previos calculations 
        if not args  and isinstance ( x , float ) and self._prev :
            if       isinstance ( self._xmin , float ) and ( self._prev[0] - x ) < abs ( self._xmin - x ) :
                xmn = self._prev[0]
                dlt = self._prev[1]
            elif not isinstance ( self._xmin , float ) :                        
                xmn = self._prev[0]
                dlt = self._prev[1]
                
        ## use scipy
        result = self._integrate_ ( xmn , x , *args )
        result += dlt 
        
        if not args and isinstance ( x , float ) :
            self._prev =  x , result 
            
        return result 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.logger.line import line 
    logger.info ( __file__  + '\n' + line  ) 
    logger.info ( 80*'*'   )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 


# =============================================================================
# The END 
# =============================================================================
