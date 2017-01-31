#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  integral.py
#  Simple wrapper over scipy integration in a spirit of derivative.py 
#  - In case scipy is not available, it provides a reasonable replacement
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
# =============================================================================
"""Simple wrapper over scipy integration in a spirit of derivative.py 

>>> func = lambda x : x*x
>>> print integral ( func , 0 , 1 )

there are also object form:

>>> func = ...
>>> integ = Integral   ( func , 0 )
>>> print integ ( 0.1 )

In case scipy is not available, it provides a reasonable replacement
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-06"
__all__     = (
    ##
    "integral"      , ## numerical integration (as function,  using scipy if possible)
    "Integral"      , ## numerical integration (as as object, using scipy if possible)
    ##
    "romberg"       , ## numerical integration using romberg's method 
    ##
    "IntegralCache" , ## numerical integration (as object, using scipy is possible)
    ) 
# =============================================================================
import ROOT, warnings
from   ostap.math.ve   import VE
from   ostap.math.base import isequal, iszero 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.integral' )
else                       : logger = getLogger ( __name__              )
# =============================================================================

# =============================================================================
## Straw-man replacement of scipy.integrate.quad when it is not available.
#  Actually it is a primitive form of Romberg's adaptive integration
#  @see https://en.wikipedia.org/wiki/Romberg's_method
#  @see https://en.wikipedia.org/wiki/Richardson_extrapolation
#  @see https://en.wikipedia.org/wiki/Numerical_integration#Adaptive_algorithms
#  @code
#  func = lamnda x : x*x 
#  v    = romberg ( func , 0 , 1 ) 
#  @endcode 
#  CPU performance is not superb, but it is numerically stable.
def romberg ( fun                ,
              x0                 ,
              x                  ,
              err      = False   , 
              epsabs   = 1.49e-8 ,
              epsrel   = 1.49e-8 ,
              limit    = 10      , ## ignored, kept to mimic consistency with 
              args     = ()      ,
              nmax     = 8       , ## steps in Richardson's extrapolation
              maxdepth = 10        ## the maxmal depth 
              ) : 
    """Straw-man replacement of scipy.integrate.quad when it is not available.
    Actually it is a primitive form of Romberg's adaptive integration
    - see https://en.wikipedia.org/wiki/Romberg's_method
    - see https://en.wikipedia.org/wiki/Richardson_extrapolation
    - see https://en.wikipedia.org/wiki/Numerical_integration#Adaptive_algorithms
    CPU performance is not superb, but it is numerically stable.

    >>> func = lamnda x : x*x 
    >>> v    = romberg ( func , 0 , 1 ) 
    """
    # =========================================================================
    ## internal recursive function
    #  Romberg's adaptive integration wirh Richardson's extrapolation
    def _romberg_ ( f , a , b , ea , er , nmax , depth = 0 ) :
        
        rp = nmax * [0.0]
        rc = nmax * [0.0]

        ## starting value for integration step: full interval 
        h     =             (  b   -     a   )

        ## here we have R(0,0) 
        rc[0] = 0.5 * h * ( f( b ) + f ( a ) )
        
        i2    = 1
        nf    = 2  ## number of function calls
        for n in range ( 1 , nmax ) :
            
            nf     += i2     ## count number of function calls 
            rp,rc   = rc,rp  ## swap rows
            
            h      *= 0.5    ## decrement the step 

            ## use simple trapezoid rule to calculate R(n,0) 
            rr      = 0.0
            i2     *= 2 
            for k in xrange ( 1 , i2 , 2 ) :
                rr += f ( a + h * k )
                
            ## here we have R(n,0) 
            rc[0] = 0.5*rp[0] + h*rr
            
            ## calculate R(n,m) using Richardson's extrapolation
            p4 = 1 
            for m in xrange ( 1 , n + 1 ) :
                p4   *= 4 
                ## here we have R(n,m)
                rc[m] = rc[m-1] + (rc[m-1]-rp[m-1])/(p4-1)

            ## inspect last two elements in the row     
            e  = rc[n  ]
            p  = rc[n-1]
            
            ## check their absolute and relative difference 
            d  = abs ( e - p )
            ae = max ( abs ( e ) , abs ( p ) , abs( rc[0] ) ) 
            
            if d <= 0.5 * max ( ea , 0.5 * ae * er ) :
                return e , d , depth , nf       ## RETURN

        # =====================================================================
        ## Need to split the interval into subintervals
        # =====================================================================

        ## first check the maximal depth 
        if 0 < maxdepth <= depth :
            warnings.warn('Romberg integration: maximal depth is reached')
            return  e, d , depth , nf 
        
        ## prepare to split the interval into 2**n subintervals 
        n2 = 2**n

        ## the absolute uncertainty needs to be adjusted 
        new_ea = ea / n2                                ## ok
        ## keep the same relative uncertainty 
        new_er = er 

        ## the final result, its uncertainty and the max-depth 
        rr , ee , dd  = 0.0 , 0.0 , 0

        ## split the region and start recursion: 
        for i in range( n2 ) :
            
            ai = a  + i*h
            bi = ai +   h
            
            ## recursion:
            #  get result, error, depth and number of function calls from subinterval
            r , e , d , n = _romberg_ ( f , ai , bi , new_ea , new_er , nmax , depth + 1 )
            
            rr += r                ## get the integral as a sum over subintervals  
            ee += abs ( e      )   ## direct summation of uncertainties 
            dd  = max ( dd , d )   ## maximal depth from all sub-intervals 
            nf += n                ## total number of function calls 

        return rr , ee , dd , nf   ## RETURN 

    ## adjust input data
    
    a  = float ( x0     ) 
    b  = float ( x      ) 
    ea = float ( epsabs )
    er = float ( epsrel )
    
    ## the same edges 
    if isequal ( a , b )  :
        return VE( 0 , 0 ) if err else 0.0

    ## crazy settings for uncertainties 
    if  ( ea <= 0 or iszero ) and ( er <=0 or iszero( er )  ) : 
        raise AttributeError("Romberg: absolute and relative tolerances can't be non-positive simultaneously")

    ## adjust maximal number of steps in Richardson's extrapolation 
    nmax        = max ( 2 , nmax )

    ## ensure that function returnd double values: 
    func        = lambda x : float ( fun ( x , *args ) )

    ## finally use Romberg's integration 
    v, e, d, nc = _romberg_ ( func , a , b , ea , er , nmax )

    ## form the final result :
    return VE ( v , e * e ) if err else v 

# =============================================================================
try :
    from scipy import integrate
    # =============================================================================
    ## Calculate the integral (from x0 to x) for the 1D-function 
    #  @code 
    #  func = ...
    #  v    = integral ( func , 0 , 1 )
    #  @endcode 
    #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    #  @date   2014-06-06
    def integral ( fun , x0 , x      ,
                   err    = False    ,
                   **kwargs          ) :
        """Calculate the integral for the 1D-function using scipy
        
        >>> func = lambda x : x * x 
        >>> v = integral(func,0,1)
        """
        func   = lambda x : float ( fun ( x ) ) 
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = integrate.quad ( func , x0 , x , **kwargs )
            return VE( result[0] , result[1] * result[1] ) if err else result[0]
        
except ImportError :
    logger.warning ("scipy.integrate is not available, use local ``romberg''-replacement")
    ## use romberg integration as default method when scipy is not available 
    integral = romberg
# =============================================================================


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
    def __init__ ( self , func , xlow = 0 , err = False , **args ) :
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
    def _integrate_ ( self , xmn , xmx , **args ) :
        args   = args if args else self._args
        ## try : 
        return integral ( self._func , xmn , xmx , self._err , **args )
        ##except :
        ##    print 'EXCEPT' , xmn, xmx , type(xmn), type(xmx)
            
    ## Calculate the integral for the 1D-function using scipy
    def __call__ ( self , x , *args , **kwargs ) :
        """Calculate the integral for the 1D-function using scipy
        
        >>> func = ...
        >>> func_0 = Integral(func,0)
        >>> func_0 ( 10 ) 
        """
        a         = kwargs.copy()
        a['args'] = args 
        return self._integrate_ ( self._xmin , x , **a )

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
    def __init__ ( self , func , xlow = 0 , err = False , **args ) :
        """Calculate the integral for the 1D-function using scipy
        
        >>> func = ...
        >>> func_0 = Integral(func,0)
        """
        Integral.__init__ ( self , func , xlow , err , *args )
        self._prev   = None 
        
    ## Calculate the integral for the 1D-function using scipy
    def __call__ ( self , x , **args ) :
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
        result = self._integrate_ ( xmn , x , **args )
        result += dlt 
        
        if not args and isinstance ( x , float ) :
            self._prev =  x , result 
            
        return result 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
