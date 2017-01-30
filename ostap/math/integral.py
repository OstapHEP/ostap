#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
    _scipy = True 
except ImportError :
    _scipy = False 
    logger.warning ('scipy.integrate is not available')
    

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
    from scipy import integrate
    func   = lambda x : float ( fun ( x ) ) 
    result = integrate.quad ( func , x0 , x , **kwargs )
    return VE( result[0] , result[1] * result[1] ) if err else result[0]

# =============================================================================
## straw-man replacement of scipy.integrate.quad when it is not available
#  actually it is a primitive form of Romberg's adaptive intergation
#  @see https://en.wikipedia.org/wiki/Romberg's_method
def integral_romberg ( fun              ,
                       x0               ,
                       x                ,
                       err    = False   , 
                       epsabs = 1.49e-8 ,
                       epsrel = 1.49e-8 ,
                       limit  = 10      , ## ignored.... 
                       args   = ()      
                       ) : 
    """Straw-man replacement of scipy.integrate.quad when it is not available
    actually it is a primitive form of Romberg's adaptive intergation
    - see https://en.wikipedia.org/wiki/Romberg's_method
    """
    _nmax   = 6 
    def _romberg_ ( f , a , b , ea , er ) :
        
        rp = _nmax * [0.0]
        rc = _nmax * [0.0]
        
        h     =             (  b   -     a   )

        ## here we have R(0,0) 
        rp[0] = 0.5 * h * ( f( b ) + f ( a ) ) 

        for n in range ( 1 , _nmax ) :
            
            rr   = 0.0
            h   *= 0.5 
            for k in xrange ( 1 , 2**( n - 1 ) + 1 ) :
                rr += f ( a + h * ( 2 * k - 1 ) )

            ## here we have R(n,0) 
            rc[0] = 0.5*rp[0] + h*rr
            
            for m in range ( 1 , n + 1 ) :
                ## here we have R(n,m) 
                rc[m] = rc[m-1] + (rc[m-1]-rp[m-1])/(4**m-1)

            ## inspect last two elements in the row    
            e = rc[n  ]
            p = rc[n-1]
                
            d  = abs ( e - p )
            ae = max ( abs ( e ) , abs ( p ) )

            ## check the relative error 
            if ae and d <= er * ae :
                return e , 0.5*d

            ## check the absolute error 
            if        d <= ea      :
                return e , 0.5*d 

            ## swap rows   
            rp,rc = rc,rp 

        ## no result..  split the range 
            
        new_ea = ea / 2**n ## ok 
        new_er = er * 2    ## ?
        
        ## the final result and uncertainty:
        rr,ee  = 0.0, 0.0

        ## split the region and iterate 
        for i in range( 2**n ) :
            
            ai = a  + i*h
            bi = ai +   h 

            r, e = _romberg_ ( f , ai , bi , new_ea , new_er ) 
            rr += r
            ee += abs(e) ## direct summation of uncertaities 
            
        return rr , ee 
        
    v , e = _romberg_ ( lambda x : float( fun ( x  , *args ) ) ,
                        float (x0)      , float(x)             ,
                        float (epsabs)  , float(epsrel) )      
    
    return VE( v , e * e ) if err else float(v) 

# =============================================================================
## use romberg integration as default method when scipy is not available 
if not _scipy :
    integral = integral_romberg 
 
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
