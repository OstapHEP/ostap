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

>>> func1 = lambda x : x*x
>>> print integral  ( func1 , 0 , 1 )

>>> func2 = lambda x,y: x*x+y*y
>>> print integral2 ( func2 , 0 , 1 , 0 , 2 )

>>> func3 = lambda x,y,z: x*x+y*y+z*z
>>> print integral3 ( func3 , 0 , 1 , 0 , 2 , 0 , 3 )

there are also object form:

>>> func = ...
>>> integ = Integral   ( func , 0 )
>>> print integ ( 0.1 )

In case scipy is not available, it provides reasonable replacements
for 1D,2D and 3D integration:
 - 1D    integration using Romberg's method
 - 2D&3D integration using Genz&Malik's method 

"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-06"
__all__     = (
    ##
    "integral"       , ## (1D) numerical integration (as function, using scipy if possible)
    "integral2"      , ## (2D) numerical integration (as function, using scipy if possible)
    "integral3"      , ## (3D) numerical integration (as function, using scipy if possible)
    ##
    "Integral"       , ## (1D) numerical integration (as as object, using scipy if possible)
    "Integral2"      , ## (2D) numerical integration (as as object, using scipy if possible)
    "Integral3"      , ## (3D) numerical integration (as as object, using scipy if possible)
    ##
    'Integrate2D_X'  , ## partial integrationn of 2D-function over x-range  
    'Integrate2D_Y'  , ## partial integrationn of 2D-function over y-range
    ##
    'Integrate3D_X'  , ## partial integrationn of 2D-function over x-range  
    'Integrate3D_Y'  , ## partial integrationn of 2D-function over y-range
    'Integrate3D_Z'  , ## partial integrationn of 2D-function over z-range
    ##
    "romberg"        , ## (1D) numerical integration using romberg's method 
    'genzmalik2'     , ## (2D) numerical integration using Genz&Malik's method
    'genzmalik3'     , ## (3D) numerical integration using Genz&Malik's method
    ##
    "IntegralCache"  , ## (1D) numerical integration (as object, using scipy is if ssible)
    ##    
    ) 
# =============================================================================
import ROOT, warnings, math 
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
#  func = lambda x : x*x 
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

    >>> func = lambda x : x*x 
    >>> v    = romberg ( func , 0 , 1 ) 
    """
    # =========================================================================
    ## internal recursive function
    #  Romberg's adaptive integration with Richardson's extrapolation
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
    if  ( ea <= 0 or iszero ( ea ) ) and ( er <= 0 or iszero  ( er )  ) : 
        raise AttributeError("Romberg: absolute and relative tolerances can't be non-positive simultaneously")

    ## adjust maximal number of steps in Richardson's extrapolation 
    nmax        = max ( 2 , nmax )

    ## ensure that function returns double values: 
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
    def integral ( fun , xmin , xmax ,
                   args  = ()        , 
                   err   = False     ,
                   **kwargs          ) :
        """Calculate the integral for the 1D-function using scipy
        
        >>> func = lambda x : x * x 
        >>> v = integral(func,0,1)
        """
        func   = lambda x : float ( fun ( x , *args ) ) 
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = integrate.quad ( func , xmin , xmax , **kwargs )
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
    """Calculate the integral for the 1D-function
    
    >>> func  = lambda x : x * x      ## define function 
    >>> iint  = Integral ( func , 0 ) ## specify x_low 
    >>> value = iint (  10  )         ## specify x_high 
    """
    ## Calculate the integral for the 1D-function
    def __init__ ( self , func , xlow = 0 , args = () , err = False , **kwargs ) :
        """Calculate the integral for the 1D-function
        
        >>> func   = ...
        >>> func_0 = Integral(func,0)
        >>> value  = func_
        """
        self.__func   = func 
        self.__xmin   = float ( xlow ) 
        self.___err   = err
        self.__args   = args
        self.__kwargs = kwargs
        
    ## Calculate the integral for the 1D-function
    def _integrate_ ( self , xmn , xmx , args = () ) :
        args   = args if args else self.args
        return integral ( self.func   , xmn , xmx ,
                          args = args , err = self.err , **self.kwargs )
            
    ## Calculate the integral for the 1D-function 
    def __call__ ( self , x , *args ) :
        """Calculate the integral for the 1D-function 
    
        >>> func = ...
        >>> func_0 = Integral(func,0)
        >>> func_0 ( 10 ) 
        """
        return self._integrate_ ( self.xmin , x , args = args )

    @property
    def xmin ( self  ) :
        """Low-limit of inntegrtaion
        """
        return self.__xmin 

    @property
    def err ( self ) :
        """Flag to evaluate the integration uncnertainnties"""
        return self.__err

    @property
    def func ( self ) :
        """The integrand"""
        return self.__func
    
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
    """Calculate the integral for the 1D-function
    >>> func   = lambda x : x*x 
    >>> iint   = IntegralCache ( func , 0 ) ## specify x_low 
    >>> value  = iint ( 10 )                ## specify x_high
    """
    ## Calculate the integral for the 1D-function using scipy
    def __init__ ( self , func , xlow = 0 , args = () , err = False , **kwargs ) :
        """Calculate the integral for the 1D-function
        
        >>> func = ...
        >>> func_0 = Integral(func,0)
        """
        Integral.__init__ ( self , func , xlow , args , err , **kwargs )
        self._prev   = None 
        
    ## Calculate the integral for the 1D-function using scipy
    def __call__ ( self , x , **args ) :
        """Calculate the integral for the 1D-function using scipy
        
        >>> func = ...
        >>> func_0 = Integral(func,0)
        >>> func_0 ( 10 ) 
        """
        x      = float ( x ) 
        args   = args if args else self.args
        
        xmn    = self.xmin
        delta  = 0
        
        _x     = float ( x )
        
        ## there is previos calculation
        if self._prev :
            ##
            prev_args   = self._prev[0]
            prev_x      = self._prev[1]
            prev_result = self._prev[2]
            ##
            if prev_args == args : ## the same additional arguments                
                
                ## the point is good! 
                if prev_x == _x or isequal ( _x , prev_x ) :
                    return prev_result                                 ## RETURN
                
                ## old point is good 
                if abs ( prev_x - _x ) <= abs ( self.xmin - _x ) :
                    xmn   = prev_x
                    delta = prev_result
                    
        result  = self._integrate_ ( xmn , _x , args =   args )
        result += delta 
        
        ## fill cache 
        self._prev = args , _x , result 
            
        return result 

# =============================================================================
## 2D&3D integration 
# =============================================================================
_w  = ( lambda N : (12824-9120*N+400*N*N)/19683. ,
        lambda N : (980                 )/ 6561. ,
        lambda N : (1820-400*N          )/19683. ,
        lambda N : (200                 )/19683. ,
        lambda N : (6859                )/19683. ) 
_w2 = tuple ( [  w(2) for w in _w ] )
_w3 = tuple ( [  w(3) for w in _w ] )
_wp = ( lambda N : (729-950*N+50*N*N  )/ 729. ,
        lambda N : 245                 / 486. ,
        lambda N : (265-100*N         )/1458. ,
        lambda N : 25                  / 729. )
_w2p = tuple ( [  w(2) for w in _wp ] )
_w3p = tuple ( [  w(3) for w in _wp ] )
_l2  = math.sqrt ( 9/70. )
_l3  = math.sqrt ( 9/10. )
_l4  = math.sqrt ( 9/10. )
_l5  = math.sqrt ( 9/19. )
del _w,_wp
# =============================================================================
## Genz&Malik's basic rule for N=2 cubature
#  A.C. Genz, A.A. Malik, ``Remarks on algorithm 006: An adaptive algorithm for
#  numerical integration over an N-dimensional rectangular region'',
#  in Journal of Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295,
#  ISSN 0377-0427
#  @see https://doi.org/10.1016/0771-050X(80)90039-X.
#  @see http://www.sciencedirect.com/science/article/pii/0771050X8090039X)
#  @code
#  func  = lambda x,y : x*x + y * y
#  i7,i5 = _genzmalik2_ ( func , (-1,1) , (-1,1) )
#  print 'Integral %s, error: %s' % ( i7,  abs(i7-i5) )  
#  @endcode 
def _genzmalik2_ ( func , xlims , ylims , args = () ) :
    """ Genz&Malik's basic rule for N=2 cubature
    
    A.C. Genz, A.A. Malik, ``Remarks on algorithm 006: An adaptive algorithm for
    numerical integration over an N-dimensional rectangular region'',
    in Journal of Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295,
    ISSN 0377-0427
    - see https://doi.org/10.1016/0771-050X(80)90039-X.
    - see http://www.sciencedirect.com/science/article/pii/0771050X8090039X)
    
    Example:
    ========
    
    >>> func  = lambda x,y : x*x + y * y
    >>> i7,i5 = _genzmalik2_ ( func , (-1,1) , (-1,1) )
    >>> print 'Integral %s, error: %s' % ( i7,  abs(i7-i5) )  
    
    """
    xc = 0.5 * ( xlims[1] + xlims[0] )
    dx = 0.5 * ( xlims[1] - xlims[0] )
    yc = 0.5 * ( ylims[1] + ylims[0] )
    dy = 0.5 * ( ylims[1] - ylims[0] )

    _func = lambda x,y :   float ( func ( x , y , *args ) )
    
    s1  = _func ( xc , yc ) 
    
    s2  = _func ( xc + _l2 * dx , yc            ) + _func ( xc - _l2 * dx , yc              )
    s2 += _func ( xc            , yc + _l2 * dy ) + _func ( xc            , yc - _l2 * dy )
    
    s3  = _func ( xc + _l3 * dx , yc            ) + _func ( xc - _l3 * dx , yc              )
    s3 += _func ( xc            , yc + _l3 * dy ) + _func ( xc            , yc - _l3 * dy )
    
    s4  = _func ( xc + _l4 * dx , yc + _l4 * dy ) + _func ( xc + _l4 * dx , yc - _l4 * dy )
    s4 += _func ( xc - _l4 * dx , yc + _l4 * dy ) + _func ( xc - _l4 * dx , yc - _l4 * dy )

    s5  = _func ( xc + _l5 * dx , yc + _l5 * dy ) + _func ( xc + _l5 * dx , yc - _l5 * dy )
    s5 += _func ( xc - _l5 * dx , yc + _l5 * dy ) + _func ( xc - _l5 * dx , yc - _l5 * dy )

    i7  = _w2 [0] * s1 + _w2 [1] * s2 + _w2 [2] * s3 + _w2 [3] * s4 + _w2 [4] * s5/2**2 
    i5  = _w2p[0] * s1 + _w2p[1] * s2 + _w2p[2] * s3 + _w2p[3] * s4
    
    vol = 4 * dx * dy

    return i7*vol,i5*vol

# =============================================================================
## Genz&Malik's basic rule for N=3 cubature
#  A.C. Genz, A.A. Malik, ``Remarks on algorithm 006: An adaptive algorithm for
#  numerical integration over an N-dimensional rectangular region'',
#  in Journal of Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295,
#  ISSN 0377-0427
#  @see https://doi.org/10.1016/0771-050X(80)90039-X.
#  @see http://www.sciencedirect.com/science/article/pii/0771050X8090039X)
#  @code
#  func  = lambda x,y : x*x + y*y + z*z
#  i7,i5 = _genzmalik2_ ( func , (-1,1) , (-1,1) , (-1,1) )
#  print 'Integral %s, error: %s' % ( i7,  abs(i7-i5) )  
#  @endcode 
def _genzmalik3_ ( func , xlims , ylims , zlims , args =  () ) :
    """ Genz&Malik's basic rule for N=3 cubature
    
    A.C. Genz, A.A. Malik, ``Remarks on algorithm 006: An adaptive algorithm for
    numerical integration over an N-dimensional rectangular region'',
    in Journal of Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295,
    ISSN 0377-0427
    - see https://doi.org/10.1016/0771-050X(80)90039-X.
    - see http://www.sciencedirect.com/science/article/pii/0771050X8090039X)
    
    Example:
    ========
    
    >>> func  = lambda x,y,z : x*x + y*y + z*z
    >>> i7,i5 = _genzmalik2_ ( func , (-1,1) , (-1,1) , (-1,1) )
    >>> print 'Integral %s, error: %s' % ( i7,  abs(i7-i5) )  
    
    """
    xc = 0.5 * ( xlims[1] + xlims[0] )
    dx = 0.5 * ( xlims[1] - xlims[0] )
    yc = 0.5 * ( ylims[1] + ylims[0] )
    dy = 0.5 * ( ylims[1] - ylims[0] )
    zc = 0.5 * ( zlims[1] + zlims[0] )
    dz = 0.5 * ( zlims[1] - zlims[0] )

    _f = lambda x,y,z : float ( func ( x , y , z , *args ) )
    
    s1  = _f ( xc , yc , zc ) 
    
    s2  = _f ( xc + _l2 * dx , yc            , zc            )
    s2 += _f ( xc - _l2 * dx , yc            , zc            )
    s2 += _f ( xc            , yc + _l2 * dy , zc            )
    s2 += _f ( xc            , yc - _l2 * dy , zc            )
    s2 += _f ( xc            , yc            , zc + _l2 * dz )
    s2 += _f ( xc            , yc            , zc - _l2 * dz )

    s3  = _f ( xc + _l3 * dx , yc            , zc            )
    s3 += _f ( xc - _l3 * dx , yc            , zc            )
    s3 += _f ( xc            , yc + _l3 * dy , zc            )
    s3 += _f ( xc            , yc - _l3 * dy , zc            )
    s3 += _f ( xc            , yc            , zc + _l3 * dz )
    s3 += _f ( xc            , yc            , zc - _l3 * dz )

    s4  = _f ( xc + _l4 * dx , yc + _l4 * dy , zc            )
    s4 += _f ( xc + _l4 * dx , yc            , zc + _l4 * dz )
    s4 += _f ( xc            , yc + _l4 * dy , zc + _l4 * dz )
    s4 += _f ( xc + _l4 * dx , yc - _l4 * dy , zc            )
    s4 += _f ( xc + _l4 * dx , yc            , zc - _l4 * dz )
    s4 += _f ( xc            , yc + _l4 * dy , zc - _l4 * dz )
    s4 += _f ( xc - _l4 * dx , yc + _l4 * dy , zc            )
    s4 += _f ( xc - _l4 * dx , yc            , zc + _l4 * dz )
    s4 += _f ( xc            , yc - _l4 * dy , zc + _l4 * dz )
    s4 += _f ( xc - _l4 * dx , yc - _l4 * dy , zc            )
    s4 += _f ( xc - _l4 * dx , yc            , zc - _l4 * dz )
    s4 += _f ( xc            , yc - _l4 * dy , zc - _l4 * dz )

    s5  = _f ( xc + _l5 * dx , yc + _l5 * dy , zc + _l5 * dz )
    s5 += _f ( xc + _l5 * dx , yc + _l5 * dy , zc - _l5 * dz )
    s5 += _f ( xc + _l5 * dx , yc - _l5 * dy , zc + _l5 * dz )
    s5 += _f ( xc + _l5 * dx , yc - _l5 * dy , zc - _l5 * dz )
    s5 += _f ( xc - _l5 * dx , yc + _l5 * dy , zc + _l5 * dz )
    s5 += _f ( xc - _l5 * dx , yc + _l5 * dy , zc - _l5 * dz )
    s5 += _f ( xc - _l5 * dx , yc - _l5 * dy , zc + _l5 * dz )
    s5 += _f ( xc - _l5 * dx , yc - _l5 * dy , zc - _l5 * dz )
    
    i7  = _w3 [0] * s1 + _w3 [1] * s2 + _w3 [2] * s3 + _w3 [3] * s4 + _w3 [4] * s5/2**3 
    i5  = _w3p[0] * s1 + _w3p[1] * s2 + _w3p[2] * s3 + _w3p[3] * s4
    
    vol = 8 * dx * dy * dz 

    return i7*vol,i5*vol

# ============================================================================
## split 2D-region into  four smaller pieces
#  @code
#  region     =  (-1,1),( -2,5)
#  newregions = _split2_( region )
#  for nr in  newregions : print nr
#  @endcode 
def _split2_ ( xlims , ylims ) :
    """Split 2D-region into   four smaller pieces

    Example:
    ========
    
    >>> region     =  (-1,1),( -2,5)
    >>> newregions = _split2_( region )
    >>> for nr in  newregions : print nr    
    """
    xc = 0.5 * ( xlims[1] + xlims[0] )
    dx = 0.5 * ( xlims[1] - xlims[0] )
    yc = 0.5 * ( ylims[1] + ylims[0] )
    dy = 0.5 * ( ylims[1] - ylims[0] )
    
    return ( ( ( xc-dx , xc    ) , ( yc-dy , yc    ) ) ,
             ( ( xc    , xc+dx ) , ( yc-dy , yc    ) ) ,
             ( ( xc-dx , xc    ) , ( yc    , yc+dy ) ) ,
             ( ( xc    , xc+dx ) , ( yc    , yc+dy ) ) )
    

# ============================================================================
## split 3D-region into  eight smaller pieces
#  @code
#  region     =  (-1,1),( -2,5) , (-3,3)
#  newregions = _split3_( region )
#  for nr in  newregions : print nr
#  @endcode 
def _split3_ ( xlims , ylims , zlims ) :
    """Split 3D-region into eight smaller pieces

    Example:
    ========
    
    >>> region     =  (-1,1),( -2,5) , ( -3,3) 
    >>> newregions = _split3_( region )
    >>> for nr in  newregions : print nr    
    """
    xc = 0.5 * ( xlims[1] + xlims[0] )
    dx = 0.5 * ( xlims[1] - xlims[0] )
    yc = 0.5 * ( ylims[1] + ylims[0] )
    dy = 0.5 * ( ylims[1] - ylims[0] )
    zc = 0.5 * ( zlims[1] + zlims[0] )
    dz = 0.5 * ( zlims[1] - zlims[0] )
    
    return ( ( ( xc-dx , xc    ) , ( yc-dy , yc    ) , ( zc-dz , zc      ) ) ,
             ( ( xc-dx , xc    ) , ( yc-dy , yc    ) , ( zc    , zc + dz ) ) ,
             ( ( xc-dx , xc    ) , ( yc    , yc+dy ) , ( zc-dz , zc      ) ) ,
             ( ( xc-dx , xc    ) , ( yc    , yc+dy ) , ( zc    , zc + dz ) ) ,
             ( ( xc    , xc+dx ) , ( yc-dy , yc    ) , ( zc-dz , zc      ) ) ,
             ( ( xc    , xc+dx ) , ( yc-dy , yc    ) , ( zc    , zc + dz ) ) ,
             ( ( xc    , xc+dx ) , ( yc    , yc+dy ) , ( zc-dz , zc      ) ) ,
             ( ( xc    , xc+dx ) , ( yc    , yc+dy ) , ( zc    , zc + dz ) ) )
             
# =============================================================================
## Driving routine for Adaptive numerical 2D/3D integration using Genz&Malik's basic rule
# 
#  A.C. Genz, A.A. Malik, ``Remarks on algorithm 006: An adaptive algorithm for
#  numerical integration over an N-dimensional rectangular region'',
#  in Journal of Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295,
#   ISSN 0377-0427
#  @see https://doi.org/10.1016/0771-050X(80)90039-X.
#  @see http://www.sciencedirect.com/science/article/pii/0771050X8090039X)
def _genzmalik_( func , limits , basic_rule , splitter ,
                 args = () ,  epsabs = 1.5e-7 , epsrel = 1.5e-7 ) :
    """Driving routine for Adaptive numerical 2D/3D integration using Genz&Malik's basic rule
    
    A.C. Genz, A.A. Malik, ``Remarks on algorithm 006: An adaptive algorithm for
    numerical integration over an N-dimensional rectangular region'',
    in Journal of Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295,
    ISSN 0377-0427
    - see https://doi.org/10.1016/0771-050X(80)90039-X.
    - see http://www.sciencedirect.com/science/article/pii/0771050X8090039X)
    """

    i7,i5 = basic_rule ( func , *limits , args = args )
    
    nfc   = 17
    
    err     =  abs ( i7 -  i5 )
    relerr  =  0.5 *  ( abs(i7) + abs(i5) ) * epsrel 
    
    if   err <= epsabs : return i7, err, nfc , 1 
    elif err <= relerr : return i7, err, nfc , 1 
    
    stack = {}
    r     = limits 
    
    stack [ r ] = err, i7

    while 1 < 2 :
        
        serr =  0    ## (current) sum of errors
        res  =  0    ## (current) result 
        rmx  =  None ## (current) region with the maximal error 
        emx  = -1    ## (current) maximal error

        for r,entry in stack.iteritems() :

            err     = entry[0]            
            serr   += err
            res    += entry[1]

            if 0 > emx or err >= emx :
                emx = err
                rmx = r 

            if   serr <= 0.5 * epsabs : pass
            elif serr <= 0.5 * relerr : pass
            else                      : break ##  BREAK

        else :
            # cumulated uncertainty is small enough, return  
            lstack = len(stack)
            del stack 
            return res, serr, nfc , lstack

        ##  ??? 
        if not rmx : break     ## BREAK 
        
        rr = splitter ( *rmx )
        for r in  rr :
            r7,r5  = basic_rule (  func , *r , args = args )
            nfc   += 17 
            stack[ r ] = abs(r7-r5),r7
        ## remove large region from the stack 
        del stack[rmx]            

    lstack = len(stack)
    del stack 
    return res, serr, nfc , lstack
        

# =============================================================================
## Adaptive numerical 2D integration using Genz&Malik's basic rule
# 
#  A.C. Genz, A.A. Malik, ``Remarks on algorithm 006: An adaptive algorithm for
#  numerical integration over an N-dimensional rectangular region'',
#  in Journal of Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295,
#   ISSN 0377-0427
#  @see https://doi.org/10.1016/0771-050X(80)90039-X.
#  @see http://www.sciencedirect.com/science/article/pii/0771050X8090039X)
#  @code
#  func   = lambda x,y : x*x + y*y
#  r      = genzmalik2 ( func , xmin=-1 , xmax=2 , ymin=-1 , ymax=2 )
#  print 'Integral: %s ' % r 
#  @endcode 
def genzmalik2 ( func , xmin , xmax , ymin ,  ymax , args = () , err = False , epsabs = 1.5e-7 , epsrel = 1.5e-7 ) :
    """ Adaptive numerical 2D integration using Genz&Malik's basic rule
    
    A.C. Genz, A.A. Malik, ``Remarks on algorithm 006: An adaptive algorithm for
    numerical integration over an N-dimensional rectangular region'',
    in Journal of Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295,
    ISSN 0377-0427
    - see https://doi.org/10.1016/0771-050X(80)90039-X.
    - see http://www.sciencedirect.com/science/article/pii/0771050X8090039X)
    
    Example:
    ========
    
    >>> func   = lambda x,y : x*x + y*y
    >>> r      = genzmalik2 ( func , xmin=-1 , xmax=2 , ymin=-1 , ymax=2 )
    >>> print 'Integral: %s ' % r 
    
    """

    limits  = ( xmin , xmax ) , ( ymin , ymax ) 
    r,e,n,s = _genzmalik_ ( func , limits , _genzmalik2_ , _split2_ , args , abs(epsabs) , abs(epsrel) )
    
    return VE ( r , e * e )  if err else r 

# =============================================================================
## Adaptive numerical 3D integration using Genz&Malik's basic rule
# 
#  A.C. Genz, A.A. Malik, ``Remarks on algorithm 006: An adaptive algorithm for
#  numerical integration over an N-dimensional rectangular region'',
#  in Journal of Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295,
#   ISSN 0377-0427
#  @see https://doi.org/10.1016/0771-050X(80)90039-X.
#  @see http://www.sciencedirect.com/science/article/pii/0771050X8090039X)
#  @code
#  func   = lambda x,y : x*x + y*y + z*z 
#  r      = genzmalik3 ( func , xmin=-1 , xmax=2 , ymin=-1 , ymax=2 , zmin = -4, zmax = 7)
#  print 'Integral: %s ' % r 
#  @endcode 
def genzmalik3 ( func , xmin , xmax , ymin ,  ymax , zmin ,  zmax , args = () , err = False , epsabs = 1.5e-7 , epsrel = 1.5e-7 ) :
    """ Adaptive numerical 3D integration using Genz&Malik's basic rule
    
    A.C. Genz, A.A. Malik, ``Remarks on algorithm 006: An adaptive algorithm for
    numerical integration over an N-dimensional rectangular region'',
    in Journal of Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295,
    ISSN 0377-0427
    - see https://doi.org/10.1016/0771-050X(80)90039-X.
    - see http://www.sciencedirect.com/science/article/pii/0771050X8090039X)
    
    Example:
    ========
    
    >>> func   = lambda x,y : x*x + y*y + z*z 
    >>> r      = genzmalik3 ( func , xmin=-1 , xmax=2 , ymin=-1 , ymax=2 , zmin = -4, zmax = 7)
    >>> print 'Integral: %s ' % r 
    
    """

    limits  = ( xmin , xmax ) , ( ymin , ymax ) , ( zmin , zmax ) 
    r,e,n,s = _genzmalik_ ( func , limits , _genzmalik3_ , _split3_ , args , abs(epsabs) , abs(epsrel) )
    
    return VE ( r , e * e )  if err else r 



# =============================================================================
try :
    from scipy import integrate
    # =============================================================================
    ## Calculate the integral (from ) for the 2D-function 
    #  @code 
    #  func = lambda x,y : x*x + y*y
    #  v    = integral2 ( func , 0 , 1 ,  -2 , 2 )
    #  @endcode 
    #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    #  @date   2014-06-06
    def integral2 ( fun  ,
                    xmin , xmax   ,
                    ymin , ymax   ,
                    args  =  ()   ,
                    err   = False ,
                    **kwargs      ) :
        """Calculate the integral for the 2D-function using scipy
        
        >>> func = lambda x,y : x*x+y*y 
        >>> v = integral2(func,0,1,-2,2)
        """
        func   = lambda x,y : float ( fun ( x , y , *args ) ) 
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = integrate.dblquad ( func ,
                                         ymin , ymax     ,
                                         lambda x : xmin ,
                                         lambda x : xmax , **kwargs )
            return VE( result[0] , result[1] * result[1] ) if err else result[0]
        
except ImportError :
    logger.warning ("scipy.integrate is not available, use local ``genz&malik''-replacement")
    ## use Genz&Malik integration as default method when scipy is not available 
    integral2 = genzmalik2 
# =============================================================================



# =============================================================================
try :
    from scipy import integrate
    # =============================================================================
    ## Calculate the inteegral for the 3D-function 
    #  @code 
    #  func = lambda x,y,z: x*x+y*y+z*z
    #  v    = integral3 ( func , 0 , 1 ,  0, 2 , 0 , 3 )
    #  @endcode 
    #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    #  @date   2014-06-06
    def integral3 ( fun  ,
                    xmin , xmax   ,
                    ymin , ymax   ,
                    zmin , zmax   ,
                    args  =  ()   ,
                    err   = False ,
                    **kwargs      ) :
        """Calculate the integral for the 3D-function using scipy
        
        >>> func = lambda x,y,z : x*x+y*y+z*z
        >>> v = integral3(func,0,1,0,2,0,3)
        """
        func   = lambda x,y,z : float ( fun ( x , y , z , *args ) ) 
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = integrate.tplquad ( func ,
                                         zmin , zmax ,
                                         lambda z   : ymin ,
                                         lambda z   : ymax ,
                                         lambda y,z : xmin ,
                                         lambda y,z : xmax , **kwargs )
            return VE( result[0] , result[1] * result[1] ) if err else result[0]
        
except ImportError :
    logger.warning ("scipy.integrate is not available, use local ``genz&malik''-replacement")
    ## use Genz&Malik integration as default method when scipy is not available 
    integral3 = genzmalik3 
# =============================================================================


# =============================================================================
## @class Integrate2D_Y
# helper class to perform (partial) integration of 2D function
# \f$  f(x) = \int_{y_{min}}^{y_{max}} f_{2D}(x,y) dy \f$ 
# @code
# fun2d  = ... ## 2D-function
# fx     = Integral2D_Y( fun2d , ymin = 0 , ymax = 1 )
# print fx ( 1 ) 
# @endcode 
class Integrate2D_Y(Integral) :
    """Helper class to perform (partial) integration of 2D function
    
    f(x) = \int_{y_{min}}^{y_{max}} f_{2D}(x,y) dy
    
    >>> fun2d  = ... ## 2D-function
    >>> fx     = Integral2D_Y( fun2d , ymin = 0 , ymax = 1 )
    >>> print fx ( 1 )
    """

    ## construct the integration object 
    def __init__  ( self , fun2d , ymin , ymax , args = () , err = False , **kwargs ) :
        """Construct the integration object
        
        f(x) = \int_{y_{min}}^{y_{max}} f_{2D}(x,y) dy
        
        >>> fun2d  = ... ## 2D-function
        >>> fx     = Integral2D_Y( fun2d , ymin = 0 , ymax = 1 )
        >>> print fx ( 1 )
        """
        self.__ymax = ymax 
        Integral.__init__ ( self , fun2d , ymin , args , err , **kwargs )

    ##  evaluat the funcntion (perform y-integration) 
    def __call__ ( self , x  , *args ) :
        """Evaluat the function (perform y-integration)
        >>> fun2d  = ... ## 2D-function
        >>> fx     = Integral2D_Y( fun2d , ymin = 0 , ymax = 1 )
        >>> print fx ( 1 )        
        """
        
        args = args if  args else self.args

        ## create the helper function 
        _funy_ = lambda y , *_ : self.func ( x , y , *_ )
        
        return intergal ( _funy_ , self.ymin , self.ymax , 
                          args = args  , err  = self.err , **self.kwargs )
    
    @property 
    def ymin ( self ) :
        """Low integration limit (the same as xmin)"""
        return self.xmin
    
    @property 
    def ymax ( self ) :
        """High integration limit"""
        return self.__ymax


# =============================================================================
## @class Integrate2D_X
# helper class to perform (partial) integration of 2D function
# \f$  f(y) = \int_{x_{min}}^{x_{max}} f_{2D}(x,y) dx \f$ 
# @code
# fun2d  = ... ## 2D-function
# fy     = Integral2D_X( fun2d , xmin = 0 , xmax = 1 )
# print fy ( 1 ) 
# @endcode 
class Integrate2D_X(Integral) :
    """Helper class to perform (partial) integration of 2D function
    
    f(y) = \int_{x_{min}}^{x_{max}} f_{2D}(x,y) dx
    
    >>> fun2d  = ... ## 2D-function
    >>> fy     = Integral2D_X( fun2d , xmin = 0 , xmax = 1 )
    >>> print fy ( 1 )
    """

    ## construct the integration object 
    def __init__  ( self , fun2d , xmin , xmax , args = () , err = False , **kwargs ) :
        """Construct the integration object
        
        f(y) = \int_{x_{min}}^{x_{max}} f_{2D}(x,y) dx
        
        >>> fun2d  = ... ## 2D-function
        >>> fy     = Integral2D_X( fun2d , xmin = 0 , xmax = 1 )
        >>> print fy ( 1 )
        """
        self.__xmax = xmax 
        Integral.__init__ ( self , fun2d , xmin , args , err , **kwargs )

    ##  evaluat the funcntion (perform y-integration) 
    def __call__ ( self , y  , *args ) :
        """Evaluat the function (perform y-integration)
        >>> fun2d  = ... ## 2D-function
        >>> fy     = Integral2D_X( fun2d , xmin = 0 , xmax = 1 )
        >>> print fy ( 1 )        
        """
        
        args = args if  args else self.args

        ## create the helper function 
        _funx_ = lambda x , *_ : self.func ( x , y , *_ )
        
        return intergal ( _funx_ , self.xmin , self.xmax , 
                          args = args  , err  = self.err , **self.kwargs )
    
    @property 
    def xmax ( self ) :
        """High integration limit"""
        return self.__xmax


# =============================================================================
## @class Integrate3D_Y
# helper class to perform (partial) integration of 3D function
# \f$  f(x,z) = \int_{y_{min}}^{y_{max}} f_{2D}(x,y,z) dy \f$ 
# @code
# fun3d  = ... ## 2D-function
# fxy     = Integral2D_Y( fun3d , ymin = 0 , ymax = 1 )
# print fxy ( 1 , 2 ) 
# @endcode 
class Integrate3D_Y(Integrate2D_Y) :
    """Helper class to perform (partial) integration of 3D function
    
    f(x,z) = \int_{y_{min}}^{y_{max}} f_{3D}(x,y,z) dy
    
    >>> fun3d  = ... ## 3D-function
    >>> fxy     = Integral3D_Y( fun3d , ymin = 0 , ymax = 1 )
    >>> print fxy ( 1 , 2 )
    """

    ## construct the integration object 
    def __init__  ( self , fun3d , ymin , ymax , args = () , err = False , **kwargs ) :
        """Construct the integration object
        
        f(x,z) = \int_{y_{min}}^{y_{max}} f_{3D}(x,y,z) dy
        
        >>> fun3d  = ... ## 3D-function
        >>> fxy     = Integral3D_Y( fun3d , ymin = 0 , ymax = 1 )
        >>> print fxy ( 1 , 2 )
        """
        Integrate2D_Y.__init__ ( self , fun3d , ymin , ymax , args , err , **kwargs )
        
    ##  evaluat the funcntion (perform y-integration) 
    def __call__ ( self , x  , z , *args ) :
        """Evaluate the function (perform y-integration)
        
        f(x,z) = \int_{y_{min}}^{y_{max}} f_{3D}(x,y,z) dy
        
        >>> fun3d  = ... ## 3D-function
        >>> fxy     = Integral3D_Y( fun3d , ymin = 0 , ymax = 1 )
        >>> print fxy ( 1 , 2 )
        """
        args = args if  args else self.args

        ## create the helper function 
        _funy_ = lambda y , *_ : self.func ( x , y , z , *_ )
        
        return intergal ( _funy_ , self.ymin , self.ymax , 
                          args = args  , err  = self.err , **self.kwargs )


# =============================================================================
## @class Integrate3D_X
# helper class to perform (partial) integration of 3D function
# \f$  f(y,z) = \int_{x_{min}}^{x_{max}} f_{2D}(x,y,z) dx \f$ 
# @code
# fun3d  = ... ## 2D-function
# fyz     = Integral2D_X( fun3d , xmin = 0 , xmax = 1 )
# print fyz ( 1 , 2 ) 
# @endcode 
class Integrate3D_X(Integrate2D_X) :
    """Helper class to perform (partial) integration of 3D function
    
    f(y,z) = \int_{x_{min}}^{x_{max}} f_{3D}(x,y,z) dx
    
    >>> fun3d  = ... ## 3D-function
    >>> fyz     = Integral3D_X( fun3d , xmin = 0 , xmax = 1 )
    >>> print fyz ( 1 , 2 )
    """

    ## construct the integration object 
    def __init__  ( self , fun3d , xmin , xmax , args = () , err = False , **kwargs ) :
        """Construct the integration object
        
        f(y,z) = \int_{x_{min}}^{x_{max}} f_{3D}(x,y,z) dx
        
        >>> fun3d  = ... ## 3D-function
        >>> fyz     = Integral3D_X( fun3d , xmin = 0 , xmax = 1 )
        >>> print fyz ( 1 , 2 )
        """
        Integrate2D_X.__init__ ( self , fun3d , xmin , xmax , args , err , **kwargs )
        
    ##  evaluat the funcntion (perform y-integration) 
    def __call__ ( self , y  , z , *args ) :
        """Evaluate the function (perform y-integration)
        
        f(y,z) = \int_{x_{min}}^{x_{max}} f_{3D}(x,y,z) dx
        
        >>> fun3d  = ... ## 3D-function
        >>> fyz     = Integral3D_X( fun3d , xmin = 0 , xmax = 1 )
        >>> print fyz ( 1 , 2 )
        """
        args = args if  args else self.args

        ## create the helper function 
        _funx_ = lambda x , *_ : self.func ( x , y , z , *_ )
        
        return intergal ( _funx_ , self.xmin , self.xmax , 
                          args = args  , err  = self.err , **self.kwargs )
    

# =============================================================================
## @class Integrate3D_Z
# helper class to perform (partial) integration of 3D function
# \f$  f(x,y) = \int_{z_{min}}^{z_{max}} f_{2D}(x,y,z) dz \f$ 
# @code
# fun3d  = ... ## 2D-function
# fxy     = Integral2D_Z( fun3d , zmin = 0 , zmax = 1 )
# print fxy ( 1 , 2 ) 
# @endcode 
class Integrate3D_Z(Integrate2D_X) :
    """Helper class to perform (partial) integration of 3D function
    
    f(x,y) = \int_{z_{min}}^{z_{max}} f_{3D}(x,y,z) dz
    
    >>> fun3d  = ... ## 3D-function
    >>> fxy     = Integral3D_Z( fun3d , zmin = 0 , zmax = 1 )
    >>> print fxy ( 1 , 2 )
    """

    ## construct the integration object 
    def __init__  ( self , fun3d , zmin , zmax , args = () , err = False , **kwargs ) :
        """Construct the integration object
        
        f(x,y) = \int_{z_{min}}^{z_{max}} f_{3D}(x,y,z) dz
        
        >>> fun3d  = ... ## 3D-function
        >>> fxy     = Integral3D_Z( fun3d , zmin = 0 , zmax = 1 )
        >>> print fxy ( 1 , 2 )
        """
        Integrate3D_X.__init__ ( self , fun3d , xmin , xmax , args , err , **kwargs )
        
    ##  evaluat the funcntion (perform y-integration) 
    def __call__ ( self , x  , y , *args ) :
        """Evaluate the function (perform z-integration)
        
        f(x,y) = \int_{z_{min}}^{z_{max}} f_{3D}(x,y,z) dz
        
        >>> fun3d  = ... ## 3D-function
        >>> fxy     = Integral3D_Z( fun3d , zmin = 0 , zmax = 1 )
        >>> print fxy ( 1 , 2 )        
        """
        args = args if  args else self.args

        ## create the helper function 
        _funx_ = lambda z , *_ : self.func ( x , y , z , *_ )
        
        return intergal ( _funx_ , self.zmin , self.zmax , 
                          args = args  , err  = self.err , **self.kwargs )
    
    @property 
    def zmin ( self ) :
        """Low integration limit"""
        return self.xmin
    @property 
    def zmax ( self ) :
        """High integration limit"""
        return self.xmax


    
# =============================================================================
## 2D-integration
# =============================================================================

# =============================================================================
## @class Integral2
#  Calculate the integral (from x0 to x and y0 to  y ) for the 2D-function 
#  @code 
#  func  = lambda x,y : x*x+y*y 
#  iint  = Integral2 ( func , 0 , 0 ) ## specify x_low, y_low 
#  value = iiint (  10  , 20 )        ## specify x_high, y_high  
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class Integral2(Integral) :
    """Calculate the integral for the 2D-function
    
    >>> func  = lambda x,y : x*x+y*y  ## define function 
    >>> iint  = Integral ( func , 0 , 0 ) ## specify x_low,  y_low  
    >>> value = iint (  10  , 20 )        ## specify x_high, y_high
    """
    ## Calculate the integral for the 2D-function
    def __init__ ( self , func ,
                   xlow = 0  ,
                   ylow = 0  ,
                   args = () , err = False , **kwargs ) :
        """Calculate the integral for the 2D-function
        
        >>> func   = ...
        >>> func_0 = Integral(func,0,0)
        >>> value  = func_0(10,20) 
        """
        Integral.__init__ (   self , func , xlow , args , err , **kwargs )
        self._ymin   = float ( ylow ) 
        
    ## Calculate the integral for the 2D-function 
    def _integrate_ ( self , xmn , xmx , ymn , ymx , args = () ) :
        args = args if  args else self._args 
        return integral2 ( self.func ,
                           xmn  , xmx ,
                           ymn  , ymx ,
                           args       , 
                           self.err  , **self.kwargs )
    
    ## Calculate the integral for the 2D-function 
    def __call__ ( self , x , y , *args ) :
        """Calculate the integral for the 2D-function
        
        >>> func = ...
        >>> func_0 = Integral(func,0,0)
        >>> func_0 ( 10 ,  20 ) 
        """
        return self._integrate_ ( self.xmin , x ,
                                  self._ymin , y , args = args )


# =============================================================================
## @class Integral3
#  Calculate the integral (from x0 to x, y0 to  y, z0  to   z ) for the 3D-function 
#  @code 
#  func  = lambda x,y : x*x+y*y+z*z 
#  iint  = Integral2 ( func , 0 , 0 , 0 ) ## specify x_low, y_low, z_low 
#  value = iiint (  10  , 20 , 30 )        ## specify x_high, y_high , z_high
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class Integral3(Integral2) :
    """Calculate the integral for the 3D-function
    
    >>> func  = lambda x,y : x*x+y*y+z*z  ## define function 
    >>> iint  = Integral ( func , 0 , 0 , 0 ) ## specify x_low,  y_low, z_low
    >>> value = iint (  10  , 20 , 30 )       ## specify x_high, y_high, z_high
    """
    ## Calculate the integral for the 3D-function 
    def __init__ ( self , func ,
                   xlow = 0  ,
                   ylow = 0  ,
                   zlow = 0  ,
                   args = () , err = False , **kwargs ) :
        """Calculate the integral for the 3D-function 
        
        >>> func   = ...
        >>> func_0 = Integral(func,0,0,0)
        >>> value  = func_0(10,20,30) 
        """
        Integral2.__init__ ( self , func , xlow , ylow , args , err  , **kwargs )
        self._zmin   = float ( zlow )
        
    ## Calculate the integral for the 3D-function 
    def _integrate_ ( self ,
                      xmn , xmx ,
                      ymn , ymx ,
                      zmn , zmx ,
                      args = () ) :
        args = args if  args else self._args 
        return integral3 ( self.func  ,
                           xmn , xmx  ,
                           ymn , ymx  ,
                           zmn , zmx  ,
                           args       , 
                           self.err  , **self.kwargs )
    
    ## Calculate the integral for the 3D-function 
    def __call__ ( self , x , y , z , *args ) :
        """Calculate the integral for the 3D-function
        
        >>> func = ...
        >>> func_0 = Integral(func,0,0,0)
        >>> func_0 ( 10 ,  20 ,   30 ) 
        """
        return self._integrate_ ( self.xmin , x ,
                                  self.ymin , y ,
                                  self.zmin , z , args = args )




# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
