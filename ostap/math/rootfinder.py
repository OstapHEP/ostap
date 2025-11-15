#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/rootfinder.py
#  Module with some useful utilities for root finding
#  - a kind of replacement for Brent's method when scipy is not accessible
#
#  The main entry point is a function <code>findroot</code>.
#  It finds a bracketed root of the given function for given interval.
#  - it is either a thin wrapper for <code>scipy.optimize.brentq</code> or
#  - a homemade replacement when <code>scipy</code> is not accesible
#
#  A few addtitional trivial  functions:
#
#  - <code>halley_newton</code>      : single step of Halley/Newton's   iteration
#  - <code>steffensen</code>         : single step of Steffensen's      iteration
#  - <code>inverse_parabolic</code>  : single step of inverse parabolic interpolation 
#  - <code>inverse_cubic</code>      : single step of inverse cubic     interpolation 
#  - <code>inverse_polynomial</code> : single step of inverse polynomial (1,2,3) interpolation 
#  - <code>secant</code>             : single step of secant/regular falsi/inverse linear
#  - <code>regular_falsi</code>      : ditto 
#  - <code>inverse_linear</code>     : ditto 
#  - <code>bisection</code>          : single step of bisection
#  - <code>aitken_delta2</code>      : aitken delta2 acceleration process
#
#  @see https://en.wikipedia.org/wiki/Halley%27s_method
#  @see https://en.wikipedia.org/wiki/Newton%27s_method
#  @see https://en.wikipedia.org/wiki/Steffensen%27s_method
#  @see https://en.wikipedia.org/wiki/Secant_method
#  @see https://en.wikipedia.org/wiki/False_position_method
#  @see https://en.wikipedia.org/wiki/Inverse_quadratic_interpolation
#  @see https://en.wikipedia.org/wiki/Aitken%27s_delta-squared_process
#  @see https://en.wikipedia.org/wiki/Bisection_method
#
#  @author Vanya Belyaev@itep.ru
#  @date   2018-06-15
# =============================================================================
""" Module with some useful utilities for root finding
- a kind of home-made replacement for Brent's method when scipy is not accessible

The main entry point is a function <code>findroot</code>.
It finds a bracketed root of the given function for given interval.
- it is either a thin wrapper for <code>scipy.optimize.brentq</code> or
- a homemade replacement when <code>scipy</code> is not accesible

A few addtitional trivial  functions:

- halley_newton      : single step of Halley/Newton's   iteration
- steffensen         : single step of Steffensen's      iteration
- inverse_parabolic  : single step of inverse parabolic interpolation 
- inverse_cubic      : single step of inverse cubic     interpolation 
- inverse_polynomial : single step of inverse polynomial (1,2,3) interpolation 
- secant             : single step of secant/regular falsi/inverse linear
- regular_falsi      : ditto 
- inverse_linear     : ditto
- bisection          : single step of bisection
- aitken_delta2      : aitken delta-2 acceleration process

"""
# =============================================================================
__version__ = "$Revision:$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = (
    ## 
    'findroot'           , ## find roots for given function in the interval
    'solve'              , ## solve  f(x)=C equation 
    ## homemade stuff 
    'find_root'          , ## local homemade rootfinder
    'RootFinder'         , ## the actual root-finder 
    'Point'              , ## helper class, namedtuple: (x,fun(x))
    ## trivial helper functions
    'halley_newton'      , ## single step of Halley/Newton's iteration 
    'steffensen'         , ## single step of Steffensen's    iteration 
    'inverse_parabolic'  , ## single step of inverse parabolic interpolation 
    'inverse_cubic'      , ## single step of inverse cubic     interpolation 
    'inverse_polynomial' , ## single step of inverse polynomial-(1,2,3) interpolation 
    'secant'             , ## single step of secant/regular falsi/inverse linear
    'bisection'          , ## single step of bisection
    'aitken_delta2'      , ## aitken delta2 acceleration process   
)
# =============================================================================
from   ostap.core.ostap_types import num_types  
from   ostap.math.base        import samesign, iszero, isequal, isfinite, signum    
from   ostap.utils.basic      import counted
from   ostap.stats.counters   import EffCounter
from   ostap.logger.pretty    import pretty_float, fmt_pretty_values 
from   ostap.logger.colorized import attention, allright  
import ostap.logger.table     as     T 
import math, sys, collections
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.rootfinder' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
_xtol = 2.e-12
_rtol = sys.float_info.epsilon * 6 
# =============================================================================
## helper class to keef (x,funxc)) information together 
Point = collections.namedtuple('Point', ( 'x','fx' ) )
# =============================================================================
## Single step of Halley/Newton's root-finding iteration process 
#  @see https://en.wikipedia.org/wiki/Halley%27s_method
#  @see https://en.wikipedia.org/wiki/Newton%27s_method
#  @code 
#  fun   = lambda x : sin ( x )
#  deriv = lambda x : cos ( x )
#  x = 0.5
#  for i in  range ( 10 ) :
#      x = halley_newton ( fun , x , deriv )
#      print i, x
#  @endcode 
def halley_newton ( fun            ,      ## the function 
                    x              ,      ## x 
                    deriv1         , * ,  ## the first derivative 
                    deriv2  = None ,      ## the second derivative
                    fx      = None ) :    ## (optional) value of fun(x)
    """ Single step of Newton/Halley's algorithms

    Parameters
    -----------
    
    fun    : the function 
    x      : the initial guess for the root positon
    deriv  : the first derivative  (for Newton and Halley's method) 
    deriv2 : the second derivative (for Halley's method)
    fx     : (optional) the value of function at the guess point

    Returns
    -------
    The next approximation to the root or `None`
    
    Example
    -------
    
    >>> fun   = lambda x : sin ( x )
    >>> deriv = lambda x : cos ( x )
    >>> x = 0.5
    >>> for i in  range ( 10 ) :
    ...      x = halley_newton ( fun , x , deriv )
    ...      print i, x

    References
    ----------
    - see https://en.wikipedia.org/wiki/Halley%27s_method
    - see https://en.wikipedia.org/wiki/Newton%27s_method     
   
    """
    
    ## get funtion value if not specifiv
    assert callable ( fun ) , "Invalid function `fun'"
    
    fx = float ( fun ( x ) ) if fx is None else fx  
    if not fx : return x                ## RETURN 
 
    assert callable ( deriv1 ), "Invalid first derivative `deriv1'"
    ## the  first derivative:
    d1 = float ( deriv1 ( x ) )
    if d1  : rn = fx / d1
    else   : return None                ## ERROR HERE!  
    
    ## make corrections: Halley's  steps
    if deriv2 :  
        assert callable ( deriv2 ) , "Invalid secon derivative `deriv2'"
        ## the second derivative
        d2 = float ( deriv2 ( x ) )
        if d2 : rn /= ( 1.0  - 0.5 * rn * d2 / d1 ) ## Halley's correction 
        
    return x - rn  ## Newton/Halley's iteration

# =============================================================================
## Single step of Steffensen's method
#  @see https://en.wikipedia.org/wiki/Steffensen%27s_method
#  @code
#  fun   = lambda x : -0.5*sin ( x )
#  x = 0.5
#  for i in  range ( 10 ) :
#      x = steffensen ( fun , x )
#      print i, x
#  @endcode 
def steffensen ( fun  ,  x  , fx = None ) :
    """ Single step of Steffensen's method

    Parameters
    ----------
    
    fun  : the function itself
    x    : the initial guess for the root
    fx   : (optional) function value at `x`, `fx=fun(x)`
  
    Returns
    -------
    The next approximation to the root or `None`

    Example
    -------
    
    >>> fun   = lambda x : -0.5*sin ( x )
    >>> x = 0.5
    >>> for i in  range ( 10 ) :
    ...      x = steffensen ( fun , x )
    ...      print i, x

    
    References
    ----------
    - https://en.wikipedia.org/wiki/Steffensen%27s_method

    """
    
    ## get funtion value if not specifiv
    assert callable ( fun ) , "Invalid function `fun'"

    print ( 'STEFFENSON HERE/1' , x , fx , fun ( x ) )
    
    fx = fun ( x ) if fx is None else fx  
    if not fx : return x                ## RETURN 

    print ( 'STEFFENSON HERE/2' , x , fx , fun ( x ) ) 

    gx = ( fun ( x + fx ) - fx ) / fx
    if not gx : return None             ## RETURN 
    
    print ( 'STEFFENSON HERE/3' , x , fx , gx , x - fx / gx ) 

    return x - fx / gx

# ============================================================================
## Single step of STFA method: Imporved regular falsi + steffenson-like
#  @see Xinyuan Wu, Zuhe Shen, Jianlin Xia, "An improved regula falsi method
#  with quadratic convergence of both diameter and point for enclosing
#  simple zeros of nonlinear equations", Applied Mathematics and Computation,
#  144, 2 2003
#  @see https://doi.org/10.1016/S0096-3003(02)00414-9},
#  @see https://www.sciencedirect.com/science/article/pii/S0096300302004149},
def STFA ( fun , a , b , c = None ) :
    """ Single step of STFA method:
    - Imporoved regular falsi + steffenson-like
    - see Xinyuan Wu, Zuhe Shen, Jianlin Xia, "An improved regula falsi method
    with quadratic convergence of both diameter and point for enclosing
    simple zeros of nonlinear equations", Applied Mathematics and Computation, 144, 2 2003
    - see https://doi.org/10.1016/S0096-3003(02)00414-9},
    - see https://www.sciencedirect.com/science/article/pii/S0096300302004149},
    """
    ## check brackets 
    if   0 == a.fx or iszero ( a.fx ) : return a , a , b
    elif 0 == b.fx or iszero ( b.fx ) : return b , a , b
    ## 
    ## (1) regular falsi step
    if c is None or not a.x < c.x < b.x : 
        c = ( a.x * b.fx - b.x * a.fx ) / ( b.fx - a.fx )
        c = Point ( c , fun ( c ) )
        ## 
    if   0 == c.fx or iszero ( c.fx ) : return c , a , b 
    elif samesign ( a.fx , c.fx )     : abar , bbar = c , b
    elif samesign ( b.fx , c.fx )     : abat , bbar = a , c
    else                              : return c , a , b 
    ##
    

# ===========================================================================
## Single step of Regular failsi-Bisection-Parabolic method
#
#  @see Alojz Suhadolnik, "Combined bracketing methods for solving nonlinear equations", ,
#       Applied Mathematics Letters}, 25, 11 (2012)
#  @see https://doi.org/10.1016/j.aml.2012.02.006},
#  @see https://www.sciencedirect.com/science/article/pii/S0893965912000778},
#
#  @see Somkid Intep, "A review of bracketing methods for finding zeros of nonlinear functions".
#       Applied Mathematical Sciences, Vol. 12, 2018, no. 3, 137-146
#  @see https://doi.org/10.12988/ams.2018.811
def RBP ( fun , a , b , c = None ) : 
    """ Single step of Regular failsi-Bisection-Parabolic method
    
    - see Alojz Suhadolnik, "Combined bracketing methods for solving nonlinear equations", ,
           Applied Mathematics Letters}, 25, 11 (2012)
    - see https://doi.org/10.1016/j.aml.2012.02.006},
    - see https://www.sciencedirect.com/science/article/pii/S0893965912000778},
    
    - see Somkid Intep, "A review of bracketing methods for finding zeros of nonlinear functions".
       Applied Mathematical Sciences, Vol. 12, 2018, no. 3, 137-146
    - see https://doi.org/10.12988/ams.2018.811
    """
    ## check brackets 
    if   0 == a.fx or iszero ( a.fx ) : return a , a , b
    elif 0 == b.fx or iszero ( b.fx ) : return b , a , b
    ## 
    ## (1) regular falsi step
    if c is None or not a.x < c.x < b.x : 
        c = ( a.x * b.fx - b.x * a.fx ) / ( b.fx - a.fx )
        if not a.x < c < b.x : c = 0.5 * ( a.x + b.x ) 
        c = Point ( c , fun ( c ) )
        if 0 == c.fx or iszero ( c.fx ) : return c , a , b     ## RETURN 
        print ( 'C is :' , c , a.x < c.x < b.x )
        
    ab = a.x - b.x
    ac = a.x - c.x
    bc = b.x - c.x 
    
    A  = ( a.fx - c.fx ) / ac
    A += ( c.fx - b.fx ) / bc 
    A /= ab 
    
    B  = ( c.fx - a.fx ) * bc / ac 
    B -= ( c.fx - b.fx ) * ac / bc   
    B /= ab
    
    C = c.fx

    pp = lambda x : A * x * x + B * x + C 
    
    print ( 'PARABOLA?' , pp ( a.x       ) , pp ( b.x      ) , pp ( c.x )       )
    print ( 'PARABOLA!' , pp ( a.x - c.x ) , pp ( b.x -c.x ) , pp ( c.x - c.x ) )
    
    D = B * B - 4 * A * C

    print ( 'HERE-3/' , A , B , C , D ) 
    if D < 0                       : return None , a , b ## RETURN 

    Q = B + signum ( B ) * math.sqrt ( D )
    p = c.x - 2 * C / Q
    
    print ( 'HERE-3/' , A , B , C , D , p ) 
    
    if not a.x < p < b.x            : return None , a , b  ## RETURN 
    ## 
    p = Point ( p , fun ( p ) )

    ## zero is found ? 
    if 0 == p.fx or iszero ( p.fx ) : return p , a , b     ## RETURN 

    ## update the interval
    
    if samesign ( a.fx , p.fx ) :
        a , b  = p , c if samesign ( b.fx , c.fx ) else b 
    else :
        b , a  = p , c if samesign ( a.fx , c.fx ) else a 

    dab = ( a.fx - b.fx ) / ( b.x - a.x )

    if 0.1 < abs ( dab ) < 10 :
        ## regular falsi 
        c = ( a.x * b.fx - b.x * a.fx ) / ( b.fx - a.fx )
    else :
        ## bisection 
        c = 0.5 * ( a.x + b.x )
        
    return Point ( c , fun ( c ) ) , a , b            ## RETURN 

# ============================================================================
## Single step of Muller's emthod
#  @see https://en.wikipedia.org/wiki/Muller%27s_method
def muller ( x0 , x1 , x2 , *other ) :
    """ Single step of Muller's method
    - see https://en.wikipedia.org/wiki/Muller%27s_method
    """
    points  = ( x0 , x1 , x2 ) + other
    x0 , x1 , x2 = points [ -3 : ]

    h0 = x1.x - x0.x
    h1 = x2.x - x1.x

    if not h0      : return None
    if not h1      : return None
    if not h1 + h1 : return None

    print ( 'MULLER 2' , h0 , h1 )
    
    d0 = ( x1.fx - x0.fx ) / h0 
    d1 = ( x2.fx - x1.fx ) / h1

    ## parabola's coefficients: 
    a  = ( d1 - d0 ) / ( h1 + h0 )
    b  = a * h1 + d1
    c  = x2.fx

    pp = lambda x : a * x * x + b * x + c

    print ( 'PARABOLA?' , pp ( x0.x  - x2.x ) , pp ( x1.x - x2.x ) , pp ( 0 ) )
    
    ## Parabola's discriminant 
    D  = b * b - 4 * a * c
    
    if D < 0 : return None  ## RETURN
    
    Q = b + signum ( b ) * math.sqrt ( D )
    
    if not Q : return None 
    
    dx = -2 * c / Q

    return x2.x + dx
     
# ============================================================================
## make an inverse linear interpolation, aka "secant", aka "regular falsi"
#  @see https://en.wikipedia.org/wiki/Secant_method
#  @see https://en.wikipedia.org/wiki/False_position_method
#  @code
#  xa , xb = ...
#  fun     = ...
#  r = inverse_linear ( Point ( xa , fun ( xa ) ) ,
#                       Point ( xb , fun ( xb ) ) ) 
#  r = secant         ( Point ( xa , fun ( xa ) ) ,
#                       Point ( xb , fun ( xb ) ) ) 
#  r = regular_falsi  ( Point ( xa , fun ( xa ) ) ,
#                       Point ( xb , fun ( xb ) ) ) 
#  @endcode 
def inverse_linear ( a , b , *other ) :
    """ Make a linear interpolation, aka 'secant', aka 'regular falsi'

    Parameters
    ----------
    a :  `xa,f(xa)`-point
    b :  `xb,f(xb)`-point

    Returns
    -------
    (Inverse) linear/secant/regular falsi approximation for the root or `None`

    Example
    -------
    
    >>> xa , xb = ...
    >>> fun     = ...
    >>> r = inverse_linear ( Point ( xa , fun ( xa ) ) ,
    ...                      Point ( xb , fun ( xb ) ) ) 
    >>> r = secant         ( Point ( xa , fun ( xa ) ) ,
    ...                      Point ( xb , fun ( xb ) ) ) 
    >>> r = regular_falsi  ( Point ( xa , fun ( xa ) ) ,
    ...                      Point ( xb , fun ( xb ) ) ) 

    References
    ----------
    - https://en.wikipedia.org/wiki/Secant_method
    - https://en.wikipedia.org/wiki/False_position_method
  
    """
    points = ( a , b ) + other
    b , a = points [ -2 : ]
    
    x0 , f0 = a.x , a.fx
    x1 , f1 = b.x , b.fx

    if   not f0 : return x0
    elif not f1 : return x1 
    
    if f0 == f1 or isequal ( f0 , f1 ) : return None

    return ( x0 * f1  - x1 * f0 ) / ( f1 - f0 ) 

secant        = inverse_linear 
regular_falsi = inverse_linear

# =============================================================================
## Inverse parabolic interppolation
#  @see https://en.wikipedia.org/wiki/Inverse_quadratic_interpolation
#  @code
#  xa , xb , xc = ...
#  fun          = ...
#  r = inverse_parabolic ( Point ( xa , fun ( xa ) ) ,
#                          Point ( xb , fun ( xb ) ) ,
#                          Point ( xc , fun ( xc ) ) )
#  @endcode
def inverse_parabolic ( a , b , c , *other ) :
    """ Inverse parabolic interpolation via the points:
    - ( xa , fun ( xa ) ), ( xb , fun ( xb ) ), ( xc , fun ( xc ) )

    Parameters
    ----------
    a :  `xa,f(xa)`-point
    b :  `xb,f(xb)`-point
    c :  `xb,f(xc)`-point

    Returns
    -------
    Inverse parabolic approximation for the root or `None`

    Example
    -------
    
    >>> xa , xb , xc = ...
    >>> fun          = ...
    >>> r = inverse_parabolic ( Point ( xa , fun ( xa ) ) ,
    ...                         Point ( xb , fun ( xb ) ) ,
    ...                         Point ( xc , fun ( xc ) ) )

    References
    ----------
    - https://en.wikipedia.org/wiki/Inverse_quadratic_interpolation
  
    """
    points = ( a , b , c ) + other
    c , b , a = points [ -3 : ]
    
    x0 , f0 = a.x , a.fx
    if not f0 : return x0 
    
    x1 , f1 = b.x , b.fx
    if not f1 : return x1
    
    x2 , f2 = c.x , c.fx
    if not f2 : return x2 
    
    if   f0 == f1 or isequal ( f0 , f1 ) : return inverse_linear ( a , c )  
    elif f0 == f2 or isequal ( f0 , f2 ) : return inverse_linear ( a , b ) 
    elif f1 == f2 or isequal ( f1 , f2 ) : return inverse_linear ( a , b ) 
    
    f01 = 1.0 / ( f0 - f1 ) ;  f10 = - f01
    f02 = 1.0 / ( f0 - f2 ) ;  f20 = - f02
    f12 = 1.0 / ( f1 - f2 ) ;  f21 = - f12

    xx  = x0 * f1 * f2 * f01 * f02
    xx += x1 * f0 * f2 * f10 * f12
    xx += x2 * f0 * f1 * f20 * f21
    
    return xx

# =============================================================================
## Inverse cubic interppolaiton
#  @code
#  xa , xb , xc , xd = ...
#  fun          = ...
#  r = inverse_cubic ( Point ( xa , fun ( xa ) ) ,
#                      Point ( xb , fun ( xb ) ) ,
#                      Point ( xc , fun ( xc ) ) ,
#                      Point ( xd , fun ( xd ) ) )
#  @endcode
def inverse_cubic ( a , b , c , d , *other ) :
    """ Inverse cubic interpolaton via the last four approximations to the root

    Parameters
    ----------
    
    a :  `xa,f(xa)`-point
    b :  `xb,f(xb)`-point
    c :  `xb,f(xc)`-point
    d :  `xd,f(xd)`-point

    Returns
    -------
    Inverse cubic approximation for the root or `None`

    Example
    -------
    
    >>> xa , xb , xc , xd = ...
    >>> fun               = ...
    >>> r = inverse_cubic ( Point ( xa , fun ( xa ) ) ,
    ...                     Point ( xb , fun ( xb ) ) ,
    ...                     Point ( xc , fun ( xc ) ) ,
    ...                     Point ( xd , fun ( xd ) ) )

    
    """
    points = ( a , b , c , b ) + other
    d , c , b , a = points [ -4 : ]
    
    x0 , f0 = a.x , a.fx
    if not f0 : return x0 
    
    x1 , f1 = b.x , b.fx
    if not f1 : return x1 
    
    x2 , f2 = c.x , c.fx
    if not f2 : return x2 
    
    x3 , f3 = d.x , d.fx
    if not f3 : return x3 
    
    ## switch to inverse parabolic if some function values coincide 
    if   f0 == f1 or isequal ( f0 , f1 ) : return inverse_parabolic ( a , c , d ) 
    elif f0 == f2 or isequal ( f0 , f2 ) : return inverse_parabolic ( a , b , d )
    elif f0 == f3 or isequal ( f0 , f3 ) : return inverse_parabolic ( a , b , c )
    elif f1 == f2 or isequal ( f1 , f2 ) : return inverse_parabolic ( a , b , d ) 
    elif f1 == f3 or isequal ( f1 , f3 ) : return inverse_parabolic ( a , c , d ) 
    elif f2 == f3 or isequal ( f2 , f3 ) : return inverse_parabolic ( a , b , d ) 

    f01 = 1.0 / ( f0 - f1 ) ;  f10 = - f01
    f02 = 1.0 / ( f0 - f2 ) ;  f20 = - f02
    f03 = 1.0 / ( f0 - f3 ) ;  f30 = - f03
    f12 = 1.0 / ( f1 - f2 ) ;  f21 = - f12
    f13 = 1.0 / ( f1 - f3 ) ;  f31 = - f13
    f23 = 1.0 / ( f2 - f3 ) ;  f32 = - f23

    xx  = -x0 * f1 * f2 * f3 * f01 * f02 * f03 
    xx += -x1 * f0 * f2 * f3 * f10 * f12 * f13 
    xx += -x2 * f0 * f2 * f3 * f20 * f21 * f23 
    xx += -x3 * f0 * f1 * f2 * f30 * f31 * f32
    
    return xx

# ============================================================================
## make an inverse polynomial interpolation 
#  - Actually switch to cubic/parabolic or linear interpolation
#  @code
#  xa , xb , ... , xn = ...
#  fun     = ...
#  r = inverse_polynomial ( Point ( xa , fun ( xa ) ) ,
#                           Point ( xb , fun ( xb ) ) ,
#                           ...                       ,
#                           Point ( xn , fun ( xn ) ) )
#  @endcode 
def inverse_polynomial ( a , b , *other ) :
    """ Make an inverse polynomial interpolation
    - Actually switch to cubic/parabolic or linear interpolation

    Parameters
    ----------
    
    a     :  `xa,f(xa)`-point
    b     :  `xb,f(xb)`-point
    other : other points   

    Returns
    -------
    
    Inverse polynomial approximation for the root or `None`

    Example
    --------
    
    >>> xa , xb , ... , xn = ...
    >>> fun                = ...
    >>> r = inverse_polynomial ( Point ( xa , fun ( xa ) ) ,
    ...                      Point ( xb , fun ( xb ) ) ,
    ...                      ...                       ,
    ...                      Point ( xn , fun ( xn ) ) )
    
    """
    
    if   4 <= 2 + len ( other ) : return inverse_cubic     ( a , b , *other )
    ## elif 3 <= 2 + len ( other ) : return inverse_parabolic ( a , b , *other )
    
    return inverse_linear ( a , b , *other)

# =============================================================================
## trivial bisection method
#  @code
#
#  xa , xb = ...
#  fun     = ...
#  a  , b  = Point ( xa  , fun ( xa ) ) , Point ( xb  , fun ( xb ) )
#  a , n   = bisection ( fun , a , b )
#
#  @endcode 
#  @see https://en.wikipedia.org/wiki/Bisection_method
def bisection ( fun , a  , b ) :
    """ Trivial bisection method

    Parameters
    ----------
    fun   : the function
    a     :  `xa,f(xa)`-point
    b     :  `xb,f(xb)`-point

    Returns
    -------
    New interval with located root

    Example
    -------
    
    xa , xb = ...
    fun     = ...
    a  , b  = Point ( xa  , fun ( xa ) ) , Point ( xb  , fun ( xb ) )
    a , n   = bisection ( fun , a , b )

    References
    ----------
    - see https://en.wikipedia.org/wiki/Bisection_method
    """
    
    a , b = ( a , b ) if a.x <= b.x else ( b , a )

    c  = 0.5 * ( a.x + b.x )
    fc = fun   ( c ) 
    c  = Point ( c , fc )
    
    if   not fc or iszero ( fc ) : return c , c ## RETURN 
    elif samesign ( a.fx , fc )  : return c , b ## RETURN
    elif samesign ( b.fx , fc )  : return a , c ## RETURN
    
    return c , c 

# =============================================================================
## try to find next approximation using Aitken delta-squared process
#  @code
#
#  xcurr , xm1 , xm2 = ...
#  r = aitken_delta2 ( xcurr , xm1 , xm2 )
#
#  @endcode 
#  @see https://en.wikipedia.org/wiki/Aitken%27s_delta-squared_process
def aitken_delta2 ( x0 , x1 , x2 , *others  ) :
    """ Try to find next approximation to the root using Aitken delta-squared process.
    
    It often helps hear the multiple root, where all sophisticated schemes result
    in very  slow (linear) convergency 
    
    Parameters
    -----------
    
    xl  :  the last root estimate              (`n`)
    xl1 :  the previous root estimate          (`n-1`)
    xl2 :  the previous-previous root estimate (`n-2`) 
    
    Return
    ------
    
    Improvement in the root position or `None`

    Example
    -------
    
    >>> xcurr , xm1 , xm2 = ...
    >>> r = aitken_delta2 ( xcurr , xm1 , xm2 ) 

    References
    ----------
    - https://en.wikipedia.org/wiki/Aitken%27s_delta-squared_process
    """
    points = ( x0 , x1 , x2 ) + others 
    xn , xn1 , xn2 = [ p.x for p in points [ -3 : ] ] 

    dd  = ( xn - xn1 ) + ( xn2 - xn1 )
    if not dd or iszero ( dd ) : return None
    return ( xn * xn2 - xn1 * xn1 )  / dd 

# ========================================================================================
## @class RootResults
#  Helper class to keep results of root-finding procedure
#  - It is very similar (almost clone) of corresponding class
#  <code>RootResults</code> from <code>scipy.optimize.zeros</code>
class RootResults ( object ):
    """ Helper class to keep results of root-finding procedure
     - It is very similar (almost clone) of corresponding class
    `RootResults` from scipy.optimize.zeros
    """
    __slots__ = ( 'status'            ,
                  'flag'              ,
                  'function_calls'    ,
                  'derivative1_calls' ,
                  'derivative2_calls' ,                  
                  'iterations'        , 
                  'root'              ,
                  'value'             ,
                  'bracket'           , 
                  'counters'          )
    
    CONVERGED = 'converged'
    SIGNERR   = 'sign error'
    CONVERR   = 'convergence error'
    flag_map  = { 0 : CONVERGED, -1 : SIGNERR , -2 : CONVERR }
    
    def __init__  ( self                     ,
                    root                     ,
                    iterations               ,
                    function_calls           , * ,
                    derivative1_calls = 0    ,
                    derivative2_calls = 0    ,                                      
                    flag              = 0    , 
                    value             = None ,
                    bracket           = ()   , 
                    counters          = {}   ) : 
        
        self.root              = root
        self.iterations        = iterations
        self.function_calls    = function_calls
        self.derivative1_calls = derivative1_calls
        self.derivative2_calls = derivative2_calls
        self.flag              = flag
        self.status            = self.flag_map.get ( flag , '*UNKNOWN*' )
        self.value             = value
        self.counters          = counters
        self.bracket           = bracket 
    # ====================================================================
    ## format it as the table 
    def table   ( self             ,
                  precision = 6    ,
                  width     = 8    ,
                  title     = ''   ,
                  prefix    = ''   ,
                  style     = None ) :
        """ Format it as a table
        """
        header = '' , 'value' , '' 
        rows   = [ header ]

        r , e  = pretty_float ( self.root , precision = precision , width = width , with_sign = True )
        row    = 'Root' , r , '10^%+d' % e if e else ''
        rows.append ( row )

        row    = 'Status' , attention ( self.status ) if self.flag else allright ( self.status ) 
        rows.append ( row )

        if self.bracket :
            a , b = self.bracket
            fmtv , expo = fmt_pretty_values ( a , b ,
                                              precision = precision , 
                                              width     = width     ,
                                              with_sign = True      )
            if expo : 
                fa = fmtv % ( a / 10**expo )
                fb = fmtv % ( b / 10**expo )
                v  = '[%s,%s]' % ( fa , fb )
                row = 'Bracket' , v , '10^%+d' % expo
            else :
                fa = fmtv % ( a )
                fb = fmtv % ( b )
                v  = '[%s,%s]' % ( fa , fb )
                row = 'Bracket' , v , '' 
            rows.append ( row ) 

        if not self.value is None :            
            r , e  = pretty_float ( self.value , precision = precision , width = width , with_sign = True )
            row    = 'Function value' , r , '10^%+d' % e if e else ''            
            rows.append ( row )

        if 0 < self.iterations : 
            row = 'Iterations' , '%d' % self.iterations
            rows.append ( row )

        if 0 < self.function_calls :
            row = 'Function calls' , '%d' % self.function_calls
            rows.append ( row )
            
        if 0 < self.derivative1_calls :
            row = 'Derivative1 calls' , '%d' % self.derivative1_calls
            rows.append ( row )
            
        if 0 < self.derivative2_calls :
            row = 'Derivative2 calls' , '%d' % self.derivative2_calls
            rows.append ( row )

        if self.counters :
            row = 'Methods' , ' success / failure'
            rows.append ( row ) 
            for key,cnt  in self.counters.items() :
                success = cnt.success
                failure = cnt.failure 
                row     = key , '%7d /%- 7d' % ( success , failure )
                rows.append ( row ) 

        rows  = T.remove_empty_columns ( rows )
        title = title if title else 'Root-finder result'
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lcc' , style = style )
    
    __repr__ = table
    __str__  = table
        
# ===========================================================================================
## @class RootFinder
#  The main root-finding engine.
#
#  For each step it use sequence of well-defined "actions" The result of each "action"
#  is considered to be "successfull" if the obtained root approximation falls into the current
#  bracketing interval, otherwise it is considered as "failure"
# 
#  In the successfull case the obtained root it stored into list, and the bracketing
#  interval is correspodingly updated.

#  We start iterations from adding the current brackets into collection of roots
# 
#  The sequential list of actions is:
#  - 'Halley'     :                             if the second derivative is provided, make Halley's step
#  - 'Newton'     : if previous step is failed  and the first derivative is provdied, make Newton step
#  - 'Steffensen' : if previous step is failed                                      , make Steffensen's step
#  - 'Secant'     : if previous step is failed (it should not, but numerical failures can occur)
#  - 'Bisection'  : if previous step is failed 
#  - 'Inversenterpolation' : at this point we alwasy have at least two or more  approximations
#    for the root stored into the list (including initial brackets), so an inverse polynomial interpolation is performed: 
#     -- if              there are >=4 points, try with cubic
#     -- if it fails and there are >=3 points, try with parabolic
#     -- if it fails,                          try with secant   (and this one can't fail!)    
#  - 'Aitken'     : if there are >=3 approximations, try to use Aitken delta2 acceleratiion scheme  
#  - 'Bisection'  : if the updated bracketing interval is not shrinked significantly for the current step, a bisection step is forced 
#
#  After the step the convergency is checked:
#  - either the value of the  function at the currect root  approximation is small enough
#  - or distance between two subsequent approximations is very small 
#  - or the updated bracketing interval is very small
#
#  Note, that:
#  - similar to Brent's and TOMS74 algortihms, due to the 'Bisection' action the algorithm converges ``almost always''.
#  - due to intermediate 'Aitken' step, a slow linear or sublinear convergency (e.g. in vicinity of multiple root) is often drastically improved
# 
#  @see https://en.wikipedia.org/wiki/Halley%27s_method
#  @see https://en.wikipedia.org/wiki/Newton%27s_method
#  @see https://en.wikipedia.org/wiki/Steffensen%27s_method
#  @see https://en.wikipedia.org/wiki/Secant_method
#  @see https://en.wikipedia.org/wiki/False_position_method
#  @see https://en.wikipedia.org/wiki/Inverse_quadratic_interpolation
#  @see https://en.wikipedia.org/wiki/Aitken%27s_delta-squared_process
#  @see https://en.wikipedia.org/wiki/Bisection_method
#
#  The algorithm is largely inspired  by the Brent's method and newer TOMS748 algorithm from <code>Boost.math rootfinding</code> 
#  @see https://en.wikipedia.org/wiki/Brent%27s_method
#  @see http://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf
#  @see https://www.boost.org/doc/libs/1_67_0/boost/math/tools/toms748_solve.hpp
#  @see https://www.boost.org/doc/libs/release/libs/math/doc/html/root_finding.html
#  @see https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/roots_noderiv/bracket_solve.html
#  @see https://www.boost.org/doc/libs/release/boost/math/tools/toms748_solve.hpp
class RootFinder(object) :
    """ The main root-finding engine.
    
    For each step it use sequence of well-defined `actions'.
    The result of each `action' is considered to be `successfull',
    if the obtained root approximation falls into the current
    bracketing interval, otherwise it is considered as `failure'.
    
    In the successfull case the obtained root it stored into list, and the bracketing
    interval is correspodingly updated.
    
    We start iterations from adding the current brackets into collection of roots
    
    The sequential list of `actions' is:
    -  'Halley'     :                              if the second derivative is provided, make Halley's step
    -  'Newton'     : if previous step is failed  and the first  derivative is provided, make Newton step
    -  'Steffensen' : if previous step is failed                                       , make Steffensen's step
    -  'Secant'     : if previous step is failed (it should not, but numerical failures can occur)
    -  'Bisection'  : if previous step is failed 
    -  'InversInterpolation' : at this point we alwasy have at least two or more  approximations
    for the root stored into the list (including initial brackets), so an inverse polynomial interpolation is performed: 
    -- if              there are >=4 points, try with cubic
    -- if it fails and there are >=3 points, try with parabolic
    -- if it fails,                          try with secant   (and this one can't fail!)    
    -  'Aitken'     : if there are >=3 approximations, try to use Aitken delta2 acceleratiion scheme  
    -  'Bisection'  : if the updated bracketing interval is not shrinked significantly for the current step, a bisection step is forced 
    
    After the step the convergency is checked and the proicess is converged  if:
    - either the value of the function at the currect root  approximation is small enough
    - or the distance between two subsequent approximations is very small 
    - or the updated bracketing interval is very small
    
    Note, that:
    - similar to Brent's and TOMS748 algortihms, due to the 'Bisection' action the algorithm converges ``almost always''.
    - due to intermediate 'Aitken' step, a slow linear (or sublinear) convergency
      (e.g. in vicinity of multiple root) is often drastically improved
    
    The algorithm is largely inspired  by the Brent's method and newer TOMS748 algorithm from Boost.math root finding library

    References
    ----------
    
    - see https://en.wikipedia.org/wiki/Halley%27s_method
    - see https://en.wikipedia.org/wiki/Newton%27s_method
    - see https://en.wikipedia.org/wiki/Steffensen%27s_method
    - see https://en.wikipedia.org/wiki/Secant_method
    - see https://en.wikipedia.org/wiki/False_position_method
    - see https://en.wikipedia.org/wiki/Inverse_quadratic_interpolation
    - see https://en.wikipedia.org/wiki/Aitken%27s_delta-squared_process
    - see https://en.wikipedia.org/wiki/Bisection_method
    
    - see https://en.wikipedia.org/wiki/Brent%27s_method
    - see http://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf
    - see https://www.boost.org/doc/libs/1_67_0/boost/math/tools/toms748_solve.hpp
    - see https://www.boost.org/doc/libs/release/libs/math/doc/html/root_finding.html
    - see https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/roots_noderiv/bracket_solve.html
    - see https://www.boost.org/doc/libs/release/boost/math/tools/toms748_solve.hpp

    Performance
    -----------

    - `RootFinder` is able to find roots for many cases, where `scipy.optimize.brentq` fails,
    in particular for the multiple roots, e.g.

    >>> N = 11
    >>> fun  = lambda x : 1.0 * ( x - 1.0 )**N

    for this function the local algorithm succesfully finds the root:
    
    >>> print 'Plain :' , find_root ( fun , 0.5 , 10 , full_output = True , disp = True )

    ... Plain : (0.9999997973130655,             converged: True
    ...                 flag: 'converged'
    ...       function_calls: 68
    ...           iterations: 14
    ...                 root: 0.9999997973130655
    ...           Steffensen: (14, 2)
    ...               Secant: (12, 12)
    ... InverseInterpolation: (14, 6)
    ...               Aitken: (11, 3)
    ...            Bisection: (14, 14))

    while the standard Brent's algorithm fails:
    
    >>> import scipy.optimize as SP
    >>> print 'Brent :' , SP.brentq ( fun , 0.5 , 10 , full_output = True , disp = True )

    ... Traceback (most recent call last):
    ... <...> 
    ... r = _zeros._brentq(f,a,b,xtol,rtol,maxiter,args,full_output,disp)
    ... RuntimeError: Failed to converge after 100 iterations.

    - for many cases is makes less iterations and in spite of the larger number of calls
    per iteration, very often the total number of functon calls is somewhat smaller: 
    
    >>> N = 7
    >>> fun  = lambda x :   1.0 * ( x - 1.0 )**N
    >>> print 'Plain :' , find_root ( fun , 0.5 , 10 , full_output = True , disp = True )
    
    ... Plain : (0.9999964237209732,             converged: True
    ...                 flag: 'converged'
    ...       function_calls: 59                   ## <--- HERE 
    ...           iterations: 12                   ## <--- HERE 
    ...                 root: 0.9999964237209732
    ...           Steffensen: (12, 3)
    ...               Secant: (9, 9)
    ... InverseInterpolation: (12, 4)
    ...               Aitken: (8, 4)
    ...            Bisection: (12, 12))

    >>> print 'Brent :' , SP.brentq ( fun , 0.5 , 10 , full_output = True , disp = True )
    
    ... Brent : (0.9999999999356982,       converged: True
    ...           flag: 'converged'
    ... function_calls: 102           ## <--- HERE
    ...     iterations: 100           ## <---- HERE
    ... root: 0.9999999999356982)
    
    - however, for the trivial functions overall CPU performance is factor of 2 worse
    than the CPU performance of the standard Brent's algorithm from `scipy.optimize`
    Likely it is a direct sequence of C-implementation of `scipy.optimize.brentq`
    versus pure python implementation of this root-finder for this version of `scipy`:
    `scipy.version.version='0.18.1'`
    For complicated functions with expensive evaluation, the CPU performance
    is well comparable or even better that the CPU performance of `scipy.optimize.brentq`
    
    """
    # =========================================================================
    def __init__ ( self                ,
                   fun                 , * , 
                   deriv1      = None  , 
                   deriv2      = None  ,
                   maxiter     = 500   ,
                   args        = ()    , ## optional positional arguments 
                   kwargs      = {}    , ## optional keyword argumenst 
                   xtol        = _xtol ,
                   rtol        = _rtol ,
                   full_output = False ,
                   disp        = True  ) :
        """ Create the root-finder

        Parameters
        ----------
        fun         : the function
        deriv1      : (optionaly) the 1st derivative
        deriv2      : (optionaly) the 2nd derivative
        maxiter     : maximal number of iterations
        full_output : If `full_output` is `False`, the root is returned.
                      If `full_output` is `True`, the return value is `(x, r)`, where `x` is the root, and `r` is a `RootResults` object.
        disp        : If `True`, raise `RuntimeError` if the algorithm didn’t converge.
        """

        assert callable ( fun ) , "Invalid function `fun'"
        
        self.__full_output = full_output

        self.__args   = args   if args   else ()
        self.__kwargs = kwargs if kwargs else {} 
        
        has_args = True if self.args or self.kwargs else False 
        if has_args : ffun = lambda x : fun ( x, *self.args, **self.kwargs )
        else        : ffun = fun 
        
        fderiv1 = None
        fderiv2 = None 
        
        if deriv1 : 
            assert callable ( deriv1 ) , "Invalid first derivative `deriv1'"
            if has_args : fderiv1 = lambda x : deriv1 ( x , *self.args , **self.kwargs )
            else        : fderiv1 = deriv1
            if deriv2 : 
                assert callable( deriv2 ) , "Invalid second derivative`deriv2'"
                if has_args : fderiv2 = lambda x : deriv2 ( x , *self.args , **self.kwargs )
                else        : fderiv2 = deriv2
     
     
        self.__fun     = counted ( ffun    ) if self.__full_output             else ffun  
        self.__deriv1  = counted ( fderiv1 ) if self.__full_output and fderiv1 else fderiv1 
        self.__deriv2  = counted ( fderiv2 ) if self.__full_output and fderiv2 else fderiv2

        self.__xtol    = max ( xtol , _xtol )   
        self.__rtol    = max ( rtol , _rtol )
        
        self.__maxiter = maxiter if 1 <= maxiter else 500
        
        self.__disp    = True    if disp        else False
        
        ## keep several previous approximations for inverse polynomial, Muller or Aitken
        self.__roots   = [] 

        self.__counters = collections.defaultdict ( EffCounter )
        
    @property
    def args ( self ) :
        """`args`: optional positional arguments for fun/deriv1/deriv2
        """
        return self.__args 
    @property
    def kwargs ( self ) :
        """`kwargs`: optional keyword arguments for fun/deriv1/deriv2
        """
        return self.__kwargs 
    
    # =========================================================================
    ## Return result/root with optional report 
    def __result ( self       ,
                   root      , * ,
                   iteration = 0    ,
                   value     = None ,
                   bracket   = ()   , 
                   flag      = 0    ) :
        """ Return result/root with optional report 
        """
        if not self.__full_output : return root

        return root , RootResults ( root ,
                                    iterations        = iteration        ,
                                    function_calls    = self.__fun.calls ,
                                    derivative1_calls = self.__deriv1.calls if self.__deriv1 else 0 ,
                                    derivative2_calls = self.__deriv2.calls if self.__deriv2 else 0 ,
                                    value             = value           , 
                                    flag              = flag            ,
                                    bracket           = bracket         , 
                                    counters          = self.__counters )

    # =========================================================================
    ## solve the function and find the root 
    def find ( self , a , b , guess = None ) :
        """ Find a root
        
        Parameters
        ----------
        
        a     : low  edge of bracketing interval
        b     : high edge of bracketing interval
        guess : guess on the root position 
        
        Returns
        -------
        root or (root,report)  
        
        Example
        -------
        
        >>> rf = RootFinder ( ....   )
        >>> r  = rf.find    ( 0 , 10 )
        """

        ## swap them if needed 
        a , b = ( a , b ) if a < b else ( b , a )

        ## check the root at left edge 
        fa = self.__fun ( a )
        if 0 == fa or iszero ( fa ) : return self.__result ( a , value = fa , bracket = ( a , b ) ) ## RETURN

        ## check the root at right edge 
        fb = self.__fun ( b )
        if 0 == fb or iszero ( fb ) : return self.__result ( b , value = fb , bracket = ( a , b  ) ) ## RETURN

        ## correct bracketing interval ?
        if samesign ( fa ,  fb ) :
            if self.__disp :
                raise RuntimeError ('RootFinder.solve: invalid bracketing interval (%s,%s)/(%s,%s)' % ( a , fa , b , fb ) )
            else           :
                return self.__result ( 0.5 * ( a  + b ) , flag = -1 , bracket = ( a , b ) )  ## RETURN
            
        a = Point ( a ,  fa )
        b = Point ( b ,  fb )

        ## initialize 
        self.__roots = [ a , b ] 
        
        ## Guess is not specified or invalid 
        if guess is None or not a.x <= guess <= b.x :
            
            ## use the secant as the approximation 
            guess = secant ( a , b )
            ok    = not guess is None and isfinite ( guess ) and a.x <= guess <= b.x
            self.__counters [ 'secant' ] += ok
            if not ok :
                guess = 0.5 * ( a.x + b.x )
                self.__counters [ 'bisection' ] += 1

        ## if valid numerical value 
        assert isinstance ( guess , num_types ) and a.x <= guess <= b.x , "Invalid `guess' %s" % guess 
        
        ## calcualate fun(guess) 
        guess = Point ( guess , self.__fun ( guess ) )

        ## check the function value and shrink the interval, if/when  possible 
        if 0 == guess.fx or iszero ( guess.fx  ) : return self.__result ( guess.x , value = guess.fx , bracket = ( a.x , b.x ) ) ##  RETURN
        elif samesign ( a.fx , guess.fx ) : a , b = guess , b 
        elif samesign ( b.fx , guess.fx ) : a , b = a     , guess
        else                        : return self.__result ( guess.x , value = guess.fx , bracket = ( a.x , b.x ) ) ##   RETURN
        
        ## initialize Aitken and Inverse Interpolation structures 
        self.__roots  = [ a , b ] 
        self.__aitken = [ a , b ]   
        
        print ( 'FIND/0' , a , b )

        x0 = guess.x 
        for i in range ( self.__maxiter ) :

            ## get new x and new bracketing interval 
            x , a , b = self.__make_step ( x0 , a , b )

            print ( 'FIND/i/1' , i , x ,  a , b )
            
            ## zero is found ?
            if not x.fx or iszero ( x.fx ) : return self.__result ( x.x , iteration = i + 1 , value = x.fx , bracket = ( a.x , b.x ) ) 
        
            print ( 'FIND/i/2' , i , x ,  a , b )

            ## dx1 = self.__xtol + self.__rtol * max ( abs ( x.x ) , abs ( x0 ) )            
            ## if abs ( x.x - x0 ) * 4 < dx1 : return self.__result ( x.x , iteration = i + 1 , value = x.fx )
            
            print ( 'FIND/i/3' , i , x ,  a , b )

            ## bracketing interval is already very small 
            dx2 = self.__xtol + self.__rtol * abs ( a.x - b.x )
            if abs ( b.x - a.x   ) * 4 < dx2 : return self.__result ( x.x , iteration = i + 1 , value = x.fx , bracket = ( a.x , b.x ) ) 
            
            print ( 'FIND/i/4' , i , x ,  a , b )

            x0 = x.x
            
        message = 'RootFinder.solve: no convergency for %d iterations' % i 
        if self.__disp : raise RuntimeError ( message )
        else           : logger.warning     ( message ) 
            
        return self.__result ( x.x , iteration = i , flag = -2 , value = x.fx , bracket = ( a.x , b.x ) ) 

    # =========================================================================
    def __make_step ( self , x0 , a , b ) :

        ## do not keep too many approximations
        
        self.__roots = self.__roots [ -5 : ]
        print  ( 'MAKE STEP' , self.__roots )
        
        ## 
        old_len =  b.x - a.x
        
        x   = None
        fx0 = None

        # =====================================================================
        ## (1,2) try Halley's and/or Newton methods 
        if x is None and self.__deriv1 :            
        
            if fx0 is None               : fx0 = self.__fun ( x0 )
            if not fx0 or iszero ( fx0 ) : return x0 , a , b 
            ## 
            xx = halley_newton ( self.__fun    ,
                                 x0            ,
                                 self.__deriv1 , 
                                 deriv2 = self.__deriv2 , fx = fx0 )
            ## OK ?
            ok = not xx is None and isfinite ( xx ) and a.x <= xx <= b.x 
            if self.__deriv2 : self.__counters [ 'halley' ] += ok  
            else             : self.__counters [ 'newton' ] += ok  
            if ok : 
                x = Point ( xx , self.__fun ( xx ) ) 
                self.__roots.append ( x  )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN                
                print ( 'Newton:' , x , a , b , self.__roots ) 

        """
        # =====================================================================
        ## (3) try  Steffensen's method
        if x is None or not a.x <= x.x <= b.x : 
            ##
            if fx0 is None               : fx0 = self.__fun ( x0 ) 
            if not fx0 or iszero ( fx0 ) : return x0 , a , b
            ##
            if a.x < x0 + fx0 < b.x : 
                ##
                xx = steffensen ( self.__fun , x0 , fx = fx0 )
                ## OK ?
                ok = not xx is None and isfinite ( xx ) and a.x <= xx <= b.x
                if ok : 
                    x = Point ( xx , self.__fun ( xx ) ) 
                    self.__roots.append ( x  )
                    if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                    elif samesign ( a.fx , x.fx )     : a , b = x , b  
                    elif samesign ( b.fx , x.fx )     : a , b = a , x  
                    else                              : return x , a , b ## RETURN                    
                    print ( 'Steffensen' , x , a , b )
            else :
                print ( 'No steffensen', a.x , x0 + fx0 , b.x ) 
        """
        
        # =====================================================================
        ## (4) Try Muller's method
        if 3 <= len ( self.__roots ) and ( x is None or not a.x < x.x <= b.x ) :
            print ( 'BEFORE MULLER' , self.__roots) 
            xx = muller ( *self.__roots )
            ok = not xx is None and isfinite ( xx ) and a.x <= xx <= b.x
            self.__counters [ 'muller' ] += ok
            if ok :
                x = Point ( xx , self.__fun ( xx ) ) 
                self.__roots.append ( x  )
                print ( 'MULLER: added ' , x )
                if   0 == x.fx or iszero ( x.fx )     : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )         : a , b = x , b  
                elif samesign ( b.fx , x.fx )         : a , b = a , x  
                else                                  : return x , a , b ## RETURN
                print ( 'MULLER' , x , a , b , self.__roots ) 
                                            
        # =====================================================================
        ## (5) Inverse polynomial interpolation: we always have some points here
        if False : ## ( x is None or not a.x <= x.x <= b.x ) and 2 <= len ( self.__roots ) :
                  
            xx = inverse_polynomial ( *self.__roots )
            ok = not xx is None and isfinite ( xx ) and a.x <= xx <= b.x
            self.__counters [ 'inverse' ] += ok
            if ok : 
                x = Point ( xx , self.__fun ( xx ) ) 
                self.__roots.append ( x  )
                if   0 == x.fx or iszero ( x.fx )     : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )         : a , b = x , b  
                elif samesign ( b.fx , x.fx )         : a , b = a , x  
                else                                  : return x , a , b ## RETURN
                print ( 'INVERSE' , x , a , b , self.__roots ) 

        # =====================================================================
        ## (6) secant is "almost never" fails 
        if x is None or not a.x <= x.x <= b.x :
            
            xx = secant ( a , b )
            ## OK ?
            ok = not xx is None and isfinite ( xx ) and a.x <= xx <= b.x 
            self.__counters [ 'secant' ] += ok
            if ok : 
                x = Point ( xx , self.__fun ( xx ) ) 
                self.__roots.append ( x  )
                print ( 'SECANT: added ' , x )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN                 
                print ( 'SECANT' , x , a , b , self.__roots )

        # =====================================================================
        ## (7) Try Aitken-delta2 scheme if we have enough approximations
        if 3 <= len ( self.__aitken ) : 

            xx = aitken_delta2 ( *self.__roots )
            ok = not xx is None and isfinite ( xx ) and a.x <= xx <= b.x
            ## accept Aitken only if it improves the estimate 
            if ok :
                fxx = self.__fun ( xx )
                if not x is None :
                    ok  = abs ( fxx ) <= abs ( x.fx )                
            print ( 'AITKEN/2' , self.__roots , xx , ok )            
            self.__counters [ 'aitken' ] += ok
            if ok :
                ## x = Point ( xx , self.__fun ( xx ) ) 
                x = Point ( xx , fxx ) 
                self.__roots.append ( x  )
                print ( 'AITKEN: added' , x ) 
                ## self.__aitken .append ( xx )            
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x
                else                              : return x , a , b ## RETURN                
                print ( 'Aitken:' , x , a , b , self.__roots )

        # =====================================================================        
        ## (8) Force bisection as "ultima ratio regum"
        if x is None or not a.x <= x.x <= b.x :
            
            a , b = bisection ( self.__fun , a , b )
            ok = not samesign ( a.fx , b.fx ) 
            self.__counters [ 'bisection' ] += ok
            if ok : 
                xx = 0.5 * ( a.x + b.x  )  
                x  = Point ( xx , self.__fun ( xx ) )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x
                else                              : return x , a , b ## RETURN                
                print ( 'BISECTION/1' , x , a , b )

        """
        ## treat the end point intervals 
        delta = b.x - a.x
        dd    = delta / 10
        
        ax    = a.x + dd
        bx    = b.x - dd 
        if   x.x < ax  :            
            ## (a) agressive split at the left  edge
            nb = Point  ( ax , self.__fun ( ax ) )
            if 0 == nb.fx or iszero ( nb.fx ) : return nb , a , b
            if samesign ( nb.fx , b.fx      ) :
                b = nb 
                print ( 'AGRESSIVE SPLIT LEFT' , x , nb , a , b , b.x - a.x ) 
        elif bx < x.x :
            ## (b) agressive split at the right edge
            na = Point  ( bx , self.__fun ( bx ) )
            if 0 == na.fx or iszero ( na.fx ) : return na , a , b 
            if samesign ( na.fx , a.fx      ) :
                a = na 
                print ( 'AGRESSIVE SPLIT RIGHT' , x , na , a , b , b.x - a.x ) 
        
        delta = b.x - a.x
        if 3 * delta > old_len :
            
            a , b = bisection ( self.__fun , a , b )
            ok = not samesign ( a.fx , b.fx ) and a.x <= x.x <= b.x 
            self.__counters [ 'bisection' ] += ok
            if not ok :                
                xx = 0.5 * ( a.x + b.x  )  
                x  = Point ( xx , self.__fun ( xx ) )
                ## self.__inverse.append ( x  )
                ## self.__aitken .append ( xx )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x
                else                              : return x , a , b ## RETURN                
                print ( 'BISECTION/2' , x , a , b )
        """
                   
        return x , a , b 

# ===========================================================================================
## The main function for root-finding.
#
#  For each step it use sequence of well-defined "actions" The result of each "action"
#  is considered to be "successfull" if the obtained root approximation falls into the current
#  bracketing interval, otherwise it is considered as "failure"
# 
#  In the successfull case the obtained root it stored into list, and the bracketing
#  interval is correspodingly updated.

#  We start iterations from adding the current brackets into collection of roots
# 
#  The sequential list of actions is:
#  - 'Halley'     :                             if the second derivative is provided, make Halley's step
#  - 'Newton'     : if previous step is failed  and the first  derivative is provdied, make Newton step
#  - 'Steffensen' : if previous step is failed                                      , make Steffensen's step
#  - 'Inversenterpolation' : at this point we alwasy have at least two or more  approximations
#    for the root stored into the list (including initial brackets), so an inverse polynomial interpolation is performed: 
#     -- if              there are >=4 points, try with cubic
#     -- if it fails and there are >=3 points, try with parabolic
#     -- if it fails,                          try with secant   (and this one can't fail!)    
#  - 'Aitken'     : if there are >=3 approximations, try to use Aitken delta2 acceleratiion scheme  
#  - 'Bisection'  : if the updated bracketing interval is not shrinked significantly for the current step, a bisection step is forced 
#
#  After the step the convergency is checked:
#  - either the valeu of the  function   at the currect root  approximation is small enough
#  - or distance between two subsequent approximations is very small 
#  - or the updated bracketing interval is very small
#
#  Note, that:
#  - similar to Brent's and TOMS74 algortihms, due to the 'Bisection' action the algorithm converges ``almost always''.
#  - due to intermediate 'Aitken' step, a slow linear convergency (e.g. in vicinity of multiple root) is often drastically improved
# 
#  @see https://en.wikipedia.org/wiki/Halley%27s_method
#  @see https://en.wikipedia.org/wiki/Newton%27s_method
#  @see https://en.wikipedia.org/wiki/Steffensen%27s_method
#  @see https://en.wikipedia.org/wiki/Secant_method
#  @see https://en.wikipedia.org/wiki/False_position_method
#  @see https://en.wikipedia.org/wiki/Inverse_quadratic_interpolation
#  @see https://en.wikipedia.org/wiki/Aitken%27s_delta-squared_process
#  @see https://en.wikipedia.org/wiki/Bisection_method
#
#  The algorithm is largely inspired  by the Brent's method and newer TOMS748 algorithm from <code>Boost.math rootfinding</code> 
#  @see https://en.wikipedia.org/wiki/Brent%27s_method
#  @see http://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf
#  @see https://www.boost.org/doc/libs/1_67_0/boost/math/tools/toms748_solve.hpp
#  @see https://www.boost.org/doc/libs/release/libs/math/doc/html/root_finding.html
#  @see https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/roots_noderiv/bracket_solve.html
#  @see https://www.boost.org/doc/libs/release/boost/math/tools/toms748_solve.hpp
def find_root ( fun                 ,     ## the function 
                a                   ,     ## low-edge of bracketing 
                b                   , * , ## up-edge nbracketing               
                args        = ()    ,     ## additional positioal arguments for function (&derivatives) call  
                kwargs      = {}    ,     ## additional keyword   arguments for function (&derivatives) call   
                maxiter     = 100   ,
                xtol        = _xtol ,
                rtol        = _rtol ,
                full_output = False ,
                disp        = True  , 
                deriv1      = None  , 
                deriv2      = None  ) :
    """ The main function for root-finding

    Parameters
    ----------
    f           : the function
    a           : low edge of bracketing interval
    b           : high edge of bracketing interval
    args        : extra  positioanl argument for a function call
    kwargs      : extra  keyword  arguments for functon call 
    xtol        : absolute tolerance
    rtol        : relative tolerance
    full_output : If `full_output` is `False`, the root is returned.
                  If `full_output` is `True`, the return value is `(x, r)`, where `x` is the root, and `r` is a `RootResults` object.
    disp        : If `True`, raise `RuntimeError` if the algorithm didn’t converge.
    deriv1      : (optional) the first derivative 
    deriv2      : (optional) the second derivative

    Example
    -------
    
    >>> fun = ...
    >>> r = find_root ( fun , -10 , 10 )
    
    >>> fun  = ...
    >>> der1 = ...  ##  th first derivative 
    >>> r = find_root ( fun , -10 , 10 , deriv1 = der1 )

    >>> fun  = ...
    >>> der1 = ...  ##  the first derivative 
    >>> der2 = ...  ##  the second derivative 
    >>> r = find_root ( fun , -10 , 10 , deriv1 = der1 , deriv2 = der2 )


    For each step it use sequence of well-defined `actions'.
    The result of each ``action'' is considered to be `successfull',
    if the obtained root approximation falls into the current
    bracketing interval, otherwise it is considered as `failure'
    
    In the successfull case the obtained root it stored into list, and the bracketing
    interval is correspodingly updated.
    
    We start iterations from adding the current brackets into collection of roots
    
    The sequential list of `actions' is:
    - 'Halley'     :                             if the second derivative is provided, make Halley's step
    - 'Newton'     : if previous step is failed  and the first  derivative is provided, make Newton step
    - 'Steffensen' : if previous step is failed                                      , make Steffensen's step
    - 'Inversenterpolation' : at this point we alwasy have at least two or more  approximations
    for the root stored into the list (including initial brackets), so an inverse polynomial interpolation is performed: 
    -- if              there are >=4 points, try with cubic
    -- if it fails and there are >=3 points, try with parabolic
    -- if it fails,                          try with secant   (and this one can't fail!)    
    - 'Aitken'     : if there are >=3 approximations, try to use Aitken delta2 acceleratiion scheme  
    - 'Bisection'  : if the updated bracketing interval is not shrinked significantly for the current step, a bisection step is forced 
    
    After the step the convergency is checked and the process is converged  if:
    - either the value of the function at the currect root approximation is small enough
    - or the distance between two subsequent approximations is very small 
    - or the updated bracketing interval is very small
    
    Note, that:
    - similar to Brent's and TOMS748 algortihms, due to the 'Bisection' action the algorithm converges ``almost always''.
    - due to intermediate 'Aitken' step, a slow linear convergency (e.g. in vicinity of multiple root) is often drastically improved
    
    
    The algorithm is largely inspired  by the Brent's method and newer TOMS748 algorithm from Boost.math root finding library

    References
    ----------
    
    - see https://en.wikipedia.org/wiki/Halley%27s_method
    - see https://en.wikipedia.org/wiki/Newton%27s_method
    - see https://en.wikipedia.org/wiki/Steffensen%27s_method
    - see https://en.wikipedia.org/wiki/Secant_method
    - see https://en.wikipedia.org/wiki/False_position_method
    - see https://en.wikipedia.org/wiki/Inverse_quadratic_interpolation
    - see https://en.wikipedia.org/wiki/Aitken%27s_delta-squared_process
    - see https://en.wikipedia.org/wiki/Bisection_method
    
    - see https://en.wikipedia.org/wiki/Brent%27s_method
    - see http://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf
    - see https://www.boost.org/doc/libs/1_67_0/boost/math/tools/toms748_solve.hpp
    - see https://www.boost.org/doc/libs/release/libs/math/doc/html/root_finding.html
    - see https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/roots_noderiv/bracket_solve.html
    - see https://www.boost.org/doc/libs/release/boost/math/tools/toms748_solve.hpp

    """

    solver = RootFinder ( fun         = fun         ,
                          deriv1      = deriv1      ,
                          deriv2      = deriv2      ,
                          args        = args        , 
                          kwargs      = kwargs      ,
                          maxiter     = maxiter     ,
                          xtol        = xtol        ,
                          rtol        = rtol        ,
                          full_output = full_output , 
                          disp        = disp        )
    
    return solver.find ( a , b )

# =========================================================================
from scipy.optimize import brentq as scipy_brentq
# =========================================================================
## a tiny wrapper for `scipy.optimize.brentq`
def findroot_scipy ( fun                 , ## the function
                     a                   , ## low-edge  of bracheting interval 
                     b                   , ## high edge of bracketing interval 
                     args        = ()    ,     ## additional positioal arguments for function (&derivatives) call  
                     kwargs      = {}    ,     ## additional keyword   arguments for function (&derivatives) call   
                     maxiter     = 100   ,
                     xtol        = _xtol ,
                     rtol        = _rtol ,
                     full_output = False ,
                     disp        = True  , **other ) :
    
    if other :
        logger.warning ( "Extra arguments [%s] are ignored!" % ( ','.join ( k for k in other ) ) )

    ## wrap keyword arguemnts 
    if kwargs : the_fun = lambda x, *a : fun  ( x , *a, **kwargs )
    else      : the_fun = fun
    
    return scipy_brentq ( the_fun ,
                          a       ,
                          b       ,
                          args        = args        ,
                          maxiter     = maxiter     ,
                          xtol        = xtol        ,
                          rtol        = rtol        ,
                          full_output = full_output ,
                          disp        = disp        ) 


## as default, use scipy version 
findroot = findroot_scipy

# =============================================================================
## solve equation \f$ f(x)=C \f$
#  @code
#  fun = ...
#  C   = 
#  x   = sp_solve ( fun ,  0.0 , 1.0 , C ) 
#  @endcode 
def sp_solve ( fun           ,
               xmin          ,
               xmax          , * , 
               C = 0         ,
               solver = findroot , **other ) :
    """ Solve equation fun(x)=C
    >>> fun = ...
    >>> C   = ... 
    >>> x   = sp_solve ( fun , 0 , 1 , C  ) 
    """
    
    if iszero ( C ) : the_fun = fun 
    else            : the_fun = lambda x , *a , **kw : fun ( x , *a , **kw ) - C
    ## 
    return solver ( the_fun , xmin , xmax , **other )

# =============================================================================
## solve equation \f$ f(x)=C \f$
#  @code
#  fun = ...
#  C   = 
#  x   = sp_solve ( fun ,  0.0 , 1.0 , C = C ) 
#  @endcode 
solve = sp_solve 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if findroot is find_root : 
        logger.warning ("scipy.optimize.brentq is not available, use local `find_root'-replacement")

# =============================================================================
##                                                                      The END 
# =============================================================================
