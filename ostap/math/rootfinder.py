#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/rootfinder.py
#  Module with some useful utilities for root finding
#  - a kind of replacement for Brent's method if/when scipy is not accessible
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
    'secant'             , ## single step of secant/regular falsi/inverse linear
    'bisection'          , ## single step of bisection
    'aitken_delta2'      , ## aitken delta2 acceleration process   
)
# =============================================================================
from   ostap.core.ostap_types import num_types  
from   ostap.math.base        import samesign, iszero, isequal, isfinite, signum, Ostap     
from   ostap.utils.basic      import counted
from   ostap.stats.counters   import SE, EffCounter
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
    
    ## get function value if not specified 
    assert callable ( fun ) , "Invalid function `fun'"

    fx = fun ( x ) if fx is None else fx  
    if 0 == fx or iszero ( fx ) : return x                ## RETURN 

    if isequal ( x + fx , x )   : return None             ## RETURN
    
    gx = fun ( x + fx ) / fx - 1
    if 0 == gx or iszero ( gx ) : return None
    
    if not -1 < gx < 0 : return None 

    return x - fx / gx

# ============================================================================
## Single step of STFA method: Imporved regular falsi + steffenson-like
#  @see Xinyuan Wu, Zuhe Shen, Jianlin Xia, "An improved regula falsi method
#  with quadratic convergence of both diameter and point for enclosing
#  simple zeros of nonlinear equations", Applied Mathematics and Computation,
#  144, 2 2003
#  @see https://doi.org/10.1016/S0096-3003(02)00414-9
#  @see https://www.sciencedirect.com/science/article/pii/S0096300302004149
def STFA ( fun , a , b , x = None ) :
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
    d = b.x - a.x 
    if x is None or not a.x < x.x < b.x :        
        x = 0.5 * ( a.x + b.x )
        x = Point ( x , fun ( x ) )
        if 0 == x.fx or iszero ( x.fx ) : return x , a , b

    ## (1) regular falsi step
    c = ( a.x * b.fx - b.x * a.fx ) / ( b.fx - a.fx )
    c = Point ( c , fun ( c ) )
    
    if 0 == c.fx or iszero ( c.fx )     : return c , a , b
    elif samesign ( a.fx , c.fx )       : abar , bbar = c , b
    elif samesign ( b.fx , c.fx )       : abar , bbar = a , c
    else                                : return c , a , b

    ##
    if  c.fx == x.fx or isequal ( c.fx , x.fx ) : return None , abar , bbar

    mu   = ( b.x - a.x ) / ( b.fx - a.fx )
    cbar = x.x - mu * x.fx * x.fx / ( x.fx - c.fx )    
    cbar = Point ( cbar , fun ( cbar ) )

    if abar.x <= cbar.x <= bbar.x :
        x = cbar        
        if   samesign ( abar.fx , cbar.fx ) : a , b = cbar , bbar
        elif samesign ( bbar.fx , cbar.fx ) : a , b = abar , cbar 
    else : 
        x , a , b = c , a , b 

    return x , a , b 

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

    D = B * B - 4 * A * C
    if D < 0                       : return None , a , b ## RETURN 

    Q = B + signum ( B ) * math.sqrt ( D )
    p = c.x - 2 * C / Q
    
    if not a.x < p < b.x            : return None , a , b  ## RETURN 
    ## 
    p = Point ( p , fun ( p ) )

    ## zero is found ? 
    if 0 == p.fx or iszero ( p.fx ) : return p , a , b     ## RETURN 

    ## update the interval
    
    if isequal ( c.x , p.x ) :
        
        ## we have only three points here, new interval is trivial         
        a , b = ( c , b ) if samesign ( a.fx , c.fx ) else ( a , c )

    else :
        
        ## we have four points (three intervals), choose proper new interval  
        
        if   c.x < p.x : a , r1 , r2 , b = a , c , p , b
        elif c.x > p.x : a , r1 , r2 , b = a , p , c , b
        
        if   not samesign ( a.fx  , r1.fx ) : a , b = a  , r1
        elif not samesign ( r1.fx , r2.fx ) : a , b = r1 , r2 
        elif not samesign ( r2.fx ,  b.fx ) : a , b = r2 ,  b 

    dab = ( a.fx - b.fx ) / ( b.x - a.x )
    if 0.1 < abs ( dab ) < 10 :
        ## regular falsi 
        c = ( a.x * b.fx - b.x * a.fx ) / ( b.fx - a.fx )
    else :
        ## bisection 
        c = 0.5 * ( a.x + b.x )
        
    return Point ( c , fun ( c ) ) , a , b            ## RETURN 

# ============================================================================
## Single step of Muller's method
#  @see https://en.wikipedia.org/wiki/Muller%27s_method
def muller ( x0 , x1 , x2 , *other ) :
    """ Single step of Muller's method
    - see https://en.wikipedia.org/wiki/Muller%27s_method
    """
    points  = ( x0 , x1 , x2 ) + other
    x0 , x1 , x2 = points [ -3 : ]

    if x1.x == x0.x or isequal ( x1.x , x0.x ) : return None
    if x2.x == x1.x or isequal ( x2.x , x1.x ) : return None
    if x2.x == x0.x or isequal ( x2.x , x0.x ) : return None
    
    h0 = x1.x - x0.x
    h1 = x2.x - x1.x
    
    d0 = ( x1.fx - x0.fx ) / h0 
    d1 = ( x2.fx - x1.fx ) / h1

    ## parabola's coefficients: 
    a  = ( d1 - d0 ) / ( h1 + h0 )
    b  = a * h1 + d1
    c  = x2.fx
    
    ## Parabola's discriminant 
    D  = b * b - 4 * a * c
    
    if D < 0 : return None  ## RETURN
    
    Q = b + signum ( b ) * math.sqrt ( D )
    if iszero ( Q ) : return None 
    
    dx = -2 * c / Q

    return x2.x + dx

# =============================================================================
# secant method/regular falsi 
# - fallback to bisection if the regular falsi astimate is too close to the edges
def secant ( a , b ) :
    """ Secant/regular falsi 
    - fallback to bisection if the regular falsi astimate is too close to the edges
    """
    
    if   0 == a.fx or iszero ( a.fx ) : return a.x
    elif 0 == b.fx or iszero ( b.fx ) : return b.x
    ## 
    afx = abs ( a.fx )
    bfx = abs ( b.fx )
    dfx = 1.e-4 * ( afx + bfx )
    ##
    ## use bisection here...
    if   afx < dfx : return 0.5 * ( a.x + b.x )  ## ATTENTION! 
    elif bfx < dfx : return 0.5 * ( a.x + b.x )  ## ATTENTION! 
    ## regular falsi 
 
    ## bisection 
    if   isequal ( a.x  , b.x  ) : return 0.5 * ( a.x + b.x )
    elif isequal ( a.fx , b.fx ) : return 0.5 * ( a.x + b.x ) 
    
    ## regular falsi     
    return ( a.x * b.fx - b.x * a.fx ) / ( b.fx - a.fx )

regular_falsi = secant

# =============================================================================
## Inverse parabolic interpolation
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
    ## get last three points 
    c , b , a = ( ( a , b , c ) + other ) [-3:]
    
    if   a.x  == b.x  or isequal ( a.x  , b.x  ) : return secant ( a , c )
    elif a.x  == c.x  or isequal ( a.x  , c.x  ) : return secant ( a , b )
    elif b.x  == c.x  or isequal ( b.x  , c.x  ) : return secant ( a , b )    
    elif a.fx == b.fx or isequal ( a.fx , b.fx ) : return secant ( a , c )
    elif a.fx == c.fx or isequal ( a.fx , c.fx ) : return secant ( a , b )
    elif b.fx == b.fx or isequal ( b.fx , b.fx ) : return secant ( a , b )
    
    x0 , f0 = a.x , a.fx
    x1 , f1 = b.x , b.fx
    x2 , f2 = c.x , c.fx
    
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
    d , c , b , a = ( ( a , b , c , b ) + other ) [ -4 : ]
    
    if   a.x  == b.x  or isequal ( a.x  , b.x  ) : return inverse_parabolic ( a , c , d  )
    elif a.x  == c.x  or isequal ( a.x  , c.x  ) : return inverse_parabolic ( a , b , d  )
    elif a.x  == d.x  or isequal ( a.x  , d.x  ) : return inverse_parabolic ( a , b , c  )
    elif b.x  == c.x  or isequal ( b.x  , c.x  ) : return inverse_parabolic ( a , b , d  )
    elif b.x  == d.x  or isequal ( b.x  , d.x  ) : return inverse_parabolic ( a , b , c  )
    elif c.x  == d.x  or isequal ( c.x  , d.x  ) : return inverse_parabolic ( a , b , c  )
    elif a.fx == b.fx or isequal ( a.fx , b.fx ) : return inverse_parabolic ( a , c , d  )
    elif a.fx == c.fx or isequal ( a.fx , c.fx ) : return inverse_parabolic ( a , b , d  )
    elif a.fx == d.fx or isequal ( a.fx , d.fx ) : return inverse_parabolic ( a , b , c  )
    elif b.fx == c.fx or isequal ( b.fx , c.fx ) : return inverse_parabolic ( a , b , d  )
    elif b.fx == d.fx or isequal ( b.fx , d.fx ) : return inverse_parabolic ( a , b , c  )
    elif c.fx == d.fx or isequal ( c.fx , d.fx ) : return inverse_parabolic ( a , b , c  )
    
    x0 , f0 = a.x , a.fx
    x1 , f1 = b.x , b.fx
    x2 , f2 = c.x , c.fx
    x3 , f3 = d.x , d.fx
    
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
    
    if   0 == fc or iszero ( fc ) : return c , a ## RETURN 
    elif samesign ( a.fx , c.fx ) : return c , b ## RETURN
    elif samesign ( b.fx , c.fx ) : return a , c ## RETURN
    
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
    x0 , x1 , x2 = [ p.x for p in ( ( x0 , x1 , x2 ) + others ) [ -3 : ] ]
    
    dd = x2 + x0 - 2 * x1
    if not dd or iszero ( dd ) : return
    
    return ( x2 * x0 - x1 * x1 ) / dd 
    ## return x0 - ( x1 - x0 ) * ( x1 - x0 ) / dd 

# ========================================================================================
## A bit simpified version of single step of TOMS748
def toms748 ( fun , a , b , c = None ) :
    """ A bit simpified version of single step of TOMS748
    """
    if c is None or not a.x < c.x < b.x :
        c  = secant ( a , b )
        ok = not c is None and a.x < c < b.x
        if not ok : c = 0.5 * ( a.x + b.x )
        c = Point ( c , fun ( c ) )
        ## c  is a root ? 
        if 0 == c.fx or iszero ( c.fx )   : return c , a , b    ## RETURN 

    d = inverse_parabolic ( a , b , c )
    if d is None or not a.x < d < b.x     : return None , a , b ## RETURN
    d = Point ( d , fun ( d ) )
    
    ## d is a root ? 
    if 0 == d.fx or iszero ( d.fx )       : return d , a , b    ## RETURN 

    ## 
    e = inverse_cubic ( a , b , c , d )
    if e is None or not a.x < e < b.x :
        s = a if abs ( a.fx ) <= abs ( b.fx ) else b          
        e = inverse_parabolic ( s , c  , d )
        if e is None or not a.x < e < b.x  : return None , a , b ## RETURN
        
    e = Point ( e , fun ( e ) )
    
    ## e is a root ? 
    if 0 == e.fx or iszero ( e.fx )        : return e , a , b    ## RETURN 

    ## (MANUAL) re-bracketing 
    x1 , x2 , x3 = sorted ( ( c , d , e ) )
    if   not samesign ( a .fx , x1.fx ) : a , b = a  , x1
    elif not samesign ( x1.fx , x2.fx ) : a , b = x1 , x2
    elif not samesign ( x2.fx , x3.fx ) : a , b = x2 , x3
    elif not samesign ( x3.fx ,  b.fx ) : a , b = x3 ,  b 

    if   a.x <= e.x <= b.x : return e , a , b
    elif a.x <= d.x <= b.x : return d , a , b
    elif a.x <= c.x <= b.x : return c , a , b
    
    f = secant ( a , b )
    f = Point  ( f , fun ( f ) )
    
    return f , a , b 

# ============================================================================================
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
                  'counters'          ,
                  'shrinks'           )
    
    CONVERGED   = 'converged'
    SIGNERR     = 'sign error'
    CONVERR     = 'convergence error'
    ZEROFOUND   = 'converged/zero found'
    STABLEPOINT = 'converged/stable point'
    BRACKETED   = 'converged/bracketed'
    
    flag_map  = {  0 : CONVERGED , -1 : SIGNERR     , -2 : CONVERR   ,
                   1 : ZEROFOUND ,  2 : STABLEPOINT ,  3 : BRACKETED }
    
    def __init__  ( self                     ,
                    root                     ,
                    iterations               ,
                    function_calls           , * ,
                    derivative1_calls = 0    ,
                    derivative2_calls = 0    ,                                      
                    flag              = 0    , 
                    value             = None ,
                    bracket           = ()   , 
                    counters          = {}   ,
                    shrinks           = {}   ) : 
        
        self.root              = root
        self.iterations        = iterations
        self.function_calls    = function_calls
        self.derivative1_calls = derivative1_calls
        self.derivative2_calls = derivative2_calls
        self.flag              = flag
        self.status            = self.flag_map.get ( flag , '*UNKNOWN*' )
        self.value             = value
        self.counters          = counters
        self.shrinks           = shrinks
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
                   
        row    = 'Status' , attention ( self.status ) if self.flag < 0 else allright ( self.status ) 
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
                v  = '[%s,%s]' % ( fa , fb)
                row = 'Bracket' , v , '10^%+d' % expo
            else :
                fa = fmtv % ( a )
                fb = fmtv % ( b )
                v  = '[%s,%s]' % ( fa , fb )
                row = 'Bracket' , v , '' 
            rows.append ( row )
            
            c = b - a 
            fc , ec  = pretty_float ( c , precision = precision , width = width , with_sign = True )
            if ec : row = 'Lenght' , fc , '10^%+d' % ec 
            else  : row = 'Length' , fc  
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
            row = 'Methods:' , ' success / failure'
            rows.append ( row ) 
            for key,cnt  in self.counters.items() :
                success = cnt.success
                failure = cnt.failure 
                row     = key , '%7d /%- 7d' % ( success , failure )
                rows.append ( row )
                
        if self.shrinks :
            row = 'Shrinks:' , '' 
            rows.append ( row ) 
            for key,cnt  in self.shrinks.items() :
                smean       = cnt.mean() 
                value, expo = smean.pretty_print ( precision   = 2  ,
                                                   width       = 4  ,
                                                   parentheses = False     )
                row     = key , value , '10^%+d' % expo if expo else '' 
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
                assert callable( deriv2 ) , "Invalid second derivative `deriv2'"
                if has_args : fderiv2 = lambda x : deriv2 ( x , *self.args , **self.kwargs )
                else        : fderiv2 = deriv2
     
        self.__fun     = counted ( ffun    ) if self.__full_output             else ffun  
        self.__deriv1  = counted ( fderiv1 ) if self.__full_output and fderiv1 else fderiv1 
        self.__deriv2  = counted ( fderiv2 ) if self.__full_output and fderiv2 else fderiv2

        self.__xtol    = max ( xtol , _xtol )   
        self.__rtol    = max ( rtol , _rtol )

        self.__iteratuon = 0 
        self.__maxiter   = maxiter if 1 <= maxiter else 500
        
        self.__disp      = True    if disp        else False
        
        ## keep several previous approximations for inverse polynomial, Muller or Aitken
        self.__roots   = [] 

        self.__counters = collections.defaultdict ( EffCounter )
        self.__shrinks  = collections.defaultdict ( SE )
        
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

    @property
    def full_output ( self ) :
        """`full_output` : provide detailed information
        """
        return self.__full_output
    
    # =========================================================================
    ## Return result/root with optional report 
    def __result ( self       ,
                   root      , * ,
                   value     = None ,
                   bracket   = ()   , 
                   flag      = 0    ) :
        """ Return result/root with optional report 
        """
        if not self.__full_output : return root

        return root , RootResults ( root ,
                                    iterations        = self.__iteration ,
                                    function_calls    = self.__fun.calls ,
                                    derivative1_calls = self.__deriv1.calls if self.__deriv1 else 0 ,
                                    derivative2_calls = self.__deriv2.calls if self.__deriv2 else 0 ,
                                    value             = value           , 
                                    flag              = flag            ,
                                    bracket           = bracket         , 
                                    counters          = self.__counters , 
                                    shrinks           = self.__shrinks  )
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

        self.__counters .clear()        
        self.__shrinks  .clear()
        self.__iteration  = 0 

    
        ## swap them if needed 
        a , b = ( a , b ) if a < b else ( b , a )

        ## check the root at left edge 
        fa = self.__fun ( a )
        if 0 == fa or iszero ( fa ) :
            return self.__result ( a , value = fa , bracket = ( a , b ) , flag = 1 )

        ## check the root at right edge 
        fb = self.__fun ( b )
        if 0 == fb or iszero ( fb ) :
            return self.__result ( b , value = fb , bracket = ( a , b ) , flag = 1 )

        ## correct bracketing interval ?
        if samesign ( fa ,  fb ) :
            if self.__disp :
                raise RuntimeError ('RootFinder.solve: invalid bracketing interval (%s,%s)/(%s,%s)' % ( a , fa , b , fb ) )
            else           :
                return self.__result ( 0.5 * ( a + b ) , flag = -1   , bracket = ( a , b ) )  ## RETURN
            
        a = Point ( a ,  fa )
        b = Point ( b ,  fb )

        ## initialize 
        self.__roots = [ a , b ] 
        
        ## Guess is not specified or invalid

        assert guess is None or isinstance ( guess , num_types ) , \
            "Invaild `guess` type %s" % typename ( guess )

        ## check the provided guess 
        if guess is None or not a.x <= guess <= b.x :
            
            ## use the secant as the 1st approximation 
            guess = secant ( a , b )
            ok    = not guess is None and isfinite ( guess ) and a.x <= guess <= b.x
            if self.__full_output : self.__counters [ 'secant' ] += ok
            if not ok :
                guess = 0.5 * ( a.x + b.x )
                if self.__full_output : self.__counters [ 'bisection' ] += 1

        ## calcualate fun(guess) 
        guess = Point ( guess , self.__fun ( guess ) )

        ## check the function value and shrink the interval, if/when  possible 
        if 0 == guess.fx or iszero ( guess.fx  ) :
            return self.__result ( guess.x , value = guess.fx , bracket = ( a.x , b.x ) , flag = 1 ) ##  RETURN
        elif samesign ( a.fx , guess.fx ) : a , b = guess , b 
        elif samesign ( b.fx , guess.fx ) : a , b = a     , guess
        else :
            return self.__result ( guess.x , value = guess.fx , bracket = ( a.x , b.x ) , flag = 1 ) ##  RETURN
        
        ## initialize Aitken and Inverse Interpolation structures 
        self.__roots     = [] 
        
        ## 
        self.__iteration = 0
        self.__iteration = 0

        zeroes = [ guess ]
        while self.__iteration <= self.__maxiter :

            self.__iteration += 1 

            ## get new x and new bracketing interval 
            x , a , b = self.__make_step ( a , b , guess )
        
            ## zero is found ?
            if not x.fx or iszero ( x.fx ) :
                return self.__result ( x.x , value = x.fx , bracket = ( a.x , b.x ) , flag = 1  )
            
            ## no difference with respect to the previous estimate 
            if 2 < self.__iteration :
                dx = abs ( x.x - guess.x )
                if dx < 0.5 * max ( self.__xtol , self.__rtol * dx )  :
                   return self.__result ( x.x , value = x.fx , bracket = ( a.x , b.x ) , flag = 2 )

            ## Bracketing interval is already very small 
            dab = abs ( b.x - a.x )
            if dab < 0.5 * max ( self.__xtol , self.__rtol * dab )  :
                return self.__result ( x.x , value = x.fx , bracket = ( a.x , b.x ) , flag = 3 )
            
            ## next guess: 
            guess = x 
            
        message = 'RootFinder.solve: no convergency for %d iterations' % self.__iteration 
        if self.__disp : raise RuntimeError ( message )
        else           : logger.warning     ( message ) 
            
        return self.__result ( x.x , flag = -2 , value = x.fx , bracket = ( a.x , b.x ) ) 

    # =========================================================================
    def __make_step ( self , a , b , x = None ) :
        
        ## do not keep too many approximations

        self.__roots = self.__roots [ -5 : ]
        
        ## the lenght of the interval brackeing interval 
        the_len = b.x - a.x

        ## if guess is not specified or invalid 
        if x is None or not a.x <= x.x <= b.x : 
            ## (0.1) make a secant step
            xx = secant ( a , b )
            ok = not xx is None and isfinite ( xx ) and a.x <= xx <= b.x
            if self.__full_output : self.__counters [ 'secant' ] += ok
            if not ok :
                ## (0.2) 
                xx =  0.5 * ( a.x + b.x )            
                if self.__full_output : self.__counters [ 'bisection' ] += ok
            x = Point ( xx , self.__fun ( xx ) )

        updated = False 
        # =====================================================================
        ## (1&2) try Halley and/or Newton methods
        if self.__deriv1 :            
            ## 
            xx = halley_newton ( self.__fun    ,
                                 x.x           ,
                                 self.__deriv1 , 
                                 deriv2 = self.__deriv2 , fx = x.fx )
            ## OK ?
            ok = not xx is None and isfinite ( xx ) and a.x <= xx <= b.x
            if self.__full_output : 
                if self.__deriv2 : self.__counters [ 'halley' ] += ok  
                else             : self.__counters [ 'newton' ] += ok  
            if ok :
                cur_len = b.x - a.x 
                x = Point ( xx , self.__fun ( xx ) ) 
                if not x in self.__roots : self.__roots.append ( x  )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN                 
                updated = True
                if self.__full_output : 
                    shrink  = ( b.x - a.x ) / cur_len
                    if self.__deriv2 : self.__shrinks [ 'halley' ] += shrink
                    else             : self.__shrinks [ 'newton' ] += shrink
                
        # =====================================================================
        ## (3)  SFTA: "falsi-setffensen"
        # =====================================================================
        if not updated or not a.x <= x.x <= b.x :
            cur_len = b.x - a.x 
            xx , a , b = STFA ( self.__fun , a , b , x )
            ok = not xx is None and a.x <= xx.x <= b.x
            if self.__full_output : self.__counters [ 'STFA' ] += ok              
            if ok :
                x = xx
                if not x in self.__roots : self.__roots.append ( x  )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN                 
                updated = True
                if self.__full_output : 
                    shrink  = ( b.x - a.x ) / cur_len
                    self.__shrinks [ 'STFA' ] += shrink
        
        # =====================================================================
        ## (4) try  Steffensen's method only if the  interval is already shrinked
        if not updated                            and \
           a.x <= x.x <= b.x                      and \
           a.x <= x.x + x.fx <= b.x               and \
           abs ( x.fx ) < 0.1 * abs ( b.x - a.x ) and \
           not isequal ( x.x + x.fx , x.x ) :         
            ##
            xx = steffensen ( self.__fun , x.x , fx = x.fx )
            ## OK ?
            ok = not xx is None and isfinite ( xx ) and a.x <= xx <= b.x
            if self.__full_output : self.__counters [ 'steffensen' ] += ok
            if ok :
                cur_len = b.x - a.x 
                x = Point ( xx , self.__fun ( xx ) ) 
                if not x in self.__roots : self.__roots.append ( x  )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN                    
                updated = True
                if self.__full_output : 
                    shrink  = ( b.x - a.x ) / cur_len 
                    self.__shrinks [ 'Steffensen' ] += shrink 
                
        # =====================================================================
        ## (5) Try TOMS748  method
        if True : ## not updated or not a.x <= x.x <= b.x :
            cur_len = b.x - a.x 
            xx , a , b = toms748 ( self.__fun , a , b , x )
            ok = not xx is None and a.x <= xx.x <= b.x
            if self.__full_output : self.__counters [ 'TOMS748' ] += ok 
            if ok :
                x = xx 
                if not x in self.__roots : self.__roots.append ( x  )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN                    
                updated = True                
                if self.__full_output :                 
                    shrink  = ( b.x - a.x ) / cur_len
                    self.__shrinks [ 'TOMS748' ] += shrink                 
                
        # =====================================================================
        ## (6) Try RBP method
        if not updated or not a.x <= x.x <= b.x :
            cur_len = b.x - a.x 
            xx , a , b  = RBP ( self.__fun , a , b , x )
            ok = not xx is None and a.x <= xx.x <= b.x
            if self.__full_output : self.__counters [ 'RBP' ] += ok
            if ok :
                x = xx
                ## if not x in self.__roots : self.__roots.append ( x  )
                if   0 == x.fx or iszero ( x.fx )     : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )         : a , b = x , b  
                elif samesign ( b.fx , x.fx )         : a , b = a , x  
                else                                  : return x , a , b ## RETURN
                updated = True
                if self.__full_output : 
                    shrink  = ( b.x - a.x ) / cur_len 
                    self.__shrinks [ 'RBP' ] += shrink 

        # =====================================================================
        ## (7) Try Muller's method
        if False : ## 3 <= len ( self.__roots )  : ## and ( not updated or not a.x <= x.x <= b.x ) :
            xx = muller ( *self.__roots )
            ok = not xx is None and isfinite ( xx ) and a.x <= xx <= b.x
            if self.__full_output : self.__counters [ 'muller' ] += ok
            if ok :
                cur_len = b.x - a.x 
                x = Point ( xx , self.__fun ( xx ) ) 
                if not x in self.__roots : self.__roots.append ( x  )
                if   0 == x.fx or iszero ( x.fx )     : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )         : a , b = x , b  
                elif samesign ( b.fx , x.fx )         : a , b = a , x  
                else                                  : return x , a , b ## RETURN
                updated = True
                if self.__full_output :                     
                    shrink  = ( b.x - a.x ) / cur_len 
                    self.__shrinks [ 'muller' ] += shrink

        # =====================================================================
        ## (8) secant is "almost never" fails 
        if not updated or not a.x <= x.x <= b.x :
            xx = secant ( a , b )
            ## OK ?
            ok = not xx is None and isfinite ( xx ) and a.x <= xx <= b.x 
            if self.__full_output : self.__counters [ 'secant/R' ] += ok
            if ok :
                cur_len = b.x - a.x 
                x = Point ( xx , self.__fun ( xx ) ) 
                ## if not x in self.__roots : self.__roots.append ( x  )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN                 
                updated = True
                if self.__full_output :                     
                    shrink  = ( b.x - a.x ) / cur_len 
                    self.__shrinks [ 'secant' ] += shrink 

        # =====================================================================
        ## (9)  check inverse cubic:        
        if False : ## 4 <= len ( self.__roots ) :
            xx = inverse_cubic ( *self.__roots )
            ok = not xx is None and isfinite ( xx ) and a.x <= xx <= b.x
            if self.__full_output : self.__counters [ 'inverse_cubic' ] += ok
            if ok :
                cur_len = b.x - a.x 
                x = Point ( xx , self.__fun ( xx ) ) 
                if not x in self.__roots : self.__roots.append ( x  )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN                 
                updated = True
                if self.__full_output :                     
                    shrink  = ( b.x - a.x ) / cur_len 
                    self.__shrinks [ 'invese_cubic' ] += shrink                     
            
        # =====================================================================
        ## (10) Try Aitken-delta2 scheme if we have enough approximations
        if False : ## 3 <= len ( self.__roots ) :
            xx = aitken_delta2 ( *self.__roots )
            ok = not xx is None and isfinite ( xx ) and a.x <= xx <= b.x
            ## accept Aitken only in case it really improves the root estimate 
            if ok :
                fxx = self.__fun ( xx )
                ok  = not x is None and abs ( fxx ) <= abs ( x.fx )
            if self.__full_output : self.__counters [ 'aitken' ] += ok
            if ok :
                cur_len = b.x - a.x 
                x = Point ( xx , fxx ) 
                ## if not x in self.__roots : self.__roots.append ( x  )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x
                else                              : return x , a , b ## RETURN                
                updated = True
                if self.__full_output :                     
                    shrink  = ( b.x - a.x ) / cur_len 
                    self.__shrinks [ 'aitken' ] += shrink 

        # =====================================================================        
        ## (11) Force bisection as "the ultima ratio regum"
        new_len = b.x - a.x
        if not updated or not a.x <= x.x <= b.x or 3 * new_len > the_len :
            cur_len = b.x - a.x 
            a , b = bisection ( self.__fun , a , b )
            ok = not samesign ( a.fx , b.fx ) 
            if self.__full_output : self.__counters [ 'bisection' ] += ok
            if ok : 
                if not a.x <= x.x <= b.x :
                    x = secant ( a , b )  
                    x = Point  ( x , self.__fun ( x ) )
                ## 
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x
                else                              : return x , a , b ## RETURN
                ## 
                if self.__full_output :                     
                    shrink  = ( b.x - a.x ) / cur_len 
                    self.__shrinks [ 'bisection' ] += shrink 

        ## final... 
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
                guess       = None  , 
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

    return solver.find ( a , b , guess = guess )

# =========================================================================
from scipy.optimize import brentq as scipy_brentq
# =========================================================================
## a tiny wrapper for `scipy.optimize.brentq`
def findroot_scipy ( fun                 ,     ## the function
                     a                   ,     ## low-edge  of bracheting interval 
                     b                   , * , ## high edge of bracketing interval
                     args        = ()    ,     ## additional positioal arguments for function (&derivatives) call  
                     kwargs      = {}    ,     ## additional keyword   arguments for function (&derivatives) call   
                     maxiter     = 100   ,
                     xtol        = _xtol ,
                     rtol        = _rtol ,
                     full_output = False ,
                     disp        = True  , **other ) :
    
    if other : logger.warning ( "Extra arguments [%s] are ignored!" % ( ','.join ( k for k in other ) ) )
    
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

# =============================================================================
findroot_brent = findroot_scipy
findroot_ostap = find_root 
# =============================================================================
## as default, use scipy version 
findroot       = findroot_scipy


## a tiny wrapper for `scipy.optimize.brentq`
def findroot_ostap2 ( fun                 ,     ## the function
                      a                   ,     ## low-edge  of bracheting interval 
                      b                   , * , ## high edge of bracketing interval
                      args        = ()    ,     ## additional positioal arguments for function (&derivatives) call  
                      kwargs      = {}    ,     ## additional keyword   arguments for function (&derivatives) call   
                      maxiter     = 100   ,
                      xtol        = _xtol ,
                      rtol        = _rtol ,
                      full_output = False ,
                      disp        = True  , 
                      deriv1      = None  , 
                      deriv2      = None  , **other ) :
    
    if other : logger.warning ( "Extra arguments [%s] are ignored!" % ( ','.join ( k for k in other ) ) )
    
    ## wrap keyword arguemnts 
    if args or kwargs : the_fun = lambda x : fun  ( x , *args , **kwargs )
    else              : the_fun = fun
                        
    ##  get c++ RootFinder 
    rf = Ostap.Math.RootFinder ( maxiter , -1.0 , xtol , rtol )

    af = float ( a ) 
    bf = float ( b )
     
    import ctypes 
    aa = ctypes.c_double ( min ( af , bf ) ) 
    bb = ctypes.c_double ( max ( af , bf ) ) 
    
    ## use half sum as the first approximaiton to the root 
    rr = ctypes.c_double ( 0.5  * ( af  + bf ) ) 
    
    if deriv1 and callable ( deriv1 ) :
        
        if args or kwargs : the_der1 = lambda x : deriv1 ( x , *args , ** kwargs )
        else              : the_der1 = deriv1
        
        if deriv2 and callable ( deriv2 ) :
            
            if args or kwargs : the_der2 = lambda x : deriv2 ( x , *args , ** kwargs )
            else              : the_der2 = deriv1  
             
            sc = rf.root ( the_fun  , rr , aa , bb , the_der1 , the_der2 )
        else :
            sc = rf.root ( the_fun , rr , aa , bb , the_der1 )
    else : 
        sc = rf.root ( the_fun , rr , aa , bb )

    r = rr.value
    a = aa.value 
    b = bb.value 
    
    if not sc.isSuccess() :
        message = 'RootFinder.root: Status code %s for %d fun-calls' % ( sc , rf.ncalls() ) 
        if    disp : raise RuntimeError ( message )
        else       : logger.warning     ( message ) 
      
    if not full_output : return r 
    
    return r , RootResults ( r ,
                             iterations     = rf.ncalls ()  ,
                             function_calls = rf.ncalls ()  ,
                             bracket        = ( a , b )     ) 


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

# =============================================================================
##                                                                      The END 
# =============================================================================
