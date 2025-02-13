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
    'solve'              , ## solev f(x)=C equation 
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
from   ostap.math.base            import isequal, iszero , samesign
import sys, collections, warnings 
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.rootfinder' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
_xtol = 1.e-14
_rtol = sys.float_info.epsilon*2 
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
def halley_newton ( fun            ,    ## the function 
                    x              ,    ## x 
                    deriv1         ,    ## the first derivative 
                    deriv2  = None ,    ## the second derivative
                    fx      = None ,    ## value of fun(x)
                    args    = ()   ) :  ## additional arguments for function calls 
    """Single step of Newton/Halley's algorithm

    Parameters
    -----------
    
    fun    : the function 
    x      : the initial guess for the root positon
    deriv  : the first derivative  (for Newton and Halley's method) 
    deriv2 : the second derivative (for Halley's method)
    fx     : (optional) the value of funnction at the guess point
    args   : (optional) parameters for function/derivative calls 

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
    
    ## newton corrections
    d1 = float ( deriv1 ( x , *args ) ) 
    fx = float ( fun    ( x , *args ) ) if fx is None else fx 
    
    if d1  : rn = fx / d1
    else   : return None                               ## error here! 
    
    ## make corrections: Halley's  steps
    if deriv2 :  
        d2 = float ( deriv2 ( x , *args ) )
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
def steffensen ( fun  ,  x  , fx = None , args = () ) :
    """Single step of Steffensen's method

    Parameters
    ----------
    
    fun  : the function itself
    x    : the initial guess for the root
    fx   : (optional) function value at `x`, `fx=fun(x)`
    args : (optional) parameters for function calls

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
    
    if fx is None : fx = float ( fun ( x , *args ) ) ## reuse if already calculated
    if fx : 
        gx = ( fun ( x + fx , *args ) - fx ) / fx
        if gx : return x - fx / gx


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
def inverse_linear ( a , b ) :
    """Make a linear interpolation, aka ''secant'', aka ''regular falsi''

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
    
    x0 , f0 = a.x , a.fx
    x1 , f1 = b.x , b.fx

    if f0 == f1 or isequal ( f0 , f1 ) : return None

    return  ( x0 * f1  - x1 * f0 ) / ( f1 - f0 ) 


secant        = inverse_linear 
regular_falsi = inverse_linear



# =============================================================================
## Inverse parabolic interppolaiton
#  @see https://en.wikipedia.org/wiki/Inverse_quadratic_interpolation
#  @code
#  xa , xb , xc = ...
#  fun          = ...
#  r = inverse_parabolic ( Point ( xa , fun ( xa ) ) ,
#                          Point ( xb , fun ( xb ) ) ,
#                          Point ( xc , fun ( xc ) ) )
#  @endcode
def inverse_parabolic ( a , b , c ) :
    """Inverse parabolic interpolation via the points:
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

    
    x0 , f0 = a.x , a.fx
    x1 , f1 = b.x , b.fx
    x2 , f2 = c.x , c.fx
    
    if f0 == f1 or isequal ( f0 , f1 ) : return inverse_linear ( a , c )  
    if f0 == f2 or isequal ( f0 , f2 ) : return inverse_linear ( a , b ) 
    if f1 == f2 or isequal ( f1 , f2 ) : return inverse_linear ( a , b ) 
    
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
def inverse_cubic ( a , b , c , d ) :
    """Inverse cubic interpolaton via the last four approximations to the root

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

    x0 , f0 = a.x , a.fx
    x1 , f1 = b.x , b.fx
    x2 , f2 = c.x , c.fx
    x3 , f3 = d.x , d.fx

    ## switch to inverse parabolic if some function values coincide 
    if f0 == f1 or isequal ( f0 , f1 ) : return inverse_parabolic ( a , c , d ) 
    if f0 == f2 or isequal ( f0 , f2 ) : return inverse_parabolic ( a , b , d )
    if f0 == f3 or isequal ( f0 , f3 ) : return inverse_parabolic ( a , b , c )
    if f1 == f2 or isequal ( f1 , f2 ) : return inverse_parabolic ( a , b , d ) 
    if f1 == f3 or isequal ( f1 , f3 ) : return inverse_parabolic ( a , b , c ) 
    if f2 == f3 or isequal ( f2 , f3 ) : return inverse_parabolic ( a , b , c ) 

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
    """Make an inverse polynomial interpolation
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
    
    if   2 <= len ( other ) : return inverse_cubic     ( a , b , *other[:2] )
    elif 1 <= len ( other ) : return inverse_parabolic ( a , b , *other[:1] )
    
    return inverse_linear ( a , b )

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
def bisection ( fun , a  , b , args = () ) :
    """Trivial bisection method

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
    
    a , b = ( a, b ) if a.x < b.x else ( b , a )

    c = 0.5 * ( a.x + b.x )
    c = Point ( c , fun ( c , *args ) )
    
    if 0 == c.fx or iszero ( c.fx ) : return c , c ## RETURN 
    elif samesign ( a.fx , c.fx )   : return c , b ## RETURN
    elif samesign ( b.fx , c.fx )   : return a , c ## RETURN
    
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
def aitken_delta2 ( xl , xl1 , xl2 , *others ) :
    """Try to find next approximation to the root using Aitken delta-squared process.
    
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
    
    dd =  ( xl -  xl1 ) - ( xl1 - xl2 )
    
    if not dd : return None

    return xl - ( xl - xl1 ) **2 / dd 
    

# =============================================================================
## create 'counted' function to know number of function calls
#  @code
#
#  fun = ...
#  func = counted ( fun ) ## use as function
#
#  @endcode
def counted ( f ):
    """create 'counted' function to knon number of function calls

    Example
    -------
    
    >>> fun = ...
    >>> func = counted ( fun ) ## use as function

    """
    def wrapped ( *args, **kwargs ):
        wrapped.calls += 1
        return f( *args , **kwargs )
    wrapped.calls = 0
    return wrapped

# ========================================================================================
## @class RootResults
#  Helper class to keep results of root-finding preocedure
# 
#  It is very similar (almost clone) of corresponding class
#  <code>RootResults</code> from <code>scipy.optimize.zeros</code>
class RootResults ( object ):

    __slots__ = ( 'converged'            ,
                  'flag'                 ,
                  'function_calls'       ,
                  'iterations'           , 
                  'root'                 ,
                  ##
                  'derivative1_calls'    ,
                  'derivative2_calls'    ,
                  ##
                  'HalleyNewton'         ,                  
                  'Steffensen'           ,
                  'Secant'               ,
                  'InverseInterpolation' ,
                  'Aitken'               ,
                  'Bisection'            )

    CONVERGED = 'converged'
    SIGNERR   = 'sign error'
    CONVERR   = 'convergence error'
    flag_map  = { 0 : CONVERGED, -1 : SIGNERR , -2 : CONVERR }
    
    def __init__  ( self , root , iterations , function_calls , flag = 0  ):
        
        self.root             = root
        self.iterations       = iterations
        self.function_calls   = function_calls
        self.converged = flag == 0
        try:
            self.flag = self.flag_map[flag]
        except KeyError:
            self.flag = 'unknown error %d' % (flag,)
            
    def __repr__(self):
        
        m   = max ( map ( len , self.__slots__ ) ) + 1
        return '\n'.join ( [ a.rjust ( m ) + ': ' + repr ( getattr ( self , a ) )
                             for a in self.__slots__ if hasattr ( self , a ) ] )
    
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
#  - 'Hailey'     :                             if the second derivative is provided, make Halley's step
#  - 'Newton'     : if previous step is failed  and the first  derivative is provdied, make Newton step
#  - 'Steffensen' : if previous step is failed                                      , make Steffensen's step
#  - 'Secant'     : if previous step is failed (it should not, but numerical failures occurs)
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
class RootFinder(object) :
    """ The main root-finding engine.
    
    For each step it use sequence of well-defined ``actions''.
    The result of each ``action'' is considered to be ``successfull'',
    if the obtained root approximation falls into the current
    bracketing interval, otherwise it is considered as ``failure''
    
    In the successfull case the obtained root it stored into list, and the bracketing
    interval is correspodingly updated.
    
    We start iterations from adding the current brackets into collection of roots
    
    The sequential list of ``actions'' is:
    -  'Hailey'     :                             if the second derivative is provided, make Halley's step
    -  'Newton'     : if previous step is failed  and the first  derivative is provdied, make Newton step
    -  'Steffensen' : if previous step is failed                                      , make Steffensen's step
    -  'Secant'     : if previous step is failed (it should not, but numerical failures occurs)
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
                   fun                 ,
                   deriv1      = None  , 
                   deriv2      = None  ,
                   maxiter     = 200   ,
                   args        = ()    , 
                   xtol        = _xtol ,
                   rtol        = _xtol ,
                   full_output = False ,
                   disp        = True  ) :
        """Create the root-finder

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
        
        self.__full_output = full_output

        if fun    : ffun    = lambda  x , *a : float ( fun    ( x , *a ) )
        else      : ffun    = None            
        if deriv1 : fderiv1 = lambda  x , *a : float ( deriv1 ( x , *a ) )
        else      : fderiv1 = None 
        if deriv2 : fderiv2 = lambda  x , *a : float ( deriv2 ( x , *a ) ) 
        else      : fderiv2 = None 

        self.__fun     = counted ( ffun    ) if full_output             else ffun
        self.__deriv1  = counted ( fderiv1 ) if full_output and fderiv1 else fderiv1 
        self.__deriv2  = counted ( fderiv2 ) if full_output and fderiv2 else fderiv2
        
        self.__xtol    = xtol    if 0 < xtol    else _xtol 
        self.__rtol    = rtol    if 0 < rtol    else _rtol 
        self.__maxiter = maxiter if 0 < maxiter else   200
        self.__disp    = True    if disp        else False
        
        self.__try_halley         = 0
        self.__try_steffensen     = 0
        self.__try_secant         = 0
        self.__try_inverse        = 0
        self.__try_aitken         = 0
        self.__try_bisection      = 0
        
        self.__success_halley     = 0
        self.__success_steffensen = 0
        self.__success_secant     = 0
        self.__success_inverse    = 0
        self.__success_aitken     = 0
        self.__success_bisection  = 0
        self.__args               = args
        
        ## keep the last three roots for Aitken scheme 
        self.__aitken  = []
        
        ## keep the last three/four approximations for inverse cubic/parabolic interpolation
        self.__inverse = [] 

    # ========================================================================= 
    def __make_report ( self , root , i = 0 , flag = 0 ) :
        """Prepare a report
        """

        result = RootResults ( root , i , self.__fun.calls )
        
        if self.__deriv1 and hasattr ( self.__deriv1 , 'calls' ) :
            result.derivative1_calls = self.__deriv1.calls
            
        if self.__deriv2 and hasattr ( self.__deriv2 , 'calls' ) :
            result.derivative2_calls = self.__deriv2.calls
        
        if self.__try_halley     :
            result.HalleyNewton         = self.__try_halley     , self.__success_halley
        if self.__try_steffensen :
            result.Steffensen           = self.__try_steffensen , self.__success_steffensen
        if self.__try_secant     :
            result.Secant               = self.__try_secant     , self.__success_secant
        if self.__try_inverse    :
            result.InverseInterpolation = self.__try_inverse    , self.__success_inverse
        if self.__try_aitken     :
            result.Aitken               = self.__try_aitken     , self.__success_aitken 
        if self.__try_bisection  :
            result.Bisection            = self.__try_bisection  , self.__success_bisection

        return root , result

    # =========================================================================
    ## solve the function and find the root 
    def find ( self , a , b , guess = None ) :
        """Find a root
        
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
        fa = self.__fun ( a , *self.__args )
        if 0 == fa or iszero ( fa ) :
            if self.__full_output : return self.__report ( a ) ## RETURN
            return a                                           ## RETURN
        
        ## check the root at right edge 
        fb = self.__fun ( b , *self.__args )
        if 0 == fb or iszero ( fb ) :
            if self.__full_output : return self.__report ( b ) ## RETURN
            return b                                           ## RETURN

        ## correct bracketing interval ?
        if samesign ( fa ,  fb ) :
            if self.__disp :
                raise RuntimeError ('RootFinder.solve: invalid bracketing interval (%s,%s)/(%s,%s)' % ( a , fa , b , fb ) )
            else           :
                return self.__report ( 0.5 * ( a  + b ) , 0 , -1 )  ## RETURN
            
        a = Point ( a ,  fa )
        b = Point ( b ,  fb )
        
        ## if guess if not specified or outside the interval
        if guess is None or not a.x <= guess <= b.x :
            
            ## use the secant as the approximation 
            guess  = ( a.x * b.fx - b.x * a.fx ) / ( b.fx - a.fx  )
            
            l =  abs ( b.x - a.x )

            if   abs ( a.x - guess ) < 0.01 * l : guess = a.x + 0.10 * l 
            elif abs ( b.x - guess ) < 0.01 * l : guess = b.x - 0.10 * l
            
            fg  = self.__fun ( guess , *self.__args )
            g   = Point ( guess , fg  ) 
            if 0 == fg or iszero ( fg ) : return self.__result ( guess )     ##  RETURN
            elif samesign ( a.fx , fg ) : a , b = g , b 
            elif samesign ( b.fx , fg ) : a , b = a , g
            else                        : return self.__result ( guess )     ##   RETURN
            
        ## initialize Aitken and Inverse Interpolation structures 
        self.__inverse = [                        a   , b   ] 
        self.__aitken  = [ 0.5  * ( a.x + b.x ) , a.x , b.x ]  
            
        for i in range ( self.__maxiter ) :

            x , a , b = self.__make_step ( guess , a , b )

            ## zero is found ?
            if 0 == x.fx or iszero ( x.fx ) : return self.__result ( x.x , i + 1 ) 

            ## no change with respect to the previous step
            dx1 = self.__xtol + self.__rtol * max ( abs ( x.x ) , abs ( guess ) )            
            if abs ( x.x - guess ) * 4 < dx1 : return self.__result ( x.x , i + 1 ) 
            
            ## bracketing interval is already very small 
            dx2 = self.__xtol + self.__rtol * max ( abs ( a.x ) , abs ( b.x ) )
            if abs ( b.x - a.x   ) * 4 < dx2 : return self.__result ( x.x , i + 1 ) 
            
            guess = x.x
            
        if self.__disp :
            raise RuntimeError ( 'RootFinder.solve: no convergency for %d iterations' % i ) 
            
        return self.__result ( x.x , i  , -2 ) 

    # =========================================================================
    def __result    ( self , x , i = 0 , flag = 0 ) :
        if self.__full_output : return self.__make_report ( x , i , flag ) 
        return x
    
    # =========================================================================
    def __make_step ( self , x0 , a , b ) :

        ## do not keep too many approximations 
        if 10 < len ( self.__inverse ) : del self.__inverse[10:]
        
        old_len =  b.x - a.x
        
        x   = None

        fx0 = None

        # =====================================================================
        ## (1,2) try Halley's and/or Newton methods 
        if x is None and self.__deriv1 :            
            self.__try_halley += 1
            if fx0 is None : fx0 = self.__fun ( x0 , *self.__args )
            xx = halley_newton ( self.__fun    , x0 ,
                                 self.__deriv1 , self.__deriv2 , fx = fx0 , args = self.__args )
            if not xx is None and a.x <= xx <= b.x : 
                self.__success_halley += 1 
                x = Point ( xx , self.__fun ( xx , *self.__args ) ) 
                self.__inverse.insert ( 0 , x  )
                self.__aitken .insert ( 0 , xx )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN
                
        # =====================================================================
        ## (3) try  Steffensen's method
        if x is None and self.__try_steffensen - self.__success_steffensen < 2 :            
            self.__try_steffensen +=1
            if fx0 is None : fx0 = self.__fun ( x0 , *self.__args ) 
            xx = steffensen ( self.__fun , x0 , fx = fx0 , args = self.__args )
            if not xx is None and a.x <= xx <= b.x :                
                self.__success_steffensen +=1                
                x = Point ( xx , self.__fun ( xx , *self.__args ) ) 
                self.__inverse.insert ( 0 , x  )
                self.__aitken .insert ( 0 , xx )                
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN
                
        # =====================================================================
        ## (4) secant is "almost never" fails 
        if x is None :            
            self.__try_secant +=1 
            xx = inverse_linear ( a , b )            
            if not xx is None and a.x <= xx <= b.x :
                self.__success_secant +=1 
                x = Point ( xx , self.__fun ( xx , *self.__args ) ) 
                self.__inverse.insert ( 0 , x  )
                self.__aitken .insert ( 0 , xx )                
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN
                
        # =====================================================================
        ## (5) Inverse polynomial interpolation: we always have some points here
        xx = inverse_polynomial ( *self.__inverse )
        self.__try_inverse += 1 
        if not xx is None and a.x <= xx <= b.x : 
            self.__success_inverse += 1 
            x = Point ( xx , self.__fun ( xx , *self.__args ) ) 
            self.__inverse.insert ( 0 , x  )
            self.__aitken .insert ( 0 , xx )            
            if   0 == x.fx or iszero ( x.fx )     : return x , a , b ## RETURN
            elif samesign ( a.fx , x.fx )         : a , b = x , b  
            elif samesign ( b.fx , x.fx )         : a , b = a , x  
            else                                  : return x , a , b ## RETURN

        # =====================================================================
        ## (6) Try Aitken-delta2 scheme if we have enough approximations
        if 3 <= len ( self.__aitken ) : 
            xx = aitken_delta2 ( *self.__aitken )
            self.__try_aitken += 1 
            if not xx is None and a.x <= xx <= b.x :
                self.__success_aitken += 1 
                x = Point ( xx , self.__fun ( xx ) ) 
                self.__inverse.insert ( 0 , x  )
                del self.__aitken[:] ## ATTENTION: clear Aitken container!                 
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x
                else                              : return x , a , b ## RETURN
               
        # =====================================================================
        ## (7) Force bisection if interval is not very small 
        new_len = b.x - a.x
        if 5 * new_len > old_len and new_len > self.__xtol :
            self.__try_bisection +=1 
            a , b = bisection ( self.__fun , a , b , args = self.__args ) 
            self.__success_bisection +=1 
            new_len = b.x - a.x 

        # =====================================================================    
        ## (8) use bisection as the ultima ratio regum 
        if x is None or not a.x <= x.x <= b.x : 
            self.__try_bisection     +=1            
            xx = 0.5  * ( a.x + b.x )            
            self.__success_bisection +=1             
            x  = Point ( xx , self.__fun ( xx , *self.__args ) ) 
            self.__inverse.insert ( 0 , x  )
            self.__aitken .insert ( 0 , xx )            
            if   0 == x.fx or iszero ( x.fx )     : return x , a , b ## RETURN
            elif samesign ( a.fx , x.fx )         : a , b = x , b  
            elif samesign ( b.fx , x.fx )         : a , b = a , x
            else                                  : return x , a , b ## RETURN
   
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
def find_root ( f                   , ## the function 
                a                   , ## low-edge of bracketing 
                b                   ,                
                args        = ()    ,   
                maxiter     = 100   ,
                xtol        = _xtol ,
                rtol        = _xtol ,
                full_output = False ,
                disp        = True  , 
                deriv1      = None  , 
                deriv2      = None  ) :
    """The main function for root-finding

    Parameters
    ----------
    f           : the function
    a           : low edge of bracketing interval
    b           : high edge of bracketing interval
    args        : extra argument for a function call, function is invoked as `fun(x,*args)`
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


    For each step it use sequence of well-defined ``actions''.
    The result of each ``action'' is considered to be ``successfull'',
    if the obtained root approximation falls into the current
    bracketing interval, otherwise it is considered as ``failure''
    
    In the successfull case the obtained root it stored into list, and the bracketing
    interval is correspodingly updated.
    
    We start iterations from adding the current brackets into collection of roots
    
    The sequential list of ``actions'' is:
    - 'Hailey'     :                             if the second derivative is provided, make Halley's step
    - 'Newton'     : if previous step is failed  and the first  derivative is provdied, make Newton step
    - 'Steffensen' : if previous step is failed                                      , make Steffensen's step
    - 'Inversenterpolation' : at this point we alwasy have at least two or more  approximations
    for the root stored into the list (including initial brackets), so an inverse polynomial interpolation is performed: 
    -- if              there are >=4 points, try with cubic
    -- if it fails and there are >=3 points, try with parabolic
    -- if it fails,                          try with secant   (and this one can't fail!)    
    - 'Aitken'     : if there are >=3 approximations, try to use Aitken delta2 acceleratiion scheme  
    - 'Bisection'  : if the updated bracketing interval is not shrinked significantly for the current step, a bisection step is forced 
    
    After the step the convergency is checked and the proicess is converged  if:
    - either the value of the function at the currect root  approximation is small enough
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

    solver = RootFinder ( fun         = f           ,
                          deriv1      = deriv1      ,
                          deriv2      = deriv2      ,
                          maxiter     = maxiter     ,
                          xtol        = xtol        ,
                          rtol        = rtol        ,
                          full_output = full_output , 
                          disp        = True        )
    
    return solver.find ( a , b )

# =============================================================================
try : 
    # =========================================================================
    with warnings.catch_warnings():
        warnings.simplefilter ( "ignore" )
        from scipy.optimize import brentq as scipy_brentq 
        findroot = scipy_brentq
    # =========================================================================
except ImportError :
    # =========================================================================
    ## logger.warning ("scipy.optimize.brentq is not available, use local ``find_root''-replacement")
    findroot = find_root
    # =========================================================================

# =============================================================================
## solve equation \f$ f(x)=C \f$
#  @code
#  fun = ...
#  C   = 
#  x   = solve ( fun ,  0.0 , 1.0 , C ) 
#  @endcode 
def sp_solve ( fun , xmin ,  xmax , C = 0 , args = () ) :
    """ Solve equation fun(x)=C
    >>> fun = ...
    >>> C   = ... 
    >>> x   = solve ( fun , 0 , 1 , C  ) 
    """
    ##
    if iszero ( C ) :
        return findroot ( fun , xmin ,  xmax , args = args )
    ##
    func = lambda x , *a : fun(x,*a)-C
    return findroot ( func , xmin ,  xmax , args = args )

# =============================================================================
##
solve = sp_solve 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if findroot is find_root : 
        logger.warning ("scipy.optimize.brentq is not available, use local ``find_root''-replacement")

# =============================================================================
##                                                                      The END 
# =============================================================================
