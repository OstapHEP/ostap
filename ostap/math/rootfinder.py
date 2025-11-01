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
from   ostap.math.base    import samesign , iszero , isequal  , isfinite   
from   ostap.utils.basic  import counted 
import sys, collections
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
    fx     : (optional) the value of funnction at the guess point

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
def steffensen ( fun  ,  x  , fx = None , args = () ) :
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
    
    fx = float ( fun ( x ) ) if fx is None else fx  
    if not fx : return x                ## RETURN 
 
    gx = ( fun ( x + fx ) - fx ) / fx
    if not gx : return None             ## RETURN 
    
    return x - fx / gx

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
    
    x0 , f0 = a.x , a.fx
    x1 , f1 = b.x , b.fx

    if   not f0 : return x0
    elif not f1 : return x1 
    
    if f0 == f1 or isequal ( f0 , f1 ) : return None

    return ( x0 * f1  - x1 * f0 ) / ( f1 - f0 ) 

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
def inverse_cubic ( a , b , c , d ) :
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
def aitken_delta2 ( xl , xl1 , xl2 , *others ) :
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
    
    a  = x1 - xl1 
    dd = a  - ( xl1 - xl2 )
    if not dd : return None

    return xl - a *  a / dd 
    
# ========================================================================================
## @class RootResults
#  Helper class to keep results of root-finding procedure
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
    -  'Newton'     : if previous step is failed  and the first  derivative is provided, make Newton step
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
     
     
        self.__fun    = counted ( ffun    ) if self.__full_output             else ffun  
        self.__deriv1 = counted ( fderiv1 ) if self.__full_output and fderiv1 else fderiv1 
        self.__deriv2 = counted ( fderiv2 ) if self.__full_output and fderiv2 else fderiv2

        self.__xtol    = max ( xtol , _xtol )   
        self.__rtol    = max ( rtol , _rtol )
        
        self.__maxiter = maxiter if 1 <= maxiter else 500
        
        self.__disp    = True    if disp        else False
        
            
        ## keep the last three root estimates for the Aitken scheme 
        self.__aitken  = []
        
        ## keep the last three/four approximations for inverse cubic/parabolic interpolation
        self.__inverse = [] 

        self.__total   = collections.defaultdict ( int )
        self.__success = collections.defaultdict ( int )
        
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
    def __make_report ( self , root , i = 0 , flag = 0 ) :
        """ Prepare a report
        """

        result = RootResults ( root , i , self.__fun.calls )
        
        if self.__deriv1 and hasattr ( self.__deriv1 , 'calls' ) :
            result.derivative1_calls = self.__deriv1.calls
            
        if self.__deriv2 and hasattr ( self.__deriv2 , 'calls' ) :
            result.derivative2_calls = self.__deriv2.calls
        
        return root , result

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
        if 0 == fa or iszero ( fa ) :
            if self.__full_output : return self.__report ( a ) ## RETURN
            return a                                           ## RETURN
        
        ## check the root at right edge 
        fb = self.__fun ( b )
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

        ## initialize Aitken and Inverse Interpolation structures 
        self.__inverse = [ a   , b   ] 
        self.__aitken  = [ a.x , b.x ]   

        ## if guess is not specified or outside the interval
        if guess is None or not a.x < guess < b.x :
            
            self.__total [ 'secant' ] += 1 
            ## use the secant as the approximation 
            guess = secant ( a , b ) 
            if guess is None  : guess = 0.5 * ( a.x + b.x )
            else :  self.__success [ 'secant' ]
                   
            fg  = self.__fun ( guess )
            g   = Point ( guess , fg  )
            
            self.__inverse.append ( g   ) 
            self.__aitken .append ( g.x )
            
            if 0 == fg or iszero ( fg ) : return self.__result ( guess ) ##  RETURN
            elif samesign ( a.fx , fg ) : a , b = g , b 
            elif samesign ( b.fx , fg ) : a , b = a , g
            else                        : return self.__result ( guess ) ##   RETURN
            
        ## initialize Aitken and Inverse Interpolation structures 
        self.__inverse = [ a   , b           ] 
        self.__aitken  = [ a.x , b.x , guess ]   

        for i in range ( self.__maxiter ) :

            ## get newx and bew bracketing innterval 
            x , a , b = self.__make_step ( guess , a , b )

            ## zero is found ?
            if not x.fx or iszero ( x.fx ) : return self.__result ( x.x , i + 1 ) 

            ## no change with respect to the previous step
            dx1 = self.__xtol + self.__rtol * max ( abs ( x.x ) , abs ( guess ) )            
            if abs ( x.x - guess ) * 4 < dx1 : return self.__result ( x.x , i + 1 ) 
            
            ## bracketing interval is already very small 
            dx2 = self.__xtol + self.__rtol * abs ( a.x - b.x )
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
        self.__inverse = self.__inverse [ -4 : ]
        
        ## 
        old_len =  b.x - a.x
        
        x   = None
        fx0 = None

        # =====================================================================
        ## (1,2) try Halley's and/or Newton methods 
        if x is None and self.__deriv1 :            
        
            if fx0 is None : fx0 = self.__fun ( x0 )
            if not fx0 or iszero ( fx0 ) : return x0 , a , b 
            
            self.__total [ 'newton' ] += 1 
            xx = halley_newton ( self.__fun , x0 , self.__deriv1 , 
                                 deriv2 = self.__deriv2 , fx = fx0 )
            if   xx is None          : pass
            elif not isfinite ( xx ) : pass 
            elif a.x <= xx <= b.x : 
                self.__success [ 'newton' ] += 1 
                x = Point ( xx , self.__fun ( xx ) ) 
                self.__inverse.append ( x  )  
                self.__aitken .append ( xx )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN
                
        # =====================================================================
        ## (3) try  Steffensen's method
        if x is None :            
            
            if fx0 is None : fx0 = self.__fun ( x0 ) 
            if not fx0 or iszero ( fx0 ) : return x0 , a , b 
            
            self.__total [ 'steffensen' ] += 1 
            xx = steffensen ( self.__fun , x0 , fx = fx0 )
            if   xx is None          : pass
            elif not isfinite ( xx ) : pass 
            elif a.x <= xx <= b.x :                
                self.__success [ 'steffenson' ] += 1                
                x = Point ( xx , self.__fun ( xx ) ) 
                self.__inverse.append ( x  )
                self.__aitken .append ( xx )                
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN
                
        # =====================================================================
        ## (4) secant is "almost never" fails 
        if x is None :            
            
            self.__total [ 'secant' ] += 1 
            xx = inverse_linear ( a , b )
            if   xx is None         : pass 
            elif not isfinite( xx ) : pass 
            elif a.x <= xx <= b.x :
                self.__success [ 'secant' ] += 1 
                x = Point ( xx , self.__fun ( xx ) ) 
                self.__inverse.append ( x  )
                self.__aitken .append ( xx )                
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x  
                else                              : return x , a , b ## RETURN
                 
        # =====================================================================
        ## (5) Inverse polynomial interpolation: we always have some points here
        if x is None and 2 <= len ( self.__inverse ) :
                  
            self.__total [ 'inverse' ] += 1 
            xx = inverse_polynomial ( *self.__inverse [-4:]  ) 
            if   xx is None          : pass 
            elif not isfinite ( xx ) : pass
            elif a.x <= xx <= b.x    :
                self.__success [ 'inverse' ]  += 1 
                x = Point ( xx , self.__fun ( xx ) ) 
                self.__inverse.append ( x  )
                self.__aitken .append ( xx )            
                if   0 == x.fx or iszero ( x.fx )     : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )         : a , b = x , b  
                elif samesign ( b.fx , x.fx )         : a , b = a , x  
                else                                  : return x , a , b ## RETURN

        # =====================================================================
        ## (6) Try Aitken-delta2 scheme if we have enough approximations
        if x is None and 3 <= len ( self.__aitken ) : 
            
            self.__total [ 'aitken' ] += 1 
            xx = aitken_delta2 ( *self.__aitken )
            if   xx is None         : pass 
            elif not isfinite( xx ) : pass
            elif a.x <= xx <= b.x :
                self.__success ['aitken'] += 1 
                x = Point ( xx , self.__fun ( xx ) ) 
                self.__inverse.append ( x  )
                if   0 == x.fx or iszero ( x.fx ) : return x , a , b ## RETURN
                elif samesign ( a.fx , x.fx )     : a , b = x , b  
                elif samesign ( b.fx , x.fx )     : a , b = a , x
                else                              : return x , a , b ## RETURN
                
        # =====================================================================
        ## (7) Force bisection as ultima ratio regum or if theinterval is not reduced  
        new_len = b.x - a.x
        if x is None or not a.x <= x.x <= b.x or ( 3 * new_len > old_len ) :
            self.__total   [ 'bisection' ] += 1 
            a , b = bisection ( self.__fun , a , b ) 
            self.__success [ 'bisection' ] +=1 
            if x is None or not a.x <= x.x <= b.x :
                xx = 0.5 * ( a.x + b.x  )  
                x  = Point ( xx , self.__fun ( xx ) ) 
                
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
                args        = ()    ,   
                kwargs      = {}    , 
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
    - 'Hailey'     :                             if the second derivative is provided, make Halley's step
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
                          disp        = True        )
    
    return solver.find ( a , b )

# =========================================================================
# =========================================================================
try : # ===================================================================
    # =====================================================================
    from scipy.optimize import brentq as scipy_brentq  
    findroot = scipy_brentq
    # =====================================================================
except ImportError : # ====================================================
    # =====================================================================
    ## logger.warning ("scipy.optimize.brentq is not available, use local `find_root'-replacement")
    findroot = find_root
    
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
        logger.warning ("scipy.optimize.brentq is not available, use local `find_root'-replacement")

# =============================================================================
##                                                                      The END 
# =============================================================================
