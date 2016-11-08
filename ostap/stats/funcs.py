#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file  LHCbMath/deriv.py 
#  Simple adaptive numerical differentiation (for pyroot/PyRoUts/Ostap)
#  R. De Levie, "An improved numerical approximation for the first derivative"
#  @see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
#
#  .. and also very simple wrapper to numerical integation using scipy
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
#  
#                    $Revision$
#  Last modification $Date$
#  by                $Author$
# =============================================================================
"""Simple adaptive numerical differentiation (for pyroot/PyRoUts/Ostap/...)

>>> func = lambda x : x*x
>>> print derivative ( func , 1 )

... and also very simple wrapper to numerical integation using scipy

>>> func = lambda x : x*x
>>> print integral ( func , 0 , 1 )

there are also object form:

>>> func = ...
>>> deriv = Derivative ( func )
>>> integ = Integral   ( func , 0 )
>>> print deriv ( 0.1 ) , integ ( 0.1 )  
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-06"
__all__     = (
    ##
    "derivative"    , ## numerical differentiation (as function)
    "integral"      , ## numerical integration     (as function,  using scipy)
    "Derivative"    , ## numerical differentiation (as object) 
    "Integral"      , ## numerical integration     (as as object, using scipy)
    ##
    "IntegralCache" , ## numerical integration     (as object, using scipy and cache)
    ##
    'EvalVE'        , ## evaluate the function taking argument's unicertainty 
    #
    ## stat-quantities   (based on generic SciPy-actions)
    "Moment"        , ## calculate N-th moment of functions/distribitions, etc (scipy)
    "CentralMoment" , ## calculate N-th central moment of functions/distribitions
    "Mean"          , ## calculate "mean"     for functions/distribitions, etc (scipy)
    "Variance"      , ## calculate "variance" for functions/distribitions, etc (scipy)
    "RMS"           , ## calculate "RMS"      for functions/distribitions, etc (scipy)
    "Skewness"      , ## calculate "skewness" for functions/distribitions, etc (scipy)
    "Median"        , ## calculate "median"   for functions/distribitions, etc (scipy)
    "Quantile"      , ## calculate "quantile" for functions/distribitions, etc (scipy)
    "Mode"          , ## calculate "mode"     for functions/distribitions, etc (scipy)
    "Width"         , ## calculate "width"    for functions/distribitions, etc (scipy)
    "CL_symm"       , ## calcualte symmetrical confidence intervals            (scipy)
    "CL_asymm"      , ## calcualte asymmetrical confidence intervals           (scipy)
    ##
    ## stat-quantities   (based on generic SciPy-actions)
    "moment"        , ## calculate N-th moment of functions/distribitions, etc (scipy)
    "central_moment", ## calculate N-th moment of functions/distribitions, etc (scipy)
    "mean"          , ## calculate "mean"     for functions/distribitions, etc (scipy)
    "variance"      , ## calculate "variance" for functions/distribitions, etc (scipy)
    "rms"           , ## calculate "RMS"      for functions/distribitions, etc (scipy)
    "skewness"      , ## calculate "skeness"  for functions/distribitions, etc (scipy)
    "kurtosis"      , ## calculate "kurtosis" for functions/distribitions, etc (scipy)
    "median"        , ## calculate "median"   for functions/distribitions, etc (scipy)
    "quantile"      , ## calculate "quantile" for functions/distribitions, etc (scipy)
    "mode"          , ## calculate "mode"     for functions/distribitions, etc (scipy)
    "width"         , ## calculate "width"    for functions/distribitions, etc (scipy)
    "cl_symm"       , ## calculate symmetrical confidence intervals            (scipy)
    "cl_asymm"      , ## calculate asymmetrical confidence intervals           (scipy)
    ##
    "is_equal"      , ## equality of two floating point numbers  
    "is_zero"       , ## is floaitng number close to zero? 
    ) 
# =============================================================================
from sys import float_info
_eps_   = float_info.epsilon
if not 0.75 < _eps_ * 2**52 < 1.25 :
    import warnings
    warnings.warn ('"epsilon" in not in the expected range!Math could be suboptimal')
# =============================================================================
from LHCbMath.Types import cpp 
VE = cpp.Gaudi.Math.ValueWithError
_mULPs_       = 1000 
_next_double_ = cpp.Gaudi.Math.next_double
def _delta_ ( x , ulps = _mULPs_ ) :
    n1 = _next_double_ ( x ,  ulps )
    n2 = _next_double_ ( x , -ulps )
    return max ( abs ( n1 - x ) , abs ( n2 - x ) )
## is floating point bumber close to zero ?
is_zero        = cpp.LHCb.Math.Zero     ('double')()
## are two floating point number close enough?
is_equal       = cpp.LHCb.Math.Equal_To ('double')()


# =============================================================================
## calculate 1st (and optionally 3rd)  derivative with the given step
#  - f'      is calculated as O(h**2)
#  - f^(III) is calculated as O(h**2)
def _h2_ ( func , x , h , der = False ) :
    """Calculate 1st (and optionally 3rd) derivative 
    - f'      is calculated as O(h**2)
    - f^(III) is calculated as O(h**2)
    """
    fm1 = func ( x - 1 * h )
    fp1 = func ( x + 1 * h )

    d1   = 1.0 * ( fp1 - fm1  ) / (2 * h    )
    if not der : return d1 
    
    fm2 = func ( x - 2 * h )    
    fp2 = func ( x + 2 * h )

    d3  = -2.0 * ( fp1 - fm1 ) + ( fm2 - fm2 )   
    d3 /= (2 * h*h*h)

    return d1,d3

# =============================================================================
## calculate 1st (and optionally 5th) derivative with the given step
#  - f'     is calculated as O(h**4)
#  - f^(V)  is calculated as O(h**2)
def _h4_ ( func , x , h , der = False ) :
    """Calculate 1st (and optionally 5th) derivative  with O(h*h) precision
    - f'     is calculated as O(h**4)
    - f^(V)  is calculated as O(h**2)
    """    
    fm2 = func ( x - 2 * h )    
    fm1 = func ( x - 1 * h )
    fp1 = func ( x + 1 * h )
    fp2 = func ( x + 2 * h )
    
    d1  = 8.0 * ( fp1 - fm1 ) - ( fp2 - fm2 )
    d1 /= 12 * h 
    
    if not der : return d1 
    
    fm3 = func ( x - 3 * h )    
    fp3 = func ( x + 3 * h )
    
    d5  =  5.0 * ( fp1 - fm1 ) - 4 * ( fp2 - fm2 ) + ( fp3 - fm3 )
    d5 /= 2 * h**5 
    
    return d1,d5


# =============================================================================
## calculate 1st (and optionally 7th) derivative with the given step
#  - f'      is calculated as O(h**6)
#  - f^(VII) is calculated as O(h**2)
def _h6_ ( func , x , h , der = False ) :
    """Calculate 1st (and optionally 7th) derivative
    - f'      is calculated as O(h**6)
    - f^(VII) is calculated as O(h**2)
    """        
    fm3 = func ( x - 3 * h )    
    fm2 = func ( x - 2 * h )    
    fm1 = func ( x - 1 * h )
    fp1 = func ( x + 1 * h )
    fp2 = func ( x + 2 * h )
    fp3 = func ( x + 3 * h )
    
    d1  = 45.0 * ( fp1 - fm1 ) - 9 * ( fp2 - fm2 ) + 1 * ( fp3 - fm3 )
    d1 /=  60*h

    if not der : return d1
    
    fm4 = func ( x - 4 * h )    
    fp4 = func ( x + 4 * h )
    
    d7  = -14.0 * ( fp1 - fm1 ) + 14 * ( fp2 - fm2 ) - 6 * ( fp3 - fm3 ) + ( fp4 - fm4 )
    d7 /= 2* h**7
    
    return d1,d7


# =============================================================================
## calculate 1st (and optionally 9th) derivative with the given step
#  - f'     is calculated as O(h**8)
#  - f^(IX) is calculated as O(h**2)
def _h8_ ( func , x , h , der = False ) :
    """Calculate 1st (and optionally 9th) derivative
    - f'     is calculated as O(h**8)
    - f^(IX) is calculated as O(h**2)
    """            
    fm4 = func ( x - 4 * h )    
    fm3 = func ( x - 3 * h )    
    fm2 = func ( x - 2 * h )    
    fm1 = func ( x - 1 * h )
    fp1 = func ( x + 1 * h )
    fp2 = func ( x + 2 * h )
    fp3 = func ( x + 3 * h )
    fp4 = func ( x + 4 * h )
    
    d1  = 672.0 * ( fp1 - fm1 ) - 168 * ( fp2 - fm2  ) + 32 * ( fp3 - fm3 ) - 3  * ( fp4 - fm4 )
    d1 /= 840 * h  
    
    if not der : return d1 
    
    fm5 = func ( x - 5 * h )    
    fp5 = func ( x + 5 * h )
    d9  = ( -fm5 +8*fm4 -27*fm3 + 48*fm2 - 42*fm1 + 42*fp1 - 48*fp2 +27*fp3 -8*fp4 + fp5 )/(  2 * h**9)

    d9  = 42.0 * ( fp1  - fm1 ) - 48 * ( fp2  - fm2  ) + 27 * ( fp3 - fm3 ) - 8 * ( fp4  - fm4 ) + ( fp5 - fm5  )  
    d9 /=  2 * h**9
    
    return d1,d9

# =============================================================================
## calculate 1st (and optionally 11th) derivative with the given step
#  - f'     is calculated as O(h**10)
#  - f^(XI) is calculated as O(h**2)
def _h10_ ( func , x , h , der = False ) :
    """Calculate 1st (and optionally 11th) derivative
    - f'     is calculated as O(h**10)
    - f^(XI) is calculated as O(h**2)
    """
    fm5 = func ( x - 5 * h )    
    fm4 = func ( x - 4 * h )    
    fm3 = func ( x - 3 * h )    
    fm2 = func ( x - 2 * h )    
    fm1 = func ( x - 1 * h )
    fp1 = func ( x + 1 * h )
    fp2 = func ( x + 2 * h )
    fp3 = func ( x + 3 * h )
    fp4 = func ( x + 4 * h )
    fp5 = func ( x + 5 * h )
    
    d1  =  2100.0 * ( fp1 - fm1 ) -  600 * ( fp2 - fm2 ) + 150 * ( fp3 - fm3 ) - 25 * ( fp4 - fm4 ) + 2 * ( fp5 - fm5 )
    d1 /=  2520.0 * h 
    

    if not der : return d1
    
    fm6  =  func ( x - 6 * h )    
    fp6  =  func ( x + 6 * h )
    
    d11  = -132.0 * ( fp1 - fm1 ) + 165 * ( fp2 - fm2 ) - 110 * ( fp3 - fm3 ) + 44 * ( fp4 - fm4 ) - 10 * ( fp5 - fm5 ) + ( fp6 - fm6 ) 
    d11 /= 2*h**11
    
    return d1,d11

# =============================================================================
## calculate 1st (and optionally 13th) derivative with the given step
#  - f'     is calculated as O(h**12)
#  - f^(XIII) is calculated as O(h**2)
def _h12_ ( func , x , h , der = False ) :
    """Calculate 1st (and optionally 13th) derivative
    - f'     is calculated as O(h**12)
    - f^(XIII) is calculated as O(h**2)
    """
    
    fm6 = func ( x - 6 * h )    
    fm5 = func ( x - 5 * h )    
    fm4 = func ( x - 4 * h )    
    fm3 = func ( x - 3 * h )    
    fm2 = func ( x - 2 * h )    
    fm1 = func ( x - 1 * h )
    fp1 = func ( x + 1 * h )
    fp2 = func ( x + 2 * h )
    fp3 = func ( x + 3 * h )
    fp4 = func ( x + 4 * h )
    fp5 = func ( x + 5 * h )
    fp6 = func ( x + 6 * h )
    
    d1  =  23760.0 * ( fp1 - fm1 ) - 7425 * ( fp2 - fm2 ) + 2200 * ( fp3 - fm3 ) - 495 * ( fp4 - fm4 ) + 72 * ( fp5 - fm5 ) - 5 * ( fp6 - fm6 )
    d1 /=  27720.0 * h 
    
    if not der : return d1
    
    fm7  =  func ( x - 7 * h )    
    fp7  =  func ( x + 7 * h )

    d13  = 429.0 * ( fp1 - fm1 ) - 572 * ( fp2 - fm2 ) + 429 * ( fp3 - fm3 ) - 208 * ( fp4 - fm4 ) + 65 * ( fp5 - fm5 ) - 12 * ( fp6 - fm6 ) + ( fp7 - fm7 ) 
    d13 /= 2*h**13
    
    return d1,d13

# =============================================================================
## calculate 1st (and optionally 15th) derivative with the given step
#  - f'     is calculated as O(h**14)
#  - f^(XV) is calculated as O(h**2)
def _h14_ ( func , x , h , der = False ) :
    """Calculate 1st (and optionally 15th) derivative
    - f'     is calculated as O(h**12)
    - f^(XV) is calculated as O(h**2)
    """
    
    fm7 = func ( x - 7 * h )    
    fm6 = func ( x - 6 * h )    
    fm5 = func ( x - 5 * h )    
    fm4 = func ( x - 4 * h )    
    fm3 = func ( x - 3 * h )    
    fm2 = func ( x - 2 * h )    
    fm1 = func ( x - 1 * h )
    fp1 = func ( x + 1 * h )
    fp2 = func ( x + 2 * h )
    fp3 = func ( x + 3 * h )
    fp4 = func ( x + 4 * h )
    fp5 = func ( x + 5 * h )
    fp6 = func ( x + 6 * h )
    fp7 = func ( x + 7 * h )
    
    d1  = 315315.0 * ( fp1 - fm1 ) - 105105 * ( fp2 - fm2 ) + 35035 * ( fp3 - fm3 ) - 9555 * ( fp4 - fm4 ) + 1911 * ( fp5 - fm5 ) - 245 * ( fp6 - fm6 ) + 15 * ( fp7 - fm7 ) 
    d1 /= 360360.0 * h 
    
    if not der : return d1
    
    fm8  =  func ( x - 8 * h )    
    fp8  =  func ( x + 8 * h )

    d15  = -1430.0 * ( fp1 - fm1 ) + 2002 * ( fp2 - fm2 ) - 1638 * ( fp3 - fm3 ) + 910 * ( fp4 - fm4 ) - 350 * ( fp5 - fm5 ) + 90 * ( fp6 - fm6 ) - 14 * ( fp7 - fm7 ) + ( fp8 - fm8 ) 
    d15 /= 2*h**15
    
    return d1,d15

# =============================================================================
## calculate 1st (and optionally 17th) derivative with the given step
#  - f'     is calculated as O(h**16)
#  - f^(XVII) is calculated as O(h**2)
def _h16_ ( func , x , h , der = False ) :
    """Calculate 1st (and optionally 17th) derivative
    - f'       is calculated as O(h**16)
    - f^(XVII) is calculated as O(h**2)
    """
    
    fm8 = func ( x - 8 * h )    
    fm7 = func ( x - 7 * h )    
    fm6 = func ( x - 6 * h )    
    fm5 = func ( x - 5 * h )    
    fm4 = func ( x - 4 * h )    
    fm3 = func ( x - 3 * h )    
    fm2 = func ( x - 2 * h )    
    fm1 = func ( x - 1 * h )
    fp1 = func ( x + 1 * h )
    fp2 = func ( x + 2 * h )
    fp3 = func ( x + 3 * h )
    fp4 = func ( x + 4 * h )
    fp5 = func ( x + 5 * h )
    fp6 = func ( x + 6 * h )
    fp7 = func ( x + 7 * h )
    fp8 = func ( x + 8 * h )
    
    d1  = 640640.0 * ( fp1 - fm1 ) - 224224 * ( fp2 - fm2 ) + 81536 * ( fp3 - fm3 ) - 25480 * ( fp4 - fm4 ) + 6272 * ( fp5 - fm5 ) - 1120 * ( fp6 - fm6 ) + 128 * ( fp7 - fm7 ) - 7 * ( fp8 - fm8 ) 
    d1 /= 720720.0 * h 
    
    if not der : return d1
    
    fm9  =  func ( x - 9 * h )    
    fp9  =  func ( x + 9 * h )

    d17  = 4862.0 * ( fp1 - fm1 ) - 7072 * ( fp2 - fm2 ) + 6188 * ( fp3 - fm3 ) - 3808 * ( fp4 - fm4 ) + 1700 * ( fp5 - fm5 ) - 544 * ( fp6 - fm6 ) + 119 * ( fp7 - fm7 ) - 16 * ( fp8 - fm8 ) + ( fp9 - fm9 )  
    d17 /= 2*h**17
    
    return d1,d17

# =============================================================================
# The actual setup for
# =============================================================================
_funcs_ = (
    _h2_   ,  ## 0 == 1 
    _h2_   ,  ## 1 
    _h4_   ,  ## 2 
    _h6_   ,  ## 3
    _h8_   ,  ## 4 
    _h10_  ,  ## 5 
    _h12_  ,  ## 6
    _h14_  ,  ## 7
    _h16_     ## 8
    )
_numbers_ = (
    ( 0.5*10**(-10./ 3) ,
      0.5*10**(-10./ 3) ,
      0.5*10**(-10./ 5) ,
      0.5*10**(-10./ 7) ,
      0.5*10**(-10./ 9) ,
      0.5*10**(-10./11) , 
      0.5*10**(-10./13) , 
      0.5*10**(-10./15) , 
      0.5*10**(-10./17) ) , 
    (      6 , 4.5324e-17 , 5.1422e-6 , 6.0554e-6 ) , ## I=1, J= 3  3-point rule 
    (     30 , 6.0903e-17 , 8.5495e-4 , 7.4009e-4 ) , ## I=2, J= 5  5-point rule 
    (    140 , 6.9349e-17 , 7.7091e-3 , 5.8046e-4 ) , ## I=3, J= 7  7-point rule 
    (    630 , 7.4832e-17 , 2.6237e-2 , 1.8227e-2 ) , ## I=4, J= 9  9-point rule 
    (   2772 , 7.8754e-17 , 5.7292e-2 , 3.7753e-2 ) , ## I=5, J=11 11-point rule 
    (  12012 , 8.1738e-17 , 9.8468e-2 , 6.2500e-2 ) , ## I=6, J=13 13-point rule 
    (  51480 , 8.4108e-17 , 1.4656e-1 , 9.0454e-1 ) , ## I=7, J=15 15-point rule 
    ( 218790 , 8.6047e-17 , 1.9873e-1 , 1.2000e-1 ) , ## I=8, J=17 17-point rule 
    )
# =============================================================================
## Calculate the first derivative for the function
#  R. De Levie, "An improved numerical approximation for the first derivative"
#  @see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
#  @code
#  >>> fun =  lambda x : x*x
#  >>> print derivative ( fun , x = 1 ) 
#  @endcode
#  @param fun  (INPUT) the function itself
#  @param x    (INPUT) the argument
#  @param h    (INPUT) the guess for the step used in numeric differentiation
#  @param I    (INPUT) the rule to be used ("N-point rule" = 2*I+1)
#  @param err  (INPUT) calcualte the uncertainty?
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
def derivative ( fun , x , h = 0  , I = 2 , err = False ) : 
    """Calculate the first derivative for the function
    #  @code
    #  >>> fun =  lambda x : x*x
    #  >>> print derivative ( fun , x = 1 ) 
    """

    func = lambda x : float ( fun ( x ) )
    
    ## get the function value at the given point 
    f0 = func(x)

    ## adjust the rule 
    I  = min ( max ( I , 1 ) , 8 )
    J  = 2 * I + 1
    
    _dfun_ = _funcs_[I]
    delta  = _delta_ ( x )
    
    ## if the intial step is too small, choose another one 
    if abs ( h ) <  _numbers_[I][3] or abs ( h ) < delta :  
        if is_zero( x ) : h    =             _numbers_[0][I]
        else            : h    = abs ( x ) * _numbers_[I][3] 
        
    ## 1) find the estimate for first and "J"th the derivative with the given step 
    d1,dJ = _dfun_( func , x , h , True )
        
    ## find the optimal step 
    if is_zero  ( dJ ) or  ( is_zero ( f0 ) and is_zero ( x * d1 ) ) :
        if  is_zero ( x )   : hopt =             _numbers_[0][I] 
        else                : hopt = abs ( x ) * _numbers_[I][3]
    else : 
        hopt = _numbers_[I][2]*(  ( abs ( f0 ) + abs ( x * d1 ) ) / abs ( dJ ) )**( 1.0 / J )

    ## finally get the derivative 
    if not err  :  return _dfun_ ( func , x , hopt , False )

    ## estimate the uncrtainty, if needed  
    d1,dJ =  _dfun_ ( func , x , hopt , True )
    e     = _numbers_[I][1]/(J-1)/_numbers_[I][2]
    e2    = e * e * ( J * _eps_ + abs ( f0 ) + abs( x * d1 ) )**(2-2./J) * abs( dJ )**(2./J)
    return VE ( d1 , 4 * e2 )

# =============================================================================
## @class Derivative
#  Calculate the first derivative for the function
#  R. De Levie, "An improved numerical approximation for the first derivative"
#  @see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class Derivative(object) :
    """Calculate the first derivative for the function
    R. De Levie, ``An improved numerical approximation for the first derivative''
    see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
    """
    def __init__ ( self , func , h = 0 , I = 2 , err = False ) :
        """
        Calculate the first derivative for the function
        R. De Levie, ``An improved numerical approximation for the first derivative''
        see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
        
        >>> func = math.sin
        >>> deri = Derivative ( func )
        
        """
        self._func = func
        self._h    = h    
        self._I    = I
        self._err  = err

    ## evaluate the derivative 
    def __call__ ( self , x ) :
        """
        Calculate the first derivative for the function
        R. De Levie, ``An improved numerical approximation for the first derivative''
        see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
        
        >>> func  = math.sin
        >>> deriv = Derivative ( func )
        
        >>> print deriv(0.1) 
        """
        return derivative ( self._func ,
                            x          ,
                            self._h    ,
                            self._I    ,
                            self._err  )

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
## @class EvalVE
#  Evaluate the function taking into account the uncertainty in the argument
#  @code
#  import math 
#  x = VE(1,0.1**2)
#  sin = EvalVE( math.sin , lambda s : math.cos(s) )
#  print 'sin(x) = %s ' % sin(x) 
#  @endcode
#  @see LHCbMath.math_ve 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-02-23
class EvalVE(object) :
    """Evaluate the function taking into account uncertainty in the argument
    
    >>> import math 
    >>> x    = VE(1,0.1**2)
    >>> sin1 = EvalVE( math.sin , lambda s : math.cos(s) )
    >>> print 'sin1(x) = %s ' % sin1(x) 
    >>> sin2 = EvalVE( math.sin )
    >>> print 'sin2(x) = %s ' % sin2(x) 
    see also LHCbMath.math_ve
    """
    ## constructor
    def __init__ ( self , func , deriv = None  , name = '' ) :
        """ Constructor from the function, derivative and name
        >>> sin1 = EvalVE( math.sin , lambda s : math.cos(s) )
        >>> sin2 = EvalVE( math.sin )  ## numerical derivative will be used 
        """
        self._func  = func
        if    deriv : self._deriv  = deriv
        elif  hasattr ( func , 'derivative' ) :
            try :
                self._deriv = func.derivative()  ## derivative as object?
            except:
                self._deriv = func.derivative    ## derivative as function 
        elif  hasattr ( func , 'Derivative' ) :
            try :
                self._deriv = func.Derivative()  ## derivative as object?
            except:
                self._deriv = func.Derivative    ## derivative as function 
        else :
            ## use numerical differentiation
            self._deriv = Derivative(func)
            
        if   name                          : self.__name__ =  name 
        elif hasattr ( func , '__name__' ) : self.__name__ = func.__name__
        else                               : self.__name__ = 'EvalVE'
            
    ## get a value 
    def _value_ ( self , x , *args ) :
        """Evaluate a function"""
        return self._func( float( x ) , *args )
    
    ## evaluate the function 
    def __call__ ( self , x , *args ) :
        #
        ## evaluate the function 
        val  = self._value_ ( x , *args )
        #
        ## no uncertainties? 
        if   isinstance ( x , ( float , int , long ) )   : return VE ( val , 0 )
        # ignore small or invalid uncertanties 
        elif 0 >= x.cov2() or is_zero ( x.cov2() )       : return VE ( val , 0 )
        # evaluate the derivative  
        d    = self._deriv ( float ( x ) , *args )
        ## calculate the variance 
        cov2 = d * d * x.cov2()
        ## get a final result 
        return VE ( val , cov2 )
    
    def __str__ ( self ) : return self.__name__
    
# =============================================================================
## @class Moment
#  Calculate the N-th moment for the distribution 
#  @code
#   xmin,xmax = 0,math.pi 
#   mean  = Moment(1,xmin,xmax)  ## specify min/max
#   value = mean  ( math.sin )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class Moment(object) :
    """
    Calculate the N-th moment for the distribution
    
    >>> xmin,xmax = 0,math.pi 
    >>> mean  = Moment(1,xmin,xmax)  ## specify min/max
    >>> value = mean  ( math.sin )
    
    """
    ## constructor
    def __init__ ( self , N , xmin , xmax , err = False , x0 = 0 , *args ) :
        """
        Contructor 
        """
        if not isinstance ( N , ( int , long ) ) or 0 > N  :
            raise TypeError('Moment: illegal order')
        
        self._N    = N 
        self._xmin = float ( xmin ) if isinstance ( xmin , ( int , long ) ) else xmin 
        self._xmax = float ( xmax ) if isinstance ( xmax , ( int , long ) ) else xmax 
        self._x0   = x0  
        self._err  = err
        self._args = args
        self._moms = {} 

    ## make an integral 
    def _integral_ ( self , func , xmn , xmx , *args ) :
        integrator = Integral ( func , xmn , err = self._err )
        return integrator ( xmx , *args )
    
    ## calculate un-normalized 0-moment  
    def _moment0_ ( self , func , *args ) :
        """
        Calculate un-normalized 0-moment
        """
        return self._integral_ ( func , self._xmin , self._xmax , *args )
    
    ## calculate un-normalized k-moment  
    def _momentK_ ( self , k , func , mu = None , *args ) :
        """
        Calculate unnormalized k-moment
        """
        x0     = self._x0 if mu is None else mu 
        func_N = lambda x,*a : func( x , *a ) * ( ( x - x0 ) ** k  )
        return self._integral_ ( func_N , self._xmin , self._xmax , *args )
    
    ## calculate the moment 
    def __call__ ( self , func , *args ) :
        ## 
        args   = args if args else self._args
        ##
        n0  = self._moment0_ (            func ,            *args ) 
        nN  = self._momentK_ ( self._N  , func , self._x0 , *args ) 
        ##
        return nN/n0

    ## print it!
    def __str__ ( self ) :
        return "Moment(%d,%s,%s,%s,%s)" % ( self._N    ,
                                            self._xmin ,
                                            self._xmax ,
                                            self._err  ,
                                            self._x0   )                                            
                                            

# =============================================================================
## @class CentralMoment
#  Calculate the N-th central moment for the distribution 
#  @code
#   xmin,xmax = 0,math.pi 
#   mc        = CentralMoment(1,xmin,xmax)  ## specify min/max
#   value     = mome  ( math.sin )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class CentralMoment(Moment) :
    """
    Calculate the N-th central moment for the distribution
    
    >>> xmin,xmax = 0,math.pi 
    >>> mc        = CentralMoment(1,xmin,xmax)  ## specify min/max
    >>> value     = mc  ( math.sin )
    
    """
    ## constructor
    def __init__ ( self , N , xmin , xmax , err = False, *args ) :
        """
        Contructor 
        """
        Moment.__init__ ( self , N , xmin , xmax , err , 0.0 , *args ) 

    ## calculate the central moment
    def __call__ ( self , func , *args ) :
        ## 
        args   = args if args else self._args
        ##
        n0  = self._moment0_ (     func ,             *args ) 
        n1  = self._momentK_ ( 1 , func , mu = 0.0  , *args )
        ## get mean
        mu  = float(n1/n0)
        ## use it 
        nN  = self._momentK_ ( self._N  , func , mu , *args ) 
        ##
        return nN/n0

    ## print it!
    def __str__ ( self ) :
        return "CentralMoment(%d,%s,%s,%s)" % ( self._N    ,
                                                self._xmin ,
                                                self._xmax ,
                                                self._err  )
                                                
# =============================================================================
## @class Mean
#  Calculate the mean-value for the distribution 
#  @code
#   xmin,xmax = 0,math.pi 
#   mean  = Mean(xmin,xmax)  ## specify min/max
#   value = mean  ( math.sin )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class Mean(Moment) :
    """
    Calculate the N-th moment for the distribution or function 
    
    >>> xmin,xmax = 0,math.pi 
    >>> mean  = Mean ( xmin , xmax )  ## specify min/max
    >>> value = mean ( math.sin    )
    
    """
    def __init__ ( self , xmin , xmax , err = False ) :
        Moment.__init__ ( self , 1 , xmin , xmax , err )

    def __str__ ( self ) :
        return "Mean(%s,%s,%s)" % ( self._xmin ,
                                    self._xmax ,
                                    self._err  )                                            

# =============================================================================
## @class Variance
#  Calculate the variance for the distribution or function  
#  @code
#   xmin,xmax = 0,math.pi 
#   variance  = Variance ( xmin,xmax )  ## specify min/max
#   value     = variance ( math.sin  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class Variance(Mean) :
    """
    Calculate the variance for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> variance  = Variance ( xmin,xmax )  ## specify min/max
    >>> value     = variance ( math.sin  )
    """
    def __init__ ( self , xmin , xmax , err = False ) :
        Mean.__init__ ( self , xmin , xmax , err )
    ## calculate the variance 
    def __call__ ( self , func , *args ) :
        ## 
        args   = args if args else self._args
        ##
        n0 = self._moment0_ (     func ,            *args ) ## moment-0
        n1 = self._momentK_ ( 1 , func , mu = 0.0 , *args ) ## moment-1 
        ##
        mu = float(n1/n0)                        ## mean-value 
        ## central moment 
        m2 = self._momentK_ ( 2 , func , mu , *args ) 
        ##
        return m2/n0
    
    def __str__ ( self ) :
        return "Variance(%s,%s,%s)" % ( self._xmin ,
                                        self._xmax ,
                                        self._err  )                                            

# =============================================================================
## @class RMS
#  Calculate the variance for the distribution or function  
#  @code
#   xmin,xmax = 0,math.pi 
#   rms       = RMS ( xmin,xmax )  ## specify min/max
#   value     = rms ( math.sin  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class RMS(Variance) :
    """
    Calculate the RMS for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> rms       = RMS ( xmin,xmax )  ## specify min/max
    >>> value     = rms ( math.sin  )
    """
    def __init__ ( self , xmin , xmax , err = False ) :
        Variance.__init__ ( self , xmin , xmax , err )
        
    ## calculate the variance 
    def __call__ ( self , func , *args ) :
        """
        Calculate the RMS for the distribution or function          
        """
        ##
        args   = args if args else self._args
        ##
        var2 = Variance.__call__ ( self , func , *args )
        import LHCbMath.math_ve as ME 
        return ME.sqrt ( var2 ) 

    def __str__ ( self ) :
        return "RMS(%s,%s,%s)" % ( self._xmin ,
                                   self._xmax ,
                                   self._err  )                                            

# =============================================================================
## @class Skewness
#  Calculate the skewness for the distribution or function  
# \f$ k = \frac{\mu_3}{\sigma^3} \f$
#  @code
#   xmin,xmax = 0,math.pi 
#   skew     = Skewness ( xmin,xmax )  ## specify min/max
#   value     = skew ( math.sin  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-08-09
class Skewness(Variance) :
    """
    Calculate the variance for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> skew      = Skewness ( xmin,xmax )  ## specify min/max
    >>> value     = skew     ( math.sin  )
    """
    def __init__ ( self , xmin , xmax , err = False ) :
        Variance.__init__ ( self , xmin , xmax , err )
    ## calculate the variance 
    def __call__ ( self , func , *args ) :
        ## 
        args   = args if args else self._args
        ##
        n0 = self._moment0_ (     func ,          *args ) ## norm
        n1 = self._momentK_ ( 1 , func , mu = 0 , *args ) ## m1
        ## get mean-value 
        mu = float(n1/n0) ## mean-value
        ## 
        m2 = self._momentK_ ( 2 , func , mu     , *args ) ## mu2 
        m3 = self._momentK_ ( 3 , func , mu     , *args ) ## mu3 
        ##
        m2 /= n0 ## normalize 
        m3 /= n0 ## normalize
        ## 
        return m3/(m2**(3.0/2))
    
    def __str__ ( self ) :
        return "Skewness(%s,%s,%s)" % ( self._xmin ,
                                        self._xmax ,
                                        self._err  )                                            

# =============================================================================
## @class  Kurtosis
#  Calculate the (excess) kurtosis for the distribution or function
#  \f$ k = \frac{\mu_4}{\sigma_2}-3\f$
#  @code
#   xmin,xmax = 0,math.pi 
#   kurt      = Kurtosis ( xmin,xmax )  ## specify min/max
#   value     = kurt ( math.sin  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-08-09
class Kurtosis(Skewness) :
    """
    Calculate the variance for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> kurt      = Kurtosis ( xmin,xmax )  ## specify min/max
    >>> value     = kurt     ( math.sin  )
    """
    def __init__ ( self , xmin , xmax , err = False ) :
        Skewness.__init__ ( self , xmin , xmax , err )
    ## calculate the variance 
    def __call__ ( self , func , *args ) :
        ## 
        args   = args if args else self._args
        ##
        n0 = self._moment0_ (     func ,          *args ) ## norm
        n1 = self._momentK_ ( 1 , func , mu = 0 , *args ) ## m1
        ## get mean-value 
        mu = float(n1/n0) ## mean-value
        ## 
        m2 = self._momentK_ ( 2 , func , mu     , *args ) ## mu2 
        m4 = self._momentK_ ( 4 , func , mu     , *args ) ## mu3 
        ##
        m2 /= n0 ## normalize 
        m4 /= n0 ## normalize
        ## 
        return m4/(m2*m2)-3.0 
    
    def __str__ ( self ) :
        return "Kurtosis(%s,%s,%s)" % ( self._xmin ,
                                        self._xmax ,
                                        self._err  )                                            


# =============================================================================
## @class Median 
#  Calculate the median for the distribution or function  
#  @code
#  xmin,xmax = 0,math.pi 
#  median    = Median ( xmin,xmax )  ## specify min/max
#  value     = median ( math.sin  )
#  @endcode 
#  @see https://en.wikipedia.org/wiki/Median#Inequality_relating_means_and_medians
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-12
class Median(RMS) :
    """
    Calculate median for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> median    = Median ( xmin,xmax )  ## specify min/max
    >>> value     = median ( math.sin  )
    """
    def __init__ ( self , xmin , xmax ) :
        RMS.__init__ ( self , xmin , xmax , err = False )

    ## calculate he median
    def _median_ ( self , func , xmin , xmax , *args ) :
        ## need to know the integral
        iint   = IntegralCache ( func ,  xmin , False , *args )
        half   = 2.0 / iint    ( xmax ) 
        
        from scipy import optimize
        ifun   = lambda x : iint( x ) * half - 1.0

        ## @see https://en.wikipedia.org/wiki/Median#Inequality_relating_means_and_medians
        try: 
            meanv = Mean . __call__ ( self , func , *args )
            sigma = RMS  . __call__ ( self , func , *args )
            import math
            xmn   = meanv - 2 * sigma ## use 2 instead of 1 
            xmx   = meanv + 2 * sigma ## use 2 instead of 1
            #
            if isinstance ( xmin , float ) : xmn = max ( xmn , xmin ) 
            if isinstance ( xmax , float ) : xmx = min ( xmx , xmax )
            #
            result = optimize.brentq ( ifun , xmn , xmx )
        except :
            result = optimize.brentq ( ifun , xmin , xmax )
            
        return result

        
    ## calculate the median 
    def __call__ ( self , func , *args ) :
        return self._median_ ( func , self._xmin , self._xmax )

    def __str__ ( self ) :
        return "Median(%s,%s)" % ( self._xmin , self._xmax )
    
# =============================================================================
## get the quantile
#  Calculate the quantile for the distribution or function  
#  @code
#  xmin,xmax = 0,math.pi 
#  quantile  = Quantile ( 0.1 , xmin,xmax )  ## specify min/max
#  value     = quantile ( math.sin  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-12
class Quantile(Median) :
    """
    Calculate quantiles for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> quantile  = Quantile ( 0.1 , xmin,xmax )  ## specify min/max
    >>> value     = quantile ( math.sin  )
    """
    def __init__ ( self , Q , xmin , xmax ) :
        Median.__init__ ( self , xmin , xmax )
        #
        if Q < 0 : raise ArrtibuteError ( 'Quantile is invalid %s' % Q )
        if Q > 1 : raise ArrtibuteError ( 'Quantile is invalid %s' % Q )
        self._Q = float( Q ) 
        
    def __str__ ( self ) :
        return "Quantile(%s,%s,%s)" % ( self._Q , self._xmin , self._xmax )

    ## calculate the median 
    def __call__ ( self , func , *args ) :
        ##

        if    0.5 == self._Q : return Median.__call__ ( self , func , *args ) 
        elif  0.0 == self._Q : return self._xmin
        elif  1.0 == self._Q : return self._xmax

        ## need to know the integral
        iint = IntegralCache ( func, self._xmin, False , *args )
        quan = 1.0 / iint    (  self._xmax ) / self._Q 
        
        from scipy import optimize
        ifun   = lambda x : iint( x ) * quan - 1.0

        xmn = self._xmin
        xmx = self._xmax

        p   = 0.5
        l   = 0.5

        ## make some bracketing before next step 
        while ( not isinstance ( xmn , float ) ) or ( not isinstance ( xmx , float ) ) or l>0.1 :   
        
            l /= 2            
            m = self._median_ ( func , xmn , xmx , *args ) 
            
            if   self._Q < p :
                xmn   = xmn 
                xmx   = float( m ) 
                p    -= l 
            elif self._Q > p :
                xmn   = float ( m ) 
                xmx   = xmx
                p    +=  l  
            else : return m               ## RETURN 

        ## finally, calculate quantile 
        result = optimize.brentq ( ifun , xmn , xmx )
            
        return result


# =============================================================================
## @class Mode
#  Calculate the mode for the distribution or function  
#  @code
#  xmin,xmax = 0,math.pi 
#  mode      = Mode ( xmin,xmax )  ## specify min/max
#  value     = mode ( math.sin  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-12
class Mode(Median) :
    """
    Calculate the mode for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> mode      = Mode ( xmin,xmax )  ## specify min/max
    >>> value     = mode ( math.sin  )
    """
    def __init__ ( self , xmin , xmax ) :
        Median.__init__ ( self , xmin , xmax )
        
    ## calculate the mode 
    def __call__ ( self , func , *args ) :
        ##
        
        ## use mean    as intial approximation for mode 
        m1     = Mean   .__call__ ( self , func , *args )
        
        ## use median as intial approximation for mode 
        ## m2     = Median.__call__ ( self , func , *args )
        
        ## use the point intermediate between mean and median as approximation 
        ## m0     = 0.5 * ( m1 + m2 )

        m0 = m1 
        ifun = lambda x,*a : -1.0 * float( func ( x , *a ) )
        
        from scipy import optimize
        result = optimize.minimize (
            ifun                                   , 
            x0     = float ( m0 )                  ,
            bounds = [ (self._xmin , self._xmax) ] ,
            args   = args )
        
        return result.x[0]
    
    def __str__ ( self ) :
        return "Mode(%s,%s)" % ( self._xmin , self._xmax )

# =============================================================================
## @class Width
#  Calculate the full width at half heigh for the distribution or function  
#  @code
#  xmin,xmax = 0,math.pi 
#  width     = Width ( xmin,xmax )  ## specify min/max
#  x1,x2     = width ( math.sin  )
#  fwhm      = x2-x1
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-12
class Width(Mode) :
    """
    Calculate the mode for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> width     = Width ( xmin,xmax )  ## specify min/max
    >>> x1,x2     = width ( math.sin )
    >>> fwhm      = x2-x1
    """
    def __init__ ( self , xmin , xmax , height_factor = 0.5 ) :
        Mode.__init__ ( self , xmin , xmax )
        self._hfactor = height_factor
        
    ## calculate the width
    def __call__ ( self , func , *args ) :
        ##

        ## get the position of the mode
        m0  = Mode.__call__ ( self , func , *args )

        ## function  value at the maximum
        v0      = func ( m0 , *args )

        ## half height 
        vheight = 1.0 * v0 * self._hfactor

        
        ## use scipy to find solution 
        from scipy import optimize        
        ifun = lambda x,*a : float(func (x,*a))-vheight
        x1 = optimize.brentq ( ifun , self._xmin , m0         , args = args )
        x2 = optimize.brentq ( ifun , m0         , self._xmax , args = args ) 
        
        return x1,x2

    def __str__ ( self ) :
        return "Width(%s,%s,%s)" % ( self._xmin , self._xmax , self._hfactor)
    

# =============================================================================
## @class CL_symm
#  Calcualate symmetic confidence interval around x0
#  for the given function on (xmin,xmax) interval
#  function is assumed to be zero outside the interval
#  @code
#  fun = lambda x : exp( -0.5*x*x)
#  reg = CL_symm ( 0.68 , -10 , 10 )
#  print reg ( fun )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-08-03
class CL_symm(object) :
    """
    Calculate symmetic confidence interval around x0
    for the given function on (xmin,xmax) interval
    function is assumed to be zero outside the interval
    >>> fun = lambda x : exp( -0.5*x*x)
    >>> reg = CL_symm ( 0.68 , -10 , 10 )
    >>> print reg ( fun )
    """
    def __init__ ( self , prob ,  xmin , xmax , x0 = None , *args ) :
        
        if not 0.0 < prob < 1.0 :
            raise AttributeError ("Invalid value of prob/CL=%g" % prob)

        self._prob = prob 
        self._xmin = float ( xmin ) if isinstance ( xmin , ( int , long ) ) else xmin 
        self._xmax = float ( xmax ) if isinstance ( xmax , ( int , long ) ) else xmax 
        self._x0   = float ( x0   ) if isinstance ( x0   , ( int , long ) ) else x0  
        self._args = args
        
    def __str__ ( self ) :
        return "CL_sym(%s,%s,%s,%s)" % ( self._prob ,
                                         self._xmin , self._xmax , self._x0   )

    def __call__ ( self , func , *args ) :

        ## additional arguments
        args   = args if args else self._args
        
        #
        ## define integration rules
        #
        if hasattr ( func , 'integral' ) :
            _integral_ = lambda f , low , high : f.integral (      low , high , *args )
        else                             :
            _integral_ = lambda f , low , high :   integral ( f  , low , high , False , *args )

        #
        ## xmin/max
        #
        xmin,xmax = self._xmin, self._xmax
        
        #
        ## calculate x0 as "mean"-value
        #
        x0 = self._x0 
        if x0 is None :
            if hasattr ( func , 'mean' ) : x0 = func.mean()
            else                         :
                m  = Mean ( xmin , xmax , False )
                x0 = m    ( func , *args )  

        #
        ## check validity of x0
        #
        if not xmin <= x0 <= xmax :
            raise AttributeError ("Invalid x0 value %s<=%s<=%s" % ( xmin , x0 , xmax ) )

        #
        ## get the normalization
        #
        norm  = _integral_ ( func , xmin , xmax )
        if 0 >= norm :
            raise AttributeError ("Normalization integral is not positive %s" % norm )

        #
        ## Equation:  ifun(x) \equiv \int_{x0-x}^{x0+x}f(t)dt - N*prob = 0   
        #
        yval  = self._prob * norm
        def ifun ( x ) :
            if 0 >= x : return -yval 
            return _integral_ ( func , max ( xmin , x0 - x ) , min ( xmax , x0 + x )  ) - yval  
        
        ## use scipy to find solution 
        from scipy import optimize        
        s = optimize.brentq (  ifun , 0 , max ( xmax - x0 , x0 - xmin ) )
        
        return VE ( x0 , s * s )


# =============================================================================
## @class CL_asymm
#  Calcualate asymmetic confidence interval (x1,x2) 
#  for the given function on (xmin,xmax) interval,
#  such as 
#  \f$ \begin{array}{l} f(x_1)=f(x_2)            \\
#  \int_{x_1}^{x_2}f(t)dr = p \int_{x_{min}}^{x_{max}}f(t)ft \end{array}\f$
#  Assumtions on function
#  - nonnegative
#  - unimodal
#  - positive integral between (xmin,xmax)
#  - zero outside interval     (xmin,xmax)
#  @code
#  fun   = lambda x : exp( -0.5*x*x)
#  reg   = CL_asymm ( 0.68 , -10 , 10 )
#  print reg ( fun )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-08-03
class CL_asymm(object) :
    """
    Calculate asymmetic confidence interval around x0
    for the given function on (xmin,xmax) interval
    function is assumed to be zero outside the interval
    >>> fun = lambda x : exp( -0.5*x*x)
    >>> reg = CL_asymm ( 0.68 , -10 , 10 )
    >>> print reg ( fun )
    """
    def __init__ ( self , prob ,  xmin , xmax , *args ) :
        
        if not 0.0 < prob < 1.0 :
            raise AttributeError ("Invalid value of prob/CL=%g" % prob)

        self._prob = prob 
        self._xmin = float ( xmin ) if isinstance ( xmin , ( int , long ) ) else xmin 
        self._xmax = float ( xmax ) if isinstance ( xmax , ( int , long ) ) else xmax 
        self._args = args

    def __str__ ( self ) :
        return "CL_asymm(%s,%s,%s)" % ( self._prob , self._xmin , self._xmax )

    ## solve equation f(x)=a 
    def _solve_  ( self , func , fval , xmn , xmx , *args ) :

        ifun = lambda x,*a : func(x,*a) - fval

        ## use scipy to find solution 
        from scipy import optimize        
        return optimize.brentq (  ifun , xmn , xmx , args = args )

                   
    def __call__ ( self , func , *args ) :

        ## additional arguments
        args   = args if args else self._args
        
        #
        ## define integration rules
        #
        if hasattr ( func , 'integral' ) :
            _integral_ = lambda f , low , high : f.integral (      low , high , *args )
        else                             :
            _integral_ = lambda f , low , high :   integral ( f  , low , high , False , *args )

        #
        ## xmin/max
        #
        xmin,xmax = self._xmin, self._xmax
        
        #
        # calculate mode
        if hasattr ( func , 'mode' ) : xmode = func.mode()
        else :
            md    = Mode ( xmin , xmax )
            xmode = md   ( func , *args )

        if not xmin <= xmode <= xmax :
            raise AttributeError ("Invalid mode value %s<=%s<=%s" % ( xmin , xmode , xmax ) )

        #
        ## get the normalization
        #
        norm  = _integral_ ( func , xmin , xmax )
        if 0 >= norm :
            raise AttributeError ("Normalization integral is not positive %s" % norm )

        normL = _integral_ ( func , xmin  , xmode )
        normR = _integral_ ( func , xmode , xmax  )

        ## solve equation f(x)=a 
        def _solve_  ( func , fval , xmn , xmx , *args ) :
            ##
            if is_equal (  xmn , xmx   )  : return xmn
            ## 
            ifun = lambda x,*a : func(x,*a) - fval
            ## 
            fmn = ifun  ( xmn )
            if is_zero  ( ifun ( xmn ) )  : return xmn
            fmx = ifun  ( xmx )
            if is_zero  ( ifun ( xmx ) )  : return xmx 
            ##
            if 0 < fmx * fmn : ## more or less arbitrary choice 
                return xmx if abs ( fmx ) <= abs ( fmn ) else xmn 
            #
            ## use scipy to find solution 
            from scipy import optimize
            return optimize.brentq (  ifun , xmn , xmx , args = args )

        yval = self._prob * norm
        fm   = func ( xmode ) 
        def iifun ( f ) :

            if   is_zero  ( f      ) : x1 , x2 = xmin,xmax
            elif is_equal ( f , fm ) : return -yval 
            else : 
                x1 = _solve_ ( func ,  f , xmin  , xmode )
                x2 = _solve_ ( func ,  f , xmode , xmax  )

            return _integral_( func , x1 , x2 ) - yval 
            
        from scipy import optimize
        l  = optimize.brentq (  iifun , 0 , func ( xmode ) )
        x1 = _solve_ ( func ,  l , xmin  , xmode )
        x2 = _solve_ ( func ,  l , xmode , xmax  )

        return x1 , x2 

 
# =============================================================================
## calculate some statistical quantities of variable,
#  considering function to be PDF 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def sp_action ( func , actor , xmin = None , xmax = None ) :
    """
    Calculate some statistical quantities of variable, considering function to be PDF 
    """
    ##
    import numpy
    ##
    if   isinstance  ( xmin , (float,int,long) ) : xmn =  float ( xmin            ) 
    elif hasattr     ( func ,'GetXmin' )         : xmn =  float ( func.GetXmin () )
    elif hasattr     ( func ,'xmin'    )         : xmn =  float ( func.xmin    () ) 
    else                                         : xmn = -numpy.inf
    ##
    if   isinstance  ( xmax , (float,int,long) ) : xmx =  float ( xmax            )
    elif hasattr     ( func ,'GetXmax' )         : xmx =  float ( func.GetXmax () ) 
    elif hasattr     ( func ,'xmax'    )         : xmx =  float ( func.xmax    () )
    else                                         : xmx = +numpy.inf
    ##
    xmn = float ( xmn ) if isinstance ( xmn , ( int , long ) ) else xmn 
    xmx = float ( xmx ) if isinstance ( xmx , ( int , long ) ) else xmx
    #
    ## instantiate calculator and use it 
    calc = actor ( xmn , xmx )
    ##
    return calc  ( func )

# =============================================================================
## get the N-moment of variable, considering function to be PDF 
# @code 
# >>> fun  = ...
# >>> m5   = moment( fun , 5 , xmin = 10 , xmax = 50 )
# @endcode
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2015-07-11
def moment ( func , N , xmin = None , xmax = None , err = False , x0 = 0 ) :
    """ Get the moment for the distribution using scipy/numpy
    >>> fun  = ...
    >>> mom5 = moment ( fun , 5 , xmin = 10 , xmax = 50 )
    """
    ## get the functions from LHCbMath.deriv 
    actor = lambda x1,x2 : Moment ( N , x1 , x2 , err , x0 ) 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the N-th centralmoment of variable, considering function to be PDF 
# @code 
# >>> fun  = ...
# >>> m5   = central_moment( fun , 5 , xmin = 10 , xmax = 50 )
# @endcode
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2015-07-11
def central_moment ( func , N , xmin = None , xmax = None , err = False ) :
    """Get the central moment for the distribution using scipy/numpy
    >>> fun  = ...
    >>> mom5 = central_moment ( fun , 5 , xmin = 10 , xmax = 50 )
    """
    ## get the functions from LHCbMath.deriv 
    actor = lambda x1,x2 : CentralMoment ( N , x1 , x2 , err ) 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the mean value of variable, considering function to be PDF 
#  @code 
#  >>> fun = ...
#  >>> m   = mean( fun , xmin = 10 , xmax = 50 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def mean ( func , xmin = None , xmax = None , err = False ) :
    """Get the mean-value for the distribution using scipy/numpy
    >>> fun = ...
    >>> m   = mean( fun , xmin = 10 , xmax = 50 )
    """
    ## get the functions from LHCbMath.deriv 
    actor = lambda x1,x2 : Mean ( x1 , x2 , err ) 
    ## use it! 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the variance of the variable, considering function to be PDF 
#  @code 
#  >>> fun = ...
#  >>> v   = variance( fun , xmin = 10 , xmax = 50 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def variance ( func , xmin = None , xmax = None , err = False ) :
    """Get the variance for the distribution using scipy/numpy
    >>> fun = ...
    >>> v   = variance( fun , xmin = 10 , xmax = 50 )
    """
    ##
    ## get the functions from LHCbMath.deriv 
    actor = lambda x1,x2 : Variance ( x1 , x2 , err ) 
    ## use it! 
    return sp_action ( func , actor  , xmin , xmax )

# =============================================================================
## get the rms of the variable, considering function to be PDF 
#  @code 
#  >>> fun = ...
#  >>> v   = rms( fun , xmin = 10 , xmax = 50 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def rms ( func , xmin = None , xmax = None , err = False ) :
    """Get RMS for the distribution using scipy/numpy
    >>> fun = ...
    >>> v   = rms( fun , xmin = 10 , xmax = 50 )
    """
    ##
    ## get the functions from LHCbMath.deriv 
    actor = lambda x1,x2 : RMS ( x1 , x2 , err ) 
    ## use it! 
    return sp_action ( func , actor  , xmin , xmax )

# =============================================================================
## get the skewness of the variable, considering function to be PDF 
#  @code 
#  >>> fun  = ...
#  >>> skew = skewness ( fun , xmin = -10 , xmax = 10 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def skewness ( func , xmin = None , xmax = None , err = False ) :
    """Get the skewness for the distribution using scipy/numpy
    >>> fun = ...
    >>> v   = skewness ( fun , xmin = -10 , xmax = 10 )
    """
    ##
    ## get the functions from LHCbMath.deriv 
    actor = lambda x1,x2 : Skewness ( x1 , x2 , err ) 
    ## use it! 
    return sp_action ( func , actor  , xmin , xmax )


# =============================================================================
## get the excessive kurtosis of the variable, considering function to be PDF
#  \f$ k = \frac{\mu_4}{\sigma^4}-3\f$ 
#  @code 
#  >>> fun  = ...
#  >>> kurt = kurtosis ( fun , xmin = -10 , xmax = 10 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def kurtosis ( func , xmin = None , xmax = None , err = False ) :
    """Get the (exessive) kurtosis for the distribution using scipy/numpy
    >>> fun  = ...
    >>> kurt = kurtosis ( fun , xmin = 10 , xmax = 50 )
    """
    ##
    ## get the functions from LHCbMath.deriv 
    actor = lambda x1,x2 : Kurtosis ( x1 , x2 , err ) 
    ## use it! 
    return sp_action ( func , actor  , xmin , xmax )

# =============================================================================
## get the median the variable, considering function to be PDF 
#  @code 
#  >>> fun = ...
#  >>> med = median ( fun , xmin = 10 , xmax = 50 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def median ( func , xmin = None , xmax = None ) :
    """Get the median for the distribution using scipy/numpy
    >>> fun = ...
    >>> v   = median( fun , xmin = 10 , xmax = 50 )
    """
    ## get the functions from LHCbMath.deriv 
    actor = lambda x1,x2 : Median ( x1 , x2 ) 
    ##
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the quantile of variable, considering function to be PDF 
#  @code 
#  >>> fun  = ...
#  >>> quan = quantile( fun , 0.1 , xmin = 10 , xmax = 50 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def quantile ( func , Q , xmin = None , xmax = None , err = False , x0 = 0 ) :
    """Get quantile for the distribution using scipy/numpy
    >>> fun  = ...
    >>> quan = quantile ( fun , 0.1 , xmin = 10 , xmax = 50 )
    """
    ## get the functions from LHCbMath.deriv 
    actor = lambda x1,x2 : Quantile ( Q , x1 , x2 ) 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the mode the variable, considering function to be PDF 
#  @code 
#  >>> fun = ...
#  >>> m   = mode( fun ,  xmin = 10 , xmax = 50 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def mode ( func , xmin = None , xmax = None ) :
    """Get the mode for the distribution using scipy/numpy
    >>> fun = ...
    >>> v   = mode( fun ,  xmin = 10 , xmax = 50 )
    """
    ## get the functions from LHCbMath.deriv 
    ## use it! 
    ## get the functions from LHCbMath.deriv 
    actor = lambda x1,x2 : Mode  ( x1 , x2 ) 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the width, considering function to be PDF 
#  @code 
#  >>> fun   = ...
#  >>> x1,x2 = width( fun ,  xmin = 10 , xmax = 50 )
#  >>> fwhm  = x2-x1
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def width ( func , xmin = None , xmax = None , height_factor = 0.5 ) :
    """ Get the width for the distribution using scipy/numpy
    >>> fun   = ...
    >>> x1,x2 = width ( fun ,  xmin = 10 , xmax = 50 )
    >>> fwhm  = x2-x1   
    """
    ## get the functions from LHCbMath.deriv 
    ## use it! 
    ## get the functions from LHCbMath.deriv
    actor = lambda x1,x2 : Width  ( x1 , x2 , height_factor ) 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the symmetric confidence interval around x0 for (xmin,xmax) interval 
#  @code 
#  fun  = lambda x : exp( - 0.5 * x * x ) 
#  x_1  = cl_symm ( fun , 0.68 , -10 , 10 )
#  print x_1 
#  x_2  = cl_symm ( fun , 0.68 ,   0 , 10 , x0 = 0 )
#  print x_2 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-08-03
def cl_symm ( func , prob , xmin = None , xmax = None , x0 = None ) :
    """ Symmetric confidence interval around x0    
    >>> fun  = lambda x : exp( - 0.5 * x * x )
    >>> x_1  = cl_symm ( fun , 0.68 , -10 , 10 )
    >>> print x_1 
    >>> x_2  = cl_symm ( fun , 0.68 ,   0 , 10 , x0 = 0 )
    >>> print x_2 
    """
    ## get the functions
    actor = lambda x1,x2 : CL_symm ( prob , x1 , x2 , x0 = x0 ) 
    ## and use it! 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the asymmetric confidence interval around x0 for (xmin,xmax) interval 
#  @code 
#  fun  = lambda x : exp( - 0.5 * x * x ) 
#  x_1,x_2  = cl_asymm ( fun , 0.68 , -10 , 10 )
#  print x_1 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-08-03
def cl_asymm ( func , prob , xmin = None , xmax = None ) :
    """Asymmetric confidence interval around x0    
    >>> fun  = lambda x : exp( - 0.5 * x * x )
    >>> x_1,x_2  = cl_asymm ( fun , 0.68 , -10 , 10 )
    >>> print x_1,x_2 
    """
    ## get the functions
    actor = lambda x1,x2 : CL_asymm ( prob , x1 , x2 ) 
    ## and use it! 
    return sp_action ( func , actor , xmin , xmax )


# =============================================================================
if '__main__' == __name__ :
    
    print 80*'*'
    print __doc__
    print ' Author  : ' , __author__
    print ' Version : ' , __version__
    print ' Date    : ' , __date__    
    print ' Symbols : ' , __all__    
    print 80*'*'
    
    import math

    functions = [
        ( lambda x : math.cos(100.*x)  , lambda x : -100*math.sin(100.*x)                , 0   ) , 
        ( lambda x : x**3              , lambda x : 3.0*x*x                              , 1   ) , 
        ( lambda x : math.exp(x)       , lambda x : math.exp(x)                          , 5   ) ,
        ( lambda x : x**8              , lambda x : 8.0*x**7                             , 2.2 ) ,
        ( lambda x : math.tanh(10.*x)  , lambda x : 10*(1.-math.tanh(10.*x)**2)          , 0.1 ) ,
        ( lambda x : math.tan (10.*x)  , lambda x : 10*(1.+math.tan (10.*x)**2)          , 2.1 ) ,
        ( lambda x : 1./math.sin(2.*x) , lambda x : -2./math.sin(2.*x)**2*math.cos(2.*x) , 0.2 ) , 
        ( lambda x : 1.11*x            , lambda x : 1.11                                 , 0.4 ) , 
        ( lambda x : 1000./x**2        , lambda x : -2000./x**3                          , 0.8 ) , 
        ( lambda x : 1.11111           , lambda x : 0.0                                  , 0.6 ) , 
        ( lambda x : x**50             , lambda x : 50.*x**49                            , 1.3 ) , 
        ]


    for o in functions :

        fun = o[0] ## function 
        der = o[1] ## derivative 
        x   = o[2] ## argument 
        
        print 80*'*'
        
        for i in range(1,6 ) :
            
            res = derivative ( fun , x , 1.e-20 , i  , err = True )
            f1  = res 
            d   = der(x)  ## the exact value for derivative 
            if is_zero ( d ) : 
                print 'Rule=%2d' % ( 2*i+1 ) , f1 , d , (f1.value()-d)    
            else      :
                print 'Rule=%2d' % ( 2*i+1 ) , f1 , d , (f1.value()-d)/d  
                    
    print 80*'*'

    ## the function 
    func    = math.sin
    ## analysitical derivative 
    deriv_a = math.cos
    ## numerical first derivative     
    deriv_1 = Derivative ( func    , I = 5 )  
    ## numerical second derivative     
    deriv_2 = Derivative ( deriv_1 , I = 5 )  

    import random

    for i in range ( 0 , 20 ) : 

        x = random.uniform( 0 , 0.5*math.pi )
        
        fmt = "x=%10.5g f=%10.5g delta(f')= %+10.4g delta(f\")= %+10.4g "
        print fmt %  ( x       ,
                       func(x) ,
                       1-deriv_1(x)/deriv_a(x),
                       1+deriv_2(x)/func(x) )
        

    #
    ## test mean/vars
    #
    import math
    mean_ = Mean     (0, math.pi)
    print 'sin@[0,pi]                mean: %s ' % mean_ (math.sin) 

    var2 = Variance (0, math.pi)
    print 'sin@[0,pi]            variance: %s ' % var2 (math.sin) 

    med  = Median   (0, math.pi)
    print 'sin@[0,pi]              median: %s ' % med  (math.sin) 

    mode_ = Mode     (0, math.pi)
    print 'sin@[0,pi]                mode: %s ' % mode_ (math.sin) 

    def fwhm ( fun ) :
        _w    =   Width   (0, math.pi)
        x1,x2 = _w ( fun )
        return x2-x1
    print 'sin@[0,pi]                fwhm: %s ' % fwhm (math.sin) 

    rms_ = RMS     (0, math.pi)
    print 'sin@[0,pi]                 rms: %s ' % rms_  (math.sin) 

    mom5 = Moment  ( 5, 0, math.pi)
    print 'sin@[0,pi]          5th moment: %s ' % mom5 (math.sin) 

    mom5 = CentralMoment  ( 5, 0, math.pi)
    print 'sin@[0,pi] 5th central moment : %s ' % mom5 (math.sin) 

    mom1 = Moment  ( 1, 0, math.pi)
    print 'sin@[0,pi]          1st moment: %s ' % mom1 (math.sin) 

    mom1 = CentralMoment  ( 1, 0, math.pi)
    print 'sin@[0,pi] 1st central moment : %s ' % mom1 (math.sin) 

    s    = Skewness ( 0 , math.pi )
    print 'sin@[0,pi]            skewness: %s ' % s    (math.sin)
    
    k    = Kurtosis ( 0 , math.pi )
    print 'sin@[0,pi]            kurtosis: %s ' % k    (math.sin) 

    quan = Quantile ( 0.980 , 0, 10)
    print 'sin@[0,pi]      0.980-quantile: %s ' % quan (math.sin) 


    quan = Quantile ( 0.252 , 0, 10)
    print '1@[0,10]        0.252-quantile: %s ' % quan ( lambda x : 1 ) 

    print '1@[0,1]         0.501-quantile: %s ' % quantile ( lambda x : 1 , 0.501 , 0 , 1 ) 
    print '1@[0,1]         0.201-quantile: %s ' % quantile ( lambda x : 1 , 0.201 , 0 , 1 ) 
    
    print 80*'*'

    from math import exp 
    gau = lambda x : exp(-0.5*x*x)

    print 'CL(gauss,0.68,-10,10)  %s ' % cl_symm ( gau , 0.68 , -10 , 10 )
    print 'CL(gauss,0.68,0,10,0)  %s ' % cl_symm ( gau , 0.68 ,   0 , 10 , x0 = 0 )

    gau1 = lambda x : exp(-0.5*x*x) if x > 0 else 0

    print 'CLa(gauss,0.68,-10,10) (%.3f,%.3f) ' % cl_asymm ( gau  , 0.68 , -10 , 10 )
    print 'CLa(aga  ,0.68,-10,10) (%.3f,%.3f) ' % cl_asymm ( gau1 , 0.68 , -10 , 10 )

    print 'Skewness(gauss,-10,10) %s ' % skewness ( gau  , -10 , 10 )
    print 'Kurtosis(gauss,-10,10) %s ' % kurtosis ( gau  , -10 , 10 )

    print 'Skewness(agau ,-10,10) %s ' % skewness ( gau1 , -10 , 10 )
    print 'Kurtosis(agau ,-10,10) %s ' % kurtosis ( gau1 , -10 , 10 )
    
# =============================================================================
# The END 
# =============================================================================
