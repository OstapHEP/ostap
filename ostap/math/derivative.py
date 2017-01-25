#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  derivative.py 
#  Simple adaptive numerical differentiation (for pyroot/PyRoUts/Ostap)
#  R. De Levie, "An improved numerical approximation for the first derivative"
#  @see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
#
#  .. and also very simple wrapper to numerical integation using scipy
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
#  
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
    "derivative" , ## numerical differentiation (as function)
    "Derivative" , ## numerical differentiation (as object)
    ## 
    "partial"    , ## numerical partial derivatives (as function)
    "Partial"    , ## numerical partial derivatives (as object) 
    ##
    'EvalVE'     , ## evaluate the function taking argument's unicertainty 
    'Eval2VE'    , ## evaluate 2-argument function with argument's unicertainties
    ) 
# =============================================================================
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.deriv' )
else                       : logger = getLogger ( __name__           )
# =============================================================================
from sys import float_info
_eps_   = float_info.epsilon
if not 0.75 < _eps_ * 2**52 < 1.25 :
    import warnings
    warnings.warn ('"epsilon" in not in the expected range!Math could be suboptimal')
# =============================================================================
from ostap.math.base import cpp,iszero,isequal
from ostap.math.ve   import VE 
# =============================================================================
_next_double_ = cpp.Ostap.Math.next_double
_mULPs_       = 1000 
def _delta_ ( x , ulps = _mULPs_ ) :
    n1 = _next_double_ ( x ,  ulps )
    n2 = _next_double_ ( x , -ulps )
    return max ( abs ( n1 - x ) , abs ( n2 - x ) )

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
        if iszero( x )  : h    =             _numbers_[0][I]
        else            : h    = abs ( x ) * _numbers_[I][3] 
        
    ## 1) find the estimate for first and "J"th the derivative with the given step 
    d1,dJ = _dfun_( func , x , h , True )
        
    ## find the optimal step 
    if iszero   ( dJ ) or  ( iszero ( f0 ) and iszero ( x * d1 ) ) :
        if  iszero ( x )    : hopt =             _numbers_[0][I] 
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
#  @code
#  func  = math.sin
#  deriv = Derivative ( func )    
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class Derivative(object) :
    """Calculate the first derivative for the function
    R. De Levie, ``An improved numerical approximation for the first derivative''
    see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
    >>> func = math.sin
    >>> deri = Derivative ( func )        
    """
    # =========================================================================
    ## constructor 
    #  @param func   the function
    #  @param step   proposed initial step for evaluation of derivatives
    #  @param order  derivative is calcualated using 2*I(+1) point
    #  @param err    evaluate numerical uncertainties?
    def __init__ ( self , func , step = 0 , order = 2 , err = False ) :
        """Calculate the first derivative for the function
        R. De Levie, ``An improved numerical approximation for the first derivative''
        see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
        - func:   the function
        - step:   proposed initial step for evaluation of derivatives
        - order  derivative is calcuakted using 2*I+1 point
        - err    evaluate numerical uncertainties?
        >>> func = math.sin
        >>> deri = Derivative ( func )        
        """
        self._func  = func
        self._step  = float( step  ) 
        self._order = int  ( order )
        self._err   = True if err else False 
        if self._order < 0 :
            raise AttributeError("Invalid ``order''-parameter!")

    # =========================================================================
    ## evaluate the derivative
    #  @code 
    #  func  = math.sin
    #  deriv = Derivative ( func )
    #  print deriv(0.1)
    #  @endcode 
    def __call__ ( self , x ) :
        """Calculate the first derivative for the function
        R. De Levie, ``An improved numerical approximation for the first derivative''
        see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
        
        >>> func  = math.sin
        >>> deriv = Derivative ( func )
        
        >>> print deriv(0.1) 
        """
        return derivative ( self._func  ,
                            x           ,
                            self._step  ,
                            self._order ,
                            self._err   )

# =============================================================================
## Calculate the partial derivative for the function
#  @see derivative
#  R. De Levie, "An improved numerical approximation for the first derivative"
#  @see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
#  @code
#  >>> fun2 =  lambda x,y : x*x+y*y
#  >>> print partial ( 0 , fun , (1.0,2.0) ) ) 
#  @endcode
#  @param index  (INPUT) indx of the variable 
#  @param func   (INPUT) the function itself
#  @param x      (INPUT) the argument
#  @param h      (INPUT) the guess for the step used in numeric differentiation
#  @param I      (INPUT) the rule to be used ("N-point rule" = 2*I+1)
#  @param err    (INPUT) calcualte the uncertainty?
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
def partial ( index , func , x , h = 0  , I = 2 , err = False ) : 
    """Calculate the partial derivative for the function
    - index  (INPUT) indx of the variable 
    - func   (INPUT) the function itself
    - x      (INPUT) the argument
    - h      (INPUT) the guess for the step used in numeric differentiation
    - I      (INPUT) the rule to be used ("N-point rule" = 2*I+1)
    - err    (INPUT) calcualte the uncertainty?
    
    Algorithm used:
    R. De Levie, ``An improved numerical approximation for the first derivative''
    - see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
    >>> fun2 =  lambda x,y : x*x+y*y
    >>> print partial ( 0 , fun , (1.0,2.0) ) )     
    - see derivative
    """
    
    if len(x) <= index :
        raise AttributeError("Invalid argument length/index %d/%d" %  ( len(x) , index ) )
    
    _x =   [ float(a) for a in x ]
    
    ## create wrapper function 
    def _wrap ( z ) :
        _z              = _x[index] 
        _x[index] =  z
        _r = func ( *_x )
        _x[index] = _z
        return _r
    
    xi = _x[ index ]
    return derivative ( _wrap , xi , h , I , err )

# =============================================================================
## calcuate the partial derivative for the function
#  @code
#  func = lambda x, y: x * x + y * y  
#  dFdX = Partial ( 0 , func )
#  dFdY = Partial ( 1 , func )
#  x = 1
#  y = 2
#  print ' f(%f,%f)=%f    ' % ( x , y , func( x, y ) ) 
#  print ' dFdX=%f dFdY=%f' % ( dFdX(x,y), dFdY ( x, y ) ) 
#  @endcode 
#  @see Derivative
#  @see partial 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-01-25
class Partial(Derivative) :
    """Calcuate the partial derivative for the function
    >>> func = lambda x, y: x * x + y * y  
    >>> dFdX = Partial( 0 , func )
    >>> dFdY = Partial( 1 , func )
    >>> x = 1
    >>> y = 2
    >>> print ' f(%f,%f)=%f    ' % ( x , y , func( x, y ) ) 
    >>> print ' dFdX=%f dFdY=%f' % ( dFdX(x,y), dFdY ( x, y ) ) 
    """
    # =========================================================================
    ## constructor
    #  @param index index of the variable  (>=0)
    #  @param func  the function
    #  @param step  proposed initial step for derivatives 
    #  @param order derivative is evaluated using 2*I(+1) point  (>=0) 
    #  @param err   estimate the numerical uncertainty?
    def __init__ ( self        ,
                   index       ,   ## index of the variabale 
                   func        ,   ## the function 
                   step  = 0   ,   ## proposed initial step for derivatives
                   order = 2   ,   ## J=2*I(+1) is a number of points used for evaluation of derivative
                   err = False ) : ## estimate the uncertainty?
        """Calculate the partial derivative for the function
        - index  (INPUT) indx of the variable 
        - func   (INPUT) the function itself
        - step   (INPUT) the guess for the step used in numeric differentiation
        - order  (INPUT) the rule to be used ("N-point rule" = 2*I+1)
        - err    (INPUT) calcualte the uncertainty?
        >>> func = lambda x, y: x * x + y * y  
        >>> dFdX = Partial( 0 , func )
        >>> dFdY = Partial( 1 , func )
        >>> x = 1
        >>> y = 2
        >>> print ' f(%f,%f)=%f    ' % ( x , y , func( x, y ) ) 
        >>> print ' dFdX=%f dFdY=%f' % ( dFdX(x,y), dFdY ( x, y ) ) 
        """
        if isinstance ( index , (int,long) ) and 0 <= index :
            self._index = index
        else :
            raise AttributeError("Invalid variable index %s" % index)

        ## keep the function 
        self._func2 = func
        
        ## initialize the base 
        Derivative.__init__ ( self , func , step , order , err )
        
    # =========================================================================
    ## evaluate the derivative
    #  @code 
    #  func  = math.sin
    #  deriv = Derivative ( func )
    #  print deriv(0.1)
    #  @endcode 
    def __call__ ( self , *x ) :
        """Calcuate the partial derivative for the function
        >>> func = lambda x, y: x * x + y * y  
        >>> dFdX = Partial( 0 , func )
        >>> dFdY = Partial( 1 , func )
        >>> x = 1
        >>> y = 2
        >>> print ' f(%f,%f)=%f    ' % ( x , y , func( x, y ) ) 
        >>> print ' dFdX=%f dFdY=%f' % ( dFdX(x,y), dFdY ( x, y ) )
        """
        return partial ( self._index  ,
                         self._func2  ,
                         x            ,
                         self._step   ,
                         self._order  ,
                         self._err    )
        
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
            
        if   name                          : self.name =  name 
        elif hasattr ( func , '__name__' ) and '<lambda>' != func.__name__ :
            self.name = func.__name__
        else                               : self.name = 'Eval2VE'

    ## printout 
    def __str__ ( self ) : return str ( self.name )
    __repr__ = __str__
    
    ## get a value 
    def _value_ ( self , x , *args ) :
        """Evaluate a function"""
        return self._func( float( x ) , *args )

    # =========================================================================
    ## Evaluate the function taking into account uncertainty in the argument
    #  @code
    #  import math 
    #  x    = VE(1,0.1**2)
    #  sin1 = EvalVE( math.sin , lambda s : math.cos(s) )
    #  print 'sin1(x) = %s ' % sin1(x) 
    #  sin2 = EvalVE( math.sin )
    #  print 'sin2(x) = %s ' % sin2(x)
    #  @endcode
    def __call__ ( self , x , *args ) :
        """Evaluate the function taking into account uncertainty in the argument        
        >>> import math 
        >>> x    = VE(1,0.1**2)
        >>> sin1 = EvalVE( math.sin , lambda s : math.cos(s) )
        >>> print 'sin1(x) = %s ' % sin1(x) 
        >>> sin2 = EvalVE( math.sin )
        >>> print 'sin2(x) = %s ' % sin2(x)
        """
        #
        ## evaluate the function 
        val  = self._value_ ( x , *args )
        #
        ## no uncertainties? 
        if   isinstance ( x , ( float , int , long ) )   : return VE ( val , 0 )
        # ignore small or invalid uncertanties 
        elif 0 >= x.cov2() or iszero ( x.cov2()  )       : return VE ( val , 0 )
        # evaluate the derivative  
        d    = self._deriv ( float ( x ) , *args )
        ## calculate the variance 
        cov2 = d * d * x.cov2()
        ## get a final result 
        return VE ( val , cov2 )
    
# ===================================================================================
## @class Eval2VE
#  Evaluate the 2-argument function taking into account the uncertaintines
#  @code
#  func2 = lambda x,y : x*x + y*y
#  eval2 = Eval2VE ( func2 )
#  x = VE(1,0.1**2)
#  y = VE(2,0.1**2)
#  print eval2(x,y)    ## treat x,y as uncorrelated 
#  print eval2(x,y, 0) ## ditto 
#  print eval2(x,y,+1) ## treat x,y as 100% correlated 
#  print eval2(x,y,-1) ## treat x,y as 100% anti-correlated
#  @endcode
#  Partial derivatives can be provided explictely:
#  @code
#  func2 = lambda x,y : x*x + y*y
#  eval2 = Eval2VE ( func2 , dFdX = lambda x,y : 2*x , dFdY = lambda x,y : 2*y ) 
#  @endcode
#  If derivatves are not provided, numerical differentiation will be used 
#  @see EvalVE
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-02-23
class Eval2VE(object) :
    """ Evaluate the 2-argument function taking into account the uncertaintines
    >>> func2 = lambda x,y : x*x + y*y
    >>> eval2 = Eval2VE ( func2 )
    >>> x = VE(1,0.1**2)
    >>> y = VE(2,0.1**2)
    >>> print eval2(x,y)    ## treat x,y as uncorrelated 
    >>> print eval2(x,y, 0) ## ditto 
    >>> print eval2(x,y,+1) ## treat x,y as 100% correlated 
    >>> print eval2(x,y,-1) ## treat x,y as 100% anti-correlated
    Partial derivatives can be provided explictely:
    >>> func2 = lambda x,y : x*x + y*y
    >>> eval2 = Eval2VE ( func2 , dFdX = lambda x,y : 2*x , dFdY = lambda x,y : 2*y )
    If derivatves are not provided, numerical differentiation will be used 
    """
    ## constructor
    #  @param func  the 2-argument function
    #  @param dFdX  (optional) the partial derivative d(dunc)/dX
    #  @param dFdY  (optional) the partial derivative d(dunc)/dY
    #  @param name  (optional) the function name 
    def __init__ ( self , func , dFdX = None , dFdY = None , name = '' ) :
        """ Constructor from the function, optional partial derivative and name
        >>> func2 = lambda x,y : x*x + y*y
        >>> eval2 = Eval2VE ( func2 )
        >>> x = VE(1,0.1**2)
        >>> y = VE(2,0.1**2)
        >>> print eval2(x,y)    ## treat x,y as uncorrelated 
        >>> print eval2(x,y, 0) ## ditto 
        >>> print eval2(x,y,+1) ## treat x,y as 100% correlated 
        >>> print eval2(x,y,-1) ## treat x,y as 100% anti-correlated
        Partial derivatives can be provided explictely:
        >>> func2 = lambda x,y : x*x + y*y
        >>> eval2 = Eval2VE ( func2 , dFdX = lambda x,y : 2*x , dFdY = lambda x,y : 2*y )
        If derivatves are not provided, numerical differentiation will be used 
        """
        self._func  = func
        
        self._dFdX = dFdX if dFdX else Partial ( 0 , func )
        self._dFdY = dFdY if dFdY else Partial ( 1 , func )

        if not hasattr ( self._dFdX , '_step' ) : self._dFdX._step = 0 
        if not hasattr ( self._dFdY , '_step' ) : self._dFdY._step = 0 
            
        if   name                          : self.name =  name 
        elif hasattr ( func , '__name__' ) and '<lambda>' != func.__name__ :
            self.name = func.__name__
        else                               : self.name = 'Eval2VE'

    ## printout 
    def __str__ ( self ) : return str ( self.name )
    __repr__ = __str__
    
    # =========================================================================
    ## get a value 
    def _value_ ( self , x , y ) :
        """Evaluate a function"""
        return self._func( float ( x ) , float(y) )
    
    # =========================================================================
    ## get a partial derivatives,
    #  adjust the step for numerical differentiation: use the initial value of 0.1*error
    def _partial_ ( self , f , x , y , step2 ) :
        """ get a partial derivatives,
        adjust the step for numerical differentiation:
        - use the initial value of 0.1*error
        """
        
        old_step   = f._step
        ## magic number choice: 1.e-8 <= ( step ~ 0.1 * error ) 
        if 1.e-14 <= step2 :
            from math import sqrt as _sqrt_ 
            f._step = 0.1 * _sqrt_( step2 ) 
        _r = f ( x , y )
        f._step = old_step
        return _r 
    
    ## get a partial derivatives d(func)/d(x)
    def dFdX  ( self , x , y ) :
        """Get a partial derivatives d(func)/d(X)"""
        x = VE ( x )
        return self._partial_ ( self._dFdX , x.value() , y          , x.cov2() )
    
    ## get a partial derivatives d(func)/d(y)
    def dFdY  ( self , x , y ) :
        """Get a partial derivatives d(func)/d(Y)"""
        y = VE ( y )
        return self._partial_ ( self._dFdY , x         , y.value()  , y.cov2() )

    # =========================================================================
    ## evaluate the function 
    #  @code
    #  func2 = lambda x,y : x*x + y*y
    #  eval2 = Eval2VE ( func2 )
    #  x = VE(1,0.1**2)
    #  y = VE(2,0.1**2)
    #  print eval2(x,y)    ## treat x,y as uncorrelated 
    #  print eval2(x,y, 0) ## ditto 
    #  print eval2(x,y,+1) ## treat x,y as 100% correlated 
    #  print eval2(x,y,-1) ## treat x,y as 100% anti-correlated
    #  @endcode
    def __call__ ( self , x , y , cxy = 0 ) :
        """Evaluate the function 
        >>> func2 = lambda x,y : x*x + y*y
        >>> eval2 = Eval2VE ( func2 )
        >>> x = VE(1,0.1**2)
        >>> y = VE(2,0.1**2)
        >>> print eval2(x,y)    ## treat x,y as uncorrelated 
        >>> print eval2(x,y, 0) ## ditto 
        >>> print eval2(x,y,+1) ## treat x,y as 100% correlated 
        >>> print eval2(x,y,-1) ## treat x,y as 100% anti-correlated
        """
        ## evaluate the function 
        val = self._value_ ( x , y )
        #
        x   =  VE ( x )
        y   =  VE ( y )
        #
        x_plain = x.cov2() <= 0 or iszero ( x.cov2() )
        y_plain = y.cov2() <= 0 or iszero ( y.cov2() ) 
        #
        if x_plain and y_plain : return VE ( val , 0 )
        #
        ## here we need to calculate the uncertainties
        # 
        cov2 = 0.0
        #
        fx   = self.dFdX ( x         , y.value() ) if not x_plain else 0 
        fy   = self.dFdY ( x.value() , y         ) if not y_plain else 0 
        #
        if not x_plain : cov2 += fx * fx * x.cov2()            
        if not y_plain : cov2 += fy * fy * y.cov2()
        # 
        if not x_plain and not y_plain : 
            ## adjust the correlation coefficient:
            cxy   = min ( max ( -1 , cxy ) , 1 )
            if not iszero ( cxy ) :
                cov2 += 2 * cxy * fx * fy * x.error() * y.error() 

        return VE ( val , cov2 )
 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

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
        
        logger.info ( 80*'*' ) 
        
        for i in range(1,6 ) :
            
            res = derivative ( fun , x , 1.e-20 , i  , err = True )
            f1  = res 
            d   = der(x)  ## the exact value for derivative 
            if iszero ( d ) : 
                logger.info ( 'Rule=%2d %s %s %s' % ( 2*i+1 , f1 , d , (f1.value()-d)   ) ) 
            else      :
                logger.info ( 'Rule=%2d %s %s %s' % ( 2*i+1 , f1 , d , (f1.value()-d)/d ) )  
                    
    logger.info ( 80*'*' ) 

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
        logger.info (  fmt %  ( x       ,
                                func(x) ,
                                1-deriv_1(x)/deriv_a(x) ,
                                1+deriv_2(x)/func   (x) ) ) 
        
    logger.info ( 80*'*' ) 

    ## the function
    func2 = lambda x,y : x*x + y*y
    
    ## use explicit   partial derivatives 
    eval2_1 = Eval2VE( func2 , dFdX = lambda x,y : 2*x , dFdY = lambda x,y : 2*y )
    
    ## use numerical  partial derivatives 
    eval2_2 = Eval2VE( func2 )

    for x,y in [ (0,0) , (1,1) , (1,2) , (2,2) ] :

        x = VE(x,0.1**2)
        y = VE(x,0.1**2)
        
        logger.info ( 'x=%-17s, y=%-17s' % ( x , y ) )
        
        logger.info ( '  eval2_1(x,y)=   %-17s eval2_3(x,y)=   %-17s ' % ( eval2_1 ( x, y ) ,
                                                                           eval2_2 ( x, y ) ) )
        ## use correlation coefficient:
        for c in (0,-1,1) :
            logger.info ( '  eval2_1(x,y,%+2d)=%-17s eval2_3(x,y,%+2d)=%-17s ' % ( c , eval2_1 ( x, y ) ,
                                                                                   c , eval2_2 ( x, y ) ) ) 
    logger.info ( 80*'*' ) 
    
# =============================================================================
# The END 
# =============================================================================
