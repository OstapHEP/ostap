#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
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
    "Derivative"    , ## numerical differentiation (as object) 
    ##
    'EvalVE'        , ## evaluate the function taking argument's unicertainty 
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
        elif 0 >= x.cov2() or iszero ( x.cov2()  )       : return VE ( val , 0 )
        # evaluate the derivative  
        d    = self._deriv ( float ( x ) , *args )
        ## calculate the variance 
        cov2 = d * d * x.cov2()
        ## get a final result 
        return VE ( val , cov2 )
    
    def __str__ ( self ) : return self.__name__
    



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

# =============================================================================
# The END 
# =============================================================================
