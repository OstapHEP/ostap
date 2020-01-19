#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/math/derivative.py 
#  Simple adaptive numerical differentiation (for pyroot/PyRoUts/Ostap)
#  @see R. De Levie, "An improved numerical approximation for the first derivative"
#  @see https://link.springer.com/article/10.1007/s12039-009-0111-y
#  @see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
# =============================================================================
"""Simple adaptive numerical differentiation (for pyroot/PyRoUts/Ostap/...)

>>> func = lambda x : x*x
>>> print derivative ( func , 1 )


>>> func = lambda x : x*x
>>> print integral ( func , 0 , 1 )

there are also object form:

>>> func = ...
>>> deriv = Derivative ( func )
>>> integ = Integral   ( func , 0 )
>>> print deriv ( 0.1 ) , integ ( 0.1 )  

- see R. De Levie, ``An improved numerical approximation for the first derivative''
- see https://link.springer.com/article/10.1007/s12039-009-0111-y
- see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf

"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-06"
__all__     = (
    ##
    "derivative"        , ## numerical differentiation (as function)
    "Derivative"        , ## numerical differentiation (as object)
    ## 
    "partial"           , ## numerical partial derivatives (as function)
    "Partial"           , ## numerical partial derivatives (as object) 
    ##
    'EvalVE'            , ## evaluate the function taking argument's uncertainty 
    'Eval2VE'           , ## evaluate 2-argument function with argument's uncertainties
    'EvalNVE'           , ## evaluate N-argument function with argument's uncertainties
    'EvalNVEcov'        , ## evaluate N-argument function with argument's uncertainties
    'EvalNVEcor'        , ## evaluate N-argument function with argument's uncertainties
    ##
    'complex_derivative', ## evaluiate a complex derivatibe for analytical funtion
    'ComplexDerivative' , ## evaluiate a complex derivatibe for analytical funtion (as object) 
    ) 
# =============================================================================
from   builtins import range 
import ROOT, math
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.derivative' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
from sys import float_info
_eps_   = float_info.epsilon
if not 0.75 < _eps_ * 2**52 < 1.25 :
    import warnings
    warnings.warn ('"epsilon" in not in the expected range!Math could be suboptimal')
# =============================================================================
from ostap.math.base  import cpp , iszero , isequal
from ostap.core.ostap_types import num_types , is_integer
from ostap.math.ve    import VE 
# =============================================================================
_next_double_ = cpp.Ostap.Math.next_double
_mULPs_       = 1000 
def _delta_ ( x , ulps = _mULPs_ ) :
    n1 = _next_double_ ( x ,  ulps )
    n2 = _next_double_ ( x , -ulps )
    return max ( abs ( n1 - x ) , abs ( n2 - x ) )

# =============================================================================
# Four versions:

## use dot_fma from ostap            ## 17.7s
dot_fma = cpp.Ostap.Math.dot_fma  
import array
ARRAY   = lambda x : array.array ( 'd' , x )

# (2) use dot based on Kahan summation  ## 17.5s 
# dot_fma = cpp.Ostap.Math.dot_kahan
# import array 
# ARRAY =  lambda x : array.array ( 'd' , x )
# (3) use numpy variant                 ## 17.9s
# import numpy
# ARRAY   = lambda x     : numpy.array ( x , dtype=float) 
# dot_fma = lambda n,x,y : numpy.dot(x,y)
# (4) SLOWEST:  use math.fsum           ## 21s
# import numpy
# ARRAY   = lambda x     : numpy.array ( x , dtype=float)
# from math import fsum 
# dot_fma = lambda n,x,y : fsum ( ( (i*j) for i,j in zip(x,y) ) )


# ======================================================================================
## calculate  1st (and optionally Nth) derivative with the given step
#  - f'      is calcualted as O(h^(N-1))
#  - f^{(N)} is calcualted as O(h^2)
#  Simple adaptive numerical differentiation
#  @see R. De Levie, "An improved numerical approximation for the first derivative"
#  @see https://link.springer.com/article/10.1007/s12039-009-0111-y
#  @see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
class DeLevie(object) :
    """Calculate  1st (and optionally Nth) derivative with the given step
    - f'      is calcualted as O(h^(N-1))
    - f^{(N)} is calcualted as O(h^2) 
    Simple adaptive numerical differentiation 
    - see R. De Levie, ``An improved numerical approximation for the first derivative''
    - see https://link.springer.com/article/10.1007/s12039-009-0111-y
    - see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
    """
    __d1h = [
        ( ARRAY ( [      +1 ,       0                                                      ] ) ,      2 ) ,
        ( ARRAY ( [      +8 ,      -1 ,      0                                             ] ) ,     12 ) ,       
        ( ARRAY ( [      45 ,      -9 ,     +1 ,      0                                    ] ) ,     60 ) ,
        ( ARRAY ( [    +672 ,    -168 ,    +32 ,     -3 ,     0                            ] ) ,    840 ) ,
        ( ARRAY ( [   +2100 ,    -600 ,   +150 ,    -25 ,    +2 ,     0                    ] ) ,   2520 ) ,
        ( ARRAY ( [  +23760 ,   -7425 ,  +2200 ,   -495 ,   +72 ,    -5 ,    0             ] ) ,  27720 ) ,
        ( ARRAY ( [ +315315 , -105105 , +35035 ,  -9555 , +1911 ,  -245 ,  +15 ,   0       ] ) , 360360 ) ,
        ( ARRAY ( [ +640640 , -224224 , +81536 , -25480 , +6272 , -1120 , +128 ,  -7 ,  0  ] ) , 720720 ) ] 

    __d2h = [
        ( ARRAY ( [      -2 ,      +1                                                      ] ) , 2      ) , 
        ( ARRAY ( [      +5 ,      -4 ,     +1                                             ] ) , 2      ) , 
        ( ARRAY ( [     -14 ,     +14 ,     -6 ,     +1                                    ] ) , 2      ) , 
        ( ARRAY ( [     +42 ,     -48 ,    -27 ,     -8 ,    +1                            ] ) , 2      ) , 
        ( ARRAY ( [    -132 ,    +165 ,   -110 ,    +44 ,   -10 ,    +1                    ] ) , 2      ) , 
        ( ARRAY ( [    +429 ,    -572 ,   +429 ,   -208 ,   +65 ,   -12 ,   +1             ] ) , 2      ) , 
        ( ARRAY ( [   -1430 ,   +2002 ,  -1638 ,   +910 ,  -350 ,   +90 ,  -14 ,  +1       ] ) , 2      ) , 
        ( ARRAY ( [   +4862 ,   -7072 ,  +6188 ,  -3808 , +1700 ,  -544 , +119 , -16 , +1  ] ) , 2      ) ]

    ## constructor with the order parameter
    def __init__ ( self , o ) :

        assert is_integer ( o ), 'Invalid order %s' % o
        
        self.__order = o
        if   o < 0                : self.__order = 0
        elif o >= len(self.__d1h) : self.__order = len(self.__d1h) - 1 

        self.__d1          = self.__d1h[ self.__order ][0]
        self.__sf1         = self.__d1h[ self.__order ][1] 
        self.__d2          = self.__d2h[ self.__order ][0]
        self.__sf2         = self.__d2h[ self.__order ][1]
        
        ## vector of function differences 
        self.__df          = ARRAY ( ( self.__order + 2 ) * [ 0 ] ) 

    @property
    def order ( self ) :
        """Order of differentiation: 2*order+1 points witll be used"""
        return self.__order
        
    # =========================================================================
    ## calculate  1st (and optionally Nth) derivative with the given step
    #  - f'      is calcualted as O(h^(N-1))
    #  - f^{(N)} is calcualted as O(h^2)
    #  Simple adaptive numerical differentiation
    #  @see R. De Levie, "An improved numerical approximation for the first derivative"
    #  @see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
    def __call__ ( self , func , x , h , der = False ) :
        """Calculate  1st (and optionally Nth) derivative with the given step
        - f'      is calcualted as O(h^(N-1))
        - f^{(N)} is calcualted as O(h^2)
        Simple adaptive numerical differentiation
        - see R. De Levie, ``An improved numerical approximation for the first derivative''
        - see https://link.springer.com/article/10.1007/s12039-009-0111-y
        - see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
        """

        ## calculate differences 
        imax = self.__order + 2 if der else self.__order + 1
        i = 0
        while i < imax : 
            j = i + 1
            self.__df[i] = func ( x + j * h ) - func ( x - j * h )
            i += 1
            
        ## 1) calculate 1st derivative 
        result = dot_fma ( self.__order + 1 , self.__df , self.__d1 ) / ( self.__sf1 * h )            
        if not der : return result 
            
        ## 2) calculate Nth derivative 
        dd     = dot_fma ( self.__order + 2 , self.__df , self.__d2 ) / ( self.__sf2 * h**(self.__order*2+3) ) 
        
        return result, dd 
                            
# =============================================================================
# The actual setup for
# =============================================================================
_funcs_ = (
    DeLevie(0) , ## 0 == 1 
    DeLevie(0) , ## 1 
    DeLevie(1) , ## 2 
    DeLevie(2) , ## 3
    DeLevie(3) , ## 4 
    DeLevie(4) , ## 5 
    DeLevie(5) , ## 6
    DeLevie(6) , ## 7
    DeLevie(7) , ## 8
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
    (    140 , 6.9349e-17 , 7.7091e-3 , 5.8046e-3 ) , ## I=3, J= 7  7-point rule 
    (    630 , 7.4832e-17 , 2.6237e-2 , 1.8227e-2 ) , ## I=4, J= 9  9-point rule 
    (   2772 , 7.8754e-17 , 5.7292e-2 , 3.7753e-2 ) , ## I=5, J=11 11-point rule 
    (  12012 , 8.1738e-17 , 9.8468e-2 , 6.2500e-2 ) , ## I=6, J=13 13-point rule 
    (  51480 , 8.4108e-17 , 1.4656e-1 , 9.0454e-2 ) , ## I=7, J=15 15-point rule 
    ( 218790 , 8.6047e-17 , 1.9873e-1 , 1.2000e-1 ) , ## I=8, J=17 17-point rule 
    )

# =============================================================================
## Calculate the first derivative for the function
#  @see R. De Levie, "An improved numerical approximation for the first derivative"
#  @see https://link.springer.com/article/10.1007/s12039-009-0111-y
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
    - see R. De Levie, ``An improved numerical approximation for the first derivative''
    - see https://link.springer.com/article/10.1007/s12039-009-0111-y
    - see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
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

    h = max ( h , 2 * delta )
    
    ## 1) find the estimate for first and "J"th the derivative with the given step 
    d1 , dJ = _dfun_( func , x , h , True )
        
    ## find the optimal step 
    if iszero   ( dJ ) or  ( iszero ( f0 ) and iszero ( x * d1 ) ) :
        if  iszero ( x )    : hopt =             _numbers_[0][I] 
        else                : hopt = abs ( x ) * _numbers_[I][3]
    else : 
        hopt = _numbers_[I][2] * ( ( abs ( f0 ) + abs ( x * d1 ) ) / abs ( dJ ) )**( 1.0 / J )

    ## finally get the derivative 
    if not err  :  return _dfun_ ( func , x , hopt , False )

    ## estimate the uncertainty, if needed  
    d1,dJ =  _dfun_ ( func , x , hopt , True )
    
    e     =  _numbers_[I][1] / _numbers_[I][2] * J / ( J - 1 ) 
    e2    =  e * e * ( J * _eps_ + abs ( f0 ) + abs( x * d1 ) )**( 2 - 2./J ) * abs( dJ )**(2./J) 
    return VE ( d1 , 4 * e2 ) 

# =============================================================================
## @class Derivative
#  Calculate the first derivative for the function
#  @see R. De Levie, "An improved numerical approximation for the first derivative"
#  @see https://link.springer.com/article/10.1007/s12039-009-0111-y
#  @see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
#  @code
#  func  = math.sin
#  deriv = Derivative ( func )    
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class Derivative(object) :
    """Calculate the first derivative for the function
    - see R. De Levie, ``An improved numerical approximation for the first derivative''
    - see https://link.springer.com/article/10.1007/s12039-009-0111-y
    - see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
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
        - order  derivative is calculated using 2*I+1 point
        - err    evaluate numerical uncertainties?
        >>> func = math.sin
        >>> deri = Derivative ( func )        
        """
        assert is_integer ( order ) and 0 <= order, "Invalid ``order''-parameter!"
        
        self.__func  = func
        self.__step  = float( step  ) 
        self.__order = int  ( order )
        self.__err   = True if err else False 

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
        return derivative ( self.__func  ,
                            x           ,
                            self.__step  ,
                            self.__order ,
                            self.__err   )

    @property
    def func ( self ) :
        """The function to be differentiated"""
        return self.__func
    @property
    def step ( self ) :
        """proposed initial step for evaluation of derivatives"""
        return self.__step
    @property
    def order ( self ) :
        """ Derivative is calculated using 2*order+1 point"""
        return self.__order
    @property
    def err   (  self ) :
        """Evaluate numerical uncertainties?"""
        return self.__err

        
# =============================================================================
## Calculate the partial derivative for the function
#  @see R. De Levie, "An improved numerical approximation for the first derivative"
#  @see https://link.springer.com/article/10.1007/s12039-009-0111-y
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
    - see R. De Levie, ``An improved numerical approximation for the first derivative''
    - see https://link.springer.com/article/10.1007/s12039-009-0111-y
    - see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
    >>> fun2 =  lambda x,y : x*x+y*y
    >>> print partial ( 0 , fun , (1.0,2.0) ) )     
    - see derivative
    """
    
    if len(x) <= index :
        raise AttributeError("Invalid argument length/index %d/%d" %  ( len(x) , index ) )
    
    _x =  [ float(a) for a in x ]
    
    ## create wrapper function 
    def _wrap ( z ) :
        _z        = _x[index] 
        _x[index] =  z
        _r = func ( *_x )
        _x[index] = _z
        return _r
    
    x_i = _x[ index ]
    return derivative ( _wrap , x = x_i , h = h , I = I , err = err )

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
        - index  (INPUT) index of the variable 
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
        assert is_integer ( index ) and 0 <= index, \
               "Invalid variable index %s" % index 

        ## get the index 
        self.__index = index 
        
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
        return partial ( self.__index ,
                         self.func    ,
                         x            ,
                         self.step    ,
                         self.order   ,
                         self.err     )

    @property 
    def index ( self ) :
        "Index of the variable to be differentiated"
        return self.__index 
        
# =============================================================================
## @class EvalVE
#  Evaluate the function taking into account the uncertainty in the argument
#  @code
#  import math 
#  x = VE(1,0.1**2)
#  sin = EvalVE( math.sin , lambda s : math.cos(s) )
#  print 'sin(x) = %s ' % sin(x) 
#  @endcode
#  @see ostap.math.math_ve 
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
    see also ostap.math.math_ve
    """
    ## constructor
    def __init__ ( self , func , deriv = None  , name = '' ) :
        """ Constructor from the function, derivative and name
        >>> sin1 = EvalVE( math.sin , lambda s : math.cos(s) )
        >>> sin2 = EvalVE( math.sin )  ## numerical derivative will be used 
        """
        self.__func  = func
        if    deriv and callable ( deriv ) : self.__derivative = deriv
        else                               : self.__derivative = Derivative(func)
            
        if   name                          : self.__name__ =  name 
        elif hasattr ( func , '__name__' ) and '<lambda>' != func.__name__ :
            self.__name__ = func.__name__
        else                               : self.__name__ = 'Eval2VE'

    ## printout 
    def __str__ ( self ) : return str ( self.__name__ )
    __repr__ = __str__

    ## get a value 
    def func_eval ( self , x , *args ) :
        """Evaluate a function"""
        return self.__func( float ( x ) , *args )

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
        ## 1) evaluate the function 
        val  = self.func_eval ( x , *args )
        #
        ## no uncertainties? 
        if   isinstance ( x , num_types          ) : return VE ( val , 0 )
        # ignore small or invalid uncertanties 
        elif 0 >= x.cov2() or iszero ( x.cov2()  ) : return VE ( val , 0 )
        #
        ## 2) evaluate the derivative
        dfun = self.__derivative
        d    = dfun (  float ( x ) , *args ) 
        ## 3) calculate the variance 
        cov2 = d * d * x.cov2()
        ## 4) get a final result 
        return VE ( val , cov2 )

    @property
    def func       ( self ) :
        """Get the function itself"""
        return self.__func 
    @property
    def derivative ( self ):
        """Get the derivative object/function"""
        return self.__derivative
        
# ===================================================================================
## @class EvalNVE
#  Calculate the value of scalar function of N (non-correlated) arguments,
#  taking into account the argument uncertainties.
#  @code
#  fun2   = lambda x , y : x**2+y**2 
#  fun2e  = EvalNVE ( 2 , func )
#  x , y  = ...
#  result = fun2e (  
#  @endcode
class  EvalNVE(object) :
    """Calculate the value of scalar function of N-arguments,
    taking into account the argument uncertainties
    >>> fun2   = lambda x , y : x**2+y**2 
    >>> fun2e  = EvalNVE( 2 , func )e
    >>> x , y  = ...
    >>> result = fun2e ( x , y )  
    """
    # =========================================================================
    ## create the object
    #  @param N       the dimensionality of the problem
    #  @param func    the function
    #  @param partial optional vector of partial derivatives 
    def __init__ ( self , N , func , partial = () , name = '' ) :
        """
        """
        assert is_integer ( N )  and 1 < N , 'Invalid "N"=%s'       % N
        assert callable   ( func )         , 'Invalid "func" %s/%s' % ( func , type ( func ) )
        assert ( not partial ) or len ( partial ) == N, 'Invalid "partial"'
        
        self.__func = func
        self.__N    = N
        
        derivatives = []
        for i in  range ( N ) :
            
            if i < len ( partial ) and callable ( partial[i] ) : d_i = partial[i] 
            else                                               : d_i = Partial( i , func )
                
            derivatives.append ( d_i )
            
        self.__partial = tuple (  derivatives )

        if   name  : self.__name__  =  name 
        elif hasattr ( func , '__name__' ) and '<lambda>' != func.__name__ :
            self.__name__ = func.__name__
        else       : self.__name__ = 'EvalNVE'

    @property 
    def func    ( self ) :
        """The original function"""
        return  self.__func
    @property
    def partial ( self ) :
        """Get tuple of partial derivatives (functions)"""
        return self.__partial
    @property
    def N       ( self ) :
        """N -    number of arguments"""
        return self.__N
    
    ## printout 
    def __str__ ( self ) : return str ( self.__name__ )
    __repr__ = __str__

    ## calculate the gradient
    def gradient ( self , *x ) :
        n = self.__N
        assert len ( x ) == n , 'Invalid argument size'

        xx = [ VE(i).value()        for i in x           ]
        gr = [ self.partial[i](*xx) for i in range ( n ) ]
        
        return tuple ( gr ) 

    ## the main method 
    def __call__ ( self , *x ) :

        n = self.__N 
        assert len ( x ) == n , 'Invalid  argument size'

        xve = [ VE(i)        for i in  x   ] ## force everything to be VE 
        xv  = [    i.value() for i in  xve ] ## get only the values
        
        xv  = tuple ( xv ) 
        
        ## value of the function 
        value = self.__func ( *xv )
        
        ## calculate the covariance 
        cov2  = 0
        for i in range ( n ) :

            c2 = xve[i].cov2()
            if c2 <= 0 or iszero ( c2 ) : continue

            ## calculate the gradient 
            df = self.partial[i] ( *xv ) 
            if iszero ( df )            : continue 
            
            cov2 += c2 * df**2 

        return VE ( value , cov2 ) 
                

                   
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
class Eval2VE(EvalNVE) :
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

        EvalNVE.__init__ ( self , 2 ,  func , ( dFdX , dFdY ) , name )
        
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
        assert isinstance ( cxy , num_types ) and \
               ( abs ( cxy ) <= 1 or isequal ( abs ( cxy ) , 1 ) ) , \
               'Invalid correlation coefficient %s' % cxy 
        
        x   = VE ( x ) 
        y   = VE ( y )
        
        xv  = x.value()
        yv  = y.value()
        
        val = self.func ( xv , yv )

        xc2 = x.cov2()
        yc2 = x.cov2()
        
        x_plain = xc2 <= 0 or iszero ( xc2 )
        y_plain = yc2 <= 0 or iszero ( yc2 )
        
        #
        if x_plain and y_plain : return VE ( val , 0 )
        
        #
        ## here we need to calculate the uncertainties
        # 
        dx   = self.partial[0] ( xv , yv ) if not x_plain else 0 
        dy   = self.partial[1] ( xv , yv ) if not y_plain else 0 
        #
        
        cov2 = dx * dx * xc2 + dy * dy * yc2
        
        if cxy and xc2 and yc2 :
            cov2 += 2 * cxy * dx * dy * math.sqrt ( xc2 * yc2 ) 
            
        return VE ( val , cov2 )
    
# =============================================================================


# =============================================================================
# @class EvalNVEcov
# Calcualte the value of the scalar function of N (scalar) arguments
# with the given covariance matrix
# @code
# fun2   = lambda x , y : ...
# fun2e  = EvalNVEcov ( fun2 , 2 )
# cov2   = ... # covariance matrix:
# result = fun2e ( (1,5) , cov2 ) ##   
# @endcode 
class EvalNVEcov(EvalNVE) :
    """Calcualte the value of the scalar function of N (scalar) arguments 
    with the given NxN covariance matrix
    
    >>> fun2   = lambda x , y : ...
    >>> fun2e  = EvalNVEcov ( fun2 , 2 )
    >>> cov2   = ... # 2x2 symmetric covariance matrix
    >>> result = fun2e ( (1,5) , cov2 ) ##   

    """
    def __init__ ( self , func , N , partial = () , name = '' , cov2getter = None ) :
        
        EvalNVE.__init__ ( self , func = func , N = N , partial = partial , name = name  )
        if cov2getter is None :
            from ostap.math.linalg import  matrix as get_cov2
            cov2getter = get_cov2 
        self.__cov2 = cov2getter 

    def __call__ ( self , args , cov2 = None ) :
        """Calcualte the value of the scalar function of N (scalar) arguments 
        with the given NxN covariance matrix
        
        >>> fun2   = lambda x , y : ...
        >>> fun2e  = EvalNVEcov ( fun2 , 2 )
        >>> cov2   = ... # 2x2 symmetric covariance matrix
        >>> result = fun2e ( (1,5) , cov2 ) ##   
        
        """
        
        n = self.N 
        assert len ( args ) == n , 'Invalid argument size'

        ## get value of the function 
        val = float ( self.func ( *args ) ) 
        
        ## no covariance matrix is specified ?
        if not cov2 : return val
        
        c2 = 0
        g  = n * [ 0 ] ## gradient 
        
        for i in range ( n ) :

            di    = self.partial[i] ( *args )
            if iszero ( di ) : continue
            
            g [i] = di  

            cii   = self.__cov2 ( cov2 , i , i )
            c2   += di * di * cii
            
            for j in range ( i ) : 

                dj   = g [ j ]
                if iszero ( dj ) : continue
                
                cij  =     self.__cov2 ( cov2 , i , j ) 
                c2  += 2 * di * dj * cij 

        return VE ( val , c2 ) 


# =============================================================================
# @class EvalNVEcor
# Calcualte the value of the scalar function of N arguments
# with uncrtainties and correlations 
# @code
# fun2   = lambda x , y : ...
# fun2e  = EvalNVEcor ( fun2 , 2 )
# cor    = ... # correlation matrix
# x      = ... # 
# y      = ... 
# result = fun2e ( ( x,y) , cor ) ##   
# @endcode 
class EvalNVEcor(EvalNVEcov) :
    """Calcualte the value of the scalar function of N arguments
    with uncertainties and correlations 
    
    >>> fun2   = lambda x , y : ...
    >>> fun2e  = EvalNVEcor ( fun2 , 2 )
    >>> cor    = ... # 2x2 symmetric correlation matrix
    >>> x      = ...
    >>> y      = ...    
    >>> result = fun2e ( (x,y) , cor ) ##   

    """
    def __init__ ( self , func , N , partial = () , name = '' , corrgetter = None  ) :
        
        EvalNVEcov.__init__ ( self , func =  func , N = N , partial = partial , name = name  , cov2getter = corrgetter )
        if corrgetter is None :
            from ostap.math.linalg import matrix as get_corr
            corrgetter = get_corr 
        self.__corr = corrgetter 

    def __call__ ( self , args , cor = None ) :
        """Calcualte the value of the scalar function of N arguments
        with uncertainties and correlations 
        
        >>> fun2   = lambda x , y : ...
        >>> fun2e  = EvalNVEcor ( fun2 , 2 )
        >>> cor    = ... # 2x2 symmetric correlation matrix
        >>> x      = ...
        >>> y      = ...    
        >>> result = fun2e ( (x,y) , cor ) ##           
        """
        n = self.N 
        assert len ( args ) == n , 'Invalid argument size'

        ## get value of the function 
        val = self.func ( *args )
        
        c2  = 0
        
        x   = n * [ 0 ]  ## argument 
        g   = n * [ 0 ]  ## gradient 
        
        for i in range ( n ) :

            xi    = VE ( args[i] )
            x [i] = x 
            
            ci    = xi.cov2()
            if ci < 0  or iszero ( ci ) : continue
            
            di    = self.partial[i] ( *args )
            if iszero ( di )            : continue
            
            ei    = xi.error() 

            g [i] = di  
            e [i] = ei
            
            ## diagonal correlation coefficients are assumed to be 1 and ignored! 
            c2   += ci * di * di
            
            for j in range ( i ) : 

                xj   = x  [ j ]
                cj   = xj.cov2 () 
                if cj < 0  or iszero ( cj ) : continue
                dj   = d  [ j ]
                if iszero ( dj )            : continue                
                ej   = e  [ j ]

                rij  =  self.__corr ( cor , i , j ) if cor else 0 
                assert -1 <= rij <= 1 or isequal ( abs ( rij ) , 1 ) ,\
                       'Invalid correlaation coefficient (%d,%d)=%s ' % ( i , j , rij )
                
                c2  += 2.0 * di * dj * rij * ei * ej 
                
        return VE ( val , c2 ) 


# =============================================================================
## Calcualte the  derivative from analytical function \f$ \frac{{\mathrm{d}} f }{{\mathrm{d}} z }\f$.
#  The function is assumed to be analytical (Cauchy-Riemann conditions are valid).
#  The derivative of the function
#  \f$  f(z) \equiv u(x,y) + iv(x,u) \f$
#  is calculated as 
#  \f$ \frac{{\mathrm{d}} f}{{\mathrm{d}} z} \equiv \frac{{\mathrm{d}} u}{{\mathrm{d}} x} + i\frac{{\mathrm{d}} v}{{\mathrm{d}} x}\f$
#  @code
#  fun = ...
#  d   = complex_derivative ( fun , x = 1+2j )
#  @endcode
#  @param fun  (INPUT) the analytical function
#  @param z    (INPUT) the complex point
#  @param h    (INPUT) the guess for the step used in numeric differentiation
#  @param I    (INPUT) the rule to be used ("N-point rule" = 2*I+1)
#  @param err  (INPUT) calcualte the uncertainty?
#  @return the derivative d(fun)/dz  at point <code>z</code>
#  
def complex_derivative ( fun , z , h = 0 , I = 3 , err = False ,  real = True , imag = True ) :
    """Get a complex derivative
    - The function is assumed to be analytical (Cauchy-Riemann conditions are valid).
    >>> fun = ...
    >>> d   = complex_derivative ( fun , x = 1+2j ) 
    """
    
    Z = complex ( z )
    
    X = Z.real
    Y = Z.imag

    ## few altenatives to calculate the real and imaginary part
    
    if  real :
        UX =  lambda x : complex ( fun ( complex ( x , Y ) ) ).real
        ## Real part 
        re =  derivative ( UX , X , h = h , I = I , err = err )
    else :
        VY =  lambda y : complex ( fun ( complex ( X , y ) ) ).imag             
        ## Real part 
        re =  derivative ( VY , Y , h = h , I = I , err = err )

    if imag :  
        VX =  lambda x : complex ( fun ( complex ( x , Y ) ) ).imag        
        ## Imaginary part 
        im =  derivative ( VX , X , h = h , I = I , err = err )
    else :
        UY =  lambda y : complex ( fun ( complex ( X , y ) ) ).real
        ## Imaginary part 
        im = -derivative ( UY , Y , h = h , I = I , err = err )
    
    if not err : return complex ( re , im )
    
    result = complex ( re.value() , im.value() )
    error  = ( re.cov2() +  im.cov2() ) ** 0.5 
    
    return result , error 


# =============================================================================
## @class ComplexDerivative
#  Calculate the first derivative for the function
#  @see R. De Levie, "An improved numerical approximation for the first derivative"
#  @see https://link.springer.com/article/10.1007/s12039-009-0111-y
#  @see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
#  @code
#  func  = cmath.sin
#  deriv = ComplexDerivative ( func )    
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class ComplexDerivative(Derivative) :
    """Calculate the first derivative for the function
    - see R. De Levie, ``An improved numerical approximation for the first derivative''
    - see https://link.springer.com/article/10.1007/s12039-009-0111-y
    - see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
    >>> func = cmath.sin
    >>> deri = ComplexDerivative ( func )        
    """
    def __init__ ( self , func   ,
                   real  = True  ,
                   imag  = True  ,
                   step  = 0     ,
                   order = 2     ,
                   err   = False ) :
        
        Derivative.__init__ ( self , func = func , step = step , order = order , err = err )
        
        self.__real = True if real else False 
        self.__imag = True if imag else False 
               
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
        
        >>> func  = cmath.sin
        >>> deriv = ComplexDerivative ( func )
        
        >>> print deriv(0.1) 
        """
        return complex_derivative ( self.func          ,
                                    complex ( x )      ,
                                    h    = self.step   ,
                                    I    = self.order  ,
                                    err  = self.err    ,
                                    real = self.__real ,
                                    imag = self.__imag )
    
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
    deriv_1 = Derivative ( func    , order = 5 )  
    ## numerical second derivative     
    deriv_2 = Derivative ( deriv_1 , order = 5 )  

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
