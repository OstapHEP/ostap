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
    "derivative"         , ## numerical differentiation (as function)
    "Derivative"         , ## numerical differentiation (as object)
    ##
    'partial_derivative' , ## numerical partial derivatives (as function)
    "PartialDerivative"  , ## numerical partial derivatives (as object) 
    ##
    'complex_derivative' , ## evaluiate a complex derivatibe for analytical funtion
    'ComplexDerivative'  , ## evaluiate a complex derivatibe for analytical funtion (as object) 
    ## 
    'complex_derivative' , ## evaluate a complex derivatibe for analytical funtion
    'ComplexDerivative'  , ## evaluate a complex derivatibe for analytical funtion (as object) 
    ) 
# =============================================================================
from   ostap.math.base        import Ostap , iszero  , isequal
from   ostap.math.ve          import VE 
from   ostap.math.finitediffs import Rule  , the_dot , darray , delta 
from   ostap.utils.utils      import classprop 
import sys 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.delevie' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
epsilon = sys.float_info.epsilon
if not 0.75 < epsilon * 2**52 < 1.25 :
    logger.warning ('"epsilon" in not in the expected range! Math could be suboptimal')
    
# ======================================================================================
## The rules for numerical differentiation (1st derivative) 
#  - It uses central differences
#  - In addition it allows to avaluate Nth derivative 
#  @see R. De Levie, "An improved numerical approximation for the first derivative"
#  @see https://link.springer.com/article/10.1007/s12039-009-0111-y
#  @see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
class DeLevie(Rule) :
    """The intermediate base class for numerical differentiation
    - see R. De Levie, ``An improved numerical approximation for the first derivative''
    - see https://link.springer.com/article/10.1007/s12039-009-0111-y
    - see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
    """
    ## central differences for the 1st derivative 
    __cd_D1 = (
        ( darray ( [      +1 ,       0                                                      ] ) ,      2 ) , 
        ( darray ( [      +8 ,      -1 ,      0                                             ] ) ,     12 ) ,       
        ( darray ( [      45 ,      -9 ,     +1 ,      0                                    ] ) ,     60 ) ,       
        ( darray ( [    +672 ,    -168 ,    +32 ,     -3 ,     0                            ] ) ,    840 ) ,       
        ( darray ( [   +2100 ,    -600 ,   +150 ,    -25 ,    +2 ,     0                    ] ) ,   2520 ) ,       
        ( darray ( [  +23760 ,   -7425 ,  +2200 ,   -495 ,   +72 ,    -5 ,    0             ] ) ,  27720 ) ,       
        ( darray ( [ +315315 , -105105 , +35035 ,  -9555 , +1911 ,  -245 ,  +15 ,   0       ] ) , 360360 ) ,       
        ( darray ( [ +640640 , -224224 , +81536 , -25480 , +6272 , -1120 , +128 ,  -7 ,  0  ] ) , 720720 ) ) 

    ## central differences for the Nth derivative 
    __cd_DN = (
        ( darray ( [      -2 ,      +1                                                      ] ) , 2      ) , 
        ( darray ( [      +5 ,      -4 ,     +1                                             ] ) , 2      ) , 
        ( darray ( [     -14 ,     +14 ,     -6 ,     +1                                    ] ) , 2      ) , 
        ( darray ( [     +42 ,     -48 ,    -27 ,     -8 ,    +1                            ] ) , 2      ) , 
        ( darray ( [    -132 ,    +165 ,   -110 ,    +44 ,   -10 ,    +1                    ] ) , 2      ) , 
        ( darray ( [    +429 ,    -572 ,   +429 ,   -208 ,   +65 ,   -12 ,   +1             ] ) , 2      ) , 
        ( darray ( [   -1430 ,   +2002 ,  -1638 ,   +910 ,  -350 ,   +90 ,  -14 ,  +1       ] ) , 2      ) , 
        ( darray ( [   +4862 ,   -7072 ,  +6188 ,  -3808 , +1700 ,  -544 , +119 , -16 , +1  ] ) , 2      ) )

    
    ## choice for the optimal step:
    __deltas  = (
        0.5*10**(-10./ 3) ,  ## I=0 
        0.5*10**(-10./ 5) ,  ## I=1 
        0.5*10**(-10./ 7) ,  ## I=2
        0.5*10**(-10./ 9) ,  ## I=3 
        0.5*10**(-10./11) ,  ## I=4
        0.5*10**(-10./13) ,  ## I=5
        0.5*10**(-10./15) ,  ## I=6 
        0.5*10**(-10./17) )  ## I=7
    
    ## choice for the optimal step: Table 2 from DeLevie's paper 
    __table2 = (
        ##   1/bj  cj           dj          eps^(1/j) 
        (      6 , 4.5324e-17 , 5.1422e-6 , 6.0554e-6 ) , ## I=0, J= 3  3-point rule 
        (     30 , 6.0903e-17 , 8.5495e-4 , 7.4009e-4 ) , ## I=1, J= 5  5-point rule 
        (    140 , 6.9349e-17 , 7.7091e-3 , 5.8046e-3 ) , ## I=2, J= 7  7-point rule 
        (    630 , 7.4832e-17 , 2.6237e-2 , 1.8227e-2 ) , ## I=3, J= 9  9-point rule 
        (   2772 , 7.8754e-17 , 5.7292e-2 , 3.7753e-2 ) , ## I=4, J=11 11-point rule 
        (  12012 , 8.1738e-17 , 9.8468e-2 , 6.2500e-2 ) , ## I=5, J=13 13-point rule 
        (  51480 , 8.4108e-17 , 1.4656e-1 , 9.0454e-2 ) , ## I=6, J=15 15-point rule 
        ( 218790 , 8.6047e-17 , 1.9873e-1 , 1.2000e-1 ) , ## I=7, J=17 17-point rule 
        )

    @classprop 
    def IMAX ( cls ) :
        """``IMAX'' : maximal rule index"""
        return 7  

    # =========================================================================
    ## standard constructor with I-parameter
    #  number of function points to be evaluated is 
    def __init__ ( self , I , with_error = False , max_step = -1 ) :
        
        assert isinstance ( I , int ) and 0 <= I , "Invalid value of ``I''-parameter"
        
        if DeLevie.IMAX < I :
            logger.warning ("DeLevie: I/%d exceeds IMAX/%d" % ( I , DeLevie.IMAX ) )
            I = DeLevie.IMAX 

        J = I * 2 + 3
        
        ## initialize the base class 
        Rule.__init__ ( self                              ,  
                        D          =  1                   ,
                        N          = J - 1                ,
                        interval   =  ( - I - 1 , I + 1 ) ,
                        with_error = with_error           ,
                        max_step   = max_step             ) 
        
        self.__I    = I 
        self.__J    = J 
        
        self.__df   = darray ( ( self.__I + 2 ) * [ 0.0 ] )
        
    @property
    def I ( self ) :
        """``I''-parameter: actual function is evaluated in 2*I+2+(+1) point for the 1st derivative
        """
        return self.__I
    
    @property
    def J ( self ) :
        """``J''-parameter : number of points, len of stencil. 
        - derivatives are calculated up with precision f^(J) delta^(J-1)
        """
        return self.__J

    @property
    def cd_D1 ( self ) :
        """``cd_D1'' : central differences configuration for the 1st derivative"""
        return self.__cd_D1 [ self.__I ] 
    
    @property
    def cd_DN ( self ) :
        """``cd_DN'' : central differences configuration for the Nth derivative"""
        return self.__cd_DN [ self.__I ] 

    @property
    def table2 ( self ) :
        """``table2'' : Table2 from DeLevie's paper"""
        return self.__table2 [ self.__I ]  

    @property
    def deltas ( self ) :
        """``deltas'' : step sise for """
        return self.__deltas

    # =========================================================================
    ## The main method of the rule: evaluate the numerical derivative with the given step 
    #  @param func   the function <code>func(x,*args,**kwargs)</code>
    #  @param x      point where derivative is evaluated
    #  @param h      step size 
    #  @param args   addtional positional arguments
    #  @param kwargs addtional kewword arguments
    def __call__ ( self , func , x , h , args = () , kwargs = {}) :
        """The main method: evaluate the numerical derivative 
        - func   the function with signature `func(x,*args,**kwargs)`
        - x      point where derivative is evaluated
        - h      step size 
        - args   addtional positional arguments
        - kwargs addtional keyword arguments
        """
        
        if not self.with_error : 
            return self.derivatives ( False , func , x , h , args = args , kwargs = kwargs )

        ## estimate the uncertainty, if needed  
        d1 , dJ =  self.derivatives ( True  , fun  , x , h , args = args , kwargs = kwargs  )
        
        ## get the function value at the given point 
        f0    = fun ( x , *args , **kwargs )
        
        i     = self.__I 
        j     = self.__J

        c_j   = self.table2 [ 1 ]
        d_j   = self.table2 [ 2 ]
        
        e     = j * c_j / ( d_j * ( j - 1 ) )
        a     = j * epsilon + abs ( f0 ) + abs ( x * d1 )
        e2    = e * e * ( a ** ( 2 - 2./j) ) * ( abs ( dJ ) ** (2.0/j) ) 
        
        return VE ( d1 , 4 * e2 ) 

    # =========================================================================
    ## calculate  1st (and optionally Jth) derivative with the given step
    #  - f'      is calcualted as O(h^(J-1))
    #  - f^{(J)} is calcualted as O(h^2)
    #  Simple adaptive numerical differentiation
    #  @see R. De Levie, "An improved numerical approximation for the first derivative"
    #  @see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
    def derivatives ( self , both , func , x , h , args = () , kwargs = {}) :
        """Calculate  1st (and optionally Jth) derivative with the given step
        - f'      is calcualted as O(h^(J-1))
        - f^{(J)} is calcualted as O(h^2)
        Simple adaptive numerical differentiation
        - see R. De Levie, ``An improved numerical approximation for the first derivative''
        - see https://link.springer.com/article/10.1007/s12039-009-0111-y
        - see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
        """

        o = self.__I
        
        ## calculate differences 
        imax = o + 2 if both else o + 1
        i = 0
        while i < imax : 
            j = i + 1
            self.__df[i] = func ( x + j * h , *args , **kwargs ) - func ( x - j * h , *args , **kwargs )
            i += 1
            
        ## 1) calculate 1st derivative
        d1     = self.cd_D1 [ 0 ] 
        s1     = self.cd_D1 [ 1 ] 
        result = the_dot ( o + 1 , self.__df , d1 ) / ( s1 * h )
        
        if not both : return result 
        
        ## 2) calculate Nth derivative
        d2     = self.cd_DN [ 0 ]
        s2     = self.cd_DN [ 1 ]        
        dN     = the_dot ( o + 2 , self.__df , d2) / ( s2 * h**(self.__J) ) 
        
        return result, dN

    # =====================================================================================
    ## Get the optimal step (and f(x)) for numerical differentiation
    #  \f$ \delta_{opt} = d_j \left( \frac{ \left| f \right| +
    #   \left| x  f^{\prime}\right| } {\left| f^{J} \right| } \right)^{1/J} \f$ 
    #  @code
    #  deLevie = ...
    #  hopt , f0 = deLevier.optimal_step ( fun , x = 0 , h = 0.1 ) 
    #  @endcode 
    def optimal_step ( self , fun , x , h , hmax = -1  , args = () , kwargs = {}) :
        """Get the optimal step (and f(x)) for numerical differentiation
        
        >>> deLevie = ...
        >>> hopt , f0 = deLevier.optimal_step ( fun , x = 0 , h = 0.1 )
        
        """
        
        ## adjust inital step size 
        h = self.adjust_step ( x , h , hmax )

        
        ## get the function value at the given point 
        f0    = fun ( x , *args , **kwargs )
        
        dx    = delta ( x ) 
        
        xph   = x   + h
        h     = xph - x 

        h_max = hmax if dx < hmax else self.max_step
        
        ## adjust h
        if 0 < h_max < abs ( h ) :
            h = math.copysign ( h_max , h ) 

        i     = self.__I
        j     = self.__J
        
        ## if the intial step is too small, choose another one
        eps_j = self.table2 [ 3 ]
        
        if abs ( h ) <  eps_j or abs ( h ) < dx  :  
            if iszero ( x ) : h = self.deltas [ i ]
            else            : h =  ( 1 + abs ( x ) ) * eps_j ## Sect. 7 
            
        ## 1) find the estimate for first and "J"th the derivative with the given step 
        d1 , dJ = self.derivatives ( True , fun , x , h , args = args , kwargs = kwargs )
            
        ## find the optimal step 
        if iszero  ( dJ ) or ( iszero ( f0 ) and iszero ( x * d1 ) ) :
            
            if iszero ( x ) : hopt = self.deltas [ i ]
            else            : hopt = ( 1 + abs ( x ) ) * eps_j ## Sect. 7
            
        else :
            
            d_j  = self.table2 [ 2 ]
            hopt = d_j * ( ( abs ( f0 ) + abs ( x * d1 ) ) / abs ( dJ ) )** ( 1.0 / j )

        ## final adjustment 
        hopt = self.adjust_step ( x , hopt , hmax )
            
        return hopt , f0  

    # =====================================================================================
    ## Get the value of the 1st derivative using the adaptive rule with the optimal step 
    def derivative ( self , fun , x , h = 0 , args = () , kwargs = {} ) :
        """Get the value of the 1st derivative using the adaptive rule with the optimal step
        """
        
        ## get the optimal step and value f(x) 
        hopt , f0 = self.optimal_step ( fun , x , h , args = args , kwargs = kwargs )

        ## adjust it if needed 
        hopt = self.adjust_step ( x , hopt )
            
        if not self.with_error :
            return self.derivatives ( False , fun , x , hopt , args = args , kwargs = kwargs ) 
        
        ## estimate the uncertainty, if needed  
        d1 , dJ =  self.derivatives ( True  , fun , x , hopt , args = args , kwargs = kwargs  )
        
        i     = self.__I 
        j     = self.__J

        c_j   = self.table2 [ 1 ]
        d_j   = self.table2 [ 2 ]
        
        e     = j * c_j / ( d_j * ( j - 1 ) )
        a     = j * epsilon + abs ( f0 ) + abs ( x * d1 )
        e2    = e * e * ( a ** ( 2 - 2./j) ) * ( abs ( dJ ) ** (2.0/j) ) 
        
        return VE ( d1 , 4 * e2 ) 
                     
# =============================================================================
# The actual setup for adaptive numerical differentiation 
# =============================================================================
_deLevie_  = tuple ( DeLevie ( i                     ) for i in range ( DeLevie.IMAX + 1 ) )
_deLevieE_ = tuple ( DeLevie ( i , with_error = True ) for i in range ( DeLevie.IMAX + 1 ) )



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
#  @param err  (INPUT) calculate the uncertainty?
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
def derivative ( fun , x , h = 0  , I = 2 , err = False , args = ()  , kwargs = {} ) : 
    """Calculate the first derivative for the function
    #  @code
    #  >>> fun =  lambda x : x*x
    #  >>> print derivative ( fun , x = 1 ) 
    - see R. De Levie, ``An improved numerical approximation for the first derivative''
    - see https://link.springer.com/article/10.1007/s12039-009-0111-y
    - see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
    """

    assert isinstance ( I , int ) and 0 <= I , \
           "derivative: Invalid value for ``I''-parameteter!"
    
    if DeLevie.IMAX < I :
        ## logger.warning ("derivative: I/%d exceeds IMAX/%d" % ( I , DeLevie.IMAX ) )
        I = DeLevie.IMAX

    dl = _deLevieE_ if err else _deLevie_ 

    return dl [ I ].derivative ( fun , x , h , args = args  , kwargs = kwargs )

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
    @classprop 
    def IMAX ( cls ) :
        """``IMAX'' : maximal rule index"""
        return DeLevie.IMAX   
    # =========================================================================
    ## constructor 
    #  @param func   the function
    #  @param step   proposed initial step for evaluation of derivatives
    #  @param calc   actual derivateoive calcualtor 
    #  @param err    evaluate numerical uncertainties?
    def __init__ ( self               ,
                   func               ,
                   step       = 0     ,
                   calc       = 2     ,
                   with_error = False ,
                   max_step   = -1    , **kwargs ) :
        """Calculate the first derivative for the function
        R. De Levie, ``An improved numerical approximation for the first derivative''
        see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
        - func:   the function
        - step:   proposed initial step for evaluation of derivatives
        - config  derivative is calculated using 2*I+3 point
        - err    evaluate numerical uncertainties?
        >>> func = math.sin
        >>> deri = Derivative ( func )        
        """
        
        self.__func    = func
        self.__step    = float   ( step  )         
        self.__order   = int     ( calc                               )
        self.__deLevie = DeLevie ( calc ,
                                   max_step   = max_step   , 
                                   with_error = with_error , **kwargs )
        
        
    # =========================================================================
    ## evaluate the derivative
    #  @code 
    #  func  = math.sin
    #  deriv = Derivative ( func )
    #  print deriv(0.1)
    #  @endcode 
    def __call__ ( self , x , *args  , **kwargs  ) :
        """Calculate the first derivative for the function
        R. De Levie, ``An improved numerical approximation for the first derivative''
        see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
        
        >>> func  = math.sin
        >>> deriv = Derivative ( func )
        
        >>> print deriv(0.1) 
        """

        return self.__deLevie.derivative ( self.__func     ,
                                           x               ,
                                           self.__step     ,  
                                           args   = args   ,
                                           kwargs = kwargs )
    
    @property
    def func       ( self ) :
        """The function to be differentiated"""
        return self.__func
    @property
    def step       ( self ) :
        """proposed initial step for evaluation of derivatives"""
        return self.__step
    @property
    def order      ( self ) :
        """ Derivative is calculated using 2*order+1 point"""
        return self.__order
    @property
    def with_error (  self ) :
        """``with_error'' : evaluate numerical uncertainties?"""
        return self.__deLevie.with_error 
    @property
    def calculator ( self ) :
        """``calculator'' : get actual calcualtor for numeric derivatives"""
        return self.__deLevie


    
        
# =============================================================================
## Calculate the partial derivative for the function
#  @see R. De Levie, "An improved numerical approximation for the first derivative"
#  @see https://link.springer.com/article/10.1007/s12039-009-0111-y
#  @see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
#  @code
#  >>> fun2 =  lambda x,y : x*x+y*y
#  >>> print ( partial_derivative ( 0 , fun , (1.0,2.0) ) ) )
#  @endcode
#  @param index  (INPUT) indx of the variable 
#  @param func   (INPUT) the function itself
#  @param x      (INPUT) the argument
#  @param h      (INPUT) the guess for the step used in numeric differentiation
#  @param I      (INPUT) the rule to be used ("N-point rule" = 2*I+1)
#  @param err    (INPUT) calcualte the uncertainty?
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
def partial_derivative ( index , func , x , h = 0  , I = 2 , err = False ) : 
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
    >>> print ( partial_derivative ( 0 , fun , (1.0,2.0) ) ) )
    - see derivative
    """
    
    if len ( x ) <= index :
        raise AttributeError("Invalid argument length/index %d/%d" %  ( len(x) , index ) )
    
    _x =  [ float ( a ) for a in x ]
    
    ## create wrapper function 
    def _wrap ( z ) :
        _z           = _x [ index ] 
        _x [ index ] =  z
        _r = func ( *_x )
        _x [ index ] = _z
        return _r
    
    x_i = _x[ index ]
    
    return derivative ( _wrap , x = x_i , h = h , I = I , err = err )

# =============================================================================
## calcuate the partial derivative for the function
#  @code
#  func = lambda x, y: x * x + y * y  
#  dFdX = PartialDerivative ( 0 , func )
#  dFdY = PartialDerivative ( 1 , func )
#  x = 1
#  y = 2
#  print ' f(%f,%f)=%f    ' % ( x , y , func( x, y ) ) 
#  print ' dFdX=%f dFdY=%f' % ( dFdX(x,y), dFdY ( x, y ) ) 
#  @endcode 
#  @see Derivative
#  @see partial 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-01-25
class PartialDerivative(Derivative) :
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
    def __init__ ( self               ,
                   index              ,   ## index of the variable 
                   func               ,   ## the function 
                   step       = 0     ,   ## proposed initial step for derivatives
                   order      = 2     ,   ## J=2*I(+1) is a number of points used for evaluation of derivative
                   with_error = False ,   ## estimate the uncertainty?
                   max_step   = -1    , **kwargs ) :

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
        assert isinstance ( index , int ) and 0 <= index, \
               "Invalid variable index %s" % index 
        
        ## get the index 
        self.__index = index 
        
        ## initialize the base 
        Derivative.__init__ ( self , func , step , order , with_error = with_error , max_step = max_step , **kwargs )
        
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
        return partial_derivative ( self.__index    ,
                                    self.func       ,
                                    x               ,
                                    self.step       ,
                                    self.order      ,
                                    self.with_error )

    @property 
    def index ( self ) :
        "Index of the variable to be differentiated"
        return self.__index 
        
# =============================================================================
## Calculate the  derivative from analytical function \f$ \frac{{\mathrm{d}} f }{{\mathrm{d}} z }\f$.
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
#  @param I    (INPUT) the rule to be used ("N-point rule" = 2*I+3)
#  @param err  (INPUT) calcualte the uncertainty?
#  @return the derivative d(fun)/dz  at point <code>z</code>  
def complex_derivative ( fun , z , h = 0 , I = 2 , err = False ,  real = True , imag = True ) :
    """Get a complex derivative
    - The function is assumed to be analytical (Cauchy-Riemann conditions are valid).
    >>> fun = ...
    >>> d   = complex_derivative ( fun , x = 1+2j ) 
    """
    
    Z = complex ( z )
    
    X = Z.real
    Y = Z.imag

    ## few alternatives to calculate the real and imaginary part
    
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
class ComplexDerivative(object) :
    """Calculate the first derivative for the function
    - see R. De Levie, ``An improved numerical approximation for the first derivative''
    - see https://link.springer.com/article/10.1007/s12039-009-0111-y
    - see http://www.ias.ac.in/chemsci/Pdf-Sep2009/935.pdf
    >>> func = cmath.sin
    >>> deri = ComplexDerivative ( func )        
    """
    def __init__ ( self            ,
                   func            ,
                   real     = True ,   ## use real to calculate real 
                   imag     = True ,   ## use imag to calculate imag 
                   step     = 0    ,
                   calc     = 1    ,   ## 5-point rule
                   max_step = -1   ) : ## maximal allowed stap 
        
        self.__func    = func
        self.__step    = float   ( step  )

        if isinstance  ( calc , int     ) and 0 <= calc :
            
            self.__order   = int     ( calc )
            self.__deLevie = DeLevie ( calc , with_error = False , step = step , max_step = max_step )
            
        elif isinstance ( calc , Derivative ) : 

            self.__order    = calc.order
            self.__deLevie  = calc.calculator 
            
        elif isinstance ( calc, DeLevie ) :
            
            self.__order   = calc.I
            self.__deLevie = calc 

        else :
            
            raise TypeError ( "Invalid type of ``calc'': %s/%s" % ( calc , type ( calc ) ) )
        
        
        self.__real  = True if real else False 
        self.__imag  = True if imag else False

    ## calculate the complex derivative 
    def __call__ ( self , x , args = () , kwargs = {} ) : 
        """Calculate the complex derivative
        """
        
        Z = complex ( z )
        
        X = Z.real
        Y = Z.imag
        
        if  self.__real :
            
            UX =  lambda x,*p,**kw : complex ( fun ( complex ( x , Y , *p , **kw ) ) ).real
            
            ## calculate the real part 
            re = self.__deLevie.derivative ( UX , X , self.__step , args = args , kwargs = kwargs )
            
        else :
            
            VY =  lambda y,*p,**kw : complex ( fun ( complex ( X , y ) , *p , **kw ) ).imag
            
            ## calculate the real part 
            re =  self.__deLevie.derivative ( VY , Y , self.__step , args = args , kwargs = kwargs )
            
        if self.__imag :
            
            VX =  lambda x,*p,**kw : complex ( fun ( complex ( x , Y ) , *p , **kw ) ).imag
            
            ## get imaginary part 
            im =  self.__deLevie.derivative ( VX , X , self.__step , args = args , kwargs = kwargs )
            
        else :
            
            UY =  lambda y,*p,**kw : complex ( fun ( complex ( X , y ) , *p , **kw ) ).real
            
            ## get imaginary part: nute minus sign!  
            im = - self.__deLevie.derivative ( UY , Y , self.__step , args = args , kwargs = kwargs )
 
        return complex ( re , im )
    
    @property
    def func       ( self ) :
        """The function to be differentiated"""
        return self.__func
    @property
    def step       ( self ) :
        """proposed initial step for evaluation of derivatives"""
        return self.__step
    @property
    def order      ( self ) :
        """ Derivative is calculated using 2*order+3 points"""
        return self.__order
    @property
    def calculator ( self ) :
        """``calculator'' : get actual calculator for numeric derivatives"""
        return self.__deLevie
    

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if not 0.75 < epsilon * 2**52 < 1.25 :
        logger.warning ('"epsilon" in not in the expected range! Math could be suboptimal')

# =============================================================================
##                                                                      The END 
# =============================================================================
