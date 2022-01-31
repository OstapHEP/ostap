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
    'EvalVE'             , ## evaluate the function taking argument's uncertainty 
    'Eval2VE'            , ## evaluate 2-argument function with argument's uncertainties
    'EvalNVE'            , ## evaluate N-argument function with argument's uncertainties
    'EvalNVEcov'         , ## evaluate N-argument function with argument's uncertainties
    'EvalNVEcor'         , ## evaluate N-argument function with argument's uncertainties
    ##
    'Derivative1'        , ## calculate 1st derivative
    'Derivative2'        , ## calculate 2nd derivative
    'Derivative3'        , ## calculate 3rd derivative
    'Derivative4'        , ## calculate 4th derivative
    'Derivative5'        , ## calculate 5th derivative
    'Derivative6'        , ## calculate 6th derivative
    ) 
# =============================================================================
from   builtins import range 
import ROOT, math, abc, array  
# =============================================================================
from ostap.math.base        import Ostap, iszero , isequal
from ostap.core.ostap_types import num_types , is_integer
from ostap.math.ve          import VE 
from sys                    import float_info
from ostap.math.finitediffs import  ( Rule , the_dot , darray ,
                                      Derivative1 , Derivative2 , Derivative3 ,
                                      Derivative4 , Derivative5 , Derivative6 )
                                      
from ostap.math.delevie     import ( derivative         ,        Derivative ,
                                     complex_derivative , ComplexDerivative , 
                                     partial_derivative , PartialDerivative ) 

# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.derivative' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
_eps_ = float_info.epsilon
if not 0.75 < _eps_ * 2**52 < 1.25 :
    import warnings
    warnings.warn ('"epsilon" in not in the expected range! Math could be suboptimal')
    
_next_double_ = Ostap.Math.next_double
_mULPs_       = 1000 
def _delta_ ( x , ulps = _mULPs_ ) :
    n1 = _next_double_ ( x ,  ulps )
    n2 = _next_double_ ( x , -ulps )
    return max ( abs ( n1 - x ) , abs ( n2 - x ) )

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

    ## get a value of function 
    def func_eval        ( self , x , args = () , kwargs = {} ) :
        """Evaluate a function"""
        return self.__func( float ( x ) , *args , **kwargs )
    ## get a value of derivative
    def derivative_eval  ( self , x , args = () , kwargs = {} ) :
        """Evalaute the derivative"""
        return self.__derivative ( float ( x ) , *args , **kwargs )
    
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
    def __call__ ( self , x , args = () , kwargs = {} ) :
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
        val  = self.func_eval ( x , args = args , kwargs = kwargs )
        #
        ## no uncertainties? 
        if   isinstance ( x , num_types          ) : return VE ( val , 0 )
        # ignore small or invalid uncertanties 
        elif 0 >= x.cov2() or iszero ( x.cov2()  ) : return VE ( val , 0 )
        #
        ## 2) evaluate the derivative
        d    = self.derivative_eval (  float ( x ) , args = args , kwargs = kwargs ) 
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
    >>> fun2e  = EvalNVE( 2 , func )
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
            
            if i < len ( partial ) and callable ( partial [ i ] ) : d_i = partial [i] 
            else                                                  : d_i = PartialDerivative ( i , func )
                
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

    # =========================================================================
    ## calculate the gradient (as array)
    #  @code
    #  fun = ...
    #  g   = fun.gradient ( x , y , x , t ) 
    #  @endcode
    def gradient ( self , *x ) :
        """Calculate the gradient (as array)
        >>> fun = ...
        >>> g   = fun.gradient ( x , y , x , t ) 
        """
        n = self.__N
        assert len ( x ) == n , 'Invalid argument size'

        xx = tuple ( VE(i).value() for i in x )

        return darray ( self.partial[i](*xx) for i in range ( n ) ) 
        
    # =========================================================================
    ## the main method (assume that all x are independent) 
    def __call__ ( self , *x ) :
        """The main method (assume that all x are independent)
        """

        n = self.__N 
        assert len ( x ) == n , 'Invalid  argument size'

        xve = tuple ( VE(i)        for i in  x   ) ## force everything to be VE 
        xv  = tuple (    i.value() for i in  xve ) ## get only the values
        
        ## value of the function 
        value = self.__func ( *xv )
        
        ## calculate the covariance 
        cov2 = 0.0
        for i , v in enumerate ( xv ) : 

            ## get covariance 
            c2 = v.cov2()
            if c2 <= 0 or iszero ( c2 ) : continue

            ## calculate the partical derivative 
            df = self.partial [ i ] ( *xv ) 
            if iszero ( df )            : continue 

            ## update covariance for result 
            cov2 += c2 * df * df 

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
        yc2 = y.cov2()
        
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
            from ostap.math.linalg import  mgetter as  get_cov2
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
            from ostap.math.linalg import mgetter as get_corr
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
    def __init__ ( self          ,
                   func          ,
                   real   = True , ## use real to calculate real 
                   imag   = True , ## use imag to calculate imag 
                   step   = 0    ,
                   calc   = 1    ) : ## 5-point rule 
        
        self.__func    = func
        self.__step    = float   ( step  )

        if isinstance  ( calc , int     ) and 0 <= calc :
            
            self.__order   = int     ( calc )
            self.__deLevie = DeLevie ( calc )
            
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
## some rules for (non-adaptive) numerical differentiation 
# =============================================================================


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    
    ## ## the function 
    ## func    = math.sin
    ## ## analytical  derivative 
    ## deriv_a = math.cos

    ## funcs       = [ math.sin for i in range ( DeLevie.IMAX ) ] 
    ## derivatives = ( lambda x :    math.cos(x) ,
    ##                 lambda x : -1*math.sin(x) ,
    ##                 lambda x : -1*math.cos(x) ,
    ##                 lambda x :    math.sin(x) ) 

    
    ## table =  [ ['Function'] + [ 'I=%d' % i for i in range ( DeLevie.IMAX  ) ] ] 
    ## for d , D in enumerate  ( progress_bar ( derivatives ) ) :
        
    ##     funcs = [ Derivative ( funcs [i] , order = i ) for i in range ( DeLevie.IMAX ) ] 
    ##     row = [ '%2d' % (d+1) ]
        
    ##     for f in funcs : 
            
    ##         cnt1 = SE ()            
    ##         for j in range ( 100 ) :
                
    ##             x    = random.uniform ( 0 , 2 * math.pi )
    ##             res  = f     ( x   ) 
    ##             dif  = float ( res ) - D ( x ) 
    ##             cnt1 += dif 
                
    ##         mmax1 = abs ( cnt1.mean ()  *10**12 ) 
    ##         row.append ( '%7.3f ' % mmax1 )
                
    ##     table.append ( row )

    ## import ostap.logger.table as T
    ## table = T.table ( table , prefix = '# ' , alignment=9*'c' )
    ## logger.info ('Sequential differentiation: Mean ifference [10^12]\n%s' % table ) 
    

    ## for i in range ( 0 , 20 ) : 

    ##     x = random.uniform( 0 , 0.5*math.pi )
        
    ##     fmt = "x=%10.5g f=%10.5g delta(f')= %+10.4g delta(f\")= %+10.4g %s %s "
    ##     logger.info (  fmt %  ( x       ,
    ##                             func(x) ,
    ##                             1-deriv_1(x)/deriv_a(x) ,
    ##                             1+deriv_2(x)/func   (x) ,
    ##                             1+deriv_3(x)/deriv_a(x) ,
    ##                             1-deriv_4(x)/func   (x) ) ) 
        


# =============================================================================
##                                                                      The END 
# =============================================================================
