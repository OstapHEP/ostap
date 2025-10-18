#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/operations.py
#  Collection of simple classes to wrap the basic math operations for callables 
#  @code
#  fun1 = ...
#  fun2 = ...
#  ops  = Operation2( fun1 ,  fun2 , lambda a , b : a + b )
# 
#  sin2 = Mul ( math.sin , math.sin )
#  cos2 = Mul ( math.cos , math.cos )
#  fun  = Sum ( sin2 , cos2  )
#
#  p2  = Pow ( lambda x : x , lambda x : x**2  )
#  @endcode
#
#  Integration, differentiation, moments, etc are available  via
#  - ostap.math.derivative
#  - ostap.math.integral 
#  - ostap.stats.moments  
#
#  + converison to ROOT.TF(1,2,3)
#  + visualization via ROOT.TF(1,2,3)
#
#  @attention It could be coded in much simpler way,
#             but since it is extensively used for pickling,
#             we need to  avoid lambdas and non-pickable functions  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2019-02-18
# =============================================================================
""" Collection of simple classes to wrap the basic operations for callables 

>>> fun1 = ...
>>> fun2 = ...
>>> ops  = Operation2( fun1 ,  fun2 , lambda a , b : a + b )

>>> sin2 = Mul ( math.sin , math.sin )
>>> cos2 = Mul ( math.cos , math.cos )
>>> fun  = Sum ( sin2     , cos2  )

>>> p2   = Pow ( lambda x : x , lambda x : x**2  )

Integration, differentiation, moments, etc are available  via
- ostap.math.derivative
- ostap.math.integral 
- ostap.stats.moments  

+ converison to ROOT.TF(1,2,3)
+ visualization via ROOT.TF(1,2,3)

NB: It could be coded in much simpler way,
but since it is extensively used for pickling,
we need to avoid lambdas and non-pickable constructions   

"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2019-02-18"
__version__ = "Version$Revision$"
# =============================================================================
__all__     = (
    'Operation2' , ## wrap the operation
    'WrapOper2'  , ## wrap the operation
    'Constant'   , ## constant function
    'Wrapper'    , ## wrapper function
    'Compose'    , ## composition of two functions 
    'Descartes'  , ## Descartes' product of two functions
    'Sum'        , ## wrapped summation 
    'Sub'        , ## wrapped subtraction
    'Mul'        , ## wrapped multiplication  
    'Div'        , ## wrapped (true)division
    'Pow'        , ## wrapped pow-function 
    'Mod'        , ## wrapped modulo-operation
    'Max'        , ## wrapped max-operation
    'Min'        , ## wrapped min-operation
    'Or_l'       , ## wrapped OR(logical)-operation
    'And_l'      , ## wrapped AND(logical)-operation
    'Or_b'       , ## wrapped OR(bitwise)-operation
    'And_b'      , ## wrapped AND(bitwise)-operation
    'Xor_b'      , ## wrapped XOR(bitwise)-operation
    'Exp'        , ## exponent
    'Log'        , ## natural logarithm
    'Log10'      , ## decimal logarithm
    'Sin'        , ## sine 
    'Sinh'       , ## hyperbolic sine 
    'ASin'       , ## inverse sine 
    'ASinh'      , ## inverse hyperbolic sine 
    'Cos'        , ## cosine
    'Cosh'       , ## hyperbolic cosine 
    'ACos'       , ## inverse cosine 
    'ACosh'      , ## inverse hyperbolic cosine 
    'Tan'        , ## tangent
    'Tanh'       , ## hyperbolic tangent  
    'ATan'       , ## inverse tangent
    'ATanh'      , ## inverse hyperbolic tangent
    'Sech'       , ## hyperbolic secant: 1/cosh 
    'Erf'        , ## error function
    'Erfc'       , ## complementary error function
    'Gamma'      , ## gamma function
    'LogGamma'   , ## logarithm of gamma fnuction
    'iGamma'     , ## 1/Gamma 
    'Sqrt'       , ## square root
    'Cbrt'       , ## cubic root
    'Square'     , ## square 
    'Cube'       , ## cube
    'digitize'   , ## digitize the fuction to numpy-array 
    ) 
# =============================================================================
import operator, abc, math 
from   ostap.core.ostap_types import num_types, integer_types, list_types  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.operations' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
class _Fmax(object) :
    def __call__ ( self , x , y ) : return max ( x , y )
class _Fmin(object) :
    def __call__ ( self , x , y ) : return min ( x , y )
class _For (object) :
    def __call__ ( self , x , y ) : return x or  y 
class _Fand(object) :
    def __call__ ( self , x , y ) : return x and y 
# =============================================================================
#  @class Function
class Function(object) :
    __metaclass__ = abc.ABCMeta
    @abc.abstractmethod
    def __call__ ( self , *args ) :
        """ Use the function! """
        return NotImplemented
    ## simple operations 
    def __add__      ( self , other ) : return Sum ( self  , other )
    def __sub__      ( self , other ) : return Sub ( self  , other )
    def __mul__      ( self , other ) : return Mul ( self  , other )
    def __div__      ( self , other ) : return Div ( self  , other )
    def __truediv__  ( self , other ) : return Div ( self  , other )
    def __pow__      ( self , other ) : return Pow ( self  , other )
    def __mod__      ( self , other ) : return Mod ( self  , other )
    def __radd__     ( self , other ) : return Sum ( other , self  )
    def __rmul__     ( self , other ) : return Mul ( other , self  )
    def __rsub__     ( self , other ) : return Sub ( other , self  )
    def __rdiv__     ( self , other ) : return Div ( other , self  )
    def __rtruediv__ ( self , other ) : return Div ( other , self  )
    def __rpow__     ( self , other ) : return Pow ( other , self  )
    def __rmod__     ( self , other ) : return Mod ( other , self  )

    def __iadd__     ( self , other ) : return NotImplemented
    def __imul__     ( self , other ) : return NotImplemented
    def __isub__     ( self , other ) : return NotImplemented
    def __idiv__     ( self , other ) : return NotImplemented

# =============================================================================
## @class Constant
#  trivial "constant"  function
class Constant(Function) :
    """ Trivial `constant' function
    """
    def __init__ ( self , value ) :
        self.__value = value
    def __call__ ( self , *x    ) :
        return  self.__value
    def __str__ ( self ) : return "%s" % self.value          
    __repr__ = __str__ 
    @property
    def value ( self ) :
        """`value' : value of the constant"""
        return self.__value

# =============================================================================
## @class Wrapper
#  Wrap another function
class Wrapper(Function) :
    """ Trivial wrapper function
    """
    def __init__ ( self , function , name = '' ) :
        if not isinstance ( function , Wrapper ) : self.__function = function
        else                                     : self.__function = function.function        
        self.__name = name if name else str(function) 
    def __call__ ( self , *x ) :
        return  self.__function ( *x ) 
    def __str__ ( self ) : return "%s" % self.function           
    __repr__ = __str__ 
    @property
    def function ( self ) :
        """`function' : the actual function"""
        return self.__function
    @property
    def name     ( self ) :
        """`name' : function name"""
        return self.__name
    def __str__  ( self ) :
        return self.__name
    __repr__ = __str__ 

# =============================================================================
## helper class to wrap certain operation
#  @code
#  fun1  = math.sin
#  fun2  = lambda x : x**2
#  fsum  = WrapOper ( fun1 , fun2 , operator.add ) 
#  @endcode 
class WrapOper2(Function) :
    """Helper class to wrap certain operation
    >>> fun1  = math.sin
    >>> fun2  = lambda x : x**2
    >>> fsum  = WrapOper2 ( fun1 , fun2 , operator.add ) 
    """

    def __init__ ( self , a , b , operation ) :
        
        o = operation 
        assert callable ( a ) , 'Invalid type of the first  operand: %s/%s' %  ( a , type ( a ) )  
        assert callable ( b ) , 'Invalid type of the second operand: %s/%s' %  ( a , type ( b ) )
        assert callable ( o ) , 'Invalid type of operation         : %s/%s' %  ( o , type ( o ) )  

        self.__afun = a
        self.__bfun = b
        self.__oper = o 

    ## the main method 
    def __call__ ( self , *x ) :
        
        afun =  self.__afun
        bfun =  self.__bfun
        oper =  self.__oper
        
        return oper ( afun ( *x ) , bfun ( *x ) )
        
    @property
    def a ( self ) :
        """`a' - the first operand/function"""
        return self.__afun
    
    @property
    def b ( self ) :
        """`b' - the second operand/function"""
        return self.__bfun 

    @property
    def operation ( self ):
        """`operation' - the actual operation"""
        return  self.__oper
    
# =============================================================================
## @class Descartes
#  Make the Descartes product of two functions
#  @code
#  funx = math.sin
#  funy = math.sin
#  fun2 = Descartes (  funx , funy , 1 )
#  @endcode
class Descartes(Function) :
    """ Make the Descartes product of two functions
    >>> funx = math.sin
    >>> funy = math.sin
    >>> fun2 = Descartes (  funx , funy , 1 )
    """
    def __init__  ( self , a , b , N = 1 ) :
        """ Constructor from  two  functions and arity of the first function
        >>> funx = math.sin
        >>> funy = math.sin
        >>> fun2 = Descartes (  funx , funy , 1 )        
        """
        
        assert isinstance ( N , integer_types ) and 0 <= N ,\
               'Invalid arity of the first operand %s/%s' % ( N , type ( N ) )
        
        afun = a
        bfun = b
        
        ## trivial case 
        if isinstance ( a , num_types ) :
            afun = Constant ( a * 1.0 ) 
            
        ## trivial case 
        if isinstance ( b , num_types ) :
            bfun = Constant ( b * 1.0 ) 
            
        assert callable ( afun ) , 'Invalid type of the first  operand: %s/%s' %  ( a , type( a ) )  
        assert callable ( bfun ) , 'Invalid type of the second operand: %s/%s' %  ( a , type( b ) )  

        self.__afun = afun 
        self.__bfun = bfun 
        self.__N = N
        
    # =========================================================================
    def __call__ ( self , *x ) :

        afun = self.__afun
        bfun = self.__bfun
        n    = self.__N

        return afun ( *x[:n] ) * bfun ( *x[n:] )
        
    @property
    def a ( self ) :
        """`a' - the first operand/function"""
        return self.__afun
    
    @property
    def b ( self ) :
        """`b' - the second operand/function"""
        return self.__bfun 

    @property
    def N ( self ):
        """`N' - arity  of the first function"""
        return  self.__N

# =============================================================================
## @class Compose
#  trivial  composition of two functions
#  @code
#  fun1 = math.sin
#  fun2 = lambda x : x*x
#  comp = Compose (fun2 , fun1 )
#  @endcode 
class Compose(Function) :
    """ Trivial  composition of two functions
    >>> fun1 = math.sin
    >>> fun2 = lambda x : x*x
    >>> comp = Compose ( fun2 , fun1 )
    """
    def __init__ ( self , a , b , na = '' , nb = '' )  :
        
        afun = a
        bfun = b
        
        ## trivial case 
        if isinstance ( a , num_types ) : afun = Constant ( a ) 
            
        ## trivial case 
        if isinstance ( b , num_types ) : bfun = Constant ( b ) 
            
        assert callable ( afun ) , 'Invalid type of the first  operand: %s/%s' %  ( a , typename ( a ) )  
        assert callable ( bfun ) , 'Invalid type of the second operand: %s/%s' %  ( a , typename ( b ) )  
        
        self.__afun = afun 
        self.__bfun = bfun
        
        self.__na   = na if na else str ( afun ) 
        self.__nb   = nb if nb else str ( bfun ) 
        
    def __call__ ( self , *x ) :

        afun = self.__afun
        bfun = self.__bfun
        
        br   = bfun ( *x )
        if not isinstance ( br , list_types ) : br = br ,
        
        return afun ( *br )
        
    @property
    def a ( self ) :
        """`a' - the first operand/function"""
        return self.__afun
    
    @property
    def b ( self ) :
        """`b' - the second operand/function"""
        return self.__bfun

    def __str__ ( self ) :
        return "%s(%s)" % ( self.__na , self.__nb )
    __repr__ = __str__ 
        
# =============================================================================
## @class Operation2
#  helper class to  define some operations
#  @code
#  fun1 = ...
#  fun2 = ...
#  ops  = Operation2( fun1 ,  fun2 , lambda a , b : a + b )  
#  @endcode
#  - first it tries to use the native operations for callables
#  - if it fails, generic WrapOper2 is used
class Operation2(Function) :
    """ Helper base class to define the operation
    >>> fun1 = ...
    >>> fun2 = ...
    >>> ops  = Operation2( fun1 ,  fun2 , lambda a , b : a + b )      
    """
    def __init__ ( self , a ,  b , oper , symbol = '?' , prefix = '' ) : 

        afun = a
        bfun = b

        anum = isinstance ( a , num_types )
        bnum = isinstance ( b , num_types )

        funab = None
        if anum and bnum :
            afun  = Constant ( a * 1.0 ) 
            bfun  = Constant ( b * 1.0 ) 
            funab = Constant ( oper  ( a * 1.0 , b * 1.0 ) )
        elif anum :
            afun  = Constant ( a * 1.0 ) 
        elif bnum :
            bfun  = Constant ( b * 1.0 ) 
            
        assert callable ( afun ) , 'Invalid type of the first  operand: %s/%s' %  ( a , typename ( a ) )  
        assert callable ( bfun ) , 'Invalid type of the second operand: %s/%s' %  ( b , typename ( b ) )  

        self.__afun = afun
        self.__bfun = bfun
        self.__oper = oper

        ## try to use the native operations, if/when defined

        if isinstance   ( afun , Operation2 ) : afun = afun.result
        if isinstance   ( bfun , Operation2 ) : bfun = bfun.result
        
        fa = isinstance ( afun , Function   )
        fb = isinstance ( bfun , Function   )

        if fa and fb and not funab : 
            funab = WrapOper2 ( afun , bfun , oper )

        if not funab :

            native = [  ( afun , bfun ) ]
            
            if anum : native = [ ( a    , bfun ) , ( afun , bfun ) ]
            if bnum : native = [ ( afun , b    ) , ( afun , bfun ) ]

            if   isinstance ( afun , Operation2 ) and isinstance ( bfun , Operation2 ) :
                native.append ( ( afun.result , bfun.result ) )
            elif isinstance ( afun , Operation2 ) :        
                native.append ( ( afun.result , bfun        ) )
            elif isinstance ( bfun , Operation2 ) :        
                native.append ( ( afun        , bfun.result ) )

            for _a , _b in native :
                if funab is None : 
                    try :
                        funab = oper ( _a , _b )
                        if not callable ( funab ) :
                            funab = None 
                            raise TypeError ()
                        self.__shortcut = True
                        ## logger.info ('Native  operation is used %s vs %s' % ( type(a) , type(b) ) )
                    except :
                        pass                
            del native
                    
        ## use generic wrapper, if no native operations are defined 
        if not funab :
            funab = WrapOper2 ( afun , bfun , oper )
            self.__shortcut = False 
            ## logger.info ('Generic operation is used %s vs %s' % ( type(a) , type(b) ) ) 

        ## store the result 
        self.__funab  = funab 
        self.__symbol = symbol
        self.__prefix = prefix
        
    def __call__ ( self , *x ) :
        fab = self.__funab         
        return fab ( *x ) 
    
    @property
    def a ( self ) :
        """`a' - the first operand/function"""
        return self.__afun
    
    @property
    def b ( self ) :
        """`b' - the second operand/function"""
        return self.__bfun 

    @property
    def operation ( self ):
        """`operation' - the actual operation"""
        return  self.__oper
    
    @property
    def result (  self ) :
        """`result': the result of operation"""
        return self.__funab
    
    @property
    def shortcut ( self ) :
        """`shortcut' :  Has native/shortcut approach used?"""
        return  self.__shortcut

    def __str__ ( self ) :
        return "%s((%s)%s(%s))" % ( self.__prefix , self.a ,
                                    self.__symbol , self.b )
            
    __repr__ = __str__
    
# =============================================================================
## @class Sum
#  Summation operation
#  @code
#  op = Sum ( math.sin , math.cos ) 
#  @endcode 
class Sum (Operation2)  :
    """ Summation operation
    >>> op = Sum ( math.sin , math.cos )
    """
    def __init__ ( self , a , b )  :
        super(Sum,self).__init__ ( a , b , operator.add , '+' )

# =============================================================================
## @class Sub
#  Subtraction operation
#  @code
#  op = Sub ( math.sin , math.cos ) 
#  @endcode 
class Sub (Operation2)  :
    """ Subtraction operation
    >>> op = Sub ( math.sin , math.cos )
    """
    def __init__ ( self , a , b )  :
        super(Sub,self).__init__ ( a , b , operator.sub , '-' )

# =============================================================================
## @class Mul
#  Multiplication operation
#  @code
#  op = Mul ( math.sin , math.cos ) 
#  @endcode 
class Mul (Operation2)  :
    """ Multiplication operation
    >>> op = Mul ( math.sin , math.cos )
    """
    def __init__ ( self , a , b )  :
        super(Mul,self).__init__ ( a , b , operator.mul  , '*' )

# =============================================================================
## @class Div
#  Division operation
#  @code
#  op = Div ( math.sin , math.cos ) 
#  @endcode 
class Div (Operation2)  :
    """ Division operation
    >>> op = Div ( math.sin , math.cos )
    """
    def __init__ ( self , a , b )  :
        super(Div,self).__init__ ( a , b , operator.truediv , '/')

# =============================================================================
## @class Pow
#  Pow-operation
#  @code
#  op = Pow ( lambda x : x , lambda x : x**2 ) 
#  @endcode 
class Pow (Operation2)  :
    """ Pow-operation
    >>> op = Pow ( lambda x : x , lambda x : x**2  )
    """
    def __init__ ( self , a , b )  :
        super(Pow,self).__init__ ( a , b , operator.pow , '**')

# =============================================================================
## @class Square
#  square-operation
class Square(Pow)  :
    """ Square-operation
    """
    def __init__ ( self , a )  :
        super(Square,self).__init__ ( a , 2 )

# =============================================================================
## @class Cube
#  cube-operation
class Cube(Pow)  :
    """ Cube-operation
    """
    def __init__ ( self , a )  :
        super(Cube,self).__init__ ( a , 3 )

# =============================================================================
## @class Mod
#  modulo operation
#  @code
#  op = Mod ( math.sin , math.cos ) 
#  @endcode 
class Mod (Operation2)  :
    """ Modulo operation
    >>> op = Mod ( math.sin , math.cos )
    """
    def __init__ ( self , a , b )  :
        super(Mod ,self).__init__ ( a , b , operator.mod , '%')

# =============================================================================
## @class Max
#  Max-operation
#  @code
#  op = Max ( math.sin , math.cos ) 
#  @endcode 
class Max (Operation2)  :
    """ Max-operation
    >>> op = Max ( math.sin , math.cos ) 
    """
    def __init__ ( self , a , b )  :        
        super(Max,self).__init__ ( a , b , _Fmax() , ',' , 'max')

# =============================================================================
## @class Min
#  Min-operation
#  @code
#  op = Min ( math.sin  , math.cos ) 
#  @endcode 
class Min (Operation2)  :
    """ Min-operation
    >>> op = Min ( math.sin , math.cos ) 
    """
    def __init__ ( self , a , b )  :
        super(Min,self).__init__ ( a , b , _Fmin() , '' , 'min')

# =============================================================================
## @class Or_l
#  OR (logical) - operation
#  @code
#  op = LOr ( fun1 , fun2 ) 
#  @endcode 
class Or_l (Operation2)  :
    """ OR (logical) -operation
    >>> op = Or_l ( fun1 , fun2 ) 
    """
    def __init__ ( self , a , b )  :        
        super(Or_l,self).__init__ ( a , b , _For () , ' or ' )

# =============================================================================
## @class And_l
#  AND (logical) -operation
#  @code
#  op = And_l ( fun1 , fun2 ) 
#  @endcode 
class And_l (Operation2)  :
    """ AND (logical) -operation
    >>> op = And_l ( fun1 , fun2 ) 
    """
    def __init__ ( self , a , b )  :
        super(And_l,self).__init__ ( a , b , _Fand() , ' and ' )


# =============================================================================
## @class Or_b
#  OR (bitwise) - operation
#  @code
#  op = Or_b ( fun1 , fun2 ) 
#  @endcode 
class Or_b (Operation2)  :
    """ OR (bitwise) -operation
    >>> op = Or_b ( fun1 , fun2 ) 
    """
    def __init__ ( self , a , b )  :        
        super(Or_b,self).__init__ ( a , b , operator.or_ , '|' )

# =============================================================================
## @class And_b
#  AND (bitwise) -operation
#  @code
#  op = And_b ( fun1 , fun2 ) 
#  @endcode 
class And_b (Operation2)  :
    """ AND (bitwise) -operation
    >>> op = And_b ( fun1 , fun2 ) 
    """
    def __init__ ( self , a , b )  :
        super(And_b,self).__init__ ( a , b , operator.and_ , '&' )

# =============================================================================
## @class Xor_b
#  XOR (bitwise) - operation
#  @code
#  op = Xor_b ( fun1 , fun2 ) 
#  @endcode 
class Xor_b (Operation2)  :
    """ XOR (bitwise) -operation
    >>> op = Xor_b ( fun1 , fun2 ) 
    """
    def __init__ ( self , a , b )  :        
        super(Xor_b,self).__init__ ( a , b , operator.xor , '^' )

# =============================================================================
## @class Abs
#  absolute value  
class Abs(Compose) :
    """ Absolute value"""
    def __init__ ( self , func ) :
        super(Abs,self).__init__ ( abs , func , 'abs')

# =============================================================================
## @class Exp
#  simple 'exponent' 
class Exp(Compose) :
    """ Exponent for the function"""
    def __init__ ( self , func ) :
        super(Exp,self).__init__ ( math.exp , func , 'exp')


# =============================================================================
## @class Log
#  simple 'log' 
class Log(Compose) :
    """ Log for the function"""
    def __init__ ( self , func ) :
        super(Log,self).__init__ ( math.log , func , 'log')

# =============================================================================
## @class Log10
#  simple 'log10' 
class Log10(Compose) :
    """ Log10 for the function"""
    def __init__ ( self , func ) :
        super(Log10,self).__init__ ( math.log10 , func , 'log10' )

# =============================================================================
## @class Sin
#  simple 'sin' 
class Sin(Compose) :
    """ Sin for the function"""
    def __init__ ( self , func ) :
        super(Sin,self).__init__ ( math.sin , func , 'sin')

# =============================================================================
## @class Cos
#  simple 'cos' 
class Cos(Compose) :
    """ Cos for the function"""
    def __init__ ( self , func ) :
        super(Cos,self).__init__ ( math.cos , func , 'cos' )

# =============================================================================
## @class Tan
#  simple 'tan' 
class Tan(Compose) :
    """ Tan for the function"""
    def __init__ ( self , func ) :
        super(Tan,self).__init__ ( math.tan , func , 'tan' )

# =============================================================================
## @class Sinh
#  simple 'sinh' 
class Sinh(Compose) :
    """ Sinh for the function"""
    def __init__ ( self , func ) :
        super(Sinh,self).__init__ ( math.sinh , func , 'sinh' )

# =============================================================================
## @class Cosh
#  simple 'cosh' 
class Cosh(Compose) :
    """ Cos for the function"""
    def __init__ ( self , func ) :
        super(Cosh,self).__init__ ( math.cosh , func , 'cosh' )

# =============================================================================
## @class Tanh
#  simple 'tan' 
class Tanh(Compose) :
    """ Tanh for the function"""
    def __init__ ( self , func ) :
        super(Tanh,self).__init__ ( math.tanh , func , 'tanh' )

# =============================================================================
## @class ASin
#  simple 'asin' 
class ASin(Compose) :
    """ ASin for the function"""
    def __init__ ( self , func ) :
        super(ASin,self).__init__ ( math.asin , func , 'asinh' )

# =============================================================================
## @class ACos
#  simple 'acos' 
class ACos(Compose) :
    """ ACos for the function"""
    def __init__ ( self , func ) :
        super(ACos,self).__init__ ( math.acos , func , 'acos' )

# =============================================================================
## @class ATan
#  simple 'atan' 
class ATan(Compose) :
    """ ATan for the function"""
    def __init__ ( self , func ) :
        super(ATan,self).__init__ ( math.atan , func , 'atan' )

# =============================================================================
## @class ASinh
#  simple 'asinh' 
class ASinh(Compose) :
    """ ASinh for the function"""
    def __init__ ( self , func ) :
        super(ASinh,self).__init__ ( math.asinh , func , 'asinh' )

# =============================================================================
## @class ACosh
#  simple 'acosh' 
class ACosh(Compose) :
    """ ACosh for the function"""
    def __init__ ( self , func ) :
        super(ACosh,self).__init__ ( math.acosh , func , 'acosh' )

# =============================================================================
## @class ATanh
#  simple 'atanh' 
class ATanh(Compose) :
    """ ATanh for the function"""
    def __init__ ( self , func ) :
        super(ATanh,self).__init__ ( math.atanh , func , 'atanh' )

# =============================================================================
## @class Erf
#  simple 'erf' 
class Erf(Compose) :
    """ Erf for the function"""
    def __init__ ( self , func ) :
        super(Erf,self).__init__ ( math.erf , func ,  'erf' )

# =============================================================================
## @class Erfc
#  simple 'erfc' 
class Erfc(Compose) :
    """ Erfc for the function"""
    def __init__ ( self , func ) :
        super(Erfc,self).__init__ ( math.erfc , func , 'erfc' )

# =============================================================================
## @class Gamma
#  simple 'Gamma' 
class Gamma(Compose) :
    """ Gamma for the function"""
    def __init__ ( self , func ) :
        super(Gamma,self).__init__ ( math.gamma , func , 'G' )

# =============================================================================
## @class LogGamma
#  simple 'LogGamma' 
class LogGamma(Compose) :
    """ LogGamma for the function"""
    def __init__ ( self , func ) :
        super(LogGamma,self).__init__ ( math.lgamma , func , 'lnG')

# =============================================================================
## @class iGamma
#  1/Gamma
class iGamma(Compose) :
    """ 1/Gamma for the function"""
    def __init__ ( self , func ) :
        from ostap.math.math_ve import igamma as _igamma         
        super(iGamma,self).__init__ ( _igamma , func , 'iG')

# =============================================================================
## @class Sqrt
#  simple 'Sqrt' 
class Sqrt(Compose) :
    """ Square root for the function"""
    def __init__ ( self , func ) :
        super(Sqrt,self).__init__ ( math.sqrt , func , 'sqrt')

# =============================================================================
## @class Cbrt
#  Cubic root 
class Cbrt(Compose) :
    """ Cubic root for the function"""
    def __init__ ( self , func ) :
        from ostap.math.math_ve import cbrt as _cbrt 
        super(Cbrt,self).__init__ ( _cbrt , func , 'cbrt')

# =============================================================================
## @class Sech
#  1/cosh 
class Sech (Compose) :
    """ Sech (1/cosh) for the function"""
    def __init__ ( self , func ) :
        from ostap.math.math_ve import sech as _sech 
        super(Sech,self).__init__ ( _sech , func , 'sech')
        
# =============================================================================
## convert function object to TF1 
def tf1  ( fun , **kwargs ) :
    """ Convert function object to TF1     
    """
    from ostap.math.models import tf1 as _tf1 
    return _tf1 ( fun , **kwargs )

# =============================================================================
## convert function object to TF2 
def tf2  ( fun , **kwargs ) :
    """ Convert function object to TF2     
    """
    from ostap.math.models import tf2 as _tf2
    return _tf2 ( fun , **kwargs )

# =============================================================================
## convert function object to TF3 
def tf3  ( fun , **kwargs ) :
    """ Convert function object to TF3     
    """
    from ostap.math.models import tf3 as _tf3 
    return _tf3 ( fun , **kwargs )

# =============================================================================
## draw function via conversion to TF1 
def draw1D ( fun , **kwargs ) :
    """ draw function via conversion to TF1"""
    from ostap.math.models import f1_draw
    return f1_draw ( fun , **kwargs )

# =============================================================================
## draw function via conversion to TF2 
def draw3D ( fun , **kwargs ) :
    """ draw function via conversion to TF2"""
    from ostap.math.models import _f2_draw_
    return f2_draw ( fun , **kwargs )

# =============================================================================
## draw function via conversion to TF3 
def draw3D ( fun , **kwargs ) :
    """ draw function via conversion to TF3"""
    from ostap.math.models import _f3_draw_
    return f3_draw ( fun , **kwargs ) 

# =============================================================================
## digitize the function as np-array
#  @code
#  >>>  fun   = ...
#  >>> array1 =    digitize ( fun , 0 , 1 , 1000  )
#  >>> array2 =    digitize ( fun , 0 , 1 , 0.001 )
#  @endcode 
def digitize ( func , xmin , xmax , N ) :
    """ Digitize the function as np-array
    >>>  fun   = ...
    >>> array1 =    digitize ( fun , 0 , 1 , 1000  )
    >>> array2 =    digitize ( fun , 0 , 1 , 0.001 )
    """
    try : 
        import numpy as np
    except ImportError :
        logger.error ( 'NumPy cannot be used!')
        return ()
    
    from   ostap.core.ostap_types import integer_types
    
    if   isinstance   ( N , integer_types ) and 1 < N :
        x = np.linspace ( xmin * 1.0 , xmax * 1.0 , N )
    elif isinstance   ( N , float_types   ) and 0 < N < abs ( xmax - xmin ) : 
        x = np.arange   ( xmin * 1.0 , xmax * 1.0 , N )
    else :
        raise TypeError ('Invalid xmin/xmax/N: %s/%s/%s' % ( xmin ,
                                                             xmax ,
                                                             N    ) )
    ## vectorize function 
    vfunc = np.vectorize ( func)
    
    return  vfunc ( x ) 
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                     The  END
# =============================================================================

