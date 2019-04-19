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
#  @attention It could be coded in much simpler way,
#             but since it is extensively used for pickling,
#             we need to  avoid lambdas and non-pickable functions  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2019-02-18
# =============================================================================
"""Collection of simple classes to wrap the basic operations for callables 

>>> fun1 = ...
>>> fun2 = ...
>>> ops  = Operation2( fun1 ,  fun2 , lambda a , b : a + b )

>>> sin2 = Mul ( math.sin , math.sin )
>>> cos2 = Mul ( math.cos , math.cos )
>>> fun  = Sum ( sin2     , cos2  )

>>> p2  = Pow ( lambda x : x , lambda x : x**2  )

NB: It could be coded in much simpler way,
but since it is extensively used for pickling,
we need to  avoid lambdas and non-pickable functions  

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
    ) 
# =============================================================================
import operator
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
## @class Constant
#  trivial "constant"  function
class Constant(object) :
    """Trivial ``constant'' function
    """
    def __init__ ( self , value ) :
        self.__value = value
    def __call__ ( self , *x    ) :
        return  self.__value
    @property
    def value ( self ) :
        """``value'' : value of the constant"""
        return self.__value
# =============================================================================
## helper class to wrap certain operation
#  @code
#  fun1  = math.sin
#  fun2  = lambda x : x**2
#  fsum  = WrapOper ( fun1 , fun2 , operator.add ) 
#  @endcode 
class WrapOper2(object) :
    """Helper class to wrap certain operation
    >>> fun1  = math.sin
    >>> fun2  = lambda x : x**2
    >>> fsum  = WrapOpen ( fun1 , fun2 , operator.add ) 
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
        
        return oper (  afun ( *x )  , bfun ( *x ) )
        
    @property
    def a ( self ) :
        """``a'' - the first operand/function"""
        return self.__afun
    
    @property
    def b ( self ) :
        """``b'' - the second operand/function"""
        return self.__bfun 

    @property
    def operation ( self ):
        """``operation'' - the actual operation"""
        return  self.__oper
    
# =============================================================================
## @class Descartes
#  Make the Descartes product of two functions
#  @code
#  funx = math.sin
#  funy = math.sin
#  fun2 = Descartes (  funx , funy , 1 )
#  @endcode
class Descartes(object) :
    """Make the Descartes product of two functions
    >>> funx = math.sin
    >>> funy = math.sin
    >>> fun2 = Descartes (  funx , funy , 1 )
    """
    def __init__  ( self , a , b , N = 1 ) :
        """Constructor from  two  functions and arity of the first function
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
            afun = Constant ( float ( a ) ) 
            
        ## trivial case 
        if isinstance ( b , num_types ) :
            bfun = Constant ( float ( b ) ) 
            
        assert callable ( afun ) , 'Invalid type of the first  operand: %s/%s' %  ( a , type( a ) )  
        assert callable ( bfun ) , 'Invalid type of the second operand: %s/%s' %  ( a , type( b ) )  

        self.__afun = afun 
        self.__bfun = bfun 
        self.__N = N
        

    def __call__ ( self , *x ) :

        afun = self.__afun
        bfun = self.__bfun
        n    = self.__N

        return afun ( *x[:n] ) * bfun ( *x[n:] )
        
    @property
    def a ( self ) :
        """``a'' - the first operand/function"""
        return self.__afun
    
    @property
    def b ( self ) :
        """``b'' - the second operand/function"""
        return self.__bfun 

    @property
    def N ( self ):
        """``N'' - arity  of the first function"""
        return  self.__N


# =============================================================================
## @class Compose
#  trivial  composition of two functions
#  @code
#  fun1 = math.sin
#  fun2 = lambda x : x*x
#  comp = Compose (fun2 , fun1 )
#  @endcode 
class Compose(object) :
    """trivial  composition of two functions
    >>> fun1 = math.sin
    >>> fun2 = lambda x : x*x
    >>> comp = Compose ( fun2 , fun1 )
    """
    def __init__ ( self , a , b )  :
        
        afun = a
        bfun = b
        
        ## trivial case 
        if isinstance ( a , num_types ) :
            afun = Constant ( float ( a ) ) 
            
        ## trivial case 
        if isinstance ( b , num_types ) :
            bfun = Constant ( float ( b ) ) 
            
        assert callable ( afun ) , 'Invalid type of the first  operand: %s/%s' %  ( a , type( a ) )  
        assert callable ( bfun ) , 'Invalid type of the second operand: %s/%s' %  ( a , type( b ) )  
        self.__afun = afun 
        self.__bfun = bfun 
        
    def __call__ ( self , *x ) :

        afun = self.__afun
        bfun = self.__bfun
        
        br  = bfun ( *x )
        if not isinstance ( br , list_types ) : br = br ,
        
        return afun ( *br )
        
    @property
    def a ( self ) :
        """``a'' - the first operand/function"""
        return self.__afun
    
    @property
    def b ( self ) :
        """``b'' - the second operand/function"""
        return self.__bfun 

        
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
class Operation2(object) :
    """Helper base class to define the operation
    >>> fun1 = ...
    >>> fun2 = ...
    >>> ops  = Operation2( fun1 ,  fun2 , lambda a , b : a + b )      
    """
    def __init__ ( self , a ,  b , oper ) : 

        afun = a
        bfun = b

        ## trivial case 
        if isinstance ( a , num_types ) :
            afun = Constant ( float ( a ) ) 
            
        ## trivial case 
        if isinstance ( b , num_types ) :
            bfun = Constant ( float ( b ) ) 
            
        assert callable ( afun ) , 'Invalid type of the first  operand: %s/%s' %  ( a , type( a ) )  
        assert callable ( bfun ) , 'Invalid type of the second operand: %s/%s' %  ( a , type( b ) )  

        self.__afun = afun
        self.__bfun = bfun
        self.__oper = oper
        
        ## try to use the native operations, if defined  

        native = [  ( afun , bfun ) ]

        if   isinstance ( afun , Operation2 ) and isinstance ( bfun , Operation2 ) :
            native.append ( ( afun.result , bfun.result ) )
        elif isinstance ( afun , Operation2 ) :        
            native.append ( ( afun.result , bfun        ) )
        elif isinstance ( bfun , Operation2 ) :        
            native.append ( ( afun        , bfun.result ) )

        funab = None 
        for _a , _b in native :
            if funab is None : 
                try :
                    funab = self.__oper ( _a , _b )
                    if not callable ( funab ) :
                        funab = None 
                        raise TypeError ()
                    self.__shortcut = True
                    ## logger.info ('Native  operation is used %s vs %s' % ( type(a) , type(b) ) ) 
                except :
                    pass                
        del native
                    
        ## use generic wrapper, if no native operations are defined 
        if funab is None :
            funab = WrapOper2 ( afun , bfun , oper )
            self.__shortcut = False 
            ## logger.info ('Generic operation is used %s vs %s' % ( type(a) , type(b) ) ) 

        ## store the result 
        self.__funab = funab 

    def __call__ ( self , *x ) :

        fab = self.__funab 
        
        return fab ( *x ) 
    
    @property
    def a ( self ) :
        """``a'' - the first operand/function"""
        return self.__afun
    
    @property
    def b ( self ) :
        """``b'' - the second operand/function"""
        return self.__bfun 

    @property
    def operation ( self ):
        """``operation'' - the actual operation"""
        return  self.__oper
    
    @property
    def result (  self ) :
        """``result'': the result of operation"""
        return self.__funab
    
    @property
    def shortcut ( self ) :
        """``shortcut'' :  Has native/shortcut approach used?"""
        return  self.__shortcut
    
# =============================================================================
## @class Sum
#  Summation operation
#  @code
#  op = Sum ( math.sin , math.cos ) 
#  @endcode 
class Sum (Operation2)  :
    """Summation operation
    >>> op = Sum ( math.sin , math.cos )
    """
    def __init__ ( self , a , b )  :
        super(Sum,self).__init__ ( a , b , operator.add )

# =============================================================================
## @class Sub
#  Subtraction operation
#  @code
#  op = Sub ( math.sin , math.cos ) 
#  @endcode 
class Sub (Operation2)  :
    """Subtraction operation
    >>> op = Sub ( math.sin , math.cos )
    """
    def __init__ ( self , a , b )  :
        super(Sub,self).__init__ ( a , b , operator.sub )

# =============================================================================
## @class Mul
#  Multiplication operation
#  @code
#  op = Mul ( math.sin , math.cos ) 
#  @endcode 
class Mul (Operation2)  :
    """Multiplication operation
    >>> op = Mul ( math.sin , math.cos )
    """
    def __init__ ( self , a , b )  :
        super(Mul,self).__init__ ( a , b , operator.mul )

# =============================================================================
## @class Div
#  Division operation
#  @code
#  op = Div ( math.sin , math.cos ) 
#  @endcode 
class Div (Operation2)  :
    """Division operation
    >>> op = Div ( math.sin , math.cos )
    """
    def __init__ ( self , a , b )  :
        super(Div,self).__init__ ( a , b , operator.truediv )

# =============================================================================
## @class Pow
#  Pow-operation
#  @code
#  op = Pow ( lambda x : x , lambda x : x**2 ) 
#  @endcode 
class Pow (Operation2)  :
    """Pow-operation
    >>> op = Pow ( lambda x : x , lambda x : x**2  )
    """
    def __init__ ( self , a , b )  :
        super(Pow,self).__init__ ( a , b , operator.pow )

# =============================================================================
## @class Mod
#  modulo operation
#  @code
#  op = Mod ( math.sin , math.cos ) 
#  @endcode 
class Mod (Operation2)  :
    """Modulo operation
    >>> op = Mod ( math.sin , math.cos )
    """
    def __init__ ( self , a , b )  :
        super(Mod ,self).__init__ ( a , b , operator.mod )

# =============================================================================
## @class Max
#  Max-operation
#  @code
#  op = Max ( math.sin , math.cos ) 
#  @endcode 
class Max (Operation2)  :
    """Max-operation
    >>> op = Max ( math.sin , math.cos ) 
    """
    def __init__ ( self , a , b )  :        
        super(Max,self).__init__ ( a , b , _Fmax() )

# =============================================================================
## @class Min
#  Min-operation
#  @code
#  op = Min ( math.sin  , math.cos ) 
#  @endcode 
class Min (Operation2)  :
    """Min-operation
    >>> op = Min ( math.sin , math.cos ) 
    """
    def __init__ ( self , a , b )  :
        super(Min,self).__init__ ( a , b , _Fmin() )

# =============================================================================
## @class Or_l
#  OR (logical) - operation
#  @code
#  op = LOr ( fun1 , fun2 ) 
#  @endcode 
class Or_l (Operation2)  :
    """OR (logical) -operation
    >>> op = Or_l ( fun1 , fun2 ) 
    """
    def __init__ ( self , a , b )  :        
        super(Or_l,self).__init__ ( a , b , _For () )

# =============================================================================
## @class And_l
#  AND (logical) -operation
#  @code
#  op = And_l ( fun1 , fun2 ) 
#  @endcode 
class And_l (Operation2)  :
    """AND (logical) -operation
    >>> op = And_l ( fun1 , fun2 ) 
    """
    def __init__ ( self , a , b )  :
        super(And_l,self).__init__ ( a , b , _Fand() )


# =============================================================================
## @class Or_b
#  OR (bitwise) - operation
#  @code
#  op = Or_b ( fun1 , fun2 ) 
#  @endcode 
class Or_b (Operation2)  :
    """OR (bitwise) -operation
    >>> op = Or_b ( fun1 , fun2 ) 
    """
    def __init__ ( self , a , b )  :        
        super(Or_b,self).__init__ ( a , b , operator.or_ )

# =============================================================================
## @class And_b
#  AND (bitwise) -operation
#  @code
#  op = And_b ( fun1 , fun2 ) 
#  @endcode 
class And_b (Operation2)  :
    """AND (bitwise) -operation
    >>> op = And_b ( fun1 , fun2 ) 
    """
    def __init__ ( self , a , b )  :
        super(And_b,self).__init__ ( a , b , operator.and_ )

# =============================================================================
## @class Xor_b
#  XOR (bitwise) - operation
#  @code
#  op = Xor_b ( fun1 , fun2 ) 
#  @endcode 
class Xor_b (Operation2)  :
    """XOR (bitwise) -operation
    >>> op = Xor_b ( fun1 , fun2 ) 
    """
    def __init__ ( self , a , b )  :        
        super(Xor_b,self).__init__ ( a , b , operator.xor )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The  END
# =============================================================================

