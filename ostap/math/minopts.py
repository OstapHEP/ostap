#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/minopts.py
#  set of utilities to deal with ROOT::Math::MinimizerOptions
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2023-04-07
# =============================================================================
"""Set of utilities to deal with ROOT::Math::MinimizerOptions
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    ## as classes/objects 
    'PrintLevel'         , ## context mamager to deal with DefaultPrintLevel  
    'Strategy'           , ## context manager to deal with DefaultStrategy 
    'Tolerance'          , ## context manager to deal with DefaultTolerance
    'MaxFunctionCalls'   , ## context manager to deal with DefaultMaxFunctionCalls
    'MaxIterations'      , ## context manager to deal with DefaultMaxIterations
    'TypeAlgo'           , ## context manager to deal with DefaultMinimizerType/Algo
    ## as functions
    'print_level'        , ## context mamager to deal with DefaultPrintLevel  
    'strategy'           , ## context manager to deal with DefaultStrategy 
    'tolerance'          , ## context manager to deal with DefaultTolerance
    'max_calls'          , ## context manager to deal with DefaultMaxFunctionCalls
    'max_iterations'     , ## context manager to deal with DefaultMaxIterations
    'type_algo'          , ## context manager to deal with DefaultMinimizerType/Algo
    ) 
# =============================================================================
from   ostap.core.ostap_types import num_types, integer_types, string_types 
import ROOT, ctypes
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger      import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.math.minopts' )
else                       : logger = getLogger( __name__ )
# =============================================================================

RRMMO = ROOT.ROOT.Math.MinimizerOptions
                             
# =============================================================================
## @class PrintLevel
#  Context manager to play with Default Print level
#  @code
#  with PrintLevel ( 0 ) : ... 
#  @endcode 
#  @see ROOT::Math::MinimizerOptions
#  @see ROOT::Math::MinimizerOptions::DefaultPrintLevel
#  @see ROOT::Math::MinimizerOptions::SetDefaultPrintLevel
class PrintLevel(object) :
    """ Context manager to play with Default Print level
    
    >>> with PrintLevel ( 0 ) : ...
    
    - see ROOT::Math::MinimizerOptions
    - see ROOT::Math::MinimizerOptions::DefaultPrintLevel
    - see ROOT::Math::MinimizerOptions::SetDefaultPrintLevel
    """
    def __init__ ( self , level ) :
        assert level in  ( 0 , 1, 2 , 3 ) , 'Invalid print level'
        self.__level = level
    def __enter__ ( self ) :
        self.__prev = RRMMO.DefaultPrintLevel()
        RRMMO.SetDefaultPrintLevel ( self.__level )
        return self
    def __exit__ ( self , *_ ) :
        RRMMO.SetDefaultPrintLevel ( self.__prev  )

# =============================================================================
## @class Strategy
#  Context manager to play with Default Strategy
#  @code
#  with Strategy ( 0 ) : ... 
#  @endcode 
#  @see ROOT::Math::MinimizerOptions
#  @see ROOT::Math::MinimizerOptions::DefaultStrategy
#  @see ROOT::Math::MinimizerOptions::SetDefaultStrategy
class Strategy(object) :
    """Context manager to play with default strategy 

    >>> with Strategy ( 0 ) : ...

    - see ROOT::Math::MinimizerOptions
    - see ROOT::Math::MinimizerOptions::DefaultStrategt 
    - see ROOT::Math::MinimizerOptions::SetDefaultPrintLevel
    """
    def __init__ ( self , strategy ) :
        assert strategy in  ( 0 , 1 , 2 ) , 'Invalid strategy!'
        self.__strategy = strategy
    def __enter__ ( self ) :
        self.__prev = RRMMO.DefaultStrategy()
        RRMMO.SetDefaultStrategy ( self.__strategyl )
        return self
    def __exit__ ( self , *_ ) :
        RRMMO.SetDefaultStrategy ( self.__prev  )

# =============================================================================
## @class Tolerance
#  Context manager to play with Default Tolerance
#  @code
#  with Tolerance ( 0.001 ) : ... 
#  @endcode 
#  @see ROOT::Math::MinimizerOptions
#  @see ROOT::Math::MinimizerOptions::DefaultTolerance
#  @see ROOT::Math::MinimizerOptions::SetDefaultTolerance
class Tolerance(object) :
    """ Context manager to play with default strategy
    
    >>> with Tolerance ( 0.001 ) : ...

    - see ROOT::Math::MinimizerOptions
    - see ROOT::Math::MinimizerOptions::DefaultTolerance
    - see ROOT::Math::MinimizerOptions::SetDefaultTolerance
    """
    def __init__ ( self , tolerance ) :
        assert isinstance ( tolerance , num_types ) and 0 < tolerance , 'Invalid tolerance!'
        self.__tolerance = tolerance
    def __enter__ ( self ) :
        self.__prev = RRMMO.DefaultTolerance()
        RRMMO.SetDefaultTolerance ( self.__tolerance )
        return self
    def __exit__ ( self , *_ ) :
        RRMMO.SetDefaultTolerance ( self.__prev  )

# =============================================================================
## @class MaxFunctionCalls
#  Context manager to play with Default MaxFunctionCalls 
#  @code
#  with MaxFunctionCalls ( 100000 ) : ... 
#  @endcode 
#  @see ROOT::Math::MinimizerOptions
#  @see ROOT::Math::MinimizerOptions::DefaultMaxFunctionCalls
#  @see ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls
class MaxFunctionCalls(object) :
    """Context manager to play with default max functon calls 
    
    >>> with MaxFunctionCalls ( 100000 ) : ...

    - see ROOT::Math::MinimizerOptions
    - see ROOT::Math::MinimizerOptions::DefaultMaxFunctionCalls
    - see ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls
    """
    def __init__ ( self , maxcalls ) :
        assert isinstance ( maxcalls , integer_types ) and 0 < maxcalls , 'Invalid maxcalls!'
        self.__maxcalls = maxcalls
    def __enter__ ( self ) :
        self.__prev = RRMMO.DefaultMaxFunctionCalls ()
        RRMMO.SetDefaultMaxFunctionCalls ( self.__maxcalls )
        return self
    def __exit__ ( self , *_ ) :
        EEMMO.SetDefaultMaxFunctionCalls ( self.__prev  )

# =============================================================================
## @class MaxIterations
#  Context manager to play with Default MaxIterations
#  @code
#  with MaxIterations ( 1000 ) : ... 
#  @endcode 
#  @see ROOT::Math::MinimizerOptions
#  @see ROOT::Math::MinimizerOptions::DefaultMaxIterations
#  @see ROOT::Math::MinimizerOptions::SetDefaultMaxIterations
class MaxIterations(object) :
    """ Context manager to play with default max iterations
    
    >>> with MaxIterations ( 1000 ) : ...

    - see ROOT::Math::MinimizerOptions    - see ROOT::Math::MinimizerOptions
    - see ROOT::Math::MinimizerOptions::DefaultMaxIterations
    - see ROOT::Math::MinimizerOptions::SetDefaultMaxIterations
    """
    def __init__ ( self , maxiter ) :
        assert isinstance ( maxiter , integer_types ) and 0 < maxiter , 'Invalid maxiters!'
        self.__maxiter = maxiter
    def __enter__ ( self ) :
        self.__prev = RRMMO.DefaultMaxIterations ()
        RRMMO.SetDefaultMaxIterations ( self.__maxiter )
        return self
    def __exit__ ( self , *_ ) :
        RRMMO.SetDefaultMaxIterations ( self.__prev  )

# =============================================================================
## @class TypeAlgo
#  Context manager to play with Default minimized type and algorithm 
#  @code 
#  with TypeAlgo ( 'Miniit2' , 'Migrad' ) : ...
#  @endcode 
#  @see ROOT::Math::MinimizerOptions
#  @see ROOT::Math::MinimizerOptions::DefaultMinimizerType
#  @see ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo
#  @see ROOT::Math::MinimizerOptions::SetDefaultMiminimizer
class TypeAlgo(object) :
    """Context manager to play with default max iterations

    >>> with TypeAlgo ( 'Miniit2' , 'Migrad' ) : ...
    
    - see ROOT::Math::MinimizerOptions
    - see ROOT::Math::MinimizerOptions::DefaultMinimizerType
    - see ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo
    - see ROOT::Math::MinimizerOptions::SetDefaultMimimizer
    """
    def __init__ ( self , mintype , minalgo = '' ) :
        assert isinstance ( mintype , string_types ) , 'Invalid mintype!'
        assert isinstance ( minalgo , string_types ) , 'Invalid minalgo!'
        self.__typealgo = mintype, minalgo        
    def __enter__ ( self ) :
        self.__prev = RRMMO.DefaultMinimizerType() , RRMMO.DefaultMinimizerAlgo() 
        RRMMO.SetDefaultMimimizer ( *self.__typealgo )
        return self
    def __exit__ ( self , *_ ) :
        RRMMO.SetDefaultMinimizer ( *self.__prev  )




# =============================================================================
## context manager to define default print level 
def print_level ( level ) :
    """ Context manager to define default print level
    """
    return PrintLevel ( level )

# =============================================================================
## context manager to define default strategy
def strategy ( strat ) :
    """ Context manager to define default strategy
    """
    return Strategy ( strat )

# =============================================================================
## context manager to define default tolerance 
def tolerance  ( tol ) :
    """ Context manager to define default tolerance
    """
    return Tolerance ( tol )

# =============================================================================
## context manager to define default max function calls 
def max_calls  ( mx ) :
    """ Context manager to define default max function calls 
    """
    return MaxFunctionCalls ( mx )

# =============================================================================
## context manager to define default max iterations 
def max_iterations  ( mx ) :
    """ Context manager to define default max iterations  
    """
    return MaxIterations ( mx )


# =============================================================================
## context manager to define default minimizer type and algorithm 
def type_algo  ( mtype , malgo = ''  ) :
    """ Context manager to define default minimizer type and algorithm 
    """
    return TypeAlgo ( mtype , malgo  )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
