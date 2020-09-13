#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some simple but useful utilities fro GSL Error handling 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
# =============================================================================
"""Module with some simple but useful utilities for GSL error handling 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    #
    'gslIgnore'          , ## context manager to ignore GSL errore
    'gslError'           , ## context manager to print  GSL errors 
    'gslException'       , ## context manager to turn   GSL errors into C++/Python exceptions
    'GslIgnore'          , ## context manager to ignore GSL errore
    'GslError'           , ## context manager to print  GSL errors 
    'GslException'       , ## context manager to turn   GSL errors into C++/Python exceptions
    'setHandler'         , ## use ``global'' GSL handler 
    'useHandler'         , ## ditto 
    )
# =============================================================================
## helper base class/context manager  
class ErrHandler(object) :
    def __init__  ( self ) :
        self.err_handler = None 
    def __enter__ ( self ) :
        self.err_handler = self.handler ()
        return self
    def __exit__  ( self , *_ ) :
        if self.err_handler : del self.err_handler
        self.err_handler = None
# =============================================================================
## @class GslIgnore
#  Simple context manager to ignore all GSL errors
#  @code
#  with GslIgnore() :
#      ... do something 
#  @endcode 
class GslIgnore(ErrHandler) :
    """Simple context manager to ignore all GSL errors
    >>> with GslIgnore() :
    >>>    ... do something...
    """
    def __init__ ( self ) :
        from ostap.core.core import Ostap 
        self.handler = Ostap.Utils.GslIgnore
        super(GslIgnore,self).__init__()
# =============================================================================
## @class GslError
#  Simple context manager to print GSL errors to stderr 
#  @code
#  with GslError() :
#      ... do something 
#  @endcode 
class GslError(ErrHandler) :
    """Simple context manager to print GSL errors to stderr
    >>> with GslError() :
    >>>    ... do something...
    """
    def __init__ ( self ) :
        from ostap.core.core import Ostap 
        self.handler = Ostap.Utils.GslError   
        super(GslError,self).__init__()
        
# =============================================================================
## @class GslException
#  Simple context manager to turn GSL errors into C++/Python exceptions 
#  @code
#  with GslException() :
#      ... do something 
#  @endcode 
class GslException (ErrHandler) :
    """Simple context manager to turn GSL Errors into C++/Python exceptions 
    >>> with GslException() :
    >>>    ... do something...
    """
    def __init__ ( self ) : 
        from ostap.core.core import Ostap 
        self.handler = Ostap.Utils.GslException 
        super(GslException,self).__init__()

# =============================================================================
## Simple context manager to ignore all GSL errors
#  @code
#  with gslIgnore() :
#      ... do something 
#  @endcode 
def gslIgnore   () :
    """Simple context manager to ignore all GSL errors
    >>> with gslIgnore() :
    >>>    ... do something...
    """
    return GslIgnore()

# =============================================================================
## Simple context manager to print GSL errors to stderr 
#  @code
#  with gslError() :
#      ... do something 
#  @endcode 
def gslError    () :
    """Simple context manager to print GSL errors to stderr
    >>> with gslError() :
    >>>    ... do something...
    """
    return GslError()

# =============================================================================
## Simple context manager to turn GSL errors into C++/Python exceptions 
#  @code
#  with gslException() :
#      ... do something 
#  @endcode 
def gslException () :
    """Simple context manager to turn GSL Errors into C++/Python exceptions 
    >>> with glException() :
    >>>    ... do something...
    """
    return GslException()


# =============================================================================
_global_gsl_handler = [] 
def _setHandler ( handler ) :
    global _global_gsl_handler
    while _global_gsl_handler : 
        _global_gsl_handler.pop() 
    if handler: _global_gsl_handler.append ( handler ) 
    return _global_gsl_handler

# =============================================================================
## Make use ``global'' GSL handler
#  @code
#  setHandler ( None        ) ## clean up global  handlers 
#  setHandler ( 'Ignore'    ) ## ignore all GSL erorrs 
#  setHandler ( 'Error'     ) ## print GSL errors to stderr and continue
#  setHandler ( 'Exception' ) ## convert GSL errors into C++/Python exceptions 
#  setHandler ( 'Raise'     ) ## ditto 
#  setHandler ( 'Throw'     ) ## ditto 
#  @endcode 
def setHandler ( handler ) :
    """Use ``global'' GSL handler
    >>> setGlobalHandler ( None        ) ## clean up global  handlers 
    >>> setGlobalHandler ( 'Ignore'    ) ## ignore all GSL erorrs 
    >>> setGlobalHandler ( 'Error'     ) ## print GSL errors to stderr and continue
    >>> setGlobalHandler ( 'Exception' ) ## convert GSL errors into C++/Python exceptions 
    >>> setGlobalHandler ( 'Raise'     ) ## ditto 
    >>> setGlobalHandler ( 'Throw'     ) ## ditto 
    """
    #
    from   ostap.logger.logger import getLogger
    logger = getLogger( 'ostap.utils.gsl' )
    #
    from ostap.core.core import Ostap 
    #
    global _global_gls_handler
    if   not  handler : _setHandler ( handler )
    elif isinstance ( handler , str ) : 
        hl = handler.lower()
        if 'ignore' == hl  :
            _setHandler ( Ostap.Utils.GslIgnore    () ) 
            logger.debug('Global GSL error Handler: Ignore all GLS errors') 
        elif hl in ( 'error' , 'print' ) : 
            _setHandler ( Ostap.Utils.GslError     () ) 
            logger.debug('Global GSL error Handler: print all GLS errors to stderr') 
        elif hl in ( 'exception' , 'raise' , 'throw' ) : 
            _setHandler ( Ostap.Utils.GslException () )
            logger.debug('Global GSL error Handler: convert GLS errors to C++/Python exceptions')
        else :
            raise TypeError ( 'Unknown handler type %s' % handler )
    elif isinstance ( handler , Ostap.Utils.GslError ) :
        _setHandler ( handler  )
        logger.debug('Global Eror Handler: %s' % handler ) 
    elif issubclass ( handler , Ostap.Utils.GslError ) :
        h = _setHandler ( handler () )
        logger.debug('Global Eror Handler: %s' % h )         
    else : 
        raise TypeError('Unknown handler type %s' % handler )

## ditto 
useHandler = setHandler 
# =============================================================================
if '__main__' == __name__ :
    
    from   ostap.logger.logger import getLogger
    if '__main__' ==  __name__ : logger = getLogger ( 'ostap.utils.gsl' )
    else                       : logger = getLogger ( __name__          )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    setHandler ( 'Error' )
    with gslIgnore() , gslException() :
        with gslError() :
            setHandler ( 'Exception' )
            setHandler ( 'Exception' )
            setHandler ( 'Exception' )
            setHandler ( 'Exception' )
            setHandler ( 'Ignore'    )
            setHandler ( 'Error'     )
            setHandler ( 'Error'     )
            setHandler ( 'Error'     )
            setHandler ( 'Error'     )
        
    logger.info ( 'Active handlers %s' % _global_gsl_handler ) 
    del   _global_gsl_handler[:]
    logger.info ( 'Active handlers %s' % _global_gsl_handler ) 
    logger.info ( 80*'*' )
    
# =============================================================================
# The END 
# =============================================================================
    
    
