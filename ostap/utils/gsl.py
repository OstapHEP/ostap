#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some simple but useful utilities fro GSL Error handling 
#   - timing
#   - memory
#   - profiling
#   - ... 
#
#  It is recommended to install psutil module 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
#  
# =============================================================================
"""Module with some simple but useful utilities for GSL erro handling 
- timing
- memory
- profiling
- etc

It is recommended to install psutil module 
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
    )
# =============================================================================
import ROOT,cppyy 
# =============================================================================
class ErrHandler(object) :
    def __init__  ( self ) :
        self.err_handler = None
    def __exit__  ( self , *_ ) :
        if self.err_handler : del self.err_handler
        self.err_handler = None
    def __del__   ( self ) :
        if self.err_handler : del self.err_handler

# =============================================================================
## @class GslIgnore
#  Simple context manager to ignore all GSL errors
#  @code
#  with GslIgnore() :
#      ... do something 
#  @endcode 
class GslIgnore   (ErrHandler) :
    """Simple context manager to ignore all GSL errors
    >>> with GslIgnore() :
    >>>    ... do something...
    """
    def __enter__ ( self ) :
        Ostap = cppyy.gbl.Ostap
        self.err_handler = Ostap.Utils.GslIgnore    ()
        return self
    
# =============================================================================
## @class GslError
#  Simple context manager to print GSL errors to stderr 
#  @code
#  with GslError() :
#      ... do something 
#  @endcode 
class GslError    (ErrHandler) :
    """Simple context manager to print GSL errors to stderr
    >>> with GslError() :
    >>>    ... do something...
    """
    def __enter__ ( self ) :
        Ostap = cppyy.gbl.Ostap
        self.err_handler = Ostap.Utils.GslError     ()
        return self
            
# =============================================================================
## @class GslException
#  Simple context manager to turn GSL errors into C++/Python exceptions 
#  @code
#  with GslException() :
#      ... do something 
#  @endcode 
class GslException (ErrHandler) :
    """Simple context manager to turn GSL Errors into C++/Python exceptions 
    >>> with GslError() :
    >>>    ... do something...
    """
    def __enter__ ( self ) :
        Ostap = cppyy.gbl.Ostap
        self.err_handler = Ostap.Utils.GslException ()
        return self

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
    >>> with GslError() :
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
    >>> with GslError() :
    >>>    ... do something...
    """
    return GslException()

# =============================================================================
if '__main__' == __name__ :
    
    from   ostap.logger.logger import getLogger, isatty 
    if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.gsl' )
    else                       : logger = getLogger( __name__          )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    logger.info ( 80*'*' ) 
    
# =============================================================================
# The END 
# =============================================================================
    
    
