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
    'gslCount'           , ## context manager to count  GSL errore
    'gslError'           , ## context manager to print  GSL errors 
    'gslException'       , ## context manager to turn   GSL errors into C++/Python exceptions
    'GslIgnore'          , ## context manager to ignore GSL errore
    'GslError'           , ## context manager to print  GSL errors 
    'GslCount'           , ## context manager to count  GSL errors 
    'GslException'       , ## context manager to turn   GSL errors into C++/Python exceptions
    'setHandler'         , ## use ``global'' GSL handler 
    'useHandler'         , ## ditto 
    )
# =============================================================================
from ostap.core.core import Ostap 
# =============================================================================

## helper base class/context manager  
class ErrHandler(object) :
    def __init__  ( self , force = True ) :
        self.err_handler = None
        self.__force     = True if force else False 
    def __enter__ ( self ) :
        self.err_handler = self.handler ( self.force )
        return self
    def __exit__  ( self , *_ ) :
        if self.err_handler : del self.err_handler
        self.err_handler = None
    @property
    def force ( self ) :
        """``force'' : force usge of error handler?"""
        return self.__force
    
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
    def __init__ ( self , force = True ) :
        self.handler = Ostap.Utils.GslIgnore
        super(GslIgnore,self).__init__( force )

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
    def __init__ ( self , force = True ) :
        self.handler = Ostap.Utils.GslError   
        super(GslError,self).__init__(  force )

# =============================================================================
## @class GslCount
#  Simple context manager to count GSL errors to stderr 
#  @code
#  with GslCount() :
#      ... do something 
#  @endcode 
class GslCount(ErrHandler) :
    """Simple context manager to count GSL errors 
    >>> with GslCount() :
    >>>    ... do something...
    """
    def __init__ ( self , force = True ) :
        self.handler = Ostap.Utils.GslCount   
        super(GslCount,self).__init__(  force )
        
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
    def __init__ ( self , force = True ) : 
        self.handler = Ostap.Utils.GslException 
        super(GslException,self).__init__( force )

# =============================================================================
## Simple context manager to ignore all GSL errors
#  @code
#  with gslIgnore() :
#      ... do something 
#  @endcode 
def gslIgnore   ( force = True ) :
    """Simple context manager to ignore all GSL errors
    >>> with gslIgnore() :
    >>>    ... do something...
    """
    return GslIgnore ( force )

# =============================================================================
## Simple context manager to count GSL errors
#  @code
#  with gslCount() :
#      ... do something 
#  @endcode 
def gslCount   ( force = True ) :
    """Simple context manager to count GSL errors
    >>> with gslCount () :
    >>>    ... do something...
    """
    return GslCount ( force )

# =============================================================================
## Simple context manager to print GSL errors to stderr 
#  @code
#  with gslError() :
#      ... do something 
#  @endcode 
def gslError    ( force = True ) :
    """Simple context manager to print GSL errors to stderr
    >>> with gslError() :
    >>>    ... do something...
    """
    return GslError ( force )

# =============================================================================
## Simple context manager to turn GSL errors into C++/Python exceptions 
#  @code
#  with gslException() :
#      ... do something 
#  @endcode 
def gslException ( force = True ) :
    """Simple context manager to turn GSL Errors into C++/Python exceptions 
    >>> with glException() :
    >>>    ... do something...
    """
    return GslException ( force )


# =============================================================================
_global_gsl_handler = [] 
def _setHandler ( handler ) :
    
    global _global_gsl_handler
    while  _global_gsl_handler : 
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
def setHandler ( handler , force = True ) :
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
            _setHandler ( Ostap.Utils.GslIgnore    ( force ) ) 
            logger.debug('Global GSL error Handler: Ignore all GLS errors') 
        elif hl in ( 'error' , 'print' ) : 
            _setHandler ( Ostap.Utils.GslError     ( force ) ) 
            logger.debug('Global GSL error Handler: print all GLS errors to stderr') 
        elif hl in ( 'count' , 'count' ) : 
            _setHandler ( Ostap.Utils.GslCount     ( force ) ) 
            logger.debug('Global GSL error Handler: count all GLS errors') 
        elif hl in ( 'exception' , 'raise' , 'throw' ) : 
            _setHandler ( Ostap.Utils.GslException ( force ) )
            logger.debug('Global GSL error Handler: convert GLS errors to C++/Python exceptions')
        else :
            raise TypeError ( 'Unknown handler type %s' % handler )
    elif isinstance ( handler , Ostap.Utils.GslError ) :
        _setHandler ( handler  , force )
        logger.debug('Global Eror Handler: %s' % handler ) 
    elif issubclass ( handler , Ostap.Utils.GslError ) :
        h = _setHandler ( handler ( force ) )
        logger.debug('Global Eror Handler: %s' % h )         
    else : 
        raise TypeError('Unknown handler type %s' % handler )


## ditto 
useHandler = setHandler 

# =============================================================================
# catch GSL errors from C++ and print summary table at exit 
# =============================================================================
import atexit 
@atexit.register
def print_gsl_errors () :
    """Catch GSL errors from C++ and print the summary table at exit
    """
    
    gsl_cnt =  Ostap.Utils.GslCount
    if 0    == gsl_cnt.size() : return  ## No GSL errors 
    
    ## get the summary 
    table = gsl_cnt.table()
    rows  = [] 
    for tline in table :
        
        try: 
            n , code , msg , reason , file , line = tline
            code = int ( code )
            n    = int ( n    )
        except :
            logger.warning ('print_gs_errors: failure to decode line: %s, skip it!' % str ( tline ) )
            continue 
        
        row = '%4d' % n , '%3d:%s' % ( code , msg ) , reason , file , line 
        rows.append ( row ) 
        
    if rows :
        
        from   ostap.logger.logger import getLogger
        logger = getLogger ( 'ostap.utils.gsl' )
        
        rows = [  ( '#' , 'error' , 'reason' , 'file', 'line' ) ] + rows 
        
        title = 'Summary of GSL errors'
        
        import ostap.logger.table     as     T
        from   ostap.logger.colorized import attention 
        logger.error ( '%s\n%s' % ( attention ( title ) , T.table ( rows , title = title, prefix = '# ' , alignment = 'ccccl') ) ) 
        
    ## clear the errors 
    gsl_cnt.clear() 
    del gsl_cnt 
    

# =============================================================================
if '__main__' == __name__ :
    
    
    from ostap.utils.docme import docme
    from  ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.utils.gsl' )
    docme ( __name__ , logger = logger )

    setHandler ( 'Error' , force = True )
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
            
    with gslCount ( force = True ) :
        for i in range ( 20 ) :
            Ostap.Math.psi ( -1 ) 
        
    logger.info ( 'Active handlers %s' % _global_gsl_handler ) 
    del   _global_gsl_handler[:]
    logger.info ( 'Active handlers %s' % _global_gsl_handler ) 
    logger.info ( 80*'*' )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
    
    
