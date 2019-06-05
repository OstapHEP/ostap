#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file
# Trivial utility for documentation
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date   2016-12-02
# =============================================================================
"""Trivial utility for documentation purposes
>>> docme ( module )
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    'docme' , ## the only one useful symbol from this module
    )
# =============================================================================
## primitive types
from sys import version_info as python_version 
if python_version.major > 2 : 
    primitive =  int        , float , bool , str , bytes   , list , dict , set
else :
    primitive =  int , long , float , bool , str , unicode , list , dict , set
    
# =============================================================================
## simple function to allow coherent self-print for all ostap modules
#  @code
#  docme ( module )
#  @endcode
def docme( module , symbols = {} , logger = None ) :
    """Simple function to allow coherent self-print for all ostap modules
    >>> docme ( module ) 
    """
    if not logger :
        from ostap.logger.logger import getLogger
        if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.docme' )
        else                       : logger = getLogger( __name__            )
        
    if isinstance ( module , str ) :
        import sys
        mod = sys.modules.get( module , None ) 
        if not mod :
            logger.error ( "Can't get module %s" % module )
            return
        module = mod
        
    import inspect
    if not inspect.ismodule ( module ) :
        logger.error ( 'The module is *not* a module')
        return

    logger.info ( 80*'*' )

    ## the file name and the banner
    _file_ = module.__file__ if hasattr ( module , '__file__' ) else None
    from ostap import banner
    if _file_ : logger.info ( _file_ + '\n' + banner ) 
    else      : logger.info (                 banner ) 

    logger.info ( 80*'*' )

    ## own documentation 
    _doc_ = module.__doc__ if hasattr ( module , '__doc__' ) else None
    if _doc_ :
        logger.info ( _doc_ ) 
        logger.info ( 80*'*' )

    ## Author 
    _author_  = module.__author__  if hasattr ( module , '__author__' ) else None
    if _author_ :
        logger.info ( "Author  : %s" % _author_ ) 

    ## Version
    _version_ = module.__version__ if hasattr ( module , '__version__' ) else None
    if _version_ :
        logger.info ( "Version : %s" % _version_ ) 

    ## Date 
    _date_    = module.__date__    if hasattr ( module , '__date__' ) else None
    if _version_ :
        logger.info ( "Date    : %s" % _date_ ) 

    _done_ = set ()
    
    ## Public symbols 
    _all_     = module.__all__    if hasattr ( module , '__all__' ) else None
    if _all_ :
        logger.info ( "Symbols : %s" % list(_all_) ) 
        for key in _all_ :
            sym = getattr ( module , key )            
            if isinstance ( sym , primitive ) : continue
            if isinstance ( sym , (tuple, list, set) ) :
                logger.info ( "Symbol ``%s''\n# - %s" % ( key , str ( sym ) ) )                
            elif hasattr ( sym , '__doc__' ) and sym.__doc__ :
                d = sym.__doc__.replace( '\n' , '\n#' )
                logger.info ( "Symbol ``%s''\n# - %s" % ( key , d ) )
            else : 
                logger.info ( "``%s''" %        key ) 
            _done_.add ( sym ) 
        logger.info ( 80*'*' )

    _klasses_  = getattr ( module , '_decorated_classes_' , () )
    _klasses_  = set ( _klasses_ ) 
    if _klasses_ :
        logger.info ( "Decorated classes : %s" % [ getattr ( k ,'__name__' , k ) for k in _klasses_ ] )
        for sym in _klasses_ :
            if hasattr ( sym , '__doc__' ) and sym.__doc__ :
                d = sym.__doc__.replace( '\n' , '\n#' )
                logger.info ( "Decorated class ``%s''\n# - %s" % ( getattr ( sym , '__name__' , sym ) , d ) )
        
        logger.info ( 80*'*' )

    _methods_  = getattr ( module , '_new_methods_' , () )
    _methods_  = set ( _methods_ ) 
    if _methods_ :
        for sym in _methods_ :
            if hasattr ( sym , '__doc__' ) and sym.__doc__ :
                d = sym.__doc__.replace( '\n' , '\n#' )
                logger.info ( "New method\n# - %s" % d )        
        logger.info ( 80*'*' )
        
    if isinstance ( symbols , dict ) :
        
        for key in sorted ( symbols.keys () ) :
            if key.startswith('_')      : continue
            sym = symbols[key]
            if sym in _done_            : continue 
            if inspect.ismodule ( sym ) : continue 
            if hasattr ( sym , '__doc__' ) and sym.__doc__ and not isinstance ( sym , primitive ) : 
                d = sym.__doc__.replace( '\n' , '\n#' )
                logger.info ( "Symbol ``%s''\n# - %s" % ( key , d ) )
            else : 
                logger.info ( "``%s''" %        key ) 

        logger.info ( 80*'*' )

    elif isinstance ( symbols , ( list , set , tuple ) ) :

        for sym in symbols : 
            if sym in _done_ : continue
            if inspect.ismodule ( sym ) : continue 
            if hasattr ( sym , '__doc__' ) and sym.__doc__ and not isinstance ( sym , primitive ) : 
                d = sym.__doc__.replace( '\n' , '\n#' )
                logger.info ( "Symbol ``%s''\n# - %s" % ( key , d ) )
            else : 
                logger.info ( "``%s''" %        key ) 
                        
        logger.info ( 80*'*' )

# =============================================================================
if '__main__' == __name__ :

    docme ( __name__ ) 

# =============================================================================
# The END 
# =============================================================================

