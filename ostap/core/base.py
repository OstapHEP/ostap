#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/base.py
#  Some base objects for ostap 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
""" Core/base objects for ostap 
"""
# ============================================================================= 
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    ##
    'cpp'                 , ## Global C++ namespace 
    'std'                 , ## C++ namespace std 
    'Ostap'               , ## C++ namespace Ostap 
    'ROOTIgnore'          , ## control ROOT   verbosity, suppress ROOT errors
    'RooSilent'           , ## control RooFit verbosity, suppress ROOT errors
    'rootException'       , ## context manager to perform ROOT Error -> C++/Python exception    
    'RootError2Exception' , ## context manager to perform ROOT Error -> C++/Python exception
    ## 
    'rooSilent'           , ## control RooFit verbosity
    'roo_silent'          , ## control RooFit verbosity 
    'rootError'           , ## control ROOT verbosity 
    'rootWarning'         , ## control ROOT verbosity
    ## 
    'valid_pointer'       , ## valid C++ pointer?
)
# =============================================================================
from  ostap.utils.basic import NoContext 
import ROOT, cppyy 
# ============================================================================= 
## get global C++ namespace
cpp   = cppyy.gbl
# =============================================================================
## C++ namespace std 
std   = cpp.std
# =============================================================================
## C++ namespace Ostap
Ostap = cpp.Ostap 
# =============================================================================
## Simple context manager to control RooFit evrbosity
# =============================================================================
## very simple context manager to suppress RooFit printout
#
#  @code
#
#  >>> with rooSilent( 4 , False ) :
#  ...        some_RooFit_code_here()
#
#  @endcode
#  @see RooMgsService
#  @see RooMgsService::globalKillBelow
#  @see RooMgsService::silentMode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-09
class RooSilent(object) :
    """ Very simple context manager to suppress RooFit printout
    
    >>> with rooSilent( 4 , False ) :
    ...        some_RooFit_code_here ()
    
    """
    ## constructor
    #  @param level  (INPUT) print level 
    #  @param silent (print level 
    # 
    def __init__ ( self   ,
                   level  = ROOT.RooFit.ERROR ,
                   silent = True              ) :
        """ Constructor
        @param level  (INPUT) print level 
        @param silent (print level 
        
        >>> with rooSilent( ROOT.RooFit.ERROR , True  ) :
        ...        some_RooFit_code_here ()
        
        
        >>> with rooSilent( ROOT.RooFit.INFO , False  ) :
        ...        some_RooFit_code_here ()
                
        """
        #
        if level > ROOT.RooFit.FATAL : level = ROOT.RooFit.FATAL 
        if level < ROOT.RooFit.DEBUG : level = ROOT.RooFit.DEBUG 
        #
        self.__roo_level   = level 
        self.__roo_silent  = True if silent else False  

        svc = ROOT.RooMsgService.instance()
        self.__prev_level  = svc.globalKillBelow  () 
        self.__prev_silent = svc.silentMode       () 

    # =========================================================================
    ## context manager : ENTER 
    def __enter__ ( self ) :
        """ Contex manager: ENTER 
        """
        svc = ROOT.RooMsgService.instance()
        self.__prev_level  = svc.globalKillBelow  () 
        self.__prev_silent = svc.silentMode       () 
        ##
        svc.saveState           ()
        svc.setGlobalKillBelow  ( self.__roo_level  )
        svc.setSilentMode       ( self.__roo_silent )
        ## 
        return self
    
    # =========================================================================
    ## context manager: EXIT 
    def __exit__ ( self , *_ ) : 
        """ Contex manager: EXIT
        """
        svc = ROOT.RooMsgService.instance()
        svc.setSilentMode      ( self.__prev_silent )
        svc.setGlobalKillBelow ( self.__prev_level  )
        svc.restoreState       ()
        
# =============================================================================
## Very simple context manager to suppress ROOT printout
#  @code
#  >>> with ROOTIgnore( ROOT.kError + 1 ) : some_ROOT_code_here()
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-30
class ROOTIgnore ( RooSilent ) :
    """ Very simple context manager to suppress ROOT printout
    >>> with ROOTIgnore ( ROOT.kError + 1 ) : some_ROOT_code_here()
    """
    ## constructor
    #  @param level  (INPUT) print level 
    #  @param silent (print level 
    # 
    def __init__ ( self , level ) :
        """ Constructor:        
        >>> with rootError   () : some_ROOT_code_here()
        >>> with rootWarning () : some_ROOT_code_here()
        """
        #
        level = int ( level )
        if   ROOT.kUnset == level : level = ROOT.kInfo 
        elif level > ROOT.kFatal  : level = ROOT.kFatal 
        elif level < ROOT.kPrint  : level = ROOT.kPrint 
        #
        level = max ( level , ROOT.kPrint )
        level = min ( level , ROOT.kFatal )
        ## 
        self._level   = level 
        
        if   ROOT.kBreak   <= level : rlevel = ROOT.RooFit.FATAL
        elif ROOT.kError   <= level : rlevel = ROOT.RooFit.ERROR
        elif ROOT.kWarning <= level : rlevel = ROOT.RooFit.WARNING
        elif ROOT.kInfo    <= level : rlevel = ROOT.RooFit.INFO  
        else                        : rlevel = ROOT.RooFit.DEBUG 

        silent = ROOT.kWarning <= level
        
        RooSilent.__init__ ( self , rlevel , silent ) 
        
    # =========================================================================
    ## context manager: ENTER 
    def __enter__ ( self ) :
        """ The actual context manager: ENTER
        """
        ## ROOT 
        self._old = int ( ROOT.gErrorIgnoreLevel ) 
        if self._old != self._level :
            groot = ROOT.ROOT.GetROOT()
            if groot : groot.ProcessLine ( "gErrorIgnoreLevel= %d ; " % self._level ) 
        ## RooFit 
        RooSilent.__enter__ ( self )
        return self 
    
    # =========================================================================
    ## context manager: EXIT 
    def __exit__ ( self , *_ ) : 
        """ The actual context manager: EXIT
        """
        ## RooFit 
        RooSilent.__exit__ ( self , *_ )
        ## ROOT 
        if self._old != int ( ROOT.gErrorIgnoreLevel )  :
            groot = ROOT.ROOT.GetROOT()            
            if groot : groot.ProcessLine ( "gErrorIgnoreLevel= %d ; " % self._old ) 


# =============================================================================
## very simple context manager to suppress ROOT printout
#  @code
#  >>> with rootError () : some_ROOT_code_here()
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-30
def rootError   ( level = 1 ) :
    """ Very simple context manager to suppress ROOT printout
    >>> with rootError () : some_ROOT_code_here()
    """
    return ROOTIgnore ( ROOT.kError   + level )

# =============================================================================
## very simple context manager to suppress ROOT printout
#  @code
#  >>> with rootError () : some_ROOT_code_here()
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-30
def rootWarning ( level = 1 ) :
    """ Very simple context manager to suppress ROOT printout
    >>> with rootWarning () : some_ROOT_code_here()
    """
    return ROOTIgnore ( ROOT.kWarning + level )


# =============================================================================
## very simple context manager to suppress RooFit printout
#
#  @code
#
#  >>> with rooSilent( 4 , False ) :
#  ...        some_RooFit_code_here()
#
#  @endcode
#  @see RooMgsService
#  @see RooMgsService::globalKillBelow
#  @see RooMgsService::silentMode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-09
def rooSilent ( level = ROOT.RooFit.ERROR , silent = True ) :
    """ Very simple context manager to suppress RooFit printout
    >>> with rooSilent( 4 , False ) :
    ...        some_RooFit_code_here()    
    """
    return RooSilent ( level , silent ) 

# =============================================================================
## helper context manager
#  @code
#
#  >>> with roo_silent( True ) : 
#  ...        some_RooFit_code_here()
#
#  @endcode
#  @see rooSilent
#  @see NoContex
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-09
def roo_silent ( silence , *args ) :
    """ Helper context manager#
    >>> with roo_silent ( True ) : 
    ...        some_RooFit_code_here()
    """
    return rooSilent ( *args ) if silence else NoContext() 

# =============================================================================

            
# =============================================================================
## helper context manager to activate ROOT Error -> Python exception converter 
#  @see Ostap::Utils::useErrorHandler
#  @see Ostap::Utils::ErrorSentry
#  @code
#  with RootError2Exception() :
#  .... do something here 
#  @endcode 
class RootError2Exception (object) :
    """ Helper context manager to activate ROOT Error -> Python exception converter
    >>> with RootError2Exception() :
    ...      do something here 
    """
    def __init__ ( self ) :
        self.e_handler  = Ostap.Utils.useErrorHandler 
        self.m_previous = False 

    ## context manager entry point  
    def __enter__ ( self ) :    
        self.m_previous = self.e_handler ( True ) 
        return self
    
    ## context manager exit point
    def __exit__ ( self , *_ ) :    
        if self.m_previous : self.e_handler ( False ) 
        self.m_previous = False 

# =============================================================================
## helper context manager to activate ROOT Error -> Python exception converter 
#  @see Ostap::Utils::useErrorHandler
#  @see Ostap::Utils::ErrorSentry
#  @code
#  with rootException () :
#  .... do something here 
#  @endcode
def rootException () :
    """ Helper context manager to activate ROOT Error -> Python exception converter
    #
    with rootException() :
    ... do something here 
    """
    return RootError2Exception()

# =============================================================================
with ROOTIgnore ( ROOT.kError ) :
    ## valid C++ pointer ? 
    _valid_pointer_ = Ostap.Utils.valid_pointer
    ## used by ROOT/RooFit ?  
    usedRootID = Ostap.Utils.usedRootID 

# =============================================================================
## Is it a valid C++ pointer?
#  @code
#  ptr = ...
#  print 'Is the pointer valid? %s'  % valid_pointer ( prt ) 
#  @endcode 
#  @see Ostap::Utils::valid_pointer 
def valid_pointer ( obj ) :
    """ Is it a valid C++ pointer?
    - see Ostap::Utils::valid_pointer 
    >>> ptr = ...
    >>> print 'Is the C++ pointer valid? %s'  % valid_pointer ( ptr ) 
    """
    r = _valid_pointer_ ( obj )
    return True if r else False

# =============================================================================
if '__main__' == __name__ :

    from ostap.logger.logger import getLogger 
    logger = getLogger( 'ostap.core.base' )

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# ============================================================================= 

