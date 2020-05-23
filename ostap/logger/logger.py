#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright Ostap developers
# =============================================================================
"""Ostap simple logger.
Based on the logging of the Gaudi software project of CERN:
- Simple control (global and local)  over logging threshold.
 Primitive utilities for colorized logging.
"""
# =============================================================================
import logging, os, sys  
# =============================================================================
__all__ = (
    'getLogger'      , ## get (configured) logger
    'setLogging'     , ## set disable level according to MSG.Level
    'LogLevel'       , ## context manager to control output level 
    'logLevel'       , ## helper function to control output level
    'logVerbose'     , ## helper function to control output level
    'logDebug'       , ## helper function to control output level
    'logInfo'        , ## helper function to control output level
    'logWarning'     , ## helper function to control output level
    'logError'       , ## helper function to control output level
    'logColor'       , ## context manager to switch on  color logging locally  
    'logNoColor'     , ## context manager to switch off color logging locally  
    'noColor'        , ## context manager to switch off color logging locally  
    'keepColor'      , ## context manager preserve to preserve coloring
    'make_colors'    , ## force colored logging 
    'reset_colors'   , ## reset colored logging
    ##
    'ALL', 'VERBOSE', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'FATAL' ,
    )
# =============================================================================
# Message levels   (a'la Gaudi) 
ALL     = 0
VERBOSE = 1
DEBUG   = 2
INFO    = 3
WARNING = 4
ERROR   = 5
FATAL   = 6
# =============================================================================
## some manipulations with logging module
if not hasattr ( logging , 'VERBOSE' ) : logging.VERBOSE = 5
# =============================================================================
## some manipulations with logging module
if not hasattr ( logging , 'VERBOSE' ) : logging.VERBOSE = 5

# =============================================================================
## Log message with severity 'VERBOSE'
def _verbose1_(self, msg, *args, **kwargs):
    """Log 'msg % args' with severity 'VERBOSE'.
    """
    if self.isEnabledFor(logging.VERBOSE):
        self._log(logging.VERBOSE, msg, args, **kwargs)
        
# =============================================================================
## Log message with severity 'VERBOSE'
def _verbose2_(msg, *args, **kwargs):
    """Log a message with severity 'VERBOSE' on the root logger.
    """
    if not logging.root.handlers : logging.basicConfig()
    logging.root.verbose (msg, *args, **kwargs)

# =============================================================================
# add method 'verbose' to logger 
logging.Logger.verbose = _verbose1_

# =============================================================================
## add method 'verbose' to root logger 
logging.verbose        = _verbose2_
# =============================================================================
## convert MSG::Level into logging level 
def setLogging ( output_level ) :
    """Convert MSG::Level into logging level 
    """
    if   FATAL   <= output_level : logging.disable ( logging.FATAL   - 1 )
    elif ERROR   <= output_level : logging.disable ( logging.ERROR   - 1 )
    elif WARNING <= output_level : logging.disable ( logging.WARNING - 1 )
    elif INFO    <= output_level : logging.disable ( logging.INFO    - 1 )
    elif DEBUG   <= output_level : logging.disable ( logging.DEBUG   - 1 )
    elif VERBOSE <= output_level : logging.disable ( logging.VERBOSE - 1 )

# =============================================================================
## define standard logging names
logging_levels = { logging.CRITICAL : 'FATAL'   ,
                   logging.WARNING  : 'WARNING' ,
                   logging.DEBUG    : 'DEBUG  ' ,
                   logging.INFO     : 'INFO   ' ,
                   logging.ERROR    : 'ERROR  ' ,
                   logging.VERBOSE  : 'VERBOSE' }
for a in logging_levels : logging.addLevelName ( a ,  logging_levels[a]  )
# =============================================================================
logging_format      = '# %(name)-32s %(levelname)-7s %(message)s'
logging_file_format = '# %(asctime)s %(name)-32s %(levelname)-7s %(message)s'
logging_date_format = "%Y-%m-%d %H:%M:%S" 

# =============================================================================
## The basic configuration 
logging.basicConfig (
    level    = logging.INFO   ,
    format   = logging_format )

# =============================================================================
## get configured logger
#  @code
#  logger1 = getLogger ( 'LOGGER1' )
#  logger2 = getLogger ( 'LOGGER2' , level = logging.INFO )
#  @endcode 
def getLogger ( name   = 'ostap' ,
                format = ''      ,
                level  = None    ,
                stream = None    ) :
    
    """Get the proper logger
    >>> logger1 = getLogger ( 'LOGGER1' )
    >>> logger2 = getLogger ( 'LOGGER2' , level = logging.INFO )
    """
    #
    logger = logging.getLogger ( name )
    logger.propagate = True  ## ???

    ## if not logger.handlers :         
    ##     lh  = logging.StreamHandler ( stream ) if stream else logging.StreamHandler ()            
    ##     fmt = logging.Formatter     ( fmt    )
    ##     lh  . setFormatter          ( fmt    )
    ##     logger.addHandler           ( lh     )

    ## redefine log-level if needed 
    if level and level != logger.level :
        logger.setLevel ( level  )
    #
    return logger

# =============================================================================
## @class LogLevel
#  Temporarily enable/disable certain logger levels
#  @code
#  with LogLevel( logging.CRITICAL ) :
#       ...do something... 
#  @endcode
class LogLevel(object) :
    """Temporarily enable/disable certain logger levels
    >>> with LogLevel( logging.CRITICAL ) :
    ...  do something here ...
    """
    def __init__  ( self , level = logging.INFO - 1 ) :
        self.new_level = level 
        self.old_level = logging.root.manager.disable

    ## context manager: ENTER 
    def __enter__ ( self ) :
        self.old_level = logging.root.manager.disable
        logging.disable ( self.new_level )
        return self

    ## context manager: EXIT 
    def __exit__ ( self , *_ ) :        
        logging.disable ( self.old_level )
        
# =============================================================================
#  Temporarily enable/disable certain logger levels
#  @code
#  with logLevel( logging.CRITICAL ) :
#       ...do something... 
#  @endcode
def logLevel ( level = logging.INFO - 1 ) :
    """Temporarily enable/disable certain logger levels
    >>> with logLevel( logging.CRITICAL ) :
    >>>  ...do something...
    """
    return LogLevel ( level )

# =============================================================================
#  Temporarily enable/disable all loggers with level less then DEBUG 
#  @code
#  with logVerbose() :
#       ...do something... 
#  @endcode
def logVerbose () :
    """Temporarily disable all loggers with level less then INFO 
    >>> with logVerbose() :
    >>>  ...do something...
    """    
    return logLevel (  logging.VERBOSE   - 1 )

# =============================================================================
#  Temporarily enable/disable all loggers with level less then DEBUG 
#  @code
#  with logInfo() :
#       ...do something... 
#  @endcode
def logDebug   () :
    """Temporarily disable all loggers with level less then INFO 
    >>> with logDebug() :
    >>>  ...do something...
    """    
    return logLevel (  logging.DEBUG   - 1 )
# =============================================================================
#  Temporarily enable/disable all loggers with level less then INFO 
#  @code
#  with logInfo() :
#       ...do something... 
#  @endcode
def logInfo    () :
    """Temporarily disable all loggers with level less then INFO 
    >>> with logInfo() :
    >>>  ...do something...
    """
    return logLevel (  logging.INFO    - 1 )

# =============================================================================
#  Temporarily enable/disable all loggers with level less then WARNING
#  @code
#  with logInfo() :
#       ...do something... 
#  @endcode
def logWarning () : 
    """Temporarily disable all loggers with level less then WARNING
    >>> with logWarning() :
    >>>  ...do something...
    """   
    return logLevel (  logging.WARNING - 1 )

# =============================================================================
#  Temporarily enable/disable all loggers with level less then ERROR 
#  @code
#  with logError() :
#       ...do something... 
#  @endcode
def logError   () :
    """Temporarily disable all loggers with level less then ERROR
    >>> with logWarning() :
    >>>  ...do something...
    """       
    return logLevel (  logging.ERROR   - 1 )


# =============================================================================
## BASIC   colorization
# =============================================================================
from ostap.logger.colorized import ( with_colors    ,
                                     colored_string ,
                                     attention      ,
                                     allright       ,
                                     infostr        ,
                                     isatty         ,
                                     decolorize     )
# =============================================================================
__colored_logger = []
# =============================================================================
## reset colorization of logging 
def reset_colors () :
    """Reset colorization of logging 
    >>> reset_colors()
    """
    for a in logging_levels :
        logging.addLevelName ( a ,  logging_levels [ a ] )
    #
    while __colored_logger :
        __colored_logger.pop()

    from ostap.logger.colorized import set_with_colors
    set_with_colors ( False )
    return with_colors()

# =============================================================================
## make colors 
def make_colors () :
    """Colorize logging
    """

    if __colored_logger : return
    
    from ostap.logger.colorized import set_with_colors
    set_with_colors ( True )

    if not with_colors () : return
    
    def  makeName ( level , fg = None  , bg = None , blink = False , underline = False , bgb = False , fgb = False ) :

        name = logging.getLevelName ( level )
        bold = fg is None and bg is None and not uderline 
        bold = True
        return colored_string ( name , fg , bg , bold , blink , underline , fg_bright = fgb , bg_bright = bgb ) 

    from ostap.logger.colorized import RED , BLUE , YELLOW , GREEN , WHITE
    
    logging.addLevelName ( logging.CRITICAL ,  makeName ( logging.CRITICAL , fg = RED    , bg  = BLUE   , blink     = True ) )
    logging.addLevelName ( logging.WARNING  ,  makeName ( logging.WARNING  , fg = RED    , bg  = YELLOW , underline = True , bgb = True ) )
    logging.addLevelName ( logging.ERROR    ,  makeName ( logging.ERROR    , fg = YELLOW , bg  = RED    , blink     = True , bgb = True , fgb = True ) )
    logging.addLevelName ( logging.INFO     ,  makeName ( logging.INFO     , bg = BLUE   , fg  = WHITE  ) )
    logging.addLevelName ( logging.DEBUG    ,  makeName ( logging.DEBUG    , bg = GREEN  , fg  = WHITE  ) )

    __colored_logger.append ( 1 ) 
    return with_colors() 

# =============================================================================
## @class ColorLogging
#  Simple context manager to swicth on coloring
#  @code
#  with ColorLogging():
#      ... do something ... 
#  @endcode 
class ColorLogging(object) :
    """Simple context manager to swith on coloring
    
    >>> with ColorLogging() :
    ...     do something ... 
    """
    def __init__  ( self , color = True ) :
        self.color = color 
        
    def __enter__ ( self ) :
        self.with_color = with_colors() 
        if   self.color      and not self.with_color  : make_colors ()
        elif self.with_color and not self.color       : reset_colors ()
        return self
    
    def __exit__  ( self , *_ ) :
        if   self.color      and not self.with_color  : reset_colors ()
        elif self.with_color and not self.color       : make_colors ()

# =============================================================================
## simple context manager to switch on color logging 
#  @code
#  with logColor() :
#      ... do something ... 
#  @endcode 
def logColor ( color = True ) :
    """Simple context manager to switch on coloring
    
    >>> with logColor () :
    ...     do something ... 
    """
    return ColorLogging ( color )

# =============================================================================
## simple context manager to switch off color logging 
#  @code
#  with logNoColor() :
#      ... do something ... 
#  @endcode 
def logNoColor () :
    """Simple context manager to switch on coloring
    
    >>> with logNoColor () :
    ...     do something ... 
    """
    return ColorLogging ( False )

# =============================================================================
## simple context manager to switch off color logging 
#  @code
#  with noColor() :
#      ... do something ... 
#  @endcode 
def noColor () :
    """Simple context manager to switch on coloring
    
    >>> with noColor () :
    ...     do something ... 
    """
    return ColorLogging ( False )


# =============================================================================
## @class KeepColorLogging
#  Simple context manager to preserve coloring
#  @code
#  with KeepColorLogging():
#      ... do something ... 
#  @endcode 
class KeepColorLogging(object) :
    """Simple context manager to preserve coloring
    
    >>> with KeepColorLogging() :
    ...     do something ... 
    """
    def __enter__ ( self ) :
        self.with_color = with_colors() 
        return self
    
    def __exit__  ( self , *_ ) :
        if   self.with_color and not with_colors()   :  make_colors ()
        elif with_colors()   and not self.with_color : reset_colors ()

# =============================================================================
## simple context manager to preserve color logging 
#  @code
#  with keepColor() :
#      ... do something ... 
#  @endcode 
def keepColor () :
    """Simple context manager to preserve color logging 
    
    >>> with keepColor () :
    ...     do something ... 
    """
    return KeepColorLogging ()

# =============================================================================
# Actions!
# =============================================================================

## reset colors
if isatty  () : make_colors()

## define the default logging thresholds as 'INFO'
setLogging ( 3 )

# =============================================================================
# Log file?
# =============================================================================
log_file = os.getenv ( 'OSTAP_LOGFILE' , '' )
if log_file : 

    ## set buffering to be 1-line and decolorize the output   
    class LogHandler(logging.FileHandler) :
        def __init__(self, filename, mode='w', encoding=None, delay=0):
            logging.FileHandler.__init__ ( self , filename , mode , encoding, delay ) 
        def _open(self):
            """
            Open the current base file with the (original) mode and encoding.
            Return the resulting stream.
            """
            stream = open(self.baseFilename, self.mode, buffering = 1 )
            return stream
        
        def emit(self, record):
            """Emit an ddecolorize the record
            """
            lname = logging_levels.get ( record.levelno , '' )
            if not lname : lname = '%s' % record.levelno
            record.levelname = lname
            if with_colors () : record.msg = decolorize ( record.msg ) 
            return logging.FileHandler.emit ( self , record ) 
    
    loglev = os.getenv ( 'OSTAP_LOGLEVEL' , '%s' % logging.INFO )
    try :
        loglev = int ( loglev )
        if not loglev in logging_levels : loglev = logging.INFO 
    except :
        loglev = logging.INFO
    log_handler = LogHandler ( log_file , mode = 'w' )
    log_handler.setLevel ( loglev ) 
    formatter   = logging.Formatter ( logging_file_format , logging_date_format )
    log_handler.setFormatter ( formatter   ) 
    logging.root.addHandler  ( log_handler )

if log_file :

    logger = getLogger('ostap.logger.logger')
    func   = lambda : logger.info ( 'Log-file is %s' %  log_file )
    func () 
    import atexit    
    atexit.register ( func )

# =============================================================================
if __name__ == '__main__' :

    setLogging ( 0 )
    
    logger = getLogger ( 'ostap.logger.logger')

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
 
    logger.verbose  ( 'This is VERBOSE  message'  ) 
    logger.debug    ( 'This is DEBUG    message'  ) 
    logger.info     ( 'This is INFO     message'  ) 
    logger.warning  ( 'This is WARNING  message'  ) 
    logger.error    ( 'This is ERROR    message'  ) 
    logger.fatal    ( 'This is FATAL    message'  ) 
    logger.critical ( 'This is CRITICAL message'  ) 

    with logColor() : 
        
        logger.verbose  ( 'This is VERBOSE  message'  ) 
        logger.debug    ( 'This is DEBUG    message'  ) 
        logger.info     ( 'This is INFO     message'  )
        logger.warning  ( 'This is WARNING  message'  ) 
        logger.error    ( 'This is ERROR    message'  ) 
        logger.fatal    ( 'This is FATAL    message'  ) 
        logger.critical ( 'This is CRITICAL message'  )

        with noColor () : 
            logger.verbose  ( 'This is VERBOSE  message'  ) 
            logger.debug    ( 'This is DEBUG    message'  ) 
            logger.info     ( 'This is INFO     message'  )
            logger.warning  ( 'This is WARNING  message'  ) 
            logger.error    ( 'This is ERROR    message'  ) 
            logger.fatal    ( 'This is FATAL    message'  ) 
            logger.critical ( 'This is CRITICAL message'  )
            
        logger.verbose  ( 'This is VERBOSE  message'  ) 
        logger.debug    ( 'This is DEBUG    message'  ) 
        logger.info     ( 'This is INFO     message'  )
        logger.warning  ( 'This is WARNING  message'  ) 
        logger.error    ( 'This is ERROR    message'  ) 
        logger.fatal    ( 'This is FATAL    message'  ) 
        logger.critical ( 'This is CRITICAL message'  )

    with keepColor() :
        logger.verbose  ( 'This is VERBOSE  message'  ) 
        logger.debug    ( 'This is DEBUG    message'  ) 
        logger.info     ( 'This is INFO     message'  ) 
        logger.warning  ( 'This is WARNING  message'  )

        make_colors()
        
        logger.error    ( 'This is ERROR    message'  ) 
        logger.fatal    ( 'This is FATAL    message'  ) 
        logger.critical ( 'This is CRITICAL message'  ) 
        
    logger.verbose  ( 'This is VERBOSE  message'  ) 
    logger.debug    ( 'This is DEBUG    message'  ) 
    logger.info     ( 'This is INFO     message'  ) 
    logger.warning  ( 'This is WARNING  message'  ) 
    logger.error    ( 'This is ERROR    message'  ) 
    logger.fatal    ( 'This is FATAL    message'  ) 
    logger.critical ( 'This is CRITICAL message'  ) 
    logger.info ( 80*'*'  )

# =============================================================================
##                                                                      The END 
# =============================================================================
