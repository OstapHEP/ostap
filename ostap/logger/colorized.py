#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Simple colorization of strings 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
# =============================================================================
"""Simple colorization of strings"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    'colored_string'   , ## make a colored string
    'decolorize'       , ## decolorize the colored string 
    'attention'        , ## make "attention" string
    'allright'         , ## make "allright" string
    'infostr'          , ## just for information
    #
    'with_colors'      , ## Is colorization enabled?
    'set_with_colors'  , ## Enable/Disable colorization 
    )
# =============================================================================
import os, sys
# =============================================================================
# - is sys.stdout attached to terminal or not ?
# from ostap.utils.basic import isatty
def isatty ( stream = None ) :
    """Is the stream is attached to terminal?
    >>> stream = ...
    >>> if isatty( stream ) : print('Teminal!')
    >>> if isatty() : print('stdout is terminal!')
    """
    if not stream : stream = sys.stdout
    #
    try :
        return stream.isatty()
    except : pass 
    #
    try :
        return os.isatty ( stream.fileno() ) 
    except : pass
    #
    return False
# =============================================================================
## global flag to indicate if we use colors 
__with_colors__ = isatty ()   
# =============================================================================
## Is colorization enabled ? 
def with_colors() :
    """Is colorization enabled ?"""
    global __with_colors__
    return bool ( __with_colors__ ) and isatty() 
# =============================================================================
## Enable/disable colorization
def set_with_colors ( use ) :
    """Enable/disable colorization"""
    global __with_colors__
    __with_colors__ = bool ( use )
    return with_colors () 

# =============================================================================
## BASIC ASCII colors :
BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = list ( range ( 8 ) )
# =============================================================================
## provide colored string
#  @code
#  print colored_string ( 'Hello' , foreground = RED , background = YELLOW , bold = True )
#  @endcode
#  @see https://en.wikipedia.org/wiki/ANSI_escape_code#Colors
def colored_string ( what               ,
                     foreground = None  ,
                     background = None  ,
                     bold       = False ,
                     blink      = False ,
                     underline  = False ,
                     fg_bright  = False ,
                     bg_bright  = False ) :
    """
    >>> print colored_string ( 'Hello' , foreground = RED , background = YELLOW , bold = True , blink = True , underline = True )
    - see https://en.wikipedia.org/wiki/ANSI_escape_code#Colors
    """
    ## nothing to colorize or no coloring is activated
    if not what or not with_colors() or not isatty() : return what

    ## nothing to do 
    if ( foreground is None ) and ( background is None ) :
        if ( not bold ) and ( not blink ) and ( not underline ) : return what 

    RESET_SEQ = "\033[0m"

    COLOR_SEQ = "\033[%dm"
    BOLD_SEQ  = "\033[1m"    if bold      else ''
    BLINK_SEQ = "\033[5m"    if blink     else '' 
    ULINE_SEQ = "\033[4m"    if underline else ''
    
    ESCAPE    = '\033['

    keys = [] 
    
    if bold      : keys.append  ( '1' )
    if underline : keys.append  ( '4' )
    if blink     : keys.append  ( '5' )
    
    if foreground is None : pass 
    else : keys.append ( str ( (  90 if fg_bright else 30 ) + foreground % 8 ) ) 
    
    if background is None : pass 
    else : keys.append ( str ( ( 100 if bg_bright else 40 ) + background % 8 ) ) 

    if not keys : return what

    prefix = '\033[%sm' %  ( ';'.join ( k for k in keys ) )

    return '{prefix}{what}{reset}'.format ( prefix = prefix , what = what , reset = RESET_SEQ )


# =============================================================================
import re 
_decolor = re.compile ( r'\033\[\d{1,3}(;\d{1,3})*m' )
# =============================================================================
## decolorize the string: eliminate all colorization stuff  
def decolorize ( what ) :
    """Decolorize the string: eliminate all colorization stuff
    """
    return _decolor.sub ( '' , what ) if what else what 

# =====-=======================================================================
## attention!
def attention ( what ) :
    """Attention string """
    return colored_string ( what                ,
                            foreground = YELLOW ,
                            background = RED    ,
                            bold       = True   ,
                            blink      = True   ,
                            underline  = True   ,
                            fg_bright  = True   ,
                            bg_bright  = True   )

# =============================================================================
## allright 
def allright ( what ) :
    """Allright string """
    return colored_string ( what                ,
                            foreground = YELLOW ,
                            background = GREEN  ,
                            bold       = True   ,
                            blink      = False  ,
                            underline  = False  ,
                            fg_bright  = True   ,
                            bg_bright  = False  )

# ==============================================================================
## just for information 
def infostr ( what ) :
    """Just for information"""
    return colored_string ( what                ,
                            foreground = WHITE  ,
                            background = BLUE   ,
                            bold       = True   ,
                            blink      = False  ,
                            underline  = False  ,
                            bg_bright  = False  )
 
# =============================================================================
if __name__ == '__main__' :


    from ostap.logger.logger import getLogger 
    logger = getLogger ( 'ostap.logger.colorized')
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    fmt = 'colored_string(fg=%d,bg=%d,bold=%s,blink=%s,underline=%s,fg_bright=%s,bg_bright=%s)'
    
    for bg in range ( 8 ) :
        for fg in range ( 8 ) :
            for fg_bright  in ( True , False ) :
                for bg_bright  in ( True , False ) :
                    if  ( fg , fg_bright ) == ( bg , bg_bright ) : continue                    
                    for bold in ( True , False ) :
                        for blink in ( True , False ) :
                            for underline  in ( True , False ) :
                                pars = ( fg        ,
                                         bg        ,
                                         bold      ,
                                         blink     ,
                                         underline ,
                                         fg_bright ,
                                         bg_bright )
                                logger.info ( colored_string ( fmt % pars , *pars ) )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
