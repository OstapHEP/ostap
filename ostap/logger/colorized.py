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
    'attstr'           , ## just for information/attenttion
    ##
    'critical_info'    , ## `critical`    information 
    'error_info'       , ## `error`       information 
    'warning_info'     , ## `warning`     information
    'attention_info'   , ## `attention`   information
    'attention_info'   , ## `attention`   information
    'info_info'        , ## `information` information
    ##
    'with_colors'      , ## Is colorization enabled?
    'set_with_colors'  , ## Enable/Disable colorization
    ##
    'markup'           , ## some primitive markup  
    )
# =============================================================================
import os, sys, re 
# =============================================================================
# - is sys.stdout attached to terminal or not ?
# from ostap.utils.basic import isatty
def isatty ( stream = None ) :
    """ Is the stream is attached to terminal?
    >>> stream = ...
    >>> if isatty( stream ) : print('Teminal!')
    >>> if isatty() : print('stdout is terminal!')
    """
    if not stream : stream = sys.stdout
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        return stream.isatty()
        # =====================================================================
    except : pass # ===========================================================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        return os.isatty ( stream.fileno() ) 
        # =====================================================================
    except : pass # ===========================================================
    # =========================================================================
    #
    return False
# =============================================================================
## global flag to indicate if we use colors 
__with_colors__ = isatty ()   
# =============================================================================
## Is colorization enabled? 
def with_colors() :
    """ Is colorization enabled? """
    global __with_colors__
    return bool ( __with_colors__ ) and isatty() 
# =============================================================================
## Enable/disable colorization
def set_with_colors ( use ) :
    """ Enable/disable colorization"""
    global __with_colors__
    __with_colors__ = bool ( use )
    return with_colors () 

# =============================================================================
## BASIC ASCII colors :
BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = list ( range ( 8 ) )
# =============================================================================
COLOR_FMT      = "\033[%sm"  ## general format for color 

RESET_SEQ      = COLOR_FMT % 0 
BOLD_SEQ       = COLOR_FMT % 1
LIGHT_SEQ      = COLOR_FMT % 2
ITALIC_SEQ     = COLOR_FMT % 3
UNDERLINE_SEQ  = COLOR_FMT % 4
BLINK_FAST_SEQ = COLOR_FMT % 5
BLINK_SLOW_SEQ = COLOR_FMT % 6 
REVERSE_SEQ    = COLOR_FMT % 7
HIDE_SEQ       = COLOR_FMT % 8
CROSS_SEQ      = COLOR_FMT % 9
BLINK_SEQ      = BLINK_FAST_SEQ

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
    """ Attention string: yellow on red """
    return colored_string ( what                ,
                            foreground = YELLOW ,
                            background = RED    ,
                            bold       = True   ,
                            blink      = True   ,
                            underline  = True   ,
                            fg_bright  = True   ,
                            bg_bright  = True   )

# =====-=======================================================================
## attention!
def attstr ( what ) :
    """ Attention string: white on blue"""
    return colored_string ( what                ,
                            foreground = WHITE  ,
                            background = BLUE   ,
                            bold       = True   ,
                            blink      = True   ,
                            underline  = False  ,
                            fg_bright  = False  ,
                            bg_bright  = False  )

# =============================================================================
## allright 
def allright ( what ) :
    """ Allright string """
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
    """ Just for information"""
    return colored_string ( what                ,
                            foreground = WHITE  ,
                            background = BLUE   ,
                            bold       = True   ,
                            blink      = False  ,
                            underline  = False  ,
                            bg_bright  = False  ,
                            fg_bright  = False  )

# =============================================================================
## Critical information
def critical_info ( what ) :
    """ Critical information """ 
    return colored_string ( what                ,
                            foreground = RED    ,
                            background = BLUE   ,
                            bold       = True   ,
                            blink      = True   ,
                            underline  = True   ,
                            bg_bright  = True   ,
                            fg_bright  = True   )

# =============================================================================
## Error information
def error_info ( what ) :
    """ Error information """ 
    return colored_string ( what                ,
                            foreground = YELLOW ,
                            background = RED    ,                            
                            blink      = True   ,
                            bold       = True   ,
                            underline  = False  , 
                            fg_bright  = True   ,
                            bg_bright  = True   )

# ==============================================================================
## Critical information
def warning_info ( what ) :
    """ Warning information """ 
    return colored_string ( what                ,
                            foreground = RED    ,
                            background = YELLOW ,
                            bold       = False  ,
                            blink      = False  ,
                            underline  = True   ,
                            bg_bright  = True   ,
                            fg_bright  = True   )

# ==============================================================================
## Attention information
def attention_info  ( what ) :
    """ Attention information """ 
    return colored_string ( what                ,
                            foreground = WHITE  ,
                            background = BLUE   ,
                            blink      = True   ,
                            bold       = True   , 
                            underline  = False  ,
                            bg_bright  = True   ,
                            fg_bright  = True   )

# ==============================================================================
## Info information
def info_info ( what ) :
    """ Info information """ 
    return colored_string ( what                ,
                            foreground = WHITE  ,
                            background = BLUE   ,
                            blink      = False  ,
                            bold       = False  ,
                            underline  = False  , 
                            fg_bright  = False  ,
                            bg_bright  = False  )


# =============================================================================
## Apply some oversimplified markup tranformation (a'la github)
_markup_ = (
    ## undeline      :           __ text __   
    ( re.compile ( r'(?P<UL>__.+?__)'         ) , 2 , '%s%%s%s'   % ( UNDERLINE_SEQ ,              RESET_SEQ ) ) ,
    ## bold & blink  :          *** text *** 
    ( re.compile ( r'(?P<UL>\*\*\*.+?\*\*\*)' ) , 3 , '%s%s%%s%s' % ( BOLD_SEQ      , BLINK_SEQ   , RESET_SEQ ) ) , 
    ## bold & italic :           ** text ** 
    ( re.compile ( r'(?P<UL>\*\*.+?\*\*)'     ) , 2 , '%s%s%%s%s' % ( BOLD_SEQ      , ITALIC_SEQ  , RESET_SEQ ) ) ,
    ## bold          :            * text * 
    ( re.compile ( r'(?P<UL>\*.+?\*)'         ) , 1 , '%s%%s%s'   % ( BOLD_SEQ      ,               RESET_SEQ ) ) ,
    ## italic        :            _ text _ 
    ( re.compile ( r'(?P<UL>~.+?~)'           ) , 1 , '%s%%s%s'   % ( ITALIC_SEQ    ,               RESET_SEQ ) ) ,
    ## blink         :            > text < 
    ( re.compile ( r'(?P<UL>>.+?<)'           ) , 1 , '%s%%s%s'   % ( BLINK_SEQ      ,              RESET_SEQ ) ) ,
    ## cross-out/sstrike-out :   -- text __ 
    ( re.compile ( r'(?P<UL>--.+?--)'         ) , 2 , '%s%%s%s'   % ( CROSS_SEQ      ,              RESET_SEQ ) ) , 
    ## reverse       :           < text >
    ( re.compile ( r'(?P<UL>\^.+?\^)'         ) , 1 , '%s%%s%s'   % ( REVERSE_SEQ    ,              RESET_SEQ ) ) , 
)

# ================================================================================
## Very simeel markup-like tranformations
#  - <code>  __ text __  </code> underline 
#  - <code> *** text *** </code> bold and blink
#  - <code>  ** text **  </code> bold and italic
#  - <code>   * text *   </code> bold
#  - <code>   ~ text ~   </code> italic
#  - <code>   ^ text ^   </code> reverse 
#  - <code>   > text <   </code> blink
#  - <code>  -- text --  </code> cross-out/strike-out
def markup ( what ) :
    """ Very simple markup-like tranformations
    - `  __ text __  ` underline 
    - ` *** text *** ` bold and blink
    - `  ** text **  ` bold and italic
    - `   * text *   ` bold
    - `   ~ text ~   ` italic
    - `   > text <   ` blink
    - `   ^ text ^   ` reverse 
    - `  -- text --  ` cross-out/strike-out 
    """
    ## nothing to colorize
    if not what : return what
    
    ## if colorizing is isabled
    enabled = with_colors()
    
    for c , l , fmt in _markup_ :
        qq = tuple ( s for s in re.findall ( c , what ) )
        for s in qq :
            q    = s [ l : -l ]
            ## transform ( or just remove indicators)
            r    = fmt % q if enabled else q 
            what = what.replace ( s , r )
        
    return what


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

    ## markup ?
    logger.info ( markup ( 'Markup: ***BOLD&BLINK***' ) )
    logger.info ( markup ( 'Markup: **BOLD&ITALIC**'  ) )
    logger.info ( markup ( 'Markup: *BOLD*'           ) )
    logger.info ( markup ( 'Markup: >BLINK<'          ) )
    logger.info ( markup ( 'Markup: ~ITALIC~'         ) )
    logger.info ( markup ( 'Markup: __UNDERLINE__'    ) )
    logger.info ( markup ( 'Markup: ^REVERSE^'        ) )
                     
# =============================================================================
##                                                                      The END 
# =============================================================================
