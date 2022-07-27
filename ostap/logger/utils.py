#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some simple but useful utilities
#   - suppression of stdout/stderr 
#   - dumpting of stdout/stderr into file 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
#  
# =============================================================================
"""Module with some simple but useful utilities"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    #
    'tee_py'             , ## tee for Python's printouts
    'tee_cpp'            , ## tee for C++'s    printouts
    'output'             , ## redirect stdout/stderr into the file 
    'mute_py'            , ## suppress stdout/strerr Python printout 
    'silence_py'         , ## ditto 
    'mute'               , ## context manager to suppress stdout/strerr printout 
    'silence'            , ## ditto 
    'NoContext'          , ## empty context manager
    ## logging   
    'logColor'           , ## switch on  locally the colored logging
    'logNoColor'         , ## switch off locally the colored logging
    'logVerbose'         , ## redefine (locally) the logging level
    'logDebug'           , ## redefine (locally) the logging level
    'logInfo'            , ## redefine (locally) the logging level
    'logWarning'         , ## redefine (locally) the logging level
    'logError'           , ## redefine (locally) the logging level
    ##
    'pretty_float'       , ## pretty print of the floatig number 
    'pretty_ve'          , ## pretty print of the value with error 
    'pretty_2ve'         , ## pretty print of the value with asymmetric errors
    ##
    'fmt_pretty_float'   , ## format for pretty print of the floatig number 
    'fmt_pretty_ve'      , ## format for pretty print of the value with error 
    'fmt_pretty_2ve'     , ## format for pretty print of the value with asymmetric errors 
    ##
    'multicolumn'        , ## format the list of strings into multicolumn block
    'DisplayTree'        , ## display tree-like structures 
    )
# =============================================================================
import time, os, sys, math  ## attention here!!
from   ostap.logger.logger    import logVerbose,  logDebug, logInfo, logWarning, logError
from   ostap.logger.mute      import ( mute   , mute_py ,
                                       tee_py , tee_cpp ,
                                       output , silence , silence_py ,
                                       MuteC  , MutePy  ,
                                       TeeCpp , TeePy   , OutputC    )
from   ostap.utils.basic      import NoContext
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger    import getLogger, logColor, logNoColor 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.logger.utils' )
else                       : logger = getLogger( __name__ )
del getLogger 
# =============================================================================

# =============================================================================
## Format for nice printout of the floating number (string + exponent)
#  @code
#  fmt , n = fmt_pretty_float ( number ) 
#  @endcode
#  @return format for nice string and the separate exponent 
def fmt_pretty_float ( value , width = 8 , precision = 6 ) :
    """Format for nice printout of the floating number
    - return format for nice string and the separate exponent 
    >>> fmt , n = fmt_pretty_float ( number ) 
    """
    from   ostap.core.ostap_types import integer_types, num_types  

    assert isinstance ( value     , num_types     ),\
           'Invalid value parameter %s/%s'   % ( value , type ( value ) )

    ## not finite?
    from ostap.math.base import isfinite 
    if not isfinite ( value ) : return "%s" , 0
    
    assert isinstance ( width     , integer_types ) and \
           isinstance ( precision , integer_types ) and 2 <= precision < width, \
           "Invalid width/precision parameters %s/%s" % ( width , precision )
    
    v  = value 
    av = abs ( v ) 
    if   100 <= av < 1000 :
        fmt2 = '%%+%d.%df' % ( width , precision - 2 )
        return fmt2 , 0 
    elif 10  <= av < 100  :
        fmt1 = '%%+%d.%df' % ( width , precision - 1 )
        return fmt1 , 0  
    elif 1   <= av < 10   :
        fmt0 = '%%+%d.%df' % ( width , precision     )
        return fmt0 , 0  
    elif 0.1 <= av < 1    :
        fmt0 = '%%+%d.%df' % ( width , precision     )
        return fmt0 , 0  
    
    from  ostap.math.base        import frexp10, iszero    
    if iszero ( av ) :
        fmt0 = '%%+%d.%df' % ( width , precision )
        return fmt0 , 0  
    
    v_a , v_e = frexp10 ( v )
    if 0 == v_e and iszero ( v_a  ) :
        fmt0 = '%%+0%d.%df' % ( width , precision )
        return fmt0 , 0 
            
    v_a *= 10 
    v_e -=  1 

    n , r = divmod  ( v_e , 3 )
        
    ra = v_a * ( 10**r )

    fmt , p = fmt_pretty_float ( ra , width , precision )

    return fmt , p + 3 * n


# ===============================================================================
## Formats nice printout of the ValueWithError object  ( string + exponent)
#  @code
#  fmt, fmt_v , fmt_e  , n = fmt_pretty_ve ( number ) 
#  @endcode
#  @return formats nice string and the separate exponent 
def fmt_pretty_ve ( value , width = 8 , precision = 6 , parentheses = True ) :
    """Formats for nice printout of the ValueWithError object  ( string + exponent)
    - return formats for nice stirng and the separate exponent 
    >>> fmt , fmt_v , fmt_e , n = pretty_ve ( number ) 
    """
    
    from ostap.math.ve          import VE
    from ostap.core.ostap_types import integer_types 
    
    v =           value.value ()
    e = max ( 0 , value.error () )

    assert isinstance ( value , VE ), 'Invalid type for value %s' % value 

    assert isinstance ( width     , integer_types ) and \
           isinstance ( precision , integer_types ) and 2 <= precision < width, \
           "Invalid width/precision parameters %s/%s" % ( width , precision ) 
    
    av = max ( abs ( v ) ,  e )

    if   100 <= av < 1000 :
        fmtv = '%%+%d.%df' % ( width , precision - 2 )
        fmte = '%%-%d.%df' % ( width , precision - 2 )
        fmt  = '%%+%d.%df +/- %%-%d.%df' % ( width , precision - 2 , width , precision - 2 )
        if parentheses : fmt = '( ' + fmt + ' )'
        return fmt , fmtv , fmte , 0 
    elif 10  <= av < 100  :
        fmtv = '%%+%d.%df' % ( width , precision - 1  )   
        fmte = '%%-%d.%df' % ( width , precision - 1 )   
        fmt  = '%%+%d.%df +/- %%-%d.%df' % ( width , precision - 1 , width , precision - 1 )   
        if parentheses : fmt = '( ' + fmt + ' )'
        return fmt , fmtv , fmte , 0 
    elif 1   <= av < 10   :
        fmtv = '%%+%d.%df' % ( width , precision     )
        fmte = '%%-%d.%df' % ( width , precision     )
        fmt  = '%%+%d.%df +/- %%-%d.%df' % ( width , precision     , width , precision     )
        if parentheses : fmt = '( ' + fmt + ' )'        
        return fmt , fmtv , fmte , 0  
    elif 0.1 <= av < 1     :
        fmtv = '%%+%d.%df' % ( width , precision     )
        fmte = '%%-%d.%df' % ( width , precision     )
        fmt  = '%%+%d.%df +/- %%-%d.%df' % ( width , precision     , width , precision     )
        if parentheses : fmt = '( ' + fmt + ' )'
        return fmt , fmtv , fmte , 0  


    from  ostap.math.base        import frexp10, iszero 
    if iszero ( av ) :
        fmtv = '%%+%d.%df' % ( width , precision     )
        fmte = '%%-%d.%df' % ( width , precision     )
        fmt  = '%%+%d.%df +/- %%-%d.%df' % ( width , precision     , width , precision     )
        if parentheses : fmt = '( ' + fmt + ' )'        
        return fmt , fmtv , fmte , 0  

    v_a , v_e = frexp10 ( av )

    v   /= 10**(v_e-1)  
    e   /= 10**(v_e-1)
    v_e -= 1 
    
    n , r = divmod  ( v_e , 3 )    

    v *= 10**r
    e *= 10**r 

    fmt , fmtv , fmte  , p = fmt_pretty_ve ( VE ( v , e * e ) , width, precision , parentheses )    
    return fmt , fmtv , fmte , p + 3 * n 


# ===============================================================================
## Formats for nice printout of the object with asymmetric  errors   ( string + exponent)
#  @code
#  fmt , fmtv , fmte , n = fmt_pretty_2ve ( number , ehigh , elow ) 
#  @endcode
#  @return formats for nice string and the separate exponent 
def fmt_pretty_2ve ( value              ,
                     eh                 ,
                     el                 ,
                     width       = 8    ,
                     precision   = 6    ,
                     parentheses = True ) :
    
    from ostap.math.ve          import VE
    from ostap.core.ostap_types import integer_types, num_types  

    assert isinstance ( value     , num_types     ),\
           'Invalid value parameter %s/%s'   % ( value , type ( value ) )      
    assert isinstance ( eh       , num_types     ),\
           'Invalid eh    parameter %s/%s'   % ( eh    , type ( eh    ) )      
    assert isinstance ( el       , num_types     ),\
           'Invalid el    parameter %s/%s'   % ( el    , type ( el    ) )      

    v = value 
    e = max ( abs ( eh ), abs ( el ) )

    assert isinstance ( width     , integer_types ) and \
           isinstance ( precision , integer_types ) and 2 <= precision < width, \
           "Invalid width/precision parameters %s/%s" % ( width , precision ) 

    assert 0 <= eh or 0 <= el, 'Both errors cannot be negative!'
    
    if   eh < 0  and el < 0 :
        eh , el =   el , eh
    elif eh >= 0 and el < 0 :
        eh , el = eh , abs ( el )
    
    av = max ( abs ( v ) , e )

    if   100 <= av < 1000 :
        fmtv  = '%%+%d.%df' %  ( width , precision - 2 )
        fmte  = '%%-%d.%df' %  ( width , precision - 2 )
        fmt   = '%%+%d.%df +/%%-%d.%df -/%%-%d.%df' %  ( width , precision - 2 , width , precision - 2 , width , precision - 2 )
        if parentheses : fmt = '( ' + fmt + ' )'        
        
        return fmt , fmtv , fmte , 0  
    elif 10  <= av < 100  :        
        fmtv  = '%%+%d.%df' %  ( width , precision - 1 )
        fmte  = '%%-%d.%df' %  ( width , precision - 1 )
        fmt   = '%%+%d.%df +/%%-%d.%df -/%%-%d.%df' %  ( width , precision - 1 , width , precision - 1 , width , precision - 1 )
        if parentheses : fmt = '( ' + fmt + ' )'        
        return fmt , fmtv  , fmte , 0
    elif 1   <= av < 10   :
        fmtv  = '%%+%d.%df' %  ( width , precision     )
        fmte  = '%%-%d.%df' %  ( width , precision     )
        fmt   = '%%+%d.%df +/%%-0%d.%df -/%%-%d.%df' %  ( width , precision     , width , precision     , width , precision     )
        if parentheses : fmt = '( ' + fmt + ' )'        
        return fmt , fmtv , fmte , 0 
    elif 0.1 <= av < 1    :
        fmtv  = '%%+%d.%df' %  ( width , precision     )
        fmte  = '%%-%d.%df' %  ( width , precision     )
        fmt   = '%%+%d.%df +/%%-0%d.%df -/%%-%d.%df' %  ( width , precision     , width , precision     , width , precision     )
        if parentheses : fmt = '( ' + fmt + ' )'        
        return fmt , fmtv , fmte , 0 

    from  ostap.math.base        import frexp10, iszero
    if iszero ( av ) :
        fmtv  = '%%+%d.%df' %  ( width , precision     )
        fmte  = '%%-%d.%df' %  ( width , precision     )
        fmt   = '%%+%d.%df +/%%-%d.%df -/%%-%d.%df' %  ( width , precision     , width , precision     , width , precision     )
        if parentheses : fmt = '( ' + fmt + ' )'        
        return fmt , fmtv , fmte , 0 
                
    v_a , v_e = frexp10 ( av )

    v   /= 10**(v_e-1)
    eh  /= 10**(v_e-1)
    el  /= 10**(v_e-1)
    v_e -= 1

    n , r = divmod  ( v_e , 3 )    

    v  *= 10**r
    eh *= 10**r 
    el *= 10**r 

    fmt , fmtv  , fmte  , p = fmt_pretty_2ve ( v , eh , el , width , precision , parentheses )
    
    return fmt , fmtv , fmte  , p + 3 * n 


# =============================================================================
## Nice printout of the floating number (string + exponent)
#  @code
#  s , n = pretty_float ( number ) 
#  @endcode
#  @return nice stirng and the separate exponent 
def pretty_float ( value , width = 8 , precision = 6 ) :
    """Nice printout of the floating number
    - return nice string and the separate exponent 
    >>> s , n = pretty_float ( number ) 
    """
    fmt , n = fmt_pretty_float ( value , width , precision )
    return  fmt % ( value / 10**n ) , n 

# ===============================================================================
## nice printout of the ValueWithError object  ( string + exponent)
#  @code
#  s , n = pretty_ve ( number ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_ve ( value , width = 8 , precision = 6 , parentheses = True ) :
    """Nice printout of the ValueWithError object  ( string + exponent)
    - return nice stirng and the separate exponent 
    >>> s , n = pretty_ve ( number ) 
    """
    from ostap.math.ve          import VE
    value =  VE ( value )
    
    fmt , fmtv , fmte , n = fmt_pretty_ve ( value , width , precision , parentheses )
    
    v =           value.value ()   
    e = max ( 0 , value.error () ) 

    return fmt % ( v / 10**n , e / 10**n ) , n

# ===============================================================================
## nice printout of the object with asymmetric  errors   ( string + exponent)
#  @code
#  s , n = pretty_2ve ( number , ehigh , elow ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_2ve ( value              ,
                 eh                 ,
                 el                 ,
                 width       = 8    ,
                 precision   = 6    ,
                 parentheses = True ) :

    assert 0 <= eh or 0 <= el, 'Both errors cannot be negative!'
    
    if   eh < 0  and el < 0 :
        eh , el =   el , eh
    elif eh >= 0 and el < 0 :
        eh , el = eh , abs ( el )

    fmt , fmtv , fmte , n = fmt_pretty_2ve ( value , eh , el , width , precision , parentheses )
    
    return fmt  % ( value / 10**n , eh / 10**n , el / 10**n ) , n

# =============================================================================
## format list of strings into multicolumn string
#  @code 
#  >>> strings =  ....
#  >>> table   = multicolumn (  strings , indent = 2 )
#  >>> print table
#  @endcode 
def multicolumn ( lines , term_width=None , indent = 0 , pad = 2 ):
    """Format list of strings into multicolumn string
    >>> strings =  ....
    >>> table   = multicolumn (  strings , indent = 2 )
    >>> print table 
    """
    n_lines = len(lines)
    if n_lines == 0:
        return
    
    if not term_width :
        from ostap.utils.basic import terminal_size 
        h , term_width = terminal_size()
        
    col_width = max(len(line) for line in lines)
    n_cols = int((term_width + pad - indent)/(col_width + pad))
    n_cols = min(n_lines, max(1, n_cols))
    
    col_len = int(n_lines/n_cols) + (0 if n_lines % n_cols == 0 else 1)
    if (n_cols - 1) * col_len >= n_lines:
        n_cols -= 1
        
    cols = [lines[i*col_len : i*col_len + col_len] for i in range(n_cols)]
    
    rows        = list(zip(*cols))
    rows_missed = list(zip(*[col[len(rows):] for col in cols[:-1]]))
    rows.extend(rows_missed)

    result = []
    for row in rows:
        line = " "*indent + (" "*pad).join(line.ljust(col_width) for line in row)
        result.append ( line )
    return '\n'.join ( result )


# =======================================================================================
## @class DisplayTree
#  Very simple class to allow a rendering of Tree-like objects,
#  in particular the structure and the content of directories 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2020-01-029
class DisplayTree ( object ) :
    """Very simple class to allow a rendering of Tree-like objects,
    in particular the structure and the content of directories 
    """    
    prefix_item_middle   = '\033(0\x74\033(B' + '\033(0\x71\033(B' + '\033(0\x71\033(B'
    prefix_item_last     = '\033(0\x6d\033(B' + '\033(0\x71\033(B' + '\033(0\x71\033(B'
    prefix_parent_middle = '    '  
    prefix_parent_last   = '\033(0\x78\033(B' +  '   '


    def __init__  ( self            ,
                    payload         ,   ## the object 
                    parent  = None  ,   ## parent if any  
                    display = None  ,   ## display function: how to display it? 
                    isdir   = False ,   ## directory object ?
                    last    = False ) : ## last in the list? 
        
        self.__payload = payload
        self.__parent  = parent

        if   display :         ## use display function if specified 
            self.__display = display
        elif parent  :         ## otherwise pick if from the parent 
            self.__display = parent.__display
            
        self.__last    = last
        self.__isdir   = isdir        
        self.__depth   = self.parent.depth + 1 if parent else 0 

    @property
    def payload ( self ) :
        """``payload''  : the actual object, hold by the tree leaf"""
        return self.__payload
    @property
    def parent ( self ) :
        """``parent'' : the parent element in the tree"""
        return self.__parent
    @property
    def last   (  self ) :
        """``last''  : is it a alst object oin the list? """
        return self.__last
    @property
    def isdir  (  self ) :
        """``isdir''  : is it `directory-like' item? """
        return self.__isdir    
    @property
    def depth  (  self ) :
        """``depth''  : how deep is our tree? """
        return self.__depth 
    
    @property 
    def itemname ( self ) :
        """``itemname'' : how the name of the leaf is displayed?"""
        return self.__display ( self.payload , self.isdir )

    # ==========================================================================
    ##  show the (sub)tree  
    def showme ( self , prefix = '' ) :
        """Show the constructed (sub)tree  
        """        
        if self.parent is None:
            return self.itemname

        symbols = self.prefix_item_last if self.last else self.prefix_item_middle

        graph   = [ prefix + '{!s} {!s}'.format ( symbols , self.itemname ) ]

        parent = self.parent        
        while parent and parent.parent is not None:
            pp     = self.prefix_parent_middle if parent.last else self.prefix_parent_last
            graph.append  ( pp )
            parent = parent.parent

        return ''.join ( reversed ( graph ) )

    def __str__ (  self ) :
        return 'Leaf(%s,isdir=%s,last=%s)' % ( self.itemname ,
                                               self.isdir    ,
                                               self.last     ) 
    __repr__ = __str__
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap import banner
    logger.info ( __file__  + '\n' + banner )
    logger.info ( 80*'*'   )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 
    
# =============================================================================
# The END 
# =============================================================================
