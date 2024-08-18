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
    'pretty_float'       , ## pretty print of the floating number 
    'pretty_err1'        , ## pretty print of the floating number with error
    'pretty_err2'        , ## pretty print of the floating number with asymmetric error
    'pretty_errs'        , ## pretty print of the asymmetric error
    ## 
    'fmt_pretty_float'   , ## format for pretty print of the floatig number 
    'fmt_pretty_err1'    , ## format for pretty print of the value with error
    'fmt_pretty_err2'    , ## format for pretty print of the value with asymmetric error
    'fmt_pretty_errs'    , ## format for pretty print of the value with multiple errors
    ""
    'add_expo'           , ## add an exponential factor to sting representaiton of the object 
    ##
    'multicolumn'        , ## format the list of strings into multicolumn block
    'DisplayTree'        , ## display tree-like structures 
)
# =============================================================================
from   ostap.core.ostap_types import integer_types, num_types  
from   ostap.logger.logger    import logVerbose,  logDebug, logInfo, logWarning, logError
from   ostap.math.base        import isfinite, iszero, frexp10
from   ostap.math.ve          import VE
from   ostap.logger.mute      import ( mute   , mute_py ,
                                       tee_py , tee_cpp ,
                                       output , silence , silence_py ,
                                       MuteC  , MutePy  ,
                                       TeeCpp , TeePy   , OutputC    )
from   ostap.utils.basic      import NoContext
import time, os, sys, math  ## attention here!!
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger    import getLogger, logColor, logNoColor 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.logger.utils' )
else                       : logger = getLogger( __name__ )
del getLogger 
# =============================================================================
## Formats for nice printout of the object with errors   ( string + exponent)
#  @code
#  fmtv , fmte , expo = fmt_pretty_errs ( number , ( e1 , e2 , e3 ) ) 
#  @endcode
#  @return formats for nice string and the separate exponent 
def fmt_pretty_errs ( value              ,
                      errors      = ()   , 
                      width       = 6    ,
                      precision   = 4    ) : 
    """ Formats for nice printout of the object with errors  ( strings + exponent)
    >>> fmtv , fmte , expo = fmt_pretty_errs ( number , ( e1 , e2 , e3 )  ) 
    """
    
    assert isinstance ( width  , integer_types ) and \
        isinstance ( precision , integer_types ) and 2 <= precision < width, \
        "Invalid width/precision parameters: %s/%s" % ( width , precision ) 
    
    v = value
    e = max ( abs ( e ) for e in errors ) if errors else abs ( v ) 
    
    ## quantity that actually defines the format 
    av    = max ( abs ( v ) , e )

    if   100 <= av < 1000 :
        
        fmtv  = '%%+%d.%df' %  ( width , precision - 2 )
        fmte  = '%%-%d.%df' %  ( width , precision - 2 )        
        return fmtv, fmte , 0
    
    elif 10  <= av < 100  :
        
        fmtv  = '%%+%d.%df' %  ( width , precision - 1 )
        fmte  = '%%-%d.%df' %  ( width , precision - 1 )
        return fmtv  , fmte , 0
    
    elif 0.1 <= av < 10 :
        
        fmtv  = '%%+%d.%df' %  ( width , precision     )
        fmte  = '%%-%d.%df' %  ( width , precision     )
        return fmtv , fmte , 0

    if iszero ( av ) :
        fmtv  = '%%+%d.%df' %  ( width , precision     )
        fmte  = '%%-%d.%df' %  ( width , precision     )
        return fmtv , fmte , 0 

    #
    ## here we scale input data and try to get formats for scaled data
    #
    
    v_a , v_e = frexp10 ( av )
    v_ee  = v_e - 1
    n , r  = divmod  ( v_ee , 3 )    

    scale  = 10** ( r - v_ee  )    
    v     *= scale
    errs   = [ e*scale for e in errors ]  

    #
    ## get formats for scaled data
    # 
    fmtv , fmte , expo = fmt_pretty_errs ( value     = v         ,
                                           errors    = errs      ,
                                           width     = width     ,
                                           precision = precision )
    
    return fmtv , fmte , expo + 3 * n 

# =============================================================================
## Format for nice printout of the floating number (string + exponent)
#  @code
#  fmt , expo = fmt_pretty_float ( number ) 
#  @endcode
#  @return format for nice string and the separate exponent 
def fmt_pretty_float ( value         ,
                       width     = 6 ,
                       precision = 4 ) :
    """Format for nice printout of the floating number
    - return format for nice string and the separate exponent 
    >>> fmt , expo = fmt_pretty_float ( number ) 
    """

    assert isinstance ( value     , num_types     ),\
        "Invalid `value' parameter: %s"   % type ( value ) 
    
    value = float ( value )

    ## finite?
    if not isfinite ( value ) : return "%s" , 0

    ## get the format 
    fmtv , _ , expo = fmt_pretty_errs ( value , width = width , precision = precision ) 

    return fmtv, expo

# ===============================================================================
## Formats nice printout of the valuw with error
#  @code
#  fmt, fmtv , fmte  , expo = fmt_pretty_err1 ( value , error  ) 
#  @endcode
#  @return formats nice string and the separate exponent 
def fmt_pretty_err1 ( value              ,
                      error              , 
                      width       = 6    ,
                      precision   = 4    ,
                      parentheses = True ) :
    """Formats for nice printout of the vaoleu with error 
    - return formats for nice stirng and the separate exponent 
    >>> fmt , fmt_v , fmt_e , n = pretty_err1 ( value , error  ) 
    """    
    assert isinstance ( value     , num_types     ),\
        "Invalid `value' parameter: %s"   % type ( value ) 
    assert isinstance ( error     , num_types     ),\
        "Invalid `error' parameter: %s"   % type ( error ) 
    
    v = float ( value )
    e = max   ( 0.0 , float ( error  ) )
    
    ## get the formats 
    fmtv , fmte , expo = fmt_pretty_errs ( value     = v         ,
                                           errors    = ( e , )   ,
                                           width     = width     ,
                                           precision = precision ) 
    
    fmt = '%s +/- %s' % ( fmtv , fmte ) 
    if parentheses : fmt = '( ' + fmt + ' )'
    
    return fmt, fmtv, fmte , expo

# ===============================================================================
## Formats for nice printout of the object with asymmetric  errors ( string + exponent)
#  @code
#  fmt , fmtv , fmte , expo = fmt_pretty_err2 ( number , errlow  , errhigh ) 
#  @endcode
#  @return formats for nice string and the separate exponent 
def fmt_pretty_err2 ( value              ,
                      errlow             ,
                      errhigh            ,
                      width       = 6    ,
                      precision   = 4    ,
                      parentheses = True ) :
    """ Formats for nice printout of the object with asymmetric  errors ( string + exponent)
    >>> fmt , fmtv , fmte , expo = fmt_pretty_err2 ( number , elow  , ehigh ) 
    """ 
    assert isinstance ( value     , num_types     ),\
        "Invalid `value'  parameter: %s"   % type ( value  ) 
    assert isinstance ( errlow    , num_types     ),\
        "Invalid `errlow' parameter: %s"   % type ( errlow  ) 
    assert isinstance ( errhigh   , num_types     ),\
        "Invalid `errigh' parameter: %s"   % type ( errhigh ) 
   
    v  = float ( value   )
    eh = float ( errhigh )
    el = float ( errlow  )

    ## get the formats 
    fmtv , fmte , expo = fmt_pretty_errs ( value     = v           ,
                                           errors    = ( eh , el ) ,
                                           width     = width       ,
                                           precision = precision   ) 

    fmt = '%s -/%s +/%s' % ( fmtv , fmte , fmte )
    if parentheses : fmt = '( ' + fmt + ' )'
    
    return fmt, fmtv, fmte , expo


# =============================================================================
## Nice printout of the floating number (string + exponent)
#  @code
#  s , expo = pretty_float ( number ) 
#  @endcode
#  @return nice stirng and the separate exponent 
def pretty_float ( value         ,
                   width     = 6 ,
                   precision = 4 ) :
    """Nice printout of the floating number
    - return nice string and the separate exponent 
    >>> s , expo = pretty_float ( number )
    """
    assert isinstance ( value     , num_types     ),\
        "Invalid `value'  parameter: %s"   % type ( value  ) 

    v = float ( value )
    
    if   math.isinf ( v ) :
        return '+inf' if 0 < v else '-inf' , 0    
    elif math.isnan ( v ) :
        return 'NaN' , 0

    fmt , expo = fmt_pretty_float ( value     = v         ,
                                    width     = width     ,
                                    precision = precision )
    if expo :
        scale = 10 ** expo
        v /= scale
        
    return  fmt % v , expo 

# ===============================================================================
## nice printout of value with error  ( string + exponent)
#  @code
#  s , expo = pretty_err1 ( number , error ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_err1 ( value              ,
                  error              , 
                  width       = 6    ,
                  precision   = 4    ,
                  parentheses = True ) :
    """Nice printout of the value with error ( string + exponent)
    - return nice stirng and the separate exponent 
    >>> s , expo = pretty_err1 ( value , error ) 
    """
    assert isinstance ( value     , num_types     ),\
        "Invalid `value'  parameter: %s"   % type ( value  ) 
    assert isinstance ( error     , num_types     ),\
        "Invalid `error'  parameter: %s"   % type ( error  ) 
    
    v =             float ( value ) 
    e = max ( 0.0 , float ( error ) ) 
    
    finv = isfinite ( v )
    fine = isfinite ( e )
    
    if not finv and not fine  :
        fmt = '%%+%ds +/- %%-%ds' % ( width , width ) 
        if parentheses : fmt = '( ' + fmt + ' )'
        return fmt % ( v , e ) , 0 
    elif not finv  :
        fe , ne = pretty_float ( value = e , width = width , precision = precision ) 
        fmt = '%%+%ds +/- %%-%ds' % ( width , width ) 
        if parentheses : fmt = '( ' + fmt + ' )'
        return fmt % ( v , fe ) , ne 
    elif not fine  :
        fv , nv = pretty_float ( value = v , width = width , precision = precision ) 
        fmt = '%%%ds +/- %%-%ds' % ( width , width ) 
        if parentheses : fmt = '( ' + fmt + ' )'
        return fmt % ( fv , e ) , nv
    
    fmt , _ , _ , expo = fmt_pretty_err1  ( value       = v           ,
                                            error       = e           ,
                                            width       = width       ,
                                            precision   = precision   ,
                                            parentheses = parentheses )
    values = v , e
    if expo :
        scale = 10 ** expo 
        values = tuple ( v / scale for v in values ) 
    
    return fmt % values , expo 

# ===============================================================================
## nice printout of the object with asymmetric  errors  ( string + exponent)
#  @code
#  s , expo = pretty_err2 ( value , errlow , errhigh ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_err2 ( value              ,
                  errlow             ,
                  errhigh            ,
                  width       = 6    ,
                  precision   = 4    ,
                  parentheses = True ) :
    """Nice printout of the object with asymmetric  errors   ( string + exponent)
    >>> s , expo = pretty_err2 ( number , elow , ehigh ) 
    """

    assert isinstance ( value     , num_types     ),\
        "Invalid `value'  parameter: %s"   % type ( value  ) 
    assert isinstance ( errlow    , num_types     ),\
        "Invalid `errlow' parameter: %s"   % type ( errlow  ) 
    assert isinstance ( errhigh   , num_types     ),\
        "Invalid `errigh' parameter: %s"   % type ( errhigh ) 

    v     = float ( value   )
    ehigh = float ( errhigh )
    elow  = float ( errlow  )

    assert 0 <= elow or 0 <= ehigh, 'Both errlow/errhigh cannot be negative!'

    if   elow  <= 0          <= ehigh : pass
    elif ehigh <= 0          <= elow  : elow, ehigh = ehigh, elow 
    elif 0     <= elow and 0 <= ehigh : elow = -elow

    ## get the format 
    fmt , _ , _ , expo = fmt_pretty_err2 ( value       = v           ,
                                           errlow      = elow        ,
                                           errhigh     = ehigh       ,
                                           width       = width       ,
                                           precision   = precision   ,
                                           parenthesis = parentheses )

    values = value , abs ( elow ) , ehigh    
    if expo :
        scale = 10 ** expo 
        values = tuple ( v / scale for v in values ) 
    
    return fmt % values , expo 


# ===============================================================================
## nice printout of the asymmetric  errors  ( string + exponent)
#  @code
#  s , expo = pretty_errs ( errlow , errhigh ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_errs ( errlow             ,
                  errhigh            ,
                  width       = 6    ,
                  precision   = 4    ,
                  parentheses = True ) :
    """Nice printout of asymmetric  errors   ( string + exponent)
    >>> s , expo = pretty_errs ( errlow , errhigh ) 
    """
    assert isinstance ( errlow    , num_types     ),\
        "Invalid `errlow' parameter: %s"   % type ( errlow  ) 
    assert isinstance ( errhigh   , num_types     ),\
        "Invalid `errigh' parameter: %s"   % type ( errhigh ) 

    ehigh = float ( errhigh )
    elow  = float ( errlow  )

    assert 0 <= elow or 0 <= ehigh, 'Both errlow/errhigh cannot be negative!'

    if   elow  <= 0          <= ehigh : pass
    elif ehigh <= 0          <= elow  : elow, ehigh = ehigh, elow 
    elif 0     <= elow and 0 <= ehigh : elow = -elow

    ## get the format 
    _ , fmte , expo = fmt_pretty_errs ( value     = 0.0              ,
                                        errors    = ( elow , ehigh ) ,
                                        width     = width            ,
                                        precision = precision        )

    ## construct format 
    fmt = ' -/%s +/%s' % ( fmte , fmte )
    if parentheses : fmt = '( ' + fmt + ' )'

    values = value , abs ( elow ) , ehigh    
    if expo :
        scale = 10 ** expo 
        values = tuple ( v / scale for v in values ) 
    
    return fmt % values , expo 

# =============================================================================
##  Add the exponent to the sting representation of th eobejuct
#   @code
#   value = ...
#   value , expo = pretty_XXX ( value , .... )
#   result = add_expo ( value , expo ) 
#   @endcode 
def add_expo ( value , expo , fmt = '%s x 10^%+d' ) :
    """ Add the exponent to the sting representation of th eobejuct
    >>> value = ...
    >>> value , expo = pretty_XXX ( value , .... )
    >>> result = add_expo ( value , expo ) 
    """
    assert isinstance ( expo , integer_types ) ,  "Invalid type of `expo':%s" % type ( expo )
    ## 
    if not expo : return value
    ## 
    return fmt % ( value , expo )
    
# =============================================================================
## format list of strings into multicolumn string
#  @code 
#  >>> strings =  ....
#  >>> table   = multicolumn (  strings , indent = 2 )
#  >>> print table
#  @endcode 
def multicolumn ( lines , term_width = None , indent = 0 , pad = 2 ):
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

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    vv = 1./3.0 
    for i in range ( -10 , 20 ) :
        v = vv * ( 10**i )
        print ('VALUE' ,  v , pretty_float ( v ) )
        print ('FMTS'  ,  fmt_pretty_errs  (  v , ( v, v ) ) ) 
        
        
    
# =============================================================================
##                                                                     The END 
# =============================================================================
