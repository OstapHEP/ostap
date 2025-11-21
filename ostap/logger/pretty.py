#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file pretty.py
#  Helper functions for pretty prints  of some object
#   - float value 
#   - loat varules with errors 
#   - float values with asymmetric values
#   - float value swith multiple asymmetric values 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07 
# =============================================================================
""" Helper functions for pretty prints of some object
   - float value 
   - loat varules with errors 
   - float values with asymmetric values
   - float value swith multiple asymmetric values 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    ##
    ## Formats
    ##
    'fmt_pretty_values'  , ## (common) format for several floating values
    'fmt_pretty_errors'  , ## (common) format for values with uncertaintines 
    'fmt_pretty_float'   , ## (common) format for floaitng value
    'fmt_pretty_error'   , ## (common) format for floaitng value with error 
    'fmt_pretty_error2'  , ## (common) format for floaitng value with asymmetric errors
    ##
    ## Formatters
    ##
    'pretty_float'       , ## pretty representation for floating value
    'pretty_error'       , ## pretty represetnation for floating value with (symmetric) error 
    'pretty_error2'      , ## pretty representation for floating value with asymmetric errors
    'pretty_errors'      , ## pretty representation for asymmetric errors
    ## 
    'pretty_print'       , ## dispatcher for pretty prints 
    'nice_print'         , ## ready-to-usepretty prints 
    'add_expo'           , ## add an exponetaion factor for sting representaion of the object
    ## Force LaTeX format 
    'pretty_latex'       , ## dispatcher for pretty prints 
    'nice_latex'         , ## ready-to-usepretty prints 
)

# =============================================================================
from   ostap.core.ostap_types import integer_types, num_types, sequence_types
from   ostap.logger.symbols   import plus_minus, times
from   ostap.utils.basic      import typename 
import math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.logger.pretty' )
else                       : logger = getLogger( __name__                )
# =============================================================================
logger.debug ( "Helper functions for pretty prints of some object" )
# =============================================================================

# =============================================================================
## 1-formats 
# =============================================================================
## Formats for nice printout of the floating objects 
#  @code
#  fmtv , expo = fmt_pretty_values ( e1 , e2 , e3 )
#  @endcode
#  @return formats for nice string and the separate exponent 
def fmt_pretty_values ( *values             ,
                        width       = 6     ,
                        precision   = 4     ,
                        with_sign   = True  ,
                        latex       = False ) : 
    """ Formats for nice printout of the numerical objects ( strings + exponent)
    >>> fmtv , expo = fmt_pretty_values ( e1 , e2 , e3 ) 
    """
    if not values : values = (  0, )  
    assert isinstance ( width  , integer_types ) and \
        isinstance ( precision , integer_types ) and 2 <= precision < width, \
        "Invalid width/precision parameters: %s/%s" % ( width , precision ) 

    assert values and all ( isinstance ( v , num_types ) for v in values ) , "Invalid type of `values'"
    
    av = max ( abs ( float ( v )  ) for v in values )

    from ostap.math.base import iszero, frexp10
    the_format =  '%%+%d.%df' if with_sign else '%%%d.%df'
    if   100 <= av < 1000 : return the_format % ( width , precision - 2 ) , 0 
    elif 10  <= av < 100  : return the_format % ( width , precision - 1 ) , 0 
    elif 0.1 <= av < 10   : return the_format % ( width , precision     ) , 0 
    elif iszero ( av )    : return the_format % ( width , precision     ) , 0 

    ## here we scale input data and try to get formats for scaled data
    
    v_a , v_e = frexp10 ( av )
    v_ee      = v_e - 1
    n , r     = divmod  ( v_ee , 3 )    

    scale  = 10 ** ( r - v_ee  )
    
    scaled = av * scale 

    ## get formats for properly scaled data
    fmtv , expo = fmt_pretty_values ( scaled                ,
                                      width     = width     ,
                                      precision = precision ,
                                      with_sign = with_sign ,
                                      latex     = latex     )
    
    return fmtv , expo + 3 * n 

# ==================================================================================
## Formats for nice printout of the object with errors  ( string + exponent)
#  @code
#  fmtv , fmte , expo = fmt_pretty_errors ( number , ( e1 , e2 , e3 ) ) 
#  @endcode
#  @return formats for nice string and the separate exponent 
def fmt_pretty_errors ( value               ,
                        *errors             ,  
                        width       = 6     ,
                        precision   = 4     ,
                        with_sign   = True  ,
                        latex       = False ) : 
    """ Formats for nice printout of the object with errors  ( strings + exponent)
    >>> fmtv , fmte , expo = fmt_pretty_errors ( number , ( e1 , e2 , e3 )  ) 
    """
    assert isinstance ( width  , integer_types ) and \
        isinstance ( precision , integer_types ) and 2 <= precision < width, \
        "Invalid width/precision parameters: %s/%s" % ( width , precision ) 

    assert isinstance ( value , num_types ) and errors and \
        all ( isinstance ( e , num_types ) for e in errors ) , "Invalid value/errors!"
    
    from ostap.math.base import iszero, frexp10 

    v = float ( value ) 
    e = max   ( abs ( float ( q ) ) for q in errors ) if errors else abs ( v ) 
    
    ## quantity that actually defines the format 
    av = max ( abs ( v ) , e )

    vformat = '%%+%d.%df' if with_sign else '%%%d.%df'
    eformat = '%%-%d.%df'
    
    if   100 <= av < 1000 :
        
        fmtv  = vformat %  ( width , precision - 2 )
        fmte  = eformat %  ( width , precision - 2 )        
        return fmtv, fmte , 0
    
    elif 10  <= av < 100  :
        
        fmtv  = vformat  %  ( width , precision - 1 )
        fmte  = eformat  %  ( width , precision - 1 )
        return fmtv  , fmte , 0
    
    elif 0.1 <= av < 10 :
        
        fmtv  = vformat  %  ( width , precision     )
        fmte  = eformat  %  ( width , precision     )
        return fmtv , fmte , 0

    if iszero ( av ) :
        fmtv  = vformat %  ( width , precision     )
        fmte  = eformat %  ( width , precision     )
        return fmtv , fmte , 0 

    #
    ## here we scale input data and try to get formats for scaled data
    #
    
    v_a , v_e = frexp10 ( av )
    v_ee   = v_e - 1
    n , r  = divmod  ( v_ee , 3 )    

    scale  = 10** ( r - v_ee  )    
    v     *= scale
    errs   = tuple ( e * scale for e in errors ) 

    #
    ## get formats for scaled data
    # 
    fmtv , fmte , expo = fmt_pretty_errors ( v                     ,
                                             *errs                 , 
                                             width     = width     ,
                                             precision = precision ,
                                             with_sign = with_sign ,
                                             latex     = latex     )
    
    return fmtv , fmte , expo + 3 * n 

# =============================================================================
## Format for nice printout of the floating number (string + exponent)
#  @code
#  fmt , expo = fmt_pretty_float ( number ) 
#  @endcode
#  @return format for nice string and the separate exponent 
def fmt_pretty_float ( value             , * , 
                       width     = 6     ,
                       precision = 4     ,
                       with_sign = True  ,
                       latex     = False ) :
    """ Format for nice printout of the floating number
    - return format for nice string and the separate exponent 
    >>> fmt , expo = fmt_pretty_float ( number ) 
    """
    
    assert isinstance ( value , num_types ) , "Invalid `value' parameter: %s"   % typename ( value ) 
    
    value = float ( value )

    ## finite?
    if   math.isinf ( value ) : return "%s" , 0
    elif math.isnan ( value ) : return "%s" , 0

    ## get the format 
    return fmt_pretty_values ( value                 ,
                               width     = width     ,
                               precision = precision ,
                               with_sign = with_sign ,
                               latex     = latex     ) 

# ===============================================================================
## Formats nice printout of the valuw with (sinle, symmetric)  error
#  @code
#  fmt, fmtv, fmte, expo = fmt_pretty_error ( value , error  ) 
#  @endcode
#  @return formats nice string and the separate exponent 
def fmt_pretty_error ( value               ,
                       error               , * ,  
                       width       = 6     ,
                       precision   = 4     ,
                       with_sign   = True  , 
                       parentheses = True  ,
                       latex       = False ) :
    """ Formats for nice printout of the value with (single symmetric) error 
    - return formats for nice stirng and the separate exponent 
    >>> fmt , fmt_v , fmt_e , n = pretty_error ( value , error  ) 
    """    
    assert isinstance ( value , num_types ) , "Invalid `value' parameter: %s"   % typename ( value ) 
    assert isinstance ( error , num_types ) , "Invalid `error' parameter: %s"   % typename ( error ) 
    
    v = float ( value )
    e = max   ( 0.0 , float ( error  ) )
    
    ## get the formats 
    fmtv , fmte , expo = fmt_pretty_errors ( v                     ,
                                             e                     , 
                                             width     = width     ,
                                             precision = precision ,
                                             with_sign = with_sign ,
                                             latex     = latex     ) 
    
    if latex : fmt = '%s %s %s' % ( fmtv , '\\pm'     , fmte )
    else     : fmt = '%s %s %s' % ( fmtv , plus_minus , fmte )
    
    if ( expo and latex ) or parentheses : fmt = '( ' + fmt + ' )'
    
    return fmt, fmtv, fmte , expo

# ===============================================================================
## Formats for nice printout of the object with asymmetric  errors ( string + exponent)
#  @code
#  fmt , fmtv , fmte , expo = fmt_pretty_error2 ( number , errlow  , errhigh ) 
#  @endcode
#  @return formats for nice string and the separate exponent 
def fmt_pretty_error2 ( value               ,
                        errlow              ,
                        errhigh             , * , 
                        width       = 6     ,
                        precision   = 4     ,
                        with_sign   = True  , 
                        parentheses = True  ,
                        latex       = False ) :
    """ Formats for nice printout of the object with asymmetric  errors ( string + exponent)
    >>> fmt , fmtv , fmte , expo = fmt_pretty_err2 ( number , elow  , ehigh ) 
    """ 
    assert isinstance ( value     , num_types ) , "Invalid `value'  parameter: %s"   % typename ( value  ) 
    assert isinstance ( errlow    , num_types ) , "Invalid `errlow' parameter: %s"   % typename ( errlow  ) 
    assert isinstance ( errhigh   , num_types ) , "Invalid `errigh' parameter: %s"   % typename ( errhigh ) 
    
    v  = float ( value   )
    eh = float ( errhigh )
    el = float ( errlow  )

    ## get the formats 
    fmtv , fmte , expo = fmt_pretty_errors ( v , el , eh ,
                                             width     = width     ,
                                             precision = precision ,
                                             with_sign = with_sign ,
                                             latex     = latex     ) 

    if latex : fmt = '%s_{%s}^{%s}' % ( fmtv , fmte , fmte )
    else     : fmt = '%s -/%s +/%s' % ( fmtv , fmte , fmte )
    
    if ( latex and expo ) or parentheses : fmt = '( ' + fmt + ' )' 
    
    return fmt, fmtv, fmte , expo

# =============================================================================
## - Formatters 
# =============================================================================

# =============================================================================
## Nice printout of the floating number (string + exponent)
#  @code
#  s , expo = pretty_float ( number ) 
#  @endcode
#  @return nice stirng and the separate exponent 
def pretty_float ( value             , * , 
                   width     = 6     ,
                   precision = 4     ,
                   with_sign = True  ,
                   latex     = False ) :
    """ Nice printout of the floating number
    - return nice string and the separate exponent 
    >>> s , expo = pretty_float ( number )
    """
    assert isinstance ( value     , num_types ) , "Invalid `value'  parameter: %s"   % typename ( value  ) 

    v = float ( value )
    
    if   math.isinf ( v ) : return '+inf' if 0 < v else '-inf' , 0    
    elif math.isnan ( v ) : return 'NaN' , 0

    ## get format 
    fmt , expo = fmt_pretty_float ( value     = v         ,
                                    width     = width     ,
                                    precision = precision ,
                                    with_sign = with_sign ,
                                    latex     = latex     )
    if expo :
        scale = 10 ** expo
        v /= scale
        
    return  fmt % v , expo 

# ===============================================================================
## nice printout of value with (singular, symetric) error  ( string + exponent)
#  @code
#  s , expo = pretty_err1 ( number , error ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_error ( value               ,
                   error               , * , 
                   width       = 6     ,
                   precision   = 4     ,
                   parentheses = True  ,
                   latex       = False ) :
    """ Nice printout of the value with singuar, symmetric, error ( string + exponent)
    - return nice stirng and the separate exponent 
    >>> s , expo = pretty_err1 ( value , error ) 
    """
    assert isinstance ( value , num_types ) , "Invalid `value'  parameter: %s"   % typename ( value  ) 
    assert isinstance ( error , num_types ) , "Invalid `error'  parameter: %s"   % typename ( error  ) 
    
    v =             float ( value ) 
    e = max ( 0.0 , float ( error ) ) 
    
    finv = not ( math.isinf ( v ) or math.isnan ( v ) ) 
    fine = not ( math.isinf ( e ) or math.isnan ( e ) ) 

    pm = plus_minus if not latex else '\\pm'
    
    if not finv and not fine  :
        fmt = '%%+%ds %s %%-%ds' % ( width , pm , width ) 
        if latex or parentheses : fmt = '( %s )' % fmt 
        return fmt % ( v , e ) , 0 
    elif not finv  :
        fe , ne = pretty_float ( value = e , width = width , precision = precision , latex = latex ) 
        fmt = '%%+%ds %s %%-%ds' % ( width , pm , width ) 
        if ( latex and ne ) or  parentheses : fmt = '( %s )' % fmt 
        return fmt % ( v , fe ) , ne 
    elif not fine  :
        fv , nv = pretty_float ( value = v , width = width , precision = precision , latex = latex ) 
        fmt = '%%%ds %s %%-%ds' % ( width , pm  , width ) 
        if ( latex and nv ) or parentheses : fmt = '( %s )' % fmt 
        return fmt % ( fv , e ) , nv

    ## get format 
    fmt , _ , _ , expo = fmt_pretty_error ( v , e                     , 
                                            width       = width       ,
                                            precision   = precision   ,
                                            parentheses = parentheses ,
                                            latex       = latex       )
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
def pretty_error2 ( value               ,
                    errlow              ,
                    errhigh             , * , 
                    width       = 6     ,
                    precision   = 4     ,
                    parentheses = True  ,
                    latex       = False ) :
    """ Nice printout of the object with asymmetric  errors   ( string + exponent)
    >>> s , expo = pretty_err2 ( number , elow , ehigh ) 
    """

    assert isinstance ( value   , num_types ) , "Invalid `value'   parameter: %s" % typename ( value  ) 
    assert isinstance ( errlow  , num_types ) , "Invalid `errlow'  parameter: %s" % typename ( errlow  ) 
    assert isinstance ( errhigh , num_types ) , "Invalid `errhigh' parameter: %s" % typename ( errhigh ) 

    v     = float ( value   )
    ehigh = float ( errhigh )
    elow  = float ( errlow  )

    assert 0 <= elow or 0 <= ehigh, 'Both errlow/errhigh cannot be negative!'

    if   elow  <= 0          <= ehigh : pass
    elif ehigh <= 0          <= elow  : elow, ehigh = ehigh, elow 
    elif 0     <= elow and 0 <= ehigh : elow = -elow

    ## get the format 
    fmt , _ , _ , expo = fmt_pretty_error2 ( v , elow , ehigh ,
                                             width       = width       ,
                                             precision   = precision   ,
                                             parentheses = parentheses ,
                                             latex       = latex       )
    
    values = value , abs ( elow ) , ehigh    
    if expo :
        scale = 10 ** expo 
        values = tuple ( v / scale for v in values ) 
    
    return fmt % values , expo 

# ===============================================================================
## nice printout of the asymmetric  errors  ( string + exponent)
#  @code
#  s , expo = pretty_errors ( errlow , errhigh ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_errors ( errlow              ,
                    errhigh             , * , 
                    width       = 6     ,
                    precision   = 4     ,
                    parentheses = True  ,
                    latex       = False ) :
    """ Nice printout of asymmetric  errors   ( string + exponent)
    >>> s , expo = pretty_errors ( errlow , errhigh ) 
    """
    assert isinstance ( errlow  , num_types ) , "Invalid `errlow' parameter: %s" % typename ( errlow  ) 
    assert isinstance ( errhigh , num_types ) , "Invalid `errigh' parameter: %s" % typename ( errhigh ) 
    
    ehigh = float ( errhigh )
    elow  = float ( errlow  )

    assert 0 <= elow or 0 <= ehigh, 'Both errlow/errhigh cannot be negative!'

    if   elow  <= 0          <= ehigh : pass
    elif ehigh <= 0          <= elow  : elow, ehigh = ehigh, elow 
    elif 0     <= elow and 0 <= ehigh : elow = -elow

    ## get the format 
    _ , fmte , expo = fmt_pretty_errors ( 0.0 , elow , ehigh    ,
                                          width     = width     ,
                                          precision = precision ,
                                          latex     = latex     )
    
    ## construct format 
    if latex : fmt = '_{%s}^{%s}' % ( fmte , fmte )
    else     : fmt = ' -/%s +/%s' % ( fmte , fmte )
    if ( latex and expo ) or parentheses : fmt = '( ' + fmt + ' )' 

    values = value , abs ( elow ) , ehigh    
    if expo :
        scale = 10 ** expo 
        values = tuple ( v / scale for v in values ) 
    
    return fmt % values , expo 

# =============================================================================
##  Add the exponent to the string representation of the obejuct
#   @code
#   value = ...
#   value , expo = pretty_XXX ( value , .... )
#   result = add_expo ( value , expo ) 
#   @endcode 
def add_expo ( value , expo , fmt = '%%s %s 10^%%+d' % times ) :
    """ Add the exponent to the sting representation of the objeuct
    >>> value = ...
    >>> value , expo = pretty_XXX ( value , .... )
    >>> result = add_expo ( value , expo ) 
    """
    assert isinstance ( expo , integer_types ) , "Invalid type of `expo':%s" % typename ( expo )
    if not expo : return value
    return fmt % ( value , expo )

# =======================================================================
## nice printout of asymmetric errors ( string + exponent)
#  @code
#  ae = AsymErrors ( ... ) 
#  s , expo = pretty_ae (  ae ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_ae  ( errors              , * , 
                 width       = 6     ,
                 precision   = 4     ,
                 parentheses = True  ,
                 latex       = False ) :
    """ Nice printout of asymmetric errors ( string + exponent)
    - return nice stirng and the separate exponent 
    >>> s , expo = pretty_ve ( number ) 
    """
    from ostap.utils.valerrors  import pretty_ae as _pretty_ae_ 
    return _pretty_ae_ ( errors      = errors      ,
                         width       = width       ,
                         precision   = precision   ,
                         parentheses = parentehses ,
                         latex       = latex       ) 

# =======================================================================
## nice printout of the ValWithErrors object  ( string + exponent)
#  @code
#  vae = ValWithErors ( ... ) 
#  s , expo = pretty_vae (  ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_vae ( value               , * , 
                 width       = 6     ,
                 precision   = 4     ,
                 parentheses = True  ,
                 latex       = False ) :
    """ Nice printout of the ValueWithError object  ( string + exponent)
    - return nice stirng and the separate exponent 
    >>> s , expo = pretty_ve ( number ) 
    """
    from ostap.utils.valerrors  import pretty_vae as _pretty_vae_ 
    return _pretty_vae_ ( value       = value       ,
                          width       = width       ,
                          precision   = precision   ,
                          parentheses = parentheses ,
                          latex       = latex       ) 

# =======================================================================
## Pretty print for the object
#  @code
#  obj = ...
#  prnt , expo = pretty_print ( obj ) 
#  @endcode
def pretty_print ( what              , * , 
                   width     = 6     ,
                   precision = 4     ,
                   latex     = False , 
                   **kwargs          ) :
    """ Pretty print for the object
    >>> obj = ...
    >>> prnt , expo = pretty_print ( obj ) 
    """
    if  hasattr ( what , 'pretty_print' ) :
        return what.pretty_print ( width     = width     ,
                                   precision = precision ,
                                   lates     = latex     , **kwargs )
    elif isinstance ( what , num_types ) :
        return pretty_float ( width     = width     ,
                              precision = precision ,
                              latex     = latex     , **kwargs )
    
    from ostap.utils.basic import typename 
    raise TypeError( "pretty_print: invalid `what` type: %s:" % typename ( what ) )

# =======================================================================
## Nice print for the object
#  @code
#  obj  = ...
#  prnt = nice_print ( obj ) 
#  @endcode
def nice_print ( what              , * , 
                 width     = 6     ,
                 precision = 4     ,
                 latex     = False ,
                 **kwargs          ) :
    """ Nice print for the object
    >>> obj = ...
    >>> prnt = nice_print ( obj ) 
    """
    if  hasattr ( what , 'nice_print' ) :
        return what.nice_print ( width     = width     ,
                                 precision = precision ,
                                 latex     = latex     , **kwargs )

    result , expo = pretty_print ( what                  ,
                                   width     = width     ,
                                   precision = precision ,
                                   latex     = latex     , **kwargs )
    
    if   expo and latex : result = '%s %s 10^{%+d}' % ( result , '\\times' , expo )
    elif expo           : result = '%s%s10^%+d'     % ( result ,    times  , expo )
    ##
    return result


# =======================================================================
## Pretty LaTeX print for the object
#  @code
#  obj = ...
#  prnt , expo = pretty_latex ( obj ) 
#  @endcode
def pretty_latex ( what              , * , 
                   width     = 6     ,
                   precision = 4     ,
                   **kwargs          ) :
    """ Pretty LaTeX print for the object
    >>> obj = ...
    >>> prnt , expo = pretty_latex ( obj ) 
    """
    kwargs.pop ('latex', None  )
    return pretty_print ( what ,
                          width     = width     ,
                          precision = precision ,
                          latex     = True      , **kwargs )

# =======================================================================
## Nice print for the object
#  @code
#  obj  = ...
#  prnt = nice_latex ( obj ) 
#  @endcode
def nice_latex  ( what              , * , 
                  width     = 6     ,
                  precision = 4     ,
                  **kwargs          ) :
    """ Nice print for the object
    >>> obj = ...
    >>> prnt = nice_latex ( obj ) 
    """
    kwargs.pop ('latex', None  )
    return nice_print ( what ,
                        width     = width     ,
                        precision = precision ,
                        latex     = True      , **kwargs )

# =======================================================================
## nice printout of the ValWithMultiErrors object  ( string + exponent)
#  @code
#  vae = ValWithMultiErros ( ... ) 
#  s , expo = pretty_vme (  ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_vme ( value               , * , 
                 width       = 6     ,
                 precision   = 4     ,
                 parentheses = True  ,
                 latex       = False ) :
    """ Nice printout of the ValueWithMultiError object  ( string + exponent)
    - return nice stirng and the separate exponent 
    >>> s , expo = pretty_ve ( number ) 
    """
    from ostap.utils.valerrors  import pretty_vme as _pretty_vme_ 
    return _pretty_vme_ ( value       = value       ,
                          width       = width       ,
                          precision   = precision   ,
                          parentheses = parentehses ,
                          latex       = latex       ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
