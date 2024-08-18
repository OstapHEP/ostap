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
    'pretty_float'       , ## pretty print of the floating number 
    'pretty_err1'        , ## pretty print of the floating number with error
    'pretty_err2'        , ## pretty print of the floating number with asymmetric error
    'pretty_errs'        , ## pretty print of the asymmetric error
    ## 
    'pretty_ve'          , ## pretty print of VE  / ValueWithError      object
    'pretty_ae'          , ## pretty print of AE  / AsymErrors          object
    'pretty_vae'         , ## pretty print of VAE / ValWithErrors       object
    'pretty_vme'         , ## pretty print of VME / ValueWithMultiError object
    ##
    'add_expo'           , ## add an exponetaion factor for sting representaion of the object
) 
# =============================================================================
from ostap.logger.utils import  ( pretty_float    ,
                                  pretty_err1     ,
                                  pretty_err2     ,
                                  pretty_errs     ,
                                  add_expo        , 
                                  fmt_pretty_err1 )
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.logger.pretty' )
else                       : logger = getLogger( __name__                )
# =============================================================================
logger.debug ( "Helper functions for pretty prints of some object" )
# =============================================================================
## Formats f0r nice printout of the ValueWithError object  ( string + exponent)
#  @code
#  fmt , fmtv , fmte , expo = fmt_pretty_ve ( number ) 
#  @endcode
#  @return nice string and the separate exponent 
def fmt_pretty_ve ( value              ,
                    width       = 6    ,
                    precision   = 4    ,
                    parentheses = True ) :
    """ormats f0r nice printout of the ValueWithError object  ( string + exponent)
    >>> fmt , fmtv , fmte , expo = fmt_pretty_ve ( number ) 
    """
    from ostap.math.ve import VE
    assert isinstance ( value , VE ) , \
        "Invalid `value' parameter: %s" % type ( value )
    ## decode object 
    v , e =  value.value () , max ( 0 , value.error () )
    return fmt_pretty_err1 ( value       = v           ,
                             error       = e           ,
                             width       = width       ,
                             precision   = precision   ,
                             parentheses = parentheses ) 

# =============================================================================
## nice printout of the ValueWithError object  ( string + exponent)
#  @code
#  s , expo = pretty_ve ( number ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_ve ( value              ,
                width       = 6    ,
                precision   = 4    ,
                parentheses = True ) :
    """Nice printout of the ValueWithError object  ( string + exponent)
    - return nice stirng and the separate exponent 
    >>> s , expo = pretty_ve ( number ) 
    """
    from ostap.math.ve import VE
    assert isinstance ( value , VE ) , \
        "Invalid `value' parameter: %s" % type ( value )
    ## decode object 
    v , e =  value.value () , max ( 0 , value.error () )
    ## delegate 
    return pretty_err1 ( v , e , width = width , precision = precision )

# =======================================================================
## nice printout of asymmetric errors ( string + exponent)
#  @code
#  ae = AsymErrors ( ... ) 
#  s , expo = pretty_ae (  ae ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_ae  ( errors             ,
                 width       = 6    ,
                 precision   = 4    ,
                 parentheses = True ) :
    """Nice printout of asymmetric errors ( string + exponent)
    - return nice stirng and the separate exponent 
    >>> s , expo = pretty_ve ( number ) 
    """
    from ostap.utils.valerrors  import pretty_ae as _pretty_ae_ 
    return _pretty_ae_ ( errors      = errors      ,
                         width       = width       ,
                         precision   = precision   ,
                         parentheses = parentehses ) 

# =======================================================================
## nice printout of the ValWithErrors object  ( string + exponent)
#  @code
#  vae = ValWithErors ( ... ) 
#  s , expo = pretty_vae (  ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_vae ( value              ,
                 width       = 6    ,
                 precision   = 4    ,
                 parentheses = True ) :
    """Nice printout of the ValueWithError object  ( string + exponent)
    - return nice stirng and the separate exponent 
    >>> s , expo = pretty_ve ( number ) 
    """
    from ostap.utils.valerrors  import pretty_vae as _pretty_vae_ 
    return _pretty_vae_ ( value       = value       ,
                          width       = width       ,
                          precision   = precision   ,
                          parentheses = parentehses ) 

# =======================================================================
## nice printout of the ValWithMultiErrors object  ( string + exponent)
#  @code
#  vae = ValWithMultiErros ( ... ) 
#  s , expo = pretty_vme (  ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_vme ( value              ,
                 width       = 6    ,
                 precision   = 4    ,
                 parentheses = True ) :
    """Nice printout of the ValueWithError object  ( string + exponent)
    - return nice stirng and the separate exponent 
    >>> s , expo = pretty_ve ( number ) 
    """
    from ostap.utils.valerrors  import pretty_vme as _pretty_vme_ 
    return _pretty_vme_ ( value       = value       ,
                          width       = width       ,
                          precision   = precision   ,
                          parentheses = parentehses ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
