#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file
# Set of utilities for rounding according to PDG prescription
# @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
# @see section 5.3 of doi:10.1088/0954-3899/33/1/001
# Quote:
# The basic rule states that:
# - if the three highest order digits of the error lie between 100 and 354,
#   we round to two significant digits. 
# - If they lie between 355 and 949,
#   we round to one significant digit.
# - Finally, if they lie between 950 and 999,
#   we round up to 1000 and keep two significant digits.
# 
# In all cases, the central value is given with a precision that matches that of the error.
# 
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2015-07-20
#
# =============================================================================
""" Set of utilities for rounding according to PDG prescription
see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
see section 5.3 of doi:10.1088/0954-3899/33/1/001

Quote:

The basic rule states that:
 - if the three highest order digits of the error lie between 100 and 354,
   we round to two significant digits. 
 - If they lie between 355 and 949,
   we round to one significant digit.
 - Finally, if they lie between 950 and 999,
   we round up to 1000 and keep two significant digits.

In all cases, the central value is given with a precision that matches that of the error.
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2015-07-15"
__version__ = "$Revision$"
__all__ = (
    'pdg_round'    , ## round value,error-pair according to PDG prescription
    'pdg_format_'  , ## format value&error according to PDG prescription
    'pdg_format'   , ## format value+2errors according to PDG
    )
# ===============================================================================
from   ostap.math.math_base   import frexp10
from   ostap.core.ostap_types import integer_types, string_types
from   ostap.logger.pretty    import pretty_float 
from   ostap.logger.symbols   import times, plus_minus
import math, enum  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.utils.pdg_format' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
## the reference error type (in case seevral errors are specified) 
class ErrMode(enum.IntEnum):
    TOTAL      = 0  ## Use total uncertainty 
    MIN        = 1  ## Use minimal uncertainty 
    MAX        = 2  ## Use maximal uncertainty 
    MEAN       = 3  ## Use mean
    AVERAGE    = 3  ## Use mean
    QUADRATIC  = 4  ## Use quadratic (root mean square) 
    RMS        = 4  ## Use quadratic (root mean square) 

# =============================================================================
## get the `reference/representative error' from the list of uncertainties
#  @code
#  error = ref_error ( 'total'  , 0.1 , 0.2 , 0.3 ) 
#  error = ref_error ( Mode.MIN , 0.1 , 0.2 , 0.3 ) 
#  @endcode 
def ref_error ( mode , error , *errors ) :
    """ Get the `reference/representative error' from the list of uncertainties
    >>> error = ref_error ( 'total'  , 0.1 , 0.2 , 0.3 ) 
    >>> error = ref_error ( Mode.MIN , 0.1 , 0.2 , 0.3 ) 
    """

    if not errors  : return abs ( error ) 

    if not mode : mode = ErrMode.TOTAL 

    umode = mode

        
    if isinstance ( mode , string_types ) :
        
        umode = mode.upper()
        
        if   umode in ( 'MIMIMAL'    , 'MINIMUM' , 'MIN'     , 'MN'        ) : umode = ErrMode.MIN      .name 
        elif umode in ( 'MAXIMAL'    , 'MAXIMUM' , 'MAX'     , 'MX'        ) : umode = ErrMode.MAX      .name 
        elif umode in ( 'A' , 'AV'   , 'AVE'     , 'MEAN'                  ) : umode = ErrMode.AVERAGE  .name        
        elif umode in ( 'Q' , 'QUAD' , 'QUADR'   , 'QUADRAT' , 'QUADRATIC' ) : umode = ErrMode.QUADRATIC.name
        elif umode in ( 'R' , 'RMS'                                        ) : umode = ErrMode.QUADRATIC.name
        elif umode                                                           : umode = ErrMode.TOTAL    .name 
        
        assert umode in ErrMode.__members__ , 'ref_error: Unknown string mode: %s' % mode 
        
        umode = ErrMode [ umode ]
        
    elif isinstance ( mode , integer_types ) :

        for _ , v in ErrMode.__members__.items() :
            if v == mode :
                umode = v
                break
        else :
            raise ValueError("ref_error: Unknown integer mode %s" % mode )

    assert isinstance ( umode , ErrMode ), 'ref_error: Unknown mode %s' % umode

    if umode == ErrMode.TOTAL :
        
        result = error * error 
        for e in errors : result += e * e 
        return math.sqrt ( result )
    
    elif umode == ErrMode.MIN :
        
        result = abs ( error )  
        for e in errors :
            ae = abs ( e )
            if ae and 0 < ae :
                result = ae if not result else min ( result , ae )          
        return result
    
    elif umode == ErrMode.MAX :
        
        result = abs ( error )  
        for e in errors : result = max ( result , abs ( e ) )            
        return result
    
    elif umode == ErrMode.QUADRATIC or umode == ErrMode.RMS : 
        
        result = error * error 
        ne     = 1 
        for e in errors :
            result += e * e 
            ne     += 1                 
        return result if 0 == result else math.sqrt ( result / ne )

    
    ## MEAN 
    result = abs ( error )
    ne     = 1 
    for e in errors :
        result += abs ( e ) 
        ne     += 1                 
    return result  / float ( ne ) 
    

# ==============================================================================
#  Classify according to PDG prescription
#  @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
#  @see section 5.3 of doi:10.1088/0954-3899/33/1/001
#
#  The basic rule states that
#   - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
#   - If they lie between 355 and 949, we round to one significant digit.
#   - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
#
#  Cases:
#   -  -3, error is NaN
#   -  -2, error is Inf
#   -  -1, error is `not isfinite`
#   -   0, error is zero 
#   -   1, the first regular case
#   -   2, the second regular case
#   -   3, the third  regular case

#  @code
#  case , error , expo = pdg_case ( error ) 
#  @endcode 
def pdg_case ( error ) :
    """ Classify according to PDG prescription
    - see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    - see section 5.3 of doi:10.1088/0954-3899/33/1/001
    
    The basic rule states that
    - if the three highest order digits of the error lie between 100 and 354
      we round to two significant digits.
    - If they lie between 355 and 949
      we round to one significant digit.
    - Finally, if they lie between 950 and 999
      we round up to 1000 and keep two significant digits.
    
    Cases:
    -  -3, error is NaN
    -  -2, error is Inf 
    -  -1, error is nontfinite 
    -   0, error is zero 
    -   1, the first regular case
    -   2, the second regular case
    -   3, the third  regular case
    
    >>> error = ...
    >>> case , rounded_error , expo = pdg_case ( error ) 
    """
    
    if        math.isnan    ( error ) : return -3 , error , 0
    elif      math.isinf    ( error ) : return -2 , error , 0  
    elif  not math.isfinite ( error ) : return -1 , error , 0 
    
    m , expo = frexp10 ( error )
  
    am    = abs ( m )  

    ##   
    if   am <  0.100 : return 0 ,         m       ,  expo   
    elif am <= 0.354 : return 1 , round ( m , 2 ) , expo 
    elif am <= 0.949 : return 2 , round ( m , 1 ) , expo  
    #
    return 3 , 1.0 , expo  

# ==============================================================================

# ==============================================================================
## Round value and error according to PDG prescription
#  @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
#  @see section 5.3 of doi:10.1088/0954-3899/33/1/001
#
#  The basic rule states that
#   - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
#   - If they lie between 355 and 949, we round to one significant digit.
#   - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
#
#  @code
#  value, error = ...
#  value, error = pdg_round ( value , error )
#  @ndcode
def pdg_round ( value , error ) :
    """ Round value and error according to PDG prescription
    - see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    - see section 5.3 of doi:10.1088/0954-3899/33/1/001
    
    The basic rule states that
    - if the three highest order digits of the error lie between 100 and 354
      we round to two significant digits.
    - If they lie between 355 and 949
      we round to one significant digit.
    - Finally, if they lie between 950 and 999
      we round up to 1000 and keep two significant digits.
    """
   
    the_case , the_error , expo = pdg_case ( error )
    if  the_case <= 0 : return value , error 
    
    scale = 10 ** expo if expo else 1 
    
    vv    = value / scale 
    ev    = the_error
    
    if   1 == the_case : vv = round ( vv , 2 )
    elif 2 == the_case : vv = round ( vv , 1 )
    elif 3 == the_case : vv = round ( vv , 1 )
    
    return vv * scale , ev * scale 

# ============================================================================
## Format value and errors according to PDG prescription
#  @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
#  @see section 5.3 of doi:10.1088/0954-3899/33/1/001
#
#  The basic rule states that
#   - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
#   - If they lie between 355 and 949, we round to one significant digit.
#   - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
#
#  @code
#  value, error = ...
#  case , fmtv, fmte, expo = fmt_pdg ( value , error )
#  @ndcode
def fmt_pdg ( value , error ) :
    """ Format value and error according to PDG prescription
    - see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    - see section 5.3 of doi:10.1088/0954-3899/33/1/001
    
    The basic rule states that
    - if the three highest order digits of the error lie between 100 and 354
      we round to two significant digits.
    - If they lie between 355 and 949
      we round to one significant digit.
    - Finally, if they lie between 950 and 999
      we round up to 1000 and keep two significant digits.
      
    >>> value , error = ...
    >>> case , fmtv, fmte, expo = fmt_pdg  ( value , error )
    """
    
    the_case , the_error , expo = pdg_case ( error )
    if  the_case <= 0 : return the_case , '%+g', '%g' , expo
        
    scale = 10 ** expo if expo else 1 
    vv    = value / scale
    
    av    = abs ( vv ) 
    ev    = the_error

    rr = expo % 3
    
    r0 = 0 == rr
    r1 = 1 == rr
    r2 = 2 == rr
    
    if   1 == the_case :
        
        if   r1             : return the_case  , '%+.1f' , '%.1f' , expo - 1 
        elif r2 and 1 <= av : return the_case  , '%+.3f' , '%.3f' , expo + 1 
        elif r2             : return the_case  , '%+.0f' , '%.0f' , expo - 2 
        return the_case  , '%+.2f' , '%.2f' , expo 

    if   r1 and 10 <= av : return the_case  , '%+.3f' , '%.3f' , expo + 2 
    elif r1              : return the_case  , '%+.0f' , '%.0f' , expo - 1    
    elif r2              : return the_case  , '%+.2f' , '%.2f' , expo + 1 
    return the_case , '%+.1f' , '%.1f' , expo 

       
# ============================================================================
## Format value and errors according to PDG prescription
#  @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
#  @see section 5.3 of doi:10.1088/0954-3899/33/1/001
#
#  The basic rule states that
#   - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
#   - If they lie between 355 and 949, we round to one significant digit.
#   - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
#
#  @code
#  value, error = ...
#  result , expo = pdg_format_  ( value , error )
#  @ndcode     
def pdg_format_ ( value         , *errors , 
                  mode          = 'TOTAL' , 
                  latex         = False   , 
                  neglect_error = 1.1e-8  ) :
    """ Format value and error according to PDG prescription
    - see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    - see section 5.3 of doi:10.1088/0954-3899/33/1/001
    
    The basic rule states that
    - if the three highest order digits of the error lie between 100 and 354
      we round to two significant digits.
    - If they lie between 355 and 949
      we round to one significant digit.
    - Finally, if they lie between 950 and 999
      we round up to 1000 and keep two significant digits.
      
    >>> value , error = ...
    >>> result , expo = pdg_format_ ( value , error )
    """  
    sumerr2  = sum (  float ( e ) * float ( e ) for e in errors )
    if not errors or not sumerr2 or ( 0 < neglect_error and math.sqrt ( sumerr2 ) <= neglect_error * abs ( value ) ) : 
        return pretty_float ( value , latex = latex , precision = 6 )
    
    ## get the representative error 
    the_error = ref_error ( mode , *errors )
    
    pm = ' %s ' %  ( '\\pm' if latex else plus_minus ) 
        
    the_case , fmtv , fmte , expo = fmt_pdg ( value , the_error )
    if  the_case <= 0 :
        result  = ( fmtv % value ) + pm 
        result += pm.join ( ( fmte % e ) for e in errors )  
        return result , expo 
    
    scale   = 10 ** expo if expo else 1      
    result  = ( fmtv %  ( value / scale ) ) + pm 
    result += pm.join ( ( fmte % ( e / scale ) ) for e in errors )  
    
    return result , expo 
 
# ============================================================================
## Format value and errors according to PDG prescription
#  @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
#  @see section 5.3 of doi:10.1088/0954-3899/33/1/001
#
#  The basic rule states that
#   - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
#   - If they lie between 355 and 949, we round to one significant digit.
#   - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
#
#  @code
#  value, error = ...
#  result = pdg_format ( value , error )
#  @ndcode     
def pdg_format ( value         , *errors , 
                 mode          = 'TOTAL' ,
                 latex         = False   ,
                 neglect_error = 1.1e-8  ) :
    
    """ Format value and error according to PDG prescription
    - see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    - see section 5.3 of doi:10.1088/0954-3899/33/1/001
    
    The basic rule states that
    - if the three highest order digits of the error lie between 100 and 354
      we round to two significant digits.
    - If they lie between 355 and 949
      we round to one significant digit.
    - Finally, if they lie between 950 and 999
      we round up to 1000 and keep two significant digits.
      
    >>> value , error = ...
    >>> result = pdg_format ( value , error )
    """
    ##
    sumerr2  = sum (  float ( e ) * float ( e ) for e in errors )
    if not errors or not sumerr2 or ( 0 < neglect_error and math.sqrt ( sumerr2 ) <= neglect_error * abs ( value ) ) : 
        result , expo = pretty_float ( value , precision = 6 , latex = latex  )
    else :    
        result , expo  = pdg_format_ ( value, *errors , mode = mode, latex = latex , neglect_error = neglect_error )

    if   expo and latex : result = '%s %s 10^{%+d}' % ( result , '\\times' , expo )
    elif expo           : result = '%s %s 10^%+d'   % ( result ,    times  , expo )
    ##
    return result

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
