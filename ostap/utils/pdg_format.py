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
    'round_N'      , ## round floating value to N-significant digits
    'pdg_round'    , ## round value,error-pair according to PDG prescription
    'pdg_format_'  , ## format value&error according to PDG prescription
    'pdg_format'   , ## format value+2errors according to PDG
    )
# ===============================================================================
from   ostap.math.ve          import VE
from   ostap.math.base        import ( frexp10  , round_N , 
                                       isfinite , isclose , iszero, isequal )
from   ostap.core.ostap_types import integer_types, string_types
from   ostap.logger.symbols   import times, plus_minus 
import ROOT,  math, sys, enum  
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
    GEOMETRIC  = 4  ## Use geometric mean
    QUADRATIC  = 5  ## Use quadratic (root mean square) 
    RMS        = 5  ## Use quadratic (root mean square) 

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
        
    umode = mode

    if isinstance ( mode , string_types ) :
        
        umode = mode.upper()
        
        if   umode in ( 'MIMIMAL'    , 'MINIMUM' , 'MIN'     , 'MN'        ) : umode = ErrMode.MIN      .name 
        elif umode in ( 'MAXIMAL'    , 'MAXIMUM' , 'MAX'     , 'MX'        ) : umode = ErrMode.MAX      .name 
        elif umode in ( 'A' , 'AV'   , 'AVE'     , 'MEAN'                  ) : umode = ErrMode.AVERAGE  .name        
        elif umode in ( 'G' , 'GEO'  , 'GEOM'    , 'GEOMET'  , 'GEOMETRIC' ) : umode = ErrMode.GEOMETRIC.name
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
            if not e : continue  
            ae = abs ( e ) 
            result = ae if not result else min ( result , ae )          
        return result
    
    elif umode == ErrMode.MAX :
        
        result = abs ( error )  
        for e in errors : result = max ( result , abs ( e ) )            
        return result
    
    elif umode == ErrMode.GEOMETRIC :
        
        result = abs ( error )
        ne     = 1 
        for e in errors :
            ae = abs ( e )
            if 0 < ae :
                result *= ae
                ne     += 1                 
        return result if 0 == result else pow ( result , 1.0 / ne )
    
    elif umode == ErrMode.QUADRATIC :
        
        result = error * error 
        ne     = 1 
        for e in errors :
            result += e * e 
            ne     += 1                 
        return result if 0 == result else math.sqrt ( result / ne )
    
    else :
        
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
#  fmtv, fmte, expo = pdg_fmt  ( value , error )
#  @ndcode
def pdg_fmt ( value , error ) :
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
    >>> fmtv, fmte, expo = pdg_fmt  ( value , error )
    """
    
    the_case , the_error , expo = pdg_case ( error )
    if  the_case <= 0 : return the_case , '%+g', '%g' , expo
        
    scale = 10 ** expo if expo else 1 
    vv    = value / scale 
    ev    = the_error
        
    if   1 == the_case : return the_case  , '%+.2f' , '%.2f' , expo 
        
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
def pdg_format_ ( value  , *errors , 
                  mode   = 'TOTAL' , 
                  latex  = False   ) :
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
  
    if not errors :
        result = '%+g' % value
        if parentheses : result = '( %s )' % result 
        return result , 0 
    
    ## get the representative error 
    the_error = ref_error ( mode , *errors )
    
    pm = ' %s ' %  ( '\\pm' if latex else plus_minus ) 
        
    the_case , fmtv , fmte , expo = pdg_fmt ( value , the_error )
    if  the_case <= 0 :
        result  = ( fmtv % value ) + pm 
        result += pm.join ( ( fmte % e ) for e in errors )  
        return result , expo 
    
    scale = 10 ** expo if expo else 1  
    
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
def pdg_format ( value , *errors , 
                 mode  = 'TOTAL' ,
                 latex = False   ) :
    
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
    result, expo  = pdg_format_ ( value, *errors , mode = mode, latex = latex )
    if not expo : return result 
    
    if latex : return '( %s )\\times 10^{%d}' % ( result , expo )
    return '( %s ) %s 10^%+d' %  ( result , times , expo )
     
# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    N = 37

    rows = [ ( 'n1' , 'n2' , 'n3' , 'n4' , 'n5' ) ] 
    for i in range ( 1 , N ) :
        
        v  = float(i)/N
         
        v1 = VE ( v        , 0.01 * v*v             )
        v2 = VE ( v        , 0.01 * v*v/100         )
        v3 = VE ( v        , 0.01 * v*v/10000       )
        v4 = VE ( v        , 0.01 * v*v/1000000     )
        v5 = VE ( v        , 0.01 * v*v/100000000   )
        v6 = VE ( v        , 0.01 * v*v/10000000000 )
        v7 = VE ( math.inf , 0.01 * v*v/10000000000 )
        v8 = VE ( v        , math.nan        )
        
        for e in [ -100 , -50 , -10 ] + [ g for g in range ( -5 , 6 ) ] + [ 10 , 50 , 100 ] : 
            
            w1 = v1 * 10 ** e 
            w2 = v2 * 10 ** e 
            w3 = v3 * 10 ** e 
            w4 = v4 * 10 ** e 
            w5 = v5 * 10 ** e 
            w6 = v6 * 10 ** e 
            w7 = v7 * 10 ** e 
            w8 = v8 * 10 ** e 
            
            row = ( pdg_format ( w1.value() , w1.error() ) ,
                    pdg_format ( w2.value() , w2.error() ) ,
                    pdg_format ( w3.value() , w3.error() ) ,
                    pdg_format ( w4.value() , w4.error() ) ,
                    pdg_format ( w5.value() , w5.error() ) ,
                    pdg_format ( w6.value() , w6.error() ) ,
                    pdg_format ( w7.value() , w7.error() ) ,
                    pdg_format ( w8.value() , w8.error() ) ) 
            
            rows.append  ( row )
            
    import ostap.logger.table as T
    logger.info ( 'PDG roundings:\n%s ' % T.table ( rows , prefix = '# ' , alignment = 'ccccc' ) ) 
                  
    logger.info ( 80*'*' )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
