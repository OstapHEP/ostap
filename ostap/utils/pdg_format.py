#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file
# Set of utilities for rounding according to PDG prescription
# @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
# @see section 5.3 of doi:10.1088/0954-3899/33/1/001
# Quote:
# The basic rule states that if the three highest order digits of the error
# lie between 100 and 354, we round to two significant digits. If they lie between
# 355 and 949, we round to one significant digit. Finally,
# if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
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
| The basic rule states that if the three highest order digits of the error
| lie between 100 and 354, we round to two significant digits. If they lie between
| 355 and 949, we round to one significant digit. Finally,
| if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
| In all cases, the central value is given with a precision that matches that of the error.

"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2015-07-15"
__version__ = "$Revision$"
__all__ = (
    'round_N'     , ## round floating value to N-significant digits
    'pdg_round'   , ## round value,error-pair according to PDG prescription
    'pdg_format'  , ## format value&error according to PDF prescription
    'pdg_format2' , ## format value+2errors according to PDG
    'pdg_format3' , ## format value+3errors according to PDG
    )
# ===============================================================================
from   ostap.math.ve          import VE
from   ostap.math.base        import frexp10, isfinite, isclose  
from   ostap.core.ostap_types import integer_types, string_types 
import ROOT,  math, sys, enum  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.utils.pdg_format' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
class ErrMode(enum.IntEnum):
    TOTAL      = 0  ## Use total uncertainty 
    MIN        = 1  ## Use minimal uncertainty 
    MAX        = 2  ## Use maximal uncertainty 
    MEAN       = 3  ## Use mean
    AVERAGE    = 3  ## Use mean
    GEOMETRIC  = 4  ## Use geometric mean
    QUADRATIC  = 5  ## USe quadratic (root mean sqaure) 
    RMS        = 5  ## USe quadratic (root mean sqaure) 

# =============================================================================
## get the `reference/representative error' from the list of uncertainties
#  @code
#  error = ref_error ( 'total'  , 0.1 , 0.2 , 0.3 ) 
#  error = ref_error ( Mode.MIN , 0.1 , 0.2 , 0.3 ) 
#  @endcode 
def ref_error ( mode , error , *errors ) :
    """Get the `reference/representative error' from the list of uncertainties
    >>> error = ref_error ( 'total'  , 0.1 , 0.2 , 0.3 ) 
    >>> error = ref_error ( Mode.MIN , 0.1 , 0.2 , 0.3 ) 
    """

    if not errors  : return abs ( error ) 
        
    umode = mode

    if isinstance ( mode , string_types ) :
        
        umode = model.upper()
        
        if   umode in ( 'MIMIMAL'    , 'MINIMUM' , 'MIN'     , 'MN'        ) : umode = ErrMode.MIN      .name 
        elif umode in ( 'MAXIMAL'    , 'MAXIMUM' , 'MAX'     , 'MX'        ) : umode = ErrMode.MAX      .name 
        elif umode in ( 'A' , 'AV'   , 'AVE'     , 'MEAN'                  ) : umode = ErrMode.AVERAGE  .name        
        elif umode in ( 'G' , 'GEO'  , 'GEOM'    , 'GEOMET'  , 'GEOMETRIC' ) : umode = ErrMode.QUADRATIC.name
        elif umode in ( 'Q' , 'QUAD' , 'QUADR'   , 'QUADRAT' , 'QUADRATIC' ) : umode = ErrMode.QUADRATIC.name
        elif umode in ( 'R' , 'RMS'                                        ) : umode = ErrMode.QUADRATIC.name
        elif umode                                                           : umode = ErrMode.TOTAL    .name 
        
        assert umode in ErrMode.__members__ ,\
               'ref_error: Unknown string mode: %s' % mode 
        
        umode = ErrMode[umode]
        
    elif isinstance ( mode , integer_types ) :

        for k,v in ErrMode.__members__.items() :
            if v == mode :
                umode = v
                break
        else :
            raise ValueError("ref_error: Unknown integer mode %s" % mode )

    assert isinstance ( umode , ErrMode ),\
           'ref_error: Unknown mode %s' % umode

    if umode == ErrMode.TOTAL :
        
        result = error * error 
        for e in errors : result +=  e * e 
        return math.sqrt ( result )
    
    elif umode == ErrMode.MIN :
        
        result = abs ( error )  
        for e in errors : result = min ( result , abs ( e ) )            
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
        
        result = abs ( error )
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


# =============================================================================
## round to nearest integer, rounds half integers to nearest even integer 
#  It is just a simple wrapper around boost::numeric::converter 
#  @see Ostap::Math::round 
## cpp_round   = Ostap.Math.round
## cpp_round_N = Ostap.Math.round_N 

## ============================================================================
#  round value to N-digits
#  @code
#  new_value = round_N ( value , 3 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20
def round_N ( value , n ) :
    """Round value to to N-digits    
    >>> new_value = round_N ( value , 3 )    
    """
    
    assert isinstance ( n , integer_types ) and 0 <= n,\
           "round_N: invalid `n' %s (must be non-negative integer)" % n 

    v = float ( value )
    
    if 0 == v : return 0

    a , b = frexp10 ( v ) 
    
    e = b - 1
    m = a * 10
    
    ni = n - 1

    f1 = 10 ** ni 

    f2 = 1 
    if   ni < e : f2 =     ( 10 ** ( e  - ni ) ) 
    elif ni > e : f2 = 1.0/( 10 ** ( ni - e  ) )

    return round ( m * f1 ) * float ( f2 ) 

# =============================================================================
## get ``mantissa'' (1<=m<10) and exponent for radix10
#  similar for frexp, but use radix=10
#  @code
#  m,e = _frexp10_ ( value ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20
def _frexp10_ ( value ) :
    """Get ``mantissa'' (1<=m<10) and exponent for radix10
    similar for frexp, but use radix=10
    >>> m , e = _frexp10_ ( value ) 
    """
    m , e = frexp10 ( value ) 
    return m * 10 , e-1

# =============================================================================
## get three significant digits from the floating value
#  @code
#  nums = three_digits ( value ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20
def three_digits ( value ) :
    """Get three first significant digits
    
    >>> nums = three_digits ( value ) 
    """
    #
    if not 0.1<= abs ( value ) < 1 : value  = frexp10 ( float ( value ) ) [0]
    #
    return int ( round ( float ( value ) * 1000 , 0 ) )


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
#   -  -2, error is NaN
#   -  -1, error is Inf
#   -   0, error is zero 
#   -   1, the first regular case
#   -   2, the second regular case
#   -   3, the third  regular case

#  @code
#  case , rounded_error = pdg_case ( error ) 
#  @endcode 
def pdg_case ( error ) :
    """Classify according to PDG prescription
    - see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    - see section 5.3 of doi:10.1088/0954-3899/33/1/001
    
    The basic rule states that
    - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
    - If they lie between 355 and 949, we round to one significant digit.
    - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    
    Cases:
    -  -2, error is NaN
    -  -1, error is Inf
    -   0, error is zero 
    -   1, the first regular case
    -   2, the second regular case
    -   3, the third  regular case
    
    >>> error = ...
    >>> case , rounded_error = pdg_case ( error ) 
    """
    
    if    math.isnan ( error ) : return -2 , error
    elif  math.isinf ( error ) : return -1 , error
    
    ne = abs ( three_digits ( error ) )
    
    if   0 == ne    : return  0 , 0 
    elif ne  <= 354 : return  1 , round_N ( error , 2 )    
    elif ne  <= 949 : return  2 , round_N ( error , 1 )    
    else            : return  3 , round_N ( error , 1 )    
    
    
# ==============================================================================
# make a rounding according to PDG prescription
#  @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
#  @see section 5.3 of doi:10.1088/0954-3899/33/1/001
#  
#  The basic rule states that if the three highest order digits of the error
#  lie between 100 and 354, we round to two significant digits. If they lie between
#  355 and 949, we round to one significant digit. Finally,
#   if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
#  In all cases, the central value is given with a precision that matches that of the error.
#
#  @code
#  val, err , exponnet, err_case  = pdg_round__ ( value , error )
#  print ( ' Rounded value +/- error is  (%s +/- %s)' % ( val , err ) )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20 
def pdg_round__ ( value , error ) :
    """Make a rounding according to PDG prescription
    - see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    - see section 5.3 of doi:10.1088/0954-3899/33/1/001
    
    Quote: 
    The basic rule states that if the three highest order digits of the error
    lie between 100 and 354, we round to two significant digits. If they lie between
    355 and 949, we round to one significant digit. Finally,
    if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    In all cases, the central value is given with a precision that matches that of the error.
    
    >>> val , err , exponent , err_case = pdg_round__ ( value , error )    
    >>> print ( ' Rounded value +/- error is  (%s +/- %s)' % ( val , err ) ) 
    """

    ecase , err = pdg_case ( error )

    assert -2 <= ecase <= 3 ,\
           'pdg_round: invalid error case %s/%s' %  ( ecase , error ) 

    ## irregular casses :
    if ecase <= 0 or not isfinite ( value ) : 
        return value , error , 0 , ecase

    ## regular error 
    ee , be = _frexp10_ ( error )
    
    if   1 == ecase : q = be + 1 - 2
    elif 2 == ecase : q = be + 1 - 1
    elif 3 == ecase : q = be + 1 - 2

    r   = 10 ** q 
    val = round ( float ( value ) / r ) * r * 1.0 

    return val , err , q , ecase

# ===============================================================================
## make a rounding according to PDG prescription
#  @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
#  @see section 5.3 of doi:10.1088/0954-3899/33/1/001
#  
#  The basic rule states that if the three highest order digits of the error
#  lie between 100 and 354, we round to two significant digits. If they lie between
#  355 and 949, we round to one significant digit. Finally,
#   if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
#  In all cases, the central value is given with a precision that matches that of the error.
#
#  @code
#  val , err , exponent  = pdg_round_ ( value , error )
#  print ( ' Rounded value +/- error is  (%s +/- %s)' % ( val , err ) ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20 
def pdg_round_ ( value , error ) :
    """Make a rounding according to PDG prescription
    - see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    - see section 5.3 of doi:10.1088/0954-3899/33/1/001
    
    Quote: 
    The basic rule states that if the three highest order digits of the error
    lie between 100 and 354, we round to two significant digits. If they lie between
    355 and 949, we round to one significant digit. Finally,
    if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    In all cases, the central value is given with a precision that matches that of the error.
    
    >>> val, err, exponent = pdg_round_ ( value , error )    
    >>> print ( ' Rounded value +/- error is  (%s +/- %s)' % ( val , err ) )
    """
    ##
    val, err , q , c = pdg_round__ ( value , error )
    ##
    return v , e , q  

# =============================================================================
## make a rounding according to PDG prescription
#  @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
#  @see section 5.3 of doi:10.1088/0954-3899/33/1/001
#  
#  The basic rule states that
#   - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
#   - If they lie between 355 and 949, we round to one significant digit.
#   - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
#  In all cases, the central value is given with a precision that matches that of the error.
#
#  @code
#  val , err = pdg_round ( value , error )
#  print( ' Rounded value +/- error is  (%s +/- %s)' % ( val , err ) ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20 
def pdg_round ( value , error ) :
    """Make a rounding according to PDG prescription
    - see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    - see section 5.3 of doi:10.1088/0954-3899/33/1/001
    
    Quote:
    
    The basic rule states that
    - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
    - If they lie between 355 and 949, we round to one significant digit.
    - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    In all cases, the central value is given with a precision that matches that of the error.
    
    >>> val, err = pdg_round ( value , error )    
    >>> print ( ' Rounded value +/- error is  (%s +/- %s)' % ( val , err ) )
    """
    ##
    val , err , f = pdg_round_ ( value , error )
    ##
    return val , err

# =============================================================================
## Round value/error accoriding to PDG prescription and format it for print
#  @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
#  @see section 5.3 of doi:10.1088/0954-3899/33/1/001
#  
#  The basic rule states that
#   - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
#   - If they lie between 355 and 949, we round to one significant digit.
#   - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
# 
#  In all cases, the central value is given with a precision that matches that of the error.
#
#  @code
#  print ' Rounded value/error is %s ' % pdg_format ( value , error , True ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20 
def pdg_format ( value , error , latex = False ) :
    """Round value/error accoridng to PDG prescription and format it for print
    - see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    - see section 5.3 of doi:10.1088/0954-3899/33/1/001
      
    Quote:
    
    The basic rule states that
    - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
    - If they lie between 355 and 949, we round to one significant digit.
    - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    In all cases, the central value is given with a precision that matches that of the error.
    
    >>> value, error = ...
    >>> print ' Rounded value/error is %s ' % pdg_format ( value , error , True ) 
    """

    val , err , q , ecase = pdg_round__ ( value , error )  

    if ecase <= 0 :
        if not isfinite ( val ) :
            return ( '%+g \\pm %-g ' % ( val , err ) ) if latex else ( '% +g +/- %-g ' % ( val , err ) ) 
        else :
            qv , bv = _frexp10_ ( val )
            if 0 != bv : 
                if latex : return '(%+.2f \\pm %-s)\\times 10^{%d}' % ( qv , err / 10**bv , bv )
                else     : return ' %+.2f +/- %-s)*10^{%d} '        % ( qv , err / 10**bv , bv )
            else :
                if latex : return ' %+.2f \\pm %-s ' % ( qv , err )
                else     : return ' %+.2f +/- %-s '  % ( qv , err )
                
    qe , be = _frexp10_ ( error )

    a , b = divmod ( be , 3 ) 

    if   1 == ecase :

        if   0 == b :
            nd = 1
        elif 1 == b :
            nd  = 3
            a  += 1 
        elif 2 == b :
            a  += 1 
            nd = 2

    elif 2 == ecase :

        if   0 == b :
            nd = 0
            
            if 2 == a % 3 :
                nd = 3
                a  = a + 1 
                
        elif 1 == b :            
            nd = 2
            a += 1
            
        elif 2 == b :
            nd = 1
            a += 1

    elif 3 == ecase :
        
        if   0 == b :
            nd = 0
            
            if 2 == a % 3 :
                nd = 3 
                a  = a + 1
                
        elif 1 == b :
            nd = 2
            a += 1 
        elif 2 == b :
            nd = 1
            a += 1 

    if 0 == a :
        
        if latex: fmt = '(%%+.%df \\pm %%.%df)' %  ( nd , nd ) 
        else    : fmt = ' %%+.%df +/- %%.%df '  %  ( nd , nd )

        return fmt % ( val , err )

        
    if latex: fmt = '(%%+.%df \\pm %%.%df)\\times 10^{%%d}' %  ( nd , nd ) 
    else    : fmt = '(%%+.%df +/- %%.%df)*10^{%%d}'         %  ( nd , nd ) 
        
        
    scale = 1.0/10**(3*a)

    return fmt % ( val * scale , err * scale , 3 * a )

# =============================================================================
## Round value/error according to PDG prescription and format it for print
#  @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
#  @see section 5.3 of doi:10.1088/0954-3899/33/1/001
#  
#  The basic rule states that
#   - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
#   - If they lie between 355 and 949, we round to one significant digit.
#   - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
# 
#  In all cases, the central value is given with a precision that matches that of the error.
#
#  @code
#  print ' Rounded value/error is %s ' % pdg_format2 ( value , error , error2 , True ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20 
def pdg_format2( value , error1 , error2  , latex = False , mode = 'min' ) :
    """Round value/error accoridng to PDG prescription and format it for print
    @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    @see section 5.3 of doi:10.1088/0954-3899/33/1/001
      
    Quote:
    
    The basic rule states that
    - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
    - If they lie between 355 and 949, we round to one significant digit.
    - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    In all cases, the central value is given with a precision that matches that of the error.
    
    >>> value, error1, error2 = ...
    >>> print ' Rounded value/error is %s ' % pdg_format2 ( value , error1 , error2 , True ) 
    """

    error = ref_error ( mode , error1 , error2 ) 

    val , err , q , ecase = pdg_round__ ( value , error )  

    if ecase <= 0 or ( not isfinite ( error1 ) ) or ( not isfinite ( error2 ) ) :
        
        if not isfinite ( val ) : 
            return ( '%+g \\pm %-g \\pm %-g ' % ( val , error1 , error2 ) ) if latex else \
                   ( '%+g +/- %-g +/- %-g'    % ( val , error1 , error2 ) ) 
        
        else :
            qv , bv = _frexp10_ ( val )
            if 0 != bv :
                scale = 1.0 / 10**bv 
                if latex : return '(%+.2f \\pm %-s \\pm %-s)\\times 10^{%d}' % ( qv , error1 * scale , error2 * scale , bv )
                else     : return ' %+.2f +/- %-s +/ %-s)*10^{%d} '          % ( qv , error1 * scale , error2 * scale , bv )
            else : 
                if latex : return ' %+.2f \\pm %-s \\pm %-s ' % ( qv , error1 , error2 )
                else     : return ' %+.2f +/- %-s +/- %-s '   % ( qv , error1 , error2  )
                
    qe , be = _frexp10_ ( error )                
    a , b = divmod ( be , 3 ) 
        
    
    if   1 == ecase :

        err1 = round_N ( error1 , 2 ) ## if isclose ( error1 , error , rel_tol = 1.e-2 ) else err 
        err2 = round_N ( error2 , 2 ) ## if isclose ( error2 , error , rel_tol = 1.e-2 ) else err 
        
        if   0 == b :
            nd = 1
        elif 1 == b :
            nd  = 3
            a  += 1 
        elif 2 == b :
            a  += 1 
            nd = 2

    elif 2 == ecase :

        err1 = round_N ( error1 , 1 ) ## if isclose ( error1 , error , rel_tol = 1.e-2 ) else err 
        err2 = round_N ( error2 , 1 ) ## if isclose ( error2 , error , rel_tol = 1.e-2 ) else err 

        if   0 == b :
            nd = 0
            
            if 2 == a % 3 :
                nd = 3
                a  = a + 1 
                
        elif 1 == b :            
            nd = 2
            a += 1
            
        elif 2 == b :
            nd = 1
            a += 1

    elif 3 == ecase :
        
        err1 = round_N ( error1 , 2 ) ## if isclose ( error1 , error , reL_tol = 1.e-2 ) else err  
        err2 = round_N ( error2 , 2 ) ## if isclose ( error2 , error , rel_tol = 1.e-2 ) else err  

        if   0 == b :
            nd = 0
            
            if 2 == a % 3 :
                nd = 3 
                a  = a + 1
                
        elif 1 == b :
            nd = 2
            a += 1 
        elif 2 == b :
            nd = 1
            a += 1 

    if 0 == a :
        
        if latex: fmt = '(%%+.%df \\pm %%.%df \\pm %%.%df )' %  ( nd , nd , nd ) 
        else    : fmt = ' %%+.%df +/- %%.%df +/- %%.%df '    %  ( nd , nd , nd )

        return fmt % ( val , err1 , err2  )

        
    if latex: fmt = '(%%+.%df \\pm %%.%df \\pm %%.%df )\\times 10^{%%d}' %  ( nd , nd , nd ) 
    else    : fmt = '(%%+.%df +/- %%.%df +/- %%.%df)*10^{%%d}'           %  ( nd , nd , nd ) 
        
        
    scale = 1.0/10**(3*a)

    return fmt % ( val * scale , err1 * scale , err2 * scale , 3 * a )




# =============================================================================
## Round value/error according to PDG prescription and format it for print
#  @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
#  @see section 5.3 of doi:10.1088/0954-3899/33/1/001
#  
#  The basic rule states that
#   - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
#   - If they lie between 355 and 949, we round to one significant digit.
#   - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
# 
#  In all cases, the central value is given with a precision that matches that of the error.
#
#  @code
#  print ' Rounded value/error is %s ' % pdg_format2 ( value , error , error2 , True ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20 
def pdg_format3( value , error1 , error2 , error3 , latex = False , mode = 'min' ) :
    """Round value/error accoridng to PDG prescription and format it for print
    @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    @see section 5.3 of doi:10.1088/0954-3899/33/1/001
      
    Quote:
    
    The basic rule states that
    - if the three highest order digits of the error lie between 100 and 354, we round to two significant digits.
    - If they lie between 355 and 949, we round to one significant digit.
    - Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    In all cases, the central value is given with a precision that matches that of the error.
    
    >>> value, error1, error2 = ...
    >>> print ' Rounded value/error is %s ' % pdg_format2 ( value , error1 , error2 , True ) 
    """

    error = ref_error ( mode , error1 , error2 , error3 ) 

    val , err , q , ecase = pdg_round__ ( value , error )  
   
    if ecase <= 0 or  ( not isfinite ( error1 ) ) or ( not isfinite ( error2 ) ) or ( not isfinite ( error3 ) ) :
        
        if not isfinite ( val ) : 
            return ( '%+g \\pm %-g \\pm %-g \\pm %-g ' % ( val , error1 , error2 , error3 ) ) if latex else \
                   ( '%+g +/- %-g +/- %-g +/- %-g'     % ( val , error1 , error2 , error3 ) ) 
        
        else :

            qv , bv = _frexp10_ ( val )
            if 0 != bv :
                scale = 1.0 / 10**bv 
                if latex : return '(%+.2f \\pm %-s \\pm %-s)\\times 10^{%d}' % ( qv , error1 * scale , error2 * scale , error3 * scale , bv )
                else     : return ' %+.2f +/- %-s +/ %-s +/- %-s )*10^{%d} ' % ( qv , error1 * scale , error2 * scale , error3 * scale , bv )
            else :  
                if latex : return ' %+.2f \\pm %-s \\pm %-s \\pm %-s ' % ( qv , error1 , error2 , error3 )
                else     : return ' %+.2f +/- %-s +/- %-s +/- %-s '    % ( qv , error1 , error2 , error3 )



    qe , be = _frexp10_ ( error )                
    a , b = divmod ( be , 3 ) 
        

    if   1 == ecase :

        err1 = round_N ( error1 , 2 ) ## if isclose ( error1 , error , rel_tol = 1.e-2 ) else err  
        err2 = round_N ( error2 , 2 ) ## if isclose ( error2 , error , rel_tol = 1.e-2 ) else err  
        err3 = round_N ( error3 , 2 ) ## if isclose ( error3 , error , rel_tol = 1.e-2 ) else err  
        
        if   0 == b :
            nd = 1
        elif 1 == b :
            nd  = 3
            a  += 1 
        elif 2 == b :
            a  += 1 
            nd = 2

    elif 2 == ecase :

        err1 = round_N ( error1 , 1 ) ## if isclose ( error1 , error , rel_tol = 1.e-2 ) else err  
        err2 = round_N ( error2 , 1 ) ## if isclose ( error2 , error , rel_tol = 1.e-2 ) else err  
        err3 = round_N ( error3 , 1 ) ## if isclose ( error3 , error , rel_tol = 1.e-2 ) else err  

        if   0 == b :
            nd = 0
            
            if 2 == a % 3 :
                nd = 3
                a  = a + 1 
                
        elif 1 == b :            
            nd = 2
            a += 1
            
        elif 2 == b :
            nd = 1
            a += 1

    elif 3 == ecase :
        
        err1 = round_N ( error1 , 2 ) ## if isclose ( error1 , error , rel_tol = 1.e-2 ) else err  
        err2 = round_N ( error2 , 2 ) ## if isclose ( error2 , error , rel_tol = 1.e-2 ) else err  
        err3 = round_N ( error3 , 2 ) ## if isclose ( error3 , error , rel_tol = 1.e-2 ) else err  

        if   0 == b :
            nd = 0
            
            if 2 == a % 3 :
                nd = 3 
                a  = a + 1
                
        elif 1 == b :
            nd = 2
            a += 1 
        elif 2 == b :
            nd = 1
            a += 1 

    if 0 == a :
        
        if latex: fmt = '(%%+.%df \\pm %%.%df \\pm %%.%df \\pm %%.%df)' %  ( nd , nd , nd , nd ) 
        else    : fmt = ' %%+.%df +/- %%.%df +/- %%.%df +/- %%.%df '    %  ( nd , nd , nd . nd )

        return fmt % ( val , err )

        
    if latex: fmt = '(%%+.%df \\pm %%.%df \\pm %%.%df \\pm %%.%df)\\times 10^{%%d}' %  ( nd , nd , nd , nd ) 
    else    : fmt = '(%%+.%df +/- %%.%df +/- %%.%df +/- %%.%df)*10^{%%d}'          %  ( nd , nd , nd , nd ) 
        
        
    scale = 1.0/10**(3*a)

    return fmt % ( val * scale , err1 * scale , err2 * scale , err3 * scale , 3 * a )

# ====================================================================================
## make a rounding according to PDG prescription
#  @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
#  @see section 5.3 of doi:10.1088/0954-3899/33/1/001
#  
#  The basic rule states that if the three highest order digits of the error
#  lie between 100 and 354, we round to two significant digits. If they lie between
#  355 and 949, we round to one significant digit. Finally,
#   if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
#  In all cases, the central value is given with a precision that matches that of the error.
#
#  @code
#  ve  = VE( ...
#  vr  = ve.pdg()
#  print ' Rounded value with error is  %s ' % vr 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20 
def _ve_pdg_ ( ve ) :
    """Make a rounding according to PDG prescription
    @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    @see section 5.3 of doi:10.1088/0954-3899/33/1/001
      
    Quote:
    The basic rule states that if the three highest order digits of the error
    lie between 100 and 354, we round to two significant digits. If they lie between
    355 and 949, we round to one significant digit. Finally,
    if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    In all cases, the central value is given with a precision that matches that of the error.

    >>> ve  = VE( ...
    >>> vr  = ve.pdg()
    >>> print ' Rounded value with error is  %s ' % vr
    """
    #
    v , e  = pdg_round  ( ve.value() , ve.error() )
    # 
    return VE ( v , e * e )

# =============================================================================
## Round value/error accoridng to PDG prescription and format it for print
#  @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
#  @see section 5.3 of doi:10.1088/0954-3899/33/1/001
#  
#  The basic rule states that if the three highest order digits of the error
#  lie between 100 and 354, we round to two significant digits. If they lie between
#  355 and 949, we round to one significant digit. Finally,
#   if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
#  In all cases, the central value is given with a precision that matches that of the error.
#
#  @code
#  ve = VE( ... ) 
#  print ' Rounded value/error is %s ' % ve.pdg_format ()
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20 
def _ve_pdg_format_ ( ve , latex = False ) :
    """Round value/error accoridng to PDG prescription and format it for print
    @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    @see section 5.3 of doi:10.1088/0954-3899/33/1/001
      
    Quote:
    The basic rule states that if the three highest order digits of the error
    lie between 100 and 354, we round to two significant digits. If they lie between
    355 and 949, we round to one significant digit. Finally,
    if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    In all cases, the central value is given with a precision that matches that of the error.
    
    >>> ve = VE(... ) 
    >>> print ' Rounded value/error is %s ' % ve.pdg_format ()
    """
    return pdg_format ( ve.value() , ve.error() , latex ) 


# =============================================================================
## finally decorate class ValueWith Error
# =============================================================================

VE.pdg        = _ve_pdg_
VE.pdg_format = _ve_pdg_format_ 
    
# =============================================================================
## insert it to math
if not hasattr ( math , 'frexp10' ) : math.frexp10 = frexp10
if not hasattr ( math , 'round_N' ) : math.round_N = round_N


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
