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
# Version           $Revision$
# Last modification $Date$
#                by $Author$
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
    'frexp10'     , ## similar to math.frexp but with radix=10
    'round_N'     , ## round floaing value to N-significant digits
    'pdg_round'   , ## round value,error-pair according to PDG prescription
    'pdg_format'  , ## format value&error according to PDF prescription
    'pdg_format2' , ## format value+2errorr according to PDF
    'pdg_format3' , ## format value+3errors according to PDF
    )
# ===============================================================================
import ROOT,         math
from   ostap.math.ve import cpp, VE 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.utils.pdg_format' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================

# =============================================================================
## round to nearest integer, rounds half integers to nearest even integer 
#  It is just a simple wrapper around boost::numeric::converter 
#  @see Ostap::Math::round 
cpp_round   = cpp.Ostap.Math.round
cpp_frexp10 = cpp.Ostap.Math.frexp10 
cpp_round_N = cpp.Ostap.Math.round_N 

# =============================================================================
## get mantissa (0.1<=m<1) and exponent for radix10
#  similar for frexp, but use radix=10
#  @code
#  m,e = frexp10 ( value ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20
def frexp10 ( value ) :
    """
    Get the mantissa (0.1<=m<1) and exponent for radix10
    (similar for frexp, but use radix=10)
    
    >>> a,b = frexp10 ( value ) 
    """
    #
    p = cpp_frexp10 ( value )
    return p.first, p.second

# =============================================================================
## get ``mantissa'' (1<=m<10) and exponent for radix10
#  similar for frexp, but use radix=10
#  @code
#  m,e = _frexp10_ ( value ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20
def _frexp10_ ( value ) :
    """
    Get ``mantissa'' (1<=m<10) and exponent for radix10
    similar for frexp, but use radix=10
    >>> m,e = _frexp10_ ( value ) 
    """
    m,e = frexp10 ( value ) 
    return m*10,e-1

# =============================================================================
## get three significant digits from the floating value
#  @code
#  nums = _3digits_ ( value ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20
def _3digits_ ( value ) :
    """
    Get three first significant digits
    
    >>> nums = _3digits_ ( value ) 
    """
    #
    if not 0.1<= abs ( value ) < 1 : value  = frexp10 ( value ) [0]
    #
    return int ( round ( float ( value ) * 1000 , 0 ) )

## ============================================================================
#  round value to N-digits
#  @code
#  new_value = round_N ( value , 3 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20
def round_N ( value , n ) :
    """
    Round value to to N-digits
    
    >>> new_value = round_N ( value , 3 )    
    """
    return cpp_round_N ( value , n ) 

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
#  v,e,f = pdg_round_ ( value , error )
#  print ' Rounded value +/- error is  (%s +/- %s)' % ( v , e )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20 
def pdg_round_ ( value , error ) :
    """
    Make a rounding according to PDG prescription
    @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    @see section 5.3 of doi:10.1088/0954-3899/33/1/001
    
    Quote: 
    The basic rule states that if the three highest order digits of the error
    lie between 100 and 354, we round to two significant digits. If they lie between
    355 and 949, we round to one significant digit. Finally,
    if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    In all cases, the central value is given with a precision that matches that of the error.
    
    >>> v,e,f = pdg_round_ ( value , error )    
    >>> print ' Rounded value +/- error is  (%s +/- %s)' % ( v , e )
    """
    ##
    ne = _3digits_ ( error )
    ##
    if   abs ( ne ) <= 354 : pr,pe = 2,2
    elif abs ( ne ) <= 949 : pr,pe = 1,1
    else                   : pr,pe = 2,1
    ##
    ev , bv = _frexp10_ ( value )
    ee , be = _frexp10_ ( error )
    ##
    e = round_N ( error , pe )
    ##
    q  = be + 1 - pe 
    r  = 10**q
    v  = cpp_round ( value / r ) * r
    ##
    return v , e , q 

# =============================================================================
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
#  v,e = pdg_round ( value , error )
#  print ' Rounded value +/- error is  (%s +/- %s)' % ( v , e )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20 
def pdg_round ( value , error ) :
    """
    Make a rounding according to PDG prescription
    @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    @see section 5.3 of doi:10.1088/0954-3899/33/1/001
    
    Quote: 
    The basic rule states that if the three highest order digits of the error
    lie between 100 and 354, we round to two significant digits. If they lie between
    355 and 949, we round to one significant digit. Finally,
    if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    In all cases, the central value is given with a precision that matches that of the error.
    
    >>> v,e = pdg_round ( value , error )    
    >>> print ' Rounded value +/- error is  (%s +/- %s)' % ( v , e )
    """
    ##
    v,e,f = pdg_round_ ( value , error )
    ##
    return v, e 

# =============================================================================
## Round value/error accoriding to PDG prescription and format it for print
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
#  print ' Rounded value/error is %s ' % pdg_format ( value , error , True ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20 
def pdg_format ( value , error , latex = False ) :
    """
    Round value/error accoridng to PDG prescription and format it for print
    @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    @see section 5.3 of doi:10.1088/0954-3899/33/1/001
      
    Quote:
    The basic rule states that if the three highest order digits of the error
    lie between 100 and 354, we round to two significant digits. If they lie between
    355 and 949, we round to one significant digit. Finally,
    if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    In all cases, the central value is given with a precision that matches that of the error.
    
    >>> value, error = ...
    >>> print ' Rounded value/error is %s ' % pdg_format ( value , error , True ) 
    """

    ## round value, erorr and get the scaling factor
    v , e , q = pdg_round_  ( value , error )
    
    ## get scaling factors for printout 
    qv , bv   = _frexp10_ (  value  )
    n         = int       ( math.floor ( bv / 3. ) ) 


    q = abs ( q ) 
    short = ( 0 == n ) or (  1 == n and  abs(v) < 10000 ) or ( -1 == n and  abs(v) > 0.1 ) 
    
    if not short :

        q += 3 * n
        
        q  = abs ( q )
        
        if latex : fmt = "$(%%.%df \pm %%.%df)\cdot10^{%%d}$" % ( q , q )
        else     : fmt =  "(%%.%df +/- %%.%df)*10^{%%d}"      % ( q , q )

        v  /= 10**(3*n)
        e  /= 10**(3*n)
        ##
                
        return fmt % ( v , e , n*3 )
        
    else :
        
        q = abs ( q ) 

        if latex : fmt = "$%%.%df \pm %%.%df$" % ( q , q )
        else     : fmt =  "%%.%df +/- %%.%df"  % ( q , q )
        ##

        return fmt % ( v , e )

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
#  print ' Rounded value/error is %s ' % pdg_format2 ( value , error , error2 , True ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20 
def pdg_format2( value , error1 , error2  , latex = False , mode = 'total' ) :
    """
    Round value/error accoridng to PDG prescription and format it for print
    @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    @see section 5.3 of doi:10.1088/0954-3899/33/1/001
      
    Quote:
    The basic rule states that if the three highest order digits of the error
    lie between 100 and 354, we round to two significant digits. If they lie between
    355 and 949, we round to one significant digit. Finally,
    if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    In all cases, the central value is given with a precision that matches that of the error.
    
    >>> value, error1, error2 = ...
    >>> print ' Rounded value/error is %s ' % pdg_format2 ( value , error1 , error2 , True ) 
    """

    if isinstance ( mode , str ) : mode = mode.lower()
    ## 
    if   mode in ( 'min'   , 'mn'  ) : 
        error = min ( abs ( error1 ) , abs ( error2 ) )
    elif mode in ( 'max'   , 'mx'  ) : 
        error = min ( abs ( error1 ) , abs ( error2 ) )
    elif mode in ( 'total' , 'tot' ) : 
        error = math.sqrt ( 1.0 * error1 * error1 + error2 * error2 )
    else :
        ## use the default policy 
        error = math.sqrt ( 1.0 * error1 * error1 + error2 * error2 )

    ## round value, error and get the scaling factor
    v  , e  , q  = pdg_round_  ( value , error )

    v_ , e1 , q_ = pdg_round_  ( value , error1 )
    v_ , e2 , q_ = pdg_round_  ( value , error2 )

    ## get scaling factors for printout 
    qv , bv   = frexp10 (  value  )
    n         = int     ( math.floor ( bv / 3. ) ) 


    q = abs ( q ) 
    short = ( 0 == n ) or (  1 == n and  abs(v) < 10000 ) or ( -1 == n and  abs(v) > 0.1 ) 
    
    if not short :

        q += 3 * n
        
        q  = abs ( q )
        
        if latex : fmt = "$(%%.%df \pm %%.%df \pm %%.%df)\cdot10^{%%d}$" % ( q , q , q )
        else     : fmt =  "(%%.%df +/- %%.%df +/- %%.%df)*10^{%%d}"      % ( q , q , q )

        v  /= 10**(3*n)
        e1 /= 10**(3*n)
        e2 /= 10**(3*n)
        ##
                
        return fmt % ( v , e1 , e2 , n*3 )
        
    else :
        
        q = abs ( q ) 

        if latex : fmt = "$%%.%df \pm %%.%df \pm %%.%df$" % ( q , q , q )
        else     : fmt =  "%%.%df +/- %%.%df +/- %%.%df"  % ( q , q , q )
        ##

        return fmt % ( v , e1 , e2 )


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
#  print ' Rounded value/error is %s ' % pdg_format2 ( value , error , error2 , True ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20 
def pdg_format3( value , error1 , error2 , error3 , latex = False , mode = 'total' ) :
    """
    Round value/error accoridng to PDG prescription and format it for print
    @see http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf
    @see section 5.3 of doi:10.1088/0954-3899/33/1/001
      
    Quote:
    The basic rule states that if the three highest order digits of the error
    lie between 100 and 354, we round to two significant digits. If they lie between
    355 and 949, we round to one significant digit. Finally,
    if they lie between 950 and 999, we round up to 1000 and keep two significant digits.
    In all cases, the central value is given with a precision that matches that of the error.
    
    >>> value, error1, error2 = ...
    >>> print ' Rounded value/error is %s ' % pdg_format2 ( value , error1 , error2 , True ) 
    """

    if isinstance ( mode , str ) : mode = mode.lower()
    ## 
    if   mode in ( 'min'   , 'mn'  ) : 
        error = min ( abs ( error1 ) , abs ( error2 ) , abs ( error3 ) )
    elif mode in ( 'max'   , 'mx'  ) : 
        error = min ( abs ( error1 ) , abs ( error2 ) , abs ( error3 ) )
    elif mode in ( 'total' , 'tot' ) : 
        error = math.sqrt ( 1.0 * error1 * error1 + error2 * error2 + error3 * error3 )
    else :
        ## use the default policy 
        error = math.sqrt ( 1.0 * error1 * error1 + error2 * error2 + error3 * error3 )

    ## round value, error and get the scaling factor
    v  , e  , q  = pdg_round_  ( value , error )

    v_ , e1 , q_ = pdg_round_  ( value , error1 )
    v_ , e2 , q_ = pdg_round_  ( value , error2 )
    v_ , e3 , q_ = pdg_round_  ( value , error3 )

    ## get scaling factors for printout 
    qv , bv   = frexp10 (  value  )
    n         = int     ( math.floor ( bv / 3. ) ) 


    q = abs ( q ) 
    short = ( 0 == n ) or (  1 == n and  abs(v) < 10000 ) or ( -1 == n and  abs(v) > 0.1 ) 
    
    if not short :

        q += 3 * n
        
        q  = abs ( q )
        
        if latex : fmt = "$(%%.%df \pm %%.%df \pm %%.%df \pm %%.%df)\cdot10^{%%d}$" % ( q , q , q , q )
        else     : fmt =  "(%%.%df +/- %%.%df +/- %%.%df +/- %%.%df)*10^{%%d}"      % ( q , q , q , q )

        v  /= 10**(3*n)
        e1 /= 10**(3*n)
        e2 /= 10**(3*n)
        e3 /= 10**(3*n)
        ##
                
        return fmt % ( v , e1 , e2 , e3 , n*3 )
        
    else :
        
        q = abs ( q ) 

        if latex : fmt = "$%%.%df \pm %%.%df \pm %%.%df \pm %%.%df$" % ( q , q , q , q )
        else     : fmt =  "%%.%df +/- %%.%df +/- %%.%df +/- %%.%df"  % ( q , q , q , q )
        ##

        return fmt % ( v , e1 , e2 , e3 )

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
    """
    Make a rounding according to PDG prescription
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
    """
    Round value/error accoridng to PDG prescription and format it for print
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

    from ostap.logger.line import line 
    logger.info ( __file__  + '\n' + line  ) 
    logger.info ( 80*'*'   )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 
  
    v = 1./3 
    for s in ( 1.e-6 , 1.e-3 , 1.e3 , 1.e6 ) : 
        for e in ( 0.000123 , 0.000456 , 0.000900 , 0.000986 ) :
            ve = VE ( 1./3 , e*e ) * s 
            logger.info ( ' Value %s,\tround: %s,\tTeX: %s' % ( ve , ve.pdg() , ve.pdg_format( True ) ) )

    logger.info ( 80*'*' ) 
    for s in ( 1 , 10000 , 0.01 ) :
        
        logger.info ( pdg_format  ( 1.0*s/3 , 0.000012345 ) ) 
        logger.info ( pdg_format  ( 1.0*s/3 , 0.0020      ) ) 
        logger.info ( pdg_format  ( 1.0*s/3 , 0.0050      ) ) 
        logger.info ( pdg_format  ( 1.0*s/3 , 0.0099      ) ) 
        
        logger.info ( pdg_format2 ( 1.0*s/3 , 0.000012345 , 0.00500      ) ) 
        logger.info ( pdg_format2 ( 1.0*s/3 , 0.0020      , 0.0040       ) ) 
        logger.info ( pdg_format2 ( 1.0*s/3 , 0.0050      , 0.000012345  ) ) 
        logger.info ( pdg_format2 ( 1.0*s/3 , 0.0099      , 0.0020       ) )
        
        logger.info ( pdg_format3 ( 1.0*s/3 , 0.000012345 , 0.00500      , 0.0001 ) ) 
        logger.info ( pdg_format3 ( 1.0*s/3 , 0.0020      , 0.0040       , 0.0001 ) )
        logger.info ( pdg_format3 ( 1.0*s/3 , 0.0050      , 0.000012345  , 0.0001 ) ) 
        logger.info ( pdg_format3 ( 1.0*s/3 , 0.0099      , 0.0020       , 0.0001 ) ) 
        
    logger.info ( 80*'*' ) 
# =============================================================================
# The END 
# =============================================================================
