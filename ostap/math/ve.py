#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# =============================================================================
# @file ostap/math/ve.py
# Simple ``value-with-error'' concept
# @code
# a = VE (1,0.1**2)
# b = VE (2,0.2**2)
# c = a + b 
# @endcode
# @see Ostap::Math::ValueWithError
# ============================================================================= 
r""" Simple `value-with-error' concept

>>> a = VE (1,0.1**2) 
>>> b = VE (2,0.2**2)
>>> c = a + b

"""
# ============================================================================= 
__all__     = (
    'VE'  ,  # Value with error  
    'VVE' ,  # vector of values with errors
    )
# ============================================================================= 
from   ostap.math.base       import Ostap, std , iszero, isequal
from   ostap.math.random_ext import poisson as _poisson 
import ROOT
# ============================================================================= 
# logging 
# ============================================================================= 
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.ve' )
else                       : logger = getLogger ( __name__        )
# ============================================================================= 

VE        = Ostap.Math.ValueWithError
VVE       = std.vector ( VE )

VE.Vector = VVE
VE.Vector .__str__   = lambda s : str( [ i for i in s ])
VE.Vector .__repr__  = lambda s : str( [ i for i in s ])

VVVE = std.vector( VVE )
VVE.Vector = VVVE 
VVE.Vector . __str__  = lambda s : str( [ i for i in s ] )
VVE.Vector . __repr__ = lambda s : str( [ i for i in s ] )
VVE.Vector . __len__  = lambda s : s.size ()

# ============================================================================= 
## Sum the contents of the vector
def _ve_sum_ ( s ) :
    """Sum the contents of the vector.
    >>> v = ...
    >>> s = v.sum()
    """
    return Ostap.Math.sum ( s )

# ============================================================================= 
## Sum the contents of the vector
def _ve_asum_ ( s ) :
    """Sum the contents of the vector.
    >>> v = ...
    >>> s = v.abssum()
    """
    return Ostap.Math.abssum ( s )

_ve_sum_  . __doc__ += '\n' + Ostap.Math.sum    .__doc__
_ve_asum_ . __doc__ += '\n' + Ostap.Math.abssum .__doc__


for t in ( Ostap.Math.ValueWithError         ,
           Ostap.Math.Point3DWithError       ,
           Ostap.Math.Vector3DWithError      ,
           Ostap.Math.LorentzVectorWithError ) :
    
    if not hasattr ( t , '_new_str_' ) :
        t._new_str_ = t.toString
        t.__str__   = t.toString
        t.__repr__  = t.toString

from ostap.logger.symbols import plus_minus
ve_fmt = '%%+.4g %s %%-.4g' % plus_minus 
def _ve_str_ ( ve ) : return ve.toString ( ve_fmt )  
Ostap.Math.ValueWithError .__str__  = _ve_str_
Ostap.Math.ValueWithError .__repr__ = _ve_str_

# ============================================================================= 
## Calculate the "effective" background-to-signal ratio from the valeu 
#  and its uncertainty using the identity
#  \f$ \frac{\sigma(S)}{S} = \frac{\sqrt{S}}{S}\sqrt{1+\frac{B}{S}}\f$.
#  From this identity one gets
#  \f$ \left.\frac{B}{S}\right|_{\mathrm{eff}} \equiv \frac{\sigma^2(S)}{S} -1 \f$
#  @param v the value 
#  @return the effective backround-to-signal ratio or -1 
#  @code
#  v = VE( ... )
#  print 'B/S=', v.b2s() 
#  @endcode
#  @see Ostap::Math::b2s 
def _ve_b2s_ ( s )  :
    """ Get background-over-signal ratio B/S estimate from the equation:
    error(S) = 1/sqrt(S) * sqrt ( 1 + B/S).
    >>> v = ...
    >>> b2s = v.b2s() ## get B/S estimate
    - see `Ostap.Math/b2s`
    """
    return Ostap.Math.b2s ( VE ( s ) ) 

# =============================================================================
## Calculate the "effective purity" ratio using the identity
#  \f$ p_{\mathrm{eff}} = \frac{S}{S+B} = \frac{1}{1+\frac{B}{S}}\f$
#  and the effective "background-to-signal" ratio is estimated as 
#  \f$ \left.\frac{B}{S}\right|_{\mathrm{eff}} = \frac{\sigma^2(S)}{S} -1 \f$, 
#  finally one gets 
#  \f$ p_{\mathrm{eff}} \equiv \frac{S}{\sigma^2(S)}\f$
#  @see Ostap::Math::b2s 
#  @see Ostap::Math::purity 
#  @param v the value 
#  @return the effective purity or -1
#  @see Ostap.Math.purity 
def _ve_purity_ ( s ) :
    """ Calculate the ``effective purity'' ratio using the identity
    p = S/(S+B) = 1/( 1 + B/S ), 
    - and the effective ``background-to-signal'' ratio B/S is estimated as
    B/S = sigma^2(S)/S - 1
    - Finally one gets
    p = S / sigma^2(S) 
    - see `Ostap.Math.b2s` 
    - see `Ostap.Math.purity`
    """
    return Ostap.Math.purity  ( VE ( s ) ) 

# ============================================================================= 
## Get precision with ``some'' error estimate.
#  @see Ostap::Math::precision 
def _ve_precision_ ( s )  :
    """ Get precision with `some' error estimate.
    >>> v = ...
    >>> p = v.prec()
    - see `Ostap.Math.precision`
    """
    return Ostap.Math.precision ( VE ( s ) )

# =============================================================================
## Calculate the useful Figure-of-merit, aka "significance" :  S/(S+B)
#  @see Ostap::Math::FoM
def _ve_fom_  ( s ) :
    """ Calculate the useful Figure-of-merit, aka "significance" :  S/(S+B)
    - see `Ostap.Math.FoM`
    """
    return Ostap.Math.FoM ( VE ( s ) )

# =============================================================================
## Calculate the useful Figure-of-merit, aka "significance times purity" 
#  @see Ostap::Math::FoM2
def _ve_fom2_  ( s ) :
    """ Calculate the useful Figure-of-merit, aka "significance times purity" 
    - see `Ostap.Math.FoM2`
    """
    return Ostap.Math.FoM2 ( VE ( s ) )

VE . b2s        = _ve_b2s_
VE . B2S        = _ve_b2s_
VE . purity     = _ve_purity_
VE . precision  = _ve_precision_
VE . prec       = _ve_precision_
VE . FoM        = _ve_fom_
VE . FoM1       = _ve_fom_
VE . FoM2       = _ve_fom2_

_is_le_    = Ostap.Math.LessOrEqual ( 'double' )()

# =============================================================================
## Comparison of ValueWithError object with other objects
#  @attention it is comparison by value only, errors are ignored 
def _ve_lt_ ( self , other ) :
    """ Comparison of ValueWithError object with other objects
    >>> a = VE( ... )
    >>> print a < b 
    Attention: comparison by value only!
    """
    return float(self) < float(other)

# =============================================================================
## Comparison of ValueWithError object with other objects
#  @attention it is comparison by value only, errors are ignored 
def _ve_le_ ( self , other ) :
    """ Comparison of ValueWithError object with other objects
    >>> a = VE( ... )
    >>> print a <= b 
    Attention: comparison by value only!
    """
    return _is_le_ ( float(self) , float(other) )

# =============================================================================
## Comparison of ValueWithError object with other objects
#  @attention it is comparison by value only, errors are ignored 
def _ve_gt_ ( self , other ) :
    """ Comparison of ValueWithError object with other objects
    >>> a = VE( ... )
    >>> print a > b 
    Attention: comparison by value only!
    """
    return float(self) > float(other)

# =============================================================================
## Comparison of ValueWithError object with other objects
#  @attention it is comparison by value only, errors are ignored 
def _ve_ge_ ( self , other ) :
    """ Comparison of ValueWithError object with other objects
    >>> a = VE( ... )
    >>> print a >= b 
    Attention: comparison by value only!
    """
    return _is_le_ ( float(other) , float(self) )
    
VE.__lt__ = _ve_lt_
VE.__le__ = _ve_le_ 
VE.__gt__ = _ve_gt_ 
VE.__ge__ = _ve_ge_


_is_equal_ = Ostap.Math.Equal_To    ( 'double' )()
_is_zero_  = Ostap.Math.Zero        ( 'double' )()

if not hasattr ( VE ,  '__truediv__' ) : VE. __truediv__  = VE. __div__ 
if not hasattr ( VE , '__itruediv__' ) : VE.__itruediv__  = VE.__idiv__ 
if not hasattr ( VE , '__rtruediv__' ) : VE.__rtruediv__  = VE.__rdiv__ 

# =============================================================================
## Equality for ValueWithError objects
def _ve_eq_ ( self , other ) :
    """Equality for ValueWithError objects
    >>> a = VE( ... )
    >>> b = VE( ... )
    >>> print a == b 
    """
    if isinstance ( other , VE ) :
        v1 = self .value()
        v2 = other.value()
        return _is_equal_ ( v1 , v2 ) and _is_equal_ ( self.cov2() , other.cov2() )
    elif _is_zero_ ( self.cov2() )  :
        return _is_equal_ ( float ( self ) ,  float ( other )  ) 
    else :
        raise NotImplementedError ( ' Equality for %s and  %s is not implemented' % ( self , other ) ) 

# =============================================================================
## inequality for ValueWithError objects
def _ve_ne_ ( self , other ) :
    """ Inequality for ValueWithError objects
    >>> a = VE( ... )
    >>> b = VE( ... )
    >>> print a != b 
    """
    try: 
        return not self == other
    except NotImplementedError :
        raise NotImplementedError ( ' Inequality for %s and  %s is not implemented' % ( self , other ) ) 

VE . __eq__ = _ve_eq_
VE . __ne__ = _ve_ne_

# ==============================================================================
## get easy (and coherent)  way  to access min/max for
#  the value with error object: (value-n*error,value+n*error)
#  @code
#  ve = VE(2,2)
#  print ve.minmax()
#  print ve.minmax(2)
#  @endcode 
def _ve_minmax_ ( s , n = 1 ) :
    """ Get an easy and coherent way to access `min/max' for
    the value with error object:  (value-n*error,value+n*error)
    >>> ve = VE(2,2) 
    >>> print ve.minmax()
    >>> print ve.minmax(2)
    """
    v = s.value()
    e = s.error() 
    if e <= 0 : return v,v
    v1 = v + e * n
    v2 = v - e * n
    if v1 <= v2 : return v1 , v2
    return v2,v1 

VE.minmax = _ve_minmax_

# =============================================================================
## hashing for VE object
#  @code
#  v = VE   ( ... )
#  h = hash ( v   ) 
#  @endcode
def _ve_hash_ ( v ) :
    """ Hashing function for VE objecte
    >>> v = VE   ( ... )
    >>> h = hash ( v   ) 
    """
    return  hash ( ( v.value() , v.cov2() ) )

VE.__hash__ = _ve_hash_

# =============================================================================
## Get the "positive" error
#  @code
#  ve = ...
#  pos_error = ve.pos_error
#  @endcode 
def _ve_pos_error_ ( v ) :
    """ Get the "positive" error
    >>> ve = ...
    >>> pos_error = ve.pos_error
    """
    err = v.error()
    return err if 0 <= err else 0 

# =============================================================================
## Get the "negative" error
#  @code
#  ve = ...
#  neg_error = ve.neg_error
#  @endcode 
def _ve_neg_error_ ( v ) :
    """ Get the "negative" error
    >>> ve = ...
    >>> neg_error = ve.neg_error
    """
    err = v.error()
    return -err if 0 <= err else 0 

# =============================================================================
## Get "negative&positive" errors
#  @code
#  ve = ...
#  neg, pos = ve.errors () 
#  @endcode  
def _ve_errors_  ( v ) :
    """ Get `negative&positive' errors
    >>> ve = ...
    >>> neg, pos = ve.errors () 
    """
    err = v.error()
    return (-err,err) if 0 <= err else (0,0)

VE.pos_error = property ( _ve_pos_error_ , None , None , _ve_pos_error_ .__doc__ ) 
VE.neg_error = property ( _ve_neg_error_ , None , None , _ve_neg_error_ .__doc__ ) 
VE.errors    = property ( _ve_errors_    , None , None , _ve_errors_    .__doc__ ) 

# =============================================================================
from random import gauss as _gauss     
# =============================================================================
## get the (gaussian) random number according to parameters
#
#  @code
#    >>> v = ...  ## the number with error
#
#    ## get 100 random numbers 
#    >>> for i in range ( 0, 100 ) : print v.gauss()
#    
#    ## get only non-negative numbers 
#    >>> for j in range ( 0, 100 ) : print v.gauss( lambda s : s > 0 )
#
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-08-10
# 
def _ve_gauss_ ( s , accept = lambda a : True , nmax = 1000 ) : 
    """ Get the gaussian random number
    >>> v = ...  ## the number with error
    ## get 100 random numbers 
    >>> for i in range ( 0, 100 ) : print v.gauss()    
    ## get only non-negative numbers 
    >>> for j in range ( 0, 100 ) : print v.gauss( lambda s : s > 0 )
    """
    #
    if 0 >= s.cov2() or iszero ( s.cov2 () ) : return s.value() ## return
    #
    v = s.value ()
    e = s.error ()
    #
    for i in range ( nmax  ) :
        r = _gauss  ( v , e ) 
        if accept   ( r     ) : return r
        
    logger.warning("Can'n generate proper random number %s" % s )
    return v

# =============================================================================
## generate poisson random number according to parameters 
#  @code
#    >>> v = ...  ## the number with error
#
#    ## get 100 random numbers 
#    >>> for i in range ( 0, 100 ) : print v.poisson ( fluctuate = True )
#    
#    ## get only odd numbers 
#    >>> for j in range ( 0, 100 ) : print v.poisson ( fluctuate = True , accept = lambda s : 1 ==s%2 )
#
#    ## do not fluctuate the mean of poisson:    
#    >>> for j in range ( 0, 100 ) : print v.poisson ( fluctuate = False  )
#
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-08-10   
def _ve_poisson_ ( s , fluctuate , accept = lambda s : True ) :
    """ Generate poisson random number according to parameters 
    >>> v = ...  ## the number with error    
    ## get 100 random numbers 
    >>> for i in range ( 0, 100 ) : print v.poisson()    
    ## get only odd numbers 
    >>> for j in range ( 0, 100 ) : print v.poisson ( accept = lambda s : 1 ==s%2 )    
    ## do not fluctuate the mean of poisson:    
    >>> for j in range ( 0, 100 ) : print v.poisson ( fluctuate = False  )
    """
    s = VE( s ) 
    v = s.value() 
    if v < 0 and not fluctuate :
        raise TypeError ( 'Negative mean without fluctuations (1)' )
    if v < 0 and s.cov2() <= 0 :
        raise TypeError ( 'Negative mean without fluctuations (2)' )

    e = s.error() 
    if v < 0 and abs(v) > 3 * e  :
        logger.warning ( "Very inefficient mean fluctuations: %s" % s ) 

    mu = v
    if fluctuate :
        mu = s.gauss ()
        while mu < 0 :
            mu = s.gauss ()
    
    return _poisson ( mu ) 

VE.gauss   = _ve_gauss_
VE.poisson = _ve_poisson_ 





# ===================================================================================
## Formats f0r nice printout of the ValueWithError object  ( string + exponent)
#  @code
#  fmt , fmtv , fmte , expo = fmt_pretty_ve ( number ) 
#  @endcode
#  @return nice string and the separate exponent 
def fmt_pretty_ve ( value              ,
                    width       = 6    ,
                    precision   = 4    ,
                    parentheses = True ) :
    """ Formats f0r nice printout of the ValueWithError object  ( string + exponent)
    >>> fmt , fmtv , fmte , expo = fmt_pretty_ve ( number ) 
    """
    assert isinstance ( value , VE ) , \
        "Invalid `value' parameter: %s" % type ( value )
    ## decode object 
    v , e =  value.value () , max ( 0 , value.error () )
    from ostap.logger.pretty import fmt_pretty_error 
    return fmt_pretty_error ( v , e ,
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
    """ Nice printout of the ValueWithError object  ( string + exponent)
    - return nice stirng and the separate exponent 
    >>> s , expo = pretty_ve ( number ) 
    """
    assert isinstance ( value , VE ) , "Invalid `value' parameter: %s" % type ( value )
    ## decode object 
    v , e =  value.value () , max ( 0 , value.error () )
    ## delegate 
    from ostap.logger.pretty import pretty_error 
    return pretty_error ( v , e , width = width , precision = precision , parentheses = parentheses  )

# ==============================================================================
Ostap.Math.ValueWithError.pretty_print = pretty_ve 

# ==============================================================================
## factory for unpickling of <code>Ostap::Math::ValueWithError</code>
#  @see Ostap::Math::ValueWithError
def ve_factory ( value , cov2 ) :
    """ Factory for unpickling of <code>Ostap::Math::ValueWithError</code>
    - see Ostap::Math::ValueWithError
    """
    return  VE ( value , cov2 ) 
    
# =============================================================================
## reduce <code>Ostap::Math::ValueWithError</code>
#  @see Ostap::Math::ValueWithError
def ve_reduce ( v ) :
    """reduce `Ostap.Math.ValueWithError`
    - see Ostap.Math.ValueWithError
    """
    return ve_factory , ( v.value() , v.cov2() )

Ostap.Math.ValueWithError.__reduce__ = ve_reduce

# =============================================================================
## decorated classes 
_decorated_classes_  = (
    Ostap.Math.ValueWithError         ,
    Ostap.Math.ValueWithError.Vector  ,
    Ostap.Math.Point3DWithError       ,
    Ostap.Math.Vector3DWithError      ,
    Ostap.Math.LorentzVectorWithError )

# =============================================================================



# =============================================================================
## decorated methods 
_new_methods_ = (
    VE.Vector  . __str__  ,
    VE.Vector  . __repr__ , 
    VVE.Vector . __str__  ,
    VVE.Vector . __repr__ , 
    VVE.Vector . __len__  , 
    VE . b2s              , 
    VE . B2S              , 
    VE . purity           ,  
    VE . precision        , 
    VE . prec             , 
    VE . FoM              , 
    VE . FoM1             , 
    VE . FoM2             , 
    VE . __lt__           ,
    VE . __le__           , 
    VE . __gt__           , 
    VE . __ge__           , 
    VE . __eq__           , 
    VE . __ne__           , 
    VE . __hash__         , 
    VE . minmax           ,
    VE . gauss            , 
    VE . poisson          ,
    VE . __reduce__       ,
    VE . pretty_print     ,
   )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme  import docme
    docme ( __name__ , logger = logger )
    
    a = VE(100,100)
    b = VE(400,400)
    
    logger.info ( 'a=%s, b=%s' % ( a , b ) )
    
    logger.info ( 'a+b         %s' % ( a + b ) )
    logger.info ( 'a-b         %s' % ( a - b ) )
    logger.info ( 'a*b         %s' % ( a * b ) )
    logger.info ( 'a/b         %s' % ( a / b ) )
    logger.info ( 'a/(a+b)     %s' % ( a.frac ( b ) ) )
    logger.info ( '(a-b)/(a+b) %s' % ( a.asym ( b ) ) )

    logger.info ( 80*'*')
    
# =============================================================================
##                                                                      The END
# =============================================================================
