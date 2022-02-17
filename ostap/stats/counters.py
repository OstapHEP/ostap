#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  @see Ostap::StatEntity
#  @see Ostap::WStatEntity
#  @see Ostap::NStatEntity
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
"""Simple counters 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = ""
# =============================================================================
__all__     = (
    'SE'  , ## simple smart counter (1-bi histo)  C++ Ostap::StatEntity 
    'WSE' , ## simple smart counter with weight :     Ostap::WStatEntity 
    'NSE' , ## simple smart running counter     :     Ostap::NStatEntity 
    ) 
# =============================================================================
import ROOT, cppyy
from   ostap.math.ve   import Ostap, VE
from   ostap.math.base import isequal, isequalf  
# =============================================================================
_new_methods_ = []
# =============================================================================
SE    = Ostap.StatEntity 
WSE   = Ostap.WStatEntity 
NSE   = Ostap.NStatEntity 

# =============================================================================
# temporary trick, to be removed 
# =============================================================================

SE .__repr__ = lambda s : 'Stat: '+ s.toString()
SE .__str__  = lambda s : 'Stat: '+ s.toString()

# =============================================================================
# minor decoration for StatEntity 
# ============================================================================= 
if not hasattr ( SE , '_orig_sum'  ) : 
    _orig_sum    = SE.sum
    SE._orig_sum = _orig_sum
    
if not hasattr ( SE , '_orig_mean' ) : 
    _orig_mean    = SE.mean
    SE._orig_mean = _orig_mean


SE. sum     = lambda s : VE ( s._orig_sum  () , s.sum2 ()      )
SE. minmax  = lambda s :    ( s.min        () , s.max  ()      ) 
SE. mean    = lambda s : VE ( s._orig_mean () , s.meanErr()**2 )

# ==============================================================================
## Update the counter
#  @code
#  >>> cnt  = ...
#  >>> cnt += value
#  @endcode
def _se_iadd_ ( self , other ) :
    """Update the counter
    >>> cnt  = ...
    >>> cnt += value 
    """
    return self.add ( other )

# ==============================================================================
## Update the counter
#  @code
#  >>> cnt  = ...
#  >>> cnt -= value
#  @endcode
def _se_isub_ ( self , other ) :
    """Update the counter
    >>> cnt  = ...
    >>> cnt -= value 
    """
    return self.add ( -other )

# ==============================================================================
## add two counters
#  @code
#  >>> cnt1 = ...
#  >>> cnt2 = ...
#  >>> cnt  = cnt1 +  cnt2 
#  @endcode
def _se_add_ ( self , other ) :
    """Update the counter
    >>> cnt  = ...
    >>> cnt += value 
    """
    r  = SE ( self )
    r.add ( other )
    return r

# ==============================================================================
## add two counters
#  @code
#  >>> cnt1 = ...
#  >>> cnt2 = ...
#  >>> cnt  = cnt1 +  cnt2 
#  @endcode
def _wse_add_ ( self , other ) :
    """Update the counter
    >>> cnt  = ...
    >>> cnt += value 
    """
    r = WSE ( self )
    r.add ( other )
    return r
# =============================================================================

SE  .__iadd__ =  _se_iadd_
SE  .__add__  =  _se_add_
WSE .__add__  = _wse_add_


# ==============================================================================
## numbers are close enough?
def closef ( a , b , N = 100 ) :
    return ( a == b ) or isequalf ( a , b ) or abs ( a - b ) * N <= abs ( a ) + abs ( b )
    
# =============================================================================
## equal counters?
def _se_eq_ ( s1 , s2 ) :
    """Numerically equal counters?"""
    ##
    if not isinstance ( s2 , SE ) : return NotImplemented
    ## 
    ## closef   ( s1.mu2 () , s2.mu2 () ) and \
    return s1.n     () == s2.n   ()           and \
           isequalf ( s1.mu  () , s2.mu  () ) and \
           isequalf ( s1.mu2 () , s2.mu2 () ) and \
           s1.min () == s2.min ()             and \
           s1.max () == s2.max ()
# =============================================================================
## non-equal counters?
def _se_ne_ ( s1 , s2 ) :
    """Numerically non-equal counters?
    """
    if not isinstance ( s2 , SE ) : return NotImplemented 
    return not ( s1 == s2 )

SE.__eq__ = _se_eq_
SE.__ne__ = _se_ne_

_new_methods_ += [
    SE.sum       ,
    SE.mean      ,
    SE.minmax    ,
    SE.__add__   ,
    SE.__iadd__  ,    
    SE.__isub__  ,    
    SE.__repr__  ,
    SE.__str__   ,
    SE.__str__   ,
    SE.__eq__    , 
    SE.__ne__    ,
    ]

# =============================================================================
## equal counters?
def _wse_eq_ ( s1 , s2 ) :
    """Numerically equal counters?"""
    ##
    if not isinstance ( s2 , WSE ) : return NotImplemented
    ##
    ## closef     ( s1.mu2 () , s2.mu2 () ) and \
    return s1.n       () == s2.n       ()       and \
           isequalf   ( s1.mu  () , s2.mu  () ) and \
           isequalf   ( s1.mu2 () , s2.mu2 () ) and \
           s1.weights () == s2.weights ()       and \
           s1.values  () == s2.values  ()

# =============================================================================
## non-equal counters?
def _wse_ne_ ( s1 , s2 ) :
    """Numerically non-equal counters?"""
    ##
    if not isinstance ( s2 , WSE ) : return NotImplemented
    ##    
    return not ( s1 == s2 )

WSE.__eq__ = _wse_eq_
WSE.__ne__ = _wse_ne_

# =============================================================================
# minor decoration for WStatEntity 
# ============================================================================= 
if not hasattr ( WSE , '_orig_sum'  ) : 
    _orig_sum     = WSE.sum
    WSE._orig_sum = _orig_sum

if not hasattr ( WSE , '_orig_mean' ) : 
    _orig_mean_wse = WSE.mean
    WSE._orig_mean = _orig_mean_wse
    
WSE. sum     = lambda s : VE ( s._orig_sum  () , s.sum2()       )
WSE. mean    = lambda s : VE ( s._orig_mean () , s.meanErr()**2 )
WSE. minmax  = lambda s :            s.values  ().minmax() 
WSE.__repr__ = lambda s : 'WStat: '+ s.toString()
WSE.__str__  = lambda s : 'WStat: '+ s.toString()

_new_methods_ += [
    WSE.sum      ,
    WSE.mean     ,
    WSE.minmax   ,
    WSE.__add__  ,
    WSE.__repr__ ,
    WSE.__str__  ,
    WSE.__str__  ,
    WSE.__eq__   ,
    WSE.__ne__   ,
    ]

# =============================================================================
## Count iterable
#  @code
#  c = SE.count ( [1,2,3,4] )  
#  @endcode
def se_count ( values ) :
    """Count iterable
    >>> c = SE.count ( [1,2,3,4] )  
    """
    cnt = SE ()
    for v in values : cnt.add ( v ) 
    return cnt
SE.count = staticmethod ( se_count )

_new_methods_ += [
    SE.count ,
    ]

_new_methods_ = tuple ( _new_methods_ )

# =============================================================================
_decorated_classes_ = (
    SE , WSE , NSE 
    )
# =============================================================================
if '__main__' == __name__  :

    from ostap.logger.logger import getLogger
    logger = getLogger( 'ostap.stats.counters' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    import random
    
    cnt = SE() 
    from builtins import range 
    for i in range(10000) :
        cnt += random.gauss(1,1)
        
    logger.info ( 'Counter: %s' % cnt        ) 
    logger.info ( 'Mean   : %s' % cnt.mean() ) 
    logger.info ( 'RMS    : %s' % cnt.rms () ) 
    
    logger.info (80*'*')
    
# =============================================================================
##                                                                      The END 
# =============================================================================
