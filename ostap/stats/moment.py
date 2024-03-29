#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  moment.py
#  Decorate moment-counters
#  @see Ostap::Math::Moment
#  @see Ostap::Math::Moment_
#  @see Ostap::Math::WMoment
#  @see Ostap::Math::WMoment_
#  @see  Pebay, P., Terriberry, T.B., Kolla, H. et al. 
#        "Numerically stable, scalable formulas for parallel and online 
#        computation of higher-order multivariate central moments with 
#        arbitrary weights". Comput Stat 31, 1305–1325 (2016). 
#  @see https://doi.org/10.1007/s00180-015-0637-z
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-06-08  
# =============================================================================
"""Decorate moment--counters
- see Ostap::Math::Moment
- see Ostap::Math::Moment_
- see Ostap::Math::WMoment
- see Ostap::Math::WMoment_
- see  Pebay, P., Terriberry, T.B., Kolla, H. et al.
   ``Numerically stable, scalable formulas for parallel and online
     computation of higher-order multivariate central moments with 
     arbitrary weights''. Comput Stat 31, 1305–1325 (2016). 
- see https://doi.org/10.1007/s00180-015-0637-z
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2020-06-08"
__all__     = ()
# =============================================================================
from   ostap.core.ostap_types import integer_types, num_types
from   ostap.math.base        import isfinite 
from   ostap.core.core        import Ostap, VE
from   ostap.core.meta_info   import root_version_int, root_info 
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.moment' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
# new stuff: Ostap::Math::Moment_<N> 
# =============================================================================
## get a mean for the moment-counter
#  @code
#  m = ...
#  v = m.mean() 
#  @endcode
#  If order of the moment-counter exceeds 1, the uncertainty is also evaluated 
def _om_mean ( obj ) :
    """Get a mean for the moment-counter
    >>> m = ...
    >>> v = m.mean() 
    - If order of the moment-counter exceeds 1, the uncertainty is also evaluated 
    """
    o = obj.order
    assert 1 <= o , 'mean: the order must be >=1!'
    if root_info < ( 6, 18 ) :
        return obj.mu () 
    return Ostap.Math.Moments.mean ( obj )  

# =============================================================================
## get a variance for the moment-counter
#  @code
#  m = ...
#  v = m.varinace() 
#  @endcode
#  If order of the moment-counter exceeds 3, the uncertainty is also evaluated 
def _om_variance ( obj ) :
    """Get a variance for the moment-counter
    >>> m = ...
    >>> v = m.variance() 
    - If order of the moment-counter exceeds 3, the uncertainty is also evaluated 
    """
    o = obj.order
    assert 2 <= o , 'variance: the order must be >=2!'
    if root_info < ( 6 , 18 ) :
        return obj.moment ( 2 )
    return Ostap.Math.Moments.variance ( obj )

# =============================================================================
## get a skewness for the moment-counter
#  @code
#  m = ...
#  v = m.skewness () 
#  @endcode
def _om_skewness ( obj ) :
    """Get a skewness for the moment-counter
    >>> m = ...
    >>> v = m.skewness() 
    """    
    assert 3 <= obj.order , 'skewness: the order must be >=3!' 
    if root_info < ( 6 , 18 ) :
        m3 = obj.moment ( 3 )
        m2 = obj.moment ( 2 )        
        return m3 / pow ( m2 , 0.5 * 3 )
    return Ostap.Math.Moments.skewness ( obj )  

# =============================================================================
## get an excess kurtosis for the moment-counter
#  @code
#  m = ...
#  v = m.kurtosis () 
#  @endcode
def _om_kurtosis( obj ) :
    """Get an excess  kurtosis  for the moment-counter
    >>> m = ...
    >>> v = m.kurtosis () 
    """    
    assert 4 <= obj.order , 'kurtosis: the order must be >=4!'
    if root_info < ( 6 , 18 ) :
        m4 = obj.moment ( 4 )
        m2 = obj.moment ( 2 )        
        return m4 / ( m2 * m2 ) - 3.0 
    return Ostap.Math.Moments.kurtosis ( obj )  

# =============================================================================
## get unbiased 2nd order moment from the moment-counter 
#  @code
#  m = ...
#  v = m.unbiased_2nd() 
#  @endcode
def _om_u2nd ( obj ) :
    """Get an unbiased 2nd order moment from the moment-counter 
    >>> m = ...
    >>> v = m.unbiased_3nd() 
    """    
    assert 2 <= obj.order , 'unbiased 2nd moment: the order must be >=2!'
    if root_info < ( 6 , 18 ) :
        return obj.moment ( 2 )    
    return Ostap.Math.Moments.unbiased_2nd ( obj )  

# =============================================================================
## get unbiased 3rd order moment from the moment-counter 
#  @code
#  m = ...
#  v = m.unbiased_3rd() 
#  @endcode
def _om_u3rd ( obj ) :
    """Get an unbiased 3rd order moment fro the moment-counter 
    >>> m = ...
    >>> v = m.unbiased_3rd() 
    """    
    assert 3 <= obj.order , 'unbiased 3rd moment: the order must be >=3!' 
    if root_info < ( 6 , 18 ) :
        return obj.moment ( 3 )    
    return Ostap.Math.Moments.unbiased_3rd ( obj )  

# =============================================================================
## get unbiased 4th order moment from the moment-counter 
#  @code
#  m = ...
#  v = m.unbiased_4th() 
#  @endcode
def _om_u4th ( obj ) :
    """Get an unbiased 4th order moment from the moment-counter 
    >>> m = ...
    >>> v = m.unbiased_4th() 
    """    
    assert 4 <= obj.order , 'unbiased 4th moment the order must be >=4!' 
    if root_info < ( 6 , 18 ) :
        return obj.moment ( 4 )    
    return Ostap.Math.Moments.unbiased_4th ( obj )  

# =============================================================================
## get unbiased 5th order moment from the moment-counter 
#  @code
#  m = ...
#  v = m.unbiased_5th() 
#  @endcode
def _om_u5th ( obj ) :
    """Get an unbiased 5th order moment from the moment-counter 
    >>> m = ...
    >>> v = m.unbiased_5th() 
    """    
    assert 5 <= obj.order  , 'unbiased 5th moment: the order must be >=4!' 
    if root_info < ( 6 , 18 ) :
        return obj.moment ( 5 )    
    return Ostap.Math.Moments.unbiased_5th ( obj )  


pos_infinity = float('+Inf')
neg_infinity = float('-Inf')


# =============================================================================
## get central moment 
#  @code
#  m = ...
#  v = m.central_moment() 
#  v = m.moment() 
#  @endcode
if   ( 6 , 22 ) <= root_info :
    ##
    def _om_cm2 ( obj , order  ) :
        assert isinstance  ( order , integer_types ) and order <= obj.order ,\
               'central_moment: invalid order %s/%d' % ( order , obj.order )
        if 2 <= order and obj.empty() : return neg_infinity
        return obj.moment_[order]()
    ##
else :
    ##
    def _om_cm2 ( obj , order  ) :
        assert isinstance  ( order , integer_types ) and order <= obj.order ,\
               'central_moment: invalid order %s/%d' % ( order , obj.order )        
        if 2 <= order and obj.empty() : return neg_infinity
        return obj.moment ( order ) 

_om_cm2.__doc__ = \
   """Get a central moment fro the moment-counter 
   >>> m = ...
   >>> v = m.central_moment ( 3 ) ## ditto 
   >>> v = m.cmoment        ( 3 ) ## ditto 
   """

# =============================================================================
## get central moment 
#  @code
#  m = ...
#  v = m.central_moment() 
#  v = m.moment() 
#  @endcode
if   ( 6 , 22 ) <= root_info :
    ##
    def _om_cm3 ( obj , order  ) :
        assert isinstance  ( order , integer_types ) and order <= obj.order ,\
               'central_moment: invalid order %s/%d' % ( order , obj.order )
        if 2 <= order and not obj.ok () : return neg_infinity
        ##
        return obj.moment_[order]()
    ##
else :
    ##
    def _om_cm3 ( obj , order  ) :    
        assert isinstance  ( order , integer_types ) and order <= obj.order ,\
               'central_moment: invalid order %s/%d' % ( order , obj.order )
        if 2 <= order and not obj.ok () : return neg_infinity
        return obj.moment ( order ) 

_om_cm3.__doc__ = \
    """ Get a central moment fro the moment-counter 
    >>> m = ...
    >>> v = m.central_moment ( 3 ) ## ditto 
    >>> v = m.cmoment        ( 3 ) ## ditto 
    """


# =============================================================================
## get a RMS 
#  @code
#  m = ...
#  v = m.rms  () 
#  @endcode
def _om_rms  ( obj ) :
    """Get a RMS value for the moment-counter
    >>> m = ...
    >>> v = m.rms () 
    """    
    assert 2 <= obj.order  , 'rms: the order must be >=2!'
    ##
    var = obj.variance ()
    if isinstance ( var , VE ) :
        return VE ( max ( var.value() , 0 ) , var.cov2 () ) **  0.5 
    
    return  max ( var , 0 ) **  0.5 

# =============================================================================
## print object as a table
#  @code
#  m = ...
#  t = m.table()
#  @endcode  
def _om_table ( obj , title = '' , prefix = '' , standard = False ) :
    """Print object as a table
    >>> m = ...
    >>> t = m.table()
    """

    from ostap.logger.utils import pretty_float, pretty_ve 

    IM = Ostap.Math.Moments.invalid_moment() 
    
    rows  = []
    
    order  = obj.order
    size   = obj.size   ()
    
    s = obj.size() 
    if 1.e+6 < s :
        field , n = pretty_float ( s * 1.0 ) 
        row = "#entries" , '' if not n else '[10^%+d]' % n , field 
    else :                
        row = "#entries" , '' , '%d' % s     
    rows.append ( row )
                    
    if hasattr  ( obj , 'w2' ) :
            
        w2 = obj.w2 ()
        field , n = pretty_float ( w2 ) 
        row = "sum(w^2)" , '' if not n else '[10^%+d]' % n , field 
        rows.append ( row )
        
    if hasattr  ( obj , 'w' ) :
        
        w = obj.w ()
        field , n = pretty_float ( w ) 
        row = "sum(w)" , '' if not n else '[10^%+d]' % n , field 
        rows.append ( row )
        
    if hasattr  ( obj , 'nEff' ) :
            
        neff = obj.nEff()
        if neff != size and 0 < neff and isfinite ( neff ) :
            field , n = pretty_float ( neff ) 
            row = "nEff" , '' if not n else '[10^%+d]' % n , field 
            rows.append ( row )
            
    if 1 <= order and 1 <= size and obj.ok () and hasattr ( obj , 'mean' )  : 
        
        v  = obj.mean ()
        vv = float   ( v ) 
        if IM != float ( v ) and isfinite ( vv ) : 
            if isinstance ( v , VE ) : field , n = pretty_ve    ( v )
            else                     : field , n = pretty_float ( v )
            row = "mean" , '' if not n else '[10^%+d]' % n , field 
            rows.append ( row )
            
    if 2 <= order and 2 <= size and obj.ok() :

        if hasattr ( obj , 'rms' ) : 
            v  = obj.rms (    )
            vv = float   ( v  )                             
            if isfinite  ( vv ) and IM != vv and 0 <= vv :                    
                if isinstance ( v , VE ) : field , n = pretty_ve    ( v )
                else                     : field , n = pretty_float ( v )
                row = "rms"  , '' if not n else '[10^%+d]' % n , field 
                rows.append ( row )
                
        if hasattr ( obj , 'variance' ) :
            v  = obj.variance (    )
            vv = float   ( v  )                             
            if isfinite  ( vv ) and IM != vv and 0 <= vv :                    
                if isinstance ( v , VE ) : field , n = pretty_ve    ( v )
                else                     : field , n = pretty_float ( v )
                row = "variance"  , '' if not n else '[10^%+d]' % n , field 
                rows.append ( row )
                
    if 3 <= order and 3 <= size and obj.ok () and hasattr  ( obj , 'skewness' ) :
        
        v  = obj.skewness ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != vv :
            if isinstance ( v , VE ) : field , n = pretty_ve    ( v )
            else                     : field , n = pretty_float ( v )
            row = "skewness"  , '' if not n else '[10^%+d]' % n , field 
            rows.append ( row )
            
    if 4 <= order and 4 <= size and  obj.ok () and hasattr ( obj , 'kurtosis' ) : 
        
        v  = obj.kurtosis ()
        vv = float   ( v )                                     
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = pretty_ve    ( v )
            else                     : field , n = pretty_float ( v )
            row = "kurtosis"  , '' if not n else '[10^%+d]' % n , field 
            rows.append ( row )

    if 2 <= order and 2 <= size and obj.ok () and hasattr ( obj , 'unbiased_2nd' ) :
        
        v  = obj.unbiased_2nd ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = pretty_ve    ( v )
            else                     : field , n = pretty_float ( v )
            row = "M[2]/unb" , '' if not n else '[10^%+d]' % n , field 
            rows.append ( row )

    if 3 <= order and 3 <= size and obj.ok () and hasattr ( obj , 'unbiased_3rd' ) :
        
        v  = obj.unbiased_3rd ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = pretty_ve    ( v )
            else                     : field , n = pretty_float ( v )
            row = "M[3]/unb" , '' if not n else '[10^%+d]' % n , field 
            rows.append ( row )

    if 4 <= order and 4 <= size and obj.ok () and hasattr ( obj , 'unbiased_4th' ) :
        
        v  = obj.unbiased_4th ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = pretty_ve    ( v )
            else                     : field , n = pretty_float ( v )
            row = "M[4]/unb" , '' if not n else '[10^%+d]' % n , field 
            rows.append ( row )
        
    if 5 <= order and 5 <= size and obj.ok () and hasattr ( obj , 'unbiased_5th' ) : 

        v  = obj.unbiased_5th ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = pretty_ve    ( v )
            else                     : field , n = pretty_float ( v )
            row = "M[5](unb)" , '' if not n else '[10^%+d]' % n , field 
            rows.append ( row )

    fmt = '[%d]'
    if   10  <= order < 100   : fmt = '[%02d]'
    elif 100 <= order < 1000  : fmt = '[%03d]'
  
    for i in range ( 2 , order + 1 ) :
        
        if standard  :
            
            v  = obj.std_moment ( i )
            vv = float   ( v )                         
            if isfinite ( vv ) and IM != float ( v ) :                
                if isinstance ( v , VE ) : field , n = pretty_ve    ( v )
                else                     : field , n = pretty_float ( v )
                row = ( "M" + fmt + "/std" )  % i , '' if not n else '[10^%+d]' % n , field 
                rows.append ( row )

        else :

            v  = obj.central_moment ( i )
            vv = float   ( v )                         
            if isfinite ( vv ) and IM != float ( v ) :                
                if isinstance ( v , VE ) : field , n = pretty_ve    ( v )
                else                     : field , n = pretty_float ( v )
                row = ( "M" + fmt + "/raw" )  % i  , '' if not n else '[10^%+d]' % n , field 
                rows.append ( row )

        
    rows = tuple ( [ ( '' , '' , 'value') ] +  rows  )
        
    import ostap.logger.table as T

    tit = title 
    if not title :
        if   isinstance ( obj , Ostap.Math. Moment ) : tit =  'Moment_[%d]' % obj.order 
        elif isinstance ( obj , Ostap.Math.WMoment ) : tit = 'WMoment_[%d]' % obj.order 
    
    return T.table ( rows  , title = tit , prefix = prefix )

Ostap.Math.Moment.unbiased_2nd   = _om_u2nd 
Ostap.Math.Moment.unbiased_3rd   = _om_u3rd
Ostap.Math.Moment.unbiased_4th   = _om_u4th 
Ostap.Math.Moment.unbiased_5th   = _om_u5th 

Ostap.Math.Moment.mean           = _om_mean    
Ostap.Math.Moment.rms            = _om_rms 
Ostap.Math.Moment.variance       = _om_variance
Ostap.Math.Moment.skewness       = _om_skewness
Ostap.Math.Moment.kurtosis       = _om_kurtosis
Ostap.Math.Moment.cmoment        = _om_cm2
Ostap.Math.Moment.central_moment = _om_cm2
Ostap.Math.Moment.table          = _om_table

Ostap.Math.WMoment.mean           = _om_mean    
Ostap.Math.WMoment.rms            = _om_rms 
Ostap.Math.WMoment.variance       = _om_variance
Ostap.Math.WMoment.skewness       = _om_skewness
Ostap.Math.WMoment.kurtosis       = _om_kurtosis
Ostap.Math.WMoment.cmoment        = _om_cm3
Ostap.Math.WMoment.central_moment = _om_cm3
Ostap.Math.WMoment.table          = _om_table


for t in ( Ostap.Math.WMoment ,
           Ostap.Math. Moment ) :
    t.__str__  = _om_table
    t.__repr__ = _om_table

M0  = Ostap.Math.Moment_(0)
M1  = Ostap.Math.Moment_(1)
WM0 = Ostap.Math.Moment_(0)
WM1 = Ostap.Math.Moment_(1)
for m in ( M0 , M1 , WM0 , WM1 ) :
    m.mean  = _om_mean
    m.table = _om_table

if not hasattr (  M0 , 'order' ) :  M0.order = 0
if not hasattr (  M1 , 'order' ) :  M1.order = 1
if not hasattr ( WM0 , 'order' ) : WM0.order = 0
if not hasattr ( WM1 , 'order' ) : WM1.order = 1

_decorated_classes = (
    Ostap.Math.Moment     ,
    Ostap.Math.WMoment    ,
    Ostap.Math.Moment_(0) , 
    Ostap.Math.Moment_(1) ,
    Ostap.Math.Moment_(0) , 
    Ostap.Math.Moment_(1) , 
    )

_new_methods_ = (
    ##
    Ostap.Math.Moment.unbiased_2nd   , 
    Ostap.Math.Moment.unbiased_3rd   ,
    Ostap.Math.Moment.unbiased_4th   ,
    Ostap.Math.Moment.unbiased_5th   ,
    ##
    Ostap.Math.Moment.mean           ,
    Ostap.Math.Moment.rms            ,
    Ostap.Math.Moment.variance       ,  
    Ostap.Math.Moment.skewness       ,
    Ostap.Math.Moment.kurtosis       ,
    Ostap.Math.Moment.cmoment        ,
    Ostap.Math.Moment.central_moment , 
    Ostap.Math.Moment.table          ,
    ##
    Ostap.Math.WMoment.mean           ,
    Ostap.Math.WMoment.rms            ,
    Ostap.Math.WMoment.variance       ,  
    Ostap.Math.WMoment.skewness       ,
    Ostap.Math.WMoment.kurtosis       ,
    Ostap.Math.WMoment.cmoment        ,
    Ostap.Math.WMoment.central_moment , 
    Ostap.Math.WMoment.table          ,
    ##
    )
                          
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
