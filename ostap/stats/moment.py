#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  moment.py
#  Decorate moment0-counters
#  @see Ostap::Math::Moment
#  @see Ostap::Math::Moment_
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-06-08  
# =============================================================================
"""Decorate moment0-counters
- see Ostap::Math::Moment
- see Ostap::Math::Moment_
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2020-06-08"
__all__     = ()
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.moment' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
from   ostap.core.ostap_types import integer_types, num_types 
from   ostap.core.core        import Ostap, VE 
# =============================================================================
# new stuff: Ostap::Math::Moment_<N> 
# =============================================================================
## get a variance for the moment-counter
#  @code
#  m = ...
#  v = m.varinace() 
#  @encode
#  If order of the moment-counter exceeds 4, the uncertainty is also evaluated 
def _om_variance ( obj ) :
    """Get a variance for the moment-counter
    >>> m = ...
    >>> v = m.varinace() 
    - If order of the moment-counter exceeds 4, the uncertainty is also evaluated 
    """
    o = obj.order
    assert 2 <= o , 'variance: the order must be >=2!'
    return Ostap.Math.Moments.variance ( obj )  

# =============================================================================
## get a skewness for the moment-counter
#  @code
#  m = ...
#  v = m.skewness () 
#  @encode
def _om_skewness ( obj ) :
    """Get a skewness for the moment-counter
    >>> m = ...
    >>> v = m.skewness() 
    """    
    assert 3 <= obj.order , 'skewness: the order must be >=3!' 
    return Ostap.Math.Moments.skewness ( obj )  

# =============================================================================
## get an excess  kurtosis for the moment-counter
#  @code
#  m = ...
#  v = m.kurtosis () 
#  @encode
def _om_kurtosis( obj ) :
    """Get an excess  kurtosis  for the moment-counter
    >>> m = ...
    >>> v = m.kurtosis () 
    """    
    assert 4 <= obj.order , 'kurtosis: the order must be >=4!' 
    return Ostap.Math.Moments.kurtosis ( obj )  

# =============================================================================
## get unbiased 2nd order moment from the moment-counter 
#  @code
#  m = ...
#  v = m.unbiased_2nd() 
#  @encode
def _om_u2nd ( obj ) :
    """Get an unbiased 2nd order moment from the moment-counter 
    >>> m = ...
    >>> v = m.unbiased_3nd() 
    """    
    assert 2 <= obj.order , 'unbiased 2nd moment: the order must be >=2!' 
    return Ostap.Math.Moments.unbiased_2nd ( obj )  

# =============================================================================
## get unbiased 3rd order moment from the moment-counter 
#  @code
#  m = ...
#  v = m.unbiased_3rd() 
#  @encode
def _om_u3rd ( obj ) :
    """Get an unbiased 3rd order moment fro the moment-counter 
    >>> m = ...
    >>> v = m.unbiased_3rd() 
    """    
    assert 3 <= obj.order , 'unbiased 3rd moment: the order must be >=3!' 
    return Ostap.Math.Moments.unbiased_3rd ( obj )  

# =============================================================================
## get unbiased 4th order moment from the moment-counter 
#  @code
#  m = ...
#  v = m.unbiased_4th() 
#  @encode
def _om_u4th ( obj ) :
    """Get an unbiased 4th order moment from the moment-counter 
    >>> m = ...
    >>> v = m.unbiased_4th() 
    """    
    assert 4 <= obj.order , 'unbiased 4th moment the order must be >=4!' 
    return Ostap.Math.Moments.unbiased_4th ( obj )  

# =============================================================================
## get unbiased 5th order moment from the moment-counter 
#  @code
#  m = ...
#  v = m.unbiased_5th() 
#  @encode
def _om_u5th ( obj ) :
    """Get an unbiased 5th order moment from the moment-counter 
    >>> m = ...
    >>> v = m.unbiased_5th() 
    """    
    assert 5 <= obj.order  , 'unbiased 5th moment: the order must be >=4!' 
    return Ostap.Math.Moments.unbiased_5th ( obj )  

# =============================================================================
## get central moment 
#  @code
#  m = ...
#  v = m.unbiased_5th() 
#  @encode
def _om_cm ( obj , order  ) :
    """Get an unbiased 5th order moment fro the moment-counter 
    >>> m = ...
    >>> v = m.central_moment ( 3 ) ## ditto 
    >>> v = m.cmoment        ( 3 ) ## ditto 
    """
    assert isinstance  ( order , integer_types ) and 2<= order , 'Invalid order %s'% order
    assert order <= obj.order , 'central_moment: invalid order cmbiarions %s/%s' % ( order , obj.order )

    if order * 2  <= obj.order :
        T = Ostap.Math.Moments._central_moment_2 [ order , obj.order ]
        M = Ostap.Math.Moments()
        return T ( M , obj )

    return obj.moment ( order ) 

# =============================================================================
## get a mean 
#  @code
#  m = ...
#  v = m.mean () 
#  @encode
def _om_mean ( obj ) :
    """Get a mean value for the moment-counter
    >>> m = ...
    >>> v = m.mean () 
    """    
    assert 1 <= obj.order  , 'mean: the order must be >=1!'
    ##
    if  1 == obj.order or obj.size() < 2 : return obj.mu()
    ##   
    cov2 = float ( obj.variance() ) / obj.size() 
    ##
    return VE ( obj.mu () , cov2 ) 

# =============================================================================
## get a RMS 
#  @code
#  m = ...
#  v = m.rms  () 
#  @encode
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
def _om_table ( obj , title = '' , prefix = '' ) :
    """Print object as a table
    >>> m = ...
    >>> t = m.table()
    """
    rows = []
    item  = obj
    while 1 < 2 :
        
        order  = item.order
        moment = item.moment()  
        
        if   0 == order :
            row = "#"        , "%d" % item.size() 
            rows.append ( row )
        elif 1 == order :
            v = obj.mean()
            if  isinstance ( v , VE ) : field = v.toString( "%.6g +/- %.6g" )
            else                      : field = "%.6g" % v 
            row = "Mean"     , field 
            rows.append ( row )
        elif 2 == order :
            v = obj.variance ()
            if  isinstance ( v , VE ) : field = v.toString( "%.6g +/- %.6g" )
            else                      : field = "%.6g" % v            
            row = "Variance" , field 
            rows.append ( row )
            ##
            v = obj.rms  ()
            if  isinstance ( v , VE ) : field = v.toString( "%.6g +/- %.6g" )
            else                      : field = "%.6g" % v            
            row = "RMS" , field 
            rows.append ( row )
        elif 3 == order :
            v = obj.skewness ()
            if  isinstance ( v , VE ) : field = v.toString( "%.6g +/- %.6g" )
            else                      : field = "%.6g" % v            
            row = "Skewness" , field 
            rows.append ( row )
        elif 4 == order :
            v = obj.kurtosis ()
            if  isinstance ( v , VE ) : field = v.toString( "%.6g +/- %.6g" )
            else                      : field = "%.6g" % v            
            row = "Kurtosis" , field 
            rows.append ( row )
        else :
            v = obj.cmoment ( order )
            if  isinstance ( v , VE ) : field = v.toString( "%.6g +/- %.6g" )
            else                      : field = "%.6g" % v            
            row = "M<%d>" % order , field 
            rows.append ( row )

        if hasattr ( item , 'previous' ) : item = item.previous()
        else                             : break
        
        
    rows = tuple ( [ ('', 'Value') ] +  [ r for r in reversed ( rows ) ] )
        
    import ostap.logger.table as T
        
    return T.table ( rows  , title = title , prefix = prefix )

  
Ostap.Math.Moment.unbiased_2nd   = _om_u2nd 
Ostap.Math.Moment.unbiased_3rd   = _om_u3rd
Ostap.Math.Moment.unbiased_4th   = _om_u4th 
Ostap.Math.Moment.unbiased_5th   = _om_u5th 

Ostap.Math.Moment.mean           = _om_mean    
Ostap.Math.Moment.rms            = _om_rms 
Ostap.Math.Moment.variance       = _om_variance
Ostap.Math.Moment.skewness       = _om_skewness
Ostap.Math.Moment.kurtosis       = _om_kurtosis
Ostap.Math.Moment.cmoment        = _om_cm
Ostap.Math.Moment.central_moment = _om_cm
Ostap.Math.Moment.table          = _om_table

M0 = Ostap.Math.Moment_(0)
M1 = Ostap.Math.Moment_(1)
for m in ( M0 , M1 ) :
    m.mean = _om_mean
M0.order = 0
M1.order = 1

_decorated_classes = (
    Ostap.Math.Moment ,
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
    )
                          
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
