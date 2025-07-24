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
""" Decorate moment--counters
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
from   ostap.math.base        import isfinite, isequal
from   ostap.core.core        import Ostap, VE
from   ostap.logger.pretty    import pretty_float
from   ostap.logger.symbols   import times, sum_symbol 
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.moment' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
invalid_moment = Ostap.Math.Moments.invalid_moment () 
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
    """ Get a mean for the moment-counter
    >>> m = ...
    >>> v = m.mean() 
    - If order of the moment-counter exceeds 1, the uncertainty is also evaluated 
    """
    o = obj.order
    assert 1 <= o , 'mean: the order must be >=1!'
    return Ostap.Math.Moments.mean ( obj )

# ===========================================================================
## get min/max  for the moment-counter
#  @code
#  m = ...
#  mn, mx  = m.minmax() 
#  @endcode
def _om_minmax ( obj ) :
    """ Get min/.max for the moment-counter
    >>> m = ...
    >>> mn , mx  = m.menamax () 
    """
    o = obj.order
    assert 1 <= o , 'minmax: the order must be >=1!'
    return m.min () , m.max() 

# =============================================================================
## get a variance for the moment-counter
#  @code
#  m = ...
#  v = m.varinace() 
#  @endcode
#  If order of the moment-counter exceeds 3, the uncertainty is also evaluated 
def _om_variance ( obj ) :
    """ Get a variance for the moment-counter
    >>> m = ...
    >>> v = m.variance() 
    - If order of the moment-counter exceeds 3, the uncertainty is also evaluated 
    """
    o = obj.order
    assert 2 <= o , 'variance: the order must be >=2!'
    return Ostap.Math.Moments.variance ( obj )

# =============================================================================
## get a RMS 
#  @code
#  m = ...
#  v = m.rms  () 
#  @endcode
def _om_rms  ( obj ) :
    """ Get a RMS value for the moment-counter
    >>> m = ...
    >>> v = m.rms () 
    """    
    assert 2 <= obj.order  , 'rms: the order must be >=2!'
    ##
    var = Ostap.Math.Moments.variance ( obj )    
    if   isinstance ( var , VE ) :
        if 0 > var.value()  : var = VE ( 0 , var.cov2() )
    elif 0 > var : var = 0.0
    
    return  var ** 0.5 

# =============================================================================
## get a skewness for the moment-counter
#  @code
#  m = ...
#  v = m.skewness () 
#  @endcode
def _om_skewness ( obj ) :
    """ Get a skewness for the moment-counter
    >>> m = ...
    >>> v = m.skewness() 
    """    
    assert 3 <= obj.order , 'skewness: the order must be >=3!' 
    return Ostap.Math.Moments.skewness ( obj )  

# =============================================================================
## get an excess kurtosis for the moment-counter
#  @code
#  m = ...
#  v = m.kurtosis () 
#  @endcode
def _om_kurtosis ( obj ) :
    """ Get an excess  kurtosis  for the moment-counter
    >>> m = ...
    >>> v = m.kurtosis () 
    """    
    assert 4 <= obj.order , 'kurtosis: the order must be >=4!'
    return Ostap.Math.Moments.kurtosis ( obj )  

# =============================================================================
## Get 1st cumulant, well it is equal to the mean 
#  @code
#  m = ...
#  v = m.cumulant_1st() 
#  @endcode
def _om_cumulant_1st( obj ) :
    """ Get 1st unbiased cumulant, well it is equal to the mean
    
    >>> m = ...
    >>> v = m.cumulant_1st() 
    """
    return obj.mean() 

# =============================================================================
## Get 2nd cumulant, well it is equal to the variance 
#  @code
#  m = ...
#  v = m.cumulant_2nd() 
#  @endcode
def _om_cumulant_2nd ( obj ) :
    """ Get 2nd unbiased cumulant, well it is equal to the variance 
    
    >>> m = ...
    >>> v = m.cumulant_2nd() 
    """
    assert 2 <= obj.order , 'cumulant_2nd: the order must be >=2!'
    return Ostap.Math.Moments.cumulant_2nd  ( obj )  

# =============================================================================
## Get 3nd unbiased cumulant, well it is equal to the 3rd central moment  
#  @code
#  m = ...
#  v = m.cumulant_3rd() 
#  @endcode
def _om_cumulant_3rd ( obj ) :
    """ Get 3rd unbiased cumulant, well it is equal to the 3rd central moment 
    
    >>> m = ...
    >>> v = m.cumulant_3rd() 
    """
    assert 3 <= obj.order , 'cumulant_3rd: the order must be >=3!'
    return Ostap.Math.Moments.cumulant_3rd ( obj )  

# =============================================================================
## Get 4th nb unbiased cumulant,
#  @code
#  m = ...
#  v = m.cumulant_4th() 
#  @endcode
def _om_cumulant_4th ( obj ) :
    """ Get 4th unbiased cumulant, well it is equal to the 3rd central moment 
    
    >>> m = ...
    >>> v = m.cumulant_4th() 
    """
    assert 4 <= obj.order , 'cumulant_4th: the order must be >=3!'
    return Ostap.Math.Moments.cumulant_4th ( obj )  

# =============================================================================
## Get Nth cumulant,
#  @code
#  m = ...
#  v = m.cumulant(6) 
#  @endcode
def _om_cumulant_( obj , order ) :
    """ Get Nth cumulant 
    
    >>> m = ...
    >>> v = m.cumulant_4th() 
    """
    O = min ( obj.order , 10 ) 
    assert isinstance ( order , integer_types ) and 1 <= order <= O, 'cumulant: invalid order!'
    
    if   1 == order : return obj.cumulant_1st ()
    elif 2 == order : return obj.cumulant_2nd ()
    elif 3 == order : return obj.cumulant_3rd ()
    elif 4 == order : return obj.cumulant_4th ()
    ## 
    return obj.cumulant_[order]() 

# =============================================================================
## get unbiased 2nd order moment from the moment-counter 
#  @code
#  m = ...
#  v = m.unbiased_2nd() 
#  @endcode
def _om_u2nd ( obj ) :
    """ Get an unbiased 2nd order moment from the moment-counter 
    >>> m = ...
    >>> v = m.unbiased_3nd() 
    """    
    assert 2 <= obj.order , 'unbiased 2nd moment: the order must be >=2!'
    ## 
    r2 = Ostap.Math.Moments.unbiased_2nd ( obj )
    if isinstance ( r2 , VE ) : return r2
    ## 
    m2 = Ostap.Math.Moments.moment[2]( obj )
    if isinstance ( m2 , VE ) and 0 < m2.cov2() : return VE ( r2 , m2.cov2() )
    return r2

# ===========================================================================
## get unbiased 3rd order moment from the moment-counter 
#  @code
#  m = ...
#  v = m.unbiased_3rd() 
#  @endcode
def _om_u3rd ( obj ) :
    """ Get an unbiased 3rd order moment fro the moment-counter 
    >>> m = ...
    >>> v = m.unbiased_3rd() 
    """    
    assert 3 <= obj.order , 'unbiased 3rd moment: the order must be >=3!'
    ## 
    r3 = Ostap.Math.Moments.unbiased_3rd ( obj )  
    if isinstance ( r3 , VE ) : return r3
    ## 
    m3 = Ostap.Math.Moments.moment[3]( obj )
    if isinstance ( m3 , VE ) and 0 < m3.cov2() : return VE ( r3 , m3.cov2() )
    return r3

# =============================================================================
## get unbiased 4th order moment from the moment-counter 
#  @code
#  m = ...
#  v = m.unbiased_4th() 
#  @endcode
def _om_u4th ( obj ) :
    """ Get an unbiased 4th order moment from the moment-counter 
    >>> m = ...
    >>> v = m.unbiased_4th() 
    """    
    assert 4 <= obj.order , 'unbiased 4th moment the order must be >=4!'
    ## 
    r4 = Ostap.Math.Moments.unbiased_4th ( obj )
    if isinstance ( r4 , VE ) : return r4
    ## 
    m4 = Ostap.Math.Moments.moment[4]( obj )
    if isinstance ( m4 , VE ) and 0 < m4.cov2() : return VE ( r4 , m4.cov2() )
    return r4
    
# =============================================================================
## get unbiased 5th order moment from the moment-counter 
#  @code
#  m = ...
#  v = m.unbiased_5th() 
#  @endcode
def _om_u5th ( obj ) :
    """ Get an unbiased 5th order moment from the moment-counter 
    >>> m = ...
    >>> v = m.unbiased_5th() 
    """    
    assert 5 <= obj.order  , 'unbiased 5th moment: the order must be >=4!'
    ## 
    r5 = Ostap.Math.Moments.unbiased_5th ( obj )
    if isinstance ( r5 , VE ) : return r5
    ## 
    m5 = Ostap.Math.Moments.moment[5]( obj )
    if isinstance ( m5 , VE ) and 0 < m5.cov2() : return VE ( r5 , m5.cov2() )
    return r5
    
# =============================================================================
## get central moment 
#  @code
#  m = ...
#  v = m.central_moment() 
#  v = m.moment() 
#  @endcode
def _om_cm2 ( obj , order  ) :
   """ Get a central moment fro the moment-counter 
   >>> m = ...
   >>> v = m.central_moment ( 3 ) ## ditto 
   >>> v = m.cmoment        ( 3 ) ## ditto 
   """
   assert isinstance  ( order , integer_types ) and order <= obj.order ,\
    'central_moment: invalid order %s/%d' % ( order , obj.order )
   if 2 <= order and obj.empty() : return invalid_moment
   ## 
   return Ostap.Math.Moments.moment [ order ]( obj )

# =============================================================================
## get central moment 
#  @code
#  m = ...
#  v = m.central_moment() 
#  v = m.moment() 
#  @endcode
def _om_cm3 ( obj , order  ) :
    """ Get a central moment for the moment-counter 
    >>> m = ...
    >>> v = m.central_moment ( 3 ) ## ditto 
    >>> v = m.cmoment        ( 3 ) ## ditto 
    """
    assert isinstance  ( order , integer_types ) and order <= obj.order ,\
        'central_moment: invalid order %s/%d' % ( order , obj.order )
    if 2 <= order and not obj.ok () : return invalid_moment
    ##
    return Ostap.Math.Moments.moment [ order ]( obj )


# =============================================================================
## get centralized moment 
#  @code
#  m = ...
#  v = m.centralized_moment ( center = 15. ) 
#  @endcode
def _om_cm4 ( obj , order , center ) :
   """ Get the centralized  moment f
   >>> m = ...
   >>> v = m.centralized_moment ( 3 , center = 15.0 ) ## ditto 
   """
   assert isinstance  ( order , integer_types ) and order <= obj.order ,\
       'centralized_moment: invalid order %s/%d' % ( order , obj.order )
   assert isinstance  ( cener , num_types ) ,\
       'centralized_moment: invalid center %s' % center
   ## 
   return Ostap.Math.centralized_moment [ order ]( center )

# =============================================================================
## Get the standartized moment of order `order` 
def _om_std ( obj , order ) :
    """ Get the standartized moment of order `order` 
    """
    assert isinstance ( order , integer_types ) and 0 <= order <= obj.order , \
        "Invalid `order`!"

    if   0 == order : return 1            ## RETURN 
    elif 1 == order : return 0            ## RETURN 
    elif 2 == order : return 1            ## RETURN

    return Ostap.Math.Moments.std_moment [ order ] ( obj  )
    
# =============================================================================
## print object as a table
#  @code
#  m = ...
#  t = m.table()
#  @endcode  
def _om_table ( obj , * , 
                title     = ''    ,
                prefix    = ''    ,
                standard  = True  ,
                cumulants = False ,
                style     = None  ) :
    """ Print object as a table
    >>> m = ...
    >>> t = m.table()
    """

    IM = Ostap.Math.Moments.invalid_moment() 
    
    rows   = []
    
    order  = obj.order
    size   = obj.size   ()

    fmt_factor = '[%s10^%%+d]' % times 
    s = obj.size() 
    if 1.e+6 < s :
        field , n = pretty_float ( s * 1.0 ) 
        row = "#entries" , field    , '' if not n else fmt_factor % n 
    else :                
        row = "#entries" , '%d' % s , ''      
    rows.append ( row )
                    
    if hasattr  ( obj , 'w' ) :
        
        w = obj.w ()
        field , n = pretty_float ( w ) 
        row = "%sw" % sum_symbol  , field , '' if not n else fmt_factor % n
        rows.append ( row )
        
    if hasattr  ( obj , 'w2' ) :
            
        w2 = obj.w2 ()
        field , n = pretty_float ( w2 ) 
        row = "%sw^2" % sum_symbol , field , '' if not n else fmt_factor % n 
        rows.append ( row )
        
    if 1 <= size and obj.ok () and hasattr ( obj , 'wmin' ) :
        
        v  = obj.wmin ()
        vv = float   ( v ) 
        if IM != float ( v ) and isfinite ( vv ) : 
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = "weight/min" , field , '' if not n else fmt_factor % n 
            rows.append ( row )
        
    if 1 <= size and obj.ok () and hasattr ( obj , 'wmax' ) :
        
        v  = obj.wmax ()
        vv = float   ( v ) 
        if IM != float ( v ) and isfinite ( vv ) : 
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = "weight/max" , field , '' if not n else fmt_factor % n 
            rows.append ( row )

    if hasattr  ( obj , 'nEff' ) :
            
        neff = obj.nEff()
        if neff != size and 0 < neff and isfinite ( neff ) :
            field , n = pretty_float ( neff ) 
            row = "nEff" , field , '' if not n else fmt_factor % n 
            rows.append ( row )

    if 1 <= order and 1 <= size and obj.ok () and hasattr ( obj , 'min' ) :
        
        v  = obj.min ()
        vv = float   ( v ) 
        if IM != float ( v ) and isfinite ( vv ) : 
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = "value/min" , field , '' if not n else fmt_factor % n 
            rows.append ( row )
        
    if 1 <= order and 1 <= size and obj.ok () and hasattr ( obj , 'max' ) :
        
        v  = obj.max ()
        vv = float   ( v ) 
        if IM != float ( v ) and isfinite ( vv ) : 
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = "value/max" , field , '' if not n else fmt_factor % n 
            rows.append ( row )
                      
    if 1 <= order and 1 <= size and obj.ok () and hasattr ( obj , 'mean' )  : 
        
        v  = obj.mean ()
        vv = float   ( v ) 
        if IM != float ( v ) and isfinite ( vv ) : 
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = "mean" , field , '' if not n else fmt_factor % n 
            rows.append ( row )
            
    if 2 <= order and 2 <= size and obj.ok() :

        if hasattr ( obj , 'rms' ) : 
            v  = obj.rms (    )
            vv = float   ( v  )                             
            if isfinite  ( vv ) and IM != vv and 0 <= vv :                    
                if isinstance ( v , VE ) : field , n = v.pretty_print () 
                else                     : field , n = pretty_float ( v )
                row = "rms"  , field , '' if not n else fmt_factor % n
                rows.append ( row )
                
        if hasattr ( obj , 'variance' ) :
            v  = obj.variance (    )
            vv = float   ( v  )                             
            if isfinite  ( vv ) and IM != vv and 0 <= vv :                    
                if isinstance ( v , VE ) : field , n = v.pretty_print ()
                else                     : field , n = pretty_float ( v )
                row = "variance"  , field , '' if not n else fmt_factor % n 
                rows.append ( row )
                
    if 3 <= order and 3 <= size and obj.ok () and hasattr  ( obj , 'skewness' ) :
        
        v  = obj.skewness ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != vv :
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = "skewness"  , field , '' if not n else fmt_factor % n
            rows.append ( row )
            
    if 4 <= order and 4 <= size and  obj.ok () and hasattr ( obj , 'kurtosis' ) : 
        
        v  = obj.kurtosis ()
        vv = float   ( v )                                     
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = "kurtosis (excess)"  , field , '' if not n else fmt_factor % n  
            rows.append ( row )
            
    fmt = '[%d]'
    if   order <   10 : fmt = '[%d]'
    elif order <  100 : fmt = '[%02d]'
    elif order < 1000 : fmt = '[%03d]'
    
    if 2 <= order and 2 <= size and obj.ok () and hasattr ( obj , 'unbiased_2nd' ) :
        
        v  = obj.unbiased_2nd ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = v.pretty_print ()  
            else                     : field , n = pretty_float ( v )
            row = ( "M" + fmt + "/unb" ) % 2 , field , '' if not n else fmt_factor % n 
            rows.append ( row )

    if 3 <= order and 3 <= size and obj.ok () and hasattr ( obj , 'unbiased_3rd' ) :
        
        v  = obj.unbiased_3rd ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = ( "M" + fmt + "/unb" ) % 3 , field , '' if not n else fmt_factor % n 
            rows.append ( row )

    if 4 <= order and 4 <= size and obj.ok () and hasattr ( obj , 'unbiased_4th' ) :
        
        v  = obj.unbiased_4th ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = ( "M" + fmt + "/unb" ) % 4 , field , '' if not n else fmt_factor % n
            rows.append ( row )
        
    if 5 <= order and 5 <= size and obj.ok () and hasattr ( obj , 'unbiased_5th' ) : 

        v  = obj.unbiased_5th ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = ( "M" + fmt + "/unb" ) % 5 , field , '' if not n else fmt_factor % n  
            rows.append ( row )

    for i in range ( 2 , order + 1 ) :

        if standard and 2 < i :
    
            v  = obj.std_moment2 ( i )

            vv = float  ( v  )                         
            if isfinite ( vv ) and IM != float ( v ) :                
                if isinstance ( v , VE ) : field , n = v.pretty_print () 
                else                     : field , n = pretty_float ( v )
                if 4 == i : 
                    row = "kurtosis (classic)" , field , '' if not n else fmt_factor % n
                elif 5  == i : 
                    row = "hyperskewness" , field , '' if not n else fmt_factor % n
                elif 6 == i :
                    row = "hypertailedness" , field , '' if not n else fmt_factor % n
                else  :
                    row = ( "M" + fmt + "/std" )  % i , field , '' if not n else fmt_factor % n
                rows.append ( row )
                
        elif not standard :

            v  = obj.central_moment ( i )
            vv = float   ( v )                         
            if isfinite ( vv ) and IM != float ( v ) :                
                if isinstance ( v , VE ) : field , n = v.pretty_print () 
                else                     : field , n = pretty_float ( v )
                row = ( "M" + fmt + "" )  % i  , field , '' if not n else fmt_factor % n  
                rows.append ( row )

    if cumulants and 1 <= order and 1 <= size and obj.ok () and hasattr ( obj , 'cumulant_1st' ) :

        v  = obj.cumulant_1st ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = ( "K" + fmt + "/unb" ) % 1 , field , '' if not n else fmt_factor % n 
            rows.append ( row )

    if cumulants and 2 <= order and 2 <= size and obj.ok () and hasattr ( obj , 'cumulant_2nd' ) :

        v  = obj.cumulant_2nd ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = ( "K" + fmt + "/unb" ) % 2 , field , '' if not n else fmt_factor % n 
            rows.append ( row )

    if cumulants and 3 <= order and 3 <= size and obj.ok () and hasattr ( obj , 'cumulant_3rd' ) :

        v  = obj.cumulant_3rd ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = ( "K" + fmt + "/unb" ) % 3 , field , '' if not n else fmt_factor % n
            rows.append ( row )

    if cumulants and 4 <= order and 4 <= size and obj.ok () and hasattr ( obj , 'cumulant_4th' ) :

        v  = obj.cumulant_4th ()
        vv = float   ( v )                         
        if isfinite ( vv ) and IM != float ( v ) :
            if isinstance ( v , VE ) : field , n = v.pretty_print () 
            else                     : field , n = pretty_float ( v )
            row = ( "K" + fmt + "/unb" ) % 4 , field , '' if not n else fmt_factor % n  
            rows.append ( row )
            
    if cumulants and 5 <= order and obj.ok () and hasattr ( obj , 'cumulant' ) :
        kmax = min ( order , 10 )        
        if 10 <= kmax : fmt = '[%02d]'
        for k in range ( 5 , kmax + 1 ) :            
            v  = obj.cumulant ( k )
            vv = float   ( v )                         
            if isfinite ( vv ) and IM != float ( v ) :
                if isinstance ( v , VE ) : field , n = v.pretty_print () 
                else                     : field , n = pretty_float ( v )
                row = ( "K" + fmt + '' ) % k , field , '' if not n else fmt_factor  % n
                rows.append ( row )
                
    rows = tuple ( [ ( '' , 'value' , 'factor' ) ] +  rows  )
        
    import ostap.logger.table as T
    rows = T.remove_empty_columns ( rows )
    
    tit = title 
    if not title :
        if   isinstance ( obj , Ostap.Math. Moment ) : tit =  'Moment_[%d]' % obj.order 
        elif isinstance ( obj , Ostap.Math.WMoment ) : tit = 'WMoment_[%d]' % obj.order 
    
    return T.table ( rows  , title = tit , prefix = prefix , style = style )

Ostap.Math.Moment.unbiased_2nd   = _om_u2nd 
Ostap.Math.Moment.unbiased_3rd   = _om_u3rd
Ostap.Math.Moment.unbiased_4th   = _om_u4th 
Ostap.Math.Moment.unbiased_5th   = _om_u5th 

Ostap.Math.Moment.mean               = _om_mean    
Ostap.Math.Moment.rms                = _om_rms 
Ostap.Math.Moment.variance           = _om_variance
Ostap.Math.Moment.skewness           = _om_skewness
Ostap.Math.Moment.kurtosis           = _om_kurtosis
Ostap.Math.Moment.cmoment            = _om_cm2
Ostap.Math.Moment.central_moment     = _om_cm2
Ostap.Math.Moment.centralized_moment = _om_cm4
Ostap.Math.Moment.std_moment2        = _om_std
Ostap.Math.Moment.table              = _om_table

Ostap.Math.Moment.cumulant_1st   = _om_cumulant_1st ## 1st cumulant is a mens
Ostap.Math.Moment.cumulant_2nd   = _om_cumulant_2nd ## 2nd cumulant is a variance 
Ostap.Math.Moment.cumulant_3rd   = _om_cumulant_3rd ## 3nd cumulant is a 3rd central moment 
Ostap.Math.Moment.cumulant_4th   = _om_cumulant_4th ## 4th cumulant is different 
Ostap.Math.Moment.cumulant       = _om_cumulant_  

Ostap.Math.WMoment.mean           = _om_mean    
Ostap.Math.WMoment.rms            = _om_rms 
Ostap.Math.WMoment.variance       = _om_variance
Ostap.Math.WMoment.skewness       = _om_skewness
Ostap.Math.WMoment.kurtosis       = _om_kurtosis
Ostap.Math.WMoment.cmoment        = _om_cm3
Ostap.Math.WMoment.central_moment = _om_cm3
Ostap.Math.WMoment.std_moment2    = _om_std
Ostap.Math.WMoment.table          = _om_table

Ostap.Math.Moment.minmax           = _om_minmax
Ostap.Math.WMoment.minmax          = _om_minmax

for t in ( Ostap.Math.WMoment ,
           Ostap.Math. Moment ) :
    t.__str__  = _om_table
    t.__repr__ = _om_table

M0  = Ostap.Math.Moment_(0)
M1  = Ostap.Math.Moment_(1)
WM0 = Ostap.Math.WMoment_(0)
WM1 = Ostap.Math.WMoment_(1)
for m in ( M0 , M1 , WM0 , WM1 ) :
    m.mean  = _om_mean
    m.table = _om_table

if not hasattr (  M0 , 'order' ) :  M0.order = 0
if not hasattr (  M1 , 'order' ) :  M1.order = 1
if not hasattr ( WM0 , 'order' ) : WM0.order = 0
if not hasattr ( WM1 , 'order' ) : WM1.order = 1



# ============================================================
## equality for two counters  
def _mom0_eq_  ( cnt , another ) :
    """ Equality for two Ostap.Math.Momemnt_<0> counters
    """
    ## The same base type
    if not isinstance ( another , Ostap.Math.Moment ) : return NotImplemented
    ## order is defined 
    if not hasattr    ( another , 'order'           ) : return NotImplemented
    ## the same order 
    if cnt.order  !=  another.order                   : return NotImplemented
    ### the same size 
    return cnt.size() == another.size ()
# ============================================================
## equality for two counters  
def _mom1_eq_  ( cnt , another ) :
    """ Equality for two Ostap.Math.Momemnt_<1> counters
    """
    ## The same base type
    if not isinstance ( another , Ostap.Math.Moment ) : return NotImplemented
    ## order is defined 
    if not hasattr    ( another , 'order'           ) : return NotImplemented
    ## the same order 
    if cnt.order  !=  another.order                   : return NotImplemented
    return isequal ( cnt.mu() , another.mu() ) and  cnt.previous() == another.previous() 
# ============================================================
## equality for two counters  
def _momN_eq_  ( cnt , another ) :
    """ Equality for two Ostap.Math.Moment_<N> counters
    """
    ## The same base type
    if not isinstance ( another , Ostap.Math.Moment ) : return NotImplemented
    ## order is defined 
    if not hasattr    ( another , 'order'           ) : return NotImplemented
    ## the same order 
    if cnt.order  !=  another.order                   : return NotImplemented
    ## 
    return isequal ( cnt.M () , another.M () ) and cnt.previous() == another.previous() 

Ostap.Math.Moment . __eq__ = _momN_eq_ 
M0                . __eq__ = _mom0_eq_
M1                . __eq__ = _mom1_eq_

# ============================================================
## equality for two counters  
def _wmom0_eq_  ( cnt , another ) :
    """ Equality for two Ostap.Math.Momemnt_<0> counters
    """
    ## The same base type
    if not isinstance ( another , Ostap.Math.WMoment ) : return NotImplemented
    ## order is defined 
    if not hasattr    ( another , 'order'            ) : return NotImplemented
    ## the same order 
    if cnt.order  !=  another.order                    : return NotImplemented
    ### the same size 
    return \
        cnt.size() == another.size ()         and \
        isequal ( cnt.w  () , another.w  () ) and \
        isequal ( cnt.w2 () , another.w2 () ) 


# ============================================================
## equality for two counters  
def _wmom1_eq_  ( cnt , another ) :
    """ Equality for two Ostap.Math.WMomemnt_<1> counters
    """
    ## The same base type
    if not isinstance ( another , Ostap.Math.WMoment ) : return NotImplemented
    ## order is defined 
    if not hasattr    ( another , 'order'            ) : return NotImplemented
    ## the same order 
    if cnt.order  !=  another.order                    : return NotImplemented
    return isequal ( cnt.mu() , another.mu() ) and  cnt.previous() == another.previous() 
# ============================================================
## equality for two counters  
def _wmomN_eq_  ( cnt , another ) :
    """ Equality for two Ostap.Math.WMoment_<N> counters
    """
    ## The same base type
    if not isinstance ( another , Ostap.Math.WMoment ) : return NotImplemented
    ## order is defined 
    if not hasattr    ( another , 'order'            ) : return NotImplemented
    ## the same order 
    if cnt.order  !=  another.order                    : return NotImplemented
    ## 
    return isequal ( cnt.M () , another.M () ) and cnt.previous() == another.previous() 


Ostap.Math.WMoment . __eq__ = _wmomN_eq_ 
WM0                . __eq__ = _wmom0_eq_
WM1                . __eq__ = _wmom1_eq_
    

# ==========================================================
## REDUCE 
from ostap.math.reduce import root_factory 
# ==========================================================
## Redude Ostap::Math::Moment_<0>
def _mom0_reduce_ ( cnt ) :
    """ Reduce Ostap::Math::Moment_<0>
    """
    return root_factory , ( type ( cnt ), cnt.size () )
# =========================================================
## Redude Ostap::Math::Moment_<1>
def _mom1_reduce_ ( cnt ) :
    """ Reduce Ostap::Math::Moment_<1>
    """
    return root_factory , ( type ( cnt )    ,
                            cnt.previous () ,
                            cnt.mu       () ,
                            cnt.min      () ,
                            cnt.max      () )
# ========================================================
## Redude Ostap::Math::Moment_<N>
def _momN_reduce_ ( cnt ) :
    """ Reduce Ostap::Math::Moment_<N>
    """
    return root_factory , ( type ( cnt ) , cnt.previous () , cnt.moment () ) 

Ostap.Math.Moment . __reduce__ = _momN_reduce_
M0                . __reduce__ = _mom0_reduce_
M1                . __reduce__ = _mom1_reduce_

# ==========================================================
## Redude Ostap::Math::WMoment_<0>
def _wmom0_reduce_ ( cnt ) :
    """ Reduce Ostap::Math::WMoment_<0>
    """
    return root_factory , ( type ( cnt ) ,
                            cnt.size  () ,
                            cnt.w     () ,
                            cnt.w2    () ,
                            cnt.wmin  () ,
                            cnt.wmax  () )
# =========================================================
## Redude Ostap::Math::WMoment_<1>
def _wmom1_reduce_ ( cnt ) :
    """ Reduce Ostap::Math::WMoment_<1>
    """
    return root_factory , ( type ( cnt )    ,
                            cnt.previous () , 
                            cnt.mu       () ,
                            cnt.min      () ,
                            cnt.max      () )

# ========================================================
## Redude Ostap::Math::WMoment_<N>
def _wmomN_reduce_ ( cnt ) :
    """ Reduce Ostap::Math::WMoment_<N>
    """
    return root_factory , ( type ( cnt ) , cnt.previous () , cnt.moment () ) 

Ostap.Math.WMoment . __reduce__ = _wmomN_reduce_
WM0                . __reduce__ = _wmom0_reduce_
WM1                . __reduce__ = _wmom1_reduce_

# =============================================================================
## serialization of Harmonic&Geometric means
#  @see Ostap::Math::GeometricMean
#  @see Ostap::Math::HarmonicMean
#  @see Ostap::Math::WGeometricMean
#  @see Ostap::Math::WHarmonicMean
def _mn_reduce_ ( cnt ) :
    """ Serialization of Harmonic&Geometric means
    - see Ostap::Math::GeometricMean
    - see Ostap::Math::HarmonicMean
    - see Ostap::Math::WGeometricMean
    - see Ostap::Math::WHarmonicMean
    """
    return root_factory , ( type ( cnt ) , cnt.counter() )

Ostap.Math.GeometricMean   . __reduce__ = _mn_reduce_
Ostap.Math.HarmonicMean    . __reduce__ = _mn_reduce_
Ostap.Math.ArithmeticMean  . __reduce__ = _mn_reduce_
Ostap.Math.WGeometricMean  . __reduce__ = _mn_reduce_
Ostap.Math.WHarmonicMean   . __reduce__ = _mn_reduce_
Ostap.Math.WArithmeticMean . __reduce__ = _mn_reduce_


# =============================================================================
## serialization of Power means
#  @see Ostap::Math::PoweMean
#  @see Ostap::Math::WPoweMean
def _pm_reduce_ ( cnt ) :
    """ Serialization of Power means
    - see Ostap::Math::PowerMean
    - see Ostap::Math::WPowerMean
    """
    return root_factory , ( type ( cnt ) , cnt.p() , cnt.counter() )

Ostap.Math.PowerMean   . __reduce__ = _pm_reduce_
Ostap.Math.WPowerMean  . __reduce__ = _pm_reduce_

# =============================================================================
## serialization of Lehmermeans
#  @see Ostap::Math::LehmerMean
#  @see Ostap::Math::WLehmerMean
def _lm_reduce_ ( cnt ) :
    """ Serialization of Lehmer means
    - see Ostap::Math::LehmerMean
    - see Ostap::Math::WLehmerMean
    """
    return root_factory , ( type ( cnt ) , cnt.p() , cnt.counter1() , cnt.counter2() )

Ostap.Math.LehmerMean   . __reduce__ = _lm_reduce_
Ostap.Math.WLehmerMean  . __reduce__ = _lm_reduce_


# ==============================================================================
## equality for the mean
#  @see Ostap::Math::GeometricMean
#  @see Ostap::Math::HarmonicMean
#  @see Ostap::Math::ArithmeticMean
#  @see Ostap::Math::WGeometricMean
#  @see Ostap::Math::WHarmonicMean
#  @see Ostap::Math::WArithmeticMean
def _mn_eq_ ( cnt , another ) :
    """ Equality Harmonic&Geometric means
    - see Ostap::Math::GeometricMean
    - see Ostap::Math::HarmonicMean
    - see Ostap::Math::ArithmeticMean
    - see Ostap::Math::WGeometricMean
    - see Ostap::Math::WHarmonicMean
    - see Ostap::Math::WArithmeticMean
    """
    if not type ( cnt ) is type ( another ) : return NotImplemented
    return isequal ( cnt.value() , another.value() ) and \
        cnt.counter() == another.counter()

Ostap.Math.GeometricMean   . __eq__ = _mn_eq_
Ostap.Math.HarmonicMean    . __eq__ = _mn_eq_
Ostap.Math.ArithmeticMean  . __eq__ = _mn_eq_
Ostap.Math.WGeometricMean  . __eq__ = _mn_eq_
Ostap.Math.WHarmonicMean   . __eq__ = _mn_eq_
Ostap.Math.WArithmeticMean . __eq__ = _mn_eq_

_decorated_classes = (
    Ostap.Math.Moment          ,
    Ostap.Math.WMoment         ,
    Ostap.Math.Moment_(0)      , 
    Ostap.Math.Moment_(1)      ,
    Ostap.Math.Moment_(0)      , 
    Ostap.Math.Moment_(1)      ,
    ##
    Ostap.Math.GeometricMean   , 
    Ostap.Math.HarmonicMean    , 
    Ostap.Math.ArithmeticMean  , 
    Ostap.Math.WGeometricMean  ,
    Ostap.Math.WHarmonicMean   ,    
    Ostap.Math.WArithmeticMean , 
)

_new_methods_ = (
    ##
    Ostap.Math.Moment.unbiased_2nd   , 
    Ostap.Math.Moment.unbiased_3rd   ,
    Ostap.Math.Moment.unbiased_4th   ,
    Ostap.Math.Moment.unbiased_5th   ,
    ##
    Ostap.Math.Moment.cumulant_1st   , 
    Ostap.Math.Moment.cumulant_2nd   , 
    Ostap.Math.Moment.cumulant_3rd   , 
    Ostap.Math.Moment.cumulant_4th   ,
    ##
    Ostap.Math.Moment.cumulant       ,
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
    Ostap.Math.Moment.minmax          , 
    Ostap.Math.WMoment.minmax         , 
)
                          
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
