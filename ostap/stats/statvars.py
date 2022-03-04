#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file statvars.py 
#  Functions to collect statistics for trees and  datasets
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06  
# =============================================================================
"""Functions to collect statistics for trees and  datasets
- data_get_moment      - calculate the moment 
- data_moment          - get the moment            (with uncertainty)
- data_central_moment  - get the central moment    (with uncertainty)
- data_mean            - get the mean              (with uncertainty)
- data_variance        - get the variance          (with uncertainty)
- data_dispersion      - get the dispersion        (with uncertainty)
- data_rms             - get the RMS               (with uncertainty)
- data_skewness        - get the skewness          (with uncertainty)
- data_kurtosis        - get the (excess) kurtosis (with uncertainty)
- data_quantile        - get the quantile 
- data_median          - get the median
- data_quantiles       - the quantiles  
- data_interval        - get the interval 
- data_terciles        - get two terciles 
- data_quartiles       - get three quartiles 
- data_quintiles       - get four  quintiles 
- data_deciles         - get nine  deciles
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-06"
__all__     = (
    'data_get_moment'     , ##  calculate the moment 
    'data_moment'         , ## get the moment            (with uncertainty)
    'data_central_moment' , ## get the central moment    (with uncertainty)
    'data_mean'           , ## get the mean              (with uncertainty)
    'data_variance'       , ## get the variance          (with uncertainty)
    'data_dispersion'     , ## get the dispersion        (with uncertainty)
    'data_rms'            , ## get the RMS               (with uncertainty)
    'data_skewness'       , ## get the skewness          (with uncertainty)
    'data_kurtosis'       , ## get the (excess) kurtosis (with uncertainty)
    'data_quantile'       , ## get the quantile 
    'data_median'         , ## get the median
    'data_quantiles'      , ## get the quantiles  
    'data_interval'       , ## get the interval 
    'data_terciles'       , ## get two terciles 
    'data_quartiles'      , ## get three quartiles 
    'data_quintiles'      , ## get four  quintiles 
    'data_deciles'        , ## get nine  deciles
    'data_decorate'       , ## technical function to decorate the class
    )
# =============================================================================
from   builtins             import range
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.statvars' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
from   ostap.core.core import Ostap
import ostap.stats.moment 
# =============================================================================
StatVar = Ostap.StatVar 
# =============================================================================
## @var QEXACT
#  use it as threshold for exact/slow vs approximate/fast quanitle calcualtion 
QEXACT = 10000
# =============================================================================
## get the moment of order 'order' relative to 'center'
#  @code
#  data =  ...
#  print data_get_moment ( data , 3 , 0.0 , 'mass' , 'pt>1' ) 
#  print data.get_moment (        3 , 0.0 , 'mass' , 'pt>1' ) ## ditto
#  @endcode
#  @see Ostap::StatVar::get_moment
def data_get_moment ( data , order , center , expression , cuts = '' , *args ) :
    """Get the moment of order 'order' relative to 'center'
    >>> data =  ...
    >>> print data_get_moment ( data , 3 , 0.0 , 'mass' , 'pt>1' )
    >>> print data.get_moment (        3 , 0.0 , 'mass' , 'pt>1' ) ## ditto
    - see Ostap::StatVar::get_moment
    """
    assert isinstance ( order  , int ) and 0<= order , 'Invalid order  %s'  % order
    return StatVar.get_moment ( data       ,
                                order      ,
                                expression ,
                                center     ,
                                cuts       , *args )

# =============================================================================
## get the moment (with uncertainty) of order 'order' 
#  @code
#  data =  ...
#  print data_moment ( data , 3 , 'mass' , 'pt>1' ) 
#  print data.moment (        3 , 'mass' , 'pt>1' ) ## ditto
#  @endcode
#  @see Ostap::StatVar::moment
def  data_moment ( data , order , expression , cuts  = ''  , *args ) :
    """Get the moment of order 'order' relative to 'center'
    >>> data =  ...
    >>> print data_moment ( data ,  3 , 'mass' , 'pt>1' ) 
    >>> print data.moment (         3 , 'mass' , 'pt>1' ) ## ditto
    - see Ostap::StatVar::moment
    """
    assert isinstance ( order  , int ) and 0<= order , 'Invalid order  %s'  % order
    return StatVar.moment ( data       ,
                            order      ,
                            expression ,
                            cuts       , *args )

# =============================================================================
## get the central moment (with uncertainty) of order 'order' 
#  @code
#  data =  ...
#  print data_central_moment ( data , 3 , 'mass' , 'pt>1' ) 
#  print data.central_moment (        3 , 'mass' , 'pt>1' ) ## ditto
#  @endcode
#  - The moments are calculated with  O(1/n) precision 
#  - For 3rd and  4th moments,  explicit  correction is applied
#  - Uncertainty is calculated with O(1/n^2) precision
#  @see Ostap::StatVar::central_moment
def  data_central_moment ( data , order , expression , cuts  = '' , *args ) :
    """Get the moment of order 'order' relative to 'center'
    >>> data =  ...
    >>> print data_central_moment ( data , 3 , 'mass' , 'pt>1' ) 
    >>> print data.central_moment (        3 , 'mass' , 'pt>1' ) ## ditto
    - The moments are calculated with  O(1/n) precision 
    - For 3rd and  4th moments,  explicit  correction is applied
    - Uncertainty is calculated with O(1/n^2) precision
    - see Ostap::StatVar::central_moment
    """
    assert isinstance ( order  , int ) and 0<= order , 'Invalid order  %s'  % order
    return StatVar.central_moment ( data       ,
                                    order      ,
                                    expression ,
                                    cuts       , *args )

# =============================================================================
## get the  skewness (with uncertainty)
#  @code
#  data =  ...
#  print data_skewness ( data , 'mass' , 'pt>1' ) 
#  print data.skewness (        'mass' , 'pt>1' ) ## ditto
#  @endcode
#  @see Ostap::StatVar::skewness
def data_skewness ( data , expression , cuts  = '' , *args ) :
    """Get the  skewness
    >>> data =  ...
    >>> print data_skewness ( data , 'mass' , 'pt>1' ) 
    >>> print data.skewness (        'mass' , 'pt>1' ) ## ditto
    - see Ostap::StatVar::skewness
    """
    return StatVar.skewness ( data , expression , cuts , *args )

# =============================================================================
## get the (excess) kurtosis (with uncertainty)
#  @code
#  data =  ...
#  print data_kustosis ( data , 'mass' , 'pt>1' ) 
#  print data.kustosis (        'mass' , 'pt>1' ) ## ditto
#  @endcode
#  @see Ostap::StatVar::kurtosis
def data_kurtosis ( data , expression , cuts  = '' , *args ) :
    """Get the kurtosis
    >>> data =  ...
    >>> print data_kurtosis ( data , 'mass' , 'pt>1' ) 
    >>> print data.kurtosis (        'mass' , 'pt>1' ) ## ditto
    - see Ostap::StatVar::kurtosis
    """
    return StatVar.kurtosis ( data , expression , cuts , *args )

# =============================================================================
## get the  quantile 
#  @code
#  data =  ...
#  print data_quantile ( data , 0.10 , 'mass' , 'pt>1' ) 
#  print data.quantile (        0.10 , 'mass' , 'pt>1' ) ## ditto
#  @endcode
#  @see Ostap::StatVar::quantile
#  @see Ostap::StatVar::p2quantile
def  data_quantile ( data , q , expression , cuts  = '' , exact = QEXACT , *args ) :
    """Get the quantile
    >>> data =  ...
    >>> print data_quantile ( data , 0.1 , 'mass' , 'pt>1' ) 
    >>> print data.quantile (        0.1 , 'mass' , 'pt>1' ) ## ditto
    - see Ostap.StatVar.quantile
    - see Ostap.StatVar.p2quantile
    """
    assert isinstance ( q , float ) and 0 < q < 1 , 'Invalid quantile:%s' % q

    if   exact is True  :
        ##  exact slow algorithm 
        qn = StatVar.  quantile ( data , q , expression , cuts , *args )

    elif exact is False :
        ## approximate fast formula 
        qn = StatVar.p2quantile ( data , q , expression , cuts , *args )

    elif isinstance ( exact , int ) and len ( data ) <= exact  :
        ## exact slow algorithm 
        qn = StatVar.  quantile ( data , q , expression , cuts , *args )

    else :
        ## approximate fast algorithm 
        qn = StatVar.p2quantile ( data , q , expression , cuts , *args )

    return qn

# =============================================================================
## get the interval
#  @code
#  data =  ...
#  print data_interval ( data , 0.05, 0.95 , 'mass' , 'pt>1' ) ##   get 90% interval
#  print data.interval (        0.05, 0.95 , 'mass' , 'pt>1' ) ##   get 90% interval
#  @endcode
#  @see Ostap::StatVar::interval
#  @see Ostap::StatVar::p2interval
def data_interval ( data , qmin ,  qmax , expression , cuts = '' , exact = QEXACT , *args ) :
    """Get the interval 
    >>> data =  ...
    >>> print data_interval ( data , 0.05 , 0.95 , 'mass' , 'pt>1' ) ## get 90% interval
    >>> print data.interval (        0.05 , 0.95 , 'mass' , 'pt>1' ) ## get 90% interval
    - see Ostap::StatVar::interval
    """
    assert isinstance ( qmin , float ) and 0 < qmin < 1 , 'Invalid quantile-1:%s' % qmin 
    assert isinstance ( qmax , float ) and 0 < qmax < 1 , 'Invalid quantile-2:%s' % qmax

    qmin, qmax = min ( qmin ,  qmax ) , max  ( qmin, qmax  )
    
    if   exact is True  :
        ##  exact slow algorithm 
        rn = StatVar.  interval ( data , qmin , qmax  , expression , cuts , *args )

    elif exact is False :
        ## approximate fast formula 
        rn = StatVar.p2interval ( data , qmin , qmax , expression , cuts , *args )

    elif isinstance ( exact , int ) and len ( data ) <= exact  :
        ## exact slow algorithm 
        rn = StatVar.  interval ( data , qmin , qmax  , expression , cuts , *args )

    else :
        ## approximate fast formula 
        rn = StatVar.p2interval ( data , qmin , qmax , expression , cuts , *args )

    ## x
    return rn

# =============================================================================
## Get the median
#  @code
#  data =  ...
#  print data_median ( data , 'mass' , 'pt>1' )
#  print data.median (        'mass' , 'pt>1' ) ## ditto
#  @endcode 
#  @see Ostap::StatVar::quantile
def data_median ( data , expression , cuts = '' , exact = QEXACT , *args ) :
    """Get the median
    >>> data =  ...
    >>> print data_median ( data , 'mass' , 'pt>1' ) 
    >>> print data.median (        'mass' , 'pt>1' ) ##  ditto
    - see Ostap::StatVar::quantile
    """
    return data_quantile ( data , 0.5 , expression , cuts  , exact , *args ) 

# =============================================================================
## get the  quantiles 
#  @code
#  data =  ...
#  print data_quantiles ( data , [0.1,0.3,0.5] , 'mass' , 'pt>1' ) 
#  print data_quantiles ( data , 0.12  , 'mass' , 'pt>1' ) 
#  print data_quantiles ( data , 20    , 'mass' , 'pt>1' )
#  print data.quantiles (        [0.1,0.3,0.5] , 'mass' , 'pt>1' ) 
#  print data.quantiles (        0.12  , 'mass' , 'pt>1' ) 
#  print data.quantiles (        20    , 'mass' , 'pt>1' ) 
#  @endcode
#  @see Ostap::StatVar::quantile
def data_quantiles ( data , quantiles , expression , cuts  = '' , exact = QEXACT , *args ) :
    """Get the quantiles
    >>> data =  ...
    >>> print data_quantiles ( data , 0.1       , 'mass' , 'pt>1' ) ## quantile 
    >>> print data_quantiles ( data , (0.1,0.5) , 'mass' , 'pt>1' )
    >>> print data_quantiles ( data , 10        , 'mass' , 'pt>1' ) ## deciles!     
    >>> print data.quantiles (        0.1       , 'mass' , 'pt>1' ) ## quantile 
    >>> print data.quantiles (        (0.1,0.5) , 'mass' , 'pt>1' )
    >>> print data.quantiles (        10        , 'mass' , 'pt>1' ) ## deciles!     
    - see Ostap::StatVar::quantile
    """
    if   isinstance ( quantiles , float ) and 0 < quantiles < 1 : 
        quantiles = [ quantiles ]
    elif isinstance ( quantiles , int   ) and 1 < quantiles     :
        N         = quantiles 
        quantiles =  ( float ( i ) / N for i in range ( 1 , N ) )

    qq = [] 
    for q in quantiles :
        assert isinstance ( q , float ) and 0 < q < 1 , 'Invalid quantile:%s' % q
        qq.append ( q )
    qq.sort ()

    from ostap.math.base import doubles
    qqq = doubles ( qq )

    if   exact is True  :
        ##  exact slow algorithm 
        qn = StatVar.  quantiles ( data , qqq , expression , cuts , *args )

    elif exact is False :
        ## approximate fast formula 
        qn = StatVar.p2quantiles ( data , qqq , expression , cuts , *args )

    elif isinstance ( exact , int ) and len ( data ) <= exact  :
        ## exact slow algorithm 
        qn = StatVar.  quantiles ( data , qqq , expression , cuts , *args )

    else :
        ## approximate fast algorithm 
        qn = StatVar.p2quantiles ( data , qqq , expression , cuts , *args )

    return qn 

# =============================================================================
## Get the terciles 
#  @code
#  data =  ...
#  print data_terciles( data , 'mass' , 'pt>1' )
#  print data.terciles(        'mass' , 'pt>1' )   ##  ditto
#  @endcode 
#  @see Ostap::StatVar::quantile
def data_terciles ( data , expression , cuts = '' , exact = QEXACT , *args ) :
    """Get the terciles
    >>> data =  ...
    >>> print data_terciles ( data , 'mass' , 'pt>1' ) 
    >>> print data.terciles (        'mass' , 'pt>1' ) ## ditto
    - see Ostap::StatVar::quantile
    """
    return data_quantiles ( data  , 3 , expression  , cuts , exact , *args ) 

# =============================================================================
## Get the quartiles 
#  @code
#  data =  ...
#  print data_quartiles( data , 'mass' , 'pt>1' )
#  print data.quartiles(        'mass' , 'pt>1' ) ##  ditto
#  @endcode 
#  @see Ostap::StatVar::quantile
def data_quartiles ( data , expression , cuts = '' , exact = QEXACT , *args ) :
    """Get the quartiles
    >>> data =  ...
    >>> print data_quartiles ( data , 'mass' , 'pt>1' ) 
    >>> print data.quartiles (        'mass' , 'pt>1' ) ##  ditto
    - see Ostap::StatVar::quantile
    """
    return data_quantiles ( data  , 4 , expression  , cuts , exact , *args ) 

# =============================================================================
## Get the quintiles
#  @code
#  data =  ...
#  print data_quintiles( data , 'mass' , 'pt>1' )
#  print data.quintiles(        'mass' , 'pt>1' ) ## ditto
#  @endcode 
#  @see Ostap::StatVar::quantile
def data_quintiles ( data , expression , cuts = '' , exact = QEXACT , *args ) :
    """Get the quartiles
    >>> data =  ...
    >>> print data.quintiles ( 'mass' , 'pt>1' ) 
    >>> print data.quintiles ( 'mass' , 'pt>1' ) ## ditto
    - see Ostap::StatVar::quantile
    """
    return data_quantiles ( data  , 5 , expression  , cuts , exact , *args ) 

# =============================================================================
## Get the deciles 
#  @code
#  data =  ...
#  print data_deciles( data , 'mass' , 'pt>1' )
#  print data.deciles(        'mass' , 'pt>1' ) ## ditto
#  @endcode 
#  @see Ostap::StatVar::quantile
def data_deciles ( data , expression , cuts = '' , exact = QEXACT , *args ) :
    """Get the deciles
    >>> data =  ...
    >>> print data_deciles ( data , 'mass' , 'pt>1' ) 
    >>> print data.deciles (        'mass' , 'pt>1' ) ## ditto
    - see Ostap::StatVar::quantile
    """
    return data_quantiles ( data  , 10 , expression  , cuts , exact , *args ) 

# =============================================================================
## Get the mean (with uncertainty):
#  @code
#  data = ...
#  data_mean( data , 'mass*mass', 'pt>0')
#  data.mean(        'mass*mass', 'pt>0') ## ditto
#  @endcode 
#  @see Ostap::StatVar::moment
def data_mean ( data  , expression , cuts  = '' , *args ) :
    """Get the   mean (with uncertainty):
    >>> data = ...
    >>> data_mean( data , 'mass*mass', 'pt>0')
    >>> data.mean(        'mass*mass', 'pt>0') ## ditto
    """
    return data_moment ( data , 1 , expression ,  cuts , *args )

# =============================================================================
## Get the variance(with uncertainty):
#  @code
#  data = ...
#  data_variance (  data , 'mass*mass', 'pt>0')
#  data.variance (         'mass*mass', 'pt>0') ## ditto
#  @endcode 
#  @see Ostap::StatVar::central_moment
def data_variance ( data , expression , cuts  = '' , *args ) :
    """Get the variance (with uncertainty):
    >>> data = ...
    >>> data_variance ( data , 'mass*mass', 'pt>0')
    >>> data.variance (        'mass*mass', 'pt>0') ## ditto
    """
    return data_central_moment ( data , 2 , expression , cuts , *args )

# =============================================================================
## Get the dispersion(with uncertainty):
#  @code
#  data = ...
#  data_dispersion (  data , 'mass*mass', 'pt>0')
#  data.dispersion (         'mass*mass', 'pt>0') ## ditto
#  @endcode 
#  @see Ostap::StatVar::central_moment
def data_dispersion ( data , expression , cuts  = '' , *args ) :
    """Get the variance (with uncertainty):
    >>> data = ...
    >>> data_variance ( data , 'mass*mass', 'pt>0')
    >>> data.variance (        'mass*mass', 'pt>0') ## ditto
    """
    return data_variance ( data , expression , cuts , *args ) 

# =============================================================================
## Get the rms(with uncertainty):
#  @code
#  data = ...
#  data_rms(  data , 'mass*mass', 'pt>0')
#  data.rms(         'mass*mass', 'pt>0') ## ditto
#  @endcode 
#  @see Ostap::StatVar::central_moment
def data_rms ( data , expression , cuts  = '' , *args ) :
    """Get the rms (with uncertainty):
    >>> data = ...
    >>> data_rms( data , 'mass*mass', 'pt>0')
    >>> data.rms(        'mass*mass', 'pt>0') ## ditto
    """
    return data_variance ( data ,  expression , cuts , *args ) ** 0.5


data_get_moment      .__doc__ += '\n' + StatVar.get_moment     .__doc__  
data_moment          .__doc__ += '\n' + StatVar.moment         .__doc__
data_central_moment  .__doc__ += '\n' + StatVar.central_moment .__doc__ 
data_skewness        .__doc__ += '\n' + StatVar.skewness       .__doc__ 
data_kurtosis        .__doc__ += '\n' + StatVar.kurtosis       .__doc__ 
data_quantile        .__doc__ += '\n' + StatVar.quantile       .__doc__ 
data_quantiles       .__doc__ += '\n' + StatVar.quantiles      .__doc__ 
data_interval        .__doc__ += '\n' + StatVar.interval       .__doc__ 

def data_decorate ( klass ) :
    
    if hasattr ( klass , 'get_moment'     ) : klass.orig_get_moment     = klass.get_moment
    if hasattr ( klass , 'moment'         ) : klass.orig_moment         = klass.moment
    if hasattr ( klass , 'central_moment' ) : klass.orig_central_moment = klass.central_moment
    if hasattr ( klass , 'mean'           ) : klass.orig_mean           = klass.mean
    if hasattr ( klass , 'variance'       ) : klass.orig_variance       = klass.variance 
    if hasattr ( klass , 'dispersion'     ) : klass.orig_dispersion     = klass.dispersion
    if hasattr ( klass , 'rms'            ) : klass.orig_rms            = klass.rms 
    if hasattr ( klass , 'rms'            ) : klass.orig_rms            = klass.rms 
    if hasattr ( klass , 'skewness'       ) : klass.orig_skewness       = klass.skewness
    if hasattr ( klass , 'kurtosis'       ) : klass.orig_kurtosis       = klass.kurtosis
    
    if hasattr ( klass , 'quantile'       ) : klass.orig_quantile       = klass.quantile
    if hasattr ( klass , 'quantiles'      ) : klass.orig_quantiles      = klass.quantiles
    if hasattr ( klass , 'interval'       ) : klass.orig_interval       = klass.interval    
    if hasattr ( klass , 'median'         ) : klass.orig_median         = klass.median 
    if hasattr ( klass , 'terciles'       ) : klass.orig_terciles       = klass.terciles
    if hasattr ( klass , 'quartiles'      ) : klass.orig_quartiles      = klass.quartiles
    if hasattr ( klass , 'quintiles'      ) : klass.orig_quintiles      = klass.quintiles
    if hasattr ( klass , 'deciles'        ) : klass.orig_deciles        = klass.deciles

    klass.get_moment      = data_get_moment
    klass.moment          = data_moment
    klass.central_moment  = data_central_moment
    klass.mean            = data_mean
    klass.variance        = data_variance 
    klass.dispersion      = data_dispersion
    klass.rms             = data_rms
    klass.skewness        = data_skewness
    klass.kurtosis        = data_kurtosis
    
    klass.quantile        = data_quantile
    klass.quantiles       = data_quantiles
    klass.interval        = data_interval    
    klass.median          = data_median
    klass.terciles        = data_terciles
    klass.quartiles       = data_quartiles
    klass.quintiles       = data_quintiles
    klass.deciles         = data_deciles


# =============================================================================

Quantile  = StatVar.Quantile
Quantiles = StatVar.Quantiles
Interval  = StatVar.Interval 
QInterval = StatVar.QInterval 

def _q_str_  ( o ) : return "Quantile(%.5g,n=%s)"         % ( o.quantile  , o.nevents )
def _qs_str_ ( o ) : return "Quantiles(%s,n=%s)"          % ( o.quantiles , o.nevents )
def _i_str_  ( o ) : return "Interval([%.5g,%.5g])"       % ( o.low       , o.high    )
def _qi_str_ ( o ) : return "QInterval([%.5g,%.5g],n=%s)" % ( o.interval.low  ,
                                                              o.interval.high , o.nevents )
Quantile  .__str__  = _q_str_
Quantile  .__repr__ = _q_str_
Quantiles .__str__  = _qs_str_
Quantiles .__repr__ = _qs_str_
Interval  .__str__  = _i_str_
Interval  .__repr__ = _i_str_
QInterval .__str__  = _qi_str_
QInterval .__repr__ = _qi_str_


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
