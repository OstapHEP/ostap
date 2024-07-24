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
- data_get_stat        - get the momentt-based statistics 
- data_central_moment  - get the central moment    (with uncertainty)
- data_mean            - get the mean              (with uncertainty)
- data_nEff            - get the effective number of entries 
- data_sum             - get the (weigted) sum
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
- data_the_moment      - get the (central) moment   
- data_the_mean        - get mean     via moment 
- data_the_rms         - get rms      via moment 
- data_the_variance    - get variance via moment 
- data_the_skewness    - get skewness via moment 
- data_the_kurtosis    - get kurtosis via moment 
- data_harmonic_mean   - get the (weighted) harmonic  mean 
- data_geometric_mean  - get the (weighte)  geometric mean 
- data_arithmetic_mean - get the (weighted) geometric mean
- data_power_mean      - get the (weighted) power     mean 
- data_lehmer_mean     - get the (weighted) Lehmer    mean 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-06"
__all__     = (
    'data_get_moment'      , ##  calculate the moment 
    'data_moment'          , ## get the moment            (with uncertainty)
    'data_get_stat'        , ## get the momentt-based statistics 
    'data_central_moment'  , ## get the central moment    (with uncertainty)
    'data_nEff'            , ## get the effective number of entries 
    'data_sum'             , ## get the (weigted) sum
    'data_mean'            , ## get the mean              (with uncertainty)
    'data_variance'        , ## get the variance          (with uncertainty)
    'data_dispersion'      , ## get the dispersion        (with uncertainty)
    'data_rms'             , ## get the RMS               (with uncertainty)
    'data_skewness'        , ## get the skewness          (with uncertainty)
    'data_kurtosis'        , ## get the (excess) kurtosis (with uncertainty)
    'data_quantile'        , ## get the quantile 
    'data_median'          , ## get the median
    'data_quantiles'       , ## get the quantiles  
    'data_interval'        , ## get the interval 
    'data_terciles'        , ## get two terciles 
    'data_quartiles'       , ## get three quartiles 
    'data_quintiles'       , ## get four  quintiles 
    'data_deciles'         , ## get nine  deciles
    'data_the_moment'      , ## get the cental  moment             
    'data_the_mean'        , ## get mean     via moment 
    'data_the_rms'         , ## get rms      via moment 
    'data_the_variance'    , ## get variance via moment 
    'data_the_skewness'    , ## get skewness via moment 
    'data_the_kurtosis'    , ## get kurtosis via moment 
    'data_harmonic_mean'   , ## get the (weighted) harmonic   mean 
    'data_geometric_mean'  , ## get the (weighted) geometric  mean 
    'data_arithmetic_mean' , ## get the (weighted) arithmetic mean 
    'data_power_mean'      , ## get the (weighted) power      mean 
    'data_lehmer_mean'     , ## get the (weighted) Lehmer     mean 
    'data_decorate'        , ## technical function to decorate the class
    'expression_types'     , ## valid types for expressions/cuts/weights
)
# =============================================================================
from   builtins               import range
from   ostap.math.base        import isequal, iszero 
from   ostap.core.core        import Ostap, rootException, WSE   
from   ostap.core.ostap_types import ( string_types , integer_types  , 
                                       num_types    , dictlike_types )
from   ostap.trees.cuts       import expression_types, vars_and_cuts
from   ostap.utils.basic      import loop_items
import ostap.stats.moment 
import ostap.logger.table     as     T
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.statvars' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
StatVar = Ostap.StatVar 
# =============================================================================
## @var QEXACT
#  use it as threshold for exact/slow vs approximate/fast quantile calcualtion 
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
    assert isinstance ( order      , integer_types    ) and 0 <= order , 'Invalid order  %s'  % order
    assert isinstance ( center     , num_types        )                , 'Invalid center!'  
    
    assert isinstance ( expression , expression_types ) , 'Invalid type of expression!'
    assert isinstance ( cuts       , expression_types ) , 'Invalid type of cuts/weight!'

    expression = str ( expression )
    cuts       = str ( cuts       ).strip() 
    
    with rootException() : 
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
    assert isinstance ( order      , integer_types    ) and 0 <= order , 'Invalid order  %s'  % order    
    assert isinstance ( expression , expression_types ) , 'Invalid type of expression!'
    assert isinstance ( cuts       , expression_types ) , 'Invalid type of cuts/weight!'
    
    expression = str ( expression )
    cuts       = str ( cuts       ).strip() 
    
    with rootException() : 
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
    assert isinstance ( order      , integer_types    ) and 0 <= order , 'Invalid order  %s'  % order    
    assert isinstance ( expression , expression_types ) , 'Invalid type of expression!'
    assert isinstance ( cuts       , expression_types ) , 'Invalid type of cuts/weight!'
    
    expression = str ( expression )
    cuts       = str ( cuts       ).strip() 

    with rootException() : 
        return StatVar.central_moment ( data       ,
                                        order      ,
                                        expression ,
                                        cuts       , *args )

# ==============================================================================
## Get the (s)Statistic-bases statistics/counter from data
#  @code
#  statobj = Ostap.Math.MinValue()
#  data    = ...
#  result  = data.get_stat( statobj , 'x+y' , 'pt>1' ) 
#  @encode
#  @see Ostap::Math::Statistic
#  @see Ostap::Math::WStatistic
#  @see Ostap::statVar::the_moment
def data_get_stat ( data , statobj , expression , cuts = '' , *args ) :
    """Get the (W)Statistic-based statistics.counters from data  
    >>> data   = ...
    >>> stat   = Ostap.Math.HarmonicMean() 
    >>> result = data.get_stat ( stat , 'x/y+z' , '0<qq' )
    - see Ostap.Math.Moment 
    - see Ostap.Math.WMoment 
    """
    assert isinstance ( statobj , ( Ostap.Math.Statistic , Ostap.Math.WStatistic ) ) ,\
        'get_object: invalid statobj type: %s' % type ( statobj ) 

    assert isinstance ( expression , expression_types ) , 'Invalid type of expression!'
    assert isinstance ( cuts       , expression_types ) , 'Invalid type of cuts/weight!'
    
    expression = str ( expression )
    cuts       = str ( cuts       ).strip()
    
    if isinstance ( data , ROOT.TTree ) and isinstance ( statobj , Ostap.Math.Statistic ) and not cuts : 
        with rootException() :
            sc = StatVar.the_moment ( data , statobj , expression , *args )
            assert sc.isSuccess() , 'Error %s from StatVar::the_moment' % sc 
            return statobj
        
    with rootException() :
        sc = StatVar.the_moment ( data , statobj , expression , cuts , *args )
        assert sc.isSuccess() , 'Error %s from StatVar::the_moment' % sc 
        return statobj

# ==============================================================================
## Get the statistics from data
#  @code
#  data    = ...
#  result  = data_statistics ( data , 'x+y'   , 'pt>1' ) 
#  results = data_statistics ( data , 'x;y;z' , 'pt>1' ) ## result is dictionary 
#  @encode
#  @see Ostap::StatEntity
#  @see Ostap::WStatEntity
#  @see Ostap::StatVar::statVar
def data_statistics ( data , expressions , cuts = '' , *args ) :
    """Get statistics from data 
    >>> data    = ...
    >>> result  = data_statistics ( data , 'x/y+z' , '0<qq' )
    >>> results = data_statistics ( data , 'x/y;z' , '0<qq' ) ## result is dictionary
    - see Ostap.Math.StatEntity
    - see Ostap.Math.WStatEntity
    - see Ostap.StatVar.statVar
    """
    
    ## decode expressions & cuts 
    var_lst, cuts, input_string = vars_and_cuts ( expressions , cuts )

    ## only one name is specified 
    if input_string :         
        with rootException() :
            return StatVar.statVar ( data , var_lst[0] , cuts , *args )

    ## several variables are specofied
    from ostap.core.core import strings 
    names   = strings ( *var_lst )
    results = StatVar.Statistics()
    with rootException() :
        if cuts : StatVar.statVars ( data , results , names , cuts , *args )
        else    : StatVar.statVars ( data , results , names        , *args )
        assert len ( var_lst ) == len ( results ) , \
            'Invalid output from StatVar::statVars!'
        
    result = {}
    for v,r in zip ( var_lst , results ) : result [ v ] = r
    return result 


# ==============================================================================
## Get the covarince from dataxpressio
#  @code
#  data    = ...
#  result  = data_covariance ( data , 'x+y' , 'z'  , '0<u') 
#  @encode
#  @see Ostap::Math::Covariance 
#  @see Ostap::Math::ECovariance 
#  @see Ostap::StatVar::statCov
def data_covariance ( data , expression1 , expression2 , cuts = '' , *args ) :
    """Get the covariance from data 
    >>> data   = ...
    >>> result = data_covariance ( data , 'x/y+z' , 'qq' , 'u>0')
    - see Ostap.Math.Covariance
    - see Ostap.Math.WCovariance 
    - see Ostap.StatVar.statCov 
    """
    assert isinstance ( expression1 , expression_types ) , 'Invalid type of expression1!'
    assert isinstance ( expression2 , expression_types ) , 'Invalid type of expression2!'
    assert isinstance ( cuts        , expression_types ) , 'Invalid type of cuts/weight!'
    
    expression1 = str ( expression1 )
    expression2 = str ( expression2 )
    cuts        = str ( cuts        ).strip() 
    
    with rootException() :
        if isinstance ( data , ROOT.TTree ) and not cuts : 
            return StatVar.statCov ( data , expression1 , expression2 ,        *args )
        else  :
            return StatVar.statCov ( data , expression1 , expression2 , cuts , *args )


# ==============================================================================
## Get the (weighted) sum over the variable(s)
#  @code
#  data    = ...
#  result  = data_sum ( data , 'x+y'   , 'pt>1' ) 
#  results = data_sum ( data , 'x;y;z' , 'pt>1' ) ## result is dictionary 
#  @encode
#  @see Ostap::StatVar::statVar
def data_sum ( data , expressions , cuts = '' , *args ) :
    """Get (weighted) sum over the variables 
    >>> data    = ...
    >>> result  = data_sum ( data , 'x/y+z' , '0<qq' )
    >>> results = data_sum ( data , 'x/y;z' , '0<qq' ) ## result is dictionary
    - see Ostap.StatVar.statVar
    """

    result = data_statistics ( data , expressions , cuts , *args )
    
    if isinstance ( result , dictlike_type ) :
        for k , r in loop_items ( result ) : 
            result [ key ] = VE ( r.sum() , r.sum2() )
    else :
        result = VE ( result.sum() , result.sumw2() )
    
    return result 
        
# ==============================================================================
## Get the effewctoevnumebr of etries in data 
#  @code
#  data    = ...
#  result  = data_nEff ( data , 'x+y' ) 
#  @encode
#  @see Ostap::StatVar::nEff 
def data_nEff ( data , expression = '' , *args ) :
    """Get statistics from data 
    >>> data    = ...
    >>> result  = data_nEff ( data , 'x/y+z' )
    - see Ostap.StatVar.nEff
    """

    assert isinstance ( expression , expression_types ) , 'Invalid type of expression!'
    expression = str ( expression ).strip() 

    return StatVar.nEff ( expression )
        
# =============================================================================
## Get harmonic mean over the data
#  @code
#  data = ...
#  result = data_harmonic_mean ( data , 'pt' , 'eta>0' )
#  @endcode 
#  @see Ostap::Math::Moment
#  @see Ostap::Math::WMoment
#  @see Ostap::Math::HarmonicMean
#  @see Ostap::Math::WHarmonicMean
#  @see Ostap::statVar::the_moment
def data_harmonic_mean ( data , expression , cuts = '' , *args ) :
    """ Get harmonic mean over the data
    >>> data = ...
    >>> result = data_harmonic_mean ( data , 'pt' , 'eta>0' )
    - see `Ostap.Math::Moment`
    - see `Ostap.Math::WMoment`
    - see `Ostap.Math::HarmonicMean`
    - see `Ostap.Math::WHarmonicMean`
    - see `Ostap.statVar.the_moment`
    """    
    if isinstance ( data , ROOT.TTree ) and not cuts : 
        statobj = Ostap.Math.HarmonicMean()
        return data_get_stat ( data , statobj , expression , '' , *args )  
    
    statobj = Ostap.Math.WHarmonicMean()
    return data_get_stat ( data , statobj , expression , cuts , *args )  

# =============================================================================
## Get geometric mean over the data
#  @code
#  data = ...
#  result = data_geometric_mean ( data , 'pt' , 'eta>0' )
#  @endcode 
#  @see Ostap::Math::Moment
#  @see Ostap::Math::GeometricMean
#  @see Ostap::Math::WMoment
#  @see Ostap::Math::WGeometricMean
#  @see Ostap::statVar::the_moment
def data_geometric_mean ( data , expression , cuts = '' , *args ) :
    """ Get geometric mean over the data
    >>> data = ...
    >>> result = data_geometric_mean ( data , 'pt' , 'eta>0' )
    - see `Ostap.Math::Moment`
    - see `Ostap.Math::HarmonicMean`
    - see `Ostap.statVar.the_moment`
    """
    if ( not cuts ) and isinstance ( data , ROOT.TTree ) : 
        statobj = Ostap.Math.GeometricMean()
        return data_get_stat ( data , statobj , expression , '' , *args )  
    
    statobj = Ostap.Math.WGeometricMean()
    return data_get_stat ( data , statobj , expression , cuts , *args )  

# =============================================================================
## Get power mean over the data
#  @code
#  data = ...
#  result = data_power_mean ( data , 5 , 'pt' , 'eta>0' )
#  @endcode 
#  @see Ostap::Math::Moment
#  @see Ostap::Math::WMoment
#  @see Ostap::Math::PowerMean
#  @see Ostap::Math::WPowerMean
#  @see Ostap::statVar::the_moment
def data_power_mean ( data , p , expression , cuts = '' , *args ) :
    """ Get power mean over the data
    >>> data = ...
    >>> result = data_power_mean ( data , 5 , 'pt' , 'eta>0' )
    - see `Ostap.Math::Moment`
    - see `Ostap.Math::WMoment`
    - see `Ostap.Math::PowerMean`
    - see `Ostap.Math::WPowerMean`
    - see `Ostap.statVar.the_moment`
    """
    assert isinstance ( p , num_types ) , 'Invalid p-parameter type: %s' % type ( p ) 

    if    p == -1 or isequal ( p , -1. ) : return data_harmonic_mean   ( data , expression , cuts , *args ) 
    elif  p ==  0 or iszero  ( p       ) : return data_geometric_mean  ( data , expression , cuts , *args ) 
    elif  p ==  1 or isequal ( p ,  1. ) : return data_arithmetic_mean ( data , expression , cuts , *args ) 
    
    if isinstance ( data , ROOT.TTree ) and not cuts : 
        statobj = Ostap.Math.PowerMean      ( p )
        return data_get_stat ( data , statobj , expression , '' , *args )  
    
    statobj = Ostap.Math.WPowerMean         ( p )
    return data_get_stat ( data , statobj , expression , cuts , *args )

# =============================================================================
## Get arithmetic mean over the data (just for completeness)
#  @code
#  data = ...
#  result = data_arithmetic_mean ( data , 'pt' , 'eta>0' )
#  @endcode 
#  @see Ostap::Math::Moment
#  @see Ostap::Math::WMoment
#  @see Ostap::Math::ArithmeticMean
#  @see Ostap::Math::WArithmeticMean
#  @see Ostap::statVar::the_moment
def data_arithmetic_mean ( data , expression , cuts = '' , *args ) :
    """ Get power mean over the data (just for completeness)
    >>> data = ...
    >>> result = data_arithmetic_mean ( data , 5 , 'pt' , 'eta>0' )
    - see `Ostap.Math::Moment`
    - see `Ostap.Math::WMoment`
    - see `Ostap.Math::ArithmeticMean`
    - see `Ostap.Math::WArithmeticMean`
    - see `Ostap.statVar.the_moment`
    """
    
    if isinstance ( data , ROOT.TTree ) and not cuts : 
        statobj = Ostap.Math.ArithmeticMean (   )
        return data_get_stat ( data , statobj , expression , '' , *args )  
    
    statobj = Ostap.Math.WArithmeticMean (   )
    return data_get_stat ( data , statobj , expression , cuts , *args )

# =============================================================================
## Get Lehmer mean over the data
#  @code
#  data = ...
#  result = data_lehmer_mean ( data , 5 , 'pt' , 'eta>0' )
#  @endcode 
#  @see Ostap::Math::Moment
#  @see Ostap::Math::WMoment
#  @see Ostap::Math::LehmerMean
#  @see Ostap::Math::WLehmerMean
#  @see Ostap::statVar::the_moment
def data_lehmer_mean ( data , p , expression , cuts = '' , *args ) :
    """ Get Lehmer mean over the data
    >>> data = ...
    >>> result = data_lehmer_mean ( data , 5 , 'pt' , 'eta>0' )
    - see `Ostap.Math::Moment`
    - see `Ostap.Math::WMoment`
    - see `Ostap.Math::LehmerMean`
    - see `Ostap.Math::WLehmerMean`
    - see `Ostap.statVar.the_moment`
    """
    
    assert isinstance ( p , num_types ) , 'Invalid p-parameter!'
    
    if   p == 0 or iszero  ( p       ) : return data_harmonic_mean   ( data , expression , cuts , *args ) 
    elif p == 1 or isequal ( p , 1.0 ) : return data_arithmetic_mean ( data , expression , cuts , *args ) 

    if ( not cuts ) and isinstance ( data , ROOT.TTree ) :
        statobj = Ostap.Math.LehmerMean  ( p ) 
        return data_get_stat ( data , statobj , expression , '' , *args )  

    statobj = Ostap.Math.WLehmerMean     ( p )     
    return data_get_stat ( data , statobj , expression , cuts , *args )  

# =============================================================================
## Get the moments or order <code>order</code> as <code>Ostap::Math::(W)Moment_<order></code>
#  @code
#  data = ...
#  moment = data.the_momemnt ( 5 , 'x/y+z' , '0<qq' )
#  @endcode 
#  @see Ostap::Math::Moment 
#  @see Ostap::Math::WMoment
def data_the_moment ( data , order , expression , cuts = '' , *args ) :
    """Get the moments or order <code>order</code> as <code>Ostap::Math::(W)Moment_<order></code>
    >>> data = ...
    >>> moment = data.the_moment ( 5 , 'x/y+z' , '0<qq' )
    - see Ostap.Math.Moment 
    - see Ostap.Math.WMoment 
    """
    assert isinstance ( order  , int ) and 0 <= order , 'Invalid order  %s'  % order
    
    assert isinstance ( expression , expression_types ) , 'Invalid type of expression!'
    assert isinstance ( cuts       , expression_types ) , 'Invalid type of cuts/weight!'
    
    expression = str ( expression )
    cuts       = str ( cuts       ).strip() 
    
    if isinstance ( data , ROOT.TTree ) and not cuts :
        M      = Ostap.Math. Moment_(order)
        moment = M ()
        with rootException() : 
            sc = StatVar.the_moment ( data , moment , expression , *args )
            assert sc.isSuccess() , 'Error %s from StatVar::the_moment' % sc 
            return moment 

    M      = Ostap.Math.WMoment_(order)
    moment = M ()
    with rootException() :
        sc = StatVar.the_moment ( data , moment , expression , cuts , *args )
        assert sc.isSuccess() , 'Error %s from StatVar::the_moment' % sc 
        return moment 

# =============================================================================
## Get the mean as the moment
#  @code
#  data  = ...
#  value = data_the_mean ( data , 'x' ,'0<y' ).mean()  
#  @endcode
def data_the_mean ( data , expressions , cuts , errors = True , *args ) :
    """Get the mean as the moment
    >>> data  = ...
    >>> value = data_the_mean ( data , 'x' ,'0<y' ).mean()  
    """
    return data_the_moment ( data , 2 if errors else 1 , expressions , cuts = cuts , *args )
# =============================================================================
## Get the RMS as the moment
#  @code
#  data  = ...
#  value = data_the_rms ( data , 'x' ,'0<y' ).rms()  
#  @endcode
def data_the_rms ( data , expressions , cuts , errors = True , *args ) :
    """Get the mean as the moment
    >>> data  = ...
    >>> value = data_the_rmsn ( data , 'x' ,'0<y' ).rms()  
    """
    return data_the_moment ( data , 4 if errors else 2 , expressions , cuts = cuts , *args )
# =============================================================================
## Get the variance as the moment
#  @code
#  data  = ...
#  value = data_the_variance ( data , 'x' ,'0<y' ).variance()  
#  @endcode
def data_the_variance  ( data , expressions , cuts , errors = True , *args ) :
    """Get the mean as the moment
    >>> data  = ...
    >>> value = data_the_variance  ( data , 'x' ,'0<y' ).variance ()  
    """
    return data_the_moment ( data , 4 if errors else 2 , expressions , cuts = cuts , *args )
# =============================================================================
## Get the skewness as the moment
#  @code
#  data  = ...
#  value = data_the_skewness  ( data , 'x' ,'0<y' ).skewness ()  
#  @endcode
def data_the_skewness ( data , expressions , cuts , errors = True , *args ) :
    """Get the skewness as the moment
    >>> data  = ...
    >>> value = data_the_skewness  ( data , 'x' ,'0<y' ).skewness ()  
    """
    return data_the_moment ( data , 6 if errors else 3 , expressions , cuts = cuts , *args )
# =============================================================================
## Get the kurtosis as the moment
#  @code
#  data  = ...
#  value = data_the_kurtosis  ( data , 'x' ,'0<y' ).kurtosis ()  
#  @endcode
def data_the_kurtosis ( data , expressions , cuts , errors = True , *args ) :
    """Get the kurtosis as the moment
    >>> data  = ...
    >>> value = data_the_kurtosis  ( data , 'x' ,'0<y' ).kurtosis ()  
    """
    return data_the_moment ( data , 8 if errors else 4 , expressions , cuts = cuts , *args )
# =============================================================================

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
    
    assert isinstance ( expression , expression_types ) , 'Invalid type of expression!'
    assert isinstance ( cuts       , expression_types ) , 'Invalid type of cuts/weight!'
    
    expression = str ( expression )
    cuts       = str ( cuts       ).strip()
    
    with rootException() : 
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
    with rootException() : 
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

    assert isinstance ( expression , expression_types ) , 'Invalid type of expression!'
    assert isinstance ( cuts       , expression_types ) , 'Invalid type of cuts/weight!'
    
    expression = str ( expression )
    cuts       = str ( cuts       ).strip()
    
    if   exact is True  :
        ##  exact slow algorithm 
        with rootException() : 
            qn = StatVar.  quantile ( data , q , expression , cuts , *args )

    elif exact is False :
        ## approximate fast formula 
        with rootException() : 
            qn = StatVar.p2quantile ( data , q , expression , cuts , *args )
            
    elif isinstance ( exact , int ) and len ( data ) <= exact  :
        ## exact slow algorithm 
        with rootException() : 
            qn = StatVar.  quantile ( data , q , expression , cuts , *args )

    else :
        ## approximate fast algorithm 
        with rootException() : 
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

    assert isinstance ( expression , expression_types ) , 'Invalid type of expression!'
    assert isinstance ( cuts       , expression_types ) , 'Invalid type of cuts/weight!'
    
    expression = str ( expression )
    cuts       = str ( cuts       ).strip()

    qmin, qmax = min ( qmin ,  qmax ) , max  ( qmin, qmax  )
    
    if   exact is True  :
        ##  exact slow algorithm 
        with rootException() : 
            rn = StatVar.  interval ( data , qmin , qmax  , expression , cuts , *args )

    elif exact is False :
        ## approximate fast formula 
        with rootException() : 
            rn = StatVar.p2interval ( data , qmin , qmax , expression , cuts , *args )

    elif isinstance ( exact , int ) and len ( data ) <= exact  :
        ## exact slow algorithm 
        with rootException() : 
            rn = StatVar.  interval ( data , qmin , qmax  , expression , cuts , *args )

    else :
        ## approximate fast formula 
        with rootException() : 
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
    with rootException() : 
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
    assert isinstance ( expression , expression_types ) , 'Invalid type of expression!'
    assert isinstance ( cuts       , expression_types ) , 'Invalid type of cuts/weight!'
    
    expression = str ( expression )
    cuts       = str ( cuts       ).strip()

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
        with rootException() : 
            qn = StatVar.  quantiles ( data , qqq , expression , cuts , *args )

    elif exact is False :
        ## approximate fast formula 
        with rootException() : 
            qn = StatVar.p2quantiles ( data , qqq , expression , cuts , *args )
            
    elif isinstance ( exact , int ) and len ( data ) <= exact  :
        ## exact slow algorithm 
        with rootException() : 
            qn = StatVar.  quantiles ( data , qqq , expression , cuts , *args )

    else :
        ## approximate fast algorithm 
        with rootException() : 
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

# =============================================================================
## decorate certain class with some useful  methods 
def data_decorate ( klass ) :
    """Decorate certain class with some useful  methods
    """
    
    if hasattr ( klass , 'get_moment'     ) : klass.orig_get_moment     = klass.get_moment
    if hasattr ( klass , 'moment'         ) : klass.orig_moment         = klass.moment
    if hasattr ( klass , 'central_moment' ) : klass.orig_central_moment = klass.central_moment
    if hasattr ( klass , 'nEff'           ) : klass.orig_nEff           = klass.nEff
    if hasattr ( klass , 'mean'           ) : klass.orig_mean           = klass.mean
    if hasattr ( klass , 'variance'       ) : klass.orig_variance       = klass.variance 
    if hasattr ( klass , 'dispersion'     ) : klass.orig_dispersion     = klass.dispersion
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
    klass.nEff            = data_nEff 
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

    if hasattr ( klass , 'get_stats'       ) : klass.orig_get_stats       = klass.get_stats
    
    if hasattr ( klass , 'statVar'         ) : klass.orig_statVar         = klass.statVar 
    if hasattr ( klass , 'statVars'        ) : klass.orig_statVars        = klass.statVars
    if hasattr ( klass , 'sumVar'          ) : klass.orig_sumVar          = klass.sumVar 
    if hasattr ( klass , 'sumVars'         ) : klass.orig_sumVars         = klass.sumVars
    if hasattr ( klass , 'statCov'         ) : klass.orig_statCov         = klass.statCov
    
    klass.statVar  = data_statistics
    klass.statVars = data_statistics    
    klass.sumVar   = data_sum 
    klass.sumVars  = data_sum 
    klass.statCov  = data_covariance
    
    if hasattr ( klass , 'the_moment'      ) : klass.orig_the_moment      = klass.the_moment
    if hasattr ( klass , 'the_mean'        ) : klass.orig_the_mean        = klass.the_mean
    if hasattr ( klass , 'the_rms'         ) : klass.orig_the_rms         = klass.the_rms 
    if hasattr ( klass , 'the_variance'    ) : klass.orig_the_variance    = klass.the_variance
    if hasattr ( klass , 'the_skewness'    ) : klass.orig_the_skewness    = klass.the_skewness
    if hasattr ( klass , 'the_kurtosis'    ) : klass.orig_the_kurtosis    = klass.the_kurtosis
    
    if hasattr ( klass , 'harmonic_mean'   ) : klass.orig_harmonic_mean   = klass.harmonic_mean
    if hasattr ( klass , 'geometric_mean'  ) : klass.orig_geometric_mean  = klass.geometric_mean
    if hasattr ( klass , 'power_mean'      ) : klass.orig_power_mean      = klass.power_mean
    if hasattr ( klass , 'lehmer_mean'     ) : klass.orig_lehmer_mean     = klass.lehmer_mean
    if hasattr ( klass , 'arithmetic_mean' ) : klass.orig_arithmetic_mean = klass.arithmetic_mean

    klass.get_stats       = data_get_stat
    klass.the_moment      = data_the_moment
    klass.the_mean        = data_the_mean 
    klass.the_rms         = data_the_rms
    klass.the_variance    = data_the_variance
    klass.the_skewness    = data_the_skewness 
    klass.the_kurtosis    = data_the_kurtosis    
    
    klass.harmonic_mean   = data_harmonic_mean 
    klass.geometric_mean  = data_geometric_mean 
    klass.power_mean      = data_power_mean 
    klass.lehmer_mean     = data_lehmer_mean 
    klass.arithmetic_mean = data_arithmetic_mean 

    return ( klass.get_moment      , 
             klass.moment          , 
             klass.central_moment  , 
             klass.mean            , 
             klass.nEff            , 
             klass.variance        , 
             klass.dispersion      , 
             klass.rms             , 
             klass.skewness        , 
             klass.kurtosis        ,              
             klass.quantile        , 
             klass.quantiles       , 
             klass.interval        , 
             klass.median          , 
             klass.terciles        , 
             klass.quartiles       , 
             klass.quintiles       , 
             klass.deciles         , 
             klass.get_stats       ,
             klass.statVar         ,
             klass.statVars        ,
             klass.sumVar          ,
             klass.sumVars         ,
             klass.statCov         ,
             klass.the_moment      ,
             klass.the_mean        ,
             klass.the_rms         ,
             klass.the_variance    ,
             klass.the_skewness    ,
             klass.the_kurtosis    ,             
             klass.harmonic_mean   , 
             klass.geometric_mean  , 
             klass.power_mean      ,
             klass.lehmer_mean     ,
             klass.arithmetic_mean ) 

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


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
