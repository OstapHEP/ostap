#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file statvars.py 
#  Functions to collect statistics for trees and  datasets
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06  
# =============================================================================
""" Functions to collect statistics for trees and  datasets
- data_get_stat        - generic statistics 
- data_the_moment      - get colleciton of moment  
- data_moment          - get the moment  (with uncertainty) 
- data_ECDF            - get the emptipical cumulative distribution function 
- data_statistics      - get statistics ast StatEntity/WStatEntity objects 
- data_minmax          = get min/max 
- data_range           - get suitable rangess for drawing 
- data_covariance      - get the covariaces 
- data_statvector      - get the covariaces 
- data_sun             - get the sum 
- data_nEff            - get umber of effective entries 
- data_mean            - get the mean              (with uncertainty)
- data_variance        - get the variance          (with uncertainty)
- data_dispersion      - get the dispersion        (with uncertainty)
- data_rms             - get the RMS               (with uncertainty)
- data_skewness        - get the skewness          (with uncertainty)
- data_kurtosis        - get the (excess) kurtosis (with uncertainty)
- data_harmonic_mean   - get the (weighted) harmonic  mean 
- data_geometric_mean  - get the (weighte)  geometric mean 
- data_arithmetic_mean - get the (weighted) geometric mean
- data_power_mean      - get the (weighted) power     mean 
- data_lehmer_mean     - get the (weighted) Lehmer    mean 
- data_quantiles       - get the quantiles 
- data_median          - get the median
- data_terciles        - get the terciles 
- data_quartiles       - get the quartiles  
- data_quintiles       - get the quintiles 
- data_sextiles        - get the sextiles 
- data_septiles        - get the septiles 
- data_octiles         - get the octiles  
- data_deciles         - get the deciles  
- data_ventiles        - get the ventiles  
- data_percentiles     - get the percentiles  
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-06"
__all__     = (
    'data_get_stat'        , ## generic statistics 
    'data_the_moment'      , ## get colleciton of moment  
    'data_ECDF'            , ## get the emptipical cumulative distribution function 
    'data_statistic'       , ## get statistic ast StatEntity/WStatEntity objects 
    'data_minmax'          , ## get min/max 
    'data_range'           , ## get suitable rangess for drawing 
    'data_covariance'      , ## get the covariaces
    'data_statvector'      , ## get the covariaces
    'data_moment'          , ## get the moment  (with uncertainty)     
    'data_sum'             , ## get the sum 
    'data_nEff'            , ## get umber of effective entries 
    'data_mean'            , ## get the mean              (with uncertainty)
    'data_variance'        , ## get the variance          (with uncertainty)
    'data_dispersion'      , ## get the dispersion        (with uncertainty)
    'data_rms'             , ## get the RMS               (with uncertainty)
    'data_skewness'        , ## get the skewness          (with uncertainty)
    'data_kurtosis'        , ## get the kurtosis          (with uncertainty)
    'data_harmonic_mean'   , ## get the (weighted) harmonic   mean 
    'data_geometric_mean'  , ## get the (weighted) geometric  mean 
    'data_arithmetic_mean' , ## get the (weighted) arithmetic mean
    'data_power_mean'      , ## get the (weighted) power      mean 
    'data_lehmer_mean'     , ## get the (weighted) Lehmer     mean 
    'data_quantiles'       , ## get the quantiles 
    'data_median'          , ## get the median
    'data_terciles'        , ## get the terciles 
    'data_quartiles'       , ## get the quartiles  
    'data_quintiles'       , ## get the quintiles 
    'data_sextiles'        , ## get the sextiles 
    'data_septiles'        , ## get the septiles 
    'data_octiles'         , ## get the octiles  
    'data_deciles'         , ## get the deciles  
    'data_ventiles'        , ## get the ventiles  
    'data_percentiles'     , ## get the percentiles  
    ##
    'data_decorate'        , ## technical function to decorate the classes
    'expression_types'     , ## valid types for expressions/cuts/weights
)
# =============================================================================
from   ostap.math.base                import ( isequal     , iszero    ,
                                               axis_range  ,
                                               strings     , doubles   ,  
                                               all_entries , evt_range )      
from   ostap.core.core                 import Ostap, rootException, WSE, VE, std     
from   ostap.core.ostap_types          import ( string_types   , integer_types  , 
                                                num_types      , dictlike_types ,
                                                sequence_types )
from   ostap.trees.cuts                import expression_types, vars_and_cuts
from   ostap.utils.basic               import loop_items, typename, numcpu
from   ostap.utils.progress_conf       import progress_conf
import ostap.frames.frames             as     F 
import ostap.parallel.parallel_statvar as P 
import ostap.logger.table              as     T
import ostap.stats.counters 
import ostap.stats.moment 
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.statvars' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
LARGE   = 50000    ## allow frames or parallel for LARGE datasets 
# =============================================================================
MIN_ENTRIES_FOR_FRAME    = 100000
MIN_ENTRIES_FOR_PARALLEL = 200000
MIN_FILES_FOR_PARALLEL   = 2
# =============================================================================
## Good for provcssing via frames?
def good_for_frame    ( data , *args ,
                        min_events = MIN_ENTRIES_FOR_FRAME ) :
    """ Good for provcessing via frames?
    """
    if 2 < len ( args )                                    : return False
    
    ## Not TTree ?
    if not isinstance ( data , ROOT.TTree  )               : return False
    
    ## do we have at least two CPUs?
    if numcpu() < 2                                        : return False

    ## Is multithreading enabled? 
    if not ROOT.ROOT.IsImplicitMTEnabled()                 : return False 

    ## check number of events 
    first, last = evt_range ( data , *args[:2] )
    ## dataset is too small
    return 0 <= first < last and min_events <= ( last - first ) and all_entries ( data , first , last ) 

# ============================================================================
## good for parallel processing 
def good_for_parallel ( data ,
                        *args ,
                        min_events = MIN_ENTRIES_FOR_PARALLEL , 
                        min_files  = MIN_FILES_FOR_PARALLEL   ) :
    
    ## Not TTree ?
    if not isinstance ( data , ROOT.TTree  )               : return False

    ## do we have at least two CPUs?
    if numcpu() < 2                                        : return False

            
    first, last = evt_range ( data , *args[:2] )
    if 0 <= first < last and min_events < ( last - first ) : return True ## ATTENTION

    if isinstance ( data , ROOT.TChain ) : return min_files <= data.nFiles
    
    return False

# =============================================================================
StatVar = Ostap.StatVar
# =============================================================================
## reset the target 
def target_reset ( obj ) :
    """ Reset the target """
    if instance ( obj , ROOT.TH1 ) :
        obj.Reset()
        if not obj.GetSumw2() : obj.Sumw2()
    elif hasattr ( obj , 'Reset' )  : obj.Reset ()
    elif hasattr ( obj , 'reset' )  : obj.reset ()
# =============================================================================    
## Copy target 
def target_copy ( self , obj ) :
    """ Copy target """
    if instance ( obj , ROOT.TH1 ) :
        newoj = obj.Clone ()
        newobj.Reset()
        if not newobj.GetSumw2() : newobj.Sumw2()
        return newobj
    elif hasattr ( obj , 'Clone' ) : return obj.Clone()
    elif hasattr ( obj , 'clone' ) : return obj.clone()
    ## 
    T = type ( obj )
    return T ( obj )
# =============================================================================
_s1D = Ostap.Math.Statistic  , Ostap.Math.WStatistic
_s2D = Ostap.Math.Statistic2 , Ostap.Math.WStatistic2
_s3D = Ostap.Math.Statistic3 , Ostap.Math.WStatistic3
_s4D = Ostap.Math.Statistic4 , Ostap.Math.WStatistic4
# =============================================================================
## Get the (s)Statistic-bases statistic/counter from data
#  @code
#  statobj = Ostap.Math.MinValue()
#  data    = ...
#  result  = data_get_stat ( data , statobj , 'x+y' , cuts = 'pt>1' ) 
#  @encode
#  @see Ostap::Math::Statistic
#  @see Ostap::Math::WStatistic
#  @see Ostap::Math::Statistic2
#  @see Ostap::Math::WStatistic2
#  @see Ostap::Math::Statistic3
#  @see Ostap::Math::WStatistic3
#  @see Ostap::Math::Statistic4
#  @see Ostap::Math::WStatistic4
def data_get_stat ( data               ,
                    statobj            ,
                    expressions        ,
                    cuts       = ""    , *args ,  
                    cut_range  = ""    ,
                    progress   = False , 
                    use_frame  = False ,
                    parallel   = False ) :
    """ Get the (W)Statistic-based statistic/counters from data  
    >>> data   = ...
    >>> stat   = Ostap.Math.HarmonicMean() 
    >>> result = data.get_stat ( stat , 'x/y+z' , '0<qq' )
    - see Ostap.Math.Statistic
    - see Ostap.Math.WStatistic
    - see Ostap.Math.Statistic2
    - see Ostap.Math.WStatistic2
    - see Ostap.Math.Statistic3
    - see Ostap.Math.WStatistic3
    - see Ostap.Math.Statistic4
    - see Ostap.Math.WStatistic4
    """
    
    ## (1) decode expressions & cuts
    var_lst , cuts, _  = vars_and_cuts ( expressions , cuts )
    nvars = len ( var_lst )

    ## (2) check consistency
    if   1 == nvars and isinstance ( statobj , _s1D ) : pass
    elif 2 == nvars and isinstance ( statobj , _s2D ) : pass
    elif 3 == nvars and isinstance ( statobj , _s3D ) : pass
    elif 4 == nvars and isinstance ( statobj , _s4D ) : pass
    else :
        raise TypeError ( 'Inconsistent statobj %s & expression(s): %s' % \
                          ( typename ( statobj ) , str ( var_lst ) ) ) 

    ## (3) cut_range defined *only* for RooFit datasets 
    if cut_range and not isinstance ( data , ROOT.RooAbsData ) : 
        raise TypeError ( "Invalid use of `cut_range':%s" % cut_range  ) 

    ## (4) display progress ? 
    progress = progress_conf ( progress )

    ## (5) create the driver 
    sv = StatVar ( progress )

    ## (6) RooFit ?
    if isinstance ( data , ROOT.RooAbsData ) :
        
        weighted          = data.isWeighted ()
        store_errors      = weighted and data.store_errors      ()
        store_asym_errors = weighted and data.store_asym_errors ()         
        if store_errors or store_asym_errors :
            logger.warning ( "Weight uncertainties are defined, but will be ignored!" ) 

        with rootException() :
            the_args = var_lst + ( cuts , cut_range ) + args 
            sc       = sv.get_stat ( data , statobj , *the_args )
        assert sc.isSuccess() , 'Error %s from StatVar::get_stat' % sc 
        return statobj

    ## Use frame processing ?
    if   use_frame and good_for_frame ( data , *args ) : 
        return F.frame_project ( data                   ,
                                 model       = statobj  ,
                                 expressions = var_lst  ,
                                 cuts        = cuts     ,
                                 progress    = progress ,
                                 report      = progress ,
                                 lazy        = False    )
    ## Use parallel processon 
    elif parallel and good_for_parallel ( data , *args ) : 
        from ostap.parallel.parallel_stavar import parallel_get_stat        
        return parallel_get_stat ( data        ,
                                   statobj     ,
                                   expressions ,
                                   cuts        , *args         ,  
                                   progress    = progress      ,
                                   use_frame   = False         , ## NB!!
                                   chunk_size  = 2 * LARGE     ,
                                   max_files   = 1             ,
                                   silent      = not progress  ) ;
    
    assert isinstance ( data , ROOT.TTree ) , "Here data must be TTree: %s" % typename ( data ) 
    
    ## Branches to be activated
    from ostap.trees.trees import ActiveBranches
    with rootException() , ActiveBranches ( data , cuts , *var_lst ) :
        the_args = var_lst + ( cuts , ) + args         
        sc       = sv.get_stat  ( data , statobj , *the_args  )
        assert sc.isSuccess() , 'Error %s from StatVar::the_moment' % sc 
        return statobj
    
# =============================================================================
## get the moment of order 'order' relative to 'center'
#  @code
#  data =  ...
#  print data_the_moment ( data , 3 , 0.0 , 'mass' , 'pt>1' ) 
#  @endcode
#  @see Ostap::Math::Moment_
#  @see Ostap::Math::WMoment_
def data_the_moment ( data               ,
                      order              ,
                      expression         ,
                      cuts       = ''    , *args , 
                      cut_range  = ''    ,
                      as_weight  = True  , ## interpret cuts as weight 
                      progress   = False , 
                      use_frame  = False ,
                      parallel   = False ) : 
    """ Get the moment of order 'order'
    >>> data =  ...
    >>> print data_the_moment ( data , 3 , 0.0 , 'mass' , cuts = 'pt>1' )
    >>> print data.the_moment (        3 , 0.0 , 'mass' , cuts = 'pt>1' ) ## ditto
    #  @see Ostap::Math::WMoment_
    """   
    assert isinstance ( order , integer_types ) and 0 <= order , \
        'Invalid moment order: %s'  % order
    
    ## (1) decode expressions & cuts
    _ , cuts , _  = vars_and_cuts ( expression , cuts )

    if isinstance ( data , ROOT.RooAbsData ) or ( cuts and as_weight ) : 
        statobj = Ostap.Math.WMoment_[order] ()
    else :
        statobj = Ostap.Math. Moment_[order] ()

    return data_get_stat ( data                   ,
                           statobj                ,
                           expression             , 
                           cuts       , *args     ,
                           cut_range  = cut_range ,
                           progress   = progress  , 
                           use_frame  = use_frame , 
                           parallel   = parallel  )

# =============================================================================
## Is there at least one good entry ?
def data_hasEntry ( data               ,
                    cuts      = ""     , *args , 
                    cut_range = ""     ,
                    progress  = False  ,
                    use_frame = False  ,
                    parallel  = False  ) :
    """ Is there at least one good entry? 
    """
    ## (1) decode expressions & cuts
    _ , cuts , _  = vars_and_cuts ( "1" , cuts )

    ##  (2) trivial case 
    first , last = evt_range ( data , *args[:2] )
    if last <= first : return False
    
    ## (2) cut_range defined *only* for RooFit datasets 
    if cut_range and not isinstance ( data , ROOT.RooAbsData ) : 
        raise TypeError ( "Invalid use of `cut_range':%s" % cut_range  )

    ## (4) display progress ? 
    progress = progress_conf ( progress )
    
    ## (5) create the driver 
    sv = StatVar ( progress )

    ## (6) RooFit ?
    if isinstance ( data , ROOT.RooAbsData ) :
        ## trivial case 
        if not cuts and not cut_range and not data.isWeigted() : return True
        
        weighted          = data.isWeighted ()
        store_errors      = weighted and data.store_errors      ()
        store_asym_errors = weighted and data.store_asym_errors ()         
        if store_errors or store_asym_errors :
            logger.warning ( "Weight uncertainties are defined, but will be ignored!" ) 

        with rootException() :
            return sv.hasEntry ( data , cuts , cut_range , *args )
        
    ## (7) trivial case 
    if not cuts : return True 

    ## use frame ? 
    if use_frame and good_for_frame ( data , *args ) :
        return 0 < F.frame_size ( data                ,
                                  cuts     = cuts     , 
                                  progress = progress ,
                                  report   = progress ,
                                  lazy     = False    )     
    ## parallel ? 
    elif parallel and good_for_parallel ( data , *args ) :
        from ostap.parallel.parallel_statvar import parallel_size 
        return 0 < paralel_size ( data                  ,
                                  cuts      , *args     , 
                                  progress  = progress  , 
                                  use_frame = use_frame )

    assert isinstance ( data , ROOT.TTree ) , "Here data must be TTree: %s" % typename ( data ) 

    ## Branches to be activated
    from ostap.trees.trees import ActiveBranches
    with rootException() , ActiveBranches ( cuts ) :
        return sv.hasEntry ( data , cuts , *args )

# =============================================================================
## How many good entries in dataset?
def data_size ( data               ,       
                cuts      = ""     , *args , 
                cut_range = ""     ,
                progress  = False  ,
                use_frame = False  ,
                parallel =  False  ) :
    """ How many  good entries in the dataset ?
    """
    ## (1) decode expressions & cuts
    _ , cuts , _  = vars_and_cuts ( "1" , cuts )

    ##  (2) trivial case 
    first , last = evt_range ( data , *args[:2] )
    if last <= first : return 0 

    ## (2) cut_range defined *only* for RooFit datasets 
    if cut_range and not isinstance ( data , ROOT.RooAbsData ) : 
        raise TypeError ( "Invalid use of `cut_range':%s" % cut_range  )

    ## (4) display progress ? 
    progress = progress_conf ( progress )
    
    ## (5) create the driver 
    sv = StatVar ( progress )

    ## (6) RooFit ?
    if isinstance ( data , ROOT.RooAbsData ) :
        ## trivial case        
        if not cuts and not cut_range and not data.isWeighted() : return last - first 
        
        weighted          = data.isWeighted ()
        store_errors      = weighted and data.store_errors      ()
        store_asym_errors = weighted and data.store_asym_errors ()         
        if store_errors or store_asym_errors :
            logger.warning ( "Weight uncertainties are defined, but will be ignored!" ) 
            
        with rootException() :
            return sv.size ( data , cuts , cut_range , *args )

    ## (7) trivial case 
    if not cuts : return last - first 

    ## use frame ? 
    if   use_frame and good_for_frame ( data , *args ) : 
        return F.frame_size ( data                 ,
                              cuts     = cuts      , 
                              progress = progress  ,
                              report   = progress  ,
                              lazy     = False     )
    ## parallel ?
    elif parallel and good_for_parallel ( data , *args ) : 
        from ostap.parallel.parallel_statvar import parallel_size 
        return paralel_size ( data                  ,
                              cuts      , *args     ,
                              progress  = progress  , 
                              use_frame = use_frame )
    
    assert isinstance ( data , ROOT.TTree ) , "Here data must be TTree: %s" % typename ( data ) 

    ## Branches to be activated
    from ostap.trees.trees import ActiveBranches
    with rootException() , ActiveBranches ( cuts ) :
        return sv.size ( data , cuts , *args )
    
# =============================================================================

# =============================================================================
## get the moment (with uncertainty) of order 'order' 
#  @code
#  data =  ...
#  print data_moment ( data , 3 , 'mass' , 'pt>1' ) 
#  print data.moment (        3 , 'mass' , 'pt>1' ) ## ditto
#  @endcode
#  @see Ostap::StatVar::moment
def  data_moment ( data               ,
                   order              ,
                   expression         ,
                   cuts       = ''    , *args , 
                   cut_range  = ''    ,
                   as_weight  = True  , ## interpret cuts s as weiggt 
                   progress   = False , 
                   use_frame  = False ,
                   parallel   = False ) : 
    """ Get the moment of order 'order' relative to 'center'
    >>> data =  ...
    >>> print data_moment ( data ,  3 , 'mass' , 'pt>1' ) 
    >>> print data.moment (         3 , 'mass' , 'pt>1' ) ## ditto
    - see Ostap::StatVar::moment
    """
    assert isinstance ( order , integer_types ) and 0 <= order , 'Invalid order %s'  % order

    the_moment = data_the_moment ( data       ,
                                   2 * order  ,
                                   expression ,
                                   cuts       , *args     , 
                                   cut_range  = cut_range ,
                                   progress   = progress  , 
                                   as_weight  = as_weight , 
                                   use_frame  = use_frame , 
                                   parallel   = parallel  )
    ## get the actual moment 
    return the_moment.moment ( order )

# =============================================================================
## Get the empirical cumulative distributtion function 
#  @code
#  data =  ...
#  ecdf1 = data_ECDF ( data , 'mass' ) 
#  ecdf2 = data_ECDF ( data , 'mass' ,  cuts = 'PT>1') 
#  @endcode
#  @see Ostap::Math::ECDF
#  @see Ostap::Math::WECDF
def data_ECDF ( data               ,
                expression         ,
                cuts       = ''    , *args , 
                cut_range  = ''    ,
                progress   = False , 
                as_weight  = True  , ## interpret cuts as weiggt 
                use_frame  = False ,
                parallel   = False ) :     
    """ Get the empirical cumulative distribution function 
    >>> data =  ...
    >>> ecdf1 = data.ECDF ( data , 'mass' ) 
    >>> ecdf2 = data.ECDF ( data , 'mass' , cuts = 'PT>1') 
    - see `Ostap.Math.ECDF`
    - see `Ostap.Math.WECDF`
    """
    
    ## (1) decode expressions & cuts
    _ , cuts , _  = vars_and_cuts ( expression , cuts )
    
    if isinstance ( data , ROOT.RooAbsData ) or ( cuts and as_weight ) : 
        statobj = Ostap.Math.WECDF ()
    else :
        statobj = Ostap.Math. ECDF ()
        
    return data_get_stat ( data                  ,
                           statobj               ,
                           expression            , 
                           cuts                  , *args , 
                           cut_range = cut_range ,
                           progress  = progress  , 
                           use_frame = use_frame , 
                           parallel  = parallel  )

# ==============================================================================
## Get the statistic from data
#  @code
#  data    = ...
#  result  = data_statistics ( data , 'x+y'   , cuts = 'pt>1' ) 
#  results = data_statistics ( data , 'x;y;z' , cuts = 'pt>1' ) ## result is dictionary 
#  @encode
#  @see Ostap::StatEntity
#  @see Ostap::WStatEntity
#  @see Ostap::StatVar::statVar
def data_statistic ( data               , 
                     expressions        ,
                     cuts       = ''    , *args , 
                     cut_range  = ''    ,
                     progress   = False , 
                     as_weight  = True  , ## interpret cuts as weiggt 
                     use_frame  = False ,
                     parallel   = False ) :     
    """ Get statistics from data 
    >>> data    = ...
    >>> result  = data_statistics ( data , 'x/y+z' , cuts = '0<qq' )
    >>> results = data_statistics ( data , 'x/y;z' , cuts = '0<qq' ) ## result is dictionary
    - see Ostap.Math.StatEntity
    - see Ostap.Math.WStatEntity
    - see Ostap.StatVar.statVar
    """
    ## decode expressions & cuts
    var_lst, cuts, input_string = vars_and_cuts ( expressions , cuts )
    assert var_lst , "Invalid expressions!"


    ## Use frames? 
    if use_frame and good_for_frame ( data , *args ) : 
        return F.frame_statistic ( data        ,
                                   expressions ,
                                   cuts        = cuts      , 
                                   as_weight   = as_weight , 
                                   progress    = progress  ,
                                   report      = progress  ,
                                   lazy        = False     )
    ## Use parallel processing ?
    elif parallel and good_for_parallel ( data , *args ) : 
        from ostap.parallel.parallel_statvar import parallel_statistic
        return parallel_statistic ( data        ,
                                    expressions ,
                                    cuts        , *args      ,
                                    as_weight   = as_weight  , 
                                    progress    = progress   ,
                                    use_frame   = use_frame  )

    
    ##  diplay progress ? 
    progress = progress_conf ( progress )

    ## create the driver 
    sv = StatVar ( progress )
    
    if input_string :
        varname  = var_lst [ 0 ] 
        if   isinstance ( data , ROOT.RooAbsData ) :
            return sv.statVar     ( data , varname , cuts , cut_range , *args )
        elif cuts and as_weight                    :
            return sv.statVar     ( data , varname , cuts ,             *args )
        else : 
            return sv.statVar_cut ( data , varname , cuts ,             *args )
           
    ##     
    if   isinstance ( data , ROOT.RooAbsData ) :
        TCNT = Ostap.WStatEntity 
        vcnt = Ostap.StatVar.WStatVector ()
    elif cuts and as_weight                    : 
        TCNT = Ostap.WStatEntity            
        vcnt = Ostap.StatVar.WStatVector ()
    else                                       :
        TCNT = Ostap.StatEntity
        vcnt = Ostap.StatVar. StatVector ()

    ## variable names 
    vnames = strings ( var_lst ) 
    
    if isinstance ( data , ROOT.RooAbsData ) :
        
        weighted          = data.isWeighted ()
        store_errors      = weighted and data.store_errors      ()
        store_asym_errors = weighted and data.store_asym_errors ()         
        if store_errors or store_asym_errors :
            logger.warning ( "Weight uncertainties are defined, but will be ignored!" ) 
        
        with rootException() :
            sc   = sv.statVars ( data , vcnt , vnames  , cuts , cut_range , *args )
            assert sc.isSuccess() , 'Error %s from Ostap::StatVar::statVars' % sc  
            assert len ( vcnt ) == len ( vnames ) , "Mismatch in structure"
            return { name : TCNT ( cnt ) for ( name , cnt ) in zip ( var_lst , vcnt ) } 
            
    assert isinstance ( data , ROOT.TTree ) , "Invalid type for data!"
    
    from ostap.trees.trees import ActiveBranches
    with rootException() , ActiveBranches ( data  , cuts , *var_lst ) :
        sc   = sv.statVars ( data , vcnt , vnames , cuts , *args   ) 
        assert sc.isSuccess() , 'Error %s from Ostap::StatVar::statVars' % sc 
        assert len ( vcnt ) == len ( vnames ) , "Mismatch in structure"    
        return { name : TCNT ( cnt ) for ( name , cnt ) in zip ( var_lst , vcnt ) } 

# ==============================================================================
## Get the minmax from data
#  @code
#  data    = ...
#  result  = data_minmax ( data , 'x+y'   , 'pt>1' ) 
#  results = data_minmax ( data , 'x;y;z' , 'pt>1' ) ## result is dictionary 
#  @encode
#  @see Ostap::StatEntity
#  @see Ostap::WStatEntity
#  @see Ostap::StatVar::statVar
def data_minmax ( data , 
                  expressions        ,
                  cuts       = ''    , *args , 
                  cut_range  = ''    ,
                  progress   = False , 
                  use_frame  = False ,
                  parallel   = False ) :     
    """ Get min/max from data 
    >>> data    = ...
    >>> result  = data_minmax ( data , 'x/y+z' , '0<qq' )
    >>> results = data_minmax ( data , 'x/y;z' , '0<qq' ) ## result is dictionary
    - see Ostap.Math.StatEntity
    - see Ostap.Math.WStatEntity
    - see Ostap.StatVar.statVar
    """
    results = data_statistic ( data ,
                               expressions            , 
                               cuts                   , *args , 
                               cut_range = cut_range  ,
                               progress  = progress   ,
                               use_frame = use_frame  ,
                               parallel  = parallel   ) 
                                
    if isinstance ( results , dictlike_types ) :
        res = {} 
        for key, r in loop_items ( results ) :
            res  [ key ] = r.min() , r.max()
        results = res 
    else : results = results.min() , results.max()
    ## 
    return results 

# =============================================================================\
## Get suitable ranges for drawing expressions/variables
#  - In case there is no suitable range None is returned 
## @code
#  dataset = ...
#  result  = data_range ( dataset , 'sin(x)*100*y' , 'x<0' )
#  results = data_range ( dataset , 'x,y,z,t,u,v'  , 'x<0' ) ## as dictionary
#  @endcode
#  @see data_minmax
#  @see data_statistics
#  @see axis_range 
def data_range ( data               ,
                 expressions        ,
                 cuts       = ''    , *args , 
                 cut_range  = ''    ,
                 progress   = False , 
                 use_frame  = False ,
                 parallel   = False , 
                 delta      = 0.01  ) : 
    """ Get suitable ranges for drawing expressions/variables
    - In case there is no suitable range None is returned 
    >>> data = ...
    >>> result  = data_range ( data , 'sin(x)*100*y' , 'x<0' )
    >>> results = data_range ( dataset , 'x,y,z,t,u,v'  , 'x<0' ) ## as dictionary
    """
    results = data_minmax ( data                   ,
                            expressions            , 
                            cuts                   , *args , 
                            cut_range = cut_range  ,
                            progress  = progress   ,
                            use_frame = use_frame  ,
                            parallel  = parallel   ) 
    
    if isinstance ( results , dictlike_types ) :
        for k , r in loop_items ( results ) :
            mn, mx = r
            if mx < mn and ( cuts or cut_range or args ) : 
                ## recalculate without cuts and arg-ranges 
                mn , mx = data_minmax ( data , k ,        
                                        cuts      = ''         ,
                                        cut_range = ''         ,
                                        progress  = progress   ,
                                        use_frame = use_frame  ,
                                        parallel  = parallel   )  
            if mx < mn : return None               ## ATTENTION!                    
            results [ k ] = axis_range ( mn , mx , delta = delta ) 
    else :
        mn , mx = results
        if mx < mn and ( cuts or cut_range or args ) :
            mn , mx = data_minmax ( data                  ,
                                    expressions           , 
                                    progress  = progress  ,
                                    use_frame = use_frame , 
                                    parallel  = parallel  )
        if mx < mn : return None               ## ATTENTION! 
        return axis_range ( mn , mx , delta = delta ) 
      
    return results 

# ==============================================================================
## Get the covariance for expressions
#  @code
#  data    = ...
#  result  = data_covariance ( data , 'x+y' , 'z'  , '0<u') 
#  @encode
#  @see Ostap::Math::Covariance 
#  @see Ostap::Math::WCovariance 
#  @see Ostap::Math::Covariances 
#  @see Ostap::Math::WCovariances 
#  @see Ostap::StatVar::statCov
def data_covariance ( data        ,
                      expressions ,               
                      cuts        = ''    , *args , 
                      cut_range   = ''    ,
                      as_weight   = True  ,  
                      progress    = False , 
                      use_frame   = False ,
                      parallel    = False ) :
    """ Get the covariance from data 
    >>> data   = ...
    >>> result = data_covariance ( data , 'x/y+z' , 'qq' , 'u>0')
    - see Ostap.Math.Covariance
    - see Ostap.Math.WCovariance 
    - see Ostap.Math.Covariances
    - see Ostap.Math.WCovariances 
    - see Ostap.StatVar.statCov 
    """
    
    var_lst, cuts, _ =  vars_and_cuts ( expressions , cuts )
    assert 2 <= len ( var_lst ) , "At least two variables are required!"
    N = len ( var_lst ) 

    ## allow linear algebra manipulations are here 
    import ostap.math.linalg
    
    if isinstance ( data , ROOT.RooAbsData ) or ( cuts and as_weight ) :
        result = Ostap.Math.WCovariance() if 2 == N else Ostap.Math.WCovariances ( N )
    else :
        result = Ostap.Math. Covariance() if 2 == N else Ostap.Math. Covariances ( N )

    ##  diplay progress ? 
    progress = progress_conf ( progress )

    ## create the driver 
    sv = StatVar ( progress )
    
    vnames = strings ( var_lst )
    
    if isinstance ( data , ROOT.RooAbsData ) :
        
        weighted          = data.isWeighted ()
        store_errors      = weighted and data.store_errors      ()
        store_asym_errors = weighted and data.store_asym_errors ()         
        if store_errors or store_asym_errors :
            logger.warning ( "Weight uncertainties are defined, but will be ignored!" ) 

        with rootException() :
            if 2 == N : sc = sv.statCov ( data , result , var_lst[0] , var_lst[1] , cuts , cut_range , *args )
            else      : sc = sv.statCov ( data , result , vnames                  , cuts , cut_range , *args )
            assert sc.isSuccess() , 'Error %s from StatVar::statVars' % sc 
            return result 
        
    if  use_frame and good_for_frame ( data , *args )  and 2 == N : 
        return F.frame_covariance ( frame                 ,
                                    var_lst [ 0 ]         ,
                                    var_lst [ 1 ]         ,
                                    cuts      = cuts      ,
                                    as_weight = as_weight ,
                                    progress  = progress  ,
                                    lazy      = False     ) 
    
    if  parallel and good_for_parallel ( data , *args ) : 
        from ostap.parallel.parallel_statvar import parallel_covariance
        return parallel_covariance ( data        ,
                                     expressions ,
                                     cuts        , *args      , 
                                     as_weight   = as_weight  , 
                                     progress    = progress   ,
                                     use_frame   = use+_frame )
    
    assert isinstance ( data , ROOT.TTree ) , "Invalid type for data!"

    ## Branches to be activated
    from ostap.trees.trees import ActiveBranches
    with rootException() , ActiveBranches ( data , cuts , *var_lst ) :
        if 2 == N : sc = sv.statCov ( data , result , var_lst[0] , var_lst[1] , cuts , *args )
        else      : sc = sv.statCov ( data , result , vnames                  , cuts , *args )
        assert sc.isSuccess() , 'Error %s from StatVar::statVars' % sc 
        return result 

# ============================================================================
## get the effective vector of mean-values with covariances for the dataset
#  @code
#  vct = data_statvct ( data ,  'a,b,c' ) 
#  @endcode 
def data_statvector ( data        ,
                      expressions ,               
                      cuts        = ''    , *args , 
                      cut_range   = ''    ,
                      as_weight   = True  ,  
                      progress    = False , 
                      use_frame   = False ,
                      parallel    = False ) :
    """ Get the effective vector of mean-values with covariances for the dataset
    >>> vct = data_statvct ( data , 'a,y,z' , cuts = 'w' ) 
    """
    
    ## decode expressions & cuts
    var_lst , cuts, input_string = vars_and_cuts ( expressions , cuts )
    N = len ( var_lst )
    
    assert 2 <= N , "At least two variables are needed!"
        
    covs = data_covariance ( data       , var_lst   ,
                             cuts       , *args     , 
                             cut_range  = cut_range ,                             
                             as_weight  = as_weight ,  
                             progress   = progress  , 
                             use_frame  = use_frame ,
                             parallel   = parallel  )
    
    ## some linear algebra manipulations are here 
    import ostap.math.linalg

    VCT = Ostap.VectorE ( N )
    if 2 == N : 
        cov2 = Ostap.Math.covariance ( covs )
        vct  = VCT ( cov2 )
        vct.setValue ( 0 , covs.counter1().mean() ) 
        vct.setValue ( 1 , covs.counter2().mean() ) 
        return vct
    
    cov2 = covs.cov2    ()
    cov2 = cov2.smatrix () 
    vct  = VCT ( cov2 )    
    for i in range ( N ) : vct.setValue ( i , covs.counters()[i].mean() )
    return vct

# ==============================================================================
## Get the (weighted) sum over the variable(s)
#  @code
#  data    = ...
#  result  = data_sum ( data , 'x+y'   , 'pt>1' ) 
#  results = data_sum ( data , 'x;y;z' , 'pt>1' ) ## result is dictionary 
#  @encode
#  @see Ostap::StatVar::statVar
def data_sum ( data               ,
               expressions        ,
               cuts       = ''    , *args , 
               cut_range  = ''    ,
               progress   = False , 
               as_weight  = True  , ## interpret cuts as weiggt 
               use_frame  = False ,
               parallel   = False ) :     
    """ Get (weighted) sum over the variables 
    >>> data    = ...
    >>> result  = data_sum ( data , 'x/y+z' , '0<qq' )
    >>> results = data_sum ( data , 'x/y;z' , '0<qq' ) ## result is dictionary
    - see Ostap.StatVar.statVar
    """    
    result = data_statistic ( data                  ,
                              expressions           ,  
                              cuts                  , *args , 
                              cut_range = cut_range ,
                              progress  = progress  ,
                              as_weight = as_weight ,
                              use_frame = use_frame ,
                              parallel  = parallel  )
    result2 = {} 
    if isinstance ( result , dictlike_types ) :
        for key , r in loop_items ( result ) : 
            result2 [ key ] = r.sum()
    else :
        result2 = result.sum() 
    
    return result2 
        
# ==============================================================================
## Get the effective number of entries in data
#  \f$ \frac{ (\sum w)^2}{\sum w^2}  \f$ 
#  @code
#  data    = ...
#  result  = data_nEff ( data , 'x+y' ) 
#  @encode
#  @see Ostap::StatVar::nEff 
def data_nEff ( data ,
                cuts       = ''    , *args , 
                cut_range  = ''    ,
                as_weight  = True  , ## interpret cuts as weiggt 
                progress   = False , 
                use_frame  = False ,
                parallel   = False ) :         
    """ Get statistic from data 
    >>> data    = ...
    >>> result  = data_nEff ( data , 'x/y+z' )
    - see Ostap.StatVar.nEff
    """
    ## decode expressions & cuts
    var_lst, cuts, input_string = vars_and_cuts ( '1' , cuts )    
    assert 1 == len ( var_lst ) and input_string , "Invalid expressions: %s" % expression 
    expression = var_lst[0]

    stat = data_statistic ( data       ,
                            expression , 
                            cuts       , *args     , 
                            cut_range  = cut_range ,
                            as_weight  = as_weight ,
                            progress   = progress  ,
                            use_frame  = use_frame ,
                            parallel   = parallel  ) 
    return stat.nEff () 

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
def data_harmonic_mean ( data , 
                         expression         ,
                         cuts       = ''    , *args , 
                         cut_range  = ''    ,
                         progress   = False , 
                         as_weight  = True  , ## interpret cuts as weiggt 
                         use_frame  = False ,
                         parallel   = False ) :                
    """ Get harmonic mean over the data
    >>> data = ...
    >>> result = data_harmonic_mean ( data , 'pt' , 'eta>0' )
    - see `Ostap.Math::Moment`
    - see `Ostap.Math::WMoment`
    - see `Ostap.Math::HarmonicMean`
    - see `Ostap.Math::WHarmonicMean`
    - see `Ostap.statVar.the_moment`
    """    
    ## decode expressions & cuts
    var_lst, cuts, input_string = vars_and_cuts ( expression , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"
    
    ## 
    if isinstance ( data , ROOT.RooAbsData ) or ( cuts and as_weight ) :
        stat = Ostap.Math.WHarmonicMean ()
    else : 
        stat = Ostap.Math. HarmonicMean ()
        
    return data_get_stat  ( data       ,
                            stat       ,
                            expression ,
                            cuts       , *args     , 
                            cut_range  = cut_range ,
                            progress   = progress  ,
                            use_frame  = use_frame ,
                            parallel   = parallel  ) .value() 

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
def data_geometric_mean ( data , 
                          expression         ,
                          cuts       = ''    , *args , 
                          cut_range  = ''    ,
                          progress   = False , 
                          as_weight  = True  , ## interpret cuts as weiggt 
                          use_frame  = False ,
                          parallel   = False ) :                
    """ Get geometric mean over the data
    >>> data = ...
    >>> result = data_geometric_mean ( data , 'pt' , 'eta>0' )
    - see `Ostap.Math::Moment`
    - see `Ostap.Math::HarmonicMean`
    - see `Ostap.statVar.the_moment`
    """
    ## decode expressions & cuts
    var_lst, cuts, input_string = vars_and_cuts ( expression , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"
    
    ## 
    if isinstance ( data , ROOT.RooAbsData ) or ( cuts and as_weight ) :
        stat = Ostap.Math.WGeometricMean ()
    else : 
        stat = Ostap.Math. GeometricMean ()
        
    return data_get_stat  ( data       ,
                            stat       ,
                            expression ,
                            cuts       , *args     , 
                            cut_range  = cut_range ,
                            progress   = progress  ,
                            use_frame  = use_frame ,
                            parallel   = parallel  ) .value() 

    
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
def data_power_mean ( data , p ,
                      expression         ,
                      cuts       = ''    , *args , 
                      cut_range  = ''    ,
                      progress   = False , 
                      as_weight  = True  , ## interpret cuts as weiggt 
                      use_frame  = False ,
                      parallel   = False ) :                
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
    ## decode expressions & cuts
    var_lst, cuts, input_string = vars_and_cuts ( expression , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"

    if    p == -1 or isequal ( p , -1. ) :
        return data_harmonic_mean   ( data       ,
                                      expression , 
                                      cuts       , *args     , 
                                      cut_range  = cut_range ,
                                      progress   = progress  ,
                                      use_frame  = use_frame ,
                                      parallel   = parallel  ) 
    elif  p ==  0 or iszero  ( p       ) :
        return data_geometric_mean  ( data       ,
                                      expression , 
                                      cuts       , *args     ,
                                      cut_range  = cut_range ,
                                      progress   = progress  ,
                                      use_frame  = use_frame ,
                                      parallel   = parallel  ) 
    elif  p ==  1 or isequal ( p ,  1. ) :
        return data_arithmetic_mean ( data       ,
                                      expression ,
                                      cuts       , *args     ,
                                      cut_range  = cut_range ,
                                      progress   = progress  ,
                                      use_frame  = use_frame ,
                                      parallel   = parallel  ) 

    
    ## 
    if isinstance ( data , ROOT.RooAbsData ) or ( cuts and as_weight ) :
        stat = Ostap.Math.WPowerMean ( p )
    else : 
        stat = Ostap.Math. PowerMean ( p )
        
    return data_get_stat  ( data       , 
                            stat       ,
                            expression , 
                            cuts       , *args     , 
                            cut_range  = cut_range ,
                            progress   = progress  ,
                            use_frame  = use_frame ,
                            parallel   = parallel  ).value() 

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
def data_arithmetic_mean ( data , 
                           expression         ,
                           cuts       = ''    , *args , 
                           cut_range  = ''    ,
                           progress   = False , 
                           as_weight  = True  , ## interpret cuts as weiggt 
                           use_frame  = False ,
                           parallel   = False ) :                
    """ Get power mean over the data (just for completeness)
    >>> data = ...
    >>> result = data_arithmetic_mean ( data , 5 , 'pt' , 'eta>0' )
    - see `Ostap.Math::Moment`
    - see `Ostap.Math::WMoment`
    - see `Ostap.Math::ArithmeticMean`
    - see `Ostap.Math::WArithmeticMean`
    - see `Ostap.statVar.the_moment`
    """
    ##
    ## (1) decode expressions & cuts
    var_lst , cuts , _  = vars_and_cuts ( expression , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"
    
    if isinstance ( data , ROOT.RooAbsData ) or ( cuts and as_weight ) :
        stat = Ostap.Math.WArithmeticMean ()
    else : 
        stat = Ostap.Math. ArithmeticMean ()
        
    return data_get_stat  ( data                   ,
                            stat                   ,
                            expression             , 
                            cuts       , *args     , 
                            cut_range  = cut_range ,
                            progress   = progress  ,
                            use_frame  = use_frame ,
                            parallel   = parallel  ).value()  

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
def data_lehmer_mean ( data , p , 
                       expression         ,
                       cuts       = ''    , *args , 
                       cut_range  = ''    ,
                       progress   = False , 
                       as_weight  = True  , ## interpret cuts as weiggt 
                       use_frame  = False ,
                       parallel   = False ) :                    
    """ Get Lehmer mean over the data
    >>> data = ...
    >>> result = data_lehmer_mean ( data , 5 , 'pt' , 'eta>0' )
    - see `Ostap.Math::Moment`
    - see `Ostap.Math::WMoment`
    - see `Ostap.Math::LehmerMean`
    - see `Ostap.Math::WLehmerMean`
    - see `Ostap.statVar.the_moment`
    """
    
    ## (1) decode expressions & cuts
    var_lst , cuts , _  = vars_and_cuts ( expression , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"
    
    assert isinstance ( p , num_types ) , 'Invalid p-parameter!'
    
    if   p == 0 or iszero  ( p       ) :
        return data_harmonic_mean   ( data       ,
                                      expression , 
                                      cuts       , *args     , 
                                      cut_range  = cut_range ,
                                      progress   = progress  ,
                                      use_frame  = use_frame ,
                                      parallel   = parallel  )
    elif p == 1 or isequal ( p , 1.0 ) :
        return data_arithmetic_mean ( data       ,
                                      expression , 
                                      cuts       , *args     ,
                                      cut_range  = cut_range ,
                                      progress   = progress  ,
                                      use_frame  = use_frame ,
                                      parallel   = parallel  ) 
                                      

    if isinstance ( data , ROOT.RooAbsData ) or ( cuts and as_weight ) :
        stat = Ostap.Math.WLehmerMean ( p )
    else : 
        stat = Ostap.Math. LehmerMean ( p )

    return data_get_stat  ( data       ,
                            stat       ,
                            expression , 
                            cuts       , *args     ,
                            cut_range  = cut_range ,
                            progress   = progress  ,
                            use_frame  = use_frame ,
                            parallel   = parallel  ).value() 
        

# =============================================================================
## Get the mean (with uncertainty):
#  @code
#  data = ...
#  data_mean( data , 'mass*mass', 'pt>0')
#  data.mean(        'mass*mass', 'pt>0') ## ditto
#  @endcode 
#  @see Ostap::StatVar::moment
def data_mean ( data  , 
                expression         ,
                cuts       = ''    , *args , 
                cut_range  = ''    ,
                progress   = False , 
                as_weight  = True  , ## interpret cuts as weiggt 
                use_frame  = False ,
                parallel   = False ) :                                    
    """ Get the   mean (with uncertainty):
    >>> data = ...
    >>> data_mean( data , 'mass*mass', 'pt>0')
    >>> data.mean(        'mass*mass', 'pt>0') ## ditto
    """
    m2 = data_the_moment ( data , 2   ,
                           expression ,
                           cuts       , *args     ,
                           cut_range  = cut_range ,
                           progress   = progress  ,
                           use_frame  = use_frame ,
                           parallel   = parallel  ) 
    return m2.mean()

# =============================================================================
## Get the variance (with uncertainty):
#  @code
#  data = ...
#  data_mean( data , 'mass*mass', 'pt>0')
#  data.mean(        'mass*mass', 'pt>0') ## ditto
#  @endcode 
#  @see Ostap::StatVar::moment
def data_variance ( data  , 
                    expression         ,
                    cuts       = ''    , *args , 
                    cut_range  = ''    ,
                    progress   = False , 
                    as_weight  = True  , ## interpret cuts as weiggt 
                    use_frame  = False ,
                    parallel   = False ) :                                    
    """ Get the   mean (with uncertainty):
    >>> data = ...
    >>> data_mean( data , 'mass*mass', 'pt>0')
    >>> data.mean(        'mass*mass', 'pt>0') ## ditto
    """
    m4 = data_the_moment ( data   , 4 ,
                           expression , 
                           cuts       , *args     , 
                           cut_range  = cut_range ,
                           progress   = progress  ,
                           use_frame  = use_frame ,
                           parallel   = parallel  ) 
    return m4.variance ()

data_dispersion = data_variance


# =============================================================================
## Get the rms (with uncertainty):
#  @code
#  data = ...
#  data_mean( data , 'mass*mass', 'pt>0')
#  data.mean(        'mass*mass', 'pt>0') ## ditto
#  @endcode 
#  @see Ostap::StatVar::moment
def data_rms ( data               , 
               expression         ,
               cuts       = ''    , *args , 
               cut_range  = ''    ,
               progress   = False , 
               as_weight  = True  , ## interpret cuts as weiggt 
               use_frame  = False ,
               parallel   = False ) :                                    
    """ Get the   mean (with uncertainty):
    >>> data = ...
    >>> data_mean( data , 'mass*mass', 'pt>0')
    >>> data.mean(        'mass*mass', 'pt>0') ## ditto
    """
    variance  = data_variance ( data                   ,
                                expression             , 
                                cuts       , *args     , 
                                cut_range  = cut_range ,
                                progress   = progress  ,
                                use_frame  = use_frame ,
                                parallel   = parallel  ) 
    return variance ** 0.5 

# =============================================================================
## get the  skewness (with uncertainty)
#  @code
#  data =  ...
#  print data_skewness ( data , 'mass' , 'pt>1' ) 
#  print data.skewness (        'mass' , 'pt>1' ) ## ditto
#  @endcode
#  @see Ostap::StatVar::skewness
def data_skewness ( data  , 
                    expression         ,
                    cuts       = ''    , *args , 
                    cut_range  = ''    ,
                    progress   = False , 
                    as_weight  = True  , ## interpret cuts as weiggt 
                    use_frame  = False ,
                    parallel   = False ) :                                    
    """ Get the  skewness
    >>> data =  ...
    >>> print data_skewness ( data , 'mass' , 'pt>1' ) 
    >>> print data.skewness (        'mass' , 'pt>1' ) ## ditto
    - see Ostap::StatVar::skewness
    """
    result = data_the_moment ( data , 6 , expression , *args ,
                               cuts       = cuts      ,
                               cut_range  = cut_range ,
                               progress   = progress  ,
                               use_frame  = use_frame ,
                               parallel   = parallel  ) 
    return result.skewness ()

# =============================================================================
## get the (excess) kurtosis (with uncertainty)
#  @code
#  data =  ...
#  print data_kustosis ( data , 'mass' , 'pt>1' ) 
#  print data.kustosis (        'mass' , 'pt>1' ) ## ditto
#  @endcode
#  @see Ostap::StatVar::kurtosis
def data_kurtosis ( data               , 
                    expression         ,
                    cuts       = ''    , *args , 
                    cut_range  = ''    ,
                    progress   = False , 
                    as_weight  = True  , ## interpret cuts as weiggt 
                    use_frame  = False ,
                    parallel   = False ) :                                    
    """ Get the kurtosis
    >>> data =  ...
    >>> print data_kurtosis ( data , 'mass' , 'pt>1' ) 
    >>> print data.kurtosis (        'mass' , 'pt>1' ) ## ditto
    - see Ostap::StatVar::kurtosis
    """
    result = data_the_moment ( data , 8   ,
                               expression ,
                               cuts       , *args     ,
                               cut_range  = cut_range ,
                               progress   = progress  ,
                               use_frame  = use_frame ,
                               parallel   = parallel  ) 
    return result.kurtosis ()

# =============================================================================
## get the (approximate) quantiles for the data using P2-algorithm
#  @code
#  data =  ...
#  print ( data_quantiles ( data , 3 , 'mass' , 'pt>1' ) ) 
#  @endcode
#
#  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
#
#  @see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
def  data_quantiles ( data               ,
                      p                  , 
                      expression         ,
                      cuts       = ''    , *args , 
                      cut_range  = ''    ,
                      progress   = False , 
                      use_frame  = False ,
                      parallel   = False ) : 
    """ Get the (approximate) quantiles for the data using P2-algorithm 
    >>> data =  ...
    >>> data =  ...
    >>> print ( data_quantiles ( data , 3 , 'mass' , 'pt>1' ) ) 
    - see https://aakinshin.net/posts/p2-quantile-estimator-intro/
    - see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
    - see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
    - see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
    - see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
    """

    if   isinstance ( p , integer_types  ) and 1 <= p     : qq = Ostap.Math.Quantiles_[p] ()
    elif isinstance ( p , float          ) and 0 <  p < 1 : qq = Ostap.Math.Quantile  ( p  )    
    elif isinstance ( p , sequence_types ) and all ( isinstance ( v , float ) and 0 < v < 1 for v in p ) :
        qq = Ostap.Math.Quantiles ( doubles ( sorted ( float ( v ) for v in p ) ) ) 
    else :
        raise TypeError ( 'Invalid probabilities: %s/%s' % ( typename ( p ) , str ( p ) ) ) 
    
    result = data_get_stat ( data                  ,
                             qq                    ,
                             expression            ,
                             cuts       , *args    ,
                             cut_range = cut_range ,
                             progress  = progress  ,              
                             use_frame = False     ,  ## ATTENTION!
                             parallel  = False     )  ## ATTENTION!

    qq = result.quantiles()
    nq = len ( qq ) 
    return tuple ( qq [ i ] for i in range ( nq ) ) 

# =============================================================================
## get the (approximate) median for the data using P2-algorithm
#  @code
#  data =  ...
#  print data_median  ( data , 'mass' , 'pt>1' ) 
#  @endcode
#  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
#  @see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
def data_median ( data       ,
                  expression ,
                  cuts       = ""    , *args ,
                  cut_range  = ''    ,
                  progress   = False , 
                  use_frame  = False ,
                  parallel   = False ) :
    
    """ Get the (approximate) median for the data using P2-algorithm 
    >>> data =  ...
    >>> data =  ...
    >>> print data_median  ( data , 'mass' , 'pt>1' ) 
    - see https://aakinshin.net/posts/p2-quantile-estimator-intro/
    - see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
    - see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
    - see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
    - see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
    """
    
    qs = data_quantiles ( data       ,
                          0.5        ,
                          expression ,
                          cuts       , *args     ,
                          cut_range  = cut_range , 
                          progress   = progress  ,  
                          use_frame  = False     ,
                          parallel   = False     ) 
    return qs [ 1 ]

# =============================================================================
## Get the (approximate) terciles for the data using P2-algorithm
#  @code
#  data =  ...
#  print data_terciles ( data , 'mass' , 'pt>1' ) 
#  @endcode
#  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
#  @see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
def data_terciles ( data       ,
                    expression ,
                    cuts       = ""    , *args ,
                    cut_range  = ""    , 
                    progress   = False , 
                    use_frame  = False ,
                    parallel   = False ) :
    """ Get the (approximate) terciles for the data using P2-algorithm
    >>> data =  ...
    >>> print data_terciles ( data , 'mass' , 'pt>1' ) 
    """
    return data_quantiles ( data       ,
                            2          ,
                            expression ,
                            cuts       , *args     ,
                            cut_range  = cut_range , 
                            progress   = progress  ,  
                            use_frame  = False     ,
                            parallel   = False     ) 

# =============================================================================
## Get the (approximate) quartiles  for the data using P2-algorithm
#  @code
#  data =  ...
#  print data_quartiles ( data , 'mass' , 'pt>1' ) 
#  @endcode
#  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
#  @see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
def data_quartiles  ( data       ,
                      expression ,
                      cuts       = ""    , *args ,
                      cut_range  = ""    , 
                      progress   = False , 
                      use_frame  = False ,
                      parallel   = False ) :
    """ Get the (approximate) quartiles for the data using P2-algorithm
    >>> data =  ...
    >>> print data_quartiles ( data , 'mass' , 'pt>1' ) 
    """
    return data_quantiles ( data       ,
                            3          ,
                            expression ,
                            cuts       , *args     ,
                            cut_range  = cut_range , 
                            progress   = progress  ,  
                            use_frame  = False     ,
                            parallel   = False     ) 
    
# =============================================================================
## Get the (approximate) quintiles for the data using P2-algorithm
#  @code
#  data =  ...
#  print data_quintiles( data , 'mass' , 'pt>1' ) 
#  @endcode
#  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
#  @see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
def data_quintiles  ( data       ,
                      expression ,
                      cuts       = ""    , *args ,
                      cut_range  = ""    , 
                      progress   = False , 
                      use_frame  = False ,
                      parallel   = False ) :
    """ Get the (approximate) quintiles for the data using P2-algorithm
    >>> data =  ...
    >>> print data_quartiles ( data , 'mass' , 'pt>1' ) 
    """
    return data_quantiles ( data       ,
                            4          ,
                            expression ,
                            cuts       , *args     ,
                            cut_range  = cut_range , 
                            progress   = progress  ,  
                            use_frame  = False     ,
                            parallel   = False     ) 
    

# =============================================================================
## Get the (approximate) sextiles  for the data using P2-algorithm
#  @code
#  data =  ...
#  print data_quintiles( data , 'mass' , 'pt>1' ) 
#  @endcode
#  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
#  @see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
def data_sextiles   ( data       ,
                      expression ,
                      cuts       = ""    , *args ,
                      cut_range  = ""    , 
                      progress   = False , 
                      use_frame  = False ,
                      parallel   = False ) :
    """ Get the (approximate) sextile  for the data using P2-algorithm
    >>> data =  ...
    >>> print data_sextiles ( data , 'mass' , 'pt>1' ) 
    """
    return data_quantiles ( data       ,
                            5          ,
                            expression ,
                            cuts       , *args     ,
                            cut_range  = cut_range , 
                            progress   = progress  ,  
                            use_frame  = False     ,
                            parallel   = False     ) 
    

# =============================================================================
## Get the (approximate) septiles  for the data using P2-algorithm
#  @code
#  data =  ...
#  print data_septiles ( data , 'mass' , 'pt>1' ) 
#  @endcode
#  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
#  @see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
def data_septiles   ( data       ,
                      expression ,
                      cuts       = ""    , *args ,
                      cut_range  = ""    , 
                      progress   = False , 
                      use_frame  = False ,
                      parallel   = False ) :
    """ Get the (approximate) septiles for the data using P2-algorithm
    >>> data =  ...
    >>> print data_septiles ( data , 'mass' , 'pt>1' ) 
    """
    return data_quantiles ( data       ,
                            6          ,
                            expression ,
                            cuts       , *args     ,
                            cut_range  = cut_range , 
                            progress   = progress  ,  
                            use_frame  = False     ,
                            parallel   = False     ) 
    

# =============================================================================
## Get the (approximate) octiles for the data using P2-algorithm
#  @code
#  data =  ...
#  print data_octiles ( data , 'mass' , 'pt>1' ) 
#  @endcode
#  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
#  @see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
def data_octiles    ( data       ,
                      expression ,
                      cuts       = ""    , *args ,
                      cut_range  = ""    , 
                      progress   = False , 
                      use_frame  = False ,
                      parallel   = False ) :
    """ Get the (approximate) octiles for the data using P2-algorithm
    >>> data =  ...
    >>> print data_octiles ( data , 'mass' , 'pt>1' ) 
    """
    return data_quantiles ( data       ,
                            7          ,
                            expression ,
                            cuts       , *args     ,
                            cut_range  = cut_range , 
                            progress   = progress  ,  
                            use_frame  = False     ,
                            parallel   = False     ) 
    

# =============================================================================
## Get the (approximate) deciles for the data using P2-algorithm
#  @code
#  data =  ...
#  print data_deciles ( data , 'mass' , 'pt>1' ) 
#  @endcode
#  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
#  @see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
def data_deciles    ( data       ,
                      expression ,
                      cuts       = ""    , *args ,
                      cut_range  = ""    , 
                      progress   = False , 
                      use_frame  = False ,
                      parallel   = False ) :
    """ Get the (approximate) deciles for the data using P2-algorithm
    >>> data =  ...
    >>> print data_deciles ( data , 'mass' , 'pt>1' ) 
    """
    return data_quantiles ( data       ,
                            9          ,
                            expression ,
                            cuts       , *args     ,
                            cut_range  = cut_range , 
                            progress   = progress  ,  
                            use_frame  = False     ,
                            parallel   = False     ) 
    
# =============================================================================
## Get the (approximate) ventiles for the data using P2-algorithm
#  @code
#  data =  ...
#  print data_ventiles ( data , 'mass' , 'pt>1' ) 
#  @endcode
#  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
#  @see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
def data_ventiles   ( data       ,
                      expression ,
                      cuts       = ""    , *args ,
                      cut_range  = ""    , 
                      progress   = False , 
                      use_frame  = False ,
                      parallel   = False ) :
    """ Get the (approximate) ventiles for the data using P2-algorithm
    >>> data =  ...
    >>> print data_ventiles ( data , 'mass' , 'pt>1' ) 
    """
    return data_quantiles ( data       ,
                            19         ,
                            expression ,
                            cuts       , *args     ,
                            cut_range  = cut_range , 
                            progress   = progress  ,  
                            use_frame  = False     ,
                            parallel   = False     ) 
    

# =============================================================================
## Get the (approximate) percentiles for the data using P2-algorithm
#  @code
#  data =  ...
#  print data_percentiles ( data , 'mass' , 'pt>1' ) 
#  @endcode
#  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-adjusting-order/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-initialization/
#  @see https://aakinshin.net/posts/p2-quantile-estimator-rounding-issue/
#  @see https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
def data_percentiles ( data       ,
                       expression ,
                       cuts       = ""    , *args ,
                       cut_range  = ""    , 
                       progress   = False , 
                       use_frame  = False ,
                       parallel   = False ) :
    """ Get the (approximate) percentiles for the data using P2-algorithm
    >>> data =  ...
    >>> print data_percentiles ( data , 'mass' , 'pt>1' ) 
    """
    return data_quantiles ( data       ,
                            99         ,
                            expression ,
                            cuts       , *args     ,
                            cut_range  = cut_range , 
                            progress   = progress  ,  
                            use_frame  = False     ,
                            parallel   = False     ) 
    

# =============================================================================
## Get "center" for data  using Pragmastat toolkit 
#  @see https://pragmastat.dev/
def data_center ( data       ,
                  expression ,
                  cuts       = ''    , *args , 
                  cut_range  = ''    , 
                  progress   = False ,  
                  use_frame  = False ,
                  parallel   = False ) :
    """ Get "center" for data using Pragmastat toolkit 
    - see https://pragmastat.dev/
    """
    
    ## (1) decode expressions & cuts
    var_lst , cuts , _  = vars_and_cuts ( expression , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"
    
    if cuts : logger.warning ( "selection cuts      will be used only as boolean!" )
    if isinstance ( data , ROOT.RooAbsData ) and data.isWeighted () :
        logger.warning ( "Weight from dataset will be used only as boolean!" )

    npdata , weights = data_slice  ( data                   ,
                                     expression             ,
                                     cuts       , *args     ,
                                     cut_range  = cut_range ,
                                     structured = False     ,
                                     transpose  = False     ,
                                     progress   = progress  , 
                                     use_frame  = use_frame ,
                                     parallel   = parallel  )
    
    import pragmastat as PS 
    return PS.center ( npdata )

# =============================================================================
## Get "spread" for data  using Pragmastat toolkit 
#  @see https://pragmastat.dev/
def data_spread ( data       ,
                  expression ,
                  cuts       = ''    , *args , 
                  cut_range  = ''    , 
                  progress   = False ,  
                  use_frame  = False ,
                  parallel   = False ) :
    """ Get "spread" for data  using Pragmastat toolkit 
    - see https://pragmastat.dev/
    """
    
    ## (1) decode expressions & cuts
    var_lst , cuts , _  = vars_and_cuts ( expression , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"
    
    if cuts : logger.warning ( "selection cuts      will be used only as boolean!" )
    if isinstance ( data , ROOT.RooAbsData ) and data.isWeighted () :
        logger.warning ( "Weight from dataset will be used only as boolean!" )

    npdata , weights = data_slice  ( data                   ,
                                     expression             ,
                                     cuts       , *args     ,
                                     cut_range  = cut_range ,
                                     structured = False     ,
                                     transpose  = False     ,
                                     progress   = progress  , 
                                     use_frame  = use_frame ,
                                     parallel   = parallel  )
    
    import pragmastat as PS 
    return PS.spread ( npdata )

# =============================================================================
## Get "volatility" for data  using Pragmastat toolkit: spread/abs(center)
#  @see https://pragmastat.dev/
def data_volatility ( data       ,
                      expression ,
                      cuts       = ''    , *args , 
                      cut_range  = ''    , 
                      progress   = False ,  
                      use_frame  = False ,
                      parallel   = False ) :
    """ Get "volatility" for data using Pragmastat toolkit 
    - see https://pragmastat.dev/
    """
    
    ## (1) decode expressions & cuts
    var_lst , cuts , _  = vars_and_cuts ( expression , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"
    
    if cuts : logger.warning ( "selection cuts      will be used only as boolean!" )
    if isinstance ( data , ROOT.RooAbsData ) and data.isWeighted () :
        logger.warning ( "Weight from dataset will be used only as boolean!" )
        
    npdata , weights = data_slice  ( data                   ,
                                     expression             ,
                                     cuts       , *args     ,
                                     cut_range  = cut_range ,
                                     structured = False     ,
                                     transpose  = False     ,
                                     progress   = progress  , 
                                     use_frame  = use_frame ,
                                     parallel   = parallel  )
    
    import pragmastat as PS 
    return PS.volatility ( npdata )

# =============================================================================
## Get "precision" for data  using Pragmastat toolkit: 2*spread/sqrt(n) 
#  @see https://pragmastat.dev/
def data_precision ( data       ,
                     expression ,
                     cuts       = ''    , *args , 
                     cut_range  = ''    , 
                     progress   = False ,  
                     use_frame  = False ,
                     parallel   = False ) :
    """ Get "precision" for data  using Pragmastat toolkit: 2*spread/sqrt(n) 
    - see https://pragmastat.dev/
    """
    
    ## (1) decode expressions & cuts
    var_lst , cuts , _  = vars_and_cuts ( expression , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"
    
    if cuts : logger.warning ( "selection cuts      will be used only as boolean!" )
    if isinstance ( data , ROOT.RooAbsData ) and data.isWeighted () :
        logger.warning ( "Weight from dataset will be used only as boolean!" )
        
    npdata , weights = data_slice  ( data                   ,
                                     expression             ,
                                     cuts       , *args     ,
                                     cut_range  = cut_range ,
                                     structured = False     ,
                                     transpose  = False     ,
                                     progress   = progress  , 
                                     use_frame  = use_frame ,
                                     parallel   = parallel  )
    
    import pragmastat as PS 
    return PS.precision ( npdata )

# =============================================================================
## helper function to copy/clone the projection target 
#  - clone for histograms  
#  - copy constructor for other (C++) types
# 
#  Should we use here more generic (python) copy ?
def target_copy   ( t ) :
    """ Helper function to copy/clone the projection target 
    - clone for histograms  
    - copy constructor for other types 
    """
    if isinstance ( t , ROOT.TH1 ) :
        c = t.Clone ()
        if not c.GetSumw2() : c.Sumw2()
        return c
    ## copy constructor for other C++ objects
    ## Should we use here more generic (python) copy ?
    c = type ( t ) ( t )
    ## reset it! 
    return target_reset ( c )

# =============================================================================
## Helper function to reset the projection target
#  - `TH1.Reset` for histogram 
#  - `reset`     for other objects 
def target_reset  ( t ) :
    """ Helper function to reset the projection target
    - `TH1.Reset` for histogram 
    - `reset`     for other objects 
    """
    if isinstance ( t , ROOT.TH1 ) : 
        t.Reset ()
        if not t.GetSumw2() : t.Sumw2 () 
    else : t.reset ()
    ## 
    return t

# =============================================================================
## progect data (TTree or RooAbsData) into into the histogram or some other object
#  @see Ostap::Project
def data_project ( data                ,
                   target              ,
                   expressions         ,
                   cuts        = ''    , *args , 
                   cut_range   = ""    ,
                   as_weight   = ""    ,
                   progress    = True  ,
                   use_frame   = False ,
                   parallel    = False ) :
    
    """ project data (TTree or RooAbsData) into into histogram or some other object
    - see `Ostap,Project`
    """
    
    # ========================================================================
    ## (1) decode expressions & cuts
    # ========================================================================
    var_lst , cuts , _  = vars_and_cuts ( expressions , cuts )
    nvars = len ( var_lst )

    # ========================================================================
    ## profile ?
    # ========================================================================
    from ostap.histos.histos import profile_types 
    profile = isinstance ( target , profile_types )  
    
    # ========================================================================
    ## (2) check consistency
    # ========================================================================
    h1_stack = False ## stack of 1D-objects/histograms  
    if   4 == nvars and isinstance ( target , _s4D            ) : pass
    elif 3 == nvars and isinstance ( target , _s3D            ) : pass
    elif 2 == nvars and isinstance ( target , _s2D            ) : pass
    elif 1 == nvars and isinstance ( target , _s1D            ) : pass
    ## profiles 
    elif 4 == nvars and isinstance ( target , ROOT.TProfile3D ) and 3 == target.dim() : pass
    elif 3 == nvars and isinstance ( target , ROOT.TProfile2D ) and 2 == target.dim() : pass
    elif 2 == nvars and isinstance ( target , ROOT.TProfile   ) and 1 == target.dim() : pass
    ## histgrams 
    elif 3 == nvars and isinstance ( target , ROOT.TH3        ) and 3 == target.dim() : pass    
    elif 2 == nvars and isinstance ( target , ROOT.TH2        ) and 2 == target.dim() : pass
    ## stack...
    elif 1 == nvars and isinstance ( target , ROOT.TH1        ) and 1 == target.dim() : pass
    elif 1 <  nvars and 1 == target.dim() and not profile  : h1_stack = True 
    else :
        raise TypeError ( 'Target %s and expression(s): %s are inconsistent' % \
                          ( typename ( target ) , ','.join ( var_lst ) ) )
    
    # ========================================================================
    ## (3) cut_range defined *only* for RooFit datasets 
    # ========================================================================
    if cut_range and not isinstance ( data , ROOT.RooAbsData ) : 
        raise TypeError ( "Invalid use of `cut_range':%s for %s " % ( cut_range , typename ( data ) ) )
    
    # ========================================================================
    ## (4) display progress ? 
    # ========================================================================
    progress = progress_conf ( progress )

    # ========================================================================
    ## (5) create the driver 
    # ========================================================================
    pv = Ostap.Project ( progress )
    
    # ========================================================================
    ## (6) RooFit dataset ?
    # ========================================================================
    if isinstance ( data , ROOT.RooAbsData ) :
        
        # ====================================================================
        ## weight with errors ? 
        # ====================================================================
        if data.isWeighted () and Ostap.Utils.storeError ( data ) :
            if isinstance ( target , ROOT.TProfile   ) or \
               isinstance ( target , ROOT.TProfile2D ) or not isinstance ( target , ROOT.TH1 ) :
                logger.warning ( 'Dataset has weight with errors! Weight errors will be ignored!' ) 
                
        with rootException() :
            # =================================================================
            ## project several variables into the same 1D histogram 
            # =================================================================
            if h1_stack :

                target = target_reset ( target ) 
                tcopy  = target_copy  ( target ) 
                the_args = ( cuts , cut_range ) + args
                for var in var_lst :
                    sc = pv.project1 ( data , tcopy , var , *the_args )
                    assert sc.isSuccess() , 'Error %s from Ostap::Project::project1(%s,%s)' % ( sc , typename ( target ) , var )  
                    target += tcopy 
                del tcopy
                ## 
                return target                                         ## RETURN
            
            # =================================================================
            ## regular processing 
            # =================================================================
            the_args = var_lst + ( cuts , cut_range ) + args
            if   1 == nvars : sc = pv.project1 ( data , target , *the_args )
            elif 2 == nvars : sc = pv.project2 ( data , target , *the_args )
            elif 3 == nvars : sc = pv.project3 ( data , target , *the_args )            
            elif 4 == nvars : sc = pv.project4 ( data , target , *the_args )            
            assert sc.isSuccess() , \
                'Error %s from Ostap::Project::project%d(%s,%s)' % ( sc                   ,
                                                                     typename ( target )  ,
                                                                     nvars                ,
                                                                     ','.join ( var_lst ) )
            ## 
            return target                                            ## RETURN


    # =========================================================================
    ## (7) delegate processing to RDataFrame if/when  possible?
    # =========================================================================
    if  use_frame and good_for_frame ( data , *args ) and not h1_stack and not isinstance ( target , ROOT.TProfile3D ) : 
        return F.frame_project ( data                   ,
                                 model       = target   ,
                                 expressions = var_lst  ,
                                 cuts        = cuts     ,
                                 progress    = progress ,
                                 report      = progress ,
                                 lazy        = False    ) 

    # =========================================================================
    ## (8) parallel processing ?
    # =========================================================================
    if  parallel and good_for_parallel ( data , *args ) and not profile : 
        from ostap.parallel.parallel_statvars import parallel_project
        return parallel_project ( data                   ,
                                  target                 ,
                                  expressions            ,
                                  cuts       , *args     ,
                                  as_weight  = as_weight ,
                                  progress   = progress  ,
                                  use_frame  = use_frame )

    # =========================================================================
    ## (9) regular TTree processing
    # =========================================================================
    ## Branches to be activated
    from ostap.trees.trees import ActiveBranches
    with rootException() , ActiveBranches ( data , cuts , *var_lst ) :
        the_args = var_lst + ( cuts , ) + args         
        if   1 == nvars : sc = pv.project1 ( data , target , *the_args  )
        elif 2 == nvars : sc = pv.project2 ( data , target , *the_args  )
        elif 3 == nvars : sc = pv.project3 ( data , target , *the_args  )
        elif 4 == nvars : sc = pv.project4 ( data , target , *the_args  )
        assert sc.isSuccess() , 'Error %s from StatVar::project(1,2,3,4)' % sc 
        return target

# =============================================================================
## Get slice of dat s in form of Numpy array
#  @code
#  data = ...
#  arr , weight = data_slice ( data , "x,y,x" , "pt>1" ) 
#  @endif 
# =============================================================================
def data_slice ( data        ,
                 expressions ,
                 cuts        = ""    , *args , 
                 cut_range   = ""    ,
                 structured  = True  ,
                 transpose   = True  ,
                 progress    = False , 
                 use_frame   = False ,
                 parallel    = False ) :
    """ Get slice of data s in form of Numpy array
    >>> data = ...
    >>> arr , weight = data_slice ( data , "x,y,x" , "pt>1" ) 
    """
    
    ## (1) decode expressions & cuts
    var_lst , cuts, _  = vars_and_cuts ( expressions , cuts )
    
    ## (2) adjust first/last 
    first , last = evt_range ( data , *args )

    ## (3) cut_range defined *only* for RooFit datasets 
    if cut_range and not isinstance ( data , ROOT.RooAbsData ) : 
        raise TypeError ( "Invalid use of `cut_range':%s" % cut_range  ) 
    
    ## (3) RooFit ?
    if isinstance ( data , ROOT.RooAbsData ) :
        from ostap.fitting.dataset import ds_slice
        return ds_slice ( data                    ,
                          expressions             ,
                          cuts       = cuts       ,
                          cut_range  = cut_range  ,
                          structured = structured ,
                          transpose  = transpose  , 
                          first      = first      ,
                          last       = last       )

    ## Frame processing ? 
    if  use_frame and good_for_frame ( data , first , last  ) : 
        return F.frame_slice ( data                    ,
                               expression              ,
                               cuts       = cuts       ,
                               structured = structured ,
                               transpose  = transpose  , 
                               progress   = prorgess   )
    ## Parallel processing?
    elif parallel   and good_for_parallel ( data , first , last ) :
        from ostap.parallel.parallel_statvar import parallel_slice 
        return parallel_slice ( data                     ,
                                expressions              ,
                                cuts        , first      , last ,
                                structured  = structured ,
                                transpose   = transpose  , 
                                progress    = prorgess   ,
                                use_frame   = use_frame  ) 

    assert isinstance ( data , ROOT.TTree ) , "Here data must be TTree: %s" % typename ( data ) 

    from ostap.trees.trees import tree_slice
    return tree_slice ( data        ,
                        expressions ,
                        cuts        ,
                        first       = first      ,
                        last        = last       , 
                        structured  = structured ,
                        transpose   = transpose  , 
                        progress    = progress   , 
                        use_frame   = False      ,
                        parallel    = False      )

# =============================================================================
## Produce  "efficiency" histogram for boolean <code>criteriaon</c>
#  as function of valiabed listed as <code>expressions</code>
#  internally it creatd two histogram
#  - "accepted" for events accepted by (boolean) criterion
#  - "rejected" for events rekected by (boolean) criterion
#
#  @code
#  histo_1D = ...
#  data     = ...
#  eff, accepted, rejected = data_efficiency
#  ...   ( tree         ,  
#  ...     'DLL>5'      , 
#  ...     histo_1D     ,
#  ...     'PT'         ,
#  ...     cuts = 'A>2' ) 
#  @code
#
#  @code
#  histo_2D = ...
#  data     = ...
#  eff, accepted, rejected = data_efficiency
#  ...   ( tree         ,  
#  ...     'DLL>5'      , 
#  ...     histo_2D     ,
#  ...     'PT, y'      , 
#  ...     cuts = 'A>2' ) 
#  @code
#
#  @code
#  histo_3D = ...
#  data     = ...
#  eff, accepted, rejected = data_efficiency
#  ...   ( tree             ,  
#  ...     'DLL>5'          , 
#  ...     histo_3D         ,
#  ...     'PT, y, nTracks' , 
#  ...     cuts = 'A>2'     ) 
#  @code
#
#  @param tree       (INPUT)  input data
#  @param criterion  (INPUT)  (boolean) criterion
#  @param histo      (UPDATE) outptu efficiency histogram 
#  @param exressions (INPUT)  expressions for the histiogram axes
#  @param cuts       (INPUT)  (boolean) selection criteria to be applied 
#  @param weight     (INPUT)  expression to be used as weight 
#  @return triplet of histigrams: efficiency, accepted & rejected 
#
#  @attention both `cuts` and `criterion` are treated as boolean! 
#  @see tree_project  
#  @see frame_project  
def data_efficiency ( data        ,
                      criterion   ,  
                      histo       , 
                      expressions ,
                      cuts        = ''          , *args ,
                      cut_range   = ''          , 
                      weight      = ''          , 
                      use_frame   = False       , 
                      parallel    = False       , 
                      progress    = False       ) : 
    """ Produce  "efficiency" histogram for boolean <code>criteriaon</c>
    as function of valiabed listed as <code>expressions</code>
    internally it creatd two histogram
    - "acepted" for events accepted by (boolean) criterion
    - "rejected" for events rekected by (boolean) criterion 
    
    tree:       (INPUT)  input tree 
    criterion:  (INPUT)  (boolean) criterion
    histo:      (UPDATE) outptu efficiency histogram 
    exressions: (INPUT)  expressions for the histiogram axes
    cuts:       (INPUT)  (boolean) selection criteria to be applied 
    weight:     (INPUT)  expression to be used as weight 

    return triplet of histigrams: 
        - efficiency
        - distribution for accepted events
        - distribution for rejected events
        
    ATTENTION: both `cuts` and `criterion` are treated as boolean!  
    
    - see `ostap.trees.trees.tree_project` 

    >>> histo_1D = ...
    >>> tree     = ...
    >>> eff, accepted, rejected = tree_efficiency
    ...    ( tree               ,  
    ...      'DLL>5'            , ## criterion 
    ...      histo_1D           , 
    ...      'PT'               , ## axes 
    ...      cuts   = 'A>2'     , 
    ...      weight = 'sWeight' )

    >>> histo_2D = ...
    >>> tree     = ...
    >>> eff, accepted, rejected = tree_efficiency
    ...    ( tree               ,  
    ...      'DLL>5'            , ## criterion 
    ...      histo_2D           , 
    ...      'PT, y'            , ## axes 
    ...      cuts   = 'A>2'     , 
    ...      weight = 'sWeight' )

    >>> histo_3D = ...
    >>> tree     = ...
    >>> eff, accepted, rejected = tree_efficiency
    ...    ( tree               ,  
    ...      'DLL>5'            , ## criterion 
    ...      histo_3D           , 
    ...      'PT, y, nTracks'   , ## axes 
    ...      cuts   = 'A>2'     , 
    ...      weight = 'sWeight' )

    """
    
    ## (1) decode expressions & cuts
    var_lst , cuts      , _  = vars_and_cuts ( expressions , cuts      )    
    ## (2) decode expressions & criterion 
    var_lst , criterion , _  = vars_and_cuts ( var_lst     , criterion )
    ## (3) decode expressions & weight 
    var_lst , weight    , _  = vars_and_cuts ( var_lst     , weight    )

    nvars  = len ( var_lst )
    
    assert criterion, "Invalid criterion: %s" % criterion 
    assert isinstance ( histo , ROOT.TH1 ) , "Invalid `histo` type: %s" % typename ( histo )
    
    hdim = histo.GetDimension() 
    assert 1 <= hdim <= 3 and hdim == nvars , \
        "Mismatch histogram dimension/#vars: %s/%d"% ( hdim, nvars ) 
    
    ## (4) cut_range defined *only* for RooFit datasets 
    if cut_range and not isinstance ( data , ROOT.RooAbsData ) : 
        raise TypeError ( "Invalid use of `cut_range':%s" % cut_range  ) 
    
    ## (5) RooFit ?
    if isinstance ( data , ROOT.RooAbsData ) :
        
        histo.Reset() 
        if not histo.GetSumw2() : histo.Sumw2() 
        
        h_efficiency = histo.clone ()
        h_accepted   = histo.clone ()
        h_rejected   = histo.clone ()
        
        if not h_efficiency.GetSumw2() : h_efficiency.Sumw2() 
        if not h_rejected  .GetSumw2() : h_rejected  .Sumw2() 
        if not h_rejected  .GetSumw2() : h_rejected  .Sumw2() 
        
        h_accepted   .SetTitle ( "Distribution for events `accepted` by %s" % criterion ) 
        h_rejected   .SetTitle ( "Distribution for events `rejected` by %s" % criterion )
        h_efficiency .SetTitle ( "Efficiency   for criterion: %s"           % criterion )
        
        ## use cuts as boolean!     
        the_cut  = ROOT.TCut( '!!(%s)' % cuts ) if cuts else ROOT.TCut() 
        
        ## use criterion as boolean 
        accept = the_cut * ( "!!(%s)" % criterion ) ## times boolean  
        reject = the_cut * (  "!(%s)" % criterion ) ## times boolean
        
        if weight : 
            accept *= weight  
            reject *= weight 

        ## fill 1st histogram 
        h_accept = data_project ( data       , 
                                  h_accepted , 
                                  var_lst    , 
                                  cuts       = accept    , *args , 
                                  cut_range  = cut_range , 
                                  use_frame  = False     ,
                                  parallel   = False     , 
                                  progress   = progress  ) 
            
        ## fill 2nd histogram 
        h_reject = data_project ( tree       , 
                                  h_rejected , 
                                  var_lst    , 
                                  cuts       = reject    , *args , 
                                  cut_range  = cut_range ,                                       
                                  use_frame  = use_frame ,
                                  parallel   = parallel  , 
                                  progress   = progress  )
        
        ## calculate binomial efficiency and assign it to the `histo``
        h_efficiency += 1 / ( 1 + h_rejected / h_accepted )            
        histo        += h_efficiency
        return h_efficiency, h_accepted, h_rejected

    ## unpack&update  first/last evens 
    first , last = evt_range ( data , *args[:2] )
    
    if  use_frame and good_for_frame ( data , first , last ) : 
        return F.frame_efficiency ( data                    ,
                                    criterion   = criterion ,
                                    histo       = histo     ,
                                    expressions = var_lst   , 
                                    cuts        = cuts      ,
                                    weight      = weight    ,
                                    progress    = progress  ,
                                    report      = False     ,
                                    lazy        = False     )
    
    ## delegate processing to parallel machinery if requested and possible 
    elif parallel and good_for_parallel ( data , first , last ) :
        from ostap.parallel.parallel_project import parallel_efficiency 
        return parallel_efficiency ( data      ,
                                     criterion ,
                                     histo     ,
                                     var_lst   ,
                                     cuts      , first     , last ,   
                                     weight    = weight    ,
                                     use_frame = use_frame , 
                                     progress  = progress  )

    assert isinstance ( data , ROOT.TTree ) , "Here data must be TTree: %s" % typename ( data ) 
    
    from ostap.trees.trees import tree_efficiency
    return tree_efficiency ( data        ,
                             criterion   ,
                             histo       , 
                             expressions ,
                             cuts        ,
                             weight      = weight     , 
                             first       = first      ,
                             last        = last       , 
                             progress    = progress   , 
                             use_frame   = use_frame  ,
                             parallel    = parallel   )

# =============================================================================
## decorate certain class with some useful  methods 
def data_decorate ( klass ) :
    """ Decorate certain class with some useful  methods
    """
    
    if hasattr ( klass , 'get_moment'     ) : klass.orig_get_moment     = klass.get_moment
    if hasattr ( klass , 'the_moment'     ) : klass.orig_the_moment     = klass.the_moment
    
    if hasattr ( klass , 'hasEntry'       ) : klass.orig_hasEntry       = klass.hasEntry 
    if hasattr ( klass , 'has_entry'      ) : klass.orig_has_entry      = klass.has_entry 
    
    if hasattr ( klass , 'ECDF'           ) : klass.orig_ECDF           = klass.ECDF    
    if hasattr ( klass , 'moment'         ) : klass.orig_moment         = klass.moment
    if hasattr ( klass , 'nEff'           ) : klass.orig_nEff           = klass.nEff
    if hasattr ( klass , 'mean'           ) : klass.orig_mean           = klass.mean
    if hasattr ( klass , 'variance'       ) : klass.orig_variance       = klass.variance 
    if hasattr ( klass , 'dispersion'     ) : klass.orig_dispersion     = klass.dispersion
    if hasattr ( klass , 'rms'            ) : klass.orig_rms            = klass.rms 
    if hasattr ( klass , 'skewness'       ) : klass.orig_skewness       = klass.skewness
    if hasattr ( klass , 'kurtosis'       ) : klass.orig_kurtosis       = klass.kurtosis

    if hasattr ( klass , 'sumVar'         ) : klass.orig_sumVar         = klass.sumVar
    if hasattr ( klass , 'sumVars'        ) : klass.orig_sumVars        = klass.sumVars
    if hasattr ( klass , 'statVar'        ) : klass.orig_statVar        = klass.statVar
    if hasattr ( klass , 'statVars'       ) : klass.orig_statVars       = klass.statVars
    if hasattr ( klass , 'statistic'      ) : klass.orig_statistic      = klass.statistic
    if hasattr ( klass , 'statCov'        ) : klass.orig_statCov        = klass.statCovs
    if hasattr ( klass , 'statCovs'       ) : klass.orig_statCovs       = klass.statCovs

    if hasattr ( klass , 'statVct'        ) : klass.orig_statVct        = klass.statVct
    if hasattr ( klass , 'statvector'     ) : klass.orig_statvector     = klass.statvector
    
    if hasattr ( klass , 'quantiles'      ) : klass.orig_quanitles      = klass.quantiles 
    if hasattr ( klass , 'terciles'       ) : klass.orig_terciles       = klass.terciles 
    if hasattr ( klass , 'quartiles'      ) : klass.orig_quartiles      = klass.quartiles 
    if hasattr ( klass , 'quintiles'      ) : klass.orig_quintiles      = klass.quintiles
    if hasattr ( klass , 'sextiles'       ) : klass.orig_sextiles       = klass.sextiles 
    if hasattr ( klass , 'septiles'       ) : klass.orig_septiles       = klass.septiles 
    if hasattr ( klass , 'octiles'        ) : klass.orig_octiles        = klass.octiles 
    if hasattr ( klass , 'deciles'        ) : klass.orig_deciles        = klass.deciles 
    if hasattr ( klass , 'ventiles'       ) : klass.orig_ventiles       = klass.ventiles  
    if hasattr ( klass , 'percentiles'    ) : klass.orig_percentiles    = klass.percentiles 

    if hasattr ( klass , 'median'         ) : klass.orig_median         = klass.median 

    if hasattr ( klass , 'center'         ) : klass.orig_center         = klass.center
    if hasattr ( klass , 'spread'         ) : klass.orig_spread         = klass.spread
    if hasattr ( klass , 'volatility'     ) : klass.orig_volatility     = klass.volatility
    if hasattr ( klass , 'precision'      ) : klass.orig_precision      = klass.precision
    
    if hasattr ( klass , 'efficiency'     ) : klass.orig_efficiency     = klass.efficiency 

    
    klass.get_moment      = data_the_moment
    klass.the_moment      = data_the_moment

    klass.hasEntry        = data_hasEntry 
    klass.has_entry       = data_hasEntry
    
    klass.data_size       = data_size 
    klass.num_entries     = data_size 
    klass.n_entries       = data_size 

    klass.ECDF            = data_ECDF 
    klass.moment          = data_moment
    klass.nEff            = data_nEff 
    klass.mean            = data_mean
    klass.variance        = data_variance 
    klass.dispersion      = data_dispersion
    klass.rms             = data_rms
    klass.skewness        = data_skewness
    klass.kurtosis        = data_kurtosis
    
    klass.sumVar          = data_sum 
    klass.sumVars         = data_sum 

    klass.statVar         = data_statistic
    klass.statVars        = data_statistic
    klass.statistic       = data_statistic

    klass.statCov         = data_covariance
    klass.statCovs        = data_covariance
    klass.statvector      = data_statvector
    klass.statVct         = data_statvector

    klass.quantiles       = data_quantiles 

    klass.median          = data_median
    klass.terciles        = data_terciles 
    klass.quartiles       = data_quartiles 
    klass.quintiles       = data_quintiles  
    klass.sextiles        = data_sextiles 
    klass.septiles        = data_septiles 
    klass.octiles         = data_octiles 
    klass.deciles         = data_deciles  
    klass.ventiles        = data_ventiles  
    klass.percentiles     = data_percentiles 
    
    klass.center          = data_center
    klass.spread          = data_spread
    klass.volatilty       = data_volatility
    klass.precision       = data_precision
    
    klass.efficiency      = data_efficiency 

    if hasattr ( klass , 'var_minmax'      ) : klass.orig_var_minmax      = klass.var_minmax 
    if hasattr ( klass , 'var_range'       ) : klass.orig_var_range       = klass.var_range 
        
    klass.var_minmax  = data_minmax 
    klass.var_range   = data_range 
    
    if hasattr ( klass , 'harmonic_mean'   ) : klass.orig_harmonic_mean   = klass.harmonic_mean
    if hasattr ( klass , 'geometric_mean'  ) : klass.orig_geometric_mean  = klass.geometric_mean
    if hasattr ( klass , 'power_mean'      ) : klass.orig_power_mean      = klass.power_mean
    if hasattr ( klass , 'lehmer_mean'     ) : klass.orig_lehmer_mean     = klass.lehmer_mean
    if hasattr ( klass , 'arithmetic_mean' ) : klass.orig_arithmetic_mean = klass.arithmetic_mean
    
    klass.harmonic_mean   = data_harmonic_mean 
    klass.geometric_mean  = data_geometric_mean 
    klass.power_mean      = data_power_mean 
    klass.lehmer_mean     = data_lehmer_mean 
    klass.arithmetic_mean = data_arithmetic_mean 

    return ( klass.get_moment      , 
             klass.the_moment      , 
             klass.hasEntry        , 
             klass.has_entry       , 
             klass.ECDF            , 
             klass.moment          , 
             klass.nEff            , 
             klass.mean            ,             
             klass.variance        , 
             klass.dispersion      , 
             klass.rms             , 
             klass.skewness        , 
             klass.kurtosis        ,             
             klass.statVar         ,
             klass.statVars        ,
             klass.statistic       ,             
             klass.sumVar          ,
             klass.sumVars         ,             
             klass.statCov         ,
             klass.statCovs        ,
             ##
             klass.harmonic_mean   , 
             klass.geometric_mean  , 
             klass.power_mean      ,
             klass.lehmer_mean     ,
             klass.arithmetic_mean ) 

# =============================================================================

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
