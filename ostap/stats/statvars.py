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
    ##
    'data_decorate'        , ## technical function to decorate the classes
    'expression_types'     , ## valid types for expressions/cuts/weights
)
# =============================================================================
from   ostap.math.base           import ( isequal, iszero, axis_range, strings, 
                                          all_entries )      
from   ostap.core.core           import Ostap, rootException, WSE, VE, std     
from   ostap.core.ostap_types    import ( string_types , integer_types  , 
                                          num_types    , dictlike_types )
from   ostap.trees.cuts          import expression_types, vars_and_cuts
from   ostap.utils.basic         import loop_items, typename 
from   ostap.utils.progress_conf import progress_conf
import ostap.frames.frames       as     F 
import ostap.logger.table        as     T
import ostap.stats.moment 
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.statvars' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
FIRST   = Ostap.FirstEvent
LAST    = Ostap.LastEvent
# =============================================================================
LARGE   = 50000    ## allow frames or parallel for LARGE datasets 
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
#  result  = data.get_stat ( statobj , 'x+y' , 'pt>1' ) 
#  @encode
#  @see Ostap::Math::Statistic
#  @see Ostap::Math::WStatistic
#  @see Ostap::Math::Statistic2
#  @see Ostap::Math::WStatistic2
#  @see Ostap::Math::Statistic3
#  @see Ostap::Math::WStatistic3
#  @see Ostap::Math::Statistic4
#  @see Ostap::Math::WStatistic4
#  @see Ostap::statVar::the_moment
def data_get_stat ( data               ,
                    statobj            ,
                    expressions        ,
                    *args              , 
                    cuts       = ""    , 
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
        with rootException() :
            the_args = var_lst + ( cuts , cut_range ) + args 
            sc       = sv.get_stat ( data , statobj , *the_args )
        assert sc.isSuccess() , 'Error %s from StatVar::get_stat' % sc 
        return statobj

    if  use_frame and isinstance ( data , ROOT.TTee    ) and LARGE < len ( data ) :
        if all_entries ( data , *args[:2] ) and not arg[2:]  :
            return F.frame_project ( data                   ,
                                     model       = statobj  ,
                                     expressions = var_lst  ,
                                     cuts        = cuts     ,
                                     progress    = progress ,
                                     report      = progress ,
                                     lazy        = False    ) 
    
    if  parallel   and isinstance ( data , ROOT.TChain ) and LARGE < len ( data ) :
        if all_entries ( data , *args[:2] ) and not args [2:]  :
            import ostap.trees.trees
            if 1 < len ( data.files () ) :
                from ostap.parallel.parallel_stavar import parallel_get_stat
                return parallel_get_stat ( data        ,
                                           statobj     ,
                                           expressions ,
                                           cuts       = cuts          ,
                                           progress   = progress      ,
                                           use_frame  = False         , ## NB!!
                                           chunk_size = 2 * LARGE     ,
                                           maX_files  = 1             ,
                                           silent     = not progress  ) ;

    assert isinstance ( data , ROOT.TTree ) , "Invalid type for data!"
    
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
#  print data.the_moment (        3 , 0.0 , 'mass' , 'pt>1' ) ## ditto
#  @endcode
#  @see Ostap::Math::Moment_
#  @see Ostap::Math::WMoment_
def data_the_moment ( data               ,
                      order              ,
                      expression         ,
                      *args              , 
                      cuts       = ''    ,
                      as_weight  = True  , ## interpret cuts as weight 
                      cut_range  = ''    ,
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
        
    return data_get_stat ( data                  ,
                           statobj               ,
                           expression            , 
                           *args                 ,
                           cuts      = cuts      ,
                           cut_range = cut_range ,
                           progress  = progress  , 
                           use_frame = use_frame , 
                           parallel  = parallel  )

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
                   *args              , 
                   cuts       = ''    ,
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
                                   expression , *args     ,                                   
                                   cuts       = cut       ,
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
                *args              , 
                cuts       = ''    ,
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
                           *args                 ,
                           cuts      = cuts      ,
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
                     *args              , 
                     cuts       = ''    ,
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

    if use_frame and isinstance ( data , ROOT.TTree  ) and len ( data ) > LARGE :
        if all_entries ( data , *args[:2] ) and not args[2:] :
            return F.frame_statistic ( data        ,
                                       expressions ,
                                       cuts        = cuts      , 
                                       as_weight   = as_weight , 
                                       progess     = progress  ,
                                       report      = progress  ,
                                       lazy        = False     ) 

    if parallel and isinstance ( data , ROOT.TChain ) and len ( data ) > LARGE :
        if all_entries ( data , *args[:2] ) and not args [2:] : 
            import ostap.trees.trees   
            if 1 < data.nfiles () :
                from ostap.parallel.paralle_statvar import parallel_statistic
                return paralllel_statistic ( data        ,
                                             expressions ,
                                             cuts        = cuts       ,
                                             as_weight   = as_weight  , 
                                             progress    = progress   ,
                                             use_frame   = use+_frame )

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
        else                                       :
            return sv.statVar_cut ( data , varname , cuts ,             *args )
           
    ## 
    
    if   isinstance ( data , ROOT.RooAbsData ) : vcnt = Ostap.StatVar.WStatVector ()
    elif cuts and as_weight                    : vcnt = Ostap.StatVar.WStatVector ()
    else                                       : vcnt = Ostap.StatVar. StatVector ()
        
    ## variable names 
    vnames = strings ( var_lst ) 
    
    if isinstance ( data , ROOT.RooAbsData ) :
        with rootException() :
            sc   = sv.statVars ( data , vcnt , vnames  , cuts , cut_range , *args ) 
            assert sc.isSuccess() , 'Error %s from Ostap::StatVar::statVars' % sc  
            assert len ( vcnt ) == len ( vnames ) , "Mismatch in structure"
            return { name : cnt for ( name , cnt ) in zip ( var_lst , vct ) } 
            
    assert isinstance ( data , ROOT.TTree ) , "Invalid type for data!"
    
    from ostap.trees.trees import ActiveBranches
    with rootException() , ActiveBranches ( data , cuts , *var_lst ) :
        sc   = sv.statVars ( data , vcnt , vnames , cuts , *args ) 
        assert sc.isSuccess() , 'Error %s from Ostap::StatVar::statVars' % sc 
        assert len ( vcnt ) == len ( vnames ) , "Mismatch in structure"
        return { name : cnt for ( name , cnt ) in zip ( var_lst , vcnt ) } 

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
                  *args              , 
                  cuts       = ''    ,
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
                               expressions , *args    ,
                               cuts      = cuts       ,
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
                 *args              , 
                 cuts       = ''    ,
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
    results = data_minmax ( data,
                            expressions , *args ,
                            cuts      = cuts       ,
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
            mn , mx = data_minmax ( data , expressions    , 
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
                      expressions , *args ,                       
                      cuts        = ''    ,
                      cut_range   = ''    ,
                      progress    = False , 
                      as_weight   = True  ,  
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
        with rootException() :
            if 2 == N : sc = sv.statCov ( data , result , var_lst[0] , var_lst[1] , cuts , cut_range , *args )
            else      : sc = sv.statCov ( data , result , vnames                  , cuts , cut_range , *args )
            assert sc.isSuccess() , 'Error %s from StatVar::statVars' % sc 
            return result 

    if  use_frame and isinstance ( data , ROOT.TTee   ) and all_entries ( data , *args[:2] ) :
        raise NotImplementedError ( 'Not yet implemented yet!' )
    
    if  parallel  and isinstance ( data , ROOT.TChain ) and all_entries ( data , *args[:2] ) : 
        import ostap.trees.trees
        if 1 < len ( data.files () ) : raise NotImplementedError ( 'Not yet implemented yet!' )
        
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
                      expressions , *args ,                       
                      cuts        = ''    ,
                      cut_range   = ''    ,
                      progress    = False , 
                      as_weight   = True  ,  
                      use_frame   = False ,
                      parallel    = False ) :
    """ Get the effective vector of mean-values with covariances for the dataset
    >>> vct = data_statvct ( data , 'a,y,z' ,'w' ) 
    """
    
    ## decode expressions & cuts
    var_lst , cuts, input_string = vars_and_cuts ( expressions , cuts )
    N = len ( var_lst )
    
    assert 2 <= N , "At least two variables are needed!"

    covs = data_covariance ( data , var_lst , *args ,
                             cuts       = cuts      ,
                             cut_range  = cut_range ,                             
                             progress   = progress  , 
                             as_weight  = as_weight ,  
                             use_frame  = use_frame ,
                             parallel   = parallel  )

    ## HERE !!!
    ## FIX IT!!!
    
    return ## Ostap.VectorE ( N ) ( v , cov2 ) 

# ==============================================================================
## Get the (weighted) sum over the variable(s)
#  @code
#  data    = ...
#  result  = data_sum ( data , 'x+y'   , 'pt>1' ) 
#  results = data_sum ( data , 'x;y;z' , 'pt>1' ) ## result is dictionary 
#  @encode
#  @see Ostap::StatVar::statVar
def data_sum ( data        ,
               expressions        ,
               *args              , 
               cuts       = ''    ,
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
                              expressions , *args   ,
                              cuts      = cuts      ,
                              cut_range = cut_range ,
                              progress  = progress  ,
                              as_weight = as_weight ,
                              use_frame = use_frame ,
                              parallel  = parallel  )
    result2 = {} 
    if isinstance ( result , dictlike_types ) :
        for key , r in loop_items ( result ) : 
            result2 [ key ] = VE ( r.sum() , r.sum2() )
    else :
        result2 = VE ( result.sum() , result.sum2() )
    
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
                expression         ,
                *args              , 
                cuts       = ''    ,
                cut_range  = ''    ,
                progress   = False , 
                as_weight  = True  , ## interpret cuts as weiggt 
                use_frame  = False ,
                parallel   = False ) :         
    """ Get statistic from data 
    >>> data    = ...
    >>> result  = data_nEff ( data , 'x/y+z' )
    - see Ostap.StatVar.nEff
    """
    ## decode expressions & cuts
    var_lst, cuts, input_string = vars_and_cuts ( expressions , cuts )    
    assert 1 == len ( var_lst ) and input_string , "Inbvalid expressions: %s" % expression 
    expression = var_list[0]

    stat = data_statistic ( data ,
                            expression , *args     ,
                            cuts       = cuts      ,
                            cut_range  = cut_range ,
                            progress   = progress  ,
                            as_weight  = as_Weight ,
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
                         *args              , 
                         cuts       = ''    ,
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
    var_lst, cuts, input_string = vars_and_cuts ( expressions , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"
    
    ## 
    if isinstance ( ROOT.RooAbsData ) or ( cuts and as_weight ) :
        stat = Ostap.StatVar.WHarmonicMean ()
    else : 
        stat = Ostap.StatVar. HarmonicMean ()
        
    return data_get_stat  ( data , stat , var_lst , *args , 
                            cuts       = cuts      ,
                            cut_range  = cut_range ,
                            progress   = progress  ,
                            use_frame  = use_frame ,
                            parallel   = parallel  ) 

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
                          *args              , 
                          cuts       = ''    ,
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
    var_lst, cuts, input_string = vars_and_cuts ( expressions , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"
    
    ## 
    if isinstance ( data , ROOT.RooAbsData ) or ( cuts and as_weight ) :
        stat = Ostap.StatVar.WHarmonicMean ()
    else : 
        stat = Ostap.StatVar. HarmonicMean ()
        
    return fata_get_stat  ( data , stat , var_lst , *args , 
                            cuts       = cuts      ,
                            cut_range  = cut_range ,
                            progress   = progress  ,
                            use_frame  = use_frame ,
                            parallel   = parallel  ) 

    
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
                      *args              , 
                      cuts       = ''    ,
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
    var_lst, cuts, input_string = vars_and_cuts ( expressions , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"

    if    p == -1 or isequal ( p , -1. ) :
        return data_harmonic_mean   ( data , expression , *args ,
                                      cuts       = cuts      ,
                                      cut_range  = cut_range ,
                                      progress   = progress  ,
                                      use_frame  = use_frame ,
                                      parallel   = parallel  ) 
    elif  p ==  0 or iszero  ( p       ) :
        return data_geometric_mean  ( data , expression , *args , 
                                      cuts       = cuts      ,
                                      cut_range  = cut_range ,
                                      progress   = progress  ,
                                      use_frame  = use_frame ,
                                      parallel   = parallel  ) 
    elif  p ==  1 or isequal ( p ,  1. ) :
        return data_arithmetic_mean ( data , expression , *args ,
                                      cuts       = cuts      ,
                                      cut_range  = cut_range ,
                                      progress   = progress  ,
                                      use_frame  = use_frame ,
                                      parallel   = parallel  ) 

    
    ## 
    if isinstance ( data , ROOT.RooAbsData ) or ( cuts and as_weight ) :
        stat = Ostap.StatVar.WHarmonicMean ()
    else : 
        stat = Ostap.StatVar. HarmonicMean ()
        
    return data_get_stat  ( data , stat , var_lst , *args , 
                            cuts       = cuts      ,
                            cut_range  = cut_range ,
                            progress   = progress  ,
                            use_frame  = use_frame ,
                            parallel   = parallel  ) 

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
                           *args              , 
                           cuts       = ''    ,
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
    var_lst , cuts , _  = vars_and_cuts ( expressions , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"
    
    if isinstance ( ROOT.RooAbsData ) or ( cuts and as_weight ) :
        stat = Ostap.StatVar.WArithmeticMean ()
    else : 
        stat = Ostap.StatVar. ArithmeticMean ()
        
    return data_get_stat  ( data , stat , var_lst , *args , 
                            cuts       = cuts      ,
                            cut_range  = cut_range ,
                            progress   = progress  ,
                            use_frame  = use_frame ,
                            parallel   = parallel  ) 

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
                       *args              , 
                       cuts       = ''    ,
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
    var_lst , cuts , _  = vars_and_cuts ( expressions , cuts )
    assert 1 == len ( var_lst ) , "Invalid expression!"
    
    assert isinstance ( p , num_types ) , 'Invalid p-parameter!'
    
    if   p == 0 or iszero  ( p       ) :
        return data_harmonic_mean   ( data , expression , *args , 
                                      cuts       = cuts      ,
                                      cut_range  = cut_range ,
                                      progress   = progress  ,
                                      use_frame  = use_frame ,
                                      parallel   = parallel  ) 
    elif p == 1 or isequal ( p , 1.0 ) :
        return data_arithmetic_mean ( data , expression , *args ,
                                      cuts       = cuts      ,
                                      cut_range  = cut_range ,
                                      progress   = progress  ,
                                      use_frame  = use_frame ,
                                      parallel   = parallel  ) 
                                      

    if isinstance ( data , ROOT.RooAbsData ) or ( cuts and as_weight ) :
        stat = Ostap.StatVar.WArithmeticMean ()
    else : 
        stat = Ostap.StatVar. ArithmeticMean ()

    return data_get_stat  ( data , stat , var_lst , *args , 
                            cuts       = cuts      ,
                            cut_range  = cut_range ,
                            progress   = progress  ,
                            use_frame  = use_frame ,
                            parallel   = parallel  ) 
        

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
                *args              , 
                cuts       = ''    ,
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
    m2 = data_the_moment ( data , 2 , expression , *args ,
                           cuts       = cuts      ,
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
                    *args              , 
                    cuts       = ''    ,
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
    m4 = data_the_moment ( data , 4 , expression , *args ,
                           cuts       = cuts      ,
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
def data_rms ( data  , 
               expression         ,
               *args              , 
               cuts       = ''    ,
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
    
    variance  = data_variance ( data , expression , *args ,
                                cuts       = cuts      ,
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
                    *args              , 
                    cuts       = ''    ,
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
    m6 = data_the_moment ( data , 6 , expression , *args ,
                           cuts       = cuts      ,
                           cut_range  = cut_range ,
                           progress   = progress  ,
                           use_frame  = use_frame ,
                           parallel   = parallel  ) 
    return m6.skewness ()


# =============================================================================
## get the (excess) kurtosis (with uncertainty)
#  @code
#  data =  ...
#  print data_kustosis ( data , 'mass' , 'pt>1' ) 
#  print data.kustosis (        'mass' , 'pt>1' ) ## ditto
#  @endcode
#  @see Ostap::StatVar::kurtosis
def data_kurtosis ( data  , 
                    expression         ,
                    *args              , 
                    cuts       = ''    ,
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
    m8 = data_the_moment ( data , 8 , expression , *args ,
                           cuts       = cuts      ,
                           cut_range  = cut_range ,
                           progress   = progress  ,
                           use_frame  = use_frame ,
                           parallel   = parallel  ) 
    return m6.kurtosis ()

# =============================================================================
def data_project ( data        ,
                   target      ,
                   expressions ,
                   *args               , 
                   cuts        = ''    ,
                   cut_range   = ""    ,
                   as_weight   = ""    ,
                   progress    = True  ,
                   use_frame   = False ,
                   parallel    = False ) :
    
    ## (1) decode expressions & cuts
    var_lst , cuts, _  = vars_and_cuts ( expressions , cuts )
    nvars = len ( var_lst )
    
    ## (2) check consistency
    if   1 == nvars and isinstance ( target , _s1D     ) : pass
    elif 2 == nvars and isinstance ( target , _s2D     ) : pass
    elif 3 == nvars and isinstance ( target , _s3D     ) : pass
    elif 4 == nvars and isinstance ( target , _s4D     ) : pass
    elif 3 == nvars and isinstance ( target , ROOT.TH3 ) and 3 == target.dim() : pass    
    elif 2 == nvars and isinstance ( target , ROOT.TH2 ) and 2 == target.dim() : pass
    elif 1 == nvars and isinstance ( target , ROOT.TH1 ) and 1 == target.dim() : pass
    else :
        raise TypeError ( 'Inconsistent statobj %s & expression(s): %s' % \
                          ( typename ( statobj ) , str ( var_lst ) ) ) 

    ## (3) cut_range defined *only* for RooFit datasets 
    if cut_range and not isinstance ( data , ROOT.RooAbsData ) : 
        raise TypeError ( "Invalid use of `cut_range':%s" % cut_range  ) 
    
    ## (4) display progress ? 
    progress = progress_conf ( progress )

    ## (5) create the driver 
    pv = Ostap.Project ( progress )
    
    ## (6) RooFit ?
    if isinstance ( data , ROOT.RooAbsData ) :
        with rootException() :
            the_args = var_lst +  ( cuts , cut_range ) + args
            if   1 == nvars : sc = pv.project1 ( data , target , *the_args )
            elif 2 == nvars : sc = pv.project2 ( data , target , *the_args )
            elif 3 == nvars : sc = pv.project3 ( data , target , *the_args )            
            elif 4 == nvars : sc = pv.project4 ( data , target , *the_args )            
        assert sc.isSuccess() , 'Error %s from StatVar::project(1,2,3)' % sc 
        return target
    
    if  use_frame and isinstance ( data , ROOT.TTee    ) and LARGE < len ( data ) :
        if all_entries ( data , *args[:2] ) and not arg[2:]  :
            return F.frame_project ( data                   ,
                                     model       = target   ,
                                     expressions = var_lst  ,
                                     cuts        = cuts     ,
                                     progress    = progress ,
                                     report      = progress ,
                                     lazy        = False    ) 
    
    if  parallel   and isinstance ( data , ROOT.TChain ) and LARGE < len ( data ) :
        if all_entries ( data , *args[:2] ) and not args [2:]  :
            import ostap.trees.trees
            if 1 < len ( data.files () ) :
                pass
            
    ## Branches to be activated
    from ostap.trees.trees import ActiveBranches
    with rootException() , ActiveBranches ( data , cuts , *var_lst ) :
        the_args = var_lst + ( cuts , ) + args         
        if   1 == nvars : sc = pv.project1 ( data , target , *the_args  )
        elif 2 == nvars : sc = pv.project2 ( data , target , *the_args  )
        elif 3 == nvars : sc = pv.project3 ( data , target , *the_args  )
        elif 4 == nvars : sc = pv.project4 ( data , target , *the_args  )
        assert sc.isSuccess() , 'Error %s from StatVar::project(1,2,3,4)' % sc 
        return statobj
    
# =============================================================================
## decorate certain class with some useful  methods 
def data_decorate ( klass ) :
    """ Decorate certain class with some useful  methods
    """
    
    if hasattr ( klass , 'get_moment'     ) : klass.orig_get_moment     = klass.get_moment
    if hasattr ( klass , 'the_moment'     ) : klass.orig_the_moment     = klass.the_moment
    
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
    if hasattr ( klass , 'statCov'        ) : klass.orig_statCov        = klass.statCov 
    if hasattr ( klass , 'statCovs'       ) : klass.orig_statCovs       = klass.statCovs
    
    klass.get_moment      = data_the_moment
    klass.the_moment      = data_the_moment
    
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

## Quantile  = StatVar.Quantile
## Quantiles = StatVar.Quantiles
## Interval  = StatVar.Interval 
## QInterval = StatVar.QInterval 

## def _q_str_  ( o ) : return "Quantile(%.5g,n=%s)"         % ( o.quantile  , o.nevents )
## def _qs_str_ ( o ) : return "Quantiles(%s,n=%s)"          % ( o.quantiles , o.nevents )
## def _i_str_  ( o ) : return "Interval([%.5g,%.5g])"       % ( o.low       , o.high    )
## def _qi_str_ ( o ) : return "QInterval([%.5g,%.5g],n=%s)" % ( o.interval.low  ,
##                                                              o.interval.high , o.nevents )

## Quantile  .__str__  = _q_str_
## Quantile  .__repr__ = _q_str_
## Quantiles .__str__  = _qs_str_
## Quantiles .__repr__ = _qs_str_
## Interval  .__str__  = _i_str_
## Interval  .__repr__ = _i_str_
## QInterval .__str__  = _qi_str_
## QInterval .__repr__ = _qi_str_
    
# =============================================================================


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
