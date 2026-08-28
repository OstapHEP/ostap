#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof_utils.py
#  Set of utilities for goodness-of-fit studies 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2023-12-06
# =============================================================================
""" Simple utilities for Goodness-of-Fit studies 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2023-12-06"
__all__     = (
    'mean_var'           , ## mean and variance for (weighted) arrays
    'nEff'               , ## get number of effective entries
    'normalize'          , ## "normalize" variables in dataset/structured array
    'normalize_pooled'   , ## "normalize" variables in pooled dataset 
    'clip_pvalue'        , ## clip-value
    'pairwise_distances' , ## get array of all pair-wise distances between two datasets
    'nearest_neighbors'  , ## get all nearest neigbours
    'nearest_distances'  , ## get all nearest distances 
    's2u'                , ## convert structured numpy array into non-structured
    'combine_pvalues'    , ## combine p-values 
    'np2vct'             , ## numpy arrays into SVectorWithError
    ##
)
# =============================================================================
from   ostap.core.core          import Ostap, VE, SE 
from   ostap.utils.cidict       import cidict, cidict_fun
from   ostap.core.ostap_types   import string_types, num_types, numpy_floats  
from   ostap.math.math_base     import doubles, axis_range
from   ostap.math.math_ve       import significance
from   ostap.math.ve            import fmt_pretty_ve
from   ostap.stats.utils        import weight_trivial
from   ostap.stats.counters     import SE, WSE, ECDF
from   ostap.utils.basic        import ( numcpu   , num_jobs     , 
                                         typename , run_parallel ) 
from   ostap.utils.progress_bar import progress_bar
from   ostap.logger.symbols     import ( plus_minus  , greek_lower_sigma ,
                                         subscript_A , subscript_K , subscript_C , 
                                         likelihood  , script_t    , script_p    ,
                                         sub_min     , sub_max     ,
                                         sub_mean    , sub_rms     , 
                                         infinity_pos as pos_infinity_symbol     )                                         
from   ostap.logger.pretty      import pretty_float, pretty_row 
from   ostap.plotting.color     import Orange, Green, Blue
from   packaging.version        import Version
import ostap.math.linalg 
import ostap.logger.table       as     T 
import ROOT, os, sys, math, numpy, scipy, abc    
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gof_utils' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-fit studies' )#
# =============================================================================
## transform structured array to unstructured one
s2u = None 
# =============================================================================
try : # =======================================================================
    # =========================================================================  
    from numpy.lib.recfunctions import structured_to_unstructured as s2u
    # =========================================================================
except ImportError :
    # =========================================================================
    s2u = None 
# =============================================================================
## pair-wise distances between two datasets 
pairwise_distance = None # ====================================================
# =============================================================================
## (1) use cdist from scipy.spatial
# =============================================================================
try : # =======================================================================
    # =========================================================================
    from scipy.spatial.distance import cdist as _scipy_pw_distances
    ## Calculate all pair-wise distances using scipy 
    def pairwise_distances ( data1  ,
                             data2  ,
                             metric = 'sqeuclidean'   , *  ,
                             n_jobs = None            , **kwargs ) : 
        """ Calculate all pair-wise distances using scipy 
        -   scipy.spatial.distance.cdist  
        """
        return _scipy_pw_distances ( data1 , data2 , metric , **kwargs ).flatten()
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    pass
# =============================================================================
## (2) make a try to use more-efficient (?) pair-wise distances from scikit-learn
# =============================================================================
has_sklearn = True # ==========================================================
# =============================================================================
try : # =======================================================================
    # =========================================================================
    from sklearn.metrics import pairwise_distances as _sk_pw_distances    
    ## Calculate all pair-wise distances   using sklearn/scikit-learn 
    def pairwise_distances ( data1  ,
                             data2  ,
                             metric = 'sqeuclidean' , **kwargs ) : 
        """ Calculate all pair-wise distances using scikit-learn
        -  see `sklearn.metrics.pairwise_distances` 
        """
        if data1 is data2 : data2 = None 
        return _sk_pw_distances ( data1 , Y = data2 , metric = metric , **kwargs ).flatten()
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    has_sklearn = False # =====================================================
    # =========================================================================

# =============================================================================
## Nearest distances 
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import sklearn.neighbors 
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    has_sklearn = False 

# =============================================================================
## (2) Helper constants
# =============================================================================
WORKER_KEY    = 'workers' if Version ( "1.6.0" ) <= Version ( scipy.__version__ ) else 'n_jobs'
SCIPY_METRICS = {
    'minkowski' , 'euclidean' , 'l2' , 'sqeuclidean' , 
    'manhattan' , 'cityblock' , 'l1' ,
    'chebyshev' , 'infinity'
}
# =============================================================================
## (3) Nearest Neighbors implementation
# =============================================================================
def nearest_neighbors ( data ,
                        n_neighbors = 10 ,
                        metric      = 'euclidean' ,
                        p           = 2           ,
                        n_jobs      = -1           , **config ) : 
    """ Get nearest neighbors using scipy.spatial.cKDTree (fast L_p metrics) 
        or sklearn.neighbors (advanced metrics)
    """
    metric_clean = str ( metric ).lower ()
    
    # =========================================================================
    ## Use scipy.spatial.cKDTree for standard L_p metrics
    # =========================================================================
    if metric_clean in SCIPY_METRICS or isinstance ( metric , ( int , float ) ) :
        
        if   metric_clean in ( 'euclidean' , 'sqeuclidean' , 'l2' ) : p_val = 2
        elif metric_clean in ( 'manhattan' , 'cityblock'   , 'l1' ) : p_val = 1
        elif metric_clean in ( 'chebyshev' , 'infinity'           ) : p_val = numpy.inf
        else : p_val = p

        config [ WORKER_KEY ] = n_jobs
        tree = scipy.spatial.cKDTree ( data )
        _ , indices = tree.query ( data , k = n_neighbors + 1 , p = p_val , **config )
        return indices [ :, 1: ]

    # =========================================================================
    ## Use sklearn.neighbors for advanced metrics (cosine, haversine, etc.)
    # =========================================================================
    elif has_sklearn :
        
        nn = sklearn.neighbors.NearestNeighbors ( n_neighbors = n_neighbors + 1 , 
                                                  metric      = metric          , 
                                                  p           = p               , 
                                                  n_jobs      = n_jobs          , **config )
        nn.fit ( data )
        _ , indices = nn.kneighbors ( data )
        return indices [ :, 1: ]
    
    # =========================================================================
    ## Error handling
    # =========================================================================
    else : # ==================================================================
        # =====================================================================
        raise ValueError ( f"Metric '{metric}' requires sklearn, which is not installed." )
        # =====================================================================
        

# =============================================================================
## Compute distances to the k-th nearest neighbor for each point in data.
# 
#  Parameters:
#    data (array-like): Input dataset of shape (N, D).
#    k (int): Target neighbor index (default is 1 for the immediate nearest neighbor).
#       metric (str or int): Distance metric to use.
#       p (int): Power parameter for the Minkowski metric.
#         n_jobs (int): Number of parallel jobs for neighbor search (-1 uses all cores).
#         **config: Additional parameters passed to the tree search backend.
#        
#    Returns:
#        numpy.ndarray: 1D array of shape (N,) containing distances to the k-th neighbor.
def nearest_distances (  data   ,
                         k      =  1            ,
                         metric = 'euclidean'   ,
                         p      =  2            ,
                         n_jobs = -1 , **config ):
    """ Compute distances to the k-th nearest neighbor for each point in data.
    
    Parameters:
        data (array-like): Input dataset of shape (N, D).
        k (int): Target neighbor index (default is 1 for the immediate nearest neighbor).
        metric (str or int): Distance metric to use.
        p (int): Power parameter for the Minkowski metric.
        n_jobs (int): Number of parallel jobs for neighbor search (-1 uses all cores).
        **config: Additional parameters passed to the tree search backend.
        
    Returns:
        numpy.ndarray: 1D array of shape (N,) containing distances to the k-th neighbor.
    """
    metric_clean = str(metric).lower()
    
    # =========================================================================
    ## (1) Use scipy.spatial.cKDTree for fast standard L_p metrics
    # =========================================================================
    if metric_clean in SCIPY_METRICS or isinstance(metric, (int, float)):
        
        if   metric_clean in ( 'euclidean' , 'sqeuclidean' , 'l2' ) : p_val = 2
        elif metric_clean in ( 'manhattan' , 'cityblock'   , 'l1' ) : p_val = 1
        elif metric_clean in ( 'chebyshev' , 'infinity'           ) : p_val = numpy.inf
        else                                                        : p_val = p

        config [ WORKER_KEY ] = n_jobs
        tree = scipy.spatial.cKDTree( data )
        
        # Query k + 1 neighbors because the 0-th neighbor is the point itself (distance = 0.0)
        distances, _ = tree.query(data, k=k + 1, p=p_val, **config)
        
        # Extract the k-th neighbor column to guarantee a flat 1D array of shape (N,)
        return distances[:, k]

    # =========================================================================
    ## (2) Use sklearn.neighbors for advanced metrics (cosine, haversine, etc.)
    # =========================================================================
    elif has_sklearn:
        
        nn = sklearn.neighbors.NearestNeighbors ( n_neighbors = k + 1  ,
                                                  metric      = metric ,
                                                  p           = p      ,
                                                  n_jobs      = n_jobs , **config ) 
        nn.fit ( data )
        
        distances , _ = nn.kneighbors(data)
        
        # Extract the k-th neighbor column to guarantee a flat 1D array of shape (N,)
        return distances[:, k]
    
    # =========================================================================
    ## Error handling
    # =========================================================================
    else: # ===================================================================
        # =====================================================================
        raise ValueError(f"Metric '{metric}' requires scikit-learn, which is not installed.")
    
# =============================================================================
## Get the mean and variance for (1D) data array with optional (1D) weight array
#  @code
#  ds = ... ## dataset as structured array
#  mean, cov2 = mean_var ( ds ['x'] )
#  @endcode
#  - with weight 
#  @code
#  ds = ... ## dataset as structured array with weight 
#  mean, cov2 = mean_var ( ds ['x'] , ds['weight'] )
#  @endcode 
def mean_var ( data , weight = None ) :
    """ Get the mean and variance for 1D-data array with optional 1D-weight array

    >>> ds = ... ## dataset as structured array
    >>> mean, cov2 = mean_var( ds ['x'] )
    
    - with weight 
    
    >>> ds = ... ## dataset as structured array with weight 
    >>> mean, cov2 = mean_var ( ds ['x'] , ds['weight'] )
    """
    #
    if weight_trivial ( weight ) :
        mean = numpy.mean ( data , axis = 0 , dtype = float ) 
        var  = numpy.var  ( data , axis = 0 , dtype = float ) 
        return mean , var
    # 
    mean  = numpy.average (   data               , weights = weight , axis = 0 )
    var   = numpy.average ( ( data - mean ) ** 2 , weights = weight , axis = 0 )
    #
    return mean , var 

# =============================================================================
## Get the effective number of entries for 1D-array
#  \f{ n_{eff} = \frac{  \left\langle x \right\rangle^2}
#                     { \left\langle x^2 \right\rangle } \f}
def nEff ( weights ) :
    """ Get the effective number of entries for 1D-array
    n_eff = ( sum ( x )  ) ^2 / sum ( x^2 )
    """
    s1 = numpy.sum ( weights      , dtype = float )
    s2 = numpy.sum ( weights ** 2 , dtype = float )
    return s1 * s1 / s2

# =============================================================================
## Normalize/standartize several numpy arrays such that mean and rms for the pooled sample
#  are equal to 0 and 1 correspondinly
#  @code
#  ds1 = ...
#  ds2 = ...
#  ds3 = ...
#  ds1 , ds2 , ds3 = normalize_pooled ( ds1 , ds2 , ds3 ) 
#  @endcode 
def normalize_pooled ( *datasets ) :
    """ Normalize/standartize several numpy arrays such that the mean and rms for the POOLED sample
    are equal to 0 and 1 correspondinly
    
    >>> ds1 = ...
    >>> ds2 = ...
    >>> ds3 = ...
    >>> ds1 , ds2 , ds3 = normalize_pooled ( ds1 , ds2 , ds3 )
    
    """

    if not datasets : return ()

    total      = 0 
    total_mean = None
    total_std2 = None
    
    sizes      = []
    means      = []
    
    for i , data in enumerate ( datasets ) :
        
        n     = len ( data )
        
        mean = numpy.mean ( data , axis=0 , keepdims = True )
        std2 = numpy.std  ( data , axis=0 , keepdims = True ) ** 2 
        
        sizes .append ( n    ) 
        means .append ( mean )
        
        if total_mean is None : total_mean  = n * mean
        else                  : total_mean += n * mean
        
        if total_std2 is None : total_std2  = n * std2
        else                  : total_std2 += n * std2
        
    ## total number of entries in pooled sample
    total       = sum ( sizes )

    ## mean values 
    total_mean /= total

    for n , mu in zip ( sizes , means ) : total_std2 += n * ( mu - total_mean ) ** 2
    
    ## squared RMS
    total_std2 /= total

    ## mean value 
    mean = total_mean

    ## RMS 
    std  = numpy.sqrt ( total_std2 )
    
    result = []
    for i , data in enumerate ( datasets ) :
        ### shift&scale 
        result.append ( ( data - mean ) / std )

    if 1 == len ( result ) : return result [ 0 ]
    
    return tuple ( result ) 

# =============================================================================
## Get the "normalized" input datasets
#  All floating fields  are calculated as
#  \f[ x = \frac{x - \left\langle x \right\rangle}{\sigma} \f]
#  where \f$ \left\langle x \right\rangle\f$ is a mean  value
#  and \f$ \sigma \f$ is a standard deviation.
# 
#  @code
#  ds      = ... # data set as structured array
#  dsn, J  = normalize ( ds ) 
#  @endcode
#
#  - If several datasets are specified, all floating names must be the same
#  and the mean and sigma are either taken either from the first dataset,
#  if <code>first=True</code> or as combined through all datasets otherwise 
#
#  @code
#  ds1 = ... # data set as structured array
#  ds2 = ... # data set as structured array
#  ds3 = ... # data set as structured array
#  ds1n, ds2n, ds3n = normalize ( ds1 , ds2 , ds3 , first = True )
#  @endcode
#
#  - If <code>weight</code> is specified, this floating column is considered
#  as the weight
#
#  @code
#  ds  = ... # data set as structured array with weight 
#  dsn = normalize ( ds , weight = 'weight' ) 
#  @endcode
#
#  @code
#  ds1 = ... # data set as structured array without weight 
#  ds2 = ... # data set as structured array with weight 
#  datasets    = normalize ( ds1 , ds2 , weight = ( None , 'weight'  ) )
#  ds1n , ds2n = datasets 
#  @endcode
#
#  @attention Only the floating point columns are transformed! 
#  @attention Input datasets are expected to be numpy structured arrays
#
#  @code
#  ds = ... # data set as structured array
#  dsn = normalize ( ds ) 
#  @endcode
def normalize ( ds , *others , weight = () , first = True ) :
    """ Get the `normalized' input datasets
    All floating fields  are calculated as
    
    x = (x - <x>)/sigma
    
    - <x> is a mean value 
    - is a standard deviation.
    
    - If several datasets are specified, all floating names must be the same
    and the mean and sigma are either taken either from the first dataset,
    if `first=True` or as combined through all datasets, otherwise  
    
    - If `weight` is specified, this floating column is considered
    as the weight
    
    - attention Only the floating point columns are transformed! 
    - attention Input datasets are expected to be numpy structured arrays 
    """

    datasets = ( ds , ) + others
    
    nd = len ( datasets ) 
    if not weight                             : weight = nd * [ ''     ]
    elif isinstance ( weight , string_types ) : weight = nd * [ weight ]
    
    assert ( len ( weight ) == nd ) and \
           all ( ( not w ) or isinstance ( w , string_types ) for w in weight ) , \
           'Invalid specification of weight!'
        
    weight = list ( weight )
    for i , w in enumerate ( weight ) :
        if not w : weight [ i ] = '' 
    weight = tuple ( weight )

    ds     = datasets [ 0  ]
    others = datasets [ 1: ]
    
    ## collect the floating columns 
    columns = []
    w0      = weight [ 0 ] 
    for n,t in ds.dtype.fields.items () :
        if t [ 0 ] in numpy_floats  and n != w0 : columns.append ( n ) 
        
    vmeans  = [] 
    for i , c in enumerate ( columns ) :
        mean, var    = mean_var ( ds [ c ] , None if not w0 else ds [ w0 ] )
        vmeans.append ( VE ( mean , var ) )
        
    ## Number of events/effective entries 
    nevents = 1.0 * ds.shape [ 0 ] if not w0 else nEff ( ds [ w0 ] )
    
    if not first and others : 
        nevents = ds.shape[0] 
        for k , dd in enumerate ( others ) :
            
            wk = weight [ k + 1 ] 
            nn = 1.0 * dd.shape [ 0 ] if not wk else nEff ( dd [ wk ] )                
            
            for i , c in enumerate ( columns ) :
                
                mean, var = mean_var ( dd [ c ] , None if not wk else dd [ wk ] )                    
                vv = VE ( mean , var ) 
                vmeans [ i ] = Ostap.Math.two_samples ( vmeans [ i ] , nevents , vv , nn )
                
            nevents += nn
                
    result = []  
    for d in datasets :
        
        nds = d.copy ()
        for ic , c in enumerate ( columns ) :
            vv        = vmeans [ ic ]
            mean      = vv.value ()
            sigma     = vv.error ()                 
            a         = nds [ c ]
            nds [ c ] =  ( a - mean ) / sigma
            
        result.append ( nds )
            
    return tuple ( result ) 

# =============================================================================
## Normalize sample weights to match effective dataset sizes sum(|w|) = N.
def normalize_weight_to_N ( weight , N ) :
    """ Normalize sample weights using sum of absolute values 
    such that sum(|w1|) == N1 and sum(|w2|) == N2.
    """
    weight = numpy.ones ( N , dtype = numpy.float32 ) if weight_trivial ( weight ) else numpy.asarray ( weight , dtype = numpy.float32 )
    weight *= N / numpy.sum ( numpy.abs ( weight ) , dtype = numpy.float64 )
    ## 
    return weight 

# =============================================================================
## Normalize sample weights to match effective dataset sizes sum(|w|) = N.
def normalize_weights_to_N ( weight1 ,
                             weight2 ,
                             N1      ,
                             N2      ) :
    """ Normalize sample weights using sum of absolute values 
    such that sum(|w1|) == N1 and sum(|w2|) == N2.
    """
    w1_trivial = weight_trivial ( weight1 )
    w2_trivial = weight_trivial ( weight2 )

    if w1_trivial and w2_trivial :
        return None , None

    weight1 = normalize_weight_to_N ( weight1 , N1 )
    weight2 = normalize_weight_to_N ( weight2 , N2 )
    
    return weight1 , weight2 

# =============================================================================
## Short labels for various statitical estimators 
Labels = {
    'KS'  : 'Kolmogorov-Smirnov'   ,
    'K'   : 'Kuiper'               ,
    'AD'  : 'Anderson-Darling'     ,
    'CM'  : 'Cramer-von Mises'     ,
    'ZK'  : "Zhang Z%s" % subscript_K ,
    'ZA'  : "Zhang Z%s" % subscript_A ,
    'ZC'  : "Zhang Z%s" % subscript_C ,        
    'BJ'  : 'Berk-Jones'           ,        
    'NLL' : '-log%s' % likelihood  ,
    'AIC' : 'Akaike IC'            ,
    'BIC' : 'Bayesian IC'          ,    
}
# =============================================================================
## lower-case shortcuts:
Keys = {    
    'KS'  : ( 'ks'  , 'kolmogorov' , 'kolmogorovsmirnov' ) , 
    'K'   : ( 'k'   , 'kuiper'                           ) , 
    'AD'  : ( 'ad'  , 'anderson'   , 'andersondarling'   ) , 
    'CM'  : ( 'cm'  , 'cramer'     , 'cramervonmises'    ) , 
    'ZK'  : ( 'zk'  , 'zhangk'     , 'zhangzk'           ) , 
    'ZA'  : ( 'za'  , 'zhanga'     , 'zhangza'           ) , 
    'ZC'  : ( 'zc'  , 'zhangc'     , 'zhangzc'           ) , 
    'BJ'  : ( 'bj'  , 'berkjones'  , 'berk'              ) , 
    'NLL' : ( 'nll' ,              ) , 
    'AIC' : ( 'aic' , 'akaike'     ) ,  
    'BIC' : ( 'bic' , 'bayesian'   ) , 
    }
# =============================================================================
assert Labels.keys() == Keys.keys() , "Mismatch between Labels & Keys structures!"
# =============================================================================
__all_1D_names = cidict ( transform = cidict_fun )
for k,vv in Keys.items() :
    for v in vv : __all_1D_names [ v ] = k
# =============================================================================
## Get the internal  name for extrnal 1D-method
def method_1D ( mm ) :
    """ Get the internal  name for external 1D-method
    """
    return __all_1D_names.get ( mm , '' ) 
# =============================================================================        
## clip p-value 
def clip_pvalue ( pvalue , clip = 0.5 ) :
    """ Clip p-value
    """
    pv = VE ( pvalue )
    ## everything is fine, no need to clip 
    if 0 < pv.value() < 1 : return pv  
    ##
    clip = min ( 0.5 , abs ( clip ) * pv.error () )
    ## 
    if   1 <= pv.value() : pv = VE ( 1 - clip , pv.cov2() )
    elif 0 >= pv.value() : pv = VE (     clip , pv.cov2() )
    ## 
    return pv 
    
# =============================================================================
pvalue_types = num_types + ( VE , ) 
# =============================================================================
col_t_value = '%s-value'      %   script_t
col_t_mean  = '%s%s'          % ( script_t , sub_mean ) 
col_t_rms   = '%s%s'          % ( script_t , sub_rms  ) 
col_t_min   = '%s%s'          % ( script_t , sub_min  )
col_t_max   = '%s%s'          % ( script_t , sub_max  )
col_t_unit  = '%s-unit'       %   script_t 
col_p_value = '%s-value [%%]' %   script_p
col_n_sigma = '#%s'           %   greek_lower_sigma
the_header  = ( col_t_value ,
                col_t_mean  ,
                col_t_rms   ,
                col_t_min   ,
                col_t_max   ,
                col_t_unit  ,
                col_p_value ,
                col_n_sigma )
# =============================================================================
## Format the row for 2S/GoF tables
#  @code
#  tvalue = ...
#  pvalue = ...
#  ecdf   = ... 
#  header , row = format_row ( tvalue = tvalue, pvalue = pvalue , ecdf = ecdf )
#  @endcode
def format_row ( tvalue    = None ,
                 pvalue    = None ,
                 ecdf      = None ,
                 counter   = None , 
                 precision = 4    ,
                 width     = 6    ) :
    """ Format the row for the 2S/GoF tables
    >>> tvalue = ...
    >>> pvalue = ...
    >>> ecdf   = ...
    >>> header , row = format_row ( tvalue = tvalue , pvalue = pvalue  , ecdf = ecdf )
    """
    
    has_tvalue  = not tvalue  is None and isinstance ( tvalue  , num_types     ) 
    has_pvalue  = not pvalue  is None and isinstance ( pvalue  , pvalue_types  ) 
    has_ecdf    = not ecdf    is None and isinstance ( ecdf    , ECDF          ) 
    has_counter = not counter is None and isinstance ( counter , ( SE , WSE )  )

    if has_ecdf  and not has_counter  :
        counter     = ecdf.counter ()
        has_counter = True 
        
    if has_tvalue and has_pvalue and has_counter :
        
        header = the_header

        mean       = counter.mean   ()
        rms        = counter.rms    () 
        vmin, vmax = counter.minmax () 

        items      = tvalue, mean, rms , vmin , vmax
        row, unit  = pretty_row  ( *items , width = width , precision = precision )
        unit       = '[%s]' % unit if unit else unit
        
        pv         = clip_pvalue  ( pvalue ) 
        nsigma     = significance ( pv     ) ## convert  it to significance
        ## 
        if isinstance ( nsigma , VE ) and nsigma.cov2 () <= 0 : nsigma = float ( nsigma ) 

        pvalue  = pvalue * 100

        if isinstance ( pvalue , VE )   : pvalue  = '%6.2f %s %.2f' % ( pvalue.value() , plus_minus , pvalue.error () )
        else                            : pvalue  = '%6.2f'         %   float ( pvalue )
        
        if   99 < float ( nsigma      ) : nsigma = pos_infinity_symbol 
        elif isinstance ( nsigma , VE ) : nsigma = '%.2f %s %.2f'   % ( nsigma.value() , plus_minus , nsigma.error () )
        else                            : nsigma = '%.2f'           %   float ( nsigma ) 

        row = row + ( unit , pvalue , nsigma )
        
        return header , row

    elif has_tvalue and has_pvalue :
        
        header    = ( col_t_value , 
                      col_t_unit  ,
                      col_p_value ,
                      col_n_sigma ) 

        items      = tvalue, 
        row, unit  = pretty_row  ( *items , width = width , precision = precision )
        unit       = '[%s]' % unit if unit else unit
        
        pv         = clip_pvalue  ( pvalue )
        nsigma     = significance ( pv     )
        if isinstance ( nsigma , VE ) and nsigma.cov2 () <= 0 : nsigma = float ( nsigma ) 
        
        pvalue  = pvalue * 100 
        if isinstance ( pvalue , VE )   : pvalue  = '%6.2f %s %.2f' % ( pvalue.value() , plus_minus , pvalue.error () )
        else                            : pvalue  = '%6.2f'         % float ( pvalue )

        if   99 < float ( nsigma      ) : nsigma = pos_infinity_symbol 
        elif isinstance ( nsigma , VE ) : nsigma  = '%.2f %s %.2f'  % ( nsigma.value() , plus_minus , nsigma.error () )        
        else                            : nsigma  = '%.2f'          % float ( nsigma ) 
        
        row     = row + ( unit , pvalue , nsigma )
        return header , row 

    elif has_tvalue :
        
        header     = col_t_value , col_t_unit 
        items      = tvalue, 
        row, unit  = pretty_row  ( *items , width = width , precision = precision )
        unit       = '[%s]' % unit if unit else unit
        row        = row + ( unit , ) 
        return header, row

    ## no data for table 
    return () , ()

# ==========================================================================
## Get results in form of the table 
def format_table ( tvalue    = None ,
                   pvalue    = None ,
                   ecdf      = None ,
                   counter   = None ,
                   precision = 4    ,
                   width     = 6    , 
                   title     = ''   ,
                   prefix    = ''   ,
                   style     = ''   ) :
    """ Get results in form of the table 
    """
    
    header , row = format_row ( tvalue    = tvalue    ,
                                pvalue    = pvalue    ,
                                ecdf      = ecdf      ,
                                counter   = counter   , 
                                precision = precision ,
                                width     = width     ) 
    rows = [ header ]
    rows.append ( row )
    rows = T.remove_empty_columns ( rows ) 
    return T.table ( rows               ,
                     title     = title  ,
                     prefix    = prefix ,
                     alignment = 10*'c' ,
                     style     = style  )

# =============================================================================
## Draw ECDF + 2 lines when/if t-value specified
#  @code
#  ecdf   = ...
#  tvalue = ...
#  result = draw_ecdf ( ecdf , tvalue = tvalue ) 
#  @endcode 
def draw_ecdf ( ecdf          ,
                tvalue = None ,
                option = ''   , *options , **kwargs ) :
    """ Draw ECDF + 2 lines when/if t-value specified
    >>> ecdf   = ...
    >>> tvalue = ...
    >>> result = draw_ecdf ( ecdf , tvalue = tvalue ) 
    """

    has_tvalue  = isinstance ( tvalue , num_types )
    
    xmin , xmax = ecdf.xmin () , ecdf.xmax ()
    
    if has_tvalue :
        tvalue  = float ( tvalue ) 
        xmin    = min ( xmin , tvalue )
        xmax    = max ( xmax , tvalue )

    delta  = xmax - xmin
    xmin  -= 0.10 * delta 
    xmax  += 0.10 * delta 
 
    xmin , xmax = axis_range ( xmin , xmax , delta = 0.20 )

    ## some transformation  
    kw = cidict ( transform = cidict_fun , **kwargs )
        
    kw [ 'xmin'      ] = kw.pop ( 'xmin'       , xmin   ) 
    kw [ 'xmax'      ] = kw.pop ( 'xmax'       , xmax   )
    kw [ 'color'     ] = kw.pop ( 'linecolor'  , Orange )
    kw [ 'linewidth' ] = kw.pop ( 'linewidth'  , 2      )
    kw [ 'maxvalue'  ] = kw.pop ( 'maxvalue'   , 1.1    )
    kw [ 'minvalue'  ] = kw.pop ( 'minvalue'   , 1e-6   )
    kw [ 'copy'      ] = kw.pop ( 'copy'       , True   )

    result = ecdf.draw  ( option = option , *options , **kw )
    
    ## draw ECDF 
    if not  has_tvalue : return result
    
    ## vertical line 
    vline     = ROOT.TLine ( tvalue , 1e-3 , tvalue , 1 - 1e-3 )
    
    ## horisontal line 
    xmin      = kw [ 'xmin' ]
    xmax      = kw [ 'xmax' ]
    dx        = ( xmax - xmin ) / 200.0 
    e         = ecdf ( tvalue )
    hline     = ROOT.TLine ( xmin + dx , e , xmax - dx , e )

    ## 
    vline.SetLineWidth  ( 4     ) 
    vline.SetLineColor  ( Green )
    ## 
    hline.SetLineWidth  ( 2     ) 
    hline.SetLineColor  ( Blue  )  
    hline.SetLineStyle  ( 9     ) 
    ##
    ecdf.__lines = vline , hline 
    ##
    ## ROOT.SetOwnership ( vline , False ) 
    ## ROOT.SetOwnership ( hline , False )
    ## 
    vline.draw ( 'same' , copy = True )
    hline.draw ( 'same' , copy = True )
    ## 
    return result, vline, hline 

# ==============================================================================
## Convert numpy-array statistics into Ostap::SVectorWithError
#  @see Ostap::Math::SVectorWithError
def np2vct ( data , weight = None ) :
    """ Convert numpy-array statistics into `Ostap.Math.SVectorWithError`
    - see `Ostap.Math.SVectorWithError`
    """

    shape     = data.shape
    n , N     = shape
    
    w         = None if weight_trivial ( weight ) else weight
    
    from statsmodels.stats.weightstats import DescrStatsW as DSW 
    dsw       = DSW  ( data , weights = w )
    mean      = dsw.mean
    covmtrx   = dsw.cov

    ## load corresponding linear-algebra tricks:
    MN        = Ostap.Math.SymMatrix(N)
    
    ## prepare output 
    RT        = Ostap.Math.SVectorWithError [ N ]
    values    = RT.Value      () 
    covs      = RT.Covariance () 
    
    for i in range ( N  ) :
        values [ i     ] = float ( mean    [ i ]       )
        covs   [ i , i ] = float ( covmtrx [ i ] [ i ] )
        for j in range (  0 , i  ) :
            cij = float ( covmtrx [ i ] [ j ]  )
            cji = float ( covmtrx [ j ] [ i ]  )                
            cc  = 0.5 * ( cij + cji ) 
            covs [ i , j ] = cc
            covs [ j , i ] = cc  
        
    return RT ( values ,  covs ) 
        
# ==============================================================================
from scipy.stats import combine_pvalues as _combine_pvs_ # =====================
# ==============================================================================
if Version ( "1.10" ) <= Version ( scipy.__version__ ) : # =====================
    # ==========================================================================
    ## combine p-values using certain method 
    def _combine_pvalues ( data , method = "fisher" ) :
        """ Combine p-values using certain method """
        return _combine_pvs_ ( data , method = method ).pvalue
    # ==========================================================================
else : # =======================================================================
    # ==========================================================================
    ## combine p-values using certain method 
    def _combine_pvalues ( data , method = "fisher" ) :
        """ Combine p-values using certain method """
        return _combine_pvs_ ( data , method = method ) [ 1 ]
    
# ==============================================================================
## combine p-values 
#  - use toys to propagate the uncertainties 
def combine_pvalues ( pvalues , method , tol = 1.e-8 , N = 400 ) :
    """ Combine p-values
    - use toys to propagate the uncertainties 
    """
    from   scipy.stats              import combine_pvalues as _combine_pvs
    ##
    if any ( isinstance ( p , VE ) and 0 < p.cov2() and 0 <= p.value() <= 1 for p in pvalues ) :
        pvs = [ clip_pvalue ( VE ( p ) , 0.5 ) for p in pvalues ]
        ## toys
        cnt = SE()
        for i in range ( N )  :
            ## sampled 
            spv  = [ p.gauss ( accept = lambda v : 0 < v < 1 ) if 0 < p.cov2() else float ( p ) for p in pvs ]            
            cnt += _combine_pvalues ( spv , method = method )

        rms = cnt.rms()        
        return VE ( cnt.mean() , rms * rms ) 
    
    pvs = ( min ( max ( tol , float ( p ) ) , 1 - tol ) for p in pvalues )            
    return _combine_pvalues ( pvs , method = method )

# ==============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================


    
