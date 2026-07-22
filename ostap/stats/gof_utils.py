#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof_utils.py
#  Set of utulities for goodness-of-fit studies 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2023-12-06
# =============================================================================
""" Simple utulities for goodness-of-fit studies 
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
    'nearest_distances'  , ## get array of all nearest distances for the dataset 
    'nearest_neighbors'  , ## get all nearest neigbours
    's2u'                , ## convert structured numpy array into non-structured
    'combine_pvalues'    , ## combine p-values 
    
)
# =============================================================================
from   ostap.core.meta_info     import root_info 
from   ostap.core.core          import Ostap, VE, SE 
from   ostap.utils.cidict       import cidict, cidict_fun
from   ostap.core.ostap_types   import string_types, num_types, numpy_floats  
from   ostap.math.math_base     import doubles, axis_range
from   ostap.math.math_ve       import significance
from   ostap.math.ve            import fmt_pretty_ve
from   ostap.math.math_base     import pos_infinity, weight_trivial      
from   ostap.stats.counters     import SE, WSE, EffCounter
from   ostap.utils.basic        import ( numcpu   , num_jobs     , 
                                         typename , run_parallel ) 
from   ostap.utils.utils        import splitter
from   ostap.utils.memory       import memory_enough 
from   ostap.utils.progress_bar import progress_bar
from   ostap.logger.symbols     import times, plus_minus, greek_lower_sigma
from   ostap.logger.pretty      import pretty_float
from   ostap.plotting.color     import Orange, Green, Blue
from   packaging.version        import Version 
import ostap.logger.table       as     T 
import ROOT, os, sys, math, numpy, scipy   
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gof_utils' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-fit studies' )#
# =============================================================================
if  ( 6 , 32 ) <= root_info : data2vct = lambda s : s
else                        : data2vct = lambda s : doubles ( s ) 
# ============================================================================
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
        """
        return _sk_pw_distances ( data1 , data2 , metric , **kwargs ).flatten()
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    has_sklearn = False # =====================================================
    # =========================================================================
# =============================================================================
## (1) use scipy 
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import scipy # ============================================================
    import scipy.spatial # ====================================================
    # =========================================================================   
    if Version ( "1.6.0" ) <=  Version ( scipy.__version__ ) : # ==============
        # =====================================================================
        ## Get nearest distances using scipy.spatial.cKDTree 
        def nearest_distances ( data        ,
                                n_jobs = -1 , **config ) :
            """ Get nearest distances using scipy.spatial.cKDTree 
            """
            config [ 'workers' ] = n_jobs 
            tree          = scipy.spatial.cKDTree ( data )
            distances , _ = tree.query ( data , k = [ 2 ]  , **config )
            ## distances     = distances [ : , 1 ]  # DNN (Distance to Nearest Neighbour 
            return distances.flatten() 
        # ====================================================================
    else : # =================================================================
        # ====================================================================
        ## Get nearest distances using scipy.spatial.cKDTree         
        def nearest_distances ( data , n_jobs = -1 , **config ) :
            """ Get nearest distances using scipy.spatial.cKDTree 
            """
            tree          = scipy.spatial.cKDTree ( data )
            distances , _ = tree.query ( data , k = 2 , **config )
            distances     = distances [ : , 1 ]  # DNN (Distance to Nearest Neighbour 
            return distances.flatten() 
        # =====================================================================
except ImportError : # ========================================================
    # =========================================================================
    pass # ====================================================================
# =============================================================================
## (2) make a try to use NearestNeighbours from sklearn 
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import sklearn # ==========================================================
    import sklearn.neighbors # ================================================
    # =========================================================================
    ## Get nearest distances using sklearn 
    def nearest_distances ( data , n_jobs = -1 , **config ) :
        """ Get nearest distances using sklearn 
        """
        nn = sklearn.neighbors.NearestNeighbors ( n_jobs      = n_jobs ,
                                                  n_neighbors = 2      , **config )
        nn.fit ( data )
        distances ,  _  = nn.kneighbors( data )
        distances       = distances [ : , 1 ]  # DNN (Distance to Nearest Neighbour 
        return distances.flatten()
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    has_sklearn = False # =====================================================
    # =========================================================================
# =============================================================================
## (1) use scipy
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import scipy # ============================================================
    import scipy.spatial # ====================================================
    if Version ( "1.6.0" ) <= Version ( scipy.__version__  ) : # ==============
        # =====================================================================
        ## Get nearest neighbors using scipy.spatial.cKDTree 
        def nearest_neighbors ( data , n_jobs = -1 , n_neighbors = 10 , **config ) :            
            """ Get nearest neighbors using scipy.spatial.cKDTree
            """
            config [ 'workers' ] = n_jobs 
            tree        = scipy.spatial.cKDTree ( data )
            _ , indices = tree.query ( data , k = n_neighbors + 1 , **config )
            return indices [:, 1: ]
        # ====================================================================
    else : # =================================================================
        # ====================================================================
        ## Get nearest neighbors using scipy.spatial.cKDTree 
        def nearest_neighbors ( data , n_jobs = -1 , n_neighbors = 10 , **config ) :            
            """ Get nearest neighbors using scipy.spatial.cKDTree
            """
            config [ 'workers' ] = n_jobs 
            tree        = scipy.spatial.cKDTree ( data )
            _ , indices = tree.query ( data , k = n_neighbors + 1 , **config )
            return indices [:, 1: ]
    # ========================================================================
except ImportError : # =======================================================
    # ========================================================================
    pass
# ============================================================================
## (2) make try to use sklearn 
# ============================================================================
try : # ======================================================================
    # ========================================================================
    import sklearn
    import sklearn.neighbors 
    # =========================================================================
    ## Get nearest neighbors using scipy.spatial.cKDTree 
    def nearest_neighbors ( data , n_jobs = -1 , n_neighbors = 10 , **config ) :            
        """ Get nearest neighbors using scipy.spatial.cKDTree
        """
        nn = sklearn.neighbors.NearestNeighbors ( n_jobs = n_jobs , n_neighbors = n_neighbors , **config )
        nn.fit ( data )
        _ , indices  = nn.kneighbors( data )
        return indices [ :, 1: ]
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    has_sklearn = False # =====================================================
    # =========================================================================
    
  

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
## Normalize several numpy arrays such that mean and rms for the pooled sample
#  are equal to 0 and 1 correspondinly
#  @code
#  ds1 = ...
#  ds2 = ...
#  ds3 = ...
#  ds1 , ds2 , ds3 = pool_normalize ( ds1 , ds2 , ds3 ) 
#  @endcode 
def normalize_pooled ( *datasets ) :
    """ Normalize several numpy arrays such that the mean and rms for the POOLED sample
    are equal to 0 and 1 correspondinly
    
    >>> ds1 = ...
    >>> ds2 = ...
    >>> ds3 = ...
    >>> ds1 , ds2 , ds3 = np_normalize ( ds1 , ds2 , ds3 )
    
    """

    if not datasets : return ()

    ## print ( 'SKIP NORMALIZATION' ) 
    ## return datasets 
    
    total      = 0 
    total_mean = None
    total_std2 = None
    
    sizes      = []
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
## Short labels for various statitical estimators 
Labels = {
    'KS'  : 'Kolmogorov-Smirnov' ,
    'K'   : 'Kuiper'             ,
    'AD'  : 'Anderson-Darling'   ,
    'CM'  : 'Cramer-von Mises'   ,
    'ZK'  : 'Zhang/ZK'           ,
    'ZA'  : 'Zhang/ZA'           ,
    'ZC'  : 'Zhang/ZC'           ,        
    'BJ'  : 'Berk-Jones'         ,        
    'NLL' : '-log L'             ,
    'AIC' : 'Aikaike IC'         ,
    'BIC' : 'Bayesian IC'        ,
    
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
    'AIC' : ( 'aic' , 'aikaike'    ) ,  
    'BIC' : ( 'bic' , 'bayesian'   ) , 
    }
# =============================================================================
assert Labels.keys() == Keys.keys() , "Mismatch between Labels & Keys structures!"
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
## Generator of permutations for datasets (and associated weights)
#  1. pooled dataset is created
#  2. pooled dataset is split into two permuted datasets 
def make_permutations ( nToys    ,
                        data1    ,
                        data2    , *    ,
                        weight1  = None ,
                        weight2  = None ,
                        silent   = True , 
                        progress = True ) :
    """ Generator of permutations
    >>> ds1, ds2 = ...
    >>> gof      = ...
    >>> nToys    = 100 
    >>> for d1,d2,w1,w2  in make_permutations ( nToys , ds1 , ds2 ) :
    ...
    """
    assert isinstance ( nToys , int ) and 1 <= nToys , "make_permutations: Invalid `nToys`: %s" % nToys 

    n1     = len ( data1 )
    n2     = len ( data2 )
    pooled = numpy.concatenate ( [ data1 , data2 ] )

    w1_trivial = weight_trivial ( weight1 )
    w2_trivial = weight_trivial ( weight2 )
    
    if w1_trivial and w2_trivial : weights = None 
    else :
        if w1_trivial : weight1 = numpy.ones ( n1 )
        if w2_trivial : weight2 = numpy.ones ( n2 )
        weights = numpy.concatenate ( [ weight1 , weight2 ] , dtype = pooled.dtype ).reshape ( ( n1 + n2 , 1 ) )
        pooled  = numpy.hstack      ( [ pooled       , weights      ] )

    
    ## run permutations 
    for i in progress_bar ( nToys , silent = silent and not progress , description = 'Permutations:' ) :
        
        numpy.random.shuffle ( pooled )
        
        ds1 = pooled [    : n1 ]
        ds2 = pooled [ n1 :    ]
        
        if weights is None :
            ds1 , w1 = ds1 , None
            ds2 , w2 = ds2 , None
        else               :                
            ds1 , w1 = ds1 [ : , : -1 ] , ds1 [ : , -1 ]
            ds2 , w2 = ds2 [ : , : -1 ] , ds2 [ : , -1 ]
            
        yield ds1 , ds2 , w1 , w2 
            
    del pooled
    del weights
                    
# =============================================================================
## @class PERMUTATOR
#  Helper class that allow to run permutation test in parallel 
class PERMUTATOR(object) :
    """ Helper class that allow to run permutation test in parallel 
    """
    def __init__ ( self    ,
                   gof     ,
                   t_value ,
                   ds1     ,
                   ds2     ,
                   weight1 = None ,
                   weight2 = None ) :
        
        self.gof     = gof
        self.ds1     = ds1
        self.ds2     = ds2

        w1_trivial = weight_trivial ( weight1 )
        w2_trivial = weight_trivial ( weight2 )
        
        if w1_trivial and w2_trivial :
            weight1 , weight2 = None, None 
        else :
            if w1_trivial : weight1 = numpy.ones ( len ( self.ds1 ) , dtype = float )
            if w2_trivial : weight2 = numpy.ones ( len ( self.ds2 ) , dtype = float )
        
        self.weight1 = weight1
        self.weight2 = weight2 
        self.t_value = t_value
        self.ecdf    = None
        
    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object
        """
        return { 'gof'        : self.gof      ,
                 'ds1'        : self.ds1      ,   
                 'ds2'        : self.ds2      ,   
                 'weight1'    : self.weight1  ,
                 'weight2'    : self.weight2  ,   
                 't_value'    : self.t_value  , 
                 'ecdf'       : self.ecdf     }
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object
        """
        self.gof     = state.pop ( 'gof'        )
        self.ds1     = state.pop ( 'ds1'        )
        self.ds2     = state.pop ( 'ds2'        )        
        self.weight1 = state.pop ( 'weight1'    )
        self.weight2 = state.pop ( 'weight2'    )
        self.t_value = state.pop ( 't_value'    )
        self.ecdf    = state.pop ( 'ecdf'       )
        
    # =========================================================================
    ## run N-permutations 
    def __call__ ( self , N , silent = True , progress = False ) :
        
        counter, tvalues = self.run_toys ( N = N , silent = silent , progress = progress )
        
        if not self.ecdf : self.ecdf = Ostap.Math.ECDF ( tvalues , True )
        else             : self.ecdf.add ( data2vct ( tvalues )  )
        
        return counter 

    # =========================================================================
    ## run N-toys 
    def run_toys ( self, N , silent = True , progress = False ) :
        """ Run N-toys
        """        
        counter = EffCounter()
        tvalues = []

        for data1 , data2 , weight1 , weight2 in make_permutations ( N                       ,
                                                                     self.ds1                ,
                                                                     self.ds2                ,
                                                                     weight1  = self.weight1 ,
                                                                     weight2  = self.weight2 , 
                                                                     progress = progress     ,   
                                                                     silent   = silent       ) :

            ## attention: normalize = False 
            tv       = self.gof.tvalue ( data1     ,
                                         data2     ,
                                         weight1   = weight1 ,
                                         weight2   = weight2 ,
                                         normalize = False   )
            tv       = float ( tv )
            tvalues.append   ( tv  )
            counter += bool  ( self.t_value < tv  )

        return counter, tuple ( tvalues )

    # =========================================================================
    ## Run NN-permutations in parallel using the default WorkManager
    def run ( self , nToys , silent = False , progress  = True ) :
        """ Run permutations in parallel using WorkManager
        """
        me       = math.ceil ( memory_enough() ) + 1 
        njobs    = min ( 2 * numcpu () + 3 , me ) 
        the_list = [ n for n in splitter ( nToys , njobs ) ] 
        njobs    = len ( the_list ) 
        
        if not silent :
            logger.info ( 'GoF-permutations: #%d parallel subjobs to be used with WorkManager' % njobs )
            
        counter = EffCounter()
        tvalues = () 
        ## 
        ## use *BARE* interface here
        from ostap.parallel.parallel import WorkManager
        with WorkManager ( silent = silent ) as manager : 
            for result in manager.iexecute ( self.run_toys ,
                                             the_list      ,
                                             progress      = progress        ,
                                             njobs         = njobs           ,
                                             description   = 'Permutations:' ) :
                cnt , tvals = result 
                counter += cnt
                tvalues += tvals 
        ##
        if not self.ecdf : self.ecdf = Ostap.Math.ECDF ( tvalues , True )
        else             : self.ecdf.add    ( data2vct ( tvalues )     )
        ##
        return counter

# =============================================================================
## @class TOYS
#  Helper class to run toys for Goodness-of-Fit studies 
class TOYS(object) :
    """ Helper class that allow to run toys in parallel 
    """
    def __init__ ( self    , 
                   gof     , *        , 
                   t_value            ,
                   pdf                ,
                   Ndata              , 
                   sample     = False ,
                   parameters = {}    ) :
        
        self.gof        = gof
        self.pdf        = pdf
        self.Ndata      = Ndata 
        self.t_value    = t_value
        self.sample     = gof.sample 
        self.silent     = gof.silent

        if parameters : self.parameters = parameters
        else          : self.parameters = pdf.params() 
                    
        self.__ecdf    = None 

    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object
        """
        self.pdf.load_params ( self.parameters , silent = True )        
        return { 'gof'        : self.gof        ,
                 'pdf'        : self.pdf        ,
                 'Ndata'      : self.Ndata      ,
                 't_value'    : self.t_value    , 
                 'sample'     : self.sample     , 
                 'silent'     : self.silent     ,                  
                 'parameters' : self.parameters ,
                 'ecdf'       : self.ecdf       } 

    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object
        """
        self.gof         = state.pop ( 'gof'        )
        self.pdf         = state.pop ( 'pdf'        )
        self.Ndata       = state.pop ( 'Ndata'      )
        self.t_value     = state.pop ( 't_value'    )
        self.sample      = state.pop ( 'sample'     )
        self.silent      = state.pop ( 'silent'     )
        self.parameters  = state.pop ( 'parameters' )
        self.__ecdf      = state.pop ( 'ecdf'       )

        ## (1) re-load parameters 
        self.pdf.load_params ( self.parameters , silent = True )


    # =========================================================================
    ## run N-toys 
    def __call__ ( self , nToys , progress = True  ) :
        """ Run N-toys
        """
        counter , ecdf = self.run_toys ( nToys = nToys , progress = progress )
        return counter 
    
    # =========================================================================
    ## run N-toys 
    def run_toys ( self , nToys , progress = False ) :
        """ Run N-toys
        """
        ROOT.gRandom                     .SetSeed () 
        ROOT.RooRandom.randomGenerator() .SetSeed ()

        counter = EffCounter ()
        tvalues = [] 
        for i in progress_bar ( nToys , description = "Toys:" , silent = not progress ) : 

            ## for consistency
            self.pdf.load_params ( self.parameters , silent = True )
            
            dset     = self.pdf.generate ( self.Ndata , sample = self.sample )
            tv       = self.gof ( self.pdf , dset )
            counter += bool ( self.t_value > tv   ) ## NOTE THE SIGN HERE!

            tvalues.append ( tv ) 

            if isinstance  ( dset , ROOT.RooDataSet ) :
                dset.clear ()
                ROOT.SetOwnership ( dset , True )
                del dset

        tvalues = tuple ( tvalues ) 
        if not self.ecdf : self.__ecdf = Ostap.Math.ECDF ( tvalues  , True ) 
        else             : self.ecdf.add ( data2vct ( tvalues )  )

        return counter, self.ecdf 

    # =========================================================================
    ## Run N-toys in parallel using WorkManager
    def run ( self , nToys , silent = False , progress = True ) :
        """ Run toys in parallel using WorkManager
        """
        ##
        assert isinstance ( nToys , int ) and 1 <= nToys , "Invalid nToys: %s" % nToys
        ##
        ## how many processes fits into memory ?
        me       = math.ceil ( memory_enough() ) + 1 
        njobs    = min ( 2 * numcpu () + 3 , me ) 
        the_list = [ n for n in splitter ( nToys , njobs ) ] 
        njobs    = len ( the_list )
        
        if not silent :
            logger.info ( 'GoF-toys: #%d parallel subjobs to be used' % njobs )
            ##
            
        counter = EffCounter()
        tvalues = ()
        ##        
        ## use *BARE* interface here 
        from ostap.parallel.parallel import WorkManager
        with WorkManager ( silent = silent ) as manager :            
            for result in manager.iexecute ( self.run_toys ,
                                             the_list      ,
                                             progress      = progress   ,
                                             njobs         = njobs      ,
                                             description   = 'Toys:'    ) :

                cnt , ecdf = result
                
                counter += cnt
                
                if not self.__ecdf    : self.__ecdf    =  ecdf 
                else                  : self.__ecdf.add ( ecdf )
                
        return counter 

    @property 
    def ecdf ( self ) :
        """`ecdf` : empirical CDF for t-values from toys/pseudoexperiments 
        """
        return self.__ecdf
    
# =============================================================================
pvalue_types = num_types + ( VE , ) 
# =============================================================================
## Format the row for GoF tables
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
    """ Format the row for the GoF tables
    >>> tvalue = ...
    >>> pvalue = ...
    >>> ecdf   = ...
    >>> header , row = format_row ( tvalue = tvalue , pvalue = pvalue  , ecdf = ecdf )
    """
    
    has_tvalue  = not tvalue  is None and isinstance ( tvalue  , num_types       ) 
    has_pvalue  = not pvalue  is None and isinstance ( pvalue  , pvalue_types    ) 
    has_ecdf    = not ecdf    is None and isinstance ( ecdf    , Ostap.Math.ECDF ) 
    has_counter = not counter is None and isinstance ( counter , ( SE , WSE )    )

    if has_ecdf  and not has_counter  :
        counter     = ecdf.counter ()
        has_counter = True 
        
    if has_tvalue and has_pvalue and has_counter :
        
        header = ( 't-value'    ,
                   't-mean'     ,
                   't-rms'      ,
                   't-min/max'  ,                
                   '%s[..]' % times , 'p-value [%]' , '#%s' % greek_lower_sigma ) 
        
        mean       = counter.mean   ()
        rms        = counter.rms    () 
        vmin, vmax = counter.minmax () 
        
        mxv        = max ( abs ( tvalue       ) ,
                           abs ( mean.value() ) ,
                           mean.error()         , rms ,
                           abs ( vmin )  , abs ( vmax ) ) 
        
        fmt, fmtv , fmte , expo = fmt_pretty_ve ( VE ( mxv ,  mean.cov2() ) ,
                                                  precision   = precision   ,
                                                  width       = width       , 
                                                  parentheses = False       )
        
        if expo : scale = 10**expo
        else    : scale = 1
    
        vs  = tvalue / scale
        vm  = mean   / scale
        vr  = rms    / scale
        vmn = vmin   / scale
        vmx = vmax   / scale
        
        fmt2 = '%s/%s' % ( fmtv , fmtv ) 

        ##
        
        pv      = clip_pvalue  ( pvalue ) 
        nsigma  = significance ( pv     ) ## convert  it to significance
        ## 
        if isinstance ( nsigma , VE ) and nsigma.cov2 () <= 0 : nsigma = float ( nsigma ) 
        if 50 <= float ( nsigma ) : nsigma = pos_infinity 

        pvalue  = pvalue * 100

        if isinstance ( pvalue , VE ) : pvalue  = '%5.2f %s %.2f' % ( pvalue.value() , plus_minus , pvalue.error () )
        else                          : pvalue  = '%5.2f'         % float ( pvalue ) 
        if isinstance ( nsigma , VE ) : nsigma  = '%.2f %s %.2f'  % ( nsigma.value() , plus_minus , nsigma.error () )
        else                          : nsigma  = '%.2f'          % float ( nsigma ) 
        
        row = ( fmtv  % vs ,
                fmt   % ( vm.value() , vm.error() ) ,
                fmtv  % vr                          ,
                fmt2  % ( vmn , vmx )               ,
                ( '%s10^%+d' %  ( times , expo )  if expo else '' ) , pvalue , nsigma )
        
        return header , row

    elif has_tvalue and has_pvalue :
        
        header = ( 't-value'  , '%s[..]' % times , 'p-value [%]' , '#%s' % greek_lower_sigma ) 

        pv        = clip_pvalue  ( pvalue )
        nsigma    = significance ( pv     )
        if isinstance ( nsigma , VE ) and nsigma.cov2 () <= 0 : nsigma = float ( nsigma ) 
        if 50 <= float ( nsigma ) : nsigma = pos_infinity 
        
        tv , expo = pretty_float ( tvalue , precision = precision , width = width )

        pvalue  = pvalue * 100 
        if isinstance ( pvalue , VE ) : pvalue  = '%5.2f %s %.2f' % ( pvalue.value() , plus_minus , pvalue.error () )
        else                          : pvalue  = '%5.2f'         % float ( pvalue ) 
        if isinstance ( nsigma , VE ) : nsigma  = '%.2f %s %.2f'  % ( nsigma.value() , plus_minus , nsigma.error () )        
        else                          : nsigma  = '%.2f'          % float ( nsigma ) 
        
        row     = tv , '%s10^%+d' % ( times , expo ) if expo else '' , pvalue , nsigma

        return header , row 

    elif has_tvalue :
        
        header    = ( 't-value'  , '%s[..]' % times )          
        tv , expo = pretty_float ( tvalue , precision = precision , width = width )
        row       = tv , '%s10^%+d' % ( times , expo ) if expo else ''        
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


# =============================================================================
## use LigthGBM ?
#  - there is some mess with lightgbm&narwhals installation 
def useLightGBM () :
    """ Use LigthGBM ?
    - there is soem mess with ligthgbm&narwhals installation 
    """
    # ============================================================================
    try : # ======================================================================
        # ========================================================================
        import lightgbm
        logger.info ( 'LightGBM version : %s' % lightgbm.__version__ ) 
        if Version ( lightgbm.__version__ ) <  Version ( "4.7.0"  ) : return True
        import narwhals
        logger.info ( 'Narwhals version : %s' % narwhals.__version__ ) 
        return Version ( "2.0" ) <= Version ( narwhals.__version__ )
        # ========================================================================
    except ImportError : # =======================================================
        # ========================================================================
        return False 
    
# ===============================================================================
## use XGBoost ?
def useXGBoost () : 
    """ Use XGBoost
    """
    # ==========================================================================
    try : # ====================================================================
        # ======================================================================
        import xgboost        
        return Version ( "1.0" ) <= Version ( xgboost.__version__ )
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False 

# ===============================================================================
## use CatBoost ?
def useCatBoost () : 
    """ Use CatBoost
    """
    from ostap.core.cpu_info import HAS_AVX2
    if not HAS_AVX2 : return  False 
    # ==========================================================================
    try : # ====================================================================
        # ======================================================================
        import catboost
        return True 
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================


    
