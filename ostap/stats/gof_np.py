#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof_np.py
#  Set of utilities for goodness-of-fit studies for multidimensional fits
#  @see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
#  @see https://doi.org/10.1088/1748-0221/5/09/P09004
#  @see http://arxiv.org/abs/arXiv:1003.1768 
#  @author Artem Egorychev Artem.Egorychev@cern.ch 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-09-16
# =============================================================================
""" Simple utilities for goodness-of-fit studies for multidimensional fits 
- see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
- see https://doi.org/10.1088/1748-0221/5/09/P09004
- see http://arxiv.org/abs/arXiv:1003.1768
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-29"
__all__     = (
    'GoFnp'           , ## A base class for numpy-related family of methods to probe goodness-of-fit
    ##
    'MIXnp'           , ## Mixed samples                 Goodness-of-Fit method 
    'PPDnp'           , ## Point-to-Point Dissimilarity  Goodness-of-Fit method 
    'DNNnp'           , ## Distance-to-Nearest-Neighbor  Goodness-of-Fit method
    ##
    'KullbackLeibler' , ## Very crude estimator based on Kullback-Leibler's divergency 
    'Jeffrey'         , ## Very crude estimator based on Jeffrey's divergency 
    'JensenShannon'   , ## Very crude estimator based on Jensen-Shannon divergency 
    'Mahalanobis'     , ## Very crude estimator based on Mahalanobis' distance
    'Hotelling'       , ## Very crude estimator based on Hotelling's distance
    'Bhattacharyya'   , ## Very crude estimator based on Bhattacharyya's divergency 
    'Wasserstein'     , ## Very crude estimator based on Wasserstein's divergency 
    'Hellinger'       , ## Very crude estimator based on Hellinger's divergency
    ## 
    'Chi2'            , ## binned 1,2&3D chi2
    ## 
)
# =============================================================================
from   ostap.core.ostap_types   import string_types, num_types 
from   ostap.core.core          import Ostap, hID  
from   ostap.utils.utils        import split_n_range
from   ostap.utils.core         import typename
from   ostap.utils.basic        import numcpu, num_jobs, run_parallel 
from   ostap.utils.config       import Config
from   ostap.stats.gof          import AGoFnp
from   ostap.stats.gofnp        import GoFnp
from   ostap.stats.utils        import ( weight_trivial     ,
                                         check_all          , 
                                         valid_data_shape   ,
                                         num_features       ,
                                         num_samples        ,
                                         np2vct             ) 
from   ostap.stats.gof_utils    import ( normalize_pooled   ,
                                         pairwise_distances ,
                                         nearest_neighbors  , 
                                         nearest_distances  , 
                                         draw_ecdf          , s2u )
import ostap.logger.symbols    as      S
import ostap.math.math_base           
import ROOT, numpy
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gof_np' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-fit studies for multidimensional fits' )
# =============================================================================
method_MIX = ( S.to_script ( 'MIX' ) + ' ' ) if S.show else 'MIX'
method_PPD = ( S.to_script ( 'PPD' ) + ' ' ) if S.show else 'PPD'
method_DNN = ( S.to_script ( 'DNN' ) + ' ' ) if S.show else 'DNN'
method_DNN = ( S.to_script ( 'DNN' ) + ' ' ) if S.show else 'DNN'
method_KL  = ( S.to_script ( 'KL'  ) + ' ' ) if S.show else 'Kullback-Leibler'
method_J   = ( S.to_script ( 'J'   ) + ' ' ) if S.show else 'Jeffrey'
method_JS  = ( S.to_script ( 'JS'  ) + ' ' ) if S.show else 'Jensen-Shannon'
method_M   = ( S.to_script ( 'M'   ) + ' ' ) if S.show else 'Mahalanobis'
method_T2  =   S.T2                          if S.show else 'Hotelling' 
method_B   = ( S.to_script ( 'B'   ) + ' ' ) if S.show else 'Bhattacharyya'
method_W2  =   S.W2                          if S.show else 'Wasserstein' 
method_H   = ( S.to_script ( 'H'   ) + ' ' ) if S.show else 'Hellinger'

# ============================================================================
## define configuration for psi-function for PPD method
#   - distance type of <code>cdist</code>
#   - transformation function for `pairwise_distance' soutput
#   - increasing function ?
#   @code
#   distance_type , transform, increasing = psi_conf ( 'linear' )
#   @endcode
def psi_conf ( psi , scale = 1.0 ) :
    """ Define configuration for psi-function for PPD method
    """

    if   psi in ( 'euclidean'   , 'linear'   ) :                           ## psi = x 
        return 'euclidean'      , None                                   , True 
    elif psi in ( 'sqeuclidean' , 'squared'  ) :                           ## psi = x**2 
        return 'sqeuclidean'    , None                                   , True 
    elif psi in ( 'inverse'     , 'coulomb'  ) :                           ## psi = 1/x 
        return 'euclidean'      , lambda x : -1.0 / ( x [ 0 < x ] )      , True  
    elif psi in ( 'inverse2'    , 'coulomb2' ) :                           ## psi = 1/x**2 
        return 'sqeuclidean'    , lambda x : -1.0 / ( x [ 0 < x ] )      , True  
    elif psi in ( 'log'         , 'logarithm'    ) :                       ## psi = log(x)
        return 'sqeuclidean'    , lambda x :   numpy.log ( x [ 0 < x ] ) , True  
    elif psi in ( 'gauss'       , 'gaussian'     ) :                       ## psi = exp (-x*x/0.5)
        return 'sqeuclidean'    , lambda x :  -numpy.exp ( scale * x   ) , True
    elif isinstance ( psi , string_types ) :
        return psi , None , True         

    raise TypeError ( "Unknown `psi':%s" % psi ) 

# =============================================================================
## @class MIXnp
#  Implementation of `Mixed Sample' method for probing the Goodness-Of-Fit
#  @see M.Williams, "How good are your fits?
#       Unbinned multivariate goodness-of-fit tests in high energy physics"
#  @see https://doi.org/10.1088/1748-0221/5/09/P09004
#  @see http://arxiv.org/abs/arXiv:1003.1768
#
#  M.Williams writes:
#     The method <...> is easy to use and conceptually it is easy to understand.
#     It is excellent at rejecting large localized discrepancies but fairly poor
#     at rejecting small omnipresent ones.  The p-values can be calculated analytically.
#     This method would make a nice addition to the high energy physics g.o.f. toolkit.
class MIXnp(GoFnp) :
    """ Implementation of `Mixed Sample' for probing the Goodness-Of-Fit
    - see M.Williams, "How good are your fits?
       Unbinned multivariate goodness-of-fit tests in high energy physics"
    - see https://doi.org/10.1088/1748-0221/5/09/P09004
    - see http://arxiv.org/abs/arXiv:1003.1768
    
    M.Williams writes:
    The method <...> is easy to use and conceptually it is easy to understand.
    It is excellent at rejecting large localized discrepancies but fairly poor
    at rejecting small omnipresent ones.  The p-values can be calculated analytically.
    This method would make a nice addition to the high energy physics g.o.f. toolkit.
    """
    
    def __init__ ( self ,
                   nToys       = 1000  ,
                   parallel    = True  , 
                   n_neighbors = 10    , **params ) : 
        
        ## Attention!
        assert isinstance ( n_neighbors , int ) and 2 <= n_neighbors , \
            "Invalid `n_neighbors`: %s" % n_neighbors 

        ## store it
        self._k_max = n_neighbors 
        ##
        
        ## switch off parallel processing 
        if parallel and not run_parallel ( parallel ) :
            logger.info ( '%s: (internal) parallel processing is OFF' % typename ( self ) )
            parallel = False
            
        ## (re)define n_jobs     
        params [ 'n_jobs' ] = 1 if parallel else num_jobs ( params , numcpu () - 1 )

        ## initialize the base 
        GoFnp.__init__ ( self                      , 
                         nToys        = nToys      ,
                         parallel     = parallel   , 
                         method       = method_MIX ,
                         normalize    = True       , 
                         n_neighbors  = self.k_max , **params )

    # =========================================================================
    ## Are weights supported by this estimator?
    @property
    def weights_supported ( self ) :
        """`weights_supported` : Are weights supported by this estimator?
        """
        return False 
    
    # =========================================================================
    ## Good for two-samples comparison?
    #  Can this estimator be used for comparison of two samples?
    @property 
    def two_samples ( self ) :
        """`two_samples`: Can this estimator be used for comparison of two samples?
        """
        return True 
    
    # =========================================================================
    ## k_max` : number fo nearest neighbors to test
    @property
    def k_max ( self ) :
        """`k_max` : number of nearest neighbors to test
        """
        return self._k_max 

    # =========================================================================
    # calculate t-value for (non-structured) 2D arrays
    def tvalue ( self      , 
                 data1     , 
                 data2     , *    , 
                 weight1   = None , 
                 weight2   = None ,
                 normalize = True ) :
        """ Calculate t-value for (non-structured) 2D arrays
        """
        ##
        ## weights 
        w1_trivial = weight_trivial ( weight1 )
        w2_trivial = weight_trivial ( weight2 )

        if not self.weights_supported :
            if not w1_trivial : raise ValueError (  "weight1 must be *trivial*" )
            if not w2_trivial : raise ValueError (  "weight2 must be *trivial*" )

        shape1 = data1.shape
        shape2 = data2.shape
        assert 2 == len ( shape1 ) and 2 == len ( shape2 ) and shape1 [ 1 ]  == shape2 [ 1 ] , \
            "Invalid arrays: %s , %s" % ( shape1 , shape2  )

        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )
        ## normalize
        if normalize and self.normalize : uds1, uds2  = self.normalize_pooled ( uds1 , uds2  ) 

        ## check validity and consitency of ALL input parameters 
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) )
        
        ## 
        n1 = len ( uds1 ) 
        n2 = len ( uds2 ) 
        
        data   = numpy.vstack ( [ uds1  , uds2 ]         )
        labels = numpy.array  ( [ 1 ] * n1 + [ 0 ] * n2  )
        
        ## get the nearest neighbors indices 
        actual_neighbors = nearest_neighbors ( data , **self.params )
        
        source_labels    = labels [ : , numpy.newaxis ] # (N, 1)
        neighbor_labels  = labels [ actual_neighbors  ] # (N, K)

        ## I(i, k) = 1 if neighbor has the same label/source 
        I_ik = ( source_labels == neighbor_labels ) . astype ( int )
        
        result = numpy.sum ( I_ik ) / ( 1.0 * self.k_max * ( n1 + n2 ) )
        return float ( result ) 
    
# =============================================================================
## @class PPDnp
#  Implementation of concrete method "Point-To-Point Dissimilarity"
#  for probing of Goodness-Of-Fit
#  @see M.Williams, "How good are your fits?
#       Unbinned multivariate goodness-of-fit tests in high energy physics"
#  @see https://doi.org/10.1088/1748-0221/5/09/P09004
#  @see http://arxiv.org/abs/arXiv:1003.1768
#
#  M.Williams writes: 
#    The method <...> has excellent rejection power for both large localized
#    discrepancies and small omnipresent ones.  Determining the p-value
#    requires re-sampling the data (using the permutation test) which uses
#    a relatively large amount of processing time.
#    The method is not as easy to understand conceptually as some of
#    the other methods <..> .
#    These downsides are not enough to out-way its excellent performance;
#    this is a very powerful g.o.f. tool.
#
#    @attention It is rather sow to large number of events 
class PPDnp(GoFnp) : 
    """ Implementation of concrete method "Point-To-Point Dissimilarity"
    for probing of Goodness-Of-Fit
    - see M.Williams, "How good are your fits? 
                       Unbinned multivariate goodness-of-fit tests in high energy physics"
    - see https://doi.org/10.1088/1748-0221/5/09/P09004    
    - see http://arxiv.org/abs/arXiv:1003.1768 

    M.Williams writes: 
    ... The method <...> has excellent rejection power for both large localized
    ... discrepancies and small omnipresent ones.  Determining the p-value
    ... requires re-sampling the data (using the permutation test) which uses
    ... a relatively large amount of processing time.
    ... The method is not as easy to understand conceptually as some of
    ... the other methods <..> .
    ... These downsides are not enough to out-way its excellent performance;
    ... this is a very powerful g.o.f. tool.

    - ATTENTION: It is rather sow to large number of events! 
    """
    def __init__ ( self                   ,
                   mc2mc     = False      ,
                   nToys     = 1000       ,
                   psi       = 'gaussian' ,
                   sigma     = 0.10       ,
                   parallel  = True       , 
                   maxsize   = 10000000   , **params ) :
        
        ## switch off parallel processing 
        if parallel and not run_parallel ( parallel ) :
            logger.info ( '%s: (internal) parallel processing is OFF' % typename ( self ) )
            parallel = False            
        params [ 'n_jobs' ]  = 1 if parallel else num_jobs ( params , numcpu() - 1 )


        self.__mc2mc     = True if mc2mc else False
        self.__transform = None
        self.__sigma     = sigma
        self.__psi       = psi
        assert isinstance ( maxsize , int ) and 0 < maxsize , "Invalid `maxsize' : %s" % maxsize
        
        self.__maxsize   = max ( maxsize , 100000  )
        
        ## check validity of `psi`
        scale = -0.5 / ( self.sigma ** 2 ) 
        self.__distance_type , _ , _ = psi_conf ( psi , scale )

        GoFnp.__init__ ( self                   ,
                         nToys     = nToys      ,
                         parallel  = parallel   , 
                         normalize = True       ,
                         method    = method_PPD , **params )
                
    # ==================================================================================
    @property
    def config ( self ) :
        """`config` : get all configuration parameters"""
        conf = super().config 
        conf [ 'mc2mc'            ] = self.mc2mc
        conf [ 'psi'              ] = self.__psi        
        conf [ 'sigma'            ] = self.sigma
        conf [ 'maxsize'          ] = self.__maxsize 
        return conf 
            
    # =========================================================================
    ## Are weights supported by this estimator?
    @property
    def weights_supported ( self ) :
        """`weights_supported` : Are weights supported by this estimator?
        """
        return False

    # =========================================================================
    ## Good for two-samples comparison?
    #  Can this estimator be used for comparison of two samples?
    @property 
    def two_samples ( self ) :
        """`two_samples`: Can this estimator be used for comparison of two samples?
        """
        return True 

    # =========================================================================
    @property
    def mc2mc ( self ) :
        """`mc2mc` : add mc <-->mc distances to the T-value ?
        - when size of the second data set is significantly larger, 
        `mc2mc = False` can be used to speedup calculations 
        """
        return self.__mc2mc
    @property
    def psi  ( self )  :
        """`psi` : psi-function to be used for distance calculation"""
        return self.__psi
    @property
    def sigma ( self ) :
        """`sigma` : `sigma` parameter for gaussian-type of `psi`"""
        return self.__sigma
        
    # =========================================================================
    ## Calculate `sum-of-(transformed)-distances' between all elements in data1 & data2
    def sum_distances ( self, data1 , data2 ) :
        """ Calculate `sum-of-(transformed)-distances' between all elements in data1 & data2
        """
        n1     = len ( data1 )
        n2     = len ( data2 )
        ## if too many distances, process them in chunks
        nnmax  = self.__maxsize 
        if 0 < nnmax < n1 * n2 :
            # ================================================================
            if n1 > n2 : ## swap datasets 
                data1 , data2 = data2 , data1
                n1    , n2    = n2    , n1
            # =================================================================
            result = 0.0
            nsplit = ( n1 * n2 ) // nnmax  + 2
            ## split the second (larger) dataset into `nsplit` parts
            for f , l in split_n_range ( 0 , n2 , nsplit ) :
                result += self.sum_distances ( data1 , data2 [ f : l ] )
            return result 
        ##
        ## how to build distances?
        scale = -0.5 / ( self.sigma ** 2 ) 
        distance_type , transform , _ = psi_conf ( self.psi , scale )
        ##
        ## calculate all pair-wise distances
        distances = pairwise_distances ( data1 , data2 , metric = distance_type , **self.params )        
        distances = distances [ distances > 0 ]
        if transform : distances  = transform ( distances )        
        ##
        return numpy.sum ( distances )
    
    # =========================================================================
    ## Calculate the t-value for 2D arrays
    def tvalue ( self      , 
                 data1     , 
                 data2     , *    , 
                 weight1   = None , 
                 weight2   = None ,
                 normalize = True ) :
        """ Calculate t-value for (non-structured) 2D arrays
        """
        ##
        if not self.weights_supported :
            if not weight_trivial ( weight1 ) : raise ValueError  ( "weight1 must be *trivial*" ) 
            if not weight_trivial ( weight2 ) : raise ValueError  ( "weight2 must be *trivial*" )
            weight1 = None
            weight2 = None

        ## unpack data 
        uds1 , uds2 = self.unpack ( data1 , data2 ) 
        
        shape1 = uds1.shape 
        shape2 = uds2.shape 
        if 1 == len ( shape1 ) : uds1 = uds1.reshape ( -1 , shape1 [ 0 ] ) 
        if 1 == len ( shape2 ) : uds2 = uds2.reshape ( -1 , shape2 [ 0 ] ) 
        
        ## check validity and consitency of input parameters 
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) )

        ## normalize
        if normalize and self.normalize :
            uds1 , uds2 = self.normalize_pooled ( uds1 , uds2 ) 
            
        n1 = len ( uds1 ) 
        n2 = len ( uds2 )
        
        ## calculate sums of distances, Eq (3.7) 
        result  = self.sum_distances ( uds1 , uds1 ) / ( n1 * ( n1 - 1 ) )

        result -= self.sum_distances ( uds1 , uds2 ) / ( n1 *   n2       )

        ## add the distances from the second dataset? 
        if self.mc2mc : result += self.sum_distances ( uds2 , uds2 ) / ( n2 * ( n2 - 1 ) )
        
        ## 
        result = float ( result )
        self.t_value = result
        ## 

        return result

    # =========================================================================    
    ## Calculate t-value for two (structured) datasets 
    #  @code
    #  ppd   = ...
    #  data1 = ... ## the first  data set 
    #  data2 = ... ## the second data set
    #  t = ppd ( data1 , data2 , normalize = False ) 
    #  @endcode
    def __call__ ( self     , 
                  data1     , 
                  data2     , * , 
                  weight1   = None , 
                  weight2   = None , 
                  normalize = True ) :
        """ Calculate T-value for two data sets 
        >>> ppd   = ...
        >>> data1 = ... ## the first  data set 
        >>> data2 = ... ## the second data set
        >>> t = ppd ( data1 , data1 , normalize = False ) 
        >>> t = ppd ( data1 , data1 , normalize = True  ) 
        """        
        ## unpack data is if needed 
        uds1 , uds2 = self.unpack ( data1 , data2 ) 
        ## 
        shape1 = uds1.shape 
        shape2 = uds2.shape 
        if 1 == len ( shape1 ) : uds1 = uds1.reshape ( -1 , shape1 [ 0 ] ) 
        if 1 == len ( shape2 ) : uds2 = uds2.reshape ( -1 , shape2 [ 0 ] ) 
        ##
        
        ## check validity and consitency of input parameters 
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) )
        
        ## normalize if requested 
        if normalize and self.normalize :
            uds1 , uds2 = self.normalize_pooled ( uds1 , uds2 )            
        ##        
        return self.tvalue ( uds1      ,
                             uds2      ,
                             weight1   = weight1 ,
                             weight2   = weight2 ,
                             normalize = False   )
    
# =============================================================================
## @class DNNnp
#  Distance-to-Nearest-Neighour GoF-method 
#  @see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
#  @see https://doi.org/10.1088/1748-0221/5/09/P09004
#  @see http://arxiv.org/abs/arXiv:1003.1768
#
#  M.Williams writes: 
#    The method <..> is easy to use, requires very little processing time and is
#    conceptually fairly easy to understand; however it is not very powerful.
#    The U-statistic it defines does provide a useful easy-to-visualize diagnostic
#    tool (especially for very high dimensional analyses), but its quantitative
#    usefulness as a g.o.f. test is limited.  
class DNNnp(GoFnp) : 
    """ Distance-to-Nearest-Neighor GoF-method 
    - see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
    - see https://doi.org/10.1088/1748-0221/5/09/P09004    
    - see http://arxiv.org/abs/arXiv:1003.1768 

    M.Williams writes: 
    ... The method <..> is easy to use, requires very little processing time and is
    ... conceptually fairly easy to understand; however it is not very powerful.
    ... The U-statistic it defines does provide a useful easy-to-visualize diagnostic
    ... tool (especially for very high dimensional analyses), but its quantitative
    ... usefulness as a g.o.f. test is limited.  
    """
    def __init__ ( self            ,
                   histo    = None ,
                   nToys    = 1000 ,
                   parallel = True , **params ) :

        ## switch off parallel processing 
        if parallel and not run_parallel ( parallel ) :
            logger.info ( '%s: (internal) parallel processing is OFF' % typename ( self ) )
            parallel = False
            
        params [ 'n_jobs'    ] = 1 if parallel else num_jobs ( params , numcpu() - 1 )
        params [ 'normalize' ] = True 
        
        if 'metric' in params : params.pop ( 'metric' )
        if 'p'      in params : params.pop ( 'p'      )
        
        self.__histo = None 
        if   isinstance ( histo , ROOT.TH1 ) and 1 == histo.GetDimension () :
            self.__histo = histo
        elif isinstance ( histo , int      ) and 1 < histo :
            self.__histo = ROOT.TH1D ( hID () , 'U-values' , histo , -0.05 , 1.05 ) 

        GoFnp.__init__ ( self                     ,
                         nToys       = nToys      ,
                         parallel    = parallel   , 
                         method      = method_DNN , **params )
        
    # ==================================================================================
    @property
    def config ( self ) :
        """`config` : get all (configuration) parameters"""
        conf = super().config 
        conf [ 'histo' ] = self.__histo
        return conf 
            
    # =========================================================================
    ## Are weights supported by this estimator?
    @property
    def weights_supported ( self ) :
        """`weights_supported` : Are weights supported by this estimator?
        """
        return False

    # =========================================================================
    ## Good for two-samples comparison?
    #  Can this estimator be used for comparison of two samples?
    @property 
    def two_samples ( self ) :
        """`two_samples`: Can this estimator be used for comparison of two samples?
        """
        return False 
    
    @property
    def histo ( self ) :
        """`histo` : the histogram with distribution of U-values"""
        return self.__histo
    
    # =========================================================================
    ## Calculate the t-value
    #  @see Eqs. (3.16) in M.Williams' paper
    #  @param data1 actual data (as unstructured array)
    #  @param vpdf  array of PDF values  
    def tvalue ( self      ,
                 data      ,
                 vpdf      , * ,
                 weight1   = None ,
                 weight2   = None ,
                 normalize = True ) :
        """ Calculate the t-value
        - see Eqs. (3.14)&(3.15) in M/.Williams' paper
        data1 : actual data (as unstructured array)
        vpdf  : array of PDF values  
        """
        ## 
        if not weight_trivial ( weight2 ) : raise TypeError ( "DNNnp: weight2 must be *trivial*" )

        w1_trivial = weight_trivial ( weight1 ) 
        if not w1_trivial and not self.weight_supported :
            raise TypeError ( "DNNnp: weight1 is provides but not supported!" )
        
        ## unpack if needed 
        uds1 , uds2 = self.unpack ( data , vpdf )
        
        ## reshape it if needed 
        shape1 = uds1.shape 
        shape2 = uds2 .shape
        if 1 == len ( shape1 ) : uds1 = uds1.reshape ( -1 , shape1 [ 0 ] ) 
        
        if not valid_data_shape  ( uds1 ) : raise TypeError ( "DNNnp: invalid uds1 shape!"     )
        if not valid_data_shape  ( uds2 ) : raise TypeError ( "DNNnp: invalid uds2 shape!"     )
        if not 1 == num_features ( uds2 ) : raise TypeError ( "DNNnp: invalid #features(uds2)" )
        if num_samples ( uds1 ) != num_samples ( uds2 ) :
            raise TypeError ( "DNNnp: invalid #samples!" )
        
        uds2 = uds2.ravel() 
        
        ## # of points & dimensionality of the problem
        N , D = shape1
                
        ## normalize
        jacobian = 1.0  
        if normalize and self.normalize :
            jacobian = numpy.prod ( numpy.std  ( uds1 , axis = 0 , keepdims = True ) ) 
            uds1     = normalize_pooled ( uds1  )     

        """
        from sklearn.neighbors import NearestNeighbors   
        nn = NearestNeighbors ( **self.params )        
        nn.fit ( uds1 )        
        distances ,  _  = nn.kneighbors( uds1 )
        distances       = distances [ : , 1]  # DNN (Distance to Nearest Neighbor)
        """
        
        ## distances = nearest_neighbors ( uds1 , **self.params )
        
        distances = nearest_distances ( uds1 , **self.params ) 
        
        if  1 != D : distances = distances ** D
    
        ## volume of the ball in D-dimensions 
        VD = 1.0 * Ostap.Math.NBallVolume_ [ D ].unit_volume 
        
        ## total weight 
        WT = 1.0 * N if w1_trivial else numpy.sum ( weight1 )
        
        ## Collect all multiplicative factors 
        factor = - WT * VD * jacobian
        
        ## get u-values 
        ## expected weight in sphere
        uvalues = factor * distances 
        if not w1_trivial : uvalues /= weight1 
        uvalues = 1.0 - numpy.exp ( uvalues )

        delta   = 1.e-7
        uvalues = numpy.clip ( uvalues , delta , 1.0 - delta )

        ## fill the histogram of u-values (if defined)
        if self.__histo :
            self.__histo.Reset ()
            for u in uvalues : self.__histo.Fill ( u ) 

        # =======================================================================
        ## t-value as Gemini AI suggests: (modified Anderson-Darling criteria)
        result = - numpy.mean ( numpy.log ( uvalues ) + numpy.log ( 1.0 - uvalues ) )
        result = float ( result )
        self.t_value = result
        return result
    
    # ===========================================================================
    ## Calculate the t-value
    #  @see Eqs. (3.16) in M.Williams' paper
    #  @param data1 actual data (as structured array)
    #  @param vpdf  array of PDF values  
    def __call__ ( self      ,
                   data1     ,
                   vpdf      , * ,
                   weight1   = None  ,
                   weight2   = None  , 
                   normalize = True  ) :
        """" Calculate the t-value
        - see Eqs. (3.16) in M.Williams' paper
        data1: actual data (as structured array)
        vpdf : array of PDF values  
        """
        ## 
        if not weight_trivial ( weight2 ) : raise TypeError ( "DNNnp: weight2 must be *trivial*" )        
        ## 
        uds1 , uds2 = self.unpack ( data1 , vpdf ) 
        ## 
        return self.tvalue ( uds1               ,
                             uds2                ,
                             weight1   = weight1 ,
                             weight2   = weight2 ,                             
                             normalize = True    )
    
    # ============================================================================
    ## p-value is not really defined here 
    # 
    #  M.Williams writes:
    #  `Because of this I do not think p-value are worth calculating`
    # 
    #  - However, one always can run straightforward pseudoexperiments 
    def pvalue ( self , *args , **kwargs ) :
        """ p-value is not defined..
        
        M.Williams writes:
        ... `Because of this I do not think p-value are worth calculating`
        
        However, one always can run straightforward pseudoexperiments 
        """        
        raise NotImplementedError( "p-value is not defined for DNNnp!" )

# =============================================================================
## @class DistanceTest
#  Base class for Two-Samples/GoF tests based on
#  the gloabal shape parameters ("distances")
#  @attention All of them are  *VERY* crude "estimators"
class DistanceTest(GoFnp) :
    """Base class for Two-Sampels/GoF tests based on
    global shape parameters ("distances")
    - attention All of them are  *VERY* crude "estimators"
    """    
    def __init__ ( self        , *     , 
                   nToys       = 1000  ,
                   parallel    = True  , 
                   method      = "<UNSPECIFIED>" , 
                   check_vct   = True  ,
                   normalize   = True  , **params ) :         
        
        ## check data-vectors ?
        self.__check_vct    = True if check_vct else False
        self.__invalid_data = 0 
        ## initialize the base 
        GoFnp.__init__ ( self                            , 
                         nToys        = nToys            ,
                         parallel     = parallel         , ## parallel is true 
                         method       = method           ,
                         normalize    = True             , ## normailzation is FORCED
                         check_vct    = self.__check_vct , **params )
        
    # =========================================================================
    @property
    def check_vct ( self ) :
        """`check_vct` : check data vectors ? """
        return self.__check_vct

    # =========================================================================
    @property
    def invalid_data ( self ) :
        """`invaild_data` : number of cases where data-vector is invalid
        """
        return self.__invalid_data

    # ==================================================================================
    @property
    def config ( self ) :
        """`config` : get all configuration parameters"""
        conf = {} 
        conf.update ( super().config )
        if self.check_vct : conf [ 'invalid_data'] = self.invalid_data
        return conf
    
    # =========================================================================
    ## Are weights supported by this estimator?
    @property
    def weights_supported ( self ) :
        """`weights_supported` : Are weights supported by this estimator?
        """
        return True
    
    # =========================================================================
    ## Good for two-samples comparison?
    #  Can this estimator be used for comparison of two samples?
    @property 
    def two_samples ( self ) :
        """`two_samples`: Can this estimator be used for comparison of two samples?
        """
        return True 

    # =========================================================================
    ## Convert numpy-array statistics into Ostap::SVectorWithError
    #  @see Ostap::Math::SVectorWithError
    def np2vstat ( self , data , weight = None ) :
        """ Convert numpy-array statistics into `Ostap.Math.SVectorWithError`
        - see `Ostap.Math.SVectorWithError`
        """
        vct = np2vct ( data , weight = weight )
        if self.check_vct and not vct.valid () :
            self.__invalid_data += 1 
            if not self.silent :
                logger.warning ( '%s: data-vector is not valid #%d' % ( typename ( self ) , self.__invalid_data ) ) 
        return vct 
                                    
# ============================================================================
## @class KullbackLeibler 
#  Use (asymmetric) KullbackLeibler divergency to discriminiate the dataset
#  @attention it is *VERY* crude "estimator"
class KullbackLeibler(DistanceTest) :
    """ Use (asymmetric) Kullback-Leibler divergency to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self  , **params ) :         
        
        ## initialize the base 
        super() .__init__ ( method = method_KL , **params )

    # =========================================================================
    # calculate t-value for (non-structured) 2D arrays
    def tvalue ( self      , 
                 data1     , 
                 data2     , *    , 
                 weight1   = None , 
                 weight2   = None ,
                 normalize = True ) :
        """ Calculate t-value for (non-structured) 2D arrays
        """
        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## check everything
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) ) 

        ## normalize
        if normalize and self.normalize :
            uds1, uds2  = self.normalize_pooled ( uds1 , uds2  ) 
            
        v1 = self.np2vstat ( uds1 , weight1 )
        v2 = self.np2vstat ( uds2 , weight2 )
        
        return v1.kullback_leibler ( v2 ) 

# ============================================================================
## @class Jeffrey
#  Use Jeffrey' (aka symmetric Kullback-Leibler divergency  to discriminiate the datasets 
#  @attention it is *VERY* crude "estimator"
class Jeffrey(DistanceTest) :
    """ Use Jeffrey's (aka symmetric Kullback-Leibler) divergency to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self , **params ) :         
        
        ## initialize the base 
        super() .__init__ ( method = method_J , **params )

    # =========================================================================
    # calculate t-value for (non-structured) 2D arrays
    def tvalue ( self      , 
                 data1     , 
                 data2     , *    , 
                 weight1   = None , 
                 weight2   = None ,
                 normalize = True ) :
        """ Calculate t-value for (non-structured) 2D arrays
        """
        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## check everything
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) ) 

        ## normalize
        if normalize and self.normalize :
            uds1, uds2  = self.normalize_pooled ( uds1 , uds2  ) 
            
        v1 = self.np2vstat ( uds1 , weight1 )
        v2 = self.np2vstat ( uds2 , weight2 )
        
        return v1.jeffrey ( v2 ) 

# ============================================================================
## @class JensenShannon
#  Use Jensen-Shannon', another form of  symmetric Kullback-Leibler divergency  to discriminiate the datasets 
#  @attention it is *VERY* crude "estimator"
class JensenShannon(DistanceTest) :
    """ Use Jensen-Shannon's - another form of symmetric Kullback-Leibler) divergency to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self , **params ) :         
        
        ## initialize the base 
        super() .__init__ ( method = method_JS , **params )

    # =========================================================================
    # calculate t-value for (non-structured) 2D arrays
    def tvalue ( self      , 
                 data1     , 
                 data2     , *    , 
                 weight1   = None , 
                 weight2   = None ,
                 normalize = True ) :
        """ Calculate t-value for (non-structured) 2D arrays
        """
        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## check everything
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) ) 

        ## normalize
        if normalize and self.normalize :
            uds1, uds2  = self.normalize_pooled ( uds1 , uds2  ) 
            
        v1 = self.np2vstat ( uds1 , weight1 )
        v2 = self.np2vstat ( uds2 , weight2 )

        n1 = num_samples ( uds1 ) if weight_trivial ( weight1 ) else float ( numpy.sum ( weight1 ) ) 
        n2 = num_samples ( uds2 ) if weight_trivial ( weight2 ) else float ( numpy.sum ( weight2 ) ) 
        
        return v1.jensen_shannon ( v2 , n1 , n2 ) 
    
# =============================================================================
## @class Mahalanobis
#  Use Mahalanobis distance to discriminiate the dataset
#  @attention it is *VERY* crude "estimator"
class Mahalanobis(DistanceTest) :
    """ Use Mahalanobis distance to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self , **params ) :         
        
        ## initialize the base 
        super().__init__ ( method = method_M , **params )

    # =========================================================================
    # calculate t-value for (non-structured) 2D arrays
    def tvalue ( self      , 
                 data1     , 
                 data2     , *    , 
                 weight1   = None , 
                 weight2   = None ,
                 normalize = True ) :
        """ Calculate t-value for (non-structured) 2D arrays
        """
        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## check everything
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) ) 
        
        ## normalize
        if normalize and self.normalize :
            uds1, uds2  = self.normalize_pooled ( uds1 , uds2  ) 
            
        v1 = self.np2vstat ( uds1 , weight1 )
        v2 = self.np2vstat ( uds2 , weight2 )
        
        n1 = num_samples ( uds1 ) if weight_trivial ( weight1 ) else float ( numpy.sum ( weight1 ) ) 
        n2 = num_samples ( uds2 ) if weight_trivial ( weight2 ) else float ( numpy.sum ( weight2 ) ) 
        
        return v1.mahalanobis ( v2 , n1 , n2 )
    
# ============================================================================
## @class Hotelling  
#  Use Hotelling's t-squared statistics to discriminiate the datasets
#  @attention it is *VERY* crude "estimator"
class Hotelling(DistanceTest) :
    """ Use Hotelling's t-squared statistics to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self        , **params ) : 

        ## initialize the base 
        super() .__init__ ( method = method_T2 , **params )
        
    # =========================================================================
    # calculate t-value for (non-structured) 2D arrays
    def tvalue ( self      , 
                 data1     , 
                 data2     , *    , 
                 weight1   = None , 
                 weight2   = None ,
                 normalize = True ) :
        """ Calculate t-value for (non-structured) 2D arrays
        """
        ##
        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## check everything
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) ) 
        
        ## normalize
        if normalize and self.normalize :
            uds1, uds2  = self.normalize_pooled ( uds1 , uds2  ) 

        v1 = self.np2vstat ( uds1 , weight1 )
        v2 = self.np2vstat ( uds2 , weight2 )
        
        n1 = num_samples ( uds1 ) if weight_trivial ( weight1 ) else float ( numpy.sum ( weight1 ) ) 
        n2 = num_samples ( uds2 ) if weight_trivial ( weight2 ) else float ( numpy.sum ( weight2 ) ) 
                        
        return v1.hotelling ( v2 , n1 , n2 ) 

# ============================================================================
## @class Bhattacharyya 
#  Use Bhattacharyya' statistics to discriminiate the datasets
#  @attention it is *VERY* crude "estimator"
class Bhattacharyya(DistanceTest) :
    """ Use Bhatatcharyya' statistics to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self        , **params ) :
        ## initialize the base 
        super() .__init__ ( method = method_B  , **params )
        
    # =========================================================================
    # calculate t-value for (non-structured) 2D arrays
    def tvalue ( self      , 
                 data1     , 
                 data2     , *    , 
                 weight1   = None , 
                 weight2   = None ,
                 normalize = True ) :
        """ Calculate t-value for (non-structured) 2D arrays
        """
        ##
        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## check everything
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) ) 
        
        ## normalize
        if normalize and self.normalize :
            uds1, uds2  = self.normalize_pooled ( uds1 , uds2  ) 

        v1 = self.np2vstat ( uds1 , weight1 )
        v2 = self.np2vstat ( uds2 , weight2 )
        
        n1 = num_samples ( uds1 ) if weight_trivial ( weight1 ) else float ( numpy.sum ( weight1 ) ) 
        n2 = num_samples ( uds2 ) if weight_trivial ( weight2 ) else float ( numpy.sum ( weight2 ) ) 
                        
        return v1.bhattacharyya ( v2 , n1 , n2 ) 

# ============================================================================
## @class Wasserstein
#  Use Wasserstein' statistics to discriminiate the datasets
#  @attention it is *VERY* crude "estimator"
class Wasserstein(DistanceTest) :
    """ Use Wassertein' statistics to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self        , **params ) :
        ## initialize the base 
        super() .__init__ ( method = method_W2 , **params )
        
    # =========================================================================
    # calculate t-value for (non-structured) 2D arrays
    def tvalue ( self      , 
                 data1     , 
                 data2     , *    , 
                 weight1   = None , 
                 weight2   = None ,
                 normalize = True ) :
        """ Calculate t-value for (non-structured) 2D arrays
        """
        ##
        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## check everything
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) ) 

        ## normalize
        if normalize and self.normalize :
            uds1, uds2  = self.normalize_pooled ( uds1 , uds2  ) 

        v1 = self.np2vstat ( uds1 , weight1 )
        v2 = self.np2vstat ( uds2 , weight2 )

        n1 = num_samples ( uds1 ) if weight_trivial ( weight1 ) else float ( numpy.sum ( weight1 ) ) 
        n2 = num_samples ( uds2 ) if weight_trivial ( weight2 ) else float ( numpy.sum ( weight2 ) ) 
                
        return v1.wasserstein ( v2 , n1 , n2 ) 
    

# ============================================================================
## @class Hellinger
#  Use Wasserstein' statistics to discriminiate the datasets
#  @attention it is *VERY* crude "estimator"
class Hellinger(DistanceTest) :
    """ Use Hellinger' statistics to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self        , **params ) :
        ## initialize the base 
        super() .__init__ ( method = method_H , **params )
        
    # =========================================================================
    # calculate t-value for (non-structured) 2D arrays
    def tvalue ( self      , 
                 data1     , 
                 data2     , *    , 
                 weight1   = None , 
                 weight2   = None ,
                 normalize = True ) :
        """ Calculate t-value for (non-structured) 2D arrays
        """
        ##
        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## check everything
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) ) 

        ## normalize
        if normalize and self.normalize :
            uds1, uds2  = self.normalize_pooled ( uds1 , uds2  ) 

        v1 = self.np2vstat ( uds1 , weight1 )
        v2 = self.np2vstat ( uds2 , weight2 )

        n1 = num_samples ( uds1 ) if weight_trivial ( weight1 ) else float ( numpy.sum ( weight1 ) ) 
        n2 = num_samples ( uds2 ) if weight_trivial ( weight2 ) else float ( numpy.sum ( weight2 ) ) 
                
        return v1.hellinger ( v2 , n1 , n2 ) 


## import Chi2 
from   ostap.stats.gofchi2      import Chi2 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
