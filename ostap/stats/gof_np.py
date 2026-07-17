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
    'DNNnp'           , ## Distance-to-Nearest-Neighbour Goodness-of-Fit method
    ##
    'Mahalanobis'     , ## Very crude estiamtor based on Mahalanobis' disatnce
    'KullbackLeibler' , ## Very crude estimator based on Kullback-Leibler's divergency 
    'Hotelling'       , ## Very crude estimator based on Hotelling's distance 
)
# =============================================================================
from   ostap.core.ostap_types   import string_types, num_types 
from   ostap.stats.gof_utils    import PERMUTATOR
from   ostap.core.core          import SE, VE, Ostap, hID  
from   ostap.stats.counters     import EffCounter
from   ostap.utils.progress_bar import progress_bar
from   ostap.utils.utils        import split_n_range
from   ostap.utils.basic        import numcpu, typename  
from   ostap.stats.gof          import AGoFnp
from   ostap.stats.gof_utils    import ( run_parallel       ,
                                         num_jobs           , 
                                         weight_trivial     ,
                                         normalize_pooled   ,
                                         pairwise_distances ,
                                         nearest_distances  ,
                                         nearest_neighbors  , 
                                         draw_ecdf          , s2u ) 
from   ostap.utils.memory       import memory, memory_enough
from   ostap.math.math_ve       import gauss_cdf 
import ostap.math.math_base           
import ROOT, os, abc, numpy, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gof_np' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-fit studies for multidimensional fits' )
# =============================================================================
## @class GoFnp 
#  A base class for numpy-related family of methods to probe goodness-of-fit
class GoFnp (AGoFnp) :
    """ A base class for numpy-related family of methods to probe goodness-of-fit
    """
    def __init__ ( self               ,
                   nToys     = 0      ,
                   silent    = False  , 
                   parallel  = False  ,
                   method    = 'GoF'  ,
                   progress  = True   ,
                   normalize = True   , **params ) : 

        assert isinstance ( nToys , int ) and 0 <= nToys  , \
            "Invalid number of permulations/toys:%s" % nToys
        
        self.__nToys    = nToys
        ## 
        self.__silent    = True if silent    else False
        self.__parallel  = True if parallel  else False
        self.__progress  = True if progress  else False
        self.__normalize = True if normalize else False 
        ## 
        self.__method    = method
        self.__params    = params

        ## Empirical CDF for t-value distribution from permutations/toys"""
        self.__ecdf      = None
        self.__counter   = None
        self.__tvalue    = None
        self.__pvalue    = None
                
        if self.__parallel :
            mratio = memory_enough () / numcpu () 
            if mratio < 1 :
                logger.warning ( 'Available/Used memory ratio: %.1f; switch-off parallel processing' % mratio )                
                ## self.__parallel = False
                
    @property
    def params ( self ) :
        """`params` : configuration of underlying classifier"""
        return self.__params

    @property
    def normalize ( self ) :
        """`normalize` : scale and shift both datasets to have mean = 0 and rms=1 for each column of pooled data"""
        return self.__normalize
    
    # ==================================================================================
    @property
    def config ( self ) :
        """`config` : get all configuration parameters"""
        conf = {} 
        conf.update ( self.params )
        conf [ 'nToys'            ] = self.nToys
        conf [ 'silent'           ] = self.silent
        conf [ 'progress'         ] = self.progress
        conf [ 'parallel'         ] = self.parallel
        conf [ 'normalize'        ] = self.normalize
        conf [ 'method'           ] = self.method  
        conf [ 'weight_supported' ] = self.weights_supported
        return conf 
    
    # =========================================================================
    ## self-print get the configuration 
    def table (  self , prefix = '# ') : 
        """ print configuration """
        from ostap.logger.utils import map2table_ex
        title = "%s configuration " % typename ( self )
        return map2table_ex ( self.config , 
                              header      = ( 'Parameter' , 'type' , 'value' ) ,
                              ailgnment   = 'rcw'  , 
                              prefix      = prefix ,
                              title       = title  )
    
    def __str__  ( self ) : return self.table ( prefix = '' ) 
    def __repr__ ( self ) : return self.__str__ ()
    
    # =========================================================================
    @property 
    def nToys ( self ) :
        """`nToys` : number of permutations/toys used for permutation/toys test"""
        return self.__nToys
    # =========================================================================
    @property
    def silent  ( self ) :
        """`silent` : silent processing?"""
        return self.__silent
    # ========================================================================
    @property
    def parallel ( self ) :
        """`parallel` : parallel processing where/when/if possible?"""
        return self.__parallel
    # ========================================================================
    @property
    def progress ( self ) :
        """`progress` : show progress bar?"""
        return self.__progress 
    # ========================================================================
    @property
    def method ( self ) :
        """`method` : the actual GoF method """
        return self.__method
    
    # =======================================================================
    ##  Unpack data ( consvert from structured to unstructured arrays)
    #   @code
    #   gof   = 
    #   ds1   , ds2  = ...
    #   data1 , dat2 = gof.unpack ( ds1 , ds2 ) 
    #   @endcode 
    def unpack ( self , ds1 , ds2 ) :
        """ Unpack data ( convert from structured to unstructured arrays)
        >>> gof   = 
        >>> ds1   , ds2  = ...
        >>> data1 , dat2 = gof.unpack ( ds1 , ds2 ) 
        """
        ## 
        ## transform ?  
        structured1 = True if ds1.dtype.fields else False
        structured2 = True if ds2.dtype.fields else False
        ##
        ## convert to unstructured datasets 
        data1 = s2u ( ds1 , copy = False ) if structured1 else ds1
        data2 = s2u ( ds2 , copy = False ) if structured2 else ds2
        ##
        return data1 , data2 
                       
    # =======================================================================
    ## Calculate T-value for two (structured) datasets 
    #  @code
    #  adval  = ...
    #  data1 = ... ## the first  data set 
    #  data2 = ... ## the second data set
    #  t = adval ( data1 , data1 , normalize = False ) 
    #  t = adval ( data1 , data1 , normalize = True  ) 
    #  @endcode
    def __call__ ( self              ,
                   data1             ,
                   data2             , * ,
                   weight1   = None  ,
                   weight2   = None  ,
                   normalize = True  ) :
        
        """ Calculate T-value for two (STRUCTURED) data sets 
        >>> adval = ...
        >>> data1 = ... ## the first  data set 
        >>> data2 = ... ## the second data set 
        >>> t = adval ( data1 , data1 , normalize = False ) 
        >>> t = adval ( data1 , data1 , normalize = True  ) 
        """        
        
        if not self.weights_supported :
            assert weight_trivial ( weight1 ) , "weight1 must be *trivial*"
            assert weight_trivial ( weight2 ) , "weight2 must be *trivial*"
            weight1 = None
            weight2 = None

        ## transform ?  
        uds1 , uds2 = self.unpack ( data1 , data2 ) 
            
        ## normalize
        if normalize and self.normalize : uds1 , uds2 = normalize_pooled ( uds1 , uds2 ) 
        
        return self.tvalue ( uds1      ,
                             uds2      ,
                             weight1   = weight1 ,
                             weight2   = weight2 ,
                             normalize = False   )
    
    # =========================================================================
    ## Calculate the t & p-values
    #  @code
    #  gof = ...
    #  data1 , data2 = ...
    #  t , p = gof.pvalue ( data1 , data2 , normalize = False ) 
    #  @endcode 
    def pvalue ( self             , 
                 data1            ,
                 data2            , * ,
                 tvalue    = None , 
                 weight1   = None ,
                 weight2   = None ) : 
                
        """ Calculate the t & p-values
        >>> gof  = ...
        >>> data1 , data2 = ...
        >>> t   , p = gof.pvalue ( ds1 , ds2 , normalize = True ) 
        """
        
        if not self.weights_supported :
            assert weight_trivial ( weight1 ) , "weight1 must be *trivial*"
            assert weight_trivial ( weight2 ) , "weight2 must be *trivial*"
            weight1 = None
            weight2 = None

        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 ) 
        
        ## normalize ? 
        if self.normalize : uds1 , uds2 = normalize_pooled ( uds1 , uds2 ) 

        ### calculate t-value
        t_value    = tvalue if not tvalue is None else self.tvalue ( uds1      ,
                                                                     uds2      ,
                                                                     weight1   = weight1 ,
                                                                     weight2   = weight2 ,
                                                                    normalize = False   )
        
        ## use permutations to get the p-value 
        permutator = PERMUTATOR ( self              ,
                                  t_value           ,
                                  uds1              ,
                                  uds2              ,
                                  weight1 = weight1 ,
                                  weight2 = weight2 )
        
        if self.parallel and permutator.run : counter = permutator.run ( self.nToys , progress = self.progress )            
        else                                : counter = permutator     ( self.nToys , progress = self.progress )

        # ==================================================
        ## @see Phipson, Belinda; Smyth, Gordon K (2010).
        #       "Permutation p-values should never be zero:
        #       calculating exact p-values when permutations are randomly drawn".
        #       Statistical Applications in Genetics and Molecular Biology. 9 (1) 39.
        #       arXiv:1603.05766. doi:10.2202/1544-6115.1585. PMID 21044043. S2CID 10735784.
        # counter += True

        ## get the t-value distribution from permutator: ECDF & COUNTER 
                
        ## get the efficiency/p-value from the counter
        p_value      = counter.eff

        self.ecdf    = permutator.ecdf
        self.counter = counter
        self.t_value = t_value
        self.p_value = p_value 
        
        return t_value , p_value
    
    @property
    def ecdf ( self ) :
        """`ecdf` : empirical CDF for t-value distribution from permutations"""
        return self.__ecdf
    @ecdf.setter
    def ecdf ( self , value ) :
        assert value is None or isinstance ( value , Ostap.Math.ECDF ) , \
            "Invalid type for ECDF: %s" % typename ( value )
        self.__ecdf = value 

    @property
    def counter ( self ) :
        """`counter` : get the efficiency counter from toys"""
        return self.__counter
    @counter.setter
    def counter ( self , value ) :
        assert value is None or isinstance ( value , EffCounter ) , \
            "Invalid counter type %s" % typename ( value ) 
        self.__counter = value
        
    ## access the calculated t-value 
    @property 
    def t_value ( self ) :
        """`t_value` : get the calculated t-value """
        return self.__tvalue
    @t_value.setter
    def t_value ( self , value ) :
        assert value is None or isinstance ( value , num_types ) , \
            "Invalid t-value type %s" % typename ( value ) 
        self.__tvalue = value

    # ========================================================================
    ## access the calculated p-value 
    @property 
    def p_value ( self ) :
        """`p_value` : get the calculated p-value """
        return self.__pvalue
    @p_value.setter
    def p_value ( self , value ) :
        assert value is None or isinstance ( value , num_types ) or isinstance ( value , VE ) , \
            "Invalid p-value type %s" % typename ( value ) 
        self.__pvalue = value

    # =========================================================================
    ## Get results in a form of the table 
    def report ( self           ,
                 tvalue  = None ,
                 pvalue  = None ,
                 ecdf    = None ,
                 counter = None ,
                 title   = ''   ,
                 prefix  = ''   ,
                 style   = None ) :
        """ Get results in a for of the table 
        """
        return super().report ( tvalue  = tvalue  if not tvalue  is None else self.__tvalue  ,
                                pvalue  = pvalue  if not pvalue  is None else self.__pvalue  ,
                                ecdf    = ecdf    if not ecdf    is None else self.__ecdf    ,
                                counter = counter if not counter is None else self.__counter ,
                                title   = title  if title else '%s GoF-report' % typename ( self ) , 
                                prefix  = prefix ,
                                style   = style  )

    # ========================================================================
    ## Get results in form of the row in the table
    #  @code
    #  gof = ...
    #  header , row = gof.the_row ( ... ) 
    #  @endcode `
    def the_row ( self             ,
                  tvalue    = None ,
                  pvalue    = None ,
                  ecdf      = None ,
                  counter   = None ,                  
                  precision = 4    ,
                  width     = 6    ) :         
        """ Get results in form of the table 
        >>> gof = ...
        >>> header , row = gof.the_row ( ... ) 
        """
        return super().the_row ( tvalue  = tvalue  if not tvalue  is None else self.__tvalue  ,
                                 pvalue  = pvalue  if not pvalue  is None else self.__pvalue  ,
                                 ecdf    = ecdf    if not ecdf    is None else self.__ecdf    ,
                                 counter = counter if not counter is None else self.__counter ,
                                 title   = title  if title else '%s GoF-report' % typename ( self ) , 
                                 prefix  = prefix ,
                                 style   = style  )
    

    # =========================================================================
    ## Draw the empirical CDF from permutations or toys  
    def draw  ( self , option = '' , * , tvalue = None , **kwargs ) :
        """ Draw empirical CDF from permutations or toys 
        """
        ## 
    
        ecdf = self.ecdf 
        if not ecdf : return ecdf 
        ## 
        tvalue     = self.t_value if tvalue is None else tvalue 
        has_tvalue = isinstance ( tvalue , num_types ) 
        ##
        if not has_tvalue : return draw_ecdf (  ecdf , tvalue = None   , option = option , **kwargs )
        result , vline , hline =   draw_ecdf (  ecdf , tvalue = tvalue , option = option , **kwargs )
        ## 
        self._vline = vline 
        self._hline = hline 
        ##
        return result, vline, hline   
    
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
                   parallel    = False , 
                   silent      = False ,
                   progress    = True  ,
                   n_neighbors = 10    , **params ) : 
        
        # =================================================================================
        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "Parallel processing is switched OFF!" ) 
            parallel = False 
            
        ## Attention!
        assert isinstance ( n_neighbors , int ) and 2 <= n_neighbors , \
            "Invalid `n_neighbors`: %s" % n_neighbors 

        ## store it
        self._k_max = n_neighbors 
        ## 
        n_jobs = num_jobs ( params , numcpu() - 1 )
        if parallel : n_jobs = 1 

        ## initialize the base 
        GoFnp.__init__ ( self                          , 
                         nToys        = nToys          ,
                         parallel     = parallel       , 
                         silent       = silent         ,
                         progress     = progress       ,                         
                         method       = 'Mixed Sample' ,
                         normalize    = True           , 
                         n_neighbors  = self.k_max     ,
                         n_jobs       = n_jobs         , **params )

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

        shape1 = data1.shape
        shape2 = data2.shape
        assert 2 == len ( shape1 ) and 2 == len ( shape2 ) and shape1 [ 1 ]  == shape2 [ 1 ] , \
            "Invalid arrays: %s , %s" % ( shape1 , shape2  )

        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## normalize
        if normalize and self.normalize :
            uds1, uds2  = normalize_pooled ( uds1 , uds2  ) 
            
        ## 
        n1 = len ( uds1 ) 
        n2 = len ( uds2 ) 
                
        data       = numpy.vstack       ( [ uds1  , uds2 ]    )
        labels     = numpy.array        ( [ 1 ] * n1 + [ 0 ] * n2  )

        ## weights 
        w1_trivial = weight_trivial ( weight1 )
        w2_trivial = weight_trivial ( weight2 )
        
        ## combine weights, if needed 
        if   w1_trivial and w2_trivial : weights = None
        else : 
            if weight1 is None : weight1 = numpy.ones ( n1 )
            if weight2 is None : weight2 = numpy.ones ( n2 )            
            weights = numpy.concatenate ( [ weight1 , weight2 ] )

        ##

        
        ## from sklearn.neighbors import NearestNeighbors
        ## nn = NearestNeighbors ( **self.params )
        ## nn.fit ( data )
        ## _  , indices      = nn.kneighbors ( data )        
        ## actual_neighbors = indices[ : , 1: ]
        
        actual_neighbors = nearest_neighbors ( data , **self.params )
        
        
        source_labels    = labels [ : , numpy.newaxis ] # (N, 1)
        neighbor_labels  = labels [ actual_neighbors  ] # (N, K)

        # I(i, k) = 1 
        
        I_ik = ( source_labels == neighbor_labels ) . astype ( int )
        
        if weights is None :
            result = numpy.sum ( I_ik ) / ( 1.0 * self.k_max * ( n1 + n2 ) )
            return float ( result ) 

        w_i = weights [ :, numpy.newaxis ]  ## (N, 1)
        w_k = weights [ actual_neighbors ]  ## (N, K)
        
        pair_weights = w_i * w_k            ## (N, K)
        
        weighted_numerator = numpy.sum ( I_ik * pair_weights)
        total_weight_sum   = numpy.sum ( pair_weights )
        
        result = weighted_numerator / total_weight_sum

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
    
    """
    def __init__ ( self                   ,
                   mc2mc     = False      ,
                   nToys     = 1000       ,
                   psi       = 'gaussian' ,
                   sigma     = 0.10       ,
                   parallel  = False      , 
                   silent    = False      ,
                   progress  = True       , 
                   maxsize   = 1000000    , **params ) :


        n_jobs = num_jobs ( params , numcpu() - 1 )
        if parallel : n_jobs = 1 
        
        GoFnp.__init__ ( self                 ,
                         nToys     = nToys    ,
                         parallel  = parallel , 
                         silent    = silent   ,
                         progress  = progress , ## ATTENTION!                          
                         normalize = True     ,
                         n_jobs    = n_jobs   , 
                         method    = 'Point-to-Point Dissimilarity' , **params )
        
        self.__mc2mc     = True if mc2mc else False
        self.__transform = None
        self.__sigma     = sigma
        self.__psi       = psi
        assert isinstance ( maxsize , int ) and 0 < maxsize , "Invalid `maxsize' : %s" % maxsize
        
        self.__maxsize   = max ( maxsize , 100000  )
        
        ## check validity of `psi`
        scale = -0.5 / ( self.sigma ** 2 ) 
        self.__distance_type , _ , _ = psi_conf ( psi , scale )

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
            assert weight_trivial ( weight1 ) , "weight1 must be *trivial*"
            assert weight_trivial ( weight2 ) , "weight2 must be *trivial*"
            weight1 = None
            weight2 = None

        ## unpack data 
        uds1 , uds2 = self.unpack ( data1 , data2 ) 

        shape1 = uds1.shape 
        shape2 = uds2.shape 
        if 1 == len ( shape1 ) : uds1.reshape ( -1 , shape1 [ 0 ] ) 
        if 1 == len ( shape2 ) : uds2.reshape ( -1 , shape2 [ 0 ] ) 
        
        ## normalize
        if normalize and self.normalize :
            uds1 , uds2 = normalize_pooled ( uds1 , uds2 ) 
            
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
        if 1 == len ( shape1 ) : uds1.reshape ( -1 , shape1 [ 0 ] ) 
        if 1 == len ( shape2 ) : uds2.reshape ( -1 , shape2 [ 0 ] ) 
        ## 
        ## normalize
        if normalize and self.normalize : uds1 , uds2 = normalize_pooled ( uds1 , uds2 )            
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
    def __init__ ( self              ,
                   histo    = None   ,
                   nToys    = 1000   ,
                   parallel = False  , 
                   silent   = False  ,
                   progress = True   , **params ) :
        
        n_jobs = num_jobs ( params , numcpu() - 1 )
        if parallel : n_jobs = 1 

        ## 
        if 'metric' in params : params.pop ( 'metric' )
        if 'p'      in params : params.pop ( 'p'      )
        
        GoFnp.__init__ ( self                      ,
                         nToys       = nToys       ,
                         parallel    = parallel    , 
                         silent      = silent      ,
                         progress    = progress    ,
                         normalize   = True        ,                          
                         method      = 'Distance-to-Nearest-Neighbor' ,
                         n_jobs      = n_jobs      , **params )
        

        self.__histo = None 
        if   isinstance ( histo , ROOT.TH1 ) :
            self.__histo = histo
        elif isinstance ( histo , int      ) and 1 < histo :
            self.__histo = ROOT.TH1D ( hID() , 'U-values' , histo , -0.1 , 1.1 ) 

    # =========================================================================
    ## Are weights  are supported by this estimators?
    @property 
    def weights_supported ( self ) :
        """ Are weights are supported by this estimators?
        """
        return True 

    # =========================================================================
    ## Good for two-samples comparison?
    #  Can this estimator be used for comparison of two samples?
    @property 
    def two_samples ( self ) :
        """`two_samples`: Can this estimator be used for comparison of two samples?
        """
        return False 

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
        assert weight_trivial ( weight2 ) , "weight2 must be *trivial*"

        ## 
        w1_trivial = weight_trivial ( weight1 ) 

        ## unpack if needed 
        uds1 , uds2 = self.unpack ( data , vpdf ) 
        
        ## reshape it if needed 
        shape1 = uds1.shape 
        if 1 == len ( shape1 ) : uds1.reshape ( -1 , shape1 [ 0 ] ) 
        
        shape2 = uds2 .shape
        assert 2 == len ( shape1 ) and 1 == len ( shape2 ) and len ( uds1 ) == len ( uds2 ) , \
            "Invalid arrays: %s , %s" % ( shape1 , shape2 )

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
        
        distances = nearest_distances ( uds1 , **self.params ) 

        
        if  1 != D : distances = distances ** D
    
        ## volume of the ball in D-dimensions 
        VD = 1.0 * Ostap.Math.NBallVolume_ [ D ].unit_volume 
        
        ## total weight 
        WT = 1.0 * N if w1_trivial else numpy.sum ( weight1 )
        
        ## Collect all multiplicative factors 
        factor = - WT * VD * jacobian
        
        ## get u-values 
        ##  expected weight in sphere
        uvalues = factor * distances 
        if not w1_trivial : uvalues /= weight1 
        uvalues = 1.0 - numpy.exp ( uvalues )

        delta   = 1.e-8 
        uvalues = numpy.clip ( uvalues , delta , 1.0 - delta )

        ## fill the histogram of u-values (if defined)
        if self.__histo :
            self.__histo.Reset () 
            for u in uvalues : self.__histo.Fill ( u ) 
            
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
        assert weight_trivial ( weight2 ) , "weight2 must be trivial!"
        ## 
        uds1 , uds2 = self.unpack ( data1 , vpdf ) 
        ## 
        shape1 = uds1.shape 
        if 1 == len ( shape1 ) : uds1.reshape ( -1 , shape1 [ 0 ] ) 
        ##        
        return self.tvalue ( uds1               ,
                             uds2                ,
                             weight1   = weight1 ,
                             weight2   = weight2 ,                             
                             normalize = True    )
    
    '''
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
    '''
    @property
    def histo ( self ) :
        """`histo` : the histogram with distribution of U-values"""
        return self.__histo

# =============================================================================
## @class Mahalanobis
#  Use Mahalanobis distance to discriminiate the dataset
#  @attention it is *VERY* crude "estimator"
class Mahalanobis(GoFnp) :
    """ Use Mahalanobis distance to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self ,
                   nToys       = 400   ,
                   parallel    = False , 
                   silent      = False ,
                   progress    = True  , **params ) :         
        
        ## initialize the base 
        GoFnp.__init__ ( self                          , 
                         nToys        = nToys          ,
                         parallel     = parallel       , 
                         silent       = silent         ,
                         progress     = progress       ,                         
                         method       = 'Mahalanobis'  ,
                         normalize    = True           , **params )

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
    ## Convert numpy-array statictics into Ostap::SVectorWithError
    #  @see Ostap::Math::SVectorWithError
    def np2vstat ( self , data , weight = None ) :
        """ Convert numpy-array statictics into `Ostap.Math.SVectorWithError`
        - see `Ostap.Math.SVectorWithError`
        """
        
        shape     = data.shape
        n , N     = shape 

        w_trivial = weight_trivial ( weight )
        w         = None if w_trivial else weight 
        wsum      = n    if w_trivial else numpy.sum ( w )
        
        mean      = numpy.average ( data , axis = 0 , weights = w )
        centered  = data - mean
        weighted  = centered if w_trivial else centered * w [ : , numpy.newaxis ]
        covmtrx   = numpy.dot ( centered.T , weighted ) / wsum 

        ## prepare output 
        RT        = Ostap.Math.SVectorWithError [ N ]
        values    = RT.Value      () 
        covs      = RT.Covariance () 

        for i in range ( N  ) :
            values [ i ] = float ( mean [ i ] )
            for j in range ( i , N ) :
                cij = float ( covmtrx [ i ] [ j ]  )
                if i == j : covs [ i , j ] = cij
                else : 
                    cij = 0.5 * ( cij + float ( covmtrx [ j ] [ i ]  ) )
                    covs [ i, j ] = cij
                    covs [ j, i ] = cij 
                    
        return RT ( values ,  covs ) 
                        
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

        shape1 = data1.shape
        shape2 = data2.shape
        assert 2 == len ( shape1 ) and 2 == len ( shape2 ) and shape1 [ 1 ]  == shape2 [ 1 ] , \
            "Invalid arrays: %s , %s" % ( shape1 , shape2  )

        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## normalize
        if normalize and self.normalize : uds1, uds2  = normalize_pooled ( uds1 , uds2  ) 
            
        v1 = self.np2vstat ( uds1 , weight1 )
        v2 = self.np2vstat ( uds2 , weight2 )
        
        return v1.mahalanobis ( v2 )
    
# ============================================================================
## @class KullbackLeibler 
#  Use KullbackLeibler divergency to discriminiate the dataset
#  @attention it is *VERY* crude "estimator"
class KullbackLeibler(Mahalanobis) :
    """ Use Kullback-Leibler divergency to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self ,
                   nToys       = 400   ,
                   parallel    = False , 
                   silent      = False ,
                   progress    = True  ,
                   symmetric   = True  , **params ) :         
        
        self.__symmetric = True if symmetric else False
        
        ## initialize the base 
        GoFnp.__init__ ( self                            , 
                         nToys        = nToys            ,
                         parallel     = parallel         , 
                         silent       = silent           ,
                         progress     = progress         ,                           
                         method       = 'Kullback-Leibler/symmetric' if self.__symmetric else 'Kullback-Leibler',
                         symmetric    = self.__symmetric , 
                         normalize    = True             , **params )
                                        
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

        shape1 = data1.shape
        shape2 = data2.shape
        assert 2 == len ( shape1 ) and 2 == len ( shape2 ) and shape1 [ 1 ]  == shape2 [ 1 ] , \
            "Invalid arrays: %s , %s" % ( shape1 , shape2  )

        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## normalize
        if normalize and self.normalize :
            uds1, uds2  = normalize_pooled ( uds1 , uds2  ) 
            
        v1 = self.np2vstat ( uds1 , weight1 )
        v2 = self.np2vstat ( uds2 , weight2 )
        
        return v1.kullback_leibler ( v2 ) if self.__symmetric else v1.asymmetric_kullback_leibler ( v2 )

# ============================================================================
## @class Hotelling  
#  Use Hotelling's tsquared statistics to discriminiate the datasets
#  @attention it is *VERY* crude "estimator"
class Hotelling(Mahalanobis) :
    """ UseHotelling's t-squared statistics to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self ,
                   nToys       = 400   ,
                   parallel    = False , 
                   silent      = False ,
                   progress    = True  , **params ) :         
        
        ## initialize the base 
        GoFnp.__init__ ( self                            , 
                         nToys        = nToys            ,
                         parallel     = parallel         , 
                         silent       = silent           ,
                         progress     = progress         ,                           
                         method       = 'Hotelling'      , 
                         normalize    = True             , **params )
                                        
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

        shape1 = data1.shape
        shape2 = data2.shape
        assert 2 == len ( shape1 ) and 2 == len ( shape2 ) and shape1 [ 1 ]  == shape2 [ 1 ] , \
            "Invalid arrays: %s , %s" % ( shape1 , shape2  )

        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## normalize
        if normalize and self.normalize :
            uds1, uds2  = normalize_pooled ( uds1 , uds2  ) 

        w1_trivial = weight_trivial ( weight1 )
        w2_trivial = weight_trivial ( weight2 )

        nw1 = len ( uds1 )
        nw2 = len ( uds2 )
        
        if not w1_trivial :
            sumw  = numpy.sum ( weight1 )
            sumw2 = numpy.sum ( weight1 * weight1 )
            print ( 'W1-NOTRIVIAL' , typename ( self ) ,  sumw , sumw2 ) 
            nw1   = math.floor ( float ( sumw * sumw / sumw2) )
            
        if not w2_trivial :
            sumw  = numpy.sum ( weight2 )
            sumw2 = numpy.sum ( weight2 * weight2 )
            print ( 'W2-NOTRIVIAL' , typename ( self ) ,  sumw , sumw2 ) 
            nw2   = math.floor ( float ( sumw * sumw / sumw2 ) )
        
        v1 = self.np2vstat ( uds1 , weight1 )
        v2 = self.np2vstat ( uds2 , weight2 )

        return Ostap.Math.hotelling ( v1 , int ( nw1 ) , v2 , int ( nw2 ) )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
