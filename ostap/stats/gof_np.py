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
    'Mahalanobis'     , ## Very crude estimator based on Mahalanobis' distance
    'KullbackLeibler' , ## Very crude estimator based on Kullback-Leibler's divergency 
    'Hotelling'       , ## Very crude estimator based on Hotelling's distance 
)
# =============================================================================
from   ostap.core.ostap_types   import string_types, num_types 
from   ostap.core.core          import SE, VE, Ostap, hID  
from   ostap.stats.counters     import EffCounter, ECDF
from   ostap.utils.utils        import split_n_range
from   ostap.utils.core         import typename
from   ostap.utils.basic        import numcpu
from   ostap.utils.config       import Config
from   ostap.stats.gof          import AGoFnp
from   ostap.stats.utils        import ( weight_trivial     ,
                                         check_all          , 
                                         valid_data_shape   ,
                                         num_features       ,
                                         num_samples        ) 
from   ostap.stats.gof_utils    import ( run_parallel       ,
                                         num_jobs           , 
                                         normalize_pooled   ,
                                         pairwise_distances ,
                                         nearest_neighbors  , 
                                         nearest_distances  , 
                                         draw_ecdf          , s2u , np2vct ) 
from   ostap.utils.memory       import memory, memory_enough
from   ostap.math.math_ve       import gauss_cdf
from   ostap.logger.symbols     import ( symmetry  as symmetry_symbol  ,
                                         asymmetry as asymmetry_symbol )
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
class GoFnp (AGoFnp,Config) :
    """ A base class for numpy-related family of methods to probe goodness-of-fit
    """
    def __init__ ( self               , * , 
                   nToys     = 100    ,
                   silent    = False  , 
                   parallel  = False  ,
                   method    = 'GoF'  ,
                   progress  = True   ,
                   normalize = True   , **params ) : 
        
        if not isinstance ( nToys   , int ) : raise TypeError  ( "Invalid type  for `nToys`  : %s" % typename ( nToys   ) )
        if not 0 <= nToys                   : raise ValueError ( "Invalid value for `nToys`  : %d" %            nToys     )

        self.__nToys    = nToys
        ## 
        self.__parallel  = True if parallel  else False
        self.__progress  = True if progress  else False
        self.__normalize = True if normalize else False 
        ## 
        self.__method    = method

        ## Empirical CDF for t-value distribution from permutations/toys"""
        self.__ecdf      = None
        self.__counter   = None
        self.__tvalue    = None
        self.__pvalue    = None
                
        if self.__parallel :
            mratio = memory_enough ()  
            if mratio <= 1 :
                logger.warning ( 'Available/used memory ratio: %.1f; switch-off parallel processing' % mratio )                
                ## self.__parallel = False

        ## initiailze the base
        Config.__init__ ( self , silent = silent , **params ) 
                
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
        conf [ 'nToys'                      ] = self.nToys
        conf [ 'progress'                   ] = self.progress
        conf [ 'parallel'                   ] = self.parallel
        conf [ 'normalize'                  ] = self.normalize
        conf [ 'method'                     ] = self.method  
        conf [ 'weights_supported'          ] = self.weights_supported
        conf [ 'negative_weights_supported' ] = self.negative_weights_supported
        conf [ 'two_samples'                ] = self.two_samples 
        return conf 
    
    # =========================================================================
    @property 
    def nToys ( self ) :
        """`nToys` : number of permutations/toys used for permutation/toys test"""
        return self.__nToys

    # =========================================================================
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
            if not weight_trivial ( weight1 ) : raise ValueError ( "weight1 must be *trivial*" ) 
            if not weight_trivial ( weight2 ) : raise ValueError ( "weight2 must be *trivial*" ) 
            weight1 = None
            weight2 = None

        ## transform ?  
        uds1 , uds2 = self.unpack ( data1 , data2 ) 
        
        ## check vaildity and consitency of input parameters 
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) )

        ## normalize
        if normalize and self.normalize :
            uds1 , uds2 = self.normalize_pooled ( uds1 , uds2 ) 
        
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
    def pvalue ( self           , 
                 data1          ,
                 data2          , * ,
                 tvalue  = None , 
                 weight1 = None ,
                 weight2 = None ) : 
                
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
    
        ## check vaildity and consitency of input parameters 
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) ) 

        ## normalize ? 
        if self.normalize :
            uds1 , uds2 = self.normalize_pooled ( uds1 , uds2 ) 

        ## calculate t-value if not specified 
        t_value    = tvalue if not tvalue is None else self.tvalue ( uds1      ,
                                                                     uds2      ,
                                                                     weight1   = weight1 ,
                                                                     weight2   = weight2 ,
                                                                     normalize = False   )        
        ## use permutations to get the p-value
        from ostap.stats.pvalue import PERMUTATOR as RESAMPLER 
        ## from ostap.stats.pvalue import BOOTSTRAPPER as RESAMPLER 
        resampler = RESAMPLER ( self                   ,
                                t_value                , 
                                uds1                   ,
                                uds2                   ,
                                weight1 = weight1      ,
                                weight2 = weight2      ) ;
        
        if self.parallel and resampler.run : counter , _ = resampler.run ( self.nToys , progress = self.progress )            
        else                               : counter , _ = resampler     ( self.nToys , progress = self.progress )

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

        self.ecdf    = resampler.ecdf
        self.counter = counter
        
        self.t_value = t_value 
        self.p_value = p_value 
        
        return self.t_value , self.p_value
    
    @property
    def ecdf ( self ) :
        """`ecdf` : empirical CDF for t-value distribution from permutations"""
        return self.__ecdf
    @ecdf.setter
    def ecdf ( self , value ) :
        assert value is None or isinstance ( value , ECDF ) , \
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
                                title   = title  if title else '%s GoF-report [#%d]' %  ( typename ( self ) , self.nToys ) ,
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
                                 counter = counter if not counter is None else self.__counter )
    

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

    # ========================================================================
    ## Normalize/standartize two numpy arrays such that mean and rms for the pooled sample
    #  are equal to 0 and 1 correspondinly
    #  @code
    #  ds1 = ...
    #  ds2 = ...
    #  ds1 , ds2 = gof.normalize_pooled  ( ds1 , ds2 ) 
    #  @endcode 
    def normalize_pooled ( self , ds1 , ds2 ) :
        """ Normalize/standartize two numpy arrays such that the mean and rms for the POOLED sample
        are equal to 0 and 1 correspondinly
        
        >>> ds1 = ...
        >>> ds2 = ...
        >>> ds1 , ds2 = gof.normalize_pooled  ( ds1 , ds2 )        
        """
        return normalize_pooled ( ds1 , ds2 )
        
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
        
        ## Attention!
        assert isinstance ( n_neighbors , int ) and 2 <= n_neighbors , \
            "Invalid `n_neighbors`: %s" % n_neighbors 

        ## store it
        self._k_max = n_neighbors 
        ##
        
        # =================================================================================
        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "%s: Parallel processing is switched OFF!" % typename ( self ) ) 
            parallel = False 
                    
        n_jobs = 1 if parallel else num_jobs ( params , numcpu() - 1 )

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
        return False 

    # =========================================================================
    ## Are weights supported by this GoF estimator?
    @property 
    def negative_weights_supported ( self ) :
        """`negative_weghts_supported`: Are weights supported by this estimator?
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

        ## check validity and consitency of input parameters 
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
                   parallel  = False      , 
                   silent    = False      ,
                   progress  = True       , 
                   maxsize   = 10000000   , **params ) :

        
        # =================================================================================
        if parallel and not run_parallel ( parallel ) :
            logger.warning ( "%s: Parallel processing is switched OFF!" % typename ( self ) ) 
            parallel = False 
                    
        n_jobs = 1 if parallel else num_jobs ( params , numcpu() - 1 )
        
        self.__mc2mc     = True if mc2mc else False
        self.__transform = None
        self.__sigma     = sigma
        self.__psi       = psi
        assert isinstance ( maxsize , int ) and 0 < maxsize , "Invalid `maxsize' : %s" % maxsize
        
        self.__maxsize   = max ( maxsize , 100000  )
        
        ## check validity of `psi`
        scale = -0.5 / ( self.sigma ** 2 ) 
        self.__distance_type , _ , _ = psi_conf ( psi , scale )

        GoFnp.__init__ ( self                 ,
                         nToys     = nToys    ,
                         parallel  = parallel , 
                         silent    = silent   ,
                         progress  = progress , ## ATTENTION!                          
                         normalize = True     ,
                         n_jobs    = n_jobs   , 
                         method    = 'Point-to-Point Dissimilarity' , **params )
                
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
    ## Are weights supported by this GoF estimator?
    @property 
    def negative_weights_supported ( self ) :
        """`negative_weghts_supported`: Are weights supported by this estimator?
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
    def __init__ ( self              ,
                   histo    = None   ,
                   nToys    = 1000   ,
                   parallel = False  , 
                   silent   = True   ,
                   progress = True   , **params ) :
        
        n_jobs = 1 if parallel else num_jobs ( params , numcpu() - 1 )

        if 'metric' in params : params.pop ( 'metric' )
        if 'p'      in params : params.pop ( 'p'      )
        
        self.__histo = None 
        if   isinstance ( histo , ROOT.TH1 ) :
            self.__histo = histo
        elif isinstance ( histo , int      ) and 1 < histo :
            self.__histo = ROOT.TH1D ( hID () , 'U-values' , histo , -0.05 , 1.05 ) 

        GoFnp.__init__ ( self                      ,
                         nToys       = nToys       ,
                         parallel    = parallel    , 
                         silent      = silent      ,
                         progress    = progress    ,
                         normalize   = True        ,                          
                         method      = 'Distance-to-Nearest-Neighbor' ,
                         n_jobs      = n_jobs      , **params )
        
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
    ## Are weights supported by this GoF estimator?
    @property 
    def negative_weights_supported ( self ) :
        """`negative_weghts_supported`: Are weights supported by this estimator?
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
## @class Mahalanobis
#  Use Mahalanobis distance to discriminiate the dataset
#  @attention it is *VERY* crude "estimator"
class Mahalanobis(GoFnp) :
    """ Use Mahalanobis distance to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self        , *     , 
                   nToys       = 400   ,
                   parallel    = False , 
                   silent      = False ,
                   method      = "Mahalanobis" , 
                   normalize   = True  , 
                   progress    = True  , **params ) :         
        
        ## initialize the base 
        GoFnp.__init__ ( self                          , 
                         nToys        = nToys          ,
                         parallel     = parallel       , 
                         silent       = silent         ,
                         progress     = progress       ,                         
                         method       = method         ,
                         normalize    = normalize      , **params )

    # =========================================================================
    ## Are weights supported by this estimator?
    @property
    def weights_supported ( self ) :
        """`weights_supported` : Are weights supported by this estimator?
        """
        return True
    
    # =========================================================================
    ## Are weights supported by this GoF estimator?
    @property 
    def negative_weights_supported ( self ) :
        """`negative_weghts_supported`: Are weights supported by this estimator?
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
        return np2vct ( data , weight = weight )
                        
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
        
        return v1.mahalanobis ( v2 )
    
# ============================================================================
## @class KullbackLeibler 
#  Use KullbackLeibler divergency to discriminiate the dataset
#  @attention it is *VERY* crude "estimator"
class KullbackLeibler(Mahalanobis) :
    """ Use Kullback-Leibler divergency to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self        ,  *    , 
                   nToys       = 400  ,
                   parallel    = False , 
                   silent      = False ,
                   progress    = True  ,
                   symmetric   = True  , **params ) :         
        
        self.__symmetric = True if symmetric else False

        method = 'Kullback-Leibler/%s' % ( symmetry_symbol if self.symmetric else asymmetry_symbol )
        
        ## initialize the base 
        super() .__init__ ( nToys        = nToys     ,
                            parallel     = parallel  , 
                            silent       = silent    ,
                            progress     = progress  ,                           
                            method       = method    , 
                            symmetric    = symmetric , 
                            normalize    = True      , **params )

    # =========================================================================
    @property
    def symmetric ( self ) :
        """`symmetric` : Symemtric Kullback-Leibler distance?"""
        return self.__symmetric 
        
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
        
        return v1.kullback_leibler ( v2 ) if self.symmetric else v1.asymmetric_kullback_leibler ( v2 )

# ============================================================================
## @class Hotelling  
#  Use Hotelling's t-squared statistics to discriminiate the datasets
#  @attention it is *VERY* crude "estimator"
class Hotelling(Mahalanobis) :
    """ Use Hotelling's t-squared statistics to discriminiate the dataset
    - attention it is *VERY* crude "estimator"
    """    
    def __init__ ( self        , *     ,
                   nToys       = 400   ,
                   parallel    = False , 
                   silent      = False ,
                   progress    = True  , **params ) :         

        ## initialize the base 
        super().__init__ ( nToys        = nToys            ,
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
        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## check everything
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) ) 
        
        ## normalize
        if normalize and self.normalize :
            uds1, uds2  = self.normalize_pooled ( uds1 , uds2  ) 

        w1_trivial = weight_trivial ( weight1 )
        w2_trivial = weight_trivial ( weight2 )

        nw1 = len ( uds1 )
        nw2 = len ( uds2 )
        
        if not w1_trivial :
            
            sumw  = numpy.sum ( weight1 )
            sumw2 = numpy.sum ( weight1 * weight1 )
            nw1   = math.floor ( float ( sumw * sumw / sumw2 ) )
            
        if not w2_trivial :
            
            sumw  = numpy.sum ( weight2 )
            sumw2 = numpy.sum ( weight2 * weight2 )
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
