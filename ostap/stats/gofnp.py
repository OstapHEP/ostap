#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gofnp.py
#  Base class for implementation of Two-Samples/Goddness-of-Fit methos
#  @date   2024-09-16
# =============================================================================
""" Base class for implementatio Two-Samples/Goddness-of-Fit methos
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-29"
__all__     = (
    'GoFnp'           , ## A base class for numpy-related family of methods to probe goodness-of-fit
)
# =============================================================================
from   ostap.core.ostap_types   import string_types, num_types 
from   ostap.core.core          import SE, VE, Ostap, hID  
from   ostap.stats.counters     import EffCounter, ECDF
from   ostap.utils.utils        import split_n_range
from   ostap.utils.core         import typename
from   ostap.utils.basic        import numcpu, num_jobs, run_parallel 
from   ostap.utils.config       import Config
from   ostap.stats.gof          import AGoFnp
from   ostap.stats.utils        import ( weight_trivial     ,
                                         check_all          , 
                                         valid_data_shape   ,
                                         num_features       ,
                                         num_samples        ,
                                         np2vct             ) 
from   ostap.stats.gof_utils    import ( run_parallel       ,
                                         num_jobs           , 
                                         normalize_pooled   ,
                                         pairwise_distances ,
                                         nearest_neighbors  , 
                                         nearest_distances  , 
                                         draw_ecdf          , s2u ) 
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
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gofnp' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Base class for implementation of Two-Samples/Goddness-of-Fit methos' ) 
# =============================================================================
## @class GoFnp 
#  A base class for numpy-related family of methods to probe goodness-of-fit
class GoFnp (AGoFnp,Config) :
    """ A base class for numpy-related family of methods to probe goodness-of-fit
    """
    def __init__ ( self               , * , 
                   nToys     = 100    ,
                   silent    = True   , 
                   parallel  = True   ,
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
            if not weight_trivial ( weight1 ) : raise ValueError ( "weight1 must be *trivial*" ) 
            if not weight_trivial ( weight2 ) : raise ValueError ( "weight2 must be *trivial*" )
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
        
        if self.parallel and resampler.run : counter , _ = resampler.run ( self.nToys , progress = self.progress , silent = self.silent )            
        else                               : counter , _ = resampler     ( self.nToys , progress = self.progress , silent = self.silent )

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

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
        
