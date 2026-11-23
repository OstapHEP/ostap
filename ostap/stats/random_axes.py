#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/random_axes.py
#  Implement Random Axes to probe
#  the difference between two (weighted) datasets
#  
# The theoretical justification for projecting multidimensional datasets 
# onto random axes lies in the Cramér–Wold Theorem.
# 
# Theorem Statement:
#
#  A multivariate probability distribution is uniquely determined by 
#  the set of all its one-dimensional linear projections. 
#  In other words, two random vectors $X$ and $Y$ have the same joint 
#  distribution if and only if their scalar projections 
#  $v^T X$ and $v^T Y$ have the same distribution
#  for all possible unit vectors $v \in \mathbb{R}^D$.
# 
#  @see Cramér, H., & Wold, H. (1936). 
#       "Some Theorems on Distribution Functions". 
#       Journal of the London Mathematical Society, 
#        s1-11(4), 290–294.Read on Wiley Online Library
#  @see  https://doi.org/10.1112/jlms/s1-11.4.290
#
#  Cramér–Wold Theorem on Wikipedia
#  @see https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93Wold_theorem
#  
#  @see Cuesta-Albertos, J. A., Escudero, L., & Matrán, C. (2007). 
#       "A-tests: Multivariate goodness-of-fit tests based on projections Test", 16(1), 21–42 
#  @see https://doi.org/10.1016/j.csda.2006.09.007
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-04
# =============================================================================
""" Implement Random Axes to probe
    the difference between two (weighted) datasets 
    
To estimate the p-value, permutations are used.

The theoretical justification for projecting multidimensional datasets 
onto random axes lies in the Cramér–Wold Theorem.
 
Theorem Statement:

A multivariate probability distribution is uniquely determined by 
the set of all its one-dimensional linear projections. 
In other words, two random vectors $X$ and $Y$ have the same joint 
distribution if and only if their scalar projections 
$v^T X$ and $v^T Y$ have the same distribution
for all possible unit vectors

- see Cramér, H., & Wold, H. (1936). 
        "Some Theorems on Distribution Functions". 
        Journal of the London Mathematical Society, 
         s1-11(4), 290–294.Read on Wiley Online Library
- see  https://doi.org/10.1112/jlms/s1-11.4.290

- Cramér–Wold Theorem on Wikipedia
- see https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93Wold_theorem

- see Cuesta-Albertos, J. A., Escudero, L., & Matrán, C. (2007). 
      "A-tests: Multivariate goodness-of-fit tests based on projections Test", 16(1), 21–42 
- see https://doi.org/10.1016/j.csda.2006.09.007

"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2026-11-15"
__all__     = (
    'RandomAxes' , ## Random Axes Cramer-Wold Two-Samples Tests 
) 
# =============================================================================
from   ostap.core.ostap_types import sequence_types 
from   ostap.utils.core       import typename 
from   ostap.utils.basic      import classprop, LoopGC
from   ostap.stats.pvalue     import make_permutations, permutation_symbol 
from   ostap.utils.config     import Config
from   ostap.stats.gof        import AGoFnp
from   ostap.stats.counters   import SE, EffCounter
from   ostap.stats.pvalues    import combine_methods , combine_pvalues
import ostap.stats.utils      as     U
import ostap.logger.symbols   as     S 
import numpy
# =============================================================================
# Logging setup
# =============================================================================
from ostap.logger.logger import getLogger
logger = getLogger ( 'ostap.stats.random_axes' )
# =============================================================================
logger.debug ( 'Implement Random Axes/Cramer-Wold Two-Samples Test' )
# =============================================================================
# Memory threshold for dense projection matrices (200 MB)
MAX_PROJECTION_BYTES = 200 * 1024 * 1024
# =============================================================================
## @class RandomAxes
# Two-Sample Test based on the Cramer-Wold theorem
#
# Theorem Statement:
#
#  A multivariate probability distribution is uniquely determined by 
#  the set of all its one-dimensional linear projections. 
#  In other words, two random vectors $X$ and $Y$ have the same joint 
#  distribution if and only if their scalar projections 
#  $v^T X$ and $v^T Y$ have the same distribution
#  for all possible unit vectors $v \in \mathbb{R}^D$.
# 
#  @see Cramér, H., & Wold, H. (1936). 
#       "Some Theorems on Distribution Functions". 
#       Journal of the London Mathematical Society, 
#        s1-11(4), 290–294.Read on Wiley Online Library
#  @see  https://doi.org/10.1112/jlms/s1-11.4.290
#
#  Cramér–Wold Theorem on Wikipedia
#  @see https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93Wold_theorem
#  
#  @see Cuesta-Albertos, J. A., Escudero, L., & Matrán, C. (2007). 
#       "A-tests: Multivariate goodness-of-fit tests based on projections Test", 16(1), 21–42 
#  @see https://doi.org/10.1016/j.csda.2006.09.007
#
class RandomAxes(AGoFnp,Config) :
    """ Two sample Test based on the Cramer-Wold theorem
    
    Theorem Statement:
    
    A multivariate probability distribution is uniquely determined by 
    the set of all its one-dimensional linear projections. 
    In other words, two random vectors $X$ and $Y$ have the same joint 
    distribution if and only if their scalar projections 
    $v^T X$ and $v^T Y$ have the same distribution
    for all possible unit vectors
    
    - see Cramér, H., & Wold, H. (1936). 
        "Some Theorems on Distribution Functions". 
        Journal of the London Mathematical Society, 
         s1-11(4), 290–294.Read on Wiley Online Library
    - see  https://doi.org/10.1112/jlms/s1-11.4.290

    - Cramér–Wold Theorem on Wikipedia
    - see https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93Wold_theorem

    - see Cuesta-Albertos, J. A., Escudero, L., & Matrán, C. (2007). 
      "A-tests: Multivariate goodness-of-fit tests based on projections Test", 16(1), 21–42 
    - see https://doi.org/10.1016/j.csda.2006.09.007

    """
    def __init__ ( self                 , * ,
                   comparators          ,   
                   nToys        =  100  ,
                   nAxes        =   50  , 
                   silent       = True  ,
                   progress     = True  , 
                   random_state = None  , **params ) :

        if not isinstance ( nToys   , int ) : raise TypeError  ( "Invalid type  for `nToys` : %s" % typename ( nToys   ) )
        if not 0 <= nToys                   : raise ValueError ( "Invalid value for `nToys` : %d" %            nToys     )
        if not isinstance ( nAxes   , int ) : raise TypeError  ( "Invalid type  for `nAxes` : %s" % typename ( nAxes   ) )
        if not 2 <= nAxes                   : raise ValueError ( "Invalid value for `nAxes` : %d" %            nAxes     )
        
        self.__nToys        = nToys
        self.__nAxes        = nAxes 
        self.__random_state = random_state 
        
        # Lazy cache of unit projection axes of shape (nAxes, nDim)
        self.__axes         = None   
        self.__rng          = numpy.random.default_rng ( random_state )

        if not isinstance ( comparators , sequence_types ) : comparators = [ comparators ] 
            
        self.__comparators  = tuple ( c for c in comparators if isinstance ( c , AGoFnp ) and  c.two_samples )
        if not self.comparators : raise TypeError ( "No valid comparators are provided!")

        self.__weights_supported = all ( c.weights_supported  for c in self.comparators )
        
        self.__t_counters = tuple ( EffCounter () for _ in range ( self.nComps ) )
                
        ##
        self.__progress  = True if progress  else False
        self.__counter   = EffCounter  () 
        self.__pvalue    = None 
        
        Config.__init__ ( self, 
                          silent       = silent          , 
                          progress     = self.__progress , 
                          random_state = random_state    ) 
        
    # =========================================================================
    ## Symbol/tag associated with this estimator
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_RA

    # =========================================================================
    ## Actual method implemented by the estimator
    @property
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_RA
            
    # =========================================================================
    ## Are weights supported by this GoF estimator?
    @property
    def weights_supported ( self ) :
        """`weghts_supported`: Are weights supported by this estimator?
        """
        return self.__weights_supported

    # =========================================================================
    ## Good for two-samples comparison?
    #  Can this estimator be used for comparison of two samples?
    @property
    def two_samples ( self ) :
        """`two_samples`: Can this estimator be used for comparison of two samples?
        """
        return True 

    # =========================================================================
    ## max-dimensionality
    @property
    def max_dimensionality ( self ) :
        """`max_dimensionality` : maximal dimensionality of dataset
        """
        return -1 

    @property
    def nToys  ( self ) : 
        """`nToys`: number of toys/permutations to get the p-value
        """
        return self.__nToys
    
    @property
    def nAxes  ( self ) : 
        """`nAxes`: number of random axes/projections to probe
        """
        return self.__nAxes
    
    @property
    def nComps ( self ) :
        """`nComps` : number of comparators for projected datasets
        """
        return len ( self.__comparators )
    
    @property
    def random_state ( self ) :
        """`random_state` : random state for random-axes and permuttions
        """
        return self.__random_state
    
    @property
    def axes ( self ) :
        """`axes` : set of random axes/projection to probe or None  
        """
        return self.__axes 
    
    @property
    def comparators ( self ) :
        """`comparators` : comparators to compate projected datasets
        """
        return self.__comparators 
        
    @property
    def counter ( self ) : 
        """`counter`: global counter for p-value
        """
        return self.__counter
    
    @property
    def progress ( self ) : 
        """`progress`: show progress-bar? 
        """
        return self.__progress 

    @property
    def t_counters ( self ) : 
        """`t_counters`: counter of t-values for each comparator
        """
        return self.__t_counters

    @property
    def impornace_features ( self ) :
        """`importance_features` : fictive empty dictionary fo importance features
        - importance features are not supported for this comparator
        """
        return {}

    # =========================================================================
    ## Decide whether dense projection matrices fit under the memory budget
    def use_dense_projections ( self   ,
                                n1     ,
                                n2     ,
                                n_axes ,
                                dtype  = numpy.float64 ) :
        """ Check if (N1 + N2) * nAxes projection data stays below MAX_PROJECTION_BYTES.
        """
        elem_size   = numpy.dtype ( dtype ).itemsize
        total_bytes = ( n1 + n2 ) * n_axes * elem_size
        return total_bytes < MAX_PROJECTION_BYTES

    # =========================================================================
    ## Internal helper to compute t-matrix directly from precomputed projections
    def t_matrix_from_proj ( self    ,
                             proj1   ,
                             proj2   , * ,
                             weight1 = None ,
                             weight2 = None ) :
        """ Compute t-matrix from already projected 1D arrays per axis.
        """
        n_axes   = self.nAxes
        n_comps  = self.nComps
        t_matrix = numpy.empty ( ( n_comps , n_axes ) , dtype = numpy.float64 )
        
        for i in range ( n_axes ) :
            u = proj1 [ : , i ]
            v = proj2 [ : , i ]
            for c_idx , comp in enumerate ( self.comparators ) :
                t_matrix [ c_idx , i ] = comp.tvalue ( data1     = u       , 
                                                       data2     = v       , 
                                                       weight1   = weight1 , 
                                                       weight2   = weight2 , 
                                                       normalize = False   )
        return t_matrix

    # =========================================================================
    ## Internal helper to compute t-matrix by projecting onto each axis on the fly
    def t_matrix_on_the_fly ( self    ,
                              data1   ,
                              data2   , * ,
                              weight1 = None ,
                              weight2 = None ) :
        """ Project samples onto each axis one by one to avoid storing full (N, nAxes) matrices.
        """
        n_axes   = self.nAxes
        n_comps  = self.nComps
        t_matrix = numpy.empty ( ( n_comps , n_axes ) , dtype = numpy.float64 )

        for i in range ( n_axes ) :
            axis = self.axes [ i ]
            # Compute 1D projections only for the current axis and release them right after
            u    = numpy.dot ( data1 , axis )
            v    = numpy.dot ( data2 , axis )
            for c_idx , comp in enumerate ( self.comparators ) :
                t_matrix [ c_idx , i ] = comp.tvalue ( data1     = u       ,
                                                       data2     = v       ,
                                                       weight1   = weight1 ,
                                                       weight2   = weight2 ,
                                                       normalize = False   )
            del u , v

        return t_matrix
    
    # =========================================================================
    ## Public helper to compute t-matrix from precomputed projections
    def t_matrix ( self    ,
                   proj1   ,
                   proj2   , * ,
                   weight1 = None ,
                   weight2 = None ) :
        """ Compute t-matrix from already projected data 
        """
        return self.t_matrix_from_proj ( proj1     = proj1   ,
                                         proj2     = proj2   ,
                                         weight1   = weight1 ,
                                         weight2   = weight2 )
    
    # =========================================================================
    ## Calculate the matrix (#comparators x #axes) of t-values
    #  @code
    #  gof       = ...
    #  ds1 , ds2 = ...
    #  t_matrix  = gof.tvalue ( ds1 , ds2 , normalize = True )
    #  @endcode
    def tvalue ( self      ,
                 data1     ,
                 data2     , *     ,
                 weight1   = None  ,
                 weight2   = None  , 
                 normalize = True  ,
                 importace = False ) :
        """ Calculate the t-value
        >>> gof = ...
        >>> ds1 , ds2 = ...
        >>> t   = gof.tvalue ( ds1 , ds2 , normalize = True )
        """
        U.check_all ( data1   = data1    , 
                      data2   = data2    , 
                      weight1 = weight1  , 
                      weight2 = weight2  , 
                      where   = typename ( self ) )
        
        if not self.weights_supported : 
            if not U.weight_trivial ( weight1 ) : raise ValueError  ( "Weights are not supprted!")
            if not U.weight_trivial ( weight2 ) : raise ValueError  ( "Weights are not supprted!")
           
        n_dim  = U.num_features ( data1 )
        n_axes = self.nAxes     
        
        if self.axes is None or self.axes.shape != ( n_axes , n_dim ) :
            self.__axes = self.generate_axes ( n_dim , n_axes )
        
        ## normalize input data
        if normalize : 
            data1 , data2 = U.normalize_pooled  ( data1   = data1   ,
                                                  data2   = data2   ,
                                                  weight1 = weight1 ,
                                                  weight2 = weight2 )

        proj_dtype = numpy.result_type ( numpy.asarray ( data1 ).dtype ,
                                         numpy.asarray ( data2 ).dtype ,
                                         self.axes.dtype )

        # Use dense projection matrices only when they stay within the 500 MB budget
        if self.use_dense_projections ( len ( data1 ) , len ( data2 ) , n_axes , proj_dtype ) :
            proj1 = numpy.dot ( data1 , self.axes.T )
            proj2 = numpy.dot ( data2 , self.axes.T )
            res   = self.t_matrix_from_proj ( proj1     = proj1   ,
                                              proj2     = proj2   ,
                                              weight1   = weight1 ,
                                              weight2   = weight2 )
            del proj1 , proj2
            return res

        # Fall back to axis-by-axis projection for large datasets / many axes
        return self.t_matrix_on_the_fly ( data1     = data1   ,
                                          data2     = data2   ,
                                          weight1   = weight1 ,
                                          weight2   = weight2 )

    
    __call__ = tvalue
    
    # =========================================================================
    ## Calculate the t (matrix) & p (scalar) values using permutation test
    #  @code
    #  gof = ...
    #  ds1 , ds2 = ...
    #  t , p = gof.pvalue ( ds1 , ds2 )
    #  @endcode
    def pvalue ( self     ,
                 data1    ,
                 data2    , *    ,
                 weight1  = None ,
                 weight2  = None ,
                 tvalue   = None ) :
        """ Calculate the t (matrix) & p (scalar) values using permutation test
        >>> gof = ...
        >>> ds1 , ds2 = ...
        >>> t , p = gof.pvalue ( ds1 , ds2 , normalize = True )
        """
        
        ## reset/re-create the counters 
        self.__counter.reset()  
        for c in self.__t_counters : c.reset()
        self.__pvalue    = None 
        
        if not self.weights_supported : 
            if not U.weight_trivial ( weight1 ) : raise ValueError  ( "Weights are not supported!")
            if not U.weight_trivial ( weight2 ) : raise ValueError  ( "Weights are not supported!")
            weight1 , weight2 = None , None

        nc = self.nComps
        na = self.nAxes 
        
        U.check_all ( data1   = data1    , 
                      data2   = data2    , 
                      weight1 = weight1  , 
                      weight2 = weight2  , 
                      where   = typename ( self ) )
        
        data1 , data2 = U.normalize_pooled  ( data1   = data1   ,
                                              data2   = data2   ,
                                              weight1 = weight1 ,
                                              weight2 = weight2 )
        
        if tvalue is None or tvalue.shape != ( nc , na ) :
            tvalue = self.tvalue (  data1     = data1   ,  
                                    data2     = data2   , 
                                    weight1   = weight1 , 
                                    weight2   = weight2 , 
                                    normalize = False   )
        
        tv_max      = numpy.max ( tvalue , axis = 1 )
        description = permutation_symbol  + ' ' + S.fixed_width ( self.symbol , 4 ) + ':'
        
        proj_dtype  = numpy.result_type ( numpy.asarray ( data1 ).dtype ,
                                          numpy.asarray ( data2 ).dtype ,
                                          self.axes.dtype )
        
        dense_proj  = self.use_dense_projections ( len ( data1 ) , len ( data2 ) , na , proj_dtype )

        if dense_proj :
            # Precompute full projection matrices once when memory budget allows
            proj1 = numpy.dot ( data1 , self.axes.T )
            proj2 = numpy.dot ( data2 , self.axes.T )
            perm_data1 , perm_data2 = proj1 , proj2
        else :
            # Keep original feature matrices and project per axis inside each toy
            proj1  = proj2 = None
            perm_data1 , perm_data2 = data1 , data2

        ## make permutations with LoopGC over toy stream
        for ds1, ds2, w1 , w2 in LoopGC ( make_permutations ( self.nToys  ,
                                                              data1       = perm_data1      , 
                                                              data2       = perm_data2      , 
                                                              weight1     = weight1         , 
                                                              weight2     = weight2         , 
                                                              progress    = self.progress     , 
                                                              random_seed = self.random_state ,
                                                              description = description       ) ) : 
            
            ## get new matrix of t-values 
            if dense_proj :
                tmatrix = self.t_matrix_from_proj ( proj1   = ds1 ,
                                                    proj2   = ds2 , 
                                                    weight1 = w1  , 
                                                    weight2 = w2  )
            else :
                tmatrix = self.t_matrix_on_the_fly ( data1   = ds1 ,
                                                     data2   = ds2 ,
                                                     weight1 = w1  ,
                                                     weight2 = w2  )
                
            tm_max = numpy.max ( tmatrix , axis = 1 )

            exceeds = tv_max <= tm_max

            ## update counters 
            for exc, counter in zip ( exceeds, self.t_counters ) : counter += exc

            ## any ? all ? 
            self.__counter += numpy.any ( exceeds )

            ## clean up per-toy temporaries
            del tmatrix , tm_max , exceeds , ds1 , ds2

        if dense_proj :
            del proj1 , proj2
                        
        pvalue        = self.__counter.eff 
        self.__pvalue = pvalue 
        
        return tvalue , pvalue 
    
    # ========================================================================
    ## generate <code>n_axes</code>  random axes of dimension  <code>n_dim</code>
    def generate_axes ( self , n_dim  , n_axes ) :
        """  Generate `n_axes` random axes on dimension `dim` each  
        """
        axes   = self.__rng.normal  ( loc = 0.0 , scale = 1.0 , size = ( n_axes , n_dim ) )
        norms  = numpy.linalg.norm ( axes , axis = 1 , keepdims = True )
        norms [ norms == 0.0 ] = 1.0
        axes  /= norms
        ## 
        return axes

    # ==========================================================================
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
        
        pvalue = self.__pvalue if pvalue is None else pvalue 
        return super().the_row ( pvalue    = pvalue    , 
                                 precision = precision ,      
                                 width     = width     )
 
    # ==========================================================================
    ## Get results in form of the table 
    def report ( self             ,
                 tvalue    = None ,
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

        rows = []
        header =  ( '' , 'Method' )
        for cmp , pv in zip ( self.comparators ,
                              self.t_counters  ) : 
        
            hdr , row = cmp.the_row ( pvalue    = pv.eff      ,
                                      precision = precision   ,
                                      width     = width       )
            
            row = ( cmp.symbol , cmp.method ) + row
            
            rows.append ( row )

        ## combination methods
        crows = [] 
        pvals = tuple ( p.eff for p in self.t_counters )
        for method in combine_methods :
            pval   = combine_pvalues ( pvals  , method )
            pv     = U.clip_pvalue   ( pval   )
            nsigma = U.significance  ( pv     )
            pval   = U.format_pvalue ( pval   ) 
            nsigma = U.format_nsigma ( nsigma )
            row    = pval , nsigma
            m      = U.method_name ( method.title() )
            row    = ( S.o_plus , m )  + row
            crows.append ( row ) 
            
        _ , row = self.the_row ( pvalue    = self.__pvalue , 
                                 precision = precision     ,
                                 width     = width         ) 
        
        first = header + hdr 
        last  = ( self.symbol, self.method ) + row        
        rows  = [ first ] + \
            sorted ( rows  , key = lambda row : row [ 1 ] ) +  \
            sorted ( crows , key = lambda row : row [ 1 ] ) + [ last ]
        
        import ostap.logger.table as T
        rows  = T.remove_empty_columns ( rows )
        title = title if title else '%s Random Axes/%s #axes=%d #comps=%d' % ( self.symbol, self.method , self.nAxes, self.nComps ) 
        return T.table ( rows , title = title , prefix = prefix , alignment = 'llcccccccc' , style = style )
        
    def __str__  ( self ) : return self.report ()
    def __repr__ ( self ) : return self.report ()

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END
# =============================================================================
