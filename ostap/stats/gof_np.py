#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof_np.py
#  Set of utulities for goodness-of-fit studies for multidimensional fits
#  @see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
#  @see https://doi.org/10.1088/1748-0221/5/09/P09004
#  @see http://arxiv.org/abs/arXiv:1003.1768 
#  @author Artem Egorychev Artem.Egorychev@cern.ch 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-09-16
# =============================================================================
""" Simple utulities for goodness-of-fit studies for multidimensionao fits 
- see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
- see https://doi.org/10.1088/1748-0221/5/09/P09004
- see http://arxiv.org/abs/arXiv:1003.1768
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-29"
__all__     = (
    'PPDNP' , ## Point-to-Point Dissimilarity Goodness-of-fit method 
)
# =============================================================================
from   ostap.core.ostap_types   import string_types 
from   ostap.stats.gof          import normalize2 
from   ostap.core.core          import SE, VE
from   ostap.utils.progress_bar import progress_bar
from   ostap.utils.utils        import split_n_range, splitter 
from   ostap.utils.basic        import numcpu
from   ostap.stats.counters     import EffCounter 
import os, abc, warnings  
# =============================================================================
try : # =======================================================================
    # =========================================================================
    with warnings.catch_warnings(): 
        warnings.simplefilter ( "ignore" , category = UserWarning  )
        import numpy                                                  as np   
        import scipy                                                  as sp  
        from numpy.lib.recfunctions import structured_to_unstructured as s2u
        from scipy.spatial.distance import cdist                      as cdist 
    # =========================================================================
except ImportError :
    # =========================================================================
    np    = None
    sp    = None    
    s2u   = None
    cdist = None     
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gofnd' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-fit studies for multidimensional fits' )
# =============================================================================
## @class AGoFNP
#  An absract base class for numpy-related family of methods to probe goodness-of fit
#  There are two abstract methods
#  - <code>__call__</code> to evaluate t-value for two datasets 
#  - <code>pvalue</code> to evaluate (t,p)-vaues for two datasets
#  @code
#  gof = ...
#  ds1, ds2 = ...
#  t    = god        ( ds1 , ds2 , normalize = True )
#  t,p  = god.pvalue ( ds1 , ds2 , normalize = True )
#  @endcode
class AGoFNP(object) :
    """ An absract base class for numpy-related family of methods to probe goodness-of fit
    There are two abstract methods
    - `__call__` to evaluate t-value for two datasets 
    - `pvalue` to evaluate (t,p)-vaues for two datasets
    
    >>> gof = ...
    >>> ds1, ds2 = ...
    >>> t    = gof ( ds1 , ds2 , normalize = True )
    >>> t,p  = gof ( ds1 , ds2 , normalize = True )
    """
    # =========================================================================
    ## Calculate T-value for two datasets 
    #  @code
    #  data1 = ... ## the first  data set 
    #  data2 = ... ## the second data set
    #  gof   = ...
    #  t = gof ( data1 , data1 , normalize = False ) 
    #  t = gof ( data1 , data1 , normalize = True  ) 
    #  @endcode
    @abc.abstractmethod
    def __call__ ( self , data1 , data2 , normalize = True ) :
        """ Calculate T-value for two data sets 
        >>> gof    = ...
        >>> data1 = ... ## the first  data set 
        >>> data2 = ... ## the second data set
        >>> t = gof ( data1 , data1 , normalize = False ) 
        >>> t = gof ( data1 , data1 , normalize = True  ) 
        """
        return NotImplemented 
    # =========================================================================
    ## Calculate the t & p-values
    #  @code
    #  gof = ...
    #  ds1 , ds2 = ...
    #  t , p = gof.pvalue ( ds1 , ds2 , normalize = True ) 
    #  @endcode 
    @abc.abstractmethod
    def pvalue ( self , data1 , data2 , normalize = True ) :
        """ Calculate the t & p-values
        >>> gof = ...
        >>> ds1 , ds2 = ...
        >>> t , p = gof.pvalue ( ds1 , ds2 , normalize = True ) 
        """
        return NotImplemented

# =================================================================================
## @class GoFNP 
#  A base class for numpy-relarted family of methods to probe goodness-of fit
class GoFNP (object) :
    """ A base class for numpy-relarted family of methods to probe goodness-of fit
    """
    def __init__ ( self           ,
                   Nperm  = 0     ,
                   silent = False ) : 
        
        assert isinstance ( Nperm , int ) and 0 <= Nperm  , \
            "Invalid number of permulations:%s" % Nperm
        
        self.__Nperm = Nperm
        self.__silent = True if silent else False
        
    # ==========================================================================
    ## Normalize two data sets, such that each variable in pooled set
    #  has a mean of zero and variance of one 
    #  @code
    #  ds1 = ... ## the 1st structured array 
    #  ds2 = ... ## the 2nd structured array
    #  gof = ... ##    #  nd1, nd2 = god.normalize ( ds1 , ds2 ) 
    #  @encode
    def normalize ( self , data1 , data2 ) :
        """ Normalize two data sets, such that each variable in pooled set
        has a  mean of zero and variance of one 
        
        >>> ds1 = ... ## the 1st structured array 
        >>> ds2 = ... ## the 2nd structured array
        >>> gof = ... ##
        >>> nd1, nd2 = god.normalize ( ds1 , ds2 ) 
        """ 
        return normalize2  ( [ data1 , data2 ] , first = False ) 
    # =========================================================================
    ## Generator of permutations
    #  @code
    #  ds1, ds2 = ...
    #  gof      = ...
    #  for d1,d2 in gof.permutations ( 100 , ds1 , ds2 ) :
    #  ...  
    #  @endcode 
    def permutations ( self , Nperm , data1 , data2 ) :
        """ Generator of permutations
        >>> ds1, ds2 = ...
        >>> gof      = ...
        >>> for d1,d2 in gof.permutations ( 100 , ds1 , ds2 ) :
         ...
        """
        n1     = len ( data1 )
        pooled = np.concatenate ( [ data1 , data2 ] )
        ## 
        for i in progress_bar ( Nperm , silent = self.silent ) : 
            np.random.shuffle ( pooled )
            yield pooled [ : n1 ] , pooled [ n1 : ]

        del pooled 
    # =========================================================================
    @property 
    def Nperm ( self ) :
        """`Nperm` : number of permutations used for permutation test"""
        return self.__Nperm
    # ========================================================================
    @property
    def silent  ( self ) :
        """`silent` : silent processing?"""
        return self.__silent
    
# ============================================================================
## define configurtaion for psi-function for PPD method
#   - distance type of <code>cdist</code>
#   - transformation funciton for cdisct output
#   - increasing function ?
#   @code
#   distance_type , transform, increasing = psi_conf ( 'linear' )
#   @endcode
def psi_conf ( psi , scale = 1.0 ) :
    """ Define configuration for psi-function for PPD method"""

    if   psi in ( 'euclidean'  , 'linear'  ) :                 ## psi = x 
        return 'euclidean'   , None                           , True 
    elif psi in ( 'euclidean2' , 'sqeuclidean', 'squared' ) :  ## psi = x**2
        return 'sqeuclidean' , None                           , True  
    elif psi in ( 'inverse'    , 'coulomb' ) :                 ## psi = 1/x 
        return 'euclidean'   , lambda x : 1.0/x               , False 
    elif psi in ( 'log'   , 'logarithm'    ) :                 ## psi = -log(x)
        return 'sqeuclidean' , lambda x : -np.log ( x )       , False 
    elif psi in ( 'gauss' , 'gaussian'     ) :                 ## psi = exp (-x*x/0.5)
        return 'sqeuclidean' , lambda x :  np.exp ( scale*x ) , False 
    elif isinstance ( psi , string_types ) :
        return psi , None , True         

    raise TypeError ( "Unknown `psi':%s" % psi ) 

try :
    import joblib as jl
except ImportError :
    jl = None

# =====================================================================
## @class PERMUTATOR
#  Helper class that allow to run permutattion test in parallel 
class PERMUTATOR(object) :
    """ Helper class that allow to run permutation test in parallel 
    """
    def __init__ ( self, gof, t_value , ds1 , ds2 ) :
        self.gof     = gof
        self.ds1     = ds1
        self.ds2     = ds2
        self.t_value = t_value
        
    # =========================================================================
    ## run N-permutations 
    def __call__ ( self , N , silent = True ) :
        np.random.seed()
        n1      = len ( self.ds1 )
        pooled  = np.concatenate ( [ self.ds1 , self.ds2 ] )
        counter = EffCounter()
        for i in progress_bar ( N , silent = silent ) : 
            np.random.shuffle ( pooled )            
            tv       = self.gof.t_value ( pooled [ : n1 ] , pooled [ n1: ] )
            counter += bool ( self.t_value < tv  )
            
        del pooled
        return counter
    
# =============================================================================
if jl : # =====================================================================
    # =========================================================================
    ## Run NN-permutations in parallel using joblib 
    def joblib_run ( self , NN , silent = True ) :
        """ Run NN-permutations in parallel using joblib """
        nj    = 4 * numcpu () + 4
        lst   = splitter ( NN , nj )
        ## 
        conf  = { 'n_jobs' : -1 , 'verbose' : 0 }
        if    '1.3.0' <= jl.__version__ < '1.4.0' : conf [ 'return_as' ] = 'generator'           
        elif  '1.4.0' <= jl.__version__           : conf [ 'return_as' ] = 'unordered_generator'
        ##
        input   = ( jl.delayed (self)( N ) for N in lst )
        results = jl.Parallel ( **conf ) ( input )
        counter = EffCounter()
        for r in progress_bar ( results , max_value = nj , silent = silent ) :
            counter += r
        # 
        return counter 
    
    PERMUTATOR.run = joblib_run
    # =========================================================================
else : # ======================================================================
    # =========================================================================
    ## Run NN-permutations in parallel using WorkManager
    def pp_run ( self , NN , silent = True ) :
        """ Run NN-permutations in parallel using WorkManager"""
        nj    = 4 * numcpu () + 4
        lst   = splitter ( NN , nj )
        ##
        from ostap.parallel.parallel import WorkManager
        manager = WorkManager ( silent = silent )
        counter = EffCounter()
        ## 
        ## use the bare interface 
        for result in manager.iexecute ( self , lst , progress = not silent  , njobs = nj ) :
            counter += result 
        # 
        return counter 
    
    PERMUTATOR.run = pp_run
    # =========================================================================
    
# =============================================================================
## @class PPD
#  Implementation of concrete method "Point-To-Point Dissimilarity"
#  for probing of Goodness-Of-Fit
#  @see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
#  @see https://doi.org/10.1088/1748-0221/5/09/P09004
#  @see http://arxiv.org/abs/arXiv:1003.1768 
class PPDNP(AGoFNP,GoFNP) : 
    """ Implementation of concrete method "Point-To-Point Dissimilarity"
    for probing of Goodness-Of-Fit
    - see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
    - see https://doi.org/10.1088/1748-0221/5/09/P09004
    - see http://arxiv.org/abs/arXiv:1003.1768 
    """
    def __init__ ( self                   ,
                   mc2mc     = False      ,
                   Nperm     = 1000       ,
                   psi       = 'gaussian' ,
                   sigma     = 0.05       ,
                   silent    = False      ) :
        
        GoFNP.__init__ ( self , Nperm = Nperm , silent = silent )
        self.__mc2mc     = True if mc2mc else False
        self.__transform = None
        self.__sigma     = sigma
        self.__psi       = psi
        
        ## check validity of `psi`
        scale = -0.5/(self.sigma**2) 
        self.__distance_type , _ , self.__increasing = psi_conf ( psi , scale )
        
    # =========================================================================
    @property
    def mc2mc ( self ) :
        """`mc2mc` : add mc<-->mc distances to the T-value 
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
        nnmax  = 1000000
        if n1 * n2 > nnmax  :
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
        scale = -0.5/(self.sigma**2) 
        distance_type , transform , _ = psi_conf ( self.psi , scale )
        ##
        ## calculate all pair-wise distances
        distances = cdist ( data1 , data2 , distance_type ) .flatten () ## data <-> data
        distances = distances [ distances > 0 ] 
        if transform : distances = transform ( distances )
        ## 
        return np.sum ( distances )         
    # =========================================================================
    # calculate t-value for (non-structured) 2D arrays
    def t_value ( self , ds1 , ds2 ) :
        """ Calculate t-value for (non-structured) 2D arrays
        """
        ##
        sh1 = ds1.shape
        sh2 = ds2.shape
        assert 2 == len ( sh1 ) and 2 == len ( sh2 ) and sh1[1] == sh2[1] , \
            "Invalid arrays: %s , %s" % ( sh1 , sh2 )
        
        n1 = len ( ds1 ) 
        n2 = len ( ds2 ) 
        ##

        ## calculate sums of distances, Eq (3.7) 
        result  = self.sum_distances ( ds1 , ds1 ) / ( n1 * ( n1 - 1 ) )
        result -= self.sum_distances ( ds1 , ds2 ) / ( n1 * n2 )
        if self.mc2mc :
            ## add the distances from the second dataset? 
            result += self.sum_distances ( ds2 , ds2 ) / ( n2 * ( n2 - 1 ) )
        ## 
        return result            
    # =========================================================================    
    ## Calculate T-value for two datasets 
    #  @code
    #  ppd   = ...
    #  data1 = ... ## the first  data set 
    #  data2 = ... ## the second data set
    #  t = ppd ( data1 , data1 , normalize = False ) 
    #  t = ppd ( data1 , data1 , normalize = True  ) 
    #  @endcode
    def __call__ ( self , data1 , data2 , normalize = True ) :
        """ Calculate T-value for two data sets 
        >>> ppd   = ...
        >>> data1 = ... ## the first  data set 
        >>> data2 = ... ## the second data set
        >>> t = ppd ( data1 , data1 , normalize = False ) 
        >>> t = ppd ( data1 , data1 , normalize = True  ) 
        """
        
        ## transform/normalize ? 
        if normalize : ds1 , ds2 = self.normalize ( data1 , data2 )
        else         : ds1 , ds2 = data1 , data2 

        ## convert to unstructured datasets 
        uds1    = s2u ( data1 , copy = False )
        uds2    = s2u ( data2 , copy = False )
        
        return self.t_value ( uds1 , uds2 )

    # =========================================================================
    ## Calculate the t & p-values
    #  @code
    #  ppd = ...
    #  ds1 , ds2 = ...
    #  t , p = ppd.pvalue ( ds1 , ds2 , normalize = True ) 
    #  @endcode 
    def pvalue ( self , data1 , data2 , normalize = True ) :
        """ Calculate the t & p-values
        >>> ppd = ...
        >>> ds1 , ds2 = ...
        >>> t , p = ppd.pvalue ( ds1 , ds2 , normalize = True ) 
        """
        
        ## transform/normalize ? 
        if normalize : ds1 , ds2 = self.normalize ( data1 , data2 )
        else         : ds1 , ds2 = data1 , data2 

        ## convert to unstructured datasets 
        uds1 = s2u ( ds1 , copy = False )
        uds2 = s2u ( ds2 , copy = False )

        t_value    = self.t_value ( uds1 , uds2 )

        permutator = PERMUTATOR ( self , t_value , uds1 , uds2 )
        
        if permutator.run : eff = permutator.run ( self.Nperm , silent = self.silent )            
        else              : eff = permutator     ( self.Nperm , silent = self.silent )
        
        p = eff.eff  
            
        if self.__increasing : p = 1 - p
        
        return t_value , p 
            
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if not np    : logger.warning ('Numpy  is not available') 
    if not sp    : logger.warning ('Scipy  is not available') 
    if not s2u   : logger.warning ('s2u    is not available') 
    if not cdist : logger.warning ('cdist  is not available')
    if not jl    : logger.warning ('joblib is not available')
    
# =============================================================================
##                                                                      The END 
# =============================================================================
