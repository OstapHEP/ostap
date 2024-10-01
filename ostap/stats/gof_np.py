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
try : # =======================================================================
    # =========================================================================
    import numpy                                                    as np
    from   numpy.lib.recfunctions import structured_to_unstructured as s2u
    # =========================================================================
except ImportError :
    # =========================================================================
    np    = None
    s2u   = None
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import scipy                                                    as sp 
    from   scipy.spatial.distance import cdist                      as cdist 
    # =========================================================================
except ImportError :
    # =========================================================================
    sp    = None
    cdist = None 
# =============================================================================
from   ostap.core.ostap_types   import string_types 
from   ostap.stats.gof          import normalize2 
from   ostap.core.core          import SE, VE
from   ostap.utils.progress_bar import progress_bar 
import abc 
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
    #  for d1,d2 in gof.permulations ( 100 , ds1 , ds2 ) :
    #  ...  
    #  @endcode 
    def permutations ( self , data1 , data2 ) :
        """ Generator of permutations
        >>> ds1, ds2 = ...
        >>> gof      = ...
        >>> for d1,d2 in gof.permulations ( 100 , ds1 , ds2 ) :
         ...
        """
        n1     = len ( data1 )
        pooled = np.concatenate ( [ data1 , data2 ] )

        for i in progress_bar ( self.Nperm , silent = self.silent ) : 
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
                   Nperm     = 100        ,
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

        n1 = len ( ds1 )
        n2 = len ( ds2 )
        
        ## convert to unstructured datasets 
        uds1  = s2u ( data1 , copy = False )
        uds2  = s2u ( data2 , copy = False )
        
        ## how to build distances?
        scale = -0.5/(self.sigma**2) 
        distance_type , transform , _ = psi_conf ( self.psi , scale )
        
        ## list of all distances between points in the first dataset 
        dist_11 = cdist ( uds1 , uds1 , distance_type ) .flatten () ## data <-> data
        dist_11 = dist_11 [ dist_11 > 0 ] 
        if transform : dist_11  = transform ( dist_11 )

        ## result = np.sum ( dist_11 ) / ( n1 * ( n1 - 1 ) )
        result = np.mean ( dist_11 ) 
        del dist_11
        
        ## list of all distances between points in the 1st and 2nd datasets  
        dist_12 = cdist ( uds1 , uds2 , distance_type ) .flatten ()  ## data <-> mc
        dist_12 = dist_12 [ dist_12 > 0 ] 
        if transform : dist_12  = transform ( dist_12 )

        ## result -= np.sum ( dist_12 ) / ( n1 * n2 )
        result -= np.mean ( dist_12 )
        del dist_12 

        ## add the distances from the second dataset? 
        if self.mc2mc :
            ## list of all distances between points in the 2nd dataset 
            dist_22 = cdist ( uds2 , uds2 , distance_type ) . flatten() ## mc<-> mc 
            dist_22 = dist_22 [ dist_22 > 0 ] 
            if transform : dist_22  = transform ( dist_22 )
            ## result += np.sum ( dist_22 ) / ( n2 * ( n2 - 1 ) )
            result += np.mean ( dist_22 )
            del dist_22

        return result    
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

        ## calculate t-value
        tvalue = self 
        t = tvalue ( ds1 , ds2 , normalize = False )
        
        pv = SE () 
        ## start the permutation test
        for d1 , d2 in self.permutations ( ds1 , ds2 ) :
            ti  = tvalue ( d1 , d2 , normalize = False )
            pv += float  ( t < ti ) 

        p = VE ( pv.eff () , pv.effErr() ** 2 )
        
        if self.__increasing : p = 1 - p

        return t , p 
            
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if not np    : logger.warning ('Numpy is not available') 
    if not sp    : logger.warning ('Scipy is not available') 
    if not s2u   : logger.warning ('s2u   is not available') 
    if not cdist : logger.warning ('cdist is not available') 
# =============================================================================
##                                                                      The END 
# =============================================================================
