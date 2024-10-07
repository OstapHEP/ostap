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
""" Simple utulities for goodness-of-fit studies for multidimensional fits 
- see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
- see https://doi.org/10.1088/1748-0221/5/09/P09004
- see http://arxiv.org/abs/arXiv:1003.1768
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-29"
__all__     = (
    'PPDnp' , ## Point-to-Point Dissimilarity  Goodness-of-fit method 
    'DNNnp' , ## Distance-to-Nearest-Neighbour Goodness-of-fit method 
)
# =============================================================================
from   ostap.core.ostap_types   import string_types 
from   ostap.stats.gof_utils    import normalize2, PERMUTATOR 
from   ostap.core.core          import SE, VE, Ostap, hID  
from   ostap.utils.progress_bar import progress_bar
from   ostap.utils.utils        import split_n_range
from   ostap.stats.gof          import AGoFnp
import os, abc, warnings, ROOT   
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
## @class GoFnp 
#  A base class for numpy-relarted family of methods to probe goodness-of fit
class GoFnp (AGoFnp) :
    """ A base class for numpy-relarted family of methods to probe goodness-of fit
    """
    def __init__ ( self              ,
                   nToys    = 0      ,
                   silent   = False  , 
                   parallel = True   ) : 
        
        assert isinstance ( nToys , int ) and 0 <= nToys  , \
            "Invalid number of permulations/toys:%s" % nToys
        
        self.__nToys    = nToys
        self.__silent   = True if silent   else False
        self.__parallel = True if parallel else False 
    # ==========================================================================
    ## Normalize two data sets, such that each variable in pooled set
    #  has a mean of zero and variance of one 
    #  @code
    #  ds1 = ... ## the 1st structured array 
    #  ds2 = ... ## the 2nd structured array
    #  gof = ... ##    #  nd1, nd2 = god.normalize ( ds1 , ds2 ) 
    #  @encode
    def normalize ( self , *datasets ) :
        """ Normalize two data sets, such that each variable in pooled set
        has a  mean of zero and variance of one 
        
        >>> ds1 = ... ## the 1st structured array 
        >>> ds2 = ... ## the 2nd structured array
        >>> gof = ... ##
        >>> nd1, nd2 = god.normalize ( ds1 , ds2 ) 
        """ 
        return normalize2  ( datasets  , first = False ) 
    # =========================================================================
    ## Generator of permutations
    #  @code
    #  ds1, ds2 = ...
    #  gof      = ...
    #  for d1,d2 in gof.permutations ( 100 , ds1 , ds2 ) :
    #  ...  
    #  @endcode 
    def permutations ( self , nToys , data1 , data2 ) :
        """ Generator of permutations
        >>> ds1, ds2 = ...
        >>> gof      = ...
        >>> for d1,d2 in gof.permutations ( 100 , ds1 , ds2 ) :
         ...
        """
        n1     = len ( data1 )
        pooled = np.concatenate ( [ data1 , data2 ] )
        ## 
        for i in progress_bar ( nToys , silent = self.silent ) : 
            np.random.shuffle ( pooled )
            yield pooled [ : n1 ] , pooled [ n1 : ]

        del pooled 
    # =========================================================================
    @property 
    def nToys ( self ) :
        """`nToys` : number of permutations/toys used for permutation/toys test"""
        return self.__nToys
    # ========================================================================
    @property
    def silent  ( self ) :
        """`silent` : silent processing?"""
        return self.__silent
    # =======================================================================
    @property
    def parallel ( self ) :
        """`parallel` : parallel processing where/when/if possible?"""
        return self.__parallel
    # =======================================================================
    
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
class PPDnp(GoFnp) : 
    """ Implementation of concrete method "Point-To-Point Dissimilarity"
    for probing of Goodness-Of-Fit
    - see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
    - see https://doi.org/10.1088/1748-0221/5/09/P09004    
    - see http://arxiv.org/abs/arXiv:1003.1768 
    """
    def __init__ ( self                   ,
                   mc2mc     = False      ,
                   nToys     = 1000       ,
                   psi       = 'gaussian' ,
                   sigma     = 0.05       ,
                   parallel  = True       , 
                   silent    = False      ,
                   maxsize   = 1000000    ) :
        
        GoFnp.__init__ ( self                ,
                         nToys    = nToys    ,
                         parallel = parallel , 
                         silent   = silent   )
        
        self.__mc2mc     = True if mc2mc else False
        self.__transform = None
        self.__sigma     = sigma
        self.__psi       = psi
        assert isinstance ( maxsize , int ) and 0 < maxsize , "Invalid ``maxsize'' : %s" % maxsize
        
        self.__maxsize   = max ( maxsize , 100000  )
        
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
        nnmax  = self.__maxsize 
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

        ## For 1D-arrays add a fictive second dimension to please `cdist`-function
        if 1 == uds1.shape [ 1 ] : uds1 = np.c_[ uds1 , np.zeros ( len ( uds1 ) ) ] 
        if 1 == uds2.shape [ 1 ] : uds2 = np.c_[ uds2 , np.zeros ( len ( uds2 ) ) ] 
        
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

        ## For 1D-arrays add a fictive second dimension to please `cdist`-function
        if 1 == uds1.shape [ 1 ] : uds1 = np.c_[ uds1 , np.zeros ( len ( uds1 ) ) ] 
        if 1 == uds2.shape [ 1 ] : uds2 = np.c_[ uds2 , np.zeros ( len ( uds2 ) ) ] 
        
        t_value    = self.t_value ( uds1 , uds2 )

        permutator = PERMUTATOR ( self , t_value , uds1 , uds2 )
        
        if self.parallel and permutator.run :
            counter = permutator.run ( self.nToys , silent = self.silent )            
        else :
            counter = permutator     ( self.nToys , silent = self.silent )

            
        p_value = counter.eff
        
        if self.__increasing : p_value = 1 - p_value
        
        return t_value , p_value 


# =============================================================================
## @class DNNnp
#  Distance-to-Nearest-Neighour GoF-method 
#  @see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
#  @see https://doi.org/10.1088/1748-0221/5/09/P09004
#  @see http://arxiv.org/abs/arXiv:1003.1768 
class DNNnp(GoFnp) : 
    """ Distance-to-Nearest-Neighour GoF-method 
    - see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
    - see https://doi.org/10.1088/1748-0221/5/09/P09004    
    - see http://arxiv.org/abs/arXiv:1003.1768 
    """
    
    def __init__ ( self              ,
                   histo    = None   ,
                   nToys    = 1000   ,
                   parallel = True   , 
                   silent   = False  ) :

        GoFnp.__init__ ( self                ,
                         nToys    = nToys    ,
                         parallel = parallel , 
                         silent   = silent   )
        
        self.__histo = None 
        if   isinstance ( histo , ROOT.TH1 ) :
            self.__histo = histo
        elif isinstance ( histo , int  ) and 1 < histo :
            self.__histo = ROOT.TH1D ( hID() , 'U-values' , histo , 0, 1 )
            
    # =========================================================================
    ## Calculate the t-value
    #  @see Eqs. (3.16) in M.Williams' paper
    #  @param data1 actual data (as unstructured array)
    #  @param vpdf  array of PDF values  
    def t_value ( self , ds1 , vpdf ) :
        """ Calculate the t-value
        - see Eqs. (3.14)&(3.15) in M/.Williams' paper
        data1 : actual data (as unstructured array)
        vpdf  : array of PDF values  
        """
        
        sh1 = ds1 .shape
        sh2 = vpdf.shape
        assert 2 == len ( sh1 ) and 1 == len ( sh2 ) and len ( ds1 ) == len ( vpdf ) , \
            "Invalid arrays: %s , %s" % ( sh1 , sh2 )
        
        tree        = sp.spatial.KDTree ( ds1 )
        conf = { 'k' : [2] }
        if '1.6.0' <= sp.__version__ : conf [ 'workers'] = -1          
        uvalues , _ = tree.query ( ds1 , **conf )

        uvalues     = uvalues.flatten () 

        del tree

        ## dimension of the problem (it must be set in __call__)
        D = self.__D
        C = - Ostap.Math.nball_volume ( D ) * len ( ds1 )  ## constant 
        if 1 != D : uvalues = uvalues ** D
        
        uvalues *= C
        uvalues  = np.exp ( uvalues * vpdf ) 
        
        ## fill the histogram if defined 
        if self.__histo : 
            for u in uvalues : self.__histo.Fill ( float ( u ) ) 
        
        uvalues  = np.sort ( uvalues )

        n        = len ( uvalues )
        aux      = np.linspace ( 1 , n , n ) / n 
        uvalues -= aux

        return   np.sum ( uvalues ** 2 ) 
        
    # ===========================================================================
    ## Calculate the t-value
    #  @see Eqs. (3.16) in M.Williams' paper
    #  @param data1 actual data (as structured array)
    #  @param vpdf  array of PDF values  
    def __call__ ( self ,  data1 , vpdf , normalize = False ) :
        """" Calculate the t-value
        - see Eqs. (3.16) in M.Williams' paper
        data1: actual data (as structured array)
        vpdf : array of PDF values  
        """
        
        ds1 = data1

        ## convert to unstructured dataset 
        uds1    = s2u ( ds1  , copy = False )
        self.__D = uds1.shape[1] 
        
        ## For 1D-arrays add a fictive second dimension to please `cKDTree`-structure
        if 1 == uds1.shape [ 1 ] : uds1 = np.c_[ uds1 , np.zeros ( len ( uds1 ) ) ]
        
        return self.t_value ( uds1 , vpdf  )

    # ============================================================================
    ## p-value is not defined..
    #  M.Williams writes:
    #  ``
    def pvalue ( self , *args , **kwargs ) :
        raise NotImplementedError( "p-value os not defined for DDDNP!" )
    
    @property
    def histo ( self ) :
        """`histo` : the histogram with distribution of U-values"""
        return self.__histo
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    if not np    : logger.warning ( 'Numpy  is not available' ) 
    if not sp    : logger.warning ( 'Scipy  is not available' ) 
    if not s2u   : logger.warning ( 's2u    is not available' ) 
    if not cdist : logger.warning ( 'cdist  is not available' )
 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
