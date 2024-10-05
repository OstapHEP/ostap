#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof_1d.py
#  Set of utulities for goodness-of-fit studies for multidimensional fits
#  @see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
#  @see https://doi.org/10.1088/1748-0221/5/09/P09004
#  @see http://arxiv.org/abs/arXiv:1003.1768
# 
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
    'PPD' , ## Point-to-Point Dissimilarity Goodness-of-fit method 
)
# =============================================================================
from   ostap.core.core        import Ostap 
from   ostap.stats.gof_np     import PPDNP 
from   ostap.fitting.ds2numpy import ds2numpy 
import ROOT, abc 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gofnd' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-fit studies for multidimensional fits' )
# =============================================================================
## @class AGoF
#  An absract base class for family of methods to probe goodness-of fit
#  There are two abstract methods
#  - <code>__call__</code> to evaluate t-value 
#  - <code>pvalue</code> to evaluate (t,p)-values for two datasets
#  @code
#  gof  = ...
#  pdf  = ...
#  data = ...
#  pdf.fitTo ( data ) 
#  t    = gof        ( pdf , data )
#  t,p  = gof.pvalue ( pdf , ds2 )
#  @endcode
class AGoF(object) :
    """ An absract base class for numpy-related family of methods to probe goodness-of fit
    There are two abstract methods
    - `__call__` to evaluate t-value for the fit 
    - `pvalue` to evaluate (t,p)-vaues for two datasets
    
    >>> gof  = ...
    >>> pdf  = 
    >>> data = ...
    >>> pdf.fitTo ( data , ... ) 
    >>> t    = gof        ( pdf , data )
    >>> t ,p = gof.pvalue ( pdf , data )
    """
    # =========================================================================
    ## Calculate T-value for Goodness-of-Git 
    #  @code
    #  gof   = ...
    #  pdf   = ...  
    #  data  = ... 
    #  t = gof ( pdf , data ) 
    #  @endcode
    @abc.abstractmethod
    def __call__ ( self , pdf , data ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> gof   = ...
        >>> pdf   = ... 
        >>> data  = ... 
        >>> t = gof ( pdf , data ) 
        """
        return NotImplemented 
    # =========================================================================
    ## Calculate the t & p-values
    #  @code
    #  gof  = ...
    #  pdf  = ...
    #  data = ... 
    #  t , p = gof.pvalue ( pdf , data )
    #  @endcode 
    @abc.abstractmethod
    def pvalue ( self , pdf , data ) :
        """ Calculate the t & p-values
        >>> gof  = ...
        >>> pdf  = ... 
        >>> data = ... 
        >>> t , p = gof.pvalue ( pdf , data ) 
        """
        return NotImplemented

# =================================================================================
## @class GoFnD
#  A base class for numpy-relarted family of methods to probe goodness-of fit
class GoFnD (object) :
    """ A base class for numpy-relarted family of methods to probe goodness-of fit
    """
    def __init__ ( self , mcFactor = 10 ) :

        assert isinstance ( mcFactor, int ) and 1 <= mcFactor , \
            "Invalid `mcFactor':%s" % mcFactor

        self.__mcFactor = mcFactor 

    @property
    def mcFactor ( self ) :
        """`mcFactor` : scale-factor for the size of MC-dataset"""
        return self.__mcFactor

    # ==========================================================================
    ## Generate MC dataset from PDF according to model data
    def generate ( self ,  pdf , data , sample = False ) :
        """ Generate MC dataset from PDF according to data 
        """
        nEvents = len ( data ) * self.mcFactor
        return pdf.generate ( nEvents = nEvents , varset = data , sample = sample )
    
    # ==========================================================================
    ## Transform a (pdf,data) pair into (data_np, mc_np) pair 
    #  @code
    #  gof  = ...
    #  pdf  = ...
    #  data = ... ## as ROOT.RooAbsData 
    #  ds1 , ds2 = gof.transform ( pdf , data ) 
    #  @endcode
    def transform  ( self , pdf , data ) :
        """ Transform a (pdf,data) pair into (data_np, mc_np) pair 
        >>> gof  = ...
        >>> pdf  = ...
        >>> data = ... ## as ROOT.RooAbsData
        >>> ds1, ds2 = gof.transform ( pdf , data ) 
        """
        
        data1 = data
        data2 = self.generate ( pdf , data , sample = False )
        
        vs1  = data1.get()
        vs2  = data2.get()
        vlst = set() 
        for v in vs1 :
            if v in vs2 : vlst.add ( v.name )
        vlst = tuple ( vlst ) 
        
        ds1 = ds2numpy ( data1 , vlst )
        ds2 = ds2numpy ( data2 , vlst )
        
        ## delete 
        if isinstance ( data2  , ROOT.RooDataSet ) :
            data2 = Ostap.MoreRooFit.delete_data ( data2 )
            
        del data2
        
        return ds1, ds2
                   
# =============================================================================
## @class PPD
#  Implementation of concrete method "Point-To-Point Dissimilarity"
#  for probing of Goodness-Of-Fit
#  @see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
#  @see https://doi.org/10.1088/1748-0221/5/09/P09004
#  @see http://arxiv.org/abs/arXiv:1003.1768 
#  Important parameters:
#  - mcFactor : (int)  the size of mc-dataset is `mcFactor` times size od real data
#  - mc2mc    : (bool) should distances witng (large) mc-dataste be accounted? 
#  - Nperm    : (int)  number of permutations 
#  - psi      : type of psi/distance function 
#  - sigma    : sigma scale (used for psi=`gaussian`)
class PPD(AGoF,GoFnD) : 
    """ Implementation of concrete method "Point-To-Point Dissimilarity" for probing of Goodness-Of-Fit
    - see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
    - see https://doi.org/10.1088/1748-0221/5/09/P09004
    - see http://arxiv.org/abs/arXiv:1003.1768 

    Important parameters:
    
    - mc2mc    : (bool)  should distances within (large) mc-dataste be accounted? 
    - Nperm    : (int)   number of permutations 
    - psi      : (str)   type of psi/distance function 
    - sigma    : (float) sigma scale (used for psi=`gaussian`) 
    - mcFactor : (int)   the size of mc-dataset is `mcFactor` times size od real data
    """
    # =========================================================================
    ## create the estimator
    #  @param mc2mc    : (bool) should distances within (large) mc-dataste be accounted? 
    #  @param Nperm    : (int)  number of permutations 
    #  @param psi      : type of psi/distance function 
    #  @param sigma    : sigma scale (used for psi=`gaussian`)
    #  @param mcFactor : (int)  the size of mc-dataset is `mcFactor` times size od real data    
    def __init__ ( self                   ,
                   mc2mc     = False      ,
                   Nperm     = 1000       ,
                   psi       = 'gaussian' ,
                   sigma     = 0.05       ,
                   mcFactor  = 10         ) : 
        """ Create the Point-to-Point Dssimilaroity estimator
        
        Parameters  

        - mcFactor : (int)  the size of mc-dataset is `mcFactor` times size od real data
        - mc2mc    : (bool) should distances within (large) mc-dataste be accounted? 
        - Nperm    : (int)  number of permutations 
        - psi      : type of psi/distance function 
        - sigma    : sigma scale (used for psi=`gaussian`)
        """
        ## initialize the base 
        GoFnD.__init__ ( self , mcFactor = mcFactor )

        ## the actual worker: estimator for two datasets 
        self.__ppd = PPDNP ( mc2mc = mc2mc ,
                             Nperm = Nperm ,
                             psi   = psi   ,
                             sigma = sigma )
        
    @property
    def ppd ( self ) :
        """`ppd` : Point-To-Point Dissimilarity calculator for two datasets """
        return self.__ppd 

    # =========================================================================
    ## Calculate T-value for Goodness-of-Git 
    #  @code
    #  ppd   = ...
    #  pdf   = ...  
    #  data  = ... 
    #  t = ppd ( pdf , data ) 
    #  @endcode
    def __call__ ( self , pdf , data ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> ppd   = ...
        >>> pdf   = ... 
        >>> data  = ... 
        >>> t = ppd ( pdf , data ) 
        """

        ds1, ds2 = self.transform ( pdf , data ) 
        
        ## estimate t-value 
        return self.ppd ( ds1 , ds2 , normalize = True )

    # =========================================================================
    ## Calculate the t & p-values
    #  @code
    #  ppd  = ...
    #  pdf  = ...
    #  data = ... 
    #  t , p = ppd.pvalue ( pdf , data )
    #  @endcode 
    def pvalue ( self , pdf , data ) :
        """ Calculate the t & p-values
        >>> ppd  = ...
        >>> pdf  = ... 
        >>> data = ... 
        >>> t , p = ppd.pvalue ( pdf , data ) 
        """

        ds1, ds2 = self.transform ( pdf , data ) 
        
        ## estimate t&p-values 
        return self.ppd.pvalue ( ds1 , ds2 , normalize = True )
    # =========================================================================
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
