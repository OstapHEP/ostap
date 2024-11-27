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
    'PPD' , ## Point-to-Point Dissimilarity  Goodness-of-fit method 
    'DNN' , ## Distance-to-Nearest-Neighbour Goodness-Of-Fit method 
)
# =============================================================================
from   ostap.stats.gof          import AGoF  
from   ostap.core.core          import Ostap 
from   ostap.fitting.ds2numpy   import ds2numpy
from   ostap.stats.counters     import EffCounter 
from   ostap.utils.progress_bar import progress_bar
from   ostap.utils.utils        import random_name
from   ostap.stats.gof_utils    import TOYS 
import ostap.stats.gof_np       as     GNP
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gofnd' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-fit studies for multidimensional fits' )
# =============================================================================
## @class GoFnD
#  A base class for numpy-relarted family of methods to probe goodness-of fit
class GoF(AGoF) :
    """ A base class for numpy-related family of methods to probe goodness-of fit
    """
    def __init__ ( self             ,
                   gof              , ## actual GoF-evaluator 
                   mcFactor = 10    ,
                   sample   = False ) : 
        
        assert isinstance ( mcFactor, int ) and 1 <= mcFactor , \
            "Invalid `mcFactor':%s" % mcFactor
        self.__gof      = gof
        self.__mcFactor = mcFactor 
        self.__sample   = True if sample else False
        
    @property
    def gof      ( self ) :
        """`gof` : actual Goodness-of-Fit evaluator"""
        return self.__gof 
    
    @property
    def mcFactor ( self ) :
        """`mcFactor` : scale-factor for the size of MC-dataset"""
        return self.__mcFactor
    
    @property
    def nToys ( self )  :
        """`nToys` : number of toys/permutations for toys/permutations test"""
        return self.gof.nToys

    @property
    def sample ( self ) :
        """`sample` : sample numbef of events for generation step?"""
        return self.__sample
    
    @property
    def silent ( self ) :
        """`silent` : silent processing?"""
        return self.gof.silent

    @property
    def parallel ( self ) :
        """`parallel` : parallel processing where/when/if possible?"""
        return self.gof.parallel
    
    # =======================================================================
    ## Generate MC dataset from PDF according to model data
    def generate ( self ,  pdf , data ) :
        """ Generate MC dataset from PDF according to data 
        """
        nEvents = len ( data ) * self.mcFactor
        return pdf.generate ( nEvents = nEvents , varset = data , sample = self.sample )
    
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
        data2 = self.generate ( pdf , data )
        
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
#  - Ntoys    : (int)  number of permutations/toys 
#  - psi      : type of psi/distance function 
#  - sigma    : sigma scale (used for psi=`gaussian`)
class PPD(GoF) : 
    """ Implementation of concrete method "Point-To-Point Dissimilarity" for probing of Goodness-Of-Fit
    - see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
    - see https://doi.org/10.1088/1748-0221/5/09/P09004
    - see http://arxiv.org/abs/arXiv:1003.1768 

    Important parameters:
    
    - mc2mc    : (bool)  should distances within (large) mc-dataste be accounted? 
    - Ntoys   : (int)   number of permutations/toys 
    - psi      : (str)   type of psi/distance function 
    - sigma    : (float) sigma scale (used for psi=`gaussian`) 
    - mcFactor : (int)   the size of mc-dataset is `mcFactor` times size od real data
    """
    # =========================================================================
    ## create the estimator
    #  @param mc2mc    : (bool) should distances within (large) mc-dataste be accounted? 
    #  @param Ntoys    : (int)  number of permutations/toys 
    #  @param psi      : type of psi/distance function 
    #  @param sigma    : sigma scale (used for psi=`gaussian`)
    #  @param mcFactor : (int)  the size of mc-dataset is `mcFactor` times size od real data    
    def __init__ ( self                   ,
                   mc2mc     = False      ,
                   nToys     = 1000       ,
                   psi       = 'gaussian' ,
                   sigma     = 0.05       ,
                   silent    = False      ,
                   parallel  = False      , 
                   mcFactor  = 10         ) : 
        """ Create the Point-to-Point Dissimilarity estimator
        
        Parameters  

        - mcFactor : (int)  the size of mc-dataset is `mcFactor` times size od real data
        - mc2mc    : (bool) should distances within (large) mc-dataste be accounted? 
        - Ntoys    : (int)  number of permutations/toys 
        - psi      : type of psi/distance function 
        - sigma    : sigma scale (used for psi=`gaussian`)
        """
        
        GoF.__init__ ( self ,
                       gof = GNP.PPDnp ( mc2mc    = mc2mc    ,
                                         nToys    = nToys    ,
                                         psi      = psi      ,
                                         sigma    = sigma    ,
                                         parallel = parallel ,
                                         silent   = silent   ) ,                          
                         mcFactor = mcFactor )        
    @property
    def ppd ( self ) :
        """`ppd` : Point-To-Point Dissimilarity calculator for two datasets """
        return self.gof 
    
    # =========================================================================
    ## Calculate T-value for Goodness-of-Git 
    #  @code
    #  ppd    = ...
    #  pdf    = ...  
    #  data   = ... 
    #  tvalue = ppd ( pdf , data ) 
    #  @endcode
    def __call__ ( self , pdf , data ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> ppd    = ...
        >>> pdf    = ... 
        >>> data   = ... 
        >>> tvalue = ppd ( pdf , data ) 
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
class DNN(GoF) : 
    """ Implementation of concrete method "Distance-to-Nearest Neighbour" for probing of Goodness-Of-Fit
    - see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
    - see https://doi.org/10.1088/1748-0221/5/09/P09004
    - see http://arxiv.org/abs/arXiv:1003.1768 
    """

    def __init__ ( self             ,
                   histo    = 100   ,
                   nToys    = 100   ,
                   sample   = False ,
                   parallel = False , 
                   silent   = False ) :
        
        ## initialize the base 
        GoF.__init__ ( self             ,
                       mcFactor = 1     , 
                       gof      = GNP.DNNnp ( histo    = histo    ,
                                              nToys    = nToys    ,
                                              parallel = parallel , 
                                              silent   = silent ) ,
                       sample   = sample )
        
        self.__histo  = None 
        
    @property
    def dnn ( self ) :
        """`ppd` : Distance-to-Nearest-Neighbour calculator for numpy sata"""
        return self.gof  

    # =========================================================================
    ## Transform pdf&data
    def transform ( self , pdf , data ) :
        """ Transform pdf&data """

        vs1  = pdf .vars
        vs2  = data.get()
        vlst = set() 
        for v in vs1 :
            if v in vs2 : vlst.add ( v.name )
        vlst = tuple ( vlst ) 
        
        pdf_name = random_name ( prefix = 'pdf_' , suffix = '_value' )
        while pdf_name in data : pdf_name = random_name ( prefix = 'pdf_' , suffix = '_value' )            

        ## convert the input data into structured array with the additional column with PDF-values 
        ds1  = ds2numpy ( data , vlst , more_vars = { pdf_name : pdf } )
        
        ## split structured array into "normal" array of PDF-values & (structured) array of data 
        vpdf = ds1 [ pdf_name      ]         
        ds1  = ds1 [ list ( vlst ) ]

        return ds1 , vpdf 
        
    # =========================================================================
    ## Calculate T-value for Goodness-of-Git 
    #  @code
    #  dnn     = ...
    #  pdf     = ...  
    #  data    = ... 
    #  tvalue  = dnn ( pdf , data ) 
    #  @endcode
    def __call__ ( self , pdf , data ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> dnn    = ...
        >>> pdf    = ... 
        >>> data   = ... 
        >>> tvalue = dnn ( pdf , data ) 
        """

        ds1 , vpdf = self.transform ( pdf , data ) 

        ## call underlying method 
        return self.dnn ( ds1 , vpdf , normalize = True )

    # =========================================================================
    ## Calculate the t & p-values
    #  @code
    #  dnn  = ...
    #  pdf  = ...
    #  data = ... 
    #  t , p = dnn.pvalue ( pdf , data )
    #  @endcode 
    def pvalue ( self , pdf , data ) :
        """ Calculate the t & p-values
        >>> dnn  = ...
        >>> pdf  = ... 
        >>> data = ... 
        >>> t , p = dnn.pvalue ( pdf , data ) 
        """

        ds1 , vpdf = self.transform ( pdf , data ) 

        ## call underlying method 
        t_value = self.dnn ( ds1 , vpdf )

        ## save the histogram
        if   self.__histo   : pass
        elif self.dnn.histo :
            self.__histo = self.dnn.histo.clone() 
            self.dnn.histo.Reset()

        ## prepare toys
        toys = TOYS ( self , t_value , pdf = pdf , Ndata = len ( data ) , sample = self.sample )
        
        if self.parallel : counter = toys.run ( self.nToys , silent = self.silent )            
        else             : counter = toys     ( self.nToys , silent = self.silent )            
        
        p_value = 1 - counter.eff
        return t_value, p_value 

    @property
    def histo ( self ) :
        """`histo` : histogram with U-value distribution"""
        return self.__histo
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
