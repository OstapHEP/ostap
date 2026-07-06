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
    'PPD'               , ## Point-to-Point Dissimilarity  Goodness-of-fit method 
    'DNN'               , ## Distance-to-Nearest-Neighbour Goodness-Of-Fit method
    'USTAT'             , ## Alternative implementation of DNN method
    'NLL'               , ## Use -log L as GoF estimator 
    'AikaikeIC'         , ## Use Aikaike IC  as GoF estimator 
    'BayesianIC'        , ## Use Bayesian IC  as GoF estimator
    'ADVAL_LightGBM'    , ## Use Adversarial Validation as GoF estimator 
    'ADVAL_HistoGBoost' , ## Use Adversarial Validation as GoF estimator 
    'ADVAL_GBoost'      , ## Use Adversarial Validation as GoF estimator 
    'ADVAL_XGBoost'     , ## Use Adversarial Validation as GoF estimator 
)
# =============================================================================
from   ostap.core.ostap_types   import num_types, integer_types, sized_types
from   ostap.core.cpu_info      import HAS_AVX2 
from   ostap.utils.basic        import typename  
from   ostap.stats.gof          import AGoF
from   ostap.core.core          import Ostap, VE 
from   ostap.math.math_base     import axis_range
from   ostap.fitting.ds2numpy   import ds2numpy
from   ostap.stats.counters     import EffCounter 
from   ostap.utils.progress_bar import progress_bar
from   ostap.utils.utils        import random_name
from   ostap.stats.gof_utils    import TOYS
from   ostap.stats.ustat        import USTAT
from   ostap.plotting.color     import Navy, DarkGreen
from   ostap.stats.gof_utils    import format_row, draw_ecdf  
import ostap.stats.gof_np       as     GNP
import ostap.logger.table       as     T 
import ROOT, numpy 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gofnd' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-fit studies for multidimensional fits' )
# =============================================================================
## @class GoF
#  A base class for family of methods to probe goodness-of-fit
class GoF(AGoF) :
    """ A base class for family of methods to probe goodness-of-fit
    """
    def __init__ ( self             ,
                   gof              , ## actual GoF-evaluator
                   mcFactor = 20    ,
                   sample   = False ) : 
        
        assert isinstance ( mcFactor, int ) and 1 <= mcFactor , "Invalid `mcFactor':%s" % mcFactor
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
        """`nToys` : number of toys/permutations for toys (or permutations) test"""
        return self.gof.nToys

    @property
    def sample ( self ) :
        """`sample` : sample number of events for generation step?"""
        return self.__sample
    
    @property
    def silent ( self ) :
        """`silent` : silent processing?"""
        return self.gof.silent

    @property
    def progress ( self ) :
        """`progress` : show progress bar?"""
        return self.gof.progress 
    
    @property
    def parallel ( self ) :
        """`parallel` : parallel processing where/when/if possible?"""
        return self.gof.parallel
    
    @property
    def method ( self ) :
        """`method` : the actual GoF-method 
        """
        return self.gof.method
        
    # =========================================================================
    ## Calculate t-value for Goodness-of-Fit test
    #  @code
    #  gof   = ...
    #  pdf   = ...  
    #  data  = ... 
    #  t_value = gof ( pdf , data ) 
    #  @endcode
    def __call__ ( self , pdf , data ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> gof   = ...
        >>> pdf   = ... 
        >>> data  = ... 
        >>> t_value = gof ( pdf , data ) 
        """
        return self.tvalue ( pdf , data ) 

    # =======================================================================
    ## Generate MC dataset from PDF according to model data
    def generate ( self ,  pdf , data ) :
        """ Generate MC dataset from PDF according to data 
        """
        nEvents = len ( data ) * self.mcFactor
        dset    = pdf.generate ( nEvents = nEvents , varset = data , sample = self.sample )        
        assert dset or self.sample , "generate: failrue to produce non-empty dataset!"
        while not dset : dset = pdf.generate ( nEvents = nEvents , varset = data , sample = self.sample )        
        return dset 
            
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
        if isinstance ( data2  , ROOT.RooDataSet ) : data2.clear()            
        del data2 
        
        return ds1 , ds2

    # =========================================================================
    ## Transform a (pdf,wdata) pair into (data_np, mc_np, w_np) triplet 
    #  @code
    #  gof  = ...
    #  pdf  = ...
    #  data = ... ## as *WEIGHTED* ROOT.RooAbsData 
    #  ds1 , ds2 , w = gof.wtransform ( pdf , data ) 
    #  @endcode
    def wtransform  ( self , pdf , data ) :
        """ Transform a (pdf,data) pair into (data_np, mc_np) pair 
        >>> gof  = ...
        >>> pdf  = ...
        >>> data = ... ## as ROOT.RooAbsData
        >>> ds1, ds2 , w1 = gof.wtransform ( pdf , data ) 
        """

        if not data.isWeighted() :
            ds1 , ds2 = self.transform ( pdf , data )
            w1  = None 
            return ds1 , ds2 , w1 

        
        data1 = data
        data2 = self.generate ( pdf , data )

        vars1 = data1.get() 
        vars2 = data2.get()
        
        wname = data1.wname
        
        var_lst = tuple ( sorted ( v.name for v in var1  if v in var2 and v.name != wname ) )  
                   
        ds1 , w1 = data1.slice ( var_lst , silent = silent , structured = True , weight_name = data.wname() )
        ds2 , _  = data2.slice ( var_lst , silent = silent , structured = True )

        ## scale the weights such that
        if w1 is None              : pass
        elif numpy.all ( w1 == 1 ) : pass 
        else                       :
            wsum = numpy.sum ( w1 )
            assert 0 < wsum , "Sum of weights is non-positive: %g" % wsum
            w1   = w1 * ( len ( w1 ) * 1.0 / wsum )
            
        ## delete 
        if isinstance ( data2  , ROOT.RooDataSet ) : data2.clear()            
        del data2 
            
        return ds1 , ds2 , w1 

    # =========================================================================
    ## Draw the empirical CDF from permutations or toys  
    def draw  ( self , tvalue = None , option = '' , *args , **kwargs ) :
        """ Draw empirical CDF from permutations or toys 
        """
        ## 
        ecdf = self.ecdf 
        if not ecdf : return ecdf 
        ## 
        has_tvalue = not tvalue is None and isinstance ( tvalue , num_types )
        ## 
        if not has_tvalue : return draw_ecdf (  ecdf , tvalue = None  , option = option , **kwargs )
        ## 
        result, vline, hline = draw_ecdf (  ecdf , tvalue = tvalue , option = option , **kwargs )
        ## 
        self._vline = vline 
        self._hline = hline 
        ##
        return result, vline, hline   

# =============================================================================
## @class PPD
#  Implementation of concrete method "Point-To-Point Dissimilarity"
#  for probing of Goodness-Of-Fit
#  @see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
#  @see https://doi.org/10.1088/1748-0221/5/09/P09004
#  @see http://arxiv.org/abs/arXiv:1003.1768 
#  Important parameters:
#  - mcFactor : (int)  the size of mc-dataset is `mcFactor` times size of real data
#  - mc2mc    : (bool) should distances witng (large) mc-dataste be accounted? 
#  - nToys    : (int)  number of permutations/toys 
#  - psi      : type of psi/distance function 
#  - sigma    : sigma scale (used for psi=`gaussian`)
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
class PPD(GoF) : 
    """ Implementation of concrete method "Point-To-Point Dissimilarity" for probing of Goodness-Of-Fit
    - see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
    - see https://doi.org/10.1088/1748-0221/5/09/P09004
    - see http://arxiv.org/abs/arXiv:1003.1768 

    Important parameters:
    
    - mc2mc    : (bool)  should distances within (large) mc-dataste be accounted? 
    - nToys    : (int)   number of permutations/toys 
    - psi      : (str)   type of psi/distance function 
    - sigma    : (float) sigma scale (used for psi=`gaussian`) 
    - mcFactor : (int)   the size of mc-dataset is `mcFactor` times size of real data

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
    # =========================================================================
    ## create the estimator
    #  @param mc2mc    : (bool) should distances within (large) mc-dataset be accounted for? 
    #  @param nToys    : (int)  number of permutations/toys 
    #  @param psi      : type of psi/distance function 
    #  @param sigma    : sigma scale (used for psi=`gaussian`)
    #  @param mcFactor : (int)  the size of mc-dataset is `mcFactor` times size of real data    
    def __init__ ( self                   ,
                   mc2mc     = False      ,
                   nToys     = 1000       ,
                   psi       = 'gaussian' ,
                   sigma     = 0.10       ,
                   silent    = False      ,
                   parallel  = False      ,
                   mcFactor  = 10         ) : 
        """ Create the Point-to-Point Dissimilarity estimator
        
        Parameters  

        - mcFactor : (int)  the size of mc-dataset is `mcFactor` times size of real data
        - mc2mc    : (bool) should distances within (large) mc-dataset be accounted for? 
        - nToys    : (int)  number of permutations/toys 
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

        self.__ecdf   = None
        self.__tvalue = None
        self.__pvalue = None
        
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
    def tvalue ( self , pdf , data ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> ppd    = ...
        >>> pdf    = ... 
        >>> data   = ... 
        >>> tvalue = ppd.tvalue ( pdf , data ) 
        """
        ds1, ds2 = self.transform ( pdf , data )         
        ## estimate the t-value 
        t = self.ppd ( ds1 , ds2 , normalize = True )
        ## 
        self.__tvalue = float ( t ) 
        return self.__tvalue 

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
        return self.tvalue ( pdf , data ) 
    
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
        self.__tvalue , self.__pvalue = self.ppd.pvalue ( ds1 , ds2 , normalize = True )        
        return self.__tvalue , self.__pvalue

    # =========================================================================
    @property
    def ecdf ( self ) :
        """`ecdf` : empirical CDF for t-value distribution"""
        return self.ppd.ecdf

    # =========================================================================
    ## Get results in a form of the table 
    def table ( self , title = '' , prefix = '' ) :
        """ Get results in a for of the table 
        """
        return self.the_table ( tvalue = self.__tvalue ,
                                pvalue = self.__pvalue ,
                                ecdf   = self.__ecdf   , 
                                title  = title ,
                                prefix = prefx ) 
    __str__  = table
    __repr__ = table
    
# =============================================================================
## @class DNN
#  Distance-to-Nearest-Neighour GoF-method 
#  @see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
#  @see https://doi.org/10.1088/1748-0221/5/09/P09004
#  @see http://arxiv.org/abs/arXiv:1003.1768
#
#  - M.Williams writes: 
#    The method <..> is easy to use, requires very little processing time and is
#    conceptually fairly easy to understand; however it is not very powerful.
#    The U-statistic it defines does provide a useful easy-to-visualize diagnostic
#    tool (especially for very high dimensional analyses), but its quantitative
#    usefulness as a g.o.f. test is limited.  
class DNN(GoF) : 
    """ Implementation of concrete method "Distance-to-Nearest Neighbour" for probing of Goodness-Of-Fit
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
    def __init__ ( self             ,
                   histo    = 100   ,
                   nToys    = 500   ,
                   sample   = False ,
                   parallel = False , 
                   silent   = False ,
                   progress = True  ) : 
        
        ## initialize the base 
        GoF.__init__ ( self             ,
                       mcFactor = 1     , 
                       gof      = GNP.DNNnp ( histo    = histo    ,
                                              nToys    = nToys    ,
                                              parallel = parallel , 
                                              silent   = silent   ,
                                              progress = progress ) ,
                       sample   = sample )
        
        self.__ecdf   = None 
        self.__histo  = None 
        self.__tvalue = None
        self.__pvalue = None
        
    @property
    def dnn ( self ) :
        """`dnn` : Distance-to-Nearest-Neighbour calculator for numpy data
        """
        return self.gof  

    # =========================================================================
    ## (internal method) Transform pdf&data
    def transform ( self , pdf , data ) :
        """ (internal method) Transform pdf&data
        """

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
    def tvalue ( self , pdf , data ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> dnn    = ...
        >>> pdf    = ... 
        >>> data   = ... 
        >>> tvalue = dnn ( pdf , data ) 
        """
        ds1 , vpdf = self.transform ( pdf , data ) 
        ## call underlying method 
        t  = self.dnn ( ds1 , vpdf , normalize = False )
        self.__tvalue = float ( t ) 
        return self.__tvalue 

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

        ## save the histogram (do not override it!)
        if   self.__histo   : pass
        elif self.dnn.histo :
            self.__histo = self.dnn.histo.clone() 
            self.dnn.histo.Reset()

        ## prepare toys
        toys = TOYS ( self                   ,
                      t_value = t_value      ,
                      pdf     = pdf          ,
                      Ndata   = len ( data ) ,
                      sample  = self.sample  )

        if self.parallel : counter = toys.run ( self.nToys , progress = self.progress , silent = self.silent ) 
        else             : counter = toys     ( self.nToys , progress = self.progress )
        
        ## get ECDF from toys
        self.__ecdf = toys.ecdf
        
        p_value       = 1 - counter.eff
        
        self.__tvalue = t_value 
        self.__pvalue = p_value 

        return t_value, p_value
    
    @property
    def histo ( self ) :
        """`histo` : histogram with U-value distribution"""
        return self.__histo
    
    @property
    def ecdf ( self ) :
        """`ecdf` : empirical CDF for t-value distribution"""
        return self.__ecdf

    # =========================================================================
    ## Get results in a form of the table 
    def table ( self , title = '' , prefix = '' ) :
        """ Get results in a form of the table 
        """
        return self.the_table ( tvalue = self.__tvalue ,
                                pvalue = self.__pvalue ,
                                ecdf   = self.ecdf     , 
                                title  = title  ,
                                prefix = prefix ) 
    __str__  = table
    __repr__ = table

    
# =============================================================================
## @class NLL
#  Trivial `estimator' for -log N as `fit-quality'
#  - @attention the absolute value of -log L is *not* a true GoF-estimator
#  toys are needed!
class NLL(AGoF) :
    """  Trivial `estimator' for -log N as fit-quality.
    - attention the absolute value of -log L is not a true GoF-estimator
    - toys are needed!
    """
    def __init__ ( self             ,
                   fitresult        , * , 
                   fitconf  = { 'silent' : True  , 'draw' : False } ,
                   nToys    = 100   , 
                   sample   = False ,
                   parallel = False ,
                   silent   = True  ,
                   progress = True  ) :
        
        assert fitresult and isinstance ( fitresult , ROOT.RooFitResult ) , \
            "Invalid fit-result %s" % typename ( fitresult ) 

        assert isinstance ( nToys , integer_types ) and 1 <= nToys , \
            "Invalid nToys %s" % typename ( nToys )

        self.__sample    = True if sample   else False 
        self.__parallel  = True if parallel else False 
        self.__silent    = True if silent   else False 
        self.__progress  = True if progress else False 
        self.__fitconf   = fitconf
        self.__nToys     = nToys 
        self.__fitresult = fitresult 
        self.__ecdf      = None
    
    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object
        """
        return { 'sample'    : self.sample    ,
                 'parallel'  : self.parallel  ,
                 'silent'    : self.silent    ,
                 'progress'  : self.progress  ,
                 'fitconf'   : self.fitconf   ,
                 'fitresult' : self.fitresult ,
                 'nToys'     : self.nToys     ,
                 'ecdf'      : self.ecdf      }
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object
        """
        self.__sample    = state.pop ( 'sample'    )
        self.__parallel  = state.pop ( 'parallel'  )
        self.__silent    = state.pop ( 'silent'    )
        self.__progress  = state.pop ( 'progress'  )        
        self.__fitconf   = state.pop ( 'fitconf'   )
        self.__fitresult = state.pop ( 'fitresult' )        
        self.__nToys     = state.pop ( 'nToys'     )
        self.__ecdf      = state.pop ( 'ecdf'      )
    
    # =========================================================================
    @property
    def sample    ( self ) :
        """`sample` : sample number of events for toys?"""
        return self.__sample
    
    @property
    def silent    ( self ) :
        """`silent` : silent prpcessing ?"""
        return self.__silent
    
    @property
    def progress    ( self ) :
        """`progress` : show progress bar?"""
        return self.__progress
    
    @property
    def parallel  ( self ) :
        """`parallel` : parallel processing?"""
        return self.__parallel
    
    @property
    def fitconf   ( self ) :
        """`fitconf` : configuraton of fit (arguments for `PDF.fitTo`)"""
        return self.__fitconf
    
    @property
    def fitresult ( self ) :
        """`fitresult` : initially specified fit result"""
        return self.__fitresult
        
    @property
    def nToys     ( self ) :
        """`nToys: number of pseudoexperiments for p-value"""
        return self.__nToys  

    @property
    def ecdf ( self ) :
        """`ecdf: empirical CDF for t-values"""
        return self.__ecdf 

    # =========================================================================
    ## Calculate t-value for Goodness-of-Fit test
    #  @code
    #  gof   = ...
    #  pdf   = ...  
    #  data  = ... 
    #  t_value = gof ( pdf , data ) 
    #  @endcode
    def __call__ ( self , pdf , data ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> gof   = ...
        >>> pdf   = ... 
        >>> data  = ... 
        >>> t_value = gof ( pdf , data ) 
        """
        return self.tvalue ( pdf , data ) 
    
    # =========================================================================
    ## get the actual t-value from the fit-result 
    def t_value   ( self , fitresult ) :
        return fitresult.minNll () 
        
    # ===========================================================================
    ## get the t-value 
    def tvalue ( self, pdf , data ) :
        """ Get t-value: -log L here 
        """
        ## to ensure consistency:
        pdf.load_params ( self.fitresult , silent = True )

        r , _ = pdf.fitTo ( data , **self.fitconf )
        
        tv = self.t_value ( r ) 
        del r
        
        return tv

    # =========================================================================
    ## Calculate the t & p-values
    #  @code
    #  gof  = ...
    #  pdf  = ...
    #  data = ... 
    #  t_value , p_value = gof.pvalue ( pdf , data )
    #  @endcode 
    def pvalue ( self , pdf , data ) :
        """ Calculate the t & p-values
        >>> gof  = ...
        >>> pdf  = ... 
        >>> data = ... 
        >>> t_value , p_value = gof.pvalue ( pdf , data ) 
        """

        t_value   = self.t_value ( self.fitresult ) 

        ## prepare toys
        toys = TOYS ( self                        ,
                      t_value    = t_value        ,
                      pdf        = pdf            ,
                      Ndata      = len ( data )   ,
                      sample     = self.sample    ,
                      parameters = self.fitresult )
        

        if self.parallel : counter = toys.run ( self.nToys , progress = self.progress , silent = self.silent ) 
        else             : counter = toys     ( self.nToys , progress = self.progress )
                 
        ## get ECDF from toys
        self.__ecdf = toys.ecdf
        
        p_value       = 1 - counter.eff
        
        self.__tvalue = t_value 
        self.__pvalue = p_value 

        return t_value, p_value
    
    # =========================================================================
    ## Draw the empirical CDF from toys  
    def draw  ( self , option = '' , *args , **kwargs ) :
        """ Draw empirical CDF from toys 
        """
        ## 
        ecdf = self.ecdf 
        if not ecdf : return ecdf 
        ##
        t_value = self.t_value ( self.fitresult ) 
        ## 
        result, vline, hline = draw_ecdf (  ecdf , tvalue = t_value , option = option , **kwargs )
        ## 
        self._vline = vline 
        self._hline = hline 
        ##
        ecdf.__lines = vline , hline 
        ## 
        return result, vline, hline   

# =============================================================================
## @class AikaikeIC
#  Use Aikaike information criterion as `estimate` for Goodness-of-Fit
class AikaikeIC(NLL) :
    """  Trivial `estimator' for Aikaike informatio criterion as fit-quality.
    - attention: the absolute value of -log L is not a true GoF-estimator
    - toys are needed!
    """
    # =========================================================================
    ## get the actual t-value from the fit-result 
    def t_value   ( self , fitresult ) : return fitresult.aic () 

# =============================================================================
## @class BayesianIC
#  Use Bayesian information criterion as `estimate` for Goodness-of-Fit
class BayesianIC(NLL) :
    """  Trivial `estimator' for Bayesian informatio criterion as fit-quality.
    - attention: the absolute value of -log L is not a true GoF-estimator
    - toys are needed!
    """
    # =========================================================================
    def __init__ ( self             ,
                   fitresult        ,
                   data             , * , 
                   fitconf  = { 'silent' : True  , 'draw' : False } ,
                   nToys    = 100   , 
                   sample   = False ,
                   parallel = False ,
                   silent   = False ,
                   progress = True  ) :

        if   isinstance ( data , ROOT.RooDataSet ) : data = len ( data )
        elif isinstance ( data , sized_types     ) : data = len ( data )
        
        assert data and isinstance ( data , integer_types ) and 1 <= data , \
            "Invalid data %s" % typename ( data )
        
        NLL.__init__ ( self      ,
                       fitresult = fitresult , 
                       fitconf   = fitconf   ,
                       nToys     = nToys     ,  
                       sample    = sample    , 
                       parallel  = parallel  , 
                       silent    = silent    ,
                       progress  = progress  ) 

        self.__ndata = data 

    ## get the actual t-value from the fit-result 
    def t_value   ( self , fitresult ) : return fitresult.bic ( self.__ndata ) 

    # ===========================================================================
    ## get the t-value 
    def tvalue ( self, pdf , data ) :
        """ Get t-value: -log L here 
        """
        ## to ensure consistency:
        pdf.load_params ( self.fitresult , silent = True )

        r , _ = pdf.fitTo ( data , **self.fitconf )
        
        tv = r.bic ( data ) 
        del r
        
        return tv

    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object
        """
        state = NLL.__getstate__ ( self )
        state [ 'ndata' ] = self.__ndata
        return state
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object
        """
        self.__ndata = state.pop ( 'ndata' ) 
        NLL.__setstate__ ( self , state )
    
# =============================================================================

# =============================================================================
## @class ADVAL_LightGBM
#  Use "adversatial validation" method to estimate the Goodness-of-Fit
#  Actually we'll compare the dataset (possible weighted) and MC-dataset generated from PDF
#  @see ADVAL_LGBM 
class ADVAL_LightGBM(GoF) : 
    """ Implementation of concrete method LightGBM-based Adversation Validation for probing of Goodness-Of-Fit
    -   t-value is defined via AUC
    -   p-value if defined via permutations 
    Important parameters:
    
    - mcFactor : (int)   the size of mc-dataset is `mcFactor` times size of real data
    - nToys    : (int)   number of permutations/toys 
    - n_splits : (int)   number of splits for cross-validation 
    
    """
    # =========================================================================
    ## create the estimator
    #  @param mcFactor : (int)  the size of mc-dataset is `mcFactor` times size of real data    
    #  @param nToys    : (int)  number of permutations/toys 
    #  @param n_splits : (int) number of splits for cross-validation 
    def __init__ ( self               ,
                   mcFactor   = 20    , 
                   nToys      = 400   ,
                   parallel   = False ,
                   silent     = False ,
                   progress   = True  ,
                   n_splits   = 8     ,
                   ADVAL_TYPE = None  , **params ) : 
    
        """ Create the Adversarial Validation estimator 

        Parameters  

        - mcFactor : (int) the size of mc-dataset is `mcFactor` times size of real data
        - nToys    : (int) number of permutations/toys 
        - n_splits : (int) number of splits for cross-validation 
        """
        
        if ADVAL_TYPE is None : 
            from ostap.stats.adversarial_validation import ADVAL_LGBM as ADVAL_TYPE
            
        from ostap.stats.adversarial_validation import ADVAL_base
        assert issubclass ( ADVAL_TYPE , ADVAL_base ) , "Invalid ADVAL_TYPE %s" % ADVAL_TYPE 
            
        GoF.__init__ ( self ,
                       gof = ADVAL_TYPE ( nToys    = nToys    ,
                                          parallel = parallel ,
                                          silent   = silent   , 
                                          progress = progress ,
                                          n_splits = n_splits , **params ) ,                       
                       mcFactor = mcFactor )
        
        self.__ecdf   = None
        self.__tvalue = None
        self.__pvalue = None

    
    # =========================================================================
    ## Calculate T-value for Goodness-of-Fit 
    #  @code
    #  ppd    = ...
    #  pdf    = ...  
    #  data   = ... 
    #  tvalue = ppd ( pdf , data ) 
    #  @endcode
    def tvalue ( self , pdf , data ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> ppd    = ...
        >>> pdf    = ... 
        >>> data   = ... 
        >>> tvalue = ppd.tvalue ( pdf , data ) 
        """
        ds1 , ds2 , w1 = self.wtransform ( pdf , data )         
        ## estimate the t-value 
        t = self.gof ( ds1              ,
                       ds2              ,
                       weight1   = w1   ,
                       weight2   = None ,
                       normalize = True )
        ## 
        self.__tvalue = float ( t ) 
        return self.__tvalue 

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
        ds1 , ds2 , w1 = self.wtransform ( pdf , data )         
        ## estimate t&p-values 
        self.__tvalue , self.__pvalue = self.gof.pvalue ( ds1              ,
                                                          ds2              ,
                                                          weight1   = w1   , 
                                                          weight2   = None , 
                                                          normalize = True )        
        return self.__tvalue , self.__pvalue

    # =========================================================================
    @property
    def ecdf ( self ) :
        """`ecdf` : empirical CDF for t-value distribution"""
        return self.gof.ecdf

    # =========================================================================
    ## Get results in a form of the table 
    def table ( self , title = '' , prefix = '' ) :
        """ Get results in a for of the table 
        """
        return self.the_table ( tvalue = self.__tvalue ,
                                pvalue = self.__pvalue ,
                                ecdf   = self.__ecdf   , 
                                title  = title ,
                                prefix = prefx ) 
    __str__  = table
    __repr__ = table

# =============================================================================
## @class ADVAL_HistoGBoost
#  Use "Adversarial Validation" method to estimate the Goodness-of-Fit
#  Actually we'll compare the dataset (possible weighted) and MC-dataset generated from PDF
#  @see ADVAL_HGBC 
class ADVAL_HistoGBoost(ADVAL_LightGBM) : 
    """ Implementation of concrete method XGBoost-based Adversation Validation for probing of Goodness-Of-Fit
    -   t-value is defined via AUC
    -   p-value if defined via permutations
    
    Important parameters:
    
    - mcFactor : (int)   the size of mc-dataset is `mcFactor` times size of real data
    - nToys    : (int)   number of permutations/toys 
    - n_splits : (int)   number of splits for cross-validation 
    
    """
    # =========================================================================
    ## create the estimator
    #  @param mcFactor : (int)  the size of mc-dataset is `mcFactor` times size of real data    
    #  @param nToys    : (int)  number of permutations/toys 
    #  @param n_splits : (int) number of splits for cross-validation 
    def __init__ ( self               ,
                   mcFactor  = 20     , 
                   nToys     = 400    ,
                   parallel  = False  ,
                   silent    = False  ,
                   progress  = True   ,
                   n_splits  = 8      , **params ) : 
    
        """ Create the Adversarial Validation estimator 

        Parameters  

        - mcFactor : (int) the size of mc-dataset is `mcFactor` times size of real data
        - nToys    : (int) number of permutations/toys 
        - n_splits : (int) number of splits for cross-validation 
        """
        from ostap.stats.adversarial_validation import ADVAL_HGBC as ADVAL_TYPE 
        ADVAL_LightGBM.__init__ ( self ,
                                  mcFactor   = mcFactor ,
                                  nToys      = nToys    ,
                                  parallel   = parallel ,
                                  silent     = silent   ,
                                  n_splits   = n_splits ,
                                  ADVAL_TYPE = ADVAL_TYPE , **params )

# =============================================================================
## @class ADVAL_GBoost
#  Use "Adversarial Validation" method to estimate the Goodness-of-Fit
#  Actually we'll compare the dataset (possible weighted) and MC-dataset generated from PDF
#  @see ADVAL_GBC 
class ADVAL_GBoost(ADVAL_LightGBM) : 
    """ Implementation of concrete method XGBoost-based Adversation Validation for probing of Goodness-Of-Fit
    -   t-value is defined via AUC
    -   p-value if defined via permutations
    
    Important parameters:
    
    - mcFactor : (int)   the size of mc-dataset is `mcFactor` times size of real data
    - nToys    : (int)   number of permutations/toys 
    - n_splits : (int)   number of splits for cross-validation 
    
    """
    # =========================================================================
    ## create the estimator
    #  @param mcFactor : (int)  the size of mc-dataset is `mcFactor` times size of real data    
    #  @param nToys    : (int)  number of permutations/toys 
    #  @param n_splits : (int) number of splits for cross-validation 
    def __init__ ( self               ,
                   mcFactor  = 20     , 
                   nToys     = 400    ,
                   parallel  = False  ,
                   silent    = False  ,
                   progress  = True   ,
                   n_splits  = 8      , **params ) : 
    
        """ Create the Adversarial Validation estimator 

        Parameters  

        - mcFactor : (int) the size of mc-dataset is `mcFactor` times size of real data
        - nToys    : (int) number of permutations/toys 
        - n_splits : (int) number of splits for cross-validation 
        """
        from ostap.stats.adversarial_validation import ADVAL_GBC as ADVAL_TYPE 
        ADVAL_LightGBM.__init__ ( self ,
                                  mcFactor   = mcFactor ,
                                  nToys      = nToys    ,
                                  parallel   = parallel ,
                                  silent     = silent   ,
                                  n_splits   = n_splits ,
                                  ADVAL_TYPE = ADVAL_TYPE , **params )
        

# =============================================================================
## @class ADVAL_XGBoost
#  Use "Adversarial Validation" method to estimate the Goodness-of-Fit
#  Actually we'll compare the dataset (possible weighted) and MC-dataset generated from PDF
#  @see ADVAL_XGB 
class ADVAL_XGBoost(ADVAL_LightGBM) : 
    """ Implementation of concrete method XGBoost-based Adversation Validation for probing of Goodness-Of-Fit
    -   t-value is defined via AUC
    -   p-value if defined via permutations
    
    Important parameters:
    
    - mcFactor : (int)   the size of mc-dataset is `mcFactor` times size of real data
    - nToys    : (int)   number of permutations/toys 
    - n_splits : (int)   number of splits for cross-validation 
    
    """
    # =========================================================================
    ## create the estimator
    #  @param mcFactor : (int)  the size of mc-dataset is `mcFactor` times size of real data    
    #  @param nToys    : (int)  number of permutations/toys 
    #  @param n_splits : (int) number of splits for cross-validation 
    def __init__ ( self               ,
                   mcFactor  = 20     , 
                   nToys     = 400    ,
                   parallel  = False  ,
                   silent    = False  ,
                   progress  = True   ,
                   n_splits  = 8      , **params ) : 
    
        """ Create the Adversarial Validation estimator 

        Parameters  

        - mcFactor : (int) the size of mc-dataset is `mcFactor` times size of real data
        - nToys    : (int) number of permutations/toys 
        - n_splits : (int) number of splits for cross-validation 
        """
        from ostap.stats.adversarial_validation import ADVAL_XGB as ADVAL_TYPE 
        ADVAL_LightGBM.__init__ ( self ,
                                  mcFactor   = mcFactor ,
                                  nToys      = nToys    ,
                                  parallel   = parallel ,
                                  silent     = silent   ,
                                  n_splits   = n_splits ,
                                  ADVAL_TYPE = ADVAL_TYPE , **params )
        
# ==============================================================================
if HAS_AVX2 : 
    # ==========================================================================
    ## @class ADVAL_CatBoost
    #  Use "Adversarial Validation" method to estimate the Goodness-of-Fit
    #  Actually we'll compare the dataset (possible weighted) and MC-dataset generated from PDF
    #  @see ADVAL_CATB 
    class ADVAL_CatBoost(ADVAL_LightGBM) : 
        """ Implementation of concrete method CatBoost-based Adversation Validation for probing of Goodness-Of-Fit
        -   t-value is defined via AUC
        -   p-value if defined via permutations
        
        Important parameters:
        
        - mcFactor : (int)   the size of mc-dataset is `mcFactor` times size of real data
        - nToys    : (int)   number of permutations/toys 
        - n_splits : (int)   number of splits for cross-validation 
        
        """
        # =========================================================================
        ## create the estimator
        #  @param mcFactor : (int)  the size of mc-dataset is `mcFactor` times size of real data    
        #  @param nToys    : (int)  number of permutations/toys 
        #  @param n_splits : (int) number of splits for cross-validation 
        def __init__ ( self               ,
                       mcFactor  = 20     , 
                       nToys     = 400    ,
                       parallel  = False  ,
                       silent    = False  ,
                       progress  = True   ,
                       n_splits  = 8      , **params ) : 
            
            """ Create the Adversarial Validation estimator 

            Parameters  
            
            - mcFactor : (int) the size of mc-dataset is `mcFactor` times size of real data
            - nToys    : (int) number of permutations/toys 
            - n_splits : (int) number of splits for cross-validation            
            """
            from ostap.stats.adversarial_validation import ADVAL_CATB as ADVAL_TYPE 
            ADVAL_LightGBM.__init__ ( self ,
                                      mcFactor   = mcFactor ,
                                      nToys      = nToys    ,
                                      parallel   = parallel ,
                                      silent     = silent   ,
                                      n_splits   = n_splits ,
                                      ADVAL_TYPE = ADVAL_TYPE , **params )
            
            __all__ += ( 'ADVAL_CatBoost' , )


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
