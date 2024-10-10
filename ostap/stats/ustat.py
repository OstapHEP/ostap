#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================
# @file ostap/stats/ustat.py
#
# Helper module to get ``U-statistics'' useful for ``Goodnes-Of-Fit'' tests
#
# This is a simple translation of the original C++ lines written by Greig Cowan into python
#
# Usage is fairly trivial:
#
#  @code
# 
#   >>> pdf  = ...               ## pdf
#   >>> data = ...               ## dataset
#   >>> pdf.fitTo ( data , ... ) ## fit it!
#
#   >>> import ostap.stats.ustat as uStat
#
#   >>> r,histo = uStat.uPlot ( pdf , data ) 
#   >>> print ( r )              ## print fit results
#   >>> histo.Draw()             ## plot the results  
#
#  @endcode
#
# It is a numpy-free version of `Distance-t0-Nearedt-Neighbour" method
# described in M.Williams' paper
# @see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
# @see https://doi.org/10.1088/1748-0221/5/09/P09004
# @see http://arxiv.org/abs/arXiv:1003.1768
# @author Vanya Belyaev Ivan.Belyaev@cern.ch
# @date 2011-09-21
#
# ============================================================================
""" `U-statistics' useful for `Goodness-Of-Fit' tests

This is a simple translation of the original C++ lines written by Greig Cowan into Python
It is a numpy-free version of `Distance-t0-Nearedt-Neighbour" method described in M/williamns' paper
- see M.Williams, `How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics'
- see https://doi.org/10.1088/1748-0221/5/09/P09004
- see http://arxiv.org/abs/arXiv:1003.1768

Usage is fairly trivial:

   >>> pdf  = ...               ## pdf
   >>> data = ...               ## dataset
   >>> pdf.fitTo( data , ... )  ## fit it!

   >>> import ostap.stats.ustat as uStat

   >>> r,histo = uStat.uPlot ( pdf , data ) 
   >>> print ( r  )             ## print fit results
   >>> histo.Draw()             ## plot the results     
"""
# ============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2010-09-21"
__version__ = "$Revision:$"
# ============================================================================
__all__     = (
    "uPlot" ,  ## make  plot of U-statistics 
    "uCalc" ,  ## calculate U-statistics
    "uToys" ,  ## calculate p-value using toys
    "USTAT" ,  ## the same but using  a common GoF interface 
    )
# ============================================================================
from   ostap.core.ostap_types   import integer_types 
from   ostap.core.core          import Ostap, hID
from   ostap.stats.counters     import EffCounter
from   ostap.stats.gof          import AGoF
from   ostap.stats.gof_utils    import TOYS  
from   ostap.utils.progress_bar import progress_bar
import ostap.histos.histos
import ROOT, ctypes
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.ustat' )
else                       : logger = getLogger ( __name__      )
# =============================================================================
##  calculate U-statistics
#   @param pdf    (input) PDF
#   @param args   (input) arguments/variables
#   @param data   (input) dataset 
#   @param histo  (input) the histogram to be filled 
#   @author Vanya Belyaev Ivan.Belyaev@cern.ch
#   @see Analysis::UStat
#   @see Analysis::UStat::calculate
#   @date 2011-09-21
def uCalc ( pdf            ,
            data           ,
            args   = None  , 
            histo  = None  ,
            silent = False )  :
    """ Calculate U-statistics 
    """
    
    if not isinstance ( pdf , ROOT.RooAbsPdf ) or not pdf :
        from   ostap.fitting.pdfbasic import APDF1 
        assert pdf and isinstance ( pdf , APDF1 ) , "Invalid type of `pdf'!"
        pdf = pdf.pdf
        
    if not args  : args  = pdf.getObservables ( data )  
    if not histo : histo = ROOT.nullptr

    ##
    tStat = ctypes.c_double ( -1 )
    from ostap.utils.progress_conf import progress_conf
    if silent : 
        sc = Ostap.UStat.calculate ( pdf   ,
                                     data  ,
                                     tStat ,
                                     histo ,
                                     args  )
    else :
        from ostap.utils.progress_conf import progress_conf
        sc = Ostap.UStat.calculate ( progress_conf ( description = 'Entries:') ,
                                     pdf              ,
                                     data             ,
                                     tStat            ,
                                     histo            ,
                                     args             )
        
        
    if sc.isFailure() :
        logger.error ( "Error from Ostap::UStat::Calculate %s" % sc )

    if not histo : histo = None 
    
    tStat = float ( tStat.value  ) 
    return tStat, histo  
    
# =============================================================================
##  make the plot of U-statistics
#
#   @code
#
#    >>> pdf  = ...               ## pdf
#    >>> data = ...               ## dataset
#    >>> pdf.fitTo( data , ... )  ## fit it!
#    >>> vars = ...               ## get variables
#    
#    >>> import ostap.stats.ustat as uStat
#    
#    >>> t , res , histo = uStat.uPlot ( pdf , data ) 
#    >>> print ( res )            ## print fit results
#    >>> histo.Draw()             ## plot the results
#
#   @endcode
#
#   @param pdf    (input) PDF
#   @param args   (input) arguments/variables 
#   @param data   (input) dataset 
#   @param bins   (input) number  of bins in histogram 
#   @param silent (input) keep the silence 
def uPlot ( pdf            ,
            data           ,
            histo  = None  ,
            args   = None  ,
            silent = False ) :
    """ Make the plot of U-statistics 
    
    >>> pdf  = ...               ## pdf
    >>> data = ...               ## dataset
    >>> pdf.fitTo( data , ... )  ## fit it!
    
    >>> import ostap.stats.ustat as uStat
    
    >>> t, res , histo = uStat.uPlot ( pdf , data ) 
    >>> print ( res )            ## print fit results
    >>> histo.Draw()             ## plot the results  
    """

    if   isinstance ( histo , ROOT.TH1 )                    : pass
    elif isinstance ( histo , integer_types ) and 1 < histo :
        histo = ROOT.TH1F ( hID () ,'U-statistics', histo , 0 , 1 )
        histo.Sumw2()         
    elif not histo :   
        nEntries = float ( len ( data ) )
        bins     = 10 
        for nbins in ( 1000 , 500 ,
                       200  , 100 ,
                       50   , 40  ,
                       25   , 20  ,
                       16   , 10  ,
                       8    , 5   ) :
            if nEntries / float ( nbins ) < 100 : continue  
            bins = nbins
            break
        histo = ROOT.TH1F ( hID () ,'U-statistics', bins , 0 , 1 )
        histo.Sumw2()         
    else :
        raise TypeError ( "Invalid type of ``histo'':%s" % type ( histo ).__name__  )
    
    histo.SetMinimum ( 0 )
    tStat , hh = uCalc ( pdf       ,
                         data      ,
                         args      ,
                         histo     ,
                         silent    )    

    res  = histo.Fit         ( 'pol0' , 'SLQ0+' )
    func = histo.GetFunction ( 'pol0' )
    if func :
        func.SetLineWidth ( 3 )
        func.SetLineColor ( 2 )
        func.ResetBit     ( 1 << 9 )
        
    return float ( tStat ) , histo , res  


# ===========================================================================
## get p-value for GoF using toys
#  @code
#  pdf  = ...
#  data = ...
#  t_value , p_value , histo = uToys ( pdf , data , nToys = 1000 )
#  @endcode 
def uToys ( pdf            ,
            data           ,
            nToys  = 1000  ,
            histo  = None  , 
            args   = None  ,
            silent = False ,
            sample = True  ) :
    """ Get p-value for GoF using toys
    >>> pdf  = ...
    >>> data = ...
    >>> t_value , p_value , histo = uToys ( pdf , data , nToys = 1000 )
    """
    t_value, histo, _ = uPlot ( pdf , data , histo = histo , args = args , silent = silent )

    N       = len ( data )
    counter = EffCounter()
    
    for i in progress_bar ( nToys , silent = silent , description = 'Toys: ' ) :
        
        ds = pdf.generate ( N , sample = sample )
        ti, _ = uCalc ( pdf , ds , silent = True )
        counter += bool ( t_value < ti  )

        if isinstance ( ds , ROOT.RooDataSet ) :
            ds = Ostap.MoreRooFit.delete_data ( ds )
        del ds
        
    p_value = 1 - counter.eff
    return t_value, p_value, histo 

# ===========================================================================
## @class USTAT
# Goodness-of-Fit estimator for Distance-to-Nearest-Neighbour GoF test
# It is a numpy-free version of `Distance-t0-Nearedt-Neighbour" method
# described in M.Williams' paper
# @see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
# @see https://doi.org/10.1088/1748-0221/5/09/P09004
# @see http://arxiv.org/abs/arXiv:1003.1768
# @see Ostap::
class USTAT(AGoF) :
    """ Goodness-of-Fit estimator for Distance-to-Nearest-Neighbour GoF test
    It is a numpy-free version of `Distance-t0-Nearedt-Neighbour" method described in M.Williams' paper
    - see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
    - see https://doi.org/10.1088/1748-0221/5/09/P09004
    - see http://arxiv.org/abs/arXiv:1003.1768
    - see `Ostap.UStat`
    """
    def __init__ ( self             ,
                   histo    = None  ,
                   nToys    = 1000  ,
                   sample   = False ,
                   parallel = True  , 
                   silent   = False ) :
        
        self.__silent   = True if silent else False
        self.__sample   = True if sample else False 
        self.__histo    = None 
        self.__bins     = histo
        self.__nToys    = nToys 
        self.__parallel = True if parallel else False 
    # =========================================================================
    ## Calculate T-value for Goodness-of-Git 
    #  @code
    #  ustat   = ...
    #  pdf     = ...  
    #  data    = ... 
    #  tvalue  = ustat ( pdf , data ) 
    #  @endcode
    def __call__ ( self , pdf , data ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> ustat  = ...
        >>> pdf    = ... 
        >>> data   = ... 
        >>> tvalue = ustat ( pdf , data ) 
        """

        ## get t-value 
        t_value , histo , _ = uPlot ( pdf  ,
                                      data , 
                                      histo  = self.__bins , 
                                      silent = self.silent ) 

        if histo and not self.__histo :
            self.__histo = histo.clone()
            
        return t_value

    # =========================================================================
    ## Calculate the t & p-values
    #  @code
    #  ustat = ...
    #  pdf   = ...
    #  data  = ... 
    #  t , p = ustat.pvalue ( pdf , data )
    #  @endcode 
    def pvalue ( self , pdf , data ) :
        """ Calculate the t & p-values
        >>> ustat  = ...
        >>> pdf   = ... 
        >>> data  = ... 
        >>> t , p = ustat.pvalue ( pdf , data ) 
        """

        estimator = self 
        t_value   = estimator ( pdf , data )        

        ## prepare toys
        toys = TOYS ( self , t_value , pdf = pdf , Ndata = len ( data ) , sample = self.sample )
        
        silent = self.silent
        self.__silent = True 
        if self.parallel :
            counter = toys.run ( self.nToys , silent = silent )
        else :
            counter = toys     ( self.nToys , silent = self.silent )            

        self.__silent = silent 
                
        p_value = 1 - counter.eff
        return t_value, p_value 

    @property
    def histo ( self )  :
        """`histo` : get histogram with U-values distribution"""
        return self.__histo

    @property
    def nToys ( self )  :
        """`nToys` : number of toys fro p-value evaluation"""
        return self.__nToys
    
    @property
    def silent ( self ) :
        """`silent` : silent processing?"""
        return self.__silent
    
    @property
    def sample ( self ) :
        """`sample` : sample number of events for toys?"""
        return self.__sample

    @property
    def parallel ( self ) :
        """`parallel` : parallel processing where/when/if possible?"""
        return self.__parallel
    

# ===========================================================================

if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# ===========================================================================
##                                                                    The END 
# ===========================================================================
