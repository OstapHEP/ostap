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
# @see M.Williams, "How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics"
# @see https://doi.org/10.1088/1748-0221/5/09/P09004
# @see http://arxiv.org/abs/arXiv:1003.1768
# @author Vanya Belyaev Ivan.Belyaev@cern.ch
# @date 2011-09-21
#
# ============================================================================
""" `U-statistics' useful for `Goodness-Of-Fit' tests

This is a simple translation of the original C++ lines written by Greig Cowan into Python

- see M.Williams, ``How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics''
- see https://doi.org/10.1088/1748-0221/5/09/P09004
- see http://arxiv.org/abs/arXiv:1003.1768

Usage is fairly trivial:

   >>> pdf  = ...               ## pdf
   >>> data = ...               ## dataset
   >>> pdf.fitTo( data , ... )  ## fit it!

   >>> import ostap.stats.ustat as uStat

   >>> r,histo = uStat.uPlot ( pdf , data ) 
   >>> print r                  ## print fit results
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
    )
# ============================================================================
from   ostap.core.core import Ostap, hID
import ostap.histos.histos
import ROOT, math, ctypes
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.ustat' )
else                       : logger = getLogger ( __name__      )
# =============================================================================

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
    """Calculate U-statistics 
    """
    
    if not isinstance ( pdf , ROOT.RooAbsPdf ) or not pdf :
        from   ostap.fitting.pdfbasic import APDF1 
        assert pdf and isinstance ( pdf , APDF1 ) , "Invalid type of `pdf'!"
        pdf = pdf.pdf
        
    if not args  : args  = pdf.getObservables ( data )  
    if not histo : histo = ROOT.nullptr

    ##
    tStat = ctypes.c_double (-1)
    sc    = Ostap.UStat.calculate ( pdf   ,
                                    data  ,
                                    tStat ,
                                    histo ,
                                    args  )
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
#   @param bins   (input) bumbef of bins in histogram 
#   @param silent (input) keep the silence 
def uPlot ( pdf            ,
            data           ,
            bins   = None  ,
            args   = None  ,
            silent = False ) :
    """Make the plot of U-statistics 
    
    >>> pdf  = ...               ## pdf
    >>> data = ...               ## dataset
    >>> pdf.fitTo( data , ... )  ## fit it!
    
    >>> import ostap.stats.ustat as uStat
    
    >>> t, res , histo = uStat.uPlot ( pdf , data ) 
    >>> print ( res )            ## print fit results
    >>> histo.Draw()             ## plot the results  
    """
    
    if not bins or bins <= 0 :
        nEntries = float(data.numEntries())
        bins     = 10 
        for nbins in ( 1000 , 500 ,
                       200  , 100 ,
                       50   , 40  ,
                       25   , 20  ,
                       16   , 10  ,
                       8    , 5   ) :
            if nEntries/float(nbins) < 100 : continue  
            bins = nbins
            break

    histo = ROOT.TH1F ( hID () ,'U-statistics', bins , 0 , 1 )
    histo.Sumw2      (   )
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

if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# ===========================================================================
##                                                                    The END 
# ===========================================================================
