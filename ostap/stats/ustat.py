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
#   >>> pdf.fitTo( data , ... )  ## fit it!
#
#   >>> import ostap.stats.ustat as uStat
#
#   >>> r,histo = uStat.uPlot ( pdf , data ) 
#   >>> print r                  ## print fit results
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
""" ``U-statistics'' useful for ``Goodness-Of-Fit'' tests

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
from   __future__        import print_function
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2010-09-21"
__version__ = "$Revision$"
# ============================================================================
__all__     = (
    "uPlot" ,  ## make  plot of U-statistics 
    "uDist" ,  ## calculate  U-statistics 
    "uCalc" ,  ## calclulate the distance between two data points 
    )
# ============================================================================
from   ostap.core.core import cpp, Ostap, hID
import ostap.histos.histos
import ROOT, math, ctypes
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.ustat' )
else                       : logger = getLogger ( __name__      )
# =============================================================================
## calculate the distance between two data points 
#  @author Vanya Belyaev Ivan.Belyaev@cern.ch
#  @date 2011-09-21
def uDist ( x , y ) :
    """Calculate the distance between two data points 
    """
    
    ix = Ostap.Utils.Iterator ( x )
    iy = Ostap.Utils.Iterator ( y )
    
    dist = 0.0
    
    xv = ix.Next()
    yv = iy.Next()
    while xv and yv :

        if not hasattr ( xv  , 'getVal' ) : break 
        if not hasattr ( yv  , 'getVal' ) : break
        
        d     = xv.getVal()-yv.getVal()
        dist += d*d
        xv = ix.Next()
        yv = iy.Next()

    del ix
    del iy

    return math.sqrt( dist )

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
            args           , 
            data           ,
            histo          ,
            silent = False )  :
    """Calculate U-statistics 
    """
    import sys
    
    tStat = ctypes.c_double (-1)
    sc    = Ostap.UStat.calculate ( pdf   ,
                                    data  ,
                                    histo ,
                                    tStat ,
                                    args  )
    tStat = float ( tStat.value  ) 
    return histo, tStat 
    
    numEntries = data.numEntries ()
    dim        = args.getSize    ()
    data_clone = data.Clone      ()

    from ostap.logger.progress_bar import progress_bar 
    for i in progress_bar ( xrange ( numEntries ) ) : 

        event_x = data_clone.get ( i ) 
        event_i = event_x.selectCommon ( args )
        
        ## fill args and evaluate PDF 
        for a in args : a.setVal ( event_i.getRealValue ( a.GetName () ) )        
        pdfValue = pdf.getVal ( args )
        
        small_v = 1.e+100 
        small_j = 0
        for j in xrange ( 0 , numEntries ) :
            
            if j == i  : continue

            event_y = data.get( j )
            event_j = event_y.selectCommon ( args ) 

            dist = uDist ( event_i , event_j )
            
            if 0 == j or dist < small_v :
                small_v = dist
                small_j = j 
                
        value = 0 
        if   1 == dim :
            value  = small_v
            value *= numEntries * pdfValue 
        elif 2 == dim :
            value  = small_v**2
            value *= numEntries * pdfValue
            value *= math.pi 
        else :
            logger.error ( ' Not-implemented (yet) %s '  % dim )
            continue 

        value = math.exp ( -1 * value )

        histo.Fill ( value ) 

    del data_clone 
    return histo

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
#    >>> r,histo = uStat.uPlot ( pdf , data ) 
#    >>> print r                  ## print fit results
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
    
    >>> r,histo = uStat.uPlot ( pdf , data ) 
    >>> print r                  ## print fit results
    >>> histo.Draw()             ## plot the results  
    """

    if not bins or bins <= 0 :
        nEntries = float(data.numEntries())
        bins = 10 
        for nbins in ( 50  ,
                       40  ,
                       25  ,
                       20  ,
                       16  ,
                       10  ,
                       8   ,
                       5   ) :
            if nEntries/nbins < 100 : continue  
            bins = nbins
            break
        print('#bins %s' % bins)
        
    histo = ROOT.TH1F ( hID () ,'U-statistics', bins , 0 , 1 )
    histo.Sumw2      (   )
    histo.SetMinimum ( 0 )

    if not args : args = pdf.getObservables ( data )
    
    h,tStat = uCalc ( pdf       ,
                      args      ,
                      data      ,
                      histo     ,
                      silent    )    
    
    res  = histo.Fit         ( 'pol0' , 'SLQ0+' )
    func = histo.GetFunction ( 'pol0' )
    if func :
        func.SetLineWidth ( 3 )
        func.SetLineColor ( 2 )
        func.ResetBit     ( 1 << 9 )
        
    return res , histo, float(tStat)

# ===========================================================================

if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# ===========================================================================
# The END 
# ===========================================================================
