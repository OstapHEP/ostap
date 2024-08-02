#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Simple utility to creat ROC-like curves/gpaphs
#  - signal effciency vs background efficiency
#  - signal effciency vs background retention 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-07 
# =============================================================================
""" Simple utility to creat ROC-like curves/gpaphs
  - signal effciency vs background efficiency
  - signal effciency vs background retention 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-08-02"
__all__     = (
    'makeGraph'   , # make graph from primitive data
    ) 
# =============================================================================
from   ostap.core.ostap_types import string_types 
import ostap.histos.histos
import ostap.histos.graphs 
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.histos.roc' )
else                       : logger = getLogger( __name__           )
# =============================================================================
## Build the ROC-curve from signal and background disctributuions
#  @param signal    (histogram) of signal     distribution
#  @param backgrund (histogram) of background distribution
#  @return ROC-curve
#  @code
#  hsignal = ... ## signal distribution 
#  hbkg    = ... ## background distribution
#  roc     = roc_curve ( signal         = hsig ,
#                        backgrund      = hbkg ,
#                        increasing     = True , ## "keep valeus less than cut value"
#                        show_sinal     = 'efficiency' ,
#                        show_backgrund = 'rejection'  )
#  ## get AUC
#  import ostap.math,integral as I
#  auc  = I.integral ( roc , xmin = 0 , xmax = 1.0 ) 
#  @endcode 
def roc_curve ( signal                         , 
                background                     ,
                increasing                     ,                
                show_signal     = 'efficiency' , 
                show_background = 'rejection'  ) :

    """Build the ROC-curve from signal and background disctributuions
    - signal    : (histogram) of signal     distribution
    - backgrund : (histogram) of background distribution

    >>> hsignal = ... ## signal distribution 
    >>> hbkg    = ... ## background distribution
    >>> roc     = roc_curve ( signal         = hsig ,
    ...                       backgrund      = hbkg ,
    ...                       increasing     = True , ## "keep vaues that are less than the cut value"
    ...                       show_sinal     = 'efficiency' ,
    ...                       show_backgrund = 'rejection'  )
      
    >>> import ostap.math,integral as I
    >>> auc  = I.integral ( roc , xmin = 0 , xmax = 1.0 ) 
    
    """
    assert isinstance ( signal     , ROOT.TH1 ) and 1 == signal    .dim() , \
        "Invalid `signal' type: %s" % type ( signal )
    assert isinstance ( background , ROOT.TH1 ) and 1 == background.dim() , \
        "Invalid `background' type: %s" % type ( background )

    sig_fun = show_signal
    if callable ( sig_fun ) : pass
    else :
        assert isinstance ( sig_fun , string_types ) , "Invalid type of `show_signal' %s" % type ( sig_fun )
        sig_fun = str(sig_fun).strip().lower()
        if   sig_fun in ( 'e' , 'eff' , 'effic'   , 'efficiency'               ) : sig_fun = lambda s : s
        elif sig_fun in ( 'r' , 'rej' , 'rejec'   , 'reject'   , 'rejection'   ) : sig_fun = lambda s : 1-s
        else :
            raise TypeError ("Unknown `show_signal' :%s" % show_signal ) 

    bkg_fun = show_background
    if callable ( bkg_fun ) : pass
    else :
        assert isinstance ( bkg_fun , string_types ) , "Invalid type of `show_background' %s" % type ( bkg_fun )
        bkg_fun = str(bkg_fun).strip().lower()
        if   bkg_fun in ( 'e' , 'eff' , 'effic'   , 'efficiency'               ) : bkg_fun = lambda s : s
        elif bkg_fun in ( 'r' , 'rej' , 'rejec'   , 'reject'   , 'rejection'   ) : bkg_fun = lambda s : 1-s
        else :
            raise TypeError ("Unknown `show_background' :%s" % show_backgrund ) 

                                                                       
    hs = signal
    hb = background 
    
    a## signal efficiency
    hse = hs.eff ( increasing = increasing )
    
    ## background efficiency
    hbe = hb.eff ( increasing = increasing )

    ## output graph: ROC curve
    np    = len ( hse ) 
    graph = ROOT.TGraphErrors ( np )

    ## loop over signal efficiency 
    for i , xs , es in hse.items() :

        ## backrgiund efficiency
        eb = hbe ( xs.value() )

        ## tarnsform if requested:
        
        es = sig_fun ( es ) 
        eb = bkg_fun ( eb )
        
        graph [ i - 1 ] = es , eb

    return graph 
    

## ============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
