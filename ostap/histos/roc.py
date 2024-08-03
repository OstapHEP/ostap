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
    'roc_curve'   , # Make ROC curve form signal & backgrund distributions 
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
## symbols to indicate the efficiency
_effs = ( 'e' , 'eff' , 'effs'   , 'effic'     , 'efficiency' )
## symbols to indicate the rejection 
_rejs = ( 'r' , 'rej' , 'reject' , 'rejection' )
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
#                        cut_low        = True , ## "keep valeus less than cut value"
#                        show_sinal     = 'efficiency' ,
#                        show_backgrund = 'rejection'  )
# 
#  import ostap.math,integral as I
#  auc  = I.integral ( roc , xmin = 0 , xmax = 1.0 ) 
#  @endcode 
def roc_curve ( signal                         , 
                background                     ,
                cut_low                        ,                
                show_signal     = 'efficiency' , 
                show_background = 'rejection'  ) :

    """Build the ROC-curve from signal and background disctributuions
    - signal    : (histogram) of signal     distribution
    - backgrund : (histogram) of background distribution
    
    >>> hsignal = ... ## signal distribution 
    >>> hbkg    = ... ## background distribution
    >>> roc     = roc_curve ( signal          = hsig ,
    ...                       backgrund       = hbkg ,
    ...                       cut_low         = True , ## "keep vaues that are less than the cut value"
    ...                       show_sinal      = 'efficiency' ,
    ...                       show_background = 'rejection'  )
      
    >>> import ostap.math,integral as I
    >>> auc  = I.integral ( roc , xmin = 0 , xmax = 1.0 ) 
    
    """
    assert isinstance ( signal     , ROOT.TH1 ) and 1 == signal    .dim() , \
        "Invalid `signal' type: %s" % type ( signal )
    assert isinstance ( background , ROOT.TH1 ) and 1 == background.dim() , \
        "Invalid `background' type: %s" % type ( background )
    
    def _fun_ ( obj ) :
        if callable ( obj ) : return obj
        assert isinstance ( obj , string_types ) , 'Invalid type: %s' % type ( obj )
        obj  = str ( obj ).strip ().lower () 
        if   obj in _effs : return lambda e : e
        elif obj in _rejs : return lambda e : 1.0-e
        raise TypeError ( 'Invalid object: %s' % obj )

    ## transformations :
    
    sig_fun = _fun_ ( show_signal     )
    bkg_fun = _fun_ ( show_background )

                                                                       
    hs = signal
    hb = background 
    
    ## signal efficiency histogram 
    hse = hs.eff ( cut_low = cut_low  )
    
    ## background efficiency histogram 
    hbe = hb.eff ( cut_low = cut_low )
    
    ## output graph: ROC curve
    np    = len ( hse ) 
    graph = ROOT.TGraphErrors ( np )

    ## loop over signal efficiency 
    for i , xs , es in hse.items() :

        ## background efficiency
        eb = hbe ( xs.value() )

        ## transform if requested:
        
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
