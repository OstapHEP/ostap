#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file ostap/histos/tests/test_histos_roc.py
# Test module for `roc_curve' 
# ============================================================================= 
"""Test module for ostap/histos/roc.py
 - Test module for `roc_curve' 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import 
# ============================================================================= 
from   ostap.core.core        import hID 
from   ostap.histos.roc       import roc_curve
from   ostap.math.integral    import integral 
from   ostap.plotting.canvas  import use_canvas 
from   ostap.utils.root_utils import batch_env 
import ROOT, random, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'ostap.test_histos_roc' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( "Test for `roc_curve'")
# =============================================================================
## set batch form environment 
batch_env ( logger )
# =============================================================================
def test_roc () :

    logger = getLogger ( 'test_roc' )

    logger.info ( 'Test ROC-curve') 

    hs = ROOT.TH1D ( hID() , "signal     distribution" , 15 , -5 ,5 )
    hb = ROOT.TH1D ( hID() , "background distribution" , 15 , -5 ,5 )

    hs.red  ( fill = True , opacity = 0.35 )
    hs.blue ( fill = True , opacity = 0.35 )

    for i in range ( 100 ) :
        hs.Fill ( random.gauss ( +1. , 1. ) ) 
        hb.Fill ( random.gauss ( -1. , 1. ) ) 

    with use_canvas  ( 'EFF-curves/graphs' , wait = 3 ) :
        
        hse = hs.eff ( cut_low = True  )
        hbe = hb.eff ( cut_low = True  )

        hse.red  ()
        hbe.blue ()
                
        grs = hs.eff_graph ( cut_low = True )
        grb = hb.eff_graph ( cut_low = True )
        
        grs.red  ()
        grb.blue ()
        
        hse.draw ()
        hbe.draw ('same')
        grs.draw ('pl')
        grb.draw ('pl')

        logger.info ( 'First bin   efficiencies S:%s B:%s' % \
                      ( hse[ 1].toString( '%+.3f +/- %-.3f' ) ,
                        hbe[ 1].toString( '%+.3f +/- %-.3f' ) ) )
        logger.info ( 'Last  bin   efficiencies S:%s B:%s' % \
                      ( hse[-1].toString( '%+.3f +/- %-.3f' ) ,
                        hbe[-1].toString( '%+.3f +/- %-.3f' ) ) )
        logger.info ( 'First point efficiencies S:%s B:%s' % \
                      ( grs[ 0][1].toString( '%+.3f +/- %-.3f' ) ,
                        grb[ 0][1].toString( '%+.3f +/- %-.3f' ) ) )
        logger.info ( 'Last  point efficiencies S:%s B:%s' % \
                      ( grs[-1][1].toString( '%+.3f +/- %-.3f' ) ,
                        grb[-1][1].toString( '%+.3f +/- %-.3f' ) ) )
                                
    ## make ROC curve
    roc = roc_curve ( signal     = hs   ,
                      background = hb   ,
                      cut_low    = True )
    
    roc.green ()
    with use_canvas ( 'ROC-curve' , wait = 3 ) :
        roc.draw ( 'apc' , minvalue = 0, maxvalue = 1 )

        ## area uner the curve 
        auc = integral ( roc , xmin = 0 , xmax = 1 )
        logger.info( 'Areas under the ROC curve (AUC)= %.3f' % auc )
        

# =============================================================================
if '__main__' == __name__ :

    test_roc ()


# =============================================================================
##                                                                      The END 
# =============================================================================
