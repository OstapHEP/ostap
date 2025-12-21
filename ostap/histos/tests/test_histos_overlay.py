#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developpers.
# =============================================================================
# @file ostap/histos/tests/test_histos_overlay.py
# Test `overlay` functions 
# =============================================================================
""" Test `overlay` functions 
"""
# =============================================================================
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# =============================================================================
from   ostap.core.core        import hID
from   ostap.histos.histos    import histos_overlay
from   ostap.plotting.canvas  import use_canvas
from   ostap.utils.root_utils import batch_env
import ostap.histos.graphs 
import ROOT, random
# =============================================================================
# logging
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ :
    logger = getLogger ( 'ostap.test_histos_overlay' )
else :
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( 'Test for basic operations with histograms')
# =============================================================================
## set batch form environment
batch_env ( logger )
# =============================================================================


def test_overlay () :
    
    h1 = ROOT.TH1D ( hID() , '' ,  100 , -10 , 10 )
    h2 = ROOT.TH1D ( hID() , '' ,  50  , -2  , 10 )
    h3 = ROOT.TH1D ( hID() , '' , 120  , -5  , 5  )
    h4 = ROOT.TH1D ( hID() , '' , 120  , -5  , 15 )

    for i in range ( 10000) :
        x1 = random.gauss   ( -1 , 2   )
        x2 = random.gauss   (  1 , 2   )
        x3 = random.gauss   (  0 , 1   )
        x4 = random.gauss   ( -2 , 2   )
        w2 = random.uniform (  0 , 10  )
        h1.Fill  ( x1 )
        h2.Fill  ( x2 , w2 )
        h3.Fill  ( x3 )
        h4.Fill  ( x4 )
    
    h1.red()
    h2.blue()
    h3.green()
    h4.yellow() 

    for  left_log in ( True , False  ) :
        for right_log in ( True , False ) : 
            for clipx in ( True , False ) : 
                with use_canvas ( "overlay 4 Log left/right=%s/%s clipx=%s" %  ( left_log , right_log , clipx ) , wait = 2  ) : 
                    lefts  = h2 , h1 
                    rights = h3 , h4 
                    histos_overlay ( lefts , rights , right_log = right_log , left_log = left_log , clipx = clipx )
                    
                    
# =============================================================================
if '__main__' == __name__ :

    test_overlay  ()

# =============================================================================
##                                                                      The END
# =============================================================================                