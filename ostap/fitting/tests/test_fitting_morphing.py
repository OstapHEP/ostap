#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_transform.py
# Test module 
# - It tests transform-PDF  
# ============================================================================= 
""" Test module 
- It tests transform-PDF 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ROOT, random
import ostap.fitting.roofit
import ostap.fitting.models        as     Models
from   ostap.fitting.morphing_pdf  import Morphing1D_pdf 
from   ostap.utils.utils           import vrange
from   builtins                    import range
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_morphing' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

mass = ROOT.RooRealVar ( 'mass' , 'some mass' , 0 , 20 ) 
h1   = ROOT.TH1D ('h' ,'' , 200 , 0 , 20 )
N    = 10000
for i in range ( N ) :
    h1.Fill ( random.gauss ( 10 , 2.5 ) ) 

def test_morphing1 () :

    pdf1 = Models.Gauss_pdf ( 'G1' , xvar = mass , mean = 10 , sigma = 1 )
    pdf2 = Models.Gauss_pdf ( 'G2' , xvar = mass , mean = 10 , sigma = 2 )
    pdf3 = Models.Gauss_pdf ( 'G3' , xvar = mass , mean = 10 , sigma = 3 )
    pdf4 = Models.Gauss_pdf ( 'G4' , xvar = mass , mean = 10 , sigma = 4 )

    pdf  = Morphing1D_pdf ( 'M' , { 1.0 : pdf1 ,
                                    2.0 : pdf2 ,
                                    3.0 : pdf3 ,
                                    4.0 : pdf4 , } , xvar =  mass )
    
    for mu in vrange ( 1 , 3 , 6 ) :
        pdf.mu = mu
        logger.info ( 'Mu= %s' % mu ) 
        pdf.draw()

    r , f = pdf.fitHisto ( h1 , draw = True , nbins = 100 , silent = True )
    logger.info ( 'Morphing: \n%s' % r.table ( prefix = "# " ) ) 
    
# =============================================================================
if '__main__' == __name__ :

    test_morphing1   () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
