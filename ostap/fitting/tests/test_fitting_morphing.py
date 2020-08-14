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
from   ostap.fitting.morphing_pdf  import Morphing1D_pdf, Morphing2D_pdf 
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

# ============================================================================
def test_morphing1 () :

    pdf1 = Models.Gauss_pdf ( 'G1' , xvar = mass , mean = 10 , sigma = 1 )
    pdf2 = Models.Gauss_pdf ( 'G2' , xvar = mass , mean = 10 , sigma = 2 )
    pdf3 = Models.Gauss_pdf ( 'G3' , xvar = mass , mean = 10 , sigma = 3 )
    pdf4 = Models.Gauss_pdf ( 'G4' , xvar = mass , mean = 10 , sigma = 4 )

    pdf  = Morphing1D_pdf ( 'M1' , { 1.0 : pdf1 ,
                                     2.0 : pdf2 ,
                                     3.0 : pdf3 ,
                                     4.0 : pdf4 , } , xvar =  mass )
    
    for mu in vrange ( 1 , 3 , 6 ) :
        pdf.mu = mu
        logger.info ( 'Mu= %s' % mu ) 
        pdf.draw()

    r , f = pdf.fitHisto ( h1 , draw = True , nbins = 100 , silent = True )
    logger.info ( 'Morphing: \n%s' % r.table ( prefix = "# " ) ) 

# ============================================================================
## def test_morphing2 () :
if 1 < 2 :
    
    pdf11 = Models.Gauss_pdf ( 'G11' , xvar = mass , mean =  8 , sigma = 1 )
    pdf12 = Models.Gauss_pdf ( 'G12' , xvar = mass , mean = 10 , sigma = 1 )
    pdf13 = Models.Gauss_pdf ( 'G13' , xvar = mass , mean = 12 , sigma = 1 )
    pdf21 = Models.Gauss_pdf ( 'G21' , xvar = mass , mean =  8 , sigma = 2 )
    pdf22 = Models.Gauss_pdf ( 'G22' , xvar = mass , mean = 10 , sigma = 2 )
    pdf23 = Models.Gauss_pdf ( 'G23' , xvar = mass , mean = 12 , sigma = 2 )
    pdf31 = Models.Gauss_pdf ( 'G31' , xvar = mass , mean =  8 , sigma = 3 )
    pdf32 = Models.Gauss_pdf ( 'G32' , xvar = mass , mean = 10 , sigma = 3 )
    pdf33 = Models.Gauss_pdf ( 'G33' , xvar = mass , mean = 12 , sigma = 3 )    
    pdf41 = Models.Gauss_pdf ( 'G41' , xvar = mass , mean =  8 , sigma = 4 )
    pdf42 = Models.Gauss_pdf ( 'G42' , xvar = mass , mean = 10 , sigma = 4 )
    pdf43 = Models.Gauss_pdf ( 'G43' , xvar = mass , mean = 12 , sigma = 4 )    
    
    pdf  = Morphing2D_pdf ( 'M2' , { ( 8,1) : pdf11 ,
                                     (10,1) : pdf12 ,
                                     (12,1) : pdf13 ,
                                     ( 8,2) : pdf21 ,
                                     (10,2) : pdf22 ,
                                     (12,2) : pdf23 ,                                     
                                     ( 8,3) : pdf31 ,
                                     (10,3) : pdf32 ,
                                     (12,3) : pdf33 ,
                                     ( 8,4) : pdf41 ,
                                     (10,4) : pdf42 ,
                                     (12,4) : pdf43 } , 
                            xvar =  mass ,
                            setting = ROOT.RooMomentMorphND.Linear ) 
    
    r , f = pdf.fitHisto ( h1 , draw = False , silent = True )
    r , f = pdf.fitHisto ( h1 , draw = True  , nbins = 100 , silent = True )
    logger.info ( 'Morphing: \n%s' % r.table ( prefix = "# " ) ) 
    
# =============================================================================
if '__main__' == __name__ :

    pass 
    ## test_morphing1   () 
    ## test_morphing2   () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
