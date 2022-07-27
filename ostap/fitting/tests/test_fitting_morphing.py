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
import ostap.fitting.roofit
import ostap.fitting.models        as     Models
from   ostap.core.meta_info        import root_info 
from   ostap.fitting.morphing_pdf  import MorphingN1_pdf, MorphingN2_pdf 
from   ostap.utils.utils           import vrange
from   builtins                    import range
from   ostap.utils.timing          import timing 
from   ostap.plotting.canvas       import use_canvas
from   ostap.utils.utils           import wait 
import ROOT, random
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

    logger = getLogger ('test_morphing1')    
    if root_info < ( 6 , 23 )  or ( 6 , 27 ) <= root_info :
        logger.warning( 'Test is disabled for ROOT version %s' % ROOT.gROOT.GetVersion() )
        return

    shapes = {}

    mean        = 10 
    sigma_range = 0.5 , 5.0
    
    for i , sigma in  enumerate ( vrange ( *sigma_range , 10 ) ) :

        gauss = Models.Gauss_pdf ( 'G1_%d' % i ,
                                   xvar  = mass ,
                                   mean  = ROOT.RooFit.RooConst ( mean  ) ,
                                   sigma = ROOT.RooFit.RooConst ( sigma ) )
        shapes [ sigma ] = gauss
        
    ## create morphing PDF 
    pdf  = MorphingN1_pdf ( 'M1' , shapes , xvar =  mass )
        
    for mu in vrange ( *sigma_range , 6 ) :
        pdf.mu = mu
        logger.info ( 'Mu= %s' % mu ) 
        with wait ( 0.2 ) , use_canvas ( 'test_morphing1' ) :
            pdf.draw()

    r , f = pdf.fitHisto ( h1 , draw = False , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_morphing1' ) :
        r , f = pdf.fitHisto ( h1 , draw = True , nbins = 100 , silent = True )
        logger.info ( 'Morphing: \n%s' % r.table ( prefix = "# " ) ) 

# ============================================================================
def test_morphing2 () :
    
    logger = getLogger ('test_morphing3')
    
    if root_info < ( 6 , 23 )  or ( 6 , 27 ) <= root_info :
        logger.warning( 'Test is disabled for ROOT version %s' % ROOT.gROOT.GetVersion() )
        return


    shapes = {}

    sigma       = 2.5
    mean_range  = 8.0 , 12.0
    
    for j , mean in  enumerate ( vrange ( *mean_range , 10 ) ) :
        gauss = Models.Gauss_pdf ( 'G2_%d' % j  , 
                                   xvar  = mass ,
                                   mean  = ROOT.RooFit.RooConst ( mean  ) ,
                                   sigma = ROOT.RooFit.RooConst ( sigma ) )
        shapes [ mean ] = gauss
        
    ## create morphing PDF 
    pdf  = MorphingN1_pdf ( 'M2'    , shapes , 
                            xvar    =  mass ,
                            setting = ROOT.RooMomentMorphND.Linear )
    
    r , f = pdf.fitHisto ( h1 , draw = False , silent = True )
    r , f = pdf.fitHisto ( h1 , draw = False , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_morphing2' ) :
        r , f = pdf.fitHisto ( h1 , draw = True  , nbins = 100 , silent = True )
        logger.info ( 'Morphing: \n%s' % r.table ( prefix = "# " ) ) 
    
# ============================================================================
def test_morphing3 () :
    
    logger = getLogger ('test_morphing3')
    
    if root_info < ( 6 , 23 )  or ( 6 , 27 ) <= root_info :
        logger.warning( 'Test is disabled for ROOT version %s' % ROOT.gROOT.GetVersion() )
        return


    shapes = {}
    
    sigma_range = 0.5 ,  5.0
    mean_range  = 8.0 , 12.0
    
    for i , sigma in  enumerate ( vrange ( *sigma_range , 10 ) ) :
        for j , mean in  enumerate ( vrange ( *mean_range , 10 ) ) :
            gauss = Models.Gauss_pdf ( 'G3_%d_%d' % ( i , j ) , 
                                       xvar  = mass ,
                                       mean  = ROOT.RooFit.RooConst ( mean  ) ,
                                       sigma = ROOT.RooFit.RooConst ( sigma ) )
            shapes [ mean , sigma ] = gauss

    ## create morphing PDF 
    pdf  = MorphingN2_pdf ( 'M3'    , shapes , 
                            xvar    =  mass ,
                            setting = ROOT.RooMomentMorphND.Linear )
            
    r , f = pdf.fitHisto ( h1 , draw = False , silent = True )
    r , f = pdf.fitHisto ( h1 , draw = False , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_morphing3' ) :
        r , f = pdf.fitHisto ( h1 , draw = True  , nbins = 100 , silent = True )
        logger.info ( 'Morphing: \n%s' % r.table ( prefix = "# " ) ) 
    

# =============================================================================
if '__main__' == __name__ :

    with timing ("Morphing1" , logger ) :  
        test_morphing1   () 
    with timing ("Morphing2" , logger ) :  
        test_morphing2   () 
    with timing ("Morphing3" , logger ) :  
        test_morphing3   () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
