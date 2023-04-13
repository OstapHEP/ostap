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
from   ostap.utils.utils           import vrange
from   builtins                    import range
from   ostap.utils.timing          import timing 
from   ostap.plotting.canvas       import use_canvas
from   ostap.utils.utils           import wait 
from   ostap.fitting.morphing_pdf  import ( MorphingN1_pdf ,
                                            MorphingN2_pdf  ,
                                            LinearMorph_pdf ) 
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
mass.setBins ( 10000  ,'cache' )

h1   = ROOT.TH1D ('h' ,'' , 200 , 0 , 20 )
N    = 10000
for i in range ( N ) :
    h1.Fill ( random.gauss ( 10 , 2.5 ) ) 

# ============================================================================
def test_morphingL () :

    logger = getLogger ('test_morphingL')    
    if root_info < ( 6 , 23 )  : 
        logger.warning( 'Test is disabled for ROOT version %s' % ROOT.gROOT.GetVersion() )
        return
    
    pdf1 = Models.Gauss_pdf ( 'GL1'  ,
                              xvar  = mass ,
                              mean  = ROOT.RooFit.RooConst ( 10 ) ,
                              sigma = ROOT.RooFit.RooConst ( 1  ) )
    pdf2 = Models.Flat1D    ( xvar = mass )
    
    ## create morphing PDF 
    pdf  = LinearMorph_pdf ( 'ML' , pdf1 , pdf2 ) 

    amin , amax = pdf.alpha.minmax ()     
    for alpha in vrange ( amin , amax  , 10 ) :
        pdf.alpha = alpha
        logger.info ( 'alpha= %s' % alpha ) 
        with wait ( 0.2 ) , use_canvas ( 'test_morphingL' ) :
            pdf.draw()
            
    r , f = pdf.fitHisto ( h1 , draw = False , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_morphingL' ) :
        r , f = pdf.fitHisto ( h1 , draw = True , nbins = 100 , silent = True )
        logger.info ( 'Morphing: \n%s' % r.table ( prefix = "# " ) ) 

# ============================================================================
def test_morphing1 () :

    logger = getLogger ('test_morphing1')    
    if root_info < ( 6 , 23 ) : 
        logger.warning( 'Test is disabled for ROOT version %s' % ROOT.gROOT.GetVersion() )
        return

    shapes = {}

    mean        = 10 
    smin , smax = 0.5 , 5.0
    
    for i , sigma in  enumerate ( vrange ( smin , smax , 10 ) ) :

        gauss = Models.Gauss_pdf ( 'G1_%d' % i ,
                                   xvar  = mass ,
                                   mean  = ROOT.RooFit.RooConst ( mean  ) ,
                                   sigma = ROOT.RooFit.RooConst ( sigma ) )
        shapes [ sigma ] = gauss
        
    ## create morphing PDF 
    pdf  = MorphingN1_pdf ( 'M1' , shapes , xvar =  mass )
        
    for mu in vrange ( smin , smax , 6 ) :
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
    
    if root_info < ( 6 , 23 ) : 
        logger.warning( 'Test is disabled for ROOT version %s' % ROOT.gROOT.GetVersion() )
        return


    shapes = {}

    sigma       = 2.5
    mmin , mmax = 8.0 , 12.0
    
    for j , mean in  enumerate ( vrange ( mmin , mmax , 10 ) ) :
        gauss = Models.Gauss_pdf ( 'G2_%d' % j  , 
                                   xvar  = mass ,
                                   mean  = ROOT.RooFit.RooConst ( mean  ) ,
                                   sigma = ROOT.RooFit.RooConst ( sigma ) )
        shapes [ mean ] = gauss
        
    ## create morphing PDF 
    pdf  = MorphingN1_pdf ( 'M2', shapes , xvar    =  mass ) 
    
    r , f = pdf.fitHisto ( h1 , draw = False , silent = True )
    r , f = pdf.fitHisto ( h1 , draw = False , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_morphing2' ) :
        r , f = pdf.fitHisto ( h1 , draw = True  , nbins = 100 , silent = True )
        logger.info ( 'Morphing: \n%s' % r.table ( prefix = "# " ) ) 
    
# ============================================================================
def test_morphing3 () :
    
    logger = getLogger ('test_morphing3')
    
    if root_info < ( 6 , 23 ) : 
        logger.warning( 'Test is disabled for ROOT version %s' % ROOT.gROOT.GetVersion() )
        return


    shapes = {}
    
    smin , smax = 0.5 ,  5.0
    mmin , mmax = 8.0 , 12.0

    
    for i , sigma in  enumerate ( vrange ( smin , smax , 10 ) ) :
        for j , mean in  enumerate ( vrange ( mmin , mmax , 10 ) ) :
            gauss = Models.Gauss_pdf ( 'G3_%d_%d' % ( i , j ) , 
                                       xvar  = mass ,
                                       mean  = ROOT.RooFit.RooConst ( mean  ) ,
                                       sigma = ROOT.RooFit.RooConst ( sigma ) )
            shapes [ mean , sigma ] = gauss

    ## create morphing PDF 
    pdf  = MorphingN2_pdf ( 'M3' , shapes , xvar    =  mass )
            
    r , f = pdf.fitHisto ( h1 , draw = False , silent = True )
    r , f = pdf.fitHisto ( h1 , draw = False , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_morphing3' ) :
        r , f = pdf.fitHisto ( h1 , draw = True  , nbins = 100 , silent = True )
        logger.info ( 'Morphing: \n%s' % r.table ( prefix = "# " ) ) 
    

# =============================================================================
if '__main__' == __name__ :

    with timing ("MorphingL" , logger ) :  
        test_morphingL   () 
    with timing ("Morphing1" , logger ) :  
        test_morphing1   () 
    with timing ("Morphing2" , logger ) :  
        test_morphing2   () 
    with timing ("Morphing3" , logger ) :  
        test_morphing3   () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
