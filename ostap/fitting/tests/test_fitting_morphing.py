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
from   ostap.core.meta_info        import root_info 
from   ostap.utils.ranges          import vrange
from   ostap.utils.timing          import timing 
from   ostap.plotting.canvas       import use_canvas
from   ostap.utils.root_utils      import batch_env 
from   ostap.core.core             import roo_silent 
from   ostap.fitting.morphing_pdf  import ( MorphingN1_pdf ,
                                            MorphingN2_pdf  ,
                                            LinearMorph_pdf ) 
import ostap.fitting.models        as     Models
import ostap.fitting.roofit
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
## set batch from environment 
batch_env ( logger )
# =============================================================================
#
mass = ROOT.RooRealVar ( 'mass' , 'some mass' , 0 , 20 )
mass.setBins ( 10000  ,'cache' )
varset = ROOT.RooArgSet ( mass ) 
ds   = ROOT.RooDataSet ( 'ds' , '' , varset ) 

h1   = ROOT.TH1D ('h' ,'' , 200 , 0 , 20 )
N    = 10000
for i in range ( N ) :
    v = random.gauss ( 10 , 2.5 ) 
    h1.Fill ( v )
    if v in mass :
        mass.setVal ( v )
        ds.add ( varset )  
    
conf = { 'refit'     :   5 }

if (6,27) <= root_info : conf [ 'maxcalls' ] = 1000000
if (6,29) <= root_info : 
    conf [ 'minimizer'] = 'Minuit','migrad'
    conf [ 'hesse'    ] =  True 

# ============================================================================
def test_morphingL () :

    logger = getLogger ('test_morphingL')    
 
    
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
        with use_canvas ( 'test_morphingL' , wait = 0.5 ) :
            pdf.draw()
            
    with use_canvas ( 'test_morphingL' , wait = 1 ) :
        r , f = pdf.fitHisto ( h1 , draw = True , nbins = 100 , silent = True , **conf )
        logger.info ( 'Morphing: \n%s' % r.table ( prefix = "# " ) ) 

# ============================================================================
def test_morphing_N1s() :

    logger = getLogger ('test_morphing_N1s')    
  

    shapes = {}

    mean        = 10 
    smin , smax = 0.5 , 5.0
    
    for i , sigma in  enumerate ( vrange ( smin , smax , 20 ) ) :
        gauss = Models.Gauss_pdf ( 'G1_%d' % i ,
                                   xvar  = mass ,
                                   mean  = ROOT.RooFit.RooConst ( mean  ) ,
                                   sigma = ROOT.RooFit.RooConst ( sigma ) )
        shapes [ sigma ] = gauss
        
    ## create morphing PDF 
    pdf  = MorphingN1_pdf ( 'N1s' , shapes , xvar =  mass )
            
    with use_canvas ( 'test_morphing_N1s' , wait = 1 ) :
        r , f = pdf.fitHisto ( h1 , draw = True , nbins = 100 , silent = True , **conf )
        logger.info ( 'Morphing: \n%s' % r.table ( prefix = "# " ) ) 


# ============================================================================
def test_morphing_N1m () :
    
    logger = getLogger ('test_morphing_N1m')
    


    shapes = {}

    sigma       = 2.5
    mmin , mmax = 5.0 , 15.0
    
    for j , mean in  enumerate ( vrange ( mmin , mmax , 20 ) ) :
        gauss = Models.Gauss_pdf ( 'G2_%d' % j  , 
                                   xvar  = mass ,
                                   mean  = ROOT.RooFit.RooConst ( mean  ) ,
                                   sigma = ROOT.RooFit.RooConst ( sigma ) )
        shapes [ mean ] = gauss
        
    ## create morphing PDF 
    pdf  = MorphingN1_pdf ( 'Nm1', shapes , xvar    =  mass ) 
    
    with use_canvas ( 'test_morphing_N1m' , wait = 1 ) :
        r , f = pdf.fitHisto ( h1 , draw = True  , nbins = 100 , silent = True , **conf )
        logger.info ( 'Morphing: \n%s' % r.table ( prefix = "# " ) ) 
    
# ============================================================================
def test_morphing_N2 () :
    
    logger = getLogger ('test_morphing_N2')
    
  

    shapes = {}
    
    mmin , mmax = 8.0 , 12.0
    smin , smax = 0.5 ,  5.5
    
    for i , mean in  enumerate ( vrange ( mmin , mmax , 10 ) ) :
        for j , sigma in  enumerate ( vrange ( smin , smax , 10 ) ) :
            gauss = Models.Gauss_pdf ( 'G3_%d_%d' % ( i , j ) , 
                                       xvar  = mass ,
                                       mean  = ROOT.RooFit.RooConst ( mean  ) ,
                                       sigma = ROOT.RooFit.RooConst ( sigma ) )
            shapes [ mean , sigma ] = gauss
            ## logger.info ( 'Add [%2d,%2d] = [%6.3f,%6.2f] component ' % ( i , j , mean , sigma ) )
            
    ## create morphing PDF
    
    pdf  = MorphingN2_pdf ( 'N2' , shapes , xvar =  mass )
            
    with use_canvas ( 'test_morphing_N2' , wait = 1) :
        r , f = pdf.fitHisto ( h1 , draw = True , nbins = 100 , silent = True , **conf )
        logger.info ( 'Morphing: \n%s' % r.table ( prefix = "# " ) )

# =============================================================================
if '__main__' == __name__ :

    with timing ("MorphingL" , logger ) :  
        test_morphingL      () 
    with timing ("Morphing_N1s" , logger ) :  
        test_morphing_N1s   () 
    with timing ("Morphing_N1m" , logger ) :  
        test_morphing_N1m   () 
    with timing ("Morphing_N2" , logger ) :  
        test_morphing_N2    () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
