#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_shapes2.py
# ============================================================================= 
""" Test module for ostap/fitting/shapes2.py
- It tests various ``signal-like''/``peak-like'' shapes 
"""
# ============================================================================= 
from   __future__               import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.meta_info     import root_info
from   ostap.core.pyrouts    import hID, dsID, Ostap 
import ostap.fitting.models  as     Models
from   ostap.plotting.canvas import use_canvas
from   ostap.utils.timing    import timing 
from   ostap.utils.gsl       import gslCount 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_shapes2' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

xvar   = ROOT.RooRealVar ( 'x' , 'x-variable' , 0 , 10 )
yvar   = ROOT.RooRealVar ( 'y' , 'y-variable' , 0 ,  5 )

varset = ROOT.RooArgSet  ( xvar , yvar )
ds1    = ROOT.RooDataSet ( dsID()  , 'dataset' , varset )
ds2    = ROOT.RooDataSet ( dsID()  , 'dataset' , varset )


N1 = 300
N2 = 1000

## prepare some more or less random data
for i in range ( N1 ) :     
    xv = random.expovariate ( 1.0/5 )
    while not xv in xvar : xv = random.expovariate ( 1.0/5 )
    yv = random.expovariate ( 1.0/3 )
    while not yv in yvar : yv = random.expovariate ( 1.0/3 )

    xvar.setVal ( xv )
    yvar.setVal ( yv )
    ds1.add ( varset )
    
for i in range ( 3*N1 ) :
    xv = -100
    yv = -100
    while ( not xv in xvar ) or ( not yv in yvar ) : 
        x1 = random.gauss ( 0 , 3 )
        x2 = random.gauss ( 0 , 1 )
        xv = 5 + x1 + x2 
        yv = 3 + x1 - x2

    xvar.setVal ( xv )
    yvar.setVal ( yv )
    ds1.add ( varset )

## prepare some more or less random data
for i in range ( N2 ) :     
    xv = random.expovariate ( 1.0/5 )
    while not xv in xvar : xv = random.expovariate ( 1.0/5 )
    yv = random.expovariate ( 1.0/3 )
    while not yv in yvar : yv = random.expovariate ( 1.0/3 )

    xvar.setVal ( xv )
    yvar.setVal ( yv )
    ds2.add ( varset )
    
for i in range ( 3*N2  ) :
    xv = -100
    yv = -100
    while ( not xv in xvar ) or ( not yv in yvar ) : 
        x1 = random.gauss ( 0 , 3 )
        x2 = random.gauss ( 0 , 1 )
        xv = 5 + x1 + x2 
        yv = 3 + x1 - x2

    xvar.setVal ( xv )
    yvar.setVal ( yv )
    ds2.add ( varset )


ds1, ds2 = ds2, ds1 
## ds2 = ds1

flat   = Models.Flat2D( xvar = xvar , yvar = yvar , name = 'Flat2D' )    
pdfs   = []
frames = []

# ===========================================================================================
# (1) rely on histograms 
# ===========================================================================================
def test_shapes2_histos () :

    logger = getLogger('test_shapes2_histos')
    
    pdfs1 = [] 
    for nx,ny  in [   ( 10 , 10 ) , ( 15 , 15 ) , ( 20 , 20 ) , ( 22 , 22 ) ] :
        
        h2 = ROOT.TH2F ( hID() , 'histograms' ,
                         nx , xvar.getMin() , xvar.getMax() ,
                         ny , yvar.getMin() , yvar.getMax() )
        
        ds1.project ( h2 , 'y : x' )
        
        ## rely on RooHistPdf 
        pdf = Models.H2D_pdf       ( 'R2D_%02d%02d' % ( nx , ny ) , h2 , xvar = xvar , yvar = yvar , order = 2 )
        pdfs1.append ( pdf )
        
        ## rely on Histo2D_pdf 
        shape = Ostap.Math.Histo2D ( h2 , 3 , 3 )          
        pdf   = Models.Histo2D_pdf ( 'H2D_%02d%02d' % ( nx , ny ) , histo = shape , xvar = xvar , yvar = yvar )
        pdfs1.append ( pdf )


    fit_pdfs = [ Models.Sum2D ( [ p , flat ] ) for p in pdfs1 ]

    for pdf in fit_pdfs :
        
        with gslCount() , timing ( pdf.name , logger = logger ) : 
            
            pdf.fractions = 0.99 ,   
            r , _ = pdf.fitTo ( ds2 , silent = True )
            
            pdf.fractions = 0.99 , 
            r , _ = pdf.fitTo ( ds2 , silent = True )
            
            logger.info ( '%s:\n%s' % ( pdf.name , r.table ( title = pdf.name , prefix = '# ' ) ) )
            
            with use_canvas ( 'x-view: %s' %  pdf.name )  :
                fx    = pdf.draw1 ( ds2 , nbins = 50 )
                frames.append ( fx )
                
            with use_canvas ( 'y-view: %s' %  pdf.name )  :
                fy    = pdf.draw2 ( ds2 , nbins = 50 )
                frames.append ( fy )
                    
# =========================================================================================
## (2) use non-parametric polynomials and convert them to the shapes
# =========================================================================================
def test_shapes2_poly () :

    logger = getLogger('test_shapes2_poly')


    if root_info < (6,18) :
        logger.warning ( "Test is disabled for ROOT version %s" % str ( root_info ) )
        return 

    pdfs2 = []
    
    for nx,ny  in [  ( 5 , 5 ) ,  ( 7 , 7 ) , ( 10 , 10 ) , ( 12 , 12 ) , ( 15 , 15 ) ] :
        
        p2 = Ostap.Math.LegendreSum2 ( nx , ny ,  xvar.getMin() , xvar.getMax() , yvar.getMin() , yvar.getMax() )
        
        ds1.project ( p2 , 'y : x' )
        
        ## rely on Shape2D_pdf
        shape = p2 
        pdf   = Models.Shape2D_pdf ( 'P2D_%02d%02d' %  (nx ,ny ) , shape = shape , xvar = xvar , yvar = yvar )
        pdfs2.append ( pdf )

    fit_pdfs = [ Models.Sum2D ( [ p , flat ] ) for p in pdfs2 ]

    for pdf in fit_pdfs :
        
        with gslCount() , timing ( pdf.name , logger = logger ) :
            
            pdf.fractions = 0.99 ,   
            r , _ = pdf.fitTo ( ds2 , silent = True )
            
            pdf.fractions = 0.99 , 
            r , _ = pdf.fitTo ( ds2 , silent = True )
            
            logger.info ( '%s:\n%s' % ( pdf.name , r.table ( title = pdf.name , prefix = '# ' ) ) )
            
            with use_canvas ( 'x-view: %s' %  pdf.name )  :
                fx    = pdf.draw1 ( ds2 , nbins = 50 )
                frames.append ( fx )
                
            with use_canvas ( 'y-view: %s' %  pdf.name )  :
                fy    = pdf.draw2 ( ds2 , nbins = 50 )
                frames.append ( fy )
                

# =========================================================================================
## (3) Use RooKeysPDF 
# =========================================================================================
def test_shapes2_keys () :

    logger = getLogger('test_shapes2_keys')

    pdfs3 = [] 

    pdf = Models.RooKeys2D_pdf ( 'RKeys' , xvar = xvar , yvar = yvar , data = ds1 )
    pdfs3.append ( pdf ) 
    
    fit_pdfs = [ Models.Sum2D ( [ p , flat ] ) for p in pdfs3 ]
    
    for pdf in fit_pdfs :
        
        with gslCount() , timing ( pdf.name , logger = logger ) :
            
            pdf.fractions = 0.99 ,   
            r , _ = pdf.fitTo ( ds2 , silent = True )
            
            pdf.fractions = 0.99 , 
            r , _ = pdf.fitTo ( ds2 , silent = True )
            
            logger.info ( '%s:\n%s' % ( pdf.name , r.table ( title = pdf.name , prefix = '# ' ) ) )
            
            with use_canvas ( 'x-view: %s' %  pdf.name )  :
                fx    = pdf.draw1 ( ds2 , nbins = 50 )
                frames.append ( fx )
                
            with use_canvas ( 'y-view: %s' %  pdf.name )  :
                fy    = pdf.draw2 ( ds2 , nbins = 50 )
                frames.append ( fy )
                
    
# =============================================================================
if '__main__' == __name__ :
    
    test_shapes2_histos  () 
    test_shapes2_poly    () 
    test_shapes2_keys    () 
    
# =============================================================================
##                                                                      The END 
# ============================================================================= 

