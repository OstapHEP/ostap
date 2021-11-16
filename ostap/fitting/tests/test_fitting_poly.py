#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_poly.py
# Test module for ostap/fitting/models.py
# - It tests various ``background-like'' smooth polynomial shapes
# ============================================================================= 
""" Test module for ostap/fitting/models.py
- It tests various ``background-like'' smooth polynomial shapes
"""
# ============================================================================= 
from   __future__               import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import 
# ============================================================================= 
import ROOT, random
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import cpp, VE, dsID
from   ostap.logger.utils       import rooSilent 
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'ostap/fitting/tests/test_fitting_poly' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( 'Test for some polynomials models from ostap')
# =============================================================================
## make simple test 
mass     = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 2 , 10 )
x        = mass 

## book very simple data set
varset   = ROOT.RooArgSet  ( mass )


events = 10000

logger.debug('Make a test data using Gamma-Distribution')
m_gamma0 = Models.GammaDist_pdf( 'GD0' , x )
m_gamma0.k    .setVal( 2 )
m_gamma0.theta.setVal( 1 )

dataset = m_gamma0.pdf.generate ( varset , events ) 

 
models = set()

# =============================================================================
## Test  Poly(4)-Distribution
# =============================================================================
def test_poly4 () :

    logger = getLogger ( 'test_poly4' ) 
    logger.info("Test  Poly(4)-Distribution")
    model = Models.PolyPos_pdf('P4'  , x , power = 3 )
    
    result,f  = model.fitTo ( dataset , silent = True , refit = 5 )
    with wait ( 1 ) , use_canvas ( 'test_poly4' ) : 
        model.draw ( dataset )        
        
    title = 'Positive polynomial'
    logger.info('%s \n%s' % ( title , result.table ( title = title , prefix = '# ' ) ) ) 
            
    models.add ( model ) 

# =============================================================================
## Test  monotonic Poly(5)-Distribution
# =============================================================================
def test_monotonic5 () :

    logger = getLogger ( 'test_monotonic5' )
    logger.info("Test  monotonic Poly(5)-Distribution") 
    model = Models.Monotonic_pdf( 'M4'  , x , power = 4, increasing = False  )
    
    result,f  = model.fitTo ( dataset , silent = True , refit = 5 )  
    with wait ( 1 ) , use_canvas ( 'test_monotonic5' ) : 
        model.draw ( dataset )        
        
    title = 'Positive decreasing polynomial'
    logger.info('%s \n%s' % ( title , result.table ( title = title , prefix = '# ' ) ) ) 

    models.add ( model ) 

# =============================================================================
## Test  convex Poly(4)-Distribution
# =============================================================================
def test_convex4() :
    
    logger = getLogger ( 'test_convex4' )

    logger.info("Test  convex Poly(4)-Distribution")
    model = Models.Convex_pdf('C4'  , x , power = 4 , increasing = False , convex = True  )
    
    result,f  = model.fitTo ( dataset , silent = True , refit = 5 )  
    with wait ( 1 ) , use_canvas ( 'test_convex4' ) : 
        model.draw ( dataset )        
        
    title = 'Positive decreasing convex polynomial'
    logger.info('%s \n%s' % ( title , result.table ( title = title , prefix = '# ' ) ) ) 
            
    models.add ( model ) 

# =============================================================================
## Test  Poly(2)*Expo -Distribution
# =============================================================================
def test_expopoly2() :

    logger = getLogger ( 'test_expopoly2' )

    logger.info("Test  Poly(2)*Expo -Distribution")
    model = Models.Bkg_pdf('P2e'  , x , power = 1 )
    model.tau.fix(-1.25)

    result,f  = model.fitTo ( dataset , silent = True )
    model.tau.release() 
    result,f  = model.fitTo ( dataset , silent = True , refit = 5 )    
    with wait ( 1 ) , use_canvas ( 'test_expopoly2' ) : 
        model.draw ( dataset )        
    
    title = 'Expo*polynom'
    logger.info('%s \n%s' % ( title , result.table ( title = title , prefix = '# ' ) ) ) 
            
    models.add ( model ) 

# =============================================================================
## Test positive spline: order 3 with 2 inner knots
# =============================================================================
def test_pspline () : 
    
    logger = getLogger ( 'test_pspline' )

    logger.info ("Test positive spline of order 3 with 1 inner knot ")
    ## define spline 
    spline = cpp.Ostap.Math.PositiveSpline( x.xmin() , x.xmax() , 1 , 3 )
    ## build the model
    model  = Models.PSpline_pdf ( 'S3' , x , spline )

    ## fit it! 
    result,f  = model.fitTo ( dataset , silent = True , refit = 5  )
    with wait ( 1 ) , use_canvas ( 'test_pspline' ) : 
        model.draw ( dataset )        
        
    title = 'Positive spline'
    logger.info('%s \n%s' % ( title , result.table ( title = title , prefix = '# ' ) ) ) 
            
    models.add ( model ) 

# =============================================================================
## Test positive decreasing: order 3 with 2 inner knots ")
# =============================================================================
def test_mspline () :
    
    logger = getLogger ( 'test_mspline' )

    logger.info ("Test positive decreasing spline of order 3 with 1 inner knot ")
    spline = cpp.Ostap.Math.MonotonicSpline( x.xmin() , x.xmax() , 1 , 3 , False )
    model  = Models.MSpline_pdf ( 'M3' , x , spline )

    ## fit it! 
    result,f  = model.fitTo ( dataset , silent = True , refit = 5 )
    with wait ( 1 ) , use_canvas ( 'test_mspline' ) : 
        model.draw ( dataset )        
        
    title = 'Positive descreasing spline'
    logger.info('%s\n%s' % ( title , result.table ( title = title , prefix = '# ' ) ) )

    models.add ( model ) 

# =============================================================================
## Test positive decreasing convex: order 3 with 2 inner knots 
# =============================================================================
def test_cspline () :
    
    logger = getLogger ( 'test_cspline' )
    logger.info ("Test positive decreasing convex spline of order 3 with 1 inner knot")
    spline = cpp.Ostap.Math.ConvexSpline( x.xmin() , x.xmax() , 1 , 3 , increasing = False , convex = True )
    model  = Models.CSpline_pdf ( 'C3' , x , spline )

    ## fit it! 
    result,f  = model.fitTo ( dataset , silent = True , refit = 5 )
    with wait ( 1 ) , use_canvas ( 'test_cspline' ) : 
        model.draw ( dataset )        

    title = 'Positive descreasing convex spline'
    logger.info('%s \n%s' % ( title , result.table ( title = title , prefix = '# ' ) ) ) 
            
    models.add ( model ) 


# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' )
    logger.info('Saving all objects into DBASE')    
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( name = 'Save everything to DBASE'), DBASE.tmpdb() as db : 
        db['mass,vars'] = mass, varset
        db['dataset'  ] = dataset
        db['models'   ] = models
        db.ls() 

# =============================================================================
if '__main__' == __name__ :

    with timing ( "Poly4" , logger ) : 
        test_poly4       () ## Polynomial (4)
    with timing ( "Monotonic5" , logger ) : 
        test_monotonic5 () ## Monotonic polynomial (5)
    with timing ( "Convex4"   , logger ) : 
        test_convex4     () ## Convex polynomial (4)
    with timing ( "ExpoP2"    , logger ) : 
        test_expopoly2   () ## Exponent times positive polynomial (2)
    with timing ( "p-Spline"  , logger ) :         
        test_pspline     () ## Positive spline of order 3 with two knots 
    with timing ( "m-Spline"  , logger ) :         
        test_mspline     () ## Positive monotonic spline of order 3 with two knots 
    with timing ( "c-Spline"  , logger ) :         
        test_cspline     () ## Positive monotonic convex spline of order 3 with two knots 

    ## check finally that everything is serializeable:
    with timing ( "Save to DB"  , logger ) :         
        test_db          ()          

# =============================================================================
##                                                                      The END 
# ============================================================================= 
