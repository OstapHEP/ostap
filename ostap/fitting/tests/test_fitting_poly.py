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
from   __future__        import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import 
# ============================================================================= 
import ROOT, random
import ostap.fitting.roofit 
import ostap.fitting.models as     Models 
from   ostap.core.core      import cpp, VE, dsID
from   ostap.logger.utils   import rooSilent 
from   ostap.utils.timing   import timing 
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
    
    logger.info("Test  Poly(4)-Distribution")
    model = Models.PolyPos_pdf('P4'  , x , power = 4 )
    
    result,f  = model.fitTo ( dataset , silent = True )  
    model.draw ( dataset )        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        for phi in model.phis : 
            logger.info ( "\tPoly4:        phi= %-17s " % phi.ve() ) 
            
    models.add ( model ) 

# =============================================================================
## Test  monotonic Poly(4)-Distribution
# =============================================================================
def test_monotonic4 () :
    
    logger.info("Test  monotonic Poly(4)-Distribution")
    model = Models.Monotonic_pdf('M4'  , x , power = 4 , increasing = False  )
    
    with rooSilent() : 
        result,f  = model.fitTo ( dataset )  
        result,f  = model.fitTo ( dataset )  
        model.draw ( dataset )        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        for phi in model.phis : 
            logger.info ( "\tMonotonic4:  phi= %-17s " % phi.ve() ) 
            
    models.add ( model ) 

# =============================================================================
## Test  convex Poly(4)-Distribution
# =============================================================================
def test_convex4() :
    
    logger.info("Test  convex Poly(4)-Distribution")
    model = Models.Convex_pdf('C4'  , x , power = 4 , increasing = False , convex = True  )
    
    result,f  = model.fitTo ( dataset , silent = True )  
    model.draw ( dataset )        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        for phi in model.phis : 
            logger.info ( "\tConvex4:      phi= %-17s " % phi.ve() ) 
            
    models.add ( model ) 

# =============================================================================
## Test  Poly(2)*Expo -Distribution
# =============================================================================
def test_expopoly2() : 
    logger.info("Test  Poly(2)*Expo -Distribution")
    model = Models.Bkg_pdf('P2e'  , x , power = 2 )
    model.tau.fix(-1.25)

    result,f  = model.fitTo ( dataset , silent = True )
    model.tau.release() 
    result,f  = model.fitTo ( dataset , silent = True )  
    model.draw ( dataset )        
    
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        logger.info ( "\tTau:          tau= %-17s " %  result( model.tau ) [0] ) 
        for phi in model.phis : 
            logger.info ( "\tExpoPoly2:    phi= %-17s " % phi.ve() ) 
            
    models.add ( model ) 

# =============================================================================
## Test positive spline: order 3 with 2 inner knots
# =============================================================================
def test_pspline () : 
    
    logger.info ("Test positive spline of order 3 with 2 inner knots ")
    ## define spline 
    spline = cpp.Ostap.Math.PositiveSpline( x.xmin() , x.xmax() , 2 , 3 )
    ## build the model
    model  = Models.PSpline_pdf ( 'S3' , x , spline )

    ## fit it! 
    result,f  = model.fitTo ( dataset , silent = True )
    model.draw ( dataset )        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        for phi in model.phis : 
            logger.info ( "\tpSpline:      phi= %-17s " % phi.ve() ) 
            
    models.add ( model ) 

# =============================================================================
## Test positive decreasing: order 3 with 2 inner knots ")
# =============================================================================
def test_mspline () :
    
    logger.info ("Test positive decreasing spline of order 3 with 2 inner knots ")
    spline = cpp.Ostap.Math.MonotonicSpline( x.xmin() , x.xmax() , 2 , 3 , False )
    model  = Models.MSpline_pdf ( 'S3' , x , spline )

    ## fit it! 
    result,f  = model.fitTo ( dataset , silent = True )
    model.draw ( dataset )        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        for phi in model.phis : 
            logger.info ( "\tmSpline:      phi= %-17s " % phi.ve() ) 
            
    models.add ( model ) 

# =============================================================================
## Test positive decreasing convex: order 3 with 2 inner knots 
# =============================================================================
def test_cspline () :
    
    logger.info ("Test positive decreasing convex spline of order 3 with 2 inner knots")
    spline = cpp.Ostap.Math.ConvexSpline( x.xmin() , x.xmax() , 2 , 3 , False , True )
    model  = Models.CSpline_pdf ( 'C3' , x , spline )

    ## fit it! 
    result,f  = model.fitTo ( dataset , silent = True )
    model.draw ( dataset )        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        for phi in model.phis : 
            logger.info ( "\tcSpline:      phi= %-17s " % phi.ve() ) 
            
    models.add ( model ) 


# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    from ostap.core.meta_info import root_version_int 
    if root_version_int >= 62200 :
        logger.warning("test_db: test is disabled for ROOT version %s" % root_version_int )
        return 
    
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
    with timing ( "Monotonic4" , logger ) : 
        test_monotonic4 () ## Monotonic polynomial (4)
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
