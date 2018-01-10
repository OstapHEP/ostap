#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_models_2D.py
# Test module for ostap/fitting/models_2d.py
# - It tests various 2D-non-factrorizeable models 
# ============================================================================= 
""" Test module for ostap/fitting/models_2d.py
- It tests various 2D-non-factrorizeable models 
"""
# ============================================================================= 

import ROOT, random
import ostap.fitting.roofit 
import ostap.fitting.models as     Models 
from   ostap.core.core      import Ostap, VE, dsID
from   ostap.logger.utils   import rooSilent 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_models_2D' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make simple test mass 
m_x     = ROOT.RooRealVar ( 'mass_x' , 'Some test mass(X)' , 3 , 3.2 )
m_y     = ROOT.RooRealVar ( 'mass_y' , 'Some test mass(Y)' , 3 , 3.2 )

## book very simple data set
varset  = ROOT.RooArgSet  ( m_x , m_y )
dataset = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset )  

## fill it with uniform (x,y) 
m = VE(3.100,0.015**2)
w = VE(3.100,0.100**2) 
for i in xrange(0,5000) :
    
    m_x.setVal ( random.uniform ( *m_x.minmax() )  
                 m_y.setVal ( random.uniform ( *m_y.minmax() )
    
    dataset.add ( varset  )

logger.info ( 'DataSet: %s' % dataset )

models = set() 


# =============================================================================
## Positive polynomial in X and Y 
# =============================================================================
def test_polypos2D() :
    
    logger.info ('Test PolyPos2D_pdf: positive polynomial in 2D' )
    model = Models.PolyPos2D_pdf ( 'P2D'  ,
                                   m_x    ,
                                   m_y    ,
                                   nx = 2 ,
                                   ny = 2 )
    
    with rooSilent() : 
        result, f = model.fitTo ( dataset ) 

    result, f = model.fitTo ( dataset , silent = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :
        logger.info ( 'Bernstein Coefficients:\n%s' % model.pars() )  
            
       
    models.add ( model ) 
    
# =============================================================================
## Positive *SYMMETRIC* polynomial in X and Y 
# =============================================================================
def test_polypossym2D() :
    
    logger.info ('Test PolyPos2Dsym_pdf: Symmetric positive polynomial' )
    model = Models.PolyPos2Dsym_pdf ( 'P2Ds ',
                                      m_x    ,
                                      m_y    ,
                                      n = 2  )

    with rooSilent() : 
        result, f = model.fitTo ( dataset ) 

    result, f = model.fitTo ( dataset , silent = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :
        logger.info ( 'Bernstein Coefficients:\n%s' % model.pars() )  
            
       
    models.add ( model ) 

# =============================================================================
## product of phase space factors, modulated by positive polynomial in X and Y 
# =============================================================================
def test_pspol2D() : 
    
    logger.info ('Test PSPol2D_pdf: product of phase space factors, modulated by positive polynomial in X and Y' )
    
## "fictive phase space"
    psx   = Ostap.Math.PhaseSpaceNL ( 0 , 10 , 2 , 10 )
    psy   = Ostap.Math.PhaseSpaceNL ( 0 , 10 , 2 , 10 )
    model = Models.PSPol2D_pdf ( 'PS2D',
                                 m_x    , m_y  ,
                                 psx    , psy  , 
                                 nx = 2 , ny = 2 )
    
    with rooSilent() : 
        result, f = model.fitTo ( dataset ) 
        
    result, f = model.fitTo ( dataset , silent = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :
        logger.info ( 'Bernstein Coefficients:\n%s' % model.pars() )  
            
       
    models.add ( model ) 
 
# =============================================================================
## *SYMMETRIC* product of phase space factors, modulated by positive polynomial in X and Y 
# =============================================================================
def test_pspolsym2D() : 
    logger.info ('Test PSPol2Dsym_pdf: *SYMMETRIC* product of phase space factors, modulated by positive polynomial in X and Y ')
    
    ## "fictive phase space"
    ps    = Ostap.Math.PhaseSpaceNL ( 0 , 10 , 2 , 10 )
    model = Models.PSPol2Dsym_pdf ( 'PS2Ds',
                                    m_x   , m_y  ,
                                    ps    , n = 2 )
    
    with rooSilent() : 
        result, f = model.fitTo ( dataset ) 
        
    result, f = model.fitTo ( dataset , silent = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :
        logger.info ( 'Bernstein Coefficients:\n%s' % model.pars() )  
            
       
    models.add ( model )

# =============================================================================
## exponential times phase space factor, modulated by positive polynomial in X and Y 
# =============================================================================
def test_expopspol2D() :
    
    logger.info ('Test ExpoPSPol2D_pdf: Exponential times phase space factor, modulated by positive polynomial in X and Y ')
    
    ## "fictive phase space"
    psy   = Ostap.Math.PhaseSpaceNL ( 0 , 10 , 2 , 10 )
    model = Models.ExpoPSPol2D_pdf ( 'EPS',
                                     m_x    , m_y  ,
                                     psy    , 
                                     nx = 2 , ny = 2 )
    
    with rooSilent() : 
        result, f = model.fitTo ( dataset ) 
        
    result, f = model.fitTo ( dataset , silent = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :
        logger.info ( 'Bernstein Coefficients:\n%s' % model.pars() )  
            
       
    models.add ( model )

# =============================================================================
## exponential times exponential modulated by positive polynomial in X and Y 
# =============================================================================
def test_expopol2D() :
    
    logger.info ('Test ExpoPol2D_pdf: Exponential times exponential modulated by positive polynomial in X and Y ')
## "fictive phase space"
    model = Models.ExpoPol2D_pdf ( 'EP',
                                   m_x    , m_y  ,
                                   nx = 2 , ny = 2 )
    

    with rooSilent() : 
        result, f = model.fitTo ( dataset ) 
        
    result, f = model.fitTo ( dataset , silent = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :
        logger.info ( 'Bernstein Coefficients:\n%s' % model.pars() )  
            
       
    models.add ( model )

# =============================================================================
## symmetric exponential times exponential modulated by positive polynomial in X and Y 
# =============================================================================
def test_expopolsym2D() : 
    logger.info ('Test ExpoPol2Dsym_pdf: symmetric exponential times exponential modulated by positive polynomial in X and Y')

    model = Models.ExpoPol2Dsym_pdf ( 'EPs',
                                      m_x   ,
                                      m_y  ,
                                      n = 2 )
    
    with rooSilent() : 
        result, f = model.fitTo ( dataset ) 
        
    result, f = model.fitTo ( dataset , silent = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :
        logger.info ( 'Bernstein Coefficients:\n%s' % model.pars() )  
            
       
    models.add ( model )

# =============================================================================
## 2D-spline 
# =============================================================================
def test_spline2D() :
    
    logger.info ('Test Spline2D_pdf : 2D-spline')
    
    s1 = Ostap.Math.BSpline  ( m_x.xmin(), m_x.xmax() , 1 , 2 ) 
    s2 = Ostap.Math.BSpline  ( m_y.xmin(), m_y.xmax() , 1 , 2 ) 
    s3 = Ostap.Math.Spline2D ( s1 , s2 )
    
    model = Models.Spline2D_pdf ( 'S2D' , m_x , m_y, s3 )
    
    model = Models.ExpoPol2Dsym_pdf ( 'EPs',
                                      m_x   ,
                                      m_y  ,
                                      n = 2 )
    
    with rooSilent() : 
        result, f = model.fitTo ( dataset ) 
        
    result, f = model.fitTo ( dataset , silent = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
       
    models.add ( model )

# =============================================================================
## Symetric 2D-spline 
# =============================================================================
def test_splinesym2D() :

    logger.info ('Test Spline2Dsym_pdf: Symetric 2D-spline')
    
    ss    = Ostap.Math.BSpline     ( m_x.xmin(), m_x.xmax() , 1 , 2 ) 
    ss3   = Ostap.Math.Spline2DSym ( ss)
    
    model = Models.Spline2Dsym_pdf ( 'SS2D' , m_x , m_y, ss3 )

    with rooSilent() : 
        result, f = model.fitTo ( dataset ) 
        
    result, f = model.fitTo ( dataset , silent = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
       
    models.add ( model )


# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :
    logger.info('Saving all objects into DBASE')
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( name = 'Save everything to DBASE'), DBASE.tmpdb() as db : 
        db['m_x'     ] = m_x
        db['m_y'     ] = m_y
        db['vars'    ] = varset
        db['dataset' ] = dataset
        db['models'  ] = models
        db.ls() 
    
# =============================================================================
if '__main__' == __name__ :

    test_polypos2D    () 
    test_polypossym2D () 
    test_pspol2D      () 
    test_pspolsym2D   () 
    test_expopspol2D  () 
    test_expopol2D    () 
    test_expopolsym2D () 
    test_spline2D     () 
    test_splinesym2D  () 
    
    ## check finally that everything is serializeable:
    test_db ()          
    
# =============================================================================
# The END 
# ============================================================================= 
