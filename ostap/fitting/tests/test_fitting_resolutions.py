#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_resolutions.py
# Test module for ostap/fitting/resolution.py
# - It tests various resoltuion models
# ============================================================================= 
"""Test module for ostap/fitting/resolution.py
- It tests various resoltuion models
"""
# ============================================================================= 
from   __future__        import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ROOT, random
import ostap.fitting.roofit 
from   ostap.core.core       import VE, dsID
try:
    from builtins import range
except ImportError:
    from __builtin__ import range
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_resolutions' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make simple test mass 
mass     = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , -3 , 3 )

## book very simple data set
varset0  = ROOT.RooArgSet  ( mass )
dataset0 = ROOT.RooDataSet ( dsID() , 'Test Data set-0' , varset0 )  

mmin = mass.getMin()
mmax = mass.getMax()

## fill it 
m1 = VE(0,0.1**2)
for i in range(0,20000) :
    mass.setVal  ( m1.gauss () )
    dataset0.add ( varset0     )
    
m2 = VE(0,0.2**2)
for i in range(0,10000) :
    mass.setVal  ( m2.gauss () )
    dataset0.add ( varset0     )
    
logger.info ('DATASET %s' % dataset0 )

models = set() 


# =============================================================================
## Single gauss
# =============================================================================
def test_gauss () :
    logger.info ('Test ResoGauss:  single Gaussian resolution model' )
    from   ostap.fitting.resolution import ResoGauss 
    reso_gauss = ResoGauss( 'G1' , mass ,  0.1 )
    reso_gauss.sigma.release()

    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso_gauss. fitTo ( dataset0 )
        result, frame = reso_gauss. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoGauss:  RMS         %s ' % reso_gauss.rms    () )  
        logger.info ( 'ResoGauss:  FWHM        %s ' % reso_gauss.fwhm   () )  
        logger.info ( 'ResoGauss:  Resolution: %s ' % result.sigma_G1.ve() )  

    models.add ( reso_gauss )
    
# =============================================================================
## Double gauss
# =============================================================================
def test_2gauss () :
    logger.info ('Test ResoGauss2:  double Gaussian resolution model' )
    from   ostap.fitting.resolution import ResoGauss2
    reso_2gauss = ResoGauss2( 'G2' , mass ,  0.1 )
    reso_2gauss.sigma   .release()
    reso_2gauss.scale   .release()
    reso_2gauss.fraction.release()
    
    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso_2gauss. fitTo ( dataset0 )
        result, frame = reso_2gauss. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoGauss2: RMS        %s ' % reso_2gauss.rms          () )  
        logger.info ( 'ResoGauss2: FWHM       %s ' % reso_2gauss.fwhm         () )  
        logger.info ( 'ResoGauss2: Resolution %s ' % result.sigma_G2       .ve() )  
        logger.info ( 'ResoGauss2: Scale      %s ' % result.SigmaScale_G2  .ve() )  
        logger.info ( 'ResoGauss2: Fraction   %s ' % result.CoreFraction_G2.ve() )  
        
    models.add ( reso_2gauss )

        
# =============================================================================
## Symmetric Apolonious
# =============================================================================
def test_apo2 () :
    
    logger.info ('Test ResoApo2:  symmetric Apolonios resolution model' )
    from   ostap.fitting.resolution import ResoApo2
    reso_apo2 = ResoApo2( 'A2' , mass ,  0.1 )
    reso_apo2.sigma .release()
    reso_apo2.beta  .release()
    
    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso_apo2. fitTo ( dataset0 )
        result, frame = reso_apo2. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoApo2:   RMS        %s ' % reso_apo2.rms          () )  
        logger.info ( 'ResoApo2:   FWHM       %s ' % reso_apo2.fwhm         () )  
        logger.info ( 'ResoApo2:   Resolution %s ' % result.sigma_A2     .ve() )  
        logger.info ( 'ResoApo2:   Beta       %s ' % result.ResoBeta_A2  .ve() )  
    
    models.add ( reso_apo2)

        
# =============================================================================
## Symmetric double-sided Crystal Ball 
# =============================================================================
def test_cb2 () :
    
    logger.info ('Test ResoCB2: symmetric double-sided Crystal Ball resolution model' )
    from   ostap.fitting.resolution import ResoCB2
    reso_cb2 = ResoCB2( 'CB' , mass ,  0.1 )
    reso_cb2.sigma .release()
    reso_cb2.alpha .release()
    reso_cb2.n     .release()
    
    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso_cb2. fitTo ( dataset0 )
        result, frame = reso_cb2. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoCB2:    RMS        %s ' % reso_cb2.rms           () )  
        logger.info ( 'ResoCB2:    FWHM       %s ' % reso_cb2.fwhm          () )  
        logger.info ( 'ResoCB2:    Resolution %s ' % result.sigma_CB     .ve() )  
        logger.info ( 'ResoCB2:    Alpha      %s ' % result.ResoAlpha_CB .ve() )  
        logger.info ( 'ResoCB2:    N          %s ' % result.ResoN_CB     .ve() )  
        
    models.add ( reso_cb2)


# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :
    logger.info('Saving all objects into DBASE')
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( name = 'Save everything to DBASE'), DBASE.tmpdb() as db : 
        db['mass,vars'] = mass, varset0
        db['dataset'  ] = dataset0
        db['models'   ] = models
        db.ls() 
        
# =============================================================================
if '__main__' == __name__ :

    test_gauss  () ## single Gaussian resolution model
    test_2gauss () ## double Gaussian resolution model
    test_apo2   () ## double Gaussian resolution model
    test_cb2    () ## double Gaussian resolution model
    
    ## check finally that everything is serializeable:
    test_db ()          

# =============================================================================
# The END 
# ============================================================================= 
