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
from   __future__               import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ROOT, random
import ostap.fitting.roofit 
from   ostap.core.core          import VE, dsID
from   builtins                 import range
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
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
    
logger.info ('DATASET\n%s' % dataset0.table ( prefix = "# " ) ) 

models = set() 


# =============================================================================
## Single gauss
# =============================================================================
def test_gauss () :

    logger = getLogger ( 'test_gauss' )
    
    logger.info ('Test ResoGauss:  single Gaussian resolution model' )
    from   ostap.fitting.resolution import ResoGauss 
    reso_gauss = ResoGauss( 'G1' , mass ,  0.1 )
    reso_gauss.sigma.release()

    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso_gauss. fitTo ( dataset0 )
        with wait ( 1 ) , use_canvas ( 'test_gauss' ) : 
            result, frame = reso_gauss. fitTo ( dataset0 , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoGauss:  RMS         %s ' % reso_gauss.rms    () )  
        logger.info ( 'ResoGauss:  FWHM        %s ' % reso_gauss.fwhm   () )
        logger.info ( "ResoGauss: fit results\n%s" % result.table ( title = 'Gaussian resolution model' , prefix = '# ' ) )

    models.add ( reso_gauss )
    
# =============================================================================
## Double gauss
# =============================================================================
def test_2gauss () :
    
    logger = getLogger ( 'test_2gauss' )

    logger.info ('Test ResoGauss2:  double Gaussian resolution model' )
    from   ostap.fitting.resolution import ResoGauss2
    reso_2gauss = ResoGauss2( 'G2' , mass ,  0.1 )
    reso_2gauss.sigma   .release()
    reso_2gauss.scale   .release()
    reso_2gauss.fraction.release()
    
    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso_2gauss. fitTo ( dataset0 )
        with wait ( 1 ) , use_canvas ( 'test_2gauss' ) : 
            result, frame = reso_2gauss. fitTo ( dataset0 , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoGauss2: RMS        %s ' % reso_2gauss.rms          () )  
        logger.info ( 'ResoGauss2: FWHM       %s ' % reso_2gauss.fwhm         () )
        logger.info ( "ResoGauss2: fit results\n%s" % result.table ( title = 'double-Gaussian resolution model' , prefix = '# ' ) )
        
    models.add ( reso_2gauss )

        
# =============================================================================
## Symmetric Apollonios
# =============================================================================
def test_apo2 () :
    
    logger = getLogger ( 'test_apo2' )

    logger.info ('Test ResoApo2:  symmetric Apollonios resolution model' )
    from   ostap.fitting.resolution import ResoApo2
    reso_apo2 = ResoApo2( 'A2' , mass ,  0.1 )
    reso_apo2.sigma .release()
    reso_apo2.beta  .release()
    
    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso_apo2. fitTo ( dataset0 )
        with wait ( 1 ) , use_canvas ( 'test_apo2' ) : 
            result, frame = reso_apo2. fitTo ( dataset0 , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoApo2:   RMS        %s ' % reso_apo2.rms          () )  
        logger.info ( 'ResoApo2:   FWHM       %s ' % reso_apo2.fwhm         () )
        logger.info ( "ResoApo2:   fit results\n%s" % result.table ( title = 'symmetric Apollonios resolution model' , prefix = '# ' ) )
    
    models.add ( reso_apo2)

        
# =============================================================================
## Symmetric double-sided Crystal Ball 
# =============================================================================
def test_cb2 () :
    
    logger = getLogger ( 'test_cb2' )

    logger.info ('Test ResoCB2: symmetric double-sided Crystal Ball resolution model' )
    from   ostap.fitting.resolution import ResoCB2
    reso_cb2 = ResoCB2( 'CB' , mass ,  0.1 )
    reso_cb2.sigma .release()
    reso_cb2.alpha .release()
    reso_cb2.n     .release()
    
    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso_cb2. fitTo ( dataset0 )
        with wait ( 1 ) , use_canvas ( 'test_cb2' ) : 
            result, frame = reso_cb2. fitTo ( dataset0 , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoCB2:    RMS        %s ' % reso_cb2.rms           () )  
        logger.info ( 'ResoCB2:    FWHM       %s ' % reso_cb2.fwhm          () )
        logger.info ( "ResoCB2:    fit results\n%s" % result.table ( title = 'symmetric double-sided Crystal Ball resolution model' , prefix = '# ' ) )
        
    models.add ( reso_cb2)


# =============================================================================
## Hyperbolic secant 
# =============================================================================
def test_sech () :
    
    logger = getLogger ( 'test_sech' )

    logger.info ('Test ResoSech: hyperbolic secant resolution model' )
    from   ostap.fitting.resolution import ResoSech
    reso_sech = ResoSech ( 'Sech' , mass ,  0.1 )
    reso_sech.sigma .release()
    
    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso_sech. fitTo ( dataset0 )
        with wait ( 1 ) , use_canvas ( 'test_sech' ) : 
            result, frame = reso_sech. fitTo ( dataset0 , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoSech:   RMS        %s ' % reso_sech.rms          () )  
        logger.info ( 'ResoSech:   FWHM       %s ' % reso_sech.fwhm         () )
        logger.info ( "ResoSech:   fit results\n%s" % result.table ( title = 'hyperbolic secant/sech resolution model' , prefix = '# ' ) )
        
    models.add ( reso_sech )

# =============================================================================
## Logistic 
# =============================================================================
def test_logistic () :
    
    logger = getLogger ( 'test_logistic' )

    logger.info ('Test ResoLogistic: logistic (sech-squared) resolution model' )
    from   ostap.fitting.resolution import ResoLogistic
    reso_log = ResoLogistic ( 'Logistic' , mass ,  0.1 )
    reso_log.sigma .release()
    
    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso_log. fitTo ( dataset0 )
        result, frame = reso_log. fitTo ( dataset0 , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoLog :   RMS        %s ' % reso_log.rms          () )  
        logger.info ( 'ResoLog :   FWHM       %s ' % reso_log.fwhm         () )
        logger.info ( "ResoLog :   fit results\n%s" % result.table ( title = 'Logistic/sech-squared resolution model' , prefix = '# ' ) )
        
    models.add ( reso_log )


# =============================================================================
## symmetric Bukin
# =============================================================================
def test_bukin () :
    
    logger = getLogger ( 'test_bukin' )

    logger.info ('Test ResoBukin: symmetric Bukin resolution model' )
    from   ostap.fitting.resolution import ResoBukin
    reso = ResoBukin ( 'Bukin' , mass ,  0.1 ,rho = (0, 0 , 10 ) )
    reso.sigma .release()
    reso.rho   .release()
    
    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso. fitTo ( dataset0 )
        with wait ( 1 ) , use_canvas ( 'test_bukin' ) : 
            result, frame = reso. fitTo ( dataset0 , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoBukin:  RMS        %s ' % reso.rms          () )  
        logger.info ( 'ResoBukin:  FWHM       %s ' % reso.fwhm         () )
        logger.info ( "ResoBukin: fit results\n%s" % result.table ( title = 'symmetric Bukin resolution model' , prefix = '# ' ) )
        
    models.add ( reso )

# =============================================================================
## symmetric Johnson's SU 
# =============================================================================
def test_johnsonSU () :
    
    logger.info ('Test JohnsonSU: symmetric JohnsonSU  resolution model' )
    from   ostap.fitting.resolution import ResoJohnsonSU 
    reso = ResoJohnsonSU ( 'JohnsonSU' , mass ,
                           delta = ( 1.7 , 1.e-6 , 1000 ) ,
                           lambd = ( 0.2 , 1.e-6 , 1000 ) ) 
    
    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso. fitTo ( dataset0 )
        with wait ( 1 ) , use_canvas ( 'test_johnsonSU' ) : 
            result, frame = reso. fitTo ( dataset0 , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoJohnsonSU:  RMS        %s ' % reso.rms          () )  
        logger.info ( 'ResoJohnsonSU:  FWHM       %s ' % reso.fwhm         () )
        logger.info ( "ResoJohnsonSU: fit results\n%s" % result.table ( title = 'symmetric JohnsonSU resolution model' , prefix = '# ' ) )
        
    models.add ( reso )



# =============================================================================
## Sinh-Asinh
# =============================================================================
def test_sinhasinh () :
    
    logger = getLogger ( 'test_sinhasinh' )

    logger.info ('Test SinhAsinh: symmetric SinhAsinh resolution model' )
    from   ostap.fitting.resolution import ResoSinhAsinh
    reso = ResoSinhAsinh ( 'SinhAsinh' , mass ,  delta = ( 0.7 , 1.e-5 , 1000 ) )
    
    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso. fitTo ( dataset0 )
        with wait ( 1 ) , use_canvas ( 'test_sinhasinh' ) : 
            result, frame = reso. fitTo ( dataset0 , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoSinhAsinh:  RMS        %s ' % reso.rms          () )  
        logger.info ( 'ResoSinhAsinh:  FWHM       %s ' % reso.fwhm         () )
        logger.info ( "ResoSinhAsinh: fit results\n%s" % result.table ( title = 'symmetric SinhAsinh resolution model' , prefix = '# ' ) )
        
    models.add ( reso )


# =============================================================================
## Hyperbolic 
# =============================================================================
def test_hyperbolic () :
    
    logger = getLogger ( 'test_hyperbolic' )

    logger.info ('Test Hyperbolic: symmetric Hyperbolic resolution model' )
    from   ostap.fitting.resolution import ResoHyperbolic
    reso = ResoHyperbolic ( 'Hyperbolic' , mass ,  zeta = ( 1.0 , 1.e-5 , 1.e+5 ) )
    
    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso. fitTo ( dataset0 )
        with wait ( 1 ) , use_canvas ( 'test_hyperbolic' ) : 
            result, frame = reso. fitTo ( dataset0 , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoHyperbolic:  RMS        %s ' % reso.rms          () )  
        logger.info ( 'ResoHyperbolic:  FWHM       %s ' % reso.fwhm         () )
        logger.info ( "ResoHyperbolic: fit results\n%s" % result.table ( title = 'symmetric Hyperbolic resolution model' , prefix = '# ' ) )
        
    models.add ( reso )


# =============================================================================
## Generalized Hyperbolic 
# =============================================================================
def test_genhyperbolic () :
    
    logger = getLogger ( 'test_genhyperbolic' )

    logger.info ('Test Hyperbolic: symmetric generalised Hyperbolic resolution model' )
    from   ostap.fitting.resolution import ResoGenHyperbolic
    reso = ResoGenHyperbolic ( 'GenHyperbolic' , mass ,  zeta = ( 1.0 , 1.e-5 , 1.e+5 ) , lambd =  (-100,100) )
    
    from   ostap.logger.utils   import rooSilent
    with rooSilent() : 
        result, frame = reso. fitTo ( dataset0 )
        with wait ( 1 ) , use_canvas ( 'test_genhyperbolic' ) : 
            result, frame = reso. fitTo ( dataset0 , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        logger.info ( 'ResoGenHyperbolic:  RMS        %s ' % reso.rms          () )  
        logger.info ( 'ResoGenHyperbolic:  FWHM       %s ' % reso.fwhm         () )
        logger.info ( "ResoGenHyperbolic: fit results\n%s" % result.table (
            title = 'symmetric generalised Hyperbolic resolution model' , prefix = '# ' ) )
        
    models.add ( reso )



# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' )

    logger.info('Saving all objects into DBASE')
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db : 
        db['mass,vars'] = mass, varset0
        db['dataset'  ] = dataset0
        for m in models : db['model %s' % m.name ] = m
        db['models'   ] = models
        db.ls() 
        
# =============================================================================
if '__main__' == __name__ :

    with timing ("Gauss"     , logger ) :  
        test_gauss      () ## single Gaussian resolution model
        
    with timing ("2-Gauss"   , logger ) :  
        test_2gauss     () ## double Gaussian resolution model
        
    with timing ("Apo2"      , logger ) :  
        test_apo2       () ## symmetric Apollonios resoltuion model
        
    with timing ("CB2"       , logger ) :  
        test_cb2        () ## double-sided Crystal Ball resoltuion model
        
    with timing ("Sech"      , logger ) :  
        test_sech       () ## hyperbolic secant resolution model

    with timing ("Logistic"  , logger ) :  
        test_logistic   () ## logistic resolution model
        
    with timing ("Bukin"     , logger ) :  
        test_bukin      () ## Bukin resolution model
    
    with timing ("SinhAsinh" , logger ) :  
        test_sinhasinh  () ## SinhAsinh resolution model

    with timing ("JohnsonSU" , logger ) :  
        test_johnsonSU  () ## JohnsonSU resolution model

    with timing ("Hyperbolic" , logger ) :  
        test_hyperbolic  () ## Hyperbolic resolution model
    
    with timing ("GenHyperbolic" , logger ) :  
        test_genhyperbolic  () ## generalised Hyperbolic resolution model
    
    ## check finally that everything is serializeable:
    with timing ("Save to DB"    , logger ) :  
        test_db ()          


# =============================================================================
##                                                                      The END 
# ============================================================================= 
