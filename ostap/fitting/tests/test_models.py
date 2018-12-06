#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_models.py
# Test module for ostap/fitting/models.py
# - It tests various ``signal-like''/``peak-like'' shapes
# ============================================================================= 
""" Test module for ostap/fitting/models.py
- It tests various ``signal-like''/``peak-like'' shapes 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ROOT, random
import ostap.fitting.roofit 
import ostap.fitting.models as     Models 
from   ostap.core.core      import cpp, VE, dsID
from   ostap.logger.utils   import rooSilent 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_models' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make simple test mass 
mass     = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 3.0 , 3.2 )

## book very simple data set
varset0  = ROOT.RooArgSet  ( mass )
dataset0 = ROOT.RooDataSet ( dsID() , 'Test Data set-0' , varset0 )  

mmin,mmax = mass.minmax()

## fill it 
m = VE(3.100,0.015**2)
for i in xrange(0,5000) :
    mass.value = m.gauss () 
    dataset0.add ( varset0 )

for i in xrange(0,500) :
    mass.value = random.uniform ( mmin , mmax ) 
    dataset0.add ( varset0   )

logger.info ('DATASET %s' % dataset0 )

models = set() 

## signal component 
signal_gauss = Models.Gauss_pdf ( name  = 'Gauss'    , ## the name 
                                  xvar  = mass       , ## the variable 
                                  mean  = m.value () , ## mean value (fixed)
                                  sigma = m.error () ) ## sigma      (fixed)

# =============================================================================
## gauss PDF
# =============================================================================
def test_gauss() :
    
    logger.info ('Test Gauss_pdf:  simple Gaussian signal' )
    
    ## release the sigma of signal:
    signal_gauss.mean.release()
    
    ## simple fit with gaussian only 
    result = signal_gauss . fitTo ( dataset0 , silent = True )
    
    ## construct composite model: signal + background 
    model_gauss = Models.Fit1D( signal     = signal_gauss ,
                                background = Models.Bkg_pdf ('BkgGauss', xvar = mass , power = 0 ) )
    model_gauss.background.tau.fix(0)
    
    with rooSilent() : 
        result, frame = model_gauss . fitTo ( dataset0 )
        result, frame = model_gauss . fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print result 
    else :     
        logger.info( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) )
        logger.info( 'Mean   & Sigma      are: %-28s & %-28s ' % ( result ( 'mean_Gauss')[0] , result( 'sigma_Gauss' )[0] ) )

    models.add ( model_gauss )
    signal_gauss.mean.fix ( m.value() )

# =============================================================================
## CrystalBall PDF
# =============================================================================
def test_crystalball () :
    
    logger.info ('Test CrystalBall_pdf: Crystal Ball  function' )
    
    ## composite model: signal + background 
    model_cb = Models.Fit1D (
        signal     = Models.CrystalBall_pdf ( name  = 'CB'     , ## the name 
                                              xvar  = mass     , ## the variable   
                                              alpha = (2,1,5)  , ## tail parameter
                                              n     = (3,1,9)  , ## tail parameter 
                                              sigma = signal_gauss.sigma ,   ## reuse sigma from gauss
                                              mean  = signal_gauss.mean  ) , ## reuse mean  from gauss 
        background = Models.Bkg_pdf ('BkgCB', xvar = mass , power = 0 )
        ) 
    model_cb.signal.n.fix(8) 
    with rooSilent() : 
        result, frame = model_cb. fitTo ( dataset0 )
        model_cb.signal.alpha.release()
        result, frame = model_cb. fitTo ( dataset0 )
        result, frame = model_cb. fitTo ( dataset0 )
    
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else : 
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean   & Sigma      are: %-28s & %-28s ' % ( result ( 'mean_Gauss')[0] , result( 'sigma_Gauss' )[0] ) ) 
        logger.info ( 'Alpha  & n          are: %-28s & %-28s ' % ( result ( model_cb.signal.alpha ) [ 0 ] , result ( model_cb.signal.n ) [ 0 ] ) )
        
    models.add ( model_cb )
                      

# =============================================================================
## right side CrystalBall PDF
# =============================================================================
def test_crystalball_RS () : 
    logger.info ('Test CrystalBallRS_pdf:  right-side Crystal Ball function' )
    model_cbrs = Models.Fit1D (
        signal = Models.CrystalBallRS_pdf ( name  = 'CBRS' , 
                                            xvar  = mass               ,
                                            sigma = signal_gauss.sigma ,
                                            alpha = (1.5, 0.5 , 3.0)   ,
                                            n     = (5,1,10)           , 
                                            mean  = signal_gauss.mean  ) ,
        background = Models.Bkg_pdf ('BkgCBRS', xvar = mass , power = 0 )
        )
    
    model_cbrs.S.value  = 5000
    model_cbrs.B.value  =  500
    
    with rooSilent() : 
        result, frame = model_cbrs. fitTo ( dataset0 )
        model_cbrs.signal.alpha.release()
        result, frame = model_cbrs. fitTo ( dataset0 )
        result, frame = model_cbrs. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result
    else : 
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean   & Sigma      are: %-28s & %-28s ' % ( result ( 'mean_Gauss')[0] , result( 'sigma_Gauss' )[0] ) ) 
        
    models.add ( model_cbrs  )

# =============================================================================
## double sided CrystalBall PDF
# =============================================================================
def test_crystalball_DS () :
    
    logger.info ('Test CrystalBallDS_pdf: double-sided Crystal Ball function' )
    model_cbds = Models.Fit1D (
        signal = Models.CB2_pdf ( name   = 'CB2'              , 
                                  xvar   = mass               ,
                                  nL     = 10                 , 
                                  nR     = 10                 , 
                                  alphaL = (1.5,0.5,3)        , 
                                  alphaR = (1.4,0.5,3)        , 
                                  sigma  = signal_gauss.sigma ,  
                                  mean   = signal_gauss.mean  ) ,
        background = Models.Bkg_pdf ('BkgCBDS', xvar = mass , power = 0 )
        )

    model_cbds.S.value  = 5000
    model_cbds.B.value  =  500
    
    with rooSilent() : 
        result, frame = model_cbds. fitTo ( dataset0 )
        model_cbds.signal.aL.release()
        model_cbds.signal.aR.release()
        result, frame = model_cbds. fitTo ( dataset0 )
        model_cbds.signal.aL.fix(1.5) 
        model_cbds.signal.aR.fix(1.5)    
        result, frame = model_cbds. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean   & Sigma      are: %-28s & %-28s ' % ( result ( 'mean_Gauss')[0] , result( 'sigma_Gauss' )[0] ) ) 

    models.add ( model_cbds  )

# =============================================================================
## Needham PDF
# =============================================================================
def test_needham() :
    
    logger.info ('Test Needham_pdf: Crystal Ball with alpha=f(sigma)' )
    model_matt = Models.Fit1D (
        signal = Models.Needham_pdf ( name  = 'Matt'             , 
                                      xvar  = mass               ,
                                      sigma = signal_gauss.sigma ,  
                                      mean  = signal_gauss.mean  ) ,
        background = Models.Bkg_pdf ('BkgCBDS', xvar = mass , power = 0 )
        )
    
    with rooSilent() : 
        result, frame = model_matt. fitTo ( dataset0 )
        model_matt.signal.mean .release()
        model_matt.signal.sigma.release()
        result, frame = model_matt. fitTo ( dataset0 )
        result, frame = model_matt. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else : 
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean   & Sigma      are: %-28s & %-28s ' % ( result ( 'mean_Gauss')[0] , result( 'sigma_Gauss' )[0] ) ) 
        
    models.add ( model_matt  )


# ==========================================================================
## Apolonios
# ==========================================================================
def test_apolonios () :
    
    logger.info ('Test Apolonios_pdf: Modified gaussian with power-law and exponential tails' ) 
    model_apolonios = Models.Fit1D (
        signal = Models.Apolonios_pdf ( name  = 'APO', 
                                        xvar  = mass ,
                                        mean  = signal_gauss.mean    ,
                                        sigma = signal_gauss.sigma ,
                                        b     =  1 ,
                                        n     = 10 ,
                                        alpha =  3 ) ,
        background = Models.Bkg_pdf ('BkgAPO', xvar = mass , power = 0 )
        )
    
    model_apolonios.S.setVal(5000)
    model_apolonios.B.setVal( 500)
    
    with rooSilent() : 
        result, frame = model_apolonios. fitTo ( dataset0 )
        result, frame = model_apolonios. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else : 
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean   & Sigma      are: %-28s & %-28s ' % ( result ( 'mean_Gauss')[0] , result( 'sigma_Gauss' )[0] ) ) 

        
    models.add ( model_apolonios )

# ==========================================================================
## Apolonios2
# ==========================================================================
def test_apolonios2() :
    
    logger.info ('Test Apolonios2_pdf: modified Gaussian with exponential tails' ) 
    model_apolonios2 = Models.Fit1D (
        signal = Models.Apolonios2_pdf ( name = 'AP2' , 
                                         xvar      = mass ,
                                         mean      = signal_gauss.mean  ,
                                         sigma     = signal_gauss.sigma ,
                                         beta      =  ( 0.5 , 2 )       ,
                                         asymmetry = 0 ) ,
        background = Models.Bkg_pdf ('BkgAPO2', xvar = mass , power = 0 )
        )
    
    model_apolonios2.signal.mean.fix( m.value() )    
    model_apolonios2.S.value  = 5000
    model_apolonios2.B.value  =  500
    model_apolonios2.signal.sigma.release() 
    
    with rooSilent() :
        result, frame = model_apolonios2. fitTo ( dataset0 )
        model_apolonios2.signal.asym.release ()
        result, frame = model_apolonios2. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else : 
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean   & Sigma      are: %-28s & %-28s ' % ( result ( 'mean_Gauss')[0] , result( 'sigma_Gauss' )[0] ) ) 

        
    models.add ( model_apolonios2 )


# =============================================================================
## Bifurcated gauss PDF
# =============================================================================
def test_bifurcated () :
    
    logger.info ('Test BifurcatedGauss_pdf: Bifurcated Gaussian' )
    
    signal_bifurcated = Models.BifurcatedGauss_pdf ( name = 'BfGau' ,
                                                     mean  = signal_gauss.mean  ,
                                                     sigma = signal_gauss.sigma ,
                                                     xvar  = mass    )
    signal_bifurcated . asym  . setVal ( 0          )
    
    model_bifurcated = Models.Fit1D(
        signal     = signal_bifurcated       ,
        background = Models.Bkg_pdf ('BkgBFG', xvar  = mass , power = 0 )) 
    
    model_bifurcated.B.setVal (  500 )
    model_bifurcated.S.setVal ( 6000 )
    with rooSilent() : 
        result, frame = model_bifurcated . fitTo ( dataset0 )
        result, frame = model_bifurcated . fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :     
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean   & Sigma      are: %-28s & %-28s ' % ( result ( 'mean_Gauss')[0] , result( 'sigma_Gauss' )[0] ) )
        
    models.add ( model_bifurcated  )


# =============================================================================
## Bifurcated gauss PDF
# =============================================================================
def test_2gauss () :
    
    logger.info ('Test DoubleGauss_pdf: Double Gaussian' )
    
    signal_2gauss = Models.DoubleGauss_pdf ( name = 'Gau2' ,
                                             mean  = signal_gauss.mean  ,
                                             sigma = signal_gauss.sigma ,
                                             xvar  = mass               ,
                                             fraction = 0.9             ,
                                             scale    = 1.2             )
    
    model_2gauss = Models.Fit1D(
        signal     = signal_2gauss      ,
        background = Models.Bkg_pdf ('Bkg22G', xvar = mass , power = 0 )) 
    
    model_2gauss.B.setVal (  500 )
    model_2gauss.S.setVal ( 6000 )
    with rooSilent() : 
        result, frame = model_2gauss. fitTo ( dataset0 )
        signal_2gauss.fraction.release() 
        result, frame = model_2gauss. fitTo ( dataset0 )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :     
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean   & Sigma      are: %-28s & %-28s ' % ( result ( 'mean_Gauss')[0] , result( 'sigma_Gauss' )[0] ) )
        
    models.add ( model_2gauss  )


# =============================================================================
## GenGaussV1
# =============================================================================
def test_gengauss_v1 () :
    logger.info ('Test GenGaussV1_pdf: Generalized Gaussian V1' ) 
    model_gauss_gv1 = Models.Fit1D (
        signal = Models.GenGaussV1_pdf ( name = 'Gv1' , 
                                         xvar = mass  ,
                                         mean = signal_gauss.mean ) ,
        background = Models.Bkg_pdf ('BkgGGV1', xvar = mass , power = 0 )) 
    
    model_gauss_gv1.signal.beta .fix(2)
    model_gauss_gv1.signal.mean .fix( m.value() ) 
    model_gauss_gv1.S.setVal(5000)
    model_gauss_gv1.B.setVal( 500)
    
    with rooSilent() : 
        result, frame = model_gauss_gv1. fitTo ( dataset0 )
        model_gauss_gv1.signal.alpha.release()
        result, frame = model_gauss_gv1. fitTo ( dataset0 )
        model_gauss_gv1.signal.mean .release() 
        result, frame = model_gauss_gv1. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :     
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        

    models.add ( model_gauss_gv1  )

# =============================================================================
## GenGaussV2
# =============================================================================
def test_gengauss_v2 () : 
    logger.info ('Test GenGaussV2_pdf: Generalized Gaussian function V2' ) 
    model_gauss_gv2 = Models.Fit1D (
        signal = Models.GenGaussV2_pdf ( name = 'Gv2' , 
                                         xvar = mass  ,
                                         mean = signal_gauss.mean ) ,
        background = Models.Bkg_pdf ('BkgGGV2', xvar = mass , power = 0 )) 
    
    model_gauss_gv2.signal.kappa.fix(0)
    
    with rooSilent() : 
        result, frame = model_gauss_gv2. fitTo ( dataset0 )
        model_gauss_gv2.signal.mean.release() 
        model_gauss_gv2.S.setVal(5000)
        model_gauss_gv2.B.setVal( 500)
        result, frame = model_gauss_gv2. fitTo ( dataset0 )
        ##model_gauss_gv2.signal.kappa.release() 
        model_gauss_gv2.S.setVal(5000)
        model_gauss_gv2.B.setVal( 500)
        result, frame = model_gauss_gv2. fitTo ( dataset0 )
        result, frame = model_gauss_gv2. fitTo ( dataset0 )
        model_gauss_gv2.S.setVal(5000)
        model_gauss_gv2.B.setVal( 500)
        result, frame = model_gauss_gv2. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :     
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        
    models.add ( model_gauss_gv2  )

# =============================================================================
## SkewGauss
# =============================================================================
def test_skewgauss() :
    
    logger.info ('Test SkewGauss_pdf: Skew Gaussian function' ) 
    model_gauss_skew = Models.Fit1D (
        signal = Models.SkewGauss_pdf ( name = 'GSk' , 
                                        xvar = mass  , mean = signal_gauss.mean ) ,
        background = Models.Bkg_pdf ('BkgSkG', xvar = mass , power = 0 )) 
    
    model_gauss_skew.signal.alpha.fix(0)
    model_gauss_skew.S.setVal(5000)
    model_gauss_skew.B.setVal( 500)
    
    with rooSilent() : 
        result, frame = model_gauss_skew. fitTo ( dataset0 )
        result, frame = model_gauss_skew. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :     
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 

    models.add ( model_gauss_skew  )

# =============================================================================
## QGauss
# =============================================================================
def test_qgauss () :
    logger.info ('Test QGaussian_pdf: q-Gaussian' ) 
    model_qgauss = Models.Fit1D (
        signal = Models.QGaussian_pdf ( name = 'qG'  , 
                                        xvar = mass  ,
                                        q    = (1,0.7,1.2), 
                                        mean = signal_gauss.mean   ,
                                        scale = signal_gauss.sigma ) ,
        background = Models.Bkg_pdf ('BkgQG', xvar = mass , power = 0 )) 
    
    s = model_qgauss.signal
    s.scale = 0.015
    
    with rooSilent() : 
        result, frame = model_qgauss. fitTo ( dataset0 )
        model_qgauss.signal.scale.release()
        result, frame = model_qgauss. fitTo ( dataset0 )
        model_qgauss.signal.q .release() 
        result, frame = model_qgauss. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :     
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Scale  & Q          are: %-28s & %-28s ' % ( result ( s.scale     )[0] , result( s.q           )[0] ) ) 
        

    models.add ( model_qgauss )


# =============================================================================
## Bukin
# =============================================================================
def test_bukin() :
    
    logger.info ('Test Bukin_pdf: Bukin function: skew gaussian core + exponenial/gaussian  tails' ) 
    model_bukin = Models.Fit1D (
        signal = Models.Bukin_pdf ( name  = 'Bukin' ,
                                    xvar  = mass    ,
                                    xi    = 0    ,
                                    rhoL  = 0    ,
                                    rhoR  = 0    , 
                                    mean  = signal_gauss.mean  , 
                                    sigma = signal_gauss.sigma ) ,
        background = Models.Bkg_pdf ('BkgBK', xvar = mass , power = 0 )) 
    
    model_bukin.signal.mean .fix  ( m.value() )
    model_bukin.signal.sigma.fix  ( m.error() )
    model_bukin.S.setVal(5000)
    model_bukin.B.setVal( 500)
    
    with rooSilent() : 
        result, frame = model_bukin. fitTo ( dataset0 )
        model_bukin.signal.xi  .release()     
        result, frame = model_bukin. fitTo ( dataset0 )
        model_bukin.signal.rhoL.release()     
        result, frame = model_bukin. fitTo ( dataset0 )
        model_bukin.signal.rhoR.release()     
        result, frame = model_bukin. fitTo ( dataset0 )
        model_bukin.signal.mean .release() 
        model_bukin.signal.sigma.release() 
        result, frame = model_bukin. fitTo ( dataset0 )
    
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :     
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean   & Sigma      are: %-28s & %-28s ' % ( result ( 'mean_Gauss')[0] , result( 'sigma_Gauss' )[0] ) )
        
    models.add ( model_bukin  )

# =============================================================================
## StudentT
# =============================================================================
def test_studentT () :
    
    logger.info ('Test StudentT_pdf: Student-t distribution' ) 
    model_student = Models.Fit1D (
        signal = Models.StudentT_pdf ( name = 'ST' , 
                                       xvar = mass ,
                                       mean = signal_gauss.mean ) ,
        background = Models.Bkg_pdf ('BkgST', xvar = mass , power = 0 )) 
    
    model_student.signal.n    .setVal(20)
    model_student.signal.sigma.setVal(0.013)
    model_student.S.setVal(5000)
    model_student.B.setVal( 500)
    
    with rooSilent() : 
        result, frame = model_student. fitTo ( dataset0 )
        result, frame = model_student. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :     
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean                is : %-28s '         %   result ( 'mean_Gauss')[0]  )
        
    models.add ( model_student  )

# =============================================================================
## Bifurcated StudentT
# =============================================================================
def test_bifstudentT(): 
    logger.info ('Test bifurcated StudentT_pdf: bifurcated Student-t' ) 
    model = Models.Fit1D (
        signal = Models.BifurcatedStudentT_pdf ( name  = 'BfST' , 
                                                 xvar  = mass   ,
                                                 nL    = 25     ,
                                                 nR    = 25     ,                                                 
                                                 mean  = signal_gauss.mean   , 
                                                 sigma = signal_gauss.sigma  ) ,
        background = Models.Bkg_pdf ('BkgST2', xvar = mass , power = 0 )) 
    
    signal = model.signal 
    model.S.setVal(5000)
    model.B.setVal( 500)
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 )
        signal.nL   .release()
        signal.nR   .release()
        signal.sigma.release()
        result, frame = model. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :     
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean                 is: %-28s ' %  result ( signal.mean  )[0] )
        logger.info ( 'Sigma                is: %-28s ' %  result ( signal.sigma )[0] )
        logger.info ( 'Asymmetry            is: %-28s ' %  result ( signal.asym  )[0] )
        logger.info ( 'n(L)                 is: %-28s ' %  result ( signal.nL    )[0] )
        logger.info ( 'n(R)                 is: %-28s ' %  result ( signal.nR    )[0] )
        
    models.add ( model )

# =============================================================================
## Test  SinhAsinh-Distribution
# =============================================================================
def test_sinhasinh() :
    
    logger.info("Test  SinhAsinh-Distribution")
    model = Models.Fit1D (
        signal = Models.SinhAsinh_pdf( 'SASH'                   ,
                                       xvar = mass              , 
                                       mean = signal_gauss.mean ) ,
        background = Models.Bkg_pdf ('BkgSAS', xvar = mass , power = 0 )) 
    
    signal = model.signal
    
    model.S.setVal(5000)
    model.B.setVal( 500)
    
    signal.mu      .setVal (  3.10  )
    signal.sigma   .setVal (  0.015 ) 
    signal.epsilon .setVal (  0.021 ) 
    signal.delta   .setVal (  1.0   ) 

    with rooSilent() : 
        result,f  = model.fitTo ( dataset0 )  
        result,f  = model.fitTo ( dataset0 )  
        signal.delta.release()
        result,f  = model.fitTo ( dataset0 )  
        signal.epsilon.release()
        result,f  = model.fitTo ( dataset0 )  
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print result
    else :
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mu                   is: %-28s ' %  result ( signal.mu      )[0] )
        logger.info ( 'Sigma                is: %-28s ' %  result ( signal.sigma   )[0] )
        logger.info ( 'Epsilon              is: %-28s ' %  result ( signal.epsilon )[0] )
        logger.info ( 'delta                is: %-28s ' %  result ( signal.delta   )[0] )

    models.add ( model )


# =============================================================================
## Test  JohnsonSU-Distribution
# =============================================================================
def test_johnsonSU () :
    
    logger.info("Test  JohnsonSU-Distribution")
    model = Models.Fit1D (
        signal = Models.JohnsonSU_pdf( 'JSU'                    ,
                                       xvar = mass              , 
                                       xi   = signal_gauss.mean ) ,
        background = Models.Bkg_pdf ('BkgJSU', xvar = mass , power = 0 )) 
    
    signal = model.signal
    
    model.S.setVal(5000)
    model.B.setVal( 500)
    
    with rooSilent() : 
        result,f  = model.fitTo ( dataset0 )  
        result,f  = model.fitTo ( dataset0 )  
        signal.lambd .release()
        signal.delta.release()
        result,f  = model.fitTo ( dataset0 )  
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print result
    else :
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Xi                   is: %-28s ' %  result ( signal.xi     )[0] )
        logger.info ( 'Lambda               is: %-28s ' %  result ( signal.lambd  )[0] )
        logger.info ( 'Delta                is: %-28s ' %  result ( signal.delta  )[0] )
        logger.info ( 'Gamma                is: %-28s ' %  result ( signal.gamma  )[0] )

    models.add ( model )

# =============================================================================
## Test  ATLAS
# =============================================================================
def test_atlas () :
    
    logger.info("Test  ATLAS: Modified Gaussian, used by ATLAS/Zeus")
    model = Models.Fit1D (
        signal = Models.Atlas_pdf( 'ATLAS'                  ,
                                   xvar = mass              , 
                                   mean = signal_gauss.mean ) ,
        background = Models.Bkg_pdf ('BkgATLAS', xvar = mass , power = 0 )) 
    
    signal = model.signal
    model.S.setVal(5000)
    model.B.setVal( 500)
    
    with rooSilent() : 
        result,f  = model.fitTo ( dataset0 )  
        result,f  = model.fitTo ( dataset0 )  
        signal.mean  .release()
        signal.sigma .release()
        result,f  = model.fitTo ( dataset0 )  
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print result
    else :
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean                 is: %-28s ' %  result ( signal.mean  )[0] )
        logger.info ( 'Sigma                is: %-28s ' %  result ( signal.sigma )[0] )

    models.add ( model )

# =============================================================================
## Test  SECH
# =============================================================================
def test_sech() :
    
    logger.info("Test  SECH:  Sech(1/cosh) distribution")
    model = Models.Fit1D (
        signal = Models.Sech_pdf( 'SECH'                    ,
                                  xvar = mass              , 
                                  mean = signal_gauss.mean ) ,
        background = Models.Bkg_pdf ('BkgSECH', xvar = mass , power = 0 )) 
    
    signal = model.signal
    model.S.setVal(5000)
    model.B.setVal( 500)
    
    with rooSilent() : 
        result,f  = model.fitTo ( dataset0 )  
        result,f  = model.fitTo ( dataset0 )  
        signal.mean  .release()
        signal.sigma .release()
        result,f  = model.fitTo ( dataset0 )  
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print result
    else :
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean                 is: %-28s ' %  result ( signal.mean  )[0] )
        logger.info ( 'Sigma                is: %-28s ' %  result ( signal.sigma )[0] )

    models.add ( model )

# =============================================================================
## Test  LOGISTIC
# =============================================================================
def test_logistic () :
    
    logger.info("Test  LOGISTIC: Logistic distribution")
    model = Models.Fit1D (
        signal = Models.Logistic_pdf( 'LOGI'                    ,
                                      xvar = mass              , 
                                      mean = signal_gauss.mean ) ,
        background = Models.Bkg_pdf ('BkgLOGI', xvar = mass , power = 0 )) 
    
    signal = model.signal
    model.S.setVal(5000)
    model.B.setVal( 500)
    
    with rooSilent() : 
        result,f  = model.fitTo ( dataset0 )  
        result,f  = model.fitTo ( dataset0 )  
        signal.mean  .release()
        signal.sigma .release()
        result,f  = model.fitTo ( dataset0 )  
    
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print result
    else :
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean                 is: %-28s ' %  result ( signal.mean  )[0] )
        logger.info ( 'Sigma                is: %-28s ' %  result ( signal.sigma )[0] )

    models.add ( model )
# =============================================================================
## Voigt
# =============================================================================
def test_voigt () :
    
    logger.info ('Test Voigt_pdf: Breit-Wigner convoluted with Gauss' )
    model = Models.Fit1D (
        signal = Models.Voigt_pdf ( 'V' , 
                                    xvar  = mass                ,
                                    mean  = signal_gauss.mean   , 
                                    sigma = signal_gauss.sigma  ) , 
        background = Models.Bkg_pdf ('BkgV', xvar = mass , power = 0 )) 
    
    signal = model.signal
    signal.sigma.fix ( m.error() )
    signal.gamma.fix ( 0.002     )
    signal.mean .fix() 
    
    model.B.setVal( 500)
    model.S.setVal(5000)
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 )
        result, frame = model. fitTo ( dataset0 )
        signal.sigma.release() 
        result, frame = model. fitTo ( dataset0 )
    
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print result
    else :
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean                 is: %-28s ' %  result ( signal.mean  )[0] )
        logger.info ( 'Sigma                is: %-28s ' %  result ( signal.sigma )[0] )
        logger.info ( 'Gamma                is: %-28s ' %  result ( signal.gamma )[0] )

    models.add ( model )


# =============================================================================
## PseudoVoigt
# =============================================================================
def test_pvoigt () :
    
    logger.info ('Test PSeudoVoigt_pdf: fast approximation to Voigt profile' )
    model = Models.Fit1D (
        signal = Models.PseudoVoigt_pdf ( 'PV' , 
                                          xvar  = mass                ,
                                          mean  = signal_gauss.mean   , 
                                          sigma = signal_gauss.sigma  ) , 
        background = Models.Bkg_pdf ('BkgV', xvar = mass , power = 0 )) 
    
    signal = model.signal
    signal.mean .fix() 
    signal.sigma.fix ( m.error() )
    signal.gamma.fix ( 0.002     )
    
    model.B.setVal( 500)
    model.S.setVal(5000)
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 )
        result, frame = model. fitTo ( dataset0 )
        model.signal.sigma.release() 
        result, frame = model. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print result
    else :
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean                 is: %-28s ' %  result ( signal.mean  )[0] )
        logger.info ( 'Sigma                is: %-28s ' %  result ( signal.sigma )[0] )
        logger.info ( 'Gamma                is: %-28s ' %  result ( signal.gamma )[0] )

    models.add ( model )


# =============================================================================
## Breit-Wigner
# =============================================================================
def test_bw () :
    
    logger.info ('Test BreitWigner_pdf' )
    
    ff = cpp.Ostap.Math.FormFactors.BlattWeisskopf( 1 , 3.5 ) ## formfactor 
    bw = cpp.Ostap.Math.BreitWigner (
    m.value() ,
    m.error() ,
    0.150     , ## m1 
    0.150     , ## m2 
    1         , ## orbital momentum
    ff          ## formfactor 
    )
    
    model = Models.Fit1D (
        signal = Models.BreitWigner_pdf
        ( name        = 'BW'              ,
          breitwigner = bw                ,     
          xvar        = mass              ,
          mean        = signal_gauss.mean ) ,
        background = Models.Bkg_pdf ('BkgBW', xvar = mass , power = 0 )) 

    signal = model.signal 
    model.S.setVal(5000)
    model.B.setVal(500)
    
    signal.mean.fix ( m.value() )
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 )
        signal.mean .release()
        signal.gamma.release()
        result, frame = model. fitTo ( dataset0 )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print result
    else :
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean                 is: %-28s ' %  result ( signal.mean  )[0] )
        logger.info ( 'Gamma                is: %-28s ' %  result ( signal.gamma )[0] )

    models.add ( model )



# =============================================================================
## Slash
# =============================================================================
def test_slash(): 
    logger.info ('Test Slash shape' ) 
    model = Models.Fit1D (
        signal = Models.Slash_pdf ( name  = 'Slash' , 
                                    xvar  = mass   ,
                                    mean  = signal_gauss.mean   , 
                                    scale = signal_gauss.sigma  ) ,
        background = None  )   
    
    signal = model.signal
    signal.scale.release() 
    signal.mean.fix()
    
    model.S.setVal(5000)
    model.B.setVal( 500)
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 )
        result, frame = model. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :     
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean                 is: %-28s ' %  result ( signal.mean  )[0] )
        logger.info ( 'Scale                is: %-28s ' %  result ( signal.scale )[0] )
        
    models.add ( model )


# =============================================================================
## Asymmetric Laplace 
# =============================================================================
def test_laplace(): 
    logger.info ('Test Asymmetric Laplace shape' ) 
    model = Models.Fit1D (
        signal = Models.AsymmetricLaplace_pdf ( name  = 'AL' , 
                                                xvar  = mass   ,
                                                mean  = signal_gauss.mean   , 
                                                slope = signal_gauss.sigma  ) ,
        background = None  )
    
    signal = model.signal
    signal.slope.release() 
    signal.mean.fix()
    
    model.S.setVal(5000)
    model.B.setVal( 500)
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 )
        result, frame = model. fitTo ( dataset0 )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print result 
    else :     
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean                 is: %-28s ' %  result ( signal.mean  )[0] )
        logger.info ( 'Slope                is: %-28s ' %  result ( signal.slope )[0] )
        logger.info ( 'Asymmetry            is: %-28s ' %  result ( signal.asym  )[0] )
        
    models.add ( model )

# =============================================================================
## Test Rasing cosine 
# =============================================================================
def test_rasingcosine () :
    
    logger.info("Test RaisingCosine")
    model = Models.Fit1D (
        signal = Models.RaisingCosine_pdf( 'RC'                     ,
                                           xvar = mass              , 
                                           mean = signal_gauss.mean ) ,
        background = 1 ) 
    
    signal = model.signal
    model.S.setVal(5000)
    model.B.setVal( 500)
    
    with rooSilent() : 
        result,f  = model.fitTo ( dataset0 )  
        result,f  = model.fitTo ( dataset0 )  
        signal.mean  .release()
        signal.scale .release()
        result,f  = model.fitTo ( dataset0 )  
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )
        print result
    else :
        logger.info ( 'Signal & Background are: %-28s & %-28s ' % ( result ( 'S'         )[0] , result( 'B'           )[0] ) ) 
        logger.info ( 'Mean                 is: %-28s ' %  result ( signal.mean  )[0] )
        logger.info ( 'Scale                is: %-28s ' %  result ( signal.scale )[0] )

    models.add ( model )


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

    test_gauss          () ## simple Gaussian PDF                       + background 
    test_crystalball    () ## Crystal Ball                              + background
    test_crystalball_RS () ## right-side Crystal Ball                   + background  
    test_crystalball_DS () ## double side Crystal Ball                  + background 
    test_needham        () ## Needham function (CB with alpha=f(sigma)) + background 
    test_apolonios      () ## Apolonios function                        + background 
    test_apolonios2     () ## modified Apolonios function               + background 
    test_bifurcated     () ## bifurcated Gaussian function              + background 
    test_2gauss         () ## double     Gaussian function              + background 
    test_gengauss_v1    () ## generalized Gaussian function V1          + background 
    test_gengauss_v2    () ## generalized Gaussian function V2          + background 
    test_skewgauss      () ## skew gaussian                             + background 
    test_qgauss         () ## q-Gaussian function                       + background 
    test_bukin          () ## Bukin - skew Gaussian core with exponential tails  + background 
    test_studentT       () ## Student-t shape                           + background 
    test_bifstudentT    () ## Bifurcated Student-t shape                + background 
    test_sinhasinh      () ## Sinh-Asinh distribution                   + background 
    test_johnsonSU      () ## Johnson-SU distribution                   + background 
    test_atlas          () ## Modified Gaussian used by ATLAS/Zeus      + background 
    test_sech           () ## Sech (1/cosh)  distribution               + background 
    test_logistic       () ## Logistic distribution                     + background 
    test_voigt          () ## Voigt profile                             + background 
    test_pvoigt         () ## Pseudo-Voigt(approximation to Voigt)      + background 
    test_bw             () ## Breit-Wigner(+resolution)                 + background 
    test_slash          () ## Slash-function                            + background 
    test_rasingcosine   () ## Raising Cosine                            + background 
    test_laplace        () ## Laplace-function                            + background 
    
    ## check finally that everything is serializeable:
    test_db ()          

# =============================================================================
# The END 
# ============================================================================= 
