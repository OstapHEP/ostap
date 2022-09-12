#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_models.py
# Test module for ostap/fitting/models.py
# - It tests various ``signal-like''/``peak-like'' shapes
# ============================================================================= 
""" Test module for ostap/fitting/models.py
- It tests various ``signal-like''/``peak-like'' shapes 
"""
# ============================================================================= 
from   __future__               import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import cpp, VE, dsID, rooSilent 
from   ostap.utils.timing       import timing
from   builtins                 import range
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait
from   ostap.fitting.background import make_bkg
from   ostap.logger.colorized   import attention
import ostap.logger.table       as     T
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_models' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

## make simple test mass 
mass     = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 2.9 , 3.3 )

## book very simple data set
varset0  = ROOT.RooArgSet  ( mass )

dataset0 = ROOT.RooDataSet ( dsID() , 'Test Data set-0' , varset0 )  

mmin , mmax = mass.minmax()

NS = 10000
NB =  1000

## fill it 
m  = VE(3.100,  0.015**2)
m1 = VE(3.100,2*0.015**2)
NS1 = int ( 0.75 * NS )
NS2 = NS  - NS1 
for i in range(0,NS1) :
    mass.value = m.gauss () 
    dataset0.add ( varset0 )
for i in range(0,NS2) :
    mass.value = m1.gauss () 
    dataset0.add ( varset0 )

for i in range(0,NB) :
    mass.value = random.uniform ( mmin , mmax ) 
    dataset0.add ( varset0   )

logger.info ('DATASET\n%s' % dataset0 )

models = set () 

## signal component 
signal_gauss = Models.Gauss_pdf ( name  = 'Gauss'    , ## the name 
                                  xvar  = mass       , ## the variable 
                                  mean  = m.value () , ## mean value (fixed)
                                  sigma = m.error () ) ## sigma      (fixed)

background = make_bkg ( 0 , 'Bkg' , xvar = mass , logger = logger ) 
## construct composite model: signal + background 
model_gauss = Models.Fit1D(
    signal     = signal_gauss ,
    background = background   
    )

signal_gauss.sigma.setMin (  0.1 * m.error () )
signal_gauss.sigma.setMax ( 10.0 * m.error () )

S = model_gauss.S
B = model_gauss.B


stats   = {}
results = [] 
# =============================================================================
def make_print ( pdf , fitresult , title , logger = logger ) :

    table  = fitresult.table ( title  = title , prefix = '# ' ) 
    
    if 0 != fitresult.status() or 3 != fitresult.covQual() :
        logger.warning ('%s: fit result\n%s' % ( title , table ) )
    else :
        logger.info    ('%s: fit result\n%s' % ( title , table ) )

    signal = pdf.signal 

    mean, mode, median, midpoint, rms , fwhm, skewness , kurtosis = 8 * ( attention ( '<error>' ) , )
    
    model = signal
    
    try    :
        mean     =  "%+.3g" % model.get_mean     () 
    except :
        pass
    
    try    :
        mode     =  "%+.3g" % model .mode       () 
    except :
        pass 

    try    :
        median   =  "%+.3g" % model .median     () 
    except :
        pass
    
    try    :
        midpoint =  "%+.3g" % model .mid_point  () 
    except :
        pass

    try    :
        rms      =  "%+.3g" % model .rms        () 
    except :
        pass
    
    try    :
        fwhm     =  "%+.3g" % model .fwhm       () 
    except :
        pass
    
    try    :
        skewness =  "%+.3g" % model .skewness   () 
    except :
        pass

    try    :
        kurtosis =  "%+.3g" % model .kurtosis   () 
    except :
        pass
    
    row   = model.name  , mean , mode , median , midpoint , rms , fwhm , skewness , kurtosis

    stats [ model.name ] = row

    with wait ( 1 ), use_canvas ( title ) : 
        pdf.draw (  dataset0 )
    
# =============================================================================
## gauss PDF
# =============================================================================
def test_gauss() :

    logger = getLogger ( 'test_gauss' )
    
    logger.info ('Test Gauss_pdf:  simple Gaussian signal' )

    model  = model_gauss
    signal = model.signal
    
    ## release the sigma of signal:
    signal.mean.release()
    
    ## simple fit with gaussian only 
    result = signal. fitTo ( dataset0 , silent = True )
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , 'Simple Gaussian model' , logger )

    models.add ( model )
    results.append ( result  )

# =============================================================================
## CrystalBall PDF
# =============================================================================
def test_crystalball () :
    
    logger = getLogger ( 'test_crystalball' )
    logger.info ('Test CrystalBall_pdf: Crystal Ball  function' )
    
    ## composite model: signal + background 
    model = Models.Fit1D (
        signal     = Models.CrystalBall_pdf ( name  = 'CB'     , ## the name 
                                              xvar  = mass     , ## the variable   
                                              alpha = (2,1,5)  , ## tail parameter
                                              n     = (3,1,9)  , ## tail parameter 
                                              sigma = signal_gauss.sigma ,   ## reuse sigma from gauss
                                              mean  = signal_gauss.mean  ) , ## reuse mean  from gauss 
        background = background   ,
        S = S , B = B 
        )
    
    model.signal.n.fix(8)
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        model.signal.alpha.release()
        model.signal.n    .release()
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , 'Crystal Ball model' , logger )
                      
    models.add ( model )
    results.append ( result  )

# =============================================================================
## right side CrystalBall PDF
# =============================================================================
def test_crystalball_RS () :
    
    logger = getLogger ( 'test_crystalball_RS' )

    logger.info ('Test CrystalBallRS_pdf:  right-side Crystal Ball function' )
    model = Models.Fit1D (
        signal = Models.CrystalBallRS_pdf ( name  = 'CBRS' , 
                                            xvar  = mass               ,
                                            sigma = signal_gauss.sigma ,
                                            alpha = (1.5, 0.5 , 3.0)   ,
                                            n     = (5,1,10)           , 
                                            mean  = signal_gauss.mean  ) ,
        background = background   ,
        S = S , B = B 
        )
    
    model.S = NS 
    model.B = NB
    
    model.signal.n.fix(8)

    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        model.signal.alpha.release()
        model.signal.n    .release()
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , '(Right-side) Crystal Ball model' , logger )
    
    models.add ( model )
    results.append ( result  )

# =============================================================================
## double sided CrystalBall PDF
# =============================================================================
def test_crystalball_DS () :
    
    logger = getLogger ( 'test_crystalball_DS' )
    logger.info ('Test CrystalBallDS_pdf: double-sided Crystal Ball function' )
    model = Models.Fit1D (
        signal = Models.CB2_pdf ( name   = 'CB2'              , 
                                  xvar   = mass               ,
                                  nL     = 10                 , 
                                  nR     = 10                 , 
                                  alphaL = (1.5,0.5,3)        , 
                                  alphaR = (1.4,0.5,3)        , 
                                  sigma  = signal_gauss.sigma ,  
                                  mean   = signal_gauss.mean  ) ,
        background = background   ,
        S = S , B = B 
        )
    
    model.S = NS 
    model.B = NB
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 )
        model.signal.aL.release()
        model.signal.aR.release()
        result, frame = model. fitTo ( dataset0 )
        model.signal.aL.fix(1.5) 
        model.signal.aR.fix(1.5)    
        result, frame = model. fitTo ( dataset0 )
        
    make_print ( model, result , 'Double-sided Crystal Ball model' , logger )

    models.add ( model )
    results.append ( result  )

# =============================================================================
## Needham PDF
# =============================================================================
def test_needham() :

    logger = getLogger ( 'test_needham' )
  
    logger.info ('Test Needham_pdf: Crystal Ball with alpha=f(sigma)' )
    model = Models.Fit1D (
        signal = Models.Needham_pdf ( name  = 'Matt'             , 
                                      xvar  = mass               ,
                                      sigma = signal_gauss.sigma ,  
                                      mean  = signal_gauss.mean  ) ,
        background = background   ,
        S = S , B = B 
        )
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        model.signal.mean .release()
        model.signal.sigma.release()
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )

    make_print ( model , result , 'Needham model' , logger )

    models.add ( model )
    results.append ( result  )

# ==========================================================================
## Apollonios
# ==========================================================================
def test_apollonios () :

    logger = getLogger ( 'test_apollonios' )

    logger.info ('Test Apollonios_pdf: Modified gaussian with power-law and exponential tails' ) 
    model = Models.Fit1D (
        signal = Models.Apollonios_pdf ( name  = 'APO', 
                                         xvar  = mass ,
                                         mean  = signal_gauss.mean    ,
                                         sigma = signal_gauss.sigma ,
                                         b     =  1 ,
                                         n     = 10 ,
                                         alpha =  3 ) ,
        background = background   ,
        S = S , B = B 
        )
    
    model.S = NS 
    model.B = NB 
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , 'Apollonios model' , logger )

    models.add ( model )
    results.append ( result  )

# ==========================================================================
## Apollonios2
# ==========================================================================
def test_apollonios2() :
    
    logger = getLogger ( 'test_apollonios2' )

    logger.info ('Test Apollonios2_pdf: modified Gaussian with exponential tails' ) 
    model = Models.Fit1D (
        signal = Models.Apollonios2_pdf ( name = 'AP2' , 
                                          xvar      = mass ,
                                          mean      = signal_gauss.mean  ,
                                         sigma     = signal_gauss.sigma ,
                                         beta      =  ( 0.5 , 2 )       ,
                                         asymmetry = 0 ) ,
        background = background   ,
        S = S , B = B 
        )
    
    model.signal.mean.fix( m.value() )    
    model.S = NS 
    model.B = NB
    model.signal.sigma.release() 
    
    with rooSilent() :
        result, frame = model. fitTo ( dataset0 , silent = True )
        model.signal.asym.release ()
        model.signal.mean.release ()
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model, result , 'Apollonios2 model' , logger )
        
    models.add ( model )
    results.append ( result  )

# =============================================================================
## Bifurcated gauss PDF
# =============================================================================
def test_bifurcated () :
    
    logger = getLogger ( 'test_bifurcated' )
    logger.info ('Test BifurcatedGauss_pdf: Bifurcated Gaussian' )
    
    signal = Models.BifurcatedGauss_pdf ( name = 'BfGau' ,
                                          mean  = signal_gauss.mean  ,
                                          sigma = signal_gauss.sigma ,
                                          xvar  = mass    )
    signal.asym = 0 
    
    model = Models.Fit1D(
        signal     = signal     ,
        background = background ,
        S = S , B = B 
        ) 
    
    model.B = NB 
    model.S = NS
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , 'Bifurcated Gaussian model' , logger )

    models.add ( model )
    results.append ( result  )

# =============================================================================
## Double gauss PDF
# =============================================================================
def test_2gauss () :
    
    logger = getLogger ( 'test_2gauss' )
    logger.info ('Test DoubleGauss_pdf: Double Gaussian' )
    
    signal = Models.DoubleGauss_pdf ( name = 'Gau2' ,
                                             mean  = signal_gauss.mean  ,
                                             sigma = signal_gauss.sigma ,
                                             xvar  = mass               ,
                                             fraction = 0.9             ,
                                             scale    = 1.2             )
    
    model = Models.Fit1D(
        signal     = signal      ,
        background = background  ,
        S = S , B = B 
        )
    
    model.B = NB 
    model.S =NS
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.fraction.release() 
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , 'Double Gaussian model' , logger )
        
    models.add ( model )
    results.append ( result  )

# =============================================================================
## GenGaussV1
# =============================================================================
def test_gengauss_v1 () :
    
    logger = getLogger ( 'test_gengauss_v1' )

    logger.info ('Test GenGaussV1_pdf: Generalized Gaussian V1' ) 
    model = Models.Fit1D (
        signal = Models.GenGaussV1_pdf ( name = 'Gv1' , 
                                         xvar = mass  ,
                                         mean = signal_gauss.mean ) ,
        background = background   ,
        S = S , B = B 
        ) 
    
    model.signal.beta .fix(2)
    model.signal.mean .fix( m.value() ) 
    model.S = NS 
    model.B = NB 
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        model.signal.alpha.release()
        result, frame = model. fitTo ( dataset0 , silent = True )
        model.signal.beta .release() 
        model.signal.mean .release() 
        result, frame = model . fitTo ( dataset0 , silent = True )
        result, frame = model . fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , 'Generalized Gaussian V1 model' , logger )

    models.add ( model )
    results.append ( result  )

# =============================================================================
## GenGaussV2
# =============================================================================
def test_gengauss_v2 () : 

    logger = getLogger ( 'test_gengauss_v2' )
    
    logger.info ('Test GenGaussV2_pdf: Generalized Gaussian function V2' ) 
    model = Models.Fit1D (
        signal = Models.GenGaussV2_pdf ( name = 'Gv2' , 
                                         xvar = mass  ,
                                         mean = signal_gauss.mean ) ,
        background = background   ,
        S = S , B = B 
        ) 
    
    model.signal.kappa.fix(0)
    
    model.S = NS 
    model.B = NB 
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        model.signal.kappa.release() 
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )        

    make_print ( model , result , 'Generalized Gaussian V2 model' , logger )
    
    models.add ( model )
    results.append ( result  )

# =============================================================================
## SkewGauss
# =============================================================================
def test_skewgauss() :
    
    logger = getLogger ( 'test_skewgauss' )
    
    logger.info ('Test SkewGauss_pdf: Skew Gaussian function' ) 
    model = Models.Fit1D (
        signal = Models.SkewGauss_pdf ( name = 'GSk' , 
                                        xvar = mass  , mean = signal_gauss.mean ) ,
        background = background   ,
        S = S , B = B 
        ) 

    signal = model.signal 
    model.S = NS 
    model.B = NB 
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.mean .release() 
        signal.sigma.release() 
        result, frame = model. fitTo ( dataset0 , silent = True )

    make_print ( model , result , 'Skew  Gaussian model' , logger )

    models.add ( model )
    results.append ( result  )

# =============================================================================
## qGauss
# =============================================================================
def test_qgauss () :
    
    logger = getLogger ( 'test_qgauss' )

    logger.info ('Test QGaussian_pdf: q-Gaussian' ) 
    model = Models.Fit1D (
        signal = Models.QGaussian_pdf ( name = 'qG'  , 
                                        xvar = mass  ,
                                        q    = (1,0.7,1.2), 
                                        mean = signal_gauss.mean   ,
                                        scale = signal_gauss.sigma ) ,
        background = background   ,
        S = S , B = B 
        ) 

    signal = model.signal
    model.S = NS 
    model.B = NB 

    signal.scale = 0.015
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.scale.release()
        signal.mean.release()
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.q .release() 
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , 'q-Gaussian model' , logger )

    models.add ( model )
    results.append ( result  )


# =============================================================================
## kGauss
# =============================================================================
def test_kgauss () :
    
    logger = getLogger ( 'test_kgauss' )

    logger.info ('Test KGaussian_pdf: k-Gaussian' ) 
    model = Models.Fit1D (
        signal = Models.KGaussian_pdf ( name  = 'kG'  , 
                                        xvar  = mass  ,
                                        kappa = (0, -10, 10 ), 
                                        mean  = signal_gauss.mean   ,
                                        scale = signal_gauss.sigma ) ,
        background = background   ,
        S = S , B = B 
        ) 

    signal = model.signal
    model.S = NS 
    model.B = NB 

    signal.scale = 0.015
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.scale.release()
        signal.mean.release()
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.kappa.release() 
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , 'k-Gaussian model' , logger )

    models.add ( model )
    results.append ( result  )

# =============================================================================
## Novosibirs
# =============================================================================
def test_novosibirsk () :
    
    logger = getLogger ( 'test_novosibirsk' )

    logger.info ('Test Bukin_pdf: Novisibirsk function: asymmetric tail' ) 
    model = Models.Fit1D (
        signal = Models.Novosibirsk_pdf ( name  = 'Novosibirsk' ,
                                          xvar  = mass    ,
                                          tau   = ( 0 , -1 , 1 ) ,
                                          mean  = signal_gauss.mean  , 
                                          sigma = signal_gauss.sigma ) ,
        background = background   ,
        S = S , B = B 
        )
    
    model.signal.mean .fix  ( m.value() )
    model.signal.sigma.fix  ( m.error() )
    model.S = NS 
    model.B = NB 
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        model.signal.tau .release()     
        model.signal.mean .release() 
        model.signal.sigma.release() 
        result , frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , 'Novisibirsk model' , logger )        

    models.add ( model )
    results.append ( result  )

# =============================================================================
## Bukin
# =============================================================================
def test_bukin() :
    
    logger = getLogger ( 'test_bukin' )

    logger.info ('Test Bukin_pdf: Bukin function: skew gaussian core + exponenial/gaussian  tails' ) 
    model = Models.Fit1D (
        signal = Models.Bukin_pdf ( name  = 'Bukin' ,
                                    xvar  = mass    ,
                                    xi    = 0    ,
                                    rhoL  = 0    ,
                                    rhoR  = 0    , 
                                    mean  = signal_gauss.mean  , 
                                    sigma = signal_gauss.sigma ) ,
        background = background   ,
        S = S , B = B 
        )
    
    model.signal.mean .fix  ( m.value() )
    model.signal.sigma.fix  ( m.error() )
    model.S = NS 
    model.B = NB 
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        model.signal.xi  .release()     
        result, frame = model. fitTo ( dataset0 , silent = True )
        model.signal.rhoL.release()     
        result, frame = model. fitTo ( dataset0 , silent = True )
        model.signal.rhoR.release()     
        result, frame = model. fitTo ( dataset0 , silent = True )
        model.signal.mean .release() 
        model.signal.sigma.release() 
        result, frame = model. fitTo ( dataset0 , silent = True )

    make_print ( model , result , 'Bukin (modified Novosibirsk) model' , logger )        

    models.add ( model )
    results.append ( result  )

# =============================================================================
## StudentT
# =============================================================================
def test_studentT () :

    logger = getLogger ( 'test_studentT' )
       
    logger.info ('Test StudentT_pdf: Student-t distribution' ) 
    model = Models.Fit1D (
        signal = Models.StudentT_pdf ( name = 'ST' , 
                                       xvar = mass ,
                                       mean = signal_gauss.mean ) ,
        background = background   ,
        S = S , B = B 
        ) 
    
    model.signal.n      = 20 
    model.signal.sigma = 0.015 

    model.S = NS 
    model.B = NB 

    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Student's t-distribution" , logger )
    
    models.add ( model )
    results.append ( result  )

# =============================================================================
## Bifurcated StudentT
# =============================================================================
def test_bifstudentT():

    logger = getLogger ( 'test_bifstudentT' )

    logger.info ('Test bifurcated StudentT_pdf: bifurcated Student-t' ) 
    model = Models.Fit1D (
        signal = Models.BifurcatedStudentT_pdf ( name  = 'BfST' , 
                                                 xvar  = mass   ,
                                                 nL    = 25     ,
                                                 nR    = 25     ,                                                 
                                                 mean  = signal_gauss.mean   , 
                                                 sigma = signal_gauss.sigma  ) ,
        background = background   ,
        S = S , B = B 
        ) 
    
    signal = model.signal 
    model.S = NS 
    model.B = NB
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.nL   .release()
        signal.nR   .release()
        signal.sigma.release()
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Bifurcated Student's t-distribution" , logger )        

    models.add ( model )
    results.append ( result  )

# ==========================================================================
## PearsonIV 
# ==========================================================================
def test_PearsonIV () :
    
    logger = getLogger ( 'test_PearsonIV' )
       
    
    logger.info ('Test PearsonIV_pdf: asymmetric pdf ' ) 
    model = Models.Fit1D (
        signal = Models.PearsonIV_pdf ( name = 'PIV' , 
                                        xvar      = mass                ,
                                        mu        = signal_gauss.mean   ,
                                        varsigma  = signal_gauss.sigma  ,
                                        n         = ( 10 ,  1   , 500 ) , 
                                        kappa     = ( 0  , -100 , 100 ) ) ,
        background = background   ,
        S = S , B = B ,
        )


    signal = model.signal 
    model.S = NS 
    model.B = NB

    model.signal.kappa.fix (   0 )
    model.signal.n.fix     ()
    signal.mu    .fix      ( 3.1 )  
    with rooSilent() :
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.sigma .release ()
        signal.n     .release () 
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.kappa .release ()
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Pearson Type IV distribution" , logger )        

    models.add ( model )
    results.append ( result  )

# =============================================================================
## Test  SinhAsinh-Distribution
# =============================================================================
def test_sinhasinh() :

    logger = getLogger ( 'test_sinhasinh' ) 

    logger.info("Test  SinhAsinh-Distribution")
    model = Models.Fit1D (
        signal = Models.SinhAsinh_pdf( 'SASH'                   ,
                                       xvar = mass              , 
                                       mean = signal_gauss.mean ) ,
        background = background   ,
        S = S , B = B 
        )
    
    signal = model.signal
    
    model.S = NS 
    model.B = NB

    
    signal.mu      = 3.10  
    signal.sigma   = 0.015  
    signal.epsilon = 0.021 
    signal.delta   = 1.0   

    with rooSilent() : 
        result,f  = model.fitTo ( dataset0 , silent = True )  
        result,f  = model.fitTo ( dataset0 , silent = True )  
        signal.delta.release()
        result,f  = model.fitTo ( dataset0 , silent = True )  
        signal.epsilon.release()
        result,f  = model.fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Sinh-asinh model" , logger )        

    models.add ( model )
    results.append ( result  )

# =============================================================================
## Test  JohnsonSU-Distribution
# =============================================================================
def test_johnsonSU () :

    logger = getLogger ( 'test_johnsonSU' ) 
     
    logger.info("Test  JohnsonSU-Distribution")
    model = Models.Fit1D (
        signal = Models.JohnsonSU_pdf( 'JSU'                    ,
                                       xvar = mass              , 
                                       xi   = signal_gauss.mean ) ,
        background = background   ,
        S = S , B = B 
        )
    
    signal = model.signal
    
    model.S = NS 
    model.B = NB
    
    with rooSilent() : 
        result,f  = model.fitTo ( dataset0 , silent = True )  
        signal.lambd .release()
        signal.delta.release()
        result,f  = model.fitTo ( dataset0 , silent = True )  
        signal.gamma.release()
        signal.mean .release()
        result,f  = model.fitTo ( dataset0 , silent = True )
        
    make_print ( model, result , "Johnson's SU model" , logger )        

    models.add ( model )
    results.append ( result  )

# =============================================================================
## Test  ATLAS
# =============================================================================
def test_atlas () :
    
    logger = getLogger ( 'test_atlas' ) 
        

    logger.info("Test  ATLAS: Modified Gaussian, used by ATLAS/Zeus")
    model = Models.Fit1D (
        signal = Models.Atlas_pdf( 'ATLAS'                  ,
                                   xvar = mass              , 
                                   mean = signal_gauss.mean ) ,
        background = background   ,
        S = S , B = B 
        )
    
    signal = model.signal

    model.S = NS 
    model.B = NB
    
    with rooSilent() : 
        result,f  = model.fitTo ( dataset0 , silent = True )  
        result,f  = model.fitTo ( dataset0 , silent = True )  
        signal.mean  .release()
        signal.sigma .release()
        result,f  = model.fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "ATLAS/ZEUS model" , logger )        

    models.add ( model )
    results.append ( result  )
       
# ==========================================================================
## Das
# ==========================================================================
def test_das_1 () :
    
    logger = getLogger ( 'test_das_1' )
       
    
    logger.info ('Test Das_pdf: Das pdf with two tails' ) 
    model = Models.Fit1D (
        signal = Models.Das_pdf ( name = 'Das1' , 
                                  xvar      = mass               ,
                                  mu        = signal_gauss.mean  ,
                                  sigma     = signal_gauss.sigma ,
                                  kL        = ( 1 , 0.1 , 10.0 ) ,
                                  kR        = ( 1 , 0.1 , 10.0 ) ) ,
        background = background   ,
        S = S , B = B ,
        )

    signal = model.signal 
    signal.mean .fix ( m.value() )
    signal.sigma.fix ( m.error() )

    model.S = NS 
    model.B = NB
    
    with rooSilent() :
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.mean .release() 
        signal.sigma.release() 
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Das model (1) " , logger )        

    models.add ( model )
    results.append ( result  )

# ==========================================================================
## Das
# ==========================================================================
def test_das_2 () :
    
    logger = getLogger ( 'test_das_2' )
       
    logger.info ('Test Das_pdf: Das pdf with tail and asymmetry' ) 
    model = Models.Fit1D (
        signal = Models.Das_pdf ( name = 'Das2' , 
                                  xvar      = mass               ,
                                  mu        = signal_gauss.mean  ,
                                  sigma     = signal_gauss.sigma ,
                                  k         = ( 1 , 0.1 , 10.0 ) ,
                                  kappa     = ( 0  , -1 , 1    ) ) ,
        background = background   ,
        S = S , B = B ,
        )

    signal = model.signal 
    signal.mean .fix ( m.value() )
    signal.sigma.fix ( m.error() )

    model.S = NS 
    model.B = NB
    
    with rooSilent() :
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.mean .release() 
        signal.sigma.release() 
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Das model (2) " , logger )        

    models.add ( model )
    results.append ( result  )

# ==========================================================================
## Hat
# ==========================================================================
def test_hat  () :
    
    logger = getLogger ( 'test_hat' )
       
    logger.info ('Test Hat_pdf: smooth simmetricfinite function' ) 
    model = Models.Fit1D (
        signal = Models.Hat_pdf ( name     = 'Hat'              , 
                                  xvar     = mass               ,
                                  mean     = signal_gauss.mean  ,
                                  varsigma = signal_gauss.sigma ) ,
        background = background   ,
        S = S , B = B ,
        )

    signal = model.signal 

    signal.sigma.release() 
    signal.mean .release() 

    model.S = NS 
    model.B = NB
    
    with rooSilent() :
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Hat model" , logger )        

    models.add ( model )

    signal.mean.release() 
    results.append ( result  )

# ==========================================================================
## Up
# ==========================================================================
def test_up  () :
    
    logger = getLogger ( 'test_up' )
       
    logger.info ('Test Up_pdf: smooth atomic finite function' ) 
    model = Models.Fit1D (
        signal = Models.Up_pdf ( name     = 'Up' , 
                                 xvar     = mass               ,
                                 mean     = signal_gauss.mean  ,
                                 varsigma = signal_gauss.sigma ) ,
        background = background   ,
        S = S , B = B ,
        )
    
    signal = model.signal 
    
    signal.mean .release() 
    
    signal.sigma = 0.057  
    signal.mean .fix ( m.value() ) 
    
    model.S = NS 
    model.B = NB
    
    with rooSilent() :
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Up model" , logger )        

    models.add ( model )

    signal.mean.release() 
    results.append ( result  )

# ==========================================================================
## FupN
# ==========================================================================
def test_fupn  () :
## if 1 < 2 :
    
    logger = getLogger ( 'test_fupN' )

    for N in range ( 8 , 9 ) :
        
        logger.info ('Test Fup%d_pdf: smooth atomic finite function' % N ) 
        model = Models.Fit1D (
            signal = Models.FupN_pdf ( name     = 'Fup%d' % N        , 
                                       xvar     = mass               ,
                                       N        = N                  ,
                                       mean     = signal_gauss.mean  ,
                                       varsigma = signal_gauss.sigma ) ,
            background = background   ,
            S = S , B = B ,
            )
        
        signal = model.signal 
        
        
        signal.sigma = 0.051  
        signal.mean .fix ( m.value() ) 
        
        model.S = NS 
        model.B = NB
        
        with rooSilent() :
            result, frame = model. fitTo ( dataset0 , silent = True )
            result, frame = model. fitTo ( dataset0 , silent = True )
            
        make_print ( model , result , "Fup%d model" % N  , logger )        
        
        models.add ( model )
        signal.mean.release() 
        results.append ( result  )

# =============================================================================
## Test  SECH
# =============================================================================
def test_sech() :
    
    logger = getLogger ( 'test_sech' )
    
    logger.info("Test  SECH:  Sech(1/cosh) distribution")
    model = Models.Fit1D (
        signal = Models.Sech_pdf( 'Sech'                   ,
                                  xvar = mass              , 
                                  mean = signal_gauss.mean ) ,
        background = background   ,
        S = S , B = B 
        )
    
    signal = model.signal
    model.S = NS 
    model.B = NB
    
    with rooSilent() : 
        result,f  = model.fitTo ( dataset0 , silent = True )  
        result,f  = model.fitTo ( dataset0 , silent = True )  
        signal.mean  .release()
        signal.sigma .release()
        result,f  = model.fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Sech model" , logger )        

    models.add ( model )
    results.append ( result  )

# =============================================================================
## Test  LOGISTIC
# =============================================================================
def test_logistic () :
    
    logger = getLogger ( 'test_logistic' )
 
    logger.info("Test  LOGISTIC: Logistic distribution")
    model = Models.Fit1D (
        signal = Models.Logistic_pdf( 'Logistic'                    ,
                                      xvar = mass              , 
                                      mean = signal_gauss.mean ) ,
        background = background   ,
        S = S , B = B 
        )
    
    signal = model.signal
    model.S = NS 
    model.B = NB
    
    with rooSilent() : 
        result,f  = model.fitTo ( dataset0 , silent = True )  
        result,f  = model.fitTo ( dataset0 , silent = True )  
        signal.mean  .release()
        signal.sigma .release()
        result,f  = model.fitTo ( dataset0 , silent = True )

    make_print ( model , result , "Logistic  model" , logger )        

    models.add ( model )
    results.append ( result  )

# =============================================================================
## Test  LOSEV
# =============================================================================
def test_losev() :
    
    logger = getLogger ( 'test_losev' )
    
    logger.info("Test  Losev: asymmetric hyperbilic secant distribution")
    model = Models.Fit1D (
        signal = Models.Losev_pdf( 'LOSEV'                  ,
                                   xvar = mass              , 
                                   mean = signal_gauss.mean ) ,
        background = background   ,
        S = S , B = B 
        ) 
    
    signal = model.signal
    model.S = NS 
    model.B = NB
    
    with rooSilent() : 
        result,f  = model.fitTo ( dataset0 , silent = True )  
        result,f  = model.fitTo ( dataset0 , silent = True )  
        signal.mean  .release()
        signal.alpha .release()
        signal.beta  .release()
        result,f  = model.fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Losev  model" , logger )        

    models.add ( model )
    results.append ( result  )

# =============================================================================
## Slash
# =============================================================================
def test_slash():

    logger = getLogger ( 'test_slash' )
    
    logger.info ('Test Slash shape' ) 
    model = Models.Fit1D (
        signal = Models.Slash_pdf ( name  = 'Slash' , 
                                    xvar  = mass   ,
                                    mean  = signal_gauss.mean   , 
                                    scale = signal_gauss.sigma  ) ,
        background = background   ,
        S = S , B = B 
        )
    
    signal = model.signal
    signal.scale.release() 
    signal.mean.fix()
    
    model.S = NS 
    model.B = NB

    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.mean.release ()
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Slash model" , logger )
    
    models.add ( model )
    results.append ( result  )


# =============================================================================
## Test Rasing cosine 
# =============================================================================
def test_raisngcosine () :
    
    logger = getLogger ( 'test_raisngcosine' )
    
    logger.info("Test RaisingCosine")
    model = Models.Fit1D (
        signal = Models.RaisingCosine_pdf( 'RCos'                     ,
                                           xvar = mass              , 
                                           mean = signal_gauss.mean ) ,
        background = background   ,
        S = S , B = B 
        )
    
    signal = model.signal
    model.S = NS 
    model.B = NB
    
    with rooSilent() : 
        result,f  = model.fitTo ( dataset0 , silent = True )  
        result,f  = model.fitTo ( dataset0 , silent = True )  
        signal.mean  .release()
        signal.scale .release()
        result,f  = model.fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Rising Cosine model" , logger )        

    models.add ( model )
    results.append ( result  )

# =============================================================================
## Asymmetric Laplace 
# =============================================================================
def test_laplace():
    
    logger = getLogger ( 'test_laplace' )

    logger.info ('Test Asymmetric Laplace shape' ) 
    model = Models.Fit1D (
        signal = Models.AsymmetricLaplace_pdf ( name  = 'Laplace' , 
                                                xvar  = mass   ,
                                                mean  = signal_gauss.mean   , 
                                                slope = signal_gauss.sigma  ) ,
        background = background   ,
        S = S , B = B 
        )
    
    signal = model.signal
    signal.slope.release() 
    signal.mean.fix()
    
    model.S = NS 
    model.B = NB
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.slope.release() 
        signal.mean .release()
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Laplace model" , logger )        

    models.add ( model )
    results.append ( result  )

# ==========================================================================
## ExGauss
# ==========================================================================
def test_exgauss () :
    
    logger = getLogger ( 'test_ExGauss' )
       
    
    logger.info ('Test ExGauss_pdf: single tail pdf ' ) 
    model = Models.Fit1D (
        signal = Models.ExGauss_pdf ( name = 'ExGauss' , 
                                      xvar      = mass               ,
                                      mu        = signal_gauss.mean  ,
                                      varsigma  = signal_gauss.sigma ,
                                      k         = ( 1 , -10 , 10.0 ) ) ,
        background = background   ,
        S = S , B = B ,
        )
    

    signal = model.signal
    model.S = NS 
    model.B = NB

    with rooSilent() :
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.mu    .release()
        signal.sigma.release()        
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "ExGauss model" , logger )        

    models.add ( model )
    results.append ( result  )

# ==========================================================================
## NormalLaplace 
# ==========================================================================
def test_normlapl () :
    
    logger = getLogger ( 'test_NormalLaplace' )
       
    
    logger.info ('Test NormalLaplace_pdf: two-sided tails pdf ' ) 
    model = Models.Fit1D (
        signal = Models.NormalLaplace_pdf ( name = 'NormalLaplace' , 
                                            xvar      = mass               ,
                                            mu        = signal_gauss.mean  ,
                                            varsigma  = signal_gauss.sigma ,
                                            kL        = ( 1 , -10 , 10.0 ) , 
                                            kR        = ( 1 , -10 , 10.0 ) ) ,
        background = background   ,
        S = S , B = B ,
        )
    
    signal_gauss.mean .fix ( m.value() )
    signal_gauss.sigma.fix ( m.error() )

    model.S = NS 
    model.B = NB
    
    with rooSilent() :
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal_gauss.mean .release() 
        signal_gauss.sigma.release() 
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Normal Laplace model" , logger )        

    models.add ( model )
    results.append ( result  )
        
# ==========================================================================
## Hyperbolic
# ==========================================================================
def test_hyperbolic() :

    logger = getLogger ( 'test_hyperbolic' )
       
    
    logger.info ('Test Hyperbolic_pdf: hyperbolic distribution (exponential tails)' ) 
    model = Models.Fit1D (
        signal = Models.Hyperbolic_pdf ( name = 'HB' , 
                                         xvar      = mass               ,
                                         mu        = signal_gauss.mean  ,
                                         sigma     = signal_gauss.sigma ,
                                         zeta      = ( 1   , 1.e-6 , 1e+6  ) ,
                                         kappa     = ( 0   , -1    ,   1 ) ) ,
        background = background   ,
        S = S , B = B ,
        )
    
    signal_gauss.mean .fix ( m.value() )
    signal_gauss.sigma.fix ( m.error() )

    signal  = model.signal 
    model.S = NS 
    model.B = NB
    
    with rooSilent() :
        signal.kappa.fix ( 0 )    
        result, frame = model . fitTo ( dataset0 , silent = True )
        signal.kappa.release ()
        signal.zeta .release ()
        signal.mean .release ()
        signal.sigma.release ()
        result, frame = model . fitTo ( dataset0 , silent = True )
        result, frame = model . fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Hyperbolic model" , logger )        

    models.add ( model )
    results.append ( result  )

# ==========================================================================
## Generalised Hyperbolic
# ==========================================================================
def test_genhyperbolic() :

    logger = getLogger ( 'test_genhyperbolic' )
       
    
    logger.info ('Test GenHyperbolic_pdf: generalised hyperbolic distribution (exponential tails)' ) 
    model = Models.Fit1D (
        signal = Models.GenHyperbolic_pdf ( name = 'GHB' , 
                                            xvar      = mass               ,
                                            mu        = signal_gauss.mean  ,
                                            sigma     = signal_gauss.sigma ,
                                            zeta      = ( 1   , 1.e-6 , 1e+6  ) ,
                                            lambd     = ( -2 , -10 , 10       ) , 
                                            kappa     = ( 0   , -1 ,   1 ) ) ,
        background = background   ,
        S = S , B = B ,
        )
    
    signal_gauss.mean .fix ( m.value() )
    signal_gauss.sigma.fix ( m.error() )

    signal  = model.signal 
    model.S = NS 
    model.B = NB

    from ostap.utils.gsl import gslCount
    with rooSilent() , gslCount() :
        signal.kappa.fix ( 0 )    
        result, frame = model . fitTo ( dataset0 , silent = True )
        signal.kappa.release ()
        signal.zeta .release ()
        signal.sigma.release ()
        signal.mean .release ()
        result, frame = model . fitTo ( dataset0 , silent = True )
        result, frame = model . fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Generalized Hyperbolic model" , logger )
    
    models.add ( model )
    results.append ( result  )

# ==========================================================================
## Hypatia 
# ==========================================================================
def test_hypatia () :
## if 1 < 2 :
    
    logger = getLogger ( 'test_hypatia' )
       
    
    logger.info ('Test Hypatia_pdf: Hypatia pdf' ) 
    model = Models.Fit1D (
        signal = Models.Hypatia_pdf ( name      = 'Hypatia' , 
                                      xvar      = mass               ,
                                      mu        = signal_gauss.mean  ,
                                      sigma     = signal_gauss.sigma ,
                                      zeta      = ( 1   , 1.e-6 , 1e+6 ) ,
                                      lambd     = ( -2 , -10 , 10      ) , 
                                      kappa     = ( 0   , -1 ,   1 ) ,
                                      sigma0    = 0.005     ) , 
        background = background   ,
        S = S , B = B ,
        )

    signal_gauss.mean .fix ( m.value() )
    signal_gauss.sigma.fix ( m.error() )

    signal  = model.signal 
    model.S = NS 
    model.B = NB

    from ostap.utils.gsl import gslCount
    with rooSilent() , gslCount() :
        signal.kappa.fix ( 0 )    
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.kappa.release ()
        signal.zeta .release ()
        signal.mean .release ()
        signal.sigma.release ()
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Hypatia model" , logger )        
    
    models.add ( model )
    results.append ( result  )

# =============================================================================
## Voigt
# =============================================================================
def test_voigt () :
        
    logger = getLogger ( 'test_voigt' )
    
    logger.info ('Test Voigt_pdf: Breit-Wigner convoluted with Gauss' )
    model = Models.Fit1D (
        signal = Models.Voigt_pdf ( 'Voigt' , 
                                    xvar  = mass                ,
                                    m0    = signal_gauss.mean   , 
                                    sigma = signal_gauss.sigma  ) , 
        background = background   ,
        S = S , B = B 
        )
    
    signal = model.signal
    
    signal.sigma.fix ( m.error() )
    signal.gamma.fix ( 0.002     )
    signal.m0   .fix() 
    
    model.S = NS 
    model.B = NB
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.m0   .release() 
        signal.sigma.release() 
        signal.gamma.release() 
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Voigt model" , logger )        
    
    models.add ( model )
    results.append ( result  )

# =============================================================================
## PseudoVoigt
# =============================================================================
def test_pvoigt () :
    
    logger = getLogger ( 'test_pvoigt' )
 
    logger.info ('Test PSeudoVoigt_pdf: fast approximation to Voigt profile' )
    model = Models.Fit1D (
        signal = Models.PseudoVoigt_pdf ( 'pVoigt' , 
                                          xvar  = mass                ,
                                          m0    = signal_gauss.mean   , 
                                          sigma = signal_gauss.sigma  ) , 
        background = background   ,
        S = S , B = B 
        )
    
    signal = model.signal
    
    model.S = NS 
    model.B = NB

    signal.m0   .fix () 
    signal.sigma.fix ( m.error() )
    signal.gamma.fix ( 0.002     )

    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.m0   .release() 
        signal.sigma.release() 
        result, frame = model. fitTo ( dataset0 , silent = True )
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "pseudo-Voigt model" , logger )        

    models.add ( model )
    results.append ( result  )

# =============================================================================
## Breit-Wigner
# =============================================================================
def test_bw () :
    
    logger = getLogger ( 'test_bw' )

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
          m0          = signal_gauss.mean ) ,
        background = background   ,
        S = S , B = B 
        )

    signal = model.signal 
    model.S = NS 
    model.B = NB
    
    signal.mean.fix ( m.value() )
    
    with rooSilent() : 
        result, frame = model. fitTo ( dataset0 , silent = True )
        signal.m0   .release()
        signal.gamma.release()
        result, frame = model. fitTo ( dataset0 , silent = True )
        
    make_print ( model , result , "Breit-Wigner model" , logger )        

    models.add ( model )
    results.append ( result  )

## # =============================================================================
## ## Breit-Wigner with interference
## # =============================================================================
## def test_bwi () :
## ## if 1 < 2 :     
##     logger.info ('Test BWI_pdf' )
    
##     ff = cpp.Ostap.Math.FormFactors.BlattWeisskopf( 1 , 3.5 ) ## formfactor 
##     bw = cpp.Ostap.Math.BreitWigner (
##     m.value() ,
##     m.error() ,
##     0.150     , ## m1 
##     0.150     , ## m2 
##     1         , ## orbital momentum
##     ff          ## formfactor 
##     )
    
##     model = Models.Fit1D (
##         signal = Models.BWI_pdf
##         ( name        = 'BWI'             ,
##           breitwigner = bw                ,     
##           xvar        = mass              ,
##           mean        = signal_gauss.mean ,
##           bkg         = -1                ) ,
##         suffix     = '_a' , 
##         background = Models.PolyPos_pdf ('BkgBWI', xvar = mass , power = 1 ) ,
##         )
    
##     signal = model.signal 
##     model.S.setVal ( NS )
##     model.B.setVal ( NB )
    
##     signal.mean.fix ( m.value() )
    
##     with rooSilent() : 
##         result, frame = model. fitTo ( dataset0 )
##         signal.mean .release()
##         signal.gamma.release()
##         result, frame = model. fitTo ( dataset0 )
##         model.draw (  dataset0 )

##     if 0 != result.status() or 3 != result.covQual() :
##         logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual()  ) )

##     logger.info ( "Breit-Wigner function\n%s" % result.table ( prefix = "# " ) )

##     models.add ( model )

# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        db['mass'] = mass 
        db['vars'] = varset0 
        db['dataset'  ] = dataset0
        for m in models :
            db['model:' + m.name ] = m
            db['roo_tot:%s' % m.name ] = m.pdf
            db['roo_sig:%s' % m.name ] = m.signal    .pdf
            db['roo_bkg:%s' % m.name ] = m.background.pdf
        db['models'   ] = models
        for r in results :
            db ['result:%s' % r.name ] = r
        db['results'   ] = results 
        db.ls() 

# ==============================================================================
## dump all models
# ==============================================================================
def dump_models () :
    
    header =  'Model'   , \
             'mean'     , 'mode' , 'midpoint' , 'median' , \
             'rms'      , 'fwhm' ,  \
             'skewness' , 'kurtosis'
    
    rows = [ header ] 
    for m in sorted ( stats  ) :
        rows.append ( stats [ m ]  )

    table = T.table ( rows , title = "Model's features" ,  prefix = '# ' )
    logger.info ( 'Features of models\n%s' % table )
  

# =============================================================================
if '__main__' == __name__ :

    ## simple Gaussian PDF                       + background
    with timing ('test_gauss'          , logger ) :
        test_gauss          () 
        
    ## Crystal Ball                              + background
    with timing ('test_crystalball'    , logger ) :
        test_crystalball    () 
        
    ## right-side Crystal Ball                   + background
    with timing ('test_crystalball_RS' , logger ) :
        test_crystalball_RS () 

    ## double side Crystal Ball                  + background
    with timing ('test_crystalball_DS' , logger ) :
        test_crystalball_DS () 

    ## Needham function (CB with alpha=f(sigma)) + background 
    with timing ('test_needham'        , logger ) :
        test_needham        ()

    ## Apollonios function                       + background 
    with timing ('test_apollonios'     , logger ) :
        test_apollonios     () 

    ## modified Apollonios function              + background
    with timing ('test_apolloniois2'   , logger ) :
        test_apollonios2    () 

    ## bifurcated Gaussian function              + background
    with timing ('test_bifurcated'     , logger ) :
        test_bifurcated     () 

    ## double     Gaussian function              + background
    with timing ('test_2gauss'         , logger ) :
        test_2gauss         () 

    ## generalized Gaussian function V1          + background
    with timing ('test_gengauss_v1'    , logger ) :
        test_gengauss_v1    () 
        
    ## generalized Gaussian function V2          + background
    with timing ('test_gengauss_v2'    , logger ) :
        test_gengauss_v2    () 
    
    ## skew gaussian                             + background
    with timing ('test_skewgauss'      , logger ) :
        test_skewgauss      () 
        
    ## q-Gaussian function                       + background
    with timing ('test_qgauss'         , logger ) :
        test_qgauss         () 
        
    ## k-Gaussian function                       + background
    with timing ('test_kgauss'         , logger ) :
        test_kgauss         () 

    ## Novisibirsk + background
    with timing ('test_novosibirsk'          , logger ) :
        test_novosibirsk     ()

    ## Bukin - skew Gaussian core with exponential tails  + background         
    with timing ('test_bukun'          , logger ) :
        test_bukin          ()

    ## Student-t shape                           + background 
    with timing ('test_studentT'       , logger ) :
        test_studentT       () 
    
    ## Bifurcated Student-t shape                + background
    with timing ('test_bifstudentT'    , logger ) :
        test_bifstudentT    ()
        
    ## PearsonIV                                      + background
    with timing ('test_PearsonIV'          , logger ) :
        test_PearsonIV () 

    ## Sinh-Asinh distribution                   + background
    with timing ('test_sinhasinh'      , logger ) :
        test_sinhasinh      () 
    
    ## Johnson-SU distribution                   + background 
    with timing ('test_johnsonSU'      , logger ) :
        test_johnsonSU      () 

    ## Modified Gaussian used by ATLAS/Zeus      + background 
    with timing ('test_atlas'          , logger ) :
        test_atlas          () 
        
    ## Das/1                                       + background 
    with timing ('test_das_1'             , logger ) :
        test_das_1            () 

    ## Das/2                                       + background 
    with timing ('test_das_2'             , logger ) :
        test_das_2            () 

    ## Hat                                       + background 
    with timing ('test_hat'             , logger ) :
        test_hat            () 

    ## Up                                       + background 
    with timing ('test_up'             , logger ) :
        test_up            () 

    ## FupN                                       + background 
    with timing ('test_fupN'             , logger ) :
        test_fupn          () 

    ## Sech (1/cosh)  distribution               + background    
    with timing ('test_sech'           , logger )  : 
        test_sech           ()

    ## Logistic distribution                     + background 
    with timing ('test_logistic'       , logger ) :
        test_logistic       () 

    ## Asymmetric hyperbilic secant distribution + background         
    with timing ('test_losev'          , logger ) :
        test_losev          () 

    ## Slash-function                            + background 
    with timing ('test_slash'          , logger ) :
        test_slash          () 

    ## Raising Cosine                            + background 
    with timing ('test_raisngcosine'   , logger ) :
        test_raisngcosine   () 

    ## Laplace-function                            + background 
    with timing ('test_laplace'        , logger ) :
        test_laplace        () 

    ## ExGauss                                       + background 
    with timing ('test_ExGauss'             , logger ) :
        test_exgauss            () 

    ## Normal Laplace                                       + background
    with timing ('test_NormalLaplas'        , logger ) :
        test_normlapl           () 

    ## Hyperbolic                                 + background 
    with timing ('test_hyperbolic'     , logger ) :
        test_hyperbolic        () 
        
    ## Generalised Hyperbolic                      + background 
    with timing ('test_genhyperbolic'     , logger ) :
        test_genhyperbolic     () 
        
    ## ## Hypatia                                     + background 
    ## with timing ('test_hypatia'           , logger ) :
    ##     test_hypatia           ()

    ## Voigt profile                             + background
    with timing ('test_voigt'          , logger ) :
        test_voigt          () 
    
    ## Pseudo-Voigt(approximation to Voigt)      + background
    with timing ('test_pvoigt'         , logger ) :
        test_pvoigt         () 
        
    ## Breit-Wigner(+resolution)                 + background 
    with timing ('test_bw'             , logger ) :
        test_bw             () 

    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()

    dump_models () 
         
# =============================================================================
##                                                                      The END 
# ============================================================================= 
