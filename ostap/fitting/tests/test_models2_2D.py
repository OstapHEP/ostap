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
from   ostap.core.core      import Ostap, std, VE, dsID
from   ostap.logger.utils   import rooSilent 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_models2_2D' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make simple test mass 
m_x     = ROOT.RooRealVar ( 'mass_x' , 'Some test mass(X)' , 3 , 3.2 )
m_y     = ROOT.RooRealVar ( 'mass_y' , 'Some test mass(Y)' , 3 , 3.2 )

## book very simple data set
varset  = ROOT.RooArgSet  ( m_x , m_y )
dataset = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset )  

m = VE(3.100,0.015**2)

N_ss =  5000
N_sb =  2500
N_bs =  2500
N_bb =  5000

random.seed(0)

## fill it : 5000 events  Gauss * Gauss 
for i in xrange(0,N_ss) : 
    m_x.setVal  ( m.gauss() )
    m_y.setVal  ( m.gauss() )
    dataset.add ( varset  )


## fill it : 2500 events  Gauss * const  
for i in xrange(0,N_sb) : 
    m_x.setVal  ( m.gauss() )
    m_y.setVal  ( random.uniform ( m_y.getMin() , m_y.getMax() ) ) 
    dataset.add ( varset  )

## fill it : 2500 events  const * Gauss
for i in xrange(0,N_bs) : 
    m_x.setVal  ( random.uniform ( m_x.getMin() , m_x.getMax() ) ) 
    m_y.setVal  ( m.gauss() )
    dataset.add ( varset  )

## fill it : 5000 events  const * const
for i in xrange(0,N_bb) :

    m_x.setVal ( random.uniform ( m_x.getMin() , m_x.getMax() ) ) 
    m_y.setVal ( random.uniform ( m_y.getMin() , m_y.getMax() ) ) 
    dataset.add ( varset  )


logger.info ( 'Dataset:%s ' % dataset ) 


models = set()
# =============================================================================


# =============================================================================
## gauss as signal, const as background 
# =============================================================================
def test_const () :
    
    logger.info ('Simplest (factorized) fit model:  ( Gauss + const ) x ( Gauss + const ) ' )
    model   = Models.Fit2D (
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y )
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    model.bkg2   .tau  .fix ( 0 )
    model.bkgA   .tau  .fix ( 0 )
    model.bkgB   .tau  .fix ( 0 )
    
    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )
        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' % result ( model.ss ) [0] )
        logger.info ('S1xB2 : %20s' % result ( model.sb ) [0] )
        logger.info ('B1xS2 : %20s' % result ( model.bs ) [0] )
        logger.info ('B1xB2 : %20s' % result ( model.bb ) [0] )

    models.add ( model ) 

# =============================================================================
## gauss as signal, second order polynomial as background 
# =============================================================================
def test_p2xp2 () : 
    logger.info ('Simple (factorized) fit model:  ( Gauss + P2 ) (x) ( Gauss + P2 ) ' )
    model   = Models.Fit2D (
        suffix   = '_2' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        power1   = 2 , 
        power2   = 2 ,
        powerA   = 2 ,
        powerB   = 2 , 
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    model.bkg2   .tau  .fix ( 0 )
    model.bkgA   .tau  .fix ( 0 )
    model.bkgB   .tau  .fix ( 0 )
    
    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )
        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' % result ( model.ss ) [0] )
        logger.info ('S1xB2 : %20s' % result ( model.sb ) [0] )
        logger.info ('B1xS2 : %20s' % result ( model.bs ) [0] )
        logger.info ('B1xB2 : %20s' % result ( model.bb ) [0] )

    models.add ( model ) 


# =============================================================================
## gauss as signal, 1st order polynomial as background + non-factorizeable BB
# =============================================================================
def test_p1xp1_BB () : 
    logger.info ('Simplest non-factorized fit model:  ( Gauss + P1 ) (x) ( Gauss + P1 ) + BB' )
    model   = Models.Fit2D (
        suffix   = '_3' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        power1   = 1 , 
        power2   = 1 ,
        bkg2D    = Models.PolyPos2D_pdf ( 'P2D' , m_x , m_y , nx = 2 , ny = 2 ) 
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    model.bkg2   .tau  .fix ( 0 )
    
    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )
        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' % result ( model.ss ) [0] )
        logger.info ('S1xB2 : %20s' % result ( model.sb ) [0] )
        logger.info ('B1xS2 : %20s' % result ( model.bs ) [0] )
        logger.info ('B1xB2 : %20s' % result ( model.bb ) [0] )

    models.add ( model ) 


# =============================================================================
## gauss as signal, 1st order polynomial as background 
# =============================================================================
def test_p1xp1_BBs () :
    
    logger.info ('Non-factorized symmetric background fit model:  ( Gauss + P1 ) (x) ( Gauss + P1 ) + BBsym' )
    model   = Models.Fit2D (
        suffix   = '_4' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        power1   = 1 , 
        power2   = 1 ,
        bkg2D    = Models.PolyPos2Dsym_pdf ( 'P2Ds' , m_x , m_y , n = 2 ) 
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    model.bkg2   .tau  .fix ( 0 )
    
    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )
        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' % result ( model.ss ) [0] )
        logger.info ('S1xB2 : %20s' % result ( model.sb ) [0] )
        logger.info ('B1xS2 : %20s' % result ( model.bs ) [0] )
        logger.info ('B1xB2 : %20s' % result ( model.bb ) [0] )

    models.add ( model ) 

# =============================================================================
## gauss as signal, 1st order polynomial as background 
# =============================================================================
def test_p1xp1_BBss () :
    logger.info ('Symmetrised fit model with non-factorized symmetric background:  ( Gauss + P1 ) (x) ( Gauss + P1 ) + BBsym' )
    sb      = ROOT.RooRealVar('sb','SB',0,10000)
    model   = Models.Fit2D (
        suffix   = '_5' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        power1   = 1  , 
        power2   = 1  ,
        bkg2D    = Models.PolyPos2Dsym_pdf ( 'P2Ds' , m_x , m_y , n = 2 ) ,
        sb       = sb ,
        bs       = sb 
        )

    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    model.bkg2   .tau  .fix ( 0 )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )
        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' % result ( model.ss ) [0] )
        logger.info ('S1xB2 : %20s' % result ( model.sb ) [0] )
        logger.info ('B1xS2 : %20s' % result ( model.bs ) [0] )
        logger.info ('B1xB2 : %20s' % result ( model.bb ) [0] )

    models.add ( model ) 

# =============================================================================
## gauss as signal, 1st order polynomial as background 
# =============================================================================
def test_p1xp1_BBsym () :
    logger.info ('Symmetric non-factorized fit model:  ( Gauss + P1 ) (x) ( Gauss + P1 ) + BBsym' )
    sb      = ROOT.RooRealVar('sb','SB',0,10000)
    model   = Models.Fit2DSym (
        suffix   = '_6' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        bkg1     = 1 , 
        bkg2D    = Models.PolyPos2Dsym_pdf ( 'P2D5' , m_x , m_y , n = 2 ) ,
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    model.bkg2   .tau  .fix ( 0 )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' %   result ( model.ss ) [0]   )
        logger.info ('S1xB2 : %20s' % ( result ( model.sb ) [0]/2 ) ) 
        logger.info ('B1xS2 : %20s' % ( result ( model.bs ) [0]/2 ) ) 
        logger.info ('B1xB2 : %20s' %   result ( model.bb ) [0]   )

    models.add ( model ) 


# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_pbxpb_BB  () :
    logger.info ('Non-factorizeable background component:  ( Gauss + expo*P1 ) (x) ( Gauss + expo*P1 ) + (expo*P1)**2')
    model   = Models.Fit2D (
        suffix   = '_7' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        power1   = 1 , 
        power2   = 1 ,
        bkg2D    = Models.ExpoPol2D_pdf ( 'P2D7' , m_x , m_y , nx = 1 , ny = 1 ) 
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    model.bkg2   .tau  .fix ( 0 )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' %   result ( model.ss ) [0]   )
        logger.info ('S1xB2 : %20s' %   result ( model.sb ) [0]   ) 
        logger.info ('B1xS2 : %20s' %   result ( model.bs ) [0]   ) 
        logger.info ('B1xB2 : %20s' %   result ( model.bb ) [0]   )

    models.add ( model ) 

    
# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_pbxpb_BBs () :
    logger.info ('Non-factorizeable background component:  ( Gauss + expo*P1 ) (x) ( Gauss + expo*P1 ) + Sym(expo*P1)**2')
    model   = Models.Fit2D (
        suffix   = '_8' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        power1   = 1 , 
        power2   = 1 ,
        bkg2D    = Models.ExpoPol2Dsym_pdf ( 'P2D8' , m_x , m_y , n = 1 ) 
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    model.bkg2   .tau  .fix ( 0 )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' %   result ( model.ss ) [0]   )
        logger.info ('S1xB2 : %20s' %   result ( model.sb ) [0]   ) 
        logger.info ('B1xS2 : %20s' %   result ( model.bs ) [0]   ) 
        logger.info ('B1xB2 : %20s' %   result ( model.bb ) [0]   )

    models.add ( model ) 
 
# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_pbxpb_BBsym () :
    logger.info ('Symmetric fit model with non-factorizeable background component:  ( Gauss + expo*P1 ) (x) ( Gauss + expo*P1 ) + Sym(expo*P1)**2')
    model   = Models.Fit2DSym (
        suffix   = '_9' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        bkg1     = 1 , 
        bkg2D    = Models.ExpoPol2Dsym_pdf ( 'P2D9' , m_x , m_y , n = 1 ) 
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    
    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' %   result ( model.ss ) [0]   )
        logger.info ('S1xB2 : %20s' %   ( result ( model.sb ) [0]  / 2 ) )
        logger.info ('B1xS2 : %20s' %   ( result ( model.bs ) [0]  / 2 ) )
        logger.info ('B1xB2 : %20s' %   result ( model.bb ) [0]   )

    models.add ( model ) 
                     
# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_psxps_BBs () :
        
    logger.info ('Non-factorizeable symmetric background component:  ( Gauss + expo*P1 ) (x) ( Gauss + expo*P1 ) + (PS*P1)**2')
    PS      = Ostap.Math.PhaseSpaceNL( 1.0  , 5.0 , 2 , 5 )
    model   = Models.Fit2D (
        suffix   = '_11' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        power1   = 1 , 
        power2   = 1 ,
        bkg2D    = Models.PSPol2Dsym_pdf ( 'P2D11' , m_x , m_y , ps = PS , n = 1 ) 
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    model.bkg2   .tau  .fix ( 0 )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' %   result ( model.ss ) [0]   )
        logger.info ('S1xB2 : %20s' %   result ( model.sb ) [0]   ) 
        logger.info ('B1xS2 : %20s' %   result ( model.bs ) [0]   ) 
        logger.info ('B1xB2 : %20s' %   result ( model.bb ) [0]   )

    models.add ( model ) 

# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_psxps_BBsym () :
    logger.info ('Simmetric fit model with non-factorizeable background component:  ( Gauss + expo*P1 ) (x) ( Gauss + expo*P1 ) + (PS*P1)**2')
    PS      = Ostap.Math.PhaseSpaceNL( 1.0  , 5.0 , 2 , 5 )
    model   = Models.Fit2DSym (
        suffix   = '_12' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        bkg1     = 1 , 
        bkg2D    = Models.PSPol2Dsym_pdf ( 'P2D12' , m_x , m_y , ps = PS , n = 1 ) 
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    model.bkg2   .tau  .fix ( 0 )
    
    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' %   result ( model.ss ) [0]   )
        logger.info ('S1xB2 : %20s' %   ( result ( model.sb ) [0]  / 2 ) )
        logger.info ('B1xS2 : %20s' %   ( result ( model.bs ) [0]  / 2 ) )
        logger.info ('B1xB2 : %20s' %   result ( model.bb ) [0]   )

    models.add ( model ) 
    
# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_model_13 () :
    
    logger.info ('Non-factorizeable fit component:  ( Gauss + expo*P1 ) (x) ( Gauss + expo*P1 ) + (Expo*PS)**2')
    PS      = Ostap.Math.PhaseSpaceNL( 1.0  , 5.0 , 2 , 5 )
    model   = Models.Fit2D (
        suffix   = '_13' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        power1   = 1 , 
        power2   = 1 ,
        bkg2D    = Models.ExpoPSPol2D_pdf ( 'P2D13' , m_x , m_y , psy = PS , nx = 1 , ny = 1 ) 
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    model.bkg2   .tau  .fix ( 0 )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' %   result ( model.ss ) [0]   )
        logger.info ('S1xB2 : %20s' %   result ( model.sb ) [0]   ) 
        logger.info ('B1xS2 : %20s' %   result ( model.bs ) [0]   ) 
        logger.info ('B1xB2 : %20s' %   result ( model.bb ) [0]   )

    models.add ( model ) 



    
knots = std.vector('double')()
knots.push_back(      m_x.xmin()              )
knots.push_back( 0.5*(m_x.xmin()+m_x.xmax() ) )
knots.push_back(                 m_x.xmax()   )
spline1 = Ostap.Math.BSpline( knots , 2 )

# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_model_14 () :
    
    logger.info ('Non-factorazeable background component (spline):  ( Gauss + expo*P1 ) (x) ( Gauss + expo*P1 ) + Spline2D')
    SPLINE  = Ostap.Math.Spline2D ( spline1 , spline1 ) 
    model   = Models.Fit2D (
        suffix   = '_14' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        power1   = 1 , 
        power2   = 1 ,
        bkg2D    = Models.Spline2D_pdf ( 'P2D14' , m_x , m_y , spline = SPLINE ) 
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    model.bkg2   .tau  .fix ( 0 )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' %   result ( model.ss ) [0]   )
        logger.info ('S1xB2 : %20s' %   result ( model.sb ) [0]   ) 
        logger.info ('B1xS2 : %20s' %   result ( model.bs ) [0]   ) 
        logger.info ('B1xB2 : %20s' %   result ( model.bb ) [0]   )

    models.add ( model )
    
# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_model_15 () :
    
    logger.info ('Non-factorized symmetric background component (spline):  ( Gauss + expo*P1 ) (x) ( Gauss + expo*P1 ) + Spline2Dsym')
    SPLINES = Ostap.Math.Spline2DSym ( spline1 ) 
    model   = Models.Fit2D (
        suffix   = '_15' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        power1   = 1 , 
        power2   = 1 ,
        bkg2D    = Models.Spline2Dsym_pdf ( 'P2D15' , m_x , m_y , spline = SPLINES ) 
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )
    model.bkg2   .tau  .fix ( 0 )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' %   result ( model.ss ) [0]   )
        logger.info ('S1xB2 : %20s' %   result ( model.sb ) [0]   ) 
        logger.info ('B1xS2 : %20s' %   result ( model.bs ) [0]   ) 
        logger.info ('B1xB2 : %20s' %   result ( model.bb ) [0]   )

    models.add ( model )

   
# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_model_16() : 
    logger.info ('Symmetric fit model with non-factorazeable symmetric (spline) backgrond component:  ( Gauss + expo*P1 ) (x) ( Gauss + expo*P1 ) + Spline2Dsym')
    SPLINES = Ostap.Math.Spline2DSym ( spline1 ) 
    model   = Models.Fit2DSym (
        suffix   = '_16' , 
        signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
        signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
        bkg1     = 1 , 
        bkg2D    = Models.Spline2Dsym_pdf ( 'P2D16' , m_x , m_y , spline = SPLINES ) 
        )
    
    model.signal1.sigma.fix ( m.error () )
    model.signal2.sigma.fix ( m.error () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.signal1.mean .fix ( m.value () )
    model.signal2.mean .fix ( m.value () )
    model.bkg1   .tau  .fix ( 0 )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal1.sigma.release () 
        model.signal2.sigma.release ()
        model.signal1.mean .release () 
        model.signal2.mean .release () 
        result, frame = model. fitTo ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print result 
    else :
        logger.info ('S1xS2 : %20s' %   result ( model.ss ) [0]      )
        logger.info ('S1xB2 : %20s' % ( result ( model.sb ) [0] /2 ) ) 
        logger.info ('B1xS2 : %20s' % ( result ( model.bs ) [0] /2 ) ) 
        logger.info ('B1xB2 : %20s' %   result ( model.bb ) [0]      )

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
        db['models'  ] = models
        db['dataset' ] = dataset
        db.ls()
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.timing import timing
    ## with timing ('test_const'       ) : test_const       ()          
    ## with timing ('test_p2xp2'       ) : test_p2xp2       ()          
    ## with timing ('test_p1xp1_BB'    ) : test_p1xp1_BB    ()          
    ## with timing ('test_p1xp1_BBss'  ) : test_p1xp1_BBss  ()          
    ## with timing ('test_p1xp1_BBsym' ) : test_p1xp1_BBsym ()          
    ## with timing ('test_pbxpb_BB'    ) : test_pbxpb_BB    ()          
    ## with timing ('test_pbxpb_BBs'   ) : test_pbxpb_BBs   ()          
    ## with timing ('test_pbxpb_BBsym' ) : test_pbxpb_BBsym ()          
    ## with timing ('test_psxps_BBs'   ) : test_psxps_BBs   ()          
    ## with timing ('test_psxps_BBsym' ) : test_psxps_BBsym ()          
    ## with timing ('test_model_13'    ) : test_model_13    ()          
    ## with timing ('test_model_14'    ) : test_model_14    ()          
    with timing ('test_model_15'    ) : test_model_15    ()          
    with timing ('test_model_16'    ) : test_model_16    ()          

    ## check finally that everything is serializeable:
    test_db ()          
    
# =============================================================================
# The END 
# =============================================================================
