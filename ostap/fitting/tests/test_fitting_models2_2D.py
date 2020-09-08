#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_models2_2D.py
# Test module for ostap/fitting/models_2d.py
# - It tests various 2D-non-factrorizeable models 
# ============================================================================= 
""" Test module for ostap/fitting/models_2d.py
- It tests various 2D-non-factrorizeable models 
"""
# ============================================================================= 
from   __future__        import print_function
# ============================================================================= 
import ROOT, random
import ostap.fitting.roofit 
import ostap.fitting.models as     Models 
from   ostap.core.core      import Ostap, std, VE, dsID
from   ostap.logger.utils   import rooSilent 
import ostap.io.zipshelve   as     DBASE
from   ostap.utils.timing   import timing 
from   builtins             import range
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_models2_2D' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
root_version = ROOT.ROOT.GetROOT().GetVersionInt()
# =============================================================================
## make simple test mass 
m_x     = ROOT.RooRealVar ( 'mass_x' , 'Some test mass(X)' , 3 , 3.2 )
m_y     = ROOT.RooRealVar ( 'mass_y' , 'Some test mass(Y)' , 3 , 3.2 )

## book very simple data set
varset  = ROOT.RooArgSet  ( m_x , m_y )
dataset = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset )  

m = VE(3.100,0.015**2)

N_ss = 5000
N_sb =  500
N_bs =  500
N_bb = 1000

random.seed(0)

## fill it : 5000 events  Gauss * Gauss 
for i in range(0,N_ss) : 
    m_x.value = m.gauss() 
    m_y.value = m.gauss() 
    dataset.add ( varset  )


## fill it : 2500 events  Gauss * const  
for i in range(0,N_sb) : 
    m_x.value = m.gauss() 
    m_y.value = random.uniform ( *m_y.minmax() )  
    dataset.add ( varset  )

## fill it : 2500 events  const * Gauss
for i in range(0,N_bs) : 
    m_x.value = random.uniform ( *m_x.minmax() ) 
    m_y.value = m.gauss() 
    dataset.add ( varset  )

## fill it : 5000 events  const * const
for i in range(0,N_bb) :

    m_x.value  = random.uniform ( *m_x.minmax() )
    m_y.value  = random.uniform ( *m_y.minmax() )
    dataset.add ( varset  )

logger.info ( 'Dataset:%s ' % dataset ) 


models = set()
# =============================================================================


signal1  = Models.Gauss_pdf ( 'Gx'  , xvar = m_x ) 
signal2  = Models.Gauss_pdf ( 'Gy'  , xvar = m_y )
signal2s = signal1.clone ( name = 'GyS' , xvar = m_y )

signal1.mean  = m.value ()
signal1.sigma = m.error ()
signal2.mean  = m.value ()
signal2.sigma = m.error ()


# =============================================================================
## gauss as signal, const as background 
# =============================================================================
def test_const () :
    
    logger.info ('Simplest (factorized) fit model:  ( Gauss + const ) x ( Gauss + const ) ' )
    model   = Models.Fit2D (
        signal_x = signal1  ,
        signal_y = signal2s ,
        )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal_x.sigma.release () 
        model.signal_y.sigma.release ()
        model.signal_x.mean .release () 
        model.signal_y.mean .release () 
        result, frame = model. fitTo ( dataset )
        model. draw1 ( dataset )
        model. draw2 ( dataset )

        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        logger.info ('S1xS2 : %20s' % result ( model.SS ) [0] )
        logger.info ('S1xB2 : %20s' % result ( model.SB ) [0] )
        logger.info ('B1xS2 : %20s' % result ( model.BS ) [0] )
        logger.info ('B1xB2 : %20s' % result ( model.BB ) [0] )

    models.add ( model ) 

# =============================================================================
## gauss as signal, second order polynomial as background 
# =============================================================================
def test_p2xp2 () : 
    logger.info ('Simple (factorized) fit model:  ( Gauss + P1 ) (x) ( Gauss + P1 ) ' )
    model   = Models.Fit2D (
        suffix   = '_2' , 
        signal_x = signal1  ,
        signal_y = signal2s ,
        bkg_1x     = -1 , 
        bkg_1y     = -1 ,
        bkg_2x     = -1 ,
        bkg_2y     = -1 , 
        )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal_x.sigma.release () 
        model.signal_y.sigma.release ()
        model.signal_x.mean .release () 
        model.signal_y.mean .release () 
        result, frame = model. fitTo ( dataset )
        model. draw1 ( dataset )
        model. draw2 ( dataset )
        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        logger.info ('S1xS2 : %20s' % result ( model.SS ) [0] )
        logger.info ('S1xB2 : %20s' % result ( model.SB ) [0] )
        logger.info ('B1xS2 : %20s' % result ( model.BS ) [0] )
        logger.info ('B1xB2 : %20s' % result ( model.BB ) [0] )

    models.add ( model ) 


# =============================================================================
## gauss as signal, 1st order polynomial as background + non-factorizeable BB
# =============================================================================
def test_p1xp1_BB () : 
    logger.info ('Simplest non-factorized fit model:  ( Gauss + P1 ) (x) ( Gauss + P1 ) + BB' )
    model   = Models.Fit2D (
        suffix   = '_3' , 
        signal_x = signal1  ,
        signal_y = signal2s ,
        bkg_1x     = -1 , 
        bkg_1y     = -1 ,
        bkg_2D     = Models.PolyPos2D_pdf ( 'P2D' , m_x , m_y , nx = 2 , ny = 2 ) 
        )
    
    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal_x.sigma.release () 
        model.signal_y.sigma.release ()
        model.signal_x.mean .release () 
        model.signal_y.mean .release () 
        result, frame = model. fitTo ( dataset )
        model. draw1 ( dataset )
        model. draw2 ( dataset )
        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        logger.info ('S1xS2 : %20s' % result ( model.SS ) [0] )
        logger.info ('S1xB2 : %20s' % result ( model.SB ) [0] )
        logger.info ('B1xS2 : %20s' % result ( model.BS ) [0] )
        logger.info ('B1xB2 : %20s' % result ( model.BB ) [0] )

    models.add ( model ) 


# =============================================================================
## gauss as signal, 1st order polynomial as background 
# =============================================================================
def test_p1xp1_BBs () :
    
    logger.info ('Non-factorized symmetric background fit model:  ( Gauss + P1 ) (x) ( Gauss + P1 ) + BBsym' )
    model   = Models.Fit2D (
        suffix   = '_4' , 
        signal_x = signal1  ,
        signal_y = signal2s ,
        bkg_1x     = -1 , 
        bkg_1y     = -1 ,
        bkg_2D    = Models.PolyPos2Dsym_pdf ( 'P2Ds' , m_x , m_y , n = 2 ) 
        )
    
    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal_x.sigma.release () 
        model.signal_y.sigma.release ()
        model.signal_x.mean .release () 
        model.signal_y.mean .release () 
        result, frame = model. fitTo ( dataset )
        model. draw1 ( dataset )
        model. draw2 ( dataset )
        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        logger.info ('S1xS2 : %20s' % result ( model.SS ) [0] )
        logger.info ('S1xB2 : %20s' % result ( model.SB ) [0] )
        logger.info ('B1xS2 : %20s' % result ( model.BS ) [0] )
        logger.info ('B1xB2 : %20s' % result ( model.BB ) [0] )

    models.add ( model ) 

# =============================================================================
## gauss as signal, 1st order polynomial as background 
# =============================================================================
def test_p1xp1_BBss () :
    
    logger.info ('Symmetrised fit model with non-factorized symmetric background:  ( Gauss + P1 ) (x) ( Gauss + P1 ) + BBsym' )
    sb      = ROOT.RooRealVar('sb','SB',2500 , 0,10000)
    model   = Models.Fit2D (
        suffix   = '_5' , 
        signal_x = signal1  ,
        signal_y = signal2s ,
        bkg_1x     = -1     , 
        bkg_2D    = Models.PolyPos2Dsym_pdf ( 'P2Ds' , m_x , m_y , n = 1 ) ,
        sb       = sb ,
        bs       = sb 
        )

    model.SS = N_ss
    model.BB = N_bb
    model.SB = 2500

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal_x.sigma.release () 
        model.signal_y.sigma.release ()
        model.signal_x.mean .release () 
        model.signal_y.mean .release () 
        result, frame = model. fitTo ( dataset )
        model. draw1 ( dataset )
        model. draw2 ( dataset )
        
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        logger.info ('S1xS2 : %20s' % result ( model.SS ) [0] )
        logger.info ('S1xB2 : %20s' % result ( model.SB ) [0] )
        logger.info ('B1xS2 : %20s' % result ( model.BS ) [0] )
        logger.info ('B1xB2 : %20s' % result ( model.BB ) [0] )

    models.add ( model ) 

# =============================================================================
## gauss as signal, 1st order polynomial as background 
# =============================================================================
def test_p1xp1_BBsym () :
    logger.info ('Symmetric non-factorized fit model:  ( Gauss + P1 ) (x) ( Gauss + P1 ) + BBsym' )
    sb      = ROOT.RooRealVar('sb','SB',0,10000)
    model   = Models.Fit2DSym (
        suffix   = '_6' , 
        signal_x = signal1  ,
        signal_y = signal2s ,
        bkg_1x     = -1 , 
        bkg_2D    = Models.PolyPos2Dsym_pdf ( 'P2D5' , m_x , m_y , n = 2 ) ,
        )
    
    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal_x.sigma.release () 
        model.signal_y.sigma.release ()
        model.signal_x.mean .release () 
        model.signal_y.mean .release () 
        result, frame = model. fitTo ( dataset )
        model. draw1 ( dataset )
        model. draw2 ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        logger.info ('S1xS2 : %20s' %   result ( model.SS ) [0]      )
        logger.info ('S1xB2 : %20s' % ( result ( model.SB ) [0] /2 ) ) 
        logger.info ('B1xS2 : %20s' % ( result ( model.BS ) [0] /2 ) ) 
        logger.info ('B1xB2 : %20s' %   result ( model.BB ) [0]      )

    models.add ( model ) 


# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_pbxpb_BB  () :
    logger.info ('Non-factorizeable background component:  ( Gauss + expo*P1 ) (x) ( Gauss + expo*P1 ) + (expo*P1)**2')
    model   = Models.Fit2D (
        suffix   = '_7' , 
        signal_x = signal1  ,
        signal_y = signal2s ,
        bkg_1x     = 1 , 
        bkg_1y     = 1 ,
        bkg_2D    = Models.ExpoPol2D_pdf ( 'P2D7' , m_x , m_y , nx = 1 , ny = 1 ) 
        )
    
    model.bkg_1x   .tau  .fix ( 0 )
    model.bkg_1y   .tau  .fix ( 0 )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal_x.sigma.release () 
        model.signal_y.sigma.release ()
        model.signal_x.mean .release () 
        model.signal_y.mean .release () 
        result, frame = model. fitTo ( dataset )
        model. draw1 ( dataset )
        model. draw2 ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print(result)
    else :
        logger.info ('S1xS2 : %20s' % result ( model.SS ) [0]    )
        logger.info ('S1xB2 : %20s' % result ( model.SB ) [0]    )
        logger.info ('B1xS2 : %20s' % result ( model.BS ) [0]    )
        logger.info ('B1xB2 : %20s' % result ( model.BB ) [0]    )

    models.add ( model ) 

    
# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_pbxpb_BBs () :
    logger.info ('Non-factorizeable background component:  ( Gauss + expo*P1 ) (x) ( Gauss + expo*P1 ) + Sym(expo*P1)**2')
    model   = Models.Fit2D (
        suffix   = '_8' , 
        signal_x = signal1  ,
        signal_y = signal2s ,
        bkg_1x     = 1 , 
        bkg_1y     = 1 ,
        bkg_2D    = Models.ExpoPol2Dsym_pdf ( 'P2D8' , m_x , m_y , n = 1 ) 
        )

    model.bkg_1x   .tau  .fix ( 0 )
    model.bkg_1y   .tau  .fix ( 0 )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal_x.sigma.release () 
        model.signal_y.sigma.release ()
        model.signal_x.mean .release () 
        model.signal_y.mean .release () 
        result, frame = model. fitTo ( dataset )
        model. draw1 ( dataset )
        model. draw2 ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print(result)
    else :

        logger.info ('S1xS2 : %20s' % result ( model.SS ) [0]    )
        logger.info ('S1xB2 : %20s' % result ( model.SB ) [0]    )
        logger.info ('B1xS2 : %20s' % result ( model.BS ) [0]    )
        logger.info ('B1xB2 : %20s' % result ( model.BB ) [0]    )


    models.add ( model ) 
 
# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_pbxpb_BBsym () :
    logger.info ('Symmetric fit model with non-factorizeable background component:  ( Gauss + P1 ) (x) ( Gauss + P1 ) + Sym(expo*P1)**2')
    model   = Models.Fit2DSym (
        suffix   = '_9' , 
        signal_x = signal1  ,
        signal_y = signal2s ,
        bkg_1x     = -1   , 
        bkg_2D    = Models.ExpoPol2Dsym_pdf ( 'P2D9' , m_x , m_y , n = 1 ) 
        )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal_x.sigma.release () 
        model.signal_y.sigma.release ()
        model.signal_x.mean .release () 
        model.signal_y.mean .release () 
        result, frame = model. fitTo ( dataset )
        model. draw1 ( dataset )
        model. draw2 ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print(result)
    else :

        logger.info ('S1xS2 : %20s' %   result ( model.SS ) [0]       )
        logger.info ('S1xB2 : %20s' % ( result ( model.SB ) [0] /2  ) )
        logger.info ('B1xS2 : %20s' % ( result ( model.BS ) [0] /2  ) )
        logger.info ('B1xB2 : %20s' %   result ( model.BB ) [0]       )

    models.add ( model ) 
                     
# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_psxps_BBs () :
        
    logger.info ('Non-factorizeable symmetric background component:  ( Gauss + P1 ) (x) ( Gauss + P1 ) + (PS*P1)**2')
    PS      = Ostap.Math.PhaseSpaceNL( 1.0  , 5.0 , 2 , 5 )
    model   = Models.Fit2D (
        suffix   = '_11' , 
        signal_x = signal1  ,
        signal_y = signal2s ,
        bkg_1x     = -1 , 
        bkg_1y     = -1 ,
        bkg_2D    = Models.PSPol2Dsym_pdf ( 'P2D11' , m_x , m_y , ps = PS , n = 1 ) 
        )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal_x.sigma.release () 
        model.signal_y.sigma.release ()
        model.signal_x.mean .release () 
        model.signal_y.mean .release () 
        result, frame = model. fitTo ( dataset )
        model. draw1 ( dataset )
        model. draw2 ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print(result)
    else :

        logger.info ('S1xS2 : %20s' % result ( model.SS ) [0]     )
        logger.info ('S1xB2 : %20s' % result ( model.SB ) [0]     )
        logger.info ('B1xS2 : %20s' % result ( model.BS ) [0]     )
        logger.info ('B1xB2 : %20s' % result ( model.BB ) [0]     )

    models.add ( model ) 

# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_psxps_BBsym () :
    logger.info ('Simmetric fit model with non-factorizeable background component:  ( Gauss + P1 ) (x) ( Gauss + P1 ) + (PS*P1)**2')
    PS      = Ostap.Math.PhaseSpaceNL( 1.0  , 5.0 , 2 , 5 )
    model   = Models.Fit2DSym (
        suffix   = '_12' , 
        signal_x = signal1  ,
        signal_y = signal2s ,
        bkg_1x     = -1 , 
        bkg_2D    = Models.PSPol2Dsym_pdf ( 'P2D12' , m_x , m_y , ps = PS , n = 1 ) 
        )
    
    
    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal_x.sigma.release () 
        model.signal_y.sigma.release ()
        model.signal_x.mean .release () 
        model.signal_y.mean .release () 
        result, frame = model. fitTo ( dataset )
        model. draw1 ( dataset )
        model. draw2 ( dataset )

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d '
                       % ( result.status() , result.covQual()  ) )
        print(result)
    else :

        logger.info ('S1xS2 : %20s' %   result ( model.SS ) [0]       )
        logger.info ('S1xB2 : %20s' % ( result ( model.SB ) [0]  /2 ) )
        logger.info ('B1xS2 : %20s' % ( result ( model.BS ) [0]  /2 ) ) 
        logger.info ('B1xB2 : %20s' %   result ( model.BB ) [0]       )

    models.add ( model ) 
    
 
# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    from ostap.core.meta_info import root_version_int
    if root_version_int >= 62200 :
        logger.warning("test_db: test is disabled for ROOT version %s" % root_version_int )
        return
    
    with timing ( 'Save everything to DBASE' , logger ), DBASE.tmpdb() as db : 
        db['m_x'     ] = m_x
        db['m_y'     ] = m_y
        db['vars'    ] = varset
        for m in models : db[ 'model:' + m.name ] = m 
        db['models'  ] = models
        db['dataset' ] = dataset
        db.ls()
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.timing import timing
    with timing ('test_const'       , logger ) : test_const       ()          
    with timing ('test_p2xp2'       , logger ) : test_p2xp2       ()          
    with timing ('test_p1xp1_BB'    , logger ) : test_p1xp1_BB    ()          
    with timing ('test_p1xp1_BBss'  , logger ) : test_p1xp1_BBss  ()          
    with timing ('test_p1xp1_BBsym' , logger ) : test_p1xp1_BBsym ()          
    with timing ('test_pbxpb_BB'    , logger ) : test_pbxpb_BB    ()          
    with timing ('test_pbxpb_BBs'   , logger ) : test_pbxpb_BBs   ()          
    with timing ('test_pbxpb_BBsym' , logger ) : test_pbxpb_BBsym ()          
    with timing ('test_psxps_BBs'   , logger ) : test_psxps_BBs   ()          
    with timing ('test_psxps_BBsym' , logger ) : test_psxps_BBsym ()          

    with timing ('Save to DB' ) : test_db ()          
    
# =============================================================================
##                                                                      The END 
# =============================================================================
