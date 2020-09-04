#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_models6_2D.py
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
    logger = getLogger ( 'test_fitting_models6_2D' )
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


    
knots = std.vector('double')()
knots.push_back (      m_x.xmin()              )
knots.push_back ( 0.5*(m_x.xmin()+m_x.xmax() ) )
knots.push_back (                 m_x.xmax()   )
spline1 = Ostap.Math.BSpline ( knots , 2 )

   
# =============================================================================
## gauss as signal, expo times 1st order polynomial as background 
# =============================================================================
def test_model_16() : 
    logger.info ('Symmetric fit model with non-factorazeable symmetric (spline) background component:  ( Gauss + expo*P1 ) (x) ( Gauss + expo*P1 ) + Spline2Dsym')
    SPLINES = Ostap.Math.PositiveSpline2DSym ( spline1 ) 
    model   = Models.Fit2DSym (
        suffix   = '_16' , 
        signal_x = signal1  ,
        signal_y = signal2s ,
        bkg_1x   = 1 , 
        bkg_2D   = Models.Spline2Dsym_pdf ( 'P2D16' , m_x , m_y , spline = SPLINES ) 
        )

    ## fit with fixed mass and sigma
    with rooSilent() : 
        result, frame = model. fitTo ( dataset )
        model.signal_x.sigma.release () 
        model.signal_y.sigma.release ()
        model.signal_x.mean .release () 
        model.signal_y.mean .release () 
        result, frame = model. fitTo ( dataset )
        model.draw1 ( dataset )        
        model.draw2 ( dataset )

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
## check that everything is serializable
# =============================================================================
def test_db() :

    from ostap.core.meta_info import root_version_int 
    if root_version_int >= 62200 :
        logger.warning("test_db: test is disabled for ROOT verison %s" % root_version_int )
        return 
    
    logger.info('Saving all objects into DBASE')
    with timing('Save everything to DBASE', logger ), DBASE.tmpdb() as db : 
        db['m_x'     ] = m_x
        db['m_y'     ] = m_y
        db['vars'    ] = varset
        for m in models : db['model:' + m.name ] = m
        db['models'  ] = models
        db['dataset' ] = dataset
        db.ls()
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.timing import timing

    with timing ('test_model_16'   , logger ) : test_model_16    ()          
    
    ## check finally that everything is serializeable:
    with timing ( 'save to DB'     , logger ) : test_db ()          
    
# =============================================================================
##                                                                      The END 
# =============================================================================
