#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_spectra.py
# Test module for ostap/fitting/spectra.py
# ============================================================================= 
""" Test module for ostap/fitting/pectra.py
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.core          import cpp, VE, dsID, rooSilent 
from   ostap.utils.timing       import timing
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait, batch_env 
from   ostap.logger.colorized   import attention
from   ostap.fitting.models     import ( Hagedorn_pdf  ,
                                         Tsallis_pdf   ,
                                         QGSM_pdf      ,
                                         GammaDist_pdf ) 
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_spectra' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
#

models = set() 
GeV    = 1.0
MeV    = GeV/1000

## fitting variable 
pt     = ROOT.RooRealVar ( 'pT' , 'Pt' , 0 , 10 * GeV )

## particle mass 
mass = 139 * MeV
mass =   3 * GeV 

## inverse temperature 
beta   = 1 / ( 200 * MeV ) 


## create the model 
model = Tsallis_pdf ( name = 'T0'   ,
                      xvar = pt     ,
                      m0   = mass ,
                      n    = 50     ,  
                      T    = 1/beta )
    
dataset = model.generate ( 100000 )

models.add ( model )

# =============================================================================
## test Hagedorn model
def test_hagedorn ( ) :
    """Test Hagedorn model"""
    

    logger = getLogger ( 'test_hagedorn' )
    logger.info ( 'Test Hagedorn model'  ) 
    
    ## create the model 
    model = Hagedorn_pdf ( name = 'H',
                           xvar = pt ,
                           m0   = ( mass , mass / 100 , mass * 100 ) ,
                           beta = ( beta   , beta   / 100 , beta   * 100 ) )
    

    with rooSilent() : 
        fitresult, _ = model.fitTo ( dataset , silent = True ) 
        fitresult, _ = model.fitTo ( dataset , silent = True ) 
        
    title = 'Hagedorn model'
    table = fitresult.table ( title  = title , prefix = '# ' ) 

    if 0 != fitresult.status() or 3 != fitresult.covQual() :
        logger.warning ('%s: fit result\n%s' % ( title , table ) )
    else :
        logger.info    ('%s: fit result\n%s' % ( title , table ) )

    with wait ( 2 ) , use_canvas ( 'Hagedorn model' ) :
        model.draw ( dataset ) 

    models.add ( model )

# =============================================================================
## test Tsallis model
def test_tsallis ( ) :
    """Test Tsallis model"""
    

    logger = getLogger ( 'test_tsallis' )
    logger.info ( 'Test Tsallis model'  ) 
    
    ## create the model 
    model = Tsallis_pdf ( name = 'T',
                          xvar = pt ,
                          m0   = ( mass , mass / 100 , mass * 100 ) ,
                          n    = ( 50 , 2  , 1.e+6  ) , 
                          T    = ( 1/beta   , (1/beta)   / 100 , (1/beta) * 100 ) )
    
    
    with rooSilent() : 
        fitresult, _ = model.fitTo ( dataset , silent = True ) 
        fitresult, _ = model.fitTo ( dataset , silent = True ) 
        
    title = 'Tsallis model'
    table = fitresult.table ( title  = title , prefix = '# ' ) 

    if 0 != fitresult.status() or 3 != fitresult.covQual() :
        logger.warning ('%s: fit result\n%s' % ( title , table ) )
    else :
        logger.info    ('%s: fit result\n%s' % ( title , table ) )

    with wait ( 2 ) , use_canvas ( 'Tsallis model' ) :
        model.draw ( dataset ) 

    models.add ( model )


# =============================================================================
## test QGSM model
def test_qgsm  ( ) :
    """Test QGSM model"""
    

    logger = getLogger ( 'test_qgsm' )
    logger.info ( 'Test QGSM model'  ) 
    
    ## create the model 
    model = QGSM_pdf ( name = 'Q',
                       xvar = pt ,
                       m0   = ( mass , mass / 100 , mass * 100 ) ,
                       b    = ( 1 , 1.e-6 , 1.e+6  ) )
    
    with rooSilent() : 
        fitresult, _ = model.fitTo ( dataset , silent = True ) 
        fitresult, _ = model.fitTo ( dataset , silent = True ) 
        
    title = 'QGSM model'
    table = fitresult.table ( title  = title , prefix = '# ' ) 

    if 0 != fitresult.status() or 3 != fitresult.covQual() :
        logger.warning ('%s: fit result\n%s' % ( title , table ) )
    else :
        logger.info    ('%s: fit result\n%s' % ( title , table ) )

    with wait ( 2 ) , use_canvas ( 'QGSM model' ) :
        model.draw ( dataset ) 

    models.add ( model )

# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        db['pt'       ] = pt
        db['dataset'  ] = dataset
        for m in models :
            db['model:' + m.name     ] = m
            db['roo_tot:%s' % m.name ] = m.pdf
        db['models'   ] = models
        db.ls() 

# =============================================================================
if '__main__' == __name__ :

    ## Hagehorn 
    with timing ('test_hagedorn' , logger ) :
        test_hagedorn ()
        
    ## Tsallis 
    with timing ('test_tsallis' , logger ) :
        test_tsallis  ()

    ## QGSM 
    with timing ('test_qgsm' , logger ) :
        test_qgsm     ()
        
    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()

# =============================================================================
##                                                                      The END 
# ============================================================================= 
