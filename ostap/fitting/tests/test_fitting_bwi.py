#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_bwi.py
# test for some Breit-Wigner models 
# ============================================================================= 
""" Test module for Breit-Wigner wirth interferenc e 
"""
# ============================================================================= 
from   __future__        import print_function
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models  as     Models 
from   ostap.core.core       import Ostap
import ostap.io.zipshelve    as     DBASE
from   ostap.utils.timing    import timing
from   ostap.utils.utils     import wait
from   ostap.plotting.canvas import use_canvas 
from   builtins              import range
import ROOT, time 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_bwi' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

m1    = 1.0 
m2    = 1.0
m0    = 5.0
g0    = 1.0 

## create BreitWigner function 
breitwigner = Ostap.Math.BreitWigner ( m0 , g0 , m1 , m2  )
phasespace  = Ostap.Math.PhaseSpace2 ( m1 , m2 )

## mass observable 
mass       = ROOT.RooRealVar ( 'mass' , 'mass-observable' , 0 , 15 )  

## create BW function
bw        = Models.BreitWigner_pdf ( 'BW'  ,
                                     m0     = ( m0 , m0 - 1 , m0 + 1 ) , 
                                     gamma  = ( g0 , g0/2   , g0 * 2 ) , 
                                     breitwigner = breitwigner , xvar = mass )

## create BWI function
bwi        = Models.BWI_pdf        ( 'BWI' ,
                                     m0    = bw.m0    ,
                                     gamma = bw.gamma ,                                     
                                     breitwigner = breitwigner , xvar = mass )
bwi.magnitude.phis =  0.05 , 
bwi.phase    .phis = -0.7  , 
bwi.scale1         =  0.3
bwi.scale2         = 11.0 

## model without interference
bkg1   = Models.PSLeftExpoPol_pdf ( 'B1' , xvar = mass       , 
                                    phasespace  = phasespace ,
                                    scale       = ROOT.RooFit.RooConst ( 1 ) ,
                                    tau         = ROOT.RooFit.RooConst ( 0 ) ,
                                    power       = 5 )
model1   = Models.Fit1D ( signal = bw , background = bkg1 , suffix = '1' ) 
model1.S = 100
model1.B =  50

## model with interference
bkg2   = Models.PSLeftExpoPol_pdf ( 'B2' , xvar = mass       ,
                                    phasespace  = phasespace ,
                                    scale       = ROOT.RooFit.RooConst ( 1 ) ,
                                    tau         = ROOT.RooFit.RooConst ( 0 ) ,
                                    power       = 1 )
model2   = Models.Fit1D ( signal = bwi , background = bkg2 , suffix = '2' ) 
model2.S = 630
model2.B = 330

ds1 = model1.generate ( 1000 )
ds2 = model2.generate ( 1000 )

models  = set()
results = [] 

models.add  ( model1 ) 
models.add  ( model2 ) 
models.add  ( model1.signal ) 
models.add  ( model2.signal ) 
models.add  ( model1.background ) 
models.add  ( model2.background ) 

# ==============================================================================
with wait ( 1 ), use_canvas ( 'Fit ds1 with model1' ) : 
    r11 , f = model1.fitTo ( ds1 , silent = True )
    r11 , f = model1.fitTo ( ds1 , silent = True , draw = True )
    title = 'Fit ds1 with model1'
    logger.info ( '%s\n%s' % ( title , r11.table ( title , prefix = '# ' ) ) )
    results.append ( r11 ) 
    
with wait ( 1 ), use_canvas ( 'Fit ds2 with model2' ) : 
    r22 , f = model2.fitTo ( ds2 , silent = True )
    r22 , f = model2.fitTo ( ds2 , silent = True , draw = True )
    title = 'Fit ds2 with model2'
    logger.info ( '%s\n%s' % ( title , r22.table ( title , prefix = '# ' ) ) ) 
    results.append ( r22 ) 

with wait ( 1 ), use_canvas ( 'Fit ds1 with model2' ) : 
    r12 , f = model2.fitTo ( ds1 , silent = True )
    r12 , f = model2.fitTo ( ds1 , silent = True , draw = True )
    title = 'Fit ds1 with model2'
    logger.info ( '%s\n%s' % ( title , r12.table ( title , prefix = '# ' ) ) ) 
    results.append ( r12 ) 

with wait ( 1 ), use_canvas ( 'Fit ds2 with model1' ) : 
    r21 , f = model1.fitTo ( ds2 , silent = True )
    r21 , f = model1.fitTo ( ds2 , silent = True , draw = True )
    title = 'Fit ds2 with model1'
    logger.info ( '%s\n%s' % ( title , r21.table ( title , prefix = '# ' ) ) ) 
    results.append ( r21 ) 

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
        db['dataset1' ] = ds1
        db['dataset2' ] = ds2
        for m in models :
            db['model:' + m.name ] = m
            db['roo:%s' % m.name ] = m.pdf
        db['models'   ] = models
        for i, r in enumerate ( results ) :
            db ['result:%s' % r.name ] = r
        db['results'   ] = results 
        db.ls() 


# =============================================================================
if '__main__' == __name__ :
    
    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()

# =============================================================================
##                                                                      The END 
# =============================================================================
