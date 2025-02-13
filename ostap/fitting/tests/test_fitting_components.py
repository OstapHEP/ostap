#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_components.py
# Test module 
# - It tests various multicomponents models 
# ============================================================================= 
""" Test module 
- It tests various multicomponents models 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.core      import cpp, VE, dsID, rooSilent 
from   ostap.utils.timing   import timing 
from   ostap.utils.utils    import batch_env 
import ostap.fitting.models as     Models 
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_components' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch form environment 
batch_env ( logger )
# =============================================================================

## make simple test mass 
mass    = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 0 , 15 )

## book very simple data set
varset  = ROOT.RooArgSet  ( mass )
dataset = ROOT.RooDataSet ( dsID() , 'Test Data set' , varset )  

mmin, mmax = mass.minmax()

m1 = VE(  4   , 0.125**2)
m2 = VE(  7.5 , 0.200**2)
m3 = VE( 11   , 0.125**2)

N1 = 5000
N2 = 1000
N3 = 1000

## fill it with three gausissians, 5k events each
for i in range ( N1 ) :
    for m in (m1,m2,m3) : 
        mass.value = m.gauss () 
        dataset.add ( varset     )


## add 5k events of uniform background
for i in range ( N2 ) :
    mass.value = random.uniform ( *mass.minmax() )
    dataset.add ( varset   )

## make background less trivial:
w1 = VE ( 6.0, 3**2 )
w2 = VE ( 9.0, 3**2 )
n1 = VE ( 3.0, 1**2 )
n2 = VE (12.0, 1**2 )

for i in range ( N3 ) :
    for w in ( w1 , w2 , n1 , n2 ) :
        v = w.gauss()
        if v in mass :
            mass.value =  v 
            dataset.add(varset)            

logger.info ('Dataset:\n%s' % dataset.table ( prefix = '# ' ) )   

NT = 3.0 * N1 + N2 + 4 * N3 

## various fit components

signal_1 = Models.Gauss_pdf ( 'G1'  , xvar = mass , mean = m1.value() , sigma = m1.error() )  
signal_2 = Models.Gauss_pdf ( 'G2'  , xvar = mass , mean = m2.value() , sigma = m2.error() ) 
signal_3 = Models.Gauss_pdf ( 'G3'  , xvar = mass , mean = m3.value() , sigma = m3.error() )

wide_1   = Models.Gauss_pdf ( 'GW1' , xvar = mass , mean = w1.value() , sigma = w1.error() )
wide_2   = Models.Gauss_pdf ( 'GW2' , xvar = mass , mean = w2.value() , sigma = w2.error() )

narrow_1 = Models.Gauss_pdf ( 'GN1' , xvar = mass , mean = n1.value() , sigma = n1.error() )
narrow_2 = Models.Gauss_pdf ( 'GN2' , xvar = mass , mean = n2.value() , sigma = n2.error() )


for s in ( signal_1 , signal_2  , signal_3 ) :
    s.mean  .fix()
    s.sigma.fix()

for s in ( wide_1 , wide_2 ) :
    s.mean  .fix()
    s.sigma.fix()

for s in ( narrow_1 , narrow_2 ) :
    s.mean  .fix()
    s.sigma.fix()
    
models = set()

# =============================================================================
## Test     extended multi-component fit'
def test_extended1 () :
    
    logger.info ('Test     extended multi-component fit')
    
    model = Models.Fit1D (
        name             = 'E1'                    , 
        signal           = signal_1                , 
        othersignals     = [ signal_2 , signal_3 ] ,
        background       = Models.Bkg_pdf ('P1' , xvar = mass , power = 0 , tau = 0 ) ,
        otherbackgrounds = [ wide_1   , wide_2   ] ,
        others           = [ narrow_1 , narrow_2 ] , 
        suffix           = '_a'
        )

    ## signals
    model.S = N1 , N1 , N1 
    
    ## backgrounds 
    model.B = N2 , N3 , N3 
    
    ## "components"
    model.C = N3 , N3 
    
    r, f = model.fitTo ( dataset , draw = False , silent = True )
    r, f = model.fitTo ( dataset , draw = True  , silent = True )
    
    logger.info ( 'Model %s Fit result\n%s' % ( model.name , r.table ( prefix = '# ') ) ) 

    models.add ( model ) 

# =============================================================================
## Test     extended combined multi-component fit'
def test_extended2 () :
    
    logger.info ('Test     extended combined multi-component fit')
    
    model = Models.Fit1D (
        name                = 'E2'                    , 
        signal              = signal_1                , 
        othersignals        = [ signal_2 , signal_3 ] ,
        background          = Models.Bkg_pdf ('P2' , xvar = mass , power = 0 , tau = 0 ) ,
        otherbackgrounds    = [ wide_1   , wide_2   ] ,
        others              = [ narrow_1 , narrow_2 ] ,
        combine_signals     =  True  , ## ATTENTION!
        combine_others      =  True  , ## ATTENTION! 
        combine_backgrounds =  True  , ## ATTENTION!
        suffix              = '_b'
        )

    model.S  = N1 + N1 + N1 
    model.B  = N2 + N3 + N3 
    model.C  = N3 + N3

    model.fS = N1 * 1.0 / ( N1 + N1 + N1 ) , N1 * 1.0 / ( N1 + N1 )
    model.fB = N2 * 1.0 / ( N2 + N3 + N3 ) , N3 * 1.0 / ( N3 + N3 )
    model.fC = N3 * 1.0 / ( N3 + N3 ) 
    
    r, f = model.fitTo ( dataset , draw = False , silent = True )        
    r, f = model.fitTo ( dataset , draw = True  , silent = True )
    
    logger.info ( 'Model %s Fit result\n%s' % ( model.name , r.table ( prefix = '# ' ) ) ) 

    models.add ( model ) 

# ==============================================================================
## Test non-extended multi-component fit'
def test_nonextended1 () :
    
    logger.info ('Test non-extended multi-component fit')
    
    model = Models.Fit1D (
        name                = 'N1'                    , 
        signal              = signal_1                , 
        othersignals        = [ signal_2 , signal_3 ] ,
        background          = Models.Bkg_pdf ('P3' , xvar = mass , power = 0 , tau = 0 ) ,
        otherbackgrounds    = [ wide_1 , wide_2 ] ,
        others              = [ narrow_1 , narrow_2 ] , 
        extended            = False , ## ATTENTION! 
        suffix              = '_c'
        )

    model.F = 0.2 , 0.25 , 0.36 , 0.65 , 0.05 , 0.25 , 0.50 

    r, f = model.fitTo ( dataset , draw = False , silent = True )
    r, f = model.fitTo ( dataset , draw = True  , silent = True )
    
    logger.info ( 'Model %s Fit result\n%s' % ( model.name , r.table ( prefix = '# ' ) ) ) 

    models.add ( model ) 

# ==============================================================================
## Test non-extended combined multi-component fit
def test_nonextended2 () :
    
    logger.info ('Test non-extended combined multi-component fit')
    
    model = Models.Fit1D (
        name                = 'N2'                    , 
        signal              = signal_1                , 
        othersignals        = [ signal_2 , signal_3 ] ,
        background          = Models.Bkg_pdf ('P4' , xvar = mass , power = 0 , tau = 0 ) ,
        otherbackgrounds    = [ wide_1 , wide_2 ] ,
        others              = [ narrow_1 , narrow_2   ] , 
        combine_signals     =  True  , ## ATTENTION!
        combine_others      =  True  , ## ATTENTION! 
        combine_backgrounds =  True  , ## ATTENTION!   
        extended            = False  , ## ATTENTION!
        suffix              = '_d'
        )

    model.F = N1 * 3.0 / ( 3 * N1 + N2 + 4 * N3 ) , ( N2 + N3 + N3 ) * 1.0/ ( N2 + 4 * N3 )
    
    model.fS = N1 * 1.0 / ( N1 + N1 + N1 ) , N1 * 1.0 / ( N1 + N1 )
    model.fB = N2 * 1.0 / ( N2 + N3 + N3 ) , N3 * 1.0 / ( N3 + N3 )
    model.fC = N3 * 1.0 / ( N3 + N3 ) 

    r, f = model.fitTo ( dataset , draw = False , silent = True )        
    r, f = model.fitTo ( dataset , draw = True  , silent = True )        
        
    logger.info ( 'Model %s Fit result\n%s' % ( model.name , r.table ( prefix = '# ' ) ) ) 

    models.add ( model ) 

# ==============================================================================
## Test non-extended multi-component non-recursive fit'
def test_nonextended3 () :
    
    logger.info ('Test non-extended non-recursive multi-component fit')
    
    model = Models.Fit1D (
        name                = 'N3'                    , 
        signal              = signal_1                , 
        othersignals        = [ signal_2 , signal_3 ] ,
        background          = Models.Bkg_pdf ('P5' , xvar = mass , power = 0 , tau = 0 ) ,
        otherbackgrounds    = [ wide_1 , wide_2 ] ,
        others              = [ narrow_1 , narrow_2  ] , 
        extended            = False , ## ATTENTION!
        recursive           = False , ## ATTENTION!
        suffix              = "_e"
        )
    

    with rooSilent() :
        
        model.F = N1 / NT , N1 / NT , N1 / NT , \
                  N2 / NT , N3 / NT , N3 / NT , \
                  N3 / NT , N3 / NT 
        
        r, f = model.fitTo ( dataset , draw = False , silent = True )
        r, f = model.fitTo ( dataset , draw = True , silent = True )
        
    logger.info ( 'Model %s Fit result\n%s' % ( model.name , r.table ( prefix = '# ' ) ) ) 

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
        db['mass'] = mass 
        db['vars'] = varset
        db['dataset'  ] = dataset
        for m in models :
            db['model:' + m.name ] = m
            db['roo_tot:%s' % m.name ] = m.pdf
            for i,s in enumerate ( m.signals ) :
                db['roo_sig%d:%s' % ( i , m.name ) ] = s
            for i, b in enumerate ( m.backgrounds ) : 
                db['roo_bkg%d:%s' % ( i , m.name ) ] = s
            for a in m.alist1 : 
                db['cmp:%s/%s' % ( m.name , a.name ) ] = a
        db['models'   ] = models
        db.ls() 

# =============================================================================
if '__main__' == __name__ :

    with timing ( "Extended1"     , logger ) : 
        test_extended1    ()
        
    with timing ( "Extended2"     , logger ) : 
        test_extended2    ()
        
    with timing ( "non-Extendd1" , logger ) : 
        test_nonextended1 ()
        
    with timing ( "non-Extended2" , logger ) : 
        test_nonextended2 () 
        
    with timing ( "non-Extended3" , logger ) : 
        test_nonextended3 () 

    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()

# =============================================================================
##                                                                      The END 
# =============================================================================
