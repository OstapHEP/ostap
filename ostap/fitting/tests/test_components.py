#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_components.py
# Test module for ostap/fitting/models.py
# - It tests various multicomponents models 
# ============================================================================= 
""" Test module for ostap/fitting/models.py
- It tests various multicomponents models 
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
    logger = getLogger ( 'test_components' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

## make simple test mass 
mass    = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 0 , 10 )

## book very simple data set
varset  = ROOT.RooArgSet  ( mass )
dataset = ROOT.RooDataSet ( dsID() , 'Test Data set' , varset )  

mmin, mmax = mass.minmax()

m1 = VE(3,0.500**2)
m2 = VE(5,0.250**2)
m3 = VE(7,0.125**2)

## fill it with three gausissians, 5k events each
for i in xrange(0,5000) :
    for m in (m1,m2,m3) : 
        mass.value = m.gauss () 
        dataset.add ( varset     )


## add 5k events of uniform background
for i in xrange(0,5000) :
    mass.value = random.uniform ( *mass.minmax() )
    dataset.add ( varset   )

## make background less trivial:
w1 = VE(4.0, 4**2 )
w2 = VE(6.0, 4**2 )
n1 = VE(2.0, 1**2 )
n2 = VE(4.0, 1**2 )

for i in xrange(1000) :
    for w in ( w1 , w2 , n1 , n2 ) :
        v = w.gauss()
        if v in mass :
            mass.value =  v 
            dataset.add(varset)            

logger.info ('Dataset: %s' % dataset )  



## various fit components

signal_1 = Models.Gauss_pdf ( 'G1' , xvar = mass  , mean = m1.value() , sigma = m1.error() )  
signal_2 = Models.Gauss_pdf ( 'G2' , xvar = mass  , mean = m2.value() , sigma = m2.error() ) 
signal_3 = Models.Gauss_pdf ( 'G3' , xvar = mass  , mean = m3.value() , sigma = m3.error() )

wide_1   = Models.Gauss_pdf ( 'GW1', xvar = mass  , mean = 4.0  , sigma = 4 )
wide_2   = Models.Gauss_pdf ( 'GW2', xvar = mass  , mean = 6.0  , sigma = 4 )

narrow_1 = Models.Gauss_pdf ( 'GN1' , xvar = mass , mean = 2.0  , sigma =  1 )
narrow_2 = Models.Gauss_pdf ( 'GN2' , xvar = mass , mean = 4.0  , sigma =  1 )


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
        )
 
    with rooSilent() : 
    ## signals
        model.S[0].fix    ( 5000 )
        model.S[1].fix    ( 5000 )
        model.S[2].fix    ( 5000 )
    
        ## backgrounds 
        model.B[0].fix    ( 5000 )
        model.B[1].fix    ( 1000 )
        model.B[2].setVal ( 1000 )
        
        ## "components"
        model.C[0].value = 1000
        model.C[1].value =  500
        
        r, f = model.fitTo ( dataset , draw = False )

        model.S[0].release() 
        model.S[1].release()
        model.S[2].release()
        
        model.B[0].release()

        r, f = model.fitTo ( dataset , draw = False , silent = True )
        
        model.B[1].release()
        model.C[0].release() 
        
    r, f = model.fitTo ( dataset , draw = False , silent = True )
        
    logger.info ( 'Model %s Fit result \n#%s ' % ( model.name , r ) ) 


# =============================================================================
## Test     extended combined multi-component fit'
def test_extended2 () :
    
    logger.info ('Test     extended combined multi-component fit')
    
    model = Models.Fit1D (
        name                = 'E2'                    , 
        signal              = signal_1                , 
        othersignals        = [ signal_2 , signal_3 ] ,
        background          = Models.Bkg_pdf ('P1' , xvar = mass , power = 0 , tau = 0 ) ,
        otherbackgrounds    = [ wide_1   , wide_2   ] ,
        others              = [ narrow_1 , narrow_2 ] ,
        combine_signals     =  True  , ## ATTENTION!
        combine_others      =  True  , ## ATTENTION! 
        combine_backgrounds =  True  , ## ATTENTION!   
        )
    

    with rooSilent() :
        model.S.fix ( 15000 )
        model.B.fix (  7000 )
        model.C.fix (  2000 )
        r, f = model.fitTo ( dataset , draw = False , silent = True )
        
    model.S.release() 
    model.B.release() 
    model.C.release() 
    r, f = model.fitTo ( dataset , draw = False , silent = True )
    
    logger.info ( 'Model %s Fit result \n#%s ' % ( model.name , r ) ) 

# ==============================================================================
## Test non-extended multi-component fit'
def test_nonextended1 () :
    
    logger.info ('Test non-extended multi-component fit')
    
    model = Models.Fit1D (
        name                = 'N1'                    , 
        signal              = signal_1                , 
        othersignals        = [ signal_2 , signal_3 ] ,
        background          = Models.Bkg_pdf ('P2' , xvar = mass , power = 0 , tau = 0 ) ,
        otherbackgrounds    = [ wide_1 , wide_2 ] ,
        others              = [ narrow_1 , narrow_2 ] , 
        extended            = False  ## ATTENTION! 
        )

    model.F[0].setVal ( 0.20 )
    model.F[1].setVal ( 0.25 )
    model.F[2].setVal ( 0.36 )

    model.F[3].setVal ( 0.65 )
    model.F[4].setVal ( 0.05 )
    model.F[5].setVal ( 0.25 )
    model.F[6].setVal ( 0.50 )
    
    r, f = model.fitTo ( dataset , draw = False , silent = True )
    logger.info ( 'Model %s Fit result \n#%s ' % ( model.name , r ) ) 

# ==============================================================================
## Test non-extended combined multi-component fit
def test_nonextended2 () :
    
    logger.info ('Test non-extended combined multi-component fit')
    
    model = Models.Fit1D (
        name                = 'N2'                    , 
        signal              = signal_1                , 
        othersignals        = [ signal_2 , signal_3 ] ,
        background          = Models.Bkg_pdf ('P3' , xvar = mass , power = 0 , tau = 0 ) ,
        otherbackgrounds    = [ wide_1 , wide_2 ] ,
        others              = [ narrow_1 , narrow_2   ] , 
        combine_signals     =  True  , ## ATTENTION!
        combine_others      =  True  , ## ATTENTION! 
        combine_backgrounds =  True  , ## ATTENTION!   
        extended            = False  , ## ATTENTION!
        )
                  ## fB_1    3.3333e-01    9.4219e-01 +/-  7.64e-02  <none>
                  ## fB_2    3.3333e-01    4.6419e-04 +/-  6.25e-01  <none>
                  ## fC_1    5.0000e-01    4.9921e-01 +/-  3.95e-02  <none>
                  ## fS_1    3.3333e-01    3.3002e-01 +/-  5.35e-03  <none>
                  ## fS_2    3.3333e-01    4.9896e-01 +/-  5.60e-03  <none>
                  ##  f_1    3.3333e-01    6.3647e-01 +/-  6.05e-03  <none>
                  ##  f_2    3.3333e-01    7.3854e-01 +/-  1.71e-02  <none>

    with rooSilent() :
        model.F[0] .fix( 0.63 ) 
        model.F[1] .fix( 0.74 )
        
        model.fB[0].setVal( 0.95 )
        model.fB[1].setVal( 0.01 )
        model.fC[0].fix   ( 0.50 )
        model.fS[0].fix   ( 0.33 )
        model.fS[1].fix   ( 0.50 )
        r, f = model.fitTo ( dataset , draw = False , silent = True )
        
        model.F[0]. release() 
        model.F[1] .release() 
        
        model.fC[0].release() 
        model.fS[0].release() 
        model.fS[1].release()
        
        r, f = model.fitTo ( dataset , draw = False , silent = True )        
        
    logger.info ( 'Model %s Fit result \n#%s ' % ( model.name , r ) ) 

# ==============================================================================
## Test non-extended multi-component non-recursive fit'
def test_nonextended3 () :
    
    logger.info ('Test non-extended non-recursive multi-component fit')
    
    model = Models.Fit1D (
        name                = 'N3'                    , 
        signal              = signal_1                , 
        othersignals        = [ signal_2 , signal_3 ] ,
        background          = Models.Bkg_pdf ('P2' , xvar = mass , power = 0 , tau = 0 ) ,
        otherbackgrounds    = [ wide_1 , wide_2 ] ,
        others              = [ narrow_1 , narrow_2  ] , 
        extended            = False , ## ATTENTION!
        recursive           = False   ## ATTENTION! 
        )
    

    with rooSilent() :
        
        model.F[0].fix    ( 0.20 )
        model.F[1].fix    ( 0.20 )
        model.F[2].fix    ( 0.20 )
        model.F[3].fix    ( 0.25 )
        
        model.F[4].setVal ( 0.01  )
        model.F[5].setVal ( 0.02  )
        model.F[6].setVal ( 0.02  )
        
        r, f = model.fitTo ( dataset , draw = False , silent = True )
        
        model.F[0].release() 
        model.F[1].release() 
        model.F[2].release()
        model.F[3].release()  
        
        r, f = model.fitTo ( dataset , draw = False , silent = True )
        
    logger.info ( 'Model %s Fit result \n#%s ' % ( model.name , r ) ) 
    
# =============================================================================
if '__main__' == __name__ :

    test_extended1    () 
    test_extended2    () 
    test_nonextended1 () 
    test_nonextended2 () 
    test_nonextended3 () 
    
# =============================================================================
# The END 
# =============================================================================
