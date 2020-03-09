#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_efficiency.py
# Test module for ostap/fitting/efficiency.py
# ============================================================================= 
""" Test module for ostap/fitting/efficiency.py
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ROOT, random, math 
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import cpp, VE, dsID, Ostap
from   ostap.logger.utils       import rooSilent
from   ostap.fitting.efficiency import Efficiency1D
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_efficiency' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make
x           = ROOT.RooRealVar ( 'x',  'test' , 0 , 10 )
xmin , xmax = x.minmax()

acc = ROOT.RooCategory( 'cut','cut')
acc.defineType('accept',1)
acc.defineType('reject',0)
varset  = ROOT.RooArgSet  ( x , acc )
ds      = ROOT.RooDataSet ( dsID() , 'test data' ,  varset )

eff0       = Models.Monotonic_pdf ( 'E0' , xvar = x , power = 2 , increasing = True )
eff0.phis  = 3.1415/2 , 3.1415/2 
margin     = 1.25 
emax       = margin * eff0 ( x.getMax() ) 

for i in range ( 10000 ) :
    
    xv = random.uniform ( xmin , xmax )
    
    x.setVal ( xv )
    
    ev = random.uniform ( 0 , emax )
    
    if eff0( xv )  > ev : acc.setIndex(1)
    else                : acc.setIndex(0) 
    
    ds.add ( varset ) 

np     = 20
dx     = (xmax-xmin)/np 
points = [ dx * i for i in range ( np + 1 ) ]

# =============================================================================
# use some PDF to parameterize efficiciency
def test_pdf () :


    effPdf = Models.Monotonic_pdf ( 'P6' , xvar = x , power = 3 , increasing = True )

    maxe   = margin * effPdf ( xmax )
    
    s0     = min ( 1.0 / emax , 1.0 / maxe ) 
    scale  = ROOT.RooRealVar ( 'scaleX' , 'scaleX' , s0 , 0.2 * s0 , 5.0 * s0  )
    
    eff2   = Efficiency1D ( 'E2' , effPdf , cut = acc  , scale = scale )
    
    r2     = eff2.fitTo ( ds )
    f2     = eff2.draw  ( ds )
    
    print (r2)
    
    for p in points :
        e  = eff2 ( p , error = True )
        ev = e.value()
        e0 = eff0 ( p ) / emax  
        print (' Point/Eff %4.1f %s%% (%.2f%%)'   % ( p , (100*e).toString ( '(%5.2f+-%4.2f)' ) ,  e0 * 100 ) )


# =============================================================================
# use some functions  to parameterize efficiciency
def test_vars1 () :

    from ostap.fitting.roofuncs import BernsteinPoly as BP 
    
    f      = BP ( 'G' , xvar = x , power = 4 )
    f.pars = 0.2 , 0.2 , 0.2 , 0.2 
        
    eff2   = Efficiency1D ( 'E3' , f.fun , cut = acc  , xvar = x )
    
    r2     = eff2.fitTo ( ds )
    f2     = eff2.draw  ( ds )
    
    print (r2)
    
    for p in points :
        e  = eff2 ( p , error = True )
        ev = e.value()
        e0 = eff0 ( p ) / emax  
        print (' Point/Eff %4.1f %s%% (%.2f%%)'   % ( p , (100*e).toString ( '(%5.2f+-%4.2f)' ) ,  e0 * 100 ) )


# =============================================================================
# use some functions  to parameterize efficiciency
def test_vars2 () :
    
    from ostap.fitting.roofuncs import MonotonicPoly as MP 

    f      = MP ( 'G' , xvar = x , increasing = True , power = 4 )
    f.pars = 0.6 , 0.8 , -0.1 , -0.6
    f.a    = 0.06
    f.b    = 2.72
    f.a.release ()
    f.b.release ()

    eff2   = Efficiency1D ( 'E4' , f , cut = acc  , xvar = x )
    
    r2     = eff2.fitTo ( ds )
    f2     = eff2.draw  ( ds )
    
    print (r2)
    
    for p in points :
        e  = eff2 ( p , error = True )
        ev = e.value()
        e0 = eff0 ( p ) / emax  
        print (' Point/Eff %4.1f %s%% (%.2f%%)'   % ( p , (100*e).toString ( '(%5.2f+-%4.2f)' ) ,  e0 * 100 ) )



    
# =============================================================================
if '__main__' == __name__ :


    test_pdf    ()
    test_vars1  ()
    test_vars2  ()

# =============================================================================
##                                                                      The END 
# ============================================================================= 
