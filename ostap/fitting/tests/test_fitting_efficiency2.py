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
from   ostap.core.core          import cpp, VE, dsID
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

A  = 0.1
B  = 0.8
C  = 0.5
X0 = 6.0

def eff0 ( x ) :
    return A + 0.5 * B * ( 1.0 + math.tanh ( C * 1.0 * ( x - X0 ) ) ) 

emax       = 1.0  

for i in range ( 50000 ) :
    
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
# use RooFormulaVar to parameterise efficiency:
def test_formula () :
    
    a       = ROOT.RooRealVar    ( 'a'  , 'a'  , 0.05 , 0.0  , 0.2 )
    b       = ROOT.RooRealVar    ( 'b'  , 'b'  , 0.50 , 0.1  , 0.9 )
    c       = ROOT.RooRealVar    ( 'c'  , 'c'  , 0.05 , 0.01 , 5   )
    x0      = ROOT.RooRealVar    ( 'x0' , 'x0' , 4    , 1    , 9   )
    effFunc = ROOT.RooFormulaVar ( "effFunc" , "a+0.5*b*(1+tanh((x-x0)*c))" , ROOT.RooArgList ( x , a , b , c , x0 ) )
    
    eff1 = Efficiency1D( 'E1' , effFunc , cut  = acc , xvar = x )
    r1 = eff1.fitTo ( ds )
    r1 = eff1.fitTo ( ds )
    f1 = eff1.draw  ( ds )
    print(r1)
    for p in points :
        print(' Point/Eff %4.1f %s%% %.2f%%'   % ( p , (100*eff1 ( p , error = True )).toString ( '(%5.2f+-%4.2f)' ) , 100*eff0(p) ))

    
# =============================================================================
if '__main__' == __name__ :
    
    test_formula     ()



# =============================================================================
# The END 
# ============================================================================= 
