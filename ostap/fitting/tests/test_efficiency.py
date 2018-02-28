#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_efficiency.py
# Test module for ostap/fitting/efficiency.py
# ============================================================================= 
""" Test module for ostap/fitting/efficiency.py
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
    logger = getLogger ( 'test_efficiency' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make
x   = ROOT.RooRealVar ( 'x',  'test' , 0 , 1 )
acc = ROOT.RooCategory( 'cut','cut')
acc.defineType('accept',1)
acc.defineType('reject',0)
varset  = ROOT.RooArgSet ( x , acc )
ds      = ROOT.RooDataSet( dsID() , 'test data' ,  varset )

for i in range ( 10000 ) :
    v = random.gauss ( 3 , 2 ) / 10 
    ## v = random.uniform ( 0 , 1  ) 
    if not 0  <  v < 1 : continue 
    x.setVal ( v ) 
    acc.setIndex(0)
    ds.add ( varset )
    
for i in range ( 10000 ) :
    v = random.gauss ( 7 , 2 ) / 10
    ## v = random.uniform ( 0 , 1  ) 
    if not  0  <  v < 1 : continue
    x.setVal ( v ) 
    acc.setIndex(1)
    ds.add ( varset )

from ostap.fitting.efficiency import Efficiency1D

# =============================================================================
# use RooFormulaVar to parameterise efficiency:
def test_formula () :
    
    a       = ROOT.RooRealVar('a','a',0.0001,1.e-6,0.95)
    b       = ROOT.RooRealVar('b','b',  1,0,2)
    c       = ROOT.RooRealVar('c','c',  1,0.01,10)
    effFunc = ROOT.RooFormulaVar ("effFunc","(1-a)+a*cos((x-b)/c)",ROOT.RooArgList(x,a,b,c))
        
    eff1 = Efficiency1D( 'E1' , effFunc , cut  = acc , xvar = x )
    r1 = eff1.fitTo ( ds )
    f1 = eff1.draw  ( ds )
    
    
# =============================================================================
# use some PDF to parameterize efficiciency
def test_pdf () : 
    effPdf = Models.PolyPos_pdf ( 'B' , xvar = x , power = 6 )
    
    eff2 = Efficiency1D( 'E2' , effPdf , cut = acc  )
    r2 = eff2.fitTo ( ds )
    f2 = eff2.draw  (  ds )
    
# =============================================================================
if '__main__' == __name__ :
    
    test_formula ()
    test_pdf     ()
    

# =============================================================================
# The END 
# ============================================================================= 
