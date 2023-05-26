#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fixes_formula.py
# Test module for Ostap::FormulaVar and friends
# @see https://github.com/root-project/root/commit/a470a3d85e8b5c93917d3e84c39e9d5f0066da97
# ============================================================================= 
""" Test module for Ostap.FormulaVar & friends 
- see https://github.com/root-project/root/commit/a470a3d85e8b5c93917d3e84c39e9d5f0066da97
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# =============================================================================
from   ostap.core.meta_info         import root_version 
from   ostap.core.core              import Ostap
import ostap.fitting.roofit 
import ostap.fitting.roocollections 
from   ostap.fitting.variables      import make_formula
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fixes_formula' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

def test_formula() :

    a    = ROOT.RooRealVar     ( 'a' , '' , 0 , 100 )
    b    = ROOT.RooRealVar     ( 'b' , '' , 0 , 100 )
    c    = ROOT.RooRealVar     ( 'c' , '' , 0 , 100 )
    
    vlst = ROOT.RooArgList     ( a , b , c ) 
    
    f1   = Ostap.FormulaVar    ( 'f1'  , 'a+1'   , vlst , True  ) 
    f2   = Ostap.FormulaVar    ( 'f2'  , 'a+b'   , vlst , True  ) 
    f3   = Ostap.FormulaVar    ( 'f3'  , 'a+b+c' , vlst , True  ) 
    f4   = Ostap.FormulaVar    ( 'f4'  , 'c+a'   , vlst , True  ) 
    
    f11  = Ostap.FormulaVar    ( 'f11' , 'a+1'   , vlst , False ) 
    f21  = Ostap.FormulaVar    ( 'f21' , 'a+b'   , vlst , False ) 
    f31  = Ostap.FormulaVar    ( 'f31' , 'a+b+c' , vlst , False ) 
    f41  = Ostap.FormulaVar    ( 'f41' , 'c+a'   , vlst , False ) 
    
    f12  = make_formula        ( 'f12' , 'a+1'   , vlst ) 
    f22  = make_formula        ( 'f22' , 'a+b'   , vlst ) 
    f32  = make_formula        ( 'f32' , 'a+b+c' , vlst ) 
    f42  = make_formula        ( 'f42' , 'c+a'   , vlst ) 

    used1 = Ostap.usedVariables ( 'a+1' , vlst )    
    assert 1 == len ( used1 ) , \
           'Invalid "used1" stuff %s for ROOT %s' % ( used1 , root_version )
    
    used2 = Ostap.usedVariables ( 'a+b' , vlst )    
    assert 2 == len ( used2 ) , \
           'Invalid "used2" stuff %s for ROOT %s' % ( used2 , root_version )
    
    used3 = Ostap.usedVariables ( 'a+b+c' , vlst )    
    assert 3 == len ( used3 ) , \
           'Invalid "used3" stuff %s for ROOT %s' % ( used3 , root_version )

    
    
# =============================================================================
if '__main__' == __name__ :

    test_formula() 

# =============================================================================
##                                                                      The END 
# ============================================================================= 
