#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_efficiency2.py
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
from   ostap.utils.utils        import timing
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_efficiency2' )
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

def eff  ( x , a , b ,  c , x0 ) :
    return a + 0.5 * b * ( 1.0 + math.tanh ( c * 1.0 * ( x - x0 ) ) ) 

def eff0 ( x ) : return eff ( x , A , B , C , X0 ) 

emax       = 1.0  

for i in range ( 1000 ) :
    
    xv = random.uniform ( xmin , xmax )
    
    x.setVal ( xv )
    
    ev = random.uniform ( 0 , emax )
    
    if eff0( xv )  > ev : acc.setIndex(1)
    else                : acc.setIndex(0) 
    
    ds.add ( varset ) 

np     = 20
dx     = (xmax-xmin)/np 
points = [ dx * i for i in range ( np + 1 ) ]

a       = ROOT.RooRealVar    ( 'a'  , 'a'  , 0.09 , 0.0  , 0.2 )
b       = ROOT.RooRealVar    ( 'b'  , 'b'  , 0.50 , 0.1  , 0.9 )
c       = ROOT.RooRealVar    ( 'c'  , 'c'  , 0.05 , 0.01 , 5   )
x0      = ROOT.RooRealVar    ( 'x0' , 'x0' , 4    , 1    , 9   )

vars    = ROOT.RooArgList ( x , a , b , c , x0 ) 

#x0.fix ( X0 )
#c .fix ( C  )
#b .fix ( B  )
    
# =============================================================================
# use RooFormulaVar to parameterise efficiency:
def test_formula () :
    
    effFunc = ROOT.RooFormulaVar ( "effFunc" , "a+0.5*b*(1+tanh((x-x0)*c))" , vars )
    
    eff1 = Efficiency1D( 'E1' , effFunc , cut  = acc , xvar = x )
    r1 = eff1.fitTo ( ds )
    r1 = eff1.fitTo ( ds )
    f1 = eff1.draw  ( ds , nbins = 20 )
    print(r1)
    for p in points :
        print(' Point/Eff %4.1f %s%% %.2f%%'   % ( p , (100*eff1 ( p , error = True )).toString ( '(%5.2f+-%4.2f)' ) , 100*eff0(p) ))

# =============================================================================
def test_pyvar () :

    from ostap.core.meta_info import root_version_int 
    if root_version_int >= 62200 :
        logger.warning("test_pyvar: test is disabled for ROOT verison %s" % root_version_int )
        return 

    from ostap.fitting.pyvar import PyVAR
    class MyVar(PyVAR) :

        def  evaluate ( self ) :

            vlist = self.varlist

            _x  = float ( vlist[0] )
            _a  = float ( vlist[1] ) 
            _b  = float ( vlist[2] ) 
            _c  = float ( vlist[3] ) 
            _x0 = float ( vlist[4] ) 

            r = eff ( _x , _a , _b , _c , _x0 )
            
            return r 
        
    myEff2 = MyVar ( 'myEff2' , vars = vars , title = 'title' )
    
    eff2 = Efficiency1D( 'E2' , myEff2.var , cut  = acc , xvar = x )

    r2 = eff2.fitTo ( ds )
    r2 = eff2.fitTo ( ds )
    f2 = eff2.draw  ( ds , nbins = 20 )
    print(r2)
    for p in points :
        print(' Point/Eff %4.1f %s%% %.2f%%'   % ( p , (100*eff2 ( p , error = True )).toString ( '(%5.2f+-%4.2f)' ) , 100*eff0(p) ))

# =============================================================================
def test_pyvar2 () :

    from ostap.core.meta_info import root_version_int 
    if root_version_int >= 62200 :
        logger.warning("test_pyvar2: test is disabled for ROOT verison %s" % root_version_int )
        return 

    from ostap.fitting.pyvar import PyVAR2

    
    myEff3 = PyVAR2 ( name = 'myEff3' , vars = vars , function  = eff )
    
    eff3 = Efficiency1D( 'E3' , myEff3.var , cut  = acc , xvar = x )

    r2 = eff3.fitTo ( ds )
    r2 = eff3.fitTo ( ds )
    f2 = eff3.draw  ( ds , nbins = 20 )
    print(r2)
    for p in points :
        print(' Point/Eff %4.1f %s%% %.2f%%'   % ( p , (100*eff3 ( p , error = True )).toString ( '(%5.2f+-%4.2f)' ) , 100*eff0(p) ))

    
# =============================================================================
if '__main__' == __name__ :
    
    with timing ('RooFormulaVar ', logger ) : test_formula ()
    
    with timing ('PyVAR         ', logger ) : test_pyvar   ()
    with timing ('PyVAR2        ', logger ) : test_pyvar2  ()


# =============================================================================
##                                                                      The END 
# ============================================================================= 
