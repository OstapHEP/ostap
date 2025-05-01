#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_roofuncs.py
# Test module for ostap/fitting/roofuncs.py
# - It tests basic FUN methods 
# ============================================================================= 
""" Test module for ostap/fitting/fnubasic.py
- It tests basic FUN methods 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.utils.timing       import timing
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.root_utils   import batch_env 
from   ostap.fitting.funbasic   import Id, Fun1D
from   ostap.math.math_ve       import *
import ostap.fitting.roofuncs   as     RF 
import ROOT, random, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_roofuncs' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
    
x = ROOT.RooRealVar ( 'x' , 'x-variable' , 1.e-2 , 5 )
y = ROOT.RooRealVar ( 'y' , 'y-variable' , 0     , 1 )

funs = set() 

# =============================================================================
## test basic operationsw with functions
def test_roofuncs_1 () :
    """Test basic operations with functions
    """
    
    logger = getLogger ( 'test_roofuncs_1' ) 
    logger.info  ( "Test basic operations with functions" ) 
    
    power = 6
    
    B  = RF.BernsteinPoly        ( 'Bernstein'   , xvar = x , power = power )
    I  = RF.MonotonicPoly        ( 'BernsteinI'  , xvar = x , power = power , increasing = True  )
    D  = RF.MonotonicPoly        ( 'DernsteinD'  , xvar = x , power = power , increasing = False )
    IX = RF.ConvexPoly           ( 'ConvexIX'    , xvar = x , power = power , increasing = True  , convex = True  )
    IV = RF.ConvexPoly           ( 'ConvexIV'    , xvar = x , power = power , increasing = True  , convex = False )
    DX = RF.ConvexPoly           ( 'ConvexDX'    , xvar = x , power = power , increasing = False , convex = True  )
    DV = RF.ConvexPoly           ( 'ConvexDV'    , xvar = x , power = power , increasing = False , convex = False )
    CX = RF.ConvexOnlyPoly       ( 'ConvexOnlyX' , xvar = x , power = power , convex = True  )
    CV = RF.ConvexOnlyPoly       ( 'ConvexOnlyV' , xvar = x , power = power , convex = False )
    RT = RF.RationalFun          ( 'Rational'    , xvar = x , n     = power , d = 2 )
    RB = RF.RationalBernsteinFun ( 'RationalBernstein'  , xvar = x , p = 3 , q = 3 )

    for p in ( B , I , D , IX , IV , DX , DV , CX , CV , RT , RB ) :
        p.pars = [ random.uniform ( -5 , 5 ) for i in range ( power + 2 ) ]
        with use_canvas ( 'test_roofuncs:%s' % p.name  , wait = 1 ) : p.draw ()
        funs.add  ( p )
        
# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db : 
        db['xvar' ] = x 
        db['yvar' ] = y 
        for f in funs :
            db['fun:' + f.name ] = f
            db['roo:' + f.name ] = f.fun
        db['funs' ] = funs
        db.ls() 

# =============================================================================
if '__main__' == __name__ :

    with timing ('test_roofuncs_1'          , logger ) :
        test_roofuncs_1          () 

    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()


         
# =============================================================================
##                                                                      The END 
# ============================================================================= 
