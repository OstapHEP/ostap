#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file ostap/histos/tests/test_histos_parameterisation3.py
# Test module for ostap/histos/param.py
# - It tests parameterisation of histograms 
# ============================================================================= 
""" Test module for ostap/histos/param.py
- It tests parameterisations of histograms 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import 
# ============================================================================= 
import ROOT, random, ostap.histos.param, ostap.histos.histos, ostap.fitting.funcs
from   builtins import range
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'ostap.test_histos_parameterisation3' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( 'Test for histogram parameterisation')
# =============================================================================
try :
    import scipy
except ImportError :
    scipy = None 
    
# =============================================================================
from ostap.histos.param import legendre_sum, chebyshev_sum
from ostap.core.core    import hID, fID 
from ostap.utils.timing import timing

h1   = ROOT.TH1F ( hID() , 'histogram' , 100, 0 , 1 ) ; h1.Sumw2() 
h2   = ROOT.TH1F ( hID() , 'histogram' , 100, 0 , 1 ) ; h2.Sumw2() 
h3   = ROOT.TH1F ( hID() , 'histogram' , 100, 0 , 1 ) ; h3.Sumw2() 
h4   = ROOT.TH1F ( hID() , 'histogram' , 100, 0 , 1 ) ; h4.Sumw2() 
h5   = ROOT.TH1F ( hID() , 'histogram' , 100, 0 , 1 ) ; h5.Sumw2() 
h6   = ROOT.TH1F ( hID() , 'histogram' , 100, 0 , 1 ) ; h6.Sumw2() 

f1   = ROOT.TF1  ( fID() , '(x-1)**2'         , 0 , 1 )
f2   = ROOT.TF1  ( fID() , 'x**2'             , 0 , 1 )
f3   = ROOT.TF1  ( fID() , '1-(x-1)**2'       , 0 , 1 )
f4   = ROOT.TF1  ( fID() , '1-x**2'           , 0 , 1 )
f5   = ROOT.TF1  ( fID() , '4*(x-0.5)**2'     , 0 , 1 )
f6   = ROOT.TF1  ( fID() , '1-4*(x-0.5)**2'   , 0 , 1 )

f_2  = ROOT.TF2  ( fID() , 'x*x+y*y'     , -1 , 1 , 0 , 2          )
f_3  = ROOT.TF3  ( fID() , 'x*x+y*y+z*z' , -1 , 1 , 0 , 2 , -1 , 2 )

h_2  = ROOT.TH2F ( hID() , '' , 50 , -1 , 1 , 50 , 0 , 2 )
h_3  = ROOT.TH3F ( hID() , '' , 20 , -1 , 1 , 20 , 0 , 2 , 20 , -1 , 2 ) 

h_2 += f_2
h_3 += f_3

entries = 100000

## random.seed(10) 
for i in range ( 0 , entries ) :
    h1.Fill ( f1.GetRandom() )
    h2.Fill ( f2.GetRandom() )
    h3.Fill ( f3.GetRandom() )
    h4.Fill ( f4.GetRandom() )
    h5.Fill ( f5.GetRandom() )
    h6.Fill ( f6.GetRandom() )
    
# h1 - decreasing convex
# h2 - increasing convex
# h3 - increasing concave
# h4 - decreasing concave 
# h5 - non-monotonic convex     (symmetric)
# h6 - non-monotonic concave    (symmetric)


## make a quadratic difference between two functions 
def _diff2_ ( fun1 , fun2 , xmin , xmax ) :

    _fun1_  = lambda x : float(fun1(x))**2 
    _fun2_  = lambda x : float(fun2(x))**2 
    _fund_  = lambda x : (float(fun1(x))-float(fun2(x)))**2 
                              
    from ostap.math.integral import integral as _integral 
    import warnings
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        d1 = _integral ( _fun1_ , xmin , xmax )
        d2 = _integral ( _fun2_ , xmin , xmax )
        dd = _integral ( _fund_ , xmin , xmax )
        
    import math
    return math.sqrt(dd/(d1*d2))

## make a quadratic difference between histogram and function 
def diff1 ( result, histo ) :

    _fun0  = result[4]
    _norm  = result[3]
    
    _fun1  = lambda x : _fun0(x)*_norm 
    _fun2  = lambda x : float(histo(x))
        
    return _diff2_ ( _fun1 , _fun2 , histo.xmin() , histo.xmax() )

    
# =============================================================================
def test_positive_3 () :
    
    with timing ( 'Positive [4]' , logger ) :
        rB1 = h1.pdf_positive   ( 4 , silent = True , draw = True )
        rB2 = h2.pdf_positive   ( 4 , silent = True , draw = True )
        rB3 = h3.pdf_positive   ( 4 , silent = True , draw = True )
        rB4 = h4.pdf_positive   ( 4 , silent = True , draw = True )
        rB5 = h5.pdf_positive   ( 4 , silent = True , draw = True )
        rB6 = h6.pdf_positive   ( 4 , silent = True , draw = True )        
        logger.info ( 'Positive [4]: diff   %s ' %  [ diff1(*p) for p in [ (rB1 , h1) ,
                                                                           (rB2 , h2) ,
                                                                           (rB3 , h3) ,
                                                                           (rB4 , h4) ,
                                                                           (rB5 , h5) ,
                                                                           (rB6 , h6) ] ] )

# =============================================================================
def test_monotonic_3 () :
    
    with timing ( 'Monotonic[4]' , logger ) :
        rB1 = h1.pdf_decreasing ( 4 , silent = True , draw = True )
        rB2 = h2.pdf_increasing ( 4 , silent = True , draw = True )
        rB3 = h3.pdf_increasing ( 4 , silent = True , draw = True )
        rB4 = h4.pdf_decreasing ( 4 , silent = True , draw = True )
        logger.info ( 'Monotonic[4]: diff   %s ' %  [ diff1(*p) for p in [ (rB1 , h1) ,
                                                                           (rB2 , h2) ,
                                                                           (rB3 , h3) ,
                                                                           (rB4 , h4) ] ] )


# =============================================================================
def test_convex_3 () :
    
    with timing ( 'Convex   [4]' , logger ) :
        rB1 = h1.pdf_convex_decreasing  ( 4 , silent = True , draw = True )
        rB2 = h2.pdf_convex_increasing  ( 4 , silent = True , draw = True )
        rB3 = h3.pdf_concave_increasing ( 4 , silent = True , draw = True )
        rB4 = h4.pdf_concave_decreasing ( 4 , silent = True , draw = True )
        logger.info ( 'Convex   [4]: diff   %s ' %  [ diff1(*p) for p in [ (rB1 , h1) ,
                                                                           (rB2 , h2) ,
                                                                           (rB3 , h3) ,
                                                                           (rB4 , h4) ] ] )


# =============================================================================
def test_convexonly_3 () :
    
    with timing ( 'ConvexP  [4]' , logger ) :
        rB1 = h1.pdf_convexpoly  ( 4 , silent = True , draw = True )
        rB2 = h2.pdf_convexpoly  ( 4 , silent = True , draw = True )
        rB3 = h3.pdf_concavepoly ( 4 , silent = True , draw = True )
        rB4 = h4.pdf_concavepoly ( 4 , silent = True , draw = True )
        rB5 = h5.pdf_convexpoly  ( 4 , silent = True , draw = True )
        rB6 = h6.pdf_concavepoly ( 4 , silent = True , draw = True )
        
        logger.info ( 'ConvexP  [4]: diff   %s ' %  [ diff1(*p) for p in [ (rB1 , h1) ,
                                                                           (rB2 , h2) ,
                                                                           (rB3 , h3) ,
                                                                           (rB4 , h4) ,
                                                                           (rB5 , h5) ,
                                                                           (rB6 , h6) ] ] )
        
# =============================================================================
def test_positive_spline_3 () :
    
    with timing ('P-spline[3,2]' , logger ) :
        
        rC1 = h1.pdf_pSpline ( (3,2) , silent = True , draw = True )
        rC2 = h2.pdf_pSpline ( (3,2) , silent = True , draw = True )
        rC3 = h3.pdf_pSpline ( (3,2) , silent = True , draw = True )
        rC4 = h4.pdf_pSpline ( (3,2) , silent = True , draw = True )
        rC5 = h5.pdf_pSpline ( (3,2) , silent = True , draw = True )
        rC6 = h6.pdf_pSpline ( (3,2) , silent = True , draw = True )
        
        logger.info ( 'P-spline[3,2]: diff      %s ' %  [ diff1(*p) for p in [ (rC1 , h1) ,
                                                                               (rC2 , h2) ,
                                                                               (rC3 , h3) ,
                                                                               (rC4 , h4) ,
                                                                               (rC5 , h5) ,
                                                                               (rC6 , h6) ] ] )
        
# =============================================================================
def test_monotonic_spline_3 () :
    
    with timing ('M-spline[2,2]' , logger ) :
        
        rC1 = h1.pdf_mSpline ( ( 2 , 2 , False ) , silent = True , draw = True )
        rC2 = h2.pdf_mSpline ( ( 2 , 2 , True  ) , silent = True , draw = True )
        rC3 = h3.pdf_mSpline ( ( 2 , 2 , True  ) , silent = True , draw = True )
        rC4 = h4.pdf_mSpline ( ( 2 , 2 , False ) , silent = True , draw = True )
        
        logger.info ( 'M-spline[2,2]: diff      %s ' %  [ diff1(*p) for p in [ (rC1 , h1) ,
                                                                               (rC2 , h2) ,
                                                                               (rC3 , h3) ,
                                                                               (rC4 , h4) ] ] )

# =============================================================================
def test_convex_spline_3 () :
    
    with timing ('C-spline[2,2]' , logger ) :

        rC1 = h1.pdf_cSpline ( ( 2 , 2 , False , True  ) , silent = True , draw = True )
        rC2 = h2.pdf_cSpline ( ( 2 , 2 , True  , True  ) , silent = True , draw = True )
        rC3 = h3.pdf_cSpline ( ( 2 , 2 , True  , False ) , silent = True , draw = True )
        rC4 = h4.pdf_cSpline ( ( 2 , 2 , False , False ) , silent = True , draw = True )
        
        logger.info ( 'C-spline[2,2]: diff      %s ' %  [ diff1(*p) for p in [ (rC1 , h1) ,
                                                                               (rC2 , h2) ,
                                                                               (rC3 , h3) ,
                                                                               (rC4 , h4) ] ] )

# =============================================================================
def test_convexonly_spline_3 () :
    
    with timing ('C-spline[2,2]' , logger ) :

        rC1 = h1.pdf_convexSpline  ( ( 2 , 2 ) , silent = True , draw = True )
        rC2 = h2.pdf_convexSpline  ( ( 2 , 2 ) , silent = True , draw = True )
        rC3 = h3.pdf_concaveSpline ( ( 2 , 2 ) , silent = True , draw = True )
        rC4 = h4.pdf_concaveSpline ( ( 2 , 2 ) , silent = True , draw = True )
        rC5 = h5.pdf_convexSpline  ( ( 2 , 2 ) , silent = True , draw = True )
        rC6 = h6.pdf_concaveSpline ( ( 2 , 2 ) , silent = True , draw = True )
        
        logger.info ( 'C-spline[2,2]: diff      %s ' %  [ diff1(*p) for p in [ (rC1 , h1) ,
                                                                               (rC2 , h2) ,
                                                                               (rC3 , h3) ,
                                                                               (rC4 , h4) ,
                                                                               (rC5 , h5) ,
                                                                               (rC6 , h6) ] ] )


# =============================================================================
def test_convex_only_spline () :
    
    with timing ('Convex/Concave-spline[2,2]' , logger ) :
        
        rC1 = h1.convexSpline  ( degree = 2 , knots = 2 ) 
        rC2 = h2.convexSpline  ( degree = 2 , knots = 2 ) 
        rC3 = h3.concaveSpline ( degree = 2 , knots = 2 ) 
        rC4 = h4.concaveSpline ( degree = 2 , knots = 2 ) 
        rC5 = h1.convexSpline  ( degree = 2 , knots = 2 ) 
        rC6 = h1.concaveSpline ( degree = 2 , knots = 2 ) 
        
        logger.info ( 'Convex/Concave-spline[2,2]: diff      %s ' %  [ diff2(*p) for p in [ (rC1 , h1) ,
                                                                                (rC2 , h2) ,
                                                                                (rC3 , h3) ,
                                                                                (rC4 , h4) ,
                                                                                (rC5 , h5) ,
                                                                                (rC6 , h6) ] ] )


# =============================================================================
if '__main__' == __name__ :
    
    logger.info ( 100*'*')
    logger.info ( 'Parameterizations techniques using RooFit')
    logger.info ( 100*'*')
    
    test_positive_3          ()
    test_monotonic_3         ()
    test_convex_3            ()
    test_convexonly_3        ()

    test_positive_spline_3   () 
    test_monotonic_spline_3  () 
    test_convex_spline_3     () 
    test_convexonly_spline_3 () 


# =============================================================================
##                                                                      The END 
# =============================================================================
