#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file test_histo_parameterisation.py
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
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'ostap/histos/tests/test_histo_parameterisation' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( 'Test for histogram parameterisation')
# =============================================================================
use_scipy = False 
try :
    import scipy
    use_scipy = True 
except ImportError :
    use_scipy = False 
    
# =============================================================================
from ostap.histos.param import legendre_sum, chebyshev_sum
from ostap.core.core    import hID
from ostap.utils.timing import timing

h1 = ROOT.TH1F( hID() , 'histogram' , 100, 0 , 1 ) ; h1.Sumw2() 
h2 = ROOT.TH1F( hID() , 'histogram' , 100, 0 , 1 ) ; h2.Sumw2() 
h3 = ROOT.TH1F( hID() , 'histogram' , 100, 0 , 1 ) ; h3.Sumw2() 
h4 = ROOT.TH1F( hID() , 'histogram' , 100, 0 , 1 ) ; h4.Sumw2() 
h5 = ROOT.TH1F( hID() , 'histogram' , 100, 0 , 1 ) ; h5.Sumw2() 
h6 = ROOT.TH1F( hID() , 'histogram' , 100, 0 , 1 ) ; h6.Sumw2() 

f1 = ROOT.TF1('f3','(x-1)**2'         ,0,1)
f2 = ROOT.TF1('f4','x**2'             ,0,1)
f3 = ROOT.TF1('f3','1-(x-1)**2'       ,0,1)
f4 = ROOT.TF1('f4','1-x**2'           ,0,1)
f5 = ROOT.TF1('f5','4*(x-0.5)**2'     ,0,1)
f6 = ROOT.TF1('f5','1-4*(x-0.5)**2'   ,0,1)

entries = 100000

## random.seed(10) 
for i in xrange(0,entries) :
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
def diff1 ( func , histo ) :

    _fun1  = lambda x : func(x)
    _fun2  = lambda x : float(histo(x))
        
    return _diff2_ ( _fun1 , _fun2 , histo.xmin() , histo.xmax() )

## make a quadratic difference between histogram and function 
def diff2 ( func , histo ) :

    _f     =  func[2]
    _n     =  float(func[-1])

    _fun1  = lambda x : _n*_f(x)
    _fun2  = lambda x : float(histo(x))
        
    return _diff2_ ( _fun1 , _fun2 , histo.xmin() , histo.xmax() )

## make a quadratic difference between histogram and function 
def diff3 ( func , histo ) :

    _f     =  func[2]
    _n     =  float(func[3])

    _fun1  = lambda x : _n*_f(x)
    _fun2  = lambda x : float(histo(x))
        
    return _diff2_ ( _fun1 , _fun2 , histo.xmin() , histo.xmax() )

# =============================================================================
def test_bernstein_sum() :

    with timing ( 'Bezier[4]' ) :
        rB1 = h1.bernstein_sum ( 4 )
        rB2 = h2.bernstein_sum ( 4 )
        rB3 = h3.bernstein_sum ( 4 )
        rB4 = h4.bernstein_sum ( 4 )
        rB5 = h5.bernstein_sum ( 4 )
        rB6 = h6.bernstein_sum ( 4 )
        logger.info ( 'Bezier[4]: diff      %s ' %  [ diff1(*p) for p in [ (rB1 , h1) ,
                                                                           (rB2 , h2) ,
                                                                           (rB3 , h3) ,
                                                                           (rB4 , h4) ,
                                                                           (rB5 , h5) ,
                                                                           (rB6 , h6) ] ] )

# =============================================================================
def test_bernsteineven_sum() :
            
    
    with timing ( 'BezierEven[2]' , logger ) :
        rB1 = h1.bernsteineven_sum ( 2 )
        rB2 = h2.bernsteineven_sum ( 2 )
        rB3 = h3.bernsteineven_sum ( 2 )
        rB4 = h4.bernsteineven_sum ( 2 )
        rB5 = h5.bernsteineven_sum ( 2 )
        rB6 = h6.bernsteineven_sum ( 2 )
        logger.info ( 'BezierEven[2]: diff  %s ' %  [ diff1(*p) for p in [ (rB1 , h1) ,
                                                                           (rB2 , h2) ,
                                                                           (rB3 , h3) ,
                                                                           (rB4 , h4) ,
                                                                           (rB5 , h5) ,
                                                                           (rB6 , h6) ] ] )

# =============================================================================
def test_legendre_sum() :
    
    with timing ( 'Legendre[4]' , logger ) :
        rB1 = h1.legendre_sum ( 4 )
        rB2 = h2.legendre_sum ( 4 )
        rB3 = h3.legendre_sum ( 4 )
        rB4 = h4.legendre_sum ( 4 )
        rB5 = h5.legendre_sum ( 4 )
        rB6 = h6.legendre_sum ( 4 )
        logger.info ( 'Legendre[4]: diff    %s ' %  [ diff1(*p) for p in [ (rB1 , h1) ,
                                                                           (rB2 , h2) ,
                                                                           (rB3 , h3) ,
                                                                           (rB4 , h4) ,
                                                                           (rB5 , h5) ,
                                                                           (rB6 , h6) ] ] )
        
# =============================================================================
def test_chebyshev_sum() :
    
    with timing ( 'Chebyshev[4]' , logger ) :
        rB1 = h1.chebyshev_sum ( 4 )
        rB2 = h2.chebyshev_sum ( 4 )
        rB3 = h3.chebyshev_sum ( 4 )
        rB4 = h4.chebyshev_sum ( 4 )
        rB5 = h5.chebyshev_sum ( 4 )
        rB6 = h6.chebyshev_sum ( 4 )
        logger.info ( 'Chebyshev[4]: diff   %s ' %  [ diff1(*p) for p in [ (rB1 , h1) ,
                                                                           (rB2 , h2) ,
                                                                           (rB3 , h3) ,
                                                                           (rB4 , h4) ,
                                                                           (rB5 , h5) ,
                                                                           (rB6 , h6) ] ] )

# =============================================================================
def test_fourier_sum() :
    
    with timing ( 'Fourier[12]' , logger ) :
        rB1 = h1.fourier_sum ( 12 )
        rB2 = h2.fourier_sum ( 12 )
        rB3 = h3.fourier_sum ( 12 )
        rB4 = h4.fourier_sum ( 12 )
        rB5 = h5.fourier_sum ( 12 )
        rB6 = h6.fourier_sum ( 12 )
        logger.info ( 'Fourier[12]: diff    %s ' %  [ diff1(*p) for p in [ (rB1 , h1) ,
                                                                           (rB2 , h2) ,
                                                                           (rB3 , h3) ,
                                                                           (rB4 , h4) ,
                                                                           (rB5 , h5) ,
                                                                           (rB6 , h6) ] ] )

# =============================================================================
def test_cosine_sum() :
    
    if not use_scipy :
        logger.warning("No scipy is avilable, skip 'cosine_sum' test")
        return

    with timing ( 'Cosine[12]' , logger ) :
        rB1 = h1.cosine_sum ( 12 )
        rB2 = h2.cosine_sum ( 12 )
        rB3 = h3.cosine_sum ( 12 )
        rB4 = h4.cosine_sum ( 12 )
        rB5 = h5.cosine_sum ( 12 )
        rB6 = h6.cosine_sum ( 12 )
        logger.info ( 'Cosine[12]: diff     %s ' %  [ diff1(*p) for p in [ (rB1 , h1) ,
                                                                           (rB2 , h2) ,
                                                                           (rB3 , h3) ,
                                                                           (rB4 , h4) ,
                                                                           (rB5 , h5) ,
                                                                           (rB6 , h6) ] ] )

logger.info ( 100*'*')
logger.info ( 'Parameterizations techniques using ROOT::TH1::Fit (in general slow)')
logger.info ( 100*'*')

# =============================================================================
def test_bernstein() :
    
    with timing ( 'Bernstein[4]' , logger ) :
        rB1 = h1.bernstein      ( 4 )
        rB2 = h2.bernstein      ( 4 )
        rB3 = h3.bernstein      ( 4 )
        rB4 = h4.bernstein      ( 4 )
        rB5 = h5.bernstein      ( 4 )
        rB6 = h6.bernstein      ( 4 )        
        logger.info ( 'Bernstein[4]: diff   %s ' %  [ diff2(*p) for p in [ (rB1 , h1) ,
                                                                           (rB2 , h2) ,
                                                                           (rB3 , h3) ,
                                                                           (rB4 , h4) ,
                                                                           (rB5 , h5) ,
                                                                           (rB6 , h6) ] ] )

# =============================================================================
def test_chebyshev() :
    
    with timing ( 'Chebyshev[4]' , logger ) : 
        rC1 = h1.chebyshev  ( 4 )
        rC2 = h2.chebyshev  ( 4 )
        rC3 = h3.chebyshev  ( 4 )
        rC4 = h4.chebyshev  ( 4 )
        rC5 = h5.chebyshev  ( 4 )
        rC6 = h6.chebyshev  ( 4 )
        logger.info ( 'Chebyshev[4]: diff   %s ' %  [ diff2(*p) for p in [ (rC1 , h1) ,
                                                                           (rC2 , h2) ,
                                                                           (rC3 , h3) ,
                                                                           (rC4 , h4) ,
                                                                           (rC5 , h5) ,
                                                                           (rC6 , h6) ] ] )

# =============================================================================
def test_legendre() :
    
    with timing ( 'Legendre[4]' , logger ) :
        rL1 = h1.legendre   ( 4 )
        rL2 = h2.legendre   ( 4 )
        rL3 = h3.legendre   ( 4 )
        rL4 = h4.legendre   ( 4 )
        rL5 = h5.legendre   ( 4 )
        rL6 = h6.legendre   ( 4 )        
        logger.info ( 'Legendre[4]: diff    %s ' %  [ diff2(*p) for p in [ (rL1 , h1) ,
                                                                           (rL2 , h2) ,
                                                                           (rL3 , h3) ,
                                                                           (rL4 , h4) ,
                                                                           (rL5 , h5) ,
                                                                           (rL6 , h6) ] ] )

# =============================================================================
def test_monomial() :
    
    with timing ( 'Monomial[4]' , logger ) : 
        rP1 = h1.polynomial ( 4 )
        rP2 = h2.polynomial ( 4 )
        rP3 = h3.polynomial ( 4 )
        rP4 = h4.polynomial ( 4 )
        rP5 = h5.polynomial ( 4 )
        rP6 = h6.polynomial ( 4 )
        logger.info ( 'Monomial[4]: diff    %s ' %  [ diff2(*p) for p in [ (rP1 , h1) ,
                                                                           (rP2 , h2) ,
                                                                           (rP3 , h3) ,
                                                                           (rP4 , h4) ,
                                                                           (rP5 , h5) ,
                                                                           (rP6 , h6) ] ] )

# =============================================================================
def test_fourier () : 
    with timing ( 'Fourier[8]' , logger ) :
        rF1 = h1.fourier    ( 8 )
        rF2 = h2.fourier    ( 8 )
        rF3 = h3.fourier    ( 8 )
        rF4 = h4.fourier    ( 8 )
        rF5 = h5.fourier    ( 8 )
        rF6 = h6.fourier    ( 8 )
        logger.info ( 'Fourier[8] : diff    %s ' %  [ diff2(*p) for p in [ (rF1 , h1) ,
                                                                           (rF2 , h2) ,
                                                                           (rF3 , h3) ,
                                                                           (rF4 , h4) ,
                                                                           (rF5 , h5) ,
                                                                           (rF6 , h6) ] ] )
# =============================================================================
def test_cosine() :
    
    if not use_scipy :
        logger.warning("No scipy is avilable, skip 'cosine' test")
        return
    
    with timing ( 'Cosine[8]' , logger ) :
        rC1 = h1.cosine     ( 8 )
        rC2 = h2.cosine     ( 8 )
        rC3 = h3.cosine     ( 8 )
        rC4 = h4.cosine     ( 8 ) 
        rC5 = h5.cosine     ( 8 )
        rC6 = h6.cosine     ( 8 ) 
        logger.info ( 'Cosine[8]: diff      %s ' %  [ diff2(*p) for p in [ (rC1 , h1) ,
                                                                           (rC2 , h2) ,
                                                                           (rC3 , h3) ,
                                                                           (rC4 , h4) ,
                                                                           (rC5 , h5) ,
                                                                           (rC6 , h6) ] ] )

# =============================================================================
if '__main__' == __name__ :
    
    logger.info ( 100*'*')
    logger.info ( 'Parameterizations techniques using histogram values only (in general fast)')
    logger.info ( 100*'*')
    
    test_bernstein_sum     ()
    test_bernsteineven_sum ()
    test_legendre_sum      ()
    test_chebyshev_sum     ()
    test_fourier_sum       ()
    test_cosine_sum        ()

    logger.info ( 100*'*')
    logger.info ( 'Parameterizations techniques using ROOT::TH1::Fit (could be slow)')
    logger.info ( 100*'*')
    
    test_bernstein         ()
    test_legendre          ()
    test_chebyshev         ()
    test_monomial          ()
    test_fourier           ()
    test_cosine            ()
    
# =============================================================================
# The END 
# =============================================================================
