#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file ostap/histos/tests/test_histos_parameterisation2.py
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
    logger = getLogger ( 'ostap.test_histos_parameterisation2' )
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

f_2  = ROOT.TF2 ( fID() , 'x*x+y*y'     , -1 , 1 , 0 , 2          )
f_3  = ROOT.TF3 ( fID() , 'x*x+y*y+z*z' , -1 , 1 , 0 , 2 , -1 , 2 )

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
## qudratic difference between 2D-functions 
def diff_2 ( fun1 , fun2 , xmin , xmax , ymin , ymax , N = 100000 ) :

    d = 0
    for i in range ( N ) :
        x   = random.uniform ( xmin , xmax )
        y   = random.uniform ( ymin , ymax )
        f1  = fun1 ( x , y )
        f2  = fun2 ( x , y )
        d  += ( f1 - f2 ) * ( f1 - f2 ) 

    d /= N
    
    return d

# =============================================================================
## qudratic difference between 3D-functions 
def diff_3 ( fun1 , fun2 , xmin , xmax , ymin , ymax , zmin , zmax , N = 1000000 ) :

    d = 0
    for i in range ( N ) :
        x   = random.uniform ( xmin , xmax )
        y   = random.uniform ( ymin , ymax )
        z   = random.uniform ( zmin , zmax )
        f1  = fun1 ( x , y , z )
        f2  = fun2 ( x , y , z )
        d  +=  ( f1 - f2 ) * ( f1 - f2 ) 

    d /= N
    
    return d
    
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
def test_monotonic() :
    
    with timing ( 'Monotonic[4]' , logger ) :
        rL1 = h1.monotonic ( 4 , increasing = False )
        rL2 = h2.monotonic ( 4 , increasing = True  )
        rL3 = h3.monotonic ( 4 , increasing = True  )
        rL4 = h4.monotonic ( 4 , increasing = False )
        logger.info ( 'Monotonic[4]: diff    %s ' %  [ diff2(*p) for p in [ (rL1 , h1) ,
                                                                            (rL2 , h2) ,
                                                                            (rL3 , h3) ,
                                                                            (rL4 , h4) ] ] )

# =============================================================================
def test_convex () :
    
    with timing ( 'Convex[4]' , logger ) :
        rL1 = h1.convex ( 4 , increasing = False , convex = True  )
        rL2 = h2.convex ( 4 , increasing = True  , convex = True  )
        rL3 = h3.convex ( 4 , increasing = True  , convex = False )
        rL4 = h4.convex ( 4 , increasing = False , convex = False )
        logger.info ( 'Convex[4]: diff    %s ' %  [ diff2(*p) for p in [ (rL1 , h1) ,
                                                                            (rL2 , h2) ,
                                                                            (rL3 , h3) ,
                                                                            (rL4 , h4) ] ] )

# =============================================================================
def test_convex_poly () :
    
    with timing ( 'Convex/ConcavePoly[4]' , logger ) :
        rL1 = h1.convexpoly  ( 4 )
        rL2 = h2.convexpoly  ( 4 )
        rL3 = h3.concavepoly ( 4 )
        rL4 = h4.concavepoly ( 4 )
        rL5 = h5.convexpoly  ( 4 )
        rL6 = h6.concavepoly ( 4 )
        logger.info ( 'Convex/ConcavePoly[4]: diff    %s ' %  [ diff2(*p) for p in [ (rL1 , h1) ,
                                                                                     (rL2 , h2) ,
                                                                                     (rL3 , h3) ,
                                                                                     (rL4 , h4) ,
                                                                                     (rL5 , h5) ,
                                                                                     (rL6 , h6) ] ] )

# =============================================================================
def test_fourier () : 
    with timing ( 'Fourier[4]' , logger ) :
        # rF1 = h1.fourier    ( 8 )
        # rF2 = h2.fourier    ( 8 )
        # rF3 = h3.fourier    ( 8 )
        # rF4 = h4.fourier    ( 8 )
        rF5 = h5.fourier    ( 4 )
        rF6 = h6.fourier    ( 4 )
        logger.info ( 'Fourier[4] : diff    %s ' %  [ diff2(*p) for p in [ # (rF1 , h1) ,
                                                                           # (rF2 , h2) ,
                                                                           # (rF3 , h3) ,
                                                                           # (rF4 , h4) ,
                                                                           (rF5 , h5) ,
                                                                           (rF6 , h6) ] ] )
# =============================================================================
def test_cosine() :
    
    if not use_scipy :
        logger.warning("No scipy is avilable, skip 'cosine' test")
        return
    
    with timing ( 'Cosine[6]' , logger ) :
        rC1 = h1.cosine     ( 6 )
        rC2 = h2.cosine     ( 6 )
        rC3 = h3.cosine     ( 6 )
        rC4 = h4.cosine     ( 6 ) 
        rC5 = h5.cosine     ( 6 )
        rC6 = h6.cosine     ( 6 ) 
        logger.info ( 'Cosine[6]: diff      %s ' %  [ diff2(*p) for p in [ (rC1 , h1) ,
                                                                           (rC2 , h2) ,
                                                                           (rC3 , h3) ,
                                                                           (rC4 , h4) ,
                                                                           (rC5 , h5) ,
                                                                           (rC6 , h6) ] ] )
# =============================================================================
def test_generic_spline () :
    
    with timing ('B-spline[1,4]' , logger ) :
        
        rC1 = h1.bSpline ( degree = 1 , knots = 4 ) 
        rC2 = h2.bSpline ( degree = 1 , knots = 4 ) 
        rC3 = h3.bSpline ( degree = 1 , knots = 4 ) 
        rC4 = h4.bSpline ( degree = 1 , knots = 4 ) 
        rC5 = h5.bSpline ( degree = 1 , knots = 4 ) 
        rC6 = h6.bSpline ( degree = 1 , knots = 4 )
        
        logger.info ( 'B-spline[1,4]: diff      %s ' %  [ diff2(*p) for p in [ (rC1 , h1) ,
                                                                               (rC2 , h2) ,
                                                                               (rC3 , h3) ,
                                                                               (rC4 , h4) ,
                                                                               (rC5 , h5) ,
                                                                               (rC6 , h6) ] ] )
        
# =============================================================================
def test_positive_spline () :
    
    with timing ('P-spline[2,2]' , logger ) :
        
        rC1 = h1.pSpline ( degree = 2 , knots = 2 ) 
        rC2 = h2.pSpline ( degree = 2 , knots = 2 ) 
        rC3 = h3.pSpline ( degree = 2 , knots = 2 ) 
        rC4 = h4.pSpline ( degree = 2 , knots = 2 ) 
        rC5 = h5.pSpline ( degree = 2 , knots = 2 ) 
        rC6 = h6.pSpline ( degree = 2 , knots = 2 )
        
        logger.info ( 'P-spline[2,2]: diff      %s ' %  [ diff2(*p) for p in [ (rC1 , h1) ,
                                                                               (rC2 , h2) ,
                                                                               (rC3 , h3) ,
                                                                               (rC4 , h4) ,
                                                                               (rC5 , h5) ,
                                                                               (rC6 , h6) ] ] )
        
# =============================================================================
def test_monotonic_spline () :
    
    with timing ('M-spline[2,2]' , logger ) :
        
        rC1 = h1.mSpline ( degree = 2 , knots = 2 , increasing = False ) 
        rC2 = h2.mSpline ( degree = 2 , knots = 2 , increasing = True  ) 
        rC3 = h3.mSpline ( degree = 2 , knots = 2 , increasing = True  ) 
        rC4 = h4.mSpline ( degree = 2 , knots = 2 , increasing = False ) 
        
        logger.info ( 'M-spline[2,2]: diff      %s ' %  [ diff2(*p) for p in [ (rC1 , h1) ,
                                                                               (rC2 , h2) ,
                                                                               (rC3 , h3) ,
                                                                               (rC4 , h4) ] ] )


# =============================================================================
def test_convex_spline () :
    
    with timing ('C-spline[2,2]' , logger ) :
        
        rC1 = h1.cSpline ( degree = 2 , knots = 2 , increasing = False , convex = True  ) 
        rC2 = h2.cSpline ( degree = 2 , knots = 2 , increasing = True  , convex = True  ) 
        rC3 = h3.cSpline ( degree = 2 , knots = 2 , increasing = True  , convex = False ) 
        rC4 = h4.cSpline ( degree = 2 , knots = 2 , increasing = False , convex = False ) 
        
        logger.info ( 'C-spline[2,2]: diff      %s ' %  [ diff2(*p) for p in [ (rC1 , h1) ,
                                                                               (rC2 , h2) ,
                                                                               (rC3 , h3) ,
                                                                               (rC4 , h4) ] ] )


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
def test_legendre_fast () :
    
    with timing ( 'LegendreFast[6]' , logger ) :
        rL1 = h1.legendre_fast ( 6 )
        rL2 = h2.legendre_fast ( 6 )
        rL3 = h3.legendre_fast ( 6 )
        rL4 = h4.legendre_fast ( 6 )
        rL5 = h5.legendre_fast ( 6 )
        rL6 = h6.legendre_fast ( 6 )        
        logger.info ( 'LegendreFast[6]: diff %s ' %  [ diff1(*p) for p in [ (rL1 , h1) ,
                                                                            (rL2 , h2) ,
                                                                            (rL3 , h3) ,
                                                                            (rL4 , h4) ,
                                                                            (rL5 , h5) ,
                                                                            (rL6 , h6) ] ] )
        

# =============================================================================
def test_legendre2_fast () :
    
    with timing ( 'Legendre2D' , logger ) :
        
        r22 = h_2.legendre_fast ( 2 , 2 )
        r32 = h_2.legendre_fast ( 3 , 2 )
        r42 = h_2.legendre_fast ( 4 , 2 )
        r52 = h_2.legendre_fast ( 5 , 2 )

        r33 = h_2.legendre_fast ( 3 , 3 )
        r43 = h_2.legendre_fast ( 4 , 3 )
        r53 = h_2.legendre_fast ( 5 , 3 )
        
    limits = h_2.xmin() , h_2.xmax() , h_2.ymin() , h_2.ymax()
    logger.info ( 'Legendre2D: diff %s ' %  [ diff_2(f,p,*limits) for f,p in [ (r22 , f_2) ,
                                                                               (r32 , f_2) ,
                                                                               (r42 , f_2) ,
                                                                               (r52 , f_2) ,
                                                                               (r33 , f_2) ,
                                                                               (r43 , f_2) ,
                                                                               (r53 , f_2) ] ] )
    

# =============================================================================
def test_legendre3_fast () :
    
    with timing ( 'Legendre3D' , logger ) :
        
        r222 = h_3.legendre_fast ( 2 , 2 , 2 )  
        r333 = h_3.legendre_fast ( 3 , 3 , 3 )
        r432 = h_3.legendre_fast ( 4 , 3 , 2 ) 
        r535 = h_3.legendre_fast ( 5 , 3 , 5 )
        
    limits = h_3.xmin() , h_3.xmax() , h_3.ymin() , h_3.ymax() , h_3.zmin() , h_3.zmax()
    logger.info ( 'Legendre3D: diff %s ' %  [ diff_3(f,p,*limits) for f,p in [ ( r222 , f_3) ,
                                                                               ( r333 , f_3) ,
                                                                               ( r432 , f_3) ,
                                                                               ( r535 , f_3) ] ] )
    
# =============================================================================
if '__main__' == __name__ :
    
    logger.info ( 100*'*')
    logger.info ( 'Parameterizations techniques using ROOT::TH1::Fit (could be slow)')
    logger.info ( 100*'*')
    
    ## test_bernstein         ()
    ## test_legendre          ()
    ## test_chebyshev         ()
    ## test_monomial          ()

    ## test_monotonic         () 
    ## test_convex            () 
    ## test_convex_poly       ()
    
    ## test_fourier           ()
    ## test_cosine            ()

    ## test_generic_spline       () 
    ## test_positive_spline      () 
    ## test_monotonic_spline     () 
    ## test_convex_spline        () 
    ## test_convex_only_spline   () 
    
    test_legendre_fast       ()
    test_legendre2_fast      ()
    test_legendre3_fast      ()

# =============================================================================
# The END 
# =============================================================================
