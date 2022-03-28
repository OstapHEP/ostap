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
import ROOT, random, time
from   builtins              import range
# =============================================================================
import ostap.histos.param
import ostap.histos.histos
import ostap.fitting.funcs
from   ostap.plotting.canvas import use_canvas
from   ostap.utils.utils     import wait 
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
try :
    import numpy 
    import scipy
    use_scipy = True 
except ImportError :
    use_scipy = False 
    logger.warning ("Numpy/scipy-dependent are disables!")

# =============================================================================
from ostap.histos.param import legendre_sum, chebyshev_sum
from ostap.core.core    import hID, fID 
from ostap.utils.timing import timing

h1   = ROOT.TH1F ( hID () , 'decreasing convex ' , 100 , 0 , 1 ) ; h1.Sumw2 () 
h2   = ROOT.TH1F ( hID () , 'increasing convex ' , 100 , 0 , 1 ) ; h2.Sumw2 () 
h3   = ROOT.TH1F ( hID () , 'increasing concave' , 100 , 0 , 1 ) ; h3.Sumw2 () 
h4   = ROOT.TH1F ( hID () , 'decreasing concave' , 100 , 0 , 1 ) ; h4.Sumw2 () 
h5   = ROOT.TH1F ( hID () , 'symmetric  convex ' , 100 , 0 , 1 ) ; h5.Sumw2 () 
h6   = ROOT.TH1F ( hID () , 'symmetric  concave' , 100 , 0 , 1 ) ; h6.Sumw2 () 

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

## all histograms 
histos = h1 , h2 , h3 , h4 , h5 , h6

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
    return "%.4e" % math.sqrt(dd/math.sqrt(d1*d2)) 

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

    logger =   getLogger("test_bernstein")    
    with timing ( 'Bernstein [4]' , logger ) :
        params = [ h.bernstein ( 4 ) for h in  histos ]

    for h , f in zip ( histos , params ) :
        with wait ( 2 ) ,  use_canvas ( 'test_bernstein %s' % h.GetTitle()  ) : 
            h    .draw()
            f.tf1.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )
            
# =============================================================================
def test_chebyshev() :
    
    logger =   getLogger("test_chebyshev")    
    with timing ( 'Chebyshev [4]' , logger ) :
        params = [ h.chebyshev ( 4 ) for h in  histos ]

    for h , f in zip ( histos , params ) :
        with wait ( 2 ) ,  use_canvas ( 'test_chebyshev %s' % h.GetTitle()  ) : 
            h    .draw()
            f.tf1.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )
                
# =============================================================================
def test_legendre() :

    logger =   getLogger("test_legendre")    
    with timing ( 'Legendre [4]' , logger ) :
        params = [ h.legendre ( 4 ) for h in  histos ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 2 ) ,  use_canvas ( 'test_legendre %s' % h.GetTitle()  ) : 
            h    .draw()
            f.tf1.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )

# =============================================================================
def test_monomial() :

    logger =   getLogger("test_monomial")
    with timing ( 'Monomial [4]' , logger ) :
        params = [ h.polynomial ( 4 ) for h in  histos ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 2 ) ,  use_canvas ( 'test_monomial %s' % h.GetTitle()  ) : 
            h    .draw()
            f.tf1.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )

# =============================================================================
def test_positive() :

    logger =   getLogger("test_positive")
    with timing ( 'Positive [4]' , logger ) :
        params = [ h.positive ( 4 ) for h in  histos ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 2 ) ,  use_canvas ( 'test_positive %s' % h.GetTitle()  ) : 
            h    .draw()
            f.tf1.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )


# =============================================================================
def test_monotonic() :
    
    logger =   getLogger("test_monotonic")
    with timing ( 'Monotonic [4]' , logger ) :
        params = [ h1.monotonic ( 4 , increasing = False ) , 
                   h2.monotonic ( 4 , increasing = True  ) , 
                   h3.monotonic ( 4 , increasing = True  ) , 
                   h4.monotonic ( 4 , increasing = False ) ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 2 ) ,  use_canvas ( 'test_monotonic %s' % h.GetTitle()  ) : 
            h    .draw()
            f.tf1.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )
            
# =============================================================================
def test_convex () :
    
    logger =   getLogger("test_convex")
    with timing ( 'Convex [4]' , logger ) :
        params = [ h1.convex ( 4 , increasing = False , convex = True  ) , 
                   h2.convex ( 4 , increasing = True  , convex = True  ) , 
                   h3.convex ( 4 , increasing = True  , convex = False ) , 
                   h4.convex ( 4 , increasing = False , convex = False ) ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 2 ) ,  use_canvas ( 'test_convex %s' % h.GetTitle()  ) : 
            h    .draw()
            f.tf1.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )
            
# =============================================================================
def test_convex_poly () :

    logger =   getLogger("test_convex_poly")
    with timing ( 'Convex/ConcavePoly [4]' , logger ) :
        params = [ h1.convexpoly  ( 4 ) , 
                   h2.convexpoly  ( 4 ) ,
                   h3.concavepoly ( 4 ) ,
                   h4.concavepoly ( 4 ) ,
                   h5.convexpoly  ( 4 ) , 
                   h6.concavepoly ( 4 ) ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 2 ) ,  use_canvas ( 'test_convex_poly %s' % h.GetTitle()  ) : 
            h    .draw()
            f.tf1.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )

# =============================================================================
def test_fourier () :
    
    logger =   getLogger("test_fourier")
    
    if not use_scipy :
        logger.warning("No numpy/scipy is avilable, skip 'fourier' test")
        return

    with timing ( 'Fourier [4]' , logger ) :
        params = [ h5.fourier ( 4 ) , 
                   h6.fourier ( 6 ) ] 
        
    for h , f in zip  ( ( h5 , h6 ) , params ) :
        with wait ( 2 ) , use_canvas ( 'test_fourier: %s' % h.GetTitle () ) : 
            h    .draw  ()
            f.tf1.draw  ('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )

# =============================================================================
def test_cosine() :
    
    logger =   getLogger("test_cosine")
    if not use_scipy :
        logger.warning("No numpy/scipy is avilable, skip 'cosine' test")
        return
    
    with timing ( 'Cosine [4]' , logger ) :
        params = [ h.cosine ( 4 ) for h in  histos ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 2 ) ,  use_canvas ( 'test_cosine: %s' % h.GetTitle()  ) : 
            h    .draw()
            f.tf1.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )


# =============================================================================
def test_positive_spline () :
    
    logger =   getLogger("test_positive_spline")
    with timing ('P-spline [2,3]' , logger ) :
        params = [ h.pSpline ( degree = 2 , knots = 3 ) for h in  histos ]

    for h , f in zip ( histos , params ) :
        with wait ( 2 ) ,  use_canvas ( 'test_positive_spline: %s' % h.GetTitle()  ) : 
            h    .draw()
            f.tf1.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )
        
# =============================================================================
def test_monotonic_spline () :
    
    logger =   getLogger("test_monotonic_spline")
    with timing ('M-spline [2,3]' , logger ) :
        params = [ h1.mSpline ( degree = 2 , knots = 3 , increasing = False ) , 
                   h2.mSpline ( degree = 2 , knots = 3 , increasing = True  ) ,
                   h3.mSpline ( degree = 2 , knots = 3 , increasing = True  ) , 
                   h4.mSpline ( degree = 2 , knots = 3 , increasing = False ) ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 2 ) ,  use_canvas ( 'test_monotonic_spline %s' % h.GetTitle()  ) : 
            h    .draw()
            f.tf1.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )


# =============================================================================
def test_convex_spline () :
    
    logger =   getLogger("test_convex_spline")
    with timing ('C-spline [2,3]' , logger ) :
        params = [ h1.cSpline ( degree = 2 , knots = 3 , increasing = False , convex = True  ) , 
                   h2.cSpline ( degree = 2 , knots = 3 , increasing = True  , convex = True  ) ,
                   h3.cSpline ( degree = 2 , knots = 3 , increasing = True  , convex = False ) ,
                   h4.cSpline ( degree = 2 , knots = 3 , increasing = False , convex = False ) ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 2 ) ,  use_canvas ( 'test_convex_spline %s' % h.GetTitle()  ) : 
            h    .draw()
            f.tf1.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )


# =============================================================================
def test_convex_only_spline () :
    
    logger =   getLogger("test_convex_only_spline")
    with timing ('Convex/Concave-spline [2,3]' , logger ) :
        params = [ h1.convexSpline  ( degree = 2 , knots = 3 ) , 
                   h2.convexSpline  ( degree = 2 , knots = 3 ) , 
                   h3.concaveSpline ( degree = 2 , knots = 3 ) ,
                   h4.concaveSpline ( degree = 2 , knots = 3 ) ,
                   h5.convexSpline  ( degree = 2 , knots = 3 ) ,
                   h6.concaveSpline ( degree = 2 , knots = 3 ) ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 2 ) ,  use_canvas ( 'test_convex_only_spline %s' % h.GetTitle()  ) : 
            h    .draw()
            f.tf1.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff2 ( f , h ) ) )

# =============================================================================
#  Legendre fast
# =============================================================================
def test_legendre_fast () :
    
    logger =   getLogger("test_legenedre_fast")
    with timing ( 'Legendre-fast [6]' , logger ) :
        rL1 = h1.legendre_fast ( 6 )
        rL2 = h2.legendre_fast ( 6 )
        rL3 = h3.legendre_fast ( 6 )
        rL4 = h4.legendre_fast ( 6 )
        rL5 = h5.legendre_fast ( 6 )
        rL6 = h6.legendre_fast ( 6 )        
    logger.info ( 'difference %s ' %  [ diff1(*p) for p in [ (rL1 , h1) ,
                                                             (rL2 , h2) ,
                                                             (rL3 , h3) ,
                                                             (rL4 , h4) ,
                                                             (rL5 , h5) ,
                                                             (rL6 , h6) ] ] )
    
# =============================================================================
def test_legendre2_fast () :
    
    logger =   getLogger("test_legendre2_fast")
    with timing ( 'Legendre-2D' , logger ) :
        r22 = h_2.legendre_fast ( 2 , 2 )
        r32 = h_2.legendre_fast ( 3 , 2 )
        r42 = h_2.legendre_fast ( 4 , 2 )
        r52 = h_2.legendre_fast ( 5 , 2 )

        r33 = h_2.legendre_fast ( 3 , 3 )
        r43 = h_2.legendre_fast ( 4 , 3 )
        r53 = h_2.legendre_fast ( 5 , 3 )
        
    limits = h_2.xmin() , h_2.xmax() , h_2.ymin() , h_2.ymax()
    logger.info ( 'difference %s ' %  [ diff_2(f,p,*limits) for f,p in [ (r22 , f_2) ,
                                                                         (r32 , f_2) ,
                                                                         (r42 , f_2) ,
                                                                         (r52 , f_2) ,
                                                                         (r33 , f_2) ,
                                                                         (r43 , f_2) ,
                                                                         (r53 , f_2) ] ] )
    

# =============================================================================
def test_legendre3_fast () :
    
    logger =   getLogger("test_legendre3_fast")
    with timing ( 'Legendre-3D' , logger ) :
        
        r222 = h_3.legendre_fast ( 2 , 2 , 2 )  
        r333 = h_3.legendre_fast ( 3 , 3 , 3 )
        r432 = h_3.legendre_fast ( 4 , 3 , 2 ) 
        r535 = h_3.legendre_fast ( 5 , 3 , 5 )
        
    limits = h_3.xmin() , h_3.xmax() , h_3.ymin() , h_3.ymax() , h_3.zmin() , h_3.zmax()
    logger.info ( 'difference %s ' %  [ diff_3(f,p,*limits) for f,p in [ ( r222 , f_3) ,
                                                                         ( r333 , f_3) ,
                                                                         ( r432 , f_3) ,
                                                                         ( r535 , f_3) ] ] )
    
# =============================================================================
if '__main__' == __name__ :
    
    logger.info ( 100*'*')
    logger.info ( 'Parameterizations techniques using ROOT::TH1::Fit (could be slow)')
    logger.info ( 100*'*')
    
    test_bernstein              ()
    test_legendre               ()
    test_chebyshev              ()
    test_monomial               ()

    test_positive               ()
    test_monotonic              () 
    test_convex                 () 
    test_convex_poly            ()

    test_fourier                ()
    test_cosine                 ()
    
    test_positive_spline        ()
    test_monotonic_spline       ()
    test_convex_spline          ()
    test_convex_only_spline     ()
    

    test_legendre_fast          ()
    test_legendre2_fast         ()
    test_legendre3_fast         ()

# =============================================================================
##                                                                      The END 
# =============================================================================
