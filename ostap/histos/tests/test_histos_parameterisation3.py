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
from   ostap.histos.param       import legendre_sum, chebyshev_sum
from   ostap.core.core          import hID, fID 
from   ostap.utils.timing       import timing
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import batch_env, wait  
import ostap.histos.param
import ostap.histos.histos
import ostap.fitting.funcs
import ROOT, random, time
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
## set batch form environment 
batch_env ( logger )
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import scipy
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    scipy = None 
    
# =============================================================================

h1   = ROOT.TH1F ( hID () , 'decreasing convex ' , 100 , 0 , 1 ) ; h1.Sumw2 () 
h2   = ROOT.TH1F ( hID () , 'increasing convex ' , 100 , 0 , 1 ) ; h2.Sumw2 () 
h3   = ROOT.TH1F ( hID () , 'increasing concave' , 100 , 0 , 1 ) ; h3.Sumw2 () 
h4   = ROOT.TH1F ( hID () , 'decreasing concave' , 100 , 0 , 1 ) ; h4.Sumw2 () 
h5   = ROOT.TH1F ( hID () , 'symmetric  convex ' , 100 , 0 , 1 ) ; h5.Sumw2 () 
h6   = ROOT.TH1F ( hID () , 'symmetric  concave' , 100 , 0 , 1 ) ; h6.Sumw2 () 

f1   = ROOT.TF1  ( fID () , '(x-1)**2'         , 0 , 1 )
f2   = ROOT.TF1  ( fID () , 'x**2'             , 0 , 1 )
f3   = ROOT.TF1  ( fID () , '1-(x-1)**2'       , 0 , 1 )
f4   = ROOT.TF1  ( fID () , '1-x**2'           , 0 , 1 )
f5   = ROOT.TF1  ( fID () , '4*(x-0.5)**2'     , 0 , 1 )
f6   = ROOT.TF1  ( fID () , '1-4*(x-0.5)**2'   , 0 , 1 )

f_2  = ROOT.TF2  ( fID () , 'x*x+y*y'     , -1 , 1 , 0 , 2          )
f_3  = ROOT.TF3  ( fID () , 'x*x+y*y+z*z' , -1 , 1 , 0 , 2 , -1 , 2 )

h_2  = ROOT.TH2F ( hID () , '' , 50 , -1 , 1 , 50 , 0 , 2 )
h_3  = ROOT.TH3F ( hID () , '' , 20 , -1 , 1 , 20 , 0 , 2 , 20 , -1 , 2 ) 

h_2 += f_2
h_3 += f_3

entries = 1000000

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
    return "%.4e" % math.sqrt(dd/(d1*d2))

## make a quadratic difference between histogram and function 
def diff1 ( result, histo ) :

    _fun0  = result[4]
    _norm  = result[3]
    
    _fun1  = lambda x : _fun0(x)*_norm 
    _fun2  = lambda x : float(histo(x))
        
    return _diff2_ ( _fun1 , _fun2 , histo.xmin() , histo.xmax() )

    
# =============================================================================
def test_positive_pdf () :

    logger = getLogger("test_positive_pdf")
    with timing ( 'Positive [6]' , logger ) :
        params = [ h.pdf_positive ( 6 , silent = True , draw = True ) for h in histos ]
            
    for h , f in zip ( histos , params ) :
        with wait ( 1 ) ,  use_canvas ( 'test_positive_pdf: ' + h.GetTitle() )  :
            f.plot.draw () 
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )

# =============================================================================
def test_monotonic_pdf () :
    
    logger = getLogger("test_monotonic_pdf")
    with timing ( 'Monotonic[6]' , logger ) :
        params = [ h1.pdf_decreasing ( 6 , silent = True , draw = True ) , 
                   h2.pdf_increasing ( 6 , silent = True , draw = True ) ,
                   h3.pdf_increasing ( 6 , silent = True , draw = True ) , 
                   h4.pdf_decreasing ( 6 , silent = True , draw = True ) ]

    for h , f in zip ( histos , params ) :
        with wait ( 1 ) ,  use_canvas ( 'test_monotonic_pdf: ' + h.GetTitle() )  :
            f.plot.draw () 
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )

# =============================================================================
def test_convex_pdf () :
    
    logger = getLogger("test_convex_pdf")
    with timing ( 'Convex   [6]' , logger ) :
        params = [ h1.pdf_convex_decreasing  ( 6 , silent = True , draw = True ) , 
                   h2.pdf_convex_increasing  ( 6 , silent = True , draw = True ) , 
                   h3.pdf_concave_increasing ( 6 , silent = True , draw = True ) , 
                   h4.pdf_concave_decreasing ( 6 , silent = True , draw = True ) ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 1 ) ,  use_canvas ( 'test_convex_pdf: ' + h.GetTitle() )  :
            f.plot.draw () 
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )


# =============================================================================
def test_convexonly_pdf () :
    
    logger = getLogger("test_convexonly_pdf")
    with timing ( 'ConvexP  [4]' , logger ) :
        params = [ h1.pdf_convexpoly  ( 4 , silent = True , draw = True ) , 
                   h2.pdf_convexpoly  ( 4 , silent = True , draw = True ) ,
                   h3.pdf_concavepoly ( 4 , silent = True , draw = True ) ,
                   h4.pdf_concavepoly ( 4 , silent = True , draw = True ) ,
                   h5.pdf_convexpoly  ( 4 , silent = True , draw = True ) ,
                   h6.pdf_concavepoly ( 4 , silent = True , draw = True ) ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 1 ) ,  use_canvas ( 'test_convexonly_pdf: ' + h.GetTitle() )  :
            f.plot.draw () 
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )
        
# =============================================================================
def test_positive_spline_pdf () :
    
    logger = getLogger("test_positive_spline_pdf")
    with timing ('P-spline [3,2]' , logger ) :
        params = [ h.pdf_pSpline ( (5,2) , silent = True , draw = True ) for h in histos ] 

    for h , f in zip ( histos , params ) :
        with wait ( 1 ) ,  use_canvas ( 'test_positive_spline_pdf: ' + h.GetTitle() )  :
            f.plot.draw () 
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )
            
# =============================================================================
def test_monotonic_spline_pdf () :
    
    logger = getLogger("test_monotonic_spline_pdf")
    with timing ('M-spline [2,2]' , logger ) :
        params = [ h1.pdf_mSpline ( ( 5 , 2 , False ) , silent = True , draw = True ) , 
                   h2.pdf_mSpline ( ( 5 , 2 , True  ) , silent = True , draw = True ) , 
                   h3.pdf_mSpline ( ( 5 , 2 , True  ) , silent = True , draw = True ) , 
                   h4.pdf_mSpline ( ( 5 , 2 , False ) , silent = True , draw = True ) ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 1 ) ,  use_canvas ( 'test_monotonic_spline_pdf: ' + h.GetTitle() )  :
            f.plot.draw () 
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )
            
# =============================================================================
def test_convex_spline_pdf () :
    
    logger = getLogger("test_convex_spline_pdf")
    with timing ('C-spline [2,2]' , logger ) :

        params = [ h1.pdf_cSpline ( ( 5 , 2 , False , True  ) , silent = True , draw = True ) , 
                   h2.pdf_cSpline ( ( 5 , 2 , True  , True  ) , silent = True , draw = True ) , 
                   h3.pdf_cSpline ( ( 5 , 2 , True  , False ) , silent = True , draw = True ) , 
                   h4.pdf_cSpline ( ( 5 , 2 , False , False ) , silent = True , draw = True ) ] 
        
    for h , f in zip ( histos , params ) :
        with wait ( 1 ) ,  use_canvas ( 'test_convex_spline_pdf: ' + h.GetTitle() )  :
            f.plot.draw () 
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )

# =============================================================================
def test_convexonly_spline_pdf () :
    
    logger = getLogger("test_convex_spline_pdf")
    with timing ('C-spline [2,2]' , logger ) :
        params = [ h1.pdf_convexSpline  ( ( 5 , 2 ) , silent = True , draw = True ) ,
                   h2.pdf_convexSpline  ( ( 5 , 2 ) , silent = True , draw = True ) ,
                   h3.pdf_concaveSpline ( ( 5 , 2 ) , silent = True , draw = True ) ,
                   h4.pdf_concaveSpline ( ( 5 , 2 ) , silent = True , draw = True ) ,
                   h5.pdf_convexSpline  ( ( 5 , 2 ) , silent = True , draw = True ) ,
                   h6.pdf_concaveSpline ( ( 5 , 2 ) , silent = True , draw = True ) ]
        
    for h , f in zip ( histos , params ) :
        with wait ( 1 ) ,  use_canvas ( 'test_convexonly_spline_pdf: ' + h.GetTitle() )  :
            f.plot.draw () 
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )



# =============================================================================
def test_rational_pdf () :
    
    logger = getLogger("test_rational_pdf")
    p      = 2 
    with timing ('Rational [%s]' % p  , logger ) :
        for h in histos :
            with use_canvas ( 'test_rational_pdf: ' + h.GetTitle() , wait = 1 )  :
                h.draw() 
                for d in range ( 1 , p + 3 ) :
                    f = h.pdf_rational ( p , d , silent = True , draw = True )
                    f.plot.draw('same', color = d + 1 ) 
                    logger.info ( "%-25s : [%d] difference %s" %  ( h.title , d , diff1 ( f , h ) ) )

# =============================================================================
if '__main__' == __name__ :
    
    logger.info ( 100*'*')
    logger.info ( 'Parameterizations techniques using RooFit')
    logger.info ( 100*'*')
    
    test_positive_pdf          ()
    test_monotonic_pdf         ()
    test_convex_pdf            ()
    test_convexonly_pdf        ()
    
    test_rational_pdf          ()

    test_positive_spline_pdf   () 
    test_monotonic_spline_pdf  () 
    test_convex_spline_pdf     () 
    test_convexonly_spline_pdf () 


# =============================================================================
##                                                                      The END 
# =============================================================================
