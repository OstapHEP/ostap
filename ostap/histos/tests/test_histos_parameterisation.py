#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file ostap/histos/tests/test_histos_parameterisation.py
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
from   builtins import range
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'ostap.test_histos_parameterisation' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
import ostap.histos.param
import ostap.histos.histos
import ostap.fitting.funcs
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
from ostap.core.core    import hID , fID 
from ostap.utils.timing import timing

h1 = ROOT.TH1F ( hID () , 'decreasing convex ' , 100 , 0 , 1 ) ; h1.Sumw2 () 
h2 = ROOT.TH1F ( hID () , 'increasing convex ' , 100 , 0 , 1 ) ; h2.Sumw2 () 
h3 = ROOT.TH1F ( hID () , 'increasing concave' , 100 , 0 , 1 ) ; h3.Sumw2 () 
h4 = ROOT.TH1F ( hID () , 'decreasing concave' , 100 , 0 , 1 ) ; h4.Sumw2 () 
h5 = ROOT.TH1F ( hID () , 'symmetric  convex ' , 100 , 0 , 1 ) ; h5.Sumw2 () 
h6 = ROOT.TH1F ( hID () , 'symmetric  concave' , 100 , 0 , 1 ) ; h6.Sumw2 () 

f1 = ROOT.TF1  ( fID ()  , '(x-1)**2'       , 0 , 1 )
f2 = ROOT.TF1  ( fID ()  , 'x**2'           , 0 , 1 )
f3 = ROOT.TF1  ( fID ()  , '1-(x-1)**2'     , 0 , 1 )
f4 = ROOT.TF1  ( fID ()  , '1-x**2'         , 0 , 1 )
f5 = ROOT.TF1  ( fID ()  , '4*(x-0.5)**2'   , 0 , 1 )
f6 = ROOT.TF1  ( fID ()  , '1-4*(x-0.5)**2' , 0 , 1 )

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
    return "%.4e" % math.sqrt(dd/(d1*d2))

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

    logger = getLogger("test_bernstein_sum")

    with timing ( 'Bernstein-sum[6]' , logger ) :
        params  = [ h.bernstein_sum ( 6 ) for h in histos ]
    
    for h , f in zip ( histos , params ) :
        h.draw()
        f.draw('same')
        logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )
        time.sleep  (1) 

# =============================================================================
def test_bernsteineven_sum() :
            
    logger = getLogger("test_bernsteineven_sum")

    with timing ( 'Bernstein-(even)-sum[6]' , logger ) :
        params  = [ h.bernsteineven_sum ( 6 ) for h in histos ]
    
    for h , f in zip ( histos , params ) :
        h.draw()
        f.draw('same')
        logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )
        time.sleep  (1) 

# =============================================================================
def test_legendre_sum() :
    
    logger = getLogger("test_legendre_sum")

    with timing ( 'Legendre-sum[6]' , logger ) :
        params  = [ h.legendre_sum ( 6 ) for h in histos ]

    for h , f in zip ( histos , params ) :
        h.draw()
        f.draw('same')
        logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )
        time.sleep  (1) 
        
# =============================================================================
def test_chebyshev_sum() :
    
    logger = getLogger("test_chebyshev_sum")
    
    with timing ( 'Chebyshev-sum[6]' , logger ) :
        params  = [ h.chebyshev_sum ( 6 ) for h in histos ]

    for h , f in zip ( histos , params ) :
        h.draw()
        f.draw('same')
        logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )
        time.sleep  (1) 

# =============================================================================
def test_fourier_sum() :
            
    logger = getLogger("test_fourier_sum")

    with timing ( 'Fourier-sum[16]' , logger ) :
        params  = [ h.fourier_sum ( 16 ) for h in histos ]

    for h , f in zip ( histos , params ) :
        h.draw()
        f.draw('same')
        logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )
        time.sleep  (1) 

# =============================================================================
def test_cosine_sum() :

    logger = getLogger("test_cosine_sum")
    if not use_scipy :
        logger.warning("No scipy is avilable, skip the test test")
        return

    with timing ( 'Cosine-sum[16]' , logger ) :
        params  = [ h.cosine_sum ( 16 ) for h in histos ]
        
    for h , f in zip ( histos , params ) :
        h.draw()
        f.draw('same')
        logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )
        time.sleep  (1) 
        
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
    
# =============================================================================
##                                                                      The END 
# =============================================================================
