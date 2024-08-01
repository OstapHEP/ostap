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
from   builtins                 import range
import ostap.histos.param
import ostap.histos.histos
import ostap.fitting.funcs
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait
from   ostap.histos.param       import legendre_sum, chebyshev_sum
from   ostap.core.core          import hID , fID, SE, Ostap  
from   ostap.utils.timing       import timing
from   ostap.utils.progress_bar import progress_bar 
import ostap.logger.table       as     T 
import ROOT, random, time
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'ostap.test_histos_parameterisation' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( 'Test for histogram parameterisation')
# =============================================================================
use_scipy = False 
try :
    import numpy 
    import scipy
    use_scipy = True 
except ImportError :
    use_scipy = False 
    logger.warning ("Numpy/scipy-dependent test are disabled!")
    
# =============================================================================

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

# =============================================================================
## 2D&2D histograms
h2D = ROOT.TH2F ( hID() , '2D histogram' ,
                  40 , -10 , 10 ,
                  40 , -10 , 10 )
h3D = ROOT.TH3F ( hID() , '3D histogram' ,
                  20 , -10 , 10 ,
                  20 , -10 , 10 ,
                  20 , -10 , 10 )

## NMAX = 2000000
## NN   = 0
## while NN < NMAX :
    
##     x = random.gauss ( 0 , 6 )
##     if not -10 < x < 10 : continue
    
##     y = random.gauss ( 0 , 6 )
##     if not -10 < y < 10 : continue

##     z = random.gauss ( 0 , 6 )
##     if not -10 < z < 10 : continue

##     NN += 1
    
##     h2D.Fill ( x , y     )
##     h2D.Fill ( y , x     )    
##     h2D.Fill ( x , z     )
##     h2D.Fill ( z , x     )    
##     h2D.Fill ( y , z     )
##     h2D.Fill ( z , y     )
    
##     h3D.Fill ( x , y , z )    
##     h3D.Fill ( x , z , y )    
##     h3D.Fill ( y , x , z )
##     h3D.Fill ( y , z , x )
##     h3D.Fill ( z , x , y )
##     h3D.Fill ( z , y , x )


fun2D = lambda x,y   : ( 220.0 - x * x ) * ( 160.0 - y * y )
fun3D = lambda x,y,z : ( 220.0 - x * x ) * ( 160.0 - y * y ) * ( 120.0 - z * z ) 


h2D  += fun2D
h3D  += fun3D 
    
# =============================================================================
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

# =============================================================================
## make a quadratic difference between histogram and function 
def diff1 ( func , histo ) :

    _fun1  = lambda x : func(x)
    _fun2  = lambda x : float(histo(x))
        
    return _diff2_ ( _fun1 , _fun2 , histo.xmin() , histo.xmax() )

# =============================================================================
## make a quadratic difference between histogram and function 
def diff2 ( func , histo ) :

    _f     =  func[2]
    _n     =  float(func[-1])

    _fun1  = lambda x : _n*_f(x)
    _fun2  = lambda x : float(histo(x))
        
    return _diff2_ ( _fun1 , _fun2 , histo.xmin() , histo.xmax() )

# =============================================================================
## make a quadratic difference between histogram and function 
def diff3 ( func , histo ) :

    _f     =  func[2]
    _n     =  float(func[3])

    _fun1  = lambda x : _n*_f(x)
    _fun2  = lambda x : float(histo(x))
        
    return _diff2_ ( _fun1 , _fun2 , histo.xmin() , histo.xmax() )


# =============================================================================
def test_bernstein_sum_orig() :

    logger = getLogger("test_bernstein_sum_orig")

    with timing ( 'Bernstein-sum-orig[6]' , logger ) :
        params  = [ h.bernstein_sum_orig ( 6 ) for h in histos ]

    for h , f in zip ( histos , params ) :
        with use_canvas ( 'test_bernstein_sum_orig: ' + h.GetTitle() , wait = 1 )  :
            h.draw()
            f.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )


# =============================================================================
def test_bernstein_sum_fill() :

    logger = getLogger("test_bernstein_sum_fill")

    with timing ( 'Bernstein-sum-fill[6]' , logger ) :
        params  = [ h.bernstein_sum_fill ( 6 ) for h in histos ]

    for h , f in zip ( histos , params ) :
        with use_canvas ( 'test_bernstein_sum_fill: ' + h.GetTitle() , wait = 1 )  :
            h.draw()
            f.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )


# =============================================================================
def test_bernstein_sum() :

    logger = getLogger("test_bernstein_sum")

    with timing ( 'Bernstein-sum[6]' , logger ) :
        params  = [ h.bernstein_sum ( 6 ) for h in histos ]

    for h , f in zip ( histos , params ) :
        with use_canvas ( 'test_bernstein_sum: ' + h.GetTitle() , wait = 1 )  :
            h.draw()
            f.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )


# =============================================================================
def test_bernsteineven_sum() :
            
    logger = getLogger("test_bernsteineven_sum")

    with timing ( 'Bernstein-(even)-sum[6]' , logger ) :
        params  = [ h.bernsteineven_sum ( 6 ) for h in histos[4:] ]
    
    for h , f in zip ( histos[4:] , params ) :
        with use_canvas ( 'test_bernsteineven_sum: ' + h.GetTitle() , wait = 1 )  :
            h.draw()
            f.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )

# =============================================================================
def test_legendre_sum_orig () :
    
    logger = getLogger("test_legendre_sum_orig")

    with timing ( 'Legendre-sum-orig[6]' , logger ) :
        params  = [ h.legendre_sum_orig ( 6 ) for h in histos ]

    for h , f in zip ( histos , params ) :
        with use_canvas ( 'test_legendre_sum_orig: ' + h.GetTitle() , wait = 1 )  :
            h.draw()
            f.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )

# =============================================================================
def test_legendre_sum_fill () :
    
    logger = getLogger("test_legendre_sum_fill ")

    with timing ( 'Legendre-sum-fill[6]' , logger ) :
        params  = [ h.legendre_sum_fill ( 6 ) for h in histos ]

    for h , f in zip ( histos , params ) :
        with use_canvas ( 'test_legendre_sum_fill: ' + h.GetTitle() , wait = 1 )  :
            h.draw()
            f.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )

# =============================================================================
def test_legendre_sum() :
    
    logger = getLogger("test_legendre_sum")

    with timing ( 'Legendre-sum[6]' , logger ) :
        params  = [ h.legendre_sum ( 6 ) for h in histos ]

    for h , f in zip ( histos , params ) :
        with use_canvas ( 'test_legendre_sum: ' + h.GetTitle() , wait = 1 )  :
            h.draw()
            f.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )


# =============================================================================
def test_chebyshev_sum_orig () :
    
    logger = getLogger("test_chebyshev_sum_orig")
    
    with timing ( 'Chebyshev-sum-orig[6]' , logger ) :
        params  = [ h.chebyshev_sum_orig ( 8 ) for h in histos ]

    for h , f in zip ( histos , params ) :
        with use_canvas ( 'test_chebyshev_sum_orig: ' + h.GetTitle() , wait = 1 )  :
            h.draw()
            f.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )

# =============================================================================
def test_chebyshev_sum_fill () :
    
    logger = getLogger("test_chebyshev_sum_fill")
    
    with timing ( 'Chebyshev-sum-fill[6]' , logger ) :
        params  = [ h.chebyshev_sum_fill ( 8 ) for h in histos ]

    for h , f in zip ( histos , params ) :
        with use_canvas ( 'test_chebyshev_sum_fill: ' + h.GetTitle() , wait = 1 )  :
            h.draw()
            f.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )

# =============================================================================
def test_chebyshev_sum() :
    
    logger = getLogger("test_chebyshev_sum")
    
    with timing ( 'Chebyshev-sum[6]' , logger ) :
        params  = [ h.chebyshev_sum ( 8 ) for h in histos ]

    for h , f in zip ( histos , params ) :
        with use_canvas ( 'test_chebyshev_sum: ' + h.GetTitle() , wait = 1 )  :
            h.draw()
            f.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )

# =============================================================================
def test_rational_fun () :
    
    logger = getLogger("test_rational_fun")

    n  = 8 
    with timing ( 'Rational_fun[%s]' % n  , logger ) :
        for h in histos :
            with use_canvas ( 'test_rational_fun: ' + h.GetTitle() , wait = 1 )  :
                h.draw()
                for d in range ( 0 , n + 1 ) :
                    f = h.rational_fun ( n , d ) 
                    f.draw ( 'same' , color = d + 1 )
                    logger.info ( "%-25s : [%d] difference %s" %  ( h.title , d , diff1 ( f , h ) ) )
                    
# =============================================================================
def test_fourier_sum() :
    
    logger = getLogger("test_fourier_sum")
    if not use_scipy :
        logger.warning("No numpy/scipy is avilable, skip the test test")
        return

    hh = [ h for h in histos if hasattr ( h , 'fourier_sum' ) ]
    with timing ( 'Fourier-sum[10]' , logger ) :
        params  = [ h.fourier_sum ( 10 ) for h in hh ]

    for h , f in zip ( hh , params ) :
        with use_canvas ( 'test_fourier_sum: ' + h.GetTitle() , wait = 1 )  :
            h.draw()
            f.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )

# =============================================================================
def test_cosine_sum() :

    logger = getLogger("test_cosine_sum")
    if not use_scipy :
        logger.warning("No numpy/scipy is avilable, skip the test test")
        return

    hh = [ h for h in histos if hasattr ( h , 'cosine_sum' ) ]
    with timing ( 'Cosine-sum[10]' , logger ) :
        params  = [ h.cosine_sum ( 10 ) for h in hh ]
        
    for h , f in zip ( hh , params ) :
        with use_canvas ( 'test_cosine_sum: ' + h.GetTitle() , wait = 1 )  :
            h.draw()
            f.draw('same')
            logger.info ( "%-25s : difference %s" %  ( h.title , diff1 ( f , h ) ) )
                        
# =============================================================================
def test_legendre_2D () :

    logger = getLogger("test_legendre_2D")

    lsums = {}
    lsums[ ('*','*') ] = fun2D 
    
    N1    = 2
    N2    = 11
    with timing ( '2D Legendre parameterisation' , logger ) :                   
        for nx in range ( N1 , N2 ) :
            for ny in range ( nx , N2 ) :
                lsum = h2D.legendre ( nx , ny )
                lsums [ ( nx , ny ) ] = lsum 
            
    def rdif ( x , y ) :
        fx = float ( x )
        fy = float ( y )
        return ( fx - fy ) / ( abs ( fx ) + abs ( fy ) )

    scale = 1.e+2
    rows  = [ ('nx' , 'ny' , 'mean [%] ' , 'rms [%]' , 'min [%]' , 'max [%]')  ]
    
    for n in progress_bar ( lsums ) :
        lsum    = lsums [ n ]
        nx , ny = n
        cnt     = SE()

        for i in range ( 1000 ) :
            
            x  = random.uniform ( -10 , 10 )
            y  = random.uniform ( -10 , 10 )
            
            vl   = lsum ( x  , y  )
            vh   = h2D  ( x  , y  ) 
            cnt += rdif ( vl , vh ) * scale 

            vl   = lsum ( y  , x  )
            vh   = h2D  ( y  , x  ) 
            cnt += rdif ( vl , vh ) * scale 

        row = '%s' % nx , '%s' % ny , \
              '%+.3g' % ( float ( cnt.mean() ) ) , \
              '%.3g'  % ( float ( cnt.rms () ) ) , \
              '%+.3g' % ( float ( cnt.min () ) ) , \
              '%+.3g' % ( float ( cnt.max () ) ) 
        rows.append ( row )

    title = '2D Legendre parameterisations'
    table = T.table ( rows, title = title , prefix = '# ' , alignment = 'rrllrl' )
    logger.info ( '%s\n%s' % ( title , table ) ) 

# =============================================================================
def test_bernstein_2D () :

    logger = getLogger("test_bernstein_2D")

    bsums = {}
    bsums[ ('*','*') ] = fun2D 
    
    N1    = 2
    N2    = 11  
    with timing ( '2D Bernstein parameterisation' , logger ) :                   
        for nx in range ( N1 , N2 ) :
            for ny in range ( nx , N2 ) :
                bsum = h2D.bernstein ( nx , ny )
                bsums [ ( nx , ny ) ] = bsum 
                
    def rdif ( x , y ) :
        fx = float ( x )
        fy = float ( y )
        return ( fx - fy ) / ( abs ( fx ) + abs ( fy ) )

    scale = 1.e+2
    rows  = [ ('nx' , 'ny' , 'mean [%] ' , 'rms [%]' , 'min/max [%s]')  ]
    
    for n in progress_bar ( bsums ) :
        bsum    = bsums [ n ]
        nx , ny = n
        cnt     = SE()

        for i in range ( 1000 ) :
            
            x  = random.uniform ( -10 , 10 )
            y  = random.uniform ( -10 , 10 )
            
            vl   = bsum ( x  , y  )
            vh   = h2D  ( x  , y  ) 
            cnt += rdif ( vl , vh ) * scale 

            vl   = bsum ( y  , x  )
            vh   = h2D  ( y  , x  ) 
            cnt += rdif ( vl , vh ) * scale 

        row = '%s' % nx , '%s' % ny , \
              '%+.3g' % ( float ( cnt.mean() ) ) , \
              '%.3g'  % ( float ( cnt.rms () ) ) , \
              '%+.3g' % ( float ( cnt.min () ) ) , \
              '%+.3g' % ( float ( cnt.max () ) ) 
        rows.append ( row )
        
    title = '2D Bernstein parameterisations'
    table = T.table ( rows, title = title , prefix = '# ' , alignment = 'rrllrl' )
    logger.info ( '%s\n%s' % ( title , table ) ) 

# =============================================================================
def test_legendre_3D () :

    logger = getLogger("test_legendre_3D")

    lsums = {}

    lsums[ ('*','*','*') ] = fun3D 

    
    N1    = 2
    N2    = 7
    with timing ( '3D Legendre parameterisation' , logger ) :                   
        for nx in range ( N1 , N2 ) :
            for ny in range ( nx , N2 ) :
                for nz in range ( ny , N2 ) :
                    lsum = h3D.legendre ( nx , ny , nz )
                    lsums [ ( nx , ny , nz ) ] = lsum 
            
    def rdif ( x , y ) :
        fx = float ( x )
        fy = float ( y )
        return ( fx - fy ) / ( abs ( fx ) + abs ( fy ) )

    scale = 1.e+2
    rows  = [ ('nx' , 'ny' , 'nz' , 'mean [%] ' , 'rms [%]' , 'min [%s]' , 'max [%]')  ]
    
    for n in progress_bar ( lsums ) :
        lsum         = lsums [ n ]
        nx , ny , nz = n
        cnt          = SE()

        for i in range ( 1000 ) :
            
            x  = random.uniform ( -10 , 10 )
            y  = random.uniform ( -10 , 10 )
            z  = random.uniform ( -10 , 10 )
            
            vl   = lsum ( x  , y  , z )
            vh   = h3D  ( x  , y  , z ) 
            cnt += rdif ( vl , vh ) * scale 

            vl   = lsum ( y  , x  , z )
            vh   = h3D  ( y  , x  , z ) 
            cnt += rdif ( vl , vh ) * scale 

            vl   = lsum ( y  , z  , x )
            vh   = h3D  ( y  , z  , x ) 
            cnt += rdif ( vl , vh ) * scale 

            vl   = lsum ( z  , y  , x )
            vh   = h3D  ( z  , y  , x ) 
            cnt += rdif ( vl , vh ) * scale 

            vl   = lsum ( x  , z  , y )
            vh   = h3D  ( x  , z  , y ) 
            cnt += rdif ( vl , vh ) * scale 

            vl   = lsum ( z  , x  , y )
            vh   = h3D  ( z  , x  , y ) 
            cnt += rdif ( vl , vh ) * scale 


        row = '%s' % nx , '%s' % ny , '%s' % nz  , \
              '%+.3g' % ( float ( cnt.mean() ) ) , \
              '%.3g'  % ( float ( cnt.rms () ) ) , \
              '%+.3g' % ( float ( cnt.min () ) ) , \
              '%+.3g' % ( float ( cnt.max () ) ) 
        rows.append ( row )

    title = '3D Legendre parameterisations'
    table = T.table ( rows, title = title , prefix = '# ' , alignment = 'rrrllrl' )
    logger.info ( '%s\n%s' % ( title , table ) ) 


# =============================================================================
def test_bernstein_3D () :

    logger = getLogger("test_bernstein_3D")

    bsums = {}
    
    bsums[ ('*','*','*') ] = fun3D 

    N1    = 2
    N2    = 6
    with timing ( '3D Bernstein parameterisation' , logger ) :                   
        for nx in range ( N1 , N2 ) :
            for ny in range ( nx , N2 ) :
                for nz in range ( ny , N2 ) :
                    bsum = h3D.bernstein ( nx , ny , nz )
                    bsums [ ( nx , ny , nz ) ] = bsum 
                    
    def rdif ( x , y ) :
        fx = float ( x )
        fy = float ( y )
        return ( fx - fy ) / ( abs ( fx ) + abs ( fy ) )

    scale = 1.e+2
    rows  = [ ('nx' , 'ny' , 'nz' , 'mean [%] ' , 'rms [%]' , 'min [%s]' , 'max [%]')  ]
    
    for n in progress_bar ( bsums ) :
        bsum         = bsums [ n ]
        nx , ny , nz = n
        cnt          = SE()

        for i in range ( 1000 ) :
            
            x  = random.uniform ( -10 , 10 )
            y  = random.uniform ( -10 , 10 )
            z  = random.uniform ( -10 , 10 )
            
            vl   = bsum ( x  , y  , z )
            vh   = h3D  ( x  , y  , z ) 
            cnt += rdif ( vl , vh ) * scale 

            vl   = bsum ( y  , x  , z )
            vh   = h3D  ( y  , x  , z ) 
            cnt += rdif ( vl , vh ) * scale 

            vl   = bsum ( y  , z  , x )
            vh   = h3D  ( y  , z  , x ) 
            cnt += rdif ( vl , vh ) * scale 

            vl   = bsum ( z  , y  , x )
            vh   = h3D  ( z  , y  , x ) 
            cnt += rdif ( vl , vh ) * scale 

            vl   = bsum ( x  , z  , y )
            vh   = h3D  ( x  , z  , y ) 
            cnt += rdif ( vl , vh ) * scale 

            vl   = bsum ( z  , x  , y )
            vh   = h3D  ( z  , x  , y ) 
            cnt += rdif ( vl , vh ) * scale 


        row = '%s' % nx , '%s' % ny , '%s' % nz  , \
              '%+.3g' % ( float ( cnt.mean() ) ) , \
              '%.3g'  % ( float ( cnt.rms () ) ) , \
              '%+.3g' % ( float ( cnt.min () ) ) , \
              '%+.3g' % ( float ( cnt.max () ) ) 
        rows.append ( row )

    title = '3D Bernstein parameterisations'
    table = T.table ( rows, title = title , prefix = '# ' , alignment = 'rrrllrl' )
    logger.info ( '%s\n%s' % ( title , table ) ) 


# =============================================================================
if '__main__' == __name__ :
    
    logger.info ( 100*'*')
    logger.info ( 'Parameterizations techniques using histogram values only (in general fast)')
    logger.info ( 100*'*')

    test_bernstein_sum_orig ()
    test_bernstein_sum_fill ()
    test_bernstein_sum      ()

    test_bernsteineven_sum  ()
    
    test_legendre_sum_orig  ()
    test_legendre_sum_fill  ()
    test_legendre_sum       ()
    
    test_chebyshev_sum_orig ()
    test_chebyshev_sum_fill ()
    test_chebyshev_sum      ()
    
    test_rational_fun       ()

    test_fourier_sum        ()
    test_cosine_sum         ()

    test_legendre_2D        () 
    test_bernstein_2D       () 

    test_legendre_3D        () 
    test_bernstein_3D       () 
    """
    
# =============================================================================
##                                                                      The END 
# =============================================================================
