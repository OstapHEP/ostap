#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_funbasic.py
# Test module for ostap/fitting/funbasic.py
# - It tests basic FUN methods 
# ============================================================================= 
""" Test module for ostap/fitting/funbasic.py
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
import ROOT, random, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_funbasic' )
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
def test_funbasic_1 () :
    """Test basic operationsw with functions
    """
    
    logger = getLogger ( 'test_funbasic_1' ) 
    logger.info  ( "Test basic operationsw with functions" ) 

    Zero = Fun1D ( 0 , x , 'zero' )
    One  = Fun1D ( 1 , x , 'one'  )                                
    X    = Id   ( x )

    
    with use_canvas ( 'Zero' , wait = 1 ) :  Zero.draw ()
    with use_canvas ( 'One'  , wait = 1 ) :  One .draw ()
    with use_canvas ( 'X'    , wait = 1 ) :  X   .draw ()
        
    x2 = X*X
    
    with use_canvas ( 'x2'    , wait = 1 ) :  x2.draw ()

    pp = ( X - 2 ) ** 2 - 2 

    with use_canvas ( '(x-2)**2-2'  , wait = 1 ) :  pp.draw ()
    
    aa = abs ( sin ( x2 ) ) 
    with use_canvas ( 'abs(sin(x2))' , wait = 1 ) : aa.draw ()

    f1 = minv ( x2 ,  2.0 )
    with use_canvas ( 'min(x2,2)'    , wait = 1 ) : f1.draw ()

    f2 = maxv ( 1 - X ,  x2 )
    with use_canvas ( 'max(x2,1-x)'   , wait = 1 ) : f2.draw ()

    ix = 1/X
    with use_canvas ( '1/x'           , wait = 1 ) : ix.draw ()

    ex = 1-X
    with use_canvas ( '1-x'           , wait = 1 ) : ex.draw ()

    wx = 3+X
    with use_canvas ( '3+x'           , wait = 1 ) : wx.draw ()

    gf = exp ( -0.5 * ( X - 2.5) ** 2 )
    with use_canvas ( 'gauss'         , wait = 1 ) : gf.draw ()
    

    for ff in ( abs   ,
                exp   , erf    , erfc   , 
                sin   , cos    , tan    , atan , 
                sinh  , cosh   , tanh   ,
                log   , log10  , sech   , 
                gamma , lgamma , igamma ) :

        
        f = ff ( X )
        with use_canvas ( f.name  , wait = 1 ) : f.draw ()
        funs.add ( f ) 
                
    funs.add ( Zero )
    funs.add ( One  )
    funs.add ( X    )
    funs.add ( x2   )
    funs.add ( pp   )
    funs.add ( aa   )
    funs.add ( f1   )
    funs.add ( f2   )
    funs.add ( ix   )
    funs.add ( ex   )
    funs.add ( wx   )
    funs.add ( gf   )


# =============================================================================
## Test functionality from F1AUX mixin 
def test_funbasic_2 () :
    """Test functionality from F1AUX mixin 
    """
    
    logger = getLogger ( 'test_funbasic_2' ) 
    logger.info  ( "Test functionality from F1AUX mixin" )  

    
    Y = Id    ( y )
    
    mean  = 0.5
    sigma = 0.1

    G = exp   ( ( -0.5 / sigma**2 ) * ( Y - mean ) ** 2 ) / ( sigma * (2*math.pi)**0.5 )
    U = Fun1D ( 1 , y )
    
    with use_canvas ( 'Y' , wait = 1 ) : Y.draw ()
    with use_canvas ( 'U' , wait = 1 ) : U.draw ()
    with use_canvas ( 'G' , wait = 1 ) : G.draw ()
    

    header = 'Quantity' , 'uniform' , 'gauss'
    rows   =  [ header ]

    row    = 'get_mean'  , '%+.4g' % U.get_mean()  , '%+.4g' % G.get_mean () 
    rows.append ( row )

    row    = 'rms'       , '%+.4g' % U.rms()       , '%+.4g' % G.rms      () 
    rows.append ( row )

    row    = 'fwhm'      , ''                      , '%+.4g' % G.fwhm     () 
    rows.append ( row )

    row    = 'skewness'  , '%+.4g' % U.skewness()  , '%+.4g' % G.skewness () 
    rows.append ( row )

    row    = 'kurtosis'  , '%+.4g' % U.kurtosis()  , '%+.4g' % G.kurtosis () 
    rows.append ( row )

    row    = 'mode'      , ''                      , '%+.4g' % G.mode     () 
    rows.append ( row )

    row    = 'median'    , '%+.4g' % U.median  ()  , '%+.4g' % G.median   () 
    rows.append ( row )
    
    for i in range ( 7 ) : 
        row    = 'moment(%d)' % i , '%+.4g' % U.moment  ( i )  , '%+.4g' % G.moment ( i ) 
        rows.append ( row )

    for i in range ( 7 ) : 
        row    = 'central_moment(%d)' % i , \
                 '%+.4g' % U.central_moment ( i )  , \
                 '%+.4g' % G.central_moment ( i ) 
        rows.append ( row )

    for i in range ( 7 ) : 
        row    = 'std_moment(%d)' % i , \
                 '%+.4g' % U.std_moment ( i )  , \
                 '%+.4g' % G.std_moment ( i ) 
        rows.append ( row )
        
    row    = 'quantile(0.1)' , '%+.4g' % U.quantile(0.1) , '%+.4g' % G.quantile ( 0.1 )  
    rows.append ( row )

    row    = 'quantile(0.5)' , '%+.4g' % U.quantile(0.5) , '%+.4g' % G.quantile ( 0.5 )  
    rows.append ( row )

    row    = 'quantile(0.9)' , '%+.4g' % U.quantile(0.9) , '%+.4g' % G.quantile ( 0.9 )  
    rows.append ( row )


    for k in ( True , False ) :
        for i in range ( 7 ) : 
            row    = 'roo_moment(%d,central=%s)' %  (i, k)  , \
                     '%+.4g' % U.roo_moment ( i , central = k )  , \
                     '%+.4g' % G.roo_moment ( i , central = k ) 
            rows.append ( row )
            
    row    = 'roo_mean'  , '%+.4g' % U.roo_mean()  , '%+.4g' % G.roo_mean () 
    rows.append ( row )

    row    = 'roo_variance'  , '%+.4g' % U.roo_variance ()  , '%+.4g' % G.roo_variance  () 
    rows.append ( row )

    row    = 'roo_rms'      , '%+.4g' % U.roo_rms       ()  , '%+.4g' % G.roo_rms       () 
    rows.append ( row )

    row    = 'roo_skewness' , '%+.4g' % U.roo_skewness  ()  , '%+.4g' % G.roo_skewness  () 
    rows.append ( row )

    row    = 'roo_kurtosis' , '%+.4g' % U.roo_kurtosis  ()  , '%+.4g' % G.roo_kurtosis  () 
    rows.append ( row )

    for p in ( -0.2 , 0 , 0.2 ) :
        row    = 'derivative %g' % p , '%+.4g' % U.derivative ( p )  , '%+.4g' % G.derivative  ( p ) 
        rows.append ( row )

    xmin = 0.5
    xmax = 1.0 
    row  = 'integral(%g,%s)' %  ( xmin , xmax )  , \
           '%+.4g' % U.integral ( xmin , xmax )  , \
           '%+.4g' % G.integral ( xmin , xmax )  
    rows.append ( row )
    
        
    title = 'Test methdos from F1AUX moixin'
    import ostap.logger.table       as     T 
    table = T.table ( rows , title = title ,  prefix = '# ' , alignment = 'lll')
    logger.info ( '%s\n%s' % ( title , table ) ) 

    funs.add ( Y )
    funs.add ( U )
    funs.add ( G )
    
    
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

    with timing ('test_funbasic_1'          , logger ) :
        test_funbasic_1          () 

    with timing ('test_funbasic_2'          , logger ) :
        test_funbasic_2          () 

    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()


         
# =============================================================================
##                                                                      The END 
# ============================================================================= 
