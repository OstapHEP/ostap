#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_efficiency.py
# Test module for ostap/fitting/efficiency.py
# ============================================================================= 
""" Test module for ostap/fitting/efficiency.py
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import cpp, VE, dsID, Ostap, rooSilent 
from   ostap.fitting.efficiency import Efficiency1D
from   ostap.fitting.variables  import FIXVAR 
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
import ROOT, random, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_efficiency' )
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

eff0       = Models.Monotonic_pdf ( 'E0' , xvar = x , power = 3 , increasing = True )
eff0.phis  = 3.1415/1 , 3.1415/2 , 3.1415/3  
margin     = 1.25 
emax       = margin * eff0 ( x.getMax() ) 

N = 20000

for i in range ( N ) :
    
    xv = random.uniform ( xmin , xmax )
    
    x.setVal ( xv )
    
    ev = random.uniform ( 0 , emax )
    
    if eff0( xv )  > ev : acc.setIndex(1)
    else                : acc.setIndex(0) 
    
    ds.add ( varset ) 

np     = 20
dx     = (xmax-xmin)/np 
points = [ dx * i for i in range ( np + 1 ) ]



funs = set () 
# =================================================================================
## make comparison table 
def make_table ( func , title , prefix = "# ") :

    rows = [ ( 'x' , 'fitted eff [%]' , 'true eff [%]' , 'delta [%]' ) ]
    for p in points :

        e1  = 100 * func (  p , error = True ) 
        e2  = 100 * eff0 ( p ) / emax 
        d   = e1 - e2    
        row = "%4.2f" % p , \
              "%s"    %  e1.toString ( '(%5.2f+-%4.2f)'  ) ,\
              "%.2f"  %  e2  ,\
              "%s"    %  d .toString ( '(%5.2f+-%4.2f)'  )
        
        rows.append ( row )
    from ostap.logger.table import table
    return table ( rows , title = title , prefix = prefix ) 


# =============================================================================
# use some PDF to parameterize efficienct
def test_pdf () :
    
    logger = getLogger ( 'test_pdf' )

    effPdf = Models.Monotonic_pdf ( 'P4' , xvar = x , power = 4 , increasing = True )

    maxe   = margin * effPdf ( xmax )
    
    s0     = min ( 1.0 / emax , 1.0 / maxe ) 
    scale  = ROOT.RooRealVar ( 'scaleX' , 'scaleX' , s0 , 0.2 * s0 , 5.0 * s0  )
    
    eff2   = Efficiency1D ( 'EPDF' , effPdf , cut = acc  , scale = scale )
    
    r2     = eff2.fitTo ( ds )
    
    logger.info ( "Fit result using-Monotonic_pdf \n%s" % r2.table ( prefix = "# ") )
    logger.info ( "Compare with true efficiency (using Monotonic_pdf)\n%s" % make_table (
        eff2 , title = 'using Monotonic_pdf') )
    
    with wait ( 2 ) , use_canvas ( 'test_pdf' ) : 
        f2     = eff2.draw  ( ds , nbins = 25 )


    funs.add ( effPdf )
    funs.add ( eff2   )
    

# =============================================================================
# use some functions  to parameterize efficiciency
def test_vars1 () :
    
    from ostap.fitting.roofuncs import BernsteinPoly as BP 
    
    logger = getLogger ( 'test_vars1' )

    power = 3 
    f     = BP ( 'B3'           ,
                 xvar  = x     ,
                 power = power ,
                 pars  = ( power + 1 ) * [ ( 0.2 , 0 , 1 ) ] )

    for p in f.pars :
        p.setMin ( 0 )
        p.setMax ( 1 )
        p.setVal ( 0.2)
        p.release()
        
    ## f.release_par() 
    
    eff2   = Efficiency1D ( 'Ev1' , f.fun , cut = acc  , xvar = x )
    
    r2     = eff2.fitTo ( ds )

    logger.info ( "Fit result using-BernsteinPoly \n%s" % r2.table ( prefix = "# ") )
    logger.info ( "Compare with true efficiency (using BernsteinPoly)\n%s" % make_table (
        eff2 , title = 'using BernsteinPoly') )
    
    with wait ( 2 ) , use_canvas ( 'test_vars1' ) : 
        f2     = eff2.draw  ( ds , nbins = 25 )

    funs.add ( f    ) 
    funs.add ( eff2 ) 

# =============================================================================
# use some functions  to parameterize efficiciency
def test_vars2 () :
    
    logger = getLogger ( 'test_vars2' )

    from ostap.fitting.roofuncs import MonotonicPoly as MP 

    f      = MP ( 'M4' , xvar = x , increasing = True , power = 4 )
    f.pars = 0.6 , 0.8 , -0.1 , -0.6
    f.a    = 0.06
    f.b    = 2.72
    f.a.release ()
    f.b.release ()

    eff2   = Efficiency1D ( 'Ev2' , f , cut = acc  , xvar = x )
    
    r2     = eff2.fitTo ( ds )
    
    logger.info ( "Fit result using-MonotonicPoly \n%s" % r2.table ( prefix = "# ") )
    logger.info ( "Compare with true efficiency (using MonotonicPoly)\n%s" % make_table (
        eff2 , title = 'using MonotonicPoly') )
    
    with wait ( 2 ) , use_canvas ( 'test_vars2' ) : 
        f2     = eff2.draw  ( ds , nbins = 25 )
        
    funs.add ( f    ) 
    funs.add ( eff2 ) 

# =============================================================================
# use some functions  to parameterize efficiciency
def test_vars3 () :
    
    logger = getLogger ( 'test_vars3' )

    from ostap.fitting.roofuncs import BSplineFun as BS

    for power in range ( 4 ) :
        
        f      = BS ( 'BS%s' % power  , xvar = x ,
                      knots = ( 0 , 3 , 7 , 10 ) ,
                      power = power ,
                      pars  = 12 * [ ( 0.2 , 0 , 1 ) ] )
        
        eff2   = Efficiency1D ( 'Ev3_%s' % power , f , cut = acc  , xvar = x )
        
        r2     = eff2.fitTo ( ds )
        
        logger.info ( "Fit result using BSpline(%d) \n%s" %  ( power , r2.table ( prefix = "# ") ) )
        logger.info ( "Compare with true efficiency (using BSpine(%d))\n%s" % ( power , make_table (
            eff2 , title = 'using MonotonicPoly') ) ) 
        
        with wait ( 2 ) , use_canvas ( 'test_var3_%s' % power  ) : 
            f2     = eff2.draw  ( ds , nbins = 25 )
        
        funs.add ( f      )
        funs.add ( eff2   )

# =============================================================================
# use some functions  to parameterize efficiciency
def test_vars4 () :
## if 1 < 2 :     
    logger = getLogger ( 'test_vars4' )

    a  = ROOT.RooRealVar  ( 'A', 'a' , 0.05  ,   0   , 1   )
    b  = ROOT.RooRealVar  ( 'B', 'b' , 0.02  , -0.05 , 0.1 )
    c  = ROOT.RooRealVar  ( 'C', 'c' , 0.005 ,   0   , 0.1 )

    import ostap.fitting.roofuncs as     R
    from   ostap.fitting.funbasic import Fun1D 
    X   = Fun1D ( x , xvar = x , name = 'X' )
    
    F      = a + b * X + c * X**2
    F      = a         + c * X**2
    
    eff2   = Efficiency1D ( 'E5' , F , cut = acc  , xvar = x )
    
    r2     = eff2.fitTo ( ds )

    logger.info ( "Fit result using-Fun1D \n%s" % r2.table ( prefix = "# ") )
    logger.info ( "Compare with true efficiency (using Fun1D)\n%s" % make_table (
        eff2 , title = 'using Fnu1D') )

    
    with wait ( 2 ) , use_canvas ( 'test_vars4' ) : 
        f2     = eff2.draw  ( ds , nbins = 25 )
    
    funs.add ( X    )
    funs.add ( F    )
    funs.add ( eff2 )


# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        db['x'        ] = x  
        db['varset'   ] = varset
        db['dataset'  ] = ds
        print ('before loop' )
        for m in funs :
            db['funs:'     + m.name ] = m
            if hasattr ( m , 'fun' ) :
                db['fun/F:' + m.name ] = m.fun
            if hasattr ( m , 'pdf' ) :
                db['pdf/F:' + m.name ] = m.pdf
            if hasattr ( m , 'eff_fun' ) :
                db['eff_fun/F:' + m.name ] = m.eff_fun
            if hasattr ( m , 'eff_pdf' ) :
                db['eff_pdf/F:' + m.name ] = m.eff_pdf
        db['funs'   ] = funs
        db.ls() 

# =============================================================================
if '__main__' == __name__ :
    
    with timing ("PDF"   , logger ) :  
        test_pdf   ()
        
    with timing ("Vars1" , logger ) :  
        test_vars1 ()
        
    with timing ("Vars2" , logger ) :        
        test_vars2 ()
        
    with timing ("Vars3" , logger ) :        
        test_vars3 ()

    with timing ("Vars4" , logger ) :        
       test_vars4 ()

    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()

    pass


# =============================================================================
##                                                                      The END 
# ============================================================================= 
