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
from   ostap.fitting.roofit     import FIXVAR 
from   ostap.core.core          import cpp, VE, dsID, Ostap, rooSilent 
from   ostap.fitting.efficiency import Efficiency1D
from   ostap.fitting.variables  import FIXVAR 
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait, batch_env 
from   ostap.logger.logger      import attention 
import ostap.fitting.models     as     Models 
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
## set batch from environment 
batch_env ( logger )
# =============================================================================

## make
x           = ROOT.RooRealVar ( 'x',  'test' , 0 , 10 )
xmin , xmax = x.minmax()

acc = ROOT.RooCategory( 'cut','cut')
acc.defineType('accept',1)
acc.defineType('reject',0)

eff0       = Models.Monotonic_pdf ( 'E0' , xvar = x , power = 3 , increasing = True )
eff0.phis  = 3.1415/1 , 3.1415/2 , 3.1415/3

margin     = 1.25 
emax       = margin * eff0 ( x.getMax() ) 

conf = { 'minimizer' : ('Minuit','migrad') , 'maxcalls' : 1000000 , 'refit' : 5 }


N = 20000

varset = ROOT.RooArgSet  ( x , acc )
ds     = ROOT.RooDataSet ( dsID() , 'test data' ,  varset )

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
        dd  = e1 - e2
        d   = dd .toString ( '(%5.2f+-%4.2f)' )
        
        if 5 < abs ( dd.value() ) : d = attention ( d )
            
        row = "%4.2f" % p , \
              "%s"    %  e1.toString ( '(%5.2f+-%4.2f)'  ) ,\
              "%.2f"  %  e2  ,\
              d
        
        rows.append ( row )
        
    from ostap.logger.table import table
    return table ( rows , title = title , prefix = prefix ) 


# =============================================================================
# use some PDF to parameterize efficienct
def test_eff_PDF () :
    
    logger = getLogger ( 'test_eff_PDF' )

    ## for power in range ( 3 , 4 ) :
    for power in range ( 1 , 2 ) :

        effPdf = Models.Monotonic_pdf ( 'M%d' % power , xvar = x , power = power , increasing = True )
        ## effPdf = Models.Bkg_pdf ( 'MGK' , xvar = x , power = 0 ) 
        
        maxe   = margin * effPdf ( xmax )
        
        s0     = 1.0 / effPdf ( 9.99 )
        
        ## smx    = min ( 1.0 / maxe , s0 * 100 ) 
        ## smn    = s0 / 100 
        scale  = ROOT.RooRealVar ( 'scale_M%d' % power , 'scaleX' , s0 , s0/100 , s0 * 2 ) 

        scale.setVal ( 0.40 / effPdf ( 4.0 ) ) 
        
        eff2   = Efficiency1D ( 'EPDF_M%d' % power , effPdf , cut = acc  , scale = scale )

        with FIXVAR ( effPdf.phis ) : 
            r2     = eff2.fitTo ( ds , **conf     )
            f2     = eff2.draw  ( ds , nbins = 25 )
            
        r2 = eff2.fitTo ( ds , **conf )
            
        title = 'Fit result using Monotonic_pdf[%d]' % power         
        logger.info ( "%s\n%s" % ( title , r2.table ( title = title , prefix = "# " ) ) )
        
        title = "using Monotonic_pdf[%d]" % power         
        logger.info ( "Compare with true efficiency (%s)\n%s" % ( title , make_table (
            eff2 , title = title ) ) ) 
        
        with wait ( 2 ) , use_canvas ( 'test_eff_PDF_M%d' % power) : 
            f2     = eff2.draw  ( ds , nbins = 25 )
            
        funs.add ( effPdf )
        funs.add ( eff2   )
        
# =============================================================================
# use some functions  to parameterize efficiciency
def test_eff_BP () :
    
    from ostap.fitting.roofuncs import BernsteinPoly as BP 
    
    logger = getLogger ( 'test_eff_BP' )

    for power in range ( 1 , 4 ) :
        
        f     = BP ( 'BP_%s' % power  ,
                     xvar  = x        ,
                     power = power    ,
                     pars  = ( power + 1 ) * [ ( 0.2 , 0 , 1 ) ] )
        
        for p in f.pars :
            p.setMin ( 0 )
        p.setMax ( 1 )
        p.setVal ( 0.2)
        p.release()
                
        eff2   = Efficiency1D ( 'Ef1_BP%d' % power , f.fun , cut = acc  , xvar = x )
        
        r2     = eff2.fitTo ( ds , refit = 5 )
        
        title = "Fit result using BernsteinPoly[%d]" % power 
        logger.info ( "%s\n%s" % ( title  , r2.table ( title = title , prefix = "# " ) ) )
        
        title = "using BernsteinPoly[%d]" % power 
        logger.info ( "Compare with true efficiency (%s)\n%s" % ( title , make_table (
            eff2 , title = title ) ) ) 
        
        with wait ( 2 ) , use_canvas ( 'test_eff_BP_PB%s' % power  ) : 
            f2     = eff2.draw  ( ds , nbins = 25 )
            
        funs.add ( f    ) 
        funs.add ( eff2 ) 
            
# =============================================================================
# use some functions  to parameterize efficiciency
def test_eff_MP () :
    
    logger = getLogger ( 'test_eff_MP' )

    from ostap.fitting.roofuncs import MonotonicPoly as MP 

    for power in ( 1 , 5 ) : 
        f      = MP ( 'MP%d' % power  , xvar = x , increasing = True , power = power )
        f.pars = 0.6 , 0.8 , -0.1 , -0.6
        f.a    = 0.06
        f.b    = 2.72
        f.a.release ()
        f.b.release ()
        
        eff2   = Efficiency1D ( 'Eff_MP%d'% power  , f , cut = acc  , xvar = x )
        
        r2     = eff2.fitTo ( ds , **conf )
        
        title = "Fit result using MonotonicPoly[%d]" % power 
        logger.info ( "%s\n%s" % ( title , r2.table ( title = title , prefix = "# " ) ) )
        title = "using MonotonicPoly[%d]" % power         
        logger.info ( "Compare with true efficiency (%s)\n%s" % ( title , make_table (
            eff2 , title = title ) ) ) 
        
        with wait ( 2 ) , use_canvas ( 'test_eff_MP_MP%d' % power ) : 
            f2     = eff2.draw  ( ds , nbins = 25 )
            
        funs.add ( f    ) 
        funs.add ( eff2 ) 
            
# =============================================================================
# use some functions  to parameterize efficiciency
def test_eff_BS () :
    
    logger = getLogger ( 'test_eff_BS' )

    from ostap.fitting.roofuncs import BSplineFun as BS

    for power in range ( 5 ) :
        
        f      = BS ( 'BSP_%s' % power  , xvar = x ,
                      knots = ( 0 , 3 , 7 , 10 ) ,
                      power = power ,
                      pars  = 12 * [ ( 0.2 , 0 , 1 ) ] )
        
        eff2   = Efficiency1D ( 'Eff_BSP_%s' % power , f , cut = acc  , xvar = x )
        
        r2     = eff2.fitTo ( ds , **conf )
        
        title = "Fit result using BSpline[%d]" % power 
        logger.info ( "%s \n%s" %  ( title, r2.table ( title = title , prefix = "# " ) ) )
        title = "using BSpline[%d]" % power         
        logger.info ( "Compare with true efficiency (%s)\n%s" % ( title , make_table (
            eff2 , title = title ) ) ) 
        
        with wait ( 2 ) , use_canvas ( 'test_eff_BS_%s' % power  ) : 
            f2     = eff2.draw  ( ds , nbins = 25 )
        
        funs.add ( f      )
        funs.add ( eff2   )


# =============================================================================
# use some functions  to parameterize efficiciency
def test_eff_RF () :
    
    logger = getLogger ( 'test_eff_RF' )

    from ostap.fitting.roofuncs import RationalFun as RF 

    n = 4
    
    for d in range ( 0 , n + 1 ) :
        
        f      = RF ( 'RF_%d%d'% ( n , d )  , xvar = x , n = n , d = d )

        for p in f.pars :
            p.setVal (  0.5  ) 
            p.setMin ( -0.01 )
            p.setMax (  1.01 )
        
        eff2   = Efficiency1D ( 'Eff_ER%d%d' % ( n , d ) , f.fun , cut = acc  , xvar = x )
        
        r2     = eff2.fitTo ( ds , **conf  )
        
        title = "Fit result using RationalFun[%d/%s]" % ( f.n , f.d )
        logger.info ( "%s\n%s" % ( title , r2.table ( title , prefix = "# " ) ) )
        
        title = "using RationalFun[%d/%s]" % ( f.n , f.d )        
        logger.info ( "Compare with true efficiency (%s)\n%s" % ( title , make_table (
            eff2 , title = title ) ) ) 
        
        with use_canvas ( 'test_eff_RF_[%s/%d]' % ( f.n , f.d )  , wait = 2 ) : 
            f2     = eff2.draw  ( ds , nbins = 25 )
            
        funs.add ( f    ) 
        funs.add ( eff2 )

# =============================================================================
# use some functions  to parameterize efficiciency
def test_eff_RB () :
    
    logger = getLogger ( 'test_eff_RB' )

    from ostap.fitting.roofuncs import RationalBernsteinFun as RB 

    p = 1
    
    for q in range ( 1 , 4 ) :
        
        f      = RB ( 'RB%d%d'% ( p , q )  , xvar = x , p = p , q = q )
        for v in f.pars [ : p + 1 ] :
            v.setVal (  0.02 ) 
            v.setMin (  0    )
            v.setMax (  100  )

                
        eff2   = Efficiency1D ( 'Eff_BR%d%d' % ( p , q ) , f.fun , cut = acc  , xvar = x )

        ## pre-fit with fixed denominator 
        with FIXVAR ( f.qpars ) : r2     = eff2.fitTo ( ds )

        ## fit 
        r2  = eff2.fitTo ( ds )
        
        title = "Fit result using RationalBernstein[%d/%d]" % ( f.p , f.q )
        logger.info ( "%s \n%s" % ( title , r2.table ( title = title , prefix = "# ") ) )
        
        title = "using RationalBernstein[%d/%d]" % ( f.p , f.q )
        logger.info ( "Compare with true efficiency (%s)\n%s" % ( title , make_table (
            eff2 , title = title ) ) ) 
        
        with use_canvas ( 'test_eff_RB_[%s/%d]' % ( p , q )  , wait = 2 ) : 
            f2     = eff2.draw  ( ds , nbins = 25 )
            
        funs.add ( f    ) 
        funs.add ( eff2 )


# =============================================================================
# use some functions  to parameterize efficiciency
def test_eff_FUN () :

    logger = getLogger ( 'test_eff_FUN' )

    a  = ROOT.RooRealVar  ( 'A', 'a' , 0.05  ,   0   , 1   )
    b  = ROOT.RooRealVar  ( 'B', 'b' , 0.02  , -0.05 , 0.1 )
    c  = ROOT.RooRealVar  ( 'C', 'c' , 0.005 ,   0   , 0.1 )

    import ostap.fitting.roofuncs as     R
    from   ostap.fitting.funbasic import Fun1D 
    X   = Fun1D ( x , xvar = x , name = 'X' )
    
    F      = a + b * X + c * X**2
    F      = a         + c * X**2
    
    eff2   = Efficiency1D ( 'E5' , F , cut = acc  , xvar = x )
    
    r2     = eff2.fitTo ( ds , **conf )

    logger.info ( "Fit result using-Fun1D \n%s" % r2.table ( prefix = "# ") )
    logger.info ( "Compare with true efficiency (using Fun1D)\n%s" % make_table (
        eff2 , title = 'using Fnu1D') )

    
    with wait ( 2 ) , use_canvas ( 'test_eff_FUN' ) : 
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

    """
    with timing ("test_eff_pdf"   , logger ) :  
        test_eff_PDF   ()

    """
    with timing ("test_eff_BP" , logger ) :  
        test_eff_BP    ()

    
    with timing ("test_eff_MP" , logger ) :        
        test_eff_MP  ()

    with timing ("test_eff_BS" , logger ) :        
        test_eff_BS ()

    with timing ("test_eff_RF" , logger ) :        
        test_eff_RF ()
        
    with timing ("test_eff_RB" , logger ) :        
        test_eff_RB()
        
    with timing ("test_eff_FUN" , logger ) :        
       test_eff_FUN ()
       
    
    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()


# =============================================================================
##                                                                      The END 
# ============================================================================= 
