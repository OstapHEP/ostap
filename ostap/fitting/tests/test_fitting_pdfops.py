#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_pdfops.py
# Test module for ostap/fitting/pdf_ops.py
# - It tests basic FUN methods 
# ============================================================================= 
""" Test module for ostap/fitting/pdfops.py
- It tests basic PDF methods 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.utils.timing        import timing
from   ostap.plotting.canvas     import use_canvas
from   ostap.utils.utils         import wait, batch_env 
from   ostap.fitting.models      import Gauss_pdf 
from   ostap.fitting.resolution  import ResoGauss 
from   ostap.math.math_ve        import *
import ROOT, random, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_pdfops' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
#

x = ROOT.RooRealVar ( 'x' , 'x-variable' , 0  , 10 )
y = ROOT.RooRealVar ( 'y' , 'y-variable' , 10 , 20 )
z = ROOT.RooRealVar ( 'z' , 'z-variable' , 20 , 30 )

pdfs = set() 

# =============================================================================
def test_pdfops_1 ( ) :
    
    logger = getLogger ( 'test_pdfops_1' ) 

    gx = Gauss_pdf ( 'GX' , xvar = x , mean = (2 , 1, 3) , sigma = (2,1,3) )
    gy = Gauss_pdf ( 'GY' , xvar = y , mean = (15,14,16) , sigma = (2,1,3) )
    gz = Gauss_pdf ( 'GZ' , xvar = z , mean = (27,26,28) , sigma = (2,1,3) )

    pdfs.add ( gx ) 
    pdfs.add ( gy ) 
    pdfs.add ( gz ) 

    logger.info ( 'gx       : %s' % gx )
    logger.info ( 'gy       : %s' % gy )
    logger.info ( 'gz       : %s' % gz )
    
    with wait ( 1 ), use_canvas ( 'test_pdfops_1:%s' % gx.name ) :  gx.draw ()
    with wait ( 1 ), use_canvas ( 'test_pdfops_1:%s' % gy.name ) :  gx.draw ()
    with wait ( 1 ), use_canvas ( 'test_pdfops_1:%s' % gz.name ) :  gx.draw ()

    gxy = gx*gy
    gyz = gy*gz
    gxz = gx*gz
    
    pdfs.add ( gxy ) 
    pdfs.add ( gyz ) 
    pdfs.add ( gxz ) 

    logger.info ( 'gx*gy    : %s' % gxy )
    logger.info ( 'gy*gz    : %s' % gyz )
    logger.info ( 'gx*gz    : %s' % gxz )

    for g in ( gxy , gyz , gxz ) :

        with wait ( 1 ) , use_canvas ( 'test_pdfops_1:%s.draw1' % g.name ) :  g.draw1 ()
        with wait ( 1 ) , use_canvas ( 'test_pdfops_1:%s.draw2' % g.name ) :  g.draw2 ()
        
    gxyz = gx*gy*gz
    gxyy = gx*gy*gy
    gxyx = gx*gy*gx
    gxxx = gx*gx*gx

    logger.info ( 'gx*gy*gz : %s' % gxyz )
    logger.info ( 'gx*gy*gy : %s' % gxyy )
    logger.info ( 'gx*gy*gx : %s' % gxyx )
    logger.info ( 'gx*gx*gx : %s' % gxxx )

    
    with wait ( 1 ) , use_canvas ( 'test_pdfops_1:%s.draw1' % gxyz.name ) :  gxyz.draw1 ()
    with wait ( 1 ) , use_canvas ( 'test_pdfops_1:%s.draw2' % gxyz.name ) :  gxyz.draw2 ()
    with wait ( 1 ) , use_canvas ( 'test_pdfops_1:%s.draw3' % gxyz.name ) :  gxyz.draw3 ()

    with wait ( 1 ) , use_canvas ( 'test_pdfops_1:%s.draw1' % gxyy.name ) :  gxyy.draw1 ()
    with wait ( 1 ) , use_canvas ( 'test_pdfops_1:%s.draw2' % gxyy.name ) :  gxyy.draw2 ()

    with wait ( 1 ) , use_canvas ( 'test_pdfops_1:%s.draw1' % gxyx.name ) :  gxyx.draw1 ()
    with wait ( 1 ) , use_canvas ( 'test_pdfops_1:%s.draw2' % gxyx.name ) :  gxyx.draw2 ()

    with wait ( 1 ) , use_canvas ( 'test_pdfops_1:%s.draw'  % gxxx.name ) :  gxxx.draw  ()

    pdfs.add ( gxyz ) 
    pdfs.add ( gxyy ) 
    pdfs.add ( gxyx ) 
    pdfs.add ( gxxx ) 

# =============================================================================
def test_pdfops_2 ( ) :
    
    logger = getLogger ( 'test_pdfops_2' ) 

    ## narrow gaussian 
    g0 = Gauss_pdf ( 'G1' , xvar = x , mean = (5,4,6) , sigma = (0.2,0.1,3) )

    ## resoltuion object
    reso = ResoGauss ( 'R' , xvar = x , sigma = 1 )

    with wait ( 1 ) , use_canvas ( 'test_pdfops_2/0: original' ) :  g0.draw ()

    g1 = g0 % reso
    with wait ( 1 ) , use_canvas ( 'test_pdfops_2/1: %s' % g1.name ) :  g1.draw ()

    g2 = g0 % 0.5
    with wait ( 1 ) , use_canvas ( 'test_pdfops_2/2: %s' % g2.name ) :  g2.draw ()

    ## g3 = g0 % ( 0.5 , 0.1 , 1.0 ) 
    ## with wait ( 1 ) , use_canvas ( 'test_pdfops_2/3: %s' % g3.name ) :  g3.draw ()

    ## g4 = ( 0.5 , 0.1 , 1.0 ) % g0  
    ## with wait ( 1 ) , use_canvas ( 'test_pdfops_2/4: %s' % g4.name ) :  g4.draw ()

        
    pdfs.add ( g0   ) 
    pdfs.add ( reso ) 
    pdfs.add ( g1   ) 
    pdfs.add ( g2   ) 
    ## pdfs.add ( g3   ) 
    ## pdfs.add ( g4   ) 

    
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
        db['zvar' ] = z 
        for m in pdfs :
            db['model:' + m.name ] = m
            db['roo_tot:%s' % m.name ] = m.pdf
            for i,s in enumerate ( m.signals ) :
                db['roo_sig%d:%s' % ( i , m.name ) ] = s
            for i, b in enumerate ( m.backgrounds ) : 
                db['roo_bkg%d:%s' % ( i , m.name ) ] = s
            for a in m.alist1 : 
                db['cmp:%s/%s' % ( m.name , a.name ) ] = a
        db['pdfs' ] = pdfs
        db.ls() 

# =============================================================================
if '__main__' == __name__ :

    with timing ('test_pdfops_1'          , logger ) :
        test_pdfops_1          () 

    with timing ('test_pdfops_2'          , logger ) :
        test_pdfops_2          () 

    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()


         
# =============================================================================
##                                                                      The END 
# ============================================================================= 
