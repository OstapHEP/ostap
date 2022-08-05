#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_constraints.py
# Test module for consraubned fits 
# ============================================================================= 
""" Test module for constrained fits 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import VE
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait
from   ostap.utils.timing       import timing
import ROOT, random, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_constraints' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make
x           = ROOT.RooRealVar ( 'x',  'test' , 0 , 10 )
signal      = Models.Gauss_pdf ( 'Gauss' ,
                                 xvar  = x                ,
                                 mean  = ( 5 , 2    , 8 ) ,
                                 sigma = ( 1 , 0.01 , 5 ) ) 
background  = Models.Bkg_pdf   ( 'Bkg' , xvar = x , power = 1 )
background.tau.fix(0)

model       = Models.Fit1D ( signal = signal , background = background )
model.S     = 100
model.B     = 1000

data        = model.generate ( 4*1100 )

models = set() 
# =============================================================================
# use some functions  to parameterize efficiciency
def test_constraint () :

    logger = getLogger ( 'test_constraint' )


    r0 , f0 = model.fitTo ( data , silent = True )
    
    with wait ( 1 ) , use_canvas ( "Unconstrained fit" ) : 
        r0 , f0 = model.fitTo ( data , silent = True , draw = True , nbins = 100 ,
                                minos = ( 'S', 'mean_Gauss' , 'sigma_Gauss' ) )  
        logger.info ( "Uncontrained fit\n%s" % r0.table ( prefix = '# ' ) )


    ## prepare asymmetric constraints 
    a0 = model.soft_constraint2 ( signal.sigma , 1 , -0.7  , 0.05 ) 
    a1 = model.soft_constraint2 ( signal.mean  , 5 , -0.05 , 5    )
    
    with wait ( 1 ) , use_canvas ( "Constrained (asymmetric) fit" ) : 
        r1 , f1 = model.fitTo ( data , silent = True ,
                                constraints = ( a0 , a1 ) )        
        r1 , f1 = model.fitTo ( data , silent = True , draw = True , nbins = 100 ,                                
                                constraints = ( a0 , a1 ) ,        
                                minos = ( 'S', 'mean_Gauss' , 'sigma_Gauss' ) )  
        logger.info ( "Constrained (asymmetric) fit\n%s" % r1.table ( prefix = '# ' ) )


    ## prepare symmetric constrains 
    s0 = model.soft_constraint ( signal.sigma , VE ( 1 , 0.15**2 ) )
    s1 = model.soft_constraint ( signal.mean  , VE ( 5 , 0.10**2 ) )
    
    with wait ( 1 ) , use_canvas ( "Constrained (symmetric) fit" ) : 
        r2 , f2 = model.fitTo ( data , silent = True ,
                                constraints = ( s0 , s1 ) )        
        r2 , f2 = model.fitTo ( data , silent = True , draw = True , nbins = 100 ,                                
                                constraints = ( s0 , s1 ) ,        
                                minos = ( 'S', 'mean_Gauss' , 'sigma_Gauss' ) )  
        logger.info ( "Symmetric constrained fit\n%s" % r2.table ( prefix = '# ' ) )

    models.add ( model ) 

# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        for m in models :
            db['model:' + m.name ] = m
            db['roo_tot:%s' % m.name ] = m.pdf
            for i,s in enumerate ( m.signals ) :
                db['roo_sig%d:%s' % ( i , m.name ) ] = s
            for i, b in enumerate ( m.backgrounds ) : 
                db['roo_bkg%d:%s' % ( i , m.name ) ] = s
            for a in m.alist1 : 
                db['cmp:%s/%s' % ( m.name , a.name ) ] = a
        db['models'   ] = models
        db.ls() 
  

    
# =============================================================================
if '__main__' == __name__ :


    test_constraint ()
        
    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()

# =============================================================================
##                                                                      The END 
# ============================================================================= 
