#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_models7_2D.py
# ============================================================================= 
""" Test module for ostap/fitting/models7_2d.py
"""
# ============================================================================= 
from   __future__               import print_function
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import Ostap, std, VE, dsID
from   ostap.logger.utils       import rooSilent 
import ostap.io.zipshelve       as     DBASE
from   ostap.utils.timing       import timing 
from   builtins                 import range
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait
from   ostap.utils.timing import timing
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_models7_2D' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
m_x     = ROOT.RooRealVar ( 'mass_x' , 'Some test mass(X)' , 0 , 10 )
m_y     = ROOT.RooRealVar ( 'mass_y' , 'Some test mass(Y)' , 0 , 10 )

models = set()

def test_gauss2D() :

    logger = getLogger("test_gauss2D")

    logger.info ( "Test Gauss2D_pdf" )

    model = Models.Gauss2D_pdf ( 'G2' ,
                                 xvar   = m_x   ,
                                 yvar   = m_y   ,
                                 muX    = ( 5   , 1   , 9  ) ,
                                 muY    = ( 5   , 1   , 9  ) ,
                                 sigmaX = ( 0.5 , 0.1 , 5  ) ,
                                 sigmaY = ( 2.5 , 0.1 , 6  ) ,
                                 theta  = ( 1.5 , -1  , 10 ) )
    dataset = model.generate(5000)
    
    result , f = model.fitTo ( dataset , silent = True )
    result , f = model.fitTo ( dataset , silent = True )
    
    logger.info( 'Simple 2D-Gaussian model\n%s' % result.table ( prefix = "# " ) )
    
    with wait ( 2 ), use_canvas ( 'test_gauss2D_x' ) : 
        model.draw1 (  dataset  )
    with wait ( 2 ), use_canvas ( 'test_gauss2D_y' ) : 
        model.draw2 (  dataset  )
        
    models.add ( model )
        
    
 # =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info('Saving all objects into DBASE')
    with timing('Save everything to DBASE', logger ), DBASE.tmpdb() as db : 
        db['m_x'     ] = m_x
        db['m_y'     ] = m_y
        for m in models : db['model:' + m.name ] = m
        db['models'  ] = models
        db.ls()
        
   
# =============================================================================
if '__main__' == __name__ :
    
    with timing ('test_gauss2D'   , logger ) :
        test_gauss2D   ()
            
    ## check finally that everything is serializeable:
    with timing ( 'save to DB'     , logger ) : test_db ()          
    
# =============================================================================
##                                                                      The END 
# =============================================================================
