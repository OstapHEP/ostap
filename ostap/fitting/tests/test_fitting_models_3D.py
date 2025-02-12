#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_models_3D.py
# ============================================================================= 
""" Test module for ostap/fitting/models_3D.py
"""
# ============================================================================= 
from   ostap.core.core          import Ostap, std, VE, dsID, rooSilent 
from   ostap.utils.timing       import timing 
from   builtins                 import range
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait, batch_env 
from   ostap.utils.timing       import timing
import ostap.io.zipshelve       as     DBASE
import ostap.fitting.models     as     Models 
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_models_3D' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
m_x     = ROOT.RooRealVar ( 'mass_x' , 'Some test mass(X)' , 0 , 10 )
m_y     = ROOT.RooRealVar ( 'mass_y' , 'Some test mass(Y)' , 0 , 10 )
m_z     = ROOT.RooRealVar ( 'mass_z' , 'Some test mass(Z)' , 0 , 10 )

models = set()

def test_gauss3D() :

    logger = getLogger("test_gauss3D")

    logger.info ( "Test Gauss3D_pdf" )

    model = Models.Gauss3D_pdf ( 'G3' ,
                                 xvar   = m_x   ,
                                 yvar   = m_y   ,
                                 zvar   = m_z   ,
                                 muX    = ( 5   , 1   , 9  ) ,
                                 muY    = ( 5   , 1   , 9  ) ,
                                 muZ    = ( 5   , 1   , 9  ) ,
                                 sigmaX = ( 0.5 , 0.1 , 5  ) ,
                                 sigmaY = ( 1.5 , 0.1 , 6  ) ,
                                 sigmaZ = ( 2.5 , 0.1 , 6  ) ,
                                 phi    = ( 0.5 , -1  , 10 ) ,
                                 theta  = ( 1.0 , -1  , 10 ) ,
                                 psi    = ( 1.3 , -1  , 10 ) )
    
    dataset = model.generate(1000)
    
    result , _ = model.fitTo ( dataset , silent = True )
    result , _ = model.fitTo ( dataset , silent = True )
    
    logger.info( 'Simple 3D-Gaussian model\n%s' % result.table ( prefix = "# " ) )
    
    with wait ( 2 ), use_canvas ( 'test_gauss3D_x' ) : 
        model.draw1 (  dataset  )
    with wait ( 2 ), use_canvas ( 'test_gauss3D_y' ) : 
        model.draw2 (  dataset  )
    with wait ( 2 ), use_canvas ( 'test_gauss3D_z' ) : 
        model.draw3 (  dataset  )
        
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
        db['m_y'     ] = m_z
        for m in models :
            db['model:' + m.name ] = m
            db['roo_tot:%s' % m.name ] = m.pdf
            for i,s in enumerate ( m.signals ) :
                db['roo_sig%d:%s' % ( i , m.name ) ] = s
            for i, b in enumerate ( m.backgrounds ) : 
                db['roo_bkg%d:%s' % ( i , m.name ) ] = s
            for a in m.alist1 : 
                db['cmp:%s/%s' % ( m.name , a.name ) ] = a
        db['models'  ] = models
        db.ls()
        
   
# =============================================================================
if '__main__' == __name__ :
    
    with timing ('test_gauss3D'   , logger ) :
        test_gauss3D   ()
            
    ## check finally that everything is serializeable:
    with timing ( 'save to DB'     , logger ) : test_db ()          
    
# =============================================================================
##                                                                      The END 
# =============================================================================
