#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_convolution.py
# Test module for ostap/fitting/convolution.py
# ============================================================================= 
""" Test module for ostap/fitting/convolution.py
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.core          import cpp, VE, dsID, rooSilent
from   ostap.utils.timing       import timing, wait 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.root_utils   import batch_env 
import ostap.fitting.models     as     Models 
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_convolution' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch form environment 
batch_env ( logger )
# =============================================================================

## make
x = ROOT.RooRealVar ( 'x',  'test' , 1 , 10 )
models = set()

# =============================================================================
## Asymmetric Laplace 
# =============================================================================
def test_laplace(): 

    logger = getLogger ( 'test_laplace' )
    
    logger.info ('Test Asymmetric Laplace shape' )
    laplace = Models.AsymmetricLaplace_pdf ( name  = 'AL', 
                                             xvar  = x   ,
                                             mean  = 5   , 
                                             slope = 1   )
    
    from ostap.fitting.convolution import  Convolution_pdf

    ## constant resolution  
    laplace_1 = Convolution_pdf ( name = 'L1' , pdf = laplace, resolution = 0.75 )

    ## resolution PDF
    from ostap.fitting.resolution import ResoApo2
    rAp = ResoApo2 ( 'A' , x , 0.75  )
    
    ## resolution as PDF 
    laplace_2 = Convolution_pdf ( name = 'L2' , pdf = laplace, resolution = rAp )


    ## resolution as PDF via operator 
    laplace_3  = laplace % rAp

    with use_canvas ( 'test_laplace' ) : 
        f  = laplace  .draw ( silent = True )        
        f1 = laplace_1.draw ( silent = True )        
        f2 = laplace_2.draw()
        f3 = laplace_3.draw()
        
        with wait ( 2 ) :
            f .draw()
            f1.draw('same')
            f2.draw('same')
            f3.draw('same')
        
    logger.info ( 'Convolution/1:\n%s' % laplace_1.cnv.table ( title = 'Laplace/1' , prefix = '# ' ) ) 
    logger.info ( 'Convolution/2:\n%s' % laplace_2.cnv.table ( title = 'Laplace/2' , prefix = '# ' ) ) 
    logger.info ( 'Convolution/3:\n%s' % laplace_3.cnv.table ( title = 'Laplace/3' , prefix = '# ' ) ) 

    models.add ( laplace   )
    models.add ( laplace_1 )
    models.add ( laplace_2 )
    models.add ( laplace_3 )
    
# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_laplace' )

    logger.info('Saving all objects into DBASE')
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing('Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        for  i , m in enumerate ( models ) :
            db['model/%-2d: %s' % ( i , m.name ) ] = m 
            db['roo/%-2d: %s'   % ( i , m.name ) ] = m.pdf
        db['models'   ] = models
        db.ls() 


# =============================================================================
if '__main__' == __name__ :

    ## Laplace-function + background
    with timing('Laplace'    , logger ) :
        test_laplace () 

    
    ## check finally that everything is serializeable:
    with timing('Save to DB' , logger ) :
        test_db ()          

# =============================================================================
##                                                                      The END 
# ============================================================================= 
