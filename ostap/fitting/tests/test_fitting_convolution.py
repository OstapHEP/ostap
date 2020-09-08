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
import ROOT, random
import ostap.fitting.roofit 
import ostap.fitting.models as     Models 
from   ostap.core.core      import cpp, VE, dsID
from   ostap.logger.utils   import rooSilent
from   ostap.utils.timing   import timing
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_convolution' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make
x = ROOT.RooRealVar ( 'x',  'test' , 1 , 10 )
models = set()

# =============================================================================
## Asymmetric Laplace 
# =============================================================================
def test_laplace(): 

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

    laplace.draw( silent = True )

    laplace_1.draw( silent = True )

    laplace_2.draw()

    models.add ( laplace  )
    models.add ( laplace_1 )
    models.add ( laplace_2 )
    
# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    from ostap.core.meta_info import root_version_int 
    if root_version_int >= 62200 :
        logger.warning("test_db: test is disabled for ROOT version %s" % root_version_int )
        return 

    
    logger.info('Saving all objects into DBASE')
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing('Save everything to DBASE', logger ), DBASE.tmpdb() as db : 
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
