#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_models.py
# Test module for ostap/fitting/models.py
# - It tests various ``signal-like''/``peak-like'' shapes
# ============================================================================= 
""" Test module for ostap/fitting/models.py
- It tests various ``signal-like''/``peak-like'' shapes 
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
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_models' )
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

    laplace.draw()
    laplace_1.draw()
    laplace_2.draw()

    models.add ( laplace  )
    models.add ( laplace_1 )
    models.add ( laplace_2 )
    
# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :
    logger.info('Saving all objects into DBASE')
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( name = 'Save everything to DBASE'), DBASE.tmpdb() as db : 
        db['models'   ] = models
        db.ls() 
        
# =============================================================================
if '__main__' == __name__ :

    test_laplace        () ## Laplace-function                            + background 
    
    ## check finally that everything is serializeable:
    test_db ()          

# =============================================================================
# The END 
# ============================================================================= 
