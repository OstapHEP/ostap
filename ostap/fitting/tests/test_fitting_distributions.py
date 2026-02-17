#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_distributions.py
# Test module for ostap/fitting/distributions.py
# ============================================================================= 
""" Test module for ostap/fitting/distributions.py
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.utils.timing          import timing 
from   ostap.plotting.canvas       import use_canvas
from   ostap.utils.root_utils      import batch_env 

from   ostap.utils.basic           import typename 
import ostap.fitting.models        as     Models 
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_distributions' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================


x = ROOT.RooRealVar ( 'x' , 'x-variable' ,  0  , 10 ) 
y = ROOT.RooRealVar ( 'y' , 'y-variable' , -5 ,   5 ) 


the_plots  = set ()
the_models = set ()

# =============================================================================
def _check_the_model_ ( model_type ) :


    model  = model_type ( xvar = x )
    logger = getLogger("test_%s" % model.name )
    
    with use_canvas ( 'test %s' % model.name , wait = 3 ) :
        logger.info ( 'Test the model %s:\n%s' % ( typename ( model ) , model ) )
        plot = model.draw()
    
    the_models.add ( model )
    the_plots .add ( plot  ) 

# =============================================================================
def test_models () :

    from   ostap.fitting.distributions import models , spectra
    models = list ( models )
    logger.info ( 'Test distributions' )
    for model in models :
        if not model in spectra : _check_the_model_ ( model ) 
        
# =============================================================================
def test_Tsallis () :

    logger = getLogger("test_Tsallis")
    
    model  = Models.Tsallis_pdf ( name = 'Tsallis' , 
                                  xvar = x         ,
                                  m0 = ROOT.RooFit.RooConst ( 0.135 ) ,
                                  n  = 10   ,
                                  T  = 0.2  )
    
    with use_canvas ( 'Tsallis_pdf' ) :
        plot = model.draw()
        
    the_models.add ( model )
    the_plots .add ( plot  ) 

# =============================================================================
def test_QGSM () :

    logger = getLogger("test_QGSM")
    
    model  = Models.QGSM_pdf ( name = 'QGSM'   , 
                               xvar = x  ,
                               m0 = ROOT.RooFit.RooConst ( 0.135 ) ,
                               b  = 10  )
    
    with use_canvas ( 'QGSM_pdf' ) :
        plot = model.draw()
        
    the_models.add ( model )
    the_plots .add ( plot  ) 

# =============================================================================
def test_Hagedorn () :

    logger = getLogger("test_Hagedorn")
    
    model  = Models.Hagedorn_pdf ( name = 'Hagedorn' , 
                                   xvar = x ,
                                   m0   = 1   ,
                                   beta = 10  )
    
    with use_canvas ( 'QGSM_pdf' ) :
        plot = model.draw()
        
    the_models.add ( model )
    the_plots .add ( plot  ) 





# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        for m in the_models :
            db['model:' + m.name ] = m
            db['roo:%s' % m.name ] = m.pdf
        db['models'   ] = the_models
        db['plots'    ] = the_plots 
        db.ls() 

# =============================================================================
if '__main__' == __name__ :
    
    test_models       ()
    
    test_Tsallis      ()
    test_QGSM         ()
    test_Hagedorn     ()

    ## check that everything is serializeable 
    test_db ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
