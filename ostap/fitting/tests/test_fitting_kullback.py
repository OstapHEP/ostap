#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_kullback.py
# simpel tetx module 
# ============================================================================= 
""" Simple test module 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   builtins                 import range
from   ostap.core.core          import SE, hID  
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models
from   ostap.fitting.fithelpers import SETPARS
from   ostap.utils.progress_bar import progress_bar 
from   ostap.utils.timing       import timing
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait
import ostap.histos.histos  
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_kullback' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

x        = ROOT.RooRealVar ( 'x' , 'x-observable' , -10 , 10 )

signal  = Models.Gauss_pdf ( 'G0' , xvar = x ,
                             mean  = ( 0 ,-0.5 , 0.5) ,
                             sigma = ( 1 , 0.1 , 2.5) )
model   = Models.Fit1D ( signal = signal , background = -1 )

NS = 700
NB = 300

model.S = NS 
model.B = NB 

with SETPARS ( model ) : 
    dataset0 = model.generate ( NS + NB )
    r0 , f0  = model.fitTo ( dataset0 , silent = True , draw = True )   
    logger.info ( 'Fit result\n%s' % r0 )
    
def test_kullback () :
    
    with SETPARS ( model ) : 
        r , f = model.fitTo ( dataset0 , silent = True ) 

    cnt = SE()
    hh  = ROOT.TH1F ( hID() , 'Kullback Leibler divergency' , 10 , 0 , 100 )
    
    for i in progress_bar ( range ( 100 ) ) :
        
        ds = model.generate ( NS + NB )
        
        with SETPARS ( model ) :
            
            r , f = model.fitTo ( ds , silent = True ) 
            kl = r.kullback_leibler ( r0 )
            cnt += kl
            hh.Fill ( kl ) 
                
        ds.clear()
        del ds
            
    with wait ( 2 ) , use_canvas ( 'test_kullback' ) :

        hh.draw()
        logger.info ( 'Counter %s' % cnt )
        x
# =============================================================================
if '__main__' == __name__ :


    test_kullback()
    
# =============================================================================
##                                                                      The END
# =============================================================================
