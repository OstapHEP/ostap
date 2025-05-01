#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_transform.py
# Test module 
# - It tests transform-PDF  
# ============================================================================= 
""" Test module 
- It tests transform-PDF 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.meta_info        import root_info 
from   ostap.core.core             import cpp, VE, dsID
from   ostap.fitting.transform_pdf import TrPDF
from   ostap.fitting.funbasic      import Fun1D
from   ostap.utils.timing          import timing 
from   ostap.plotting.canvas       import use_canvas
from   ostap.utils.root_utils      import batch_env 
import ostap.fitting.models        as     Models
import ostap.fitting.roofit
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_transform' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

# =============================================================================
def  test_transform () :

    logger = getLogger ( 'test_transform' )
    
    if not (6,20) <= root_info :
        logger.info ("Not for this version of ROOT!")
        return 

    lmin = 0
    lmax = 5
    
    x  = ROOT.RooRealVar('x','', 10.**lmin , 10.**lmax )
    mean  =  1000
    sigma =  1000 

    ## book very simple data set
    varset  = ROOT.RooArgSet  ( x  )
    dataset = ROOT.RooDataSet ( dsID() , 'Test Data set' , varset )  
    while len ( dataset ) < 5000 : 
        v = random.gauss ( mean , sigma )
        if v in x :
            x.value = v 
            dataset.add ( varset  )
            
    dataset.add_var ( 'lx' , 'log10(x)' ) 
    dataset.lx.setMin ( lmin ) 
    dataset.lx.setMax ( lmax )
    
    ## original PDF 
    gauss = Models.Gauss_pdf ( 'G' , xvar = x ,
                           mean  = ( mean  , mean  / 2 , mean  * 2 ) , 
                           sigma = ( sigma , sigma / 2 , sigma * 2 ) )
    
    with use_canvas ( 'test_transform/x' , wait = 1 ) : 
        r1 , f1  = gauss.fitTo  ( dataset , draw = True  , silent = True )
        logger.info ( 'Fit x:\n%s' % r1.table() ) 
    
    lx = dataset.lx
    LX = Fun1D (  lx , lx )
    NX = 10 ** LX
    
    ## transformed PDF 
    tgauss  = TrPDF ( pdf = gauss , new_var = NX )
    
    with use_canvas ( 'test_transform/log10' , wait = 1 ): 
        r2 , f2 = tgauss.fitTo  ( dataset , draw = True  , silent = True )
        logger.info ( 'Fit log10(x):\n%s' % r2.table() ) 

    
# =============================================================================
if '__main__' == __name__ :

    with timing ("transform" , logger ) :  
        test_transform    () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
