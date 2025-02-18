#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_weights.py
# - It tests various issues for  weighted datasets and minuit 
# ============================================================================= 
""" It tests various issues for  weighted datasets and minuit 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.core       import cpp, VE, dsID, rooSilent 
from   ostap.utils.timing    import timing
from   ostap.plotting.canvas import use_canvas
from   ostap.utils.utils     import wait, batch_env  
import ostap.fitting.models  as     Models 
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_minuit_weighted' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

# =============================================================================
def test_minuit_weighted () : 

    logger = getLogger ( "test_minuit_weighted" )
    
    ## make simple test variable
    xvar    = ROOT.RooRealVar ( 'test_x' , 'Some test x' , 0 , 10 )
    yvar    = ROOT.RooRealVar ( 'test_y' , 'Some test y' , 0 , 1  )
    ## book very simple data set
    varset  = ROOT.RooArgSet  ( xvar , yvar )
    dataset = ROOT.RooDataSet ( dsID() , 'Test Data set-0' , varset )  
    
    NS = 10000
    NB = 10000
    
    xmin , xmax = xvar.minmax()
    ymin , ymax = yvar.minmax()
    
    ## fill it 
    x = VE(5.0,1**2)
    for i in range(0,NS) :
        xvar.value = x.gauss ()
        y = -1
        while not ymin <= y <= ymax : y = random.expovariate ( 1/0.1 ) 
        yvar.value = y 
        dataset.add ( varset )
    
    for i in range(0,NB) :
        xvar.value = random.uniform ( xmin , xmax ) 
        y = -1
        while not ymin <= y <= ymax : y = random.expovariate ( 1/0.9 ) 
        yvar.value = y 
        dataset.add ( varset   )

    logger.info ('DATASET\n%s' % dataset.table( prefix = '# ' ) ) 

    signal  = Models.Gauss_pdf ('G' , xvar = xvar , mean = ( 4, 6 ) , sigma =   ( 0.5 , 1.5 ) )
    model   = Models.Fit1D ( signal = signal , background = 0 )
    model.S = NS 
    model.B = NB 

    with use_canvas ( "test_minuit_weighted" ) , wait ( 2 ) : 
        r0 , _  = model.fitTo ( dataset , silent = True , draw = False )
        r0 , f0 = model.fitTo ( dataset , silent = True , draw = True  , nbins = 100 )
    
    logger.info ('Fit result\n%s' % r0.table ()  )
    
    ## make splot
    model.sPlot ( dataset )
    dsw =  dataset.makeWeighted ( "S_sw" )
    
    logger.info ('Weighted DATASET\n%s' % dsw.table ( prefix = '# ' ) ) 
    
    ## fitting function for the  exponential stuff
    expf = Models.Bkg_pdf ( 'E' , xvar = yvar , power = 0 )

    rw , fw = expf.fitTo ( dsw , silent = True , draw = False )
    rw , fw = expf.fitTo ( dsw , silent = True , draw = True  , nbins = 50 )

    logger.always ( str ( type ( rw ) ) )
    
    logger.info ('Weighted fit result (incorrect)\n%s' % rw.table ( title = 'incorrect')  )
    
    r2 , f2 = expf.fitTo ( dsw , silent = True , sumw2 = True , draw = False )
    r2 , f2 = expf.fitTo ( dsw , silent = True , sumw2 = True , draw = True  , nbins = 50 )
    
    ## r2.Print('vvv')
    
    logger.info ('Weighted fit result (SumW2=True)\n%s' % r2.table ( 'SumW2=True')  )
    
    with use_canvas ( "test_minuit_weighted" ) , wait ( 2 ) : 
        ra , fa = expf.fitTo ( dsw , asymptotic = True , draw = False , silent = True )
        ra , fa = expf.fitTo ( dsw , asymptotic = True , draw = True  , silent = True , nbins = 50 )
        logger.info ('Weighted fit result (asymptotic=True)\n%s' % ra.table ( 'Asymptotic=True')  )
        
    m = expf.minuit ( dataset = dsw , scale = False , silent = True )
    m.migrad ( 5 )
    rm = m.save() 
    
    logger.info ('MiNUIT fit result (incorrect) \n%s' % rm.table ( title = 'incorrect')  )
    
    ms = expf.minuit ( dataset = dsw , scale = True , silent = True )
    ms.migrad ( 5 )
    rs = ms.save() 
    
    logger.info ('MiNUIT (scale=True) fit result\n%s' % rs.table ( title = 'scaled')  )
    
    
# =============================================================================
if '__main__' == __name__ :

    with timing ("minuit weighted" , logger ) :
        test_minuit_weighted ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
