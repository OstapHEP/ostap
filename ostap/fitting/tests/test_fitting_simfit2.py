#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_simfit2.py
# Test module for ostap/fitting/simfit.py
# - It tests the most simple "Simultaneous fit"
#
# Simultaneous fit for two different 1D-distributions/ranges/observables 
# ============================================================================= 
""" Test module for ostap/fitting/simfit.py
- It tests the most simple ``Simultaneous fit''

Simultaneous fit for two different 1D-distributions/ranges/observables 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   builtins    import range 
from   ostap.core.core          import dsID, rooSilent 
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_simfit2' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make simple test mass 
mass1    = ROOT.RooRealVar ( 'test_mass1' , 'Some test mass' ,  0 , 10 )
mass2    = ROOT.RooRealVar ( 'test_mass2' , 'Some test mass' , 10 , 20 )

## book very simple data set:
varset1  = ROOT.RooArgSet  ( mass1 )
dataset1 = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset1 )  

## book very simple data set:
varset2  = ROOT.RooArgSet  ( mass2 )
dataset2 = ROOT.RooDataSet ( dsID() , 'Test Data set-2' , varset2 )  

## high statistic, low-background "control channel"
mean1  = 5.0
sigma1 = 1.0
NS1    =  10000
NB1    =   1000

for i in range ( NS1 )  :
    v1 = random.gauss( mean1 , sigma1 )
    if v1 in mass1 :
        mass1.setVal ( v1 )
        dataset1.add ( varset1 )

for i in range ( NB1 ) :
    v1 = random.uniform ( 0 , 10  )
    if v1 in mass1 :
        mass1.setVal  ( v1 )
        dataset1.add ( varset1 )
        
## low statistic, high-background "control channel"
NS2     =    500
NB2     =  10000     
mean2   = mean1  + 10.0
sigma2  = sigma1 *  0.5

for i in range(NS2)  :
    v2 = random.gauss( mean2 , sigma2 )
    if v2 in mass2 :
        mass2.setVal ( v2 )
        dataset2.add ( varset2 )

for i in range (NB2 ) :
    v2 = random.uniform ( 10 , 20  )
    if v2 in mass2 :
        mass2.setVal ( v2 )
        dataset2.add ( varset2 )

# =============================================================================
def test_simfit2 ( ) :

    logger = getLogger ( 'test_simfit2' ) 
    # =========================================================================    
    signal1  = Models.Gauss_pdf ( 'G1'                 ,
                                  xvar  = mass1        ,
                                  mean  = (1.5 , 6.5 ) ,
                                  sigma = (0.1 , 2.5 ) )
    
    model1   = Models.Fit1D ( suffix = 'M1' , signal = signal1 ,  background = -1 )
    model1.S = NS1
    model1.B = NB1 
    
    
    mean2    = signal1.vars_add      ( signal1.mean  , 10.0 )
    sigma2   = signal1.vars_multiply ( signal1.sigma , 0.5  )
    
    signal2  = Models.Gauss_pdf ( 'G2'            ,
                                  xvar  = mass2   ,
                                  mean  = mean2   ,
                                  sigma = sigma2  )
    
    model2  = Models.Fit1D ( suffix = 'M2' , signal = signal2 ,  background = -1  )
    model2.S = NS2
    model2.B = NB2 

    with use_canvas ( 'test_simfit2' ) : 
        # =========================================================================
        ## fit 1
        with wait ( 1 ) : 
            r1 , f1 = model1.fitTo ( dataset1 , draw = True , nbins = 50 , silent = True )        
            title = 'Results of fit to dataset1'
            logger.info ( '%s\n%s' % ( title , r1.table ( title = title , prefix = '# ' ) ) )
        ## fit 2
        with wait ( 1 ) : 
            r2 , f2 = model2.fitTo ( dataset2 , draw = True , nbins = 50 , silent = True )
            title = 'Results of fit to dataset2'
            logger.info ( '%s\n%s' % ( title , r2.table ( title = title , prefix = '# ' ) ) )
        # =========================================================================
    
    ## combine data
    sample  = ROOT.RooCategory ('sample','sample'  , 'A' , 'B' )
    
    ## combine datasets
    from ostap.fitting.simfit import combined_data 
    vars    = ROOT.RooArgSet ( mass1 , mass2 )
    dataset = combined_data  ( sample , vars , { 'A' : dataset1 , 'B' : dataset2 } )
    
    
    ## combine PDFs
    model_sim  = Models.SimFit (
        sample , { 'A' : model1  , 'B' : model2 } , name = 'X'
        )
    
    # =========================================================================
    r , f = model_sim.fitTo ( dataset , silent = True )
    r , f = model_sim.fitTo ( dataset , silent = True )
    
    title = 'Results of simultaneous fit'
    logger.info ( '%s\n%s' % ( title , r.table ( title = title , prefix = '# ' ) ) )

    with use_canvas ( 'test_simfit2' ) : 
        with wait ( 1 ) : fA    = model_sim.draw ( 'A' , dataset , nbins = 50 )
        with wait ( 1 ) : fB    = model_sim.draw ( 'B' , dataset , nbins = 50 )
        
    ## fNLL  = model_sim.draw_nll ( 'SM2' , dataset , range =   (0,1000)  )
    ## significance 
    ## wilks = model_sim.wilks    ( 'SM2' , dataset  )
    


# =============================================================================
if '__main__' == __name__ :

    with timing ("simfit-2", logger) :  
        test_simfit2 () 

# =============================================================================
##                                                                      The END 
# =============================================================================
