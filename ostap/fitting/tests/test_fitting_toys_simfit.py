#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_toys_simfit.py
# Test toys+simfit
# ============================================================================= 
""" Test toys + simfit
"""
# ============================================================================= 
from builtins    import range 
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import dsID, rooSilent
import ostap.fitting.toys       as     Toys
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_simfit1' )
else : 
    logger = getLogger ( __name__ )

        
# =============================================================================
def test_toys_simfit_1 () :

    logger = getLogger ( 'test_toys_simfit_1' )
    
    ## make simple test mass 
    mass     = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 0 , 5 )
    
    ## book very simple data set:
    varset1  = ROOT.RooArgSet  ( mass )
    dataset1 = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset1 )  
    
    ## book very simple data set:
    varset2  = ROOT.RooArgSet  ( mass )
    dataset2 = ROOT.RooDataSet ( dsID() , 'Test Data set-2' , varset2 )
    
    ## high statistic, low-background "control channel"
    mean1  = 2.0
    sigma1 = 0.50
    NS1    =   1000
    NB1    =    250
    
    for i in range ( NS1 )  :
        v1 = random.gauss( mean1 , sigma1 )
        if v1 in mass :
            mass.setVal ( v1 )
            dataset1.add ( varset1 )
            
    for i in range ( NB1 ) :
        v1 = random.uniform ( 0 , 5  )
        if v1 in mass :
            mass.setVal ( v1 )
            dataset1.add ( varset1 )
            
    ## low statistic, high-background "control channel"
    NS2     =    250
    NB2     =   1000     
    mean2   = mean1  + 1.0
    sigma2  = sigma1 * 0.5
    
    for i in range(NS2)  :
        v2 = random.gauss( mean2 , sigma2 )
        if v2 in mass :
            mass.setVal ( v2 )
            dataset2.add ( varset2 )
    for i in range ( NB2 ) :
        v2 = random.uniform ( 0 , 5  )
        if v2 in mass :
            mass.setVal ( v2 )
            dataset2.add ( varset2 )

            
    signal1  = Models.Gauss_pdf ( 'G1'                 ,
                                  xvar  = mass         ,
                                  mean  = (0.5 , 2.5 ) ,
                                  sigma = (0.1 , 1.0 ) )
    
    model1   = Models.Fit1D ( suffix = 'M1' , signal = signal1 ,  background = -1 )
    model1.S = NS1
    model1.B = NB1 
    
    mean2    = signal1.vars_add      ( signal1.mean  , 1.0  )
    sigma2   = signal1.vars_multiply ( signal1.sigma , 0.5  )
    
    signal2  = Models.Gauss_pdf ( 'G2'            ,
                                  xvar  = mass    ,
                                  mean  = mean2   ,
                                  sigma = sigma2  )
    
    model2  = Models.Fit1D ( suffix = 'M2' , signal = signal2 ,  background = model1.background )
    model2.S = NS2
    model2.B = NB2 

    with use_canvas ( 'test_toys_simfit_1' ) : 
        # =========================================================================
        ## fit 1 
        with wait ( 1 ) :
            r1 , f1 = model1.fitTo ( dataset1 , draw = True , nbins = 50 , silent = True )
        
        ## fit 2
        with wait ( 1 ) :
            r2 , f2 = model2.fitTo ( dataset2 , draw = True , nbins = 50 , silent = True )
        # =========================================================================
    
    ## combine data    
    sample  = ROOT.RooCategory ('sample','sample'  , 'A' , 'B' )
    
    ## combine datasets
    from ostap.fitting.simfit import combined_data 
    vars    = ROOT.RooArgSet ( mass )
    dataset = combined_data  ( sample , vars , { 'A' : dataset1 , 'B' : dataset2 } )
    
    ## combine PDFs
    model_sim  = Models.SimFit (
        sample , { 'A' : model1  , 'B' : model2 } , name = 'X'
        )
    
    # =========================================================================
    r , f = model_sim.fitTo ( dataset , silent = True )
    r , f = model_sim.fitTo ( dataset , silent = True )

    
    with use_canvas ( 'test_toys_simfit_1' ) : 
        with wait ( 1 ) : fA = model_sim.draw ( 'A' , dataset , nbins = 50 )
        with wait ( 1 ) : fB = model_sim.draw ( 'B' , dataset , nbins = 50 )
        
    logger.info ( 'Fit  results are: %s ' % r.table ( prefix = "# " ) ) 

    ## Make toys
    results , stats = Toys.make_toys  (
        pdf         = model_sim ,
        nToys       = 100       ,
        data        = [ mass ]  , 
        gen_config  = { 'nEvents' : ( NS1 + NB1 , NS2 + NB2 ) , 'sample' : True } ,
        fit_config  = { 'silent'  : True } ,
        init_pars   = { 'mean_G1'  : mean1   , 'BM1' : NB1 , 'BM2' : NB2 , 
                        'sigma_G1' : sigma1  , 'SM1' : NS1 , 'SM2' : NS2 ,
                        'phi0_Bkg_FitG1_M1'  :  0
                        } ,
        silent      = True , 
        progress    = True )

# =============================================================================
if '__main__' == __name__ :

    with timing ( "simfit-1" , logger ) :
        test_toys_simfit_1 () 

# =============================================================================
##                                                                      The END 
# =============================================================================
