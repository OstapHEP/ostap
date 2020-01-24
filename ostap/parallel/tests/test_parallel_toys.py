#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_toys.py
# Test module for ostap/fitting/toys.py
# - make some fitting toys 
# ============================================================================= 
""" Test module for ostap/fitting/toys.py
- make some fitting toys 
"""
# ============================================================================= 
from   __future__        import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ROOT, random
import ostap.fitting.roofit 
import ostap.fitting.models as     Models 
from   ostap.core.core      import cpp, VE, dsID
from   ostap.logger.utils   import rooSilent
from   builtins             import range
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_toys' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
import ROOT, time 
from   ostap.core.pyrouts           import hID 
import ostap.fitting.models         as     Models
import ostap.parallel.parallel_toys as     Toys

mass      = ROOT.RooRealVar ( 'mass' , '', 0 , 1 )  
gen_gauss = Models.Gauss_pdf ( 'GG' , xvar = mass )
fit_gauss = Models.Gauss_pdf ( 'FG' , xvar = mass )
gen_gauss.mean  = 0.4
gen_gauss.sigma = 0.1

# ==============================================================================
## Perform toy-study for possible fit bias and correct uncertainty evaluation
#  - generate <code>nToys</code> pseudoexperiments with some PDF <code>pdf</code>
#  - fit teach experiment with the same PDF
#  - store  fit results
#  - calculate staistics of pulls
#  - fill distributions for fit results
#  - fill distribution of pulls 
def test_parallel_toys ( ) :
    """Perform toys-study for possible fit bias and correct uncertainty evaluation
    - generate `nToys` pseudoexperiments with some PDF `pdf`
    - fit teach experiment with the same PDF
    - store  fit results
    - calculate statistics of pulls
    - fill distributions of fit results
    - fill distributions of pulls 
    """
    
    results , stats = Toys.parallel_toys  (
        pdf         = gen_gauss ,
        nToys       = 1000      ,
        nSplit      = 20        ,
        data        = [ mass ]  , 
        gen_config  = { 'nEvents' : 200  } ,
        fit_config  = { 'silent'  : True } ,
        init_pars   = { 'mean_GG' : 0.4 , 'sigma_GG' : 0.1 } ,
        silent      = True  , 
        progress    = False )

    for p in stats :
        logger.info (  "Toys: %-20s : %s" % (  p, stats [ p ] ) )

    ## make histos:
    
    h_mean       = ROOT.TH1F ( hID() , 'mean of Gauss ' , 100 ,  0    ,  0.80 )
    h_sigma      = ROOT.TH1F ( hID() , 'sigma of Gauss' , 100 ,  0.05 ,  0.15 )
    
    for r in results [ 'mean_GG'  ] : h_mean  .Fill ( r ) 
    for r in results [ 'sigma_GG' ] : h_sigma .Fill ( r )

    for h in ( h_mean , h_sigma ) :
        
        h.draw()
        logger.info ( "%s  :\n%s"  % ( h.GetTitle() , h.dump ( 30 , 10 ) ) )
        time.sleep ( 1 )

# =============================================================================
## Perform toy-study for possible fit bias and correct uncertainty evaluation
#  - generate <code>nToys</code> pseudoexperiments with some PDF <code>gen_pdf</code>
#  - fit teach experiment with the PDF <code>fit_pdf</code>
#  - store  fit results
#  - fill distributions for fit results
def test_parallel_toys2 ( ) :    
    """Perform toys-study for possible fit bias and correct uncertainty evaluation
    - generate `nToys` pseudoexperiments with some PDF `gen_pdf`
    - fit teach experiment with the PDF `fit_pdf`
    - store  fit results
    - fill distributions of fit results
    """    
            
    results , stats = Toys.parallel_toys2 (
        gen_pdf     = gen_gauss ,
        fit_pdf     = fit_gauss ,
        nToys       = 1000      ,
        nSplit      = 20        ,
        data        = [ mass ]  , 
        gen_config  = { 'nEvents' : 200  } ,
        fit_config  = { 'silent'  : True } ,
        gen_pars    = { 'mean_GG' : 0.4 , 'sigma_GG' : 0.1 } ,
        fit_pars    = { 'mean_GF' : 0.4 , 'sigma_GF' : 0.1 } ,
        silent      = True  , 
        progress    = False )
    
    for p in stats :
        logger.info (  "Toys: %-20s : %s" % (  p, stats [ p ] ) )

    ## make histos
        
    h_mean       = ROOT.TH1F ( hID() , 'mean of Gauss ' , 50 ,  0    ,  0.80 )
    h_sigma      = ROOT.TH1F ( hID() , 'sigma of Gauss' , 50 ,  0.05 ,  0.15 )

    for r in results ['mean_FG'  ] : h_mean .Fill ( r ) 
    for r in results ['sigma_FG' ] : h_sigma.Fill ( r )

    for h in ( h_mean , h_sigma ) :
        
        h.draw()
        logger.info ( "%s  :\n%s"  % ( h.GetTitle() , h.dump ( 30 , 10 ) ) )
        time.sleep ( 1 )

# =============================================================================
## Perform toy-study for significance of the signal 
#  - generate <code>nToys</code> pseudoexperiments using background-only hypothesis 
#  - fit teach experiment with "signal+background" hypothesis
#  - store  fit results
#  - fill distributions for fit results
def test_parallel_significance_toys ( ) :
    """Perform toy-study for significance of the signal 
    - generate `nToys` pseudoexperiments using background-only hypothesis 
    - fit each experiment with signal+background hypothesis
    - store  fit results
    - fill distributions for fit results
    """
    
    ## only background hypothesis
    bkg_only = Models.Bkg_pdf    ( "BKG" , xvar =  mass , power = 0 , tau = 0      )

    
    signal   = Models.Gauss_pdf  ( 'S'   , xvar = mass , mean = 0.5 , sigma = 0.1 )

    signal.mean .fix ( 0.4 ) 
    signal.sigma.fix ( 0.1 ) 


    ## signal + background hypothesis
    model    = Models.Fit1D      ( signal = signal , background = 1 )
    model.background.tau.fix ( 0 )
    
    results , stats = Toys.parallel_toys2 (
        gen_pdf     = bkg_only  ,
        fit_pdf     = model     ,
        nToys       = 1000      ,
        nSplit      = 20        , 
        data        = [ mass ]  , 
        gen_config  = { 'nEvents'  : 100 , 'sample' : True } ,
        fit_config  = { 'silent'   : True } ,
        gen_pars    = { 'tau_BKG'  : 0.   } , ## initial values for generation 
        fit_pars    = { 'B' : 100         ,
                        'S' : 10          ,
                        'phi0_Bkg_S': 0.0 } , ## initial fit values for parameters 
        silent      = True , 
        progress    = True )

    for p in stats :
        logger.info (  "Toys: %-20s : %s" % (  p, stats [ p ] ) )

    h_S      = ROOT.TH1F ( hID() , '#S' , 60 ,  0 , 60 )
    
    for r in results ['S'  ] : h_S .Fill ( r )
    
    for h in ( h_S ,  ) :
        
        h.draw()
        logger.info ( "%s  :\n%s"  % ( h.GetTitle() , h.dump ( 30 , 10 ) ) )
        time.sleep  ( 1 )


# =============================================================================
if '__main__' == __name__ :

    test_parallel_toys  () 
    test_parallel_toys2 () 
    test_parallel_significance_toys ( ) 


# =============================================================================
##                                                                      The END 
# =============================================================================
