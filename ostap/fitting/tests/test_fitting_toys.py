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
## from   __future__        import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.pyrouts       import hID
import ostap.fitting.models     as     Models
import ostap.fitting.toys       as     Toys
import ostap.histos.histos
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
from   ostap.fitting.toys       import pull_var
import ROOT, time, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_toys' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
nominal_mean  = 0.4
nominal_sigma = 0.1

mass      = ROOT.RooRealVar ( 'mass' , '', 0 , 1 )  
gen_gauss = Models.Gauss_pdf ( 'GG' , xvar = mass )
fit_gauss = Models.Gauss_pdf ( 'FG' , xvar = mass )
gen_gauss.mean  = nominal_mean
gen_gauss.sigma = nominal_sigma 
model     = Models.Fit1D ( signal     = gen_gauss , background = 'flat'    )


# ==============================================================================
## Perform toy-study for possible fit bias and correct uncertainty evaluation
#  - generate <code>nToys</code> pseudoexperiments with some PDF <code>pdf</code>
#  - fit each experiment with the same PDF
#  - store  fit results
#  - calculate staistics of pulls
#  - fill distributions for fit results
#  - fill distribution of pulls 
def test_toys ( ) :
    """Perform toys-study for possible fit bias and correct uncertainty evaluation
    - generate `nToys` pseudoexperiments with some PDF `pdf`
    - fit teach experiment with the same PDF
    - store  fit results
    - calculate statistics of pulls
    - fill distributions of fit results
    - fill distributions of pulls 
    """

    logger = getLogger ( 'test_toys' )
    
    more_vars   = { 'pull:mean_GG'  : pull_var ( 'mean_GG'  , nominal_mean  ) ,  
                    'pull:sigma_GG' : pull_var ( 'sigma_GG' , nominal_sigma ) }
    
    with timing ( 'Toys      analysis' , logger = logger )  :        
        results , stats = Toys.make_toys  (
            pdf         = gen_gauss ,
            nToys       = 1000      ,
            data        = [ mass ]  ,
            gen_config  = { 'nEvents'  : 200  , 'sample'  : True } ,
            fit_config  = { 'silent'   : True , 'refit'   : 5   } ,
            init_pars   = { 'mean_GG'  : nominal_mean   ,
                            'sigma_GG' : nominal_sigma  } ,
            more_vars   = more_vars , 
            silent      = True , 
            progress    = True )

    ## make histos:
    
    h_mean       = ROOT.TH1F ( 'h1' , 'mean  of Gauss ' , 100 ,  0                   , nominal_mean  * 2 )   
    h_sigma      = ROOT.TH1F ( 'h2' , 'sigma of Gauss'  , 100 ,  0.5 * nominal_sigma , nominal_sigma * 2 )
    
    for r in results [ 'mean_GG'  ] : h_mean  .Fill ( r ) 
    for r in results [ 'sigma_GG' ] : h_sigma .Fill ( r )
    
    for h in ( h_mean , h_sigma ) :
        with use_canvas ( 'test_toys  %s' % h.title  , wait = 1 ) : h.draw()

# =============================================================================
## Perform toy-study for possible fit bias and correct uncertainty evaluation
#  - generate <code>nToys</code> pseudoexperiments with some PDF <code>gen_pdf</code>
#  - fit teach experiment with the PDF <code>fit_pdf</code>
#  - store  fit results
#  - fill distributions for fit results
def test_toys2 ( ) :
    """Perform toys-study for possible fit bias and correct uncertainty evaluation
    - generate `nToys` pseudoexperiments with some PDF `gen_pdf`
    - fit teach experiment with the PDF `fit_pdf`
    - store  fit results
    - fill distributions of fit results
    """    

    logger = getLogger ( 'test_toys2' )
    
    more_vars   = { 'pull:mean_FG'  : pull_var ( 'mean_FG' , nominal_mean  ) ,  
                    'pull:sigma_FG' : pull_var ( 'sigma_FG', nominal_sigma ) }
    
    with timing ( 'Toys2     analysis' , logger = logger )  :        
        results , stats = Toys.make_toys2 (
            gen_pdf     = gen_gauss ,
            fit_pdf     = fit_gauss ,
            nToys       = 1000      ,
            data        = [ mass ]  , 
            gen_config  = { 'nEvents' : 200  , 'sample' : True } ,
            fit_config  = { 'silent'  : True } ,
            gen_pars    = { 'mean_GG' : nominal_mean , 'sigma_GG' : nominal_sigma  } ,
            fit_pars    = { 'mean_FG' : nominal_mean , 'sigma_FG' : nominal_sigma } ,
            more_vars   = more_vars , 
            silent      = True      , 
            progress    = True       )
        
        
    ## make histos
        
    h_mean       = ROOT.TH1F ( 'h1' , 'mean  of Gauss ' , 100 ,  0                   , nominal_mean  * 2 )   
    h_sigma      = ROOT.TH1F ( 'h2' , 'sigma of Gauss'  , 100 ,  0.5 * nominal_sigma , nominal_sigma * 2 )

    for r in results ['mean_FG'  ] : h_mean .Fill ( r ) 
    for r in results ['sigma_FG' ] : h_sigma.Fill ( r )

    for h in ( h_mean , h_sigma ) :
        with use_canvas ( 'test_toys2 %s' % h.title  , wait = 1 ) : h.draw()


# ==============================================================================
## Perform toy-study for Jackknife 
def test_jackknife ( ) :
    """Perform toys-study for Jackknife
    """

    logger = getLogger ( 'test_parallel_jackknife' ) 
    
    with timing ( 'Jackknife analysis' , logger = logger )  :        
        model.S = 1000
        model.B = 100
        
        data    = model.generate ( 1100 ) 
        
        Toys.make_jackknife ( pdf        = model ,
                              data       = data  ,
                              fit_config = { 'silent' : True } ,                     
                              progress   = True  ,
                              silent     = True  ) 

# ==============================================================================
## Perform toy-study for Bootstrap 
def test_bootstrap ( ) :
    """Perform toys-study for Bootstrap
    """

    logger = getLogger ( 'test_parallel_bootstrap' ) 

    with timing ( 'Bootstrap analysis' , logger = logger )  :        

        model.S = 1000
        model.B = 100
        
        data    = model.generate ( 1100 ) 
        
        Toys.make_bootstrap ( pdf        = model ,
                              size       = 100   , 
                              data       = data  ,
                              fit_config = { 'silent' : True } ,                     
                              progress   = True  ,
                              silent     = True  ) 
        
# =============================================================================
## Perform toy-study for significance of the signal 
#  - generate <code>nToys</code> pseudoexperiments using background-only hypothesis 
#  - fit teach experiment with "signal+background" hypothesis
#  - store  fit results
#  - fill distributions for fit results
def test_significance ( ) :
    """Perform toy-study for significance of the signal 
    - generate `nToys` pseudoexperiments using background-only hypothesis 
    - fit each experiment with signal+background hypothesis
    - store  fit results
    - fill distributions for fit results
    """
    
    logger = getLogger ( 'test_significance' )

    with timing ( 'Significance analysis' , logger = logger )  :        
        
        ## only background hypothesis
        bkg_only = Models.Bkg_pdf    ( "BKG" , xvar = mass , power = 0 , tau = 0      )        
        signal   = Models.Gauss_pdf  ( 'S'   , xvar = mass , mean = 0.5 , sigma = 0.1 )
        signal.mean .fix ( 0.4 ) 
        signal.sigma.fix ( 0.1 ) 
        
        ## signal + background hypothesis
        model    = Models.Fit1D      ( signal = signal , background = 1 )
        model.background.tau.fix ( 0 )
        
        results , stats = Toys.make_toys2 (
            gen_pdf     = bkg_only  ,
            fit_pdf     = model     ,
            nToys       = 1000      ,
            data        = [ mass ]  , 
            gen_config  = { 'nEvents'  : 100 , 'sample' : True } ,
            fit_config  = { 'silent'   : True } ,
            gen_pars    = { 'tau_BKG'  : 0.   } , ## initial values for generation 
            fit_pars    = { 'B' : 100         ,
                            'S' : 10          ,
                            'phi0_Bkg_S': 0.0 } , ## initial fit values for parameters 
            silent      = True , 
            progress    = True )
        
    h_S      = ROOT.TH1F ( hID() , '#S' , 60 ,  0 , 60 )
        
    for r in results ['S'  ] : h_S .Fill ( r )
    
    for h in ( h_S ,  ) :
        with use_canvas ( 'test_significance  %s' % h.title  , wait = 1 ) : h.draw()
        

# =============================================================================
if '__main__' == __name__ :
    
    test_toys         ( ) 
    test_toys2        ( )
    test_jackknife    ( )
    test_bootstrap    ( )    
    test_significance ( ) 
    

# =============================================================================
##                                                                      The END 
# =============================================================================
