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
from   ostap.core.pyrouts   import hID 
import ostap.fitting.models as     Models
import ostap.fitting.toys   as     Toys

mass      = ROOT.RooRealVar ( 'mass' , '', 0 , 1 )  
gen_gauss = Models.Gauss_pdf ( 'GG' , xvar = mass )
fit_gauss = Models.Gauss_pdf ( 'FG' , xvar = mass )
gen_gauss.mean  = 0.4
gen_gauss.sigma = 0.1

def test_toys ( ) :
    
    results , stats = Toys.make_toys  (
        pdf         = gen_gauss ,
        nToys       = 1000      ,
        data        = [ mass ]  , 
        gen_config  = { 'nEvents' : 200  } ,
        fit_config  = { 'silent'  : True } ,
        init_pars   = { 'mean_GG' : 0.4 , 'sigma_GG' : 0.1 } ,
        silent      = True , 
        progress    = True )

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
def test_toys2 ( ) :
                
    results , stats = Toys.make_toys2 (
        gen_pdf     = gen_gauss ,
        fit_pdf     = fit_gauss ,
        nToys       = 1000      ,
        data        = [ mass ]  , 
        gen_config  = { 'nEvents' : 200  } ,
        fit_config  = { 'silent'  : True } ,
        gen_pars    = { 'mean_GG' : 0.4 , 'sigma_GG' : 0.1 } ,
        fit_pars    = { 'mean_GF' : 0.4 , 'sigma_GF' : 0.1 } ,
        silent      = True , 
        progress    = True )
    
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
if '__main__' == __name__ :

    test_toys  () 
    test_toys2 () 


# =============================================================================
##                                                                      The END 
# =============================================================================
