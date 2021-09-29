#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_jackknife.py
# Test module for ostap/fitting/toys.py
# - run Jackknife  analysis
# ============================================================================= 
""" Test module for ostap/fitting/toys.py
-  run Jackknife  analysis
"""
# ============================================================================= 
## from   __future__        import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ROOT, time, random, math
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_jackknife' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
from   ostap.core.pyrouts       import hID, VE 
import ostap.fitting.models     as     Models
import ostap.fitting.toys       as     Toys
import ostap.histos.histos
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
# =============================================================================
mass        = ROOT.RooRealVar  ( 'mass' , '', 0 , 1 )  
gauss       = Models.Gauss_pdf ( 'G'    , xvar = mass )
gauss.mean  = 0.4
gauss.sigma = 0.1
## model       = Models.Fit1D ( signal = gauss ) 
## model.S     = 100
## model.B     = 100 
model = gauss

# ==============================================================================
## Perform jackknife study for possible fit bias and correct uncertainty evaluation 
def test_jackknife ( ) :
    """Perform toys-study for possible fit bias and correct uncertainty evaluation
    - generate `nToys` pseudoexperiments with some PDF `pdf`
    - fit teach experiment with the same PDF
    - store  fit results
    - calculate statistics of pulls
    - fill distributions of fit results
    - fill distributions of pulls 
    """

    logger = getLogger ( 'test_jackknife' )

    N = 400 
    dataset = model.generate ( N , sample = False )
    
    more_vars   = { 'vm' : lambda  r, *_ : ( r.mean_G - 0.4 ) / 0.1       ,
                    'vs' : lambda  r, *_ :   r.sigma_G        / 0.1 - 1   ,
                    'vr' : lambda  r, *_ :   r.sigma_G * 1    / r.mean_G  } 
    
     ## prefit the whole dataset
    with use_canvas ( 'test_jackknife' ) : 
        res , f = model.fitTo ( dataset , draw = True , nbins = 100 , silent = True , refit = 5 )

    ## start Jackknife process 
    results , stats = Toys.make_jackknife  (
        pdf         = model      ,
        data        = dataset    , 
        fit_config  = { 'silent' : True , 'refit'   : 5   } ,
        fit_pars    = { 'mean_G' : 0.4  , 'sigma_G' : 0.1 } ,
        more_vars   = more_vars  , 
        silent      = True , 
        progress    = True )

    ## fit the whole sample 
    with use_canvas ( 'test_jackknife' ) : 
        res , f = model.fitTo ( dataset , draw = True  , nbins = 100 , silent = True , refit = 5 )
    
    logger.info  ('Fit results:\n%s' % res.table ( title     = 'Fit results' ,
                                                   prefix    = '# '          ,
                                                   more_vars = more_vars     )  ) 
    ## print the final table
    Toys.print_jackknife ( res   ,
                           stats ,
                           morevars = dict ( ( k , more_vars [ k ]( res , model ) ) for k in more_vars ) ,
                           logger   = logger )

    time.sleep ( 2 ) 


# =============================================================================
if '__main__' == __name__ :

    test_jackknife ( ) 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
