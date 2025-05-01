#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_simfit7.py
# Test module for ostap/fitting/simfit.py
# - It tests the most simple "Simultaneous fit" versus constrained fit 
# ============================================================================= 
""" Test module for ostap/fitting/simfit.py
- It tests the most simple `Simultaneous fit' versus constrained fit 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.core          import dsID, rooSilent 
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.root_utils   import batch_env 
from   ostap.fitting.simfit     import combined_data
import ostap.fitting.models     as     Models 
import ostap.io.zipshelve       as     DBASE
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_simfit7' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
#
## make simple test mass 
mass      = ROOT.RooRealVar ( 'mass'     , 'Some test mass' ,  0 , 20 )

signal_MC = Models.CB2_pdf ( 'MC' ,
                             xvar   = mass                 ,
                             alphaL = ( 1.75 , 1.5 , 2.0 ) ,  
                             alphaR = ( 1.75 , 1.5 , 2.0 ) ,  
                             nL     = ( 5    , 2   , 100 ) ,  
                             nR     = ( 5    , 2   , 100 ) ,  
                             mean   = ( 10 , 3     , 18  )    ,
                             sigma  = ( 1 , 0.5 , 2.0 )    )

## scale factor for resolution 
sfactor     = ROOT.RooRealVar ( 'sfactor' , 'Resolution scale factor' , 1.05 , 0.50 , 1.50 )
## resolution on DATA
sigma_DATA  = signal_MC.vars_multiply ( sfactor , signal_MC.sigma )
signal_DATA = Models.CB2_pdf ( 'DATA' ,
                               xvar   = mass           , 
                               mean   = signal_MC.mean ,
                               alphaL = signal_MC.aL   ,
                               alphaR = signal_MC.aR   ,
                               nL     = signal_MC.nL   ,
                               nR     = signal_MC.nR   ,                               
                               sigma  = sigma_DATA     )

model_MC   = Models.Fit1D ( signal = signal_MC   , background = 'flat' , suffix = 'MC'   )
model_DATA = Models.Fit1D ( signal = signal_DATA , background = -1     , suffix = 'DATA' )

model_MC.S   = 1000
model_MC.B   =   10

model_DATA.S = 1000
model_DATA.B = 1000

ds_MC   = model_MC  .generate ( 100000 )
ds_DATA = model_DATA.generate (  16000 )

# =============================================================================
## (0) (separate) pre-fit MC 
rMC, _  = model_MC  .fitTo ( ds_MC , silent = True )
logger.info  ( 'Pre-fit MC\n%s'   % rMC.table ( prefix = '# ' ) )
## create MC-constraint 
constr_MC  = model_MC.soft_multivar_constraint ( [ 'aL_MC', 'aR_MC' , 'mean_MC'  , 
                                                   'nL_MC', 'nR_MC' , 'sigma_MC' ] , rMC ) 
# =============================================================================
## (1) simultaneous fit of MC and DATA
# =============================================================================
logger.info ( 'Make simultaneous DATA/MC fit to get scale-factor' ) 
category = ROOT.RooCategory ( 'sample' , 'sample' , 'MC' , 'DATA' )
vars     = ROOT.RooArgSet   ( mass  )
cds      = combined_data ( category ,
                           vars     ,  
                           datasets = { 'MC'   : ds_MC   ,
                                        'DATA' : ds_DATA } )
logger.info ( 'Combined dataset/1:\n%s' % cds.table ( prefix = '# ' ) )                            
model_sim  = Models.SimFit (
    category , { 'DATA' : model_DATA ,
                 'MC'   : model_MC   } , name = 'X' )

rSIM , _ = model_sim.fitTo ( cds , silent = True )
rSIM , _ = model_sim.fitTo ( cds , silent = True )
logger.info ( 'Simultaneous fit result:\n%s' % rSIM.table ( prefix = "# " ) )
with use_canvas ( 'test_simfit7: simultaneous fit/DATA' ) :
    fDATA = model_sim.draw ( 'DATA'  , cds )
with use_canvas ( 'test_simfit7: simultaneous fit/MC' ) :
    fMC   = model_sim.draw ( 'MC'    , cds )


# =============================================================================
## (2) Perform DATA fit with MC constraints to get scale-factor
# =============================================================================
logger.info ( 'Perform DATA fit with MC constraints to get scale-factor' ) 

with use_canvas ( 'test_simfit7: constrained fit/DATA' ) :
    rDATA , fDATA = model_DATA.fitTo ( ds_DATA , silent = True , constraints = constr_MC ) 
    rDATA , fDATA = model_DATA.fitTo ( ds_DATA , silent = True , constraints = constr_MC , draw = True ) 
    logger.info ( 'Fit with constraits result:\n%s' % rDATA.table ( prefix = "# " ) )

    
logger.info ( "Scale factor from 'simultaneous' fit: %s" % ( rSIM .sfactor * 1 ) )
logger.info ( "Scale factor from 'constrained'  fit: %s" % ( rDATA.sfactor * 1 ) )


# =========================================================================
## test creation of dataset
# =========================================================================
ds_gen = model_sim.generate ( nEvents = { 'MC'   : len ( ds_MC   ) ,
                                          'DATA' : len ( ds_DATA ) } ,
                              varset  = vars  )

rg , f = model_sim.fitTo ( ds_gen , silent = True )
rg , f = model_sim.fitTo ( ds_gen , silent = True )

title = 'Results of simultaneous fit to generated dataset'
logger.info ( '%s\n%s' % ( title , rg.table ( title = title , prefix = '# ' ) ) )



    
# =============================================================================
##                                                                      The END 
# =============================================================================
