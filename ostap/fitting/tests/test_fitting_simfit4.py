#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_simfit4.py
# Test module for ostap/fitting/simfit.py
# - It tests the most simple "Simultaneous fit"
#
#  Simultaneous fit of 1D and 2D-distributions 
# ============================================================================= 
""" Test module for ostap/fitting/simfit.py
- It tests the most simple ``Simultaneous fit''

Simultaneous fit of two 1D distributions  
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.fitting.fithelpers import FitHelper 
from   ostap.core.core          import dsID, rooSilent
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.ranges       import vrange
from   ostap.utils.root_utils   import batch_env
import ostap.fitting.models     as     Models
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_simfit4' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
#
## make simple test mass 
mass1    = ROOT.RooRealVar ( 'test_mass1' , 'Some test mass' ,  0  , 10 )
mass2    = ROOT.RooRealVar ( 'test_mass2' , 'Some test mass' ,  5  , 15 )

## book very simple data set:
varset1  = ROOT.RooArgSet  ( mass1 )
dataset1 = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset1 )  

## book very simple data set:
varset2  = ROOT.RooArgSet  ( mass2 )
dataset2 = ROOT.RooDataSet ( dsID() , 'Test Data set-2' , varset2 )  

## high statistic, low-background "control channel"
mean1  = 5.0
sigma1 = 0.5
NS1    = 20000
NB1    =  2000

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
        
## low statistic, high-background "signal channel"
NS2     =    500
NB2     =  10000     
mean2   =  mean1  + 5.0 
sigma2  =  sigma1 * 2.0  

for i in range(NS2)  :
    v2 = random.gauss ( mean2 , sigma2 )
    if v2 in mass2 :
        mass2.setVal ( v2 )
        dataset2.add ( varset2 )

for i in range (NB2 ) :
    v2 = random.uniform (  0 , 20  )
    if v2 in mass2 : 
        mass2.setVal ( v2 )
        dataset2.add ( varset2 )

graphs = [] 
# =============================================================================
def test_simfit4() :

    logger = getLogger( 'test_simfit4' ) 
    # =========================================================================
    VARS  = FitHelper ()
    
    ## MC/DATA ratio for signal widths
    Rs   = ROOT.RooRealVar ( 'Rs' , 'MC/DATA ratio for signal widths' , 1.0 , 0.2 , 3.0 )
    
    ## sigma for normalization peak: sigma from MC multiplied by a factor 
    sigma_N  = VARS.vars_multiply ( Rs , sigma1 , name = 'sigma_N' )
    
    ## sigma for signal peak: sigma from MC multiplied by the same factor 
    sigma_S  = VARS.vars_multiply ( Rs , sigma2 , name = 'sigma_S' )
    
    ## Normalization peak:
    signal_N = Models.Gauss_pdf ( 'GN'               ,
                                  xvar  = mass1      ,
                                  mean  = ( 3 , 8 )  ,
                                  sigma = sigma_N    )
    
    ## mean of low-statistics signal is constrained to the mean of high statistic signal
    mean_S = VARS.vars_sum ( signal_N.mean , mean2 - mean1 ) 
    
    ## low statistic signal :
    signal_S = Models.Gauss_pdf ( 'GS'           ,
                                  xvar  = mass2   ,
                                  mean  = mean_S  ,
                                  sigma = sigma_S )
    
    
    # =========================================================================
    ## 3rd order positive polynomial 
    bkg_N = Models.make_bkg ( -3 , 'Bkg_N' , mass1  )
    ## model for fitting of normalization channel 
    model_N = Models.Fit1D ( suffix     = 'N'       ,
                             signal     =  signal_N ,
                             background =  bkg_N    )
    model_N.S = NS1 
    model_N.B = NB1 
    
    
    ## 3rd order positive polynomial 
    bkg_S = Models.make_bkg ( -3 , 'Bkg_S' , mass2  )
    ## model for fitting of low-statistic signal channel 
    model_S = Models.Fit1D ( suffix     = 'S'       ,
                             signal     =  signal_S ,
                             background =  bkg_S    )
    
    model_S.S = NS2 
    model_S.B = NB2 

    # =========================================================================
    ## make fit for thew normalization channel only 
    # =========================================================================
    with use_canvas ( 'test_simfit4: norm' ) : 
        rN , fN = model_N.fitTo ( dataset1 , draw = None ,              silent = True )
        rN , fN = model_N.fitTo ( dataset1 , draw = None ,              silent = True )
        rN , fN = model_N.fitTo ( dataset1 , draw = True , nbins = 50 , silent = True )
        graphs.append ( fN )
    
    title = 'Fit to high-statistic normalisation sample'
    logger.info ( 'Fit results for normalization sample only: %s' % rN.table ( title = title ,
                                                                               prefix = '# ' ) )
    
    # =========================================================================
    ## make fit for the low-statistic signal channel only 
    # =========================================================================
    with use_canvas ( 'test_simfit4: signal' ) : 
        rS , fS = model_S.fitTo ( dataset2 , draw = None ,              silent = True )
        rS , fS = model_S.fitTo ( dataset2 , draw = None ,              silent = True )
        rS , fS = model_S.fitTo ( dataset2 , draw = True , nbins = 50 , silent = True )
        graphs.append ( fS )
        
    with use_canvas ( 'test_simfit4: profile ' ) :
        grs = model_S.graph_profile ( 'S_S' , vrange ( 50 , 1200 , 100 ) , dataset2 , draw = True )
        grs.draw('apl')
        graphs.append ( grs )
        
    title = 'Fit to low-statistics sample'
    logger.info ( 'Fit results for low-statistic signal only:\n%s' % rS.table ( title  = title ,
                                                                                prefix = '# '  )  )
    
    # =========================================================================
    ## combine data for simultaneous  fit
    # =========================================================================
    
    sample  = ROOT.RooCategory ('sample','sample'  , 'S' , 'N' )
    
    ## combine datasets
    from ostap.fitting.simfit import combined_data 
    vars    = ROOT.RooArgSet ( mass1 , mass2 )
    dataset = combined_data  ( sample , vars , { 'N' : dataset1 , 'S' : dataset2 } )
    
    ## combine PDFs
    model_sim  = Models.SimFit (
        sample , { 'S' : model_S  , 'N' : model_N } , name = 'X'
        )
    
    # =========================================================================
    rC , fC = model_sim.fitTo ( dataset , silent = True )
    rC , fC = model_sim.fitTo ( dataset , silent = True )
    rC , fC = model_sim.fitTo ( dataset , silent = True )
    rC , fC = model_sim.fitTo ( dataset , silent = True )
    
    with use_canvas ( 'test_simfit4: N' ) : fN  = model_sim.draw ( 'N'   , dataset , nbins = 50 )
    with use_canvas ( 'test_simfit4: S' ) : fS  = model_sim.draw ( 'S'   , dataset , nbins = 50 )
    
    with use_canvas ( 'test_simfit4: profile-sim' , wait = 3 ) :       
       grs = model_sim.graph_profile ( 'S_S' , vrange ( 100 , 1000 , 50 ) , dataset , draw = True )
       grs.draw('apl')
       graphs.append ( grs )
       title = 'Simultaneous fit'
       logger.info ( 'Combined fit  results are:\n%s ' % rC.table ( title  = title ,
                                                                    prefix = '#'   )  )
    # =========================================================================
    ## test creation of dataset
    # =========================================================================
    ds_gen = model_sim.generate ( nEvents = { 'N' : len ( dataset1 ) ,
                                              'S' : len ( dataset2 ) } ,
                                  varset  = vars  )

    rg , f = model_sim.fitTo ( ds_gen , silent = True )
    rg , f = model_sim.fitTo ( ds_gen , silent = True )
    
    title = 'Results of simultaneous fit to generated dataset'
    logger.info ( '%s\n%s' % ( title , rg.table ( title = title , prefix = '# ' ) ) )

    
# =============================================================================
if '__main__' == __name__ :

    with timing ( "simfit-4" , logger ) : 
        test_simfit4() 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
