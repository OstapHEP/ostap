#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_roostats4.py
# Test module for RooStats 
# ============================================================================= 
""" Test module for RooStats 
"""
# ============================================================================= 
from   __future__               import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   builtins                 import range
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models
from   ostap.fitting.variables  import SETVAR, FIXVAR  
from   ostap.core.core          import cpp, VE, dsID, hID , rooSilent, Ostap 
from   ostap.utils.timing       import timing
from   ostap.utils.utils        import vrange
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait
import ostap.logger.table       as     T
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_roostats4' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================


mass      = ROOT.RooRealVar   ('mass','mass-variable'        , 0 , 10 )

signal1   = Models.Gauss_pdf    ( 'G1',
                                  xvar  = mass                ,
                                  mean  = ( 4   , 1.0 , 9.0 ) ,
                                  sigma = ( 0.5 , 0.1 , 1.5 ) )

signal2   = Models.Gauss_pdf   ( 'G2',
                                  xvar  = mass                 ,
                                 mean  = signal1.vars_add( signal1.mean , 2.0 ) ,                                
                                 sigma = signal1.sigma       )

## S1/S2 
f1         = ROOT.RooRealVar ( 'f1'        , 'N1/N2 ratio' , 0.05    , 0 , 1 )
## Ratio of efficinecies
eff_ratio  = ROOT.RooRealVar ( 'eff_ratio' , 'ratio of efficiencies' ,  1.0 , 0.2 , 3.0 )                             
## constraint for efficieency ratio 
eff_constraint = signal1.soft_constraint ( eff_ratio , VE ( 1.0 , 0.05**2 ) )
## S2 
N2         = ROOT.RooRealVar ( 'S2'        , 'S2-signal'   , 1000 , 0 , 100000 )
## S1
f1r        = signal1.vars_multiply ( f1   , eff_ratio )
N1         = signal1.vars_multiply ( f1r , N2 ) 
##
model   = Models.Fit1D ( signals = ( signal1 , signal2 ) , background = 'e-' , S = ( N1 , N2 ) , suffix = 'RD' )
model.B = 1000
model.background.tau = -0.25 

mcsignal = Models.Gauss_pdf    ( 'MC',
                                 xvar  = mass           ,
                                 mean  = signal1.vars_add( signal1.mean , 1.0 ) ,                                
                                 sigma = signal1.sigma  )



        
# ==========================================================================
## data sets
# ==========================================================================
N_DATA   = 200
rd_data  = model   .generate ( N_DATA )

N_MC     = 10 * N_DATA 
mc_data  = mcsignal.generate ( N_MC )

# ===========================================================================
## Create sim-fit data and model 
# ===========================================================================
## combine data
sample  = ROOT.RooCategory ('sample','sample'  , 'Data' , 'MC' )

## combine datasets
from ostap.fitting.simfit import combined_data 
vars    = ROOT.RooArgSet ( mass )
dataset = combined_data  ( sample , vars , { 'Data' : rd_data , 'MC' : mc_data } )

## combine PDFs
model_sim  = Models.SimFit (
    sample , { 'Data' : model  , 'MC' : mcsignal } , name = 'X'
    )
  
# ===========================================================================
## Fit MC data only
# ===========================================================================
rMC , fMC = mcsignal.fitTo ( mc_data , silent = True )
rMC , fMC = mcsignal.fitTo ( mc_data , silent = True )
with use_canvas ( 'MC-sample' ) :  fMC = model.draw  ( rd_data , nbins = 100 )
title = 'Results of the MC-fit'
logger.info ( '%s\n%s' % ( title , rMC.table ( title = title , prefix = '# ' ) ) )

# ===========================================================================
## Fit RD data only
# ===========================================================================
rRD , fRF = model.fitTo ( rd_data , silent = True , constraints = eff_constraint )
rRD , fRD = model.fitTo ( rd_data , silent = True , constraints = eff_constraint )
with use_canvas ( 'RD-sample' ) :  fRD = model.draw  ( rd_data , nbins = 100 )
title = 'Results of the RD-fit'
logger.info ( '%s\n%s' % ( title , rRD.table ( title = title , prefix = '# ' ) ) )


## # ===========================================================================
## ## SimFit 
## # ===========================================================================
## rSIM , f = model_sim.fitTo ( dataset , silent = True )
## rSIM , f = model_sim.fitTo ( dataset , silent = True )
## with use_canvas ( 'SimFit: RD-sample' ) :  fRDsim = model_sim.draw  ( 'Data' , dataset , nbins = 100 ) 
## with use_canvas ( 'SimFit: MC-sample' ) :  fMCsim = model_sim.draw  ( 'MC'   , dataset , nbins = 100 ) 
## title = 'Results of Sim-fit'
## logger.info ( '%s\n%s' % ( title , rSIM.table ( title = title , prefix = '# ' ) ) )


# ============================================================================-
## Get the upper limit limit for small signal at fixed mass
#  - Asymptotic Calculator is used 
def test_limit_ac_1 () :
    """ Get the upper limit at given point for small signal at fixed mass
    - Asymptotic Calculator is used 
    """
    
    logger = getLogger("test_limit_ac_1")
    
    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )
    
    the_model = model  .clone ( name = 'M1' )
    the_data  = rd_data.clone()

    constraints = eff_constraint ,
    
    logger.info ( 'Dataset is\n%s' % the_data.table ( prefix = '# ' ) ) 
    rS , frame = the_model.fitTo ( the_data , nbins = 100 , constraints = eff_constraint )
    with use_canvas ( 'test_limit_ac_1:data' ) : the_model.draw ( the_data , nbins = 100 )
    title = 'Results of the fit'
    logger.info ( '%s\n%s' % ( title , rS.table ( title = title , prefix = '# ' ) ) )
    
    
    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf         = the_model       ,
                             poi         = f1              , ## parameter of interest 
                             dataset     = the_data        ,
                             constraints = constraints     ,
                             name        = 'S+B'           ,
                             snapshot    = rS              )
    
    with FIXVAR ( f1 ) :
        f1.setVal ( 0 ) 
        rB , _ = the_model.fitTo ( the_data , constraints = eff_constraint , silent = True )
        
    ## create ModelConfig  for 'B-only' model
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = f1                 , ## parameter of interest 
                             dataset     = the_data           ,
                             constraints = constraints        ,
                             workspace   = model_sb.workspace , 
                             name        = 'B-only'           ,
                             snapshot    = rB                 )
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b .name , model_b .table ( prefix = '# ' ) ) )
    
    from   ostap.core.core import rootError
    
    with timing ( "Using Asymptotic Calculator" , logger = logger ) , rooSilent () , rootError() :
        ## create the calculator 
        ac  = AsymptoticCalculator ( model_b              ,
                                     model_sb             ,
                                     dataset   = the_data ,
                                     asimov    = False    )
        ac.calculator.SetOneSided ( True ) 
        ## create Hypo Test inverter 
        hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
        
        ## make a scan 
        hti .scan_with_progress ( vrange ( 0 , 0.5 , 100 )  ) ## scan it!
        
    ## visualize the scan results 
    with use_canvas ( 'test_limit_ac_1: HypoTestInverter plot (asymptotic)' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.3f' % hti.upper_limit )

# ============================================================================-
## Get the upper limit limit for small signal at fixed mass
#  - Asymptotic Calculator is used 
def test_limit_ac_2 () :
    """Get the upper limit at given point for small signal at fixed mass
    - Asymptotic Calculator is used 
    """
    
    logger = getLogger("test_limit_ac_2")
    
    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )
    
    the_model   = model  .clone ( name = 'M2' )
    the_data    = rd_data.clone()

    constraint  = mcsignal.soft_multivar_constraint ( ( 'mean_G1' , 'sigma_G1' ) , rMC ) 

    constraints = constraint , eff_constraint 

    logger.info ( 'Dataset is\n%s' % the_data.table ( prefix = '# ' ) ) 
    rS , frame = the_model.fitTo ( the_data , nbins = 100 , constraints = constraints )
    
    with use_canvas ( 'test_limit_ac_2:data' ) : the_model.draw ( the_data , nbins = 100 )
    title = 'Results of the fit'
    logger.info ( '%s\n%s' % ( title , rS.table ( title = title , prefix = '# ' ) ) )
        
    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf         = the_model   ,
                             poi         = f1          , ## parameter of interest 
                             dataset     = the_data    ,
                             constraints = constraints ,
                             name        = 'S+B'       ,
                             snapshot    = rS          )
    
    
    with FIXVAR ( f1 ) :
        f1.setVal ( 0 ) 
        rB , _ = the_model.fitTo ( the_data , constraints = constraints , silent = True )
 
    ## create ModelConfig  for 'B-only' model
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = f1                 , ## parameter of interest 
                             dataset     = the_data           ,
                             workspace   = model_sb.workspace ,
                             constraints = constraints        ,                               
                             name        = 'B-only'           ,
                             snapshot    = rB                 )
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b .name , model_b .table ( prefix = '# ' ) ) )
    
    from   ostap.core.core import rootError
    
    with timing ( "Using Asymptotic Calculator" , logger = logger ) , rooSilent () , rootError() :
        ## create the calculator 
        ac  = AsymptoticCalculator ( model_b              ,
                                     model_sb             ,
                                     dataset   = the_data ,
                                     asimov    = False    ) 
        ac.calculator.SetOneSided ( True ) 

        ## create Hypo Test inverter 
        hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
        
        ## make a scan 
        hti .scan_with_progress ( vrange ( 0 , 0.5 , 100 )  ) ## scan it!
        
    ## visualize the scan results 
    with use_canvas ( 'test_limit_ac_2: HypoTestInverter plot (asymptotic)' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.3f' % hti.upper_limit )

# ============================================================================-
## Get the upper limit limit for small signal at fixed mass
#  - Frequentists Calculator is used 
def test_limit_fc_1 () :
    """ Get the upper limit at given point for small signal at fixed mass
    - Friquentists Calculator is used 
    """
    
    logger = getLogger("test_limit_fc_1")
    
    logger.info ( "Test Point limits with RooStats (Frequentists Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             FrequentistCalculator  ,
                                             HypoTestInverter      )
    
    the_model = model  .clone ( name = 'M3' )
    the_data  = rd_data.clone()

    constraint  = mcsignal.soft_multivar_constraint ( ( 'mean_G1' , 'sigma_G1' ) , rMC ) 
    constraints = constraint , eff_constraint 
    
    logger.info ( 'Dataset is\n%s' % the_data.table ( prefix = '# ' ) ) 
    rS , frame = the_model.fitTo ( the_data , nbins = 100 , constraints = constraints , silent = True )
    with use_canvas ( 'test_limit_ac_2:data' ) : the_model.draw ( the_data , nbins = 100 )
    title = 'Results of the fit'
    logger.info ( '%s\n%s' % ( title , rS.table ( title = title , prefix = '# ' ) ) )
    
    
    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf         = the_model   ,
                             poi         = f1          , ## parameter of interest 
                             dataset     = the_data    ,
                             constraints = constraints ,  
                             name        = 'S+B'       ,
                             ## snapshot    = f1          )
                             snapshot    = rS          )

    with FIXVAR ( f1 ) :
        f1.setVal ( 0 ) 
        rB , _ = the_model.fitTo ( the_data , constraints = constraints , silent = True )

    ## create ModelConfig  for 'B-only' model
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = f1                 , ## parameter of interest 
                             dataset     = the_data           ,
                             workspace   = model_sb.workspace ,
                             constraints = constraints        ,                               
                             name        = 'B-only'           ,
                             ## snapshot    = rB                 )
                             snapshot    = f1                 )

    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b .name , model_b .table ( prefix = '# ' ) ) )
    
    from   ostap.core.core import rootError
    
    with timing ( "Using Frequentists Calculator" , logger = logger ) , rooSilent () , rootError() :

        ## create the calculator 
        fc  = FrequentistCalculator ( model_b               ,
                                      model_sb              ,
                                      dataset    = the_data ,
                                      ntoys_null = 100      ,
                                      ntoys_alt  = 100      ) 

        ## create Hypo Test inverter 
        hti = HypoTestInverter ( fc ,  0.90 , use_CLs = True , verbose = False )
        
        ## make a scan 
        hti .scan_with_progress ( vrange ( 1.e-3 , 0.4 , 20 )  ) ## scan it!
        
    ## visualize the scan results 
    with use_canvas ( 'test_limit_fc_1: HypoTestInverter plot (frequentist)' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (frequentist)  = %.3f' % hti.upper_limit )

# =============================================================================
if '__main__' == __name__ :

    test_limit_ac_1 () 
    test_limit_ac_2 ()    
    test_limit_fc_1 () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
