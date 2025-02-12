#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_resolutions3.py
# Test module for ostap/fitting/resolution3.py
# - It demosntaret simfit for data and mc with resolution
# ============================================================================= 
"""Test module for ostap/fitting/resolution3.py
- It tests various asymmetric resolution models
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.pyrouts       import Ostap
from   ostap.fitting.simfit     import combined_data 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait, batch_env  
import ostap.fitting.models     as     Models 
import ROOT, random
# ============================================================================= 
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ :
    logger = getLogger ( 'ostap.test_fitting_resolution3' )
else : 
    logger = getLogger ( __name__ )
# ============================================================================= 
## set batch from environment 
batch_env ( logger )
# =============================================================================


GeV = 1.0
MeV = 0.001 * GeV

m_phi = 1.019 * GeV
g_phi = 4.3   * MeV
m_K   = 497   * MeV 

dm   = ROOT.RooRealVar('dm'  ,'delta-mass'  , -25 * MeV  , 25    * MeV         ) 
mass = ROOT.RooRealVar('dmKK','mass(KK)'    ,  2 * m_K  , 2 * m_K + 100 * MeV  )

# ===============================================================================
## 1a:  prepare MC dataset
# ===============================================================================

## resolution model (slightly asymmetric) 
gen_reso = Models.ResoCB2 ( 'R0'             ,
                            xvar  = dm       ,
                            sigma = 1  * MeV ,  
                            alpha = 1.5      ,
                            n     = 10       ,
                            kappaN = 0.05    , ## small asymmetry 
                            kappaA = 0.05    ) ## small asymmetry 

## MC - dataset
ds_MC = gen_reso.generate ( 100000 ) 

# ===============================================================================
## 1b: prepare DATA dataset
# ===============================================================================
gen_reso_data = Models.ResoCB2 ( 'R'                  ,
                                 xvar   = mass        ,
                                 sigma  = 1.15  * MeV ,  ## 15% wider 
                                 alpha  = 1.5         ,  ## slightly different alpha 
                                 n      = 10          ,  ## slightly different N 
                                 kappaN = 0.05        ,  ## small asymmetry 
                                 kappaA = 0.05        )  ## small asymmetry

ff      = Ostap.Math.FormFactors.BlattWeisskopf ( 1 , 3.5/GeV ) 
phi     = Ostap.Math.BreitWigner ( m_phi  , g_phi , m_K , m_K , 1 , ff ) 
gen_phi = Models.BreitWigner_pdf ( 'BW0' , xvar = mass ,
                                   breitwigner  = phi  ,
                                   m0           = ( m_phi  , m_phi - 10 * MeV , m_phi + 10 * MeV ) ,
                                   gamma        = ( g_phi  , 1 * MeV          , g_phi + 10 * MeV ) )

gen_cnv_conf = { 'resolution' : gen_reso_data  , 
                 'nbins'      : 5000           ,
                 'buffer'     : 0.25           , 
                 'bufstrat'   : 2              } 

## Signals:
signal      = Models.Convolution_pdf   ( gen_phi , **gen_cnv_conf )
ps          = Ostap.Math.PhaseSpace2   ( m_K , m_K ) 
bkg         = Models.PSLeftExpoPol_pdf ( 'B0' , xvar = mass , phasespace = ps , power = 1 , tau = 0 , scale = 1 )

## final model for generation of DATA dataset 
model       = Models.Fit1D( signal = signal , background = bkg , suffix = '_GEN' )  
model.S     = 1000
model.B     =  100

#
ds_DATA     = model.generate ( 50000 ) 


# ===============================================================================
## 2a: Create fit resolution model (MC)
# ===============================================================================
## resolution model (slightly asymmetric) 
reso_MC = Models.ResoCB2 ( 'RMC'                                      ,
                           xvar  = dm                                 ,
                           sigma =  ( 1 * MeV , 0.1 * MeV , 5  *MeV ) ,  
                           alpha =  ( 1.5  , 0.5   , 5.0  )           ,
                           n     =  ( 10   , 0     , 50   )           ,
                           kappaN = ( 0.05 , -0.95 , 0.95 )           , ## small asymmetry 
                           kappaA = ( 0.05 , -0.95 , 0.95 )           ) ## small asymmetry 
model_MC = Models.Fit1D ( signal = reso_MC , background  = None , suffix = '_MC' )

## fit MC data with resolution models
reso_MC.kappaN.fix()
reso_MC.kappaA.fix()
rMC , fMC = model_MC.fitTo ( ds_MC , draw = True , silent = True , nbins = 100 )
reso_MC.kappaN.release ()
reso_MC.kappaA.release ()
rMC , fMC = model_MC.fitTo ( ds_MC , draw = True , silent = True , nbins = 100 )
title = 'Fit to MC dataset to get MC resoluion'
logger.info ( '%s\n%s' % ( title , rMC.table ( title = title , prefix = '# ' ) ) ) 


# ===============================================================================
## 2b: Create fit resolution model (DATA)
# ===============================================================================
reso_DATA  = Models.ResoCB2 ( 'RDATA'                     ,
                              xvar  = mass                ,
                              sigma = reso_MC.sigma       ,  
                              alpha = reso_MC.alpha       ,
                              fudge = ( 1.1 , 0.1 , 5.0 ) , ## NB: fudge-factor  
                              n     =  reso_MC.n          ,
                              kappaN = reso_MC.kappaN     , ## small asymmetry 
                              kappaA = reso_MC.kappaA     ) ## small asymmetry 

phi_DATA = Models.BreitWigner_pdf ( 'BWDATA' , xvar = mass ,
                                    breitwigner  = phi  ,
                                    m0           = ( m_phi  , m_phi - 10 * MeV , m_phi + 10 * MeV ) ,
                                    gamma        = ( g_phi  , 1 * MeV          , g_phi + 10 * MeV ) )
phi_DATA.gamma.fix ( g_phi )
 
cnv_conf = { 'resolution' : reso_DATA , 
             'nbins'      : 5000      ,
             'buffer'     : 0.25      , 
             'bufstrat'   : 2         } 

## Signals:
signal      = Models.Convolution_pdf   ( phi_DATA , **cnv_conf )

## final model for generation of DATA dataset 
model_DATA       = Models.Fit1D( signal = signal , background = bkg , suffix = '_DATA' )  
model_DATA.S     = 0.9 * len ( ds_DATA ) 
model_DATA.B     = 0.1 * len ( ds_DATA )

phi_DATA.m0   .fix ( m_phi )
reso_MC.sigma .fix()

## ## fit data 
## rDATA , fDATA = model_DATA.fitTo ( ds_DATA , draw = True , silent = True , nbins = 100 )
## rDATA , fDATA = model_DATA.fitTo ( ds_DATA , draw = True , silent = True , nbins = 100 )
## title = 'Fit to DATA with floating resolution '
## logger.info ( '%s\n%s' % ( title , rDATA.table ( title = title , prefix = '# ' ) ) )


# ===============================================================================
## 3a: Construction for simultaneous fit 
# ===============================================================================

## combine data
sample  = ROOT.RooCategory ('sample','sample'  , 'DATA' , 'MC' )

## combine datasets
vars    = ROOT.RooArgSet ( mass , dm )
dataset = combined_data  ( sample , vars , { 'DATA' : ds_DATA , 'MC' : ds_MC } )

## combine PDFs
model_sim  = Models.SimFit ( sample , { 'DATA' : model_DATA , 'MC' : model_MC } , name = 'X' )

## fit it!
rSIM , fSIM = model_sim.fitTo ( dataset , silent = True )
rSIM , fSIM = model_sim.fitTo ( dataset , silent = True )
title = 'Simultaneous fit(1) to DATA&MC'
logger.info ( '%s\n%s' % ( title , rSIM.table ( title = title , prefix = '# ' ) ) ) 


phi_DATA.m0   .release() 
reso_MC.sigma .release() 

## fit it!
rSIM , fSIM = model_sim.fitTo ( dataset , silent = True )
rSIM , fSIM = model_sim.fitTo ( dataset , silent = True )
title = 'Simultaneous fit(2) to DATA&MC'
logger.info ( '%s\n%s' % ( title , rSIM.table ( title = title , prefix = '# ' ) ) ) 

with use_canvas ( 'test_resoltuion3' ) : 
    with wait ( 3 ) : fMC   = model_sim.draw ( 'MC'   , dataset , nbins = 100 )
    with wait ( 3 ) : fDATA = model_sim.draw ( 'DATA' , dataset , nbins = 100 )


# =============================================================================
##                                                                      The END 
# =============================================================================
