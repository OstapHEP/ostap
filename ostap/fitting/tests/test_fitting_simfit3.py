#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_simfit3.py
# Test module for ostap/fitting/simfit.py
# - It tests the most simple "Simultaneous fit"
#
#  Simultaneous fit of 1D and 2D-distributions 
# ============================================================================= 
""" Test module for ostap/fitting/simfit.py
- It tests the most simple ``Simultaneous fit''

Simultaneous fit of 1D and 2D-distributions 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
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
    logger = getLogger ( 'test_fitting_simfit3' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
#
## make simple test mass 
mass1    = ROOT.RooRealVar ( 'test_mass1' , 'Some test mass' ,  0 , 10 )
mass2    = ROOT.RooRealVar ( 'test_mass2' , 'Some test mass' ,  0 , 20 )
mass3    = ROOT.RooRealVar ( 'test_mass3' , 'Some test mass' ,  0 , 20 )

## book very simple data set:
varset1  = ROOT.RooArgSet  ( mass1 )
dataset1 = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset1 )  

## book very simple data set:
varset2  = ROOT.RooArgSet  ( mass2 , mass3 )
dataset2 = ROOT.RooDataSet ( dsID() , 'Test Data set-2' , varset2 )  

## high statistic, low-background "control channel"
mean1  = 5.0
sigma1 = 0.75
NS1    =  20000
NB1    =   2000

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
NC2     =    500
mean2   = mean1  + 10.0
sigma2  = sigma1 *  2.0

for i in range(NS2)  :
    v2 = random.gauss ( mean2 , sigma2 )
    v3 = random.gauss ( mean2 , sigma2 )
    if v2 in mass2 and v3 in mass3 :
        mass2.setVal ( v2 )
        mass3.setVal ( v3 )
        dataset2.add ( varset2 )

for i in range (NB2 ) :
    v2 = random.uniform (  0 , 20  )
    v3 = random.uniform (  0 , 20  )
    if v2 in mass2 and v3 in mass3 : 
        mass2.setVal ( v2 )
        mass3.setVal ( v3 )
        dataset2.add ( varset2 )

for i in range (NC2) :
     v2 = random.gauss   ( mean2 , sigma2 )
     v3 = random.uniform (  0 , 20  )
     if v2 in mass2 and v3 in mass3 :
         mass2.setVal ( v2 )
         mass3.setVal ( v3 )
         dataset2.add ( varset2 )
     v2 = random.uniform (  0 , 20  )
     v3 = random.gauss   ( mean2 , sigma2 )
     if v2 in mass2 and v3 in mass3 :
         mass2.setVal ( v2 )
         mass3.setVal ( v3 )
         dataset2.add ( varset2 )

         
graphs = []
# =============================================================================
def test_simfit3() : 

    logger = getLogger ( 'test_simfit3' ) 
    ## low statistic, high-background "signal channel"
    signal2  = Models.Gauss_pdf ( 'G2'                  ,
                                  xvar  = mass2         ,
                                  mean  = ( 15  , 20  ) ,
                                  sigma = ( 0.1 , 5.0 ) )
    
    signal3  = Models.Gauss_pdf ( 'G3'                  ,
                                  xvar  = mass3         ,
                                  mean  = signal2.mean  ,
                                  sigma = signal2.sigma ) 
    
    model_S   = Models.Fit2DSym  ( suffix   = 'M2'    ,
                                   yvar     = mass3   , 
                                   signal_x = signal2 ,
                                   signal_y = signal3 ,
                                   bkg_1x   = None    ,
                                   bkg_2x   = None    )
    model_S.BB = NB2
    model_S.SS = NS2
    model_S.SB = NC2  
    model_S.SB.setMin ( -500 ) 
    
    ## fit 2
    rS , fx = model_S.fitTo ( dataset2 , draw = None , nbins = 50 , silent = True )
    rS , fx = model_S.fitTo ( dataset2 , draw = None , nbins = 50 , silent = True )

    title = 'Results of fit to signal sample only'
    logger.info ( '%s\n%s' % ( title , rS.table ( title = title , prefix = '# ' ) ) )
    
    with use_canvas ( 'test_simfit3 : X' ) : 
        rS , fx = model_S.fitTo ( dataset2 , draw = 'X'  , nbins = 50 , silent = True )
    with use_canvas ( 'test_simfit3 : Y' ) : 
        rS , fy = model_S.fitTo ( dataset2 , draw = 'Y'  , nbins = 50 , silent = True )
        graphs.append ( fx )
        graphs.append ( fy )
    with use_canvas ( 'test_simfit3: profile ' ) :
        grs = model_S.graph_profile ( 'SSM2' , vrange ( 100 , 1000 , 50 ) , dataset2 , draw = True )
        grs.draw('apl')
        graphs.append ( grs )
                                       
    # =========================================================================
    ## high statistic, low-background "control/normalization channel"
    mean_N   = signal2.vars_subtract ( signal2.mean  , 10.0 )
    sigma_N  = signal2.vars_multiply ( signal2.sigma ,  0.5 )
    
    signal_N = Models.Gauss_pdf ( 'G1'            ,
                                  xvar  = mass1   ,
                                  mean  = mean_N  ,
                                  sigma = sigma_N )
    
    model_N   = Models.Fit1D ( suffix = 'M1' , signal = signal_N ,  background = -1 )
    model_N.S = NS1
    model_N.B = NB1 
    
    ## fit 1 
    rN , fN = model_N.fitTo ( dataset1 , draw = False , nbins = 50 , silent = True )
    with use_canvas ( 'test_simfit3' , wait = 1 ) :
        rN , fN = model_N.fitTo ( dataset1 , draw = True  , nbins = 50 , silent = True )
        title = 'Results of fit to normalization sample only'
        logger.info ( '%s\n%s' % ( title , rN.table ( title = title , prefix = '# ' ) ) )
        
    # =========================================================================
    ## combine data
    
    sample  = ROOT.RooCategory ('sample','sample'  , 'S' , 'N' )
    
    ## combine datasets
    from ostap.fitting.simfit import combined_data 
    vars    = ROOT.RooArgSet ( mass1 , mass2 , mass3 )
    dataset = combined_data  ( sample , vars , { 'N' : dataset1 , 'S' : dataset2 } )
    
    ## combine PDFs
    model_sim  = Models.SimFit (
        sample , { 'S' : model_S  , 'N' : model_N } , name = 'X'
        )
    
    # =========================================================================
    rC , fC = model_sim.fitTo ( dataset , silent = True )
    rC , fC = model_sim.fitTo ( dataset , silent = True )
    
    title = 'Results of simultaneous fit'
    logger.info ( '%s\n%s' % ( title , rC.table ( title = title , prefix = '# ' ) ) )

    with use_canvas ( 'test_simfit3:simfit N'   ) : fN  = model_sim.draw ( 'N'   , dataset , nbins = 50 )
    with use_canvas ( 'test_simfit3:simfit S/x' ) : fSx = model_sim.draw ( 'S/x' , dataset , nbins = 50 )
    with use_canvas ( 'test_simfit3:simfit S/y' ) : fSy = model_sim.draw ( 'S/y' , dataset , nbins = 50 )
    graphs.append ( fN  )
    graphs.append ( fSx )
    graphs.append ( fSy )
    
    logger.info ( ' Value |        Simple fit         |    Combined fit ' )
    logger.info ( ' #N    | %25s | %-25s '  % ( rS.SSM2     * 1 , rC.SSM2     * 1 ) )
    logger.info ( ' mean  | %25s | %-25s '  % ( rS.mean_G2  * 1 , rC.mean_G2  * 1 ) )
    logger.info ( ' sigma | %25s | %-25s '  % ( rS.sigma_G2 * 1 , rC.sigma_G2 * 1 ) )
    
    with use_canvas ( 'test_simfit3:  profile-sim' , wait = 3 ) :
        
        grs = model_sim.graph_profile ( 'SSM2' , vrange ( 100 , 1000 , 50 ) , dataset , draw = True )
        grs.draw('apl')
        graphs.append ( grs )


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

    with timing ("simfit-3", logger ) :
        test_simfit3 () 


# =============================================================================
##                                                                      The END 
# =============================================================================
