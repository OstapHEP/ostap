#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_simfit5.py
# Test module for ostap/fitting/simfit.py
# - It tests the most simple "Simultaneous fit"
# ============================================================================= 
""" Test module for ostap/fitting/simfit.py
- It tests the most simple ``Simultaneous fit''
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.core          import dsID, rooSilent 
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.timing       import wait 
from   ostap.utils.root_utils   import batch_env 
from   ostap.fitting.simfit     import combined_data
from   ostap.fitting.fithelpers import H1D_dset
import ostap.fitting.models     as     Models 
import ostap.io.zipshelve       as     DBASE
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_simfit5' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
#
## make simple test mass 
mass     = ROOT.RooRealVar ( 'mass' , 'Some test mass' ,  0 , 10 )
## and three mass-differences 
dm1      = ROOT.RooRealVar ( 'dm1'  , 'delta-m1' ,  -10 , 10 )
dm2      = ROOT.RooRealVar ( 'dm2'  , 'delta-m2' ,  -10 , 10 )
dm3      = ROOT.RooRealVar ( 'dm3'  , 'delta-m3' ,  -10 , 10 )

## historgams 
hdm1     = dm1.histo ( 10000 )
hdm2     = dm2.histo ( 10000 )
hdm3     = dm3.histo ( 10000 )

# =============================================================================
## book very simple data set for ``DATA''
# =============================================================================
varset   = ROOT.RooArgSet  ( mass )
dataset  = ROOT.RooDataSet ( dsID() , 'Test Data set' , varset )  

mean1    = 3
sigma1   = 0.5
NS1      = 250
for i in range ( NS1 ) :
    v1 = random.gauss ( mean1 , sigma1 ) 
    if v1 in mass :
        mass.setVal ( v1     )
        dataset.add ( varset )
        
mean2    = 5
sigma2   = 0.5
NS2      = 250
for i in range ( NS2 ) :
    v2 = random.gauss ( mean2 , sigma2 ) 
    if v2 in mass :
        mass.setVal ( v2     )
        dataset.add ( varset )

mean3    = 7
sigma3   = 0.5
NS3      = 250
for i in range ( NS3 ) :
    v3 = random.gauss ( mean3 , sigma3 ) 
    if v3 in mass :
        mass.setVal ( v3     )
        dataset.add ( varset )

        
NB = 2000 
for i in range ( NB ) :
    v3 = random.uniform ( 0 , 10  )
    if v3 in mass :
        mass.setVal ( v3     )
        dataset.add ( varset )
        
vset1 = ROOT.RooArgSet  ( dm1 )
dset1 = ROOT.RooDataSet ( dsID () , 'Test data set 1: resolutuon for dm1' , vset1 )
for i in range ( 50000 ) :
    v1 = random.gauss ( 0 , sigma1 ) 
    if v1 in dm1 :
        hdm1.Fill  ( v1    )  
        dm1.setVal ( v1    )
        dset1.add  ( vset1 )
            
vset2 = ROOT.RooArgSet  ( dm2 )
dset2 = ROOT.RooDataSet ( dsID () , 'Test data set 2: resolutuon for dm2' , vset2 )
for i in range ( 50000 ) :
    v2 = random.gauss ( 0 , sigma2 ) 
    if v2 in dm2 :
        hdm2.Fill  ( v2    )  
        dm2.setVal ( v2    )
        dset2.add  ( vset2 )

vset3 = ROOT.RooArgSet  ( dm3 )
dset3 = ROOT.RooDataSet ( dsID () , 'Test data set 3: resolutuon for dm3' , vset3 )
for i in range ( 50000 ) :
    v3 = random.gauss ( 0 , sigma3 ) 
    if v3 in dm3 :
        hdm3.Fill  ( v3    )  
        dm3.setVal ( v3    )
        dset3.add  ( vset3 )

category = ROOT.RooCategory ( 'sample' , 'sample' , 'data' , 'dm1' , 'dm2' , 'dm3' )
vars     = ROOT.RooArgSet   ( mass, dm1 , dm2 , dm3 )

cdataset1 = combined_data ( category ,
                            vars     ,  
                            datasets = { 'data' : dataset ,
                                         'dm1'  : dset1   ,
                                         'dm2'  : dset2   ,
                                         'dm3'  : dset3   } )
logger.info ( 'Combined dataset/1:\n%s' % cdataset1.table ( prefix = '# ' ) )                            


# ======================================================================
## make try to use binned dataset  (as weighted datasets)
# ======================================================================

dsdm1  = H1D_dset ( hdm1 , xaxis = dm1 , weighted = True , skip_zero = True )
dsdm2  = H1D_dset ( hdm2 , xaxis = dm2 , weighted = True , skip_zero = True )
dsdm3  = H1D_dset ( hdm3 , xaxis = dm3 , weighted = True , skip_zero = True )

dataset.addVar( 'h1weight' , '1.0' ) 
cdataset2 = combined_data ( category ,
                            vars     ,  
                            datasets = { 'data' : dataset, ## dataset.makeWeighted ( 'h1weight' ) ,
                                         'dm1'  : dsdm1.dset ,
                                         'dm2'  : dsdm2.dset ,
                                         'dm3'  : dsdm3.dset } )
                            ## args = ( ROOT.RooFit.WeightVar( dsdm1.wname ) , ) )
                            

logger.info ( 'Combined dataset/2:\n%s' % cdataset2.table ( prefix = '# ' ) )                            


# =============================================================================
def test_simfit5() : 
    
    logger = getLogger ( 'test_simfit5' )

    # =========================================================================
    ## resolution for the left peak 
    # =========================================================================
    reso1 = Models.ResoGauss ( 'RG1'       ,
                               xvar  = dm1 ,
                               sigma = ( sigma1 , sigma1 / 2 , sigma1 * 2 ) ,
                               mean  = ( 0 , -1 , 1 ) )
    
    # fit DM1 dataset for resolution 
    r1  , f1  = reso1.fitTo ( dset1 , silent = True )
    r1  , f1  = reso1.fitTo ( dset1 , silent = True , draw = True , nbins = 100 )
    
    # =========================================================================
    ## resolution for the central peak 
    # =========================================================================
    reso2 = Models.ResoGauss ( 'RG2'       ,
                               xvar  = dm2 ,
                               sigma = ( sigma2 , sigma2 / 2 , sigma2 * 2 ) ,
                               mean  = ( 0 , -1 , 1 ) )

    # fit DM2 dataset for resolution 
    r2 , f2 = reso2.fitTo ( dset2 , silent = True )
    r2 , f2 = reso2.fitTo ( dset2 , silent = True , draw = True , nbins = 100 )

    # =========================================================================
    ## resolution for the right peak 
    # =========================================================================
    reso3 = Models.ResoGauss ( 'RG3'       ,
                               xvar  = dm3 ,
                               sigma = ( sigma3 , sigma3 / 2 , sigma3 * 2 ) ,
                               mean  = ( 0 , -1 , 1 ) )

    # fit DM2 dataset for resolution 
    r3 , f3 = reso3.fitTo ( dset3 , silent = True )
    r3 , f3 = reso3.fitTo ( dset3 , silent = True , draw = True , nbins = 100 )

    # =========================================================================
    # model for data fit 
    # =========================================================================
    signal1 = Models.Gauss_pdf ( 'G1' ,
                                 xvar  = mass        ,
                                 sigma = reso1.sigma , 
                                 mean  = ( mean1 , mean1 - 1 , mean1 + 1 ) )
    signal2 = Models.Gauss_pdf ( 'G2' ,
                                 xvar  = mass        ,
                                 sigma = reso2.sigma , 
                                 mean  = ( mean2 , mean2 - 1 , mean2 + 1 ) )
    signal3 = Models.Gauss_pdf ( 'G3' ,
                                 xvar  = mass        ,
                                 sigma = reso3.sigma , 
                                 mean  = ( mean3 , mean3 - 1 , mean3 + 1 ) )
    
    background = Models.PolyPos_pdf ( 'Bkg' , xvar = mass , power = 1 )
    
    model = Models.Fit1D ( signals    = [ signal1 , signal2 , signal3 ] ,
                           background = background            )
    model.S = NS1 , NS2 , NS3 
    model.B = NB

    # =========================================================================
    ## Fit data (without resolution samples)
    # =========================================================================
    r0 , f0 = model.fitTo ( dataset , silent = True )
    r0 , f0 = model.fitTo ( dataset , silent = True , draw = True , nbins = 100 )

    logger.info ( 'Fit result (only data):\n%s' % r0.table ( prefix = "# " ) )
    
    ## combine PDFs
    model_sim  = Models.SimFit (
        category , { 'data' : model ,
                     'dm1'  : reso1 ,
                     'dm2'  : reso2 ,
                     'dm3'  : reso3 } , name = 'X' )
    
    # =========================================================================
    ## Simultanegous fit data with resolution samples)
    # =========================================================================
    r_1 , f_1 = model_sim.fitTo ( cdataset1 , silent = True )
    r_1 , f_1 = model_sim.fitTo ( cdataset1 , silent = True )

    logger.info ( 'Simultaneous fit result (unbinned reso):\n%s' % r_1.table ( prefix = "# " ) )
    
    with use_canvas ( 'test_simfit5: unbinned resolutions' ) , wait ( 1 ) :
        with wait ( 1 ) : fdm1 = model_sim.draw ( 'dm1'  , cdataset1 , nbins = 100 )
        with wait ( 1 ) : fdm2 = model_sim.draw ( 'dm2'  , cdataset1 , nbins = 100 )
        with wait ( 1 ) : fdm3 = model_sim.draw ( 'dm3'  , cdataset1 , nbins = 100 )
        with wait ( 1 ) : fd   = model_sim.draw ( 'data' , cdataset1 , nbins = 100 )
        

    # =========================================================================
    ## Simultanegous fit data with resolution samples)
    # =========================================================================
    r_2 , f_2 = model_sim.fitTo ( cdataset2 , silent = True , sumw2 = True )
    r_2 , f_2 = model_sim.fitTo ( cdataset2 , silent = True , sumw2 = True )

    logger.info ( 'Simultaneous fit result (binned reso):\n%s' % r_2.table ( prefix = "# " ) )
    
    with use_canvas ( 'test_simfit5: binned resolutions' ) , wait ( 1 ) :
        with wait ( 1 ) : fdm1 = model_sim.draw ( 'dm1'  , cdataset2 , nbins = 100 )
        with wait ( 1 ) : fdm2 = model_sim.draw ( 'dm2'  , cdataset2 , nbins = 100 )
        with wait ( 1 ) : fdm3 = model_sim.draw ( 'dm3'  , cdataset2 , nbins = 100 )
        with wait ( 1 ) : fd   = model_sim.draw ( 'data' , cdataset2 , nbins = 100 )

    # =========================================================================
    ## test creation of dataset
    # =========================================================================
    ds_gen = model_sim.generate ( nEvents = { 'data' : len ( dataset ) ,
                                              'dm1'  : len ( dset1   ) ,
                                              'dm2'  : len ( dset2   ) ,
                                              'dm3'  : len ( dset3   ) } , 
                                  varset  = vars  )

    rg , f = model_sim.fitTo ( ds_gen , silent = True )
    rg , f = model_sim.fitTo ( ds_gen , silent = True )
    
    title = 'Results of simultaneous fit to generated dataset'
    logger.info ( '%s\n%s' % ( title , rg.table ( title = title , prefix = '# ' ) ) )

    
    ## try to serialize everything
    logger.info('Saving all objects into DBASE')
    with timing ('Save everything to DBASE' , logger ), DBASE.tmpdb() as db : 

        db ['reso1'  ] = reso1
        db ['signal1'] = signal1
        db ['signal2'] = signal2 
        db ['signal3'] = signal3 
        db ['model'  ] = model 
        
        db ['r1']    = r1
        db ['f1']    = f1

        db ['r2']    = r2
        db ['f2']    = f2

        db ['r3']    = r3
        db ['f3']    = f3
        
        db ['r0']    = r0
        db ['f0']    = f0
        
        db ['r_1']   = r_1
        db ['f_1']   = f_1
        
        db ['r_2']   = r_2
        db ['f_2']   = f_2
        
        db ['fdm1']  = fdm1
        db ['fdm2']  = fdm2
        db ['fdm3']  = fdm3
        db ['fd']    = fd

# =============================================================================
if '__main__' == __name__ :

    with timing ("simfit-5", logger ) :
        test_simfit5 () 

# =============================================================================
##                                                                      The END 
# =============================================================================
