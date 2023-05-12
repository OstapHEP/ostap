#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_simfit6.py
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
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models
from   ostap.core.meta_info     import root_info 
from   ostap.core.core          import hID , dsID, rooSilent 
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait
from   ostap.fitting.simfit     import combined_data, combined_hdata
from   ostap.fitting.fithelpers import H1D_dset
import ostap.io.zipshelve       as     DBASE
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_simfit6' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================


# =============================================================================
def test_simfit6() :
    
    logger = getLogger ( 'test_simfit6' ) 
    if  (6,18) <= root_info < (6,20) :
        logger.warning ( "Test is disabled for ROOT verison %s" % str ( root_info ) )
        return 

    ## make simple test mass 
    mass     = ROOT.RooRealVar ( 'mass' , 'Some test mass (for splot)' ,  0 , 10 )
    m1       = ROOT.RooRealVar ( 'm1'   , 'another mass'               ,  0 ,  5 )
    ## and three mass-differences 
    dm1      = ROOT.RooRealVar ( 'dm1'  , 'delta-m1' ,  -1 , 1 )
    
    # =============================================================================
    ## prepare dataset for splot
    # =============================================================================
    signal     = Models.Gauss_pdf ( 'SG' , xvar = mass , mean  = 5 , sigma =  0.5 )
    bkg        = Models.Bkg_pdf   ( 'BG' , xvar = mass , power = 1 , tau   = -0.2 )
    
    ## "signal"
    signal_cmp = Models.Fit1D  (
        signal     = Models.Apollonios2_pdf ( 'A'             ,
                                              xvar      = m1  ,
                                              mean      = 2.5 ,
                                              sigma     = 0.2 ,
                                              beta      = 0.7 ,
                                              asymmetry = 0.1 ) , 
        background = 0   , suffix = '_Scmp' )
    signal_cmp.S = 1000
    signal_cmp.B =  500 
    
    # =============================================================================
    ## model for fitting MC 
    model_mc = signal_cmp.signal.clone ( name = 'SMC' )
    
    ## mc_dataset
    mc_dset = model_mc.generate ( 100000 )
    hmc     = m1.histo ( 500 )
    mc_dset.project ( hmc , model_mc.xvar.name )
    
    model_mc.fitTo ( hmc , silent = True )
    rmch , fmch = model_mc.fitTo ( hmc , silent = True , draw = True , nbins = 100 )
    title   = 'Fit for MC data (histo)'
    logger.info ( '%s:\n%s' % ( title , rmch.table ( title = title , prefix = '# ' ) ) )
    
    ## create the constrains according to MC
    mc_constraint = model_mc.soft_multivar_constraint ( ( model_mc.mean  ,
                                                          model_mc.sigma ,
                                                          model_mc.beta  ,
                                                          model_mc.kappa ) , rmch )
    
    # =============================================================================
    
    ## "signal"
    bkg_cmp = Models.Bkg_pdf ( 'BC' , xvar = m1 , power = 1 , tau = -1 )
    
    cmps = signal * signal_cmp
    cmpb = bkg    * bkg_cmp
    
    dss = cmps.generate ( 8000 )
    dsb = cmpb.generate ( 3000 )
    
    ds = dss + dsb 
    
    # ===========================================================================
    ## model for sPlot
    # ===========================================================================
    smodel = Models.Fit1D ( signal = signal , background = bkg , suffix = '_spl' )
    smodel.fitTo ( ds , silent = True )
    rs , fs = smodel.fitTo ( ds , silent = True , draw = True , nbins = 100 )
    title = 'Fit for sPlot'
    logger.info ( '%s:|n%s' % ( title , rs.table ( title = title , prefix = '# ' ) ) )
    smodel.sPlot ( ds )
    
    # ============================================================================
    ## background subtracted (weighted) datatset
    # ============================================================================
    weight = 'S_spl_sw'
    dsw = ds.makeWeighted ( weight )
    
    signal_cmp.fitTo ( dsw , silent = True , sumw2 = True )
    r1 , f1 = signal_cmp.fitTo ( dsw , silent = True , sumw2 = True , draw = True , nbins = 100 )
    title   = 'Fit for background subtracted dataset'
    logger.info ( '%s:\n%s' % ( title , r1.table ( title = title , prefix = '# ' ) ) )
    
    ## make a first histogram  (weighted) 
    hdata1 = m1.histo ( 200 )
    dsw.project ( hdata1 , 'm1' )
    
    ## fit data histo
    rhist , fhist = signal_cmp.fitTo ( hdata1 , silent = True ) 
    rhist , fhist = signal_cmp.fitTo ( hdata1 , silent = True ) 
    rhist , fhist = signal_cmp.fitTo ( hdata1 , silent = True , draw = True , nbins = 100 ) 
    title   = 'Fit histogram'
    logger.info ( '%s:\n%s' % ( title , rhist.table ( title = title , prefix = '# ' ) ) )
    
    ## fit data (unbinned) 
    rdsw , fdsw = signal_cmp.fitTo ( dsw , sumw2 = True , silent = True ) 
    rdsw , fdsw = signal_cmp.fitTo ( dsw , sumw2 = True , silent = True ) 
    rdsw , fdsw = signal_cmp.fitTo ( dsw , sumw2 = True , silent = True , draw = True , nbins = 100 ) 
    title   = 'Fit unbinned dataset'
    logger.info ( '%s:\n%s' % ( title , rdsw.table ( title = title , prefix = '# ' ) ) )
    
    ## fit data histo with constrains 
    rhistc , fhistc = signal_cmp.fitTo ( hdata1 , silent = True , constraints = ( mc_constraint, ) ) 
    rhistc , fhistc = signal_cmp.fitTo ( hdata1 , silent = True , constraints = ( mc_constraint, ) ) 
    rhistc , fhistc = signal_cmp.fitTo ( hdata1 , silent = True , constraints = ( mc_constraint, ) , draw = True , nbins = 100 ) 
    title   = 'Fit histogram with constraint'
    logger.info ( '%s:\n%s' % ( title , rhistc.table ( title = title , prefix = '# ' ) ) )
    
    ## fit data (unbinned) with constrains 
    rdswc , fdswc = signal_cmp.fitTo ( dsw , sumw2 = True , silent = True , constraints = ( mc_constraint, ) ) 
    rdswc , fdswc = signal_cmp.fitTo ( dsw , sumw2 = True , silent = True , constraints = ( mc_constraint, ) ) 
    rdswc , fdswc = signal_cmp.fitTo ( dsw , sumw2 = True , silent = True , constraints = ( mc_constraint, ) , draw = True , nbins = 100 ) 
    title   = 'Fit unbinned with constraint'
    logger.info ( '%s:\n%s' % ( title , rdswc.table ( title = title , prefix = '# ' ) ) )
    
    ## tru to serialize everything
    logger.info('Saving all objects into DBASE')
    with timing ('Save everything to DBASE' , logger ), DBASE.tmpdb() as db : 
        db['mass'          ] = mass
        db['m1'            ] = m1
        db['signal'        ] = signal
        db['bkg'           ] = bkg 
        db['signal_cmp'    ] = signal_cmp 
        db['model_mc'      ] = model_mc 
        db['hmc'           ] = hmc 
        db['hdata1'        ] = hdata1

        db['rmch'          ] = rmch 
        db['fmch'          ] = fmch 

        db['mc_constraint' ] = mc_constraint

        db['bkg_cmp'       ] = bkg_cmp 
        db['dss'           ] = dss
        db['dsb'           ] = dsb
        db['ds'            ] = ds
        db['dsw'           ] = dsw
        db['mc_dset'       ] = mc_dset
        
        db['smodel'        ] = smodel

        db['rs'            ] = rs
        db['fs'            ] = fs

        db['r1'            ] = r1
        db['f1'            ] = f1

        db['rhist'         ] = rhist
        db['fhist'         ] = fhist

        db['rdsw'          ] = rdsw
        db['fdsw'          ] = fdsw

        db['rhistc'        ] = rhistc
        db['fhistc'        ] = fhistc

        db['rdswc'         ] = rdswc
        db['fdswc'         ] = fdswc
        
        db.ls()

# =============================================================================
if '__main__' == __name__ :

    with timing ("simfit-6", logger ) :
        test_simfit6 () 



# =============================================================================
##                                                                      The END 
# =============================================================================
