#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_splot.py
# Test module 
# - It tests sPlotting 
# ============================================================================= 
""" Test module 
- It tests sPlotting
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import cpp, VE, dsID, rooSilent 
from   builtins                 import range
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_splot' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

def test_splot () :

    logger = getLogger ( 'test_splot' )
    
    ## make simple test mass 
    mass    = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 0 , 10 )
    tau     = ROOT.RooRealVar ( 'test_tau'  , 'Some test tau'  , 0 , 10 )
    
    ## book very simple data set
    varset  = ROOT.RooArgSet  ( mass , tau )
    dataset = ROOT.RooDataSet ( dsID() , 'Test Data set' , varset )  
    
    mmin, mmax = mass.minmax()
    tmin, tmax = tau.minmax ()
    
    m0   = VE(3,0.2**2)
    taus = 6
    taub = 0.4
    NS   = 5000
    NB   = 5000
    
    ## generate signal 
    for i in range(0,NS) :
        m =  m0.gauss()
        while not mmin < m < mmax : m =  m0.gauss() 
        t = random.expovariate ( 1./ taus )
        while not tmin < t < tmax : t = random.expovariate ( 1./ taus )
        mass.value = m
        tau.value  = t 
        dataset.add ( varset )
        
    ## generate background
    for i in range(0,NB) :
        m = random.expovariate ( 1./ 5  )
        while not mmin < m < mmax : m = random.expovariate ( 1./ 3 )
        t = random.expovariate ( 1./ taub )
        while not tmin < t < tmax : t = random.expovariate ( 1./ taus )
        mass.value = m
        tau.value  = t 
        dataset.add ( varset )

    logger.info ( "Original dataset\n%s" % dataset.table ( prefix = '# ' ) )

    signal = Models.Gauss_pdf ( 'G'  , xvar = mass  , mean = ( 3 , 2, 4 ) , sigma = ( 0.2 , 0.1 , 0.5 ) )

    model = Models.Fit1D ( signal = signal , background = 'expo-' )
    model.S = NS
    model.B = NB

    model.background.tau.fix ( -0.2 )
    signal.mean .fix ( m0.value () )
    signal.sigma.fix ( m0.error () )
    model.fitTo ( dataset , silent = True )
    signal.mean .release()
    model.fitTo ( dataset , silent = True )
    signal.sigma.release() 
    model.fitTo ( dataset , silent = True ) 
    model.background.tau.release  () 
    model.fitTo ( dataset , silent = True )
    
    with use_canvas ( 'test_splot' ) : 
        r , f = model.fitTo ( dataset , silent = True , draw = True , nbins = 50 )
    logger.info ( "Mass fit : fit results\n%s" % r.table ( title = 'Mass fit' , prefix = '# ' ) )

    ## make splot analysis
    model.sPlot ( dataset )
    logger.info ( "Dataset after sPlot\n%s" % dataset.table ( prefix = '# ' ) )

    ## make signal-weighted dataset
    ds_signal = dataset.makeWeighted  ( 'S_sw' )

    ## make background-weighted dataset
    ds_bkg    = dataset.makeWeighted  ( 'B_sw' )
    
    logger.info ( "Signal-weighted dataset\n%s" % ds_signal.table ( prefix = '# ' ) )
    logger.info ( "Background-weighted dataset\n%s" % ds_bkg.table ( prefix = '# ' ) )
    
    ##  make exponential fits for signal and background samples 
    TS = Models.Bkg_pdf  ( 'TS' , xvar = tau , power = 0 )
    TB = Models.Bkg_pdf  ( 'TB' , xvar = tau , power = 0 )
    
    TS.fitTo ( ds_signal , silent = True , sumw2 = True )
    TS.fitTo ( ds_signal , silent = True , sumw2 = True )
    with use_canvas ( 'test_splot' ) : 
        rS, f = TS.fitTo ( ds_signal , silent = True , draw  = True , nbins = 100 , sumw2 = True )
    logger.info ( "Tau/signal fit : fit results\n%s" % rS.table ( title = 'Tau signal fit' , prefix = '# ' ) )


    TB.fitTo ( ds_bkg , silent = True , sumw2 = True )
    TB.fitTo ( ds_bkg , silent = True , sumw2 = True )
    with use_canvas ( 'test_splot' ) : 
        rB, f = TB.fitTo ( ds_bkg , silent = True , draw  = True , nbins = 100 , sumw2 = True )
    logger.info ( "Tau/bkg fit : fit results\n%s" % rB.table ( title = 'Tau bkg fit' , prefix = '# ' ) )

    
    logger.info ( "Tau/signal : %28s vs %s" % ( abs ( 1.0 / rS.tau_TS ) , taus ) ) 
    logger.info ( "Tau/bkg    : %28s vs %s" % ( abs ( 1.0 / rB.tau_TB ) , taub ) ) 
        
# =============================================================================

if '__main__' == __name__ :

    with timing ( "sPlot" , logger ) :  
        test_splot        () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
