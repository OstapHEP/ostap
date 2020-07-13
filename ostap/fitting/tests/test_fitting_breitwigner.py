#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_breitwigner.py
# test for some Breit-Wigner models 
# ============================================================================= 
""" Test module for soem BReit-Wigner models
"""
# ============================================================================= 
from   __future__        import print_function
# ============================================================================= 
import ROOT, random
import ostap.fitting.roofit 
import ostap.fitting.models as     Models 
from   ostap.core.core      import Ostap, std, VE, dsID
from   ostap.logger.utils   import rooSilent 
import ostap.io.zipshelve   as     DBASE
from   ostap.utils.timing   import timing 
from   builtins             import range
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_breitwigner' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make simple test mass

GeV   = 1.0
MeV   = 0.001 * GeV

m_pi  = 139 *  MeV
m_rho = 770 * MeV
g_rho = 150 * MeV



def test_breitwigner_rho () : 

    ## 1) P-wave Breit-Wigner with Jackson's formfactor 
    bw1 = Ostap.Math.Rho0 ( m_rho ,   g_rho , m_pi )
    
    ## 2) P-wave Breit-Wigner with Blatt-Weisskopf formfactor
    ff  = Ostap.Math.FormFactors.BlattWeisskopf ( 1 , 3.5 / GeV )
    ch2 = Ostap.Math.Channel ( g_rho , m_pi , m_pi , 1 , ff ) 
    bw2 = Ostap.Math.BreitWigner ( m_rho , ch2 )
    
    ## 3) Gounaris-Sakurai lineshape
    ch2 = Ostap.Math.ChannelGS   ( g_rho , m_pi ) 
    bw3 = Ostap.Math.BreitWigner ( m_rho , ch2  )

    bw1.draw (          xmin =  200 * MeV , xmax = 1.6 * GeV , linecolor = 2 )
    bw2.draw ( 'same' , xmin =  200 * MeV , xmax = 1.6 * GeV , linecolor = 4 )
    bw3.draw ( 'same' , xmin =  200 * MeV , xmax = 1.6 * GeV , linecolor = 8 )
    
    mass    = ROOT.RooRealVar  ('mass' , 'm(pipi)' , 200 * MeV , 1.6 * GeV ) 

    model1 = Models.BreitWigner_pdf ( 'BW1' , bw1 , xvar = mass , mean = m_rho , gamma = g_rho ) 
    model2 = Models.BreitWigner_pdf ( 'BW2' , bw2 ,
                                      xvar  = model1.mass  ,
                                      mean  = model1.mean  ,
                                      gamma = model1.gamma )
    model3 = Models.BreitWigner_pdf ( 'BW3' , bw3 ,
                                      xvar  = model1.mass  ,
                                      mean  = model1.mean  ,
                                      gamma = model1.gamma )     
    
    f1 = model1.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 2 ) , ) )
    f2 = model2.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 4 ) , ) )
    f3 = model3.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 8 ) , ) )
    
    f1.draw ()
    f2.draw ( 'same' )
    f3.draw ( 'same' )
    
    

# =============================================================================
if '__main__' == __name__ :

    test_breitwigner_rho ()
    
    
# =============================================================================
##                                                                      The END  
# =============================================================================

