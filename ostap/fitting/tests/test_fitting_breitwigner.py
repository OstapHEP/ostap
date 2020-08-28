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

m_phi = 1019.46 * MeV
g_phi =   4.249 * MeV
m_K   = 493.677 * MeV

m_Bs  = 5.366 * GeV
m_X   = 3.872 * GeV 

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

    ## 3) P-wave Breit-Wigner with no formfactors 
    bw4 = Ostap.Math.BreitWigner ( m_rho , g_rho , m_pi , m_pi , 1 )

    bw1.draw (          xmin =  200 * MeV , xmax = 1.6 * GeV , linecolor = 2 )
    bw2.draw ( 'same' , xmin =  200 * MeV , xmax = 1.6 * GeV , linecolor = 4 )
    bw3.draw ( 'same' , xmin =  200 * MeV , xmax = 1.6 * GeV , linecolor = 8 )
    bw4.draw ( 'same' , xmin =  200 * MeV , xmax = 1.6 * GeV , linecolor = 5 )
    
    mass    = ROOT.RooRealVar  ('mass' , 'm(pipi)' , 200 * MeV , 1.6 * GeV ) 

    model1 = Models.BreitWigner_pdf ( 'BW1' , bw1 , xvar = mass ,  m0 = m_rho , gamma = g_rho ) 
    model2 = Models.BreitWigner_pdf ( 'BW2' , bw2 ,
                                      xvar  = model1.mass  ,
                                      m0    = model1.mean  ,
                                      gamma = model1.gamma )
    model3 = Models.BreitWigner_pdf ( 'BW3' , bw3 ,
                                      xvar  = model1.mass  ,
                                      m0    = model1.mean  ,
                                      gamma = model1.gamma )     
    
    f1 = model1.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 2 ) , ) )
    f2 = model2.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 4 ) , ) )
    f3 = model3.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 8 ) , ) )
    
    f1.draw ()
    f2.draw ( 'same' )
    f3.draw ( 'same' )
    

# =============================================================================
def test_breitwigner_phi () : 

    ## 1) P-wave Breit-Wigner with Jackson's formfactor 
    bw1 = Ostap.Math.Phi0 ( m_phi , g_phi , m_K )

    ## 2) P-wave Breit-Wigner with Blatt-Weisskopf formfactor
    ff  = Ostap.Math.FormFactors.BlattWeisskopf ( 1 , 3.5 / GeV )
    ch2 = Ostap.Math.Channel ( g_phi , m_K , m_K , 1 , ff ) 
    bw2 = Ostap.Math.BreitWigner ( m_phi , ch2 )

    ## 3) P-wave Breit-Wigner with no formfactors 
    bw3 = Ostap.Math.BreitWigner ( m_phi , g_phi , m_K , m_K , 1 )

    bw1.draw (          xmin = 0.95 * GeV , xmax = 1.5 * GeV , linecolor = 2 )
    bw2.draw ( 'same' , xmin = 0.95 * GeV , xmax = 1.5 * GeV , linecolor = 4 )
    bw3.draw ( 'same' , xmin = 0.95 * GeV , xmax = 1.5 * GeV , linecolor = 5 )

    logger.info ("bw1 fraction %s" % ( bw1.integral ( 1.1 , 1.5 ) / bw1.integral ( 0.9 , 1.1 ) ) )
    logger.info ("bw2 fraction %s" % ( bw2.integral ( 1.1 , 1.5 ) / bw2.integral ( 0.9 , 1.1 ) ) )
    logger.info ("bw3 fraction %s" % ( bw3.integral ( 1.1 , 1.5 ) / bw3.integral ( 0.9 , 1.1 ) ) )

# =============================================================================
def test_breitwigner_phi_ps () : 
##  if 1 < 2 :
    
    ## 1) P-wave Breit-Wigner with Jackson's formfactor 
    bw1 = Ostap.Math.Phi0 ( m_phi , g_phi , m_K )

    ## 2) P-wave Breit-Wigner with Blatt-Weisskopf formfactor
    ff  = Ostap.Math.FormFactors.BlattWeisskopf ( 1 , 3.5 / GeV )
    ch2 = Ostap.Math.Channel ( g_phi , m_K , m_K , 1 , ff ) 
    bw2 = Ostap.Math.BreitWigner ( m_phi , ch2 )

    ## 3) P-wave Breit-Wigner with no formfactors 
    bw3 = Ostap.Math.BreitWigner ( m_phi , g_phi , m_K , m_K , 1 )

    ps = Ostap.Math.PhaseSpaceNL ( 2 * m_K , m_Bs - m_X , 0 , 3 - 2 ) 

    f1 = Ostap.Math.BWPS ( bw1 , ps , True , True )
    f2 = Ostap.Math.BWPS ( bw2 , ps , True , True )
    f3 = Ostap.Math.BWPS ( bw3 , ps , True , True )
    
    f1.draw (          linecolor = 2 )
    f2.draw ( 'same' , linecolor = 4 )
    f3.draw ( 'same' , linecolor = 5 )

    logger.info (" f1 fraction %s" % ( f1.integral ( 1.1 , 1.5 ) / f1.integral ( 0.9 , 1.1 ) ) )
    logger.info (" f2 fraction %s" % ( f2.integral ( 1.1 , 1.5 ) / f2.integral ( 0.9 , 1.1 ) ) )
    logger.info (" f3 fraction %s" % ( f3.integral ( 1.1 , 1.5 ) / f3.integral ( 0.9 , 1.1 ) ) )


# =============================================================================
if '__main__' == __name__ :

    ## pass 
    test_breitwigner_rho    ()
    test_breitwigner_phi    ()
    test_breitwigner_phi_ps ()
    
    
# =============================================================================
##                                                                      The END  
# =============================================================================

