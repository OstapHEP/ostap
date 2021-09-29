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
import ROOT, time 
import ostap.fitting.roofit 
import ostap.fitting.models  as     Models 
from   ostap.core.core       import Ostap, std, VE, dsID
from   ostap.logger.utils    import rooSilent 
import ostap.io.zipshelve    as     DBASE
from   ostap.utils.timing    import timing
from   ostap.plotting.canvas import use_canvas 
from   builtins              import range
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

GeV      = 1.0
MeV      = 0.001 * GeV

m_pi     = 139 *  MeV
m_rho    = 770 * MeV
g_rho    = 150 * MeV

m_etap   = 958 *  MeV 

m_phi    = 1019.46 * MeV
g_phi    =   4.249 * MeV
m_K      = 493.677 * MeV

m_Bs     = 5.366 * GeV
m_X      = 3.872 * GeV 

m0_f0    = 980   * MeV
m0g1_f0  = 0.165 * GeV**2 * 4.21 
g2og1_f0 = 1/4.21 

models   = set() 

# =============================================================================
## Different rho0 parameterizations 
def test_breitwigner_rho () :
    """Different rho0 parameterizations
    """
    
    logger  = getLogger ( "test_breitwigner_rho" )
    logger.info ( "Rho0 shapes" ) 
                  
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


    with use_canvas ( 'test_breitwigner_rho' ) : 
        bw2.draw (          xmin =  200 * MeV , xmax = 1.6 * GeV , linecolor = 4 )
        bw1.draw ( 'same' , xmin =  200 * MeV , xmax = 1.6 * GeV , linecolor = 2 )
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
    
    with use_canvas ( 'test_breitwigner_rho' ) : 
        f1 = model1.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 2 ) , ) )
        f2 = model2.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 4 ) , ) )
        f3 = model3.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 8 ) , ) )
        
        f2.draw ()
        f1.draw ( 'same' )
        f3.draw ( 'same' )
        
    models.add ( bw1    )
    models.add ( bw2    )
    models.add ( bw3    )
    models.add ( model1 )
    models.add ( model2 )
    models.add ( model3 )
    
    time.sleep ( 2 ) 

# =============================================================================
## Different phi0 parameterizations 
def test_breitwigner_phi () : 
    """Different phi0 parameterizations
    """
    
    logger  = getLogger ( "test_breitwigner_phi" )
    logger.info ( "Phi0 shapes" ) 

    ## 1) P-wave Breit-Wigner with Jackson's formfactor 
    bw1 = Ostap.Math.Phi0 ( m_phi , g_phi , m_K )

    ## 2) P-wave Breit-Wigner with Blatt-Weisskopf formfactor
    ff  = Ostap.Math.FormFactors.BlattWeisskopf ( 1 , 3.5 / GeV )
    ch2 = Ostap.Math.Channel ( g_phi , m_K , m_K , 1 , ff ) 
    bw2 = Ostap.Math.BreitWigner ( m_phi , ch2 )

    ## 3) P-wave Breit-Wigner with no formfactors 
    bw3 = Ostap.Math.BreitWigner ( m_phi , g_phi , m_K , m_K , 1 )

    with use_canvas ( 'test_breitwigner_phi' ) : 
        bw1.draw (          xmin = 0.95 * GeV , xmax = 1.5 * GeV , linecolor = 2 )
        bw2.draw ( 'same' , xmin = 0.95 * GeV , xmax = 1.5 * GeV , linecolor = 4 )
        bw3.draw ( 'same' , xmin = 0.95 * GeV , xmax = 1.5 * GeV , linecolor = 5 )
        
    logger.info ("bw1 fraction %.4f" % ( bw1.integral ( 1.1 * GeV , 1.5 * GeV ) / bw1.integral ( 0.9 * GeV , 1.1 * GeV ) ) )
    logger.info ("bw2 fraction %.4f" % ( bw2.integral ( 1.1 * GeV , 1.5 * GeV ) / bw2.integral ( 0.9 * GeV , 1.1 * GeV ) ) )
    logger.info ("bw3 fraction %.4f" % ( bw3.integral ( 1.1 * GeV , 1.5 * GeV ) / bw3.integral ( 0.9 * GeV , 1.1 * GeV ) ) )

    models.add ( bw1    )
    models.add ( bw2    )
    models.add ( bw3    )
    
    time.sleep ( 2 ) 

# =============================================================================
## Phi  shapes woth    phase space   corrections 
def test_breitwigner_phi_ps () : 
    """Phi  shapes woth    phase space   corrections 
    """
    
    logger  = getLogger ( "test_breitwigner_phi_ps" )
    logger.info ( "Phi0 shapes with phase space corrections" ) 

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


    logger.info (" f1 fraction %.4f" % (  f1.integral ( 1.1 * GeV , 1.5 * GeV ) /  f1.integral ( 0.9 * GeV , 1.1 * GeV ) ) )

    logger.info (" f2 fraction %.4f" % (  f2.integral ( 1.1 * GeV , 1.5 * GeV ) /  f2.integral ( 0.9 * GeV , 1.1 * GeV ) ) )
    
    logger.info (" f3 fraction %.4f" % (  f3.integral ( 1.1 * GeV , 1.5 * GeV ) /  f3.integral ( 0.9 * GeV , 1.1 * GeV ) ) )

    mass    = ROOT.RooRealVar  ('mass' , 'm(KK)' , 0.96 * GeV , 1.5 * GeV ) 
    
    p1 = Models.BWPS_pdf ( 'P1' , f1 , xvar = mass , m0 = m_phi  , gamma = g_phi    )
    p2 = Models.BWPS_pdf ( 'P2' , f2 , xvar = mass , m0 = p1.m0  , gamma = p1.gamma )
    p3 = Models.BWPS_pdf ( 'P3' , f3 , xvar = mass , m0 = p1.m0  , gamma = p1.gamma )
    
    with use_canvas ( 'test_breitwigner_phi_ps' ) : 
        fr1 = p1.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 2 ) , ) ) 
        fr2 = p2.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 4 ) , ) ) 
        fr3 = p3.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 8 ) , ) ) 
        
    with use_canvas ( 'test_breitwigner_phi_ps' ) : 
        fr1.draw ()
        fr2.draw ('same')
        fr3.draw ('same')
        
    ## flatte
    flatte = Ostap.Math.Flatte ( m0_f0 , m0g1_f0 , g2og1_f0 , m_K , m_K , m_pi , m_pi , 0.0 )
    flatte.draw ( xmin = 960 * MeV , xmax = 1.07 * GeV ) 
    
    ##  pdf 
    f0_980 = Models.Flatte_pdf ( 'F0' ,
                                 flatte = flatte   ,  
                                 xvar   = mass     ,
                                 m0     = ( 980 * MeV , 950 * MeV , 1000 * MeV ) ,
                                 m0g1   = m0g1_f0  ,
                                 g2og1  = g2og1_f0 ,
                                 gamma0 = 0        )
    f0_980.m0    .fix ( m0_f0 )
    f0_980.gamma0.fix ( 0     )
    
    flatte_ps = Ostap.Math.BWPS ( f0_980.pdf.flatte () , ps , True , True )
    f0_ps     = Models.FlattePS_pdf ( 'FP' ,
                                      flatte = flatte_ps      ,
                                      xvar   = mass           ,
                                      m0     = f0_980.m0      ,
                                      m0g1   = f0_980.m0g1    ,
                                      g2og1  = f0_980.g2og1   ,
                                      gamma0 = f0_980.gamma0  )
    
    with use_canvas ( 'test_breitwigner_phi_ps' ) : 
        fr4 = f0_980.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 6 ) , ) ) 
        fr5 = f0_ps .draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 7 ) , ) ) 

    with use_canvas ( 'test_breitwigner_phi_ps' ) : 
        fr5.draw (     )
        fr4.draw ('same')
        
    with use_canvas ( 'test_breitwigner_phi_ps' ) : 
        fr1.draw()
        fr2.draw('same')
        fr3.draw('same')
        fr4.draw('same')
        fr5.draw('same')
        
    models.add ( bw1    )
    models.add ( bw2    )
    models.add ( bw3    )

    models.add ( flatte )
    models.add ( flatte_ps )

    models.add ( f1     )
    models.add ( f2     )
    models.add ( f3     )
    models.add ( p1     )
    models.add ( p2     )
    models.add ( p3     )
    models.add ( f0_980 )
    models.add ( f0_ps  )

    time.sleep ( 2 ) 

# =============================================================================
## Rho0 shape from eta'  decays
def test_breitwigner_rho_more () : 
    """Rho0 shape from eta'  decays
    """
    
    logger  = getLogger ( "test_breitwigner_rho_more" )
    logger.info ( "Rho0 shape from eta'  decays" ) 

    ## Rho-profile with Gounaris-Sakurai lineshape
    ch4 = Ostap.Math.ChannelGS   ( g_rho , m_pi ) 
    bw4 = Ostap.Math.BreitWigner ( m_rho , ch4  )

    ## Rho-profile from eta' decays
    bw5 = Ostap.Math.BW3L        ( bw4 , m_etap , m_pi , m_pi , 0 , 1 )
    
    
    mass    = ROOT.RooRealVar  ('mass' , 'm(pipi)' , 200 * MeV , 1.6 * GeV ) 

    model4 = Models.BreitWigner_pdf ( 'BW4' , bw4 , xvar = mass ,  m0 = m_rho , gamma = g_rho ) 
    model5 = Models.BW3L_pdf        ( 'BW5' , bw5 , xvar = mass ,  m0 = m_rho , gamma = g_rho ) 

    with use_canvas ( 'test_breitwigner_rho_more' ) : 
        f4 = model4.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 4 ) , ) )
        f5 = model5.draw ( total_fit_options = ( ROOT.RooFit.LineColor ( 8 ) , ) )
        f5.draw ()
        f4.draw ( 'same' )

    models.add ( bw4    )
    models.add ( bw5    )
    models.add ( model4 )
    models.add ( model5 )

    time.sleep ( 2 ) 

# =============================================================================
## check that everything is serializable
def test_db() :
    """check that everything is serializable
    """

    
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db : 
        for i, m in enumerate ( models ) :
            db ['model/%2d: %s' % ( i , type ( m ).__name__  ) ] = m
        db['models'   ] = models
        db.ls() 

# =============================================================================
if '__main__' == __name__ :

 
    test_breitwigner_rho      ()        
    test_breitwigner_phi      ()       
    test_breitwigner_phi_ps   ()
    test_breitwigner_rho_more ()

    ## check finally that everything is serializeable:
    test_db ()

# =============================================================================
##                                                                      The END  
# =============================================================================

