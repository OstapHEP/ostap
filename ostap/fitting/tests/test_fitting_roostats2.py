#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_roostats.py
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
from   ostap.core.meta_info     import root_info 
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
    logger = getLogger ( 'test_fitting_roostats2' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

# =============================================================================
## more complicated/sophisticated  stuff
# =============================================================================

mass    = ROOT.RooRealVar   ('mass','mass-variable', 0 , 10 )
signal  = Models.Gauss_pdf  ( 'Gauss',
                              xvar  = mass                 ,
                              mean  = ( 2.5 , 0.5  , 9.5 ) ,
                              sigma = ( 0.3 , 0.01 , 3.0 ) )
signal.mean .fix()
signal.sigma.fix()
model   = Models.Fit1D ( signal = signal , background = 'e-' )
model.background.tau = -0.25 
model.S = 55
model.B = 1000
model.S.setMax(200)

data    = model.generate ( 55 + 1000 )
data1   = data.clone()
data2   = data.clone()
data3   = data.clone()
data4   = data.clone()
data5   = data.clone()
data6   = data.clone()

# ============================================================================-
## Get the upper limit limit for small signal at fixed mass
#  - resolution is fixed
#  - Asymptotic Calculator is used 
def test_point_limit_ac () :
    """Get the upper limit at given point for small signal at fixed mass
    - resoltuion is fixed 
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_point_limit_ac")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calcultor)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )

    the_model = model.clone ( name = 'M1' )

    with use_canvas ( 'test_point_limit_ac' ) : 
        logger.info ( 'Dataset is\n%s' % data1.table ( prefix = '# ' ) ) 
        rr , frame = the_model.fitTo ( data1 , draw = True , nbins = 50 )

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf       = the_model   ,
                             poi       = the_model.S , ## parameter of interest 
                             dataset   = data1       ,
                             name      = 'S+B'       )
    
    model_sb.snapshot = the_model.S ## ATTENTION! 
    
    ## create ModelConfig  for 'B-only' model
    model_b  = ModelConfig ( pdf       = the_model          ,
                             poi       = the_model.S        , ## parameter of interest 
                             dataset   = data1              ,
                             workspace = model_sb.workspace , 
                             name      = 'B-only'           )
    
    the_model.S = 0 
    model_b.snapshot = the_model.S  ## ATTENTION! 
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )
    
    model_sb.mc.Print('vvv')
    model_b .mc.Print('vvv')


    with timing ( "Using Asymptotic Calculator" , logger = logger ) :
        ## create the calculator 
        ac  = AsymptoticCalculator ( model_b           ,
                                     model_sb          ,
                                     dataset   = data  ,
                                     asimov    = False , 
                                     one_sided = True  )
        
        ## create Hypo Test inverter 
        hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
        
        ## make a scan 
        hti .scan ( vrange ( 0 , 150 , 50 )  ) ## scan it!
        
    ## visualize the scan results 
    with use_canvas ( 'test_pointLimit: HypoTestInverter plot (asymptotic)' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.1f' % hti.upper_limit )


# ============================================================================-
## Get the upper limit limit for small signal at fixed mass
#  - resolution is fixed
#  - Frequestist Calculator is used 
def test_point_limit_fc  () :
    """Get the upper limit at given point for small signal at fixed mass
    - resolution is fixed 
    - Frequestist Calculator is used 
    """

    logger = getLogger("test_point_limit_fc")

    logger.info ( "Test Point limits with RooStats using Frequestist Calculator" )

    ## if root_info < (6,24) :
    ##    logger.info ( 'Test is disabled for ROOT version %s' % str ( root_info ) )
    ##    return 
    
    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             FrequentistCalculator ,
                                             HypoTestInverter      )

    the_model = model.clone ( name = 'M2' )

    with use_canvas ( 'test_point_limit_fc' ) : 
        logger.info ( 'Dataset is\n%s' % data2.table ( prefix = '# ' ) ) 
        rr , frame = the_model.fitTo ( data2 , draw = True , nbins = 50 )

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf       = the_model   ,
                             poi       = the_model.S , ## parameter of interest 
                             dataset   = data2       ,
                             name      = 'S+B'       )
    
    model_sb.snapshot = the_model.S ## ATTENTION! 
    
    ## create ModelConfig  for 'B-only' model
    model_b  = ModelConfig ( pdf       = the_model          ,
                             poi       = the_model.S        , ## parameter of interest 
                             dataset   = data2              ,
                             workspace = model_sb.workspace , 
                             name      = 'B-only'           )
    
    the_model.S = 0 
    model_b.snapshot = the_model.S  ## ATTENTION! 
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )
    
    model_sb.mc.Print('vvv')
    model_b .mc.Print('vvv')

    ## with Frequentist calculator
    with timing ( "Using Frequentist Calculator" , logger = logger ) :
        
        ## create the calculator 
        fc  = FrequentistCalculator ( model_b          ,
                                      model_sb         ,
                                      dataset   = data ,
                                      ntoys_null = 50  ,
                                      ntoys_alt  = 50  ,
                                      ) 
        
        
        ## create Hypo Test inverter 
        hti = HypoTestInverter ( fc ,  0.90 , use_CLs = True , verbose = False )
        
        ## make a scan 
        hti .scan_with_progress  ( vrange ( 0 , 100 , 20 )  ) ## scan it!
 
    ## visualize the scan results 
    with use_canvas ( 'test_pointLimit: HypoTestInverter plot (frequentist)' , wait = 2 ) :
        plot = hti .plot
        plot .draw('LCb 2CL')    
        logger.info ( '90%%CL upper limit (frequentist) = %.1f' % hti.upper_limit )

# ============================================================================-
## Get the upper limit limit for small signal at fixed mass
#  - resolution is fixed
#  - Hybrid Calculator is used 
def test_point_limit_hc  () :
    """Get the upper limit at given point for small signal at fixed mass
    - resolution is fixed 
    - Hybrid Calculator is used 
    """

    logger = getLogger("test_point_limit_hc")

    logger.info ( "Test Point limits with RooStats using Hybrid Calculator" )

    ## if root_info < (6,24) :
    ##    logger.info ( 'Test is disabled for ROOT version %s' % str ( root_info ) )
    ##    return 

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             HybridCalculator      ,
                                             HypoTestInverter      )

    the_model = model.clone ( name = 'M3' )

    with use_canvas ( 'test_point_limit_hc' ) : 
        logger.info ( 'Dataset is\n%s' % data3.table ( prefix = '# ' ) ) 
        rr , frame = the_model.fitTo ( data3 , draw = True , nbins = 50 )

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf       = the_model   ,
                             poi       = the_model.S , ## parameter of interest 
                             dataset   = data3       ,
                             name      = 'S+B'       )
    
    model_sb.snapshot = the_model.S ## ATTENTION! 
    
    ## create ModelConfig  for 'B-only' model
    model_b  = ModelConfig ( pdf       = the_model          ,
                             poi       = the_model.S        , ## parameter of interest 
                             dataset   = data3              ,
                             workspace = model_sb.workspace , 
                             name      = 'B-only'           )
    
    the_model.S = 0 
    model_b.snapshot = the_model.S  ## ATTENTION! 
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )
    
    model_sb.mc.Print('vvv')
    model_b .mc.Print('vvv')

    ## with Hybrid calculator
    with timing ( "Using Hybrid Calculator" , logger = logger ) :
        
        ## create the calculator 
        hc  = HybridCalculator ( model_b          ,
                                 model_sb         ,
                                 dataset   = data ,
                                 ntoys_null = 50  ,
                                 ntoys_alt  = 50  ,
                                 ) 
        
        
        ## create Hypo Test inverter 
        hti = HypoTestInverter ( hc ,  0.90 , use_CLs = True , verbose = False )
        
        ## make a scan 
        hti .scan_with_progress  ( vrange ( 0 , 100 , 20 )  ) ## scan it!
 
    ## visualize the scan results 
    with use_canvas ( 'test_pointLimit: HypoTestInverter plot (hybrid)' , wait = 2 ) :
        plot = hti .plot
        plot .draw('LCb 2CL')    
        logger.info ( '90%%CL upper limit (hybrid) = %.1f' % hti.upper_limit )


# ============================================================================-
## Get the upper limit limit for small signal at fixed mass
#  - resolution is fixed
#  - Profile Likelihoood Calculator is used 
def test_point_limit_pl () :
    """Get the upper limit at given point for small signal at fixed mass
    - resoltuion is fixed 
    - Profile-Likelihood Calculator is used 
    """

    logger = getLogger("test_point_limit_pl")

    logger.info ( "Test Point limits with RooStats (Profile Likeloihood Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig                 ,
                                             ProfileLikelihoodCalculator ,
                                             HypoTestInverter            )

    the_model = model.clone ( name = 'M1' )

    with use_canvas ( 'test_point_limit_ac' ) : 
        logger.info ( 'Dataset is\n%s' % data4.table ( prefix = '# ' ) ) 
        rr , frame = the_model.fitTo ( data4 , draw = True , nbins = 50 )

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf       = the_model   ,
                             poi       = the_model.S , ## parameter of interest 
                             dataset   = data4       ,
                             name      = 'S+B'       )
    
    model_sb.snapshot = the_model.S ## ATTENTION! 
    
    ## create ModelConfig  for 'B-only' model
    model_b  = ModelConfig ( pdf       = the_model          ,
                             poi       = the_model.S        , ## parameter of interest 
                             dataset   = data4              ,
                             workspace = model_sb.workspace , 
                             name      = 'B-only'           )
    
    the_model.S = 0 
    model_b.snapshot = the_model.S  ## ATTENTION! 
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )
    
    model_sb.mc.Print('vvv')
    model_b .mc.Print('vvv')

    with timing ( "Using Profile Likelihood Calculator" , logger = logger ) :
        ## create the calculator 
        pl  = ProfileLikelihoodCalculator ( model_sb            ,
                                            dataset     = data  ,
                                            null_params = { the_model.S : 0.0 } )

        res = pl.calculator.GetHypoTest ()
        res.Print('vvv')
        
                        
    ##     ## create Hypo Test inverter 
    ##     hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )


        
    ##     ## make a scan 
    ##     hti .scan ( vrange ( 0 , 150 , 50 )  ) ## scan it!
        
    ## ## visualize the scan results 
    ## with use_canvas ( 'test_pointLimit: HypoTestInverter plot (asymptotic)' , wait = 2 ) :
    ##     plot = hti.plot
    ##     plot .draw('LCb 2CL')                    
    ##     logger.info ( '90%%CL upper limit (asymptotic)  = %.1f' % hti.upper_limit )


# =============================================================================
## Get the upper limit limit for small signal at fixed mass 
#  - resolution is known with some finite precision 
def test_point_limit2 () :
    """Get the upperlimit limit for small signal at fixed mass 
    - resolution is known with some finite precision 
    """
    
    logger = getLogger("test_point_limit2")

    logger.info ( "Test Point limits with RooStats (constrained resolution)" )

    from   ostap.fitting.roostats   import ( ModelConfig          ,
                                             AsymptoticCalculator ,
                                             HypoTestInverter     )

    sigma = ROOT.RooRealVar( 'sigma_Gau1', 'sigma of Gaussian' , 0.3 , 0.1 , 2 )

    the_signal = signal.clone ( sigma = sigma , name = 'S1' )

    ## create "soft" constraint for sigma 
    sigma_constraint = the_signal.soft_constraint ( sigma , VE ( 0.3 , 0.01**2 ) ) 
    
    the_model = model.clone ( name = 'M5' , signal = the_signal , signals = () )

    ## all constraints 
    constraints = sigma_constraint, 

    with use_canvas ( 'test_point_limit2' ) : 
        logger.info ( 'Dataset is\n%s' % data5.table ( prefix = '# ' ) ) 
        rr , frame = the_model.fitTo ( data5 , draw = True , nbins = 50 , constraints = constraints )
        
    model_sb = ModelConfig ( pdf         = the_model        ,
                             poi         = the_model.S      , ## parameter of interest 
                             dataset     = data5            ,
                             constraints = constraints      ,   
                             name        = 'S+B'            )
    
    model_sb.snapshot = the_model.S ## ATTENTION! 
    
    
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = the_model.S        , ## parameter of interest 
                             dataset     = data5              ,
                             workspace   = model_sb.workspace , 
                             constraints = constraints        ,   
                             name        = 'B-only'           )
    
    the_model.S = 0 
    model_b.snapshot = the_model.S  ## ATTENTION! 
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )

    model_sb.mc.Print('vvv')
    model_b .mc.Print('vvv')
    
    ac  = AsymptoticCalculator ( model_b           ,
                                 model_sb          ,
                                 dataset   = data5 , 
                                 one_sided = True  )     
    hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
    
    hti.scan ( vrange ( 0 , 150 , 50 )  ) ## scan it!

    with use_canvas ( 'test_point_limit2: HypoTestInverter plot' , wait = 2 ) :
        plot = hti.plot
        plot.draw('LCb 2CL')    
        logger.info ( '90%%CL upper limit = %.1f' % hti.upper_limit )


# ============================================================================-
## Get the upper limit limit for small signal at fixed mass 
#  - resolution is known with some finite precision 
#  - efficiency is known with some finite precision 
def test_point_limit3 () :
    """Get the upper limit limit for small soignal at fixed mass 
    - resolution is known with some finite precision 
    - efficiency is known with some finite precision 
    """
    logger = getLogger("test_point_limit3")

    logger.info ( "Test Point limits with RooStats (constrained resolution&efficiency)" )

    from   ostap.fitting.roostats   import ( ModelConfig          ,
                                             AsymptoticCalculator ,
                                             HypoTestInverter     )

    sigma = ROOT.RooRealVar( 'sigma_Gau2', 'sigma of Gaussian' , 0.3 , 0.1 , 2 )

    the_signal = signal.clone ( sigma = sigma , name = 'S2' )

    ## create "soft" constraint for sigma: (0.30+/-0.01)  
    sigma_constraint = the_signal.soft_constraint ( sigma , VE ( 0.3 , 0.01**2 ) ) 

    ## efficiency 
    eff = ROOT.RooRealVar( 'eff', 'efficiency'                        , 0.9 , 0 , 1)

    ## true (efficiency corrected) signal yield 
    NS  = ROOT.RooRealVar( 'NS' , 'efficiency corrected signal yield' , 0 , 200 )

    ## create "soft" constraint for efficiency: (90+/-3)%  
    eff_constraint = the_signal.soft_constraint ( eff , VE ( 0.9 , 0.02**2 ) ) 

    ## raw/visible signal yield 
    raw_S = the_signal.vars_multiply ( NS , eff , 'raw_S' , 'raw/observed signal yeild' )

    the_model = model.clone ( name = 'M6' , signal = the_signal , signals = () , S = raw_S )

    constraints = sigma_constraint, eff_constraint 

    with use_canvas ( 'test_point_limit3' ) : 
        logger.info ( 'Dataset is\n%s' % data6.table ( prefix = '# ' ) ) 
        rr , frame = the_model.fitTo ( data6 , draw = True , nbins = 50 , constraints = constraints )
        
    model_sb = ModelConfig ( pdf         = the_model        ,
                             poi         = NS               , ## parameter of interest 
                             dataset     = data6            ,
                             constraints = constraints      ,   
                             name        = 'S+B'            )
    
    model_sb.snapshot = NS ## ATTENTION! 
    
    
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = NS                 , ## parameter of interest 
                             dataset     = data6              ,
                             workspace   = model_sb.workspace , 
                             constraints = constraints        ,   
                             name        = 'B-only'           )
    
    NS.setVal ( 0 )
    model_b.snapshot = NS   ## ATTENTION! 
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )

    model_sb.mc.Print('vvv')
    model_b .mc.Print('vvv')
    
    ac  = AsymptoticCalculator ( model_b           ,
                                 model_sb          ,
                                 dataset   = data6 , 
                                 one_sided = True  )     
    hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
    
    hti.scan ( vrange ( 0 , 150 , 50 )  ) ## scan it!

    with use_canvas ( 'test_point_limit3: HypoTestInverter plot' , wait = 2 ) :
        plot = hti.plot
        plot.draw('LCb 2CL')    
        logger.info ( '90%%CL upper limit = %.1f' % hti.upper_limit )


    
# =============================================================================
if '__main__' == __name__ :

    from ostap.core.core    import rooSilent
    
    with rooSilent ( ) : 
    
        test_point_limit_ac ()        
        test_point_limit_fc ()
        test_point_limit_hc ()
        
        ## test_point_limit_pl ()
        
        test_point_limit2   ()
        test_point_limit3   ()

# =============================================================================
##                                                                      The END 
# =============================================================================
