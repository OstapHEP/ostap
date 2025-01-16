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
N_S     = 20
N_B     = 1000
model.S = N_S 
model.B = N_B
model.S.setMax(100+3*N_S)

data    = model.generate ( N_S + N_B )

summary = [ ('method' , '90%CL' , 'time [s]') ]
plots   = []

# ============================================================================-
## Get the upper limit limit for small signal at fixed mass
#  - Resolution is fixed
#  - Asymptotic Calculator is used 
def test_point_limit_ac1() :
    """ Get the upper limit at given point for small signal at fixed mass
    - Resoltuion is fixed 
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_point_limit_ac1")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )

    the_model = model.clone ( name = 'M1' )
    
    with use_canvas ( 'test_point_limit_ac1' ) : 
        rS , _ = the_model.fitTo ( data , draw = True , nbins = 50 )

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf       = the_model   ,
                             poi       = the_model.S , ## parameter of interest 
                             dataset   = data        ,
                             name      = 'S+B'       ,
                             snapshot  = rS          ) ## ATTENTION! 
    
    with FIXVAR ( the_model.S ) :
        the_model.S = 0 
        rB , _ = the_model.fitTo ( data )

    ## create ModelConfig  for 'B-only' model
    model_b  = ModelConfig ( pdf       = the_model          ,
                             poi       = the_model.S        , ## parameter of interest 
                             dataset   = data               ,
                             workspace = model_sb.workspace , 
                             name      = 'B-only'           ,
                             snapshot  = rB                 )
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )
    
    with timing ( "Using Asymptotic Calculator" , logger = logger ) as timer :
        ## create the calculator 
        ac  = AsymptoticCalculator ( model_b           ,
                                     model_sb          ,
                                     dataset   = data  ,
                                     asimov    = False )
        ac.calculator.SetOneSided ( True ) 
        
        ## create Hypo Test inverter 
        hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
        
        ## make a scan 
        hti .scan_with_progress ( vrange ( 0 , 150 , 150 )  ) ## scan it!
        
    ## visualize the scan results 
    with use_canvas ( 'test_point_limit_ac1: HypoTestInverter plot (asymptotic)' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.1f' % hti.upper_limit )

    row = 'Asymptotic (Asimov=False)' , '%.1f' % hti.upper_limit , '%.1f' % timer.delta 
    summary.append ( row  )
    plots  .append ( plot )
        
    ## check the dataset
    stat = data.statVar('mass')
    if stat.rms() <= 0 :
        logger.error   ( 'Calculator destroyed input dataset!') 
        logger.error   ( 'Dataset is\n%s'          % data.table ( prefix = '# ' ) ) 

# ============================================================================-
## Get the upper limit limit for small signal at fixed mass
#  - Resolution is fixed
#  - Asymptotic Calculator is used 
def test_point_limit_ac2() :
    """ Get the upper limit at given point for small signal at fixed mass
    - Resolution is fixed 
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_point_limit_ac2")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator with Asimov dataset)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )

    the_model = model.clone ( name = 'M2' )
    
    with use_canvas ( 'test_point_limit_ac2' ) : 
        rS , _ = the_model.fitTo ( data , draw = True , nbins = 50 )
        
    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf       = the_model   ,
                             poi       = the_model.S , ## parameter of interest 
                             dataset   = data        ,
                             name      = 'S+B'       ,
                             snapshot  = rS          ) ## ATTENTION!
    
    with FIXVAR ( the_model.S ) :
        the_model.S = 0 
        rB , _ = the_model.fitTo ( data )
    
    ## create ModelConfig  for 'B-only' model
    model_b  = ModelConfig ( pdf       = the_model          ,
                             poi       = the_model.S        , ## parameter of interest 
                             dataset   = data               ,
                             workspace = model_sb.workspace , 
                             name      = 'B-only'           ,
                             snapshot  = rB                 )
  
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )
    
    with timing ( "Using Asymptotic+Asimov Calculator" , logger = logger ) as timer :
        ## create the calculator 
        ac  = AsymptoticCalculator ( model_b           ,
                                     model_sb          ,
                                     dataset   = data  ,
                                     asimov    = True  )
        ac.calculator.SetOneSided ( True ) 
        
        ## create Hypo Test inverter 
        hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
        
        ## make a scan 
        hti .scan_with_progress ( vrange ( 0 , 150 , 150 )  ) ## scan it!
        
    ## visualize the scan results 
    with use_canvas ( 'test_point_limit_ac2: HypoTestInverter plot (asymptotic/asimov)' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic/asimov)  = %.1f' % hti.upper_limit )

        row = 'Asymptotic (Asimov=True)' , '%.1f' % hti.upper_limit, '%.1f' % timer.delta 
        summary.append ( row  )
        plots  .append ( plot )
        
    ## check the dataset
    stat = data.statVar('mass')
    if stat.rms() <= 0 :
        logger.error   ( 'Calculator destroyed input dataset!') 
        logger.error   ( 'Dataset is\n%s'          % data.table ( prefix = '# ' ) ) 

# ============================================================================-
## Get the upper limit limit for small signal at fixed mass
#  - Resolution is fixed
#  - Frequestist Calculator is used 
def test_point_limit_fc  () :
    """ Get the upper limit at given point for small signal at fixed mass
    - Resolution is fixed 
    - Frequestist Calculator is used 
    """

    logger = getLogger("test_point_limit_fc")

    logger.info ( "Test Point limits with RooStats using Frequestist Calculator" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             FrequentistCalculator ,
                                             HypoTestInverter      )

    the_model = model.clone ( name = 'M3' )
    
    with use_canvas ( 'test_point_limit_fc' ) : 
        rS , _ = the_model.fitTo ( data , draw = True , nbins = 50 )

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf       = the_model   ,
                             poi       = the_model.S , ## parameter of interest 
                             dataset   = data        ,
                             name      = 'S+B'       ,
                             snapshot  = rS          ) 
    
    with FIXVAR ( the_model.S ) :
        the_model.S = 0 
        rB , _ = the_model.fitTo ( data )
        
    ## create ModelConfig  for 'B-only' model
    model_b  = ModelConfig ( pdf       = the_model          ,
                             poi       = the_model.S        , ## parameter of interest 
                             dataset   = data               ,
                             workspace = model_sb.workspace , 
                             name      = 'B-only'           ,
                             ### snapshot  = rB                 )
                             snapshot  = the_model.S ) 
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )
    
    ## with Frequentist calculator
    with timing ( "Using Frequentist Calculator" , logger = logger ) as timer :
        
        ## create the calculator 
        fc  = FrequentistCalculator ( model_b            ,
                                      model_sb           ,
                                      dataset    = data  ,
                                      ntoys_null = 100   ,
                                      ntoys_alt  = 100   ,
                                      ) 
        
        ## create Hypo Test inverter 
        hti = HypoTestInverter ( fc ,  0.90 , use_CLs = True , verbose = False )
        
        ## make a scan 
        hti .scan_with_progress  ( vrange ( 0.1 , 100 , 10 )  ) ## scan it!
 
    ## visualize the scan results 
    with use_canvas ( 'test_point_limits_fc: HypoTestInverter plot (frequentist)' , wait = 2 ) :
        plot = hti .plot
        plot .draw('LCb 2CL')    
        logger.info ( '90%%CL upper limit (frequentist) = %.1f' % hti.upper_limit )

        row = 'Frequentist' , '%.1f' % hti.upper_limit, '%.1f' % timer.delta
        summary.append ( row  )
        plots  .append ( plot )

    ## check the dataset
    stat = data.statVar('mass')
    if stat.rms() <= 0 :
        logger.error   ( 'Calculator destroyed input dataset!') 
        logger.error   ( 'Dataset is\n%s'          % data.table ( prefix = '# ' ) ) 

# ============================================================================-
## Get the upper limit limit for small signal at fixed mass
#  - Resolution is fixed
#  - Hybrid Calculator is used 
def test_point_limit_hc  () :
    """ Get the upper limit at given point for small signal at fixed mass
    - Resolution is fixed 
    - Hybrid Calculator is used 
    """

    logger = getLogger("test_point_limit_hc")

    logger.info ( "Test Point limits with RooStats using Hybrid Calculator" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             HybridCalculator      ,
                                             HypoTestInverter      )

    the_model = model.clone ( name = 'M4' )
    
    with use_canvas ( 'test_point_limit_hc' ) : 
        rS , _ = the_model.fitTo ( data , draw = True , nbins = 50 )

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf       = the_model   ,
                             poi       = the_model.S , ## parameter of interest 
                             dataset   = data        ,
                             name      = 'S+B'       , 
                             snapshot  = rS          ) ## attention
    
    with FIXVAR ( the_model.S ) :
        the_model.S = 0 
        rB , _ = the_model.fitTo ( data )
        
    ## create ModelConfig  for 'B-only' model
    model_b  = ModelConfig ( pdf       = the_model          ,
                             poi       = the_model.S        , ## parameter of interest 
                             dataset   = data               ,
                             workspace = model_sb.workspace , 
                             name      = 'B-only'           ,
                             ### snapshot  = rB                 )
                             snapshot  = the_model.S ) 
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )
    
    ## with Hybrid calculator
    with timing ( "Using Hybrid Calculator" , logger = logger ) as timer :
        
        ## create the calculator 
        hc  = HybridCalculator ( model_b            ,
                                 model_sb           ,
                                 dataset    = data  ,
                                 ntoys_null = 100   ,
                                 ntoys_alt  = 100   ,
                                 ) 
        
        
        ## create Hypo Test inverter 
        hti = HypoTestInverter ( hc ,  0.90 , use_CLs = True , verbose = False )
        
        ## make a scan 
        hti .scan_with_progress  ( vrange ( 0.1 , 100 , 10 )  ) ## scan it!
 
    ## visualize the scan results 
    with use_canvas ( 'test_point_limit_hc: HypoTestInverter plot (hybrid)' , wait = 2 ) :
        plot = hti .plot
        plot .draw('LCb 2CL')    
        logger.info ( '90%%CL upper limit (hybrid) = %.1f' % hti.upper_limit )

        row = 'Hybrid' , '%.1f' % hti.upper_limit, '%.1f' % timer.delta
        summary.append ( row  )
        plots  .append ( plot )

    ## check the dataset
    stat = data.statVar('mass')
    if stat.rms() <= 0 :
        logger.error   ( 'Calculator destroyed input dataset!') 
        logger.error   ( 'Dataset is\n%s'          % data.table ( prefix = '# ' ) ) 

# ============================================================================-
## Get the upper limit limit for small signal at fixed mass
#  - resolution is fixed
#  - Profile Likelihoood Calculator is used 
def test_point_limit_pl () :
    """ Get the upper limit at given point for small signal at fixed mass
    - resolution is fixed 
    - Profile-Likelihood Calculator is used 
    """

    logger = getLogger("test_point_limit_pl")

    logger.info ( "Test Point limits with RooStats (Profile Likelihood Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig                 ,
                                             ProfileLikelihoodCalculator ,
                                             HypoTestInverter            )

    the_model = model.clone ( name = 'M5' )
    
    with use_canvas ( 'test_point_limit_ac' ) : 
        rS , frame = the_model.fitTo ( data , draw = True , nbins = 50 )

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf       = the_model   ,
                             poi       = the_model.S , ## parameter of interest 
                             dataset   = data        ,
                             name      = 'S+B'       ,
                             snapshot  = rS          )
    
    with FIXVAR ( the_model.S ) :
        the_model.S = 0 
        rB , _ = the_model.fitTo ( data )
        
    ## create ModelConfig  for 'B-only' model
    model_b  = ModelConfig ( pdf       = the_model          ,
                             poi       = the_model.S        , ## parameter of interest 
                             dataset   = data               ,
                             workspace = model_sb.workspace , 
                             name      = 'B-only'           ,
                             snapshot  = rB                 )
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )
    
    with timing ( "Using Profile Likelihood Calculator" , logger = logger ) as timer :
        ## create the calculator 
        pl  = ProfileLikelihoodCalculator ( model_sb            ,
                                            dataset     = data  ,
                                            null_params = { the_model.S : 0.0 } )

        """
        res = pl.calculator.GetHypoTest ()
        
        ## create Hypo Test inverter 
        hti = HypoTestInverter ( pl ,  0.90 , use_CLs = True , verbose = False )

        ## make a scan 
        hti .scan_with_progress ( vrange ( 0 , 150 , 100 )  ) ## scan it!
        
        ## visualize the scan results 
        with use_canvas ( 'test_point_limi_pl: HypoTestInverter plot (ProfileLikelihood)' , wait = 2 ) :
            plot = hti.plot
            plot .draw('LCb 2CL')                    
            logger.info ( '90%%CL upper limit (profile_likelihood)  = %.1f' % hti.upper_limit )
            
        row = 'ProfileLikelihood' , '%.1f' % hti.upper_limit , '%.1f' % timer.delta 
        summary.append ( row  )
        """
        
# =============================================================================
## Get the upper limit limit for small signal at fixed mass 
#  - resolution is known with some finite precision 
def test_point_limit2 () :
    """ Get the upper limit limit for small signal at fixed mass 
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
    
    the_model = model.clone ( name = 'M6' , signal = the_signal , signals = () )

    ## all constraints 
    constraints = sigma_constraint, 

    with use_canvas ( 'test_point_limit2' ) : 
        rS , _ = the_model.fitTo ( data , draw = True , nbins = 50 , constraints = constraints )
        
    model_sb = ModelConfig ( pdf         = the_model        ,
                             poi         = the_model.S      , ## parameter of interest 
                             dataset     = data             ,
                             constraints = constraints      ,   
                             name        = 'S+B'            ,
                             snapshot    = rS               )
    
    with FIXVAR ( the_model.S ) :
        the_model.S = 0 
        rB , _ = the_model.fitTo ( data , constraints = constraints ) 

    model_b  = ModelConfig ( pdf         = the_model          ,
                             ## poi         = the_model.S        , ## parameter of interest
                             poi         = model_sb.poi       ,
                             dataset     = data               ,
                             workspace   = model_sb.workspace , 
                             constraints = constraints        ,   
                             name        = 'B-only'           ,
                             snapshot    = rB                 )
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )

    with timing ( "Using Asymptotic Calculator with 1 constraint" , logger = logger ) as timer :
        
        ac  = AsymptoticCalculator ( model_b           ,
                                     model_sb          ,
                                     dataset   = data  )
        ac.calculator.SetOneSided ( True ) 

        hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
    
        hti.scan_with_progress ( vrange ( 0 , 150 , 150 )  ) ## scan it!

    with use_canvas ( 'test_point_limit2: HypoTestInverter plot' , wait = 2 ) :
        plot = hti.plot
        plot.draw('LCb 2CL')    
        logger.info ( '90%%CL upper limit = %.1f' % hti.upper_limit )
        plots  .append ( plot )

    row = 'Asymptotic (1 constraint)' , '%.1f' % hti.upper_limit , '%.1f' % timer.delta 
    summary.append ( row  )

    ## check the dataset
    stat = data.statVar('mass')
    if stat.rms() <= 0 :
        logger.error   ( 'Calculator destroyed input dataset!') 
        logger.error   ( 'Dataset is\n%s'          % data.table ( prefix = '# ' ) ) 

# ============================================================================-
## Get the upper limit limit for small signal at fixed mass 
#  - Resolution is known with some finite precision 
#  - Efficiency is known with some finite precision 
def test_point_limit3 () :
    """ Get the upper limit limit for small signal at fixed mass 
    - Resolution is known with some finite precision 
    - Efficiency is known with some finite precision 
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

    ## create "soft" constraint for efficiency: (90+/-2)%  
    eff_constraint = the_signal.soft_constraint ( eff , VE ( 0.9 , 0.02**2 ) ) 

    ## raw/visible signal yield 
    raw_S = the_signal.vars_multiply ( NS , eff , 'raw_S' , 'raw/observed signal yeild' )

    the_model = model.clone ( name = 'M7' , signal = the_signal , signals = () , S = raw_S )

    constraints = sigma_constraint, eff_constraint 

    with use_canvas ( 'test_point_limit3' ) : 
        rS , _ = the_model.fitTo ( data , draw = True , nbins = 50 , constraints = constraints )
        
    model_sb = ModelConfig ( pdf         = the_model        ,
                             poi         = NS               , ## parameter of interest 
                             dataset     = data             ,
                             constraints = constraints      ,   
                             name        = 'S+B'            ,
                             snapshot    = rS               )
    
    with FIXVAR ( NS ) :
        NS.setVal ( 0 ) 
        rB , _ = the_model.fitTo ( data , constraints = constraints )
     
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = NS                 , ## parameter of interest 
                             dataset     = data               ,
                             workspace   = model_sb.workspace , 
                             constraints = constraints        ,   
                             name        = 'B-only'           ,
                             snapshot    = rB                 )
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )

    with timing ( "Using Asymptotic Calculator with 2 constraints" , logger = logger ) as timer :

        ac  = AsymptoticCalculator ( model_b           ,
                                     model_sb          ,
                                     dataset   = data  )
        ac.calculator.SetOneSided ( True )
                
        hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
        
        hti.scan_with_progress ( vrange ( 0 , 150 , 150 )  ) ## scan it!

    with use_canvas ( 'test_point_limit3: HypoTestInverter plot' , wait = 2 ) :
        plot = hti.plot
        plot.draw('LCb 2CL')    
        logger.info ( '90%%CL upper limit = %.1f' % hti.upper_limit )
        plots  .append ( plot )
        row = 'Asymptotic (2 constraints)' , '%.1f' % hti.upper_limit , '%.1f' % timer.delta 
        summary.append ( row  )

    ## check the dataset
    stat = data.statVar('mass')
    if stat.rms() <= 0 :
        logger.error   ( 'Calculator destroyed input dataset!') 
        logger.error   ( 'Dataset is\n%s'          % data.table ( prefix = '# ' ) ) 

# =============================================================================
if '__main__' == __name__ :

    from ostap.core.core    import rooSilent
    
    with rooSilent ( ) : 

        test_point_limit_ac1 ()
        test_point_limit_ac2 ()
        
        test_point_limit_fc  ()
        test_point_limit_hc  ()
        
        test_point_limit_pl  ()

        test_point_limit2    ()
        test_point_limit3    ()
        
    import ostap.logger.table as T
    title = 'Summary of 90%CL Upper Limits'
    table = T.table ( summary , title = title , prefix = '# ' , alignment = 'lcr' )
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
