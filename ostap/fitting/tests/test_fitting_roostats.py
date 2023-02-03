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
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models
from   ostap.fitting.variables  import SETVAR 
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
    logger = getLogger ( 'test_fitting_roostats' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================



# ============================================================================
def test_intervals () :

    logger = getLogger("test_intervals ")

    logger.info ( "Test Intervals with RooStats" )

    from   ostap.fitting.roostats   import ( ProfileLikelihoodInterval ,
                                             FeldmanCousinsInterval    ,
                                             BayesianInterval          ,
                                             MCMCInterval              )


    x       = ROOT.RooRealVar('x','x-variable', -10 , 10 )
    sigma   = ROOT.RooFit.RooConst( 1 ) 
    pdf     = Models.Gauss_pdf    ( 'G' , x , mean = ( 0 , -1 , 1 ) , sigma  = sigma  ) 
    
    dataset = pdf.generate ( 100 )
    
    with use_canvas ( "test_intervals" ) : 
        result , frame = pdf.fitTo ( dataset , draw = True , nbins = 20 )

    rows = [  ( 'Method' , 'Interval [90%]' , 'Lower [95%]' , 'Upper[95%]' , 'CPU' ) ]
    
    with timing ( "Profile Likelihood", logger ) as t :
        pli = ProfileLikelihoodInterval ( pdf     = pdf      ,
                                          poi     = pdf.mean , 
                                          dataset = dataset  )
        low, high = pli.interval    ( 0.90 )
        upper     = pli.upper_limit ( 0.95 )
        lower     = pli.lower_limit ( 0.95 )
        
    row = 'Profile likelihood' , \
          '[%-+.3f,%+.3f]' %  ( low , high ) , \
          '%+.3f' %  lower , '%+.3f' %  upper , '%.1f' % t.delta
    rows.append ( row )

    with use_canvas  ( 'test_intervals: Profile Likelihood plot ' , wait = 1 ) :
        pli_plot = pli.plot ()
        if pli_plot : pli_plot.draw () 

    with use_canvas ( "test_limits" ) : 
        result , frame = pdf.fitTo ( dataset , draw = True , nbins = 20 )
        
    with timing ( 'Feldman-Cousins [~30",silent]'   , logger ) as t :
        with rooSilent ( ROOT.RooFit.WARNING ) : 
            fci = FeldmanCousinsInterval ( pdf     = pdf      ,
                                           poi     = pdf.mean , 
                                           dataset = dataset  )  
            low, high = fci.interval    ( 0.90 )
            upper     = fci.upper_limit ( 0.95 )
            lower     = fci.lower_limit ( 0.95 )
            
    row = 'Feldman-Cousins' , \
          '[%-+.3f,%+.3f]' %  ( low , high ) , \
          '%+.3f' %  lower , '%+.3f' %  upper , '%.1f' % t.delta
    rows.append ( row )

    with use_canvas  ( 'test_intervals: Felfman-Cousins plot ' , wait = 1 ) :
        fci_plot = fci.plot ()
        if fci_plot : fci_plot.draw ('ap') 

    with use_canvas ( "test_limits" ) : 
        result , frame = pdf.fitTo ( dataset , draw = True , nbins = 20 )
        
    with timing ( 'Bayesian [~60"]' , logger ) as t :        
        bci = BayesianInterval ( pdf     = pdf      ,
                                 poi     = pdf.mean , 
                                 dataset = dataset  ) 
        low, high = bci.interval    ( 0.90 )
        upper     = bci.upper_limit ( 0.95 )
        lower     = bci.lower_limit ( 0.95 )
        
    row = 'Bayesian' , \
          '[%-+.3f,%+.3f]' %  ( low , high ) , \
          '%+.3f' %  lower , '%+.3f' %  upper , '%.1f' % t.delta
    rows.append ( row )
    
    with use_canvas  ( 'test_intervals: Bayesian  plot ' , wait = 1 ) :
        bci_plot = bci.plot ()
        if bci_plot : bci_plot.draw () 

    with timing ( "Markov Chain MC (~6')" , logger ) as t :
        mci = MCMCInterval ( pdf     = pdf      ,
                             poi     = pdf.mean , 
                             dataset = dataset  )          
        low, high = mci.interval    ( 0.90 )
        upper     = mci.upper_limit ( 0.95 )
        lower     = mci.lower_limit ( 0.95 )


    row = 'Markov Chain MC' , \
          '[%-+.3f,%+.3f]' %  ( low , high ) , \
          '%+.3f' %  lower , '%+.3f' %  upper , '%.1f' % t.delta
    rows.append ( row )

    with use_canvas  ( 'test_intervals: Markov Chain MC plot ' , wait = 1 ) :
        mci_plot = mci.plot ()
        if mci_plot : mci_plot.draw () 

    with use_canvas  ( 'test_intervals: Visualise CL intervals' , wait = 2 ) as cnv :
        
        cnv.Divide(2,2)
        
        cnv.cd(1) 
        pli_plot = pli.plot ()
        if pli_plot : pli_plot.draw () 
        
        cnv.cd(2) 
        fci_plot = fci.plot ()
        if fci_plot : fci_plot.draw ('ap') 
        
        cnv.cd(3) 
        bci_plot = bci.plot ()
        if bci_plot : bci_plot.draw () 
        
        cnv.cd(4) 
        mci_plot = mci.plot ()
        if mci_plot : mci_plot.draw () 
        
        title = 'Intervals & Limits'
        table = T.table ( rows , title = title , prefix = '# ' , alignment = 'rcccr' )
        logger.info ( '%s\n%s' % ( title , table ) )

# =============================================================================
## more complicated/sophisticated  stuff
# =============================================================================

mass    = ROOT.RooRealVar   ('mass','mass-variable', 0 , 10 )
signal  = Models.Gauss_pdf  ( 'Gauss',
                              xvar  = mass                 ,
                              mean  = ( 2.5  , 0.5 , 9.5 ) ,
                              sigma = ( 0.01 , 0.3 , 3.0 ) )
signal.mean .fix()
signal.sigma.fix()
model   = Models.Fit1D ( signal = signal , background = 'e-' )
model.background.tau = -0.25 
model.S = 55
model.B = 1000
model.S.setMax(200)

data    = model.generate ( 55 + 1000 )

# ============================================================================-
## Get the limit at given point 
def test_point_limit () :
    """Get the limit at given point
    """

    logger = getLogger("test_point_limit")

    logger.info ( "Test Point limits with RooStats" )

    from   ostap.fitting.roostats   import ( ModelConfig          ,
                                             AsymptoticCalculator ,
                                             HypoTestInverter     )

    the_model = model.clone ( name = 'M1' )

    with use_canvas ( 'test_point_limit' ) : 
        rr , frame = the_model.fitTo ( data , draw = True , nbins = 50 )

        
    model_sb = ModelConfig ( pdf       = the_model   ,
                             poi       = the_model.S , ## parameter of interest 
                             dataset   = data        ,
                             name      = 'S+B'       )
    
    model_sb.mc.Print('vvv') 
    model_sb.snapshot = the_model.S ## ATTENTION! 
    
    
    model_b  = ModelConfig ( pdf       = the_model          ,
                             poi       = the_model.S        , ## parameter of interest 
                             dataset   = data               ,
                             workspace = model_sb.workspace , 
                             name      = 'B-only'           )
    
    the_model.S = 0 
    model_b.snapshot = the_model.S  ## ATTENTION! 
    
    model_sb.mc.Print('vvv')
    model_b .mc.Print('vvv')
    
    ac  = AsymptoticCalculator ( model_b          ,
                                 model_sb         ,
                                 dataset   = data , 
                                 one_sided = True )     
    hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
    
    hti.scan ( vrange ( 0 , 150 , 50 )  ) ## scan it!

    with use_canvas ( 'test_pointLimit: HypoTestInverter plot' , wait = 2 ) :
        plot = hti.plot
        plot.draw('LCb 2CL')    
        logger.info ( '90%%CL upper limit = %.1f' % hti.upper_limit )


# ============================================================================-
## Get the limit at given point
#  - resolution is known with soem finite precvision 
def test_point_limit2 () :
    """Get the limit at given point
    - resolution is known with soem finite precvision 
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
    
    the_model = model.clone ( name = 'M2' , signal = the_signal , signals = () )

    constraints = sigma_constraint, 

    with use_canvas ( 'test_point_limit2' ) : 
        rr , frame = the_model.fitTo ( data , draw = True , nbins = 50 , constraints = constraints )
        
    model_sb = ModelConfig ( pdf         = the_model        ,
                             poi         = the_model.S      , ## parameter of interest 
                             dataset     = data             ,
                             constraints = constraints      ,   
                             name        = 'S+B'            )
    
    model_sb.mc.Print('vvv') 
    model_sb.snapshot = the_model.S ## ATTENTION! 
    
    
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = the_model.S        , ## parameter of interest 
                             dataset     = data               ,
                             workspace   = model_sb.workspace , 
                             constraints = constraints        ,   
                             name        = 'B-only'           )
    
    the_model.S = 0 
    model_b.snapshot = the_model.S  ## ATTENTION! 
    
    model_sb.mc.Print('vvv')
    model_b .mc.Print('vvv')
    
    ac  = AsymptoticCalculator ( model_b          ,
                                 model_sb         ,
                                 dataset   = data , 
                                 one_sided = True )     
    hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
    
    hti.scan ( vrange ( 0 , 150 , 50 )  ) ## scan it!

    with use_canvas ( 'test_point_limit2: HypoTestInverter plot' , wait = 2 ) :
        plot = hti.plot
        plot.draw('LCb 2CL')    
        logger.info ( '90%%CL upper limit = %.1f' % hti.upper_limit )


# ============================================================================-
## Get the limit at given point
#  - resolution is known with some finite precvision 
#  - efficiency is known with some finite precision 
def test_point_limit3 () :
    """Get the limit at given point
    - resolution is known with some finite precvision 
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

    the_model = model.clone ( name = 'M3' , signal = the_signal , signals = () , S = raw_S )

    constraints = sigma_constraint, eff_constraint 

    with use_canvas ( 'test_point_limit3' ) : 
        rr , frame = the_model.fitTo ( data , draw = True , nbins = 50 , constraints = constraints )
        
    model_sb = ModelConfig ( pdf         = the_model        ,
                             poi         = NS               , ## parameter of interest 
                             dataset     = data             ,
                             constraints = constraints      ,   
                             name        = 'S+B'            )
    
    model_sb.mc.Print('vvv') 
    model_sb.snapshot = NS ## ATTENTION! 
    
    
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = NS                 , ## parameter of interest 
                             dataset     = data               ,
                             workspace   = model_sb.workspace , 
                             constraints = constraints        ,   
                             name        = 'B-only'           )
    
    NS.setVal ( 0 )
    model_b.snapshot = NS   ## ATTENTION! 
    
    model_sb.mc.Print('vvv')
    model_b .mc.Print('vvv')
    
    ac  = AsymptoticCalculator ( model_b          ,
                                 model_sb         ,
                                 dataset   = data , 
                                 one_sided = True )     
    hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
    
    hti.scan ( vrange ( 0 , 150 , 50 )  ) ## scan it!

    with use_canvas ( 'test_point_limit3: HypoTestInverter plot' , wait = 2 ) :
        plot = hti.plot
        plot.draw('LCb 2CL')    
        logger.info ( '90%%CL upper limit = %.1f' % hti.upper_limit )


# ============================================================================-
## Get the limit for different peak masses 
#  - resolution is known with some finite precvision 
#  - efficiency is known with some finite precision 
def test_scan_limit1 () :
    """Get the limit for different peak masses 
    - resolution is known with some finite precvision 
    - efficiency is known with some finite precision 
    """
    logger = getLogger("test_scan_limit1")

    logger.info ( "Scan limits with RooStats (constrained resolution&efficiency)" )

    from   ostap.fitting.roostats   import ( ModelConfig          ,
                                             AsymptoticCalculator ,
                                             HypoTestInverter     , 
                                             BrasilBand           )

    sigma = ROOT.RooRealVar( 'sigma_Gau3', 'sigma of Gaussian' , 0.3 , 0.1 , 2 )

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

    the_model = model.clone ( name = 'M3' , signal = the_signal , signals = () , S = raw_S )

    constraints = sigma_constraint, eff_constraint 

    with use_canvas ( 'test_scan_limit1' ) : 
        rr , frame = the_model.fitTo ( data , draw = True , nbins = 50 , constraints = constraints )
        
    model_sb = ModelConfig ( pdf         = the_model        ,
                             poi         = NS               , ## parameter of interest 
                             dataset     = data             ,
                             constraints = constraints      ,   
                             name        = 'S+B'            )
    
    model_sb.mc.Print('vvv') 
    model_sb.snapshot = NS ## ATTENTION! 
    
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = NS                 , ## parameter of interest 
                             dataset     = data               ,
                             workspace   = model_sb.workspace , 
                             constraints = constraints        ,   
                             name        = 'B-only'           )
    
    NS.setVal ( 0 )
    model_b.snapshot = NS   ## ATTENTION! 
    
    model_sb.mc.Print('vvv')
    model_b .mc.Print('vvv')

    
    ## scan peak positions
    position = model_b.var ( the_signal.mean )

    ## prepare Brasil-plot 
    bplot = BrasilBand ( sigmas = (1,2,3) ) 
    
    rows = [ ( 'm0' , '90%UL' ) ]
    
    for m0 in vrange ( 1 , 9 , 16 ) : 
        
        with SETVAR ( position ) :
            
            position.setVal ( m0  )
            
            ac  = AsymptoticCalculator ( model_b          ,
                                         model_sb         ,
                                         dataset   = data , 
                                         one_sided = True )
            ac.calculator.SetPrintLevel ( -1 )
            
            hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
            
            hti.scan ( vrange ( 0 , 150 , 50 )  ) ## scan it!
            
            interval = hti.interval
            limit    = hti.upper_limit

            ## fill Brasil plot 
            bplot.fill ( m0 , limit , interval ) 

            row = '%.1f' % m0 , '%.1f' % limit
            rows.append ( row ) 
                
    with use_canvas ( 'test_scan_limit1: Brasil plot' , wait = 2 ) :
        
        bplot.plot  .draw ( 'a' )
        bplot.legend.draw (     )

        title = '90% Upper limits'
        table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lr' )
        logger.info ( '%s:\n%s' % ( title , table ) )
        

# ============================================================================-
## Get the limit for different peak masses 
#  - mass-dependent resolution is known with some finite precvision
#  - mass-dependent efficiency is known with some finite precision 
def test_scan_limit2 () :
    """Get the limit for different peak masses 
    - mass-dependent resolution is known with some finite precvision 
    - mass-depondent efficiency is known with some finite precision 
    """
    logger = getLogger("test_scan_limit2")

    logger.info ( "Scan limits with RooStats (constrained mass-dependent resolution&efficiency)" )


    from ostap.fitting.roofuncs import BernsteinPoly as BP 

    ## resoltuion depends on the peak position 
    sigma_fun = BP ( 'SigmaFun' , signal.mean , power = 2 , pars = ( 1 , 4 , 1 ) )
    sigma_err = 0.01
    
    ## normalize at m0=2.5 to be 0.3
    vn = sigma_fun ( 2.5 )
    for p in sigma_fun.pars :
        p.setVal ( float ( p ) * 0.3 / vn )
        p.fix()
        
    ## efficiency depends on the peak positiob 
    eff_fun   = BP ( 'EffFun' , signal.mean , power = 1 , pars = ( 0.95 , 0.35 ) )
    eff_err   = 0.03 

    for p in eff_fun.pars : p.fix()

    with use_canvas ( 'test_scan_limit2: Resolurion depends on peak position' ) : sigma_fun.draw ()
    with use_canvas ( 'test_scan_limit2: Efficiency depends on peak position' ) : eff_fun  .draw ()

    from   ostap.fitting.roostats   import ( ModelConfig          ,
                                             AsymptoticCalculator ,
                                             HypoTestInverter     , 
                                             BrasilBand           )

    sigma = ROOT.RooRealVar( 'sigma', 'sigma of Gaussian' , 0.3 , 0.1 , 2 )

    the_signal = signal.clone ( sigma = sigma , name = 'S4' )

    ## create "soft" constraint for sigma: (0.30+/-0.01)  
    sigma_constraint = the_signal.soft_constraint ( sigma , value = sigma_fun.fun , error = sigma_err ) 

    ## efficiency 
    eff = ROOT.RooRealVar( 'eff', 'efficiency'                        , 0.9 , 0 , 1)

    ## true (efficiency corrected) signal yield 
    NS  = ROOT.RooRealVar( 'NS' , 'efficiency corrected signal yield' , 0 , 200 )

    ## create "soft" constraint for efficiency: (90+/-3)%  
    eff_constraint = the_signal.soft_constraint ( eff , value = eff_fun.fun , error = eff_err ) 

    ## raw/visible signal yield 
    raw_S = the_signal.vars_multiply ( NS , eff , 'raw_S' , 'raw/observed signal yeild' )

    the_model = model.clone ( name = 'M4' , signal = the_signal , signals = () , S = raw_S )

    constraints = sigma_constraint, eff_constraint 

    with use_canvas ( 'test_scan_limit2' ) : 
        rr , frame = the_model.fitTo ( data , draw = True , nbins = 50 , constraints = constraints )
        
    model_sb = ModelConfig ( pdf         = the_model        ,
                             poi         = NS               , ## parameter of interest 
                             dataset     = data             ,
                             constraints = constraints      ,   
                             name        = 'S+B'            )
    
    model_sb.mc.Print('vvv') 
    model_sb.snapshot = NS ## ATTENTION! 
    
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = NS                 , ## parameter of interest 
                             dataset     = data               ,
                             workspace   = model_sb.workspace , 
                             constraints = constraints        ,   
                             name        = 'B-only'           )
    
    NS.setVal ( 0 )
    model_b.snapshot = NS   ## ATTENTION! 
    
    model_sb.mc.Print('vvv')
    model_b .mc.Print('vvv')

    
    ## scan peak positions
    position = model_b.var ( the_signal.mean )

    ## prepare Brasil-plot 
    bplot = BrasilBand ( sigmas = (1,2,3) ) 
    
    rows = [ ( 'm0' , '90%UL' ) ]
    
    for m0 in vrange ( 1 , 9 , 16 ) : 
        
        with SETVAR ( position ) :
            
            position.setVal ( m0  )
            
            ac  = AsymptoticCalculator ( model_b          ,
                                         model_sb         ,
                                         dataset   = data , 
                                         one_sided = True )
            ac.calculator.SetPrintLevel ( -1 )
            
            hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
            
            hti.scan ( vrange ( 0 , 150 , 50 )  ) ## scan it!
            
            interval = hti.interval
            limit    = hti.upper_limit

            ## fill Brasil plot 
            bplot.fill ( m0 , limit , interval ) 

            row = '%.1f' % m0 , '%.1f' % limit
            rows.append ( row ) 
                
    with use_canvas ( 'test_scan_limit2: Brasil plot' , wait = 2 ) :
        bplot.plot  .draw ( 'a' )
        bplot.legend.draw (     )

        title = '90% Upper limits'
        table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lr' )
        logger.info ( '%s:\n%s' % ( title , table ) )
        

# =============================================================================
if '__main__' == __name__ :

    from ostap.core.core    import rooSilent
    with rooSilent ( ) : 
        
        ## test_intervals    ()
        
        ## test_point_limit  () 
        ## test_point_limit2 ()
        
        test_scan_limit1 ()
        
        ## test_scan_limit2 () 
        

# =============================================================================
##                                                                      The END 
# =============================================================================
