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
    logger = getLogger ( 'test_fitting_roostats3' )
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

# ============================================================================-
## Scan the positon of the peak and get the limit for each peak posiiton peak
#  - resolution is known with some finite precision 
#  - efficiency is known with some finite precision 
def test_scan_limit1 () :
    """Scan the positon of the peak and get the limit for each peak posiiton peak
    - resolution is known with some finite precision 
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

    the_model = model.clone ( name = 'M7' , signal = the_signal , signals = () , S = raw_S )

    constraints = sigma_constraint, eff_constraint 

    with use_canvas ( 'test_scan_limit1' ) : 
        rr , frame = the_model.fitTo ( data , draw = True , nbins = 50 , constraints = constraints )
        
    model_sb = ModelConfig ( pdf         = the_model        ,
                             poi         = NS               , ## parameter of interest 
                             dataset     = data             ,
                             constraints = constraints      ,   
                             name        = 'S+B'            )
    
    model_sb.snapshot = NS ## ATTENTION! 
    
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = NS                 , ## parameter of interest 
                             dataset     = data               ,
                             workspace   = model_sb.workspace , 
                             constraints = constraints        ,   
                             name        = 'B-only'           )
    
    NS.setVal ( 0 )
    model_b.snapshot = NS   ## ATTENTION! 
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )

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
## Scan the positon of the peak and get the limit for each peak posiiton peak
#  - mass-dependent resolution is known with some finite precision
#  - mass-dependent efficiency is known with some finite precision 
def test_scan_limit2 () :
    """ Scan the positon of the peak and get the limit for each peak posiiton peak
    - mass-dependent resolution is known with some finite precision 
    - mass-depondent efficiency is known with some finite precision 
    """
    logger = getLogger("test_scan_limit2")

    logger.info ( "Scan limits with RooStats (constrained mass-dependent resolution&efficiency)" )

    from ostap.fitting.roofuncs import BernsteinPoly as BP 

    ## resolution depends on the peak position
    sigma_fun = BP ( 'SigmaFun' , signal.mean , power = 2 , pars = ( 1 , 4 , 1 ) )
    sigma_err = 0.01    
    ## normalize at m0 = 2.5 to be 0.3
    vn = sigma_fun ( 2.5 )
    for p in sigma_fun.pars :
        p.setVal ( float ( p ) * 0.3 / vn )
        p.fix()
    
    ## efficiency depends on the peak position
    eff_fun   = BP ( 'EffFun' , signal.mean , power = 1 , pars = ( 0.95 , 0.35 ) )
    eff_err   = 0.03 
    for p in eff_fun.pars : p.fix()

    
    with use_canvas ( 'test_scan_limit2: Resolution depends on the peak position' ) : sigma_fun.draw ()
    with use_canvas ( 'test_scan_limit2: Efficiency depends on the peak position' ) : eff_fun  .draw ()

    from   ostap.fitting.roostats   import ( ModelConfig          ,
                                             AsymptoticCalculator ,
                                             HypoTestInverter     , 
                                             BrasilBand           )
    
    ## resolution 
    sigma      = ROOT.RooRealVar( 'sigma', 'sigma of Gaussian' , 0.3 , 0.1 , 2 )

    ## signal 
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

    ## fit models 
    the_model = model.clone ( name = 'M8' , signal = the_signal , signals = () , S = raw_S )

    ## collect constraints 
    constraints = sigma_constraint, eff_constraint 

    ## fit with "S+B" model 
    with use_canvas ( 'test_scan_limit2: S+B'    ) : 
        r_sb , frame = the_model.fitTo    ( data , draw = True , nbins = 50 , constraints = constraints )

    ## fit with "B-only" model
    with use_canvas ( 'test_scan_limit2: B-only' ) :
        with FIXVAR ( NS ) :
            NS.setVal(0) 
            r_b , frame = the_model.fitTo ( data , draw = True , nbins = 50 , constraints = constraints )

    ## Create ModelConfig for "S+B" model 
    model_sb = ModelConfig ( pdf         = the_model        ,
                             poi         = NS               , ## parameter of interest 
                             dataset     = data             ,
                             constraints = constraints      ,   
                             name        = 'S+B'            ,
                             snapshot    = r_sb             ) ## ATTENTION HERE!
    
    ## model_sb.snapshot = NS ## ATTENTION! 
    ## model_sb.snapshot = r_sb ## ATTENTION! 
    
    ## Create ModelConfig for "B-only" model 
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = NS                 , ## parameter of interest 
                             dataset     = data               ,
                             workspace   = model_sb.workspace , 
                             constraints = constraints        ,   
                             name        = 'B-only'           ,
                             snapshot    = r_b                ) ## ATTENTION! 
    
    # NS.setVal ( 0 )
    # model_b.snapshot = NS   ## ATTENTION! 
    ## model_b.snapshot = r_b    ## ATTENTION! 
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )

    model_sb.mc.Print('vvv')
    model_b .mc.Print('vvv')
    
    ## peak positions  (as in workspace)
    position = model_b.var ( the_signal.mean )

    ## prepare Brasil-plot 
    bplot = BrasilBand ( sigmas = (1,2,3) ) 
    
    rows = [ ( 'm0' , '90%UL' ) ]

    ## start the scan: 
    for m0 in vrange ( 1 , 9 , 16 ) : 
        
        with SETVAR ( position ) :
            
            position.setVal ( m0  )

            ## 1) create the calculator 
            ac  = AsymptoticCalculator ( model_b          ,
                                         model_sb         ,
                                         dataset   = data , 
                                         one_sided = True )
            ac.calculator.SetPrintLevel ( -1 )

            ## 2) create the Hypo Test inverter 
            hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )

            ## define the scan 
            hti.scan ( vrange ( 0 , 150 , 50 )  ) ## scan it!

            ## get the obtained interval & upper limit 
            interval = hti.interval
            limit    = hti.upper_limit

            ## fill Brasil plot 
            bplot.fill ( m0 , limit , interval ) 

            ## fill a row in table 
            row = '%.1f' % m0 , '%.1f' % limit
            rows.append ( row ) 

    ## visualize the Brasin plot 
    with use_canvas ( 'test_scan_limit2: Brasil plot' , wait = 2 ) :
        bplot.plot  .draw ( 'a' )
        bplot.legend.draw (     )

        title = '90% Upper limits'
        table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lr' )
        logger.info ( '%s:\n%s' % ( title , table ) )
        
    print ( model_b.table() )

    
# =============================================================================
if '__main__' == __name__ :

    from ostap.core.core    import rooSilent
    
    with rooSilent ( ) : 
    
        test_scan_limit1    ()
        test_scan_limit2    () 

# =============================================================================
##                                                                      The END 
# =============================================================================
