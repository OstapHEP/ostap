#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_roostats3.py
# Test module for RooStats 
# ============================================================================= 
""" Test module for RooStats 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.meta_info     import root_info 
from   ostap.fitting.variables  import SETVAR, FIXVAR  
from   ostap.core.core          import cpp, VE, dsID, hID , rooSilent, Ostap 
from   ostap.utils.timing       import timing
from   ostap.utils.ranges       import vrange
from   ostap.utils.progress_bar import progress_bar
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.root_utils   import batch_env 
import ostap.fitting.models     as     Models
import ostap.logger.table       as     T
import ostap.fitting.roofit 
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
## set batch from environment 
batch_env ( logger )
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

graphs = [] 
# ============================================================================-
## Scan the position of the peak and get the limit for each peak posiiton peak
#  - resolution is known with some finite precision 
#  - efficiency is known with some finite precision 
def test_scan_limit1 () :
    """Scan the position of the peak and get the limit for each peak posiiton peak
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
        rS , _  = the_model.fitTo ( data1 , draw = True , nbins = 50 , constraints = constraints )
        
    model_sb = ModelConfig ( pdf         = the_model          ,
                             poi         = NS                 , ## parameter of interest 
                             dataset     = data1              ,
                             constraints = constraints        ,   
                             name        = 'S+B'              ,
                             snapshot    = rS                 )

    with FIXVAR ( NS ) :
        NS.setVal ( 0 ) 
        rB , _ = the_model.fitTo ( data , constraints = constraints )
        
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = NS                 , ## parameter of interest 
                             dataset     = data1              ,
                             workspace   = model_sb.workspace , 
                             constraints = constraints        ,   
                             name        = 'B-only'           ,
                             snapshot    = rB                 )
    
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )

    ## scan peak positions
    position = model_b.var ( the_signal.mean )

    ## prepare Brasil-plot 
    bplot = BrasilBand ( sigmas = (1,2,3) ) 
    
    rows = [ ( 'm0' , '90%UL' ) ]
    
    for m0 in progress_bar ( vrange ( 1 , 9 , 25 ) ) : 
        
        with SETVAR ( position ) :
            
            position.setVal ( m0  )
            
            ac  = AsymptoticCalculator ( model_b           ,
                                         model_sb          ,
                                         dataset   = data1 )
            
            ac.calculator.SetOneSided ( True )
            
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
        
        graphs.append ( bplot ) 
# ============================================================================-
## Scan the position of the peak and get the limit for each peak posiiton peak
#  - mass-dependent resolution is known with some finite precision
#  - mass-dependent efficiency is known with some finite precision 
def test_scan_limit2 () :
    """ Scan the positon of the peak and get the limit for each peak posiiton peak
    - mass-dependent resolution is known with some finite precision 
    - mass-depondent efficiency is known with some finite precision 
    """
    logger = getLogger("test_scan_limit2")

    logger.info ( "Scan limits with RooStats (constrained mass-dependent resolution&efficiency)" )


    if ( 6 , 28 , 0 ) <= root_info < ( 6 , 28 , 10 ) :
        logger.warning ( "Test is disabled foe this version of ROOT" )
        return
    
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
        r_sb , frame = the_model.fitTo    ( data3 , draw = True , nbins = 50 , constraints = constraints )

    ## fit with "B-only" model
    with use_canvas ( 'test_scan_limit2: B-only' ) :
        with FIXVAR ( NS ) :
            NS.setVal(0) 
            r_b , frame = the_model.fitTo ( data3 , draw = True , nbins = 50 , constraints = constraints )

    ## Create ModelConfig for "S+B" model 
    model_sb = ModelConfig ( pdf         = the_model        ,
                             poi         = NS               , ## parameter of interest 
                             dataset     = data3            ,
                             constraints = constraints      ,   
                             name        = 'S+B'            ,
                             snapshot    = r_sb             ) ## ATTENTION HERE!
    
    ## Create ModelConfig for "B-only" model 
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = NS                 , ## parameter of interest 
                             dataset     = data3              ,
                             workspace   = model_sb.workspace , 
                             constraints = constraints        ,   
                             name        = 'B-only'           ,
                             snapshot    = r_b                ) ## ATTENTION! 
    
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )
    
    ## peak positions  (as in workspace)
    position = model_b.var ( the_signal.mean )

    ## prepare Brasil-plot 
    bplot = BrasilBand ( sigmas =  ( 1 , 2 , 3 ) ) 
    
    rows = [ ( 'm0' , '90%UL' ) ]

    ## start the scan: 
    for m0 in progress_bar ( vrange ( 1 , 9 , 25 ) ) : 
        
        with SETVAR ( position ) :
            
            position.setVal ( m0  )

            ## 1) create the calculator 
            ac  = AsymptoticCalculator ( model_b           ,
                                         model_sb          ,
                                         dataset   = data3 )
            
            ac.calculator.SetOneSided   ( True )
            ac.calculator.SetPrintLevel ( -1   )

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

    model_b.ws.Print('vvv')
    
    ## visualize the Brasil plot 
    with use_canvas ( 'test_scan_limit2: Brasil plot' , wait = 2 ) :
        bplot.plot  .draw ( 'a' )
        bplot.legend.draw (     )

        title = '90% Upper limit'
        table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lr' )
        logger.info ( '%s:\n%s' % ( title , table ) )

        graphs.append ( bplot ) 

# ============================================================================-
## Scan the position of the peak and get the p0 
#  - resolution is known with some finite precision 
#  - efficiency is known with some finite precision
#  @thanks to Dima Golubkov 
def test_scan_p0_1 () :
    """ Scan the position of the peak and get p0 for each peak position peak
    - resolution is known with some finite precision 
    - efficiency is known with some finite precision 
    - thanks to Dima Golubkov 
    """
    logger = getLogger("test_scan_p0_1")

    logger.info ( "Scan p0 with RooStats (constrained resolution&efficiency)" )
    
    from   ostap.fitting.roostats   import ( ModelConfig          ,
                                             AsymptoticCalculator ,
                                             P0Plot               )

    sigma = ROOT.RooRealVar( 'sigma_Gau3', 'sigma of Gaussian' , 0.3 , 0.1 , 2 )

    the_signal = signal.clone ( sigma = sigma , name = 'S6' )

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

    the_model = model.clone ( name = 'M8' , signal = the_signal , signals = () , S = raw_S )

    constraints = sigma_constraint, eff_constraint 

    with use_canvas ( 'test_scan_p0_1' ) : 
        rr , frame = the_model.fitTo ( data4 , draw = True , nbins = 50 , constraints = constraints )
        
    model_sb = ModelConfig ( pdf         = the_model        ,
                             poi         = NS               , ## parameter of interest 
                             dataset     = data4            ,
                             constraints = constraints      ,   
                             name        = 'S+B'            )
    
    model_sb.snapshot = NS ## ATTENTION! 
    
    model_b  = ModelConfig ( pdf         = the_model          ,
                             poi         = NS                 , ## parameter of interest 
                             dataset     = data4              ,
                             workspace   = model_sb.workspace , 
                             constraints = constraints        ,   
                             name        = 'B-only'           )
    
    NS.setVal ( 0 )
    model_b.snapshot = NS   ## ATTENTION! 
    
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )

    ## scan peak positions
    position = model_b.var ( the_signal.mean )

    ## get p0-plot 
    plot = P0Plot()
    
    for m0 in progress_bar ( vrange ( 1 , 9 , 200  ) ) : 
        
        with SETVAR ( position ) :
            
            position.setVal ( m0  )
            
            ac  = AsymptoticCalculator ( model_sb          ,
                                         model_b           ,
                                         dataset   = data4 , 
                                         silent    = True  )

            ## attention!!! 
            ac.calculator.SetOneSidedDiscovery ( True )

            ## shortcut:
            # 
            plot.fill ( m0 , ac )

            ## full machinery:
            ## 
            ## ht    = ac.hypo_test 
            ## p0    = ht.     NullPValue () 
            ## p0alt = ht.AlternatePValue ()
            ## p0exp = ROOT.RooStats.AsymptoticCalculator.GetExpectedPValues (  p0 , p0alt , 0 , False ) 
            ## ## add values to the plot
            ## plot.fill ( m0 , p0 , p0exp )            
            ## del ht
            ## del ac
            
    ## visualize the P0 plot 
    with use_canvas  ( 'test_scan_p0_1: p0-plot'     , wait = 3 ) as cnv :
        plot.p0.draw ( 'ac' )
        cnv.SetLogy  ( True )
        plot.p0_expected.draw (  'c' )
        cnv.Update   ()
        
    with use_canvas ( 'test_scan_p0_1: #sigma-plot' , wait = 3 ) :
        plot.sigmas.draw ( 'ac' )

    title = 'p0-value scan'
    table = plot.table ( title = title , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )
        
    graphs.append ( plot ) 

# =============================================================================
if '__main__' == __name__ :

    from ostap.core.core    import rooSilent
    
    with rooSilent ( ) : 
    
        ## test_scan_limit1    ()
        
        ## test_scan_limit2    ()
        
        test_scan_p0_1      ()

# =============================================================================
##                                                                      The END 
# =============================================================================
