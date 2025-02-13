#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_roostats5.py
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
from   ostap.utils.utils        import vrange, lrange
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait, batch_env 
from   ostap.fitting.simfit     import combined_data
import ostap.fitting.models     as     Models
import ostap.logger.table       as     T
import ostap.fitting.roofit 
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_roostats6' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

# =============================================================================
## more complicated/sophisticated  stuff
# =============================================================================

mass   = ROOT.RooRealVar   ( 'mass' , 'mass-variable', 0 , 10 )
vars   = ROOT.RooArgSet    ( mass   ) 
signal = Models.Gauss_pdf  ( 'Gauss',
                             xvar  = mass                 ,
                             mean  = ( 5.0 , 1.5  , 8.5 ) ,
                             sigma = ( 0.5 , 0.01 , 3.0 ) )

signal.sigma.fix()
signal.mean .fix()

model   = Models.Fit1D ( signal = signal , background = 'p1' , name = 'A' , suffix = 'A' )
model.S.setMin ( -30 )
model.S = 0 
model.B = 200 



dataset = model.generate ( 200 , sample = False )
r , f = model.fitTo ( dataset , draw = True , quiet = True ) 

err = 0.1  
reff       = ROOT.RooRealVar  ( 'reff'  , 'efficiency' , 1.0 , 1.e-5 , 2 )
ceff1      = ROOT.RooGaussian ( 'ceff1' , 'efficiency constraint', reff ,
                                ROOT.RooFit.RooConst ( 1.0  ) ,
                                ROOT.RooFit.RooConst ( err  ) )
reff_mean  = ROOT.RooRealVar ( 'ref_mean' , '' , 1.0 , 1.e-5 , 2 )
reff_sigma = ROOT.RooRealVar ( 'ref_mean' , '' , err , 1.e-6 , 1 )
ceff2      = ROOT.RooGaussian ( 'ceff2' , 'efficiency constraint',
                                reff_mean , reff , reff_sigma ) 
reff_mean .fix() 
reff_sigma.fix() 


SS       = ROOT.RooRealVar      ( 'SS' , 'signal' , 0 , -30 , 100 )
NS       = signal.vars_multiply (  SS , reff ) 

summary = [ ( 'method' , '90%CL' , 'time [s]') ]
plots   = [] 

values  = vrange ( 0 , 30 , 30 ) 
# ============================================================================-
## Get the upper limit limit for small signal
#  - Asymptotic Calculator is used 
def test_limit_1() :
    """ Get the upper limit at given point for small signal
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_limit_1")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )

    the_model = model.clone ( name = 'X1' )
    data      = dataset 

    logger.info ('Fit S+B') 
    with use_canvas ( 'test_limit_1' ) :
        rS , _ = the_model.fitTo ( data , draw = True , silent = True )
                   
    POI = the_model.S  

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf         = the_model   ,
                             poi         = POI         , ## parameter of interest 
                             dataset     = data        ,
                             name        = 'S+B'       ,
                             snapshot    = POI         ) ## ATTENTION! 
    
    with FIXVAR ( POI ) :
        POI.setVal ( 0 )
        logger.info ( 'Fit B-only' )        
        rB , _ = the_model.fitTo ( data , silent = True )

        # create ModelConfig  for 'B-only' model
        model_b  = ModelConfig ( pdf         = the_model          ,
                                 poi         = POI                , ## parameter of interest
                                 dataset     = data               ,
                                 workspace   = model_sb.workspace , 
                                 name        = 'B-only'           ,
                                 snapshot    = POI                )
      
       
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
        hti .scan_with_progress ( values  ) #can it!
          
    ## visualize the scan results 
    with use_canvas ( 'test_point_limit_ac1: No efficiency' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.2f' % hti.upper_limit )

    row = 'No efficiency' , '%.3f' % hti.upper_limit , '%.1f' % timer.delta 
    summary.append ( row  )
    plots  .append ( plot )

# ============================================================================-
## Get the upper limit limit for small signal
#  - Asymptotic Calculator is used 
def test_limit_2() :
    """ Get the upper limit at given point for small signal
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_point_limit_2")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )

    the_model = model.clone  ( S = NS , name = 'X2' )
    data      = dataset
    
    reff.setVal ( 1 ) 
    reff.fix()

    POI = SS
    
    logger.info ('Fit S+B') 
    with use_canvas ( 'test_limit_2' ) :
        rS , _ = the_model.fitTo ( data , draw = True , silent = True )

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf         = the_model   ,
                             poi         = POI         , ## parameter of interest
                             dataset     = data        ,
                             name        = 'S+B'       ,
                             snapshot    = POI         ) ## ATTENTION! 
    

    with FIXVAR ( POI ) :
        POI.setVal ( 0 )
        rB , _ = the_model.fitTo ( data , silent = True )

        # create ModelConfig  for 'B-only' model
        model_b  = ModelConfig ( pdf         = the_model          ,
                                 poi         = POI                , ## parameter of interest
                                 dataset     = data               ,
                                 workspace   = model_sb.workspace , 
                                 name        = 'B-only'           ,
                                 snapshot    = POI                )
      
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
        hti .scan_with_progress ( values ) #can it!
          
    ## visualize the scan results 
    with use_canvas ( 'test_point_limit_ac2: Fixed efficiency' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.2f' % hti.upper_limit )

    row = 'Fixed efficiency' , '%.3f' % hti.upper_limit , '%.1f' % timer.delta 
    summary.append ( row  )
    plots  .append ( plot )


# ============================================================================-
## Get the upper limit limit for small signal
#  - Asymptotic Calculator is used 
def test_limit_3() :
    """ Get the upper limit at given point for small signal
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_3")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )

    the_model = model.clone  ( S = NS , name = 'X3' )
    data      = dataset
    
    reff.setVal ( 1 ) 
    reff.fix() 
    constraints = ceff1 , 

    logger.info ('Fit S+B') 
    with use_canvas ( 'test_limit_3' ) :
        rS , _ = the_model.fitTo ( data , draw = True , silent = True )

    POI = SS 

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf         = the_model   ,
                             poi         = POI         , ## parameter of interest
                             dataset     = data        ,
                             constraints = constraints ,                             
                             name        = 'S+B'       ,
                             snapshot    = POI         ) ## ATTENTION! 

    with FIXVAR ( POI ) :
        POI.setVal ( 0 )
        rB , _ = the_model.fitTo ( data , silent = True , constraints = constraints )

        # create ModelConfig  for 'B-only' model
        model_b  = ModelConfig ( pdf         = the_model          ,
                                 poi         = POI                , ## parameter of interest                                 
                                 dataset     = data               ,
                                 constraints = constraints        ,
                                 workspace   = model_sb.workspace ,                                 
                                 name        = 'B-only'           ,
                                 snapshot    = POI                )
      
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
        hti .scan_with_progress ( values ) #can it!
          
    ## visualize the scan results 
    with use_canvas ( 'test_point_limit_ac3: Fixed&constrained efficiency' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.2f' % hti.upper_limit )

    row = 'Fixed&constrained' , '%.3f' % hti.upper_limit , '%.1f' % timer.delta 
    summary.append ( row  )
    plots  .append ( plot )


# ============================================================================-
## Get the upper limit limit for small signal
#  - Asymptotic Calculator is used 
def test_limit_4() :
    """ Get the upper limit at given point for small signal
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_point_4")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )


    the_model = model.clone  ( S = NS , name = 'X4' )
    data      = dataset
    
    reff.setVal ( 1 ) 
    reff.release() 
    constraints = ceff1 , 

    logger.info ('Fit S+B') 
    with use_canvas ( 'test_limit_3' ) :
        rS , _ = the_model.fitTo ( data , draw = True , silent = True )

    POI = SS 

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf         = the_model   ,
                             poi         = POI         , ## parameter of interest
                             dataset     = data        ,
                             constraints = constraints , 
                             name        = 'S+B'       ,
                             snapshot    = POI         ) ## ATTENTION! 
    

    with FIXVAR ( POI ) :
        POI.setVal ( 0 )
        rB , _ = the_model.fitTo ( data , silent  = True , constraints = constraints )

        # create ModelConfig  for 'B-only' model
        model_b  = ModelConfig ( pdf         = the_model          ,
                                 poi         = POI                , ## parameter of interest                                 
                                 dataset     = data               ,
                                 constraints = constraints        ,                                  
                                 workspace   = model_sb.workspace , 
                                 name        = 'B-only'           ,
                                 snapshot    = POI                )
      
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
        hti .scan_with_progress ( values ) #can it!
          
    ## visualize the scan results 
    with use_canvas ( 'test_point_limit_ac4: Constrained efficinecy#1' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.2f' % hti.upper_limit )

    row = 'Constrained efficiency/1' , '%.3f' % hti.upper_limit , '%.1f' % timer.delta 
    summary.append ( row  )
    plots  .append ( plot )

# ============================================================================-
## Get the upper limit limit for small signal
#  - Asymptotic Calculator is used 
def test_limit_5() :
    """ Get the upper limit at given point for small signal
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_limit_5")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )


    the_model = model.clone  ( S = NS , name = 'X5' )
    data      = dataset
    
    reff.setVal ( 1 ) 
    reff.release() 
    constraints = ceff2 , 
    gobs = ROOT.RooArgSet ( reff_mean , reff_sigma ) 

    logger.info ('Fit S+B') 
    with use_canvas ( 'test_limit_3' ) :
        rS , _ = the_model.fitTo ( data , draw = True , silent = True )

    POI = SS 

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf                = the_model   ,
                             poi                = POI         , ## parameter of interest
                             dataset            = data        ,
                             global_observables = gobs        , 
                             constraints        = constraints , 
                             name               = 'S+B'       ,
                             snapshot           = POI         ) ## ATTENTION! 


    with FIXVAR ( POI ) :
        POI.setVal ( 0 )
        rB , _ = the_model.fitTo ( data , silent  = True , constraints = constraints )

        # create ModelConfig  for 'B-only' model
        model_b  = ModelConfig ( pdf         = the_model          ,
                                 poi         = POI                , ## parameter of interest                                 
                                 dataset     = data               ,
                                 constraints        = constraints ,
                                 global_observables = gobs              , 
                                 workspace          = model_sb.workspace , 
                                 name               = 'B-only'           ,
                                 snapshot           = POI                )
      
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
        hti .scan_with_progress ( values ) #can it!
          
    ## visualize the scan results 
    with use_canvas ( 'test_point_limit_ac5: Constrained efficiency#2' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.2f' % hti.upper_limit )

    row = 'Constrained efficiency/2' , '%.3f' % hti.upper_limit , '%.1f' % timer.delta 
    summary.append ( row  )
    plots  .append ( plot )

    
# =============================================================================
if '__main__' == __name__ :

    from ostap.core.core    import rooSilent
    
    with rooSilent ( ) : 

        test_limit_1 () ## trivial  
        test_limit_2 () ## add efficiency 
        test_limit_3 () ## add efficiency 
        test_limit_4 () ## add efficiency 
        test_limit_5 () ## add efficiency 
 
        
    import ostap.logger.table as T
    title = 'Summary of 90%CL Upper Limits'
    table = T.table ( summary , title = title , prefix = '# ' , alignment = 'lcr' )
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
