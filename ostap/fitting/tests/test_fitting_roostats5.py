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
from   ostap.utils.utils        import vrange, lrange
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait
from   ostap.fitting.simfit     import combined_data
import ostap.logger.table       as     T
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_roostats5' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

# =============================================================================
## more complicated/sophisticated  stuff
# =============================================================================

mass    = ROOT.RooRealVar   ('mass','mass-variable', 0 , 10 )
vars    = ROOT.RooArgSet    ( mass ) 
signal1 = Models.Gauss_pdf  ( 'Gauss1',
                              xvar  = mass                 ,
                              ## mean  = ( 5.0 , 1.5  , 8.5 ) ,
                              ## sigma = ( 0.5 , 0.01 , 3.0 ) )
                              mean  = ROOT.RooFit.RooConst ( 5 ) ,
                              sigma = ROOT.RooFit.RooConst ( 1 ) ) 
signal2 = Models.Gauss_pdf  ( 'Gauss2',
                              xvar  = mass                 ,
                              mean  = signal1.mean         ,
                              sigma = signal1.sigma        )         
           
## signal1.sigma.fix()
## signal1.mean .fix()

NS_A     = 200
NB_A     = 50

model1   = Models.Fit1D ( signal = signal1 , background = 'p1' , name = 'A' , suffix = '_A')
model1.S = NS_A 
model1.B = NB_A

fBA      = ROOT.RooRealVar ( 'fBA'  , 'SB/SA-ratio', 0.0 , -0.10 , 1.0 )
SB       = signal2.vars_multiply   ( model1.S , fBA , name  = 'SB' )
model2   = Models.Fit1D            ( signal = signal2 , background = 'p1' , S = SB , name = 'B' , suffix = '_B' )

NS_B     = float ( SB )  
NB_B     = 5

model2.B = NB_B 

NA   = NS_A + NB_A  
NB   = NS_B + NB_B
NTOT = NA   + NB 
 
## combine PDFs
sample     = ROOT.RooCategory ( 'sample' ,'sample'  , 'A' , 'B' )
allvars    = ROOT.RooArgSet   (  mass    , sample ) 
model_sim  = Models.SimFit    (  sample  , { 'A' : model1  , 'B' : model2 } , name = 'X' )

with FIXVAR ( fBA  ) :
    fBA.setVal ( 0 ) 
    ds_A   = model1           .generate ( NA , sample = False )
    print ( ds_A ) 
    ds_B   = model2.background.generate ( NB , sample = False )
    print ( ds_B ) 
    dataset = combined_data ( sample , allvars,  { 'A' : ds_A , 'B'  : ds_B } )
    

model_sim.fitTo ( dataset , quiet = True ) 


reff       = ROOT.RooRealVar ( 'reff' , 'ratio of efficiencies' , 1.0 , 0 , 10 )
ceff1      = model2.soft_constraint ( reff , VE ( 1.0, 0.001**2 ) )


reff_mean  = ROOT.RooRealVar  ( 'reff_mean'  , 'mean  reff' , 1.0  , 0     , 10 )
reff_mean .fix()

reff_sigma = ROOT.RooRealVar  ( 'reff_sigma' , 'sigma reff' , 0.001 , 1.e-6 , 0.5 )
## reff_sigma = ROOT.RooFit.RooConst (  0.001 ) 

reff_sigma.fix()

ceff2      = ROOT.RooGaussian ( 'CG2' , 'CG2' , reff_mean , reff , reff_sigma )

    
summary = [ ( 'method' , '90%CL' , 'time [s]') ]
plots   = [] 

values  = lrange ( 1.e-5 , 0.05 , 30 )
# ============================================================================-
## Get the upper limit limit for small signal
#  - Asymptotic Calculator is used 
def test_point_limit_1() :
    """ Get the upper limit at given point for small signal
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_point_limit_1")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )

    the_model = model_sim.clone ( name = 'X1' )
    data      = dataset 

    logger.info ('Fit S+B') 
    rS , _ = the_model.fitTo ( data , quiet = True ) 
    with use_canvas ( 'test_point_limit_1/A' ) : the_model.draw  ( 'A' , data )  
    with use_canvas ( 'test_point_limit_1/B' ) : the_model.draw  ( 'B' , data )  
                   
    POI = fBA 

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf         = the_model   ,
                             poi         = POI         , ## parameter of interest 
                             dataset     = data        ,
                             name        = 'S+B'       ,
                             snapshot    = POI         ) ## ATTENTION! 
    
    with FIXVAR ( POI ) :
        POI.setVal ( 0 )
        logger.info ( 'Fit B-only' )        
        rB , _ = the_model.fitTo ( data , quiet = True )

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
def test_point_limit_2() :
    """ Get the upper limit at given point for small signal
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_point_limit_2")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )


    S3          = signal2.vars_multiply ( SB , reff , name = 'S3' )
    model2_new  = model2.clone  ( S = S3 , name = 'M3' )

    reff.fix()
    
    the_model   = Models.SimFit (  sample  , { 'A' : model1  , 'B' : model2_new } , name = 'X3' )

    data        = dataset 
    
    rS , _ = the_model.fitTo ( data , silent = True ) 
    with use_canvas ( 'test_point_limit_2/A' ) : the_model.draw  ( 'A' , data )  
    with use_canvas ( 'test_point_limit_2/B' ) : the_model.draw  ( 'B' , data )  
                   
    POI = fBA 

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
def test_point_limit_3() :
    """ Get the upper limit at given point for small signal
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_point_limit_3")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )


    S3          = signal2.vars_multiply ( SB , reff , name = 'S4' )
    model2_new  = model2.clone  ( S = S3 , name = 'M4' )

    reff.fix()
    constraints = ceff1 , 

    the_model   = Models.SimFit (  sample  , { 'A' : model1  , 'B' : model2_new } , name = 'X4' )

    data        = dataset 
    
    rS , _ = the_model.fitTo ( data , silent = True , constraints = constraints ) 
    with use_canvas ( 'test_point_limit_3/A' ) : the_model.draw  ( 'A' , data )  
    with use_canvas ( 'test_point_limit_3/B' ) : the_model.draw  ( 'B' , data )  
                   
    POI = fBA 

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
def test_point_limit_4() :
    """ Get the upper limit at given point for small signal
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_point_limit_4")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )


    S3          = signal2.vars_multiply ( SB , reff , name = 'S5' )
    model2_new  = model2.clone  ( S = S3 , name = 'M5' )

    reff.release()
    constraints = ceff1 , 
    
    the_model   = Models.SimFit (  sample  , { 'A' : model1  , 'B' : model2_new } , name = 'X5' )

    data        = dataset 
    
    rS , _ = the_model.fitTo ( data , silent  = True , constraints = constraints ) 
    with use_canvas ( 'test_point_limit_4/A' ) : the_model.draw  ( 'A' , data )  
    with use_canvas ( 'test_point_limit_4/B' ) : the_model.draw  ( 'B' , data )  
                   
    POI = fBA 

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf         = the_model   ,
                             poi         = POI         , ## parameter of interest
                             dataset     = data        ,
                             constraints = constraints , 
                             name        = 'S+B'       ,
                             snapshot    = POI         ) ## ATTENTION! 
    

    with FIXVAR ( POI ) :
        POI.setVal ( 0 )
        rB , _ = the_model.fitTo ( data , silent  = True , constraints = True )

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
def test_point_limit_5() :
    """ Get the upper limit at given point for small signal
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_point_limit_5")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )


    S3          = signal2.vars_multiply ( SB , reff , name = 'S6' )
    model2_new  = model2.clone  ( S = S3 , name = 'M6' )

    reff.release()
    constraints = ceff2 , 

    gobs = ROOT.RooArgSet ( reff_mean , reff_sigma ) 

    the_model   = Models.SimFit (  sample  , { 'A' : model1  , 'B' : model2_new } , name = 'X6' )

    data        = dataset 
    
    rS , _ = the_model.fitTo ( data , silent = True , constraints = constraints ) 
    with use_canvas ( 'test_point_limit_5/A' ) : the_model.draw  ( 'A' , data )  
    with use_canvas ( 'test_point_limit_5/B' ) : the_model.draw  ( 'B' , data )  
                   
    POI = fBA 

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

        test_point_limit_1 () ## trivial  
        test_point_limit_2 () ## add efficiency 
        test_point_limit_3 () ## add efficiency 
        test_point_limit_4 () ## add efficiency 
        test_point_limit_5 () ## add efficiency 
 
        
    import ostap.logger.table as T
    title = 'Summary of 90%CL Upper Limits'
    table = T.table ( summary , title = title , prefix = '# ' , alignment = 'lcr' )
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
