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
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.meta_info     import root_info 
from   ostap.fitting.variables  import SETVAR, FIXVAR  
from   ostap.core.core          import cpp, VE, rooSilent, Ostap 
from   ostap.utils.timing       import timing
from   ostap.utils.utils        import vrange, lrange
from   ostap.plotting.canvas    import use_canvas
from   ostap.fitting.simfit     import combined_data
from   ostap.utils.progress_bar import progress_bar
from   ostap.utils.utils        import batch_env 
import ostap.fitting.models     as     Models
import ostap.logger.table       as     T
import ostap.fitting.roofit 
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
## set batch from environment 
batch_env ( logger )
# =============================================================================

# =============================================================================
## more complicated/sophisticated  stuff
# =============================================================================

mass    = ROOT.RooRealVar   ('mass','mass-variable', 0 , 10 )
vars    = ROOT.RooArgSet    ( mass ) 
signal1 = Models.Gauss_pdf  ( 'Gauss1',
                              xvar  = mass                 ,
                              mean  = ( 5.0 , 1.5  , 8.5 ) ,
                              sigma = ( 0.5 , 0.01 , 3.0 ) )
signal2 = Models.Gauss_pdf  ( 'Gauss2',
                              xvar  = mass                 ,
                              mean  = signal1.mean         ,
                              sigma = signal1.sigma        )         
           

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
    
ds_A.name    = 'ds_A'
ds_B.name    = 'ds_B'
dataset.name = 'ds_combined'

model_sim.fitTo ( dataset , quiet = True ) 


err        = 0.3
reff       = ROOT.RooRealVar ( 'reff' , 'ratio of efficiencies' , 1.0 , 0 , 10 )
ceff1      = model2.soft_constraint ( reff , VE ( 1.0, err**2 ) )


reff_mean  = ROOT.RooRealVar  ( 'reff_mean'  , 'mean  reff' , 1.0  , 0     , 10 )
reff_mean .fix()

reff_sigma = ROOT.RooRealVar  ( 'reff_sigma' , 'sigma reff' , err , 1.e-6 , 0.5 )

reff_sigma.fix()

ceff2      = ROOT.RooGaussian ( 'CG2' , 'CG2' , reff_mean , reff , reff_sigma )

    
summary = [ ( 'method' , '90%CL' , 'time [s]') ]
plots   = [] 

values  = lrange ( 1.e-5 , 0.08 , 50 )
# ============================================================================-
## Get the upper limit limit for small signal
#  - Asymptotic Calculator is used 
def test_point_limit_1() :
    """ Get the upper limit at given point for small signal
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_point_limit_1")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator, No efficiency)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )

    the_model = model_sim.clone ( name = 'X1' )
    data      = dataset 

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

    model_sb.ws.Import ( data )
    model_sb.ws.Import ( ds_A )
    model_sb.ws.Import ( ds_B )

    logger.info  ( 'ModelConfig:\n%s' % ( model_b.mc.table ( prefix = '# ' ) ) ) 
    logger.info  ( 'Workspace:\n%s'   % ( model_b.ws.table ( prefix = '# ' ) ) ) 

    print ( model_sb.table()                 )
    print ( model_sb.table( style = 'local') )

    print ( model_sb.mc.table()                 )
    print ( model_sb.mc.table( style = 'local') )

    print ( model_sb.ws.table()                 )
    print ( model_sb.ws.table( style = 'local') )
    
    with timing ( "Using Asymptotic Calculator, No efficiency" , logger = logger ) as timer :
        ## create the calculator 
        ac  = AsymptoticCalculator ( model_b           ,
                                     model_sb          ,
                                     dataset   = data  ,
                                     asimov    = False )
        ac.calculator.SetOneSided ( True ) 
        
        ## create Hypo Test inverter 
        hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
        
        ## make a scan 
        hti .scan_with_progress ( values  ) # scan it!
          
    ## visualize the scan results 
    with use_canvas ( 'test_point_limit_1: No efficiency' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.4f' % hti.upper_limit )

    row = 'No efficiency' , '%.4f' % hti.upper_limit , '%.2f' % timer.delta 
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

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator, Fixed efficiency)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )


    S3          = signal2.vars_multiply ( SB , reff , name = 'S3' )
    model2_new  = model2.clone  ( S = S3 , name = 'M3' )

    reff.fix()
    
    the_model   = Models.SimFit (  sample  , { 'A' : model1  , 'B' : model2_new } , name = 'X3' )

    data        = dataset 
    
    fit_conf = { 'silent' : True }
    
    rS , _ = the_model.fitTo ( data , **fit_conf ) 
    ## with use_canvas ( 'test_point_limit_2/A' ) : the_model.draw  ( 'A' , data )  
    ## with use_canvas ( 'test_point_limit_2/B' ) : the_model.draw  ( 'B' , data )  
                   
    POI = fBA 

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf         = the_model   ,
                             poi         = POI         , ## parameter of interest
                             dataset     = data        ,
                             name        = 'S+B'       ,
                             snapshot    = POI         ) ## ATTENTION! 
    

    with FIXVAR ( POI ) :
        POI.setVal ( 0 )
        rB , _ = the_model.fitTo ( data , **fit_conf )

        # create ModelConfig  for 'B-only' model
        model_b  = ModelConfig ( pdf         = the_model          ,
                                 poi         = POI                , ## parameter of interest
                                 dataset     = data               ,
                                 workspace   = model_sb.workspace , 
                                 name        = 'B-only'           ,
                                 snapshot    = POI                )
      
    logger.info ( 'Model config %s\n%s'  % ( model_sb.name , model_sb.table ( prefix = '# ' ) ) ) 
    logger.info ( 'Model config %s\n%s'  % ( model_b.name  , model_b .table ( prefix = '# ' ) ) )
    
    with timing ( "Using Asymptotic Calculator,Fixed efficiency" , logger = logger ) as timer :
        ## create the calculator 
        ac  = AsymptoticCalculator ( model_b           ,
                                     model_sb          ,
                                     dataset   = data  ,
                                     asimov    = False )
        ac.calculator.SetOneSided ( True ) 
        
        ## create Hypo Test inverter 
        hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
        
        ## make a scan 
        hti .scan_with_progress ( values ) # scan it!
          
    ## visualize the scan results 
    with use_canvas ( 'test_point_limit_ac2: Fixed efficiency' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.4f' % hti.upper_limit )

    row = 'Fixed efficiency' , '%.4f' % hti.upper_limit , '%.2f' % timer.delta 
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

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator,Fixed&Constrained efficiency)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )


    S3          = signal2.vars_multiply ( SB , reff , name = 'S4' )
    model2_new  = model2.clone  ( S = S3 , name = 'M4' )

    reff.fix()
    constraints = ceff1

    the_model   = Models.SimFit (  sample  , { 'A' : model1  , 'B' : model2_new } , name = 'X4' )

    data        = dataset 
    
    fit_conf = { 'silent' : True , 'constraints' : constraints }

    rS , _ = the_model.fitTo ( data , **fit_conf ) 
    ## with use_canvas ( 'test_point_limit_3/A' ) : the_model.draw  ( 'A' , data )  
    ## with use_canvas ( 'test_point_limit_3/B' ) : the_model.draw  ( 'B' , data )  
                   
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
        rB , _ = the_model.fitTo ( data , **fit_conf )

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
    
    with timing ( "Using Asymptotic Calculator,Fixed&Constrained efficiency" , logger = logger ) as timer :
        ## create the calculator 
        ac  = AsymptoticCalculator ( model_b           ,
                                     model_sb          ,
                                     dataset   = data  ,
                                     asimov    = False )
        ac.calculator.SetOneSided ( True ) 
        
        ## create Hypo Test inverter 
        hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
        
        ## make a scan 
        hti .scan_with_progress ( values ) # scan it!
          
    ## visualize the scan results 
    with use_canvas ( 'test_point_limit_3: Fixed&Constrained efficiency' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.4f' % hti.upper_limit )

    row = 'Fixed&Constrained efficiency' , '%.4f' % hti.upper_limit , '%.2f' % timer.delta 
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

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator, Constrained efficiency)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )

    S3          = signal2.vars_multiply ( SB , reff , name = 'S5' )
    model2_new  = model2.clone  ( S = S3 , name = 'M5' )

    reff.release()
    constraints = ceff1
    
    the_model   = Models.SimFit (  sample  , { 'A' : model1  , 'B' : model2_new } , name = 'X5' )

    data        = dataset 
    
    fit_conf = { 'silent' : True , 'constraints' : constraints }

    rS , _ = the_model.fitTo ( data , silent  = True , constraints = constraints ) 
    ## with use_canvas ( 'test_point_limit_4/A' ) : the_model.draw  ( 'A' , data )  
    ## with use_canvas ( 'test_point_limit_4/B' ) : the_model.draw  ( 'B' , data )  
                   
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
        rB , _ = the_model.fitTo ( data , **fit_conf )

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
        hti .scan_with_progress ( values ) # scan it!
          
    ## visualize the scan results 
    with use_canvas ( 'test_point_limit_4: Constrained efficinecy#1' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.4f' % hti.upper_limit )

    row = 'Constrained efficiency/1' , '%.4f' % hti.upper_limit , '%.2f' % timer.delta 
    summary.append ( row  )
    plots  .append ( plot )

        

# ============================================================================-
## Get the upper limit limit for small signal
#  - Asymptotic Calculator is used
#  - Varinat with Global observables 
def test_point_limit_5() :
    """ Get the upper limit at given point for small signal
    - Asymptotic Calculator is used 
    - Varinat with Global observables 
    """

    logger = getLogger("test_point_limit_5")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator, Consrained efficienicy&GlobalObservables)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )


    S3          = signal2.vars_multiply ( SB , reff , name = 'S6' )
    model2_new  = model2.clone  ( S = S3 , name = 'M6' )

    reff.release()
    constraints = ceff2

    global_observables  = reff_mean 

    the_model   = Models.SimFit (  sample  , { 'A' : model1  , 'B' : model2_new } , name = 'X6' )

    data        = dataset 
    
    fit_conf = { 'silent' : True , 'constraints' : constraints }
    
    rS , _ = the_model.fitTo ( data , **fit_conf ) 
    ## with use_canvas ( 'test_point_limit_5/A' ) : the_model.draw  ( 'A' , data )  
    ## with use_canvas ( 'test_point_limit_5/B' ) : the_model.draw  ( 'B' , data )  
                   
    POI = fBA 

    ## create ModelConfig  for 'S+B' model
    model_sb = ModelConfig ( pdf                = the_model          ,
                             poi                = POI                , ## parameter of interest
                             dataset            = data               ,
                             global_observables = global_observables , 
                             constraints        = constraints        , 
                             name               = 'S+B'              ,
                             snapshot           = POI                ) ## ATTENTION! 

    with FIXVAR ( POI ) :
        POI.setVal ( 0 )
        rB , _ = the_model.fitTo ( data , **fit_conf )

        # create ModelConfig  for 'B-only' model
        model_b  = ModelConfig ( pdf                = the_model          ,
                                 poi                = POI                , ## parameter of interest                                 
                                 dataset            = data               ,
                                 constraints        = constraints        ,
                                 global_observables = global_observables , 
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
        hti .scan_with_progress ( values ) # scan it!
          
    ## visualize the scan results 
    with use_canvas ( 'test_point_limit_5: Constrained efficiency#2' , wait = 2 ) :
        plot = hti.plot
        plot .draw('LCb 2CL')                    
        logger.info ( '90%%CL upper limit (asymptotic)  = %.4f' % hti.upper_limit )

    row = 'Constrained efficiency/2' , '%.4f' % hti.upper_limit , '%.2f' % timer.delta 
    summary.append ( row  )
    plots  .append ( plot )


# ============================================================================-
## Get the upper limit limit for small signal
#  - Asymptotic Calculator is used 
def test_point_limit_6() :
    """ Get the upper limit at given point for small signal
    - Asymptotic Calculator is used 
    """

    logger = getLogger("test_point_limit_6")

    logger.info ( "Test Point limits with RooStats (Asymptotic Calculator, Constrained efficiencies)" )

    from   ostap.fitting.roostats   import ( ModelConfig           ,
                                             AsymptoticCalculator  ,
                                             HypoTestInverter      )


    S3          = signal2.vars_multiply ( SB , reff , name = 'S7' )
    model2_new  = model2.clone  ( S = S3 , name = 'M7' )

    reff.release()
    
    the_model   = Models.SimFit (  sample  , { 'A' : model1  , 'B' : model2_new } , name = 'X5' )

    data        = dataset 

    for eff in progress_bar ( ( VE ( 1.0 , 0.01 **2 ) , 
                                VE ( 1.0 , 0.05 **2 ) ,
                                VE ( 1.0 , 0.10 **2 ) ,
                                VE ( 1.0 , 0.15 **2 ) ,
                                VE ( 1.0 , 0.20 **2 ) ,
                                VE ( 1.0 , 0.25 **2 ) ,
                                VE ( 1.0 , 0.30 **2 ) ,
                                VE ( 1.0 , 0.35 **2 ) ,
                                VE ( 1.0 , 0.40 **2 ) ,
                                VE ( 1.0 , 0.45 **2 ) ,
                                VE ( 1.0 , 0.50 **2 ) ) ) :
        
        ceff = model2.soft_constraint ( reff , eff )
        constraints = ceff 
                
        fit_conf = { 'silent' : True , 'constraints' : constraints }
        
        rS , _ = the_model.fitTo ( data , **fit_conf ) 
        
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
            rB , _ = the_model.fitTo ( data , **fit_conf )
            
            # create ModelConfig  for 'B-only' model
            model_b  = ModelConfig ( pdf         = the_model          ,
                                     poi         = POI                , ## parameter of interest                                 
                                     dataset     = data               ,
                                     constraints = constraints        ,                                  
                                     workspace   = model_sb.workspace , 
                                     name        = 'B-only'           ,
                                     snapshot    = POI                )
            
        with timing ( "Using Asymptotic Calculator" , logger = None ) as timer :
            ## create the calculator 
            ac  = AsymptoticCalculator ( model_b           ,
                                         model_sb          ,
                                         dataset   = data  ,
                                         asimov    = False )
            ac.calculator.SetOneSided ( True ) 
            
            ## create Hypo Test inverter 
            hti = HypoTestInverter ( ac ,  0.90 , use_CLs = True , verbose = False )
            
            ## make a scan 
            hti .scan ( values ) # scan it!
            
        row = 'Constrained efficiency %s' % eff.toString ( '%+.3f +/- %.3f' )  , '%.4f' % hti.upper_limit , '%.2f' % timer.delta 
        summary.append ( row  )
            
# =============================================================================
if '__main__' == __name__ :

    from ostap.core.core    import rooSilent
    
    with rooSilent ( ) : 

        test_point_limit_1 () ## trivial  
        # test_point_limit_2 () ## add efficiency 
        # test_point_limit_3 () ## add efficiency 
        # test_point_limit_4 () ## add efficiency 
        # test_point_limit_5 () ## add efficiency 
        # test_point_limit_6 () ## scan efficiencies 
 
        
    import ostap.logger.table as T
    title = 'Summary of 90%CL Upper Limits'
    table = T.table ( summary , title = title , prefix = '# ' , alignment = 'lcr' )
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
