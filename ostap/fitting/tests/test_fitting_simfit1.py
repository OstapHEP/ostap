#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_simfit1.py
# Test module for ostap/fitting/simfit.py
# - It tests the most simple "Simultaneous fit"
# ============================================================================= 
""" Test module for ostap/fitting/simfit.py
- It tests the most simple ``Simultaneous fit'' :
Simultaneous fit of two 1D-distributions
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.core          import dsID, rooSilent
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait, batch_env  
import ostap.io.zipshelve       as     DBASE
import ostap.fitting.models     as     Models 
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_simfit1' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
#
## make simple test mass 
mass     = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 0 , 5  )
xyz      = ROOT.RooRealVar ( 'test_xyz'  , 'Some test xyz'  , 0 , 10 )

## book very simple data set:
varset1  = ROOT.RooArgSet  ( mass , xyz )
dataset1 = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset1 )  

## book very simple data set:
varset2  = ROOT.RooArgSet  ( mass , xyz )
dataset2 = ROOT.RooDataSet ( dsID() , 'Test Data set-2' , varset2 )  

## high statistic, low-background "control channel"
mean1  = 2.0
sigma1 = 0.50
NS1    =  10000
NB1    =   1000

for i in range ( NS1 )  :
    
    v1 = random.gauss ( mean1 , sigma1 )
    while not v1 in mass :
        v1 = random.gauss ( mean1 , sigma1 )
        
    vv = random.triangular( 0 , 10, 1 ) 
    while not vv in xyz : 
        vv = random.triangular ( 0, 10 , 1 ) 

    mass.setVal  ( v1 )
    xyz .setVal  ( vv )
    
    dataset1.add ( varset1 )

for i in range ( NB1 ) :
    v1 = random.uniform ( 0 , 5  )
    while not v1 in mass :
        v1 = random.uniform ( 0 , 5  )
        
    vv = random.triangular( 0 , 10, 9 ) 
    while not vv in xyz : 
        vv = random.triangular ( 0, 10 , 9 ) 

    mass.setVal  ( v1 )
    xyz .setVal  ( vv )
    
    dataset1.add ( varset1 )
        
## low statistic, high-background "signal channel"
NS2     =    500
NB2     =  10000     
mean2   = mean1  + 1.0
sigma2  = sigma1 * 0.5

for i in range ( NS2 )  :
    v2 = random.gauss ( mean2 , sigma2 )
    while not v2 in mass : 
        v2 = random.gauss ( mean2 , sigma2 )
        
    vv = random.expovariate ( 1/2.0 )
    while not vv in xyz : 
        vv = random.expovariate ( 1/2.0 )

    mass.setVal  ( v2 )
    xyz .setVal  ( vv )
    dataset2.add ( varset2 )

for i in range (NB2 ) :
    v2 = random.uniform ( 0 , 5  )
    while not v2 in mass : 
        v2 = random.uniform ( 0 , 5  )
        
    vv = random.expovariate ( 1/2.0 )
    while not vv in xyz : 
        vv = random.expovariate ( 1/2.0 )

    mass.setVal  ( v2      )
    xyz .setVal  ( 10 - vv )
    dataset2.add ( varset2 )


sample  = ROOT.RooCategory ('sample','sample'  , 'A' , 'B' )

models  = set()
results = []
graphs  = [] 
# =============================================================================
def test_simfit1 () :
    
    logger = getLogger( 'test_simfit1' )
    
    signal1  = Models.Gauss_pdf ( 'G1'                 ,
                                  xvar  = mass         ,
                                  mean  = (0.5 , 2.5 ) ,
                                  sigma = (0.1 , 1.0 ) )
    
    model1   = Models.Fit1D ( suffix = 'M1' , signal = signal1 ,  background = -1 )
    model1.S = NS1
    model1.B = NB1 
    
    
    mean2    = signal1.vars_add      ( signal1.mean  , 1.0  )
    sigma2   = signal1.vars_multiply ( signal1.sigma , 0.5  )
    
    signal2  = Models.Gauss_pdf ( 'G2'            ,
                                  xvar  = mass    ,
                                  mean  = mean2   ,
                                  sigma = sigma2  )
    
    model2  = Models.Fit1D ( suffix = 'M2' , signal = signal2 ,  background = model1.background )
    model2.S = NS2
    model2.B = NB2 
    
    with use_canvas ( 'test_simfit1: fit dataset1' , wait = 2 ) : 
        # =========================================================================
        ## fit 1
        r1 , f1 = model1.fitTo ( dataset1 , draw = True , nbins = 50 , silent = True )
        title = 'Results of fit to dataset1'
        logger.info ( '%s\n%s' % ( title , r1.table ( title = title , prefix = '# ' ) ) )
        
    with use_canvas ( 'test_simfit1: fit dataset2' , wait = 2 ) : 
        ## fit 2
        r2 , f2 = model2.fitTo ( dataset2 , draw = True , nbins = 50 , silent = True )
        title = 'Results of fit to dataset2'
        logger.info ( '%s\n%s' % ( title , r2.table ( title = title , prefix = '# ' ) ) )
        # =========================================================================

    ## combine data
        
    ## combine datasets
    from ostap.fitting.simfit import combined_data 
    vars    = ROOT.RooArgSet ( mass , xyz )
    dataset = combined_data  ( sample , vars , { 'A' : dataset1 , 'B' : dataset2 } )
    
    ## combine PDFs
    model_sim  = Models.SimFit (
        sample , { 'A' : model1  , 'B' : model2 } , name = 'X'
        )
    
    # =========================================================================
    r , f = model_sim.fitTo ( dataset , silent = True )
    r , f = model_sim.fitTo ( dataset , silent = True )

    title = 'Results of simultaneous fit'
    logger.info ( '%s\n%s' % ( title , r.table ( title = title , prefix = '# ' ) ) )
    
    with use_canvas ( 'test_simfit1: fit both datasets & draw A' , wait = 2 ) :        
        fA = model_sim.draw ( 'A' , dataset , nbins = 50 )
        graphs.append ( fA )
    with use_canvas ( 'test_simfit1: fit both datasets & draw B' , wait = 2 ) :        
        fB = model_sim.draw ( 'B' , dataset , nbins = 50 )            
        graphs.append ( fB )
    with use_canvas ( 'test_simfit1: graph-NLL for S_M2' , wait = 2 ) :
        from ostap.utils.utils import vrange 
        grs = model_sim.graph_nll     ( 'S_M2' , vrange ( 0 , 1000 , 100 ) , dataset , draw = True )
        grs.draw('apl')
        graphs.append ( grs )
    with use_canvas ( 'test_simfit1: graph-profile for S_M2' , wait = 2 ) :
        from ostap.utils.utils import vrange 
        grs = model_sim.graph_profile ( 'S_M2' , vrange ( 0 , 1000 , 50 ) , dataset , draw = True )
        grs.draw('apl')
        graphs.append ( grs )

    with use_canvas ( 'test_simfit1: residual  for S_M2' , wait = 2 ) :
        _ , residual , _  = model_sim.draw ( 'A' , dataset , nbins = 50 , residual = 'P' )
        residual.draw()
        graphs.append ( residual )
        
    with use_canvas ( 'test_simfit1: pull for S_M2' , wait = 2 ) :
        _ , _ , pull = model_sim.draw ( 'A' , dataset , nbins = 50 , pull     = 'P' )
        pull.draw()
        graphs.append ( pull )

    models.add ( model1        )
    models.add ( model2        )
    models.add ( model_sim     )
    ## models.add ( model_sim.pdf )
    
    ## for k in model_sim.categories : models.add ( model_sim.categories[k] )
    ## for k in model_sim.drawpdfs   : models.add ( model_sim.drawpdfs  [k] )
        
    results.append ( r1 ) 
    results.append ( r2 ) 
    results.append ( r  ) 

    # ========================================================================
    logger.info ( 'Make sPlot-analysis')
    model_sim.sPlot ( dataset ) 

    with use_canvas ( 'test_simfit1: sPlot/xyz for A (signal)'     , wait = 1 ) :
        dataset.draw ( 'test_xyz' , '(sample==0)*S_M1_sw' )

    with use_canvas ( 'test_simfit1: sPlot/xyz for A (background)' , wait = 1 ) :
        dataset.draw ( 'test_xyz' , '(sample==0)*B_M1_sw' )

    with use_canvas ( 'test_simfit1: sPlot/xyz for B (signal)'     , wait = 1 ) :
        dataset.draw ( 'test_xyz' , '(sample==1)*S_M2_sw' )

    with use_canvas ( 'test_simfit1: sPlot/xyz for B (background)' , wait = 1 ) :
        dataset.draw ( 'test_xyz' , '(sample==1)*B_M2_sw' )


    # =========================================================================
    ## test creation of dataset
    # =========================================================================
    ds_gen = model_sim.generate ( nEvents = { 'A' : len ( dataset1 ) ,
                                              'B' : len ( dataset2 ) } ,
                                  varset  = vars  )

    rg , f = model_sim.fitTo ( ds_gen , silent = True )
    rg , f = model_sim.fitTo ( ds_gen , silent = True )
    
    title = 'Results of simultaneous fit to generated dataset'
    logger.info ( '%s\n%s' % ( title , rg.table ( title = title , prefix = '# ' ) ) )
   
# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :
    
    logger = getLogger ( 'test_db' ) 

    logger.info('Saving all objects into DBASE')
    with timing('Save everything to DBASE' , logger ), DBASE.tmpdb() as db : 
        db['mass'     ] = mass
        db['vars1'    ] = varset1
        db['vars2'    ] = varset2
        db['dataset1' ] = dataset1
        db['dataset2' ] = dataset2
        db['sample'   ] = sample 
        for m in models :
            logger.info ( 'save %s '% m )  
            db [ 'model:'     + m.name ] = m
        db['models'  ] = models
        for r in results : db['result' + r.name ] = r 
        db['results' ] = results
        db['graphs'  ] = graphs
        db.ls()

# =============================================================================
if '__main__' == __name__ :

    with timing( "simfit-1" ,   logger ) :  
       test_simfit1 ()
        
    ## check finally that everything is serializeable:
    with timing ('Save to DB:'     , logger ) :
       test_db ()          
 
# =============================================================================
##                                                                      The END 
# =============================================================================
