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

Simultannepous fit of two 1D-distributions

"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ostap.io.zipshelve       as     DBASE
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   builtins                 import range 
from   ostap.core.core          import dsID, rooSilent
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
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
## make simple test mass 
mass     = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 0 , 5 )

## book very simple data set:
varset1  = ROOT.RooArgSet  ( mass )
dataset1 = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset1 )  

## book very simple data set:
varset2  = ROOT.RooArgSet  ( mass )
dataset2 = ROOT.RooDataSet ( dsID() , 'Test Data set-2' , varset2 )  

## high statistic, low-background "control channel"
mean1  = 2.0
sigma1 = 0.50
NS1    =  10000
NB1    =   1000

for i in range ( NS1 )  :
    v1 = random.gauss( mean1 , sigma1 )
    if v1 in mass :
        mass.setVal ( v1 )
        dataset1.add ( varset1 )

for i in range ( NB1 ) :
    v1 = random.uniform ( 0 , 5  )
    if v1 in mass :
        mass.setVal ( v1 )
        dataset1.add ( varset1 )
        
## low statistic, high-background "control channel"
NS2     =    500
NB2     =  10000     
mean2   = mean1  + 1.0
sigma2  = sigma1 * 0.5

for i in range(NS2)  :
    v2 = random.gauss( mean2 , sigma2 )
    if v2 in mass :
        mass.setVal ( v2 )
        dataset2.add ( varset2 )

for i in range (NB2 ) :
    v2 = random.uniform ( 0 , 5  )
    if v2 in mass :
        mass.setVal ( v2 )
        dataset2.add ( varset2 )

sample  = ROOT.RooCategory ('sample','sample'  , 'A' , 'B' )

models  = set()
results = [] 
# =============================================================================
def test_simfit1 () :
## if 1 < 2 :
    
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
    
    with use_canvas ( 'test_simfit1' ) : 
        # =========================================================================
        ## fit 1
        with wait ( 1 ) : 
            r1 , f1 = model1.fitTo ( dataset1 , draw = True , nbins = 50 , silent = True )
            title = 'Results of fit to dataset1'
            logger.info ( '%s\n%s' % ( title , r1.table ( title = title , prefix = '# ' ) ) )
        
        ## fit 2
        with wait ( 1 ) : 
            r2 , f2 = model2.fitTo ( dataset2 , draw = True , nbins = 50 , silent = True )
            title = 'Results of fit to dataset2'
            logger.info ( '%s\n%s' % ( title , r2.table ( title = title , prefix = '# ' ) ) )
        # =========================================================================

    ## combine data
    
    
    ## combine datasets
    from ostap.fitting.simfit import combined_data 
    vars    = ROOT.RooArgSet ( mass )
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
    
    with use_canvas ( 'test_simfit1' ) :
        with wait ( 1 ) : 
            fA = model_sim.draw ( 'A' , dataset , nbins = 50 )
        with wait ( 1 ) : 
            fB = model_sim.draw ( 'B' , dataset , nbins = 50 )            

    models.add ( model1        )
    models.add ( model2        )
    models.add ( model_sim     )
    models.add ( model_sim.pdf )
    
    for k in model_sim.categories : models.add ( model_sim.categories[k] )
    for k in model_sim.drawpdfs   : models.add ( model_sim.drawpdfs  [k] )
        
    results.append ( r1 ) 
    results.append ( r2 ) 
    results.append ( r  ) 
    
    
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
            db['model:'     + m.name ] = m
            db['roo_tot:%s' % m.name ] = m.pdf
        db['models'  ] = models
        for r in results : db['result' + r.name ] = r 
        db['results' ] = results
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
