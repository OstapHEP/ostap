#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_gof_simfit.py
# - Goodness-of-fit machenery for "Simultaneous fit"
# ============================================================================= 
""" Test module for ostap/fitting/simfit.py
 - Goodness-of-fit machenery for "Simultaneous fit"
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# =============================================================================
from   ostap.core.meta_info     import root_info 
from   ostap.core.core          import dsID, rooSilent
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.root_utils   import batch_env  
from   ostap.fitting.simfit     import combined_data 
from   ostap.stats.gof_simfit   import GoFSimFit, GoFSimFitToys 
import ostap.io.zipshelve       as     DBASE
import ostap.fitting.models     as     Models 
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_stats_gof_simfit' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
#
## make simple test mass 
mass     = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 0 , 5  )

## book very simple data set:
varset1  = ROOT.RooArgSet  ( mass  )
dataset1 = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset1 )  

## book very simple data set:
varset2  = ROOT.RooArgSet  ( mass )
dataset2 = ROOT.RooDataSet ( dsID() , 'Test Data set-2' , varset2 )  

## high statistic, low-background "control channel"
mean1  = 2.0
sigma1 = 0.50
NS1    = 10000
NB1    = 1000

for i in range ( NS1 )  :
    
    v1 = random.gauss ( mean1 , sigma1 )
    while not v1 in mass :
        v1 = random.gauss ( mean1 , sigma1 )
        
    mass.setVal  ( v1 )
    
    dataset1.add ( varset1 )

for i in range ( NB1 ) :
    v1 = random.uniform ( 0 , 5  )
    while not v1 in mass :
        v1 = random.uniform ( 0 , 5  )
        

    mass.setVal  ( v1 )
    
    dataset1.add ( varset1 )

# ============================================================================
## low statistic, high-background "signal channel"
NS2     =    500
NB2     =  10000     
mean2   = mean1  + 1.0
sigma2  = sigma1 * 0.5

for i in range ( NS2 )  :
    v2 = random.gauss ( mean2 , sigma2 )
    while not v2 in mass : 
        v2 = random.gauss ( mean2 , sigma2 )
        
    mass.setVal  ( v2 )
    dataset2.add ( varset2 )

for i in range (NB2 ) :
    v2 = random.uniform ( 0 , 5  )
    while not v2 in mass : 
        v2 = random.uniform ( 0 , 5  )
        
    mass.setVal  ( v2      )
    dataset2.add ( varset2 )

sample  = ROOT.RooCategory ( 'sample' , 'sample'  , 'A' , 'B' )


## combined datasets
vars    = ROOT.RooArgSet ( mass )
dataset = combined_data  ( sample , vars , { 'A' : dataset1 , 'B' : dataset2 } )

models  = set()
results = []
graphs  = []

signal1  = Models.Gauss_pdf ( 'G1'                 ,
                              xvar  = mass         ,
                              mean  = (0.5 , 2.5 ) ,
                              sigma = (0.1 , 1.0 ) )

model1   = Models.Fit1D ( suffix = 'MA' , signal = signal1 ,  background = 'flat' )
model1.S = NS1
model1.B = NB1 

mean2    = signal1.vars_add      ( signal1.mean  , 1.0  )
sigma2   = signal1.vars_multiply ( signal1.sigma , 0.5  )

signal2  = Models.Gauss_pdf ( 'G2'            ,
                              xvar  = mass    ,
                              mean  = mean2   ,
                              sigma = sigma2  )

model2  = Models.Fit1D ( suffix = 'MB' , signal = signal2 ,  background = model1.background )
model2.S = NS2
model2.B = NB2 

## combine PDFs
model_sim  = Models.SimFit (
    sample , { 'A' : model1  , 'B' : model2 } , name = 'X'
)

# =============================================================================
def test_simfit1 () :
    
    logger = getLogger( 'test_gof_simfit1' )
    
    with use_canvas ( 'test_gof_simfit1: fit dataset1' , wait = 2 ) : 
        # =========================================================================
        ## fit 1
        r1 , f1 = model1.fitTo ( dataset1 , draw = True , nbins = 50 , silent = True )
        title = 'Results of fit to dataset1'
        logger.info ( '%s\n%s' % ( title , r1.table ( title = title , prefix = '# ' ) ) )
        
    with use_canvas ( 'test_gof_simfit1: fit dataset2' , wait = 2 ) : 
        ## fit 2
        r2 , f2 = model2.fitTo ( dataset2 , draw = True , nbins = 50 , silent = True )
        title = 'Results of fit to dataset2'
        logger.info ( '%s\n%s' % ( title , r2.table ( title = title , prefix = '# ' ) ) )
        # =========================================================================

    graphs.append ( f1 ) 
    graphs.append ( f2 ) 
    # =========================================================================
    r , f = model_sim.fitTo ( dataset , silent = True )
    r , f = model_sim.fitTo ( dataset , silent = True )

    title = 'Results of simultaneous fit'
    logger.info ( '%s\n%s' % ( title , r.table ( title = title , prefix = '# ' ) ) )
    
    with use_canvas ( 'test_gof_simfit1: fit both datasets & draw A' , wait = 2 ) :        
        fA = model_sim.draw ( 'A' , dataset , nbins = 50 )
    with use_canvas ( 'test_gof_simfit1: fit both datasets & draw B' , wait = 2 ) :        
        fB = model_sim.draw ( 'B' , dataset , nbins = 50 )
        
    models.add ( model1        )
    models.add ( model2        )
    models.add ( model_sim     )

    graphs.append  ( fA ) 
    graphs.append  ( fB ) 
        
    results.append ( r1 ) 
    results.append ( r2 ) 
    results.append ( r  ) 

    """
    ## GOF machinery
    
    gof = GoFSimFit ( model_sim      ,
                      dataset        ,
                      parameters = r )

    with use_canvas ( 'test_gof_simfit1: GoF-A' , wait = 2 ) : gof.draw ( 'A' )
    with use_canvas ( 'test_gof_simfit1: GoF-B' , wait = 2 ) : gof.draw ( 'B' )
        
    title = 'GoF for 1D SimFit' 
    logger.info ( '%s:\n%s' % ( title , gof.table ( title = title , prefix = '# ' ) ) )

    toys = GoFSimFitToys ( gof )
    toys.run ( 500 , silent = False , parallel = True )

    title = 'GoF for 1D SimFit %s toys)' % toys.nToys  
    logger.info ( '%s:\n%s' % ( title , toys.table ( title = title , prefix = '# ' ) ) )

    ## 
    for sample, g in toys.gofs.items () : 
        for k in g.estimators :
            with use_canvas ( 'test_gof_simfit1: GoF-%s %s' % ( sample , k ) , wait = 1 ) :
                toys .draw ( sample , k )
    """

    from   ostap.stats.gof_simfit   import PPDSimFit
    gof_ppd = PPDSimFit ( model_sim          , 
                          dataset            ,
                          parameters = r     ,                          
                          mcFactor   = 20    , 
                          nToys      = 1000  ,
                          sigma      = 0.5   ,
                          silent     = False )

    print ( 'PPD:', gof_ppd.tvalues() )
    print ( 'PPD:', gof_ppd.pvalues() )
    
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
            logger.info ( 'save %s\n%s '% ( m.name , m )  )  
            db [ 'model:'     + m.name ] = m
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
