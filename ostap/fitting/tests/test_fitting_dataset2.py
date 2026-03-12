#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_dataset.py
# - It tests some decorators fore dataset
# ============================================================================= 
"""  It tests some jacknife/bootstrapping for datatsets 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.core          import dsID, hID, Ostap 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.progress_bar import progress_bar
from   ostap.utils.memory       import memory 
from   ostap.utils.root_utils   import batch_env 
import ostap.logger.table       as     T
import ostap.fitting.roofit 
import ostap.trees.trees   
import ostap.histos.histos 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_dataset2' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch form environment 
batch_env ( logger )
# =============================================================================


evt     = ROOT.RooRealVar ( 'Evt'    , '#event'        , 0 , 1000000 )
run     = ROOT.RooRealVar ( 'Run'    , '#run'          , 0 , 1000000 )
mass    = ROOT.RooRealVar ( 'Mass'   , 'mass-variable' , 0 , 100     )
pt1     = ROOT.RooRealVar ( 'Pt1'    , 'pt1-variable'  , 0 )
pt2     = ROOT.RooRealVar ( 'Pt2'    , 'pt2-variable'  , 0 )
weight  = ROOT.RooRealVar ( 'Weight' , 'some weight'   , -10 , 10 )

varset  = ROOT.RooArgSet  ( evt , run , mass , pt1 , pt2 , weight )
dataset = ROOT.RooDataSet ( dsID () , 'Test Data set-0' , varset )  

NR  = 100
NE  = 100

with memory ( 'Create initial dataset' , logger = logger ) as dm0 :
    
    for r in progress_bar ( range ( NR ) )  :
        
        run.setVal ( r )
        for e in range ( NE ) :
            
            evt .setVal   ( e )
            mass.setVal   ( random.uniform      ( 0  , 10  ) )
            weight.setVal ( random.gauss        ( 1  , 0.1 ) )
            pt1.setVal     ( random.expovariate ( 1.0/25   ) ) 
            pt2.setVal     ( random.expovariate ( 1.0/15   ) ) 
            dataset.add ( varset )
            
            if ( 1 <= r <= 2  ) and 3 <=  e <= 3 :
                
                mass.setVal ( random.uniform ( 0 , 10 ) )
                dataset.add ( varset )

# =============================================================================
## (1) print datasets
# =============================================================================
logger.info ( 'Print initial dataset:\n%s' % dataset .table ( prefix = '# ' ) )



from ostap.math.base import std
data_ptr = std.unique_ptr(ROOT.RooAbsData)

# =============================================================================
## (2)  Bootstrapping 
# =============================================================================
with memory ( 'Bootstrapping' , logger  = logger ) as dm1 :

    size = 200
    for i , ds in progress_bar ( enumerate ( dataset.bootstrap ( size  , extended = True ) ) ,
                                 max_value   = size         , 
                                 description = "Bootstrap:" ) :
        
        ## ds = Ostap.MoreRooFit.delete_data ( ds )
        ## ds.Delete()
        ds = data_ptr ( ds )
        del ds 


# =============================================================================
## (3) Jackknife  
# =============================================================================
with memory ( 'Jackknife' , logger  = logger ) as dm2 :

    first , last = 0, 200 
    for i , ds in progress_bar ( enumerate ( dataset.jackknife ( first , last ) ) ,
                                 max_value   = 100          , 
                                 description = "Jackknife:" ) :

        ## ds = Ostap.MoreRooFit.delete_data ( ds )
        ## ds.Delete ()
        ds = data_ptr ( ds ) 
        del ds 

# ==============================================================================
## delete it! 
with memory ( 'Delete Initial dataset' , logger = logger ) as dm3 :
    
    ## dataset = Ostap.MoreRooFit.delete_data ( dataset )
    ## del dataset 
    dataset = data_ptr ( dataset )
    del dataset 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
