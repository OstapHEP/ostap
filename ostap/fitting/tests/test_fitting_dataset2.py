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
from   ostap.logger.symbols     import delta_symbol, ram 
from   ostap.fitting.dataset import data_ptr
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

NR  = 1000
NE  =  500

with memory ( 'Create initial dataset' , logger = logger ) as dm0 :
    
    for r in progress_bar ( range ( NR ) , description = "Create:")  :
        
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



# =============================================================================
## (2)  Bootstrapping 
# =============================================================================
with memory ( 'Bootstrapping' , logger  = logger ) as dm1 :

    size = 10
    for i , ds in progress_bar ( enumerate ( dataset.bootstrap ( size  , extended = True , wrap = True ) ) ,
                                 max_value   = size         , 
                                 description = "Bootstrap:" ) :

        del ds 

# =============================================================================
## (3) Jackknife  
# =============================================================================
with memory ( 'Jackknife' , logger  = logger ) as dm2 :

    first , last = 0, 10
    for i , ds in progress_bar ( enumerate ( dataset.jackknife ( first , last , wrap = True ) ) ,
                                 max_value   = 100          , 
                                 description = "Jackknife:" ) :
        del ds 

# ==============================================================================
## delete it! 
with memory ( 'Delete Initial dataset' , logger = logger ) as dm3 :
    
    ## dataset = Ostap.MoreRooFit.delete_data ( dataset )
    ## del dataset 
    ## dataset = data_ptr ( dataset )
    del dataset 


# ==============================================================================
rows = [ ( "Action" , "%s %s memory [MB]" % ( delta_symbol , ram ) ) ] 

row  = "Create"    , '%+.2f' % dm0.delta
rows.append ( row ) 

row  = "Bootstrap" , '%+.2f' % dm1.delta
rows.append ( row ) 

row  = "Jackknife" , '%+.2f' % dm2.delta
rows.append ( row ) 

row  = "Delete"    , '%+.2f' % dm3.delta
rows.append ( row ) 

title = "Memory usage"
table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lc' )
logger.info ( '%s:\n%s' % ( title , table ) )


if True :
    raise TypeError ( " MEMORY %+.3f %+.3f %+.3f " % ( dm0.delta , dm1.delta , dm2.delta ) )
                 
assert dm1.delta < max ( 1 , 2 * dm0.delta ) , "There is some memory leak in Bootstrap: %s vs %s " % ( dm1.delta , dm0.delta )
assert dm2.delta < max ( 1 , 2 * dm0.delta ) , "There is some memory leak in Jackknife: %s vs %s " % ( dm1.delta , dm0.delta )
assert dm3.delta < max ( 1 , 2 * dm0.delta ) , "Memory is not released %s" % dm3.delta

# =============================================================================
##                                                                      The END 
# =============================================================================
