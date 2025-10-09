#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_duplicates.py
# - It tests the removal of duplicated entries 
# ============================================================================= 
""" It tests the removal of duplicated entries 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.core          import dsID
from   ostap.utils.timing       import timing
from   ostap.utils.root_utils   import batch_env 
import ostap.logger.table       as     T
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_duplicates' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment  
batch_env ( logger )
# =============================================================================

evt     = ROOT.RooRealVar ( 'Evt'    , '#event'        , 0 , 1000000 )
run     = ROOT.RooRealVar ( 'Run'    , '#run'          , 0 , 1000000 )
mass    = ROOT.RooRealVar ( 'Mass'   , 'mass-variable' , 0 , 100     )
weight  = ROOT.RooRealVar ( 'Weight' , 'some weight'   , 0 , 2       )

varset  = ROOT.RooArgSet  ( evt , run , mass , weight )
dataset = ROOT.RooDataSet ( dsID () , 'Test Data set-0' , varset )  

for r in range ( 5 ) :
    
    run.setVal ( r )
    
    for e in range ( 5 ) :
        
        evt .setVal   ( e )
        mass.setVal   ( random.uniform ( 0 , 10 ) )
        weight.setVal ( random.uniform ( 0  , 2  ) )
        
        dataset.add ( varset )
        
        if ( 1 <= r <= 2  ) and 3 <=  e <= 3 :
            
            mass.setVal ( random.uniform ( 0 , 10 ) )
            dataset.add ( varset )

rows = [ ( '#' , 'Run' , 'Evt' , 'Mass' ) ] 
for i, item in enumerate ( dataset ) :
    e , _ = item 
    row = '%3d' % i , \
          '%3d' % float ( e.Run  ) , \
          '%3d' % float ( e.Evt  ) , \
          '%3d' % float ( e.Mass )
    rows.append ( row ) 
title = 'Original dataset'
logger.info ( '%s:\n%s' % ( title , T.table ( rows ,title = title , prefix = '# ' ) ) ) 

# =============================================================================
## (1) Eliminate (intermal) duplicates 
# =============================================================================
event_tag = 'Run' , 'Evt'

ds0 = dataset.make_unique ( event_tag , choice = 'first' )

rows = [ ( '#' , 'Run' , 'Evt' , 'Mass' ) ] 
for i, item in enumerate ( ds0 ) :
    e , _ = item 
    row = '%3d' % i , \
          '%3d' % float ( e.Run  ) , \
          '%3d' % float ( e.Evt  ) , \
          '%3d' % float ( e.Mass )
    rows.append ( row ) 
title = 'Unique dataset'
logger.info ( '%s:\n%s' % ( title , T.table ( rows ,title = title , prefix = '# ' ) ) ) 

dscl = dataset.emptyClone()
for dups in dataset.duplicates ( event_tag ) :
    for i in dups :
        entry , weight = dataset [ i ] 
        dscl.add ( entry  )

rows = [ ( '#' , 'Run' , 'Evt' , 'Mass' ) ] 
for i, item in enumerate ( dscl ) :
    e, w = item     
    row = '%3d' % i , \
          '%3d' % float ( e.Run  ) , \
          '%3d' % float ( e.Evt  ) , \
          '%3d' % float ( e.Mass )
    rows.append ( row ) 
title = 'Dataset with duplicates'
logger.info ( '%s:\n%s' % ( title , T.table ( rows ,title = title , prefix = '# ' ) ) ) 

# =============================================================================
## (2) Shared entries ("external duplicates") 
# =============================================================================
varset2  = ROOT.RooArgSet  ( evt , run , mass  )
dataset2 = ROOT.RooDataSet ( dsID () , 'Test Data set-2' , varset2 )  

for r in range ( 1000 ) :
    
    run.setVal ( r )
    
    for e in range ( 1000 ) :
        
        evt .setVal   ( e )
        mass.setVal   ( random.uniform ( 0 , 10 ) )
        
        dataset2.add ( varset )

        if  random.uniform ( 0 , 1 ) < 0.1  :

            mass.setVal ( random.uniform ( 0 , 10 ) )
            dataset2.add ( varset )

## Two dataset that have some common entries 
ds1 = dataset2.sample ( 0.50 )
ds2 = dataset2.sample ( 0.50 )

title = 'Dataset-1  (with shared entries)'
logger.info ( '%s:\n%s' % ( title , ds1.table ( prefix = '# ' ) ) ) 


ds1_shared = ds1.shared_data ( ds2 , event_tag , shared = True  , progress = True ) 
ds1_own    = ds1.shared_data ( ds2 , event_tag , shared = False , progress = True ) 

title = 'Dataset-1 (no   shared entries)'
logger.info ( '%s:\n%s' % ( title , ds1_own   .table ( prefix = '# ' ) ) ) 
title = 'Dataset-2 (only shared entries)'
logger.info ( '%s:\n%s' % ( title , ds1_shared.table ( prefix = '# ' ) ) ) 


n0 = len ( ds1 )
n1 = len ( ds1_shared )
n2 = len ( ds1_own    )

logger.info ( 'Sizes: %s %s %s : %s' % ( n0 , n1  , n2 , n0 - n1 - n2 ) ) 



# =============================================================================
##                                                                      The END 
# =============================================================================
