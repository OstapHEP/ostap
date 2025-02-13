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
from   ostap.utils.utils        import batch_env 
import ostap.logger.table       as     T
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_dusplicates' )
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
##                                                                      The END 
# =============================================================================
