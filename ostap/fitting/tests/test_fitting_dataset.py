#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_dataset.py
# - It tests some decorators fore dataset
# ============================================================================= 
"""  - It tests some decorators fore dataset
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   builtins                 import range
from   ostap.core.core          import dsID
import ostap.fitting.roofit 
import ostap.logger.table       as     T
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_dataset' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

evt     = ROOT.RooRealVar ( 'Evt'    , '#event'        , 0 , 1000000 )
run     = ROOT.RooRealVar ( 'Run'    , '#run'          , 0 , 1000000 )
mass    = ROOT.RooRealVar ( 'Mass'   , 'mass-variable' , 0 , 100     )
weight  = ROOT.RooRealVar ( 'Weight' , 'some weight'   , -10 , 10 )

varset  = ROOT.RooArgSet  ( evt , run , mass , weight )
dataset = ROOT.RooDataSet ( dsID () , 'Test Data set-0' , varset )  

for r in range ( 100 ) :
    
    run.setVal ( r )
    for e in range ( 100 ) :
        
        evt .setVal   ( e )
        mass.setVal   ( random.uniform ( 0  , 10   ) )
        weight.setVal ( random.gauss   ( 1  , 0.1  ) )
        
        dataset.add ( varset )
        
        if ( 1 <= r <= 2  ) and 3 <=  e <= 3 :
            
            mass.setVal ( random.uniform ( 0 , 10 ) )
            dataset.add ( varset )

# =========================================================================================
weighted = dataset.makeWeighted ( 'Weight' )

# =========================================================================================
## (1) print datasets
# =========================================================================================
logger.info ( 'Print unweighted dataset:\n%s' % dataset .table ( prefix = '# ' ) )
logger.info ( 'Print   weighted dataset:\n%s' % weighted.table ( prefix = '# ' ) )

# =========================================================================================
## (2) loop over some subset of entries 
# =========================================================================================
logger.info  ('Loop over unweighted dataset') 
for index, entry , weight in dataset .loop ( '(Evt<10) && (Mass<10)' , first = 100 , last = 1000 , progress = True ) : pass

logger.info  ('Loop over   weighted dataset') 
for index, entry , weight in weighted.loop ( '(Evt<10) && (Mass<10)' , first = 100 , last = 1000 , progress = True ) : pass 

# =========================================================================================
## (3) loop and get certain information as arrays/rows 
# =========================================================================================
logger.info  ('Loop over unweighted dataset') 
for index, row , weight in dataset  .rows ( 'Mass, 2*Mass, Mass/2' , '(Evt<5) && (Mass<5)' , first = 100 , last = 1000 , progress = False  ) :
    print ( index, row , weight ) 
logger.info  ('Loop over  weighted dataset') 
for index, row , weight in weighted .rows ( 'Mass, 2*Mass, Mass/2' , '(Evt<5) && (Mass<5)' , first = 100 , last = 1000 , progress = False  ) :
    print ( index, row , weight ) 


# =============================================================================
##                                                                      The END 
# =============================================================================
