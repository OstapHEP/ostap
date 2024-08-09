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
pt      = ROOT.RooRealVar ( 'Pt'     , 'pt-variable'   , 0 , 100     )
weight  = ROOT.RooRealVar ( 'Weight' , 'some weight'   , -10 , 10 )

varset  = ROOT.RooArgSet  ( evt , run , mass , pt , weight )
dataset = ROOT.RooDataSet ( dsID () , 'Test Data set-0' , varset )  

for r in range ( 100 ) :
    
    run.setVal ( r )
    for e in range ( 100 ) :
        
        evt .setVal   ( e )
        mass.setVal   ( random.uniform ( 0  , 10  ) )
        weight.setVal ( random.gauss   ( 1  , 0.1 ) )
        pt.setVal     ( random.uniform ( 1  , 99  ) ) 
        dataset.add ( varset )
        
        if ( 1 <= r <= 2  ) and 3 <=  e <= 3 :
            
            mass.setVal ( random.uniform ( 0 , 10 ) )
            dataset.add ( varset )

# =========================================================================================
weighted = dataset.makeWeighted ( 'Weight' )

# =========================================================================================
## (1) print datasets
# =========================================================================================
logger.info ( 'Print         unweighted dataset:\n%s' % dataset .table ( prefix = '# ' ) )
logger.info ( 'Print           weighted dataset:\n%s' % weighted.table ( prefix = '# ' ) )

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

# =========================================================================================
## (4) subset 
# =========================================================================================
ss1 = dataset  [ 1:500:10 ]
ss2 = weighted [ 1:500:10 ]
logger.info ( 'Print small   unweighted dataset:\n%s' % ss1.table ( prefix = '# ' ) )
logger.info ( 'Print small     weighted dataset:\n%s' % ss2.table ( prefix = '# ' ) )

# =========================================================================================
## (5) subset/sample (unique) 
# =========================================================================================
ss1 = dataset  . sample ( 100 )
ss2 = weighted . sample ( 100 ) 
logger.info ( 'Print small   unweighted sample:\n%s' % ss1.table ( prefix = '# ' ) )
logger.info ( 'Print small     weighted sample:\n%s' % ss2.table ( prefix = '# ' ) )

# =========================================================================================
## (6) subset/sample (wih r  
# =========================================================================================
ss1 = dataset  . choice ( 100 )
ss2 = weighted . choice ( 100 ) 
logger.info ( 'Print small   unweighted sample:\n%s' % ss1.table ( prefix = '# ' ) )
logger.info ( 'Print small     weighted sample:\n%s' % ss2.table ( prefix = '# ' ) )

# =========================================================================================
## (7) shuffle 
# =========================================================================================
ss1 = ss1. shuffle ()
ss2 = ss2. shuffle () 
logger.info ( 'Print shuffle unweighted sample:\n%s' % ss1.table ( prefix = '# ' ) )
logger.info ( 'Print shuffle weighted sample:\n%s' % ss2.table ( prefix = '# ' ) )




# =============================================================================
##                                                                      The END 
# =============================================================================
