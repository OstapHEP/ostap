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
from   __future__               import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import cpp, VE, dsID, rooSilent 
from   ostap.utils.timing       import timing
from   builtins                 import range
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait
from   ostap.fitting.background import make_bkg
from   ostap.logger.colorized   import attention
import ostap.logger.table       as     T
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


evt     = ROOT.RooRealVar ( 'Evt'  , '#event'        , 0 , 1000000 )
run     = ROOT.RooRealVar ( 'Run'  , '#run'          , 0 , 1000000 )
mass    = ROOT.RooRealVar ( 'Mass' , 'mass-variable' , 0 , 100     )

varset  = ROOT.RooArgSet  ( evt , run , mass )
dataset = ROOT.RooDataSet ( dsID() , 'Test Data set-0' , varset )  

for r in range ( 5 ) :
    run.setVal ( r ) 
    for e in range ( 5 ) :
        
        evt .setVal ( e )
        mass.setVal ( random.uniform ( 0 , 10 ) )
        
        dataset.add ( varset )
        
        if ( 1 <= r <= 2  ) and 3 <=  e <= 3 :
            
            mass.setVal ( random.uniform ( 0 , 10 ) )
            dataset.add ( varset )
            
            ## mass.setVal ( random.uniform ( 0 , 10 ) )
            ## dataset.add ( varset ) 


rows = [ ( '#' , 'Run' , 'Evt' , 'Mass' ) ] 
for i, e in enumerate ( dataset ) :

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
for i, e in enumerate ( ds0 ) :

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
        dscl.add ( dataset[i] )


rows = [ ( '#' , 'Run' , 'Evt' , 'Mass' ) ] 
for i, e in enumerate ( dscl ) :

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
