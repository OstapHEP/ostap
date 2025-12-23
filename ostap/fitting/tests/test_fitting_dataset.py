#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_dataset.py
# - It tests some decorators fore dataset
# ============================================================================= 
"""  - It tests some decorators for dataset
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.core          import dsID, hID, Ostap 
from   ostap.plotting.canvas    import use_canvas 
from   ostap.utils.root_utils   import batch_env
from   ostap.utils.basic        import typename 
import ostap.logger.table       as     T
import ostap.fitting.roofit
import ostap.fitting.ds2numpy 
import ostap.math.models   
import ostap.trees.trees   
import ostap.histos.histos 
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
## set batch from environment 
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

for r in range ( 100 ) :
    
    run.setVal ( r )
    for e in range ( 100 ) :
        
        evt .setVal   ( e )
        mass.setVal   ( random.uniform      ( 0  , 10  ) )
        weight.setVal ( random.gauss        ( 1  , 0.1 ) )
        pt1.setVal     ( random.expovariate ( 1.0/15   ) ) 
        pt2.setVal     ( random.expovariate ( 1.0/25   ) ) 
        dataset.add ( varset )
        
        if ( 1 <= r <= 2  ) and 3 <=  e <= 3 :
            
            mass.setVal ( random.uniform ( 0 , 10 ) )
            dataset.add ( varset )
            
# =============================================================================
weighted = dataset.makeWeighted ( 'Weight' )

# =============================================================================
## (1) print datasets
# =============================================================================
logger.info ( 'Print         unweighted dataset:\n%s' % dataset .table ( prefix = '# ' ) )
logger.info ( 'Print           weighted dataset:\n%s' % weighted.table ( prefix = '# ' ) )

# =============================================================================
## (2) loop over some subset of entries 
# =============================================================================
logger.info  ('Loop over unweighted dataset') 
for index, entry , weight in dataset .loop ( '(Evt<10) && (Mass<10)' , first = 100 , last = 1000 , progress = True ) : pass

logger.info  ('Loop over   weighted dataset') 
for index, entry , weight in weighted.loop ( '(Evt<10) && (Mass<10)' , first = 100 , last = 1000 , progress = True ) : pass 

# =============================================================================
## (3) loop and get certain information as arrays/rows 
# =============================================================================
logger.info  ('Loop over unweighted dataset') 
for index, row , weight in dataset  .rows ( 'Mass, 2*Mass, Mass/2' , '(Evt<5) && (Mass<5)' , first = 100 , last = 1000 , progress = False  ) :
    print ( index, row , weight )
    
logger.info  ('Loop over  weighted dataset') 
for index, row , weight in weighted .rows ( 'Mass, 2*Mass, Mass/2' , '(Evt<5) && (Mass<5)' , first = 100 , last = 1000 , progress = False  ) :
    print ( index, row , weight ) 

# =============================================================================
## (4) subset 
# =============================================================================
ss1 = dataset  [ 1:500:10 ]
ws1 = weighted [ 1:500:10 ]
logger.info ( 'Print small   unweighted dataset:\n%s' % ws1.table ( prefix = '# ' ) )
logger.info ( 'Print small     weighted dataset:\n%s' % ws1.table ( prefix = '# ' ) )
 
# =============================================================================
## (5) subset/sample (unique) 
# =============================================================================
ss2 = dataset  . sample ( 100 )
ws2 = weighted . sample ( 100 ) 
logger.info ( 'Print small/1 unweighted sample:\n%s' % ss2.table ( prefix = '# ' ) )
logger.info ( 'Print small/1   weighted sample:\n%s' % ws2.table ( prefix = '# ' ) )

# =============================================================================
## (6) subset/sample (wih replacement)   
# =============================================================================
ss3 = dataset  . choice ( 100 )
ws3 = weighted . choice ( 100 ) 
logger.info ( 'Print small/2 unweighted sample:\n%s' % ss3.table ( prefix = '# ' ) )
logger.info ( 'Print small/2   weighted sample:\n%s' % ws3.table ( prefix = '# ' ) )

# =============================================================================
## (7) shuffle 
# =============================================================================
ss4 = ss3. shuffle ()
ws4 = ws3. shuffle () 
logger.info ( 'Print shuffle unweighted sample:\n%s' % ss4.table ( prefix = '# ' ) )
logger.info ( 'Print shuffle   weighted sample:\n%s' % ws4.table ( prefix = '# ' ) )

# =============================================================================
## (8) subset
# =============================================================================
ss5 = dataset .subset ( 'Mass,Pt1,Pt2' , cuts = '(Evt<5) && (Mass<5)' , first = 100 , last = 1000 )
ws5 = weighted.subset ( 'Mass,Pt1,Pt2' , cuts = '(Evt<5) && (Mass<5)' , first = 100 , last = 1000 )
logger.info ( 'Print subset unweighted sample:\n%s' % ss5.table ( prefix = '# ' ) )
logger.info ( 'Print subset   weighted sample:\n%s' % ws5.table ( prefix = '# ' ) )

# =============================================================================
## (9) symmetrize
# =============================================================================
ss6 = dataset .symmetrize ( [ 'Pt1' , 'Pt2' ] ) 
ws6 = weighted.symmetrize ( [ 'Pt1' , 'Pt2' ] ) 
logger.info ( 'Print symmetrized unweighted sample:\n%s' % ss6.table ( prefix = '# ' ) )
logger.info ( 'Print symmetrized   weighted sample:\n%s' % ws6.table ( prefix = '# ' ) )

# =============================================================================
## (10) remove entry 
# =============================================================================
ss7 = ss6 - 3 
ws7 = ws6 - 3 
logger.info ( 'Print remove #3 unweighted sample:\n%s' % ss7.table ( prefix = '# ' ) )
logger.info ( 'Print remove #3   weighted sample:\n%s' % ws7.table ( prefix = '# ' ) )

# =============================================================================
## (11) remove variable 
# =============================================================================
ss8 = ss7 - 'Pt1' 
ws8 = ws7 - 'Pt1'
logger.info ( 'Print remove #Pt1 unweighted sample:\n%s' % ss8.table ( prefix = '# ' ) )
logger.info ( 'Print remove #Pt1   weighted sample:\n%s' % ws8.table ( prefix = '# ' ) )

# =============================================================================
## (12) remove variables 
# =============================================================================
ss9 = ss7 - 'Pt1,Pt2' 
ws9 = ws7 - 'Pt1,Pt2'
logger.info ( 'Print remove #Pt  unweighted sample:\n%s' % ss9.table ( prefix = '# ' ) )
logger.info ( 'Print remove #Pt    weighted sample:\n%s' % ws9.table ( prefix = '# ' ) )

# =============================================================================
## (13) remove variables 
# =============================================================================
ss10 = ss7 - ('Pt1','Pt2')
ws10 = ws7 - ('Pt1','Pt2')
logger.info ( 'Print remove #Pt1,Pt2 unweighted sample:\n%s' % ss10.table ( prefix = '# ' ) )
logger.info ( 'Print remove #Pt1,Pt2   weighted sample:\n%s' % ws10.table ( prefix = '# ' ) )

# =============================================================================
## (14) project 
# =============================================================================
hd = ROOT.TH1D ( hID() , 'Pt projection' , 50 , 0 , 100 )
hw = ROOT.TH1D ( hID() , 'Pt projection' , 50 , 0 , 100 )

hd = dataset .project ( hd  , 'Pt1' , cuts = 'Mass<5' )
hw = weighted.project ( hw  , 'Pt1' , cuts = 'Mass<5' )

hd.draw()
hw.draw() 

# =============================================================================
## (15) project 
# =============================================================================

hd = dataset .project ( hd  , 'Pt1,Pt2' , cuts = 'Mass<5' )
hw = weighted.project ( hw  , 'Pt1,Pt2' , cuts = 'Mass<5' )

hd.draw()
hw.draw() 

# =============================================================================
## (16) parameterise (==project) 
# =============================================================================

ld = Ostap.Math.LegendreSum ( 12 , 0 , 100 )
lw = Ostap.Math.LegendreSum ( 12 , 0 , 100 )

ld = dataset .project ( ld  , 'Pt1' , cuts = 'Mass<5' )
lw = weighted.project ( lw  , 'Pt1' , cuts = 'Mass<5' )

# =============================================================================
## (17) plot 
# =============================================================================
with use_canvas ( "test_fitting_dataset: dataset.draw" , wait = 2 ) :
    
    hd = dataset.draw ( 'Pt1' , cuts = 'Mass<5' , color = 2 , xmin = 0 , xmax = 100 )

    hd.draw ( 'same' , color = 2 )
    ld.draw ( 'same' , color = 2 , width = 3 )
    
with use_canvas ( "test_fitting_dataset: weighted.draw" , wait  =2 ) :

    hw = weighted.draw ( 'Pt1' , cuts = 'Mass<5' , color = 4 , xmin = 0 , xmax = 100 )
    
    hw.draw ( 'same' , color = 4 )
    lw.draw ( 'same' , color = 4 , width = 3 )

# =============================================================================
## (18) conversion to TTree/TChain 
# =============================================================================
## chd = dataset .ds2tree().chain 
## chw = weighted.ds2tree().chain 
chd = dataset .asTree()
chw = weighted.asTree()

logger.info ( 'Print tree/chain unweighted sample:\n%s' % chd.table ( prefix = '# ' ) )
logger.info ( 'Print tree/chain   weighted sample:\n%s' % chw.table ( prefix = '# ' ) )


hd1 = hd.clone()
hd2 = hd.clone()

hd  = dataset .project ( hd  , 'Pt1,Pt2' , cuts = 'Mass<5' )
hd1 = dataset .project ( hd1 , 'Pt1'     , cuts = 'Mass<5' )
hd2 = dataset .project ( hd2 , 'Pt2'     , cuts = 'Mass<5' )

hd.draw(color = 2)
hd1.draw( 'same' , color = 4)
hd2.draw( 'same' , color = 8)

# =============================================================================
## (19) conversion to numpy: 'to_numpy'
# =============================================================================
ds = dataset  [:50] 
ws = weighted [:50]

vars = 'Mass' , 'Pt1' , 'Pt2'

rr = ds.to_numpy ( vars )
rw = ws.to_numpy ( vars )

logger.info ( 'Using `ROOT.RooDataset.to_numpy` method' ) 
logger.info ( 'to_numpy unweighted sample: type=%s names: %s' % ( typename ( rr ) , ','.join ( sorted( str ( v ) for v in rr ) ) ) ) 
logger.info ( 'to_numpy   weighted sample: type=%s names: %s' % ( typename ( rw ) , ','.join ( sorted( str ( v ) for v in rw ) ) ) ) 


# =============================================================================
## (20) conversion to numpy: 'tonumpy'
# =============================================================================
ds = dataset  [:15] 
ws = weighted [:15]

vars = 'Mass' , 'Pt1' 

rr = ds.tonumpy ( vars )
rw = ws.tonumpy ( vars )

logger.info ( 'Using `tonumpy/ds2numpy` method/function' ) 
logger.info ( 'tonumpy unweighted sample: type=%s shape:%s names: %s\n%s' % ( typename ( rr ) ,
                                                                              rr.shape        , 
                                                                              ','.join ( rr.dtype.names ) , rr ) )
logger.info ( 'tonumpy  weighted sample: type=%s shape:%s names: %s\n%s' % ( typename ( rr ) ,
                                                                             rr.shape        , 
                                                                             ','.join ( rr.dtype.names ) , rr ) )

# =============================================================================
## (21) conversion to numpy: 'tonumpy'
# =============================================================================
ds = dataset  [:15] 
ws = weighted [:15]

vars = 'Mass' , 'Pt1' 

rr , rrw = ds.tonumpy ( vars , weight_split = True ) 
rw , rww = ws.tonumpy ( vars , weight_split = True ) 

logger.info ( 'Using `tonumpy/ds2numpy(weight_split=True)` method/function' ) 
logger.info ( 'tonumpy unweighted sample: type=%s shape:%s names: %s\n%s' % ( typename ( rr ) ,
                                                                              rr.shape        , 
                                                                              ','.join ( rr.dtype.names ) , rr ) )
logger.info ( 'tonumpy  weighted sample: type=%s shape:%s names: %s\n%s' % ( typename ( rr ) ,
                                                                             rr.shape        , 
                                                                             ','.join ( rr.dtype.names ) , rr ) )

# =============================================================================
## (22) conversion to numpy: 'tonumpy'
# =============================================================================
ds = dataset  [:15] 
ws = weighted [:15]

vars = 'Mass' , 'Pt1' 

rr = ds.tonumpy ( vars , structured = False )
rw = ws.tonumpy ( vars , structured = False )
logger.info ( 'Using `tonumpy/ds2numpy(structured=False)` method/function' ) 
logger.info ( 'tonumpy unweighted sample: type=%s shape:%s\n%s' % ( typename ( rr ) , rr.shape , rr ) )
logger.info ( 'tonumpy   weighted sample: type=%s shape:%s\n%s' % ( typename ( rr ) , rr.shape , rr ) )

# =============================================================================
## (23) conversion to numpy: 'tonumpy'
# =============================================================================
ds = dataset  [:15] 
ws = weighted [:15]

vars = 'Mass' , 'Pt1' 

rr , rrw = ds.tonumpy ( vars , structured = False , weight_split = True ) 
rw , rww = ws.tonumpy ( vars , structured = False , weight_split = True ) 
logger.info ( 'Using `tonumpy/ds2numpy(structured=False,weight_split=True)` method/function' ) 
logger.info ( 'tonumpy unweighted sample: type=%s shape:%s\n%s' % ( typename ( rr ) , rr.shape , rr ) )
logger.info ( 'tonumpy   weighted sample: type=%s shape:%s\n%s' % ( typename ( rr ) , rr.shape , rr ) )

# =============================================================================
## (24) conversion to numpy: 'slice'
# =============================================================================
ds = dataset  [:10] 
ws = weighted [:10]

vars = 'Mass' , 'Pt1' 

rr , rrw = ds.slice ( vars ) 
rw , rww = ws.slice ( vars )

logger.info ( 'Using `slice` method/function' ) 
logger.info ( 'slice   unweighted sample: type=%s shape:%s\n%s' % ( typename ( rr ) , rr.shape , rr ) )
logger.info ( 'slice     weighted sample: type=%s shape:%s\n%s' % ( typename ( rr ) , rr.shape , rr ) )

# =============================================================================
## (25) conversion to numpy: 'slice'
# =============================================================================
ds = dataset  [:10] 
ws = weighted [:10]

vars = 'Mass' , 'Pt1' 

rr , rrw = ds.slice ( vars , structured = False ) 
rw , rww = ws.slice ( vars , structured = False )

logger.info ( 'Using `slice(structured=False)` method/function' ) 
logger.info ( 'slice   unweighted sample: type=%s shape:%s\n%s' % ( typename ( rr ) , rr.shape , rr ) )
logger.info ( 'slice     weighted sample: type=%s shape:%s\n%s' % ( typename ( rr ) , rr.shape , rr ) )

# =============================================================================
## (26) conversion to numpy: 'slice'
# =============================================================================
ds = dataset  [:10] 
ws = weighted [:10]

vars = 'Mass' , 'Pt1' 

rr , rrw = ds.slice ( vars , structured = False , transpose = True ) 
rw , rww = ws.slice ( vars , structured = False , transpose = True )

logger.info ( 'Using `slice(structured=False,transpose=True)` method/function' ) 
logger.info ( 'slice   unweighted sample: type=%s shape:%s\n%s' % ( typename ( rr ) , rr.shape , rr ) )
logger.info ( 'slice     weighted sample: type=%s shape:%s\n%s' % ( typename ( rr ) , rr.shape , rr ) )


# =============================================================================
##                                                                       The END 
# =============================================================================
