#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_plotting_fitdraw.py
# Test module for fit-draw machinery 
# ============================================================================= 
""" Test module for fit-draw machinery 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.utils.timing       import timing
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.root_utils   import batch_env 
from   ostap.fitting.background import make_bkg
from   ostap.logger.colorized   import attention
import ostap.logger.table       as     T
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_plotting_fitdraw' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
batch_env ( logger ) 
# =============================================================================

x = ROOT.RooRealVar  ( 'x' , 'x-obsevable' , 0 , 10 )

## "Signals"
s1 = Models.Gauss_pdf ( 'S1' , xvar = x , mean = 3 , sigma = 0.30 ) 
s2 = Models.Gauss_pdf ( 'S2' , xvar = x , mean = 7 , sigma = 0.30 )

s1.mean .fix()
s2.mean .fix()
s1.sigma.fix()
s2.sigma.fix()

## "Components" 
c1 = Models.Gauss_pdf ( 'C1' , xvar = x , mean = 1 , sigma = 0.70 ) 
c2 = Models.Gauss_pdf ( 'C2' , xvar = x , mean = 5 , sigma = 0.70 )
c3 = Models.Gauss_pdf ( 'C3' , xvar = x , mean = 9 , sigma = 0.70 )

c1.mean .fix()
c2.mean .fix()
c3.mean .fix()
c1.sigma.fix()
c2.sigma.fix()
c3.sigma.fix()

## Backgrounds"
b1 = Models.GammaDist_pdf ( 'B1' , xvar = x , k = 2.5   , theta = 2 )
b1.k    .fix()
b1.theta.fix() 
b2 = Models.Bkg_pdf       ( 'B2' , xvar = x , power = 1 )
b2.tau.fix()

model0 = Models.Fit1D (
    signals     = [s1,s2]    ,  
    others      = [c1,c2,c3] ,
    backgrounds = [b1,b2]    ,
    suffix      = 'M0'       )

model0.S = 1000 , 500
model0.C =  500 , 500 , 500
model0.B = 1000 , 1000


data = model0.generate ( 10000 ) 

models  = set()
results = set()
plots   = set()


model1 = Models.Fit1D (
    signals             = [s1,s2]    ,  
    others              = [c1,c2,c3] ,
    backgrounds         = [b1,b2]    ,
    combine_signals     = True       ,
    combine_others      = True       ,
    combine_backgrounds = True       ,    
    suffix              = 'M1'       )

model1.S  = 3000 
model1.C  = 3000 
model1.B  = 4000
model1.fS = 0.5 
model1.fC = 0.333 , 0.500 
model1.fB = 0.5 

# =============================================================================
def test_fitdraw1 () :

    logger = getLogger( 'test_fitdraw1' )

    with use_canvas ( 'test_fitfraw1' ) :
        _ , _  = model0.fitTo ( data , silent = True ) 
        _ , _  = model0.fitTo ( data , silent = True ) 
        r , f  = model0.fitTo ( data , silent = True , draw = True , nbins = 100  ) 
        title = 'fitdraw1'
        logger.info ( '%s:\n%s' % ( title , r.table ( title = title , prefix = '# ' ) ) )
        
    models .add ( model0 )
    results.add ( r  )
    plots  .add ( f  )

# =============================================================================
def test_fitdraw2 () :

    logger = getLogger( 'test_fitdraw2' )

    with use_canvas ( 'test_fitfraw2/default' ) :
        _ , _  = model1.fitTo ( data , silent = True ) 
        _ , _  = model1.fitTo ( data , silent = True ) 
        r , f  = model1.fitTo ( data , silent = True , draw = True , nbins = 100  ) 
        title = 'fitdraw2'
        logger.info ( '%s:\n%s' % ( title , r.table ( title = title , prefix = '# ' ) ) )
        
    models .add ( model1 )
    results.add ( r  )
    plots  .add ( f  )

    with use_canvas ( 'test_fitfraw2/combined-only-1' ) :
        f1 = model1.draw ( dataset = data ,
                           nbins   = 100  ,
                           ##
                           draw_combined_signal     = True  ,   
                           draw_signals             = False , 
                           draw_combined_component  = True  ,  
                           draw_components          = False , 
                           draw_combined_background = True  ,  
                           draw_backgrounds         = False ) 
        plots  .add ( f1 )

    import ostap.plotting.fit_draw  as FD 
    
    with use_canvas ( 'test_fitfraw2/combined-only-2' ) :
        f2 = model1.draw ( dataset = data ,
                           nbins   = 100  ,
                           ##
                           combined_signal_style     = FD.default_signal_style     , 
                           combined_background_style = FD.default_background_style , 
                           combined_component_style  = FD.default_component_style  , 
                           ##
                           draw_combined_signal     = True  ,   
                           draw_signals             = False , 
                           draw_combined_component  = True  ,  
                           draw_components          = False , 
                           draw_combined_background = True  ,  
                           draw_backgrounds         = False ) 
        plots  .add ( f2 )

        

# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        db['x'        ] = x  
        db['dataset'  ] = data
        for m in models :
            db['model:'     + m.name   ] = m
            db['roo_tot:%s' % m.name ] = m.pdf
            for s in m.signals              : db['signal           : %s' % s.name ] = s
            for c in m.components           : db['component        : %s' % c.name ] = c
            for b in m.backgrounds          : db['background       : %s' % b.name ] = b
            for s in m.combined_signals     : db['signal(comb)     : %s' % s.name ] = s
            for c in m.combined_components  : db['component(comb)  : %s' % c.name ] = c
            for b in m.combined_backgrounds : db['background(comb) : %s' % b.name ] = b
                
        db['models'   ] = models
        for r in results :
            db ['result:%s' % r.name ] = r
        db['results'   ] = results
        for p in plots :
            db [ ' plot:%s' % p.name ] = p  
        db['plots'     ] = plots 
        db.ls() 

        

# =============================================================================
if '__main__' == __name__ :

    with timing ('test_fitdraw1'       , logger ) :
        test_fitdraw1 () 
    
    with timing ('test_fitdraw2'       , logger ) :
        test_fitdraw2 () 
    
    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()


        

# =============================================================================
##                                                                      The END 
# =============================================================================
