#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_components.py
# Test module 
# - It tests various multicomponents models 
# ============================================================================= 
""" Test module 
- It tests various multicomponents models 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import cpp, VE, dsID, rooSilent 
from   builtins                 import range
from   ostap.fitting.background import make_bkg 
from   ostap.plotting.canvas    import use_canvas 
from   ostap.utils.timing       import timing 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_components2_2d' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

## make simple test mass 
m_x     = ROOT.RooRealVar ( 'mass_x' , 'Some test mass(X)' , 0 , 10 )
m_y     = ROOT.RooRealVar ( 'mass_y' , 'Some test mass(Y)' , 0 , 10 )

## book very simple data set
varset  = ROOT.RooArgSet  ( m_x , m_y )
dataset = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset )  



m1 = VE(3,0.10**2)
m2 = VE(7,0.20**2)

## fill it with three gausissians, 5k events each
N_ss = 5000
N_sb =  500
N_bs =  500
N_bb = 1000


random.seed(0)

## fill it : 5000 events  Gauss * Gauss *Gauss
for m in (m1,m2) : 
    for i in range(0,N_ss) : 
        m_x.value = m.gauss() 
        m_y.value = m.gauss() 
        dataset.add ( varset  )


## fill it : 500 events  Gauss * const * Gauss  
    for i in range(0,N_sb) : 
        m_x.value = m.gauss() 
        m_y.value = random.uniform ( *m_y.minmax() )  

        dataset.add ( varset  )

## fill it : 500 events  const * Gauss * Gauss
    for i in range(0,N_bs) : 
        m_x.value = random.uniform ( *m_x.minmax() ) 
        m_y.value = m.gauss() 
        dataset.add ( varset  )

## fill it : 1000 events  const * const *Gauss
    for i in range(0,N_bb) :

        m_x.value  = random.uniform ( *m_x.minmax() )
        m_y.value  = random.uniform ( *m_y.minmax() )
        dataset.add ( varset  )


logger.info ('Dataset: %s' % dataset )  



## various fit components
signal_x1 = Models.Gauss_pdf ( 'G1x'  , xvar = m_x  , mean = m1.value() , sigma = m1.error() )  
signal_y1 = signal_x1.clone ( name='G1y'  , xvar = m_y   ) 

signal_x2 = Models.Gauss_pdf ( name='G2x'  , xvar = m_x  , mean = m2.value() , sigma = m2.error() )  
signal_y2 = signal_x2.clone ( name='G2y'  , xvar = m_y   ) 

bkg_x= make_bkg ( -1      , 'Bx' , m_x )
bkg_y= bkg_x.clone ( name= 'By' , xvar =m_y )

# S(x)*S(y) component 
ss_cmp=signal_x2*signal_y2

# S(x)*B(y) component 
sb=signal_x2*bkg_y

# B(x)*S(y) component 
bs= bkg_x*signal_y2

# B(x)*S(y)+ S(x)*B(y) component 
sb_cmp=sb+bs

# B(x)*B(y) component 
bb_cmp=bkg_x*bkg_y


models = set()

# =============================================================================
## Test  multi-component  3d fit'
def test_components_2dSymfit () :
    
    logger = getLogger ( 'test_components_2dSymfit' ) 
    logger.info ('Test  multi-component  2d Sym fit')
    
    model = Models.Fit2DSym (
        name    = 'fit_comp', 
        signal_x    = signal_x1, 
        signal_y    = signal_y1,
        bkg_1x  = bkg_x ,
        components=[ss_cmp,sb_cmp,bb_cmp]
        )
    
    with rooSilent() : 
        ## components
        model.SS.fix ( 5000 )
        model.SB.fix ( 1000 )
        model.BB.fix ( 1000 )

        model.C[0].fix ( 5000 )
        model.C[1].fix ( 1000 )
        model.C[2].fix ( 1000 )
    
                
        r = model.fitTo ( dataset , ncpu=8 )

        model.SS.release (  )
        model.SB.release (  )
        model.BB.release (  )

        
        model.C[0].release (  )
        model.C[1].release (  )
        model.C[2].release (  )
        
        r = model.fitTo ( dataset , ncpu=8 )

        with use_canvas ( 'test_components_2dSymfit' ) :
            model.draw1 (   dataset )
            model.draw2 (   dataset )
        
    logger.info ( 'Model %s Fit result \n#%s ' % ( model.name , r ) ) 


    models.add ( model )
    
# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        db['m_x,m_u'] = m_x, m_y 
        db['vars'] = varset
        db['dataset'  ] = dataset
        for m in models :
            db['model:' + m.name ] = m
            db['roo_tot:%s' % m.name ] = m.pdf
            for i,s in enumerate ( m.signals ) :
                db['roo_sig%d:%s' % ( i , m.name ) ] = s
            for i, b in enumerate ( m.backgrounds ) : 
                db['roo_bkg%d:%s' % ( i , m.name ) ] = s
            for a in m.alist1 : 
                db['cmp:%s/%s' % ( m.name , a.name ) ] = a
        db['models'   ] = models
        db.ls() 
    
# =============================================================================
if '__main__' == __name__ :


    with  timing( "2dSym-fit" ,   logger ) : 
        test_components_2dSymfit    () 

    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()

# =============================================================================
##                                                                      The END 
# =============================================================================
