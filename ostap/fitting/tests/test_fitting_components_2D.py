#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_components_2D.py
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
from   ostap.core.core          import cpp, VE, dsID, rooSilent 
from   ostap.utils.timing       import timing
from   ostap.fitting.background import make_bkg 
from   ostap.plotting.canvas    import use_canvas 
from   ostap.utils.utils        import batch_env 
import ostap.fitting.models     as     Models 
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_components_2D' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch form environment 
batch_env ( logger )
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
N_ss = 5000*2
N_sb =  500*2
N_bs =  500*2
N_bb = 1000*2


random.seed(0)

## fill it : N_ss events  Gauss * Gauss *Gauss
for m in (m1,m2) : 
    for i in range(0,N_ss) : 
        m_x.value = m.gauss() 
        m_y.value = m.gauss() 
        dataset.add ( varset  )


## fill it : N_sb  events  Gauss * const * Gauss  
    for i in range(0,N_sb) : 
        m_x.value = m.gauss() 
        m_y.value = random.uniform ( *m_y.minmax() )  

        dataset.add ( varset  )

## fill it : N_bs  events  const * Gauss * Gauss
    for i in range(0,N_bs) : 
        m_x.value = random.uniform ( *m_x.minmax() ) 
        m_y.value = m.gauss() 
        dataset.add ( varset  )

## fill it : N_bb events  const * const *Gauss
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
sb_cmp=signal_x2*bkg_y

# B(x)*S(y) component 
bs_cmp= bkg_x*signal_y2

# B(x)*B(y) component 
bb_cmp=bkg_x*bkg_y

models = set() 

# =============================================================================
## Test  multi-component  3d fit'
def test_components_2D () :

    logger = getLogger  ('test_components_2D' )
    
    logger.info ('Test  multi-component  3d fit')
    
    model = Models.Fit2D (
        name    = 'fit_comp', 
        signal_x    = signal_x1, 
        signal_y    = signal_y1,
        bkg_1x  = bkg_x ,
        bkg_1y  = bkg_y ,
        components=[ss_cmp,sb_cmp,bs_cmp,bb_cmp]
        )
    
    with rooSilent() : 
        ## components
        model.SS.fix ( N_ss  )
        model.SB.fix ( N_sb  )
        model.SB.fix ( N_bs  )
        model.BB.fix ( N_bb  )

        model.C[0].fix ( 5000 )
        model.C[1].fix ( 500  )
        model.C[2].fix ( 500  )
        model.C[3].fix ( 1000 )
    
                
        r = model.fitTo ( dataset)

        model.SS.release (  )
        model.SB.release (  )
        model.BS.release (  )
        model.BB.release (  )

        
        model.C[0].release (  )
        model.C[1].release (  )
        model.C[2].release (  )
        model.C[3].release (  )
        
        r = model.fitTo ( dataset )
        r = model.fitTo ( dataset )
        r = model.fitTo ( dataset )

        with use_canvas ( 'test_components_2D' ) : 
            
            model.draw1 ( dataset )
            model.draw2 ( dataset )

    logger.info ( 'Model %s Fit result \n#%s ' % ( model.name , r ) ) 


    models.add  ( model )

# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
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

    test_components_2D    ()
    
    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()
  

# =============================================================================
##                                                                      The END 
# =============================================================================
