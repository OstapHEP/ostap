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
import ROOT, random
import ostap.fitting.roofit 
import ostap.fitting.models as     Models 
from   ostap.core.core      import cpp, VE, dsID
from   ostap.logger.utils   import rooSilent
from   builtins             import range
from ostap.fitting.background import make_bkg 

# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_components' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

## make simple test mass 
m_x     = ROOT.RooRealVar ( 'mass_x' , 'Some test mass(X)' , 0 , 10 )
m_y     = ROOT.RooRealVar ( 'mass_y' , 'Some test mass(Y)' , 0 , 10 )
m_z     = ROOT.RooRealVar ( 'mass_z' , 'Some test mass(z)' , 0 , 10 )

## book very simple data set
varset  = ROOT.RooArgSet  ( m_x , m_y,m_z )
dataset = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset )  



m1 = VE(3,0.10**2)
m2 = VE(7,0.20**2)

## fill it with three gausissians, 5k events each
N_sss = 5000
N_ssb =  500
N_sbs =  500
N_sbb = 1000

N_bss = 500
N_bsb =  100
N_bbs =  100
N_bbb = 250

random.seed(0)

## fill it : 5000 events  Gauss * Gauss *Gauss
for m in (m1,m2) : 
    for i in range(0,N_sss) : 
        m_x.value = m.gauss() 
        m_y.value = m.gauss() 
        m_z.value = m.gauss() 
        dataset.add ( varset  )


## fill it : 500 events  Gauss * const * Gauss  
    for i in range(0,N_ssb) : 
        m_x.value = m.gauss() 
        m_y.value = random.uniform ( *m_y.minmax() )  
        m_z.value = m.gauss() 

        dataset.add ( varset  )

## fill it : 500 events  const * Gauss * Gauss
    for i in range(0,N_sbs) : 
        m_x.value = random.uniform ( *m_x.minmax() ) 
        m_y.value = m.gauss() 
        m_z.value = m.gauss() 
        dataset.add ( varset  )

## fill it : 1000 events  const * const *Gauss
    for i in range(0,N_sbb) :

        m_x.value  = random.uniform ( *m_x.minmax() )
        m_y.value  = random.uniform ( *m_y.minmax() )
        m_z.value = m.gauss() 
        dataset.add ( varset  )
## fill it : 500 events    Gauss * Gauss * const 
    for i in range(0,N_bss) : 
        m_x.value = m.gauss() 
        m_y.value = m.gauss() 
        m_z.value = random.uniform ( *m_z.minmax() )
        dataset.add ( varset  )
## fill it : 100 events  Gauss * const * const  
    for i in range(0,N_bsb) : 
        m_x.value = m.gauss() 
        m_y.value = random.uniform ( *m_y.minmax() )  
        m_z.value = random.uniform ( *m_y.minmax() )
        dataset.add ( varset  )

## fill it : 100 events  const * Gauss * const
    for i in range(0,N_bbs) : 
        m_x.value = random.uniform ( *m_x.minmax() ) 
        m_y.value = m.gauss() 
        m_z.value = random.uniform ( *m_y.minmax() )
        dataset.add ( varset  )

## fill it : 250 events  const * const * const
    for i in range(0,N_bbb) :
        m_x.value  = random.uniform ( *m_x.minmax() )
        m_y.value  = random.uniform ( *m_y.minmax() )
        m_z.value = random.uniform ( *m_y.minmax() )
        dataset.add ( varset  )


logger.info ('Dataset: %s' % dataset )  



## various fit components
signal_x1 = Models.Gauss_pdf ( 'G1x'  , xvar = m_x  , mean = m1.value() , sigma = m1.error() )  
signal_y1 = signal_x1.clone ( name='G1y'  , xvar = m_y   ) 
signal_z1 = signal_x1.clone ( name='G1z'  , xvar = m_z   )

signal_x2 = Models.Gauss_pdf ( name='G2x'  , xvar = m_x  , mean = m2.value() , sigma = m2.error() )  
signal_y2 = signal_x2.clone ( name='G2y'  , xvar = m_y   ) 
signal_z2 = signal_x2.clone ( name='G2z'  , xvar = m_z   )

bkg_x= make_bkg ( -1      , 'Bx' , m_x )
bkg_y= bkg_x.clone ( name= 'By' , xvar =m_y )
bkg_z= bkg_x.clone ( name='Bz' , xvar =m_z )

# S(x)*S(y) component 
ss_cmp=signal_x2*signal_y2

# S(x)*B(y) component 
sb_cmp=signal_x2*bkg_y

# B(x)*S(y) component 
bs_cmp= bkg_x*signal_y2

# B(x)*B(y) component 
bb_cmp=bkg_x*bkg_y

# S(x)*S(y)*S(z) component
sss_cmp=ss_cmp*signal_z2

# S(x)*B(y)*S(z)+ B(x)*S(y)*S(z)+ S(x)*S(y)*B(z) component
ssb=sb_cmp*signal_z2
sbs=bs_cmp*signal_z2
bss=ss_cmp*bkg_z
ssb_=ssb+sbs

ssb_cmp=ssb_+bss

# S(x)*B(y)*B(z)+ B(x)*S(y)*B(z)+ B(x)*B(y)*S(z) component
sbb=bb_cmp*signal_z2
bsb=sb_cmp*bkg_z
bbs=bs_cmp*bkg_z
sbb_=sbb+bsb

sbb_cmp=bsb+bbs

# B(x)*B(y)*B(z) component
bbb_cmp=bb_cmp*bkg_z






def test_comp_3dSymfit () :
    
    logger.info ('Test  multi-component  3d Sym fit')
    
    model = Models.Fit3DSym (
        name    = 'fitSym_comp', 
        signal_x    = signal_x1, 
        signal_y    = signal_y1,
        signal_z    = signal_z1,
        bkg_1x  = bkg_x ,
        bkg_2x  = 'clone' ,        
        components=[sss_cmp,ssb_cmp,sbb_cmp,bbb_cmp]
        )
    
    with rooSilent() : 
        ## components
        model.SSS.fix ( 5000 )
        model.SSB.fix ( 2000 )
        model.SBB.fix ( 1200 )
        model.BBB.fix ( 250 )

        model.C[0].fix ( 5000 )
        model.C[1].fix ( 2000 )
        model.C[2].fix ( 1000 )
        model.C[3].fix ( 250 )
    
                
        r = model.fitTo ( dataset , ncpu=8 )

        model.SSS.release (  )
        model.SSB.release (  )
        model.SBB.release (  )
        model.BBB.release (  )

        
        model.C[0].release (  )
        model.C[1].release (  )
        model.C[2].release (  )
        model.C[3].release (  )
        
        r = model.fitTo ( dataset , ncpu=8 )

    logger.info ( 'Model %s Fit result \n#%s ' % ( model.name , r ) ) 

    
# =============================================================================
if '__main__' == __name__ :

    test_comp_3dSymfit ()
# =============================================================================
# The END 
# =============================================================================
