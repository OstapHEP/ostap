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
from   ostap.core.core          import cpp, VE, dsID
from   ostap.logger.utils       import rooSilent
from   builtins                 import range
from   ostap.fitting.background import make_bkg 
from   ostap.utils.timing       import timing 
from   ostap.core.meta_info     import root_info 
from   ostap.plotting.canvas    import use_canvas 
import ROOT, random 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_components2_3D' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

# =============================================================================
def test_fitting_components2_3D () :

    logger = getLogger( 'test_fitting_components2_3D' )
    
    ## make simple test mass 
    m_x     = ROOT.RooRealVar ( 'mass_x' , 'Some test mass(X)' , 0 , 10 )
    m_y     = ROOT.RooRealVar ( 'mass_y' , 'Some test mass(Y)' , 0 , 10 )
    m_z     = ROOT.RooRealVar ( 'mass_z' , 'Some test mass(z)' , 0 , 10 )

    ## book very simple data set
    varset  = ROOT.RooArgSet  ( m_x , m_y,m_z )
    dataset = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset )  
    
    
    m1 = VE ( 3 , 0.10**2 )
    m2 = VE ( 7 , 0.20**2 )
    
    ## fill it with three gausissians, 5k events each
    N_sss = 5000
    N_ssb =  500
    N_sbs =  N_ssb
    N_bss =  N_ssb 
    
    N_bbs =  100
    N_bsb =  N_bbs
    N_sbb =  N_bbs 

    N_bbb =  500 
    
    random.seed(0)
    ## fill it : 5000 events  Gauss * Gauss *Gauss
    for m in (m1,m2) :

        ##  S x S x S 
        for i in range ( N_sss ) :
            
            m_x.value = m.gauss() 
            m_y.value = m.gauss() 
            m_z.value = m.gauss()
            
            dataset.add ( varset  )
            
        ## S x S x B 
        for i in range(0,N_bss) :
            
            m_x.value = m.gauss() 
            m_y.value = m.gauss() 
            m_z.value = random.uniform ( *m_z.minmax() )
            
            dataset.add ( varset  )
            
        ## S x B x S 
        for i in range ( N_ssb ) :
            
            m_x.value = m.gauss() 
            m_y.value = random.uniform ( *m_y.minmax() )  
            m_z.value = m.gauss() 
            
            dataset.add ( varset  )

        ## B x S x S 
        for i in range ( N_sbs ) : 
            m_x.value = random.uniform ( *m_x.minmax() ) 
            m_y.value = m.gauss() 
            m_z.value = m.gauss() 
            dataset.add ( varset  )

            
        ## B x B X S  
        for i in range ( N_sbb ) :
            
            m_x.value  = random.uniform ( *m_x.minmax() )
            m_y.value  = random.uniform ( *m_y.minmax() )
            m_z.value = m.gauss() 
            dataset.add ( varset  )

        ## B  x S x B 
        for i in range ( N_bsb ) : 

            m_x.value = random.uniform ( *m_x.minmax() ) 
            m_y.value = m.gauss() 
            m_z.value = random.uniform ( *m_y.minmax() )
            dataset.add ( varset  )

        ## S x B x B 
        for i in range ( N_sbb ) :
            
            m_x.value = m.gauss() 
            m_y.value = random.uniform ( *m_y.minmax() )  
            m_z.value = random.uniform ( *m_y.minmax() )
            
            dataset.add ( varset  )

        ## B x B x B 
        for i in range ( N_bbb ) :
            
            m_x.value  = random.uniform ( *m_x.minmax() )
            m_y.value  = random.uniform ( *m_y.minmax() )
            m_z.value = random.uniform ( *m_y.minmax() )
            
            dataset.add ( varset  )
            
    logger.info ('Dataset:\n%s' % dataset.table ( prefix = '# ') )   


    
    
    ## various fit components
    signal_x1 = Models.Gauss_pdf ( 'G1x'  , xvar = m_x  , mean = m1.value() , sigma = m1.error() )  
    signal_y1 = signal_x1.clone  ( name='G1y'  , xvar = m_y   ) 
    signal_z1 = signal_x1.clone  ( name='G1z'  , xvar = m_z   )
    
    bkg_x     = make_bkg    ( -1   , 'Bx' , m_x )
    bkg_y     = bkg_x.clone ( name = 'By' , xvar = m_y )
    bkg_z     = bkg_x.clone ( name = 'Bz' , xvar = m_z )
    
    signal_x2 = Models.Gauss_pdf ( name='G2x'  , xvar = m_x  , mean = m2.value() , sigma = m2.error() )  
    signal_y2 = signal_x2.clone  ( name='G2y'  , xvar = m_y   ) 
    signal_z2 = signal_x2.clone  ( name='G2z'  , xvar = m_z   )



    # S(x) * S(y) component 
    ss_cmp = signal_x2 * signal_y2
    
    # S(x) * B(y) component 
    sb_cmp = signal_x2 * bkg_y

    # B(x) * S(y) component 
    bs_cmp = bkg_x * signal_y2
    
    # B(x) * B(y) component 
    bb_cmp = bkg_x * bkg_y

    # S(x)*S(y)*S(z) component
    sss_cmp = ss_cmp * signal_z2
    
    # S(x) * B(y) * S(z) + B(x) * S(y) * S(z) + S(x) * S(y) * B(z) component
    ssb_    = ss_cmp * bkg_z
    sbs_    = sb_cmp * signal_z2
    bss_    = bs_cmp * signal_z2

    ssb_cmp = ssb_ + sbs_ + bss_ 

    # S(x)*B(y)*B(z) + B(x)*S(y)*B(z) + B(x)*B(y)*S(z) component
    sbb_ = sb_cmp * bkg_z
    bsb_ = bs_cmp * bkg_z
    bbs_ = bs_cmp * signal_z2 

    bbs_cmp = sbb_ + bsb_ + bsb_ 
    
    model = Models.Fit3DSym (
        name        = 'fitSym_comp', 
        signal_x    = signal_x1, 
        signal_y    = signal_y1,
        signal_z    = signal_z1,
        bkg_1x      = bkg_x ,
        bkg_2x      = 'clone' ,
        suffix      = '3D'    ,
        components  = [ sss_cmp , ssb_cmp , bbs_cmp  ]
        )


    ## components
    
    model.SSS = N_sss 
    model.SSB = N_ssb * 3 
    model.SBB = N_sbb * 3 
    model.BBB = N_bbb 
    
    model.C  = N_sss, N_ssb * 3 , N_sbb * 3
    
    r = model.fitTo ( dataset , silent = True  )
    r = model.fitTo ( dataset , silent = True  )

    with use_canvas ( 'test_fitting_components2_3D' ) : 
        model.draw1 ( dataset )
        model.draw2 ( dataset )
        model.draw3 ( dataset )
    
    logger.info ( 'Model %s Fit result\n%s ' % ( model.name , r.table ( prefix = '# ' ) ) ) 

    
# =============================================================================
if '__main__' == __name__ :

    with timing ( "3dSym-fit" , logger ) :
        test_fitting_components2_3D () 
        
# =============================================================================
##                                                                      The END 
# =============================================================================
