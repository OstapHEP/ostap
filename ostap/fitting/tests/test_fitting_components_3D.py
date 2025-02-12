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
from   builtins                 import range
from   ostap.core.core          import cpp, VE, dsID, rooSilent 
from   ostap.fitting.background import make_bkg 
from   ostap.utils.timing       import timing
from   ostap.core.meta_info     import root_info 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait, batch_env 
import ostap.fitting.models     as     Models 
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_components_3D ' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch form environment 
batch_env ( logger )
# =============================================================================

models = set()

# =============================================================================
def test_fitting_components_3D () :

    logger = getLogger ( 'test_fitting_components_3D' ) 
    ## make simple test mass 
    m_x     = ROOT.RooRealVar ( 'mass_x' , 'Some test mass(X)' , 0 , 10 )
    m_y     = ROOT.RooRealVar ( 'mass_y' , 'Some test mass(Y)' , 0 , 10 )
    m_z     = ROOT.RooRealVar ( 'mass_z' , 'Some test mass(z)' , 0 , 10 )
    
    ## book very simple data set
    varset  = ROOT.RooArgSet  ( m_x , m_y,m_z )
    dataset = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset )  
    
    
    m1 = VE ( 3 , 0.10**2 )
    m2 = VE ( 7 , 0.20**2 )
    
    ## fill it
    
    N_sss = 5000
    N_ssb =  500
    N_sbs =  500
    N_bss  = 1000
    
    N_bbs = 500
    N_bsb = 100
    N_sbb = 100
    
    N_bbb = 250

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

    
    
    ## various fit components for main fit model :
    
    signal_x1 = Models.Gauss_pdf ( 'G1x'  , xvar = m_x  , mean = m1.value() , sigma = m1.error() )  
    signal_y1 = signal_x1.clone ( name='G1y'  , xvar = m_y   ) 
    signal_z1 = signal_x1.clone ( name='G1z'  , xvar = m_z   )
    
    bkg_x     = make_bkg ( -1      , 'Bx' , m_x        )
    bkg_y     = bkg_x.clone ( name = 'By' , xvar = m_y )
    bkg_z     = bkg_x.clone ( name = 'Bz' , xvar = m_z )

    
    ## construct other components 
    
    signal_x2 = Models.Gauss_pdf ( name='G2x'  , xvar = m_x  , mean = m2.value() , sigma = m2.error() )  
    signal_y2 = signal_x2.clone  ( name='G2y'  , xvar = m_y   ) 
    signal_z2 = signal_x2.clone  ( name='G2z'  , xvar = m_z   )

    
    ## S(x) * S(y) components  
    ss_cmp = signal_x2 * signal_y2

    ## S(x) * B(y) component 
    sb_cmp = signal_x2 * bkg_y

    ## B(x) * S(y) component 
    bs_cmp = bkg_x * signal_y2

    ## B(x) * B(y) component 
    bb_cmp = bkg_x * bkg_y

    ## S(x) * S(y) * S(z) component
    sss_cmp = ss_cmp * signal_z2
    
    ## S(x) * S(y) * B(z) component
    ssb_cmp = ss_cmp * bkg_z

    ## S(x) * B(y) * S(z) component
    sbs_cmp = sb_cmp * signal_z2

    ## S(x) * B(y) * S(z) component
    bss_cmp = bs_cmp * signal_z2
    
    
    ## S(x) * B(y) * B(z) component
    sbb_cmp = sb_cmp * bkg_z 

    ## B(x) * S(y) * B(z) component
    bsb_cmp = bs_cmp * bkg_z 

    ## B(x) * B(y) * S(z) component
    bbs_cmp = bs_cmp * signal_z2

    model = Models.Fit3D (
        name        = 'fit_comp', 
        signal_x    = signal_x1 , 
        signal_y    = signal_y1 ,
        signal_z    = signal_z1 ,
        bkg_1x      = bkg_x     ,
        bkg_1y      = bkg_y     ,
        bkg_1z      = bkg_z     ,
        suffix      = "_1"      ,
        ##
        components = [ sss_cmp ,
                       ssb_cmp ,
                       sbs_cmp ,
                       bss_cmp ,
                       sbb_cmp ,
                       bsb_cmp ,
                       bbs_cmp ] )
    
    
    model.SSS = N_sss
    
    model.SSB = N_ssb
    model.SBS = N_sbs
    model.BSS = N_bss
    
    model.BBS = N_bbs 
    model.BSB = N_bsb
    model.SBB = N_sbb
    
    model.BBB = N_bbb
    
    model.C   = N_sss , N_ssb , N_sbs , N_sbb , N_bss , N_sbs , N_bbs 
    
    r , _ = model.fitTo ( dataset , silent = True )
    r , _ = model.fitTo ( dataset , silent = True )
    

    with use_canvas ( 'test_fitting_components_3D' ) : 
        with wait ( 2 ) : model.draw1 (  dataset )
        with wait ( 2 ) : model.draw2 (  dataset )
        with wait ( 2 ) : model.draw3 (  dataset )
        
    logger.info ( 'Model %s Fit result\n%s ' % ( model.name , r.table (prefix = '# ') ) ) 


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

    
    test_fitting_components_3D () 
            
    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()

# =============================================================================
##                                                                      The END 
# =============================================================================
