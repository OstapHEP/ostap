#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_in_range_2d.py
# Test module 
# - It tests various multicomponents models 
# ============================================================================= 
""" Test module 
- It tests various multicomponents models 
"""
# ============================================================================= 
from   builtins                 import range
from   ostap.core.core          import cpp, VE, dsID, rooSilent 
from   ostap.fitting.background import make_bkg
from   ostap.core.meta_info     import root_info 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait, batch_env  
import ostap.fitting.models     as     Models 
import ostap.fitting.roofit 
import ROOT, random, time, sys, datetime  
# ============================================================================= 
from   ostap.logger.logger      import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'test_fitting_in_range_2d' )
else                       : logger = getLogger ( __name__  )
# ============================================================================= 
## set batch from environment 
batch_env ( logger )
# =============================================================================

# ============================================================================= 
def test_fitting_in_range_2d () :
    
    logger = getLogger ( 'test_fitting_in_range_2d' ) 
    ## make simple test mass 
    m_x     = ROOT.RooRealVar ( 'm_x' , 'Some test mass(X)' , 0 , 5 )
    m_y     = ROOT.RooRealVar ( 'm_y' , 'Some test mass(Y)' , 6 , 10 )
    
    ## book very simple data set
    varset  = ROOT.RooArgSet  ( m_x , m_y)
    dataset = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset )  
    
    m1 = VE ( 3 , 0.10**2 )
    m2 = VE ( 7 , 0.10**2 )
    
    ## fill it with three gausissians, 5k events each
    N_ss = 5000
    N_sb = 1000
    N_bs = 500
    N_bb = 100
    
    random.seed(0)
    
    ## S x S
    for i in range (N_ss ) :
        m_x.value = m1.gauss()
        m_y.value = m2.gauss() 
        dataset.add ( varset  )
		
    ## S x B 
    for i in range ( N_sb ) :
        m_x.value = m1.gauss()
        m_y.value = random.uniform ( *m_y.minmax() )
        dataset.add ( varset  )
	
    ## B x S 
    for i in range ( N_bs ) :
        m_x.value = random.uniform ( *m_x.minmax() )
        m_y.value = m2.gauss()
        dataset.add ( varset  )
	
    ## B x B 
    for i in range ( N_bb ) :
        m_x.value  = random.uniform ( *m_x.minmax() )
        m_y.value  = random.uniform ( *m_y.minmax() )
        dataset.add ( varset  )
	
    logger.info ('Dataset:\n%s' % dataset.table ( prefix = '# ') )   
    
    ## various fit components
    signal_x1 = Models.Gauss_pdf ( 'G1x'  , xvar = m_x  , mean = m1.value() , sigma = m1.error() )  
    signal_y1 = Models.Gauss_pdf ( name='G1y'  , xvar = m_y  , mean = m2.value() , sigma = m2.error()  ) 
    
    bkg_x= make_bkg ( -1 , 'Bx' , m_x )
    bkg_y= make_bkg ( -1 , name= 'By' , xvar =m_y )
    
    ## build fit model 
    model = Models.Fit2D ( name        = 'fit_comp', 
			   signal_x    = signal_x1, 
			   signal_y    = signal_y1,
			   bkg_1x      = bkg_x ,
			   bkg_1y      = bkg_y )
    model.SS = N_ss
    model.SB = N_sb
    model.BS = N_bs
    model.BB = N_bb
    
    r = model.fitTo ( dataset , silent = True )
    r = model.fitTo ( dataset , silent = True )
    
    dataset.m_y.setRange ( 'fit' , 8,10. )
    model.yvar.setRange  ( 'fit' , 8,10. )
    
    with use_canvas ( 'test_fitting_in_range_2d' ) : 		
        with wait ( 2 ) : model.draw1(dataset,nbins=200,in_range=(6,8))
        with wait ( 2 ) : model.draw1(dataset,nbins=200, in_range='fit')
	
    dataset.m_x.setRange ( 'fit2' , 0,2.5 )
    model.xvar.setRange  ( 'fit2' , 0,2.5 )
	
    with use_canvas ( 'test_fitting_in_range_2d' ) :
        with wait ( 2 ) : model.draw2(dataset,nbins=200, in_range=(2.5,5))
        with wait ( 2 ) : model.draw2(dataset,nbins=200, in_range='fit2')
	
# =============================================================================
if '__main__' == __name__ :

    test_fitting_in_range_2d ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
		

