#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_in_range_3d.py
# Test module 
# - It tests various multicomponents models 
# ============================================================================= 
""" Test module 
- It tests various multicomponents models 
"""
# ============================================================================= 
from   ostap.core.pyrouts       import *
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import cpp, VE, dsID
from   ostap.logger.utils       import rooSilent
from   builtins                 import range
from   ostap.fitting.background import make_bkg 
from   ostap.core.meta_info     import root_info 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
import ROOT, sys, random, time
# ============================================================================= 
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'test_fitting_in_range_3d' )
else                       : logger = getLogger ( __name__  )
# ============================================================================= 

# ============================================================================= 
def test_fitting_in_range_3d () :
        
	logger = getLogger ( 'test_fitting_in_range_3d' )
	
        ## make simple test mass 
	m_x     = ROOT.RooRealVar ( 'm_x' , 'Some test mass(X)' , 0 , 5 )
	m_y     = ROOT.RooRealVar ( 'm_y' , 'Some test mass(Y)' , 6 , 10 )
	m_z     = ROOT.RooRealVar ( 'm_z' , 'Some test mass(z)' , 10 , 15 )
	
        ## book very simple data set
	varset  = ROOT.RooArgSet  ( m_x , m_y,m_z )
	dataset = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset )  
	
	
	m1 = VE (  3 , 0.10**2 )
	m2 = VE (  7 , 0.10**2 )
	m3 = VE ( 12 , 0.10**2 )
	
	N_sss = 5000
	N_ssb = 5000
	N_sbs = 5000
	N_sbb = 1000
	
	N_bss = 500
	N_bsb =  100
	N_bbs =  100
	N_bbb = 250

	random.seed(0)

	## S x S x S 
	for i in range ( N_sss ) : 
		m_x.value = m1.gauss() 
		m_y.value = m2.gauss() 
		m_z.value = m3.gauss() 
		dataset.add ( varset  )
		
	## S x S x B 
	for i in range ( N_ssb ):
		m_x.value = m1.gauss() 
		m_y.value = m2.gauss() 
		m_z.value = random.uniform ( *m_z.minmax() )
		dataset.add ( varset  )
		
	## S x B x S 
	for i in range ( N_sbs ) : 
		m_x.value = m1.gauss() 
		m_y.value = random.uniform ( *m_y.minmax() )  
		m_z.value = m3.gauss() 		
		dataset.add ( varset  )
		
	## B x S x S
	for i in range ( N_bss ) : 
		m_x.value = random.uniform ( *m_x.minmax() ) 
		m_y.value = m2.gauss() 
		m_z.value = m3.gauss() 
		dataset.add ( varset  )
		
		
	## S x B x B
	for i in range ( N_sbb ) : 
		m_x.value = m1.gauss() 
		m_y.value = random.uniform ( *m_y.minmax() )  
		m_z.value = random.uniform ( *m_z.minmax() )
		dataset.add ( varset  )
		
	## B x S x B   
	for i in range ( N_bsb ) : 
		m_x.value = random.uniform ( *m_x.minmax() ) 
		m_y.value = m2.gauss() 
		m_z.value = random.uniform ( *m_z.minmax() )
		dataset.add ( varset  )


	## B x B x S 
	for i in range ( N_bbs ) : 
		m_x.value = random.uniform ( *m_x.minmax() ) 
		m_y.value = random.uniform ( *m_y.minmax() )  
		m_z.value = m3.gauss() 
		dataset.add ( varset  )
		
	## B x B x B 
	for i in range ( N_bbb ) : 
		m_x.value = random.uniform ( *m_x.minmax() ) 
		m_y.value = random.uniform ( *m_y.minmax() )  
		m_z.value = random.uniform ( *m_z.minmax() )
		dataset.add ( varset  )
		
	logger.info ('Dataset:\n%s' % dataset.table ( prefix = '# ' ) )  

	
	signal_x1 = Models.Gauss_pdf ( 'G1x'  , xvar = m_x  , mean = m1.value() , sigma = m1.error() )  
	signal_y1 = Models.Gauss_pdf ( name='G1y'  , xvar = m_y  , mean = m2.value() , sigma = m2.error()  ) 
	signal_z1 = Models.Gauss_pdf ( name='G1z'  , xvar = m_z  , mean = m3.value() , sigma = m3.error()  )
	
	bkg_x= make_bkg ( -1 , 'Bx' , m_x )
	bkg_y= make_bkg ( -1 , name= 'By' , xvar =m_y )
	bkg_z= make_bkg ( -1 ,name='Bz' , xvar =m_z )
	
	
	model = Models.Fit3D (
		name        = 'fit3'    , 
		signal_x    = signal_x1 , 
		signal_y    = signal_y1 ,
		signal_z    = signal_z1 ,
		bkg_1x  = bkg_x ,
		bkg_1y  = bkg_y ,
		bkg_1z  = bkg_z ,
		)

	model.SSS = N_sss
	model.SSB = N_ssb
	model.SBS = N_sbs
	model.BSS = N_bss
	model.SBB = N_sbb
	model.BSB = N_bsb
	model.BBS = N_bbs
	model.BBB = N_bbb

	r = model.fitTo ( dataset , silent = True )
	r = model.fitTo ( dataset , silent = True )

	t = 1.0
        
	with use_canvas ( 'test_fitting_in_range_3d (no-ranges) ' ) :
		with wait ( t ) : model.draw1 (dataset,nbins=200)
		with wait ( t ) : model.draw2 (dataset,nbins=200)
		with wait ( t ) : model.draw3 (dataset,nbins=200)

	dataset.m_y.setRange ( 'fit' , 6,8. )
	model.yvar.setRange ( 'fit' , 6,8. )
                
	
	with use_canvas ( 'test_fitting_in_range_3d/1' ) :
		with wait ( t ) : model.draw1(dataset,nbins=200, in_range3=(11,12),in_range2=(8,10))
		with wait ( t ) : model.draw1(dataset,nbins=200, in_range3=(11,12),in_range2='fit')
		with wait ( t ) : model.draw1(dataset,nbins=200, in_range3=(11,12))
		with wait ( t ) : model.draw1(dataset,nbins=200, in_range2='fit')
	
	dataset.m_x.setRange ( 'fit2' , 2.5,3. )
	model.xvar.setRange ( 'fit2' , 2.5,3. )
	
	with use_canvas ( 'test_fitting_in_range_3d/2' ) :		
		with wait ( t ) : model.draw2(dataset,nbins=200, in_range3=(11,12),in_range1=(0,3))
		with wait ( t ) : model.draw2(dataset,nbins=200, in_range3=(11,12),in_range1='fit2')
		with wait ( t ) : model.draw2(dataset,nbins=200, in_range3=(11,12))
		with wait ( t ) : model.draw2(dataset,nbins=200, in_range1='fit2')
		
	dataset.m_x.setRange ( 'fit3' , 2.5,3. )
	model.xvar.setRange ( 'fit3' , 2.5,3. )
	
	with use_canvas ( 'test_fitting_in_range_3d/3' ) :		
		with wait ( t ) : model.draw3(dataset,nbins=200, in_range2=(6,8),in_range1=(0,3))
		with wait ( t ) : model.draw3(dataset,nbins=200, in_range2=(6,8),in_range1='fit3')
		with wait ( t ) : model.draw3(dataset,nbins=200, in_range2=(6,8))
		with wait ( t ) : model.draw3(dataset,nbins=200, in_range1='fit3')
				
# =============================================================================
if '__main__' == __name__ :

	test_fitting_in_range_3d ()

# =============================================================================
##                                                                      The END 
# =============================================================================
		
