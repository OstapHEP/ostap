#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file test_histo_compare.py
# Test module for ostap/histos/compare.py
# - It tests comparision of 1D-histograms 
# ============================================================================= 
"""Test module for ostap/histos/compare.py
- It tests comparision of 1D-histograms
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import 
# ============================================================================= 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_histo_compare' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( 'Test for 1D-histogram compare')
# =============================================================================
from   ostap.math.ve        import VE 
from   ostap.core.core      import hID 
from   ostap.histos.histos  import h1_axis 
import ostap.histos.compare 
#
## histos for gaussian distributions
# 

h1g   = ROOT.TH1D ( hID() , '' ,  40 , -5 , 5 ) ; h1g.Sumw2() 
h2g   = ROOT.TH1D ( hID() , '' ,  40 , -5 , 5 ) ; h2g.Sumw2() 
h3g   = ROOT.TH1D ( hID() , '' ,  20 , -5 , 5 ) ; h3g.Sumw2() 

bins  = [ -5 ]
## 
random.seed(10) 
for i in range(0, 15 ) : bins.append ( random.uniform ( -5 , 5 ) )
bins += [  5 ]
bins.sort()

h4g   = h1_axis ( bins ) 

#
## histos for uniform distributions
# 
h1u = h1g.clone()
h2u = h2g.clone()
h3u = h3g.clone()
h4u = h4g.clone()

#
## histos for exponential distributions
# 
h1e = h1g.clone()
h2e = h2g.clone()
h3e = h3g.clone()
h4e = h4g.clone()

## get the value 
v   = VE(1,1.75**2)

random.seed(10) 
for i in range ( 0, 50000 ) :

    g1 = v.gauss()
    g2 = v.gauss()
    g3 = g2          ## the same as g2 
    g4 = v.gauss()
    
    h1g.Fill ( g1 ) 
    h2g.Fill ( g2 ) 
    h3g.Fill ( g3 ) 
    h4g.Fill ( g4 ) 
    
    u1 = random.uniform ( -5 , 5 )
    u2 = random.uniform ( -5 , 5 )
    u3 = u2          #3 the same as u2 
    u4 = random.uniform ( -5 , 5 )
    
    h1u.Fill ( u1 ) 
    h2u.Fill ( u2 ) 
    h3u.Fill ( u3 ) 
    h4u.Fill ( u4 ) 

    e1 = -1 * random.expovariate( -0.5 )  -5 
    e2 = -1 * random.expovariate( -0.5 )  -5 
    e3 = e2  ## the same as e2 
    e4 = -1 * random.expovariate( -0.5 )  -5 
    if not -5 < e1 < 5 : continue 
    if not -5 < e2 < 5 : continue 
    if not -5 < e3 < 5 : continue 
    if not -5 < e4 < 5 : continue 

    h1e.Fill ( e1 ) 
    h2e.Fill ( e2 ) 
    h3e.Fill ( e3 ) 
    h4e.Fill ( e4 ) 

h5g = h4g.rescale_bins(1)
h5u = h4u.rescale_bins(1)
h5e = h4e.rescale_bins(1)

## compare two histograms 
def compare ( h1 , h2 , rescale = False ) :
    
    r1  = h1.cmp_fit       ( h2 , opts = 'WL0Q' , rescale = rescale )
    if r1 : logger.info    ( 'h1 vs h2 : fit probability is %.5f%% ' % ( r1.Prob()*100 ) )
    else  : logger.warning ( 'h1 vs h2 : fit problems ')
    
    r2  = h2.cmp_fit       ( h1 , opts = 'WL0Q' , rescale = rescale )
    if r2 : logger.info    ( 'h2 vs h1 : fit probability is %.5f%% ' % ( r2.Prob()*100 ) )
    else  : logger.warning ( 'h2 vs h1 : fit problems ')

    ct  = h1.cmp_cos      ( h2 , rescale = rescale ) 
    logger.info           ( 'h1 vs h2 : cos(theta)      is %s ' % ct  )
    
    dd1 = h1.cmp_dist     ( h2 , rescale = rescale ) 
    logger.info           ( 'h1 vs h2 : distance        is %s ' % dd1 )
    
    dd2 = h1.cmp_dist2    ( h2 , rescale = rescale ) 
    logger.info           ( 'h1 vs h2 : distance2       is %s ' % dd2 )
    
    h1.cmp_prnt ( h2 , 'h1 vs h2 ' , rescale = rescale ) 

# =============================================================================
## compare gaussians 
def test_compare_gaussians() : 
    logger.info( 'Compare gaussians    (1) and (2)' )
    compare ( h1g , h2g )
    logger.info( 'Compare gaussians    (1) and (3)' )
    compare ( h1g , h3g )
    logger.info( 'Compare gaussians    (1) and (4)' )
    compare ( h1g , h4g )
    logger.info( 'Compare gaussians    (1) and (4) with rescale' )
    compare ( h1g , h4g , True )
    logger.info( 'Compare gaussians    (1) and (5)' )
    compare ( h1g , h5g )
    logger.info( 'Compare gaussians    (2) and (3) : should be the same!' )
    compare ( h2g , h3g )
    logger.info( 'Compare gaussians    (2) and (4)' )
    compare ( h2g , h4g )
    logger.info( 'Compare gaussians    (2) and (4) with rescale' )
    compare ( h2g , h4g , True )
    logger.info( 'Compare gaussians    (2) and (5)' )
    compare ( h2g , h5g )
    logger.info( 'Compare gaussians    (3) and (4)' ) 
    compare ( h3g , h4g )
    logger.info( 'Compare gaussians    (3) and (4) with rescale' ) 
    compare ( h3g , h4g , True )
    logger.info( 'Compare gaussians    (3) and (5)' ) 
    compare ( h3g , h5g )
    logger.info( 'Compare gaussians    (4) and (5)' ) 
    compare ( h4g , h5g )

def test_compare_uniforms () :
    logger.info( 'Compare uniforms     (1) and (2)' )
    compare ( h1u , h2u )
    logger.info( 'Compare uniforms     (1) and (3)' )
    compare ( h1u , h3u )
    logger.info( 'Compare uniforms     (1) and (4)' )
    compare ( h1u , h4u )
    logger.info( 'Compare uniforms     (1) and (4) with rescale' )
    compare ( h1u , h4u , True )
    logger.info( 'Compare uniforms     (1) and (5)' )
    compare ( h1u , h5u )
    logger.info( 'Compare uniforms     (2) and (3) : should be the same!' )
    compare ( h2u , h3u )
    logger.info( 'Compare uniforms     (2) and (4)' )
    compare ( h2u , h4u )
    logger.info( 'Compare uniforms     (2) and (4) with rescale' )
    compare ( h2u , h4u , True )
    logger.info( 'Compare uniforms     (2) and (5)' )
    compare ( h2u , h4u )
    logger.info( 'Compare uniforms     (3) and (4)' ) 
    compare ( h3u , h4u )
    logger.info( 'Compare uniforms     (3) and (4) with rescale;' ) 
    compare ( h3u , h4u , True )
    logger.info( 'Compare uniforms     (3) and (5)' ) 
    compare ( h3u , h5u )
    logger.info( 'Compare uniforms     (4) and (5)' ) 
    compare ( h4u , h5u )

def test_compare_exponentials () :
    logger.info( 'Compare exponentials (1) and (2)' )
    compare ( h1e , h2e )
    logger.info( 'Compare exponentials (1) and (3)' )
    compare ( h1e , h3e )
    logger.info( 'Compare exponentials (1) and (4)' )
    compare ( h1e , h4e )
    logger.info( 'Compare exponentials (1) and (4) with rescale' )
    compare ( h1e , h4e , True )
    logger.info( 'Compare exponentials (1) and (5)' )
    compare ( h1e , h5e )
    logger.info( 'Compare exponentials (2) and (3) : should be the same!' )
    compare ( h2e , h3e )
    logger.info( 'Compare exponentials (2) and (4)' )
    compare ( h2e , h4e )
    logger.info( 'Compare exponentials (2) and (4) with rescale' )
    compare ( h2e , h4e , True )
    logger.info( 'Compare exponentials (2) and (5)' )
    compare ( h2e , h5e )
    logger.info( 'Compare exponentials (3) and (4)' ) 
    compare ( h3e , h4e )
    logger.info( 'Compare exponentials (3) and (4) with rescale' ) 
    compare ( h3e , h4e , True )
    logger.info( 'Compare exponentials (3) and (5)' ) 
    compare ( h3e , h5e )
    logger.info( 'Compare exponentials (4) and (5)' ) 
    compare ( h4e , h5e )

def test_compare_gauss_vs_uniform() :     
    _ig = 0 
    for ig in ( h1g , h2g , h3g , h4g , h5g ) :
        _ig += 1
        _iu  = 0 
        for iu in ( h1u , h2u , h3u , h4u , h5u ) :
            _iu += 1 
            logger.info( 'Compare gaussian  (%d) and uniform     (%d)' % ( _ig , _iu ) )
            compare ( ig , iu )
            
def test_compare_gauss_vs_exponent () :     
    _ig = 0 
    for ig in ( h1g , h2g , h3g , h4g , h5g ) :
        _ig += 1
        _ie  = 0 
        for ie in ( h1e , h2e , h3e , h4e , h5e ) :
            _ie += 1 
            logger.info( 'Compare gaussian  (%d) and exponent    (%d)' % ( _ig , _ie ) )
            compare ( ig , ie )

def test_compare_uniform_vs_exponent () :     
    _iu = 0 
    for iu in ( h1u , h2u , h3u , h4u , h5u ) :
        _iu += 1
        _ie  = 0 
        for ie in ( h1e , h2e , h3e , h4e , h5e ) :
            _ie += 1 
            logger.info( 'Compare uniform   (%d) and exponent    (%d)' % ( _iu , _ie ) )
            compare ( iu , ie )
            
# =============================================================================
if '__main__' == __name__ :
    
    test_compare_gaussians           ()
    test_compare_uniforms            ()
    test_compare_exponentials        ()
    test_compare_gauss_vs_uniform    ()
    test_compare_gauss_vs_exponent   ()
    test_compare_uniform_vs_exponent ()
    
    pass

# =============================================================================
# The END 
# =============================================================================

