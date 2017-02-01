#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file test_histo_interpolation.py
# Test module for ostap/histos/histos.py
# - It tests interpolation for 1,2&3D histograms 
# ============================================================================= 
""" Test module for ostap/histos/histos.py
- It tests interpolation for 1,2&3D histograms 
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
    logger = getLogger ( 'test_histo_interpolation' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( 'Test for 1,2&3D-histogram interpolation')
# =============================================================================

import ROOT, random, ostap.histos.histos  
from   ostap.core.core      import hID, SE 

h1 = ROOT.TH1D ( hID() , '', 5 , 0 , 1 )
h2 = ROOT.TH2D ( hID() , '', 5 , 0 , 1 , 5 , 0 , 1 ) 
h3 = ROOT.TH3D ( hID() , '', 5 , 0 , 1 , 5 , 0 , 1 , 5 , 0 , 1 ) 

fun1 = lambda x         : 0.1+x**3
fun2 = lambda x, y      : 0.1+x**3+y**2+x*y
fun3 = lambda x, y , z  : 0.1+x**3+y**2+x*y*z+z*x 

h1 += fun1
h2 += fun2
h3 += fun3

from collections import defaultdict 

##  test interpolation for 1D-historgams 
def test_1D() :
    
    N1 = 2000
    cnts_1D = defaultdict(SE) 
    for i in range(N1) :   
        x  = random.uniform(0,1)
        vf = fun1  ( x )
        for t in ( 0 , 1 , 2 , 3 ) :
            for edges in ( False , True ) :
                for extrapolate in ( False , ) :
                    key = t , edges , extrapolate
                    vh  = h1 (  x , interpolate = t , edges = edges  , extrapolate = extrapolate ) 
                    cnts_1D[key] += (vh.value()/vf-1.0)*100
                    
    logger.info('1D interpolation (wide bins): ') 
    keys = cnts_1D.keys()
    keys.sort()
    for k in keys :
        v = cnts_1D[k]
        b = v.mean()
        r = v.rms ()
        logger.info ( 'Key:%-20s bias:%15s%% RMS:%.2f%%' % ( k , b.toString("(%.2f+-%.2f)") , r ) ) 
        

# =============================================================================
##  test interpolation for 2D-historgams 
def test_2D() :

    N2 = 10
    cnts_2D = defaultdict(SE) 
    for i in range(N2) :   
        x  = random.uniform(0,1)
        for j in range(N2) :
            y  = random.uniform(0,1)        
            vf = fun2  ( x , y )        
            for tx in ( 0 , 1 , 2 , 3 ) :
                for ty in ( 0 , 1 , 2 , 3 ) :                
                    for edges in ( False , True ) :
                        for extrapolate in ( False , ) :
                            key = tx, ty , edges , extrapolate                        
                            vh  = h2 ( x , y , interpolate = (tx , ty)  ,
                                       edges = edges  , extrapolate = extrapolate  )
                            cnts_2D[key] += (vh.value()/vf-1.0)*100
                            
                            
    logger.info('2D interpolation (wide bnis): ') 
    keys = cnts_2D.keys()
    keys.sort()
    for k in keys :
        v = cnts_2D[k]
        b = v.mean()
        r = v.rms ()
        logger.info ( 'Key:%-20s bias:%15s%% RMS:%.2f%%' % ( k , b.toString("(%.2f+-%.2f)") , r ) ) 
        
# =============================================================================
##  test interpolation for 3D-historgams 
def test_3D() :

    cnts_3D = defaultdict(SE)
    
    N3 = 10
    for i in range(N3) :   
        x  = random.uniform(0,1)
        for j in range(N3) :
            y  = random.uniform(0,1)
            for k in range(N3) :
                z  = random.uniform(0,1)
                vf = fun3  ( x , y , z )        
                for tx in ( 0 , 1 , 2 , 3 ) :
                    for ty in ( 0 , 1 , 2 , 3 ) :
                        for tz in ( 0 , 1 , 2 , 3 ) :                        
                            for edges in ( False , True ) :
                                for extrapolate in ( False , ) :
                                    key = tx, ty , tz,  edges , extrapolate                              
                                    vh  = h3 ( x , y , z , interpolate= ( tx , ty , tz )  ,
                                               edges  = edges , extrapolate = extrapolate ) 
                                    cnts_3D[key] += (vh.value()/vf-1.0)*100
                                    
                
    logger.info('3D interpolation (wide bins): ') 
    keys = cnts_3D.keys()
    keys.sort()
    for k in keys :
        v = cnts_3D[k]
        b = v.mean()
        r = v.rms ()
        logger.info ( 'Key:%-20s bias:%15s%% RMS:%.2f%%' % ( k , b.toString("(%.2f+-%.2f)") , r ) ) 
        
        
# =============================================================================
if '__main__' == __name__ :

    test_1D() ## test interpolaiton for 1D-histograms
    test_2D() ## test interpolaiton for 2D-histograms
    test_3D() ## test interpolaiton for 3D-histograms
    
# =============================================================================
# The END 
# =============================================================================
