#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/stats/tests/test_stats_blue.py 
# Test module BLUE Best Linear Unbiased Estimator
# @see ostap/stats/combine.py
# Copyright (c) Ostap developpers.
# ============================================================================= 
""" Test module BLUE Best Linear Unbiased Estimator
- see ostap/stat/combine.py.
"""
# =============================================================================
import ROOT
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'tests_stats_blue' )
else                       : logger = getLogger ( __name__           )
# ==============================================================================
from   ostap.stats.combine import Combine, Ostap, VE 

def test_stats_blue1 () :
    
    m1 =  VE( 5366.779, 0.069 **2 )
    m2 =  VE( 5367.122, 0.089 **2 ) 
    
    ## fit model systematics, assumed to be uncorrelated 
    syst1 = Ostap.Math.SymMatrix2x2()
    syst1 [0,0] = 0.134**2
    syst1 [1,1] = 0.091**2
    
    ## momentum scaling, 100% correlated  
    syst2 = Ostap.Math.SymMatrix2x2()
    syst2 [0,0] = 0.315**2
    syst2 [1,1] = 0.132**2
    syst2 [0,1] = 0.132*0.315
    
    ## energy loss, 100% correlated  
    syst3 = Ostap.Math.SymMatrix2x2()
    syst3 [0,0] = 0.030**2
    syst3 [1,1] = 0.015**2
    syst3 [0,1] = 0.030*0.015
    
    ## kaon mass, 100% correlated  
    syst4 = Ostap.Math.SymMatrix2x2()
    syst4 [0,0] = 0.026**2
    syst4 [1,1] = 0.026**2
    syst4 [0,1] = 0.026*0.026
    
    cmb = Combine( [m1,m2] , syst1 + syst2 + syst3 + syst4 )
    
    r           = cmb.result        
    stat , syst = cmb.errors 
    
    for i,c  in enumerate ( ( m1 , m2 ) ) :
        logger.info ( 'Component %d: %.3f +/- %.3f ' % ( i + 1 , c.value() , c.error() ) )

    logger.info ( 'Combined   : %.3f +/- %.3f = (%.3f +/- %.3f +/- %.3f) ' % ( r.value() , r.error()  , r.value() , stat, syst ) ) 
    logger.info ( 'Weights    : %s' %  cmb.weights ) 

def test_stats_blue2 () :
    
    m1 =  VE ( 5366.779 , 0.069**2 ) ## Bs -> J/psi K*K*
    m2 =  VE ( 5367.122 , 0.089**2 ) ## Bs -> ( psi' -> J/psi pi+ pi-) (phi->K+ K-)
    m3 =  VE ( 5366.83  , 0.25 **2 ) ## Bs -> chi_c1 -> ( J/psi gamma ) K+ K- 
    m4 =  VE ( 5367.08  , 0.38 **2 ) 
    m5 =  VE ( 5366.90  , 0.22 **2 ) 
    m6 =  VE ( 5366.77  , 0.24 **2 )
    
    COV2 = Ostap.Math.SymMatrix6x6
    
    ## momentum scaling for 2010, uncorrelated 
    syst0  = COV2()
    syst0 [4,4] = 0.28**2 
    
    ## fit model systematics, assumed to be uncorrelated 
    syst1 = COV2 ()
    syst1 [0,0] = 0.134**2
    syst1 [1,1] = 0.091**2
    syst1 [2,2] = (0.01+0.00)**2
    syst1 [3,3] = (0.02+0.00)**2
    syst1 [4,4] = (0.01+0.02)**2
    syst1 [5,5] = (0.02+0.00)**2    
    
    ## momentum scaling, 100% correlated  
    syst2 = COV2()
    syst2 [0,0] = 0.315**2
    syst2 [1,1] = 0.132**2
    syst2 [2,2] = 0.26 **2
    syst2 [3,3] = 0.12 **2
    syst2 [4,4] = 0.   **2
    syst2 [5,5] = 0.124**2
    
    
    ## energy loss, 100% correlated  
    syst3 = COV2()
    syst3 [0,0] = 0.030**2
    syst3 [1,1] = 0.015**2
    syst3 [2,2] = 0.020**2
    syst3 [3,3] = 0.060**2
    syst3 [4,4] = 0.030**2
    syst3 [5,5] = 0.030**2
    
    
    ## kaon mass, 100% correlated  
    syst4 = COV2()
    syst4 [0,0] = 0.026**2
    syst4 [1,1] = 0.026**2
    syst4 [2,2] = 0.020**2
    syst4 [3,3] = 0.060**2
    syst4 [4,4] = 0.020**2
    syst4 [5,5] = 0.000**2
    
    ## they are 100% correlated
    for CV in ( syst2 , syst3 , syst4 ) :
        for i in range(6) :
            cii = CV (i,i)
            for j in  range ( i + 1 , 6 ) :
                cjj = CV (j,j)
                CV  [ i , j ] = ( cii * cjj )**0.5 
                
                
    cmb = Combine( [m1,m2,m3,m4,m5,m6] , syst0 + syst1 + syst2 + syst3 + syst4 )
    
    r          = cmb.result        
    stat, syst = cmb.errors 
    
    
    for i,c  in enumerate ( ( m1 , m2 , m3 , m4 , m5 , m6 ) ) :
        logger.info ( 'Component %d: %.3f +/- %.3f ' % ( i + 1 , c.value() , c.error() ) )
        
    logger.info ( 'Combined   : %.3f +/- %.3f = (%.3f +/- %.3f +/- %.3f) ' % ( r.value() , r.error()  , r.value() , stat, syst ) ) 
    logger.info ( 'Weights    : %s' %  cmb.weights ) 
        

# =============================================================================
if '__main__' == __name__ :

    test_stats_blue1()
    test_stats_blue2()
# =============================================================================
##                                                                      The END 
# ============================================================================= 

