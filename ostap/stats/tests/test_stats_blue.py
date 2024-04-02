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
from   ostap.stats.combine import Combine, Ostap, VE, covMatrix  
import ROOT
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'tests_stats_blue' )
else                       : logger = getLogger ( __name__           )
# ==============================================================================

# =============================================================================
def test_stats_blue1 () :

    ## values with their stat uncertainties
    
    m1 =  VE ( 5366.779, 0.069 **2 )
    m2 =  VE ( 5367.122, 0.089 **2 ) 
    
    ## fit model systematics, assumed to be uncorrelated, add also psi' mass  
    syst1 = covMatrix ( False , 0.134 , ( 0.091 ** 2 + 0.01 **2 ) ** 0.5  )
    
    ## momentum scaling systematic, 100% correlated  
    syst2 = covMatrix ( True , 0.315 , 0.132 ) 
    
    ## energy loss systematic, 100% correlated  
    syst3 = covMatrix ( True , 0.030 , 0.015 ) 
    
    ## kaon mass systematic,   100% correlated  
    syst4 = covMatrix ( True , 0.029 , 0.029 ) 

    cmb = Combine( [m1,m2] , syst1 + syst2 + syst3 + syst4 )
    
    r           = cmb.result        
    stat , syst = cmb.errors 
    
    for i,c  in enumerate ( ( m1 , m2 ) ) :
        logger.info ( 'Component %d: %.3f +/- %.3f ' % ( i + 1 , c.value() , c.error() ) )

    logger.info ( 'Combined   : %.3f +/- %.3f = (%.3f +/- %.3f +/- %.3f) ' % ( r.value() , r.error()  , r.value() , stat, syst ) ) 
    logger.info ( 'Weights    : %s'       %  cmb.weights ) 
    logger.info ( 'chi2/ndf   : %.2f/%d'  %  ( cmb.chi2 , cmb.D ) ) 
    logger.info ( 'p-value    : %.2f[%%]' %  ( cmb.pvalue * 100 ) ) 

# =============================================================================
def test_stats_blue2 () :
    
    m0 =  VE ( 5366.779 , 0.069**2 ) ## Bs -> J/psi K*K*
    m1 =  VE ( 5367.122 , 0.089**2 ) ## Bs -> ( psi' -> J/psi pi+ pi-) (phi->K+ K-)    
    m2 =  VE ( 5366.83  , 0.25 **2 ) ## Bs -> chi_c1 -> ( J/psi gamma ) K+ K- 
    m3 =  VE ( 5367.08  , 0.38 **2 ) ## Bs -> J/psi  phi phi 
    m4 =  VE ( 5366.90  , 0.22 **2 ) ## Bs -> J/psi phi
    m5 =  VE ( 5366.85  , 0.19 **2 ) ## Bs -> J/psi p pbar
    
    COV2 = Ostap.Math.SymMatrix(6)
    
    ## charmonuium masses, uncorrelated  
    syst   = covMatrix ( False , 0.0 , 0.01 , 0.04 , 0 , 0 , 0 ) 

    ## momentum scaling for 2010, uncorrelated 
    syst0  = covMatrix ( False , 0    , 0    , 0 , 0 , 0.22 , 0 )
    
    ## fit model systematics, assumed to be uncorrelated 
    syst1  = covMatrix ( False , 0.134 , 0.091 , 0.01 + 0.0 , 0.02 + 0.0 , 0.02 + 0.01 , 0.02 + 0.0 ) 
    
    ## momentum scaling, 100% correlated  
    syst2 = covMatrix  ( True , 0.315 , 0.132 , 0.26 , 0.12 , 0 , 0.124 ) 
    
    ## energy loss, 100% correlated  
    syst3 = covMatrix  ( True , 0.030 , 0.015 , 0.020 , 0.060 , 0.030 , 0.030 )
    
    ## kaon mass, 100% correlated  
    syst4 = covMatrix  ( True , 0.029 , 0.029 , 0.020 , 0.060 , 0.020 , 0.000 )

    cmb   = Combine( [m0,m1,m2,m3,m4,m5] , syst + syst0 + syst1 + syst2 + syst3 + syst4 )
    
    r          = cmb.result        
    stat, syst = cmb.errors 
    
    
    for i,c  in enumerate ( ( m0 , m1 , m2 , m3 , m4 , m5 ) ) :
        logger.info ( 'Component %d: %.3f +/- %.3f ' % ( i + 1 , c.value() , c.error() ) )
        
    logger.info ( 'Combined   : %.3f +/- %.3f = (%.3f +/- %.3f +/- %.3f) ' % ( r.value() , r.error()  , r.value() , stat, syst ) ) 
    logger.info ( 'Weights    : %s'   %  cmb.weights ) 
    logger.info ( 'chi2/ndf   : %.2f/%d'  %  ( cmb.chi2 , cmb.D ) ) 
    logger.info ( 'p-value    : %.2f[%%]' %  ( cmb.pvalue * 100 ) ) 
        

# =============================================================================
if '__main__' == __name__ :

    test_stats_blue1 ()
    test_stats_blue2 ()
    
# =============================================================================
##                                                                      The END 
# ============================================================================= 

