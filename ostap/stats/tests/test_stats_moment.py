#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/stats/tests/test_stats_moment.py
# Test module for ostap/stat/moment.py.
# Copyright (c) Ostap developpers.
# ============================================================================= 
""" Test module for ostap/stat/moment.py.
"""
# =============================================================================
import ostap.stats.moment
from   ostap.core.core import Ostap
import ROOT,random
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_stats_moment' )
else                       : logger = getLogger ( __name__            )
# ============================================================================= 
root_version = ROOT.gROOT.GetVersionInt() 
# ============================================================================= 
def test_moment1() :

    
    m0  = Ostap.Math.Moment_( 0)()
    m1  = Ostap.Math.Moment_( 1)()
    m2  = Ostap.Math.Moment_( 2)()
    m3  = Ostap.Math.Moment_( 3)()
    m4  = Ostap.Math.Moment_( 4)()
    m5  = Ostap.Math.Moment_( 5)()
    m10 = Ostap.Math.Moment_(10)()
    m20 = Ostap.Math.Moment_(20)()

    m = m0, m1, m2 , m3, m4, m5 , m10, m20
    
    for i in range ( 10000 ) :
        v = random.gauss ( 0 , 1 ) 
        for i in m : i += v

    if root_version < 61800 :
        logger.warning ( "It does not work for ancient ROOT versions")
        return

    for i in m :
        t = 'Moment counter-%2d' % i.order 
        logger.info ( "%s\n%s" % ( t , i.table ( title = t , prefix = '# ' ) ) )

    for i in m :
        t = 'Moment counter-%2d' % i.order 
        logger                  .info (' %s : size     %s ' % ( t , i.size () ) )
        if 1 <= i.order : logger.info (' %s : mean     %s ' % ( t , i.mean         () ) )
        if 2 <= i.order : logger.info (' %s : rms      %s ' % ( t , i.rms          () ) )
        if 2 <= i.order : logger.info (' %s : variance %s ' % ( t , i.variance     () ) )
        if 2 <= i.order : logger.info (' %s : u2nd     %s ' % ( t , i.unbiased_2nd () ) )
        if 3 <= i.order : logger.info (' %s : u3rd     %s ' % ( t , i.unbiased_3rd () ) )
        if 3 <= i.order : logger.info (' %s : u3rd     %s ' % ( t , i.unbiased_3rd () ) )
        if 4 <= i.order : logger.info (' %s : u4th     %s ' % ( t , i.unbiased_4th () ) )
        if 5 <= i.order : logger.info (' %s : u5th     %s ' % ( t , i.unbiased_5th () ) )
        if 3 <= i.order : logger.info (' %s : skewness %s ' % ( t , i.skewness     () ) )
        if 4 <= i.order : logger.info (' %s : kurtosis %s ' % ( t , i.kurtosis     () ) )

def test_moment2() :

    if root_version < 61800 :
        logger.warning ( "does not work for ancient ROOT versions")
        return

    m0  = Ostap.Math.WMoment_( 0)()
    m1  = Ostap.Math.WMoment_( 1)()
    m2  = Ostap.Math.WMoment_( 2)()
    m3  = Ostap.Math.WMoment_( 3)()
    m4  = Ostap.Math.WMoment_( 4)()
    m5  = Ostap.Math.WMoment_( 5)()
    m10 = Ostap.Math.WMoment_(10)()
    m20 = Ostap.Math.WMoment_(20)()

    m = m0, m1, m2 , m3, m4, m5 , m10, m20
    
    for i in range ( 10000 ) :
        v = random.gauss ( 0 , 1   ) 
        w = random.gauss ( 1 , 0.1 ) 
        for i in m : i.add ( v , w )  

    if root_version < 61800 :
        logger.warning ( "It does not work for ancient ROOT versions")
        return

    for i in m :
        t = 'Moment counter-%2d' % i.order 
        logger.info ( "%s\n%s" % ( t , i.table ( title = t , prefix = '# ' ) ) )

    for i in m :
        t = 'Moment counter-%2d' % i.order 
        logger                  .info (' %s : size     %s ' % ( t , i.size     () ) )
        logger                  .info (' %s : n_eff    %s ' % ( t , i.nEff     () ) )
        logger                  .info (' %s : sum(w)   %s ' % ( t , i.w        () ) )
        logger                  .info (' %s : sum(w^2) %s ' % ( t , i.w2       () ) )        
        if 1 <= i.order : logger.info (' %s : mean     %s ' % ( t , i.mean     () ) )
        if 2 <= i.order : logger.info (' %s : rms      %s ' % ( t , i.rms      () ) )
        if 2 <= i.order : logger.info (' %s : variance %s ' % ( t , i.variance () ) )
        if 3 <= i.order : logger.info (' %s : skewness %s ' % ( t , i.skewness () ) )
        if 4 <= i.order : logger.info (' %s : kurtosis %s ' % ( t , i.kurtosis () ) )


# =============================================================================
if '__main__' == __name__ :

    test_moment1()
    test_moment2()
        
# =============================================================================
##                                                                      The END 
# ============================================================================= 

