#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_fill.py
# Test module for filling dataset
# - make some fitting toys 
# ============================================================================= 
""" Test module for filling datasets 
- fill some datasets 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.stats.counters   import SE, WSE
from   ostap.utils.root_utils import batch_env 
import ostap.logger.table     as     T
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_stats_counters' )
else : 
    logger = getLogger ( __name__              )
# =============================================================================
batch_env ( logger ) 
# =============================================================================

threshold = 1.e-14 
# =============================================================================
def test_stats_counters_1 () :
    
    logger = getLogger("tests_stats_counters_1")
    
    data = [ random.uniform ( -1 , 1 ) for i in range ( 10000 ) ]
    
    counters = []
    for i in range ( 1000 ) :
        random.shuffle ( data )
        c = SE.count   ( data )
        counters.append  ( c ) 
        
    ## mu values
    c_mu  = SE.count ( ( c.mu  () for c in counters ) )
    c_mu2 = SE.count ( ( c.mu2 () for c in counters ) )

    rms_mu  = c_mu .rms()
    rms_mu2 = c_mu2.rms()    
    if threshold <= rms_mu  : logger.error ( "RMS for `mu'  : %.3g" % rms_mu  )
    else                    : logger.info  ( "RMS for `mu'  : %.3g" % rms_mu  )
    if threshold <= rms_mu2 : logger.error ( "RMS for `mu2' : %.3g" % rms_mu2 )
    else                    : logger.info  ( "RMS for `mu2' : %.3g" % rms_mu2 )

# =============================================================================
def test_stats_counters_2 () :
    
    logger = getLogger("tests_stats_counters_2")
    
    data = [ random.uniform ( -1 , 1 ) for i in range ( 20000 ) ]
    L    = len ( data )
    
    counters = []
    for i in range ( 1000 ) :
        
        i1 = random.randrange (      L )
        i2 = random.randrange ( i1 , L )
        i3 = random.randrange ( i2 , L )
        i4 = random.randrange ( i3 , L )

        c  = SE.count ( data[  :i1] ) + \
             SE.count ( data[i1:i2] ) + \
             SE.count ( data[i2:i3] ) + \
             SE.count ( data[i3:i4] ) + \
             SE.count ( data[i4:  ] ) 

        counters.append  ( c ) 
        
    ## mu values
    c_mu  = SE.count ( ( c.mu  () for c in counters ) )
    c_mu2 = SE.count ( ( c.mu2 () for c in counters ) )

    rms_mu  = c_mu .rms()
    rms_mu2 = c_mu2.rms()
    
    if threshold <= rms_mu  : logger.error ( "RMS for ``mu''  : %.3g" % rms_mu  )
    else                    : logger.info  ( "RMS for ``mu''  : %.3g" % rms_mu  )
    if threshold <= rms_mu2 : logger.error ( "RMS for ``mu2'' : %.3g" % rms_mu2 )
    else                    : logger.info  ( "RMS for ``mu2'' : %.3g" % rms_mu2 )

# =============================================================================
def test_stats_counters_3 () :
    
    logger = getLogger("tests_stats_counters_3")
    
    data = [ ( random.uniform ( -1 , 1 ) , random.uniform ( 0 , 2 ) )   for i in range ( 10000 ) ]
    
    counters = []
    for i in range ( 1000 ) :
        random.shuffle ( data )
        c = WSE ()
        for v , w in data : c.add ( v , w ) 
        counters.append  ( c ) 
        
    c_mu  = SE.count ( ( c.mu  () for c in counters ) )
    c_mu2 = SE.count ( ( c.mu2 () for c in counters ) )

    rms_mu  = c_mu .rms()
    rms_mu2 = c_mu2.rms()
    
    if threshold <= rms_mu  : logger.error ( "RMS for ``mu''  : %.3g" % rms_mu  )
    else                    : logger.info  ( "RMS for ``mu''  : %.3g" % rms_mu  )
    if threshold <= rms_mu2 : logger.error ( "RMS for ``mu2'' : %.3g" % rms_mu2 )
    else                    : logger.info  ( "RMS for ``mu2'' : %.3g" % rms_mu2 )
    
# =============================================================================
def test_stats_counters_4 () :
    
    logger = getLogger("tests_stats_counters_4")

    data = [ ( random.uniform ( -1 , 1 ) , random.uniform ( 0 , 2 ) )   for i in range ( 10000 ) ]
    L    = len ( data )
    
    counters = []
    for i in range ( 1000 ) :
        
        i1 = random.randrange (      L )
        i2 = random.randrange ( i1 , L )
        i3 = random.randrange ( i2 , L )
        i4 = random.randrange ( i3 , L )

        c1 = WSE ()
        c2 = WSE ()
        c3 = WSE ()
        c4 = WSE ()
        c5 = WSE ()
        
        for v , w in data [  :i1] : c1.add ( v , w ) 
        for v , w in data [i1:i2] : c2.add ( v , w ) 
        for v , w in data [i2:i3] : c3.add ( v , w ) 
        for v , w in data [i3:i4] : c4.add ( v , w ) 
        for v , w in data [i4:  ] : c5.add ( v , w ) 

        c = c1 + c2 + c3 + c4 + c5 
        counters.append  ( c ) 

        
    c_mu  = SE.count ( ( c.mu  () for c in counters ) )
    c_mu2 = SE.count ( ( c.mu2 () for c in counters ) )

    rms_mu  = c_mu .rms()
    rms_mu2 = c_mu2.rms()
    
    if threshold <= rms_mu  : logger.error ( "RMS for ``mu''  : %.3g" % rms_mu  )
    else                    : logger.info  ( "RMS for ``mu''  : %.3g" % rms_mu  )
    if threshold <= rms_mu2 : logger.error ( "RMS for ``mu2'' : %.3g" % rms_mu2 )
    else                    : logger.info  ( "RMS for ``mu2'' : %.3g" % rms_mu2 )




# =============================================================================
if "__main__" == __name__ :

    test_stats_counters_1 ()
    test_stats_counters_2 ()
    test_stats_counters_3 ()
    test_stats_counters_4 ()
    
    
    

# =============================================================================
##                                                                      The END 
# =============================================================================
