#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_cb.py
#  @see Ostap::Math::CrystalBall
#  @see Ostap::Math::CrystalBallRightside 
#  @see Ostap::Math::CrystalBallDoubleSided 
# ============================================================================= 
""" Test module for CrystallBall functions 
- see Ostap::Math::CrystalBall
- see Ostap::Math::CrystalBallRightside 
- see Ostap::Math::CrystalBallDoubleSided 
"""
# ============================================================================= 
from   ostap.core.core          import Ostap, SE
from   ostap.utils.basic        import typename 
from   ostap.utils.root_utils   import batch_env
import ostap.math.models
import random
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_cb'   )
else                       : logger = getLogger ( __name__                   )
# ============================================================================= 
## set batch from environment 
batch_env ( logger )
# =============================================================================
import ostap.math.integral      as I 

m0    = 1
sigma = 0.5

xmin = m0 - 6.1 * sigma 
xmax = m0 + 6.1 * sigma

def the_tst ( fun , logger = logger ) :

    cnt1 = SE()
    cnt2 = SE()
    
    while cnt1.nEntries () < 100000 : 
        
        x1 = random.uniform ( xmin , xmax )
        x2 = random.uniform ( x1   , xmax )
        
        ## skip very long intervals  
        if x2 - x1 > 1 * sigma  : continue
            
        i1 = fun.integral   ( x1   , x2   )
        i2 = I.integral     ( fun  , x1 , x2 )

        d  = i1 - i2
        r  = i1 / i2
        
        cnt1 += d
        cnt2 += ( r - 1 ) 
        
    delta = 1.e-5

    mnmx1 = cnt1.max() - cnt1.min()
    mnmx2 = cnt2.max() - cnt2.min()

    tn = typename ( fun ) 
    if   delta > mnmx1 : logger.info    ( 'DIFFERENCE %s\n%s'       % ( tn ,         cnt1 ) ) 
    else               : logger.warning ( 'DIFFERENCE %s %.4g \n%s' % ( tn , mnmx1 , cnt1 ) )
    
    if   delta > mnmx2 : logger.info    ( 'RATIO-1    %s\n%s'       % ( tn , cnt2         ) ) 
    else               : logger.warning ( 'RATIO-1    %s %.4g \n%s' % ( tn , mnmx2 , cnt2 ) )


# ===============================================================================
def test_cb () :

    logger = getLogger ( 'test_cb' )

    cb = Ostap.Math.CrystalBall ( m0 , sigma , 1.5 , 0.1  ) 
    the_tst ( cb , logger )


# ===============================================================================
def test_cbrs () :
    
    logger = getLogger ( 'test_cbrs' )
    
    cb = Ostap.Math.CrystalBallRightSide  ( m0 , sigma , 1.5 , 0.5 ) 
    
    return the_tst ( cb , logger )

# ===============================================================================
def test_cbds () :

    logger = getLogger ( 'test_cbds' )

    cb = Ostap.Math.CrystalBallDoubleSided  ( m0 , sigma , 1.5 , 0.1 , 1.5 , 4.2 ) 

    return the_tst ( cb , logger )

# =============================================================================
if '__main__' == __name__ :
    
    test_cb   ()
    test_cbrs ()
    test_cbds ()

# =============================================================================
##                                                                      The END 
# ============================================================================= 
