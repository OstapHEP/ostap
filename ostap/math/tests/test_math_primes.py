#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_primes.py
#  Test module for the file ostap/math/primes.py
# ============================================================================= 
""" Test module for ostap/math/primes.py

It tests local implementation of prime numbers sieve 
"""
# ============================================================================= 
from __future__ import print_function
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_math_primes'   )
else                       : logger = getLogger ( __name__             )
# ============================================================================= 
from ostap.math.primes import primes, Primes 

def test_primes1 () :

    p = primes ( 1000 )
    logger.info  ('Primes:\n%s'  % p )

def test_primes2() :

    p = Primes ( 10000 )

    sump = 0 
    for n in p :  sump+= n
    
    logger.info ('Sum of primes          : %s' % sump )
    
    sump = 0 
    for n in p.range( 100 , 1000 )  :  sump+= n
    logger.info ('Sum of primes [0:1000] : %s' % sump )


    logger.info ('primes[100]            : %s' % p[100]     )
    logger.info ('choice                 : %s' % p.choice() )
    logger.info ('choice[0:100]          : %s' % p.choice(0,100) )

    
# =============================================================================
if '__main__' == __name__ :

    test_primes1 ()
    test_primes2 ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
