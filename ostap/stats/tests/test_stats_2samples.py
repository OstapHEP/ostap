#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==============================================================================
# @file ostap/stats/tests/test_stats_gof1d.py 
# Test Goddness-of-fits for 1D fits 
# Copyright (c) Ostap developpers.
# ==============================================================================
""" Test Goddness-of-fits for 1D fits 
"""
# ==============================================================================
from   ostap.stats.twosamples import TSTest, TSToys 
from   ostap.core.meta_info   import python_info 
from   ostap.plotting.canvas  import use_canvas
from   ostap.logger.pretty    import pretty_float
from   ostap.utils.utils      import vrange 
from   ostap.math.math_ve     import significance
import ostap.logger.table     as     T
import ROOT, random, array    
# ==============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'tests_stats_2samples' )
else                       : logger = getLogger ( __name__            )
# ==============================================================================
N1  = 200
N2  = 200

ds1 = array.array ( 'd' , ( random.gauss ( 0   , 1   ) for i in range ( N1 ) ) ) ## Gaussian 
ds2 = array.array ( 'd' , ( random.gauss ( 0   , 1   ) for i in range ( N2 ) ) ) ## the same  Gaussian 
ds3 = array.array ( 'd' , ( random.gauss ( 0.5 , 1   ) for i in range ( N2 ) ) ) ## different Gaussian 
ds4 = array.array ( 'd' , ( random.gauss ( 0   , 1.3 ) for i in range ( N2 ) ) ) ## different Gaussian 

# =============================================================================
def test_2samples_good () :
    
    logger = getLogger ( 'test_2samples_good' ) 

    test   = TSTest ( ds1 , ds2 )
    title  = 'Test good samples'
    logger.info ( '%s:\n%s' % ( title , test.table ( title = title , prefix = '# ' ) ) )

    toys   = TSToys ( ds1 , ds2 )
    title  = 'Toys good samples'
    logger.info ( '%s:\n%s' % ( title , toys.table ( title = title , prefix = '# ' ) ) )

    with use_canvas ( 'test_2samples_good/test' , wait = 3 ) : test.draw() 
    with use_canvas ( 'test_2samples_good/toys' , wait = 3 ) : toys.draw() 

    
# =============================================================================
def test_2samples_same () :
    
    logger = getLogger ( 'test_2samples_same' ) 
    
    test   = TSTest ( ds1 , ds1 )                   
    title  = 'Test same samples'
    logger.info ( '%s:\n%s' % ( title , test.table ( title = title , prefix = '# ' ) ) )

    toys   = TSToys ( ds1 , ds1 )
    title  = 'Toys same samples'
    logger.info ( '%s:\n%s' % ( title , toys.table ( title = title , prefix = '# ' ) ) )

    with use_canvas ( 'test_2samples_same/test' , wait = 3 ) : test.draw() 
    with use_canvas ( 'test_2samples_same/toys' , wait = 3 ) : toys.draw() 
    
# =============================================================================
def test_2samples_bad1 () :
    
    logger = getLogger ( 'test_2samples_bad1'  ) 
    
    test   = TSTest ( ds1 , ds3 )
    title  = 'Test bad/1 samples'
    logger.info ( '%s:\n%s' % ( title , test.table ( title = title , prefix = '# ' ) ) )

    toys   = TSToys ( ds1 , ds3 )
    title  = 'Toys bad/1 samples'
    logger.info ( '%s:\n%s' % ( title , toys.table ( title = title , prefix = '# ' ) ) )
    
    with use_canvas ( 'test_2samples_bad1/test' , wait = 3 ) : test.draw() 
    with use_canvas ( 'test_2samples_bad1/toys' , wait = 3 ) : toys.draw() 
    
# =============================================================================
def test_2samples_bad2 () :
    
    logger = getLogger ( 'test_2samples_bad2'  ) 
    
    test   = TSTest ( ds1 , ds4  )
    title  = 'Test bad/2 samples'
    logger.info ( '%s:\n%s' % ( title , test.table ( title = title , prefix = '# ' ) ) )

    toys   = TSToys ( ds1 , ds4 )
    title  = 'Toys bad/2 samples'
    logger.info ( '%s:\n%s' % ( title , toys.table ( title = title , prefix = '# ' ) ) )
    
    with use_canvas ( 'test_2samples_bad2/test' , wait = 3 ) : test.draw() 
    with use_canvas ( 'test_2samples_bad2/toys' , wait = 3 ) : toys.draw() 

# ===============================================================================
if '__main__' == __name__ :

    test_2samples_good ()
    test_2samples_same ()
    test_2samples_bad1 ()
    test_2samples_bad2 ()
    
# ===============================================================================
##                                                                        The END 
# ===============================================================================
                   
