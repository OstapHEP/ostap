#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_callable.py
#  @see Ostap::Functions::PyCallable
#  @see Ostap::Functions::PyCallable2
#  @see Ostap::Functions::PyCallable3
# ============================================================================= 
""" Test module for pythonkic/C++ callable wrappers 
- see Ostap.Functions.PyCallable
- see Ostap.Functions.PyCallable2
- see Ostap.Fumctions.PyCallable3
"""
# =============================================================================
from   ostap.core.pyrouts       import Ostap, SE  
from   ostap.utils.timing       import timing
from   ostap.utils.progress_bar import progress_bar
from   ostap.math.make_fun      import make_fun1, make_fun2 , make_fun3
from   ostap.utils.memory       import memory
from   ostap.utils.root_utils   import batch_env 
import ostap.logger.table       as     T
import ROOT, math, random  
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_math_callable' )
else                       : logger = getLogger ( __name__             )
# ============================================================================= 
 ## set batch from environment 
batch_env ( logger )
# =============================================================================
#
   
## 1D function
def my_fun1 ( x ) :
    return math.exp ( x )

## 2D function
def my_fun2 ( x , y ) :
    return x + y 

## 3D function
def my_fun3 ( x , y , z ) :
    return x + y + z  

# ==========================================================================
## Test 1D callable 
def test_callable1 () :
    """Test 1D callable"""
    logger = getLogger ( "test_callable1" )
    logger.info ( 'Test Ostap.functions.PyCallable' ) 
    the_fun = Ostap.Functions.PyCallable ( my_fun1 , True )
    tag     = "PyCallable"

    s = 0 
    with timing ( tag , logger = logger ) , memory ( tag , logger = logger ) :
        for i in range ( 1000000 ) :
            s += the_fun ( random.uniform ( -1 , 1 ) )

# ==========================================================================
## Test 2D callable 
def test_callable2 () :
    """Test 2D callable"""
    logger = getLogger ( "test_callable2" )
    logger.info ( 'Test Ostap.functions.PyCallable2' ) 

    the_fun = Ostap.Functions.PyCallable2 ( my_fun2 , True )
    tag     = "PyCallable3"

    s = 0 
    with timing ( tag , logger = logger ) , memory ( tag , logger = logger ) :
        for i in range ( 1000000 ) :
            s += the_fun ( random.uniform ( -1 , 1 ) ,
                           random.uniform ( -1 , 1 ) )
                           

# ==========================================================================
## Test 3d callable 
def test_callable3 () :
    """Test 3D callable"""
    logger = getLogger ( "test_callable3" )
    logger.info ( 'Test Ostap.functions.PyCallable3' ) 

    the_fun = Ostap.Functions.PyCallable3 ( my_fun3 , True )

    tag = "PyCallable3"

    s   = 0    
    with timing ( tag , logger = logger ) , memory ( tag , logger = logger ) :
        for i in range ( 1000000 ) :
            s += the_fun ( random.uniform ( -1 , 1 ) ,
                           random.uniform ( -1 , 1 ) ,
                           random.uniform ( -1 , 1 ) )
            
# ==========================================================================
## Test 1D callable 
def test_callable_f1 () :
    """Test 1D callable"""
    logger = getLogger ( "test_callable_f1" )
    the_fun = make_fun1 ( my_fun1  )
    logger.info ( 'Test make_fun1: %s' % type ( the_fun )  ) 
    tag     = "make_fun1"

    s = 0 
    with timing ( tag , logger = logger ) , memory ( tag , logger = logger ) :
        for i in range ( 1000000) :
            s += the_fun ( random.uniform ( -1 , 1 ) )
            
# ==========================================================================
## Test 2D callable 
def test_callable_f2 () :
    """Test 2D callable"""
    logger  = getLogger ( "test_callable_f2" )    
    the_fun = make_fun2 ( my_fun2  )
    logger.info ( 'Test make_fun2: %s' % type ( the_fun )  ) 
    tag     = "make_fun2"

    s = 0 
    with timing ( tag , logger = logger ) , memory ( tag , logger = logger ) :
        for i in range ( 1000000 ) :
            s += the_fun ( random.uniform ( -1 , 1 ) ,
                           random.uniform ( -1 , 1 ) )
                           

# ==========================================================================
## Test 3d callable 
def test_callable_f3() :
    """Test 3D callable"""
    logger  = getLogger ( "test_callable_f3" )

    the_fun = make_fun3 ( my_fun3  )
    logger.info ( 'Test make_fun3: %s' % type ( the_fun )  ) 
    tag     = "make_fun3"

    s   = 0    
    with timing ( tag , logger = logger ) , memory ( tag , logger = logger ) :
        for i in range ( 1000000 ) :
            s += the_fun ( random.uniform ( -1 , 1 ) ,
                           random.uniform ( -1 , 1 ) ,
                           random.uniform ( -1 , 1 ) )
            
# =============================================================================
if '__main__' == __name__ :

    test_callable1   () 
    test_callable2   () 
    test_callable3   () 

    test_callable_f1 () 
    test_callable_f2 () 
    test_callable_f3 () 

# =============================================================================
##                                                                      The END 
# ============================================================================= 
