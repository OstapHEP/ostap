#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_base.py
#  Test module for the file ostap/math/base.py
# ============================================================================= 
""" Test module for ostap/math/base.py
"""
# ============================================================================= 
from   ostap.logger.pretty    import pretty_float
from   ostap.logger.colorized import allright, attention
from   ostap.utils.utils      import batch_env 
import ostap.math.base        as     MB
import ostap.logger.table     as     T
import random
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_math_base' ) 
else                       : logger = getLogger ( __name__         )
# ============================================================================= 
## set batch from environment 
batch_env ( logger )
# =============================================================================
#

# =============================================================================
## test frexp10 
def test_frexp10 () :

    logger = getLogger ( 'test_frexp10' )
    
    rows = [ ( 'Value' , 'exp' , 'm/py' , 'exp/py' , 'm/C++' , 'exp/C++' , 'OK?' ) ]

    ok = True 
    for i in range ( 100 ) :
        
        mexp   = random.randrange ( -100 , 101 )
        val    = random.uniform   ( -1   , 1   )
        value  = val * 10**mexp
        m1, e1 = MB.frexp10     ( value )
        m2, e2 = MB.cpp_frexp10 ( value )
        
        v, n = pretty_float ( value )
         
        if MB.isequal ( m1 , m2 ) and MB.isequal ( e1 , e2 ) :
            good = '' 
        else :
            ok = False 
            good = attention ( 'No!' )
        
        row = v , '%+3d' % n if n else  ''              , \
            '%+.4f' % ( 10 * m1 ) , '%+3d' % ( e1 - 1 ) if e1 != 1 else '' , \
            '%+.4f' % ( 10 * m2 ) , '%+3d' % ( e2 - 1 ) if e2 != 1 else '' ,  good
        
        rows.append ( row )
        
    title = 'Test for frexp10'
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lrlrlrc' )
    if ok : logger.info  ( '%s:\n%s' % ( title , table ) )
    else  : logger.error ( '%s:\n%s' % ( title , table ) )
    
# =============================================================================
if '__main__' == __name__ :

    test_frexp10() 
    
# =============================================================================
##                                                                      The END 
# ============================================================================= 
