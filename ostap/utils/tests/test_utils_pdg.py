#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/utils/tests/test_utils_pdg.py
#  Test module for the file ostap/math/ve.py
# ============================================================================= 
""" Test module for ostap/math/ve.py
"""
# =============================================================================
from   ostap.logger.symbols   import times , ellipsis
import ostap.utils.pdg_format as     PDG
import ostap.logger.table     as     T 
# =============================================================================     
# logging 
# ============================================================================= 
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_utils_pdg' )
else                       : logger = getLogger ( __name__  )
# =============================================================================

def test_pdg ():

    ss = tuple ( [ 1.0/9  * s for s in range ( 9 ) ] + [ 0.96666666666666666666666666 ] ) 
    
    v0 = 1.0 / 3

    rows = [ ( 'Value' , 'Error', 'PDG' , 'PDG' , '%s[%s]' % ( times , ellipsis ) , 'LaTeX' ) ]
    
    for i in range ( -6 , 8 ) :

        v  = v0 * ( 10 ** i )
        
        pv = '%.6g' % v

        for j in range ( -4 , 5 ) :
            
            for s in ss :

    
                sss = s * ( 10 ** j )
            
                ps  = '%.6g' % sss

                r1 , expo = PDG.pdg_format_ ( v , sss ) 
                r         = PDG.pdg_format  ( v , sss )  
                l         = PDG.pdg_format  ( v , sss , latex = True )  
            
                row = pv , ps , r , r1 , '10^%+d' % expo if expo else '' , l 

                rows.append ( row ) 

    rows  = T.remove_empty_columns ( rows )
    title = 'Test for PDG format'
    table = T.table ( rows , title = title , prefix = '# ' ) 
    logger.info ( '%s:\n%s' % ( title , table ) )
        
# =============================================================================
if '__main__' == __name__ :

    test_pdg ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
