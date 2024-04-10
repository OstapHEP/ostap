#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/logger/tests/test_logger_table.py
# Copyright (c) Ostap developpers.
# =============================================================================
""" Test module for tables 
"""
# =============================================================================
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_logger_table'  )
else                       : logger = getLogger ( __name__             )
# =============================================================================
import ostap.logger.table as T




# =============================================================================
## Test tables
def test_tables () :
    "Test tables"

    
    table1 = [
        ( 'Name'  , 'Occupation'        , 'Note' ) ,
        ( 'Alice' , '?'                 , '---'  ) ,
        ( 'Bob'   , 'unemployed person' , ''     )
        ]


    comment = 20 * 'This is a very long comment. '
    
    table2 = [
        ( 'Name'  , 'Occupation'        , 'Note'   , 'Comment' ) ,
        ( 'Alice' , 'unknown'           , comment  , ' '       ) ,
        ( 'Bob'   , 'unemployed person' , ''       , comment   )
        ]

    styles = ( 'local'  , 'ascii'    ,  'single' , 'porcelain' ,
               'github' , 'markdown' ,  'double' , 'default'   )
    
        
    tabfuns  = T.table , T.the_table
    tabtypes = 'table' , 'the_table'
    
    for tabfun , tabtyp in zip ( tabfuns , tabtypes ) :

        for style in styles : 

            title = "Test '%s' style: '%s'" % ( tabtyp , style ) 
        
            t = tabfun ( table1 , title = title , prefix = '#  ' , style = style )
            logger.info ( 'Table %s\n%s' % ( title , t ) )


    for tabfun , tabtyp in zip ( tabfuns , tabtypes ) :
                
        title = "Test '%s' for wrapped columns" % tabtyp
        t = tabfun ( table2 , title = title , prefix = '#  ' , alignment = 'llww')
        logger.info ( 'Table %s\n%s' % ( title , t ) )
        
# =============================================================================
if __name__ == '__main__' :

    test_tables () 


# =============================================================================
##                                                                      The END 
# =============================================================================


