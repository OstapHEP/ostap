#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/io/tests/test_io_check.py
# Check the availabel backends 
# =============================================================================
""" Check the available backends 
"""
# =============================================================================
from   ostap.core.meta_info   import python_info 
from   ostap.logger.colorized import attention 
import ostap.io.dbase         as     DB
import ostap.io.sqlitedict 
import ostap.logger.table     as     T

# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_io_check' )
else                       : logger = getLogger ( __name__         )
# =============================================================================
dbases = []

# =============================================================================
def test_io_check() :
    """ Check available  dbase backends
    """
    
    rows = [   ( 'Backend' , 'Module' , 'Location' ) ] 

    if DB.use_berkeleydb and DB.berkeleydb : 
        row = 'BerkeleyDB' , DB.berkeleydb.__name__ , DB.berkeleydb.__file__ ,
    else :
        row = attention ( 'BerkeleyDB' ) , attention ( '---' )   , attention ( '---' ) 
    rows.append ( row )

    if  python_info < ( 3 , 10 ) :
        if DB.use_bsddb3 and DB.bsddb3 : 
            row = 'BSDDB3'  , DB.bsddb3.__name__ , DB.bdsdb3.__file__ ,
        else :
            row = attention ( 'BSDDB3' ) , attention ( '---' )   , attention ( '---' ) 
        rows.append ( row )

    if DB.use_lmdb and DB.lmdb : 
        row = 'LMDB' , DB.lmdb.__name__ , DB.lmdb.__file__ ,
    else :
        row = attention ( 'LMDB' ) , attention ( '---' )   , attention ( '---' ) 
    rows.append ( row )

    row = 'SqliteDict' , ostap.io.sqlitedict.__name__ , ostap.io.sqlitedict.__file__ 
    rows.append ( row )
    
    if DB.db_gnu : 
        row = 'GNU DB' , DB.db_gnu.__name__ , DB.db_gnu.__file__ ,
    else :
        row = attention ( 'GNU DB' ) , attention ( '---' )   , attention ( '---' ) 

    rows.append ( row )

    if DB.db_dbm : 
        row = 'NDBM' , DB.db_dbm.__name__ , DB.db_dbm.__file__ ,
    else :
        row = attention ( 'NDBM' ) , attention ( '---' )   , attention ( '---' ) 

    rows.append ( row )
    
    if DB.db_dumb : 
        row = 'Dumb DB' , DB.db_dumb.__name__ , DB.db_dumb.__file__ ,
    else :
        row = attention ( 'Dumb DB' ) , attention ( '---' )   , attention ( '---' ) 

    rows.append ( row )

    title = 'Status of DB backends'
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lll' ) 
    logger.info ( '%s:\n%s' % ( title , table ) )
    
        
# =============================================================================
if '__main__' == __name__ :

    test_io_check()


# =============================================================================
##                                                                      The END 
# =============================================================================
