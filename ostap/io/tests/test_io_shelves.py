#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/io/tests/test_io_shelves.py
# Test module for data storages, i.e. modules
# @see rootshelve.py
# @see sqliteshelve.py
# @see zipshelve.py
# Copyright (c) Ostap developpers.
# =============================================================================
"""
Test module for data storages, i.e. modules
  - ostap/io/rootshelve.py
  - ostap/io/sqliteshelve.py
  - ostap/io/zipshelve.py
"""
# =============================================================================
import os, tempfile, random
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_io_shelves' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
from   ostap.math.base       import iszero
from   ostap.core.pyrouts    import VE
from   ostap.utils.timing    import timing 
from   sys                   import version_info as python_version 

import ostap.io.zipshelve    as     zipshelve
import ostap.io.sqliteshelve as     sqliteshelve
import ostap.io.rootshelve   as     rootshelve

# =============================================================================

db_sql_name  = tempfile.mktemp ( suffix = '.msql' ) 
db_zip_name  = tempfile.mktemp ( suffix = '.zdb'  )
db_root_name = tempfile.mktemp ( suffix = '.root' )

bins    = 1000
data    = {}
h1      = ROOT.TH1D('h1','1D-histogram',bins,-5,5) ; h1.Sumw2() 
m1      = VE(1,2)
for i in range ( 0, 100000) : h1.Fill( m1.gauss() )

bins    = 50
h2      = ROOT.TH2D('h2','2D-histogram',bins,-5,5,bins,-5,5) ; h2.Sumw2() 
for i in range ( 0, 100000) : h2.Fill( m1.gauss() , m1.gauss() )

data [ 'histo-1D' ] = h1
data [ 'histo-2D' ] = h2
data [ 'both'     ] = (123 , h1 , {'a':2}, h2,'comment')

def test_shelves():
    
    db_sql  = sqliteshelve.open    ( db_sql_name  , 'c' )
    db_zip  = zipshelve.open       ( db_zip_name  , 'c' )
    db_root = rootshelve.open      ( db_root_name , 'c' )
    
        
    for k in data :
        db_sql  [ k ] = data[k]
        db_zip  [ k ] = data[k]
        db_root [ k ] = data[k]

    logger.info('SQLiteShelve keys: %s' % list ( db_sql .keys() ) ) 
    logger.info('ZipShelve    keys: %s' % list ( db_zip .keys() ) )
    logger.info('RootShelve   keys: %s' % list ( db_root.keys() ) )

    db_sql  .close() 
    db_zip  .close()
    db_root .close()
    
    logger.info('SQLiteShelve size: %d ' % os.path.getsize ( db_sql_name  ) )
    logger.info('ZipShelve    size: %d ' % os.path.getsize ( db_zip_name  ) )    
    logger.info('RootShelve   size: %d ' % os.path.getsize ( db_root_name ) )
    
    db_sql  = sqliteshelve.open    ( db_sql_name  , 'r' )
    db_zip  = zipshelve.open       ( db_zip_name  , 'r' )
    db_root = rootshelve.open      ( db_root_name , 'r' )

    logger.info('SQLiteShelve keys: %s' % list ( db_sql .keys() ) ) 
    logger.info('ZipShelve    keys: %s' % list ( db_zip .keys() ) )
    logger.info('RootShelve   keys: %s' % list ( db_root.keys() ) )

    
    with timing ( 'h2-read/SQL'  ) : h2_sql  = db_sql  [ 'histo-2D']
    with timing ( 'h2_read/ZIP'  ) : h2_zip  = db_zip  [ 'histo-2D']
    with timing ( 'h2_read/ROOT' ) : h2_root = db_root [ 'histo-2D']

    with timing ( 'tu-read/SQL'  ) : tu_sql  = db_sql  [ 'both'    ]
    with timing ( 'tu_read/ZIP'  ) : tu_zip  = db_zip  [ 'both'    ] 
    with timing ( 'tu_read/ROOT' ) : tu_root = db_root [ 'both'    ]

    with timing ( 'h1-read/SQL'  ) : h1_sql  = db_sql  [ 'histo-1D']
    with timing ( 'h1-read/ZIP'  ) : h1_zip  = db_zip  [ 'histo-1D']
    with timing ( 'h1-read/ROOT' ) : h1_root = db_root [ 'histo-1D']

    for i in h1_sql : 
        v = h1_sql  [i] - h1_zip [i] 
        if not iszero ( v.value() ) :
            logger.error('Large difference for 1D histogram(1)!')
        v = h1_sql  [i] - h1     [i] 
        if not iszero ( v.value() ) :
            logger.error('Large difference for 1D histogram(2)!')
        v = h1_root [i] - h1     [i] 
        if not iszero ( v.value() ) :
            logger.error('Large difference for 1D histogram(3)!')

    
    for i in h2_sql : 
        v = h2_sql  [i] - h2_zip[i] 
        if not iszero ( v.value() ) :
            logger.error('Large difference for 2D histogram(1)!')
        v = h2_sql  [i] - h2    [i] 
        if not iszero ( v.value() ) :
            logger.error('Large difference for 2D histogram(2)!')
        v = h2_root [i] - h2    [i] 
        if not iszero ( v.value() ) :
            logger.error('Large difference for 2D histogram(3)!')


    h1tq = tu_sql [1]
    h1tz = tu_zip [1]
    h1tr = tu_root[1]
    
    with timing('Close SQL'  ) : db_sql .close() 
    with timing('Close ZIP'  ) : db_zip .close()
    with timing('Close ROOT' ) : db_root.close()

    with timing('Remove SQL' ) : os.remove ( db_sql_name  )
    with timing('Remove ZIP' ) : os.remove ( db_zip_name  )
    with timing('Remove ROOT') : os.remove ( db_root_name )


    for dbase in ( sqliteshelve.tmpdb() ,
                   zipshelve   .tmpdb() ,
                   rootshelve  .tmpdb() ) :

        with dbase as db :
            
            db [ 'h1'] = h1
            db [ 'h2'] = h2
            db.ls()
    
# =============================================================================
if '__main__' == __name__ :    
    test_shelves()

    del h1 , h2

# =============================================================================
# The END
# =============================================================================
