#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/io/tests/test_io_shelves.py
# Test module for data storages, i.e. modules
# @see rootshelve.py
# @see sqliteshelve.py
# @see zipshelve.py
# @see bz2shelve.py
# @see lzshelve.py
# Copyright (c) Ostap developpers.
# =============================================================================
""" Test module for data storages, i.e. modules
  - ostap/io/rootshelve.py
  - ostap/io/sqliteshelve.py
  - ostap/io/zipshelve.py
  - ostap/io/bzshelve.py  
  - ostap/io/lzshelve.py  (python3 only)
"""
# =============================================================================
import os, random
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
import ostap.utils.cleanup   as     CU 
import ostap.io.zipshelve    as     zipshelve
import ostap.io.bz2shelve    as     bz2shelve
if  2 < python_version.major :
    import ostap.io.lzshelve as     lzshelve
else :
    lzshelve = None 
import ostap.io.sqliteshelve as     sqliteshelve
import ostap.io.rootshelve   as     rootshelve

# =============================================================================
 
db_sql_name  = CU.CleanUp.tempfile ( suffix = '.sqldb'  )  
db_zip_name  = CU.CleanUp.tempfile ( suffix = '.zipdb'  ) 
db_bz2_name  = CU.CleanUp.tempfile ( suffix = '.bz2db'  )
db_root_name = CU.CleanUp.tempfile ( suffix = '.root'   )
db_lz_name   = CU.CleanUp.tempfile ( suffix = '.lzmadb' )

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
data [ 'both'     ] = (123 , h1 , {'a':2}, h2,'comment',())
data [ 'histos'   ] = {}
for i in range ( 5000 ) : 
    ht = 'histo#%d' % i
    hh = ROOT.TH1D ( ht , '' , 500 , 0 , 100 )
    for j in range ( 200 ) :
        hh.Fill ( random.gauss ( 50 , 10) )
    data['histos'][ht] = hh 

def test_shelves():
    
    db_sql  = sqliteshelve.open ( db_sql_name  , 'c' )
    db_zip  = zipshelve.open    ( db_zip_name  , 'c' )
    db_bz2  = bz2shelve.open    ( db_bz2_name  , 'c' )
    db_root = rootshelve.open   ( db_root_name , 'c' )
    
    if lzshelve : db_lz = lzshelve.open ( db_lz_name , 'c' )
    else        : db_ls = None 
        
        
    for k in data :
        db_sql  [ k ] = data[k]
        db_zip  [ k ] = data[k]
        db_bz2  [ k ] = data[k]
        if lzshelve :
            db_lz  [ k ] = data[k]
        db_root [ k ] = data[k]
        
    logger.info('SQLiteShelve #keys: %s' % len ( list ( db_sql .keys() ) ) ) 
    logger.info('ZipShelve    #keys: %s' % len ( db_zip .keys() ) )
    logger.info('Bz2Shelve    #keys: %s' % len ( db_bz2 .keys() ) )
    logger.info('RootShelve   #keys: %s' % len ( db_root.keys() ) )
    if lzshelve :
        logger.info('LzShelve     #keys: %s' % len ( db_lz .keys() ) )

    db_sql  .close() 
    db_zip  .close()
    db_bz2  .close()
    db_root .close()
    if lzshelve : db_lz .close()
    
    logger.info('SQLiteShelve size: %d ' % os.path.getsize ( db_sql_name  ) )
    logger.info('ZipShelve    size: %d ' % os.path.getsize ( db_zip_name  ) )    
    logger.info('Bz2Shelve    size: %d ' % os.path.getsize ( db_bz2_name  ) )    
    logger.info('RootShelve   size: %d ' % os.path.getsize ( db_root_name ) )
    if lzshelve :
        logger.info('LzShelve     size: %d ' % os.path.getsize ( db_lz_name    ) )    
    
    db_sql  = sqliteshelve.open    ( db_sql_name  , 'r' )
    db_zip  = zipshelve.open       ( db_zip_name  , 'r' )
    db_bz2  = bz2shelve.open       ( db_bz2_name  , 'r' )
    if lzshelve :
        db_lz  = lzshelve.open     ( db_lz_name   , 'r' )
    db_root = rootshelve.open      ( db_root_name , 'r' )

    logger.info('SQLiteShelve #keys: %s' % len ( list ( db_sql .keys() ) ) ) 
    logger.info('ZipShelve    #keys: %s' % len ( db_zip .keys() ) )
    logger.info('Bz2Shelve    #keys: %s' % len ( db_bz2 .keys() ) )
    if lzshelve :
        logger.info('LzShelve     #keys: %s' % len ( db_lz  .keys() ) )
    logger.info('RootShelve   #keys: %s' % len ( db_root.keys() ) )

    
    with timing ( 'h2-read/SQL'  ) : h2_sql  = db_sql  [ 'histo-2D']
    with timing ( 'h2_read/ZIP'  ) : h2_zip  = db_zip  [ 'histo-2D']
    with timing ( 'h2_read/BZ2'  ) : h2_bz2  = db_bz2  [ 'histo-2D']
    if lzshelve :
        with timing ( 'h2_read/LZ'  ) :
            h2_lz = db_lz  [ 'histo-2D']
    with timing ( 'h2_read/ROOT' ) : h2_root = db_root [ 'histo-2D']

    with timing ( 'tu-read/SQL'  ) : tu_sql  = db_sql  [ 'both'    ]
    with timing ( 'tu_read/ZIP'  ) : tu_zip  = db_zip  [ 'both'    ] 
    with timing ( 'tu_read/BZ2'  ) : tu_bz2  = db_bz2  [ 'both'    ] 
    if lzshelve :
        with timing ( 'tu_read/LZ'   ) :
            tu_lz   = db_lz   [ 'both'    ] 
    with timing ( 'tu_read/ROOT' ) : tu_root = db_root [ 'both'    ]

    with timing ( 'h1-read/SQL'  ) : h1_sql  = db_sql  [ 'histo-1D']
    with timing ( 'h1-read/ZIP'  ) : h1_zip  = db_zip  [ 'histo-1D']
    with timing ( 'h1-read/BZ2'  ) : h1_bz2  = db_bz2  [ 'histo-1D']
    if lzshelve : 
        with timing ( 'h1-read/LZ'   ) :
            h1_lz   = db_lz   [ 'histo-1D']
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
        v = h1_bz2  [i] - h1     [i] 
        if not iszero ( v.value() ) :
            logger.error('Large difference for 1D histogram(4)!')
        if lzshelve :
            v = h1_lz  [i] - h1     [i] 
            if not iszero ( v.value() ) :
                logger.error('Large difference for 1D histogram(5)!')
                
    
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
        v = h2_bz2  [i] - h2    [i] 
        if not iszero ( v.value() ) :
            logger.error('Large difference for 2D histogram(4)!')
        if lzshelve :
            v = h2_lz  [i] - h2    [i] 
            if not iszero ( v.value() ) :
                logger.error('Large difference for 2D histogram(5)!')
            

    h1tq = tu_sql [1]
    h1tz = tu_zip [1]
    h1tr = tu_root[1]
    
    with timing('Close SQL'  ) : db_sql .close() 
    with timing('Close ZIP'  ) : db_zip .close()
    with timing('Close BZ2'  ) : db_bz2 .close()
    if lzshelve : 
        with timing('Close LZ'   ) : db_lz  .close()
    with timing('Close ROOT' ) : db_root.close()


    for dbase in ( sqliteshelve . tmpdb () ,
                   zipshelve    . tmpdb () ,
                   bz2shelve    . tmpdb () ,
                   ## lzshelve     . tmpdb () ,
                   rootshelve   . tmpdb () ) :

        with timing () :
            
            with dbase as db :
                
                db [ 'h1'    ] = h1
                db [ 'h2'    ] = h2
                db [ 'data'  ] = data
                db [ 'histos'] = data['histos']
                db.ls()
    
# =============================================================================
if '__main__' == __name__ :    
    test_shelves()

    del h1 , h2

# =============================================================================
# The END
# =============================================================================
