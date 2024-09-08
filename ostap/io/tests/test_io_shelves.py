#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/io/tests/test_io_shelves.py
# Test module for data storages, i.e. modules
# @see rootshelve.py
# @see zipshelve.py
# @see bz2shelve.py
# @see lzshelve.py
# @see zstshelve.py
# Copyright (c) Ostap developpers.
# =============================================================================
""" Test module for data storages, i.e. modules
  - ostap/io/rootshelve.py
  - ostap/io/zipshelve.py
  - ostap/io/bzshelve.py  
  - ostap/io/lzshelve.py   (python3 only)
  - ostap/io/zstshelve.py  (python3 only)
"""
# =============================================================================
## import sys
## sys.modules['dbhash'] = None
## sys.modules['bsddb' ] = None
## sys.modules['gdbm'  ] = None
## ##  sys.modules['dbm'   ] = None
# =============================================================================
from   ostap.math.base       import iszero
from   ostap.core.pyrouts    import VE, hID 
from   ostap.utils.timing    import timing
from   ostap.utils.utils     import random_name 
from   sys                   import version_info as python_version
import ostap.utils.cleanup   as     CU
import ostap.io.zipshelve    as     zipshelve
import ostap.io.bz2shelve    as     bz2shelve
import ostap.io.rootshelve   as     rootshelve
import ostap.logger.table    as     T 
import ROOT, os, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_io_shelves' )
else                       : logger = getLogger ( __name__          )
# =============================================================================
if  ( 3 , 3 ) <= python_version :
    import ostap.io.lzshelve  as lzshelve
    if not lzshelve.lzma : lzshelve = None 
else :
    lzshelve = None
# =============================================================================
if  ( 3 , 6 ) <= python_version :    
    import ostap.io.zstshelve as zstshelve 
    if not zstshelve.zst : zstshelve = None 
    try : 
        import zstandard
    except ImportError :
        zstshelve = None
else :
    zstshelve = None 
# =============================================================================

bins    = 1000
data    = {}
h1      = ROOT.TH1D( 'h1' ,'1D-histogram',bins,-5,5) ; h1.Sumw2() 
m1      = VE(1,2)
for i in range ( 0, 100000) : h1.Fill( m1.gauss() )

bins    = 50
h2      = ROOT.TH2D('h2','2D-histogram',bins,-5,5,bins,-5,5) ; h2.Sumw2() 
for i in range ( 0, 100000) : h2.Fill( m1.gauss() , m1.gauss() )

data [ 'histo-1D' ] = h1
data [ 'histo-2D' ] = h2
data [ 'both'     ] = (123 , h1 , {'a':2}, h2,'comment',())

data [ 'histos'   ] = {}
for i in range ( 1000 ) :
    hn = hID() 
    hh = ROOT.TH1D ( hn ,'1D-histogram' , 100 , -5, 5 )
    for j in range ( 5000 ) :  hh.Fill ( random.gauss ( 0 , 1 ) ) 
    if i < 100 : data [ hn         ] = hh
    data [ 'histos'   ] [ hn ] = hh
    
# =============================================================================
def test_shelves1():

    logger = getLogger ('Test shelves')
    logger.info ( 'Test varioouts shelves' ) 

    names  = {} 
    dbases = {}

    names      [ 'db_sql'  ] = CU.CleanUp.tempfile ( suffix = '.sqlsh'  )
    names      [ 'db_zip'  ] = CU.CleanUp.tempfile ( suffix = '.zipsh'  )
    names      [ 'db_bz2'  ] = CU.CleanUp.tempfile ( suffix = '.bz2sh'  )
    names      [ 'db_root' ] = CU.CleanUp.tempfile ( suffix = '.root'   )  
    
    dbases     [ 'db_zip'  ] = zipshelve .open ( names [ 'db_zip'  ] , 'c' )
    dbases     [ 'db_sql'  ] = zipshelve .open ( names [ 'db_sql'  ] , 'c' , dbtype = 'sqlite3' )
    dbases     [ 'db_bz2'  ] = bz2shelve .open ( names [ 'db_bz2'  ] , 'c' )
    dbases     [ 'db_root' ] = rootshelve.open ( names [ 'db_root' ] , 'c' )
    
    if lzshelve :
        names  [ 'db_lzma' ] = CU.CleanUp.tempfile ( suffix = '.lzsh'  ) 
        dbases [ 'db_lzma' ] = lzshelve  .open ( names ['db_lzma'] , 'c' )
        
    if zstshelve :
        names  [ 'db_zstd' ] = CU.CleanUp.tempfile ( suffix = '.zstsh'  )
        dbases [ 'db_zstd' ] = zstshelve .open ( names ['db_zstd'] , 'c' )

    # ===================================================================================
    ## test writing 
    # ===================================================================================

    rows = [  ( 'DBASE' , 'dbtype' , 'CPU [ms]' ) ] 
    for dbname in dbases :
        db = dbases [ dbname ]
        with timing ( 'Write %10s' % dbname ) as tm : 
            for key in data :
                db [ key ] = data [ key ]    
        logger.info ( 'DB %-25s #keys: %d' % ( dbname , len ( db  ) ) )
        row = dbname , db.dbtype , '%.1f' % ( tm.delta * 1000 )
        rows.append ( row ) 
        db.close()
        
    title = 'Write DB'
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'llr' )
    logger.info ( '%s:\n%s' % ( title , table ) ) 

    # ===================================================================================
    ## test reading
    # ===================================================================================
    
    dbases = {}

    ## open for reading 
    logger.info ( 'Reopen databases as read-only' ) 
    
    dbases     [ 'db_zip'  ] = zipshelve .open ( names [ 'db_zip'  ] , 'r' )
    dbases     [ 'db_sql'  ] = zipshelve .open ( names [ 'db_sql'  ] , 'r' , dbtype = 'sqlite3' )
    dbases     [ 'db_bz2'  ] = bz2shelve .open ( names [ 'db_bz2'  ] , 'r' )
    dbases     [ 'db_root' ] = rootshelve.open ( names [ 'db_root' ] , 'r' )
    
    if lzshelve :
        dbases [ 'db_lzma' ] = lzshelve  .open ( names ['db_lzma'] , 'r' )
    if zstshelve :
        dbases [ 'db_zstd' ] = zstshelve .open ( names ['db_zstd'] , 'r' )
                  
    rows = [  ( 'DBASE' , 'dbtype' , 'CPU [ms]' ) ] 
    for dbname in dbases :
        db = dbases [ dbname ]
        with timing ( 'Read %10s' % dbname ) as tm : 
            for key in db : value = db [ key ] 
        logger.info ( 'DB %-25s #keys: %d' % ( dbname , len ( db  ) ) )
        row = dbname , db.dbtype , '%.1f' % ( tm.delta * 1000 )
        rows.append ( row ) 
        db.close()
        
    title = 'Read DB'
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'llr' )
    logger.info ( '%s:\n%s' % ( title , table ) ) 
    
    # ===================================================================================
    ## test cloning 
    # ===================================================================================
    
    backends  = [
        'lmdb'      , 
        'berkleydb' ,
        'berkley'   ,
        'bsddb3'    ,
        'sqlite'    ,
        'sqlite3'   ,
        'dbm.gnu'   ,
        'dbm.ndbm'  ,
        'dbm.dumb'  ,
        'dbhash'    ,
        'dbm'       ,        
        'gdbm'      ,        
        'dumbdbm'   ,
        'std'       ,
    ]

    clones = {}
    
    with zipshelve .open ( names [ 'db_zip'  ] , 'r' ) as original : 

        for dbtype in backends :
            
            nz = CU.CleanUp.tempfile ( suffix = '.zip'    )
            nt = CU.CleanUp.tempfile ( suffix = '.tar'    )

            tag1 = 'db_zip/%s' % dbtype 
            nn = CU.CleanUp.tempfile ( suffix = '.zipsh'  )
            clones [ tag1 ] = original.clone ( dbname = nn , dbtype = dbtype )
            
            tagz = 'db_zip/%s/zip' % dbtype 
            nz = CU.CleanUp.tempfile ( suffix = '.zip'  )
            clones [ tagz ] = original.clone ( dbname = nz , dbtype = dbtype ) 

            tagt = 'db_zip/%s/tar' % dbtype 
            nt = CU.CleanUp.tempfile ( suffix = '.tar'  )
            clones [ tagt ] = original.clone ( dbname = nt , dbtype = dbtype ) 

    rows = [  ( 'DBASE' , 'dbtype' , 'CPU [ms]' ) ] 
    for dbname in clones  :
        db = clones [ dbname ]
        with timing ( 'Read %10s' % dbname ) as tm : 
            for key in db : value = db [ key ] 
        logger.info ( 'DB %-25s #keys: %d' % ( dbname , len ( db  ) ) )
        row = dbname , db.dbtype , '%.1f' % ( tm.delta * 1000 )
        rows.append ( row ) 
        db.close()
        
    title = 'Read clones '
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'llr' )
    logger.info ( '%s:\n%s' % ( title , table ) ) 
    
# =============================================================================

def test_shelves2 () : 

    shelves = [
        zipshelve ,
        bz2shelve , 
    ]
    if lzshelve  : shelves.append ( lzshelve  )
    if zstshelve : shelves.append ( zstshelve )

    backends  = [
        'lmdb'      , 
        'berkleydb' ,
        'berkley'   ,
        'bdsdb2'    ,
        'sqlite'    ,
        'sqlite3'   ,
        ''
    ]
    
    for sh in shelves :        
        for b in backends :            
            with sh.tmpdb ( dbtype = b ) as db :
                
                db ['one'] = 1
                db ['two'] = 2                
                db.ls() 
    
# =============================================================================
if '__main__' == __name__ :
    
    test_shelves1 ()
    test_shelves2 ()

# =============================================================================
##                                                                      The END
# =============================================================================
