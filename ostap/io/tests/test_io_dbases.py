#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/io/tests/test_io_dbases.py
# Test for different databases  
# =============================================================================
""" Test module for different databases 
"""
# =============================================================================
## import sys
## sys.modules['dbhash'] = None
## sys.modules['bsddb' ] = None
## sys.modules['gdbm'  ] = None
## ##  sys.modules['dbm'   ] = None
# =============================================================================
from   ostap.utils.cleanup      import CleanUp
from   ostap.utils.progress_bar import progress_bar
from   ostap.utils.timing       import timing 
import ROOT, sys, pickle, random  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_io_dbases' )
else                       : logger = getLogger ( __name__         )
# =============================================================================
characters = "abcdefghijklmnopqrstuvwxyz0123456789_"
# =============================================================================
## generate random name 
def the_name ():
    """ Generate random name"""
    return ''.join ( random.choices ( characters , k = 8 ) )
# =============================================================================

dbases = []

# ============================================================================
try :
    from ostap.io.lmdbdict import LmdbDict
    item = 'LmdbDict' , CleanUp.tempdir ( prefix = 'ostap-LMDB-' ) , LmdbDict
    dbases.append ( item )
except ImportError:
    logger.warning ( 'LmdbDict is not accessible!' )

# ============================================================================
try :
    from ostap.io.dbase import berkeleydb_open 
    item = 'BerkeleyDB' , CleanUp.tempfile ( prefix = 'ostap-BerkeleyDB-' , suffix = '.db' ) , berkeleydb_open 
    dbases.append ( item )
except ImportError:
    logger.warning ( 'BerkeleyDB is not accessible!' )

# ============================================================================
try :
    from ostap.io.dbase import bdsdb3_open 
    item = 'BSDDB3' , CleanUp.tempfile ( prefix = 'ostap-BSDDB3-' , suffix = '.db'  ) , bsbdb3_open 
    dbases.append ( item )
except ImportError:
    logger.warning ( 'bsddb3 is not accessible!' )
   
# ============================================================================
try :
    from ostap.io.sqlitedict import SqliteDict
    item = 'SqliteDict' , CleanUp.tempfile ( prefix = 'ostap-SqliteDB-' , suffix = '.sql') , SqliteDict
    dbases.append ( item )
except ImportError:
    logger.warning ( 'SQliteDict is not accessible!' )

# ============================================================================
if ( 3,0 ) <= sys.version_info :
    # ========================================================================
    try :
        from dbm.gnu import open as _gnu_open
        def gnu_open ( filename , flag = 'r' ) : return _gnu_open ( filename , flag )        
        item = 'dbm.gnu' , CleanUp.tempfile ( prefix = 'ostap-GNUDB-' , suffix = '.db' ) , gnu_open
        dbases.append ( item )
    except ImportError:
        logger.warning ( 'dbm.gnuis not accessible!' )
    # ======================================================================
    try :
        from dbm.ndbm import open as _ndbm_open
        def ndbm_open ( filename , flag = 'r' ) : return _ndbm_open ( filename , flag )        
        item = 'dbm.ndbm' , CleanUp.tempfile ( prefix = 'ostap-NDB-' , suffix = '.db' ) , ndbm_open
        dbases.append ( item )
    except ImportError:
        logger.warning ( 'dbm.ndbm is not accessible!' )
    # ======================================================================
    try :
        from dbm.dumb import open as dumb_open
        ## def ndbm_open ( filename , flag = 'r' ) : return _ndbm_open ( filename , flag )        
        item = 'dbm.dumb' , CleanUp.tempfile ( prefix = 'ostap-DumbDB-' , suffix = '.db' ) , dumb_open
        dbases.append ( item )
    except ImportError:
        logger.warning ( 'dbm.dumb is not accessible!' )

else : ## python2
    
    # ======================================================================
    try :
        from anydbm.dbm import open as _dbm_open
        def dbm_open ( filename , flag = 'r' ) : return _dbm_open ( filename , flag )        
        item = 'anydbm.dbm' , CleanUp.tempfile ( prefix = 'ostap-DBM-' , suffix = '.db' ) , dbm_open
        dbases.append ( item )
    except ImportError:
        logger.warning ( 'anydbm.dbm accessible!' )
    # ======================================================================
    try :
        from anydbm.gdbm import open as _gdbm_open
        def gdbm_open ( filename , flag = 'r' ) : return _gdbm_open ( filename , flag )        
        item = 'anydbm.gdbm' , CleanUp.tempfile ( prefix = 'ostap-GDBM-' , suffix = '.db' ) , gdbm_open
        dbases.append ( item )
    except ImportError:
        logger.warning ( 'anydbm.gdbm accessible!' )
    # ======================================================================
    try :
        from anydbm.dbhash import open as _hash_open
        def hash_open ( filename , flag = 'r' ) : return _hash_open ( filename , flag )        
        item = 'anydbm.dbhash' , CleanUp.tempfile ( prefix = 'ostap-DBHASH-' , suffix = '.db' ) , hash_open
        dbases.append ( item )
    except ImportError:
        logger.warning ( 'anydbm.dbhash accessible!' )
    # ======================================================================


data = {
    'string'  : 'string'         ,
    'int'     : 1                , 
    'float'   : 1.0              ,
    'list'    : [1,2,3]          ,
    'set'     : set ( [1,2,3] )  ,
    'tuple'   : ('a', 'b', 'c' ) ,
    'histo'   : ROOT.TH1D ( 'h1' , '' , 100 , 0, 1 )
}


if ( 3,0 ) <= sys.version_info : 
    def the_key ( key ) :  return key.encode ( 'utf-8' )
else :
    def the_key ( key ) :  return key


dbases.reverse() 

hh = ROOT.TH1D ( 'histp' , 'title' , 100 , 0 , 10 )
for i in range ( 1000 ) : hh.Fill ( random.gauss ( 5 , 1 ) ) 
# =============================================================================
## test databases 
def test_dbases () :
    """ Test databases """

    times = [ ( 'DB' , 'write/read' , 'time [ms]' ) ]
    
    for item in dbases :
        
        tag, filename, dbase = item

        
        with timing () as tm : 
            db = dbase ( filename , flag = 'n'  )         
            for datum in data :
                
                key   = the_key ( datum ) 
                value = pickle.dumps ( datum ) 
                db [ key ] = value 
                
            ##for i in range ( 1000 ) :
            ##    kk  = the_name ()
            ##    key = the_key ( kk ) 
            ##    db [ key ] = pickle.dumps ( hh ) 
                
            if hasattr ( db , 'sync' ) : db.sync  ()
            db.close ()

        row = tag , 'write' , '%.1f' % ( tm.delta * 1000 ) 
        times.append ( row )
        
    rows = [  ( 'DB' , 'key' , 'type' , 'size' ) ]
    for item in dbases :
        
        tag, filename, dbase = item

        with timing () as tm :
            
            db = dbase ( filename , flag = 'r'  ) 
            for k in db.keys()  :
                
                key   = k
                vv    = db [ key ]
                value = pickle.loads ( vv ) 
                
                row = tag , str ( key ) , type ( value ) .__name__ , '%d' % len ( vv )
                rows.append ( row )
                
            db.close()

        row = tag , 'read' , '%.1f' % ( tm.delta * 1000 ) 
        times.append ( row )
        
    import ostap.logger.table as T
    
    title = 'Test databases'
    table = T.table ( rows ,title = title , prefix = '# ' , alignment = 'llcc' ) 
    logger.info ( '%s:\n%s' % ( title , table ) ) 

    title = 'Test CPUY'
    table = T.table ( times ,title = title , prefix = '# ' , alignment = 'lc' ) 
    logger.info ( '%s:\n%s' % ( title , table ) ) 
    

# =============================================================================
if '__main__' == __name__ :

    test_dbases ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
