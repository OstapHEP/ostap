#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/io/dbase.py
# Helper module to use databases
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date   2020-05-16
# =============================================================================
""" Helper module to use databases
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2020-05-16"
__version__ = "$Revision:$" 
# =============================================================================
__all__ = (
    'whichdb'              , ## guess the database type
    'isdbase'              , ## ditto 
    'dbopen'               , ## open database
    'Item'                 , ## item: named tuple (creation time & payload)
    'TmpDB'                , ## mixin for temporary database
    'available_backends'   , ## List of available DBASE backends 
    'preferrable_backends' , ## List of available DBASE backends 
    )
# =============================================================================
import sys, os, collections
import ostap.io.shelve_ext
# =============================================================================
from   ostap.logger.logger  import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.dbase' )
else                      : logger = getLogger ( __name__         )
# =============================================================================
try : # =======================================================================
    import sqlite3 
    from ostap.io.sqlitedict import SQLiteDict, issqlite3
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    SQLiteDict, issqlite3 = None , None 
# =============================================================================
## named tuple to DB-item: (time, payload)
Item = collections.namedtuple ( 'Item', ( 'time' , 'payload' ) )
# =============================================================================
ordered_dict = dict
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import dbm.gnu  as db_gnu
    db_gnu_types = 'dbm.gnu'  , 'gdbm' , 'gnu' , 'dbm'
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    db_gnu = None
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import dbm.ndbm as db_dbm
    db_ndbm_types =  'dbm.ndbm' , 'ndbm' 
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    db_dbm = None
# =============================================================================
import dbm.dumb as db_dumb
db_dumb_types = 'dbm.dumb' , 'dumbdbm' , 'dumb' 
db_std_types  = 'std' , 'stddb' , 'standard' 
# =============================================================================
## Check sqlite3 
# =============================================================================
if ( 3 , 13 ) <= sys.version_info : # =========================================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import dbm.sqlite3
        def sqlite3_open ( filename     ,
                           flag = 'c'   ,
                           mode = 0o660 , **kwargs ) :
            ## 
            return dbm.sqlite3.open ( filename    ,
                                      flag = flag ,
                                      mode = mode )
        use_sqlite3 = True 
        # =====================================================================
    except ImportError : # ====================================================
        # =====================================================================       
        use_sqlite3 = False
    # =========================================================================
else : # ======================================================================
    # =========================================================================
    use_sqlite3 = False
# =============================================================================
if use_sqlite3 or SQLiteDict :
    db_sqlite3_types = 'dbm.sqlite3' , 'sqlite3' , 'sqlite'  , 'sql3' , 'sql'
# =============================================================================
## Check for Berkeley DB
# =============================================================================
use_berkeleydb = False
# =============================================================================
## make a try to use berkeleydb
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import berkeleydb
    use_berkeleydb       = True    
    berkeleydb_open_mode = {
        'r' : berkeleydb.db.DB_RDONLY ,
        'w' : 0                       ,
        'c' : berkeleydb.db.DB_CREATE , 
        'n' : berkeleydb.db.DB_CREATE
    }
    ## open Berkeley DB 
    def berkeleydb_open ( filename                          ,
                          flag     = 'c'                    ,
                          mode     = 0o660                  ,
                          filetype = berkeleydb.db.DB_HASH  ,
                          dbenv    = None                   ,
                          dbname   = None                   ,
                          decode   = lambda s : s           ,
                          encode   = lambda s : s           ) :            
        """ Open Berkeley DB
        """
        assert flag in berkeleydb_open_mode, \
            "berkeleydb_open: invalid open mode %s" % flag
        
        if 'n' == flag and os.path.exists ( filename ) and os.path.isfile ( filename ) : 
            try    : os.remove ( filename )
            except : pass
            
        db = berkeleydb.db.DB ( dbenv )
        db.open ( filename , dbname , filetype , berkeleydb_open_mode [ flag ]  , mode )
        
        return db
    
    db_berkeley_types = 'berkeleydb' , 'berkeley' , 'berkeley-db' , 'dbhash'
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    berkeleydb      = None
    use_berkeleydb  = False 
# =============================================================================
## Check for BSDDB3 
# =============================================================================
use_bsddb3     = False
# =============================================================================
## make a try for dbddb3 
if sys.version_info < ( 3 , 10 ) :
    # ==a=======================================================================
    try : # ===================================================================
        # =====================================================================
        import bsddb3
        ## open bsddb3 database 
        def bsddb3_open ( filename        ,
                          flag    = 'c'   ,
                          mode    = 0o660 , **kwargs ) :
            """ Open `bsddb3` database
            """
            return bsddb3.hashopen ( filename , flag = flag  , mode = mode , **kwargs )        
        use_bsddb3  = True
        db_bsddb_types = 'bsddb3' , 'berkeleydb' , 'berkeley' , 'berkeley-db' , 'dbhash' 
        # =====================================================================
    except ImportError  : # ===================================================
        # =====================================================================
        bsddb3      = None 
        use_bsddb3  = False 
# =============================================================================
import dbm                    as std_db
std_whichdb = std_db.whichdb
# =============================================================================
## available DB types: 
available_backends = set ( db_std_types ) 
if use_berkeleydb : available_backends.update ( db_berkeley_types ) 
if use_sqlite3    : available_backends.update ( db_sqlite3_types  ) 
if db_dumb        : available_backends.update ( db_dumb_types     ) 
if db_gnu         : available_backends.update ( db_gnu_types      ) 
if db_dbm         : available_backends.update ( db_ndbm_types     ) 
if use_bsddb3     : available_backends.update ( db_bsddb3_types   ) 
available_backends = frozenset ( available_backends )
# ============================================================================
## check the environment variable for the preferred DBASE backends 
from ostap.utils.env import get_env, OSTAP_DBTYPES
preferrable_backends = tuple ( t for t in get_env ( OSTAP_DBTYPES , '' ).lower().split ( os.pathsep ) if t in available_backends ) 
# =============================================================================

# =============================================================================
##  Guess which db package to use to open a db file.
#  
#   Return values:
#  - None if the database file can't be read;
#  - empty string if the file can be read but can't be recognized
#  - the name of the dbm submodule (e.g. "ndbm" or "gnu") if recognized.
#   
# Importing the given module may still fail, and opening the
# database using that module may still fail.
# 
#  - Actually it is a bit extended  form of <code>dbm.whichdb</code>
#   that accounts for  <code>bsddb3</code> and <code>sqlite3</code>
def whichdb ( filename  ) :
    """ Guess which db package to use to open a db file.
    
    Return values:
    
    - None if the database file can't be read;
    - empty string if the file can be read but can't be recognized
    - the name of the dbm submodule (e.g. 'ndbm' or 'gnu') if recognized.
    
    Importing the given module may still fail, and opening the
    database using that module may still fail.
    
    - Actually it is a bit extended  form of `dbm.whichdb`
    that accounts for `bsddb3` and `sqlite3`
    """

    ## use the standard function 
    tst = std_whichdb ( filename  )

    ## dbase is identified 
    if tst : return tst 

    ## non-existing DB  ? 
    if tst is None     : return tst
    
    ## sqlite3 ?
    if issqlite3 ( filename ) : return 'sqlite3'

    import io , struct

    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        with io.open ( filename  ,'rb' ) as f :
            # Read the start of the file -- the magic number
            s16 = f.read(16)
        # =====================================================================
    except OSError : # ========================================================
        # =====================================================================
        return None

    ## TKRZW: Hash, Skip & Tree  
    if s16[:9] == b'TkrzwHDB\n': return 'TkrzwHDB' ## Hash
    if s16[:9] == b'TkrzwSDB\n': return 'TkrzwSDB' ## Skip 
    if s16[:4] == b'TDB\n'     : return 'TDB'      ## Tree 

    s = s16[:4]
    # Return "" if not at least 4 bytes
    if len ( s ) != 4:
        return ""

    ## ROOT file 
    if s == b'root'  :
        return 'root'

    # Convert to 4-byte int in native byte order -- return "" if impossible
    # =========================================================================
    try: # ====================================================================
        # =====================================================================
        ( magic, ) = struct.unpack("=l", s)
        # =====================================================================
    except struct.error: # ====================================================
        # =====================================================================
        return ""

    # =========================================================================
    # Check for GNU dbm
    if magic in (0x13579ace, 0x13579acd, 0x13579acf):
        return "dbm.gnu"

    # =========================================================================
    # Check for old Berkeley db hash file format v2
    if magic in ( 0x00061561 , 0x61150600 ):
        return "bsddb185"
    
    # =========================================================================
    # Later versions of Berkeley db hash file have a 12-byte pad in
    # front of the file type
    # =========================================================================
    try: # ====================================================================
        # =====================================================================
        ( magic , ) = struct.unpack("=l", s16[-4:])
        # =====================================================================
    except struct.error: # ====================================================
        # =====================================================================
        return ""

    # =========================================================================
    # Check for BSD hash
    if magic in ( 0x00061561 , 0x61150600 ):
        return "berkeleydb" if use_berkeleydb else "bsddb3"

    ## unknown 
    return ""

# ============================================================================
## another name for <code>whichdb</code>
isdbase = whichdb

# ============================================================================
## Open or create database at path given by *file*.
# 
#  Optional argument *flag* can be 'r' (default) for read-only access, 'w'
#  for read-write access of an existing database, 'c' for read-write access
#  to a new or existing database, and 'n' for read-write access to a new
#  database.
# 
#  Note: 'r' and 'w' fail if the database doesn't exist; 'c' creates it
#  only if it doesn't exist; and 'n' always creates a new database.
# 
#  - Actually it is a bit extended  form of <code>dbm.open</code>, that
#    accounts for <code>bsbdb3</code>, <code>sqlite3</code> 
def dbopen ( file               ,
             flag       = 'r'   ,
             mode       = 0o666 ,
             concurrent = True  ,
             dbtype     = ()    , ## preferred db-type as list of preferences  
             **kwargs           ) :
    """ Open or create database at path given by *file*.
    
    Optional argument *flag* can be 'r' (default) for read-only access, 'w'
    for read-write access of an existing database, 'c' for read-write access
    to a new or existing database, and 'n' for read-write access to a new
    database.
    
    Note: 'r' and 'w' fail if the database doesn't exist; 'c' creates it
    only if it doesn't exist; and 'n' always creates a new database.
    
    - Actually it is a bit extended  form of `dbm.open` that  accounts for `bsddb3`,`sqlite3`,'berkeleydb' 
    """

    if 'n' in flag and os.path.exists ( file ) and os.path.isfile ( file ) :
        os.unlink ( file )
        
    check = whichdb ( file ) if 'n' not in flag  else None

    if 'c' in flag and '' == check :
        check = None 
        if os.path.exists ( file ) and os.path.isfile ( file ) :
            # =====================================================================================
            try : # ===============================================================================
                os.unlink ( file )
                # =================================================================================
            except OSError : # ====================================================================
                # =================================================================================
                logger.warning ( "dbopen: unable to unlink file:`%s'" % file , exc_info = True ) 

    # 'n' flag is specified  or dbase does not exist and c flag is specified 
    if 'n' in flag or ( check is None and 'c' in flag ) : 
        
        if isinstance ( dbtype , str ) : db_types = dbtype.lower() ,  
        elif not dbtype                : db_types = () 
        else                           : db_types = tuple ( db.lower() for db in dbtype ) 

        
        def _extra_args_ ( **kw ) :
            if not kw : return
            from ostap.logger.utils import print_args
            title  = 'dbopen: Unused %d arguments' % len ( kw )
            table  = print_args ( prefix = '# ' , **kw )
            logger.warning ( '%s:\n%s' % ( title , table ) )
            
        the_dbs = dbtype if dbtype else preferrable_backends 
        
        ## check the preferred database type:
        for db in the_dbs :

            if   use_berkeleydb and ( db in db_berkeley_types or not db ) : 
                return berkeleydb_open ( file            , flag , mode , **kwargs )            
            elif SQLiteDict     and ( db in db_sqlite3_types  or not db ) : ## NB!!
                return SQLiteDict      ( filename = file , flag = flag , **kwargs )                        
            elif use_sqlite3    and ( db in db_sqlite3_types  or not db ) : ## NB!!
                _extra_args_ ( **kwargs )
                return sqlite3_open ( file , flag , mode )
            elif use_bsddb3     and db in db_bsddb3_types :  
                return bsddb3_open     ( file            , flag , mode , **kwargs ) 
            elif db_gnu  and db in db_gnu_types  : 
                _extra_args_ ( **kwargs )
                return db_gnu.open ( file , flag , mode )
            elif db_dbm  and db in db_ndbm_types : 
                _extra_args_ ( **kwargs )
                return db_dbm.open ( file , flag , mode )
            elif db in db_dumb_types : 
                _extra_args_ ( **kwargs )
                return db_dumb.open ( file , flag , mode )
            
            elif db in  db_std_types or not db :                                     ## NB !! 
                _extra_args_ ( **kwargs )
                return std_db.open ( file , flag , mode )

            elif db == 'root' : 
                from ostap.io.root_file import root_open                
                return root_open ( file , mode = flag , **kwargs ) 

        if the_dbs : logger.warning  ( 'DB-type hints not used: [%s]' %  (  ','.join ( db for fn in the_dbs ) ) ) 

        ## require the concurrent database 
        if concurrent :
            ## 
            if   use_berkeleydb : return berkeleydb_open ( file            , flag , mode , **kwargs ) 
            elif SQLiteDict     : return SQLiteDict      ( filename = file , flag = flag , **kwargs )
            elif use_sqlite3    :
                _extra_args_ ( **kwargs )
                return sqlite3_open ( file , flag , mode )
            elif use_bsddb3     : return bsddb3_open     ( file            , flag , mode , **kwargs )
            ## 
            logger.warning ( "No concurrent DBASE-backend is available" )
            
        _extra_args_ ( **kwargs )
        return std_db.open ( file , flag , mode ) 

    if use_berkeleydb and check in db_berkeley_types :
        return berkeleydb_open ( file , flag , mode , **kwargs ) 

    if SQLiteDict     and check in db_sqlite3_types :
        return SQLiteDict ( filename = file , flag = flag , **kwargs )

    if use_bsddb3     and check in db_bsddb3_types :
        return bsddb3.hashopen ( file , flag , mode , **kwargs ) 
    
    _extra_args_ ( **kwargs )

    ## dbm.sqlite3 
    if use_sqlite3    and check in db_sqlite3_types :
        return sqlite3_open ( file , flag , mode )

    if check == 'root' : 
        from ostap.io.root_file import root_open                
        return root_open ( file , mode = flag , **kwargs ) 

    ## as a last resort - use the standard stuff 
    return std_db.open ( file , flag , mode )  
    
# =============================================================================
## get disk size of data-base-like object
#  @code
#  num, size = dbsize ( 'mydb' ) 
#  @endcode  
def dbsize  ( filename  ) :
    """ Get disk  size of data-base-like object
    >>> num, size = dbsize ( 'mydb' ) 
    """
    size = 0
    num  = 0

    if os.path.exists ( filename  ) and os.path.isfile ( filename   ) :        
        size += os.path.getsize ( filename  )
        num  += 1
        
    for suffix in ( '.db'  ,
                    '.dir' , '.pag' ,
                    '.bak' , '.dir' , '.dat' ) :
        nfile = filename + suffix 
        if os.path.exists (  nfile ) and os.path.isfile ( nfile ) :
            size += os.path.getsize ( nfile  )
            num  += 1
            
    return num, size 

# ============================================================================
## Expected DB file names for the given basename 
def dbfiles ( dbtype , basename ) :
    """ Expected DB file names for the given basename
    """
    if   dbtype in ( 'dbm.ndbm'    , ) : 
        return '%s.pag' % basename , '%s.dir' % basename , 
    elif dbtype in ( 'dbm.dumb'    , ) : 
        return '%s.dat' % basename , '%s.dir' % basename , 
    elif dbtype in ( 'dbm.sqlite3' , ) : 
        return basename , ## '%s-wal' % basename , '%s-shm' % basename , 
    else :
        return basename ,

# ============================================================================
## @class TmpDB
#  Mixin class for temporary databases
# - remove : remove the temporary file immediately (just after `close')
# - keep   : keep the file and do not delete it
class TmpDB(object) :
    """ Mixin class for temporary databases
    - remove : remove the temporary file immediately (just after `close')
    - keep   : keep the file and do not delete it
    """
    def __init__ ( self            ,
                   suffix          ,
                   remove  = True  ,
                   keep    = False ) :
        
        self.__keep = True if keep  else False
        
        ## create temporary file name 
        import ostap.utils.cleanup as CU 
        fname = CU.CleanUp.tempfile ( prefix = 'ostap-tmpdb-' ,
                                      suffix = suffix         ,
                                      keep   = self.keep      )
        
        self.__tmp_name = fname        
        self.__remove   = True if ( remove and not self.keep ) else False 
        
    @property
    def tmp_name ( self ) :
        """`tmp_name' : get the generated temporary file name
        """
        return self.__tmp_name
    
    @property
    def remove ( self ) :
        """`remove':  remove the temporary file immediately (just after `close'),
        otherwise remove it at the shutdown
        """
        return self.__remove
    
    @property
    def keep   ( self )  :
        """`keep': keep the file and do not delete it
        """
        return self.__keep 
    
    ## remove the file 
    def clean  ( self ) :
        """ remove the file
        """
        fname = self.nominal_dbname 
        if self.remove and os.path.exists ( fname ) :
            import ostap.utils.cleanup as CU
            CU.CleanUp.remove_file ( fname ) 
            
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    logger.info  ('Available DB-backends are:' )
    if use_berkeleydb : logger.info ( ' - BerkeleyDB : %s' % str ( berkeleydb  ) ) 
    if SQLiteDict     : logger.info ( ' - SQLITEDICT : %s' % str ( ostap.io.sqlitedict ) )
    if use_sqlite3    : logger.info ( ' - SQLITE3    : %s' % str ( dbm.sqlite3 ) ) 
    if db_gnu         : logger.info ( ' - GNU DB     : %s' % str ( db_gnu      ) ) 
    if db_dbm         : logger.info ( ' - NDBM       : %s' % str ( db_dbm      ) ) 
    if db_dumb        : logger.info ( ' - DUMB       : %s' % str ( db_dumb     ) ) 
    if use_bsddb3     : logger.info ( ' - BSDDB3     : %s' % str ( bsddb3      ) ) 

    logger.info  ('Available   DB-backends are: %s' % ','.join ( b for b in available_backends   ) )
    logger.info  ('Preferrable DB-backends are: %s' % ','.join ( b for b in preferrable_backends ) )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
        
        
        
        
    
    
    
    
