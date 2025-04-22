#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/io/dbase.py
# 
# Helper module to use databases
#
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
    'whichdb'    , ## guess database type
    'isdbase'    , ## ditto 
    'dbopen'     , ## open database
    'Item'       , ## item: named tuple (time,payload)
    'TmpDB'      , ## mixing for tempoirary database
    
    )
# =============================================================================
import sys, os, collections
from   ostap.logger.logger  import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.dbase' )
else                      : logger = getLogger ( __name__         )
# =============================================================================
from ostap.io.sqlitedict import SqliteDict, issqlite3
# =============================================================================
## named tuple to DB-item: (time, payload)
Item = collections.namedtuple ( 'Item', ( 'time' , 'payload' ) )
# =============================================================================
import dbm                    as std_db
std_whichdb = std_db.whichdb
# =============================================================================
ordered_dict = dict
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import dbm.gnu  as db_gnu
    # =========================================================================
except ImportError :
    # =========================================================================
    db_gnu = None
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import dbm.ndbm as db_dbm
    # =========================================================================
except ImportError :
    # =========================================================================
    db_dbm = None
# =============================================================================
db_hash = None
# =============================================================================
import dbm.dumb as db_dumb
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
    use_berkeleydb   = True
    
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
            """ Open `bsddb3` database """
            return bsddb3.hashopen ( filename , flag = flag  , mode = mode , **kwargs )        
        use_bsddb3  = True
        # =====================================================================
    except ImportError  : # ===================================================
        # =====================================================================
        bsddb3      = None 
        use_bsddb3  = False 

# =============================================================================
## make a try to use LMDB
use_lmdb = False
# =============================================================================
## make a try for LMDB 
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import lmdb 
    from ostap.io.lmdbdict import LmdbDict, islmdb 
    use_lmdb = True
    # =========================================================================
except ImportError  : # =======================================================
    # =========================================================================
    lmdb     = None 
    use_lmdb = False 

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

    ## make a try with LMDB 
    if use_lmdb and os.path.exists  ( filename ) and os.path.isdir ( filename ) :
        if islmdb ( filename ) : return 'lmdb'
    
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

    ## TKRZW: Hash 
    if s16[:9] == b'TkrzwHDB\n': return 'TkrzwHDB'

    s = s16[:4]
    # Return "" if not at least 4 bytes
    if len ( s ) != 4:
        return ""

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

    # Check for GNU dbm
    if magic in (0x13579ace, 0x13579acd, 0x13579acf):
        return "dbm.gnu"

    # Check for old Berkeley db hash file format v2
    if magic in ( 0x00061561 , 0x61150600 ):
        return "bsddb185"

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
#    accounts for <code>bsbdb3</code>, <code>sqlite3</code> and <code>lmdbdict</code>
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
    
    - Actually it is a bit extended  form of `dbm.open` that  accounts for `bsddb3`,`sqlite3`,'berkeleydb' and `lmdb`
    """

    if 'n' in flag and os.path.exists ( file ) and os.path.isfile ( file ) :
        os.unlink ( file )
        
    check = whichdb ( file ) if 'n' not in flag  else None

    if 'c' in flag and '' == check :
        check = None 
        if os.path.exists ( file ) and os.path.isfile ( file ) : os.unlink ( file ) 

    # 'n' flag is specified  or dbase does not exist and c flag is specified 
    if 'n' in flag or ( check is None and 'c' in flag ) : 
        
        if isinstance ( dbtype , str ) : db_types = dbtype.lower() ,  
        elif not dbtype                : db_types = () 
        else                           : db_types = tuple ( db.lower() for db in dbtype ) 

        if kwargs : message = 'Ignore extra %d arguments:%s' % ( len ( kwargs ) , [ k for k in kwargs ] )
        else      : message = '' 
        
        ## check the preferred database type:
        for db in db_types :
                        
            if db             in ( 'sqlite3' , 'sqlite'  , 'sql'                      , '' ) : ## NB!! 
                return SqliteDict      ( filename = file , flag = flag , **kwargs )            
            elif use_berkeleydb and db in ( 'berkeleydb' , 'berkeley' , 'berkeley-db' , '' ) :  ## NB!!
                return berkeleydb_open ( file            , flag , mode , **kwargs ) 
            elif use_bsddb3     and 'bsddb3'     == db :
                return bsddb3_open     ( file            , flag , mode , **kwargs ) 
            elif use_lmdb       and 'lmdb'       == db :
                return LmdbDict        ( path     = file , flag = flag , **kwargs )
            elif db_gnu  and db in ( 'dbm.gnu'  , ) :
                if kwargs : logger.warning ( message ) 
                return db_gnu.open ( file , flag , mode )
            elif db_dbm  and db in ( 'dbm.ndbm' , ) :
                if kwargs : logger.warning ( message  ) 
                return db_dbm.open ( file , flag , mode )
            elif db_hash and db in ( 'dbhash' , ) :
                if kwargs : logger.warning ( message ) 
                return db_hash.open ( file , flag , mode )
            elif db in ( 'dbm.dumb' , 'dumbdbm' , 'dumb' ) :
                if kwargs : logger.warning ( message ) 
                return db_dumb.open ( file , flag , mode )
            elif db in  ( 'std' , 'standard' ) or not db :                                     ## NB !! 
                if kwargs : logger.warning ( message )                 
                return std_db.open ( file , flag , mode )

        if db_types :
            logger.warning  ( 'DB-type hints not used: [%s]' %  (  ','.join ( db for fn in db_types ) ) ) 
        
        if concurrent and use_berkeleydb :
            return berkeleydb_open ( file , flag , mode , **kwargs ) 

        if concurrent and use_bsddb3     :
            return bsddb3_open     ( file , flag , mode , **kwargs ) 

        if concurrent :
            return SqliteDict      ( filename = file , flag = flag , **kwargs )

        if kwargs : logger.warning ( message ) 
        return std_db.open ( file , flag , mode ) 

    if use_berkeleydb and check in ( 'berkeleydb' , 'bsddb3' , 'dbhash' ) :
        return berkeleydb_open ( file , flag , mode , **kwargs ) 

    if use_bsddb3     and check in ( 'berkeleydb' , 'bsddb3' , 'bsddb' , 'dbhash' , 'bsddb185' ) :
        return bsddb3.hashopen ( file , flag , mode , **kwargs ) 

    if use_lmdb       and check in ( 'lmdb' , ) :
        return LmdbDict    ( path    = file , flag = flag , **kwargs )
    
    if check in ( 'sqlite3' , 'sqlite' ) :
        return SqliteDict ( filename = file , flag = flag , **kwargs )

    if kwargs : logger.warning ( 'Ignore extra %d arguments:%s' % ( len ( kwargs ) , [ k for k in kwargs ] ) ) 

    ## as a lasty resort - use the standard stuff 
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
    if   dbtype in ( 'dbm.ndbm' , ) : 
        return '%s.pag' % basename , '%s.dir' % basename , 
    elif dbtype in ( 'dbm.dumb' , ) : 
        return '%s.dat' % basename , '%s.dir' % basename , 
    elif dbtype in ( 'lmdb', ) : 
        return ( os.path.join ( basename , ''         ) , ## directory 
                 os.path.join ( basename , 'data.mdb' ) ,
                 os.path.join ( basename , 'lock.mdb' ) )
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
    if use_berkeleydb : logger.info ( ' - BerkeleyDB : %s' % str ( berkeleydb ) ) 
    if use_bsddb3     : logger.info ( ' - BSDDB3     : %s' % str ( bsddb3     ) ) 
    if use_lmdb       : logger.info ( ' - LMDB       : %s' % str ( lmdb       ) ) 
    if db_gnu         : logger.info ( ' - GNU dbase  : %s' % str ( db_gnu     ) ) 
    if db_dbm         : logger.info ( ' - NDBM       : %s' % str ( db_dbm     ) ) 
    if db_dumb        : logger.info ( ' - DUMB       : %s' % str ( db_dumb    ) ) 
    if db_hash        : logger.info ( ' - DBHASH     : %s' % str ( db_hash    ) ) 
        
# =============================================================================
##                                                                      The END 
# =============================================================================
        
        
        
        
    
    
    
    
