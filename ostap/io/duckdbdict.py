#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Simplistic lightweighted  class that provides  wrapping to duckdb data base
# @see https://duckdb.org/
# @author Vanya BELYAEV 2024-06-17
# To provide some concurrency control, the `fasteners` module is used with explicit Read/Write locks.
# @see https://pypi.org/project/fasteners/
# It results in some sizeable performance penalty, but provides reasonaly safety, 
# The penalty is negligible for large dbase&and large records, 
# but can be sizeable for small dbase and small records.
# =============================================================================
"""
Simplistic lightweighted  class that provides  wrapping to duckdb data base
- see https://duckdb.org/

To provide some concurrency control, the `fasteners` module is used with explicit Read/Write locks.
- see https://pypi.org/project/fasteners/
It results in some (sizeable) performance penalty, but provides reasonaly safety, 
The penalty is negligible for large dbase&and large records, 
#but can be sizeable for small dbase and small records.
"""
# =============================================================================
__all__ = (
    'DuckDBDict' , ## duckdb-persistent dictionary
    'isduckdb'   , ## is it a duckdb file?
)
# =============================================================================
from   collections.abc         import MutableMapping
from   pathlib                 import Path
from   ostap.core.ostap_types  import dictlike_types 
import ostap.utils.cleanup     as     CU
import os, io, shutil, hashlib 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.io.duckdbdict' )
else                       : logger = getLogger( __name__     )
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import duckdb
    # =========================================================================
except ImportError: # =========================================================
    # =========================================================================
    logger.warning ("No `duckdb` module is available!")
    duckdb = None
# =============================================================================
RWFileLock = None
# =============================================================================
if duckdb : # =================================================================  
    # =========================================================================    
    try : # ===================================================================
        # =====================================================================
        from fasteners import InterProcessReaderWriterLock as RWFileLock
        # =====================================================================
    except ImportError : # ====================================================
        logger.warning ("No `fasteners` module is available!")
       
# =============================================================================
## Is it a duckdb file?
#  @code
#  ok = isduckdb ( 'mydbase' ) 
#  @endcode
def isduckdb ( filename ) :
    """ Is it a duckdb file?
    >>> ok = isduckdb ( 'mydbase' ) 
    """
    if not   os.path.exists  ( filename ) : return False
    if not   os.path.isfile  ( filename ) : return False
    if 100 > os.path.getsize ( filename ) : return False
    
    with io.open ( filename  , 'rb' ) as f :
        hdr = f.read(12)
        return hdr[8:12] == b'DUCK'
    
    return False

# ============================================================================
## directory for lock-file for duckdb database
lock_dir = Path ( CU.CleanUp.tempdir ( prefix = 'ostap-duckdb-locks' , date = False) ) 

# ============================================================================
## Lock-file for duckdb database
def lock_file ( filename ) :
    """ Lock-file for duckdb database
    >>> lock = lock_file ( 'mydbase' ) 
    """
    assert filename, "Invalid file name %s" % filename
    path      = Path ( filename ).resolve() 
    name      = path.name 
    hash_path = hashlib.md5(str(path).encode('utf-8')).hexdigest()
    return str ( lock_dir / f"{hash_path}.{name}.rwlock" )

# =============================================================================
## @class DuckDBDict
#  Persistent dictionary based on duckdb
#  @attention both keys and values are bytestrings!
#  @code
#
#  mydict = DuckDBDict ( 'some.db' , 'c' )
#
#  mydict [ key ] = value
#
#  for k   in mydict          : print ( k )
#  for k   in mydict.keys()   : print ( k )
#  for v   in mydict.values() : print ( v )
#  for k,v in mydict.items () : print ( k , v )
#
#  val  = mdict[ key]
#
#  del mydict [ key ]
#
# @endcode 
class DuckDBDict ( MutableMapping ) :
    """ Persistent dictionary based on duckdb 
    - both keys and values are bytestrings! 

    >>> mydict = DuckDBDict ( 'some.db' , 'c' )
    >>> mydict [ key ] = value
    
    >>> for k   in mydict          : print ( k )
    >>> for k   in mydict.keys()   : print ( k )
    >>> for v   in mydict.values() : print ( v )
    >>> for k,v in mydict.items () : print ( k , v )

    >>> val  = mdict[ key]

    >>> del mydict [ key ]
    
    """
    FLAGS      = 'rwcn'
   
    def __init__ ( self                  ,
                   filename              ,
                   flag        = 'r'     , 
                   tablename   = 'ostap' , 
                   keyencoding = 'utf-8' , **config )  :
        
        ## check presence of the tkrzwz module 
        assert duckdb     , "`duckdb` module is not available!"
        assert RWFileLock , "`fasteners` module is not available!"
        
        assert 1 == len ( flag ) and isinstance ( flag , str ) and flag.lower() in self.FLAGS , "Invalid `flag`:%s" % flag
        
        flag = flag.lower() 

        if 'n' == flag and filename and os.path.exists ( filename ) :
            if   os.path.isfile ( filename ) :
                try    : os.unlink ( filename )
                except : pass 
            elif os.path.isdir  ( filename ) :
                try    : shutil.rmtree ( filename )
                except : pass                    
            if os.path.exists ( filename ) :
                logger.warning ( 'path exists!' )
         
        self.__filename    = filename
        self.__config      = config
        self.__tablename   = tablename 
        self.__flag        = flag 
        self.__keyencoding = keyencoding   
        
        ## Read-Write Lock 
        lf = lock_file ( self.filename )
        self.__lock = RWFileLock( lf )
        
        ## check if the table exists 
        if not 'r' in self.flag :
            table_exists = False 
            if os.path.exists ( self.filename ) :
                with self.__lock.read_lock () : 
                    with duckdb.connect ( self.filename , read_only = True ) as conn:
                        query = "SELECT EXISTS (SELECT 1 FROM duckdb_tables WHERE table_name = ?)"
                        table_exists = conn.execute(query, [ self.tablename]).fetchone()[0]
            if not table_exists : 
                with self.__lock.write_lock():
                    with duckdb.connect(self.__filename, read_only=False) as conn:
                        conn.execute(f"""CREATE TABLE IF NOT EXISTS {self.tablename} (
                                        key VARCHAR PRIMARY KEY, 
                                        value BLOB )
                                     """)
                    
        self.__SQL_GET  =  f"SELECT value FROM {self.tablename} WHERE key = ?"
        self.__SQL_SET  =  f"""INSERT INTO {self.tablename} (key, value) VALUES (?, ?) 
                               ON CONFLICT (key) DO UPDATE SET value = EXCLUDED.value"""
        self.__SQL_DEL  =  f"DELETE FROM {self.tablename} WHERE key = ?"
        self.__SQL_IN   =  f"SELECT 1 FROM {self.tablename} WHERE key = ?"
        self.__SQL_ITER =  f"SELECT key, value FROM {self.tablename}"
        self.__SQL_LEN  =  f"SELECT COUNT(*) FROM {self.tablename}" 
        
    @property
    def filename ( self ) : 
        """Get the filename of the duckdb database."""
        return self.__filename
    
    @property
    def flag ( self ) :
        """Get the flag used to open the duckdb database."""
        return self.__flag
    
    @property
    def tablename ( self ) :
        """Get the name of the table used in the duckdb database."""
        return self.__tablename
        
    @property
    def keyencoding ( self ) :
        """Get the key encoding used for the duckdb database."""
        return self.__keyencoding

    ## Get data from DB
    def get ( self , key , default = None ) :
        """ Get data from DB 
        """
        with self.__lock.read_lock():
            with duckdb.connect(self.__filename, read_only=True) as conn:
                result = conn.execute( self.__SQL_GET , [key] ).fetchone()
                
        return default if result is None else self.decode_value ( result [ 0 ] ) 
            
    ## Get data from DB 
    def __getitem__ ( self , key ) :
        """ Get data from DB 
        """
        value = self.get ( key , default = None ) 
        if value is None : raise KeyError ( key )
        return value  
    
    ## Add data from DB 
    def __setitem__( self , key , value):
        """ Add data to DB """ 
        ## key = key.encode ('utf8')
        assert not 'r' in self.flag, "Database is opened as read-only!"
        with self.__lock.write_lock():
            with duckdb.connect(self.__filename, read_only=False) as conn:
                conn.execute( self.__SQL_SET , [ key , self.encode_value ( value ) ] )
          
    ## 
    def __delitem__ ( self , key ) :
        """ delete data from DB
        """
        assert not 'r' in self.flag, "Database is opened as read-only!"
        if key not in self : raise KeyError ( key )
        with self.__lock.write_lock():
            with duckdb.connect ( self.filename , read_only = False) as conn:
                conn.execute( self.__SQL_DEL , [ key ] )   
    
    def __contains__ ( self , key ):
        with self.__lock.read_lock():
           with duckdb.connect ( self.filename , read_only = True ) as conn:
                return conn.execute ( self.__SQL_IN , [ key ] ).fetchone() is not None 
                
    def __len__( self ):
        with self.__lock.read_lock():
            with duckdb.connect ( self.filename , read_only = True ) as conn:
                return conn.execute ( self.__SQL_LEN ).fetchone() [ 0 ]

    def __iter__ ( self ):
        with self.__lock.read_lock():
            with duckdb.connect ( self.filename , read_only = True ) as conn:
                cursor = conn.execute ( self.__SQL_ITER )
                keys   = [ row [ 0 ].encode ( self.keyencoding ) for row in cursor.fetchall() ]

        return iter ( keys )
        
    # =========================================================================
    ## Bulk update database 
    def update ( self , other = None , **kwargs ) :     
        """ Bulk update database 
        >>> db = ...
        >>> db.update ( other )
        >>> db.update ( key1 = value1 , key2 = value2 ) 
        """
        if     other is None        : pass 
        if not other and not kwargs : return 
        ##
        assert not 'r' in self.flag, "Database is opened as read-only!"
        with self.__lock.write_lock():
            with duckdb.connect(self.__filename, read_only=False) as conn :
                
                if isinstance ( other , dictlike_types ) : 
                    for key, value in other.items() :
                        conn.execute ( self.__SQL_SET , [ key , self.encode_value ( value ) ] )
                elif hasattr ( other , 'keys' ) : 
                    for key in other.keys() :
                        value = other [ key ]
                        conn.execute ( self.__SQL_SET , [ key , self.encode_value ( value ) ] )
                else :  
                    for key, value in other :
                        conn.execute ( self.__SQL_SET , [ key , self.encode_value ( value ) ] )  
                
                for key, value in kwargs.items () :
                    conn.execute ( self.__SQL_SET , [ key , self.encode_value ( value ) ] )     
                        
    # =========================================================================
    ## encode value: input -> dbase 
    def encode_value ( self , value ) : 
        """ encode value: 
        - input -> dbase 
        """
        return value

    # =========================================================================
    ## decode value: dbase -> ouput (no-op) 
    def decode_value ( self , value ) :
        """ decode value: (no-op) 
        - dbase -> output 
        """
        return value 
    
    # =========================================================================
    ## Close database : NO-OP
    def close ( self ) : 
        """ Close database : NO-OP
        """
        pass 
    
    # ========================================================================
    ## sync database: NO-OP 
    def sync ( self ) :
        """ Sync database
        >>> db = ...
        >>> db.sync () 
        """
        pass

    # =========================================================================
    ## context manager: ENTER 
    def __enter__ ( self ) :
        """ Context manager: ENTER 
        """
        return self
    
    # =========================================================================
    ## context manager: EXIT  
    def __exit__ ( self , *args ) :
        """ Context manager: EXIT 
        """        
        self.close()

# ==============================================================================
## disable imports 
if not duckdb  : __all__ = () 

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if not duckdb: logger.warning ( "No `duckdb` module is found!")

# =============================================================================
#                                                                       The END 
# =============================================================================
