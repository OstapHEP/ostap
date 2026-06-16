#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Simplistic lightweighted  class that provides  wrapping to duckdb data base
# @see https://duckdb.org/
"""
Simplistic lightweighted  class that provides  wrapping to duckdb data base
- see https://duckdb.org/
"""
# =============================================================================
__all__ = (
    'DuckDBDict' , ## duckdb-persistent dictionary
    'isduckdb'   , ## is it a duckdb file?
)
# =============================================================================
from   collections.abc import MutableMapping
import os, io, shutil 
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

# =============================================================================
## @class Connect
#  Helper class to implement "read-only" access to database 
class Connect ( object ) :
    """ Helper class to implement `read-only' access to database 
    """
    def __init__ ( self , filename , flag , **config ) :
        
        self.__filename = filename
        self.__flag     = flag
        self.__config   = config
        self.__connect  = None

    # =========================================================================
    ## context manager ENTER
    def __enter__ ( self ) :

        self.__connect = duckdb.connect ( self.filename , 
                                         read_only = 'r'in self.flag , 
                                         **self.config )
                    
        return self.connect

    # ========================================================================
    ## context manager EXIT 
    def __exit__ ( self , *_ ) :

        if self.connect :
            self.__connect.close()
            self.__connect = None
            
    @property
    def filename ( self ) :
        """`filename' : the original filename"""
        return self.__filename

    @property
    def flag ( self ) :
        """`flag' : the access flag"""
        return self.__flag

    @property
    def config( self ) :
        """`connfig` : addtional keyword arguments"""
        return self.__config
    
    @property
    def connect ( self ) :
        """`connect' : the actual connection"""
        return self.__connect
    
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
   
    def __init__ ( self                ,
                   filename            ,
                   flag      = 'r'     , 
                   tablename = 'ostap' , **config )  :
        
        ## check presence of the tkrzwz module 
        assert duckdb , "`duckdb` module is not available!"
        
        assert 1 == len ( flag ) and isinstance ( flag , str ) and flag.lower() in sefl.FLAGS , "Invalid `flag`:%s" % flag
        
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
        
        ## check it! 
        with Connect ( self.filename , self.flag , **config ) :
            logger.debug ( "Opening DuckDB table %r in %s" % ( tablename , filename ) )
        
        self.__filename  = filename
        self.__config    = config
        self.__tablename = tablename 
        self.__flag      = flag 
        
 
    @property
    def db  ( self ) :
        """`db` : the actual Tkrzw database ( same ad `dbm`) """
        return self.db 

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
    ## Get data from DB
    def get ( self , key , default = None ) :
        """ Get data from DB
        >>> db = ...
        >>> value = db [ key ]
        """
        assert not self.db is None , "TkrzwDB is invalid!"
        value = self.db.Get ( key )
        return default if value is None else self.decode_value ( value ) 
    
    # =========================================================================
    ## Get data from DB
    def __getitem__( self , key ) :
        """ Get data from DB
        >>> db = ...
        >>> value = db [ key ]
        """
        value = self.get ( key ) 
        if value is None: raise KeyError( 'Invalid key:%s' % key )        
        return value 

    # =========================================================================
    ## Add data to DB
    def __setitem__ ( self , key , value ) :
        """ Add data to DB 
        >>> db = ...
        >>> db [ key ] = value 
        """
        assert not self.db is None , "TkrzwDB is invalid!"
        vv = self.encode_value ( value )
        self.db.Set ( key , vv  )  

    # =========================================================================
    ## iterate over key,value pairs 
    def items ( self ) :
        """ Iterate over key, value pairs
        >>> db = ...
        >>> for key, value in db.items () : ...
        """
        assert not self.db is None , "TkrzwDB is invalid!"
        for k, value in self.db :
            yield key, self.decode_value ( value )

    # =========================================================================
    ## iterate over key,value pairs 
    def iteritems ( self ) :
        """ Iterate over key, value pairs
        >>> db = ...
        >>> for key, value in db.iteritems () : ...
        """
        return self.items() 

    # =========================================================================
    ## iterator over keys 
    def keys ( self ) :
        """ Iterator over keys
        >>>  db = ...
        >>>  for key in db.keys() : .... 
        """
        assert not self.db is None , "TkrzwDB is invalid!"
        for k, _ in self.db : yield key

    # ========================================================================
    ## Iterator over values
    def values ( self ) :
        """ Iterator over values
        >>> db = ...
        >>> for value in db.values() : ...
        """
        assert not self.db is None , "TkrzwDB is invalid!"
        for _ , value in self.db : yield self.decode_value ( value )

    # =========================================================================
    ## check the key in database 
    def __contains__ ( self , key ) :
        """ Check the key in database """
        assert not self.db is None , "TkrzwDB is invalid!"
        value = self.get ( key ) 
        return value is not None

    # =========================================================================
    ## iterator over keys 
    def __iter__ ( self ) :
        """ Iterator over keys
        >>> db = ...
        >>> for key in db : ... 
        """
        return self.keys ()

    # =========================================================================
    ## number of stored leentskeys in the database 
    def __len__ ( self ) :
        """ Number of stored elements in the database
        >>> db = ...
        >>> nu m= len ( db ) 
        """
        assert not self.db is None , "TkrzwDB is invalid!"
        return len ( self.db ) 

    # =======================================================================--
    ## pop the certain key from the database 
    def pop ( self, key , default = None ) :
        """ pop the certain key from database
        >>> db = ...
        >>> value = db.pop ( key ) 
        """
        assert not self.db is None , "TkrzwDB is invalid!"
        sc , value = self.db.RemoveAndGet( key )
        if sc.IsOK () and not value is None : return self.decode_value ( value )
        return default 

    # =========================================================================
    ## delete certain key from the database 
    def __delitem__( self , key ) :
        """ Delete certain key from the database
        >>> db = ...
        >>> del db[ key ]
        """
        assert not self.db is None , "TkrzwDB is invalid!"
        sc = self.db.Remove ( key )
        
    # ========================================================================
    ## sync database 
    def sync ( self ) :
        """ Sync database
        >>> db = ...
        >>> db.sync () 
        """
        assert not self.db is None , "TkrzwDB is invalid!"
        self.db.Synchronize ( True )

    # =========================================================================
    ## close the database 
    def close ( self ) :
        """ Close database 
        >>> db = ...
        >>> db.close() 
        """
        if self.db is None : return
        
        if self.db.IsWritable () and self.db.ShouldBeRebuild() :
            self.db.Rebuild()
            
        self.db.Close ()
        self.__db = None 

    # =========================================================================
    ## context manager: ENTER 
    def __enter__ ( self ) :
        """ Context manager: ENTER """
        return self
    
    # =========================================================================
    ## context manager: EXIT  
    def __exit__ ( self , *args ) :
        """ Context manager: EXIT """        
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
