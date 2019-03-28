#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file sqliteshelve.py
# 
# This is SQLite-version of shelve database with compressed content 
# 
# Keeping the same interface and functionlity as shelve data base,
# SQLiteShelve allows much more compact file size through
# the compression of the content
#
# The actual code has been taken from <c>sqlitedict</c> code
#  by Radim Rehurek <radimrehurek@seznam.cz>
# Hacked together from:
#  * http://code.activestate.com/recipes/576638-draft-for-an-sqlite3-based-dbm/
#  * http://code.activestate.com/recipes/526618/   ( see Google...)
#
# The compression (with zlib) is added atop of original code 
# 
# Create new DB:
#
# @code
#
# >>> import sqliteshelve as DBASE     ## import the SQLiteShelve module 
# >>> db = DBASE.open ('a_db', 'n')    ## create new DB
# ...
# >>> abcde = ...
# >>> db['some_key'] =  abcde                 ## add information to DB
# ...
# >>> db.close()
#
# @endcode 
#
# Access to DB in read-only mode :
#
# @code
#
# >>> import sqliteshelve as DBASE     ## import the SQLiteShelve module 
# >>> db = DBASE.open ('a_db' , 'r' ) ## access existing dbase in read-only mode
# ...
# >>> for key in db : print key
# ...
# >>> abcd = db['some_key']
#
# @endcode 
#
# Access existing DB in update mode :
#
# @code
#
# >>> import sqliteshelve as DBASE     ## import the SQLiteShelve module 
# >>> db = DBASE.open ('a_db' )      ## access existing dbase in update mode
# ...
# >>> for key in db : print key
# ...
# >>> abcd = db['some_key']
#
# @endcode 
#
# 
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date   2010-04-30
# 
# =============================================================================
"""This is SQLite-version of shelve database.

Keeping the same interface and functionlity as shelve data base,
SQLiteShelve allows much more compact file size through
the compression of the content

The actual code has been inspired by <c>sqlitedict</c> ( see Google...)
The compression (with zlib) is added atop of original code 
 
Create new DB:

# >>> import sqliteshelve as DBASE            ## import the SQLiteShelve module 
# >>> db = DBASE.open ('a_db', 'n')    ## create new DB
# ...
# >>> abcde = ...
# >>> db['some_key'] =  abcde                 ## add information to DB
# ...
# >>> db.close()

Access to DB in read-only mode :

# >>> import sqliteshelve as DBASE            ## import the SQLiteShelve module 
# >>> db = DBASE.open ('a_db' , 'r' )  ## access existing dbase in read-only mode
# ...
# >>> for key in db : print key
# ...
# >>> abcd = db['some_key']


 Access existing DB in update mode :

# >>> import sqliteshelve as DBASE            ## import the SQLiteShelve module 
# >>> db = DBASE.open ('a_db' )       ## access existing dbase in update mode
# ...
# >>> for key in db : print key
# ...
# >>> abcd = db['some_key']
"""
# =============================================================================
from   __future__        import print_function
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-19"
__version__ = "$Revision$" 
# =============================================================================
__all__ = (
    'SQLiteShelf' ,   ## The DB-itself
    'open'        ,   ## helper function to hide the actual DB
    'tmpdb'       ,   ## helper function to create the temporary database 
    )
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.sqliteshelve' )
else                      : logger = getLogger ( __name__ )
# =============================================================================
from   ostap.io.sqlitedict import SqliteDict
import sys, zlib 
import sqlite3
# =============================================================================
try:
    from cPickle   import Pickler, Unpickler , HIGHEST_PROTOCOL
except ImportError:
    from  pickle   import Pickler, Unpickler , HIGHEST_PROTOCOL
# =============================================================================
PROTOCOL = 2 
# =============================================================================
try : 
    from io        import             BytesIO
except ImportError : 
    from cStringIO import StringIO as BytesIO         
# =============================================================================
_modes_ = {
    # =========================================================================
    # 'r'	Open existing database for reading only
    # 'w'	Open existing database for reading and writing
    # 'c'	Open database for reading and writing, creating it if it doesn’t exist
    # 'n'	Always create a new, empty database, open for reading and writing
    # =========================================================================
    'n'        : 'n' ,
    'c'        : 'c' ,
    'r'        : 'r' ,
    'u'        : 'w' ,
    'w'        : 'w' ,
    'a'        : 'w' ,
    ##
    '+'        : 'w' ,        
    'w+'       : 'w' ,
    'rw'       : 'w' ,
    'new'      : 'n' ,
    'create'   : 'c' ,
    'recreate' : 'n' ,
    'read'     : 'r' ,        
    'write'    : 'w' ,        
    'update'   : 'w' ,        
    'append'   : 'w' ,        
    }
# =============================================================================
## @class SQLiteShelf
#  SQLite-based ``shelve-like'' database with compressed content.
#
#  Keeping the same interface and functionlity as shelve data base,
#  SQLiteShelve allows much more compact file size through
#  the compression of the content
#
#  The actual code has been taken from <c>sqlitedict</c> code
#  by Radim Rehurek <radimrehurek@seznam.cz>
#  Hacked together from:
#  * http://code.activestate.com/recipes/576638-draft-for-an-sqlite3-based-dbm/
#  * http://code.activestate.com/recipes/526618/   ( see Google...)
#
# The compression (with zlib) is added atop of original code
#
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date   2010-04-30
class SQLiteShelf(SqliteDict):
    """ SQLite-based ``shelve-like'' database with compressed content.

    Keeping the same interface and functionlity as shelve data base,
    SQLiteShelve allows much more compact file size through
    the compression of the content
    
    The actual code has been inspired by <c>sqlitedict</c> ( see Google...)
    The compression (with zlib) is added atop of original code 
    
    Create new DB:
    
    >>> import SQLiteShelve                     ## import the SQLiteShelve module 
    >>> db = SQLiteShelve.open ('a_db', 'n')    ## create new DB
    ...
    >>> abcde = ...
    >>> db['some_key'] =  abcde                 ## add information to DB
    ...
    >>> db.close()
    
    Access to DB in read-only mode :
    
    >>> import SQLiteShelve                     ## import the SQLiteShelve module 
    >>> db = SQLiteShelve.open ('a_db' , 'r' )  ## access existing dbase in read-only mode
    ...
    >>> for key in db : print key
    ...
    >>> abcd = db['some_key']
    
    Access existing DB in update mode :
    
    >>> import SQLiteShelve                    ## import the SQLiteShelve module 
    >>> db = SQLiteShelve.open ('a_db' )       ## access existing dbase in update mode
    ...
    >>> for key in db : print key
    ...
    >>> abcd = db['some_key']
    """
    def __init__ ( self                       ,
                   filename       = None      ,
                   mode           = 'c'       ,
                   tablename      = 'ostap'   ,
                   writeback      = True      , ## original name: "autocommit"
                   compress_level = zlib.Z_BEST_COMPRESSION , 
                   journal_mode   = "DELETE"  ,
                   protocol       = PROTOCOL  ) :
        """Initialize a thread-safe sqlite-backed dictionary.
        The dictionary will be a table ``tablename`` in database file
        ``filename``. A single file (=database) may contain multiple tables.
        
        If no ``filename`` is given, a random file in temp will be used
        (and deleted from temp once the dict is closed/deleted).
        
        If you enable ``writeback/autocommit`` changes will be committed
        after each operation (more inefficient but safer).
        Otherwise, changes are committed on
        ``self.commit()``,
        ``self.clear()`` and
        ``self.close()``.
        
        Set ``journal_mode`` to ``OFF``
        if you're experiencing sqlite I/O problems
        or if you need performance and don't care about crash-consistency.
        
        The `mode` parameter:
        - 'c': default mode, open for read/write, creating the db/table if necessary.
        - 'w': open for r/w, but drop `tablename` contents first (start with empty table)
        - 'n': create a new database (erasing any existing tables, not just `tablename`!).
        
        Modes: %s 
        """ % _modes_ 
        
        ## the mode 
        mode = _modes_.get( mode.lower() , '' )
        if not mode :
            logger.warning("Unknown opening mode '%s', replace with 'c'")
            mode = 'c'
            
        if not filename is None :
            import os 
            filename  = os.path.expandvars ( filename )
            filename  = os.path.expanduser ( filename )
            filename  = os.path.expandvars ( filename )
            filename  = os.path.abspath    ( filename )
            
        self.__filename = filename 
        
        SqliteDict.__init__ ( self                        ,
                              filename     = filename     ,
                              tablename    = tablename    ,
                              flag         = mode         ,
                              autocommit   = writeback    ,
                              journal_mode = journal_mode )
        
        self.__compression = compress_level 
        self.__protocol    = protocol

    @property
    def compression ( self ) :
        """The  compression level from zlib"""
        return self.__compression
    
    @property
    def protocol    ( self ) :
        """The pickling protocol"""
        return self.__protocol
    
    ## list the avilable keys 
    def __dir ( self , pattern = '' ) :
        """ List the avilable keys (patterns included).
        Pattern matching is performed accoriding to
        fnmatch/glob/shell rules [it is not regex!] 
        
        >>> db = ...
        >>> db.ls() ## all keys
        >>> db.ls ('*MC*')
        
        """
        keys_ = self.keys()
        if pattern :
            import fnmatch
            _keys = [ k for k in keys_ if fnmatch.fnmatchcase ( k , pattern ) ]
            keys_ = _keys
        #
        for key in sorted(keys_) : print(key)
        
    ## list the avilable keys 
    def ls    ( self , pattern = '' ) :
        """List the avilable keys (patterns included).
        Pattern matching is performed accoriding to
        fnmatch/glob/shell rules [it is not regex!] 

        >>> db = ...
        >>> db.ls() ## all keys
        >>> db.ls ('*MC*')        
        
        """
        return self.__dir( pattern )

# =============================================================================
## ``get-and-uncompress-item'' from dbase 
def _zip_getitem (self, key):
    """ ``get-and-uncompress-item'' from dbase 
    """
    GET_ITEM = 'SELECT value FROM %s WHERE key = ?' % self.tablename
    item = self.conn.select_one(GET_ITEM, (key,))
    if item is None: raise KeyError(key)

    f     = BytesIO ( zlib.decompress ( item[0] ) ) 
    value = Unpickler(f).load()

    return value

# =============================================================================
## ``set-and-compress-item'' to dbase 
def _zip_setitem ( self , key , value ) :
    """ ``set-and-compress-item'' to dbase 
    """
    ADD_ITEM = 'REPLACE INTO %s (key, value) VALUES (?,?)' % self.tablename
    
    f     = BytesIO()
    p     = Pickler(f, self.protocol )
    p.dump(value)
    blob  = f.getvalue() 
    zblob = zlib.compress ( blob , self.compression ) 
    self.conn.execute(ADD_ITEM, (key, sqlite3.Binary( zblob ) ) )

                      
SQLiteShelf.__setitem__ = _zip_setitem
SQLiteShelf.__getitem__ = _zip_getitem

def _sql_enter_ ( self      ) : return self
def _sql_exit_  ( self , *_ ) :
    try :
        os.close()
    except : pass

SQLiteShelf.__enter__ = _sql_enter_
SQLiteShelf.__exit__  = _sql_exit_ 

# =============================================================================
## add an object into data base
#  @code
#  dbase  = ...
#  object = ...
#  dbase.ls() 
#  object >> dbase ## add object into dbase 
#  dbase.ls() 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-06-04
def _db_rrshift_ ( db , obj ) :
    """Add an object into data base
    
    dbase  = ...
    object = ...
    dbase.ls() 
    object >> dbase ## add object into dbase 
    dbase.ls() 
    """
    if   hasattr ( obj , 'GetName' ) : name = obj.GetName()
    elif hasattr ( obj , 'name'    ) : name = obj.name   ()
    else : name =  obj.__class__.__name__
    #
    db [ name ] = obj
    
SQLiteShelf.__rrshift__ = _db_rrshift_

# =============================================================================
## open new SQLiteShelve data base
#  @code
#  import SQLiteShleve as DBASE
#  db = DBASE.open( 'data.msql')
#  db['a'] = ...
#  @endcode
#  @see SQLiteShelf
def open(*args, **kwargs):
    """See documentation of the SQLiteShelf class.
    >>> import SQLiteShleve as DBASE
    >>> db = DBASE.open('data.msql','c')
    >>> db['a'] = ...
    """
    return SQLiteShelf(*args, **kwargs)

# =============================================================================
## @class TmpSQLiteShelf
#  TEMPORARY SQLite-based ``shelve-like'' database with compressed content. 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2015-10-31
#  @see SQLiteShelf
class TmpSQLiteShelf(SQLiteShelf):
    """TEMPORARY SQLite-based ``shelve-like'' database with compressed content. 
    see SQLiteShelf
    """
    def __init__ ( self                                     ,
                   tablename      = 'ostap'                 ,
                   compress_level = zlib.Z_BEST_COMPRESSION , 
                   journal_mode   = "DELETE"                ,
                   protocol       = HIGHEST_PROTOCOL        ) :
        
        SQLiteShelf.__init__ ( self            ,
                               None            ,
                               'c'             ,
                               tablename       ,
                               True            , ## False , ## writeback/autocommit
                               compress_level  ,
                               journal_mode    ,
                               protocol        ) 
        
# =============================================================================
## open new TEMPORARY SQLiteShelve data base
# @code
# import SQLiteShleve as DBASE
# db = DBASE.tmpdb()
# db['a'] = ...
# @endcode
# @see TmpSQLiteShelf
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date   2015-10-31
def tmpdb(*args, **kwargs):
    """See documentation of the TmpSQLiteShelf class.
    >>> import SQLiteShleve as DBASE
    >>> db = DBASE.tmpdb()
    >>> db['a'] = ...
    """
    return TmpSQLiteShelf( *args , **kwargs )

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
