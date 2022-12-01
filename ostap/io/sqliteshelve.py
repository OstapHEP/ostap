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
# >>> db = DBASE.open ('a_db' , 'r' )  ## access existing dbase in read-only mode
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
# >>> db = DBASE.open ('a_db' )        ## access existing dbase in update mode
# ...
# >>> for key in db : print key
# ...
# >>> abcd = db['some_key']
#
# @endcode 
#
# @attention: When one tries to read the database with pickled ROOT object using newer
# version of ROOT, one could get a ROOT read error,
# in case of evoltuion in ROOT streamers for some  classes, e.g. <code>ROOT.TH1D</code>>
# @code 
# Error in <TBufferFile::ReadClassBuffer>: Could not find the StreamerInfo for version 2 of the class TH1D, object skipped at offset 19
# Error in <TBufferFile::CheckByteCount>: object of class TH1D read too few bytes: 2 instead of 878
# @endcode
# The solution is simple and described in  file ostap.io.dump_root
# @see ostap.io.dump_root
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

 Attention: When one tries to read the database with pickled ROOT object using newer
 version of ROOT, one could get a ROOT read error,
 in case of evoltuion in ROOT streamers for some  classes, e.g. ROOT.TH1D
 > Error in <TBufferFile::ReadClassBuffer>: Could not find the StreamerInfo for version 2 of the class TH1D, object skipped at offset 19
 > Error in <TBufferFile::CheckByteCount>: object of class TH1D read too few bytes: 2 instead of 878
 The solution is simple and described in  file ostap.io.dump_root
 - see ostap.io.dump_root
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
import os , sys , zlib, collections, datetime   
import sqlite3
from   ostap.io.sqlitedict  import SqliteDict
from   ostap.io.dbase       import Item, TmpDB 
from   ostap.core.meta_info import meta_info 
from   ostap.io.pklprotocol import PROTOCOL, HIGHEST_PROTOCOL, DEFAULT_PROTOCOL
# =============================================================================
try:
    from cPickle   import Pickler, Unpickler
except ImportError:
    from  pickle   import Pickler, Unpickler
# =============================================================================
try : 
    from io        import             BytesIO
except ImportError : 
    from cStringIO import StringIO as BytesIO         
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.sqliteshelve' )
else                      : logger = getLogger ( __name__ )
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
                   protocol       = PROTOCOL  ,
                   timeout        = 30        ) :
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
            
        if not 0 <= protocol <= HIGHEST_PROTOCOL :
            logger.warning ("Invalid protocol:%s" % protocol )
            protocol = PROTOCOL 

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
                              journal_mode = journal_mode ,
                              timeout      = timeout      )
        
        self.__compresslevel = compress_level 
        self.__protocol      = protocol
        self.__sizes         = {}

        if self.flag in (  'w' , 'n' ) :
            dct  = collections.OrderedDict() 
            dct  [ 'Created by'                  ] = meta_info.User
            dct  [ 'Created at'                  ] = datetime.datetime.now ().strftime( '%Y-%m-%d %H:%M:%S' )  
            dct  [ 'Created with Ostap version'  ] = meta_info.Ostap
            dct  [ 'Created with Python version' ] = meta_info.Python
            dct  [ 'Created with ROOT version'   ] = meta_info.ROOT 
            dct  [ 'Pickle protocol'             ] = protocol 
            dct  [ 'Compress level'              ] = self.__compresslevel 
            self [ '__metainfo__' ] = dct


    @property
    def  writeback  ( self ):
        """``writeback'' : the same as ``autocommit'':
        If one  enables `autocommit`, changes will be committed after each operation
        (more inefficient but safer).
        Otherwise, changes are committed on `self.commit()`,
        `self.clear()` and `self.close()`."""        
        return self.autocommit

    @property
    def compression ( self ) :
        """The  compression level from zlib"""
        return self.__compresslevel
    
    @property
    def compresslevel ( self ) :
        "``compress level'' : compression level"
        return self.__compresslevel 

    @property
    def protocol    ( self ) :
        """The pickling protocol"""
        return self.__protocol
    
    # =========================================================================
    ## iterator over good keys 
    def ikeys ( self , pattern = '' ) :
        """Iterator over avilable keys (patterns included).
        Pattern matching is performed accoriding to
        fnmatch/glob/shell rules [it is not regex!] 
        
        >>> db = ...
        >>> for k in db.ikeys('*MC*') : print(k)
        
        """
        keys_ = self.keys()
        
        if not pattern :
            good = lambda s,p : True
        else :
            import fnmatch
            good = lambda s,p : fnmatch.fnmatchcase ( k , p )
        
        for k in sorted ( keys_ ) :
            if good ( k , pattern ) : yield k

    ## list the avilable keys 
    def __dir ( self , pattern = '' ) :
        """ List the avilable keys (patterns included).
        Pattern matching is performed accoriding to
        fnmatch/glob/shell rules [it is not regex!] 
        
        >>> db = ...
        >>> db.ls() ## all keys
        >>> db.ls ('*MC*')
        
        """
        for k in self.ikeys ( patters ) : print ( key )
        
    # =========================================================================
    ## list the avilable keys 
    def ls    ( self , pattern = '' , load = True ) :
        """List the available keys (patterns included).
        Pattern matching is performed accoriding to
        fnmatch/glob/shell rules [it is not regex!] 

        >>> db = ...
        >>> db.ls() ## all keys
        >>> db.ls ('*MC*')        
        
        """
        n  = os.path.basename ( self.filename )
        ap = os.path.abspath  ( self.filename ) 
        
        try :
            fs = os.path.getsize ( self.filename )
        except :
            fs = -1
            
        if    fs < 0            : size = "???"
        elif  fs < 1024         : size = str(fs)
        elif  fs < 1024  * 1024 :
            size = '%.2fkB' % ( float ( fs ) / 1024 )
        elif  fs < 1024  * 1024 * 1024 :
            size = '%.2fMB' % ( float ( fs ) / ( 1024 * 1024 ) )
        else :
            size = '%.2fGB' % ( float ( fs ) / ( 1024 * 1024 * 1024 ) )
            
                        
        keys = [] 
        for k in self.ikeys ( pattern ): keys.append ( k )
        keys.sort()
        if keys : mlen = max ( [ len(k) for k in keys] ) + 2 
        else    : mlen = 2 
        fmt = ' --> %%-%ds : %%s' % mlen

        table = [ ( 'Key' , 'type' , '   size   ', ' created/modified ') ]

        meta = self.get( '__metainfo__' , {} )
        for k in meta :
            row = "META:%s" % k , '' , '' , str ( meta[k] )
            table.append ( row  ) 
        
        for k in keys :
            size = '' 
            ss   =   self.__sizes.get ( k , -1 )
            if    ss < 0    : size = '' 
            elif  ss < 1024 : size = '%7d   ' % ss 
            elif  ss < 1024 * 1024 :
                size = '%7.2f kB' %  ( float ( ss ) / 1024 )
            elif  ss < 1024 * 1024 * 1024 :
                size = '%7.2f MB' %  ( float ( ss ) / ( 1024 * 1024 ) )
            else :
                size = '%7.2f GB' %  ( float ( ss ) / ( 1024 * 1024 * 1024 ) )
                
            ot    = type ( self [ k ] )
            otype = ot.__cppname__ if hasattr ( ot , '__cppname__' ) else ot.__name__

            rawitem = self.__get_raw_item__ ( k )
            if isinstance ( rawitem , Item ) : timetag = rawitem.time
            else                             : timetag = '' 

            row = '{:15}'.format ( k ) , '{:15}'.format ( otype ) , size , timetag 
            table.append ( row )

        import ostap.logger.table as T
        t      = self.__class__.__name__
        title  = '%s:%s' % ( t  , n )
        maxlen = 0
        for row in table :
            rowlen = 0 
            for i in row : rowlen += len ( i )
            maxlen = max ( maxlen, rowlen ) 
        if maxlen + 3 <= len ( title ) :
            title = '<.>' + title [ -maxlen : ] 
        table = T.table ( table , title = title , prefix = '# ' )
        ll    = getLogger ( n )
        line  = 'Database %s:%s #keys: %d size: %s' % ( t , ap , len ( self ) , size )
        ll.info (  '%s\n%s' %  ( line , table ) )
 
   
    # =========================================================================
    ## clone the database into new one
    #  @code
    #  db  = ...
    #  ndb = db.clone ( 'new_file.db' )
    #  @endcode
    def clone ( self , new_name , keys = () ) : 
        """ Clone the database into new one
        >>> old_db = ...
        >>> new_db = new_db.clone ( 'new_file.db' )
        """
        new_db = SQLiteShelf ( new_name                           ,
                               mode           = 'n'               ,
                               tablename      = self.tablename    ,
                               writeback      = self.writeback    ,                                
                               protocol       = self.protocol     ,
                               compress_level = self.compression  , 
                               journal_mode   = self.journal_mode )
        
        ## copy the content
        if keys :
            for key in self.keys() :
                if key in keys     : new_db [ key ] = self [ key ]
        else : 
            for key in self.keys() : new_db [ key ] = self [ key ]
        
        new_db.sync ()  
        return new_db 

    # =============================================================================
    ## ``get-and-uncompress-item'' from dbase 
    def __get_raw_item__ ( self, key ):
        """ ``get-and-uncompress-item'' from dbase 
        """
        GET_ITEM = 'SELECT value FROM %s WHERE key = ?' % self.tablename
        item = self.conn.select_one ( GET_ITEM , ( key , ) )
        if item is None: raise KeyError(key)
        
        self.__sizes [ key ] = len ( item[0] ) 
        f     = BytesIO ( zlib.decompress ( item [ 0 ] ) ) 
        value = Unpickler ( f ) . load ()
        
        return value

    # =============================================================================
    ## ``get-and-uncompress-item'' from dbase 
    def __getitem__ ( self, key ):
        """ ``get-and-uncompress-item'' from dbase 
        """

        value  = self.__get_raw_item__ ( key )
        if isinstance ( value , Item ) : value = value.payload 
        return value

    # =============================================================================
    ## ``set-and-compress-item'' to dbase 
    def __setitem__ ( self , key , value ) :
        """ ``set-and-compress-item'' to dbase 
        """
        ADD_ITEM = 'REPLACE INTO %s (key, value) VALUES (?,?)' % self.tablename

        item = Item ( datetime.datetime.now ().strftime( '%Y-%m-%d %H:%M:%S' ) , value )

        f     = BytesIO ()
        p     = Pickler ( f , self.protocol )
        p.dump ( item )
        blob  = f.getvalue ( ) 
        zblob = zlib.compress ( blob , self.compression )
        
        self.__sizes [ key ] = len ( zblob )
        self.conn.execute ( ADD_ITEM, ( key , sqlite3.Binary ( zblob ) ) )


    # =========================================================================
    ## close and compress (if needed)
    def close ( self ) :
        """ Close the file (and compress it if required) 
        """

        if self.flag == 'c' : 
            dct = self.get ( '__metainfo__' , {} )
            dct  [ 'Updated at'                  ] = datetime.datetime.now().strftime( '%Y-%m-%d %H:%M:%S' )   
            dct  [ 'Updated by'                  ] = meta_info.User 
            dct  [ 'Updated with Ostap version'  ] = meta_info.Ostap 
            dct  [ 'Updated with Python version' ] = meta_info.Python 
            dct  [ 'Updated with ROOT version'   ] = meta_info.ROOT   
            self [ '__metainfo__' ] = dct

        return SqliteDict.close ( self )
        
    ## context manager
    def __enter__ ( self      ) :
        """Context manager"""
        return self

    ## context manager: close at exit
    def __exit__  ( self , *_ ) :
        """Context manager: close at exit
        """
        try :
            self.close()
        except :
            pass


# =============================================================================
## open new SQLiteShelve data base
#  @code
#  import SQLiteShleve as DBASE
#  db = DBASE.open( 'data.msql')
#  db['a'] = ...
#  @endcode
#  @see SQLiteShelf
def open ( *args , **kwargs ):
    """See documentation of the SQLiteShelf class.
    >>> import SQLiteShleve as DBASE
    >>> db = DBASE.open('data.msql','c')
    >>> db['a'] = ...
    """
    return SQLiteShelf ( *args , **kwargs )

# =============================================================================
## @class TmpSQLiteShelf
#  TEMPORARY SQLite-based ``shelve-like'' database with compressed content. 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2015-10-31
#  @see SQLiteShelf
class TmpSQLiteShelf(SQLiteShelf,TmpDB):
    """TEMPORARY SQLite-based ``shelve-like'' database with compressed content. 
    see SQLiteShelf
    """
    def __init__ ( self                                     ,
                   tablename      = 'ostap'                 ,
                   compress_level = zlib.Z_BEST_COMPRESSION , 
                   journal_mode   = "DELETE"                ,
                   protocol       = HIGHEST_PROTOCOL        ,
                   remove         = True                    ,
                   keep           = False                   ) :
        
        ## initialize the base: generate the name 
        TmpDB.__init__ ( self , suffix = '.lzdb' , remove = remove , keep = keep ) 
        
        ## open DB  
        SQLiteShelf.__init__ ( self            ,
                               self.tmp_name   ,
                               'c'             ,
                               tablename       ,
                               True            , ## False , ## writeback/autocommit
                               compress_level  ,
                               journal_mode    ,
                               protocol        ) 
            
    ## close and delete the file 
    def close ( self )  :
        ## close the shelve file
        SQLiteShelf.close ( self )
        ## delete the file
        TmpDB  .clean     ( self ) 
     
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
##                                                                      The END 
# =============================================================================
