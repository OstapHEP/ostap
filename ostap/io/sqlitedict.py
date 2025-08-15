#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# =============================================================================
# copied from https://github.com/RaRe-Technologies/sqlitedict
# copied from https://pypi.org/project/sqlitedict/
# =============================================================================
# 
# This code was inspired by:
#  * http://code.activestate.com/recipes/576638-draft-for-an-sqlite3-based-dbm/
#  * http://code.activestate.com/recipes/526618/
"""
A lightweight wrapper around Python's sqlite3 database, with a dict-like interface
and multi-thread access support::

>>> mydict = SqliteDict('some.db', autocommit=True) # the mapping will be persisted to file `some.db`
>>> mydict['some_key'] = any_picklable_object
>>> print mydict['some_key']
>>> print len(mydict) # etc... all dict functions work

Pickle is used internally to serialize the values. Keys are strings.

If you don't use autocommit (default is no autocommit for performance), then
don't forget to call `mydict.commit()` when done with a transaction.

"""
# =============================================================================
__all__ = (
    'SqliteDict' , ## sqlite3-persistent dictionary
    'issqlite3'  , ## is it a  sqlite3 file?
    )
# =============================================================================
from   collections import UserDict      as     DictClass
from   threading                        import Thread
from   queue import Queue
from   ostap.io.pickling                 import dumps, loads, HIGHEST_PROTOCOL as PICKLE_PROTOCOL
import sys, os, io, sqlite3, traceback, logging 
# =============================================================================# 
from ostap.logger.logger import getLogger
logger = getLogger( __name__ )
# =============================================================================
## Is it a sqlite3 file?
#  @code
#  ok = issqlite3 ( 'mydbase' ) 
#  @endcode
#  @see https://stackoverflow.com/questions/12932607/how-to-check-if-a-sqlite3-database-exists-in-python
def issqlite3 ( filename ) :
    """ Is it a sqlite3 file?
    >>> ok = issqlite3 ( 'mydbase' ) 
    - see https://stackoverflow.com/questions/12932607/how-to-check-if-a-sqlite3-database-exists-in-python
    """
    if not   os.path.exists  ( filename ) : return False
    if not   os.path.isfile  ( filename ) : return False
    if 100 > os.path.getsize ( filename ) : return False
    
    with io.open ( filename  , 'rb' ) as f :
        hdr = f.read(16)
        return hdr[:16] == b'SQLite format 3\x00'
    
    return False 

# =============================================================================
def reraise ( tp , value , tb = None):
    if value is None : value = tp ()
    if value.__traceback__ is not tb :
        raise value.with_traceback ( tb )
    raise value

# =============================================================================
## @class Connect
#  Helper class to implement "read-only" access to database 
class Connect ( object ) :
    """ Helper class to implement `read-only' access to database 
    """
    def __init__ ( self , filename , flag , **kwargs ) :
        
        self.__filename = filename
        self.__flag     = flag
        self.__kwargs   = kwargs
        self.__connect  = None
        self.__fd       = None

    # =========================================================================
    ## context manager ENTER
    def __enter__ ( self ) :

        if 'r' in self.flag :
            filename  = 'file:%s?mode=ro' % self.filename 
            self.__connect = sqlite3.connect ( filename , uri = True , **self.kwargs )
        else :
            self.__connect = sqlite3.connect ( self.filename         , **self.kwargs )
            
        return self.connect

    # ========================================================================
    ## context manager EXIT 
    def __exit__ ( self , *_ ) :

        if self.connect :
            self.__connect.close()
            self.__connect = None
        if not self.__fd is None  :
            os.close ( self.__fd )
            
    @property
    def filename ( self ) :
        """`filename' : the original filename"""
        return self.__filename

    @property
    def flag ( self ) :
        """`flag' : the access flag"""
        return self.__flag

    @property
    def kwargs ( self ) :
        """`kwargs` : addtional keyword arguments"""
        return self.__kwargs
    
    @property
    def connect ( self ) :
        """`connect' : the actual connection"""
        return self.__connect
    
# =============================================================================
## @class SqliteDict
#  A downgraded version of the original `SqliteDict`  (no pickling!)
#  - Keys are strings or bytes
#  - Values are bytes 
class SqliteDict(DictClass):
    """ A downgraded version of the original SqliteDict  (no pickling!)
      - Keys are strings or bytes
      - Values are bytes 
    """
    VALID_FLAGS =  ( 'c', 'r', 'w', 'n' ) 

    def __init__ ( self                    ,
                   filename     = None     ,
                   flag         = 'c'      ,
                   tablename    = 'Dict'   ,
                   autocommit   = False    ,
                   ## journal_mode = "DELETE" ,
                   journal_mode = "OFF"    ,
                   timeout      = 30       ) :
        """ Initialize a thread-safe sqlite-backed dictionary. The dictionary will
        be a table `tablename` in database file `filename`. A single file (=database)
        may contain multiple tables.

        If no `filename` is given, a random file in temp will be used (and deleted
        from temp once the dict is closed/deleted).

        If you enable `autocommit`, changes will be committed after each operation
        (more inefficient but safer). Otherwise, changes are committed on `self.commit()`,
        `self.clear()` and `self.close()`.

        Set `journal_mode` to 'OFF' if you're experiencing sqlite I/O problems
        or if you need performance and don't care about crash-consistency.

        The `flag` parameter. Exactly one of:
          'c': default mode, open for read/write, creating the db/table if necessary.
          'w': open for r/w, but drop `tablename` contents first (start with empty table)
          'r': open as read-only
          'n': create a new database (erasing any existing tables, not just `tablename`!).

        """
        tablename = 'Dict'
        
        ## temporary DB ?
        self.in_temp = filename is None
        if self.in_temp:
            import ostap.utils.cleanup as CU 
            filename = CU.CleanUp.tempfile ( prefix = 'ostap-tmpdb-' , suffix = '.sqldb' )

        ## valid flag ? 
        self.flag = flag
        if self.flag not in SqliteDict.VALID_FLAGS:
            raise RuntimeError ( "Unrecognized flag: %s" % flag )

        if self.flag == 'n':
            if os.path.exists ( filename ) and os.path.isfile ( fielname ) :
                try : 
                    os.remove ( filename )
                except :
                    pass
                
        dirname = os.path.dirname(filename)
        if dirname:
            if not os.path.exists ( dirname ):
                raise RuntimeError ( 'Error! The directory does not exist, %s' % dirname )
            
        self.filename = filename
        if '"' in tablename:
            raise ValueError ( 'Invalid tablename %r' % tablename)
        self.tablename    = tablename
        self.autocommit   = autocommit
        self.journal_mode = journal_mode
        self.timeout      = timeout 

        ## check it! 
        with Connect ( self.filename , self.flag , timeout = self.timeout ) :
            pass  

        logger.debug ( "opening Sqlite table %r in %s" % ( tablename , filename ) )
        
        MAKE_TABLE = 'CREATE TABLE IF NOT EXISTS "%s" (key TEXT PRIMARY KEY, value BLOB)' % self.tablename
        self.conn = self._new_conn()
        self.conn.execute(MAKE_TABLE)
        self.conn.commit()
        if flag == 'w' :
            self.clear()

        ## read-only DB 
        self.readonly =  'r' == flag

    #  ========================================================================
    ## make new connection 
    def _new_conn ( self ) :
        return SqliteMultithread ( self.filename                    ,
                                   self.flag                        ,
                                   autocommit   = self.autocommit   ,
                                   journal_mode = self.journal_mode ,
                                   timeout      = self.timeout      )

    # =========================================================================
    ## Context manager: ENTER 
    def __enter__(self):
        """ Context manager:ENTER
        """
        if self.conn is None: self.conn = self._new_conn()
        return self

    # =========================================================================
    ## Context manager: EXIT 
    def __exit__(self, *exc_info):
        """ Context manager: EXIT 
        """
        self.sync  ()  ## ADDED 
        self.close ()
        
    def __str__(self):
        return "SqliteDict(%s)" % (self.filename)
    
    def __repr__(self):
        return str(self)  # no need of something complex

    # =========================================================================
    ## number of keys/rows 
    def __len__ ( self ) :
        # `select count (*)` is super slow in sqlite (does a linear scan!!)
        # As a result, len() is very slow too once the table size grows beyond trivial.
        # We could keep the total count of rows ourselves, by means of triggers,
        # but that seems too complicated and would slow down normal operation
        # (insert/delete etc).
        GET_LEN = 'SELECT COUNT(*) FROM "%s"' % self.tablename
        rows = self.conn.select_one(GET_LEN)[0]
        return rows if rows is not None else 0

    # =========================================================================
    def __bool__(self):
        # No elements is False, otherwise True
        GET_MAX = 'SELECT MAX(ROWID) FROM "%s"' % self.tablename
        m = self.conn.select_one(GET_MAX)[0]
        # Explicit better than implicit and bla bla
        return True if m is not None else False

    # =========================================================================
    __nonzero__ = __bool__ 
    # =========================================================================
    def tables ( self ) :
        """ Get list of tables in DBASE"""
        GET_TABLES = "SELECT name FROM sqlite_master WHERE type='table';"
        tables = [] 
        for table in self.conn.select ( GET_TABLES ) :
            tables.append ( table[0]  ) 
        return tuple ( tables )

    # =========================================================================
    ## encode value: input -> dbase 
    def encode_value ( self , value ) : 
        """ encode value: 
        - input -> dbase 
        """
        return sqlite3.Binary ( value ) 

    # =========================================================================
    ## decode value: dbase -> ouput (no-op) 
    def decode_value ( self , value ) :
        """ decode value: (no-op) 
        - dbase -> output 
        """
        return value 

    # =========================================================================
    ## iterator over keys 
    def iterkeys   ( self ):
        """ Iterator over keys
        >>> db = ...
        >>> for k in db.iterkeys() : ... 
        """
        GET_KEYS = 'SELECT key FROM "%s" ORDER BY rowid' % self.tablename
        for key in self.conn.select ( GET_KEYS ) :
            yield key [ 0 ]

    # =========================================================================
    ## iterator over values 
    def itervalues ( self ) :
        """ Iterator over values 
        >>> db = ...
        >>> for v in db.itervalues () : ... 
        """
        GET_VALUES = 'SELECT value FROM "%s" ORDER BY rowid' % self.tablename
        for value in self.conn.select(GET_VALUES):
            yield self.decode_value ( value [ 0 ] )

    # =========================================================================
    ## iterator over (key,value) pairs 
    def iteritems  ( self ) :
        """ Iterator over (key,value) pairs 
        >>> db = ...
        >>> for k,v in db.iteritems () : ... 
        """
        GET_ITEMS = 'SELECT key, value FROM "%s" ORDER BY rowid' % self.tablename
        for key, value in self.conn.select ( GET_ITEMS ) :
            yield  key , self.decode_value ( value )

    # =========================================================================
    ## iterator over keys 
    def keys ( self ) :
        """ Iterator over keys
        >>> db = ...
        >>> for k in db.keys() : ... 
        """
        return self.iterkeys() 

    # =========================================================================
    ## iterator over values 
    def values ( self ) :
        """ Iterator over values 
        >>> db = ...
        >>> for v in db.values () : ... 
        """        
        return self.itervalues() 

    # =========================================================================
    ## iterator over *key,value) pairs 
    def items(self):
        """ Iterator over (key,value) pairs 
        >>> db = ...
        >>> for k,v in db.items () : ... 
        """
        return self.iteritems() 
    
    # =========================================================================
    ## The key in dbase? 
    def __contains__(self, key):
        """ Is the key in the dbase? """
        HAS_ITEM = 'SELECT 1 FROM "%s" WHERE key = ?' % self.tablename
        return self.conn.select_one ( HAS_ITEM , ( key , ) ) is not None
    
    # =========================================================================
    ## get value form DB 
    def get ( self , key , default = None ) :
        """ Get value from dbase 
        >>> db = ...
        >>> value = db.get ( key , 42 ) 
        """
        GET_ITEM = 'SELECT value FROM "%s" WHERE key = ?' % self.tablename
        item = self.conn.select_one(GET_ITEM, (key,))
        if item is None : return default 
        return self.decode_value  ( item [ 0 ] )

    # =========================================================================
    ## get value from DB 
    def __getitem__ ( self , key ) :
        """ Get value from dbase 
        >>> db = ...
        >>> value = db [ key ] 
        """
        GET_ITEM = 'SELECT value FROM "%s" WHERE key = ?' % self.tablename
        item = self.conn.select_one(GET_ITEM, (key,))
        if item is None: raise KeyError(key)        
        return self.decode_value  ( item [ 0 ] )

    # =========================================================================
    ## write value to dbase
    def __setitem__(self, key, value):
        """ Write value to dbase 
        >>> db = ...
        >>> db [ key ] = value 
        """
        if self.readonly :
            raise RuntimeError ( 'Refusing to write to read-only SqliteDict' )
        ADD_ITEM = 'REPLACE INTO "%s" (key, value) VALUES (?,?)' % self.tablename        
        self.conn.execute(ADD_ITEM, ( key , self.encode_value ( value ) ) )
        
    # =========================================================================
    ## delete the element from dbase 
    def __delitem__(self, key):
        """ Delete the element fomr dbase 
        >>> db = ...
        >>> del db [ key ] 
        """
        if self.readonly :
            raise RuntimeError ( 'Refusing to delete from read-only SqliteDict' )
        if key not in self: raise KeyError ( key )        
        DEL_ITEM = 'DELETE FROM "%s" WHERE key = ?' % self.tablename
        self.conn.execute ( DEL_ITEM , ( key , ) )

    # =========================================================================
    ## update many elements in a single go 
    def update ( self , items = () , **kwargs  ) :
        """ Update many elements in a single go
        >>> db = ...
        >>> db.update ( ... ) 
        """
        if self.readonly :
            raise RuntimeError('Refusing to update read-only SqliteDict')
        
        if   hasattr ( items , 'items'     ) : items = items.items     ()
        elif hasattr ( items , 'iteritems' ) : items = items.iteritems ()
        ##
        items = [ item for item in items ]
        ## 
        if  kwargs : items += [  item for item in kwargs.items () ]
        ## 
        items = [ ( k , self.encode_value ( v ) ) for k, v in items ]
        ## 
        UPDATE_ITEMS = 'REPLACE INTO "%s" (key, value) VALUES (?, ?)' % self.tablename
        self.conn.executemany ( UPDATE_ITEMS, items)
        
    # =========================================================================
    ## iterate over keys 
    def __iter__ ( self ) :
        """ Iterat eover keys 
        >>> db = ...
        >>> for key in db : ... 
        """
        return self.iterkeys()

    # =========================================================================
    ## clear database 
    def clear ( self ):
        """ Clear database 
        >>> db = ....
        >>> db.clear() 
        """
        if self.readonly :
            raise RuntimeError('Refusing to clear read-only SqliteDict')

        CLEAR_ALL = 'DELETE FROM "%s";' % self.tablename  # avoid VACUUM, as it gives "OperationalError: database schema has changed"
        self.conn.commit ()
        self.conn.execute(CLEAR_ALL)
        self.conn.commit ()
        
    @staticmethod
    def get_tablenames(filename):
        """ Get the names of the tables in an sqlite db as a list"""
        if not os.path.isfile(filename):
            raise IOError('file %s does not exist' % (filename))
        GET_TABLENAMES = 'SELECT name FROM sqlite_master WHERE type="table"'
        
        with Connect ( filename , 'r' , timeout = self.timeout ) as conn:
            cursor = conn.execute ( GET_TABLENAMES )
            res    = cursor.fetchall()

        return [ name [ 0 ] for name in res ]

    # =========================================================================
    ## persissts all data on dist 
    def commit ( self , blocking = True ):
        """ Persist all data to disk.
        
        When `blocking` is False, the commit command is queued, but the data is
        not guaranteed persisted (default implication when autocommit=True).
        """
        if self.conn is not None: self.conn.commit ( blocking )

    # =========================================================================
    ## persissts all data on dist 
    sync = commit

    # =========================================================================
    ## close database 
    def close ( self           ,
                do_log = True  ,
                force  = False ) :
        
        if do_log: logger.debug ( "Closing %s" % self)
        if not self.conn is None:
            if self.conn.autocommit and not force:
                # typically calls to commit are non-blocking when autocommit is
                # used.  However, we need to block on close() to ensure any
                # awaiting exceptions are handled and that all data is
                # persisted to disk before returning.
                self.conn.commit ( blocking = True )
            self.conn.close ( force = force )
            self.conn = None

        # =====================================================================
        ## remove temporary DBASE
        # =====================================================================
        if self.in_temp: # ====================================================
            # =================================================================
            try: # ============================================================
                # =============================================================
                os.remove ( self.filename )
                # =============================================================
            except: # =========================================================
                # =============================================================
                pass

    # =========================================================================
    ## terminate session and delete dbase 
    def terminate ( self ) :
        """ Close session and Delete the underlying database file
        - Use with care.
        """
        if self.readonly :
            raise RuntimeError('Refusing to terminate read-only SqliteDict')
        
        self.close()

        if self.filename == ':memory:': return

        logger.debug ( "deleting %s" % self.filename )
        # =====================================================================
        try: # ================================================================
            # =================================================================
            if os.path.exists ( self.filename ) and os.path.isfile ( self.filename ) :
                os.remove ( self.filename )
        except ( OSError , IOError ) : # ======================================
            # =================================================================
            logger.exception("failed to delete %s" % (self.filename))

    # =========================================================================
    ## close dbase 
    def __del__ ( self ) :
        try:
            self.close ( do_log = False ,  force = True )
        except Exception:
            # prevent error log flood in case of multiple SqliteDicts
            # closed after connection lost (exceptions are always ignored
            # in __del__ method.
            pass
        

# =============================================================================
# @class SqliteMultithread
class SqliteMultithread(Thread):
    """ Wrap sqlite connection in a way that allows concurrent requests from multiple threads.

    This is done by internally queueing the requests and processing them sequentially
    in a separate thread (in the same order they arrived).

    """
    def __init__ ( self         ,
                   filename     ,
                   flag         ,
                   autocommit   ,
                   journal_mode ,
                   timeout = 5  ) :
        
        super ( SqliteMultithread , self) .__init__()
        self.filename     = filename
        self.flag         = flag 
        self.autocommit   = autocommit
        self.journal_mode = journal_mode
        self.timeout      = timeout
        
        # use request queue of unlimited size
        self.reqs     = Queue()
        
        ## self.setDaemon(True)  # python2.5-compatible
        self.daemon    = True
        
        self.exception = None
        self.log       = logging.getLogger('sqlitedict.SqliteMultithread')
        self.start()


    # =========================================================================
    def run ( self ) :
        
        if self.autocommit :
            connect = Connect ( self.filename , self.flag , timeout = self.timeout , isolation_level = None, check_same_thread = False )
        else:
            connect = Connect ( self.filename , self.flag , timeout = self.timeout , check_same_thread = False)
            
        with connect as conn :
            return self.run__ (  conn )

    # =========================================================================
    def run__ ( self , conn ):
        
        if self.autocommit:
            conn = sqlite3.connect ( self.filename , isolation_level = None , check_same_thread = False )
        else:
            conn = sqlite3.connect ( self.filename , check_same_thread = False )
        
        conn.execute('PRAGMA journal_mode = %s' % self.journal_mode)
        conn.text_factory = str
        cursor = conn.cursor()
        conn.commit()
        cursor.execute('PRAGMA synchronous=OFF')

        res = None
        while True:
            req, arg, res, outer_stack = self.reqs.get()
            if req == '--close--':
                assert res, ('--close-- without return queue', res)
                break
            elif req == '--commit--':
                conn.commit()
                if res:
                    res.put('--no more--')
            else:
                try:
                    cursor.execute(req, arg)
                except Exception as err:
                    self.exception = (e_type, e_value, e_tb) = sys.exc_info()
                    inner_stack = traceback.extract_stack()

                    # An exception occurred in our thread, but we may not
                    # immediately able to throw it in our calling thread, if it has
                    # no return `res` queue: log as level ERROR both the inner and
                    # outer exception immediately.
                    #
                    # Any iteration of res.get() or any next call will detect the
                    # inner exception and re-raise it in the calling Thread; though
                    # it may be confusing to see an exception for an unrelated
                    # statement, an ERROR log statement from the 'sqlitedict.*'
                    # namespace contains the original outer stack location.
                    self.log.error('Inner exception:')
                    for item in traceback.format_list(inner_stack):
                        self.log.error(item)
                    self.log.error('')  # deliniate traceback & exception w/blank line
                    for item in traceback.format_exception_only(e_type, e_value):
                        self.log.error(item)

                    self.log.error('')  # exception & outer stack w/blank line
                    self.log.error('Outer stack:')
                    for item in traceback.format_list(outer_stack):
                        self.log.error(item)
                    self.log.error('Exception will be re-raised at next call.')

                if res:
                    for rec in cursor:
                        res.put(rec)
                    res.put('--no more--')

                if self.autocommit:
                    conn.commit()

        self.log.debug('received: %s, send: --no more--', req)
        conn.close()
        res.put('--no more--')

    # =========================================================================
    def check_raise_error(self):
        """ Check for and raise exception for any previous sqlite query.

        For the `execute*` family of method calls, such calls are non-blocking and any
        exception raised in the thread cannot be handled by the calling Thread (usually
        MainThread).  This method is called on `close`, and prior to any subsequent
        calls to the `execute*` methods to check for and raise an exception in a
        previous call to the MainThread.
        """
        if self.exception:
            e_type, e_value, e_tb = self.exception

            # clear self.exception, if the caller decides to handle such
            # exception, we should not repeatedly re-raise it.
            self.exception = None

            self.log.error('An exception occurred from a previous statement, view '
                           'the logging namespace "sqlitedict" for outer stack.')

            # The third argument to raise is the traceback object, and it is
            # substituted instead of the current location as the place where
            # the exception occurred, this is so that when using debuggers such
            # as `pdb', or simply evaluating the naturally raised traceback, we
            # retain the original (inner) location of where the exception
            # occurred.
            reraise(e_type, e_value, e_tb)

    # =========================================================================
    def execute(self, req, arg=None, res=None):
        """ `execute` calls are non-blocking: just queue up the request and return immediately.
        """
        self.check_raise_error()

        # NOTE: This might be a lot of information to pump into an input
        # queue, affecting performance.  I've also seen earlier versions of
        # jython take a severe performance impact for throwing exceptions
        # so often.
        stack = traceback.extract_stack()[:-1]
        self.reqs.put((req, arg or tuple(), res, stack))

    def executemany(self, req, items):
        for item in items:
            self.execute(req, item)
        self.check_raise_error()

    # =========================================================================
    def select ( self , req , arg = None ) :
        """ Unlike sqlite's native select, this select doesn't handle iteration efficiently.

        The result of `select` starts filling up with values as soon as the
        request is dequeued, and although you can iterate over the result normally
        (`for res in self.select(): ...`), the entire result will be in memory.
        """
        res = Queue()  # results of the select will appear as items in this queue
        self.execute(req, arg, res)
        while True:
            rec = res.get()
            self.check_raise_error()
            if rec == '--no more--':
                break
            yield rec

    def select_one(self, req, arg=None):
        """ Return only the first row of the SELECT, or None if there are no matching rows."""
        try:
            return next(iter(self.select(req, arg)))
        except StopIteration:
            return None

    def commit(self, blocking=True):
        if blocking:
            # by default, we await completion of commit() unless
            # blocking=False.  This ensures any available exceptions for any
            # previous statement are thrown before returning, and that the
            # data has actually persisted to disk!
            self.select_one('--commit--')
        else:
            # otherwise, we fire and forget as usual.
            self.execute('--commit--')

    def close(self, force=False):
        if force:
            # If a SqliteDict is being killed or garbage-collected, then select_one()
            # could hang forever because run() might already have exited and therefore
            # can't process the request. Instead, push the close command to the requests
            # queue directly. If run() is still alive, it will exit gracefully. If not,
            # then there's nothing we can do anyway.
            self.reqs.put(('--close--', None, Queue(), None))
        else:
            # we abuse 'select' to "iter" over a "--close--" statement so that we
            # can confirm the completion of close before joining the thread and
            # returning (by semaphore '--no more--'
            self.select_one('--close--')
            self.join()

                        
# =============================================================================
#                                                                       The END 
# =============================================================================
