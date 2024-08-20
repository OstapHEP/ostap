#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file lmdbdict.py
# 
# A lightweight wrapper around `lmdb` database, with a dict-like interface
#
# @code
#
# mydict = LmdbDict ( 'some.db' , 'c' )
#
# mydict [ key ] = value
#
# for k   in mydict          : print ( k )
# for k   in mydict.keys()   : print ( k )
# for v   in mydict.values() : print ( v )
# for k,v in mydict.items () : print ( k , v )
#
# val  = mdict[ key]
#
# del mydict [ key ]
#
# @endcode 
#
# @attention both keys and valeus are bytestrings!
# 
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2024-08-19
# =============================================================================
""" A lightweight wrapper around `lmdb` database, with a dict-like interface

    >>> mydict = LmdbDict ( 'some.db' , 'c' )
    >>> mydict [ key ] = value

    >>> for k   in mydict          : print ( k )
    >>> for k   in mydict.keys()   : print ( k )
    >>> for v   in mydict.values() : print ( v )
    >>> for k,v in mydict.items () : print ( k , v )

    >>> val  = mdict[ key]
    >>> del mydict [ key ] 

- both keys and values are bytestrings!

"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2024-08-19"
__version__ = "$Revision$"
__all__     = (
    'lmdb'     ,
    'LmdbDict' 
)
# =============================================================================
import sys, os
# =============================================================================
if ( 3 , 5 ) <= sys.version_info : from collections.abc import MutableMapping
else                             : from collections     import MutableMapping
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.io.lmdbdict' )
else                       : logger = getLogger( __name__     )
# =============================================================================
lmdb = None 
if ( 3 , 7 ) <= sys.version_info :
    try :
        import lmdb
    except ImportError:
        logger.warning ("No `lmdb` module is available!")
        pass
    
# =============================================================================
## Is it a lmdb subdir?
#  @code
#  ok = islmdb( 'my_dbase' ) 
#  @endcode
def islmdb ( path ) :
    """  Is it a lmdb subdir?
    >>> ok = islmdb( 'mydbase' ) 
    """
    if not   os.path.exists  ( path ) : return False
    if not   os.path.isdir   ( path ) : return False
    data = os.path.join ( path , 'data.mdb')
    if not   os.path.exists  ( data ) : return False
    if not   os.path.isfile  ( data ) : return False
    
    try :
        if lmdb :
            db = lmdb.open ( path , readonly = True , create = False )
            return True 
    except:
        pass

    return True 

# =============================================================================
## @class LmdbDict
#  Persistent dictionary based on LMDB
#  @attention both keys and valeus are bytestrings!
#  @code
#
#  mydict = LmdbDict ( 'some.db' , 'c' )
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
class LmdbDict ( MutableMapping ) :
    """ Persistent dictionary based on LMDB 
    - both keys and values are bytestrings! 

    >>> mydict = LmdbDict ( 'some.db' , 'c' )
    >>> mydict [ key ] = value
    
    >>> for k   in mydict          : print ( k )
    >>> for k   in mydict.keys()   : print ( k )
    >>> for v   in mydict.values() : print ( v )
    >>> for k,v in mydict.items () : print ( k , v )

    >>> val  = mdict[ key]

    >>> del mydict [ key ]
    
    """
    def __init__ ( self              ,
                   path              ,
                   flag      = 'r'   ,
                   autogrow  = 14    , **kwargs )  :

        ## check presence of the lmdb module 
        assert lmdb , "`lmdb` module is not available!"

        conf = { 'map_size'  : 2**24 ,
                 'map_async' : True  ,
                 'max_dbs'   : 1     , 
                 'mode'      : 0o755 }
        
        conf.update ( kwargs )
        
        if   'r' == flag :
            ## Open existing database for reading only (default)
            env = lmdb.open ( path , readonly = True  , create = False , **conf )
        elif 'w' == flag :
            ## Open existing database for reading and writing
            env = lmdb.open ( path , readonly = False , create = False , **conf )
        elif 'c' == flag :
            # Open database for reading and writing, creating it if it doesn't exist
            env = lmdb.open ( path , readonly = False , create = True  , **conf )
        elif 'n' == flag :
            ## remove_lmdbm(file)
            env = lmdb.open ( path , readonly = False , create = True , **conf  )
        else :
            raise ValueError ( "Invalid flag: %s" % flag )

        assert isinstance ( autogrow , int ) and  0 <= autogrow , "Invalid `autogrow` parameter: %s" % autogrow
        
        self.__autogrow = min ( autogrow , 32 )        
        self.__db   = env 
        
    @property
    def db ( self ) :
        """`db` : the actual LMDB database/Environment"""
        return self.__db 

    @property
    def autogrow ( self ) :
        """`autogrow' : allow autimatic resize of the database if/when needed of  1<times<32"""
        return self.__autogrow
    
    @property
    def map_size ( self ) :
        return self.db.info()["map_size"]
    @map_size.setter
    def map_size ( self, value ) :
        self.db.set_mapsize ( value )

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
        with self.db.begin ( write = False ) as txn:
            value = txn.get ( key )
        return default if value is None else self.decode_value ( value ) 

    # =========================================================================
    ## Get data from DB
    def __getitem__( self , key ) :
        """ Get data from DB
        >>> db = ...
        >>> value = db [ key ]
        """
        with self.db.begin ( write = False ) as txn:
            value = txn.get ( key )            
        if value is None: raise KeyError( 'Invalid key:%s' % key )        
        return self.decode_value ( value ) 

    # =========================================================================
    ## Add data to DB
    def __setitem__ ( self , key , value ) :
        """ Add data to DB 
        >>> db = ...
        >>> db [ key ] = value 
        """
        nn = max ( 1 , self.autogrow )
        vv = self.encode_value ( value ) 
        for i in range ( nn ):
            try:
                with self.db.begin ( write = True ) as txn:
                    txn.put ( key , vv )
                    return
            except lmdb.MapFullError :
                if not self.autogrow : raise
                new_map_size  = self.map_size * 2
                self.map_size = new_map_size
                logger.info ( 'DBASE is resized to %s' % self.map_size )
                
        raise lmdb.MapFullError ("LmdbDict: Failure to autogrow!") 

    # =========================================================================
    ## iterate over key,value pairs 
    def items ( self ) :
        """ Iterate over key, value pairs
        >>> db = ...
        >>> for key, value in db.items () : ...
        """
        with self.db.begin ( write = False ) as txn :
            for key , value in txn.cursor().iternext ( keys   = True ,
                                                       values = True ) :                
                yield  key , self.decode_value ( value ) 
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
        with self.db.begin ( write = False ) as txn:
            for key in txn.cursor().iternext ( keys   = True  ,
                                               values = False ):
                yield key 

    # ========================================================================
    ## Iterator over values
    def values ( self ) :
        """ Iterator over values
        >>> db = ...
        >>> for value in db.values() : ...
        """
        with self.db.begin( write = False ) as txn:
            for value in txn.cursor().iternext ( keys   = False ,
                                                 values = True  ) :
                yield self.decode_value ( value ) 

    # =========================================================================
    ## check the key in database 
    def __contains__ ( self , key ) :
        """ Check the key in database """
        with self.db.begin ( write = False ) as txn:
            value = txn.get ( key )
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
        with self.db.begin ( write = False ) as txn:
            return txn.stat()["entries"]

    # =======================================================================--
    ## pop the certain key from the database 
    def pop ( self, key , default = None ) :
        """ pop the certain key from database
        >>> db = ...
        >>> value = db.pop ( key ) 
        """
        with self.db.begin ( write = True ) as txn:
            value = txn.pop ( key )
        if value is None: return default
        return value

    # =========================================================================
    ## delete certain key from the database 
    def __delitem__( self , key ) :
        """ Delete certain key from the database
        >>> db = ...
        >>> value = db.pop ( key ) 
        """
        with self.db.begin ( write = True ) as txn:
            txn.delete ( key )

    # ========================================================================
    ## sync database 
    def sync ( self ) :
        """ Sync database
        >>> db = ...
        >>> db.sync () 
        """
        self.db.sync()

    # =========================================================================
    ## close the database 
    def close ( self ) :
        """ Close database 
        >>> db = ...
        >>> db.close() 
        """
        self.db.close ()

    # =========================================================================
    ## context manager: ENTER 
    def __enter__ ( self ) :
        """ Context manager: ENTER"""
        return self
    
    # =========================================================================
    ## context manager: EXIT  
    def __exit__ ( self , *args ) :
        """ Context manager: EXIT"""        
        self.close()

# ==============================================================================
## disable imports 
if not lmdb : __all__ = () 

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if not lmdb : logger.warning ( "No `lmdb` module is found!")
        
# =============================================================================
#                                                                       The END 
# =============================================================================

