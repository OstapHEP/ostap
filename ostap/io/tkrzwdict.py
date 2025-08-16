#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Simplistic lightweighted  class that provides  wrapping to tkrzw dbase
# @see https://dbmx.net/tkrzw/
"""
Simplistic lightweighted  class that provides  wrapping to tkrzw dbase
- see https://dbmx.net/tkrzw/
"""
# =============================================================================
__all__ = (
    'TrkzwDict'  , ## tkrzw-persistent dictionary
    'istkrzw'   , ## is it a tkrzw file?
)
# =============================================================================
from   collections.abc import MutableMapping
import os, io, shutil 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.io.tkrzwdict' )
else                       : logger = getLogger( __name__     )
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import tkrzw
    # =========================================================================
except ImportError: # =========================================================
    # =========================================================================
    logger.warning ("No `tkrzw` module is available!")
    tkrzw= None

# =============================================================================
## Is it a tkzw file?
#  @code
#  ok = istkrzw ( 'mydbase' ) 
#  @endcode
def istkrzw ( filename ) :
    """ Is it a tkrzw file?
    >>> ok = istkrzw  ( 'mydbase' ) 
    """
    if not   os.path.exists  ( filename ) : return False
    if not   os.path.isfile  ( filename ) : return False
    if 100 > os.path.getsize ( filename ) : return False
    
    with io.open ( filename  , 'rb' ) as f :
        hdr = f.read(16)
        if   s16 [:9] == b'TkrzwHDB\n': return True  ## Hash 
        elif s16 [:9] == b'TkrzwSDB\n': return True  ## Skip
        elif s16 [:4] == b'TDB\n'     : return True  ## Tree 

    return False

# =============================================================================
## @class TkrzwDict
#  Persistent dictionary based on TKRZW
#  @attention both keys and values are bytestrings!
#  @code
#
#  mydict = TkrzwDict ( 'some.db' , 'c' )
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
class TkrzwDict ( MutableMapping ) :
    """ Persistent dictionary based on TKRZW 
    - both keys and values are bytestrings! 

    >>> mydict = TkrzwDict ( 'some.db' , 'c' )
    >>> mydict [ key ] = value
    
    >>> for k   in mydict          : print ( k )
    >>> for k   in mydict.keys()   : print ( k )
    >>> for v   in mydict.values() : print ( v )
    >>> for k,v in mydict.items () : print ( k , v )

    >>> val  = mdict[ key]

    >>> del mydict [ key ]
    
    """
    FLAGS      = 'rwcn'
    EXTENSIONS = ( 'tkh'  , 'tkt'     , 'tks' )
    CLASSES    = ( 'hash' , 'hashdbm' ,
                   'tree' , 'treedbm' ,
                   'skip' , 'skipdbm' )
    
    def __init__ ( self         ,
                   path         ,
                   flag   = 'r' , **kwargs )  :
        
        ## check presence of the tkrzwz module 
        assert trkzw , "`tkrzw` module is not available!"
        
        assert 1 == len ( flag ) and isinstance ( flag , str ) and flag.lower() in sefl.FLAGS , "Invalid `flag`:%s" % flag
        
        flag = flag.lower() 
        

        conf = { 'truncate'   : False ,
                 'concurrent' : True  ,
                 'no_create'  : True  , 
                 'no_wait'    : True  , 
                 'no_lock'    : True  }
        
        conf.update ( kwargs )

        dbtype = conf.pop ( 'dbm' ) 

        
        writable  = 'r' != flag  

        if   'r' == flag :
            conf [ 'no_create' ] = True 
        elif 'w' == flag :
            conf [ 'no_create' ] = False 
        elif 'c' == flag :
            conf [ 'no_create' ] = False 
        elif 'n' == flag :
            conf [ 'no_create' ] = False 
            conf [ 'truncate'  ] = True            
        else :
            raise ValueError ( "Invalid flag: %s" % flag )
        
        if 'n' == flag and path and os.path.exists ( path ) :
            if   os.path.isfile ( path ) :
                try    : os.unlink ( path )
                except : pass 
            elif os.path.isdir  ( path ) :
                try    : shutil.rmtree ( path )
                except : pass                    
            if os.path.exists ( path ) :
                logger.warning ( 'path exists!' )

                
        self.__conf = conf
        self.__dbm  = tkrzw.DBM()
        
        sc = self.__dbm.Open ( path , writeable , **conf )
        
        self.__path     = path
        self.__conf     = conf

        
    @property
    def db ( self ) :
        """`db` : the actual Tkrzw database ( same ad `dbm`) """
        return self.__db
    
    @property
    def dbm ( self ) :
        """`dbm` : the actual Tkrzw database ( same ad `db`) """
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
        vv = self.encode_value ( value )
        self.db.Set ( key , vv  )  

    # =========================================================================
    ## iterate over key,value pairs 
    def items ( self ) :
        """ Iterate over key, value pairs
        >>> db = ...
        >>> for key, value in db.items () : ...
        """
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
        for k, _ in self.db : yield key

    # ========================================================================
    ## Iterator over values
    def values ( self ) :
        """ Iterator over values
        >>> db = ...
        >>> for value in db.values() : ...
        """
        for _ , value in self.db : yield self.decode_value ( value )

    # =========================================================================
    ## check the key in database 
    def __contains__ ( self , key ) :
        """ Check the key in database """
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
        return len ( self.db ) 

    # =======================================================================--
    ## pop the certain key from the database 
    def pop ( self, key , default = None ) :
        """ pop the certain key from database
        >>> db = ...
        >>> value = db.pop ( key ) 
        """
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
        sc = self.db.Remove ( key )
        
    # ========================================================================
    ## sync database 
    def sync ( self ) :
        """ Sync database
        >>> db = ...
        >>> db.sync () 
        """
        self.db.Synchronize ( True )

    # =========================================================================
    ## close the database 
    def close ( self ) :
        """ Close database 
        >>> db = ...
        >>> db.close() 
        """
        if not self.db : return
        
        if self.db.IsWritable () and self.db.ShouldBeRebuild() :
            self.db.Rebuild()
            
        self.db.Close ()
        del self.__db
        
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
if not tkrzw  : __all__ = () 

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if not tkrzw : logger.warning ( "No `tkrzw` module is found!")


# =============================================================================
#                                                                       The END 
# =============================================================================
