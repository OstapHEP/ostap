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
    'dbopen'     , ## open database
    'Item'       , ## item: named tuple (time,payload)
    'use_bsddb3' , ## make use of bsbdb3 ?  
    )
# =============================================================================
import sys, os, collections
from   ostap.logger.logger  import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.dbase' )
else                      : logger = getLogger ( __name__         )
# =============================================================================
## named tuple to DB-item: (time, payload)
Item = collections.namedtuple ( 'Item', ( 'time' , 'payload' ) )
# =============================================================================
use_bsddb3  = False

# =============================================================================
## python2 : bsddb is a part of Python
if 2 == sys.version_info.major : 

    import anydbm 
    from whichdb              import whichdb   as _whichdb

    from ostap.io.sqlitedict  import issqlite3 
    from ostap.io.sqlitedict  import SqliteDict
    
    # =====================================================================
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
    #  - Actually it is a bit extended  form of <code>whichdb.whichdb</code>
    #   that accounts for  <code>sqlite3</code>
    def whichdb ( filename  ) :
        """Guess which db package to use to open a db file.
        
        Return values:
        
        - None if the database file can't be read;
        - empty string if the file can be read but can't be recognized
        - the name of the dbm submodule (e.g. 'ndbm' or 'gnu') if recognized.
        
        Importing the given module may still fail, and opening the
        database using that module may still fail.
        
        - Actually it is a bit extended  form of `whichdb.whichdb`
        that accounts for `sqlite3`
        """

        tst = _whichdb ( filename )
        if tst or tst is None          : return tst

        if issqlite3 ( filename ) : return "sqlite3"

        return  tst 

    # =====================================================================
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
    #    accounts for <code>sqlite3</code>
    def dbopen ( file , flag = 'r' , mode=0o666 , **kwargs ):
        """Open or create database at path given by *file*.
        
        Optional argument *flag* can be 'r' (default) for read-only access, 'w'
        for read-write access of an existing database, 'c' for read-write access
        to a new or existing database, and 'n' for read-write access to a new
        database.
        
        Note: 'r' and 'w' fail if the database doesn't exist; 'c' creates it
        only if it doesn't exist; and 'n' always creates a new database.
        
        - Actually it is a bit extended  form of `dbm.open` that  accounts for `sqlite3`
        """
        
        result = whichdb ( file ) if 'n' not in flag  else None

        db = None 
        if result is None :
            
            # db doesn't exist or 'n' flag was specified to create a new db    
            if 'c' in flag or 'n' in flag:
                
                # file doesn't exist and the new flag was used so use bsddb3                     
                db = anydbm.open ( file , flag , mode )
                
            else :
                raise anydbm.error[0] ( "db file '%s' doesn't exist; use 'c' or 'n' flag to create a new db" % file )
        
        elif result == 'sqlite3' :
            
            db = SqliteDict ( filename = file , flag = flag , **kwargs )

        ## use ANYDBM 
        if db is None : db = anydbm.open ( file , flag , mode )
            
        logger.debug ("Open DBASE %s of type %s/%s" % ( file , whichdb ( file ) , type ( db ) ) ) 
        return db 
        

else :                              ## 3.3 <= python

    
    ## for python3 <code>bsddb</code> is not a part of the standard library
    ##  make a try to use <code>bsddb3</code>

    if sys.version_info < (3,3) :
        
        bsddb3     = None
        use_bsddb3 = None
        
    else  :

        try :
            import bsddb3
            use_bsddb3 = True 
        except ImportError :
            bsddb3     = None
            use_bsddb3 = False

        
    from ostap.io.sqlitedict  import issqlite3 
    from ostap.io.sqlitedict  import SqliteDict

    if bsddb3 and use_bsddb3 :
        
        ## <code>bsddb3</code> is available, try to use it as a defauld database 

        import dbm, io, struct 
        
        # =====================================================================
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
        #   that accounnt for  <code>bsddb3</code> and <code>sqlite3</code>
        def whichdb ( filename  ) :
            """Guess which db package to use to open a db file.

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
            tst = dbm.whichdb ( filename  )
            
            ## identified or non-existing DB  ? 
            if tst or tst is None : return tst

            ## non-identified DB  
            
            ## check for bsddb magic numbers (from python2)
            try : 
                with io.open ( filename  ,'rb' ) as f :
                    # Read the start of the file -- the magic number
                    s16 = f.read(16)
            except OSError :
                return None
            
            s = s16[0:4]
            
            # Return "" if not at least 4 bytes
            if len(s) != 4:
                return ""
            
            # Convert to 4-byte int in native byte order -- return "" if impossible
            try:
                ( magic, ) = struct.unpack("=l", s)
            except struct.error:
                return ""
            
            # Check for GNU dbm
            if magic in (0x13579ace, 0x13579acd, 0x13579acf):
                return "dbm.gnu"
            
            # Check for old Berkeley db hash file format v2
            if magic in (0x00061561, 0x61150600):
                return "bsddb185"
            
            # Later versions of Berkeley db hash file have a 12-byte pad in
            # front of the file type
            try:
                (magic,) = struct.unpack("=l", s16[-4:])
            except struct.error:
                return ""
            
            # Check for BSD hash
            if magic in (0x00061561, 0x61150600):
                return "bsddb3"

            if issqlite3 ( filename ) : return 'sqlite3'
                
            # Unknown
            return ""

        # =====================================================================
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
        #    accounts for <code>bsbdb3</code> and <code>sqlite3</code>
        def dbopen ( file , flag = 'r' , mode=0o666 , **kwargs ):
            """Open or create database at path given by *file*.
            
            Optional argument *flag* can be 'r' (default) for read-only access, 'w'
            for read-write access of an existing database, 'c' for read-write access
            to a new or existing database, and 'n' for read-write access to a new
            database.
            
            Note: 'r' and 'w' fail if the database doesn't exist; 'c' creates it
            only if it doesn't exist; and 'n' always creates a new database.
            
            - Actually it is a bit extended  form of `dbm.open` that  accounts for `bsddb3` and `sqlite3`
            """
            
            result = whichdb ( file ) if 'n' not in flag  else None

            db = None 
            if result is None :
                
                # db doesn't exist or 'n' flag was specified to create a new db
                if 'c' in flag or 'n' in flag:
                    
                    # file doesn't exist and the new flag was used so use bsddb3 
                    db = bsddb3.hashopen ( file , flag , mode ) 

                else :
                    
                    raise dbm.error[0] ( "db file '%s' doesn't exist; use 'c' or 'n' flag to create a new db" % file )
            
            elif result in ( 'bsddb' , 'dbhash' , 'bsddb3' , 'bsddb185' ) :
                
                db = bsddb3.hashopen ( file , flag , mode ) 

            elif result == 'sqlite3' :
                
                db = SqliteDict ( filename = file , flag = flag , *kwargs )

            ## use DBM
            if db is None : db = dbm.open ( file , flag , mode )  

            logger.debug ("Open DBASE %s of type %s/%s" % ( file , whichdb ( file ) , type ( db ) ) ) 
            return db 

    else :
        
        import dbm
        
        # =====================================================================
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
        #   that accounnt for <code>sqlite3</code?
        def whichdb ( filename  ) :
            """Guess which db package to use to open a db file.
            
            Return values:
            
            - None if the database file can't be read;
            - empty string if the file can be read but can't be recognized
            - the name of the dbm submodule (e.g. 'ndbm' or 'gnu') if recognized.
            
            Importing the given module may still fail, and opening the
            database using that module may still fail.
            
            - Actually it is a bit extended  form of `dbm.whichdb`
            that accounts for `sqlite3`
            """
            
            tst = dbm.whichdb ( filename )
            if tst or tst is None : return tst
            
            if issqlite3 ( filename ) : return "sqlite3"
        
            return  tst 
        
        # =====================================================================
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
        #    accounts for <code>sqlite3</code>
        def dbopen ( file , flag = 'r' , mode=0o666 , **kwargs ):
            """Open or create database at path given by *file*.
            
            Optional argument *flag* can be 'r' (default) for read-only access, 'w'
            for read-write access of an existing database, 'c' for read-write access
            to a new or existing database, and 'n' for read-write access to a new
            database.
            
            Note: 'r' and 'w' fail if the database doesn't exist; 'c' creates it
            only if it doesn't exist; and 'n' always creates a new database.
            
            - Actually it is a bit extended  form of `dbm.open` that  accounts for `sqlite3`
            """
            
            result = whichdb ( file ) if 'n' not in flag  else None

            db = None 
            if result is None :
                
                # db doesn't exist or 'n' flag was specified to create a new db
                
                if 'c' in flag or 'n' in flag:
                    
                    # file doesn't exist and the new flag was used so use bsddb3                     
                    db = dbm.open ( file , flag , mode )  

                else :
                    
                    raise dbm.error[0] ( "db file '%s' doesn't exist; use 'c' or 'n' flag to create a new db" % file )
            
            elif result == 'sqlite3' :
                
                db = SqliteDict ( filename = name , flag = flag , **kwargs )

            ##  use DBM  
            if db is None : db = dbm.open ( file , flag , mode )
        
            logger.debug ("Open DBASE %s of type %s/%s" % ( file , whichdb ( file ) , type ( db ) ) ) 
            return db 

        
# =============================================================================
## get disk size of data-base-like object
#  @code
#  num, size = dbsize ( 'mydb' ) 
#  @endcode  
def dbsize  ( filename  ) :
    """Get disk  size of data-base=like object
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

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
        
        
        
        
    
    
    
    
