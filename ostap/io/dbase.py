#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file compress_dbase.py
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
    'whichdb'  , ## guess database type  
    'dbopen'   , ## open database 
    )
# =============================================================================
import sys
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.compress_shelve' )
else                      : logger = getLogger ( __name__                   )
# =============================================================================
if  sys.version_info.major < 3 :   ## PYTHON2 

    ## for python2 <code>bdsdb</code> is part of the standard library 
    from anydbm  import open    as dbopen
    from whichdb import whichdb 

else :                              ## PYTHON3   

    
    ## for python3 <code>bdsdb</code> is not a part of the standard library
    ##  make a try to use <code>bdsdb3</code>
    
    
    try :
        import bdsdb3
    except ImportError :
        bsddb3 = None

    if bdsdb3 :
        ## <code>bdsdb3</code> is available, try to use it as a defauld database 

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
        #   that accounnt for  <code>bdsdb3</code>
        def whichdb ( filename  ) :
            """Guess which db package to use to open a db file.

            Return values:
            
            - None if the database file can't be read;
            - empty string if the file can be read but can't be recognized
            - the name of the dbm submodule (e.g. 'ndbm' or 'gnu') if recognized.
            
            Importing the given module may still fail, and opening the
            database using that module may still fail.
            
            - Actually it is a bit extended  form of `dbm.whichdb`
            that accounts for `bdsdb3`
            """
            
            ## use the standard function 
            tst = dbm.wichdb ( filename  )
            
            ## identified or non-existing DB  ? 
            if tst or tsts is None : return tst

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
                return "dbhash"
            
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
        #    accounts for <code>bsbdb3</code>
        def dbopen ( file , flag = 'r' , mode=0o666 ):
            """Open or create database at path given by *file*.
            
            Optional argument *flag* can be 'r' (default) for read-only access, 'w'
            for read-write access of an existing database, 'c' for read-write access
            to a new or existing database, and 'n' for read-write access to a new
            database.
            
            Note: 'r' and 'w' fail if the database doesn't exist; 'c' creates it
            only if it doesn't exist; and 'n' always creates a new database.
            
            - Actually it is a bit extended  form of `dbm.open` that  accounts for `bdsdb3`
            """
            
            result = whichdb ( file  ) if 'n' not in flag  else None
            
            if result is None :
                
                # db doesn't exist or 'n' flag was specified to create a new db
                
                if 'c' in flag or 'n' in flag:
                    
                    # file doesn't exist and the new flag was used so use bdsdb3 
                    
                    return bdsdb3.hasopen ( flag , mode ) 
                
                raise dbm.error[0] ( "db file '%s' doesn't exist; use 'c' or 'n' flag to create a new db" % file )
            
            elif result in ( 'bdsdb' , 'dbhash' , 'bdsdb3' , 'dbddb185' ) :
                
                return bdsdb3.hasopen ( flag , mode ) 
            
            return dbm.open ( file , flag , mode )  

        
    else :
        
        from dbm import open     as dbopen
        from dbm import whichdb
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
        
        
        
        
    
    
    
    
