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
from   ostap.io.zipshelve import ZipShelf
from   ostap.io.dbase     import TmpDB, ordered_dict 
import zlib
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.sqliteshelve' )
else                      : logger = getLogger ( __name__ )
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
class SQLiteShelf(ZipShelf):
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
    def __init__(
            self                                    ,
            filename                                ,
            mode         = 'c'                      ,
            tablename    = 'ostap'                  ,
            journal_mode = 'DELETE'                 ,              ## 'OFF' can an alternative 
            timeout      = 30                       , **kwargs ) : ## other properties 
        
        ## initialize base class 
        ZipShelf.__init__ ( self                        ,
                            filename                    ,
                            mode         = mode         ,
                            dbtype       = 'sqlite3'    , ## ATTENTION!!! 
                            ##                   
                            tablename    = tablename    ,
                            journal_mode = journal_mode , 
                            timeout      = timeout      , **kwargs )
        
# =============================================================================
## open new SQLiteShelve data base
#  @code
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
    def __init__ ( self            ,
                   remove  = True  ,
                   keep    = False , **kwargs ) :
        
        ## initialize the base: generate the name 
        TmpDB.__init__ ( self , suffix = '.sqldb' , remove = remove , keep = keep ) 
        
        ## open DB  
        SQLiteShelf.__init__ ( self                   ,
                               self.tmp_name          ,
                               mode         = 'c'     ,
                               writeback    = False   , **kwargs )

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
def tmpdb ( *args , **kwargs ):
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
