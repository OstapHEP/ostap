#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file bz2shelve.py
# 
# This is ``bzip2''-version of shelve database.
# 
# Keeping the same interface and functionlity as shelve data base,
# ZipShelf allows much more compact file size through the on-flight
# compression of the content
#
# The actual code has been inspired by <c>zipshelve</c> ( see Google...)
#
# However is contains several new features:
# 
#  - Optionally it is possible to perform the compression
#    of the whole data base, that can be rathe useful fo data base
#    with large amout of keys 
#
# The module has been developed and used with great success in
# ``Kali, framework for fine calibration of LHCb Electormagnetic Calorimeter''
#
# Create new DB:
#
# @code
#
# >>> import bz2shelve  
# >>> db = bz2shelve.open ('a_db', 'n')    ## create new DB
# ...
# >>> abcde = ...
# >>> db['some_key'] =  abcde              ## add information to DB
# ...
# >>> db.close()
#
# @endcode 
#
# Access to DB in read-only mode :
#
# @code
#
# >>> import bz2shelve  
# >>> db = bz2shelve.open ('a_db' , 'r' ) ## access existing dbase in read-only mode
# ...
# >>> for key in db : print(key)
# ...
# >>> abcd = db['some_key']
#
# @endcode 
#
# Access existing DB in update mode :
#
# @code
#
# >>> import bz2helve 
# >>> db = bz2shelve.open ('a_db' )    ## access existing dbase in update mode
# ...
# >>> for key in db : print key
# ...
# >>> abcd = db['some_key']
#
# @endcode 
#
# @attention: In case DB-name has extension ``.bz2'', the whole data base
#             will be ``bzip2''-ed ". 
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
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2010-04-30
# =============================================================================
"""This is ``bzip2''-version of shelve database.

Keeping the same interface and functionlity as shelve data base,
Bz2Shelf allows much more compact file size through the on-flight
compression of the content

The actual code has been inspired by zipshelve ( see Google...)

However is contains several new features:

 - Optionally it is possible to perform the compression
   of the whole data base, that can be rathe useful fo data base
   with large amout of keys 
   
The module has been developed and used with great success in
 ``Kali, framework for fine calibration of LHCb Electormagnetic Calorimeter''

 Create new DB:

 >>> import bz2shelve as DBASE  
 >>> db = DBASE.open ('a_db', 'n')    ## create new DB
 ...
 >>> abcde = ...
 >>> db['some_key'] =  abcde          ## add information to DB
 ...
 >>> db.close()

 Access to DB in read-only mode :

 >>> import bz2shelve as DBASE  
 >>> db = DBASE.open ('a_db' , 'r' ) ## access existing dbase in read-only mode
 ...
 >>> for key in db : print(key)
 ...
 >>> abcd = db['some_key']

 Access existing DB in update mode :

 >>> import bz2shelve as DBASE
 >>> db = DBASE.open ('a_db' )    ## access existing dbase in update mode
 ...
 >>> for key in db : print(key)
 ...
 >>> abcd = db['some_key']
 
 In case DB-name has extension ``.bz2'', the whole data base will be ``bzip2''-ed
 
 Attention: When one tries to read the database with pickled ROOT object using newer
 version of ROOT, one could get a ROOT read error,
 in case of evoltuion in ROOT streamers for some  classes, e.g. ROOT.TH1D
 > Error in <TBufferFile::ReadClassBuffer>: Could not find the StreamerInfo for version 2 of the class TH1D, object skipped at offset 19
 > Error in <TBufferFile::CheckByteCount>: object of class TH1D read too few bytes: 2 instead of 878
 The solution is simple and described in  file ostap.io.dump_root
 - see ostap.io.dump_root

"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2010-04-30"
__version__ = "$Revision$" 
# =============================================================================
__all__ = (
    'Bz2Shelf' ,   ## The DB-itself
    'open'     ,   ## helper function to hide the actual DB
    'tmpdb'    ,   ## create TEMPORARY data base 
    )
# =============================================================================
from sys import version_info as python_version 
# =============================================================================
import os, sys, shelve, shutil
import bz2         ## use bz2 to compress DB-content 
from   ostap.io.compress_shelve import CompressShelf, ENCODING, PROTOCOL, HIGHEST_PROTOCOL
from   ostap.io.dbase           import TmpDB 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.bz2shelve' )
else                      : logger = getLogger ( __name__             )
# =============================================================================
logger.debug ( "Simple generic (c)Pickle-based ``bzip2''-database"    )
# =============================================================================
## @class Bz2Shelf
#  ``Bzip2''-version of ``shelve''-database
#    Modes: 
#    - 'r' Open existing database for reading only
#    - 'w' Open existing database for reading and writing
#    - 'c' Open database for reading and writing, creating if it does not  exist (default)
#    - 'n' Always create a new, empty database, open for reading and writing
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
class Bz2Shelf(CompressShelf):
    """Zipped-version of ``shelve''-database
    Modes: 
    - 'r'  Open existing database for reading only
    - 'w'  Open existing database for reading and writing
    - 'c'  Open database for reading and writing, creating if it does not exist
    - 'n'  Always create a new, empty database, open for reading and writing
    """ 
    ## the known "standard" extensions: 
    extensions = '.tbz' , '.tbz2' , '.bz2' 
    ## 
    def __init__(
            self                                   ,
            filename                               ,
            mode        = 'c'                      ,
            dbtype      = ''                       , 
            protocol    = PROTOCOL                 , 
            compress    = 9                        ,
            writeback   = False                    ,
            silent      = False                    ,
            keyencoding = ENCODING                 , **kwargs ) :
        
        ## save arguments for pickling....
        self.__init_args = ( filename  ,
                             mode      ,
                             dbtype    , 
                             protocol  ,
                             compress  ,
                             writeback ,
                             silent    )

        ## initialize the base class 
        CompressShelf.__init__ ( self                      ,
                                 filename                  ,
                                 mode        = mode        ,
                                 dbtype      = dbtype      , 
                                 protocol    = protocol    ,
                                 compress    = compress    , 
                                 writeback   = writeback   ,
                                 silent      = silent      ,
                                 keyencoding = keyencoding , **kwargs ) 
        
    ## needed for proper (un)pickling 
    def __getinitargs__ ( self ) :
        """for proper (un_pickling"""
        return self.__init_args

    ## needed for proper (un)pickling 
    def __getstate__ ( self ) :
        """for proper (un)pickling"""
        self.sync() 
        return {}
    
    ## needed for proper (un)pickling 
    def __setstate__ ( self , dct ) :
        """for proper (un)pickling"""
        pass
    
    # =========================================================================
    ## compress (bz2) the file into temporary location, keep original
    def compress_files   ( self , files ) :
        """Compress (bz2) the file into temporary location, keep original
        """
        output = self.tempfile()
        
        import tarfile
        with tarfile.open ( output , 'w:bz2' ) as tfile :
            for file in files  :
                _ , name = os.path.split ( file )
                tfile.add ( file , name  )
        ##
        return output 

    # =========================================================================
    ## uncompress (bz2) the file into temporary location, keep original
    def uncompress_file ( self , filein ) :
        """Uncompress (bz2) the file into temporary location, keep original
        """
        items  = []
        tmpdir = self.tempdir ()
        
        ## 2) try compressed-tarfile 
        import tarfile
        if tarfile.is_tarfile ( filein ) : 
            with tarfile.open ( filein  , 'r:*' ) as tfile :
                for item in tfile  :
                    tfile.extract ( item , path = tmpdir )
                    items.append  ( os.path.join ( tmpdir , item.name ) )
                items.sort() 
                return tuple ( items )
            
        import tempfile , io   
        fd , fileout = tempfile.mkstemp ( prefix = 'ostap-tmp-' , suffix = '-bz2db' )
        with bz2.open ( filein  , 'rb' ) as fin : 
            with io.open ( fileout , 'wb' ) as fout : 
                shutil.copyfileobj ( fin , fout )                
                return fileout , 
            
    # ==========================================================================
    ## compress (bzip2)  the item  using <code>bz2.compress</code>
    def compress_item ( self , value ) :
        """Compress (zip) the item using ``bz2.compress''
        - see bz2.compress
        """
        return bz2.compress ( self.pickle ( value ) , self.compresslevel )

    # =========================================================================
    ## uncompres (bzip2) the item using <code>bz2.decompress</code>
    def uncompress_item ( self , value ) :
        """Uncompress (bzip2) the item using ``bz2.decompress''
        -  see bz2.decompress
        """        
        return self.unpickle ( bz2.decompress ( value ) ) 

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
        new_db = Bz2Shelf ( new_name                         ,
                            mode        =  'c'               ,
                            dbtype      = self.dbtype        , 
                            protocol    = self.protocol      ,
                            compress    = self.compresslevel , 
                            writeback   = self.writeback     ,
                            silent      = self.silent        ,
                            keyencoding = self.keyencoding   )
        
        ## copy the content
        copy = keys if keys else self.keys()
        for key in copy : new_db [ key ] = self [ key ]
        new_db.sync ()  
        return new_db 
                         
# =============================================================================
## helper function to access Bz2Shelve data base
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
def open ( filename                  ,
           mode          = 'c'       ,
           dbtype        = ''        , 
           protocol      = PROTOCOL  ,
           compress      = 9         , 
           writeback     = False     ,
           silent        = True      ,
           keyencoding   = ENCODING  , **kwargs ) :
    
    """Open a persistent dictionary for reading and writing.
    
    The filename parameter is the base filename for the underlying
    database.  As a side-effect, an extension may be added to the
    filename and more than one file may be created.  The optional flag
    parameter has the same interpretation as the flag parameter of
    anydbm.open(). The optional protocol parameter specifies the
    version of the pickle protocol (0, 1, or 2).
    
    See the module's __doc__ string for an overview of the interface.
    """
    
    return Bz2Shelf ( filename                  ,
                      mode        = mode        ,
                      dbtype      = dbtype      ,  
                      protocol    = protocol    ,
                      compress    = compress    ,
                      writeback   = writeback   ,
                      silent      = silent      ,
                      keyencoding = keyencoding , **kwargs )

# =============================================================================
## @class TmpBz2Shelf
#  TEMPORARY bzip2-version of ``shelve''-database
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2015-10-31
class TmpBz2Shelf(Bz2Shelf,TmpDB):
    """
    TEMPORARY ``bzip2''-version of ``shelve''-database     
    """    
    def __init__( self                           ,
                  dbtype      = ''               , 
                  protocol    = HIGHEST_PROTOCOL , 
                  compress    = 9                ,
                  silent      = False            ,
                  keyencoding = ENCODING         ,
                  remove      = True             ,
                  keep        = False            , **kwargs ) :

        ## initialize the base: generate the name 
        TmpDB.__init__ ( self , suffix = '.bz2db' , remove = remove , keep = keep ) 
        
        ## open DB 
        Bz2Shelf.__init__ ( self                      ,  
                            self.tmp_name             ,
                            dbtype      = dbtype      , 
                            mode        = 'c'         ,
                            protocol    = protocol    ,
                            compress    = compress    , 
                            writeback   = False       , ## writeback 
                            silent      = silent      ,
                            keyencoding = keyencoding , **kwargs ) 
        
    ## close and delete the file 
    def close ( self )  :
        ## close the shelve file
        Bz2Shelf.close ( self )
        ## delete the file
        TmpDB.clean    ( self )
            
# =============================================================================
## helper function to open TEMPORARY ZipShelve data base#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
def tmpdb ( dbtype      = ''               ,
            protocol    = HIGHEST_PROTOCOL ,
            compress    = 9                , 
            silent      = True             ,
            keyencoding = ENCODING         ,
            remove      = True             ,   ## immediate remove 
            keep        = False            , **kwargs ) : ## keep it 
    """Open a TEMPORARY persistent dictionary for reading and writing.
    
    The optional protocol parameter specifies the
    version of the pickle protocol (0, 1, or 2).
    
    See the module's __doc__ string for an overview of the interface.
    """
    return TmpBz2Shelf ( dbtype      = dbtype      ,
                         protocol    = protocol    ,
                         compress    = compress    ,
                         silent      = silent      ,
                         keyencoding = keyencoding ,
                         remove      = remove      ,
                         keep        = keep        , **kwargs ) 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
