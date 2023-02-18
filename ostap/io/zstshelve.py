#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file zstshelve.py
# 
# This is `zstandard'-version of shelve database.
# 
# Keeping the same interface and functionlity as shelve data base,
# ZstShelf allows much more compact file size through the on-flight
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
# `Kali, framework for fine calibration of LHCb Electormagnetic Calorimeter'
#
# Create new DB:
#
# @code
#
# >>> import zstshelve  
# >>> db = zstshelve.open ('a_db', 'n')    ## create new DB
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
# >>> import zstshelve  
# >>> db = zstshelve.open ('a_db' , 'r' ) ## access existing dbase in read-only mode
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
# >>> import zstshelve 
# >>> db = zstshelve.open ('a_db' )    ## access existing dbase in update mode
# ...
# >>> for key in db : print key
# ...
# >>> abcd = db['some_key']
#
# @endcode 
#
# @attention: In case DB-name has extension `.zst'  the whole data base
#             will be `ZST'-ed ". 
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
"""This is `zstandard'-version of shelve database.

Keeping the same interface and functionlity as shelve data base,
ZstShelf allows much more compact file size through the on-flight
compression of the content

The actual code has been inspired by zipshelve ( see Google...)

However is contains several new features:

 - Optionally it is possible to perform the compression
   of the whole data base, that can be rathe useful fo data base
   with large amout of keys 
   
The module has been developed and used with great success in
 `Kali, framework for fine calibration of LHCb Electormagnetic Calorimeter'

 Create new DB:

 >>> import zsthelve as DBASE  
 >>> db = DBASE.open ('a_db', 'n')    ## create new DB
 ...
 >>> abcde = ...
 >>> db['some_key'] =  abcde          ## add information to DB
 ...
 >>> db.close()

 Access to DB in read-only mode :

 >>> import zstshelve as DBASE  
 >>> db = DBASE.open ('a_db' , 'r' ) ## access existing dbase in read-only mode
 ...
 >>> for key in db : print(key)
 ...
 >>> abcd = db['some_key']

 Access existing DB in update mode :

 >>> import zstshelve as DBASE
 >>> db = DBASE.open ('a_db' )    ## access existing dbase in update mode
 ...
 >>> for key in db : print(key)
 ...
 >>> abcd = db['some_key']
 
 In case DB-name has extension `.zst' the whole data base will be `ZST'-ed

 Attention: When one tries to read the database with pickled ROOT object using newer
 version of ROOT, one could get a ROOT read error,
 in case of evolution in ROOT streamers for some  classes, e.g. ROOT.TH1D
 > Error in <TBufferFile::ReadClassBuffer>: Could not find the StreamerInfo for version 2 of the class TH1D, object skipped at offset 19
 > Error in <TBufferFile::CheckByteCount>: object of class TH1D read too few bytes: 2 instead of 878
 The solution is simple and described in  file ostap.io.dump_root
 - see ostap.io.dump_root
 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2023-02-17"
__version__ = "$Revision:$" 
# =============================================================================
__all__ = ()
# =============================================================================
from sys import version_info as python_version 
# =============================================================================
from   ostap.io.compress_shelve import CompressShelf, ENCODING, PROTOCOL, HIGHEST_PROTOCOL
from   ostap.io.dbase           import TmpDB 
import os, sys, shelve, shutil
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.zstshelve' )
else                      : logger = getLogger ( __name__             )
# =============================================================================
logger.debug ( "Simple generic (c)Pickle-based `ZST'-database"    )
# =============================================================================
if ( 3 , 6 )<= python_version :
    try :
        import zstandard as zst
    except ImportError :
        logger.error ( "Cannot import 'zstandard', zstshelve is disabled" )
        zst = None
else :
    logger.error ( 'zstshelve is disabled for python %s' % str ( python_version ) )
    zst = None

# =============================================================================
if zst : 
    # =============================================================================
    ## @class ZstShelf
    #  `ZST'-version of `shelve'-database
    #    Modes: 
    #    - 'r' Open existing database for reading only
    #    - 'w' Open existing database for reading and writing
    #    - 'c' Open database for reading and writing, creating if it does not  exist (default)
    #    - 'n' Always create a new, empty database, open for reading and writing
    #  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
    #  @date   2010-04-30
    class ZstShelf(CompressShelf):
        """ZST-version of `shelve'-database
        Modes: 
        - 'r'  Open existing database for reading only
        - 'w'  Open existing database for reading and writing
        - 'c'  Open database for reading and writing, creating if it does not exist
        - 'n'  Always create a new, empty database, open for reading and writing
        """ 
        ## the known "standard" extensions: 
        extensions =  '.zst', '.zstd' 
        ## 
        def __init__( self                                   ,
                      filename                               ,
                      mode        = 'c'                      , 
                      protocol    = PROTOCOL                 , 
                      compress    = 3                        , ## level in Zstandard 
                      writeback   = False                    ,
                      silent      = False                    ,
                      keyencoding = 'utf-8'                  ,
                      threads     = -1                       ) :
            
            ## save arguments for pickling....
            self.__init_args = ( filename  ,
                                 mode      ,
                                 protocol  ,
                                 compress  ,
                                 writeback ,
                                 silent    ,
                                 threads   )
            
            self.__threads      = threads 
            self.__compressor   = zst.ZstdCompressor   ( level              = compress ,
                                                         threads            = threads  ,
                                                         write_checksum     = True     ,
                                                         write_content_size = True     ) 
            self.__decompressor = zst.ZstdDecompressor ( )
            
            
            ## initialize the base class 
            CompressShelf.__init__ ( self        ,
                                     filename    ,
                                     mode        ,
                                     protocol    ,
                                     compress    , 
                                     writeback   ,
                                     silent      ,
                                     keyencoding ) 
            
            
        @property
        def threads ( self ) :
            """'threads' : now many (C)-thread can be used for compression/decompression?"""
            return self.__threads
        @property
        def compressor ( self ) :
            """'compressor' : get the actual compressor object"""
            return self.__compressor
        @property
        def decompressor ( self ) :
            """'decompressor' : get the actual decompressor object"""
            return self.__decompressor
        
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
        ## compress (zstandard) the file into temporary location, keep original
        def compress_files ( self , files ) :
            """Compress (zstandard) the file into temporary location, keep original
            """
            output = self.tempfile()
            import tarfile
            with tarfile.open ( output , 'x:gz' ) as tfile :
                for file in files  :
                    _ , name = os.path.split ( file )
                    tfile.add ( file , name  )
                return output
            
        # =========================================================================
        ## uncompress (zstandard) the file into temporary location, keep original
        def uncompress_file ( self , filein ) :
            """Uncompress (zstandard) the file into temporary location, keep original
            """
            
            items  = []
            tmpdir = self.tempdir ()
            
            ## 1) try compressed-tarfile 
            import tarfile
            if tarfile.is_tarfile ( filein ) : 
                with tarfile.open ( filein  , 'r:*' ) as tfile :
                    for item in tfile  :
                        tfile.extract ( item , path = tmpdir )
                        items.append  ( os.path.join ( tmpdir , item.name ) )
                        items.sort() 
                    return tuple ( items )
                
            ## 2) single zst-file 
            import tempfile , io   
            fd , fileout = tempfile.mkstemp ( prefix = 'ostap-tmp-' , suffix = '-zstdb' )
            
            with io.open ( filein , 'rb' ) as fin :
                with io.open ( fileout , 'wb' ) as fout :
                    self.decompressor.copy_stream ( fin , fout )
                    return fileout ,
            
        # ==========================================================================
        ## compress (ZST)  the item  using compressor 
        def compress_item ( self , value ) :
            """Compress (ZST) the item using compressor 
            """
            return self.compressor.compress (  self.pickle ( value ) )
    
        # =========================================================================
        ## uncompres (ZST) the item using decompressor 
        def uncompress_item ( self , value ) :
            """Uncompress (ZST) the item using decompressor 
            """        
            return self.unpickle ( self.decompressor.decompress ( value ) ) 
        
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
            new_db = ZstShelf ( new_name                         ,
                                mode        =  'c'               ,
                                protocol    = self.protocol      ,
                                compress    = self.compresslevel , 
                                writeback   = self.writeback     ,
                                silent      = self.silent        ,
                                keyencoding = self.keyencoding   ,
                                threads     = self.threads       )
            
            ## copy the content
            copy = keys if keys else self.keys()
            for key in copy : new_db [ key ] = self [ key ]            
            new_db.sync ()  
            return new_db 
        
    # =============================================================================
    ## helper function to access ZstShelve data base
    #  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
    #  @date   2023-02-17
    def open ( filename                            ,
               mode          = 'c'                 ,
               protocol      = PROTOCOL            ,
               compresslevel = 3                   , 
               writeback     = False               ,
               silent        = True                ,
               keyencoding   = ENCODING            ,
               threads       = -1                  ) :
        
        """Open a persistent dictionary for reading and writing.
        
        The filename parameter is the base filename for the underlying
        database.  As a side-effect, an extension may be added to the
        filename and more than one file may be created.  The optional flag
        parameter has the same interpretation as the flag parameter of
        anydbm.open(). The optional protocol parameter specifies the
        version of the pickle protocol (0, 1, or 2).
        
        See the module's __doc__ string for an overview of the interface.
        """
        
        return ZstShelf ( filename      ,
                          mode          ,
                          protocol      ,
                          compresslevel ,
                          writeback     ,
                          silent        ,
                          keyencoding   ,
                          threads       )
    
    # =============================================================================
    ## @class TmpZstShelf
    #  TEMPORARY zst-version of `shelve'-database
    #  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
    #  @date   2023-02-17
    class TmpZstShelf(ZstShelf,TmpDB):
        """
        TEMPORARY `ZST'-version of `shelve'-database     
        """    
        def __init__( self                              ,
                      protocol    = HIGHEST_PROTOCOL    , 
                      compress    = 3                   ,
                      silent      = False               ,
                      keyencoding = ENCODING            , 
                      remove      = True                ,
                      keep        = False               ,
                      threads     = -1                  ) :
            
            ## initialize the base: generate the name 
            TmpDB.__init__ ( self , suffix = '.zstdb' , remove = remove , keep = keep ) 
            
            ## open DB 
            ZstShelf.__init__ ( self          ,  
                                self.tmp_name ,
                                'c'           ,
                                protocol      ,
                                compress      , 
                                False         , ## writeback 
                                silent        ,
                                keyencoding   ,
                                threads       ) 
            
        ## close and delete the file 
        def close ( self )  :
            ## close the shelve file
            ZstShelf.close ( self )
            ## delete the file
            TmpDB  .clean  ( self ) 
            
    # =============================================================================
    ## helper function to open TEMPORARY ZstShelve data base#
    #  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
    #  @date   2023-02-17
    def tmpdb ( protocol      = HIGHEST_PROTOCOL    ,
                compresslevel = 3                   , 
                silent        = True                ,
                keyencoding   = ENCODING            ,
                remove        = True                ,  ## immediate remove 
                keep          = False               ,  ## keep it
                threads       = -1                  ) : 
        """Open a TEMPORARY persistent dictionary for reading and writing.
        
        The optional protocol parameter specifies the
        version of the pickle protocol (0, 1, or 2).
        
        See the module's __doc__ string for an overview of the interface.
        """
        return TmpZstShelf ( protocol      ,
                             compresslevel ,
                             silent        ,
                             keyencoding   ,
                             remove        ,
                             keep          ,
                             threads       ) 

    # ==========================================================================
    __all__ = (
        'ZstShelf'    , ## database
        'open'        , ## open the database 
        'tmpdb'       , ## open the temporary database 
        )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
