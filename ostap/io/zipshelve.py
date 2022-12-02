#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file zipshelve.py
# 
# This is zip-version of shelve database.
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
# >>> import zipshelve  
# >>> db = zipshelve.open ('a_db', 'n')    ## create new DB
# ...
# >>> abcde = ...
# >>> db['some_key'] =  abcde              ## add information to DB
# ...
# >>> db.close()
# @endcode 
#
# Access to DB in read-only mode :
#
# @code
# >>> import zipshelve  
# >>> db = zipshelve.open ('a_db' , 'r' )    ## access existing dbase in read-only mode
# ...
# >>> for key in db : print(key)
# ...
# >>> abcd = db['some_key']
# @endcode 
#
# Access existing DB in update mode :
#
# @code
# >>> import ziphelve  
# >>> db = zipshelve.open ('a_db' )    ## access existing dbase in update mode
# ...
# >>> for key in db : print key
# ...
# >>> abcd = db['some_key']
# @endcode 
#
# @attention: In case DB-name has extension ``.gz'', the whole data base
#             will be ``gzip''-ed. 
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
# 
# =============================================================================
"""This is zip-version of shelve database.

Keeping the same interface and functionlity as shelve data base,
ZipShelf allows much more compact file size through the on-flight
compression of the content

The actual code has been inspired by zipshelve ( see Google...)

However is contains several new features:

 - Optionally it is possible to perform the compression
   of the whole data base, that can be rathe useful fo data base
   with large amout of keys 
   
The module has been developed and used with great success in
 ``Kali, framework for fine calibration of LHCb Electormagnetic Calorimeter''

 Create new DB:

 >>> import zipshelve as DBASE  
 >>> db = DBASE.open ('a_db', 'n')      ## create new DB
 ...
 >>> abcde = ...
 >>> db['some_key'] =  abcde            ## add information to DB
 ...
 >>> db.close()

 Access to DB in read-only mode :

 >>> import zipshelve as DBASE  
 >>> db = DBASE.open ('a_db' , 'r' )    ## access existing dbase in read-only mode
 ...
 >>> for key in db : print(key)
 ...
 >>> abcd = db['some_key']

 Access existing DB in update mode :

 >>> import zipshelve as DBASE 
 >>> db = DBASE.open ('a_db' )          ## access existing dbase in update mode
 ...
 >>> for key in db : print(key)
 ...
 >>> abcd = db['some_key']
 
 In case DB-name has extension ``.gz'', the whole data base will be ``gzip''-ed
 
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
    'ZipShelf' ,   ## The DB-itself
    'open'     ,   ## helper function to hide the actual DB
    'tmpdb'    ,   ## create TEMPORARY data base 
    )
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.zipshelve' )
else                      : logger = getLogger ( __name__             )
# =============================================================================
logger.debug ( "Simple generic (c)Pickle-based ``zipped''-database"   )
# =============================================================================
from sys import version_info as python_version 
# =============================================================================
## try:
##     from cPickle   import Pickler, Unpickler
## except ImportError:
##     from  pickle   import Pickler, Unpickler
## # =============================================================================
## if 2 < python_version.major :
##     from io import BytesIO
## else : 
##     try:
##         from cStringIO import StringIO as BytesIO 
##     except ImportError:
##         from  StringIO import StringIO as BytesIO    
# ==============================================================================
import os, sys, shelve, shutil
import zlib ## use zlib to compress DB-content 
from   ostap.io.compress_shelve import CompressShelf, ENCODING, PROTOCOL, HIGHEST_PROTOCOL
from   ostap.io.dbase           import TmpDB 
# =============================================================================
## @class ZipShelf
#  Zipped-version of ``shelve''-database
#    Modes: 
#    - 'r' Open existing database for reading only
#    - 'w' Open existing database for reading and writing
#    - 'c' Open database for reading and writing, creating if it does not  exist (default)
#    - 'n' Always create a new, empty database, open for reading and writing
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
class ZipShelf(CompressShelf):
    """Zipped-version of ``shelve''-database
    Modes: 
    - 'r'  Open existing database for reading only
    - 'w'  Open existing database for reading and writing
    - 'c'  Open database for reading and writing, creating if it does not exist
    - 'n'  Always create a new, empty database, open for reading and writing
    """ 
    ## the known "standard" extensions: 
    extensions = '.zip' , '.tgz' , '.gz'  
    ## 
    def __init__(
        self                                   ,
        filename                               ,
        mode        = 'c'                      , 
        protocol    = PROTOCOL                 , 
        compress    = zlib.Z_BEST_COMPRESSION  ,
        writeback   = False                    ,
        silent      = False                    ,
        keyencoding = ENCODING                 ) :


        ## save arguments for pickling....
        self.__init_args = ( filename  ,
                             mode      ,
                             protocol  ,
                             compress  ,
                             writeback ,
                             silent    )

        ## initialize the base class 
        CompressShelf.__init__ ( self        ,
                                 filename    ,
                                 mode        ,
                                 protocol    ,
                                 compress    , 
                                 writeback   ,
                                 silent      ,
                                 keyencoding ) 
        
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
    ## compress the file into temporary location, keep original
    def compress_files   ( self , files ) :
        """Compress the files into the temporary location, keep original
        """
        output = self.tempfile ()
        
        import zipfile 
        with zipfile.ZipFile( output , 'w' , allowZip64 = True ) as zfile :
            for file in files :
                _ , name = os.path.split ( file )
                zfile.write ( file  , name  )
                
        return output 

    # =========================================================================
    ## uncompress (gunzip) the file into temporary location, keep original
    #  @code
    #  db    = ...
    #  files = db.uncompress_file ( input_cmpressed_file )   
    #  @endcode 
    def uncompress_file ( self , filein ) :
        """Uncompress (gunzip) the file into temporary location, keep original
        >>> db    = ...
        >>> files = db.uncompress_file ( input_cmpressed_file )   
        """
        items  = []
        tmpdir = self.tempdir ()
        
        ## 1) try zipfile 
        import zipfile
        if zipfile.is_zipfile ( filein ) :
            with zipfile.ZipFile ( filein , 'r' , allowZip64 = True ) as zfile :
                for item in zfile.filelist :
                    zfile.extract ( item , path = tmpdir )
                    items.append  ( os.path.join ( tmpdir , item.filename ) )
            items.sort() 
            return tuple  ( items ) 
                    
        ## 2) try compressed-tarfile 
        import tarfile
        if tarfile.is_tarfile ( filein ) :
            with tarfile.open ( filein  , 'r:*' ) as tfile :
                for item in tfile  :
                    tfile.extract ( item , path = tmpdir )
                    items.append  ( os.path.join ( tmpdir , item.name ) )
            items.sort() 
            return tuple ( items ) 

        ## 3) try old good gzipped (single) file
        import gzip , io, tempfile
        fd , fileout = tempfile.mkstemp ( prefix = 'ostap-tmp-' , suffix = '-db' )
        with gzip.open ( filein  , 'rb' ) as fin : 
            with io.open ( fileout , 'wb' ) as fout : 
                shutil.copyfileobj ( fin , fout )            
                return fileout , 
            
    # ==========================================================================
    ## compress (zip)  the item  using <code>zlib.compress</code>
    def compress_item ( self , value ) :
        """Compress (zip) the item using ``zlib.compress''
        - see zlib.compress
        """
        ## f = BytesIO ()
        ## p = Pickler ( f , self.protocol )
        ## p.dump ( value )
        ## return zlib.compress ( f.getvalue() , self.compresslevel )
        return zlib.compress ( self.pickle ( value ) , self.compresslevel )
        
    # =========================================================================
    ## uncompress (unzip) the item using <code>zlib.decompress</code>
    def uncompress_item ( self , value ) :
        """Uncompress (nuzip) the item using ``zlib.decompress''
        -  see zlib.decompress
        """        
        ## f = BytesIO ( zlib.decompress ( value ) )
        ## return Unpickler ( f ) . load ( )
        return self.unpickle ( zlib.decompress ( value ) ) 

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
        new_db = ZipShelf ( new_name                         ,
                            mode        =  'c'               ,
                            protocol    = self.protocol      ,
                            compress    = self.compresslevel , 
                            writeback   = self.writeback     ,
                            silent      = self.silent        ,
                            keyencoding = self.keyencoding   )
        
        ## copy the content
        if keys :
            for key in self.keys() :
                if key in keys     : new_db [ key ] = self [ key ]
        else : 
            for key in self.keys() : new_db [ key ] = self [ key ]
            
        new_db.sync ()  
        return new_db 
    
# =============================================================================
## helper function to access ZipShelve data base
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
def open ( filename                                 ,
           mode          = 'c'                      ,
           protocol      = PROTOCOL                 ,
           compresslevel = zlib.Z_BEST_COMPRESSION  , 
           writeback     = False                    ,
           silent        = True                     ,
           keyencoding   = ENCODING                 ) :
    
    """Open a persistent dictionary for reading and writing.
    
    The filename parameter is the base filename for the underlying
    database.  As a side-effect, an extension may be added to the
    filename and more than one file may be created.  The optional flag
    parameter has the same interpretation as the flag parameter of
    anydbm.open(). The optional protocol parameter specifies the
    version of the pickle protocol (0, 1, or 2).
    
    See the module's __doc__ string for an overview of the interface.
    """
    
    return ZipShelf ( filename      ,
                      mode          ,
                      protocol      ,
                      compresslevel ,
                      writeback     ,
                      silent        ,
                      keyencoding   )

# =============================================================================
## @class TmpZipShelf
#  TEMPORARY Zipped-version of ``shelve''-database
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2015-10-31
class TmpZipShelf(ZipShelf,TmpDB):
    """TEMPORARY Zipped-version of ``shelve''-database     
    """    
    def __init__( self                                   ,
                  protocol    = HIGHEST_PROTOCOL         , 
                  compress    = zlib.Z_BEST_COMPRESSION  ,
                  silent      = False                    ,
                  keyencoding = ENCODING                 , 
                  remove      = True                     ,   ## immediate remove 
                  keep        = False                    ) : ## keep it 
        
        ## initialize the base: generate the name 
        TmpDB.__init__ ( self , suffix = '.zdb' , remove = remove , keep = keep ) 
         
        ZipShelf.__init__ ( self          ,  
                            self.tmp_name ,
                            'c'           ,
                            protocol      ,
                            compress      , 
                            False         , ## writeback 
                            silent        ,
                            keyencoding   )
        
    ## close and delete the file 
    def close ( self )  :
        ZipShelf.close ( self )
        TmpDB.clean    ( self ) 
            
            
# =============================================================================
## helper function to open TEMPORARY ZipShelve data base#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
def tmpdb ( protocol      = HIGHEST_PROTOCOL        ,
            compresslevel = zlib.Z_BEST_COMPRESSION , 
            silent        = True                    ,
            keyencoding   = ENCODING                ,
            remove        = True                    ,    ## immediate remove 
            keep          = False                   ) :  ## keep it 
    """Open a TEMPORARY persistent dictionary for reading and writing.
    
    The optional protocol parameter specifies the
    version of the pickle protocol (0, 1, or 2).
    
    See the module's __doc__ string for an overview of the interface.
    """
    return TmpZipShelf ( protocol      ,
                         compresslevel ,
                         silent        ,
                         keyencoding   ,
                         remove        ,
                         keep          ) 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
