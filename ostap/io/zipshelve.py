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
try:
    from cPickle   import Pickler, Unpickler, HIGHEST_PROTOCOL
except ImportError:
    from  pickle   import Pickler, Unpickler, HIGHEST_PROTOCOL
# =============================================================================
## to be compatible between  Python2 and Python3 
PROTOCOL = 2
ENCODING = 'utf-8'
# =============================================================================
if python_version.major > 2  :
    from io import BytesIO
else : 
    try:
        from cStringIO import StringIO as BytesIO 
    except ImportError:
        from  StringIO import StringIO as BytesIO    
# ==============================================================================
import os, sys
import zlib        ## use zlib to compress DB-content 
import shelve      ## 
import shutil
from ostap.io.compress_shelve import CompressShelf
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
    extensions = '.gz' , 
    ## 
    def __init__(
        self                                   ,
        filename                               ,
        mode        = 'c'                      , 
        protocol    = PROTOCOL                 , 
        compress    = zlib.Z_BEST_COMPRESSION  ,
        writeback   = False                    ,
        silent      = False                    ,
        keyencoding = 'utf-8'                  ) :


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
    ## compress (gzip) the file into temporary location, keep original
    def compress_file   ( self , filein ) :
        """Compress (gzip) the file into temporary location, keep original
        """
        import tempfile , gzip , io 
        fd , fileout = tempfile.mkstemp ( prefix = 'tmp-' , suffix = '-db.gz' )
        with io.open ( filein , 'rb' ) as fin :
            with gzip.open ( fileout , 'wb') as fout : 
                shutil.copyfileobj ( fin , fout )                
                return fileout 

    # =========================================================================
    ## uncompress (gunzip) the file into temporary location, keep original
    def uncompress_file ( self , filein ) :
        """Uncompress (gunzip) the file into temporary location, keep original
        """
        import tempfile, gzip , io   
        fd , fileout = tempfile.mkstemp ( prefix = 'tmp-' , suffix = '-db' )
        with gzip.open ( filein  , 'rb' ) as fin : 
            with io.open ( fileout , 'wb' ) as fout : 
                shutil.copyfileobj ( fin , fout )            
                return fileout
    
    # ==========================================================================
    ## compress (zip)  the item  using <code>zlib.compress</code>
    def compress_item ( self , value ) :
        """Compress (zip) the item using ``zlib.compress''
        - see zlib.compress
        """
        f = BytesIO ()
        p = Pickler ( f , self.protocol )
        p.dump ( value )
        return zlib.compress ( f.getvalue() , self.compresslevel )

    # =========================================================================
    ## uncompres (unzip) the item using <code>zlib.decompress</code>
    def uncompress_item ( self , value ) :
        """Uncompress (nuzip) the item using ``zlib.decompress''
        -  see zlib.decompress
        """        
        f = BytesIO ( zlib.decompress ( value ) )
        return Unpickler ( f ) . load ( )

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
class TmpZipShelf(ZipShelf):
    """TEMPORARY Zipped-version of ``shelve''-database     
    """    
    def __init__(
        self                                   ,
        protocol    = HIGHEST_PROTOCOL         , 
        compress    = zlib.Z_BEST_COMPRESSION  ,
        silent      = False                    ,
        keyencoding = ENCODING                 ) :

        ## create temporary file name 
        import tempfile
        filename = tempfile.mktemp  ( prefix = 'tmpdb-' , suffix = '.zdb' )
        
        ZipShelf.__init__ ( self        ,  
                            filename    ,
                            'c'         ,
                            protocol    ,
                            compress    , 
                            False       , ## writeback 
                            silent      ,
                            keyencoding ) 
        
    ## close and delete the file 
    def close ( self )  :
        ## close the shelve file
        fname = self.filename 
        ZipShelf.close ( self )
        ## delete the file 
        if os.path.exists ( fname ) :
            try :
                os.unlink ( fname )
            except : 
                pass
            
# =============================================================================
## helper function to open TEMPORARY ZipShelve data base#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
def tmpdb ( protocol      = HIGHEST_PROTOCOL        ,
            compresslevel = zlib.Z_BEST_COMPRESSION , 
            silent        = True                    ,
            keyencoding   = ENCODING                ) :
    """Open a TEMPORARY persistent dictionary for reading and writing.
    
    The optional protocol parameter specifies the
    version of the pickle protocol (0, 1, or 2).
    
    See the module's __doc__ string for an overview of the interface.
    """
    return TmpZipShelf ( protocol      ,
                         compresslevel ,
                         silent        ,
                         keyencoding   ) 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
# =============================================================================
# The END 
# =============================================================================
