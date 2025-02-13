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
__all__ = (
    'ZstShelf'    , ## database
    'TmpZstShelf' , ## database
    'open'        , ## open the database 
    'tmpdb'       , ## open the temporary database 
)
# =============================================================================
from   ostap.io.compress_shelve import CompressShelf, HIGHEST_PROTOCOL
from   ostap.io.dbase           import TmpDB 
import shelve
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.zstshelve' )
else                      : logger = getLogger ( __name__             )
# =============================================================================
logger.debug ( "Simple generic (c)Pickle-based `ZST'-database"    )
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import zstandard as zst # =================================================
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    logger.error ( "Cannot import 'zstandard', zstshelve is disabled" )
    zst = None

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
    """ ZST-version of `shelve'-database
     Modes: 
    - 'r'  Open existing database for reading only
    - 'w'  Open existing database for reading and writing
    - 'c'  Open database for reading and writing, creating if it does not exist
    - 'n'  Always create a new, empty database, open for reading and writing
    """
    ## 
    def __init__( self                                   ,
                  dbname                                 ,
                  mode        = 'c'                      ,
                  compress    = 22                       , ## level in Zstandard 
                  threads     = -1                       , **kwargs ) :
        
        assert zst , "`zstandard` module is not available!"
        assert 1 <= compress <= 22 , 'Invalid `compress` for `zstandard`-compression: %s' % compress
        
        self.__threads      = threads 
        self.__compressor   = zst.ZstdCompressor   ( level              = compress ,
                                                     threads            = threads  ,
                                                     write_checksum     = True     ,
                                                     write_content_size = True     ) 
        self.__decompressor = zst.ZstdDecompressor ( )
        
        ## initialize the base class 
        CompressShelf.__init__ ( self                    ,
                                 dbname                  ,
                                 mode         = mode     ,
                                 compress     = compress ,
                                 compresstype = 'zstd'   , **kwargs ) 
        
        conf = { 'threads' : threads } 
        self.kwargs.update ( conf )
    
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
    
    # ==========================================================================
    ## compress (ZST)  the item  using compressor 
    def compress_item ( self , value ) :
        """ Compress (ZST) the item using compressor 
        """
        return self.compressor.compress (  self.pickle ( value ) )
    
    # =========================================================================
    ## uncompres (ZST) the item using decompressor 
    def uncompress_item ( self , value ) :
        """ Uncompress (ZST) the item using decompressor 
        """        
        return self.unpickle ( self.decompressor.decompress ( value ) ) 

# =============================================================================
## helper function to access ZstShelve data base
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2023-02-17
def open ( dbname , mode        = 'c' , **kwargs ) : 
    """ Open a persistent dictionary for reading and writing.
    
    The filename parameter is the base filename for the underlying
    database.  As a side-effect, an extension may be added to the
    filename and more than one file may be created.  The optional flag
    parameter has the same interpretation as the flag parameter of
    anydbm.open(). The optional protocol parameter specifies the
    version of the pickle protocol (0, 1, or 2).
        
    See the module's __doc__ string for an overview of the interface.
    """        
    return ZstShelf ( dbname , mode = mode  , **kwargs )
    
# =============================================================================
## @class TmpZstShelf
#  TEMPORARY zst-version of `shelve'-database
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2023-02-17
class TmpZstShelf(ZstShelf,TmpDB):
    """ TEMPORARY `ZST'-version of `shelve'-database     
    """    
    def __init__( self                              ,
                  protocol    = HIGHEST_PROTOCOL    , 
                  remove      = True                ,
                  keep        = False               , **kwargs ) :
        
        ## initialize the base: generate the name 
        TmpDB.__init__ ( self , suffix = '.zstdb' , remove = remove , keep = keep ) 
        
        ## open DB 
        ZstShelf.__init__ ( self                        ,  
                            self.tmp_name               ,
                            mode        = 'c'           ,
                            protocol    = protocol      ,
                            writeback   = False         , **kwargs  ) 
        
        conf = { 'remove' : remove , 'keep' : keep }
        self.kwargs.update ( conf )
        
    # =============================================================================
    ## close and delete the file 
    def close ( self )  :
        ## close the shelve file
        ZstShelf.close ( self )
        ## delete the file
        TmpDB   .clean ( self ) 

# =============================================================================
## helper function to open TEMPORARY ZstShelve data base#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2023-02-17
def tmpdb ( **kwargs ) :
    """ Open a TEMPORARY persistent dictionary for reading and writing.
        
    The optional protocol parameter specifies the
    version of the pickle protocol (0, 1, or 2).
    
    See the module's __doc__ string for an overview of the interface.
    """
    return TmpZstShelf ( **kwargs ) 

# ==========================================================================
if  not zst : __all__ = ()

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
