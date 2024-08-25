#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file lzshelve.py
# 
# This is `LZMA'-version of shelve database.
# 
# Keeping the same interface and functionlity as shelve data base,
# LzmaShelf allows much more compact file size through the on-flight
# compression of the content
#
# The actual code has been inspired by <c>zipshelve</c> ( see Google...)
#
# However is contains several new features:
# 
#  - Optionally it is possible to perform the compression
#    of the whole data base, that can be rather useful for data base
#    with very large amout of keys 
#
# The module has been developed and used with great success in
# ``Kali, framework for fine calibration of LHCb Electormagnetic Calorimeter''
#
# Create new DB:
#
# @code
#
# >>> import lzshelve  
# >>> db = lzshelve.open ('a_db', 'n')    ## create new DB
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
# >>> import lzshelve  
# >>> db = lzshelve.open ('a_db' , 'r' ) ## access existing dbase in read-only mode
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
# >>> import lzshelve 
# >>> db = bz2shelve.open ('a_db' )    ## access existing dbase in update mode
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
"""This is ``LZMA''-version of shelve database.

Keeping the same interface and functionlity as shelve data base,
LzmaShelf allows much more compact file size through the on-flight
compression of the content

The actual code has been inspired by zipshelve ( see Google...)

However is contains several new features:

 - Optionally it is possible to perform the compression
   of the whole data base, that can be rathe useful fo data base
   with large amout of keys 
   
The module has been developed and used with great success in
 ``Kali, framework for fine calibration of LHCb Electormagnetic Calorimeter''

 Create new DB:

 >>> import lzshelve as DBASE  
 >>> db = DBASE.open ('a_db', 'n')    ## create new DB
 ...
 >>> abcde = ...
 >>> db['some_key'] =  abcde          ## add information to DB
 ...
 >>> db.close()

 Access to DB in read-only mode :

 >>> import lzshelve as DBASE  
 >>> db = DBASE.open ('a_db' , 'r' ) ## access existing dbase in read-only mode
 ...
 >>> for key in db : print(key)
 ...
 >>> abcd = db['some_key']

 Access existing DB in update mode :

 >>> import lzshelve as DBASE
 >>> db = DBASE.open ('a_db' )    ## access existing dbase in update mode
 ...
 >>> for key in db : print(key)
 ...
 >>> abcd = db['some_key']
 
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
__version__ = "$Revision:$" 
# =============================================================================
__all__ = (
    'LzShelf'     , ## database
    'TMpLzShelf'  , ## database
    'open'        , ## open the database 
    'tmpdb'       , ## open the temporary database 
)
# =============================================================================
from sys import version_info as python_version 
# =============================================================================
import os, sys, shelve, shutil 
from   ostap.io.compress_shelve import CompressShelf, HIGHEST_PROTOCOL
from   ostap.io.dbase           import TmpDB 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.lzshelve' )
else                      : logger = getLogger ( __name__             )
# =============================================================================
logger.debug ( "Simple generic (c)Pickle-based ``LZMA''-database"    )
# =============================================================================
if ( 3 , 3 ) <= python_version :
    try :
        import lzma
    except ImportError :
        logger.error ( "Cannot import 'lzma', lzshelve is disabled" )
        lzma = None
else :
    logger.error ( 'lzshelve is disabled for python %s' % str ( python_version ) )
    lzma = None

# =============================================================================
## @class LzShelf
#  ``LZMA''-version of ``shelve''-database
#    Modes: 
#    - 'r' Open existing database for reading only
#    - 'w' Open existing database for reading and writing
#    - 'c' Open database for reading and writing, creating if it does not  exist (default)
#    - 'n' Always create a new, empty database, open for reading and writing
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
class LzShelf(CompressShelf):
    """ LZMA-version of ``shelve''-database
    Modes: 
    - 'r'  Open existing database for reading only
    - 'w'  Open existing database for reading and writing
    - 'c'  Open database for reading and writing, creating if it does not exist
    - 'n'  Always create a new, empty database, open for reading and writing
    """ 
    ## 
    def __init__( self                              ,
                  dbname                            ,
                  mode        = 'c'                 ,
                  compress    = lzma.PRESET_DEFAULT , **kwargs ) :

        assert lzma, "`lzma` module is not available!"
        
        ## initialize the base class 
        CompressShelf.__init__ ( self                    ,
                                 dbname                  ,
                                 mode         = mode     ,
                                 compress     = compress ,
                                 compresstype = 'lzma'   , **kwargs ) 

        self.taropts = 'x:xz'
            
    # ==========================================================================
    ## compress (LZMA)  the item  using <code>lzma.compress</code>
    def compress_item ( self , value ) :
        """ Compress (LZMA) the item using ``bz2.compress''
        - see lzma.compress
        """
        return lzma.compress ( self.pickle ( value ) , preset = self.compresslevel )
    
    # =========================================================================
    ## uncompres (LZMA) the item using <code>lzma.decompress</code>
    def uncompress_item ( self , value ) :
        """ Uncompress (LZMA) the item using ``lzma.decompress''
        -  see lzma.decompress
        """        
        return self.unpickle ( lzma.decompress ( value ) )
    
# =============================================================================
## helper function to access LzShelve data base
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
def open ( dbname , mode   = 'c' , **kwargs ) :
    """ Open a persistent dictionary for reading and writing.
        
    The filename parameter is the base filename for the underlying
    database.  As a side-effect, an extension may be added to the
    filename and more than one file may be created.  The optional flag
    parameter has the same interpretation as the flag parameter of
    anydbm.open(). The optional protocol parameter specifies the
    version of the pickle protocol (0, 1, or 2).
    
    See the module's __doc__ string for an overview of the interface.
    """    
    return LzShelf ( dbname , mode = mode , **kwargs )
    
# =============================================================================
## @class TmpLzShelf
#  TEMPORARY lzma-version of ``shelve''-database
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2015-10-31
class TmpLzShelf(LzShelf,TmpDB):
    """ TEMPORARY ``LZMA''-version of ``shelve''-database     
    """    
    def __init__( self                              ,
                  protocol    = HIGHEST_PROTOCOL    , 
                  remove      = True                ,
                  keep        = False               , **kwargs ) :
    
        ## initialize the base: generate the name 
        TmpDB.__init__ ( self , suffix = '.lzdb' , remove = remove , keep = keep ) 
        
        ## open DB 
        LzShelf.__init__ ( self                      ,  
                           self.tmp_name             ,
                           mode        = 'c'         ,
                           protocol    = protocol    ,
                           writeback   = False       , **kwargs ) 
        
        conf = { 'remove' : remove , 'keep' : keep }
        self.kwargs.update ( conf )
        
    ## close and delete the file 
    def close ( self )  :
        ## close the shelve file
        LzShelf.close ( self )
        ## delete the file
        TmpDB  .clean ( self ) 
            
# =============================================================================
## helper function to open TEMPORARY ZipShelve data base#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
def tmpdb ( **kwargs ) :
    """ Open a TEMPORARY persistent dictionary for reading and writing.
        
    The optional protocol parameter specifies the
    version of the pickle protocol (0, 1, or 2).
        
    See the module's __doc__ string for an overview of the interface.
    """
    return TmpLzShelf ( **kwargs )

# ==========================================================================
if not lzma : __all__ = ()

# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
