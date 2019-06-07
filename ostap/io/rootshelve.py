#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file rootshelve.py
# 
# This is shelve-like database with ROOT.TFile as internal storage 
#
# @see zipshelve
# @see sqliteshelve
#
#
# Create new DB:
#
# @code
#
# >>> import rootshelve as DBASE   ## import the RootShelve module 
# >>> db = DBASE.open ('a_db', 'n')    ## create new DB
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
# >>> import rootshelve  as DBASE ## import the ZipShelve module 
# >>> db = DBASE.open ('a_db' , 'r' )    ## access existing dbase in read-only mode
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
# >>> import rootshelve as DBASE  ## import the ZipShelve module 
# >>> db = DBASE.open ('a_db' )    ## access existing dbase in update mode
# ...
# >>> for key in db : print(key)
# ...
# >>> abcd = db['some_key']
#
# @endcode 
#
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2015-07-31
# 
# =============================================================================
""" This is ROOT-based version of shelve database.

 Create new DB:

 >>> import rootshelve as DBASE ## import the ZipShelve module 
 >>> db = DBASE.open ('a_db', 'n')    ## create new DB
 ...
 >>> abcde = ...
 >>> db['some_key'] =  abcde              ## add information to DB
 ...
 >>> db.close()

 Access to DB in read-only mode :

 >>> import rootshelve as DBASE  ## import the ZipShelve module 
 >>> db = DBASE.open ('a_db' , 'r' )    ## access existing dbase in read-only mode
 ...
 >>> for key in db : print(key)
 ...
 >>> abcd = db['some_key']

 Access existing DB in update mode :

 >>> import rootshelve as DBASE   ## import the RootShelve module 
 >>> db = DBASE.open ('a_db' )    ## access existing dbase in update mode
 ...
 >>> for key in db : print(key)
 ...
 >>> abcd = db['some_key']
 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2015-07-31"
__version__ = "$Revision$" 
# =============================================================================
__all__ = (
    'RootShelf'     , ## The DB-itself
    'RootOnlyShelf' , ## "data base" for ROOT-only objects
    'open'          , ## helper function to hide the actual DB
    'tmpdb'         , ## helper function to create TEMPORARY  RootShelve database 
    )
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.rootshelve' )
else                      : logger = getLogger ( __name__              )
# =============================================================================
import ROOT, shelve, zlib
import ostap.io.root_file  
from   sys import version_info as python_version 
logger.debug ( "Simple generic ROOT-based shelve-like-database" )
# =============================================================================
try : 
    from cPickle import Pickler, Unpickler, HIGHEST_PROTOCOL
except ImportError : 
    from  pickle import Pickler, Unpickler, HIGHEST_PROTOCOL
# =============================================================================    
try :
    from io     import             BytesIO 
except ImportError : 
    from shelve import StringIO as BytesIO
# =============================================================================
PROTOCOL = 2
# =============================================================================
## @class RootOnlyShelf
#  Plain vanilla DBASE for ROOT-object (only)
#  essentially it is nothing more than just shelve-like interface for ROOT-files
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2015-07-31
#  @attention It CRUCIALLY depends on the proper TFile-decorations
#                          from ostap.io.root_file module
#  @code
#  db = RootOnlyShelf('mydb.root','c')
#  h1 = ...
#  db ['histogram'] = h1
#  db.ls()
#  @endcode 
#  @see Ostap.TFileDeco
class RootOnlyShelf(shelve.Shelf):
    """Plain vanilla DBASE for ROOT-object (only)
    Essentially it is nothing more than just shelve-like
    interface for ROOT-files
    Attention: It CRUCIALLY depends on the proper
    TFile-decorations from ostap.io.root_file module
    
    >>> db = RooOnlyShelf('mydb.root','c')
    >>> h1 = ...
    >>> db ['histogram'] = h1
    >>> db.ls()
    """ 
    ## constructors 
    #  @attention it depends on proper TFile-decorations in ostap.io.root_file module
    def __init__( self                  ,
                  filename              ,
                  mode                  ,
                  writeback   = False   ,
                  args        = ()      ) :
        """ Create Root-only database 
        >>> db = RooOnlyShelf('mydb.root','c')
        >>> h1 = ...
        """
        self.__filename = filename 
        from ostap.io.root_file import ROOTCWD, open_mode  
        with ROOTCWD() : ## NB: preserve current directory in ROOT!
            rfile = ROOT.TFile.Open ( filename   , open_mode ( mode ) , *args  )
            shelve.Shelf.__init__ ( self , rfile , writeback )

    @property
    def filename    ( self       ) :
        """``filename'' : the file name for root-database"""
        return self.__filename
    
    def __enter__   ( self       ) : return self 
    def __exit__    ( self , *_  ) : self.close ()


# =============================================================================
## need to disable endcode/decode for the keys 
if python_version.major > 2 :
    
    def _ros_iter_     ( self ):
        for k in self.dict.keys() : yield k
    def _ros_contains_ ( self , key ):
        return key in self.dict
    def _ros_get_      ( self , key , default = None ) :
        if key in self.dict       : return self[key]
        return default
    def _ros_delitem_  ( self , key ) :
        del self.dict[key]
        try:
            del self.cache[key]
        except KeyError:
            pass

    RootOnlyShelf.__iter__     = _ros_iter_ 
    RootOnlyShelf.__contains__ = _ros_contains_
    RootOnlyShelf.__ros_get__  = _ros_get_
    RootOnlyShelf.__delitem__  = _ros_delitem_
    
# =============================================================================
## get item from ROOT-file
#  @code
#  obj = db['A/B/C/histo']
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2015-07-31 
def _root_getitem_ ( self , key ) :
    """Get the item from ROOT-file
    >>> obj = db['A/B/C/histo']
    """
    try:
        value = self.cache[key]
    except KeyError:
        value = self.dict[ key ] 
        if self.writeback:
            self.cache[key] = value
    return value
    
# =============================================================================
## put item into ROOT-file 
#  @code
#  db['A/B/C/histo'] = obj
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2015-07-31 
def _root_setitem_ ( self , key , value ) :
    """ Put item into ROOT-file
    >>> db['A/B/C/histo'] = obj
    """
    if self.writeback:
        self.cache[key] = value
    self.dict[ key ] = value 

RootOnlyShelf.__getitem__ = _root_getitem_
RootOnlyShelf.__setitem__ = _root_setitem_

# =============================================================================
## add an object into data base
#  @code
#  dbase  = ...
#  object = ...
#  dbase.ls() 
#  object >> dbase ## add object into dbase 
#  dbase.ls() 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-06-04
def _db_rrshift_ ( dbase , obj ) :
    """Add an object into data base
    
    dbase  = ...
    object = ...
    dbase.ls() 
    object >> dbase ## add object into dbase 
    dbase.ls() 
    """
    if   hasattr ( obj , 'GetName' ) : name = obj.GetName()
    elif hasattr ( obj , 'name'    ) : name = obj.name   ()
    else : name = obj.__class__.__name__
    #
    dbase [ name ] = obj
    
RootOnlyShelf.__rrshift__ = _db_rrshift_

# =============================================================================
## @class RootShelf
#  The actual class for ROOT-based shelve-like data base
#  it implement shelve-interface with underlying ROOT-file as storage
#  - ROOT-objects are stored directly in the ROOT-file,
#  - other objects are pickled and stored via ROOT.TObjString
#  @code
#  db = RootShelf( 'mydb.root' , 'c' )
#  db['histo'] = h1
#  db['tuple'] = ('a',1,h1) 
#  @endcode
#  @see RootOnlyShelf 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2015-07-31 
class RootShelf(RootOnlyShelf):
    """ The actual class for ROOT-based shelve-like data base
    it implement shelve-interface with underlyinog ROOT-fiel storage
    - ROOT-object are store ddirectly in the ROOT-file,
    - other objects are pickled and stored in ROOT.TObjString
    
    >>> db = RootShelf( 'mydb.root' , 'c' )
    >>> db['histo'] = h1
    >>> db['tuple'] = ('a',1,h1)
    """
    def __init__( self                                ,
                  filename                            ,
                  mode                                ,
                  writeback = False                   ,
                  protocol  = PROTOCOL                , ## pickling protocol
                  compress  = zlib.Z_BEST_COMPRESSION , ## compression level 
                  args      = ()       ):
        RootOnlyShelf.__init__ ( self , filename , mode , writeback , args = args )
        self.__protocol      = protocol
        self.__compresslevel = compress 
    @property
    def protocol ( self ) :
        """``protocol'' : pickle protocol"""
        return self.__protocol
    @property
    def compresslevel ( self ) :
        """``compresslevel'' : zlib compression level
        """
        return self.__compresslevel
        
# =============================================================================
##  get object (unpickle if needed)  from dbase
#   @code
#   obj = db['A/B/C']
#   @endcode
#   @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#   @date   2015-07-31 
def _pickled_getitem_ ( self , key ):
    """ Get object (unpickle if needed)  from dbase
    >>> obj = db['A/B/C']
    """
    ##
    try:    
        value = self.cache[key]
    except KeyError:
        value = self.dict[key]
        ## blob ?
        from  ostap.core.core import  Ostap
        if isinstance ( value , Ostap.BLOB ) :
            ## unpack it!
            z     = Ostap.blob_to_bytes ( value )
            u     = zlib.decompress ( z )
            ## unpickle it! 
            f     = BytesIO ( u )
            value = Unpickler(f).load()
            del z , u , f            
        if self.writeback:
            self.cache[key] = value
    return value

# =============================================================================
##  Add object (pickle if needed)  to dbase
#   @code
#   db['A/B/C'] = obj
#   @endcode
#   @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#   @date   2015-07-31 
def _pickled_setitem_ ( self , key , value ) :
    """ Add object (pickle if needed)  to dbase
    >>> db['A/B/C'] = obj
    """
    ##
    if self.writeback:
        self.cache [ key ] = value

    ## not TObject? pickle it and convert to Ostap.BLOB
    if not isinstance  ( value , ROOT.TObject ) :
        ## (1) pickle it 
        f = BytesIO    ( )
        p = Pickler    ( f , self.protocol )
        p.dump ( value )
        ## (2) zip it
        z      = zlib.compress ( f.getvalue() , self.compresslevel )
        ## (3) put it into  BLOB 
        from  ostap.core.core import  Ostap
        blob   = Ostap.BLOB            ( key      ) 
        status = Ostap.blob_from_bytes ( blob , z )
        value  = blob 
        del z , f , p 
        
    ## finally use ROOT 
    self.dict[key] = value
    
RootShelf.__getitem__ = _pickled_getitem_
RootShelf.__setitem__ = _pickled_setitem_

RootShelf.__rrshift__ = _db_rrshift_
# =============================================================================
## helper function to open RootShelve data base
#  @code
#  import RootShelve as DBASE
#  db = DBASE.open ( 'mydb.root' , 'c' )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
def open ( filename              ,
           mode          = 'c'   ,
           writeback     = False , *args ) : 
    """
    Helper function to open RootShelve data base
    >>> import RootShelve as DBASE
    >>> db = DBASE.open ( 'mydb.root' , 'c' )
    """    
    return RootShelf ( filename  ,
                       mode      ,
                       writeback , * args )

# =============================================================================
## @class TmpRootShelf
#  TEMPORARY The actual class for ROOT-based shelve-like data base
#  it implements shelve-interface with underlying ROOT-file as a storage
#  - ROOT-objects are stored directly in the ROOT-file,
#  - other objects are pickled and stored in ROOT.TObjString
#  @code
#  db = TmmRootShelf()
#  db['histo'] = h1
#  db['tuple'] = ('a',1,h1) 
#  @endcode
#  @see RootShelf 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2015-07-31 
class TmpRootShelf(RootShelf):
    """The actual class for TEMPORARY ROOT-based shelve-like data base
    it implement shelve-intergase with underlyinog ROOT-fiel storage
    - ROOT-object are stored directly in the ROOT-file,
    - other objects are pickled and stored via  ROOT.TObjString
    see RootShelf
    """
    def __init__( self, *args ):
        
        ## create temporary file name 
        import tempfile
        filename = tempfile.mktemp  ( suffix = '.root' )
        
        RootShelf.__init__ ( self                                   ,
                             filename                               ,
                             mode      = 'n'                        ,
                             writeback = False                      ,
                             protocol  = HIGHEST_PROTOCOL           ,
                             compress  = zlib.Z_DEFAULT_COMPRESSION ,
                             args      = args                       )
        
    ## close and delete the file 
    def close ( self )  :
        ## close the shelve file
        fname = self.filename 
        ## super(TmpRootShelf,self).close ()
        shelve.Shelf.close ( self )
        ## delete the file
        import os
        if os.path.exists ( fname ) :
            try :
                os.unlink ( fname )
            except :
                pass
            
# =============================================================================
## helper function to open RootShelve data base
#  @code
#  import RootShelve as DBASE
#  db = DBASE.open ( 'mydb.root' , 'c' )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
def tmpdb ( *args ) : 
    """ Helper function to open TEMPPORARY RootShelve data base
    >>> import RootShelve as DBASE
    >>> db = DBASE.tmpdb()
    """    
    return TmpRootShelf ( *args )


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
