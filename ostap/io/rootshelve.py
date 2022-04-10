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
й >>> db = DBASE.open ('a_db' )    ## access existing dbase in update mode
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
import ROOT, shelve, zlib, os 
import ostap.io.root_file
from   sys import version_info as python_version 
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
from   ostap.io.dbase import TmpDB 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.rootshelve' )
else                      : logger = getLogger ( __name__              )
logger.debug ( "Simple generic ROOT-based shelve-like-database" )
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
            rfile = ROOT.TFile.Open ( filename   , open_mode ( mode ) , *args )
            if not rfile or not rfile.IsOpen() :
                raise IOError("Can't open ``%s'' with mode ``%s''" % ( filename , mode ) )
            shelve.Shelf.__init__ ( self , rfile , writeback )
            
        self.nominal_dbname = filename
        
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
        new_db = RootOnlyShelf ( new_name                         ,
                                 mode        =  'c'               ,
                                 writeback   = self.writeback     )
        ## copy the content
        if keys :
            for key in self.keys () :
                if key in keys      : new_db [ key ] = self [ key ]
        else : 
            for key in self.keys () : new_db [ key ] = self [ key ]
         
        new_db.sync ()  
        return new_db 

    # =========================================================================
    ##  Iterator over avilable keys (patterns included).
    #   Pattern matching is performed accoriding to
    #   fnmatch/glob/shell rules (default) or regex 
    #   @code  
    #   db = ...
    #   for k in db.ikeys('*MC*') : print(k)
    #   @endcode  
    def ikeys ( self , pattern = '' , regex = False ) :
        """Iterator over avilable keys (patterns included).
        Pattern matching is performed according to
        fnmatch/glob/shell rules (default) or regex 
        
        >>> db = ...
        >>> for k in db.ikeys('*MC*') : print(k)
        
        """
        keys_ = self.keys()
        
        if not pattern :
            good = lambda k : True
        elif regex : 
            import re
            re_cmp = re.compile ( pattern ) 
            good = lambda k  : re_cmp.match ( k )
        else :
            import fnmatch
            good = lambda s : fnmatch.fnmatchcase ( k , pattern  )
        
        keys_ = self.keys()
        for k in sorted ( keys_ ) :
            if good ( k ) : yield k

    @property
    def filename    ( self       ) :
        """``filename'' : the file name for root-database"""
        return self.__filename
    
    def __enter__   ( self       ) : return self 
    def __exit__    ( self , *_  ) : self.close ()


    # =============================================================================
    ## get the object from the ROOT file 
    def get ( self , key , defval = None ) :
        """Get the object from the ROOT file
        """
        return self.dict.get ( key , defval ) 

    # =============================================================================
    ## get item from ROOT-file
    #  @code
    #  obj = db['A/B/C/histo']
    #  @endcode 
    #  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
    #  @date   2015-07-31 
    def __getitem__ ( self , key ) :
        """Get the item from ROOT-file
        >>> obj = db['A/B/C/histo']
        """
        try:
            value = self.cache [ key ]
        except KeyError:
            value = self.dict [ key ] 
            if self.writeback:
                self.cache [ key ] = value
        return value
    
    # =============================================================================
    ## put item into ROOT-file 
    #  @code
    #  db['A/B/C/histo'] = obj
    #  @endcode 
    #  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
    #  @date   2015-07-31 
    def __setitem__ ( self , key , value ) :
        """ Put item into ROOT-file
        >>> db ['A/B/C/histo'] = obj
        """
        if self.writeback : self.cache [ key ] = value
        self.dict [ key ] = value 

    ## close the database 
    def close ( self ) :
        """Close the database
        """
        shelve.Shelf.close ( self )
        
    # =========================================================================
    ## list the available keys 
    def ls    ( self ) :
        """List the available keys as table .
        
        >>> db = ...
        >>> db.ls() ## all keys
        
        """
        return self.dict.ls_table( prefix = "# ")


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
        self.__sizes         = {}
        
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
        new_db = RootShelf ( new_name                         ,
                             mode        =  'c'               ,
                             protocol    = self.protocol      ,
                             compress    = self.compresslevel )
        
        ## copy the content
        if keys :
            for key in self.keys() :
                if key in keys     : new_db [ key ] = self [ key ]
        else : 
            for key in self.keys() : new_db [ key ] = self [ key ]
         
        new_db.sync ()  
        return new_db 

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
    def __getitem__ ( self , key ):
        """ Get object (unpickle if needed)  from dbase
        >>> obj = db['A/B/C']
        """

        try:    
            value = self.cache [ key ]
        except KeyError:

            ## value = self.dict [ key ]
            tkey , value = self.dict.get_key_object ( key )
            self.__sizes [ key ] = tkey.GetNbytes()
            
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
    def __setitem__ ( self , key , value ) :
        """ Add object (pickle if needed)  to dbase
        >>> db['A/B/C'] = obj
        """
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
            self.__sizes [ key ] = len ( z ) 
            ## (3) put it into  BLOB 
            from  ostap.core.core import  Ostap
            blob   = Ostap.BLOB            ( key      ) 
            status = Ostap.blob_from_bytes ( blob , z )
            value  = blob 
            del z , f , p 
        
        ## finally use ROOT 
        self.dict [ key ] = value

    ## close the database 
    def close ( self ) :
        """Close the database
        """
        RootOnlyShelf.close ( self )
        
    # =========================================================================
    ## list the available keys 
    def ls    ( self , pattern = '' , load = True ) :
        """List the available keys (patterns included).
        Pattern matching is performed accoriding to
        fnmatch/glob/shell rules [it is not regex!] 

        >>> db = ...
        >>> db.ls() ## all keys
        >>> db.ls ('*MC*')        
        
        """
        n  = os.path.basename ( self.filename )
        ap = os.path.abspath  ( self.filename ) 
        
        try :
            fs = os.path.getsize ( self.filename )
        except :
            fs = -1
            
        if    fs < 0            : size = "???"
        elif  fs < 1024         : size = str(fs)
        elif  fs < 1024  * 1024 :
            size = '%.2fkB' % ( float ( fs ) / 1024 )
        elif  fs < 1024  * 1024 * 1024 :
            size = '%.2fMB' % ( float ( fs ) / ( 1024 * 1024 ) )
        else :
            size = '%.2fGB' % ( float ( fs ) / ( 1024 * 1024 * 1024 ) )
            
                        
        keys = [] 
        for k in self.ikeys ( pattern ): keys.append ( k )
        keys.sort()
        if keys : mlen = max ( [ len(k) for k in keys] ) + 2 
        else    : mlen = 2 
        fmt = ' --> %%-%ds : %%s' % mlen

        table = [ ( 'Key' , 'type' , '   size   ') ] 
        for k in keys :
            size = '' 
            ss   =   self.__sizes.get ( k , -1 )
            if    ss < 0    : size = ''  
            elif  ss < 1024 : size = '%7d   ' % ss 
            elif  ss < 1024 * 1024 :
                size = '%7.2f kB' %  ( float ( ss ) / 1024 )
            elif  ss < 1024 * 1024 * 1024 :
                size = '%7.2f MB' %  ( float ( ss ) / ( 1024 * 1024 ) )
            else :
                size = '%7.2f GB' %  ( float ( ss ) / ( 1024 * 1024 * 1024 ) )
                
            ot    = type ( self [ k ] )
            otype = ot.__cppname__ if hasattr ( ot , '__cppname__' ) else ot.__name__ 
            row = '{:15}'.format ( k ) , '{:15}'.format ( otype ) , size 
            table.append ( row )

        import ostap.logger.table as T
        t      = self.__class__.__name__
        title  = '%s:%s' % ( t  , n )
        maxlen = 0
        for row in table :
            rowlen = 0 
            for i in row : rowlen += len ( i )
            maxlen = max ( maxlen, rowlen ) 
        if maxlen + 3 <= len ( title ) :
            title = '<.>' + title [ -maxlen : ] 
        table = T.table ( table , title = title , prefix = '# ' )
        ll    = getLogger ( n )
        line  = 'Database %s:%s #keys: %d size: %s' % ( t , ap , len ( self ) , size )
        ll.info (  '%s\n%s' %  ( line , table ) )


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
           writeback     = False ,
           root_only     = False , *args ) : 
    """
    Helper function to open RootShelve data base
    >>> import RootShelve as DBASE
    >>> db = DBASE.open ( 'mydb.root' , 'c' )
    """
    db_type = RootOnlyShelf if root_only else RootShelf
    return db_type ( filename  ,
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
class TmpRootShelf(RootShelf,TmpDB):
    """The actual class for TEMPORARY ROOT-based shelve-like data base
    it implement shelve-intergase with underlyinog ROOT-fiel storage
    - ROOT-object are stored directly in the ROOT-file,
    - other objects are pickled and stored via  ROOT.TObjString
    see RootShelf
    """
    def __init__( self,
                  protocol  = HIGHEST_PROTOCOL           ,
                  compress  = zlib.Z_DEFAULT_COMPRESSION ,
                  remove    = True                       , ## immediate remove 
                  keep      = False                      , ## keep it 
                  *args                                  ):
        
        ## initialize the base: generate the name 
        TmpDB.__init__ ( self , suffix = '.root' , remove = remove , keep = keep ) 
        
        RootShelf.__init__ ( self                 ,
                             self.tmp_name        ,
                             mode      = 'n'      ,
                             writeback = False    ,
                             protocol  = protocol ,
                             compress  = compress ,
                             args      = args     )
        
    ## close and delete the file 
    def close ( self )  :
        RootShelf.close ( self )
        TmpDB.clean     ( self ) 
            
# =============================================================================
## helper function to open RootShelve data base
#  @code
#  import RootShelve as DBASE
#  db = DBASE.open ( 'mydb.root' , 'c' )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
def tmpdb ( protocol  = HIGHEST_PROTOCOL           ,
            compress  = zlib.Z_DEFAULT_COMPRESSION ,
            remove    = True                       ,            ## immediate remove 
            keep      = False                      ,
            *args                                  ) : ## keep it 
    """ Helper function to open TEMPPORARY RootShelve data base
    >>> import RootShelve as DBASE
    >>> db = DBASE.tmpdb()
    """    
    return TmpRootShelf ( protocol = protocol , 
                          compress = compress ,
                          remove   = remove   ,
                          keep     = keep     , *args ) 



# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
