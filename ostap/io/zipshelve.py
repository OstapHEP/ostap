#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
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
#
# Create new DB:
#
# @code
#
# >>> import zipshelve  ## import the ZipShelve module 
# >>> db = zipShelve.open ('a_db', 'n')    ## create new DB
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
# >>> import zipshelve  ## import the ZipShelve module 
# >>> db = zipShelve.open ('a_db' , 'r' )    ## access existing dbase in read-only mode
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
# >>> import ZipShelve  ## import the ZipShelve module 
# >>> db = zipShelve.open ('a_db' )    ## access existing dbase in update mode
# ...
# >>> for key in db : print key
# ...
# >>> abcd = db['some_key']
#
# @endcode 
#
# @attention: In case DB-name has extention "gz", the whole data base
#             will be gzipped. 
#
#
# 
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2010-04-30
# 
#                    $Revision$
#  Last Modification $Date$
#                 by $Author$
# =============================================================================
""" This is zip-version of shelve database.

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

 >>> import zipshelve as DBASE  ## import the ZipShelve module 
 >>> db = DBASE.open ('a_db', 'n')    ## create new DB
 ...
 >>> abcde = ...
 >>> db['some_key'] =  abcde              ## add information to DB
 ...
 >>> db.close()

 Access to DB in read-only mode :

 >>> import zipshelve as DBASE  ## import the ZipShelve module 
 >>> db = DBASE.open ('a_db' , 'r' )    ## access existing dbase in read-only mode
 ...
 >>> for key in db : print key
 ...
 >>> abcd = db['some_key']

 Access existing DB in update mode :

 >>> import zipshelve as DBASE ## import the ZipShelve module 
 >>> db = DBASE.open ('a_db' )    ## access existing dbase in update mode
 ...
 >>> for key in db : print key
 ...
 >>> abcd = db['some_key']
 
 In case DB-name has extension 'gz', the whole data base will be gzipped

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
try:
    from cPickle   import Pickler, Unpickler, HIGHEST_PROTOCOL
except ImportError:
    from  pickle   import Pickler, Unpickler, HIGHEST_PROTOCOL 
# =============================================================================
try:
    from cStringIO import StringIO
except ImportError:
    from  StringIO import StringIO
# ==============================================================================
import os
import zlib        ## use zlib to compress DB-content 
import shelve      ## 
import shutil
# =============================================================================
_modes_ = {
    # =========================================================================
    # 'r'	Open existing database for reading only
    # 'w'	Open existing database for reading and writing
    # 'c'	Open database for reading and writing, creating it if it doesn’t exist
    # 'n'	Always create a new, empty database, open for reading and writing
    # =========================================================================
    'n'        : 'n' ,
    'c'        : 'c' ,
    'r'        : 'r' ,
    'u'        : 'w' ,
    'w'        : 'w' ,
    'a'        : 'w' ,
    ##
    '+'        : 'w' ,        
    'w+'       : 'w' ,
    'rw'       : 'w' ,
    'new'      : 'n' ,
    'create'   : 'c' ,
    'recreate' : 'n' ,
    'read'     : 'r' ,        
    'write'    : 'w' ,        
    'update'   : 'w' ,        
    'append'   : 'w' ,        
    }
_dbases = []
# =============================================================================
## Zipped-version of ``shelve''-database
#    Modes: 
#    - 'r' Open existing database for reading only
#    - 'w' Open existing database for reading and writing
#    - 'c' Open database for reading and writing, creating it if it doesn’t exist (default)
#    - 'n' Always create a new, empty database, open for reading and writing
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
class ZipShelf(shelve.Shelf):
    """ Zipped-version of ``shelve''-database
    Modes: 
    - 'r'  Open existing database for reading only
    - 'w'  Open existing database for reading and writing
    - 'c'  Open database for reading and writing, creating it if it doesn’t exist
    - 'n'  Always create a new, empty database, open for reading and writing
    Modes: %s 
    # =========================================================================
    """ % _modes_
    def __init__(
        self                                   ,
        filename                               ,
        mode      = 'c'                        , 
        protocol  = HIGHEST_PROTOCOL           , 
        compress  = zlib.Z_BEST_COMPRESSION    ,
        writeback = False                      ,
        silent    = False                      ) :

        ## the mode 
        mode = _modes_.get( mode.lower() , '' )
        if not mode :
            logger.warning("Unknown opening mode '%s', replace with 'c'")
            mode = 'c'
 
        #
        ## expand the actual file name 
        filename  = os.path.expandvars ( filename )
        filename  = os.path.expanduser ( filename )
        filename  = os.path.expandvars ( filename )
        filename  = os.path.expandvars ( filename )
        
        self.__gzip          = False 
        self.__filename      = filename
        self.__remove        = False
        self.__silent        = silent
        self.__opened        = False 

        if not self.__silent :
            logger.info ( 'Open DB: %s' % filename ) 
        
        if filename.endswith( '.gz' ) :
            
            if os.path.exists ( filename ) and 'r' == mode :
                ## gunzip into temporary location
                filename_ = self._gunzip ( filename ) 
                if not os.path.exists ( filename_ ) :
                    raise TypeError ( "Unable to gunzip properly: %s" % filename )
                if not self.__silent : 
                    size1 = os.path.getsize ( filename  ) 
                    size2 = os.path.getsize ( filename_ )
                    logger.info("GZIP uncompression %s: %.1f%%" %  ( filename , (size2*100.0)/size1 ) )
                filename        = filename_ 
                self.__filename = filename_
                self.__remove   = True
            elif os.path.exists ( filename ) and 'r' != mode :
                ## unzip in place
                filename_     = filename[:-3]
                # remove existing file (if needed) 
                if os.path.exists ( filename_ ) : os.remove ( filename_ )
                size1 = os.path.getsize ( filename  ) 
                # gunzip in place 
                self.__in_place_gunzip  ( filename ) 
                ##
                if not os.path.exists ( filename_ ) :
                    raise TypeError ( "Unable to gunzip properly: %s" % filename )
                if not self.__silent : 
                    size2 = os.path.getsize ( filename_ )
                    logger.info("GZIP uncompression %s: %.1f%%" %  ( filename , (size2*100.0)/size1 ) ) 
                filename        = filename_ 
                self.__gzip     = True 
                self.__filename = filename_
                self.__remove   = False
            else : 
                ## 
                filename        = filename[:-3]
                self.__gzip     = True 
                self.__filename = filename

        import anydbm
        shelve.Shelf.__init__ (
            self                                   ,
            anydbm.open ( self.__filename , mode ) ,
            protocol                               ,
            writeback                              )
        
        self.compresslevel = compress
        self.__opened      = True

        ## keep in the list of known/opened databases 
        #_dbases.append ( self )

    def filename ( self ) : return self.__filename
    def opened   ( self ) : return self.__opened

    ## valid, opened DB 
    def __nonzero__ ( self ) :
        return self.opened() and not isinstance ( self.dict , shelve._ClosedDict ) and not self.dict is None 
    
    ## destructor 
    def __del__ ( self ) :
        """ Destructor 
        """
        ## close if opened 
        if self.opened () : self.close()  

    ## iterator over good keys 
    def ikeys ( self , pattern = '' ) :
        """Iterator over avilable keys (patterns included).
        Pattern matching is performed accoriding to
        fnmatch/glob/shell rules [it is not regex!] 
        
        >>> db = ...
        >>> for k in db.ikeys('*MC*') : print k 
        
        """
        keys_ = self.keys()
        keys_.sort()
        if not pattern :
            good = lambda s,p : True
        else :
            import fnmatch
            good = lambda s,p : fnmatch.fnmatchcase ( k , p )
        
        for k in keys_ :
            if good ( k , pattern ) : yield k

    ## list the avilable keys 
    def ls    ( self , pattern = '' ) :
        """List the avilable keys (patterns included).
        Pattern matching is performed accoriding to
        fnmatch/glob/shell rules [it is not regex!] 

        >>> db = ...
        >>> db.ls() ## all keys
        >>> db.ls ('*MC*')        
        
        """
        for k in self.ikeys( pattern ): print k
        
    ## close and gzip (if needed)
    def close ( self ) :
        """ Close the file (and gzip it if required) 
        """
        if not self.opened() : return 
        ##
        shelve.Shelf.close ( self )
        self.__opened = False  
        ##
        if self.__remove and os.path.exists ( self.__filename ) :
            if not self.__silent :
                logger.info( 'REMOVE: ', self.__filename )
            os.remove ( self.__filename )
        ##
        if self.__gzip and os.path.exists ( self.__filename ) :
            # get the initial size 
            size1 = os.path.getsize ( self.__filename )
            # gzip the file
            self.__in_place_gzip ( self.__filename ) 
            #
            if not os.path.exists ( self.__filename + '.gz' ) :
                logger.warning( 'Unable to compress the file %s ' % self.__filename  )
            size2 = os.path.getsize( self.__filename + '.gz' )
            if not self.__silent : 
                logger.info( 'GZIP compression %s: %.1f%%' % ( self.__filename, (size2*100.0)/size1 ) ) 

        ## remove from list of known databases 
        if self in _dbases :
            _dbases.remove ( self )
            
    ## gzip the file (``in-place'') 
    def __in_place_gzip   ( self , filein ) :
        """Gzip the file ``in-place''
        
        It is better to use here ``os.system'' or ``popen''-family,
        but it does not work properly for multiprocessing environemnt        
        """
        if os.path.exists ( filein + '.gz' ) : os.remove ( filein + '.gz' )
        #
        # gzip the file 
        fileout = self._gzip ( filein )
        # 
        if os.path.exists ( fileout ) :
            # rename the temporary file 
            shutil.move ( fileout , filein + '.gz' )
            #
            import time
            time.sleep( 3 )
            #
            # remove the original
            os.remove   ( filein ) 

    ## gunzip the file (``in-place'') 
    def __in_place_gunzip   ( self , filein ) :
        """Gunzip the file ``in-place''
        
        It is better to use here ``os.system'' or ``popen''-family,
        but unfortunately it does not work properly for multithreaded environemnt        
        """
        #
        filename = filein[:-3]            
        if os.path.exists ( filename ) : os.remove ( filename )
        #
        # gunzip the file 
        fileout = self._gunzip ( filein )
        #        
        if os.path.exists ( fileout ) :
            # rename the temporary file 
            shutil.move ( fileout , filename )
            #
            import time
            time.sleep( 3 )
            #
            # remove the original
            os.remove   ( filein ) 

    ## gzip the file into temporary location, keep original
    def _gzip   ( self , filein ) :
        """ Gzip the file into temporary location, keep original
        """
        if not os.path.exists  ( filein  ) :
            raise NameError ( "GZIP: non existing file: " + filein )
        #
        fin  =      file ( filein  , 'r' )
        #
        import tempfile 
        fd,fileout = tempfile.mkstemp ( prefix = 'tmp_' , suffix = '_zdb.gz' )
        #
        import gzip 
        fout = gzip.open ( fileout , 'w' )
        #
        try : 
            for all in fin : fout.write ( all )
        finally:
            fout.close()
            fin .close()   
            import time
            time.sleep( 3 ) 
        return fileout
        
    ## gzip the file into temporary location, keep original
    def _gunzip ( self , filein ) :
        """Gunzip the file into temporary location, keep original
        """
        if not os.path.exists  ( filein  ) :
            raise NameError ( "GUNZIP: non existing file: " + filein )
        #
        import gzip 
        fin  = gzip.open ( filein  , 'r' )
        #
        import tempfile 
        fd,fileout = tempfile.mkstemp ( prefix = 'tmp_' , suffix = '_zdb' )
        fout = file ( fileout , 'w' )
        #
        try : 
            for all in fin : fout.write ( all )
        finally: 
            fout.close()
            fin .close()
            import time
            time.sleep( 3 ) 
        return fileout

    #
    ## some context manager functionality
    # 
    def __enter__ ( self      ) : return self 
    def __exit__  ( self , *_ ) : self.close ()

    ## print .
    def __repr__ ( self ) :
        
        kls = self.__class__
        if kls.__module__.startswith('__') and kls.__module__.endswith('__') :
            kname = kls.__name__
        else :
            kname = '%s.%s' % (  kls.__module__ , kls.__name__ ) 
        
        if   self and len(self) :
            return "%s('%s'): %d object(s)" % ( kname , self.filename(), len(self) ) 
        elif self :
            return "%s('%s'): empty"        % ( kname , self.filename() ) 
        return "Invalid/Closed %s('%s')"    % ( kname , self.filename() )
    
    __str__ = __repr__
    

# =============================================================================
## ``get-and-uncompress-item'' from dbase 
def _zip_getitem (self, key):
    """ ``get-and-uncompress-item'' from dbase 
    """
    try:
        value = self.cache[key]
    except KeyError:
        f = StringIO(zlib.decompress(self.dict[key]))
        value = Unpickler(f).load()
        if self.writeback:
            self.cache[key] = value
    return value

# =============================================================================
## ``set-and-compress-item'' to dbase 
def _zip_setitem ( self , key , value ) :
    """``set-and-compress-item'' to dbase 
    """
    if self.writeback:
        self.cache[key] = value
    f = StringIO()
    p = Pickler(f, self._protocol)
    p.dump(value)
    self.dict[key] = zlib.compress( f.getvalue(), self.compresslevel)

ZipShelf.__getitem__ = _zip_getitem
ZipShelf.__setitem__ = _zip_setitem


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
def _db_rrshift_ ( db , object ) :
    """Add an object into data base
    
    dbase  = ...
    object = ...
    dbase.ls() 
    object >> dbase ## add object into dbase 
    dbase.ls() 
    """
    if   hasattr ( object , 'GetName' ) : name = object.GetName()
    elif hasattr ( object , 'name'    ) : name = object.name   ()
    else : name = object.__class__.__name__
    #
    db [ name] = object 
    
ZipShelf.__rrshift__ = _db_rrshift_

# =============================================================================
## helper finction to access ZipShelve data base#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
def open ( filename                                   ,
           mode          = 'c'                        ,
           protocol      = HIGHEST_PROTOCOL           ,
           compresslevel = zlib.Z_BEST_COMPRESSION    , 
           writeback     = False                      ,
           silent        = True                       ) : 
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
                      silent        )



# =============================================================================
## TEMPORARY Zipped-version of ``shelve''-database
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2015-10-31
class TmpZipShelf(ZipShelf):
    """
    TEMPORARY Zipped-version of ``shelve''-database     
    """    
    def __init__(
        self                                   ,
        protocol  = HIGHEST_PROTOCOL           , 
        compress  = zlib.Z_BEST_COMPRESSION    ,
        silent    = False                      ) :

        ## create temporary file name 
        import tempfile
        filename = tempfile.mktemp  ( suffix = '.zdb' )
        
        ZipShelf.__init__ ( self     ,  
                            filename ,
                            'n'      ,
                            protocol ,
                            compress , 
                            False    , ## writeback 
                            silent   ) 
        
    ## close and delete the file 
    def close ( self )  :
        ## close the shelve file
        fname = self.filename() 
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
def tmpdb ( protocol      = HIGHEST_PROTOCOL           ,
            compresslevel = zlib.Z_BEST_COMPRESSION    , 
            silent        = True                       ) : 
    """Open a TEMPORARY persistent dictionary for reading and writing.
    
    The optional protocol parameter specifies the
    version of the pickle protocol (0, 1, or 2).
    
    See the module's __doc__ string for an overview of the interface.
    """
    return TmpZipShelf ( protocol      ,
                         compresslevel ,
                         silent        )
    

# ============================================================================
logger.debug ( "Simple generic (c)Pickle-based ``zipped''-database")
# =============================================================================
## a bit more decorations for shelve  (optional)
import ostap.io.shelve_ext


# =============================================================================
# some gymnastic to close all DBs 
# =============================================================================
import atexit
def _close_dbs_ () :
    while _dbases :
        db = _dbases.pop()
        if db :
            logger.info ('Close ZipShelve database "%s"' % db.filename() ) 
            db.close()
        del db
        
atexit.register ( _close_dbs_ ) 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.logger.line import line 
    logger.info ( __file__  + '\n' + line  ) 
    logger.info ( 80*'*'   )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 
    
# =============================================================================
# The END 
# =============================================================================
