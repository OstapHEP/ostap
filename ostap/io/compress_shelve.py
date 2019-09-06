#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file compress_shelve.py
# 
# Abstract helper class for "compressed" version of of shelve database.
# 
# Keeping the same interface and functionlity as shelve data base,
# it allows much more compact file size through the on-flight
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
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2010-04-30
# 
# =============================================================================
""" Abstract helper class for ``compressed'' version of of shelve database.

Keeping the same interface and functionlity as shelve data base,
it allows much more compact file size through the on-flight
compression of the content

The actual code has been inspired by zipshelve ( see Google...)

However is contains several new features:

 - Optionally it is possible to perform the compression
   of the whole data base, that can be rathe useful fo data base
   with large amout of keys 
   
The module has been developed and used with great success in
 ``Kali, framework for fine calibration of LHCb Electormagnetic Calorimeter''

"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2010-04-30"
__version__ = "$Revision$" 
# =============================================================================
__all__ = (
    'CompressShelf' ,   ## The DB-itself
    )
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.compress_shelve' )
else                      : logger = getLogger ( __name__                   )
# =============================================================================
logger.debug ( "Absract class for  (c)Pickle-based ``compressed''-database")
# =============================================================================
## to be compatible between  Python2 and Python3 
PROTOCOL = 2
ENCODING = 'utf-8'
# ==============================================================================
import os, sys, abc, shelve, shutil , time 
from  sys import version_info as python_version
try                : import anydbm as dbm
except ImportError : import           dbm
# =============================================================================
_modes_ = {
    # =========================================================================
    # 'r'	Open existing database for reading only
    # 'w'	Open existing database for reading and writing
    # 'c'	Open database for reading and writing, creating it if it doesn't exist
    # 'n'	Always create a new, empty database, open for reading and writing
    # =========================================================================
    'n'        : 'c' ,
    'c'        : 'c' ,
    'r'        : 'r' ,
    'u'        : 'w' ,
    'w'        : 'w' ,
    'a'        : 'w' ,
    ##
    '+'        : 'w' ,        
    'w+'       : 'w' ,
    'rw'       : 'w' ,
    'new'      : 'c' ,
    'create'   : 'c' ,
    'recreate' : 'n' ,
    'read'     : 'r' ,        
    'write'    : 'w' ,        
    'update'   : 'w' ,        
    'append'   : 'w' ,        
    }
# =============================================================================
## @class CompressShelf
#  Abstract class for ``compressed'' version of ``shelve''-database.
#  It has four abstract methods:
#  - <code>compress_item</code>
#  - <code>uncompress_item</code>
#  - <code>compress_file</code>
#  - <code>uncompress_file</code>
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
class CompressShelf(shelve.Shelf,object):
    """ ``Compressed''-version of ``shelve''-database
    It has four abstract methods:
    - compress_item
    - uncompress_item
    - compress_file
    - uncompress_file
    """
    __metaclass__ = abc.ABCMeta
    extensions    = ()
    
    def __init__(
        self                   ,
        filename               ,
        mode        = 'c'      , 
        protocol    = PROTOCOL , 
        writeback   = False    ,
        silent      = False    ,
        keyencoding = 'utf-8'  ) :


        ## the mode 
        mode = _modes_.get( mode.lower() , '' )
        if not mode :
            logger.warning ( "Unknown opening mode '%s', replace with 'c'")
            mode = 'c'

        ## expand the actual file name 
        filename  = os.path.expandvars ( filename )
        filename  = os.path.expanduser ( filename )
        filename  = os.path.expandvars ( filename )
        filename  = os.path.expandvars ( filename )
        
        self.__compress      = False 
        self.__filename      = filename
        self.__remove        = False
        self.__silent        = silent
        self.__opened        = False 

        if not self.__silent :
            logger.info ( 'Open DB: %s' % filename ) 

        ## filename without extension and the extension  itself 
        fname , ext = os.path.splitext ( filename )

        self.__extension = ext

        ## predefined extension?
        if ext.lower() in self.extensions :

            fexists = os.path.exists ( filename )
            
            if  fexists and 'r' == mode :
                
                ## uncompress into the temporary location
                filename_ = self.uncompress_file ( filename ) 
                if not os.path.exists ( filename_ ) :
                    raise TypeError ( "Unable to uncompress properly: %s" % filename )
                if not self.__silent : 
                    size1 = os.path.getsize ( filename  ) 
                    size2 = os.path.getsize ( filename_ )
                    logger.info ( "Uncompression %s: %.1f%%" %  ( filename , ( size2 * 100.0 ) / size1 ) )
                filename        = filename_ 
                self.__filename = filename_
                self.__remove   = True
                
            elif fexists and 'r' != mode :
                
                ## uncompress in place (name witout extension)
                filename_     = fname  
                # remove existing file (if needed) 
                if os.path.exists ( filename_ ) : os.remove ( filename_ )
                size1 = os.path.getsize    ( filename ) 
                # gunzip in place 
                self.__in_place_uncompress ( filename ) 
                ##
                if not os.path.exists ( filename_ ) :
                    raise TypeError ( "Unable to uncompress properly: %s" % filename )
                if not self.__silent : 
                    size2 = os.path.getsize ( filename_ )
                    logger.info ( "Uncompression %s: %.1f%%" %  ( filename , ( size2 * 100.0 ) / size1 ) ) 
                filename        = filename_ 
                self.__compress = True 
                self.__filename = filename_
                self.__remove   = False
                
            else :
                
                ## 
                filename        = fname 
                self.__compress = True 
                self.__filename = filename

        if python_version.major >= 3 :
            
            shelve.Shelf.__init__ (
                self                              ,
                dbm.open ( self.filename , mode ) ,
                protocol                          ,
                writeback                         ,
                keyencoding                       )
        else :
            
            shelve.Shelf.__init__ (
                self                              ,
                dbm.open ( self.filename , mode ) ,
                protocol                          ,
                writeback                         ) 
            self.keyencoding = keyencoding
                    
        self.__opened        = True
        self.__mode          = mode
        
    @property
    def protocol ( self ) :
        """``protocol'' : pickling protocol used in the shelve"""
        return self._protocol
    
    @property 
    def filename ( self ) :
        "``filename'' :   the actual file name for database"
        return self.__filename

    @property 
    def opened   ( self ) :
        "``open'' : is data base opened?"
        return self.__opened

    @property
    def mode    ( self ) :
        "``mode'' : the actual open-mode for the database"
        return self.__mode
    
    @property
    def protocol( self ) :
        "``protocol'' : pickling protocol"
        return self._protocol

    # =========================================================================
    ## valid, opened DB 
    def __nonzero__ ( self ) :
        return self.opened \
               and not isinstance ( self.dict , shelve._ClosedDict ) \
               and not self.dict is None 

    # =========================================================================
    ## valid, opened DB 
    def __bool__    (  self ) : return self.__nonzero__ ()

    # =========================================================================
    ## destructor 
    def __del__ ( self ) :
        """ Destructor 
        """
        ## close if opened 
        if self.opened : self.close ()  

    # =========================================================================
    ## iterator over good keys 
    def ikeys ( self , pattern = '' ) :
        """Iterator over avilable keys (patterns included).
        Pattern matching is performed accoriding to
        fnmatch/glob/shell rules [it is not regex!] 
        
        >>> db = ...
        >>> for k in db.ikeys('*MC*') : print(k)
        
        """
        keys_ = self.keys()
        
        if not pattern :
            good = lambda s,p : True
        else :
            import fnmatch
            good = lambda s,p : fnmatch.fnmatchcase ( k , p )
        
        for k in sorted ( keys_ ) :
            if good ( k , pattern ) : yield k

    # =========================================================================
    ## list the avilable keys 
    def ls    ( self , pattern = '' ) :
        """List the available keys (patterns included).
        Pattern matching is performed accoriding to
        fnmatch/glob/shell rules [it is not regex!] 

        >>> db = ...
        >>> db.ls() ## all keys
        >>> db.ls ('*MC*')        
        
        """
        n  = os.path.basename ( self.filename )
        ap = os.path.abspath  ( self.filename ) 
        ll = getLogger ( n )
        
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
            
        ll.info ( 'Database: %s #keys: %d size: %s' % ( ap , len ( self ) , size ) )
                
        keys = [] 
        for k in self.ikeys ( pattern ): keys.append ( k )
        keys.sort()
        if keys : mlen = max ( [ len(k) for k in keys] ) + 2 
        else    : mlen = 2 
        fmt = ' --> %%-%ds : %%s' % mlen 
        for k in keys :
            ss = len ( self.dict[k] ) ##  compressed size 
            if    ss < 1024 : size = '%8d' % ss 
            elif  ss < 1024  * 1024 :
                size = '%8.3f kB' %  ( float ( ss ) / 1024 )
            elif  ss < 1024  * 1024 * 1024 :
                size = '%8.3f MB' %  ( float ( ss ) / ( 1024 * 1024 ) )
            else :
                size = '%8.3f GB' %  ( float ( ss ) / ( 1024 * 1024 * 1024 ) )
            
            ll.info ( fmt  % ( k , size ) ) 
        
        
    # =========================================================================
    ## close and compress (if needed)
    def close ( self ) :
        """ Close the file (and compress it if required) 
        """
        if not self.opened : return 
        ##
        shelve.Shelf.close ( self )
        self.__opened = False  
        ##
        if self.__remove and os.path.exists ( self.__filename ) :
            if not self.__silent :
                logger.info( 'REMOVE: ', self.__filename )
            os.remove ( self.__filename )
        ##
        if self.__compress and os.path.exists ( self.__filename ) :
            # get the initial size 
            size1 = os.path.getsize  ( self.__filename )
            # compress the file
            self.__in_place_compress ( self.__filename ) 
            #
            cname = self.__filename + self.__extension 
            if not os.path.exists ( cname ) :
                logger.warning( 'Unable to compress the file %s ' % self.__filename  )
            size2 = os.path.getsize ( cname )
            if not self.__silent : 
                logger.info( 'Compression %s: %.1f%%' % ( self.__filename, ( size2 * 100.0 ) / size1 ) )

    # =========================================================================
    ## compress the file (``in-place'') 
    def __in_place_compress   ( self , filein ) :
        """Compress the file ``in-place''        
        - It is better to use here ``os.system'' or ``popen''-family,
        but it does not work properly for multiprocessing environemnt        
        """
        cfname  =  filein + self.__extension
        # compress the file 
        fileout = self.compress_file ( filein )
        # 
        if os.path.exists ( fileout ) :
            # rename the temporary file 
            shutil.move   ( fileout , cfname )
            time.sleep    ( 2 )
            # remove the original
            os.remove     ( filein )
            
    # =========================================================================
    ## uncompress the file (``in-place'') 
    def __in_place_uncompress   ( self , filein ) :
        """Uncompress the file ``in-place''        
        - It is better to use here ``os.system'' or ``popen''-family,
        but unfortunately it does not work properly for multithreaded environemnt        
        """
        cuname , ext = os.path.splitext ( filein )
        if ( not ext ) or  ( ext not in self.extensions ) : 
            logger.error('Unknown extension for %s' % filein )    
        # uncompress the file 
        fileout = self.uncompress_file ( filein )
        if os.path.exists ( fileout ) :
            # rename the temporary file 
            shutil.move   ( fileout , cuname )
            time.sleep    ( 2 )
            # remove the original
            os.remove     ( filein ) 

    # =========================================================================
    ## some context manager functionality : enter 
    def __enter__ ( self      ) : return self
    
    # =========================================================================
    ## some context manager functionality : exit 
    def __exit__  ( self , *_ ) : self.close ()

    # =========================================================================
    ## print it! 
    def __repr__ ( self ) :
        
        kls = self.__class__
        if kls.__module__.startswith('__') and kls.__module__.endswith('__') :
            kname = kls.__name__
        else :
            kname = '%s.%s' % (  kls.__module__ , kls.__name__ ) 
        
        if   self and len ( self ) :
            return "%s('%s'): %d object(s)" % ( kname , self.filename , len(self) ) 
        elif self :
            return "%s('%s'): empty"        % ( kname , self.filename ) 
        return "Invalid/Closed %s('%s')"    % ( kname , self.filename )
    
    __str__ = __repr__

    # =========================================================================
    ## ``get-and-uncompress-item'' from dbase
    #  @code
    #  value =   dbase ['item']
    #  @endcode 
    def __getitem__  ( self , key ) : 
        """ ``get-and-uncompress-item'' from dbase
        >>> value = dbase['item'] 
        """
        try:            
            value = self.cache [ key ]            
        except KeyError:            
            value = self.uncompress_item ( self.dict [ key.encode ( self.keyencoding ) ] )            
            if self.writeback : self.cache [ key ] = value            
        return value       

    # =========================================================================
    ## ``set-and-compress-item'' to dbase
    #  @code
    #  dbase ['item'] = value 
    #  @endcode 
    def __setitem__  ( self , key , value ) :
        """ ``get-and-uncompress-item'' from dbase 
        >>> dbase['item'] = value 
        """
        if self.writeback : self.cache [ key ] = value
        self.dict [ key.encode ( self.keyencoding ) ] = self.compress_item ( value ) 

    # =========================================================================
    @abc.abstractmethod
    def compress_item   ( self , value ) :
        """Compress the value  using the certain compressing engine"""
        return None

    # =========================================================================
    @abc.abstractmethod
    def uncompress_item ( self , value ) :
        """Uncompress the value  using the certain compressing engine"""
        return None

    # =========================================================================
    ## Compress the file into temporary location, keep original
    @abc.abstractmethod 
    def compress_file   ( self , filein ) :
        """Compress the file into temporary location, keep the original """
        return None 

    ## Uncompress the file into temporary location, keep original
    @abc.abstractmethod 
    def uncompress_file ( self , filein ) :
        """Uncompress the file into temporary location, keep the original"""
        return None


# ============================================================================
## a bit more decorations for shelve  (optional)
import ostap.io.shelve_ext

    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
