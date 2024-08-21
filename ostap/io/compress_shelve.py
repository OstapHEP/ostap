#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file compress_shelve.py
# 
# Abstract helper class for "compressed" version of the shelve database.
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
#    of the whole data base, that can be rather useful for data base
#    with large amout of keys 
#
# The module has been developed and used with great success in
# `Kali, framework for fine calibration of LHCb Electormagnetic Calorimeter'
#
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2010-04-30
# =============================================================================
""" Abstract helper class for `compressed' version of shelve database.

Keeping the same interface and functionlity as shelve data base,
it allows much more compact file size through the on-flight
compression of the content

The actual code has been inspired by zipshelve ( see Google...)

However is contains several new features:

 - Optionally it is possible to perform the compression
   of the whole data base, that can be rathe useful fo data base
   with large amout of keys 
   
The module has been developed and used with great success in
 `Kali, framework for fine calibration of LHCb Electormagnetic Calorimeter'

"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2010-04-30"
__version__ = "$Revision$" 
# =============================================================================
__all__ = (
    'CompressShelf'    ,   ## The DB-itself
    'PROTOCOL'         ,
    'DEFAULT_PROTOCOL' ,
    'HIGHEST_PROTOCOL' ,
    'ENCODING'       
    )
# =============================================================================
import os, abc, shelve, shutil, glob, time, datetime
from   sys                  import version_info           as python_version
from   ostap.io.dbase       import dbopen, whichdb, Item, ordered_dict  
from   ostap.core.meta_info import meta_info 
from   ostap.io.pickling    import ( Pickler , Unpickler, BytesIO,
                                     PROTOCOL,
                                     HIGHEST_PROTOCOL, DEFAULT_PROTOCOL ) 

import pickle 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.compress_shelve' )
else                      : logger = getLogger ( __name__                   )
# =============================================================================
## encoding 
ENCODING = 'utf-8'
# =============================================================================
_modes_ = {
    # =========================================================================
    # 'r'	Open existing database for reading only
    # 'w'	Open existing database for reading and writing
    # 'c'	Open database for reading and writing, creating it if it doesn't exist
    # 'n'	Always create a new, empty database, open for reading and writing
    # =========================================================================
    'n'        : 'n' ,
    'c'        : 'c' ,
    'r'        : 'r' ,
    'w'        : 'w' ,
    'u'        : 'w' , ## update
    'a'        : 'w' , ## update/appends 
    ##
    '+'        : 'w' , ## update/append        
    'w+'       : 'w' , ## update/append
    'rw'       : 'w' , ## read/write 
    'new'      : 'n' , ## new
    'create'   : 'c' , ## create if needed  
    'recreate' : 'n' , ## new 
    'read'     : 'r' , ## read-only        
    'write'    : 'w' , ## write/update/append        
    'update'   : 'w' , ## write/update/append               
    'append'   : 'w' , ## write/update/append               
    }

# =============================================================================
## @class CompressShelf
#  Abstract class for `compressed' version of `shelve'-database.
#  It has four abstract methods:
#  - <code>compress_item</code>
#  - <code>uncompress_item</code>
#  - <code>compress_file</code>
#  - <code>uncompress_file</code>
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2010-04-30
class CompressShelf(shelve.Shelf,object):
    """ `Compressed'-version of `shelve'-database
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
            dbname                 ,
            mode        = 'c'      ,
            dbtype      = ''       , ## preferred data base 
            protocol    = PROTOCOL , ## pickle protocol 
            compress    = 0        , 
            writeback   = False    ,
            silent      = True     ,
            keyencoding = ENCODING , **kwargs ) :
        
        ## the mode 
        mode = _modes_.get( mode.lower() , '' )
        if not mode :
            logger.warning ( "Unknown opening mode '%s', replace with 'c'")
            mode = 'c'

        if not 0 <= protocol <= HIGHEST_PROTOCOL :
            logger.warning ("Invalid pickle protocol:%s" % protocol )
            protocol = PROTOCOL 

        ## expand the actual file name 
        dbname  = os.path.expandvars ( dbname )
        dbname  = os.path.expanduser ( dbname )
        dbname  = os.path.expandvars ( dbname )
        dbname  = os.path.expandvars ( dbname )
        
        self.__compresslevel = compress

        self.__nominal_dbname  = dbname 
        self.__actual_dbname   = dbname 
        self.__compress        = () 
        self.__remove          = () 
        self.__silent          = silent
        self.__opened          = False 
        self.__files           = ()
        self.__dbtype          = dbtype
        
        if not self.__silent :
            logger.info ( 'Open DB: %s' % dbname ) 

        ## filename without extension and the extension  itself 
        fname , ext = os.path.splitext ( dbname )

        self.__extension = ext

        ## predefined extension?
        if ext.lower() in self.extensions :

            fexists = os.path.exists ( dbname )
            
            if  fexists and 'r' == mode :
                
                ## uncompress into the temporary location 
                tfiles   = self.uncompress_file       ( dbname )
                filename = self.dbase_name            ( tfiles )
                
                self.__remove        = tfiles 
                self.__compress      = ()                                
                self.__actual_dbname = filename
                self.__files         = tfiles
                
            elif fexists and 'r' != mode :
                
                ## uncompress locally 
                tfiles   = self.__in_place_uncompress ( dbname )
                filename = self.dbase_name            ( tfiles )

                self.__compress      = tfiles 
                self.__remove        = tfiles 
                self.__files         = tfiles
                self.__actual_dbname = filename

            else :
                
                ## 
                filename             = fname 
                self.__compress      = True 
                self.__remove        = ()
                self.__actual_dbname = filename

        afiles = tuple ( [ self.dbname + suffix for suffix in (  '' , ',db' , '.dir' , '.pag' , '.dat' ) ] )
        ofiles = set ( [ i for i in glob.iglob  ( self.dbname + '*' ) if i in afiles ] ) 

        
        self.__opened = False

        ## actual database 
        dbase = dbopen ( self.dbname , flag = mode , dbtype = dbtype , **kwargs )
        conf  = { 'protocol' : protocol , 'writeback' : writeback }

        if  3 <= python_version.major  : conf [ 'keyencoding'] = keyencoding
        else                           : self.keyencoding      = keyencoding
        
        shelve.Shelf.__init__ ( self   , dbase  , **conf )


        self.__opened        = True
        self.__mode          = mode        

        ### self.sync  ()

        self.__dbtype        =  whichdb ( self.dbname )
        if hasattr ( self.dict , 'reopen' ) : self.dict.reopen() 
                             
        nfiles = set ( [ i for i in glob.iglob  ( self.dbname + '*' ) if i in afiles ] ) - ofiles 

        if not self.__files :
            
            files = []
            f     = self.dbname
            db    = self.dbtype           
            if os.path.exists ( f ) and os.path.isfile ( f ) and \
                   db in ( 'dbm.gnu' , 'gdbm' , 'dbhash' , 'bsddb185' , 'bsddb' , 'bsddb3' , 'sqlite3' , 'berkeleydb') :  
                files.append  ( f )
            elif   f + '.db'  in nfiles  and db in ( 'dbm.ndmb' , 'dbm' ) :
                files.append  ( f + '.db'  )
            elif   f + '.pag' in nfiles and f  + '.dir' in nfiles and db in ( 'dbm.ndbm' , 'dbm'     ) :
                files.append  ( f + '.pag' )
                files.append  ( f + '.dir' )
            elif   f + '.dat' in nfiles and f  + '.dir' in nfiles and db in ( 'dbm.dumb' , 'dumbdbm' ) :
                files.append  ( f + '.dat' )
                files.append  ( f + '.dir' )
            elif   f + '.pag' in nfiles                           and db in ( 'dbm.ndbm' , 'dbm'     ) :
                files.append  ( f + '.pag' )
            elif   f + '.db'  in ofiles  and db in ( 'dbm.ndmb' , 'dbm' ) :
                files.append  ( f + '.db'  )
            elif   f + '.pag' in ofiles and f  + '.dir' in ofiles and db in ( 'dbm.ndbm' , 'dbm'     ) :
                files.append  ( f + '.pag' )
                files.append  ( f + '.dir' )
            elif   f + '.dat' in ofiles and f  + '.dir' in ofiles and db in ( 'dbm.dumb' , 'dumbdbm' ) :
                files.append  ( f + '.dat' )
                files.append  ( f + '.dir' )                
            elif   f + '.pag' in ofiles                           and db in ( 'dbm.dumb' , 'dumbdbm' ) :
                files.append  ( f + '.pag' )
            ## else  :
            ##    logger.error ( 'Cannot find DB for %s|%s' % ( self.dbname , self.dbtype ) ) 

            files.sort ()
            self.__files = tuple  ( files  )
        
        if  self.__compress is True :
            self.__compress = self.files 
            self.__remove   = self.files 

            
        if    'sqlite3' == self.dbtype and self.mode in ( 'w' , 'n' ) : write = True
        elif  'sqlite3' != self.dbtype and self.mode in ( 'c' , 'n' ) : write = True
        else                                                          : write = False

        if write :
            dct  = ordered_dict()  
            dct  [ 'Created by'                  ] = meta_info.User
            dct  [ 'Created at'                  ] = datetime.datetime.now ().strftime( '%Y-%m-%d %H:%M:%S' )  
            dct  [ 'Created with Ostap version'  ] = meta_info.Ostap
            dct  [ 'Created with Python version' ] = meta_info.Python
            dct  [ 'Created with ROOT version'   ] = meta_info.ROOT 
            dct  [ 'Pickle protocol'             ] = protocol 
            dct  [ 'Compress level'              ] = self.__compresslevel 
            self [ '__metainfo__'                ] = dct

            
        if not self.silent :
            self.ls ()
            ff = [ os.path.basename ( f ) for f in self.files ]
            ff = ff [0] if 1 == len ( ff ) else ff              
            logger.info ( 'DB files are %s|%s' % ( ff, self.dbtype ) )

        self.sync ()
        
    @property
    def dbtype   ( self ) :
        """`dbtype'  : the underlying type of database"""
        return self.__dbtype
        
    @property
    def protocol ( self ) :
        """`protocol' : pickling protocol used in the shelve"""
        return self._protocol
    
    @property
    def compression ( self ) :
        "`compression' : compression level"
        return self.__compresslevel
    
    @property
    def compresslevel ( self ) :
        "`compress level' : compression level"
        return self.__compresslevel 

    @property 
    def dbname ( self ) :
        "`dbname' :   the actual name for the database"
        return self.__actual_dbname

    @property 
    def nominal_dbname ( self ) :
        "`nominal_dbname' :   the actual name for the database"
        return self.__nominal_dbname

    @property 
    def opened   ( self ) :
        "`open' : is data base opened?"
        return self.__opened

    @property
    def mode    ( self ) :
        "`mode' : the actual open-mode for the database"
        return self.__mode
    
    @property
    def silent ( self ) :
        "`silent' : silent actions?"
        return self.__silent 

    @property
    def protocol( self ) :
        "`protocol' : pickling protocol"
        return self._protocol

    @property
    def files  ( self ) :
        """`files' : the files assocated with the database"""
        return self.__files

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
    ##  Iterator over available keys (patterns included).
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

        if not pattern :
            good = lambda k : True
        elif regex : 
            import re
            re_cmp = re.compile ( pattern ) 
            good = lambda k  : re_cmp.match ( k )
        else :
            import fnmatch
            good = lambda s : fnmatch.fnmatchcase ( k , pattern  )


        for key in self.keys () :
            if good ( key ) : yield key

    # =========================================================================
    ##  Build the table of content(patterns included).
    #   Pattern matching is performed according to
    #  fnmatch/glob/shell rules [it is not regex!] 
    #  @code  
    #  db = ...
    #  print ( db.table() ) 
    #  print ( db.table ('*MC*') ) 
    #  @endcode 
    def table ( self , pattern = '' , load = True , prefix = '' ) :
        """Build table of content (patterns included).
        Pattern matching is performed according to
        fnmatch/glob/shell rules [it is not regex!] 

        >>> db = ...
        >>> print ( db.table () ) ## all keys
        >>> print ( db.table ('*MC*') ) 
        
        """
        n  = os.path.basename ( self.dbname         )
        ap = os.path.abspath  ( self.dbname         ) 
        ap = os.path.abspath  ( self.nominal_dbname ) 
        
        self.sync() 
        fs = float ( self.disk_size()  ) 
            
        if    fs <= 0                   : size = "???"
        elif  fs <  1024                : size = '%3d'    %  int ( fs )
        elif  fs <  1024  * 1024        : size = '%.2fkB' % ( fs / 1024 )
        elif  fs <  1024  * 1024 * 1024 : size = '%.2fMB' % ( fs / ( 1024 * 1024 ) )
        else                            : size = '%.2fGB' % ( fs / ( 1024 * 1024 * 1024 ) )

        keys = [] 
        for k in self.ikeys ( pattern ): keys.append ( k )
        keys.sort()

        if keys : mlen = max ( [ len ( k ) for k in keys ] ) + 2 
        else    : mlen = 2 
        fmt = ' --> %%-%ds : %%s' % mlen

        table = [ ( 'Key' , 'type',  '   size   ' , ' created/modified') ]

        meta = self.get( '__metainfo__' , ordered_dict() )
        for k in meta :
            row = "META:%s" % k , '' , '' , str ( meta[k] )
            table.append ( row  ) 
        
        for k in keys :

            ss = len ( self.dict [ k.encode ( self.keyencoding ) ] )
            
            if    ss < 1024 : size = '%7d   ' % ss 
            elif  ss < 1024 * 1024 :
                size = '%7.2f kB' %  ( float ( ss ) / 1024 )
            elif  ss < 1024 * 1024 * 1024 :
                size = '%7.2f MB' %  ( float ( ss ) / ( 1024 * 1024 ) )
            else :
                size = '%7.2f GB' %  ( float ( ss ) / ( 1024 * 1024 * 1024 ) )

            rawitem = self.__get_raw_item__ ( k )
            if isinstance ( rawitem , Item ) :
                timetag = rawitem.time
                if isinstance ( timetag , float ) and 0 < timetag :
                    timetag = datetime.datetime.fromtimestamp ( timetag  )
                if isinstance ( timetag , datetime.datetime ) :
                    timetag = timetag.strftime( '%Y-%m-%d %H:%M:%S' )
                kobj    = rawitem.payload 
            else                             :
                timetag = ''
                kobj    = rawitem  

            ot   = type ( kobj )
            if   hasattr ( ot , '__cpp_name__' ) : otype = ot.__cpp_name__ 
            elif hasattr ( ot , '__cppname__'  ) : otype = ot.__cppname__
            else                                 : otype = ot.__name__ 

            row = '{:15}'.format ( k ) , '{:15}'.format ( otype ) , size  , timetag 
            table.append ( row )
            
        import ostap.logger.table as T
        t      = type( self ).__name__
        title  = '%s:%s' % ( t  , n )
        maxlen = 0
        for row in table :
            rowlen = 0 
            for i in row : rowlen += len ( i )
            maxlen = max ( maxlen, rowlen ) 
        if maxlen + 3 <= len ( title ) :
            title = '<.>' + title [ -maxlen : ] 

        return T.table ( table , title = title , prefix = prefix , alignment = 'llcc')

    # =========================================================================
    ## List the available keys (patterns included).
    #   Pattern matching is performed according to
    #  fnmatch/glob/shell rules [it is not regex!] 
    #  @code  
    #  db = ...
    #  db.ls() ## all keys
    #  db.ls ('*MC*')        
    #  @endcode 
    def ls    ( self , pattern = '' , load = True , prefix = '# ' , logger = None ) :
        """List the available keys (patterns included).
        Pattern matching is performed according to
        fnmatch/glob/shell rules [it is not regex!] 

        >>> db = ...
        >>> db.ls() ## all keys
        >>> db.ls ('*MC*')        
        
        """

        n  = os.path.basename ( self.dbname         )
        ap = os.path.abspath  ( self.dbname         ) 
        ap = os.path.abspath  ( self.nominal_dbname ) 

        self.sync() 
        fs = self.disk_size()  
            
        if    fs <= 0            : size = "???"
        elif  fs <  1024         : size = str ( fs ) 
        elif  fs <  1024  * 1024 :
            size = '%.2fkB' % ( float ( fs ) / 1024 )
        elif  fs <  1024  * 1024 * 1024 :
            size = '%.2fMB' % ( float ( fs ) / ( 1024 * 1024 ) )
        else :
            size = '%.2fGB' % ( float ( fs ) / ( 1024 * 1024 * 1024 ) )

        t     = type( self ).__name__
        tab   = self.table ( pattern = pattern , load = load , prefix = prefix )
        
        ll    = logger if logger else getLogger ( n )
        line  = 'Database %s:%s|%s #keys: %d size: %s' % ( t , ap , self.dbtype , len ( self ) , size )
        ll.info (  '%s\n%s' %  ( line , tab ) )
                
    # =========================================================================
    ## close and compress (if needed)
    def close ( self ) :
        """ Close the file (and compress it if required) 
        """
        
        print ('CLOSE', self.opened )
        
        if not self.opened : return 
        ##
        if    'sqlite3' == self.dbtype and 'c' == self.mode : write = True
        elif  'sqlite3' != self.dbtype and 'e' == self.mode : write = True
        else                                                : write = False 
        if write : 
            dct  = self.get ( '__metainfo__' , ordered_dict () )
            dct  [ 'Updated at'                  ] = datetime.datetime.now().strftime( '%Y-%m-%d %H:%M:%S' )   
            dct  [ 'Updated by'                  ] = meta_info.User 
            dct  [ 'Updated with Ostap version'  ] = meta_info.Ostap 
            dct  [ 'Updated with Python version' ] = meta_info.Python 
            dct  [ 'Updated with ROOT version'   ] = meta_info.ROOT
            
            ## self [ '__metainfo__'                ] = dct

        if not self.silent : self.ls ()
        
        shelve.Shelf.close ( self )
        self.__opened = False  
        ##
        if self.__compress :
            self.__in_place_compress ( self.__compress ) 

        ##  remove the intermediate files 
        for f in self.__remove :
            if os.path.exists ( f ) :
                try  :
                    os.remove ( f )
                except OSError :
                    pass
                
    # =========================================================================
    ## compress the files (`in-place') 
    def __in_place_compress   ( self , files  ) :
        """Compress the file `in-place'        
        - It is better to use here `os.system' or `popen'-family,
        but it does not work properly for multiprocessing environemnt        
        """
        output   = self.nominal_dbname 
        out , _  = os.path.split ( output )
        outdir   = out if out else  '.'
        assert os.access ( outdir , os.W_OK  ) ,\
               'The directory "%s" is not writeable!' % os.abspath ( outdir )
        # 
        # compress the file 
        compressed = self.compress_files ( files )
        # remove input files  
        for f in files  :
            try :
                os.remove  ( f  )
            except OSError :
                pass
            
        shutil.move ( compressed , output )
            
    # =========================================================================
    ## uncompress the file (`in-place') 
    def __in_place_uncompress   ( self , filein ) :
        """Uncompress the file `in-place'        
        - It is better to use here `os.system' or `popen'-family,
        but unfortunately it does not work properly for multithreaded environemnt        
        """
        _ , ext = os.path.splitext ( filein )
        if ( not ext ) or  ( ext not in self.extensions ) : 
            logger.error ( 'Unknown extension for %s' % filein )
        ##
        out , _  = os.path.split   ( filein )
        outdir   = out if out else  '.'
        assert os.access ( outdir , os.W_OK  ) ,\
               'The directory "%s" is not writeable!' % os.abspath ( outdir )
        
        ## uncompress the file
        tfiles = self.uncompress_file ( filein )
        # remove the original
        os.remove     ( filein ) 
        ofiles      = [] 
        for f in tfiles :
            _ , ff  = os.path.split ( f        )
            ff      = os.path.join  ( out , ff )  
            shutil.move   ( f   , ff   )
            ofiles.append ( ff )
            
        return tuple  ( ofiles ) 
    
    # =========================================================================
    ## Context manager functionality : enter 
    def __enter__ ( self      ) : return self
    
    # =========================================================================
    ## Context manager functionality : exit 
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
            return "%s('%s')|%s: %d object(s)" % ( kname , self.dbname , self.dbtype  , len ( self ) ) 
        elif self :
            return "%s('%s')|%s: empty"        % ( kname , self.dbname , self.dbtype  )
        
        return "Invalid/Closed %s('%s')"       % ( kname , self.dbname )
    
    __str__ = __repr__

    # =========================================================================
    ## `get-and-uncompress-item' from dbase
    #  @code
    #  value =   dbase ['item']
    #  @endcode 
    def __get_raw_item__ ( self , key ) : 
        """ `get-and-uncompress-item' from dbase
        >>> value = dbase['item'] 
        """
        try:            
            value = self.cache [ key ]
        except KeyError:            
            value = self.uncompress_item ( self.dict [ key.encode ( self.keyencoding ) ] ) 
            if self.writeback : self.cache [ key ] = value
            
        return value
    
    # =========================================================================
    ## `get-and-uncompress-item' from dbase
    #  @code
    #  value =   dbase.get ( 'item' , default ) 
    #  @endcode 
    def get  ( self , key , default = None ) : 
        """ `get-and-uncompress-item' from dbase
        >>> value =   dbase.get ( 'item' , default ) 
        """
        try : 
            value = self.__get_raw_item__ ( key )
            if isinstance ( value , Item ) : value = value.payload
            return value 
        except KeyError :
            return default

    # =========================================================================
    ## `get-and-uncompress-item' from dbase
    #  @code
    #  value =   dbase ['item']
    #  @endcode 
    def __getitem__  ( self , key ) : 
        """ `get-and-uncompress-item' from dbase
        >>> value = dbase['item'] 
        """
        value = self.__get_raw_item__ ( key )        
        if isinstance ( value , Item ) : value = value.payload 
        return value
    
    # =========================================================================
    ## `set-and-compress-item' to dbase
    #  @code
    #  dbase ['item'] = value 
    #  @endcode 
    def __setitem__  ( self , key , value ) :
        """ `get-and-uncompress-item' from dbase 
        >>> dbase['item'] = value 
        """
        item = Item ( time.time()  , value )
        if self.writeback : self.cache [ key ] = item
        ##
        self.dict [ key.encode ( self.keyencoding ) ] = self.compress_item ( item )

    # =========================================================================
    ##  get the disk size of the db
    #   @code
    #   db = ...
    #   ds = db.disk_size()   
    #   @endcode  
    def disk_size ( self  ) :
        """Get the disk size of the db
        >>> db = ...
        >>> ds = db.disk_size()   
        """
        size  = 0 
        for f in self.files  :
            if os.path.exists ( f  ) and os.path.isfile ( f ) :
                size  += os.path.getsize ( f )
        return size 
        
    # =========================================================================
    ## Create the temporary directory
    #  The directory will be cleaned-up and deleted at-exit.
    @classmethod
    def tempdir ( cls , suffix = '=db-dir' , prefix = 'ostap-compress-shelve-dir-' , date = True  ) :
        """ Create the temporary directory
        The directory will be cleaned-up and deleted at-exit.
        """
        from ostap.utils.cleanup import CleanUp as CU
        return CU.tempdir ( suffix = suffix , prefix = prefix, date = date ) 

    # =========================================================================
    ## Ccreate the name for the temproary file 
    #  The file will be deleted at-axit 
    @classmethod
    def tempfile ( cls , suffix = '-db' , prefix = 'ostap-compress-shelve-' , dir = None , date = True  ) :
        """ Create the name for the temporary file 
        The file will be deleted at-axit
        """
        from ostap.utils.cleanup import CleanUp as CU
        return CU.tempfile ( suffix = suffix , prefix = prefix, dir = dir , date = date ) 

    # ========================================================================
    ## guess the name of the database from the list of (uncompressed files)
    @classmethod
    def dbase_name ( cls , files ) :
        """ Guess the name of the database from the list of (uncompressed files)
        """
        
        exts = set ( [ os.path.splitext ( f )[1] for f in files ] )
        
        if   1 == len ( files ) and whichdb ( files [ 0 ] ) : return files [ 0 ]
        elif 1 == len ( files ) and '.db' in exts :
            f , _ = os.path.splitext ( files [0] )
            if whichdb ( f  ) : return f
        elif 2 <= len ( files  ) and  '.dir' in exts and '.pag' in exts :
            f , _ = os.path.splitext ( files [0] )
            if whichdb ( f  ) : return f
        elif 1 <= len ( files  )                     and '.pag' in exts :
            f , _ = os.path.splitext ( files [0] )
            if whichdb ( f  ) : return f
        elif 2 <= len ( files  ) and  '.dir' in exts and '.dat' in exts :
            f , _ = os.path.splitext ( files [0] )
            if whichdb ( f  ) : return f

        raise TypeErrro ('Cannot identify the database name: %s' % str ( files ) )  

    # =========================================================================
    ## Pickle/serialize compressed data 
    def pickle ( self , value ) :
        """Pickle/serialize compressed data"""
        f = BytesIO ()
        p = Pickler ( f , self.protocol )
        p.dump ( value )
        return f.getvalue()
    
    # =========================================================================
    ## Unpickle/deserialize the uncompressed data
    def unpickle ( self , value ) :
        """Unpickle/deserialize uncompressed data"""
        f = BytesIO ( value )
        return Unpickler ( f ) . load ( )
    
    # =========================================================================
    @abc.abstractmethod
    def compress_item   ( self , value ) :
        """Compress the value  using the certain compressing engine"""
        return NotImplemented

    # =========================================================================
    @abc.abstractmethod
    def uncompress_item ( self , value ) :
        """Uncompress the value  using the certain compressing engine"""
        return NotImplemented

    # =========================================================================
    ## Compress the files into temporary location, keep original
    @abc.abstractmethod 
    def compress_files  ( self , files  ) :
        """Compress the files into temporary location, keep the original """
        return NotImplemented

    # =========================================================================
    ## Uncompress the file into temporary location, keep original
    @abc.abstractmethod 
    def uncompress_file ( self , filein ) :
        """Uncompress the file into temporary location, keep the original"""
        return NotImplemented

    # =========================================================================
    ## clone the database into new one
    #  @code
    #  db  = ...
    #  ndb = db.clone ( 'new_file.db' )
    #  @endcode
    @abc.abstractmethod
    def clone ( self , filename , keys = () ) :
        """ Clone the database into new one
        >>> old_db = ...
        >>> new_db = new_db.clone ( 'new_file.db' )
        """
        return NotImplemented

# ============================================================================
## a bit more decorations for shelve  (optional)
import ostap.io.shelve_ext

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
