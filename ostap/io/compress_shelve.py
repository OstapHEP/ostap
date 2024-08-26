#!/usr/bin/env python
#*- coding: utf-8 -*-
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
#    of the whole data base, that can be useful for data base
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
 `Kali, framework for fine calibration of LHCb Electromagnetic Calorimeter'

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
import os, abc, shelve, shutil, glob, time, datetime, zipfile, tarfile 
from   sys                  import version_info           as python_version
from   ostap.io.dbase       import dbopen , whichdb, Item, ordered_dict, dbfiles   
from   ostap.core.meta_info import meta_info 
from   ostap.io.pickling    import ( Pickler , Unpickler, BytesIO,
                                     PROTOCOL,
                                     HIGHEST_PROTOCOL, DEFAULT_PROTOCOL ) 
from   ostap.utils.cleanup  import CUBase
from   ostap.utils.utils    import file_size
from   ostap.utils.basic    import writeable 
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
class CompressShelf (shelve.Shelf,CUBase) :
    """ `Compressed' - version of `shelve'-database
    It has four abstract methods:
    - compress_item
    - uncompress_item
    - compress_file
    - uncompress_file
    """
    __metaclass__ = abc.ABCMeta
    ZIP_EXTS      = ( '.zip' , '.zipdb' , '.dbzip' , '.zdb' , '.dbz' ) ## whole DB is in zip-archive 
    TAR_EXTS      = ( '.tar' , '.tardb' , '.dbtar' , '.tdb' , '.dbt' ) ## whole DB is in tar-archive 

    def __init__(
            self                   ,
            dbname                 ,
            mode        = 'c'      ,
            compress    = 0        , 
            protocol    = PROTOCOL , ## pickle protocol 
            dbtype      = ''       , ## preferred data base type 
            writeback   = False    ,
            silent      = True     ,
            keyencoding = ENCODING , **kwargs ) :
        
        ## the mode 
        mode = _modes_.get( mode.lower() , '' )
        if not mode :
            logger.warning ( "Unknown opening mode '%s', replace with 'c'")
            mode = 'c'

        if not 0 <= protocol <= HIGHEST_PROTOCOL :
            logger.warning ("Invalid pickle protocol:%s, replace with:%s" % ( protocol , PROTOCOL ) ) 
            protocol = PROTOCOL 
            
        self.__opened          = False        
        
        ## expand the actual file name
        dbname_ = dbname
        dbname  = self.name_expand ( dbname ) ## from CUBase 
        
        self.__compresstype    = kwargs.pop ( 'compresstype' , '???' ) 
        self.__compresslevel   = compress
        self.__silent          = silent
        self.__protocol        = protocol
        self.__dbtype          = dbtype
        
        ## all arguments for constructor 
        self.__kwargs = {
            'dbname'       : dbname              ,
            'mode'         :  mode               ,
            'compress'     :  self.compresslevel ,
            'protocol'     : protocol            , ## from shelve.Shelf  
            'writeback'    : writeback           , ## from shelf.Shelf 
            'dbtype'       : self.dbtype         , ## preferred dbtype 
            'silent'       : self.silent         ,
            'keyencoding'  : keyencoding         ,   
        }
        self.__kwargs.update ( kwargs )

        self.__mode            = mode        
        self.__nominal_dbname  = dbname

        ## generic case 
        self.__actual_dbname   = dbname        
        self.__compress        = () 
        self.__remove          = () 
        self.__files           = ()

        if not self.__silent : logger.info ( 'Open DB: %s' % dbname ) 

        ## filename without extension and the extension
        fname , ext = os.path.splitext ( dbname )

        extlow = ext.lower ()
        
        exists = os.path.exists            ( dbname ) 
        fexist = exists and os.path.isfile ( dbname )
        
        ## use zip-compression of whole db 
        zip = extlow in self.ZIP_EXTS 
        
        ## use tar-compression of whole db 
        tar = extlow in self.TAR_EXTS 
        
        zipf   = zip and fexist and zipfile.is_zipfile ( dbname )
        tarf   = tar and fexist and tarfile.is_tarfile ( dbname )
                
        self.__zip = zip
        self.__tar = tar 
    
        if 'r' == mode and ( tarf or zipf ) : 
            
            ## uncompress into the temporary location
            tmpdir   = self.tempdir () 
            tfiles   = self.uncompress_file ( dbname , where = tmpdir )

            self.__compress      = ()                                
            self.__remove        = tfiles 
            self.__files         = tfiles
            self.__actual_dbname = os.path.join ( tmpdir , os.path.basename ( fname ) ) 

        elif 'r' != mode and ( tarf or zipf ) :

            ## uncompress locally
            outdir , _ = os.path.split        ( os.path.abspath ( dbname ) ) 
            tfiles     = self.uncompress_file ( dbname , where = outdir  )

            self.__compress      = tuple ( sorted ( tfiles ) ) 
            self.__remove        = tfiles 
            self.__files         = tfiles
            self.__actual_dbname = os.path.join ( outdir , os.path.basename ( fname ) ) 

        elif ( zip or tar ) : 
            
            filename             = fname 
            self.__compress      = True 
            self.__remove        = ()
            self.__actual_dbname = filename

        else : 
            
            filename             = dbname 
            self.__compress      = False  
            self.__remove        = ()
            self.__actual_dbname = filename

        # =====================================================================
        the_path = lambda s : os.path.normpath ( os.path.abspath ( s ) )
        ## all files  before dbopen
        ofiles  = set ( [ the_path ( i )  for i in glob.iglob  ( self.dbname + '*'  ) ] ) 
        ofiles |= set ( [ the_path ( i ) for i in glob.iglob  ( self.dbname + '/*' ) ] )         
        
        self.__ofiles = tuple ( sorted ( ofiles ) ) 
        
        # =====================================================================
        ## open/create the actual underlying database 
        # =====================================================================
        self.__opened = False        
        dbase = dbopen ( self.dbname , flag = mode , dbtype = dbtype , **kwargs )
        self.__opened        = True

        # =======================================================================
        ## all files after dbopen 
        pfiles  = set ( [ the_path ( i ) for i in glob.iglob  ( self.dbname + '*'  ) ] ) 
        pfiles |= set ( [ the_path ( i ) for i in glob.iglob  ( self.dbname + '/*' ) ] )
        ## new files 
        nfiles = pfiles - ofiles
        
        self.__pfiles = tuple ( sorted ( pfiles ) ) 
        self.__nfiles = tuple ( sorted ( nfiles ) ) 

        # =====================================================================
        ## initialize the base class
        # =====================================================================
        conf  = { 'protocol' : protocol , 'writeback' : writeback }
        if  3 <= python_version.major  : conf [ 'keyencoding'] = keyencoding
        else                           : self.keyencoding      = keyencoding        
        shelve.Shelf.__init__ ( self   , dbase  , **conf )

        # ======================================================================
        ## actual type of underlying database 
        self.__dbtype        =  whichdb ( self.dbname ) ## actual dbtype 
        
        # ======================================================================
        ## expected files for the given DB type 
        efiles = set ( the_path ( f ) for f in dbfiles ( self.dbtype , self.dbname ) ) 

        if  ( ofiles | efiles ) != pfiles :
            logger.warning ( 'Some missing or unexpected files' )
            
        files1 = pfiles & efiles ## expected and found 
        files2 = efiles - pfiles ## expected but not found 
        files3 = nfiles - efiles

        if files2 : logger.warning ( 'Expected but not found %s/%d: %s' % ( self.dbtype , len ( files2 ) , list ( files2 ) ) )            
        if files3 : logger.warning ( 'New but not expected   %s/%d: %s' % ( self.dbtype , len ( files3 ) , list ( files3 ) ) ) 

        # =====================================================================
        ## list of files (expected) 
        self.__files = tuple ( sorted ( efiles ) ) 
        
        if  self.__compress is True :
            self.__compress = self.files 
            self.__remove   = self.files 

        if 'n' == self.mode or ( 'c' == mode and not '__metainfo__' in self ) : 
            dct  = ordered_dict()  
            dct  [ 'Created by'                  ] = meta_info.User
            dct  [ 'Created at'                  ] = datetime.datetime.now ().strftime( '%Y-%m-%d %H:%M:%S' )  
            dct  [ 'Created with Ostap version'  ] = meta_info.Ostap
            dct  [ 'Created with Python version' ] = meta_info.Python
            dct  [ 'Created with ROOT version'   ] = meta_info.ROOT 
            dct  [ 'Pickle protocol'             ] = protocol
            dct  [ 'Compress level'              ] = self.compresslevel
            dct  [ 'Compress type'               ] = self.compresstype 
            dct  [ 'Underlying dbase type'       ] = self.dbtype 
            self [ '__metainfo__'                ] = dct

        if 'r' != self.mode and 'n' != self.mode :            
            dct  = self.get ( '__metainfo__' , ordered_dict () )
            if self.protocol      != dct.get ( 'Pickle protocol'       , 0  ) : dct [ 'Pickle protocol'       ] = self.protocol            
            if self.compresslevel != dct.get ( 'Compress level'        , 0  ) : dct [ 'Compress level'        ] = self.compresslevel
            if self.compresstype  != dct.get ( 'Compress type'         , '' ) : dct [ 'Compress type'         ] = self.compresstype
            if self.dbtype        != dct.get ( 'Underlying dbase type' , '' ) : dct [ 'Underlying dbase type' ] = self.dbtype
            self [ '__metainfo__' ] = dct

        if not self.silent :
            self.ls ()
            ff = [ os.path.basename ( f ) for f in self.files ]
            ff = ff [ 0 ] if 1 == len ( ff ) else ff              
            logger.info ( 'DB files are %s|%s' % ( ff , self.dbtype ) )

        self.sync ()

        self.__taropts = 'x:gz' if (3,0)<= python_version else 'w:gz'
	

    @property
    def dbtype   ( self ) :
        """`dbtype'  : the underlying type of database"""
        return self.__dbtype
        
    @property
    def protocol ( self ) :
        """`protocol' : pickling protocol used in the shelve"""
        return self.__protocol
    
    @property
    def compression ( self ) :
        """`compression' : compression level"""
        return self.__compresslevel
    
    @property
    def compresslevel ( self ) :
        """`compress level' : compression level"""
        return self.__compresslevel 

    @property
    def compresstype ( self ) :
        """`compresstype' : type of compression"""
        return self.__compresstype
    
    @property 
    def dbname ( self ) :
        """`dbname' :   the actual name for the database"""
        return self.__actual_dbname

    @property 
    def nominal_dbname ( self ) :
        """`nominal_dbname' :   the actual name for the database"""
        return self.__nominal_dbname

    @property 
    def opened   ( self ) :
        """`opened' : is data base opened?"""        
        return self.__opened
    
    @property
    def closed ( self  ) :
        """`closed` : is database closed (==not-opened)? """
        return not self.opened 

    @property
    def mode    ( self ) :
        """`mode' : the actual open-mode for the database"""
        return self.__mode
    
    @property
    def silent ( self ) :
        """`silent' : silent actions?"""
        return self.__silent 

    @property
    def protocol( self ) :
        "`protocol' : pickling protocol"
        return self.__protocol

    @property
    def files  ( self ) :
        """`files' : the files assocated with the database"""
        return self.__files

    @property
    def kwargs ( self )  :
        """`kwargs` : all constructor arguments"""
        return self.__kwargs

    @property
    def zip ( self ) :
        """`zip` : use zip-archive """
        return self.__zip

    @property
    def tar ( self ) :
        """`tar` : use tar-archive """
        return self.__tar

    @property
    def taropts ( self ) :
        """`taropts` : options to open tar-fiel for write"""
        return self.__taropts
    @taropts.setter
    def taropts ( self , value ) :
        if python_version < (3,0) and value and 'x' == value[0] : value = 'w' + value[1:] 
        self.__taropts = value
        
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

        if not pattern : good = lambda k : True
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
    #  Pattern matching is performed according to
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
        if not self.opened : return 
        ##
        if  not self.silent : logger.info   ( 'Closing database %s' % self.dbname )
        else                : logger.debug  ( 'Closing database %s' % self.dbname )
        # 
        if 'r' != self.mode and 'n' != self.mode :
            dct  = self.get ( '__metainfo__' , ordered_dict () )
            dct  [ 'Updated at'                  ] = datetime.datetime.now().strftime( '%Y-%m-%d %H:%M:%S' )   
            dct  [ 'Updated by'                  ] = meta_info.User
            dct  [ 'Updated with Ostap version'  ] = meta_info.Ostap 
            dct  [ 'Updated with Python version' ] = meta_info.Python 
            dct  [ 'Updated with ROOT version'   ] = meta_info.ROOT
            ##
            if self.protocol != dct.get ( 'Pickle protocol'            , 0  ) : 
                dct  [ 'Pickle protocol'             ] = self.protocol                        
            if self.compresslevel != dct.get ( 'Compress level'        , 0  ) : 
                dct  [ 'Compress level'              ] = self.compresslevel
            ## 
            if self.compresstype  != dct.get ( 'Compress type'         , '' ) : 
                dct  [ 'Compress type'               ] = self.compresstype
            ## 
            if self.dbtype        != dct.get ( 'Underlying dbase type' , '' ) : 
                dct  [ 'Underlying dbase type'       ] = self.dbtype
            ## 
            self [ '__metainfo__'                ] = dct

        if not self.silent : self.ls ()
        
        shelve.Shelf.close ( self )
        self.__opened = False
        
        ##
        if self.__compress :
            self.compress_files ( self.__compress , self.nominal_dbname ) 
            
        ##  remove the intermediate files 
        for f in self.__remove : self.remove_file ( f )

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
        return file_size ( *self.files )
        
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
        """ Uncompress the value  using the certain compressing engine"""
        return NotImplemented

    # =========================================================================
    ## Compress the files into specified location
    def compress_files  ( self , files  , output  ) :
        """ Compress the files into the specified location
        """
        if not self.silent : logger.info ( 'Compress %s into %s' % ( files , output ) )

        fdir   = '' 
        fdirs  = [ f for f in files if os.path.isdir ( f ) ]
        if fdirs :
            a , b = os.path.split ( fdirs  [ 0 ] )
            if not b :
                a , b = os.path.split ( a )
            fdir = os.path.join ( b , '' ) 

        if self.zip :
            with zipfile.ZipFile ( output , 'w' , allowZip64 = True ) as zfile :
                for f in sorted ( files ) :
                    _ , name = os.path.split ( f  )
                    if fdir : name = os.path.join ( fdir , name ) 
                    zfile.write ( f  , name  )        
                if not self.silent :
                    logger.info ( "Zip-file `%s` content:" % output )  
                    for f in zfile.infolist() :
                        logger.info ( '%s' % f )                    
                return output
        
        elif self.tar :
            
            with tarfile.open ( output , self.taropts ) as tfile :
                for f in sorted ( files ) :
                    _ , name = os.path.split ( f )
                    if fdir : name = os.path.join ( fdir , name )                         
                    tfile.add ( f  , name  )
                if not self.silent :
                    logger.info ( "Tar-file `%s` content:" % output )
                    tfile.list()
                return output
            
        return None  
        
    # =========================================================================
    ## Uncompress the file into specified location, keep original
    def uncompress_file ( self , filein , where ) :
        """ Uncompress the file into specofed location, keep the original """

        if not self.silent : logger.info ( 'Uncompress %s into %s' % ( filein , where ) ) 

        assert os.path.exists ( filein ) and os.path.isfile ( filein ) , \
            "Non existing/invalid file:`%s'" % filein

        assert os.path.exists ( where ) and os.path.isdir ( where ) and writeable ( where ) ,\
            "Invalid/nonwriteable  directory:`%s'" % where 

        items  = []
        
        ## 1) zip-archive ? 
        if zipfile.is_zipfile ( filein ) :
            with zipfile.ZipFile ( filein , 'r' , allowZip64 = True ) as zfile :
                if not self.silent :
                    logger.info ( "Zip-file `%s` content:" % filein )  
                    for f in zfile.infolist() :
                        logger.info ( '%s' % f ) 
                for item in zfile.filelist :
                    zfile.extract ( item , path = where )
                    name = ( os.path.join ( where , item.filename ) )
                    name = os.path.normpath  ( name ) 
                    items.append  ( name )
            return tuple  ( sorted ( items ) ) 
        
        ## 2) compressed-tar archive ?
        if tarfile.is_tarfile ( filein ) :
            with tarfile.open ( filein  , 'r:*' ) as tfile :
                if not self.silent :
                    logger.info ( "Tar-file `%s` content:" % filein )
                    tfile.list()
                for item in tfile  :
                    tfile.extract ( item , path = where  )
                    name = ( os.path.join ( where , item.name ) )
                    name = os.path.normpath  ( name ) 
                    items.append  ( name )
            return tuple ( sorted ( items ) ) 

        return () 

    # =========================================================================
    ## copy the database into new one
    #  @code
    #  db  = ...
    #  ndb = db.copy ( 'new_file.db' , copykeys , **kwargs  )
    #  @endcode
    def copy ( self , dbname , copykeys = () , **kwargs ) :
        """ Clone the database into new one
        >>> old_db = ...
        >>> new_db = new_db.clone ( 'new_file.db' )
        """
        klass = type ( self )
        conf = {}
        conf.update ( self.__kwargs ) ## use argument 
        conf.update ( kwargs        ) ## redefine them
        
        conf [ 'mode'   ] = 'c'
        conf [ 'dbname' ] = dbname
        
        ## create newdb
        newdb = klass ( **conf )
        
        ## copy the required keys: 
        for key in copykeys : newdb [ key  ] = self [ key ] 

        newdb.sync()        
        return newdb 
    
    # =========================================================================
    ## clone the database into new one (copy all keys) 
    #  @code
    #  db  = ...
    #  ndb = db.clone ( 'new_file.db' )
    #  @endcode
    def clone ( self , dbname , **kwargs ) :
        """ Clone the database into new one (copy all keys) 
        >>> old_db = ...
        >>> new_db = new_db.clone ( 'new_file.db' )
        """
        return self.copy ( dbname , copykeys = self.keys ()  , **kwargs ) 


# ============================================================================
## a bit more decorations for shelve  (optional)
import ostap.io.shelve_ext


# =============================================================================
if python_version < ( 3 , 2 ) :
    
    def _arxiv_enter_ ( s ) : return s
    def _arxiv_exit_  ( s , *_) : s.close() 
    
    import zipfile
    if not hasattr ( zipfile.ZipFile , '__enter__' ) : zipfile.ZipFile . __enter__ = _arxiv_enter_
    if not hasattr ( zipfile.ZipFile , '__exit__'  ) : zipfile.ZipFile . __exit__  = _arxiv_exit_
    
    import tarfile
    if not hasattr ( tarfile.TarFile , '__enter__' ) : tarfile.TarFile . __enter__ = _arxiv_enter_
    if not hasattr ( tarfile.TarFile , '__exit__'  ) : tarfile.TarFile . __exit__  = _arxiv_exit_
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
