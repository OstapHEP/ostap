#/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/utils/cleanup.py
#  Module with some simple but useful utilities
#  to deal  with temporary files and directories
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10  
# =============================================================================
"""Module with some simple but useful utilities
to deal  with temporary files and directories 
Module with some simple but useful utilities for
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
__all__     = (
    'CleanUp'  , ## Simple (base) class to deal with temporary files and directories
    'TempFile' , ## Simple placeholder for temporary file  
    )
# =============================================================================
import os, tempfile, datetime, weakref, re    
from   sys import version_info as python_version 
# =============================================================================
from   ostap.core.ostap_types import string_types
from   ostap.utils.basic      import make_dir, writeable, whoami, mtime 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.cleanup' )
else                       : logger = getLogger( __name__              )
del getLogger
# =============================================================================
date_format = "%Y-%b-%d"
re_format   = re.compile ( r'(?P<date>(\d{4}-\S{3}-\d{2}))-(?P<pid>(\d{2,12}))-.*' )
start       = datetime.datetime.now()
time_delta  = datetime.timedelta ( days = 4 ) ## from Friday to Tuesday :-)

# =============================================================================
user = whoami ()
# =============================================================================            
## temporary directory for <code>tempfile</code> module
base_tmp_dir = None
for_cleanup  = False 
# =============================================================================

## 1) check the environment variable OSTAP_TMPDIR 
if not base_tmp_dir :
    from ostap.utils.basic import get_env as ostap_getenv  
    ## base_tmp_dir = os.environ.get  ( 'OSTAP_TMP_DIR' , None )
    base_tmp_dir = ostap_getenv ( 'OSTAP_TMP_DIR' , None )
    ##
    if base_tmp_dir and not os.path.exists ( base_tmp_dir ) :
        base_tmp_dir = make_dir  ( base_tmp_dir ) 
    if base_tmp_dir and not writeable ( base_tmp_dir ) :
        logger.warning ("Directory `%s' is not writeable!" % base_tmp_dir )
        base_tmp_dir = None
        
## 2) get from configuration file 
if not base_tmp_dir :
    ## 2) check the configuration file 
    import ostap.core.config as OCC 
    base_tmp_dir = OCC.general.get ( 'TMP_DIR' , None )
    del OCC

    if base_tmp_dir and not os.path.exists ( base_tmp_dir ) :
        base_tmp_dir = make_dir ( base_tmp_dir ) 
    if base_tmp_dir and not writeable ( base_tmp_dir ) :
        logger.warning ('Directory ``%s'' is not writeable!' % base_tmp_dir )
        base_tmp_dir = None

# ===============================================================================
## local storage of temporary pid-dependent temporary directories 
base_tmp_pid_dirs = {}
dir_prefix        = 'ostap-session-'
dir_user_prefix   = 'ostap-session-%s-' % user 
# ===========================================================================
## create the base temporary directory
def make_base_tmp_dir () :
    """Create the base temporary directory
    """
    
    prefix = '%s' % dir_prefix 
    
    td = tempfile.gettempdir()
    if user and not user in td : prefix = dir_user_prefix
    
    prefix  = "%s%s-%d-"   %  ( prefix , start.strftime ( date_format ) , os.getpid () )

    t = tempfile.mkdtemp ( prefix = prefix ) 
    return t


# ===============================================================================
## get the process-dependent name of the temporary directory 
def tmp_dir ( pid = None ) :
    """get the process-dependent name of the temporary directory
    """
    if base_tmp_dir :
        return base_tmp_dir 
        
    if not pid : pid = os.getpid()

    if not pid in base_tmp_pid_dirs :
        piddir = make_base_tmp_dir ()
        base_tmp_pid_dirs [ pid ] = piddir
        return piddir 
            
    return base_tmp_pid_dirs [ pid ]

# ===============================================================================
## Context manager to define/redefine temporary directory for <code>tempfile</code> module
class UseTmpDir ( object ) :
    """Context manager to define/redefine TmpDir for the tempfile module
    """
    def __init__   ( self , temp_dir = None ) :
        
        self.__tmp_dir = temp_dir if ( temp_dir is None or writeable ( temp_dir ) ) else tmp_dir ( os.getpid () )  
        self.previous  = None
        
    def __enter__  ( self ) :
        
        self.previous    = tempfile.tempdir        
        if  self.tmp_dir is None or writeable ( self.tmp_dir ) : 
            tempfile.tempdir = self.tmp_dir
            
        return self.__tmp_dir
    
    def __exit__   ( self , *_ ) :
        if self.previous :
            tempfile.tempdir = self.previous
            
    @property
    def tmp_dir ( self ) :
        return self.__tmp_dir 


# =============================================================================
## @class CleanUp
#  Simple (base) class to get temporary files
#  and directories and to remove them at exit
class  CleanUp(object) :
    """Simple (base) class to get temporary files
    and directories and to remove them at exit
    """
    _tmpfiles  = set ()
    _tmpdirs   = set ()
    _protected = set ()
    _failed    = set ()
    
    ## @attention ensure that important attributes are available even before __init__
    def __new__( cls, *args, **kwargs):
        if  python_version.major > 2 : obj = super(CleanUp, cls).__new__( cls )
        else                         : obj = super(CleanUp, cls).__new__( cls , *args , **kwargs )
        ## define the local trash 
        obj.__trash   = set()
        if (3,4) <= python_version : obj.__cleaner = weakref.finalize ( obj , obj._clean_trash_ , obj.__trash )
        else                       : obj.__cleaner = None 
        return obj
 
    @property
    def tmpdir ( self ) :
        """`tmpdir' : return the name of temporary managed directory
        - the managed directory will be cleaned-up and deleted at-exit
        >>> o    = CleanUp()
        >>> tdir = o.tmpdir 
        """
        tdir = CleanUp.tempdir () 
        return tdir
    
    @property
    def tmpdirs  ( self ) :
        """`tmpdirs' - list of currently registered managed temporary directories"""
        return tuple ( self._tmpdirs )
    
    @tmpdirs.setter
    def tmpdirs ( self, values ) :
        if isinstance ( values , str ) : values = [ values ]
        for o in values :
            if o and isinstance ( o ,  str ) :
                self._tmpdirs.add ( o )
                logger.verbose ( 'temporary directory     added %s' % o )
                
    @property 
    def tmpfiles ( self ) :
        """`tempfiles' : list of registered managed temporary files"""
        return list ( self._tmpfiles )
    
    @tmpfiles.setter
    def tmpfiles ( self , other ) :
        if isinstance ( other , str ) : other = [ other ]
        for o in other :
            if o and isinstance ( o , str ) : 
                self._tmpfiles.add ( o )
                logger.verbose ( 'temporary file          added %s' % o )

    @property
    def trash ( self ) :
        """`trash' : local trash (files, directory), to be removed at deletion of the instance"""
        return self.__trash
    
    @classmethod
    def _clean_trash_ ( cls , trash ) :
        """Helper class method to clean the local trash"""
        while trash : cls.remove ( trash.pop () )
        
    ## delete all local trash (files, directory)
    def __del__ ( self ) :
        """Delete all local trash (files, directory)
        """
        if not self.__trash : return
        if ( not self.__cleaner ) or self.__cleaner.detach () :
            self._clean_trash_ ( self.__trash  )             

    # =========================================================================
    @staticmethod
    def tempdir ( suffix = '' , prefix = 'ostap-tmp-dir-' , date = True ) :
        """Get the name of the newly created temporary directory.
        The directory will be cleaned-up and deleted at-exit
        >>> dirname = CleanUp.tempdir() 
        """
        with UseTmpDir ( tmp_dir () ) :
            if date :
                if prefix and prefix.endswith('-') :
                    prefix = "%s%s-"  % ( prefix , start.strftime ( date_format ) )
                else :
                    prefix = "%s-%s-" % ( prefix , start.strftime ( date_format ) )
                    
            td = tempfile.gettempdir()
            if user and not user in td : prefix = '%s%s-' % ( prefix , user )
            
            tmp = tempfile.mkdtemp ( suffix = suffix , prefix = prefix )
            
            CleanUp._tmpdirs.add ( tmp )
            logger.verbose ( 'temporary directory requested %s' % tmp   )
            return tmp        

    # =========================================================================
    @staticmethod
    def get_temp_file ( suffix = '' , prefix = 'ostap-tmp-' , dir = None , date = True ) :
        """ Generate the name for the temporary file.
        - the method should be  avoided in favour of `CleanUp.tempfile`
        >>> fname = CleanUp.get_temp_file () 
        """
        with UseTmpDir ( tmp_dir () ) :
            if date :
                if prefix and prefix.endswith('-') :
                    prefix = "%s%s-"  % ( prefix , start.strftime ( date_format ) )
                else :
                    prefix = "%s-%s-" % ( prefix , start.strftime ( date_format ) )

            td = tempfile.gettempdir()
            if user and not user in td : prefix = '%s%s-' % ( prefix , user )

            with tempfile.NamedTemporaryFile ( suffix = suffix ,
                                               prefix = prefix ,
                                               dir    = dir    , 
                                               delete = False  ) as tfile :
                fname = tfile.name
                
            os.unlink ( fname )            
            assert not os.path.exists ( fname )
            logger.verbose  ( 'temporary file      requested %s' % fname )
            return fname

    # =========================================================================
    @staticmethod
    def tempfile ( suffix = '' , prefix = 'ostap-tmp-' , dir = None , date = True , keep = False ) :
        """ Get the name of the temporary file.
        - The file will be deleted at-exit
        >>> fname = CleanUp.tempfile() 
        """
        fname = CleanUp.get_temp_file ( suffix = suffix , prefix = prefix ,
                                        dir    = dir    , date   = date   )
        assert not os.path.exists  ( fname ) 
        CleanUp._tmpfiles.add ( fname )
        if keep : CleanUp._protected.add ( fname )
        return fname

    # =========================================================================
    @staticmethod
    def protect_file ( fname ) :
        """ Protect the temporary from removal"""        
        if os.path.exists ( fname ) and os.path.isfile ( fname ) :
            CleanUp._protected.add ( fname ) 
            logger.verbose ( 'the file is protected: %s ' % fname )
            
    # =========================================================================
    @staticmethod
    def remove_file ( fname ) :
        """ Remove the (temporary) file
        """
        if os.path.exists ( fname ) and os.path.isfile ( fname ) :

            if fname in CleanUp._protected :
                logger.verbose ( 'do not remove the protected file : %s ' % fname )
                return False

            logger.verbose ( 'remove the file : %s' % fname )
            try    : os.remove ( fname )
            except : pass
            
        if os.path.exists  ( fname ) and os.path.isfile ( fname ) :
            CleanUp._failed.add ( fname ) 
            logger.error   ( 'failed to remove file : %s' %  fname  )
            return False 
        return True 

    # =========================================================================
    @staticmethod
    def remove_dir ( fdir ) :
        """Remove the (temporary) directory
        """
        if os.path.exists ( fdir ) and os.path.isdir ( fdir ) :
            logger.verbose ( 'remove the dir : %s' % fdir  )
            ## 1: collect all files & subdirectories 
            for root, subdirs, files in os.walk ( fdir  , topdown = False ):
                ## 2: remove all files 
                for ff in files : CleanUp.remove_file ( os.path.join ( root , ff  ) ) 
                ## 3: remove all directories 
                for dd in subdirs :
                    dd = os.path.join ( root , dd )
                    logger.verbose ( 'remove subdirectory %s in the directory %s ' % ( dd , fdir ) )
                    try    : os.rmdir   ( dd  )
                    except : pass
                    if os.path.exists ( dd ) and os.path.isdir ( dd )   :
                        CleanUp._failed.add ( dd  ) 
                        logger.error ( 'failed to remove %s in the directory %s ' % ( dd , fdir ) )                        
            ## 4: finally remove the root
            try    :
                os.rmdir ( fdir )
                return True 
            except :
                pass
        if os.path.exists ( fdir ) and os.path.isdir ( fdir ) :
            CleanUp._failed.add ( fdir ) 
            logger.error ( 'failed to  remove : %s' % fdir  )
            return False        
        return True
        
    # =========================================================================
    @staticmethod
    def remove ( fname ) :
        """ Remove temporary object(file or directory) (if any) 
        """
        if os.path.exists ( fname ) :
            if   os.path.isdir  ( fname ) : return CleanUp.remove_dir  ( fname )
            elif os.path.isfile ( fname ) : return CleanUp.remove_file ( fname )
            
# ============================================================================
## @class CUBase
#  A ligth version of the CleanUp to bve used as helper base class
class CUBase(object) :
    """ A ligth version of the CleanUp to be used as helper base class
    """
    # ==========================================================================
    ## expand the (file/dir)name 
    @classmethod
    def name_expand ( cls , name ) :
        """ Expand the (file/dir)name """
        thename  = os.path.expandvars ( name    )
        thename  = os.path.expanduser ( thename )
        thename  = os.path.expandvars ( thename )
        thename  = os.path.expanduser ( thename )
        return     os.path.expandvars ( thename )
    # =========================================================================
    ## remove the file 
    @classmethod 
    def remove_file ( cls , filename ) :
        """ Remove the file """
        return CleanUp.remove_file ( filename ) 
    # =========================================================================
    ## remove the directory 
    @classmethod 
    def remove_dir ( cls , dirname ) :
        """ Remove the directory """
        return CleanUp.remove_dir ( dirname ) 
    # =========================================================================
    ## remove the file or directory 
    @classmethod 
    def remove      ( cls , name ) :
        """ Remove the directory """
        return CleanUp.remove ( name ) 
    
    # =========================================================================
    ## Create the temporary directory
    #  The directory will be cleaned-up and deleted at-exit.
    @classmethod
    def tempdir ( cls , suffix = '-TMP-dir' , prefix = 'ostap-TMP-dir-' , date = True  ) :
        """ Create the temporary directory
        The directory will be cleaned-up and deleted at-exit.
        """
        return CleanUp.tempdir ( suffix = suffix , prefix = prefix, date = date )
    
    # =========================================================================
    ## Ccreate the name for the temporary file 
    #  The file will be deleted at-axit 
    @classmethod
    def tempfile ( cls , suffix = '-tmp' , prefix = 'ostap-TMP-' , dir = None , date = True  ) :
        """ Create the name for the temporary file 
        The file will be deleted at-axit
        """
        return CleanUp.tempfile ( suffix = suffix , prefix = prefix, dir = dir , date = date ) 

# ============================================================================
## Context manager to cleanup PID-dependent directories 
class CleanUpPID(object) :
    """Context manager to cleanup PID-dependent directories
    """
    def __enter__ ( self ) :
        
        pid = os.getpid () 
        self.__piddir = tmp_dir ( pid ) 
        if not pid in base_tmp_pid_dirs : 
            self.__piddir = None
            
    def __exit__ ( self , *_ ) :
        if self.__piddir :
            CleanUp.remove_dir ( self.__piddir )
            
    @property
    def piddir ( self ) :
        """`piddir' : PID-dependent temporary directory"""
        return self.__piddir
    
# =============================================================================
if base_tmp_dir and for_cleanup :
    CleanUp().tmpdirs = base_tmp_dir 
    
# =============================================================================
import atexit
@atexit.register
def _cleanup_ () :
    
    ## 1. clean up the files 
    tmp_files  = CleanUp._tmpfiles
    logger.debug ( 'remove temporary files: %s' % list ( tmp_files ) )
    while tmp_files :
        f = tmp_files.pop()
        CleanUp.remove_file ( f )

    ## 2. clean up  the directories 
    tmp_dirs = CleanUp._tmpdirs
    logger.debug ( 'remove temporary directories: %s' % list ( tmp_dirs ) )
    hdirs    = set()
    hdirs.add ( tempfile.gettempdir() ) 
    while tmp_dirs :
        f = tmp_dirs.pop()
        CleanUp.remove_dir ( f )
        dn = os.path.dirname ( f ) 
        if os.path.exists ( dn ) and os.path.isdir ( dn ) : hdirs.add ( dn )

    ## 3.remove base directories
    global base_tmp_pid_dirs

    for k in base_tmp_pid_dirs :
        d = base_tmp_pid_dirs [ k ]
        CleanUp.remove_dir ( d )
        dn = os.path.dirname ( f ) 
        if os.path.exists ( dn ) and os.path.isdir ( dn ) : hdirs.add ( dn )
        
    base_tmp_pid_dirs = {}

    ## 4. remove base tmp directory 
    if for_cleanup and base_tmp_dir :
        CleanUp.remove_dir ( base_tmp_dir  )

    ## 5. remove old/ancient directories 
    ddirs   = set() 
    for hd in hdirs :
        if not os.path.exists ( hd ) : continue
        if not os.path.isdir  ( hd ) : continue
        for entry in os.listdir ( hd ) :
            if   entry.startswith ( dir_user_prefix ) : e = entry.replace ( dir_user_prefix , '' )
            elif entry.startswith ( dir_prefix      ) : e = entry.replace ( dir_prefix      , '' )
            else                                      : continue  
            match = re_format.search ( e )
            if match :
                date  = match.group ( 'date' )
                pid   = match.group ( 'pid'  )
                
                try :
                    date = datetime.datetime.strptime ( date , r'%Y-%b-%d')
                    pid  = int ( pid ) 
                except :
                    continue

                dd = os.path.join ( hd , entry )
                if os.path.exists ( dd ) and os.path.isdir ( dd ) :
                    mdate = mtime ( dd , subdirs = False )                    
                    if date <= mdate and mdate + time_delta < start :
                        entry = start - mdate , dd 
                        ddirs.add ( entry )
                        
    if ddirs :
        ddirs = list ( ddirs )
        ddirs.sort ()
        rows  = []
        for d , hd in ddirs :
            if os.path.exists ( hd ) and os.path.isdir ( hd ) :
                CleanUp.remove_dir ( hd )
                if not os.path.exists ( hd ) : 
                    row =  '%s'% hd 
                    if     365 <= d.days    : row = row , '%d  year'    % ( d.days    / 365  ) 
                    elif    30 <= d.days    : row = row , '%2d months'  % ( d.days    /  30  ) 
                    elif     7 <= d.days    : row = row , '%1d weeks'   % ( d.days    /   7  )
                    elif     1 <= d.days    : row = row , '%1d days'    % ( d.days           )
                    elif  3600 <= d.seconds : row = row , '%d hours'    % ( d.seconds / 3600 )
                    elif    60 <= d.seconds : row = row , '%d minutes'  % ( d.seconds /   60 )
                    else :                    row = row , '%d seconds'  % ( d.seconds        )                    
                    rows.append ( row ) 
        if rows :
            rows =  [ ( 'directory' , 'age' ) ]  + rows 
            title = 'Remove ancient directories' 
            import ostap.logger.table as T
            table = T.table ( rows , title = title , prefix = "# " , alignment = 'll' )
            logger.info ( '%s:\n%s' % ( title , table ) ) 
            
        
    if CleanUp._protected :
        title = 'Kept temporary files'
        rows = [] 
        for fname in CleanUp._protected :
            if os.path.exists ( fname ) and os.path.isfile ( fname ) :
                row = '%s' % fname
                rows.append ( row )
        if rows : 
            rows.sort()
            rows = [  ('%d' % i , f ) for i, f in enumerate ( rows , start = 1 ) ]
            rows = [ ( '#', title ) ] + rows 
            import ostap.logger.table as T
            table = T.table ( rows , title = title , prefix = "# " , alignment = 'rl' )
            logger.info ( '%s:\n%s' % ( title , table ) ) 

    if CleanUp._failed :
        title = 'Not removed directories/files'
        rows  = []

        def alldirs ( path ) :
            a , b = os.path.split ( path )
            yield a
            while a and b :
                a , b = os.path.split ( a )
                yield a

        pdirs = set() 
        for pp in CleanUp._protected :
            pdirs |= set ( ( p for p in alldirs ( pp ) ) ) 
            
        for fname in CleanUp._failed :
            if os.path.exists ( fname ) and os.path.isdir  ( fname ) and not fname in pdirs :
                row = '%s' % fname
                rows.append ( row )
        for fname in CleanUp._failed :
            if fname in CleanUp._protected : continue 
            if os.path.exists ( fname ) and os.path.isfile ( fname ) :
                row = '%s' % fname
                rows.append ( row )
                
        if rows : 
            rows.sort()
            rows = [  ('%d' % i , f ) for i, f in enumerate ( rows , start = 1 ) ]
            rows = [ ( '#', title ) ] + rows 
            import ostap.logger.table as T
            title = 'Not removed directories/files'
            table = T.table ( rows , title = title  , prefix = "# " , alignment = 'rl' )
            logger.warning ( '%s\n%s' % ( title , table ) ) 

# =============================================================================
## @class TempFile
#  Base class-placeholder for the temporary  file 
class TempFile(object) :
    """Base class-placeholder for the temporary  file 
    """
    def __init__ ( self , suffix = '' , prefix = '' , dir = None ) :
        
        self.__filename  = CleanUp.tempfile ( suffix = suffix ,
                                             prefix = prefix ,
                                             dir    = dir    )
        if (3,4) <= python_version :  
            self.__finalizer = weakref.finalize ( self                   ,
                                                  self._remove_the_file_ ,
                                                  self.__filename        )
        else :
            self.__finalizer = None   
            
    @property
    def filename ( self ) :
        """`filename' : the actual name of temporary file"""
        return self.__filename

    # =========================================================================
    ## context  manager enter: no action 
    def __enter__ ( self      ) :
        """ Context  manager enter: no action"""
        return self

    # =========================================================================
    ## Context manager  exit: delete the file 
    def __exit__  ( self , *_ ) :
        """ Context manager  exit: delete the file
        """
        self._remove_the_file_ ( self.__filename ) 
        self.__filename = ''
        
    ## delete the file 
    def __del__  ( self ) :
        """ Delete the temporary file
        """
        if not self.__filename : return 
        if ( not self.__finalizer ) or self.__finalizer.detach () :
            self._remove_the_file_ ( self.__filename )             
            self.__filename = ''
            
    # =========================================================================
    @classmethod
    def _remove_the_file_ ( cls , name ) :
        if name and os.path.exists ( name ) :
            CleanUp.remove_file ( name )
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    logger.info ( 80*'*' ) 
    # =========================================================================

  
    # =========================================================================    
    logger.info ( 80*'*' ) 
            
# =============================================================================
##                                                                      The END 
# =============================================================================
