#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some simple but useful utilities for 
#   - timing
#   - memory
#   - profiling
#   - ... 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
#  
# =============================================================================
""" Module with some simple but useful utilities for Ostap
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    'with_ipython'         , ## do we run IPython ? 
    'isatty'               , ## is the stream ``isatty'' ?
    'terminal_size'        , ## get the size of terminal cosole
    'writeable'            , ## good writeable directory?
    'whoami'               , ## who am I? 
    'commonpath'           , ## common path(prefix) for list of files
    'copy_file'            , ## copy fiel creating inetremadiet directorys if needed 
    'NoContext'            , ## empty context manager
    'mtime'                , ## last modication/creation time for the path (dir or file)
    'loop_items'           , ## loop over dictionary items 
    'items_loop'           , ## ditto
    ##
    'numcpu'               , ## number of cores/CPUs
    ##
    'typename'             , ## the typename of the object
    'prntrf'               , ## very specific printer of functions 
    ##
    'zip_longest'          , ## itertools.(i)zip.longest
    ##
    'file_size'            , ## get cumulative size of files/directories 
    ##
    'num_fds'              , ## get number of opened file descriptors 
    'get_open_fds'         , ## get list of opened file descriptors
    ##
    'file_info'            , ## very simple infrmation for the file
    ##
    'isfunction'           , ## is it a function (or lambda) ?
    'islambda'             , ## is it a lambda?
    'ismethod'             , ## is it a method?    
    ##
    # =========================================================================
) # ===========================================================================
# =============================================================================
from   ostap.core.meta_info import python_info, whoami  
from   itertools            import zip_longest
import sys, os, datetime, shutil 
# =============================================================================
## is sys.stdout attached to terminal or not  ?
#  @code
#  stream = ...
#  if isatty( stream ) : print('Teminal!')
#  @endcode 
def isatty ( stream = None ) :
    """ Is the stream is attached to terminal?
    >>> stream = ...
    >>> if isatty( stream ) : print('Teminal!')
    >>> if isatty() : print('stdout is terminal!')
    """
    if not stream : stream = sys.stdout
    # ==========================================================================
    if hasattr ( stream , 'isatty' ) : 
        try    : return stream.isatty()
        except : pass
    # ==========================================================================     
    if hasattr ( stream , 'fileno' ) :
        # ======================================================================
        try    : return os.isatty ( stream.fileno () ) 
        except : pass
    ## 
    return False

# ==============================================================================
## does the atream support unicode? 
def has_unicode ( stream = None ) :
    """ Does the stream support unicode?
    """
    if stream is None : stream = sys.stdout
    encoding  = getattr ( stream , 'encoding' , '' )
    if not encoding : return False 
    return encoding.lower().startswith ( 'utf' )
    
# =============================================================================
## helper function that allows to detect running ipython
def with_ipython()  :
    """ Helper function that allows to detect running ipython"""
    try :
        return __IPYTHON__
    except NameError :
        return False

# ============================================================================
fallback      = 80 , 50
# ============================================================================
def terminal_size ( fallback = fallback ) :
    """ Get the terminal console size (use shutil.get_terminal_size)
    >>> width, height = terminal_size () 
    """
    return shutil.get_terminal_size ( fallback ) 

# ===============================================================================
## make directory
#  @code
#  path = ...
#  make_dir( path )
#  @endcode 
def make_dir ( bdir ) :
    """ Make new directory 
    >>> path = ...
    >>> make_dir ( path )
    """
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        if bdir : # ===========================================================
            os.mkdir ( bdir )
            # =================================================================
            if os.path.exists ( bdir ) and os.path.isdir ( bdir ) :
                return os.path.abspath ( bdir)
        # =====================================================================
    except OSError : # ========================================================
        # =====================================================================
        pass
    
    return ''

# =============================================================================
## is this directory writeable?
#  @code
#  my_dir = ...
#  if wrietable ( my_dir ) : ...
#  @endcode
def writeable ( adir ) :
    """ Is this directory is writeable?
    >>> my_dir = ...
    >>> if writeable ( my_dir ) : ...
    """
    if adir and os.path.exists ( adir ) and os.path.isdir ( adir ) :
        # =====================================================================
        import tempfile
        # =====================================================================
        try : # ===============================================================
            # =================================================================
            with tempfile.TemporaryFile ( dir = adir ) : pass
            return True
        except : # ============================================================
            # =================================================================
            return False

    return False    

# =============================================================================
## get a common path(prefix) for list of paths 
commonpath = os.path.commonpath

# =============================================================================
## make directoreis  
make_dirs = os.makedirs

# =============================================================================
## @class NoContext
#  Fake empty context manager to be used as empty placeholder
#  @code
#  with NoContext() :
#  ...  do_something() 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  date 2013-01-12
class NoContext(object) :
    """ Fake (empty) context manager to be used as empty placeholder
    >>> with NoContext() :
    ...         do_something() 
    """
    def __init__  ( self , *args , **kwargs ) : pass
    ## context manager
    def __enter__ ( self         ) : return self 
    ## context manager 
    def __exit__  ( self , *args ) : pass  

# =============================================================================
## loop over dictoribnaty items
def loop_items ( dct ) :
    """ Iterate over the dictionary items
    >>> d = { 'a' : ...   , 'b' : ... , }
    >>> for e in   loop_items ( d ) : print (e) 
    """
    for item in dct.items () : yield item

# =============================================================================
## Iterate over the dictionary items
items_loop = loop_items 

# =============================================================================
## Get the modification time for the path (including subdirectories)
#  @code
#  path = ...
#  mdate = mtime ( path ) ## check file/directory
#  mdate = mtime ( path , subdirs = True ) ## check file/directory including all subdirectoreis
#  @endcode
#  @attention for <cpde>subdirs=True</code> it could be very slow!
def mtime ( path , subdirs = False ) :
    """ Get the last modification time for the path (including subdirectories)
    >>> path = ...
    >>> mdate = mtime ( path ) ## check file/directory
    >>> mdate = mtime ( path , subdirs = True ) ## check file/directory including all subdirectoreis
    - attention: for `subdirs=True` it could be very slow! 
    """
    assert os.path.exists ( path ) , "mtime: the path `%s' does not exist!" % path

    ## get the time of modification/creation 
    _mtime_ = lambda p : max ( os.path.getmtime ( p ) , os.path.getctime ( p ) ) 

    ## own/root  
    tt = _mtime_ ( path ) 
    
    if subdirs and os.path.isdir ( path ) :
        
        for root, dirs, files in os.walk ( path , topdown = False ) :
            
            if os.path.exists ( root ) : tt  = max ( tt , _mtime_ ( root ) ) 

            for d in dirs  :
                entry = os.path.join ( root , d )
                if os.path.exists ( entry ) : tt = max ( tt , _mtime_ ( entry ) )
                
            for f in files :
                entry = os.path.join ( root , f )
                if os.path.exists ( entry ) : tt = max ( tt , _mtime_ ( entry ) ) 

    return datetime.datetime.fromtimestamp ( tt )

# =========================================================================
## copy source file into destination, creating intermediate directories
#  @see https://stackoverflow.com/questions/2793789/create-destination-path-for-shutil-copy-files/49615070 
def copy_file ( source           ,
                destination      ,
                progress = False ) :
    """ Copy source file into destination, creating intermediate directories
    - see https://stackoverflow.com/questions/2793789/create-destination-path-for-shutil-copy-files/49615070 
    """
    assert os.path.exists ( source ) and os.path.isfile ( source ), \
           "copy_file: `source' %s does not exist!" % source 
    
    destination = os.path.abspath  ( destination )    
    destination = os.path.normpath ( destination )
    destination = os.path.realpath ( destination )
    
    if os.path.exists ( destination ) and os.path.isdir ( destination ) :
        destination = os.path.join ( destination , os.path.basename ( source ) )
        
    make_dirs ( os.path.dirname ( destination ) , exist_ok = True )
    
    if not progress : 
        import shutil 
        return shutil.copy2 ( source , destination )
    else :
        from ostap.utils.utils import copy_with_progress
        return copy_with_progress ( source , destination )

# =========================================================================
## Sync/copy source file into destination, creating intermediate directories, using 'rsync -a'
#  @see https://stackoverflow.com/questions/2793789/create-destination-path-for-shutil-copy-files/49615070 
def sync_file ( source              ,
                destination         ,
                progress    = False ) :
    """ Sync/Copy source file into destination, creating intermediate directories
    - see https://stackoverflow.com/questions/2793789/create-destination-path-for-shutil-copy-files/49615070 
    """
    from   ostap.utils.utils import which
    rsync = which ( 'rsync' )
    if not rsync :
        return copy_file ( source = source , destination = destination , progress = progress )
    
    assert os.path.exists ( source ) and os.path.isfile ( source ), \
           "sync_file: `source' %s does not exist!" % source 
    
    destination = os.path.abspath  ( destination )    
    destination = os.path.normpath ( destination )
    destination = os.path.realpath ( destination )
    
    if os.path.exists ( destination ) and os.path.isdir ( destination ) :
        destination = os.path.join ( destination , os.path.basename ( source ) )
        
    make_dirs ( os.path.dirname ( destination ) , exist_ok = True )

    import subprocess, shlex
    
    if progress : command = 'rsync --progress -a %s %s ' % ( source , destination ) 
    else        : command = 'rsync            -a %s %s ' % ( source , destination )
    
    subprocess.check_call ( shlex.split ( command )  )
    
    if not os.path.exists ( destination ) :
        logger.warning ( "copy_files: no expected output '%s'" % nf ) 
        return ''
    
    return destination

# ============================================================================
def __the_function () : pass
__fun_type = type ( __the_function )
# =============================================================================
## very specific printer of objhect
#  - o defiend special print for fumnctins  
def prntrf ( o ) :
    """ very specific printer of objhect
      - o defiend special print for fumnctins  
    """
    if callable ( o ) :
        func_doc = getattr ( o , 'func_doc' , '' )
        if func_doc : return func_doc
        if type ( o ) is __fun_type :                
            return getattr ( o , '__qualname__' , getattr ( o  , '__name__' , 'FUNCTION' ) ) 
    return str ( o )

# =============================================================================
## Get the type name
#  @code
#  obj = ...
#  print ( 'Object type name is %s' % typename ( obj ) ) 
#  @endcode 
def typename ( o ) :
    """ Get the type name
    >>> obj = ...
    >>> print ( 'Object type name is %s' % typename ( obj ) )
    """
    if callable ( o ) :
        to = type ( o ) 
        if to is __fun_type :
            if '<lambda>' == to.__name__ : return 'lambda'
            return getattr ( to , '__qualname__' , getattr ( to , '__name__' ) )
        
    tname = getattr ( o , '__cpp_name__'  ,\
                      getattr ( o , '__qualname__' ,\
                                getattr ( o , '__name__' , '' ) ) )
    if tname : return tname
    to = type ( o )
    return getattr ( to , '__cpp_name__'  ,\
                     getattr ( to , '__qualname__' ,\
                               getattr ( to , '__name__' ) ) )
    

# ==============================================================================
## get the total  size of files/directories
#  @code
#  size = file_size ( 'a.f' , 'b.f'  'c.dir' ) 
#  @endfcode
def file_size ( *files ) :
    """ Get the total  size of files/directories
    >>> size = file_size ( 'a.f' , 'b.f'  'c.dir' ) 
    """
    size = 0
    for name in files :
        if not os.path.exists ( name ) : continue 
        elif   os.path.islink ( name ) : continue 
        elif   os.path.isfile ( name ) : size += os.path.getsize ( name )
        elif   os.path.isdir  ( name ) :
            for dirpath , dirnames , filenames in os.walk ( name ) :
                for f in filenames:
                    fp = os.path.join ( dirpath , f )
                    if not os.path.islink ( fp ):
                        size += os.path.getsize ( fp )
    return size

# =============================================================================
## get all open file descriptors
#  The actual code is copied from http://stackoverflow.com/a/13624412
def get_open_fds():
    """ Get all open file descriptors    
    The actual code is copied from http://stackoverflow.com/a/13624412
    """
    #
    import resource
    import fcntl
    #
    fds = []
    soft , hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    for fd in range ( 0 , soft ) :
        # =====================================================================
        try: # ================================================================
            # =================================================================
            flags = fcntl.fcntl(fd, fcntl.F_GETFD)
            # =================================================================
        except IOError: # =====================================================
            # =================================================================
            continue
        fds.append ( fd )
    return tuple ( fds ) 

# =============================================================================
try : # =======================================================================
    # =========================================================================
    import psutil
    ## get number of opened file descriptors 
    def num_fds () :
        """ Get number of popened file descriptors """        
        p = psutil.Process() 
        return p.num_fds()
    # =========================================================================
except ImportError :
    # =========================================================================
    ## get number of opened file descriptors 
    def num_fds () :
        """ Get number of popened file descriptors"""        
        return len ( get_open_fds () )
    # =========================================================================
    
# =============================================================================
## get the actual file name form file descriptor 
#  The actual code is copied from http://stackoverflow.com/a/13624412
#  @warning: it is likely to be "Linux-only" function
def get_file_names_from_file_number ( fds ) :
    """ Get the actual file name from file descriptor 
    The actual code is copied from http://stackoverflow.com/a/13624412 
    """
    names = []
    for fd in fds:
        names.append(os.readlink('/proc/self/fd/%d' % fd))
    return names

# =============================================================================
## Get number of cores/CPUs
if ( 3 , 13 ) <= python_info : from os import process_cpu_count as cpu_count 
else                         : from os import         cpu_count 

# =============================================================================
## Get number of CPUs     
#  - it uses the function `cpu_count` from `%s` module  
#  - it reads OSTAP_NCPUS envrironment variable 
#  - it checks `General.NCPUS` section of global config
def numcpu () :
    """ Get number of CPUs (non-negative integer number)
    - it uses the function `cpu_count` from `%s` module  
    - it reads OSTAP_NCPUS envrironment variable 
    - it checks `General.NCPUS` section of global config
    """
    nn = cpu_count() 
    ## (1) check the environment variable 
    from   ostap.utils.env   import get_env, OSTAP_NCPUS
    # ========================================================================
    try  : # =================================================================
        nv = int ( get_env ( 'OSTAP_NCPUS' , -1 ) )
        if 1 <= nv : nn = min ( nv , nn )
        # ====================================================================
    except : # ===============================================================
        # ====================================================================
        pass
    # ========================================================================
    ## (2) check the global setting 
    import ostap.core.config as OCC
    ng = OCC.general.getint ( 'ncpus' , -1 )
    if 1 <= ng : nn = min ( nn , ng )
    ## 
    return max ( 1 , nn  ) 

# =============================================================================
## get some file info for the given path
#  - used for multiprocessing of TTree/Tchain
def file_info ( fname ) :
    """ Get some file info for the given path
    - used for multiprocessing of TTree/Tchain    
    """
    p , s , f = fname.partition ( '://' )
    if p and s : return 'Protocol'
    if os.path.exists ( fname ) and os.path.isfile ( fname ) and os.access ( fname , os.R_OK ) :
        s = os.stat ( fname )
        return ( s.st_mode  ,
                 s.st_size  , 
                 s.st_uid   ,
                 s.st_gid   ,
                 s.st_atime ,
                 s.st_mtime ,
                 s.st_ctime ) 
    return 'Invalid'

# =============================================================================
from inspect import ismethod
from types   import FunctionType, LambdaType
# =============================================================================
## is it a function (or lambda) ?
#  @code
#  obj = ...
#  print ( 'function?' , isfunction ( obj ) ) 
#  @endcode 
def isfunction ( func ) :
    """ Is it a function (or lambda) ?
    """
    return isinstance ( func , ( FunctionType , LambdaType ) )
# =============================================================================
## is it a lambda ?
#  @code
#  obj = ...
#  print ( 'lambda?' , islambda ( obj ) ) 
#  @endcode 
def isfunction ( func ) :
    """ Is it a lambda?
    """
    return isinstance ( func , LambdaType )

# =============================================================================
if __name__ == '__main__' :

    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.utils.basic' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
