#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some useful io/file-related utilities 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
# =============================================================================
""" Module with some useful io/file-related utilities
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
__all__     = (
    'make_dir'     , ## make *NEW* directory
    'make_dirs'    , ## recursively nae new directories
    'commonpath'   , ## get a common path(prefix) for the list of paths 
    'mtime'        , ## Get the modification time for the path (including subdirectories
    'copy_file'    , ## copy source file into destination, creating intermediate directories
    'sync_file'    , ## Sync/copy source file into destination, creating intermediate directories, using 'rsync -a'
    'file_size'    , ## get the total  size of files/directories
    'file_info'    , ## get some file info for the given path
    'get_open_fds' , ## get all open file descriptors
    'num_fds'      , ## get number of opened file descriptors
    'get_file_names_from_file_number' , ## get the actual file name from the file descriptors 
)
# =============================================================================
import os, datetime  
# =============================================================================
## make *NEW* directory
#  @code
#  path = ...
#  make_dir( path )
#  @endcode
#  @return name of created directory (on success) 
def make_dir ( bdir ) :
    """ Make *NEW* directory 
    - return name of created directory (on success) 
    >>> path = ...
    >>> make_dir ( path )
    """
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        if bdir : # ===========================================================
            # =================================================================
            os.mkdir ( bdir )
            # =================================================================
            if os.path.exists ( bdir ) and os.path.isdir ( bdir ) :
                return os.path.abspath ( bdir)
        # =====================================================================
    except OSError : # ========================================================
        # =====================================================================
        return ''

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
## get a common path(prefix) for the list of paths 
commonpath = os.path.commonpath

# =============================================================================
## make directories recursively 
make_dirs = os.makedirs
import datetime
import os

# =============================================================================
## Get the modification time for the path (including subdirectories)
#  @code
#  path = ...
#  mdate = mtime ( path ) ## check file/directory
#  mdate = mtime ( path , subdirs = True ) ## check file/directory including all subdirectories
#  @code
#  @attention for <code>subdirs=True</code> it could be slow on huge directories!
def mtime ( path , subdirs = False ) :
    """ Get the last modification time for the path (including subdirectories)

    >>> path = ...
    >>> mdate = mtime ( path ) ## check file/directory
    >>> mdate = mtime ( path , subdirs = True ) ## check file/directory including all subdirectories
    - attention: for `subdirs=True` it could be slow on huge directories!
    """
    assert os.path.exists ( path ) , f"mtime: the path '{path}' does not exist!"

    # =======================================================================
    # Helper function to safely get mtime/ctime from stat result or path
    def _safe_mtime(entry_or_path):
        """ Helper function to safely get mtime/ctime from stat result or path"""
        # ==================================================================
        try: # =============================================================
            # ==============================================================
            if isinstance ( entry_or_path ,  os.DirEntry ) :
                st = entry_or_path.stat(follow_symlinks=False)
            else:
                st = os.stat(entry_or_path, follow_symlinks=False)
            return max(st.st_mtime, st.st_ctime)
            # =============================================================
        except ( OSError , PermissionError ) : # ==========================
            # =============================================================
            return 0.0

    # Own/root timestamp
    tt = _safe_mtime ( path )

    if subdirs and os.path.isdir ( path ):
        # Fast filesystem traversal via os.walk
        for root, dirs, files in os.walk(path, topdown=False):
            # 1. Process files
            for f in files:
                f_path = os.path.join(root, f)
                tt = max(tt, _safe_mtime(f_path))

            # 2. Process subdirectories
            for d in dirs:
                d_path = os.path.join(root, d)
                tt = max(tt, _safe_mtime(d_path))

    return datetime.datetime.fromtimestamp(tt)

# =============================================================================
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
        
    os.makedirs ( os.path.dirname ( destination ) , exist_ok = True )
    
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
    from    ostap.utils.utils import which
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
        
    os.makedirs ( os.path.dirname ( destination ) , exist_ok = True )

    import subprocess, shlex
    
    if progress : command = 'rsync --progress -a %s %s ' % ( source , destination ) 
    else        : command = 'rsync            -a %s %s ' % ( source , destination )
    
    subprocess.check_call ( shlex.split ( command )  )
    
    if not os.path.exists ( destination ) :
        logger.warning ( "copy_files: no expected output '%s'" % nf ) 
        return ''
    
    return destination


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
        """ Get number of opened file descriptors
        """        
        p = psutil.Process() 
        return p.num_fds()
    # =========================================================================
except ImportError :
    # =========================================================================
    ## get number of opened file descriptors 
    def num_fds () :
        """ Get number of popened file descriptors
        """        
        return len ( get_open_fds () )
    # =========================================================================

# =============================================================================
## get the actual file name from the file descriptors 
#  The actual code is copied from http://stackoverflow.com/a/13624412
#  @warning: it is likely to be "Linux-only" function
def get_file_names_from_file_number ( *fds ) :
    """ Get the actual file name from the file descriptor 
    The actual code is copied from http://stackoverflow.com/a/13624412 
    """
    return tuple ( os.readlink ( '/proc/self/fd/%d' % fd ) for fd in fds ) 
    
# =============================================================================
if __name__ == '__main__' :

    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.io.utils' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
