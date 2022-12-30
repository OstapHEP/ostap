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
"""Module with some simple but useful utilities for Ostap
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    'with_ipython'  , ## do we run IPython ? 
    'isatty'        , ## is the stream ``isatty'' ?
    'terminal_size' , ## get the size of terminal cosole
    'writeable'     , ## good writeable directory?
    'whoami'        , ## who am I? 
    'commonpath'    , ## common path(prefix) for list of files
    'copy_file'     , ## copy fiel creating inetremadiet directorys if needed 
    'NoContext'     , ## empty context manager
    ##
    'loop_items'    , ## loop over dictionary items 
    'items_loop'    , ## ditto
    ##
    'has_env'       , ## case-insensitive check for environment variable   
    'get_env'       , ## case-insensitive access to environment variable   
    )
# =============================================================================
import sys, os 
# =============================================================================
## is sys.stdout attached to terminal or not  ?
#  @code
#  stream = ...
#  if isatty( stream ) : print('Teminal!')
#  @endcode 
def isatty ( stream = None ) :
    """Is the stream is attached to terminal?
    >>> stream = ...
    >>> if isatty( stream ) : print('Teminal!')
    >>> if isatty() : print('stdout is terminal!')
    """
    if not stream : stream = sys.stdout
    #
    try :
        return sys.stdout.isatty()
    except : pass 
    #
    try :
        return os.isatty ( sys.stdout.fileno() ) 
    except : pass
    #
    return False

# =============================================================================
## Who am i ?
#  @cdoe
#  print ( 'I am :' % whoami() ) 
#  @endcode 
def whoami () :
    """ Who am i ?
    >>> print ( "I am ", whoami() ) 
    """
    try :
        return os.getlogin()
    except :
        pass
    
    import getpass
    return getpass.getuser() 
    
# =============================================================================
## helper function that allows to detect running ipython
def with_ipython()  :
    """Helper function that allows to detect running ipython"""
    try :
        return __IPYTHON__
    except NameError :
        return False

# ============================================================================
## Get the terminal console size
def terminal_size():
    """Get the terminal console size
    >>> height , width = terminal_size() 
    """
    try :
        import fcntl, termios, struct
        th, tw, hp, wp = struct.unpack(
            'HHHH',fcntl.ioctl(0, termios.TIOCGWINSZ,
                               struct.pack('HHHH', 0, 0, 0, 0)))
        return th , tw  
    except :
        return 50 , 128

    
# ===============================================================================
## make directory
#  @code
#  path = ...
#  make_dir( path )
#  @endcode 
def make_dir ( bdir ) :
    """Make new directory 
    >>> path = ...
    >>> make_dir ( path )
    """
    try :

        if bdir : 
            os.mkdir ( bdir )
            if os.path.exists ( bdir ) and os.path.isdir ( bdir ) :
                return os.path.abspath ( bdir)
            
    except OSError :
        
        pass
    
    return ''

# =============================================================================
## is this directory writeable?
#  @code
#  my_dir = ...
#  if wrietable ( my_dir ) : ...
#  @endcode
def writeable ( adir ) :
    """Is this directory is writeable?
    >>> my_dir = ...
    >>> if writeable ( my_dir ) : ...
    """
    if adir and os.path.exists ( adir ) and os.path.isdir ( adir ) :
        
        import tempfile 
        try :
            with tempfile.TemporaryFile ( dir = adir ) : pass
            return True 
        except :
            return False

    return False    


# =============================================================================
if (3,5) <= sys.version_info :

    ## get a common path(prefix) for list of paths 
    commonpath = os.path.commonpath
    
else :

    ## get a common path(prefix) for list of paths
    # Fix for a  buggy <code>os.path.commonprefix</code>
    #  @see https://www.rosettacode.org/wiki/Find_common_directory_path#Python
    def commonpath ( paths ) :
        """Get a common path(prefix) for list of paths
        Fix a buggy `os.path.commonprefix`
        - https://www.rosettacode.org/wiki/Find_common_directory_path#Python
        """
        return os.path.dirname ( os.path.commonprefix ( paths ) ) 


# =============================================================================
if (3,2) <= sys.version_info :

    make_dirs = os.makedirs
    
else :
    
    def make_dirs ( name , mode = 0o777 , exist_ok = False ):
        """makedirs(path [, mode=0o777])
        
        Super-mkdir; create a leaf directory and all intermediate ones.
        Works like mkdir, except that any intermediate path segment (not
        just the rightmost) will be created if it does not exist.  This is
        recursive.
        
        """
        head, tail = os.path.split(name)
        if not tail:
            head, tail = os.path.split(head)
            
        if head and tail and not os.path.exists(head):
            
            try:
                
                make_dirs ( head , mode , exist_ok = exist_ok ) ## RECURSION!
                
            except OSError as e :
                
                # be happy if someone already created the path
                if e.errno != errno.EEXIST:
                    raise
                
            if tail == os.curdir:           # xxx/newdir/. exists if xxx/newdir exists
                return

        try :
            
            os.mkdir ( name , mode )
            
        except OSError:
            # Cannot rely on checking for EEXIST, since the operating system
            # could give priority to other errors like EACCES or EROFS
            if not exist_ok or not os.path.isdir ( name ) :
                raise

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
    """Fake (empty) context manager to be used as empty placeholder
    >>> with NoContext() :
    ...         do_something() 
    """
    def __init__  ( self , *args , **kwargs ) : pass
    ## context manager
    def __enter__ ( self         ) : return self 
    ## context manager 
    def __exit__  ( self , *args ) : pass  


# =============================================================================
if   ( 3 , 0 ) <= sys.version_info : items_loop = lambda d : d.items     () 
else                               : items_loop = lambda d : d.iteritems ()
# =============================================================================
def loop_items ( d ) :
    """Return  iterator over the dictionary   items
    >>> d = { 'a' : ...   , 'b' : ... , }
    >>> for e in   loop_items ( d ) : print (e) 
    """
    return items_loop ( d ) 

# =============================================================================
## case-insensitive check for existence of the environment variable
#  - case-insensitive
#  - space ignored
#  - underline ignored
#  @code
#  has_env ( 'Ostap_Table_Style' ) 
#  @endcode
def has_env ( variable ) :
    """
    Case-insensitive check for existence of the environment variable
    
    - case-insensitive
    - space ignored
    - underline ignored
    
    >>> has_env ( 'Ostap_Table_Style' )
    
    """
    transform = lambda v : v.replace(' ','').replace('_','').lower()
    new_var   = transform ( variable )
    for key, _ in items_loop ( os.environ ) :
        if transform ( key ) == new_var  : return True
    return False

# =============================================================================
## case-insensitive access for the environment variable
#  - case-insensitive
#  - space ignored
#  - underline ignored
#  @code
#  var = get_env ( 'Ostap_Table_Style' , '' ) 
#  @endcode
#  In ambiguous case warning message is printed and the last value is returned 
def get_env ( variable , default ) :
    """
    Case-insensitive access for the environment variable
    
    - case-insensitive
    - space ignored
    - underline ignored
    
    >>> var = has_env ( 'Ostap_Table_Style' , 'empty' )
    
    In ambiguous case  warning message is printed and the last value is returned 
    """
    transform = lambda v : v.replace(' ','').replace('_','').lower()
    new_var   = transform ( variable )
    found     = [] 
    for key, value in items_loop ( os.environ ) :
        if transform ( key ) == new_var  :
            item = key, value
            found.append  ( item )
            
    if not found : return default

    if 1 < len ( found )  :
        rows = [ ( 'Variable' , 'value' ) ]
        for k, v in found  :
            row = '%s' % k , '%s' % v
            rows.append ( row )
        title  = "'%s' matches" % variable
        import ostap.logger.table as T 
        table = T.table ( rows,  title = title  , prefix = '# ' , alignment = 'll' )
        from ostap.logger.logger import getLogger
        logger = getLogger ( 'ostap.get_env' ) 
        title2 = "Found %s matches for '%s', the last is taken" % ( len ( found ) , variable )  
        logger.warning ( '%s\n%s' %  ( title2 , table ) ) 
    
    return found [ -1] [ 1 ] 

# =========================================================================
## copy source file into destination, creating intermediate directories
#  @see https://stackoverflow.com/questions/2793789/create-destination-path-for-shutil-copy-files/49615070 
def copy_file ( source , destination , progress = False ) :
    """Copy source file into destination, creating intermediate directories
    - see https://stackoverflow.com/questions/2793789/create-destination-path-for-shutil-copy-files/49615070 
    """
    assert os.path.exists ( source ) and os.path.isfile ( source ), \
           "copy_file: ``source'' %s does not exist!" % source 
    
    destination = os.path.abspath  ( destination )    
    destination = os.path.normpath ( destination )
    destination = os.path.realpath ( destination )
    
    if os.path.exists ( destination ) and os.path.isdir ( destination ) :
        destination = os.path.join ( destination , os.path.basename ( source ) )
        
    make_dirs ( os.path.dirname ( destination ) , exist_ok = True)
    
    if not progress : 
        import shutil 
        return shutil.copy2 ( source , destination )
    else :
        from ostap.utils.utils import copy_with_progress
        return copy_with_progress ( source , destination )
    
# =============================================================================
if __name__ == '__main__' :

    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.utils.basic' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
