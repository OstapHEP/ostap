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
    'mtime'         , ## last modication/creation time for the path (dir or file)
    ##
    'loop_items'    , ## loop over dictionary items 
    'items_loop'    , ## ditto
    ##
    'has_env'       , ## case-insensitive check for environment variable   
    'get_env'       , ## case-insensitive access to environment variable
    ##
    'numcpu'        , ## number of cores/CPUs
    ##
    'typename'      , ## the typename of the object
    'prntrf'        , ## very specific printer of functions 
    ##
    'zip_longest'   , ## itertools.(i)zip.longest
    ##
    'var_separators'       , ## separators form the split  
    'split_string_respect' , ## split the string  according to separators 
    'split_string'         , ## split the string  according to separators
    # =========================================================================
) # ===========================================================================
# =============================================================================
import sys, os, datetime
# =============================================================================
if (3,0) <= sys.version_info : from itertools import  zip_longest
else                         : from itertools import izip_longest as zip_longest 
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
    ## 
    if hasattr ( stream , 'isatty' ) : 
        try    : return stream.isatty()
        except : pass
    ##     
    if hasattr ( stream , 'fileno' ) : 
        try    : return os.isatty ( stream.fileno () ) 
        except : pass
    ## 
    return False

# =============================================================================
## Who am I ?
#  @cdoe
#  print ( 'I am :' % whoami() ) 
#  @endcode 
def whoami () :
    """ Who am I ?
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
    """ Helper function that allows to detect running ipython"""
    try :
        return __IPYTHON__
    except NameError :
        return False

# ============================================================================
fallback      = 80 , 50
# ============================================================================
terminal_size = None
# ============================================================================
if ( 3 , 9 ) <= sys.version_info : # =========================================
    # ========================================================================
    try : #===================================================================
        # ====================================================================
        from terminaltables3.terminal_io import terminal_size as tt_terminal_size
        def terminal_size ( fallback = fallback ) :
            """ Get the terminal console size (use terminaltables3)
            >>> width, height = terminal_size () 
            """
            return tt_terminal_size () 
        # ====================================================================
    except ImportError : # ===================================================
        # ====================================================================
        pass
# ============================================================================
if not terminal_size : # =====================================================
    # ========================================================================
    try : #===================================================================
        # ====================================================================
        from terminaltables.terminal_io import terminal_size as tt_terminal_size
        def terminal_size ( fallback = fallback ) :
            """ Get the terminal console size (use terminaltables)
            >>> width, height = terminal_size () 
            """
            return tt_terminal_size () 
        # ====================================================================
    except ImportError : # ===================================================
        # ====================================================================
        pass
# ============================================================================
if not terminal_size and ( 3 , 3 ) <= sys.version_info : # ===================
    # ========================================================================
    import shutil 
    def terminal_size ( fallback = fallback ) :
        """ Get the terminal console size (use shutil.get_terminal_size)
        >>> width, height = terminal_size () 
        """
        return shutil.get_terminal_size ( fallback ) 
    # ========================================================================
elif not terminal_size  : # ====================================================
    # ========================================================================
    ## Get the terminal console size
    def terminal_size ( fallback = fallback ):
        """ Get the terminal console size (use local version) 
            >>> width, height = terminal_size () 
            """
        # ====================================================================
        try : # ==============================================================
            # ================================================================
            import fcntl, termios, struct
            th, tw, hp, wp = struct.unpack(
                'HHHH',fcntl.ioctl(0, termios.TIOCGWINSZ,
                                   struct.pack('HHHH', 0, 0, 0, 0)))
            return tw  , th
            # ================================================================
        except : # ============================================================
            # =================================================================
            return fallback
        # =================================================================

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
    """ Is this directory is writeable?
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

    # ========================================================================
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
    # =========================================================================    
    make_dirs = os.makedirs
    # =========================================================================        
else :
    # =========================================================================    
    def make_dirs ( name , mode = 0o777 , exist_ok = False ):
        """ makedirs(path [, mode=0o777])
        
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
if   ( 3 , 0 ) <= sys.version_info : # ========================================
    # =========================================================================
    def loop_items ( dct ) :
        """ Iterate over the dictionary items
        >>> d = { 'a' : ...   , 'b' : ... , }
        >>> for e in   loop_items ( d ) : print (e) 
        """
        for item in dct.items () : yield item
    # =========================================================================
else : # ======================================================================
    # =========================================================================
    def loop_items ( dct ) :
        """ Iterate over the dictionary items
        >>> d = { 'a' : ...   , 'b' : ... , }
        >>> for e in   loop_items ( d ) : print (e) 
        """
        for item in dct.iteritems () : yield item
    # =========================================================================
        
# =============================================================================
## Iterate over the dictionary items
items_loop = loop_items 

# =============================================================================
## case-insensitive check for existence of the environment variable
#  - case-insensitive
#  - space ignored
#  - underline ignored
#  @code
#  has_env ( 'Ostap_Table_Style' ) 
#  @endcode
def has_env ( variable ) :
    """ Case-insensitive check for existence of the environment variable
    
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
def get_env ( variable , default , silent = False ) :
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

    if 1 < len ( found )  and not silent :
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
    
# =============================================================================
if ( 3 , 4 ) <= sys.version_info : 
    # =========================================================================
    ## Get number of cores/CPUs
    from os              import cpu_count as _numcpu
    # =========================================================================
else :
    # =========================================================================
    ## Get number of cores/CPUs    
    from multiprocessing import cpu_count as _numcpu
    # =========================================================================

# =============================================================================
## defalt separators for the string expressions
var_separators = ',:;'
## rx_separators  = re.compile ( r'[,:;]\s*(?![^()]*\))' )
## rx_separators  = re.compile ( '[ ,:;](?!(?:[^(]*\([^)]*\))*[^()]*\))')
# =============================================================================
## mark for double columns: double column is a special for C++ namespaces 
dc_mark = '_SSPPLLIITT_'
# =============================================================================
## split string using separators and respecting the (),[] and {} groups.
#  - group can be nested
def split_string_respect  ( text , separators = var_separators , strip = True ) :
    """ Split string using separators and respecting the (),[] and {} groups.
    - groups can be nested
    """
    protected = False 
    if ':' in separators and '::' in text :
        text      = text.replace ( '::' , dc_mark ) 
        protected = True 
    
    flag1  = 0
    flag2  = 0
    flag3  = 0
    item   = ''
    items  = []
    for c in text:
        if   c == '(' : flag1 += 1
        elif c == ')' : flag1 -= 1
        elif c == '[' : flag2 += 1
        elif c == ']' : flag2 -= 1
        elif c == '{' : flag3 += 1
        elif c == '}' : flag3 -= 1
        elif 0 == flag1 and 0 == flag2 and 0 == flag3 and c in separators :
            items .append ( item )
            item = ''
            continue
        item += c
        
    if item : items.append ( item  )

    if protected :
        nlst = []
        for item in items :
            if dc_mark in item : nlst.append ( item.replace ( dc_mark , '::' ) )
            else               : nlst.appenf ( item ) 
        items = nlst 
                              
    ## strip items if required 
    if strip : items = [ item.strip() for item in items ] 
    ## remove empty items 
    return tuple ( item for item in items if item  )

# =============================================================================
## split string using separators:
#  @code
#  split_string ( ' a b cde,fg;jq', ',;:' )
#  @endcode
def split_string ( line                            ,
                   separators     = var_separators ,
                   strip          = False          ,
                   respect_groups = False          ) :
    """ Split the string using separators
    >>> split_string ( ' a b cde,fg;jq', ',;:' )
    """
    if respect_groups :
        return split_string_respect ( line                    ,
                                      separators = separators ,
                                      strip      = strip      )
    ##
    protected = False 
    if ':' in separators and '::' in line :
        line      = line.replace ( '::' , dc_mark ) 
        protected = True 
    
    items = [ line ]
    for s in separators :
        result = []
        for item in items :
            if s in item : result += item.split ( s )
            else         : result.append ( item ) 
        items = result

    if protected :
        nlst = []
        for item in items :
            if dc_mark in item : nlst.append ( item.replace ( dc_mark , '::' ) )
            else               : nlst.appenf ( item ) 
        items = nlst 
        
    ## strip items if required 
    if strip : items = [ i.strip() for i in items ] 
    ## remove empty items 
    return tuple ( item for item in items if item )
    
# =============================================================================
## Get number of CPUs     
def numcpu () :
    """ Get number of CPUs (non-negative integer number)
    - it uses the function `cpu_count` from `%s` module  
    """
    nc = _numcpu ()
    return nc if nc and 0 < nc else 0

numcpu.__doc__ = numcpu.__doc__ % _numcpu.__module__ \
    +  '\n' + _numcpu.__doc__

# =============================================================================
if __name__ == '__main__' :

    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.utils.basic' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
