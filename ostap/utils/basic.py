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
    ## 
    'with_ipython'         , ## do we run IPython ?
    'interactive'          , ## interactive processing?
    ##
    'isatty'               , ## is the stream ``isatty'' ?
    'terminal_size'        , ## get the size of terminal cosole
    'whoami'               , ## who am I?
    ##
    'NoContext'            , ## empty context manager
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
    'isfunction'           , ## is it a function (or lambda) ?
    'islambda'             , ## is it a lambda?
    'ismethod'             , ## is it a method?    
    ## 
    'counted'              , ## helper to count function calls 
    'memoize'              , ## lightweigth cache
    ##
    # =========================================================================
) # ===========================================================================
# =============================================================================
from   ostap.core.meta_info import python_info, whoami  
from   itertools            import zip_longest
import sys, os, datetime, shutil, functools
# =============================================================================
## Interactive processing ?
#  @see https://stackoverflow.com/questions/2356399/tell-if-python-is-in-interactive-mode
def interactive () :
    """ Interactive processing ?
    - see https://stackoverflow.com/questions/2356399/tell-if-python-is-in-interactive-mode
    """
    return hasattr ( sys , 'ps1' )
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
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        return __IPYTHON__
        # =====================================================================
    except NameError : # ======================================================
        # =====================================================================
        return False

# ============================================================================
fallback      = 80 , 50
# ============================================================================
def terminal_size ( fallback = fallback ) :
    """ Get the terminal console size (use shutil.get_terminal_size)
    >>> width, height = terminal_size () 
    """
    return shutil.get_terminal_size ( fallback ) 

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
## make directories
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
## loop over dictionary items
def loop_items ( dct ) :
    """ Iterate over the dictionary items
    >>> d = { 'a' : ...   , 'b' : ... , }
    >>> for e in   loop_items ( d ) : print (e) 
    """
    for item in dct.items () : yield item

# =============================================================================
## Iterate over the dictionary items
items_loop = loop_items 

# ============================================================================
def __the_function () : pass
__fun_type = type ( __the_function )
# =============================================================================
## very specific printer of object
#  - o defiend special print for functins  
def prntrf ( o ) :
    """ very specific printer of object
      - o defined special print for functins  
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
## Get number of cores/CPUs
if ( 3 , 13 ) <= python_info : from os import process_cpu_count as cpu_count 
else                         : from os import         cpu_count 
# =============================================================================
## Get number of CPUs     
#  - it uses the function `cpu_count` from `%s` module  
#  - it reads OSTAP_NCPUS environment variable 
#  - it checks `General.NCPUS` setting for global config
def numcpu () :
    """ Get number of CPUs (non-negative integer number)
    - it uses the function `cpu_count` from `%s` module  
    - it reads OSTAP_NCPUS envrironment variable 
    - it checks `General.NCPUS` section of global config
    """
    # ========================================================================
    ## (1) check the system 
    nn = cpu_count () 
    # ========================================================================
    ## (2) Check the global Ostap configuration: 
    import ostap.core.config as config 
    nc = config.ncpus
    if 1 <= nc : nn = min ( nn , nc )
    ## 
    return max ( 1 , nn  ) 

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
def islambda  ( func ) :
    """ Is it a lambda?
    """
    return isinstance ( func , LambdaType )

# =============================================================================
## create 'counted' function to know number of function calls
#  @code
#  fun = ...
#  func = counted ( fun ) ## use as function
#
#  # alternatively use it as decorator:
#  @counted
#  def fun2 ( ...  ) : return ...
#  @endcode
def counted ( fun ):
    """ Create 'counted' function to know number of function calls

    Example
    -------

    >>> fun = ...
    >>> func = counted ( fun ) ## use as function

    >>> @counted
    >>> def fun2 ( ...  ) : return ...
    """
    def wrapped ( *fargs, **kwargs ):
        wrapped.calls += 1
        return fun ( *fargs , **kwargs )
    wrapped.calls = 0
    return wrapped


# =============================================================================
if   ( 3 , 9 ) <= python_info : # =============================================
    # =========================================================================
    memoize = functools.cache
    # =========================================================================
else : # ======================================================================
    # =========================================================================
    ## Simple lightweight unbounded cache
    def memoize ( user_function ):
        """ Simple lightweight unbounded cache
        - see `functools.lru_cache`
        """
        return functools.lru_cache(maxsize=None)(user_function)

# =============================================================================
## Print/format warning message in one line
#  @see warnings.WarnigMessage 
def wm_print ( wm , with_category = True ) : 
    """ Print/format warning message in one line
    - see warnings.WarningMessage 
    """
    if with_category:
        msg = "%%s,category=%s,file=%%s,line#=%%d"
        msg = msg % wm._categrory_name
    else : msg = "%s,file=%s,line#=%d"
    ## 
    fname = wm.filename 
    if 60 < len ( fname ) : fname = os.path.basename ( fname ) 
    msg = msg % ( wm.message , fname , wm.lineno )
    msg = msg.replace ( '.,' , ',' )
    msg = msg.replace ( '\n' , ' ' )
    while '  ' in msg : msg = msg.replace ( '  ' , ' ' )    
    if wm.line : msg += '%s' % wm.line    
    return msg

# =============================================================================
if __name__ == '__main__' :

    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.utils.basic' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
