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
    ##
    'numcpu'               , ## number of cores/CPUs
    ##
    'prntrf'               , ## very specific printer of functions 
    ##
    'zip_longest'          , ## itertools.(i)zip.longest
    ##
    'counted'              , ## helper to count function calls 
    'memoize'              , ## lightweight cache
    ##
    # =========================================================================
) # ===========================================================================
# =============================================================================
from   ostap.core.meta_info import python_info, whoami
from   ostap.utils.core     import typename, isfunction
from   itertools            import zip_longest
import multiprocessing      as     mp_
import sys, os, datetime, shutil, functools, cppyy, threading 
# =============================================================================
## current process 
current_process = mp_.current_process()
## main process ? 
main_process    = current_process and current_process.name == 'MainProcess'  
# ==============================================================================
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

# ============================================================================
## very specific printer of object
#  - o defined special print for functinos   
def prntrf ( o ) :
    """ very specific printer of object
      - o defined special print for functions  
    """
    if callable ( o ) :
        func_doc = getattr ( o , '__doc__' , '' )
        if func_doc : return func_doc
        if isfunction ( o ) :            
            return getattr ( o , '__qualname__' , '' ) or getattr ( o  , '__name__' , 'FUNCTION' ) 
    return str ( o )
     
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

# ==============================================================================
## key-word arguments related to number of cores/threads to be used 
njobs_kwords = ( 'num_threads'  , 'num_thread'  ,
                 'numthreads'   , 'numthread'   ,
                 'n_threads'    , 'n_thread'    ,
                 'nthreads'     , 'nthread'     ,
                 'num_jobs'     , 'num_job'     ,
                 'numjobs'      , 'numjob'      ,
                 'n_jobs'       , 'n_job'       ,
                 'njobs'        , 'njob'        ,
                 'thread_count' , 'threadcount' )
# ==============================================================================
## get the value of "n_jobs/num_threads/thread_count" parameter 
def num_jobs ( params , defval = -2 ) :
    """ Get the value of "n_jobs/num_threads/thread_count" parameter
    """
    nj = params.pop ( njobs_kwords [ 0 ] , defval   )
    for kw in njobs_kwords[1:] : nj = params.pop ( kw , nj )
    return nj

# ==============================================================================
## allow parallel run ? 
#  - check "OMP_NUM_THREADS"
#  - check "MKL_NUM_THREADS"
#  - check "OPENBLAS_NUM_THREADS"
def run_parallel ( parallel ) : 
    """ Allow parallel run
    - check "OMP_NUM_THREADS"
    - check "MKL_NUM_THREADS"
    - check "OPENBLAS_NUM_THREADS"
    """
    ##
    if not parallel : return False
    ## 
    if   '1' != os.environ.get ( "OMP_NUM_THREADS"      , "" ).strip () : return False 
    elif '1' != os.environ.get ( "MKL_NUM_THREADS"      , "" ).strip () : return False 
    elif '1' != os.environ.get ( "OPENBLAS_NUM_THREADS" , "" ).strip () : return False
    ## 
    return True 

# =============================================================================
if __name__ == '__main__' :

    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.utils.basic' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
