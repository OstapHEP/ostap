#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some simple but useful utilities for 
#   - timing
#   - memory
#   - profiling
#   - ... 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
#  
# =============================================================================
"""Module with some simple but useful utilities for
- timing
- memory
- profiling
- etc
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    #
    'virtualMemory'      , ## context manager to count virtual memory increase 
    'memory'             , ## ditto
    'timing'             , ## context manager to count time 
    'timer'              , ## ditto
    'profiler'           , ## context manager to perform profiling
    #
    'Profiler'           , ## context manager to perform profiling 
    ##
    'takeIt'             , ## take and later delete ...
    'isatty'             , ## is the stream ``isatty'' ?
    'with_ipython'       , ## do we run IPython?
    ##
    'batch'              , ## context manager to keep/force certain ROOT ``batch''-mode
    ##
    'keepCanvas'         , ## context manager to keep the current ROOT canvas
    'invisibleCanvas'    , ## context manager to use the invisible current ROOT canvas
    ##
    'keepArgs'           , ## context manager to keep sys.argv
    ##
    'keepCWD'            , ## context manager to keep current working directory 
    ##
    'implicitMT'         , ## context manager to enable/disable implicit MT in ROOT 
    ##
    'Batch'              , ## context manager to keep  ROOT ``batch''-mode
    ##
    'KeepCanvas'         , ## context manager to keep the current ROOT canvas
    'InvisibleCanvas'    , ## context manager to use the invisible current ROOT canvas
    ##
    'KeepArgs'           , ## context manager to keep sys.argv
    ##
    'Wait'               , ## conitext manager to wait soem tiem bvefore and/or after action
    ## 
    'wait'               , ## conitext manager to wait soem tiem bvefore and/or after action 
    ##
    'ImplicitMT'         , ## context manager to enable/disable implicit MT in ROOT 
    ##
    'counted'            , ## decorator to create 'counted'-function
    ##
    'cmd_exists'         , ## check the existence of the certain command/executable
    ##
    'which'              , ## which command (from shutil)
    ##
    'gen_password'       , ## generate password/secret
    ##
    'vrange'             , ## helper loop over values between xmin and xmax
    ## 
    'log_range'          , ## helper loop over values between xmin and xmax in log
    ## 
    'lrange'             , ## helper loop over values between xmin and xmax in log
    'prange'             , ## helper loop over values between xmin and xmax in log
    'prange'             , ## helper loop over values between xmin and xmax in log
    ##
    'crange'             , ## helper loop over values between xmin and xmax using Chebyshev nodes 
    ##
    'split_range'        , ## helper generator to split large range into smaller chunks
    'split_n_range'      , ## helper generator to split large range into smaller chunks of approximately same size 
    ##
    'chunked'            , ## break *iterable* into chunks of length *n*:
    'divide'             , ## divide the elements from *iterable* into *n* parts
    'grouper'            , ## collect data into fixed-length chunks or blocks"
    ##
    'make_iterable'      , ## create infinite or finite iterable 
    ##
    'checksum_files'     , ## get SHA512 sum for sequence of files
    ##
    'balanced'           , ## Simple utility to check balanced parenthesis/brackets, etc...
    ##
    'random_name'        , ## get some random name
    'short_hash_name'    , ## get some short hash name
    ##
    'choices'            , ## `random.choices` function
    ## 
    'memoize'            , ## Simple lightweight unbounded cache
    'absproperty'        , ## abstract property decorator
    'classprop'          , ## class property decorator
    'numcalls'           , ## decoratro for #ncalls  
    ##
    'hadd'               , ## merge ROOT files using command `hadd`
    'num_fds'            , ## get number of opened file descriptors 
    'get_open_fds'       , ## get list of opened file descriptors
    ##
    'slow'               , ## "slow" looping with delays at each step
    ##
    'CallThem'           , ## convert sequence of callables into singel callable
    'has_symbol'         , ## Has any of the symbols?
    'is_formula'         , ## Is this string expression represent math formula?
    ##
    'RandomSeed'         , ## context manager to set/keep seed for the PYTOHN random generator
    'randomSeed'         , ## context manager to set/keep seed for the PYTHON random generator
    'random_seed'        , ## context manager to set/keep seed for the PYTHON random generator
    ## 
    'RootRandomSeed'     , ## context manager to set/keep seed for the ROOT random generator
    'rootRandomSeed'     , ## context manager to set/keep seed for the ROOT random generator
    'root_random_seed'   , ## context manager to set/keep seed for the ROOT random
    ## 
    'RooRandomSeed'      , ## context manager to set/keep seed for the RooFit random generator
    'rooRandomSeed'      , ## context manager to set/keep seed for the RooFit random generator
    'roo_random_seed'    , ## context manager to set/keep seed for the RooFit random generator
    ## 
    )

# =============================================================================
from   builtins               import range
from   itertools              import repeat, chain, islice 
from   sys                    import version_info  as python_version 
## timing stuff
from   ostap.utils.timing     import timing, timer
## other useful stuff 
from   ostap.utils.basic      import isatty, with_ipython, NoContext
from   ostap.core.ostap_types import ( integer_types  , num_types ,
                                       string_types   ,
                                       dictlike_types , listlike_types )
## ... and more useful stuff 
from   ostap.utils.memory     import memory, virtualMemory, Memory 
import ROOT, time, os , sys, math, time, functools, abc, array, random, datetime  ## attention here!!
# =============================================================================
try :
    from string import ascii_letters, digits 
except ImportError :
    from string import letters as ascii_letters
    from string import digits
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.utils' )
else                       : logger = getLogger( __name__            )
del getLogger
# =============================================================================
## symbols for name generation 
all_symbols = ascii_letters + digits 
# =============================================================================
## @class Profiler
#  Very simple profiler, based on cProfile module
#  @see https://docs.python.org/2/library/profile.html
#  @code
#  with profiler() :
#      ...  some code here ... 
#  with profiler('output.file') :
#      ...  some code here ... 
#  @endcode
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2016-07-25                     
class Profiler(object) :
    """Very simple profiler, based on cProfile module
    - see https://docs.python.org/2/library/profile.html
    
    with profiler() :
    #
    # ...  some code here ...
    #

    with profiler( 'output.file' ) :
    #
    # ...  some code here ...
    # 
    """
    def __init__  ( self , fname = '' )  :
        self.fname = fname
        
    ## enter the context
    def __enter__ ( self ) :
        import cProfile as profile
        self._profile = profile.Profile()
        self._profile.enable()
        return self
    
    ## exit the context
    def __exit__ ( self , *_ ) :
        ## end of profiling 
        self._profile.disable()
        
        import pstats
        if self.fname :
            try :
                with open ( self.fname , 'w' ) as out :
                    stat = pstats.Stats( self._profile , stream = out ).sort_stats( 'cumulative' )
                    stat.print_stats()
                del self._profile 
                return 
            except : pass
            
        ## show on screen 
        stat = pstats.Stats( self._profile ).sort_stats( 'cumulative' )
        stat.print_stats()
        del self._profile 
                
# =============================================================================
## Very simple profiler, based on cProfile module
#  @see https://docs.python.org/2/library/profile.html
#  @code
#  with profiler() :
#      ...  some code here ... 
#  @endcode
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2016-07-25                     
def profiler( name = '' ) :
    """Very simple profiler, based on cProfile module
    - see https://docs.python.org/2/library/profile.html
    
    with profiler() :
    #
    # ...  some code here ...
    # 
    """
    return Profiler ( name )
            
# =============================================================================
## @class TakeIt
#  Take some object, keep it and delete at the exit
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  date 2014-08-03    
class TakeIt(object):
    """Take some object, keep it and delete at the exit
    
    >>> ds = dataset.reduce('pt>1')
    >>> with takeIt ( ds ) :
    ...
    
    """
    def __init__  ( self , other ) :
        self.other = other
        
    def __enter__ ( self ) :
        ROOT.SetOwnership ( self.other , True )
        return self.other
    
    def __exit__  ( self , *_ ) :

        o = self.other

        ## delete it! 
        del self.other
        
        if o and hasattr ( o , 'reset'  ) : o.reset  ()
        if o and hasattr ( o , 'Reset'  ) : o.Reset  ()
        if o and hasattr ( o , 'Delete' ) : o.Delete ()
        
        if o : del o
                        
    
# =============================================================================
## Take some object, keep it and delete at the exit
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  date 2014-08-03    
def takeIt (  other ):
    """Take some object, keep it and delete at the exit
    >>> ds = dataset.reduce('pt>1')
    >>> with takeIt ( ds ) :
    ...    
    """
    return TakeIt ( other ) 

    
# =============================================================================
## get all open file descriptors
#  The actual code is copied from http://stackoverflow.com/a/13624412
def get_open_fds():
    """Get all open file descriptors    
    The actual code is copied from http://stackoverflow.com/a/13624412
    """
    #
    import resource
    import fcntl
    #
    fds = []
    soft , hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    for fd in range ( 0 , soft ) :
        try:
            flags = fcntl.fcntl(fd, fcntl.F_GETFD)
        except IOError:
            continue
        fds.append ( fd )
    return tuple ( fds ) 


# =============================================================================
try :
    # =========================================================================
    import psutil
    ## get number of opened file descriptors 
    def num_fds () :
        "Get number of popened file descriptors"        
        p = psutil.Process() 
        return p.num_fds()
    # =========================================================================
except ImportError :
    # =========================================================================
    ## get number of opened file descriptors 
    def num_fds () :
        "Get number of popened file descriptors"        
        return len ( get_open_fds () )
    # =========================================================================
    
# =============================================================================
## get the actual file name form file descriptor 
#  The actual code is copied from http://stackoverflow.com/a/13624412
#  @warning: it is likely to be "Linux-only" function
def get_file_names_from_file_number(fds):
    """Get the actual file name from file descriptor 
    The actual code is copied from http://stackoverflow.com/a/13624412 
    """
    names = []
    for fd in fds:
        names.append(os.readlink('/proc/self/fd/%d' % fd))
    return names

# =============================================================================
## context manager to keep ROOT ``batch'' state
#  @code
#  with Batch() :
#  ... do something here 
#  @endcode 
class Batch(object) :
    """Context manager to keep ROOT ``batch'' state
    >>> with Batch() :
    ... do something here 
    """
    def __init__  ( self , batch = True ) :
        self.__batch = batch 
    ## contex manahger: ENTER
    def __enter__ ( self ) :
        import ROOT
        groot = ROOT.ROOT.GetROOT()
        self.old_state = groot.IsBatch()
        if self.old_state != self.__batch : groot.SetBatch ( self.__batch ) 
        return self
    ## contex manager: EXIT
    def __exit__  ( self , *_ ) :
        import ROOT
        groot = ROOT.ROOT.GetROOT()
        if self.old_state != groot.IsBatch() : groot.SetBatch( self.old_state ) 

# =============================================================================
## context manager to keep ROOT ``batch'' state
#  @code
#  with batch() :
#  ... do something here 
#  @endcode 
def batch( batch = True ) :
    """Context manager to keep ROOT ``batch'' state
    >>> with batch() :
    ... do something here 
    """
    return Batch ( batch )

# =============================================================================
## context manager to keep the current working directory
#  @code
#  with KeepCWD ( new_dir ) :
#    ....
#  @endcode
#  - No action if no directory is specified 
class KeepCWD(object) :
    """context manager to keep the current working directory
    >>> with KeepCWD( new_dir ) :
    ...
    - No action if no directory is specified 
    """
    def __init__ ( self , new_dir = '' ) :
        
        self.__old_dir = os.getcwd ()
        self.__new_dir = new_dir

    ## ENTER : context mamager 
    def __enter__ (  self ) :
        
        self.__old_dir = os.getcwd()
        
        if   self.new_dir :
            os.chdir ( self.new_dir )
            
        return self
        
    ## EXIT : context mamager 
    def __exit__ ( self , *_ ) :
        
        if os.path.exists ( self.old_dir ) and os.path.isdir ( self.old_dir ) :
            os.chdir ( self.old_dir )
            
    @property
    def old_dir ( self ) :
        """``old_dir'' : old working directory"""
        return self.__old_dir

    @property
    def new_dir ( self ) :
        """``new_dir'' : new current working directory"""
        return self.__new_dir 

    
# =============================================================================
## context manager to keep the current working directory
#  @code
#  with keepCWD ( new_dir ) :
#    ....
#  @endcode 
#  - No action if no directory is specified 
def keepCWD ( new_dir = '' ) :
    """Context manager to keep the current working directory
    >>> with keepCWD( new_dir ) :
    ...
    - No action if no directory is specified 
    """
    return KeepCWD ( new_dir ) 


# =============================================================================
## context manager that invokes <code>time.sleep</code> before and after action
#  @code
#  with Wait ( after = 5 , before = 0 ) :
#  ...
#  @endcode
class Wait(object):
    """Context manager that invokes <code>time.sleep</code> before and after action
    >>> with Wait ( after = 5 , before = 0 ) :
    >>> ...
    """
    def __init__ ( self , after = 0 , before = 0 ) :
        self.__after  = after
        self.__before = before 

    def __enter__ ( self ) :
        if 0 < self.__before : time.sleep  ( self.__before ) 
    def __exit__ ( self , *_ ) :
        if 0 < self.__after  : time.sleep  ( self.__after  ) 
    @property
    def before ( self ) :
        """``before'': wait some time before the action"""
        return self.__before    
    @property
    def after  ( self ) :
        """``after'': wait some time after the action"""
        return self.__after

# =============================================================================
## context manager that invokes <code>time.sleep</code> before and after action
#  @code
#  with wait ( after = 5 , before = 0 ) :
#  ...
#  @endcode
def wait ( after = 0 , before = 0 ) :
    """Context manager that invokes <code>time.sleep</code> before and after action
    >>> with wait ( after = 5 , before = 0 ) :
    >>> ...
    """    
    return Wait ( after = after , before = before )


# =============================================================================
## "slow" looping over some iterable with delay for each step
# @code
# for i in slow ( range ( 5 ) , wait = 0.1 ) :
# ... 
# @endcode 
def slow ( iterable , wait = 0 ) :
    """ `slow' looping over some iterable with delay for each eatep 
    >>> for i in slow ( range ( 5 ) , wait = 0.1 ) :
    >>> ... 
    """
    if isinstance ( iterable , int ) and 0 <= iterable :
        iterable = range ( iterable )
        
    for r in iterable :
        if 0 < wait : time.sleep ( wait )
        yield r 
        
    
# =============================================================================
## @class KeepCanvas
#  helper class to keep the current canvas
#  @code
#  with KeepCanvas() :
#  ... do something here 
#  @endcode 
class KeepCanvas(Wait) :
    """Helper class to keep the current canvas
    >>> with KeepCanvas() :
    ... do something here 
    """
    def __init__ ( self , wait = 0 ) :
        Wait.__init__ ( self , after = wait ) 
        self.__old_canvas  = None 
    def __enter__ ( self ) :
        Wait.__enter__ ( self ) 
        import ROOT
        cnv  = ROOT.gPad.GetCanvas() if ROOT.gPad else None 
        self.__old_canvas = cnv if cnv else None 
    def __exit__  ( self , *_ ) :
        Wait.__exit__ ( self , *_ )         
        if self.__old_canvas:
            self.__old_canvas.cd()
        self.__old_canvas = None             
    @property
    def old_canvas ( self ) :
        """``old_canvas'': canvas to be preserved"""
        return self.__old_canvas
    
# =============================================================================
#  Keep the current canvas
#  @code
#  with keepCanvas() :
#  ... do something here 
#  @endcode
def keepCanvas() :
    """Keep the current canvas
    >>> with keepCanvas() :
    ... do something here
    """
    return KeepCanvas()


# =============================================================================
## @class InvisibleCanvas
#  Use context ``invisible canvas''
#  @code
#  with InvisibleCanvas() :
#  ... do somehing here 
#  @endcode
class InvisibleCanvas(KeepCanvas) :
    """Use context ``invisible canvas''
    >>> with InvisibleCanvas() :
    ... do something here 
    """
    ## context manager: ENTER 
    def __enter__ ( self ) :
        ## start from keeping the current canvas 
        KeepCanvas.__enter__ ( self )
        ## create new canvas in batch mode 
        with Batch( True ) : 
            import ROOT 
            self.batch_canvas = ROOT.TCanvas()
            self.batch_canvas.cd ()
            return self.canvas

    ## context manager: EXIT
    def __exit__ ( self , *_ ) :
        if self.batch_canvas :
            self.batch_canvas.Close() 
            del self.batch_canvas             
        KeepCanvas.__exit__ ( self , *_ )

# =============================================================================
## Use context ``invisible canvas''
#  @code
#  with invisibleCanvas() :
#  ... do something here 
#  @endcode
def invisibleCanvas() :
    """ Use context ``invisible canvas''
    >>> with invisibleCanvas() :
    ... do something here 
    """
    return InvisibleCanvas() 


# =============================================================================
## @class KeepArgs
#  context manager to keep/preserve sys.argv
#  @code
#  with KeepArgs() :
#    ...  
#  @endcode 
class KeepArgs(object) :
    """Context manager to keep/preserve sys.argv
    >>> with KeepArgs() :
    ...  
    """
    ## context manager  ENTER 
    def __enter__ ( self ) :
        import sys, copy
        self._args = copy.deepcopy( sys.argv )
        return self
    ## context manager  EXIT
    def __exit__ ( self , *_ ) :
        import sys, copy
        sys.argv = copy.deepcopy ( self._args )
        del self._args 


# =============================================================================
## context manager to keep/preserve sys.argv
#  @code
#  with keepArgs() :
#    ...  
#  @endcode
def keepArgs() :
    """Context manager to keep/preserve sys.argv
    >>> with keepArgs() :
    ...  
    """
    return KeepArgs()

# =============================================================================
## EnableImplicitMT
#  Context manager to enable/disable implicit MT in ROOT 
#  @see ROOT::EnableImplicitMT 
#  @see ROOT::DisableImplicitMT 
#  @see ROOT::IsImplicitMTEnabled
#  @code
#  with ImplicitMT( True ) :
#  ...
#  @endcode
class ImplicitMT(object) :
    """Context manager to enable/disable implicit MT in ROOT 
    >>> with ImplicitMT( True ) :
        ...
    - see ROOT::EnableImplicitMT 
    - see ROOT::DisableImplicitMT 
    - see ROOT::IsImplicitMTEnabled
    """
    def __init__  ( self , enable = True ) :

        if   isinstance ( enable , bool ) : 
            self.__enable   =        enable
            self.__nthreads =        0
        elif isinstance ( enable , int  ) and 0 <= enable : 
            self.__enable   = bool ( enable ) 
            self.__nthreads =        enable 
        else :
            raise  TypeError ( "ImplicitMT: invalid ``enable'' flag :%s/%s" % ( enable , type ( enable ) ) )

    @property
    def enable   ( self ) : return self.__enable
    @property
    def nthreads ( self ) : return self.__nthreads
    
    ## Context manager: ENTER 
    def __enter__ ( self ) :
            
        self.__initial = ROOT.ROOT. IsImplicitMTEnabled ()
        
        if bool ( self.__initial ) == bool ( self.enable ) : pass 
        elif self.enable : ROOT.ROOT.EnableImplicitMT  ( self.__nthreads )
        else             : ROOT.ROOT.DisableImplicitMT ()

        return self
    
    ## Context manager: EXIT
    def __exit__ ( self , *_ ) :

        _current = ROOT.ROOT.IsImplicitMTEnabled()

        if   _current == self.__initial : pass
        elif _current                   : ROOT.ROOT.DisableImplicitMT ()
        else                            : ROOT.ROOT.EnableImplicitMT  ()
            

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
def counted ( f ):
    """create 'counted' function to know number of function calls

    Example
    -------
    
    >>> fun = ...
    >>> func = counted ( fun ) ## use as function

    >>> @counted
    >>> def fun2 ( ...  ) : return ...
    """
    def wrapped ( *fargs, **kwargs ):
        wrapped.calls += 1
        return f( *fargs , **kwargs )
    wrapped.calls = 0
    return wrapped

# =============================================================================
## Context manager to enable/disable implicit MT in ROOT 
#  @see ROOT::EnableImplicitMT 
#  @see ROOT::DisableImplicitMT 
#  @see ROOT::IsImplicitMTEnabled
#  @code
#  with implicitMT( True ) :
#  ...
#  @endcode
def implicitMT ( enable = True ) :
    """Context manager to enable/disable implicit MT in ROOT 
    >>> with implicitMT( True ) :
        ...
    - see ROOT::EnableImplicitMT 
    - see ROOT::DisableImplicitMT 
    - see ROOT::IsImplicitMTEnabled
    """
    return ImplicitMT ( enable ) 

# =============================================================================
## Return the path to an executable which would be run if the given <code>cmd</code> was called.
#  If no <code>cmd</code> would be called, return <code>None</code>.
#  - <code>mode</code> is a permission mask passed to <code>os.access()</code>,
#    by default determining if the file exists and executable.
#  - When no <code>path</code> is specified, the results of <code>os.environ()</code> are used,
#    returning either the <code>“PATH”</code> value or a fallback of <code>os.defpath</code>.
#  - copied from <code>shutil</cdde> module
def local_which ( cmd, mode=os.F_OK | os.X_OK, path=None):
    """Given a command, mode, and a PATH string, return the path which
    conforms to the given mode on the PATH, or None if there is no such
    file.

    `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
    of os.environ.get("PATH"), or can be overridden with a custom search
    path.

    """
    # Check that a given file can be accessed with the correct mode.
    # Additionally check that `file` is not a directory, as on Windows
    # directories pass the os.access check.
    def _access_check(fn, mode):
        return (os.path.exists(fn) and os.access(fn, mode)
                and not os.path.isdir(fn))

    # If we're given a path with a directory part, look it up directly rather
    # than referring to PATH directories. This includes checking relative to the
    # current directory, e.g. ./script
    if os.path.dirname(cmd):
        if _access_check(cmd, mode):
            return cmd
        return None

    if path is None:
        path = os.environ.get("PATH", os.defpath)
    if not path:
        return None
    path = path.split(os.pathsep)

    if sys.platform == "win32":
        # The current directory takes precedence on Windows.
        if not os.curdir in path:
            path.insert(0, os.curdir)

        # PATHEXT is necessary to check on Windows.
        pathext = os.environ.get("PATHEXT", "").split(os.pathsep)
        # See if the given file matches any of the expected path extensions.
        # This will allow us to short circuit when given "python.exe".
        # If it does match, only test that one, otherwise we have to try
        # others.
        if any(cmd.lower().endswith(ext.lower()) for ext in pathext):
            files = [cmd]
        else:
            files = [cmd + ext for ext in pathext]
    else:
        # On other platforms you don't have things like PATHEXT to tell you
        # what file suffixes are executable, so just pass on cmd as-is.
        files = [cmd]

    seen = set()
    for dir in path:
        normdir = os.path.normcase(dir)
        if not normdir in seen:
            seen.add(normdir)
            for thefile in files:
                name = os.path.join(dir, thefile)
                if _access_check(name, mode):
                    return name
    return None

# =============================================================================

try :
    from shutil import which
except ImportError :    
    which  = local_which 

# =============================================================================
## get the command
#  @code
#  >>> if cmd_exists ( 'epstopdf' ) : ... 
#  @endcode 
#  @see https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def cmd_exists ( command ) :
    """Check the existence of certain command/executable
    >>> if cmd_exists ( 'epstopdf' ) : ...
    
    """
    return which ( command ) is not None

# =============================================================================
## @class VRange
#  Helper looper over the values between vmin and vmax :
#  @code
#  for v in VRange ( vmin = 0 , vmax = 5 , n = 100 ) :
#  ... print ( v ) 
#  @endcode 
class VRange(object) :
    """Helper looper over the values between vmin and vmax :
    >>> for v in VRange ( vmin = 0 , vmax = 5 , n = 100 ) :
    >>> ... print ( v ) 
    """
    def __init__ ( self , vmin , vmax , n = 100 , edges = True ) :
        
        assert isinstance ( n , integer_types ) and 0 < n,\
               'VRange: invalid N=%s/%s' % ( n  , type ( n ) ) 
        
        self.__vmin  = min ( vmin , vmax )
        self.__vmax  = max ( vmin , vmax ) 
        self.__n     = n
        self.__edges = True if edges else False  

    @property
    def edges ( self ) :
        """`edges`: include edges?"""
        return self.__edges
    
    @property
    def vmin  ( self ) :
        """`vmin' : minimal value"""
        return self.__vmin
    @property
    def vmax  ( self ) :
        """`vmax' : maximal value"""
        return self.__vmax    
    @property
    def n ( self )  :
        """`n' : number of steps"""
        return self.__n

    def __len__     ( self ) :

        n = self.__n
        e = self.__edges 
        return n + 1 if e else n - 1 

    def __iter__    ( self ) :
        
        n  = self.n 
        fn = 1.0 / float ( n ) 
        e  = self.edges

        vmn = self.vmin
        vmx = self.vmax
        
        if e : yield vmn

        for i in range ( 1 , n ) :
            f2 = i * fn
            f1 = 1 - f2
            yield vmn * f1 + f2 * vmx
            
        if e : yield vmx
                
# =============================================================================
## loop over values between xmin and xmax 
#  @code
#  for x in vrange ( xmin , xmax , 200 ) :
#         print (x) 
#  @endcode
def vrange ( vmin , vmax , n = 100 , edges = True ) :
    """ Loop  over range of values between xmin and xmax 
    >>> for v in vrange ( vmin , vmax , 200 ) :
    ...                print (v) 
    """
    return VRange ( vmin , vmax , n , edges )


# =============================================================================
## @class LRange
#  Helper looper over the values between vmin and vmax using log-steps 
#  @code
#  for v in LRange ( vmin = 1 , vmax = 5 , n = 100 ) :
#  ... print ( v ) 
#  @endcode 
class LRange(VRange) :
    """Helper looper over the values between vmin and vmax using log-steps
    >>> for v in LRange ( vmin = 1 , vmax = 5 , n = 100 ) :
    >>> ... print ( v ) 
    """
    def __init__ ( self , vmin , vmax , n = 100 , edges = True ) :
        
        assert 0 < vmin  and 0 < vmax,\
           'LRange: invalid  non-positive vmin/ymax values: %s/%s' %  ( vmin , vmax )

        super ( LRange , self ).__init__ ( vmin , vmax , n , edges ) 

        self.__lmin = math.log10 ( self.vmin )
        self.__lmax = math.log10 ( self.vmax )
                
    @property
    def lmin  ( self ) :
        """``lmin'' : log10(minimal value)"""
        return self.__lmin
    @property
    def lmax  ( self ) :
        """``lmax'' : log10(maximal value)"""
        return self.__lmax    

    def __iter__    ( self ) :

        n  = self.n 
        fn = 1.0 / float ( n )
        e  = self.edges

        lmn = self.__lmin
        lmx = self.__lmax
        
        if e : yield self.vmin
        
        for i in range ( 1 , n  ) :
            #
            f2 = i * fn
            f1 = 1 - f2
            yield 10.0 ** ( lmn * f1 + f2 * lmx ) 
            
        if e : yield self.vmax

# =============================================================================
## loop over values between xmin and xmax in log-scale 
#  @code
#  for x in log_range ( xmin , xmax , 200 ) :
#         print (x) 
#  @endcode
def log_range ( vmin , vmax , n = 100 , edges = True ) :
    """Loop over values between xmin and xmax in log-scale 
    >>> for x in log_range ( xmin , xmax , 200 ) :
    >>>      print (x) 
    """
    return LRange ( vmin , vmax , n , edges )

# =============================================================================
## loop over values between xmin and xmax in log-scale 
#  @code
#  for v in lrange ( vmin , vmax , 200 ) : ## ditto 
#         print (v) 
#  @endcode
def lrange ( vmin , vmax , n = 100 , edges = True ) :
    """:oop over values between vmin and vmax in log-scale 
    >>> for v in lrange ( vmin , vmax , 200 ) :  ## ditto 
    >>>      print (v) 
    """
    return LRange ( vmin , vmax , n , edges )

# =============================================================================
## @class CRange
#  Generate sequence of numbers between vmin and vmax according to Chebyshev nodes
#  It can be useful for e.g. interpolation nodes
#  @code
#  for c in CRange(-1,1,10) : print ( c ) 
#  @endcode 
class CRange(VRange):
    """ Generate sequence of numbers between vmin and vmax accrording to Chebyshev nodes
    It can be useful for e.g. interpolation nodes
    >>> for c in CRange(-1,1,10) : print ( c ) 
    """
    __nodes = {}
    
    ## number of nodes 
    def __len__     ( self ) : return self.n 
    ## 
    def __iter__    ( self ) :

        n     = self.n
        
        if not n in self.__nodes :
            nodes = array.array ( 'd' , n * [ 0.0 ] )
            n2    = n // 2 
            pn    = math.pi / ( 2 * n )
            for k in range ( n2 ) :
                vv  = math.cos ( ( 2 * k + 1 ) * pn )
                nodes [     k     ] =  vv
                nodes [ n - k - 1 ] = -vv
            self.__nodes [ n ] = nodes 

        nodes = self.__nodes [ n ]        
        mid   = 0.5 * ( self.vmin + self.vmax )
        scale = 0.5 * ( self.vmin - self.vmax )              
        for x in nodes :             
            yield mid + scale * x

# =============================================================================
### Generate sequence of numbers between vmin and vmax accrording to Chebyshev nodes
#  It can be useful for e.g. interpolation nodes
#  @code
#  for c in crange(-1,1,10) : print ( c ) 
#  @endcode 
def crange ( vmin , vmax , n = 10 ) :
    """ Generate sequence of numbers between vmin and vmax accrording to Chebyshev nodes
    It can be usefor for e.g. interpolation nodes
    >>> for c in crange(-1,1,10) : print ( c ) 
    """
    return CRange ( vmin , vmax , n )

# =============================================================================
## @class RRange
#  Helper looper over the random values between vmin and vmax
#  @code
#  for v in RRange ( vmin = 1 , vmax = 5 , n = 100 , edges = True  ) :
#  ... print ( v ) 
#  @endcode 
class RRange(VRange) :
    """Helper looper over the values between vmin and vmax using log-steps
    >>> for v in RRange ( vmin = 1 , vmax = 5 , n = 100 , edges = True ) :
    >>> ... print ( v ) 
    """
    def __iter__    ( self ) :
        
        n   = self.n
        vmn = self.vmin
        vmx = self.vmax
        e   = self.edges

        if e : yield vmn
        
        for i in range ( 1 , n ) :
            yield random.uniform ( vmn, vmx )

        if e : yield vmx
            
# =============================================================================
## Generate sequence of random numbers between vmin and vmax
#  @code
#  for c in rrange(-1,1,10) : print ( c ) 
#  @endcode 
def rrange ( vmin , vmax , n = 10 , edges = True  ) :
    """ Generate random sequence of numbers between vmin and vmax 
    >>> for c in crange(-1,1,10, edges = True ) : print ( c ) 
    """
    return RRange ( vmin , vmax , n , edges )

# =============================================================================
## @class PRange
#  Helper looper over the values between vmin and vmax with non=unifomr power-law) distributed poinnts 
#  @code
#  for v in PRange ( vmin = 1 , vmax = 5 , n = 100 , power = 2 ) , edges = True :
#  ... print ( v ) 
#  @endcode 
class PRange(VRange) :
    """Helper looper over the values between vmin and vmax with non-unifomr  (power-low) points 
    >>> for v in PRange ( vmin = 1 , vmax = 5 , n = 100 , power = 2 , edges = True ) :
    >>> ... print ( v ) 
    """
    def __init__ ( self , vmin , vmax , n = 100 , power = 2 , edges = True ) :
        
        assert isinstance ( power , num_types ) and 0 < power , 'PRange: Invalid powr %s' % power  

        super ( PRange , self ).__init__ ( vmin , vmax , n , edges ) 

        self.__power = power

    @property
    def power ( self ) :
        """`power`: the power-low exponent """
        return self.__power
    
    def __iter__    ( self ) :

        n     = self.n 
        fn    = 1.0 / float ( n ) 
        e     = self.edges
        
        p     = self.power
        vmn   = self.vmin
        vmx   = self.vmax
        delta = vmx - vmn
        
        if e : yield vmn

        for i in range ( 1 , n ) :
            x = i * fn
            yield vmn + delta * ( x ** p ) 

        if e : yield vmx


# =============================================================================
## Loop over sequence of non-uniformly distributed  (power-low) unmpers
#  @code
#  for c in prange(-1,1,10, power = 2 , edges = True ) : print ( c ) 
#  @endcode 
def prange ( vmin , vmax , n = 10 , power = 2 , edges = True  ) :
    """ Loop over the sequence of non-unifrmy distribited (power-low) numbers 
    >>> for c in prange(-1,1,10, power = 2 , edges = True  ) : print ( c ) 
    """
    return PRange ( vmin , vmax , n , power = power  , edges = edges )


# =============================================================================
## split range into smaller chunks:
#  @code
#  for i in SplitRange ( 0 , 10000 , 200 ) :
#     for j in range (*i) :
#          ... 
#  @endcode 
class SplitRange(object) :
    """Split range into smaller chunks:
    >>> for i in SplitRange ( 0 , 10000 , 200 ) :
    >>>     for j in range (*i) :
    >>>         ... 
    """
    def __init__ ( self , low , high , num ) :

        self.__low  = low
        self.__high = high 
        self.__num  = num

        self.__size = 0 
        if low < high and 1 <= num :
            self.__size , r = divmod ( self.__high - self.__low , self.__num )
            if r : self.__size += 1 

    def __iter__ ( self ) :
        
        if 1 <= self.__num : 
            while self.__low < self.__high :
                next  = min ( self.__low + self.__num , self.__high ) 
                yield self.__low , next
                self.__low = next 
                
    def __len__ ( self ) :
        return self.__size 
    

# =============================================================================
## split range into smaller chunks:
#  @code
#  for i in split_range ( 0 , 10000 , 200 ) :
#     for j in range (*i) :
#          ... 
#  @endcode 
def split_range ( low , high , num ) :
    """Split range into smaller chunks:
    >>> for i in split_range ( 0 , 10000 , 200 ) :
    >>>     for j in range (*i) :
    >>>         ... 
    """
    return SplitRange ( low , high , num )

# =============================================================================
## split range into num smaller chunks of approximate size 
#  @code 
#  for i in split_n_range ( 0 , 10000 , 200 ) :
#     for j in range (*i) :
#          ... 
#  @endcode 
def split_n_range ( low , high , num ) :
    """Split range into `num` smaller chunks of approximate size 
    >>> for i in split_n_range ( 0 , 10000 , 200 ) :
    >>>     for j in range (*i) :
    >>>         ... 
    """

    if   high <= low or num < 1            : pass 
    elif 1 == num                          : yield low , high
    elif low < high  and high <= num + low : yield low , high
    else : 
        
        nn   = high - low
        newn = nn // num
        for i in range ( 0 , num - 1 ) :
            nl = i  * newn
            nh = nl + newn
            yield low + nl , low + nh 
        yield low + num * newn - newn , high
        
# =============================================================================
if (3,6) <= sys.version_info :
    
    choices = random.choices
    
else :
    
    def choices ( population , weights = None , cum_weights = None , k = 1 ) :
        """ Simple variant of `random.choice`
        """
        assert weights is None and cum_weights is None,\
               "choices: Neither ``weigths'' nor ``cum_weights'' are supported!"
        
        return [ random.choice ( population ) for i in range ( k ) ] 
    
# ========================================================================================
## Generate some random name of given name
#  @code
#  name = random_name ( 5 ) 
#  @endcode 
def random_name ( size ) :
    """Generate some random name of given name 
    >>> name = random_name ( 5 )
    """
    assert 1 <= size , 'random_name: invalid size!'

    first = random.choice  ( ascii_letters ) 
    if 1 == size : return first
    
    return first  + ''.join ( choices ( sll_symbols , k = size - 1 ) ) 

# ========================================================================================
## generate some pseudo-random 6-symbol name from provided hash sources 
def short_hash_name ( size , name , *names ) :
    """generate some pseudo-random 6-symbol name from provided hash sources
    """

    size = max ( min ( size , 8 ) , 4 ) 

    h = size , hash ( tuple ( ord ( i ) for i in name ) )
    h = hash ( h ) 
                   
    for n in names :

        h = h , hash ( tuple ( ord ( i ) for i in n ) )
        h = hash ( h ) 
        
    h = abs ( h ) % ( 2 ** ( 4 * size ) )
    
    return ( '%%0%dx' % size ) % h 
    
# =============================================================================
## Generate the random string, that can be used as password or secret word
#  @code
#  password = gen_password () 
#  @endcode 
def gen_password ( size = 12 ) :
    """Generate the random string, that can be used as password or secret word
    >>> password = gen_password () 
    """
    import random
    ## save random state 
    state = random.getstate ()
    ## reset the random seed
    random.seed ()
    ## generate the password 
    result = ''.join ( choices ( all_symbols , k = size ) ) 
    ## restore the random state 
    random.setstate ( state )
    ## 
    return result

# =============================================================================

try :
    
    from more_itertools import chunked, divide 
    
except ImportError :
    
    from itertools import islice
    from functools import partial
    
    # =========================================================================
    ## Return first *n* items of the iterable as a list
    #  @code 
    #  take(3, range(10))  ## [0, 1, 2]
    #  take(5, range(3))   ## [0, 1, 2]
    #  @endcode
    #
    #  The function is copied from <code>more_itertools</code> 
    def take(n, iterable):
        """Return first *n* items of the iterable as a list.
        
        >>> take(3, range(10))
        [0, 1, 2]
        >>> take(5, range(3))
        [0, 1, 2]
        
        Effectively a short replacement for ``next`` based iterator consumption
        when you want more than one item, but less than the whole iterator.
        
        - the function is copied from `more_itertools`
        """
        return list(islice(iterable, n))
    
    # =========================================================================
    ## Break *iterable* into lists of length *n*:
    #  @code
    #  list(chunked([1, 2, 3, 4, 5, 6], 3)) ## [[1, 2, 3], [4, 5, 6]]
    #  @endcode
    #  If the length of *iterable* is not evenly divisible by *n*, the last
    #  returned list will be shorter:
    #  @code 
    #  list(chunked([1, 2, 3, 4, 5, 6, 7, 8], 3)) ## [[1, 2, 3], [4, 5, 6], [7, 8]]
    #  @endcode 
    #  <code>chunked</code> is useful for splitting up a computation on a large number
    #  of keys into batches, to be pickled and sent off to worker processes. One
    #  example is operations on rows in MySQL, which does not implement
    #  server-side cursors properly and would otherwise load the entire dataset
    #  into RAM on the client.
    # 
    #  The function is copied from <code>more_itertools</code>
    def chunked  ( iterable , n ):
        """Break *iterable* into lists of length *n*:
        
        >>> list(chunked([1, 2, 3, 4, 5, 6], 3))
        [[1, 2, 3], [4, 5, 6]]
        
        If the length of *iterable* is not evenly divisible by *n*, the last
        returned list will be shorter:
        
        >>> list(chunked([1, 2, 3, 4, 5, 6, 7, 8], 3))
        [[1, 2, 3], [4, 5, 6], [7, 8]]
        
        To use a fill-in value instead, see the :func:`grouper` recipe.
        
        :func:`chunked` is useful for splitting up a computation on a large number
        of keys into batches, to be pickled and sent off to worker processes. One
        example is operations on rows in MySQL, which does not implement
        server-side cursors properly and would otherwise load the entire dataset
        into RAM on the client.
        
        - the function is copied from `more_itertools`
        """
        return iter ( partial ( take , n , iter ( iterable ) ) , [] )

    # =========================================================================
    ## Divide the elements from *iterable* into *n* parts, maintaining order.
    #  @code 
    #  >>> group_1, group_2 = divide(2, [1, 2, 3, 4, 5, 6])
    #  >>> list(group_1)
    #  ...    [1, 2, 3]
    #  >>> list(group_2)
    #  ... [4, 5, 6]
    #  @endcode
    #  If the length of *iterable* is not evenly divisible by *n*, then the
    #  length of the returned iterables will not be identical:
    #  @code 
    #  >>> children = divide(3, [1, 2, 3, 4, 5, 6, 7])
    #  >>> [list(c) for c in children]
    #  ... [[1, 2, 3], [4, 5], [6, 7]]
    #  @endcode
    # 
    # If the length of the iterable is smaller than n, then the last returned
    # iterables will be empty:
    # @code
    # >>> children = divide(5, [1, 2, 3])
    # >>> [list(c) for c in children]
    # ... [[1], [2], [3], [], []]
    # @endcode
    # 
    # This function will exhaust the iterable before returning and may require
    # significant storage. If order is not important, see :func:`distribute`,
    # which does not first pull the iterable into memory.
    #
    # The function is copied from <code>more_itertools</code>
    def divide ( n , iterable):
        """Divide the elements from *iterable* into *n* parts, maintaining
        order.
        
        >>> group_1, group_2 = divide(2, [1, 2, 3, 4, 5, 6])
        >>> list(group_1)
        [1, 2, 3]
        >>> list(group_2)
        [4, 5, 6]
        
        If the length of *iterable* is not evenly divisible by *n*, then the
        length of the returned iterables will not be identical:
        
        >>> children = divide(3, [1, 2, 3, 4, 5, 6, 7])
        >>> [list(c) for c in children]
        [[1, 2, 3], [4, 5], [6, 7]]
        
        If the length of the iterable is smaller than n, then the last returned
        iterables will be empty:
        
        >>> children = divide(5, [1, 2, 3])
        >>> [list(c) for c in children]
        [[1], [2], [3], [], []]
        
        This function will exhaust the iterable before returning and may require
        significant storage. If order is not important, see :func:`distribute`,
        which does not first pull the iterable into memory.
        
        - the function is copied from `more_itertools`
        """
        if n < 1:
            raise ValueError('n must be at least 1')

        seq = tuple(iterable)
        q, r = divmod(len(seq), n)
        
        ret = []
        for i in range(n):
            start = (i * q) + (i if i < r else r)
            stop = ((i + 1) * q) + (i + 1 if i + 1 < r else r)
            ret.append(iter(seq[start:stop]))
            
        return ret

# =============================================================================
if ( 3 , 0 ) <= python_version : from itertools import zip_longest
else                           : from itertools import izip_longest as zip_longest

# =============================================================================
## Collect data into fixed-length chunks or blocks"
def grouper ( iterable , n , fillvalue = None ):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    nargs =  [ iter ( iterable ) ]  * n
    return zip_longest ( *nargs , fillvalue = fillvalue )

# =============================================================================
## Create (infinite or finite) iterable from other iterable or non-iterable
#  @code
#  for a,b in zip ( 'abcde' , make_iterable ( 1 ) ) : ...
#  for a,b in zip ( 'abcde' , make_iterable ( [1,2,3] ) ) : ...
#  for a,b in zip ( 'abcde' , make_iterable ( [1,2,3] , 0 , 2 ) ) : ...
#  @endcode 
def make_iterable ( what , default = None , size = -1 ) :
    """Create infinite or finite iterable from other iterable or no-iterable
    >>> for a,b in zip ( 'abcde' , make_iterable ( 1 ) ) : ...
    >>> for a,b in zip ( 'abcde' , make_iterable ( [1,2,3] ) ) : ...
    >>> for a,b in zip ( 'abcde' , make_iterable ( [1,2,3] , 0 , 2 ) ) : ...
    """
    from   ostap.core.ostap_types import iterable_types

    if not isinstance ( what , iterable_types ) : what = what, 

    ## make infinite iterable 
    result = chain ( what , repeat ( default ) )

    ## cut it, if needed 
    return result if size < 0 else islice ( result , size )

# =============================================================================
## calculate SHA512-checksum for the files
#  @see hashlib
#  @see hashlib.sha512
#  @code
#  s =  checksum_files ( 'a.txt', 'b.bin' ) 
#  @endcode
#  Non-existing files are ignored 
#  @param files  list of filenames
#  @return checksum for these files 
def checksum_files ( *files ) :
    """Calculate SHA512-checksum for the files
    >>> s =  checksum_files ( 'a.txt', 'b.bin' ) 
    Non-existing files are ignored 
    - see `hashlib`
    - see `hashlib.sha512`
    """
    import hashlib
    hash_obj = hashlib.sha512 ()
    for fname in files :
        if os.path.exists ( fname ) and os.path.isfile ( fname ) : 
            with open ( fname , "rb" ) as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_obj.update(chunk)
                    
    return hash_obj.hexdigest()

# =============================================================================
## Simple utility to check balanced parenthesis/brackets, etc...
#  @code
#  expression = ' .... '
#  ok = balanced ( expression ) 
#  @endcode 
def  balanced ( expression , left = '([' , right = ')]' ) :
    """Simple utility to check balanced parenthesis/brackets, etc...
    >>> expression = ' .... '
    >>> ok = balanced ( expression ) 
    """
    
    assert left and len(left) == len ( right ) ,\
           'balanced: invalid left/right arguments!'
    
    stack = []
    for i in expression :
        if   i in left  : stack.append ( i )
        elif i in right :
            pos = right.index ( i )
            if stack  and  left[ pos ] == stack [ -1 ] :
                stack.pop()
            else :
                return False

    return True if not stack else False 

# =============================================================================


# ============================================================================
if   ( 3 , 9) <= sys.version_info : memoize = functools.cache 
elif ( 3 , 2) <= sys.version_info : 
    
    # =========================================================================
    ## Simple lightweight unbounded cache
    def memoize ( user_function ):
        """Simple lightweight unbounded cache"""
        return functools.lru_cache(maxsize=None)(user_function)

else :

    from ostap.utils.basic import loop_items
    
    # =========================================================================
    ## Simple lightweight unbounded cache
    class memoize(object):
        """Simple lightweight unbounded cache
        """
        def __init__ ( self, func ) :
            
            self.func  = func
            self.cache = {}
            functools.update_wrapper ( self , func ) 
            
        def __call__ ( self, *args, **kwargs ):

            if kwargs : all_args = args , tuple ( loop_items ( kwargs ) ) 
            else      : all_args = args 
                    
            if all_args in self.cache:    
                return self.cache [ all_args ]
            
            value = self.func ( *args , **kwargs )
            self.cache [ all_args ] = value
            
            return value
        
        def __repr__(self):
            return self.func.__doc__
        def __get__(self, obj, objtype):
          return functools.partial(self.__call__, obj)
      
    
# ============================================================================
## abstract prpoperty
#  @code
#  @absproperty
#  def A ( self ) : ...  
#  @endcode
if ( 3 , 3 ) <= sys.version_info  :
    
    # =========================================================================
    ## absract property decorator
    #  @code
    #  @absproperty
    #  def A ( self ) : ...  
    #  @endcode
    def absproperty ( func ) :
        """Abstract property
        @absproperty
        def A ( self ) : ...
        """
        return property ( abc.abstractmethod ( func ) ) 
        
else :

    import abc 
    # =========================================================================
    ## abstract prpoperty
    #  @code
    #  @absproperty
    #  def A ( self ) : ...  
    #  @endcode
    absproperty = abc.abstractproperty

    
# =============================================================================
if  ( 3 , 9 ) <= sys.version_info :
    
    # =========================================================================
    ## class property decorator
    #  @code
    #  @classprop
    #  def A ( cls ) : ...  
    #  @endcode
    def classprop ( func ) :
        """Class property
        @classprop
        def A ( cls ) : ...
        """
        return classmethod ( property ( func ) ) 
        
elif ( 3 , 0 ) <= sys.version_info : 

    # =========================================================================
    ## class @classproperty
    #  class property decorator  (copied and simplified from astropy)
    #  @code
    #  @classprop
    #  def A ( cls ) : ...  
    #  @endcode
    class classprop(property):
        """Class property
        @classprop
        def A ( cls ) : ...
        """        
        def __new__(cls, fget=None, doc=None):
            if fget is None:
                # Being used as a decorator--return a wrapper that implements
                # decorator syntax
                def wrapper(func):
                    return cls(func)            
                return wrapper
            return super(classprop,cls).__new__(cls)
        
        def __init__(self, fget, doc=None, ):
            fget = self._wrap_fget(fget)
            super(classprop,self).__init__(fget=fget, doc=doc)
            
            # There is a buglet in Python where self.__doc__ doesn't
            # get set properly on instances of property subclasses if
            # the doc argument was used rather than taking the docstring
            # from fget
            # Related Python issue: https://bugs.python.org/issue24766
            if doc is not None:
                self.__doc__ = doc
                
        def __get__(self, obj, objtype):
            # The base property.__get__ will just return self here;
            # instead we pass objtype through to the original wrapped
            # function (which takes the class as its sole argument)
            val = self.fget.__wrapped__(objtype)
            return val

        def getter(self, fget):
            return super(classprop,self).getter(self._wrap_fget(fget))
        
        def setter(self, fset):
            raise NotImplementedError(
                "classproperty can only be read-only; use a metaclass to "
                "implement modifiable class-level properties")
        
        def deleter(self, fdel):
            raise NotImplementedError(
                "classproperty can only be read-only; use a metaclass to "
                "implement modifiable class-level properties")
        
        @staticmethod
        def _wrap_fget(orig_fget):
            if isinstance(orig_fget, classmethod):
                orig_fget = orig_fget.__func__
                
            # Using stock functools.wraps instead of the fancier version
            # found later in this module, which is overkill for this purpose
            
            @functools.wraps(orig_fget)
            def fget(obj):
                return orig_fget(obj.__class__)
            
            return fget
else :
    
    class classprop(object):
        def __init__(self, fget):
            self.fget = fget
        def __get__(self, inst, cls):
            return self.fget(cls)

# =============================================================================
## @class CallThem
#  Use sequence of callables as single callable
#  @code
#  fun1, fun2 , func3 = ....
#  new_fun1 = CallThem ( fun1, fun2, fun3 )
#  print ( new_fun ( x ) ) 
#  new_fun2 = CallThem ( [ fun1, fun2, fun3 ]  )
#  print ( new_fun2 ( x ) ) 
#  @endcode
class CallThem(object) :
    """ Use sequence of callables as a single callable
    >>> fun1, fun2 , func3 = ....
    >>> new_fun1 = CallThem ( fun1, fun2, fun3 )
    >>> print ( new_fun1 ( x ) )
    >>> new_fun2 = CallThem ( [ fun1, fun2, fun3 ] )
    >>> print ( new_fun2 ( x ) )
    """
    # ========================================================================
    def __init__ ( self , *callables ) :
        
        from ostap.core.ostap_types import sequence_types
        if 1 == len ( callables ) and \
           ( not callable ( callables [ 0 ] ) ) and \
           isinstance ( callables [0] , sequence_types ) : 
            
            callables = tuple ( c for c in callables [ 0 ] )

        ## here all arguments must be callables! 
        assert all ( callable ( c ) for c in callables ) , \
               'All parameters must be callables!'
        
        self.__callables = callables
        
    ## call all callables! 
    def __call__ ( self , *args , **kwargs ) :
        """call all callables!"""
        return tuple ( c ( *args , **kwargs ) for c in self.__callables ) 

    @property
    def callables ( self ) :
        """'callables' : get all callables"""
        return self.__callables 

    def __str__ ( self ) :
        return ','.join ( str ( c ) for c in self.__callables )
    __repr__ = __str__

# =============================================================================
## @class NumCalls
#  Count a number of  times a callable object is invoked
class NumCalls (object):
    """Count a number of  times a callable object is invoked"""
    def __init__ ( self , func ) :
        self.__func  = func
        self.__count = 0
        functools.update_wrapper ( self, func ) 
    def __call__ ( self, *cargs , **kwargs ) :
        self.__count +=1
        return self.__func ( *cargs , **kwargs )
    @property
    def count ( self ) :
        """'count': number of times the function was invoked"""
        return self.__count
    
# ==============================================================================        
#  Count a number of  times a callable object is invoked
numcalls = NumCalls


# ==============================================================================
## @class RandomSeed 
#  Context manager to set/keep seed/state of PYTHON random generator
#  @code
#  import random 
#  with RandomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, random.uniform (0,1) )
#  @endcode
#  Seed can be None, int,float,string, bytes, .
class RandomSeed(object) :
    """Context manager to set/keep seed/state of PYTHON random generator
    >>> import random 
    >>> with RandomSeed ( 'this is seed' ) :
    >>> ... for i in range ( 10 ) : print ( i, random.uniform (0,1) )
    - Seed can be None, int,float,string, bytes, .
    """
    def __init__ ( self , seed = None )  :
        
        self.__seed  = seed
        self.__state = None  

    def __enter__  ( self ) :
        
        self.__state = random.getstate ()
        random.seed ( self.__seed )
        return self
    
    def __exit__   ( self , *_ ) :
        random.setstate ( self.__state )
        
    @property
    def seed ( self ) :
        """'seed' : the seed for random
        """
        return self.__seed 
    
# ==============================================================================
## Context manager to set/keep seed/state of random generator
#  @code
#  import random 
#  with randomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, random.uniform (0,1) )
#  @endcode
#  Seed can be None, int,float,string, bytes, ... 
def randomSeed ( seed = None ) :
    """Context manager to set/keep seed/state of random generator
    >>> import random 
    >>> with randomSeed ( 'this is seed' ) :
    >>> ... for i in rnage ( 10 ) : print ( i, random.uniform (0,1) )
    - Seed can be None, int,float,string, bytes, ...
    """
    return RandomSeed ( seed )

# ==============================================================================
## Context manager to set/keep seed/state of random generator
#  @code
#  import random 
#  with random_seed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, random.uniform (0,1) )
#  @endcode
#  Seed can be None, int,float,string, bytes, ...
random_seed = randomSeed     

# ==============================================================================
## @class RootRandomSeed 
#  Context manager to set/keep seed/state of ROOT random generator
#  @code
#  import random 
#  with RootRandomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, ROOT.gRandom.Rndm () )
#  @endcode
#  Seed can be None, int,float,string, bytes, .
class RootRandomSeed(object) :
    """Context manager to set/keep seed/state of ROOT random generator
    >>> import random 
    >>> with RootRandomSeed ( 'this is seed' ) :
    >>> ... for i in range ( 10 ) : print ( i, ROOT.gRandom.Rndm() )
    - Seed can be None, int,float,string, bytes or any hashable object 
    """
    def __init__ ( self , seed = None )  :
        
        if   seed is None              : self.__seed = 0
        elif isinstance ( seed , int ) : self.__seed = seed
        else                           : self.__seed = hash ( seed ) 
        
        self.__old  = 0   
        
    def __enter__  ( self ) :

        if ROOT.gRandom :
            self.__old = ROOT.gRandom.GetSeed()
            ROOT.gRandom.SetSeed ( self.__seed ) 
            
        return self
    
    def __exit__   ( self , *_ ) :
        if ROOT.gRandom :
            ROOT.gRandom.SetSeed ( self.__old ) 
        
    @property
    def seed ( self ) :
        """'seed' : the actual seed for random
        """
        return self.__seed 


# ==============================================================================
## Context manager to set/keep seed/state of ROOT random generator
#  @code
#  import random 
#  with rootRandomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, ROOT.gRandom.Rndm () )
#  @endcode
#  Seed can be None, int,float,string, bytes, .
def rootRandomSeed ( seed = None ) :
    """Context manager to set/keep seed/state of ROOT random generator
    >>> import random 
    >>> with rootRandomSeed ( 'this is seed' ) :
    >>> ... for i in range ( 10 ) : print ( i, ROOT.gRandom.Rndm() )
    - Seed can be None, int,float,string, bytes or any hashable object 
    """
    return RootRandomSeed ( seed ) 

# ==============================================================================
## Context manager to set/keep seed/state of ROOT random generator
#  @code
#  import random 
#  with rootRandomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, ROOT.gRandom.Rndm () )
#  @endcode
#  Seed can be None, int,float,string, bytes or any hashable object 
root_random_seed = rootRandomSeed


# ==============================================================================
## @class RooRandomSeed 
#  Context manager to set/keep seed/state of RooFit random generator
#  @code
#  import random 
#  with RooRandomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, ROOT.RooRandom.uniform () )
#  @endcode
#  Seed can be None, int,float,string, bytes, .
class RooRandomSeed(object) :
    """Context manager to set/keep seed/state of ROOT random generator
    >>> import random 
    >>> with RooRandomSeed ( 'this is seed' ) :
    >>> ... for i in range ( 10 ) : print ( i, ROOT.RooRandom.uniform() )
    - Seed can be None, int,float,string, bytes or any hashable object
    """

    def __init__ ( self , seed = None )  :
        
        if   seed is None              : self.__seed = 0
        elif isinstance ( seed , int ) : self.__seed = seed
        else                           : self.__seed = hash ( seed ) 
        
        self.__old  = 0   
        

    def __enter__  ( self ) :
        
        rg = ROOT.RooRandom.randomGenerator() 
        if rg :
            self.__old = rg.GetSeed()
            rg.SetSeed ( self.__seed ) 
            
        return self
    
    def __exit__   ( self , *_ ) :
        
        rg = ROOT.RooRandom.randomGenerator() 
        if rg : rg.SetSeed ( self.__old ) 
        
    @property
    def seed ( self ) :
        """'seed' : the actual seed for random
        """
        return self.__seed 

# ==============================================================================
## Context manager to set/keep seed/state of RooFit random generator
#  @code
#  import random 
#  with rooRandomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, ROOT.RooRandom.uniform () )
#  @endcode
#  Seed can be None, int,float,string, bytes, .
def rooRandomSeed ( seed = None ) :
    """Context manager to set/keep seed/state of RooFit random generator
    >>> import random 
    >>> with rooRandomSeed ( 'this is seed' ) :
    >>> ... for i in range ( 10 ) : print ( i, ROOT.RooRandom.uniform() )
    - Seed can be None, int,float,string, bytes or any hashable object 
    """
    return RooRandomSeed ( seed ) 

# ==============================================================================
## Context manager to set/keep seed/state of RooFit random generator
#  @code
#  import random 
#  with roo_random_seed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, ROOT.RooRandom.uniform () )
#  @endcode
#  Seed can be None, int,float,string, bytes or any hashable object 
roo_random_seed = rooRandomSeed


# ==============================================================================
## Copy file with the progress
#  @code
#  copy_with_progress ( 'inputfilename.ext' , 'outputfilename.ext' ) 
#  @endcode
def copy_with_progress ( source  , destination ) :
    """Copy file with progress
    >>> copy_with_progress ( 'inputfilename.ext' , 'outputfilename.ext' ) 
    """
    assert os.path.exists ( source ) and os.path.isfile ( source ), \
           "copy_with_progress: ``source'' %s does nto exist!" % source
    
    total = os.stat ( source ) . st_size
    BLOCK = 512 * 1024

    destination = os.path.abspath  ( destination )    
    destination = os.path.normpath ( destination )
    destination = os.path.realpath ( destination )
    if os.path.exists ( destination ) and os.path.isdir ( destination ) :
        destination = os.path.join ( destination , os.path.basename ( source ) )
        
    from ostap.utils.progress_bar import ProgressBar 
    read = 0
    
    with ProgressBar ( total , silent = total < 3 * BLOCK )  as pbar : 
        with open ( source , 'rb' ) as fin :
            with open ( destination , 'wb' ) as fout :
                while True :
                    
                    block = fin.read ( BLOCK )
                    fout.write ( block )
                    
                    read  += len ( block )
                    pbar.update_amount ( read )
                    if not block : break             ## BREAK
                
    assert os.path.exists ( destination ) and \
           os.path.isfile ( destination ) and \
           os.stat ( destination ).st_size == total, \
           "Invalid ``destination'' %s " % destination
    
    return os.path.realpath ( destination )


# =========================================================================
## Has any of the symbols?
def has_symbol ( expression , symbols ) :
    """Has any of the symbols?"""
    return any ( s in expression for s in symbols ) 
# =========================================================================
## markers of the math formulas 
math_symbols = ' +-*/=><()[]^%&|'
## Is this string expression represend math formula?
def is_formula ( expr , symbols = math_symbols ) :
    """Is this string expression represend math formula?"""
    return has_symbol ( expr.strip() , symbols ) 

# =========================================================================
## merge all files using <code>hadd</code> script from ROOT
#  @param output  name of the output merged file, if None,
#                 the temporary name will be generated,
#                 that will be deleted at the end of the session
#  @param opts   options for command <code>hadd</code>
#  @return the name of the merged file
# OPTIONS:
# -a                                   Append to the output
# -k                                   Skip corrupt or non-existent files, do not exit
# -T                                   Do not merge Trees
# -O                                   Re-optimize basket size when merging TTree
# -v                                   Explicitly set the verbosity level: 0 request no output, 99 is the default
# -j                                   Parallelize the execution in multiple processes
# -dbg                                 Parallelize the execution in multiple processes in debug mode (Does not delete partial files stored inside working directory)
# -d                                   Carry out the partial multiprocess execution in the specified directory
# -n                                   Open at most 'maxopenedfiles' at once (use 0 to request to use the system maximum)
# -cachesize                           Resize the prefetching cache use to speed up I/O operations(use 0 to disable)
# -experimental-io-features            Used with an argument provided, enables the corresponding experimental feature for output trees
# -f                                   Gives the ability to specify the compression level of the target file(by default 4) 
# -fk                                  Sets the target file to contain the baskets with the same compression
#                                      as the input files (unless -O is specified). Compresses the meta data
#                                      using the compression level specified in the first input or the
#                                      compression setting after fk (for example 206 when using -fk206)
# -ff                                  The compression level use is the one specified in the first input
# -f0                                  Do not compress the target file
# -f6                                  Use compression level 6. (See TFile::SetCompressionSettings for the support range of value.)  
def hadd ( files , output = None , dir = None , opts = "-ff -O" ) :
    """Merge all files using <code>hadd</code> script from ROOT
    - `output`  name of the output merged file
    - `opts`   options for command <code>hadd</code>
    It returns the name of the merged file
    
    If no output file name is specified, the temporary name
    will be generate and the temporary file will be deleted
    at the end of the session
    
    OPTIONS:
    # -a                                   Append to the output
    # -k                                   Skip corrupt or non-existent files, do not exit
    # -T                                   Do not merge Trees
    # -O                                   Re-optimize basket size when merging TTree
    # -v                                   Explicitly set the verbosity level: 0 request no output, 99 is the default
    # -j                                   Parallelize the execution in multiple processes
    # -dbg                                 Parallelize the execution in multiple processes in debug mode (Does not delete partial files stored inside working directory)
    # -d                                   Carry out the partial multiprocess execution in the specified directory
    # -n                                   Open at most 'maxopenedfiles' at once (use 0 to request to use the system maximum)
    # -cachesize                           Resize the prefetching cache use to speed up I/O operations(use 0 to disable)
    # -experimental-io-features            Used with an argument provided, enables the corresponding experimental feature for output trees
    # -f                                   Gives the ability to specify the compression level of the target file(by default 4) 
    # -fk                                  Sets the target file to contain the baskets with the same compression
    #                                      as the input files (unless -O is specified). Compresses the meta data
    #                                      using the compression level specified in the first input or the
    #                                      compression setting after fk (for example 206 when using -fk206)
    # -ff                                  The compression level use is the one specified in the first input
    # -f0                                  Do not compress the target file
    # -f6                                  Use compression level 6. (See TFile::SetCompressionSettings for the support range of value.)                            
    """

    if isinstance ( files , string_types ) : files = [ files ]

    import glob
    all_files = []
    for p in files :
        all_files += [ f for f in glob.iglob ( p ) ]

    all_files = list  ( set ( all_files ) )
    all_files.sort()
    all_files = tuple ( all_files )
    
    if not output :
        import ostap.utils.cleanup as CU
        suffix = '.root'
        if files :
            base , suffix = os.path.splitext  ( all_files [0] )
            if base : base    = os.path.basename ( base   )
            if base : suffix  = '-%s%s' % ( base , suffix )            
        output = CU.CleanUp.tempfile ( prefix = 'ostap-hadd-merged-' , suffix = suffix , dir = dir )
            
                            
    cargs    = [ 'hadd' ] + opts.split() + [ output ] + [ f for f in all_files ]

    import subprocess
    subprocess.check_call ( cargs )
    
    if os.path.exists ( output ) and os.path.isfile ( output ) :
        return output 
    
    raise IOError ( "The output file %s does not exist!" % output )


# =============================================================================
def hadd2 ( args ) :

    if   isinstance ( args , string_types     ) : return hadd (  args )
    elif isinstance ( args , listlike_types   ) \
         and all ( ( isinstance ( i , string_types ) and i ) for i in args ) :
        return hadd ( args )    
    elif isinstance ( args , dictlike_types   ) :
        return hadd ( **args )

    return hadd  ( *args ) 


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    logger.info ( 80*'*' ) 
            
# =============================================================================
##                                                                      The END 
# =============================================================================
