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
    'rootException'      , ## context manager to perform ROOT Error -> C++/Python exception
    #
    'Profiler'           , ## context manager to perform profiling 
    'RootError2Exception', ## context manager to perform ROOT Error -> C++/Python exception
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
    ##
    'chunked'            , ## break *iterable* into chunks of length *n*:
    'divide'             , ## divide the elements from *iterable* into *n* parts
    'grouper'            , ## collect data into fixed-length chunks or blocks"
    )

# =============================================================================
import ROOT, time, os , sys, math ## attention here!!
from   builtins            import range
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.utils' )
else                       : logger = getLogger( __name__            )
del getLogger
# =============================================================================
from sys                    import version_info  as python_version 
## timing stuff
from ostap.utils.timing     import timing, timer
## other useful stuff 
from ostap.utils.basic      import isatty, with_ipython
from ostap.core.ostap_types import integer_types 
## ... and more useful stuff 
from ostap.utils.memory     import memory, virtualMemory, Memory 
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
    
    def __exit__  ( self , *args ) :

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
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    for fd in range(0, soft):
        try:
            flags = fcntl.fcntl(fd, fcntl.F_GETFD)
        except IOError:
            continue
        fds.append(fd)
    return fds

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
## helper context manager to activate ROOT Error -> Python exception converter 
#  @see Ostap::Utils::useErrorHandler
#  @see Ostap::Utils::ErrorSentry
#  @code
#  with RootError2Exception() :
#  .... do something here 
#  @endcode 
class RootError2Exception (object) :
    """Helper context manager to activate ROOT Error -> Python exception converter
    #
    with RootError2Exception() :
    ... do something here 
    """
    def __init__ ( self ) :
        import ROOT,cppyy 
        Ostap = cppyy.gbl.Ostap
        self.e_handler  = Ostap.Utils.useErrorHandler 
        self.m_previous = False 

    ## context manager entry point  
    def __enter__ ( self ) :    
        self.m_previous = self.e_handler ( True ) 
        return self
    
    ## context manager exit point
    def __exit__ ( self , *_ ) :    
        if self.m_previous : self.e_handler ( False ) 
        self.m_previous = False 

    def __del__ ( self ) :
        if self.m_previous : self.e_handler ( False ) 
        

# =============================================================================
## helper context manager to activate ROOT Error -> Python exception converter 
#  @see Ostap::Utils::useErrorHandler
#  @see Ostap::Utils::ErrorSentry
#  @code
#  with rootException () :
#  .... do something here 
#  @endcode
def rootException () :
    """Helper context manager to activate ROOT Error -> Python exception converter
    #
    with rootException() :
    ... do something here 
    """
    return RootError2Exception()


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
class KeepCWD(object) :
    """context manager to keep the current working directory
    >>> with KeepCWD( new_dir ) :
    ...
    """
    def __init__ ( self , new_dir = '' ) :
        self.cwd     = os.getcwd() 
        self.new_dir = new_dir
        
    def __enter__ (  self ) :
        self.cwd = os.getcwd() 
        if self.new_cdir and os.path.exists ( self.new_dir ) and os.path.isdir ( self.new_dir ) :
            os.chdir ( self.new_dir )
            return self
        
    def __exit__ ( self , *_ ) :
        if os.path.exists ( self.cwd ) and os.path.isdir ( self.cwd ) :
            os.chdir ( os.cwd ) 
            
# =============================================================================
## context manager to keep the current working directory
#  @code
#  with keepCWD ( new_dir ) :
#    ....
#  @endcode 
def keepCWD ( new_dir = '' ) :
    """context manager to keep the current working directory
    >>> with keepCWD( new_dir ) :
    ...
    """
    return KeepCWD (  new_dir ) 
        
# =============================================================================
## @class KeepCanvas
#  helper class to keep the current canvas
#  @code
#  with KeepCanvas() :
#  ... do something here 
#  @endcode 
class KeepCanvas(object) :
    """Helper class to keep the current canvas
    >>> with KeepCanvas() :
    ... do something here 
    """
    def __enter__ ( self ) :
        import ROOT 
        self.canvas = ROOT.gPad
    def __exit__  ( self , *_ ) :
        if self.canvas:
            self.canvas.cd()

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
    def wrapped ( *args, **kwargs ):
        wrapped.calls += 1
        return f( *args , **kwargs )
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

## return any( os.access ( os.path.join ( path , command  ) , os.X_OK ) for path in os.environ["PATH"].split(os.pathsep) )


# =============================================================================
## loop over values between xmin and xmax 
#  @code
#  for x in vrange ( xmin , xmax , 200 ) :
#         print (x) 
#  @endcode
def vrange ( xmin , xmax , n = 100 ) :
    """ Loop  over range of values between xmin and xmax 
    >>> for x in vrange ( xmin , xmax , 200 ) :
    ...                print (x) 
    """
    assert isinstance ( n , integer_types ) and 0 < n,\
           'vrange: invalid N=%s/%s' % ( n  , type ( n ) ) 

    fn = 1.0 / float ( n ) 
    for i in range ( n + 1 ) :
        #
        if   0 == i : yield xmin
        elif n == i : yield xmax
        else        :
            f2 = i * fn
            f1 = 1 - f2
            yield xmin * f1 + f2 * xmax 

# =============================================================================
## loop over values between xmin and xmax in log-scale 
#  @code
#  for x in log_range ( xmin , xmax , 200 ) :
#         print (x) 
#  @endcode
def log_range ( xmin , xmax , n = 100 ) :
    """:oop over values between xmin and xmax in log-scale 
    >>> for x in log_range ( xmin , xmax , 200 ) :
    >>>      print (x) 
    """
    assert 0 < xmin  and 0 < xmax,\
           'log_range: invalid  xmin/xmax values: %s/%s' %  ( xmin , xmax )

    ## loop
    for x in vrange ( math.log10 ( xmin ) , math.log10 ( xmax ) , n ) :
        yield 10.0**x 

# =============================================================================
## loop over values between xmin and xmax in log-scale 
#  @code
#  for x in lrange ( xmin , xmax , 200 ) : ## ditto 
#         print (x) 
#  @endcode
def lrange ( xmin , xmax , n = 100 ) :
    """:oop over values between xmin and xmax in log-scale 
    >>> for x in lrange ( xmin , xmax , 200 ) :  ## ditto 
    >>>      print (x) 
    """
    for x in log_range ( xmin , xmax , n ) : yield x 
        


# =============================================================================
## Generate the random string, that can be used as password or secret word
#  @code
#  password = gen_password () 
#  @endcode 
def gen_password ( len = 12 ) :
    """Generate the random string, that can be used as password or secret word
    >>> password = gen_password () 
    """
    import random , string
    symbols = string.ascii_letters + string.digits
    ## save random state 
    state = random.getstate ()
    ## reset the random seed
    random.seed ()
    ## generate the password 
    result = ''.join ( random.choice ( symbols ) for i in range ( len ) )
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
    def chunked(iterable, n):
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
        return iter(partial(take, n, iter(iterable)), [])

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
        
if ( 3 , 0 ) <= python_version :
    from itertools import zip_longest
else :
    from itertools import izip_longest as zip_longest

# =============================================================================
## Collect data into fixed-length chunks or blocks"
def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    logger.info ( 80*'*' ) 
            
# =============================================================================
##                                                                      The END 
# =============================================================================
