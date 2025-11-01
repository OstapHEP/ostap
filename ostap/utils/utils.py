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
""" Module with some simple but useful utilities for
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
    'isatty'             , ## is the stream ``isatty'' ?
    'with_ipython'       , ## do we run IPython?
    ##
    'keepArgs'           , ## context manager to keep sys.argv
    ##
    'keepCWD'            , ## context manager to keep current working directory 
    ##
    'KeepArgs'           , ## context manager to keep sys.argv
    ##
    'counted'            , ## decorator to create 'counted'-functioniterable'      , ## create infinite or finite iterable 
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
    'numcalls'           , ## decorator for #ncalls
    ##
    'slow'               , ## "slow" looping with delays at each step
    ##
    'CallThem'           , ## convert sequence of callables into single callable
    'AttrGetter'         , ## helper class to have pickeable `operator.attrgetter`
    ##
    'Singleton'          , ## Metaclass for the singleton 
    ##
    )
# =============================================================================
from   itertools              import repeat, chain, islice
from   ostap.core.meta_info   import python_info 
from   ostap.utils.timing     import timing, timer
from   ostap.utils.basic      import ( isatty   , with_ipython , 
                                      NoContext , zip_longest  , 
                                      counted   , memoize      )  
from   ostap.core.ostap_types import ( integer_types  , num_types ,
                                       string_types   ,
                                       dictlike_types , listlike_types )
from   ostap.utils.memory     import memory, virtualMemory, Memory
import ROOT, time, os , sys, math, time, functools, abc, array, random, datetime, operator  ## attention here!!
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.utils' )
else                       : logger = getLogger( __name__            )
del getLogger
# =============================================================================
try : # =======================================================================
    # =========================================================================
    from string import ascii_letters, digits
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    from string import letters as ascii_letters
    from string import digits
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
    """ Very simple profiler, based on cProfile module
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
## context manager to keep the current working directory
#  @code
#  with KeepCWD ( new_dir ) :
#    ....
#  @endcode
#  - No action if no directory is specified 
class KeepCWD(object) :
    """ Context manager to keep the current working directory
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
    """ Context manager to keep the current working directory
    >>> with keepCWD( new_dir ) :
    ...
    - No action if no directory is specified 
    """
    return KeepCWD ( new_dir ) 

# =============================================================================
## "slow" looping over some iterable with delay for each step
# @code
# for i in slow ( range ( 5 ) , wait = 0.1 ) :
# ... 
# @endcode 
def slow ( iterable , wait = 0 ) :
    """ `slow' looping over some iterable with delay for each step 
    >>> for i in slow ( range ( 5 ) , wait = 0.1 ) :
    >>> ... 
    """
    if isinstance ( iterable , int ) and 0 <= iterable : iterable = range ( iterable )
    ## 
    for r in iterable :
        yield r
        if 0 < wait : time.sleep ( wait )
        
# =============================================================================
## @class KeepArgs
#  context manager to keep/preserve sys.argv
#  @code
#  with KeepArgs() :
#    ...  
#  @endcode 
class KeepArgs(object) :
    """ Context manager to keep/preserve sys.argv
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
    """ Context manager to keep/preserve sys.argv
    >>> with keepArgs() :
    ...  
    """
    return KeepArgs()

# =============================================================================
## Return the path to an executable which would be run if the given <code>cmd</code> was called.
#  If no <code>cmd</code> would be called, return <code>None</code>.
#  - <code>mode</code> is a permission mask passed to <code>os.access()</code>,
#    by default determining if the file exists and executable.
#  - When no <code>path</code> is specified, the results of <code>os.environ()</code> are used,
#    returning either the <code>“PATH”</code> value or a fallback of <code>os.defpath</code>.
#  - copied from <code>shutil</cdde> module
def local_which ( cmd, mode=os.F_OK | os.X_OK, path=None):
    """ Given a command, mode, and a PATH string, return the path which
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
try : # =======================================================================
    # =========================================================================
    from shutil import which
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    which  = local_which 

# =============================================================================
## get the command
#  @code
#  >>> if cmd_exists ( 'epstopdf' ) : ... 
#  @endcode 
#  @see https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def cmd_exists ( command ) :
    """ Check the existence of certain command/executable
    >>> if cmd_exists ( 'epstopdf' ) : ...
    
    """
    return which ( command ) is not None


# =============================================================================
## split range into smaller chunks:
#  @code
#  for i in SplitRange ( 0 , 10000 , 200 ) :
#     for j in range (*i) :
#          ... 
#  @endcode 
class SplitRange(object) :
    """ Split range into smaller chunks:
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
    """ Split range into smaller chunks:
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
    """ Split range into `num` smaller chunks of approximate size 
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
        
# ======================================================================================
## use accumulae form itertools
from itertools import accumulate
# ====================================================================================
## use choiced from random 
choices = random.choices
# ====================================================================================
## Generate some random name of given name
#  @code
#  name = random_name ( 5 ) 
#  @endcode 
def random_name ( size = 6 , prefix = '' , suffix = '' ) :
    """ Generate some random name of given name 
    >>> name = random_name ( 5 )
    """
    assert 1 <= size , 'random_name: invalid size!'
    ## 
    first = random.choice  ( ascii_letters ) 
    if 1 == size : return prefix + first + suffix 
    ## 
    return prefix + first  + ''.join ( random.choices ( all_symbols , k = size - 1 ) ) + suffix 

# ========================================================================================
## generate some pseudo-random 6-symbol name from provided hash sources 
def short_hash_name ( size , name , *names ) :
    """ Generate some pseudo-random 6-symbol name from provided hash sources
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
    """ Generate the random string, that can be used as password or secret word
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
try : # =======================================================================
    # =========================================================================
    from more_itertools import chunked, divide 
    # =========================================================================    
except ImportError : # ========================================================
    # =========================================================================    
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
        """ Return first *n* items of the iterable as a list.
        
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
        """ Break *iterable* into lists of length *n*:
        
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
        """ Divide the elements from *iterable* into *n* parts, maintaining
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
## Collect data into fixed-length chunks or blocks"
def grouper ( iterable , n , fillvalue = None ):
    """ Collect data into fixed-length chunks or blocks """ 
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
    """ Create infinite or finite iterable from other iterable or no-iterable
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
    """ Calculate SHA512-checksum for the files
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
    """ Simple utility to check balanced parenthesis/brackets, etc...
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
## The simplest splitter of N-object into n groups :
#  @code
#  for i in spliter ( 10 , 3 ) : print ( i ) 
#  @endcode
def splitter ( N , n ) :
    """ The simplest splitter of N-objects into n groups
    >>> for i in spliter ( 10 , 3 ) : print ( i ) 
    """
    assert isinstance ( N , int ) and 0 <= N , "Invalid N"
    assert isinstance ( n , int ) and 1 <= n , "Invalid n"
    if   not N  : pass
    elif 1 == n : yield N
    else :
        a , b = divmod ( N , n )
        if   not a :
            for i in range ( b ) : yield 1
        elif not b :
            for i in range ( n ) : yield a
        else :
            for i in range ( b     ) : yield a + 1
            for i in range ( b , n ) : yield a 
# =========================================================================
## absract property decorator
#  @code
#  @absproperty
#  def A ( self ) : ...  
#  @endcode
def absproperty ( func ) :
    """ Abstract property
    @absproperty
    def A ( self ) : ...
    """
    return property ( abc.abstractmethod ( func ) ) 
    
# =============================================================================
if  ( 3 , 12 ) <= python_info : # =============================================
    # =========================================================================
    ## class property decorator
    #  @see https://stackoverflow.com/questions/76249636/class-properties-in-python-3-11
    class classprop :
        """ Class property decorator
        - see https://stackoverflow.com/questions/76249636/class-properties-in-python-3-11
        """
        def __init__(self, func):
            self.fget = func
        def __get__(self, instance, owner):
            return self.fget(owner)
    # =========================================================================
elif  ( 3 , 9 ) <= python_info : # ============================================
    # =========================================================================
    ## class property decorator
    #  @code
    #  @classprop
    #  def A ( cls ) : ...  
    #  @endcode
    def classprop ( func ) :
        """ Class property
        @classprop
        def A ( cls ) : ...
        """
        return classmethod ( property ( func ) ) 
    # =========================================================================
else : 
    # =========================================================================
    ## class @classproperty
    #  class property decorator  (copied and simplified from astropy)
    #  @code
    #  @classprop
    #  def A ( cls ) : ...  
    #  @endcode
    class classprop(property):
        """ Class property
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

# ============================================================================
## @class AttrGetter
#  simple class to bypass <code>operator.attrgetter</code> that
#  has some problem with serialization for multiprocessing
class AttrGetter(object):
    """ Simple class to bypass `operator.attrgetter` that
    has some problem with serializaton for multiprocessing
    """
    def __init__ ( self , *attributes ) :
        self.__attributes = attributes 
    def __call__ ( self , obj ) :
        getter = operator.attrgetter( *self.__attributes )
        return getter ( obj )

    @property 
    def attributes ( self ) :
        """`attributes': the actual attributes
        """
        return self.__attributes
    # print attributes 
    def __str__ ( self ) : return ','.join ( self.__attributes )
    __repr__ = __str__ 

# =============================================================================
## @class NumCalls
#  Count a number of  times a callable object is invoked
class NumCalls (object):
    """ Count a number of  times a callable object is invoked"""
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


## @class Singleton
#  Simple metaclass for the singleton
#  @code
#  class MyClass(metaclass = Singleton ) :
#      ....
#  a = MyClass()
#  b = MyClass()
#  print ( a is b ) 
#  endcode 
class Singleton(type) :
    """ Simple metaclass for the singleton:
    >>> class MyClass(metaclass = Singleton ) :
        ...
    >>> a = MyClass()
    >>> b = MyClass()
    >>> print ( a is b ) 
    """
    _INSTANCES = {}
    def __call__ ( klass , *args , **kwargs):
        """ Singleton """
        cls = klass._INSTANCES.get ( klass , None )
        if cls is None :
            ## create it:  
            cls = super(Singleton, klass).__call__ ( *args , **kwargs )
            klass._INSTANCES [ klass ] = cls
        else:
            ## Re-initialize it 
            cls.__init__ ( *args , **kwargs )

        return cls 


# ==============================================================================
## Copy file with the progress
#  @code
#  copy_with_progress ( 'inputfilename.ext' , 'outputfilename.ext' ) 
#  @endcode
def copy_with_progress ( source  , destination ) :
    """ Copy file with progress
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


# =============================================================================

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    logger.info ( 80*'*' ) 
            
# =============================================================================
##                                                                      The END 
# =============================================================================
