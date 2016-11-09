#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file
#  Module with some simple but useful utilities
#   - timing
#   - memory
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
#  
# =============================================================================
"""Module with some simple but useful utilities"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    #
    'virtualMemory'  , ## context manager to count virtual memory increase 
    'memory'         , ## ditto
    'clocks'         , ## context manager to count clocks 
    'timing'         , ## context manager to count time 
    'timer'          , ## ditto
    'profiler'       , ## context manager to perform profiling
    #
    'Profiler'       , ## context manager to perform profiling
    ##
    'takeIt'         , ## take and later delete ...
    )
# =============================================================================
import ROOT,cppyy, time, os,sys ## attention here!!
cpp = cppyy.gbl
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.utils' )
else                       : logger = getLogger( __name__            )
del getLogger 
# =============================================================================
## @class Memory
#  Simple context manager to measure the virtual memory increase
#
#  @see System::virtualMemory
#  @code
#
#  with Memory() :
#     <whatever action is>
#     at the exit it prints the chaneg in virtual memory 
#  @endcode
#
# Or:
#
#  @code
#
#  with Memory() as Q :
#     <whatever action is>
#     at the exit it prints the chaneg in virtual memory
#
#  print Q.delta 
# 
#  @endcode
#
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-02-10
class Memory(object):
    """Simple class to evaluate the change in virtual memory
    to be used as context manager:
    
    >>> with Memory('here...') :
    ...     <whatever action is>
    at the exit it prints the change in virtual memory
    
    >>> with Memory('here...') as M :
    >>> <whatever action is>
    at the exit it prints the change in virtual memory
    
    >>> delta = M.delta    
    """
    _logger = logger 
    def __init__  ( self , name = '' , logger = None , format = 'Memory %-18s %.1fMB') :
        self.name   = name
        self.logger = logger if logger else self._logger 
        self.format = format
    def __enter__ ( self ) :
        self.memory = cpp.System.virtualMemory()
        return self 
    def __exit__  ( self, *_ ) :
        self.delta  = cpp.System.virtualMemory() - self.memory
        try :
            message = self.format          % ( self.name , self.delta ) 
        except TypeError :
            message = 'Memory %-18s %.1fMB'% ( self.name , self.delta )

        self.logger.info ( message )
 
# ============================================================================
## create the context manager to monitor the virtual memory increase  
def virtualMemory ( name = '' ) :
    """Create the context manager to monitor the virtual memory increase:
    
    >>> with memory('here...') :
    ...   <whatever action is>
    at the exit it prints the change in virtual memory
          
    >>> with memory('here...') as m :
    ...   <whatever action is>
    at the exit it prints the change in virtual memory
    
    >>> delta = m.delta    
    """
    return Memory( name )

## ditto 
memory = virtualMemory  ## ditto

# =============================================================================
## @class Clock
#  Smple context manager to measure the clock counts
#
#  @code
#
#  with Clock() :
#     <whatever action is>
#     at the exit it prints the clock counts 
#  @endcode
#
# Or:
#
#  @code
#
#  with Clock() as c :
#     <whatever action is>
#     at the exit it prints the clock counts 
#
#  print c.delta 
# 
#  @endcode
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-02-10                    
class Clock(object):
    """Simple context manager to measure the clock counts
    >>> with Clock() :
    ...  <whatever action is>
    at the exit it prints the clock counts 
    
    >>> with Clock() as c :
    ...  <whatever action is>
    at the exit it prints the clock counts 
    
    >>> print c.delta 
    """
    _logger = logger 
    def __init__  ( self , name = '' , logger = None , format = 'Clocks %-18s %s') :
        self.name   = name
        self.logger = logger if logger else self._logger 
        self.format = format
    def __enter__ ( self ) :
        self.clock = time.clock()
        return self 
    def __exit__  ( self, *_ ) :
        self.delta = time.clock() - self.clock
        try :
            message = self.format       % ( self.name , self.delta ) 
        except TypeError :
            message = 'Clocks %-18s %s' % ( self.name , self.delta )

        self.logger.info ( message )
        
# =============================================================================
## @class Timer
#  Simple context manager to measure the time 
#  @code
#
#  with Timer() :
#     <whatever action is>
#     at the exit it prints the time 
#  @endcode
#
# Or:
#
#  @code
#
#  with Timer() as t :
#     <whatever action is>
#     at the exit it prints the clock counts 
#
#  print ct.delta 
# 
#  @endcode
#
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-02-10
#
class Timer(object):
    """Simple context manager to measure the time
    
    >>> with Timer() :
    ...  <whatever action is>
    at the exit it prints the time 
    
    Or:
    
    >>> with Timer() as t :
    ...  <whatever action is>
    at the exit it prints the clock counts 
    
    >>> print ct.delta 
    """
    _logger = logger 
    def __init__  ( self , name = '' , logger = None , format = 'Timing %-18s %.3f' ) :
        self.name   = name
        self.logger = logger if logger else self._logger 
        self.format = format
    def __enter__ ( self ) :
        self.time = time.time()
        return self 
    def __exit__  ( self, *_ ) :
        self.delta = time.time() - self.time
        
        try :
            message = self.format       % ( self.name , self.delta ) 
        except TypeError :
            message = 'Timing %-18s %s' % ( self.name , self.delta )

        self.logger.info ( message )
            
# =============================================================================
## Simple context manager to measure the clock counts
#
#  @code
#
#  with clocks () :
#     <whatever action is>
#     at the exit it prints the clock counts
#
#  @endcode
#
# Or:
#
#  @code
#
#  with clocks () as c :
#     <whatever action is>
#     at the exist it prints the clock counts 
#
#  print c.delta 
# 
#  @endcode
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-02-10                    
def clocks ( name = '' ) :
    """Simple context manager to measure the clock counts 
    
    >>> with clocks () :
    ...   <whatever action is>
    at the exit it prints the clock counts 
    
    >>> with clocks () as c :
    ...   <whatever action is>
    at the exit it prints the clock counts 
    
    >>>print c.delta
    """
    return Clock ( name )

# =============================================================================
## Simple context manager to measure the time
#
#  @code
#
#  with timer () :
#     <whatever action is>
#     at the exit it prints the time
#
#  @endcode
#
# Or: 
#
#  @code
#
#  with timer () as t :
#     <whatever action is>
#     at the exit it prints the clock counts 
#
#  print t.delta 
# 
#  @endcode
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-02-10                    
def timing ( name = '' , logger = None ) :
    """Simple context manager to measure the clock counts 
    
    >>> with timing () :
    ...   <whatever action is>
    at the exit it prints the clock counts 
    
    >>> with timing () as c :
    ...   <whatever action is>
    at the exit it prints the clock counts 
    
    >>> print c.delta
    """
    return Timer ( name , logger )

# =============================================================================
## ditto 
timer = timing   # ditto


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
    """  Very simple profiler, based on cProfile module
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
    """ Very simple profiler, based on cProfile module
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
#  @endocode 
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
## @class TakeIt#
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
    """ Take some object, keep it and delete at the exit
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
    """ Get the actual file name from file descriptor 
    The actual code is copied from http://stackoverflow.com/a/13624412 
    """
    names = []
    for fd in fds:
        names.append(os.readlink('/proc/self/fd/%d' % fd))
    return names


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.logger.line import line 
    logger.info ( __file__  + '\n' + line  ) 
    logger.info ( 80*'*'   )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 
    
# =============================================================================
# The END 
# =============================================================================
