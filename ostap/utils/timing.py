
__all__     = [
    'clocks'         , ## context manager to count clocks 
    'timing'         , ## context manager to count time 
    'timer'          , ## ditto
    ]

from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils' )
else                       : logger = getLogger( __name__ )
del getLogger 
from ostap.logger.logger import  logVerbose,  logDebug, logInfo, logWarning, logError

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

## ditto 
timer = timing   # ditto
