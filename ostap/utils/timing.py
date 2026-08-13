#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file
# Set of useful utilities
# =============================================================================
"""Set of useful utilisties for timing
"""
# =============================================================================
__all__     = (
    ## 
    'timing'         , ## context manager to count time 
    'timer'          , ## ditto
    'Timer'          , ## context manager to count time 
    ##
    'Wait'           , ## conitext manager to wait soem tiem bvefore and/or after action
    'wait'           , ## conitext manager to wait soem tiem bvefore and/or after action 
    ##
   )
# =============================================================================
from   ostap.logger.symbols import clock as clock_symbol
import time 
# =============================================================================
from   ostap.logger.logger  import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.utils.timing' )
else                       : logger = getLogger ( __name__             )
_logger_t = type ( logger )
del getLogger
# =============================================================================
clock_symbol = clock_symbol + ' ' if clock_symbol else '' 
# =============================================================================
## format for timing prints 
TIMING_FORMAT = 'Timing %-20s | WALL: %.3fs | CPU: %.3fs | CPU Load %.1f%%' 
# =============================================================================
## @class Timer
#  Simple context manager to measure the time 
#  @code
#  with Timer() :
#     ... whatever action is here 
#     ... at the exit it prints the time 
#  @endcode
# Or:
#  @code
#  with Timer() as t :
#     ... whatever action is here 
#     ... at the exit it prints the clock counts 
#  print t.delta
#  @endcode
#  Print to specified logger :
#  @code
#  logger = ...
#  wtih Timer ( logger = logger ) : ... 
#  @endcode 
#  Print to specified logger level  :
#  @code
#  logger = ...
#  wtih Timer ( logger = logger.error  ) : ... 
#  @endcode
#  no print at all
#  @code
#  logger = ...
#  wtih Timer ( logger = lambda *s : '' ) : ... 
#  @endcode
class Timer(object):
    """ Simple context manager to measure the time
    
    >>> with Timer() :
    ...  whatever action is
    at the exit it prints the time 
    
    Or:
    
    >>> with Timer() as t :
    ...  whatever action is
    at the exit it prints the clock counts 
    
    >>> print t.delta

    - Print to specified logger :
    >>> logger = ...
    >>> with Timer ( logger = logger ) : ... 

    - Print to specified logger level  :
    >>> logger = ...
    >>> with Timer ( logger = logger.error  ) : ...
    
    - no print at all
    >>> logger = ...
    >>> with Timer ( logger = lambda *s : '' ) : ... 
    """
    __logger = logger.info
    ##
    def __init__  ( self                   ,
                    name   = ''            , * , 
                    logger = ''            ,
                    format = TIMING_FORMAT ,
                    start  = ''                   ) : 
        
        self.__name   = name

        if   logger and isinstance ( logger , _logger_t ) : self.logger = logger.info
        elif logger and callable   ( logger )             : self.logger = logger
        elif logger is None                               : self.logger = None 
        else                                              : self.logger = self.__logger 

        
        if not format : format = TIMING_FORMAT 
        
        self.format        = format
        
        if   start              : self.start_message = start
        elif name and not start : self.start_message = 'Start  %s' % name
        else                    : self.start_message = ''
        
        self.__start_wall = 0
        self.__start_cpu  = 0
        self.__delta_wall = 0
        self.__delta_cpu  = 0
        
    # ==========================================================================
    ## context manager: ENTER 
    def __enter__ ( self ) :
        """ Context manager: ENTER
        """
        if self.start_message and self.logger :
            self.logger ( clock_symbol + self.start_message )
            
        self.__delta_wall = 0
        
        self.__delta_cpu  = 0        
        self.__start_wall = time.perf_counter ()
        self.__start_cpu  = time.process_time () 
        
        return self
    
    # ==========================================================================
    ## context manager: EXIT 
    def __exit__  ( self, *_ ) :
        """ Context manager: EXIT
        """

        delta_wall = time.perf_counter () - self.__start_wall
        delta_cpu  = time.process_time () - self.__start_cpu
        
        # Calculate CPU utilization (can exceed 100% with multi-threading)
        cpu_util = delta_cpu / delta_wall * 100 if 0 < delta_wall else 0 

        self.__delta_wall = delta_wall 
        self.__delta_cpu  = delta_cpu
        

        # =====================================================================
        if self.logger : # ====================================================
            # =================================================================
            try : # ===========================================================
                # =============================================================
                message = self.format % ( self.__name , delta_wall , delta_cpu , cpu_util )
                # =============================================================
            except TypeError : # ==============================================
                # =============================================================
                format  = TIMING_FORMAT 
                message = format % ( self.__name , delta_wall , delta_cpu , cpu_util  )

            ## print it! 
            self.logger ( clock_symbol + message )
            
    @property
    def name ( self ) :
        """`name' : Timer name"""
        return self.__name

    @property
    def start_cpu  ( selt ) :
        """`start_cpu`: start CPU timer"""
        return self.__start_cpu
    
    @property
    def start_wall ( selt ) :
        """`start_wall`: start WALL timer"""
        return self.__start_wall 
    
    @property
    def delta_cpu ( self ) :
        """`delta_cpu' : stop - start for CPU timer """
        return max ( 0 , self.__delta_cpu ) 
    
    @property
    def delta_wall ( self ) :
        """`delta_wall' : stop - start for WALL timer """
        return max ( 0 , self.__delta_wall ) 

    @property
    def start ( self ) :
        """`start` : start Timer(s)"""
        return 0.5 * ( self.__start_wall + self.__start__cpu )

    @property
    def delta  ( self ) :
        """`delta` : stop - start for CPU timer"""
        return self.delta_cpu 

    # =========================================================================
    ## Calculate CPU utilization [%] (can exceed 100% with multi-threading)
    @property
    def cpu_utilization ( self ) :
        """`cpu_utlization` : CPU utilization [%] (can exceed 100% with multi-threading)
        """
        if self.__delta_wall <= 0 : return 0
        cpu_util = self.__delta_cpu / self.__delta_wall * 100
        return cpu_util

# =============================================================================
## Simple context manager to measure the time
#
#  @code
#  with timer () :
#     ... whatever action is here 
#     ... at the exit it prints the time
#  @endcode
#
# Or: 
#
#  @code
#  with timer () as t :
#     ... whatever action is here 
#     ... at the exit it prints the clock counts 
#  print t.delta 
#  @endcode
def timing ( name   = ''   , 
             logger = None , * , 
             format = 'Timing %-18s | WALL: %.3fs | CPU: %.3fs | CPU Load %.1f%%' , 
             **kwargs  ) :
    """ Simple context manager to measure the clock counts 
    
    >>> with timing () :
    ...   whatever action is here
    at the exit it prints the clock counts 
    
    >>> with timing () as c :
    ...   whatever action is here 
    at the exit it prints the clock counts 
    
    >>> print c.delta
    """
    return Timer ( name = name , logger = logger , format = format , **kwargs ) 

## ditto 
timer = timing   # ditto


# =============================================================================
## context manager that invokes <code>time.sleep</code> before and after action
#  @code
#  with Wait ( after = 5 , before = 0 ) :
#  ...
#  @endcode
class Wait(object):
    """ Context manager that invokes <code>time.sleep</code> before and after action
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
    """ Context manager that invokes <code>time.sleep</code> before and after action
    >>> with wait ( after = 5 , before = 0 ) :
    >>> ...
    """    
    return Wait ( after = after , before = before )

# =============================================================================
if '__main__' == __name__ :
    
    with timer ( logger = logger ) :
        from ostap.utils.docme import docme
        docme ( __name__ , logger = logger )
        
# =============================================================================
##                                                                     The END 
# =============================================================================
