#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================
## @file
#  Module with some simple but useful utilities for memory profiling 
#  - It is recommended to install psutil module
#  @see http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler/
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
# ==========================================================================
""" Module with some simple but useful utilities for momory profiling 
- It is recommended to install psutil module (e.g. from pip)
see https://github.com/giampaolo/psutil
see http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler/
"""
# ==========================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# ==========================================================================
__all__     = (
    'virtualMemory'  , ## context manager to count virtual memory increase 
    'memory'         , ## ditto
    'Memory'         , ## ditto
    'memory_usage'   , ## report current memory usage 
    'delta_ram'      , ## symbol for delta-RAM
    'PeakMemory'     , ## context manager to check the peak memory
    'peak_memory'    , ## ditto 
    )
# ==========================================================================
from   ostap.logger.symbols import ( delta_symbol , sum_symbol ,
                                     ram      as ram_symbol    , 
                                     enough   as enough_symbol ,
                                     mountain as peak_symbol   ) 
import os, psutil, time, threading  
# ==========================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.memory' )
else                       : logger = getLogger( __name__             )
_logger_t = type ( logger )
del getLogger
# ==========================================================================
ram_symbol = ram_symbol + ' ' if ram_symbol else ''
delta_ram  = '%s%s' % ( delta_symbol , ram_symbol )
peak_ram   = '%s%s' % ( peak_symbol  , ram_symbol ) 
# ==========================================================================
## Fast & Cross-Platform Memory Usage Reporter
# ==========================================================================
## report current memory usage (in MB)
#  @attention it is ps-based version, it can be slow...  :-(
#  @code
#  print memory_usage() 
#  @endcode
def _memory_usage ( proc = None ) :
    """ Fallback Linux-only procfs RSS memory reader (in MB).
    """
    # =====================================================================
    try : # ===============================================================
        # =================================================================
        if not proc : proc = f"/proc/{os.getpid()}/stat"
        with open ( proc , "r" ) as p :
            line = p.readline ()
            # Field 23 in /proc/[pid]/stat contains RSS in pages (typically 4096 bytes each)
            pages = int ( line.split () [ 23 ] )
            return ( pages * 4096 ) / ( 1024 * 1024 )
        # =================================================================
    except Exception : # ==================================================
        # =================================================================
        return -1.0

# =========================================================================
# Cache the current Process instance once at import to minimize overhead
PSUTIL_CURRENT_PROCESS = psutil.Process ( os.getpid () )
# =========================================================================
## Report current process RSS memory usage in MB.
#  Uses a cached psutil Process object for maximum speed.
def memory_usage ( include_children = True  ) :
    """ Report current process RSS memory usage in MB.
    Uses a cached psutil Process object for maximum speed.
    
    >>> print ( memory_usage () )
    """
    # =====================================================================
    try : # ===============================================================
        # =================================================================
        # Fetch RSS (Resident Set Size) in bytes directly from kernel
        mem_bytes = PSUTIL_CURRENT_PROCESS.memory_info ().rss
        # =================================================================
        # Check for sub-processes 
        if include_children : # ===========================================
            # =============================================================
            for child in PSUTIL_CURRENT_PROCESS.children ( recursive = True ) :
                # =========================================================
                try : # ===================================================
                    # =====================================================
                    mem_bytes += child.memory_info ().rss
                # =====================================================
                except ( psutil.NoSuchProcess , psutil.AccessDenied ) : # =
                    # =====================================================
                    pass
                
        return mem_bytes / ( 1024 * 1024 )
        # ==================================================================
    except Exception : # ===================================================
        # ==================================================================
        return _memory_usage ()
    
# ============================================================================
## Get the available memory (in MB) 
def memory_available () :
    """ Get the available memory (in MB) 
    """
    if not psutil : return -1
    vm = psutil.virtual_memory()
    return vm.available / float ( 2 ** 20 )
# =============================================================================
## Get the ratio of available memory to used
def memory_enough () :
    """ Get the ratio of available memory to used
    """
    ma = memory_available () 
    mu = memory_usage     ()
    return ma * 1.0 / mu if 0 < mu and 0 < ma else -1 

# ============================================================================
## format for memory prints 
MEMORY_FORMAT = 'Memory %%-20s %s=%%+.3fMB %s=%%.3fGB %s%%.1f' % ( delta_ram , sum_symbol , enough_symbol )
# =============================================================================
## @class Memory
#  Simple context manager to measure the virtual memory increase
#
#  @see System::virtualMemory
#  @code
#
#  with Memory() :
#     whatever action is here 
#     at the exit it prints the change in virtual memory 
#  @endcode
#
# Or:
#
#  @code
#
#  with Memory() as Q :
#     whatever action is here 
#     at the exit it prints the change in virtual memory
#
#  print Q.delta 
# 
#  @endcode
#
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-02-10
class Memory(object):
    """ Simple class to evaluate the change in virtual memory
    to be used as context manager:
    
    >>> with Memory('here...') :
    ...     whatever action is here
    at the exit it prints the change in virtual memory
    
    >>> with Memory('here...') as M :
    ...     whatever action is here
    at the exit it prints the change in virtual memory
    
    >>> delta = M.delta    
    """
    _logger  = logger     
    
    def __init__  ( self ,
                    name   = ''   , * , 
                    logger = ''   ,
                    format = MEMORY_FORMAT ) :
        
        
        if   logger and isinstance ( logger , _logger_t ) : self.logger = logger.info
        elif logger and callable   ( logger )             : self.logger = logger
        elif logger is None                               : self.logger = None 
        else                                              : self.logger = self._logger 

        self.name   = name

        if not format : format = MEMORY_FORMAT

        self.format = format
        self._proc  = None
        
        self.__memory = self.current 
        self.__delta  = None

    # =======================================================================
    ## Context manager: ENTER 
    def __enter__ ( self ) :
        """ Context manager: ENTER """ 
        self.__memory = self.current 
        self.__delta  = None  
        return self

    # =========================================================================
    ## Context maanger: EXIT 
    def __exit__  ( self, *_ ) :
        """ Context manager: EXIT """ 
        current      = self.current 
        self.__delta = current - self.memory
        total        = current / 1024         ## in GB
        enough       = max ( 0 , memory_enough () ) 
        
        # =====================================================================
        if self.logger  : # ===================================================
            # =================================================================
            try : # ===========================================================
                # =============================================================
                message = self.format % ( self.name , self.delta , total , enough )
                # =============================================================
            except TypeError : # ==============================================
                # =============================================================
                format  = MEMORY_FORMAT 
                message = format % ( self.name , self.delta , total , enough )
                # =============================================================                
            self.logger ( ram_symbol + message )

    # =========================================================================
    ## delta-memory
    @property 
    def delta ( self ) :
        """`delta`: delta-memory [MB] between enter&exit or None
        """
        return self.__delta 

    # =========================================================================
    ## current memory
    @property
    def current ( self ) :
        """`current` : current memory usage [MB]
        """
        return memory_usage ( True )
    
    # =========================================================================
    ## memory at the start of measurement
    @property
    def memory ( self ) :
        """`memory` : initial/start memory usage [MB] 
        """
        return self.__memory 
        
# ============================================================================
## create the context manager to monitor the virtual memory increase  
def virtualMemory ( name = '' , **kwargs ) :
    """ Create the context manager to monitor the virtual memory increase:
    
    >>> with memory('here...') :
    ...   whatever action is here
    at the exit it prints the change in virtual memory
          
    >>> with memory('here...') as m :
    ...   whatever action is here
    at the exit it prints the change in virtual memory
    
    >>> delta = m.delta    
    """
    return Memory( name , **kwargs  )

## ditto 
memory = virtualMemory  ## ditto

# =============================================================================
## format for memory prints 
PEAK_MEMORY_FORMAT = 'Peak memory %%-20s %s=%%.3fGB %s=%%.3fGB' % ( peak_symbol , sum_symbol )
# =============================================================================
## @class PeakMemory
#  Context manager reprting the peak memory consumption
#  @attention : there is some overehead! 
#  @code
#  with PeakMemory ( "my code block" ) :
#  ... 
#  @endcode
class PeakMemory ( object ) :
    """ Context manager reporting the peak memory consumption
    - attention : there is some overhead! 
    >>> with PeakMemory ( "my code block" ) :
    ...  <some code here> 
    """
    _logger  = logger     
    
    def __init__ ( self ,
                   name     = ''                 , * ,
                   logger   = ''                 ,
                   format   = PEAK_MEMORY_FORMAT ,                    
                   interval = 0.005 ) :
        
        if   logger and isinstance ( logger , _logger_t ) : self.logger = logger.info
        elif logger and callable   ( logger )             : self.logger = logger
        elif logger is None                               : self.logger = None 
        else                                              : self.logger = self._logger 
        
        if not format  : format = PEAK_MEMORY_FORMAT 
        
        self.name             = name
        self.format           = format 
        self.__interval       = interval
        self.__process        = PSUTIL_CURRENT_PROCESS 
        self.__peak_bytes     = 0
        self.__peak_GB        = 0.0
        self.__stop_event     = threading.Event ()
        self.__monitor_thread = None
    
    # =========================================================================
    def _poll_memory ( self ) :
        # =====================================================================
        while not self.__stop_event.is_set () : # =============================
            # =================================================================
            try : # ===========================================================
                # =============================================================
                current_bytes = self.__process.memory_info ().rss
                if current_bytes > self.__peak_bytes :
                    self.__peak_bytes = current_bytes
                # =============================================================
            except Exception : # ==============================================
                # =============================================================
                pass
            time.sleep ( self.__interval )

    # =========================================================================
    ## Context manager: ENTER
    def __enter__ ( self ) :
        """ Context manger: ENTER """
        self.__peak_bytes = self.__process.memory_info ().rss
        self.__stop_event.clear ()
        self.__monitor_thread = threading.Thread ( target = self._poll_memory , daemon = True )
        self.__monitor_thread.start ()
        return self

    # =========================================================================
    ## Context manger: EXIT 
    def __exit__ ( self , *_ ) :
        """ Context manger: EXIT """
        
        self.__stop_event.set ()
        if self.__monitor_thread is not None :
            self.__monitor_thread.join ()
            
        self.__peak_MB = self.__peak_bytes / ( 1024 * 1024 )
        total = self.current               /   1024 

        peak_GB = self.__peak_MB / 1024 
        # =====================================================================
        if self.logger : # ====================================================
            # =================================================================
            try : # ===========================================================
                # =============================================================
                message = self.format % ( self.name , peak_GB , total )
                # =============================================================
            except TypeError : # ==============================================
                # =============================================================
                format  = PEAK_MEMORY_FORMAT 
                message = format % ( self.name , peak_GB , total )
                # =============================================================                
            self.logger ( peak_ram + message )

        return False

    # =========================================================================
    ## current memory
    @property
    def current ( self ) :
        """`current` : current memory usage [MB]
        """
        return memory_usage ( True )

    # =========================================================================
    def peak ( self ) :
        """`peak` : peak memory usage (in MB)"""
        return self.__peak_MB 


# ===============================================================================
## Context manager repщrting the peak memory consumption
#  @attention : there is some overehead! 
#  @code
#  with peak_memory ( "my code block" ) :
#  ... 
#  @endcode
def peak_memory ( name = '' , **kwargs ) : 
    """ Context manager reporting the peak memory consumption
    - attention : there is some overehead! 
    >>> with peak_memory ( "my code block" ) :
    ...  <some code here> 
    """
    return PeakMemory ( name , **kwargs )

# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    with peak_memory() , memory(), memory() , memory()  :
        logger.info ( 80*'*' ) 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
