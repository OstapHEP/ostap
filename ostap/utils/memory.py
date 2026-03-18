#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some simple but useful utilities for memory profiling 
#  - It is recommended to install psutil module
#  @see http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler/
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
# =============================================================================
""" Module with some simple but useful utilities for momory profiling 
- It is recommended to install psutil module (e.g. from pip)
see https://github.com/giampaolo/psutil
see http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler/
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    'virtualMemory'  , ## context manager to count virtual memory increase 
    'memory'         , ## ditto
    'Memory'         , ## ditto
    'memory_usage'   , ## report current memory usage 
    'delta_ram'      , ## symbol for delta-RAM
    )
# =============================================================================
from   ostap.logger.symbols import ram          as ram_symbol
from   ostap.logger.symbols import delta_symbol, sum_symbol
import os
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.memory' )
else                       : logger = getLogger( __name__             )
del getLogger
# =============================================================================
ram_symbol = ram_symbol + ' ' if ram_symbol else ''
# =============================================================================
delta_ram  = '%s%s' % ( delta_symbol , ram_symbol ) 
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import psutil 
    # =========================================================================
    ## report current memory usage (in MB)
    #  ps-based version, slow...  :-(
    #  @code
    #  print memory_usage() 
    #  @endcode
    def memory_usage ( *args ):
        """ Report current memory usage (in MB)
        (psutil-based version, likely fast and efficient)
        - see help(psutil)
        >>> print memory_usage()
        """
        process = psutil.Process ( os.getpid () )
        mem     = process.memory_info()[0] / float( 2 ** 20 )
        for p in process.children ( recursive = True ) :
            mem += p.memory_info()[0] / float( 2 ** 20 )
        return mem
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================    
    psutil = None 
    # =========================================================================
    ## report current memory usage (in MB)
    #  @attention it is ps-based version, it can be slow...  :-(
    #  @code
    #  print memory_usage() 
    #  @endcode    
    def memory_usage ( proc = None ) :
        """ Report current memory usage (in MB)
        """
        if not proc : 
            import os 
            proc = '/proc/%d/stat' % os.getpid()
        # ====================================================================
        try : # ==============================================================
            # ================================================================
            with open ( proc , 'r' ) as p :
                for l in  p : return int ( l.split(' ')[22] )/1024./1024
            # ================================================================
        except: # ============================================================
            # ================================================================
            return -1
# ============================================================================
        
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
    return ma * 1.0 / mu

# =============================================================================
## @class Memory
#  Simple context manager to measure the virtual memory increase
#
#  @see System::virtualMemory
#  @code
#
#  with Memory() :
#     whatever action is here 
#     at the exit it prints the chaneg in virtual memory 
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
    _printed = False
    
    def __init__  ( self , name = '' , logger = None , format = 'Memory %%-30s %s=%%+.2fMB %s=%%.2fGB' % ( delta_ram , sum_symbol ) ) :
        self.name   = name
        self.logger = logger if logger else self._logger 
        self.format = format
        self._proc  = None 
        if not psutil :
            ## if not self._printed :
            ##    self.logger.warning('Memory:"psutil" module is not available, "/proc/[pid]/stat"-based replacement is in use')
            ##    self.__class__._printed = True
            self.proc = '/procs/%d/stat' %  os.getpid()
            
        self.__memory = self.current 
        self.__delta  = None

    def __enter__ ( self ) :
        
        self.__memory  = self.current 
        self.__delta = None  
        return self
    
    def __exit__  ( self, *_ ) :

        current      = self.current 
        self.__delta = current - self.memory
        total        = current / 1024         ## in GB 
        # ====================================================================
        try : # ==============================================================
            # ================================================================
            message = self.format % ( self.name , self.delta , total )
            # ================================================================
        except TypeError : # =================================================
            # ================================================================
            message = 'Memory %-30s %s=%+.2fMB %s=%.2fGB'% ( self.name , delta_ram , self.delta , sum_symbol , total  )
            # ================================================================
        self.logger.info ( ram_symbol + message )

    ## delta-memory
    @property 
    def delta ( self ) :
        """`delta`: delta-memory [MB] between enter&exit or None
        """
        return self.__delta 

    ## current memory
    @property
    def current ( self ) :
        """`current` : current memory usage [MB]
        """
        return memory_usage ( self._proc )

    ## memory at the start of measurement
    @property
    def memory ( self ) :
        """`memory` : initial/start memory usage [MB] 
        """
        return self.__memory 
        
# ============================================================================
## create the context manager to monitor the virtual memory increase  
def virtualMemory ( name = '' , logger = None ) :
    """ Create the context manager to monitor the virtual memory increase:
    
    >>> with memory('here...') :
    ...   whatever action is here
    at the exit it prints the change in virtual memory
          
    >>> with memory('here...') as m :
    ...   whatever action is here
    at the exit it prints the change in virtual memory
    
    >>> delta = m.delta    
    """
    return Memory( name , logger = logger )

## ditto 
memory = virtualMemory  ## ditto

# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    with memory(), memory() , memory()  :
        logger.info ( 80*'*' ) 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
