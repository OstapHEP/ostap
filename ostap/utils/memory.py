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
    'virtualMemory'  , # context manager to count virtual memory increase 
    'memory'         , # ditto
    'Memory'         , # ditto
    'memory_usage'   , # report current memory usage 
    )
# =============================================================================
import os  
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.memory' )
else                       : logger = getLogger( __name__             )
del getLogger
# =============================================================================
from sys import version_info as python_version 
if   2 < python_version.major : LONG = int
else                          : LONG = long 
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
        try : 
            with open ( proc , 'r' ) as p :
                for l in  p : return LONG ( l.split(' ')[22] )/1024./1024
        except:
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
#     at the exit it prints the chaneg in virtual memory
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
    from ostap.logger.logger import getLogger 
    _logger = getLogger( 'ostap.utils.memory' )
    del getLogger
    
    _printed = False
    def __init__  ( self , name = '' , logger = None , format = 'Memory %-18s %+.1fMB/[%.2fGB]') :
        self.name   = name
        self.logger = logger if logger else self._logger 
        self.format = format
        self._proc  = None 
        if not psutil :
            ## if not self._printed :
            ##    self.logger.warning('Memory:"psutil" module is not available, "/proc/[pid]/stat"-based replacement is in use')
            ##    self.__class__._printed = True
            self.proc = '/procs/%d/stat' %  os.getpid()            
    def __enter__ ( self ) :
        self.memory = memory_usage ( self._proc )
        return self 
    def __exit__  ( self, *_ ) :

        current     = memory_usage ( self._proc )
        self.delta  = current - self.memory
        try :
            message = self.format                    % ( self.name , self.delta , current / 1024. ) 
        except TypeError :
            message = 'Memory %-18s %+.1fMB/[%.2fGB]'% ( self.name , self.delta , current / 1024. )

        self.logger.info ( message )
 
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
