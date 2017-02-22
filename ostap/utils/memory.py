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
"""Module with some simple but useful utilities for momory profiling 
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
import os , sys ## attention here!!
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.memory' )
else                       : logger = getLogger( __name__             )
del getLogger
# =============================================================================
_psutil = False 

try :
    import psutil 
    _psutil = True
    
    # =========================================================================
    ## report current memory usage (in MB)
    #  ps-based version, slow...  :-(
    #  @code
    #  print memory_usage() 
    #  @endcode
    def memory_usage ():
        """Report current memory usage (in MB)
        (psutil-based version, fast and efficient)
        - see help(psutil)
        >>> print memory_usage()
        """
        process = psutil.Process(os.getpid())
        mem     = process.get_memory_info()[0] / float(2 ** 20)
        return mem
    
except ImportError :
    _psutil = False

    # 
    # import resource
    # def memory_usage ():
    #     rusage_denom = 1024.
    #     if sys.platform == 'darwin':
    #         # ... it seems that in OSX the output is different units ...
    #         rusage_denom = rusage_denom * rusage_denom
    #     mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    #     return mem

    # =========================================================================
    ## report current memory usage (in MB)
    #  @attention it is ps-based version, slow...  :-(
    #  @code
    #  print memory_usage() 
    #  @endcode    
    def memory_usage () :
        """Report current memory usage (in MB)
        (ATTENTION: ps-based version, very slow 
        >>> print memory_usage()
        """
        import subprocess
        out = subprocess.Popen( [ 'ps', 'v', '-p', str(os.getpid())],
                                stdout=subprocess.PIPE).communicate()[0].split(b'\n')
        vsz_index = out[0].split().index(b'RSS')
        mem = float(out[1].split()[vsz_index]) / 1024
        
        return mem

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
    """Simple class to evaluate the change in virtual memory
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
    _logger = getLogger( 'ostap.utils.utils' )
    del getLogger
    
    _printed = False 
    def __init__  ( self , name = '' , logger = None , format = 'Delta(memory) %-18s %.2fMB') :
        self.name   = name
        self.logger = logger if logger else self._logger 
        self.format = format
        if not _psutil :
            if not self._printed :
                self.logger.warning('Memory:"psutil" module is not available, "ps"-based replacement is in use')
                self.__class__._printed = True
                
    def __enter__ ( self ) :
        self.memory = memory_usage ()
        return self 
    def __exit__  ( self, *_ ) :
        self.delta  = memory_usage () - self.memory
        try :
            message = self.format          % ( self.name , self.delta ) 
        except TypeError :
            message = 'Delta(memory) %-18s %.2fMB'% ( self.name , self.delta )

        self.logger.info ( message )
 
# ============================================================================
## create the context manager to monitor the virtual memory increase  
def virtualMemory ( name = '' ) :
    """Create the context manager to monitor the virtual memory increase:
    
    >>> with memory('here...') :
    ...   whatever action is here
    at the exit it prints the change in virtual memory
          
    >>> with memory('here...') as m :
    ...   whatever action is here
    at the exit it prints the change in virtual memory
    
    >>> delta = m.delta    
    """
    return Memory( name )

## ditto 
memory = virtualMemory  ## ditto

# =============================================================================
if '__main__' == __name__ :

        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    with memory(), memory() , memory()  :
        logger.info ( 80*'*' ) 
    
# =============================================================================
# The END 
# =============================================================================
