#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some simple but useful utilities for 
#   - timing
#   - memory
#   - profiling
#   - ... 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
#  
# =============================================================================
"""Module with some simple but useful utilities for Ostap
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    'with_ipython'  , ## do we run IPython ? 
    'isatty'        , ## is the stream ``isatty'' ?
    'terminal_size' , ## get the size of terminal cosole 
    )
# =============================================================================
import sys,os 
# =============================================================================
## is sys.stdout attached to terminal or not  ?
#  @code
#  stream = ...
#  if isatty( stream ) : print('Teminal!')
#  @endcode 
def isatty ( stream = None ) :
    """Is the stream is attached to terminal?
    >>> stream = ...
    >>> if isatty( stream ) : print('Teminal!')
    >>> if isatty() : print('stdout is terminal!')
    """
    if not stream : stream = sys.stdout
    #
    try :
        return sys.stdout.isatty()
    except : pass 
    #
    try :
        return os.isatty ( sys.stdout.fileno() ) 
    except : pass
    #
    return False

# =============================================================================
## helper function that allows to detect running ipython
def with_ipython()  :
    """Helper function that allows to detect running ipython"""
    try :
        return __IPYTHON__
    except NameError :
        return False

# ============================================================================
## Get the terminal console size
def terminal_size():
    """Get the terminal console size
    >>> height , width = terminal_size() 
    """
    try :
        import fcntl, termios, struct
        th, tw, hp, wp = struct.unpack(
            'HHHH',fcntl.ioctl(0, termios.TIOCGWINSZ,
                               struct.pack('HHHH', 0, 0, 0, 0)))
        return th , tw  
    except :
        return 20 , 80
    



# =============================================================================
if __name__ == '__main__' :

    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.utils.basic' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================

