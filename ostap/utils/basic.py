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
"""Module with some simple but useful utilities for
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    'with_ipython' , ## do we run IPython ? 
    'isatty'       , ## is the stream ``isatty'' ? 
    )
# =============================================================================
import sys,os 
# =============================================================================
## is sys.stdout attached to terminal or not  ?
#  @code
#  stream = ...
#  if isatty( stream ) : print 'Teminal!'
#  @endcode 
def isatty ( stream = None ) :
    """Is the stream is attached to terminal?
    >>> stream = ...
    >>> if isatty( stream ) : print 'Teminal!'
    >>> if isatty() : print 'stdout is terminal!'
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

# =============================================================================
if __name__ == '__main__' :

    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.utils.basic' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================

