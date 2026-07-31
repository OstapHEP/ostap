#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/utils/singleton,py
#  Very simple singleton-metaclass
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2024-10-10
# =============================================================================
""" Very simple singleton-metaclass
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-10-10"
# =============================================================================
__all__     = (
    ##
    'Singleton'          , ## Metaclass for the singleton 
    ##
)
# =============================================================================
import abc,threading
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.utils.singletor' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
## @class Singleton
#  Simple metaclass for the singleton
#  @code
#  class MyClass(metaclass = Singleton ) :
#      ....
#  a = MyClass()
#  b = MyClass()
#  print ( a is b ) 
#  endcode 
class Singleton(abc.ABCMeta) :
    """ Simple metaclass for the singleton:
    >>> class MyClass(metaclass = Singleton ) :
        ...
    >>> a = MyClass()
    >>> b = MyClass()
    >>> print ( a is b ) 
    """
    _INSTANCES = {}
    _LOCK      = threading.Lock()
    # ==========================================================================
    def __call__ ( klass , *args , **kwargs ) :
        """ Singleton """
        ## already created ? 
        if klass not in klass._INSTANCES:
            with klass._LOCK :
                if klass not in klass._INSTANCES:
                    instance = super().__call__ ( *args , **kwargs )
                    klass._INSTANCES [ klass ] = instance
                    return instance
                
        return klass._INSTANCES [ klass ]


# =============================================================================
if '__main__' == __name__ :
    # =========================================================================
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    logger.info ( 80*'*' ) 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
