#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/clsgetter.py
#  Property-getter for classmethod 
#  @code
#  class A(object) : pass
#  A.x = classgetter( lambda cls : 'QUQU! )
#  @endcode
#  or as decorator
#  @code
#  class A(object) :
#      @classgetter
#      def x ( cls ) :
#         return "QUQU!"
#  @endcode
#
#  @see https://stackoverflow.com/questions/128573/using-property-on-classmethods/13624858#13624858
#
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
"""Property-getter for classmethod 
>>> class A(object) : pass
>>> A.x = classgetter( lambda cls : 'QUQU! )
-   or as decorator
>>> class A(object) :
...     @classgetter
...     def x( cls ) :
...         return 'QUQU!'
- see https://stackoverflow.com/questions/128573/using-property-on-classmethods/13624858#13624858
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2020-04-24"
__version__ = ""
# =============================================================================
__all__     = (
    'classgetter' , ## classmethof proeprty-=getter 
    )
# =============================================================================
from ostap.math.base        import isequal , iszero , std, Ostap
from ostap.core.ostap_types import num_types, integer_types
# =============================================================================

# ============================================================================
## Property-getter for classmethod 
#  @code
#  class A(object) : pass
#  A.x = classgetter( lambda cls : 'QUQU! )
#  @endcode
#  or as decorator
#  @code
#  class A(object) :
#      @classgetter
#      def x ( cls ) :
#         return "QUQU!"
#  @endcode 
#  @see https://stackoverflow.com/questions/128573/using-property-on-classmethods/13624858#13624858
class classgetter(object):
    """Property-getter for classmethod 
    >>> class A(object) : pass
    >>> A.x = classgetter( lambda cls : 'QUQU! )
    -   or as decorator
    >>> class A(object) :
    ...     @classgetter
    ...     def x( cls ) :
    ...         return 'QUQU!'
    - see https://stackoverflow.com/questions/128573/using-property-on-classmethods/13624858#13624858
    """
    def __init__(self, fget):
        self.fget = fget
    def __get__(self, owner_self, owner_cls):
        return self.fget(owner_cls)
    
# ============================================================================
if '__main__' == __name__  :
    
    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.utils.clsgetter' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
