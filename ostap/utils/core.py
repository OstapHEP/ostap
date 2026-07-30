#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with simple core utilities 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2026-07-30
# =============================================================================
""" Module with some simple core utilities for Ostap 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2026-07-20"
# =============================================================================
__all__     = (
    ##
    'is_cpp_class' , ## true for pure C++ classes
    'typename'     , ## the typename of the object
    ##
    'isfunction'   , ## is it a function (or lambda) ?
    'islambda'     , ## is it a lambda?
    'ismethod'     , ## is it a method? (taked frmo inspecr) 
    # =========================================================================
) # ===========================================================================
# =============================================================================
from inspect import ismethod
from types   import FunctionType, LambdaType
import cppyy, functools
# =============================================================================
## Determine whether `cls` is a pure C++ class/struct registered in cppyy/ROOT.
#  
#  @param cls: Object/type to inspect.
#  @return True if `cls` is a native C++ class, False otherwise.
def is_cpp_class ( cls ) -> bool :
    """ Determine whether `cls` is a pure C++ class/struct registered in cppyy/ROOT.
    
    :param cls: Object/type to inspect.
    :return: True if `cls` is a native C++ class, False otherwise.
    """
    if not isinstance ( cls, type ) : return False

    # =======================================================================
    # 1. In ROOT >= 6.26, __cpp_name__ is a string property containing the full C++ path
    # (e.g., "Ostap::Kinematics::Dalitz" or "TH1F")
    # =======================================================================
    cpp_name = getattr ( cls, '__cpp_name__'  , None )
    if not cpp_name or not isinstance ( cpp_name , str ) : return False

    # =======================================================================
    # Strip leading global scope operator '::' if present
    # =======================================================================
    clean_name = cpp_name[2:] if cpp_name.startswith('::') else cpp_name

    # =======================================================================
    # 2. Resolve the full namespace path inside C++ global scope (cppyy.gbl)
    # =======================================================================
    try: # ==================================================================
        # ===================================================================
        resolved_cls = functools.reduce ( getattr , clean_name.split ( '::' ), cppyy.gbl )
        # ===================================================================
        # Identity match holds ONLY for pure C++ classes.
        # For a Python subclass, `resolved_cls` points to the C++ parent class instead.
        # ==================================================================
        return resolved_cls is cls
        # ==================================================================
    except  ( AttributeError , TypeError ) : # =============================
        # ==================================================================
        return False

# ==========================================================================
## Get the type name
#  @code
#  obj = ...
#  print ( 'Object type name is %s' % typename ( obj ) ) 
#  @endcode 
def typename ( o ) :
    """ Get the type name
    >>> obj = ...
    >>> print ( 'Object type name is %s' % typename ( obj ) )
    """
    if callable ( o ) :
        to = type ( o ) 
        if to is __fun_type :
            if '<lambda>' == to.__name__ : return 'lambda'
            return getattr ( to , '__qualname__' , getattr ( to , '__name__' ) )
        
    to = type ( o )
    if is_cpp_class ( to )  :        
        tname = ( getattr ( to , '__cpp_name__' , None ) or                  
                  getattr ( to , '__qualname__' , None ) or 
                  getattr ( to , '__name__'     , None ) )
        if tname : return tname
        
    return getattr ( to , '__qualname__' , None ) or getattr ( to , '__name__' ) 
    
# =============================================================================
## is it a function (or lambda) ?
#  @code
#  obj = ...
#  print ( 'function?' , isfunction ( obj ) ) 
#  @endcode 
def isfunction ( func ) :
    """ Is it a function (or lambda) ?
    """
    return isinstance ( func , ( FunctionType , LambdaType ) )
# =============================================================================
## is it a lambda ?
#  @code
#  obj = ...
#  print ( 'lambda?' , islambda ( obj ) ) 
#  @endcode 
def islambda  ( func ) :
    """ Is it a lambda?
    """
    return isinstance ( func , LambdaType )

# =============================================================================
if __name__ == '__main__' :

    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.utils.core' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
