#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/ostap_types.py
#  core types for ostap 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Core types for ostap 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    ##
    'num_types'       , ## numerical types (real)
    'integer_types'   , ## integer types
    'string_types'    , ## string  types
    'str_types'       , ## string  types (only string!)
    'listlike_types'  , ## list-like types (list, tuple, set, collections.Iterable)
    'list_types'      , ## list types (list, tuple)
    'dict_types'      , ## dict types 
    'dictlike_types'  , ## dict-like types 
    'long_type'       , ## long-type
    'sequence_types'  , ## sequence types
    'iterable_types'  , ## iterable 
    'sized_types'     , ## sized types
    'path_types'      , ## path-like types
    'ordered_dict'    , ## normal/ordered dict 
    ##
    'is_integer'      , ## is a value of int-like type?
    'is_number'       , ## is a value of numeric  type?
    'is_good_number'  , ## is a value of numeric  type and not-NaN,no-Inf ?
    ##
    'is_string'       , ## is a value of str-type? 
    'is_string_like'  , ## is a value of string-like type? 
    'is_list'         , ## is a value of list/tuple type?
    'is_list_like'    , ## is a value of list-like type?
    ##
    'all_integers'    , ## all argumets of integer types?
    'all_numerics'    , ## all argumets of numeric types?
    'all_strings'     , ## all argumets of string  types?
    )
# =============================================================================
from   collections.abc import Collection, Sequence, Iterable, Mapping, Sized, Generator   
import array, sys, os. math  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.core.ostap_types' )
else                       : logger = getLogger( __name__     )
# =============================================================================
logger.debug ( 'Core objects/classes/types for Ostap')
# =============================================================================
long           = int
string_types   = bytes , str 
integer_types  = int   ,
long_type      = int
# =============================================================================
iterable_types  = Iterable ,
num_types       = integer_types + ( float , ) 
str_types       = str ,
list_types      = list , tuple
listlike_types  = list_types + ( set , Sequence , array.array )
# =============================================================================
if sys.warnoptions or os.environ.get ( 'OSTAP_CMAKE_TEST', False ) :
    import warnings 
    with warnings.catch_warnings():
        warnings.simplefilter ( "ignore" , category = DeprecationWarning )        
        warnings.simplefilter ( "ignore" , category = UserWarning        )        
        import cppyy
else :    
    import cppyy
# ==============================================================================
std = cppyy.gbl.std 
if hasattr ( std , 'string'      ) : string_types += ( std.string      , )
if hasattr ( std , 'string_view' ) : string_types += ( std.string_view , )
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import numpy as np
    listlike_types  = listlike_types + ( np.ndarray , )
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    pass 
# =============================================================================
dict_types      = dict ,
dictlike_types  = dict ,  Mapping  
sequence_types  = listlike_types + ( Sequence , Collection , Iterable , Generator )
sized_types     = Sized ,
path_types      = string_types
# =============================================================================
path_types = string_types + ( os.PathLike , )
# =============================================================================
## sometimes we need to ensure that dictionary is ordered 
ordered_dict = dict

# =============================================================================
## Is this number of a proper integer?
def is_integer ( v ) :
    """ Is this number of a proper integer?"""
    return isinstance ( v , integer_types ) 

# =============================================================================
## Is this number of a proper numeric type
def is_number  ( v ) :
    """ Is this number of a proper numeric?"""
    return isinstance ( v , num_types   ) 

# =============================================================================
## good numeric value 
def is_good_number  ( v ) :
    """ Is numeric type and good value?"""
    return isinstance ( v , num_types   ) and math.isfinite ( v ) 

# =============================================================================
## is  value of str-type?
def is_string ( v ) :
    """ Is  value of str-type?"""
    return isinstance ( v , str_types )

# =============================================================================
## is  value of string-like -type?
def is_string_like ( v ) :
    """Is  value of string-like type?"""
    return isinstance ( v , string_types )

# =============================================================================
## is list type?
def is_list  ( v ) :
    """ Is value of list type (list ot tuple)"""
    return isinstance ( v , list_type )

# =============================================================================
## is list-like type?
def is_list_like  ( v ) :
    """ Is value of list-like type (list but not a string-like!)"""
    return isinstance ( v , listlike_type ) and not isinstance ( v , string_types )

# =============================================================================
## are all values of integer-line types?
def all_integers ( *args ) :
    """ Are all value of integer-line types?
    """
    return all ( isinstance ( a , integer_types ) for a in args )

# =============================================================================
## are all values of numeric types?
def all_numerics ( *args ) :
    """ Are all values of numeric types?
    """
    return all ( isinstance ( a , num_types ) for a in args )

# =============================================================================
## are all values of string types?
def all_strings ( *args ) :
    """ Are all values of string types?
    """
    return all ( isinstance ( a , string_types ) for a in args )
    
# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    
# =============================================================================
##                                                                      The END 
# =============================================================================
