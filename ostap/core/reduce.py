#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/reduce.py
#  Module with some fundamental utilities for serialization&deserialiation 
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = (
    'root_factory'       , ## a simple factory for generic deserialization
    'root_store_factory' , ## a simple factory for generic deserialization 
)
# =============================================================================
from functools          import reduce as ft_reduce 
import cppyy, copyreg
# =============================================================================
## Trivial factory for deserialization of generic objects
def root_factory ( klass , *params ) :
    """ Trivial factory for (de)serialization of generic objects
    """
    return klass ( *params )
# =============================================================================
## Factory for deserialization of generic objects
#  - dirty trick is here
#  @attention it stores the constructor parameters as local attributes
def root_store_factory ( klass , *params ) :
    """ Factory for deserialization of generic object
    - dirty trick is here
    - attention: it stores the constructor parameters as local attributes
    """
    ## create the objects 
    obj = root_factory ( klass , *params )
    ## keep arguments with the newly created obnject  
    obj.__store = params    ## Attention - keep arguments with newly created object!
    return obj 

# =============================================================================
## de-seriailze cpp-metaclass
#  Restores a C++ class/type from its fully qualified __scopedname__ 
#  using a multi-tier fallback strategy.
def cpptype_factory ( the_name ) :
    """ de-seriailze cpp-metaclass
    Restores a C++ class/type from its fully qualified __scopedname__ 
    using a multi-tier fallback strategy.
    """

    print ( 'RECONSTRUCT' , the_name ) 

    # =========================================================================
    ## Remove leading global namespace specifier (e.g., '::Physics::Particle' -> 'Physics::Particle')
    # =========================================================================
    clean_name = the_name.removeprefix ( "::" )

    # =========================================================================
    ## Tier 1: Direct lookup in cppyy.gbl (fastest if cppyy natively parses '::')
    # =========================================================================    
    the_type = getattr ( cppyy.gbl, clean_name , None )
    if not the_type is None : return the_type 

    # =========================================================================
    ## Tier 2: Traversal through nested namespaces (e.g., cppyy.gbl.Physics.Particle)
    # =========================================================================    
    try: # ====================================================================
        # =====================================================================
        return ft_reduce ( getattr, clean_name.split("::") , cppyy.gbl )
        # =====================================================================
    except AttributeError: # ==================================================
        # =====================================================================
        pass

    # =========================================================================
    ## Tier 3: Heavy-duty JIT evaluation (handles complex templated types like std::vector<...>)
    # =========================================================================    
    try: # ====================================================================
        return eval ( f"cppyy.gbl.{clean_name}" )
         # ====================================================================
    except Exception as err: # ================================================
        # =====================================================================
        raise TypeError ( f"Failed to restore C++ type for '{the_name}': {err}" )

# =============================================================================
## Serialize cpp-metaclass
#  - try __scopedname__
#  - try __cpp_name__
#  - Fallback:  __name__
def cpptype_reduce ( meta_type , protocol = None ) :
    """ Serialize cpp-metaclass
    """
    print ( 'REDUCE' , meta_type )
    
    the_name                   = getattr ( meta_type , '__scopedname__'   , None )
    if not the_name : the_name = getattr ( meta_type , '__cpp_name__'     , None )
    if not the_name : the_name = getattr ( meta_type , '__display_name__' , None )
    if not the_name : the_name = getattr ( meta_type , '__qualname__'     , None )
    if not the_name : the_name = getattr ( meta_type , '__name__'         )
    return cpptype_factory , ( the_name , )

# =============================================================================
## Registration for C++ Metaclass Serialization
# =============================================================================
import ctypes 

# ============================================================================
try : # ======================================================================
    # ========================================================================
    # Anchor on a basic C++ class to traverse the MRO up to the root C++ metaclass
    # Index [-3] points to the base cppyy.CPPScope/CPPMetaclass before (type, object)
    CPP_META = type ( cppyy.gbl.TObject ).__mro__[-3]
    # =======================================================================
except Exception : # ========================================================
    # =======================================================================
    CPP_META = None # =======================================================
    # =======================================================================

# ===========================================================================
if CPP_META : # =============================================================
    # =======================================================================
    # Register via copyreg for standard pickle support - it doesn not WORK! 
    # =======================================================================
    import copyreg 
    copyreg.pickle ( CPP_META , cpptype_reduce )

    # =======================================================================
    # Register in dispatch_table to intercept class-level pickle.dumps(CPP_CLASS)
    # ======================================================================
    copyreg.dispatch_table [ CPP_META ] = cpptype_reduce

# =============================================================================
if '__main__' == __name__ :

    # =========================================================================
    from   ostap.logger.logger import getLogger 
    logger = getLogger ( 'ostap.core.reduce' )
  
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if not CPP_META : logger.warning ( "CPP_META is invalid!" ) 
    
# =============================================================================
##                                                                      The END 
# =============================================================================

