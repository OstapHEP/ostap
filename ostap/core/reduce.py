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
    'root_factory'              , ## a simple factory for generic deserialization
    'root_store_factory'        , ## a simple factory for generic deserialization
    'is_cpp_class'              , ## pure C++ class? 
    'is_python_subclass_of_cpp' , ## python subclass of  C++ class? 
)
# =============================================================================
from   ostap.utils.core  import is_cpp_class 
import sys, cppyy, copyreg, functools, importlib, pickle   
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
    ## keep arguments with the newly created object  
    obj.__store = params    ## Attention - keep arguments with newly created object!
    return obj 

# ==========================================================================
## Determine whether `cls` is a Python-defined class inheriting from a C++ class.
#  @param  cls Object/type to inspect.
#  @return True if `cls` is a Python subclass of a C++ type, False otherwise.
def is_python_subclass_of_cpp ( cls ) -> bool:
    """ Determine whether `cls` is a Python-defined class inheriting from a C++ class.
    
    :param cls: Object/type to inspect.
    :return: True if `cls` is a Python subclass of a C++ type, False otherwise.
    """
    if not isinstance ( cls , type) : return False

    # A pure C++ class cannot be its own Python subclass
    if is_cpp_class ( cls ) : return False

    # Check if any ancestor class in the MRO (excluding self) is a pure C++ class
    return any ( is_cpp_class ( base ) for base in cls.__mro__[1:] )

# ==============================================================================
if 'win32' == sys.platform : # =================================================
    # ==========================================================================
    SHARED_LIB_EXTS = ( '.dll' ,  )
    # ==========================================================================
    def is_shared_lib ( module_name: str) -> bool :
        """ Check if the module string represents a C++ dynamic library."""
        return module_name.endswith ( SHARED_LIB_EXTS )
else :
    # ==========================================================================
    SHARED_LIB_EXTS = ( ".so", ".dylib")
    def is_shared_lib ( module_name: str) -> bool :
        """ Check if the module string represents a C++ dynamic library."""
        return module_name.startswith ( "lib" ) or module_name.endswith ( SHARED_LIB_EXTS )

# =============================================================================
## de-seriailze cpp-metaclass
#  Restores any class object (C++ or Python) from its name and optional module scope.
#  
#  @param the_name: Fully qualified class name (e.g., 'Ostap::Kinematics::Dalitz', 'MyPyClass')
#  @param module: Optional module or dynamic library name (e.g., 'ROOT', 'cppyy.gbl', 'libOstap')
#  @return Resolved class object.
def cpptype_factory(the_name, module=None):
    """ Restores any class object (C++ or Python) from its name and optional module scope.

    :param the_name: Fully qualified class name (e.g., 'Ostap::Kinematics::Dalitz', 'MyPyClass')
    :param module: Optional module or dynamic library name (e.g., 'ROOT', 'cppyy.gbl', 'libOstap')
    :return: Resolved class object.
    """
    if not the_name or not isinstance(the_name, str):
        raise pickle.UnpicklingError(f"Invalid type name for cpptype_factory: {the_name!r}")

    # =========================================================================
    # Remove leading global C++ namespace specifier (e.g., '::ROOT::Math' -> 'ROOT::Math')
    clean_name = the_name.removeprefix("::")

    # =========================================================================
    ## Tier 1: Try resolving via Shared Library or Python Module first
    # =========================================================================
    if module: # ==============================================================
        # =====================================================================
        if is_shared_lib(module):
            # ROOT.gSystem.Load returns:
            #  0 : Library successfully loaded
            #  1 : Library was already loaded in the current session
            # <0 : Error loading library
            status = ROOT.gSystem.Load(module)
            # If library loading fails, fall through to attempt Cling/cppyy resolution
            if status < 0 : pass
            # =================================================================
        else: # ===============================================================
            # =================================================================
            try: # ============================================================
                # =============================================================
                mod_obj = importlib.import_module(module)
                # Split attributes by '::' or '.' to support nested modules/classes
                delimiters = "::" if "::" in clean_name else "."
                return functools.reduce(getattr, clean_name.split(delimiters), mod_obj)
                # =============================================================
            except (ImportError, AttributeError): # ===========================
                # =============================================================
                pass

    # =========================================================================
    ## Tier 2: Direct C++ Scope Resolution in cppyy.gbl
    # =========================================================================
    try: # ====================================================================
        # =====================================================================
        the_type = getattr ( cppyy.gbl , clean_name , None )
        if the_type is not None: return the_type
        # =====================================================================
    except Exception: # =======================================================
        # =====================================================================        
        pass

    # =========================================================================
    ## Tier 3: Hierarchical C++ Namespace Traversal (e.g., cppyy.gbl.Ostap.Kinematics)
    # =========================================================================
    try: # ====================================================================
        # =====================================================================
        return functools.reduce ( getattr , clean_name.split("::") , cppyy.gbl )
        # =====================================================================
    except (AttributeError, TypeError): # =====================================
        # =====================================================================
        pass

    # =========================================================================
    ## Tier 4: Dynamic Cling JIT Lookup (for templated types like std::vector<...>)
    # =========================================================================
    try: # ====================================================================
        # =====================================================================
        return cppyy.gbl.__getattr__( clean_name )
        # =====================================================================
    except Exception as err: # ================================================
        # =====================================================================
        raise pickle.UnpicklingError ( f"Failed to restore class '{the_name}' (module={module!r}): {err}" ) from err

# =============================================================================
_INTERNAL_NAMES = 'cppyy.gbl' , 'ROOT' , 'ROOT._facade' , '' 
# =============================================================================
## Serializes class types preserving both full class name and module/namespace scope.
#  Correctly handles C++ types, pure Python classes, and Python subclasses of C++ types.
def cpptype_reduce(meta_type, protocol=None):
    """ Serializes class types preserving both full class name and module/namespace scope.
    Correctly handles C++ types, pure Python classes, and Python subclasses of C++ types.
    """
    # =========================================================================
    # 1. Native C++ Class
    # =========================================================================
    if is_cpp_class(meta_type):
        the_name = ( getattr ( meta_type , '__cpp_name__', None ) or
                     getattr ( meta_type , '__qualname__', None ) or
                     meta_type.__name__  )
        return cpptype_factory, ( the_name , None )
    
    # =========================================================================
    # 2. Python Class or Python Subclass of C++ Type (e.g. MyGauss1)
    # =========================================================================
    # Force reading from __dict__ to bypass C++ dispatcher overriding __qualname__
    the_name = ( meta_type.__dict__.get ( '__qualname__' ) or 
                 meta_type.__dict__.get ( '__name__'     ) or 
                 getattr ( meta_type , '__name__', str ( meta_type ) ) )

    # ==========================================================================
    # Safely extract module attributes
    mod_attr = getattr ( meta_type , '__module__'    , None )
    raw_mod  = meta_type.__dict__.get( '__module__'  , mod_attr )

    mod_name = raw_mod

    # ==========================================================================
    # Clean up ROOT internal facades, cppyy pseudo-modules, and empty strings
    if ( not raw_mod or raw_mod in _INTERNAL_NAMES or 'cppyy_internal' in raw_mod ):
        mod_obj  = sys.modules.get ( mod_attr ) if mod_attr else None
        mod_name = getattr ( mod_obj , '__name__', '__main__' )
        if mod_name and 'cppyy_internal' in mod_name : mod_name = '__main__'

    return cpptype_factory, (the_name, mod_name)

# ============================================================================
## Registration for C++ Metaclass Serialization
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

