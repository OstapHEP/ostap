#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file pickling.py
# Helper module to define pickling for various databases 
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2010-04-30
# =============================================================================
""" Helper module to define pickling for various databases
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2010-04-30"
__version__ = "$Revision$" 
# =============================================================================
__all__ = (
    'DEFAULT_PROTOCOL'    ,
    'HIGHEST_PROTOCOL'    ,
    'PROTOCOL'            ,
    'Pickler'             , 
    'Unpickler'           ,
    'BytesIO'             ,
    'dumps'               ,
    'loads'               ,
    ## 
    'pickle_dependencies' , ## inspect pickeld stream 
)
# =============================================================================
import pickle 
from   io                import BytesIO 
import ostap.core.config as     config
from   ostap.core.reduce import cpptype_reduce, CPP_META 
import sys, array, pickletools, re, cppyy 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.pickling' )
else                      : logger = getLogger ( __name__               )
# =============================================================================
DEFAULT_PROTOCOL = pickle.DEFAULT_PROTOCOL
HIGHEST_PROTOCOL = pickle.HIGHEST_PROTOCOL

# =============================================================================
PROTOCOL         = config.protocol
if   HIGHEST_PROTOCOL < PROTOCOL : PROTOCOL = HIGHEST_PROTOCOL
elif 0                > PROTOCOL : PROTOCOL = DEFAULT_PROTOCOL
# =============================================================================

# =============================================================================
## Regular pickle stuff
# =============================================================================

dump             = pickle.dump
load             = pickle.load
dumps            = pickle.dumps
loads            = pickle.loads
Unpickler        = pickle.Unpickler
Pickler          = pickle.Pickler

    
# =============================================================================
if CPP_META : # ===============================================================
    # =========================================================================
    ## @class CppAwarePickler
    #  Custom Pickler that intercepts C++ classes/metaclasses at the C level
    #  before CPython's pickle.c attempts to serialize them via standard save_type.
    class CppAwarePickler ( pickle.Pickler ) :
        """ Custom Pickler that intercepts C++ classes/metaclasses at the C level
        before CPython's _pickle.c attempts to serialize them via standard save_type.
        """
        ## Intercept C++ classes/metaclasses via fast C-level type check
        def reducer_override(self, obj ):
            """ Intercept C++ classes/metaclasses via fast C-level type check"""
            if CPP_META is not None and ( isinstance ( obj , CPP_META ) or obj is CPP_META ):
                return cpptype_reduce ( obj ) 
            
            # Fallback to default CPython serialization behavior for all other objects
            return NotImplemented

    # =========================================================================
    ## Transparent wrapper for pickle.dump supporting C++ metaclasses.
    def _cpp_dump ( obj             , file ,
                    protocol        = None , * ,
                    fix_imports     = True ,
                    buffer_callback = None ) :
        """Transparent wrapper for pickle.dump supporting C++ metaclasses."""
        pickler = CppAwarePickler ( file                              , 
                                    protocol        = protocol        , 
                                    fix_imports     = fix_imports     , 
                                   buffer_callback = buffer_callback )
        pickler.dump ( obj )
        
        
    # =========================================================================
    ## Transparent wrapper for pickle.dumps supporting C++ metaclasses
    def _cpp_dumps ( obj ,
                     protocol        = None , *,
                     fix_imports     = True ,
                     buffer_callback = None ) :
        """ Transparent wrapper for pickle.dumps supporting C++ metaclasses."""
        buf     = BytesIO ()
        _cpp_dump ( obj             , buf             ,
                    protocol        = protocol        ,
                    fix_imports     = fix_imports     ,
                    buffer_callback = buffer_callback )
        return buf.getvalue()
    
        
    
    Pickler = CppAwarePickler
    dump    = _cpp_dump
    dumps   = _cpp_dumps
    
else :

    logger.warning ( "CPP_META is invalid!" ) 

# =============================================================================
## @class XVersionUnpickler
#  (Deserialization - ALWAYS AVAILABLE)
class XVersionUnpickler(pickle.Unpickler):
    """ Cross-version Unpickler mapping all legacy and modern PyROOT internal 
    deserialization modules and C++ types. Always accessible.
    """
    def find_class(self, module, name):
        if ( module.startswith ( 'ROOT.libROOTPythonizations') or 
             module.startswith ( 'libROOTPythonizations') or 
             module.startswith ( 'libROOTTPython') or 
             module in ('ROOT.pyroot', 'PyROOT', 'cppyy.gbl')):
            
            # A. Try modern ROOT internal module first
            try:
                import ROOT.libROOTPythonizations as root_pyons
                if hasattr(root_pyons, name):
                    return getattr(root_pyons, name)
            except (ImportError, AttributeError):
                pass
            
            # B. Fallback to ROOT module
            try:
                import ROOT
                if hasattr(ROOT, name):
                    return getattr(ROOT, name)
            except ImportError:
                pass

            # C. Fallback to cppyy
            try:
                if hasattr(cppyy.gbl, name):
                    return getattr(cppyy.gbl, name)
            except Exception:
                pass

        return super().find_class(module, name)


# =============================================================================
## Transparent wrapper for pickle.load with ROOT cross-version support
def _xv_load ( file , *,
                 fix_imports = True     ,
                 encoding    = "ASCII"  , 
                 errors      = "strict" ,
                 buffers     = None     ) :
    """ Transparent wrapper for pickle.load with ROOT cross-version support.
    """
    unpickler = XVersionUnpickler ( file ,
                                    fix_imports = fix_imports ,
                                    encoding    = encoding    ,
                                    errors      = errors      ,
                                    buffers     = buffers     )
    return unpickler.load()


# =============================================================================
## Transparent wrapper for pickle.loads with ROOT cross-version support."""
def _xv_loads ( bytes_object , * ,
                fix_imports  = True     ,
                encoding     = "ASCII"  ,
                errors       = "strict" ,
                buffers      = None     ) :
    """ Transparent wrapper for pickle.loads with ROOT cross-version support.
    """
    buf = BytesIO(bytes_object)
    return _xv_load ( buf ,
                      fix_imports = fix_imports ,
                      encoding    = encoding    ,
                      errors      = errors      ,
                      buffers     = buffers     )

Unpickler = XVersionUnpickler
load      = _xv_load
loads     = _xv_loads

# =============================================================================
## py-moduel patters for pickeltools 
_PY_MODULE_PATTERN        = re.compile(r'^[a-zA-Z_][a-zA-Z0-9_]*(\.[a-zA-Z_][a-zA-Z0-9_]*)+$')
_CPP_GLOBAL_CLASS_PATTERN = re.compile(r'^[A-Z][a-zA-Z0-9_]*$')
# Exclude common non-C++ uppercase strings found in pickles (codecs, constants, etc.)
_EXCLUDED_STRINGS = {'UTF-8', 'ASCII', 'ASCII_GENERAL_CI', 'RAW', 'NONE', 'TRUE', 'FALSE'}

# =============================================================================
## Pickle Dependency Inspector
#  Inspects raw pickle bytes without executing/loading them to extract 
#   required Python modules, functions, and C++ scopes/types.
# 
#  @param pickled_bytes : Raw bytes or bytearray of a pickled object
#  @returm return  : A dictionary containing lists of required Python modules,
#                  C++ types, and all extracted string symbols.
# @code
# pickled = ... ## pickled representation of something
# result  = pickle_dependencies ( pickled ) 
# @encode
# =============================================================================
def pickle_dependencies ( pickled_bytes ):
    """
    Inspects raw pickle bytes without executing/loading them to extract 
    required Python modules, functions, and C++ scopes/types.
    
    - param pickled_bytes: Raw bytes or bytearray of a pickled object.
    - return  A dictionary containing lists of required Python modules,
             C++ types, and all extracted string symbols.
    """
    found_strings  = set ()
    global_imports = set ()

    # =========================================================================
    ## Iterate through all pickle opcodes safely
    for opcode, arg, pos in pickletools.genops ( pickled_bytes ):
        # =====================================================================
        ## 1. Direct GLOBAL opcode imports (module_name class_or_func_name)
        if opcode.name == 'GLOBAL' and arg:
            global_imports.add(arg.split()[0])

        # =====================================================================
        ## 2. Extract string constants (module paths, C++ names, or arguments)
        elif opcode.name in ('SHORT_BINUNICODE', 'BINUNICODE', 'UNICODE') and arg:
            found_strings.add(str(arg))

    # =========================================================================
    ## Extract Python modules using pre-compiled regex (e.g., 'ostap.core.reduce')
    # =========================================================================
    py_modules = {s for s in found_strings if _PY_MODULE_PATTERN.match(s)}
    py_modules.update(global_imports)
    
    # =========================================================================
    # Extract C++ scopes/types:
    # - Namespaced C++ types (e.g., 'Ostap::Kinematics::Dalitz')
    # - Global C++ classes starting with uppercase (e.g., 'TH1F', 'TTree', 'MyGlobalClass')
    cpp_types = {
        s for s in found_strings 
        if '::' in s or (_CPP_GLOBAL_CLASS_PATTERN.match(s) and s.upper() not in _EXCLUDED_STRINGS)
    }
    
    return {
        "modules"   : tuple ( sorted ( list ( py_modules    ) ) ) ,
        "cpp_types" : tuple ( sorted ( list ( cpp_types     ) ) ) ,
        "symbols"   : tuple ( sorted ( list ( found_strings ) ) ) }


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    logger.info ( 'Pickling protocol: %s' % PROTOCOL )

    if not CPP_META : logger.warning ( "CPP_META is invalid!" ) 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
