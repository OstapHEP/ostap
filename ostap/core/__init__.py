#!/usr/env python
# -*- coding: utf-8 -*-
# =============================================================================
import sys, os
if sys.warnoptions or os.environ.get ( 'OSTAP_CMAKE_TEST', False ) :
    import warnings 
    with warnings.catch_warnings() :
        warnings.simplefilter ( "ignore" , category = DeprecationWarning )        
        warnings.simplefilter ( "ignore" , category = UserWarning        )        
        import cppyy
        
import ostap.core.meta_info 
# =============================================================================
##                                                                      The END 
# =============================================================================
