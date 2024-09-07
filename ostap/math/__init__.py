#!/usr/env python
# -*- coding: utf-8 -*-
# =============================================================================
import sys, os
if sys.warnoptions or os.environ.get ( 'OSTAP_CMAKE_TEST', False ) :
    import warnings 
    with warnings.catch_warnings():
        warnings.simplefilter ( "always" )
        import cppyy
# =============================================================================
import ostap.math.reduce 
# =============================================================================
##                                                                      The END 
# =============================================================================
