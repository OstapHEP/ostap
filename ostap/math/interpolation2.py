#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/interpolation2.py
#  Module with some useful scipy-based utilities for dealing with interpolation.
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2020-02-28
# =============================================================================
"""Useful scipy-based utilities for dealing with interpolation.
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2018-07-22"
__all__     = (
    )
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.interpolation2' )
else                       : logger = getLogger ( __name__                    )
# =============================================================================
try : # =======================================================================
    from ostap.math.sp_interpolation import SplineInterpolator
    __all__ =  'SplineInterpolator',
    # =========================================================================
except  ImportError : # =======================================================
    # =========================================================================
    pass 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
