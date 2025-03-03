#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/linalg.py
#  Few utilities to simplify linear algebra manipulations 
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
"""Few utilities to simplify linear algebra manipulations 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2009-09-12"
__version__ = ""
# =============================================================================
__all__     = (
    'mgetter'  , ## get  (i,j) element from matrix-like object
    'checkops' , ## check the allowed operations
    'LinAlgT'  , ## LinAlgenra type&decorator store 
    )
# =============================================================================
from   ostap.math.linalg2   import mgetter, checkops  
from   ostap.math.linalgt   import LinAlgT         
import ostap.math.linalgg
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.linalg2' )
else                       : logger = getLogger ( __name__             )
# =============================================================================

# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
