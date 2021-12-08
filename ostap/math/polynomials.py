#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/polynomials.py
#  Module with some useful utilities for dealing with polynonials
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = ()
# =============================================================================
import  ROOT
import  ostap.math.base 
import  ostap.math.reduce   
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.polynomials' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================


# =============================================================================
## decorated classes 
_decorated_classes_  = ()

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================

