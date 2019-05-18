#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some multiprocessing functionality for Ostap 
#  Currently it is not loaded on default, and requires manual activation
#  @see GaudiMP.Parallel
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
# =============================================================================
""" Multiprocessing functionality for Ostap
Currently it is not loaded on default, and requires manual activation
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    ) 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.kisa' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
logger.debug ( 'Multiprocessing functionality for Ostap')
# =============================================================================
import ostap.parallel.parallel as Parallel
WorkManager = Parallel.WorkManager
# =============================================================================


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
