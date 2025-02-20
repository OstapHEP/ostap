#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostath/integrator.py
#  Decoration modlel for C++ class Ostap::Math::Integrator
#  @see Ostap::Math::Integrator
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2022-12-18
# =============================================================================
""" Decoration modlel for C++ class Ostap::Math::Integrator
- see `Ostap.Math.Integrator`
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2022-12-18"
__all__     = (
    "Integrator" , ## C++ numerical integrator  
    )
# =============================================================================
from ostap.math.base      import Ostap, doubles 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.integrator' )
else                       : logger = getLogger ( __name__                )

# =============================================================================
## actual C++ interator 
Integrator = Ostap.Math.Integrator
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                     The EMD   
# =============================================================================

