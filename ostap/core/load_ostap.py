#!/usr/bin/env ipython 
# -*- coding: utf-8 -*-
# =============================================================================
## @file load_ostap.py
#  @date   2012-02-15
#  @author Vanya BELYAEV Ivan.Belyaevitep.ru
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
__logger = getLogger ( 'ostap.core.load_ostap' )
# =============================================================================
import ROOT, sys
# =============================================================================
## 1) load Ostap-style file
# =============================================================================
import ostap.plotting.style 
# =============================================================================
# The Heart 
# =============================================================================
## load zillions of decorations for ROOT-objects
import ostap.core.pyrouts        ## NB: the most important line!
import ostap.io.zipshelve as DBASE

# =============================================================================
## minor decoration for default shelve module 
import ostap.io.shelve_ext

# =============================================================================
## import useful context managers
from ostap.logger.utils       import *
from ostap.utils.utils        import *
from ostap.utils.progress_bar import progress_bar
# ============================================================================= 
## prepend the path 
if '.' not in sys.path :
    __logger.debug('Prepend sys.path with $PWD')
    sys.path = ['.'] + sys.path 
    
# =============================================================================
from ostap.core.core import cpp, Ostap, VE, SE, WSE, hID
from ostap.math.base import doubles
# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = __logger )

# =============================================================================
##                                                                      The END 
# =============================================================================

