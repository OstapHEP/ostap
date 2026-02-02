#!/usr/bin/env ipython 
# -*- coding: utf-8 -*-
# =============================================================================
## @file load_ostap.py
#  @date   2012-02-15
#  @author Vanya BELYAEV Ivan.Belyaevitep.ru
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
    
# =============================================================================
from ostap.core.core     import cpp, Ostap, VE, SE, WSE, hID, fID, dsID, funID 
from ostap.math.base     import doubles
from ostap.io.root_files import ROOTCWD 
# =============================================================================
if '__main__' == __name__ :

    # =============================================================================
    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.core.load_ostap' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = __logger )

# =============================================================================
##                                                                      The END 
# =============================================================================

