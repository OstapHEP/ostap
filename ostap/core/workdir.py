#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-03-10
# =============================================================================
"""
""" 
# =============================================================================
__author__  = "Vanya BELYAEV  Ivan.Belyaev@itep.ru"
__date__    = "2014-03-10"
__version__ = "$Revision$"
__all__     = (
    'workdir' , ## workdir where selector cache, etc. can be placed  
    )
# =============================================================================
import os 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.core.workdir'  )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## import ostap.core.config as _config 
workdir = os.environ.get('OSTAP_DIR') or os.environ.get('OSTAPDIR') or '$HOME/.ostap'

workdir = os.path.expandvars ( workdir )
workdir = os.path.expanduser ( workdir )
workdir = os.path.expandvars ( workdir )
workdir = os.path.expanduser ( workdir )

if not os.path.exists ( workdir ):
    try : 
        os.mkdir    ( workdir )
        wdir    = os.path.join(workdir, "cache")
        os.mkdir ( wdir    )
        logger.debug ( 'Create working cache directory: %s' % wdir )
    except :
        pass 

# =============================================================================
if '__main__' == __name__ : 
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
