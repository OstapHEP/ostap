#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/default_config.py
#  Default configuration of ostap
#  Ostap parses the following configuration files :
#  - <code>$OSTAPDIR/.ostaprc</code>
#  - <code>$HOME/.ostaprc</code>
#  - <code>~/.ostaprc</code>
#  - <code>- ~/.config/ostap/.ostaprc<.code>
#  - <code>- $HOME/.config/ostap/.ostaprc<.code>
#  - <code>.ostaprc</code>
#  - <code>$OSTAP_CONFIG</code>
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-05-19
# =============================================================================
"""Default configuration of ostap
"""
# =============================================================================
__version__  = "$Revision$"
__author__   = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__     =  "2019-05-19"
__all__      = (
    'quiet'        , ## quiet   processing?
    'verbose'      , ## verbose processing?
    'config_files' , ## configuration files to read 
    )
# =============================================================================
quiet        =  False
verbose      =  False
web          = 'off'

## configuration files to read 
config_files = (
    u'$OSTAPDIR/.ostaprc'           , ## .ostaprc from central directory 
    u'$HOME/.ostaprc'               , ## .ostaprc from home directory 
    u'~/.ostaprc'                   , ## .ostaprc from home directory 
    u'$HOME/.config/ostap/.ostaprc' , ## .ostaprc from config directory 
    u'~/.config/ostap/.ostaprc'     , ## .ostaprc from config directory 
    u'.ostaprc'                       ## .ostaprc from local directory 
    )

# =============================================================================
if '__main__' == __name__ :
    
    from   ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.core.default_config' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
