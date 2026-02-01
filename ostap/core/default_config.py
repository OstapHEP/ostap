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
""" Default configuration of ostap
"""
# =============================================================================
__version__  = "$Revision$"
__author__   = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__     =  "2019-05-19"
__all__      = (
    'verbose'       , ## verbose processing?
    'config_files'  , ## configuration files to read 
    'startup_files' , ## startup files 
    )
# =============================================================================

# =============================================================================
arg_parse    = True    ## Parse command-line arguments? 
batch        = False

silent       = False   ## silent pprocesisng ?
quiet        = False   ## quiet   processing ?
debug        = False   ## debug   processing ?
verbose      = False   ## verbose processing ?
level        = -1      ## print level 
color        = True    ## use colors ? 
show_unicode = False   ## show unicode in logfiles?
#
dump_config  = '.ostap_config.dump'

build_dir    = ''                    ## ROOT/Ostap build directory 
cache_dir    = '$HOME/.cache/ostap'  ## Cache directory 
tmp_dir      = ''                    ## TMP directory 

webdisplay   = 'off'                 ## Use web-display? 

parallel     = 'PATHOS'              ## parallel engine 
ncpus        = -1                    ## use all CPUs 
implicitMT   = True                  ## implicit multithreading 
profile      = False                 ## profile the execution? 

table_style  = 'default'             ## Table style


protocol     = ''                    ## pickling protocol 

# =============================================================================
## configuration files to read
# =============================================================================
config_files = [
    u'$OSTAPDIR/.ostaprc'           , ## .ostaprc from central directory 
    u'$HOME/.ostaprc'               , ## .ostaprc from home directory 
    u'~/.ostaprc'                   , ## .ostaprc from home directory 
    u'$HOME/.config/ostap/.ostaprc' , ## .ostaprc from config directory 
    u'~/.config/ostap/.ostaprc'     , ## .ostaprc from config directory 
    u'.ostaprc'                       ## .ostaprc from local directory 
    ]

# =============================================================================
## startup/logon files to be executed:
startup_files = [ '$HOME/.ostap.py' ,
                  '~/.ostap.py'     ,        
                  './.ostap.py'     ]

# =============================================================================
## C++ macros to load
macros        = []

# =============================================================================
## Python commands to be executed 
commands      = []
# =============================================================================
if '__main__' == __name__ :
    
    from   ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.core.default_config' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
