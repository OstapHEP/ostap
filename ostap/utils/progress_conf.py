#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file progress_conf.py
#  Default configuration for C++ progress bar
#  @see Ostap::Utils::ProgressBar
#  @see Ostap::Utils::ProgressConf
#  @date 2022-09-06
#  @author Vanya BELYAEV Ivan/Belyaev@itep.ru
# =============================================================================
""" Default configuration for C++ progress bar
- see Ostap.Utils.ProgressBar
- see Ostap.Utils.ProgressConf
"""
# =============================================================================
__version__ = "$Revision:$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2022-09-06"
__all__     = (
    'progress_conf', ## defautl configuration of the C++ progrees bar 
    )
# =============================================================================
from   ostap.core.core   import Ostap
from   ostap.utils.basic import isatty, terminal_size
from   ostap.logger.colorized   import allright
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.progress_conf')
else                       : logger = getLogger( __name__ )
# =============================================================================
## configuration of the progress bar 
twidth        = terminal_size () [ 1 ] if isatty () else 110 
progress_conf = Ostap.Utils.ProgressConf (
    twidth - 15 if 25 < twidth else 0 , ## silent if is termnal is too narrow
    allright ( '#'   )                , ## 'done' symbol 
    ' '                               , ## 'not-yet' symbol 
    allright ( ' [ ' )                , ## left 
    allright ( '] '  )                , ## right 
    True                              ) ## use the timer 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END
# =============================================================================
