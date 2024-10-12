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
    'progress_conf', ## default configuration of the C++ progrees bar 
)
# =============================================================================
from ostap.core.core        import Ostap
from ostap.utils.basic      import isatty, terminal_size
from ostap.logger.colorized import allright
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.progress_conf')
else                       : logger = getLogger( __name__ )
# =============================================================================
## configuration of the progress bar
#  @see Ostap::Utils::ProgressConf
#  @see Ostap::Utils::ProgressBar 
def progress_conf ( show = True , timer = True , description = 'Entries:' ) :
    """ Configuration of the progress bar
    - see `Ostap.Utils.ProgressConf`
    - see `Ostap.Utils.ProgressBar` 
    """
    if not show : return Ostap.Utils.ProgressConf ( 0 )
    tty    = isatty () 
    twidth = terminal_size() [1] if tty else 105
    left   = '[ '
    right  = ' ]'
    ld     = len ( description ) 
    return Ostap.Utils.ProgressConf ( 
	max ( 0 , twidth - 18 - ld )          , ## silent if is terminal is too narrow
	allright ( '#'   ) if tty else '#'    , ## 'done' symbol 
        ' '                                   , ## 'not-yet' symbol
        ## allright ( left ) if tty else left    , ## left
        ## allright ( ' [ ' ) if tty else ' [ '  , ## left
        left                               , ## left 
        ## allright ( ' ] '  ) if tty else ' ] ' , ## right
        right                                    , ## right
        description                              , ## what/descritoiont
        timer                                    , ## use the timer
        tty                                      ) ## tty ? 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END
# =============================================================================
