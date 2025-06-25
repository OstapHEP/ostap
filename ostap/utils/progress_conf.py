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
from ostap.core.ostap_types import integer_types 
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
    PC = Ostap.Utils.ProgressConf 
    if   isinstance ( show , PC )  : return show 
    elif not show                  : return PC ( 0 )
    
    tty    = isatty ()
    if isinstance ( show , integer_types ) and 40 <= show <= 512 : twidth = show 
    else : twidth = terminal_size() [ 0 ] if tty else 105
    
    done   = '#'
    notyet = ' '
    assert len ( done ) == len ( notyet ) , "Mismatch in symbol lengths"
    left   = '[ '
    right  = ' ]'
    ld     = len ( description ) 
    ll     = len ( left        )
    lr     = len ( right       )
    twidth = max ( 0 , twidth - 20 - ll - lr - ld ) // len ( done )
    if tty :
        done  = allright ( done  )
        left  = allright ( left  )
        right = allright ( right )
        
    return Ostap.Utils.ProgressConf ( 
        twidth        , ## silent if terminal is too narrow
	    done          , ## 'done' symbol 
        notyet        , ## 'not-yet' symbol
        left          , ## left 
        right         , ## right
        description   , ## what/description
        timer         , ## use the timer
        tty           ) ## tty ? 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END
# =============================================================================
