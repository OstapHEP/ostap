#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some simple but useful utilities
#   - suppression of stdout/stderr 
#   - dumpting of stdout/stderr into file 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
#  
# =============================================================================
"""Module with some simple but useful utilities"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    'checked_yes'   , 
    'checked_no'    ,
    'question_mark' ,
    'hand_ok'       ,
    'squared_ok'    ,
    'thumb_up'      ,
    'thumb_down'    ,
) # ===========================================================================
from   ostap.utils.basic import isatty, has_unicode 
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger    import getLogger, logColor, logNoColor 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.logger.stmbols' )
else                       : logger = getLogger( __name__ )
# =============================================================================

show = isatty() and has_unicode()

checked_yes   = '\u2705'     if show else "+"
checked_no    = '\u274c'     if show else "-"
question_mark = '\u2753'     if show else "?"
hand_ok       = '\U0001f44c' if show else 'ok'
squared_ok    = '\U0001f197' if show else 'ok'
thumb_up      = '\U0001f44d' if show else '+'
thumb_down    = '\U0001f44e' if show else '-'

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

        
# =============================================================================
##                                                                     The END 
# =============================================================================
