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
    'checked_yes'      ,   
    'checked_no'       ,
    'question_mark'    ,
    'hand_ok'          ,
    'squared_ok'       ,
    'thumb_up'         ,
    'thumb_down'       ,
    'clock'            ,
    'ram'              ,
    'runner'           ,
    'finish'           , 
    'clock_ticks'      ,
    'arrow_left'       ,   
    'arrow_right'      ,  
    'arrow_rightleft'  , 
    'times'            , 
    'ditto'            , 
    'plus_minus'       , 
    'minus_plus'       ,
    'less_or_equal'    ,
    'greater_or_equal' , 
    'much_less'        , 
    'much_greater'     , 
    'equivalent'       ,
)
# ===========================================================================
from   ostap.utils.basic import isatty, has_unicode 
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger    import getLogger, logColor, logNoColor 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.logger.stmbols' )
else                       : logger = getLogger( __name__ )
# =============================================================================

show = isatty() and has_unicode()

checked_yes      = '\u2705'      if show else "+"
checked_no       = '\u274c'      if show else "-"
question_mark    = '\u2753'      if show else "?"
hand_ok          = '\U0001f44c'  if show else 'ok'
squared_ok       = '\U0001f197'  if show else 'ok'
thumb_up         = '\U0001f44d'  if show else '+'
thumb_down       = '\U0001f44e'  if show else '-'
clock            = '\U0001f550'  if show else '' 
ram              = '\U0001f40f'  if show else ''
runner           =  '\U0001f3c3' if show else ''
finish           =  '\U0001f3c1' if show else ''

arrow_left       = '\U00002190' if show else '<-'
arrow_right      = '\U00002192' if show else '->'
arrow_rightleft  = '\U00002194' if show else '<->'

arrows_all       = ''.join ( ( '\U00002190' ,'\U00002196' ,
                               '\U00002191' ,'\U00002197' ,
                               '\U00002192' ,'\U00002198' ,
                               '\U00002193' ,'\U00002199' ) ) \
                               if show else ( '<-' , '\\' , '|' , '/' , '->' , '\\' , '|' , '/' )

clock_ticks      = ''.join ( ( '\U0001f558' , '\U0001f567' ,
                               '\U0001f550' , '\U0001f55C' , 
                               '\U0001f551' , '\U0001f55D' , 
                               '\U0001f552' , '\U0001f55E' , 
                               '\U0001f553' , '\U0001f55F' , 
                               '\U0001f554' , '\U0001f560' , 
                               '\U0001f555' , '\U0001f561' , 
                               '\U0001f556' , '\U0001f562' , 
                               '\U0001f557' , '\U0001f563' , 
                               '\U0001f558' , '\U0001f564' , 
                               '\U0001f559' , '\U0001f565' , 
                               '\U0001f55A' , '\U0001f566' ) ) \
                               if show else '|/-\\'

times            = '\U00002a2f' if show else 'x'
plus_minus       = '\U000000B1' if show else '+/-'
minus_plus       = '\U00002213' if show else '-/+'
ditto            = '\U00003003' if show else '//'

less_or_equal    = '\U00002266' if show else '<='
greater_or_equal = '\U00002267' if show else '=>'
much_less        = '\U0000226A' if show else '<<'
much_greater     = '\U0000226B' if show else '>>'
equivalent       = '\U00002261' if show else '='

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

        
# =============================================================================
##                                                                     The END 
# =============================================================================
