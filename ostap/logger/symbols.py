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
""" Module with some simple but useful utilities """
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    'show'                , ## 
    'checked_yes'         ,   
    'checked_no'          ,
    'question_mark'       ,
    'hand_ok'             ,
    'squared_ok'          ,
    'thumb_up'            ,
    'thumb_down'          ,
    'clock'               ,
    'ram'                 ,
    'runner'              ,
    'finish'              ,
    ## 
    'clock_ticks'         ,
    'arrow_left'          ,   
    'arrow_right'         ,  
    'arrow_rightleft'     ,
    ##
    'times'               , 
    'ditto'               , 
    ##
    'plus_minus'          , 
    'minus_plus'          ,
    #3
    'less_or_equal'       ,
    'greater_or_equal'    , 
    'much_less'           , 
    'much_greater'        , 
    'equivalent'          ,
    'similar'             ,
    'approximate'         ,
    'not_equal'           ,
    ##
    'langle'              , 
    'rangle'              , 
    'ellipsis'            ,
    ##
    'union'               , 
    'intersection'        ,
    ##
    'tree'                ,    
    'chain'               ,    
    'branch'              , 
    'leaves'              , 
    'cabinet'             , 
    'frame'               ,       
    'histogram'           ,  
    'graph'               ,      
    'palette'             ,    
    'document'            ,   
    'tape'                ,
    'tape_cartridge'      ,
    'folder'              ,
    'light_bulb'          , 
    'weight_lifter'       ,
    ##
    'delta_symbol'        , 
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

langle           = '\U00003008' if show else '<'
rangle           = '\U00003009' if show else '>'
ellipsis         = '\U00002026' if show else '...'

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

tree             = '\U0001f334' if show else ''
chain            = '\U000026d3' if show else '' 
branch           = '\U00002E19' if show else '' 
leaves           = '\U0001F343' if show else '' 
cabinet          = '\U0001F5c4' if show else '' 
frame            = '\U0001F5BC' if show else '' 
histogram        = '\U0001F4CA' if show else ''
graph            = '\U0001F4C8' if show else '' 
palette          = '\U0001f3A8' if show else '' 
document         = '\U0001F5CE' if show else '' 
tape             = '\U00002707' if show else '' 
tape_cartridge   = '\U0001F5AD' if show else '' 
folder           = '\U0001F4C2' if show else '' 
light_bulb       = '\U0001F4A1' if show else '' 

less_or_equal    = '\U00002264' if show else '<='
greater_or_equal = '\U00002265' if show else '=>'
much_less        = '\U0000226A' if show else '<<'
much_greater     = '\U0000226B' if show else '>>'
equivalent       = '\U00002261' if show else '='
similar          = '\U0000223C' if show else '~'
approximate      = '\U00002248' if show else '~='
not_equal        = '\U00002260' if show else '!='


weight_lifter    = '\U0001F3CB' if show else '' 

union            = '\U000022C3' if show else ''
intersection     = '\U000022C2' if show else ''
exclusive_or     = '\U000022BB' if show else '^'
difference       = '\U000022BB' if show else '-'

## capital Greek Sigma 
sum_symbol       = '\U00002211'           if show else 'sum '
## lowercase Greek sigma 
rms_symbol       = '\U000003C3'           if show else 'rms '
## squared lower case Greek sigma 
dispersion_sym   = '\U000003c3\U000000B2' if show else 'D '
## squared lower case Greek sigma 
variance_sym     = '\U000003c3\U000000B2' if show else 'var '
## Delta symbol 
delta_symbol     = '\U00000394'           if show else 'delta'


# ==================================================
def the_sum  ( what ) : return '%s%s'   % ( sum_symbol , what ) 
def the_mean ( what ) : return '%s%s%s' % ( langle , what , rangle ) 
def the_rms  ( what ) : return '%s%d'   % ( rms_symbol , what ) 

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

        
# =============================================================================
##                                                                     The END 
# =============================================================================
