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
    'scissors'            ,
    'oil_drum'            , ## old drum
    'brain'               ,
    'kitchen_knife'       ,
    'axe'                 , 
    ##
    'delta_symbol'        ,
    'number'              , 
    ##
    'union'               , 
    'intersection'        , 
    'exclusive_or'        ,  
    'difference'          ,
    ## 
    'iteration'           , 
    ##    
    'labels'               
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
scissors         = '\U00002704' if show else '' ## ...02 ? 
oil_drum         = '\U0001F6E2' if show else '' 
brain            = '\U0001F9E0' if show else ''
kitchen_knife    = '\U00002796' if show else ''
axe              = '\U0001FA93' if show else ''

union            = '\U000022C3' if show else ''
intersection     = '\U000022C2' if show else ''
exclusive_or     = '\U000022BB' if show else '^'
difference       = '\U000022BB' if show else '-'

iteration        = '\U00003005' if show else ''

## indices: circled numbes from 0 to 50 (inclusive) 
indices2 = '\U000024FF' + \
    '\U0000278A\U0000278B\U0000278C\U0000278D\U0000278E\U0000278F\U00002790\U00002791\U00002792\U00002793' + \
    '\U000024EB\U000024EC\U000024ED\U000024EE\U000024EF\U000024F0\U000024F1\U000024F2\U000024F3\U000024F4' + \
    '\U00003251\U00003252\U00003253\U00003254\U00003255\U00003256\U00003257\U00003258\U00003259\U0000325A' + \
    '\U0000325B\U0000325C\U0000325D\U0000325E\U0000325F\U000032B1\U000032B2\U000032B3\U000032B4\U000032B5' + \
    '\U000032B6\U000032B7\U000032B8\U000032B9\U000032BA\U000032BB\U000032BC\U000032BD\U000032BE\U000032BF' if show else tuple ( '%s' % i for i in range ( 51 ) ) 

## indices = '\U000024EA' + \
indices = '\U0001F10B' + \
    '\U00002780\U00002781\U00002782\U00002783\U00002784\U00002785\U00002786\U00002787\U00002788\U00002789' + \
    '\U0000246A\U0000246B\U0000246C\U0000246D\U0000246E\U0000246F\U00002470\U00002471\U00002472\U00002473' + \
    '\U00003251\U00003252\U00003253\U00003254\U00003255\U00003256\U00003257\U00003258\U00003259\U0000325A' + \
    '\U0000325B\U0000325C\U0000325D\U0000325E\U0000325F\U000032B1\U000032B2\U000032B3\U000032B4\U000032B5' + \
    '\U000032B6\U000032B7\U000032B8\U000032B9\U000032BA\U000032BB\U000032BC\U000032BD\U000032BE\U000032BF' if show else tuple ( '%s' % i for i in range ( 51 ) ) 


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
## Number
number           = '\U00002116'           if show else '#'

# ==================================================
def the_sum  ( what ) : return '%s%s'   % ( sum_symbol , what ) 
def the_mean ( what ) : return '%s%s%s' % ( langle , what , rangle ) 
def the_rms  ( what ) : return '%s%d'   % ( rms_symbol , what ) 

# ==============================================================================
## Generate sequence of numerical labels
#  @code
#  for l in labels ( 10 , 'ABC' ) : ..
#  @endcode 
def labels ( N , labs = () )  :
    """ Generate sequence of numerical labels
    >>>  for l in labels ( 10 , 'ABC' ) : ..
    """ 
    assert isinstance ( N , int ) and 0 <= N , 'Invalid number of labels!'

    q = 0 
    for i , l in enumerate ( labs ) :
        if i < N :
            q += 1 
            yield l 
        else     : return

    for j in range ( q , len ( indices ) ) : 
        if j < N :
            q += 1
            c  = indices  [ j ]
            if show and q < 22 : yield c + ' '
            else               : yield c 
            
        else     : return
        
    for k in range ( q , N ) : yield '%d' % k
    
# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

        
# =============================================================================
##                                                                     The END 
# =============================================================================
