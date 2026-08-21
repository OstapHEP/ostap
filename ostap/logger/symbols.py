
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some useful Unicode symbols
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
# =============================================================================
""" Module with some useful unicode symbols
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    'approximate'         , ## ≈
    'arrow_down'          , ## ↓
    'arrow_left'          , ## ←
    'arrow_right'         , ## →
    'arrow_rightleft'     , ## ↔
    'arrow_up'            , ## ↑
    'arrows_all'          , ## ←↖↑↗→↘↓↙
    'asterisk'            , ## ✶️
    'axe'                 , ## 🪓  (wide)
    'brain'               , ## 🧠  (wide)
    'branch'              , ## ⸙
    'cabinet'             , ## 🗄️  (wide)
    'chain'               , ## ⛓️
    'checked_no'          , ## ❌  (wide)
    'checked_yes'         , ## ✅  (wide)
    'chi2'                , ## 𝛘²
    'chi2ndf'             , ## 𝛘²/ndf
    'chisq'               , ## 𝛘²
    'clock'               , ## 🕐  (wide)
    'clock_ticks'         , ## 🕜 2-колон...
    'delta_symbol'        , ## Δ
    'difference'          , ## ⊻
    'dispersion_sym'      , ## σ²
    'ditto'               , ## 〃
    'document'            , ## 🗎️
    'efficiency'          , ## ε
    'ellipsis'            , ## …
    'enough'              , ## ∃
    'equivalent'          , ## ≡
    'exclusive_or'        , ## ⊻
    'finish'              , ## 🏁  (wide)
    'folder'              , ## 📂  (wide)
    'frame'               , ## 🖼️  (wide)
    'gear'                , ## ⚙️
    'graph'               , ## 📈  (wide)
    'greater_or_equal'    , ## ≥
    'greek_lower_alpha'   , ## α
    'greek_lower_beta'    , ## β
    'greek_lower_chi'     , ## χ
    'greek_lower_delta'   , ## δ
    'greek_lower_epsilon' , ## ε
    'greek_lower_eta'     , ## η
    'greek_lower_gamma'   , ## γ
    'greek_lower_iota'    , ## ι
    'greek_lower_kappa'   , ## κ
    'greek_lower_lambda'  , ## λ
    'greek_lower_mu'      , ## μ
    'greek_lower_nu'      , ## ν
    'greek_lower_omega'   , ## ω
    'greek_lower_omicron' , ## ο
    'greek_lower_phi'     , ## φ
    'greek_lower_pi'      , ## π
    'greek_lower_psi'     , ## ψ
    'greek_lower_rho'     , ## ρ
    'greek_lower_sigma'   , ## σ
    'greek_lower_tau'     , ## τ
    'greek_lower_theta'   , ## θ
    'greek_lower_xi'      , ## ξ
    'greek_lower_ypsilon' , ## υ
    'greek_lower_zeta'    , ## ζ
    'greek_upper_alpha'   , ## Α
    'greek_upper_beta'    , ## Β
    'greek_upper_chi'     , ## Χ
    'greek_upper_delta'   , ## Δ
    'greek_upper_epsilon' , ## Ε
    'greek_upper_eta'     , ## Η
    'greek_upper_gamma'   , ## Γ
    'greek_upper_iota'    , ## Ι
    'greek_upper_kappa'   , ## Κ
    'greek_upper_lambda'  , ## Λ
    'greek_upper_mu'      , ## Μ
    'greek_upper_nu'      , ## Ν
    'greek_upper_omega'   , ## Ω
    'greek_upper_omicron' , ## Ο
    'greek_upper_phi'     , ## Φ
    'greek_upper_pi'      , ## Π
    'greek_upper_psi'     , ## Ψ
    'greek_upper_rho'     , ## Ρ
    'greek_upper_sigma'   , ## Σ
    'greek_upper_tau'     , ## Τ
    'greek_upper_theta'   , ## Θ
    'greek_upper_xi'      , ## Ξ
    'greek_upper_ypsilon' , ## Υ
    'greek_upper_zeta'    , ## Ζ
    'hammer_and_wrench'   , ## 🛠️  (wide)
    'hand_ok'             , ## 👌  (wide)
    #
    'hebrew_aleph'        , ## א
    'hebrew_bet'          , ## ב
    'hebrew_dalet'        , ## ד
    'hebrew_gimel'        , ## ג
    #
    'histogram'           , ## 📊  (wide)
    'indices'             , ## ⓿➊➋...
    'intersection'        , ## ⋂
    'iteration'           , ## 々
    'kitchen_knife'       , ## ➖
    'langle'              , ## 〈
    'leaves'              , ## 🍃  (wide)
    'less_or_equal'       , ## ≤
    'light_bulb'          , ## 💡  (wide)
    'likelihood'          , ## ℒ  (wide)
    'minus_plus'          , ## ∓
    'mountain'            , ## ⛰️  (wide)
    'much_greater'        , ## ≫
    'much_less'           , ## ≪
    'not_equal'           , ## ≠
    'number'              , ## №
    'oil_drum'            , ## 🛢️  (wide)
    'palette'             , ## 🎨  (wide)
    'permille'            , ## ‰
    'plus_minus'          , ## ±
    'question_mark'       , ## ❓  (wide)
    'ram'                 , ## 🐏  (wide)
    'rangle'              , ## 〉
    'rms_symbol'          , ## σ
    'runner'              , ## 🏃  (wide)
    'same'                , ## ≡
    'scissors'            , ## ✂️
    'script_A'            , ## 𝒜
    'script_B'            , ## ℬ
    'script_E'            , ## ℰ
    'script_F'            , ## ℱ
    'script_H'            , ## ℋ
    'script_L'            , ## ℒ  (wide)
    'script_M'            , ## ℳ
    'script_P'            , ## 𝒫
    'script_R'            , ## ℛ
    'script_l'            , ## 𝓁
    'script_map'          , ## mapping dict
    'show'                , ## bool
    'similar'             , ## ∼
    'size'                , ## ⌀
    'squared_ok'          , ## 🆗  (wide)
    'subscript_A'         , ## ᴀ
    'subscript_C'         , ## ᴄ
    'subscript_K'         , ## ᴋ
    'subscript_a'         , ## ₐ
    'subscript_c'         , ## ꜀
    'subscript_k'         , ## ₖ
    'sum_symbol'          , ## ∑
    'superscript_map'     , ## mapping dict
    'symmetry'            , ## ⌯
    'tape'                , ## ✂️
    'tape_cartridge'      , ## 🖭️
    'thumb_down'          , ## 👎  (wide)
    'thumb_up'            , ## 👍  (wide)
    'times'               , ## ⨯
    'toys'                , ## 🧸  (wide)
    'tree'                , ## 🌴  (wide)
    'union'               , ## ⋃
    'variance_sym'        , ## σ²
    'weight_lifter'       , ## 🏋️  (wide)
    'weight_scale'        , ## ⚖️  (wide)
    'weierstrass_p'       , ##  ℘
    'wrench'              , ## 🔧  (wide)
    # -------------------------------------------------------------------------
    # functions & generators
    # -------------------------------------------------------------------------
    'labels'              , ## generator
    'to_script'           , ## func
    'the_mean'            , ## func
    'the_rms'             , ## func
    'the_sum'             , ## func
)
# ===========================================================================
from   ostap.utils.basic import has_unicode
from   ostap.core.config import show_unicode
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger    import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.logger.symbols' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## show unicode symbols? 
show = show_unicode and has_unicode ()

checked_yes      = '✅ ' if show else "+"
checked_no       = '❌ ' if show else "-"
question_mark    = '❓ ' if show else "?"
hand_ok          = '👌 ' if show else 'ok'
squared_ok       = '🆗 ' if show else 'ok'
thumb_up         = '👍 ' if show else '+'
thumb_down       = '👎 ' if show else '-'
clock            = '🕐 ' if show else '' 
ram              = '🐏 ' if show else ''
runner           = '🏃 ' if show else ''
finish           = '🏁 ' if show else ''

arrow_left       = '←'  if show else '<-'
arrow_right      = '→'  if show else '->'
arrow_rightleft  = '↔'  if show else '<->'
arrow_up         = '↑'  if show else '|'
arrow_down       = '↓'  if show else '|'

arrows_all       = '←↖↑↗→↘↓↙' if show else ( '<-' , '\\' , '|' , '/' , '->' , '\\' , '|' , '/' )

langle           = '〈'  if show else '<'
rangle           = '〉'  if show else '>'
ellipsis         = '…'  if show else '...'
same             = '≡'  if show else 'same'

clock_ticks      = '🕐🕜🕑🕝🕒🕞🕓🕟🕔🕠🕕🕡🕖🕢🕗🕣🕘🕤🕙🕟🕚🕦🕛🕧' if show else '|/-\\'
clock_ticks      = tuple ( clock_ticks )

times            = '⨯'  if show else 'x'
plus_minus       = '±'  if show else '+/-'
minus_plus       = '∓'  if show else '-/+'
ditto            = '〃' if show else '-//-'

tree             = '🌴 ' if show else ''
chain            = '⛓ '  if show else '' 
branch           = '⸙'  if show else '' 
leaves           = '🍃 ' if show else '' 
cabinet          = '🗄 ' if show else '' 
frame            = '🖼 ' if show else '' 
histogram        = '📊 ' if show else ''
graph            = '📈 ' if show else '' 
palette          = '🎨 ' if show else '' 
document         = '🗎'  if show else '' 
tape             = '✂'  if show else '' 
tape_cartridge   = '🖭' if show else '' 
folder           = '📂 ' if show else '' 
light_bulb       = '💡 ' if show else '' 

less_or_equal    = '≤'   if show else '<='
greater_or_equal = '≥'   if show else '=>'
much_less        = '≪'   if show else '<<'
much_greater     = '≫'   if show else '>>'
equivalent       = '≡'   if show else '='
similar          = '∼'   if show else '~'
approximate      = '≈'   if show else '~='
not_equal        = '≠'   if show else '!='
weight_lifter    = '🏋'  if show else ''
weight_scale     = '⚖' if show else '' 

scissors         = '✂'  if show else '' 
oil_drum         = '🛢 ' if show else '' 
brain            = '🧠 ' if show else ''
kitchen_knife    = '➖' if show else ''
axe              = '🪓 ' if show else ''
gear             = '⚙ '  if show else ''
wrench           = '🔧 ' if show else ''
hammer_and_wrench= '🛠 ' if show else ''

union            = '⋃'  if show else ''
intersection     = '⋂'  if show else ''
exclusive_or     = '⊻'  if show else '^'
difference       = '⊻'  if show else '-'

iteration        = '々' if show else ''
permille         = '‰'  if show else '/1000'
size             = '⌀'  if show else 'size'

## indices: circled numbers from 0 to 50 (inclusive)
indices = '⓿➊➋➌➍➎➏➐➑➒⓫⓬⓯⓰⓱⓲⓳⓴㉑㉒㉓㉔㉕㉖㉗㉘㉙㉚㉛㉜㉝㉞㉟㊱㊲㊳㊴㊵㊶㊷㊸㊹㊺㊻㊼㊽㊾㊿' if show else tuple ( '%s' % i for i in range ( 51 ) ) 
indices = tuple ( indices )

## capital Greek Sigma 
sum_symbol       = '∑'   if show else 'sum'
## lowercase Greek sigma 
rms_symbol       = 'σ'   if show else 'rms'
## squared lower case Greek sigma 
dispersion_sym   = 'σ²'  if show else 'D'
## squared lower case Greek sigma 
variance_sym     = 'σ²'  if show else 'var'
## Delta symbol 
delta_symbol     = 'Δ'   if show else 'delta'
## Number
number           = '№'   if show else '#'
## chi2
chi2             = '𝛘²'  if show else 'chi2'
## chi2
chisq            = chi2 
## chi2/ndf 
chi2ndf          = '%s/ndf' % chi2 

## symmetric
symmetry         = '⌯'   if show else 'sym'
enough           = '∃'   if show else 'enough'
mountain         = '⛰ '  if show else 'peak'

## Lowercase Greek letters 
greek_lower_alpha     = 'α' if show else 'alpha'
greek_lower_beta      = 'β' if show else 'beta'
greek_lower_gamma     = 'γ' if show else 'gamma'
greek_lower_delta     = 'δ' if show else 'delta'
greek_lower_epsilon   = 'ε' if show else 'epsilon'
greek_lower_zeta      = 'ζ' if show else 'zeta'
greek_lower_eta       = 'η' if show else 'eta'
greek_lower_theta     = 'θ' if show else 'theta'
greek_lower_iota      = 'ι' if show else 'iota'
greek_lower_kappa     = 'κ' if show else 'kappa'
greek_lower_lambda    = 'λ' if show else 'lambda'
greek_lower_mu        = 'μ' if show else 'mu'
greek_lower_nu        = 'ν' if show else 'nu'
greek_lower_xi        = 'ξ' if show else 'xi'
greek_lower_omicron   = 'ο' if show else 'omicron'
greek_lower_pi        = 'π' if show else 'pi'
greek_lower_rho       = 'ρ' if show else 'rho'
greek_lower_sigma     = 'σ' if show else 'sigma'
greek_lower_tau       = 'τ' if show else 'tau'
greek_lower_ypsilon   = 'υ' if show else 'ypsilon'
greek_lower_phi       = 'φ' if show else 'phi'
greek_lower_chi       = 'χ' if show else 'chi'
greek_lower_psi       = 'ψ' if show else 'psi'
greek_lower_omega     = 'ω' if show else 'omega'

## Uppercase Greek letters 
greek_upper_alpha     = 'Α' if show else 'A'
greek_upper_beta      = 'Β' if show else 'B'
greek_upper_gamma     = 'Γ' if show else 'Gamma'
greek_upper_delta     = 'Δ' if show else 'Delta'
greek_upper_epsilon   = 'Ε' if show else 'E'
greek_upper_zeta      = 'Ζ' if show else 'Z'
greek_upper_eta       = 'Η' if show else 'H'
greek_upper_theta     = 'Θ' if show else 'Theta'
greek_upper_iota      = 'Ι' if show else 'I'
greek_upper_kappa     = 'К' if show else 'K'
greek_upper_lambda    = 'Λ' if show else 'Lambda'
greek_upper_mu        = 'Μ' if show else 'M'
greek_upper_nu        = 'Ν' if show else 'N'
greek_upper_xi        = 'Ξ' if show else 'Xi'
greek_upper_omicron   = 'Ο' if show else 'O'
greek_upper_pi        = 'Π' if show else 'Pi'
greek_upper_rho       = 'Ρ' if show else 'P'
greek_upper_sigma     = 'Σ' if show else 'Sigma'
greek_upper_tau       = 'Τ' if show else 'T'
greek_upper_ypsilon   = 'Υ' if show else 'Y'
greek_upper_phi       = 'Φ' if show else 'Phi'
greek_upper_chi       = 'Χ' if show else 'X'  
greek_upper_psi       = 'Ψ' if show else 'Psi'
greek_upper_omega     = 'Ω' if show else 'Omega'


## Hebrew letters (Mathematical Alphanumeric LTR symbols - won't break Bidi layout)
hebrew_aleph          = 'ℵ' if show else 'aleph'   # U+2135
hebrew_bet            = 'ℶ' if show else 'bet'     # U+2136
hebrew_gimel          = 'ℷ' if show else 'gimel'   # U+2137
hebrew_dalet          = 'ℸ' if show else 'dalet'   # U+2138


## Script / Calligraphic symbols
script_A              = '𝒜' if show else 'A'
script_B              = 'ℬ' if show else 'B'
script_E              = 'ℰ' if show else 'E'
script_F              = 'ℱ' if show else 'F'
script_H              = 'ℋ' if show else 'H'
script_L              = 'ℒ ' if show else 'L'
script_M              = 'ℳ' if show else 'M'
script_P              = '𝒫' if show else 'P'
script_R              = 'ℛ' if show else 'R'
script_l              = '𝓁' if show else 'l'

## likelihood: script L
likelihood            = script_L

subscript_a           = 'ₐ' if show else 'a'
subscript_c           = '꜀' if show else 'c'
subscript_k           = 'ₖ' if show else 'k'
subscript_A           = 'ᴀ' if show else 'A'
subscript_C           = 'ᴄ' if show else 'C'
subscript_K           = 'ᴋ' if show else 'K'

## use epsilon symbol for efficiency
efficiency       = greek_lower_epsilon if show else 'eff'

## toys = teddy bear 
toys                  = '🧸 ' if show else 'toys'

## star/convolution operator
asterisk              = '✶' if show else '*'
weierstrass_p         = '℘' if show else 'p'  # U+2118 WEIERSTRASS ELLIPTIC FUNCTION

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

# ============================================================================
# Mapping tables
superscript_map = str.maketrans({
    '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
    '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
    '-': '⁻', '+': '⁺', '.': 'ˑ', 'e': 'ᵉ', 'E': 'ᵉ'
})

# Mapping table for Script / Calligraphic chars (Latin UPPER + lower)
_latin_norm   = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
_latin_script = "𝒜ℬ𝒞𝒟ℰℱ𝒢ℋℐ𝒥𝒦ℒℳ𝒩𝒪𝒫𝒬ℛ𝒮𝒯𝒰𝒱𝒲𝒳𝒴𝒵𝒶𝒷𝒸𝒹ℯ𝒻ℊ𝒽𝒾𝒿𝓀𝓁𝓂𝓃ℴ𝓅𝓆𝓇𝓈𝓉𝓊𝓋𝓌𝓍𝓎𝓏"
script_map    = str.maketrans(_latin_norm, _latin_script) if show else str.maketrans("", "")

def to_script(text):
    """Convert standard ASCII string to Mathematical Script Unicode text
    >>> to_script("LHCb")
    'ℒℋ𝒞𝒷'
    """
    return text.translate(script_map) if show else text


# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    globs  = globals ().copy() 
    header = "Name" , "Value" , 'Width'
    import ostap.logger.table as T

    rows = [] 
    for name in __all__  :
        symb = globs.get ( name , None )
        if symb is None or not isinstance ( symb , str ) : continue        
        row = name , symb , '%d' % T.visible_width ( symb )
        rows.append ( row ) 
    
    title  = "Symbols"
    table  = T.table ( rows , title = title , prefix = '# ' , alignment = 'lcc' ) 
    logger.info ( "%s:\n%s" % ( title , table ) )

# =============================================================================
##                                                                     The END 
# =============================================================================

