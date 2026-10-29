#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some useful Unicode symbols
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
# =============================================================================
""" Module with some useful unicode symbols (with wide symbol definitions)
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
    'asterisk'            , ## ✶
    'asymmetry'           , ## ⌿
    'attention'           , ## 🚨
    'axe'                 , ## 🪓  (wide)
    ##
    'black_arrow_left'    , ## ⬅  
    'black_arrow_right'   , ## ⮕  
    'black_arrow_up'      , ## ⬆  
    'black_arrow_down'    , ## ⬇  
    ##
    'brain'               , ## 🧠  (wide)
    'branch'              , ## ⸙
    'brickwall'           , ## 🧱
    'bulls_eye'           , ## 🎯
    ## 
    'cat_face'            , ## 🐱
    'cabinet'             , ## 🗄️  (wide)
    'chain'               , ## ⛓️
    'checked_no'          , ## ❌  (wide)
    'checked_yes'         , ## ✅  (wide)
    'chi2'                , ## 𝛘²
    'chi2ndf'             , ## 𝛘²/ndf
    'chisq'               , ## 𝛘²
    'clock'               , ## 🕐  (wide)
    'clock_ticks'         , ## 🕜 2-колон...
    'clubs'               , ## ♣
    'column'              , ## ⋮   
    'concurrency'         , ## 🔀
    'conifer_tree'        , ## 🌲
    'christmas_tree'      , ## 🎄
    ## 
    'delta_symbol'        , ## Δ
    'deciduous_tree'      , ## 🌳    
    'descartes_product'   , ## ⨂
    'dice'                , ## '🎲'
    'difference'          , ## ⊻
    'dispersion_sym'      , ## σ²
    'ditto'               , ## 〃
    'document'            , ## 🗎️
    ## 
    'efficiency'          , ## ε
    'ellipsis'            , ## …
    'enough'              , ## ∃
    'equivalent'          , ## ≡
    'exclusive_or'        , ## ⊻
    ## 
    'finish'              , ## 🏁  (wide)
    'flashlight'          , ## 🔦
    'folder'              , ## 📂  (wide)
    'frame'               , ## 🖼️  (wide)
    ##
    'game_die'            , ## 🎲 
    'gear'                , ## ⚙️
    'graph'               , ## 📈  (wide)
    ##
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
    ##
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
    'hiking_boot'         , ## 🥾
    'histogram'           , ## 📊  (wide)
    #
    'hourglass'           , ## ⏳    
    ##
    'infinity'            , ## ∞
    'infinity_pos'        , ## +∞   
    'infinity_neg'        , ## -∞   
    # #   
    'indices'             , ## ⓿➊➋...
    'intersection'        , ## ⋂
    'iteration'           , ## 々
    ##
    'kitchen_knife'       , ## ➖
    ## 
    'langle'              , ## 〈
    'label'               , ## 🏷️
    'leaf'                , ## 🍁
    'leaves'              , ## 🍃  (wide)
    'less_or_equal'       , ## ≤
    'light_bulb'          , ## 💡  (wide)
    'likelihood'          , ## ℒ  (wide)
    'lightning'           , ## ⚡
    'lock_closed'         , ## 🔒
    'lock_open'           , ## 🔓
    ## 
    'matrix'              , ## ▦
    'minus_plus'          , ## ∓
    'mountain'            , ## ⛰️  (wide)
    'much_greater'        , ## ≫
    'much_less'           , ## ≪
    ## 
    'nabla'               , ## ∇
    'not_equal'           , ## ≠
    'number'              , ## №
    'oil_drum'            , ## 🛢️  (wide)
    ##
    'package'             , ## 📦    
    'palette'             , ## 🎨  (wide)
    'palm'                , ## '🌴'      
    'permille'            , ## ‰
    'permutation'         , ## ⤭
    'postal_horn'         , ## 📯
    'plus_minus'          , ## ±
    'proportion'          , ## ∷
    'quadratic_sum'       , ## ⨁ 
    'question_mark'       , ## ❓  (wide)
    'ram'                 , ## 🐏  (wide)
    'rangle'              , ## 〉
    'repeat'              , ##  🔁     
    'rms_symbol'          , ## σ
    'rocket'              , ## 🚀
    'ruler'               , ## 📏
    'runner'              , ## 🏃  (wide)
    'same'                , ## ≡
    'scale'               , ## ⚖️  (wide)
    'scissors'            , ## ✂️
    ##
    'script_A'            , ## 𝒜
    'script_B'            , ## ℬ
    'script_E'            , ## ℰ
    'script_F'            , ## ℱ
    'script_H'            , ## ℋ
    'script_L'            , ## ℒ  (wide)
    'script_M'            , ## ℳ
    'script_P'            , ## 🫫 / 🫥
    'script_R'            , ## ℛ
    'script_map'          , ## mapping dict
    'show'                , ## bool
    'shuffle'             , ## 🔀    
    'similar'             , ## ∼
    'size'                , ## ⌀
    'squared_ok'          , ## 🆗  (wide)
    'sparkles'            , ## ✨
    'subscript_A'         , ## ᴀ
    'subscript_C'         , ## ᴄ
    'subscript_K'         , ## ᴋ
    'subscript_a'         , ## ₐ
    'subscript_c'         , ## ꜀
    'subscript_k'         , ## ₖ    
    'sub_min'             , ## ₘᵢₙ
    'sub_max'             , ## ₘₐₓ
    'sub_mean'            , ## ₘₑₐₙ
    'sub_rms'             , ## ᵣₘₛ
    'sup_eff'             , ## ᵉᶠᶠ
    'sum_symbol'          , ## ∑    
    'superscript_map'     , ## mapping dict
    'symmetry'            , ## ⌯
    'tag'                 , ## ❖
    'tape'                , ## ➿
    'tape_cartridge'      , ## 🖭️
    'tensor_product'      , ## ⨂
    'thumb_down'          , ## 👎  (wide)
    'thumb_up'            , ## 👍  (wide)
    'times'               , ## × 
    'times_colored'       , ## ❌ (Красивый цветной крестик)
    'triangular_flag'     , ## 🚩
    'toys'                , ## 🧸  (wide)
    'tree'                , ## 🌳 (wide)
    'trophy'              , ## 🏆
    ##
    'union'                  , ## ⋃
    ## 
    'variance_sym'           , ## σ²
    ## 
    'warning'                , ## ⚠️
    'weight_lifter'          , ## 🏋️  (wide)
    'weight_scale'           , ## ⚖️  (wide)
    'weierstrass_p'          , ##  ℘
    ##
    'white_arrow_left'       , ## ⇦
    'white_arrow_right'      , ## ⇨
    'white_arrow_up'         , ## ⇧
    'white_arrow_down'       , ## ⇩
    'white_arrow_up_left'    , ## ⇖ 
    'white_arrow_up_right'   , ## ⇗ 
    'white_arrow_down_left'  , ## ⇙ 
    'white_arrow_down_right' , ## ⇘
    'white_arrows_all'       , 
    ##
    'wood'                , ## 🪵
    'wrench'              , ## 🔧  (wide)
    ## 
    'T2'                  , ## 🯯² (Note: handled or wrapped)
    'W2'                  , ## 🯰² 
    # -------------------------------------------------------------------------
    # functions & generators & collections
    # -------------------------------------------------------------------------
    'labels'              , ## generator
    'the_mean'            , ## func
    'the_rms'             , ## func
    'the_sum'             , ## func
    ##
    'to_script'           , ## convert text to Script 
    'to_roman'            , ## convert text to Roman 
    'to_italic'           , ## convert text to Italic
    'to_double_struck'    , ## convert text to Double-Dtruck 
    'to_fraktur'          , ## convert text to Fractur 
    'to_sans_serif'       , ## convert text to Sans-Serif  
    'to_bold'             , ## convert text to Bold 
    'to_bold_italic'      , ## convert text to Bold Italic   
    ## 
    'remove_invisible'    , ## remove variantion selectors 
    'exceptional_symbols' , ## set of exceptional wide symbols
    'visible_width'       , ## visible width of the character string 
)
# ===========================================================================
from   ostap.utils.basic import has_unicode
from   ostap.core.config import show_unicode
import sys, re, wcwidth, unicodedata
# =============================================================================
## show unicode symbols? 
show = show_unicode and has_unicode ()

checked_yes      = '✅' if show else "+"
checked_no       = '❌' if show else "-"
question_mark    = '❓' if show else "?"
hand_ok          = '👌' if show else 'ok'
squared_ok       = '🆗' if show else 'ok'
thumb_up         = '👍' if show else '+'
thumb_down       = '👎' if show else '-'
clock            = '🕐' if show else '' 
ram              = '🐏' if show else ''
runner           = '🏃' if show else ''
finish           = '🏁' if show else ''


arrow_left       = '←'  if show else '<-'
arrow_right      = '→'  if show else '->'
arrow_rightleft  = '↔'  if show else '<->'
arrow_up         = '↑'  if show else '|'
arrow_down       = '↓'  if show else '|'
arrow_up_left    = '↖'  if show else '\\'
arrow_up_right   = '↗'  if show else '/'
arrow_down_left  = '↙'  if show else '/'
arrow_down_right = '↘'  if show else '\\'

arrows_all       = ( arrow_left       , 
                     arrow_up_left    , 
                     arrow_up         , 
                     arrow_up_right   , 
                     arrow_right      , 
                     arrow_down_right , 
                     arrow_down       , 
                     arrow_down_left  ) 
if show : arrows_all = ''.join ( arrows_all ) 

black_arrow_left      = '⬅'  if show else arrow_left 
black_arrow_right     = '⮕'  if show else arrow_right 
black_arrow_up        = '⬆'  if show else arrow_up 
black_arrow_down      = '⬇'  if show else arrow_down 

white_arrow_left      = '⇦'  if show else arrow_left 
white_arrow_right     = '⇨'  if show else arrow_right 
white_arrow_up        = '⇧'  if show else arrow_up 
white_arrow_down      = '⇩'  if show else arrow_down 
white_arrow_up_left   = '⇖'  if show else arrow_up_left 
white_arrow_up_right  = '⇗'  if show else arrow_up_right 
white_arrow_down_left = '⇙'  if show else arrow_down_left 
white_arrow_down_right= '⇘'  if show else arrow_down_right 

white_arrows_all       = ( white_arrow_left       , 
                           white_arrow_up_left    , 
                           white_arrow_up         , 
                           white_arrow_up_right   , 
                           white_arrow_right      , 
                           white_arrow_down_right , 
                           white_arrow_down       , 
                           white_arrow_down_left  ) 
if show : white_arrows_all = ''.join ( white_arrows_all ) 


langle           = '〈'  if show else '<'
rangle           = '〉'  if show else '>'
ellipsis         = '…'  if show else '...'
same             = '≡'  if show else 'same'

clock_ticks      = '🕐🕜🕑🕝🕒🕞🕓🕟🕔🕠🕕🕡🕖🕢🕗🕣🕘🕤🕙🕟🕚🕦🕛🕧' if show else '|/-\\'
clock_ticks      = tuple ( clock_ticks )

times            = '×'  if show else 'x'   
times_colored    = '❌' if show else times 
plus_minus       = '±'  if show else '+/-'
minus_plus       = '∓'  if show else '-/+'
ditto            = '〃' if show else '-//-'

palm              = '🌴'  if show else ''       
conifer_tree      = '🌲'  if show else '' 
christmas_tree    = '🎄'  if show else ''      
deciduous_tree    = '🌳'  if show else ''      
tree              = deciduous_tree 

chain            = '⛓️'   if show else ''        # Добавлен вариационный селектор для стабильности
branch           = '⸙'   if show else ''
leaf             = '🍁'  if show else '' 
leaves           = '🍃'  if show else '' 
cabinet          = '🗄️'  if show else '' 
frame            = '🖼️'  if show else '' 
histogram        = '📊'  if show else ''
graph            = '📈'  if show else '' 
palette          = '🎨'  if show else '' 
document         = '🗎️'  if show else '' 
tape             = '➿'  if show else '' 
tape_cartridge   = '🖭️' if show else '' 
folder           = '📂' if show else '' 
light_bulb       = '💡' if show else '' 

less_or_equal    = '≤'   if show else '<='
greater_or_equal = '≥'   if show else '=>'
much_less        = '≪'   if show else '<<'
much_greater     = '≫'   if show else '>>'
equivalent       = '≡'   if show else '='
similar          = '∼'   if show else '~'
approximate      = '≈'   if show else '~='
not_equal        = '≠'   if show else '!='
weight_lifter    = '🏋️'  if show else ''

scale            = '⚖️'  if show else '' 
weight_scale     = scale  

scissors         = '✂️' if show else '' 
oil_drum         = '🛢️' if show else ''        # Добавлен селектор
brain            = '🧠' if show else ''
kitchen_knife    = '🔪' if show else '➖'      # Заменено на красивый нож вместо минуса
axe              = '🪓' if show else ''
gear             = '⚙️'  if show else ''        # Добавлен селектор
wrench           = '🔧' if show else ''
hammer_and_wrench= '🛠️'  if show else ''
package          = '📦' if show else ''     

union            = '⋃'  if show else ''
intersection     = '⋂'  if show else ''
exclusive_or     = '⊻'  if show else '^'
difference       = '⊻'  if show else '-'

iteration        = '々' if show else ''
permille         = '‰'  if show else '/1000'
size             = '⌀'  if show else 'size'

## indices: circled numbers from 0 to 50 (inclusive)
indices = '⓪①②③④⑤⑥⑦⑧⑨⑩⑪⑫⑬⑭⑮⑯⑰⑱⑲⑳㉑㉒㉓㉔㉕㉖㉗㉘㉙㉚㉛㉜㉝㉞㉟㊱㊲㊳㊴㊵㊶㊷㊸㊹㊺㊻㊼㊽㊾㊿' if show else tuple ( '%s' % i for i in range ( 51 ) ) 
indices = tuple ( indices )

## capital Greek Sigma 
sum_symbol       = '∑'   if show else 'sum'

tensor_product    = '⨂'  if show else 'x'
descartes_product =  tensor_product 
quadratic_sum     = '⨁'  if show else '(+)'

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

nabla            = '∇'   if show else ''
proportion       = '∷'   if show else ''
clubs            = '♣'  if show else ''
postal_horn      = '📯'  if show else '' 

## symmetric

symmetry         = '⌯'   if show else 'sym'
asymmetry        = '⌿'   if show else 'asym'

enough           = '∃'   if show else 'enough'
mountain         = '⛰️'   if show else 'peak'   # Добавлен селектор

concurrency      = '🔀'  if show  else '||'    
matrix           = '▦'  if show  else ''
permutation      = '⤭'   if show  else ''
column           = '⋮'    if show  else ''

hiking_boot      = '🥾' if show  else 'bootstrap'
game_die         = '🎲' if show  else ''
shuffle          = '🔀' if show  else '' 

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
greek_lower_phi       = 'φ' if show else 'phi'
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
greek_upper_eta       = 'Η' if show else 'H'
greek_upper_zeta      = 'Ζ' if show else 'Z'
greek_upper_theta     = 'Θ' if show else 'Theta'
greek_upper_iota      = 'Ι' if show else 'I'
greek_upper_kappa     = 'К' if show else 'K'
greek_upper_lambda    = 'Λ' if show else 'Lambda'
greek_upper_mu        = 'Μ' if show else 'M'
greek_upper_nu        = 'Ν' if show else 'N'
greek_upper_xi        = 'Ξ' if show else 'Xi'
greek_upper_omicron   = 'Ο' if show else 'O'
greek_upper_phi       = 'Φ' if show else 'Phi'
greek_upper_pi        = 'Π' if show else 'Pi'
greek_upper_rho       = 'Ρ' if show else 'P'
greek_upper_sigma     = 'Σ' if show else 'Sigma'
greek_upper_tau       = 'Τ' if show else 'T'
greek_upper_ypsilon   = 'Υ' if show else 'Y'
greek_upper_phi       = 'Φ' if show else 'Phi'
greek_upper_chi       = 'Χ' if show else 'X'  
greek_upper_psi       = 'Ψ' if show else 'Psi'
greek_upper_omega     = 'Ω' if show else 'Omega'


## Hebrew letters
hebrew_aleph          = 'ℵ' if show else 'aleph'   
hebrew_bet            = 'ℶ' if show else 'bet'     
hebrew_gimel          = 'ℷ' if show else 'gimel'   
hebrew_dalet          = 'ℸ' if show else 'dalet'   


## Script / Calligraphic symbols
script_A              = '𝒜' if show else 'A'
script_B              = 'ℬ' if show else 'B'
script_E              = 'ℰ'  if show else 'E'
script_F              = 'ℱ' if show else 'F'
script_H              = 'ℋ' if show else 'H'
script_L              = 'ℒ' if show else 'L'
script_M              = 'ℳ' if show else 'M'
script_P              = '𝒫'  if show else 'P'
script_R              = 'ℛ' if show else 'R'
script_l              = '𝓁'  if show else 'l'

likelihood            = script_L

subscript_a           = 'ₐ' if show else 'a'
subscript_c           = '꜀' if show else 'c'
subscript_k           = 'ₖ' if show else 'k'
subscript_A           = 'ᴀ' if show else 'A'
subscript_C           = 'ᴄ' if show else 'C'
subscript_K           = 'ᴋ' if show else 'K'

efficiency       = greek_lower_epsilon if show else 'eff'

toys                  = '🧸'  if show else 'toys'
repeat                = '🔁'  if show else ''     

asterisk              = '✶' if show else '*'
weierstrass_p         = '℘' if show else 'p'  

infinity              = "\N{INFINITY}" if show else 'inf'
infinity_pos          = "+%s" % infinity 
infinity_neg          = "-%s" % infinity

trophy                = '🏆'  if show else ''
sparkles              = '✨'  if show else ''
bulls_eye             = '🎯'  if show else '' 
triangular_flag       = '🚩'  if show else ''

cat_face              = '🐱'  if show else ''
rocket                = '🚀'  if show else ''
flashlight            = '🔦'  if show else '' 
brickwall             = '🧱'  if show else '' 
wood                  = '🪵'  if show else '' 
ruler                 = '📏'  if show else '' 

hourglass             = '⏳'  if show else ''
lightning             = '⚡'  if show else ''
lock_closed           = '🔒'  if show else ''
lock_open             = '🔓'  if show else ''
warning               = '⚠️'  if show else ''
attention             = '🚨'  if show else '' 

dice                  = '🎲'  if show else ''

# ==================================================
def the_sum  ( what ) : return '%s%s'   % ( sum_symbol , what ) 
def the_mean ( what ) : return '%s%s%s' % ( langle , what , rangle ) 
def the_rms  ( what ) : return '%s%d'   % ( rms_symbol , what ) 

# ==============================================================================
def labels ( N , labs = () )  :
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

    
label = '🏷️'  if show else ''
tag   = '❖'  if show else '' 

superscript_map = str.maketrans({
    '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
    '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
    '-': '⁻', '+': '⁺', '.': 'ˑ', 'e': 'ᵉ', 'E': 'ᵉ'
})

# ==============================================================================
# Text Styling Utilities (Unicode-based)
# ==============================================================================

# Base Latin alphabet for mapping
latin_norm            = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
# Character sets for various Unicode typography styles
script_letters        = "𝒜ℬ𝒞𝒟ℰℱ𝒢ℋℐ𝒥𝒦ℒℳ𝒩𝒪𝒫𝒬ℛ𝒮𝒯𝒰𝒱𝒲𝒳𝒴𝒵𝒶𝒷𝒸𝒹ℯ𝒻ℊ𝒽𝒾𝒿𝓀𝓁𝓂𝓃ℴ𝓅𝓆𝓇𝓈𝓉𝓊𝓋𝓌𝓍𝓎𝒵"
roman_letters         = "𝐀𝐁𝐂𝐃𝐄𝐅𝙶𝐇𝐈𝙹𝐊𝐋𝙼𝙽𝙾𝙿𝐐𝐑𝚂𝐓𝚄𝚅𝚆𝚇𝚈𝚉𝐚𝐛𝐜𝐝𝐞𝐟𝐠𝐡𝐢𝐣𝚔𝚕𝚖𝚗𝐨𝐩𝕢𝐫𝐬𝐭𝐮𝐯𝐰𝐱𝐲𝐳"
italic_letters        = "𝐴𝐵𝐶𝐷𝐸𝐹𝐺𝐻𝐼𝐽𝐾𝐿𝑀𝑁𝑂𝑃𝑄𝑅𝑆𝑇𝑈𝑉𝑊𝑋𝑌𝑍𝑎𝑏𝑐𝑑𝑒𝑓𝑔ℎ𝑖𝑗𝑘𝑙𝑚𝑛𝑜𝑝𝑞𝑟𝑠𝑡𝑢𝑣𝑤𝑥𝑦𝑧"
double_struck_letters = "𝔸𝔹ℂ𝔻𝔼𝔽𝔾ℍ𝕀𝕁𝕂𝕃𝕄ℕ𝕆ℙℚℝ𝕊𝕋𝕌𝕍𝕎𝕏𝕐ℤ𝕒𝕓𝕔𝕕𝕖𝔣𝕘𝕙𝕚𝕛𝕜𝕝𝕞𝕟𝕠𝕡𝕢𝕣𝕤𝕥𝕦𝚟𝕨𝕩𝚢𝕫"
fraktur_letters       = "𝔄𝔅ℭ𝔇𝔈𝔉𝔊ℌℑ𝔍𝔎𝔏𝔐𝔑𝔒𝔓𝔔ℜ𝔖𝔗𝘘𝔙𝔚𝔛𝔜ℨ𝔞𝔟𝔠𝔡𝔢𝔣𝔤𝔥𝔦𝔧𝔨𝔩𝔪𝔫𝔬𝔭𝔮𝔯𝔰𝔱𝔲𝔳𝔴𝔵𝔶𝔷"
sans_serif_letters    = "𝗔𝗕𝗖𝗗𝗘𝗙𝗚𝗛𝗜𝗝𝗞𝗟𝗠𝗡𝗢𝗣𝗤𝗥𝗦𝗧𝗨𝗩𝗪𝗫𝗬𝗭𝗮𝗯𝗰𝗱𝗲𝗳𝗴𝗵𝗶𝗷𝗸𝚕𝚖𝗻𝗼𝗽𝗾𝗿𝘀𝘁𝚞𝘃𝘄𝘅𝘆𝐳"
bold_italic_letters   = "".join ( chr ( 0x1D468 + i ) for i in range ( 26 ) ) + "".join ( chr ( 0x1D482 + i ) for i in range ( 26 ) )
bold_letters          = "".join ( chr ( 0x1D400 + i ) for i in range ( 26 ) ) + "".join ( chr ( 0x1D41A + i ) for i in range ( 26 ) )

# Translation maps (automatically disabled if unicode support is off)
script_map        = str.maketrans ( latin_norm  , script_letters        ) if show else str.maketrans ( "" , "" )
roman_map         = str.maketrans ( latin_norm  , roman_letters         ) if show else str.maketrans ( "" , "" )
italic_map        = str.maketrans ( latin_norm  , italic_letters        ) if show else str.maketrans ( "" , "" )
double_struck_map = str.maketrans ( latin_norm  , double_struck_letters ) if show else str.maketrans ( "" , "" )
fraktur_map       = str.maketrans ( latin_norm  , fraktur_letters       ) if show else str.maketrans ( "" , "" )
sans_serif_map    = str.maketrans ( latin_norm  , sans_serif_letters    ) if show else str.maketrans ( "" , "" )
bold_italic_map   = str.maketrans ( latin_norm  , bold_italic_letters   ) if show else str.maketrans ( "" , "" ) 
bold_map          = str.maketrans ( latin_norm  , bold_letters          ) if show else str.maketrans ( "" , "" )

# Formatting helper functions
def to_script        ( text ) : return text.translate ( script_map        ) if show else text
def to_roman         ( text ) : return text.translate ( roman_map         ) if show else text
def to_italic        ( text ) : return text.translate ( italic_map        ) if show else text
def to_double_struck ( text ) : return text.translate ( double_struck_map ) if show else text
def to_fraktur       ( text ) : return text.translate ( fraktur_map       ) if show else text
def to_sans_serif    ( text ) : return text.translate ( sans_serif_map    ) if show else text
def to_bold_italic   ( text ) : return text.translate ( bold_italic_map   ) if show else text 
def to_bold          ( text ) : return text.translate ( bold_map          ) if show else text


sub_min  = 'ₘᵢₙ'  if show else '-min'
sub_max  = 'ₘₐₓ'  if show else '-max'
sub_rms  = 'ᵣₘₛ'  if show else '-rms'
sub_mean = 'ₘₑₐₙ' if show else '-mean'
sup_eff  = 'ᵉᶠᶠ'  if show else '-eff'

T2              = '𝒯 ²' if show else 'T2'
W2              = '𝒲 ²' if show else 'W2'


# =============================================================================
## Single comprehensive regex for cleaning:
#  1. ANSI escape color codes
#  2. Variation selectors (U+FE00 - U+FE0F)
#  3. Zero-width spaces and bidi directional marks (U+200B - U+200F)
# =============================================================================
INVISIBLE_CHARS_RE = re.compile(
    r'('
    r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])'  # ANSI escape sequences (colors, styles)
    r'|[\uFE00-\uFE0F]'                        # Variation Selectors (VS1 - VS16)
    r'|[\u200B-\u200F]'                        # Zero-width spaces & Bidi marks
    r')'
)
# =============================================================================
## Remove all invisible formatting, zero-width, and variation characters
# =============================================================================
def remove_invisible ( text ) :
    """ Remove invisible control characters, variation selectors, 
    zero-width spaces, and directional marks from the text.
    """
    return '' if not text else INVISIBLE_CHARS_RE.sub ( '' , text ) 

# =============================================================================
## Unified dictionary of wide/multi-column symbols and their exact terminal widths.
#  Includes standard wide icons (width 2) and hardcoded overrides for composite emojis.
# =============================================================================
exceptional_symbols = {
    cabinet           : 2 ,
    chain             : 2 , 
    document          : 2 , 
    frame             : 2 , 
    gear              : 2 , 
    hammer_and_wrench : 2 , 
    label             : 2 , 
    mountain          : 2 , 
    oil_drum          : 2 , 
    scissors          : 2 , 
    tape_cartridge    : 2 ,
    warning           : 2 , 
    weight_lifter     : 2 ,
    weight_scale      : 2 , 
}
# =============================================================================
# Normalize keys by stripping variation selectors so they match cleaned text reliably
exceptional_symbols = { remove_invisible ( k ) : v  for k , v in exceptional_symbols.items() if k }

## conver to frozendict if/when possible
if ( 3 , 15 ) <= sys.version_info :  
    exceptional_symbols = frozendict ( exceptional_symbols )
    
# =============================================================================
## Visible width of the string 
#   Calculate the true visible width of a string in a terminal,
#   using single-pass clean_text sanitization and direct dictionary lookup 
#   for wide symbols.
def visible_width ( what ) :
    """ Calculate the true visible width of a string in a terminal,
    using single-pass clean_text sanitization and direct dictionary lookup 
    for wide symbols.
    """
    if not what : return 0
    
    # Sanitize string ( strips ANSI colors, variation selectors, zero-width & bidi chars)
    text = remove_invisible ( what )
    if not text: return 0
    
    width = 0
    
    # Iterate through each clean character directly
    for char in text:
        # Check if this single character has a custom width override in our dictionary
        if char in exceptional_symbols :
            width += exceptional_symbols [ char ]
            continue
        
        # Fallback to standard wcwidth / East Asian Width
        w = wcwidth.wcwidth ( char )
        if 0 <= w : width += w
        else:
            ea = unicodedata.east_asian_width ( char )
            if ea in ('W', 'F') : width += 2 
            else                : width += 1 
            
    return width


# =============================================================================
if '__main__' == __name__ :

    # =========================================================================
    from   ostap.logger.logger    import getLogger
    logger = getLogger ( 'ostap.logger.symbols' ) 
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    globs  = globals ().copy() 
    header = "Name" , "Value" , 'Width'
    import ostap.logger.table as T

    rows = [ header ] 
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
