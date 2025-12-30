#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/plotting/color.py
#  Simple modification/extension of class TColor
#  @see TColor
#  @date   2012-02-15
#  @author Vanya BELYAEV Ivan.Belyaevitep.ru
# =============================================================================
""" Simple modification/extension of class TColor
"""
# =============================================================================
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = "2014-10-19"
__version__ = '$Revision$'
__all__     = ()
# =============================================================================
from ostap.core.meta_info import root_info
from ostap.utils.cidict   import cidict, cidict_fun 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.plotting.color' )
else                       : logger = getLogger( __name__ )
# =============================================================================
import ROOT 
# =============================================================================
## get RGB for the color
#  @code
#  color = ...
#  r, g, b = color.rgb() 
#  @endcode
#  @see TColor::GetRed
#  @see TColor::GetGreen
#  @see TColor::GetBlue
def color_rgb ( c ) :
    """ Get RGB for the color
    >>> color = ...
    >>> r, g, b = color.rgb() 
    - see ROOT.TColor.GetRed
    - see ROOT.TColor.GetGreen
    - see ROOT.TColor.GetBlue
    """
    return c.GetRed() , c.GetGreen() , c.GetBlue()

# =============================================================================
## get RGBA for the color
#  @code
#  color = ...
#  r, g, b, a  = color.rgba() 
#  @endcode
#  @see TColor::GetRed
#  @see TColor::GetGreen
#  @see TColor::GetBlue
#  @see TColor::GetAlpha
def color_rgba ( c ) :
    """ Get RGBA fot the color
    >>> color = ...
    >>> r, g, b, a = color.rgba() 
    - see ROOT.TColor.GetRed
    - see ROOT.TColor.GetGreen
    - see ROOT.TColor.GetBlue
    - see ROOT.TColor.GetAlpha
    """
    return c.GetRed() , c.GetGreen() , c.GetBlue(), c.GetAlpha() 

# =============================================================================
## get HLS for the color
#  @code
#  color = ...
#  h, l, s = color.hls() 
#  @endcode
#  @see TColor::GetHue
#  @see TColor::GetLight
#  @see TColor::GetSaturation
def color_hls ( c ) :
    """ Get HLS fot the color
    >>> color = ...
    >>> h, l, s = color.hls() 
    - see ROOT.TColor.GetHue
    - see ROOT.TColor.GetLight
    - see ROOT.TColor.GetSaturation
    """
    return c.GetHue() , c.GetLight() , c.GetSaturation()

# =============================================================================
## Get CMYK for the color 
#  @code
#  color = ...
#  c, m , y ,  k  = color.cmyk() 
#  @endcode
#  @see TColor::GetRed
#  @see TColor::GetGreen
#  @see TColor::GetBlue
#  @see https://www.rapidtables.com/convert/color/rgb-to-cmyk.html
def color_cmyk ( c ) :
    """ Get CMYK for the color
    >>> color = ...
    >>> h, l, s = color.hls() 
    - see ROOT.TColor.GetRed
    - see ROOT.TColor.GetGreen
    - see ROOT.TColor.GetBlue
    - see https://www.rapidtables.com/convert/color/rgb-to-cmyk.html
    """
    r , g , b = c.rgb ()

    r    = int ( r * 255 ) / 255
    g    = int ( g * 255 ) / 255
    b    = int ( b * 255 ) / 255
    
    cmax = max ( r , g , b )
    
    if not cmax : return 0 , 0 , 0 , 1   ## BLACK!
    
    k = 1 - cmax 
    
    c = ( 1 - r - k ) / ( 1 - k )
    m = ( 1 - g - k ) / ( 1 - k ) 
    y = ( 1 - b - k ) / ( 1 - k )

    return c, m, y, k 

ROOT.TColor.rgb  = color_rgb
ROOT.TColor.rgba = color_rgba
ROOT.TColor.hls  = color_hls
ROOT.TColor.cmyk = color_cmyk 

ROOT.TColor.RGB  = color_rgb
ROOT.TColor.RGBZ = color_rgba
ROOT.TColor.HLS  = color_hls
ROOT.TColor.CMYK = color_cmyk 

# =============================================================================
# gefine the color for latex <code>xcolor</code>
def color_rgb_latex ( c , name ) :
    """ Define the color for latex <code>xcolor</code>
    """
    r , g , b = c.rgb() 
    return "\\definecolor{%s}[rgb]{%5.3f,%5.3f.%5.3f}" % ( name , r , g , b )

# =============================================================================
# gefine the color for latex <code>xcolor</code>
def color_hls_latex ( c , name ) :
    """ Define the color for latex <code>xcolor</code>
    """
    h , l , s = c.hls() 
    return "\\definecolor{%s}[hls]{%d,%5.3f.%5.3f}" % ( name , int ( h )  , l , s )

ROOT.TColor.rgb_latex = color_rgb_latex
ROOT.TColor.hls_latex = color_hls_latex

# =============================================================================
# Create set of colors with given RGB
# @code
# color = create_color ( R , G , B , 0.1 , 'the_name' ) 
# @endcode 
def create_color ( R , G , B , alpha = 1 , name = '' ) :
    """ Create set of colors with given RGB
    >>> color = create_color ( R , G , B , 'my_color') 
    """
    if isinstance ( R , int ) and 0 <= R <= 255 and \
       isinstance ( G , int ) and 0 <= G <= 255 and \
       isinstance ( B , int ) and 0 <= B <= 255 : pass
    
    elif  0 <= R <= 1 and  0 <= G <= 1 and  0 <= B <= 1 :
        
        return create_color ( int ( R * 255 ) ,
                              int ( G * 255 ) ,
                              int ( B  *255 ) , alpha , name )
    
    else :
        
        raise TypeError ( 'ceate_color: Invalid `RGB`=(%s,%s,%s)' %  ( R , G , B ) )

    ## valid alpha? 
    assert 0 <= alpha <= 1 , "create_color: Invalid `alpha`=%s" % alpha
    
    ## color exist ? 
    if ( 6 , 32 ) <= root_info : color = ROOT.TColor.GetColor ( R , G , B , alpha )
    else                       : color = ROOT.TColor.GetColor ( R , G , B )
    
    if 0 <= color :
        groot = ROOT.GetROOT()
        if groot : 
            color = groot.GetColor ( color )
            if color and name and color.GetName() != name : color.SetName ( name )
            if color : return color 
    ## create new color 
    color = ROOT.TColor.GetFirstFreeColorIndex()
    if color < 0 : color = ROOT.TColor.GetFreeColorIndex()
    if not name :
        name = 'COLOR_%X%X%X'    % ( R , G , B )
        if 1 != alpha : name = '%s/%d%%' % ( name , int ( alpha * 100 ) )
    ## finally create the color 
    color = ROOT.TColor ( color , R , G , B , name , alpha ) 
    return color

# ==============================================================================
## print TColor 
def _color_str_ ( color ) :
    """ Print TColor Object
    """
    r, g , b , alpha = color.rgba ()
    n = color.GetNumber() 
    name        = color.GetName() 
    if not name : name = 'Color_%d' % n
    return '%s #%d RGB=%s' % ( name , n , color.AsHexString () )

ROOT.TColor.__str__  = _color_str_
ROOT.TColor.__repr__ = _color_str_

# =============================================================================
## reduce ROOT.TColor object  
def _color_reduce_ ( color ) :
    """ Reduce the `ROOT.TColor` object 
    """
    return create_color , ( int ( color.GetRed   () * 255 ) ,
                            int ( color.GetGreen () * 255 ) ,
                            int ( color.GetBlue  () * 255 ) ,
                            color.GetAlpha () ,
                            color.GetName  () )

ROOT.TColor.__reduce__ = _color_reduce_

# =============================================================================
# Create colors with given RGB and name 
# @code
# color = make_color ( R , G , B , 'the_name' ) 
# @endcode 
def make_color ( R , G , B , * , name = '' , alpha = 1 ) :
    """ Create set of colors with given RGB
    >>> color = make_color ( R , G , B , 'my_color') 
    """
    color = create_color ( R , G , B , alpha , name )
    return color.GetNumber()

# =============================================================================
## Web colors : https://en.wikipedia.org/wiki/Web_colors
# =============================================================================
## pink colors
# =============================================================================
MediumVioletRed	     = make_color ( 199,  21, 133 , name = "MediumVioletRed" )
DeepPink	     = make_color (  55,  20, 147 , name = "DeepPink"        )
PaleVioletRed        = make_color ( 219, 112, 147 , name = "PaleVioletRed"   )
HotPink	             = make_color ( 255, 105, 180 , name = "HotPink"         )
LightPink            = make_color ( 255, 182, 193 , name = "LightPink"       ) 
Pink	             = make_color ( 255, 192, 203 , name = "Pink"            )
# =============================================================================
## Purple, violet, and magenta colors
# =============================================================================
Indigo	             = make_color (  75,   0, 130 , name = "Indigo"          )
Purple	             = make_color ( 128,   0, 128 , name = "Purple"          ) 
DarkMagenta	     = make_color ( 139,   0, 139 , name = "DasrkMagenta"    ) 
DarkViolet	     = make_color ( 148,   0, 211 , name = "DarkViolet"      ) 
DarkSlateBlue	     = make_color (  72,  61, 139 , name = "DarkSlateBlue"   ) 
BlueViolet	     = make_color ( 138,  43, 226 , name = "BlueViolet"      ) 
DarkOrchid           = make_color ( 153,  50, 204 , name = "DarkOrchid"      ) 
Fuchsia	             = make_color ( 255,   0, 255 , name = "Fuchsia"         ) 
Magenta	             = make_color ( 255,   0, 255 , name = "Magenta"         ) 
SlateBlue	     = make_color ( 106,  90, 205 , name = "SlateBlue"       ) 
MediumSlateBlue	     = make_color ( 123, 104, 238 , name = "MediumSlateBlue" ) 
MediumOrchid	     = make_color ( 186,  85, 211 , name = "MediumOrchid"    )
MediumPurple         = make_color ( 147, 112, 219 , name = "MediumPurple"    ) 
Orchid	             = make_color ( 218, 112, 214 , name = "Orhchid"         ) 
Violet	             = make_color ( 238, 130, 238 , name = "Violet"          ) 
Plum	             = make_color ( 221, 160, 221 , name = "Plum"            ) 
Thistle	             = make_color ( 216, 191, 216 , name = "Thistle"         ) 
Lavender	     = make_color ( 230, 230, 250 , name = "Lavender"        ) 
# ===============================================================================
## Green colors
# ===============================================================================
DarkGreen	     = make_color (   0, 100,   0 , name = "DarkGreen"         )
Green	             = make_color (   0, 128,   0 , name = "Green"             ) 
DarkOliveGreen	     = make_color (  85, 107,  47 , name = "DarkOliveGreen"    ) 
ForestGreen	     = make_color (  34, 139,  34 , name = "ForestGreen"       ) 
SeaGreen	     = make_color (  46, 139,  87 , name = "SeaGreen"          ) 
Olive	             = make_color ( 128, 128,   0 , name = "Olive"             ) 
OliveDrab	     = make_color ( 107, 142,  35 , name = "OliveDrab"         ) 
MediumSeaGreen	     = make_color (  60, 179, 113 , name = "MediumSeaGreen"    )
LimeGreen	     = make_color (  50, 205,  50 , name = "LimeGreen"         ) 
Lime	             = make_color (   0, 255,   0 , name = "Lime"              ) 
SpringGreen	     = make_color (   0, 255, 127 , name = "SpringGreen"       ) 
MediumSpringGreen    = make_color (   0, 250, 154 , name = "MediumSpringGreen" ) 
DarkSeaGreen	     = make_color ( 143, 188, 143 , name = "DarkSeeGreen"      )  
MediumAquamarine     = make_color ( 102, 205, 170 , name = "MediumAquamarine"  ) 
YellowGreen	     = make_color ( 154, 205,  50 , name = "YellowGreen"       ) 
LawnGreen	     = make_color ( 124, 252,   0 , name = "LawnGreen"         )
Chartreuse	     = make_color ( 127, 255,   0 , name = "Chartreuse"        ) 
LightGreen	     = make_color ( 144, 238, 144 , name = "LightGreen"        ) 
GreenYellow	     = make_color ( 173, 255,  47 , name = "GreenYellow"       ) 
PaleGreen	     = make_color ( 152, 251, 152 , name = "PaleGreen"         )
# ===============================================================================
## Red colors
# ===============================================================================
DarkRed	             = make_color ( 139,   0,   0 , name = "DarkRed"           )
Red	             = make_color ( 255,   0,   0 , name = "Red"               ) 
Firebrick	     = make_color ( 178,  34,  34 , name = "Firebrick"         ) 
Crimson	             = make_color ( 220,  20,  60 , name = "Crimson"           )
IndianRed	     = make_color ( 205,  92,  92 , name = "IndianRed"         ) 
LightCoral	     = make_color ( 240, 128, 128 , name = "LightCoral"        )
Salmon	             = make_color ( 250, 128, 114 , name = "Salmon"            ) 
DarkSalmon	     = make_color ( 233, 150, 122 , name = "DarkSalmon"        )
LightSalmon	     = make_color ( 255, 160, 122 , name = "LightSalmon"       )
# ============================================================================
## Orange colors
# ============================================================================
OrangeRed	     = make_color ( 255,  69,   0 , name = "OrangeRed"         ) 
Tomato	             = make_color ( 255,  99,  71 , name = "Tomato"            ) 
DarkOrange	     = make_color ( 255, 140,   0 , name = "DarkOrange"        )
Coral	             = make_color ( 255, 127,  80 , name = "Coral"             ) 
Orange	             = make_color ( 255, 165,   0 , name = "Orange"            ) 
# ===========================================================================
## Yellow colors
# ===========================================================================
DarkKhaki	     = make_color ( 189, 183, 107 , name = "DarkKhaki"            ) 
Gold	             = make_color ( 255, 215,   0 , name = "Gold"                 ) 
Khaki	             = make_color ( 240, 230, 140 , name = "Khaki"                ) 
PeachPuff	     = make_color ( 255, 218, 185 , name = "PeachPuff"            ) 
Yellow	             = make_color ( 255, 255, 0   , name = "Yellow"               ) 
PaleGoldenrod	     = make_color ( 238, 232, 170 , name = "PaleGoldenrod"        ) 
Moccasin	     = make_color ( 255, 228, 181 , name = "Mocassin"             ) 
PapayaWhip	     = make_color ( 255, 239, 213 , name = "PapayaWhip"           ) 
LightGoldenrodYellow = make_color ( 250, 250, 210 , name = "LightGoldenrodYellow" ) 
LemonChiffon         = make_color ( 255, 250, 205 , name = "LemonChiffon"         ) 
LightYellow	     = make_color ( 255, 255, 224 , name = "LightYellow"          )
# ===========================================================================
## Blue colors
# ===========================================================================
MidnightBlue	     = make_color (  25,  25, 112 , name = "MignightBlue"         ) 
Navy	             = make_color (   0,   0, 128 , name = "Navy"                 )
DarkBlue	     = make_color (   0,   0, 139 , name = "DarkBlue"             )  
MediumBlue	     = make_color (   0,   0, 205 , name = "MediumBlue"           ) 
Blue	             = make_color (   0,   0, 255 , name = "Blue"                 ) 
RoyalBlue	     = make_color (  65, 105, 225 , name = "RoyalBlue"            ) 
SteelBlue	     = make_color (  70, 130, 180 , name = "SteelBlue"            ) 
DodgerBlue	     = make_color (  30, 144, 255 , name = "DodgenBlue"           ) 
DeepSkyBlue	     = make_color (   0, 191, 255 , name = "DeepSkyBlue"          ) 
CornflowerBlue	     = make_color ( 100, 149, 237 , name = "CornflowerBlue"       ) 
SkyBlue	             = make_color ( 135, 206, 235 , name = "SkyBlue"              ) 
LightSkyBlue	     = make_color ( 135, 206, 250 , name = "LightSkyBlue"         ) 
LightSteelBlue	     = make_color ( 176, 196, 222 , name = "LightSteelBlue"       ) 
LightBlue	     = make_color ( 173, 216, 230 , name = "LightBlue"            ) 
PowderBlue	     = make_color ( 176, 224, 230 , name = "PowdrBlue"            )
# ============================================================================
## White colors
# ============================================================================
MistyRose	     = make_color ( 255, 228, 225 , name = "MistyRose"     ) 
AntiqueWhite	     = make_color ( 250, 235, 215 , name = "AntiqueWhite"  ) 
Linen	             = make_color ( 250, 240, 230 , name = "Linen"         ) 
Beige	             = make_color ( 245, 245, 220 , name = "Beige"         ) 
WhiteSmoke	     = make_color ( 245, 245, 245 , name = "WhiteSmoke"    ) 
LavenderBlush	     = make_color ( 255, 240, 245 , name = "LavenderBlush" ) 
OldLace	             = make_color ( 253, 245, 230 , name = "OldLace"       ) 
AliceBlue	     = make_color ( 240, 248, 255 , name = "AliceBlue"     ) 
Seashell	     = make_color ( 255, 245, 238 , name = "Seashell"      )
GhostWhite	     = make_color ( 248, 248, 255 , name = "GhostWhite"    ) 
Honeydew	     = make_color ( 240, 255, 240 , name = "Honewdew"      ) 
FloralWhite	     = make_color ( 255, 250, 240 , name = "FloralWhite"   ) 
Azure	             = make_color ( 240, 255, 255 , name = "Azure"         ) 
MintCream	     = make_color ( 245, 255, 250 , name = "MintCream"     ) 
Snow	             = make_color ( 255, 250, 250 , name = "Snow"          ) 
Ivory	             = make_color ( 255, 255, 240 , name = "Ivory"         ) 
White	             = make_color ( 255, 255, 255 , name = "White"         )
# =============================================================================
## Brown colors
# =============================================================================
Maroon	             = make_color ( 128,   0,   0 , name = "Maroon"         )  
Brown	             = make_color ( 165,  42,  42 , name = "Brown"          )
SaddleBrown	     = make_color ( 139,  69,  19 , name = "SaddleBrown"    )
Sienna	             = make_color ( 160,  82,  45 , name = "Sienna"         )  
Chocolate	     = make_color ( 210, 105,  30 , name = "Chocolate"      ) 
DarkGoldenrod	     = make_color ( 184, 134,  11 , name = "DarkGoldenrod"  ) 
Peru	             = make_color ( 205, 133,  63 , name = "Peru"           )
RosyBrown	     = make_color ( 188, 143, 143 , name = "RosyBrown"      ) 
Goldenrod	     = make_color ( 218, 165,  32 , name = "GoldenRod"      ) 
SandyBrown	     = make_color ( 244, 164,  96 , name = "SandyBrown"     )
Tan	             = make_color ( 210, 180, 140 , name = "Tan"            )
Burlywood	     = make_color ( 222, 184, 135 , name = "Burlywood"      )      
Wheat	             = make_color ( 245, 222, 179 , name = "Wheat"          ) 
NavajoWhite	     = make_color ( 255, 222, 173 , name = "NavajoWhite"    ) 
Bisque	             = make_color ( 255, 228, 196 , name = "Binque"         )
BlanchedAlmond	     = make_color ( 255, 235, 205 , name = "BlanchedAlmond" )
Cornsilk	     = make_color ( 255, 248, 220 , name = "Cornsilk"       )
# =============================================================================
## Cyan colors
# =============================================================================
Teal	             = make_color (   0, 128, 128 , name = "Teal"            )  
DarkCyan	     = make_color (   0, 139, 139 , name = "DarkCyan"        )  
LightSeaGreen	     = make_color (  32, 178, 170 , name = "LightSeaGreen"   )  
CadetBlue	     = make_color (  95, 158, 160 , name = "CaderBlue"       )  
DarkTurquoise	     = make_color (   0, 206, 209 , name = "DarkTurquoise"   ) 
MediumTurquoise	     = make_color (  72, 209, 204 , name = "MediumRurquoise" ) 
Turquoise	     = make_color (  64, 224, 208 , name = "Torquise"        )
Aqua	             = make_color (   0, 255, 255 , name = "Aqua"            )
Cyan	             = make_color (   0, 255, 255 , name = "Cyan"            )
Aquamarine	     = make_color ( 127, 255, 212 , name = "Aquamarine"      )
PaleTurquoise	     = make_color ( 175, 238, 238 , name = "PaleTurquoise"   ) 
LightCyan	     = make_color ( 224, 255, 255 , name = "LightCyan"       ) 
# =============================================================================
## Gray and black colors
# =============================================================================
Black	             = make_color (   0,   0,   0 , name = "Black"         ) 
DarkSlateGray	     = make_color (  47,  79,  79 , name = "DarkSlateGrey" ) 
DimGray	             = make_color ( 105, 105, 105 , name = "DimGrey"       ) 
SlateGray	     = make_color ( 112, 128, 144 , name = "SlateGrey"     ) 
Gray	             = make_color ( 128, 128, 128 , name = "Gray"          ) 
LightSlateGray	     = make_color ( 119, 136, 153 , name = "LightGray"     ) 
DarkGray	     = make_color ( 169, 169, 169 , name = "DarkGrey"      ) 
Silver	             = make_color ( 192, 192, 192 , name = "Silver"        ) 
LightGray	     = make_color ( 211, 211, 211 , name = "LightGray"     ) 
Gainsboro	     = make_color ( 220, 220, 220 , name = "Gainsboro"     ) 

# =============================================================================
##  color collections
colors_pink    = ( MediumVioletRed      , 
                   DeepPink	        , 
                   PaleVioletRed        , 
                   HotPink	        ,
                   LightPink            , 
                   Pink	                )
colors_purple  = ( Indigo	        , 
                   Purple	        ,
                   DarkMagenta	        , 
                   DarkViolet	        , 
                   DarkSlateBlue        , 
                   BlueViolet	        , 
                   DarkOrchid           , 
                   Fuchsia	        , 
                   Magenta	        , 
                   SlateBlue	        , 
                   MediumSlateBlue      , 
                   MediumOrchid	        , 
                   MediumPurple         , 
                   Orchid	        , 
                   Violet	        , 
                   Plum	                , 
                   Thistle	        , 
                   Lavender	        )
colors_violet  = colors_purple
colors_magenta = colors_purple
colors_green   = ( DarkGreen	        , 
                   Green	        , 
                   DarkOliveGreen       ,
                   ForestGreen	        ,
                   SeaGreen	        ,
                   Olive	        ,
                   OliveDrab	        , 
                   MediumSeaGreen       , 
                   LimeGreen	        , 
                   Lime	                ,
                   SpringGreen	        , 
                   MediumSpringGreen    ,
                   DarkSeaGreen	        , 
                   MediumAquamarine     ,
                   YellowGreen	        ,
                   LawnGreen	        ,
                   Chartreuse	        ,
                   LightGreen	        ,
                   GreenYellow	        ,
                   PaleGreen	        )
colors_red     = ( DarkRed	        , 
                   Red	                , 
                   Firebrick	        , 
                   Crimson	        , 
                   IndianRed	        , 
                   LightCoral	        , 
                   Salmon	        , 
                   DarkSalmon	        , 
                   LightSalmon	        )
colors_orange  = ( OrangeRed	        , 
                   Tomato	        , 
                   DarkOrange	        , 
                   Coral	        , 
                   Orange	        )
colors_yellow  = ( DarkKhaki	        , 
                   Gold	                , 
                   Khaki	        , 
                   PeachPuff	        , 
                   Yellow	        , 
                   PaleGoldenrod        , 
                   Moccasin	        , 
                   PapayaWhip	        , 
                   LightGoldenrodYellow , 
                   LemonChiffon         , 
                   LightYellow	        )
colors_blue    = ( MidnightBlue	        , 
                   Navy	                , 
                   DarkBlue	        , 
                   MediumBlue	        , 
                   Blue	                , 
                   RoyalBlue	        ,
                   SteelBlue	        , 
                   DodgerBlue	        , 
                   DeepSkyBlue	        , 
                   CornflowerBlue	, 
                   SkyBlue	        , 
                   LightSkyBlue	        , 
                   LightSteelBlue	, 
                   LightBlue	        , 
                   PowderBlue	        )
colors_white  = ( MistyRose	        , 
                  AntiqueWhite	        , 
                  Linen	                , 
                  Beige	                , 
                  WhiteSmoke	        , 
                  LavenderBlush	        , 
                  OldLace	        , 
                  AliceBlue	        , 
                  Seashell	        , 
                  GhostWhite	        , 
                  Honeydew	        , 
                  FloralWhite	        , 
                  Azure	                , 
                  MintCream	        , 
                  Snow	                , 
                  Ivory	                ,  
                  White	                )
colors_brown  = ( Maroon	        , 
                  Brown	                , 
                  SaddleBrown	        , 
                  Sienna	        , 
                  Chocolate	        , 
                  DarkGoldenrod	        , 
                  Peru	                , 
                  RosyBrown	        , 
                  Goldenrod	        , 
                  SandyBrown	        , 
                  Tan	                , 
                  Burlywood	        , 
                  Wheat	                , 
                  NavajoWhite	        , 
                  Bisque	        , 
                  BlanchedAlmond	,
                  Cornsilk	        )
colors_cyan   = ( Teal	                , 
                  DarkCyan	        , 
                  LightSeaGreen	        , 
                  CadetBlue	        , 
                  DarkTurquoise	        , 
                  MediumTurquoise	, 
                  Turquoise	        , 
                  Aqua	                , 
                  Cyan	                , 
                  Aquamarine	        , 
                  PaleTurquoise	        , 
                  LightCyan	        )
colors_gray   = ( Black	                , 
                  DarkSlateGray	        , 
                  DimGray	        , 
                  SlateGray	        , 
                  Gray	                , 
                  LightSlateGray	, 
                  DarkGray	        , 
                  Silver	        , 
                  LightGray	        , 
                  Gainsboro	        )
colors_black = colors_gray
# =============================================================================
## dictionatry with color lines 
# =============================================================================
colors = cidict ( transform = cidict_fun )
colors [ 'pink'    ] = colors_pink
colors [ 'red'     ] = colors_red 
colors [ 'purple'  ] = colors_purple
colors [ 'violet'  ] = colors_violet
colors [ 'magenta' ] = colors_magenta
colors [ 'green'   ] = colors_green 
colors [ 'orange'  ] = colors_orange
colors [ 'yellow'  ] = colors_yellow 
colors [ 'blue'    ] = colors_blue 
colors [ 'white'   ] = colors_white
colors [ 'cyan'    ] = colors_cyan
colors [ 'gray'    ] = colors_gray
colors [ 'black'   ] = colors_black
keys = sorted ( colors.keys () )
for key in keys :
    key_dark   = '%s_dark'   % key
    key_bright = '%s_bright' % key
    the_colors = colors [ key ]  
    colors [ key_dark   ] = tuple ( ROOT.TColor.GetColorDark   ( c ) for c in the_colors )
    colors [ key_bright ] = tuple ( ROOT.TColor.GetColorBright ( c ) for c in the_colors )
    
# =============================================================================

_decorated_classes_  = (
    ROOT.TColor  , 
    )

_new_methods_        = (
    ROOT.TColor.rgb       , 
    ROOT.TColor.rgba      , 
    ROOT.TColor.hls       , 
    ROOT.TColor.rgb_latex , 
    ROOT.TColor.hls_latex , 
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
