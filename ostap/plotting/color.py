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
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.plotting.color' )
else                       : logger = getLogger( __name__ )
# =============================================================================
import ROOT 
# =============================================================================

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
    """ Get RGB fot the color
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


ROOT.TColor.rgb       = color_rgb
ROOT.TColor.rgba      = color_rgba
ROOT.TColor.hls       = color_hls

# =============================================================================
# gefine the color for latex <code>xcolor</code>
def color_rgb_latex ( c , name ) :
    """Define the color for latex <code>xcolor</code>
    """
    r , g , b = c.rgb() 
    return "\\definecolor{%s}[rgb]{%5.3f,%5.3f.%5.3f}" % ( name , r , g , b )

# =============================================================================
# gefine the color for latex <code>xcolor</code>
def color_hls_latex ( c , name ) :
    """Define the color for latex <code>xcolor</code>
    """
    h , l , s = c.hls() 
    return "\\definecolor{%s}[hls]{%d,%5.3f.%5.3f}" % ( name , int ( h )  , l , s )


ROOT.TColor.rgb_latex = color_rgb_latex
ROOT.TColor.hls_latex = color_hls_latex

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
