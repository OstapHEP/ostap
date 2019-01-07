#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file fit_draw.py
#  Default drawing options
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Default drawing options"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'lineWidth'             , ## line width
    'lineColor'             , ## line color 
    'lineStyle'             , ## line style 
    ##
    'keys'                  , ## predefines keys for draw-options
    'draw_options'          , ## pickup draw-options from the dictionary
    ##
    'data_options'          , ## draw options for data 
    'data_options_nobars'   , ## draw options for data without bars 
    'signal_options'        , ## draw options for "signal"     component(s) 
    'background_options'    , ## draw options for "background" component(s)
    'crossterm1_options'    , ## draw options for "crossterm1" component(s)
    'crossterm2_options'    , ## draw options for "crossterm2" component(s)    
    'component_options'     , ## draw options for "other"      component(s)
    'total_fit_options'     , ## draw options for the total fit curve
    ##
    'base_signal_color'     , ## base color for "signal"     component(s)
    'base_background_color' , ## base color for "background" component(s)
    'base_crossterm1_color' , ## base color for "crossterm1" component(s)
    'base_crossterm2_color' , ## base color for "crossterm2" component(s)
    'base_component_color'  , ## base color for "other"      component(s)
    ##
    'Color'                 , ## simple class to define color for the commponents 
    )
# =============================================================================
import ROOT
from   ostap.core.types import integer_types, list_types 
# =============================================================================
from   ostap.logger.logger  import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.plotting.fit_draw' )
else                       : logger = getLogger ( __name__         )
# =============================================================================
def lineWidth ( w ) : return ROOT.RooFit.LineWidth ( w )
def lineStyle ( s ) : return ROOT.RooFit.LineStyle ( s )
def lineColor ( c ) : return ROOT.RooFit.LineColor ( c )
#
## the list of predefined "draw"-keys 
keys = (
    'data_options'          ,
    ##
    'background_options'    ,
    'background_color'      ,
    ##
    'crossterm1_options'    ,
    'crossterm1_color'      ,
    ##
    'crossterm2_options'    ,
    'crossterm2_color'      ,
    ##
    'component_options'     ,
    'component_color'       ,
    ##
    'signal_options'        ,
    'signal_color'          ,
    ##
    'draw_axis_title'       , ## draw the titles for the axes ?
    'draw_options'          , ## generic ROOT draw options, e.g. 'same'
    )
# =============================================================================
## get draw options:
#  Collect predefined keys from  the dictioary
#  @code
#  def somefunc (  ... , **kwargs ) :
#      draw_opts = draw_options ( kwargs )
#  @endcode 
def draw_options ( **kwargs ) :
    """Collect predefined keys from  the dictioary
    >>> def somefunc (  ... , **kwargs ) :
    ...      draw_opts = draw_options ( kwargs )
    """
    options = {}
    for k,v in kwargs.iteritems() :
        if k.lower() in keys : options[ k.lower() ] = v
        if k.lower() in ( 'draw' , 'draw_option' , 'draw_options' ) :
            if isinstance ( v , dict ) : options.update ( v ) 
    return options

# =============================================================================
##  simpel class to represent the line colro for the component 
class Color(object) :
    def __init__ ( self                 ,
                   colors = None        ,
                   base   = ROOT.kBlue  ,
                   step   = 1           ) :

        self.__color_func = None
        self.__color_list = None
        
        if   isinstance ( base , integer_types   ) :
            bc = ROOT.gROOT.GetColor ( base )
            if not bc : logger.warning('Non-existing color %s is specified' % base ) 
        elif isinstance ( base , ROOT.TColor ) :
            bi = base.GetNumber()
            bc = ROOT.gROOT.GetColor( bi )
            if not bc : logger.warning('Non-existing color %s is specified' % bi )
            base = bi

        self.__base = base
        self.__step = step 
        
        if isinstance  ( colors , ROOT.TColor ) :
            bi = colors.GetNumber()
            bc = ROOT.gROOT.GetColor( bi )
            if not bc : logger.warning('Non-existing color %s is specified' % bi )
            colors = bi 
            
        if isinstance  ( colors , integer_types ) and 0 <= colors :
            bc = ROOT.gROOT.GetColor ( colors )
            if not bc : logger.warning('Non-existing color %s is specified' % colors )
            the_color = colors 
            colors = lambda i : the_color 
            
        if not colors :            
            colors = lambda i : base + i * step
            self.__color_func = colors            
        elif isinstance ( colors , Color ) :
            self.__base       = colors.__base
            self.__step       = colors.__step
            self.__color_func = colors.__color_func
            self.__color_list = colors.__color_list
        elif colors and isinstance ( colors , list_types ) :            
            self.__color_list = tuple ( [ i for i in colors ] ) 
        elif callable ( colors ) :
            
            try :
                i0 , i1 = colors ( 0 ) , colors ( 1 )
                assert isinstance ( i0 , integer_types ) , 'Invalid color(0) type'
                assert isinstance ( i1 , integer_types ) , 'Invalid color(1) type'                
                self.__color_func = colors            
            except :
                pass
                        
        assert self.__color_func or self.__color_list , 'Invalid color sequence!'
        
    ## get the color index for the certain fit component 
    def __call__ ( self , cmp ) :

        assert ( cmp , integer_types  ) and 0 <= cmp , 'Inavalid component index'
        if self.__color_func :
            color_func = self.__color_func
            return color_func ( cmp )

        colors = self.__color_list 
        nc = len ( colors )
        if cmp < nc :
            return colors [ cmp ] 
        else        :
            i = cmp - nc
            return colors[-1] + ( i + 1 ) * self.__step 
        
# =============================================================================
## plain, default
data_options           = ()

## suppress small bars at the end of error bars 
data_options_nobars     = ( ROOT.RooFit.MarkerStyle ( 20   ) ,
                            ROOT.RooFit.DrawOption  ( "zp" ) )

## signal:          thin dotted line
signal_options          = lineWidth ( 1 ) , lineStyle ( 1  )

## 1D background:   thin long-dashed line
background_options      = lineWidth ( 1 ) , lineStyle ( 7  )

## "component":     thin dash-dotted line
component_options       = lineWidth ( 1 ) , lineStyle ( ROOT.kDashDotted )

## total fit curve: thick red orange line 
total_fit_options       = lineWidth ( 3 ) , lineColor ( ROOT.kOrange + 1 ) , lineStyle ( 1 ) 

## background:  thin short-dashed line
background2D_options    = lineWidth ( 1 ) , lineStyle ( ROOT.kDashed     ) 
base_background2D_color = ROOT.kBlue 

crossterm1_options      = lineWidth ( 1 ) , lineStyle ( 7 )  
crossterm2_options      = lineWidth ( 1 ) , lineStyle ( 9 )  

signal_color      = Color ( ( ROOT.kRed          , ROOT.kRed     -  7 , ROOT.kRed     + 2 ,
                              ROOT.kRed     - 10 , ROOT.kRed     +  3 , ROOT.kRed     - 9 ,
                              ROOT.kRed     -  2 , ROOT.kRed     -  9 , ROOT.kRed     + 4 ) )
background_color  = Color ( ( ROOT.kBlue         , ROOT.kBlue    -  9 , ROOT.kBlue    + 3 ,
                              ROOT.kBlue    -  2 , ROOT.kBlue    - 10 , ROOT.kBlue    + 4 ) )
component_color   = Color ( ( ROOT.kGreen   +  3 , ROOT.kCyan         , ROOT.kMagenta     ,
                              ROOT.kOrange  -  6 , ROOT.kViolet  +  1 , ROOT.kYellow  + 1 ) )
crossterm1_color  = Color ( ( ROOT.kMagenta +  3 , ROOT.kMagenta - 10 , ROOT.kMagenta + 3 ,
                              ROOT.kMagenta -  6 , ROOT.kMagenta +  4 , ROOT.kMagenta - 5 ) )
crossterm2_color  = Color ( ( ROOT.kGreen   +  1 , ROOT.kGreen   -  1 , ROOT.kMagenta - 8 ,
                              ROOT.kGreen   +  2 , ROOT.kGreen   - 10 , ROOT.kGreen   - 4 ) )


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
