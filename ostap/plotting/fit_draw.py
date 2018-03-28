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
    )
# =============================================================================
import ROOT
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
    'base_background_color' ,
    ##
    'crossterm1_options'    ,
    'base_crossterm1_color' ,
    ##
    'crossterm2_options'    ,
    'base_crossterm2_color' ,
    ##
    'component_options'     ,
    'base_component_color'  ,
    ##
    'signal_options'        ,
    'base_signal_color'     ,
    ##
    'background_step_color' , 
    'crossterm1_step_color' ,
    'crossterm2_step_color' ,
    'component_step_color'  ,
    'signal_step_color'     ,
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

base_signal_color       = ROOT.kRed
base_background_color   = ROOT.kBlue
base_component_color    = ROOT.kMagenta

## background:  thin short-dashed line
background2D_options    = lineWidth ( 1 ) , lineStyle ( ROOT.kDashed     ) 
base_background2D_color = ROOT.kBlue 

crossterm1_options      = lineWidth ( 1 ) , lineStyle ( 7 )  
crossterm2_options      = lineWidth ( 1 ) , lineStyle ( 9 )  

base_crossterm1_color   = ROOT.kMagenta
base_crossterm2_color   = ROOT.kGreen+1 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
