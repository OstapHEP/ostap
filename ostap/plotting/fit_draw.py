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
    'signal_style'          , ## style for "signal"     component(s)
    'background_style'      , ## style for "background" component(s)
    'crossterm1_style'      , ## style for "crossterm1" component(s)
    'crossterm2_style'      , ## style for "crossterm2" component(s)
    'component_style'       , ## style for "other"      component(s)
    ##
    'Color'                 , ## simple class to define color for the commponents 
    )
# =============================================================================
import ROOT
from   ostap.core.types import integer_types, list_types 
import ostap.plotting.style  
import ostap.plotting.canvas
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
    'background_style'      ,
    ##
    'crossterm1_options'    ,
    'crossterm1_style'      ,
    ##
    'crossterm2_options'    ,
    'crossterm2_style'      ,
    ##
    'component_options'     ,
    'component_style'       ,
    ##
    'signal_options'        ,
    'signal_style'          ,
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
signal_options          = lineWidth ( 1 ) ,

## 1D background:   thin long-dashed line
background_options      = lineWidth ( 1 ) ,

## "component":     thin dash-dotted line
component_options       = lineWidth ( 1 ) ,

## total fit curve: thick red orange line 
total_fit_options       = lineWidth ( 3 ) , lineColor ( ROOT.kOrange + 1 ) , 

## background:  thin short-dashed line
background2D_options    = lineWidth ( 1 ) , 
crossterm1_options      = lineWidth ( 1 ) , 
crossterm2_options      = lineWidth ( 1 ) , 


# ===========================================================================================
## @class Style
#  Store the drawinng style for the component 
class Style(object):

    def __init__ ( self                     ,
                   linecolor = ROOT.kBlack  ,
                   linestyle = 1            , 
                   fillcolor = None         ,
                   fillstyle = None         , *args ) :

        options = []
        
        if isinstance ( linecolor , ROOT.TColor ) : linecolor = linecolor.GetNumber()
        if isinstance ( fillcolor , ROOT.TColor ) : fillcolor = fillcolor.GetNumber()

        for a in args :
            if isinstance  ( a , ROOT.RooCmdArg ) :
                if   'LineColor' == a.GetName () : linecolor = a
                elif 'LineStyle' == a.GetName () : linestyle = a
                elif 'FillColor' == a.GetName () : fillcolor = a
                elif 'FillStyle' == a.GetName () : fillstyle = a

        if isinstance ( linecolor , integer_types ) and 0 < linecolor :
            options.append ( ROOT.RooFit.LineColor ( linecolor ) )
        elif isinstance ( linecolor , ROOT.RooCmdArg ) and 'LineColor' == linecolor.GetName () : 
            options.append ( linecolor )
            
        if isinstance ( linestyle , integer_types ) and 0 < linestyle :
            options.append ( ROOT.RooFit.LineStyle ( linestyle ) )
        elif isinstance ( linestyle , ROOT.RooCmdArg ) and 'LineStyle' == linestyle.GetName () : 
            options.append ( linestyle  )

        _fillopt = False 
        if isinstance ( fillcolor , integer_types ) and 0 < fillcolor :
            fc = ROOT.TColor.GetColorTransparent ( fillcolor , 0.25 ) ## transparency!! 
            options.append ( ROOT.RooFit.FillColor ( fc ) )
            _fillopt = True 
        elif isinstance ( fillcolor , ROOT.RooCmdArg ) and 'FillColor' == fillcolor.GetName () : 
            _fillopt = True 
            options.append ( fillcolor  )
            
        if isinstance ( fillstyle , integer_types ) and 0 < fillstyle :
            options.append ( ROOT.RooFit.FillStyle ( fillstyle ) )
            _fillopt = True 
        elif isinstance ( fillstyle , ROOT.RooCmdArg )  and 'FillStyle' == fillstyle.GetName () : 
            options.append ( fillstyle  )
            _fillopt = True 

        if _fillopt :
            options.append ( ROOT.RooFit.VLines     (      ) )
            options.append ( ROOT.RooFit.DrawOption ( "FL" ) )

        self.__options = tuple ( options )
        
    @property 
    def options ( self ) :
        """``options'' : get the constructed list of options"""
        return self.__options

    
class Styles(object) :

    def __init__ ( self            ,
                   styles     = [] ,
                   linecolors = [ i for i in range ( 1 , 51 ) ] ,
                   linestyles = [ i for i in range ( 1 , 16 ) ] ,
                   fillcolors = [ i for i in range ( 1 , 51 ) ] ,
                   fillstyles = [] ) :

        _styles = []
        
        for s in styles :
            ss = None 
            if   isinstance ( s , Style          ) : ss = s 
            elif isinstance ( s , dict           ) : ss = Style ( **s )
            elif isinstance ( s , list_types     ) : ss = Style (  *s )
            elif isinstance ( s , integer_types  ) : ss = Style ( linecolor = s  )
            elif isinstance ( s , ROOT.TColor    ) : ss = Style ( linecolor = s  )
            elif isinstance ( s , ROOT.RooCmdArg ) : ss = Style ( None , None , None  , None  , s  ) 
            if ss : _styles.append ( ss ) 

        self.__styles     = tuple ( _styles )
        
        self.__linecolors = linecolors
        self.__linestyles = linestyles
        self.__fillcolors = fillcolors
        self.__fillstyles = fillstyles

        
    def __get_it ( self , index , objects , defval = None ) :

        if isinstance ( objects , list_types ) and objects :
            l = len  ( objects ) 
            if index < l : return objects [ index     ] 
            elif objects : return objects [ index % l ] 
        
        if   isinstance ( objects , integer_types   ) : return objects
        elif isinstance ( objects , ROOT.RooCmdArg  ) : return objects
        elif callable   ( objects                   ) : return objects ( index ) 
        return defval 
    
    @property
    def styles ( self ) :
        return self.__styles

    def linecolor ( self , index ) :
        s = self.__get_it ( index  , self.__linecolors  )
        if   isinstance  ( s , integer_types  ) : return ROOT.RooFit.LineColor  ( s )
        elif isinstance  ( s , ROOT.RooCmdArg ) : return s
        return None

    def linestyle ( self , index ) :
        s = self.__get_it ( index  , self.__linestyles  )
        if   isinstance  ( s , integer_types  ) : return ROOT.RooFit.LineStyle  ( s )
        elif isinstance  ( s , ROOT.RooCmdArg ) : return s
        return None

    def fillcolor ( self , index ) :
        s = self.__get_it ( index  , self.__fillcolors  )
        if   isinstance  ( s , integer_types  ) : return ROOT.RooFit.FillColor  ( s )
        elif isinstance  ( s , ROOT.RooCmdArg ) : return s
        return None

    def fillstyle ( self , index ) :
        s = self.__get_it ( index  , self.__fillstyles  )
        if   isinstance  ( s , integer_types  ) : return ROOT.RooFit.FillStyle  ( s )
        elif isinstance  ( s , ROOT.RooCmdArg ) : return s
        return None
    
    def __call__ ( self , index ) :

        if self.styles :
            l = len ( self.styles ) 
            return self.styles[ index % l ].options

        l = len ( self.styles )
        lc = self.linecolor ( index - l ) 
        ls = self.linestyle ( index - l ) 
        fc = self.fillcolor ( index - l ) 
        fs = self.fillstyle ( index - l ) 

        ss = Style ( lc  , ls , fc , fs )
        
        return ss.options
    
signal_style = Styles ( (
    Style( ROOT.kRed      , 1 , ROOT.kRed      , 1001 ) ,
    Style( ROOT.kRed -  7 , 1 , ROOT.kRed -  7 , 3395 ) ,
    Style( ROOT.kRed +  2 , 1 , ROOT.kRed +  2 , 3305 ) ,
    Style( ROOT.kRed - 10 , 1 , ROOT.kRed - 10 , 3345 ) ,    
    Style( ROOT.kRed +  3 , 1 , ROOT.kRed +  3 , 3354 ) ,    
) ) 
background_style = Styles ( (
    Style( ROOT.kBlue         ,   7 ) ,
    Style( ROOT.kBlue    -  9 ,  11 ) ,
    Style( ROOT.kBlue    +  3 ,  12 ) ,
    Style( ROOT.kBlue    -  2 ,  13 ) ,
    Style( ROOT.kBlue    - 10 ,  14 ) ,
) )
component_style = Styles ( (
    Style( 8 , 1 , 8 , 3354 ) ,
    Style( 6 , 1 , 6 , 3345 ) ,
    Style( 7 , 1 , 7 , 3395 ) ,
    Style( 9 , 1 , 9 , 3305 ) ,
    Style( ROOT.kPink   , 1 , ROOT.kPink   , 3315 ) ,
    Style( ROOT.kOrange , 1 , ROOT.kOrange , 3351 ) ,
    Style( 5 , 1 , 5 , 1001 ) ,
) ) 
crossterm1_style = Styles ( (
    Style( ROOT.kMagenta +  3 ,  11 , None , 3345 ) ,
    Style( ROOT.kMagenta - 10 ,  12 , None , 3354 ) ,
    Style( ROOT.kMagenta +  3 ,  13 , None , 3395 ) ,
    Style( ROOT.kMagenta -  3 ,  14 , None , 3305 ) ,
) ) 
crossterm2_style = Styles ( (
    Style( ROOT.kGreen   +  1 ,  14 ) ,
    Style( ROOT.kGreen   -  1 ,  13 ) ,
    Style( ROOT.kGreen   - 10 ,  12 ) ,
    Style( ROOT.kGreen   +  1 ,  11 ) ,
) )


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
