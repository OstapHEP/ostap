#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/plotting/fit_draw.py
#  The default drawing options for Ostap/RooFit.
#
#  The module defines the drawing options for Ostap/RooFit.
#  The options are split into two categories:
#   - <code>options</code>: RooFit optioms that are the same for the given component
#             class, e.g. options that are common for all "signal" or "backround: components 
#   - <code>styles</code>: RooFit options that are the same for the given component
#             class, e.g. options that are common for all "signal" or "backround: components 
#
#  The major options and  styles  are :
#   - data_options'                 : data options to draw the data points 
#   - signal_options'               : draw options for "signal"       component(s)
#   - background_options'           : draw options for "background"   component(s)
#   - background2D_options'         : draw options for "background2D" component(s)
#   - crossterm1_options'           : draw options for "crossterm1"   component(s)
#   - crossterm2_options'           : draw options for "crossterm2"   component(s)    
#   - component_options'            : draw options for "other"        component(s)
#   - total_fit_options'            : draw options for the total fit curve
#
#  The drawing can be done via the explicit usage of the options
#  @code
#  pdf.fitTo ( ... )
#  pdf.draw  ( ... , signal_options = ... , signal_style = ... ) 
#  @endcode
#
#  Also one can specify the drawing options for PDF :
#  @code
#  pdf.draw_options['data_options'] = ...
#  pdf.draw_options['signal_style'] = ...
#  pdf.fitTo ( ... )
#  pdf.draw  ( ... ) 
#  @endcode
#
#  The options can be specified also via the <code>[Fit Draw]</code> section in
#  the configuration file: 
#  @code
#  ``
#  [Fit Draw]
#  Data_options = ROOT.RooFit.MarkerStyle ( 20   ) ,
#                 ROOT.RooFit.DrawOption  ( "zp" ) 
#  signal_style = Style ( linecolor = ROOT.kRed ,
#                         linewidth = 2         ,
#                         fillcolor = ROOT.kRed ,
#                         fillstyle = 1001      )
#  default_background_style = Line ( linecolor = ROOT.kBlue      , linestyle =  7 ) ,
#                             Line ( linecolor = ROOT.kBlue -  9 , linestyle = 11 ) ,
#                             Line ( linecolor = ROOT.kBlue +  3 , linestyle = 12 ) ,
#                             Line ( linecolor = ROOT.kBlue -  2 , linestyle = 13 ) ,
#                             Line ( linecolor = ROOT.kBlue - 10 , linestyle = 14 ) ,
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Default drawing options for Ostap/RooFit 

The module defines the drawing options for Ostap/RooFit.
The options are split into two categories:
- options : RooFit optioms that are the same for the given component
  class, e.g. options that are common for all 'signal' or 'backround' components 
- styles  : RooFit options that are the same for the given component
class, e.g. options that are common for all 'signal' or 'backround' components 

The major options and  styles  are :
- `data_options'                 : data options to draw the data points 
- `signal_options'               : draw options for `signal '      component(s)
- `background_options'           : draw options for `background'   component(s)
- `background2D_options'         : draw options for `background2D' component(s)
- `crossterm1_options'           : draw options for `crossterm1'   component(s)
- `crossterm2_options'           : draw options for `crossterm2'   component(s)    
- `component_options'            : draw options for `other'        component(s)
- `total_fit_options'            : draw options for the total fit curve

The drawing can be done via the explicit usage of the options

>>> pdf.fitTo ( ... )
>>> pdf.draw  ( ... , signal_options = ... , signal_style = ... ) 

Also one can specify the drawing options for PDF :

>>> pdf.draw_options['data_options'] = ...
>>> pdf.draw_options['signal_style'] = ...
>>> pdf.fitTo ( ... )
>>> pdf.draw  ( ... ) 

The options can be specified also via the [Fit Draw] section in the configuration file: 

``
[Fit Draw]
Data_options = ROOT.RooFit.MarkerStyle ( 20   ) ,
               ROOT.RooFit.DrawOption  ( 'zp' ) 
Signal_style = Style ( linecolor = ROOT.kRed ,
                       linewidth = 2         ,
                       fillcolor = ROOT.kRed ,
                       fillstyle = 1001      )
''
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'lineWidth'                    , ## line width
    'lineColor'                    , ## line color 
    'lineStyle'                    , ## line style 
    ##
    'keys'                         , ## predefines keys for draw-options
    'draw_options'                 , ## pickup draw-options from the dictionary
    ##
    'data_options_plain'           , ## draw options for data 
    'data_options_nobars'          , ## draw options for data without bars
    'data_options_small'           , ## draw options for data without bars, small markers
    ##
    'data_options'                 , ## data options to draw the data points 
    'signal_options'               , ## draw options for "signal"       component(s)
    'background_options'           , ## draw options for "background"   component(s)
    'background2D_options'         , ## draw options for "background2D" component(s)
    'crossterm1_options'           , ## draw options for "crossterm1"   component(s)
    'crossterm2_options'           , ## draw options for "crossterm2"   component(s)    
    'component_options'            , ## draw options for "other"        component(s)
    'total_fit_options'            , ## draw options for the total fit curve
    'curve_options'                , ## draw options for the curve 
    ##
    'signal_style'                 , ## style for "signal"        component(s)
    'background_style'             , ## style for "background"    component(s)
    'background2D_style'           , ## style for "background-2D" component(s)
    'crossterm1_style'             , ## style for "crossterm1"    component(s)
    'crossterm2_style'             , ## style for "crossterm2"    component(s)
    'component_style'              , ## style for "other"         component(s)
    ##
    'combined_signal_options'      , ## draw options for combined "signal"       component(s)
    'combined_background_options'  , ## draw options for combined "background"   component(s)
    'combined_component_options'   , ## draw options for "other"        component(s)
    ##
    'combined_signal_style'        , ## style for combined "signal"        component(s)
    'combined_background_style'    , ## style for combined "background"    component(s)
    'combined_component_style'     , ## style for combined "other"         component(s)
    ##
    'default_data_options'         , ## defautl data options
    'default_signal_options'       , ## default options for the signal component
    'default_background_options'   , ## draw options for "background" component(s)
    'default_background2D_options' , ## draw options for "background" component(s)
    'default_crossterm1_options'   , ## draw options for "crossterm1" component(s)
    'default_crossterm2_options'   , ## draw options for "crossterm2" component(s)
    'default_component_options'    , ## draw options for "other"      component(s)
    'default_total_fit_options'    , ## draw options for the total fit curve    
    'default_curve_options'        , ## draw options for the curve    
    ##
    'default_signal_style'         , ## style for "signal"        component(s)
    'default_background_style'     , ## style for "background"    component(s)
    'default_background2D_style'   , ## style for "background-2D" component(s)
    'default_crossterm1_style'     , ## style for "crossterm1"    component(s)
    'default_crossterm2_style'     , ## style for "crossterm2"    component(s)
    'default_component_style'      , ## style for "other"         component(s)
    ##
    'Style'                        , ## helper class to define the style for the component
    'Line'                         , ## helper class to define the style for the component
    'Area'                         , ## helper class to define the style for the component
    'Styles'                       , ## helper class to define the style for the component
    )
# =============================================================================
import ROOT
from   ostap.core.ostap_types import integer_types, list_types
from   ostap.core.core  import items_loop
import ostap.plotting.style  
import ostap.plotting.canvas
import ostap.fitting.roocmdarg   
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
    'data_options'                ,
    ##
    'background_options'          ,
    'background_style'            ,
    ##
    'background2D_options'        ,
    'background2D_style'          ,
    ##
    'crossterm1_options'          ,
    'crossterm1_style'            ,
    ##
    'crossterm2_options'          ,
    'crossterm2_style'            ,
    ##
    'component_options'           ,
    'component_style'             ,
    ##
    'signal_options'              ,
    'signal_style'                ,
    ##
    'draw_axis_title'             , ## draw the titles for the axes ?
    'draw_options'                , ## generic ROOT draw options, e.g. 'same'
    ##
    'total_fit_options'           ,
    'curve_options'               ,
    ##
    'combined_background_options' ,
    'combined_background_style'   ,
    ##
    'combined_signal_options'     ,
    'combined_signal_style'       ,
    ##
    'combined_component_options'  ,
    'combined_component_style'    ,
    ## 
    'draw_combined_signal'        ,
    'draw_combined_background'    ,
    'draw_combined_component'     ,
    ##
    )
# =============================================================================
## get draw options:
#  Collect predefined keys from  the dictionary
#  @code
#  def somefunc (  ... , **kwargs ) :
#      draw_opts = draw_options ( kwargs )
#  @endcode 
def draw_options ( **kwargs ) :
    """Collect predefined keys from  the dictionary
    >>> def somefunc (  ... , **kwargs ) :
    ...      draw_opts = draw_options ( kwargs )
    """
    options = {}
    for k,v in items_loop ( kwargs ) :
        if k.lower() in keys : options[ k.lower() ] = v
        if k.lower() in ( 'draw' , 'draw_option' , 'draw_options' ) :
            if isinstance ( v , dict ) : options.update ( v ) 
    return options


# =============================================================================
## @class Style
#  Store the drawing style for the component.
#  Known attributes :
#  - LineColor
#  - LineStyle
#  - LineWidth
#  - FillColor
#  - FillStyle
#
#  The specification  is case-insensitive and underscore-blind 
#  @code
#  s = Style(  line_color = 4 , lineWidth = 3 , LineStyle =4 ) 
#  s = Style( ROOT.RooFit.LineColor ( 4 ) )
#  s = Style(  FillColor = 4 , fill_style = 1004  ) 
#  @endcode 
class Style(object):
    """Store the drawing style for the component
    - LineColor
    - LineStyle
    - LineWidth
    - FillColor
    - FillStyle
    
    The specification  is case-insensitive and underscore-blind
    
    >>> s = Style (  line_color = 4 , lineWidth = 3 , LineStyle =4 ) 
    >>> s = Style ( ROOT.RooFit.LineColor ( 4 ) )
    >>> s = Style (  FillColor = 4 , fill_style = 1004  ) 
    
    """
    def __init__ ( self  , *args , **kwargs ) :
        
        from ostap.utils.cidict import cidict
        kw = cidict ( transform = lambda k : k.lower().replace('_','') , **kwargs )

        linecolor    = kw.pop ( 'linecolor'    , ROOT.kBlack )
        linestyle    = kw.pop ( 'linestyle'    , 1           )
        linewidth    = kw.pop ( 'linewidth'    , 1           )
        fillcolor    = kw.pop ( 'fillcolor'    , None        )
        fillstyle    = kw.pop ( 'fillstyle'    , None        )
        transparency = kw.pop ( 'transparency' , 0.35        )
        
        if kw : logger.warning ("Style: Unknown arguments: %s" % kw.keys() )
        
        options = []

        self.__linecolor = None
        self.__linestyle = None
        self.__linewidth = None
        self.__fillcolor = None
        self.__fillstyle = None
        
        for i , a in enumerate ( args ) :
            
            if isinstance  ( a , ROOT.RooCmdArg ) :
                
                if   'LineColor' == a.GetName () : linecolor = a
                elif 'LineStyle' == a.GetName () : linestyle = a
                elif 'FillColor' == a.GetName () : fillcolor = a
                elif 'FillStyle' == a.GetName () : fillstyle = a
                else : logger.warning('Style:Unknown argument: %s ' % a )
                
            elif i == 0 and isinstance ( a , ROOT.TColor )                 : linecolor    = a
            elif i == 0 and isinstance ( a , int         ) and 0 <= a      : linecolor    = a
            elif i == 1 and isinstance ( a , int         ) and 0 <= a      : linestyle    = a
            elif i == 2 and isinstance ( a , int         ) and 0 <= a      : linewidth    = a
            elif i == 3 and isinstance ( a , ROOT.TColor )                 : fillcolor    = a
            elif i == 3 and isinstance ( a , int         ) and 0 <= a      : fillcolor    = a
            elif i == 4 and isinstance ( a , int         ) and 0 <= a      : fillstyle    = a
            elif i == 5 and isinstance ( a , float       ) and 0 <  a <  1 : transparency = a
            else :
                logger.warning('Style:Unknown argument #%d : %s/%s '  % ( i , a , type ( a ) ) )
                
        if isinstance ( linecolor , ROOT.TColor ) : linecolor = linecolor.GetNumber()
        if isinstance ( fillcolor , ROOT.TColor ) : fillcolor = fillcolor.GetNumber()

        if   isinstance ( linecolor , integer_types ) and 0 < linecolor :
            self.__linecolor = linecolor 
            options.append ( ROOT.RooFit.LineColor ( linecolor ) )
        elif isinstance ( linecolor , ROOT.RooCmdArg ) and 'LineColor' == linecolor.GetName () :
            self.__linecolor = linecolor.getInt(0) 
            options.append ( linecolor )
            
        if   isinstance ( linestyle , integer_types ) and 0 < linestyle :
            self.__linestyle  = linestyle 
            options.append ( ROOT.RooFit.LineStyle ( linestyle ) ) 
        elif isinstance ( linestyle , ROOT.RooCmdArg ) and 'LineStyle' == linestyle.GetName () :
            self.__linestyle = linestyle.getInt(0) 
            options.append ( linestyle  )

        if   isinstance ( linewidth , integer_types ) and 0 < linewidth :
            self.__linewidth = linewidth 
            options.append ( ROOT.RooFit.LineWidth ( linewidth ) ) 
        elif isinstance ( linewidth , ROOT.RooCmdArg ) and 'LineWidth' == linewidth.GetName () :
            self.__linewidth = linewidth.getInt(0) 
            options.append ( linewidth )

        _fillopt = False
        if   isinstance ( fillcolor , integer_types ) and 0 < fillcolor :
            self.__fillcolor = fillcolor 
            cc = ROOT.gROOT.GetColor ( fillcolor )
            if cc and 1.0 == cc.GetAlpha () :
                ## add tranparency
                fillcolor = ROOT.TColor.GetColorTransparent ( fillcolor , transparency ) 
            options.append ( ROOT.RooFit.FillColor ( fillcolor ) )
            _fillopt = True 
        elif isinstance ( fillcolor , ROOT.RooCmdArg ) and 'FillColor' == fillcolor.GetName () :
            self.__fillcolor = fillcolor.getInt(0)            
            options.append ( fillcolor  )               
            _fillopt = True 
        
        if   isinstance ( fillstyle , integer_types ) and 0 < fillstyle :
            self.__fillstyle = fillstyle 
            options.append ( ROOT.RooFit.FillStyle ( fillstyle ) )
            _fillopt = True 
        elif isinstance ( fillstyle , ROOT.RooCmdArg ) and 'FillStyle' == fillstyle.GetName () :
            self.__fillstyle = fillstyle.getInt(0)            
            options.append ( fillstyle  )
            _fillopt = True 

        if _fillopt :
            options.append ( ROOT.RooFit.VLines     (      ) )
            options.append ( ROOT.RooFit.DrawOption ( "FL" ) )

        self.__options = tuple ( options )

    def __str__ ( self ) :
        return "Style(LineColor=%s,LineStyle=%s,LineWidth=%s,FillColor=%s,FillStyle=%s)" % (
            self.linecolor ,
            self.linestyle ,
            self.linewidth ,
            self.fillcolor ,
            self.fillstyle )
    __repr__ = __str__

    @property 
    def options ( self ) :
        """``options'' : get the constructed list of options"""
        return self.__options

    @property
    def linecolor ( self ) : return self.__linecolor
    @property
    def linestyle ( self ) : return self.__linestyle
    @property
    def linewidth ( self ) : return self.__linewidth
    @property
    def fillcolor ( self ) : return self.__fillcolor
    @property
    def fillstyle ( self ) : return self.__fillstyle
        
# =============================================================================
## @class Area 
#  Define style for the area
#  @code
#  l = Area ( ROOT.kRed )
#  l = Area ( ROOT.kRed , fillstyle = 3315 ) 
#  @endcode
class Area (Style) :
    """Define style for the ares
    >>> l = Area ( ROOT.kRed )
    >>> l = Area ( ROOT.kRed , fillstyle = 3315 ) 
    """
    def __init__ ( self , color = ROOT.kRed , style = 1001 , **kwargs ) :

        from ostap.utils.cidict import cidict
        kw = cidict ( transform = lambda k : k.lower().replace('_','') , **kwargs )
        
        fillcolor = kw.pop ( 'fill_color' , color )
        fillstyle = kw.pop ( 'fill_style' , style )
        
        linecolor = kw.pop ( 'line_color' , None  ) ## ignore 
        linestyle = kw.pop ( 'line_style' , None  ) ## ignore 
        linewidth = kw.pop ( 'line_width' , None  ) ## ignore 
        
        Style.__init__( self                  , 
                        fillcolor = fillcolor ,
                        fillstyle = fillstyle , **kw)
        
# =============================================================================
## @class Line
#  Define style for the line
#  @code
#  l = Line ( ROOT.kRed )
#  l = Line ( ROOT.kRed , linestyle = 2 ) 
#  l = Line ( ROOT.kRed , linewidth = 3 ) 
#  @endcode
class Line (Style) :
    """Define the line style:
    >>> l = Line ( ROOT.kRed )
    >>> l = Line ( ROOT.kRed , linestyle = 2 ) 
    >>> l = Line ( ROOT.kRed , linewidth = 3 )     
    """
    def __init__ ( self                 ,             
                   color = ROOT.kBlack  ,
                   style = 1            , 
                   width = 1            , **kwargs ) :

        
        from ostap.utils.cidict import cidict
        kw = cidict ( transform = lambda k : k.lower().replace('_','') , **kwargs )
        
        fillcolor = kw.pop ( 'fill_color' , None  ) ## ignore 
        fillstyle = kw.pop ( 'fill_style' , None  ) ## ignore 
        
        linecolor = kw.pop ( 'line_color' , color ) ## ignore 
        linestyle = kw.pop ( 'line_style' , style ) ## ignore 
        linewidth = kw.pop ( 'line_width' , width ) ## ignore 

        Style.__init__( self                  , 
                        linecolor = linecolor ,
                        linestyle = linestyle ,
                        linewidth = linewidth , **kw  )

# =======================================================================================
## list of styles 
class Styles(object) :
    """List of styles
    """
    def __init__ ( self            ,
                   styles     = [] ,
                   linecolors = [ i for i in range ( 1 , 51 ) ] ,
                   linestyles = [ i for i in range ( 1 , 16 ) ] ,
                   fillcolors = [ i for i in range ( 1 , 51 ) ] ,
                   fillstyles = [] , 
                   linewidths = [] ) :

        _styles = []
        
        for s in styles :
            ss = None 
            if   isinstance ( s , Style          ) : ss = s 
            elif isinstance ( s , dict           ) : ss = Style ( **s )
            elif isinstance ( s , list_types     ) : ss = Style (  *s )
            elif isinstance ( s , integer_types  ) : ss = Style ( linecolor = s  )
            elif isinstance ( s , ROOT.TColor    ) : ss = Style ( linecolor = s  )
            elif isinstance ( s , ROOT.RooCmdArg ) : ss = Style ( None ,
                                                                  None ,
                                                                  None ,
                                                                  None ,
                                                                  None , s ) 
            if ss : _styles.append ( ss ) 

        self.__styles     = tuple ( _styles )
        
        self.__linecolors = linecolors
        self.__linestyles = linestyles
        self.__linewidths = linewidths
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

    def linewidth ( self , index ) :
        s = self.__get_it ( index  , self.__linewidths  )
        if   isinstance  ( s , integer_types  ) : return ROOT.RooFit.LineWidth  ( s )
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

    def __getitem__ ( self , index ) :
        
        if self.styles :
            l = len ( self.styles ) 
            return self.styles[ index % l ]

        l  = len ( self.styles )
        lc = self.linecolor ( index - l ) 
        ls = self.linestyle ( index - l ) 
        fc = self.fillcolor ( index - l ) 
        fs = self.fillstyle ( index - l ) 
        lw = self.linewidth ( index - l ) 

        ss = Style ( linecolor = lc , linestyle = ls , linewidth = lw ,
                     fillcolor = fc , fillstyle = fs )
        
        return ss

    ## the effective length 
    def __len__  ( self ) :
        """The effective length
        """
        return max ( len ( self.__styles     ) , 
                     len ( self.__linecolors ) , 
                     len ( self.__linestyles ) , 
                     len ( self.__linewidths ) , 
                     len ( self.__fillcolors ) , 
                     len ( self.__fillstyles ) )

    ## iteration over the defined styles 
    def __iter__ ( self ) :
        """Iteration over the defined styles
        """
        for i in range ( len ( self ) ) :
            yield self[ i ] 

    ## get the style for the given index 
    def __call__ ( self , index ) :
        """Get the style for the given index
        """
        ss = self[index]
        return ss.options

    def __str__ ( self ) :

        return str ( ( s for s in self )  )
            

# =============================================================================
## default drawing precision
precision = 1.e-4

# =============================================================================
## plain, default
data_options_plain              = ()

## suppress small bars at the end of error bars 
data_options_nobars             = ROOT.RooFit.MarkerStyle ( 20   ) , \
                                  ROOT.RooFit.DrawOption  ( "zp" ) , \
                                  ROOT.RooFit.MarkerSize  ( 0.75 ) 

## sutable for large number of bins: small markers 
data_options_small              = ROOT.RooFit.MarkerStyle ( 20   ) , \
                                  ROOT.RooFit.DrawOption  ( "zp" ) , \
                                  ROOT.RooFit.MarkerSize  ( 0.5  ) 

## default options 
default_data_options            = data_options_nobars

## default signal options 
default_signal_options          = ROOT.RooFit.Precision ( precision ) , 

default_background_options      = ROOT.RooFit.Precision ( precision ) , 

## default option for "components" 
default_component_options       = ROOT.RooFit.Precision ( precision ) , 

## default total fit curve : thick red orange line 
default_total_fit_options       = ( lineWidth ( 3 )                     ,
                                    lineColor ( ROOT.kOrange + 1 )      ,
                                    ROOT.RooFit.Precision ( precision ) )

## default curve : thick red orange line 
default_curve_options           = default_total_fit_options 

## default optios for cross-terms 
default_crossterm1_options      = ROOT.RooFit.Precision ( precision  ) ,   

## default optios for cross-terms 
default_crossterm2_options      = ROOT.RooFit.Precision ( precision  ) , 

## background:  thin short-dashed line
default_background2D_options    = default_background_options

# ==============================================================================
## @var draw_combined_signals
#  Should one draw the combined signals    (if any?)
draw_combined_signal            = False
# =============================================================================
## @var draw_combined_backgrounds    
#  Should one draw the combined signals    (if any?)
draw_combined_background        = True
# =============================================================================
## @var draw_combined_components
#  Should one draw the combined components (if any?)
draw_combined_component        = False 
# =============================================================================

# =============================================================================
## get the options from configuration parser 
# =============================================================================
def  get_options ( config , option , default ) :
    
    if not option in config : return default

    opts = config.get ( option , fallback = '()' )
    opts = opts.strip   ( ) 
    opts = opts.replace ('\n', ' ' ) 
    try : 
        options = eval( opts , globals() )
        if isinstance ( options , ROOT.RooCmdArg ) : options = options ,
        if isinstance ( options , list           ) : options = tuple ( options )
        if any ( not isinstance ( o , ROOT.RooCmdArg ) for o in options ) :
            raise TypeError('Invalid type  for %s : %s' % ( option , opts ) ) 
    except :
        logger.error("Can't parse/eval option %s : %s" % ( option , opts ) ) 
        options = default 

    return options
    
# =============================================================================
## The actual drawing options
# =============================================================================
import ostap.core.config as CONFIG 
data_options       = get_options (
    CONFIG.fit_draw , 'data_options'                , default_data_options         )
signal_options     = get_options (
    CONFIG.fit_draw , 'signal_options'              , default_signal_options       )
background_options = get_options (
    CONFIG.fit_draw , 'background_options'          , default_background_options   )
component_options  = get_options (
    CONFIG.fit_draw , 'component_options'           , default_component_options    )
crossterm1_options = get_options (
    CONFIG.fit_draw , 'crossterm1_options'          , default_crossterm1_options   )
crossterm2_options = get_options (
    CONFIG.fit_draw , 'crossterm2_options'          , default_crossterm2_options   )
total_fit_options  = get_options (
    CONFIG.fit_draw , 'total_fit_options'           , default_total_fit_options    )
curve_options      = get_options (
    CONFIG.fit_draw , 'curve_options'               , default_curve_options        )
background2D_options = get_options (
    CONFIG.fit_draw , 'background2D_options'        , default_background2D_options )
## combined stuff 
combined_signal_options     = get_options (
    CONFIG.fit_draw , 'combined_signal_options'     , default_signal_options       )
combined_background_options = get_options ( 
    CONFIG.fit_draw , 'combined_background_options' , default_background_options   )
combined_component_options  = get_options (
    CONFIG.fit_draw , 'combined_component_options'  , default_component_options    )


# =============================================================================
## Drawing styles
# =============================================================================

default_signal_style  = (
    Style ( linecolor = ROOT.kRed      , linewidth = 2 , fillcolor = ROOT.kRed     , fillstyle = 1001 ) ,
    Style ( linecolor = ROOT.kBlue     , linewidth = 2 , fillcolor = ROOT.kBlue    , fillstyle = 1001 ) ,
    Style ( linecolor =    8           , linewidth = 2 , fillcolor =    8          , fillstyle = 1001 ) ,    
    Style ( linecolor = ROOT.kMagenta  , linewidth = 2 , fillcolor = ROOT.kMagenta , fillstyle = 1001 ) ,
    Style ( linecolor = ROOT.kCyan     , linewidth = 2 , fillcolor = ROOT.kCyan    , fillstyle = 1001 ) ,
    Style ( linecolor = ROOT.kYellow   , linewidth = 2 , fillcolor = ROOT.kYellow  , fillstyle = 1001 ) ,
    ) 

default_background_style = (
    Line  ( linecolor = ROOT.kBlue         , linestyle = 14 ) ,
    Line  ( linecolor = ROOT.kBlue    -  9 , linestyle = 11 ) ,
    Line  ( linecolor = ROOT.kBlue    +  3 , linestyle = 12 ) ,
    Line  ( linecolor = ROOT.kBlue    -  2 , linestyle = 13 ) ,
    Line  ( linecolor = ROOT.kBlue    - 10 , linestyle =  9 ) ,
    )

default_combined_signal_style     = (
    Line  ( linecolor = ROOT.kMagenta , linestyle =  3 , linewidth = 2 ) ,
    )
default_combined_background_style = (
    Line  ( linecolor = 8             , linestyle = 15 , linewidth = 2 ) ,
    )
default_combined_component_style  = (
    Line  ( linecolor = ROOT.kYellow  , linestyle =  5 , linewidth = 2 ) ,
    )

default_component_style  = (
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3345 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3354 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3305 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3395 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3422 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3477 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3544 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3590 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3509 ) ,
    )

default_crossterm1_style = (
    Line  ( linecolor = ROOT.kMagenta +  3 , linestyle = 11 ) ,
    Line  ( linecolor = ROOT.kMagenta - 10 , linestyle = 12 ) , 
    Line  ( linecolor = ROOT.kMagenta +  3 , linestyle = 13 ) , 
    Line  ( linecolor = ROOT.kMagenta -  3 , linestyle = 14 ) ,
    )

default_crossterm2_style = (
    Line  ( linecolor = ROOT.kGreen   +  1 , linestyle = 14 ) ,
    Line  ( linecolor = ROOT.kGreen   -  1 , linestyle = 13 ) ,
    Line  ( linecolor = ROOT.kGreen   - 10 , linestyle = 12 ) ,
    Line  ( linecolor = ROOT.kGreen   +  1 , linestyle = 11 ) ,
    ) 

default_background2D_style = default_background_style 


# =============================================================================
## get the options/styles from configurtaion parser 
# =============================================================================
def  get_style ( config , style , default ) :
    
    if not style in config : return default

    opts = config.get ( style , fallback = '()' )
    opts = opts.strip   ( ) 
    opts = opts.replace ('\n', ' ' ) 
    try : 
        options = eval( opts , globals() )
        if isinstance ( options , Style ) : opts = opts ,
        if isinstance ( options , list  ) : opts = tuple ( options )
        if any ( not isinstance ( o , Style ) for o in options ) :
            raise TypeError('Invalid style type for %s : %s' % ( style , opts ) ) 
    except :
        logger.error("Can't parse %s : %s" % (  style , opts ) ) 
        options = default 

    return options

# =============================================================================
## the actual styles 
# =============================================================================
signal_style       = get_style (
    CONFIG.fit_draw , 'signal_style'       , default_signal_style       )
background_style   = get_style (
    CONFIG.fit_draw , 'background_style'   , default_background_style   )
component_style    = get_style (
    CONFIG.fit_draw , 'component_style'    , default_component_style    )
crossterm1_style   = get_style (
    CONFIG.fit_draw , 'crossterm1_style'   , default_crossterm1_style   )
crossterm2_style   = get_style (
    CONFIG.fit_draw , 'crossterm2_style'   , default_crossterm2_style   )
background2D_style = get_style (
    CONFIG.fit_draw , 'background2D_style' , default_background2D_style )
# ==============================================================================
combined_signal_style     = get_style (
    CONFIG.fit_draw , 'combined_signal_style'      , default_combined_signal_style     )
combined_background_style = get_style (
    CONFIG.fit_draw , 'combined_background_style'  , default_combined_background_style )
combined_component_style  = get_style (
    CONFIG.fit_draw , 'combined_component_style'   , default_combined_component_style  )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
