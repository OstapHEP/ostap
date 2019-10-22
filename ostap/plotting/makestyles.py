#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/plotting/makestyles.py
# Helper utilities to deal with ROOT styles 
# =============================================================================
"""Helper utilities to deal with ROOT styles 
"""
# =============================================================================
import ROOT
import ostap.plotting.color 
__all__ = (
    'StyleStore'       , ## the storage/dictionary of created/known styles
    'dump_style'       , ## dump a style into dicitonary
    'set_style'        , ## configure style from the dictionary
    'make_styles'      , ## create the styles from the configuration
    'make_ostap_style' , ## create Ostap-like style
    )
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.plotting.makestyles' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## font 
ostap_font       = 132     ## Times-Roman 
## line thickness
ostap_line_width =   1
## define style for text
ostap_label = ROOT.TText   (             )
ostap_label . SetTextFont  ( ostap_font  )
ostap_label . SetTextColor (  1          )
ostap_label . SetTextSize  (  0.04       )
ostap_label . SetTextAlign ( 12          )
## define style of latex text
ostap_latex = ROOT.TLatex   () 
ostap_latex . SetTextFont   ( ostap_font )
ostap_latex . SetTextColor  ( 1          )
ostap_latex . SetTextSize   ( 0.04       )
ostap_latex . SetTextAlign  ( 12         )
# =============================================================================
## @class StyleStore
#  Store for all created/cofigures styles
class StyleStore(object) :
    """Store for all created/cofigures styles
    """
    __styles = {}
    @classmethod
    def styles ( kls ) :
        return kls.__styles 
    
# =============================================================================
## get the essential methods of class ROOT.TStyle
#  @code
#  getters, setters, special = style_methods() 
#  @endcode
def style_methods () :
    """Get the essential methods of class ROOT.TStyle
    >>> getters, setters, special = style_methods() 
    """
    # The style getters 
    _getters  = set   ( [ i[3:] for i in dir ( ROOT.TStyle ) if i.startswith ( 'Get' ) ] ) 
    _getters_ = set   ( [ i[3:] for i in dir ( ROOT.TNamed ) if i.startswith ( 'Get' ) ] )
    _getters  = list  ( _getters - _getters_ )
    _getters  . sort  ( )
    _getters  = tuple ( _getters ) 
    
    # The style setters  
    _setters  = set   ( [ i[3:] for i in dir ( ROOT.TStyle ) if i.startswith ( 'Get' ) ] ) 
    _setters_ = set   ( [ i[3:] for i in dir ( ROOT.TNamed ) if i.startswith ( 'Set' ) ] )
    _setters  = list  ( _setters - _setters_ )
    _setters  . sort  ( )
    _setters  = tuple ( _setters )

    # 
    _setters_int   = set()
    _setters_float = set()
    _setters_str   = set()
    
    _style = ROOT.TStyle('tmp_style','Helper style')
    for s in _setters :

        fun = getattr ( _style , 'Set' + s , None )
        if not fun : continue
        
        try :
            fun ( 0.0 )
            _setters_float.add ( s )
            continue
        except :
            pass

        try :
            fun ( 1 )
            _setters_int.add   ( s )
            continue
        except :
            pass

        try :
            fun ( '' )
            _setters_str.add   ( s )
            continue
        except :
            pass
        
    del _style
    
    # special methods
    _special = (
        ## very special 
        'LineStyleString' ,
        'AttDate'         ,
        'PaperSize'       ,
        'ColorPalette'    , 
        ## not so special 
        'AxisColor'       ,
        'TickLength'      ,
        'Ndivisions'      ,                    
        'LabelColor'      , 'LabelFont'   , 'LabelOffset' , 'LabelSize' ,
        'TitleColor'      , 'TitleFont'   , 'TitleOffset' , 'TitleSize' )
    #
    return _getters , ( _setters_float , _setters_int , _setters_str ) , _special

# =============================================================================
##  th especial methods 
style_getters , style_setters, style_special = style_methods () 
# =============================================================================
## dump the style to the dictionary
#  @code
#  style = ...
#  conf  = dump_style ( style )
#  conf  = style.dump () ## ditto
#  conf  = style.get  () ## ditto
#  @endcode 
def dump_style ( style ) :
    """Dump the style to the dictionary
    >>> style = ...
    >>> conf  = dump_style ( style )
    >>> conf  = style.dump () ## ditto
    >>> conf  = style.get  () ## ditto    
    """
    config = {} 
    ## regular attributes 
    for g in style_getters :
        
        if g in style_special : continue
        
        fun = getattr ( style , 'Get' + g , None )
        if not fun :  return
        config [ g ] = fun ()

    ## half-special attributes 
    for attr in ( 'AxisColor'   ,
                  'TickLength'  ,
                  'Ndivisions'  ,                    
                  'LabelColor'  , 'LabelFont'   , 'LabelOffset' , 'LabelSize' ,
                  'TitleColor'  , 'TitleFont'   , 'TitleOffset' , 'TitleSize' ) :
        
        for axis in ( 'X' , 'Y' , 'Z' ) :
            
            fun = getattr ( style , 'Get' + attr ) 
            config [ '%s_%s'  %  ( attr , axis ) ] = fun ( axis )

    ## very special attribute
    import array
    x = array.array('f',[0] )
    y = array.array('f',[0] )
    style.GetPaperSize ( x , y )
    config ['PaperSize_X' ] = x[0]
    config ['PaperSize_Y' ] = y[0]

    ## very special attribute
    for i in range(31) :
        l = style.GetLineStyleString(i)
        l = l.strip()
        if l : config [ 'LineStyleString_%s' % i ] = l 
        
    return config


# =============================================================================
## Set the style from the configuration dictionary
#  @code
#  config = ...
#  style  = ...
#  set_style ( style , config )
#  style.set ( config ) ##   ditto 
#  @endcode
def set_style ( style , config ) :
    """Set the style from the configurtaion dictionary
    >>> config = ...
    >>> style  = ...
    >>> set_style ( style , config )
    >>> style.set ( config ) ## ditto 
    """
    
    for attr in style_setters [0] :
        
        if not attr in config : continue

        try : 
            value  = float ( config [ attr ] )
            setter = getattr ( style , 'Set' + attr )
            setter ( value )
            logger.debug  ("Set (float) attribute %s/%s/%s " %  ( attr , config[attr] , value ) ) 
        except :
            logger.warning("Can't set (float) attribute %s/%s, skip " %  ( attr , config[attr] ) ) 
            pass 

    for attr in style_setters [1] :  
        
        if not attr in config : continue

        try  :
            value  = int ( config [ attr ] )
            setter = getattr ( style , 'Set' + attr )
            setter ( value )
            logger.debug  ("Set (int)   attribute %s/%s/%s " %  ( attr , config[attr] , value ) ) 
        except:
            logger.warning("Can't set (int)   attribute %s/%s, skip " %  ( attr , config[attr] ) ) 
            pass 
            
    for attr in style_setters [2] :
        
        if not attr in config : continue

        try : 
            value  = config [ attr ] 
            setter = getattr ( style , 'Set' + attr )
            setter ( value )
            logger.debug  ("Set (str)   attribute %s/%s/%s " %  ( attr , config[attr] , value ) ) 
        except :
            logger.warning("Can't set (str)   attribute %s/%s, skip " %  ( attr , config[attr] ) ) 
            pass 

    ## special attributes  
    for axis in ( 'X' , 'Y' , 'Z' ) :

        key = 'AxisColor_%s'     % axis
        try :
            if key in config : style.SetAxisColor     ( int   ( config [ key ] ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 

        key = 'TickLength_%s'    % axis
        try : 
            if key in config : style.SetTickLength    ( float ( config [ key ] ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
            

        key = 'Ndivisions_%s'    % axis
        try : 
            if key in config : style.SetNdivisions    ( int   ( config [ key ] ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'LabelColor_%s'    % axis 
        try : 
            if key in config : style.SetLabelColor    ( int   ( config [ key ] ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'LabelFont_%s'     % axis 
        try : 
            if key in config : style.SetLabelFont     ( int   ( config [ key ] ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'LabelOffset_%s'   % axis 
        try : 
            if key in config : style.SetLabelOffset   ( float ( config [ key ] ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'LabelSize_%s'     % axis 
        try : 
            if key in config : style.SetLabelSize     ( float ( config [ key ] ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'TitleColor_%s'    % axis 
        try : 
            if key in config : style.SetTitleColor    ( int   ( config [ key ] ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'TitleFont_%s'     % axis 
        try : 
            if key in config : style.SetTitleFont     ( int   ( config [ key ] ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'TitleOffset_%s'   % axis 
        try : 
            if key in config : style.SetTitleOffset   ( float ( config [ key ] ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'TitleSize_%s'     % axis 
        try : 
            if key in config : style.SetTitleSize     ( float ( config [ key ] ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
    ## very special attribute 
    if 'PaperSize_X' in config and 'PaperSize_Y' in config :
        key = 'PaperSize/1'
        try :
            style.SetPaperSize ( float ( config ['PaperSize_X']  ) ,
                                 float ( config ['PaperSize_Y']  ) )
        except :
            logger.warning ( "Can't set attribute %s" % key )         
    elif 'PaperSize' in config :        
        key = 'PaperSize/2'
        try :
            style.SetPaperSize ( int   ( config ['PaperSize'] ) )            
        except :
            logger.warning ( "Can't set attribute %s" % key )         

    ## one more very special attribute
    for i in range ( 31 ) :
        k = 'LineStyleString_%s' % i
        if k in config :
            style.SetLineStyleString ( i , config[k].strip() ) 
            
    return style

# ============================================================================
ROOT.TStyle.dump = dump_style 
ROOT.TStyle.get  = dump_style  
ROOT.TStyle.set  =  set_style  


# =============================================================================
## Parse the configuration and create
#  all the styles according to configuration 
def make_styles ( config = None ) :
    """Parse the configuration and create
    all the styles according to configuration 
    """
    
    if config is None : 
        import ostap.core.config as _CONFIG
        config = _CONFIG.config
        
    for key in config :
        
        if not key.upper().startswith('STYLE') :  continue
        section   = config [ key ]
        s , c , n = key.partition (':')
        if not c : continue
        ## the style name 
        name        = n.strip ( )
        description = section.get        ( 'description' , fallback = 'The style %s' % name )
        ok          = section.getboolean ( 'ostaplike'   , fallback = False )

        ## create ostap-like style 
        if ok : make_ostap_style ( name , description , section ) 
        else :
            ## generic style 
            logger.info ( 'Create Generic style  %s/%s' % ( name , description ) )             
            style       = ROOT.TStyle ( name , description )
            set_style ( style , section )
            if name in StyleStore.styles() :
                logger.warning ( "The configuration %s replaced" % name  ) 
            StyleStore.styles().update ( { name : style } ) 
            
# ==============================================================================
def get_float ( config , name , default ) :
    try :
        if hasattr ( config , 'getfloat' ) :
            value = config.getfloat ( name , fallback = default )
        else : value = config.get   ( name , default ) 
        return float ( value )
    except :
        return default 

# =============================================================================
def get_int    ( config , name , default ) :
    
    try :
        if hasattr ( config , 'getint') :         
            value = config.getint ( name , fallback = default )
        else : value = config.get ( name , default ) 
        return int ( value )
    except :
        return default 

# =============================================================================
def get_str    ( config , name , default ) :
    
    try :
        if hasattr ( config , 'getint') :         
            value = config.get ( name , fallback = default )
        else : value = config.get ( name , default ) 
        return str ( value )
    except :
        return default 

# ============================================================================
## make Ostap-like style
def make_ostap_style ( name                           ,
                       description = 'The Style'      , 
                       config      = {}               ,
                       colz        = False            ,
                       scale       = 1.0              , 
                       font        = ostap_font       ,
                       line_width  = ostap_line_width ) :

    description = config.get ( 'description' , 'The Style' )
    
    conf  = {}
    conf.update ( config )

    conf [ 'FrameBorderMode'  ] = get_int ( config , 'FrameBorderMode'  , 0 )
    conf [ 'CanvasBorderMode' ] = get_int ( config , 'CanvasBorderMode' , 0 ) 
    conf [ 'PadBorderMode'    ] = get_int ( config , 'PadBorderMode'    , 0 ) 
    
        
    conf [ 'PadColor'         ] = get_int ( config , 'PadColor'         , 0 )
    conf [ 'CanvasColor'      ] = get_int ( config , 'CanvasColor'      , 0 )
    conf [ 'StatColor'        ] = get_int ( config , 'StatColor'        , 0 )

    
    if 'PaperSize_X' in config  or 'PaperSize_Y' in config :
        conf ['PaperSize_X' ] = get_float ( config , 'PaperSize_X' , 20 )
        conf ['PaperSize_Y' ] = get_float ( config , 'PaperSize_Y' , 26 )
    else :
        
        a = str (  config.get ( 'PaperSize' ) ).upper()         
        if   'A4'     in a :  conf [ 'PaperSize' ] = ROOT.TStyle.kA4      
        elif 'US'     in a :  conf [ 'PaperSize' ] = ROOT.TStyle.kUSletter
        elif 'LETTER' in a :  conf [ 'PaperSize' ] = ROOT.TStyle.kUSletter 
        else :  conf ['PaperSize'   ] = get_int ( config , 'PaperSize' , ROOT.TStyle.kA4 )
            
    conf [ 'PadTopMargin'      ] = get_float ( config , 'PadTopMargin'    , 0.04                   ) 
    conf [ 'PadRightMargin'    ] = get_float ( config , 'PadRightMargin'  , 0.14 if colz else 0.04 ) 
    conf [ 'PadLeftMargin'     ] = get_float ( config , 'PadLeftMargin'   , 0.14                   ) 
    conf [ 'PadBottomMargin'   ] = get_float ( config , 'PadBottomMargin' , 0.14                   ) 
    
    conf [ 'TextFont'          ] = get_int   ( config , 'TextFont'        , font         ) 
    conf [ 'TextSize'          ] = get_float ( config , 'FontSize'        , 0.08 * scale ) 
    
    conf [ 'LabelFont_X'       ] = get_int   ( config , 'LabelFont_X' , font )
    conf [ 'LabelFont_Y'       ] = get_int   ( config , 'LabelFont_Y' , font )
    conf [ 'LabelFont_Z'       ] = get_int   ( config , 'LabelFont_Z' , font )
    
    conf [ 'LabelSize_X'       ] = get_float ( config , 'LabelSize_X' , 0.05 * scale ) 
    conf [ 'LabelSize_Y'       ] = get_float ( config , 'LabelSize_Y' , 0.05 * scale ) 
    conf [ 'LabelSize_Z'       ] = get_float ( config , 'LabelSize_Z' , 0.05 * scale ) 
    
    conf [ 'TitleFont_X'       ] = get_int   ( config , 'TitleFont_X' , font ) 
    conf [ 'TitleFont_Y'       ] = get_int   ( config , 'TitleFont_Y' , font ) 
    conf [ 'TitleFont_Z'       ] = get_int   ( config , 'TitleFont_Z' , font ) 
    
    conf [ 'TitleSize_X'       ] = get_float ( config , 'TitleSize_X' , -1 ) 
    conf [ 'TitleSize_Y'       ] = get_float ( config , 'TitleSize_Y' , 0.05 * scale  ) 
    conf [ 'TitleSize_Z'       ] = get_float ( config , 'TitleSize_Z' , 0.05 * scale  ) 
    
    conf [ 'LineWidth'         ] = get_int   ( config , 'LineWidth'     , line_width  )
    conf [ 'FrameWidth'        ] = get_int   ( config , 'FrameWidth'    , line_width  ) 
    conf [ 'HistLineWidth'     ] = get_int   ( config , 'HistLineWidth' , line_width  ) 
    conf [ 'FuncWidth'         ] = get_int   ( config , 'FuncWidth'     , line_width  ) 
    conf [ 'GridWidth'         ] = get_int   ( config , 'FuncWidth'     , line_width  ) 
    
    conf [ 'MarkerStyle'       ] = get_int   ( config , 'MarkerStyle'   , 20  ) 
    conf [ 'MarkerSize'        ] = get_float ( config , 'MarkerSize'    , 1.2 )
    
    conf [ 'LabelOffset'       ] = get_float ( config , 'LabelOffset'   , 0.015 ) 


    conf [ 'StatFormat'        ] = get_str   ( config , 'StatFormat'    , '6.3g')
    
    conf [ 'OptTitle'          ] = get_int   ( config , 'OptTitle'      , 0    )
    conf [ 'OptFit'            ] = get_int   ( config , 'OptFit'        , 0    )
    conf [ 'OptStat'           ] = get_int   ( config , 'OptStat'       , 0    )

    conf [ 'LegendFont'        ] = get_int   ( config , 'LegendFont'    , font )

    
    ## size of small lines at the end of error bars
    conf [ 'EndErrorsSize'     ] = get_float ( config , 'EndErrorsSize'  , 3 ) 

    ## statistics box
    conf [ 'StatBorderSize'    ] = get_int   ( config , 'StatBorderSize' , 0            )
    conf [ 'StatFont'          ] = get_int   ( config , 'StatFont'       , font         ) 
    conf [ 'StatFontSize'      ] = get_float ( config , 'StatFontSize'   , 0.05 * scale ) 

    conf [ 'StatX'             ] = get_float ( config , 'StatX'          , 0.9  ) 
    conf [ 'StatY'             ] = get_float ( config , 'StatY'          , 0.9  ) 
    conf [ 'StatW'             ] = get_float ( config , 'StatW'          , 0.25 ) 
    conf [ 'StatH'             ] = get_float ( config , 'StatH'          , 0.14 )
    
    conf [ 'PadTickX'          ] = get_int   ( config , 'PadTickX'       , 1 ) 
    conf [ 'PadTickY'          ] = get_int   ( config , 'PadTickY'       , 1 ) 
    
    conf [ 'Ndivisions_X'      ] = get_int   ( config , 'Ndivisions_X'   , 505 )
    conf [ 'Ndivisions_Y'      ] = get_int   ( config , 'Ndivisions_Y'   , 510 )
    conf [ 'Ndivisions_Z'      ] = get_int   ( config , 'Ndivisions_Z'   , 510 )
    
    ##  dark-body radiator pallete
    conf [ 'Palette'           ] = get_int   ( config , 'Paletter' , ROOT.kDarkBodyRadiator )
    conf [ 'NumberContours'    ] = get_int   ( config , 'NumberContours' , 255 )


    conf [ 'LineStyleString_2'  ] = "12 12"
    conf [ 'LineStyleString_11' ] = "76 24"
    conf [ 'LineStyleString_12' ] = "60 16 8 16"
    conf [ 'LineStyleString_13' ] = "168 32"
    conf [ 'LineStyleString_14' ] = "32  32"
    conf [ 'LineStyleString_15' ] = "80  20"

    ## create the style 
    style       = ROOT.TStyle ( name , description )
    set_style    ( style , conf ) 
    logger.debug ('Create Ostap   style %s' % style.GetName() )

    if name in StyleStore.styles() :
        logger.info ( "The configuration %s replaced" % name  ) 
    StyleStore.styles().update ( { name : style } ) 

    if  name.startswith('Style')  :
        nname  = name[5:]
        if nname in StyleStore.styles() :
            logger.info ( "The configuration %s replaced" % nname  ) 
        StyleStore.styles().update ( { nname : style } )
        
    return style

# =============================================================================
## read the configuration files and create the styles 
make_styles () 
    
# =============================================================================
_decorated_classes_ = (
    ROOT.TStyle ,
    )
_new_methods_   = (
    #
    ROOT.TStyle.dump ,
    ROOT.TStyle.get  ,
    ROOT.TStyle.set  ,
    )
# =============================================================================
if '__main__' == __name__ :
         
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
