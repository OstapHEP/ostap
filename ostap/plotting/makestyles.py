#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/plotting/makestyles.py
# Helper utilities to deal with ROOT styles 
# =============================================================================
"""Helper utilities to deal with ROOT styles 
"""
# =============================================================================
import ostap.plotting.color
from   ostap.utils.cidict import cidict
from   ostap.core.core    import cidict_fun
from   ostap.utils.utils  import classprop 
import ROOT, ctypes 
# =============================================================================
__all__ = (
    'StyleStore'       , ## the storage/dictionary of created/known styles
    'dump_style'       , ## dump a style into dicitonary
    'set_style'        , ## configure style from the dictionary
    'make_styles'      , ## create the styles from the configuration
    'make_ostap_style' , ## create Ostap-like style
    'canvas_width'     , ## (default) canvas width
    'canvas_height'    , ## (default) canvas height
    'margin_top'       , ## (default) top margin 
    'margin_bottom'    , ## (default) top margin 
    'margin_left'      , ## (default) top margin 
    'margin_right'     , ## (default) top margin 
    )
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.plotting.makestyles' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## Canvas size 
margin_top    = 0.05
margin_bottom = 0.12
margin_left   = 0.12
margin_right  = 0.07
# ==============================================================================
import ostap.core.config as OCC
width_  = OCC.canvas.get ( 'Width'       , fallback = '1000' )
try               : width_  = int  ( width_ ) 
except ValueError : width_  = 1000
canvas_width  = width_  if 20 <= width_  else 1000
# 
height_ = OCC.canvas.get ( 'Height'      , fallback =  '800' )
try               : height_ = int  ( height_ ) 
except ValueError : height_ = 800
canvas_height = height_ if 20 <= height_ else  800
##
w2h  =  float ( canvas_width  ) / float ( canvas_height )
h2w  =  float ( canvas_height ) / float ( canvas_width  )

tmp_    = OCC.canvas.get ( 'MarginRight'   , fallback = '%f'  % margin_right       )
try               : tmp_  = float ( tmp_ ) 
except ValueError : tmp_  = margin_right
if   0 <      tmp_ < 1            : margin_right = tmp_
elif 0 < -1 * tmp_ < canvas_width : margin_right = abs ( 1.0 * tmp_ / canvas_width )
##
tmp_    = OCC.canvas.get ( 'MarginBottom'   , fallback = '%f'  % margin_bottom ) 
try               : tmp_  = float ( tmp_ ) 
except ValueError : tmp_  = margin_bottom
if   0 <      tmp_ < 1             : margin_bottom = tmp_
elif 0 < -1 * tmp_ < canvas_height : margin_bottom = abs ( 1.0 * tmp_ / canvas_height )  

tmp_    = OCC.canvas.get ( 'MarginLeft' , fallback = '%f' %  margin_left  )
try               : tmp_  = float ( tmp_ ) 
except ValueError : tmp_  = margin_left 
if   0 <      tmp_ < 1            : margin_left = tmp_
elif 0 < -1 * tmp_ < canvas_width : margin_left = abs ( 1.0 * tmp_ / canvas_width )  
#

tmp_    = OCC.canvas.get ( 'MarginTop'   , fallback = '%f'  % margin_top ) 
try               : tmp_  = float ( tmp_ ) 
except ValueError : tmp_  = margin_top 
if   0 <      tmp_ < 1             : margin_top = tmp_
elif 0 < -1 * tmp_ < canvas_height : margin_top = abs ( 1.0 * tmp_ / canvas_weight )  
#

logger.verbose ( 'Default Canvas parameters: ')
logger.verbose ( '  canvas_height : %s ' % canvas_height )
logger.verbose ( '  canvas_width  : %s ' % canvas_width  )
logger.verbose ( '  margin_top    : %s ' % margin_top    )
logger.verbose ( '  margin_bottom : %s ' % margin_bottom )
logger.verbose ( '  margin_left   : %s ' % margin_left   )
logger.verbose ( '  margin_right  : %s ' % margin_right  )

# =============================================================================
## default font 
ostap_font       = 132     ## Times-Roman 
## line thickness
ostap_line_width =   1
# =============================================================================
tmp_    = OCC.canvas.get ( 'TextFont'   , fallback = '%d'  % ostap_font ) 
try               : tmp_  = int ( tmp_ ) 
except ValueError : tmp_  = ostap_font
ostap_font = tmp_
# =============================================================================
tmp_    = OCC.canvas.get ( 'LineWidth'   , fallback = '%d'  % ostap_line_width ) 
try               : tmp_  = int ( tmp_ ) 
except ValueError : tmp_  = ostap_line_width
ostap_line_width = tmp_
# ==============================================================================
## define style for text
ostap_label = ROOT.TText    (             )
ostap_label . SetTextFont   ( ostap_font  )
ostap_label . SetTextColor  (  1          )
ostap_label . SetTextSize   (  0.04       )
ostap_label . SetTextAlign  ( 12          )
## define style of latex text
ostap_latex = ROOT.TLatex   () 
ostap_latex . SetTextFont   ( ostap_font  )
ostap_latex . SetTextColor  ( 1           )
ostap_latex . SetTextSize   ( 0.04        )
ostap_latex . SetTextAlign  ( 12          )
# ==============================================================================
## @class StyleStore
#  Store for all created/configured styles
class StyleStore(object) :
    """Store for all created/cofigures styles
    """
    __styles = {}
    @classprop 
    def styles ( kls ) :
        return kls.__styles 

# =============================================================================
## get the ROOT style by name
#  @code
#  style = root_style ( 'MyStyle' )
#  if style : ...
#  @endcode
def root_style ( name ) :
    """ Get the ROOT style by name
    >>> style = root_style ( 'MyStyle' )
    >>> if style : ...
    """
    groot      = ROOT.ROOT.GetROOT ()        
    style_list = groot.GetListOfStyles()
    for style in style_list :
        if style and style.GetName () == name  : return style
    return None 
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
        'LineStyleString'    ,
        'AttDate'            ,
        'PaperSize'          ,
        'ColorPalette'       ,
        'MarkerLineWidth'    , 
        'MarkerStyleBase'    ,
        'ExponentOffset'     , ## new special stuff (August 2024)
        'TextSizePercent'    , ## new special stuff (August 2024)
        ## not so special 
        'AxisColor'       ,
        'TickLength'      ,
        'Ndivisions'      ,                    
        'LabelColor'      , 'LabelFont'   , 'LabelOffset' , 'LabelSize' ,
        'TitleColor'      , 'TitleFont'   , 'TitleOffset' , 'TitleSize' )
    #
    return _getters , ( _setters_float , _setters_int , _setters_str ) , _special

# =============================================================================
##  the special methods 
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
        if not fun :  continue 
        config [ g ] = fun ()

    ## half-special attributes 
    for attr in ( 'AxisColor'   ,
                  'TickLength'  ,
                  'Ndivisions'  ,                    
                  'LabelColor'  , 'LabelFont'   , 'LabelOffset' , 'LabelSize' ,
                  'TitleColor'  , 'TitleFont'   , 'TitleOffset' , 'TitleSize' ) :
        
        for axis in ( 'X' , 'Y' , 'Z' ) :
            
            fun = getattr ( style , 'Get' + attr , None) 
            if fun : config [ '%s_%s'  %  ( attr , axis ) ] = fun ( axis )

    ## very special attribute
    x = ctypes.c_float()
    y = ctypes.c_float()
    style.GetPaperSize ( x , y )
    config ['PaperSize_X' ] = x.value 
    config ['PaperSize_Y' ] = y.value 

    ## special... 
    config [ 'TGaxisMaxDigits' ] = ROOT.TGaxis.GetMaxDigits()
    
    ## very special attribute
    for i in range(31) :
        l = style.GetLineStyleString(i)
        l = l.strip()
        if l : config [ 'LineStyleString_%s' % i ] = l 
        
    return config

# =============================================================================
## Dump the style as a table 
def table_style ( style , prefix = '' , title = '' ) : 
    """Dump the style as a table"""

    conf = dump_style ( style )
    
    for i in range ( 31 ) :
        key = 'LineStyleString_%s' % i  
        fmt = style.GetLineStyleString ( i )
        if fmt : conf [ key ] = fmt 

    table = [ ( '#' , 'Parameter' , 'value' ) ] 
    for i, key in enumerate ( sorted ( conf ) , start = 1 ) : 
        
        value = conf[ key ] 
        row   = '%3d' % i , key , '%s' % value 
        table.append ( row )

    title = title if title else 'Style %s/%s' % ( style.GetName() , style.GetTitle () )
    import ostap.logger.table as T
    return T.table ( table , title = title , prefix = prefix , alignment = 'll' )

ROOT.TStyle.table = table_style
# =============================================================================
## extr aattribvutes 
extra = 'NumberOfColors' , 'showeditor', 'showeventstatus', 'showtoolbar', 'basestyle', 'ostaplike'
# =============================================================================
## Set the style from the configuration dictionary
#  @code
#  config = ...
#  style  = ...
#  set_style ( style , config )
#  style.set ( config ) ##   ditto 
#  @endcode
def set_style ( style , config , base_style = '' , **kwargs ) :
    """Set the style from the configurtaion dictionary
    >>> config = ...
    >>> style  = ...
    >>> set_style ( style , config )
    >>> style.set ( config ) ## ditto 
    """

    base_style  = base_style if base_style else config.get ( 'base_style' , base_style )
    
    base_style  = root_style ( base_style ) if base_style else None
    base_config = dump_style ( base_style ) if base_style else {}  
            
    conf = cidict ( transform = cidict_fun )    
    conf.update ( base_config ) ## the base configuration 
    conf.update ( config      ) ## 
    conf.update ( kwargs      ) 
          
    changed = {}
    
    for attr in style_setters  [ 0 ] :

        if not attr in conf : continue

        try :
            
            value     = float   ( conf.pop ( attr ) )
            setter    = getattr ( style , 'Set' + attr )

            old_value = None 
            try : 
                if attr in style_getters :
                    getter    = getattr ( style , 'Get' + attr )
                    old_value = getter () 
            except :
                pass
            
            setter ( value )
            
            if not old_value is None :
                changed [ attr ] = old_value 
                    
            logger.debug  ("Set (float) attribute `%s' to be %s " %  ( attr , value ) )
            
        except :
            
            logger.warning("Can't set (float) attribute %s, skip " %  attr  )
            

    for attr in style_setters [ 1 ] :  
        
        if not attr in conf : continue

        try  :
            
            value  = int ( conf.pop ( attr )  )
            setter = getattr ( style , 'Set' + attr )
            
            old_value = None 
            try : 
                if attr in style_getters :
                    getter    = getattr ( style , 'Get' + attr )
                    old_value = getter () 
            except :
                pass

            setter ( value )

            if not old_value is None :
                changed [ attr ] = old_value 
            
            logger.debug  ("Set (int)   attribute `%s' to %s " %  ( attr , value ) )
            
        except:
            
            logger.warning("Can't set (int)   attribute %s/%s, skip " %  attr )
        
    for attr in style_setters [ 2 ] :
        
        if not attr in conf : continue

        try :
            
            value  = conf.pop  ( attr ) 
            setter = getattr ( style , 'Set' + attr )

            old_value = None 
            try : 
                if attr in style_getters :
                    getter    = getattr ( style , 'Get' + attr )
                    old_value = getter () 
            except :
                pass

            setter ( value )

            if not old_value is None :
                changed [ attr ] = old_value 
                        
            logger.debug  ("Set (str)   attribute `%s; to `%s'" %  ( attr , value ) )
            
        except :
            
            logger.warning("Can't set (str)   attribute `%s', skip " % attr  ) 


    ## half-special attributes 
    for attr in ( 'AxisColor'   ,
                  'TickLength'  ,
                  'Ndivisions'  ,                    
                  'LabelColor'  , 'LabelFont'   , 'LabelOffset' , 'LabelSize' ,
                  'TitleColor'  , 'TitleFont'   , 'TitleOffset' , 'TitleSize' ) :

        if attr in conf :
            
            x_attr = '%s_X' % attr
            y_attr = '%s_Y' % attr
            z_attr = '%s_Z' % attr

            if ( not x_attr in conf ) and \
               ( not y_attr in cong ) and \
               ( not z_attr in conf ) :
                
                value = conf.pop ( attr )
                conf [ x_attr ] = value
                conf [ y_attr ] = value
                conf [ z_attr ] = value
                

    ## special attributes  
    for axis in ( 'X' , 'Y' , 'Z' ) :

        key = 'AxisColor_%s'     % axis
        try :
            if key in conf :
                changed [ key ] = style.GetAxisColor( axis )
                style.SetAxisColor     ( int   ( conf.pop ( key ) ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 

        key = 'TickLength_%s'    % axis
        try : 
            if key in conf :
                changed [ key ] = style.GetTickLength ( axis )                
                style.SetTickLength    ( float ( conf.pop ( key ) ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 

        key = 'Ndivisions_%s'    % axis
        try : 
            if key in conf  :
                changed [ key ] = style.GetNdivisions ( axis )                                
                style.SetNdivisions    ( int   ( conf.pop ( key ) ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'LabelColor_%s'    % axis 
        try : 
            if key in conf :
                changed [ key ] = style.GetLabelColor ( axis )                                                
                style.SetLabelColor    ( int   ( conf.pop ( key ) ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 

        
        key = 'LabelFont_%s'     % axis 
        try : 
            if key in conf :
                changed [ key ] = style.GetLabelFont ( axis )                                                
                style.SetLabelFont     ( int   ( conf.pop ( key )  ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'LabelOffset_%s'   % axis 
        try : 
            if key in conf :
                changed [ key ] = style.GetLabelOffset ( axis )                                                
                style.SetLabelOffset   ( float ( conf.pop ( key ) ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'LabelSize_%s'     % axis 
        try : 
            if key in conf :
                changed [ key ] = style.GetLabelSize ( axis )                                                
                style.SetLabelSize     ( float ( conf.pop ( key ) ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'TitleColor_%s'    % axis 
        try : 
            if key in conf :
                changed [ key ] = style.GetTitleColor ( axis )                                                
                style.SetTitleColor    ( int   ( conf.pop ( key ) ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'TitleFont_%s'     % axis 
        try : 
            if key in conf :
                changed [ key ] = style.GetTitleFont ( axis )                                                
                style.SetTitleFont     ( int   ( conf.pop ( key ) ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'TitleOffset_%s'   % axis 
        try : 
            if key in conf :
                changed [ key ] = style.GetTitleOffset ( axis )                                                
                style.SetTitleOffset   ( float ( conf.pop ( key ) ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 
        
        key = 'TitleSize_%s'     % axis 
        try : 
            if key in conf :
                changed [ key ] = style.GetTitleSize ( axis )                                                
                style.SetTitleSize     ( float ( conf.pop ( key ) ) , axis )
        except :
            logger.warning ( "Can't set attribute %s" % key ) 

    ## very special attribute 
    if 'PaperSize_X' in conf and 'PaperSize_Y' in conf :
        key = 'PaperSize/1'
        try :
            x = ctypes.c_float()
            y = ctypes.c_float()
            style.GetPaperSize ( x , y )
            changed [ 'PaperSize_X' ] = x.value 
            changed [ 'PaperSize_Y' ] = y.value 
            style.SetPaperSize ( float ( conf.pop ( 'PaperSize_X' ) ) ,
                                 float ( conf.pop ( 'PaperSize_Y' ) ) )
        except :
            logger.warning ( "Can't set attribute %s" % key )
            
    elif 'PaperSize' in conf :        
        key = 'PaperSize/2'
        try :            
            x = ctypes.c_float()
            y = ctypes.c_float()
            style.GetPaperSize ( x , y )
            changed [ 'PaperSize_X' ] = x.value 
            changed [ 'PaperSize_Y' ] = y.value 
            style.SetPaperSize ( int   ( conf.pop ( 'PaperSize' ) ) )            
        except :
            logger.warning ( "Can't set attribute %s" % key )         

    ## one more very special attribute
    for i in range ( 31 ) :
        k = 'LineStyleString_%s' % i
        if k in conf :
            changed [ key ] = style.GetLineStyleString ( i )                         
            style.SetLineStyleString ( i , conf.pop ( k ) .strip() ) 

    if 'TGaxisMaxDigits' in conf : 
        ROOT.TGaxis.SetMaxDigits ( conf.pop ( 'TGaxisMaxDigits' ) )
        
    if 'palette' in conf :
        style.SetPalette ( conf.pop ( 'palette' ) )

    ## Extra attribtes 
    for e in extra :
        if e in conf : conf.pop ( e ) 
        
    if conf :
        logger.warning ( "set_style: unprocessed parameters: %s" % list ( conf.keys() ) )
        
    return changed 

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
        ok          = section.getboolean ( 'ostap_like'  , fallback =  True )

        ## create ostap-like style 
        if ok : style = make_ostap_style ( name , description , section )
        else  :
            ## generic style
            style = root_style ( name )
            if not style : 
                logger.debug ( 'Create new generic style  %s/%s' % ( name , description ) )             
                style       = ROOT.TStyle ( name , description )
                
            set_style ( style , section )
            
        if name in StyleStore.styles :
            logger.warning ( "The configuration %s replaced" % name  ) 
        StyleStore.styles.update ( { name : style } )
        
        if  name.startswith('Style')  :
            nname  = name[5:]
            if nname in StyleStore.styles :
                logger.info ( "The configuration %s replaced" % nname  ) 
                StyleStore.styles.update ( { nname : style } )
        
            
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
def get_bool    ( config , name , default ) :
    
    try :
        if hasattr ( config , 'getboolean') :         
            value = config.getboolean ( name , fallback = default )
        else : value = config.get ( name , default ) 
        return bool ( value )
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
def make_ostap_style ( name                      ,
                       description = 'The Style' ,   
                       config      = {}          ,
                       base_style  = ''          , **kwargs ) :
    

    kw = cidict ( transform = cidict_fun )
    kw.update   ( kwargs ) 
    
    colz       = kw.pop ( 'colz'        , get_bool  ( config , 'colz'       , False            ) )
    scale      = kw.pop ( 'scale'       , get_float ( config , 'scale'      , 1.0              ) )
    font       = kw.pop ( 'font'        , get_int   ( config , 'font'       , ostap_font       ) )
    line_width = kw.pop ( 'line_width'  , get_int   ( config , 'line_width' , ostap_line_width ) )
    
    if kw : logger.warning ("make_ostap_style: unprocessed keys: %s" % kw ) 
    
    description = config.get ( 'description' , description )
    
    conf  = {}
    conf.update ( config )  ## own configuration 
    
    conf [ 'AxisColor_X'       ] = get_int   ( config , 'AxisColor_X'         , 1   )
    conf [ 'AxisColor_Y'       ] = get_int   ( config , 'AxisColor_Y'         , 1   )
    conf [ 'AxisColor_Z'       ] = get_int   ( config , 'AxisColor_Z'         , 1   )

    conf [ 'BarOffset'         ] = get_float ( config , 'BarOffset'           , 0.0 )
    conf [ 'BarWidth'          ] = get_float ( config , 'BarWidth'            , 1.0 )
    

    conf [ 'CanvasBorderMode'  ] = get_int   ( config , 'CanvasBorderMode'    , 0      ) 
    conf [ 'CanvasBorderSize'  ] = get_int   ( config , 'CanvasBorderSize'    , 2      ) 
    conf [ 'CanvasColor'       ] = get_int   ( config , 'CanvasColor'         , 0      ) 
    conf [ 'CanvasDefH'        ] = get_int   ( config , 'CanvasDefH'          , canvas_height ) 
    conf [ 'CanvasDefW'        ] = get_int   ( config , 'CanvasDefW'          , canvas_width  ) 
    conf [ 'CanvasDefX'        ] = get_int   ( config , 'CanvasDefX'          , 10     ) 
    conf [ 'CanvasDefY'        ] = get_int   ( config , 'CanvasDefY'          , 10     ) 

    conf [ 'DateX'             ] = get_float ( config , 'DateX'               , 0.01   ) 
    conf [ 'DateY'             ] = get_float ( config , 'DateY'               , 0.01   ) 

    conf [ 'DrawBorder'        ] = get_int   ( config , 'DrawBorder'          , 0      ) 

    conf [ 'EndErrorSize'      ] = get_float ( config , 'EndErrorSize'        , 2.0    * scale )
    conf [ 'ErrorX'            ] = get_float ( config , 'ErrorX'              , 0.5    )

    conf [ 'FitFormat'         ] = get_str   ( config , 'FitFormat'           , '5.4g' ) 
    
    conf [ 'FrameBorderMode'   ] = get_int   ( config , 'FrameBorderMode'     , 0      )
    conf [ 'FrameBorderSize'   ] = get_int   ( config , 'FrameBorderSize'     , 1      )
    conf [ 'FrameFillColor'    ] = get_int   ( config , 'FrameFillColor'      , 0      )
    conf [ 'FrameFillStyle'    ] = get_int   ( config , 'FrameFillStyle'      , 1001   )
    conf [ 'FrameLineColor'    ] = get_int   ( config , 'FrameLineColor'      , 1      )
    conf [ 'FrameLineStyle'    ] = get_int   ( config , 'FrameLineStyle'      , 1      )
    conf [ 'FrameLineWidth'    ] = get_int   ( config , 'FrameLineWidth'      , line_width  )

    conf [ 'FuncColor'         ] = get_int   ( config , 'FuncColor'           , 2 )
    conf [ 'FuncStyle'         ] = get_int   ( config , 'FuncStyle'           , 1 )
    conf [ 'FuncWidth'         ] = get_int   ( config , 'FuncWidth'           , line_width )
    
    conf [ 'GridColor'         ] = get_int   ( config , 'GridColor'           , 1 )
    conf [ 'GridStyle'         ] = get_int   ( config , 'GridStyle'           , 3 )
    conf [ 'GridWidth'         ] = get_int   ( config , 'GridWidth'           , 1 )

    conf [ 'HatchesLineWidth'  ] = get_int   ( config , 'HatchesLineWidth'    , 1     )
    conf [ 'HatchesSpacing'    ] = get_float ( config , 'HatchesSpacing'      , 1.0   )

    conf [ 'HistFillColor'     ] = get_int   ( config , 'HistFillColor'       , 0     )
    conf [ 'HistFillStyle'     ] = get_int   ( config , 'HistFillStyle'       , 1001  )
    conf [ 'HistLineColor'     ] = get_int   ( config , 'HistLineColor'       , 1     )
    conf [ 'HistLineStyle'     ] = get_int   ( config , 'HistLineStyle'       , 1     )
    conf [ 'HistLineWidth'     ] = get_int   ( config , 'HistLineStyle'       , line_width )

    conf [ 'HistMinimumZero'   ] = get_bool  ( config , 'HistMinimumZero'     , False )
    conf [ 'HistTopMargin'     ] = get_float ( config , 'HistTopMargin'       , 0.05  )

    conf [ 'JoinLinePS'        ] = get_int   ( config , 'JoinLinePS'          , 0.    )

    conf [ 'LabelColor_X'      ] = get_int   ( config , 'LabelColor_X'        , 1     )
    conf [ 'LabelColor_Y'      ] = get_int   ( config , 'LabelColor_Y'        , 1     )
    conf [ 'LabelColor_Z'      ] = get_int   ( config , 'LabelColor_Z'        , 1     )

    conf [ 'LabelFont_X'       ] = get_int   ( config , 'LabelFont_X'         , font  )
    conf [ 'LabelFont_Y'       ] = get_int   ( config , 'LabelFont_Y'         , font  )
    conf [ 'LabelFont_Z'       ] = get_int   ( config , 'LabelFont_Z'         , font  )

    conf [ 'LabelOffset_X'     ] = get_float ( config , 'LabelOffset_X'       , 0.015 )
    conf [ 'LabelOffset_Y'     ] = get_float ( config , 'LabelOffset_Y'       , 0.005 )
    conf [ 'LabelOffset_Z'     ] = get_float ( config , 'LabelOffset_Z'       , 0.005 ) 

    conf [ 'LabelSize_X'       ] = get_float ( config , 'LabelSize_X'         , 0.05  * scale )
    conf [ 'LabelSize_Y'       ] = get_float ( config , 'LabelSize_Y'         , 0.05  * scale )
    conf [ 'LabelSize_Z'       ] = get_float ( config , 'LabelSize_Z'         , 0.05  * scale )
     
    conf [ 'LegendBorderSize'  ] = get_int   ( config , 'LegendBorderSize'    , 4    )
    conf [ 'LegendFillColor'   ] = get_int   ( config , 'LegendFillColor'     , 0    )
    conf [ 'LegendFont'        ] = get_int   ( config , 'LegendFont'          , font )
    conf [ 'LegendTextSize'    ] = get_float ( config , 'LegendTextSize'      , 0.0  )

    conf [ 'LegoInnerR'        ] = get_float ( config , 'LegoInnerR'          , 0.5  )
    conf [ 'LineScalePS'       ] = get_float ( config , 'LineScalePS'         , 3.0  )

    for i in range ( 31 ) :
        key = 'LineStyleString_%d' % i
        fmt = get_str   ( config , key , '' )
        if fmt : conf [ key ] = fmt
        
    if not 'LineStyleString_2'  in conf : conf [ 'LineStyleString_2'  ] = " 12 12"
    if not 'LineStyleString_11' in conf : conf [ 'LineStyleString_11' ] = " 76 24"
    if not 'LineStyleString_12' in conf : conf [ 'LineStyleString_12' ] = " 60 16 8 16"
    if not 'LineStyleString_13' in conf : conf [ 'LineStyleString_13' ] = "168 32"
    if not 'LineStyleString_14' in conf : conf [ 'LineStyleString_14' ] = " 32 32"
    if not 'LineStyleString_15' in conf : conf [ 'LineStyleString_15' ] = " 80 20"
    if not 'LineStyleString_16' in conf : conf [ 'LineStyleString_16' ] = " 40 10"

    conf [ 'Ndivisions_X'      ] = get_int   ( config , 'Ndivisions_X'       , 505  )
    conf [ 'Ndivisions_Y'      ] = get_int   ( config , 'Ndivisions_Y'       , 510  )
    conf [ 'Ndivisions_Z'      ] = get_int   ( config , 'Ndivisions_Z'       , 510  )

    conf [ 'NumberContours'    ] = get_int   ( config , 'NumberContours'     , 127  )
    
    conf [ 'NumberOfColors'    ] = get_int   ( config , 'NumberOfColors'     , 255  )

    conf [ 'OptDate'           ] = get_int   ( config , 'OptDate'            , 0    )
    conf [ 'OptFile'           ] = get_int   ( config , 'OptFile'            , 0    )
    conf [ 'OptFit'            ] = get_int   ( config , 'OptFit'             , 0    )
    conf [ 'OptLogx'           ] = get_int   ( config , 'OptLogx'            , 0    )
    conf [ 'OptLogy'           ] = get_int   ( config , 'OptLogy'            , 0    )
    conf [ 'OptLogz'           ] = get_int   ( config , 'OptLogz'            , 0    )
    conf [ 'OptStat'           ] = get_int   ( config , 'OptStat'            , 0    )
    conf [ 'OptTitle'          ] = get_int   ( config , 'OptTitle'           , 0    )

    conf [ 'PadBorderMode'     ] = get_int   ( config , 'PadBorderMode'      , 0     ) 
    conf [ 'PadBorderSize'     ] = get_int   ( config , 'PadBorderSize'      , 2     ) 
    conf [ 'PadBottomMargin'   ] = get_float ( config , 'PadBottomMargin'    , margin_bottom )    
    conf [ 'PadColor'          ] = get_int   ( config , 'PadColor'           , 0     )
    conf [ 'PadGridX'          ] = get_bool  ( config , 'PadGridX'           , False )
    conf [ 'PadGridY'          ] = get_bool  ( config , 'PadGridY'           , False )
    conf [ 'PadLeftMargin'     ] = get_float ( config , 'PadLeftMargin'      , margin_left   )    
    conf [ 'PadRightMargin'    ] = get_float ( config , 'PadRightMargin'     , margin_left if colz else margin_right  )    
    conf [ 'PadTickX'          ] = get_int   ( config , 'PadTickX'           , 1     )
    conf [ 'PadTickY'          ] = get_int   ( config , 'PadTickY'           , 1     )
    conf [ 'PadTopMargin'      ] = get_float ( config , 'PadTopMargin'       , margin_top    )    

    conf [ 'PaintTextFormat'   ] = get_str   ( config , 'PaintTextFormat'    , 'g' )    

    if 'PaperSize_X' in config  or 'PaperSize_Y' in config :
        
        conf ['PaperSize_X' ] = get_float ( config , 'PaperSize_X' , 20 )
        conf ['PaperSize_Y' ] = get_float ( config , 'PaperSize_Y' , 26 )
        
    elif 'PaperSize' in config :
        
        a = str (  config.get ( 'PaperSize' ) ).upper()         
        if   'A4'     in a :  conf [ 'PaperSize' ] = ROOT.TStyle.kA4      
        elif ''       in a :  conf [ 'PaperSize' ] = ROOT.TStyle.kUSLetter
        elif 'LETTER' in a :  conf [ 'PaperSize' ] = ROOT.TStyle.kUSLetter 
        else               :  conf [ 'PaperSize' ] = get_int ( config , 'PaperSize' , ROOT.TStyle.kA4 )

    conf [ 'ScreenFactor'      ] = get_float ( config , 'ScreenFactor'        , 1.0  )
    
    ## conf [ 'ShowEditor'        ] = get_int   ( config , 'ShowEditor'          , 0      ) 
    ## conf [ 'ShowEventStatus'   ] = get_int   ( config , 'ShowEventStatus'     , 0      ) 
    ## conf [ 'ShowToolBar'       ] = get_int   ( config , 'ShowToolBar'         , 0      ) 

    conf [ 'StatBorderSize'    ] = get_int   ( config , 'StatBorderSize'      , 0      ) 
    conf [ 'StatColor'         ] = get_int   ( config , 'StatColor'           , 0      ) 
    conf [ 'StatFont'          ] = get_int   ( config , 'StatFont'            , font   ) 
    conf [ 'StatFontSize'      ] = get_float ( config , 'StatFontSize'        , 0.05   * scale ) 
    conf [ 'StatFormat'        ] = get_str   ( config , 'StatFormat'          , '6.3g' ) 
    conf [ 'StatH'             ] = get_float ( config , 'StatH'               , 0.15   ) ## ??? 
    conf [ 'StatStyle'         ] = get_int   ( config , 'StatStyle'           , 1001   )
    conf [ 'StatW'             ] = get_float ( config , 'StatW'               , 0.25   ) 
    conf [ 'StatX'             ] = get_float ( config , 'StatX'               , 0.9    ) 
    conf [ 'StatY'             ] = get_float ( config , 'StatY'               , 0.9    ) 

    conf [ 'StripDecimals'     ] = get_int   ( config , 'StripDecimals'       , 1      ) 


    conf [ 'TickLength_X'      ] = get_float ( config , 'TickLength_X'        , 0.03   )
    conf [ 'TickLength_Y'      ] = get_float ( config , 'TickLength_Y'        , 0.03   )
    conf [ 'TickLength_Z'      ] = get_float ( config , 'TickLength_Z'        , 0.03   )

    conf [ 'TimeOffset'        ] = get_float ( config , 'TimeOffset'          , 788918400.0 ) ## ??

    conf [ 'TitleAlign'        ] = get_int   ( config , 'TitleAlign'          , 13     ) 
    conf [ 'TitleBorderSize'   ] = get_int   ( config , 'TitleBorderSize'     , 2      ) 
    
    conf [ 'TitleColor_X'      ] = get_int   ( config , 'TitleColor_X'        , 1      ) 
    conf [ 'TitleColor_Y'      ] = get_int   ( config , 'TitleColor_Y'        , 1      ) 
    conf [ 'TitleColor_Z'      ] = get_int   ( config , 'TitleColor_Z'        , 1      ) 

    conf [ 'TitleFillColor'    ] = get_int   ( config , 'TitleFillColor'      , 19     ) 

    conf [ 'TitleFont_X'       ] = get_int   ( config , 'TitleFont_X'         , font   ) 
    conf [ 'TitleFont_Y'       ] = get_int   ( config , 'TitleFont_Y'         , font   ) 
    conf [ 'TitleFont_Z'       ] = get_int   ( config , 'TitleFont_Z'         , font   ) 

    conf [ 'TitleFontSize'     ] = get_int   ( config , 'TitleFontSize'       , 0.0    * scale )
    
    conf [ 'TitleH'            ] = get_float ( config , 'TitleH'              , 0.0    )
    
    conf [ 'TitleOffset_X'     ] = get_float ( config , 'TitleOffset_X'       , 1.0    ) 
    conf [ 'TitleOffset_Y'     ] = get_float ( config , 'TitleOffset_Y'       , 0.0    ) ## NB!!
    conf [ 'TitleOffset_Z'     ] = get_float ( config , 'TitleOffset_Z'       , 1.0    ) 

    conf [ 'TitlePS'           ] = get_str   ( config , 'TitlePS'             , ''     ) 

    conf [ 'TitleSize_X'       ] = get_float ( config , 'TitleSize_X'         , -1.0   * scale ) 
    conf [ 'TitleSize_Y'       ] = get_float ( config , 'TitleSize_Y'         ,  0.05  * scale ) 
    conf [ 'TitleSize_Z'       ] = get_float ( config , 'TitleSize_Z'         ,  0.05  * scale ) 
    
    conf [ 'TitleStyle'        ] = get_int   ( config , 'TitleStyle'          , 1001   ) 
    conf [ 'TitleTextColor'    ] = get_int   ( config , 'TitleTextColor'      , 1      ) 
    
    conf [ 'TitleW'            ] = get_float ( config , 'TitleW'              ,  0.0   )
    conf [ 'TitleX'            ] = get_float ( config , 'TitleX'              ,  0.01  )
    conf [ 'TitleXOffset'      ] = get_float ( config , 'TitleXOffset'        ,  1.0   )
    conf [ 'TitleXSize'        ] = get_float ( config , 'TitleXSize'          , -1.0   * scale )
    conf [ 'TitleY'            ] = get_float ( config , 'TitleY'              ,  0.99  * scale )
    conf [ 'TitleYOffset'      ] = get_float ( config , 'TitleYOffset'        ,  0.0   * scale )
    conf [ 'TitleYSize'        ] = get_float ( config , 'TitleYSize'          ,  0.05  * scale )

    ##
    ## Line attributes 
    ##
    
    conf [ 'LineColor'         ] = get_int   ( config , 'LineWidth'           , 1           )
    conf [ 'LineStyle'         ] = get_int   ( config , 'LineStyle'           , 1           )
    conf [ 'LineWidth'         ] = get_int   ( config , 'LineWidth'           , line_width  )

    ##
    ## Fill attributes 
    ##
    conf [ 'FillColor'         ] = get_int   ( config , 'FillWidth'           , 19          )
    conf [ 'FillStyle'         ] = get_int   ( config , 'FillStyle'           , 1001        )
    
    ##
    ## Marker attributes
    ##
    
    conf [ 'MarkerColor'       ] = get_int   ( config , 'MarkerColor'         , 1   ) 
    conf [ 'MarkerStyle'       ] = get_int   ( config , 'MarkerStyle'         , 20  ) 
    conf [ 'MarkerSize'        ] = get_float ( config , 'MarkerSize'          , 1.0 * ( 1 + ( scale - 1 ) * 0.5 ) )

    ##
    ## Text attributes
    ##
    
    conf [ 'TextAlign'         ] = get_int   ( config , 'TextAlign'           , 11   )
    conf [ 'TextAngle'         ] = get_float ( config , 'TextAngle'           , 0.0  )
    conf [ 'TextColor'         ] = get_int   ( config , 'TextColor'           , 1    )
    conf [ 'TextFont'          ] = get_int   ( config , 'TextFont'            , font )
    conf [ 'TextSize'          ] = get_float ( config , 'TextSize'            , 0.08 * scale )

    ## maximal number of digits for the axis labels 
    conf [ 'TGaxisMaxDigits'   ] = 3

    ## create the style
    style = root_style ( name )
    if not style : 
        logger.debug ( "Create new Ostap style `%s'" % name )
        style = ROOT.TStyle ( name , description )
        
    set_style  ( style , conf , base_style = base_style ) 

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
##                                                                      The END 
# =============================================================================
