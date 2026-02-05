#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/core.py
#  Core objects for ostap 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
""" Core objects for ostap 
"""
# ============================================================================= 
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'cpp'                 , ## global C++ namespace
    'std'                 , ## C++ namespace std
    'Ostap'               , ## C++ namespace Ostap
    'rootID'              , ## global identifier for ROOT objects
    'funcID'              , ## global identifier for ROOT functions 
    'funID'               , ## global identifier for ROOT functions 
    'fID'                 , ## global identifier for ROOT functions 
    'histoID'             , ## global identifier for ROOT histograms 
    'hID'                 , ## global identifier for ROOT histograms 
    'grID'                , ## global identifier for ROOT graphs 
    'dsID'                , ## global identifier for ROOT/RooFit datasets
    'usedRootID'          , ## check if the  name/identified already used by ROOT/RooFit  
    ##
    'VE'                  , ## shortcut for Gaudi::Math::ValuewithError
    'SE'                  , ## shortcut for StatEntity
    'WSE'                 , ## shortcut for Gaudi::Math::WStatEntity 
    ##
    'binomEff'            , ## binomial efficiency  
    'binomEff2'           , ## binomial efficiency
    'zechEff'             , ## binomial efficiency: Zech's recipe 
    'wilsonEff'           , ## binomial efficiency: Wilson 
    'agrestiCoullEff'     , ## binomial efficiency: Agresti-Coull
    ##
    'iszero'              , ## comparison with zero  for doubles  
    'isequal'             , ## comparison for doubles 
    'isint'               , ## Is float value actually int  ? 
    'islong'              , ## Is float value actually long ?
    'inrange'             , ## Is float walue in range ?  
    ##
    'natural_entry'       , ## natural entry?  @see Ostap::Math::natural_entry 
    'natural_number'      , ## natural nunber? @see Ostap::Math::natural_number
    ##
    'valid_pointer'       , ## Is it a valid C++ pointer?
    'root_enum'           , ## Get enum from ROOT by name 
    ##
    'strings'             , ## construct std::vector<std::string>
    ##
    'StatusCode'          , ## status code
    'SUCCESS'             , ## status code SUCCESS 
    'FAILURE'             , ## status code FAILURE
    ##
    'is_sorted'           , ## check that list is sorted
    ## 
    'RooSilent'           , ## control RooFit verbosity
    'ROOTIgnore'          , ## control ROOT verbosity, suppress ROOT errors
    'rooSilent'           , ## control RooFit verbosity
    'roo_silent'          , ## control RooFit verbosity 
    'rootError'           , ## control ROOT verbosity 
    'rootWarning'         , ## control ROOT verbosity
    ## 
    'rootException'       , ## context manager to perform ROOT Error -> C++/Python exception
    'RootError2Exception' , ## context manager to perform ROOT Error -> C++/Python exception
    ##
)
# =============================================================================
from   ostap.core.meta_info   import root_info
from   ostap.core.base        import ( ROOTIgnore          ,
                                       RooSilent           ,
                                       rootException       ,
                                       RootError2Exception ,
                                       valid_pointer       ,
                                       rooSilent           ,
                                       roo_silent          ,                                        
                                       rootError           ,
                                       rootWarning         ) 
from   ostap.math.base        import ( Ostap    , std      , cpp ,
                                       iszero   , isequal  ,
                                       isint    , islong   ,
                                       inrange  , strings  , 
                                       natural_number      ,
                                       natural_entry       )
from   ostap.math.ve          import VE
from   ostap.stats.counters   import SE , WSE 
from   ostap.core.ostap_types import integer_types, sequence_types, string_types
from   ostap.utils.basic      import NoContext, loop_items, typename
import ostap.core.ostap_setup 
import ostap.plotting.color   
import ROOT, cppyy, math, sys, os, re  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.core.core' )
else                       : logger = getLogger( __name__     )
# =============================================================================
logger.debug ( 'Core objects/classes/functions for Ostap')
# =============================================================================
binomEff        = Ostap.Math.binomEff
binomEff2       = Ostap.Math.binomEff2
zechEff         = Ostap.Math.zechEff
wilsonEff       = Ostap.Math.wilsonEff
agrestiCoullEff = Ostap.Math.agrestiCoullEff
# =============================================================================
## global identifier for ROOT objects 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @see Ostap::Utils::rootID 
#  @date   2011-06-07
def rootID ( prefix = 'root_' , suffix = '' ) :
    """ Construct the unique ROOT-id 
    """
    if prefix : prefix = prefix.replace ( ' ' , '_' )
    else      : prefix = 'root_'
    if suffix : suffix = suffix.replace ( ' ' , '_' )
    ## 
    if            not prefix.endswith   ( '_' ) : prefix = prefix + '_'
    if suffix and not suffix.startswith ( '_' ) : suffix = '_' + suffix 
    ## 
    return  Ostap.Utils.rootID ( prefix , suffix )

# =================================================================================
## used by ROOT/RooFit ?  
usedRootID = Ostap.Utils.usedRootID 

# =============================================================================
## global ROOT identified for function objects 
def funcID  ( prefix = 'f_'  , suffix = '' ) : return rootID  ( prefix , suffix )
## global ROOT identified for function objects 
def funID   ( prefix = 'f_'  , suffix = '' ) : return funcID  ( prefix , suffix )
## global ROOT identified for function objects 
def fID     ( prefix = 'f_'  , suffix = '' ) : return funcID  ( prefix , suffix )
## global ROOT identified for histogram objects 
def histoID ( prefix = 'h_'  , suffix = '' ) : return rootID  ( prefix , suffix )
## global ROOT identified for histogram objects 
def histID  ( prefix = 'h_'  , suffix = '' ) : return histoID ( prefix , suffix )
## global ROOT identified for histogram objects 
def hID     ( prefix = 'h_'  , suffix = '' ) : return histoID ( prefix , suffix )
## global ROOT identified for dataset objects 
def dsID    ( prefix = 'ds_' , suffix = '' ) : return rootID  ( prefix , suffix )
## global ROOT identified for graphs objects 
def grID    ( prefix = 'gr_' , suffix = '' ) : return rootID  ( prefix , suffix )

# ==================================================================================
## get current directory in ROOT
#  @code
#  d = cwd()
#  print d
#  @endcode 
def cwd() :
    """ Get current directory in ROOT
    >>> d = cdw() 
    """
    groot = ROOT.ROOT.GetROOT ()
    return groot.CurrentDirectory()

# =================================== ===============================================
## get current directory in ROOT
#  @code
#  print pwd() 
#  @endcode 
def pwd() :
    """ Get current directory in ROOT
    >>> print pwd() 
    """
    groot = ROOT.ROOT.GetROOT ()
    return groot.CurrentDirectory().GetPath() 

# =============================================================================
_FAILURE = Ostap.StatusCode.FAILURE 
## printout of status code
def _sc_print_ ( sc ) :
    """ Print the Status Code
    >>> st = ...
    >>> print st
    """
    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = list ( range ( 8 ) )
    ##
    from ostap.logger.colorized import colored_string
    if   sc.isSuccess     () : return colored_string ( 'SUCCESS'     , WHITE , GREEN  , True ) 
    elif sc.isRecoverable () : return colored_string ( 'RECOVERABLE' , RED   , YELLOW , True ) 
    elif _FAILURE != sc.code  () :
        return colored_string ( 'FAILURE[%d]' % sc.code() , YELLOW , RED   , True )
    return colored_string ( 'FAILURE' , YELLOW , RED , True ) 
        
StatusCode = Ostap.StatusCode 
StatusCode .__repr__ = _sc_print_
StatusCode .__str__  = _sc_print_

SUCCESS     = StatusCode ( Ostap.StatusCode.SUCCESS     )
FAILURE     = StatusCode ( Ostap.StatusCode.FAILURE     )

# =============================================================================
## Get value of enum form ROOT by name
#  @code
#  red = root_enum ( 'Red' ) 
#  @endcode 
def root_enum ( name , default = None ) :
    """ Get value of enum from ROOT by name
    >>> red   = root_enum ( 'Red'   )
    >>> blue  = root_enum ( 'Blue'  )    
    >>> green = root_enum ( 'Green' )    
    """    
    assert isinstance ( name , string_types ) , \
           'Invalid type for enum name %s' % type ( name )
           
    return getattr ( ROOT, 'k' + name , default )

# =============================================================================
# predefined ROOT colors
# @code
#  enum EColor { kWhite =0,   kBlack =1,   kGray=920,
#                kRed   =632, kGreen =416, kBlue=600, kYellow=400, kMagenta=616, kCyan=432,
#                kOrange=800, kSpring=820, kTeal=840, kAzure =860, kViolet =880, kPink=900 };
# @endcode 
__root_colors = {
    'white'   : 0   ,
    'black'   : 1   ,
    'gray'    : 920 ,
    'red'     : 632 ,
    'green'   : 416 ,
    'blue'    : 600 ,
    'yellow'  : 400 ,
    'magenta' : 616 ,
    'cyan'    : 432 ,
    'orange'  : 800 ,
    'spring'  : 820 ,
    'teal'    : 840 ,
    'azure'   : 860 ,
    'violet'  : 880 ,
    'pink'    : 900 
    }
# =============================================================================
## get predefiend ROOT color by name
#  @code
#  c1 = root_color ( 'Red'  )
#  c1 = root_color ( 'kRed' )
#  c1 = root_color ( 'RED'  )
#  c1 = root_color ( 'KRED' )
#  @endcode
def root_colors ( color ) :
    """ Get predefiend ROOT color by name
    >>> c1 = root_color ( 'Red'  )
    >>> c1 = root_color ( 'kRed' )
    >>> c1 = root_color ( 'RED'  )
    >>> c1 = root_color ( 'KRED' )
    """
    if not isinstance ( color , string_types ) : return color 
    c = color.lower().replace('_','')
    if c and c[0] == 'k' : c = c[1:]
    return __root_colors.get ( c , color )
    
# =============================================================================
## Convert color name into color index
def check_color ( color ) :
    """ Convert color name into color index 
    """
    if isinstance ( color , string_types ) :
        
        clow = color.lower() 
        if   clow in ( 'black'   ,     ) : return 1 
        elif clow in ( 'red'     , 'r' ) : return 2 
        elif clow in ( 'blue'    , 'b' ) : return 4 
        elif clow in ( 'green'   , 'g' ) : return 8 ## dark green here 
        elif clow in ( 'yellow'  , 'y' ) : return 5 
        elif clow in ( 'cyan'    , 'c' ) : return 6 
        elif clow in ( 'magenta' , 'm' ) : return 7                

        ## check predefined ROOT color name
        c = root_colors ( color)
        if isinstance   ( c , integer_type ) and 0 < c : return int ( c ) 
        
        ## get the color from ROOT by color name
        c = root_enum   ( clow.capitalize() , None )
        if c is None    : logger.error ('Unknown color:"%s"' % color ) 
        elif isinstance ( c , integer_types ) and 0 < c :
            return int ( c )                               ## RETURN 

    return color

# =============================================================================
## alowed arguments for draw-function 
draw_args = frozenset ( [ 'linecolor'   , 'markercolor' , 'fillcolor'   , 'color' , 
                          'linestyle'   , 'linewidth'   , 'width'       ,
                          'markerstyle' , 'markersize'  , 'marker'      ,
                          'fillstyle'   , 'fill'        , 'opacity'     , 'opaque' ,
                          'minimum'     , 'minimal'     , 'min'         , 'minval' , 'minvalue' ,
                          'maximum'     , 'maximal'     , 'max'         , 'maxval' , 'maxvalue' ,
                          'labelsize'   , 'labelfont'   , 'labelscale'  ,
                          'xaxislabeloffset' ,
                          'yaxislabeloffset' ,
                          'zaxislabeloffset' ,
                          'copy'  ,
                          'logx'  ,
                          'logy'  , 
                          'gridx' ,
                          'gridy' ] )
## ============================================================================
## remove "draw-args" from dictionary of arguments 
def remove_draw_args ( kwargs  ) :
    """ remove "draw-args" from dictionary of arguments """
    nargs = {}
    from ostap.utils.cidict import cidict_fun 
    for k, v in loop_items ( kwargs ) :
        key = cidict_fun ( k )
        if key in draw_args : continue
        nargs [ k ] = v
    return nargs 
# =============================================================================
## Silent draw 
# =============================================================================
if not hasattr ( ROOT.TObject , 'draw' ) :

    ## save old method
    ROOT.TObject._old_draw_ = ROOT.TObject.Draw
    ##  new draw method: silent draw
    def _TO_draw_ ( obj , option = '' , *options , **kwargs ) :
        """ (silent) Draw of ROOT object, optionally setting the drawing attributes when applicable:
        
        - LineColor
        - LineStyle
        - LineWidth
        - MarkerColor
        - MarkerStyle
        - MarkerSize
        - FillColor
        - FillStyle
        - Min/Minimal
        - Max/Maximal

        - see ROOT.TAttLine, ROOT.TAttMarker and ROOT.TAttFill

        (case-insensitive and underscore-blind)
        
        >>> obj
        >>> obj.draw ()  ##
        >>> obj.draw ( linecolor   = 2 ) 
        >>> obj.draw ( LineColor   = 2 , fill_style = 1003 ) 
        >>> obj.draw ( markerStyle = 22  )

        >>> obj.draw ( minimal     = -1  )
        >>> obj.draw ( max         = 100 )
        
        - Combined FillStyle/FillColor options
        
        >> obj.draw ( fill = 29   ) ## Color if    0 < fill < 1000 
        >> obj.draw ( fill = 2045 ) ## Style if 1000 < fill 

        - LogX/logY options :

        >>> obj.draw ( logx = True )
        >>> obj.draw ( logy = True )

        - Opacity/Opaque for TAttFill 

        >>> obg.draw ( opacity = 0.45 ) 
        >>> obg.draw ( opaque  = 0.45 ) 
        
        `args` are treated as "additional" options if they are of string types 
        otherwise are ignored and warning message is issued  

        """
        if   isinstance ( option , string_types   ) : pass
        elif isinstance ( option , sequence_types ) and all ( isinstance ( o , string_types ) for o in option ) :
            option = ' '.join ( o for o in option )
        else :
            logger.warning ( "draw: Invalid type of `option`: %s " % typename ( option ) ) 
            
        ## treat positional arguments as addtional options :
        if   options and all ( isinstance ( o , string_types ) for o in options ) :
            option  = option + ( ' '.join ( o for o in options ) )
            options = ()
        elif options :
            logger.warning ( "draw: Invalid type of `options` : (%s)" % ( ','.join ( typename ( o ) for o in options ) ) ) 

            
        from   ostap.utils.cidict import cidict, cidict_fun
        import ostap.plotting.style 
        kw = cidict ( transform = cidict_fun , **kwargs )

        ## Global color (Line, Marker, Fill) 
        if 'Color' in kw :
            color = check_color ( kw.pop('Color' ) ) 
            if not 'LineColor'   in kw and hasattr ( obj , 'SetLineColor'   ) :
                obj.SetLineColor   ( color )
            if not 'MarkerColor' in kw and hasattr ( obj , 'SetMarkerColor' ) :
                obj.SetMarkerColor ( color )
            if not 'FillColor'   in kw and hasattr ( obj , 'SetFillColor'   ) :
                obj.SetFillColor   ( color )
                
        ## Line
        if isinstance ( obj , ROOT.TAttLine ) :
            if   'LineColor' in kw : obj.SetLineColor ( check_color ( kw.pop('LineColor' ) ) ) 
            if   'LineStyle' in kw : obj.SetLineStyle ( kw.pop('LineStyle' ) )
            if   'LineWidth' in kw : obj.SetLineWidth ( kw.pop('LineWidth' ) ) 
            elif 'Width'     in kw : obj.SetLineWidth ( kw.pop('Width'     ) ) 
                                                    
        ## Marker
        if isinstance ( obj , ROOT.TAttMarker ) :               
            if   'MarkerColor' in kw : obj.SetMarkerColor ( check_color ( kw.pop('MarkerColor' ) ) ) 
            if   'MarkerSize'  in kw : obj.SetMarkerSize  ( kw.pop('MarkerSize'  ) )
            if   'MarkerStyle' in kw : obj.SetMarkerStyle ( kw.pop('MarkerStyle' ) )
            elif 'Marker'      in kw : obj.SetMarkerStyle ( kw.pop('Marker'      ) )

        ## Area
        if isinstance ( obj , ROOT.TAttFill  ) :               
            fcolor = False 
            fstyle = False
            if 'FillColor'   in kw :                
                fcolor = True                            
                obj.SetFillColor   ( check_color ( kw.pop('FillColor' ) ) ) 
            if 'FillStyle'   in kw :
                fstyle  = True  
                obj.SetFillStyle   ( kw.pop('FillStyle' ) )
                
            if 'Fill' in kw :

                fill = kw.pop ( 'Fill' )
            
                if   fill is True  and ( not fstyle ) and has_style : obj.SetFillStyle ( 1001 )
                elif fill is False and ( not fstyle ) and has_style : obj.SetFillStyle (    0 )
                elif isinstance ( fill , integer_types ) :
                    if      0 < fill < 1000 and not fcolor :
                        obj.SetFillColor ( fill )
                        if not fstyle : obj.SetFillStyle ( 1001 )
                    elif 1000 < fill        and not fstyle :                        
                        obj.SetFillStyle ( fill  )
            
            if ( ( 'Opacity' in kw ) or ( 'Opaque' in kw ) ) and 1001 == obj.GetFillStyle () :
                fc = obj.GetFillColor()
                if fc :
                    if   'Opacity' in kw : obj.SetFillColorAlpha ( fc , kw.pop ('Opacity' ) ) 
                    elif 'Opaque'  in kw : obj.SetFillColorAlpha ( fc , kw.pop ('Opaque'  ) ) 
                    
        ## Min/max values  

        if   'Minimum'     in kw and hasattr ( obj , 'SetMinimum' ) :
            obj.SetMinimum     ( kw.pop ( 'Minimum'  ) )
        elif 'Minimal'     in kw and hasattr ( obj , 'SetMinimum' ) :
            obj.SetMinimum     ( kw.pop ( 'Minimal'  ) )
        elif 'Min'         in kw and hasattr ( obj , 'SetMinimum' ) :
            obj.SetMinimum     ( kw.pop ( 'Min'      ) )
        elif 'MinVal'      in kw and hasattr ( obj , 'SetMinimum' ) :
            obj.SetMinimum     ( kw.pop ( 'MinVal'   ) )
        elif 'MinValue'    in kw and hasattr ( obj , 'SetMinimum' ) :
            obj.SetMinimum     ( kw.pop ( 'MinValue' ) )
            
        if   'Maximum'     in kw and hasattr ( obj , 'SetMaximum' ) :
            obj.SetMaximum     ( kw.pop ( 'Maximum'  ) )
        elif 'Maximal'     in kw and hasattr ( obj , 'SetMaximum' ) :
            obj.SetMaximum     ( kw.pop ( 'Maximal'  ) )
        elif 'Max'         in kw and hasattr ( obj , 'SetMaximum' ) :
            obj.SetMaximum     ( kw.pop ( 'Max'      ) )
        elif 'MaxVal'      in kw and hasattr ( obj , 'SetMaximum' ) :
            obj.SetMaximum     ( kw.pop ( 'MaxVal'   ) )
        elif 'MaxValue'    in kw and hasattr ( obj , 'SetMaximum' ) :
            obj.SetMaximum     ( kw.pop ( 'MaxValue' ) )

        if 'LabelSize' in kw or 'LabelFont' in kw or 'LabelScale' in kw  :

            axis = [] 
            if hasattr ( obj , 'GetXaxis'     ) :
                xa = obj.GetXaxis()
                if xa : axis.append ( xa )
            if hasattr ( obj , 'GetYaxis'     ) :
                ya = obj.GetYaxis()
                if ya : axis.append ( ya ) 
            if hasattr ( obj , 'GetZaxis'     ) :
                za = obj.GetZaxis()
                if za : axis.append ( za ) 

            if axis and 'LabelSize' in kw :
                ls = kw.pop ( 'LabelSize' ) 
                for a in axis : a.SetLabelSize ( ls ) 

            if axis and 'LabelFont' in kw :
                lf = kw.pop ( 'LabelFont' ) 
                for a in axis : a.SetLabelFont ( lf )                

            if axis and 'LabelScale' in kw :
                ls = kw.pop ( 'LabelScale' ) 
                for a in axis : a.SetLabelSize  ( ls * a.GetLabelSize () ) 

        if 'XaxisLabelOffset' in kw and hasattr ( obj , 'GetXaxis' ) :
            xa = obj.GetXaxis() 
            if xa : xa.SetLabelOffset ( kw.pop  ( 'XaxisLabelOffset' ) ) 

        if 'YaxisLabelOffset' in kw and hasattr ( obj , 'GetYaxis' ) :
            ya = obj.GetYaxis() 
            if ya : ya.SetLabelOffset ( kw.pop  ( 'YaxisLabelOffset' ) ) 
                                        
        if 'ZaxisLabelOffset' in kw and hasattr ( obj , 'GetZaxis' ) :
            za = obj.GetZaxis() 
            if za : za.SetLabelOffset ( kw.pop  ( 'ZaxisLabelOffset' ) ) 
                                        
        # =====================================================================
        ## finally, draw it!        
        with rootWarning() , rooSilent ( 2 )  :
            
            copy = kw.pop ( 'copy' , False ) and hasattr ( obj , 'DrawCopy' )

            if copy : result = obj.DrawCopy ( option ) ## *args )
            else    : result = obj.Draw     ( option ) ## *args )
                
            result = obj

        # =====================================================================
        ## update the pad/canvas
        # =====================================================================
        try : # ===============================================================
            # =================================================================
            pad = Ostap.Utils.get_pad() 
            if pad : Ostap.Utils.pad_update  ( pad )\
            # =================================================================
        except: # =============================================================
            # =================================================================
            logger.error ( 'Exception from get_pad/pad_update!' , exc_info = True )
            
        # =====================================================================
        ## Lin/log scale ?
        if 'LogY' in kw or 'LogX' in kw :            
            pad = Ostap.Utils.get_pad() 
            if pad and 'LogX' in kw : pad.SetLogx ( kw.pop ( 'LogX' ) )
            if pad and 'LogY' in kw : pad.SetLogy ( kw.pop ( 'LogY' ) )
            
        # =====================================================================
        ## GridX/GridY  ?
        if 'GridY' in kw or 'GridX' in kw :            
            pad = Ostap.Utils.get_pad() 
            if pad and 'GridX' in kw : pad.SetGridx ( kw.pop ( 'GridX' ) )
            if pad and 'GridY' in kw : pad.SetGridy ( kw.pop ( 'GridY' ) )

        # =====================================================================
        # =====================================================================
        if kw or options : # =====================================================
            # =================================================================
            from ostap.logger.utils import print_args
            title  = kw.pop ( 'title'  , '' )
            prefix = kw.pop ( 'prefix' , '' )
            if   kw and options : title  = 'draw: Unused %d+%d arguments' % ( len ( options ) , len ( kw ) )
            elif kw             : title  = 'draw: Unused %d arguments'    %   len ( kw      )
            elif args           : title  = 'draw: Unused %d arguments'    %   len ( options )            
            prefix = '# '
            table  = print_args ( *options , title  = title , prefix = prefix , **kw )
            logger.warning ( '%s:\n%s' % ( title , table ) )
            
        return result 

    ROOT.TObject.draw       = _TO_draw_
    ROOT.TObject.draw_ostap = _TO_draw_

# =============================================================================
## Set/Set name/title 
# =============================================================================
def _tn_name_get_  ( self )         : return self.GetName()
def _tn_name_set_  ( self , value ) : self.SetName( value )
_tn_name_doc_  = "`name' of the object using GetName/SetName"
def _tn_title_get_ ( self )         : return self.GetTitle()
def _tn_title_set_ ( self , value ) : self.SetTitle( value )
_tn_title_doc_ = "`title' of the object using GetTitle/SetTitle"
ROOT.TNamed.name  = property ( _tn_name_get_  ,  _tn_name_set_  , None , _tn_name_doc_  ) 
ROOT.TNamed.title = property ( _tn_title_get_ ,  _tn_title_set_ , None , _tn_title_doc_ ) 

# =============================================================================
## Get the full path of the named object
#  @code
#  obj  = ...
#  path = obj.path
#  path = obj.fullpath
#  path = obj.full_path
#  @endcode 
def tnamed_path ( obj ) :
    """ Get the full path of the named object
    >>> obj  = ...
    >>> path = obj.path
    >>> path = obj.full_path
    >>> path = obj.fullpath
    """
    if not valid_pointer ( obj ) : return '<INVALID-OBJECT>'
    if isinstance ( obj , ROOT.TDirectory ) :
        path = obj.GetPath ()
        _ , _ , p = path.rpartition(':/')
        if p : return  '/'.join ( [ p , obj.GetName () ] )
        return obj.GetName() 

    ## 
    if not hasattr ( obj ,  'GetDirectory' ) : return obj.GetName()
    
    rdir = obj.GetDirectory() 
    if not valid_pointer ( rdir ) : return obj.GetName()
        
    path = rdir.GetPath() 
    h , s , p = path.rpartition(':/')
    if p : return  '/'.join ( [ p , obj.GetName () ] )
    return obj.GetName() 

ROOT.TNamed.path      = property ( tnamed_path , None , None , tnamed_path.__doc__ ) 
ROOT.TNamed.fullpath  = property ( tnamed_path , None , None , tnamed_path.__doc__ )  
ROOT.TNamed.full_path = property ( tnamed_path , None , None , tnamed_path.__doc__ ) 

# =============================================================================
## check that list is sorted 
def is_sorted ( lst ) :
    """ Check that list is sorted
    """
    l = len ( lst )  
    return all ( lst [ i ] <= lst [ i + 1 ] for i in range ( l - 1 ) ) if lst else True 

# =============================================================================
## "Contains" method for <code>TCollection</code> object
#   @code
#   collection = ... 
#   'name' in collection
#   obj    in collection 
#   @endcode 
#   @see TCollection
#   @see TCollection::FindObject
def _rtc_contains_ ( lst , item ) :
    """ Contains method for `TCollection` object
    >>> collection = ...
    >>> 'name' in collection
    >>> obj    in collection 
    - see TCollection
    - see TCollection.FindObject
    """
    return lst.FindObject ( item ) 

# =============================================================================
## "Contains" method for <code>TSeqCollection</code> object
#   @code
#   lst = ... 
#   'name' in lst 
#   obj    in lst
#   5      in lst 
#   @endcode 
#   @see TSeqCollection
#   @see TSeqCollection::FindObject
def _rtl_contains_ ( lst , item ) :
    """ Contains method for `TSeqCollection` object
    >>> lst  = ...
    >>> 'name' in lst 
    >>> obj    in lst
    >>> 5 in lst 
    - see TSeqCollection
    - see TSeqCollection.FindObject
    """
    if isinstance ( item , integer_types ) :
        return 0 <= item < len ( lst ) 
        
    return _rtc_contains_ ( lst , item ) 

ROOT.TCollection   . __contains__ = _rtc_contains_
ROOT.TSeqCollection. __contains__ = _rtl_contains_



# =============================================================================
## "Get-item" method for <code>TCollection</code> object
#   @code
#   collection = ... 
#   collection['name'] 
#   @endcode 
#   @see TCollection
#   @see TCollection::FindObject
def _rtc_getitem_ ( lst , item ) :
    """``Get-item'' method for `TCollection` object
    >>> collection = ...
    >>> collection['name'] 
    - see TCollection
    - see TCollection.FindObject
    """
    
    obj = lst.FindObject ( item )
    if not valid_pointer ( obj ) :
        raise KeyError("No such key is found: %s" % str ( item ) )
    return obj 
    
# =============================================================================
## "Get-item " method for <code>TSeqCollection</code> object
#   @code
#   lst = ... 
#   lst['name']
#   lst[ 5    ]
#   lst[ -3   ]
#   lst[ 5:10 ]
#   lst[ 5:   ]
#   @endcode 
#   @see TSeqCollection
#   @see TSeqCollection::FindObject
#   @see TSeqCollection::At
def _rtl_getitem_ ( lst , item ) :
    """``Get-item'' method for `TSeqCollection` object
    >>> lst  = ...
    >>> lst['lst'] 
    >>> lst[ 5    ]
    >>> lst[ -3   ]
    >>> lst[ 5:10 ]
    >>> lst[ 5:   ]
    - see TSeqCollection
    - see TSeqCollection.FindObject
    - see TSeqCollection.At
    """

    if isinstance ( item , slice ) :

        result  = type ( lst ) ()
        for index in item.indices ( len ( lst ) ) :
            result.Add ( lst.At ( index ) )

        return result

    elif isinstance ( item , integer_types ) :
        
        len_lst = len ( lst ) 
        if item < 0 : item += len ( lst ) ## allow negative indices! 
        
        if not 0 <= item < len_lst :
            raise KeyError("No such key is found: %s" % item )

        return lst.At ( item )
            
    return _rtc_getitem_ ( lst , item ) 

ROOT.TCollection    . __getitem__ = _rtc_getitem_
ROOT.TSeqCollection . __getitem__ = _rtl_getitem_

# =============================================================================
## "Get" method for <code>TCollection</code> object
#   @code
#   collection = ... 
#   collection.get('name') 
#   @endcode 
#   @see TCollection
#   @see TCollection::FindObject
def _rtc_get_ ( lst , item , default = None ) :
    """``Get'' method for `TCollection` object
    >>> collection = ...
    >>> collection.get( 'name' , defaultvalue ) 
    - see TCollection
    - see TCollection.FindObject
    """
    obj = lst.FindObject ( item )
    return obj if valid_pointer ( obj ) else default 

# =============================================================================
## "Get" method for <code>TSeqCollection</code> object
#   @code
#   lst = ... 
#   lst.get('name') 
#   @endcode 
#   @see TSeqCollection
#   @see TSeqCollection::FindObject
#   @see TSeqCollection::At 
def _rtl_get_ ( lst , item , default = None ) :
    """``Get'' method for `TSeqCollection` object
    >>> lst = ...
    >>> lst.get( 'name' , defaultvalue ) 
    - see TSeqCollection
    - see TSeqCollection.FindObject
    - see TSeqCollection.At 
    """
    if isinstance ( item , integer_types ) :
        
        len_lst = len ( lst )
        if item < 0 : item += len_lst  ## allow negative indices
        
        if not 0 <= item < len_lst  : return default

        return lst.At ( item ) 
    
    return _rtc_get_ ( lst , item , default ) 

ROOT.TCollection    . get = _rtc_get_
ROOT.TSeqCollection . get = _rtl_get_

# =============================================================================
## Add/append element (or elements) to <code>TCollection</code> container
#  @code
#  collectiom   = ...
#  item         = ...
#  collection  += item
#  items        = ...
#  collection  += items
#  collection  += (item1, item2, item3) 
#  @endcode 
#  @see TCollection::Add 
#  @see TCollection::AddAll 
def _rtc_iadd_ ( self , item ) :
    """ Add/append element or elements to `ROOT.TCollection` container
    >>> collectiom   = ...
    >>> item         = ...
    >>> collection  += item
    >>> items        = ...
    >>> collection  += items
    >>> collection  += (item1, item2, item3)
    - see ROOT.TCollection.Add 
    - see ROOT.TCollection.AddAll 
    """
    
    if   isinstance ( item , ROOT.TCollection ) : self.AddAll ( item )
    elif isinstance ( item , ROOT.TObject     ) : self.Add    ( item ) 
    elif isinstance ( item , sequence_types   ) :
        for i in item : self += i
    else : 
        return NotImplemented

    return self 


ROOT.TCollection. __iadd__ = _rtc_iadd_ 
    


# =============================================================================
## define the build/cache&tmp directories for ROOT&Ostap 


# ==============================================================================
## make the default style 
import ostap.plotting.style 
# =============================================================================
_decorated_classes_ = (
    ROOT.TObject        ,
    ROOT.TNamed         ,
    ROOT.TDirectory     ,
    ROOT.TCollection    , 
    ROOT.TSeqCollection , 
    )

_new_methods_       = (
    #
    ROOT.TObject.draw  ,
    ROOT.TNamed .name  ,
    ROOT.TNamed .title ,
    ROOT.TNamed .path  ,
    #
    ROOT.TCollection    . __contains__ ,
    ROOT.TCollection    . __getitem__  ,
    ROOT.TCollection    . get          ,
    ROOT.TSeqCollection . __contains__ ,
    ROOT.TSeqCollection . __getitem__  ,
    ROOT.TSeqCollection . get          ,
    #
    ROOT.TCollection    . __iadd__     ,
    )

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
