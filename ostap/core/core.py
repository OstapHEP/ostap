#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/core.py
#  Core objects for ostap 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Core objects for ostap 
"""
# ============================================================================= 
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'cpp'                 , ## global C++ namespace
    'std'                 , ## C++ namespace std
    'Ostap'               , ## C++ namespace Ostap
    'ROOTCWD'             , ## context manager to keep/preserve ROOT current directory
    'rootID'              , ## global identifier for ROOT objects
    'funcID'              , ## global identifier for ROOT functions 
    'funID'               , ## global identifier for ROOT functions 
    'fID'                 , ## global identifier for ROOT functions 
    'histoID'             , ## global identifier for ROOT histograms 
    'hID'                 , ## global identifier for ROOT histograms 
    'grID'                , ## global identifier for ROOT graphs 
    'dsID'                , ## global identifier for ROOT/RooFit datasets
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
    'in_test'             , ## Are we in CMAKE-test regime?
    )
# =============================================================================
from   ostap.core.meta_info   import root_info
from   ostap.math.base        import ( Ostap    , std     , cpp ,  
                                       iszero   , isequal ,
                                       isint    , islong  ,
                                       inrange  , strings , 
                                       natural_number     ,
                                       natural_entry      ,
                                       ROOTIgnore         )
from   ostap.math.ve          import VE
from   ostap.stats.counters   import SE , WSE 
from   ostap.core.ostap_types import integer_types, sequence_types, string_types
from   ostap.utils.basic      import NoContext, loop_items         
import ROOT, cppyy, math, sys, os, re  
# =============================================================================
## ROOT.ROOT.EnableThreadSafety()
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
## @class ROOTCWD
#  context manager to preserve current directory (rather confusing stuff in ROOT)
#  @code
#  groot = ROOT.ROOT.GetROOT()
#  print groot.CurrentDirectory() 
#  with ROOTCWD() :
#     print groot.CurrentDirectory() 
#     rfile = ROOT.TFile( 'test.root' , 'recreate' )
#     print groot.CurrentDirectory() 
#  print groot.CurrentDirectory() 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
if root_info < ( 6, 29 ) : # ==================================================
    # =========================================================================
    class ROOTCWD(object) :
        """ Context manager to preserve current directory
        (rather confusing stuff in ROOT) 
        >>> print the_ROOT.CurrentDirectory() 
        >>> with ROOTCWD() :
        ...     print the_ROOT.CurrentDirectory() 
        ...     rfile = ROOT.TFile( 'test.root' , 'recreate' )
        ...     print the_ROOT.CurrentDirectory() 
        ... print the_ROOT.CurrentDirectory() 
        """
        def __init__ ( self ) :
            self._dir = None
            
        def __del__ ( self ) :
            self._dir = None 
            del self._dir
            
        ## context manager ENTER 
        def __enter__ ( self ) :
            "Save current working directory"
            self._dir = None
            
            ## ROOT::TDirectory::TContext appears in ROOT 6/23/01
            groot     = ROOT.ROOT.GetROOT ()
            if groot :
                cwd = groot.CurrentDirectory()
                cwd = cwd.load() ## resolve std::atomic 
                if cwd : self._dir = cwd
                
            return self

        ## context manager EXIT 
        def __exit__  ( self , *_ ) :
            "Make the previous directory current again"
            
            if self._dir :
                
                odir = self._dir
                
                self._dir = None
                
                fdir = odir.GetFile () if isinstance ( odir , ROOT.TDirectoryFile ) else None
                
                if fdir and not fdir.IsOpen () :
                    
                    groot = ROOT.ROOT.GetROOT ()
                    groot.cd ()
                    
                else :
                    
                    odir.cd()

            self._dir = None 
            del self._dir
    # ========================================================================
elif root_info < ( 6 , 37 ) : # ==============================================
    # ========================================================================
    class ROOTCWD(object) :
        """ Context manager to preserve current directory
        (rather confusing stuff in ROOT) 
        >>> print the_ROOT.CurrentDirectory() 
        >>> with ROOTCWD() :
        ...     print the_ROOT.CurrentDirectory() 
        ...     rfile = ROOT.TFile( 'test.root' , 'recreate' )
        ...     print the_ROOT.CurrentDirectory() 
        ... print the_ROOT.CurrentDirectory() 
        """

        def __enter__ ( self ) :            
            self._cntx = ROOT.TDirectory.TContext()
            self._cntx.__enter__()

        def __exit__ ( self , *_ ) :
            result = self._cntx.__exit__ ( *_ )
            groot     = ROOT.ROOT.GetROOT ()
            if groot :
                cwd = groot.CurrentDirectory().load()
                if valid_pointer ( cwd ) :
                    if isinstance ( cwd , ROOT.TDirectoryFile ) :
                        fdir = cwd.GetFile ()
                        if valid_pointer ( fdir ) and not fdir.IsOpen() :
                            groot = ROOT.ROOT.GetROOT ()
                            groot.cd ()
                            
            return result

    # ========================================================================
else : # =================================================
    # ========================================================================
    class ROOTCWD(object) :
        """ Context manager to preserve current directory
        (rather confusing stuff in ROOT) 
        >>> print the_ROOT.CurrentDirectory() 
        >>> with ROOTCWD() :
        ...     print the_ROOT.CurrentDirectory() 
        ...     rfile = ROOT.TFile( 'test.root' , 'recreate' )
        ...     print the_ROOT.CurrentDirectory() 
        ... print the_ROOT.CurrentDirectory() 
        """

        def __enter__ ( self ) :            
            self._cntx = ROOT.TDirectory.TContext()
            self._cntx.__enter__()

        def __exit__ ( self , *_ ) :
            result = self._cntx.__exit__ ( *_ )
            ## recheck the current directory 
            cwd = ROOT.TDirectory.CurrentDirectory().load()
            if valid_pointer ( cwd ) and isinstance ( cwd , ROOT.TDirectoryFile ) :
                fdir = cwd.GetFile ()
                if valid_pointer ( fdir ) and not fdir.IsOpen() :
                    groot = ROOT.ROOT.GetROOT ()
                    groot.cd ()                            
            return result
    
    # ========================================================================

# =============================================================================
## global identifier for ROOT objects 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @see Ostap::Utils::rootID 
#  @date   2011-06-07
def rootID ( prefix = 'o_' ) :
    """ Construct the unique ROOT-id 
    """
    return  Ostap.Utils.rootID ( prefix )

# =============================================================================
## global ROOT identified for function objects 
def funcID  ( prefix = 'f_'  ) : return rootID  ( prefix )
## global ROOT identified for function objects 
def funID   ( prefix = 'f_'  ) : return funcID  ( prefix )
## global ROOT identified for function objects 
def fID     ( prefix = 'f_'  ) : return funcID  ( prefix )
## global ROOT identified for histogram objects 
def histoID ( prefix = 'h_'  ) : return rootID  ( prefix )
## global ROOT identified for histogram objects 
def histID  ( prefix = 'h_'  ) : return histoID ( prefix )
## global ROOT identified for histogram objects 
def hID     ( prefix = 'h_'  ) : return histoID ( prefix )
## global ROOT identified for dataset objects 
def dsID    ( prefix = 'ds_' ) : return rootID  ( prefix )
## global ROOT identified for graphs objects 
def grID    ( prefix = 'gr_' ) : return rootID  ( prefix )

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

_valid_pointer_ = Ostap.Utils.valid_pointer
# =============================================================================
## Is it a valid C++ pointer?
#  @code
#  ptr = ...
#  print 'Is the pointer valid? %s'  % valid_pointer ( prt ) 
#  @endcode 
#  @see Ostap::Utils::valid_pointer 
def valid_pointer ( obj ) :
    """ Is it a valid C++ pointer?
    - see Ostap::Utils::valid_pointer 
    >>> ptr = ...
    >>> print 'Is the C++ pointer valid? %s'  % valid_pointer ( ptr ) 
    """
    r = _valid_pointer_ ( obj )
    return True if r else False

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
## alowed argument for draw-function 
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
                          'copy' ,
                          'logx' ,
                          'logy' ] )
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
    def _TO_draw_ ( obj , option = '', *args , **kwargs ) :
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

        """

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

            if copy : result = obj.DrawCopy ( option , *args )
            else    : result = obj.Draw     ( option , *args )
                
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
            logger.error ( 'Exception here!!' , exc_info = True )
            
        # =====================================================================
        ## Lin/log scale ?
        if 'LogY' in kw or 'LogX' in kw :            
            pad = Ostap.Utils.get_pad() 
            if pad and 'LogX' in kw : pad.SetLogx ( kw.pop ( 'LogX' ) )
            if pad and 'LogY' in kw : pad.SetLogy ( kw.pop ( 'LogY' ) )
                
        # =====================================================================
        if kw : # =============================================================
            # =================================================================
            import ostap.logger.table as T 
            rows = [ ( 'Argument' , 'Value' ) ]
            for k , v in loop_items ( kw ) :
                row = k , str ( v )
                rows.append ( row )
            title = 'draw: %d unused arguments' % len ( kw ) 
            table = T.table ( rows , title = 'Unused arguments' , prefix = '# ' , alignment = 'll' )    
            logger.warning ( '%s\n%s' % ( title , table ) )
            
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
## `topdir': get the top directory for the given directory
#  @code
#  rdir = ...
#  tdir = tdit.top_dir 
#  tdir = tdit.topdir 
#  @endcode 
def top_dir ( rdir ) :    
    """`topdir': get the top directory for the given directory/object
    >>> rdir = ...
    >>> tdir = tdit.top_dir 
    >>> tdir = tdit.topdir 
    """
    
    if not rdir : return None

    with ROOTCWD()  :

        top = rdir 
        if   isinstance ( rdir , ROOT.TDirectory ) : top = rdir
        elif hasattr    ( rdir , 'GetDirectory'  ) : top = rdir.GetDirectory()
        else                                       : return None 
            
        while top : ## and hasattr ( top , 'GetMotherDir' ) : 
            moth = top.GetMotherDir()
            if not moth : return top  
            top  = moth
        else :
            return None 

# ROOT.TDirectory.top_dir = property ( top_dir , None , None , top_dir . __doc__ )
# ROOT.TDirectory.topdir  = property ( top_dir , None , None , top_dir . __doc__ )
ROOT.TNamed.top_dir = property ( top_dir , None , None , top_dir . __doc__ )
ROOT.TNamed.topdir  = property ( top_dir , None , None , top_dir . __doc__ )

# =============================================================================
## valid TDirectory?
#  - check valid C++ TDirectory pointer 
#  - for file directories check validity of the file
#  @code
#  odir = ...
#  if odir : ...
#  @endcode
def _rd_valid_ ( rdir ) :
    """ Valid TDirectory ?
    - check valid C++ TDirectory pointer 
    - for file directories check validity of the file 
    >>> odir = ...
    >>> if odir : ...
    """
    # =========================================================================
    ## check validity of C++ pointer 
    if not valid_pointer ( rdir ) : return False

    # ========================================================================
    ## for the file directories check also the validity of the file
    if isinstance ( rdir , ROOT.TDirectoryFile ) : # =========================
        # ====================================================================
        fdir = rdir.GetFile()
        if not valid_pointer ( fdir ) or ( not fdir.IsOpen () ) or fdir.IsZombie () :
            return False 
        
    return True 
        
ROOT.TDirectory.__bool__     = _rd_valid_
ROOT.TDirectory.__nonzero__  = _rd_valid_

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
## very simple context manager to suppress RooFit printout
#
#  @code
#
#  >>> with rooSilent( 4 , False ) :
#  ...        some_RooFit_code_here()
#
#  @endcode
#  @see RooMgsService
#  @see RooMgsService::globalKillBelow
#  @see RooMgsService::silentMode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-09
class RooSilent(object) :
    """ Very simple context manager to suppress RooFit printout
    
    >>> with rooSilent( 4 , False ) :
    ...        some_RooFit_code_here ()
    
    """
    ## constructor
    #  @param level  (INPUT) print level 
    #  @param silent (print level 
    # 
    def __init__ ( self , level = ROOT.RooFit.ERROR , silent = True ) :
        """ Constructor
        @param level  (INPUT) print level 
        @param silent (print level 
        
        >>> with rooSilent( ROOT.RooFit.ERROR , True  ) :
        ...        some_RooFit_code_here ()
        
        
        >>> with rooSilent( ROOT.RooFit.INFO , False  ) :
        ...        some_RooFit_code_here ()
        
        
        """
        #
        import ROOT
        #
        if level > ROOT.RooFit.FATAL : level = ROOT.RooFit.FATAL 
        if level < ROOT.RooFit.DEBUG : level = ROOT.RooFit.DEBUG 
        #
        self._level  = level 
        self._silent = True if silent else False  
        self._svc    = ROOT.RooMsgService.instance()
        
    ## context manager
    def __enter__ ( self ) :

        self._prev_level  = self._svc.globalKillBelow  () 
        self._prev_silent = self._svc.silentMode       () 
        
        self._svc.setGlobalKillBelow  ( self._level      )
        self._svc.setSilentMode       ( self._silent     )
        
        return self
    
    ## context manager 
    def __exit__ ( self , *_ ) : 
            
        self._svc.setSilentMode      ( self._prev_silent )
        self._svc.setGlobalKillBelow ( self._prev_level  )

# =============================================================================
## very simple context manager to suppress ROOT printout
#  @code
#  >>> with rootError () : some_ROOT_code_here()
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-30
def rootError   ( level = 1 ) :
    """ Very simple context manager to suppress ROOT printout
    >>> with rootError () : some_ROOT_code_here()
    """
    return ROOTIgnore ( ROOT.kError   + level )

# =============================================================================
## very simple context manager to suppress ROOT printout
#  @code
#  >>> with rootError () : some_ROOT_code_here()
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-30
def rootWarning ( level = 1 ) :
    """ Very simple context manager to suppress ROOT printout
    >>> with rootWarning () : some_ROOT_code_here()
    """
    return ROOTIgnore ( ROOT.kWarning + level )


# =============================================================================
## very simple context manager to suppress RooFit printout
#
#  @code
#
#  >>> with rooSilent( 4 , False ) :
#  ...        some_RooFit_code_here()
#
#  @endcode
#  @see RooMgsService
#  @see RooMgsService::globalKillBelow
#  @see RooMgsService::silentMode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-09
def rooSilent ( level = ROOT.RooFit.ERROR , silent = True ) :
    """ Very simple context manager to suppress RooFit printout
    >>> with rooSilent( 4 , False ) :
    ...        some_RooFit_code_here()    
    """
    return RooSilent ( level , silent ) 

# =============================================================================
## helper context manager
#  @code
#
#  >>> with roo_silent( True ) : 
#  ...        some_RooFit_code_here()
#
#  @endcode
#  @see rooSilent
#  @see NoContex
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-09
def roo_silent ( silence , *args ) :
    """ Helper context manager#
    >>> with roo_silent ( True ) : 
    ...        some_RooFit_code_here()
    """
    return rooSilent ( *args ) if silence else NoContext() 

# =============================================================================


# =============================================================================
## helper context manager to activate ROOT Error -> Python exception converter 
#  @see Ostap::Utils::useErrorHandler
#  @see Ostap::Utils::ErrorSentry
#  @code
#  with RootError2Exception() :
#  .... do something here 
#  @endcode 
class RootError2Exception (object) :
    """ Helper context manager to activate ROOT Error -> Python exception converter
    #
    with RootError2Exception() :
    ... do something here 
    """
    def __init__ ( self ) :
        self.e_handler  = Ostap.Utils.useErrorHandler 
        self.m_previous = False 

    ## context manager entry point  
    def __enter__ ( self ) :    
        self.m_previous = self.e_handler ( True ) 
        return self
    
    ## context manager exit point
    def __exit__ ( self , *_ ) :    
        if self.m_previous : self.e_handler ( False ) 
        self.m_previous = False 

# =============================================================================
## helper context manager to activate ROOT Error -> Python exception converter 
#  @see Ostap::Utils::useErrorHandler
#  @see Ostap::Utils::ErrorSentry
#  @code
#  with rootException () :
#  .... do something here 
#  @endcode
def rootException () :
    """ Helper context manager to activate ROOT Error -> Python exception converter
    #
    with rootException() :
    ... do something here 
    """
    return RootError2Exception()

# =============================================================================
## define the build directory for ROOT 
import ostap.core.build_dir

# =============================================================================
## come general configuration 
import ostap.core.config as _OCC

## if _OCC.general.getboolean ( 'ThreadSafety' , fallback = False )  :
##     logger.debug ("Thread safety is enabled via 'ROOT::ROOT::EnableThreadSAfety' call") 
##     ROOT.ROOT.EnableThreadSafety()

if _OCC.general.getboolean ( 'ImplicitMT' , fallback = False )  :
    if not ROOT.ROOT.IsImplicitMTEnabled() : 
        logger.debug ("Implicit MT is enabled")
        ROOT.ROOT.EnableImplicitMT  ()
else :
    if ROOT.ROOT.IsImplicitMTEnabled() : 
        logger.debug ("Implicit MT is disabled")
        ROOT.ROOT.DisableImplicitMT ()

# =============================================================================
## Use/Force BATCH processing ? 
if not ROOT.ROOT.GetROOT().IsBatch() :
    if any ( k.lower() in ( '-b' , '--batch' ) for k in sys.argv ) : 
        ROOT.ROOT.GetROOT().SetBatch ( True )
        logger.attention ( "BATCH processing is activated (command line) " )
    else : 
        from ostap.utils.root_utils import batch_env
        batch_env ( logger )

# =============================================================================
## Are we in CMAKE-test regime?
def in_test () :
    """ Are we in CMAKE-test regime?"""
    return os.environ.get ( 'OSTAP_CMAKE_TEST' , False )

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
