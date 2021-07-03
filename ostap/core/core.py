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
from   __future__        import print_function
# ============================================================================= 
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'cpp'              ,  ## global C++ namespace
    'std'              ,  ## C++ namespace std
    'Ostap'            ,  ## C++ namespace Ostap
    'ROOTCWD'          ,  ## context manager to keep/preserve ROOT current directory
    'rootID'           ,  ## global identifier for ROOT objects
    'funcID'           ,  ## global identifier for ROOT functions 
    'funID'            ,  ## global identifier for ROOT functions 
    'fID'              ,  ## global identifier for ROOT functions 
    'histoID'          ,  ## global identifier for ROOT histograms 
    'hID'              ,  ## global identifier for ROOT histograms 
    'grID'             ,  ## global identifier for ROOT graphs 
    'dsID'             ,  ## global identifier for ROOT/RooFit datasets
    ##
    'VE'               ,  ## shortcut for Gaudi::Math::ValuewithError
    'SE'               ,  ## shortcut for StatEntity
    'WSE'              ,  ## shortcut for Gaudi::Math::WStatEntity 
    ##
    'binomEff'         ,  ## binomial efficiency  
    'binomEff2'        ,  ## binomial efficiency
    'zechEff'          ,  ## binomial efficiency: Zech's recipe 
    'wilsonEff'        ,  ## binomial efficiency: Wilson 
    'agrestiCoullEff'  ,  ## binomial efficiency: Agresti-Coull
    ##
    'iszero'           ,  ## comparison with zero  for doubles  
    'isequal'          ,  ## comparison for doubles 
    'isint'            ,  ## Is float value actually int  ? 
    'islong'           ,  ## Is float value actually long ?
    'inrange'          ,  ## Is float walue in range ?  
    ##
    'natural_entry'    ,  ## natural entry?  @see Ostap::Math::natural_entry 
    'natural_number'   ,  ## natural nunber? @see Ostap::Math::natural_number
    ##
    'valid_pointer'    ,  ## Is it a valid C++ pointer?
    ##
    'strings'          , ## construct std::vector<std::string>
    'split_string'     , ## split the string  according to separators 
    ##
    'StatusCode'       , ## status code
    'SUCCESS'          , ## status code SUCCESS 
    'FAILURE'          , ## status code FAILURE
    ##
    'loop_items'       , ## loop over dictionary items 
    'items_loop'       , ## ditto
    ##
    'is_sorted'        , ## check that list is sorted 
    )
# =============================================================================
import math, sys, os 
import ROOT, cppyy
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.core.core' )
else                       : logger = getLogger( __name__     )
# =============================================================================
logger.debug ( 'Core objects/classes/functions for Ostap')
# =============================================================================
## ROOT.ROOT.EnableThreadSafety()

from   ostap.math.base      import ( Ostap    , std     , cpp ,  
                                     iszero   , isequal ,
                                     isint    , islong  ,
                                     inrange  , strings , 
                                     natural_number     ,
                                     natural_entry      )

from   sys                  import version_info  as python_version 
from   ostap.math.ve        import VE
from   ostap.stats.counters import SE , WSE 
from   builtins             import range
from   ostap.core.meta_info import root_info  
# =============================================================================

## @var global ROOT/gROOT object 
binomEff        = Ostap.Math.binomEff
binomEff2       = Ostap.Math.binomEff2
zechEff         = Ostap.Math.zechEff
wilsonEff       = Ostap.Math.wilsonEff
agrestiCoullEff = Ostap.Math.agrestiCoullEff
# =============================================================================
## helper function for case-insensitive dictionary with ignorance of underscores and blanks 
cidict_fun = lambda k : k.lower().replace('_','').replace(' ','')
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
class ROOTCWD(object) :
    """Context manager to preserve current directory
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
        del self._dir

    ## context manager ENTER 
    def __enter__ ( self ) :
        "Save current working directory"
        self._dir = None

        ## ROOT::TDirectory::TContext appears in ROOT 6/23/01
            
        groot     = ROOT.ROOT.GetROOT ()
        if groot :
            cwd = groot.CurrentDirectory()
            if root_info < ( 6, 23, 1 ) : pass
            else                        : cwd = cwd.load() 
            if cwd : self._dir = cwd
                           
            ## self._dir = ROOT.TDirectory.TContext()
            ## pass
        
        return self

    ## context manager EXIT 
    def __exit__  ( self , *_ ) :
        "Make the previous directory current again"

        ## ROOT::TDirectory::TContext appears in ROOT 6/23/01
        if self._dir :

            odir = self._dir
                
            self._dir = None
            
            fdir = odir.GetFile () if isinstance ( odir , ROOT.TDirectoryFile ) else None
            
            if fdir and not fdir.IsOpen () :
                
                groot = ROOT.ROOT.GetROOT ()
                groot.cd ()
                
            else :
                odir.cd()

            del self._dir
            
        self._dir = None
            
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
    """Print the Status Code
    >>> st = ...
    >>> print st
    """
    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = list ( range ( 8 ) )
    ##
    from ostap.logger.colorized import colored_string
    if   sc.isSuccess     () : return colored_string ( 'SUCCESS'     , WHITE , GREEN  , True ) 
    elif sc.isRecoverable () : return colored_string ( 'RECOVERABLE' , RED   , YELLOW , True ) 
    elif _FAILURE != sc.getCode  () :
        return colored_string ( 'FAILURE[%d]' % sc.getCode() , YELLOW , RED   , True )
    return colored_string ( 'FAILURE' , YELLOW , RED , True ) 
        
StatusCode = Ostap.StatusCode 
StatusCode .__repr__ = _sc_print_
StatusCode .__str__  = _sc_print_

SUCCESS = StatusCode ( Ostap.StatusCode.SUCCESS )
FAILURE = StatusCode ( Ostap.StatusCode.FAILURE )

_valid_pointer_ = Ostap.Utils.valid_pointer
# =============================================================================
## Is it a valid C++ pointer?
#  @code
#  ptr = ...
#  print 'Is the pointer valid? %s'  % valid_pointer ( prt ) 
#  @endcode 
#  @see Ostap::Utils::valid_pointer 
def valid_pointer ( obj ) :
    """Is it a valid C++ pointer?
    - see Ostap::Utils::valid_pointer 
    >>> ptr = ...
    >>> print 'Is the C++ pointer valid? %s'  % valid_pointer ( ptr ) 
    """
    r = _valid_pointer_ ( obj )
    return True if r else False

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

        """
        
        from ostap.utils.cidict import cidict
        kw = cidict ( transform = cidict_fun , **kwargs )
        
        ## Line
        
        if 'LineColor'  in kw and hasattr ( obj , 'SetLineColor' ) :
            obj.SetLineColor   ( kw.pop('LineColor' ) )
        if 'LineStyle'  in kw and hasattr ( obj , 'SetLineStyle' ) :
            obj.SetLineStyle   ( kw.pop('LineStyle' ) )
        if 'LineWidth'  in kw and hasattr ( obj , 'SetLineWidth' ) :
            obj.SetLineWidth   ( kw.pop('LineWidth' ) )

        ## Marker
            
        if 'MarkerColor' in kw and hasattr ( obj , 'SetMarkerColor' ) :
            obj.SetMarkerColor ( kw.pop('MarkerColor' ) )
        if 'MarkerStyle' in kw and hasattr ( obj , 'SetMarkerStyle' ) :
            obj.SetMarkerStyle ( kw.pop('MarkerStyle' ) )
        if 'MarkerSize'  in kw and hasattr ( obj , 'SetMarkerSize'  ) :
            obj.SetMarkerSize  ( kw.pop('MarkerSize'  ) )

        ## Area
            
        if 'FillColor'   in kw and hasattr ( obj , 'SetFillColor' ) :
            obj.SetFillColor   ( kw.pop('FillColor' ) )
        if 'FillStyle'   in kw and hasattr ( obj , 'SetFillStyle' ) :
            obj.SetFillStyle   ( kw.pop('FillStyle' ) )

        ## Min/max values  
            
        if   'Minimum'     in kw and hasattr ( obj , 'SetMinimum' ) :
            obj.SetMinimum     ( kw.pop ( 'Minimum' ) )
        elif 'Min'         in kw and hasattr ( obj , 'SetMinimum' ) :
            obj.SetMinimum     ( kw.pop ( 'Min'     ) )
        if   'Maximum'     in kw and hasattr ( obj , 'SetMaximum' ) :
            obj.SetMaximum     ( kw.pop ( 'Maximum' ) )
        elif 'Max'         in kw and hasattr ( obj , 'SetMaximum' ) :
            obj.SetMaximum     ( kw.pop ( 'Max'     ) )


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
                                        
                
        if kw : logger.warning('draw: unknown attributes: %s' % kw.keys() )
            
        from ostap.logger.utils import  rootWarning,   rooSilent 
        with rootWarning() , rooSilent ( 2 )  :
            result = obj.Draw( option , *args )
            
        groot = ROOT.ROOT.GetROOT ()
        pad   = groot.GetSelectedPad()
        if   pad and not ROOT.gPad :
            c = pad.GetCanvas()
            if c : c.Update()
        elif ROOT.gPad :
            c = ROOT.gPad
            if c : c.Update()
            
        return result 

    ROOT.TObject.draw       = _TO_draw_
    ROOT.TObject.draw_ostap = _TO_draw_

# =============================================================================
## Set/Set name/title 
# =============================================================================
def _tn_name_get_ ( self )         : return self.GetName()
def _tn_name_set_ ( self , value ) : self.SetName( value )
_tn_name_doc_ = "``name'' of the object using GetName/SetName"
def _tn_title_get_ ( self )         : return self.GetTitle()
def _tn_title_set_ ( self , value ) : self.SetTitle( value )
_tn_title_doc_ = "``title'' of the object using GetTitle/SetTitle"
ROOT.TNamed.name  = property ( _tn_name_get_  ,  _tn_name_set_  , None , _tn_name_doc_  ) 
ROOT.TNamed.title = property ( _tn_title_get_ ,  _tn_title_set_ , None , _tn_title_doc_ ) 

# =============================================================================
## Get the full path of the named object
#  @code
#  obj  = ...
#  path = obj.path
#  @endcode 
def _tn_path_ ( obj ) :
    """Get the full path of the named object
    >>> obj  = ...
    >>> path = obj.path
    """
    if not valid_pointer ( obj ) : return "<INVALID-OBJECT>"
    d = obj.GetDirectory() if hasattr ( obj , 'GetDirectory' ) else None 
    if not d : return obj.GetName()
    dp = d.GetPath()
    h , s , p = dp.rpartition(':/')
    if p : return  '/'.join ( [ p , obj.GetName () ] )
    return obj.GetName() 

ROOT.TNamed.path = property ( _tn_path_ , None , None , None  ) 

# =============================================================================
## split string using separators:
#  @code
#  split_string ( ' a b cde,fg;jq', ',;:' )
#  @endcode
def split_string ( line , separators = ',;:' ) :
    """Split the string using separators
    >>> split_string ( ' a b cde,fg;jq', ',;:' )
    """
    
    if ' ' in separators : items = line.split()
    else                 : items = [ line ]
    
    for s in separators :
        result = []
        for item in items :
            if s in item : result += item.split ( s )
            else         : result.append ( item ) 
        items = result

    ## remove empty items 
    while '' in items: items.remove ( '' )

    return items 

# =============================================================================
if python_version.major > 2 : items_loop = lambda d : d.items     () 
else                        : items_loop = lambda d : d.iteritems ()
# =============================================================================
def loop_items ( d ) :
    """Return  iterator over the dictionary   items
    >>> d = { 'a' : ...   , 'b' : ... , }
    >>> for e in   loop_items ( d ) : print (e) 
    """
    return items_loop ( d ) 

# =============================================================================
## valid TDirectory?
#  - check valid C++ TDirectory pointer 
#  - for file directories check validity of the file
#  @code
#  odir = ...
#  if odir : ...
#  @endcode
def _rd_valid_ ( rdir ) :
    """Valid TDirectory ?
    - check valid C++ TDirectory pointer 
    - for file directories check validity of the file 
    >>> odir = ...
    >>> if odir : ...
    """

    ## check validity of C++ pointer 
    if not valid_pointer ( rdir ) : return False

    ## for the file directories check the validity of the file 
    if isinstance ( rdir , ROOT.TDirectoryFile ) :
        fdir = rdir.GetFile()
        if not valid_pointer ( fdir ) or not fdir.IsOpen () or fdir.IsZombie () :
            return False 

    return True 
        
ROOT.TDirectory.__bool__     = _rd_valid_
ROOT.TDirectory.__nonzero__  = _rd_valid_


# =============================================================================
## check that list is sorted 
def is_sorted ( lst ) :
    """Check that list is sorted
    """
    l = len ( lst )  
    return all ( lst [ i ] <= lst [ i + 1 ] for i in range ( l - 1 ) ) if lst else True 


# =============================================================================
## define the build directory for ROOT 
import ostap.core.build_dir 

# =============================================================================
_decorated_classes_ = (
    ROOT.TObject    ,
    ROOT.TNamed     ,
    ROOT.TDirectory , 
    )

_new_methods_       = (
    #
    ROOT.TObject.draw  ,
    ROOT.TNamed .name  ,
    ROOT.TNamed .title ,
    ROOT.TNamed .path  ,
    )

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
