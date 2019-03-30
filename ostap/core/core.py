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
    'natural_entry'    ,  ## natual entry?   @see Gaudi::Math::natural_entry 
    'natural_number'   ,  ## natual numnber? @see Gaudi::Math::natural_number
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
    )
# =============================================================================
import ROOT, cppyy, math, sys
cpp = cppyy.gbl
std = cpp.std
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.core' )
else                       : logger = getLogger( __name__     )
# =============================================================================
logger.debug ( 'Core objects/classes/functions for Ostap')
# =============================================================================
from ostap.math.base      import ( Ostap    ,
                                   iszero   , isequal ,
                                   isint    , islong  ,
                                   inrange  , strings , 
                                   natural_number     ,
                                   natural_entry      )

from sys                  import version_info  as python_version 
from ostap.math.ve        import VE
from ostap.stats.counters import SE , WSE 
from builtins             import range 
#
binomEff        = Ostap.Math.binomEff
binomEff2       = Ostap.Math.binomEff2
zechEff         = Ostap.Math.zechEff
wilsonEff       = Ostap.Math.wilsonEff
agrestiCoullEff = Ostap.Math.agrestiCoullEff

# =============================================================================
## @class ROOTCWD
#  context manager to preserve current directory (rather confusing stuff in ROOT)
#  @code
#  print ROOT.gROOT.CurrentDirectory() 
#  with ROOTCWD() :
#     print ROOT.gROOT.CurrentDirectory() 
#     rfile = ROOT.TFile( 'test.root' , 'recreate' )
#     print ROOT.gROOT.CurrentDirectory() 
#  print ROOT.gROOT.CurrentDirectory() 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
class ROOTCWD(object) :
    """Context manager to preserve current directory
    (rather confusing stuff in ROOT) 
    >>> print ROOT.gROOT.CurrentDirectory() 
    >>> with ROOTCWD() :
    ...     print ROOT.gROOT.CurrentDirectory() 
    ...     rfile = ROOT.TFile( 'test.root' , 'recreate' )
    ...     print ROOT.gROOT.CurrentDirectory() 
    ... print ROOT.gROOT.CurrentDirectory() 
    """
    def __init__ ( self ) :
        self._dir = None
        
    ## context manager ENTER 
    def __enter__ ( self ) :
        "Save current working directory"
        self._dir = ROOT.gROOT.CurrentDirectory()
        return self 
        
    ## context manager EXIT 
    def __exit__  ( self , *_ ) :
        "Make the previous directory current again"
        if self._dir :

            odir = self._dir
            
            ## is is a directory in the file?
            fdir = odir.GetFile ()
            if not fdir : odir.cd()
            else        :
                ## check that fiel is still Open 
                if fdir.IsOpen() : odir.cd()
                
            
            #if  isinstance ( self._dir , ROOT.TFile ) :
            #    if self._dir.IsOpen() : self._dir.cd()
                
                
        self._dir = None 
            
# =============================================================================
## global identifier for ROOT objects 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def rootID ( prefix = 'o_' ) :
    """ Construct the unique ROOT-id 
    """
    _fun = lambda i : prefix + '%d'% i
    
    _root_ID = 1000
    ## 
    with ROOTCWD() : ## keep the current working directory:
        
        _id  = _fun ( _root_ID )
        grd  = ROOT.gROOT
        cwd  = grd.CurrentDirectory()
        while grd.FindObject    ( _id ) or \
              grd.FindObjectAny ( _id ) or \
              cwd.FindObject    ( _id ) or \
              cwd.FindObjectAny ( _id ) : 
                
            _root_ID += 10 
            _id       = _fun ( _root_ID ) 
            
    return _id                 ## RETURN
# =============================================================================
## global ROOT identified for function objects 
def funcID  () : return rootID  ( 'f_' )
## global ROOT identified for function objects 
def funID   () : return funcID  ( )
## global ROOT identified for function objects 
def fID     () : return funcID  ( )
## global ROOT identified for histogram objects 
def histoID () : return rootID  ( 'h_' )
## global ROOT identified for histogram objects 
def histID  () : return histoID ( )
## global ROOT identified for histogram objects 
def hID     () : return histoID ( )
## global ROOT identified for dataset objects 
def dsID    () : return rootID  ( 'ds_' )
## global ROOT identified for graphs objects 
def grID    () : return rootID  ( 'gr_' )

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
    return ROOT.gROOT.CurrentDirectory()

# =================================== ===============================================
## get current directory in ROOT
#  @code
#  print pwd() 
#  @endcode 
def pwd() :
    """ Get current directory in ROOT
    >>> print pwd() 
    """
    return ROOT.gROOT.CurrentDirectory().GetPath() 

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
    from ostap.logger.logger   import colored_string
    if   sc.isSuccess     () : return colored_string( 'SUCCESS'     , WHITE , GREEN  , True ) 
    elif sc.isRecoverable () : return colored_string( 'RECOVERABLE' , RED   , YELLOW , True ) 
    elif _FAILURE != sc.getCode  () :
        return colored_string('FAILURE[%d]' % sc.getCode() , YELLOW , RED   , True )
    return colored_string('FAILURE' , YELLOW , RED , True ) 
        
StatusCode = Ostap.StatusCode 
StatusCode .__repr__ = _sc_print_
StatusCode .__str__  = _sc_print_

SUCCESS = StatusCode(Ostap.StatusCode.SUCCESS)
FAILURE = StatusCode(Ostap.StatusCode.FAILURE)

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
    def _to_draw_ ( obj , *args , **kwargs ) :
        """ (silent) Draw of ROOT object
        >>> obj
        >>> obj.Draw()  ##
        >>> obj.draw()  ## ditto
        """
        from ostap.logger.utils import  rootWarning
        with rootWarning() :
            return obj.Draw( *args , **kwargs )

    ROOT.TObject.draw = _to_draw_

# =============================================================================
## Set/Set name/title 
# ==========================================================================
def _tn_name_get_ ( self )         : return self.GetName()
def _tn_name_set_ ( self , value ) : self.SetName( value )
_tn_name_doc_ = "``name'' of the object using GetName/SetName"
def _tn_title_get_ ( self )         : return self.GetTitle()
def _tn_title_set_ ( self , value ) : self.SetTitle( value )
_tn_title_doc_ = "``title'' of the object using GetTitle/SetTitle"
ROOT.TNamed.name  = property ( _tn_name_get_  ,  _tn_name_set_  , None , _tn_name_doc_  ) 
ROOT.TNamed.title = property ( _tn_title_get_ ,  _tn_title_set_ , None , _tn_title_doc_ ) 

# =============================================================================
## split string using separators: blanks,
#  @code
#  split_string ( ' a b cde,fg;jq', ',;:' )
#  @endcode
def split_string ( line , separators = ',;:' ) :
    """Split the string using separators
    >>> split_string ( ' a b cde,fg;jq', ',;:' )
    """
    items = line.split()
    for s in separators :
        result = []
        for item in items :
            if s in item : result += item.split(s)
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
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
