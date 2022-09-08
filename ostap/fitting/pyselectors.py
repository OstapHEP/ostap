#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/fitting/pyselectors.py
# 
# Helper module to fix a problems in communication of
# <c>TTree/TChain::Process</c> and <c>TPySelector</c>.
# In PyROOT some of original C++ methods are disable.
# The module provides the "recovery" for missing methods
#
# @see TTree::Process
# @see TChain::Process
# @see TPySelector
#
# @code
#
# from ostap.fitting.pyselectors import Selector
#
# class MySelector ( Selector ) :
#
#      def __init__ ( self ) :
#
#           return Selector ( self )
#
#      def Process  ( self , entry ) :
#           # == getting the next entry from the tree
#           if self.GetEntry ( entry ) <= 0 : return 0             ## RETURN 
#        
#           # == for more convenience
#           tree = self.tree
#
#           # apply trivial "acceptance" cuts 
#           if not 2 <= tree.y   <=  4.5   : return 0               ## RETURN
#           if not 1 <= tree.pt  <=  10.0  : return 0               ## RETURN
#
#           return 1
# 
# selector = MySelector()
#
# chain = ...
# chain.process ( selector )  ## NB: note lowercase "process" here !!!
#
# @endcode
#
# The module has been developed and used with great success in
# "Kali, framework for fine calibration of LHCb Electromagnetic Calorimeter"
#
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date   2010-04-30
# 
# =============================================================================
"""Helper module to fix a problems in communication of
TTree/TChain.Process and TPySelector.

In PyROOT some of original C++ methods are disable.
The module provides the 'recovery' for missing methods

# from ostap.fitting.pyselectors import Selector
#
# class MySelector ( Selector ) :
#
#    def __init__ ( self ) :
#       return Selector ( self )
#    def Process  ( self , entry ) :
#        # == getting the next entry from the tree
#        if self.GetEntry ( entry ) <= 0 : return 0             ## RETURN 
#    
#        # == for more convenience
#        tree=self.fChain
#    
#        # apply trivial 'acceptance' cuts 
#        if not 2 <= tree.y   <=  4.5   : return 0               ## RETURN
#        if not 1 <= tree.pt  <=  10.0  : return 0               ## RETURN
#    
#        return 1
#
#
# selector = MySelector()
# 
# chain = ...
#
# chain.process ( selector )  ## NB: note lowercase 'process' here !!!
#

The module has been developed and used with great success in
'Kali, framework for fine calibration of LHCb Electormagnetic Calorimeter'

More complicated (and more useful) cases are covered by
SelectorWithCuts and SelectorWithVars classes 

"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2010-04-30"
__version__ = "$Revision$" 
# =============================================================================
__all__ = (
    'Selector'         ,        ## The "fixed" TPySelector
    'SelectorWithCuts' ,        ## The "fixed" TPySelector with TTree-formula 
    'SelectorWithVars' ,        ## Generic selctor to fill RooDataSet from TTree/TChain
    'Variable'         ,        ## helper class to define variable
    'DataSet_NEW_FILL' ,        ## Is new efficient filing machinery activated? 
    'SelectorWithVarsCached'    ## Generic selector with cache   
)
# =============================================================================
from   ostap.core.meta_info     import root_info 
from   ostap.core.core          import ( cpp  , Ostap , items_loop ,
                                         dsID , valid_pointer , binomEff ) 
from   ostap.core.ostap_types   import num_types, string_types, integer_types
from   ostap.core.meta_info     import old_PyROOT
import ostap.fitting.roofit 
from   ostap.utils.progress_bar import ProgressBar
from   ostap.math.reduce        import root_factory 
from   ostap.trees.trees        import Chain 
import ROOT, cppyy, math, sys
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger      import getLogger
from   ostap.logger.colorized   import attention, allright
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.pyselectors' )
else                       : logger = getLogger ( __name__          )
# =============================================================================
_new_methods_ = []
# =============================================================================
## @class Selector
#  Useful intermediate class for implementation of (py)selectors 
#  @see Ostap::Selector
#  @see TPySelector
#  @see TSelector
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-09
class Selector ( Ostap.Selector ) :
    """Useful intermediate class for implementation of (py)selectors     
    """
    ## constructor 
    def __init__ ( self , tree = None , silence = False  ) :

        if isinstance ( tree , Chain ) : tree = tree.chain 
        if tree is None : tree = ROOT.nullptr
        
        if old_PyROOT : super (Selector, self).__init__ ( self , tree ) 
        else          : super (Selector, self).__init__ (        tree ) 
        """Standart constructor
        """
    @property 
    def tree ( self ) :
        """'tree' : get the actual TTree pointer"""
        return self.get_tree() 
    
    # =========================================================================
    ## the actual major method
    #  it needs to be redefined 
    def process_entry ( self ) :
        """The actual major method
        - it needs to be redefined
        """
        raise NotImplementedError ("Selector: process_entry is not implemented!")
        return True
    
    ## reduce the object 
    def __reduce__ ( self ) :
        """Reduce the object"""
        tree = self.tree 
        if tree : return root_factory , ( type ( self ) , Chain ( tree ) ) 
        else    : return root_factory , ( type ( self ) , ) 
        
# =============================================================================
if old_PyROOT :

    # =========================================================================
    ## Start master processing
    #  @see TPySelector::Begin
    #  @see Ostap::Selector::Begin
    def Selector_Begin      ( self , tree = ROOT.nullptr  ) :
        """Start master processing
        - see ROOT.TPySelector.Begin
        - see Ostap.Selector.Begin
        """
        pass
    # =========================================================================
    ## Start slave processing 
    #  @see TPySelector::SlaveBegin
    #  @see Ostap::Selector::SlaveBegin
    def Selector_SlaveBegin ( self , tree ) :
        """Start master processing
        - see ROOT.TPySelector.SlaveBegin
        - see Ostap.Selector.SlaveBegin
        """
        assert valid_pointer ( tree ) , 'SlaveBegin: invalid TTree*'
        pass
    # =========================================================================
    ## Initialize
    #  - Invoked for <code>2<=Version</code>
    #  @see TPySelector::Init
    #  @see Ostap::Selector::Init 
    def Selector_Init  ( self , tree ) :
        """Initialize 
        - see ROOT.TPySelector.Init
        - see Ostap.Selector.Init
        """
        pass
    # =========================================================================
    ## Terminate slave processing 
    #  @see TPySelector::SlaveTerminate 
    #  @see Ostap::Selector::SlaveTerminiate 
    def Selector_SlaveTerminate ( self ) :
        """Terminate slave processing
        - see ROOT.TPySelector.SlaveTerminiate 
        - see Ostap.Selector.SlaveTerminate
        """
        pass
    # =========================================================================
    ## Terminate master processing 
    #  @see TPySelector::Terminate 
    #  @see Ostap::Selector::Terminiate 
    def Selector_Terminate ( self ) :
        """Terminate master processing
        - see ROOT.TPySelector.Terminiate 
        - see Ostap.Selector.Terminate
        """
        pass
    
    # =========================================================================
    ## Version
    #  @see TPySelector::Version
    #  @see Ostap::Selector::Version
    def Selector_Version ( self ) :
        """Version 
        - see ROOT.TPySelector.Version
        - see Ostap.Selector.Version
        """
        return 3 
    # ========================================================================= 
    ## The major method
    # -  note that  call for check for GetEntry is already called!
    # -  No need to redefine this method 
    # @see TPySelector::Process 
    # @see TPySelector::GetEntry 
    # @see Ostap::Selector::Process
    def Selector_Process        ( self , entry    ) :
        """The major method
        -  note that  call for check for GetEntry is already called!
        -  No need to redefine this method  
        -  see Ostap::Selector::Process
        """
        return self.process_entry ()

    Selector.Begin          = Selector_Begin
    Selector.SlaveBegin     = Selector_SlaveBegin
    Selector.Init           = Selector_Init
    Selector.Version        = Selector_Version
    Selector.Terminate      = Selector_Terminate
    Selector.SlaveTerminate = Selector_SlaveTerminate
    Selector.Process        = Selector_Process 

# =============================================================================
## @class SelectorWithCuts
#  Efficient selector that runs only for "good"-events  
#  @see Ostap::SelectorWithCuts
#  @see Ostap::Selector
#  @see TPySelector
#  @see TSelector
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-09
class SelectorWithCuts (Ostap.SelectorWithCuts) :
    """Efficient selector that runs only for 'good'-events  
    """
    ## constructor 
    def __init__ ( self             , 
                   selection        ,
                   silence = False  ,
                   tree    = None   ,
                   logger  = logger ) :
        """ Standart constructor
        """
        if   isinstance ( tree , Chain ) : tree = tree.chain 
        elif tree is None : tree = ROOT.nullptr
        
        self.__silence = silence
        self.__logger  = logger
        
        ## initialize the base
        self.__selection = str ( selection ).strip()
        
        
        if   old_PyROOT : super ( SelectorWithCuts , self ).__init__ ( self , self.selection , tree )
        else            : super ( SelectorWithCuts , self ).__init__ (        self.selection , tree )
        
        if self.cuts () and not self.silence :
            self.logger.info ( 'SelectorWithCuts: %s' % self.cuts() )

    @property
    def silence ( self ) :
        """'silence'  : silent processing?"""
        return self.__silence 
    
    @property
    def selection ( self ) :
        """'selection' -  selection to be used to preprocess TTree/TChain"""
        return self.__selection

    @property
    def tree ( self ) :
        """'tree' : get the tree/chain  itself"""
        return self.get_tree()  

    @property
    def logger ( self ) :
        """'logger' : get the local logger object"""
        return self.__logger 
    
    # =========================================================================
    ## the actual major method to process "good" entry
    #  it needs to be redefined
    #  @returns <code>True</code> for successfull processing 
    def process_entry ( self ) :
        """The actual major method
        - it needs to be redefined!
        - returns True for successfull processing 
        """
        raise NotImplementedError ("SelectorWithCuts: process_entry is not implemented!")
        return True

    ## reduce of the object 
    def __reduce__ ( self ) :
        """Reduce the object"""
        tree = self.tree
        if tree : return root_factory , ( type ( self )  ,
                                          self.selection ,
                                          self.silence   ,                                
                                          Chain ( tree ) )
        else    : return root_factory , ( type ( self )  ,
                                          self.selection ,
                                          self.silence   ) 
            
# =============================================================================
if old_PyROOT :

    # =========================================================================
    ## The major method for <cdoe>SelectorWithCuts</code>: do NOT redefine it!
    #  @see TPySelector::Process
    #  @see Ostap::Selector::Process
    #  @see Ostap::SelectorWithCuts::Process
    def SelectorWithCuts_Process        ( self , entry ) :
        """The major method for `SelectorWithCuts`: do NOT redefine it!
         - see TPySelector::Process
         - see Ostap::Selector::Process
         - see Ostap::SelectorWithCuts::Process
         """
        return self.process_entry () 
    
    SelectorWithCuts.Begin          = Selector_Begin
    SelectorWithCuts.SlaveBegin     = Selector_SlaveBegin
    SelectorWithCuts.Init           = Selector_Init  
    SelectorWithCuts.Version        = Selector_Version
    SelectorWithCuts.Terminate      = Selector_Terminate
    SelectorWithCuts.SlaveTerminate = Selector_SlaveTerminate
    SelectorWithCuts.Process        = SelectorWithCuts_Process 

    
# =============================================================================
_maxv =  0.99 * sys.float_info.max
_minv = -0.99 * sys.float_info.max
# =============================================================================
## @class Variable
#  Helper   structure to manage/keep/create the variable for   SelectorWithVars
#  @see SelectorWithVars
class Variable(object) :
    """Helper   structure to manage/keep/create the variable for   SelectorWithVars
    
    - Get a variable 'my_name1' from the tree/chain:
    
    >>> v = Variable ( 'my_name1' , 'my_description1' , -100 , 100 ) ]
    
    - Get a variable 'my_name' from the tree/chain using the explicit accessor function, making some on-fly transforomation:
    
    >>> v = Variable ( 'my_name2' , 'my_description2' , -100 , 100 , lambda s : s.my_name2/1000 ) ]
    
    - Use less    trivial expressions:
    
    >>> v = Variable ( 'my_name3' , 'my_description3' , -1  , 2 , lambda s : s.var1+s.var2 ) ]
    
    - Any callable that gets TChain/Tree and evaluates to double( useful case - e.g. it could be TMVAReader)
    
    >>> def myvar ( chain ) : ...
    >>> v = Variable ( 'my_name4' , 'my_description4' , accessor = myvar )  ]

    - Use already booked variables:

    >>> v5 = ROOT.RooRealVar( .... ) 
    >>> v  = Variable ( v5 , accessor = lambda s : s.var5 ) ]

    - Use already booked variables:
    
    >>> v6 = ROOT.RooRealVal( 'my_name6' )
    >>> variables += [  Variable ( v6 ) ] ## get variable 'my_name6'

    """
    def __init__ ( self                ,
                   var                 ,
                   description = ''    ,
                   vmin        = _minv , 
                   vmax        = _maxv , 
                   accessor    = None  ) :
        
        assert var and isinstance ( var ,  ( str , ROOT.RooRealVar ) ) , \
               "Variable: invalid type for 'var'  %s/%s"      % ( var         , type (  var        ) ) 
        assert accessor is None or callable ( accessor ) or isinstance ( accessor , str )  , \
               "Variable: illegal type for 'accessor' %s/%s"  % ( accessor    , type ( accessor    ) )

        if isinstance ( var , str ) :

            var = var.strip() 
            assert isinstance ( description , str ) , \
                   "Variable: illegal type for 'description'"     % ( description , type ( description ) ) 
            assert isinstance ( vmin , num_types ) , \
                   "Variable: illegal type for 'vmin' %s/%s"      % ( vmin        , type ( vmin        ) )
            assert isinstance ( vmax , num_types ) , \
                   "Variable: illegal type for 'vmax' %s/%s"      % ( vmax        , type ( vmax        ) )
            assert vmin < vmax, \
                   "Variable: invalid 'minmax' range (%g,%g)"     % ( vmin , vmax ) 

            if    description                   : pass
            elif  isinstance ( accessor , str ) : description = accessor
            else                                : description = "'%s'-variable" % var
            
            description = description.replace('\n',' ')

            ## create the variable
            var = ROOT.RooRealVar ( var , description , vmin , vmax )
            
        ## redefine min/max  
        self.__vmin , self.__vmax = var.minmax()                    
        
        self.__formula = None
        if accessor is None :
            varname        = var.GetName()
            accessor       = varname
            self.__formula = varname 
            
        if isinstance  ( accessor , str ) :
            accessor = accessor.strip() 
            if accessor  : 
                from ostap.trees.funcs import FuncFormula as FuncVar
                self.__formula  = accessor 
                accessor = FuncVar ( accessor ) 

        assert callable ( accessor ), \
               'Invalid accessor function! %s/%s' % ( accessor , type ( accessor ) )
        
        self.__var         = var
        self.__minmax      = var.minmax()
        self.__accessor    = accessor
        self.__checked     = False
        
    @property
    def var         ( self ) :
        """'var' : the variable itself/ROOT.RooRealVar"""
        return self.__var
    @property
    def name        ( self ) :
        """'name' : the name of the variable"""
        return self.__var.GetName() 
    @property
    def description ( self ) :
        """'description' : the variable description/title"""
        return self.__var.GetTitle()
    @property
    def minmax      (  self ) :
        """'minmax' : the range for the variable """
        return self.__minmax    
    @property
    def vmin ( self ) :
        """'vmin' : minimal value (if set) """
        if not self.__vmin is None : return self.__vmin
        mnmx = self.minmax
        if mnmx : return mnmx  [0]
        return None    
    @property
    def vmax ( self ) :
        """'vmax' : maximal value (if set) """
        if not self.__vmax is None : return self.__vmax
        mnmx = self.minmax
        if mnmx : return mnmx [1]
        return None    
    @property
    def accessor    ( self ) :
        """'accessor' - the actual callable to get the value from TTree/TChain"""
        return self.__accessor
    @property
    def formula ( self ) :
        """'formula' - formula for this variable (when applicable)?"""
        return self.__formula
    
    @property
    def trivial ( self ) :
        """'trivial' - is this variable 'trivial'?"""
        return self.__formula ## and self.__formula == self.name

    @property
    def really_trivial ( self ) :
        """'really_trivial' - is it really trivial enough for RooFit?"""
        triv = self.trivial 
        return triv and  ( not '[' in triv ) and ( not ']' in triv )
    
    @property
    def checked ( self ) :
        """'checked' : checked?
        """
        return self.__checked

    ## reduce the object
    def __reduce__ ( self ) :
        """Reduce the object"""
        return root_factory , ( type ( self )    ,
                                self.var.name    ,
                                self.description , 
                                self.vmin        ,
                                self.vmax        ,
                                self.formula if self.formula else self.accessor )

# =============================================================================
## @class Variables
#  Helper structure to manage/keep/check list of variables
class Variables(object) :
    """Helper structure to manage/keep/check list of variables
    """
    def __init__ ( self , variables ) :

        self.__variables = []
        self.__triv_vars = True

        names = set()
        vset  = ROOT.RooArgSet()
        
        for v in variables :

            if   isinstance   ( v , str              ) : vv = Variable (   v ) 
            elif isinstance   ( v , ROOT.RooAbsReal  ) : vv = Variable (   v )
            elif isinstance   ( v , ( tuple , list ) ) : vv = Variable (  *v )
            elif isinstance   ( v , dict             ) : vv = Variable ( **v )
            elif isinstance   ( v , Variable         ) : vv = v  
            
            assert isinstance ( vv , Variable ), 'Invalid variable %s/%s' % ( vv , type ( vv ) )

            name = vv.name
            if name in names :
                logger.error ( "Variable %s is already in the list, skip it" % name ) 
                continue

            if name in vset  : 
                logger.error ( 'Variable  %s is already defined! skip it' % name )
                continue

            if vv.var in vset : 
                logger.error ( 'Variable  %s is already defined! skip it' % name )
                continue
            
            names.add ( name )
            
            self.__variables.append ( vv     )
            vset.add    ( vv.var )
            #
            if   vv.trivial and vv.name == vv.formula : pass
            elif vv.really_trivial                    : pass
            elif vv.formula                           : pass
            else                                      :
                self.__triv_vars = False

        self.__variables = tuple ( self.__variables ) 

    @property
    def variables ( self ) :
        """'variables' :  the actual list/tuple of variables"""
        return self.__variables

    @property
    def trivial_vars( self ) :
        """'trivial_vars' : are all variables 'trivial' (suitable for fast-processing)?"""
        return self.__triv_vars

    ## number of variables 
    def __len__ ( self ) : 
        """Number of variables"""
        return len ( self.__variables )
    ## get item
    def __getitem__ ( self , item ) :
        """Get certain item"""
        return self.__variables [ item ]
    ## iterator 
    def __iter__    ( self ) :
        """Make an iteration over variables"""
        return iter ( self.__variables )

    # =========================================================================
    ## check if variable (or index) is already in the list 
    def __contains__ ( self , item ) :
        """Check if variable (or index) is already in the list
        """
        
        if isinstance ( item , int ) :
            return 0 < item < len ( self.__variables )

        if isinstance ( item , Variable ) :
            return ( item.name in self ) and ( item.var in self ) 
        
        if isinstance ( item , string_types ) :
            name = str ( item )
            for v in self.__variables :
                if name == v.name : return True
            return False
        
        if isinstance ( item , ROOT.RooAbsArg ) :
            return item in self.__varset
        
        return False

    ## reduce the object 
    def __reduce__ ( self ) :
        ## redcue the object 
        return root_factory , ( type ( self ) , self.variables  ) 
  
# =============================================================================
## Is this expression corresponds to a valid RooFit formula?
def valid_formula ( expression , varset ) :
    """Is this expression corresponds to a valid RooFit formula?
    """

    if isinstance ( expression , ROOT.TCut ) : expression =  str ( expression )
    
    expression = expression.strip()
    if not expression : return True
    
    if isinstance  ( varset  , ROOT.RooArgSet ) :
        vlst = ROOT.RooArgList()
        for v in varset : vlst.add ( v  )
        result = valid_formula ( expression , vlst )
        del vlst
        return result

    assert isinstance ( varset  , ROOT.RooArgList ), 'Invalid type %s' % type (varset)

    
    ## from ostap.logger.utils import rooSilent, rootErro
    ## with rooSilent ( ROOT.RooFit.FATAL + 1 , True ) :
    ## with rootError( ROOT.kError + 1 ) :
    
    from ostap.logger.utils import mute 
    with mute ( True , True ) :
        return Ostap.validFormula ( expression , varset )
            
# ==============================================================================
## @class SelStat
#  Helper class to keep the statististics for SelectorWithVars 
class SelStat(object) :
    """Helper class to keep the statististics for SelectorWithVars
    """
    def __init__ ( self ,  total = 0 , processed = 0 , skipped = 0 ) :

        assert 0<= total and 0 <= processed and 0 <= skipped , \
               "Invalid counters: %s/%s/%s" % ( total, processed , skipped)
        
        self.__total     = total
        self.__processed = processed 
        self.__skipped   = skipped  
        
    @property
    def total ( self ) :
        """'total'   : total number of events"""
        return self.__total
    @total.setter
    def total ( self , value ) :
        assert isinstance ( value , integer_types ) and 0 <= value ,\
               "Invalid value for 'total'"
        self.__total = value
            
    @property
    def processed ( self ) :
        """'processed'   : number of processed events"""
        return self.__processed
    @processed.setter
    def processed ( self , value ) :
        assert isinstance ( value , integer_types ) and 0 <= value ,\
               "Invalid value for 'processed'"
        self.__processed = value
        
    @property
    def skipped ( self ) :
        """'skipped'   : number of skipped events (e.g. due to variabel ranges)"""
        return self.__skipped
    @skipped.setter
    def skipped ( self , value ) :
        assert isinstance ( value , integer_types ) and 0 <= value ,\
               "Invalid value for 'skipped'"
        self.__skipped = value

    def __str__ ( self ) :
        return 'SelStat(total=%s,processed=%s,skipped=%s)' % ( self.__total     ,
                                                               self.__processed ,
                                                               self.__skipped   )
    __repr__ = __str__

    ## reduce the object 
    def __reduce__ ( self ) :
        """Reduce the object"""
        return root_factory , ( type ( self )   ,
                                self.total      ,
                                self.processed  ,
                                self.skipped    )
    
# ==============================================================================
## Define generic selector to fill RooDataSet from TChain
#
#  @code
#  variables = [ ... ]
#  ## add a variable 'my_name1' from the tree/chain 
#  variables += [ 
#   #  name       descriptor           min-value , max-value  
#   Variable ( 'my_name1' , 'my_description1' , low       , high     )
#   ]
# 
#  ## get a variable 'my_name' from the tree/chain with accessor function, e.g. rescale it on-fligh
#  variables += [ 
#   #  name       descriptor           min-value , max-value , access function   
#   Variable ( 'my_name2' , 'my_description2' , low       , high      , lambda s : s.my_name2/1000 )
#   ]
#   
#  ## get  less trivial expression
#  variables += [ 
#   #  name       descriptor           min-value , max-value , access function   
#   Variable ( 'my_name3' , 'my_description3' , low       , high      , lambda s : s.var1+s.var2 ) 
#  ]
#
#  ## any function that gets TChain/Tree entry and evaluates to double.
#  #  e.g. it could be TMVAReader
#  def myvar ( chain ) : ....
#  variables += [ 
#   #  name       descriptor           min-value , max-value , access function   
#   Variable ( 'my_name4' , 'my_description4' , low       , high      , myvar ) 
#  ]
#
#  ## add already booked variables: 
#  v5 = ROOT.RooRealVal( 'my_name5' )
#  variables += [  Variable ( v5 , accessor = lambda s : s.var5 ) ]
#
#  ## add already booked variables: 
#  v6 = ROOT.RooRealVal( 'my_name6' )
#  variables += [  Variable ( v6                     ) ] ## get variable "my_name6"
#
#  ## finally create selector
#  # 
#  selector = SelectorWithVars (
#             variables                             ,
#             selection = " chi2vx<30 && pt>2*GeV " ,  ## filtering
#            )
#  chain = ...
#  chain.process ( selector )
#  dataset = selector.dataset
#  @endcode
#  @date   2014-03-02
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  - thanks to  Alexander BARANOV 
class SelectorWithVars(SelectorWithCuts) :
    """Create and fill the basic dataset for RooFit
    
    - Define the list of 'variables' for selector:
    
    >>> variables = [ ... ]
    
    Add a variable 'my_name1' from the tree/chain:
    
    >>> variables += [ # name       descriptor         min-value , max-value  
    ...    Variable ( 'my_name1' , 'my_description1' , low       , high     ) ]
    
    Get a variable 'my_name' from the tree/chain using the accessor function, e.g. rescale it on-fligh:
    
    >>> variables += [ #  name       descriptor        min-value , max-value , access function   
    ...    Variable ( 'my_name2' , 'my_description2' , low       , high      , lambda s : s.my_name2/1000 ) ]
    
    Use less trivial expression:
    
    >>> variables += [ #  name       descriptor        min-value , max-value , access function   
    ...    Variable ( 'my_name3' , 'my_description3' , low       , high      , lambda s : s.var1+s.var2 ) ]
    
    Any callable that gets TChain/Tree and evaluates to double.
    ( useful case - e.g. it could be TMVAReader)
    
    >>> def myvar ( chain ) : ...
    >>> variables += [ #  name       descriptor        min-value , max-value , access function   
    ...    Variable ( 'my_name4' , 'my_description4' , low       , high      , myvar )  ]

    Use already booked variables:
    
    >>> v5 = ROOT.RooRealVal( 'my_name5' )
    >>> variables += [  Variable ( v5 , accessor = lambda s : s.var5 ) ]

    Add already booked variables:
    
    >>> v6 = ROOT.RooRealVal( 'my_name6' )
    >>> variables += [  Variable ( v6 ) ] ## get variable 'my_name6'

    - Finally create selector

    >>> selector = SelectorWithVars (
    ...       variables                             ,
    ...       selection = ' chi2vx<30 && pt>2*GeV ' ) ## filtering

    - Use selector to fill RooDataSet 
    >>> tree  = ...
    >>> chain.process ( selector )

    - Get dataset  from the selector 
    >>> dataset = selector.data   
    """
    ## constructor 
    def __init__ ( self                            ,
                   variables                       ,  ## list of variables  
                   selection                       ,  ## Tree-selection 
                   cuts          = None            ,  ## python function  
                   roo_cuts      = ''              ,  ## selection based on dataset variables (RooFormula)  
                   name          = ''              ,
                   fullname      = ''              ,
                   silence       = False           ,
                   tree          = ROOT.nullptr    ,
                   logger        = logger          ) :

        
        if not name : name = dsID ()
            
        if not fullname : fullname = name 

        self.__name     = name
        self.__fullname = fullname 
        #
        ## create the logger 
        #
        assert 0 < len(variables) , "Empty list of variables"
        #
        ## instantiate the base class
        # 
        SelectorWithCuts.__init__ ( self ,
                                    selection = selection ,
                                    silence   = silence   ,
                                    tree      = tree      ,
                                    logger    = logger    ) ## initialize the base
        
        self.__cuts      = cuts
        self.__variables = [] 
        self.__varset    = ROOT.RooArgSet()
        
        ## vvars  = set ()
        ## vnames = set ()
        ## for v in variables :

        ##     vv = v 
        ##     if   isinstance ( v , str              ) : vv = Variable (   v ) 
        ##     elif isinstance ( v , ROOT.RooAbsReal  ) : vv = Variable (   v )
        ##     elif isinstance ( v , ( tuple , list ) ) : vv = Variable (  *v )
        ##     elif isinstance ( v , dict             ) : vv = Variable ( **v )
        ##     elif isinstance ( v , Variable         ) : vv = v  

        ##     assert isinstance  ( vv , Variable ), 'Invalid variable %s/%s' % ( vv , type ( vv ) )

        ##     if vv.name in self.__varset :
        ##         logger.error ( 'Variable  %s is already defined! skip it' % vv.name )
        ##         continue 

        ##     self.__variables.append ( vv     )
        ##     self.__varset   .add    ( vv.var )
        ##     #
        ##     if   vv.trivial and vv.name == vv.formula : pass
        ##     elif vv.really_trivial                    : pass
        ##     elif vv.formula                           : pass
        ##     else                                      :
        ##         self.__triv_vars = False
        ##     #
        ##     vvars.add ( vv )
            
        ## self.__variables = tuple( self.__variables ) 

        self.__variables = Variables ( variables ) 
        for v in self.__variables : self.__varset.add ( v.var )
            
        assert 1 <= len ( self.__variables ), \
               'Invalid setting of variablesl!'
        
        self.__triv_sel  = valid_formula ( selection , self.varset ) 
        triv_cuts        = not cuts
        
        self.__trivial = self.trivial_vars and self.__triv_sel and triv_cuts
        if not self.silence :
            tv = allright ( 'True' ) if self.trivial_vars else attention ( 'False' )
            ts = allright ( 'True' ) if self.__triv_sel   else attention ( 'False' )
            tc = allright ( 'True' ) if triv_cuts         else attention ( 'False' )
            self.logger.info ( "Suitable for fast processing: variables:%s, selection:%s, py-cuts:%s" % ( tv , ts , tc ) )
            
        if not self.silence: 
            nl = 0
            dl = 0 
            for v in self.__variables :
                nl = max ( nl , len( v.name        ) ) 
                dl = max ( dl , len( v.description ) )                 
            dl = max ( dl , len ( 'Description' ) + 2 ) 
            nl = max ( nl , len ( 'Variable'    ) + 2 ) 

            fmt_name = '%%%ds' % nl
            fmt_desc = '%%%ds' % dl
            fmt_min  = '%+11.3g'
            fmt_max  = '%-+11.3g'
            fmt_triv = '%4s' 
            
            header = ( ( '{:^%d}' % nl ).format ( 'Variable'    ) ,
                       ( '{:^%d}' % dl ).format ( 'Description' ) ,
                       ( '{:>11}'      ).format ( 'min '        ) ,
                       ( '{:<11}'      ).format ( ' max'        ) ,
                       ( '{:^8}'       ).format ( 'Trivial?'    ) )
            
            table_data =  [ header ]
            for v in  self.__variables :
                triv = allright ( '    {:<4}'.format('yes') ) if v.really_trivial else attention ( '    {:<4}'.format ( 'no' ) )            
                triv = allright ( 'Yes' ) if v.really_trivial else attention ( 'No' )            
                table_data.append ( ( fmt_name % v.name        ,
                                      fmt_desc % v.description ,
                                      fmt_min  % v.minmax [0]  ,
                                      fmt_max  % v.minmax [1]  , triv ) )
                
            if fullname != self.name : 
                title = '%s("%s","%s") %d variables' % ( 'RooDataSet', self.name , fullname , len ( self.varset ) )
            else :
                title = '%s("%s") %d variables'      % ( 'RooDataSet', self.name ,            len ( self.varset ) )
                        
            import ostap.logger.table as T
            t  = T.table (  table_data , title , '# ' )
            self.logger.info ( "Booked dataset: %s\n%s" % ( title , t ) ) 

            
        ## Book dataset
        self.__data = ROOT.RooDataSet (
            ##
            self.name ,
            fullname  , 
            ##
            self.varset
            )

        ## selection using RooFit machinery 
        self.__roo_cuts    = ''
        self.__roo_formula = None

        if roo_cuts :
            varlist = self.__data.varlist() 
            roosel  = Ostap.FormulaVar    ( roo_cuts, varlist , False )          
            assert roosel .ok() , "SelectorWithVars: invalid 'roo_cuts': %s" % roo_cuts
            del    roosel
            used    = Ostap.usedVariables ( roo_cuts , varlist  )            
            roosel  = Ostap.FormulaVar    ( roo_cuts , used     , True )
            self.__roo_formula = roosel 
            self.__roo_cuts    = roo_cuts 
    
        ## it is still very puzzling for me: should this line be here at all??
        ROOT.SetOwnership ( self.__data  , False )
        
        self.__progress = None 
        from collections import defaultdict
        self.__skip     = defaultdict(int)
        self.__notifier = None
        self.__stat     = SelStat() 
        self.__last     = -1
        
    @property 
    def name ( self ) :
        """'name'  : the name of selector/dataset"""
        return self.__name
    
    @property 
    def fullname ( self ) :
        """'fullname' : the fullname of selector/dataset"""
        return self.__fullname 
    
    @property 
    def data ( self ) :
        """'data'  : the dataset"""
        return self.__data
    @data.setter
    def data ( self , dataset ) :
        assert isinstance ( dataset , ROOT.RooAbsData ), \
               "Incorrect type of data %s/%s " % ( dataset ,   type ( dataset ) )
        self.logger.debug ("Selector(%s), add dataset %s" % (  self.__name , dataset ) )
        self.__data = dataset 

    @property 
    def variables ( self ) :
        """'variables' : the list/tuple of variables (cleared in Terminate)"""
        return self.__variables
    
    @property
    def varset ( self ) :
        """'varset' : the structure of RooDataSet"""
        return self.__varset
    
    @property
    def morecuts ( self ) :
        """'morecuts' : additional cust to be applied in selection"""
        return self.__cuts

    @property
    def roo_cuts( self ) :
        """'roo-cuts' : addtional selection/cuts based on RooFit machinery (RooFormula)"""
        return self.__roo_cuts

    @property
    def trivial_vars( self ) :
        """'trivial_vars' : are all variables 'trivial' (suitable for fast-processing)?"""
        return self.__variables.trivial_vars

    @property
    def really_trivial ( self ) :
        """'really_trivial' : is a set of variables really trivial (for RooFit)?"""
        for v in self.variables :
            if not v.really_trivial : return False
        return  not '[' in self.selection and not ']' in self.selection 

    @property
    def trivial_sel( self ) :
        """'trivial_sel' : is the selection 'trivial' (suitable for fast-processing)?"""
        return self.__triv_sel
    
    @property
    def trivial ( self ) :
        """'trivial' : Are variables/selection/cuts 'trivial' (suitable for fast-processing)?"""
        return self.__trivial

    @property
    def skip ( self ) :
        """'skip' : dictionary of skept entries"""
        return self.__skip
    
    @property
    def skipped ( self ) :
        """'skipped': total number of skipped entries"""
        return self.stat.skipped
    
    @property
    def processed  ( self ) :
        """'processed' : number of processeed events (after cuts)"""
        return self.stat.processed 
    
    @property
    def total  ( self ) :
        """`total' : total number of processeed events (before cuts)"""
        return self.stat.total

    @property
    def stat ( self ) :
        """'stat' : Total/processed/skipped events"""
        return self.__stat    
    @stat.setter
    def stat ( self , value  ) :
        assert isinstance ( value , SelStat ) , 'Invalid "value":%s' % str ( value )
        self.__stat = value 

    ## get the dataset 
    def dataset   ( self  ) :
        """ Get the data-set """ 
        return self.__data

    ## # =========================================================================
    ## ## the only one actually important method 
    ## def Process ( self, entry ):
    ##     """ Fill data set 
    ##     """
    ##     #
    ##     ## == getting the next entry from the tree
    ##     #
    ##     if self.GetEntry ( entry ) <=  0 : return 0             ## RETURN 
    ##     #
    
    # =========================================================================
    ## the only one actually important method 
    def process_entry ( self ):
        """ Fill data set 
        """

        if not self.__progress and not self.silence :
            tree            = self.tree
            self.stat.total = tree.GetEntries()
            self.logger.info ( "Selector(%s): processing TChain('%s') #entries: %d" % ( self.name , tree.GetName() , self.total ) )
            ## decoration:
            from ostap.utils.progress_bar import ProgressBar
            self.__progress = ProgressBar ( max_value = self.total   ,
                                            silent    = self.silence )
        if not self.silence :
            
            evnt  = self.event ()
            total = self.total
            
            step  = 10000 if total <= 1 else divmod ( total , 350 ) [0]
            
            if evnt    < 20 or 0 < total < evnt + 20 or self.__last + step <= evnt :
                self.__progress.update_amount ( evnt )
                self.__last = evnt 
            
        self.stat.processed += 1
        
        #
        ## == for more convenience
        #
        bamboo = self.tree 
        ##        
        return self.fill ( bamboo )
        
    # =========================================================================
    ## fill it! 
    def fill ( self , bamboo ) :
        """The  actual processing for the given 'bamboo'
        Note that   this method is independent on TTree/TChain and can be used directy
        One just needs to  ensure that:
        - 'accessor functions' for the variables and 'cuts' agree with the type of 'bamboo'
        """
        
        ## apply cuts (if needed) 
        if self.__cuts and not self. __cuts ( bamboo )  : return 0 

        ## loop over all variables
        for v in self.__variables :

            var       = v.var                ## The variable
            vmin,vmax = v.minmax             ## min/max range 
            vfun      = v.accessor           ## accessor function

            ## use the accessor function 
            value     = vfun ( bamboo )
            if not vmin <= value <= vmax :   ## MUST BE IN RANGE!
                self.__skip[v.name] += 1     ## SKIP EVENT
                self.stat.skipped   += 1     ## SKIP EVENT 
                return 0                     ## RETURN 

            var.setVal ( value ) 

        ## no roo-cuts are specified or roo-cuts are satisfied 
        if ( not self.__roo_formula ) or self.__roo_formula.getVal() : 
            self.__data .add ( self.varset )
            
        return 1 

    # =========================================================================
    ## 'callable' interface 
    def __call__ ( self ,  entry ) :
        """'callable' interface to Selector
        """
        return self.fill ( entry ) 

    # =========================================================================
    ## clone this selector
    #  @code
    #  sel = ...
    #  new_sel = sel.clone ( name = 'QUQU' , fullname = 'FullName  ) 
    #  @endcode
    def clone ( self , **kwargs ) :
        """Clone the selector
        >>> sel = ...
        >>> new_sel = sel.clone ( name = 'QUQU' , fullname = 'FullName  ) 
        """
        
        kw = {}
        kw [ 'variables' ] = self.variables
        kw [ 'selection' ] = self.selection 
        kw [ 'roo_cuts'  ] = self.roo_cuts
        kw [ 'cuts'      ] = self.morecuts
        kw [ 'silence'   ] = self.silence

        kw.update ( kwargs )
        return SelectorWithVars ( **kw ) 

    # =========================================================================
    ## termination 
    def Terminate ( self  ) :
        #
        ## Aborted? 
        if   0 != self.GetAbort() :
            self.logger.fatal('Selector(%s): process has been aborted!' % self.__name )

            self.__data = None 
            del self.__varset
            del self.__variables
            self.__varset      = ()
            self.__variables   = ()
            self.__roo_formula = None
            
            return  ## RETURN
        
        if self.roo_cuts :
            del self.__roo_formula
            self.__roo_formula = None
            
        ## get total number of input events from base class 
        self.stat.total = self.event()
        
        if not self.silence :
            skipped = 'Skipped:%d' % self.skipped
            skipped = '/' + attention ( skipped ) if self.skipped else ''
            cuts    = allright ( '"%s"' % self.cuts () ) if self.trivial_sel else attention ( '"%s"'  % self.cuts() ) 
            self.logger.info (
                'Selector(%s): Events Total:%d/Processed:%d%s CUTS: %s' % (
                self.__name    ,
                self.total     ,
                self.processed ,
                skipped        , 
                cuts           ) )            
            
        if self.__data and not self.silence :
            vars = []
            for v in self.__variables :
                s    = self.__data.statVar( v.name )
                mnmx = s.minmax ()
                mean = s.mean   ()
                rms  = s.rms    ()
                r    = ( v.name        ,                       ## 0 
                         v.description ,                       ## 1 
                         ('%+.5g' % mean.value() ).strip() ,   ## 2
                         ('%.5g'  % rms          ).strip() ,   ## 3 
                         ('%+.5g' % mnmx[0]      ).strip() ,   ## 4
                         ('%+.5g' % mnmx[1]      ).strip() )   ## 5
                s = self.__skip [ v.name] 
                if s : skip = '%-d' % s
                else : skip = '' 
                r +=  skip,                                    ## 6 
                vars.append ( r )

            vars.sort()
            
            name_l  = len ( 'Variable'    ) + 2 
            desc_l  = len ( 'Description' ) + 2 
            mean_l  = len ( 'mean' ) + 2 
            rms_l   = len ( 'rms'  ) + 2
            min_l   = len ( 'min'  ) + 2 
            max_l   = len ( 'max'  ) + 2 
            skip_l  = len ( 'Skip' ) 
            for v in vars :
                name_l = max ( name_l , len ( v [ 0 ] ) )
                desc_l = max ( desc_l , len ( v [ 1 ] ) )
                mean_l = max ( mean_l , len ( v [ 2 ] ) )
                rms_l  = max ( rms_l  , len ( v [ 3 ] ) )
                min_l  = max ( min_l  , len ( v [ 4 ] ) )
                max_l  = max ( max_l  , len ( v [ 5 ] ) )
                skip_l = max ( skip_l , len ( v [ 6 ] ) )

            sep      = '# -%s+%s+%s+%s+%s-' % ( ( name_l       + 2 ) * '-' ,
                                                ( desc_l       + 2 ) * '-' ,
                                                ( mean_l+rms_l + 5 ) * '-' ,
                                                ( min_l +max_l + 5 ) * '-' ,
                                                ( skip_l       + 2 ) * '-' )
            fmt = '#   %%%ds | %%-%ds | %%%ds / %%-%ds | %%%ds / %%-%ds | %%-%ds   '  % (
                name_l ,
                desc_l ,
                mean_l ,
                rms_l  ,
                min_l  ,
                max_l  ,
                skip_l
                )
            
            report  = 'Dataset(%s) created:' % self.__name
            report += ' ' + allright ( '%s entries, %s variables' %  ( len ( self.__data ) , len ( self.variables ) ) )
            if self.trivial_vars : report += ' Vars:' + allright  ('trivial'       ) + ';'
            else                 : report += ' Vars:' + attention ('non-trivial'   ) + ';'
            if self.trivial_sel  : report += ' Cuts:' + allright  ('trivial'       ) + ';'
            else                 : report += ' Cuts:' + attention ('non-trivial'   ) + ';'
            if not self.__cuts   : report += ' '      + allright  ( 'no py-cuts'   )  
            else                 : report += ' '      + attention ( 'with py-cuts' )

            
            fmt_name = '%%-%ds' % name_l 
            fmt_desc = '%%-%ds' % desc_l
            fmt_mean = '%%%ds'  % mean_l
            fmt_rms  = '%%-%ds' % rms_l
            fmt_min  = '%%%ds'  % min_l
            fmt_max  = '%%-%ds' % max_l
            fmt_skip = '%%-%ds' % skip_l
            
            header = ( ( '{:^%d}' % name_l ).format ( 'Variable'    ) ,
                       ( '{:^%d}' % desc_l ).format ( 'Description' ) ,
                       ( '{:^%d}' % mean_l ).format ( 'mean'        ) ,
                       ( '{:^%d}' % rms_l  ).format ( 'rms'         ) ,
                       ( '{:^%d}' % min_l  ).format ( 'min'         ) ,
                       ( '{:^%d}' % max_l  ).format ( 'max'         ) ,
                       ( '{:^%d}' % skip_l ).format ( 'skip'        ) )
            
            table_data = [  header ]
            for v in vars :
                table_data.append ( ( fmt_name % v [ 0 ] ,
                                      fmt_desc % v [ 1 ] ,
                                      fmt_mean % v [ 2 ] ,
                                      fmt_rms  % v [ 3 ] ,
                                      fmt_min  % v [ 4 ] ,
                                      fmt_max  % v [ 5 ] ,
                                      attention ( fmt_skip % v [ 6 ] ) if v[6]  else v[6] ) ) 
            
            import ostap.logger.table as T
            title = report 
            t  = T.table ( table_data , title , '# ')
            self.logger.info ( title + '\n' + t )
            
            
        if not self.__data or not len ( self.__data ) :
            skip = 0
            for k,v in items_loop ( self.__skip ) : skip += v 
            self.logger.warning("Selector(%s): empty dataset! Total:%s/Processed:%s/Skipped:%d"
                                  % ( self.__name  , self.total , self.processed , skip ) ) 
            
        ## attention: delete these

        del self.__varset
        del self.__variables
        
        self.__varset     =  ()
        self.__variables  =  ()
        
    # =========================================================================
    ## Initialize the selector
    #  @see Ostap::SelectorWithCuts::Init 
    def Init    ( self, tree ) :
        """Initialize the selector
        - see Ostap::SelectorWithCuts::Init 
        """

        ## reset the formula 
        self.reset_formula ( tree )
        if valid_pointer ( tree ) :
            assert self.ok(), 'Init: formula is invalid!'
            
        if self.__progress and not self.silence :
            self.__progress.update_amount ( self.event () )


    # ===========================================================
    ## Start master processing
    #  @see Ostap::SelectorWithCuts::Begin
    def Begin          ( self , tree = ROOT.nullptr ) :
        """Start master processing
        - see Ostap::SelectorWithCuts::Begin
        """

        ## reset the formula 
        self.reset_formula ( tree ) 

        # reset Roo-formula (if specified) 
        if self.roo_cuts :
            roo_cuts = self.roo_cuts 
            del self.__roo_formula 
            varlist = self.__data.varlist() 
            roosel  = Ostap.FormulaVar    ( roo_cuts, varlist , False )          
            assert roosel .ok() , "SelectorWithVars: invalid 'roo_cuts': %s" % roo_cuts
            del    roosel
            used    = Ostap.usedVariables ( roo_cuts , varlist  )            
            roosel  = Ostap.FormulaVar    ( roo_cuts , used     , True )
            self.__roo_formula = roosel 
            
        ## take care on the progress bar 
        if self.__progress and not self.silence :
            self.__progress.update_amount ( self.event () )
            
        if valid_pointer ( tree ) :
            assert self.ok(), 'Begin: formula is invalid!'

    # =========================================================================
    ## Start slave processing
    #  @see Ostap::SelectorWithCuts::SlaveBegin
    def SlaveBegin     ( self , tree ) :
        """Start slave processing
        - see Ostap::SelectorWithCuts::SlaveBegin
        """
        #
        assert valid_pointer ( tree ),  'SlaveBegin:TTree is invalid!'
        #
        ## reset the formula 
        self.reset_formula  ( tree  ) 
        #
        assert self.ok () , 'SlaveBegin::Formula is invalid!'
        #
        #
        # reset Roo-formula (if specified) 
        if self.roo_cuts :
            roo_cuts = self.roo_cuts 
            del self.__roo_formula 
            varlist = self.__data.varlist() 
            roosel  = Ostap.FormulaVar    ( roo_cuts, varlist , False )          
            assert roosel .ok() , "SelectorWithVars: invalid 'roo_cuts': %s" % roo_cuts
            del    roosel
            used    = Ostap.usedVariables ( roo_cuts , varlist  )            
            roosel  = Ostap.FormulaVar    ( roo_cuts , used     , True )
            self.__roo_formula = roosel 
            
        ## take care on the progress bar 
        if not self.__progress and not self.silence :
            tree            = self.tree
            self.stat.total = tree.GetEntries()
            self.logger.info ( "Selector(%s): processing TChain('%s') #entries: %d" % ( self.name , tree.GetName() , self.total ) )
            ## decoration:
            from ostap.utils.progress_bar import ProgressBar
            self.__progress = ProgressBar ( max_value = self.total   ,
                                            silent    = self.silence )
            
        if self.__progress and not self.silence :
            self.__progress.update_amount ( self.event () )
        #
        self.stat.total = tree.GetEntries()
        #
        if self.__notifier :
            self.__notifier.exit()
            del self.__notifier
        
        self.__notifier = Ostap.Utils.Notifier( tree )
        for v in self.__variables :
            if isinstance ( v.accessor , ROOT.TObject ) :
                self.__notifier.add  ( v.accessor ) 
        
    # =========================================================================
    ## Notify  (e.g. another TTree in the chain
    #  @see Ostap::SelectorWithCuts::Notify 
    def Notify         ( self ) :
        """Notify  (e.g. another TTree in the chain
        - see Ostap::SelectorWithCuts::Notify
        """
        #
        result = True 
        if self.formula() : result = self.formula().Notify()
        # 
        if self.__progress and not self.silence :
            self.__progress.update_amount ( self.event () )
            
        return result 

    # ========================================================================
    ## Terminate slave processing
    #  @see Ostap::SelectorWithCuts::SlaveTerminate
    def SlaveTerminate ( self               ) :
        """Terminate slave processing
        - see Ostap::SelectorWithCuts::SlaveTerminate
        """
        
        if self.roo_cuts :
            del self.__roo_formula
            self.__roo_formula = None
            
        if self.__progress and not self.silence :
            self.__progress.update_amount ( self.event () )
            self.__progress.end () 
            
        if self.__notifier :
            self.__notifier.exit()
            self.__notifier = None  

    # =========================================================================
    ## reduce the object 
    def __reduce__ ( self ) :
        """Reduce the object"""
        tree = self.tree
        if tree : return root_factory , ( type ( self )  ,
                                          self.variables ,
                                          self.selection ,                                
                                          self.cuts      ,
                                          self.roo_cuts  ,
                                          self.name      ,
                                          self.fullname  ,
                                          self.silence   ,
                                          Chain ( tree ) )
        else    : return root_factory , ( type ( self )  ,
                                          self.variables ,
                                          self.selection ,                                
                                          self.cuts      ,
                                          self.roo_cuts  ,
                                          self.name      ,
                                          self.fullname  ,
                                          self.silence   )
    
# =============================================================================
import os
from   ostap.core.workdir import workdir
from   ostap.io.zipshelve import ZipShelf 
# =============================================================================

# ==============================================================================
## @class  SelectorWithVarsCached
#  Generic selector which loads already loaded datasets from cache
#  @date   2014-07-02
#  @author Sasha Baranov a.baranov@cern.ch
class SelectorWithVarsCached(SelectorWithVars) :
    """Create and fill the basic dataset for RooFit. Or just load it from cache.
    """
    ## constructor 
    def __init__ ( self                           ,
                   variables                      ,  ## list of variables  
                   selection                      ,  ## Tree-selection 
                   files                          ,  ## List of files
                   cuts         = None            ,  ## Tree-based cuts 
                   roo_cuts     = ''              ,  ## RooFit-based cuts 
                   name         = ''              ,
                   fullname     = ''              ,
                   tree         = ROOT.nullptr    ,
                   logger       = logger          ) : 

        SelectorWithVars.__init__( self ,
                                   variables = variables ,
                                   selection = selection ,
                                   cuts      = cuts      ,
                                   roo_cuts  = roo_cuts  ,
                                   name      = name      ,
                                   fullname  = fullname  ,
                                   silence   = silence   ,
                                   tree      = tree      ,
                                   logger    = logger    )

        # Try load from cache
        self._loaded_from_cache = False
        self._cachepath         = self._compute_cache_name()

        if self._loadcache():
            self.Process = lambda entry: 1
            self._loaded_from_cache = True

    def _l_internals(self, lmbd):
        " Returns str of lambda expression internals"
        if not hasattr(lmbd, "__code__"):
            return ""

        code = lmbd.__code__
        return str((code.co_code, code.co_consts))

    def _get_files_str(self):
        "Returns string with used files and their modification times"
        ret = ""
        for filename in sorted(self.__filelist):
            ret += filename
            ret += str(int(os.stat(filename).st_mtime))
        return ret

    def _compute_cache_name(self):
        " Computes dataset cache path(SHA512) "
        h =  self._get_files_str()
        h += self._l_internals(self.morecuts)

        for var , vdesc , vmin , vmax , vfun in self._variables:
            h += var.GetName() + vdesc + str(vmin) + str(vmax)
            h += self._l_internals(vfun)

        h += str(self.selection)

        import hashlib 
        filename = hashlib.sha512(h).hexdigest() + ".shelve"
        return os.path.join(workdir, "cache", filename)

    def _loadcache(self):
        " Loads cache from ZipShelf. Returns `True` on success, `False` othervise"
        if not os.path.exists(self._cachepath):
            return False

        with ZipShelf(self._cachepath, 'r' ) as db : 
            self.data = db['data']

        return True

    def _savecache(self):
        
        " Saves cache to file. "
        with ZipShelf(self._cachepath) as db : 
            db['data'] = self.data

        return True
    
    #
    def Terminate ( self  ) :
        if not self._loaded_from_cache:
            SelectorWithVars.Terminate(self)

            if len(self.data):
                self._savecache()

        else:
            self._logger.info('Loaded from cache!')

        return 1

# =============================================================================
## Is new (very efficient) fill machinery activated? 
DataSet_NEW_FILL = ( 6 , 26 ) <= root_info 

# =============================================================================
## Create RooDataset from the tree
#  @code 
#  tree = ...
#  ds   = tree.make_dataset_old ( [ 'px , 'py' , 'pz' ] ) 
#  @endcode
def make_dataset_old ( tree              ,
                       variables         , ## varibales 
                       selection = ''    , ## TTree selection 
                       roo_cuts  = ''    , ## Roo-Fit selection  
                       name      = ''    , 
                       title     = ''    ,
                       silent    = False ) :
    """Create the dataset from the tree
    >>> tree = ...
    >>> ds = tree.make_dataset ( [ 'px , 'py' , 'pz' ] ) 
    """

    if DataSet_NEW_FILL :
        if not silent : logger.info ( "Switch to more efficient function 'make_dataset'" )
        return make_dataset ( tree                  ,
                              variables = variables ,
                              seletion  = selection ,
                              roo_cuts  = roo_cuts  ,
                              name      = name      ,
                              title     = title     ,
                              silent    = silent    )
    
    import ostap.trees.cuts
    import ostap.fitting.roofit

    variables = Variables ( variables )
    if not variables.trivial_vars or ( '[' in selection ) or ( '[' in roo_cuts ) :
        if not silent : logger.info ( "Variables/selection are not trivial, switch to 'fill_dataset'" )
        return fill_dataset ( tree                  ,
                              variables = variables ,
                              selection = selection ,
                              roo_cuts  = roo_cuts  ,
                              name      = name      ,
                              title     = title     ,
                              shortcut  = True      ,
                              use_frame = True      , ## IMPORTANT!
                              silent    = silent    ) 
    
    
    varset   = ROOT.RooArgSet()
    vars     = set()

    formulas  = []
    
    selection = str ( selection ) if isinstance ( selection , ROOT.TCut ) else selection  
    selection = selection.strip() if isinstance ( selection , str       ) else selection 

    branches = set ( tree.branches () )
    leaves   = set ( tree.leaves   () ) 
    
    limits    = []

    for vv in variables :
        
        if vv.name != vv.formula :            
            if vv.name in branches :
                logger.error ( "Variable %s exists as branch! skip it" % vv.name )
                continue
            if vv.name in leaves   :
                logger.error ( "Variable %s exists as leaf! skip it"   % vv.name )
                continue

        if vv.trivial and vv.name == vv.formula : 
            
            assert hasattr  ( tree , vv.name ) , "Tree/Chain has no branch '%s'" % vv.name
            
            varset.add  ( vv.var )
            vars.add    ( vv )
            
        elif vv.formula :
            
            formulas.append ( vv )
            continue

        else :
            
            logger.error("Do not know how to treat the variable %s, skip it" % vv.name )
            continue 
        
        mn , mx = vv.minmax
        if _minv < mn : limits.append ( "(%.16g <= %s)" % ( mn      , vv.name ) ) 
        if _maxv > mx : limits.append ( "(%s <= %.16g)" % ( vv.name , mx      ) )


    ## 
    limits = ROOT.TCut (' && '.join ( limits ) ) if limits else ROOT.TCut()
    cuts   = ROOT.TCut ( limits )
    
    ## extended varset
    stor    = set() 
    varsete = ROOT.RooArgSet()
    for v in varset : varsete.add ( v )

    expressions = [ f.formula for f in formulas ]

    if selection :
        cuts = cuts & ROOT.TCut ( selection ) 
        expressions.append      ( selection ) 

    if expressions :

        tt = None 
        if isinstance ( tree , ROOT.TChain ) :
            nf = len ( tree.files() )
            for i in range ( nf ) :
                tt = tree[i]
                if tt : break 
        if not tt : tt = tree

        lvars = tt.the_variables ( *expressions )

        assert not lvars is None , 'Unable to get the basic variables for %s' % expressions
        if not silent :
            logger.info ("make_dataset: temporary variables to be added %s" % str ( lvars ) ) 
        for lname in lvars :
            if not lname in varsete :
                v = Variable ( lname )
                varsete.add  ( v.var )
                stor.add ( v )
                
    if not name : name = dsID () 
    if not title and tree.GetName() != tree.GetTitle  :
        title = tree.GetTitle ()

    total     = len ( tree )
    processed = tree.statVar ( '1' , selection ).nEntries()
    skipped   = tree.statVar ( '1' , cuts      ).nEntries() 
 
    stat = SelStat ( total , processed , processed - skipped )

    f1 = len ( varsete ) * 1.0 / len ( tree.branches() )
    f2 = skipped         * 1.0 / max ( 1 , total ) 
    if f1   < 0.20 and f2 < 0.10 and not silent :
        logger.warning ( "Only tiny fraction of variables (% 4.1f%%) and data (% 4.1F%%) is requested: prefiltering can speedup process" % ( f1 * 100 , f2 * 100 ) ) 
    elif f1 < 0.20 and not silent :
        logger.warning ( "Only tiny fraction of variables (% 4.1f%%) is requested: prefiltering can speedup process" % ( f1 * 100  ) )
    elif f2 < 0.10 and not silent : 
        logger.warning ( "Only tiny fraction of data (% 4.1F%%) is requested: prefiltering can speedup process" % ( f2 * 100 ) ) 
        
    from ostap.core.core    import rooSilent, rootError, NoContext
    from ostap.utils.timing import timing
    TIMING = timing if not silent else NoContext

    with TIMING ( 'Fill RooDataSet' , logger = logger ) : 
        with rooSilent ( ROOT.RooFit.ERROR  , True ) :
            with rootError ( ROOT.kWarning ) :
                ds      = ROOT.RooDataSet ( name  , title , tree , varsete , str ( cuts ) )
                varsete = ds.get()
                
    if not silent :
        logger.debug( "make_dataset: Initial dataset\n%s" % ds.table ( prefix = "# " ) ) 
    
    ## add complex expressions 
    if formulas or roo_cuts :
        # a
        vset = ds.get()
        vlst = ROOT.RooArgList()
        for v in vset : vlst.add ( v )

        fcols = ROOT.RooArgList() 

        ffs   = []
        fcuts = []
        for f in formulas :

            fv = Ostap.FormulaVar ( f.name , f.description , f.formula , vlst , False )
            assert fv.ok() , 'Invalid formula: %s' % f.formula 
            ffs.append ( fv )
            fcols.add  ( fv )
            mn , mx = f.minmax            
            if _minv < mn : fcuts.append ( "(%.16g <= %s)" % ( mn      , fv.name ) )
            if _maxv > mx : fcuts.append ( "(%s <= %.16g)" % ( fv.name , mx      ) )

        with rooSilent ( ROOT.RooFit.ERROR + 1 , True ) :
            with rootError( ROOT.kError ) :
                ## ## it causes some strange behaviour
                ## ds.addColumns ( fcols )
                ## ## it is ok: 
                for  f in fcols : ds.addColumn ( f ) 
                del fcols
                del ffs 

        ## insert Roo-Fit cuts here!
        if roo_cuts : fcuts.append ( roo_cuts )  ## ATTENTION: insert Roo-Fit cuts here!
        
        ##  apply cuts (if any) for the  complex expressions 
        if fcuts :
            fcuts = [ '(%s)' % f for f in fcuts ]
            fcuts = ' && '.join ( fcuts )
            _vars = ds.get()
            with rooSilent ( ROOT.RooFit.ERROR  + 1  , True ) :
                with rootError( ROOT.kError + 1 ) :
                    ds1 = ROOT.RooDataSet ( dsID() , ds.title , ds , _vars , fcuts ) 
                    ds.clear()
                    del ds
                    ds = ds1
                    varsete = ds.get()
            if not silent :
                logger.info ( "make_dataset: dataset after (f)cuts\n%s" % ds.table ( prefix = "# " ) )
                
        nvars = ROOT.RooArgSet()
        for v in varset   : nvars.add ( v     )
        for v in formulas : nvars.add ( v.var )
        varset  = nvars 
        varsete = ds.get() 

    if formulas and not silent : 
        logger.info ( "make_dataset: dataset with expressions\n%s" % ds.table ( prefix = "# " ) ) 

    ##  remove all temporary variables  
    if len ( varset ) != len ( varsete ) :
        vs  = ROOT.RooArgSet()
        vrm = set ( )
        for v in ds.get() :
            if  v in varset : vs .add ( v      )
            else            : vrm.add ( v.name )
        if vrm and not silent :
            logger.info  ("make_dataset: temporary variables to be removed %s" % str ( tuple ( vrm) ) )
        ds1 = ds.reduce ( vs , '' )
        ds.clear()
        del ds
        ds = ds1
        
    if not silent : 
        skipped = 'Skipped:%d' % stat.skipped 
        skipped = '/' + attention ( skipped ) if stat.skipped else ''
        table   = ds.table () 
        report = 'make_dataset: Events Total:%d/Processed:%s%s CUTS:"%s" dataset\n%s' % (
            stat.total     ,
            stat.processed ,
            skipped        ,
            selection      ,
            ds             )
        logger.info (  report.replace ( '\n','\n# ' ) )
        
    return ds , stat 

ROOT.TTree.make_dataset_old = make_dataset_old




# =============================================================================
## Create RooDataset from the tree using Tree->Frame->Dataset transformation 
#  @code 
#  tree = ...
#  ds   = tree.make_dataset2 ( [ 'px , 'py' , 'pz' ] ) 
#  @endcode
def make_dataset ( tree              ,
                   variables         , ## varibales 
                   selection = ''    , ## TTree selection 
                   roo_cuts  = ''    , ## Roo-Fit selection  
                   name      = ''    , 
                   title     = ''    ,
                   silent    = False ) :
    """Create the dataset from the tree via intermediate Frame 
    >>> tree = ...
    >>> ds = tree.make_dataset ( [ 'px , 'py' , 'pz' ] ) 
    """
    if not title : title = 'Dataset from %s' % tree.GetName()
        
    if not DataSet_NEW_FILL :
        return make_dataset_old ( tree      = tree      ,
                                  variables = variables ,
                                  selection = selection ,
                                  roo_cuts  = roo_cuts  ,
                                  name      = name      ,
                                  title     = title     ,
                                  silent    = silent    )
    
    if not silent :
        logger.info  ("Use new efficient RooDataSet filling machinery")
        
    import ostap.trees.cuts
    import ostap.fitting.roofit

    from ostap.frames.frames import DataFrame, report_print, frame_columns   
    frame = DataFrame ( tree )

    total = len ( tree )
    
    if not silent :
        pb = frame.ProgressBar ( total )
                    
    columns = set ( tree.branches() ) | set ( tree.leaves() )

    scuts  = [] 
    limits = []

    vars = Variables ( variables ) 
    for v in vars :
        
        mn , mx = v.minmax
        
        if _minv < mn :
            lcut = "(%.16g <= %s)" % ( mn     , v.name )
            if silent : scuts.append  (   lcut )
            else      : limits.append ( ( lcut , 'RANGE(%s,low)'  % v.name ) )
            
        if _maxv > mx :
            hcut = "(%s <= %.16g)" % ( v.name , mx     )
            if silent : scuts.append ( hcut )
            else      : limits.append ( ( hcut , 'RANGE(%s,high)' % v.name ) )
            
        if v.name == v.formula and v.name in columns : continue

        ## define new variable 
        frame = frame.Define ( v.name , v.formula )

        columns |= set ( frame_columns ( frame  ) )
        
    ## define 'range/limit' cuts: 
    for c , f in limits :
        frame = frame.Filter ( c , f )
        
    if scuts :
        acut  = ' && '.join ( ( "(%s)" % c for c in scuts ) )
        frame = frame.Filter ( acut      , 'RANGES' )
        
    if selection :
        frame = frame.Filter ( selection , 'SELECTION' )

    if roo_cuts  :
        frame = frame.Filter ( roo_cuts  , 'ROO-CUTS'  )

    report   = frame.Report()
    
    varset   = ROOT.RooArgSet()
    for v in vars : varset.add ( v.var )
    
    name     = dsID() 
    title    = "Data set from DataFrame"
    
    rds = frame.Book (
        ROOT.std.move ( ROOT.RooDataSetHelper( name , title , varset ) ) ,
        tuple ( v.name for v in vars ) )
    
    ds   = rds.GetValue()


    selcuts  = None
    roocuts  = None
    lastcut  = None
    ncuts    = 0 
    for c in report :
        ncuts += 1 
        if   c.GetName() == 'SELECTION' : selcuts = c.GetAll () , c.GetPass ()
        elif c.GetName() == 'ROO-CUTS'  : roocuts = c.GetAll () , c.GetPass ()
        else                            : lastcut = c.GetAll () , c.GetPass () 

    if not silent and ncuts : 
        title = 'Tree->DataFrame->RooDataSet transformation'
        logger.info ( title + '\n%s' % report_print ( report , title , '# ' ) )
    
    if selcuts and roocuts :
        skipped   = total - selcuts [ 0 ] 
        processed = roocuts [ 1 ] 
    elif selcuts :
        skipped   = total - selcuts [ 0 ] 
        processed = selcuts [ 1 ] 
    elif roocuts :
        skipped   = total - roocuts [ 0 ] 
        processed =         roocuts [ 1 ] 
    elif lastcut :
        if 0 == lastcut [1]  :
            skipped   = total
            processed = 0 
        else :
            skipped   = total - lastcut[1] 
            processed = lastcut[1]
    else : 
        skipped   = 0 
        processed = total
        
        
    return ds , SelStat ( total , processed , skipped ) 

ROOT.TTree.make_dataset = make_dataset

        
# =============================================================================
## define the helper function for proper decoration of ROOT.TTree/TChain
#
# @code
# from ostap.fitting.pyselectors import SelectorWithVars 
# selector = SelectorWithVars ( ... ) 
# chain    = ...
# chain.fill_dataset2 ( selector )  ## NB: note lowercase "process" here !!!
# @endcode
#
# The module has been developed and used with great success in
# "Kali, framework for fine cailbration of LHCb Electormagnetic Calorimeter"
#
# @see Ostap::Selector 
# @see Ostap::Process 
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date   2010-04-30
def fill_dataset2 ( self              ,
                    selector          ,
                    nevents   = -1    ,
                    first     =  0    ,
                    shortcut  = True  ,
                    silent    = False ,
                    use_frame = 50000 ) :
    """ 'Process' the tree/chain with proper TPySelector :
    
    >>> from ostap.fitting.pyselectors import SelectorWithCVars     
    >>> selector = SelectorWithVars ( ... ) 
    >>> chain    = ...
    >>> chain.fill_dataset2 ( selector )  ## NB: note lowercase 'process' here !!!    
    """

    ## if use_frame and root_info < (6,15) :
    ##    logger.warning ( 'Processing via DataFrame is disabled for %s' % str ( root_info ) ) 
    ##    use_frame = False 

    ## process all events? 
    all = 0 == first and ( 0 > nevents or len ( self ) <= nevents )

    if all and shortcut and isinstance ( self , ROOT.TTree ) and isinstance ( selector , SelectorWithVars ) :
        
        if ( not selector.morecuts )  and \
               selector.trivial_vars  and \
               ( DataSet_NEW_FILL or selector.really_trivial ) :
            
            ## if selector.really_trivial and not selector.morecuts ) and \
            ##    ( not '[' in selector.selection ) : 
            
            if not silent : logger.info ( "Make try to use the *SHORTCUT*!" )
            variables     = selector.variables 
            ds , stat     = self.make_dataset ( variables = variables          ,
                                                selection = selector.selection ,
                                                roo_cuts  = selector.roo_cuts  ,
                                                silent    = silent             )
            selector.data = ds
            selector.stat = stat
            
            return ds , stat  
        
    # =========================================================================
    ## If the length is large and selection is not empty,
    #  try to pre-filter using RDataFrame machinery into  temporary file 
    #  @see ROOT.RDataFrame.Filter
    #  @see ROOT.RDataFrame.Snapshot
    #  It can be very efficient is selection/filtering cuts are harsh enough.    
    if all and 0 < use_frame and isinstance ( self , ROOT.TTree ) and use_frame <= len ( self ) :
        
        if isinstance ( selector , SelectorWithVars ) and selector.selection :

            if not silent : logger.info ( "Make try to use the intermediate DataFrame!" )

            selection = selector.selection
            
            from ostap.utils.utils    import ImplicitMT 
            from ostap.frames.frames  import DataFrame 
            
            total  = len ( self )

            frame  = DataFrame ( self , enable = True )

            frame_main = frame
            
            if not silent :
                pb = frame.ProgressBar ( len ( self ) )
                
            columns = set ( frame.columns() ) 
            
            vars   = list ( selector.variables )
            tvars  = [ v for v in vars if     v.trivial ] ## trivial vars 
            vars_  = [ v for v in vars if not v.trivial ] ## non-trivial vars
            nvars  = []
            scuts  = []
            rcuts  = selector.roo_cuts
            
            dvars  = [ v.name for v in vars if v.name != v.formula and v.name in columns ]
            assert not dvars, "Can't redefine existing variables: %s" % dvars  
            
            ranges = []

            ## loop over trivial vars :
            for v in tvars :
                
                mn , mx = v.minmax
                if  v.name == v.formula :  ## really trivial variable 
                    nvars.append ( v )
                else :                    ## almost trivial: formula exists  
                    newv  = Variable  ( v.var                       ,
                                        description = v.description ,
                                        vmin        = v.vmin        ,
                                        vmax        = v.vmax        ,
                                        accessor    = v.name        ) ## make it trivial!
                    nvars.append ( newv )
                    logger.debug  ( 'PROCESS: define %s as %s ' % ( v.name , v.formula ) )                    
                    frame = frame.Define ( v.name , v.formula )  ## define new variable  for the frame
                    
                if _minv < mn :
                    lcut = "(%.16g <= %s)" % ( mn     , v.name )
                    if silent : scuts.append  (   lcut )
                    else      : ranges.append ( ( lcut , 'RANGE(%s,low)'  % v.name ) )
                        
                if _maxv > mx :
                    hcut = "(%s <= %.16g)" % ( v.name , mx     )
                    if silent : scuts.append ( hcut )
                    else      : ranges.append ( ( hcut , 'RANGE(%s,high)' % v.name ) )

            frame  = frame.Filter ( selection , 'SELECTION' )

            for  c , f in ranges : 
                frame = frame.Filter ( c  , f )
                logger.debug  ( 'PROCESS: add cut %s ' % c )

            if scuts :
                acut  = ' && '.join ( ( "(%s)" % c for c in scuts ) )
                frame = frame.Filter ( acut , 'RANGES' )
                logger.debug  ( 'PROCESS: add cut %s ' % acut )

            if rcuts :
                frame = frame.Filter ( rcuts , 'ROO-CUTS' )
                
            from ostap.utils.cleanup import TempFile
            with TempFile ( suffix = '.root' , prefix = 'ostap-frame-' ) as tf :
                if not silent : logger.info ( 'Prepare snapshot/loop over the tree %s' % tf.filename )
                report   = frame . Report ()

                if vars_  :
                    ## If some variables are non-trivial: dump the whole tree 
                    snapshot = frame . Snapshot ( 'tree' , tf.filename )
                else      :
                    ## otherwise dump only needed variables 
                    bvars  = self.the_variables( [ v.formula for v in tvars ] )
                    avars  = list ( bvars ) + [ v.name for v in nvars if not v in bvars ] 
                    avars  = list ( set ( avars ) )
                    avars.sort() 
                    logger.debug ('PROCESS: dump only %s' % list ( avars ) )
                    from ostap.core.core import strings as _strings
                    avars    = _strings ( avars ) 
                    snapshot = frame . Snapshot ( 'tree' , tf.filename , avars )

                if not silent :
                    from ostap.frames.frames import report_print
                    title =  'Tree -> Frame -> Tree filter/transformation '

                    ## reduce in number of branches 
                    otv = len ( self.branches()     )
                    with ROOT.TFile.Open ( tf.filename  , 'read' ) as tt : 
                        ntv = len ( tt.tree.branches () )
                    eff = binomEff ( ntv , otv ) * 100
                    fmt_eff =  '%4.1g +- %-4.1g'
                    eff = fmt_eff % ( eff.value() , eff.error() )
                    eff = '{:^20}'.format ( eff ) 
                    rows = [ ( '#BRANCHES' , '%s' % otv , '%s' % ntv , eff , '' ) ]             
                    logger.info ( title + '\n%s' % report_print ( report , title , '# ' , more_rows = rows ) )
                    
                total_0 = -1
                for c in report :
                    if  total_0 < 0 : total_0 = c.GetAll()
                    else            : break
                        
                if not silent :
                    s = -1 
                    try :
                        s = os.path.getsize ( tf.filename ) 
                    except :
                        s = -1
                    if 0 < s : 
                        logger.info ( 'Snapshot at %s %.3g[MB]' % ( tf.filename , float ( s ) / 2**20 ) )  
                    else :
                        logger.info ( 'Snapshot at %s '         %   tf.filename ) 

                            
                import ostap.io.root_file

                nn = frame.Count().GetValue()
                if nn <= 0 :
                    logger.warning ( 'Selection frame result is empty' )
                elif not silent :
                    logger.info    ( 'Selected frame %d entries' % nn  )
                
                with ROOT.TFile.Open ( tf.filename  , 'read' ) as tt : 
                    tree         = tt.tree
                    if not silent :
                        import ostap.logger.table as  T
                        title = 'Filtered frame/tree for further processing'
                        t = tree.table ( prefix = '# ' )
                        ## t = str ( tree )
                        ## t = T.add_prefix ( t , '# ')
                        logger.info ( 'Filtered frame/tree for the futher processing:\n%s' % t )
                        
                    new_selector = SelectorWithVars ( nvars + vars_ ,
                                                      ''            , ## no cuts here! 
                                                      cuts     = selector.morecuts ,
                                                      roo_cuts = selector.roo_cuts , 
                                                      name     = selector.name     ,
                                                      fullname = selector.fullname ,
                                                      silence  = selector.silence  )
                    
                    if not silent : logger.info ( 'Redirect to (re)processing' )
                    
                    result       = fill_dataset2 ( tree                ,
                                                   new_selector        ,
                                                   nevents   = -1      ,
                                                   first     =  0      ,
                                                   shortcut  = True    ,
                                                   silent    = silent  ,
                                                   use_frame = -1      )
                    
                    selector.data           = new_selector.data
                    selector.stat.total     = total_0
                    selector.stat.processes = new_selector.stat.processed 
                    selector.stat.skiped    = new_selector.stat.skipped 
                    del new_selector

                    return result

    # =========================================================================
    ## Standard processing: no tricks, no shortcuts, ...
    # =========================================================================
    
    import ostap.fitting.roofit

    nevents = nevents if 0 <= nevents else ROOT.TChain.kMaxEntries
    if isinstance ( self , ROOT.TTree ) and isinstance ( selector , SelectorWithVars ) :
        if not silent : logger.info ( "No shortcuts&frame tricks posisble: use plain Selector" )
        args   =  () if all else ( nevents , first )
        result = Ostap.Utils.process ( self , selector , *args )
        if result < 0   : logger.error ("TTree::Process: result is %s" % result )
        elif not silent : logger.info  ("TTree::Process: result is %s" % result )  
        return selector.data, selector.stat
    
    if isinstance ( self , ROOT.TTree ) and isinstance ( selector , ROOT.TSelector ) :
        if not silent : logger.info ( "No shortcuts&frame tricks posisble: use plain Selector" )
        args =  () if all else ( nevents , first )
        result = Ostap.Utils.process ( self , selector , *args )
        if result < 0   : logger.error ("TTree::Process: result is %s" % result )
        elif not silent : logger.info  ("TTree::Process: result is %s" % result )  
        return result 

    ## RooDataSet is here:
    if not silent : logger.info ( "Process RooDataSet!" )
    
    assert not self.isWeighted() , \
           'Processing of the weighted dataset is not possible (yet?)'
    
    store = self.store()
    if store and store.tree() :
        tree = store.tree()
        return fill_dataset2 ( tree      ,
                               selector  ,
                               nevents   = nevents   ,
                               first     = first     ,
                               shortcut  = shortcut  ,
                               silent    = silent    ,
                               use_frame = use_frame )

    from ostap.fitting.roofit import useStorage
    
    with useStorage() :
        
        logger.info ('Prepare the cloned dataset with TTree-storage type')
        from ostap.core.core import dsID            
        cloned = self.Clone ( dsID() )        
        result = fill_dataset2 ( cloned    ,
                                 selector  ,
                                 nevents   = nevents   ,
                                 first     = first     ,
                                 shortcut  = shortcut  ,
                                 silent    = silent    ,
                                 use_frame = use_frame )
        cloned.reset()
        del cloned
        return result

fill_dataset2. __doc__ += '\n' + Ostap.Utils.process.__doc__

ROOT.TTree.fill_dataset2 = fill_dataset2

# =============================================================================
## Create RooDataset from the tree
#  @code 
#  tree = ...
#  ds   = tree.fill_dataset1 ( [ 'px , 'py' , 'pz' ] ) 
#  @endcode
def fill_dataset1 ( tree                 ,
                    variables            ,  ## list of variables 
                    selection    = ''    ,  ## TTree-cuts 
                    roo_cuts     = ''    ,  ## RooFit cuts 
                    cuts         = None  ,  ## python callable 
                    name         = ''    ,
                    title        = ''    ,
                    shortcut     = True  ,
                    use_frame    = 50000 , 
                    silent       = False ) :
    """Create the dataset from the tree
    >>> tree = ...
    >>> ds = tree.fill_dataset1 ( [ 'px , 'py' , 'pz' ] ) 
    """
    selector = SelectorWithVars ( variables ,
                                  selection = selection ,
                                  cuts      = cuts      , 
                                  roo_cuts  = roo_cuts  ,
                                  name      = name      , 
                                  silence   = silent    ) 
    tree.fill_dataset2 ( selector , silent = silent , shortcut  = shortcut , use_frame = use_frame )
    data = selector.data
    stat = selector.stat
    del selector
    
    return data , stat 

ROOT.TTree.fill_dataset1 = fill_dataset1

# =============================================================================
## Create RooDataset from the tree
#  @code 
#  tree = ...
#  ds   = tree.fill_dataset ( [ 'px , 'py' , 'pz' ] ) 
#  @endcode
def fill_dataset ( tree      ,
                   variables , **kwargs ) :
    """Create the dataset from the tree
    >>> tree = ...
    >>> ds = tree.fill_dataset ( [ 'px , 'py' , 'pz' ] )
    - see `ROOT.TTree.fill_dataset1`
    - see `ROOT.TTree.fill_dataset2`
    """
    if isinstance ( variables , SelectorWithVars ) :
        selector = variables 
        return fill_dataset2 ( tree , selector , **kwargs )
    return fill_dataset1 ( tree , variables , **kwargs )

fill_dataset.__doc__ += '\n' + fill_dataset2.__doc__
fill_dataset.__doc__ += '\n' + fill_dataset2.__doc__

ROOT.TTree.fill_dataset = fill_dataset


# =============================================================================
## define the helper function for proper decoration of ROOT.TTree/TChain
#
# @code
#
# from ostap.fitting.pyselectors import Selector
#
# class MySelector ( Selector ) :
#
#      def __init__ ( self ) :
#
#           return Selector ( self )
#
#      def Process  ( self , entry ) :
#          # == getting the next entry from the tree
#           if self.GetEntry ( entry ) <= 0 : return 0             ## RETURN 
#        
#           # == for more convenience
#           tree = self.tree 
#
#           # apply trivial "acceptance" cuts 
#           if not 2 <= tree.y   <=  4.5   : return 0               ## RETURN
#           if not 1 <= tree.pt  <=  10.0  : return 0               ## RETURN
#
#           return 1
# 
# selector = MySelector()
#
# chain = ...
# chain.process ( selector )  ## NB: note lowercase "process" here !!!
#
# @endcode
#
# The module has been developed and used with great success in
# "Kali, framework for fine cailbration of LHCb Electormagnetic Calorimeter"
#
# @see Ostap::Selector 
# @see Ostap::Process 
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date   2010-04-30
def _process_ ( self , selector , nevents = -1 , first = 0 , **kwargs ) :
    """ 'Process' the tree/chain with proper TPySelector :
    
    >>> from ostap.fitting.pyselectors import Selector    
    >>> class MySelector ( Selector ) : ... 
    ...
    >>> selector = MySelector()    
    >>> chain = ...
    >>> chain.process ( selector )  ## NB: note lowercase 'process' here !!!    
    """

    if isinstance ( self , ROOT.TTree ) and isinstance ( selector , SelectorWithVars ) :
        return fill_dataset2 ( self              ,
                               selector          ,
                               nevents = nevents ,
                               first   = first   , **kwargs ) 
    
    # =========================================================================
    ## Standard processing: no tricks, no shortcuts, ...
    # =========================================================================
    
    import ostap.fitting.roofit
    
    nevents = nevents if 0 <= nevents else ROOT.TChain.kMaxEntries
    args    =  () if all else ( nevents , first )
    
    return Ostap.Utils.process ( self , selector , *args ) 


_process_. __doc__ += '\n' + Ostap.Utils.process.__doc__

# =============================================================================
## finally: decorate TTree/TChain
for t in ( ROOT.TTree      ,
           ROOT.TChain     ,
           ROOT.RooAbsData ) : t.process  = _process_ 


_new_methods_ = [
    ROOT.TTree.make_dataset  ,
    ROOT.TTree.fill_dataset  ,
    ROOT.TTree.fill_dataset1 ,
    ROOT.TTree.fill_dataset2 ,    
    ROOT.TTree.process       ,
    ROOT.RooAbsData.process  , ## senseless :-( 
    ]

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
         
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================

