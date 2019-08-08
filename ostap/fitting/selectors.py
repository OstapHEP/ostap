#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/fitting/selectors.py
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
# from ostap.fitting.selectors import Selector
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
#           tree=self.fChain
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
# ``Kali, framework for fine calibration of LHCb Electormagnetic Calorimeter''
#
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date   2010-04-30
# 
# =============================================================================
"""Helper module to fix a problems in communication of
TTree/TChain.Process and TPySelector.

In PyROOT some of original C++ methods are disable.
The module provides the 'recovery' for missing methods

# from ostap.fitting.selectors import Selector
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
# chain.process ( selector )  ## NB: note lowercase ``process'' here !!!
#

The module has been developed and used with great success in
``Kali, framework for fine calibration of LHCb Electormagnetic Calorimeter''

More complicated (and more useful) cases are covered by
SelectorWithCuts and SelectorWithVars classes 

"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2010-04-30"
__version__ = "$Revision$" 
# =============================================================================
__all__ = (
    'Selector'         ,        ## The ``fixed'' TPySelector
    'Selector2'        ,        ## The ``fixed'' TPySelector
    'SelectorWithCuts' ,        ## The ``fixed'' TPySelector with TTree-formula 
    'SelectorWithVars' ,        ## Generic selctor to fill RooDataSet from TTree/TChain
    'Variable'         ,        ## helper class to define variable 
    'SelectorWithVarsCached'    ## Generic selector with cache   
)
# =============================================================================
import ROOT, cppyy, math, sys 
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger, attention, allright
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.selectors' )
else                       : logger = getLogger ( __name__          )
# =============================================================================
from   ostap.core.core        import cpp, Ostap, items_loop 
from   ostap.core.ostap_types import num_types, string_types 
import ostap.fitting.roofit 
# =============================================================================
## C++ Selector 
Selector  = Ostap.Selector
# =============================================================================
## @class Selector2
#  Useful intermediate class for implementation of (py)selectors 
#  @see Ostap::Selector
#  @see TPySelector
#  @see TSelector
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-09
class Selector2 ( Ostap.Selector) :
    """Useful intermediate class for implementation of (py)selectors     
    """
    ## constructor 
    def __init__ ( self ) :
        """Standart constructor
        """
        ## initialize the base 
        Ostap.Selector.__init__ ( self , None , self )
        
    def Notify         ( self               ) : return True
    def Terminate      ( self               ) : pass 
    def SlaveTerminate ( self               ) : pass 
    def Begin          ( self , tree = None ) : pass
    def SlaveBegin     ( self , tree        ) : pass 
    def Init           ( self , tree        ) : pass
    ## the major method
    def Process        ( self , entry       ) :
        """The major method 
        """
        # load data 
        if self.GetEntry ( entry ) < 0 : return 0

        #
        ## put your code here 
        #
        
        return 1
    
# =============================================================================
## @class SelectorWithCuts
#  Efficient selector that runs only for ``good''-events  
#  @see Ostap::SelectorWithCuts
#  @see Ostap::Selector
#  @see TPySelector
#  @see TSelector
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-09
class SelectorWithCuts (Ostap.SelectorWithCuts) :
    """Efficient selector that runs only for ``good''-events  
    """
    ## constructor 
    def __init__ ( self , selection ) :
        """ Standart constructor
        """
        ## initialize the base
        self.__selection = str ( selection ).strip()  
        Ostap.SelectorWithCuts.__init__ ( self , self.selection , None , self )
        if self.cuts() : logger.info ( 'SelectorWithCuts: %s' % self.cuts() )
            
    @property
    def selection ( self ) :
        """``selection'' -  selection to be used to preprocess TTree/TChain"""
        return self.__selection
    
    def Notify         ( self               ) : return True
    def Terminate      ( self               ) : pass 
    def SlaveTerminate ( self               ) : pass 
    def Begin          ( self , tree = None ) : pass
    def SlaveBegin     ( self , tree        ) :
        if not self.ok() or not self.formula() :
            raise RuntimeError ("SlaveBegin:Invalid Formula %s " % self.cuts() )        
    def Init           ( self , tree        ) :
        if not self.ok() or not self.formula() :
            raise RuntimeError ("Init:      Invalid Formula %s " % self.cuts() )
        

    ## the major method
    def Process        ( self , entry       ) :
        """The major method 
        """
        # load data 
        if self.GetEntry ( entry ) < 0 : return 0

        #
        ## put your code here 
        #
        
        return 1
    
# =============================================================================
_maxv =  0.95 * sys.float_info.max
_minv = -0.95 * sys.float_info.max
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
               "Variable: invalid type for ``var''  %s/%s"      % ( var         , type (  var        ) ) 
        assert accessor is None or callable ( accessor ) or isinstance ( accessor , str )  , \
               "Variable: illegal type for ``accessor'' %s/%s"  % ( accessor    , type ( accessor    ) )

        self.__vmin = None
        self.__vmax = None
        
        if isinstance ( var , str ) :

            var = var.strip() 
            assert isinstance ( description , str ) , \
                   "Variable: illegal type for ``description''"     % ( description , type ( description ) ) 
            assert isinstance ( vmin , num_types ) , \
                   "Variable: illegal type for ``vmin'' %s/%s"      % ( vmin        , type ( vmin        ) )
            assert isinstance ( vmax , num_types ) , \
                   "Variable: illegal type for ``vmax'' %s/%s"      % ( vmax        , type ( vmax        ) )
            assert vmin < vmax, \
                   "Variable: invalid ``minmax'' range (%g,%g)"     % ( vmin , vmax ) 

            if    description                   : pass
            elif  isinstance ( accessor , str ) : description = accessor
            else                                : description = "``%s''-variable" % var
            
            description = description.replace('\n',' ')

            ## create the variable
            self.__vmin = vmin 
            self.__vmax = vmax 
            var = ROOT.RooRealVar ( var , description , vmin , vmax ) 

        if accessor is None :
            varname  = var.GetName()
            accessor = varname
            
        self.__formula = None
        if isinstance  ( accessor , str ) :
            accessor = accessor.strip() 
            if accessor  : 
                from ostap.trees.funcs import FormulaFunc as FuncVar
                self.__formula  = accessor 
                accessor = FuncVar ( accessor ) 

        assert callable ( accessor ), \
               'Invalid accessor function! %s/%s' % ( accessor , type ( accessor ) )
        
        self.__var         = var
        self.__minmax      = var.minmax()
        self.__accessor    = accessor
        
    @property
    def var         ( self ) :
        """``var'' - the variable itself/ROOT.RooRealVar"""
        return self.__var
    @property
    def name        ( self ) :
        """``name'' - the name of the variable"""
        return self.__var.GetName() 
    @property
    def description ( self ) :
        """``description'' - the variable description/title"""
        return self.__var.GetTitle()
    @property
    def minmax      (  self ) :
        """``minmax'' - the range for the variable """
        return self.__minmax    
    @property
    def vmin ( self ) :
        """``vmin'' : minimal value (if set) """
        if not self.__vmin is None : return self.__vmin
        mnmx = self.minmax
        if mnmx : return mnmx  [0]
        return None    
    @property
    def vmax ( self ) :
        """``vmax'' : maximal value (if set) """
        if not self.__vmax is None : return self.__vmax
        mnmx = self.minmax
        if mnmx : return mnmx [1]
        return None    
    @property
    def accessor    ( self ) :
        """``accessor'' - the actual callable to get the value from TTree/TChain"""
        return self.__accessor
    @property
    def formula ( self ) :
        """``formula'' - formula for this variable (when applicable)?"""
        return self.__formula
    
    @property
    def trivial ( self ) :
        """``trivial'' - is this variable ``trivial''?"""
        return self.__formula ## and self.__formula == self.name

    @property
    def really_trivial ( self ) :
        """``really_trivial'' - is it really trivia enough for RooFit?"""
        triv = self.trivial 
        return triv and  ( not '[' in triv ) and ( not ']' in triv )


# =============================================================================
## is this expression corresponds to a valid RooFit formula?
def valid_formula ( expression , varset ) :

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
    from ostap.logger.utils import rooSilent, rootError  
    with rooSilent ( ROOT.RooFit.FATAL + 1 , True ) :
        with rootError( ROOT.kError + 1 ) :
            _f = ROOT.RooFormulaVar( '' , expression , varset )
            fok = _f.ok ()
            del _f
            
    return fok
            
# ==============================================================================
## Define generic selector to fill RooDataSet from TChain
#
#  @code
# 
#  variables = [ ... ]
#
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
#
#  ## add already booked variables: 
#  v5 = ROOT.RooRealVal( 'my_name5' )
#  variables += [  Variable ( v5 , accessor = lambda s : s.var5 ) ]
#
#  
#  ## add already booked variables: 
#  v6 = ROOT.RooRealVal( 'my_name6' )
#  variables += [  Variable ( v6                     ) ] ## get variable "my_name6"
#
#
#  #
#  ## finally create selector
#  # 
#  selector = SelectorWithVars (
#             variables                             ,
#             selection = " chi2vx<30 && pt>2*GeV " ,  ## filtering
#            )
#  chain = ...
#  chain.process ( selector )
#  dataset = selector.dataset
# 
#  @endcode
#  @date   2014-03-02
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  - thanks to  Alexander BARANOV 
class SelectorWithVars(SelectorWithCuts) :
    """Create and fill the basic dataset for RooFit
    
    - Define the list of ``variables'' for selector:
    
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
    def __init__ ( self                           ,
                   variables                      ,  ## list of variables  
                   selection                      ,  ## Tree-selection 
                   cuts         = None            ,
                   name         = ''              ,
                   fullname     = ''              ,
                   silence      = False           ) :
        
        if not     name :
            from   ostap.core.core import dsID 
            name = dsID()
            
        if not fullname : fullname = name 

        self.__name     = name
        self.__fullname = fullname 
        #
        ## create the logger 
        #
        from ostap.logger.logger  import getLogger
        self.__logger = logger ## getLogger ( fullname ) 
        #
        self.__silence  = silence

        ##
        assert 0 < len(variables) , "Empty list of variables"
        #
        ## instantiate the base class
        # 
        SelectorWithCuts.__init__ ( self , selection ) ## initialize the base

        self.__cuts      = cuts
        self.__variables = [] 
        self.__varset    = ROOT.RooArgSet() 

        self.__triv_vars = True
        vvars = set() 
        for v in variables :

            vv = v 
            if   isinstance ( v , str              ) : vv = Variable (   v ) 
            elif isinstance ( v , ROOT.RooAbsReal  ) : vv = Variable (   v )
            elif isinstance ( v , ( tuple , list ) ) : vv = Variable (  *v )
            elif isinstance ( v , dict             ) : vv = Variable ( **v )
            elif isinstance ( v , Variable         ) : vv = v  

            assert isinstance  ( vv , Variable ), 'Invalid variable %s/%s' % ( vv , type ( vv ) )

            self.__variables.append ( vv     )
            self.__varset   .add    ( vv.var )
            #
            if   vv.trivial and vv.name == vv.formula : pass
            elif vv.really_trivial                    : pass
            elif vv.formula                           : pass
            else                                      :
                self.__triv_vars = False
            #
            vvars.add ( vv ) 
            
        self.__variables = tuple( self.__variables ) 

        self.__triv_sel  = valid_formula ( selection , self.__varset ) 
        triv_cuts        = not cuts
        
        self.__trivial = self.__triv_vars and self.__triv_sel and triv_cuts
        if not silence :
            tv = allright ( 'True' ) if self.__triv_vars else attention ( 'False' )
            ts = allright ( 'True' ) if self.__triv_sel  else attention ( 'False' )
            tc = allright ( 'True' ) if triv_cuts        else attention ( 'False' )
            self.__logger.info ( "Suitable for fast processing: variables:%s, selection:%s, py-cuts:%s" % ( tv , ts , tc ) )
            
        if not self.__silence: 
            nl = 0
            dl = 0 
            for v in self.__variables :
                nl = max ( nl , len( v.name        ) ) 
                dl = max ( dl , len( v.description ) )                 
            dl = max ( dl , len ( 'Description' ) + 2 ) 
            nl = max ( nl , len ( 'Variable'    ) + 2 ) 
        
            line1    = '\n# | %%%ds | %%-%ds |         min / max         | Trivial? | ' % ( nl , dl ) 
            line2    = '\n# | %%%ds | %%-%ds | %%+11.3g / %%-+11.3g | %%s | ' % ( nl , dl )         
            the_line = 'Booked %d  variables:' % len ( self.variables ) 
            sep      = '\n# +%s+%s+%s+%s+' % ( (nl+2)*'-' , (dl+2)*'-' , 27*'-', 10*'-' )
            the_line += sep 
            the_line += line1 % ( 'Variable' , 'Description' )
            the_line += sep
            for v in self.__variables :
                trivial = allright ('True') + 4* ' ' if v.trivial else attention ( 'False' ) + 3 * ' '
                    
                fmt = line2 % ( v.name        , 
                                v.description ,
                                v.minmax[0]   ,
                                v.minmax[1]   ,
                                trivial       )
                the_line += fmt
            the_line += sep 
            self.__logger.info ( the_line )
            
        ## Book dataset
        self.__data = ROOT.RooDataSet (
            ##
            self.name ,
            fullname  , 
            ##
            self.__varset
            )
        
        #
        ## it is still very puzzling for me: should this line be here at all??
        ROOT.SetOwnership ( self.__data  , False )
        
        self.__progress = None 
        from collections import defaultdict
        self.__skip     = defaultdict(int)
        self.__notifier = None
        self.__stat    = [ 0 , 0 , 0 ] 

    @property 
    def name ( self ) :
        """``name''  - the name of selector/dataset"""
        return self.__name
    
    @property 
    def fullname ( self ) :
        """``fullname''  - the fullname of selector/dataset"""
        return self.__fullname 
    
    @property 
    def data ( self ) :
        """``data''  - the dataset"""
        return self.__data
    @data.setter
    def data ( self , dataset ) :
        assert isinstance ( dataset , ROOT.RooAbsData ), \
               "Incorrect type of data %s/%s " % ( dataset ,   type ( dataset ) )
        self.__logger.debug ("Selector(%s), add dataset %s" % (  self.__name , dataset ) )
        self.__data = dataset 

    @property 
    def variables ( self ) :
        """``variables'' - the list/tuple of variables (cleared in Terminate)"""
        return self.__variables

    @property
    def varset ( self ) :
        """``varset'' : the structure of RooDataSet"""
        return self.__varset
    
    @property
    def morecuts ( self ) :
        """``morecuts'' -   additional cust to be applied in selection"""
        return self.__cuts

    @property
    def trivial_vars( self ) :
        """``trivial_vars'' : are all variables ``trivial'' (suitable for fast-processing)?"""
        return self.__triv_vars

    @property
    def really_trivial ( self ) :
        """``really_trivial'' : is a set of variables really trivial (for RooFit)?"""
        for v in self.variables :
            if not v.really_trivial : return False
        return True
    
    @property
    def trivial_sel( self ) :
        """``trivial_sel'' : is the selection ``trivial'' (suitable for fast-processing)?"""
        return self.__triv_sel
    
    @property
    def trivial ( self ) :
        """``trivial'' : Are variables/selection/cuts ``trivial'' (suitable for fast-processing)?"""
        return self.__trivial

    @property
    def skip ( self ) :
        """``skip'' : dictionary of skept entries"""
        return self.__skip
    
    @property
    def skipped ( self ) :
        """``skipped'' : total number of skept entries"""
        return self.__stat[2]
    
    @property
    def processed  ( self ) :
        """``processed'' : number of processeed events (after cuts)"""
        return self.__stat[1]
    
    @property
    def total  ( self ) :
        """``total'' : total number of processeed events (before cuts)"""
        return self.__stat[0]

    @property
    def stat ( self ) :
        """``stat'' : Total/processed/skipped events"""
        return tuple(self.__stat)
    @stat.setter
    def stat ( self , value  ) :
        assert 2<= len(value), 'Invalid "value":%s' % str ( value )
        self.__stat[0] = value[0]
        self.__stat[1] = value[1]
        self.__stat[2] = value[2]

    ## silent processing ?
    @property
    def silence (  self ) : return self.__silence
    
    ## get the dataset 
    def dataset   ( self  ) :
        """ Get the data-set """ 
        return self.__data

    # =========================================================================
    ## the only one actually important method 
    def Process ( self, entry ):
        """ Fill data set 
        """
        #
        ## == getting the next entry from the tree
        #
        if self.GetEntry ( entry ) <=  0 : return 0             ## RETURN 
        #
        
        if not self.__progress and not self.__silence :
            self.__stat[0] =  self.fChain.GetEntries()
            self.__logger.info ( "Selector(%s): processing TChain('%s') #entries: %d" % ( self.name , self.fChain.GetName() , self.total ) )
            ## decoration:
            from ostap.utils.progress_bar import ProgressBar
            self.__progress = ProgressBar ( max_value = self.total     ,
                                            silent    = self.__silence )
            
        if not self.__silence :
            if 0 == self.processed % 1000 or 0 == entry % 1000 or 0 == self.event() % 1000 : 
                self.__progress.update_amount ( self.event () )
                
        self.__stat[1] += 1
        
        #
        ## == for more convenience
        #
        bamboo = self.fChain

        return  self.fill ( bamboo )

    # =========================================================================
    ## fill it! 
    def fill ( self , bamboo ) :
        """The  actual processing for the given ``bamboo''
        Note that   this method is independent on TTree/TChain and can be used directy
        One just needs to  ensure that:
        - 'accessor functions' for the variables and 'cuts' agree with the type of ``bamboo''
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
                self.__stat[2]      += 1     ## SKIP EVENT 
                return 0                     ## RETURN 

            var.setVal ( value ) 


        self.__data .add ( self.__varset )
        
        return 1 

    # =========================================================================
    ## ``callable'' interface 
    def __call__ ( self ,  entry ) :
        """``callable'' interface to Selector
        """
        return self.fill ( entry ) 

    ## termination 
    def Terminate ( self  ) :
        #
        if self.__progress :
            self.__progress.end() 
        #
        ## Aborted? 
        if   0 != self.GetAbort() :
            self.__logger.fatal('Selector(%s): process has been aborted!' % self.__name )

            self.__data = None 
            del self.__varset
            del self.__variables
            self.__varset     =  ()
            self.__variables  =  ()
            
            return  ## RETURN

        ##get total number of input events from base class 
        self.__stat[0] = self.event()
        
        if not self.__silence :
            skipped = 'Skipped:%d' % self.skipped
            skipped = '/' + attention ( skipped ) if self.skipped else ''
            cuts    = allright ( '"%s"' % self.cuts () ) if self.trivial_sel else attention ( '"%s"'  % self.cuts() ) 
            self.__logger.info (
                'Selector(%s): Events Total:%d/Processed:%d%s CUTS: %s' % (
                self.__name    ,
                self.total     ,
                self.processed ,
                skipped        , 
                cuts           ) )            
            self.__logger.info ( 'Selector(%s): dataset created:%s' %  ( self.__name ,  self.__data ) )
            
        if self.__data and not self.__silence :
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
                name_l = max ( name_l , len ( v[0] ) )
                desc_l = max ( desc_l , len ( v[1] ) )
                mean_l = max ( mean_l , len ( v[2] ) )
                rms_l  = max ( rms_l  , len ( v[3] ) )
                min_l  = max ( min_l  , len ( v[4] ) )
                max_l  = max ( max_l  , len ( v[5] ) )
                skip_l = max ( skip_l , len ( v[6] ) )

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
                                                            
            header  = fmt % ( 'Variable'    ,
                              'Description' ,
                              'mean' ,
                              'rms'  ,
                              'min'  ,
                              'max'  ,
                              'skip' )
            report += '\n' + sep
            report += '\n' + header
            report += '\n' + sep            
            for v in vars :
                line    =  fmt % ( v[0] , v[1] , v[2] , v[3] , v[4] , v[5] , attention ( v[6] ) )
                report += '\n' + line  
            report += '\n' + sep
            self.__logger.info ( report ) 
        
        if not len ( self.__data ) :
            skip = 0
            for k,v in items_loop ( self.__skip ) : skip += v 
            self.__logger.warning("Selector(%s): empty dataset! Total:%s/Processed:%s/Skipped:%d"
                                  % ( self.__name  , self.total , self.processed , skip ) ) 
            
        ## attention: delete these

        del self.__varset
        del self.__variables
        
        self.__varset     =  ()
        self.__variables  =  ()

    def Init    ( self, chain ) :
        # 
        result = SelectorWithCuts.Init ( self , chain ) 

        if self.__progress and not self.__silence :
            self.__progress.update_amount ( self.event () )
        #
        return result 

    def Begin          ( self , tree = None ) :
        ## 
        result = SelectorWithCuts.Begin ( self , tree )

        if self.__progress and not self.__silence :
            self.__progress.update_amount ( self.event () )

        return result
    #
    def SlaveBegin     ( self , tree        ) :
        #
        result = SelectorWithCuts.SlaveBegin ( self , tree )
        #
        if self.__progress and not self.__silence :
            self.__progress.update_amount ( self.event () )
        #
        self.__stat[0] =  tree.GetEntries()
        #
        if self.__notifier :
            self.__notifier.exit()
            del self.__notifier
        
        self.__notifier = Ostap.Utils.Notifier( tree )
        for v in self.__variables :
            if isinstance ( v.accessor , ROOT.TObject ) :
                self.__notifier.add  ( v.accessor ) 
        
        return result 
    #
    def Notify         ( self ) :
        #
        result = SelectorWithCuts.Notify ( self )
        if self.__progress and not self.__silence :
            self.__progress.update_amount ( self.event () )

        return result 
            
    def SlaveTerminate ( self               ) :
        # 
        result = SelectorWithCuts.SlaveTerminate ( self )

        if self.__progress and not self.__silence :
            self.__progress.update_amount ( self.event () )

        if self.__notifier :
            self.__notifier.exit()
            self.__notifier = None  
            
        return result 


# =============================================================================
import os
from   ostap.core.workdir import workdir
from   ostap.io.zipshelve import ZipShelf 
# =============================================================================

# ==============================================================================
## Generic selector which loads already loaded datasets from cache
#
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
                   cuts         = None            ,
                   name         = ''              ,
                   fullname     = ''              ) : 

        SelectorWithVars.__init__(self, variables, selection, cuts, name, fullname)

        # Try load from cache
        self._loaded_from_cache = False
        self._cachepath = self._compute_cache_name()

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
## Create thew dataset from the tree
#  @code 
#  tree = ...
#  ds = tree.make_dataset ( [ 'px , 'py' , 'pz' ] ) 
#  @endcode
def make_dataset ( tree , variables , selection = '' , name = '' , title = '' , silent = False ) :
    """Create the dataset from the tree
    >>> tree = ...
    >>> ds = tree.make_dataset ( [ 'px , 'py' , 'pz' ] ) 
    """
    import ostap.trees.cuts
    import ostap.fitting.roofit

    varset   = ROOT.RooArgSet()
    vars     = set()

    formulas  = []
    
    selection = str ( selection ) if isinstance ( selection , ROOT.TCut ) else selection  
    selection = selection.strip() if isinstance ( selection , str       ) else selection 
    
    cuts = [ selection ] if selection else [] 
    for v in variables :

        if   isinstance  ( v , str              ) : vv = Variable (   v )
        elif isinstance  ( v , ROOT.RooRealVar  ) : vv = Variable (   v )
        elif isinstance  ( v , ( tuple , list ) ) : vv = Variable (  *v )
        elif isinstance  ( v , dict             ) : vv = Variable ( **v )
        elif isinstance  ( v , Variable         ) : vv = v 
        else :
            logger.error("Do not know how to treat the variable %s/%s, skip it" % ( v , type ( v ) ) )
            continue

        if vv.trivial and vv.name == vv.formula : 
            
            assert hasattr  ( tree , vv.name ) , "Tree/Chain has no branch ``%s''" % vv.name
            
            varset.add  ( vv.var )
            vars.add ( vv )
            
        elif vv.formula :
            
            formulas.append ( vv )
            continue
        
        else :
            
            logger.error("Do not know how to treat the variable %s, skip it" % vv.name )
            continue 
        
        mn , mx = vv.minmax
        if _minv < mn : cuts.append ( "(%.16g <= %s)" % ( mn      , vv.name ) ) 
        if _maxv > mx : cuts.append ( "(%s <= %.16g)" % ( vv.name , mx      ) )

    ## 
    cuts = ROOT.TCut(' && '.join(cuts) ) if cuts else ROOT.TCut()

    ## extended varset
    stor    = set() 
    varsete = ROOT.RooArgSet()
    for v in varset : varsete.add ( v )

    expressions = [ f.formula for f in formulas ]
    if selection : expressions.append ( selection ) 
        
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
        for lname in lvars :
            if not lname in varsete :
                v = Variable ( lname )
                varsete.add  ( v.var )
                stor.add ( v )
                
    if not name :
        from ostap.core.core import dsID 
        name = '%s_%s' % ( dsID() , tree.GetName() )
    if not title : title = '%s/%s' % ( name , tree.GetTitle() )

    total     = len ( tree )
    processed = tree.statVar ( '1' , selection    ).nEntries()
    skipped   = tree.statVar ( '1' , str ( cuts ) ).nEntries() 

    stat = total, processed , processed - skipped

    from ostap.logger.utils import rooSilent, rootError
    if not silent : logger.info ( "Start to fill the dataset") 
    with rooSilent ( ROOT.RooFit.ERROR  , True ) :
        with rootError( ROOT.kWarning ) :
            ds = ROOT.RooDataSet ( name  , title , tree , varsete , str( cuts ) )
            varsete = ds.get()
            
    ## add complex expressions 
    if formulas :
        # a
        vset   = ds.get()
        vlst = ROOT.RooArgList()
        for v in vset : vlst.add ( v )

        fcols = ROOT.RooArgList() 
        
        ffs   = []
        fcuts = [] 
        for f in formulas :            
            fv = ROOT.RooFormulaVar ( f.name , f.description , f.formula , vlst )
            assert fv.ok() , 'Invalid formula: %s' % f.formula 
            ffs.append ( fv )
            fcols.add  ( fv )
            mn , mx = f.minmax            
            if _minv < mn : fcuts.append ( "(%.16g <= %s)" % ( mn      , fv.name ) )
            if _maxv > mx : fcuts.append ( "(%s <= %.16g)" % ( fv.name , mx      ) )

        ds.addColumns ( fcols )
        ##  apply cuts (if any) for the  complex expressions 
        if fcuts :
            fcuts = [ '(%s)' % f for f in fcuts ]
            fcuts = ' && '.join ( fcuts )
            _vars = ds.get()
            ds1 = ROOT.RooDataSet ( dsID() , ds.title , ds , _vars , fcuts ) 
            ds.clear()
            del ds
            ds = ds1
            varsete = ds.get()
            
        nvars = ROOT.RooArgSet()
        for v in varset   : nvars.add ( v     )
        for v in formulas : nvars.add ( v.var )
        varset  = nvars 
        varsete = ds.get() 
        
    ##  remove all temporary variables  
    if len ( varset ) != len ( varsete ) :
        ds1 = ds.reduce ( varset )
        ds.clear()
        del ds
        ds = ds1
        
    if not silent : 
        skipped = 'Skipped:%d' % stat[2]
        skipped = '/' + attention ( skipped ) if stat[2] else '' 
        logger.info (
            'make_dataset: Events Total:%d/Processed:%s%s CUTS: "%s"\n# %s' % (
            stat[0] ,
            stat[1] ,
            skipped ,
            selection  , ds ) )            
        
    return ds , stat 

ROOT.TTree.make_dataset = make_dataset
        

# =============================================================================
## define the helper function for proper decoration of ROOT.TTree/TChain
#
# @code
#
# from ostap.fitting.selectors import Selector
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
#           tree=self.fChain
#
#           # apply trivial "acceptance" cuts 
#           if not 2 <= tree.y   <=  4.5   : return 0               ## RETURN
#           if not 1 <= tree.pt  <=  10.0  : return 0               ## RETURN
#
#           return 1
#
# 
# selector = MySelector()
#
# chain = ...
# chain.process ( selector )  ## NB: note lowercase "process" here !!!
#
# @endcode
#
# The module has been developed and used with great success in
# ``Kali, framework for fine cailbration of LHCb Electormagnetic Calorimeter''
#
# @see Ostap::Selector 
# @see Ostap::Process 
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date   2010-04-30
#
def _process_ ( self , selector , nevents = -1 , first = 0 , shortcut = True , silent = False  , use_frame = 50000 ) :
    """ ``Process'' the tree/chain with proper TPySelector :
    
    >>> from ostap.fitting.selectors import Selector    
    >>> class MySelector ( Selector ) : ... 
    ...
    >>> selector = MySelector()    
    >>> chain = ...
    >>> chain.process ( selector )  ## NB: note lowercase ``process'' here !!!    
    """

    ## process all events? 
    all = 0 == first and ( 0 > nevents or len ( self ) <= nevents )

    if all and shortcut and isinstance ( self , ROOT.TTree ) and isinstance ( selector , SelectorWithVars ) :
        if selector.really_trivial and not selector.morecuts :
            if not silent : logger.info ( "Make try to use the SHORTCUT!" )
            ds , stat  = self.make_dataset ( variables = selector.variables , selection = selector.selection , silent = silent )
            selector.data = ds
            selector.stat = stat 
            return 1
        
    # =========================================================================
    ## If the length is large and selection is not empty,
    #  try to pre-filter using RDataFrame machinery into  temporary file 
    #  @see ROOT.RDataFrame.Filter
    #  @see ROOT.RDataFrame.Snapshot
    #  It can be very efficient is selection/filtering cuts are harsh enough.    
    if all and 0 < use_frame and isinstance ( self , ROOT.TTree ) and use_frame <= len ( self ) :    
        if isinstance ( selector , SelectorWithVars ) and selector.selection :
            if not silent : logger.info ( "Make try to use intermediate FRAME!" )
            
            from ostap.utils.utils   import ImplicitMT 
            import ostap.frames.frames
            
            total  = len ( self )
            
            frame  = Ostap.DataFrame ( self , enable = True )
            
            frame  = frame.Filter ( selector.selection , 'SELECTION' )
            
            vars   = list ( selector.variables )
            tvars  = [ v for v in vars if     v.trivial ] ## trivial vars 
            vars_  = [ v for v in vars if not v.trivial ] ## non-trivial vars
            nvars  = []
            scuts  = []

            for v in tvars :
                mn , mx = v.minmax
                if v.name == v.formula :                    
                    nvars.append ( v )                    
                else :                    
                    newv  = Variable  ( v.var ,
                                        description = v.description ,
                                        vmin        = v.vmin        ,
                                        vmax        = v.vmax        ,
                                        accessor    = v.name        ) ## make it trivial!
                    nvars.append ( newv )
                    logger.debug  ( 'PROCESS: define %s as %s ' % ( v.name , v.formula ) )
                    frame = frame.Define ( v.name , v.formula )
                    
                if _minv < mn :
                    lcut = "(%.16g <= %s)" % ( mn     , v.name )
                    if silent : scuts.append ( lcut )
                    else      : frame = frame.Filter ( lcut , 'RANGE(%s,low)'  % v.name )
                    logger.debug  ( 'PROCESS: add cut %s ' % lcut )
                    
                if _maxv > mx :
                    hcut = "(%s <= %.16g)" % ( v.name , mx     )
                    if silent : scuts.append ( hcut )
                    else      : frame = frame.Filter ( hcut , 'RANGE(%s,high)' % v.name )
                    logger.debug  ( 'PROCESS: add cut %s ' % hcut )

            if scuts :
                acut = ' && '.join ( ( "(%s)" % c for c in scuts ) )
                frame = frame.Filter ( acut , 'RANGES' )
                logger.debug  ( 'PROCESS: add cut %s ' % acut )
                
            from ostap.utils.cleanup import TempFile
            with TempFile ( suffix = '.root' , prefix = 'frame_' ) as tf :
                if not silent : logger.info ( 'Prepare snapshot/loop over the tree   %s' % tf.filename )
                report   = frame . Report ()

                if vars_  :
                    ## is some variables are non-trivial: dump the whole tree 
                    snapshot = frame . Snapshot ( 'tree' , tf.filename )
                else      :
                    ## otherwise sump only needed variables 
                    bvars = self.the_variables( [ v.formula for v in tvars ] )
                    avars = set ( [ v.name for v in nvars ] )
                    for b in bvars : avars.add ( b )
                    avars = tuple ( avars ) 
                    from ostap.core.core import strings as _strings 
                    logger.debug  ('PROCESS: dump only %s' % list ( avars ) )
                    snapshot = frame . Snapshot ( 'tree' , tf.filename , _strings ( *avars ) )

                if not silent :
                    from ostap.frames.frames import report_prnt
                    txt = report_prnt ( report , 'Tree -> Frame -> Tree filter-transformation: ' )
                    logger.info ( txt )
                    
                if not silent : logger.info ( 'Write %s' % tf.filename  ) 
                import ostap.io.root_file 
                with ROOT.TFile.Open ( tf.filename  , 'read' ) as tt : 
                    tree         = tt.tree
                    new_selector = SelectorWithVars ( nvars + vars_ ,
                                                      ''            , ## no cuts here! 
                                                      cuts     = selector.morecuts ,
                                                      name     = selector.name     ,
                                                      fullname = selector.fullname ,
                                                      silence  = selector.silence  )
                    
                    if not silent : logger.info ( 'redirect to other processing' )
                    
                    result       = _process_ ( tree , new_selector ,
                                               nevents   = -1      ,
                                               first     =  0      ,
                                               shortcut  = True    ,
                                               silent    = silent  ,
                                               use_frame = -1      )
                    
                    selector.data = new_selector.data
                    selector.stat = new_selector.stat
                    del new_selector
                    
                    return result
                
    import ostap.fitting.roofit
    
    nevents = nevents if 0 <= nevents else ROOT.TChain.kMaxEntries
    if   isinstance ( self , ROOT.TTree ) :
        args =  () if all else ( nevents , first)        
        return Ostap.Process.process ( self , selector , *args ) 

    ## RooDataSet is here:
    
    assert not self.isWeighted() , \
           'Processing of the weighted dataset is not possible (yet?)'
    
    store = self.store()
    if store and store.tree() :
        tree = store.tree()
        return _process_ ( tree     ,  selector  ,
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
        
        result = _process_ ( cloned , selector ,
                             nevents   = nevents   ,
                             first     = first     ,
                             shortcut  = shortcut  ,
                             silent    = silent    ,
                             use_frame = use_frame )
        cloned.reset()
        del cloned
        return result
        
_process_. __doc__ += '\n' + Ostap.Process.process.__doc__

# =============================================================================
## finally: decorate TTree/TChain
for t in ( ROOT.TTree , ROOT.TChain , ROOT.RooAbsData ) : t.process  = _process_ 


# =============================================================================
if '__main__' == __name__ :
         
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================

