#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file selectors.py
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
    'SelectorWithVars' ,        ## Generic selctor to fill RooDataSet form TTree/TChain
    'Variable'         ,        ## helper class to define variable 
    'SelectorWithVarsCached'    ## Generic selector with cache   
)
# =============================================================================
import ROOT, cppyy, math, sys 
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.selectors' )
else                       : logger = getLogger ( __name__          )
# =============================================================================
from   ostap.core.core     import cpp, Ostap
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
        self.__selection = selection 
        Ostap.SelectorWithCuts.__init__ ( self , selection , None , self )
        if not self.cuts() :
            raise RuntimeError ("__init__:  Invalid Formula %s " % self.cuts() )
        
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
    
_maxv =  0.95 * sys.float_info.max
_minv = -0.95 * sys.float_info.max
##  is this a a trivial variable?
##_non_trivial = ' +-*/=%?()[]|&^'
##def trivial ( var ) :
##    for s in _non_trivial :
##        if 0 <= var.find ( s ) : return False
##    return True 
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

        if isinstance ( var , str ) :

            var = var.strip() 
            assert isinstance ( description , str ) , \
                   "Variable: illegal type for ``description''"     % ( description , type ( desctiption ) ) 
            assert isinstance ( vmin ,  ( int,long,float) ) , \
                   "Variable: illegal type for ``vmin'' %s/%s"      % ( vmin        , type ( vmin        ) )
            assert isinstance ( vmax ,  ( int,long,float) ) , \
                   "Variable: illegal type for ``vmax'' %s/%s"      % ( vmax        , type ( vmax        ) )
            assert vmin < vmax, \
                   "Variable: invalid ``minmax'' range (%g,%g)"     % ( vmin , vmax ) 
            
            description = description if description else "``%s''-variable" % var
            description = description.replace('\n',' ')

            ## create the variable
            var = ROOT.RooRealVar ( var , description , vmin , vmax ) 

        if accessor is None :
            varname  = var.GetName()
            accessor = varname
            
        self.__formula = False
        if isinstance  ( accessor , str ) :
            accessor = accessor.strip() 
            if accessor  : 
                from ostap.trees.funcs import FormulaFunc as FuncVar
                self.__formula  =  accessor 
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
    def accessor    ( self ) :
        """``accessor'' - the actual callable to get the value from TTree/TChain"""
        return self.__accessor
    @property
    def formula ( self ) :
        """``formula'' - formual for this variable (if applicable)?"""
        return self.__formula
    @property
    def trivial ( self ) :
        """``trivial'' - is this variable ``trivial''?"""
        return self.__formula and self.__formula == self.name

# =============================================================================
## is this expression corresponds to valid formula?
def valid_formula ( expression , varset ) :

    expression = expression.strip()
    if not expression : return True
    
    if isinstance  ( varset  , ROOT.RooArgSet ) :
        vlst = ROOT.RooArgList()
        for v in varset : vlst.add ( v  )
        result = valid_formula ( expression , vlst )
        del vlst
        return result

    assert isinstance ( varset  , ROOT.RooArgList ), 'nvaild type %s' % type (varset)
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

        self.__name = name 
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
        for v in variables :

            if   isinstance ( v , str              ) : v  = Variable (   v ) 
            elif isinstance ( v , ROOT.RooAbsReal  ) : v  = Variable (   v )
            elif isinstance ( v , ( tuple , list ) ) : v  = Variable (  *v )
            elif isinstance ( v , dict             ) : v  = Variable ( **v )

            assert isinstance  ( v , Variable ), 'Invalid variable %s/%s' % ( v , type(v) )

            self.__variables.append ( v     )
            self.__varset   .add    ( v.var )
            if not v.trivial : self.__triv_vars = False

        self.__variables = tuple( self.__variables ) 

        self.__triv_sel = valid_formula ( selection , self.__varset ) 
        triv_cuts      = not cuts
        
        self.__trivial = self.__triv_vars and self.__triv_sel and triv_cuts
        if not silence :
            self.__logger.info ( "Suitable for fast processing variables: %s, selection: %s, cuts: %s" % ( self.__triv_vars ,
                                                                                                           self.__triv_sel  ,
                                                                                                           triv_cuts      ) )
        if not self.__silence: 
            nl = 0
            dl = 0 
            for v in self.__variables :
                nl = max ( nl , len( v.name        ) ) 
                dl = max ( dl , len( v.description ) )                 
            dl = max ( dl , len ( 'Description' ) + 2 ) 
            nl = max ( nl , len ( 'Variable'    ) + 2 ) 
        
            line1    = '\n# | %%%ds | %%-%ds |         min / max         | Trivial? | ' % ( nl , dl ) 
            line2    = '\n# | %%%ds | %%-%ds | %%+11.3g / %%-+11.3g | %%-8s | ' % ( nl , dl )         
            the_line = 'Booked %d  variables:' % len ( self.variables ) 
            sep      = '\n# +%s+%s+%s+%s+' % ( (nl+2)*'-' , (dl+2)*'-' , 27*'-', 10*'-' )
            the_line += sep 
            the_line += line1 % ( 'Variable' , 'Description' )
            the_line += sep 
            for v in self.__variables :
                fmt = line2 % ( v.name        , 
                                v.description ,
                                v.minmax[0]   ,
                                v.minmax[1]   ,
                                v.trivial     )
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
        """``trivial_vars'' : are all variables ``trivial'' (sutable for fast-processing)?"""
        return self.__triv_vars
    
    @property
    def trivial_sel( self ) :
        """``trivial_sel'' : is the selection ``trivial'' (sutable for fast-processing)?"""
        return self.__triv_sel
    @property
    def trivial ( self ) :
        """``trivial'' : Are variables/selection/cuts ``trivial'' (sutable for fast-processing)?"""
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
        
        
    ## get the dataset 
    def dataset   ( self  ) :
        """ Get the data-set """ 
        return self.__data
 
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

        #
        ## apply cuts (if needed) 
        # 
        if self.__cuts and not self. __cuts ( bamboo )  : return 0 

        #
        ## loop over all variables
        # 
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
            from ostap.logger.logger import attention 
            skipped = 'Skipped:%d' % self.skipped 
            if self.skipped : skipped = attention ( skipped ) 
            self.__logger.info (
                'Selector(%s): Events Total:%d/Processed:%d/%s CUTS: "%s"' % (
                self.__name    ,
                self.total     ,
                self.processed ,
                skipped        , 
                self.cuts ()  ) )            
            self.__logger.info ( 'Selector(%s): dataset created:%s' %  ( self.__name ,  self.__data ) )
            
        if self.__data and not self.__silence :
            from ostap.logger.logger import attention 
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

                sep      = '# +%s+%s+%s+%s+%s+' % ( ( name_l       + 2 ) * '-' ,
                                                    ( desc_l       + 2 ) * '-' ,
                                                    ( mean_l+rms_l + 5 ) * '-' ,
                                                    ( min_l +max_l + 5 ) * '-' ,
                                                    ( skip_l       + 2 ) * '-' )
            fmt = '# | %%%ds | %%-%ds | %%%ds / %%-%ds | %%%ds / %%-%ds | %%-%ds | '  % (
                name_l ,
                desc_l ,
                mean_l ,
                rms_l  ,
                min_l  ,
                max_l  ,
                skip_l )
            
            report  = 'Dataset(%s) created: %s entries, %s variables ' %  ( self.__name ,
                                                                        len ( self.__data    ) ,
                                                                        len ( self.variables ) )
            
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
            for k,v in self.__skip.iteritems() : skip += v 
            self.__logger.warning("Selector(%s): empty dataset! Processed/Total/Skept %d/%d/%d"
                                  % ( self.__name  , self.__events , self.total , skip ) ) 
            
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
## Create the dataset from the tree
#  @code 
#  tree = ...
#  ds = tree.make_dataset ( [ 'px , 'py' , 'pz' ] ) 
#  @endcode
def _make_dataset_ ( tree , variables , selection = '' , name = '' , title = '' ) :
    """Create the dataset from the tree
    >>> tree = ...
    >>> ds = tree.make_dataset ( [ 'px , 'py' , 'pz' ] ) 
    """
    import ostap.trees.cuts
    import ostap.fitting.roofit

    cuts   = ROOT.TCut ( selection )
    varset = ROOT.RooArgSet()
    for v in variables :

        if   isinstance  ( v  , str              ) : v = Variable (   v )
        elif isinstance  ( v  , ( tuple , list ) ) : v = Variable (  *v )
        elif isinstance  ( v  , dict             ) : v = Variable ( **v )

        assert isinstance  ( v , Variable  )  , "Can't create Variable from %s/%s" % ( v , type ( v ) )
        assert v.trivial                      , "Variable %s is not ``trivial''"   % v.name   
        assert hasattr     ( tree , v.name )  , "Tree/Chain has no branch ``%s''"  % v.name

        varset.add  ( v.var )
        mn , mx = v.minmax
        if _minv < mn : cuts &= "%.16g <= %s"  % ( mn      , v.name   )
        if _maxv > mx : cuts &= "%s  <= %.16g" % ( v.name  , mx       )

    if not name :
        from ostap.core.core import dsID 
        name = '%s_%s' % ( dsID() , tree.GetName() )
    if not title : title = '%s/%s' % ( name , tree.GetTitle() )


    total     = len ( tree )
    processed = tree.statVar ( '1' , selection    ).nEntries()
    skept     = tree.statVar ( '1' , str ( cuts ) ).nEntries() 

    stat = total, processed , processed - skept

    from ostap.logger.utils import rooSilent, rootError  
    with rooSilent ( ROOT.RooFit.ERROR  , True ) :
        with rootError( ROOT.kWarning ) :
            ds = ROOT.RooDataSet ( name  , title , tree , varset , str( cuts ) )

    print cuts 
    ##ll = len(ds) 
    ##assert ll + stat[2] == stat[1], 'Invalid statustics; %s/%s' % ( ll , stat )
    ##print 'DS:' , len(ds), stat 
    return ds , stat 

ROOT.TTree.make_dataset = _make_dataset_ 
        
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
def _process_ ( self , selector , nevents = -1 , first = 0 , shortcut = True , silent = False  ) :
    """ ``Process'' the tree/chain with proper TPySelector :
    
    >>> from ostap.fitting.selectors import Selector    
    >>> class MySelector ( Selector ) : ... 
    ...
    >>> selector = MySelector()    
    >>> chain = ...
    >>> chain.process ( selector )  ## NB: note lowercase ``process'' here !!!    
    """

    print  'PROCESS(0)'
    ## process all events? 
    all = 0 == first and ( 0 > nevents or len ( self ) <= nevents )

    if all and shortcut and isinstance ( self , ROOT.TTree ) and isinstance ( selector , SelectorWithVars ) and selector.trivial :
        print  'PROCESS (2)'
        if not silent : logger.info ( "Make try to use the shortcut!" )
        ds , stat  = self.make_dataset( variables = selector.variables , selection = selector.selection )
        selector.data = ds
        selector.stat = stat 
        if not silent :
            from ostap.logger.logger import attention 
            skipped = 'Skipped:%d' % stat[2]
            if stat[2] : skipped = attention ( skipped )
            logger.info (
                'Selector(%s): Events Total:%d/Processed:%s/%s CUTS: "%s"\n# %s' % (
                selector.name    ,
                stat[0]          ,
                stat[1]          ,
                skipped          ,
                selector.cuts()  , ds ) )            
        return 1 
                            
    print  'PROCESS(3)'

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
        return _process_ ( tree , selector , events ,   silent = silent )

    from ostap.fitting.roofit import useStorage
    with useStorage() :
        logger.info ('Prepare the cloned dataset with TTree-storage type')
        from ostap.core.core import dsID            
        cloned = self.Clone ( dsID() )
    result = _process_ ( cloned , selector , events , silent = silent )
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

