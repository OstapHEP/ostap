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
from   ostap.logger.logger    import getLogger
from   ostap.logger.colorized import attention, allright
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.selectors' )
else                       : logger = getLogger ( __name__          )
# =============================================================================
from   ostap.core.core        import cpp, Ostap, items_loop 
from   ostap.core.ostap_types import num_types, string_types, integer_types  
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
    def __init__ ( self , selection , silence = False ) :
        """ Standart constructor
        """
        self.__silence = silence

        ## initialize the base
        self.__selection = str ( selection ).strip()  
        Ostap.SelectorWithCuts.__init__ ( self , self.selection , 0 , self )
        
        if self.cuts () and not self.silence :
            logger.info ( 'SelectorWithCuts: %s' % self.cuts() )
            
    @property
    def silence ( self ) :
        """``silence''  : silent processing?"""
        return self.__silence 
    
    @property
    def selection ( self ) :
        """``selection'' -  selection to be used to preprocess TTree/TChain"""
        return self.__selection
    
    def Notify         ( self               ) : return True
    def Terminate      ( self               ) : pass 
    def SlaveTerminate ( self               ) : pass 
    def Begin          ( self , tree = None ) : pass
    def SlaveBegin     ( self , tree        ) :
        if not self.ok () :
            raise RuntimeError ( "SlaveBegin:Invalid Formula %s " % self.cuts () )        
    def Init           ( self , tree        ) :
        if not self.ok () :
            raise RuntimeError ( "Init:      Invalid Formula %s " % self.cuts () )
        
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
                from ostap.trees.funcs import FuncFormula as FuncVar
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
            _f  = Ostap.FormulaVar( expression , varset , False )
            fok = _f.ok ()
            del _f
            
    return fok
            
# ==============================================================================
## @class SelStat
#  Helper class to keep the statististics for SelectorWithVars 
class SelStat(object) :
    """Helper class to keep the statististics for SelectorWithVars
    """
    def __init__ ( self ,  total = 0 , processed = 0 , skipped = 0 ) :
        self.__total     = total
        self.__processed = processed 
        self.__skipped   = skipped  

    @property
    def total ( self ) :
        """``total''   : total number of events"""
        return self.__total
    @total.setter
    def total ( self , value ) :
        assert isinstance ( value , integer_types ) and 0 <= value ,\
               "Invalid value for ``total''"
        self.__total = value
            
    @property
    def processed ( self ) :
        """``processed''   : number of processed events"""
        return self.__processed
    @processed.setter
    def processed ( self , value ) :
        assert isinstance ( value , integer_types ) and 0 <= value ,\
               "Invalid value for ``processed''"
        self.__processed = value
        
    @property
    def skipped ( self ) :
        """``skipped''   : number of skipped events (e.g. due to variabel ranges)"""
        return self.__skipped
    @skipped.setter
    def skipped ( self , value ) :
        assert isinstance ( value , integer_types ) and 0 <= value ,\
               "Invalid value for ``skipped''"
        self.__skipped = value

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
        assert 0 < len(variables) , "Empty list of variables"
        #
        ## instantiate the base class
        # 
        SelectorWithCuts.__init__ ( self , selection , silence ) ## initialize the base

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
        if not self.silence :
            tv = allright ( 'True' ) if self.__triv_vars else attention ( 'False' )
            ts = allright ( 'True' ) if self.__triv_sel  else attention ( 'False' )
            tc = allright ( 'True' ) if triv_cuts        else attention ( 'False' )
            self.__logger.info ( "Suitable for fast processing: variables:%s, selection:%s, py-cuts:%s" % ( tv , ts , tc ) )
            
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
                title = '%s("%s","%s") %d variables' % ( 'RooDataSet', self.name , fullname , len ( self.__varset ) )
            else :
                title = '%s("%s") %d variables'      % ( 'RooDataSet', self.name ,            len ( self.__varset ) )
                        
            import ostap.logger.table as T
            t  = T.table (  table_data , title , '# ' )
            self.__logger.info ( "Booked dataset: %s\n%s" % ( title , t ) ) 

            
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
        self.__stat     = SelStat() 

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
        return  not '[' in self.selection and not ']' in self.selection 

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
        """``skipped'' : total number of skipped entries"""
        return self.stat.skipped
    
    @property
    def processed  ( self ) :
        """``processed'' : number of processeed events (after cuts)"""
        return self.stat.processed 
    
    @property
    def total  ( self ) :
        """``total'' : total number of processeed events (before cuts)"""
        return self.stat.total

    @property
    def stat ( self ) :
        """``stat'' : Total/processed/skipped events"""
        return self.__stat    
    @stat.setter
    def stat ( self , value  ) :
        assert isinstance ( value , SelStat ) , 'Invalid "value":%s' % str ( value )
        self.__stat = value 

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
        
        if not self.__progress and not self.silence :
            self.stat.total =  self.fChain.GetEntries()
            self.__logger.info ( "Selector(%s): processing TChain('%s') #entries: %d" % ( self.name , self.fChain.GetName() , self.total ) )
            ## decoration:
            from ostap.utils.progress_bar import ProgressBar
            self.__progress = ProgressBar ( max_value = self.total   ,
                                            silent    = self.silence )
            
        if not self.silence :
            if 0 == self.processed % 1000 or 0 == entry % 1000 or 0 == self.event() % 1000 : 
                self.__progress.update_amount ( self.event () )
                
        self.stat.processed += 1
        
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
                self.stat.skipped   += 1     ## SKIP EVENT 
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
        kw [ 'cuts'      ] = self.morecuts
        kw [ 'silence'   ] = self.silence

        kw.update ( kwargs )
        return SelectorWithVars ( **kw ) 


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
        self.stat.total = self.event()
        
        if not self.silence :
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
            self.__logger.info ( title + '\n' + t )
            
            
        if not self.__data or not len ( self.__data ) :
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

        if self.__progress and not self.silence :
            self.__progress.update_amount ( self.event () )
        #
        return result 

    def Begin          ( self , tree = None ) :
        ## 
        result = SelectorWithCuts.Begin ( self , tree )

        if self.__progress and not self.silence :
            self.__progress.update_amount ( self.event () )

        return result
    #
    def SlaveBegin     ( self , tree        ) :
        #
        result = SelectorWithCuts.SlaveBegin ( self , tree )
        #
        if self.__progress and not self.silence :
            self.__progress.update_amount ( self.event () )
        #
        self.stat.total =  tree.GetEntries()
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
        if self.__progress and not self.silence :
            self.__progress.update_amount ( self.event () )

        return result 
            
    def SlaveTerminate ( self               ) :
        # 
        result = SelectorWithCuts.SlaveTerminate ( self )

        if self.__progress and not self.silence :
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
## @class    SelectorWithVarsCached
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
## Create RooDataset from the tree
#  @code 
#  tree = ...
#  ds   = tree.make_dataset ( [ 'px , 'py' , 'pz' ] ) 
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
        name = dsID () 
    if not title and tree.GetName() != tree.GetTitle  :
        title = tree.GetTitle ()

    total     = len ( tree )
    processed = tree.statVar ( '1' , selection    ).nEntries()
    skipped   = tree.statVar ( '1' , str ( cuts ) ).nEntries() 

    stat = SelStat ( total , processed , processed - skipped )

    from ostap.logger.utils import rooSilent, rootError
    from ostap.utils.timing import timing
        
    if not silent : logger.info ( "Start to fill the dataset") 
    with rooSilent ( ROOT.RooFit.ERROR  , True ) :
        with rootError( ROOT.kWarning ) :
            ds = ROOT.RooDataSet ( name  , title , tree , varsete , str( cuts ) )
            varsete = ds.get()
                
    ## add complex expressions 
    if formulas :
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

ROOT.TTree.make_dataset = make_dataset
        
# =============================================================================
## Create RooDataset from the tree
#  @code 
#  tree = ...
#  ds   = tree.fill_dataset2 ( [ 'px , 'py' , 'pz' ] ) 
#  @endcode
def fill_dataset ( tree                 ,
                   variables            ,
                   selection    = ''    ,
                   name         = ''    ,
                   title        = ''    ,
                   shortcut     = True  ,
                   use_frame    = 50000 , 
                   silent       = False ) :
    """Create the dataset from the tree
    >>> tree = ...
    >>> ds = tree.fill_dataset ( [ 'px , 'py' , 'pz' ] ) 
    """
    selector = SelectorWithVars ( variables , selection , silence = silent ) 
    tree.process ( selector , silent = silent , shortcut  = shortcut , use_frame = use_frame )
    data = selector.data
    stat = selector.stat
    del selector 
    return data , stat 

ROOT.TTree.fill_dataset = fill_dataset

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
        if selector.really_trivial and not selector.morecuts and not '[' in selector.selection : 
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
            
            if not silent : logger.info ( "Make try to use the intermediate TFrame!" )
            
            from ostap.utils.utils   import ImplicitMT 
            import ostap.frames.frames
            
            total  = len ( self )

            frame  = Ostap.DataFrame ( self , enable = True )
            if not silent :
                pb = frame.ProgressBar ( len ( self ) )
            
            vars   = list ( selector.variables )
            tvars  = [ v for v in vars if     v.trivial ] ## trivial vars 
            vars_  = [ v for v in vars if not v.trivial ] ## non-trivial vars
            nvars  = []
            scuts  = []

            ranges = []

            ## loop over trivial vars :
            for v in tvars :
                
                mn , mx = v.minmax
                if v.name == v.formula :  ## really trivial variable 
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

            frame  = frame.Filter ( selector.selection , 'SELECTION' )            

            for  c , f in ranges : 
                frame = frame.Filter ( c  , f )
                logger.debug  ( 'PROCESS: add cut %s ' % c )

            if scuts :
                acut = ' && '.join ( ( "(%s)" % c for c in scuts ) )
                frame = frame.Filter ( acut , 'RANGES' )
                logger.debug  ( 'PROCESS: add cut %s ' % acut )

            from ostap.utils.cleanup import TempFile
            with TempFile ( suffix = '.root' , prefix = 'frame-' ) as tf :
                if not silent : logger.info ( 'Prepare snapshot/loop over the tree %s' % tf.filename )
                report   = frame . Report ()

                if vars_  :
                    ## is some variables are non-trivial: dump the whole tree 
                    snapshot = frame . Snapshot ( 'tree' , tf.filename )
                else      :
                    ## otherwise dump only needed variables 
                    bvars = self.the_variables( [ v.formula for v in tvars ] )
                    avars = list ( bvars ) + [ v.name for v in nvars if not v in bvars ] 
                    from ostap.core.core import strings as _strings 
                    logger.debug ('PROCESS: dump only %s' % list ( avars ) )
                    snapshot = frame . Snapshot ( 'tree' , tf.filename , _strings ( *avars ) )

                if not silent :
                    from ostap.frames.frames import report_print
                    title =  'Tree -> Frame -> Tree filter/transformation '
                    logger.info ( title + '\n%s' % report_print ( report , title , '# ') )
                    
                total_0 = -1
                for c in report :
                    if  total_0 < 0 : total_0 = c.GetAll()
                    else            : break
                        
                if not silent : logger.info ( 'Write %s' % tf.filename  ) 
                import ostap.io.root_file

                assert frame.Count().GetValue()>0 , 'Selection result is empty'
                
                with ROOT.TFile.Open ( tf.filename  , 'read' ) as tt : 
                    tree         = tt.tree
                    if not silent :
                        import ostap.logger.table as  T
                        t = str ( tree )
                        t = T.add_prefix ( t , '# ')
                        logger.info ( 'Filtered frame/tree for the futher processing:\n%s' % t )
                        
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

                    selector.stat.total     = total_0
                    selector.stat.processes = new_selector.stat.processed 
                    selector.stat.skiped    = new_selector.stat.skipped 
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

