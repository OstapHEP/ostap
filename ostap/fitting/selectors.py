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
def _process_ ( self , selector , events = -1 , silent = False  ) :
    """ ``Process'' the tree/chain with proper TPySelector :
    
    >>> from ostap.fitting.selectors import Selector    
    >>> class MySelector ( Selector ) : ... 
    ...
    >>> selector = MySelector()    
    >>> chain = ...
    >>> chain.process ( selector )  ## NB: note lowercase ``process'' here !!!    
    """
    import ostap.fitting.roofit 
    
    if   isinstance ( self , ROOT.TTree ) :


        if 0 < events : return Ostap.Process.process ( self , selector , events  )
        else          : return Ostap.Process.process ( self , selector )

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
        assert accessor is None or callable ( accessor ) or isinstance ( str )  , \
               "Variable: illegal type for ``accessor'' %s/%s"  % ( accessor    , type ( accessor    ) )

        if isinstance ( var , str ) :

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
            
        if isinstance  ( accessor , str ) and accessor :
            from ostap.trees.funcs import FormulaFunc as FuncVar
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
                   cuts         = lambda s : True ,
                   name         = ''              ,
                   fullname     = ''              ,
                   silence      = False           ) :
        
        if not     name :
            from   ostap.core.core import dsID 
            name = dsID()
            
        if not fullname : fullname = "%s/%s " % ( __name__ , name )

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

        #
        ## keep the cuts
        # 
        self.__cuts       = cuts
        
        #
        ## variables
        # 
        self.__varset     = ROOT.RooArgSet()
        self.__variables  = []
        
        #
        ## add the variables one by one 
        #
        for v in variables :
            self.__addVariable ( v )

        self.__variables = tuple ( self.__variables ) 

            
        if not self.__silence: 
            nl = 0
            dl = 0 
            for v in self.__variables :
                nl = max ( nl , len( v.name        ) ) 
                dl = max ( dl , len( v.description ) )                 
            dl = max ( dl , len ( 'Description' ) + 2 ) 
            nl = max ( nl , len ( 'Variable'    ) + 2 ) 
        
            line1    = '\n# | %%%ds | %%-%ds |         min / max         |' % ( nl , dl ) 
            line2    = '\n# | %%%ds | %%-%ds | %%+11.3g / %%-+11.3g |' % ( nl , dl )         
            the_line = 'Booked %d  variables:' % len ( self.variables ) 
            sep      = '\n# +%s+%s+%s+' % ( (nl+2)*'-' , (dl+2)*'-' , 27*'-' )
            the_line += sep 
            the_line += line1 % ( 'Variable' , 'Description' )
            the_line += sep 
            for v in self.__variables :
                fmt = line2 % ( v.name        , 
                                v.description ,
                                v.minmax[0]   ,
                                v.minmax[1]   )
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
        
        self.__events   = 0
        self.__progress = None 
        self.__total    = 1
        from collections import defaultdict
        self.__skip     = defaultdict(int)

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
    def morecuts ( self ) :
        """``morecuts'' -   additional cust ot be applied in selection"""
        return self.__cuts
    
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
            self.__total =  self.fChain.GetEntries()
            self.__logger.info ( "Selector(%s): processing TChain('%s') #entries: %d" % ( self.name , self.fChain.GetName() , self.__total ) )
            ## decoration:
            from ostap.utils.progress_bar import ProgressBar
            self.__progress = ProgressBar ( max_value = self.__total   ,
                                            silent    = self.__silence )
            
        if not self.__silence :
            if 0 == self.__events % 1000 or 0 == entry % 1000 : 
                self.__progress.update_amount ( self.event () )
                
        self.__events += 1
        
        #
        ## == for more convenience
        #
        bamboo = self.fChain

        #
        ## apply cuts (if needed) 
        # 
        if not self. __cuts ( bamboo )  : return 0 

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
                return 0                     ## RETURN 

            var.setVal ( value ) 


        self.__data .add ( self.__varset )
        
        return 1 

    ## add declared variable to RooDataSet 
    def __addVariable ( self , variable  ) :
        """Add declared variable to RooDataSet 
        """
        if   isinstance ( variable , str          ) : variable = Variable (   variable )
        elif isinstance ( variable , (tuple,list) ) : variable = Variable (  *variable )
        elif isinstance ( variable , dict         ) : variable = Variable ( **variable )

        assert isinstance ( variable , Variable ), \
               "Invalid type of ``variable'' %s/%s" % (  variable , type ( variable ) )
        
        self.__varset.add       ( variable.var ) 
        self.__variables.append ( variable     )
        if not self.__silence: 
            self.__logger.verbose ( "Selector(%s): add variable name/desc/(min,max): ``%s''/``%s''/(%.3g,%.3g)" % (
                self.name             , 
                variable.name         ,
                variable.description  , 
                variable.minmax[0]    , 
                variable.minmax[1]    ) )
            
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
            
        if not self.__silence :
            skip = 0
            for k,v in self.__skip.iteritems() : skip += v             
            self.__logger.info (
                'Selector(%s): Events Processed/Total/Skept %d/%d/%d CUTS: "%s"' % (
                self.__name    ,
                self.__events  ,
                self.__total   ,
                skip           , 
                self.cuts () ) )            
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
                         ('%-.5g' % mean.error() ).strip() ,   ## 3 
                         ('%.5g'  % rms          ).strip() ,   ## 4 
                         ('%+.5g' % mnmx[0]      ).strip() ,   ## 5
                         ('%+.5g' % mnmx[1]      ).strip() ,   ## 6 
                         ('%-d'   % self.__skip[ v.name]   ) ) ## 7  
                vars.append ( r )

            vars.sort()
            
            name_l  = len ( 'Variable'    ) + 2 
            desc_l  = len ( 'Description' ) + 2 
            mean_l  = len ( 'mean' ) + 2 
            err_l   = len ( 'err'  ) + 2 
            rms_l   = len ( 'rms'  ) + 2
            min_l   = len ( 'min'  ) + 2 
            max_l   = len ( 'max'  ) + 2 
            skip_l  = len ( 'Skip' ) 
            for v in vars :
                name_l = max ( name_l , len ( v[0] ) )
                desc_l = max ( desc_l , len ( v[1] ) )
                mean_l = max ( mean_l , len ( v[2] ) )
                err_l  = max ( err_l  , len ( v[3] ) )
                rms_l  = max ( rms_l  , len ( v[4] ) )
                min_l  = max ( min_l  , len ( v[5] ) )
                max_l  = max ( max_l  , len ( v[6] ) )
                skip_l = max ( skip_l , len ( v[7] ) )

            sep      = '# +%s+%s+%s+%s+%s+%s+' % ( ( name_l       + 2 ) * '-' ,
                                                   ( desc_l       + 2 ) * '-' ,
                                                   ( mean_l+err_l + 6 ) * '-' ,
                                                   ( rms_l        + 2 ) * '-' ,
                                                   ( min_l +max_l + 5 ) * '-' ,
                                                   ( skip_l       + 2 ) * '-' )
            fmt = '# | %%%ds | %%-%ds | %%%ds +- %%-%ds | %%%ds | %%%ds / %%-%ds | %%-%ds | '  % (
                name_l ,
                desc_l ,
                mean_l ,
                err_l  ,
                rms_l  ,
                min_l  ,
                max_l  ,
                skip_l )
            
            report  = 'Dataset created: %s entries, %s variables ' %  ( len ( self.__data    ) ,
                                                                        len ( self.variables ) )
            
            header  = fmt % ( 'Variable'    ,
                              'Description' ,
                              'mean' ,
                              'err'  ,
                              'rms'  ,
                              'min'  ,
                              'max'  ,
                              'skip' )
            report += '\n' + sep
            report += '\n' + header
            report += '\n' + sep            
            for v in vars :
                line    =  fmt % ( v[0] , v[1] , v[2] , v[3] , v[4] , v[5] , v[6] , v[7] )
                report += '\n' + line  
            report += '\n' + sep
            self.__logger.info ( report ) 
        
        if not len ( self.__data ) :
            skip = 0
            for k,v in self.__skip.iteritems() : skip += v 
            self.__logger.warning("Selector(%s): empty dataset! Processed/Total/Skept %d/%d/%d"
                                  % ( self.__name  , self.__events , self.__total , skip ) ) 
            
        ## attention: delete these

        del self.__varset
        del self.__variables
        
        self.__varset     =  ()
        self.__variables  =  ()

    def Init    ( self, chain ) :
        # 
        if self.__progress and not self.__silence :
            self.__progress.update_amount ( self.event () )
        #
        return SelectorWithCuts.Init ( self , chain ) 

    def Begin          ( self , tree = None ) :
        ## 
        if self.__progress and not self.__silence :
            self.__progress.update_amount ( self.event () )
    #
    def SlaveBegin     ( self , tree        ) :
        # 
        if self.__progress and not self.__silence :
            self.__progress.update_amount ( self.event () )
    #
    def Notify         ( self ) :
        #
        if self.__progress and not self.__silence :
            self.__progress.update_amount ( self.event () )
            
    def SlaveTerminate ( self               ) :
        # 
        if self.__progress and not self.__silence :
            self.__progress.update_amount ( self.event () )
        #

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
                   cuts         = lambda s : True ,
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
if '__main__' == __name__ :
    
     
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================

