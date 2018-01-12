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
# from Ostap.Selectors import Selector
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

# from Ostap.Selectors import Selector
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
'SelectorWithVarsCached'    ## Generic selector with cache   
)
# =============================================================================
import ROOT, cppyy, math
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.selectors' )
else                       : logger = getLogger ( __name__          )
# =============================================================================
from ostap.core.core import cpp, Ostap 
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
        Ostap.SelectorWithCuts.__init__ ( self , selection , None , self )
        if not self.cuts() :
            raise RuntimeError ("__init__:  Invalid Formula %s " % self.cuts() )
        
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
# from Ostap.PySelector import Selector
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
def _process_ ( self , selector , *args ) :
    """ ``Process'' the tree/chain with proper TPySelector :
    
    >>> from Ostap.Selectors import Selector    
    >>> class MySelector ( Selector ) : ... 
    ...
    >>> selector = MySelector()    
    >>> chain = ...
    >>> chain.process ( selector )  ## NB: note lowercase ``process'' here !!!    
    """
    import ostap.fitting.roofit
    return Ostap.Process.process ( self , selector , *args )

_process_. __doc__ += '\n' + Ostap.Process.process.__doc__


# =============================================================================
## finally: decorate TTree/TChain
for t in ( ROOT.TTree , ROOT.TChain ) : t.process  = _process_ 

# =============================================================================
## helper function to decode information about the variable
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2010-04-30
def makeEntry( var , *args ) :
    """Helper function to decode information about the variable
    >>> e = makeEntry( 'a' , 'description' , -1000 , 1000 ,  lambda s : s.pt/2 )
    >>> e = makeEntry( 'a' , 'description' , -1000 ,         lambda s : s.pt/2 )
    >>> e = makeEntry( 'a' , 'description'                   lambda s : s.pt/2 )
    >>> e = makeEntry( 'a' , 'description' , -1000 , 1000 )
    >>> e = makeEntry( 'a' , 'description' , -1000 )
    >>> e = makeEntry( 'a' , 'description' )
    """

    if   isinstance ( var , str ) :      ## just the name of variable

        vname = var            ## name

        a = list ( args )

        ## accessor function
        if a and callable ( a[-1] )    : vfun  = a.pop(-1)
        else                           : vfun  = lambda s : getattr ( s , vname )

        if   not a                     : vdesc = 'Variable: %s'   % vname
        elif isinstance ( a[0] , str ) : vdesc = a.pop(0)
        else                           : vdesc = 'Variable:   %s' % vname

        import sys
        _f  =  sys.float_info
        _fs = float,int,long

        if   not a                     : vmin = - _f.max
        elif isinstance ( a[0] , _fs ) : vmin =    float(a.pop(0))
        else : raise AttributeError (" Invalid min-value specification %s %s " % ( var , list ( args ) ) )

        if   not a                     : vmax =   _f.max
        elif isinstance ( a[0] , _fs ) : vmax =    float(a.pop(0))
        else : raise AttributeError (" Invalid max-value specification %s %s " % ( var , list ( args ) ) )

        if vmax <= vmin :
            raise AttributeError (" Invalid min/max specification %s %s "      % ( var , list ( args ) ) )

        var = ROOT.RooRealVar ( vname , vdesc , vmin , vmax )

    elif isinstance ( var , ROOT.RooRealVar ) : # variable itself

        vname = var.GetName  () ## name
        vdesc = var.GetTitle () ## description
        vmin  = var.getMin   () ## min-value
        vmax  = var.getMax   () ## max-value
        #

        a = list (  args )

        ## accessor function
        if a and callable ( a[-1] )    : vfun  = a.pop(-1)
        else                           : vfun  = lambda s : getattr ( s , vname )

    else :

        raise AttributeError ( 'Invalid variable description  %s %s' % (  var , list   ( args ) ) )

    ## unprocessed aruments ?
    if a : raise AttributeError ( "Can't recognize all arguments %s  %s "  % ( var , list ( args ) ) )

    ## finally the entry
    return var , vdesc , vmin , vmax , vfun

# =============================================================================
## helper class to decode/keep information about the variable in c
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date   2010-04-30
class VEntry(object) :
    """ Helper class to decode/keep infomration about the variable 
    """
    def __init__ ( self , var , *args ) :
        """Add declared variable to RooDataSet 
        """
        if   isinstance ( var , str ) :      ## just the name of variable   
            
            self.vname = var            ## name 
            self.vdesc = args[0]        ## description 
            self.vmin  = args[1]        ## min-value 
            self.vmax  = args[2]        ## max-value 
            #
            ## accessor function
            #
            if 3 < len ( args ) : self.vfun = args[3]
            else                : self.vfun = lambda s : getattr( s , self.vname )
            # 
            self.var = ROOT.RooRealVar ( self.vname , self.vdesc , self.vmin , self.vmax )

        elif isinstance ( var , ROOT.RooRealVar ) : # variable itself 

            self.vname = var.GetName  () ## name 
            self.vdesc = var.GetTitle () ## description
            self.vmin  = var.getMin   () ## min-value 
            self.vmax  = var.getMax   () ## max-value 
            #
            ## accessor function
            #
            if 0 < len ( args ) : self.vfun = args[0]
            else                : self.vfun = lambda s : getattr( s , self.vname )

        else :

            self._logger.error   ( 'Invalid variable description!' )
            raise AttributeError,  'Invalid variable description!'

        ## finally the entry
        self.entry = ( self.var , self.vdesc , self.vmin , self.vmax , self.vfun )
         
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
#   ( 'my_name1' , 'my_description1' , low       , high     )
#   ]
# 
#  ## get a variable 'my_name' from the tree/chain with accessor function, e.g. rescale it on-fligh
#  variables += [ 
#   #  name       descriptor           min-value , max-value , access function   
#   ( 'my_name2' , 'my_description2' , low       , high      , lambda s : s.my_name2/1000 )
#   ]
#   
#  ## get  less trivial expression
#  variables += [ 
#   #  name       descriptor           min-value , max-value , access function   
#   ( 'my_name3' , 'my_description3' , low       , high      , lambda s : s.var1+s.var2 ) 
#  ]
#
#  ## any function that gets TChain/Tree entry and evaluates to double.
#  #  e.g. it could be TMVAReader
#  def myvar ( chain ) : ....
#  variables += [ 
#   #  name       descriptor           min-value , max-value , access function   
#   ( 'my_name4' , 'my_description4' , low       , high      , myvar ) 
#  ]
#
#
#  ## add already booked variables: 
#  v5 = ROOT.RooRealVal( 'my_name5' )
#  variables += [  ( v5 , lambda s : s.var5 ) ]
#
#  
#  ## add already booked variables: 
#  v6 = ROOT.RooRealVal( 'my_name6' )
#  variables += [  ( v6                     ) ] ## get variable "my_name6"
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
#
#  @date   2014-03-02
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  - thanks to  Alexander BARANOV 
class SelectorWithVars(SelectorWithCuts) :
    """ Create and fill the basic dataset for RooFit
    # 
    #  variables = [ ... ]
    #
    #  ## add a variable 'my_name1' from the tree/chain 
    #  variables += [ 
    #   #  name       descriptor           min-value , max-value  
    #   ( 'my_name1' , 'my_description1' , low       , high     )
    #   ]
    # 
    #  ## get a variable 'my_name' from the tree/chain with accessor function,
    #  ## e.g. rescale it on-fligh
    #  variables += [ 
    #   #  name       descriptor           min-value , max-value , access function   
    #   ( 'my_name2' , 'my_description2' , low       , high      , lambda s : s.my_name2/1000 )
    #   ]
    #   
    #  ## get  less trivial expression
    #  variables += [ 
    #   #  name       descriptor           min-value , max-value , access function   
    #   ( 'my_name3' , 'my_description3' , low       , high      , lambda s : s.var1+s.var2 ) 
    #  ]
    #
    #  ## any function that gets Tchain/Tree and avaluated to double.
    #  #  e.g. it coudl be TMVAReader
    #  def myvar ( chain ) : ....
    #  variables += [ 
    #   #  name       descriptor           min-value , max-value , access function   
    #   ( 'my_name4' , 'my_description4' , low       , high      , myvar ) 
    #  ]
    #
    #
    #  ## add already booked variables: 
    #  v5 = ROOT.RooRealVal( 'my_name5' )
    #  variables += [  ( v5 , lambda s : s.var5 ) ]
    #
    #  
    #  ## add already booked variables: 
    #  v6 = ROOT.RooRealVal( 'my_name6' )
    #  variables += [  ( v6                     ) ] ## get variable 'my_name6'
    #
    #
    #  #
    #  ## finally create selector
    #  # 
    #  selector = SelectorWithVars (
    #             variables                             ,
    #             selection = ' chi2vx<30 && pt>2*GeV ' ,  ## filtering
    #            )
    #  chain = ...
    #  chain.process ( selector )
    #  dataset = selector.dataset
    # 
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

        #
        ## create the logger 
        #
        from  ostap.logger.logger           import getLogger
        self._logger = getLogger ( fullname ) 
        #

        self._events   = 0
        self._progress = None 
        self._total    = 1
        self._skip     = 0 
        self._silence  = silence

        #
        ## instantiate the base class
        # 
        SelectorWithCuts.__init__ ( self , selection ) ## initialize the base

        #
        ## keep the cuts
        # 
        self._cuts   = cuts

        #
        ## variables
        # 
        self.varset      = ROOT.RooArgSet()
        self._variables  = []
        
        #
        ## add the variables one by one 
        #
        for v in variables : self.addVariable ( *v )

        #
        ## Book dataset
        # 
        self.data    = ROOT.RooDataSet (
            ##
            name      ,
            fullname  , 
            ##
            self.varset
            )
        
        #
        ## it is still very puzzling for me: should this line be here at all??
        ROOT.SetOwnership ( self.data  , False )
        

    ## ## delete the selector, try to clear and delete the dataset 
    ## def __del__    ( self  )  :
    ##     #
    ##     if hasattr ( self , 'data' ) and self.data : 
    ##         self.data.Clear()
    ##         self.data.reset()
    ##         #
    ##         del self.data
        
    ## get the dataset 
    def dataset   ( self  ) :
        """ Get the data-set """ 
        return self.data
    
    ## the only one actually important method 
    def Process ( self, entry ):
        """
        Fills data set 
        """
        #
        ## == getting the next entry from the tree
        #
        if self.GetEntry ( entry ) <=  0 : return 0             ## RETURN 
        #
        
        if not self._progress and not self._silence :
            self._total =  self.fChain.GetEntries()
            self._logger.info ( "Processing TChain('%s') #entries: %d" % ( self.fChain.GetName() , self._total ) )
            ## decoration:
            from ostap.utils.progress_bar import ProgressBar
            self._progress = ProgressBar ( max_value = self._total   ,
                                           silent    = self._silence )
            
        if not self._silence :
            if 0 == self._events % 1000 or 0 == entry % 1000 : 
                self._progress.update_amount ( self.event () )
            
        self._events += 1
        #
        ## == for more convenience
        #
        bamboo = self.fChain

        #
        ## apply cuts (if needed) 
        # 
        if not self . _cuts ( bamboo )  : return 0 

        #
        ## loop over all varibales
        # 
        for v in self._variables :

            var    =  v[0]  ## variable 
            vmin   =  v[2]  ## min-value 
            vmax   =  v[3]  ## max-value 
            vfun   =  v[4]  ## accessor-function 

            value  = vfun ( bamboo )
            if not vmin <= value <= vmax :       ## MUST BE IN RANGE!
                self._skip += 1 
                return 0                         ## RETURN 

            var.setVal ( value ) 


        self.data .add ( self.varset )
        
        return 1 

    ## add declared variable to RooDataSet 
    def addVariable ( self , var , *args ) :
        """
        Add decared variable to RooDataSet 
        """
        
        entry = makeEntry  ( var , *args )

        self.varset.add        ( entry[0] )
        self._variables.append ( entry    )
        if not self._silence:
            self._logger.info ( 'Add variable name/desc/min/max %s/%s/%.2g/%.2g' % ( entry[0].GetName() ,
                                                                                     entry[1] ,
                                                                                     entry[2] ,
                                                                                     entry[3] ) )
     #
    def Terminate ( self  ) :
        #
        if self._progress :
            self._progress.end() 
        #
        if not self._silence : 
            self._logger.info (
                'Events Processed/Total/Skept %d/%d/%d\nCUTS: "%s"' % (
                self._events ,
                self._total  ,
                self._skip   , 
                self.cuts () ) ) 
            self.data.Print('v')
            
        if not len ( self.data ) :
            self._logger.warning("Empty dataset!")
        ##
        if 0 != self.GetAbort() :
            self._logger.error('Process has been aborted!')

        ##
        logger.debug('Terminate: DELETE all variables')
        del self._variables
        
    # 
    def Init    ( self, chain ) :
        # 
        if self._progress and not self._silence :
            self._progress.update_amount ( self.event () )
        #
        return SelectorWithCuts.Init ( self , chain ) 

    def Begin          ( self , tree = None ) :
        ## 
        if self._progress and not self._silence :
            self._progress.update_amount ( self.event () )
    #
    def SlaveBegin     ( self , tree        ) :
        # 
        if self._progress and not self._silence :
            self._progress.update_amount ( self.event () )
    #
    def Notify         ( self ) :
        #
        if self._progress and not self._silence :
            self._progress.update_amount ( self.event () )
            
    def SlaveTerminate ( self               ) :
        # 
        if self._progress and not self._silence :
            self._progress.update_amount ( self.event () )
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

        self.__selection = selection
        self.__cuts      = cuts
        self.__filelist  = files

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
        h += self._l_internals(self.__cuts)

        for var , vdesc , vmin , vmax , vfun in self._variables:
            h += var.GetName() + vdesc + str(vmin) + str(vmax)
            h += self._l_internals(vfun)

        h += str(self.__selection)

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

