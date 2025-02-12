#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/trees/funcs.py
#  Module with helper base classes ...
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Decoration of Tree/Chain objects for efficient use in python"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'FuncTree'          , ## helper base class for 'TTree-function'
    'FuncData'          , ## helper base class for 'RooAbsData-function'
    'PyTreeFunction'    , ## 'TTree-function' that uses python function/callable 
    'PyDataFunction'    , ## 'Data-function' that uses python function/callable
    "pyfun_tree"        , ## ditto, but as fnuction 
    "pyfun_data"        , ## ditto, but as fnuction 
    'FuncFormula'       , ## simple wrapper over TTreeFormula/Ostap::Formula  
    'FuncRooFormula'    , ## simple wrapper over RooFormulaVar  
    'FuncTH1'           , ## TH1-based Tree-function 
    'FuncTH2'           , ## TH2-based Tree-function 
    'FuncTH3'           , ## TH3-based Tree-function 
    ) 
# =============================================================================
from   ostap.core.meta_info    import old_PyROOT, root_info  
from   ostap.core.core         import Ostap, valid_pointer
from   ostap.utils.basic       import typename, prntrf 
import ostap.trees.treereduce 
import ROOT, abc 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.trees.funcs' )
else                       : logger = getLogger( __name__             )
# =============================================================================
## @class FuncTree
#  Helper class to implement "TTree-function"
#  @see Ostap::Functions::PyFuncTree
class FuncTree(Ostap.Functions.PyFuncTree) :
    """ Helper class to implement ``TTree-function'' in python 
    """
    def __init__ ( self , tree = None , clone = None ) :
        ## initialize the base class
        if tree is None : tree = ROOT.nullptr
        ##
        assert isinstance ( tree  , ROOT.TTree ) or not tree , \
            "FuncTree: inbalid 'tree'  type %s" % typename ( tree ) 
        assert not clone or isinstance ( clone , FuncTree ) , \
            "FuncTree: Invalid 'clone' type:%s" % typename ( clone )
        ## 
        if clone        : super (FuncTree,self).__init__ ( clone )
        elif old_PyROOT : super (FuncTree,self).__init__ ( self , tree )
        else            : super (FuncTree,self).__init__ (        tree )
        
    @property
    def the_tree ( self ) :
        """``the_tree'' : the actual pointer to the ROOT.TTree"""
        return self.tree ()

    ## the main method, "double-abstract" one
    @abc.abstractmethod 
    def evaluate ( self ) : raise NotImplementedError("%s:  method `evaluate' is not implemented!" % typename ( self ) ) 
    
    ## tricky clone method,  "double-abstract" one
    @abc.abstractmethod 
    def clone ( self , name = ""  ) : raise NotImplementedError("%s:  method `clone' is not implemented!" % typename ( self ) ) 

    ## we need to serialize the instance, "double-abstrac" one 
    @abc.abstractmethod
    def __reduce__ ( self ) :  raise NotImplementedError("%s:  method `__reduce__' is not implemented!" % typename ( self ) ) 

    
    def __str__  ( self ) : return typename ( self )
    def __repr__ ( self ) : return typename ( self )

# =============================================================================
## @class FuncData
#  Helper class to implement "RooAbsData-function"
#  @see Ostap::Functions::PyFuncData
class FuncData(Ostap.Functions.PyFuncData) :
    """ Helper class to implement ``TTree-function''
    """
    def __init__ ( self , data = None , clone = None ) :
        ## initialize the base class
        if  data is None : data = ROOT.nullptr
        ## 
        assert isinstance ( data , ROOT.RooAbsData ) or not data , \
            "FuncData: invalid 'data'  type:%s" % typename ( data ) 
        assert not clone  or isinstance ( clone , FuncData ) , \
            "FuncData: Invalid 'clone' type: %s" % typename ( clone )
        ## 
        if clone        : super (FuncData,self).__init__ ( clone )            
        elif old_PyROOT : super (FuncData,self).__init__ ( self , data )
        else            : super (FuncData,self).__init__ (        data )

    def __str__  ( self ) : return typename ( self )
    def __repr__ ( self ) : return typename ( self )
    
    @property
    def the_data ( self ) :
        """``the_data'' : the actual pointer to the ROOT.RooAbsData"""
        return self.data ()
    
    ## the main method, "double-abstract" one
    @abc.abstractmethod 
    def evaluate ( self ) : raise NotImplementedError("%s:  method `evaluate' is not implemented!" % typename ( self ) ) 
    
    ## tricky clone method,  "double-abstract" one
    @abc.abstractmethod 
    def clone ( self , name = ""  ) : raise NotImplementedError("%s:  method `clone' is not implemented!" % typename ( self ) ) 

    ## we need to serialize the instance, "double-abstrac" one 
    @abc.abstractmethod
    def __reduce__ ( self ) :  raise NotImplementedError("%s:  method `__reduce__' is not implemented!" % typename ( self ) ) 

    def __str__  ( self ) : return typename ( self )+'QUQU2'
    def __repr__ ( self ) : return typename ( self )+'QUQU2'

# =============================================================================
## @class PyTreeFunction
#  The concrete implementation on of <code>PyFunTree</code>, <code>FuncTree</code>
#  that delegates the evaluation to provided function/callable object
#  @see ostap.trees.funcs.FunTree
#  @see Ostap.Functions.PyFuncTree
#  @see Ostap.IFuncTree 
class PyTreeFunction(FuncTree) :
    """ The concrete implementation on of PyFunTree, FuncTree,
    that delegates the evaluation to provided function/callable object
    - see ostap.trees.funcs.FunTree
    - see Ostap.Functions.PyFuncTree
    - see Ostap.IFuncTree 
    """
    store  = set()
    
    def __init__ ( self , the_function , tree = None , clone = None ) :
        """ Constructor from the function/callable and (optional) tree"""
        if tree is None : tree = ROOT.nullptr
        ##
        assert not clone or isinstance ( clone , PyTreeFunction ) , \
            "PyTreeFunction: Invalid 'clone' type!" % typename ( clone )
        assert clone or callable  ( the_function ), \
            'PyTreeFunction: Invalid callable: %s'  % typename ( the_function ) 
        ## 
        super(PyTreeFunction,self).__init__ ( tree = tree , clone = clone )
        ## 
        if clone : self.__function = clone.the_function
        else     : self.__function = the_function
        
    @property
    def the_function ( self ) :
        """``the_function'' : the actual function/callable object"""
        return self.__function

    ## clone function 
    def clone ( self , name = "" ) :
        """ Clone it! """
        cloned = PyTreeFunction ( the_function = self.the_function , tree = self.the_tree , clone = self )
        ROOT.SetOwnership ( cloned , False )
        return cloned 

    ## REDUCE 
    def __reduce__  ( self ) : return ptf_factory , ( self.the_function , ) 

    # ==============================================================================
    ## The only one method, delegate to the function    
    def evaluate  ( self ) :
        tree = self.the_tree
        return self.__function ( tree )

    def __str__  ( self ) :
        return '%s(%s/%s)' %  ( 'PyTreeFunction' , typename ( self.the_function ) , prntrf ( self.the_function ) ) 
    __repr__ = __str__

# =============================================================================
## @class PyDataFunction
#  The concrete implementation on of <code>PyFuncData</code>, <code>FuncData</code>
#  that delegates the evaluation to provided function/callable object
#  @see ostap.trees.funcs.FuncData
#  @see Ostap.Functions.PyFuncData
#  @see Ostap.IFuncData 
class PyDataFunction(FuncData) :
    """The concrete implementation on of PyFunData, FuncData,
    that delegates the evaluation to provided function/callable object
    - see ostap.trees.funcs.FuncData
    - see Ostap.Functions.PyFuncData
    - see Ostap.IFuncData 
    """
    def __init__ ( self , the_function , data = None , clone = None  ) :        
        if data is None : data = ROOT.nullptr

        assert not clone or isinstance ( clone , PyDataFunction ) , \
            "PyDataFunction: Invalid 'clone' type!" % typename ( clone )
        assert clone or callable  ( the_function ), \
            'PyDataFunction: Invalid callable %s' % typename ( the_function )
        ## 
        super(PyDataFunction,self).__init__ ( data = data , clone = clone )
        ## 
        if clone : self.__function = clone.the_function
        else     : self.__function = the_function
        
    ## clone function 
    def clone ( self , name = "" ) :
        """ Clone it! """
        cloned = PyDataFunction ( th_function = self.__function , data = self.the_data , clone = self )
        ROOT.SetOwnership ( cloned , False )
        return cloned 

    @property
    def the_function ( self ) :
        """``the_function'' : the actual function/callable object"""
        return self.__function
    
    ## the only one method, delegate to the function  
    def evaluate  ( self ) :
        data = self.the_data 
        return self.__function ( data  )

    def __reduce__  ( self ) : return pdf_factory , ( self.the_function , ) 

    def __str__  ( self ) : return '%s(%s/%s)' %  ( 'PyDataFunction' , typename ( self.the_function ) , prntrf ( self.the_function ) ) 
    __repr__ = __str__

# =============================================================================
## Helper function to (de)serialize PyTreeFunction 
def ptf_factory ( *args ) :
    """ Helper f to (de)serialize PyTreeFunction"""
    return  PyTreeFunction ( *args  )
# =============================================================================
## Helper function to (de)serialize PyFataFunction 
def pdf_factory ( *args ) :
    """ Helper function to (de)serialize PyFataFunction"""
    return  PyDataFunction ( *args  )
    
# =============================================================================
## create the Ostap.ITreeFunc obejct from python function/callable
#  @code
#  def ququ ( tree ) : return tree.pz/tree.pt
#  fun = pyfun_tree ( ququ ) 
#  @endcode
#  @see Ostap::IFuncTree
#  @see Ostap::Functions::PyFuncTree 
#  @see ostap.funcs.FuncTree 
#  @see ostap.funcs.PyTreeFuction
def pyfun_tree ( function , tree = None ) :
    """Create the Ostap.ITreeFunc obejct from python function/callable
    >>> def ququ ( tree ) : return tree.pz/tree.pt
    >>> fun = pyfun_tree ( ququ )
    - see Ostap::IFuncTree
    - see Ostap::Functions::PyFuncTree 
    - see ostap.trees.funcs.FuncTree 
    - see ostap.trees.funcs.PyTreeFuction    
    """
    return PyTreeFunction ( function , tree ) 
    
# =============================================================================
## create the Ostap.IDataFunc obejct from python function/callable
#  @code
#  def ququ ( data ) : return data.pz/data.pt
#  fun = pyfun_data ( ququ ) 
#  @endcode
#  @see Ostap::IFuncData
#  @see Ostap::Functions::PyFuncData
#  @see ostap.funcs.FuncData 
#  @see ostap.funcs.PyDataFuction
def pyfun_data ( function , data = None ) :
    """Create the Ostap.IDataFunc obejct from python function/callable
    >>> def ququ ( data ) : return data.pz/data.pt
    >>> fun = pyfun_data ( ququ )
    - see Ostap::IFuncData
    - see Ostap::Functions::PyFuncData 
    - see ostap.trees.funcs.FuncData
    - see ostap.trees.funcs.PyDataFuction    
    """
    return PyDataFunction ( function , tree ) 
    

# =============================================================================
## print for FuncFormula 
# =============================================================================

def _pff_str_ ( f ) : return f.expression()
def _pfe_str_ ( f ) : return f.expression()
def _pfr_str_ ( f ) : return f.expression()

Ostap.Functions.FuncFormula   .__str__  = _pff_str_ 
Ostap.Functions.FuncFormula   .__repr__ = _pff_str_ 
Ostap.Functions.FuncRooFormula.__str__  = _pfr_str_ 
Ostap.Functions.FuncRooFormula.__repr__ = _pfr_str_ 
Ostap.Functions.Expression    .__str__  = _pfe_str_ 
Ostap.Functions.Expression    .__repr__ = _pfe_str_ 
# =============================================================================

# =============================================================================
## the ITreeFunc based on TTree formula expression 
FuncFormula    = Ostap.Functions.FuncFormula

## the IDataFunc based on RooFormula machinery 
FuncRooFormula = Ostap.Functions.FuncRooFormula

FuncTH1        = Ostap.Functions.FuncTH1
FuncTH2        = Ostap.Functions.FuncTH2
FuncTH3        = Ostap.Functions.FuncTH3


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
