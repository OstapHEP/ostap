#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with helper base classes ...
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Decoration of Tree/Chain objects for efficient use in python"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    ## 'FuncTree'       , ## helper base class for 'TTree-function'
    ## 'FuncData'       , ## helper base class for 'RooAbsData-function'
    'FormulaFunc'    , ## simple wrapper over TTreeFormula/Ostap::Formula  
    'RooFormulaFunc' , ## simple wrapper over RooFormulaVar  
    ) 
# =============================================================================
import ROOT
from   ostap.core.core import Ostap, valid_pointer
# =============================================================================
## ## @class FuncTree
## #  Helper class to implement "TTree-function"
## #  @see Ostap::Functions::PyFuncTree
## class FuncTree(Ostap.Functions.PyFuncTree) :
##     """Helper class to implement ``TTree-function''
##     """
##     def __init__ ( self , tree = None ) :
##         ## initialize the base class 
##         Ostap.Functions.PyFuncTree.__init__ ( self , data  ) 
        
##     ## the main method 
##     def evaluate ( self ) :
##         tree =  self.tree
##         assert valid_pointer ( tree ) , 'Invalid TTree object'
##         return -1

## # =============================================================================
## ## @class FuncData
## #  Helper class to implement "RooAbsData-function"
## #  @see Ostap::Functions::PyFuncData
## class FuncData(Ostap.Functions.PyFuncData) :
##     """Helper class to implement ``TTree-function''
##     """
##     def __init__ ( self , data = None ) :
##         ## initialize the base class 
##         Ostap.Functions.PyFuncData.__init__ ( self , data  ) 
        
##     ## the main method 
##     def evaluate ( self ) :
##         data =  self.data
##         assert valid_pointer ( data ), 'Invalid RooAbsData object'
##         return -1
    
# ==============================================================================
## @class FormulaFunc
#  simple  class that uses Ostap::Formula/TTreeFormula to get ``TTree-function''
#  @see Ostap::Formula
#  @see TTreeFormula
#  @see Ostap::Functions::FuncFormula 
class FormulaFunc(Ostap.Functions.FuncFormula) :
    def __init__ (  self,  expression , tree = None  , name = '' ) :
        Ostap.Functions.FuncFormula.__init__ ( self , expression , tree , name ) 

# ==============================================================================
## @class RooFormulaFunc
#  simple  class that uses RooFormulaVar to get ``RooAbsData-function''
#  @see RooFormulaVar
#  @see Ostap::Functions::FuncRooFormula 
class RooFormulaFunc(Ostap.Functions.FuncFormula) :
    def __init__ (  self,  expression , tree = None  , name = '' ) :
        Ostap.Functions.FuncRooFormula.__init__ ( self , expression , tree , name ) 
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
