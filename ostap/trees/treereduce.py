#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/treereduce.py
#  Reduction/serialization/deserialization for some Ostap/ROOT classes  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
""" Reduction/serialization/deserialization for some Ostap/ROOT classes
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2022-09-09"
__all__     = () 
# =============================================================================
from   ostap.core.core     import Ostap
from   ostap.math.reduce   import root_factory
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.treereduce' )
else                       : logger = getLogger( __name__ )
# =============================================================================


# =============================================================================
## reduce Ostap::Functions::FunTH1 
def _ofth1_reduce_ ( fun ):
    """Reduce Ostap.Functions.FuncTH1"""
    return root_factory , ( type ( fun ) ,
                            fun.histo () ,
                            fun.x     () )
# =============================================================================
## reduce Ostap::Functions::FunTH3 
def _ofth2_reduce_ ( fun ):
    """Reduce Ostap.Functions.FuncTH2"""
    return root_factory , ( type ( fun ) ,
                            fun.histo () ,
                            fun.x     () , 
                            fun.y     () )
# =============================================================================
## reduce Ostap::Functions::FunTH3
def _ofth3_reduce_ ( fun ):
    """Reduce Ostap.Functions.FuncTH3"""
    return root_factory , ( type ( fun ) ,
                            fun.histo () ,
                            fun.x     () ,
                            fun.y     () ,
                            fun.z     () )
# =============================================================================
## reduce Ostap::Functions::FuncFormula
def _offf_reduce_ ( fun ):
    """Reduce Ostap.Functions.FuncFormula"""
    return root_factory , ( type ( fun ) , fun.expression() )

Ostap.Functions.FuncTH1     . __reduce__ = _ofth1_reduce_
Ostap.Functions.FuncTH2     . __reduce__ = _ofth2_reduce_
Ostap.Functions.FuncTH3     . __reduce__ = _ofth3_reduce_
Ostap.Functions.FuncFormula . __reduce__ =  _offf_reduce_

_decorated_classes = (
    Ostap.Functions.FuncTH1     ,
    Ostap.Functions.FuncTH2     ,
    Ostap.Functions.FuncTH3     ,
    Ostap.Functions.FuncFormula ,
    )

_new_methods_ = ( 
    Ostap.Functions.FuncTH1     . __reduce__ , 
    Ostap.Functions.FuncTH2     . __reduce__ , 
    Ostap.Functions.FuncTH3     . __reduce__ , 
    Ostap.Functions.FuncFormula . __reduce__ ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
