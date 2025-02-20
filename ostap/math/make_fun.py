#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/make_fun.py
#  Trivial utilities to convert python function/callable into C++ function
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2023-02-14
# =============================================================================
""" Trivial utilities to convert python function/callable into C++ function
#"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2023-02-14"
__all__     = (
    'make_fun1'  , ## convert 1D python function/callable into C++ function
    'make_fun2'  , ## convert 2D python function/callable into C++ function
    'make_fun3'  , ## convert 3D python function/callable into C++ function    
    )
# =============================================================================
from ostap.core.core      import Ostap
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.make_fun' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## convert python 1D-function into C++ functor 
def make_fun1 ( fun , forcepc = False ) :
    """ Convert python 1D-function into C++ functor"""
    assert callable ( fun ) , 'Function must be callable!'
    ##
    if forcepc :
        PC  = Ostap.Functions.PyCallable
        return fun if isinstance ( fun , PC ) else PC ( fun , True )
    ##
    return Ostap.Math.Apply  ( fun )
## convert python 2D-function into C++ functor 
def make_fun2 ( fun , forcepc = False ) :
    """ Convert python 2D-function into C++ functor"""
    assert callable ( fun ) , 'Function must be callable!'
    ##
    if forcepc :
        PC2 = Ostap.Functions.PyCallable2
        return fun if isinstance ( fun , PC2 ) else PC2 ( fun , True )
    ##
    return Ostap.Math.Apply2 ( fun )
## convert python 3D-function into C++ functor 
def make_fun3 ( fun , forcepc = False ) :
    """ Convert python 3D-function into C++ functor"""
    assert callable ( fun ) , 'Function must be callable!'
    ##
    if forcepc :
        PC3 = Ostap.Functions.PyCallable3
        return fun if isinstance ( fun , PC3 ) else PC3 ( fun , True )
    ##
    return Ostap.Math.Apply3 ( fun )
    # =========================================================================

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
