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
from   ostap.utils.basic   import typename 
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
    """ Reduce Ostap.Functions.FuncTH1"""
    return root_factory , ( type ( fun ) ,
                            fun.histo () ,
                            fun.x     () )
# =============================================================================
## reduce Ostap::Functions::FunTH3 
def _ofth2_reduce_ ( fun ):
    """ Reduce Ostap.Functions.FuncTH2"""
    return root_factory , ( type ( fun ) ,
                            fun.histo () ,
                            fun.x     () , 
                            fun.y     () )
# =============================================================================
## reduce Ostap::Functions::FunTH3
def _ofth3_reduce_ ( fun ):
    """ Reduce Ostap.Functions.FuncTH3"""
    return root_factory , ( type ( fun ) ,
                            fun.histo () ,
                            fun.x     () ,
                            fun.y     () ,
                            fun.z     () )
# =============================================================================
## reduce Ostap::Functions::FuncFormula
def _offf_reduce_ ( fun ):
    """Reduce Ostap.Functions.FuncFormula & Ostap.Functions.Expression
    """
    return root_factory , ( type ( fun ) , fun.expression() )

Ostap.Functions.FuncTH1        . __reduce__ = _ofth1_reduce_
Ostap.Functions.FuncTH2        . __reduce__ = _ofth2_reduce_
Ostap.Functions.FuncTH3        . __reduce__ = _ofth3_reduce_
Ostap.Functions.FuncFormula    . __reduce__ =  _offf_reduce_
Ostap.Functions.Expression     . __reduce__ =  _offf_reduce_

for o in ( Ostap.IFuncTree                ,
           Ostap.IFuncData                ,
           #
           Ostap.Functions.PyFuncTree     ,
           Ostap.Functions.PyFuncData     ,
           #
           Ostap.Functions.Func1D         ,
           Ostap.Functions.Func2D         ,
           Ostap.Functions.Func3D         ,
           #
           Ostap.Functions.FuncTH1        ,
           Ostap.Functions.FuncTH2        ,
           Ostap.Functions.FuncTH3        ,
           #
           Ostap.Functions.FuncFormula    , 
           Ostap.Functions.Expression     ,  
           #           
           Ostap.Functions.RooTreeFun     ,  
           # 
           Ostap.Functions.FuncRooFormula ,
           #
           Ostap.Functions.FuncRoo1D      ,
           Ostap.Functions.FuncRoo1D      ,
           Ostap.Functions.FuncRoo1D      )  : 
           # 
    o.__str__  = lambda s : typename ( s )
    o.__repr__ = lambda s : typename ( s )
    
def _q_str_  ( s ) :  return "%s(%s)"       % ( typename ( s ) , s.expression() )
def _q_str3_ ( s ) :  return "%s(%s,%s,%s)" % ( typename ( s ) , s.x() , s.y() , s.z() )
def _q_str2_ ( s ) :  return "%s(%s,%s)"    % ( typename ( s ) , s.x() , s.y() )
def _q_str1_ ( s ) :  return "%s(%s)"       % ( typename ( s ) , s.x() )

for o in ( Ostap.Functions.FuncFormula    ,
           Ostap.Functions.FuncRooFormula ,
           Ostap.Functions.Expression     ) :
    
    o.__str__  = _q_str_
    o.__repr__ = _q_str_


for o in ( Ostap.Functions.Func1D     ,
           Ostap.Functions.FuncTH1    ,
           Ostap.Functions.FuncRoo1D  ,
           Ostap.Functions.FuncRooTH1 ) :
    
    o.__str__  = _q_str1_
    o.__repr__ = _q_str1_

for o in ( Ostap.Functions.Func2D     ,
           Ostap.Functions.FuncTH2    ,
           Ostap.Functions.FuncRoo2D  ,
           Ostap.Functions.FuncRooTH2 ) :
    
    o.__str__  = _q_str2_
    o.__repr__ = _q_str2_

for o in ( Ostap.Functions.Func3D     ,
           Ostap.Functions.FuncTH3    ,
           Ostap.Functions.FuncRoo3D  ,
           Ostap.Functions.FuncRooTH3 ) :
    
    o.__str__  = _q_str3_
    o.__repr__ = _q_str3_

_decorated_classes = ( Ostap.IFuncTree                ,
                       Ostap.IFuncData                ,
                       #
                       Ostap.Functions.PyFuncTree     ,
                       Ostap.Functions.PyFuncData     ,
                       #
                       Ostap.Functions.Func1D         ,
                       Ostap.Functions.Func2D         ,
                       Ostap.Functions.Func3D         ,
                       #
                       Ostap.Functions.FuncTH1        ,
                       Ostap.Functions.FuncTH2        ,
                       Ostap.Functions.FuncTH3        ,
                       #
                       Ostap.Functions.FuncFormula    , 
                       Ostap.Functions.Expression     ,  
                       #           
                       Ostap.Functions.RooTreeFun     ,  
                       # 
                       Ostap.Functions.FuncRooFormula ,
                       #
                       Ostap.Functions.FuncRoo1D      ,
                       Ostap.Functions.FuncRoo1D      ,
                       Ostap.Functions.FuncRoo1D      )

_new_methods_ = ( Ostap.Functions.FuncTH1        . __reduce__ , 
                  Ostap.Functions.FuncTH2        . __reduce__ , 
                  Ostap.Functions.FuncTH3        . __reduce__ , 
                  Ostap.Functions.FuncFormula    . __reduce__ ,
                  Ostap.Functions.Expression     . __reduce__ )    + \
                  tuple ( t.__str__  for t in _decorated_classes ) + \
                  tuple ( t.__repr__ for t in _decorated_classes ) 


import pickle 
def _no_reduce_ ( t ) :
    logger.always ( "I AM REDUCE %s" % type ( t )  )
    raise pickle.PicklingError ( 'Class %%s is nonpickleable!' % typename ( t ) )

notpickleable_types = ( Ostap.Functions.Func1D , 
                        Ostap.Functions.Func2D ,
                        Ostap.Functions.Func3D )

## for k in  notpickleable_types :
##    print ( ' had reduce?' , k , typename  ( k ) , hasattr ( k , '__reduce__' ) ) 
##    k.__reduce__ = _no_reduce_

## import ostap.io.pickling as OP
## checker = OP.PickleChecker()
## checker.add_nonpickleable ( *nonpickleable_types )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
