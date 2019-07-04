#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/trees/param.py
#  1,2,3&4D parameterization of datf  from TTree
#  @see Ostap::DataParam 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-07-03
# =============================================================================
""" 1,2,3&4D parameterization of datf  from TTree
- see Ostap.DataParam 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (  
  ) 
# =============================================================================
import ROOT
from   ostap.core.core   import Ostap
import ostap.math.models
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.trees.param' )
else                       : logger = getLogger( __name__ )
# =============================================================================
_large = 2**64 -1 
# =============================================================================
## parameterize 1D unbinned ddistribution from TTree in terms of Legendre sum
#  @code
#  l = LegendreSum ( 5 , -1.0 , 1.0 )
#  tree = ...
#  l.parameterize ( tree , 'X' , 'Y>0' ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-07-03
def _l1_parameterize_ ( l1   ,
                        tree ,
                        var  ,
                        cut = '' , first = 0 , last = _large ) :
    """Parameterize 1D unbinned ddistribution from TTree in terms of Legendre sum
    
    >>> l = LegendreSum ( 5 , -1.0 , 1.0 )
    >>> tree = ...
    >>> l.parameterize ( tree , 'X' , 'Y>0' ) 
    """
    return Ostap.DataParam.parameterize ( tree        , l1    ,
                                          var         ,
                                          str ( cut ) , first , last )

# =============================================================================
## parameterize 2D unbinned ddistribution from TTree in terms of Legendre sum
#  @code
#  l = LegendreSum2 ( 5 , 4 , -1.0 , 1.0 , 0.0 , 1.0 )
#  tree = ...
#  l.parameterize ( tree , 'X' , 'Z' , 'Y>0' ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-07-03
def _l2_parameterize_ ( l2   ,
                        tree ,
                        xvar ,
                        yvar ,
                        cut = '' , first = 0 , last = _large ) :
    """Parameterize 2D unbinned ddistribution from TTree in terms of Legendre sum
    
    >>> l = LegendreSum3 ( 5 , 3 , -1.0 , 1.0 , 0.0 , 1.0  )
    >>> tree = ...
    >>> l.parameterize ( tree , 'X' , 'Y' , 'Z>0' ) 
    """
    return Ostap.DataParam.parameterize ( tree        , l2    ,
                                          xvar        , yvar  ,
                                          str ( cut ) , first , last )

# =============================================================================
## parameterize 3D unbinned ddistribution from TTree in terms of Legendre sum
#  @code
#  l = LegendreSum3 ( 5 , 4 , 2  , -1.0 , 1.0 , 0.0 , 1.0 , -2.0 , 2.0 )
#  tree = ...
#  l.parameterize ( tree , 'X' , 'Y' ,  'Z' , 'T>0' ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-07-03
def _l3_parameterize_ ( l3   ,
                        tree ,
                        xvar ,
                        yvar ,
                        zvar ,
                        cut = '' , first = 0 , last = _large ) :
    """Parameterize 3D unbinned ddistribution from TTree in terms of Legendre sum
    
    >>> l = LegendreSum3 ( 5 , 3 , 2 , -1.0 , 1.0 , 0.0 , 1.0  , 0.0 , 5.0 )
    >>> tree = ...
    >>> l.parameterize ( tree , 'X' , 'Y' , 'Z' , 'T>0' ) 
    """
    return Ostap.DataParam.parameterize ( tree        , l3    ,
                                          xvar        , yvar  , zvar ,
                                          str ( cut ) , first , last )
                                          
# =============================================================================
## parameterize 4D unbinned ddistribuition from TTree in terms of Legendre sum
#  @code
#  l = LegendreSum4 ( 5 , 4 , 2 , 1  , -1.0 , 1.0 , 0.0 , 1.0 , -2.0 , 2.0 , 0. , 10.0  )
#  tree = ...
#  l.parameterize ( tree , 'X' , 'Y' ,  'Z' , 'U' , 'q>0' ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-07-03
def _l4_parameterize_ ( l4   ,
                        tree ,
                        xvar ,
                        yvar ,
                        zvar ,
                        uvar ,
                        cut = '' , first = 0 , last = _large ) :
    """Parameterize 4D unbinned ddistribuition from TTree in terms of Legendre sum
    
    >>> l = LegendreSum4 ( 5 , 3 , 2 , 2 , -1.0 , 1.0 , 0.0 , 1.0  , 0.0 , 5.0 , 0.0 , 1.0 )
    >>> tree = ...
    >>> l.parameterize ( tree , 'X' , 'Y' , 'Z' , 'U' , 'q>0' ) 
    """
    return Ostap.DataParam.parameterize ( tree        , l4    ,
                                          xvar        , yvar  , zvar , uvar ,
                                          str ( cut ) , first , last )

Ostap.Math.LegendreSum .parameterize = _l1_parameterize_
Ostap.Math.LegendreSum2.parameterize = _l2_parameterize_
Ostap.Math.LegendreSum3.parameterize = _l3_parameterize_
Ostap.Math.LegendreSum4.parameterize = _l4_parameterize_

# =============================================================================
## parameterize 1,2,3&4D unbinned distributions in terms of Legeendre series
#  @see Ostap.Math.LegendreSum
#  @see Ostap.Math.LegendreSum2
#  @see Ostap.Math.LegendreSum3
#  @see Ostap.Math.LegendreSum4
#  @code
#  tree  = ....
#  l2    = Ostap.Math.LegendreSum2 ( 5 , -1.0 , 1.0 )
#  tree.parameterize ( l2 , 'X' , 'Y' , cut = 'T>0' , first = 0 , last = 100000 ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-07-03
def _rt_parameterize_ ( tree , func , *args , **kwargs ) :
    """Parameterize 1,2,3&4D unbinned distributions in terms of Legeendre series
    >>> tree  = ....
    >>> l2    = Ostap.Math.LegendreSum2 ( 5 , -1.0 , 1.0 )
    >>> tree.parameterize ( l2 , 'X' , 'Y' , cut = 'T>0' , first = 0 , last = 100000 )
    
    - see Ostap.Math.LegendreSum
    - see Ostap.Math.LegendreSum2
    - see Ostap.Math.LegendreSum3
    - see Ostap.Math.LegendreSum4
    """
    assert isinstance ( func , ( Ostap.Math.LegendreSum  ,
                                 Ostap.Math.LegendreSum2 ,
                                 Ostap.Math.LegendreSum3 ,
                                 Ostap.Math.LegendreSum4 ) ) ,\
                                 'Invalid function object %s/%s' % ( func , type  ( func) )
    
    return func.parameterize ( tree , *args , **kwargs )

    
ROOT.TTree.parameterize = _rt_parameterize_


_decorated_classes_ = (
    ROOT.TTree             ,
    Ostap.Math.LegendreSum  ,
    Ostap.Math.LegendreSum2 ,
    Ostap.Math.LegendreSum3 , 
    Ostap.Math.LegendreSum4 , 
    )

_new_methods_       = (
    Ostap.Math.LegendreSum .parameterize , 
    Ostap.Math.LegendreSum2.parameterize ,
    Ostap.Math.LegendreSum3.parameterize , 
    Ostap.Math.LegendreSum4.parameterize , 
    ROOT.TTree.parameterize ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
