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
    'param_types_1D' , ## list of valid 1D-types for parameterisation
    'param_types_2D' , ## list of valid 1D-types for parameterisation
    'param_types_3D' , ## list of valid 1D-types for parameterisation
    'param_types_4D' , ## list of valid 1D-types for parameterisation
  ) 
# =============================================================================
from   ostap.core.core   import Ostap
import ostap.math.models
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.trees.param' )
else                       : logger = getLogger( __name__ )
# =============================================================================
_large = 2**64 -1 
# =============================================================================
## valid non-histogram projection types 
# ============================================================================-
param_types_1D = Ostap.Math.LegendreSum  , Ostap.Math.Bernstein   , Ostap.Math.ChebyshevSum
param_types_2D = Ostap.Math.LegendreSum2 , Ostap.Math.Bernstein2D , 
param_types_3D = Ostap.Math.LegendreSum3 , Ostap.Math.Bernstein3D , 
param_types_4D = Ostap.Math.LegendreSum4 , 
# =============================================================================
## parameterize 1D unbinned distribution from TTree in terms of
#  Legendre/chebyshev/Bernstein sums
#  @code
#  l = LegendreSum ( 5 , -1.0 , 1.0 )
#  tree = ...
#  l.parameterize ( tree , 'X' , 'Y>0' )
#
#  c = ChebyshevSum ( 5 , -1.0 , 1.0 )
#  tree = ...
#  l.parameterize ( tree , 'X' , 'Y>0' )
#
#  b = Bernstein  ( 5 , -1.0 , 1.0 )
#  tree = ...
#  l.parameterize ( tree , 'X' , 'Y>0' ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-07-03
def _l1_parameterize_ ( l1   ,
                        tree ,
                        var  ,
                        cut = '' , first = 0 , last = _large ) :
    """Parameterize 1D unbinned distribution from TTree in terms of Legendre or Chebyshev sum
    
    >>> l = LegendreSum ( 5 , -1.0 , 1.0 )
    >>> tree = ...
    >>> l.parameterize ( tree , 'X' , 'Y>0' )
    
    >>> c = ChebyshevSum ( 5 , -1.0 , 1.0 )
    >>> tree = ...
    >>> c.parameterize ( tree , 'X' , 'Y>0' )
    
    >>> b = Bernstein  ( 5 , -1.0 , 1.0 )
    >>> tree = ...
    >>> b.parameterize ( tree , 'X' , 'Y>0' ) 
    """
    return Ostap.DataParam.parameterize ( tree        , l1    ,
                                          var         ,
                                          str ( cut ) , first , last )

# =============================================================================
## parameterize 2D unbinned distribution from TTree in terms of Legendre sum
#  @code
#  tree = ...
#  l = LegendreSum2 ( 5 , 4 , -1.0 , 1.0 , 0.0 , 1.0 )
#  l.parameterize ( tree , 'X' , 'Z' , 'Y>0' ) 
#
#  b = Bernsteinn2D ( 5 , 4 , -1.0 , 1.0 , 0.0 , 1.0 )
#  b.parameterize ( tree , 'X' , 'Z' , 'Y>0' ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-07-03
def _l2_parameterize_ ( l2   ,
                        tree ,
                        xvar ,
                        yvar ,
                        cut = '' , first = 0 , last = _large ) :
    """Parameterize 2D unbinned distribution from TTree in terms of Legendre sum
    
    >>> tree = ...
    >>> l = LegendreSum2 ( 5 , 3 , -1.0 , 1.0 , 0.0 , 1.0  )
    >>> l.parameterize (  tree , 'X' , 'Y' , 'Z>0' ) 
    >>> b = BErnstein2D  ( 5 , 3 , -1.0 , 1.0 , 0.0 , 1.0  )
    >>> b.parameterize ( tree , 'X' , 'Y' , 'Z>0' ) 
    """
    return Ostap.DataParam.parameterize ( tree        , l2    ,
                                          xvar        , yvar  ,
                                          str ( cut ) , first , last )

# =============================================================================
## parameterize 3D unbinned distribution from TTree in terms of Legendre sum
#  @code
#  tree = ...
#  l = LegendreSum3 ( 5 , 4 , 2  , -1.0 , 1.0 , 0.0 , 1.0 , -2.0 , 2.0 )
#  l.parameterize ( tree , 'X' , 'Y' ,  'Z' , 'T>0' ) 
#  b = Bernsteinn3D ( 5 , 4 , 2  , -1.0 , 1.0 , 0.0 , 1.0 , -2.0 , 2.0 )
#  b.parameterize ( tree , 'X' , 'Y' ,  'Z' , 'T>0' ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-07-03
def _l3_parameterize_ ( l3   ,
                        tree ,
                        xvar ,
                        yvar ,
                        zvar ,
                        cut = '' , first = 0 , last = _large ) :
    """Parameterize 3D unbinned distribution from TTree in terms of Legendre sum
    
    >>> tree = ...
    >>> l = LegendreSum3 ( 5 , 3 , 2 , -1.0 , 1.0 , 0.0 , 1.0  , 0.0 , 5.0 )
    >>> l.parameterize ( tree , 'X' , 'Y' , 'Z' , 'T>0' ) 
    >>> b = Bernsteinn3D ( 5 , 3 , 2 , -1.0 , 1.0 , 0.0 , 1.0  , 0.0 , 5.0 )
    >>> b.parameterize ( tree , 'X' , 'Y' , 'Z' , 'T>0' ) 
    """
    return Ostap.DataParam.parameterize ( tree        , l3    ,
                                          xvar        , yvar  , zvar ,
                                          str ( cut ) , first , last )
                                          
# =============================================================================
## parameterize 4D unbinned distribution from TTree in terms of Legendre sum
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
    """Parameterize 4D unbinned distribuition from TTree in terms of Legendre sum
    
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

Ostap.Math.ChebyshevSum.parameterize = _l1_parameterize_
Ostap.Math.Bernstein   .parameterize = _l1_parameterize_
Ostap.Math.Bernstein2D .parameterize = _l2_parameterize_
Ostap.Math.Bernstein3D .parameterize = _l3_parameterize_

# =============================================================================
## parameterize 1,2,3&4D unbinned distributions in terms of Legendre series
#  @see Ostap.Math.LegendreSum
#  @see Ostap.Math.LegendreSum2
#  @see Ostap.Math.LegendreSum3
#  @see Ostap.Math.LegendreSum4
#  @see Ostap.Math.ChebyshevSum
#  @see Ostap.Math.Bernstein
#  @see Ostap.Math.Bernstein2D
#  @see Ostap.Math.Bernstein3D
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
    - see Ostap.Math.ChebyshevSum
    - see Ostap.Math.Bernstein
    - see Ostap.Math.Bernstein2D
    - see Ostap.Math.Bernstein3D
    """
    assert isinstance ( func , ( Ostap.Math.LegendreSum  ,
                                 Ostap.Math.LegendreSum2 ,
                                 Ostap.Math.LegendreSum3 ,
                                 Ostap.Math.LegendreSum4 ,
                                 Ostap.Math.ChebyshevSum ,
                                 Ostap.Math.Bernstein    ,
                                 Ostap.Math.Bernstein2D  ,
                                 Ostap.Math.Bernstein3D  ) ) ,\
                                 'Invalid function object %s/%s' % ( func , type  ( func) )
    
    return func.parameterize ( tree , *args , **kwargs )

    
ROOT.TTree.parameterize = _rt_parameterize_

_decorated_classes_ = (
    ROOT.TTree              ,
    Ostap.Math.LegendreSum  ,
    Ostap.Math.LegendreSum2 ,
    Ostap.Math.LegendreSum3 , 
    Ostap.Math.LegendreSum4 , 
    Ostap.Math.ChebyshevSum ,
    Ostap.Math.Bernstein    ,
    Ostap.Math.Bernstein2D  ,
    Ostap.Math.Bernstein3D  ,
    )

_new_methods_       = (
    Ostap.Math.LegendreSum .parameterize , 
    Ostap.Math.LegendreSum2.parameterize ,
    Ostap.Math.LegendreSum3.parameterize , 
    Ostap.Math.LegendreSum4.parameterize , 
    Ostap.Math.ChebyshevSum.parameterize , 
    Ostap.Math.Bernstein   .parameterize , 
    Ostap.Math.Bernstein2D .parameterize , 
    Ostap.Math.Bernstein3D .parameterize , 
    ROOT.TTree.parameterize ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
