#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/covtransform.py
#  Transformation of covariance matrices
#  for \f$ y = y ( x ) \f$, it gets
#  \f$  C_y = J C_x J^\mathrm{T} \f$,
#  where \f$ J = \left( \frac{\partial y }{\partial x } \right) \f$
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2020-05-14
# =============================================================================
""" Transformation of covariand matrices 
- for y = y ( x ) it gets  C(y) = J C(x) J^T,
- where J is Jacobi matrix 
"""
# =============================================================================
from   __future__  import print_function
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = ""
# =============================================================================
__all__     = (
    'transform'  , ## transfrom covarinance matrix 
    )
# =============================================================================
from   builtins  import range 
import ROOT
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger   import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.covtransform' )
else                       : logger = getLogger ( __name__                  )
# =============================================================================
from   ostap.core.core       import Ostap, VE 
from   ostap.math.derivative import Partial
import ostap.math.linalg
# =============================================================================


# =============================================================================
## Transform the covariance nmatrix C at point X to the variables Y(X)
#  for \f$ y = y ( x ) \f$, it gets
#  \f[  C(y)= J C(x) J^\mathrm{T} \f],
#  where \f$ J = \left( \frac{\partial y }{\partial x } \right) \f$
#  @code
#  X = 1 , 2
#  C = Ostap.SymMatrix(2)()
#  C[ 0 , 0 ] = 0.20
#  C[ 1 , 1 ] = 0.05
#  C[ 1 , 1 ] = 0.30
#  r   = lambda x , y : (x*x+y*y)**2
#  phi = lambda x , y : math.atan2 ( y , x )
#  C_polar = transform ( C , X , r , phi ) 
#  @endcode
#  @param C  "old" covatiance matrix
#  @param X  "old" varibales  (arary iof values)
#  @param Y  "new" variables  (array of callables)
#  @return covarinance matrix for variables Y  
def transform ( C , X , *Y ) :
    """ Transform the covariance nmatrix C at point X to the variables Y(X)
    >>> X = 1 , 2
    >>> C = Ostap.SymMatrix(2)()
    >>> C [ 0 , 0 ] = 0.20
    >>> C [ 0 , 1 ] = 0.05
    >>> C [ 1 , 1 ] = 0.30
    >>> r   = lambda x , y : (x*x+y*y)**2
    >>> phi = lambda x , y : math.atan2 ( y , x )
    >>> C_polar = transform ( C , X , r , phi ) 
    """

    ny = len ( Y ) 
    assert 1 <= ny , 'Invalid size of Y!'

    nx = len ( X ) 
    if C is None and 1 <= nx : 
        C = Ostap.SymMatrix ( nx ) ()
        for i , x  in enumerate ( X ) :
            xx = VE ( x ) 
            C [ i, i ] = xx.cov2 ()

    shape = C.shape
    assert shape [ 0 ] == shape[1] and shape[0] == nx  , 'Invalid shape of matrix C!' 


    CC  = Ostap.SymMatrix  ( nx ) ()
    for i in range  ( CC.kRows ) :
        CC [ i , i ] = C  ( i , i ) 
        for j in range  ( i + 1  , CC.kCols ) :
            v = 0.5 * ( C ( i , j ) + C ( j , i ) )
            CC [ i , j ]  = v
            
    XX = Ostap.Vector ( nx  ) ()
    for i , x  in enumerate ( X ) : 
        XX [ i ] = float ( x )

    ## get vector-with-errors 
    XX = Ostap.VectorE ( nx ) ( XX , CC )

    R  = XX.transform ( *Y )
    
    return R.cov2()  

# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
