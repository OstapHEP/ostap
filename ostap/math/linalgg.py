#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/linalgg.py
#  Few utilities to simplify linear algebra manipulations using GSL 
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
"""Few utilities to simplify linear algebra manipulations unisg GSL 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = ""
# =============================================================================
__all__     = ( )
# =============================================================================
from   ostap.math.base   import Ostap
import ostap.math.linalg as     LA 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.linalgg' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
Matrix      = Ostap.GSL.Matrix
Vector      = Ostap.GSL.Vector
Permutation = Ostap.GSL.Permutation

# =============================================================================
def _m_iadd_ ( m , value ) :
    """ Add someting to matrtx 
    """
    if isinstance ( value , Matrix ) :
        if m.nRows() != value.nRows () : return NotImplemented
        if m.nCols() != value.nCols () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                               : return NotImplemented 
                    
    m.iadd ( value )
    
    return m

# =============================================================================
## Subtract from the matrix
def _m_isub_ ( m , value ) :
    """ Subtract from the matrix
    """
    if isinstance ( value , Matrix ) :
        if m.nRows() != value.nRows () : return NotImplemented
        if m.nCols() != value.nCols () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                               : return NotImplemented 
                    
    m.isub ( value )
    
    return m 


# =============================================================================
## Multiply by the matrox 
def _m_imul_ ( m , value ) :
    """ Multiply by the matrix 
    """
    if isinstance ( value , Matrix ) :
        if m.nCols() != value.nRows () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                               : return NotImplemented 
                    
    m.imul ( value )
    
    return m 

# =============================================================================
## scale/divdie by the matrix 
def _m_imul_ ( m , value ) :
    """ Multiply by the matrix 
    """
    if isinstance ( value , num_types ) : value = float ( value ) 
    else                                : return NotImplemented 
                    
    m.idiv ( value )
    
    return m 

# ==============================================================================
## multiply matrix 
def _m_mul_ ( m ,  value ) :
    
    if   isinstance ( value , Matrix ) :
        if m.nCols() != value.nRows () : return NotImplemented
    elif isinstance ( value , Vector ) :
        if m.nCols() != value.size  () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                               : return NotImplemented 
    
    return m.multiply ( value )
    

# ==============================================================================
## Convert GSL matrix to SMatrix 
def _to_smatrix_ ( mtrx ) :
    """ Convert GSL matrix to SMatrix 
    """
    ## get dimension 
    nr , nr = mtrx.nRows() , mtrx.nCols()
    ## get the result 
    result  = Ostap.Math.Matrix ( nr , nc )() 
    ## fill it!
    for i in range ( nr ) :
        for j in range ( nc ) :
            result [ i , j ] = mtrx ( i , j )
    return result

Matrix.asSMatrix  = _to_smatrix_
Matrix.to_SMatrix = _to_smatrix_
Matrix.as_SMatrix = _to_smatrix_
Matrix.to_smatrix = _to_smatrix_
Matrix.as_smatrix = _to_smatrix_

# ==============================================================================
## Convert GSL matrix to symmetric S-matrix 
def _to_ssymmatrix_ ( mtrx ) :
    """ Convert GSL matrix to (symmetric) SMatrix 
    """
    ## get dimension 
    nr , nr = mtrx.nRows() , mtrx.nCols()
    assert nr and nr == nc , "Impossible to crrate symmetruc from rectangular matrix!"
    
    ## get the result 
    result  = Ostap.Math.SymMatrix ( nr  )() 
    ## fill it!
    for i in range ( nr ) :
        result [ i , i ] = mtrx ( i , i ) 
        for j in range ( i + 1  , nr ) :
            result [ i , j ] =  0.5 * ( mtrx ( i , j ) + mtrx ( j , i ) )
            
    return result

Matrix.asSymSMatrix  = _to_ssymmatrix_
Matrix.to_SymSMatrix = _to_ssymmatrix_
Matrix.as_SymSMatrix = _to_ssymmatrix_
Matrix.to_symsmatrix = _to_ssymmatrix_
Matrix.as_sysmmatrix = _to_ssymmatrix_

# ==============================================================================
## Convert GSL matrix to TMatrix 
def _to_tmatrix_ ( mtrx ) :
    """ Convert GSL matrix to TMatrix 
    """
    ## get dimension 
    nr , nr = mtrx.nRows() , mtrx.nCols()
    ## get the result 
    result  = Ostap.Math.TMatrixD ( nr , nc ) 
    ## fill it!
    for i in range ( nr ) :
        for j in range ( nc ) :
            result [ i , j ] = mtrx ( i , j )
    ## 
    return result

Matrix.asTMatrix  = _to_tmatrix_
Matrix.to_TMatrix = _to_tmatrix_
Matrix.as_TMatrix = _to_tmatrix_
Matrix.to_tmatrix = _to_tmatrix_
Matrix.as_tmatrix = _to_tmatrix_

# ==============================================================================
## Convert GSL matrix to symmetric 
def _to_symtmatrix_ ( mtrx ) :
    """ Convert GSL matrix to (symmetric) TMatrix 
    """
    ## get dimension 
    nr , nr = mtrx.nRows() , mtrx.nCols()
    assert nr and nr == nc , "Impossible to crrate symmetruc from rectangular matrix!"
    
    ## get the result 
    result  = Ostap.Math.TMatrixSymD ( nr  )
    ## fill it!
    for i in range ( nr ) :
        result [ i , i ] = mtrx ( i , i ) 
        for j in range ( i + 1  , nr ) :
            rij = 0.5 * ( mtrx ( i , j ) + mtrx ( j , i ) ) 
            result [ i , j ] = rij
            result [ j , i ] = rij
            
    return result

Matrix.asSymTMatrix  = _to_symtmatrix_
Matrix.to_SymTMatrix = _to_symtmatrix_
Matrix.as_SymTMatrix = _to_symtmatrix_
Matrix.to_symtmatrix = _to_symtmatrix_
Matrix.as_symtmatrix = _to_symtmatrix_

Matrix.GetNrows      = Matrix.nRows
Matrix.GetNcols      = Matrix.nCols

Matrix.kRows = property ( Matrix.nRows , None , None , "`kRows` : number of rows "    )
Matrix.kCols = property ( Matrix.nCols , None , None , "`kCols` : number of columns " )


def _m_pretty_print_ ( mtrx , **kwargs ) :    
    return LA.LinAlgT.M_PRETTY ( mtrx , **kwargs )
def _m_str_          ( mtrx , **kwargs ) :    
    return LA.LinAlgT.M_STR    ( mtrx , **kwargs )

Matrix.pretty_print  = _m_pretty_print_
Matrix.table         = _m_str_ 
Matrix.__str__       = _m_str_ 
Matrix.__repr__      = _m_str_ 

def _v_pretty_print_ ( mtrx , **kwargs ) :    
    return LA.LinAlgT.V_PRETTY ( mtrx , **kwargs )
def _v_str_          ( mtrx , **kwargs ) :    
    return LA.LinAlgT.V_STR    ( mtrx , **kwargs )

Vector.pretty_print  = _v_pretty_print_
Vector.table         = _v_str_ 
Vector.__str__       = _v_str_ 
Vector.__repr__      = _v_str_ 

Permutation.pretty_print  = _v_pretty_print_
Permutation.table         = _v_str_ 
Permutation.__str__       = _v_str_ 
Permutation.__repr__      = _v_str_ 

_new_methods_ = (
    ##
    Matrix.asSMatrix         ,
    Matrix.to_SMatrix        ,
    Matrix.as_SMatrix        ,
    Matrix.to_smatrix        ,
    Matrix.as_smatrix        ,
    ## 
    Matrix.asSymSMatrix      , 
    Matrix.to_SymSMatrix     ,
    Matrix.as_SymSMatrix     ,
    Matrix.to_symsmatrix     ,
    Matrix.as_sysmmatrix     ,
    ##
    Matrix.asTMatrix         ,
    Matrix.to_TMatrix        ,
    Matrix.as_TMatrix        ,
    Matrix.to_tmatrix        ,
    Matrix.as_tmatrix        ,
    ## 
    Matrix.asSymTMatrix      ,
    Matrix.to_SymTMatrix     ,
    Matrix.as_SymTMatrix     ,
    Matrix.to_symtmatrix     ,
    Matrix.as_symtmatrix     ,
    ##
    Matrix.kRows              , 
    Matrix.kCols              ,
    Matrix.GetNrows           , 
    Matrix.GetNcols           , 
    ##
    Matrix.pretty_print       , 
    Matrix.table              , 
    Matrix.__str__            , 
    Matrix.__repr__           , 
    ##
    Vector.pretty_print       , 
    Vector.table              , 
    Vector.__str__            , 
    Vector.__repr__           , 
    ##
    Permutation.pretty_print  , 
    Permutation.table         , 
    Permutation.__str__       , 
    Permutation.__repr__      , 
    ## 
)
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    

# =============================================================================
##                                                                      The END
# =============================================================================
