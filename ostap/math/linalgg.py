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
from   ostap.core.ostap_types import num_types 
from   ostap.math.base        import Ostap
import ostap.math.linalg      as     LA 
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
## matrix += value 
def _m_iadd_ ( m , value ) :
    """ matrix += value
    """
    if isinstance ( value , Matrix ) :
        if m.nRows() != value.nRows () : return NotImplemented
        if m.nCols() != value.nCols () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                               : return NotImplemented 
    ## 
    m.iadd ( value )    
    return m

# =============================================================================
## matrix -= value 
def _m_isub_ ( m , value ) :
    """ matrix -= value 
    """
    if isinstance ( value , Matrix ) :
        if m.nRows() != value.nRows () : return NotImplemented
        if m.nCols() != value.nCols () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                               : return NotImplemented 
    ## 
    m.isub ( value )    
    return m 


# =============================================================================
## matrix *= value 
def _m_imul_ ( m , value ) :
    """ matrix *= vale 
    """
    if isinstance ( value , Matrix ) :
        if m.nCols() != value.nRows () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                               : return NotImplemented 
    ## 
    m.imul ( value )    
    return m 

# =============================================================================
## matrix /= value 
def _m_idiv_ ( m , value ) :
    """ matrix /= value 
    """
    if isinstance ( value , num_types ) : value = float ( value ) 
    else                                : return NotImplemented 
    ## 
    m.idiv ( value )
    return m 

# ============================================================================
## matrix + value 
def _m_add_ ( m , value ) :
    """ matrix + value 
    """
    if   isinstance ( value , Matrix ) :
        if m.nRows() != value.nRows () : return NotImplemented
        if m.nCols() != value.nCols () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                               : return NotImplemented 
    # 
    result  = Matrix ( m )
    result += value
    #
    return result 

# ============================================================================
## matrix - value 
def _m_sub_ ( m , value ) :
    """ Matrix - value 
    """
    if   isinstance ( value , Matrix ) :
        if m.nRows() != value.nRows () : return NotImplemented
        if m.nCols() != value.nCols () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                               : return NotImplemented 
    # 
    result  = Matrix ( m )
    result -= value
    #
    return result 

# ==============================================================================
## matrix * value 
def _m_mul_ ( m ,  value ) :
    """ matrix * value 
    """
    if   isinstance ( value , Matrix ) :
        if m.nCols() != value.nRows () : return NotImplemented
    elif isinstance ( value , Vector ) :
        if m.nCols() != value.size  () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                               : return NotImplemented 
    ## 
    result = m.multiply ( value )
    return result

# ==============================================================================
## matrix / value 
def _m_div_ ( m ,  value ) :
    """ matrix / value 
    """
    if isinstance ( value , num_types ) : value = float ( value ) 
    else                               : return NotImplemented 
    ## 
    result  = Matrix ( m )
    result /= value 
    return result 

# =============================================================================
## value + matrix 
def _m_radd_ ( m , value ) :
    """ value + matrix 
    """
    if isinstance ( value , num_types ) : value = float ( value ) 
    else                                : return NotImplemented
    ## 
    return m + value

# =============================================================================
## value - matrix 
def _m_rsub_ ( m , value ) :
    """ value - matrix 
    """
    if isinstance ( value , num_types ) : value = float ( value ) 
    else                                : return NotImplemented
    ## 
    return ( m * -1 ) + value 

# =============================================================================
## value * matrix 
def _m_rmul_ ( m , value ) :
    """ value * matrix 
    """
    if isinstance ( value , num_types ) : value = float ( value ) 
    else                                : return NotImplemented
    ## 
    return m * value

Matrix.__iadd__     = _m_iadd_ 
Matrix.__isub__     = _m_isub_ 
Matrix.__imul__     = _m_imul_ 
Matrix.__idiv__     = _m_idiv_ 
Matrix.__itruediv__ = _m_idiv_ 

Matrix.__add__      = _m_add_ 
Matrix.__sub__      = _m_sub_ 
Matrix.__mul__      = _m_mul_ 
Matrix.__div__      = _m_div_ 
Matrix.__truediv__  = _m_div_ 

Matrix.__radd__     = _m_radd_ 
Matrix.__rsub__     = _m_rsub_ 
Matrix.__rmul__     = _m_rmul_ 

# ================================================================================
## vector += value 
def _v_iadd_ ( v , value ) :
    """ vector += value 
    """
    if isinstance ( value , Vector ) :
        if v.size () != value.size () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                              : return NotImplemented
    ##
    v.iadd ( value )
    return v 

# ================================================================================
## vector -= value 
def _v_isub_ ( v , value ) :
    """ vector -= value 
    """
    if isinstance ( value , Vector ) :
        if v.size () != value.size () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                              : return NotImplemented
    ##
    v.isub ( value )
    return v 

# ================================================================================
## vector *= value 
def _v_imul_ ( v , value ) :
    """ vector *= value 
    """
    if isinstance ( value , Matrix ) :
        if v.size () != value.nRows() : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                              : return NotImplemented
    ##
    v.imul ( value )
    return v 

# ================================================================================
## vector /= value 
def _v_idiv_ ( v , value ) :
    """ vector /= value 
    """
    if isinstance ( value , num_types ) : value = float ( value ) 
    else                              : return NotImplemented
    ##
    v.idiv ( value )
    return v 

# =============================================================================
## vector + value 
def _v_add_ ( v , value ) :
    """ vector + value 
    """
    if isinstance ( value , Vector ) :
        if v.size () != value.size () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                              : return NotImplemented
    ##
    result  = Vector ( v )
    result += value 
    return result 

# =============================================================================
## vector - value 
def _v_sub_ ( v , value ) :
    """ vector - value 
    """
    if isinstance ( value , Vector ) :
        if v.size () != value.size () : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                              : return NotImplemented
    ##
    result  = Vector ( v )
    result -= value 
    return result 

# =============================================================================
## vector * value 
def _v_mul_ ( v , value ) :
    """ vector * value 
    """
    if isinstance ( value , Matrix ) :
        if v.size () != value.nRows() : return NotImplemented
    elif isinstance ( value , num_types ) : value = float ( value ) 
    else                              : return NotImplemented
    ##
    result  = Vector ( v )
    result *= value 
    return result 

# =============================================================================
## vector / value 
def _v_div_ ( v , value ) :
    """ vector / value 
    """
    if isinstance ( value , num_types ) : value = float ( value ) 
    else                              : return NotImplemented
    ##
    result  = Vector ( v )
    result /= value 
    return result 

# =============================================================================
## value + vector
def _v_radd_ ( v , value ) :
    """ value + vector  
    """
    if isinstance ( value , num_types ) : value = float ( value ) 
    else                                : return NotImplemented
    ## 
    return v + value

# =============================================================================
## value - vector  
def _v_rsub_ ( v , value ) :
    """ value - vector  
    """
    if isinstance ( value , num_types ) : value = float ( value ) 
    else                                : return NotImplemented
    ## 
    return ( v * -1 ) + value 

# =============================================================================
## value * vector 
def _v_rmul_ ( v , value ) :
    """ value * vector 
    """
    if isinstance ( value , num_types ) : value = float ( value ) 
    else                                : return NotImplemented
    ## 
    return v * value

Vector.__iadd__      = _v_iadd_ 
Vector.__isub__      = _v_isub_ 
Vector.__imul__      = _v_imul_ 
Vector.__idiv__      = _v_idiv_ 
Vector.__itruediv__  = _v_idiv_ 

Vector.__add__       = _v_add_ 
Vector.__sub__       = _v_sub_ 
Vector.__mul__       = _v_mul_ 
Vector.__div__       = _v_div_ 
Vector.__truediv__   = _v_div_ 

Vector.__radd__      = _v_radd_ 
Vector.__rsub__      = _v_rsub_ 
Vector.__rmul__      = _v_rmul_ 


# =============================================================================
## permutation * matrix 
def _p_mul_ ( p , value ) :
    """ permutation * matrix 
    """
    if isinstance ( value , Matrix ) :
        if p.size() != value.nRows() : return NotImplemenetd
    else                             : return NotImplemented 
    return p.apply ( value ) 

Permutation.__mul__ = _p_mul_

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

def _p_pretty_print_ ( mtrx , **kwargs ) :    
    return LA.LinAlgT.P_PRETTY ( mtrx , **kwargs )
def _p_str_          ( mtrx , **kwargs ) :    
    return LA.LinAlgT.P_STR    ( mtrx , **kwargs )

Permutation.pretty_print  = _p_pretty_print_
Permutation.table         = _p_str_ 
Permutation.__str__       = _p_str_ 
Permutation.__repr__      = _p_str_ 

# =============================================================================
# True LinearAlebra stuff
# =============================================================================
## Get (P)LU decomposition of matrix into P,L,U , sch  as \f$  PA = LU \f$, where 
#  - P is permutation
#  - L is lower triangular matrix
#  - U is upper triangular matrix with all diagonal elements equal to 1  
#  @code
#  A = ...
#  P, L, U = A.PLU() 
#  @endcode
def _m_PLU_ ( A ) : 
    """ Get (P)LU decomposition of matrix into P,L,U , sch  as \f$  PA = LU \f$, where 
      - P is permutation
      - L is lower triangular matrix
    - U is upper triangular matrix with all diagonal elements equal to 1  
    >>> A = ...
    >>> P, L, U = A.PLU() 
    """
    M, N = A.nRows() , A.nCols ()
    K    = min ( M , N ) 
    L    = Matrix ( M , K )
    U    = Matrix ( K , N )
    P    = Ostap.GSL.PLU ( A , L , U )
    ##
    return P, L , U 

Matrix.PLU  = _m_PLU_

_new_methods_ = (
    ##
    Matrix.__iadd__          , 
    Matrix.__isub__          , 
    Matrix.__imul__          , 
    Matrix.__idiv__          , 
    Matrix.__itruediv__      , 
    ## 
    Matrix.__add__           , 
    Matrix.__sub__           , 
    Matrix.__mul__           , 
    Matrix.__div__           , 
    Matrix.__truediv__       , 
    ## 
    Matrix.__radd__          , 
    Matrix.__rsub__          , 
    Matrix.__rmul__          ,
    ##
    Vector.__iadd__          , 
    Vector.__isub__          , 
    Vector.__imul__          , 
    Vector.__idiv__          , 
    Vector.__itruediv__      , 
    ## 
    Vector.__add__           , 
    Vector.__sub__           , 
    Vector.__mul__           , 
    Vector.__div__           , 
    Vector.__truediv__       , 
    ## 
    Vector.__radd__          , 
    Vector.__rsub__          , 
    Vector.__rmul__          , 
    ##
    Permutation.__mul__      , 
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
    Matrix.PLU                , 
)
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    

# =============================================================================
##                                                                      The END
# =============================================================================
