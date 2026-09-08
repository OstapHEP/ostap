#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/linalgg.py
#  Few utilities to simplify linear algebra manipulations using GSL
#  - easy-to-use wrappers for vector and matrix classes
#  - easy vector & matrix manipulations and basic operations
# 
#  @see https://www.gnu.org/software/gsl/doc/html/vectors.html
#
#
#  Basic Linear Algebra:
#  - (P)LU       decomposition of general (rectangular) matrix \f$ PA = LU      \f$
#  - (P)QR       decomposition of general (rectangular) matrix \f$ AP = QR      \f$
#  - LQ          decomposition of general (rectangular) matrix \f$ A  = LQ      \f$
#  - QL          decomposition of general (reclangular) matrix \f$ A  = QL      \f$
#  - COD         decomposition of general (rectangular) matrix \f$ AP = QRZ^T   \f$
#  - SVD         decomposition of general (rectangular) matrix \f$ AP = U S V^T \f$
#  - Cholesky    decomposition of symmetric positive-definite matrix \f$ A = L L^T    \f$
#  - Cholesky    decomposition of symmetric positive-definite matrix \f$ PSASP^T = L D L^T    \f$
#  - Tridiagonal decomposition of symmetric matrix     \f$ A = Q D_3 Q^T \f$
#  - Hessenberg  decomposition of square matrix        \f$ A = U H U^T   \f$
#  - Bidiagonalization of general (rectangular) matrix \f$ A = U B_2 U^T \f$
#  - Polar       decomposition of square matrix        \f$ A = UP        \f$
#  - Schur'      decomposition of square matrix        \f$ A = Z S Z^T \f$
#
#  @see https://www.gnu.org/software/gsl/doc/html/linalg.html
#  @see https://www.gnu.org/software/gsl/doc/html/linalg.html#lu-decomposition
#  @see https://www.gnu.org/software/gsl/doc/html/linalg.html#qr-decomposition-with-column-pivoting
#  @see https://www.gnu.org/software/gsl/doc/html/linalg.html#lq-decomposition
#  @see https://www.gnu.org/software/gsl/doc/html/linalg.html#ql-decomposition
#  @see https://www.gnu.org/software/gsl/doc/html/linalg.html#complete-orthogonal-decomposition
#  @see https://www.gnu.org/software/gsl/doc/html/linalg.html#singular-value-decomposition
#  @see https://www.gnu.org/software/gsl/doc/html/linalg.html#cholesky-decomposition
#  @see https://www.gnu.org/software/gsl/doc/html/linalg.html#pivoted-cholesky-decomposition
#  @see https://www.gnu.org/software/gsl/doc/html/linalg.html#tridiagonal-decomposition-of-real-symmetric-matrices
#  @see https://www.gnu.org/software/gsl/doc/html/linalg.html#hessenberg-decomposition-of-real-matrices
#  @see https://www.gnu.org/software/gsl/doc/html/linalg.html#bidiagonalization
#  @see https://www.gnu.org/software/gsl/doc/html/eigen.html#real-generalized-nonsymmetric-eigensystems

#  @see Ostap::Math::GSL
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
""" Few utilities to simplify linear algebra manipulations using GSL

  - easy-to-use wrappers for vector and matrix classes
  - easy vector & matrix manipulations and basic operations 

see https://www.gnu.org/software/gsl/doc/html/vectors.html

  Basic Linear Algebra:
  - (P)LU       decomposition of general (rectangular) matrix         PA = LU     
  - (P)QR       decomposition of general (rectangular) matrix         AP = QR     
  - LQ          decomposition of general (rectangular) matrix         A  = LQ     
  - QL          decomposition of general (reclangular) matrix         A  = QL      
  - COD         decomposition of general (rectangular) matrix         AP = QRZ^T   
  - SVD         decomposition of general (rectangular) matrix         AP = U S V^T 
  - Cholesky    decomposition of symmetric positive-definite matrix   A = L L^T    
  - Cholesky    decomposition of symmetric positive-definite matrix   PSASP^T = L D L^T  
  - Tridiagonal decomposition of symmetric matrix                     A = Q D_3 Q^T 
  - Hessenberg  decomposition of square matrix                        A = U H U^T  
  - Bidiagonalization of general (rectangular) matrix                 A = U B_2 U^T 
  - Polar       decomposition of square matrix                        A = UP        
  - Schur'      decomposition of square matrix                        A = Z S Z^T 

see https://www.gnu.org/software/gsl/doc/html/linalg.html
see https://www.gnu.org/software/gsl/doc/html/linalg.html#lu-decomposition
see https://www.gnu.org/software/gsl/doc/html/linalg.html#qr-decomposition-with-column-pivoting
see https://www.gnu.org/software/gsl/doc/html/linalg.html#lq-decomposition
see https://www.gnu.org/software/gsl/doc/html/linalg.html#ql-decomposition
see https://www.gnu.org/software/gsl/doc/html/linalg.html#complete-orthogonal-decomposition
see https://www.gnu.org/software/gsl/doc/html/linalg.html#singular-value-decomposition
see https://www.gnu.org/software/gsl/doc/html/linalg.html#cholesky-decomposition
see https://www.gnu.org/software/gsl/doc/html/linalg.html#pivoted-cholesky-decomposition
see https://www.gnu.org/software/gsl/doc/html/linalg.html#tridiagonal-decomposition-of-real-symmetric-matrices
see https://www.gnu.org/software/gsl/doc/html/linalg.html#hessenberg-decomposition-of-real-matrices
see https://www.gnu.org/software/gsl/doc/html/linalg.html#bidiagonalization
see https://www.gnu.org/software/gsl/doc/html/eigen.html#real-generalized-nonsymmetric-eigensystems

"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = ""
# =============================================================================
__all__     = ( 
    'Matrix'      , 
    'Vector'      , 
    'Permutation' 
)
# =============================================================================
from   ostap.core.ostap_types import num_types 
from   ostap.math.math_base   import Ostap
from   ostap.utils.gsl        import gsl_info 
import ostap.math.linalg      as     LA 
import ctypes   
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.linalgg' )
else                       : logger = getLogger ( __name__             )
# =============================================================================

Matrix      = Ostap.Math.GSL.Matrix
Vector      = Ostap.Math.GSL.Vector
Permutation = Ostap.Math.GSL.Permutation
Zero        = Matrix.Zero 
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
    elif isinstance ( value , Permutation  ) :
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

Matrix.__imatmul__  = _m_imul_ 
Matrix.__matmul__   = _m_mul_ 

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
        if p.size() != value.nRows() : return NotImplemented
    else                             : return NotImplemented 
    return p.apply ( value ) 

Permutation.__mul__ = _p_mul_

# ==============================================================================
## Convert GSL matrix to SMatrix 
def _to_smatrix_ ( mtrx ) :
    """ Convert GSL matrix to SMatrix 
    """
    ## get dimension 
    nr , nc = mtrx.nRows() , mtrx.nCols()
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
    nr , nc = mtrx.nRows() , mtrx.nCols()
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
    nr , nc = mtrx.nRows() , mtrx.nCols()
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
    nr , nc = mtrx.nRows() , mtrx.nCols()
    if not nr or  nr != nc : raise TypeError ( "Impossible to create symmetric from rectangular matrix!" )     
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
_m_shape_    = lambda m : ( m.nRows(), m.nCols() ) 
Matrix.shape = property ( _m_shape_ , None , None , "`shape` : shape f matrix: (#rows,#columns)" )

def _m_pretty_print_ ( mtrx , **kwargs ) : return LA.LinAlgT.M_PRETTY ( mtrx , **kwargs )
def _m_str_          ( mtrx , **kwargs ) : return LA.LinAlgT.M_STR    ( mtrx , **kwargs )

Matrix.pretty_print  = _m_pretty_print_
Matrix.table         = _m_str_ 
Matrix.__str__       = _m_str_ 
Matrix.__repr__      = _m_str_ 

def _v_pretty_print_ ( mtrx , **kwargs ) : return LA.LinAlgT.V_PRETTY ( mtrx , **kwargs )
def _v_str_          ( mtrx , **kwargs ) : return LA.LinAlgT.V_STR    ( mtrx , **kwargs )

Vector.pretty_print  = _v_pretty_print_
Vector.table         = _v_str_ 
Vector.__str__       = _v_str_ 
Vector.__repr__      = _v_str_ 

def _p_pretty_print_ ( mtrx , **kwargs ) : return LA.LinAlgT.P_PRETTY ( mtrx , **kwargs )
def _p_str_          ( mtrx , **kwargs ) : return LA.LinAlgT.P_STR    ( mtrx , **kwargs )

Permutation.pretty_print  = _p_pretty_print_
Permutation.table         = _p_str_ 
Permutation.__str__       = _p_str_ 
Permutation.__repr__      = _p_str_ 

# =============================================================================
# True Linear Algebra stuff
# =============================================================================

# =============================================================================
## Get (P)LU decomposition of matrix into P,L,U , such  as \f$  PA = LU \f$, where 
#  - P is permutation
#  - L is lower triangular matrix
#  - U is upper triangular matrix with all diagonal elements equal to 1  
#  @code
#  A = ...
#  P, L, U = A.PLU() 
#  @endcode
def _m_PLU_ ( A ) : 
    """ Get (P)LU decomposition of matrix into P,L,U , such  as \f$  PA = LU \f$, where 
      - P is permutation
      - L is lower triangular matrix
    - U is upper triangular matrix with all diagonal elements equal to 1  
    >>> A = ...
    >>> P, L, U = A.PLU() 
    """
    M, N = A.nRows() , A.nCols ()
    K    = min ( M , N ) 
    L    = Matrix ( M , K , Zero () )
    U    = Matrix ( K , N , Zero () )
    P    = Permutation ( M ) 
    sc   = Ostap.Math.GSL.PLU ( A , P , L , U )
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL::PLU" % sc )        
    return P , L , U

# ===============================================================================
## Get (P)QR decompositoon with column pivoting such as  \f$ AP = QR\f$
#  - A is input MxN matrix 
#  - P is permutation (NxN) 
#  - Q is orthogonal MxM matrix
#  - R is right triangular MxN matrix
def _m_PQR_ ( A ) :
    """ Get QR decompositoon with column pivoting such as  AP = QR
    - A is input MxN matrix 
    - P is permutation (NxN) 
    - Q is orthogonal MxM matrix
    - R is right triangular MxN matrix
    
    >>> A = ...
    >>> P, Q, R = A.PQR() 
    """
    M, N = A.nRows() , A.nCols ()
    Q    = Matrix ( M , M )
    R    = Matrix ( M , N , Zero () )
    P    = Permutation ( N ) 
    sc   = Ostap.Math.GSL.PQR( A , P , Q , R )
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL::PQR" % sc )        
    ##
    return P, Q , R 

# ===============================================================================
## Get (P)QR decompositoon with column pivoting such as  \f$ AP = QR\f$
#  - A is input MxN matrix 
#  - P is permutation (NxN) 
#  - Q is orthogonal MxM matrix
#  - R is right triangular MxN matrix
#  - r is the reciprocal condition number of R
def _m_PQRr_ ( A ) :
    """ Get QR decompositoon with column pivoting such as  AP = QR
    - A is input MxN matrix 
    - P is permutation (NxN) 
    - Q is orthogonal MxM matrix
    - R is right triangular MxN matrix
    - r is the reciprocal condition number of R
    >>> A = ...
    >>> P, Q, R, r = A.PQRr () 
    """
    M, N = A.nRows() , A.nCols ()
    Q    = Matrix ( M , M )
    R    = Matrix ( M , N , Zero () )
    P    = Permutation ( N ) 
    r    = ctypes.c_c_double( -1.0 )
    sc   = Ostap.Math.GSL.PQRr ( A , P , Q , R , r )
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL::PQR" % sc )        
    ##
    return P, Q , R, r.value 

# ===============================================================================
## Get LQ decomposition such as  \f$ A = LQ\f$
#  - A is input MxN matrix 
#  - L is lower trapezoidal  MxN matrix
#  - Q is orthogonal NxN matrix
def _m_LQ_ ( A ) :
    """ Get LQ decomposition such as  A = LQ
    - A is input MxN matrix 
    - L is lower trapezoidal  MxN matrix
    - Q is orthogonal NxN matrix    
    >>> A = ...
    >>> L, Q = A.LQ() 
    """
    M, N = A.nRows() , A.nCols ()
    L    = Matrix ( M , N , Zero ()  )
    Q    = Matrix ( N , N )
    sc   = Ostap.Math.GSL.LQ ( A , L , Q )
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL::LQ" % sc )        
    return L , Q 

# ===============================================================================
## Get QL decomposition such as  \f$ A = QL \f$
#  - A is input MxN matrix 
#  - Q is orthogonal MxM matrix
#  - L is lower trapezoidal  MxN matrix
def _m_QL_ ( A ) :
    """ Get QL decomposition with column piviting such as  A = QL
    - A is input MxN matrix 
    - Q is orthogonal MxM matrix    
    - L is lower trapezoidal  MxN matrix
    >>> A = ...
    >>> Q, L = A.QL() 
    """
    M, N = A.nRows() , A.nCols ()
    Q    = Matrix ( M , M )
    L    = Matrix ( M , N , Zero() )
    sc   = Ostap.Math.GSL.QL ( A , Q , L )
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL::QL" % sc )        
    return Q , L

# ===============================================================================
##  COD - Complete Orthogonal Decomposion
#   \f$ AP = Q R Z^T \f$ 
#  - A input MxN matrix 
#  - P is permutation matrix 
#  - Q is MxM orthogonal matrix 
#  - Z is NxN orthogonal matrix 
#  - R is 2x2 block matrix with top-left blobck being right triangular matrix and
#    other blocks are zeroes   
def _m_COD_ ( A ) :
    """ COD - Complete Orthogonal Decomposion: AP = Q R Z^T 
    - A input MxN matrix 
    - P is permutation matrix 
    - Q is MxM orthogonal matrix 
    - Z is NxN orthogonal matrix 
    - R is 2x2 block matrix with top-left blobck being right triangular matrix and
    other blocks are zeroes   
    >>> A = ...
    >>> P , Q , R , Z = A.COD() 
    """
    M, N = A.nRows() , A.nCols ()
    Q    = Matrix      ( M , M )
    R    = Matrix      ( M , N , Zero () )
    Z    = Matrix      ( N , N )
    P    = Permutation ( N ) 
    sc   = Ostap.Math.GSL.COD ( A , P , Q , R , Z )
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL:COD" % sc )            
    return P , Q , R , Z 

# ===============================================================================
## SVD : singular Value Decomposition  \f$ A = U S V^T\f$
#   - A input MxN matrix 
#   - K = min ( M , N ) : 
#   - U MxK orthogonal matrix 
#   - S KxK Diagonal matrix of singular values 
#   - V NxK orthogonal matrix 
#   @param golub (input) use Golub or Jacobi algorithm 
#   @return vector of singular values 
#  -  Jacobi algorithm is more precise  and Golub algorithm is more CPU efficient 
def _m_SVD_ ( A , golub = True ) :
    """ SVD : singular Value Decomposition  \f$ A = U S V^T\f$
    - A input MxN matrix 
    - K = min ( M , N ) : 
    - U MxK orthogonal matrix 
    - S KxK Diagonal matrix of singular values 
    - V NxK orthogonal matrix 
    >>> A = ...
    >>> S , U , V = A.SVD() 
    """
    M, N = A.nRows() , A.nCols ()
    U    = Matrix ( M , N )
    V    = Matrix ( N , N )
    S    = Vector ( min ( M , N ) ) 
    sc   = Ostap.Math.GSL.SVD ( A , S , U , V , True if golub else False )
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL::SVD" % sc )        
    ##
    return S , U , V 

# ===============================================================================
##  LLT: Cholesky decomposition of the square matrix A: \f$ A = L L^T \f$
#  Only lower triangular part of A is used, the upper part is ignored.
#  - A input MxM matrix
#  - L is lower triangular matrix
def _m_LLT_ ( A ) :
    """ LLT: Cholesky decomposition of the square matrix A: A = L L^T
    - A input MxM matrix
    - L is lower triangular matrix
    """
    M, N = A.nRows() , A.nCols ()
    assert M == N , "LLT decomposition is defined only for square matrices!"
    L  = Matrix ( M , M )
    sc = Ostap.Math.GSL.LLT ( A , L )
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL::LLT" % sc )
    ##
    return L


# ===============================================================================
##  LDLT: Cholesky decomposition of the square matrix A: \f$ PSASP^T = L D L^T \f$
#  Only lower triangular part of A is used, the upper part is ignored.
#  - A input MxM matrix
#  - S is scale vector ("diagonal matrix" )
#  - P is permutation  
#  - L is lower triangular matrix
#  - D is vector ("diagonal matrix" )
#  @code
#  A = ...
#  S , P , L , D = A.LDLT()
#  @endcode
def _m_LDLT_ ( A ) :
    """ LLT: Cholesky decomposition of the square matrix A: PSASP^T = L D L^T
    - A input MxM matrix
    - S is scale vector ('diagonal matrix')
    - P is permutation  
    - L is lower triangular matrix
    - D is vector ('diagonal matrix')    
    >>> A  = ...
    >>> S , P , L , D = A.LDLT()
    """
    M, N = A.nRows() , A.nCols ()
    assert M == N , "LDLT decomposition is defined only for square matrices!"
    S  = Vector      ( M )
    P  = Permutation ( M )
    L  = Matrix ( M , M )
    D  = Vector      ( M )
    ## 
    sc = Ostap.Math.GSL.LDLT ( A , S , P , L , D  )
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL::LDLT" % sc )
    ##
    return S , P , L , D 

# ===============================================================================
## D3 : decomposition of symmetric matrix \f$ A = Q D_3 Q^T \f$, where
#  - \f$ Q \f$ is orthogonal matrix
#  - \f$ D_3\fF is symmetric  trigiagonal matrix 
#  @param A (INPUT) input matrix A
#  @param Q (OUTPUT/UPDATE) orthogonal matrix Q
#  @param d (OUTPUT/UPDATE) main diagonal of symmetric  matrix \f$ D_3 \f$
#  @param s (OUTPUT/UPDATE) sub-diagonal of symmetric  matrix \f$ D_3 \f$
#  @return status code
def _m_D3_ ( A ) :
    """ D3 : decomposition of symmetric matrix:  A = Q D_3 Q^T , where
    - Q   is orthogonal matrix
    - D_3 is symmetric  trigiagonal matrix

    Output :
    - Q : orthogonal matrix Q
    - D : symmetric tridiagonal matrix 
    """
    M, N = A.nRows() , A.nCols ()
    assert M == N , "P3 decomposition is defined only for square matrices!"
    assert 2 <= M , "P3 decomposition is defined for 2x2 or larger matrices!" 
    Q = Matrix ( M , M )
    D = Matrix ( M , M )
    ##
    sc = Ostap.Math.GSL.D3 ( A , Q , D )     
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL::D3" % sc )
    ## 
    return Q , D 

# ================================================================================
## Hessenberg decomposition of square matrix \f$ A = U H U^T \f$, where
#  - \f$ U \f$ is orthogonal 
#  - \f$ H \f$ is Hessenberg' matrix: \f$ H(i,i)=0 \f$ for \f$ i > j + 1 \f$
#  @param A (INPUT) input matrix A
#  @param Q (OUTPUT/UPDATE) orthogonal matrix Q
#  @param H (OUTPUT/UPDATE) Hessenberg matrix 
def _m_UHUT_ ( A ) :
    """ Hessenberg decomposition of square matrix: A = U H Q^T, where
    - U   is orthogonal 
    - H  is Hessenberg' matrix:  H(i,i)=0 for i > j + 1 
    >>> A = ...
    >>> U , H = A.UHUT ()
    """
    M, N = A.nRows() , A.nCols ()
    assert M == N , "Hessenberg decomposition is defined only for square matrices!"
    assert 2 <= M , "Hessenberg decomposition is defined for 2x2 or larger matrices!" 
    U = Matrix ( M , M )
    H = Matrix ( M , M )
    ##
    sc = Ostap.Math.GSL.UHUT ( A , U , H )     
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL::UHUT" % sc )
    ## 
    return U , H 
    
# ================================================================================
## Bidiagonalization of of general matrix \f$ A = U B_2 V^T \f$ , where
# - \f$ A \f$ is \f$ M \times N \f$ matrix
# - \f$ U \f$ is \f$ M \times K \f$ orthogonal matrix 
# - \f$ B \f$ is \f$ K \times K \f$ square biadiagonal matrix 
# - \f$ V \f$ is \f$ N \times K \f$ orthogonal matrix
# - K = min ( M , N )
# 
#  @param A (INPUT) input matrix A
#  @param U (OUTPUT/UPDATE) orthogonal matrix U
#  @param B (OUTPUT/UPDATE) bidiagonal matrix B
#  @param V (OUTPUT/UPDATE) orthogonal matrix V
def _m_UBVT_  ( A ) :
    """ Bidiagonalization of of general matrix: A = U B V^T , where
    - A  is  M times N  matrix
    - U  is  M times K  orthogonal matrix
    - B  is  K times K  square biadiagonal matrix 
    - V  is  N times K  orthogonal matrix
    - K = min ( M , N ) 
    >>> A = ...
    >>> U , B , V = A.UBVT()
    """
    
    M, N = A.nRows() , A.nCols ()
    assert 2 <=  min ( M , N )  , "Bidiagonalization is defined for 2x2 or larger matrices!"

    K = min ( M , N )
    
    U = Matrix ( M , K )
    B = Matrix ( K , K )
    V = Matrix ( N , K )
    
    sc = Ostap.Math.GSL.UBVT ( A , U , B , V ) 
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL::UBVT" % sc )
    ## 
    return U , B , V 

# ===============================================================================
## Polar decompositon of the square matrix A: \f$ A = UP \f$
#  - U is orthogonal 
#  - P is positive semi-definitive 
def _m_POLAR_ ( A ) :
    """ Polar decomposition of the square matrix A: A = UP
    - U is orthogonal 
    - P is positive semi-definitive 
    >>> A = ...,
    >>> U , P = A.POLAR() 
    """
    M, N = A.nRows() , A.nCols ()
    assert M == N , "Polar decomposition is defined only for square matrices!"
    U  = Matrix ( M , M )
    P  = Matrix ( M , M )
    sc = Ostap.Math.GSL.POLAR ( A , U , P )
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL::POLAR!" % sc )
    ## 
    return U , P 

# ===============================================================================
## Schur decompositon of the square matrix A: \f$ A = Z T z^t \f$
#  - Z is orthogonal 
#  - T is a Schur form
#  @see https://www.gnu.org/software/gsl/doc/html/eigen.html#real-generalized-nonsymmetric-eigensystems
def _m_SCHUR_ ( A ) :
    """ Schur decomposition of the square matrix A: A = Z T Z^T
    - Z is orthogonal 
    - T is a Schur forms  
    >>> A = ...,
    >>> Z , T = A.SCHUR ()
    
    see https://www.gnu.org/software/gsl/doc/html/eigen.html#real-generalized-nonsymmetric-eigensystems    
    """
    M, N = A.nRows() , A.nCols ()
    assert M == N , "Schur decomposition is defined only for square matrices!"
    Z  = Matrix ( M , M )
    T  = Matrix ( M , M )
    sc = Ostap.Math.GSL.SCHUR ( A , Z , T )
    if sc.isFailure () : raise ValueError ( "Error code %s from Ostap::Math::GSL::SCHUR" % sc )    
    return Z , T  

Matrix.PLU       = _m_PLU_
Matrix.PQR       = _m_PQR_ 
Matrix.PQRr      = _m_PQRr_
Matrix.LQ        = _m_LQ_ 
Matrix.QL        = _m_QL_ 
Matrix.COD       = _m_COD_
Matrix.SVD       = _m_SVD_
Matrix.LLT       = _m_LLT_
Matrix.LDLT      = _m_LDLT_
Matrix.D3        = _m_D3_
Matrix.UHUT      = _m_UHUT_
Matrix.UBVT      = _m_UBVT_

Matrix.POLAR     = _m_POLAR_ 
Matrix.SCHUR     = _m_SCHUR_  

Matrix.t         = Matrix.T
Matrix.transpose = Matrix.T

_new_methods_ = () 

if  ( 2 , 7 ) <= gsl_info : 
    Matrix.QL = _m_QL_ 
    _new_methods_ = Matrix.QL ,  

_new_methods_ += (
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
    Matrix.t                  , 
    Matrix.transpose          , 
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
    #
    Matrix.PQR                , 
    Matrix.PQRr               ,
    #
    Matrix.LQ                 ,
    Matrix.QL                 ,
    ## 
    Matrix.COD                ,
    Matrix.SVD                ,
    ## 
    Matrix.LLT                ,
    Matrix.LDLT               ,
    Matrix.D3                 ,
    Matrix.UHUT               ,
    Matrix.UBVT               ,
    ## 
    Matrix.SCHUR              , 
    Matrix.POLAR              ,
    ##
)

# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    

# =============================================================================
##                                                                      The END
# =============================================================================
