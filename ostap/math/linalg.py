#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/linalg.py
#  Few utilities to simplify linear algebra manipulations 
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
"""Few utilities to simplify linear algebra manipulations 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = ""
# =============================================================================
__all__     = () ## nothing to be imported !
# =============================================================================
import ROOT, cppyy 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.linalg' )
else                       : logger = getLogger ( __name__            )
# =============================================================================
from ostap.math.base  import isequal,iszero
from ostap.core.types import num_types, is_integer

cpp   = cppyy.gbl

## C++ namespace std
std   = cpp.std

## C++ namespace Ostap 
Ostap = cpp.Ostap

## ROOT::Math namespace
_RM = ROOT.ROOT.Math

a   = Ostap.Math.ValueWithError()
# =============================================================================
## try to pickup the vector
@staticmethod
def _vector_ ( i , typ = 'double' ) :
    """Pick up the vector of corresponding size
    >>> V3   = Ostap.Math.Vector(3)
    >>> vct  = V3 ()
    """
    assert is_integer ( i ) and 0 <= i , 'Invalid vector size %s' % i
    v = _RM.SVector ( typ , i )
    return deco_vector ( v ) 

# =============================================================================
## try to pickup the matrix
@staticmethod
def _matrix_ ( i , j , typ = 'double' ) :
    """Pick up the matrix of corresponding size
    >>> M3x4   = Ostap.Math.Matrix(3,4)
    >>> matrix = M3x4 ()    
    """
    assert is_integer ( i ) and 0 <= i , 'Invalid matrix size (%s,%s)' % ( i , j )
    assert is_integer ( j ) and 0 <= j , 'Invalid matrix size (%s,%s)' % ( i , j )
    m = _RM.SMatrix ( "%s,%d,%d" % ( typ , i , j ) )
    return deco_matrix( m )  

# =============================================================================
## try to pickup the symmeric matrix
@staticmethod
def _sym_matrix_ ( i , typ = 'double' ) :
    """Pick up the symmetric matrix of corresponding size
    >>> SymM3  = Ostap.Math.SymMatrix(3)
    >>> matrix = SymM3 ()
    """
    assert is_integer ( i ) and 0 <= i , 'Invalid matrix size %s' %  i 
    m = _RM.SMatrix('%s,%d,%d,ROOT::Math::MatRepSym<%s,%d>' %  ( typ , i , i , typ , i ) )
    return deco_symmatrix  ( m ) 

Ostap.Vector         =     _vector_
Ostap.Math.Vector    =     _vector_
Ostap.Matrix         =     _matrix_
Ostap.Math.Matrix    =     _matrix_
Ostap.SymMatrix      = _sym_matrix_
Ostap.Math.SymMatrix = _sym_matrix_

# =============================================================================
## self-printout of S-vectors
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
def _v_str_ ( self , fmt = ' %g' ) :
    """Self-printout of SVectors: (...)
    """
    index  = 0
    result = ''
    while index < self.kSize :
        if 0 != index : result += ', '
        result += fmt % self.At( index )
        index  += 1
    return "( " + result + ' )'

# =============================================================================
## iterator for SVector
#  @code
#  vct = ...
#  for i in vct : print i 
#  @endcode
def _v_iter_ ( self ) :
    """Iterator for SVector
    >>> vct = ...
    >>> for i in vct : print i 
    """
    for i in range(self.kSize) :
        yield self(i)
        
# =============================================================================
## iterator for SVector
#  @code
#  vct = ...
#  for i,v in vct.iteritems() : print i,v 
#  @endcode
def _v_iteritems_ ( self ) :
    """Iterator for SVector
    >>> vct = ...
    >>> for i,v in vct.iteritems() : print i,v 
    """
    for i in range(self.kSize) :
        yield i,self(i)

# =============================================================================
## convert vector to numpy array:
#  @code
#  vct = ...
#  na = vct.to_numpy() 
#  @endcode
def _v_to_numpy_ ( self ) :
    """Convert vector to numpy array:
    >>> vct = ...
    >>> na = vct.to_numpy() 
    """
    import numpy
    s  = self.kSize
    na = numpy.empty( s )
    for i in range(s) : na[i] = self[i] 
    return na

# =============================================================================
## convert vector to plain array:
#  @code
#  vct = ...
#  aa  = vct.to_array() 
#  @endcode
def _v_to_array_ ( self ) :
    """Convert vector to plain array:
    >>> vct = ...
    >>> na  = vct.to_array() 
    """
    import array
    return array.array('f', self)

# =============================================================================
## the multiplication operators 
_mult_ops_ = {}
## get the proper multiplication operator 
def _get_mult_op_ ( klass1 , klass2 ) :
    """Get the proper multiplication operator
    """
    t   = klass1 , klass2
    ops = _mult_ops_.get( t , None )
    if ops : return ops                   ## RETURN  

    ## try to load the operators 
    try :
        ops = Ostap.Math.MultiplyOp ( klass1 , klass2 )
        _mult_ops_ [ t ] = ops
        return ops                       ## RETURN 
    except TypeError:
        return None                      ## RETURN
    
    return None                          ## RETURN


# =============================================================================
## equailty operators 
_eq_ops_ = {}
## get the proper equality operator 
def _get_eq_op_ ( klass1 , klass2 ) :
    """Get the proper equality operator
    """
    t   = klass1 , klass2
    ops = _eq_ops_.get( t , None )
    if ops : return ops                   ## RETURN  

    ## try to load the operators 
    try :
        ops = Ostap.Math.EqualityOp ( klass1 , klass2 )
        _eq_ops_ [ t ] = ops
        return ops                       ## RETURN 
    except TypeError:
        return None                      ## RETURN
    
    return None                          ## RETURN


# =============================================================================
## helper function for Linear Algebra: multiplications
#  multiplication of vectors, matrices, constants 
#  @code
#  vector1 = ...
#  vector2 = ...
#  matrix1 = ...
#  matrix2 = ...
#  print vector1 * vector2 
#  print vector1 * matrix1
#  print matrix1 * vector2
#  print matrix1 * matrix2
#  print vector1 * 2
#  print matrix1 * 2
#  @endcode
def _linalg_mul_ ( a  , b ) :
    """Multiplication of vectors, matrices, etc
    >>> vector1 = ...
    >>> vector2 = ...
    >>> matrix1 = ...
    >>> matrix2 = ...
    >>> print vector1 * vector2 
    >>> print vector1 * matrix1
    >>> print matrix1 * vector2
    >>> print matrix1 * matrix2
    >>> print vector1 * 2
    >>> print matrix1 * 2
    """
    ## simple cases: multiply by a constant 
    if isinstance ( b , num_types ) :
        b  = float( b )
        v  = a.__class__( a )
        v *= b
        return v
    
    ## get the proper operator 
    ops   = _get_mult_op_ ( a.__class__ , b.__class__ )
    if not ops : return NotImplemented 
    return ops.multiply ( a , b )

# =============================================================================
## helper function for "right" multiplication (Linear Algebra)
#  "right multiplication" for a constant
#  @code
#  vector = ...
#  matrix = ...
#  print 2 * vector
#  print 2 * matrix 
#  @endcode 
def _linalg_rmul_ ( a , b ) :
    """``right multiplication'' for a constant
    >>> vector = ...
    >>> matrix = ...
    >>> print 2 * vector
    >>> print 2 * matrix
    """
    if isinstance ( b , num_types ) :
        b  = float( b )
        v  = a.__class__( a )
        v *= b
        return v
    return NotImplemented 

# =============================================================================
## helper function for Linear Algebra divisions 
#  Division by a constant
#  @code
#  vector = ...
#  matrix = ...
#  print vector / 2 
#  print matrix / 2 
#  @endcode 
def _linalg_div_ ( a  , b ) :
    """Division by a constant
    >>> vector = ...
    >>> matrix = ...
    >>> print vector / 2 
    >>> print matrix / 2 
    """
    if isinstance ( b , num_types ) :
        b  = float( b )
        v  = a.__class__( a )
        v /= b
        return v
    return NotImplemented

# =============================================================================
## "cross-product" of two vectors to get a matrix
#  @code 
#  vector1 = ...
#  vector2 = ...
#  matrix =  vector1.cross ( vector2 )
#  @endcode 
def _vector_cross_ ( a, b ) :
    """Cross-product of two vectors to get a matrix
    >>> vector1 = ...
    >>> vector2 = ...
    >>> matrix =  vector1.cross ( vector2 ) 
    """
    ## get the proper operator 
    ops   = _get_mult_op_ ( a.__class__ , b.__class__ )
    if not ops : return NotImplemented 
    return ops.cross ( a , b )


# =============================================================================
## equality of vectors
#  @code
#  vector1 = ...
#  vector2 = ...
#  print vector1 == vector2
#  print vector1 == ( 0, 2, 3 )
#  print vector1 == [ 0, 2, 3 ]
#  @endcode 
def _vector_eq_ ( a , b ) :
    """Equality for vectors
    >>> vector1 = ...
    >>> vector2 = ...
    >>> print vector1 == vector2
    >>> print vector1 == ( 0, 2, 3 )
    >>> print vector1 == [ 0, 2, 3 ]
    """
    if         a    is      b          : return True
    elif not hasattr ( b , '__len__' ) : return False 
    elif  len ( a ) != len ( b )       : return False        
    #
    ops = _get_eq_op_ ( a.__class__ , b.__class__ )
    if ops : return ops.equal ( a , b )  ## RETURN
    ## compare elements  
    for i in range ( len( a ) ) :
        if not isequal ( a[i] , b[i] ) : return False
        
    return True 
    
# =============================================================================
## equality of matrices
#  @code
#  matrix1 = ...
#  matrix2 = ...
#  print matrix1 == matrix2 
#  @endcode 
def _matrix_eq_ ( a , b ) :
    """Equality for matrices
    >>> matrix1 = ...
    >>> matrix2 = ...
    >>> print matrix1 == matrix2 
    """
    if  a is b : return True
    
    try :
        if   a.kRows != b.kRows : return False
        elif a.kCols != b.kCols : return False
    except :
        pass
        
    ops = _get_eq_op_ ( a.__class__ , b.__class__ )
    if not ops : return NotImplemented
    return ops.equal ( a , b )

# =============================================================================
## decorate vector 
def deco_vector ( t ) :

    if not hasattr ( t , '_decorated' ) :

        t ._old_str_    = t.__str__
        t ._old_repr_   = t.__repr__
        
        t ._old_add_    = t.__add__
        t ._old_radd_   = t.__radd__
        t ._old_mul_    = t.__mul__
        t ._old_rmul_   = t.__rmul__
        t ._old_sub_    = t.__sub__
        t ._old_rsub_   = t.__rsub__
        t ._old_div_    = t.__div__

        _operations   = Ostap.Math.VctrOps( t )
        
        t.__add__       = lambda a,b : _operations.add  ( a , b )
        t.__sub__       = lambda a,b : _operations.sub  ( a , b )
        t.__radd__      = lambda a,b : _operations.add  ( a , b )
        t.__rsub__      = lambda a,b : _operations.rsub ( a , b )
        
        t.__mul__       = _linalg_mul_    
        t.__rmul__      = _linalg_rmul_    
        t.__div__       = _linalg_div_    
        
        t.__eq__        = _vector_eq_    
        t.__neq__       = lambda a,b : not ( a == b )
        t.__neg__       = lambda s : s*(-1)
        
        t.cross         = _vector_cross_        
        t.__rdiv__      = lambda s,*a :  NotImplemented 

        t. _new_str_    = _v_str_
        t. __str__      = _v_str_
        t. __repr__     = _v_str_

        t. __len__      = lambda s : s.kSize 
        t. __contains__ = lambda s, i : 0<=i<s.kSize

        t. __iter__     = _v_iter_        
        t. iteritems    = _v_iteritems_

        t.to_numpy      = _v_to_numpy_
        t.to_array      = _v_to_array_
        
        t. _decorated = True
        
    return t


Ostap.Vector2             = Ostap.Vector(2)
Ostap.Vector3             = Ostap.Vector(3)
Ostap.Vector4             = Ostap.Vector(4)
Ostap.Vector5             = Ostap.Vector(5)
Ostap.Vector6             = Ostap.Vector(6)
Ostap.Vector8             = Ostap.Vector(8)

Ostap.Math.Vector2        = Ostap.Vector2
Ostap.Math.Vector3        = Ostap.Vector3
Ostap.Math.Vector4        = Ostap.Vector4
Ostap.Math.Vector5        = Ostap.Vector5
Ostap.Math.Vector6        = Ostap.Vector6
Ostap.Math.Vector8        = Ostap.Vector8

## vectors of vectors
Ostap.Vectors2            = std.vector ( Ostap.Vector2 )
Ostap.Vectors3            = std.vector ( Ostap.Vector3 )
Ostap.Vectors4            = std.vector ( Ostap.Vector4 )
Ostap.Math.Vectors2       = Ostap.Vectors2
Ostap.Math.Vectors3       = Ostap.Vectors3
Ostap.Math.Vectors4       = Ostap.Vectors4

# ============================================================================
## self-printout of matrices
def _mg_str_ ( self , fmt = ' %+11.4g') :
    """Self-printout of matrices
    >>> matrix = ...
    >>> print matrix 
    """
    _rows = self.kRows
    _cols = self.kCols
    _line = ''
    for _irow in range ( 0 , _rows ) :
        _line += ' |'
        for _icol in range ( 0 , _cols ) :
            _line += fmt % self( _irow , _icol )
        _line += ' |'
        if ( _rows - 1 )  != _irow : _line += '\n'
    return _line

# =============================================================================
## self-printout of symmetrical matrices
def _ms_str_ ( self , fmt = ' %+11.4g' , width = 12 ) :
    """Self-printout of symmetrical matrices
    >>> matrix = ...
    >>> print matrix 
    """
    _rows = self.kRows
    _cols = self.kCols
    _line = ''
    for _irow in range ( 0 , _rows ) :
        _line += ' |'
        for _icol in range ( 0 , _cols  ) :
            if _icol < _irow : _line += width*' '
            else             : _line += fmt % self( _irow , _icol )
        _line += ' |'
        if ( _rows - 1 )  != _irow : _line += '\n'
    return _line

# =============================================================================
## get the correlation matrix
def _ms_corr_ ( self ) :
    """Get the correlation matrix
    >>> mtrx = ...
    >>> corr = mtrx.correlations()
    """
    from math import sqrt

    _t    = type ( self )
    _c    = _t   ()
    _rows = self.kRows
    _ok1  = True 
    _ok2  = True 
    for i in range ( 0 , _rows ) :
        
        ii  = self  ( i , i )
        if 0 > ii or iszero ( ii ) :
            _ok1  = False 
            _nan = float('nan')
            for j in range ( i , _rows ) : _c [ i , j ] = _nan 
            continue
        
        sii = sqrt ( ii )
        _c[ i , i ] = 1
        
        for j in range ( i + 1 , _rows ) :            
            jj  = self ( j , j ) 
            sjj = sqrt ( jj )
            ij  = self ( i , j )
            eij = ij / ( sii * sjj )
            if  1 < abs ( eij ) : _ok2 = False  
            _c [ i , j ] = eij 
            
    if not _ok1 : logger.error ( "correlations: zero or negative diagonal element" ) 
    if not _ok2 : logger.error ( "correlations: invalid non-diagonal element"      ) 
        
    return _c

# =============================================================================
## "getter"
def _m_get_ ( o , i , j ) :

    try :
        return o ( i , j )
    except :
        pass
    
    try :
        return o [ i , j ]
    except :
        pass

    return o [ i ][ j ]
        
# =============================================================================
## add some matrix-like object to the matrix
#  @code
#  m = ...
#  o = ...
#  m.add_to  ( o ) 
#  @endcode
def _mg_increment_ ( m , o ) : 
    """ Add some ``matrix-like'' object to the matrix
    >>> m = ...
    >>> o = ...
    >>> m.increment  ( o ) 
    """
    for i in range ( m.kRows ) :
        for j in range ( m.kCols ) :
            m[i,j] = m(i,j) + _m_get_ ( o , i , j )
                    
    return m


# =============================================================================
## add some matrix-like object to symmetric matrix, preserving the symmetry 
#  @code
#  m = ...
#  o = ...
#  m.add_to  ( o ) 
#  @endcode
def _ms_increment_ ( m , o ) : 
    """ Add some ``matrix-like'' object to the matrix
    >>> m = ...
    >>> o = ...
    >>> m.increment  ( o ) 
    """
    for i in range ( m.kRows ) :
        for j in range ( i , m.kCols ) :
            m[i,j] = m(i,j) + 0.5 * ( _m_get_ ( o , i , j ) + _m_get_ ( o , j , i ) )  
                    
    return m

# =============================================================================
## iterator for SMatrix
#  @code
#  matrix = ...
#  for i in matrix : print i 
#  @endcode
def _m_iter_ ( self ) :
    """Iterator for SMatrix
    >>> matrix = ...
    >>> for i in matrix : print i 
    """
    for i in range(self.kRows) :
        for j in range(self.kCols) :
            yield self(i,j)

# =============================================================================
## iterator for SMatrix
#  @code
#  matrix = ...
#  for i,j,v in matrix.iteritems() : print i,j,v 
#  @endcode
def _m_iteritems_ ( self ) :
    """Iterator for SMatrix
    >>> matrix = ...
    >>> for i,j,v in matrix.iteritems() : print i,j,v
    """
    for i in range(self.kRows) :
        for j in range(self.kCols) :
            yield i,j,self(i,j)

# =============================================================================
## convert matrix into numpy.matrix
#  @code
#  mtrx = ...
#  m    = mtrx.to_numpy() 
#  @endcode 
def _m_to_numpy_ ( self ) :
    """Convert matrix into numpy.matrix
    >>> mtrx = ...
    >>> m    = mtrx.to_numpy() 
    """
    import numpy
    n = numpy.empty( [ self.kRows, self.kCols ] )
    m = numpy.matrix ( n )
    for i in range(self.kRows) :
        for j in range(self.kCols) :
            m [i,j] = self(i,j)
    return m

# =============================================================================
## construct ``similarity'' with ``vector-like'' object
#  @code
#  m = ...
#  v = ...
#  m.sim ( v )
#  @endcode 
def  _ms_sim_ ( m , v ) :
    """ construct ``similarity'' with ``vector-like'' object
    >>> m = ...
    >>> v = ...
    >>> m.sim ( v )
    """
    sim = 0.0
    for i in range(m.kRows) :
        for j in range(m.kCols) :
            sim += v[i]*m(i,j)*v[j] 
    return sim 


# =============================================================================
## get the eigenvalues for symmetric matrices :
def _eigen_1_ ( self , sorted = True ) :
    """Get the eigenvalues for symmetric matrices :
    >>> mtrx = ...
    >>> values = mtrx.eigenValues ( sorted = True )
    """
    return Ostap.Math.EigenSystems.eigenValues ( self , sorted )

# =============================================================================
## get the eigevectors for symmetric matrices :
def _eigen_2_ ( self , sorted = True ) :
    """Get the eigevectors for symmetric matrices :
    >>> mtrx = ...
    >>> values, vectors = mtrx.eigenVectors( sorted = True )
    """
    if   2 == self.kCols :
        _values  = Ostap.Vector2  ()
        _vectors = Ostap.Vectors2 ()
    elif 3 == self.kCols :
        _values  = Ostap.Vector3  ()
        _vectors = Ostap.Vectors3 ()
    elif 4 == self.kCols :
        _values  = Ostap.Vector4  ()
        _vectors = Ostap.Vectors4 ()
    else :
        raise AttributeError, "Not implemented for dimention: %s" % self.kCols

    st = Ostap.Math.EigenSystems.eigenVectors ( self , _values , _vectors , sorted )
    if st.isFailure () :
        print 'EigenVectors: Failure from EigenSystems' , st

    return ( _values , _vectors )


_eigen_1_ .__doc__ += '\n' +  Ostap.Math.EigenSystems.eigenValues  . __doc__
_eigen_2_ .__doc__ += '\n' +  Ostap.Math.EigenSystems.eigenVectors . __doc__

            
# =============================================================================
##  decorate the matrix  type 
def deco_matrix ( m  ) :
    
    if not hasattr ( m , '_decorated' ) :

        ##  save 'old method'
        m. _old_str_   = m . __str__
        m. _old_repr_  = m . __repr__

        m. _new_str_   = _mg_str_
        m. __repr__    = _mg_str_
        m. __str__     = _mg_str_
        
        m ._old_add_    = m.__add__
        m ._old_radd_   = m.__radd__
        m ._old_mul_    = m.__mul__
        m ._old_rmul_   = m.__rmul__
        m ._old_sub_    = m.__sub__
        m ._old_rsub_   = m.__rsub__
        m ._old_div_    = m.__div__

        _operations     = Ostap.Math.MtrxOps( m )
        
        m.__add__       = lambda a,b : _operations.add  ( a , b )
        m.__sub__       = lambda a,b : _operations.sub  ( a , b )
        m.__radd__      = lambda a,b : _operations.add  ( a , b )
        m.__rsub__      = lambda a,b : _operations.rsub ( a , b )
        
        m.__mul__       = lambda a,b : _linalg_mul_     ( a , b )
        m.__rmul__      = lambda a,b : _linalg_rmul_    ( a , b )
        m.__div__       = lambda a,b : _linalg_div_     ( a , b )
        
        m.__rdiv__      = lambda s,*a :  NotImplemented 

        m.__eq__        = lambda a,b : _matrix_eq_      ( a , b )
        m.__neq__       = lambda a,b : not ( a == b ) 
        m.__neg__       = lambda s   : s*(-1)

        
        m._increment_  = _mg_increment_
        m. increment   = _mg_increment_

        m.__iter__     = _m_iter_
        m.iteritems    = _m_iteritems_
        m.__contains__ = lambda s,ij : 0<=ij[0]<s.kRows and 0<=ij[1]<s.kCols

        m.to_numpy     = _m_to_numpy_ 
        m._decorated   = True
        
    return m

# =============================================================================
##  decorate the symmetrix matrix  type 
def deco_symmatrix ( m ) :

    if not hasattr ( m , '_decorated' ) :

        ##  save 'old method'
        m. _old_str_   = m . __str__
        m. _old_repr_  = m . __repr__

        m.correlations = _ms_corr_
        m. _new_str_   = _ms_str_
        m. __repr__    = _ms_str_
        m. __str__     = _ms_str_

        m ._old_add_    = m.__add__
        m ._old_radd_   = m.__radd__
        m ._old_mul_    = m.__mul__
        m ._old_rmul_   = m.__rmul__
        m ._old_sub_    = m.__sub__
        m ._old_rsub_   = m.__rsub__
        m ._old_div_    = m.__div__

        _operations     = Ostap.Math.MtrxOps( m )
        
        m.__add__       = lambda a,b : _operations.add  ( a , b )
        m.__sub__       = lambda a,b : _operations.sub  ( a , b )
        m.__radd__      = lambda a,b : _operations.add  ( a , b )
        m.__rsub__      = lambda a,b : _operations.rsub ( a , b )
        
        m.__mul__       = _linalg_mul_  
        m.__rmul__      = _linalg_rmul_  
        m.__div__       = _linalg_div_    
        
        m.__rdiv__      = lambda s,*a :  NotImplemented 

        m.__eq__        = _matrix_eq_    
        m.__neq__       = lambda a,b : not ( a == b ) 
        m.__neg__       = lambda s   : s*(-1) 
        
        m._increment_  = _ms_increment_
        m.increment    = _ms_increment_

        m.__iter__     = _m_iter_
        m.iteritems    = _m_iteritems_
        m.__contains__ = lambda s,ij : 0<=ij[0]<s.kRows and 0<=ij[1]<s.kCols
        
        m._sim_        = _ms_sim_
        m.sim          = _ms_sim_

        if m.kCols < 4 : 
            m.eigenValues  = _eigen_1_
            m.eigenVectors = _eigen_2_
        
        m.to_numpy     = _m_to_numpy_
        m._decorated   = True
        
    return m

Ostap.SymMatrix2x2        = Ostap.SymMatrix(2)
Ostap.SymMatrix3x3        = Ostap.SymMatrix(3)
Ostap.SymMatrix4x4        = Ostap.SymMatrix(4)
Ostap.SymMatrix5x5        = Ostap.SymMatrix(5)
Ostap.SymMatrix6x6        = Ostap.SymMatrix(6)
Ostap.SymMatrix7x7        = Ostap.SymMatrix(7)
Ostap.SymMatrix8x8        = Ostap.SymMatrix(8)
Ostap.SymMatrix9x9        = Ostap.SymMatrix(9)


Ostap.Math.SymMatrix2x2   = Ostap.SymMatrix2x2
Ostap.Math.SymMatrix3x3   = Ostap.SymMatrix3x3
Ostap.Math.SymMatrix4x4   = Ostap.SymMatrix4x4
Ostap.Math.SymMatrix5x5   = Ostap.SymMatrix5x5
Ostap.Math.SymMatrix6x6   = Ostap.SymMatrix6x6
Ostap.Math.SymMatrix7x7   = Ostap.SymMatrix7x7
Ostap.Math.SymMatrix8x8   = Ostap.SymMatrix8x8
Ostap.Math.SymMatrix9x9   = Ostap.SymMatrix9x9

for i in range(11) :
    
    t0 = Ostap.Vector ( i )
    deco_vector    ( t0 )
    
    t1 = Ostap.SymMatrix(i)
    deco_symmatrix ( t1 )
    
    for j in range(11) : 
        t2 = Ostap.Matrix(i,j)
        deco_matrix ( t2 ) 

# =============================================================================
_decorated_classes_ = (
    Ostap.Vector(2)     ,
    Ostap.Matrix(2,3)   , 
    Ostap.SymMatrix(2)  , 
    )

_new_methods_ = (
    Ostap.Vector    , 
    Ostap.Matrix    , 
    Ostap.SymMatrix , 
    )
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    logger.info('TEST vectors: ')
    
    LA3 = Ostap.Math.Vector3
    l1  = LA3(0,1,2)
    l2  = LA3(3,4,5)
    
    logger.info ( 'l1 , l2 : %s %s '  % ( l1 , l2  ) )
    logger.info ( 'l1 + l2 : %s    '  % ( l1 + l2  ) )
    logger.info ( 'l1 - l2 : %s    '  % ( l1 - l2  ) )
    logger.info ( 'l1 * l2 : %s    '  % ( l1 * l2  ) )
    logger.info ( 'l1 *  2 : %s    '  % ( l1 *  2  ) )
    logger.info ( ' 2 * l2 : %s    '  % ( 2  * l2  ) )
    logger.info ( 'l1 +  2 : %s    '  % ( l1 +  2  ) )
    logger.info ( ' 2 + l2 : %s    '  % ( 2  + l2  ) )
    logger.info ( 'l1 -  2 : %s    '  % ( l1 -  2  ) )
    logger.info ( ' 2 - l2 : %s    '  % ( 2  - l2  ) )
    logger.info ( 'l1 /  2 : %s    '  % ( l1 /  2  ) )
    
    l1 /= 2 
    logger.info ( 'l1 /= 2 : %s    '  % l1 )
    l1 *= 2 
    logger.info ( 'l1 *= 2 : %s    '  % l1 )
    l1 += 2 
    logger.info ( 'l1 += 2 : %s    '  % l1 )
    l1 -= 2 
    logger.info ( 'l1 -= 2 : %s    '  % l1 )

    logger.info('TEST matrices: ')
    
    m22 = Ostap.Math.Matrix(2,2) ()
    m23 = Ostap.Math.Matrix(2,3) ()
    s22 = Ostap.Math.SymMatrix(2)()
    
    l2  = Ostap.Math.Vector(2)()
    l3  = Ostap.Math.Vector(3)()
    
    l2[0]    = 1
    l2[1]    = 2
    
    l3[0]    = 1
    l3[1]    = 2
    l3[1]    = 3
    
    logger.info ( 'l2 , l3 : %s %s '  % ( l2 , l3  ) )
    
    m22[0,0] = 1
    m22[0,1] = 1
    m22[1,1] = 1
    
    m23[0,0] = 1
    m23[1,1] = 1
    m23[0,2] = 1
    
    s22[0,0] = 2
    s22[1,0] = 1
    s22[1,1] = 3
    
    logger.info ( 'm22\n%s'   % m22     ) 
    logger.info ( 's22\n%s'   % s22     ) 
    logger.info ( 'm23\n%s'   % m23     ) 
    logger.info ( 'm22/3\n%s' % (m22/3) ) 
    logger.info ( 'm23*3\n%s' % (m23*3) ) 

    logger.info ( 'm22 * m23 :\n%s' % ( m22 * m23 ) ) 
    logger.info ( 'm22 *  l2 : %s ' % ( m22 * l2  ) ) 
    logger.info ( 'l2  * m22 : %s ' % ( l2  * m22 ) ) 
    logger.info ( 'm23 *  l3 : %s ' % ( m23 * l3  ) ) 
    logger.info ( 'l2  * m23 : %s ' % ( l2  * m23 ) )
    
    logger.info ( 'm22 * s22 + 2 * m22 :\n%s ' %  ( m22*s22 + 2*m22  ) )
    logger.info ( 'm22 == m22*1.0 : %s ' % (  m22 == m22 * 1.0 ) )
    logger.info ( 'm22 != m22*1.1 : %s ' % (  m22 != m22 * 1.1 ) )
    logger.info ( 'm23 == m23*1.0 : %s ' % (  m23 == m23 * 1.0 ) )
    logger.info ( 'm23 != m23*1.1 : %s ' % (  m23 != m23 * 1.1 ) )
    logger.info ( 'l1  == l1 *1.0 : %s ' % (  l1  == l1  * 1.0 ) )
    logger.info ( 'l1  != l1 *1.1 : %s ' % (  l1  != l1  * 1.1 ) )
    logger.info ( 's22 == s22*1.0 : %s ' % (  s22 == s22 * 1.0 ) )
    logger.info ( 's22 != s22*1.1 : %s ' % (  s22 != s22 * 1.1 ) )
    
    logger.info ( ' l1 == (0,1,2) : %s ' % (  l1 == ( 0 , 1 , 2 ) ) )
    logger.info ( ' l1 == [0,1,2] : %s ' % (  l1 == [ 0 , 1 , 2 ] ) )
    

    m22[0,0] = 1
    m22[0,1] = 2
    m22[1,0] = 2
    m22[1,1] = 3
    
    s22[0,0] = 1
    s22[0,1] = 2
    s22[1,1] = 3
    
    logger.info ( ' m22 == s22     : %s ' % ( m22 == s22       ) )
    logger.info ( ' m22 == s22*1.0 : %s ' % ( m22 == s22 * 1.0 ) )
    logger.info ( ' m22 != s22*1.1 : %s ' % ( m22 != s22 * 1.1 ) )

# =============================================================================
# The END 
# =============================================================================
