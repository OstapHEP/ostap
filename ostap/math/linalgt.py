#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/linalgt.py
#  Few utilities to simplify linear algebra manipulations 
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
"""Few utilities to simplify linear algebra manipulations 
"""
# =============================================================================
from   __future__        import print_function
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = ""
# =============================================================================
__all__     = ( )
# =============================================================================
from   builtins import range 
import ROOT, re 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.linalgt' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
from ostap.math.base        import isequal , iszero , std, Ostap
from ostap.core.ostap_types import num_types, integer_types
from ostap.utils.clsgetter  import classgetter 
# =============================================================================
import ostap.math.linalg2   as LA 

# =============================================================================
revct = re.compile ( r'TVector<(?P<TYPE>[^,>]+)'      )
remtx = re.compile ( r'TMatrix[^<]*<(?P<TYPE>[^,>]+)' )
# =============================================================================

# =============================================================================
## @class LinAlgT
#  Collection of operations with TMatrixT/TVectorT
class LinAlgT(LA.LinAlg) :
    """Collection of operations with TMatrixT/TVectorT
    """
    
    # ========================================================================
    ## get number of rows for TMatrix
    #  @code
    #  mtrx = ...
    #  nrows = mtrx.kRows 
    #  @endcode  
    @staticmethod
    def T_KROWS ( mtrx ) :
        """Get number of rows for TMatrix
        >>> mtrx = ...
        >>> nrows = mtrx.kRows 
        """
        assert mtrx.IsValid () , 'Matrix is not valid!'
        return mtrx.GetNrows()
    
    # ========================================================================
    ## get number of columns for TMatrix
    #  @code
    #  mtrx = ...
    #  ncols = mtrx.kCols 
    #  @endcode  
    @staticmethod
    def T_KCOLS ( mtrx ) :
        """Get number of rows for TMatrix
        >>> mtrx = ...
        >>> ncols = mtrx.kCols 
        """
        assert mtrx.IsValid() , 'Matrix is not valid!'
        return mtrx.GetNcols()

    # ========================================================================
    ## get number ofelemets in TVectorT 
    #  @code
    #  tvec = ...
    #  size = mtrx.kSize  
    #  @endcode  
    @staticmethod
    def T_KSIZE ( tvec ) :
        """Get number of elememnt in TVectorT 
        >>> tvec = ...
        >>> tvec = tvec.kSize 
        """
        assert tvec.IsValid() , 'TVector is not valid!'
        return tvec.GetNrows()

    # =========================================================================
    ## get row from the matrix
    #  @code
    #  mtrx = ...
    #  row2 = mtrx.row ( 2 ) 
    #  @endcode
    @staticmethod
    def T_ROW ( mtrx  , index ) :
        """Get row from the matrix
        >>> mtrx = ...
        >>> row2 = mtrx.row ( 2 ) 
        """

        if not 0 <= index < mtrx.kRows :
            raise IndexError ( 'Invalid row index %s' % index )

        row = mtrx.Vector( mtrx.kCols )()  
        for i in range ( mtrx.kCols ) :
            row[i] = mtrx ( index , i )
            
        return row
        
    # =========================================================================
    ## get column from the matrix
    #  @code
    #  mtrx = ...
    #  c2 = mtrx.column ( 2 ) 
    #  @endcode
    @staticmethod
    def T_COLUMN ( mtrx  , index ) :
        """Get column from the matrix
        >>> mtrx = ...
        >>> c2 = mtrx.column ( 2 ) 
        """

        if not 0 <= index < mtrx.kCols :
            raise IndexError ( 'Invalid column index %s' % index )

        col = mtrx.Vector( mtrx.kRows )()  
        for i in range ( mtrx.kCols ) :
            col[i] = mtrx ( i , index )
            
        return col

    # =========================================================================
    ##  "Symmetrized" setter for TMatrixTSym 
    @staticmethod
    def TS_SETITEM ( mtrx , index  , value ) :
        """``Symmetrized'' setter for TMatrixTSym
        """

        i , j = index
        if not 0 <= i < mtrx.GetNrows() : raise IndexError ( 'Invalid row    index %s' % i  )
        if not 0 <= j < mtrx.GetNcols() : raise IndexError ( 'Invalid column index %s' % j  )
        
        mtrx._old_setitem_ ( ( i , j ) , value )
        if i != j :
            mtrx._old_setitem_ ( ( j , i ) , value )

    # =========================================================================
    ## convert matrix/vector-like-object to TMatrixT/TVectorT
    # @code
    # obj = ...
    # res = LinAlg.T_LA ( obj )
    # @endcode  
    @staticmethod
    def T_LA ( obj ) :
        """Convert matrix/vector-like object to TMatrixT/TVectorT object
        >>> obj = ...
        >>> res = LinAlg.T_LA ( obj )
        """
        s =  obj.shape
        assert 1 <= len ( s ) <= 2 , 'Invalid shape of the object %s' % s

        ## vector 
        if 1 == len ( s ) :
            l = s[0]            
            v = Ostap.TVectorD ( l )
            for i , p in enumerate  ( obj ) : v [ i ] = p
            return v

        ##matrix 
        mget  = LinAlgT.mgetter
        nr , nc = s
        m = Ostap.TMatrixD ( nr , nc )
        for i in range ( nr ) :
            for j in range ( nc ) :
                m [ i, j ] = mget ( obj , i, j )
            
        return m 


    # =========================================================================
    ## Decorate TVector 
    @staticmethod
    def deco_vector ( t ) :
        """ Decorate SVector
        """
        
        if t in LinAlgT.decorated_vectors : return t 

        LinAlgT.decorated_vectors.add ( t )

        LinAlgT.backup ( t )

        t. __add__      = LinAlgT. ADD
        t.__radd__      = LinAlgT.RADD
        t.__iadd__      = LinAlgT.IADD
        
        t. __sub__      = LinAlgT. SUB 
        t.__rsub__      = LinAlgT.RSUB 
        t.__isub__      = LinAlgT.ISUB 

        t. __mul__      = LinAlgT. MUL        
        t.__rmul__      = LinAlgT.RMUL 
        t.__imul__      = LinAlgT.IMUL 

        t. __matmul__   = LinAlgT. MUL ## Py3
        t.__rmatmul__   = LinAlgT.RMUL ## Py3 
        t.__imatmul__   = LinAlgT.IMUL ## Py3 

        t. __div__      = LinAlgT. DIV 
        t.__idiv__      = LinAlgT.IDIV 
        t. __truediv__  = LinAlgT. DIV 
        t.__itruediv__  = LinAlgT.IDIV 

        t.dot           = LinAlgT.DOT   
        t.cross         = LinAlgT.CROSS 
        ## t.Sim           = LinAlg.SIM
        ## t.sim           = LinAlg.SIM
        

        t. __eq__       = lambda a , b : LinAlgT.EQ ( a , b )
        t.__ne__        = lambda a , b : LinAlgT.NE ( a , b ) 
        t.__neg__       = lambda s : s*(-1)

        
        t. __str__      = LinAlgT.V_STR
        t. __repr__     = LinAlgT.V_STR

        t. __len__      = lambda s : s.GetNrows ()  
        t. __contains__ = lambda s, i : 0<=i<len(s)

        t. __iter__     = LinAlgT.V_ITER      
        t. iteritems    = LinAlgT.V_ITEMS
        t.     items    = LinAlgT.V_ITEMS

        t.to_array      = LinAlgT.V_ARRAY  ## plain array.array
        
        if LinAlgT.with_numpy : 
            t.to_numpy  = LinAlgT.V_NUMPY ## numpy array
            
            
        t.shape         = property ( LinAlgT.V_SHAPE , None , None )
        t.kSize         = property ( LinAlgT.T_KSIZE , None , None ) 

        s = revct.search ( t.__name__ )
        if s :
            stype = s.group('TYPE')
            t.Scalar     = classgetter ( lambda cls : stype ) 
            t.Element    = classgetter ( lambda cls : stype ) 
            t.value_type = classgetter ( lambda cls : stype ) 
            
        return t
        
    # =========================================================================
    ## decorate TMatrixT class
    @staticmethod
    def deco_tmatrix ( m ) :
        """Decorate TMatrixT class
        """
        if m in LinAlgT.decorated_matrices : return m

        LinAlgT.backup ( m  )
        
        m.kRows         = property ( LinAlgT.T_KROWS , None , None ) 
        m.kCols         = property ( LinAlgT.T_KCOLS , None , None ) 
        m.shape         = property ( LinAlgT.M_SHAPE , None , None )
        
        m.__str__       = LinAlgT.M_STR 
        m.__repr__      = LinAlgT.M_STR 

        m. __add__      = LinAlgT. ADD
        m.__radd__      = LinAlgT.RADD
        m.__iadd__      = LinAlgT.IADD
        
        m. __sub__      = LinAlgT. SUB 
        m.__rsub__      = LinAlgT.RSUB 
        m.__isub__      = LinAlgT.ISUB 

        m. __mul__      = LinAlgT. MUL 
        m.__rmul__      = LinAlgT.RMUL 
        m.__imul__      = LinAlgT.IMUL 

        m. __matmul__   = LinAlgT. MUL ## Py3
        m.__rmatmul__   = LinAlgT.RMUL ## Py3 
        m.__imatmul__   = LinAlgT.IMUL ## Py3 

        m. __div__      = LinAlgT. DIV 
        m.__idiv__      = LinAlgT.IDIV 
        m. __truediv__  = LinAlgT. DIV 
        m.__itruediv__  = LinAlgT.IDIV 

        m. __eq__       = lambda a , b : LinAlgT.EQ ( a , b )
        m. __ne__       = lambda a , b : LinAlgT.NE ( a , b )
        m.__neg__       = lambda s : s*(-1)

        m. __iter__     = LinAlgT.M_ITER      
        m. iteritems    = LinAlgT.M_ITEMS
        m.     items    = LinAlgT.M_ITEMS

        m.row           = LinAlgT.T_ROW
        m.column        = LinAlgT.T_COLUMN
        m.rows          = LinAlgT.M_ROWS
        m.columns       = LinAlgT.M_COLUMNS
        
        m.__pow__       = LinAlgT.M_POW 
        m.sym           = LinAlgT.M_SYM
        m.asym          = LinAlgT.M_ASYM 

        if LinAlgT.with_numpy : 
            m.to_numpy  = LinAlgT.M_NUMPY ## numpy array
            
        s = remtx.search ( m.__name__ )
        if s :
            stype = s.group('TYPE')
            m.Scalar     = classgetter ( lambda cls : stype ) 
            m.Element    = classgetter ( lambda cls : stype ) 
            m.value_type = classgetter ( lambda cls : stype )

            ## add vector type 
            m.Vector = ROOT.TVectorT ( m.Element ) 
        
        return m

    
    # =========================================================================
    ## decorate TMatrixTSym class
    @staticmethod
    def deco_tsymmatrix ( m ) :
        """Decorate TMatrixTSym class
        """
        if m in LinAlgT.decorated_matrices : return m

        LinAlgT.deco_tmatrix( m ) 
        
        ## m.__str__       = LinAlgT.TS_STR
        ## m.__repr__      = LinAlgT.TS_STR

        m._old_call_    = m.__call__
        m._old_setitem_ = m.__setitem__
        
        m.__setitem__   = LinAlgT.TS_SETITEM 
        
        return m
    
Ostap.TVector          = ROOT.TVectorT    ( 'double'  )
Ostap.TVectorF         = ROOT.TVectorT    ( 'float'   )
Ostap.TVectorD         = ROOT.TVectorT    ( 'double'  )
Ostap.TMatrix          = ROOT.TMatrixT    ( 'double'  )
Ostap.TMatrixF         = ROOT.TMatrixT    ( 'float'   )
Ostap.TMatrixD         = ROOT.TMatrixT    ( 'double'  )
Ostap.TMatrixSym       = ROOT.TMatrixTSym ( 'double'  )
Ostap.TMatrixSymF      = ROOT.TMatrixTSym ( 'float'   )
Ostap.TMatrixSymD      = ROOT.TMatrixTSym ( 'double'  )
Ostap.Math.TVector     = ROOT.TVectorT    ( 'double'  )
Ostap.Math.TVectorF    = ROOT.TVectorT    ( 'float'   )
Ostap.Math.TVectorD    = ROOT.TVectorT    ( 'double'  )
Ostap.Math.TMatrix     = ROOT.TMatrixT    ( 'double'  )
Ostap.Math.TMatrixF    = ROOT.TMatrixT    ( 'float'   )
Ostap.Math.TMatrixD    = ROOT.TMatrixT    ( 'double'  )
Ostap.Math.TMatrixSym  = ROOT.TMatrixTSym ( 'double'  )
Ostap.Math.TMatrixSymF = ROOT.TMatrixTSym ( 'float'   )
Ostap.Math.TMatrixSymD = ROOT.TMatrixTSym ( 'double'  )


for t in ( 'float' , 'double' , 'long double' ) :

    m1 = ROOT.TMatrixT      ( t  )
    LinAlgT.deco_tmatrix    ( m1 ) 

    m2 = ROOT.TMatrixTSym   ( t  )
    LinAlgT.deco_tsymmatrix ( m2 ) 
    
    v  = ROOT.TVectorT      ( t  )
    LinAlgT.deco_vector     ( v  ) 
    
import atexit
atexit.register ( LinAlgT.CLEANUP ) 

# =============================================================================
_decorated_classes_ = (
    )

_decorated_classes_ = _decorated_classes_ + tuple ( LinAlgT.decorated_vectors  )
_decorated_classes_ = _decorated_classes_ + tuple ( LinAlgT.decorated_matrices )

_new_methods_ = (
    )
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    logger.info('TEST vectors: ')
    
    logger.info('Test Linaear Algebra: ')
    
    l1 = Ostap.TVector(3)
    l2 = Ostap.TVector(3)

    l1[0],l1[1],l1[2] = 0,1,2
    l2[0],l2[1],l2[2] = 3,4,5
    
    logger.info ( 'l1 , l2 : %s %s '  % ( l1 , l2  ) )
    logger.info ( 'l1 + l2 : %s    '  % ( l1 + l2  ) )
    
    logger.info ( 'l1 - l2 : %s    '  % ( l1 - l2  ) )
    logger.info ( 'l1 * l2 : %s    '  % ( l1 * l2  ) )
    logger.info ( 'l1 *  2 : %s    '  % ( l1 *  2  ) )
    logger.info ( ' 2 * l2 : %s    '  % ( 2  * l2  ) )
    logger.info ( 'l1 /  2 : %s    '  % ( l1 /  2  ) )
    
    l1 /= 2 
    logger.info ( 'l1 /= 2 : %s    '  % l1 )
    l1 *= 2 
    logger.info ( 'l1 *= 2 : %s    '  % l1 )

    logger.info('TEST matrices: ')
    
    m22 = Ostap.Math.TMatrix(2,2)
    m23 = Ostap.Math.TMatrix(2,3) 
    s22 = Ostap.Math.TMatrixSym(2)
    
    l2  = Ostap.TVector(2)
    l3  = Ostap.TVector(3)
    
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

    logger.info ( 'm22**3\n%s' % m22**3  ) 
    logger.info ( 's22**4\n%s' % s22**4  ) 

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
##                                                                      The END 
# =============================================================================
