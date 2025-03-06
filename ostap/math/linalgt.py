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
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = ""
# =============================================================================
__all__     = ( )
# =============================================================================
from   sys                    import version_info as python_version
from   ostap.math.base        import isequal , iszero , std, Ostap
from   ostap.core.ostap_types import num_types, integer_types
from   ostap.utils.clsgetter  import classgetter
from   ostap.utils.gsl        import gsl_info  
import ostap.math.linalg2     as     LA 
import ROOT, re 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.linalgt' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
revct = re.compile ( r'TVector<(?P<TYPE>[^,>]+)'      )
remtx = re.compile ( r'TMatrix[^<]*<(?P<TYPE>[^,>]+)' )
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
    ## get number of elements in TVectorT 
    #  @code
    #  tvec = ...
    #  size = mtrx.kSize  
    #  @endcode  
    @staticmethod
    def T_KSIZE ( tvec ) :
        """Get number of elements in TVectorT 
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

    # ============================================================================
    # transpose the matrix 
    @staticmethod
    def T_T ( mtrx ) :
        """ Transpose the matrix (make a copy)
        """
        newmatrix = Ostap.Math.TMatrix ( mtrx )
        newmatrix.T() 
        return newmatrix 

    # ============================================================================
    # transpose the matrix 
    @staticmethod
    def TS_T ( mtrx ) :
        """ Transpose the matrix (make a copy) 
        """
        return Ostap.Math.TMatrixSym ( mtrx )
    
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
    ## convert TVector to SVector
    @staticmethod
    def TV_SV ( obj ) :
        """Convert TVector to SVector"""
        N      = len ( obj )
        result = Ostap.SVector( N )()
        for i in range ( N ) :
            result [ i ] = obj[i]
        return result 

    # =========================================================================
    ## convert TMatrix to SMatrix
    @staticmethod
    def TM_SM ( obj ) :
        """Convert TMatrix to SMatrix"""
        nr, nc = obj.shape        
        result = Ostap.SMatrix( nr , nc )()
        for i in range ( nr ) :
            for j in range ( nc ) : 
                result [ i,j ] = obj[i,j]                
        return result 

    # =========================================================================
    ## convert TMatrixTSym to SMatrixSym
    @staticmethod
    def TMS_SMS ( obj ) :
        """Convert TMatrixTSym to SMatrixSym"""
        nr, nc = obj.shape        
        result = Ostap.SymMatrix ( nr )()
        for i in range ( nr ) :
            for j in range ( i , nr ) :                
                result [ i,j ] = 0.5 * ( obj (i,j) + obj(j,i) )                
        return result 


    # =========================================================================
    ## (P)LU decomposition
    #  @see Ostap::GSL::PLU
    #  @code
    #  matrix = ...
    #  P , L , U = matrix.PLU() 
    #  @endcode    
    @staticmethod
    def T_PLU  ( mtrx ) :
        """ Perform (P)LU decomposition of the matrix 
        >>> matrix = ...\
        >>> P, L, U = matarix.PLU() 
        - see `Ostap.GSL.PLU` 
        """
        assert mtrx.IsValid () , 'Matrix is not valid!'
        ## convert to GLS 
        A = mtrx.to_GSL()
        ## mape (P)LU decomposiiton 
        P, L, U = A.PLU ()
        ## convert BACK:
        P = Ostap.GSL.Matrix ( P )
        ## 
        P = P.to_TMatrix()
        L = L.to_TMatrix()
        U = U.to_TMatrix()
        ## 
        return P , L , U 
        
    # =========================================================================
    ## (P)QR decomposition
    #  @see Ostap::GSL::PQR
    #  @code
    #  matrix = ...
    #  P , Q , R = matrix.PQR() 
    #  @endcode    
    @staticmethod
    def T_PQR  ( mtrx ) :
        """ Perform (P)QR decomposition of the matrix 
        >>> matrix = ...\
        >>> P, Q, R = matarix.PQR() 
        - see `Ostap.GSL.PQR` 
        """ 
        assert mtrx.IsValid () , 'Matrix is not valid!'
        ## convert to GLS 
        A = mtrx.to_GSL()
        ## mape (P)LU decomposiiton 
        P, Q, R = A.PQR ()
        ## convert BACK:
        P = Ostap.GSL.Matrix ( P ).T() 
        ## 
        P = P.to_TMatrix()
        Q = Q.to_TMatrix()
        R = R.to_TMatrix()
        ## 
        return P , Q , R 

    # ==========================================================================
    ## Perform LQ decomposition of the matrix
    @staticmethod
    def T_LQ ( mtrx ) :
        """ Perfrom LQ decompositionof matrix A : A = LQ 
        - A is input MxN matrix 
        - L is lower trapezoidal  MxN matrix
        - Q is orthogonal NxN matrix    
        >>> A = ...
        >>> L, Q = A.LQ() 
        """
        assert mtrx.IsValid () , 'Matrix is not valid!'
        ## convert matrix to GSL 
        A = mtrx.to_GSL()
        ## make decomposition
        L , Q = A.LQ()
        ## convert back
        L = L.to_SMatrix()
        Q = Q.to_SMatrix()
        return L, Q

    ## Perform QL decomposition of the matrix
    @staticmethod
    def T_QL ( mtrx ) :
        """ Perfrom LQ decompositionof matrix A : A = QL 
        - A is input MxN matrix 
        - Q is orthogonal NxN matrix    
        - L is lower trapezoidal  MxN matrix
        >>> A = ...
        >>> Q , L = A.QL () 
        """
        assert mtrx.IsValid () , 'Matrix is not valid!'
        ## convert matrix to GSL 
        A = mtrx.to_GSL()
        ## make decomposition
        Q , L = A.QL ()
        ## convert back
        Q = Q.to_TMatrix()
        L = L.to_TMatrix()
        return Q , L 

    # =========================================================================
    ## COD: Complete orthogonal decomposition AP = Q R Z^T
    @staticmethod
    def T_COD ( mtrx ) :
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
        assert mtrx.IsValid () , 'Matrix is not valid!'
        ## convert matrix to GSL 
        A = mtrx.to_GSL()
        ## make decomposition
        P , Q , R , Z = A.COD()
        ## convert BACK:
        P = Ostap.GSL.Matrix ( P ).T() 
        ## 
        P = P.to_TMatrix()
        Q = Q.to_TMatrix()
        R = R.to_TMatrix()
        Z = Z.to_TMatrix()
        return P, Q, R , Z 
    
    # ===============================================================================
    ## SVD : singular Value Decomposition  \f$ A = U S V^T\f$
    def T_SVD ( mtrx , golub = True ) : 
        """ SVD : singular Value Decomposition  A = U S V^T 
        - A input MxN matrix 
        - K = min ( M , N ) : 
        - U MxK orthogonal matrix 
        - S KxK Diagonal matrix of singular values 
        - V NxK orthogonal matrix 
        >>> A = ...
        >>> S , U , V = A.SVD() 
        """
        assert mtrx.IsValid () , 'Matrix is not valid!'
        ## convert matrix to GSL 
        A = mtrx.to_GSL()
        ## make decomposition
        S , U , V = A.SVD() 
        ## convert BACK:
        N = S.size() 
        M = Ostap.TMatrixSym ( N )
        for i in range ( N ) : M [ i , i ] = S ( i )
        S = M 
        U = U.to_TMatrix()
        V = V.to_TMatrix()
        return S , U , V
    
    # ===============================================================================
    ## SCHUR : Schur Decomposition  \f$ A = Z T Z^T \f$ 
    def T_SCHUR ( mtrx ) :
        """ Schur decomposition of the square matrix A: A = Z T Z^T
        - Z is orthogonal 
        - T is Schur form 
        >>> A = ...,
        >>> Z , T = A.SCHUR() 
        """
        assert mtrx.IsValid () , 'Matrix is not valid!'
        assert mtrx.GetNrows() == mtrx.GetNcols() , \
            "Schur decomposition is defined only for square matrices!"
        ## convert matrix to GSL 
        A = mtrx.to_GSL()
        ## make decomposition
        Z , T = A.SCHUR () 
        ## convert BACK:
        Z = Z.to_TMatrix()
        T = T.to_TMatrix()
        return Z, T 

    # ===============================================================================
    ## POLAR : Polar Decomposition  \f$ A = UP \f$ 
    def T_POLAR ( mtrx ) :
        """ Polar decomposition of the square matrix A: A = UP
        - U is orthogonal 
        - P is positive semi-definitive 
        >>> A = ...,
        >>> U , P = A.POLAR() 
        """
        assert mtrx.IsValid () , 'Matrix is not valid!'
        assert mtrx.GetNrows() == mtrx.GetNcols() , \
            "Polar decomposition is defined only for square matrices!"
        ## convert matrix to GSL 
        A = mtrx.to_GSL()
        ## make decomposition
        U , P = A.POLAR() 
        ## convert BACK:
        P = P.to_TMatrix()
        U = U.to_TMatrix()
        return U , P 
    
    # ==========================================================================
    
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

        if ( 3 , 5 ) <= python_version : 
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
        
        t.__abs__       = LinAlgT.V_ABS

        t. __str__      = LinAlgT.V_STR
        t. __repr__     = LinAlgT.V_STR
        t. table        = LinAlgT.V_STR        
        t.pretty_print  = LinAlgT.V_PRETTY

        t. __len__      = lambda s : s.GetNrows ()  
        t. __contains__ = lambda s, i : 0<=i<len(s)

        t. __iter__     = LinAlgT.V_ITER      
        t. iteritems    = LinAlgT.V_ITEMS
        t.     items    = LinAlgT.V_ITEMS

        t.to_array      = LinAlgT.V_ARRAY  ## plain array.array

        if LinAlgT.with_numpy : 
            t.to_numpy  = LinAlgT.V_NUMPY ## numpy array

        ## convert TVector to SVector 
        t.svector       = LinAlgT.TV_SV 
        
        t.to_gsl        = LinAlgT.V_2GSL 
        t.as_gsl        = LinAlgT.V_2GSL 
        t.to_GSL        = LinAlgT.V_2GSL 
        t.as_GSL        = LinAlgT.V_2GSL 
        
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
        """ Decorate TMatrixT class
        """
        if m in LinAlgT.decorated_matrices : return m

        LinAlgT.backup ( m  )
        
        m.kRows         = property ( LinAlgT.T_KROWS , None , None ) 
        m.kCols         = property ( LinAlgT.T_KCOLS , None , None ) 
        m.shape         = property ( LinAlgT.M_SHAPE , None , None )
        
        m.__str__       = LinAlgT.M_STR 
        m.__repr__      = LinAlgT.M_STR 
        m.table         = LinAlgT.M_STR
        m.pretty_print  = LinAlgT.M_PRETTY

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
        
        m.t             = LinAlgT.T_T   
        m.transpose     = LinAlgT.T_T   

        m.__pow__       = LinAlgT.M_POW 
        m.sym           = LinAlgT.M_SYM
        m.asym          = LinAlgT.M_ASYM

        m.diagonal      = LinAlgT.M_DIAGONAL
        m.lnorm         = LinAlgT.M_LNORM 
        m.mnorm         = LinAlgT.M_MNORM 

        m.inverse       = LinAlgT.M_INVERSE            

        ## convert TMatrix to SMatrix 
        m.smatrix       = LinAlgT.TM_SM

        if LinAlgT.with_numpy : 
            m.to_numpy  = LinAlgT.M_NUMPY ## numpy array

        m.to_gsl        = LinAlgT.M_2GSL 
        m.as_gsl        = LinAlgT.M_2GSL 
        m.to_GSL        = LinAlgT.M_2GSL 
        m.as_GSL        = LinAlgT.M_2GSL 

        m.PLU           = LinAlgT.T_PLU 
        m.PQR           = LinAlgT.T_PQR
        m.LQ            = LinAlgT.T_LQ    
        m.COD           = LinAlgT.T_COD
        m.SVD           = LinAlgT.T_SVD
        
        if  ( 2  , 7 ) <= gsl_info :
            m.QL        = LinAlgT.T_QL
            
        m.SCHUR          = LinAlgT.T_SCHUR
        m.POLAR          = LinAlgT.T_POLAR
        
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
        """ Decorate TMatrixTSym class
        """
        if m in LinAlgT.decorated_matrices : return m

        LinAlgT.deco_tmatrix( m ) 
        
        ## convert TMatrixSym to SymMatrix 
        m.smatrix       = LinAlgT.TMS_SMS

        m.__str__       = LinAlgT.MS_STR 
        m.__repr__      = LinAlgT.MS_STR 
        m.table         = LinAlgT.MS_STR
        m.pretty_print  = LinAlgT.MS_PRETTY 

        m.t             = LinAlgT.TS_T   
        m.transpose     = LinAlgT.TS_T   


        m.Sim            = LinAlgT.SIM
        m.sim            = LinAlgT.SIM
        m.SimT           = LinAlgT.SIMT
        m.simT           = LinAlgT.SIMT

        m._old_call_     = m.__call__
        m._old_setitem_  = m.__setitem__
        
        m.__setitem__    = LinAlgT.TS_SETITEM 
        
        return m
    
Ostap.TVector          = ROOT.TVectorT    ( 'double'  )
Ostap.TVectorD         = ROOT.TVectorT    ( 'double'  )
Ostap.TMatrix          = ROOT.TMatrixT    ( 'double'  )
Ostap.TMatrixD         = ROOT.TMatrixT    ( 'double'  )
Ostap.TMatrixSym       = ROOT.TMatrixTSym ( 'double'  )
Ostap.TMatrixSymD      = ROOT.TMatrixTSym ( 'double'  )
Ostap.Math.TVector     = ROOT.TVectorT    ( 'double'  )
Ostap.Math.TVectorD    = ROOT.TVectorT    ( 'double'  )
Ostap.Math.TMatrix     = ROOT.TMatrixT    ( 'double'  )
Ostap.Math.TMatrixD    = ROOT.TMatrixT    ( 'double'  )
Ostap.Math.TMatrixSym  = ROOT.TMatrixTSym ( 'double'  )
Ostap.Math.TMatrixSymD = ROOT.TMatrixTSym ( 'double'  )


##  for t in ( 'float' , 'double' , 'long double' ) :
for t in ( 'double' , ) :

    m1 = ROOT.TMatrixT      ( t  )
    LinAlgT.deco_tmatrix    ( m1 ) 

    m2 = ROOT.TMatrixTSym   ( t  )
    LinAlgT.deco_tsymmatrix ( m2 ) 
    
    v  = ROOT.TVectorT      ( t  )
    LinAlgT.deco_vector     ( v  ) 
    

import atexit
atexit.register ( LinAlgT.CLEANUP ) 

# =============================================================================
_decorated_classes_ = ()

_decorated_classes_ = _decorated_classes_ + tuple ( LinAlgT.decorated_vectors  )
_decorated_classes_ = _decorated_classes_ + tuple ( LinAlgT.decorated_matrices )

_new_methods_ = ()
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    

# =============================================================================
##                                                                      The END 
# =============================================================================
