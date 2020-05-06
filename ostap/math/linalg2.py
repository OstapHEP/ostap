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
from   __future__  import print_function
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = ""
# =============================================================================
__all__     = (
    'mgetter'     , ## get i,j-element                  from matrix-like object
    'correlation' , ## get i,j-correlation coeffiecient from matrix-like object
    )
# =============================================================================
from   builtins    import range 
import ROOT, cppyy, re 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.linalg2' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
from ostap.math.base        import isequal , iszero , std, Ostap
from ostap.core.ostap_types import num_types, integer_types
from ostap.utils.clsgetter  import classgetter 
# =============================================================================
try :
    import numpy as np
except ImportError :
    np = None

revct = re.compile ( r'SVector<(?P<TYPE>[^,>]+)' )
remtx = re.compile ( r'SMatrix<(?P<TYPE>[^,>]+)' )

# =============================================================================
## Helper method: get  i,j element from matrix-like object
#  @code
#  mtrx  = ...
#  value = matrix ( mgetter , 1 , 2 ) 
#  @endcode 
def mgetter ( mtrx , i , j ) :
    """Helper method: get  i,j element from matrix-like object
    >>> mtrx  = ...
    >>> value = mgetter ( m , 1 , 2 ) 
    """
    
    if callable ( mtrx  ) :
        
        kRows = getattr ( mtrx , 'kRows' , None )
        if kRows is None        :
            nRows = getattr ( mtrx , 'GetNrows' , None )
            if nRows and not 0 <= i < nRows() :
                raise IndexError  ( "Invalid row index %s" % i )            
        elif not 0 <= i < kRows :
            raise IndexError  ( "Invalid row index %s" % i )
        
        kCols = getattr ( mtrx , 'kCols' , None )
        if kCols is None        :
            nCols = getattr ( mtrx , 'GetNcols' , None )
            if nCols and not 0 <= i < nCols() :
                raise IndexError  ( "Invalid column index %s" % i )
        elif not 0 <= j < kCols :
            raise IndexError  ( "Invalid column index %s" % j )

        return mtrx ( i , j )
    
    if np and isinstance ( obj , np.ndarray ) and 2 == len ( obj.shape ) :
        return obj [ i , j ]
    
    try :
        return m [ i ] [ j ]
    except :
        pass

    try :
        return m [ i , j ]
    except :
        pass
    
    return TypeError("Can't get m(%d,%d) for m=%s" % ( i , j , mtrx ) )

# =============================================================================
## Get correlation element from the matrix-like object
#  @code
#  mtrx = ...
#  corr = correlation ( mtrx , 1 , 2 ) 
#  @endcode
def correlation ( mtrx , i , j ):
    """Get the correlation element from the matrix-like object
    >>> mtrx = ...
    >>> corr = correlation ( mtrx , 1 , 2 ) 
    """
    v = matrix ( mtrx , i , j )
    
    if abs ( v ) <= 1 : return v
    elif 0 < v and isequal ( v  ,  1 ) :  return  1  
    elif 0 > v and isequal ( v  , -1 ) :  return -1  
    
    return TypeError("Can't get corr(%d,%d) for m=%s" % ( i , j , mtrx ) )

# =============================================================================
## @class Method
#  access and keep the method
class Method(object) :
    """Access and keep the method
    """
    def __init__ ( self , factory ) :
        
        self.__factory = factory
        self.__methods = {}
        
    @property
    def factory ( self ) :
        return self.__factory 

    ## get the method
    def method ( self , a , b = None ) :

        if b is None :
            args  =        a   ,
            targs = type ( a ) ,  
        elif isinstance ( b , num_types ) :
            args  =        a   ,        b 
            targs = type ( a ) , 'double'
        else :
            args  =        a   ,        b 
            targs = type ( a ) , type ( b ) 

        mm = self.__methods.get ( targs , None )
        if not mm : 
            try : 
                mm = self.__factory ( *targs )
                self.__methods[ targs ] = mm 
            except TypeError :
                mm = None 
            
        return mm
    
    def __call__ ( self , a  , b = None ) :
        return self.method ( a , b )

    def __nonzero__  ( self ) : return bool ( self.__methods ) ## or bool ( self.__factory )  
    def __bool__     ( self ) : return bool ( self.__methods ) ## or bool ( self.__factory )
    def clear        ( self ) :
        while self.__methods : self.__methods.popitem () 
# =============================================================================
## @class Method2
#  access and keep two methods 
class Method2(object) :
    """Access and keep two methods
    """

    def __init__ ( self , factory1, factory2  ) :
        
        self.__method1 = Method ( factory1 ) 
        self.__method2 = Method ( factory2 ) 
        
    @property
    def method1 ( self ) :
        return self.__method1
    @property
    def method2 ( self ) :
        return self.__method2
        
    @property
    def factory1 ( self ) :
        return self.method1.factory 
    @property
    def factory2 ( self ) :
        return self.method2.factory 

    ## get the method
    def methods ( self , a , b = None ) :

        oper  = self.__method1 ( a , b )
        if not oper : return None, None    
        check = self.__method2 ( a , b  )
        return oper, check 

    def __call__ ( self , a  , b = None ) :
        return self.methods ( a , b )
    
    def __nonzero__  ( self ) : return bool ( self.__method1 ) or bool ( self.__method2 ) 
    def __bool__     ( self ) : return bool ( self.__method1 ) or bool ( self.__method2 ) 
    def clear        ( self ) :
        self.__method1.clear() 
        self.__method2.clear() 

# ==============================================================================
## dummy function
def dummy ( *a, **b ) : return NotImplemented


# ==============================================================================
## @class LinAlg
#  collection of decorators for vectors/matrices 
class LinAlg(object) :
    """Collection of decorators for vectors/matrices
    """

    with_numpy    = np

    methods_ADD   = Method2 ( Ostap.Math.Ops.Add    , Ostap.Math.Ops.CanAdd   ) 
    methods_RADD  = Method2 ( Ostap.Math.Ops.RAdd   , Ostap.Math.Ops.CanAdd   )
    methods_IADD  = Method2 ( Ostap.Math.Ops.IAdd   , Ostap.Math.Ops.CanAdd   )

    methods_SUB   = Method2 ( Ostap.Math.Ops.Sub    , Ostap.Math.Ops.CanAdd   ) 
    methods_RSUB  = Method2 ( Ostap.Math.Ops.RSub   , Ostap.Math.Ops.CanAdd   )
    methods_ISUB  = Method2 ( Ostap.Math.Ops.ISub   , Ostap.Math.Ops.CanAdd   )

    methods_MUL   = Method2 ( Ostap.Math.Ops.Mul    , Ostap.Math.Ops.CanMul   )
    methods_RMUL  = Method2 ( Ostap.Math.Ops.RMul   , Ostap.Math.Ops.CanRMul  )
    methods_IMUL  = Method2 ( Ostap.Math.Ops.IMul   , Ostap.Math.Ops.CanIMul  )

    methods_DIV   = Method2 ( Ostap.Math.Ops.Div    , Ostap.Math.Ops.CanDiv   )
    methods_IDIV  = Method2 ( Ostap.Math.Ops.IDiv   , Ostap.Math.Ops.CanIDiv  )

    
    methods_DOT   = Method2 ( Ostap.Math.Ops.Dot    , Ostap.Math.Ops.CanDot   )
    methods_CROSS = Method2 ( Ostap.Math.Ops.Cross  , Ostap.Math.Ops.CanDot   )

    methods_SIM   = Method2 ( Ostap.Math.Ops.Sim    , Ostap.Math.Ops.CanSim   )
    methods_SIMT  = Method2 ( Ostap.Math.Ops.SimT   , Ostap.Math.Ops.CanSimT  )
    
    methods_POW   = Method2 ( Ostap.Math.Ops.Pow    , Ostap.Math.Ops.CanPow   )
    methods_SYM   = Method2 ( Ostap.Math.Ops.Sym    , Ostap.Math.Ops.CanPow   )
    methods_ASYM  = Method2 ( Ostap.Math.Ops.ASym   , Ostap.Math.Ops.CanPow   )
    
    method_EIGEN  = Method  ( Ostap.Math.Ops.Eigen  )
    method_TM     = Method  ( Ostap.Math.Ops.TM     )
    method_EQ     = Method  ( Ostap.Math.Ops.Eq     )

    known_ssymmatrices = {}
    known_smatrices    = {}
    known_svectors     = {}
    
    decorated_matrices = set ()
    decorated_vectors  = set ()

    # =========================================================================
    ## Backup useful attributes
    @staticmethod 
    def backup ( klass , attributes = ( '__add__'       , '__iadd__'     , '__radd__'    ,
                                        '__sub__'       , '__isub__'     , '__rsub__'    ,
                                        '__mul__'       , '__imul__'     , '__rmul__'    ,
                                        '__div__'       , '__idiv__'     , '__pow__'     ,
                                        '__eq__'        , '__ne__'       , '__neg__'     ,
                                        '__truediv__'   , '__itruediv__' ,
                                        '__contais__'   , '__iter__'     , 
                                        '__matmul__'    , '__rmatmul__'  , '__imatmul__' ,       
                                        '__str__'       , '__repr__'     ,
                                        '__contains__'  , '__iter__'     , '__len__'     ,
                                        'iteritems'     , 'items'        ,            
                                        'row'           , 'rows'         ,
                                        'column'        , 'columns'      ,
                                        'cross'         , 'dot'          ,
                                        'sym'           , 'asym'         , 'skew'         ,
                                        'Sim'           ,   
                                        ) ) :
        """Backup useful attributes
        """
        oa = '__old_attributes__' 
        if not hasattr  ( klass , oa ) : setattr ( klass , oa , {} )



        oatts = getattr ( klass , oa ) 
        
        for a in attributes  :
            if hasattr ( klass  , a ) :
                oatts [ a ] = getattr ( klass , a )
                
    # =========================================================================
    ## restore useful attributes
    @staticmethod 
    def restore ( klass , delete = ( 'to_numpy'      , 'to_array'       ,
                                     'tmatrix'       ,
                                     'correlations'  , 'shape'          ,
                                     'eigenValues'   , 'eigen_values'   , 
                                     'eigenVectors'  , 'eigen_vectors'  , 
                                     'sim'           , 'Sim'            ,
                                     'simT'          , 'SimT'           ,
                                     'iteritems'     , 'items'          ,                                                 
                                     'row'           , 'rows'           ,                                     
                                     'column'        , 'columns'        , 
                                     'cross'         , 'dot'            , 
                                     'sym'           , 'asym'           , 'skew'  ) ) :
        """restore useful attributesc
        """
        
        oa = '__old_attributes__' 
        if not hasattr ( klass , oa  ) : return
        
        oatts = getattr ( klass , oa )
        for k in oatts :
            setattr ( klass , k , oatts[k] ) 

        delattr ( klass , oa )
            
        for a in delete :
                if hasattr ( klass , a ) :
                    setattr ( klass , a , dummy  )
                    
    # =========================================================================
    ## cleanup the LinAlg
    @staticmethod
    def CLEANUP () :
        """Cleanup LinAlg
        """

        print ('CLEANUP-START') 

        while LinAlg.decorated_matrices :
            LinAlg.restore ( LinAlg.decorated_matrices.pop() ) 
        while LinAlg.decorated_vectors  :
            LinAlg.restore ( LinAlg.decorated_vectors .pop() ) 
        
        while LinAlg.known_ssymmatrices : LinAlg.known_ssymmatrices . popitem ()
        while LinAlg.known_smatrices    : LinAlg.known_smatrices    . popitem ()
        while LinAlg.known_svectors     : LinAlg.known_svectors     . popitem ()

        LinAlg.methods_ADD   . clear ()
        LinAlg.methods_RADD  . clear ()
        LinAlg.methods_IADD  . clear ()

        LinAlg.methods_SUB   . clear ()
        LinAlg.methods_RSUB  . clear ()
        LinAlg.methods_ISUB  . clear ()        

        LinAlg.methods_MUL   . clear ()
        LinAlg.methods_RMUL  . clear ()
        LinAlg.methods_IMUL  . clear ()        

        LinAlg.methods_SIV   . clear ()
        LinAlg.methods_IDIV  . clear ()        

        LinAlg.methods_DOT   . clear ()        
        LinAlg.methods_CROSS . clear ()        

        LinAlg.methods_SIM   . clear ()        
        LinAlg.methods_SIMT  . clear ()        

        LinAlg.methods_POW   . clear ()        
        LinAlg.methods_SYM   . clear ()        
        LinAlg.methods_ASYM  . clear ()        

        LinAlg.method_EIGEN  . clear ()        
        LinAlg.method_TM     . clear ()        
        LinAlg.method_EQ     . clear ()        

        print ('CLEANUP-END') 

    mgetter = staticmethod ( mgetter  )
        
    # =========================================================================
    ## convert vector to numpy array:
    #  @code
    #  vct = ...
    #  na = vct.to_numpy() 
    #  @endcode
    @staticmethod 
    def V_NUMPY ( vct ) :
        """Convert vector to numpy array:
        >>> vct = ...
        >>> na = vct.to_numpy() 
        """
        return np.array ( vct  , dtype = 'd' ) 
    
    # =========================================================================
    ## convert matrix into numpy.matrix
    #  @code
    #  mtrx = ...
    #  m    = mtrx.to_numpy() 
    #  @endcode 
    @staticmethod 
    def M_NUMPY ( mtrx ) :
        """Convert matrix into numpy.matrix
        >>> mtrx = ...
        >>> m    = mtrx.to_numpy() 
        """
        
        npa = np.empty ( mtrx.shape , dtype = 'd' )         
        for ij, v in mtrx.items() : npa [ ij ] = v
        return npa 

    # =========================================================================
    ## Convert matrix/vector-like object to SMatrix/SVector
    # @code
    # obj = ...
    # res = LinAlg.S_LA ( obj )
    # @endcode  
    @staticmethod
    def S_LA ( obj ) :
        """Convert matrix/vector-like object to SMatrix/SVector
        >>> obj = ...
        >>> res = LinAlg.S_LA ( obj )
        """
        s =  obj.shape
        assert 1 <= len ( s ) <= 2 , 'Invalid shape of the object %s' % s

        ## vector
        if 1 == len ( s ) : 
            
            l = s[0] 
            v = LinAlg.Vector( l )()
            for i , p in enumerate  ( obj ) : v[ i ] = p
            return v 

        ## matrix
        mget  = LinAlg.mgetter 
        nr , nc = s
        m = LinAlg.Matrix ( nr , nc )()
        for i in range ( nr ) :
            for j in range ( nc ) :
                m [ i, j ] = mget ( obj , i, j )
            
        return m 

        
    # =========================================================================
    ## addtion of matrix/vector objects
    #  @code
    #  C = A + B 
    #  @endcode 
    @staticmethod
    def ADD ( a  , b ) :
        """ Addition of vector/,atrix objects:
        >>>  C = A + B
        """
        
        if isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            s1 = a.shape
            s2 = b.shape
            if s1 != s2 : return NotImplemented
            return a.to_numpy() + b
        
        
        oper , check  = LinAlg.methods_ADD ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.add ( a, b )
            return result
        
        return NotImplemented
    
    # =========================================================================
    ## Increment  of matrix/vector objects
    #  @code
    #  A += B  
    #  @endcode 
    @staticmethod
    def IADD ( a  , b ) :
        """Increment  of matrix/vector objects
        >>> A += B  
        """

        if isinstance ( b , num_types ) : b = float( b )
        
        oper , check  = LinAlg.methods_IADD ( a , b )

        if oper and check and check.ok ( a, b ) :
            result = oper.iadd ( a, b )
            return result
        
        return NotImplemented 
                
    # =========================================================================
    ## Right-addtion of matrix/vector objects
    #  @code
    #  C = B + A  
    #  @endcode 
    @staticmethod
    def RADD ( a  , b ) :
        """ Right-addition of vector/,atrix objects:
        >>> C = B + A
        """
        if isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            s1 = a.shape
            s2 = b.shape
            if s1 != s2 : return NotImplemented
            return b + a.to_numpy()
        
        oper , check  = LinAlg.methods_RADD ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.radd ( a, b )
            return result
        
        return NotImplemented 
        
    # =========================================================================
    ## Subtraction of matrix/vector objects
    #  @code
    #  C = A - B 
    #  @endcode 
    @staticmethod
    def SUB ( a  , b ) :
        """ Subtraction of vector/,atrix objects:
        >>>  C = A - B
        """
        
        if isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            s1 = a.shape
            s2 = b.shape
            if s1 != s2 : return NotImplemented
            return a.to_numpy() - b 
        
        oper , check  = LinAlg.methods_SUB ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.sub ( a, b )
            return result
        
        return NotImplemented 

    # =========================================================================
    ## Decrement  of matrix/vector objects
    #  @code
    #  A -= B  
    #  @endcode 
    @staticmethod
    def ISUB ( a  , b ) :
        """ Decrement  of matrix/vector objects
        >>> A -= B  
        """
                
        if isinstance ( b , num_types ) : b = float( b )
        
        oper , check  = LinAlg.methods_ISUB ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.isub ( a, b )
            return result
        
        return NotImplemented 
                
    # =========================================================================
    ## Right-subtraction of matrix/vector objects
    #  @code
    #  C = B - A  
    #  @endcode 
    @staticmethod
    def RSUB ( a  , b ) :
        """ Right-subtraction of vector/matrix objects:
        >>> C = B - A
        """                
        if isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            s1 = a.shape
            s2 = b.shape
            if s1 != s2 : return NotImplemented
            return b - a.to_numpy()
        
        oper , check  = LinAlg.methods_RSUB ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.rsub ( a, b )
            return result
        
        return NotImplemented 

    # =========================================================================
    ## Multiplication/scaling of matrix/vector objects
    #  @code
    #  C = A * B 
    #  @endcode 
    @staticmethod
    def MUL ( a  , b ) :
        """  Multiplication/scaling of vector/matrix objects:
        >>>  C = A * B
        """
        
        if isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            sa = a.shape            
            sb = b.shape
            print('MUL SHAPES:', sa, sb ) 
            if sa[-1] != sb[0] :
                ##  return NotImplemented
                raise NotImplementedError ( "Cannot multiply %s/%s with %s" % ( type(a), sa , sb ) )
            return np.matmul ( a.to_numpy() , b )  
        
        oper , check  = LinAlg.methods_MUL ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.mul ( a, b )
            return result
        
        return NotImplemented 

    # =========================================================================
    ## Multiplicative increment/scaling of matrix/vector objects
    #  @code
    #  A *= B  
    #  @endcode 
    @staticmethod
    def IMUL ( a  , b ) :
        """ Multiplicative increment/scaling of matrix/vector objects
        >>> A *= B  
        """
        
        if isinstance ( b , num_types ) : b = float( b )
        
        oper , check  = LinAlg.methods_IMUL ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.imul ( a, b )
            return result
        
        return NotImplemented 
                
    # =========================================================================
    ## Right-multiplication of matrix/vector objects
    #  @code
    #  C = B * A 
    #  @endcode 
    @staticmethod
    def RMUL ( a  , b ) :
        """Right-multiplication of matrix/vector objects
        >>> C = B * A 
        """
        
        if isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            sa = a.shape            
            sb = b.shape
            print('RMUL SHAPES:', sa, sb )
            if sb[-1] != sa[0] :
                ## return NotImplemented
                raise NotImplementedError ( "Cannot multiply %s/%s with %s" % ( type(a), sa , sb ) )
            return np.matmul ( b , a.to_numpy() )  
        
        oper , check  = LinAlg.methods_RMUL ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.rmul ( a, b )
            return result
        
        return NotImplemented 
        
    # =========================================================================
    ## Division/scaling of matrix/vector objects
    #  @code
    #  C = A / B 
    #  @endcode 
    @staticmethod
    def DIV ( a  , b ) :
        """ Division/scaling of matrix/vector objects        
        >>> C = A / B 
        """
        
        if isinstance ( b , num_types ) : b = float( b )
        
        oper , check  = LinAlg.methods_DIV ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.div ( a, b )
            return result
        
        return NotImplemented 

    # =========================================================================
    ## Multiplicative decrement/scaling of matrix/vector objects
    #  @code
    #  A /= B  
    #  @endcode 
    @staticmethod
    def IDIV ( a  , b ) :
        """Multiplicative decrement/scaling of matrix/vector objects
        >>> A /= B  
        """
        
        if isinstance ( b , num_types ) : b = float( b )
        
        oper , check  = LinAlg.methods_IDIV ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.idiv ( a, b )
            return result
        
        return NotImplemented 
                
    # =========================================================================
    ## Equality for matrix/vector objects
    #  @code
    #  A == B  
    #  @endcode 
    @staticmethod
    def EQ ( a  , b ) :
        """Equality for matrix/vector objects
        >>> A == B  
        """

        if LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            s1 = a.shape            
            s2 = b.shape
            if s1 != s2 : return False 
            return np.array_equal( a.to_numpy() , b )

        oper = LinAlg.method_EQ ( a , b )
        if not oper : return NotImplemented
        
        result = oper.eq ( a, b )
        
        return result

    # =========================================================================
    ## Non-equality for matrix/vector objects
    #  @code
    #  A != B  
    #  @endcode 
    @staticmethod
    def NE ( a  , b ) :
        """Non-equality for matrix/vector objects
        >>> A != B  
        """

        result = LinAlg.EQ ( a , b )
        if result is NotImplemented : return NotImplemented
            
        return not result

    # =========================================================================
    ## Dot-product (scalar) of two vectors 
    #  @code
    #  scalar = v1.dot ( v2 ) 
    #  @endcode 
    @staticmethod
    def DOT ( a  , b ) :
        """ Dot-product (scalar) of two vectors 
        >>> scalar = v1.dot ( v2 ) 
        """
        
        if LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            s1 = a.shape            
            s2 = b.shape
            if s1 != s2  or 1 != len ( s1 ) :
                raise NotImplementedError ( "No DOT for %s/%s and %s/%s" % ( a , type ( a ) , b , type ( b ) ) )
            return np.dot ( a.to_numpy() , b )
            
        
        if isinstance ( b , num_types ) : b = float( b )
        
        oper , check  = LinAlg.methods_DOT ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.dot ( a, b )
            return result
        
        raise NotImplementedError ( "No DOT for %s/%s and %s/%s" % ( a , type ( a ) , b , type ( b ) ) )

    # =========================================================================
    ## Cross-product (D1xD2 matrix) of two vectors 
    #  @code
    #  matrix = v1.cross ( v2 ) 
    #  @endcode 
    @staticmethod
    def CROSS ( a  , b ) :
        """Cross-product (D1xD2 matrix) of two vectors 
        >>> matrix = v1.cross ( v2 ) 
        """
        
        if LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            s1 = a.shape            
            s2 = b.shape
            if s1 != s2  or 1 != len ( s1 ) :
                raise NotImplementedError ( "No CROSS for %s/%s and %s/%s" % ( a , type ( a ) , b , type ( b ) ) )
            return np.tensordot ( a.to_numpy() , b )
            
        oper , check  = LinAlg.methods_CROSS ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.cross ( a, b )
            return result
        
        raise NotImplementedError ( "No CROSS for %s/%s and %s/%s" % ( a , type ( a ) , b , type ( b ) ) )



    # =========================================================================
    ## Similarity  operation \f$ C = B A B^T \f$ 
    #  @code
    #  A = ...
    #  B = ...
    #  C = A.Sim ( B )  
    #  C = A.sim ( B )  ## ditto
    #  @endcode 
    @staticmethod
    def SIM ( a  , b ) :
        """Similarity  operation: C = B A B^T 
        >>> A = ...
        >>> B = ...
        >>> C = A.Sim ( B )   
        >>> C = A.sim ( B ) ## ditto
        """

        if isinstance ( b , num_types ) : b = float( b )
        
        oper , check  = LinAlg.methods_SIM ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.sim ( a, b )
            return result
        
        raise NotImplementedError ( "No SIM for %s/%s and %s/%s" % ( a , type ( a ) , b , type ( b ) ) )

    # =========================================================================
    ## Similarity  operation \f$ C = B^T A B  \f$ 
    #  @code
    #  A = ...
    #  B = ...
    #  C = A.SimT ( B )  
    #  C = A.simT ( B )   ## ditto 
    #  @endcode 
    @staticmethod
    def SIMT ( a  , b ) :
        """Similarity  operation C = B^T A B 
        >>> A = ...
        >>> B = ...
        >>> C = A.SimT ( B )   
        >>> C = A.simT ( B ) ## ditto 
        """

        if isinstance ( b , num_types ) : b = float( b )
        
        oper , check  = LinAlg.methods_SIMT ( a , b )
        
        if oper and check and check.ok ( a, b ) :
            result = oper.simt ( a, b )
            return result
        
        raise NotImplementedError ( "No SIMT for %s/%s and %s/%s" % ( a , type ( a ) , b , type ( b ) ) )


    # =========================================================================
    ## power function for square matrices 
    #  @code
    #  C = A ** 6 
    #  @endcode 
    @staticmethod
    def M_POW ( a  , n ) :
        """ Power function for square matrices  
        >>>  C = A ** 6 
        """

        if isinstance ( n , integer_types ) and 0 <= n : pass
        else  :
            return NotImplemented

        oper , check  = LinAlg.methods_POW ( a  )
        
        if oper and check and check.ok ( a ) :
            result = oper.pow ( a , n )
            return result
        
        return NotImplemented 

    # =========================================================================
    ## symmetric part  frot square marix
    #  @code
    #  C = A.sym() 
    #  @endcode 
    @staticmethod
    def M_SYM ( a  ) :
        """ Symmetric part of square matrices  
        >>>  C = A.sym() 
        """
        
        oper , check  = LinAlg.methods_SYM ( a  )
        
        if oper and check and check.ok ( a ) :
            result = oper.sym ( a )
            return result
        
        raise NotImplementedError ( "Cannot symmetrise %s/%s" % ( a , type(a) ) )
    
    # =========================================================================
    ## Asymmetric part  frot square marix
    #  @code
    #  C = A.asym() 
    #  C = A.skew() 
    #  @endcode 
    @staticmethod
    def M_ASYM ( a  ) :
        """ Asymmetric part of square matrices  
        >>>  C = A.asym() 
        >>>  C = A.skew() 
        """
        
        oper , check  = LinAlg.methods_ASYM ( a  )
        
        if oper and check and check.ok ( a ) :
            result = oper.asym ( a  )
            return result
        
        raise NotImplementedError ( "Cannot anti-symmetrise %s/%s" % ( a , type(a) ) )
    
    
    # =============================================================================
    ##  get matrix shape
    #   @code
    #   mtrx = ...
    #   shape = mrx.shape
    #   @endcode
    @staticmethod
    def M_SHAPE  ( mtrx ) :
        """ Get matrix shape
        >>> mtrx         = ...
        >>> nrows, ncols = mrx.shape
        """
        return mtrx.kRows, mtrx.kCols 

    # =============================================================================
    ##  get vector shape
    #   @code
    #   shape = vct.shape
    #   @endcode
    @staticmethod
    def V_SHAPE  ( mtrx ) :
        """ Get vector shape
        >>> shape = vct.shape
        """
        return mtrx.kSize,


    # =============================================================================
    ## iterator for SVector
    #  @code
    #  vct = ...
    #  for i in vct : print i 
    #  @endcode
    @staticmethod
    def V_ITER ( vct ) :
        """Iterator for SVector
        >>> vct = ...
        >>> for i in vct : print i 
        """
        s = vct.kSize 
        for i in range ( s ) : yield vct ( i )
            
    # =============================================================================
    ## iterator for SVector
    #  @code
    #  vct = ...
    #  for i,v in vct.items     () : print i,v 
    #  for i,v in vct.iteritems () : print i,v ## ditto 
    #  @endcode
    @staticmethod
    def V_ITEMS ( vct ) :
        """Iterator for SVector
        >>> vct = ...
        >>> for i,v in vct.items    () : print i,v 
        >>> for i,v in vct.iteritems() : print i,v ## ditto
        """
        s = vct.kSize 
        for i in range ( s ) :
            yield i , vct ( i )
        
    # =============================================================================
    ## self-printout of S-vectors
    #  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
    #  @date 2009-09-12
    @staticmethod
    def V_STR ( vct , fmt = ' %g' ) :
        """Self-printout of SVectors: (...)
        """

        result = '('
        for i , v in enumerate ( vct ) :
            if 0 != i : result += ', '
            result += fmt % v 
        return result + ')'

    # =============================================================================
    ## convert vector to plain array:
    #  @code
    #  vct = ...
    #  aa  = vct.to_array() 
    #  @endcode
    @staticmethod
    def V_ARRAY ( vct ) :
        """Convert vector to plain array:
        >>> vct = ...
        >>> na  = vct.to_array() 
        """
        import array
        return array.array( 'd', vct  )


    # =============================================================================
    ## iterator for SMatrix
    #  @code
    #  matrix = ...
    #  for i in matrix : print i 
    #  @endcode
    @staticmethod 
    def M_ITER ( mtrx ) :
        """Iterator for SMatrix
        >>> matrix = ...
        >>> for i in matrix : print i 
        """
        krows = mtrx.kRows
        kcols = mtrx.kCols
        for i in range ( krows ) :
            for j in range ( kcols ) :
                yield self ( i , j )
                
    # =============================================================================
    ## iterator for SMatrix
    #  @code
    #  matrix = ...
    #  for ij,v in matrix.items     () : print ij,v 
    #  for ij,v in matrix.iteritems () : print ij,v ## ditto
    #  @endcode
    @staticmethod     
    def M_ITEMS  ( mtrx ) :
        """Iterator for SMatrix
        >>> matrix = ...
        >>> for ij,v in matrix.items    () : print ij,v
        >>> for ij,v in matrix.iteritems() : print ij,v ## ditto
        """
        krows = mtrx.kRows
        kcols = mtrx.kCols
        for i in range ( krows ) :
            for j in range ( kcols ) :
                yield   ( i , j ) , mtrx ( i , j )

    # =========================================================================
    ##  Self-printout of matrices
    #   @code  
    #   matrix = ...
    #   print matrix
    #   @endcode
    @staticmethod 
    def M_STR ( mtrx , fmt = ' %+11.4g') :
        """Self-printout of matrices
        >>> matrix = ...
        >>> print matrix 
        """
        rows = mtrx.kRows
        cols = mtrx.kCols
        line = ''
        for i in range ( rows ) :
            line += ' |'
            for j in range ( cols ) :
                line += fmt % mtrx ( i , j )
            line += ' |'
            if ( rows - 1 )  != i : line += '\n'
        return line

    # =============================================================================
    ##  Self-printout of symmetric matrices
    #   @code  
    #   matrix = ...
    #   print matrix
    #   @endcode
    @staticmethod 
    def MS_STR ( mtrx , fmt = ' %+11.4g' , width = 12 ) :
        """Self-printout of symmetric matrices
        >>> matrix = ...
        >>> print matrix 
        """
        rows = mtrx.kRows
        cols = mtrx.kCols
        line = ''
        for i in range ( rows ) :
            line += ' |'
            for j in range ( cols  ) :
                if    j < i : line += width*' '
                else        : line += fmt % mtrx ( i , j )
            line += ' |'
            if ( rows - 1 ) != i : line += '\n'
        return line

    # =========================================================================
    ## get the correlation matrix
    @staticmethod 
    def MS_CORR ( mtrx ) :
        """Get the correlation matrix
        >>> mtrx = ...
        >>> corr = mtrx.correlations()
        """
        from math import sqrt
        
        _t    = type ( mtrx )
        _c    = _t   ()
        
        rows = mtrx.kRows
        ok1  = True 
        ok2  = True 
        for i in range ( rows ) :
            
            ii  = mtrx  ( i , i )
            if 0 > ii or iszero ( ii ) :
                ok1  = False 
                nan = float('nan')
                for j in range ( i , rows ) : _c [ i , j ] = nan 
                continue
        
            sii = sqrt ( ii )
            _c[ i , i ] = 1
            
            for j in range ( i + 1 ,_rows ) :            
                jj  = mtrx ( j , j ) 
                sjj = sqrt ( jj    )
                ij  = mtrx ( i , j )
                eij = ij / ( sii * sjj )
                if  1 < abs ( eij ) : ok2 = False  
                _c [ i , j ] = eij 
            
        if not ok1 : logger.error ( "correlations: zero or negative diagonal element" ) 
        if not ok2 : logger.error ( "correlations: invalid non-diagonal element"      ) 
                
        return _c


    # =========================================================================
    ## get row from the matrix
    #  @code
    #  mtrx = ...
    #  row2 = mtrx.row ( 2 ) 
    #  @endcode
    @staticmethod
    def M_ROW ( mtrx  , index ) :
        """Get row from the matrix
        >>> mtrx = ...
        >>> row2 = mtrx.row ( 2 ) 
        """

        if not 0 <= index < mtrx.kRows :
            raise IndexError ( 'Invalid row index %s' % index )

        row = mtrx.Row ( index ) 
        tr  = type ( row )
        if not tr in LinAlg.decorated_vectors : LinAlg.deco_vector ( tr )
        
        return row
        
    # =========================================================================
    ## get column from the matrix
    #  @code
    #  mtrx = ...
    #  c2 = mtrx.column ( 2 ) 
    #  @endcode
    @staticmethod
    def M_COLUMN ( mtrx  , index ) :
        """Get column from the matrix
        >>> mtrx = ...
        >>> c2 = mtrx.column ( 2 ) 
        """

        if not 0 <= index < mtrx.kCols :
            raise IndexError ( 'Invalid column index %s' % index )

        col = mtrx.Col ( index ) 
        tc  = type ( col )
        if not tc in LinAlg.decorated_vectors : LinAlg.deco_vector ( tc )
        
        return col  

    # =========================================================================
    ## Iterator over rows
    #  @code
    #  for row in mtrx.rows () : ...
    #  @endcode
    @staticmethod
    def M_ROWS ( mtrx ) :
        """Iterator over rows
        >>> for row in mtrx.rows () : ...
        """
        krows = mtrx.kRows
        for i in range  ( krows ) :
            yield mtrx.row ( i ) 
            
    # =========================================================================
    ## Iterator over columns 
    #  @code
    #  for col in mtrx.columns () : ...
    #  @endcode
    @staticmethod
    def M_COLUMNS ( mtrx ) :
        """Iterator over columns
        >>> for col in mtrx.columns () : ...
        """
        kcols = mtrx.kCols
        for i in range  ( kcols ) :
            yield mtrx.column ( i ) 
            

    # =========================================================================
    ## get the eigenvalues for symmetric matrix
    #  @code
    #  matrix = ...
    #  values = matrix.eigenValues( sorted  = True ) 
    #  @endcode
    @staticmethod
    def MS_EIGENVALUES ( mtrx , sorted = True ) :
        """Get the eigenvalues for symmetric matrices :
        >>> mtrx = ...
        >>> values = mtrx.eigenValues ( sorted = True )
        """
        
        eigen = LinAlg.method_EIGEN ( mtrx  )

        if not eigen :
            raise NotImplementedError ("EigenValues: not implemented for %s" % type ( mtrx ) )
        
        vct    = LinAlg.Vector ( mtrx.kRows )()

        st = eigen.values ( mtrx  , vct , sorted )
        assert st.isSuccess () , "Eigen-values: status code %s" % st

        return vct 

    # =========================================================================
    ## get the eigenvalues and eigen vectors for symmetric matrix
    #  @code
    #  matrix = ...
    #  values, vectors = matrix.eigenVectors( sorted  = True ) 
    #  vectors = [ vectors.column{i) for i in range ( mtrx.rCols ) ] 
    #  @endcode
    @staticmethod
    def MS_EIGENVECTORS ( mtrx , sorted = True ) :
        """Get the eigenvalues for symmetric matrices :
        >>> mtrx = ...
        >>> values, vectors = mtrx.eigenVectors ( sorted = True )
        >>> vectors = [ vectors.column{i) for i in range ( mtrx.rCols ) ] 
        """
        
        eigen = LinAlg.method_EIGEN ( mtrx )
        if not eigen :
            raise NotImplementedError ("EigenVectors: not implemented for %s" % type ( mtrx ) )

        krows = mtrx.kRows
        kcols = mtrx.kCols
        
        values  = LinAlg.Vector ( krows        ) ()
        vectors = LinAlg.Matrix ( krows, kcols ) () 
        
        st = eigen.vectors ( mtrx  , values , vectors , sorted )
        assert st.isSuccess () , "Eigen-vectors: status code %s" % st

        return values, vectors  


    # =========================================================================
    ## Convert matrix to TMatrix, vector to TVector 
    #  @code
    #  matrix = ...
    #  tm = matrix.tmatrix() 
    #  @endcode
    #  @code
    #  vector = ...
    #  tv     =  vector.tvector () 
    #  @endcode
    @staticmethod
    def M_TM ( sobj ) :
        """Convert matrix to TMatrix, vector to TVector 
        >>> matrix = ...
        >>> tm     = matrix.tmatrix() 
        >>> vector = ...
        >>> tv =   = vector.tvector() 
        """
        
        cnv = LinAlg.method_TM ( sobj )
        if not cnv :
            raise NotImplementedError ("SMatrix->TMatrix/SVector->TVector: not implemented for %s" % type ( sobj ) )

        tobj = cnv.transform ( sobj ) 
        assert tobj.IsValid() , "Smatrix->TMatrix/SVector->TVector: invaild TMatrix/TVector!" 

        return tobj

    # =========================================================================
    ## Decorate SVector 
    @staticmethod
    def deco_vector ( t ) :
        """ Decorate SVector
        """
        
        if t in LinAlg.decorated_vectors : return t 

        LinAlg.decorated_vectors.add ( t )

        LinAlg.backup ( t )

        t. __add__      = LinAlg. ADD
        t.__radd__      = LinAlg.RADD
        t.__iadd__      = LinAlg.IADD
        
        t. __sub__      = LinAlg. SUB 
        t.__rsub__      = LinAlg.RSUB 
        t.__isub__      = LinAlg.ISUB 

        t. __mul__      = LinAlg. MUL 
        t.__rmul__      = LinAlg.RMUL 
        t.__imul__      = LinAlg.IMUL 

        t. __matmul__   = LinAlg. MUL ## Py3
        t.__rmatmul__   = LinAlg.RMUL ## Py3 
        t.__imatmul__   = LinAlg.IMUL ## Py3 

        t. __div__      = LinAlg. DIV 
        t.__idiv__      = LinAlg.IDIV 
        t. __truediv__  = LinAlg. DIV 
        t.__itruediv__  = LinAlg.IDIV 

        t.dot           = LinAlg.DOT   
        t.cross         = LinAlg.CROSS 
        t.Sim           = LinAlg.SIM
        t.sim           = LinAlg.SIM
        

        t. __eq__       = lambda a , b : LinAlg.EQ ( a , b )
        t. __ne__       = lambda a , b : LinAlg.NE ( a , b )
        
        t.__neg__       = lambda s : s*(-1)

        
        t. __str__      = LinAlg.V_STR
        t. __repr__     = LinAlg.V_STR

        t. __len__      = lambda s : s.kSize 
        t. __contains__ = lambda s, i : 0<=i<s.kSize

        t. __iter__     = LinAlg.V_ITER      
        t. iteritems    = LinAlg.V_ITEMS
        t.     items    = LinAlg.V_ITEMS

        t.to_array      = LinAlg.V_ARRAY ## plain array.array 

        t.tvector       = LinAlg.M_TM 

        t.shape         = property ( LinAlg.V_SHAPE , None , None )

        s = revct.search ( t.__name__ )
        if s :
            stype = s.group('TYPE')
            t.Scalar     = classgetter ( lambda cls : stype ) 
            t.Element    = classgetter ( lambda cls : stype ) 
            t.value_type = classgetter ( lambda cls : stype ) 
            
        if LinAlg.with_numpy : 
            t.to_numpy  = LinAlg.V_NUMPY ## numpy array 
            
        return t

    # =========================================================================
    ## Decorate SMatrix 
    @staticmethod 
    def deco_matrix ( m  ) :
        """Decorate SMatrix
        """
        
        if m in LinAlg.decorated_matrices : return m 
        
        LinAlg.decorated_matrices.add ( m )

        ##  save 'old method'
        LinAlg.backup ( m )
        
        m. __add__      = lambda a , b : LinAlg. ADD ( a , b ) 
        m.__radd__      = lambda a , b : LinAlg.RADD ( a , b ) 
        m.__iadd__      = lambda a , b : LinAlg.IADD ( a , b )
        
        m. __sub__      = LinAlg. SUB 
        m.__rsub__      = LinAlg.RSUB 
        m.__isub__      = LinAlg.ISUB 

        m. __mul__      = LinAlg. MUL 
        m.__rmul__      = LinAlg.RMUL 
        m.__imul__      = LinAlg.IMUL 

        m. __matmul__   = LinAlg. MUL ## Py3
        m.__rmatmul__   = LinAlg.RMUL ## Py3 
        m.__imatmul__   = LinAlg.IMUL ## Py3 

        m. __div__      = LinAlg. DIV 
        m.__idiv__      = LinAlg.IDIV 
        m. __truediv__  = LinAlg. DIV 
        m.__itruediv__  = LinAlg.IDIV 
        m. __eq__       = lambda a , b : LinAlg.EQ ( a , b )
        m. __ne__       = lambda a , b : LinAlg.NE ( a , b )
        m.__neg__       = lambda s   : s*(-1)

        m.__str__       = LinAlg.M_STR
        m.__repr__      = LinAlg.M_STR

        m.__iter__      = LinAlg.M_ITER 
        m.iteritems     = LinAlg.M_ITEMS 
        m.    items     = LinAlg.M_ITEMS
        m.__contains__  = lambda s,ij : 0<=ij[0]<s.kRows and 0<=ij[1]<s.kCols

        m.row           = LinAlg.M_ROW
        m.column        = LinAlg.M_COLUMN
        m.rows          = LinAlg.M_ROWS
        m.columns       = LinAlg.M_COLUMNS

        m.tmatrix       = LinAlg.M_TM 

        m.shape         = property ( LinAlg.M_SHAPE , None , None )

        m.__pow__       = LinAlg.M_POW 
        m.sym           = LinAlg.M_SYM
        m.asym          = LinAlg.M_ASYM 
        m.skew          = LinAlg.M_ASYM 

        s = remtx.search ( m.__name__ )
        if s :
            stype = s.group('TYPE')
            m.Scalar     = classgetter ( lambda cls : stype ) 
            m.Element    = classgetter ( lambda cls : stype ) 
            m.value_type = classgetter ( lambda cls : stype ) 

        if LinAlg.with_numpy  : 
            m.to_numpy  = LinAlg.M_NUMPY ## numpy array 

        return m
    
    # =========================================================================
    ## Decorate symmetric SMatrix 
    @staticmethod 
    def deco_symmatrix ( m ) :
        """Decorate the symmetrix  SMatrix
        """

        if m in LinAlg.decorated_matrices : return m 

        LinAlg.deco_matrix ( m )
        
        m.__str__        = LinAlg.MS_STR
        m.__repr__       = LinAlg.MS_STR
        m.correlations   = LinAlg.MS_CORR
        
        m.Sim            = LinAlg.SIM
        m.sim            = LinAlg.SIM
        m.SimT           = LinAlg.SIMT
        m.simT           = LinAlg.SIMT
        
        m.eigenValues    = LinAlg.MS_EIGENVALUES 
        m.eigenVectors   = LinAlg.MS_EIGENVECTORS
        m.eigen_values   = LinAlg.MS_EIGENVALUES 
        m.eigen_vectors  = LinAlg.MS_EIGENVECTORS

        return m

    # =========================================================================
    ##  Pick up the vector of corresponding size
    #   @code
    #   V3   = Ostap.Math.Vector(3)
    #   vct  = V3 ()
    #   @endcode  
    @staticmethod
    def Vector ( n , t = 'double' ) :
        """Pick up the vector of corresponding size
        >>> V3   = Ostap.Math.Vector(3)
        >>> vct  = V3 ()
        """
        assert isinstance  ( n , integer_types ) and 0 <= n,\
               'Invalid length of the vector %s/%s' % ( n , type ( l ) )

        tt = n , t 
        v = LinAlg.known_svectors.get ( tt , None )
        if  v is None :
            v = ROOT.ROOT.Math.SVector ( t , n )
            LinAlg.known_svectors [ tt ] = v
            LinAlg.deco_vector    ( v )  
        ##
        return v 
    
    # =========================================================================
    ## Pick up the matrix of corresponding size
    #  @code 
    #  M3x4   = Ostap.Math.Matrix(3,4)
    #  matrix = M3x4 ()    
    #  @endcode 
    @staticmethod
    def Matrix ( k , n  , t = 'double' ) :
        """Pick up the matrix of corresponding size
        >>> M3x4   = Ostap.Math.Matrix(3,4)
        >>> matrix = M3x4 ()    
        """
        assert isinstance  ( k , integer_types ) and 0 <= k ,\
               'Invalid matrix dimension %s/%s' % ( k , type ( k ) )
        assert isinstance  ( n , integer_types ) and 0 <= n ,\
               'Invalid matrix dimension %s/%s' % ( n , type ( n ) )

        tt = k , n , t
        m = LinAlg.known_smatrices.get ( tt , None )
        if m is None :
            m = ROOT.ROOT.Math.SMatrix ( t , k , n )
            LinAlg.known_smatrices [ tt ] = m
            LinAlg.deco_matrix ( m )
            
        return m 


    # =========================================================================
    ## Pick up the symmetric matrix of corresponding size
    #  @code
    #  SymM3  = Ostap.Math.SymMatrix(3)
    #  matrix = SymM3 ()
    #  @endcode  
    @staticmethod
    def SymMatrix  ( n , t  = 'double' ) :
        """Pick up the symmetric matrix of corresponding size
        >>> SymM3  = Ostap.Math.SymMatrix(3)
        >>> matrix = SymM3 ()
        """
        assert isinstance  ( n , integer_types ) and 0 <= n ,\
               'Invalid matrix dimension %s/%s' % ( n , type ( n ) )

        tt = n , t
        m = LinAlg.known_ssymmatrices.get ( tt , None )
        if  m is None  : 
            m = ROOT.ROOT.Math.SMatrix('%s,%d,%d,ROOT::Math::MatRepSym<%s,%d>' %  ( t , n , n , t , n ) )
            LinAlg.known_ssymmatrices [ tt ] = m
            LinAlg.deco_symmatrix  ( m )

        return m
    
    # =========================================================================
    ##  is this matrix/vector type decorated properly?
    @staticmethod
    def decorated ( t ) :
        """Is this matrix/vector type decorated properly?
        """        
        return t in LinAlg.decorated_vectors or t in LinAlg.decorated_matrices
    

# =======================================================================

Ostap.Vector         =  staticmethod ( LinAlg.Vector    )
Ostap.Matrix         =  staticmethod ( LinAlg.Matrix    )
Ostap.SymMatrix      =  staticmethod ( LinAlg.SymMatrix ) 
Ostap.Math.Vector    =  staticmethod ( LinAlg.Vector    )
Ostap.Math.Matrix    =  staticmethod ( LinAlg.Matrix    )
Ostap.Math.SymMatrix =  staticmethod ( LinAlg.SymMatrix ) 



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
    t1 = Ostap.SymMatrix(i)
    for j in range(11) : t2 = Ostap.Matrix(i,j)

# =============================================================================
_decorated_classes_ = (
    Ostap.Vector(2)     ,
    Ostap.Matrix(2,3)   , 
    Ostap.SymMatrix(2)  , 
    )

_decorated_classes_ = _decorated_classes_ + tuple ( LinAlg.decorated_vectors  )
_decorated_classes_ = _decorated_classes_ + tuple ( LinAlg.decorated_matrices )

_new_methods_ = (
    Ostap.Vector    , 
    Ostap.Matrix    , 
    Ostap.SymMatrix , 
    )


import atexit
atexit.register ( LinAlg.CLEANUP ) 


# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    logger.info('Test Linaear Algebra: ')
    
    LA3 = Ostap.Vector3
    l1  = LA3(0,1,2)
    l2  = LA3(3,4,5)
    
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
    
    logger.info ( 'm22\n%s'    % m22     ) 
    logger.info ( 's22\n%s'    % s22     ) 
    logger.info ( 'm23\n%s'    % m23     ) 
    logger.info ( 'm22/3\n%s'  % (m22/3) ) 
    logger.info ( 'm23*3\n%s'  % (m23*3) ) 

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
