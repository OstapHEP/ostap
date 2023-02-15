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
    'checkops'    , ## check the allowed operations  
    'correlation' , ## get i,j-correlation coeffiecient from matrix-like object
    )
# =============================================================================
import ROOT, re, ctypes, array  
from   sys       import version_info as python_version
from   builtins  import range
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.linalg2' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
from ostap.math.base        import isequal   , iszero, std , Ostap, typename 
from ostap.core.ostap_types import num_types , integer_types
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
            except TypeError :
                return None
            
            if mm :
                self.__methods[ targs ] = mm 

        return mm.operation if mm else None 

    
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

    def __init__ ( self , factory1 , factory2  ) :

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

        operation = self.__method1 ( a , b )        
        checker   = self.__method2 ( a , b ) if operation else None 
        
        return operation, checker  

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

    methods_ADD   = Method2 ( Ostap.Math.Ops.Add    , Ostap.Math.Ops.CanAdd     ) 
    methods_RADD  = Method2 ( Ostap.Math.Ops.RAdd   , Ostap.Math.Ops.CanAdd     )
    methods_IADD  = Method2 ( Ostap.Math.Ops.IAdd   , Ostap.Math.Ops.CanAdd     )

    methods_SUB   = Method2 ( Ostap.Math.Ops.Sub    , Ostap.Math.Ops.CanAdd     )   
    methods_RSUB  = Method2 ( Ostap.Math.Ops.RSub   , Ostap.Math.Ops.CanAdd     )
    methods_ISUB  = Method2 ( Ostap.Math.Ops.ISub   , Ostap.Math.Ops.CanAdd     )

    methods_MUL   = Method2 ( Ostap.Math.Ops.Mul    , Ostap.Math.Ops.CanMul     )
    methods_RMUL  = Method2 ( Ostap.Math.Ops.RMul   , Ostap.Math.Ops.CanRMul    )
    methods_IMUL  = Method2 ( Ostap.Math.Ops.IMul   , Ostap.Math.Ops.CanIMul    )

    methods_DIV   = Method2 ( Ostap.Math.Ops.Div    , Ostap.Math.Ops.CanDiv     )
    methods_IDIV  = Method2 ( Ostap.Math.Ops.IDiv   , Ostap.Math.Ops.CanIDiv    )

    
    methods_DOT   = Method2 ( Ostap.Math.Ops.Dot    , Ostap.Math.Ops.CanDot     )
    methods_CROSS = Method2 ( Ostap.Math.Ops.Cross  , Ostap.Math.Ops.CanDot     )

    methods_SIM   = Method2 ( Ostap.Math.Ops.Sim    , Ostap.Math.Ops.CanSim     )
    methods_SIMT  = Method2 ( Ostap.Math.Ops.SimT   , Ostap.Math.Ops.CanSimT    )
    
    methods_POW   = Method2 ( Ostap.Math.Ops.Pow    , Ostap.Math.Ops.CanPow     )
    methods_SYM   = Method2 ( Ostap.Math.Ops.Sym    , Ostap.Math.Ops.CanSym     )
    methods_ASYM  = Method2 ( Ostap.Math.Ops.ASym   , Ostap.Math.Ops.CanASym    )
    
    methods_EQ    = Method2 ( Ostap.Math.Ops.Eq     , Ostap.Math.Ops.CanEq      )
    methods_INV   = Method2 ( Ostap.Math.Ops.Invert , Ostap.Math.Ops.CanInvert  )
    
    method_EIGEN  = Method  ( Ostap.Math.Ops.Eigen  )
    method_TM     = Method  ( Ostap.Math.Ops.TM     )


    ## method_EQ     = Method  ( Ostap.Math.Ops.Eq     )

    
    known_ssymmatrices = {}
    known_smatrices    = {}
    known_svectors     = {}
    known_svectorse    = {}
    
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
        """restore useful attributes
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

        while LinAlg.decorated_matrices :
            LinAlg.restore ( LinAlg.decorated_matrices.pop() ) 
        while LinAlg.decorated_vectors  :
            LinAlg.restore ( LinAlg.decorated_vectors .pop() ) 
        
        
        while LinAlg.known_ssymmatrices : LinAlg.known_ssymmatrices . popitem ()
        while LinAlg.known_smatrices    : LinAlg.known_smatrices    . popitem ()
        while LinAlg.known_svectors     : LinAlg.known_svectors     . popitem ()
        while LinAlg.known_svectorse    : LinAlg.known_svectorse    . popitem ()

        if LinAlg.methods_ADD    :
            LinAlg.methods_ADD   . clear ()
            LinAlg.methods_ADD   = None 

        if LinAlg.methods_RADD   :
            LinAlg.methods_RADD  . clear ()
            LinAlg.methods_RADD  = None 
        
        if  LinAlg.methods_IADD  :
            LinAlg.methods_IADD  . clear ()
            LinAlg.methods_IADD  = None 

        if LinAlg.methods_SUB    :
            LinAlg.methods_SUB   . clear ()
            LinAlg.methods_SUB   = None 

        if LinAlg.methods_RSUB   : 
            LinAlg.methods_RSUB  . clear ()
            LinAlg.methods_RSUB  = None 

        if LinAlg.methods_ISUB   :
            LinAlg.methods_ISUB  . clear ()        
            LinAlg.methods_ISUB  = None         

        if LinAlg.methods_MUL    : 
            LinAlg.methods_MUL   . clear ()
            LinAlg.methods_MUL   = None
            
        if LinAlg.methods_RMUL   : 
            LinAlg.methods_RMUL  . clear ()
            LinAlg.methods_RMUL  = None 

        if LinAlg.methods_IMUL   : 
            LinAlg.methods_IMUL  . clear ()
            LinAlg.methods_IMUL  = None 

        if LinAlg.methods_DIV    : 
            LinAlg.methods_DIV   . clear ()
            LinAlg.methods_DIV    = None
            
        if LinAlg.methods_IDIV   : 
            LinAlg.methods_IDIV  . clear ()
            LinAlg.methods_IDIV   = None 

        if LinAlg.methods_DOT    : 
            LinAlg.methods_DOT   . clear ()
            LinAlg.methods_DOT    = None 

        if LinAlg.methods_CROSS  : 
            LinAlg.methods_CROSS . clear ()
            LinAlg.methods_CROSS  = None 

        if LinAlg.methods_SIM    : 
            LinAlg.methods_SIM   . clear ()
            LinAlg.methods_SIM   = None 

        if LinAlg.methods_SIMT   : 
            LinAlg.methods_SIMT  . clear ()
            LinAlg.methods_SIMT  = None 

        if LinAlg.methods_POW    : 
            LinAlg.methods_POW   . clear ()
            LinAlg.methods_POW   = None 

        if LinAlg.methods_SYM    : 
            LinAlg.methods_SYM   . clear ()
            LinAlg.methods_SYM   = None 


        if LinAlg.methods_ASYM   : 
            LinAlg.methods_ASYM  . clear ()
            LinAlg.methods_ASYM  = None 

        if LinAlg.method_EIGEN   : 
            LinAlg.method_EIGEN  . clear ()
            LinAlg.method_EIGEN  = None 

        if LinAlg.method_TM      : 
            LinAlg.method_TM     . clear ()
            LinAlg.method_TM     = None 

        if LinAlg.methods_EQ      : 
            LinAlg.methods_EQ     . clear ()
            LinAlg.methods_EQ     = None
            
        if LinAlg.methods_INV     : 
            LinAlg.methods_INV    . clear ()
            LinAlg.methods_INV    = None 

        return


    mgetter = staticmethod ( mgetter  )

    # =========================================================================
    ## create matrix/vector from numpy array/matrix  
    @staticmethod
    def NUMPY ( arr ) :
        pass
    
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
        """ Addition of vector/matrix objects:
        >>>  C = A + B
        """
        
        if isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            sa = a.shape
            sb = b.shape
            if sa != sb : raise NotImplementedError ( "Cannot  add %s/%s with %s" % ( typename  ( a ) , sa , sb ) )            
            return LinAlg.ADD ( a , LinAlg.toSObject ( b ) ) 
                
        operation , check  = LinAlg.methods_ADD ( a , b )        
        if operation and check and check ( a, b ) :
            result = operation ( a, b )
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
        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            sa = a.shape
            sb = b.shape
            if sa != sb : raise NotImplementedError ( "Cannot iadd %s/%s with %s" % ( typename  ( a ) , sa , sb ) )            
            return LinAlg.IADD ( a , LinAlg.toSObject ( b ) ) 

        operation , check  = LinAlg.methods_IADD ( a , b )
        if operation and check and check ( a, b ) :
            r = operation ( a, b )
            return a 

        return NotImplemented 
                
    # =========================================================================
    ## Right-addition of matrix/vector objects
    #  @code
    #  C = B + A  
    #  @endcode 
    @staticmethod
    def RADD ( a  , b ) :
        """ Right-addition of vector/,atrix objects:
        >>> C = B + A
        """
        if   isinstance ( b , num_types ) : b = float( b )

        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            
            return NotImplemented
        
            ## s1 = a.shape
            ## s2 = b.shape            
            ## if s1 != s2 : return NotImplemented
            ## return LinAlg.toSObject ( b ) + a 
        
        operation , check  = LinAlg.methods_RADD ( a , b )        
        if operation and check and check ( a, b ) :
            result = operation ( a, b )
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
        
        if   isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :            
            sa = a.shape
            sb = b.shape
            if sa != sb : raise NotImplementedError ( "Cannot  sub %s/%s with %s" % ( typename  ( a ) , sa , sb ) )            
            return LinAlg.SUB ( a , LinAlg.toSObject ( b ) ) 
        
        operation , check  = LinAlg.methods_SUB ( a , b )        
        if operation and check and check ( a, b ) :
            result = operation ( a, b )
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
                
        if   isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :            
            sa = a.shape
            sb = b.shape
            if sa != sb : raise NotImplementedError ( "Cannot isub %s/%s with %s" % ( typename  ( a ) , sa , sb ) )            
            return LinAlg.ISUB ( a , LinAlg.toSObject ( b ) )

        operation , check  = LinAlg.methods_ISUB ( a , b )        
        if operation and check and check ( a, b ) :
            r = operation ( a, b )
            return a 
        
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
        if   isinstance ( b , num_types ) : b = float( b )

        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :

            return NotImplemented

            ## s1 = a.shape
            ## s2 = b.shape
            ## if s1 != s2 : return NotImplemented            
            ## return LinAlg.toSObject ( b ) - a

        operation , check  = LinAlg.methods_RSUB ( a , b )        
        if operation and check and check ( a, b ) :
            result = operation ( a, b )
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
        
        if   isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            sa = a.shape
            sb = b.shape
            if sa [ -1 ] != sb [ 0 ]:
                raise NotImplementedError ( "Cannot  mul %s/%s with %s" % ( typename  ( a ) , sa , sb ) )            
            return LinAlg.MUL ( a , LinAlg.toSObject ( b ) )
        
        operation , check  = LinAlg.methods_MUL ( a , b )
        if operation and check and check ( a, b ) :
            result = operation ( a, b )
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
        
        if   isinstance ( b , num_types ) : b = float( b )        
        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            sa = a.shape
            sb = b.shape
            if sa [ -1 ] != sb [ 0 ] or 2 != len ( sb ) or sb [ 0 ] != sb [ 1 ] :
                raise NotImplementedError ( "Cannot imul %s/%s with %s" % ( typename  ( a ) , sa , sb ) )            
            return LinAlg.IMUL ( a , LinAlg.toSObject ( b ) )

        operation , check  = LinAlg.methods_IMUL ( a , b )
        if operation and check and check ( a, b ) :
            r = operation ( a, b )
            return a 
        
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

        if   isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            
            return NotImplemented

            ## sa = a.shape            
            ## sb = b.shape
            ## if sb [ -1 ] != sa[0] :
            ##     raise NotImplementedError ( "Cannot multiply %s/%s with %s" % ( typename ( a ), sa , sb ) )            
            ## return LinAlg.toSObject ( b ) * a 

        operation , check  = LinAlg.methods_RMUL ( a , b )        
        if operation and check and check ( a, b ) :
            result = operation  ( a, b )
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
        
        operation ,  check  = LinAlg.methods_DIV ( a , b )        
        if operation and check and check ( a, b ) :
            result = operation ( a, b )
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
        
        operation , check  = LinAlg.methods_IDIV ( a , b )        
        if operation and check and check ( a, b ) :
            r = operation ( a, b )
            return a 
        
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

        if isinstance ( b , num_types ) : b = float( b )

        ## numpy 
        if LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            s1 = a.shape            
            s2 = b.shape
            if s1 != s2 : return NotImplemented 
            return np.array_equal ( a.to_numpy() , b )

        ## vector-like stuff 
        if isinstance ( b , ( list , tuple , array.array ) ) and 1 == len ( a.shape ) : 
            if len ( a ) != len ( b )   : return NotImplemented 
            for i , j in zip ( a , b )  :
                if i == j or isequal ( i , j ) : continue 
                return False
            return True 

        operation , check  = LinAlg.methods_EQ ( a , b )
        if operation and check and check ( a, b ) :
            result = operation ( a, b )
            return result
        
        operation , check  = LinAlg.methods_EQ ( b , a )
        if operation and check and check ( b , a ) :
            result = operation ( b , a )
            return result

        return NotImplemented 

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

        if isinstance ( b , num_types ) : b = float( b )

        ## numpy 
        if LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            s1 = a.shape            
            s2 = b.shape
            if s1 != s2 : return NotImplemented 
            return not np.array_equal ( a.to_numpy() , b )

        ## vector-like stuff 
        if isinstance ( b , ( list , tuple , array.array ) ) and 1 == len ( a.shape ) : 
            if len ( a ) != len ( b )   : return NotImplemented
            for i , j in zip ( a , b )  :
                if i != j or not isequal ( i , j ) : return True 
            return False 
        
        operation , check  = LinAlg.methods_EQ ( a , b )
        if operation and check and check ( a, b ) :
            result = operation ( a, b )
            return not result
        
        operation , check  = LinAlg.methods_EQ ( b , a )
        if operation and check and check ( b , a ) :
            result = operation ( b , a )
            return not result

        return NotImplemented 

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
                raise NotImplementedError ( "No DOT for %s/%s and %s/%s" % ( a , typename  ( a ), b , typename  ( b ) ) )
            return np.dot ( a.to_numpy() , b )
            
        
        if isinstance ( b , num_types ) : b = float( b )        
        operation , check  = LinAlg.methods_DOT ( a , b )        
        if operation and check and check ( a, b ) :
            result = operation ( a, b )
            return result
        
        raise NotImplementedError ( "No DOT for %s/%s and %s/%s" % ( a , typename ( a ) , b , typename ( b ) ) )

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
            if 1 != len ( s1 ) or 1 != len ( s1 ) :
                raise NotImplementedError ( "No CROSS for %s/%s and %s/%s" % ( a , typename ( a ) , b , typename ( b ) ) )
            return LinAlg.CROSS ( a ,  LinAlg.toSVector ( b ) ) 
            
        operation , check  = LinAlg.methods_CROSS ( a , b )
        if operation and check and check ( a, b ) :
            result = operation ( a, b )
            return result
        
        raise NotImplementedError ( "No CROSS for %s/%s and %s/%s" % ( a , typename ( a ) , b , typename ( b ) ) )

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

        if LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            sa = a.shape            
            sb = b.shape
            ## square matrix
            if 2 == len ( sa ) and sa[0] == sa [1] :
                if 1   == len ( sb ) and sb [ 0 ] == sa [ 0 ] :
                    bb = LinAlg.toSVector ( b )
                    return LinAlg.SIM ( a , bb )  
                elif 2 == len ( sb ) and sb [ 1 ] == sa [ 0 ] :  
                    bb = LinAlg.toSMatrix ( b )
                    return LinAlg.SIM ( a , bb )  
            raise NotImplementedError ( "Cannot  sim %s/%s with %s" % ( typename  ( a ) , sa , sb ) )            
        
        if isinstance ( b , num_types ) : b = float( b )
        
        operation , check  = LinAlg.methods_SIM ( a , b )        
        if operation  and check and check ( a, b ) :
            result = operation ( a, b )
            return result
        
        raise NotImplementedError ( "No SIM for %s/%s and %s/%s" % ( a , typename ( a ) , b , typename ( b ) ) )

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

        if LinAlg.with_numpy and isinstance ( b , np.ndarray ) :
            sa = a.shape            
            sb = b.shape
            ## square matrix
            if 2 == len ( sa ) and sa [ 0 ] == sa [1] :
                if 2 == len ( sb ) and sb [ 0 ] == sa [ 0 ] :  
                    bb = LinAlg.toSMatrix ( b )
                    return LinAlg.SIMT ( a , bb )      
            raise NotImplementedError ( "Cannot simT %s/%s with %s" % ( typename  ( a ) , sa , sb ) )
        
        if isinstance ( b , num_types ) : b = float( b )
        
        operation , check  = LinAlg.methods_SIMT ( a , b )
        if operation and check and check ( a, b ) :
            result = operation ( a, b )
            return result
        
        raise NotImplementedError ( "No SIMT for %s/%s and %s/%s" % ( a , typename ( a ) , b , typename ( b ) ) )


    # =========================================================================
    ## power function for square matrices 
    #  @code
    #  C = A ** 6 
    #  @endcode 
    @staticmethod
    def M_POW ( a  , n ) :
        """ Power function for the square matrices  
        >>>  C = A ** 6 
        """

        if   isinstance ( n , integer_types ) : pass 
        elif isinstance ( n , num_types     ) :
            n = float ( n )
            if Ostap.Math.isint ( n ) : n = int ( n ) 
        else :
            return NotImplemented 

        ##  square matrix ?
        square = a.kRows == a.kCols
        if not square : return NotImplemented
        
        ## 1.  1x1 matrix : exponent can be not onl integer 
        if   1 == a.kRows                         : return pow ( a ( 0 , 0 ) , n )
        ## 2.  non-integer exponent?
        elif not isinstance ( n , integer_types ) : return NotImplemented
        
        ## 3 negative integer exponent : first invert and then pow  
        if n < 0 :
            try : 
                a_inv = LinAlg.M_INVERSE ( a )
            except ValueError : ## matrix cannot be inverted
                return NotImplemented
            if -1 == n : return a_inv 
            return LinAlg.M_POW ( a_inv , abs ( n ) )

        ## 4. regular case: square matrix and non-negative integer exponent 
        operation , check  = LinAlg.methods_POW ( a  )
        if operation and check and check ( a , n ) :
            result = operation ( a , n )
            return result
        
        return NotImplemented 

    # =========================================================================
    ## symmetric part of square marix
    #  @code
    #  C = A.sym() 
    #  @endcode 
    @staticmethod
    def M_SYM ( a  ) :
        """ Symmetric part of square matrices  
        >>>  C = A.sym() 
        """
        
        operation , check  = LinAlg.methods_SYM ( a  )
        if operation and check and check ( a ) :
            result = operation ( a )
            return result
        
        raise NotImplementedError ( "Cannot symmetrise %s/%s" % ( a , typename ( a ) ) )
    
    # =========================================================================
    ## Antisymmetric/skew part of square marix
    #  @code
    #  C = A.asym() 
    #  C = A.skew() 
    #  @endcode 
    @staticmethod
    def M_ASYM ( a  ) :
        """ Antisymmetric/skew part of square matrices  
        >>>  C = A.asym() 
        >>>  C = A.skew() 
        """
        
        operation , check  = LinAlg.methods_ASYM ( a  )
        if operation and check and check ( a ) :
            result = operation ( a  )
            return result
        
        raise NotImplementedError ( "Cannot anti-symmetrise %s/%s" % ( a , typename ( a ) ) )
    
    
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
        s = len ( vct )
        for i in range ( s ) : yield vct ( i )

    # =============================================================================
    ## iterator for SVector
    #  @code
    #  vct = ...
    #  for i in vct.keys() : print i 
    #  @endcode
    @staticmethod
    def V_KEYS ( vct ) :
        """Iterator for SVector
        >>> vct = ...
        >>> for i in vct.keys() : print i 
        """
        s = len ( vct )
        for i in range ( s ) : yield i 
            
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
    ## self-printout of S-vectors-with-errors
    #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    #  @date 2020-05-15
    @staticmethod
    def VE_STR ( vct , fmt = ' %g' ) :
        """Self-printout of SVectorsWithError: (...)
        """

        v = vct.value      ()
        c = vct.covariance ()
        return "%s\n%s" % ( v , c )

    # =============================================================================
    ## Transform vector-with-errors to another variable
    #  for \f[ y = y ( x ) , C(y) = J C(x) J^T, \f]
    #  where  \f$ J = \left( \frac{\partial y}{\partial x }\right) \f$.
    #
    #  Transofrm to (r,phi)-varibales: 
    #  @code
    #  vct = ...  ##
    #  r   = lambda x,y : (x*x+y*y)**0.5
    #  phi = lambda x,y : math.atan2  ( y , x )
    #  result = vct.tranform  ( r  , phi ) 
    #  @endcode    
    #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    #  @date 2020-05-15
    @staticmethod
    def VE_TRANSFORM ( vct , *Y ) :
        """Transform vector-with-errors to another variable
         - y = y ( x )
         - C(y) = J C(x) J^T
         where J = dy/dx
         - Transofrm to (r,phi)-varibales: 
         >>> vct = ...  ##
         >>> r   = lambda x,y : (x*x+y*y)**0.5
         >>> phi = lambda x,y : math.atan2  ( y , x )
         >>> result = vct.tranform  ( r  , phi ) 
        """
        assert Y , 'Tranform: At least one function need to be specified!'

        ## variables in form of tuple 
        x = tuple ( [ v for v in vct ] )  ## values in form of tuple 
        
        y = []        
        for i , f in enumerate ( Y ) :
            assert callable ( f ) , 'Y(%d) is not calllable!' % i  
            try :
                _ = f ( *x )
            except :
                logger.error("Transform: Y(%s) can not be called with %s" % ( i , str ( x ) ) ) 
                raise
            y.append ( f )
            
        y = tuple  ( y ) 

        nx = len ( x )
        ny = len ( y )
        
        from ostap.math.derivative  import PartialDerivative as PD 
        from ostap.math.ve          import VE 

        ## Jacobi matrix  
        J  = LinAlg.Matrix ( ny , nx ) ()
        for j , f  in enumerate ( y ) : 
            for i in range ( nx ) :
                P = PD ( i , f  )            
                J [ j , i ] = P ( *x )

        ##  result 
        res = Ostap.Vector ( ny ) () 
        for i , f in enumerate ( y ) : res [ i ] = f ( *x  ) 

        ## calcaulte covariance  
        c2old = vct.covariance()
        c2new = c2old.Sim ( J )
        
        if 1 == ny : return VE ( res[0] , c2new ( 0 , 0 ) )
        
        return Ostap.VectorE ( ny ) ( res , c2new )
        
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
                yield mtrx ( i , j )

    # =============================================================================
    ## iterator for SMatrix
    #  @code
    #  matrix = ...
    #  for i in matrix.keys() : print i 
    #  @endcode
    @staticmethod 
    def M_KEYS ( mtrx ) :
        """Iterator for SMatrix
        >>> matrix = ...
        >>> for i in matrix.keys() : print i 
        """
        krows = mtrx.kRows
        kcols = mtrx.kCols
        for i in range ( krows ) :
            for j in range ( kcols ) :
                yield ( i , j )

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
            
            for j in range ( i + 1 , rows ) :            
                jj  = mtrx ( j , j ) 
                sjj = sqrt ( jj    )
                ij  = mtrx ( i , j )
                eij = ij / ( sii * sjj )
                if  1 < abs ( eij ) : ok2 = False  
                _c [ i , j ] = eij 
                _c [ j , i ] = eij 
            
        if not ok1 : logger.error ( "correlations: zero or negative diagonal element" ) 
        if not ok2 : logger.error ( "correlations: invalid non-diagonal element"      ) 
                
        return _c


    # =========================================================================
    ## Get the inverse matrix
    #  @code
    #  m     = ...
    #  m_inv = m.inverse () 
    #  @endcode
    @staticmethod 
    def M_INVERSE ( mtrx ) :
        """ Get the inverse matrix
        >>> m     = ...
        >>> m_inv = m.inverse () 
        """

        square = mtrx.kRows == mtrx.kCols 
        
        if not square :
            raise NotImplementedError ('Inversion is defined only for square matrices!') 

        ## regular case: square matrix and non-negative integer exponent 
        operation , check  = LinAlg.methods_INV ( mtrx )
        if operation and check and check ( mtrx ) :
            flag   = ctypes.c_int(0)
            result = operation ( mtrx , flag )
            if 0  != flag.value :
                raise ValueError('Matrix cannot be inverted %s (invalid/degenerated,...)' % typename ( mtrx ) )
            return result
    
        raise NotImplementedError ('Matrix inversion is not possible here %s' % typename ( mtrx ) ) 

        
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
        
        operation  = LinAlg.method_EIGEN ( mtrx  )
        if not operation  : raise NotImplementedError ("EigenValues: not implemented for %s" % typename ( mtrx ) )
        
        values = LinAlg.Vector ( mtrx.kRows )()
        st     = operation ( mtrx  , values , sorted )
        assert st.isSuccess () , "Eigen-values: status code %s" % st
        ##
        return values 

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
        
        operation = LinAlg.method_EIGEN ( mtrx )
        if not operation :
            raise NotImplementedError ("EigenVectors: not implemented for %s" % typename ( mtrx ) )
        
        krows   = mtrx.kRows
        kcols   = mtrx.kCols        
        values  = LinAlg.Vector ( krows        ) ()
        vectors = LinAlg.Matrix ( krows, kcols ) () 
        st      = operation ( mtrx  , values , vectors , sorted )
        assert st.isSuccess () , "Eigen-vectors: status code %s" % st
        ## 
        return values, vectors  

    # =========================================================================
    ## reduce SVector
    @staticmethod
    def V_REDUCE ( vct ) :
        """Reduce SVector"""
        return svct_factory, ( array.array ( 'd' , vct ) , )

    # =========================================================================
    ## reduce SMatrix 
    @staticmethod
    def M_REDUCE ( mtrx ) :
        """Reduce SMatrix
        """
        NR   = mtrx.rep_size
        a    = mtrx.Array() 
        data = array.array ( 'd' , ( a[i] for i in range ( NR ) ) ) 
        return smtrx_factory, ( mtrx.kRows, mtrx.kCols , data ) 

    # =========================================================================
    ## reduce symmetric SMatrix 
    @staticmethod
    def MS_REDUCE ( mtrx ) :
        """Reduce symmetric SMatrix
        """
        NR   = mtrx.rep_size 
        a    = mtrx.Array() 
        data = array.array ( 'd' , ( a[i] for i in range ( NR ) ) ) 
        return symmm_factory, ( mtrx.kRows , data ) 

    # =========================================================================
    ## reduce SVectorWithErrors 
    @staticmethod
    def VE_REDUCE ( vct ) :
        """reduce SVectorWithErrors"""
        return svcte_factory , ( vct.value() , vct.cov2() ) 
        
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
        operation = LinAlg.method_TM ( sobj )
        if not operation : raise NotImplementedError ("SMatrix->TMatrix/SVector->TVector: not implemented for %s" % typename ( sobj ) )
        
        tobj = operation ( sobj ) 
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
        
        if ( 3 , 5 ) <= python_version :
            
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
        t.  keys        = LinAlg.V_KEYS 
        t. ikeys        = LinAlg.V_KEYS 
        
        t.to_array      = LinAlg.V_ARRAY ## plain array.array 

        t.tvector       = LinAlg.M_TM 

        t.shape         = property ( LinAlg.V_SHAPE , None , None )

        t.__reduce__    = LinAlg.V_REDUCE 
        
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
        
        R = m.rep_type

        m. __add__      = lambda a , b : LinAlg. ADD ( a , b ) 
        m.__radd__      = lambda a , b : LinAlg.RADD ( a , b ) 
        m.__iadd__      = lambda a , b : LinAlg.IADD ( a , b )
        
        m. __sub__      = LinAlg. SUB 
        m.__rsub__      = LinAlg.RSUB 
        m.__isub__      = LinAlg.ISUB 

        m. __mul__      = LinAlg. MUL 
        m.__rmul__      = LinAlg.RMUL 
        m.__imul__      = LinAlg.IMUL 

        if ( 3 , 5 ) <= python_version : 
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
        m. keys         = LinAlg.M_KEYS  
        m.ikeys         = LinAlg.M_KEYS  
        m.__contains__  = lambda s,ij : 0<=ij[0]<s.kRows and 0<=ij[1]<s.kCols

        m.row           = LinAlg.M_ROW
        m.column        = LinAlg.M_COLUMN
        m.rows          = LinAlg.M_ROWS
        m.columns       = LinAlg.M_COLUMNS

        m.tmatrix       = LinAlg.M_TM 

        if m.kRows == m.kCols :
            m.inverse   = LinAlg.M_INVERSE            

        ## conversion to float 
        if 1 == m.kRows and 1 == m.kCols :
            m.__float__ = lambda s : s ( 0 , 0 )

        ## should be class property!! 
        m.shape         = property ( LinAlg.M_SHAPE , None , None )

        m.__pow__       = LinAlg.M_POW 
        m.sym           = LinAlg.M_SYM
        m.asym          = LinAlg.M_ASYM 
        m.skew          = LinAlg.M_ASYM 


        m.__reduce__    = LinAlg.M_REDUCE 
        m.rep_size      = classgetter ( lambda cls : cls.rep_type.kSize ) 

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
        
        R = m.rep_type
        
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

        m.__reduce__     = LinAlg.MS_REDUCE 
        m.rep_size       = classgetter ( lambda cls : cls.rep_type.kSize ) 
        
        return m

    # =========================================================================
    ## Decorate SVectorWithError 
    @staticmethod
    def deco_vectore ( t ) :
        """ Decorate SVectorWithError
        """
        
        if t in LinAlg.decorated_vectors : return t 

        LinAlg.decorated_vectors.add ( t )

        LinAlg.backup ( t )
        
        t. __str__      = LinAlg.VE_STR
        t. __repr__     = LinAlg.VE_STR

        t. __len__      = lambda s : s.kSize 
        t. __contains__ = lambda s, i : 0 <= i < s.kSize

        t. __iter__     = LinAlg.V_ITER      
        t.transform     = LinAlg.VE_TRANSFORM

        t.__reduce__    = LinAlg.VE_REDUCE 
                
        return t


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
               'Invalid length of the vector %s/%s' % ( n , typename ( l ) )

        tt = n , t 
        v  = LinAlg.known_svectors.get ( tt , None )
        if  v is None :
            v = ROOT.ROOT.Math.SVector ( t , n )
            LinAlg.known_svectors [ tt ] = v
            LinAlg.deco_vector    ( v )
            ##
            if not tt  in LinAlg.known_ssymmatrices : LinAlg.SymMatrix ( n ,     t )    
            if not tt  in LinAlg.known_svectorse    : LinAlg.  VectorE ( n ,     t )    
            tt1 = n , n , t 
            if not tt1 in LinAlg.known_smatrices    : LinAlg.   Matrix ( n , n , t )

            ## create also all smaller vectors 
            for i in range ( 2 , n ) :
                ti = i , t
                if not ti in LinAlg.known_svectors : LinAlg.Vector ( i ,  t )                
                
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
               'Invalid matrix dimension %s/%s' % ( k , typename ( k ) )
        assert isinstance  ( n , integer_types ) and 0 <= n ,\
               'Invalid matrix dimension %s/%s' % ( n , typename ( n ) )

        tt = k , n , t
        m = LinAlg.known_smatrices.get ( tt , None )
        if m is None :
            
            m = ROOT.ROOT.Math.SMatrix ( t , k , n )
            LinAlg.known_smatrices [ tt ] = m
            LinAlg.deco_matrix ( m )

            ## check k-dimension 
            tt1 = k , t
            if not tt1 in LinAlg.known_svectors     : LinAlg.   Vector ( k , t )
            if not tt1 in LinAlg.known_ssymmatrices : LinAlg.SymMatrix ( k , t )    

            ## check N-dimension
            tt2 = n , t
            if not tt2 in LinAlg.known_svectors     : LinAlg.   Vector ( n , t )
            if not tt2 in LinAlg.known_ssymmatrices : LinAlg.SymMatrix ( n , t )    

            ## check transposed matrix 
            tt3 = n , k , t
            if not tt3 in LinAlg.known_smatrices    : LinAlg.   Matrix ( n , k , t ) 

            ## create also smaller vectors 
            for i in range ( 2 , max ( k , n )  ):
                ti = i , t
                if not ti in LinAlg.known_svectors  : LinAlg.Vector ( i ,  t )                
                
            ## create also smaller matrices 
            for i  in range ( 2 , k + 1 ) :
                for j in range ( 2 , n + 1 ) :
                    tti = i , j , t
                    if not tti in LinAlg.known_smatrices :  LinAlg.   Matrix ( i , j , t ) 
                    
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
               'Invalid matrix dimension %s/%s' % ( n , typename ( n ) )

        tt = n , t
        m = LinAlg.known_ssymmatrices.get ( tt , None )
        if  m is None  : 
            m = ROOT.ROOT.Math.SMatrix('%s,%d,%d,ROOT::Math::MatRepSym<%s,%d>' %  ( t , n , n , t , n ) )
            LinAlg.known_ssymmatrices [ tt ] = m
            LinAlg.deco_symmatrix  ( m )
            ##
            if not tt  in LinAlg.known_svectors  : LinAlg.Vector ( n ,     t )

            tt1 = n , n , t 
            if not tt1 in LinAlg.known_smatrices : LinAlg.Matrix ( n , n , t )

            ## create also all smaller vectors and matrices 
            for i in range ( 2 , n ):
                ti = i , t
                if not ti in LinAlg.known_svectors     : LinAlg.   Vector ( i , t )                
                if not ti in LinAlg.known_ssymmatrices : LinAlg.SymMatrix ( i , t )

            ## create also smaller matrices 
            for i in range ( 2 , n + 1 ) :
                tv = i , t
                if not tv in LinAlg.known_svectors      : LinAlg.Vector ( i ,  t )
                for j in range ( i , n + 1 ) :
                    tti = i , j , t
                    if not tti in LinAlg.known_smatrices :  LinAlg. Matrix ( i , j , t ) 
                                    
        return m
    
    # =========================================================================
    ##  Pick up the vector-with-errors
    #   @code
    #   VE3  = Ostap.Math.VectorE(3)
    #   vctE = VE3 ()
    #   @endcode
    #   @see Ostap::Math::SVrctorWithError
    @staticmethod
    def VectorE ( n , t = 'double' ) :
        """Pick up the vector-with-error  of corresponding size
        >>> VE3  = Ostap.Math.VectorE(3)
        >>> vctE = VE3 ()
        """
        assert isinstance  ( n , integer_types ) and 0 <= n,\
               'Invalid length of the vector %s/%s' % ( n , typename ( l ) )

        tt = n , t
        v = LinAlg.known_svectorse.get ( tt , None )
        if  v is None :
            v = Ostap.Math.SVectorWithError ( n , t )
            LinAlg.known_svectorse [ tt ] = v
            LinAlg.deco_vectore    ( v )
            ##
            if not tt  in LinAlg.known_svectors     : LinAlg.   Vector ( n ,     t )    
            if not tt  in LinAlg.known_ssymmatrices : LinAlg.SymMatrix ( n ,     t )    
            tt1 = n , n , t 
            if not tt1 in LinAlg.known_smatrices    : LinAlg.   Matrix ( n , n , t )

            for i in range ( 2 , n ) :
                ti = i , t
                if not ti in LinAlg.known_svectorse  : LinAlg. VectorE ( i , t )                
                
        return v 

    # =========================================================================
    ##  is this matrix/vector type decorated properly?
    @staticmethod
    def decorated ( t ) :
        """Is this matrix/vector type decorated properly?
        """        
        return t in LinAlg.decorated_vectors or t in LinAlg.decorated_matrices
    



# =============================================================================

Ostap.Vector         =  staticmethod ( LinAlg.Vector    )
Ostap.VectorE        =  staticmethod ( LinAlg.VectorE   )
Ostap.Matrix         =  staticmethod ( LinAlg.Matrix    )
Ostap.SymMatrix      =  staticmethod ( LinAlg.SymMatrix ) 
Ostap.Math.Vector    =  staticmethod ( LinAlg.Vector    )
Ostap.Math.VectorE   =  staticmethod ( LinAlg.VectorE   ) 
Ostap.Math.Matrix    =  staticmethod ( LinAlg.Matrix    )
Ostap.Math.SymMatrix =  staticmethod ( LinAlg.SymMatrix ) 

# =============================================================================


if np :

    # =========================================================================
    ## create <code>SMatrix</code> from <code>numpy</code> array
    #  @code
    #  a = ...
    #  m = toSMatrix ( a ) 
    #  @endcode 
    def toSMatrix ( a , shape = None ) :
        """ Create `SMatrix` from 2D `numpy` array
        >>> a = ...
        >>> m = toSMatrix ( a ) 
        """
        assert isinstance ( a , np.ndarray ) and 2 == len ( a.shape ) ,\
               'toSMatrix: Invalid type/shape of input array'
        
        if shape :
            assert shape == a.shape, 'toSMatrix: Invalid shape of input array'

        mtrx  = LinAlg.Matrix ( *a.shape ) () 
        N , K = a.shape 
        for i in range ( N ) :
            for j in range ( K ) :
                index = i , j 
                mtrx [ index ] = a [ index ]
                
        return mtrx 

    # =========================================================================
    ## create <code>SVector</code> from <code>numpy</code> array
    #  @code
    #  a = ...
    #  v = toSVector ( a ) 
    #  @endcode 
    def toSVector ( a , length = None ) :
        """ Create `SVector` from 1D `numpy` array
        >>> a = ...
        >>> v = toSVector( a ) 
        """
        assert isinstance ( a , np.ndarray ) and 1 == len ( a.shape ) ,\
               'toSVector: Invalid type/shape of input array'
        
        if not length is None :
            assert shape == a.shape, 'toMatrix: Invalid shape/size of input array'

        N   = a.shape[0]
        vct = LinAlg.Vector ( N ) () 
        for i in range ( N ) :
            vct [ i ] = a [ i ]

        ## vct = LinAlg.Vector ( a.shape [ 0 ] ) ( *tuple ( a  ) )
            
        return vct

    # =========================================================================
    ## convert numpy array to SMatrix/SVecrtor
    #  @code
    #  a = ...
    #  s = toSObject ( a ) 
    #  @endcode
    def toSObject ( a ) :
        """Convert 1D or 2D numpy array to SMatrix/SVecrtor
        >>> a = ...
        >>> s = toSObject ( a ) 
        """
        assert isinstance ( a , np.ndarray ) and 1 <= len ( a.shape ) <= 2 ,\
               'toSObject: Invalid type/shape of input array'
        return toSVector ( a ) if 1 == len( a.shape ) else toSMatrix ( a )


    LinAlg.toSMatrix      = staticmethod ( toSMatrix ) 
    LinAlg.toSVector      = staticmethod ( toSVector ) 
    LinAlg.toSObject      = staticmethod ( toSObject ) 

    Ostap.toSMatrix       = staticmethod ( LinAlg.toSMatrix ) 
    Ostap.toSVector       = staticmethod ( LinAlg.toSVector )
    Ostap.toSObject       = staticmethod ( LinAlg.toSObject )
    Ostap.Math.toSMatrix  = staticmethod ( LinAlg.toSMatrix ) 
    Ostap.Math.toSVector  = staticmethod ( LinAlg.toSVector )
    Ostap.Math.toSObject  = staticmethod ( LinAlg.toSObject )



# =============================================================================
## factory for vectors
def svct_factory ( data ) :
    """Factory for vectors
    """
    N  = len ( data ) 
    VN = LinAlg.Vector( N )
    v  = VN()
    v.SetElements ( data , N )
    return v

# =============================================================================
## Factory for matrices
def smtrx_factory ( N , K , data ) :
    """ Factory for matrices
    """
    
    MNK = LinAlg.Matrix ( N , K )
    m   = MNK()
    a   = m.Array()
    for i , d in enumerate ( data ) :
        a [ i ] = d
    return m

# =============================================================================
## factory for symmetric matrices
def symmm_factory ( N , data ) :
    """Factory for symmetric matrices
    """    
    SNN = LinAlg.SymMatrix ( N )
    m   = SNN ()
    a   = m.Array()
    for i , d in enumerate ( data ) :
        a [ i ] = d
    return m

# =============================================================================
## factory for vectors with errors
def svcte_factory ( vals , cov2 ) :
    """factory for vectors with errors
    """
    N    = len ( vals )
    SVEN = LinAlg.VectorE ( N )
    return SVEN ( vals , cov2 ) 

 
# =============================================================================
_decorated_classes_ = (
    )

_decorated_classes_ = _decorated_classes_ + tuple ( LinAlg.decorated_vectors  )
_decorated_classes_ = _decorated_classes_ + tuple ( LinAlg.decorated_matrices )

_new_methods_ = (
    Ostap.Vector         , 
    Ostap.VectorE        , 
    Ostap.Matrix         , 
    Ostap.SymMatrix      , 
    Ostap.Math.Vector    , 
    Ostap.Math.VectorE   , 
    Ostap.Math.Matrix    , 
    Ostap.Math.SymMatrix , 
    )

if np :
    _new_methods_ = _new_methods_ + (
        Ostap.toSMatrix      ,
        Ostap.toSVector      ,
        Ostap.toSObject      ,
        Ostap.Math.toSMatrix ,
        Ostap.Math.toSVector ,
        Ostap.Math.toSObject ,
        )

import atexit
atexit.register ( LinAlg.CLEANUP ) 

# =============================================================================
## check what LinAlg operations are defined for these two objects
#  @code
#  obj1 = ...
#  obj2 = ...
#  checkops ( obj1 , obj2 ) 
#  @endcode 
def checkops ( a , b , logger = logger ) :
    """check what LinAlg operations are defined for these two objects    
    >>> obj1 = ...
    >>> obj2 = ...
    >>> checkops ( obj1 , obj2 ) 
    """

    rows = [ ( 'Method' , 'checker' , 'operation' , 'ok' , 'result' ) ]

    methods = ( 
        ( '+'     , LinAlg.methods_ADD   ) ,
        ( '+/r'   , LinAlg.methods_RADD  ) ,
        ( '+='    , LinAlg.methods_IADD  ) ,
        ##
        ( '-'     , LinAlg.methods_SUB   ) ,
        ( '-/r'   , LinAlg.methods_RSUB  ) ,
        ( '-='    , LinAlg.methods_ISUB  ) ,
        ##
        ( '*'     , LinAlg.methods_MUL   ) ,
        ( '*/r'   , LinAlg.methods_RMUL  ) ,
        ( '*='    , LinAlg.methods_IMUL  ) ,
        ##
        ( '/'     , LinAlg.methods_DIV   ) ,
        ( '/='    , LinAlg.methods_IDIV  ) ,
        ##
        ( '=='    , LinAlg.methods_EQ    ) ,
        ## 
        ( 'dot'   , LinAlg.methods_DOT   ) ,
        ( 'cross' , LinAlg.methods_CROSS ) ,
        ##
        ( 'sim'   , LinAlg.methods_SIM   ) ,
        ( 'simt'  , LinAlg.methods_SIMT  ) )

    
    for symbol , method in methods :
        
        operation , checker = method ( a , b )  
        
        result = ''
        if checker : ok = True if checker ( a, b ) else False
        else       : ok = ''
        result = '%s' % type ( operation ( a , b ) ).__name__ if operation and ok else '' 
        
        row = ( symbol                       ,
                'ok' if operation else '---' ,
                'ok' if checker   else '---' ,
                '%s' % ok ,
                result )
        rows.append ( row )
        
    methods2 = ( 
        ( '**'   , LinAlg.methods_POW   ) ,
        )

    methods3 = ( 
        ( 'sym'  , LinAlg.methods_SYM   ) ,
        ( 'ssym' , LinAlg.methods_ASYM  ) ,
        )
    
    
    
    import ostap.logger.table as T
    title = 'Allowed binary operations'
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lcccl' )
    logger.info ( "%s for '%s' and '%s':\n%s" % ( title , type ( a ).__name__ , type ( b ).__name__ , table ) ) 
                      
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
