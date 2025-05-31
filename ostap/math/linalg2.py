#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/linalg.py
#  Few utilities to simplify linear algebra manipulations 
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
""" Few utilities to simplify linear algebra manipulations 
"""
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
from   ostap.core.ostap_types import num_types , integer_types
from   ostap.math.base        import isequal   , iszero, std , Ostap, numpy
from   ostap.core.core        import hID 
from   ostap.utils.basic      import typename 
from   ostap.utils.clsgetter  import classgetter
from   ostap.logger.pretty    import fmt_pretty_float, fmt_pretty_error 
from   ostap.math.base        import pretty_array
from   ostap.logger.colorized import infostr
from   ostap.utils.gsl        import gsl_info
from   ostap.logger.symbols   import ditto, times, labels 
from   ostap.logger.colorized import colored_string  
import ostap.logger.table     as     T
import ROOT, math, re, ctypes, array, random 
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.linalg2' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
revct = re.compile ( r'SVector<(?P<TYPE>[^,>]+)' )
remtx = re.compile ( r'SMatrix<(?P<TYPE>[^,>]+)' )
NaN   = float('nan')
# =============================================================================
def diag ( what ) : 
    return colored_string ( what             ,  
                           foreground = 0    , 
                           background = 7    ,
                           fg_bright  = True , 
                           underline  = True )
# =============================================================================
## Helper method: get  i,j element from matrix-like object
#  @code
#  mtrx  = ...
#  value = matrix ( mgetter , 1 , 2 ) 
#  @endcode 
def mgetter ( mtrx , i , j ) :
    """ Helper method: get  (i,j)  element from matrix-like object
    >>> mtrx  = ...
    >>> value = mgetter ( mtrx , 1 , 2 ) 
    """
    
    if callable ( mtrx  ) :
        
        kRows = getattr ( mtrx , 'kRows' , None )
        if kRows is None        :
            nRows = getattr ( mtrx , 'GetNrows' , None )
            if nRows and not 0 <= i < nRows() : raise IndexError  ( "Invalid row index %s" % i )            
        elif not 0 <= i < kRows : raise IndexError  ( "Invalid row index %s" % i )
        
        kCols = getattr ( mtrx , 'kCols' , None )
        if kCols is None        :
            nCols = getattr ( mtrx , 'GetNcols' , None )
            if nCols and not 0 <= i < nCols() : raise IndexError  ( "Invalid column index %s" % i )
        elif not 0 <= j < kCols : raise IndexError  ( "Invalid column index %s" % j )

        return mtrx ( i , j )
    
    if numpy and isinstance ( obj , numpy.ndarray ) and 2 == len ( obj.shape ) :
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
    """ Get the correlation element from the matrix-like object
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
    """ Access and keep the method
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
    """ Access and keep two methods
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
    """ Collection of decorators for vectors/matrices
    """

    with_numpy    = numpy

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
        """ Restore useful attributes
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
        """ Cleanup LinAlg
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
        """ Convert vector to numpy array:
        >>> vct = ...
        >>> na = vct.to_numpy() 
        """
        return numpy.array ( vct  , dtype = 'd' ) 
    
    # =========================================================================
    ## convert matrix into numpy.matrix
    #  @code
    #  mtrx = ...
    #  m    = mtrx.to_numpy() 
    #  @endcode 
    @staticmethod 
    def M_NUMPY ( mtrx ) :
        """ Convert matrix into numpy.matrix
        >>> mtrx = ...
        >>> m    = mtrx.to_numpy() 
        """
        
        npa = numpy.empty ( mtrx.shape , dtype = 'd' )         
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
        """ Convert matrix/vector-like object to SMatrix/SVector
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
        elif LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :
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
        """ Increment  of matrix/vector objects
        >>> A += B  
        """

        if isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :
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

        elif LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :
            
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
        elif LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :            
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
        elif LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :            
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

        elif LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :

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
        """ Multiplication/scaling of vector/matrix objects:
        >>>  C = A * B
        """
        
        if   isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :
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
        elif LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :
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
        """ Right-multiplication of matrix/vector objects
        >>> C = B * A 
        """

        if   isinstance ( b , num_types ) : b = float( b )
        elif LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :
            
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
        """ Multiplicative decrement/scaling of matrix/vector objects
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
        """ Equality for matrix/vector objects
        >>> A == B  
        """

        if isinstance ( b , num_types ) : b = float( b )

        ## numpy 
        if LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :
            s1 = a.shape            
            s2 = b.shape
            if s1 != s2 : return NotImplemented 
            return numpy.array_equal ( a.to_numpy() , b )

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
        """ Non-equality for matrix/vector objects
        >>> A != B  
        """

        if isinstance ( b , num_types ) : b = float( b )

        ## numpy 
        if LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :
            s1 = a.shape            
            s2 = b.shape
            if s1 != s2 : return NotImplemented 
            return not numpy.array_equal ( a.to_numpy() , b )

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
        
        if LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :
            s1 = a.shape            
            s2 = b.shape
            if s1 != s2  or 1 != len ( s1 ) :
                raise NotImplementedError ( "No DOT for %s/%s and %s/%s" % ( a , typename  ( a ), b , typename  ( b ) ) )
            return numpy.dot ( a.to_numpy() , b )
            
        
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
        """ Cross-product (D1xD2 matrix) of two vectors 
        >>> matrix = v1.cross ( v2 ) 
        """
        
        if LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :
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
        """ Similarity  operation: C = B A B^T 
        >>> A = ...
        >>> B = ...
        >>> C = A.Sim ( B )   
        >>> C = A.sim ( B ) ## ditto
        """

        if LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :
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
    ## Similarity operation \f$ C = B^T A B  \f$ 
    #  @code
    #  A = ...
    #  B = ...
    #  C = A.SimT ( B )  
    #  C = A.simT ( B )   ## ditto 
    #  @endcode 
    @staticmethod
    def SIMT ( a  , b ) :
        """ Similarity  operation C = B^T A B 
        >>> A = ...
        >>> B = ...
        >>> C = A.SimT ( B )   
        >>> C = A.simT ( B ) ## ditto 
        """

        if LinAlg.with_numpy and isinstance ( b , numpy.ndarray ) :
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

        ## 1. zero exponent: return identity matrix 
        if   0 == n or iszero ( n ) :
            if isinstance ( a , ROOT.TObject ) : return type ( a ) ( a.kUnit , a ) 
            return type ( a ) ( ROOT.ROOT.Math.SMatrixIdentity() )        
        ## 2. unit  exponent: return copy 
        elif 1 == n or isequal ( n , 1 ) : return type ( a ) ( a )
        ## 3.  1x1 matrix : exponent can be not only integer 
        elif 1 == a.kRows                         : return pow ( a ( 0 , 0 ) , n )        
        ## 4.  non-integer exponent?
        elif not isinstance ( n , integer_types ) : return NotImplemented
        ## 5. negative integer exponent : first invert and then pow  
        elif n < 0 :
            try : 
                a_inv = LinAlg.M_INVERSE ( a )
            except ValueError : ## matrix cannot be inverted
                return NotImplemented
            if -1 == n : return a_inv 
            return LinAlg.M_POW ( a_inv , abs ( n ) )

        ## 6. regular case: square matrix and non-negative integer exponent 
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
    ## Get minimal element of the matrix
    #  @code
    #  mtrx = ...
    #  mtrx.min_element() 
    #  @endcode
    #  @see Ostap::Math::min_element 
    @staticmethod
    def M_MINELEMENT ( mtrx ) :
        """ Get minimal element of the matrix
        >>> mtrx = ...
        >>> mtrx.min_element() 
        """
        return Ostap.Math.min_element ( mtrx )

    # =============================================================================
    ## Get maximal element of the matrix
    #  @code
    #  mtrx = ...
    #  mtrx.max_element() 
    #  @endcode
    #  @see Ostap::Math::max_element 
    @staticmethod
    def M_MAXELEMENT ( mtrx ) :
        """ Get maximal element of the matrix
        >>> mtrx = ...
        >>> mtrx.max_element() 
        """
        return Ostap.Math.max_element ( mtrx )

    # =============================================================================
    ## Get element with the minimal absolute value of the matrix
    #  @code
    #  mtrx = ...
    #  mtrx.minabs_element() 
    #  @endcode
    #  @see Ostap::Math::minabs_element 
    @staticmethod
    def M_MINABSELEMENT ( mtrx ) :
        """ Get element with the minimal absolute value of the matrix
        >>> mtrx = ...
        >>> mtrx.minabs_element() 
        """
        return Ostap.Math.minabs_element ( mtrx )

    # =============================================================================
    ## Get the element with the maximal absolute value of the matrix
    #  @code
    #  mtrx = ...
    #  mtrx.maxabs_element() 
    #  @endcode
    #  @see Ostap::Math::maxabs_element 
    @staticmethod
    def M_MAXABSELEMENT ( mtrx ) :
        """ Get element with the maximal absoluye value of the matrix
        >>> mtrx = ...
        >>> mtrx.maxabs_element() 
        """
        return Ostap.Math.maxabs_element ( mtrx )
    
    # =============================================================================
    ## Get index of the minimal element of the matrix
    #  @code
    #  mtrx = ...
    #  mtrx.min_element_index () 
    #  @endcode
    #  @see Ostap::Math::ind_min_element 
    @staticmethod
    def M_MININDEX ( mtrx ) :
        """ Get index of the minimal element of the matrix
        >>> mtrx = ...
        >>> mtrx.min_element_index () 
        """
        i , j = Ostap.Math.ind_min_element ( mtrx )
        return i , j
    
    # =============================================================================
    ## Get idnex of maximal element of the matrix
    #  @code
    #  mtrx = ...
    #  mtrx.max_element_index () 
    #  @endcode
    #  @see Ostap::Math::ind_max_element
    @staticmethod
    def M_MAXINDEX ( mtrx ) :
        """ Get maximal element of the matrix
        >>> mtrx = ...
        >>> mtrx.max_element_index () 
        """
        i, j = Ostap.Math.ind_max_element ( mtrx )
        return i, j 

    # =============================================================================
    ## Get index of the element with minial absolute value of the matrix
    #  @code
    #  mtrx = ...
    #  mtrx.minabs_element_index () 
    #  @endcode
    #  @see Ostap::Math::ind_minabs_element 
    @staticmethod
    def M_MINABSINDEX ( mtrx ) :
        """ Get index of the element with minimal absoluet value of the matrix
        >>> mtrx = ...
        >>> mtrx.minabs_element_index () 
        """
        i , j = Ostap.Math.ind_minabs_element ( mtrx )
        return i , j
    
    # =============================================================================
    ## Get index of element with  maximal absolute value of the matrix
    #  @code
    #  mtrx = ...
    #  mtrx.maxabs_element_index () 
    #  @endcode
    #  @see Ostap::Math::ind_maxabs_element
    @staticmethod
    def M_MAXABSINDEX ( mtrx ) :
        """ Get index of element with maximal absolute value of the matrix
        >>> mtrx = ...
        >>> mtrx.maxabs_element_index () 
        """
        i, j = Ostap.Math.ind_maxabs_element ( mtrx )
        return i, j 

    # =============================================================================
    ## Get minimal diagonal element of the matrix
    #  @code
    #  mtrx = ...
    #  mtrx.min_diagnal() 
    #  @endcode
    #  @see Ostap::Math::min_diagonal 
    @staticmethod
    def M_MINDIAGONAL ( mtrx ) :
        """ Get minimal diagonal element of the matrix
        >>> mtrx = ...
        >>> mtrx.min_diagonal () 
        """
        return Ostap.Math.min_diagonal( mtrx )

    # =============================================================================
    ## Get maximal diagonal element of the matrix
    #  @code
    #  mtrx = ...
    #  mtrx.max_diagnal() 
    #  @endcode
    #  @see Ostap::Math::max_diagonal 
    @staticmethod
    def M_MAXDIAGONAL ( mtrx ) :
        """ Get maximal diagonal element of the matrix
        >>> mtrx = ...
        >>> mtrx.max_diagonal () 
        """
        return Ostap.Math.max_diagonal( mtrx )

    # =============================================================================
    ## Get diagonal element with minamal absoluite value  of the matrix
    #  @code
    #  mtrx = ...
    #  mtrx.minabs_diagnal() 
    #  @endcode
    #  @see Ostap::Math::minabs_diagonal 
    @staticmethod
    def M_MINABSDIAGONAL ( mtrx ) :
        """ Get diagonal element with mininal diagoinal value of the matrix
        >>> mtrx = ...
        >>> mtrx.minabs_diagonal () 
        """
        return Ostap.Math.minabs_diagonal( mtrx )

    # =============================================================================
    ## Get diagonal element with maximal absolute value of the matrix
    #  @code
    #  mtrx = ...
    #  mtrx.maxabs_diagnal() 
    #  @endcode
    #  @see Ostap::Math::maxabs_diagonal 
    @staticmethod
    def M_MAXABSDIAGONAL ( mtrx ) :
        """ Get diagonal element with maximal diagonal value of the matrix
        >>> mtrx = ...
        >>> mtrx.maxabs_diagonal () 
        """
        return Ostap.Math.maxabs_diagonal( mtrx )

    # =============================================================================
    ## Are all elements of matrix/vector finite?
    #  @code
    #  mtrx = ...
    #  mtrx.isfinite() 
    #  @endcode
    @staticmethod
    def M_ISFINITE ( mtrx ) :
        """ Are all elements of matrix/vector finite?
        >>> mtrx = ...
        >>> mtrx.isfinite() 
        """
        return Ostap.Math.isfinite ( mtrx )
        
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
    #  for i in vct : print ( i )
    #  @endcode
    @staticmethod
    def V_ITER ( vct ) :
        """ Iterator for SVector
        >>> vct = ...
        >>> for i in vct : print ( i )
        """
        s = len ( vct )
        for i in range ( s ) : yield vct ( i )

    # =============================================================================
    ## iterator for SVector
    #  @code
    #  vct = ...
    #  for i in vct.keys() : print ( i )  
    #  @endcode
    @staticmethod
    def V_KEYS ( vct ) :
        """ Iterator for SVector
        >>> vct = ...
        >>> for i in vct.keys() : print ( i ) 
        """
        s = len ( vct )
        for i in range ( s ) : yield i 
            
    # =============================================================================
    ## iterator for SVector
    #  @code
    #  vct = ...
    #  for i,v in vct.items     () : print ( i , v )  
    #  for i,v in vct.iteritems () : print ( i , v ) ## ditto 
    #  @endcode
    @staticmethod
    def V_ITEMS ( vct ) :
        """ Iterator for SVector
        >>> vct = ...
        >>> for i,v in vct.items    () : print ( i , v )  
        >>> for i,v in vct.iteritems() : print ( i , v ) ## ditto
        """
        s = vct.kSize 
        for i in range ( s ) :
            yield i , vct ( i )
 
    # =============================================================================
    ## Absolute value of a vector (sqrt from sum of squared elements) 
    #  @code
    #  vct = ...
    #  result = abs ( vct ) 
    #  @endcode
    @staticmethod
    def V_ABS ( vct ) :
        """ Absolute value of a vector (sqrt from sum of squared elements) 
        >>> vct  = ...
        >>> result = abs ( vct ) 
        """
        return sum ( v * v for v in vct ) ** 0.5 
   
    # =============================================================================
    ## self-printout of S-vectors
    #  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
    #  @date 2009-09-12
    @staticmethod
    def V_PRETTY ( vct               ,
                   row_labs  = ()    , 
                   fmt       = '%.g' ,
                   title     = ''    ,
                   prefix    = ''    ,
                   style     = ''    , 
                   width     = 6     ,
                   precision = 4     ) :
        """ Self-printout of SVectors: (...)
        >>> vct = ...
        >>> result, expo = vct.pretty_print( ... ) 
        """
        
        N = len ( vct )

        if not title : title = typename ( vct )
        
        ## the maximal element 
        maev = abs ( Ostap.Math.maxabs_element  ( vct ) )
        
        fmtv , expo = fmt_pretty_float ( value     = maev      ,
                                         width     = width     ,
                                         precision = precision )
        
        # =====================================================================
        if expo : # ===========================================================
            scale = 10 ** expo
            title = ( '[%s10^%+d] ' % ( times , expo ) ) + title 
        else    :
            scale = 1
            
        zeros =   fmtv % ( +0 ) , fmtv % ( -0 )

        if row_labs : rows = [ tuple (  v for v in labels ( N , row_labs ) ) ] 
        else        : rows = []

        row   = []
        for v in vct :
            vv = v / scale
            if    0 == v or 0 == vv            : item = ''
            elif iszero ( v ) or iszero ( vv ) : item = ' 0'
            else :
                item = fmtv % vv
                if item in zeros : item = ' 0.0'
            row.append ( item )
        rows.append ( row )

        
        table = T.table  ( rows                     ,
                           alignment       = N*'c'  ,
                           prefix          = prefix ,
                           title           = title  ,
                           style           = style  ,
                           colorize_header = True   )

        ## 
        return table, expo 


    # =============================================================================
    ## self-printout of S-vectors
    #  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
    #  @date 2009-09-12
    @staticmethod
    def P_PRETTY ( vct , fmt = '%.g' , title = '' , prefix = '' , width = 6 , precision = 4 ) :
        """ Self-printout of Permutations : (...)
        >>> vct = ...
        >>> result, expo = vct.pretty_print( ... ) 
        """
        
        N = len ( vct )

        if not title : title = typename ( vct )

        row  = [ ( '%d' % v ) for v in vct ]
            
        table = [ row ]
        table = T.table  ( table , alignment = N*'c' , prefix = prefix , title = title , colorize_header = False )
        ## 
        return table, 0 

    
    # =============================================================================
    ## self-printout of S-vectors-with-errors
    #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    #  @date 2020-05-15
    @staticmethod
    def VE_PRETTY ( vct , fmt = '' , prefix = '' , title = '' , correlations = False , width = 6 , precision = 4 ) :
        """ Self-printout of SVectorWithError: (...)
        >>> vcte = ...
        >>> result, expo = vcte.pretty_print ( ... ) 
        """
        
        values = vct.value      ()    
        cov2   = vct.covariance ()

        if not title : title = typename ( vct )
        
        N    = len ( values )
        cols = N
        rows = N
                
        ## the maximal element 
        maev = abs ( Ostap.Math.maxabs_element  ( values ) )

        if correlations :
            maec   = abs ( Ostap.Math.maxabs_diagonal ( cov2   ) )
            perror = math.sqrt ( maec )
        else :
            maec   = abs ( Ostap.Math.maxabs_element  ( cov2   ) )
            perror = maec 
            
        _ , fmtv , fmte , expo = fmt_pretty_error ( maev      ,
                                                    perror    , 
                                                    width     = width     ,
                                                    precision = precision )

        if expo :
            scale = 10 ** expo
            title = ( '[%s10^%+d] ' %( times , expo ) ) + title 
        else    : scale = 1

        if correlations :
            title  = title + ' correlations' 
            mtrx   = cov2.correlations ()
            mscale = 0.01
            fmtm   = '%+.1f%%'
            fmte   = '+/-' + fmte 
        else            :
            mtrx   = cov2
            mscale = scale
            fmtm   = fmtv

        zeros = [] 
        for f in ( fmtv , fmte , fmtm ) :
            zeros.append ( f % ( +0 ) )
            zeros.append ( f % ( -0 ) )
        zeros = tuple ( zeros )
        
        ## table = [ tuple ( [ '\\' ] + [ '%d' % i for i in range ( cols ) ] ) ] 
        table = [ tuple ( [ '\\' ] + [ l for l in labels ( cols ) ] ) ] 

        ## 1st row : values
        row    = [ infostr ( 'V' )  ] 
        for v in values :
            vv = v / scale 
            if   0 == v or 0 == vv             : item = ''
            elif iszero ( v ) or iszero ( vv ) : item = ' 0'
            else             :
                item = fmtv % vv
                if item in zeros : item = ' 0.0' 
            row.append ( item )    
        table.append ( row )

        ## for `correlation`: the 2nd row: errors 
        if correlations :
            row    = [ infostr ( 'e' )  ] 
            for i in range  ( N ) :
                v  = cov2 ( i , i )
                if   iszero ( v )  : item = '+/-0'
                elif v < 0         : item = '???'
                else :
                    vv = math.sqrt ( v ) / scale 
                    if  iszero ( vv ) :  item = '+/-0'
                    else : 
                        item = fmte % vv
                        if item in zeros : item = '+/-0' 
                row.append ( item )
            table.append ( row )
                
        ## print the matrix (covariance or correlation )
        for i in range ( rows ) :
            row = [ infostr ( '%d' % i ) ]
            for j in range ( cols ) : 
                if j < i :  item = ditto  
                else     :
                    v  = mtrx ( i , j )
                    vv = v / mscale
                    if   iszero ( v ) or iszero ( vv ) : item = ' 0'
                    else:
                        item = fmtm % vv
                        if item in zeros : item = ' 0.0'
                if item and i == j : item = diag ( item )
                row.append ( item ) 
            table.append ( row ) 

        table = T.table  ( table , alignment = 'r'+cols*'c' , prefix = prefix , title = title ) 
        return table , expo 
    
    # =============================================================================
    ## self-printout of S-vectors
    #  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
    #  @date 2009-09-12
    @staticmethod
    def V_STR ( vct , fmt = '%.g' , title = '' , prefix = '' , width = 6 , precision = 4 ) :
        """ Self-printout of SVectors: (...)
        >>> vct = ...
        >>> print ( vct ) 
        """
        N = len ( vct ) 
        if 15 < N  : return '[ ' + ( ', '.join ( fmt % v for v in vct ) ) + ' ]'

        result, _ = vct.pretty_print ( title     = title     ,
                                       prefix    = prefix    ,
                                       width     = width     ,
                                       precision = precision )
        return result

    # =============================================================================
    ## Convert vector into 1D-histogram
    #  @code
    #  vct = ...
    #  h   = vct.th1() 
    #  @endcode 
    @staticmethod
    def V_TH1 ( vct ) :
        """ Convert vector into 1D-histogram
        >>> vct = ...
        >>> h   = vct.th1() 
        """
        if hasattr ( vct , '_th1' ) : return  vct._th1 
        N     = len ( vct )
        xlow  = 0 - 0.5
        xhigh = N - 0.5
        th1   = ROOT.TH1F ( hID() , '1D-Histogram from the vector' , N , xlow , xhigh )
        for i , v in eumerate ( vct ) : th1.SetBinContent( i + 1 , v ) 
        self._th1 = th1
        return th1

    # =============================================================================
    ## Draw the vector via conversion into into 1D-histogram
    #  @code
    #  vct = ...
    #  vct.draw() 
    #  @endcode 
    @staticmethod
    def V_DRAW ( vct , opts = '' , *args , **kwargs ) :
        """ Draw the vector via conversion into 1D-histogram
        >>> vct = ...
        >>> vct.draw()  
        """
        h1 = vct.th1()
        h1.draw ( opts , *args , **kwargs )
        return vct 
    
    # =============================================================================
    ## self-printout of permutations 
    #  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
    #  @date 2009-09-12
    @staticmethod
    def P_STR ( vct , fmt = '%.g' , title = '' , prefix = '' , width = 6 , precision = 4 ) :
        """ Self-printout of permutatons (...)
        >>> vct = ...
        >>> print ( vct ) 
        """
        N = len ( vct ) 
        if 15 < N  : return '[ ' + ( ', '.join ( fmt % v for v in vct ) ) + ' ]'

        result, _ = vct.pretty_print ( title     = title     ,
                                       prefix    = prefix    ,
                                       width     = width     ,
                                       precision = precision )
        return result
            
    # =============================================================================
    ## self-printout of S-vectors-with-errors
    #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    #  @date 2020-05-15
    @staticmethod
    def VE_STR ( vct , fmt = '' , prefix = '' , title = '' , correlations = False , width = 6 , precision = 4 ) :
        """ Self-printout of SVectorWithError: (...)
        >>> vct = ...
        >>> print ( vct ) 
        """

        result, _ = vct.pretty_print ( title        = title        ,
                                       prefix       = prefix       , 
                                       correlations = correlations ,
                                       width        = width        ,
                                       precision    = precision    )
        return result 

    # =============================================================================
    ## Transform vector-with-errors to another variable
    #  for \f[ y = y ( x ) , C(y) = J C(x) J^T, \f]
    #  where  \f$ J = \left( \frac{\partial y}{\partial x }\right) \f$.
    #
    #  Transform to (r,phi)-variables 
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
        """ Transform vector-with-errors to another variable
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
        """ Convert vector to plain array:
        >>> vct = ...
        >>> na  = vct.to_array() 
        """
        return array.array( 'd', vct  )


    # =============================================================================
    ## iterator for SMatrix
    #  @code
    #  matrix = ...
    #  for i in matrix : print ( i ) 
    #  @endcode
    @staticmethod 
    def M_ITER ( mtrx ) :
        """ Iterator for SMatrix
        >>> matrix = ...
        >>> for i in matrix : print ( i )
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
    #  for i in matrix.keys() : print ( i )
    #  @endcode
    @staticmethod 
    def M_KEYS ( mtrx ) :
        """ Iterator for SMatrix
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
    #  for ij,v in matrix.items     () : print ( ij,v ) 
    #  for ij,v in matrix.iteritems () : print ( ij,v ) ## ditto
    #  @endcode
    @staticmethod     
    def M_ITEMS  ( mtrx ) :
        """ Iterator for SMatrix
        >>> matrix = ...
        >>> for ij,v in matrix.items    () : print ( ij,v ) 
        >>> for ij,v in matrix.iteritems() : print ( ij,v ) ## ditto
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
    #   result, expo = matrix.pretty_print ( ... ) 
    #   @endcode
    @staticmethod 
    def M_PRETTY ( mtrx           ,
                   row_labs  = () , ## labels for rows 
                   col_labs  = () , ## labels for rcolumns 
                   fmt       = '' ,
                   prefix    = '' ,
                   title     = '' ,
                   style     = '' , 
                   width     =  6 ,
                   precision =  4 ) :
        """ Self-printout of matrices
        >>> matrix = ...
        >>> result , expo = matrix.pretty_print ( ... )
        """
        rows = mtrx.kRows
        cols = mtrx.kCols

        if not title : title = typename ( mtrx )

        mae = abs ( Ostap.Math.maxabs_element ( mtrx ) )
        fmtv , expo = fmt_pretty_float ( mae , width = width , precision = precision )

        zeros = fmtv % ( +0.0 ) , fmtv % ( -0.0 )        
        if expo :
            scale  = 10 ** expo
            title  = ( '[%s10^%+d] ' % ( times , expo ) ) + title 
        else    :
            scale = 1 
            
        ## table = [ tuple ( [ '\\' ] + [ '%d' % i for i in range ( cols ) ] ) ]
        table = [ tuple ( [ '\\' ] + [ l for l in labels ( cols , col_labs ) ] ) ]

        for i,l in enumerate ( labels ( rows , row_labs ) ) :
            row = [ infostr ( l )  ]
            for j in range ( cols ) :                
                v     = mtrx ( i , j ) 
                value = v / scale 
                if   0 == v or 0 == value             : item = ''
                elif iszero ( v ) or iszero ( value ) : item = '0'
                else :
                    item = fmtv % value
                    if item in zeros : item = ' 0.0'
                if item and i == j : item = diag ( item )
                row.append ( item ) 
            table.append ( row )
            
        table = T.table  ( table , title = title , alignment = 'r'+cols*'c' , prefix = prefix , style = style )
        return table, expo 

    # =============================================================================
    ## Get the L(p,q) norm of matrix
    #  - \f$ 1< p \f$ 
    #  - \f$ 1< q \f$ 
    #  @code
    #  matrix = ..
    #  norm   = matrix.lnorm () 
    #  @endcode
    #  @see https://en.wikipedia.org/wiki/Matrix_norm
    @staticmethod
    def M_LNORM ( mtrx , p = 2 , q = 2 ) :
        """ Get the L(p,q) norm of matrix (1<p,q) 
        - see https://en.wikipedia.org/wiki/Matrix_norm
        >>> matrix = ..
        >>> norm  = matrix.lnorm ( p , q ) 
        """
        assert isinstance ( p , num_types ) and 1 <= p , "Invalid p : %s" % p
        assert isinstance ( q , num_types ) and 1 <= q , "Invalid q : %s" % q

        result = 0
        rows   = mtrx.kRows
        cols   = mtrx.kCols

        ## Frobenius norm 
        if 2 == p and 2 == q :
            for j in range ( cols ) :
                for i in range ( rows ) :
                    value   = mtrx ( i , j ) 
                    result += value * value
            return math.sqrt ( result )

        qop  = float ( q ) / float ( p )        
        nd   = min ( rows , cols )
        for j in range ( cols ) :
            row = 0 
            for i in range ( rows ) :
                v    = mtrx ( i , j ) 
                row += pow ( abs ( v ) , p )
            result  += pow ( row , qop )
            
        result = pow ( result , 1. / q )
        return result

    # =============================================================================
    ## Get the max-norm of matrix: \f$ max_{ij} |a_{ij}| \f$ 
    #  @code
    #  matrix = ..
    #  norm   = matrix.mnorm () 
    #  @endcode
    #  @see https://en.wikipedia.org/wiki/Matrix_norm
    @staticmethod
    def M_MNORM ( mtrx , submult = False ) :
        """ Get the max-norm of the matrix: 
        - maximal absolute value of all elements 
        >>> matrix = ..
        >>> norm  = matrix.mnorm ( ) 
        """
        
        maev = abs ( Ostap.Math.maxabs_element  ( mtrx ) )
        
        if submult :            
            rows  = mtrx.kRows
            cols  = mtrx.kCols
            maev *= math.sqrt ( rows * cols )
            
        return maev
    
    # =============================================================================
    ## Get the matrix the same type with diagonal elements, other elements are zeros
    #  @code
    #  matrix = ..
    #  diag   = matrix.diagonal() 
    #  @endcode 
    @staticmethod
    def M_DIAGONAL ( mtrx ) :
        """ Get the matrix the same typewith diagonal elements, other elements are zeros
        >>> matrix = ..
        >>> diag   = matrix.diagonal() 
        """
        newm = type ( mtrx ) () ## make new matrix of the same type 
        rows = mtrx.kRows
        cols = mtrx.kCols
        nd   = min ( rows , cols )
        ## copy diagonal elements 
        for  i in range ( nd ) : newm [ i , i ] = mtrx ( i , i )
        return newm 

    # ============================================================================
    # transpose the matrix 
    @staticmethod
    def S_T ( mtrx ) :
        """ Transpose the matrix (make a copy)
        """
        nr, nc  = mtrx.kRows , mtrx.kCols 
        return Ostap.Math.Matrix( nc , nr ) ( ROOT.Math.Transpose ( mtrx ) ) 

    # ============================================================================
    # trabnspose the matrix 
    @staticmethod
    def MS_T ( mtrx ) :
        """ Transpose the matrix (makea copy) 
        """
        nr, nc  = mtrx.kRows , mtrx.kCols 
        return Ostap.Math.SymMatrix( nc ) ( mtrx ) 
    
    # =============================================================================
    ##  Self-printout of symmetric matrices
    #   @code  
    #   matrix = ...
    #   result, exp = matrix.pretty_print ( ... ) 
    #   @endcode
    @staticmethod 
    def MS_PRETTY ( mtrx           ,
                    fmt       = '' ,
                    row_labs  = () , ## labels for rows 
                    col_labs  = () , ## labels for rcolumns 
                    prefix    = '' ,
                    title     = '' ,
                    style     = '' , 
                    width     = 6  ,
                    precision = 4  ) :
        """ Self-printout of symmetric matrices
        >>> matrix = ...
        >>> result, expo = matrix.pretty_print ( ... ) 
        """
        rows = mtrx.kRows
        cols = mtrx.kCols
        
        if not title : title = typename ( mtrx )

        mae = abs ( Ostap.Math.maxabs_element ( mtrx ) )
        fmtv , expo = fmt_pretty_float ( mae , width = width , precision = precision )

        zeros = fmtv % ( +0.0 ) , fmtv % ( -0.0 )        
        if expo :
            scale  = 10 ** expo
            title  = ( '[%s10^%+d] ' % ( times , expo ) ) + title 
        else    :
            scale = 1 
            
        ## table = [ tuple ( [ '\\' ] + [ '%d' % i for i in range ( cols ) ] ) ]
        table = [ tuple ( [ '\\' ] + [ l for l in labels ( cols , col_labs ) ] ) ]
        ## for i,l in range ( rows ) :
        for i,l in enumerate ( labels ( rows , row_labs ) ) :
            row = [ infostr ( l ) ]
            for j in range ( cols ) : 
                if j < i : item = ditto 
                else     :
                    v     = mtrx ( i , j ) 
                    value = v / scale
                    if   0 == v or 0 == value             : item = ''
                    elif iszero ( v ) or iszero ( value ) : item = ' 0'
                    else :
                        item = fmtv % value
                        if item in zeros : item = ' 0.0'
                if item and i == j : item = diag ( item ) 
                row.append ( item )                     
            table.append ( row )
            
        table = T.table  ( table , alignment = 'r'+cols*'c' , prefix = prefix , title = title , style = style )
        return table, expo 
    
    # =========================================================================
    ##  Self-printout of matrices
    #   @code  
    #   matrix = ...
    #   print matrix
    #   @endcode
    @staticmethod 
    def M_STR ( mtrx , fmt = '' , prefix = '' , title = '' , width = 6 , precision = 4 ) :
        """ Self-printout of matrices
        >>> matrix = ...
        >>> print (matrix)
        """
        result , _ = mtrx.pretty_print ( prefix    = prefix    ,
                                         title     = title     ,
                                         width     = width     ,
                                         precision = precision )
        return result

    # =============================================================================
    ## Convert matrix into 2D-histogram
    #  @code
    #  mtrx = ...
    #  h    = vct.th2() 
    #  @endcode 
    @staticmethod
    def M_TH2 ( mtrx ) :
        """ Convert matrix into 2D-histogram
        >>> mtrx = ...
        >>> h   = mtrx.th2() 
        """
        if hasattr ( mtrx , '_th2' ) : return  mtrx._th2
        
        NX    = mtrx.kCols
        NY    = mtrx.kRows
        
        xlow  = 0  - 0.5
        xhigh = NX - 0.5
        
        ylow  = 0  - 0.5
        yhigh = NY - 0.5
        
        th2   = ROOT.TH2F ( hID() ,
                            '2D-Histogram from the matrix' ,
                            NX , xlow , xhigh , 
                            NY , ylow , yhigh )
    
        for i in range ( NX ) :
            for j in range ( NY ) : 
                th2.SetBinContent( i + 1 , j + 1 , mtrx ( i , NY - j - 1 ) ) 

        mtrx._th2 = th2
        return th2
    
    # =============================================================================
    ## Draw the matrix via conversion into into 2D-histogram
    #  @code
    #  mtrx = ...
    #  mtrx.draw() 
    #  @endcode 
    @staticmethod
    def M_DRAW ( mtrx , opts = '' , *args , **kwargs ) :
        """ Draw the matrix via conversion into 2D-histogram
        >>> mtrx = ...
        >>> mtrx.draw()  
        """
        h2 = mtrx.th2()
        h2.draw ( opts , *args , **kwargs )
        return mtrx 
 
    # =============================================================================
    ##  Self-printout of symmetric matrices
    #   @code  
    #   matrix = ...
    #   print matrix
    #   @endcode
    @staticmethod 
    def MS_STR ( mtrx , fmt = '' , prefix = '' , title = '' , width = 6 , precision = 4 ) :
        """ Self-printout of symmetric matrices
        >>> matrix = ...
        >>> print(matrix)
        """
        result , _ = mtrx.pretty_print ( prefix    = prefix    ,
                                         title     = title     ,
                                         width     = width     ,
                                         precision = precision )
        return result
    
    # =========================================================================
    ## get the correlation matrix
    @staticmethod 
    def MS_CORR ( mtrx ) :
        """ Get the correlation matrix
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
                ok1 = False 
                for j in range ( i , rows ) : _c [ i , j ] = NaN 
                continue
        
            sii = sqrt ( ii )
            _c[ i , i ] = 1
            
            for j in range ( i + 1 , rows ) :            
                jj  = mtrx ( j , j ) 
                sjj = sqrt ( jj    )
                if 0 > sjj or iszero ( sjj ) :
                    ok1 = False
                    _c [ i , j ] = NaN
                    _c [ j , i ] = NaN 
                else : 
                    ij  = mtrx ( i , j )
                    eij = ij / ( sii * sjj )
                    if  1 < abs ( eij ) : ok2 = False  
            
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
        """ Get row from the matrix
        >>> mtrx = ...
        >>> row2 = mtrx.row ( 2 ) 
        """
        ## allow slightly negative indices 
        if index < 0 : index += mtrx.kRows

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
        """ Get column from the matrix
        >>> mtrx = ...
        >>> c2 = mtrx.column ( 2 ) 
        """
        
        ## allow slightly negative indices 
        if index < 0 :index += mtrx.kCols
        
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
        """ Iterator over rows
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
        """ Iterator over columns
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
        """ Get the eigenvalues for symmetric matrices :
        >>> mtrx = ...
        >>> values = mtrx.eigenValues ( sorted = True )
        """
        ## 
        operation  = LinAlg.method_EIGEN ( mtrx  )
        if not operation  : raise NotImplementedError ("EigenValues: not implemented for %s" % typename ( mtrx ) )
        ## 
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
    def MS_EIGENVECTORS ( mtrx , sorted = True , ascending = True ) :
        """ Get the eigenvalues for symmetric matrices :
        >>> mtrx = ...
        >>> values, vectors = mtrx.eigenVectors ( sorted = True )
        >>> vectors = [ vectors.column{i) for i in range ( mtrx.rCols ) ] 
        """
        ## 
        operation = LinAlg.method_EIGEN ( mtrx )
        if not operation :
            raise NotImplementedError ("EigenVectors: not implemented for %s" % typename ( mtrx ) )
        ## 
        krows   = mtrx.kRows
        kcols   = mtrx.kCols        
        values  = LinAlg.Vector ( krows        ) ()
        vectors = LinAlg.Matrix ( krows, kcols ) () 
        st      = operation ( mtrx  , values , vectors , sorted , ascending )
        assert st.isSuccess () , "Eigen-vectors: status code %s" % st
        ## 
        return values, vectors  

    # =========================================================================
    ## Iterator over (eigenvalue/eigenvector) pairs for symmetric matrix
    #  @code
    #  matrix = ...
    #  for e, v in matrix.eigenitems ( sorted  = True ) :
    #      ... 
    #  @endcode
    @staticmethod
    def MS_EIGENITEMS ( mtrx , sorted = True , ascending = True ) :
        """ Iterator over (eigenvalue,eigenvector) pairs  for symmetric matrices :
        >>> mtrx = ...
        >>> for e, v in mtrx.eigenitems ( sorted = True ) : 
        >>>      ...
        """
        ## 
        operation = LinAlg.method_EIGEN ( mtrx )
        if not operation :
            raise NotImplementedError ("EigenVectors: not implemented for %s" % typename ( mtrx ) )
        ## 
        ## eigenvalues & eigenvectors
        values, vectors = mtrx.eigenVectors ( sorted = sorted , ascending = ascending ) 
        for i in range ( vectors.kCols ) :
            yield values [ i ] , vectors.column ( i )

    # =========================================================================
    ## reduce SVector
    @staticmethod
    def V_REDUCE ( vct ) :
        """ Reduce SVector"""
        return svct_factory, ( array.array ( 'd' , vct ) , )

    # =========================================================================
    ## reduce SMatrix 
    @staticmethod
    def M_REDUCE ( mtrx ) :
        """ Reduce SMatrix
        """
        NR   = mtrx.rep_size
        a    = mtrx.Array() 
        data = array.array ( 'd' , ( a[i] for i in range ( NR ) ) ) 
        return smtrx_factory, ( mtrx.kRows, mtrx.kCols , data ) 

    # =========================================================================
    ## reduce symmetric SMatrix 
    @staticmethod
    def MS_REDUCE ( mtrx ) :
        """ Reduce symmetric SMatrix
        """
        NR   = mtrx.rep_size 
        a    = mtrx.Array() 
        data = array.array ( 'd' , ( a[i] for i in range ( NR ) ) ) 
        return symmm_factory, ( mtrx.kRows , data ) 

    # =========================================================================
    ## reduce SVectorWithErrors 
    @staticmethod
    def VE_REDUCE ( vct ) :
        """ Reduce SVectorWithErrors"""
        return svcte_factory , ( vct.value() , vct.cov2() ) 

    # =========================================================================
    @staticmethod
    def VE_RANDOM ( vct , N , use_numpy = True ) :
        """ Generate random numbers from the vector
        >>> vct = ...
        >>> for x in vct.random ( 1000 ) :
        >>> ... 
        """
        
        assert isinstance ( N , integer_types ) and 0 <= N , 'Invalid number of shoots!'
        
        v  = vct.value ()
        c2 = vct.cov2  ()
        n  = len ( v )
        
        L  = Ostap.Matrix ( n , n , vct._scalar_ )() 
        ok = Ostap.Math.cholesky ( vct , L )
        assert ok , 'random: Cholesky decomposition failed!'

        ## use numpy if/when available
        if numpy and use_numpy :
            
            v = v.to_numpy() 
            l = L.to_numpy()
            for i in range  ( N ) : 
                u = numpy.random.normal ( loc = 0 , scale = 1, size = n) 
                yield v + numpy.dot ( l , u )
                
        else :
            
            x = Ostap.Vector ( n , v._scalar_ )()
            for i in range ( N ) :                
                for i in range ( n ) : x [ i ] = random.gauss ( 0 , 1 )                
                yield v + L * x
                
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
    ## Convert SMatrix/TMatrix objects into GSL MAtrix
    @staticmethod 
    def M_2GSL ( mtrx ) :
        """ Convert SMatrix/TMatrix objects into GSL Matrix
        """
        return Ostap.GSL.matrix ( mtrx )
    
    # =========================================================================
    ## Convert SVector/TVector objects into GSL Vector
    @staticmethod 
    def V_2GSL ( vct ) :
        """ Convert SVecto/TVector objects into GSL Vector
        """
        return Ostap.GSL.vector ( vct )
    
    # =========================================================================
    ## (P)LU decomposition
    #  @see Ostap::GSL::PLU
    #  @code
    #  matrix = ...
    #  P , L , U = matrix.PLU() 
    #  @endcode    
    @staticmethod
    def S_PLU  ( mtrx ) :
        """ Perform (P)LU decomposition of the matrix 
        >>> matrix = ...\
        >>> P, L, U = matarix.PLU() 
        - see `Ostap.GSL.PLU` 
        """
        ## convert to GLS 
        A = mtrx.to_GSL()
        ## mape (P)LU decomposiiton 
        P, L, U = A.PLU ()
        ## convert BACK:
        P = Ostap.GSL.Matrix ( P )
        ## 
        P = P.to_SMatrix()
        L = L.to_SMatrix()
        U = U.to_SMatrix()
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
    def S_PQR  ( mtrx ) :
        """ Perform (P)QR decomposition of the matrix 
        >>> matrix = ...\
        >>> P, Q, R = matarix.PQR() 
        - see `Ostap.GSL.PQR` 
        """
        ## convert to GLS 
        A = mtrx.to_GSL()
        ## mape (P)LU decomposiiton 
        P, Q, R = A.PQR ()
        ## convert BACK:
        P = Ostap.GSL.Matrix ( P ).transpose() 
        ## 
        P = P.to_SMatrix()
        Q = Q.to_SMatrix()
        R = R.to_SMatrix()
        ## 
        return P , Q , R 

    # ==========================================================================
    ## Perform LQ decomposition of the matrix
    @staticmethod
    def S_LQ ( mtrx ) :
        """ Perfrom LQ decompositionof matrix A : A = LQ 
        - A is input MxN matrix 
        - L is lower trapezoidal  MxN matrix
        - Q is orthogonal NxN matrix    
        >>> A = ...
        >>> L, Q = A.LQ() 
        """
        ## convert matrix to GSL 
        A = mtrx.to_GSL()
        ## make decomposition
        L , Q = A.LQ()
        ## convert back
        L = L.to_SMatrix()
        Q = Q.to_SMatrix()
        return L, Q
    
    # ==========================================================================
    ## Perform QL decomposition of the matrix
    @staticmethod
    def S_QL ( mtrx ) :
        """ Perfrom LQ decompositionof matrix A : A = QL 
        - A is input MxN matrix 
        - Q is orthogonal NxN matrix    
        - L is lower trapezoidal  MxN matrix
        >>> A = ...
        >>> Q , L = A.QL () 
        """
        ## convert matrix to GSL 
        A = mtrx.to_GSL()
        ## make decomposition
        Q , L = A.QL ()
        ## convert back
        Q = Q.to_SMatrix()
        L = L.to_SMatrix()
        return Q , L 

    # =========================================================================
    ## COD: Complete orthogonal decomposition AP = Q R Z^T
    @staticmethod
    def S_COD ( mtrx ) :
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
        ## convert matrix to GSL 
        A = mtrx.to_GSL()
        ## make decomposition
        P , Q , R , Z = A.COD()
        ## convert BACK:
        P = Ostap.GSL.Matrix ( P ).transpose() 
        ## 
        P = P.to_SMatrix()
        Q = Q.to_SMatrix()
        R = R.to_SMatrix()
        Z = Z.to_SMatrix()
        return P, Q, R , Z 
    
    # ===============================================================================
    ## SVD : singular Value Decomposition  \f$ A = U S V^T\f$
    def S_SVD ( mtrx , golub = True ) : 
        """ SVD : singular Value Decomposition  A = U S V^T 
        - A input MxN matrix 
        - K = min ( M , N ) : 
        - U MxK orthogonal matrix 
        - S KxK Diagonal matrix of singular values 
        - V NxK orthogonal matrix 
        >>> A = ...
        >>> S , U , V = A.SVD() 
        """
        ## convert matrix to GSL 
        A = mtrx.to_GSL()
        ## make decomposition
        S , U , V = A.SVD() 
        ## convert BACK:
        N = S.size() 
        M = Ostap.Math.SymMatrix ( N )()
        for i in range ( N ) : M [ i , i ] = S ( i )
        S = M 
        U = U.to_SMatrix()
        V = V.to_SMatrix()
        return S , U , V 

    # ===============================================================================
    ## SCHUR : Schur Decomposition  \f$ A = Z T Z^T \f$ 
    def S_SCHUR ( mtrx ) :
        """ Schur  decomposition of the square matrix A: A = Z T Z^T
        - Z is orthogonal 
        - T is Shur's form 
        >>> A = ...,
        >>> Z , T = A.SCHUR () 
        """
        ## convert matrix to GSL 
        A = mtrx.to_GSL()
        ## make decomposition
        Z , T = A.SCHUR() 
        ## convert BACK:
        Z = Z.to_SMatrix()
        T = T.to_SMatrix()
        return Z , T  

    # ===============================================================================
    ## POLAR : Polar Decomposition  \f$ A = UP \f$ 
    def S_POLAR ( mtrx ) :
        """ Polar decomposition of the square matrix A: A = UP
        - U is orthogonal 
        - P is positive semi-definitive 
        >>> A = ...,
        >>> U , P = A.POLAR() 
        """
        ## convert matrix to GSL 
        A = mtrx.to_GSL()
        ## make decomposition
        U , P = A.POLAR() 
        ## convert BACK:
        P = P.to_SMatrix()
        U = U.to_SMatrix()
        return U , P 

    
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

        t.__abs__       = LinAlg.V_ABS
        
        t. __str__      = LinAlg.V_STR
        t. __repr__     = LinAlg.V_STR
        t. table        = LinAlg.V_STR        
        ## pretty printout  
        t.pretty_print  = LinAlg.V_PRETTY

        ## convert to TH1
        t. th1          = LinAlg.V_TH1
        ## Draw it via convertsion to TH!
        t. draw         = LinAlg.V_DRAW
        
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

        t.isfinite      = LinAlg.M_ISFINITE

        t.to_gsl        = LinAlg.V_2GSL 
        t.as_gsl        = LinAlg.V_2GSL 
        t.to_GSL        = LinAlg.V_2GSL 
        t.as_GSL        = LinAlg.V_2GSL 
        
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
        """ Decorate SMatrix
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
        m.table         = LinAlg.M_STR
        m.pretty_print  = LinAlg.M_PRETTY

        ## cnopversion to TH2 
        m. th2          = LinAlg.M_TH2
        ## draw ti vi conversion to TH2 
        m. draw         = LinAlg.M_DRAW 

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

        m.t             = LinAlg.S_T   
        m.transpose     = LinAlg.S_T   
        
        m.tmatrix       = LinAlg.M_TM 

        m.min_element          = LinAlg.M_MINELEMENT
        m.max_element          = LinAlg.M_MAXELEMENT
        m.minabs_element       = LinAlg.M_MINABSELEMENT
        m.maxabs_element       = LinAlg.M_MAXABSELEMENT
        
        m.min_diagonal         = LinAlg.M_MINDIAGONAL
        m.max_diagonal         = LinAlg.M_MAXDIAGONAL
        m.minabs_diagonal      = LinAlg.M_MINABSDIAGONAL
        m.maxabs_diagonal      = LinAlg.M_MAXABSDIAGONAL
        
        m.min_element_index    = LinAlg.M_MININDEX
        m.max_element_index    = LinAlg.M_MAXINDEX
        m.minabs_element_index = LinAlg.M_MINABSINDEX
        m.maxabs_element_index = LinAlg.M_MAXABSINDEX 

        m.diagonal             = LinAlg.M_DIAGONAL
        m.lnorm                = LinAlg.M_LNORM 
        m.mnorm                = LinAlg.M_MNORM 

        m.isfinite             = LinAlg.M_ISFINITE

        if m.kRows == m.kCols :
            m.inverse   = LinAlg.M_INVERSE            

        ## conversion to float 
        if 1 == m.kRows and 1 == m.kCols :
            m.__float__ = lambda s : s ( 0 , 0 )

        ## should be class property!! 
        m.shape         = property ( LinAlg.M_SHAPE , None , None )

        if m.kRows == m.kCols :
            m.__pow__       = LinAlg.M_POW 
            m.sym           = LinAlg.M_SYM
            m.asym          = LinAlg.M_ASYM 
            m.skew          = LinAlg.M_ASYM 
            
        m.__reduce__    = LinAlg.M_REDUCE 
        m.rep_size      = classgetter ( lambda cls : cls.rep_type.kSize ) 

        m.to_gsl        = LinAlg.M_2GSL 
        m.as_gsl        = LinAlg.M_2GSL 
        m.to_GSL        = LinAlg.M_2GSL 
        m.as_GSL        = LinAlg.M_2GSL 

        m.PLU           = LinAlg.S_PLU 
        m.PQR           = LinAlg.S_PQR
        m.LQ            = LinAlg.S_LQ 
        m.COD           = LinAlg.S_COD
        m.SVD           = LinAlg.S_SVD

        if  ( 2  , 7 ) <= gsl_info :
            m.QL        = LinAlg.S_QL
            
        if m.kRows == m.kCols :
            m.SCHUR     = LinAlg.S_SCHUR
            m.POLAR     = LinAlg.S_POLAR
                        
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
        """ Decorate the symmetrix  SMatrix
        """

        if m in LinAlg.decorated_matrices : return m 

        LinAlg.deco_matrix ( m )
        
        R = m.rep_type
        
        m.__str__        = LinAlg.MS_STR
        m.__repr__       = LinAlg.MS_STR
        m.table          = LinAlg.MS_STR
        m.pretty_print   = LinAlg.MS_PRETTY 

        m.correlations   = LinAlg.MS_CORR
        
        m.Sim            = LinAlg.SIM
        m.sim            = LinAlg.SIM
        m.SimT           = LinAlg.SIMT
        m.simT           = LinAlg.SIMT
        
        m.t               = LinAlg.MS_T   
        m.transpose       = LinAlg.MS_T   
        
        m.eigenValues    = LinAlg.MS_EIGENVALUES 
        m.eigenVectors   = LinAlg.MS_EIGENVECTORS
        m.eigen_values   = LinAlg.MS_EIGENVALUES 
        m.eigen_vectors  = LinAlg.MS_EIGENVECTORS
        m.eigen_items    = LinAlg.MS_EIGENITEMS 
        m.eigenitems     = LinAlg.MS_EIGENITEMS 
        m.eigenItems     = LinAlg.MS_EIGENITEMS 
        m.eigenpairs     = LinAlg.MS_EIGENITEMS 
        
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
        t. table        = LinAlg.VE_STR
        t. pretty_print = LinAlg.VE_PRETTY 

        t. __len__      = lambda s : s.kSize 
        t. __contains__ = lambda s, i : 0 <= i < s.kSize

        t. __iter__     = LinAlg.V_ITER      
        t.transform     = LinAlg.VE_TRANSFORM

        t.__reduce__    = LinAlg.VE_REDUCE 

        t.random        = LinAlg.VE_RANDOM

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
            v._scalar_ = t
            
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
            m._scalar_ = t
            
            LinAlg.known_smatrices [ tt ] = m
            LinAlg.deco_matrix     ( m )

            ## check k-dimension 
            tt1 = k , t
            if not tt1 in LinAlg.known_svectors     : LinAlg.   Vector ( k , t )
            if not tt1 in LinAlg.known_ssymmatrices : LinAlg.SymMatrix ( k , t )    

            ## check N-dimension
            tt2 = n , t
            if not tt2 in LinAlg.known_svectors     : LinAlg.   Vector ( n , t )
            if not tt2 in LinAlg.known_ssymmatrices : LinAlg.SymMatrix ( n , t )    

            ## check transposed matrix
            if n != k : 
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
            m._scalar_ = t

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
            v._scalar_ = t
                        
            LinAlg.known_svectorse [ tt ] = v
            LinAlg.deco_vectore    ( v )
            ##
            if not tt  in LinAlg.known_svectors     : LinAlg.   Vector ( n ,     t )    
            if not tt  in LinAlg.known_ssymmatrices : LinAlg.SymMatrix ( n ,     t )
            
            tt1 = n , n , t 
            if not tt1 in LinAlg.known_smatrices    : LinAlg.   Matrix ( n , n , t )

            for i in range ( 1 , n ) :
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


# =============================================================================
if numpy : # ==================================================================
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
        assert isinstance ( a , numpy.ndarray ) and 2 == len ( a.shape ) ,\
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
        assert isinstance ( a , numpy.ndarray ) and 1 == len ( a.shape ) ,\
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
        """ Convert 1D or 2D numpy array to SMatrix/SVecrtor
        >>> a = ...
        >>> s = toSObject ( a ) 
        """
        assert isinstance ( a , numpy.ndarray ) and 1 <= len ( a.shape ) <= 2 ,\
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

# =============================================================================
if numpy : # ==================================================================
    # =========================================================================
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
    """ Check what LinAlg operations are defined for these two objects    
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
        result = '%s' % typename ( operation ( a , b ) ) if operation and ok else '' 
        
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
    logger.info ( "%s for '%s' and '%s':\n%s" % ( title , typename ( a ) , typename ( b ) , table ) ) 
                      
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
