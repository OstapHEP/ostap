#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/base.py
#
#  Simple file to provide "easy" access in python for the basic ROOT::Math classes
#  @see Ostap/Point3DTypes.h
#  @see Ostap/Vector3DTypes.h
#  @see Ostap/Vector4DTypes.h
#  @see Ostap/GenericVectorTypes.h
#
#  The usage is fairly trivial:
#
#  @code
#
#  import ostap.math.base
#
#  @endcode
#
#  Important: All types are defined in corresponding C++ namespaces
#
#  @code
#
#  import ostap.math.base 
#
#  import cppyy
#  cpp   = cppyy.gbl                           ## global C++ namespace 
#  Ostap = cpp.Ostap                           ## get C++ namespace Ostap
#
#  p3 = Ostap.XYZPoint(0,1,2)               ## use C++ type Ostap::XYZPoint
#
#  @endcode
#
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
"""Simple file to provide 'easy' access in python for the basic ROOT::Math classes

  see Ostap/Point3DTypes.h
  see Ostap/Vector3DTypes.h
  see Ostap/Vector4DTypes.h
  see Ostap/GenericVectorTypes.h
  see Ostap/Line.h

  The lines and planes are decorated:
     see Ostap/GeomFun.h

  The usage is fairly trivial:

  >>> import ostap.math.base 

  Important: All types are defined in corresponding
               C++ namespaces: Ostap & Ostap::Math

  >>> Ostap = cpp.Ostap                          ## get C++ namespace Gaudi
  >>> p3 = Ostap.XYZPoint(0,1,2)                 ## use C++ type Gaudi::XYZPoint

  >>> dir( Ostap.Math )
  >>> dir( Ostap      )
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = "Version$Revision$"
# =============================================================================
__all__     = (
    'iszero'         , ## zero     for doubles 
    'isequal'        , ## equality for doubles 
    'isint'          , ## Is equal  to int ? 
    'islong'         , ## Is equal  to long?
    'signum'         , ## sign of the number 
    'samesign'       , ## two number of the same sign 
    ##
    'inrange'        , ## is double number in certain range?
    ## 
    'natural_number' , ## See Ostap::Math::natural_number  
    'natural_entry'  , ## See Ostap::Math::natural_entry
    ##
    'make_vector'    , ## construct std::vector
    'doubles'        , ## construct std::vector<double>
    'ints'           , ## construct std::vector<int>
    'uints'          , ## construct std::vector<unsigned int>
    'longs'          , ## construct std::vector<long>
    'ulongs'         , ## construct std::vector<unsigned long>
    'complexes'      , ## construct std::vector<std::complex<double>>
    'strings'        , ## construct std::vector<std::string>
    ##
    'vDoubles'       , ## std::vector<double>
    'vFloats'        , ## std::vector<float>
    'vInts'          , ## std::vector<int>
    'vLongs'         , ## std::vector<long>
    ) 
# =============================================================================
import ROOT, cppyy 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.base' )
else                       : logger = getLogger ( __name__          )
# =============================================================================
## get global C++ namespace
cpp = cppyy.gbl

## C++ namespace Gaudi
std = cpp.std

## C++ namespace Ostap
Ostap = cpp.Ostap

iszero   = Ostap.Math.Zero     ('double')()
isequal  = Ostap.Math.Equal_To ('double')()
isint    = Ostap.Math.isint 
islong   = Ostap.Math.islong

vDoubles = std.vector ( 'double' )
vFloats  = std.vector ( 'float'  )
vInts    = std.vector ( 'int'    )
vLongs   = std.vector ( 'long'   )

# =============================================================================
##  get the sign of the number 
def signum ( x ) :
    """Get the sign of the number
    >>>  signum ( -10  ) , signum(0),   signum ( +2.5 ) 
    """
    ### for integers
    from ostap.core.types import is_integer as _is_integer 
    if _is_integer ( x ) : return 0 if 0 == x else +1 if 0<x else -1
    ## for floating numbers
    return 0 if iszero ( x ) else +1 if 0 < x else -1

# =============================================================================
## the same sign ?
def samesign ( a , b ) :
    """The same sign for two numbers?
    """
    return ( 0 < a and 0 < b ) or ( 0 > a and 0 > b ) 

# =============================================================================
## natural number ?
natural_number = Ostap.Math.natural_number
## natural etry in histo-bin ? 
natural_entry  = Ostap.Math.natural_entry

# ==============================================================================
## Is float value in range?
#  @code
#  x   = 1.1
#  a,b = 1,2
#  print 'Is %s<=%s<=%s ? %s' % ( a , x , b , inrange ( a , x , b ) )
#  @endcode 
def inrange ( a , x , b ) :
    """Is float value in range?
    >>> x   = 1.1
    >>> a,b = 1,2
    >>> print 'Is x between a and b? %s' % inrange ( a , x , b ) 
    """
    _a  = float(a)
    _b  = float(b)
    _x  = float(x)
    return ( _a <= _x or isequal ( _a , _x ) ) and ( _x <= _b or isequal ( _x , _b ) ) 


# =============================================================================
## decorate some basic std::vectors
for t in ( 'int'                  ,
           'long'                 ,
           'long long'            ,
           'unsigned int'         ,
           'unsigned long'        ,
           'unsigned long long'   ,
           'float'                ,
           'double'               ,
           'std::complex<double>' , 
           'std::string'          ) :
    v = std.vector( t )
    v.asList   = lambda s :       [ i for i in s ]   ## convert vector into list
    v.toList   = v.asList
    v.asTuple  = lambda s : tuple ( s.asList() )
    v.toTuple  = v.asTuple
    v.__repr__ = lambda s : str ( [ i for i in s ] ) ## print it !
    v.__str__  = lambda s : str ( [ i for i in s ] ) ## print it !
    if not hasattr ( v , '__iter__' ) :
        def _v_iter_ ( v ) :
            _l = len(v)
            for i in  range(_l) :
                yield v[i]
        v.__iter__ = _v_iter_ 
    

# =============================================================================
## self-printout of TMaxtrix 
def _tmg_str_ ( self , fmt = ' %+11.4g') :
    """Self-printout of TMatrix
    """
    _rows = self.GetNrows()
    _cols = self.GetNcols()
    _line = ''
    for _irow in range ( 0 , _rows ) :
        _line += ' |'
        for _icol in range ( 0 , _cols ) :
            _line += fmt % self( _irow , _icol )
        _line += ' |'
        if ( _rows - 1 )  != _irow : _line += '\n'
    return _line


ROOT.TMatrix.__repr__  = _tmg_str_
ROOT.TMatrix.__str__   = _tmg_str_


    
# =============================================================================
## add something to std::vector 
def _add_to ( vct , cnv , arg1 , *args ) :
    """Add something to std::vector 
    """
    from types       import GeneratorType as GT
    from collections import Iterable      as IT

    VT = type(vct)

    VS = std.vector('std::string')
    
    ## the special treatment of vector of  strings 
    if isinstance ( vct , VS ) :
        if isinstance ( arg1 , ( str , std.string ) ) :
            vct.push_back ( arg1 )
        else                                          :
            for a in arg1 : vct.push_back ( cnv ( a ) )
    ## the first argument: iterable or generator?
    elif isinstance ( arg1 , VT ) :
        for a in arg1 : vct.push_back ( a )
    elif isinstance ( arg1 , ( GT , IT ) ) :
        for a in arg1 : vct.push_back ( cnv ( a ) )                
    else :        
        try :
            vct.push_back ( cnv ( arg1 ) )
        except TypeError :
            for a in arg1 : _add_to (    vct , cnv , a )

    ##  treat other arguments (recursively) 
    for a in args : _add_to ( vct , cnv , a )

    return vct 

# =============================================================================
##  create C++ std-vector from components
#   @code
#   a = make_vector( 'string' , str   , 1 , 2, 3 , 4 )
#   a = make_vector( 'double' , float , 1 , 2, 3 , 4 )
#   a = make_vector( 'int'    , int   , [ 1,2,3 ] , [4,5,6] )
#   @endcode 
def make_vector ( TYPE , cnv , *args ) :
    """ Create C++ std::vector<TYPE> from components
    TYPE : the internal type
    cnv  : the converter (if needed)
    arg1 : the first mandatory arugment
    args : other components
    
    :Examples:

    >>> a = make_vector ( 'string' , str   , 1 , 2, 3 , 4 )
    >>> a = make_vector ( 'double' , float , 1 , 2, 3 , 4 )
    >>> a = make_vector ( 'int'    , int   , [ 1,2,3 ] , [4,5,6] )
    
    """
    ## create new vector 
    VT  = std.vector( TYPE ) ## vector type
    vct = VT ( )             ## vector instance
    if not  args : return vct 
    ## add arguments to the vector
    return _add_to ( vct , cnv , args[0] , *args[1:] )

# =============================================================================
## construct std::vector<double> from the arguments
def doubles ( arg1 , *args ) :
    """Construct the std::vector<double> from the arguments
    >>> v1 = doubles ( 1.01 )
    >>> v2 = doubles ( 1.01 , 1.02 , 1.03  )
    >>> v3 = doubles ( [ 1.01 , 1.02 , 1.03 ] )    
    """
    return make_vector ( 'double' , float , arg1 , *args )

# =============================================================================
## construct std::vector<ints> from the arguments
def ints ( arg1 , *args ) :
    """Construct the std::vector<int> from the arguments    
    >>> v1 = ints ( 1 )
    >>> v2 = ints ( 1 , 1 , 1  )
    >>> v3 = ints ( [ 1 , 2 , 3 ] )    
    """
    return make_vector ( 'int' , int , arg1 , *args )

# =============================================================================
## construct std::vector<unsigned int> from the arguments
def uints ( arg1 , *args ) :
    """Construct the std::vector<unsigned int> from the arguments    
    >>> v1 = uints ( 1 )
    >>> v2 = uints ( 1 , 1 , 1  )
    >>> v3 = uints ( [ 1 , 2 , 3 ] )    
    """
    return make_vector ( 'unsigned int' , long , arg1 , *args )

# =============================================================================
## construct std::vector<long> from the arguments
def longs ( arg1 , *args ) :
    """Construct the std::vector<long> from the arguments    
    >>> v1 = longs ( 1 )
    >>> v2 = longs ( 1 , 1 , 1  )
    >>> v3 = longs ( [ 1 , 2 , 3 ] )    
    """
    return make_vector ( 'long' , long , arg1 , *args )

# =============================================================================
## construct std::vector<unsigned long> from the arguments
def ulongs ( arg1 , *args ) :
    """Construct the std::vector<unsigned long> from the arguments    
    >>> v1 = ulongs ( 1 )
    >>> v2 = ulongs ( 1 , 1 , 1  )
    >>> v3 = ulongs ( [ 1 , 2 , 3 ] )    
    """
    return make_vector ( 'unsigned long' , long , arg1 , *args )

# =============================================================================
## construct std::vector<std::string> from the arguments
def strings ( *args ) :
    """Construct the std::vector<string> from the arguments
    >>>  v1 = strings( 'a','b','c')
    >>>  v2 = strings( ['a','b','c'] )
    >>>  v3 = strings( ['a','b'],'c',['d','e'] )
    """

    return make_vector ( 'std::string' , std.string , *args ) 
    ##  VS = std.vector(std.string) 
    ## vs = VS()
    ## if not args :  return vs 
    ## return _add_to ( vs , std.string , args[0] , *args[1:] )

SPD = std.pair('double','double')
SPD.asTuple  = lambda s :      (s.first,s.second)
SPD.__str__  = lambda s : str( (s.first,s.second) )
SPD.__repr__ = SPD.__str__

# =============================================================================
# Improve operations with std.complex 
# =============================================================================
COMPLEX  = cpp.std.complex('double'      )
COMPLEXf = cpp.std.complex('float'       )
COMPLEXl = cpp.std.complex('long double' )
# =============================================================================
def _cmplx_to_complex_ ( s ) :
    """Convert C++ complex to Python's complex"""
    return  complex    ( s.real() , s.imag() )

# =============================================================================
def _cmplx_negate_     ( s ) :
    """Negation:
    >>> v  = ...
    >>> v1 = -v
    """
    return -complex    ( s.real() , s.imag() )

# =============================================================================
def _cmplx_abs_        ( s ) :
    """Absolute value
    >>> print abs(v) 
    """
    import math
    sr = s.real()
    si = s.imag()
    return math.sqrt( sr * sr + si * si ) 

# =============================================================================
def _cmplx_norm_       ( s ) :
    """Norm (squared absolute value)
    >>> print v.norm()
    """
    sr = s.real()
    si = s.imag()
    return sr * sr + si * si

# =============================================================================
def _cmplx_conjugate_  ( s ) :
    """Get complex conjugated
    >>> vc = v.conjugate() 
    """
    return complex     ( s.real() , -s.imag() )
    
# =============================================================================
def _cmplx_add_        ( s , o ) :
    """add complex values 
    >>> r = v + other  
    """
    return o + complex ( s.real() , s.imag() )

# =============================================================================
def _cmplx_mul_        ( s , o ) :
    """multiply  complex values 
    >>> r = v * other  
    """
    return o * complex ( s.real() , s.imag() )

# =============================================================================
def _cmplx_div_        ( s , o ) :
    """divide complex values 
    >>> r = v / other  
    """
    return (1.0/o) * complex ( s.real() , s.imag() )

# =============================================================================
def _cmplx_rdiv_       ( s , o ) :
    """divide complex values 
    >>> r = other / v 
    """
    return o       * ( 1.0 / complex ( s.real() , s.imag() ) )

# =============================================================================
def _cmplx_sub_        ( s , o ) :
    """subtract complex values 
    >>> r = v - other 
    """
    return (-o   ) + complex ( s.real() , s.imag() )

# =============================================================================
def _cmplx_rsub_       ( s , o ) :
    """subtract complex values 
    >>> r = other - v 
    """
    return   o     - complex ( s.real() , s.imag() )

# =============================================================================
def _cmplx_pow_  ( s , o ) :
    """power function 
    >>> r = v ** other  
    """
    if isinstance ( o , COMPLEX ) :
        o = complex ( o.real() , o.imag() ) 
    return complex ( s.real() , s.imag() ) ** o

# =============================================================================
def _cmplx_rpow_  ( s , o ) :
    """power function 
    >>> r = other **v  
    """
    return o ** complex ( s.real() , s.imag() )


# =============================================================================
def _cmplx_eq_    ( s , o ) :
    """equality:
    >>> r = v == other  
    """
    if isinstance ( o, COMPLEX ) :
        return s.real() == o.real() and s.imag() == o.imag()
    return complex( s.real() , s.imag() ) == o

# =============================================================================
def _cmplx_ne_    ( s , o ) :
    """non-equality:
    >>> r = v != other  
    """
    if isinstance ( o, COMPLEX ) :
        return s.real() != o.real() or  s.imag() != o.imag()
    return complex( s.real() , s.imag() ) != o 

# =============================================================================
def _cmplx_iadd_ ( s , o ) :
    x = s + o
    s.real(x.real)
    s.imag(x.imag)
    
# =============================================================================
def _cmplx_isub_ ( s , o ) :
    x = s - o
    s.real(x.real)
    s.imag(x.imag)

# =============================================================================
def _cmplx_imul_ ( s , o ) :
    x = s * o
    s.real(x.real)
    s.imag(x.imag)

# =============================================================================
def _cmplx_idiv_ ( s , o ) :
    x = s / o
    s.real(x.real)
    s.imag(x.imag)

# =============================================================
for CMPLX in ( COMPLEX , COMPLEXf , COMPLEXl ) :
    
    if not hasattr ( CMPLX , '_old_init_' ) : 
        CMPLX._old_init_  = CMPLX.__init__
        ## construct complex 
        def _cmplx_new_init_ ( s , a = 0 , *b ) :
            """Construct complex from complex or from real/imaginary parts
            """
            if not b :
                a     = complex  ( a )
                a , b = a.real , a.imag
            elif 1 != len( b ) :
                raise TypeError("Can't create std::complex!")
            else :
                b =   b[0]
            
            return s._old_init_ ( a , b )
                
        CMPLX.__init__    = _cmplx_new_init_
        
    CMPLX.__complex__ = _cmplx_to_complex_
    
    CMPLX.__add__     = _cmplx_add_
    CMPLX.__mul__     = _cmplx_mul_
    CMPLX.__div__     = _cmplx_div_
    CMPLX.__sub__     = _cmplx_sub_
    
    CMPLX.__radd__    = _cmplx_add_
    CMPLX.__rmul__    = _cmplx_mul_
    CMPLX.__rdiv__    = _cmplx_rdiv_
    CMPLX.__rsub__    = _cmplx_rsub_
    
    CMPLX.__iadd__    = _cmplx_iadd_
    CMPLX.__imul__    = _cmplx_imul_
    CMPLX.__idiv__    = _cmplx_idiv_
    CMPLX.__isub__    = _cmplx_isub_
    
    CMPLX.__repr__    = lambda s : "%s" % complex ( s.real(), s.imag() )
    CMPLX.__str__     = lambda s : "%s" % complex ( s.real(), s.imag() )
    CMPLX.__abs__     = _cmplx_abs_
    CMPLX.__pow__     = _cmplx_pow_
    CMPLX.__rpow__    = _cmplx_rpow_
    CMPLX.__neg__     = _cmplx_negate_
    
    CMPLX.__eq__      =  _cmplx_eq_
    CMPLX.__ne__      =  _cmplx_ne_
    
    if not hasattr ( CMPLX , 'cpp_conj' ) :
        CMPLX.cpp_conj = lambda s : CMPLX ( s.real() , -s.imag() )
    
    CMPLX.norm        = _cmplx_norm_
    CMPLX.conjugate   = _cmplx_conjugate_
    CMPLX.conj        = _cmplx_conjugate_
    CMPLX.to_complex  = _cmplx_to_complex_ 
    CMPLX.as_complex  = _cmplx_to_complex_ 


# =============================================================================
## construct std::vector<std::complex> from the arguments
def complexes ( arg1 , *args ) :
    """Construct the std::vector<std::complex<double>> from the arguments    
    >>> v1 = complexs( 1+2j )
    >>> v2 = complexs( 1 , 1+2j , 1  )
    >>> v3 = complexs( [ 1 , 2 , 3+3j ] )    
    """
    return make_vector( COMPLEX , COMPLEX , arg1 , *args )

# =============================================================================
## complex value ?   
def is_complex ( value ) :
    """Complex value?
    """
    return isinstance ( value , ( complex, COMPLEX , COMPLEXf  , COMPLEXl ) )

## decorated classes 
_decorated_classes_  = (
    COMPLEX ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    _v = [ l for l in dir(Ostap     ) if 0 != l.find('__') ]
    print ' dir(Ostap)      : '
    _v.sort()
    for v in _v : print v
    _v = [ l for l in dir(Ostap.Math) if 0 != l.find('__') ]
    print ' dir(Ostap.Math) : '
    _v.sort()
    for v in _v : print v


# =============================================================================
# The  END
# =============================================================================
