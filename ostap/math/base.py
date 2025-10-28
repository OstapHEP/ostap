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
""" Simple file to provide 'easy' access in python for the basic ROOT::Math classes

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
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2009-09-12"
__version__ = "Version: $Revision$"
# =============================================================================
__all__     = (
    'iszero'         , ## zero     for doubles 
    'isequal'        , ## equality for doubles 
    'isint'          , ## Is equal  to int ? 
    'islong'         , ## Is equal  to long?
    'signum'         , ## sign of the number 
    'samesign'       , ## two number of the same sign
    'isfinite'       , ## `isfinite` for float values 
    'isnan'          , ## `isfinite` for float values 
    'isclose'        , ## `isclose`  for float values
    'lround'         , ## round a value to integer/long  
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
    ##
    'frexp10'        , ## similar to math.frexp but with radix=10,
    ## 
    'cpp'            , ## C++ global  namespace 
    'Ostap'          , ## C++ namespace Ostap 
    'std'            , ## C++ namespace Ostap
    'axis_range'     , ## suitable axis range
    ## 
    'gcd'            , ## gcd-function 
    'lcm'            , ## lcm-function
    ##
    'ROOTIgnore'     , ## control ROOT verbosity, suppress ROOT errors
    ##
    'complex_types'  , ## list of complex & complex-like types
    ##
    'pos_infinity'   , ## positive infinity  
    'neg_infinity'   , ## negative infinity
    ##
    'FIRST_ENTRY'    , ## the first entryfor evetn loops
    'LAST_ENTRY'     , ## the last entry for event loops 
    'evt_range'      , ## get the actual range of entries
    'all_entries'    , ## Are all entreis required to process? 
    ##
    'numpy'          , ## numpy or None
    'scipy'          , ## scipy or None
    'np2raw'         , ## numpy array to raw C++ buffer 
    ) 
# =============================================================================
from   ostap.core.meta_info    import python_info
from   ostap.core.ostap_types  import sequence_types, sized_types 
from   collections.abc         import Sized
import ROOT, cppyy, sys, math, ctypes, array   
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.base' )
else                       : logger = getLogger ( __name__          )
# =============================================================================
## get global C++ namespace
cpp   = cppyy.gbl
# =============================================================================
## C++ namespace std 
std   = cpp.std
# =============================================================================
## C++ namespace Ostap
Ostap = cpp.Ostap 

# =============================================================================
## positive and negative infinities 
pos_infinity = float('+inf')
neg_infinity = float('-inf')

# =============================================================================
## Very simple context manager to suppress ROOT printout
#  @code
#  >>> with ROOTIgnore( ROOT.kError + 1 ) : some_ROOT_code_here()
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-30
class ROOTIgnore( object ) :
    """ Very simple context manager to suppress ROOT printout
    >>> with ROOTIgnore ( ROOT.kError + 1 ) : some_ROOT_code_here()
    """
    ## constructor
    #  @param level  (INPUT) print level 
    #  @param silent (print level 
    # 
    def __init__ ( self , level ) :
        """ Constructor:        
        >>> with rootError   () : some_ROOT_code_here()
        >>> with rootWarning () : some_ROOT_code_here()
        """
        #
        self._level = int ( level )
        
    ## context manager: ENTER 
    def __enter__ ( self ) :
        "The actual context manager: ENTER"
        self._old = int ( ROOT.gErrorIgnoreLevel ) 
        if self._old != self._level :
            groot = ROOT.ROOT.GetROOT()
            groot.ProcessLine("gErrorIgnoreLevel= %d ; " % self._level ) 
            
        return self
    
    ## context manager: EXIT 
    def __exit__ ( self , *_ ) : 
        "The actual context manager: EXIT"
        if self._old != int ( ROOT.gErrorIgnoreLevel )  :
            groot = ROOT.ROOT.GetROOT()            
            groot.ProcessLine("gErrorIgnoreLevel= %d ; " % self._old ) 
            
# =============================================================================
from ostap.logger.mute  import mute
logger.debug ("Suppress error/warnings from ROOT")
with ROOTIgnore ( ROOT.kWarning + 1 ) : 
    with mute ( True  , True ) : _ = ROOT.RooRealVar() 
    iszero   = Ostap.Math.Zero     ('double')()
    isequal  = Ostap.Math.Equal_To ('double')()
    isequalf = Ostap.Math.Equal_To ('float' )()
    isint    = Ostap.Math.isint 
    islong   = Ostap.Math.islong
    
vDoubles = std.vector ( 'double' )
vFloats  = std.vector ( 'float'  )
vInts    = std.vector ( 'int'    )
vLongs   = std.vector ( 'long'   )

# =============================================================================
## local version of <code>isfinite</code>
isfinite = math.isfinite 
## local version of <code>isnana</code>
isnan    = math.isnan
    
# =============================================================================
## local version of <code>isclose</code>
isclose = math.isclose 
            
# =============================================================================
##  get the sign of the number 
def signum ( x ) :
    """ Get the sign of the number
    >>> signum ( -10  ) , signum(0),   signum ( +2.5 ) 
    """
    ### for integers
    from ostap.core.ostap_types import is_integer as _is_integer 
    if _is_integer ( x ) : return 0 if 0 == x else +1 if 0<x else -1
    ## for floating numbers
    return 0 if iszero ( x ) else +1 if 0 < x else -1

# =============================================================================
## the same sign ?
def samesign ( a , b ) :
    """ The same sign for two numbers?
    """
    return ( 0 < a and 0 < b ) or ( 0 > a and 0 > b ) 

# =============================================================================
## natural number ?
natural_number = Ostap.Math.natural_number
# =============================================================================
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
    """ Is float value in range?
    >>> x   = 1.1
    >>> a,b = 1,2
    >>> print 'Is x between a and b? %s' % inrange ( a , x , b ) 
    """
    _a  = float ( a )
    _b  = float ( b )
    _x  = float ( x )
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
    v.asTuple  = lambda s : tuple ( i for i in s ) 
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
## Convert vector of floating types to string 
def float_vct_str ( vct , format = '%.5g' ) :
    """ Convert vector of floating types to string"""
    try :
        return '[ ' + ', '.join ( [ format % v for v in vct ] ) + ' ]'  
    except TypeError :
        pass
    return float_vct_str ( vct , format = '%.5g' )

# =============================================================================
## Convert vector of complex types to string 
def complex_vct_str ( vct , format = '%.5g%-+.5gj' ) :
    """ Convert vector of complex types to string"""
    try :
        lst = [] 
        for c in vct :
            cc   = complex ( c )
            item = format % ( cc.real , cc.imag )
            lst.append ( item )        
        return '[ ' + ', '.join ( lst ) + ' ]'  
    except TypeError :
        pass
    return complex_vct_str ( vct , format = '%.5g%-+.5gj' )

for t in ( 'float' , 'double' ):
    
    v = std.vector ( t )
    v.__repr__ = float_vct_str 
    v.__str__  = float_vct_str 

for t in ( 'std::complex<double>' , 'std::complex<float>'  ) :

    v = std.vector( t )
    v.__repr__ = complex_vct_str    
    v.__str__  = complex_vct_str    
    
# =============================================================================
## add something to <code>std::vector</code>
def _add_to ( vct , cnv , arg1 , *args ) :
    """ Add something to `std.vector`
    """
    from ostap.core.ostap_types import sequence_types
    
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
    elif isinstance ( arg1 , sequence_types ) :      ## SEQUENCE
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
    
    vct.reserve ( len ( args ) )
    
    if not  args : return vct 
    ## add arguments to the vector
    return _add_to ( vct , cnv , args[0] , *args[1:] )

# =============================================================================
## construct std::vector<double> from the arguments
def doubles ( arg1 , *args ) :
    """ Construct the std::vector<double> from the arguments
    >>> v1 = doubles ( 1.01 )
    >>> v2 = doubles ( 1.01 , 1.02 , 1.03  )
    >>> v3 = doubles ( [ 1.01 , 1.02 , 1.03 ] )    
    """
    return make_vector ( 'double' , float , arg1 , *args )

# =============================================================================
## construct std::vector<ints> from the arguments
def ints ( arg1 , *args ) :
    """ Construct the std::vector<int> from the arguments    
    >>> v1 = ints ( 1 )
    >>> v2 = ints ( 1 , 1 , 1  )
    >>> v3 = ints ( [ 1 , 2 , 3 ] )    
    """
    return make_vector ( 'int' , int , arg1 , *args )

# =============================================================================
## construct std::vector<unsigned int> from the arguments
def uints ( arg1 , *args ) :
    """ Construct the std::vector<unsigned int> from the arguments    
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

# =============================================================================
SPD = std.pair ( 'double' , 'double' )
SPD.asTuple  = lambda s :       ( s.first , s.second )
SPD.__str__  = lambda s : str(  ( s.first , s.second ) )
SPD.__repr__ = SPD.__str__
SPD.__len__  = lambda s : 2
def _spd_getitem_ ( s ,  i ) :
    """ Get item from the pair:
    >>> p = ...
    >>> p[0], p[1]
    """
    if    0 == i : return s.first  ## first 
    elif  1 == i : return s.second ## second 
    elif -1 == i : return s.second ## last 
    raise IndexError('Invalid index %s' % i )
def _spd_setitem_ ( s ,  i , value ) :
    """Set item from the pair:
    >>> p = ...
    >>> p[0] = 1
    >>> p[1] = 2
    """
    if    0 == i : s.first  = value ## first 
    elif  1 == i : s.second = value ## second 
    elif -1 == i : s.second = value ## last 
    raise IndexError('Invalid index %s' % i )
SPD.__getitem__ = _spd_getitem_
SPD.__setitem__ = _spd_setitem_
    
# =============================================================================
# Improve operations with std.complex 
# =============================================================================
COMPLEX    = cpp.std.complex ( 'double'      )
COMPLEXf   = cpp.std.complex ( 'float'       )
COMPLEXl   = cpp.std.complex ( 'long double' )

VCOMPLEX   = cpp.std.vector ( COMPLEX )
VDOUBLE    = cpp.std.vector ('double' )
VCT_TYPES  = VDOUBLE, VCOMPLEX
# =============================================================================
## list of complex types: native and X++
complex_types = complex , COMPLEX, COMPLEXf, COMPLEXl
# =============================================================================
def _real_ ( s ) : return s.real 
def _imag_ ( s ) : return s.imag 
def _cmplx_to_complex_ ( s ) :
    """ Convert C++ complex to Python's complex"""
    return  complex    ( s.real , s.imag  )

# =============================================================================
def _cmplx_negate_     ( s ) :
    """ Negation:
    >>> v  = ...
    >>> v1 = -v
    """
    return -complex    ( s )

# =============================================================================
def _cmplx_abs_        ( s ) :
    """ Absolute value
    >>> print abs(v) 
    """
    import math 
    return math.sqrt ( s.norm () )

# =============================================================================
def _cmplx_norm_       ( s ) :
    """ Norm ( squared absolute value )
    >>> print v.norm()
    """
    sr = _real_ ( s ) 
    si = _imag_ ( s ) 
    return sr * sr + si * si

# =============================================================================
def _cmplx_conjugate_  ( s ) :
    """ Get complex conjugated
    >>> vc = v.conjugate() 
    """
    return complex ( _real_ ( s )  , - _imag_  ( s )  )
    
# =============================================================================
def _cmplx_add_        ( s , o ) :
    """ Add complex values 
    >>> r = v + other  
    """
    return o + complex ( s )

# =============================================================================
def _cmplx_mul_        ( s , o ) :
    """ Multiply  complex values 
    >>> r = v * other  
    """
    return o * complex ( s  )

# =============================================================================
def _cmplx_div_        ( s , o ) :
    """ Divide complex values 
    >>> r = v / other  
    """
    return ( 1.0 / o ) * complex ( s )

# =============================================================================
def _cmplx_rdiv_       ( s , o ) :
    """ Divide complex values 
    >>> r = other / v 
    """
    return o           * ( 1.0 / complex ( s ) )

# =============================================================================
def _cmplx_sub_        ( s , o ) :
    """ Subtract complex values 
    >>> r = v - other 
    """
    return (-o   ) + complex ( s )

# =============================================================================
def _cmplx_rsub_       ( s , o ) :
    """ Subtract complex values 
    >>> r = other - v 
    """
    return   o     - complex ( s )

# =============================================================================
def _cmplx_pow_  ( s , o ) :
    """ Power function 
    >>> r = v ** other  
    """
    if isinstance ( o , COMPLEX ) : o = complex ( o ) 
    return complex ( s ) ** o

# =============================================================================
def _cmplx_rpow_  ( s , o ) :
    """ Power function 
    >>> r = other ** v  
    """
    return o ** complex ( s )

# =============================================================================
def _cmplx_eq_    ( s , o ) :
    """ Equality:
    >>> r = v == other  
    """
    if isinstance ( o, COMPLEX ) :
        return _real_ ( s ) == _real_ ( o ) and _imag_ ( s ) == _imag_ ( o )
    return complex ( s ) == o

# =============================================================================
def _cmplx_ne_    ( s , o ) :
    """ Non-equality:
    >>> r = v != other  
    """
    if isinstance ( o, COMPLEX ) :
        return _real_ ( s ) != _real_ ( o ) or _imag_ ( s ) != _imag_ ( o )
    return complex ( s ) != o 

# ==============================================================================
## Deserialize the complex numbers
def _cmplx_factory_ ( cmplxt , re , im ) :
    """ Deserialize the complex numbers
    """
    return cmplxt ( re , im )
# ==============================================================================
## reduce complex numbers 
def _cmplx_reduce_ ( c ) :
    """ Reduce complex numbers"""
    return _cmplx_factory_ , ( type ( c ) , c.real , c.imag )
    
# =========================================================================
def _cmplx_iadd_ ( s , o ) :
    x = s + o
    T = type ( s )
    t = T ( x.real , x.imag )
    s.__assign__ ( t )         
    return s
# =========================================================================
def _cmplx_isub_ ( s , o ) :
    x = s - o
    T = type ( s )
    t = T ( x.real , x.imag )
    s.__assign__ ( t )         
    return s
# =========================================================================
def _cmplx_imul_ ( s , o ) :
    x = s * o
    T = type ( s )
    t = T ( x.real , x.imag )
    s.__assign__ ( t )         
    return s
# =========================================================================
def _cmplx_idiv_ ( s , o ) :
    x = s / o
    T = type ( s )
    t = T ( x.real , x.imag )
    s.__assign__ ( t ) 
    return s
        
# =============================================================================
for CMPLX in ( COMPLEX , COMPLEXf , COMPLEXl ) :
    
    if not hasattr ( CMPLX , '_old_init_' ) : 
        CMPLX._old_init_  = CMPLX.__init__
        ## construct complex 
        def _cmplx_new_init_ ( s , a = 0 , *b ) :
            """ Construct complex from complex or from real/imaginary parts
            """
            if not b :
                a     = complex  ( a )
                a , b = a.real , a.imag
            elif 1 != len( b ) :
                raise TypeError("Can't create std::complex!")
            else :
                b =   b [ 0 ]
            
            return s._old_init_ ( a , b )
                
        CMPLX.__init__ = _cmplx_new_init_
        
    CMPLX.__complex__  = _cmplx_to_complex_
    
    CMPLX.__add__      = _cmplx_add_
    CMPLX.__mul__      = _cmplx_mul_
    CMPLX.__div__      = _cmplx_div_
    CMPLX.__sub__      = _cmplx_sub_
    CMPLX.__truediv__  = _cmplx_div_
    
    CMPLX.__radd__     = _cmplx_add_
    CMPLX.__rmul__     = _cmplx_mul_
    CMPLX.__rdiv__     = _cmplx_rdiv_
    CMPLX.__rsub__     = _cmplx_rsub_
    CMPLX.__rtruediv__ = _cmplx_rdiv_
    
    CMPLX.__iadd__     = _cmplx_iadd_
    CMPLX.__imul__     = _cmplx_imul_
    CMPLX.__idiv__     = _cmplx_idiv_
    CMPLX.__isub__     = _cmplx_isub_
    CMPLX.__itruediv__ = _cmplx_idiv_

    CMPLX.__repr__    = lambda s : "%s" % complex ( s )
    CMPLX.__str__     = lambda s : "%s" % complex ( s )
    CMPLX.__abs__     = _cmplx_abs_
    CMPLX.__pow__     = _cmplx_pow_
    CMPLX.__rpow__    = _cmplx_rpow_
    CMPLX.__neg__     = _cmplx_negate_
    
    CMPLX.__eq__      =  _cmplx_eq_
    CMPLX.__ne__      =  _cmplx_ne_
    
    if not hasattr ( CMPLX , 'cpp_conj' ) :
        CMPLX.cpp_conj = lambda s : CMPLX ( _real_ ( s ) , -_imag_ ( s ) )
    
    CMPLX.norm        = _cmplx_norm_
    CMPLX.conjugate   = _cmplx_conjugate_
    CMPLX.conj        = _cmplx_conjugate_
    CMPLX.to_complex  = _cmplx_to_complex_ 
    CMPLX.as_complex  = _cmplx_to_complex_ 

    ## make python 3 happy!
    if not hasattr ( CMPLX ,  '__truediv__' ) : CMPLX. __truediv__  = CMPLX. __div__ 
    if not hasattr ( CMPLX , '__itruediv__' ) : CMPLX.__itruediv__  = CMPLX.__idiv__ 
    if not hasattr ( CMPLX , '__rtruediv__' ) : CMPLX.__rtruediv__  = CMPLX.__rdiv__ 

    CMPLX.__reduce__   = _cmplx_reduce_ 



# =============================================================================
## construct std::vector<std::complex> from the arguments
def complexes ( arg1 , *args ) :
    """ Construct the std::vector<std::complex<double>> from the arguments    
    >>> v1 = complexs( 1+2j )
    >>> v2 = complexs( 1 , 1+2j , 1  )
    >>> v3 = complexs( [ 1 , 2 , 3+3j ] )    
    """
    return make_vector( COMPLEX , COMPLEX , arg1 , *args )

# =============================================================================
## complex value ?   
def is_complex ( value ) :
    """ Complex value?
    """
    return isinstance ( value , ( complex, COMPLEX , COMPLEXf  , COMPLEXl ) )

# =============================================================================
## decorated classes 
_decorated_classes_  = (
    COMPLEX  ,
    COMPLEXf ,
    COMPLEXl ,
    )


# =======================================================================
## nice printout of complex numbers ( string + exponent)
#  @code
#  ae = complex ( ... ) 
#  s , expo = pretty_complex (  ae ) 
#  @endcode
#  @return nice string and the separate exponent 
def pretty_complex ( value              ,
                     width       = 6    ,
                     precision   = 4    ,
                     parentheses = True ) :
    """ Nice printout of complex number ( string + exponent)
    - return nice stirng and the separate exponent 
    >>> s , expo = pretty_complex( number ) 
    """
    v  = complex ( value )
    ##
    from ostap.logger.pretty import fmt_pretty_values 
    fmtv , expo = fmt_pretty_values ( v.real () ,
                                      v.imag () ,
                                      width       = width     ,
                                      precision   = precision ,
                                      parentheses = False     )
    ## 
    fmt = '%s%sj' % ( fmtv , fmtv ) 
    if parentheses : fmt = '( ' + fmt + ' )' 
    ##
    if expo :
        scale = 10 ** expo
        v /= scale
    ## 
    return fmt % ( v.real , v.imag ) , expo

# =======================================================================
## Nice pritout of arrays of numerical values
#  @code
#  array = ...
#  result, expo = pretty_array ( array ) 
#  @endcode 
#  @return nice string and the separate exponent 
def pretty_array ( values             ,
                   width       = 6    ,
                   precision   = 4    ,
                   parentheses = True ) :
    """ Nice pritout of arrays of numerical values
    - return nice string and the separate exponent 
    >>> array = ...
    >>> result, expo = pretty_array ( array ) 
    """
    assert isinstance ( values , sequence_types ) , \
        "Invalid type of `values':%s" % type ( values )
    
    if not values : return '' if not parentheses else '[]'
    
    assert all ( isinstance ( v , num_types ) for v in values ) , \
        "Invalid content of `value': %s" % str ( values ) 

    from ostap.logger.pretty import fmt_pretty_values 
    fmtv , expo = fmt_pretty_values ( *values                 , 
                                      width       = width     ,
                                      precision   = precision )
    
    if expo :
        scale  = 10 ** expo
        result = ', '.join ( fmtv % ( float ( v ) / scale ) for v in values )
    else :
        result = ', '.join ( fmtv %   float ( v )           for v in values )
    ## 
    if parentheses : result = '[ ' + result + ' ]'
    return result, expo

# =============================================================================
## C++ version of frexp with radix 10 
cpp_frexp10 = Ostap.Math.frexp10 
# =============================================================================
## get mantissa (0.1<=m<1) and exponent for radix10
#  similar for frexp, but use radix=10
#  @code
#  m,e = frexp10 ( value ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-20
def frexp10 ( value ) :
    """ Get the mantissa (0.1<=m<1) and exponent for radix10
    (similar for frexp, but use radix=10)
    
    >>> a,b = frexp10 ( value ) 
    """

    ## a bit better treatment of near-zero numbers 
    p = cpp_frexp10 ( value )
    return p.first, p.second


    xv = abs ( value )
    if iszero ( xv )  : return  ( 0 , 0 ) 

    q  = math.floor ( math.log10 ( float ( xv ) ) )

    if    0 < q : xv /= 10**q 
    elif  0 > q : xv *= 10**abs ( q ) 
    
    if 1 <= xv :
        xv /= 10
        q  += 1
        
    return ( xv , q ) if ( 0 <= value ) else ( -xv , q ) 

# =============================================================================
## Define some "range" for the given value:
#  @code
#  value = ...
#  mn, mx = num_range ( value ) 
#  @endcode 
def num_range ( value , N = 1 ) : 
    """ Define some `range' for the given value:
    >>> value = ...
    >>> mn, mx = num_range ( value ) 
    """
    
    if   iszero ( value ) : return ( -0.5 , 0.5 )
    
    assert isinstance ( N , int ) and 1 <= N , \
        "Inavalid `N' parameter:%s" % N
    
    a , b = frexp10 ( value )         
    b  -= N 

    NN = 10 ** N
    
    if not isfinite ( a * 2* NN ):
        logger.error ( "num_range: not finite: %s %s %s %s %s" % ( value , N , a , b , NN ) )
        
    af = math.floor ( a * 2 * NN )
    ac = math.ceil  ( a * 2 * NN )
    
    if isequal ( af , ac ) :
        af, ac  = af - 0.5 , ac + 0.5
        
    xmin = af * ( 10 ** b ) * 0.5 
    xmax = ac * ( 10 ** b ) * 0.5 

    return xmin , xmax
    
# =============================================================================
## Find suitable range for histogram axis 
def axis_range ( xmin , xmax , delta = 0.02 , log = False ) :
    """ Find suitable range for histogram axis
    """
    xmn = min ( xmin , xmax )
    xmx = max ( xmin , xmax )
    
    ## 1) special case
    if iszero ( xmn ) and iszero ( xmx ) : return ( -1.0 , 1.0 )
    
    ## 2) special case 
    if isequal ( xmn , xmx ) : return num_range ( 0.5 * ( xmn + xmx ) , N = 2 ) 
    
    ## 3) special case
    if islong ( xmn - 0.5 ) and islong ( xmx + 0.5 ) :
        return math.floor ( xmn - 0.1 ) , math.ceil ( xmx + 0.1 ) 

    d = xmx - xmn

    delta = abs ( delta        )
    fr    = min ( delta , 0.9  ) 

    if iszero ( xmn ) and xmn < xmx :
        
        xmin = 0 
        xmax = xmx + delta * d 

    elif iszero ( xmx ) and xmn < xmx :
        
        xmin = xmn - delta * d 
        xmax = 0 
        
    elif xmn < 0 < xmx  :
        
        xmin = ( 1 + delta ) * xmn  
        xmax = ( 1 + delta ) * xmx
        
    elif 0 < xmn < xmx :
        
        xmin = max ( xmn * ( 1 - fr ) , xmn - delta * d )
        xmax =                          xmx + delta * d 
        
    elif xmn < xmx < 0 :
        
        xmin =                          xmn - delta * d 
        xmax = min ( xmx * ( 1 - fr ) , xmx + delta * d )
        
    else : 
    
        xmin = xmn - delta * d        
        xmax = xmx + delta * d
        
    xmin , _     = num_range ( xmin , N = 2 )
    _    , xmax  = num_range ( xmax , N = 2 )
     
    return xmin, xmax 

# =============================================================================
if   ( 3 , 9 ) <= python_info : # ========================================
    # =========================================================================
    ## Least Common Multiple.
    lcm = math.lcm
    # =========================================================================
    ## Greatest Common Divisor 
    gcd = math.gcd
    # =========================================================================
else : # ======================================================================
    # =========================================================================
    ## Greatest Common Divisor 
    gcd = math.gcd 
    ## Least Common Multiple.
    # =========================================================================
    def lcm ( a , b ) :
        """ Least Common Multiple """
        
        aa     = abs ( a )
        bb     = abs ( b )
        
        max_ab = max ( aa , bb )
        min_ab = min ( aa , bb )
        
        return min_ab * ( max_ab // gcd ( a , b ) ) 
    # =========================================================================

# =============================================================================
_round = Ostap.Math.round
# =============================================================================
## round the value
#  @see Ostap::Math::round 
def lround ( x ) :
    """ Round the value 
    - see `Ostap.Math.round`
    """
    return x if isinstance  ( x , int ) else _round ( x ) 

# ============================================================================
## The first entry for event loops
FIRST_ENTRY = Ostap.FirstEvent 
## The last entry for event loops 
LAST_ENTRY  = Ostap.LastEvent
# ============================================================================
assert isinstance ( FIRST_ENTRY , int ) , "Invalid First Entry type!"
assert isinstance ( LAST_ENTRY  , int ) , "Invalid Last  Entry type!"
assert 0 <= FIRST_ENTRY < LAST_ENTRY    , "Invalid First/Last entries!"
# ============================================================================
## Get the actual range of entries
#  @code
#  tree  = 
#  first , last = evt_range ( 100 , 1000 ) 
#  @endcode
def evt_range ( sized , first = FIRST_ENTRY , last = LAST_ENTRY  ) :
    """ Get the actual range of entries  
    >>> tree = ....
    >>> first , last = evt_range ( tree , 0 , 1000 ) 
    """
    assert isinstance ( first , int ) , "evt_range: Invalid `first' type!"
    assert isinstance ( last  , int ) , "evt_range: Invalid `last' type!"
    ##
    size = sized
    if isinstance ( sized , Sized ) : size  = len ( sized )
    assert isinstance ( size , int ) and 0 <= size , 'Invalid size!'
    ##
    if not size : return 0 , 0   ## empty range
    ##
    if first < 0 : first += size
    if last  < 0 : last  += size 
    assert 0 <= first <= last , "Invalid first/last setting!"
    ##
    if   size <= first : return size  , size 
    elif size <= last  : return first , min ( size , last )
    ##
    return first , last 
# =========================================================================
## Are all entries required to process?
#  @code
#  tree  = 
#  if not  all_events ( tree , 100 , 1000 ) :
#  @endcode
def all_entries ( sized , first = FIRST_ENTRY , last = LAST_ENTRY  ) :
    """ Get the actual range of entries  
    >>> tree = ....
    >>> if not all_entries  ( tree , 0 , 1000 ) : ...
    """
    ##
    assert isinstance ( first , int ) , "all_entries: Invalid `first' type!"
    assert isinstance ( last  , int ) , "evt_entries: Invalid `last' type!"
    ##
    size = sized
    if isinstance ( sized , Sized ) : size  = len ( sized )
    assert isinstance ( size , int ) and 0 <= size , 'Invalid size!'
    ##
    if first < 0 : first += size
    if last  < 0 : last  += size 
    assert 0 <= first <= last , "Invalid first/last setting!"
    ##
    return 0 == first and size <= last 

# =============================================================================
## Numpy & scipy 
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import numpy
    numpy_version = tuple ( int ( i ) for i in numpy.__version__.split( '.' ) )
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    numpy = None
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import scipy
    scipy_version = tuple ( int ( i ) for i in scipy.__version__.split( '.' ) )
    # ========================================================================
except ImportError :
    # ========================================================================
    scipy = None
    scipy_version = () 
# =============================================================================
np2raw = None
# =============================================================================
## Converrt numpy array into raw C++ buffer
#  @code
#  data = ...
#  raw , size = np2raw ( data ) 
#  @endcoder 
def np2raw ( data ) :
    """ Converrt numpy array into raw C++ buffer
    >>> data = ...
    >>> raw , size = np2raw ( data ) 
    """
    assert numpy and isinstance ( data , numpy.ndarray ) , \
        "np2raw: argument must be `numpy.ndarray`"
    size   = len ( data ) 
    dtype  = data.dtype
    ctype  = numpy.ctypeslib._ctype_from_dtype ( dtype  )
    buffer = data.ctypes.data_as ( ctypes.POINTER ( ctype ) )
    # 
    return buffer , size      

# =============================================================================
## Call a function of scalar argument with array-like argument
#  - a kind of straw-man vectorization
#  - function signature is `result = fun1 ( x , **kwargs )`
def vct1_call ( fun1 , x , *more , **kwargs ) :
    """ Call a function of scalar argument with array-like argument
    - a kind of straw-man vectorization
    - function signature is `result = fun1 ( x , **kwargs )`
    """    
    assert callable ( fun1 ) , "`fun2` must be callable!"
    
    ## (0) get the actual function, wrap arguments 
    if kwargs : fun = lambda x : fun1 ( x , **kwargs )
    else      : fun = fun1 
    
    ## (1)  sequence of arguments, "x"  must be scalar
    if more :
        more = ( x , ) + more 
        return tuple ( fun ( v ) for v in more )
    
    ## (2) single argument is a sequence        
    if isinstance ( x , sequence_types ) :
        gen = ( fun ( v ) for v in x )
        ## as numpy array or array.array 
        if   isinstance ( x , numpy.ndarray ) : return numpy.fromiter ( gen , dtype = float )
        elif isinstance ( x , array.array   ) : return array.array    ( 'd' , gen           ) 
        ## as simple tuple 
        return tuple ( gen )
    
    ## (3) argument is a scalar 
    return fun ( x )

# =============================================================================
## Call a function of two scalar arguments with array-like argument(s)
#  - a kind of straw-man vectorization
#  - function signature is `result = fun2 ( x , y , *args , **kwargs )`
def vct2_call ( fun2 , x , y ,  *args , **kwargs ) :
    """ Call a function of scalar argument with array-like argument
      - a kind of straw-man vectorization
      - function signature is `result = fun2  ( x , y , **kwargs )`
    """    
    assert callable ( fun2 ) , "`fun2` must be callable!"
    
    ## (0) get the actual function, wrap arguments 
    if args or kwargs : fun = lambda x, y : fun2 ( x , y , *args , **kwargs )
    else              : fun = fun2 
    
    ## (1) chech the arguments 
    xseq = isinstance ( x , sequence_types )
    yseq = isinstance ( y , sequence_types )

    ## (2) treat array-like arguments 
    if   xseq and yseq : gen = ( fun ( vx , vy ) for vx , vy in zip ( x ,          y   ) )
    elif xseq          : gen = ( fun ( vx , vy ) for vx , vy in zip ( x , repeat ( y ) ) )
    elif yseq          : gen = ( fun ( vx , vy ) for vy , vx in zip ( y , repeat ( x ) ) )
    else :
        ## (3) both arguments are scalars
        return fun ( x , y ) 
    
    ## (4) as numpy if any of arguments is numpy
    if ( xseq and isinstance ( x , numpy.ndarray ) ) or \
       ( yseq and isinstance ( y , numpy.ndarray ) ) : return numpy.fromiter ( gen , dtype = float )
    
    ## (5) as array if any of arguments is array 
    if ( xseq and isinstance ( x , array.array   ) ) or \
       ( yseq and isinstance ( y , array.array   ) ) : return array.array    ( 'd' , gen )
    
    ## (6) as simple tuple 
    return tuple ( gen ) 

# =============================================================================
## Call a function of three scalar arguments with array-like argument(s)
#  - a kind of straw-man vectorization
#  - function signature is `result = fun3 ( x , y , z , *args , **kwargs )`
def vct3_call ( fun3 , x , y , z , *args , **kwargs ) :
    """ Call a function of three scalar arguments with array-like argument(s)
    - a kind of straw-man vectorization
    - function signature is `result = fun3 ( x , y , z , *args , **kwargs )`
    """    
    assert callable ( fun3 ) , "`fun3` must be callable!"
    
    ## (0) get the actual function, wrap arguments 
    if args or kwargs : fun = lambda x, y, z : fun3 ( x , y , z , *args , **kwargs )
    else              : fun = fun3 

    ## (1) check argument types 
    xseq = isinstance ( x , sequence_types )
    yseq = isinstance ( y , sequence_types )
    zseq = isinstance ( z , sequence_types )

    ## (2) treat array-like arguments 
    if   xseq and yseq and zseq : gen = ( fun ( vx , vy , vz ) for vx , vy , vz in zip ( x , y ,          z   ) )
    elif xseq and yzeq          : gen = ( fun ( vx , vy , vz ) for vx , vy , vz in zip ( x , y , repeat ( z ) ) ) 
    elif xseq and zseq          : gen = ( fun ( vx , vy , vz ) for vx , vz , vy in zip ( x , z , repeat ( y ) ) )
    elif yseq and zseq          : gen = ( fun ( vx , vy , vz ) for vy , vz , vx in zip ( y , z , repeat ( x ) ) )
    elif xseq                   : gen = ( fun ( vx , vy , vz ) for vx , vy , vz in zip ( x , repeat ( y ) , repeat ( z ) ) )
    elif yseq                   : gen = ( fun ( vx , vy , vz ) for vy , vx , vz in zip ( y , repeat ( x ) , repeat ( z ) ) )
    elif zseq                   : gen = ( fun ( vx , vy , vz ) for vz , vx , vy in zip ( z , repeat ( x ) , repeat ( y ) ) )
    else :
        ## (3) all arguments are scalars! 
        return fun ( x , y , z ) 

    ## (4) as numpy if any of arguments is numpy
    if ( xseq and isinstance ( x , numpy.ndarray ) ) or \
       ( yseq and isinstance ( y , numpy.ndarray ) ) or \
       ( zseq and isinstance ( z , numpy.ndarray ) ) : return numpy.fromiter ( gen , dtype = float )
    
    ## (5) as array if any of arguments is array 
    if ( xseq and isinstance ( x , array.array   ) ) or \
       ( yseq and isinstance ( y , array.array   ) ) or \
       ( zseq and isinstance ( z , array.array   ) ) : return array.array    ( 'd' , gen )
    
    ## (6) as simple tuple 
    return tuple ( gen ) 

# =============================================================================
## vector call wrappers 
# =============================================================================
## decorator to enhance call-method for certain class 
def vct1_call_method ( method ) :
    """ decorator to enhace call-method for certain class
    """
    def decorated_call1 ( who , x , *more , **kwargs ) :
        fun1 = lambda x : method ( who , x , **kwargs )
        return vct1_call ( fun1 , x , *more )
    return decorated_call1

# =============================================================================
## decorator to enhance call-method for certain class 
def vct2_call_method ( method ) :
    """ decorator to enhace call-method for certain class
    """
    def decorated_call2 ( who , x , y ,  *args , **kwargs ) :
        fun2 = lambda x, y : method ( who , x , y , *args ,  **kwargs  )
        return vct2_call ( fun2 , x , y )
    return decorated_call2

# =============================================================================
## decorator to enhance call-method for certain class 
def vct3_call_method ( method ) :
    """ decorator to enhance call-method for certain class
    """
    def decorated_call3 ( who , x , y , z , *args , **kwargs ) :
        fun3 = lambda x, y, z : method ( who , x , y , z , *args , **kwargs )
        return vct3_call ( fun3 , x , y , z )
    return decorated_call3 

# =============================================================================
## imports at the end of the module to avoid the circular dependency 
# =============================================================================

import ostap.math.reduce  
import ostap.math.polynomials 
    
# =============================================================================
if not '__main__' == __name__ : # =============================================
    # =========================================================================
    from ostap.io.checker import PickleChecker as Checker 
    checker = Checker ()

    checker.add ( *complex_types )
    checker.add ( *( std.vector ( t ) for t in complex_types[1:] ) )     
    checker.add ( vDoubles , vFloats, vInts , vLongs )
    if numpy : checker.add ( numpy.ndarray ) 

# =============================================================================
if '__main__' == __name__ :


    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    _v = [ l for l in dir(Ostap     ) if 0 != l.find('__') ]
    logger.info ('dir(Ostap)      : ')
    _v.sort()
    for v in _v : logger.info ( v )
    
    _v = [ l for l in dir(Ostap.Math) if 0 != l.find('__') ]
    logger.info ('dir(Ostap.Math) : ')
    _v.sort()
    for v in _v : logger.info ( v )

    if not numpy : logger.warning ( "Numpy module is not accesible!")
    if not scipy : logger.warning ( "Scipy module is not accesible!")
    
# =============================================================================
##                                                                     The  END
# =============================================================================
