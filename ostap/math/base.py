#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#
#  Simple file to provide "easy" access in python for
#  the basic ROOT::Math classes
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
               C++ namespaces: Gaudi & Gaudi::Math

  >>> import LHCbMath.Types
  >>> from GaudiPython.Bindings import gbl as cpp ## get global C++ namespace
  >>> Gaudi = cpp.Gaudi                           ## get C++ namespace Gaudi
  >>> p3 = Gaudi.XYZPoint(0,1,2)                  ## use C++ type Gaudi::XYZPoint

  >>> dir( Gaudi.Math )
  >>> dir( Gaudi      )

  Last modification $Date$
                 by $Author$

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
    ##
    'inrange'        , ## is double number in certain range?
    ## 
    'natural_number' , ## See Ostap::Math::natural_number  
    'natural_entry'  , ## See Ostap::Math::natural_entry
    ##
    'doubles'        , ## construct std::vector<double>
    'ints'           , ## construct std::vector<int>
    'uints'          , ## construct std::vector<unsigned int>
    'longs'          , ## construct std::vector<long>
    'ulongs'         , ## construct std::vector<unsigned long>
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
for t in ( 'int'                ,
           'long'               ,
           'long long'          ,
           'unsigned int'       ,
           'unsigned long'      ,
           'unsigned long long' ,
           'float'              ,
           'double'             ) :
    v = std.vector( t )
    v.asList   = lambda s :       [ i for i in s ]   ## convert vector into list
    v.toList   = v.asList
    v.__repr__ = lambda s : str ( [ i for i in s ] ) ## print it !
    v.__str__  = lambda s : str ( [ i for i in s ] ) ## print it !


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
def _add_to ( vct , arg1 , *args ) :
    ##
    if hasattr ( arg1 , '__iter__' ) :
        for a in arg1 : vct.push_back ( a ) 
    else : vct.push_back ( arg1 ) 
    #
    for a in args : _add_to ( vct , a )
        
# =============================================================================
## construct std::vector<double> from the arguments
def doubles ( arg1 , *args ) :
    """Construct the std::vector<double> from the arguments
    >>> v1 = doubles ( 1.01 )
    >>> v2 = doubles ( 1.01 , 1.02 , 1.03  )
    >>> v3 = doubles ( [ 1.01 , 1.02 , 1.03 ] )    
    """
    ## create new vector 
    VT  = std.vector('double')
    vct = VT()
    ## add arguments to the vector 
    _add_to ( vct , arg1 , *args )
    ## 
    return vct

# =============================================================================
## construct std::vector<ints> from the arguments
def ints ( arg1 , *args ) :
    """Construct the std::vector<int> from the arguments    
    >>> v1 = ints ( 1 )
    >>> v2 = ints ( 1 , 1 , 1  )
    >>> v3 = ints ( [ 1 , 2 , 3 ] )    
    """
    ## create new vector 
    VT  = std.vector('int')
    vct = VT()
    ## add arguments to the vector 
    _add_to ( vct , arg1 , *args )
    ## 
    return vct

# =============================================================================
## construct std::vector<unsigned int> from the arguments
def uints ( arg1 , *args ) :
    """Construct the std::vector<unsigned int> from the arguments    
    >>> v1 = uints ( 1 )
    >>> v2 = uints ( 1 , 1 , 1  )
    >>> v3 = uints ( [ 1 , 2 , 3 ] )    
    """
    ## create new vector 
    VT  = std.vector('unsigned int')
    vct = VT()
    ## add arguments to the vector 
    _add_to ( vct , arg1 , *args )
    ## 
    return vct

# =============================================================================
## construct std::vector<long> from the arguments
def longs ( arg1 , *args ) :
    """Construct the std::vector<long> from the arguments    
    >>> v1 = longs ( 1 )
    >>> v2 = longs ( 1 , 1 , 1  )
    >>> v3 = longs ( [ 1 , 2 , 3 ] )    
    """
    ## create new vector 
    VT  = std.vector('long')
    vct = VT()
    ## add arguments to the vector 
    _add_to ( vct , arg1 , *args )
    ## 
    return vct

# =============================================================================
## construct std::vector<unsigned long> from the arguments
def ulongs ( arg1 , *args ) :
    """Construct the std::vector<unsigned long> from the arguments    
    >>> v1 = ulongs ( 1 )
    >>> v2 = ulongs ( 1 , 1 , 1  )
    >>> v3 = ulongs ( [ 1 , 2 , 3 ] )    
    """
    ## create new vector 
    VT  = std.vector('unsigned long')
    vct = VT()
    ## add arguments to the vector 
    _add_to ( vct , arg1 , *args )
    ## 
    return vct

SPD = std.pair('double','double')
SPD.asTuple  = lambda s : (s.first,s.second)
SPD.__str__  = lambda s : str( (s.first,s.second) )
SPD.__repr__ = SPD.__str__

# =============================================================================
# Imporve operations with std.complex 
# =============================================================================
COMPLEX = cpp.std.complex('double')

def _cmplx_to_complex_ ( s ) :
    """convert to complex"""
    return  complex    ( s.real() , s.imag() )
def _cmplx_negate_     ( s ) :
    """Negation:
    >>> v  = ...
    >>> v1 = -v
    """
    return -complex    ( s.real() , s.imag() )
def _cmplx_abs_        ( s ) :
    """Absolute value
    >>> print abs(v) 
    """
    import math
    sr = s.real()
    si = s.imag()
    return math.sqrt( sr * sr + si * si ) 
def _cmplx_norm_       ( s ) :
    """Norm (squared absolute value)
    >>> print v.norm()
    """
    sr = s.real()
    si = s.imag()
    return sr * sr + si * si
def _cmplx_conjugate_  ( s ) :
    """Get complex conjugated
    >>> vc = v.conjugate() 
    """
    return complex     ( s.real() , -s.imag() )
    
def _cmplx_add_        ( s , o ) :
    """add complex values 
    >>> r = v + other  
    """
    return o + complex ( s.real() , s.imag() )
def _cmplx_mul_        ( s , o ) :
    """multiply  complex values 
    >>> r = v * other  
    """
    return o * complex ( s.real() , s.imag() )

def _cmplx_div_        ( s , o ) :
    return (1.0/o) * complex ( s.real() , s.imag() )
    """divide complex values 
    >>> r = v / other  
    """
def _cmplx_rdiv_       ( s , o ) :
    """divide complex values 
    >>> r = other / v 
    """
    return o       * ( 1.0 / complex ( s.real() , s.imag() ) )

def _cmplx_sub_        ( s , o ) :
    """subtract complex values 
    >>> r = v - other 
    """
    return (-o   ) + complex ( s.real() , s.imag() )
def _cmplx_rsub_       ( s , o ) :
    """subtract complex values 
    >>> r = other - v 
    """
    return   o     - complex ( s.real() , s.imag() )

def _cmplx_pow_  ( s , o ) :
    """power function 
    >>> r = v ** other  
    """
    if isinstance ( o , COMPLEX ) :
        o = complex ( o.real() , o.imag() ) 
    return complex ( s.real() , s.imag() ) ** o

def _cmplx_rpow_  ( s , o ) :
    """power function 
    >>> r = other **v  
    """
    return o ** complex ( s.real() , s.imag() )


def _cmplx_eq_    ( s , o ) :
    """equality:
    >>> r = v == other  
    """
    if isinstance ( o, COMPLEX ) :
        return s.real() == o.real() and s.imag() == o.imag()
    return complex( s.real() , s.imag() ) == o

def _cmplx_ne_    ( s , o ) :
    """non-equality:
    >>> r = v != other  
    """
    if isinstance ( o, COMPLEX ) :
        return s.real() != o.real() or  s.imag() != o.imag()
    return complex( s.real() , s.imag() ) != o 
    
COMPLEX.__complex__ = _cmplx_to_complex_

COMPLEX.__add__     = _cmplx_add_
COMPLEX.__mul__     = _cmplx_mul_
COMPLEX.__div__     = _cmplx_div_
COMPLEX.__sub__     = _cmplx_sub_

COMPLEX.__radd__    = _cmplx_add_
COMPLEX.__rmul__    = _cmplx_mul_
COMPLEX.__rdiv__    = _cmplx_rdiv_
COMPLEX.__rsub__    = _cmplx_rsub_

def _cmplx_iadd_ ( s , o ) :
    x = s + o
    s.real(x.real)
    s.imag(x.imag)
    
def _cmplx_isub_ ( s , o ) :
    x = s - o
    s.real(x.real)
    s.imag(x.imag)

def _cmplx_imul_ ( s , o ) :
    x = s * o
    s.real(x.real)
    s.imag(x.imag)

def _cmplx_idiv_ ( s , o ) :
    x = s / o
    s.real(x.real)
    s.imag(x.imag)

COMPLEX.__iadd__    = _cmplx_iadd_
COMPLEX.__imul__    = _cmplx_imul_
COMPLEX.__idiv__    = _cmplx_idiv_
COMPLEX.__isub__    = _cmplx_isub_

COMPLEX.__repr__    = lambda s : "%s" % complex ( s.real(), s.imag() )
COMPLEX.__str__     = lambda s : "%s" % complex ( s.real(), s.imag() )
COMPLEX.__abs__     = _cmplx_abs_
COMPLEX.__pow__     = _cmplx_pow_
COMPLEX.__rpow__    = _cmplx_rpow_
COMPLEX.__neg__     = _cmplx_negate_

COMPLEX.__eq__      =  _cmplx_eq_
COMPLEX.__ne__      =  _cmplx_ne_

COMPLEX.norm        = _cmplx_norm_
COMPLEX.conjugate   = _cmplx_conjugate_
COMPLEX.conj        = _cmplx_conjugate_
COMPLEX.to_complex  = _cmplx_to_complex_ 
COMPLEX.as_complex  = _cmplx_to_complex_ 


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
