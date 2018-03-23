#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Set of useful "kinematic" utilities 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2009-09-12
# =============================================================================
""" Set of useful ``kinematic'' utilities 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2009-09-12"
__version__ = "Version$Revision$"
# =============================================================================
__all__     = ()
# =============================================================================
import ROOT
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.kinematic' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
from ostap.math.base import cpp , COMPLEX 

## C++ namespace Ostap 
Ostap = cpp.Ostap

## ROOT::Math namespace
_RM = ROOT.ROOT.Math
Ostap.LorentzVector        = _RM.LorentzVector ('ROOT::Math::PxPyPzE4D<double>')
Ostap.ComplexLorentzVector = _RM.LorentzVector ('ROOT::Math::PxPyPzE4D<std::complex<double> >')

Ostap.Math.LorentzVector         = Ostap.LorentzVector
Ostap.Math.ComplexLorentzVector  = Ostap.ComplexLorentzVector


V4D = Ostap.LorentzVector 
V4C = Ostap.ComplexLorentzVector

# =============================================================================
def _v4c_init_ ( self , *args ) :
    """Construct  complex lorenzt vector from another vector or from coordinates
    >>> v4  = ...
    >>> lv  = ComplexLorenztVector ( v4 )
    >>> lv1 = ComplexLorenztVector () 
    >>> lv2 = ComplexLorenztVector ( 1    ) 
    >>> lv3 = ComplexLorenztVector ( 1+2j ) 
    >>> lv4 = ComplexLorenztVector ( 1,  2, 3, 4 ) 
    >>> lv5 = ComplexLorenztVector ( 1,  2, 3, 4+3j ) 
    """
    if args and 1 == len ( args ) :
        a = args[0]
        if   isinstance ( a , V4C ) : return V4C._old_init_ ( self , a ) 
        elif isinstance ( a , V4D ) :
            return V4C._old_init_ ( self ,
                                    COMPLEX( a.X() ) ,
                                    COMPLEX( a.Y() ) ,
                                    COMPLEX( a.Z() ) ,
                                    COMPLEX( a.T() ) )

    assert len( args ) <= 4 ,'Invalid  argument length!'
    ac = [ COMPLEX ( i ) for i in args ]
    while len ( ac ) < 4 : ac.append ( COMPLEX() ) 
    return V4C._old_init_ ( self , *ac )

if not  hasattr ( V4C , '_old_init_' ) :
    V4C. _old_init_ =  V4C. __init__
    V4C. __init__   = _v4c_init_ 
    
## ============================================================================
## 4-vectors 
def _v4_iadd_ ( s , other ) :
    """Increment 4-vector with other 4-vector
    >>> p4  = ...
    >>> p4 += other 
    """
    s.SetE    ( s.E  () + other.E  () )
    s.SetPx   ( s.Px () + other.Px () )
    s.SetPy   ( s.Py () + other.Py () )
    s.SetPz   ( s.Pz () + other.Pz () )
    return s

## ============================================================================
def _v4_isub_ ( s , other ) :
    """Decrement 4-vector with other 4-vector
    >>> p4  = ...
    >>> p4 -= other 
    """
    s.SetE    ( s.E  () - other.E  () )
    s.SetPx   ( s.Px () - other.Px () )
    s.SetPy   ( s.Py () - other.Py () )
    s.SetPz   ( s.Pz () - other.Pz () )
    return s

## ============================================================================
## 4-vectors 
def _v4c_iadd_ ( s , other ) :
    """Increment 4-vector with other 4-vector
    >>> p4  = ...
    >>> p4 += other 
    """
    s.SetE    ( COMPLEX ( s.E  () + other.E  () ) )
    s.SetPx   ( COMPLEX ( s.Px () + other.Px () ) ) 
    s.SetPy   ( COMPLEX ( s.Py () + other.Py () ) ) 
    s.SetPz   ( COMPLEX ( s.Pz () + other.Pz () ) ) 
    return s

## ============================================================================
def _v4c_isub_ ( s , other ) :
    """Decrement 4-vector with other 4-vector
    >>> p4  = ...
    >>> p4 -= other 
    """
    s.SetE    ( COMPLEX ( s.E  () - other.E  () ) )
    s.SetPx   ( COMPLEX ( s.Px () - other.Px () ) ) 
    s.SetPy   ( COMPLEX ( s.Py () - other.Py () ) ) 
    s.SetPz   ( COMPLEX ( s.Pz () - other.Pz () ) ) 
    return s

## ============================================================================
def _v4_dot_   ( s , other ) :
    """``Dot''-prodcut of two 4-vectors 
    >>> p4    = ...
    >>> other = ...
    >>> print 'Q2 is ', p4.Dot ( other )
    """
    res  = s.e  () * other.e  () 
    res -= s.px () * other.px ()
    res -= s.py () * other.py ()
    res -= s.pz () * other.pz ()
    return res 


if not hasattr ( V4D , '__iadd__' ) : V4D. __iadd__ = _v4_iadd_ 
if not hasattr ( V4D , '__isub__' ) : V4D. __isub__ = _v4_isub_ 
if not hasattr ( V4D , 'Dot'      ) : V4D.Dot       = _v4_dot_
if not hasattr ( V4C , '__iadd__' ) : V4C. __iadd__ = _v4c_iadd_ 
if not hasattr ( V4C , '__isub__' ) : V4C. __isub__ = _v4c_isub_ 
if not hasattr ( V4C , 'Dot'      ) : V4C.Dot       = _v4_dot_


## ============================================================================
def _v4_mul_ ( self , other ) :
    """Multiplication/scaling of Lorentz Vectors 
    >>> vct = ...
    >>> a   = vct * 2
    
    >>> vct2 = ...
    >>> prod = vct * vct2 ## NB! 
    """
    if isinstance ( other , ( V4D , V4C ) ) : return self.Dot ( other )
    # 
    # scaling by the number :
    tmp   = self.__class__ ( self )
    tmp  *= other
    return tmp

## ============================================================================
def _v4_add_ ( self , other ) :
    """ Addition of Lorentz Vectors 
    >>> vct1 = ...
    >>> vct2 = ...
    >>> a    = vct1 + vct2
    """
    tmp   = self.__class__ ( self )
    tmp  += other
    return tmp

## ============================================================================
def _v4_sub_ ( self , other ) :
    """Subtraction of Lorentz Vectors 
    >>> vct1 = ...
    >>> vct2 = ...
    >>> a    = vct1 - vct2
    """
    tmp   = self.__class__ ( self )
    tmp  -= other
    return tmp

## ============================================================================
def _v4_div_ ( self , other ) :
    """Division/scaling of Lorentz Vectors     
    >>> vct = ...
    >>> a   = vct / 2 
    """
    tmp   = self.__class__ ( self )
    tmp  /= other
    return tmp

# =============================================================================
def _v4_pow_ ( self , e ) :
    """Squared length of the Lorentz vector
    >>> print ' mass-squared is:', p4**2 
    """
    if 2 != e : return NotImplemented
    return self.M2   ()

# =============================================================================
## Self-printout of 4D-vectors
def _v4d_str_ ( self , fmt = "[(%g,%g,%g),%g]" ) :
    """Self-printout of 4D-vectors
    >>> print p4 
    """
    return fmt % ( self.X() , self.Y( ), self.Z() , self.E() )
# =============================================================================
## Self-printout of 4D-vectors
def _v4c_str_ ( self , fmt = "[(%s,%s,%s),%s]" ) :
    """Self-printout of 4D-vectors
    >>> print p4 
    """
    return fmt % ( self.X() , self.Y( ), self.Z() , self.E() )

for _v4 in ( V4D , V4C ) : 
    _v4 . __mul__   = _v4_mul_
    _v4 . __add__   = _v4_add_
    _v4 . __sub__   = _v4_sub_
    _v4 . __div__   = _v4_div_    
    _v4 . __radd__  = lambda s,o : s+o 
    _v4 . __rmul__  = lambda s,o : s*o 
    _v4 . __pow__   = _v4_pow_

# =============================================================================
if not hasattr ( V4D , '_new_str_' ) :
    V4D . _new_str_ = _v4d_str_
    V4D . __str__   = _v4d_str_
    V4D . __repr__  = _v4d_str_
    
# =============================================================================
if not hasattr ( V4C , '_new_str_' ) :
    V4C . _new_str_ = _v4c_str_
    V4C . __str__   = _v4c_str_
    V4C . __repr__  = _v4c_str_

# =============================================================================
def _v4c_conjugate_ ( s ) :
    """Complex conjugation for 4-vector
    >>> lv = ...
    >>> cc = lv.conjugate()
    """
    return V4C ( s.X ().cpp_conj () ,
                 s.Y ().cpp_conj () ,
                 s.Z ().cpp_conj () ,
                 s.T ().cpp_conj () )

V4C.conj      = _v4c_conjugate_
V4C.conjugate = _v4c_conjugate_
V4D.conj      = lambda s : V4D ( s ) 
V4D.conjugate = lambda s : V4D ( s )  

# =============================================================================
if not hasattr ( Ostap.Math , 'Vector4' ) :
    import ostap.math.linalg
    
_V4 = Ostap.Math.Vector4
# =============================================================================
## convert LorentzVector into SVector
#  @code
#  lv = ...
#  v4 = lv.asSVector()
#  @endcode 
def _v4_as_v4_ ( self ) :
    """Convert LorentzVector into SVector
    >>> lv = ...
    >>> v4 = lv.asSVector()
    """
    v4 = _V4()
    v4[0] = self.X()
    v4[1] = self.Y()
    v4[2] = self.Z()
    v4[3] = self.T()
    return _v4 

V4D.asSVector = _v4_as_v4_ 

# =============================================================================
__euclidianNorm2 = Ostap.Math.euclidianNorm2 
__restMomentum   = Ostap.Math.restMomentum
__restEnergy     = Ostap.Math.restEnergy 
__boost          = Ostap.Math.boost
# =============================================================================
def _v4_en2_ ( v ) :
    """Get ``euclidian  norm squared'' of 4-vector
    >>>  v4 = ...
    >>>  print v4.euclidianNorm2()
    """
    return __euclidianNorm2 ( v )
# =============================================================================
def _v4_restE_ ( v , M ) :
    """Get the value of the energy in the rest-frame of particle 'M'
    >>>  v4 = ...
    >>>  M  = ... 
    >>>  print v4.restEnergy ( M )
    """
    m2 = M.M2()
    assert 0 < m2 , "RestEnergy: the reference system is not time-like: %s/%s"   % ( M , m2 ) 
    return __restEnergy ( v , M )
# =============================================================================
def _v4_restP_ ( v , M ) :
    """Get the value of the momentum in the rest-frame of particle 'M'
    >>>  v4 = ...
    >>>  M  = ... 
    >>>  print v4.restMomentum ( M )
    """
    m2 = M.M2()
    assert 0 < m2 , "RestMomentum: the reference system is not time-like: %s/%s"   % ( M , m2 ) 
    return __restMomentum ( v , M )
# =============================================================================
def _v4_boost_ ( v , M ) :
    """Boost the momentum into rest-frame of particle 'M'
    >>>  v4 = ...
    >>>  M  = ... 
    >>>  print v4.boost   ( M )
    >>>  print v4.boosted ( M )  ## ditto 
    """
    m2 = M.M2()
    assert 0 < m2 , "boost: the reference system is not time-like: %s/%s"   % ( M , m2 ) 
    return __boost ( v , M )

if not hasattr ( V4D , 'euclidianNorm2' ) : V4D.euclidianNorm2 = _v4_en2_ 
if not hasattr ( V4D , 'restEnergy'     ) : V4D.restEnergy     = _v4_restE_ 
if not hasattr ( V4D , 'restMomentum'   ) : V4D.restMomentum   = _v4_restP_ 
if not hasattr ( V4D , 'boost'          ) : V4D.boost          = _v4_boost_ 
if not hasattr ( V4D , 'boosted'        ) : V4D.boosted        = _v4_boost_ 



# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    lv0 = Ostap.LorentzVector        (2,1,0,3)

    logger.info ( 'Lorentz-vector ((px,py,pz),E)    : %s' % lv0 )

    CLV = Ostap.ComplexLorentzVector
    
    lv1 = CLV ()
    lv2 = CLV ( lv0 )
    lv3 = CLV ( 1+2j ) 
    lv4 = CLV ( 1,  2, 3, 4 ) 
    lv5 = CLV ( 1,  2, 3, 4+3j ) 

    logger.info ( 'Lorentz-vector ((px,py,pz),E)    : %s' % lv1 )
    logger.info ( 'Lorentz-vector ((px,py,pz),E)    : %s' % lv2 )
    logger.info ( 'Lorentz-vector ((px,py,pz),E)    : %s' % lv3 )
    logger.info ( 'Lorentz-vector ((px,py,pz),E)    : %s' % lv4 )
    logger.info ( 'Lorentz-vector ((px,py,pz),E)    : %s' % lv5 )


    logger.info ( 'Plus     %s' % (lv3+lv2) )
    logger.info ( 'Minus    %s' % (lv3-lv2) )
    logger.info ( 'Multiply %s' % (lv3*lv2) )
    logger.info ( 'Pow      %s' % (lv3**2 ) )
    
    
    logger.info ( 80*'*' ) 

# =============================================================================
# The END 
# =============================================================================
