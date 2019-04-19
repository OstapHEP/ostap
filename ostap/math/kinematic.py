#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/kinematic.py
#  Set of useful "kinematic" utilities 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2009-09-12
# =============================================================================
""" Set of useful ``kinematic'' utilities 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2009-09-12"
__version__ = "Version $Revision:$"
# =============================================================================
__all__     = (
    ##
    'EtaVsP'      , ## eta=eta(p)  for the fixed PT
    'EtaVsPt'     , ## eta=eta(pt) for the fixed P
    'YvsP'        , ## y=yy(p)     for the fixed PT and mass  
    'YvsPt'       , ## y=y(pt)     for the fixed P  and mass
    'EtaVsPPT'    , ## eta=eta(p,pt) 
    'PtVsPEta'    , ## pt=pt(p,eta) 
    'PvsPtEta'    , ## p= p(pt,eta)
    ##
    'kallen'      , ## Kallen ``lambda''/``triangle'' function
    'phasespace2' , ## 2-body phase space 
    'phasespace3' , ## the full 3-body phase space
    'phasespace4' , ## the full 4-body phase space
    'phasespace'  , ## the full N-body phase space
    ##
    'G'           , ## the basic universal 4-body function
    ##
    )
# =============================================================================
import ROOT, math
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.kinematic' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
from ostap.math.base  import cpp , COMPLEX 
from ostap.core.ostap_types import num_types

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
def _v4_pow_ ( self , e = 2 ) :
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
    _v4 . __mul__     = _v4_mul_
    _v4 . __add__     = _v4_add_
    _v4 . __sub__     = _v4_sub_
    _v4 . __div__     = _v4_div_    
    _v4 . __truediv__ = _v4_div_    
    _v4 . __radd__    = lambda s,o : s+o 
    _v4 . __rmul__    = lambda s,o : s*o 
    _v4 . __pow__     = _v4_pow_

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

if hasattr ( V4C , 'mag'         ) :  del V4C.mag
if hasattr ( V4C , 'mass'        ) :  del V4C.mass
if hasattr ( V4C , 'mt'          ) :  del V4C.mt
if hasattr ( V4C , 'r'           ) :  del V4C.r
if hasattr ( V4C , 'theta'       ) :  del V4C.theta
if hasattr ( V4C , 'phi'         ) :  del V4C.phi
if hasattr ( V4C , 'isSpaceLike' ) :  del V4C.isSpaceLike
if hasattr ( V4C , 'isTimeLike'  ) :  del V4C.isTimeLike
if hasattr ( V4C , 'isLightLike' ) :  del V4C.isLightLike
if hasattr ( V4C , 'M'           ) :  del V4C.M
if hasattr ( V4C , 'R'           ) :  del V4C.R
if hasattr ( V4C , 'Theta'       ) :  del V4C.Theta
if hasattr ( V4C , 'Phi'         ) :  del V4C.Phi
if hasattr ( V4C , 'Mt'          ) :  del V4C.Mt
if hasattr ( V4C , 'Et'          ) :  del V4C.Et

# =============================================================================
_acosh = math.acosh
_atanh = math.atanh
_sqrt  = math.sqrt
_cosh  = math.cosh
# =============================================================================
if not hasattr ( math , 'coth' ) :
    math.coth  = lambda x : 1.0/math.tanh (     x )
    math.coth. __doc__ = """coth(x)
    Return the   hyperbolic cotangent of x 
    """
    logger.debug ("Insert coth function into math")
# =============================================================================
if not hasattr ( math , 'acoth' ) :
    math.acoth = lambda x :     math.tanh ( 1.0/x )
    math.acoth. __doc__ = """acoth(x)
    Return the hyperbolic area cotangent of x (|x|>1) 
    """
    logger.debug ("Insert acoth function into math")

# =============================================================================

# =============================================================================
## @class EtaVsP
#  very simple function \f$ \eta = \eta(p) \f$  for the fixed transverse momentum
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-07-19
class EtaVsP(object) :
    """
    Very simple function \f$ \eta = \eta(p) \$  for the fixed transverse momentum
    """
    def __init__ ( self , pt ) :
        assert isinstance ( pt , num_types ) and 0<=pt , "PT is invalid!"
        self.pt = float ( pt )
    def __call__ ( self , p )  :
        return _acosh ( max ( p , self.pt )  / self.pt )

# =============================================================================
## @class EtaVsPt
#  very simple function \f$ \eta = \eta(p_T) \f$  for the fixed momentum
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-07-19
class EtaVsPt(object) :
    """Very simple function \f$ \eta = \eta(p_T) \$  for the fixed momentum"""
    def __init__ ( self , p  ) :
        assert isinstance ( p , num_types ) and 0<=p , "P is invalid!"
        self.p = float ( p )
    def __call__ ( self , pt )  :
        return _acosh ( self.p / min ( pt , self.p ) ) 

# =============================================================================
## rapidity as function of P for fixed pt and mass 
class YvsP(object) :
    """Rapidity as function of P for fixed pt and mass"""
    def __init__ ( self , pt , mass ) :
        assert isinstance ( pt   , num_types ) and 0<=pt   , "PT is invalid!"
        assert isinstance ( mass , num_types ) and 0<=mass , "M  is invalid!"
        self.pt2 = float ( pt   ) * float ( pt   ) 
        self.m2  = float ( mass ) * float ( mass )
        
    def __call__ ( self , p ) :
        
        p2  = p * p
        e2  = p2 + self.m2  
        pz2 = p2 - self.pt2 
        
        return _atanh ( _sqrt ( max ( pz2 , 0.0 ) / e2  ) )
        
# =============================================================================
## rapidity as function of Pt for fixed p and mass 
class YvsPt(object) :
    """Rapidity as function of Pt for fixed p and mass"""
    def __init__ ( self , p , mass ) :
        assert isinstance ( p    , num_types ) and 0<=p    , "P  is invalid!"
        assert isinstance ( mass , num_types ) and 0<=mass , "M  is invalid!"
        self.p2  = float ( p    ) * float ( p    ) 
        self.m2  = float ( mass ) * float ( mass )
        self.e2  = self.p2 + self.m2
    def __call__ ( self , pt ) :
        pt2 = pt * pt
        pz2 = self.p2 - pt2 
        return _atanh ( _sqrt ( max ( pz2 , 0.0 ) / self.e2  ) )

# =============================================================================
## helper wrapper 
class _WF1(object) :
    def __init__ ( self , obj           ) :        self.obj = obj 
    def __call__ ( self , x , pars = [] ) : return self.obj ( x[0] ) 

# =============================================================================
## convert the objects to the functions 
def _as_TF1_ ( obj , xmin , xmax ) :
    """Convert the objects to the functions"""

    from  ostap.core.core import funID
    
    fobj     = _WF1( obj ) 
    fun      = ROOT.TF1( funID() , fobj , xmin , xmax ) 
    fun._obj = fobj
    
    fun.SetNpx(500)
    
    return fun

EtaVsP  . asTF1 = _as_TF1_
EtaVsPt . asTF1 = _as_TF1_
YvsP    . asTF1 = _as_TF1_
YvsPt   . asTF1 = _as_TF1_

# =============================================================================
## eta = eta(p ,pt )"
class EtaVsPPT(object) :
    "eta = eta(p ,pt )"
    def __call__ ( self , p   , pt  ) : return _acosh ( p / pt ) 

# =============================================================================
## pt  = pt (p ,eta)
class PtVsPEta(object) :
    "pt  = pt (p ,eta)"
    def __call__ ( self , p   , eta ) : return p  / _cosh ( eta )

# =============================================================================
## p   = p  (pt,eta)
class PvsPtEta(object) :
    "p   = p  (pt,eta)"
    def __call__ ( self , pt  , eta ) : return pt * _cosh ( eta ) 

# =============================================================================
## helper wrapper 
class _WF2(object) :
    def __init__ ( self , obj           ) :        self.obj = obj 
    def __call__ ( self , x , pars = [] ) : return self.obj ( x[0] , x[1] ) 

# =============================================================================
## convert the objects to the function 
def _as_TF2_ ( obj , xmin , xmax , ymin , ymax ) :
    """Convert the objects to the function"""

    from ostap.core.core import funID

    fobj       = _WF2(obj) 
    fun        = ROOT.TF2( funID() , fobj , xmin , xmax , ymin , ymax ) 
    fun._obj   = fobj

    fun.SetNpx(250)
    fun.SetNpy(250)
    
    return fun

EtaVsPPT  . asTF2 = _as_TF2_
PtVsPEta  . asTF2 = _as_TF2_
PvsPtEta  . asTF2 = _as_TF2_

# =============================================================================
## Kallen function, aka ``lambda''/``triangle'' function 
#  @see https://en.wikipedia.org/wiki/K%C3%A4ll%C3%A9n_function
def kallen ( x , y , z ) :
    """ Kallen function, aka ``triangle'' function 
    - see https://en.wikipedia.org/wiki/K%C3%A4ll%C3%A9n_function
    """
    return x * x + y * y + z * z - 2.0 * x * y - 2.0 * y * z -  2.0 * z * x

# =============================================================================
## Calculate the two-body phase space
#  \f$ R_2 = \frac{ \pi \lambda^{1/2}( M^2 , m_1^2 , m_2^2) }{2*M^2} \f$ 
#  @code
#  M, m1  , m2 = ...
#  ps = phasespace2 ( M , m1 , m2 ) 
#  @endcode 
def phasespace2 ( M ,  m1 , m2 ) :
    r"""Calculate the full two-body phase space
    \f$ R_2 = \frac{ \pi \lambda^{1/2}( M^2 , m_1^2 , m_2^2) }{2*M^2} \f$ 
    >>> M, m1  , m2 = ...
    >>> ps = phasespace2 ( M , m1 , m2 ) 
    """    
    assert 0<M and 0<=m1 and 0<=m2, 'Invalid setting of masses!'
    
    ##
    if m1 + m2 >= M : return 0   ## RETURN!

    s =  M * M 
    import math
    return math.pi * math.sqrt ( kallen ( s , m1 * m1 , m2 * m2 ) ) / ( 2.0 * s ) 

# =============================================================================
## Calculate the three body phase space 
#  @code
#  M, m1  , m2 , m3 = ...
#  ps3 = phasespace3 ( M , m1  , m2 , m3 ) 
#  @endcode 
def phasespace3 ( M ,  m1 , m2 , m3 ) :
    """Calculate the full three body phase space:
    >>> M, m1  , m2 , m3 = ...
    >>> ps3 = phasespace3 ( M , m1  , m2 , m3 ) 
    """
    assert 0<M and 0<=m1 and 0<=m2 and 0<=m3 , 'Invalid setting of masses!'

    ##
    if m1 + m2 + m3 >= M : return 0   ## RETURN! 

    s    =  M * M
    m1_2 = m1 * m1 
    m2_2 = m2 * m2
    m3_2 = m3 * m3
    
    high = ( M  - m1 ) ** 2
    low  = ( m2 + m3 ) ** 2
    
    import math
    func = lambda x : math.sqrt ( kallen ( x , s     , m1_2 ) *
                                  kallen ( x , m2_2  , m3_2 ) ) / x
    
    from ostap.math.integral import integral 
    
    r = integral ( func ,  low , high , err = False )
    
    return ( math.pi**2 ) * r  / ( 4.0 * s ) 


# =============================================================================
## Calculate the four body phase space 
#  @code
#  M, m1  , m2 , m3 , m4 = ...
#  ps4 = phasespace4 ( M , m1  , m2 , m3 , m4 ) 
#  @endcode
#  The algorithm includes two embedded numerical integration -> could be relatively slow 
def phasespace4 ( M ,  m1 , m2 , m3 , m4 ) :
    """Calculate the full four body phase space
    >>> M, m1  , m2 , m3 , m4 = ...
    >>> ps4 = phasespace4 ( M , m1  , m2 , m3 , m4 ) 
    - The algorithm includes two embedded numerical integration -> could be relatively slow 
    """
    assert 0<M and 0<=m1 and 0<=m2 and 0<=m3 and 0<=m4 , 'Invalid setting of masses!'

    ##
    if m1 + m2 + m3 + m4 >= M : return 0   ## RETURN! 
    
    low  = m1 + m2
    high = M  - m3 - m4  
    
    func = lambda x : 2.0 * x * phasespace3 ( M , x , m3 , m4 ) * phasespace2 ( x , m1 , m2  )

    from ostap.math.integral import integral 
    
    return integral ( func ,  low , high , err = False )

# ==============================================================================
## calculate full N-body phase space
#  @code
#  M, m1 , m2 , ... , mn = ...
#  ps = phasespace ( M , m1 , m2 , ... , mn ) 
#  @endcode 
#  The algorithm includes embedded numerical integrations -> could be relatively slow 
def phasespace ( M , *args ) :
    """Calculate full  N-body phase space
    >>> M, m1 , m2 , ... , mn = ...
    >>> ps = phasespace ( M , m1 , m2 , ... , mn )
    - The algorithm includes embedded numerical integrations -> could be relatively slow 
    """
    
    assert 0 < M and 2 <= len ( args ) , 'Invalid setting of masses!'
    
    summ = 0.0
    for m in args :
        assert 0 <= m , 'Invalid setting of masses'
        summ += m
        
    if summ >= M : return 0
    
    N = len ( args )
    if   2 == N : return phasespace2 ( M , *args )
    elif 3 == N : return phasespace3 ( M , *args )
    elif 4 == N : return phasespace4 ( M , *args )

    ## split particles into two groups & (recursively) apply the splitting formula
    
    k  = N/2

    args1 = args[k:]
    args2 = args[:k]

    low  =     sum ( args1 )
    high = M - sum ( args2 ) 

    func = lambda x : 2.0 * x * phasespace ( M , x , *args2 ) *  phasespace ( x , *args1 )
    
    from ostap.math.integral import integral 
    
    return integral ( func ,  low , high , err = False )


# =============================================================================
## the basic universal 4-body function G
#  @see E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
#       London, New York, Sydney, Toronto, 1973, p.89, eq. (5.23)
#  @see https://userweb.jlab.org/~rafopar/Book/byckling_kajantie.pdf
#  E.g. physical range for 2->2 scattering process is defined as
#  \f$ G(s,t,m_2^2, m_a^2, m_b^2, m_1^2) \le 0 \f$
# or the phsyical range  for Dalitz plot is
#   \f$ G(s_2, s_1,  m_3^2, m_1^2, s , m_2^2) \le 0 \f$ 
def G ( x , y , z , u , v , w ) :
    r"""The basic universal 4-body function G
    - see E.Byckling, K.Kajantie, ``Particle kinematics'', John Wiley & Sons,
    London, New York, Sydney, Toronto, 1973 p.89, eq. (5.23)
    - see https://userweb.jlab.org/~rafopar/Book/byckling_kajantie.pdf
    
    E.g. physical range for 2->2 scattering process is defined as
    \f$ G(s,t,m_2^2, m_a^2, m_b^2 , m_1^2)     \le 0 \f$
    or the physical range  for Dalitz plot is
    \f$ G(s_2, s_1,  m_3^2, m_1^2 , s , m_2^2) \le 0 \f$ 
    """
    r1 = x * x * y + x * y * y + z * z * x + z * u * u + v * v * w + v *  w * w 
    r2 = x * y * w + x * u * v + y * z * w + y * u * w
    r3 = - x * y * ( z + u + v + w )
    r4 = - z * y * ( x + y + v + w )
    r5 = - v * w * ( x + y + z + u )
    
    return 0.0 + r1 + r2 + r3 + r4 + r5 
    

    
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
