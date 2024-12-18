#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/math_ve.py
#  Simple file to provide some wrapper function for dealing with
#  Gaudi::Math::ValueWithError objects 
#  @see Ostap::Math::ValueWithError
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-06-02
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-02"
__version__ = ""
# =============================================================================
__all__     = (
    'exp'          , 'expm1'      ,
    'log'          , 'log10'      , 'log1p'   , 
    'sqrt'         , 'cbrt'       , 'pow'     ,   
    'sin'          , 'cos'        , 'tan'     ,
    'cot'          , 'sec'        , 'csc'     , 'cas'      ,
    'sinh'         , 'cosh'       , 'tanh'    , 'sech'     ,
    'asin'         , 'acos'       , 'atan'    , 'atan2'    , 
    'asinh'        , 'acosh'      , 'atanh'   ,
    'erf'          , 'erfc'       , 'erfi'    , 'erfcx'    ,
    'sinc'         , 
    'probit'       , 'pochhammer' , 
    'gamma'        , 'tgamma'     , 'lgamma'  , 'igamma'   ,
    'psi'          , 'polygamma'  , 'digamma' , 'trigamma' ,
    'beta'         , 'lnbeta'     ,
    'gd'           , 'gd_inv'     ,
    'Ai'           , 'Bi'         ,
    'li'           , 'Li'         ,           ## integral logarithms 
    'Ei'           , 'Ein'        , 'E1'    , ## integral exponents 
    'Si'           , 'Ci'         , 'Cin'   , ## integral trigonometry
    'exp2'         , 'log2'       ,
    'bessel_J'     , 'bessel_Y'   ,           ## Bessel functions  
    'bessel_I'     , 'bessel_K'   ,           ## Bessel functions 
    'gauss_pdf'    , 'gauss_cdf'  ,           ## Gaussian PDF and CDF 
    'hypot'        , 'sigmoid'    , 'fma'   ,
    'minv'         , 'maxv'       ,
    'agm'          , 'ghm'        ,
    'bring'        , 'BR'         ,
    'am'           , 'dn'         ,            ## elliptic functions 
    'sn'           , 'cn'         , 'sc'     , ## elliptic functions
    'elliptic_K'   , 'elliptic_K' ,            ## elliptic integrals 
    'significance' , 'nsigmas'    , 'nsigma'
    )
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.math_ve' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
from   ostap.math.ve          import VE
from   ostap.math.base        import Ostap, iszero, isequal, complex_types 
from   ostap.core.ostap_types import num_types, is_integer, integer_types  
import ROOT, math, cmath 
# =============================================================================
_ln2_i = 1/math.log(2.0)                 ## useful constant 
# =============================================================================

# =============================================================================
## define `exp` function 
def exp ( x ) :
    """ `exp` function taking into account the uncertainties
    """
    fun = getattr ( x , '__exp__' , None )
    if fun : return fun()
    if isinstance ( x , complex_types ) : return cmath.exp ( x )
    return math.exp ( x )

# =============================================================================
## define `exp2` function 
def exp2 ( x ) :
    """ `exp2` function taking into account the uncertainties
    """
    fun = getattr ( x , '__exp2__' , None )
    if fun : return fun()
    return 2**x 

# =============================================================================
## define `expm1` function 
def expm1 ( x ) :
    """ `expm1` function taking into account the uncertainties
    """
    fun = getattr ( x , '__expm1__' , None )
    if fun : return fun()
    return math.expm1 ( x )

# =============================================================================
## define `log` function 
def log ( x ) :
    """ `log` function taking into account the uncertainties
    """
    fun = getattr ( x , '__log__' , None )
    if fun : return fun()
    if isinstance ( x , complex_types ) : return cmath.log ( x )
    return math.log ( x )

# =============================================================================
## define `log2' function 
def log2 ( x ) :
    """ `log2` function taking into account the uncertainties
    """
    fun = getattr ( x , '__log2__' , None )
    if fun : return fun()
    ## 
    return _ln2_i * math.log ( x ) 

# =============================================================================
## define `log10' function 
def log10 ( x ) :
    """ `log10` function taking into account the uncertainties
    """
    fun = getattr ( x , '__log10__' , None )
    if fun : return fun()
    if isinstance ( x , complex_types ) : return cmath.log10 ( x )
    return math.log10 ( x )

# =============================================================================
## define `log1p` function 
def log1p ( x ) :
    """ `log1p` function taking into account the uncertainties
    """
    fun = getattr ( x , '__log1p__' , None )
    if fun : return fun()
    return math.log1p ( x )

# =============================================================================
## define `sqrt` function 
def sqrt ( x ) :
    """ `sqrt` function taking into account the uncertainties
    """
    fun = getattr ( x , '__sqrt__' , None )
    if fun : return fun()
    if isinstance ( x , complex_types ) : return cmath.sqrt ( x )
    return math.sqrt ( x )

# =============================================================================
## define `cbrt` function 
def cbrt ( x ) :
    """ `cbrt` function taking into account the uncertainties
    """
    fun = getattr ( x , '__cbrt__' , None )
    if fun : return fun()
    return math.pow ( x , 1.0/3.0 )

# =============================================================================
## define `pow' function 
def pow ( x , y , *a ) :
    """ `pow` function taking into account the uncertainties
    """
    if   isinstance ( x , VE ) or isinstance ( y , VE ) : return x**y 
    return math.pow ( x , y , *a ) 
 
# =============================================================================
## define `sin` function 
def sin ( x ) :
    """ `sin` function taking into account the uncertainties
    """
    fun = getattr ( x , '__sin__' , None )
    if fun : return fun()
    if isinstance ( x , complex_types ) : return cmath.sin ( x )
    return math.sin ( x )

# =============================================================================
## define `cos` function 
def cos ( x ) :
    """ `cos` function taking into account the uncertainties
    """
    fun = getattr ( x , '__cos__' , None )
    if fun : return fun()
    if isinstance ( x , complex_types ) : return cmath.cos ( x )
    return math.cos ( x )

# =============================================================================
## define `tan` function 
def tan ( x ) :
    """ `tangent` function taking into account the uncertainties
    """
    fun = getattr ( x , '__tan__' , None )
    if fun : return fun()
    if isinstance ( x , complex_types ) : return cmath.tan ( x )
    return math.tan ( x )

# =============================================================================
## define `cot` function 
def cot ( x ) :
    """ `tangent` function taking into account the uncertainties
    """
    fun = getattr ( x , '__cot__' , None )
    if fun : return fun()
    return 1.0/tan ( x )

# =============================================================================
## define `sec` function 
def sec ( x ) :
    """ `secant` function, 1/cos,  taking into account the uncertainties
    """
    fun = getattr ( x , '__sec__' , None )
    if fun : return fun()
    return 1.0/cos ( x )

# =============================================================================
## define `csc` function 
def csc ( x ) :
    """ `cosecant` function, 1/sin,  taking into account the uncertainties
    """
    fun = getattr ( x , '__csc__' , None )
    if fun : return fun()
    return 1.0/cos ( x )

# =============================================================================
## define `cas` function 
def cas ( x ) :
    """ `cas` function: sin(x) + cos(s), ,  taking into account the uncertainties
    """
    fun = getattr ( x , '__cas__' , None )
    if fun : return fun()
    return sin ( x ) + cos ( x ) 

# =============================================================================
## define `sinh` function 
def sinh ( x ) :
    """ `sinh`function taking into account the uncertainties
    """
    fun = getattr ( x , '__sinh__' , None )
    if fun : return fun()
    if isinstance ( x , complex_types ) : return cmath.sinh ( x )
    return math.sinh ( x )

# =============================================================================
## define `cosh` function 
def cosh ( x ) :
    """ `cosh` function taking into account the uncertainties
    """
    fun = getattr ( x , '__cosh__' , None )
    if fun : return fun()
    if isinstance ( x , complex_types ) : return cmath.cosh ( x )
    return math.cosh ( x )

# =============================================================================
## define `tanh` function 
def tanh ( x ) :
    """ `tanh' function taking into account the uncertainties
    """
    fun = getattr ( x , '__tanh__' , None )
    if fun : return fun()
    if isinstance ( x , complex_types ) : return cmath.tanh ( x )
    return math.tanh ( x )

_erf_   = Ostap.Math.erf 
_erfc_  = Ostap.Math.erfc
_erfi_  = Ostap.Math.erfi
_erfcx_ = Ostap.Math.erfcx

# =============================================================================
## define `erf` function 
#  @see https://en.wikipedia.org/wiki/Error_function
def erf ( x ) :
    """ Error function taking into account the uncertainties
    - see https://en.wikipedia.org/wiki/Error_function
    """
    fun = getattr ( x , '__erf__' , None )
    if fun : return fun()
    return _erf_ ( x )

# =============================================================================
## define `erfc' function : complementary error function 
#  @see https://en.wikipedia.org/wiki/Error_function
def erfc ( x ) :
    """ Complemenatry error function taking into account the uncertainties
    - see https://en.wikipedia.org/wiki/Error_function
    """
    fun = getattr ( x , '__erfc__' , None )
    if fun : return fun()
    return _erfc_( x )

# =============================================================================
## define `erfcx' function 
#  @see https://en.wikipedia.org/wiki/Error_function
def erfcx ( x ) :
    """ Complemenatry scaled error function taking into account the uncertainties
    - see https://en.wikipedia.org/wiki/Error_function
    """
    fun = getattr ( x , '__erfcx__' , None )
    if fun : return fun()
    return _erfcx_( x )

# =============================================================================
## define `erfi` function 
#  @see https://en.wikipedia.org/wiki/Error_function
def erfi ( x ) :
    """ Imaginary error function taking into account the uncertainties
    - see https://en.wikipedia.org/wiki/Error_function
    """
    fun = getattr ( x , '__erfi__' , None )
    if fun : return fun()
    return _erfi_( x )


# =============================================================================
## define `asin` function 
def asin ( x ) :
    """ `asin` function taking into account the uncertainties
    """
    fun = getattr ( x , '__asin__' , None )
    if fun : return fun()
    return math.asin ( x )

# =============================================================================
## define `acos` function 
def acos ( x ) :
    """ `acos` function taking into account the uncertainties
    """
    fun = getattr ( x , '__acos__' , None )
    if fun : return fun()
    return math.acos ( x )

# =============================================================================
## define `atan` function 
def atan ( x ) :
    """ `atan` function taking into account the uncertainties
    """
    fun = getattr ( x , '__atan__' , None )
    if fun : return fun()
    return math.atan ( x )

# =============================================================================
## define `atan2` function 
def atan2 ( x , b = 1 ) :
    """ `atan2` function
    """
    fun = getattr ( x , '__atan2__' , None )
    if fun : return fun ( b )
    return math.atan2 ( x , b )


# =============================================================================
## define `asinh` function 
def asinh ( x ) :
    """ `asinh` function taking into account the uncertainties
    """
    fun = getattr ( x , '__asinh__' , None )
    if fun : return fun()
    return math.asinh ( x )

# =============================================================================
## define `acosh` function 
def acosh ( x ) :
    """ `acosh` function taking into account the uncertainties
    """
    fun = getattr ( x , '__acosh__' , None )
    if fun : return fun()
    return math.acosh ( x )

# =============================================================================
## define `atanh` function 
def atanh ( x ) :
    """ `atanh` function taking into account the uncertainties
    """
    fun = getattr ( x , '__atanh__' , None )
    if fun : return fun()
    return math.atanh ( x )


# =============================================================================
## define `tgamma` function 
def tgamma ( x ) :
    """ `tgamma` function taking into account the uncertainties
    """
    fun = getattr ( x , '__tgamma__' , None )
    if fun : return fun()
    return math.gamma ( x )

# =============================================================================
## define `gamma' function 
gamma = tgamma

# =============================================================================
## define `lgamma` function 
def lgamma ( x ) :
    """ `lgamma` function taking into account the uncertainties
    """
    fun = getattr ( x , '__lgamma__' , None )
    if fun : return fun()
    return math.lgamma ( x )

# =============================================================================
_igamma_ = Ostap.Math.igamma

# =============================================================================
## define `igamma` function
#  \f$ f(x) = \frac{1}{\Gamma(x)}\f$
#  @see https://en.wikipedia.org/wiki/Reciprocal_gamma_function
def igamma ( x ) :
    r""" `igamma` function taking into account the uncertainties
    \f[ f(x) = \frac{1}{\Gamma(x)} \f]
    - see https://en.wikipedia.org/wiki/Reciprocal_gamma_function
    """
    fun = getattr ( x , '__igamma__' , None )
    if fun : return fun()
    return _igamma_ ( x )

_psi_ = Ostap.Math.psi 
# =============================================================================
## define polygamma function
def psi ( x , n = 0  ) :
    """ Polygamma function"""
    assert isinstance ( n , integer_types ) and 0 <= n ,\
           'Invalid parameter n=%s' % n  
    fun = getattr ( x , '__psi__' , None )
    if fun : return fun ( n )
    return _psi_ ( x , n )

# =============================================================================
## define digamma function
def digamma  ( x ) :
    """ Digamma function"""
    return psi ( x , 0 )

# =============================================================================
## define trigamma function
def trigamma  ( x ) :
    """ Trigamma function"""
    return psi ( x , 1 )

# =============================================================================
## define polygamma function
def polygamma  ( x , n  ) :
    """Polygamma function"""
    return psi ( x , n )

# =============================================================================
_sinc_ = Ostap.Math.sinc  
# =============================================================================
## define sinc function  \f$ \frac{ \sin x }{x} \f$ 
def sinc  ( x ) :
    """ Sinc function: sin(x)/x
    """
    fun = getattr ( x , '__sinc__' , None )
    if fun : return fun ()
    return _sinc_ ( x )

_beta_ = Ostap.Math.beta  
# =============================================================================
## define Beta function
def beta ( x , y ) :
    """ Beta function"""    
    fun = getattr ( x , '__beta__' , None )
    if fun : return fun ( y )
    fun = getattr ( y , '__beta__' , None )
    if fun : return fun ( x )
    return _beta_ ( x , y )

_lnbeta_ = Ostap.Math.lnbeta  
# =============================================================================
## define log(Beta) function
def lnbeta ( x , y ) :
    """ log(Beta) function"""    
    fun = getattr ( x , '__lnbeta__' , None )
    if fun : return fun ( y )
    fun = getattr ( y , '__lnbeta__' , None )
    if fun : return fun ( x )
    return _lnbeta_ ( x , y )

# =============================================================================
_gd_ = Ostap.Math.gd
# =============================================================================
## define Gudermannian function
#  @see https://en.wikipedia.org/wiki/Gudermannian_function
def gd ( x ) :
    """ Gudermannian function
    - see https://en.wikipedia.org/wiki/Gudermannian_function
    """    
    fun = getattr ( x , '__gd__' , None )
    if fun : return fun ()
    return _gd_ ( x )

# =============================================================================
_gd_inv_ = Ostap.Math.gd_inv
# =============================================================================
## define inverse Gudermannian function
#  @see https://en.wikipedia.org/wiki/Gudermannian_function
def gd_inv ( x ) :
    """ Inverse Gudermannian function
    - see https://en.wikipedia.org/wiki/Gudermannian_function
    """    
    fun = getattr ( x , '__gd_inv__' , None )
    if fun : return fun ()
    return _gd_inv_ ( x )

# =============================================================================
_Ai_ = Ostap.Math.Ai
# =============================================================================
##  Airy function Ai 
#   @see https://en.wikipedia.org/wiki/Airy_function
#   @parameter x argument 
#   @return value of Airy fuinction Ai 
def Ai ( x ) :
    """ Airy function Ai 
    - see https://en.wikipedia.org/wiki/Airy_function
    """
    fun = getattr ( x , '__Ai__' , None )
    if fun : return fun ()
    return _Ai_ ( x )

# =============================================================================
_Bi_ = Ostap.Math.Bi
# =============================================================================
##  Airy function Bi 
#   @see https://en.wikipedia.org/wiki/Airy_function
#   @parameter x argument 
#   @return value of Airy fuinction Bi 
def Bi ( x ) :
    """ Airy function Ai 
    - see https://en.wikipedia.org/wiki/Airy_function
    """
    fun = getattr ( x , '__Bi__' , None )
    if fun : return fun ()
    return _Bi_ ( x )

# =============================================================================
_sech_ = Ostap.Math.sech 
# =============================================================================
## define 'sech' function 
def sech ( x ) :
    """ Sech-function:
    sech(x)=1/cosh(x)
    """
    fun = getattr ( x , '__sech__' , None )
    if fun : return fun()
    return _sech_ ( x )

# =============================================================================
_probit_ = Ostap.Math.probit  
# =============================================================================
## define `probit' function 
#  @see https://en.wikipedia.org/wiki/Probit
def probit ( x ) :
    """ Probit function taking into account the uncertainties
    - see https://en.wikipedia.org/wiki/Probit
    """
    fun = getattr ( x , '__probit__' , None )
    if fun : return fun()
    return _probit_ ( x )

# =============================================================================
## define `min' function \f$ \min (x,y) \f$ 
def minv ( x , y ) :
    """ `minv` function: min (x,y) 
    """
    fun = getattr ( x , '__minv__' , None )
    if fun : return fun ( y )
    fun = getattr ( y , '__minv__' , None )
    if fun : return fun ( x )
    return min ( x , y )

# =============================================================================
## define `max' function \f$ \max (x,y) \f$ 
def maxv ( x , y ) :
    """ `maxv` function: max (x,y) 
    """
    fun = getattr ( x , '__maxv__' , None )
    if fun : return fun ( y )
    fun = getattr ( y , '__maxv__' , None )
    if fun : return fun ( x )
    return max ( x , y )

# =============================================================================
_fma_ = Ostap.Math.fma
# =============================================================================
## evaluate fma(x,y,z) = x*y+x 
#  @param y    (INPUT) the parameter 
#  @param x    (INPUT) the parameter 
#  @param z    (INPUT) the parameter 
#  @param cxy  (INPUT) the correlation coefficient   -1<=c_xy<=1 
#  @param cxz  (INPUT) the correlation coefficient   -1<=c_xz<=1 
#  @param cyz  (INPUT) the correlation coefficient   -1<=c_yz<=1 
#  @return  fma(x,y,z)
#  @code
#  x = ...
#  y = ...
#  z = ...
#  print fma ( x , y , z ) 
#  @endcode 
#  @warning invalid and small covariances are ignored
def fma ( x , y , z , cxy = 0 , cxz = 0 , cyz = 0 ) : 
    """ Evaluate `fma(x,y,z)=x*y+z` with uncertainties
    
    >>> x = ...
    >>> y = ...
    >>> z = ...
    >>> print fma ( x , y , z )
    
    """
    _x = VE ( x )
    _y = VE ( y )
    _z = VE ( z )
    return _fma_ ( _x , _y , _z , cxy , cxz , cyz )

# =============================================================================
_hypot_ = Ostap.Math.hypot
# =============================================================================
## evaluate hypot(x,y) = sqrt(x*x+y*y)
#   \f$ \sqrt( x^2 + y^2 ) \f$
#  @param x (INPUT) the first parameter
#  @param y (INPUT) the second parameter
#  @param c (INPUT) the correlation coefficient  (-1<=c<=1)
#  @return the value of <code>hypot</code> function
#  @warning invalid and small covariances are ignored
def hypot ( x , y , c = 0 ) : 
    """ Evaluate hypot(x,y)=sqrt{x*x+y*y} with uncertainties 
    """
    fun = getattr ( x , '__hypot__' , None )
    if fun : return fun ( y )

    _x = VE ( x )
    _y = VE ( y )
    return _hypot_ ( _x , _y , c )


# =============================================================================
## evaluate sigmoid(x,y) = (1+tanh(xy))/2 
#   \f$ \frac{1+\tanh xy}{2} \f$
#  @param x (INPUT) the first parameter
#  @param y (INPUT) the second parameter
#  @return the value of <code>hypot</code> function
#  @warning invalid and small covariances are ignored
def sigmoid ( x , y = 1  ) : 
    """ Evaluate sigmoid(x,y)=(1+tanh(x*y))/2  with uncertainties 
    """
    fun = getattr ( x , '__sigmoid__' , None )
    if fun : return fun ( y )
    xy = x * y 
    return ( 1.0 + tanh ( xy ) ) * 0.5 

# =============================================================================
_pochhammer = Ostap.Math.pochhammer
# =============================================================================
## \overload calculate Pochhammer's  symbol
#  \f[ (x)^n = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f]
#  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
#  @param x (INPUT) the parameter 
#  @param n (INPUT) the parameter 
#  @return  pochhammer  symbol 
#  @warning invalid and small covariances are ignored 
#  @see Ostap::Math::rising_factorial
#  @see Ostap::Math::falling_factorial
#  @see Ostap::Math::pochhammer 
def pochhammer ( x , n ) :
    r""" calculate Pochhammer's  symbol
    \f[ (x)^n = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
    - see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
    """
    assert is_integer ( n ) and  0<= n < 2**16, \
           "pochhammer: invalid n=%s" % n
    
    fun = getattr ( x , '__pochhammer__' , None )
    if fun : return fun()

    _x = x if isinstance ( x , num_types ) else VE ( x )
    
    return _pochhammer ( _x , n ) 

# =============================================================================
_agm_ = Ostap.Math.agm
_ghm_ = Ostap.Math.ghm
# =============================================================================
## get Arithmetic-geometric mean 
#  @see https://en.wikipedia.org/wiki/Arithmetic%E2%80%93geometric_mean
#  @param x  x-value, \f$ x \ge 0 \f$   
#  @param y  y-value, \f$ y \ge 0 \f$    
#  @return Arithmetic-geometric mean 
#  @see Ostap::Math::agm
def agm ( x  , y ) :
    """ Get Arithmetic-geometric mean 
    - see https://en.wikipedia.org/wiki/Arithmetic%E2%80%93geometric_mean
    - see `Ostap.Math.agm`
    """
    return _agm_ ( x , y )
# =============================================================================
## get Geometric-harmonic  mean         
#  @see https://en.wikipedia.org/wiki/Geometric%E2%80%93harmonic_mean
#  @param x  x-value, \f$ x > 0 \f$   
#  @param y  y-value, \f$ y > 0 \f$    
#  @return Geometric-harmonic  mean         
#  @see Ostap::Math::agm
def ghm ( x  , y ) :
    """ Get Geometric-harmonic  mean         
    - see https://en.wikipedia.org/wiki/Geometric%E2%80%93harmonic_mean
    - see `Ostap.Math.ghm`
    """
    return _ghm_ ( x , y )
# =============================================================================

# =============================================================================
_bessel_Jn  = Ostap.Math.bessel_Jn 
_bessel_Jnu = Ostap.Math.bessel_Jnu
_bessel_Yn  = Ostap.Math.bessel_Yn 
_bessel_Ynu = Ostap.Math.bessel_Ynu
_bessel_In  = Ostap.Math.bessel_In 
_bessel_Inu = Ostap.Math.bessel_Inu
_bessel_Kn  = Ostap.Math.bessel_Kn 
_bessel_Knu = Ostap.Math.bessel_Knu

# =============================================================================
## Evaluate regular Bessel function \f$ J_{\nu}(x)\f$
#  @see https://en.wikipedia.org/wiki/Bessel_function
#  @see Ostap::Math::bessel_Jn
#  @see Ostap::Math::bessel_Jnu
def bessel_J ( nu , x ) : 
    """ Evaluate regular Bessel function J_nu (x)
    - see https://en.wikipedia.org/wiki/Bessel_function
    - see `Ostap.Math.bessel_Jn`
    - see `Ostap.Math.bessel_Jnu`    
    """
    if   isinstance ( nu , integer_types ) : return _bessel_Jn  ( int ( nu ) , x )
    elif isinstance ( nu , float         ) : return _bessel_Jnu (       nu   , x )
    else :
        raise TypeError ("invalid nu/index/order type!") 
    
# =============================================================================
## Evaluate irregular Bessel function \f$ Y_{\nu}(x)\f$
#  @see https://en.wikipedia.org/wiki/Bessel_function
#  @see Ostap::Math::bessel_Yn
#  @see Ostap::Math::bessel_Ynu
def bessel_Y ( nu , x ) : 
    """ Evaluate irreguar Bessel function Y_nu (x)
    - see https://en.wikipedia.org/wiki/Bessel_function
    - see `Ostap.Math.bessel_Yn`
    - see `Ostap.Math.bessel_Ynu`    
    """
    assert isinstance ( nu , num_types ) , "Invalid index/orderNu type!"

    if   isinstance ( nu , integer_types ) : return _bessel_Yn  ( int ( nu ) , x )
    elif isinstance ( nu , float         ) : return _bessel_Ynu (       nu   , x )
    else :
        raise TypeError ("invalid nu/index/order type!") 
    
# =============================================================================
## Evaluate modified Bessel function \f$ I_{\nu}(x)\f$
#  @see https://en.wikipedia.org/wiki/Bessel_function
#  @see Ostap::Math::bessel_In
#  @see Ostap::Math::bessel_Inu
def bessel_I ( nu , x ) : 
    """ Evaluate modified Bessel function I_nu (x)
    - see https://en.wikipedia.org/wiki/Bessel_function
    - see `Ostap.Math.bessel_In`
    - see `Ostap.Math.bessel_Inu`    
    """
    if   isinstance ( nu , integer_types ) : return _bessel_In  ( int ( nu ) , x )
    elif isinstance ( nu , float         ) : return _bessel_Inu (       nu   , x )
    else :
        raise TypeError ("invalid nu/index/order type!") 
    
# =============================================================================
## Evaluate modified Bessel function \f$ K_{\nu}(x)\f$
#  @see https://en.wikipedia.org/wiki/Bessel_function
#  @see Ostap::Math::bessel_Kn
#  @see Ostap::Math::bessel_Knu
def bessel_K ( nu , x ) : 
    """ Evaluate modified Bessel function K_nu (x)
    - see https://en.wikipedia.org/wiki/Bessel_function
    - see `Ostap.Math.bessel_Kn`
    - see `Ostap.Math.bessel_Knu`    
    """
    if   isinstance ( nu , integer_types ) : return _bessel_Kn  ( int ( nu ) , x )
    elif isinstance ( nu , float         ) : return _bessel_Knu (       nu   , x )
    else :
        raise TypeError ("invalid nu/index/order type!") 
    
# =============================================================================
_gauss_pdf_ = Ostap.Math.gauss_pdf
# =============================================================================
## calculate the standard gaussian PDF
#  @param x x-value
#  @param mu mu-parameter (location)
#  @param sigma sigma-parameter (width)
#  @return gaussian PDF 
def gauss_pdf ( x , mu = 0.0 , sigma = 1.0 ) :
    """ Standard gaussian PDF:
    
    >>> x,mu, sigma = ....
    >>> pdf = gauss_pdf ( x  , mu , sigma )
    """
    y =  VE ( x ) 
    return _gauss_pdf_ ( y if 0 < y.cov2() else y.value () , mu , sigma )

_gauss_cdf_ = Ostap.Math.gauss_cdf
# =============================================================================
## calculate the standard gaussian CDF
#  @param x x-value
#  @param mu mu-parameter (location)
#  @param sigma sigma-parameter (width)
#  @return gaussian CDF 
def gauss_cdf ( x , mu = 0.0 , sigma = 1.0 ) :
    """ Standard gaussian CDF:    
    >>> x,mu, sigma = ....
    >>> cdf = gauss_cdf ( x  , mu , sigma )
    """
    y =  VE ( x ) 
    return _gauss_cdf_ ( y if 0 < y.cov2() else y.value () , mu , sigma )

# ==============================================================================
_FD_ = Ostap.Math.fermi_dirac 
# =============================================================================
## Complete Fermi-Dirac integral
#  \f$ F_j(x) = \frac{1}{\Gamma(j+1)}\int^{+\infty}_{0} \frac{t^j}{ \exp ( t-x) + 1 } dt \f$ 
#  @param j parameter
#  @param x argument 
#  @return return value of complete Fermi-Dirac integral 
def fermi_dirac ( j , x ) :
    """ Complete Fermi-Dirac integral
    """
    assert isinstance ( j , integer_types ) and 0 <= j , 'Invalid value of j!'
    ## 
    return _FD_ ( j , x )

# =============================================================================
_li_ = Ostap.Math.li
# =============================================================================
## get the Logarithmic integral function 
#   \f$ li(x) \equiv \int\limits_{0}^{x} \frac{dt}{\log t} \f$ 
#   for \f$ 0 < x \f$ 
#  @see https://en.wikipedia.org/wiki/Logarithmic_integral_function
def li ( x ) :
    """ Get the Logarithmic integral function li(x)
    - see https://en.wikipedia.org/wiki/Logarithmic_integral_function
    """
    fun = getattr ( x , '__li__' , None )
    if fun : return fun ()
    return _li_ ( x ) 
# ============================================================================
_Li_ = Ostap.Math.Li
# =============================================================================
## get the Logarithmic integral function 
#   \f$ li(x) \equiv \int\limits_{0}^{x} \frac{dt}{\log t} \f$ 
#   for \f$ 0 < x \f$ 
#  @see https://en.wikipedia.org/wiki/Logarithmic_integral_function
def Li ( x ) :
    """ Get the Logarithmic integral function Li(x) = li(x) - li(2)
    - see https://en.wikipedia.org/wiki/Logarithmic_integral_function
    """
    fun = getattr ( x , '__Li__' , None )
    if fun : return fun ()
    return _Li_ ( x ) 
# ============================================================================
_Ei_ = Ostap.Math.Ei
# =============================================================================
## get the Exponential  integral function Ei(x) 
#    \f$ Ei(x) \equiv \int\limits_{-\infty}^{x} \frac{e^t}{t}dt \f$ 
# @see https://en.wikipedia.org/wiki/Exponential_integral
def Ei ( x ) :
    """ Get the Exponential  integral function Ei(x) 
     - see https://en.wikipedia.org/wiki/Exponential_integral
    """
    fun = getattr ( x , '__Ei__' , None )
    if fun : return fun ()
    return _Ei_ ( x ) 
# ============================================================================
_Ein_ = Ostap.Math.Ein
# =============================================================================
## get the Exponential  integral function Ein(x) 
#   \f$ Ein(x) \equiv \int\limits_{0}^{x} \frac{1-e^{-t}}{t}dt \f$ 
# @see https://en.wikipedia.org/wiki/Exponential_integral
def Ein ( x ) :
    """ Get the Exponential  integral function Ein(x) 
     - see https://en.wikipedia.org/wiki/Exponential_integral
    """
    fun = getattr ( x , '__Ein__' , None )
    if fun : return fun ()
    return _Ein_ ( x ) 
# ============================================================================
_E1_ = Ostap.Math.E1
# =============================================================================
## get the Exponential  integral function E1(x) 
#  \f$ E_1(x) = -\gamma - \log x + Ein(x) \f$ 
# @see https://en.wikipedia.org/wiki/Exponential_integral
def E1 ( x ) :
    """ Get the Exponential  integral function Ein(x) 
     - see https://en.wikipedia.org/wiki/Exponential_integral
    """
    fun = getattr ( x , '__E1__' , None )
    if fun : return fun ()
    return _E1_ ( x ) 
# ============================================================================
_Si_ = Ostap.Math.Si
# =============================================================================
## Get the integral sine 
#  \f$ Si(x) = \int\limits_{0}^{x} \frac{\sin t}{t}dt \f$ 
#  @see https://en.wikipedia.org/wiki/Trigonometric_integral
def Si ( x ) :
    """ Get the integral sine Si(x) 
     @see https://en.wikipedia.org/wiki/Trigonometric_integral
    """
    fun = getattr ( x , '__Si__' , None )
    if fun : return fun ()
    return _Si_ ( x ) 
# ============================================================================
_Ci_ = Ostap.Math.Ci
# =============================================================================
## Get the integral cosine 
#   \f$ Ci(x) = - \int\limits_{x}^{+\infty} \frac{\cos t}{t}dt \f$ 
#  @see https://en.wikipedia.org/wiki/Trigonometric_integral
def Ci ( x ) :
    """ Get the integral cosine Ci(x) 
     @see https://en.wikipedia.org/wiki/Trigonometric_integral
    """
    fun = getattr ( x , '__Ci__' , None )
    if fun : return fun ()
    return _Ci_ ( x ) 
# ============================================================================
_Cin_ = Ostap.Math.Cin
# =============================================================================
## Get the integral cosine 
#  \f$ Cin(x) = \int\limits_{0}^{x} \frac{1 - \cos t}{t}dt \f$ 
#  @see https://en.wikipedia.org/wiki/Trigonometric_integral
def Cin ( x ) :
    """ Get the integral cosine Cin(x) (entire function)
     @see https://en.wikipedia.org/wiki/Trigonometric_integral
    """
    fun = getattr ( x , '__Cin__' , None )
    if fun : return fun ()
    return _Cin_ ( x ) 

# =============================================================================
_BR_ = Ostap.Math.bring 
# =============================================================================
## Get the Bring radical/ultraradical, a real solution of 
#  the equation \f$ x*5 + x + a = 0 \f$
#  @see https://en.wikipedia.org/wiki/Bring_radical
#  @code
#  a = ...
#  r = bring ( a )
#  r = BR    ( a ) ## ditto
#  @endcode
def bring ( x ) :
    """ Get the Bring radical/ultraradical, a real solution of 
    the equation x*5 + x + a = 0 
     - see https://en.wikipedia.org/wiki/Bring_radical    
    >>> a = ...
    >>> r = bring ( a )
    >>> r = BR    ( a ) ## ditto
    """
    fun = getattr ( x , '__BR__' , None )
    if fun : return fun ()
    return _BR_ ( x )
# =============================================================================
BR = bring



# =============================================================================
_am_ = Ostap.Math.am
# =============================================================================
## Elliptic amplitude \f$ \mathrm{am}(u,m)=\phi\f$, where 
#  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
# @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
def am ( u , m  ) :
    """ Elliptic functon am(u,m)
    - see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions    
    """
    fun = getattr ( u , '__am__' , None )
    if fun : return fun ( m )
    return _am_ ( u , m  )

# =============================================================================
_dn_ = Ostap.Math.dn
# =============================================================================
## Elliptic delta amplitude \f$ \mathrm{sn} (u,m)=\frac{d}{du} \mathrm{am} ( u, m ) \f$, where 
#   \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
#  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
def dn ( u , m  ) :
    """ Elliptic functon dn(u,m)
    - see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions    
    """
    fun = getattr ( u , '__dn__' , None )
    if fun : return fun ( m )
    return _dn_ ( u , m  )

# =============================================================================
_sn_ = Ostap.Math.sn
# =============================================================================
## Elliptic sine amplitude \f$ \mathrm{sn} (u,m)=\sin \mathrm{am} ( u, m ) \f$, where 
#  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
#  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
def sn ( u , m  ) :
    """ Elliptic function sn(u,m)
    - see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions    
    """
    fun = getattr ( u , '__sn__' , None )
    if fun : return fun ( m )
    return _sn_ ( u , m  )

# =============================================================================
_cn_ = Ostap.Math.cn
# =============================================================================
## Elliptic cosine amplitude \f$ \mathrm{sn} (u,m)=\cos \mathrm{am} ( u, m ) \f$, where 
#  \f[ u = \int\limit_0^{\phi} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}\f]
#  @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
def cn ( u , m  ) :
    """ Elliptic function cn(u,m)
    - see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions    
    """
    fun = getattr ( u , '__cn__' , None )
    if fun : return fun ( m )
    return _cn_ ( u , m  )


# =============================================================================
_sc_ = Ostap.Math.sc
# =============================================================================
## Elliptic function sc = sn/cn 
#  \f[ \matmrm{sc}\,(u,m) = \frac{ \mathrm{sn} ( u, m) } { \mathrm{cn} ( u , m ) } \f] 
#   @see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions
def sc ( u , m  ) :
    """ Elliptic function sc(u,m) = sn / cn 
    - see https://en.wikipedia.org/wiki/Jacobi_elliptic_functions    
    """
    fun = getattr ( u , '__sc__' , None )
    if fun : return fun ( m )
    return _sc_ ( u , m  )

# =============================================================================
_elliptic_K_ = Ostap.Math.elliptic_K 
# =============================================================================
## Complete elliptic integral K(k) 
#  @see https://en.wikipedia.org/wiki/Elliptic_integral
def elliptic_K ( k ) :
    """ Complete elliptic integral K(k) 
    - see https://en.wikipedia.org/wiki/Elliptic_integral
    """
    fun = getattr ( k , '__K__' , None )
    if fun : return fun (   )
    return _elliptic_K_ ( k )

# =============================================================================
_elliptic_E_ = Ostap.Math.elliptic_E 
# =============================================================================
## Complete elliptic integral E(k) 
#  @see https://en.wikipedia.org/wiki/Elliptic_integral
def elliptic_E ( k ) :
    """ Complete elliptic integral E(k) 
    - see https://en.wikipedia.org/wiki/Elliptic_integral
    """
    fun = getattr ( k , '__E__' , None )
    if fun : return fun (   )
    return _elliptic_E_ ( k )

# =============================================================================
## FIX
#  @see https://sft.its.cern.ch/jira/browse/ROOT-6627'
_a  = VE ( 1  , 1 )
_b  = VE ( _a     )
if   isequal ( _a.error () , _b.error () ) : pass 
else :
    jira = 'https://sft.its.cern.ch/jira/browse/ROOT-6627'
    from ostap.core.meta_info import root_version 
    logger.warning ( 'The problem %s is not solved yet ( ROOT %s) ' %  ( jira , root_version ) )
    logger.warning ( 'Temporarily disable cast of VE to float' )
    del VE.__float__

# =============================================================================
## Convert p-value/probability into significance (# of sigmas)
#  @code
#  pvalue = ...
#  ns     = significance ( pvalue )
#  ns     = nsigma       ( pvalue ) ## ditto 
#  ns     = nsigmas      ( pvalue ) ## ditto 
#  @endcode 
#  @param p-value
#  @return significance as number of standard deviations 
def significance ( pvalue ) :
    """ Convert p-value/probability into significance (# of sigmas)
    """
    return probit (  ( 1.0 + ( 1.0 - pvalue ) ) / 2.0 ) 

# ============================================================================
nsigma  = significance
# ============================================================================
nsigmas = significance

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme  import docme
    docme ( __name__ , logger = logger )
    
    funcs = [ exp    , expm1   ,
              log    , log10   , log1p    ,
              sqrt   , cbrt    ,
              sin    , cos     , tan      ,
              sinh   , cosh    , tanh     , sech   ,
              asin   , acos    , atan     ,
              asinh  , acosh   , atanh    ,
              erf    , erfc    , erfi     , erfcx  ,
              sinc   ,
              probit ,
              gd     , gd_inv  , 
              Ai     , Bi      , 
              gamma  , tgamma  , lgamma   , igamma ,
              psi    , digamma , trigamma ,
              gauss_pdf ,
              gauss_cdf ] 
    
    from ostap.math.derivative import EvalVE
    funcs += [ EvalVE ( math.sin , math.cos ) ,
               EvalVE ( math.sin )            ]
    
    vars  = [ VE ( 0.001 , 0.0001**2 ) , VE(1,0) , VE(1,0.1**2) , VE(10,0.01**2) ]
    
    for v in vars :
        logger.info ( 'Var = %s ' % v )
        for f in funcs :
            logger.info ( "\t%12s\t%s = %s " % ( f.__name__ , v ,  f(v) ) )
            
    logger.info ( 80*'*')
    
# =============================================================================
##                                                                      The END
# =============================================================================
