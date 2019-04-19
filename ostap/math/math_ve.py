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
    'exp'        , 'expm1'      ,
    'log'        , 'log10'      , 'log1p'  , 
    'sqrt'       , 'cbrt'       , 'pow'    ,   
    'sin'        , 'cos'        ,  'tan'   , 
    'sinh'       , 'cosh'       , 'tanh'   , 'sech'   ,
    'asin'       , 'acos'       , 'atan'   , 
    'asinh'      , 'acosh'      , 'atanh'  ,
    'erf'        , 'erfc'       , 'erfi'   , 'erfcx'  ,
    'probit'     , 'pochhammer' , 
    'gamma'      , 'tgamma'     , 'lgamma' , 'igamma' , 
    'exp2'       , 'log2'       ,
    'gauss_pdf'  , 'gauss_cdf'  ,
    'hypot'      , 'fma'        ,
    )
# =============================================================================
import ROOT,math
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.math_ve' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
from   ostap.math.ve    import VE
from   ostap.math.base  import cpp, iszero, isequal
from   ostap.core.ostap_types import num_types, is_integer 

_ln2_i = 1/math.log(2.0)                 ## useful constant 
# =============================================================================

# =============================================================================
## define ``exp'' function 
def exp ( x ) :
    """ 'exp' function taking into account the uncertainties
    """
    fun = getattr ( x , '__exp__' , None )
    if fun : return fun()
    return math.exp ( x )

# =============================================================================
## define ``exp2'' function 
def exp2 ( x ) :
    """ 'exp2' function taking into account the uncertainties
    """
    fun = getattr ( x , '__exp2__' , None )
    if fun : return fun()
    return 2**x 

# =============================================================================
## define ``expm1'' function 
def expm1 ( x ) :
    """ 'expm1' function taking into account the uncertainties
    """
    fun = getattr ( x , '__expm1__' , None )
    if fun : return fun()
    return math.expm1 ( x )

# =============================================================================
## define ``log'' function 
def log ( x ) :
    """'log' function taking into account the uncertainties
    """
    fun = getattr ( x , '__log__' , None )
    if fun : return fun()
    return math.log ( x )

# =============================================================================
## define ``log2'' function 
def log2 ( x ) :
    """'log2' function taking into account the uncertainties
    """
    fun = getattr ( x , '__log2__' , None )
    if fun : return fun()
    ## 
    return _ln2_i * math.log ( x ) 

# =============================================================================
## define ``log10'' function 
def log10 ( x ) :
    """'log10' function taking into account the uncertainties
    """
    fun = getattr ( x , '__log10__' , None )
    if fun : return fun()
    return math.log10 ( x )

# =============================================================================
## define ``log1p'' function 
def log1p ( x ) :
    """'log1p' function taking into account the uncertainties
    """
    fun = getattr ( x , '__log1p__' , None )
    if fun : return fun()
    return math.log1p ( x )

# =============================================================================
## define ``sqrt'' function 
def sqrt ( x ) :
    """'sqrt' function taking into account the uncertainties
    """
    fun = getattr ( x , '__sqrt__' , None )
    if fun : return fun()
    return math.sqrt ( x )

# =============================================================================
## define ``cbrt'' function 
def cbrt ( x ) :
    """'cbrt' function taking into account the uncertainties
    """
    fun = getattr ( x , '__cbrt__' , None )
    if fun : return fun()
    return math.pow ( x , 1.0/3.0 )

# =============================================================================
## define ``pow'' function 
def pow ( x , y , *a ) :
    """'pow' function taking into account the uncertainties
    """
    if   isinstance ( x , VE ) or isinstance ( y , VE ) : return x**y 
    return math.pow ( x , y , *a ) 
 
# =============================================================================
## define ``sin'' function 
def sin ( x ) :
    """'Sine' function taking into account the uncertainties
    """
    fun = getattr ( x , '__sin__' , None )
    if fun : return fun()
    return math.sin ( x )

# =============================================================================
## define ``cos'' function 
def cos ( x ) :
    """'Cosine' function taking into account the uncertainties
    """
    fun = getattr ( x , '__cos__' , None )
    if fun : return fun()
    return math.cos ( x )

# =============================================================================
## define ``tan'' function 
def tan ( x ) :
    """'tangent' function taking into account the uncertainties
    """
    fun = getattr ( x , '__tan__' , None )
    if fun : return fun()
    return math.tan ( x )

# =============================================================================
## define ``sinh'' function 
def sinh ( x ) :
    """'Sinh' function taking into account the uncertainties
    """
    fun = getattr ( x , '__sinh__' , None )
    if fun : return fun()
    return math.sinh ( x )

# =============================================================================
## define ``cosh'' function 
def cosh ( x ) :
    """'Cosh' function taking into account the uncertainties
    """
    fun = getattr ( x , '__cosh__' , None )
    if fun : return fun()
    return math.cosh ( x )

# =============================================================================
## define ``tanh'' function 
def tanh ( x ) :
    """'tanh' function taking into account the uncertainties
    """
    fun = getattr ( x , '__tanh__' , None )
    if fun : return fun()
    return math.tanh ( x )

_erf_   = cpp.Ostap.Math.erf 
_erfc_  = cpp.Ostap.Math.erfc
_erfi_  = cpp.Ostap.Math.erfi
_erfcx_ = cpp.Ostap.Math.erfcx

# =============================================================================
## define ``erf'' function 
#  @see https://en.wikipedia.org/wiki/Error_function
def erf ( x ) :
    """ Error function taking into account the uncertainties
    - see https://en.wikipedia.org/wiki/Error_function
    """
    fun = getattr ( x , '__erf__' , None )
    if fun : return fun()
    return _erf_ ( x )

# =============================================================================
## define ``erfc'' function 
#  @see https://en.wikipedia.org/wiki/Error_function
def erfc ( x ) :
    """ Complemenatry error function taking into account the uncertainties
    - see https://en.wikipedia.org/wiki/Error_function
    """
    fun = getattr ( x , '__erfc__' , None )
    if fun : return fun()
    return _erfc_( x )

# =============================================================================
## define ``erfcx'' function 
#  @see https://en.wikipedia.org/wiki/Error_function
def erfcx ( x ) :
    """ Complemenatry scaled error function taking into account the uncertainties
    - see https://en.wikipedia.org/wiki/Error_function
    """
    fun = getattr ( x , '__erfcx__' , None )
    if fun : return fun()
    return _erfcx_( x )

# =============================================================================
## define ``erfi'' function 
#  @see https://en.wikipedia.org/wiki/Error_function
def erfi ( x ) :
    """ Imaginary error function taking into account the uncertainties
    - see https://en.wikipedia.org/wiki/Error_function
    """
    fun = getattr ( x , '__erfi__' , None )
    if fun : return fun()
    return _erfi_( x )


# =============================================================================
## define ``asin'' function 
def asin ( x ) :
    """'asin' function taking into account the uncertainties
    """
    fun = getattr ( x , '__asin__' , None )
    if fun : return fun()
    return math.asin ( x )

# =============================================================================
## define ``acos'' function 
def acos ( x ) :
    """'acos' function taking into account the uncertainties
    """
    fun = getattr ( x , '__acos__' , None )
    if fun : return fun()
    return math.acos ( x )

# =============================================================================
## define ``atan'' function 
def atan ( x ) :
    """'atan' function taking into account the uncertainties
    """
    fun = getattr ( x , '__atan__' , None )
    if fun : return fun()
    return math.atan ( x )

# =============================================================================
## define ``asinh'' function 
def asinh ( x ) :
    """'asinh' function taking into account the uncertainties
    """
    fun = getattr ( x , '__asinh__' , None )
    if fun : return fun()
    return math.asinh ( x )

# =============================================================================
## define ``acosh'' function 
def acosh ( x ) :
    """'acosh' function taking into account the uncertainties
    """
    fun = getattr ( x , '__acosh__' , None )
    if fun : return fun()
    return math.acosh ( x )

# =============================================================================
## define ``atanh'' function 
def atanh ( x ) :
    """'atanh' function taking into account the uncertainties
    """
    fun = getattr ( x , '__atanh__' , None )
    if fun : return fun()
    return math.atanh ( x )


# =============================================================================
## define ``tgamma'' function 
def tgamma ( x ) :
    """'tgamma' function taking into account the uncertainties
    """
    fun = getattr ( x , '__tgamma__' , None )
    if fun : return fun()
    return math.gamma ( x )

## define ``gamma'' function 
gamma = tgamma

# =============================================================================
## define ``lgamma'' function 
def lgamma ( x ) :
    """'lgamma' function taking into account the uncertainties
    """
    fun = getattr ( x , '__lgamma__' , None )
    if fun : return fun()
    return math.lgamma ( x )

_igamma_ = cpp.Ostap.Math.igamma 
# =============================================================================
## define ``igamma'' function
#  \f$ f(x) = \frac{1}{\Gamma(x)}\f$
#  @see https://en.wikipedia.org/wiki/Reciprocal_gamma_function
def igamma ( x ) :
    r"""'igamma' function taking into account the uncertainties
    \f$ f(x) = \frac{1}{\Gamma(x)}\f$
    - see https://en.wikipedia.org/wiki/Reciprocal_gamma_function
    """
    fun = getattr ( x , '__igamma__' , None )
    if fun : return fun()
    return _igamma_ ( x )

_sech_ = cpp.Ostap.Math.sech 
# =============================================================================
## define ``sech'' function 
def sech ( x ) :
    """ Sech-function:
    sech(x)=1/cosh(x)
    """
    fun = getattr ( x , '__sech__' , None )
    if fun : return fun()
    return _sech_ ( x )

_probit_ = cpp.Ostap.Math.probit  
# =============================================================================
## define ``probit'' function 
#  @see https://en.wikipedia.org/wiki/Probit
def probit ( x ) :
    """ Probit function taking into account the uncertainties
    - see https://en.wikipedia.org/wiki/Probit
    """
    fun = getattr ( x , '__probit__' , None )
    if fun : return fun()
    return _probit_ ( x )

_fma_ = cpp.Ostap.Math.fma
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
    """ Evaluate fma(x,y,z)=x*y+z with uncertainties
    >>> x = ...
    >>> y = ...
    >>> z = ...
    >>> print fma ( x , y , z ) 
    """
    _x = VE ( x )
    _y = VE ( y )
    _z = VE ( z )
    return _fma_ ( _x , _y , _z , cxy , cxz , cyz )

_hypot_ = cpp.Ostap.Math.hypot
# =============================================================================
## evaluate hypot(x,y) = sqrt(x*x+y*y)
#   \f$ \sqrt( x^2 + y^2 ) \f$
#  @param x (INPUT) the first parameter
#  @param y (INPUT) the second parameter
#  @param c (INPUT) the correlation coefficient  (-1<=c<=1)
#  @return the value of <code>hypot</code> function
#  @warning invalid and small covariances are ignored
def hypot ( x , y , c = 0 ) : 
    """ evaluate hypot(x,y)=sqrt{x*x+y*y} witn uncertainties 
    """
    _x = VE ( x )
    _y = VE ( y )
    return _hypot_ ( _x , _y , c )


_pochhammer = cpp.Ostap.Math.pochhammer
# =============================================================================
## calculate Pochhammer's  symbol
#   \f[ (x)^n = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
#  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
#  @param x (INPUT) the parameter 
#  @param n (INPUT) the parameter 
#  @return  pochhammer  symbol 
#  @warning invalid and small covariances are ignored 
#  @see Ostap::Math::rising_factorial
#  @see Ostap::Math::falling_factorial
#  @see Ostap::Math::pochhammer 
def pochhammer ( x , n ) :
    """ calculate Pochhammer's  symbol
    \f[ (x)^n = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
    @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
    """
    assert is_integer ( n ) and  0<= n < 2**16, \
           "pochhammer: invalid n=%s" % n
    
    fun = getattr ( x , '__pochhammer__' , None )
    if fun : return fun()

    _x = x if isinstance ( x , num_types ) else VE ( x )
    
    return _pochhammer ( _x , n ) 

_gauss_pdf_ = cpp.Ostap.Math.gauss_pdf
# =============================================================================
## calculate the standard gaussian PDF
#  @param x x-value
#  @param mu mu-parameter (location)
#  @param sigma sigma-parameter (width)
#  @return gaussian PDF 
def gauss_pdf( x , mu = 0.0 , sigma = 1.0 ) :
    """Standard gaussian PDF
    >>> x,mu, sigma = ....
    >>> pdf = gauss_pdf ( x  , mu , sigma )
    """
    m =       float ( mu    )
    s = abs ( float ( sigma ) )  
    y      = ( VE ( x ) - m ) / s
    ## 
    return _gauss_pdf_ ( y ) / s if  0 < y.cov2() else _gauss_pdf_ ( y.value() ) / s

_gauss_cdf_ = cpp.Ostap.Math.gauss_cdf
# =============================================================================
## calculate the standard gaussian CDF
#  @param x x-value
#  @param mu mu-parameter (location)
#  @param sigma sigma-parameter (width)
#  @return gaussian CDF 
def gauss_cdf ( x , mu = 0.0 , sigma = 1.0 ) :
    """Standard gaussian CDF
    >>> x,mu, sigma = ....
    >>> cdf = gauss_cdf ( x  , mu , sigma )
    """
    m  = float ( mu )
    s  = abs ( float ( sigma ) )
    y  = ( VE ( x ) - m ) / s
    ##
    return _gauss_cdf_ ( y ) if 0 < y.cov2() else _gauss_cdf_ ( y.value() )

        
# =============================================================================
## FIX
#  @see https://sft.its.cern.ch/jira/browse/ROOT-6627'
_a  = VE( 1  , 1 )
_b  = VE( _a     )
if     isequal ( _a.error () , _b.error () ) : pass 
else :
    jira = 'https://sft.its.cern.ch/jira/browse/ROOT-6627'
    vers = ROOT.gROOT.GetVersion() 
    logger.warning ( 'The problem %s is not solved yet ( ROOT %s) ' %  ( jira , vers ) )
    logger.warning ( 'Temporarily disable cast of VE to float' )
    del VE.__float__


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme  import docme
    docme ( __name__ , logger = logger )
    
    funcs = [ exp    , expm1  ,
              log    , log10  , log1p  ,
              sqrt   , cbrt   ,
              sin    , cos    , tan    ,
              sinh   , cosh   , tanh   , sech   ,
              asin   , acos   , atan   ,
              asinh  , acosh  , atanh  ,
              erf    , erfc   , erfi   , erfcx  ,
              probit ,
              gamma  , tgamma , lgamma , igamma ,
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
# The END
# =============================================================================
