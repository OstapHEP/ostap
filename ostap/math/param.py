#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/param.py
#  Module with utilities for parameterization of functions/histograms 
#
## (1) using histogram data only:
# 
# - as      Legendre         sum  
# - as      Chebyshev        sum  
# - as      Fourier          sum  (relies on numpy.fft.rfft & numpy.linspace)
# - as      Cosine/Fourier   sum  (relies on scipy.fftpack)
# - as      Bernstein/Bezier sum  
# - as even Bernstein/Bezier sum  
# - as rational function inspired by Floater-Hormann's rationa barycentric interpolant

# Typical usage:
#
# @code 
# histo = ...
# fun   = histo.legendre_sum(5)
# print b
# x = ...
# print 'fun(%s)=%s' % ( x , fun(x) ) 
# @endcode
# 
# A little bit more generic example: 
#
# @code
# my_fun = lambda x : x * x 
# fun    = legendre_sum( my_fun , 5, xmin = -1 , xmax = 1 ) 
# print fun
# x = ...
# print 'fun(%s)=%s' % ( x , fun(x) )
# @endcode
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date   2011-06-07
# =============================================================================
"""Module with utilities for parameterization of functions/histograms

## (1) using histogram data only:

- as      Legendre         sum  
- as      Chebyshev        sum  
- as      Fourier          sum  (relies on numpy.fft.rfft & numpy.linspace)
- as      Cosine/Fourier   sum  (relies on scipy.fftpack)
- as      Bernstein/Bezier sum  
- as even Bernstein/Bezier sum  

Typical usage:

>>> histo = ...
>>> fun   = histo.legendre_sum(5)
>>> print b
>>> x = ...
>>> print 'fun(%s)=%s' % ( x , fun(x) ) 

A little bit more generic

>>> my_fun = lambda x : x * x 
>>> fun    = legendre_sum( my_fun , 5, xmin = -1 , xmax = 1 ) 
>>> print fun
>>> x = ...
>>> print 'fun(%s)=%s' % ( x , fun(x) ) 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    ## several stand-alone function for parameteriztaion of of functions 
    'legendre_sum'      , ## Legendre         sum for the given function/object 
    'chebyshev_sum'     , ## Chebyshev        sum for the given function/object
    'bezier_sum'        , ## Bezier/Bernstein sum for the given function/object
    'bernstein_sum'     , ## - ditto -
    'beziereven_sum'    , ## even Bezier/Bernstein sum for the given function/object
    'bernsteineven_sum' , ## - ditto -
    'rational_fun'      , ## Rational (a'la Floater-Hormann's interpolant) parameterisaton 
    ) 
# =============================================================================
from   ostap.core.core         import cpp, Ostap
from   ostap.core.ostap_types  import is_integer, num_types
from   ostap.utils.utils       import crange 
import ostap.math.models
import ROOT, ctypes
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.math.param' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Some parameterization utilities')
# =============================================================================
inf_pos =  float('inf') ## positive infinity
inf_neg = -float('inf') ## negative infinity
# =============================================================================
## helper function to catch xmin/xmax from histogram/function 
def _get_xminmax_ ( func , xmin , xmax , name = 'get_xminmax') :
    """ Helper function to catch xmin/xmax from histogram/function 
    """
    ## xmin 
    if not isinstance ( xmin , num_types ) :
        if   hasattr ( func , 'xmin'     ) : xmin = func.xmin     () 
        elif hasattr ( func , 'GetXmin'  ) : xmin = func.GetXmin  ()
        elif hasattr ( func , 'GetXaxis' ) : xmin = func.GetXaxis ().GetXmin()
        elif hasattr ( func , 'GetRange' ) :
            xmn  = ctypes.c_double()
            xmx  = ctypes.c_double()
            func.GetRange(xmn,xmx)
            xmn  = float ( xmn.value )
            xmx  = float ( xmx.value )
            xmin = xmn 
        else :
            raise AttributeError( "%s: unable to catch xmin %s" % ( name , xmin ) )
        
    ## xmax 
    if not isinstance ( xmax , num_types ) :
        if   hasattr ( func , 'xmax'     ) : xmax = func.xmax     () 
        elif hasattr ( func , 'GetXmax'  ) : xmax = func.GetXmax  () 
        elif hasattr ( func , 'GetXaxis' ) : xmax = func.GetXaxis ().GetXmax()  
        elif hasattr ( func , 'GetRange' ) :
            xmn  = ctypes.c_double()
            xmx  = ctypes.c_double()
            func.GetRange ( xmn , xmx )
            xmn  = float ( xmn.value )
            xmx  = float ( xmx.value )            
            xmax = xmx
        else :
            raise AttributeError( "%s: unable to catch xmax %s" % ( name , xmax ) )

    return min ( xmin , xmax ) , max ( xmin , xmax )

# =============================================================================
## make a function representation in terms of Legendre polynomials
#  @code 
#  func = lambda x : x * x
#  fsum = legendre_sum ( func , 4 , -1 , 1 )
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::LegendreSum
#  @author Vanya Belyaev Ivan.Belyaev@iter.ru
#  It is not very CPU efficient, but stable enough...
#  @date 2015-07-26
def legendre_sum ( func , N , xmin , xmax , **kwargs ) :
    """ Make a function representation in terms of Legendre polynomials
    [It is not very CPU efficient, but stable enough...]
    >>> func = lambda x : x * x
    >>> fsum = legendre_sum ( func , 4 , -1 , 1 )
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    see Ostap.Math.LegendreSum
    """
    from copy import deepcopy

    assert is_integer ( N ) and 0 <= N, "legendre_sum: invalid N %s " % N 

    xmin,xmax = _get_xminmax_ ( func , xmin , xmax , 'legendre_sum' )

    ## prepare the result
    lsum  = Ostap.Math.LegendreSum ( N , xmin , xmax )
    
    ## the type for the basic legendre polynomials, used for integration 
    L_    = Ostap.Math.Legendre 

    ## transform x to local variable -1<t<1 
    tx    = lambda x : lsum.t ( x )

    idx   = 1.0 / ( xmax - xmin ) ## scale factor 

    from ostap.math.integral import integral as _integral 
    
    args = deepcopy ( kwargs ) if kwargs else {} 
    for n in range ( N + 1 ) :
        
        li     = L_ ( n ) 
        fun_n  = lambda x : func ( x ) * li ( tx ( x ) )        
        c_n    = _integral ( fun_n , xmin , xmax , **args ) * ( 2 * n + 1 ) * idx
            
        lsum.setPar ( n , c_n ) 
        
    return lsum

# =============================================================================
## make a function representation in terms of Chebyshev polynomials
#  @code 
#  func = lambda x : x * x
#  fsum = chebyshev_sum ( func , 4 , -1 , 1 )
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Gaudi::Math::ChebyshevSum
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2015-07-26
def chebyshev_sum ( func , N , xmin , xmax ) :
    """ Make a function representation in terms of Chebyshev polynomials
    >>> func = lambda x : x * x
    >>> fsum = chebyshev_sum ( func , 4 , -1 , 1 )
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    see Ostap::Math::ChebyshevSum
    """
    assert is_integer ( N ) and 0 <= N, "chebyshev_sum: invalid N %s " % N 

    xmin,xmax = _get_xminmax_ ( func , xmin , xmax , 'chebyshev_sum' )

    ## prepare the result
    csum = Ostap.Math.ChebyshevSum ( N , xmin , xmax )
    
    ## transform x to local variable -1<t<1 
    tx   = lambda x : csum.t ( x )

    import math 
    _cos = math.cos
    _piN = math.pi/N

    ## conversion from  -1<t<1 to  xmin<x<xmax 
    xt   = lambda t : csum.x ( t )
    
    ## precalculate function
    fk   = [ func ( xt ( _cos ( _piN * ( k + 0.5 ) ) ) ) for k in range ( 0 , N ) ]
    
    scale   = 2.0 / N ## scale factor 
    for n in range ( N + 1 ) :
        
        c_n = 0.0
        for k in range ( 0, N  ) :
            c_n += fk [ k ] * _cos ( _piN * n * ( k + 0.5 ) )
            
        if 0 == n : c_n *= 0.5
        csum.setPar ( n , c_n * scale ) 
        
    return csum

# =============================================================================
try : # =======================================================================
    # =========================================================================
    import numpy
    # =========================================================================
    ## make a function representation in terms of Fourier series
    #  @code 
    #  func = lambda x : x * x
    #  fsum = fourier_sum ( func , 4 , -1 , 1 )
    #  print fsum
    #  x = ...
    #  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    #  @endcode 
    #  @see Ostap::Math::FourierSum
    #  @author Vanya Belyaev Ivan.Belyaev@itep.ru
    #  @date 2015-07-26
    def fourier_sum ( func , N , xmin , xmax , fejer = False ) :
        """ Make a function/histogram representation in terms of Fourier series
        >>> func = lambda x : x * x
        >>> fsum = fourier_sum ( func , 4 , -1 , 1 )
        >>> print fsum
        >>> x = ...
        >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
        """
        assert is_integer ( N ) and 0 <= N , "fourier_sum: invalid N %s " % N 
        
        xmin , xmax = _get_xminmax_ ( func , xmin , xmax , 'fourier_sum' )
        
        ## 0) check the type of return value in the mid point
        xmid = 0.5 * ( xmin + xmax )
        fval = func ( xmid )
        #
        if isinstance ( fval , num_types ) : the_fun = func
        else : the_fun = lambda x : float ( func ( x ) ) 
        
        ## 1) vectorize the function
        vfunc = numpy.vectorize ( the_fun ) 
        
        ## prepare sampling 
        f_sample = 2 * N
        t, dt    = numpy.linspace ( xmin , xmax , f_sample + 2, endpoint=False , retstep=True )
        
        ## make Fast Fourier Transform 
        y  = numpy.fft.rfft ( vfunc ( t ) ) ## / t.size
        y *= 2.0 / t.size
        
        #
        ## decode the results:
        #
        
        a0 = y[0   ].real
        a  = y[1:-1].real
        b  = y[1:-1].imag
        
        #
        ## prepare the output
        #
        fsum = Ostap.Math.FourierSum ( N , xmin , xmax , fejer )
        
        #
        ## fill it!
        #
        fsum.setPar( 0, a0 )
        for i in range ( 1 , N + 1 ) :
            
            if 0 == i % 2 :
                fsum.setA ( i ,  a[i-1] )
                fsum.setB ( i , -b[i-1] )
            else          :
                fsum.setA ( i , -a[i-1] )
                fsum.setB ( i ,  b[i-1] )
                
        return fsum

    __all__ = __all__ + (
        'fourier_sum' , ## Fourier        sum for the given function/object
        )
    
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    numpy = None 
    pass

# =============================================================================
try : # =======================================================================
    # =========================================================================
    import numpy
    import scipy 
    import scipy.fftpack 
    # =========================================================================
    ## make a function representation in terms of cosine Fourier series
    #  @code 
    #  func = lambda x : x * x
    #  fsum = cosine_sum ( func , 4 , -1 , 1 )
    #  print fsum.pars()
    #  x = ...
    #  print 'FUN(%s) = %s ' % ( x , fsum ( x ) )
    #  @endcode 
    #  @see Ostap::Math::CosineSum
    #  @author Vanya Belyaev Ivan.Belyaev@itep.ru
    #  @date 2015-07-26
    def cosine_sum ( func , N , xmin , xmax , fejer = False ) :
        """ Make a function/histiogram representation in terms of Fourier series
        >>> func = lambda x : x * x
        >>> fsum = fourier_sum ( func , 4 , -1 , 1 )
        >>> print fsum
        >>> x = ...
        >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
        """
        assert is_integer ( N ) and 0 <= N , "cosine_sum: invalid N %s " % N 
        
        xmin,xmax = _get_xminmax_ ( func , xmin , xmax , 'cosine_sum' )

        #  0) check the functon type at the misdpoint
        xmid = 0.5 * ( xmin + xmax )
        fval = func ( xmid )
        
        if isinstance ( fval , num_types ) : the_fun = func
        else : the_fun = lambda x : float ( func ( x ) ) 
        
        ## 1) prepare sampling
        t     = numpy.linspace ( xmin , xmax , N + 1 , endpoint=True )
                
        ## 2) vectorize the function
        vfunc = numpy.vectorize ( the_fun )
        
        ## make cosine fourier transform 
        r = scipy.fftpack.dct ( vfunc ( t ) , 1 ) / N 
        
        #
        ## decode the results & prepare the output
        #
        csum = Ostap.Math.CosineSum ( N, xmin , xmax , fejer )
        for i in range ( 0 , N + 1 ) : csum.setPar ( i , r[i] )
        
        return csum
    
    __all__ = __all__ + (
        'cosine_sum' , ## Cosine Fourier sum for the given function/object 
        )
    
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    scipy = None
    numpy = None
    
# =============================================================================
## make a function representation in terms of Bezier sum
#  (sum over Bernstein polynomials)
#  @code 
#  func = lambda x : x * x
#  fsum = bezier_sum ( func , 2 , -1 , 1 )
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::Bernstein
#  @param func (INPUT) the function
#  @param N    (INPUT) the polynomial degree
#  @param xmin (INPUT) the low  edge of the domain
#  @param xmax (INPUT) the high edge of the domain
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2016-07-03
def bezier_sum ( func , N , xmin , xmax , **kwargs ) :
    """ Make a function/histiogram representation in terms of Bezier sum
    (sum of Bernstein Polynomials)
    >>> func = lambda x : x * x
    >>> fsum = bezier_sum ( func , 4 , 0 , 1 )
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    """
    assert is_integer ( N ) and 0 <= N , "bezier_sum: invalid N %s " % N 
    
    xmin,xmax = _get_xminmax_ ( func , xmin , xmax , 'bezier_sum' )

    ## result 
    bsum = Ostap.Math.Bernstein ( N , xmin , xmax )
 
    from ostap.math.integral import integral as _integral 

    args = {}
    
    basis = Ostap.Math.BernsteinDualBasis.basis ( N  )
    assert basis , 'cannot acquire Bernstein dual Basis!'

    ##for i, dual  in  enumerate ( range ( 0 , N + 1 ) :
    for i, dual  in  enumerate ( basis ) : 

        ## if not index in _bernstein_dual_basis_ :
        ##     ## create the dual basic function
        ##      _DUAL = Ostap.Math.BernsteinDualBasis
        ##     _dual = _DUAL ( *index )
        ##     _bernstein_dual_basis_ [ index ] = _dual
        ##  ## get the dual basic function
        ##  dual = _bernstein_dual_basis_[ index ]
        
        ## dual = Ostap.Math.BernsteinDualBasis.element ( N , i )
        ## assert dual , 'cannot acquare BernsteinDual Basis!'
            
        ## get the integration function 
        fun_i = lambda x : float ( func ( x ) ) * dual ( bsum.t( x ) )

        ## use integration 
        if kwargs : args = deepcopy ( kwargs )
        c_i = _integral ( fun_i , xmin , xmax , **args ) / ( xmax - xmin )
            
        bsum.setPar( i , c_i )
        
    return bsum


# ============================================================================-
## make a function representation in terms of Bezier sum (sum over Bernstein polynomials)
bernstein_sum = bezier_sum

# ============================================================================-
## make a function representation in terms of even Bezier sum
#  (sum over Bernstein polynomials), such as
#  \f$ f( \frac{x_{min}+x_{max}}{2} - x ) = f(  \frac{x_{min}+x_{max}}{2} + x ) \f$ 
#  @code 
#  func = lambda x : x * x
#  fsum = beziereven_sum ( func , 2 , -1 , 1 )
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::BernsteinEven
#  @param func (INPUT) the function
#  @param N    (INPUT) the degree of even polynomial 
#  @param xmin (INPUT) the low  edge of the domain
#  @param xmax (INPUT) the high edge of the domain
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2016-07-03
def beziereven_sum ( func , N , xmin , xmax , **kwargs ) :
    """Make a function/histogram representation in terms of Bezier sum
    (sum of Bernstein Polynomials)
    >>> func = lambda x : x * x
    >>> fsum = beziereven_sum ( func , 2 , -1 , 1 )
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    """
    assert is_integer ( N ) and 0 <= N , "beziereven_sum: invalid N %s " % N 
    
    xmin,xmax = _get_xminmax_ ( func , xmin , xmax , 'beziereven_sum' )

    ## result 
    bsum  = Ostap.Math.BernsteinEven ( N , xmin , xmax )
    
    npars = bsum.bernstein().npars() 
    ## constants 
    b_i   = npars * [ 0.0 ]
    
    xmid  = 0.5 * ( xmin + xmax ) 
    ## symmetric function: f(xmid-x)=f(xmid+x) 
    def _sym_func_ ( x ) :
        x1 =            x
        x2 = 2 * xmid - x
        y1 = float ( func ( x1 ) )
        y2 = float ( func ( x2 ) )
        return 0.5 * ( y1 + y2 ) 
    
    from ostap.math.integral import integral as _integral 
    
    args = {}

    ## basis = Ostap.Math.BernsteinDualBasis.element ( bsum.degree() , i ) 
    ## assert dual, 'Cannot acquare element of Bernstein dual basis!'

    for i in  range ( len ( b_i ) )  : 
        
        ## index = 2 * N + 1 , i
        index = bsum.degree () , i
        
        ## if not index in _bernstein_dual_basis_ :
        ##    ## create the dual basic function
        ##    _DUAL = Ostap.Math.BernsteinDualBasis
        ##    _dual = _DUAL ( *index )
        ##    _bernstein_dual_basis_ [ index ] = _dual
        ##  ## get the dual basic function
        ##  dual = _bernstein_dual_basis_[ index ]
        dual = Ostap.Math.BernsteinDualBasis.element ( bsum.degree() , i ) 
        assert dual, 'Cannot acquare element of Bernstein dual basis!'
        
        ## get the integration function 
        fun_i = lambda x : _sym_func_ ( x ) * dual ( bsum.t( x ) )

        ## use integration 
        if kwargs : args = deepcopy ( kwargs )
        c_i = _integral ( fun_i , xmin , xmax , **args ) / ( xmax - xmin )
            
        b_i[i] = c_i 
        
    ## fill result with symmetrized coefficients
    last = npars - 1 
    for i in bsum :
        bsum.setPar ( i , 0.5 * ( b_i [ i ] + b_i [ last - i ] ) ) 

    return bsum



# =============================================================================
## make a function representation in terms of rational function
#  based on the Floater-Hormann's rational barycentric interpolant 
#  @code 
#  func = lambda x : x * x
#  rat  = rational( func , 3 , 3  , xmin = -1 , xmax = 1 )
#  print ( rat ) 
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , rat  ( x ) ) 
#  @endcode 
#  @see Ostap::Math::Rational
#  @see Ostap::Math::FloaterHormann
#  @author Vanya Belyaev Ivan.Belyaev@iter.ru
#  @date 2023-09-18
def rational_fun  ( func , n , d  , xmin , xmax , **kwargs ) :
    """ Make a function representation in terms of rational function
    based on the Floater-Hormann's rational barycentric interpolant
    >>> func = lambda x : x * x
    >>> rat  = rational ( func , N = 4 , d = 1  , xmin = -1 , xmax = 1 )
    >>> print ( rat ) 
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , rat ( x ) ) 
    - see Ostap.Math.Rational 
    - see Ostap.Math.FloaterHormann
    """
    from copy import deepcopy
    
    assert is_integer ( n ) and 0 <= n      , "rational: invalid n %s " % N
    assert is_integer ( d ) and 0 <= d <= n , "rational: invalid d %s " % d
    
    xmin , xmax = _get_xminmax_ ( func , xmin , xmax , 'rational' )
    
    ## prepare the result
    rational = Ostap.Math.Rational( n + 1 , d , xmin , xmax )
    for i, r in enumerate ( crange ( xmin , xmax ,  n + 1 ) ) :
        rational [ i ] = float ( func ( r ) ) 

    return rational 

# =============================================================================
## make a function representation in terms of even Bezier sum (sum over Bernstein polynomials)
bernsteineven_sum = beziereven_sum


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if not 'fourier_sum' in __all__ :
        logger.warning ( "Since numpy is not available, fourier_sum is disabled") 
    if not 'cosine_sum'  in __all__ :
        logger.warning ( "Since scipy is not available, cosine_sum is disabled")
        
# =============================================================================
##                                                                      The END 
# =============================================================================
