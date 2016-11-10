#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file
#  Simple "value with error" utility
#  @see Ostap::Math::ValueWithError
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
""" Simple ``value with error'' utility
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = "Version$Revision$"
# =============================================================================
__all__     = (
    'VE'  ,  ## Value with error  
    'VVE' ,  ## vector of Values with Errors 
    ) 
# =============================================================================
import ROOT, cppyy
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.ve' )
else                       : logger = getLogger ( __name__        )
# =============================================================================

cpp       = cppyy.gbl
std       = cpp.std 
Ostap     = cpp.Ostap

VE        = Ostap.Math.ValueWithError
VVE       = std.vector ( VE )

VE.Vector = VVE
VE.Vector .__str__   = lambda s : str( [ i for i in s ])
VE.Vector .__repr__  = lambda s : str( [ i for i in s ])

VVVE = std.vector( VVE )
VVE.Vector = VVVE 
VVE.Vector . __str__  = lambda s : str( [ i for i in s ] )
VVE.Vector . __repr__ = lambda s : str( [ i for i in s ] )
VVE.Vector . __len__  = lambda s : s.size ()



## Sum the contents of the vector
def _ve_sum_ ( s ) :
    """Sum the contents of the vector
    >>> v = ...
    >>> s = v.sum()
    """
    return Ostap.Math.sum ( s )


## Sum the contents of the vector
def _ve_asum_ ( s ) :
    """Sum the contents of the vector
    >>> v = ...
    >>> s = v.abssum()
    """
    return Ostap.Math.abssum ( s )

_ve_sum_  . __doc__ += '\n' + Ostap.Math.sum    .__doc__
_ve_asum_ . __doc__ += '\n' + Ostap.Math.abssum .__doc__


Ostap.Math.SVector2WithError = Ostap.Math.SVectorWithError(2,'double')
Ostap.Math.SVector3WithError = Ostap.Math.SVectorWithError(3,'double')
Ostap.Math.SVector4WithError = Ostap.Math.SVectorWithError(4,'double')
Ostap.Math.SVector5WithError = Ostap.Math.SVectorWithError(5,'double')
Ostap.Math.SVector6WithError = Ostap.Math.SVectorWithError(6,'double')
Ostap.Math.SVector8WithError = Ostap.Math.SVectorWithError(8,'double')

Ostap.Math.SVector2WithError  . __len__ = lambda s : 2 
Ostap.Math.SVector3WithError  . __len__ = lambda s : 3 
Ostap.Math.SVector4WithError  . __len__ = lambda s : 4 
Ostap.Math.SVector5WithError  . __len__ = lambda s : 5 
Ostap.Math.SVector6WithError  . __len__ = lambda s : 6 
Ostap.Math.SVector8WithError  . __len__ = lambda s : 8 

for t in ( Ostap.Math.ValueWithError         ,
           Ostap.Math.Point3DWithError       ,
           Ostap.Math.Vector3DWithError      ,
           Ostap.Math.LorentzVectorWithError ,
           Ostap.Math.SVector2WithError      ,
           Ostap.Math.SVector3WithError      ,
           Ostap.Math.SVector4WithError      ,
           Ostap.Math.SVector5WithError      ,
           Ostap.Math.SVector6WithError      ,
           Ostap.Math.SVector8WithError      ) :
    if not hasattr ( t , '_new_str_' ) :
        t._new_str_ = t.toString
        t.__str__   = t.toString
        t.__repr__  = t.toString
        

# =============================================================================
## get the B/S estimate from the formula 
#  \f$ \sigma  = \fras{1}{S}\sqrt{1+\frac{B}{S}}\f$
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-10-15
def _b2s_ ( s )  :
    """Get B/S estimate from the equation:
    error(S) = 1/sqrt(S) * sqrt ( 1 + B/S)
    >>> v = ...
    >>> b2s = v.b2s() ## get B/S estimate
    """
    #
    c2 = s.cov2  ()
    #
    if s.value() <= 0  or c2 <= 0 : return VE(-1,0) 
    #
    return c2/s - 1 

# =============================================================================
## get the precision with some  error estimation 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-10-15
def _prec2_ ( s )  :
    """Get precision with ``some'' error estimate 
    >>> v = ...
    >>> p = v.prec ()     
    """
    if not hasattr ( s , 'value' ) :
        return _prec_ ( VE ( s , 0 ) )
    #
    c =       s.error ()
    #
    if     c <  0 or s.value() == 0  : return VE(-1,0)
    elif   c == 0                    : return VE( 0,0)
    #
    return c/abs(s) 


VE . b2s        = _b2s_
VE . prec       = _prec2_
VE . precision  = _prec2_


_is_le_    = Ostap.Math.LessOrEqual ( 'double' )()

# =============================================================================
## Comparison of ValueWithError object with other objects
#  @attention it is comparison by value only, errors are ignored 
def _ve_lt_ ( self , other ) :
    """Comparison of ValueWithError object with other objects
    >>> a = VE( ... )
    >>> print a < b 
    Attention: comparison by value only!
    """
    return float(self) < float(other)

# =============================================================================
## Comparison of ValueWithError object with other objects
#  @attention it is comparison by value only, errors are ignored 
def _ve_le_ ( self , other ) :
    """Comparison of ValueWithError object with other objects
    >>> a = VE( ... )
    >>> print a <= b 
    Attention: comparison by value only!
    """
    return _is_le_ ( float(self) , float(other) )

# =============================================================================
## Comparison of ValueWithError object with other objects
#  @attention it is comparison by value only, errors are ignored 
def _ve_gt_ ( self , other ) :
    """Comparison of ValueWithError object with other objects
    >>> a = VE( ... )
    >>> print a > b 
    Attention: comparison by value only!
    """
    return float(self) > float(other)

# =============================================================================
## Comparison of ValueWithError object with other objects
#  @attention it is comparison by value only, errors are ignored 
def _ve_ge_ ( self , other ) :
    """Comparison of ValueWithError object with other objects
    >>> a = VE( ... )
    >>> print a >= b 
    Attention: comparison by value only!
    """
    return _is_le_ ( float(other) , float(self) )
    
VE.__lt__ = _ve_lt_ 
VE.__le__ = _ve_le_ 
VE.__gt__ = _ve_gt_ 
VE.__ge__ = _ve_ge_ 


_is_equal_ = Ostap.Math.Equal_To    ( 'double' )()
_is_zero_  = Ostap.Math.Zero        ( 'double' )()


# =============================================================================
## Equality for ValueWithError objects
def _ve_eq_ ( self , other ) :
    """Equality for ValueWithError objects
    >>> a = VE( ... )
    >>> b = VE( ... )
    >>> print a == b 
    """
    if isinstance ( other , VE ) :
        v1 = self .value()
        v2 = other.value()
        return _is_equal_ ( v1 , v2 ) and _is_equal_ ( self.cov2() , other.cov2() )
    elif _is_zero_ ( self.cov2() )  :
        return _is_equal_ ( float ( self ) ,  float ( other )  ) 
    else :
        raise NotImplementedError,' Equality for %s and  %s is not implemented' % ( self, other) 

# =============================================================================
## inequality for ValueWithError objects
def _ve_ne_ ( self , other ) :
    """Inequality for ValueWithError objects
    >>> a = VE( ... )
    >>> b = VE( ... )
    >>> print a != b 
    """
    try: 
        return not self == other
    except NotImplemented,s :
        raise NotImplementedError,' Inequality for %s and  %s is not implemented' % ( self, other) 

VE . __eq__ = _ve_eq_
VE . __ne__ = _ve_ne_


# =============================================================================
try :
    from scipy.random import normal as _gauss
except ImportError :
    from random       import gauss  as _gauss
    
# =============================================================================
## get the (gaussian) random number according to parameters
#
#  @code
#    >>> v = ...  ## the number with error
#
#    ## get 100 random numbers 
#    >>> for i in range ( 0, 100 ) : print v.gauss()
#    
#    ## get only non-negative numbers 
#    >>> for j in range ( 0, 100 ) : print v.gauss( lambda s : s > 0 )
#
#  @endcode
#  @attention scipy is needed! 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-08-10
# 
def _ve_gauss_ ( s , accept = lambda a : True , nmax = 1000 ) : 
    """ Get the gaussian random number
    >>> v = ...  ## the number with error
    ## get 100 random numbers 
    >>> for i in range ( 0, 100 ) : print v.gauss()    
    ## get only non-negative numbers 
    >>> for j in range ( 0, 100 ) : print v.gauss( lambda s : s > 0 )
    
    Attention: scipy is needed! 
    """
    #
    if 0 >= s.cov2() or iszero ( s.cov2 () ) : return s.value() ## return
    #
    v = s.value ()
    e = s.error ()
    #
    for i in range( 0 , nmax ) :
        r = _gauss ( v , e ) 
        if accept ( r ) : return r

    logger.warning("Can'n generate proper random number %s" % s )
    return v  

# =============================================================================
try :
    from random                import poisson as _poisson 
except ImportError :
    from ostap.math.random_ext import poisson as _poisson 

# =============================================================================
## generate poisson random number according to parameters 
#  @code
#    >>> v = ...  ## the number with error
#
#    ## get 100 random numbers 
#    >>> for i in range ( 0, 100 ) : print v.poisson ( fluctuate = True )
#    
#    ## get only odd numbers 
#    >>> for j in range ( 0, 100 ) : print v.poisson ( fluctuate = True , accept = lambda s : 1 ==s%2 )
#
#    ## do not fluctuate the mean of poisson:    
#    >>> for j in range ( 0, 100 ) : print v.poisson ( fluctuate = False  )
#
#  @endcode
#  @attention scipy is needed! 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-08-10   
def _ve_poisson_ ( s , fluctuate , accept = lambda s : True ) :
    """Generate poisson random number according to parameters 
    >>> v = ...  ## the number with error    
    ## get 100 random numbers 
    >>> for i in range ( 0, 100 ) : print v.poisson()    
    ## get only odd numbers 
    >>> for j in range ( 0, 100 ) : print v.poisson ( accept = lambda s : 1 ==s%2 )    
    ## do not fluctuate the mean of poisson:    
    >>> for j in range ( 0, 100 ) : print v.poisson ( fluctuate = False  )

    Attention: scipy is needed! 
    """
    s = VE( s ) 
    v = s.value() 
    if v < 0 and not fluctuate :
        raise TypeError, 'Negative mean without fluctuations (1)'
    if v < 0 and s.cov2() <= 0 :
        raise TypeError, 'Negative mean without fluctuations (2)'

    e = s.error() 
    if v < 0 and abs(v) > 3 * e  :
        logger.warning ( "Very inefficient mean fluctuations: %s" % s ) 

    mu = v
    if fluctuate :
        mu = s.gauss ()
        while mu < 0 :
            mu = s.gauss ()
    
    return _poisson ( mu ) 

VE.gauss   = _ve_gauss_
VE.poisson = _ve_poisson_ 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.logger.line import line 
    logger.info ( __file__  + '\n' + line  ) 
    logger.info ( 80*'*'   )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 

    a = VE(100,100)
    b = VE(400,400)

    logger.info ( 'a=%s, b=%s' % ( a , b ) )
    
    logger.info ( 'a+b         %s' % ( a + b ) )
    logger.info ( 'a-b         %s' % ( a - b ) )
    logger.info ( 'a*b         %s' % ( a * b ) )
    logger.info ( 'a/b         %s' % ( a / b ) )
    logger.info ( 'a/(a+b)     %s' % ( a.frac ( b ) ) )
    logger.info ( '(a-b)/(a+b) %s' % ( a.asym ( b ) ) )
    
    logger.info ( 80*'*' )
    
    
# =============================================================================
# The END 
# =============================================================================
