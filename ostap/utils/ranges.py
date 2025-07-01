#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Various looping rnages:
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
#  
# =============================================================================
""" Various looping ranges 
- uniform
- logarithmic
- power-law
- Chebyshev 
- random 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    ##
    'VRange'    , ## simple looping from xmin to xmax in N-steps 
    'vrange'    , ## simple looping from xmin to xmax in N-steps
    ##
    'LRange'    , ## simple looping from xmin to xmax in N-steps in log-scale 
    'lrange'    , ## simple looping from xmin to xmax in N-steps in log-scale 
    'log_range' , ## simple looping from xmin to xmax in N-steps in log-scale
    ## 
    'CRange'    , ## simple looping from xmin to xmax in N-steps using chebshev nodes 
    'crange'    , ## simple looping from xmin to xmax in N-steps using chebshev nodes 
    ##
    ## 
    'PRange'    , ## simple looping from xmin to xmax in N-steps using non-uniform power-law poins 
    'prange'    , ## simple looping from xmin to xmax in N-steps using non-uniform power-law points 
    ##
    'RRange'    , ## simple looping from xmin to xmax in N-steps using random points 
    'rrange'    , ## simple looping from xmin to xmax in N-steps using random points 
    ##
    )
#
# =============================================================================
from   ostap.core.ostap_types import integer_types, num_types 
import math, array, random  
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.utils.utils' )
else                       : logger = getLogger ( __name__            )
# =============================================================================
## @class VRange
#  Helper looper over the values between vmin and vmax :
#  @code
#  for v in VRange ( vmin = 0 , vmax = 5 , n = 100 ) :
#  ... print ( v ) 
#  @endcode 
class VRange(object) :
    """ Helper looper over the values between vmin and vmax :
    >>> for v in VRange ( vmin = 0 , vmax = 5 , n = 100 ) :
    >>> ... print ( v ) 
    """
    def __init__ ( self , vmin , vmax , n = 100 , edges = True ) :
        
        assert isinstance ( n , integer_types ) and 0 < n,\
               'VRange: invalid N=%s/%s' % ( n  , type ( n ) ) 
        
        self.__vmin  = min ( vmin , vmax )
        self.__vmax  = max ( vmin , vmax ) 
        self.__n     = n
        self.__edges = True if edges else False  

    @property
    def edges ( self ) :
        """`edges`: include edges?"""
        return self.__edges
    
    @property
    def vmin  ( self ) :
        """`vmin' : minimal value"""
        return self.__vmin
    @property
    def vmax  ( self ) :
        """`vmax' : maximal value"""
        return self.__vmax    
    @property
    def n ( self )  :
        """`n' : number of steps"""
        return self.__n

    def __len__     ( self ) :

        n = self.__n
        e = self.__edges 
        return n + 1 if e else n - 1 

    def __iter__    ( self ) :
        
        n  = self.n 
        fn = 1.0 / float ( n ) 
        e  = self.edges

        vmn = self.vmin
        vmx = self.vmax
        
        if e : yield vmn

        for i in range ( 1 , n ) :
            f2 = i * fn
            f1 = 1 - f2
            yield vmn * f1 + f2 * vmx
            
        if e : yield vmx
                
# =============================================================================
## loop over values between xmin and xmax 
#  @code
#  for x in vrange ( xmin , xmax , 200 ) :
#         print (x) 
#  @endcode
def vrange ( vmin , vmax , n = 100 , edges = True ) :
    """ Loop  over range of values between xmin and xmax 
    >>> for v in vrange ( vmin , vmax , 200 ) :
    ...                print (v) 
    """
    return VRange ( vmin , vmax , n , edges )


# =============================================================================
## @class LRange
#  Helper looper over the values between vmin and vmax using log-steps 
#  @code
#  for v in LRange ( vmin = 1 , vmax = 5 , n = 100 ) :
#  ... print ( v ) 
#  @endcode 
class LRange(VRange) :
    """ Helper looper over the values between vmin and vmax using log-steps
    >>> for v in LRange ( vmin = 1 , vmax = 5 , n = 100 ) :
    >>> ... print ( v ) 
    """
    def __init__ ( self , vmin , vmax , n = 100 , edges = True ) :
        
        assert 0 < vmin  and 0 < vmax,\
           'LRange: invalid  non-positive vmin/ymax values: %s/%s' %  ( vmin , vmax )

        super ( LRange , self ).__init__ ( vmin , vmax , n , edges ) 

        self.__lmin = math.log10 ( self.vmin )
        self.__lmax = math.log10 ( self.vmax )
                
    @property
    def lmin  ( self ) :
        """``lmin'' : log10(minimal value)"""
        return self.__lmin
    @property
    def lmax  ( self ) :
        """``lmax'' : log10(maximal value)"""
        return self.__lmax    

    def __iter__    ( self ) :

        n  = self.n 
        fn = 1.0 / float ( n )
        e  = self.edges

        lmn = self.__lmin
        lmx = self.__lmax
        
        if e : yield self.vmin
        
        for i in range ( 1 , n  ) :
            #
            f2 = i * fn
            f1 = 1 - f2
            yield 10.0 ** ( lmn * f1 + f2 * lmx ) 
            
        if e : yield self.vmax

# =============================================================================
## loop over values between xmin and xmax in log-scale 
#  @code
#  for x in log_range ( xmin , xmax , 200 ) :
#         print (x) 
#  @endcode
def log_range ( vmin , vmax , n = 100 , edges = True ) :
    """Loop over values between xmin and xmax in log-scale 
    >>> for x in log_range ( xmin , xmax , 200 ) :
    >>>      print (x) 
    """
    return LRange ( vmin , vmax , n , edges )

# =============================================================================
## loop over values between xmin and xmax in log-scale 
#  @code
#  for v in lrange ( vmin , vmax , 200 ) : ## ditto 
#         print (v) 
#  @endcode
def lrange ( vmin , vmax , n = 100 , edges = True ) :
    """ Loop over values between vmin and vmax in log-scale 
    >>> for v in lrange ( vmin , vmax , 200 ) :  ## ditto 
    >>>      print (v) 
    """
    return LRange ( vmin , vmax , n , edges )

# =============================================================================
## @class CRange
#  Generate sequence of numbers between vmin and vmax according to Chebyshev nodes
#  It can be useful for e.g. interpolation nodes
#  @code
#  for c in CRange(-1,1,10) : print ( c ) 
#  @endcode 
class CRange(VRange):
    """ Generate sequence of numbers between vmin and vmax accrording to Chebyshev nodes
    It can be useful for e.g. interpolation nodes
    >>> for c in CRange(-1,1,10) : print ( c ) 
    """
    __nodes = {}
    
    ## number of nodes 
    def __len__     ( self ) : return self.n 
    ## 
    def __iter__    ( self ) :

        n     = self.n
        
        if not n in self.__nodes :
            nodes = array.array ( 'd' , n * [ 0.0 ] )
            n2    = n // 2 
            pn    = math.pi / ( 2 * n )
            for k in range ( n2 ) :
                vv  = math.cos ( ( 2 * k + 1 ) * pn )
                nodes [     k     ] =  vv
                nodes [ n - k - 1 ] = -vv
            self.__nodes [ n ] = nodes 

        nodes = self.__nodes [ n ]        
        mid   = 0.5 * ( self.vmin + self.vmax )
        scale = 0.5 * ( self.vmin - self.vmax )              
        for x in nodes :             
            yield mid + scale * x

# =============================================================================
### Generate sequence of numbers between vmin and vmax accrording to Chebyshev nodes
#  It can be useful for e.g. interpolation nodes
#  @code
#  for c in crange(-1,1,10) : print ( c ) 
#  @endcode 
def crange ( vmin , vmax , n = 10 ) :
    """ Generate sequence of numbers between vmin and vmax accrording to Chebyshev nodes
    It can be usefor for e.g. interpolation nodes
    >>> for c in crange(-1,1,10) : print ( c ) 
    """
    return CRange ( vmin , vmax , n )

# =============================================================================
## @class PRange
#  Helper looper over the values between vmin and vmax with nonuniform power-law) distributed poinnts 
#  @code
#  for v in PRange ( vmin = 1 , vmax = 5 , n = 100 , power = 2 ) , edges = True :
#  ... print ( v ) 
#  @endcode 
class PRange(VRange) :
    """ Helper looper over the values between vmin and vmax with non-unifomr (power-low) points 
    >>> for v in PRange ( vmin = 1 , vmax = 5 , n = 100 , power = 2 , edges = True ) :
    >>> ... print ( v ) 
    """
    def __init__ ( self , vmin , vmax , n = 100 , power = 2 , edges = True ) :
        
        assert isinstance ( power , num_types ) and 0 < power , 'PRange: Invalid powr %s' % power  

        super ( PRange , self ).__init__ ( vmin , vmax , n , edges ) 

        self.__power = power

    @property
    def power ( self ) :
        """`power`: the power-low exponent """
        return self.__power
    
    def __iter__    ( self ) :

        n     = self.n 
        fn    = 1.0 / float ( n ) 
        e     = self.edges
        
        p     = self.power
        vmn   = self.vmin
        vmx   = self.vmax
        delta = vmx - vmn
        
        if e : yield vmn

        for i in range ( 1 , n ) :
            x = i * fn
            yield vmn + delta * ( x ** p ) 

        if e : yield vmx

# =============================================================================
## Loop over sequence of non-uniformly distributed  (power-low) unmpers
#  @code
#  for c in prange(-1,1,10, power = 2 , edges = True ) : print ( c ) 
#  @endcode 
def prange ( vmin , vmax , n = 10 , power = 2 , edges = True  ) :
    """ Loop over the sequence of non-unifrmy distribited (power-low) numbers 
    >>> for c in prange(-1,1,10, power = 2 , edges = True  ) : print ( c ) 
    """
    return PRange ( vmin , vmax , n , power = power  , edges = edges )

# =============================================================================
## @class RRange
#  Helper looper over the random values between vmin and vmax
#  @code
#  for v in RRange ( vmin = 1 , vmax = 5 , n = 100 , edges = True  ) :
#  ... print ( v ) 
#  @endcode 
class RRange(VRange) :
    """ Helper looper over the values between vmin and vmax using random-steps
    >>> for v in RRange ( vmin = 1 , vmax = 5 , n = 100 , edges = True ) :
    >>> ... print ( v ) 
    """
    def __iter__    ( self ) :
        
        n   = self.n
        vmn = self.vmin
        vmx = self.vmax
        e   = self.edges

        if e : yield vmn
        
        for i in range ( 1 , n ) :
            yield random.uniform ( vmn, vmx )

        if e : yield vmx
            
# =============================================================================
## Generate sequence of random numbers between vmin and vmax
#  @code
#  for c in rrange(-1,1,10) : print ( c ) 
#  @endcode 
def rrange ( vmin , vmax , n = 10 , edges = True  ) :
    """ Generate random sequence of numbers between vmin and vmax 
    >>> for c in crange(-1,1,10, edges = True ) : print ( c ) 
    """
    return RRange ( vmin , vmax , n , edges )


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
