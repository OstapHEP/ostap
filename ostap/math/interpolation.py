#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/interpolation.py
#  Module with some useful utilities for dealing with interpolation.
#
#  In particular, it providies very efficient Barycentric Lagrange interpolation
#  - it takes O(n)   flops for initialization with Chebyshev/Lobatto or Uniform abscissas 
#  - it takes O(n^2) flops for initialization with arbitrary interpolation      abscissas
#  - it takes O(n)   flops for evaluation! It is very fast!
#
#  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
#       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
#       ISSN (print): 0036-1445
#       ISSN (online): 1095-7200
#  @see https://doi.org/10.1137/S0036144502417715
#  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
#  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
#  @see Ostap::Math::Barycentric
#  @see Ostap::Math::Newton
#
#  Straightforward Lagrange and Neville algorithm are also provided
#  - Both are rather slow: O(n^2) flops for evaluation,  while Neville is a bit faster 
#  - Lagrange algorithm is not stable numerically,  and Neville algorithm is more stable
#
#  @see Ostap::Math::Interpolation
#  @see Ostap::Math::Neville 
#  @see Ostap::Math::Lagrange 
#  @see Ostap::Math::Newton 
#  @see Ostap::Math::Barycentric
#
#  For  completeness see also:
#  - interpolation with Bersntein polynomials using on Newton-Bernstein algorithm
#  - interpolation with B-splines
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2018-07-22
# =============================================================================
"""Useful utilities for dealing with interpolation.

In particular it provides very efficient ``Barycentric Lagrange interpolation'':

- it takes O(n)   flops for initialization with Chebyshev/Lobatto or uniform abscissas 
- it takes O(n^2) flops for initialization with arbitrary interpolation      abscissas
- it takes O(n)   flops for evaluation! It is very fast! 

- see Jean-Paul Berrut and Lloyd N. Trefethen, 
...   Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
...   ISSN (print): 0036-1445
...   ISSN (online): 1095-7200
- see https://doi.org/10.1137/S0036144502417715
- see https://en.wikipedia.org/wiki/Lagrange_polynomial
- see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf

Straightforward Lagrange and Neville algorithms are also provided:
- both are rather slow: O(n^2) flops for evaluation,  while Neville is a bit faster 
- Lagrange algorithm is not stable numerically, while Neville algorithm is more stable

For completeness see also:
- interpolation with bersntein polynomials using on Newton-Bernstein algorithm
- interpolation with B-splines 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2018-07-22"
__all__     = (
    #
    'interpolate'           , ## Very  efficient Barycentric Lagrange interpolation
    'lagrange'              , ## ditto 
    'barycentric'           , ## ditto
    #
    'points'                , ## helper function to create the collection of interpolation points
    'Abscissas'             , ## collecion of interpolation abscissas 
    'Table'                 , ## the interpolation table (collection of points)
    # for completeness  
    'interpolate_bernstein' , ## Newton-Bernstein interpolation 
    'interpolate_bspline'   , ## Basic spline interpolation
    )
# =============================================================================
import  ROOT, math
from    builtins          import range
from    ostap.core.core   import cpp, Ostap
from    ostap.core.types  import is_integer
from    ostap.math.base   import iszero, isequal, doubles 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.interpolation' )
else                       : logger = getLogger ( __name__                   )
# =============================================================================

# =============================================================================
## Interpolation abscissas
# ============================================================================= 

# ============================================================================= 
## printout for the interpolation abscissas
#  @code
#  print a 
#  @endcode 
def _a_str_ ( self , nmax = 8 ) :
    """Printout for interpolation Abscissas
    >>> print a  
    """
    _n = self.n()
    if _n <= nmax : return 'Abscissas(%s)' % self.x()
    ##
    n2 = nmax/2
    s1 = ', '.join( ( str ( x ) for x in self.x() [    : n2 ] ) )
    s2 = ', '.join( ( str ( x ) for x in self.x() [ -1 :    ] ) )
    return 'Abscissas(n=%d,[%s, ... , %s])' % ( self.n() , s1 , s2 ) 
# =======================================================================================
## iterator over interpolation abscissas
#  @code
#  a =  ...
#  for x in a : print x  
#  @endcode 
def _a_iter_ ( self ) :
    """ Iterator over interpolation abscissas
    >>> a =  ...
    >>> for x in a : print x  
    """
    N = self.n()
    for i in range ( N ) :
        yield self.x ( i ) 
# =============================================================================
## Iterator over (index,x) pairs
#  a =  ...
#  for i,x in a.iteritems() : print i,x
#  @endcode 
def _a_iteritems_ ( self ) :
    """Iterator over (index,x) pairs
    >>> a =  ...
    >>> for i,x in a.iteritems() : print i,x
    """
    N = self.n()
    for i in range ( N ) :
        yield i , self.x ( i ) 
# =============================================================================
## Get the item or slice 
#  @code
#  a  = ...
#  x1 = a[1]
#  a1 = a[:10] 
#  @endcode d
def _a_getitem_ ( self , i ) :
    """Get the item or slice 
    >>> a  = ...
    >>> x1 = a[1]
    >>> a1 = a[:10]
    """
    if   isinstance ( i , int   ) :
        
        if 0 <= i < self.n () : return self.x ( i )
        raise IndexError ('Invalid key %s' % i )
    
    elif isinstance ( i , slice ) :
        
        start , stop , step = i.indices ( self.n() )
        if   1 == step : return self.slice ( start , stop )
        _x = self.x() [i] 
        if 0 <  step : return Ostap.Math.Interpolation.Abscissas ( _xi , True )
        else         : return Ostap.Math.Interpolation.Abscissas ( _xi        )
        
    raise TypeError ('Invalid key/value %s' % i )

# =============================================================================
## remove point from  the   abscissas
#  @code
#  a = ...
#  del a[2]
#  @endcode 
def _a_delitem_ (   self , i ) : 
    """Remove point from abscissas
    >>> a = ...
    >>> del a[2]
    """
    if   isinstance ( i , int   ) :
        if 0 <= i < self.n () : return self.remove ( i )
        raise IndexError ('Invalid key %s' % i )
    
    raise TypeError ('Invalid key/value %s' % i )

Ostap.Math.Interpolation.Abscissas.__str__      = _a_str_
Ostap.Math.Interpolation.Abscissas.__repr__     = _a_str_
Ostap.Math.Interpolation.Abscissas.__iter__     = _a_iter_
Ostap.Math.Interpolation.Abscissas.__len__      = lambda s : s.n()
Ostap.Math.Interpolation.Abscissas.iteritems    = _a_iteritems_ 
Ostap.Math.Interpolation.Abscissas.__getitem__  = _a_getitem_  
Ostap.Math.Interpolation.Abscissas.__delitem__  = _a_delitem_  


# =============================================================================
## Interpolation points 
# ============================================================================= 

# =============================================================================
## printout for the interpolation points
#  @code
#  print a 
#  @endcode 
def _p_str_ ( self , nmax = 7 ) :
    """Printout for interpolation Table
    >>> print a  
    """
    _n = self.n()
    if _n <= nmax :
        s = ', '.join ( ( "%s: %s" %  (x,y)  for x,y in self ) ) 
        return 'Table({%s})' % s 
    ##
    n2 = nmax/2

    s1 = ', '.join ( ( '%s: %s' % self[i] for i in  range ( n2 ) ) ) 
    s2 =               '%s: %s' % self[ self.n() -1 ]
    
    return 'Table(n=%d,{%s, ... , %s})' % ( self.n() , s1 , s2 )
                     
# =======================================================================================
## iterator over interpolation points 
#  @code
#  a =  ...
#  for x,y in a : print x,y  
#  @endcode 
def _p_iter_ ( self ) :
    """ Iterator over interpolation points 
    >>> a =  ...
    >>> for x,y in a : print x,y  
    """
    N = self.n()
    for i in range ( N ) :
        yield self.x ( i ) , self.y ( i )  
# =============================================================================
## Iterator over (index,x) pairs
#  a =  ...
#  for i,p in a.iteritems() : print i,p
#  @endcode 
def _p_iteritems_ ( self ) :
    """Iterator over (index,point) pairs
    >>> a =  ...
    >>> for i,p in a.iteritems() : print i,p
    """
    N = self.n()
    for i in range ( N ) :
        yield i , ( self.x ( i ) , self.y ( i ) )  

# =============================================================================
## Get the item or slice 
#  @code
#  a  = ...
#  x1 = a[1]
#  a1 = a[:10] 
#  @endcode d
def _p_getitem_ ( self , i ) :
    """Get the item or slice 
    >>> a  = ...
    >>> x1 = a[1]
    >>> a1 = a[:10]
    """
    if   is_integer ( i ) :
        
        if 0 <= i < self.n () :
            return self.x ( i ) , self.y ( i )
        
        raise IndexError ('Invalid key %s' % i )
    
    elif isinstance ( i , slice ) :
        
        start , stop , step = i.indices ( self.n() )
        if   1 == step : return self.slice ( start , stop )
        _x = self.x() [i] 
        _y = self.y() [i] 
        return Ostap.Math.Interpolation.Table ( _x , _y )
        
    raise TypeError ('Invalid key type/value %s' % i )

# =============================================================================
## remove point from the point collection
#  @code
#  a = ...
#  del a[2]
#  @endcode 
def _p_delitem_ (   self , i ) : 
    """Remove point from abscissas
    >>> a = ...
    >>> del a[2]
    """
    if   isinstance ( i , int   ) :
        if 0 <= i < self.n () : return self.remove ( i )
        raise IndexError ('Invalid key %s' % i )
    
    raise TypeError ('Invalid key/value %s' % i )

Ostap.Math.Interpolation.Table.__str__      = _p_str_
Ostap.Math.Interpolation.Table.__repr__     = _p_str_
Ostap.Math.Interpolation.Table.__iter__     = _p_iter_
Ostap.Math.Interpolation.Table.__len__      = lambda s : s.n()
Ostap.Math.Interpolation.Table.iteritems    = _p_iteritems_ 
Ostap.Math.Interpolation.Table.__getitem__  = _p_getitem_  
Ostap.Math.Interpolation.Table.__delitem__  = _p_delitem_  

# =============================================================================
## Weights for barycentric Lagrange interpolation
# ============================================================================= 

# =============================================================================
## printout for the interpolation points
#  @code
#  print a 
#  @endcode 
def _w_str_ ( self , nmax = 7 ) :
    """Printout for interpolation Weights 
    >>> print a  
    """
    _n = self.n()
    if _n <= nmax :
        s = ', '.join ( ( "%s: %s" %  (x,y)  for x,y in self ) ) 
        return 'Weights({%s})' % s 
    ##
    n2 = nmax/2

    s1 = ', '.join ( ( '%s: %s' % self[i] for i in  range ( n2 ) ) ) 
    s2 =               '%s: %s' % self[ self.n() -1 ]
    
    return 'Weights(n=%d,{%s, ... , %s})' % ( self.n() , s1 , s2 )
                     
# =======================================================================================
## iterator over interpolation points 
#  @code
#  a =  ...
#  for x,w in a : print x,w
#  @endcode 
def _w_iter_ ( self ) :
    """ Iterator over interpolation weights 
    >>> a =  ...
    >>> for x,w in a : print x,w  
    """
    N = self.n()
    for i in range ( N ) :
        yield self.x ( i ) , self.w ( i )  
# =============================================================================
## Iterator over (i,(x,w)) triplets 
#  a =  ...
#  for i,p in a.iteritems() : print i,p
#  @endcode 
def _w_iteritems_ ( self ) :
    """Iterator over (index,point) pairs
    >>> a =  ...
    >>> for i,p in a.iteritems() : print i,p
    """
    N = self.n()
    for i in range ( N ) :
        yield i , ( self.x ( i ) , self.w( i ) )  

# =============================================================================
## Get the item or slice 
#  @code
#  a  = ...
#  x1 = a[1]
#  a1 = a[:10] 
#  @endcode d
def _w_getitem_ ( self , i ) :
    """Get the item or slice 
    >>> a  = ...
    >>> x1 = a[1]
    >>> a1 = a[:10]
    """
    if   is_integer ( i ) :
        
        if 0 <= i < self.n () :
            return self.x ( i ) , self.w ( i )
        
        raise IndexError ('Invalid key %s' % i )
    
    elif isinstance ( i , slice ) :
        
        start , stop , step = i.indices ( self.n() )
        if   1 == step : return self.slice ( start , stop )
        _x = self.x() [i] 
        return Ostap.Math.Interpolation.Weights ( _x )
        
    raise TypeError ('Invalid key type/value %s' % i )

# =============================================================================
## remove the point from  the point collection
#  @code
#  a = ...
#  del a[2]
#  @endcode 
def _w_delitem_ (   self , i ) : 
    """Remove point from weights 
    >>> a = ...
    >>> del a[2]
    """
    if   isinstance ( i , int   ) :
        if 0 <= i < self.n () : return self.remove ( i )
        raise IndexError ('Invalid key %s' % i )
    
    raise TypeError ('Invalid key/value %s' % i )

Ostap.Math.Interpolation.Weights.__str__      = _w_str_
Ostap.Math.Interpolation.Weights.__repr__     = _w_str_
Ostap.Math.Interpolation.Weights.__iter__     = _w_iter_
Ostap.Math.Interpolation.Weights.__len__      = lambda s : s.n()
Ostap.Math.Interpolation.Weights.iteritems    = _w_iteritems_ 
Ostap.Math.Interpolation.Weights.__getitem__  = _w_getitem_  
Ostap.Math.Interpolation.Weights.__delitem__  = _w_delitem_  


# ==================================================================================
## Neville & Lagrange interpolants
# ==================================================================================
## printout for Neville interpolant 
#  @code
#  print a 
#  @endcode 
def _n_str_ ( self , nmax = 7 ) :
    """Printout for Neville interpolant
    >>> print a  
    """
    _n = self.n()
    if _n <= nmax :
        s = ', '.join ( ( "%s: %s" %  (x,y)  for x,y in self ) ) 
        return 'Neville({%s})' % s 
    ##
    n2 = nmax/2

    s1 = ', '.join ( ( '%s: %s' % self[i] for i in  range ( n2 ) ) ) 
    s2 =               '%s: %s' % self[ self.n() -1 ]
    
    return 'Neville(n=%d,{%s, ... , %s})' % ( self.n() , s1 , s2 )
# =============================================================================
## printout for Largange interpolant 
#  @code
#  print a 
#  @endcode 
def _l_str_ ( self , nmax = 7 ) :
    """Printout for Lagrange interpolant
    >>> print a  
    """
    _n = self.n()
    if _n <= nmax :
        s = ', '.join ( ( "%s: %s" %  (x,y)  for x,y in self ) ) 
        return 'Lagrange({%s})' % s 
    ##
    n2 = nmax/2

    s1 = ', '.join ( ( '%s: %s' % self[i] for i in  range ( n2 ) ) ) 
    s2 =               '%s: %s' % self[ self.n() -1 ]
    
    return 'Lagrange(n=%d,{%s, ... , %s})' % ( self.n() , s1 , s2 )


Ostap.Math.Neville .__str__      = _n_str_
Ostap.Math.Neville .__repr__     = _n_str_
Ostap.Math.Lagrange.__str__      = _l_str_
Ostap.Math.Lagrange.__repr__     = _l_str_

# ==================================================================================
## Barycentric Lagrange interpolant 
# ==================================================================================

# =============================================================================
## printout for Barycentric Largange interpolant 
#  @code
#  print a 
#  @endcode 
def _b_str_ ( self , nmax = 7 ) :
    """Printout for Lagrange interpolant
    >>> print a  
    """
    _n = self.n()
    if _n <= nmax :
        s = ', '.join ( ( "%s: %s" %  (x,y)  for x,y in self ) ) 
        return 'Barycentric({%s})' % s 
    ##
    n2 = nmax/2

    s1 = ', '.join ( ( '%s: %s' % self[i] for i in  range ( n2 ) ) ) 
    s2 =               '%s: %s' % self[ self.n() -1 ]
    
    return 'Barycentric(n=%d,{%s, ... , %s})' % ( self.n() , s1 , s2 )


# =======================================================================================
## iterator over interpolation points 
#  @code
#  a =  ...
#  for x,y in a : print x,y  
#  @endcode 
def _b_iter_ ( self ) :
    """ Iterator over interpolation points 
    >>> a =  ...
    >>> for x,y in a : print x,y  
    """
    N = self.n()
    for i in range ( N ) :
        yield self.x ( i ) , self.y ( i )  
# =============================================================================
## Iterator over (index,x) pairs
#  a =  ...
#  for i,p in a.iteritems() : print i,p
#  @endcode 
def _b_iteritems_ ( self ) :
    """Iterator over (index,point) pairs
    >>> a =  ...
    >>> for i,p in a.iteritems() : print i,p
    """
    N = self.n()
    for i in range ( N ) :
        yield i , ( self.x ( i ) , self.y ( i ) )  

# =============================================================================
## Get the item or slice 
#  @code
#  a  = ...
#  x1 = a[1]
#  a1 = a[:10] 
#  @endcode d
def _b_getitem_ ( self , i ) :
    """Get the item or slice 
    >>> a  = ...
    >>> x1 = a[1]
    >>> a1 = a[:10]
    """
    if   is_integer ( i ) :
        
        if 0 <= i < self.n () :
            return self.x ( i ) , self.y ( i )
        
        raise IndexError ('Invalid key %s' % i )
    
    elif isinstance ( i , slice ) :
        
        start , stop , step = i.indices ( self.n() )
        if   1 == step : return self.slice ( start , stop )
        _x = self.x() [i] 
        _y = self.y() [i] 
        return Ostap.Math.Interpolation.Barycentric( _x , _y )
        
    raise TypeError ('Invalid key type/value %s' % i )

# =============================================================================
## remove point from  the point collection
#  @code
#  a = ...
#  del a[2]
#  @endcode 
def _b_delitem_ (   self , i ) : 
    """Remove point from abscissas
    >>> a = ...
    >>> del a[2]
    """
    if   isinstance ( i , int   ) :
        if 0 <= i < self.n () : return self.remove ( i )
        raise IndexError ('Invalid key %s' % i )
    
    raise TypeError ('Invalid key/value %s' % i )


Ostap.Math.Barycentric.__str__      = _b_str_
Ostap.Math.Barycentric.__repr__     = _b_str_
Ostap.Math.Barycentric.__iter__     = _b_iter_
Ostap.Math.Barycentric.__len__      = lambda s : s.n()
Ostap.Math.Barycentric.iteritems    = _b_iteritems_ 
Ostap.Math.Barycentric.__getitem__  = _b_getitem_  
Ostap.Math.Barycentric.__delitem__  = _b_delitem_  



Abscissas = Ostap.Math.Interpolation.Abscissas
Table     = Ostap.Math.Interpolation.Table

# ==================================================================================
## Construct very  efficient Barycentric Largange polynomial 
#  @code
#
#  b1 = interpolate ( lambda x : x*x , [0,0.5,1,2]  )  
#  b2 = interpolate ( { 0:0 , 0.5:0.25 , 1:1 }      )  
#  b3 = interpolate ( [0,0.25,1,4] , [ 0,0.5, 1,2]  )  
#  b4 = interpolate ( lambda x : x*x ,  Abscissas( 4 , -2 , 2 , 1 ) )
#  
#  b1 = lagrange    ( lambda x : x*x , [0,0.5,1,2] )  
#  b2 = lagrange    ( { 0:0 , 0.5:0.25 , 1:1 }     )  
#  b3 = lagrange    ( [0,0.25,1,4] , [ 0,0.5, 1,2] )    
#  b4 = lagrange    ( lambda x : x * x , Abscissas( 4 , -2 , 2 , 1 ) )
#
#  b1 = barycentric ( lambda x : x*x , [0,0.5,1,2] )  
#  b2 = barycentric ( { 0:0 , 0.5:0.25 , 1:1 }     )  
#  b3 = barycentric ( [0,0.25,1,4] , [ 0,0.5, 1,2] )    
#  b4 = barycentric ( lambda x : x * x , Abscissas( 4 , -2 , 2 , 1 ) )  
#  @endcode 
#  @param func  the function or list of function values
def lagrange ( func , abscissas  = None ) :
    """Construct very efficient Barycentric Lagrange interpolant
    
    - it takes O(n)   flops for initialization with Chebyshev/Lobatto or uniform abscissas 
    - it takes O(n^2) flops for initialization with arbitrary interpolation      abscissas
    - it takes O(n)   flops for evaluation! It is very fast! 
    
    - see Jean-Paul Berrut and Lloyd N. Trefethen, 
    ...   Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
    ...   ISSN (print): 0036-1445
    ...   ISSN (online): 1095-7200
    - see https://doi.org/10.1137/S0036144502417715
    - see https://en.wikipedia.org/wiki/Lagrange_polynomial
    - see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
    
    func      : the ``function''
    abscissas : abscissas

    :Example:
    
    >>> b1 = interpolate ( lambda x : x*x , [0,0.5,1,2] )  
    >>> b2 = interpolate ( { 0:0 , 0.5:0.25 , 1:1 }     )  
    >>> b3 = interpolate ( [0,0.25,1,4] , [ 0,0.5, 1,2] )    
    >>> b4 = interpolate ( lambda x : x * x , Abscissas( 4 , -2 , 2 , 1 ) )
    
    >>> b1 = lagrange    ( lambda x : x*x , [0,0.5,1,2] )  
    >>> b2 = lagrange    ( { 0:0 , 0.5:0.25 , 1:1 }     )  
    >>> b3 = lagrange    ( [0,0.25,1,4] , [ 0,0.5, 1,2] )    
    >>> b4 = lagrange    ( lambda x : x * x , Abscissas( 4 , -2 , 2 , 1 ) )  

    >>> b1 = barycentric ( lambda x : x*x , [0,0.5,1,2] )  
    >>> b2 = barycentric ( { 0:0 , 0.5:0.25 , 1:1 }     )  
    >>> b3 = barycentric ( [0,0.25,1,4] , [ 0,0.5, 1,2] )    
    >>> b4 = barycentric ( lambda x : x * x , Abscissas( 4 , -2 , 2 , 1 ) )  
    
    """
       
    from types       import GeneratorType as GT
    from collections import Iterable      as IT
    from collections import Mapping       as MT

    ## switch on abscissas:
    _W = isinstance ( abscissas , Ostap.Math.Interpolation.Weights   )
    _A = isinstance ( abscissas , Ostap.Math.Interpolation.Abscissas )
    
    if   isinstance ( abscissas , GT ) : abscissas = [ x for x in abscissas ]
    elif abscissas is None             :
        if   isinstance ( func , MT )  :
            _new_abscissas = []
            _new_func      = []
            _keys          = func.keys()
            for _k in sorted ( _keys ) :
                _new_abscissas.append (      _k  )
                _new_func     .append ( func[_k] )
            abscissas = _new_abscissas 
            func      = _new_func
        elif isinstance ( func , Ostap.Math.Interpolation.Table ) :
            return Ostap.Math.Barycentric ( func )                   ## ready
        else  :
            raise TypeError ( "For ``None''-abscissas ``func'' must be Table or Mapping!" ) 

    ## switch on func         
    if   callable   ( func ) :        
        nfunc = func        
        if _W or _A : func = ( nfunc ( x ) for x in abscissas.x() )
        else        : func = ( nfunc ( x ) for x in abscissas     )
    elif isinstance ( func , GT ) : pass                         ## generator 
    elif isinstance ( func , IT ) : pass                         ## iterable
    else :
        raise TypeError("Can't treat ``func''=%s"  %  func )
    ##
    from ostap.math.base import doubles
    if _W or _A : return Ostap.Math.Barycentric ( abscissas , doubles ( func ) )
    ##
    return Ostap.Math.Barycentric ( doubles ( abscissas ) , doubles ( func ) )

# =============================================================================
## make new name 
interpolate = lagrange
barycentric = lagrange 

# ==================================================================================
## Construct collection of interpolated points 
#  @code
#  p1 = points ( lambda x : x*x , [0,0.5,1,2]  )  
#  p2 = points ( { 0:0 , 0.5:0.25 , 1:1 }      )  
#  b3 = points ( [0,0.25,1,4] , [ 0,0.5, 1,2]  )  
#  b4 = points ( lambda x : x * x , Abscissas ( 4 , -2 , 2 , 1 ) )  
#  @endcode 
#  @param func  the function or list of function values
def points ( func , abscissas  = None ) :
    """Construct collection of interpolation points
    
    func      : the ``function''
    abscissas : abscissas
    
    :Example:
    
    >> b1 = points ( lambda x : x*x , [0,0.5,1,2] )  
    >> b2 = points ( { 0:0 , 0.5:0.25 , 1:1 }     )  
    >> b3 = points ( [0,0.25,1,4] , [ 0,0.5, 1,2] )    
    >> b4 = points ( lambda x : x * x , Abscissas ( 4 , -2 , 2 , 1 ) )  
    
    """
       
    from types       import GeneratorType as GT
    from collections import Iterable      as IT
    from collections import Mapping       as MT

    ## switch on abscissas:
    _A = isinstance ( abscissas , Ostap.Math.Interpolation.Abscissas )
    
    if   isinstance ( abscissas , GT ) : abscissas = [ x for x in abscissas ]
    elif abscissas is None             :
        if   isinstance ( func , MT )  :
            _new_abscissas = []
            _new_func      = []
            _keys          = func.keys()
            for _k in sorted ( _keys ) :
                _new_abscissas.append (      _k  )
                _new_func     .append ( func[_k] )
            abscissas = _new_abscissas 
            func      = _new_func
        elif isinstance ( func , Ostap.Math.Interpolation.Table ) : return func
        elif hasattr ( func , 'iteritems' ) : 
            _new_abscissas = []
            _new_func      = []
            for _a,_f in func.iteritems() :
                _new_abscissas.append ( _a )
                _new_func     .append ( _f )
            abscissas = _new_abscissas 
            func      = _new_func
        elif isinstance ( func , IT ) :
            try :
                _new_abscissas = []
                _new_func      = []
                for _a , _f in func :
                    _new_abscissas.append ( _a  )
                    _new_func     .append ( _f  )
                abscissas = _new_abscissas
                func      = _new_func
            except :
                raise TypeError ( "Invalid iterable ``func'' for ``None''-abscissas!" )                 
        else  :
            raise TypeError ( "For ``None''-abscissas ``func'' must be Table or Mapping!" ) 

    ## switch on func         
    if   callable   ( func ) :
        nfunc  = func 
        func   = ( nfunc ( x ) for x in abscissas )
    elif isinstance ( func , GT ) : pass                         ## generator 
    elif isinstance ( func , IT ) : pass                         ## iterable
    else :
        raise TypeError("Can't treat ``func''=%s"  %  func )
    ##
    from ostap.math.base import doubles
    ##
    if _A : return Ostap.Math.Interpolation.Table (   abscissas   , doubles ( func ) )
    return Ostap.Math.Interpolation.Table ( doubles ( abscissas ) , doubles ( func ) )


# =============================================================================
## Bernstein & Bspline interpolation
# =============================================================================
from ostap.math.bernstein import interpolate as interpolate_bernstein
from ostap.math.bspline   import interpolate as interpolate_bspline  
# =============================================================================

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
