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
import  ROOT, math, sys 
from    builtins          import range
from    ostap.core.core   import cpp, Ostap
from    ostap.core.ostap_types  import is_integer
from    ostap.math.base   import iszero, isequal, doubles 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.interpolation' )
else                       : logger = getLogger ( __name__                   )
# =============================================================================
if (3,3) <= sys.version_info : from collections.abc  import Iterable, Mapping
else                         : from collections      import Iterable, Mapping
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
    n2 = nmax//2
    s1 = ', '.join( ( '%.3g' % x  for x in self.x() [    : n2 ] ) ) 
    s2 = ', '.join( ( '%.3g' % x  for x in self.x() [ -1 :    ] ) ) 
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
Ostap.Math.Interpolation.Abscissas.items        = _a_iteritems_ 
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
    _n = self.size ()
    if _n <= nmax :
        s = ', '.join ( ( "%.3g: %.3g" %  (x,y)  for x,y in self ) ) 
        return 'Table({%s})' % s 
    ##
    n2 = nmax/2

    s1 = ', '.join ( ( '%.3g: %.3g' % self[i] for i in  range ( n2 ) ) ) 
    s2 =               '%.3g: %.3g' % self[ self.n() -1 ]
    
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
    N = len ( self ) 
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
    N = len ( self  )
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
        
        if 0 <= i < len ( self ) :
            return self.x ( i ) , self.y ( i )
        
        raise IndexError ('Invalid key %s' % i )
    
    elif isinstance ( i , slice ) :
        
        start , stop , step = i.indices ( self.n() )
        if   1 == step : return self.slice ( start , stop )

    raise TypeError ('Invalid key type/value %s' % i )

## # =============================================================================
## ## remove point from the point collection
## #  @code
## #  a = ...
## #  del a[2]
## #  @endcode 
## def _p_delitem_ (   self , i ) : 
##     """Remove point from abscissas
##     >>> a = ...
##     >>> del a[2]
##     """
##     if   isinstance ( i , int   ) :
##         if 0 <= i < self.n () : return self.remove ( i )
##         raise IndexError ('Invalid key %s' % i )
    
##     raise TypeError ('Invalid key/value %s' % i )

Ostap.Math.Interpolation.Table.__str__      = _p_str_
Ostap.Math.Interpolation.Table.__repr__     = _p_str_
Ostap.Math.Interpolation.Table.__iter__     = _p_iter_
Ostap.Math.Interpolation.Table.__len__      = lambda s : s.size()
Ostap.Math.Interpolation.Table.iteritems    = _p_iteritems_ 
Ostap.Math.Interpolation.Table.items        = _p_iteritems_ 
Ostap.Math.Interpolation.Table.__getitem__  = _p_getitem_  
## Ostap.Math.Interpolation.Table.__delitem__  = _p_delitem_  

# ==================================================================================
## dump interpoaltion table in a form of readabel table
#  @code
#  table = ...
#  print ( table.table (prefix = '# ' )
#  @endcode 
def _p_table_ ( t , title = '' , prefix = '' , fmt = '%+.6g' ) :
    """Dump interpoaltion table in a form of readabel table
    >>> table = ...
    >>> print ( table.table ( prefix = '# ' )
    """
    n = len ( t )
    rows = [ ( '   X' , '   Y' )  ]
    for i in range ( n ) :
        row = fmt % t.x( i ) , fmt % t.y ( i )
        rows.append ( row )
            
    import ostap.logger.table as T 
    if not title : title = 'Interpolation table %s ' % t.__class__.__name__ 
    return T.table ( rows , title = title , prefix = prefix , alignment = 'll' )


Ostap.Math.Interpolation.Table.table = _p_table_ 

# ==================================================================================
## Neville, Lagrange, & Berrut  interpolants
# ==================================================================================

# ==================================================================================
## printout for interpolants 
#  @code
#  print a 
#  @endcode 
def print_interpolant ( self , name , nmax = 7 ) :
    """Printout for Neville interpolant
    >>> print a  
    """
    n =  len ( self ) 
    if n <= nmax :
        s = ', '.join ( ( "%.3g: %.3g" %  (x,y)  for x,y in self ) ) 
        return '%s({%s})' % ( name , s )  
    ##
    n2 = min ( 2 , n ) 
    s1 = ', '.join ( ( '%.3g: %.3g' % self[i] for i in  range ( n2 ) ) ) 
    s2 =               '%.3g: %.3g' % self[ n - 1 ]
    
    return '%s(n=%d,{%s, ... , %s})' % ( name , n , s1 , s2 )

Ostap.Math.Neville         .__str__   = lambda s : print_interpolant ( s , 'Neville'        , 7 )
Ostap.Math.Lagrange        .__str__   = lambda s : print_interpolant ( s , 'Lagrange'       , 7 )
Ostap.Math.Newton          .__str__   = lambda s : print_interpolant ( s , 'Newton'         , 7 )
Ostap.Math.Berrut1st       .__str__   = lambda s : print_interpolant ( s , 'Berrut1st'      , 7 )
Ostap.Math.Berrut2nd       .__str__   = lambda s : print_interpolant ( s , 'Berrut2nd'      , 7 )
Ostap.Math.FloaterHormann  .__str__   = lambda s : print_interpolant ( s , 'FloaterHormann' , 7 )
Ostap.Math.Barycentric     .__str__   = lambda s : print_interpolant ( s , 'Barycentric'    , 7 )


# ==================================================================================
## Barycentric Lagrange interpolant 
# ==================================================================================


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
    IT = Iterable
    MT = Mapping

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
    if _A : return Ostap.Math.Barycentric ( abscissas , doubles ( func ) )
    ##
    return Ostap.Math.Barycentric ( doubles ( abscissas ) , doubles ( func ) )

# =============================================================================
## make new name 
interpolate = lagrange
barycentric = lagrange 

# ==================================================================================
## Construct collection of interpolated points/interpolation table  
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
    IT = Iterable
    MT = Mapping

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
    ##
    return Ostap.Math.Interpolation.Table ( doubles ( abscissas ) , doubles ( func ) )

# =============================================================================
## Bernstein & BSpline interpolation
# =============================================================================
from ostap.math.bernstein import interpolate as interpolate_bernstein
from ostap.math.bspline   import interpolate as interpolate_bspline  
# =============================================================================

# =============================================================================
## Scipy-based spline interpolation
# =============================================================================
try :
    from ostap.math.sp_interpolation import SplineInterpolator
    __all__ =  __all__ + ( 'SplineInterpolator', ) 
except  ImportError :
    pass 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
