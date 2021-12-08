#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/interpolation.py
#  Module with some useful utilities for dealing with interpolation.
#
#  It providesl follwoing interpolations (C++ fast)
#  - Lagrange  polynomial interpolation                                 (relativelu slow)
#  - Neville   polynomial interpolartion aka Nevile-Aitken interpolaton (slow)
#  - Newton    polynimial interpolation                                 (fastest)
#  - True Barycentric polymomial interpolation                          (relatively fast)
#  - Berrut's 1st rational interpoaltion                                (fast)
#  - Berrut's 2nd rational interpoaltion                                (fast)
#  - Floater-Hormann rational interpolation                             (fast)
#  - Bernstein polynomial inetrpoaltion                                  
#
#  @see Ostap::Math::Lagrange 
#  @see Ostap::Math::Neville 
#  @see Ostap::Math::Newton  
#  @see Ostap::Math::Berycentric 
#  @see Ostap::Math::Berrut1st 
#  @see Ostap::Math::Berrut2nd 
#  @see Ostap::Math::FloaterHormann
#
#  Features 
#  - Lagrange scheme allows to calcualet also the derivative with respect to \f$y_i\f$
#  - Neville scheme allows to calculate also the derivative with respect to \f$x\f$
#  - Newton, Berrut 1st, Berrut's 2nd, true-Barycentric and Floater-Hormann are very fast,
#    but they require some tiem for initialization adn thie time can be large.
#  - All polymonial interpoaltion behaves badly for largde degrees and unform/random grids
#  - Usage of dedicated Chebyshev and/or Lobatto grids partly solved this issue, butonly partly
#  - For large nmber of pointes rational interpolants behaves more stable, particularly
#    Floater-Hormann wwith relatively small value of parameter \f$d\f$ 
#
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
#  For  completeness see also:
#  - interpolation with Bersntein polynomials using on Newton-Bernstein algorithm
#  - interpolation with B-splines
#
#  In addition purely python interpolators are provided 
#  - Berrut1st
#  - Berrut2nd
#  - Barycentric
#  - FloaterHormann
# 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2018-07-22
# =============================================================================
"""Useful utilities for dealing with interpolation.

  It providesl follwoing interpolations (C++ fast)
  - Lagrange  polynomial interpolation                                 (relativelu slow)
  - Neville   polynomial interpolartion aka Nevile-Aitken interpolaton (slow)
  - Newton    polynimial interpolation                                 (fastest)
  - True Barycentric polymomial interpolation                          (relatively fast)
  - Berrut's 1st rational interpoaltion                                (fast)
  - Berrut's 2nd rational interpoaltion                                (fast)
  - Floater-Hormann rational interpolation                             (fast)
  - Bernstein polynomial inetrpoaltion                                  
  
  - see Ostap.Math.Lagrange 
  - see Ostap.Math.Neville 
  - see Ostap.Math.Newton  
  - see Ostap.Math.Berycentric 
  - see Ostap.Math.Berrut1st 
  - see Ostap.Math.Berrut2nd 
  - see Ostap.Math.FloaterHormann

  Features 
  - Lagrange scheme allows to calcualet also the derivative with respect to \f$y_i\f$
  - Neville scheme allows to calculate also the derivative with respect to \f$x\f$
  - Newton, Berrut 1st, Berrut's 2nd, true-Barycentric and Floater-Hormann are very fast,
  but they require some tiem for initialization adn thie time can be large.
  - All polymonial interpoaltion behaves badly for largde degrees and unform/random grids
  - Usage of dedicated Chebyshev and/or Lobatto grids partly solved this issue, butonly partly
  - For large nmber of pointes rational interpolants behaves more stable, particularly
  Floater-Hormann wwith relatively small value of parameter \f$d\f$ 
  

  - see Jean-Paul Berrut and Lloyd N. Trefethen, 
  Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
  
  - see https://doi.org/10.1137/S0036144502417715
  - see https://en.wikipedia.org/wiki/Lagrange_polynomial
  - see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
  
  Straightforward Lagrange and Neville algorithms are also provided:
  - both are rather slow: O(n^2) flops for evaluation,  while Neville is a bit faster 
  - Lagrange algorithm is not stable numerically, while Neville algorithm is more stable
  
  For completeness see also:
  - interpolation with bersntein polynomials using on Newton-Bernstein algorithm
  - interpolation with B-splines

  In addition purely python interpolators are provided
  - Berrut1st
  - Berrut2nd
  - Barycentric
  - FloaterHormann
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
    # interpoaltion abscissas
    'uniform_abscissas'     , ## generator for uniform absicssas (non-optimal!)
    'chebyshev_abscissas'   , ## generator for Chebyshev absicssas 
    'lobatto_abscissas'     , ## generator for Lobatto absicssas
    # python interpolators
    'Berrut1st'             , ## rational Berrut's 1st interpolant 
    'Berrut2nd'             , ## rational Berrut's 2nd interpolant 
    'Barycentric'           , ## polynomian true Barycentric interpolant 
    'FloaterHormann'        , ## rational Floater-Hormann interpolant 
    )
# =============================================================================
import  ROOT, math, sys, abc  
from    array                  import array
from    builtins               import range
from    ostap.core.core        import cpp, Ostap
from    ostap.core.ostap_types import ( is_integer, sequence_types,
                                        integer_types , dictlike_types )  
from    ostap.math.base        import iszero, isequal, doubles 
from    ostap.utils.utils      import vrange
import  ostap.math.reduce      
# =============================================================================
from   ostap.logger.logger     import getLogger
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
    n = self.n()
    a = self.atype() 
    if n <= nmax or 0 <= a :
        if   Ostap.Math.Interpolation.Abscissas.Uniform    == a : 
            return 'Abscissas(%d,%+.4g,%+.4g,%s)' % ( n , self.xmin () , self.xmax() , 'Uniform'    )
        elif Ostap.Math.Interpolation.Abscissas.Chebyshev  == a : 
            return 'Abscissas(%d,%+.4g,%+.4g,%s)' % ( n , self.xmin () , self.xmax() , 'Chebyshev'  )
        elif Ostap.Math.Interpolation.Abscissas.Chebyshev2 == a : 
            return 'Abscissas(%d,%+.4g,%+.4g,%s)' % ( n , self.xmin () , self.xmax() , 'Chebyshev2' )
        else :
            return 'Abscissas(%d,%s)'    % ( n , self.x () ) 
            
    ##
    n2 = max ( 1 , nmax//4 ) 
    s1 = ', '.join( ( '%.3g' % x  for x in self.x() [    : n2 ] ) ) 
    s2 = ', '.join( ( '%.3g' % x  for x in self.x() [ -1 :    ] ) )
    
    return 'Abscissas(n=%d,[%s, ... , %s])'    % ( n , s1 , s2 )
        
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
## create abscissas
def _a_new_init_ ( a , arg1 , *args ) :
    """ create abscissas
    """
    
    if isinstance ( arg1 , Ostap.Math.Interpolation.Abscissas ) :
        return a._old_init_ (           arg1   ) 

    if isinstance ( arg1 , sequence_types ) :
        if not args : 
            return a._old_init_ ( doubles ( arg1 ) )
        elif 1 == len ( args ) and isinstance ( args[0] , bool ) :
            return a._old_init_ ( doubles ( arg1 ) , args[0] )
            
    if isinstance ( arg1 , integer_types ) and 0 <= arg1 :
        if   2 == len ( args ) :
            return a._old_init_ ( arg1 , args[0] , args[1] )
        elif 3 == len ( args ) and isinstance ( args[2] , integer_types ) and 0 <= args[2] :
            return a._old_init_ ( arg1 , *args  )

    return a._old_init_ ( doubles ( arg1 , *args ) )
        
# =============================================================================
Abscissas = Ostap.Math.Interpolation.Abscissas 
if not hasattr ( Abscissas , '_old_init_' ) :
    Abscissas . _old_init_ = Abscissas . __init__
    ## modified contructor 
    def __aa_new_init__ ( a , arg1 , *args ) :
        """Modified constructor for interpolaiton abscissas 
        """
        _a_new_init_ ( a , arg1 , *args )
    __aa_new_init__ .__doc__ += '\n' + _a_new_init_ .__doc__
    __aa_new_init__ .__doc__ += '\n' + Abscissas . _old_init_ . __doc__
    Abscissas . __init__  = __aa_new_init__
    
       
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
    n = self.size ()
    if n <= nmax :
        s = ', '.join ( ( "%.3g: %.3g" %  (x,y)  for x,y in self ) ) 
        return 'Table({%s})' % s 
    ##
    n2 = nmax // 3 

    s1 = ', '.join ( ( '%.3g: %.3g' % self[i] for i in  range ( n2 ) ) ) 
    s2 =               '%.3g: %.3g' % self[ n - 1 ]
    
    return 'Table(n=%d,{%s, ... , %s})' % ( n , s1 , s2 )
                     
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


# =============================================================================
## convert interpolation table to Graph
def _p_graph_ ( self ) :
    """Convert interpolation table to Graph
    """
    N  = len ( self ) 
    gr = ROOT.TGraph ( N )
    for i in range ( N ) :
        gr.SetPoint ( i , self.x(i) , self.y(i) )
    return gr 
    
    

Ostap.Math.Interpolation.Table.__str__      = _p_str_
Ostap.Math.Interpolation.Table.__repr__     = _p_str_
Ostap.Math.Interpolation.Table.__iter__     = _p_iter_
Ostap.Math.Interpolation.Table.__len__      = lambda s : s.size()
Ostap.Math.Interpolation.Table.iteritems    = _p_iteritems_ 
Ostap.Math.Interpolation.Table.items        = _p_iteritems_ 
Ostap.Math.Interpolation.Table.__getitem__  = _p_getitem_  
## Ostap.Math.Interpolation.Table.__delitem__  = _p_delitem_  

Ostap.Math.Interpolation.Table.graph        = _p_graph_


# ==================================================================================
## Uniform interpolation abscissas
#  @code
#  for a in uniform_abscissas ( 0, 1 , 10 ) :
#  ... print ( a ) 
#  @endcode 
def uniform_abscissas ( low , high , N ) :
    """Uniform interpoaltion abscissas
    for a in uniform_abscissas ( 0, 1 , 10 ) :
    ... print ( a ) 
    """
    return vrange ( low , high , max ( N , 2 ) - 1 )

# ==================================================================================
## Chebyshev interpolation abscissas - optimal for polynoimial interpolation 
#  @code
#  for a in chebyshev_abscissas ( 0, 1 , 10 ) :
#  ... print ( a ) 
#  @endcode 
def chebyshev_abscissas ( low , high , N ) : 
    """Chebyshev interpolation abscissas - optimal for polynoimial interpolation 
    >>> for a in chebyshev_abscissas ( 0, 1 , 10 ) :
    >>> ... print ( a ) 
    """
    N = max ( N , 1 ) 
    for i in range ( N ) :
        a = 1.0 * ( 2 * ( N - i ) - 1 ) * math.pi / ( 2 * N ) 
        x = math.cos ( a )
        yield 0.5 * (  ( 1 - x ) * low + ( 1 + x ) * high )


# ==================================================================================
## Lobatto interpolation abscissas - optimal for polynoimial interpolation 
#  @code
#  for a in lobatto_abscissas ( 0, 1 , 10 ) :
#  ... print ( a ) 
#  @endcode 
def lobatto_abscissas ( low , high , N ) : 
    """Lobatto interpolation abscissas - optimal for polynoimial interpolation 
    >>> for a in Lobatto_abscissas ( 0, 1 , 10 ) :
    >>> ... print ( a ) 
    """
    N = max ( N , 2 ) 
    for i in range ( N ) :
        a = 1.0 * ( N - i - 1 ) * math.pi / ( N - 1  ) 
        x = math.cos ( a )
        yield 0.5 * (  ( 1 - x ) * low + ( 1 + x ) * high )
        
# ==================================================================================
## Print interpolation table as table
#  @code
#  table= ...
#  print ( table.table() )
#  @endcode 
def _tab_print_ ( t , title = '' , prefix = '' , alignment = 'll' , xfmt = '%+.5g' , yfmt = '%+-.5g' ) :
    """Print interpolation table as table
    >>> table= ...
    >>> print ( table.table() )
    """
    rows = [ ('Abscissa' , 'Value' ) ] 
    for i in range ( t.size() ) :
        x = t.x ( i )
        y = t.y ( i )
        row = xfmt %  x, yfmt % y
        rows.append ( row )
        
    if not title : title = 'Interpolation Table' 
    import ostap.logger.table as T
    return T.table ( rows , title = title , prefix = prefix , alignment = alignment )

Ostap.Math.Interpolation.Table.table = _tab_print_ 

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
    n = len ( self )
    a = self.atype()
    if   Ostap.Math.Interpolation.Abscissas.Uniform    == a : 
        return '%s(%s[%d,%+.4g,%+.4g])' % ( name , 'Uniform'   , n , self.xmin () , self.xmax() )
    elif Ostap.Math.Interpolation.Abscissas.Chebyshev  == a : 
        return '%s(%s[%d,%+.4g,%+.4g])' % ( name , 'Chebyshev' , n , self.xmin () , self.xmax() )
    elif Ostap.Math.Interpolation.Abscissas.Chebyshev2 == a : 
        return '%s(%s[%d,%+.4g,%+.4g])' % ( name , 'Lobatto'   , n , self.xmin () , self.xmax() )
        
    if n <= nmax :
        s = ', '.join ( ( "%.3g: %.3g" %  (x,y)  for x,y in self ) ) 
        return '%s({%s})' % ( name , s )  
    ##
    n2 = min ( 2 , n ) 
    s1 = ', '.join ( ( '%.3g: %.3g' % self[i] for i in  range ( n2 ) ) ) 
    s2 =               '%.3g: %.3g' % self[ n - 1 ]
    
    return '%s(n=%d,{%s, ... , %s})' % ( name , n , s1 , s2 )

Ostap.Math.Neville         .__str__   = lambda s : print_interpolant ( s , 'Neville'                       , 7 )
Ostap.Math.Lagrange        .__str__   = lambda s : print_interpolant ( s , 'Lagrange'                      , 7 )
Ostap.Math.Newton          .__str__   = lambda s : print_interpolant ( s , 'Newton'                        , 7 )
Ostap.Math.Berrut1st       .__str__   = lambda s : print_interpolant ( s , 'Berrut1st'                     , 7 )
Ostap.Math.Berrut2nd       .__str__   = lambda s : print_interpolant ( s , 'Berrut2nd'                , 7 )
Ostap.Math.FloaterHormann  .__str__   = lambda s : print_interpolant ( s , 'FloaterHormann%d' % s.d() , 7 )
Ostap.Math.Barycentric     .__str__   = lambda s : print_interpolant ( s , 'Barycentric'                   , 7 )

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
## get all weigths from Berrut1st and Berrut2nd interpolants
#  @code
#  interpolant = ....
#  interpolant.weights()
#  @endcode 
def _b12_weights_ ( self ) :
    """Get all weigths from Berrut1st and Berrut2nd interpolants
    >>> interpolant = ....
    >>> interpolant.weights()
    """
    N = len ( self ) 
    return array ( 'd' , ( self.weight ( i ) for i in range ( N ) ) )  

for t in ( Ostap.Math.Berrut1st   , 
           Ostap.Math.Berrut2nd   ) :
    
    t.weights = _b12_weights_ 

# =============================================================================
## Get sum of weights for barycentric interpolation is muts be zero
#  @code
#  interpolant = ...
#  interpolant.sumw() 
#  @encode
def _bi_sumw_ ( self ) :
    """Get sum of weights for barycentric interpolation is muts be zero
    >>> interpolant = ...
    >>> interpolant.sumw() 
    """
    N = len ( self )
    if  0 == N : return 0 
    return sum ( self.weights ()  ) 

# =============================================================================
## Does barycentric interpolant has poles ?
#  @code
#  interpolant = ...
#  interpolant.poles() 
#  @encode
def _bi_poles_ ( self ) :
    """Does barycentric interpolant has poles ?
    >>> interpolant = ...
    >>> interpolant.poles() 
    """
    N = len ( self )
    for i in range ( N -1 ) :
        if 0 < self.weight ( i ) * self.weight ( i + 1 ) : return True 
    return False 

for t in ( Ostap.Math.Berrut1st      , 
           Ostap.Math.Berrut2nd      ,
           Ostap.Math.Barycentric    , 
           Ostap.Math.FloaterHormann ) : 
    
    t.sumw  = _bi_sumw_ 
    t.poles = _bi_poles_ 


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
from operator import mul as op_mul
from operator import add as op_add 

# =============================================================================
## Base class for barycentric-like interpolation 
class BaseInterpolant(object) :
    """Base class for barycentric-like interpolation 
    """
    
    __metaclass__ = abc.ABCMeta
  
    def __init__ ( self             ,
                   data             ,
                   scaler  = op_mul , 
                   adder   = op_add ) :

        self.__table = [] 
        if   isinstance ( data , dictlike_types ) :
            for x , v in data.items() :
                self.__table.append ( ( x , v ) )
        elif isinstance ( data , sequence_types ) :
            for x , v in data :
                self.__table.append ( ( x , v ) )
        else :            
            raise TypeError('Invalid type %s' % type ( data ) )  
        
        self.__table.sort ( key = lambda i : i[0] )
        assert self.table , 'Interpolation table must not be empty!'
        
        self.__adder  = adder
        self.__scaler = scaler 

        self.__xmin   = self.__table[ 0][0]
        self.__xmax   = self.__table[-1][0]
        
    ## make the actual interpolation 
    def __call__ ( self , x ) :
        """Make the actual interpolation""" 
    
        s1 = None
        s2 = 0

        x  = float ( x )
        
        for i, item in enumerate ( self.__table  ) :
            
            xi , yi = item
            
            if x == xi or isequal ( x , xi ) : return yi  ## RETURN
            
            ## calculate weight 
            wi = self.weight ( i ) * 1.0 / ( x - xi ) 

            if 0 == i :
                s1 = self.__scaler ( yi , wi )
            else :
                s1 = self.__adder  ( s1 , self.__scaler ( yi , wi ) )
                
            s2 += wi 

        return self.__scaler ( s1 , 1.0/s2 )

    @property
    def scaler ( self ) :
        """``scaler'' : operation ( obj , weight ) -> obj"""
        return self.scaler

    @property
    def adder ( self ) :
        """``adder'' : operation  ( obj , obj ) -> obj"""
        return self.adder 
        
    # =========================================================================
    @abc.abstractmethod
    def weight ( self , index ) :
        """Get the weigth for the given interpolation node"""
        return NotImplemented
         
    # =========================================================================
    @property
    def table ( self ) :
        """``table'' : actual interpolation table"""
        return self.__table

    # =========================================================================
    ## the length of the interpolation table
    #  - number of interpolation points 
    def __len__ ( self ) :
        """the length of the interpolation table
        - number of interpolation points
        """
        return len ( self.__table )

    # ==========================================================================
    ## array of weigths
    def weights ( self ) :
        """Get array of weights
        """
        N = len ( self ) 
        return array ( 'd' , ( self.weight ( i ) for i in range ( N ) ) )  

    # ==========================================================================
    ## sum of all weights
    #  - it must be zero for barycentric weights 
    def sumw ( self )  :
        """sum of all weights
        _ it must be zero for barycentric weights 
        """
        N = len ( self )
        if  0 == N : return 0 
        g = ( self.weight ( i ) for i in range ( N ) ) 
        return sum ( g )

    # =========================================================================
    ## Does barycentric interpolant has poles ?
    #  @code
    #  interpolant = ...
    #  interpolant.poles() 
    #  @encode
    def poles ( self ) :
        """Does barycentric interpolant has poles ?
        >>> interpolant = ...
        >>> interpolant.poles() 
        """
        N = len ( self )
        for i in range ( N -1 ) :
            if 0 < self.weight ( i ) * self.weight ( i + 1 ) : return True 
        return False 
    
    ## get the minimal value in the interpolaiton table 
    def xmin ( self ) : return self.__xmin
    ## get the maximal minimal value in the interpolaiton table 
    def xmax ( self ) : return self.__xmax

    def __str__  ( self ) :
        return "%s(%d,%.3g,%.3g)" % ( self.__class__.__name__ ,
                                      len ( self ) ,
                                      self.xmin () ,
                                      self.xmax () )
    
# =============================================================================
## Berrut's 1st barycentric rational interpolant
class Berrut1st(BaseInterpolant) :
    """Berrut's 1st barycentric rational interpolant
    """ 
    def weight ( self , index ) :
        """Get the weigth for the given interpolation node"""
        return 1.0 if ( index % 2 ) else -1.0 
    
# =============================================================================
## Berrut's 2nd barycentric rational interpolant
class Berrut2nd(BaseInterpolant) :
    """Berrut's 2nd  barycentric rational interpolant
    """ 
    def weight ( self , index ) :
        """Get the weigth for the given interpolation node"""
        
        if index == 0 : return 1.0
        N = len ( self ) 
        if index + 1 == N  :
            return 1.0 if ( N % 2 ) else -1.0 

        return 2.0 if ( index  % 2 ) else -2.0 


# =============================================================================
## true Barycentric polymnomial interpolant
class Barycentric(BaseInterpolant) :
    """True barycentric polynomial interpolant
    """
    def __init__ ( self             ,
                   data             ,
                   scaler  = op_mul , 
                   adder   = op_add ) :
        
        super(Barycentric,self).__init__ ( data , scaler , adder )

        N = len ( self )

        ws = [] 
        for i, item in enumerate ( self.table ) :

            xi , yi = item

            ww = 1.0
            for j in range ( N ) :
                if i != j :
                    xj = self.table[j][0] 
                    ww *= ( xi - xj )
                    
            ws.append ( 1.0 / ww )
            
        self.__weigths = array ( 'd' , ws ) 
        
    def weight ( self , index ) :
        """Get the weigth for the given interpolation node"""

        return self.__weigths[index]
    
    @property
    def weights ( self ) :
        """Get list of weights"""
        return self.__weights

# =============================================================================
## FloaterHormann rational interpolant
class FloaterHormann(BaseInterpolant) :
    """FloaterHormann rational interpolant
    """
    def __init__ ( self             ,
                   data             ,
                   degree  = 3      , 
                   scaler  = op_mul , 
                   adder   = op_add ) :

        assert isinstance ( degree , integer_types ) and 0 <= degree , \
               "FloaterHormann: ``degree'' must be non-negative!"
        
        super(FloaterHormann,self).__init__ ( data , scaler , adder )

        N = len ( self )
        n = max ( 0  , N - 1 ) 

        self.__degree = min ( degree , n )
        
        d = self.__degree
            
        ws = []
        
        for i in range ( N  ) :

            xi = self.table[i][0]

            ib = 0.0
            jmin = max ( i - d , 0     )
            jmax = min ( i     , n - d ) 
            for j in range ( jmin , jmax + 1 ) :

                kmin = j
                kmax = j + d

                bb = 1.0
                for k in range ( kmin , kmax + 1 ) :

                    xk = self.table[k][0]
                    if k != i : bb /= abs ( xi - xk )
                ib += bb
                
            wi =  ib * ( 1 if i % 2 else -1 ) 
            ws.append ( wi )

        self.__weights = array( 'd',  ws ) 

    def weight ( self , index ) :
        """Get the weigth for the given interpolation node"""

        return self.__weights[index]
    
    @property
    def weights ( self ) :
        """Get list of weights"""
        return self.__weights

    @property
    def degree ( self ) :
        """``degree'' : degree for FloaterHormann rational interpolant"""
        return self.__degree

    def __str__  ( self ) :
        return "%s%s(%d,%.3g,%.3g)" % ( self.__class__.__name__ ,
                                        self.degree , 
                                        len ( self ) ,
                                        self.xmin () ,
                                        self.xmax () )
    

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
