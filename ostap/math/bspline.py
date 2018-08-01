#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file bspline.py
#
#  Module with some useful utilities for dealing with BSpline
#  @see Ostap::Math::BSpline
#
# - control_polygon      : get a control polygon for BSpline
# - upper_convex_hull    : upper convex hull for BSpline
# - lower_convex_hull    : lower convex hull for BSpline
# - convex_hull          :       convex hull for BSpline
# - crossing_points      : get crossing points of control polygon with x-axis
# - solve                : solve equation B(x) = C
# - interpolate          : construct interpolating B-spline 
# - approximate          : construct approximating B-spline
# - generate&shoot       : generate random numbers 
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-12-01
# =============================================================================
""" Module with some useful utilities for dealing with BSplines
- control_polygon      : get a control polygon for BSpline
- upper_convex_hull    : upper convex hull for BSpline
- lower_convex_hull    : lower convex hull for BSpline
- convex_hull          :       convex hull for BSpline
- crossing_points      : get crossing points of control polygon with x-axis
- solve                : solve equation B(x) = C
- interpolate          : construct interpolating B-spline 
- approximate          : construct approximating B-spline
- generate&shoot       : generate random numbers 
#"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = (
    'control_polygon'      , ## get a control polygon for BSpline
    'upper_convex_hull'    , ## upper convex hull for BSpline
    'lower_convex_hull'    , ## lower convex hull for BSpline
    'convex_hull'          , ##       convex hull for BSpline
    'crossing_points'      , ## get crossing points of control polygon with x-axis
    'solve'                , ## solve equation B(x) = C
    'interpolate'          , ## spline interpolation
    'approximate'          , ## variation diminishing approximation 
    'generate'             , ## generate random numbers 
    'shoot'                , ## generate random numbers 
    )
# =============================================================================
import  ROOT, math  
from    ostap.core.core      import cpp, Ostap, funID
from    ostap.math.base      import iszero, isequal, signum, doubles
import  ostap.math.bernstein 
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.bspline' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
## short name 
BSpline = cpp.Ostap.Math.BSpline
# =============================================================================
## get control polygon for BSpline
def control_polygon ( bs )  :
    """Get control polygon for BSpline
    >>>  bspline = ...
    >>>  cp = bspline.control_polygon ()
    >>>  cp = control_polygon( bspline )  ##  ditto 
    """
    return Ostap.Math.control_polygon (  bs.bspline () ) 
# =============================================================================
## get upper convex hull for  BSpline
def upper_convex_hull ( bs ) :
    """Get upper convex hull for  BSpline
    >>> bspline = ...
    >>> upper = bspline.upper_convex_hull  ()
    >>> upper = upper_convex_hull ( bspline ) ## ditto
    - for all values of x: bspline(x) <=  upper ( x ) 
    """
    return Ostap.Math.upper_convex_hull ( bs.bspline()  )
# =============================================================================
## get lower convex hull for  BSpline
def lower_convex_hull ( bs ) :
    """Get lower convex hull for  BSpline
    >>> bspline = ...
    >>> lower = bspline.lower_convex_hull ()
    >>> lower = lower_convex_hull ( bspline )
    - for all values of x: bspline(x) >=  lower ( x ) 
    """
    return Ostap.Math.lower_convex_hull ( bs.bspline()  )
# =============================================================================
## get upper & lower convex hulls for  Bspline
def  convex_hull ( bs ) :
    """Get lower& upper convex hulls for  BSpline
    >>> bspline = ...
    >>> lower,upper = bspline.convex_hull  ()
    >>> lower,upper = convex_hull ( bspline ) ## ditto  
    - for all values of x: lower(x) <= bspline(x) <= upper ( x ) 
    """
    bs1 = bs.bspline()
    return bs1.lower_convex_hull () , bs1.upper_convex_hull ()
# =============================================================================
## get abscissas of crosssing point of the control polygon with x-axis 
#  @param  b bernstein polynomial
#  @return abscissas of crossing points of the control  polygon with x-axis 
def crossing_points  ( bp , all = False ) :
    """Get abscissas of crosssing point of the control polygon with x-axis 
    >>> bspline = ...
    >>> xps = bspline.crossing_points()
    >>> xps = crossing_points ( bspline) ## ditto 
    """
    cps = Ostap.Math.crossing_points ( bp.bspline() , all )
    return tuple ( cp for cp in cps ) 
# =============================================================================
## evaluate the spline given by the set of knots, control points and the order
#  using de Boor's algorithm
#  @see https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
#  @param x      (INPUT) the point
#  @param order  (INPUT) the order of spline 
#  @param knots  (INPUT) vector of knots
#  @param points (INPUT) vector of control points 
#  @return value of the spline in point x 
def deboor ( x      ,
             order  ,
             knots  ,
             points ) :
    """Evaluate the spline given by the set of knots, control points and the order
    using de Boor's algorithm
    - see https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
    """
    VD = cpp.std.vector('double')
    ## convert knots if needed 
    if not isinstance ( knots  , VD ) : knots   =  doubles ( knots  )
    ## convert knots if needed 
    if not isinstance ( points , VD ) : points  =  doubles ( points )
    ##
    return Ostap.Math.deboor (  x , order , knots , points ) 

# =============================================================================
## reconstruct knot vector from greville abscissas and spline degree 
def knots_from_abscissas ( abscissas , order , convert = True ) :
    """Reconstruc knot vector from greville abscissas and spline degree
    abscissas : vector of abscissas
    order     : the order/degree  of spline
    convert   : convert to tuple ?
    >>> data = [ 0,2,3,4,5,10,11]
    >>> knots = knots_from_abscissas ( data , 3 ) 
    """

    ## VD = cpp.std.vector('double')
    ## ## convert knots if needed 
    ## if not isinstance ( abscissas , VD ) : abscissas =  doubles ( abscissas )
    ## knots = Ostap.Math.knots_from_abscissas ( abscissas , order )
    ## if convert  : knots = tuple ( knots ) 
    ## return  knots

    if len(abscissas) < order  + 1 :
        raise AttributeError("Vector of abscissas is too short")
    degree = order 
    abscissas = list(abscissas)
    abscissas.sort()
    af = abscissas[ 0]
    al = abscissas[-1]

    knots = [ af ] + abscissas + [ al ]
    
    if 1 ==  order : return knots
    while len(knots) < len(abscissas) + order + 1 :
        knots = [ af ] + knots + [ al ]

    return knots 
    
    
    knots = (degree+1)*[ af ]
    N     = len ( abscissas )
    
    for i in  range(1,N) :
        st = sum ( knots [ -(degree-1) :] )
        ti = abscissas[i] * (degree  ) - st 
        knots.append ( ti )
        
    while len(knots) < N + order + 1 :
        knots.append ( al )

    abscissas = doubles ( abscissas )
    knots2 = Ostap.Math.knots_from_abscissas ( abscissas , order )
    
    return knots, knots2  
        
    
    

# =============================================================================    
## Construct the interpolation Bs-spline 
#  @code
#  b1 = interpolate ( lambda x : x*x , [0,0.5,1,2]    , 3 )  
#  b2 = interpolate ( { 0:0 , 0.5:0.25 , 1:1 } , None , 3 )  
#  b3 = interpolate ( [0,0.25,1,4] , [ 0,0.5, 1,2]    , 3 )  
#  @endcode 
#  @param func      (INPUT) the function or list of function values
#  @param abscissas (INPUT) absiccas
#  @param order     (INPUT) the order of spline
#  @return interpolation spline 
def interpolate ( func , abscissas , spline , *args ) :
    """Construct the interpolation B-spline
    
    func      : the ``function''
    abscissas : abscissas, if None/Empty,  Greville's abscissas from spline will be used
    spline    : th epsline qwill be constructed forem ``spline'' and ``args''
    
    :Example:
    
    >> b1 = interpolate ( lambda x : x*x , [0,0.5,1,2]    , 3 )  
    >> b2 = interpolate ( { 0:0 , 0.5:0.25 , 1:1 } , None , 1 )  
    >> b3 = interpolate ( [0,0.25,1,4] , [ 0,0.5, 1,2]    , 2 )  
    >> b4 = interpolate ( lambda x : x*x , None    , 3 )  
    """
    
    bs  = Ostap.Math.BSpline ( spline , *args )
    if not abscissas :
        abscissas = bs.greville_abscissas() 
        
    from types       import GeneratorType as GT
    from collections import Iterable      as IT
    from collections import Mapping       as MT
    
    if isinstance ( abscissas , GT ):
        abscissas = [ x for x in abscissas ]
        
    if   callable ( func ) and abscissas :
        func = [ func (x)  for x in abscissas ]                  ## callable 
    elif isinstance ( func , dict ) and not abscissas :          ## mapping 
        keys = func.keys()
        keys.sort()
        abscissas = [ x       for x in keys ]
        func      = [ func[x] for x in keys ]
    elif isinstance ( func , GT   ) : func = [ f for f in func ] ## generator
    elif isinstance ( func , IT   ) : pass                       ## iterable 
    else :
        raise TypeError("Can't treat ``func''=%s"  %  func )

    ##
    from ostap.math.base import doubles
    _x = doubles ( abscissas )
    _y = doubles ( func      )

    bs  = Ostap.Math.BSpline ( spline , *args )
    if len(bs) != len(_x) :
        raise TypeError("Can't interpolate:  different number of parameters  %s/%s" %  ( len(_x) ,  len(bs) ) )

    sc  = Ostap.Math.Interpolation.bspline ( _x , _y , bs )
    if sc.isFailure() :
        raise TypeError("Ostap.Math.Bspline: Can't iterpolate!%s" %  sc )
    
    return bs 

# =============================================================================    
## Construct the variation diminishing approximation 
#  @code
#  b1 = approximate ( lambda x : x*x , [0,0.5,1,2]    , 3 )  
#  b2 = interpolate ( { 0:0 , 0.5:0.25 , 1:1 } , None , 3 )  
#  b3 = interpolate ( [0,0.25,1,4] , [ 0,0.5, 1,2]    , 3 )  
#  @endcode 
#  @param func      (INPUT) the function or list of function values
#  @param abscissas (INPUT) absiccas
#  @param order     (INPUT) the order of spline
#  @return interpolation spline 
def approximate ( func , spline , *args ) : 
    """Construct the interpolation B-spline
    """
    bs = Ostap.Math.BSpline ( spline , *args )
    xv = bs.greville_abscissas()
    N  =  bs.npars()
    for i in  range(N) : bs.setPar ( i , func ( xv[i] ) )
    return bs
    
# =============================================================================
## merge duplicated roots together 
def _merge_ ( roots , dx ) :
    """Merge duplicated roots together 
    """
    merged = []
    for x in roots :
        found = False 
        n = len(merged)
        for i in range (n) :
            y = merged[i]
            if isequal ( x, y ) or isequal ( x - y + dx , dx ) :
                found = True
                merged[i] = 0.5 * ( x + y )
        if not found : merged.append ( x )

    return tuple(merged) 

from ostap.math.bernstein import _scale_

# =============================================================================
## solve equation \f$ B(x) = C \f$
#  @param bs (INPUT) b-spline
#  @param C  (INPUT) right-hand side
#  @return solutions 
def solve ( bs , C = 0 , split = 5 ) :
    """Solve equation B(x) = C
    
    bs : (input)  b-spline
    C  : (input) right-hand side
    
    :Example:  
    >>> bs = ...
    >>> print bs.solve() 
    """
    
    _bs = bs - C 
    
    if 1 >=  _bs.degree() :
        return  _bs.crossing_points ( True )

    ## scale if needed 
    _bs = _scale_ ( _bs ) 

    ## check zeroes of control polygon
    cps = _bs.crossing_points ( )
    ncp = len(cps)
    bsn = _bs.norm() 
    ##
    xmin = _bs.xmin()
    xmax = _bs.xmax()
    
    for i in  range ( 20 ) :
        
        if not cps  : return ()                         ## RETURN 
        
        inserted = False
        for x in cps :
            if  isequal   ( x , _bs.xmin() ) : continue ## CONTINUE 
            if  isequal   ( x , _bs.xmax() ) : continue ## CONTINUE
            bx = _bs ( x )
            if  iszero ( bx ) or isequal ( bx + bsn ,   bsn ) : continue
            if _bs.insert ( x ) : inserted = True 

        ## all zeroes are already inserted....
        if not inserted :
            return cps                                  ## RETURN
        
        cps = _bs.crossing_points ( False )
        if not  cps :
            ## roots have disapperead.. it could happen near multiple roots 
            cps = _bs.crossing_points ( True )

        ##  merge duplicates (if any) 
        cps = _merge_ ( cps , max ( abs ( xmin ) , abs( xmax  ) ) ) 
            
    ## if no convergency....
    ## if not cps or 1 == len(cps) or split <= 0 : return cps 
    if not cps or split <= 0 : return cps

    
    ## come  back to the original spline 
    _bs =  bs - C
    
    ## add the current/best estimates of roots 
    for x in cps : _bs.insert ( x )


    ## get internal roots only (no edges) 
    _cps = cps
    _bsn = _bs.norm()
    left_root  = False 
    right_root = False 
    if isequal ( _bs[ 0] + _bsn , _bsn ) :
        left_root  = True 
        _cps = _cps[ 1 :    ]
    if isequal ( _bs[-1] + _bsn , _bsn ) :
        right_root = True 
        _cps = _cps[   : -1 ]
    
    cpmin = _cps[ 0]
    cpmax = _cps[-1]
    
    dc = cpmax - cpmin 

    xmin = _bs.xmin()
    xmax = _bs.xmax()
    dx   = xmax - xmin

    ## 
    if 1 == len(cps) :
        
        x1 = 0.05 * xmin + 0.95 * cpmin
        x2 = 0.05 * xmax + 0.95 * cpmax
        
        b1 = _scale_ ( Ostap.Math.BSpline ( _bs , xmin , x1   ) ) 
        b2 = _scale_ ( Ostap.Math.BSpline ( _bs , x1   , x2   ) ) 
        b3 = _scale_ ( Ostap.Math.BSpline ( _bs , x2   , xmax ) ) 

        split -= 1 
        return _merge_ ( solve ( b1 , C = 0 , split = split ) +
                         solve ( b2 , C = 0 , split = split ) +
                         solve ( b3 , C = 0 , split = split ) , max ( abs ( xmin ) , abs( xmax  ) ) )

    
    dl = 0.05 * min ( cpmin -  xmin , dc )
    dr = 0.05 * min (  xmax - cpmax , dc )

    x_left  = cpmin - dl
    x_right = cpmax + dr 
    x_mid   = 0.5 * ( cpmin + cpmax ) 
    
    import random
    ntry = 50

    for i in range ( ntry ) : 
        if not isequal ( _bs ( x_left    ) + _bsn , _bsn ) : break 
        alpha   = random.uniform ( 0.05 , 1 )
        x_left  = max ( xmin , cpmin - alpha * dl )
        
    for i in range ( ntry ) :
        if not isequal ( _bs ( x_mid     ) + _bsn , _bsn ) : break 
        alpha   = random.uniform ( -1 , 1 )
        x_mid   = 0.5 * ( cpmin + cpmax) + alpha * dc / 20 
        
    for i in range ( ntry ) :
        if not isequal ( _bs ( x_right ) + _bsn , _bsn ) : break 
        alpha   = random.uniform ( 0.05 , 1 ) 
        x_right = min ( xmax , cpmax + alpha * dr )
        
    ## split!
    split -= 1

    if dc > 0.1 * dx : ##  split 
        
        bs1 = _scale_ ( Ostap.Math.BSpline ( _bs , x_left , x_mid   ) )
        bs2 = _scale_ ( Ostap.Math.BSpline ( _bs , x_mid  , x_right ) ) 
        roots = _merge_  ( solve ( bs1 , C = 0 , split = split ) +
                           solve ( bs2 , C = 0 , split = split ) , max ( abs ( xmin ) , abs( xmax  ) ) )
        
    else :             ##  trim 
        
        bst   = _scale_ ( Ostap.Math.BSpline ( _bs , x_left , x_right ) ) 
        roots = _merge_ ( solve ( bst , C = 0 , split = split )  , max ( abs ( xmin ) , abs( xmax  ) ) )
    
    if left_root  : roots = (xmin,) + roots 
    if right_root : roots =           roots + (xmax,)
    
    return roots


# =============================================================================
## generate random numbers from b-spline-distribuitions
#  @code
#  >>> func = ...
#  >>> for x in func.generate( 1000 ) : print x 
#  @endcode
def generate ( fun , num ) :
    """Generate random numbers from bspline-like distribuitions
    >>> func = ...
    >>> for x in func.generate( 1000 ) : print x 
    """
    bs  = fun.bspline() 
    xmn = bs.xmin ()
    xmx = bs.xmax ()
    ymx = max ( bs.pars() )
    i   = 0 
    from random import uniform as _uniform_
    while i < num : 
        x = _uniform_ ( xmn , xmx ) 
        y = _uniform_ (   0 , ymx )
        v = fun ( x )
        if v >= y :
            i+= 1 
            yield x
            
# =============================================================================
## generate random numbers from b-spline-distribuitions
#  @code
#  >>> func = ...
#  >>> for x in func.generate( 1000 ) : print x 
#  @endcode
def _random_generate_bspline_ ( fun , num ) :
    """Generate random numbers from bspline-like distribuitions
    >>> func = ...
    >>> for x in func.generate( 1000 ) : print x 
    """
    xmn = fun.xmin ()
    xmx = fun.xmax ()
    ymx = max ( fun.bspline().pars() )
    i   = 0 
    while i < num : 
        x = _uniform_ ( xmn , xmx ) 
        y = _uniform_ (   0 , ymx )
        v = fun ( x )
        if v >= y :
            i+= 1 
            yield x
            
# =============================================================================
## Get random number from bspline-like distribuitions
#  @code
#  >>> func = ...
#  >>> print fun.shoot() 
#  @endcode
def shoot ( fun ) :
    """Get random number from bspline-like distribuitions
    >>> func = ...
    >>> print func.shoot()  
    """
    for x in generate(  fun , 1 ) :
        return x

# =============================================================================
## (Redefine standard constructor to allow usage of python lists&tuples)
#  Lists and tuples are converted on flight to :
# - std::vector<double> 
# - std::vector<std::complex<double> >
def _new_init_ ( t ,  *args )  :
    """(Redefine standard constructor to allow usage of python lists&tuples)
    Lists and tuples are  converted on flight to :
    - std::vector<double> 
    """
    from ostap.math.base import doubles, complexes
    
    largs = list (  args )
    alen  = len  ( largs )
    
    for i in range(alen) :
        
        arg = largs[i] 
        if not isinstance ( arg ,  ( list , tuple ) ) : continue 
        
        try: 
            _arg = doubles  ( *arg  )
            largs[i] = _arg
            continue 
        except TypeError : pass
        
    targs = tuple(largs)
    ## use old constructor 
    t._old_init_ ( *targs ) 

# =============================================================================
## set, get & iterator
from ostap.math.bernstein import _p_set_par_ , _p_get_par_, _p_iter_ 
   
# =============================================================================    
for  p in ( Ostap.Math.BSpline          ,
            Ostap.Math.PositiveSpline   ,
            Ostap.Math.ConvexOnlySpline ,
            Ostap.Math.MonotonicSpline ,
            Ostap.Math.ConvexSpline     ) :
    
    p.lower_convex_hull =  lower_convex_hull
    p.upper_convex_hull =  upper_convex_hull
    p.convex_hull       =        convex_hull
    p.convex_hulls      =        convex_hull
    p.control_polygon   =    control_polygon
    p.crossing_points   =    crossing_points
    p.solve             =              solve
    p.generate          =           generate
    p.shoot             =              shoot
    
    p.__setitem__  = _p_set_par_
    p.__getitem__  = _p_get_par_
    p.__len__      = lambda s     : s.npars() 
    p.__iter__     = _p_iter_
    p.__contains__ = lambda s , i : 0<=i<len(s)

    if not hasattr ( p , '_old_init_' ) :
        p._old_init_ = p.__init__
        ## Modifed constructor to allow python lists/tuples
        def _p_new_init_ ( s ,  *args ) :
            """Modifed constructor to allow python lists/tuples
            """
            _new_init_ ( s , *args )
            
        _p_new_init_.__doc__ += '\n' +   _new_init_.__doc__ 
        _p_new_init_.__doc__ += '\n' + p._old_init_.__doc__ 
        p.__init__ = _p_new_init_ 



for  p in ( Ostap.Math.BSpline2D           ,
            Ostap.Math.BSpline2DSym        ,
            Ostap.Math.PositiveSpline2D    ,
            Ostap.Math.PositiveSpline2DSym ) :

    p.__setitem__  = _p_set_par_
    p.__getitem__  = _p_get_par_
    p.__len__      = lambda s     : s.npars() 
    p.__iter__     = _p_iter_
    p.__contains__ = lambda s , i : 0<=i<len(s)

# =============================================================================
## set parameter for polynomial/spline functions
#  @code
#  fun = ...
#  fun[1] = 10.0
#  @endcode 
def _p2_set_par_ ( o , index , value ) :
    """Set parameter for polynomial/spline function
    >>> fun = ...
    >>> fun[1]   = 10.0
    >>> fun[1,2] = 15.0
    """
    if isinstance ( index , ( int , long ) ) :                  
        n = o.npars() 
        if   index <  0 :  index += n
        if not 0 <= index < n :
            raise IndexError('[%s] index out of range [0,%d)' % ( index , n ) ) 
        return o.setPar ( index , value )

    try :        
        ix , iy = index
    except :
        raise IndexError('Invalid index %s/%s' % ( index  , type(index) ) )
    
    if o.index ( ix , iy ) not in o :
        raise IndexError('Invalid index (%s,%s)' % ( ix , iy ) )
    return o.setPar ( ix , iy , value )
    

# =============================================================================
## get parameter from polynomial/spline functions
#  @code
#  fun = ...
#  print fun[1], fun[-1], fun[3,4] 
#  @endcode 
#  Slice  semantic is also supported:
#  @code
#  fun = ...
#  print fun[:3] , fun[4:], fun[2:8] ,  fun [2::2] 
#  @endcode 
def _p2_get_par_ ( o , index ) :
    """Get parameter from polynomial/spline function
    >>> fun = ...
    >>> print fun[ 1], fun[2,4]
    >>> print fun[-1]
    Slice  semantic is also supported:
    >>> print fun[:3]
    >>> print fun[2::3]
    """
    if isinstance ( index , ( int , long ) ) :                      
        n = o.npars() 
        if  isinstance ( index , slice ) :
            return tuple ( o.par(i) for i in range( *index.indices ( n ) ) )
        #
        if  index <  0 :  index += n
        if not 0 <= index < n :
            raise IndexError('[%s] index out of range [0,%d)' % ( index , n ) ) 
        return o.par ( index )

    try:
        ix , iy = index
    except :        
        raise IndexError('Invalid index %s/%s' % ( index  , type(index) ) )
    
    if o.index ( ix , iy ) not in o :
        raise IndexError('Invalid index (%s,%s)' % ( ix , iy ) )
    return o.par ( ix , iy )
        

for  p in ( Ostap.Math.BSpline2D           ,
            Ostap.Math.BSpline2DSym        ) :
    
    p.__setitem__  = _p2_set_par_
    p.__getitem__  = _p2_get_par_

# =============================================================================
_decorated_classes_ = set( [
    Ostap.Math.BSpline              ,
    Ostap.Math.PositiveSpline       ,
    Ostap.Math.ConvexOnlySpline     ,
    Ostap.Math.MonotonicSpline     ,
    Ostap.Math.ConvexSpline         ,  
    Ostap.Math.BSpline2D            ,
    Ostap.Math.BSpline2DSym         ,
    Ostap.Math.PositiveSpline2D     ,
    Ostap.Math.PositiveSpline2DSym  ,
    ])
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
