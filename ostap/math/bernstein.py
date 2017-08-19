#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file bernstein.py
#
#  Module with some useful utilities for dealing with Bernstein polynomials
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-12-01
# =============================================================================
""" Module with some useful utilities for dealing with Bernstein polynomials
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = (
    'control_polygon'   , ## get a control polygon for Bernstein polynomial
    'upper_convex_hull' , ## upper convex hull for Bernstein polynomial
    'lower_convex_hull' , ## lower convex hull for Bernstein polynomial
    'convex_hull'       , ##       convex hull for Bernstein polynomial
    #
    'sturm_sequence'    , ## get Sturm's sequence for Bernstein polynomial
    'nroots'            , ## number of roots (using  Sturm's sequence)
    'sign_changes'      , ## number of sign changes in sequence of coefficients 
    #
    'deflate_left'      , ## deflate Berntein polynomial at x=xmin
    'deflate_right'     , ## deflate Berntein polynomial at x=xmax
    'deflate'           , ## deflate Berntein polynomial at x=x0
    #
    'crossing_points'   , ## get crossing points of control polygon with x-axis
    'left_line_hull'    , ## get the most left point of crossing of the convex hull with x-axis
    'right_line_hull'   , ## get the most right point of crossing of the convex hull with x-axis
    )
# =============================================================================
import  ROOT, math  
from    ostap.core.core import cpp, Ostap, funID
from    ostap.math.base import iszero, isequal, signum  
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.bernstein' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
## short name 
Bernstein = cpp.Ostap.Math.Bernstein
# =============================================================================
## get control polygon for Bernstein polynomial   
def control_polygon ( bp )  :
    """Get control polygon for Bernstein polynomial
    >>>  bernstein = ...
    >>>  cp = bernstein.control_polygon ()
    >>>  cp = control_polygon( bernstein )  ##  ditto 
    """
    return  Ostap.Math.control_polygon (  bp.bernstein() ) 


# =============================================================================
##  get Sturm sequence for bernstein polynomial 
##  @see https://en.wikipedia.org/wiki/Sturm%27s_theorem
def sturm_sequence ( bp ) :
    """ Get Sturm's sequence for Bernstein polynomials.
    >>> bernstein = ...
    >>> sturm_seq = bernstein.sturm_sequence () 
    >>> sturm_seq = sturm_sequence( bernstein )  ## ditto
    - see https://en.wikipedia.org/wiki/Sturm%27s_theorem
    """
    bp = bp.bernstein()
    
    ## isolate roots
    b0  = bp
    b1  = bp.derivative()
    gcd = b0.gcd ( b1 )
    
    ## square free polynom
    b0  = b0 // gcd 
    b1  = b0.derivative()
    
    ss = [ b0 , b1 ]
    
    for i in range ( bp.degree() ) :
        b0   = ss[-2]
        b1   = ss[-1] 
        a, b = divmod ( b0 , b1 ) 
        ## check norms
        n0 = b0.norm()
        n1 = b1.norm()
        na = a.norm ()
        nt = ( na * n1 + n0 ) / 10 
        if b.small ( nt ) :
            for j in range (  b.npars() ) : b.setPar( j , 0.0 ) 
        ss.append  ( -b ) 

    return tuple(ss)

# ============================================================================
## get number of unique real roots between xmin and xmax
#  @see https://en.wikipedia.org/wiki/Sturm%27s_theorem
def nroots ( bp , xmin =  None , xmax = None ) :
    """ Get the number of unique real roots between xmin and xmax
    >>> bernstein = ...
    >>> n_roots = bernstein.nroots  ( 1 , 100 )
    >>> n_roots = nroots( bernstein , 1 , 100 )  ## ditto
    - see https://en.wikipedia.org/wiki/Sturm%27s_theorem
    """
    bp     = bp.bernstein()
    ss     = bp.sturm_sequence()

    if xmin is None : xmin = bp.xmin ()
    if xmax is None : xmax = bp.xmax ()
    
    smn   =  [ b.evaluate ( xmin ) for b in ss ]
    smx   =  [ b.evaluate ( xmax ) for b in ss ]

    ncmn  = 0
    ncmx  = 0
    nss   = len(ss)
    
    for  i in range(1,nss) :
        if 0 > smn[i]*smn[i-1] : ncmn +=1
        if 0 > smx[i]*smx[i-1] : ncmx +=1

    ## print smx , smn
    ## print ncmn, ncmx
     
    return ( ncmn - ncmx ) if xmin < xmax else ( ncmx - ncmn )  

# =============================================================================
## get upper convex hull for  Bernstein polynomial
def upper_convex_hull ( bp ) :
    """Get upper convex hull for  Bernstein polynomial
    >>> bernstein = ...
    >>> upper = bernstein.upper_convex_hull  ()
    >>> upper = upper_convex_hull ( bernstein ) ## ditto
    - for all values of x: bernstein(x) <=  upper ( x ) 
    """
    return Ostap.Math.upper_convex_hull ( bp.bernstein()  )
# =============================================================================
## get lower convex hull for  Bernstein polynomial
def lower_convex_hull ( bp ) :
    """Get lower convex hull for  Bernstein polynomial
    >>> bernstein = ...
    >>> lower = bernstein.lower_convex_hull ()
    >>> lower = lower_convex_hull ( bernstein )
    - for all values of x: bernstein(x) >=  lower ( x ) 
    """
    return Ostap.Math.lower_convex_hull ( bp.bernstein()  )
# =============================================================================
## get upper & lower convex hulls for  Bernstein polynomial
def  convex_hull ( bp ) :
    """Get lower& upper convex hulls for  Bernstein polynomial
    >>> bernstein = ...
    >>> lower,upper = bernstein.convex_hull  ()
    >>> lower,upper = convex_hull ( bernstein ) ## ditto  
    - for all values of x: lower(x) <= bernstein(x) <= upper ( x ) 
    """
    bp1 = bp.bernstein()
    return bp1.lower_convex_hull () , bp1.upper_convex_hull ()

# =============================================================================
## Deflate Bernstein polynomial at  <code>x=xmin</code>
#  \f$ b(x)-b(x_{min})=(x-x_{min})*d(x)\f$      
#   @param  b  berntein polynomial to be deflated 
#   @return deflated polinomial "d"
def deflate_left ( bp ) :
    """Deflate Bernstein polynomial at  x = x min
    b(x)-b(xmin)=(x-xmin)*d(x)
    >>> berstein = ...
    >>> deflated = berstein.deflate_left ()
    >>> deflated = deflate_left( bernstein ) ## ditto
    """
    return Ostap.Math.deflate_left ( bp.bernstein() )   

# =============================================================================
## Deflate Bernstein polynomial at  <code>x=xmax</code>
#  \f$ b(x)-b(x_{max})=(x-x_{max})*d(x)\f$      
#   @param  b  berntein polynomial to be deflated 
#   @return deflated polinomial "d"
def deflate_right  ( bp ) :
    """Deflate Bernstein polynomial at  x = xmax
    b(x)-b(xmax)=(x-xmax)*d(x)
    >>> berstein = ...
    >>> deflated = berstein.deflate_right ()
    >>> deflated = deflate_right ( bernstein ) ## ditto
    """
    return Ostap.Math.deflate_right ( bp.bernstein() )   

# ===============================================================================
## get abscissas of crosssing point of the control polygon 
#  for Bernstein polynomial
#  @param  b bernstein polynomial
#  @return abscissas of crossing points of the control  polygon with x-axis 
def crossing_points  ( bp ) :
    """Get abscissas of crosssing point of the control polygon  for Bernstein polynomial
    >>> bernstein = ...
    >>> cps = bernstein.crossing_points()
    >>> cps = crossing_points( berstein ) ## ditto 
    """
    cps = Ostap.Math.crossing_points ( bp.bernstein() )
    return tuple ( cp for cp in cps ) 

# =============================================================================
## Deflate Bernstein polynomial at  <code>x=x0</code>
#  \f$ b(x)-b(x_{0})=(x-x_{0})*d(x)\f$      
#   @param  b  berntein polynomial to be deflated 
#   @return deflated polinomial "d"
def deflate ( bp , x0 ) :
    """Deflate Bernstein polynomial at  x = x0
    b(x)-b(x0)=(x-x0)*d(x)
    >>> berstein = ...
    >>> deflated = berstein.deflate ( 0.3 )
    >>> deflated = deflate ( bernstein , 0.3 ) ## ditto
    """
    return Ostap.Math.deflate ( bp.bernstein() , x0 ) 

# =============================================================================
## get the most left crossing  point of convex hull with  x-axis 
# (it is a step  towards finding the most left root, if any 
# if convex hull does not cross the x-axis, None is returned      
def left_line_hull ( bp  ) :
    """ Get the most left crossing  point of convex hull with  x-axis 
    (it is a step  towards finding the most left root, if any 
    if convex hull does not cross the x-axis, None is returned
    >>> bernstein = ...
    >>> lcp = bernstein.left_line_hull() 
    >>> lcp = left_line_hull( bernstein ) ##  ditto 
    """
    lcp = Ostap.Math.left_line_hull ( bp.bernstein() )
    ##
    if not bp.xmin() <= lcp <= bp.xmax() :
        return None
    
    return lcp 

# =============================================================================
## get the most right crossing  point of convex hull with  x-axis 
# (it is a step  towards finding the most right root, if any 
# if convex hull does not cross the x-axis, None is returned      
def right_line_hull ( bp  ) :
    """ Get the most right crossing  point of convex hull with  x-axis 
    (it is a step  towards finding the most right root, if any 
    if convex hull does not cross the x-axis, None is returned
    >>> bernstein = ...
    >>> rcp = bernstein.right_line_hull() 
    >>> rcp = right_line_hull( bernstein ) ##  ditto 
    """
    rcp = Ostap.Math.right_line_hull ( bp.bernstein() )
    ##
    if not bp.xmin() <= rcp <= bp.xmax() :
        return None
    #
    return rcp 

# =============================================================================
## get number of (strickt) sign changes in trhe sequnce of coefficients
#  for Bernstein polynomial 
#  if  N is number of sign changes, then the number of real roots R is 
#  \f$ R = N - 2K\f$, where K is non-negative integer
def sign_changes ( bp ) :
    """Get number of (strickt) sign changes in trhe sequnce of coefficients
    for Bernstein polynomial 
    - if  N is number of sign changes, then the number of real roots R is 
    \f$ R = N - 2K\f$, where K is non-negative integer
    >>> bernstein = ....
    >>> nr = bernstein.sign_changes ()
    >>> nr = sign_changes ( bersntein ) ## ditto    
    """
    return Ostap.Math.sign_changes ( bp.bernstein() ) 

# =============================================================================
## Root finding for Bernstein polnomials using Laguerre's method
## https://en.wikipedia.org/wiki/Laguerre%27s_method
def _laguerre_ ( x , f , d1 = None , d2 = None ) :
    """ Root finding for Bernstein polnomials using Laguerre's method
    - https://en.wikipedia.org/wiki/Laguerre%27s_method
    """
    fn   = f.norm() / 10 
    n    = f.degree()
    
    xmin = f.xmin ()
    xmax = f.xmax ()
    
    l    = xmax - xmin
    
    if not d1 : d1 =  f.derivative()
    if not d2 : d2 = d1.derivative()
    
    if  not f.xmin() <= x <= f.xmax() :
        x = 0.5 * ( xmin + xmax ) 

    if 1 ==  f.degree() :
        
        p0 = f [ 0 ]
        p1 = f [ 1 ]
        
        s0 = signum ( p0 )
        s1 = signum ( p1 )
        
        return ( p1 * xmin - p0 * xmax ) / ( p1 - p0) , 
        

    ## the convergency is cubic, therefore 16 is *very* larger number of iterations 
    for i in  range(16) : 
        
        if  not xmin <= x <= xmax :
            ## something goes wrong: multiple root? go to derivative
            break 

        vx = f.evaluate(x)
        if iszero ( vx ) or isequal ( fn +  vx , fn ) : return x,
        
        G =         d1 ( x ) / vx
        H = G * G - d2 ( x ) / vx
        
        d = math.sqrt ( ( n - 1 ) * ( n * H - G * G ) )
        if G < 0 : d = -d

        a = n / ( G + d ) 

        x -= a

        if isequal ( a + x , x ) : return x,  
        
    ## look for derivative
    r =  _laguerre_ ( x , d1 , d2 , d2.derivative() )
    if r : r = r[:1] + r
    return r

    
## ============================================================================
def _scale_ ( f ) :
    fn = f.norm()
    sf = math.frexp ( fn )[1]
    if 1 !=  sf : f = f.ldexp ( 1  - sf )
    return f


# ============================================================================
## Solve   equation  B(x) = C, where B(x) is  Bernstein polynomial, and
#  C is  e.g.  constant or another polynomial
#  @code
#  bernstein = ...
#  roots = bernstein.solve ()
#  roots = solve ( bernstein , C = 0 ) ## ditto
#  @endcode
def solve (  bp  , C = 0 , split = 2 ) :
    """Solve   equation  B(x) = C, where B(x) is  Bernstein polynomial, and
    C is  e.g.  constant or another polynomial
    >>> bernstein = ...
    >>> roots = bernstein.solve ()
    >>> roots = solve ( bernstein , C = 0 ) ## ditto
    """
    
    ## construct the equation  b'(x)=0:
    bp = bp.bernstein() 
    if C : bp = bp - C
    
    ## zero-degree polynomial
    if   0 == bp.degree() :

        if iszero ( bp[0] ) : return bp.xmin(),
        return () 
    
    ## linear polynomial
    elif 1 == bp.degree() :

        x0 = bp.xmin()
        x1 = bp.xmax()
        
        p0 = bp[0]
        p1 = bp[1]

        s0 = signum ( p0 )
        s1 = signum ( p1 )
        
        if   s0 ==     0 : return x0,  ## 
        elif s1 ==     0 : return x1,  ##
        elif s0 * s1 > 0 : return ()   ## no roots
        #
        return ( p1 * x0 - p0 * x1 ) / ( p1 - p0) , 

    ## make a copy & scale is needed 
    bp  = _scale_ ( bp + 0 ) 

    ##  norm of polynomial 
    bn = bp.norm()

    ## check number of roots
    nc = bp.sign_changes()
    if not nc : return ()      ## no roots !   RETURN

    ## treat separetely roots at the left and right edges
    roots = []
    while 1 <= bp.degree() and isequal ( bp[0] +  bn , bn ) :
        bp   -= bp[0]
        bp    = bp.deflate_left() 
        bp    = _scale_ ( bp ) 
        bn    = bp.norm ()
        roots.append ( bp.xmin() )        
    if roots : return tuple(roots) + bp.solve ( split = split ) 
    
    roots = []
    while 1 <= bp.degree() and isequal ( bp[-1] +  bn , bn ) :
        bp -= bp[-1]
        bp  = bp.deflate_right() 
        bp  = _scale_ ( bp ) 
        bn  = bp.norm () 
        roots.append ( bp.xmax() )        
    if roots : return bp.solve ( split = split ) + tuple(roots) 

    ## check again number of roots
    nc = bp.sign_changes()
    if not nc : return ()      ## no roots !   RETURN
    
    # =========================================================================
    ## there are  many roots in the interval 
    # =========================================================================


    if 1 < nc :

        xmin = bp.xmin ()
        xmax = bp.xmax ()
        
        lcp  = bp.left_line_hull()
        if ( not lcp is None ) and xmin <= lcp <= xmax : xmin = max ( xmin , lcp )
        
        rcp  = bp.right_line_hull()
        if ( not rcp is None ) and xmin <= rcp <= xmax : xmax = min ( xmax , rcp )

        # =====================================================================
        ## Two strategies for isolation of roots:
        #  - use control polygon 
        #  - use the derivative 
        #  For both cases, the zeros of contol=-polygon and/or derivatives
        #  are tested to be the roots of  polynomial..
        #  If they corresponds to the roots, polinomial is properly deflated and
        #  deflated roots are collected
        #  Remaining points are used to (recursively) split interval into smaller
        #  intervals with presumably smaller number of roots 
        
        # 
        if 0 < split :

            cps = bp.crossing_points()
            splits = [ xmin ]
            for xp in cps :
                if xmin < xp < xmax : splits.append ( xp )
            splits.append ( xmax ) 

            split -= 1
            
            ## use simple bisection
            if 2 == len ( splits ) :
                if xmin < xmax : splits =     xmin   , 0.5*(   xmin  +   xmax  ) ,    xmax
                else           : splits =  bp.xmin() , 0.5*(bp.xmin()+bp.xmax()) , bp.xmax()
                
        else :
            
            ## use the roots of derivative 
            dd = bp.derivative()
            rd = dd.solve ( split = split )
            if not rd :
                ## use bisection 
                rd = xmin, 0.5 * ( xmin + xmax ), xmax
            else :
                ## use roots of derivative 
                nrd = [ xmin , xmax ] 
                for  r in  rd :
                    if xmin < r < xmax :
                        found = False 
                        for  rr in nrd :
                            if isequal ( r , rr ) :
                                found = True
                                break 
                        if not found : nrd.append ( r )
                nrd.sort() 
                splits = nrd[:]

        splits = list(splits) 

        roots = []
        for s in splits :
            
            bv  = bp.evaluate ( s ) 
            bn = bp.norm() 
            while 1 <= bp.degree() and isequal ( bv + bn , bn ) :
                bp -= bv 
                bp  = bp.deflate  ( s  )
                bp  = _scale_     ( bp ) 
                bn  = bp.norm     (    )
                bv  = bp.evaluate ( s  ) 
                roots.append  ( s )
                for q in splits :
                    if  isequal ( q , s ) :
                        splits.remove ( q )
                
        if roots :
            roots += list ( bp.solve ( split = split - 1 ) )
            roots.sort()
            return tuple(roots)

        if 2 == len(splits) :
            if isequal ( bp.xmin() , splits[0] ) and isequal ( bp.xmax() , splits[1] ) :
                xmn    = splits[0]
                xmx    = splits[1]
                splits = [ xmn , 0.5*(xmn+xmx) , xmx ] 
                
        ns = len(splits)
        for i in range(1,ns) :
            xl = splits[i-1]
            xr = splits[i  ]            
            bb = _scale_ ( Bernstein ( bp , xl , xr ) )
            roots += bb.solve ( split = split - 1 )

        roots.sort()
        return  tuple ( roots )   ##  RETURN


    # =========================================================================
    ## there is exactly 1 root here
    # =========================================================================
    
    l0 = ( bp.xmax() -  bp.xmin() ) / ( bp.degree() + 1 )

    ## trivial case 
    if 1 == bp.degree() :        
        y0 =  bp[ 0]
        y1 =  bp[-1]
        x  = ( y1 * bp.xmin() -  y0 * bp.xmax() ) / ( y1 - y0)
        return x,
    
    ## make several iterations (not to much) for better isolation of root
    for i in range ( bp.degree() + 1 ) : 

        xmin = bp.xmin ()
        xmax = bp.xmax ()
        bn   = bp.norm ()

        lcp = bp.left_line_hull  ()
        if not lcp is None :
            if   isequal ( lcp , xmin ) : return  lcp,                      ## RETURN
            elif lcp <= xmax or isequal ( lcp , xmax ) : 
                bv =  bp.evaluate ( lcp )
                if  iszero ( bv ) or isequal ( bv + bn , bn ) : return lcp, ## RETURN
            xmin = max ( xmin , lcp )
                
        rcp = bp.right_line_hull () 
        if not rcp is None :
            if   isequal ( lcp , xmax ) : return  rcp,                      ## RETURN
            elif lcp >= xmin or isequal ( lcp , xmin ) : 
                bv =  bp.evaluate ( rcp )
                if  iszero ( bv ) or isequal ( bv + bn , bn ) : return rcp, ## RETURN
            xmax = min ( xmax , rcp )

        ## 
        if isequal ( xmin , xmax ) : return 0.5*(xmin+xmax),               ## RETURN

        if xmin >= xmax : break                             ## BREAK
        
        ## avoid too iterations - it decreased the precision 
        if 10 * ( xmax - xmin )  < l0  : break              ## BREAK  
        
        s0  =  signum ( bp[ 0] )
        s1  =  signum ( bp[-1] )
        
        smn =  signum ( bp.evaluate ( xmin ) )
        smx =  signum ( bp.evaluate ( xmax ) )

        ##  we have lost the root ?
        if ( s0 * s1 ) * ( smn * smx ) < 0 :  break           ## BREAK 

        ## refine the polynomial: 
        _bp = Bernstein ( bp , xmin , xmax )
        _bp = _scale_ ( _bp )
        _nc = _bp.sign_changes()

        ## we have lost the root 
        if not _nc : break                                   ## BREAK 

        bp = _bp 
        bn = bp.norm()

    # =========================================================================
    ## start the of root-polishing machinery
    # ========================================================================= 

    f = bp 
    d1 = f .derivative ()
    d2 = d1.derivative ()
    
    cps = bp.crossing_points()
    if cps : x = cps[0]
    else   : x = 0.5*(bp.xmin()+bp.xmax()) 

    ##  use Laguerre's method to get the isolated root 
    l = _laguerre_ ( x , f , d1 , d2 )

    return l 


# =============================================================================    
for  p in ( Ostap.Math.Bernstein     ,
            Ostap.Math.BernsteinEven ,
            Ostap.Math.Positive      ,
            Ostap.Math.PositiveEven  ) :
    
    p.lower_convex_hull =  lower_convex_hull
    p.upper_convex_hull =  upper_convex_hull
    p.convex_hull       =        convex_hull
    p.convex_hulls      =        convex_hull
    p.control_polygon   =    control_polygon
    p.crossing_points   =    crossing_points 
    p.left_line_hull    =     left_line_hull
    p.right_line_hull   =    right_line_hull
    #
    p.sturm_sequence    =     sturm_sequence 
    p.nroots            =             nroots
    p.sign_changes      =       sign_changes
    #
    p.deflate_left      =       deflate_left
    p.deflate_right     =      deflate_right 
    p.deflate           =            deflate
    #
    p.solve             =              solve


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
