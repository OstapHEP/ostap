#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Set of useful "geometry" utilities 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2009-09-12
# =============================================================================
""" Set of useful ``geometry'' utilities 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2009-09-12"
__version__ = "Version$Revision$"
# =============================================================================
__all__     = ()
# =============================================================================
import ROOT, cppyy 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.geometry' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
cpp = cppyy.gbl 

## C++ namespace Ostap 
Ostap = cpp.Ostap

## ROOT::Math namespace
_RM = ROOT.ROOT.Math

## Geometry vectors 
Ostap.XYZPoint            = _RM.PositionVector3D     ('ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag')
Ostap.XYZVector           = _RM.DisplacementVector3D ('ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag')
Ostap.Plane3D             = _RM.Plane3D

Ostap.Math.XYZPoint       = Ostap.XYZPoint
Ostap.Math.XYZVector      = Ostap.XYZVector
Ostap.Math.Plane3D        = Ostap.Plane3D

Ostap.Point3D             = Ostap.XYZPoint
Ostap.Math.Point3D        = Ostap.XYZPoint

Ostap.Vector3D            = Ostap.XYZVector
Ostap.Math.Vector3D       = Ostap.XYZVector

## Ostap::Math::Line 
Ostap.Math.XYZLine        = Ostap.Math.Line( Ostap.XYZPoint, Ostap.XYZVector )
Ostap.XYZLine             = Ostap.Math.XYZLine
Ostap.Line3D              = Ostap.Math.XYZLine
Ostap.Math.Line3D         = Ostap.Math.XYZLine


## ============================================================================
## some useful decoration:
## ============================================================================


## 3-vectors 
_P3D = Ostap.XYZPoint
_V3D = Ostap.XYZVector

## 3-vectors & points

## ============================================================================
def _v3_iadd_  ( s , other ) :
    """Increment 4-vector with another 3-vector
    >>> v3  = ...
    >>> v3 += other
    """
    s.SetX     ( s.X () + other.X () )
    s.SetY     ( s.Y () + other.Y () )
    s.SetZ     ( s.Z () + other.Z () )
    return s

## ============================================================================
def _v3_isub_  ( s , other ) :
    """Increment 4-vector with another 3-vector
    >>> v3  = ...
    >>> v3 -= other
    """
    s.SetX     ( s.X () - other.X () )
    s.SetY     ( s.Y () - other.Y () )
    s.SetZ     ( s.Z () - other.Z () )
    return s

## ============================================================================
def _v3_dot_   ( s , other ) :
    """``dot''-product of two 3-vectors
    >>> v3    = ...
    >>> other = ...
    >>> print 'Dot is ' , p3.Dot ( other )
    """
    res  = s.X ( ) * other.X ( )
    res += s.Y ( ) * other.Y ( )
    res += s.Z ( ) * other.Z ( )
    return res 

## ============================================================================
def _v3_cross_  ( s , other ) :
    """``cross''-product of two 3-vectors  a * ( b x c ) 
    >>> v3    = ...
    >>> other = ...
    >>> print 'Cross product is is ' , p3.cross ( other )
    """
    bx =     s . X ()
    by =     s . Y ()
    bz =     s . Z ()
    cx = other . X ()
    cy = other . Y ()
    cz = other . Z ()
    ##
    return _V3D( by * cz - bz * cy ,
                 bz * cx - bx * cz ,
                 bx * cy - by * cx )


## ============================================================================
def _v3_triple_  ( s , b , c ) :
    """``triple''-product of three 3-vectors
    >>> v1    = ...
    >>> v2    = ...
    >>> v3    = ...
    >>> print 'Triple is ' , v1.triple ( v2  , v3 ) 
    """
    return s.Dot ( b.cross ( c ) )
                 
if not hasattr ( _V3D , '__iadd__' ) : _V3D. __iadd__ = _v3_iadd_
if not hasattr ( _V3D , '__isub__' ) : _V3D. __isub__ = _v3_isub_
if not hasattr ( _V3D , 'Dot'      ) : _V3D.Dot       = _v3_dot_
if not hasattr ( _V3D , 'dot'      ) : _V3D.dot       = _v3_dot_
if not hasattr ( _V3D , 'Cross'    ) : _V3D.Cross     = _v3_cross_
if not hasattr ( _V3D , 'cross'    ) : _V3D.cross     = _v3_cross_
if not hasattr ( _V3D , 'Triple'   ) : _V3D.Triple    = _v3_triple_
if not hasattr ( _V3D , 'triple'   ) : _V3D.triple    = _v3_triple_

if not hasattr ( _P3D , '__iadd__' ) : _P3D. __iadd__ = _v3_iadd_ 
if not hasattr ( _P3D , '__isub__' ) : _P3D. __isub__ = _v3_isub_ 


# =============================================================================
def cross ( b  , c ) :
    """``cross''-product of two 3-vectors : a = b x c 
    >>> a = ... 
    >>> b = ... 
    >>> print 'Cross product is ' , cross ( b , c )
    """
    return b.cross  ( c )


# =============================================================================
def triple ( a , b  , c ) :
    """``triple''-product of three 3-vectors   : a * (b x c) 
    >>> a = ... 
    >>> b = ... 
    >>> c = ... 
    >>> print 'Triple product is ' , triple ( a , b , c )
    """
    return a.triple ( b , c )


## ============================================================================
def _p3_as_v3_ ( self ) :
    """ Conversion to 3D-vector
    >>> print p3.asV3() 
    """
    return _V3D( self.x() , self.y() , self.z() )

## ============================================================================
def _v3_as_p3_ ( self ) :
    """ Conversion to 3D-point
    >>> print v3.asP3() 
    """    
    return _P3D( self.x() , self.y() , self.z() )

_P3D. asV3 = _p3_as_v3_
_V3D. asP3 = _v3_as_p3_

# =============================================================================
if not hasattr ( Ostap.Math , 'Vector3' ) :
    import ostap.math.linalg
    
_V3 = Ostap.Math.Vector3
# =============================================================================
## convert 3D-Vector/3D-point into SVector
#  @code
#  lv = ...
#  v3 = lv.asSVector()
#  @endcode 
def _v3_as_v3_ ( self ) :
    """Convert 3D-Vector/3D-point into SVector<double,3>
    >>> lv = ...
    >>> v3 = lv.asSVector()
    """
    _v3 = _V3()
    _v3[0] = self.X()
    _v3[1] = self.Y()
    _v3[2] = self.Z()
    return _v3

_V3D.asSVector = _v3_as_v3_ 
_P3D.asSVector = _v3_as_v3_ 

# =============================================================================
def _p3_add_ ( self , other ) :
    """Addition of 3D-point and 3D-vector    
    >>> point  = ...
    >>> vector = ...
    >>> result = point + vector 
    """
    # POINT + VECTOR = POINT 
    if isinstance ( other , _V3D ) :
        tmp   = _P3D ( self )
        tmp  += other
        return tmp
    #
    return NotImplemented


# =============================================================================
def _p3_sub_ ( self , other ) :
    """Substraction of 3D-points 
    
    >>> point1 = ...
    >>> point2 = ...
    >>> vector = ...
    >>> result_point  = point1 - vector
    >>> result_vector = point1 - point2  
    """
    # POINT - VECTOR = POINT
    if   isinstance ( other , _V3D ) :
        tmp   = _P3D ( self )
        tmp  -= other
        return tmp
    # POINT - POINT = VECTOR 
    elif isinstance ( other , _P3D ) :
        tmp   = _V3D (  self.x() ,  self.y() ,  self.z() )
        ## tmp  -= other
        tmp  -= _V3D ( other.x() , other.y() , other.z() ) 
        return tmp
    #
    return NotImplemented 

# =============================================================================
def _p3_mul_ ( self , other ) :
    """Scaling of 3D-points 
    
    >>> point  = ...
    >>> result = point1 * 2 
    """
    tmp  = _P3D ( self )
    tmp *= other
    return tmp

# =============================================================================
def _p3_div_ ( self , other ) :
    """Scaling of 3D-points 
    
    >>> point  = ...
    >>> result = point1 / 2 
    """
    tmp  = _P3D ( self )
    tmp /= other
    return tmp
    
# =============================================================================
def _v3_add_ ( self , other ) :
    """Addition  of 3D-vectors
    
    >>> vector1 = ...
    >>> vector2 = ...
    >>> result_vector = vector1 + vector2
    >>> point   =
    >>> result_point  = vector1 + point 
    """
    # VECTOR + VECTOR = VECTOR 
    if   isinstance ( other , _V3D ) :
        tmp   = _V3D ( self )
        tmp  += other
        return tmp
    # VECTOR + POINT  = POINT 
    elif isinstance ( other , _P3D ) : return other + self
    #
    return NotImplemented 

# =============================================================================
def _v3_sub_ ( self , other ) :
    """Subtraction  of 3D-vectors
    
    >>> vector1 = ...
    >>> vector2 = ...
    >>> result_vector = vector1 - vector2
    """
    # VECTOR - VECTOR = VECTOR 
    if   isinstance ( other , _V3D ) :
        tmp   = _V3D ( self )
        tmp  -= other
        return tmp
    #
    return NotImplemented 

# =============================================================================
def _v3_mul_ ( self , other ) :
    """Multiplication  of 3D-vectors
    
    >>> vector1 = ...
    >>> result  = vector1 * 2 
    >>> vector2 = ...
    >>> product = vector1 * vector2
    """
    # VECTOR * VECTOR = NUMBER
    if   isinstance ( other , _V3D ) : return self.Dot ( other )
    # VECTOR * NUMBER = NUMBER 
    elif isinstance ( other , ( float , int , long ) ) :  
        tmp  = _V3D ( self )
        tmp *= other
        return tmp
    #
    return NotImplemented

# =============================================================================
def _v3_div_ ( self , other ) :
    """Scaling of 3D-vectors
    
    >>> vector = ...
    >>> result = vector1 / 2 
    """
    tmp  = _V3D ( self )
    tmp /= other
    return tmp

_P3D . __add__  = _p3_add_
_P3D . __sub__  = _p3_sub_
_P3D . __div__  = _p3_div_
_P3D . __mul__  = _p3_mul_

_V3D . __add__  = _v3_add_
_V3D . __sub__  = _v3_sub_
_V3D . __div__  = _v3_div_
_V3D . __mul__  = _v3_mul_

_P3D . __radd__ = lambda s,o : s+o 
_P3D . __rmul__ = lambda s,o : s*o 
_V3D . __radd__ = lambda s,o : s+o 
_V3D . __rmul__ = lambda s,o : s*o 


# =============================================================================
def _v3_pow_ ( self , e ) :
    """Squared length of the 3D-vector/3D-point
    >>> print 'R-squared is ', v3**2 
    """
    if 2 != e : return NotImplemented
    return self.Mag2 ()

_V3D.__pow__ = _v3_pow_
_P3D.__pow__ = _v3_pow_

# =============================================================================
## Self-printout of 3D-points and 3D-vectors
def _v3_str_ ( self , fmt = "(%g,%g,%g) ") :
    """Self-printout of 3D-points and 3D-vectors
    >>> print p3 
    """
    return fmt % ( self.X() , self.Y( ), self.Z() )

# =============================================================================
if not hasattr ( _P3D , '_new_str_' ) :
    _P3D . _new_str_ = _v3_str_
    _P3D . __str__   = _v3_str_
    _P3D . __repr__  = _v3_str_

# =============================================================================
if not hasattr ( _V3D , '_new_str_' ) :
    _V3D . _new_str_ = _v3_str_
    _V3D . __str__   = _v3_str_
    _V3D . __repr__  = _v3_str_

# =============================================================================
## Self-printout of line
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
def _l_str_ ( self ) :
    """Self-printout of line: (point, direction)
    >>> line = ... 
    >>> print line 
    """
    return "(%s,%s)" % ( self.beginPoint() , self.direction() )

if not hasattr ( Ostap.Math.XYZLine , '_new_str_' ) :
    Ostap.Math.XYZLine._new_str_ = _l_str_
    Ostap.Math.XYZLine.__str__   = _l_str_
    Ostap.Math.XYZLine.__repr__  = _l_str_


# =============================================================================
## Self-printout of 3D-plane
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
def _p_str_ ( self ) :
    """Self-printout of 3D-plane: (point, normal)
    >>> plane = ...
    >>> print plane 
    """
    return "(%s,%s)" % ( self.ProjectOntoPlane( Ostap.XYZPoint()) , self.Normal() )

if not hasattr ( Ostap.Plane3D , '_new_str_' ) :
    Ostap.Plane3D._new_str_ = _p_str_
    Ostap.Plane3D.__str__   = _p_str_
    Ostap.Plane3D.__repr__  = _p_str_

# =============================================================================
## various decorators for GeomFun.h
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
if not hasattr ( Ostap.Math , 'XYZGeomFun' ) :
    Ostap.Math.XYZGeomFun = Ostap.Math.GF(
        Ostap.XYZPoint  ,
        Ostap.XYZLine   ,
        Ostap.Plane3D
        )
    
if not hasattr ( Ostap , 'XYZGeomFun' ) :
    Ostap.XYZGeomFun = Ostap.Math.XYZGeomFun

_GeomFun = Ostap.Math.XYZGeomFun


# =============================================================================
## intersection of line and plane
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _intersect_line_and_plane_ ( line , plane ) :
    """Find the intersection of line and plane

    >>> line  = ...
    >>> plane = ...

    >>> ok, point, mu = line.intersect ( plane )

    The return value is a tuple:
    - the point
    - the parameter along the line
    - the flag (true if intersection exists)

    """
    _point = Ostap.XYZPoint(0,0,-1.e+10)
    _mu    = ROOT.Double(-1.e+10)
    _flag  = _GeomFun.intersection ( line   ,
                                     plane  ,
                                     _point ,
                                     _mu    )
    if _flag : _flag = True
    else     : _flag = False
    return (_point,_mu,_flag)

_intersect_line_and_plane_ . __doc__ += '\n' + _GeomFun.intersection . __doc__

if not hasattr ( Ostap.XYZLine , 'intersect' ) :
    Ostap.XYZLine.intersect = _intersect_line_and_plane_

# =============================================================================
## intersect two planes
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _intersect_two_planes_ ( plane , plane1 ) :
    """Find the intersection line for two planes:

    >>> plane  = ...
    >>> plane1 = ...
    >>> line, flag = plane.line(plane1)

    Return value is a tuple:

    - the intersection line
    - the flag (true if intersection exists)

    """
    _line = Ostap.XYZLine()
    _flag = _GeomFun.intersection ( plane , plane1 , _line )
    if _flag : _flag = True
    else     : _flag = False
    return (_line,_flag)

_intersect_two_planes_ . __doc__ += '\n' + _GeomFun.intersection . __doc__

if not hasattr ( Ostap.Plane3D , 'line' ) :
    Ostap.Plane3D.line = _intersect_two_planes_


# =============================================================================
## intersect three planes
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _intersect_three_planes_ ( plane , plane1 , plane2 ) :
    """Find the intersection point for three planes:

    >>> plane  = ...
    >>> plane1 = ...
    >>> plane3 = ...
    >>> point, flag = plane.point(plane1,plane2)

    Return value is a tuple:

    - the intersection point
    - the flag (true if intersection exists)

    """
    _point = Ostap.XYZPoint(0,0,-1.e+10)
    _flag = _GeomFun.intersection ( plane , plane1 , plane2 , _point )
    if _flag : _flag = True
    else     : _flag = False
    return (_point,_flag)


_intersect_three_planes_ . __doc__ += '\n' + _GeomFun.intersection . __doc__

if not hasattr ( Ostap.Plane3D , 'point' ) :
    Ostap.Plane3D.point = _intersect_three_planes_



# =============================================================================
## intersect the planes
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _intersect_the_planes_ ( plane , plane1 , plane2 = None ) :
    """Find the intersection line/point for two or three planes:

    >>> plane  = ...
    >>> plane1 = ...
    >>> line, flag = plane.intersect(plane1)

    >>> plane  = ...
    >>> plane1 = ...
    >>> plane2 = ...
    >>> point, flag = plane.intersect(plane1,plane2)

    Return value is a tuple:

    - the intersection line/point
    - the flag (true if intersection exists)

    """
    if not plane2 : return _intersect_two_planes_ ( plane , plane1 )
    return _intersect_three_planes_ ( plane , plane1 , plane2 )

_intersect_the_planes_ . __doc__ += '\n' + _GeomFun.intersection . __doc__

if not hasattr ( Ostap.Plane3D , 'intersect' ) :
    Ostap.Plane3D.intersect = _intersect_the_planes_





# =============================================================================
## calculate the impact parameter of the line & point
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _imp_par_1_ ( line , point ) :
    """Calculate the impact parameter of the line and the point
    >>> line  = ...
    >>> point = ...
    >>> ip = line.impactParameter ( point )
    """
    return _GeomFun.impactParameter ( point , line )

_imp_par_1_ . __doc__ += '\n' + _GeomFun.impactParameter . __doc__

if not hasattr ( Ostap.XYZLine , 'impactParameter' ) :
    Ostap.XYZLine.impactParameter = _imp_par_1_
if not hasattr ( Ostap.XYZLine , 'ip'              ) :
    Ostap.XYZLine.ip              = _imp_par_1_

# =============================================================================
## calculate the impact parameter of the line & point
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _imp_par_2_ ( point , line ) :
    """Calculate the impact parameter of the line and the point
    >>> point = ...
    >>> line  = ...
    >>> ip = point.impactParameter ( line )
    """
    return _GeomFun.impactParameter ( point , line )


_imp_par_2_ . __doc__ += '\n' + _GeomFun.impactParameter . __doc__

if not hasattr ( Ostap.XYZPoint , 'impactParameter' ) :
    Ostap.XYZPoint.impactParameter = _imp_par_2_
if not hasattr ( Ostap.XYZPoint , 'ip'              ) :
    Ostap.XYZPoint.ip              = _imp_par_2_


# =============================================================================
## distance between two lines
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _distance_between_two_lines_ ( line , line1 ) :
    """Find the distance between two lines :
    >>> line  = ...
    >>> line1 = ...
    >>> dist = line.distance ( line1 )
    """
    return _GeomFun.distance ( line , line1 )

_distance_between_two_lines_ . __doc__ += '\n' + _GeomFun.distance. __doc__

if not hasattr ( Ostap.XYZLine , 'distance' ) :
    Ostap.XYZLine.distance = _distance_between_two_lines_


# =============================================================================
## find the closest points for two lines
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _closest_points_ ( line , line1 ) :
    """Calculate the two closest points between two lines

    >>> line1 = ...
    >>> line2 = ...
    >>> point1 , point2 , flag = line1.closestPoints ( line2 )

    The return values is a tuple:
    - the point onthe fist line
    - the point on the second line
    - the flag (true is everything OK)
    """
    _point1 = Ostap.XYZPoint(0,0,-1.e+10)
    _point2 = Ostap.XYZPoint(0,0,-1.e+11)
    _flag   = _GeomFun.closestPoints ( line , line1 , _point1 , _point2 )
    if    _flag : _flag = True
    else        : _flag = False
    return (_point1,_point2,_flag)


_closest_points_ . __doc__ += '\n' + _GeomFun.closestPoints . __doc__

if not hasattr ( Ostap.XYZLine , 'closestPoints' ) :
   Ostap.XYZLine.closestPoints = _closest_points_


# =============================================================================
## find the closest points for two lines
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _closest_point_params_ ( line , line1 ) :
    """Calculate the parameters for two closest points between two lines

    >>> line1 = ...
    >>> line2 = ...
    >>> mu1 , mu2 , flag = line1.closestPointParams ( line2 )

    The return values is a tuple:
    - the 'mu-parameter of closest point along the first  line
    - the 'mu-parameter of closest point along the second line
    - the flag (true is everything OK)
    """
    _mu1    = ROOT.Double(-1.e+10)
    _mu2    = ROOT.Double(-1.e+11)
    _flag   = _GeomFun.closestPointParams ( line , line1 , _mu1 , _mu2 )
    if    _flag : _flag = True
    else        : _flag = False
    return (_mu1,_mu2,_flag)


_closest_point_params_ . __doc__ += '\n' + _GeomFun.closestPointParams . __doc__

if not hasattr ( Ostap.XYZLine , 'closestPointParams' ) :
    Ostap.XYZLine.closestPointParams = _closest_point_params_

# =============================================================================
## find the point on ilne closest to the given point
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _closest_point_1_ ( line , point ) :
    """Find the point on line closest to the given point
    >>> line  = ...
    >>> point = ...
    >>> ClosestPoint  = line.closestPoint ( point )
    """
    return _GeomFun.closestPoint ( point , line )

_closest_point_1_ . __doc__ += '\n' + _GeomFun.closestPoint . __doc__

if not hasattr ( Ostap.XYZLine , 'closestPoint' ) :
    Ostap.XYZLine.closestPoint = _closest_point_1_


# =============================================================================
## find the point on ilne closest to the given point
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _closest_point_2_ ( point , line ) :
    """Find the point on line closest to the given point
    >>> point = ...
    >>> line  = ...
    >>> ClosestPoint  = point.closestPoint ( line )
    """
    return _GeomFun.closestPoint ( point , line )

_closest_point_2_ . __doc__ += '\n' + _GeomFun.closestPoint . __doc__

if not hasattr ( Ostap.XYZPoint , 'closestPoint' ) :
    Ostap.XYZPoint.closestPoint = _closest_point_2_


# =============================================================================
## find the parameter along the line to the closest point to the given point
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _closest_point_param_1_ ( line , point ) :
    """Find the parameter along the line to the closest point
    >>> line  = ...
    >>> point = ...
    >>> mu = line.closestPointParam ( point )

    """
    return _GeomFun.closestPointParam ( point , line )


_closest_point_param_1_ . __doc__ += '\n' + _GeomFun.closestPointParam .__doc__

if not hasattr ( Ostap.XYZLine , 'closestPointParam' ) :
    Ostap.XYZLine.closestPointParam = _closest_point_param_1_


# =============================================================================
## find the parameter along the line to the closest point to the given point
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _closest_point_param_2_ ( point , line ) :
    """Find the parameter along the line to the closest point
    >>> point = ...
    >>> line  = ...
    >>> mu = point.closestPointParam ( line )
    """
    return _GeomFun.closestPointParam ( point , line )

_closest_point_param_2_ . __doc__ += '\n' + _GeomFun.closestPointParam .__doc__

if not hasattr ( Ostap.XYZPoint , 'closestPointParam' ) :
    Ostap.XYZPoint.closestPointParam = _closest_point_param_2_

# =============================================================================
## check if two lines are parallel
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-10-22
def _parallel_lines_ ( line , line1 ) :
    """Check if two lines are parallel:
    >>> line  = ...
    >>> line1 = ...
    >>> par   = line.parallel ( line1 )
    """
    _flag = _GeomFun.parallel ( line , line1 )
    if not _flag : return False
    return True

_parallel_lines_ . __doc__ += '\n' + _GeomFun.parallel . __doc__

if not hasattr ( Ostap.XYZLine , 'parallel' ) :
    Ostap.XYZLine.parallel = _parallel_lines_

# =============================================================================
_decorated_classes_ = set( [
    Ostap.Math.XYZPoint          , 
    Ostap.Math.XYZVector         ,
    Ostap.Math.LorentzVector     ,
    Ostap.Math.Plane3D           ,
    Ostap.Math.Point3D           , 
    Ostap.Vector3D               ,
    Ostap.XYZLine                ,
    Ostap.LorentzVector          ,
    ] )

_new_methods_ = (
    _intersect_line_and_plane_   ,
    _intersect_two_planes_       , 
    _intersect_three_planes_     , 
    _intersect_the_planes_       , 
    _imp_par_1_                  , 
    _imp_par_2_                  ,
    _distance_between_two_lines_ ,
    _closest_points_             , 
    _closest_point_params_       , 
    _closest_point_1_            , 
    _closest_point_2_            , 
    _closest_point_param_1_      , 
    _closest_point_param_2_      , 
    _parallel_lines_             , 
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    p1  = Ostap.XYZPoint      (0,1,2)
    v1  = Ostap.XYZVector     (2,1,0)
    l1  = Ostap.XYZLine       ( p1 , v1 )
    pl1 = Ostap.Plane3D       ( 1, 1, 1, 0 )

    logger.info ( '3D-point       (x,y,z)           : %s' % p1  )
    logger.info ( '3D-vector      (x,y,z)           : %s' % v1  )
    logger.info ( '3D-line        (point,direction) : %s' % l1  )
    logger.info ( '3D-plane       (point,normal)    : %s' % pl1 )

    pnt1   = Ostap.XYZPoint(-1,-2,-3)
    pnt2   = Ostap.XYZPoint( 1,-2,-3)
    line1  = Ostap.Math.XYZLine(Ostap.XYZPoint(0,1,2), Ostap.XYZVector(1,1,1)  )
    line2  = Ostap.Math.XYZLine(Ostap.XYZPoint(1,3,0), Ostap.XYZVector(1,-1,2) )
    plane1 = Ostap.Math.Plane3D ( 0 , 1, 2,  3 )
    plane2 = Ostap.Math.Plane3D ( 1 , 8, 9,  0 )
    plane3 = Ostap.Math.Plane3D ( 4 , 5, 6, -1 )
    
    logger.info ( ' line  : intersect          : %s' % ( line1.intersect ( plane1 )     ,) ) 
    logger.info ( ' plane : intersect two      : %s' % ( plane1.line(plane2)            ,) )
    logger.info ( ' plane : intersect three    : %s' % ( plane1.point(plane2,plane3)    ,) )
    logger.info ( ' plane : intersect two      : %s' % ( plane1.intersect(plane2)       ,) )
    logger.info ( ' plane : intersect three    : %s' % ( plane1.intersect(plane2,plane3),) )
    logger.info ( ' line  : impactParameter    : %s' % ( line1.impactParameter ( pnt1 ) ,) )
    logger.info ( ' point : impactParameter    : %s' % ( pnt1 .impactParameter ( line1 ),) )
    logger.info ( ' line  : distance           : %s' % ( line1.distance  (line2)        ,) )
    logger.info ( ' line  : closestPoints      : %s' % ( line1.closestPoints(line2)     ,) )
    logger.info ( ' line  : closestPointParams : %s' % ( line1.closestPointParams(line2),) ) 
    logger.info ( ' line  : closestPoint       : %s' % ( line1.closestPoint ( pnt2  )   ,) )
    logger.info ( ' point : closestPoint       : %s' % ( pnt2 .closestPoint ( line1 )   ,) )
    logger.info ( ' line  : closestPointParam  : %s' % ( line1.closestPointParam ( pnt1),) ) 
    logger.info ( ' point : closestPointParam  : %s' % ( pnt1 .closestPointParam (line1),) ) 
    logger.info ( ' line  : parallel           : %s' % ( line1.parallel  (line2)        ,) )
    
    logger.info ( 80*'*' ) 

# =============================================================================
# The END 
# =============================================================================
