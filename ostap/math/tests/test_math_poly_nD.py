#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_poly_nD.py
#  Test script for 2&3-dimnsional polynomial functions  
#  @see Ostap::Math::LegendreSum2
#  @see Ostap::Math::LegendreSum3
#  @see Ostap::Math::Bernstein2D
#  @see Ostap::Math::Bernstein3D
# ============================================================================= 
""" Test script for 2&3-dimnsional polynomial functions  
a- see Ostap::Math::LegendreSum2
- see Ostap::Math::LegendreSum3
- see Ostap::Math::Bernstein2D
- see Ostap::Math::Bernstein3D
"""
# ============================================================================= 
from   ostap.core.core          import Ostap, SE  
from   ostap.logger.colorized   import attention
from   ostap.utils.progress_bar import progress_bar 
import ostap.math.integral      as     I 
import ostap.logger.table       as     T
import ROOT, random, math 
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_poly_nD' )
else                       : logger = getLogger ( __name__                  )
# ============================================================================= 
scale = 1.e+12

# =============================================================================
## relatiev difference between two numbers  
def diff  ( a , b ) :
    """Relatiev difference between two numbers
    """
    return ( a - b ) / ( abs ( a ) + abs ( b ) ) 

# =============================================================================
## check if counter is "bad"
def check ( counter , topic , what , logger ) :
    """check if counter is 'bad'
    """

    mean = float ( counter.mean () ) / scale
    rms  = float ( counter.rms  () ) / scale
    mn   = float ( counter.min  () ) / scale
    mx   = float ( counter.max  () ) / scale

    threshold = 1.e-8 
    if abs ( mean ) > threshold :
        logger.error ( "%s : %s  counter mean: %+.4g" %  ( topic , what , mean ) ) 
    if abs ( rms  ) > threshold :
        logger.error ( "%s : %s  counter rms : %+.4g" %  ( topic , what , rms  ) ) 
    if abs ( mn   ) > threshold :
        logger.error ( "%s : %s  counter min : %+.4g" %  ( topic , what , mn   ) ) 
    if abs ( mx   ) > threshold :
        logger.error ( "%s : %s  counter max : %+.4g" %  ( topic , what , mx   ) ) 

# =============================================================================
## numercal 2D-integral 
def int2 ( f2 ,
           xmn = None , xmx = None ,
           ymn = None , ymx = None ) :
    """numercal 2D-integral"""
    if xmn is None : xmn = f2.xmin()
    if xmx is None : xmx = f2.xmax()
    if ymn is None : ymn = f2.ymin()
    if ymx is None : ymx = f2.ymax()
    fi = lambda x,y : f2 ( x , y )
    return I.integral2 ( fi , xmn , xmx , ymn, ymx )

# =============================================================================
## numerical partial 1-integral for 2D-funcntion
def int2_x ( f2 , y , xmn = None , xmx = None ) :
    """numerical partial 1-integral for 2D-funcntion"""
    if xmn is None : xmn = f2.xmin()
    if xmx is None : xmx = f2.xmax()
    fi = lambda x : f2 ( x , y )
    return I.integral ( fi , xmn , xmx )

# =============================================================================
## numerical partial 1-integral for 2D-funcntion
def int2_y ( f2 , x , ymn = None , ymx = None ) :
    """numerical partial 1-integral for 2D-funcntion"""
    if ymn is None : ymn = f2.ymin()
    if ymx is None : ymx = f2.ymax()
    fi = lambda y : f2 ( x , y )
    return I.integral ( fi , ymn , ymx )

# =============================================================================
## numercal 3D-integral 
def int3 ( f3 ,
           xmn = None , xmx = None ,
           ymn = None , ymx = None ,
           zmn = None , zmx = None ) :
    """numercal 3D-integral"""
    if xmn is None : xmn = f3.xmin()
    if xmx is None : xmx = f3.xmax()
    if ymn is None : ymn = f3.ymin()
    if ymx is None : ymx = f3.ymax()
    if zmn is None : zmn = f3.ymin()
    if zmx is None : zmx = f3.ymax()
    fi = lambda x,y,z : f3 ( x , y , z )
    return I.integral3 ( fi , xmn , xmx , ymn, ymx , zmn, zmx )

# =============================================================================
## numerical partial 1D-integral for 3D-funcntion
def int3_x ( f3 , y , z , xmn = None , xmx = None ) :
    """numerical partial 1D-integral for 3D-funcntion"""
    if xmn is None : xmn = f3.xmin()
    if xmx is None : xmx = f3.xmax()
    fi = lambda x : f3 ( x , y , z )
    return I.integral ( fi , xmn , xmx )

# =============================================================================
## numerical partial 1D-integral for 3D-funcntion
def int3_y ( f3 , x , z , ymn = None , ymx = None ) :
    """numerical partial 1D-integral for 3D-funcntion"""
    if ymn is None : ymn = f3.ymin()
    if ymx is None : ymx = f3.ymax()
    fi = lambda y : f3 ( x , y )
    return I.integral ( fi , ymn , ymx )

# =============================================================================
## numerical partial 1D-integral for 3D-funcntion
def int3_z ( f3 , x , y , zmn = None , zmx = None ) :
    """numerical partial 1D-integral for 3D-funcntion"""
    if zmn is None : zmn = f3.zmin()
    if zmx is None : zmx = f3.zmax()
    fi = lambda z : f3 ( x , y , z )
    return I.integral ( fi , zmn , zmx )

# =============================================================================
## numerical partial 2D-integral for 3D-funcntion
def int3_xy ( f3 , z , xmn = None , xmx = None , ymn = None , ymx = None ) :
    """numerical partial 2-integral for 2D-funcntion"""
    if xmn is None : xmn = f3.xmin()
    if xmx is None : xmx = f3.xmax()
    if ymn is None : ymn = f3.ymin()
    if ymx is None : ymx = f3.ymax()
    fi = lambda x,y : f3 ( x , y , z )
    return I.integral2 ( fi , xmn , xmx , ymn , ymx )

# =============================================================================
## numerical partial 2D-integral for 3D-funcntion
def int3_xz ( f3 , y , xmn = None , xmx = None , zmn = None , zmx = None ) :
    """numerical partial 2-integral for 2D-funcntion"""
    if xmn is None : xmn = f3.xmin()
    if xmx is None : xmx = f3.xmax()
    if zmn is None : zmn = f3.zmin()
    if zmx is None : zmx = f3.zmax()
    fi = lambda x,z : f3 ( x , y , z )
    return I.integral2 ( fi , xmn , xmx , zmn , zmx )

# =============================================================================
## numerical partial 2D-integral for 3D-funcntion
def int3_yz ( f3 , x , ymn = None , ymx = None , zmn = None , zmx = None ) :
    """numerical partial 2-integral for 2D-funcntion"""
    if ymn is None : ymn = f3.ymin()
    if ymx is None : ymx = f3.ymax()
    if zmn is None : zmn = f3.zmin()
    if zmx is None : zmx = f3.zmax()
    fi = lambda y,z : f3 ( x , y , z )
    return I.integral2 ( fi , ymn , ymx , zmn , zmx )

# =============================================================================
## Test 2D polynomials
#  @see Ostap::Math::LegendreSum2
#  @see Ostap::Math::Bernstein2D
def test_poly2 () :
    """ Test 2D polynomials
    - see Ostap::Math::LegendreSum2
    - see Ostap::Math::Bernstein2D
    """

    logger = getLogger ( 'test_poly2') 
    
    NX , NY , xmin, xmax, ymin, ymax = 4 , 5 , -1 , 3  , -2 , 5
    ##  NX , NY , xmin, xmax, ymin, ymax = 4 , 5 ,  0 , 1  ,  0 , 1

    ## (1) prepare polynomials
    
    l2 = Ostap.Math.LegendreSum2 ( NX , NY , xmin , xmax , ymin , ymax )
    b2 = Ostap.Math.Bernstein2D  ( NX , NY , xmin , xmax , ymin , ymax )

    xc = 0.5 * ( xmax + xmin )
    xd = 0.3 * ( xmax - xmin )
    yc = 0.5 * ( ymax + ymin )
    yd = 0.3 * ( ymax - ymin )
    for i in range ( 1000 ) :
        a = random.gauss ( 0 , xd )
        b = random.gauss ( 0 , yd )
        x = xc +       a + 1.5  * b
        y = yc + 1.5 * a -        b
        l2.fill ( x , y ) 
        b2.fill ( x , y ) 

    rows  = [ ( '' , '' ,
                'mean [%.0g]' % ( 1/scale ) ,
                'rms  [%.0g]' % ( 1/scale ) , 
                'min  [%.0g]' % ( 1/scale ) , 
                'max  [%.0g]' % ( 1/scale ) ) ]
    
    ## (2) test their equality

    cnt = SE()
    for i in progress_bar ( 100000 ) :
        
        x = random.uniform ( xmin , xmax )
        y = random.uniform ( ymin , ymax )
        
        v1 = l2 ( x , y )
        v2 = b2 ( x , y )

        cnt += diff ( v1 , v2 ) * scale 

    check ( cnt , 'Legendre vs Bernstein' , 'values' , logger )    
    row = 'Legendre  vs Bernstein' , 'values' , \
          '%+.3g'  % float ( cnt.mean() ) , \
          '%.3g'   % float ( cnt.rms () ) , \
          '%+.3g'  % float ( cnt.min () ) , \
          '%+.3g'  % float ( cnt.max () ) 
    rows.append ( row )

    ##  (3) test 2D-integrals
    cnt1, cnt2 , cnt3 = SE(), SE() , SE()  
    for i in progress_bar  ( 1000 ) :
        
        x1 = random.uniform ( xmin , xmax )
        x2 = random.uniform ( xmin , xmax )
        y1 = random.uniform ( ymin , ymax )
        y2 = random.uniform ( ymin , ymax )
        
        v1 = l2.integral ( x1 , x2 , y1 , y2 )
        v2 = b2.integral ( x1 , x2 , y1 , y2 )
        
        v3 = int2 ( l2  , x1  , x2  , y1 , y2 ) 
        v4 = int2 ( b2  , x1  , x2  , y1 , y2 ) 
        
        cnt1 += diff ( v1 , v2 ) * scale 
        cnt2 += diff ( v1 , v3 ) * scale 
        cnt3 += diff ( v2 , v4 ) * scale 
    
    check ( cnt1 , 'Legendre  vs Bernstein' , '2D integrals' , logger )    
    row = 'Legendre  vs Bernstein' , '2D integrals' , \
          '%+.3g'  % float ( cnt1.mean() ) , \
          '%.3g'   % float ( cnt1.rms () ) , \
          '%+.3g'  % float ( cnt1.min () ) , \
          '%+.3g'  % float ( cnt1.max () )
    rows.append ( row )

    check ( cnt2 , 'Legendre  vs numeric' , '2D integrals' , logger )    
    row = 'Legendre  vs numeric' , '2D integrals' , \
          '%+.3g'  % float ( cnt2.mean() ) , \
          '%.3g'   % float ( cnt2.rms () ) , \
          '%+.3g'  % float ( cnt2.min () ) , \
          '%+.3g'  % float ( cnt2.max () )
    rows.append ( row )

    check ( cnt3 , 'Bernstein vs numeric' , '2D integrals' , logger )    
    row = 'Bernstein vs numeric' , '2D integrals' , \
          '%+.3g'  % float ( cnt3.mean() ) , \
          '%.3g'   % float ( cnt3.rms () ) , \
          '%+.3g'  % float ( cnt3.min () ) , \
          '%+.3g'  % float ( cnt3.max () )
    rows.append ( row )


    ## (4) test X-integrals 
    cnt1, cnt2, cnt3, cnt4, cnt5, cnt6  = SE(), SE(), SE(), SE(), SE(), SE() 
    for i in progress_bar  ( 100 ) :
        
        x1    = random.uniform ( xmin , xmax )
        x2    = random.uniform ( xmin , xmax )
        
        y     = random.uniform ( ymin , ymax )

        vl    = l2.integralX  ( x1, x2 ) ( y )
        vb1   = b2.integralX  ( x1, x2 ) ( y ) 
        vb2   = b2.integrateX ( y , x1 , x2 ) 
        
        nl    = int2_x ( l2 , y , x1 , x2 )
        nb    = int2_x ( b2 , y , x1 , x2 )
        
        cnt1 += diff ( vl  , nl  ) * scale 
        cnt2 += diff ( vb1 , nb  ) * scale 
        cnt3 += diff ( vl  , vb1 ) * scale 
        cnt4 += diff ( vb1 , vb2 ) * scale 

        lvv1  = l2.integralX()  ( y ) 
        lvv2  = l2.integralX( l2.xmin() , l2.xmax() ) ( y )
        
        bvv1  = b2.integralX()  ( y ) 
        bvv2  = b2.integralX( b2.xmin() , b2.xmax() ) ( y )
        
        cnt5 += diff ( lvv1 , lvv2 ) * scale 
        cnt6 += diff ( bvv1 , bvv2 ) * scale  

    check ( cnt1 , 'Legendre vs numeric' , 'integralX' , logger )    
    row = 'Legendre  vs numeric' , 'integralX' , \
          '%+.3g'  % float ( cnt1.mean() ) , \
          '%.3g'   % float ( cnt1.rms () ) , \
          '%+.3g'  % float ( cnt1.min () ) , \
          '%+.3g'  % float ( cnt1.max () )
    rows.append ( row )
    
    check ( cnt2 , 'Bernstein vs numeric' , 'integralX' , logger )    
    row = 'Bernstein vs numeric' , 'integralX' , \
          '%+.3g'  % float ( cnt2.mean() ) , \
          '%.3g'   % float ( cnt2.rms () ) , \
          '%+.3g'  % float ( cnt2.min () ) , \
          '%+.3g'  % float ( cnt2.max () )
    rows.append ( row )

    check ( cnt3 , 'Legendre vs Bernstein' , 'integralX' , logger )    
    row = 'Legendre  vs Bernstein' , 'integralX' , \
          '%+.3g'  % float ( cnt3.mean() ) , \
          '%.3g'   % float ( cnt3.rms () ) , \
          '%+.3g'  % float ( cnt3.min () ) , \
          '%+.3g'  % float ( cnt3.max () )
    rows.append ( row )

    check ( cnt4 , 'Bernstein vs Bernstein' , 'integralX/integrateX' , logger )    
    row = 'Bernstein vs Bernstein' , 'integralX/integrateX' , \
          '%+.3g'  % float ( cnt4.mean() ) , \
          '%.3g'   % float ( cnt4.rms () ) , \
          '%+.3g'  % float ( cnt4.min () ) , \
          '%+.3g'  % float ( cnt4.max () )
    rows.append ( row )

    check ( cnt5 , 'Legendre  vs Legendre' , 'integralX/integralX' , logger )    
    row = 'Legendre  vs Legendre' , 'integralX/integralX' , \
          '%+.3g'  % float ( cnt5.mean() ) , \
          '%.3g'   % float ( cnt5.rms () ) , \
          '%+.3g'  % float ( cnt5.min () ) , \
          '%+.3g'  % float ( cnt5.max () )
    rows.append ( row )

    check ( cnt6 , 'Bernstein vs Bernstein' , 'integralX/integralX' , logger )    
    row = 'Bernstein vs Bernstein' , 'integralX/integralX' , \
          '%+.3g'  % float ( cnt6.mean() ) , \
          '%.3g'   % float ( cnt6.rms () ) , \
          '%+.3g'  % float ( cnt6.min () ) , \
          '%+.3g'  % float ( cnt6.max () )
    rows.append ( row )

    ## (5) test Y-integrals 
    cnt1, cnt2, cnt3, cnt4, cnt5, cnt6  = SE(), SE(), SE(), SE(), SE(), SE() 

    for i in progress_bar  ( 100 ) :
        
        x     = random.uniform ( xmin , xmax )

        y1    = random.uniform ( ymin , ymax )
        y2    = random.uniform ( ymin , ymax )

        vl    = l2.integralY  ( y1, y2 ) ( x )
        vb1   = b2.integralY  ( y1, y2 ) ( x ) 
        vb2   = b2.integrateY ( x , y1 , y2 ) 

        
        nl    = int2_y ( l2 , x , y1 , y2 )
        nb    = int2_y ( b2 , x , y1 , y2 )
        
        cnt1 += diff ( vl  , nl  ) * scale 
        cnt2 += diff ( vb1 , nb  ) * scale 
        cnt3 += diff ( vl  , vb1 ) * scale 
        cnt4 += diff ( vb1 , vb2 ) * scale
        
        lvv1  = l2.integralY()  ( x ) 
        lvv2  = l2.integralY( l2.ymin() , l2.ymax() ) ( x )
        
        bvv1  = b2.integralY()  ( x ) 
        bvv2  = b2.integralY( b2.ymin() , b2.ymax() ) ( x )
        
        cnt5 += diff ( lvv1 , lvv2 ) * scale 
        cnt6 += diff ( bvv1 , bvv2 ) * scale  

    check ( cnt1 , 'Legendre vs numeric' , 'integralY' , logger )    
    row = 'Legendre  vs numeric' , 'integralY' , \
          '%+.3g'  % float ( cnt1.mean() ) , \
          '%.3g'   % float ( cnt1.rms () ) , \
          '%+.3g'  % float ( cnt1.min () ) , \
          '%+.3g'  % float ( cnt1.max () )
    rows.append ( row )
    
    check ( cnt2 , 'Bernstein vs numeric' , 'integralY' , logger )    
    row = 'Bernstein vs numeric' , 'integralY' , \
          '%+.3g'  % float ( cnt2.mean() ) , \
          '%.3g'   % float ( cnt2.rms () ) , \
          '%+.3g'  % float ( cnt2.min () ) , \
          '%+.3g'  % float ( cnt2.max () )
    rows.append ( row )

    check ( cnt3 , 'Legendre vs Bernstein' , 'integralY' , logger )    
    row = 'Legendre  vs Bernstein' , 'integralY' , \
          '%+.3g'  % float ( cnt3.mean() ) , \
          '%.3g'   % float ( cnt3.rms () ) , \
          '%+.3g'  % float ( cnt3.min () ) , \
          '%+.3g'  % float ( cnt3.max () )
    rows.append ( row )

    check ( cnt4 , 'Bernstein vs Bernstein' , 'integralY/integrateY' , logger )    
    row = 'Bernstein vs Bernstein' , 'integralY/integrateY' , \
          '%+.3g'  % float ( cnt4.mean() ) , \
          '%.3g'   % float ( cnt4.rms () ) , \
          '%+.3g'  % float ( cnt4.min () ) , \
          '%+.3g'  % float ( cnt4.max () )
    rows.append ( row )

    check ( cnt5 , 'Legendre  vs Legendre' , 'integralY/integralY' , logger )    
    row = 'Legendre  vs Legendre' , 'integralY/integralY' , \
          '%+.3g'  % float ( cnt5.mean() ) , \
          '%.3g'   % float ( cnt5.rms () ) , \
          '%+.3g'  % float ( cnt5.min () ) , \
          '%+.3g'  % float ( cnt5.max () )
    rows.append ( row )

    check ( cnt6 , 'Bernstein vs Bernstein' , 'integralY/integralY' , logger )    
    row = 'Bernstein vs Bernstein' , 'integralY/integralY' , \
          '%+.3g'  % float ( cnt6.mean() ) , \
          '%.3g'   % float ( cnt6.rms () ) , \
          '%+.3g'  % float ( cnt6.min () ) , \
          '%+.3g'  % float ( cnt6.max () )
    rows.append ( row )
    
    
    title = '2D polynomials'
    table = T.table ( rows , title = title , prefix = '# ' , alignment ='llcccc' )
    logger.info ( '%s\n%s' %  ( title , table ) )

# =============================================================================
## Test 3D polynomials
#  @see Ostap::Math::LegendreSum3
#  @see Ostap::Math::Bernstein3D
def test_poly3 () :
    """ Test 3D polynomials
    - see Ostap::Math::LegendreSum3
    - see Ostap::Math::Bernstein3D
    """

    logger = getLogger ( 'test_poly3') 

    NX , NY , NZ   , \
       xmin , xmax , \
       ymin , ymax , \
       zmin , zmax = 3 , 2 , 4 , -1 , 3  , -2 , 5 , 0 , 3 
       
    ##  NX , NY , xmin, xmax, ymin, ymax = 4 , 5 ,  0 , 1  ,  0 , 1

    ## (1) prepare polynomials
    
    l3 = Ostap.Math.LegendreSum3 ( NX , NY , NZ , xmin , xmax , ymin , ymax , zmin, zmax )
    b3 = Ostap.Math.Bernstein3D  ( NX , NY , NZ , xmin , xmax , ymin , ymax , zmin, zmax )

    xc = 0.5 * ( xmax + xmin )
    xd = 0.3 * ( xmax - xmin )
    yc = 0.5 * ( ymax + ymin )
    yd = 0.3 * ( ymax - ymin )
    zc = 0.5 * ( zmax + zmin )
    zd = 0.5 * ( zmax - zmin )
    
    for i in range ( 1000 ) :
        a = random.gauss ( 0 , xd )
        b = random.gauss ( 0 , yd )
        c = random.gauss ( 0 , zd )
        x = xc +       a + 1.5  * b + 0.5 * c 
        y = yc + 1.5 * a -        b - 0.5 * c
        z = zc + 1.0 * a -        b +       c
        
        l3.fill ( x , y , z ) 
        b3.fill ( x , y , z ) 

    rows  = [ ( '' , '' ,
                'mean [%.0g]' % ( 1/scale ) ,
                'rms  [%.0g]' % ( 1/scale ) , 
                'min  [%.0g]' % ( 1/scale ) , 
                'max  [%.0g]' % ( 1/scale ) ) ]
    
    ## (2) test their equality

    cnt = SE()
    for i in progress_bar ( 10000 ) :
        
        x = random.uniform ( xmin , xmax )
        y = random.uniform ( ymin , ymax )
        z = random.uniform ( zmin , zmax )
        
        v1 = l3 ( x , y , z )
        v2 = b3 ( x , y , z )

        cnt += diff ( v1 , v2 ) * scale 

    check ( cnt , 'Legendre vs Bernstein' , 'values' , logger )    
    row = 'Legendre  vs Bernstein' , 'values' , \
          '%+.3g'  % float ( cnt.mean() ) , \
          '%.3g'   % float ( cnt.rms () ) , \
          '%+.3g'  % float ( cnt.min () ) , \
          '%+.3g'  % float ( cnt.max () ) 
    rows.append ( row )


    ##  (3) test 3D-integrals
    cnt1, cnt2 , cnt3 = SE(), SE() , SE()  
    for i in progress_bar  ( 100 ) :
        
        x1 = random.uniform ( xmin , xmax )
        x2 = random.uniform ( xmin , xmax )
        y1 = random.uniform ( ymin , ymax )
        y2 = random.uniform ( ymin , ymax )
        z1 = random.uniform ( zmin , zmax )
        z2 = random.uniform ( zmin , zmax )
        
        v1 = l3.integral ( x1 , x2 , y1 , y2 , z1, z2)
        v2 = b3.integral ( x1 , x2 , y1 , y2 , z1, z2)
        
        v3 = int3 ( l3  , x1  , x2  , y1 , y2 , z1 , z2  ) 
        v4 = int3 ( b3  , x1  , x2  , y1 , y2 , z1 , z2 ) 
        
        cnt1 += diff ( v1 , v2 ) * scale 
        cnt2 += diff ( v1 , v3 ) * scale 
        cnt3 += diff ( v2 , v4 ) * scale 

    check ( cnt1 , 'Legendre  vs Bernstein' , '3D integrals' , logger )    
    row = 'Legendre  vs Bernstein' , '3D integrals' , \
          '%+.3g'  % float ( cnt1.mean() ) , \
          '%.3g'   % float ( cnt1.rms () ) , \
          '%+.3g'  % float ( cnt1.min () ) , \
          '%+.3g'  % float ( cnt1.max () )
    rows.append ( row )

    check ( cnt2 , 'Legendre  vs numeric' , '3D integrals' , logger )    
    row = 'Legendre  vs numeric' , '3D integrals' , \
          '%+.3g'  % float ( cnt2.mean() ) , \
          '%.3g'   % float ( cnt2.rms () ) , \
          '%+.3g'  % float ( cnt2.min () ) , \
          '%+.3g'  % float ( cnt2.max () )
    rows.append ( row )

    check ( cnt3 , 'Bernstein vs numeric' , '2D integrals' , logger )    
    row = 'Bernstein vs numeric' , '3D integrals' , \
          '%+.3g'  % float ( cnt3.mean() ) , \
          '%.3g'   % float ( cnt3.rms () ) , \
          '%+.3g'  % float ( cnt3.min () ) , \
          '%+.3g'  % float ( cnt3.max () )
    rows.append ( row )

    ## test 2D-integrals
    
    bxy = b3.integralXY ()
    bxz = b3.integralXZ ()
    byz = b3.integralYZ ()
    
    cnt1, cnt2 , cnt3 = SE(), SE() , SE()
    for i in progress_bar  ( 100 ) :
        
        x = random.uniform ( xmin , xmax )
        y = random.uniform ( ymin , ymax )
        z = random.uniform ( zmin , zmax )

        vx = byz ( x )
        vy = bxz ( y )
        vz = bxy ( z )

        nx = int3_yz ( b3  , x )
        ny = int3_xz ( b3  , y )
        nz = int3_xy ( b3  , z )
        
        cnt1 += diff ( vx , nx ) * scale 
        cnt2 += diff ( vy , ny ) * scale 
        cnt3 += diff ( vz , nz ) * scale 


    check ( cnt1 , 'Bernstein vs numeric' , 'integralYZ' , logger )    
    row = 'Bernstein vs numeric' , 'integralYZ' , \
          '%+.3g'  % float ( cnt1.mean() ) , \
          '%.3g'   % float ( cnt1.rms () ) , \
          '%+.3g'  % float ( cnt1.min () ) , \
          '%+.3g'  % float ( cnt1.max () )
    rows.append ( row )

    check ( cnt2 , 'Bernstein vs numeric' , 'integralXZ' , logger )    
    row = 'Bernstein vs numeric' , 'integralXZ' , \
          '%+.3g'  % float ( cnt2.mean() ) , \
          '%.3g'   % float ( cnt2.rms () ) , \
          '%+.3g'  % float ( cnt2.min () ) , \
          '%+.3g'  % float ( cnt2.max () )
    rows.append ( row )

    check ( cnt3 , 'Bernstein vs numeric' , 'integralXY' , logger )    
    row = 'Bernstein vs numeric' , 'integralXY' , \
          '%+.3g'  % float ( cnt3.mean() ) , \
          '%.3g'   % float ( cnt3.rms () ) , \
          '%+.3g'  % float ( cnt3.min () ) , \
          '%+.3g'  % float ( cnt3.max () )
    rows.append ( row )


    ## test 2D-integrals
    
    cnt1, cnt2 , cnt3 = SE(), SE() , SE()
    for i in progress_bar  ( 100 ) :
        
        x = random.uniform ( xmin , xmax )
        y = random.uniform ( ymin , ymax )
        z = random.uniform ( zmin , zmax )

        x1 = random.uniform ( xmin , xmax )
        x2 = random.uniform ( xmin , xmax )
        y1 = random.uniform ( ymin , ymax )
        y2 = random.uniform ( ymin , ymax )
        z1 = random.uniform ( zmin , zmax )
        z2 = random.uniform ( zmin , zmax )

        bxy = b3.integralXY ( x1 , x2 , y1 , y2 )
        bxz = b3.integralXZ ( x1 , x2 , z1 , z2 )
        byz = b3.integralYZ ( y1 , y2 , z1 , z2 )
        
        vx = byz ( x )
        vy = bxz ( y )
        vz = bxy ( z )

        nx = int3_yz ( b3  , x , y1 , y2 , z1 , z2 )
        ny = int3_xz ( b3  , y , x1 , x2 , z1 , z2 )
        nz = int3_xy ( b3  , z , x1 , x2 , y1 , y2 )
        
        cnt1 += diff ( vx , nx ) * scale 
        cnt2 += diff ( vy , ny ) * scale 
        cnt3 += diff ( vz , nz ) * scale 


    check ( cnt1 , 'Bernstein vs numeric' , 'integralYZ/2' , logger )    
    row = 'Bernstein vs numeric' , 'integralYZ/2' , \
          '%+.3g'  % float ( cnt1.mean() ) , \
          '%.3g'   % float ( cnt1.rms () ) , \
          '%+.3g'  % float ( cnt1.min () ) , \
          '%+.3g'  % float ( cnt1.max () )
    rows.append ( row )

    check ( cnt2 , 'Bernstein vs numeric' , 'integralXZ/2' , logger )    
    row = 'Bernstein vs numeric' , 'integralXZ/2' , \
          '%+.3g'  % float ( cnt2.mean() ) , \
          '%.3g'   % float ( cnt2.rms () ) , \
          '%+.3g'  % float ( cnt2.min () ) , \
          '%+.3g'  % float ( cnt2.max () )
    rows.append ( row )

    check ( cnt3 , 'Bernstein vs numeric' , 'integralXY/2' , logger )    
    row = 'Bernstein vs numeric' , 'integralXY/2' , \
          '%+.3g'  % float ( cnt3.mean() ) , \
          '%.3g'   % float ( cnt3.rms () ) , \
          '%+.3g'  % float ( cnt3.min () ) , \
          '%+.3g'  % float ( cnt3.max () )
    rows.append ( row )


    ## test 2D-integrals
    
    cnt1, cnt2 , cnt3 = SE(), SE() , SE()
    for i in progress_bar  ( 100 ) :
        
        x = random.uniform ( xmin , xmax )
        y = random.uniform ( ymin , ymax )
        z = random.uniform ( zmin , zmax )

        x1 = random.uniform ( xmin , xmax )
        x2 = random.uniform ( xmin , xmax )
        y1 = random.uniform ( ymin , ymax )
        y2 = random.uniform ( ymin , ymax )
        z1 = random.uniform ( zmin , zmax )
        z2 = random.uniform ( zmin , zmax )

        bxy = b3.integralXY ( x1 , x2 , y1 , y2 )
        bxz = b3.integralXZ ( x1 , x2 , z1 , z2 )
        byz = b3.integralYZ ( y1 , y2 , z1 , z2 )
        
        vx = b3.integrateYZ ( x  , y1 , y2 , z1 , z2 )
        vy = b3.integrateXZ ( y  , x1 , x2 , z1 , z2 )
        vz = b3.integrateXY ( z  , x1 , x2 , y1 , y2 )

        nx = int3_yz ( b3  , x , y1 , y2 , z1 , z2 )
        ny = int3_xz ( b3  , y , x1 , x2 , z1 , z2 )
        nz = int3_xy ( b3  , z , x1 , x2 , y1 , y2 )
        
        cnt1 += diff ( vx , nx ) * scale 
        cnt2 += diff ( vy , ny ) * scale 
        cnt3 += diff ( vz , nz ) * scale 


    check ( cnt1 , 'Bernstein vs numeric' , 'integrateYZ' , logger )    
    row = 'Bernstein vs numeric' , 'integrateYZ' , \
          '%+.3g'  % float ( cnt1.mean() ) , \
          '%.3g'   % float ( cnt1.rms () ) , \
          '%+.3g'  % float ( cnt1.min () ) , \
          '%+.3g'  % float ( cnt1.max () )
    rows.append ( row )

    check ( cnt2 , 'Bernstein vs numeric' , 'integrateXZ' , logger )    
    row = 'Bernstein vs numeric' , 'integrateXZ' , \
          '%+.3g'  % float ( cnt2.mean() ) , \
          '%.3g'   % float ( cnt2.rms () ) , \
          '%+.3g'  % float ( cnt2.min () ) , \
          '%+.3g'  % float ( cnt2.max () )
    rows.append ( row )

    check ( cnt3 , 'Bernstein vs numeric' , 'integrateXY' , logger )    
    row = 'Bernstein vs numeric' , 'integrateXY' , \
          '%+.3g'  % float ( cnt3.mean() ) , \
          '%.3g'   % float ( cnt3.rms () ) , \
          '%+.3g'  % float ( cnt3.min () ) , \
          '%+.3g'  % float ( cnt3.max () )
    rows.append ( row )


    title = '3D polynomials'
    table = T.table ( rows , title = title , prefix = '# ' , alignment ='llcccc' )
    logger.info ( '%s\n%s' %  ( title , table ) )


# =============================================================================
if '__main__' == __name__ :
    
    ## test_poly2 ()
    test_poly3 ()

# =============================================================================
##                                                                      The END  
# ============================================================================= 
