#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_models_3d.py
#  Test module for the file ostap/math/models.py
# ============================================================================= 
""" Test module for 3D models
"""
# ============================================================================= 
from   __future__          import print_function
import ostap.math.models 
from   ostap.core.core     import Ostap, SE
from   ostap.math.integral import integral3
from   ostap.math.integral import Integrate3D_X , Integrate3D_Y  , Integrate3D_Z
from   ostap.math.integral import Integrate3D_XY, Integrate3D_XZ , Integrate3D_YZ
import ROOT, random  
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_math_models_3d' ) 
else                       : logger = getLogger ( __name__         )
# ============================================================================= 


# ============================================================================
def test_models ():
    
    xmin = 0
    xmax = 2
    ymin = xmin
    ymax = xmax
    zmin = xmin
    zmax = xmax
    
    funcs = [
        Ostap.Math.Bernstein3D    ( 2  , 2 , 2   ,
                                    xmin , xmax  ,
                                    ymin , ymax  ,
                                    zmin , zmax  ) ,
        Ostap.Math.Bernstein3DSym ( 2  , xmin , xmax ) ,
        Ostap.Math.Bernstein3DMix ( 2  , 2 , xmin , xmax , zmin , zmax ) ,        
        Ostap.Math.Positive3D     ( 2  , 2 , 2   ,
                                    xmin , xmax  ,
                                    ymin , ymax  ,
                                    zmin , zmax  ) ,
        Ostap.Math.Positive3DSym  ( 2  , xmin , xmax ) ,
        Ostap.Math.Positive3DMix  ( 2  , 2 , xmin , xmax , zmin , zmax ) ,        
        ]
    
    cnt1 = SE()  
    cnt2 = SE()  
    cnt3 = SE()  
    cnt4 = SE()  
    cnt5 = SE()  
    cnt6 = SE()  
    cnt7 = SE()
    
    for f in funcs :
        ## print f.xmin() , f.xmax() , f.ymin() , f.ymax(), type(f)
        for i in range(f.npars() ) :
            f.setPar( i, random.uniform ( 1 , 5 ) )

    for i in range(0,100) :
        
        if 0 == i:
            x1,x2,y1,y2,z1,z2  = xmin,xmax,ymin,ymax,zmin,zmax 
        else :
            x1 = random.uniform ( xmin , xmax/2  )
            x2 = random.uniform ( x1   , xmax    )
            y1 = random.uniform ( ymin , ymax/2  )
            y2 = random.uniform ( y1   , ymax    )
            z1 = random.uniform ( zmin , zmax/2  )
            z2 = random.uniform ( z1   , zmax    )

        for f in funcs :
            
            i1 = f.integral (       x1 , x2 , y1 , y2 , z1 , z2 )
            i2 = integral3  ( f   , x1 , x2 , y1 , y2 , z1 , z2 )

            
            r1 = (i1-i2)/(abs(i1)+abs(i2))
            assert  abs(r1) < 1.e-5 , 'I3:ERROR: difference is too large: %s (%.2f,%.2f,%.2f,%.2f,%.2f,%.2f) %s' % ( r1 , x1 , x2  , y1 , y2 ,  z1 , z2 , type(f) )  
            cnt1 += r1

            xm = 0.5 * ( x1 +  x2 )
            ym = 0.5 * ( y1 +  y2 )
            zm = 0.5 * ( z1 +  z2 )

            i1 = f.integrateX  ( ym , zm , x1 , x2 )
            IX = Integrate3D_X ( f , x1  , x2 )
            i2 = IX ( ym , zm )
            
            r2 = (i1-i2)/(abs(i1)+abs(i2))
            assert  abs(r2) < 1.e-5 , 'IX:ERROR: difference is too large: %s (%.2f,%.2f,%.2f,%.2f) %s' % ( r2 , x1 , x2  , ym ,  zm , type(f) )  
            cnt2 += r2

            i1 = f.integrateY  ( xm , zm , y1 , y2 )
            IY = Integrate3D_Y ( f , y1  , y2 )
            i2 = IY ( xm , zm )
            
            r3 = (i1-i2)/(abs(i1)+abs(i2))
            assert  abs(r3) < 1.e-5 , 'IY:ERROR: difference is too large: %s (%.2f,%.2f,%.2f,%.2f) %s' % ( r2 , xm , y1 , y2 ,  zm , type(f) )  
            cnt3 += r3

            i1 = f.integrateZ  ( xm , ym , z1 , z2 )
            IZ = Integrate3D_Z ( f , z1  , z2 )
            i2 = IZ ( xm , ym )
        
            r4 = (i1-i2)/(abs(i1)+abs(i2))
            assert  abs(r4) < 1.e-5 , 'IY:ERROR: difference is too large: %s (%.2f,%.2f,%.2f,%.2f) %s' % ( r2 , xm , ym , z1 , z2  , type(f) )  
            cnt4 += r4

            i1 = f.integrateXY  ( zm  , x1 , x2 , y1 ,  y2  )
            I5 = Integrate3D_XY ( f , x1 , x2 , y1 , y2 ) 
            i2 = I5 ( zm )
        
            r5 = (i1-i2)/(abs(i1)+abs(i2))
            assert  abs(r5) < 1.e-5 , 'I5:ERROR: difference is too large: %s (%.2f,%.2f,%.2f,%.2f,%.2f) %s' % ( r2 , x1 , x2 , y1 , y2 , zm  , type(f) )  
            cnt5 += r5
            
            i1 = f.integrateXZ  ( ym  , x1 , x2 , z1 ,  z2  )
            I6 = Integrate3D_XZ ( f , x1 , x2 , z1 , z2 ) 
            i2 = I6 ( ym )
        
            r6 = (i1-i2)/(abs(i1)+abs(i2))
            assert  abs(r6) < 1.e-5 , 'I6:ERROR: difference is too large: %s (%.2f,%.2f,%.2f,%.2f,%.2f) %s' % ( r2 , x1 , x2 , ym , z1 , z2  , type(f) )  
            cnt6 += r6

            i1 = f.integrateYZ  ( xm  , y1 , y2 , z1 ,  z2  )
            I7 = Integrate3D_YZ ( f , y1 , y2 , z1 , z2 ) 
            i2 = I7 ( xm )
        
            r7 = (i1-i2)/(abs(i1)+abs(i2))
            assert  abs(r7) < 1.e-5 , 'I7:ERROR: difference is too large: %s (%.2f,%.2f,%.2f,%.2f,%.2f) %s' % ( r2 , xm , y1 , y2 , z1 , z2  , type(f) )  
            cnt7 += r7

            
    print ( 'COUNTER(I3):' , cnt1 )
    print ( 'COUNTER(IX):' , cnt2 ) 
    print ( 'COUNTER(IY):' , cnt3 ) 
    print ( 'COUNTER(IZ):' , cnt4 ) 
    print ( 'COUNTER(I5):' , cnt5 )
    print ( 'COUNTER(I6):' , cnt6 ) 
    print ( 'COUNTER(I7):' , cnt7 )

# =============================================================================
if '__main__' == __name__ :
        
    test_models ()

    
# =============================================================================
# The END 
# ============================================================================
