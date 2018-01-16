#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
""" Test module for 2D models
"""
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_bernstein' ) 
else                       : logger = getLogger ( __name__         )
# ============================================================================= 
import ROOT, random  
import ostap.math.models 
from   ostap.core.core     import Ostap, SE
from   ostap.math.integral import integral2, Integrate2D_X, Integrate2D_Y

# ============================================================================
def test_models ():
    
    xmin = 0
    xmax = 2
    ymin = xmin
    ymax = xmax
    
    psx = Ostap.Math.PhaseSpaceNL( xmin + 0.00001 * ( xmax - xmin ) ,
                                   xmax - 0.00001 * ( xmax - xmin ) , 3 , 7 )
    psy = Ostap.Math.PhaseSpaceNL( ymin + 0.00001 * ( ymax - ymin ) ,
                                   ymax - 0.00001 * ( ymax - ymin ) , 3 , 7 )
    
    funcs = [
        Ostap.Math.Bernstein2D        ( 2  , 2 , xmin , xmax  , ymin , ymax ) ,
        Ostap.Math.Bernstein2DSym     ( 2  , xmin , xmax ) ,
        Ostap.Math.Positive2D         ( 2  , 2 , xmin , xmax , ymin , ymax ) ,
        Ostap.Math.Positive2DSym      ( 2  , xmin , xmax ) ,
        Ostap.Math.PS2DPol            ( psx  , psy , 2 , 2 , xmin , xmax , ymin , ymax ) ,
        Ostap.Math.PS2DPolSym         ( psx  , 2 , xmin , xmax ) ,
        Ostap.Math.ExpoPS2DPol        ( psx  , xmin , xmax , 2 , 2 , ymin , ymax ) ,
        Ostap.Math.Expo2DPol          ( xmin , xmax , ymin , ymax , 2  , 2 ) ,
        Ostap.Math.Expo2DPolSym       ( xmin , xmax , 2 ) ,
        ]

    cnt1 = SE()  
    cnt2 = SE()  
    cnt3 = SE()  
    for f in funcs :
        ## print f.xmin() , f.xmax() , f.ymin() , f.ymax(), type(f)
        for i in range(f.npars() ) :
            f.setPar( i, random.uniform ( 1 , 5 ) )
            if hasattr  ( f , 'setTau'  ) : f.setTau  ( random.uniform ( -2  , 2  ) )
            if hasattr  ( f , 'setTauX' ) : f.setTauX ( random.uniform ( -2  , 2  ) )
            if hasattr  ( f , 'setTauY' ) : f.setTauY ( random.uniform ( -2  , 2  ) )
            
    for i in range(0,1000) :
        
        if 0 == i:
            x1,x2,y1,y2 = xmin,xmax,ymin,ymax
        else :
            x1 = random.uniform ( xmin , xmax/2  )
            x2 = random.uniform ( x1   , xmax    )
            y1 = random.uniform ( ymin , ymax/2  )
            y2 = random.uniform ( y1   , ymax    )
            
        for f in funcs :
            
            i1 = f.integral (       x1 , x2 , y1 , y2 )
            i2 = integral2  ( f   , x1 , x2 , y1 , y2 )
            r1 = (i1-i2)/(abs(i1)+abs(i2))
            assert  abs(r1) < 1.e-5 , 'I2:ERROR: difference is too large: %s (%.2f,%.2f,%.2f,%.2f) %s' % ( r1 , x1 , x2  , y1 , y2 , type(f) )  
            cnt1 += r1

            ym = 0.5 * (  y1 + y2 )
            
            i1 = f.integrateX ( ym , x1 , x2 )
            
            IX = Integrate2D_X ( f , x1 , x2 )
            i2 = IX( ym )

            r2 = (i1-i2)/(abs(i1)+abs(i2))
            assert  abs(r2) < 1.e-5 , 'IX:ERROR: difference is too large: %s (%.2f,%.2f,%.2f) %s' % ( r1 , x1 , x2  , ym , type(f) )  
            cnt2 += r2

            xm = 0.5 * (  x1 + x2 )
            
            i1 = f.integrateY ( xm , y1 , y2 )
            
            IY = Integrate2D_Y ( f , y1 , y2 )
            i2 = IY( xm )

            r3 = (i1-i2)/(abs(i1)+abs(i2))
            assert  abs(r3) < 1.e-5 , 'IY:ERROR: difference is too large: %s (%.2f,%.2f,%.2f) %s' % ( r1 , y1 , y2  , xm , type(f) )  
            cnt3 += r3


            
            
            
    print 'COUNTER(I2):' , cnt1
    print 'COUNTER(IX):' , cnt2
    print 'COUNTER(IY):' , cnt3
            
# =============================================================================
if '__main__' == __name__ :
        
    test_models ()

    
# =============================================================================
# The END 
# =============================================================================
