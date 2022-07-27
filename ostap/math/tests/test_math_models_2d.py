#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_models_2d.py
#  Test module for the file ostap/math/models_2d.py
# ============================================================================= 
""" Test module for 2D-models
"""
# ============================================================================= 
from   __future__          import print_function
import ostap.math.models 
from   ostap.core.core     import Ostap, SE
from   ostap.math.integral import integral2, Integrate2D_X, Integrate2D_Y
import ROOT, random  
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_math_models_2d' ) 
else                       : logger = getLogger ( __name__              )
# ============================================================================= 

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
        for i in range ( f.npars () ) :
            f.setPar( i, random.uniform ( 1 , 5 ) )
            if hasattr  ( f , 'setTau'  ) : f.setTau  ( random.uniform ( -2  , 2  ) )
            if hasattr  ( f , 'setTauX' ) : f.setTauX ( random.uniform ( -2  , 2  ) )
            if hasattr  ( f , 'setTauY' ) : f.setTauY ( random.uniform ( -2  , 2  ) )
            
    for i in range(0,100) :
        
        if 0 == i:
            
            x1, x2, y1, y2 = xmin, xmax, ymin, ymax
            
        else :
            
            x1 = random.uniform ( xmin , xmax / 2  )
            x2 = random.uniform ( x1   , xmax      )
            y1 = random.uniform ( ymin , ymax / 2  )
            y2 = random.uniform ( y1   , ymax      )
            
        for f in funcs :
            
            i1 = f.integral (       x1 , x2 , y1 , y2 )
            i2 = integral2  ( f   , x1 , x2 , y1 , y2 )
            
            r1 = ( i1 - i2 ) / ( abs ( i1 ) + abs ( i2 ) )
            if 1.e-5 < abs ( r1 ) : logger.error ( 'I2:ERROR: difference is too large: %.6g (%.4g,%.4g,%.4g,%.4g) %s' % ( r1 , x1 , x2 , y1 , y2 , type(f) )  )
            cnt1 += r1
            
            ym = 0.5 * (  y1 + y2 )            
            i1 = f.integrateX  ( ym , x1 , x2 )            
            IX = Integrate2D_X ( f  , x1 , x2 )
            i2 = IX( ym )

            r2 = ( i1 - i2 ) / ( abs ( i1 ) + abs ( i2 ) )
            if 1.e-5 < abs ( r2 ) : logger.error ( 'IX:ERROR: difference is too large: %.6g (%.4g,%.4g,%.4g,%.4g) %s' % ( r2 , x1 , x2 , y1 , y2  , type(f) )  )
            cnt2 += r2

            xm = 0.5 * (  x1 + x2 )
            
            i1 = f.integrateY  ( xm , y1 , y2 )            
            IY = Integrate2D_Y ( f  , y1 , y2 )
            i2 = IY ( xm )

            r3 = ( i1 - i2 ) / ( abs ( i1 ) + abs ( i2 ) )
            if 1.e-5 < abs ( r3 ) : logger.error ( 'IY:ERROR: difference is too large: %.6g (%.4g,%.4g,%.4g,%.4g) %s' % ( r3 , x1 , x2 , y1 , y2  , type(f) )  )
            cnt3 += r3

    logger.info ( 'Counter(I2) %s' % cnt1 )
    logger.info ( 'Counter(IX) %s' % cnt2 )
    logger.info ( 'Counter(IY) %s' % cnt3 )
    
# =============================================================================
if '__main__' == __name__ :
        
    test_models ()

    
# =============================================================================
##                                                                      The END 
# =============================================================================
