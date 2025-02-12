#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_poly.py
#  Test script for polynomial functions  
#  @see Ostap::Math::LegendreSum
#  @see Ostap::Math::HermiteSum
#  @see Ostap::Math::ChebyshevSum
#  @see Ostap::Math::Bernstein
#  @see Ostap::Math::Polynomial
# ============================================================================= 
""" Test script for 2&3-dimnsional polynomial functions  
- see Ostap::Math::LegendreSum
- see Ostap::Math::HermiteSum
- see Ostap::Math::ChebyshevSum
- see Ostap::Math::Bernstein 
- see Ostap::Math::Polynomial
"""
# ============================================================================= 
from   ostap.core.core          import Ostap, SE  
from   ostap.logger.colorized   import attention
from   ostap.utils.progress_bar import progress_bar
from   ostap.utils.utils        import batch_env 
import ostap.math.integral      as     I 
import ostap.math.derivative    as     D 
import ostap.logger.table       as     T
import ostap.math.models 
import ROOT, random, math 
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_poly' )
else                       : logger = getLogger ( __name__               )
# ============================================================================= 
## set batch from environment 
batch_env ( logger )
# =============================================================================
#

xmin  = 0
xmax  = 10
N     = 5
scale = 1.e+12

# =============================================================================
## test derivative
def test_poly_derivative () :
    """Test derivative
    """
    logger = getLogger ( 'test_poly_derivative') 

    poly = (
        Ostap.Math.ChebyshevSum ( N , xmin , xmax ) ,
        Ostap.Math.LegendreSum  ( N , xmin , xmax ) ,
        Ostap.Math.HermiteSum   ( N , xmin , xmax ) ,
        Ostap.Math.Bernstein    ( N , xmin , xmax ) ,
        Ostap.Math.Polynomial   ( N , xmin , xmax ) ,
        )


    for p in poly :
        n = p.npars()
        for i in range ( n ) : p [ i ] = random.uniform ( -10 , 10 )

    ders = tuple ( [ p.derivative() for p in poly ] )

    cnts = [ SE() for i in poly ]

    dx = 0.1 * ( xmax - xmin )
    for i in progress_bar ( range ( 1000 ) ) :
        
        x = random.uniform ( xmin + dx  , xmax - dx )
        
        for i, p in enumerate ( poly ) :
            d1 = p.derivative ( x )
            d2 = ders [ i ]   ( x ) 
            d3 = D.derivative ( p , x )

            dd = max ( abs ( d1 - d2 ) / ( abs ( d1 ) + abs ( d2 ) ) ,
                       abs ( d1 - d3 ) / ( abs ( d1 ) + abs ( d3 ) ) ,
                       abs ( d2 - d3 ) / ( abs ( d2 ) + abs ( d3 ) ) )
            
            ## dd = abs ( d2 - d3 ) / ( abs ( d2 ) + abs ( d3 ) )
            

            cnts [ i ] += 2 * dd * scale 


    names = 'Chebyshev' , 'Legendre' , 'Hermite' , 'Bernstein' , 'Polynomial'

    rows = [ ( 'Poly' ,
               'mean [%.0g]' % ( 1/scale ) ,
               'rms  [%.0g]' % ( 1/scale ) ,
               'max  [%.0g]' % ( 1/scale ) ) ]
    
    for name, cnt  in zip ( names , cnts ) :
        
        row = name , \
              '%.3g' % cnt.mean () , \
              '%.3g' % cnt.rms  () , \
              '%.3g' % cnt.max  ()
        
        rows.append ( row )
        
    title = 'Derivatives (analytical&numerical)'
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lccc' )
    logger.info ( '%s\n%s' %  ( title , table ) )


# =============================================================================
## test integrals 
def test_poly_integral () :
    """Test derivative
    """
    logger = getLogger ( 'test_poly_integral') 

    poly = (
        Ostap.Math.ChebyshevSum ( N , xmin , xmax ) ,
        Ostap.Math.LegendreSum  ( N , xmin , xmax ) ,
        Ostap.Math.HermiteSum   ( N , xmin , xmax ) ,
        Ostap.Math.Bernstein    ( N , xmin , xmax ) ,
        Ostap.Math.Polynomial   ( N , xmin , xmax ) ,
        )


    for p in poly :
        n = p.npars()
        for i in range ( n ) : p [ i ] = random.uniform ( -10 , 10 )

    ints = tuple ( [ p.indefinite_integral() for p in poly ] )

    cnts = [ SE() for i in poly ]

    dx = 0.1 * ( xmax - xmin )
    for i in progress_bar ( range ( 1000 ) ) :
        
        x = random.uniform ( xmin , xmax )
        y = random.uniform ( xmin , xmax )
        
        for i, p in enumerate ( poly ) :
            d1 = p.integral ( x , y ) 
            d2 = ints [ i ] ( y ) -  ints [ i ] ( x ) 
            d3 = I.integral ( p , x , y )

            dd = max ( abs ( d1 - d2 ) / ( abs ( d1 ) + abs ( d2 ) ) ,
                       abs ( d1 - d3 ) / ( abs ( d1 ) + abs ( d3 ) ) ,
                       abs ( d2 - d3 ) / ( abs ( d2 ) + abs ( d3 ) ) )
            
            ## dd = abs ( d2 - d3 ) / ( abs ( d2 ) + abs ( d3 ) )
            

            cnts [ i ] += 2 * dd * scale 


    names = 'Chebyshev' , 'Legendre' , 'Hermite' , 'Bernstein' , 'Polynomial'

    rows = [ ( 'Poly' ,
               'mean [%.0g]' % ( 1/scale ) ,
               'rms  [%.0g]' % ( 1/scale ) ,
               'max  [%.0g]' % ( 1/scale ) ) ]
    
    for name, cnt  in zip ( names , cnts ) :
        
        row = name , \
              '%.3g' % cnt.mean () , \
              '%.3g' % cnt.rms  () , \
              '%.3g' % cnt.max  ()
        
        rows.append ( row )
        
    title = 'Integrals (analytical&numerical)'
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lccc' )
    logger.info ( '%s\n%s' %  ( title , table ) )


# =============================================================================

if '__main__' == __name__ :
    
    test_poly_derivative ()
    test_poly_integral   ()
 

# =============================================================================
##                                                                      The END  
# ============================================================================= 
