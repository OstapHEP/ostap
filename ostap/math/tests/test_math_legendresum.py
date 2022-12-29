#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_legendresum.py
#  Test script legendre sums 
#  @see Ostap::Math::LegendreSum
#  @see Ostap::Math::LegendreSum2
#  @see Ostap::Math::LegendreSum3
# ============================================================================= 
""" Test script  foe legendre sums 
- see Ostap::Math::LegendreSum
- see Ostap::Math::LegendreSum2
- see Ostap::Math::LegendreSum3
"""
# ============================================================================= 
from   ostap.core.core        import Ostap
from   ostap.logger.colorized import attention
import ostap.math.integral    as     I 
import ostap.logger.table     as     T
import ROOT, random, math 
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_legendresum' )
else                       : logger = getLogger ( __name__                      )
# ============================================================================= 


from   ostap.core.core import Ostap, VE
import ostap.core.pyrouts

LS1 = Ostap.Math.LegendreSum
LS2 = Ostap.Math.LegendreSum2
LS3 = Ostap.Math.LegendreSum3
LS4 = Ostap.Math.LegendreSum4

delta     = 1.e-10
tolerance = 1.e-6
# =============================================================================
def test_legendresum_1 ( ) :

    logger = getLogger ("test_legendresum_1")
    logger.info ( "Test for 1D splits" )
    
    xmin , xmax = 0 , 2

    rows = [ ( 'i0' ,
               'i1' , 'i2' , 'i3+i4' ,
               'i1-i2'   ,
               'i1/i2-1' ,
               'i1-(i3+i4)'  , 
               'i1/(i3+i4)-1' ) ]
    
    ls = [] 
    for i in range ( 10 ) :
        
        ls = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( i + 1 ) ] , xmin , xmax )
        
        ## whole range integral 
        i1 = ls.integral ()

        ## shoudl be numerically very similar 
        i2 = ls.integral ( xmin + delta , xmax - delta )
        i2 = ls.integral ( xmin + delta , xmax - delta )
        
        split = random.uniform ( xmin , xmax )
        
        i3 = ls.integral ( xmin   , split )
        i4 = ls.integral ( split  , xmax  )

        i0 = I.integral  ( ls , xmin = xmin , xmax = xmax )
        
        d1 =   i1 - i2
        d2 = ( i1 - i2 ) / i1
        
        d1 = ( '%+.6g' % d1 ) if abs ( d1 ) < tolerance else  attention ( '%+.6g' % d1 )
        d2 = ( '%+.6g' % d2 ) if abs ( d2 ) < tolerance else  attention ( '%+.6g' % d2 )

        d3 =   i1 - ( i3 + i4 ) 
        d4 = ( i1 - ( i3 + i4 ) ) / i1
        
        d3 = ( '%+.6g' % d3 ) if abs ( d3 ) < tolerance else  attention ( '%+.6g' % d3 )
        d4 = ( '%+.6g' % d4 ) if abs ( d4 ) < tolerance else  attention ( '%+.6g' % d4 )

        row = ( '%+.5g' % i0 ,
                '%+.5g' % i1 ,
                '%+.5g' % i2 ,
                '%+.5g' % (i3+i4)  ,
                d1 , d2 , d3 , d4 )

        rows.append ( row )
        
    
    title = '%s' % '1D split' 
    table = T.table ( rows , title = title  , prefix = '# ')
    logger.info ( 'Test %s\n%s' % ( title , table ) )

# =============================================================================
def test_legendresum_2 ( ) :

    logger = getLogger ("test_legendresum_2")
    logger.info ( "Test for 2D splits" )

    xmin , xmax = 0 , 2

    rows = [ ( 'i0' , 'i1' , 'i2' , 'i3+i4' ,
               'i1-i2'   ,
               'i1/i2-1' ,
               'i1-(i3+i4)'  , 
               'i1/(i3+i4)-1' ) ]
    
    for ix in range ( 5 ) :
        lsx = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( ix + 1 ) ] , xmin , xmax )
        for iy in range ( 6 ) :
            lsy = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( iy + 1 ) ] , xmin , xmax )

            ls = LS2 ( lsx , lsy )

            for i in range ( ls.nx() ) :
                for j in range ( ls.ny() ) :
                    ls.setPar ( i , j , random.uniform ( 1 , 5 ) ) 
                    
            lx = ls.integralX ()
            ly = ls.integralY ()

            for l in ( lx , ly ) :

                i1 = l.integral ()
                i2 = l.integral ( xmin + delta , xmax - delta )
                
                split = random.uniform ( xmin , xmax )
                
                i3 = l.integral ( xmin   , split )
                i4 = l.integral ( split  , xmax  )
                
                i0 = I.integral  ( l , xmin = xmin , xmax = xmax )

                d1 =   i1 - i2
                d2 = ( i1 - i2 ) / i1
                
                d1 = ( '%+.6g' % d1 ) if abs ( d1 ) < tolerance else  attention ( '%+.6g' % d1 )
                d2 = ( '%+.6g' % d2 ) if abs ( d2 ) < tolerance else  attention ( '%+.6g' % d2 )
                
                d3 =   i1 - ( i3 + i4 ) 
                d4 = ( i1 - ( i3 + i4 ) ) / i1
                
                d3 = ( '%+.6g' % d3 ) if abs ( d3 ) < tolerance else  attention ( '%+.6g' % d3 )
                d4 = ( '%+.6g' % d4 ) if abs ( d4 ) < tolerance else  attention ( '%+.6g' % d4 )
                
                
                row = ( '%+.5g' % i0 , '%+.5g' % i1 , 
                        '%+.5g' % i2 , '%+.5g' % (i3+i4)  ,
                        d1 , d2 , d3 , d4 )
                
                rows.append ( row )
                                
    
    title = '%s' % '2D split' 
    table = T.table ( rows , title = title  , prefix = '# ')
    logger.info ( 'Test %s\n%s' % ( title , table ) )
        

# =============================================================================
def test_legendresum_3 ( ) :

    logger = getLogger ("test_legendresum_3")
    logger.info ( "Test for 3D splits" )

    xmin , xmax = 0 , 2

    rows = [ ( 'i0' , 'i1' , 'i2' , 'i3+i4' ,
               'i1-i2'   ,
               'i1/i2-1' ,
               'i1-(i3+i4)'  , 
               'i1/(i3+i4)-1' ) ]
    
    for ix in range ( 4 ) :
        lsx = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( ix + 1 ) ] , xmin , xmax )
        for iy in range ( 4 ) :
            lsy = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( iy + 1 ) ] , xmin , xmax )
            for iz in range ( 4 ) :
                lsz = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( iz + 1 ) ] , xmin , xmax )
                
                ls = LS3 ( lsx , lsy , lsz )
                
                for i in range ( ls.nx() ) :
                    for j in range ( ls.ny() ) :
                        for k in range ( ls.nz() ) :
                            ls.setPar ( i , j , k , random.uniform ( 1 , 5 ) ) 
                            
                l2x = ls.integralX ()
                l2y = ls.integralY ()
                l2z = ls.integralZ ()
                
                for l2 in ( l2x , l2y , l2z ) :

                    lx = l2.integralX ()
                    ly = l2.integralY ()

                    for l in ( lx , lx ) :
                        
                        i1 = l.integral ()
                        i2 = l.integral ( xmin + delta , xmax - delta )
                        
                        split = random.uniform ( xmin , xmax )
                        
                        i3 = l.integral ( xmin   , split )
                        i4 = l.integral ( split  , xmax  )
                        
                        i0 = I.integral  ( l , xmin = xmin , xmax = xmax )

                        d1 =   i1 - i2
                        d2 = ( i1 - i2 ) / i1
                        
                        d1 = ( '%+.6g' % d1 ) if abs ( d1 ) < tolerance else  attention ( '%+.6g' % d1 )
                        d2 = ( '%+.6g' % d2 ) if abs ( d2 ) < tolerance else  attention ( '%+.6g' % d2 )
                        
                        d3 =   i1 - ( i3 + i4 ) 
                        d4 = ( i1 - ( i3 + i4 ) ) / i1
                        
                        d3 = ( '%+.6g' % d3 ) if abs ( d3 ) < tolerance else  attention ( '%+.6g' % d3 )
                        d4 = ( '%+.6g' % d4 ) if abs ( d4 ) < tolerance else  attention ( '%+.6g' % d4 )
                        
                        
                        row = ( '%+.5g' % i0 ,
                                '%+.5g' % i1 ,
                                '%+.5g' % i2 , '%+.5g' % (i3+i4)  ,
                                d1 , d2 , d3 , d4 )
                        
                        rows.append ( row )
        
                        
                
    title = '%s' % '3D split' 
    table = T.table ( rows , title = title  , prefix = '# ')
    logger.info ( 'Test %s\n%s' % ( title , table ) )
        


# =============================================================================
def test_legendresum_4 ( ) :

    logger = getLogger ("test_legendresum_4")
    logger.info ( "Test reductions for 2D functions" ) 
    
    xmin , xmax = 0 , 2
    
    rows = [ ( 'i1' , 'i2' , 
               'i1-i2'   ,
               'i1/i2-1' ) ]


    for ix in range ( 4 ) :
        lsx = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( ix + 1 ) ] , xmin , xmax )
        for iy in range ( 4 ) :
            lsy = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( iy + 1 ) ] , xmin , xmax )
            
            ls = LS2 ( lsx , lsy  )
                
            for i in range ( ls.nx() ) :
                for j in range ( ls.ny() ) :
                    ls.setPar ( i , j , random.uniform ( 1 , 5 ) ) 

            y1 = random.uniform   ( xmin , xmax )
            y2 = random.uniform   ( y1   , xmax )
            
            i1 = ls.integralY(y1,y2).integral()
            i2 = ls.integralX().integral(y1,y2)

            d1 =   i1 - i2
            d2 = ( i1 - i2 ) / i1
            
            d1 = ( '%+.9g' % d1 ) if abs ( d1 ) < tolerance else  attention ( '%+.9g' % d1 )
            d2 = ( '%+.9g' % d2 ) if abs ( d2 ) < tolerance else  attention ( '%+.9g' % d2 )
            
            
            row = ( '%+.5g' % i1 , '%+.5g' % i2 , d1 , d2 )
            rows.append ( row )
            
    title = '%s' % '2D reduction' 
    table = T.table ( rows , title = title  , prefix = '# ')
    logger.info ( 'Test %s\n%s' % ( title , table ) )
        

# =============================================================================
def test_legendresum_5 ( ) :

    logger = getLogger ("test_legendresum_5")
    logger.info ( "Test reductions for 3D functions" ) 

    xmin , xmax = 0 , 2
    
    rows = [ ( 'i1' , 'i2' , 
               'i1-i2'   ,
               'i1/i2-1' ) ]

    for ix in range ( 4 ) :
        lsx = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( ix + 1 ) ] , xmin , xmax )
        for iy in range ( 4 ) :
            lsy = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( iy + 1 ) ] , xmin , xmax )
            for iz in range ( 4 ) :
                lsz = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( iz + 1 ) ] , xmin , xmax )
                
                ls = LS3 ( lsx , lsy , lsz )
                
                for i in range ( ls.nx() ) :
                    for j in range ( ls.ny() ) :
                        for k in range ( ls.nz() ) :
                            ls.setPar ( i , j , k , random.uniform ( 1 , 5 ) ) 

                z1 = random.uniform   ( xmin , xmax )
                z2 = random.uniform   ( z1   , xmax )

                i1 = ls.integralZ(z1,z2).integralY().integral()
                i2 = ls.integralY().integralX().integral(z1,z2)

                d1 =   i1 - i2
                d2 = ( i1 - i2 ) / i1

                d1 = ( '%+.9g' % d1 ) if abs ( d1 ) < tolerance else  attention ( '%+.9g' % d1 )
                d2 = ( '%+.9g' % d2 ) if abs ( d2 ) < tolerance else  attention ( '%+.9g' % d2 )
                
                row = ( '%+.5g' % i1 , '%+.5g' % i2 , d1 , d2 )
                rows.append ( row )
                
    title = '%s' % '3D reduction' 
    table = T.table ( rows , title = title  , prefix = '# ')
    logger.info ( 'Test %s\n%s' % ( title , table ) )

# =============================================================================
def test_legendresum_6 ( ) :

    logger = getLogger ("test_legendresum_6")
    logger.info ( "Test reductions for 4D functions" ) 

    xmin , xmax = 0 , 2
    
    rows = [ ( 'i1' , 'i2' , 
               'i1-i2'   ,
               'i1/i2-1' ) ]

    for ix in range ( 4 ) :
        lsx = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( ix + 1 ) ] , xmin , xmax )
        for iy in range ( 4 ) :
            lsy = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( iy + 1 ) ] , xmin , xmax )
            for iz in range ( 4 ) :
                lsz = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( iz + 1 ) ] , xmin , xmax )
                for iu in range ( 4 ) :
                    lsu = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( iu + 1 ) ] , xmin , xmax )
                
                    ls = LS4 ( lsx , lsy , lsz , lsu )
                    
                    for i in range ( ls.nx() ) :
                        for j in range ( ls.ny() ) :
                            for k in range ( ls.nz() ) :
                                for m in range ( ls.nu() ) :
                                    ls.setPar ( i , j , k , m , random.uniform ( 1 , 5 ) ) 

                    u1 = random.uniform   ( xmin , xmax )
                    u2 = random.uniform   ( u1   , xmax )

                    i1 = ls.integralU(u1,u2).integralZ().integralY().integral()
                    i2 = ls.integralZ().integralY().integralX().integral(u1,u2)
                    
                    d1 =   i1 - i2
                    d2 = ( i1 - i2 ) / i1
                    
                    d1 = ( '%+.9g' % d1 ) if abs ( d1 ) < tolerance else  attention ( '%+.9g' % d1 )
                    d2 = ( '%+.9g' % d2 ) if abs ( d2 ) < tolerance else  attention ( '%+.9g' % d2 )
                    
                    row = ( '%+.5g' % i1 , '%+.5g' % i2 , d1 , d2 )
                    rows.append ( row )
                
    title = '%s' % '4D reduction' 
    table = T.table ( rows , title = title  , prefix = '# ')
    logger.info ( 'Test %s\n%s' % ( title , table ) )

    
# =============================================================================
def test_legendresum_7 ( ) :

    logger = getLogger ("test_legendresum_7")

    logger.info ( "Test 2D integrals" ) 

    xmin , xmax = 0 , 2
    
    rows = [ ( 'i0' , 'i1' , 
               'i1-i0'   ,
               'i1/i0-1' ) ]

    for ix in range ( 4 ) :
        lsx = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( ix + 1 ) ] , xmin , xmax )
        for iy in range ( 4 ) :
            lsy = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( iy + 1 ) ] , xmin , xmax )
            
            ls = LS2 ( lsx , lsy  )
                
            for i in range ( ls.nx() ) :
                for j in range ( ls.ny() ) :
                    ls.setPar ( i , j , random.uniform ( 1 , 5 ) ) 


            x1 = random.uniform   ( xmin , xmax )
            x2 = random.uniform   ( 11   , xmax )
            y1 = random.uniform   ( xmin , xmax )
            y2 = random.uniform   ( y1   , xmax )
            
            i1 = ls.integral (      x1 , x2 , y1 , y2 )
            i0 = I.integral2 ( ls , x1 , x2 , y1 , y2 , epsabs = 1.e-9 , epsrel = 1.e-9 ) 
            
            d1 =   i1 - i0
            d2 = ( i1 - i0 ) / i1
            
            d1 = ( '%+.9g' % d1 ) if abs ( d1 ) < tolerance else  attention ( '%+.9g' % d1 )
            d2 = ( '%+.9g' % d2 ) if abs ( d2 ) < tolerance else  attention ( '%+.9g' % d2 )
            
            
            row = ( '%+.5g' % i0 , '%+.5g' % i1 , d1 , d2 )
            rows.append ( row )
            
    title = '%s' % '2D integrals' 
    table = T.table ( rows , title = title  , prefix = '# ')
    logger.info ( 'Test %s\n%s' % ( title , table ) )
        

# =============================================================================
def test_legendresum_8 ( ) :

    logger = getLogger ("test_legendresum_8")
    
    logger.info ( "Test 3D integrals" ) 
                  
    xmin , xmax = 0 , 2
    
    rows = [ ( 'i0' , 'i1' , 
               'i1-i0'   ,
               'i1/i0-1' ) ]

    for ix in range ( 2 , 4 ) :
        lsx = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( ix + 1 ) ] , xmin , xmax )
        for iy in range ( 2 , 4 ) :
            lsy = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( iy + 1 ) ] , xmin , xmax )
            for iz in range ( 2 , 4 ) :
                lsz = LS1 ( [ random.uniform ( 1 , 5 ) for j in range ( iz + 1 ) ] , xmin , xmax )
                
                ls = LS3 ( lsx , lsy , lsz )
                
                for i in range ( ls.nx() ) :
                    for j in range ( ls.ny() ) :
                        for k in range ( ls.nz() ) :
                            ls.setPar ( i , j , k , random.uniform ( 1 , 5 ) ) 


                x1 = random.uniform   ( xmin , xmax )
                x2 = random.uniform   ( 11   , xmax )
                y1 = random.uniform   ( xmin , xmax )
                y2 = random.uniform   ( y1   , xmax )
                z1 = random.uniform   ( xmin , xmax )
                z2 = random.uniform   ( z1   , xmax )
                
                i1 = ls.integral (      x1 , x2 , y1 , y2 , z1 , z2 )
                i0 = I.integral3 ( ls , x1 , x2 , y1 , y2 , z1 , z2 , epsabs = 1.e-9 , epsrel = 1.e-9 ) 
                
                d1 =   i1 - i0
                d2 = ( i1 - i0 ) / i1
                
                d1 = ( '%+.9g' % d1 ) if abs ( d1 ) < tolerance else  attention ( '%+.9g' % d1 )
                d2 = ( '%+.9g' % d2 ) if abs ( d2 ) < tolerance else  attention ( '%+.9g' % d2 )
                
                
                row = ( '%+.5g' % i0 , '%+.5g' % i1 , d1 , d2 )
                rows.append ( row )

                
    title = '%s' % '3D integrals' 
    table = T.table ( rows , title = title  , prefix = '# ')
    logger.info ( 'Test %s\n%s' % ( title , table ) )

# =============================================================================
if '__main__' == __name__ :
    
    test_legendresum_1 ()
    test_legendresum_2 ()
    test_legendresum_3 ()
    test_legendresum_4 ()
    test_legendresum_5 ()
    test_legendresum_6 ()    
    test_legendresum_7 ()
    test_legendresum_8 ()

# =============================================================================
##                                                                      The END  
# ============================================================================= 
