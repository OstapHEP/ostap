#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# =============================================================================
## @file ostap/math/tests/test_math_bernstein.py
#  Test module for the file ostap/math/bernstein.py
# ============================================================================= 
""" Test module for the file ostap/math/bernstein.py
"""
# ============================================================================= 
import random  
import ostap.math.models 
import ostap.math.bernstein
from   ostap.core.core      import Ostap, SE 
from   ostap.utils.timing   import timing 
# ============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_bernstein' ) 
else                       : logger = getLogger ( __name__                    )
# ============================================================================= 
functions = set() 
# ============================================================================
##  test solution of equation  B(x) = c
def test_solve ():
    """Test solution of equation  B(x) = c 
    """
    logger = getLogger("test_solve")
    
    # 1) construct function with known roots
    
    ## roots in [0,1]
    troots  =          [  random.uniform(0.05,0.95) for i in  range ( 4 ) ]
    troots.sort() 
    
    ## roots in [2,9]
    roots   = troots + [  random.uniform(2,9) for i in  range ( 4 ) ]
    roots.sort ()
    
    ## complex roots  
    croots  = [ complex ( random.uniform(-6,-1) , random.uniform ( -5 , 5 ) ) for i in range ( 3 ) ]
    croots += [ complex ( random.uniform( 3, 6) , random.uniform ( -5 , 5 ) ) for i in range ( 3 ) ]

    ## Bernstein polynomial with known roots
    bs = Ostap.Math.Bernstein (  0 , 1 , roots ) ##  , croots )

    logger.info ( "Bernstein: %s" % bs ) 
    
    ##  find roots of Bernstein  polynomial     
    rr = bs.solve()
    rr = [ r for r in rr ]
    rr.sort()
    
    logger.info ('Roots found : [ %s]' %  ( ', '.join ( "%.6f" % r for r in rr     ) ) )    
    logger.info ('Roots true  : [ %s]' %  ( ', '.join ( "%.6f" % r for r in troots ) ) )

    if len ( rr ) != len ( troots ) :
        logger.error ( 'Mismatch in number of roots found!' )
    else :
        diff = 0.0 
        for i,r in enumerate ( rr ) :
            diff += abs ( r  - troots[i] )            
        diff /= len ( rr ) 
        logger.info ( 'Mean root distance is %.4g' % diff ) 

    functions.add ( bs ) 

# ============================================================================
##  check number of roots using Sturm' sequence 
def test_nroots ():
    """Check number of roots using Sturm' sequence 
    """

    logger = getLogger("test_nroots")
    
    # 1) construct function with known roots
    
    ## roots in [0,1]
    troots  =          [  random.uniform(0.01,0.99) for i in  range ( 6 ) ]
    troots.sort() 

    ## roots in [1,10]
    roots   = troots + [  random.uniform(1.01,9.99) for i in  range ( 4 ) ]
    roots.sort ()
    
    ## complex roots  
    croots  = [ complex ( random.uniform(-3 , -1 ) , random.gauss ( 0 , 3 ) ) for i in range ( 4 ) ]
    croots += [ complex ( random.uniform( 3 ,  5 ) , random.gauss ( 0 , 3 ) ) for i in range ( 4 ) ]
    
    ## Bernstein polynomial with known roots
    bs = Ostap.Math.Bernstein (  0 , 1 , roots , croots )


    troots.sort() 
    logger.info ('Roots true  : [%s]' %  ( ', '.join ( "%.6f" % r for r in troots )  ) )
    
    for i in  range( 20 )  :
        
        x1 = random.uniform (  0 , 0.80 )
        x2 = random.uniform ( x1 , 1.00 )

        nr = bs.nroots ( x1 , x2 )
        nt = 0
        for r in  troots :
            if x1 <= r < x2 :  nt += 1

        if nr != nt : 
            logger.error ('Roots between [%.6f, %.6f) : %d [true is %d]'  % ( x1  , x2 , nr , nt ) )
        else : 
            logger.info  ('Roots between [%.6f, %.6f) : %d [true is %d]'  % ( x1  , x2 , nr , nt ) )
        
    functions.add ( bs ) 
       
# ============================================================================
##  test Bernstein interpolation 
def test_interpolation ():
    """Test bernstein interpolation
    """

    logger = getLogger("test_interpolation")

    from math import sin,pi, sqrt
    
    fun = lambda x  : sin(2*pi*x)

    bs  =  ostap.math.bernstein.interpolate ( fun , [0] + [  random.uniform(0.01,0.99) for i in range(25) ] + [1] , 0 , 1 )
    
    from ostap.stats.counters import SE
    s = SE()
    for i in range(10000) :
        x = random.uniform ( 0 , 1 )
        vf = fun(x)
        vb = bs (x) 
        s += vf-vb
            
    logger.info ('Interpolation quality %s' % s )

    functions.add ( bs ) 


# ============================================================================
##  test polynomial divisions
def test_division ():
    """Test polynomial division
    """

    logger = getLogger("test_division")

    # 1) construct function with known roots
    
    ## roots in [0,1]
    troots  =          [  random.uniform(0.01,0.99) for i in  range(5) ]
    troots.sort() 
    ## complex roots  
    croots  = [ complex ( random.uniform(-3, 0),random.gauss(0,3) )  for i in range(3) ]

    logger.info ('Roots true  : [%s]' % ( ', '.join( "%.6f"  % r for r in troots ) ) )    

    ## Bernstein polynomial with known roots
    bs = Ostap.Math.Bernstein (  0 , 1 , troots      , croots )
    b1 = Ostap.Math.Bernstein (  0 , 1 , troots[:1]  )
    b2 = Ostap.Math.Bernstein (  0 , 1 , troots[1:3] )
    b3 = Ostap.Math.Bernstein (  0 , 1 , troots      )

    a1,r1 = divmod ( bs, b1 )
    logger.info ('Reminder, roots: %.6f [%s]' % ( r1.norm () , ', '.join ( "%.6f" % r for r in a1.solve() ) ) )
    a2,r2 = divmod ( bs, b2 )
    logger.info ('Reminder, roots: %.6f [%s]' % ( r2.norm () , ', '.join ( "%.6f" % r for r in a2.solve() ) ) )
    a3,r3 = divmod ( bs, b3 )
    logger.info ('Reminder, roots: %.6f [%s]' % ( r3.norm () , ', '.join ( "%.6f" % r for r in a3.solve() ) ) )

    functions.add ( bs ) 
    functions.add ( b1 ) 
    functions.add ( b2 ) 
    functions.add ( b3 ) 

# =============================================================================
def check_equality ( a ,  b , message = '' , tolerance = 1.e-7 ) :
    d = abs ( a - b ) 
    if d > tolerance :
        raise ValueError ( message + ' |%s-%s|=%s>%s' % ( a , b , d , tolerance ) )

# =========================================================================----
## test for elevate/reduce
def test_elevatereduce () :

    logger = getLogger("test_elevatereduce")

    BP = Ostap.Math.Bernstein
    
    b  = BP ( 5 , 0. , 2. ) ## 5th order for x in [0,2]
    for i in b :
        b[i] = random.uniform ( -10 , 10 )

    functions.add ( b ) 

    for r in (1,2,3,4,5) : 
        be = b.elevate(r)
        br = be.reduce(r)
        for i in range(50) :
            x = random.uniform ( b.xmin() , b.xmax() )
            y  = b(x)
            ye = be(x)
            yr = br(x)
            check_equality ( y  , ye  , 'Invalid elevate' , 1.e-6 )
            check_equality ( y  , yr  , 'Invalid reduce'  , 1.e-6 )

        functions.add ( be ) 
        functions.add ( br ) 

    logger.info ('Elevate/reduce  is OK' )

# ==============================================================================
## test  for polynomials 
def test_poly () :
    
    logger = getLogger("test_poly")

    BP = Ostap.Math.Bernstein
    
    # 1) create & evaluate the polynom 
    
    b  = BP ( 5 , 0. , 2. ) ## 5th order for x in [0,2]
    b += 1 
    
    for i in  range(500) :
        x = random.uniform ( b.xmin() , b.xmax() )
        y = b ( x ) 
        check_equality ( y , 1 , 'Invalid polynom value:' )
        
    logger.info ('Constant   poly is OK' )

    # evaluate the polynom
    for i in  b  : b[i] = random.uniform(-10,10)
    ymin = min ( b.pars() )
    ymax = max ( b.pars() )
    
    for i in  range(100) :
        x = random.uniform ( b.xmin() , b.xmax() )
        y = b(x)
        if not ymin <= y <= ymax : 
            raise ValueError ( 'Invalid polynom value y(%s)=%s (%s/%s)' % ( x , y , ymin , ymax ) )
 
    logger.info ('Random     poly is OK' )

    functions.add ( b ) 

# ==============================================================================
## test  for even polynomials 
def test_even () :

    
    logger = getLogger("test_even")

    BP = Ostap.Math.BernsteinEven
        
    b  = BP ( 5 , 0. , 2. ) ## 5th order for x in [0,2]
    for i in  b  : b[i] = random.uniform(-10,10)
    
    xmid = 0.5 * ( b.xmin() + b.xmax() )

    for i in range(100) :
        x1 = random.uniform ( b.xmin() , b.xmax() )
        dx = x1   - xmid  
        x2 = xmid - dx 
        y1 = b(x1)
        y2 = b(x2)
        check_equality ( y1 , y2  , 'Invalid BernsteinEven'  , 1.e-7 )

    
    logger.info ('Even       poly is OK' )

    functions.add ( b ) 

# ==============================================================================
## test  for monotonic polynomial 
def test_monotonic () :
    
    logger = getLogger("test_monotonic")

    ## 8-9) check for Monotonic 
    BPM = Ostap.Math.Monotonic
    
    b = BPM ( 5 , 0 , 2 , True )
    for i in  b :  b[i] = random.uniform ( -10 , 10 )
    
    for i in range(500) :
        x1 = random.uniform ( b.xmin() , b.xmax() )
        y1 = b(x1)
        if y1 <= 0 :
            raise ValueError ( 'Invalid Increasing y(%s)=%s' % ( x1 , y1 ) )
        x2 = random.uniform ( x1 , b.xmax() )
        y2 = b(x2)
        if y2 <= 0 :
            raise ValueError ( 'Invalid Increasing y(%s)=%s' % ( x2 , y2 ) )
        if y1 > y2 :
            raise ValueError ( 'Invalid Increasing y(%s)=%s>y(%s)=%s' % ( x1 , y1 , x2 , y1 ) )
        
    logger.info ('Increasing poly is OK' )
    
    functions.add ( b ) 

    b = BPM ( 5 , 0 , 2 , False )
    for i in  b :  b[i] = random.uniform ( -10 , 10 )
    
    for i in range(500) :
        x1 = random.uniform ( b.xmin() , b.xmax() )
        y1 = b(x1)
        if y1 <= 0 :
            raise ValueError ( 'Invalid Decreasing y(%s)=%s' % ( x1 , y1 ) )
        x2 = random.uniform ( x1 , b.xmax() )
        y2 = b(x2)
        if y2 <= 0 :
            raise ValueError ( 'Invalid Decreasing y(%s)=%s' % ( x2 , y2 ) )
        if y1 < y2 :
            raise ValueError ( 'Invalid Decreasing y(%s)=%s>y(%s)=%s' % ( x1 , y1 , x2 , y1 ) )
        
    functions.add ( b ) 

    logger.info ('Decreasing poly is OK' )

# =============================================================================
## test for Convex positive polynmomials 
def test_convex () :
    """Test for Convex positive polynmomials 
    """
    
    logger = getLogger("test_convex")
    BPC = Ostap.Math.Convex
    
    b_11 = BPC ( 5 , 0 , 2 , True  , True  )
    b_01 = BPC ( 5 , 0 , 2 , False , True  )
    b_10 = BPC ( 5 , 0 , 2 , True  , False )
    b_00 = BPC ( 5 , 0 , 2 , False , False )
    
    for b in ( b_11   , b_01 , b_10 , b_00 ) : 
        for i in b : b[i] = random.uniform ( -10 , 10 ) 
            
    functions.add ( b_11 ) 
    functions.add ( b_01 ) 
    functions.add ( b_10 ) 
    functions.add ( b_00 ) 

    for i in range(500) :
        
        ##  note:  x1 <  x2 
        x1 = random.uniform ( b.xmin() , b.xmax() )
        x2 = random.uniform ( x1       , b.xmax() )
        xm = 0.5 * ( x1 + x2 )
        
        for b in ( b_11   , b_01 , b_10 , b_00 ) :

            b1 = b ( x1 ) 
            b2 = b ( x2 ) 

            if  b1 < 0 or b2 < 0 :
                raise ValueError ( 'Invalid Convex value (b<0)' )

            if b.increasing() and b1 > b2 :
                    raise ValueError ( 'Invalid Convex Increasing' )
                
            if b.decreasing() and b1 < b2 :
                raise ValueError ( 'Invalid Convex Decreasing' )
            
            bm = b ( xm ) 
            if b.convex  () and ( b1 + b2 ) < 2 * bm :
                    raise ValueError ( 'Invalid Convex!' )                    
            if b.concave () and ( b1 + b2 ) > 2 * bm :
                    raise ValueError ( 'Invalid Concave!' )
                                
    logger.info ('Convex     poly is OK' )


# =============================================================================
## test for Convex positive polynmomials 
def test_convexonly () :
    """Test for Convex-Only positive polynmomials 
    """
    
    logger = getLogger("test_convexonly")
    BPC = Ostap.Math.ConvexOnly
    
    b_11 = BPC ( 5 , 0 , 2 , True  )
    b_01 = BPC ( 5 , 0 , 2 , False )
    
    for b in ( b_11   , b_01 ) : 
        for i in b  : b[i] = random.uniform ( -10 , 10 ) 

    functions.add ( b_11 ) 
    functions.add ( b_01 ) 

    for i in range(500) :
        
        ##  note:  x1 <  x2 
        x1 = random.uniform ( b.xmin() , b.xmax() )
        x2 = random.uniform ( x1       , b.xmax() )
        xm = 0.5 * ( x1 + x2 )
        
        for b in ( b_11   , b_01 ) :

            b1 = b ( x1 ) 
            b2 = b ( x2 ) 

            if  b1 < 0 or b2 < 0 :
                raise ValueError ( 'Invalid Convex value (b<0)' )

            bm = b ( xm ) 
            if b.convex  () and ( b1 + b2 ) < 2 * bm :
                    raise ValueError ( 'Invalid Convex!' )                    
            if b.concave () and ( b1 + b2 ) > 2 * bm :
                    raise ValueError ( 'Invalid Concave!' )
                                
    logger.info ('ConvexOnly poly is OK' )


# =============================================================================
## test integration of 
def test_integration () :
    """Test for polynomial interation
    """
    
    logger = getLogger("test_integration")
    
    BP = Ostap.Math.Bernstein
    
    b = BP ( 5 , 0 , 2  )
    for i in b  : b[i] = random.uniform ( -1 , 50 ) 

    from ostap.math.integral import romberg 

    xmax = b.xmin() + 0.9 * ( b.xmax() - b.xmin() ) 
    for i in range ( 500 ) :
        
        ##  note:  x1 <  x2 
        x1 = random.uniform ( b.xmin() ,   xmax   )
        x2 = random.uniform ( x1       , b.xmax() )

        i1 = b.integral (      x1 , x2 )
        i2 = romberg    ( b ,  x1 , x2 )

        dd = ( i1 - i2 )
        ds = ( abs ( i1 ) + abs ( i2 ) )
        if 0 < ds < abs ( dd ) * 10**6 :
            raise ValueError ( 'Invalid Integrals! %s  %s' % ( dd , ds*1.e-6 ) )

    logger.info ('Integration     is OK' )

# =============================================================================
## test transformations
def test_transformation () :
    """Test for transformation
    """
    
    logger = getLogger("test_transformation")

    
    BP = Ostap.Math.Bernstein
    MS = Ostap.Math.Polynomial   ## monomial sum 
    CS = Ostap.Math.ChebyshevSum ## chebyshev sum 
    LS = Ostap.Math.LegendreSum  ## legendre sum 
    
    
    b = BP ( 5 , 0 , 2  )

    functions.add ( b ) 

    for i in range ( 50 ) :
        
        for i in b  : b[i] = random.uniform ( -10 , 10 ) 
        
        ## Monomial sum 
        ms = MS ( b  )
        bm = BP ( ms ) 
        
        ## Chebyshev sum 
        cs = CS ( b  )
        bc = BP ( cs ) 
        
        ## Legendre  sum 
        ls = LS ( b  )
        bl = BP ( ls ) 

        functions.add ( ms ) 
        functions.add ( bm )
        
        functions.add ( cs ) 
        functions.add ( bc )
        
        functions.add ( ls ) 
        functions.add ( bl ) 

        for i in range( 100 ) :
            
            x1 = random.uniform ( b.xmin() , b.xmax() )
            
            bv = b ( x1 )
            
            check_equality ( bv , ms ( x1 ) , 'Invalid Bernstein->Monomial  transformation' )
            check_equality ( bv , bm ( x1 ) , 'Invalid Bernstein<-Monomial  transformation' )
            
            check_equality ( bv , cs ( x1 ) , 'Invalid Bernstein->Chebyshev transformation' )
            check_equality ( bv , bc ( x1 ) , 'Invalid Bernstein<-Chebyshev transformation' )
            
            check_equality ( bv , ls ( x1 ) , 'Invalid Bernstein->Legendre  transformation' )
            check_equality ( bv , bl ( x1 ) , 'Invalid Bernstein<-Legendre  transformation' )
            
    logger.info ('Transformation  is OK' )

# =============================================================================
def test_pickle () :
    logger = getLogger ( 'test_pickle'        ) 
    logger.info ( 'Check pickling/unpickling' )

    import pickle
    rows = [ ( '#', 'before' , 'after' , 'mean' , 'rms' ) ] 
    for i, f in enumerate ( functions , start = 1 ) :

        fs = pickle.loads ( pickle.dumps ( f ) )
        s  = SE () 
        for j in range ( 1000 ) :
            x = random.uniform ( f.xmin() , f.xmax() )
            s += abs ( fs ( x ) - f ( x ) )
        row = '%d' % i , f.__class__.__name__ , fs.__class__.__name__ , '%-+.4g' % s.mean() , '%-+.4g' % s.rms() 
        rows.append ( row )

    import ostap.logger.table as T
    title = "Compare before/after eserialisation"
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'rllll' ) 
    logger.info ( '%s\n%s' % ( title , table ) ) 


# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        for i , f in enumerate ( functions , start = 1 ) :
            db[ '%03d:%s' % ( i , f.__class__.__name__ ) ] = f 
        db['functions'   ] = functions 
        db.ls() 

# =============================================================================
if '__main__' == __name__ :

    test_solve          ()
    test_nroots         ()
    test_interpolation  ()
    test_division       ()
    test_elevatereduce  ()
    test_poly           ()
    test_even           ()
    test_monotonic      ()
    test_convex         () 
    test_convexonly     ()
    test_integration    ()
    test_transformation ()

    ## check finally that everything is serializeable:
    test_pickle () 
    with timing ('test_db' , logger ) :
        test_db ()
        
# =============================================================================
##                                                                      The END 
# =============================================================================
