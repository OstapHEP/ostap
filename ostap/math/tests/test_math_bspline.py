#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_bspline.py
#  Test module for the file ostap/math/bspline.py
# ============================================================================= 
""" Test module for ostap/math/bspline.py
"""
# ============================================================================= 
import random  
import ostap.math.models 
import ostap.math.bspline
from   ostap.core.core       import Ostap, SE 
from   ostap.utils.timing    import timing 
from   ostap.core.meta_info  import root_version_int 
from   ostap.plotting.canvas import use_canvas
from   ostap.utils.utils     import wait
from   ostap.math.models     import f1_draw
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_math_bspline' ) 
else                       : logger = getLogger ( __name__            )
# ============================================================================= 
functions = set()

# ============================================================================
##  test solution of equation  B(x) = c
def test_solve ():
    """Test solution of equation  B(x) = c 
    """

    # 1) construct spline with known roots
    
    ## roots in [0,1]
    troots  =          [  random.uniform(0.01,0.99) for i in  range(4) ]
    ## roots in [1,10]
    roots   = troots + [  random.uniform(1.01,9.99) for i in  range(4) ]
    ## complex roots  
    croots  = [  complex ( random.uniform(-1,3),random.gauss(0,1) )  for i in  range(4) ]

    ## Bernstein polynomial with known roots
    bp = Ostap.Math.Bernstein (  0 , 1 , roots , croots )

    ## convert Bernstein polynomial into B-spline 
    bs = Ostap.Math.BSpline   ( bp )

    ## add several internal knots ...
    for i in range(6) :
        x = random.uniform ( 0.001 , 0.999 )
        bs.insert ( x )
        
    ##  find root of B-spline:
        
    rr = bs.solve()
    logger.info ('Roots found : %s' % list( rr ) )
    
    troots.sort() 
    logger.info ('Roots true  : %s' % troots     )

    functions.add ( bs )

    with wait ( 3 ) , use_canvas ( 'test_solve' ) :
        bp.draw (          linecolor = 2 )
        bs.draw ( 'same' , linecolor = 4 )
        

# ============================================================================
##  test spline interpolation 
def test_interpolation ():
    """Test spline interpolation
    """

    if 62006 <= root_version_int :
        logger.warning ("Test_interpolation segfaults for ROOT %s" %  root_version_int ) 
        return 
    
    from math import sin,pi, sqrt

    fun = lambda x  : sin(2*pi*x)
    
    bs  =  ostap.math.bspline.interpolate ( fun , None, [0] + [  random.uniform(0.01,0.99) for i in range(30) ] + [1] , 2 )

    from ostap.stats.counters import SE
    s = SE()
    for i in range(10000) :
        x = random.uniform ( 0 , 1 )
        vf = fun(x)
        vb = bs (x) 
        s += vf-vb
            
    logger.info ('Interpolation quality %s' % s )

    functions.add ( bs )

    with wait ( 3 ) , use_canvas ( 'test_interpolation' ) :
        f1_draw ( fun , xmin = 0 , xmax = 1 , linecolor = 2 ) 
        bs.draw ( 'same' , linecolor = 4 )
        
# ============================================================================
##  test spline approxmation
def test_approximation  ():
    """Test spline approximation
    """
    from math import sin,pi, sqrt

    fun = lambda x  : sin(2*pi*x) 
    bs  =  ostap.math.bspline.approximate ( fun , [0] + [  random.uniform(0.01,0.99) for i in range(50) ] + [1] , 2 )

    from ostap.stats.counters import SE
    s = SE()
    for i in range(10000) :
        x = random.uniform ( 0 , 1 )
        vf = fun(x)
        vb = bs (x) 
        s += vf-vb
            
    logger.info ('Approximation quality %s' % s )

    functions.add ( bs )

    with wait ( 3 ) , use_canvas ( 'test_approximation' ) :
        f1_draw ( fun , xmin = 0 , xmax = 1 , linecolor = 2 ) 
        bs.draw ( 'same' , linecolor = 4 )

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

    test_solve         ()
    test_interpolation ()
    test_approximation ()
    
    ## check finally that everything is serializeable:
    test_pickle ()    
    with timing ('test_db' , logger ) :
        test_db ()
        
# =============================================================================
##                                                                      The END 
# =============================================================================
