#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
""" Test module for ostap/math/bernstein.py.
"""
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_bernstein' ) 
else                       : logger = getLogger ( __name__         )
# ============================================================================= 
import random  
import ostap.math.models 
import ostap.math.bernstein
from   ostap.core.core  import Ostap

# ============================================================================
##  test solution of equation  B(x) = c
def test_solve ():
    """Test solution of equation  B(x) = c 
    """

    # 1) construct function with known roots
    
    ## roots in [0,1]
    troots  =          [  random.uniform(0.01,0.99) for i in  range(5) ]
    ## roots in [1,10]
    roots   = troots + [  random.uniform(1.01,9.99) for i in  range(4) ]
    ## complex roots  
    croots  = [  complex ( random.uniform(-3,0),random.gauss(0,3) )  for i in  range(4) ]

    ## Bernstein polynomial with known roots
    bs = Ostap.Math.Bernstein (  0 , 1 , roots , croots )

    ##  find root of Bernstein  polynomial 
        
    rr = bs.solve()
    logger.info ('Roots found : %s' % list( rr ) )
    
    troots.sort() 
    logger.info ('Roots true  : %s' % troots     )

# ============================================================================
##  check number of roots using Sturm' sequence 
def test_nroots ():
    """Sheck number of roots using Sturm' sequence 
    """

    # 1) construct function with known roots
    
    ## roots in [0,1]
    troots  =          [  random.uniform(0.01,0.99) for i in  range(5) ]
    ## roots in [1,10]
    roots   = troots + [  random.uniform(1.01,9.99) for i in  range(4) ]
    ## complex roots  
    croots  = [  complex ( random.uniform(-3,0),random.gauss(0,3) )  for i in  range(4) ]

    ## Bernstein polynomial with known roots
    bs = Ostap.Math.Bernstein (  0 , 1 , roots , croots )


    troots.sort() 
    logger.info ('Roots true  : %s' % troots     )
    
    for i in  range(20)  :
        x1 = random.uniform (  0 , 0.70 )
        x2 = random.uniform ( x1 , 1.00 )
        x2 = bs.xmax() 
        nr = bs.nroots ( x1 , x2 )
        nc = 0
        for r in  troots :
            if x1 < r <= x2 :  nc=+1 
        logger.info ('Roots between [%s,%s) : %d/%d '  % ( x1  , x2 , nr , nc ) )
        
       
# ============================================================================
##  test Bernstein interpolation 
def test_interpolation ():
    """Test bernstein interpolation
    """
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


# ============================================================================
##  test polynomial divisions
def test_division ():
    """Test polynomial division
    """

    # 1) construct function with known roots
    
    ## roots in [0,1]
    troots  =          [  random.uniform(0.01,0.99) for i in  range(5) ]
    ## complex roots  
    croots  = [ complex ( random.uniform(-3, 0),random.gauss(0,3) )  for i in range(3) ]

    troots.sort() 
    logger.info ('Roots true  : %s' % troots     )

    ## Bernstein polynomial with known roots
    bs = Ostap.Math.Bernstein (  0 , 1 , troots      , croots )
    b1 = Ostap.Math.Bernstein (  0 , 1 , troots[:1]  )
    b2 = Ostap.Math.Bernstein (  0 , 1 , troots[1:3] )
    b3 = Ostap.Math.Bernstein (  0 , 1 , troots      )

    a1,r1 = divmod ( bs, b1 )
    logger.info ('Reminder, roots: %s %s'  %  (  r1.norm () , a1.solve() ) )
    a2,r2 = divmod ( bs, b2 )
    logger.info ('Reminder, roots: %s %s'  %  (  r2.norm () , a2.solve() ) )
    a3,r3 = divmod ( bs, b3 )
    logger.info ('Reminder, roots: %s %s'  %  (  r3.norm () , a3.solve() ) )
    
                 

# =============================================================================
if '__main__' == __name__ :
        
    test_solve         ()
    test_nroots        ()
    test_interpolation ()
    test_division      ()
    
# =============================================================================
# The END 
# =============================================================================
