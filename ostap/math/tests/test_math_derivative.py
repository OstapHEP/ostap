#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_derivative.py
#  Test module for the file ostap/math/derivative.py
# ============================================================================= 
""" Test module for ostap/math/derivative.py

It tests local implementation of numerical derivatives 
"""
# ============================================================================= 
from   math                   import exp, sin, cos, pi, tanh
from   ostap.core.core        import VE 
from   ostap.math.models      import f1_draw
from   ostap.math.derivative  import ( Derivative  , iszero      , Eval2VE     ,
                                       Derivative1 , Derivative2 , Derivative3 ,
                                       Derivative4 , Derivative5 , Derivative6 )                 
from   ostap.math.finitediffs import CentralRule, ForwardOpen, BackwardOpen  
from   ostap.stats.counters   import SE
from   ostap.utils.timing     import timing
import ostap.logger.table     as     T
from   ostap.logger.pretty    import pretty_ve
from   ostap.utils.utils      import wait, batch_env 
from   ostap.plotting.canvas  import use_canvas
import ROOT, random, math 
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_derivative' )
else                       : logger = getLogger ( __name__                     )
# =============================================================================
## set batch form environment 
batch_env ( logger )
# =============================================================================

def derivative_testing ( der_type , func , der_exact , logger , **kwargs ) :

    IMAX = der_type.IMAX
    
    cnt1  = {}
    
    for i in range ( IMAX )  :  cnt1 [i] = SE()
        
    ders = [ der_type ( func , I , **kwargs ) for I in range ( IMAX ) ]
    
    with timing ( 'derivative test' , logger = logger ) :
        
        for i in range(10000) :
            
            x = random.uniform ( 0 , pi )
            
            d_true = der_exact ( x )
            
            for I , dd in enumerate ( ders ) :
                
                delta     = dd ( x ) - d_true
                cnt1 [I] += delta

    rows = [ ('I' , 'rms' , 'min' , 'max' ) ]
    
    for i,c in cnt1.items() :
        row = '%d' % i , '% .4g' % c.rms () , '% .4g' % c.min() , '% .4g' % c.max() 
        rows.append ( row ) 
        
    table = T.table ( rows ,
                      title  = 'Test numerical derivatives' ,
                      prefix = '# ' , alignment = 'clll'     )
    
    logger.info ( 'Test numerical derivatives\n%s' % table )
    

    

# =============================================================================
def test_derivative_1 ():

    logger = getLogger ( 'test_derivative_1' )
    
    fun = math.sin
    der = math.cos

    derivative_testing ( Derivative , fun , der , logger )
    
def test_derivative_2 ():

    logger = getLogger ( 'test_derivative_2' )
    
    fun = math.sin
    der = math.cos

    derivative_testing ( Derivative1 , fun , der , logger ) 

def test_derivative_3 ():

    logger = getLogger ( 'test_derivative_2' )
    
    fun = math.sin
    der = math.cos

    derivative_testing ( Derivative1 , fun , der , logger , richardson = 2 ) 


# =============================================================================
def test_derivative_4 ():

    logger = getLogger ( 'test_derivative_4' )


    functions = (
        ( lambda x : cos(10.*x)   , lambda x : -10*sin(10.*x)                  ) , 
        ( lambda x : x**3         , lambda x : 3.0*x*x                         ) , 
        ( lambda x : exp(x)       , lambda x : exp(x)                          ) ,
        ( lambda x : x**8         , lambda x : 8.0*x**7                        ) ,
        ( lambda x : tanh(2.*x)   , lambda x : 2*(1.-tanh(2.*x)**2)            ) ,
        ( lambda x : 1.11*x       , lambda x : 1.11                            ) , 
        ( lambda x : 1.11111      , lambda x : 0.0                             ) , 
        ( lambda x : x**10        , lambda x : 10.*x**9                        ) , 
        )
    
    from ostap.core.core import SE
    counters = {} 

    from ostap.utils.progress_bar import progress_bar
    
    IMAX  = 8 
    table =  [ ['Function'] + [ 'I=%d' % i for i in range ( IMAX ) ] ] 
    for i , o in enumerate ( progress_bar ( functions ) ) :
        
        fun = o [ 0 ] ## function 
        der = o [ 1 ] ## derivative 
        
        row = [ '%2d' % (i+1) ]
        for I in range ( IMAX ) :
            
            cnt1 = SE ()            
            cnt2 = SE ()
            
            dd   = Derivative ( fun , step = 0.001 , calc = I , with_error = True )

            for j in range ( 1000 ) :
                
                x = random.uniform ( 0.05 , 1.5 )
                res  = dd ( x ) 
                dif  = float ( res ) - der ( x ) 
                cnt1 += dif 
                if res.cov2() > 0 : cnt2 += dif/res.error()  
                
            mmax1 = abs ( cnt1.max ()  *10**12 ) 
            if 2 < cnt2.nEntries() : 
                mmax2 = cnt2.max()
                row.append ( '%7.3f / %-5.2fs' % ( mmax1, mmax2 ) )
            else :
                mmax2 = 0
                row.append ( '%7.3f / %-5.2fs' % ( mmax1, mmax2 ) )
                
        table.append ( row )

    table = T.table ( table , prefix = '# ' , alignment=9*'c' )
    logger.info ('Numerical differentiation: Max difference [10^12]\n%s' % table ) 


# =============================================================================
def test_derivative_5 ():

    logger = getLogger ( 'test_derivative_5' )

    ## the function
    func2 = lambda x,y : 0.25*x*x + (x-y)*(x-y) + y + 2*x
    
    ## use explicit   partial derivatives 
    eval2_1 = Eval2VE( func2 , dFdX = lambda x,y : 2+0.5*x+2*(x-y) ,  dFdY = lambda x,y : 1-2*(x-y) )
    
    ## use numerical  partial derivatives 
    eval2_2 = Eval2VE( func2 )
    
    table = [  ( 'x' , 'y' , 'corr' , 'F(exact)' , 'F(numerical)' ) ] 
    for x , y in [ (0,0) , (1,0) , (0,1) , (1,2) , ( 1 , -1 ) , ( -1 , 1 ) , (2,1) ] :
        
        for c in  ( -1.0 , -0.5 , 0.0 , 0.5 , 1.0 ) :

            x = VE ( x , 0.1 ** 2 )
            y = VE ( y , 0.1 ** 2 )

            row = [ x.toString ( '%+5.2f +/- %-5.2f' ) ,
                    y.toString ( '%+5.2f +/- %-5.2f' ) ,
                    '%+3.1f' %  c  ]

            v1 = eval2_1 ( x , y , c )
            v2 = eval2_2 ( x , y , c )

            row.append ( v1.toString ( '%+6.3f +/- %-5.3f') )
            row.append ( v2.toString ( '%+6.3f +/- %-5.3f') )
            table.append ( tuple ( row ) ) 
            
    title = 'Error propagation for F(x,y)'
    table = T.table ( table , title = title , prefix = '# ' , alignment=5*'c' )
    logger.info ('%s\n%s' % ( title , table ) ) 


# =============================================================================
def test_derivative_6 ():
    """Function with discontiniute derivatives
    """
    
    logger = getLogger ( 'test_derivative_6' )

    ## function
    a        = 1.5 
    fun      = lambda x : abs ( math.sin ( a * x ) )  ## NB! 
    
    ## singular points 
    singular = [ i * math.pi / a  for i in range ( -20 , 21 ) ]

    derivs = (
        Derivative1 ( fun , singular = singular , max_step = 0.5 ) ,
        Derivative2 ( fun , singular = singular , max_step = 0.5 ) ,
        Derivative3 ( fun , singular = singular , max_step = 0.5 ) ,
        Derivative4 ( fun , singular = singular , max_step = 0.5 ) ,
        Derivative5 ( fun , singular = singular , max_step = 0.5 ) ,
        Derivative6 ( fun , singular = singular , max_step = 0.5 ) ,        
        )

    lines = [] 
    with wait ( 5 ) , use_canvas ( 'test_derivative_6' ) :
        xmin, xmax = -4 , 4 
        f1_draw ( fun , xmin = xmin , xmax = xmax  , linecolor=2 , linewidth = 3 , min = -12 , max = 12 )
        for i, d in enumerate ( derivs , start = 1 ) :
            color = 2 + i  
            d.draw ( 'same' , xmin = xmin  , xmax = xmax , linecolor = color , linewidth = 2 )
            v  = a**i 
            l1 = ROOT.TLine ( xmin , -v , xmax , -v )
            l2 = ROOT.TLine ( xmin ,  v , xmax ,  v )
            for l in ( l1 ,l2 ) :
                l.SetLineColor ( color ) 
                l.SetLineStyle ( 8     ) 
                lines.append ( l )
                l.draw() 
            
        

# =============================================================================
def differences_testing ( RULE , logger ) :

    a = 1.5
    fun        =   lambda x :   math.sin ( a * x ) 
    exact_ders = ( fun                     ,           ## 0
                   lambda x :   math.cos ( a * x ) * a             , ## 1
                   lambda x : - math.sin ( a * x ) * a * a         , ## 2 
                   lambda x : - math.cos ( a * x ) * pow ( a , 3 ) , ## 3  
                   lambda x :   math.sin ( a * x ) * pow ( a , 4 ) , ## 4
                   lambda x :   math.cos ( a * x ) * pow ( a , 5 ) , ## 4
                   lambda x : - math.sin ( a * x ) * pow ( a , 6 ) ) ## 4
    
    h0   = 0
    hmax = 1.2
    N    = 1000
    cd   = SE ()
    cp   = SE ()
    ch   = SE ()
    
    for D in range ( 1 , RULE.DMAX + 1 ) :

        table = [  ( 'I' , '#' , 'step' , '' , 'delta' , '' , 'pull' ) ]
        
        for I in range ( RULE.IMAX ( D ) ) :
            
            rule  = RULE ( D , I , with_error = True )
            
            for ik in range ( N ) :
                
                z      = random.uniform ( -math.pi , math.pi )
                
                h , _  = rule.optimal_step ( fun , z , h0 , hmax )
                r      = rule ( fun , z , h )
                exact  = exact_ders [ D ] ( z )
                

                dd     = float ( r ) - exact
                
                cd    += abs ( dd )
                cp    += dd / r.error()
                ch    += h 
                
                
            hh = ch.mean()
            dd = cd.mean()

            hh , nh = pretty_ve ( hh , width = 4 , precision = 3 , parentheses = False )
            dd , nd = pretty_ve ( dd , width = 4 , precision = 3 , parentheses = False )

            row = '%2d' % I , \
                  '%2d' % len( rule.stencil)  , \
                  hh , '' if not nh else '[10^%-3d]' % nh , \
                  dd , '' if not nd else '[10^%-3d]' % nd , \
                  "%-.2f" % cp.rms() 
            table.append ( row )
            

        title = 'Finite differences test for D=%s derivative' %  D  
        table = T.table ( table , title = title , prefix = '# ' , alignment = 'rrlclcr' )
        logger.info ('%s\n%s' % ( title , table ) ) 
        

# =============================================================================
def test_central_differences ( ) :
    logger = getLogger ( 'test_central_differences' )
    return differences_testing ( CentralRule , logger )

def test_forward_differences ( ) :
    logger = getLogger ( 'test_forward_differences' )
    return differences_testing ( ForwardOpen , logger )

def test_backward_differences ( ) :
    logger = getLogger ( 'test_backward_differences' )
    return differences_testing ( BackwardOpen , logger )

# =============================================================================
def topt ( d , n ) :

    fun = lambda t :  n*pow ( t , n + d ) + ( n + d ) * pow ( t , n ) - d
    from ostap.math.rootfinder import solve
    return solve ( fun , 1.e-15 , 1-1.e-15 )

# =============================================================================
## helper function to get table for forward (open) differences
def make_forward_diffs ( order ) :
    """Helper function to get tabel for forward (open) differences
    """
    
    import findiff
    from   ostap.math.base import lcm 
    from   math            import fsum, sqrt 
    from   ostap.math.finitediffs import calc_dot 
    
    def fact ( n ) :
        return 1 if 0 == n else n * fact ( n - 1 ) 


    D = order 
    for I in range ( 1 , 26 ) :
        r = findiff.coefficients(deriv=D, offsets=tuple ( k for k in range (1,I+D+1) ), symbolic=True)
        c = r['coefficients']

        denom  = lcm (  *[ i.q for i in c ] )
        coeffs = [ denom *i for i in c ]
        e = ( fsum (   ( 1.0 * k * k for k in coeffs ) ) ** 0.5 ) / ( denom * sqrt ( 12.0 ) ) 
        if all (  [ float(i) == i for i in coeffs  ] ) :

            
            n       = r['accuracy']
            d       = D 
            stencil = tuple ( range ( 1 , I + D + 1 ) )

            fun   = lambda x : x**(n+d)
            ff    = fact ( n + d )
            a     = denom * ff
            b     = calc_dot ( fun , 0.0 , 1.0 , stencil , coeffs )
            term  = a / b
            
            to    = topt ( d , n ) 

            line = ' RuleConf ( %d , darray ( %s ) , %s , %d , darray  ( range ( 1 , %s ) ) , %.10g , %.8g  , %.4f ) ,' % ( D , coeffs , denom , n , I + D + 1 , e , term , to )
            roff = r['offsets']
            
            assert roff == stencil , 'forward invalid stencil!'
            print ( line )            

# =============================================================================
## helper function to get table for backward (open) differences
def make_backward_diffs ( order ) :
    """Helper function to get table for backward (open) differences
    """
    
    import findiff
    from   ostap.math.base import lcm 
    from   math            import fsum, sqrt 
    from   ostap.math.finitediffs import calc_dot 
    
    def fact ( n ) :
        return 1 if 0 == n else n * fact ( n - 1 ) 


    D = order 
    for I in range ( 1 , 26 ) :
        r = findiff.coefficients(deriv=D, offsets=tuple ( k for k in range (-I-D-1, 0 ) ), symbolic=True)
        c = r['coefficients']

        denom  = lcm (  *[ i.q for i in c ] )
        coeffs = [ denom *i for i in c ]
        e = ( fsum (   ( 1.0 * k * k for k in coeffs ) ) ** 0.5 ) / ( denom * sqrt ( 12.0 ) ) 
        if all (  [ float(i) == i for i in coeffs  ] ) :

            
            n       = r['accuracy']
            d       = D 
            stencil = tuple ( range ( - I - D - 1 , 0 ) )

            fun   = lambda x : x**(n+d)
            ff    = fact ( n + d )
            a     = denom * ff
            b     = calc_dot ( fun , 0.0 , 1.0 , stencil , coeffs )
            term  = a / b
            
            to    = topt ( d , n ) 

            line = ' RuleConf ( %d , darray ( %s ) , %s , %d , darray  ( range ( %s , 0 ) ) , %.10g , %.8g  , %.4f ) ,' % ( D , coeffs , denom , n , -I - D - 1 , e , term , to )
            roff = r['offsets']
            
            assert roff == stencil , 'forward invalid stencil!'
            print ( line )            

    
# =============================================================================
## helper function to get table for central differences
def make_central_diffs ( order ) :
    """Helper function to get table for central differences
    """
    
    import findiff
    from   ostap.math.base        import lcm 
    from   math                   import fsum, sqrt 
    from   ostap.math.finitediffs import calc_dot 

    def fact ( n ) :
        return 1 if 0 == n else n * fact ( n - 1 ) 

    
    D = order 
    for I in range ( 0 , 26 ) :
    ## for I in range ( 0 , 10 ) :
        
        r = findiff.coefficients(deriv=D, offsets=tuple ( k for k in range (-I-D, I+D+1 ) ), symbolic=True)
        c = r['coefficients']
        ## print ( 'offsets' , r['offsets'] ) 
        denom = lcm (  *[ i.q for i in c ] )
        
        coeffs = [ denom*i for i in c ]
        e      = ( fsum (   ( 1.0 * k * k for k in coeffs ) ) ** 0.5 ) / ( denom * sqrt ( 12.0 ) )
        
        if all (  [ float(i) == i for i in coeffs ] ) :
            
            n       = r['accuracy']
            d       = D 
            stencil = tuple ( range ( -I - D , I + D + 1 ) )

            fun   = lambda x : x**(n+d)
            ff    = fact ( n + d )
            a     = denom * ff
            b     = calc_dot ( fun , 0.0 , 1.0 , stencil , coeffs )
            term  = a / b
            
            to    = topt ( d , n ) 
            
            line = ' RuleConf ( %d , darray ( %s ) , %s , %d , darray ( range ( %s , %s ) ) , %.10g , %.1f , %.4f ) , ' % ( D , coeffs , denom , n , - I - D , I + D + 1 , e , term , to )
            roff = r['offsets']
            assert roff == stencil , 'forward invalid stencil!'
            print ( line )
            

    
# =============================================================================
if '__main__' == __name__ :
    
    test_derivative_1 ()
    test_derivative_2 ()
    test_derivative_3 ()
    test_derivative_4 ()
    test_derivative_5 ()    
    test_derivative_6 ()
    
    test_central_differences  () 
    test_forward_differences  () 
    test_backward_differences () 

    pass

# =============================================================================
##                                                                      The END 
# =============================================================================
