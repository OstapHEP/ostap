#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
## @file ostap/math/tests/test_math_fourier.py
#  Test module for fourur sums 
# ============================================================================= 
""" Test module for fourier sums 
 - FouinerSum 
 - CosineSum 
 - SineSum 

"""
# ============================================================================= 
from   ostap.core.core        import VE, SE, Ostap
from   ostap.math.models      import f1_draw
from   ostap.utils.root_utils import batch_env
from   ostap.plotting.canvas  import use_canvas 
import ostap.math.derivative  as     D
import ostap.math.integral    as     I 
import ostap.logger.table     as     T
import ROOT, random, math 
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_fourier' )
else                       : logger = getLogger ( __name__                  )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

ranges = [   ( -2 , 2 ) , ( 0 , 3 ) , ( 1 , 5 )  , ( -3 , 1 ) ]

def test_fourier_sum_1 () :
    
    logger = getLogger ( 'test_fourier_sum_1' )

    N    =  7

    for xL , xH in ranges :

        pars = [ random.uniform ( -1 , 1 ) for k in range ( N ) ]

        fs   = Ostap.Math.FourierSum ( pars , xL , xH )
        
        ## python sum 
        def pysum ( x ) :
            result = 0.5 * pars [ 0 ]
            for i , p in enumerate ( pars [ 1 : ]  ) :
                k    = ( 2 + i ) // 2
                kt   = k * fs.t ( x )
                if 0 == i % 2 : result += p * math.sin ( kt )
                else          : result += p * math.cos ( kt ) 
            return result 

        ## numerical derivative
        nd = D.Derivative ( fs )

        x0 = fs.x0()   
        ## numerical integral from x0 
        ni = I.Integral ( fs , xlow = x0 ) 
        
        cnt0 = SE() ##  | value 0 python ] 
        cnt1 = SE()
        cnt2 = SE()
        cnt3 = SE()
        cnt4 = SE()
        cnt5 = SE()

        ## derivative as object 
        fd  = fs.the_derivative ()

        ## integral as object 
        fi  = fs.the_integral   ()

        
        for i in range ( 1000  ) :
            x     = random.uniform ( xL , xH )
            cnt0 += abs ( fs            ( x ) - pysum   ( x ) ) 
            cnt1 += abs ( fs.derivative ( x ) - nd      ( x ) ) 
            cnt2 += abs ( fd            ( x ) - nd      ( x ) )
            
            cnt3 += abs ( fs.integral   ( x ) - ni ( x ) )

            ff = fi ( x ) + 0.5 * fs [ 0 ] * ( x - x0 ) - fi ( x0 ) 
            cnt4 += abs ( ff                  - ni ( x ) )
            
            cnt5 += abs ( ff  - fs.integral ( x )  ) 

                    
        if cnt0.sum() < 1.e-9 : logger.info  ( 'Values        Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt0 ) )
        else                  : logger.error ( 'Values        Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt0 ) )
        if cnt1.sum() < 1.e-9 : logger.info  ( 'Derivatives/1 Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt1 ) )
        else                  : logger.error ( 'Derivatives/1 Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt1 ) )
        if cnt2.sum() < 1.e-9 : logger.info  ( 'Derivatives/2 Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt2 ) )
        else                  : logger.error ( 'Derivatives/2 Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt2 ) ) 
        if cnt3.sum() < 1.e-9 : logger.info  ( 'Integral/1    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt3 ) ) 
        else                  : logger.error ( 'Integral/1    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt3 ) ) 
        if cnt4.sum() < 1.e-9 : logger.info  ( 'Integral/2    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt4 ) )
        else                  : logger.error ( 'Integral/2    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt4 ) ) 
        if cnt5.sum() < 1.e-9 : logger.info  ( 'Integral/3    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt5 ) ) 
        else                  : logger.error ( 'Integral/3    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt5 ) ) 


# ========================================================================================
def test_cosine_sum_1 () :
    
    logger = getLogger ( 'test_cosine_sum_1' )

    N    =  5

    for xL , xH in ranges :

        pars = [ random.uniform ( -1 , 1 ) for k in range ( N ) ]

        cs   = Ostap.Math.CosineSum ( pars , xL , xH )
        
        ## python sum 
        def pysum ( x ) :
            result = 0.5 * pars [ 0 ]
            for i , p in enumerate ( pars [ 1 : ]  ) :
                k   = i + 1 
                kt  = k * cs.t(x) 
                result += p * math.cos ( kt ) 
            return result 

        ## numerical derivative
        nd = D.Derivative ( cs )

        x0 = cs.x0()   
        ## numerical integral from x0 
        ni = I.Integral ( cs , xlow = x0 ) 
        
        cnt0 = SE() ##  | value 0 python ] 
        cnt1 = SE()
        cnt2 = SE()
        cnt3 = SE()
        cnt4 = SE()
        cnt5 = SE()

        ## derivative as object 
        fd  = cs.the_derivative ()

        ## integral as object 
        fi  = cs.the_integral   ()

        for i in range ( 1000 ) :
            x     = random.uniform ( xL , xH )
            cnt0 += abs ( cs            ( x ) - pysum   ( x ) ) 
            cnt1 += abs ( cs.derivative ( x ) - nd      ( x ) ) 
            cnt2 += abs ( fd            ( x ) - nd      ( x ) )
            
            cnt3 += abs ( cs.integral   ( x ) - ni ( x ) )

            ## ff = fi ( x ) + 0.5 * cs [ 0 ] * cs.t ( x ) 
            ff = fi ( x ) + 0.5 * cs [ 0 ] * ( x - x0 ) 
            cnt4 += abs ( ff                      - ni ( x ) )
            
            cnt5 += abs ( ff                      - cs.integral ( x ) ) 

            
        if cnt0.sum() < 1.e-9 : logger.info  ( 'Values        Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt0 ) )
        else                  : logger.error ( 'Values        Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt0 ) )
        if cnt1.sum() < 1.e-9 : logger.info  ( 'Derivatives/1 Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt1 ) )
        else                  : logger.error ( 'Derivatives/1 Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt1 ) )
        if cnt2.sum() < 1.e-9 : logger.info  ( 'Derivatives/2 Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt2 ) )
        else                  : logger.error ( 'Derivatives/2 Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt2 ) ) 
        if cnt3.sum() < 1.e-9 : logger.info  ( 'Integral/1    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt3 ) ) 
        else                  : logger.error ( 'Integral/1    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt3 ) ) 
        if cnt4.sum() < 1.e-9 : logger.info  ( 'Integral/2    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt4 ) )
        else                  : logger.error ( 'Integral/2    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt4 ) ) 
        if cnt5.sum() < 1.e-9 : logger.info  ( 'Integral/3    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt5 ) ) 
        else                  : logger.error ( 'Integral/3    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt5 ) ) 


# ============================================================================
def test_sine_sum_1 () :
    
    logger = getLogger ( 'test_sine_sum_1' )

    N    =  5

    for xL , xH in ranges :

        pars = [ random.uniform ( -1 , 1 ) for k in range ( N ) ]

        ss   = Ostap.Math.SineSum ( pars , xL , xH )
        
        ## python sum 
        def pysum ( x ) :
            result = 0 
            for i , p in enumerate ( pars ) :
                k   = i + 1 
                kt  = k * ss.t(x) 
                result += p * math.sin ( kt ) 
            return result 

        ## numerical derivative
        nd = D.Derivative ( ss )

        x0 = ss.x0()   
        ## numerical integral from x0 
        ni = I.Integral ( ss , xlow = x0 ) 
        
        cnt0 = SE() ##  | value 0 python ] 
        cnt1 = SE()
        cnt2 = SE()
        cnt3 = SE()
        cnt4 = SE()
        cnt5 = SE()

        ## derivative as object 
        fd  = ss.the_derivative ()

        ## integral as object 
        fi  = ss.the_integral ()
        
        for i in range ( 1000 ) :
            x     = random.uniform ( xL , xH )
            cnt0 += abs ( ss            ( x ) - pysum   ( x ) ) 
            cnt1 += abs ( ss.derivative ( x ) - nd      ( x ) ) 
            cnt2 += abs ( fd            ( x ) - nd      ( x ) )
            
            cnt3 += abs ( ss.integral   ( x ) - ni ( x ) )

            ff = fi ( x ) 
            cnt4 += abs ( ff                   - ni ( x ) )
            
            cnt5 += abs ( ff                   - ss.integral ( x ) ) 


        if cnt0.sum() < 1.e-9 : logger.info  ( 'Values        Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt0 ) )
        else                  : logger.error ( 'Values        Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt0 ) )
        if cnt1.sum() < 1.e-9 : logger.info  ( 'Derivatives/1 Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt1 ) )
        else                  : logger.error ( 'Derivatives/1 Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt1 ) )
        if cnt2.sum() < 1.e-9 : logger.info  ( 'Derivatives/2 Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt2 ) )
        else                  : logger.error ( 'Derivatives/2 Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt2 ) ) 
        if cnt3.sum() < 1.e-9 : logger.info  ( 'Integral/1    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt3 ) ) 
        else                  : logger.error ( 'Integral/1    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt3 ) ) 
        if cnt4.sum() < 1.e-9 : logger.info  ( 'Integral/2    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt4 ) )
        else                  : logger.error ( 'Integral/2    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt4 ) ) 
        if cnt5.sum() < 1.e-9 : logger.info  ( 'Integral/3    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt5 ) ) 
        else                  : logger.error ( 'Integral/3    Range: [%+.1f,%+.1f] : diff %s' % ( xL , xH , cnt5 ) ) 


def test_cesaro_sum_1 () :
    
    logger = getLogger ( 'test_cesaro_sum_1' )

    rows = [ ( 'N' , 'K' , 'delta' ) ] 

    def fun0 ( x ) : return 0.5 * x 
               
    for N in ( 5 , 10 , 20 , 50 , 100 , 200 , 500 , 1000 ) :

        pars =  [ (-1)**(k)/(k+1) for k in range ( N + 1 ) ]
        fun  = Ostap.Math.SineSum ( pars , 0 , math.pi ) 

        delta = fun.xmax() - fun.xmin()
        
        with use_canvas ( 'N=%s' % N , wait = 2 ) : 

            fun.draw ( xmin = -math.pi , xmax = 2 * math.pi , color = 2 )
            f1_draw  ( fun0 , 'same' , xmin = fun.xmin () , xmax = fun.xmax() , color = 1 , width = 2 ) 
            
            cesaros = tuple ( fun.cesaro ( i ) for i in range ( 3 ) )

            color = 4 
            for i,c in enumerate ( cesaros ) :

                c.draw ( 'same' , color = 3 )

                df = lambda x :  ( fun0 ( x ) - c ( x ) ) ** 2 

                dd = I.integral ( df , xmin = fun.xmin () , xmax = fun.xmax() - 0.3 * delta ) 

                row = '%s' % N , '%s' % i , '%.4g' % dd 
                rows.append ( row )
                color += 1 
                
                
    title = 'Cesaro summation'
    table = T.table ( rows , title = title , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) ) 

            
            

    
        
# =============================================================================
if '__main__' == __name__ :
    
    test_fourier_sum_1 ()
    test_cosine_sum_1  ()
    test_sine_sum_1    ()
    test_cesaro_sum_1  ()

# =============================================================================
##                                                                      The END 
# =============================================================================
