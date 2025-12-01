#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_rootfinder.py
#  Test module for the file ostap/math/rotfinder.py
# ============================================================================= 
""" Test module for ostap/math/rootfinder.py.
"""
# ============================================================================= 
from   ostap.math.rootfinder    import ( find_root       , findroot , 
                                         findroot_scipy  ,
                                         findroot_ostap  ,
                                         findroot_ostap2 )  
from   ostap.logger.pretty      import fmt_pretty_values 
import ostap.logger.table       as     T  
import ostap.math.models 
import random,math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_math_rootfinder' ) 
else                       : logger = getLogger ( __name__               )
# =============================================================================


def test_root_sin () :

    fun  = lambda  x : -1.00 * math.sin ( 0.5 * ( x - 1.0 ) ) 
    der1 = lambda  x : -0.50 * math.cos ( 0.5 * ( x - 1.0 ) ) 
    der2 = lambda  x : +0.25 * math.sin ( 0.5 * ( x - 1.0 ) ) 
    
    logger.info ( 'sin/Halley:  %.15g\n%s' % find_root ( fun , -0.5 , 2.5   ,
                                                         deriv1      = der1 ,
                                                         deriv2      = der2 ,
                                                         full_output = True ) )
    logger.info ( 'sin/Newton:  %.15g\n%s' % find_root ( fun , -0.5 , 2.5   ,
                                                         deriv2      = der1 ,
                                                         full_output = True ) )
    logger.info ( 'sin/Plain :  %.15g\n%s' % find_root ( fun , -0.5 , 2.5   ,
                                                         full_output = True ) )
    logger.info ( 'sin/Brent :  %.15g\n%s' % findroot  ( fun , -0.5 , 2.5   ,
                                                         full_output = True ) )
    
def test_root_mult () :

    K = 1
    N = 2*K + 1 
    
    fun  = lambda  x : 1.0                * ( x - 1.0 ) **   N 
    der1 = lambda  x : 1.0 * N            * ( x - 1.0 ) ** ( N - 1 )
    der2 = lambda  x : 1.0 * N * ( N -1 ) * ( x - 1.0 ) ** ( N - 2 )
    
    logger.info ( 'mult/Halley:  %.15g\n%s' % find_root ( fun , -1 , 10       ,
                                                          deriv1      = der1 ,
                                                          deriv2      = der2 ,
                                                          full_output = True  ,
                                                          disp        = False ) )
    logger.info ( 'mult/Newton:  %.15g\n%s' % find_root ( fun , -1 , 10       ,
                                                          deriv1      = der1 ,
                                                          full_output = True  ,
                                                          disp        = False ) )
    logger.info ( 'mult/Plain :  %.15g\n%s' % find_root ( fun , -1 , 10       ,
                                                          full_output = True  ,
                                                          disp        = False ) )
    logger.info ( 'mult/Brent :  %.15g\n%s' % findroot  ( fun , -1 , 10       ,
                                                          full_output = True  ,
                                                          disp        = False ) )
    

def test_root_test () :
    
    inputs = [ ( 'fun#1' , lambda x : x*2 - ( 1 - x ) ** 5                ,  0.1 ,  1.0   ) ,
               ( 'fun#2' , lambda x : math.cos ( x ) - x ** 3             ,  0.1 ,  1.0   ) , 
               ( 'fun#3' , lambda x : x * math.exp ( x ) - 1              , -1.0 ,  1.0   ) , 
               ( 'fun#4' , lambda x : math.log(x)                         ,  0.5 ,  5.0   ) , 
               ( 'fun#5' , lambda x : x ** 3                              , -0.5 ,  1.0/3 ) , 
               ( 'fun#6' , lambda x : 1.0/x - math.sin(x)  + 1            , -1.3 , -0.5   ) , 
               ( 'fun#7' , lambda x : math.exp ( x * x + 7 * x - 30 ) - 1 ,  2.8 ,  3.1   ) ]
    
    
    rows = [ ( 'Function' , 'scipy'  , 'ostap/py' , 'ostap/C++', '' , '#scipy' , '#ostap/py' , '#ostap/C++') ]
    
    for tag , f , low , high in inputs : 

        r1 , o1 = findroot_scipy  ( f , low , high , full_output = True , disp = False , maxiter = 200 )
        r2 , o2 = findroot_ostap  ( f , low , high , full_output = True , disp = False , maxiter = 200 )
        r3 , o3 = findroot_ostap2 ( f , low , high , full_output = True , disp = False , maxiter = 200 )

        fmtv, expo = fmt_pretty_values ( r1  , r2 , r3 , precision = 14 , width = 16 )
        if expo : 
            scale = 10 ** expo
            row = tag, fmtv % ( r1 / scale ) , fmtv % ( r2 / scale ) , fmtv % ( r3 / scale  ) , '10^%+d' % expo 
        else :  
            row = tag, fmtv % r1 , fmtv % r2 , fmtv % r3 , '' 
     
        n1 = o1.function_calls
        n2 = o2.function_calls + o2.derivative1_calls + o2.derivative2_calls
        n3 = o3.function_calls + o3.derivative1_calls + o3.derivative2_calls
        
        row += (  '%3d' % n1  , '%3d' % n2 ,  '%3d' % n3 )
        rows.append ( row )
        
    rows  = T.remove_empty_columns ( rows )
    title = 'Compare root-finders'
    logger.info ( '%s\n%s' % ( title , T.table ( rows , title = title , prefix = '# ', alignment = 'lccccccccccc')) )
         
        
# =============================================================================
if '__main__' == __name__ :
    
    test_root_sin  () 
    test_root_mult () 
    test_root_test ()
    
    
# =============================================================================
##                                                                      The END 
# =============================================================================
