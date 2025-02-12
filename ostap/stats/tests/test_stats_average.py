#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/stats/tests/test_stats_avergae.py 
# Test averages for inconsistend data 
# Copyright (c) Ostap developpers.
# ============================================================================= 
""" # Test averages for inconsistend data 
"""
# =============================================================================
from   ostap.math.ve         import VE 
from   ostap.plotting.canvas import use_canvas 
from   ostap.stats.average   import *
from   ostap.utils.utils     import batch_env 
import ostap.logger.table    as     T 
import random, math 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'tests_stats_avegrate' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
batch_env ( logger ) 
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import warnings
    with warnings.catch_warnings():
        # =====================================================================
        warnings.simplefilter("ignore")                
        import bayesian_average as ba         
    # =========================================================================
    import numpy            as np
    import warnings  
    def ba_standard     ( *values ) :
        data   = np.array ( [ v.value() for v in values ] )
        sigma  = np.array ( [ v.error() for v in values ] )
        r,e    = ba.average ( data , sigma , mode = 'standard' )
        return VE ( r , e * e ) 

    def ba_birge        ( *values ) :
        data   = np.array ( [ v.value() for v in values ] )
        sigma  = np.array ( [ v.error() for v in values ] )
        r,e    = ba.average ( data , sigma , mode = 'birge' )
        return VE ( r , e * e )
    
    def ba_conservative ( *values ) :
        data   = np.array ( [ v.value() for v in values ] )
        sigma  = np.array ( [ v.error() for v in values ] )
        r,e    = ba.average ( data , sigma , mode = 'cons' )
        return VE ( r , e * e ) 

    def ba_jeffreys     ( *values ) :
        data   = np.array ( [ v.value() for v in values ] )
        sigma  = np.array ( [ v.error() for v in values ] )
        r,e    = ba.average ( data , sigma , mode = 'jeffreys' )
        return VE ( r , e * e ) 
    
    def ba_plot ( *values , **kwargs ) :
        data   = np.array ( [ v.value() for v in values ] )
        sigma  = np.array ( [ v.error() for v in values ] )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")                
            return ba.plot_average ( data , sigma , **kwargs )
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    ba = None 


N     = 20
mu    = 1
sigma = 0.1

data1 = []
data2 = []

for i in range ( N ) :
    d = random.gauss ( mu , sigma )
    b = random.gauss (  0 , 10 * sigma )
    
    p1 = VE ( d     , sigma * sigma )
    p2 = VE ( d + b , sigma * sigma ) ## random bias 
    data1.append ( p1 )
    data2.append ( p2 )


## add singel outlier 
data3 = data1 + [ VE ( mu + 5 * sigma , sigma * sigma / 9 ) ]

## clear bimodal distrbution 
data4 = []
for i in range ( N ) :
    
    d1 = random.gauss (   mu , sigma )
    d2 = random.gauss ( 2*mu , sigma )
    
    p1 = VE ( d1   , sigma * sigma )
    p2 = VE ( d2   , sigma * sigma )

    if 0 == i % 3 : data4.append ( p1 )
    else          : data4.append ( p2 )
    

# ===============================================================================
def test_average1 () :

    logger = getLogger ( "test_average1" )
    
    logger.info ( "Average (in)consistent data" )

    for data, title in [ ( data1 , 'data set #1' ) ,
                         ( data2 , 'data set #2' ) ,
                         ( data3 , 'data set #3' ) ,
                         ( data4 , 'data set #4' ) ] : 
        
        rows = [ ('', 'value', 'scale' ) ]
        
        v1        = standard_average     ( *data )
        row = 'Standard' , v1.toString( '%.3f +/- %-.3f' ) , ''
        rows.append ( row )                                
        
        v2, scale = birge_average        ( *data )
        row = 'Birge/scaled' , v2.toString( '%.3f +/- %-.3f' ) , '%.3f' % scale 
        rows.append ( row )                                
    
        v3        = conservative_average ( *data )
        row = 'Conservative' , v3.toString( '%.3f +/- %-.3f' ) , ''
        rows.append ( row )                                
        
        v4        = jeffreys_average     ( *data )
        row = 'Jeffreys limit' , v4.toString( '%.3f +/- %-.3f' ) , ''
        rows.append ( row )                                
        
        if ba :
            
            v6        = ba_standard    ( *data )
            row = 'BA:Standard' , v6.toString( '%.3f +/- %-.3f' ) , ''
            rows.append ( row )                                
            
            v7       = ba_birge        ( *data )
            row = 'BA:Birge/scaled' , v7.toString( '%.3f +/- %-.3f' ) , ''  
            rows.append ( row )                                
            
            v8       = ba_conservative ( *data )
            row = 'BA:Conservative' , v8.toString( '%.3f +/- %-.3f' ) , ''  
            rows.append ( row )                                
            
            v9        = ba_jeffreys    ( *data )
            row = 'BA:Jeffreys limit' , v9.toString( '%.3f +/- %-.3f' ) , ''
            rows.append ( row )                                
            
        logger.info ( '%s\n%s' % ( title , T.table ( rows , title = title , prefix = '# ' , alignment = 'lcc' ) ) ) 


# ===============================================================================
def test_average2 () :

    logger = getLogger ( "test_average" )
    
    logger.info ( "Average (in)consistent data" )

    for data, title in [ ( data1 , 'data set #1' ) ,
                         ( data2 , 'data set #2' ) ,
                         ( data3 , 'data set #3' ) ,
                         ( data4 , 'data set #4' ) ] : 
        
        rows = [ ('', 'value', 'scale' ) ]

        with use_canvas ( 'S', wait = 1 ) :
            s = Standard ( *data )
            s.graph.draw ( 'ap' )
            
        with use_canvas ( 'Standard: %s' % title      , wait = 1 ) :            
            s = Standard ( *data )
            s.draw_nll   ( linecolor = 2 , linewidth=2 )

        with use_canvas ( 'Conservative: %s' % title  , wait = 1 ) :            
            c = Conservative ( *data ) 
            c.draw_nll   (  linecolor = 4 , linewidth = 2 )

        with use_canvas ( 'Jeffreys %s' % title       , wait = 4 ) :            
            j = Jeffreys ( *data ) 
            j.draw_nll   (  linecolor = 8 , linewidth = 2 )

        if ba :
            ba_plot ( *data , jeffreys_val = True, jeffreys_like = True, plot_data=True )
            ba_plot ( *data , cons_val     = True, cons_like     = True, plot_data=True )
            ba_plot ( *data , standard_val = True, standard_like = True, plot_data=True )
            
            
# ===============================================================================
if '__main__' == __name__ :

    test_average1 ()
    test_average2 ()
    
# ===============================================================================
##                                                                        The END 
# ===============================================================================
                   
                   
                   

            
    
                   
