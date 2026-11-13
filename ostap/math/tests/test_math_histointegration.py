#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# =============================================================================
## @file ostap/math/tests/test_math_histointegration.py
# ============================================================================= 
""" Test module
"""
# ============================================================================= 
from   ostap.core.core          import Ostap, hID, SE 
from   ostap.utils.progress_bar import progress_bar
from   ostap.utils.timing       import timing
from   ostap.logger.pretty      import format_pow10 
from   ostap.utils.root_utils   import batch_env
from   ostap.utils.utils4tests  import light_tests, active_test 
import ostap.math.integral      as     I 
import ostap.math.integrator    as     II
import ostap.logger.table       as     T 
import ostap.math.models
import ostap.histos.histos
import ROOT, random 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
logger = getLogger ( 'ostap.test_math_histo_integration' ) 
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
light , weak = light_tests ()
random_freq  = 0.25

# =============================================================================
if   light : logger.warning ( "CONTINUOUS_INTEGRATION: run *LIGHT* version of tests" ) 
elif weak  : logger.warning ( "Machine is not powerful enough: reduce some tests" )

def test_histo1_integration () :

    logger = getLogger ( 'test_histo1_integration' )

    cnts = {}
    rng = ( 1 , 2 )  if ( light or weak ) else ( 0 , 4 ) 

    for k in range ( *rng ) :
        for e in ( True , False ) : 
            cnts [ k, e ]  = SE()

    N          = 1000 
    NB         = 20  
    
    xmin, xmax = -3 , 3
    
    h1    = ROOT.TH1D ( hID() , '' , NB , xmin , xmax )
    NTOYS = 10 if ( light or weak ) else 100 
    for k in progress_bar ( NTOYS ) :

        h1.reset () 
        for i in range ( N ) :
            x = random.gauss(0,1)
            while not xmin < x < xmax : x = random.gauss ( 0 , 1 )                
            h1.Fill ( x ) 

        ii = N * 1.0 * ( xmax - xmin ) / NB 
        for key in cnts :
            
            cnt   = cnts [ key ]
            hh    = Ostap.Math.Histo1D ( h1 , *key ) 
            r1    = hh.integral ()
            cnt  += ( r1 - ii ) / ii 
            
    rows = [ ('Index' ,
              'mean [%s]'    % format_pow10 ( -3 ) ,
              'rms [%s]'     % format_pow10 ( -3 ) ,
              'min/max [%s]' % format_pow10 ( -3 ) ) ]
    
    for key in sorted ( cnts ) :
        cnt = cnts [ key ]

        v    = cnt.mean() * 1.e+3
        mean = '%+.4g +/- %-.4g' % ( v.value() , v.error() ) 
        v    = cnt.rms()  * 1.e+3 
        rms  = '% .4g' % v 

        v1   = cnt.min() * 1.e+3 
        v2   = cnt.max() * 1.e+3 
        mnmx = '%+.4g/%-+4g' % ( v1 , v2 )
        
        row  = str ( key ) , mean , rms  , mnmx

        rows.append ( row ) 

    title = '1D integration'
    table = T.table ( rows , prefix = '# ' , title = title , alignment = 'cccc' )
    logger.info ( '%s:\n%s' % ( title , table ) ) 
    

# =============================================================================
def test_histo2_integration () :
    
    logger = getLogger ( 'test_histo2_integration' )

    cnts = {} 

    rng = ( 1 , 2 )  if ( light or weak ) else ( 0 , 4 ) 

    for i in range ( *rng ) :
        for j in range ( *rng ) :
            for e in ( True , False ) : 
                cnts  [ i,j,e ]  = SE()
                
    N          = 1000 if ( light or weak ) else 5000
    NX         = 10   if ( light or weak ) else 20 
    NY         = 10   if ( light or weak ) else 20
    
    xmin, xmax = -3 , 3 
    ymin, ymax = -3 , 3 

    h2  = ROOT.TH2D ( hID() , '' , NX , xmin , xmax  , NY , ymin , ymax )
    NTOYS      = 10 if ( light or weak ) else 50
    for k in progress_bar ( NTOYS ) :
        
        h2.reset()
        for i in range  ( N  ) :
            x = random.gauss(0,1)
            while not xmin < x < xmax : x = random.gauss(0,1)                
            y = random.gauss(0,1)
            while not ymin < y < ymax : y = random.gauss(0,1)                            
            h2.Fill ( x , y ) 
        
        ii = N * 1.0 * ( xmax - xmin ) * ( ymax - ymin ) / ( NX * NY ) 
            
        for key in cnts :            
            cnt   = cnts [ key ]
            hh    = Ostap.Math.Histo2D ( h2 , *key ) 
            r1    = hh.integral ()
            cnt  += ( r1 - ii ) / ii 
            
    rows = [ ('Index' ,
              'mean [%s]'    % format_pow10 ( -3 ) ,
              'rms [%s]'     % format_pow10 ( -3 ) ,
              'min/max [%s]' % format_pow10 ( -3 ) ) ]
    for key in sorted ( cnts ) :
        cnt = cnts [ key ]

        v    = cnt.mean() * 1.e+3
        mean = '%+.4g +/- %-.4g' % ( v.value() , v.error() ) 
        v    = cnt.rms()  * 1.e+3 
        rms  = '% .4g' % v 

        v1   = cnt.min() * 1.e+3 
        v2   = cnt.max() * 1.e+3 
        mnmx = '%+.4g/%-+4g' % ( v1 , v2 )
        
        row  = str ( key ) , mean , rms  , mnmx

        rows.append ( row ) 

    title = '2D integration'
    table = T.table ( rows , prefix = '# ' , title = title , alignment = 'cccc' )
    logger.info ( '%s:\n%s' % ( title , table ) ) 

# =============================================================================
def test_histo3_integration () :
    
    logger = getLogger ( 'test_histo3_integration' )

    cnts = {}

    rng = ( 1 , 2 )  if ( light or weak ) else ( 0 , 4 ) 
    
    for i in range ( *rng ) :
        for j in range ( *rng  ) :
            for k in range ( *rng  ) :
                for e in ( True , False ) : 
                    cnts [ i,j,k,e ]  = SE()
                
    N          = 2000 if ( light or weak ) else 10000
    NX         = 10   if ( light or weak ) else 20
    NY         = 10   if ( light or weak ) else 20
    NZ         = 10   if ( light or weak ) else 20
    
    xmin, xmax = -3 , 3 
    ymin, ymax = -3 , 3 
    zmin, zmax = -3 , 3 

    h3  = ROOT.TH3D ( hID() , '' ,
                      NX , xmin , xmax  ,
                      NY , ymin , ymax  ,
                      NZ , zmin , zmax )
    
    NTOYS      = 5 if ( light or weak ) else 10    
    for k in progress_bar ( NTOYS ) :
        
        h3.reset()
        for i in range  ( N ) :
            x = random.gauss(0,1)
            while not xmin < x < xmax : x = random.gauss(0,1)                
            y = random.gauss(0,1)
            while not ymin < y < ymax : y = random.gauss(0,1)                            
            z = random.gauss(0,1)
            while not zmin < z < zmax : z = random.gauss(0,1)                            
            h3.Fill ( x , y , z ) 
        
        ii = N * 1.0 * ( xmax - xmin ) * ( ymax - ymin ) *  ( zmax - zmin ) / ( NX * NY * NZ ) 

        for key in cnts :            
            cnt   = cnts [ key ]
            hh    = Ostap.Math.Histo3D ( h3 , *key ) 
            r1    = hh.integral ()
            cnt  += ( r1 - ii ) / ii 

    rows = [ ('Index' ,
              'mean [%s]'    % format_pow10 ( -3 ) ,
              'rms [%s]'     % format_pow10 ( -3 ) ,
              'min/max [%s]' % format_pow10 ( -3 ) ) ]
    for key in sorted ( cnts ) :
        cnt = cnts [ key ]

        v    = cnt.mean() * 1.e+3
        mean = '%+.4g +/- %-.4g' % ( v.value() , v.error() ) 
        v    = cnt.rms()  * 1.e+3 
        rms  = '% .4g' % v 

        v1   = cnt.min() * 1.e+3 
        v2   = cnt.max() * 1.e+3 
        mnmx = '%+.4g/%-+4g' % ( v1 , v2 )
        
        row  = str ( key ) , mean , rms  , mnmx

        rows.append ( row ) 

    title = '3D integration'
    table = T.table ( rows , prefix = '# ' , title = title , alignment = 'cccc' )
    logger.info ( '%s:\n%s' % ( title , table ) ) 
    
# =============================================================================
if '__main__' == __name__ :

    with timing  ( 'histo1_integration' , logger = logger ) :
        if active_test ( random_freq ) : test_histo1_integration()
        
    with timing  ( 'histo2_integration' , logger = logger ) : 
        if active_test ( random_freq ) : test_histo2_integration()

    with timing  ( 'histo3_integration' , logger = logger ) : 
        if active_test ( random_freq ) : test_histo3_integration()

# =============================================================================
##                                                                      The END 
# =============================================================================
