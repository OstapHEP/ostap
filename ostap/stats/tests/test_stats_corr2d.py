#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_stats_corr2d.py
# Test module for simple decorrelation
# ============================================================================= 
""" Test module for simple decorrelation 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   builtins                 import range
from   ostap.core.core          import WSE 
from   ostap.utils.cleanup      import CleanUp
from   ostap.utils.timing       import timing
from   ostap.stats.corr2d       import Corr2D
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.progress_bar import progress_bar
import ostap.trees.trees
import ostap.histos.histos 
import ROOT, os, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_stats_corr2d' )
else : 
    logger = getLogger ( __name__             )
# =============================================================================

testdata  = CleanUp.tempfile ( suffix = '.root' , prefix ='ostap-test-tools-reweight2-' )

## prepare soem random data

N1 = 10000 

def prepare_data ( ) : 
    #
    seed =  1234567890 
    random.seed ( seed ) 
    logger.info ( 'Test *RANDOM* data will be generated/seed=%s' % seed  )
    
    with ROOT.TFile.Open ( testdata ,'recreate') as test_file:
        
        test_file.cd() 
        
        tree  = ROOT.TTree ( 'DATA', 'data-tree' )
        tree .SetDirectory ( test_file ) 
        
        from array import array 
        xvar = array    ( 'f' ,  [ 0.0 ] )
        yvar = array    ( 'f' ,  [ 0.0 ] )
        zvar = array    ( 'f' ,  [ 0.0 ] )
        
        tree.Branch ( 'x' , xvar , 'x/F' )
        tree.Branch ( 'y' , yvar , 'y/F' )
        tree.Branch ( 'z' , zvar , 'z/F' )
    
        for i in  range ( N1 ) :
            
            x1 = random.gauss   (  3 , 0.5 )
            x2 = random.gauss   ( 10 , 5.0 )
            x3 = random.uniform ( -1 , 1   )
            
            xvar [ 0 ] = x1
            yvar [ 0 ] = x1 + x2  
            zvar [ 0 ] = x1 + x2 + x3   
        
            tree.Fill()

        tree.Write()
        test_file.ls()

if not os.path.exists( testdata ) :
    with timing ( "Prepare input data" , logger = logger ) :
        prepare_data ()

# ============================================================================
## test Corr2Dstuff 
def test_corr2d  ( var1 = 'x' , var2 = 'y') :
    """ Test Corr2D stuff """
    
    logger = getLogger ( 'test_corr2d(%s,%s)' % ( var1 , var2 )  )
    
    logger.info ( 'Test Corr2D(%s,%s)' % ( var1 , var2 )  )
    
    chain = ROOT.TChain ( 'DATA' )
    chain.Add ( testdata )
    
    c2d = Corr2D ( chain , var1 , var2 )

    with use_canvas ( 'test_corr2D(%s,%s)'              %  ( var1 , var2 ) , wait = 1 ) :
        chain.draw ( [ var1 , var2  ] )
    with use_canvas ( 'test_corr2D(%s,%s) decorrelated' %  ( var1 , var2 ) , wait = 1 ) :
        chain.draw ( [ c2d.decorrelated1 , c2d.decorrelated2 ] )
    with use_canvas ( 'test_corr2D(%s,%s) normalized'   %  ( var1 , var2 ) , wait = 1 ) :
        chain.draw ( [ c2d.normalized1   , c2d.normalized2   ] )
    with use_canvas ( 'test_corr2D(%s,%s) uniform'      %  ( var1 , var2 ) , wait = 1 ) :
        chain.draw ( [ c2d.uniform1      , c2d.uniform2      ] )


    logger.info ( 'var1          %s' % chain.statVar ( var1 ) )
    logger.info ( 'var2          %s' % chain.statVar ( var2 ) )
    logger.info ( 'decorrelatd11 %s' % chain.statVar ( c2d.decorrelated1 ) )
    logger.info ( 'decorrelated2 %s' % chain.statVar ( c2d.decorrelated2 ) )
    logger.info ( 'normalized1   %s' % chain.statVar ( c2d.normalized1   ) )
    logger.info ( 'normalized2   %s' % chain.statVar ( c2d.normalized2   ) )
    logger.info ( 'uniform11     %s' % chain.statVar ( c2d.uniform1      ) )
    logger.info ( 'uniform2      %s' % chain.statVar ( c2d.uniform2      ) )
        
    cnt1  = chain.statVar ( c2d.var1          )
    cnt1d = chain.statVar ( c2d.decorrelated1 )
    cnt1n = chain.statVar ( c2d.normalized1   )
    cnt1u = chain.statVar ( c2d.uniform1      )
    
    cnt2  = WSE ()
    cnt2d = WSE ()
    cnt2n = WSE ()
    cnt2u = WSE ()
    
    for i in progress_bar ( chain ) :

        x, y = getattr ( i , var1 ) , getattr ( i , var2 ) 
        
        d1 , _ = c2d.decorrelated ( x , y ) 
        cnt2d += d1
        
        n1 , _ = c2d.normalized   ( x , y ) 
        cnt2n += n1

        u1 , _ = c2d.uniform      ( x , y ) 
        cnt2u += u1

    logger.info ( 'cnt1d: %s' % cnt1d )
    logger.info ( 'cnt2d: %s' % cnt2d )
    logger.info ( 'cnt1n: %s' % cnt1n )
    logger.info ( 'cnt2n: %s' % cnt2n )
    logger.info ( 'cnt1u: %s' % cnt1u )
    logger.info ( 'cnt2u: %s' % cnt2u )
    
# ============================================================================
if '__main__' == __name__ :

    
    test_corr2d ( 'x' , 'y' )
    test_corr2d ( 'x' , 'z' )
    test_corr2d ( 'z' , 'y' )
    

# =============================================================================
##                                                                      The END 
# =============================================================================
