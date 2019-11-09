#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/frames/tests/test_frames_frames.py
# Test module for ostap/frames/frames.py.
# Copyright (c) Ostap developpers.
# ============================================================================= 
""" Test module for ostap/frames/frames.py.
"""
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_frames_frames' )
else                       : logger = getLogger ( __name__            )
# ============================================================================= 
import ROOT, os 
from ostap.frames.frames import DataFrame
from ostap.utils.cleanup import CleanUp
from ostap.trees.trees   import Tree

tmpdir = CleanUp().tmpdir
fname  = os.path.join ( tmpdir , 'test_frame.root' )

# A simple helper function to fill a test tree
def fill_tree ( tname , fname ):
    
    tdf = DataFrame        ( 1000 )
    a   = tdf.ProgressBar  ( 1000 )
    tdf.Define   ("one", "1.0"               )\
       .Define   ("b1" , "(double) tdfentry_")\
       .Define   ("b2" , "(1.0+b1)*(1.0+b1)" )\
       .Snapshot ( tname, fname )
    
# We prepare an input tree to run on
tname = "myTree"
fill_tree ( tname, fname )

frame = DataFrame ( tname        , fname        )
tree  = Tree      ( name = tname , file = fname ).chain


from ostap.utils.utils  import implicitMT
from ostap.utils.timing import timing

def test_frame0 () :
    
    for mt in  ( False , True ) :
        
        with implicitMT ( mt ) :

            logger.info ( 80*'*' ) 
            logger.info ( 'MT enabled?  %s/%s' % ( ROOT.ROOT.IsImplicitMTEnabled() , ROOT.ROOT.GetImplicitMTPoolSize() ) ) 
            
            logger.info ( 'Tree  :\n%s' % tree  ) 
            logger.info ( 'Frame :\n%s' % frame )
            
            logger.info ( 'Len:                          %30s vs %-30s '  % ( len ( frame )         , len ( tree ) )        )
            logger.info ( 'Branches:                     %30s vs %-30s '  % ( frame.branches()      , tree.branches() )     ) 
            logger.info ( 'nEff:                         %30s vs %-30s '  % ( frame.nEff    ()      , tree.nEff    () )     )
            logger.info ( 'nEff(b1):                     %30s vs %-30s '  % ( frame.nEff    ('b1')  , tree.nEff    ('b1') ) )
            logger.info ( "m(5,50,'b1','b1/(b2+1'):      %30s vs %-30s "  % ( frame.get_moment ( 5 , 50 ,  'b1' , 'b1/(b2+1)' ) ,
                                                                              tree .get_moment ( 5 , 50 ,  'b1' , 'b1/(b2+1)' ) ) ) 
            logger.info ( "m(5,'b1','b1/(b2+1'):         %30s vs %-30s "  % ( frame.moment ( 5 , 'b1' , 'b1/(b2+1)' ) ,
                                                                              tree .moment ( 5 , 'b1' , 'b1/(b2+1)' ) ) ) 
            logger.info ( "cm(5,'b1','b1/(b2+1'):        %30s vs %-30s "  % ( frame.central_moment ( 2 , 'b1' , 'b1/(b2+1)' ) ,
                                                                              tree .central_moment ( 2 , 'b1' , 'b1/(b2+1)' ) ) ) 
            logger.info ( "mean('b1','b1/(b2+1'):        %30s vs %-30s "  % ( frame.mean ( 'b1' , 'b1/(b2+1)' ) ,
                                                                              tree .mean ( 'b1' , 'b1/(b2+1)' ) ) ) 
            logger.info ( "rms ('b1','b1/(b2+1'):        %30s vs %-30s "  % ( frame.rms  ( 'b1' , 'b1/(b2+1)' ) ,
                                                                              tree .rms  ( 'b1' , 'b1/(b2+1)' ) ) ) 
            logger.info ( "skew('b1','b1/(b2+1'):        %30s vs %-30s "  % ( frame.skewness ( 'b1' , 'b1/(b2+1)' ) ,
                                                                              tree .skewness ( 'b1' , 'b1/(b2+1)' ) ) )
            logger.info ( "kurt('b1','b1/(b2+1'):        %30s vs %-30s "  % ( frame.kurtosis ( 'b1' , 'b1/(b2+1)' ) ,
                                                                              tree .kurtosis ( 'b1' , 'b1/(b2+1)' ) ) ) 
            logger.info ( "quantile(0.3, 'b1','b1<500'): %30s vs %-30s "  % ( frame.quantile ( 0.3 , 'b1' , 'b1/(b2+1)' ) ,
                                                                              tree .quantile ( 0.3 , 'b1' , 'b1/(b2+1)' ) ) )
            logger.info ( "median(0.3, 'b1','b1<500'):   %30s vs %-30s "  % ( frame.median   (       'b1' , 'b1/(b2+1)' ) ,
                                                                              tree .median   (       'b1' , 'b1/(b2+1)' ) ) )
            logger.info ( "terciles('b1','b1<500'):      %30s vs %-30s "  % ( frame.terciles (       'b1' , 'b1/(b2+1)' ) ,
                                                                              tree .terciles (       'b1' , 'b1/(b2+1)' ) ) )
            

def test_frame1 () :
    c = 0 
    for mt in  ( False , True , False , True ) :
        with implicitMT ( mt ) :
            logger.info ( 'MT enabled?  %s/%s' % ( ROOT.ROOT.IsImplicitMTEnabled() , ROOT.ROOT.GetImplicitMTPoolSize() ) ) 
            for obj in ( (tree,'Tree') , (frame,'Frame') , ( tree , 'Tree') ) :            
                with timing( obj[1] , logger ) :
                ##logger.info ('Kurtosis: %s ' % obj[0].kurtosis ( 'b1' , 'b1/(b2+1)' ) ) 
                    c += obj[0].kurtosis ( 'b1' , 'b1/(b2+1)' ) 
                                                          

def test_frame2 ( ) :

    h1 = tree .draw('b1','1/b1')
    h2 = frame.draw('b1','1/b1')
    
        
# =============================================================================
if '__main__' == __name__ :
    
    test_frame0 () 
    test_frame1 ()
    test_frame2 ()
    
    
# =============================================================================
# The END 
# =============================================================================
