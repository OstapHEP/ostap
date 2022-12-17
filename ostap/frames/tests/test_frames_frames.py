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
from   ostap.core.pyrouts    import hID 
from   ostap.utils.cleanup   import CleanUp
from   ostap.trees.trees     import Tree
from   ostap.plotting.canvas import use_canvas
from   ostap.utils.utils     import wait 
from   ostap.core.meta_info  import root_info
from   ostap.utils.utils     import implicitMT
from   ostap.utils.timing    import timing
import ROOT, os, sys
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_frames_frames' )
else                       : logger = getLogger ( __name__            )
# ============================================================================= 


if root_info < (6,16) :
    
    logger.warning ( "Tests are disabled for this version of ROOT %s" % str ( root_info ) )
    
else :
    
    from ostap.frames.frames   import DataFrame, frame_project, frame_statVar, frame_statCov, frame_progress 
    
    # A simple helper function to fill a test tree
    def fill_tree ( tname , fname ) :
        
        tdf  = DataFrame        ( 1000 )
        a    = tdf.ProgressBar  ( 1000 )
        tdf.Define   ("one", "1.0"                    )\
                   .Define   ("b1" , "(double) tdfentry_ + 1 ")\
                   .Define   ("b2" , "(1.0+b1)*(1.0+b1)"      )\
                   .Snapshot ( tname, fname )
        
    # We prepare an input tree to run on
    tname  = "myTree"
    tmpdir = CleanUp().tmpdir
    fname  = os.path.join ( tmpdir , 'test_frame.root' )
    
    fill_tree ( tname, fname )


# =============================================================================
def test_frame0 () :
    
    logger = getLogger ( 'test_frame0' ) 
    if root_info < (6,16) : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 
    
    frame = DataFrame ( tname        , fname        )
    tree  = Tree      ( name = tname , file = fname ).chain
    
    
    logger.info ( 80*'*' ) 
    logger.info ( 'MT enabled?  %s' % ROOT.ROOT.IsImplicitMTEnabled() ) 
    
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
    
    
    

# =============================================================================
def test_frame1 ( ) :

    logger = getLogger ( 'test_frame1' ) 
    if root_info < ( 6 , 16 ) : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 
    
    frame = DataFrame ( tname        , fname        )
    tree  = Tree      ( name = tname , file = fname ).chain
    
    h1 = tree .draw ( 'b1' , '1/b1' )
    h2 = frame.draw ( 'b1' , '1/b1' )

    h1 = ROOT.TH1D( hID() , '' , 100 , 0 , 1000 )
    h2 = ROOT.TH1D( hID() , '' , 100 , 0 , 1000 )

    tree.project  ( h1    , 'b1' , '1.0/b1' )
    frame_project ( frame , h2 , 'b1' , '2.0/b1' )
    
    with wait ( 3 ), use_canvas ( 'test_frame1' ) : 
        h1.red  ()
        h2.blue ()    
        h1.draw ()
        h2.draw ( 'same hist' )

# =============================================================================
def test_frame2 ( ) :
 
    logger = getLogger ( 'test_frame2' ) 
    if root_info <  ( 6 , 16 ) : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 
    
    frame = DataFrame     ( tname        , fname        )
    tree  = Tree          ( name = tname , file = fname ).chain

    s1    = tree.statVar  (         'b1' , '1/b1' )
    s2    = frame_statVar ( frame , 'b1' , '1/b1' )
    
    logger.info ('StatTree  :%s' % s1 ) 
    logger.info ('StatFrame :%s' % s2 ) 


# =============================================================================
def test_frame3 () :

    logger = getLogger ( 'test_frame3' ) 
    if root_info < ( 6 , 16 ) : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 

    tree   = Tree      ( name = tname , file = fname ).chain

    for i in range ( 5 ) :
        
        frame = DataFrame ( tname       , fname        )
        pb    = frame_progress ( frame  , len ( tree ) )
        message = 'Value: %d %s' % ( i , pb.GetValue() )
        ## sys.stdout.write('\n')        
        logger.info ( message ) 
        
# =============================================================================
def test_frame4 () :

    logger = getLogger ( 'test_frame4' ) 
    logger = getLogger ( 'test_frame0' ) 
    if root_info < (6,16) : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 
                   
    frame = DataFrame     ( tname        , fname        )
    tree  = Tree      ( name = tname , file = fname ).chain

    pb    = frame_progress ( frame  , len ( tree ) )
    frame = frame.Filter ( 'b1>100' , '#1' )
    frame = frame.Filter ( 'b1>200' , '#2' )
    frame = frame.Filter ( 'b1>300' , '#3' )
    frame = frame.Filter ( 'b1>400' , '#4' )
    frame = frame.Filter ( 'b1>500' , '#5' )
    cnt   = frame.Count  ()
    report = frame.Report()

    from ostap.frames.frames import report_print, report_as_table
    title = 'Filter summary'
    logger.info ( '%s\n%s' % ( title , report_print ( report , title = title , prefix = '# ') ) )
    
# =============================================================================
if '__main__' == __name__ :

    test_frame0 () 
    test_frame1 ()
    test_frame2 ()
    test_frame3 ()
    test_frame4 ()
    

# =============================================================================
##                                                                      The END 
# =============================================================================
