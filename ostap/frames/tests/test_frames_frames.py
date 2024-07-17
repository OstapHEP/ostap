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
from   ostap.core.pyrouts    import hID, Ostap  
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
    
    from ostap.frames.frames   import DataFrame, frame_project, \
         frame_statVar, frame_statCov, frame_progress, has_std_move  
    
    # A simple helper function to fill a test tree
    def fill_tree ( tname , fname ) :
        
        tdf  = DataFrame        ( 10000 )
        a    = tdf.ProgressBar  ( 10000 )
        tdf.Define   ("one", "1.0"                    )\
                   .Define   ("b1" , "(double) tdfentry_ + 1 ") \
                   .Define   ("b2" , "(1.0+b1)*(1.0+b1)"      ) \
                   .Define   ("b3" , "(1.0+b1)*(1.0+b1)"      ) \
                   .Define   ("b4" , "gRandom->Gaus()"        ) \
                   .Define   ("b5" , "gRandom->Gaus()"        ) \
                   .Define   ("b6" , "gRandom->Gaus()"        ) \
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

    s1    = tree.statVar  (         'b1' ) ## , '1/b1' )
    s2    = frame_statVar ( frame , 'b1' ) ## , '1/b1' )

    print ( 'FIRST s1/s2', type ( s1 ) , type ( s2 ) )
    
    logger.info ('StatTree  :%s' % s1 ) 
    logger.info ('StatFrame :%s' % s2 ) 

    s1    = tree.statVar  (         'b1' , '1/b1' )
    s2    = frame_statVar ( frame , 'b1' , '1/b1' )

    print ( 'SECONF s1/s2', type ( s1 ) , type ( s2 ) )
    
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
                   
    frame  = DataFrame     ( tname        , fname        )
    tree   = Tree      ( name = tname , file = fname ).chain

    pb     = frame_progress ( frame  , len ( tree ) )
    frame  = frame.Filter ( 'b1>100' , '#1' )
    frame  = frame.Filter ( 'b1>200' , '#2' )
    frame  = frame.Filter ( 'b1>300' , '#3' )
    frame  = frame.Filter ( 'b1>400' , '#4' )
    frame  = frame.Filter ( 'b1>500' , '#5' )
    cnt    = frame.Count  ()
    report = frame.Report()

    from ostap.frames.frames import report_print, report_as_table
    title = 'Filter summary'
    logger.info ( '%s\n%s' % ( title , report_print ( report , title = title , prefix = '# ') ) )

# =============================================================================
def test_frame5 () :

    logger = getLogger ( 'test_frame5' ) 
    if root_info < (6,16) : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 
                   
    if not has_std_move : 
        logger.warning ( "Test is disabled (no std::move)" )
        return     
    
    frame = DataFrame      ( tname        , fname        )
    tree  = Tree           ( name = tname , file = fname ).chain
    
    from ostap.frames.frames import frame_table

    ft    = frame_table ( frame )

    logger.info ( 'The frame:\n%s' % ft ) 

# =============================================================================
def test_frame6 () :

    logger = getLogger ( 'test_frame6' ) 
    if root_info < (6,16) : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 
                   
    frame = DataFrame      ( tname        , fname        )
    tree  = Tree           ( name = tname , file = fname ).chain
    
    pb    = frame_progress ( frame  , len ( tree ) )
    
    bs = Ostap.Math.Bernstein    ( 10 , -3 , 3 )
    ls = Ostap.Math.LegendreSum  ( 10 , -3 , 3 )
    cs = Ostap.Math.ChebyshevSum ( 10 , -3 , 3 )

    from ostap.frames.frames import frame_project

    h1 = ROOT.TH1D ( hID() , '' , 20 , -3, 3 )
    
    rh = frame_project  ( frame , h1 , 'b6' )

    with use_canvas ( 'test_frame6: b4' , wait = 2 ) :
        h1.draw () 


    ## report = frame.Report()
    ## from ostap.frames.frames import report_print, report_as_table
    ## title = 'Params summary'
    ## logger.info ( '%s\n%s' % ( title , report_print ( report , title = title , prefix = '# ') ) )


# =============================================================================
def test_frame7 () :

    logger = getLogger ( 'test_frame7' ) 
    if root_info < (6,16) or not has_std_move : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return
    
    if not has_std_move : 
        logger.warning ( "Test is disabled (no std::move)" )
        return     
    
    frame = DataFrame      ( tname        , fname        )
    tree  = Tree           ( name = tname , file = fname ).chain
    
    pb    = frame_progress ( frame  , len ( tree ) )

    
    bs = Ostap.Math.Bernstein    ( 10 , -3 , 3 )
    ls = Ostap.Math.LegendreSum  ( 10 , -3 , 3 )
    cs = Ostap.Math.ChebyshevSum ( 10 , -3 , 3 )

    from ostap.frames.frames import frame_param


    rb = frame_param ( frame , bs , 'b4' , '(b4>-1)*1.01' )
    rl = frame_param ( frame , ls , 'b5' )
    rc = frame_param ( frame , cs , 'b6' )

    polc = rc.GetValue()
    poll = rl.GetValue()
    polb = rb.GetValue()

    with use_canvas ( 'test_frame5: b4,b5.b6 as polynomials' , wait = 2 ) :
        poll.draw (          linecolor = 2 ) 
        polc.draw ( 'same' , linecolor = 4 ) 
        polb.draw ( 'same' , linecolor = 8 ) 


# =============================================================================
def test_frame8 () :

    logger = getLogger ( 'test_frame8' ) 
    if root_info < (6,16) or not has_std_move : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return
    
    if not has_std_move : 
        logger.warning ( "Test is disabled (no std::move)" )
        return     
    
    frame = DataFrame      ( tname        , fname        )
    tree  = Tree           ( name = tname , file = fname ).chain
    
    pb    = frame_progress ( frame  , len ( tree ) )

    from ostap.frames.frames   import frame_the_moment

    m1    = frame_the_moment ( frame , 20 , 'b6'          )  ## simple 
    m2    = frame_the_moment ( frame , 20 , 'b6' , '0<b4' )  ## cut 
    m3    = frame_the_moment ( frame , 20 , 'b6' , 'b4'   )  ## weight 

    mm1 = m1.GetValue()
    mm2 = m2.GetValue()
    mm3 = m3.GetValue()

    title = 'Simple moments'
    logger.info ( '%s:\n%s' % ( title , mm1.table ( title = title , prefix = '# ' ) ) ) 
    title = 'moments with cuts'
    logger.info ( '%s:\n%s' % ( title , mm2.table ( title = title , prefix = '# ' ) ) ) 
    title = 'moments with weights'
    logger.info ( '%s:\n%s' % ( title , mm3.table ( title = title , prefix = '# ' ) ) ) 
    
# =============================================================================
if '__main__' == __name__ :

    if (6,16) <= root_info :
        
        ## test_frame0 () 
        ## test_frame1 ()
        test_frame2 ()
        test_frame3 ()
        test_frame4 ()   
        test_frame5 ()   
        test_frame6 ()
        test_frame7 ()
        test_frame8 ()
        
    else :
        
        logger.warning ( "All tests are disabled for this version of ROOT %s" % str ( root_info ) )


# =============================================================================
##                                                                      The END 
# =============================================================================
