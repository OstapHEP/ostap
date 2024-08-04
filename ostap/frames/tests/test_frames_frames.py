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
import ostap.logger.table    as     T 
import ROOT, os, sys
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_frames_frames' )
else                       : logger = getLogger ( __name__            )
# ============================================================================= 
from ostap.frames.frames   import Frames_OK 

if not Frames_OK : 
    
    logger.warning ( "Tests are disabled for this version of ROOT %s" % str ( root_info ) )
    
else :
    
    from ostap.frames.frames   import * 
    
    # A simple helper function to fill a test tree
    def fill_tree ( tname , fname ) :

        N = 10000 
        tdf      = DataFrame      ( N       )
        tdf , pb = frame_progress ( tdf , N ) 
        tdf.Define   ("one", "1.0"                    )\
                   .Define   ("b1" , "tdfentry_ + 1000.0") \
                   .Define   ("b2" , "(1.0+b1)*(1.0+b1)"      ) \
                   .Define   ("b3" , "(1.0+b1)*(1.0+b1)"      ) \
                   .Define   ("b4" , "gRandom->Gaus()"        ) \
                   .Define   ("b5" , "gRandom->Gaus()"        ) \
                   .Define   ("b6" , "10+0.5*gRandom->Gaus()" ) \
                   .Snapshot ( tname, fname )
        
    # We prepare an input tree to run on
    tname  = "myTree"
    tmpdir = CleanUp().tmpdir
    fname  = os.path.join ( tmpdir , 'test_frame.root' )
    
    fill_tree ( tname, fname )


# =============================================================================
def test_frame0 () : 

    logger = getLogger ( 'test_frame0' ) 
    if not Frames_OK : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 
    
    frame = DataFrame ( tname        , fname        )
    tree  = Tree      ( name = tname , file = fname ).chain
    
    logger.info ( 80*'*' ) 
    logger.info ( 'MT enabled?  %s' % ROOT.ROOT.IsImplicitMTEnabled() ) 

    logger.info ( 'Tree  :\n%s' % tree  ) 
    logger.info ( 'Frame :\n%s' % frame )

    rows  = [ ( 'Parameter' , 'Tree' , 'Frame' ) ] 

    row = 'Lengh' , '%d' % len ( tree ) , '%d' % len ( frame )
    rows.append ( row )

    row = 'Branches' , str ( tree.branches() ) , str ( frame.branches() ) 
    rows.append ( row )
    
    row = 'nEff' , '%.2f' % tree.nEff()  , '%.2f' % frame.nEff() 
    rows.append ( row )
    
    row = 'nEff(b6)' , '%.2f' % tree.nEff('b6')  , '%.2f' % frame.nEff('b6')
    rows.append ( row )

    row = 'mean      (b6)' , \
        tree .mean ( 'b6' ).toString ( '%+.2f +/- %-.2f' ) , \
        frame.mean ( 'b6' ).toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )

    row = 'rms       (b6)' , \
        tree .rms  ( 'b6' ).toString ( '%+.2f +/- %-.2f' ) , \
        frame.rms  ( 'b6' ).toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )

    row = 'variance  (b6)' , \
        tree .variance ( 'b6' ).toString ( '%+.2f +/- %-.2f' ) , \
        frame.variance ( 'b6' ).toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )

    row = 'mean      (b6,1)' , \
        tree .mean ( 'b6', '1' ).toString ( '%+.2f +/- %-.2f' ) , \
        frame.mean ( 'b6', '1' ).toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )

    row = 'rms       (b6,1)' , \
        tree .rms  ( 'b6' , '1' ).toString ( '%+.2f +/- %-.2f' ) , \
        frame.rms  ( 'b6' , '1' ).toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )

    row = 'variance  (b6,1)' , \
        tree .variance ( 'b6' , '1' ).toString ( '%+.2f +/- %-.2f' ) , \
        frame.variance ( 'b6' , '1' ).toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )
    
    row = 'skewness  (b6,1)' , \
        tree .skewness ( 'b6' , '1' ).toString ( '%+.2f +/- %-.2f' ) , \
        frame.skewness ( 'b6' , '1' ).toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )

    row = 'kurtosis  (b6,1)' , \
        tree .kurtosis ( 'b6' , '1' ).toString ( '%+.2f +/- %-.2f' ) , \
        frame.kurtosis ( 'b6' , '1' ).toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )

    row = 'mean      (b6,1/b6)' , \
        tree .mean ( 'b6' , '1/b6' ).toString ( '%+.2f +/- %-.2f' ) , \
        frame.mean ( 'b6' , '1/b6' ).toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )

    row = 'rms       (b6,1/b6)' , \
        tree .rms  ( 'b6' , '1/b6' ).toString ( '%+.2f +/- %-.2f' ) , \
        frame.rms  ( 'b6' , '1/b6' ).toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )

    row = 'variance  (b6,1/b6)' , \
        tree .variance ( 'b6' , '1/b6' ).toString ( '%+.2f +/- %-.2f' ) , \
        frame.variance ( 'b6' , '1/b6' ).toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )
    
    row = 'skewness  (b6,1/b6)' , \
        tree .skewness ( 'b6' , '1/b6' ).toString ( '%+.2f +/- %-.2f' ) , \
        frame.skewness ( 'b6' , '1/b6' ).toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )

    row = 'kurtosis  (b6,1/b6)' , \
        tree .kurtosis ( 'b6' , '1/b6' ).toString ( '%+.2f +/- %-.2f' ) , \
        frame.kurtosis ( 'b6' , '1/b6' ).toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )
    
    row = 'quantile  (0.3,b6,1/b6)' , \
        '%s' % tree .quantile ( 0.3 , 'b6' , '1/b6' ) ,\
        '%s' % frame.quantile ( 0.3 , 'b6' , '1/b6' )
    rows.append ( row )

    row = 'terciles  (b6,1/b6)' , \
        '%s' % tree .terciles ( 'b6' , '1/b6' ) ,\
        '%s' % frame.terciles ( 'b6' , '1/b6' )
    rows.append ( row )
    
    row = 'median    (b6,1/b6)' , \
        '%s' % tree .median  ( 'b6' , '1/b6' ) ,\
        '%s' % frame.median  ( 'b6' , '1/b6' )
    rows.append ( row )
    
    title = 'Frame/Tree statistics' 
    logger.info ( '%s\n%s' % ( title , T.table ( rows , title = title , prefix = '# ' , alignment = 'lrl' ) ) ) 
    

# =============================================================================
def test_frame1 ( ) :

    logger = getLogger ( 'test_frame1' ) 
    if not Frames_OK : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 
    
    logger.info ( 80*'*' ) 
    logger.info ( 'MT enabled?  %s' % ROOT.ROOT.IsImplicitMTEnabled() ) 
    
    frame = DataFrame ( tname        , fname        )
    tree  = Tree      ( name = tname , file = fname ).chain
    
    with use_canvas ( 'test_frame1/draw'   , wait = 2 ) : 
        hh1 = tree .draw ( 'b1' , '1.0/b1' , color = 2                                       )
        hh2 = frame.draw ( 'b1' , '1.0/b1' , color = 4 , opts  = 'same hist' , report = True )

    
    xmnmx1 = hh1.xminmax()
    xmnmx2 = hh2.xminmax()
    
    logger.info ( "Histo-draw/tree  : #bins=%3d min/max %+.3f/%+.3f " % ( hh1.nbins() , xmnmx1[0] , xmnmx1[1] ) )
    logger.info ( "Histo-draw/frame : #bins=%3d min/max %+.3f/%+.3f " % ( hh2.nbins() , xmnmx2[0] , xmnmx2[1] ) )
                      
    h1 = ROOT.TH1D ( hID() , '' , 120 , 0 , 12000 )
    h2 = ROOT.TH1D ( hID() , '' , 120 , 0 , 12000 )

    tree.project  (         h1 , 'b1' , '1.0/b1' )
    frame_project ( frame , h2 , 'b1' , '1.0/b1' , report = True , lazy = False )
    
    with use_canvas ( 'test_frame1/project' , wait = 2 ) : 
        h1.red  ()
        h2.blue ()    
        h1.draw ()
        h2.draw ( 'same hist' )

        
# =============================================================================
def test_frame2 ( ) :
 
    logger = getLogger ( 'test_frame2' ) 
    if not Frames_OK : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 
      
    logger.info ( 80*'*' ) 
    logger.info ( 'MT enabled?  %s' % ROOT.ROOT.IsImplicitMTEnabled() ) 
    
    frame = DataFrame     ( tname        , fname        )
    tree  = Tree          ( name = tname , file = fname ).chain

    rows = [ ( 'Statistics' , 'Tree' , 'Frame' ) ]

    s1    = tree.statVar  (         'b1' ) 
    s2    = frame_statVar ( frame , 'b1' ) 

    row = 'StatVar    mean (b1)' , s1.mean().toString ( '%+.2f +/- %-.2f' ) , s2.mean().toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )

    s1    = tree.statVar  (         'b1' , '1/b1' ) 
    s2    = frame_statVar ( frame , 'b1' , '1/b1' ) 

    row = 'StatVar    mean (b1,1/b1)' , s1.mean().toString ( '%+.2f +/- %-.2f' ) , s2.mean().toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )
    
    s1    = tree.statVar  (         'b1' , 'b1>0' ) 
    s2    = frame_statVar ( frame , 'b1' , 'b1>0' ) 

    row = 'StatVar    mean (b1,b1>0)' ,  s1.mean().toString ( '%+.2f +/- %-.2f' ) , s2.mean().toString ( '%+.2f +/- %-.2f' ) 
    rows.append ( row )

    if Frames_OK :
        
        s1    = tree.arithmetic_mean      (         'b1' ) 
        s2    = frame_arithmetic_mean ( frame , 'b1' ) 
        
        row = 'Arithmetic mean (b1)' , '%+.2f' % s1.mean() , '%+.2f' % s2.value() 
        rows.append ( row )
        
        s1    = tree.arithmetic_mean  (         'b1' , '1/b1') 
        s2    = frame_arithmetic_mean ( frame , 'b1' , '1/b1') 
        
        row = 'Arithmetic mean (b1,1/b1)' , '%+.2f' % s1.mean()  , '%+.2f' % s2.value() 
        rows.append ( row )

        s1    = tree.arithmetic_mean  (         'b1' , 'b1>0') 
        s2    = frame_arithmetic_mean ( frame , 'b1' , 'b1>0') 
        row = 'Arithmetic mean (b1,b1>0)' , '%+.2f' % s1.mean() , '%+.2f' % s2.value() 
        rows.append ( row )
        
        s1    = tree.geometric_mean  (         'b1' ) 
        s2    = frame_geometric_mean ( frame , 'b1' ) 
    
        row = 'Geometric  mean (b1)' , '%+.2f' % s1.mean()  , '%+.2f' % s2.value() 
        rows.append ( row )
    
        s1    = tree.geometric_mean  (         'b1' , '1/b1') 
        s2    = frame_geometric_mean ( frame , 'b1' , '1/b1') 
        
        row = 'Geometric  mean (b1,1/b1)' , '%+.2f' % s1.mean()   , '%+.2f' % s2.value() 
        rows.append ( row )
        
        s1    = tree.geometric_mean  (         'b1' , 'b1>0') 
        s2    = frame_geometric_mean ( frame , 'b1' , 'b1>0') 
        
        row = 'Geometric  mean (b1,b1>0)' , '%+.2f' % s1.mean() , '%+.2f' % s2.value() 
        rows.append ( row )

        s1    = tree.harmonic_mean  (         'b1' ) 
        s2    = frame_harmonic_mean ( frame , 'b1' ) 
        
        row = 'Harmonic   mean (b1)' , '%+.2f' % s1.mean() , '%+.2f' %  s2.value() 
        rows.append ( row )
        
        s1    = tree.harmonic_mean  (         'b1' , '1/b1') 
        s2    = frame_harmonic_mean ( frame , 'b1' , '1/b1') 
        
        row = 'Harmonic   mean (b1,1/b1)' , '%+.2f' % s1.mean() , '%+.2f' % s2.value() 
        rows.append ( row )
        
        s1    = tree.harmonic_mean  (         'b1' , 'b1>0') 
        s2    = frame_harmonic_mean ( frame , 'b1' , 'b1>0') 
        
        row = 'Harmonic   mean (b1,b1>0)' , '%+.2f' % s1.mean() , '%+.2f' % s2.value() 
        rows.append ( row )
        
        s1    = tree.power_mean  (         2 , 'b1' ) 
        s2    = frame_power_mean ( frame , 2 , 'b1' ) 
        
        row = 'Power      mean (2,b1)' , '%+.2f' % s1.mean() , '%+.2f' % s2.value()
        rows.append ( row )
        
        s1    = tree.power_mean  (         2 , 'b1' , '1/b1') 
        s2    = frame_power_mean ( frame , 2 , 'b1' , '1/b1') 
    
        row = 'Power      mean (2,b1,1/b1)' , '%+.2f' % s1.mean()  , '%+.2f' % s2.value() 
        rows.append ( row )
        
        s1    = tree.power_mean  (         2 , 'b1' , 'b1>0') 
        s2    = frame_power_mean ( frame , 2 , 'b1' , 'b1>0') 
        
        row = 'Power      mean (2,b1,b1>0)' , '%+.2f' % s1.mean() , '%+.2f' % s2.value() 
        rows.append ( row )
        
        s1    = tree.lehmer_mean  (         3 , 'b1' ) 
        s2    = frame_lehmer_mean ( frame , 3 , 'b1' ) 
        
        row = 'Lehmer     mean (3,b1)' , '%+.2f' % s1.mean()  , '%+.2f' % s2.value() 
        rows.append ( row )
        
        s1    = tree.lehmer_mean  (         3 , 'b1' , '1/b1') 
        s2    = frame_lehmer_mean ( frame , 3 , 'b1' , '1/b1') 
        
        row = 'Lehmer     mean (3,b1,1/b1)' , '%+.2f' % s1.mean() , '%+.2f' % s2.value() 
        rows.append ( row )
        
        s1    = tree.lehmer_mean  (         3 , 'b1' , 'b1>0') 
        s2    = frame_lehmer_mean ( frame , 3 , 'b1' , 'b1>0') 
        
        row = 'Lehmer     mean (3,b1,b1>0)' , '%+.2f' % s1.mean() , '%+.2f' % s2.value() 
        rows.append ( row )
        
    title = 'Frame/Tree statistics' 
    logger.info ( '%s\n%s' % ( title , T.table ( rows , title = title , prefix = '# ' , alignment = 'lcc' ) ) ) 
        
# =============================================================================
def test_frame3 () :

    logger = getLogger ( 'test_frame3' ) 
    if not Frames_OK : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 
   
    logger.info ( 80*'*' ) 
    logger.info ( 'MT enabled?  %s' % ROOT.ROOT.IsImplicitMTEnabled() ) 

    logger.info ( 80*'*' ) 
    logger.info ( 'MT enabled?  %s' % ROOT.ROOT.IsImplicitMTEnabled() ) 

    tree   = Tree      ( name = tname , file = fname ).chain

    for i in range ( 5 ) :
        
        frame   = DataFrame ( tname       , fname        )
        fr, pb  = frame_progress ( frame  , len ( tree ) )
        message = 'Value: %d %s' % ( i , pb.GetValue() )
        logger.info ( message ) 
        
# =============================================================================
def test_frame4 () :

    logger = getLogger ( 'test_frame4' ) 
    if not Frames_OK : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 
                   
    logger.info ( 80*'*' ) 
    logger.info ( 'MT enabled?  %s' % ROOT.ROOT.IsImplicitMTEnabled() ) 

    frame  = DataFrame     ( tname        , fname        )
    tree   = Tree      ( name = tname , file = fname ).chain

    pb     = frame_progress ( frame  , len ( tree ) )
    for i,c in enumerate ( range ( 0 , 10000 , 1000 ) , start = 1 )  : 
        frame  = frame.Filter ( 'b1>%s' % c , '#cut%d' % i )
    cnt    = frame.Count ()
    report = frame.Report()

    from ostap.frames.frames import report_print, report_as_table
    title = 'Filter summary'
    logger.info ( '%s\n%s' % ( title , report_print ( report , title = title , prefix = '# ') ) )

# =============================================================================
def test_frame5 () :

    logger = getLogger ( 'test_frame5' ) 
    if not Frames_OK : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 
    
    logger.info ( 80*'*' ) 
    logger.info ( 'MT enabled?  %s' % ROOT.ROOT.IsImplicitMTEnabled() ) 

    frame = DataFrame      ( tname        , fname        )
    tree  = Tree           ( name = tname , file = fname ).chain
    
    from ostap.frames.frames import frame_table

    ft    = frame_table ( frame )

    logger.info ( 'The frame:\n%s' % ft ) 

# =============================================================================
def test_frame6 () :

    logger = getLogger ( 'test_frame6' ) 
    if not Frames_OK : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return 
                      
    logger.info ( 80*'*' ) 
    logger.info ( 'MT enabled?  %s' % ROOT.ROOT.IsImplicitMTEnabled() ) 

    frame = DataFrame      ( tname        , fname        )
    tree  = Tree           ( name = tname , file = fname ).chain
    
    pb    = frame_progress ( frame  , len ( tree ) )
    
    from ostap.frames.frames import frame_project
    
    h1 = ROOT.TH1D ( hID() , '' , 20 , -3, 3 )
    
    rh = frame_project  ( frame , h1 , 'b4' , progress = True , report = True )

    with use_canvas ( 'test_frame6: b4' , wait = 2 ) :
        h1.draw () 

# =============================================================================
def test_frame7 () :

    logger = getLogger ( 'test_frame7' ) 
    if not Frames_OK : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return
     
    logger.info ( 80*'*' ) 
    logger.info ( 'MT enabled?  %s' % ROOT.ROOT.IsImplicitMTEnabled() ) 

    frame = DataFrame      ( tname        , fname        )
    tree  = Tree           ( name = tname , file = fname ).chain
    
    pb    = frame_progress ( frame  , len ( tree ) )
    
    bs = Ostap.Math.Bernstein    ( 10 , -3 , 3 )
    ls = Ostap.Math.LegendreSum  ( 12 , -3 , 3 )
    cs = Ostap.Math.ChebyshevSum ( 15 , -3 , 3 )

    from ostap.frames.frames import frame_param

    with timing ( 'param: Bernstein [12]' , logger = logger ) : 
        rb = frame_param ( frame , bs , 'b4' , '(b4>-1)*1.01' )
    with timing ( 'param: Legendre  [15]' , logger = logger ) :     
        rl = frame_param ( frame , ls , 'b4' , '(b4>-1)*1.01' )
    with timing ( 'param: Chebyshev [18]' , logger = logger ) :             
        rc = frame_param ( frame , cs , 'b4' , '(b4>-1)*1.01' )

    with timing ( 'run  : Bernstein [12]' , logger = logger ) : polb = rb.GetValue()
    with timing ( 'run  : Legendre  [15]' , logger = logger ) : poll = rl.GetValue()
    with timing ( 'run  : Chebyshev [18]' , logger = logger ) : polc = rc.GetValue()

    with use_canvas ( 'test_frame7: b4 as polynomials [12/15/18] ' , wait = 2 ) :
        polb.draw (          linecolor = 2 , linewidth = 2 ) 
        poll.draw ( 'same' , linecolor = 4 , linewidth = 2 ) 
        polc.draw ( 'same' , linecolor = 8 , linewidth = 2 ) 

# =============================================================================
def test_frame8 () :

    logger = getLogger ( 'test_frame8' ) 
    if not Frames_OK : 
        logger.warning ( "Test is disabled for this version of ROOT %s" % str ( root_info ) )
        return
      
    logger.info ( 80*'*' ) 
    logger.info ( 'MT enabled?  %s' % ROOT.ROOT.IsImplicitMTEnabled() ) 
       
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

    if Frames_OK : 
        
        test_frame0 () 
        test_frame1 ()
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
