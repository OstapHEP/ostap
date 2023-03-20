#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file ostap/histos/tests/test_histos_histos.py
# Test module for ostap/histos/histos.py
# - It tests the basic operations with histograms  
# ============================================================================= 
"""Test module for ostap/histos/histos.py
- It tests the basic operations with histograms  
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import 
# ============================================================================= 
from   ostap.math.ve        import VE 
from   ostap.core.core      import hID 
from   ostap.histos.histos  import h1_axis, h2_axes 
from   builtins             import range
import ROOT, random, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'ostap.test_histos_histos' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( 'Test for basic operations with histograms')
# =============================================================================
# =============================================================================
## Test for very basic operations with 1D-histograms
def test_basic_1D() :

    logger.info ( 'Test for very basic operations with 1D-histograms')

    h1 = ROOT.TH1D ( hID() , '' , 10 , 0 , 1 )

    ## clone it! 
    h2 = h1.clone() 

    ## set content for certain bins:
    h1[1] = 1.55
    h1[2] = VE(1,0.2**2)

    ## get content from certain bins:
    a = h1[3]

    h1 += VE(3,0.5**2)
    
    ## add a function
    h1 += lambda x : VE(25*x*x,25*x*x)

    ## Divide by function
    h1 /= lambda x : float(1+x*x) 

    ## multiply by constant
    h1 *= 15.0

    ## loop over bins using index
    for i in h1 :
        logger.info( "bin# %2d, content %s" % (  i , h1[i] ) )

    ## loop over bins i nreverse order 
    for i in reversed( h1 ) :
        logger.info( "bin# %2d, content %s" % (  i , h1[i] ) )

    ## iterate over all 
    for item in h1.items() : 
        logger.info ( "items:  %s" % str ( item ) )

       
    ## interpolate
    for x in range(10) :
        x = random.uniform(0,1)
        logger.info ( 'h1(%.3f) = %s[d] %s[0] %s[1] %s[2] %s[3] ' % ( x     ,
                                                                      h1(x) ,
                                                                      h1(x, interpolate=0) ,
                                                                      h1(x, interpolate=1) ,
                                                                      h1(x, interpolate=2) ,
                                                                      h1(x, interpolate=3) ) )

    ## running sum
    h1_d = h1.sumv()                     ## default 
    h1_i = h1.sumv( increasing = True  ) ## ditto
    h1_r = h1.sumv( increasing = False ) ## ditto

    ## histogram efficiney of cuts
    eff_i = h1.effic( increasing = True   )
    eff_r = h1.effic( increasing = False  )

    ## efficiency of certaine value
    e1_i  = h1.efficiency( 0.3 , increasing = True  )
    e2_i  = h1.efficiency( 0.3 , increasing = False )


    ## smear the histogram
    h0    = ROOT.TH1F( hID() , '', 400 , -1 , 1 )
    for i in range( 1000 ) : h0.Fill( random.gauss(0,0.10 ) )
    hs    = h0.smear ( sigma = 0.10 )

    logger.info ( 'Original RMS %20s , smeared(0.10) %-20s' % ( h0.rms() , hs.rms() ) )  
    ## specific transformation:

    ## "precision"
    hp  = h1.precision()
    
    ## "B/S"
    hb2s = h1.b2s() 

    ## rescale histo

    hs1 = h1.scale(1)
    hs2 = h1.rescale_bins(1) 

    ## sample
    hr  = h1.sample()


    ## "figure of merit" to maximize precision:

    fm2 = h1.fom_2()
    fm2 = h1.FoM_2()

    ## rebin template:
    h_tmpl = h1_axis ( [ 0 , 0.1 , 0.5 , 0.55, 0.80 , 0.95 , 1.0 ] )

    ## rebin as "numbers"
    h1_n   = h1.rebinNumbers  ( h_tmpl )
    
    ## rebin as "function"
    h1_f   = h1.rebinFunction ( h_tmpl )

    ## slice for the histogram
    hs   = h1[4:30]

    ## accumulate
    a1 = h1.accumulate( low = 1    , high = 10 ) 
    a2 = h1.accumulate( xmin = 0.1 , xmax = 0.30 ) 

    ## shift
    hs = h1.shift ( -0.3 )

    ## other operations
    
    h2 = h1**2

    h2 = h1*h1

    h2 = h1*2

    h2 = h1/2

    h2 = 2*h1
    
    h2 = abs(h1) 
        
    from ostap.math.math_ve import sin, exp 
    h2 = sin ( h1 )
    h2 = exp ( h1 )

    
    ## shift for bins
    hs = h1 >> 4 
    hs = h1 << 5 

    ## integration (taking into account bin-width)
    i1 = h1.integrate()
    i2 = h1.integrate( lowx = 1   , highx = 15  ) 
    i2 = h1.integrate( xmin = 0.1 , xmax  = 0.6 )  

    ## statistics
    for i in range ( 0 , 5 ) :
        logger.info ( "        Moment (%d): %-20s" % ( i , h1.moment        ( i ) ) )
    ## statistics
    for i in range ( 0 , 5 ) :
        logger.info ( "Central moment (%d): %-20s" % ( i , h1.centralMoment ( i ) ) )

    logger.info (     "              Mean : %-20s" % h1.mean     () ) 
    logger.info (     "               RMS : %-20s" % h1.rms      () )  
    logger.info (     "          Skewness : %-20s" % h1.skewness () )  
    logger.info (     "          Kurtosis : %-20s" % h1.kurtosis () ) 

    logger.info (     "  Stat: %s" % h1. stat() )
    logger.info (     " WStat: %s" % h1.wstat() )
    logger.info (     " XStat: %s" % h1.xstat() )
        
    logger.info ( '  minmax  %20s' % str( h1. minmax() ) )
    logger.info ( 'x-minmax  %20s' % str( h1.xminmax() ) )
    logger.info ( 'y-minmax  %20s' % str( h1.yminmax() ) )


    hh = ROOT.TH1D ( hID() , 'Gaussian' , 500 , -5 , 5 )
    for i in range ( 1000000 ) : hh.Fill ( random.gauss ( 0 , 1 ) ) 
    
    title = 'Histogram moments'
    m     = hh.the_moment ( 24 ) 
    logger.info ( '%s:\n%s' % ( title , m.table  ( title = title , prefix = '# ' ) ) ) 
    logger.info ( 'neff       %-20s' % hh.nEff     () )
    logger.info ( 'mean       %-20s' % hh.mean     () )
    logger.info ( 'rms        %-20s' % hh.rms      () )
    logger.info ( 'skewness   %-20s' % hh.skewness () )
    logger.info ( 'kurtosis   %-20s' % hh.kurtosis () )
                     
# =============================================================================
## Test for very basic operations with 2D-histograms
def test_basic_2D   () :
    
    logger.info ( 'Test for very basic operations with 2D-histograms')

    h2 = h2_axes ( [1,2,3,4,5,6,7] , [1,2,3,4,5] )

    h2 += lambda x,y: VE(1,1)+VE(x*x+2*y*y,x*x+2*y*y)

    ## access the content
    logger.info ( "h[%1d, %1d]=%s"  % ( 2    , 3 , h2[2,3]  ) )
    logger.info ( "h[%1d][%1d]=%s"  % ( 2    , 3 , h2[2][3] ) )
    logger.info ( "h[%1d](%.2f)=%s" % ( 2    , 3.01 , h2[2](3.01   ) ) )
    logger.info ( "h(%.2f,%.2f)=%s" % ( 2.01 , 3.01 , h2(2.01,3.01 ) ) )    

    ## get the 1st X-slice 
    hx = h2.sliceX ( 1 )

    ## get the slice in 1,4,5 X-bins 
    hx = h2.sliceX ( [1,4,5] )
    
    ## get the 2nd Y-slice 
    hy = h2.sliceY ( 2 )
    
    ## get the slice in 1,3,5 Y-bins 
    hy = h2.sliceY ( [1,3,4] )

    ## projection in X 
    hx = h2.projX() 

    ## projection in Y
    hy = h2.projY() 

    logger.info ( '  minmax  %20s' % str( h2. minmax() ) )
    logger.info ( 'x-minmax  %20s' % str( h2.xminmax() ) )
    logger.info ( 'y-minmax  %20s' % str( h2.yminmax() ) )
    logger.info ( 'z-minmax  %20s' % str( h2.zminmax() ) )
                  

# =============================================================================
## Test for "efficiencies
def test_efficiency() :

    logger.info ( 'Test for "efficiencies" ')

    hA = ROOT.TH1F ( hID() , 'accepted' , 10 , 0 , 10 )
    hR = ROOT.TH1F ( hID() , 'rejected' , 10 , 0 , 10 ) 
    hT = ROOT.TH1F ( hID() , 'total'    , 10 , 0 , 10 )

    random.seed (100) 
    for i in range(10000) :

        vx       = random.uniform ( *hA.xminmax() )

        accepted = random.uniform (0,1) < 0.01
        
        if accepted : hA.Fill ( vx )
        else        : hR.Fill ( vx )
        
        hT.Fill ( vx )

    #
    ## applicbale for all cases: 
    #

    eff_0 = hA  / hT                  ## incorrect estimate of uncertainties efficiency

    
    ## (almost) correct estimate of uncertainties
    eff_1 = 1 / ( 1 + hR/hA )         ## correct estimate of uncertainties
    
    ## use Zech's prescription
    eff_2 = hA % hT                   ## use Zech's prescription, applicable for all cases
    eff_3 = hA.zechEff ( hT )         ## ditto 

    ## only for histograms with true natural entries:
    
    eff_4 = hA // hT                  ## use binomial uncertainties
    eff_5 = hA.binomEff        ( hT ) ## ditto 
    eff_6 = hA.agrestiCoullEff ( hT ) ## use Argesti-Coull's recipe 
    eff_7 = hA.wilsonEff       ( hT ) ## use Wilson's recipe 

    ## asymmetric binomial intervals
    #  only for true natural entries
    #  result in form of graph
    
    geff_8  = hA.eff_wald                    ( hR )
    geff_9  = hA.eff_wilson_score            ( hR )
    geff_10 = hA.eff_wilson_score_continuity ( hR )
    geff_11 = hA.eff_arcsin                  ( hR )
    geff_12 = hA.eff_agresti_coull           ( hR )
    geff_13 = hA.eff_jeffreys                ( hR )
    geff_14 = hA.eff_clopper_pearson         ( hR )

    for e,n in [ ( eff_0   , 'Naive'                                      ) ,
                 ( eff_1   , 'General'                                    ) ,
                 ( eff_2   , 'Zech'                                       ) ,
                 ( eff_3   , 'Zech/Operator'                              ) ,
                 ( eff_4   , 'Binomial: Simple'                           ) ,
                 ( eff_5   , 'Binomial: Simple/Operator'                  ) ,
                 ( eff_6   , 'Binomial: Agresti-Coull'                    ) ,
                 ( eff_7   , 'Binomial: Wilson'                           ) ,
                 ( geff_8  , 'Binomial: Wald interval'                    ) ,
                 ( geff_9  , 'Binomial: Wilson score interval'            ) ,
                 ( geff_10 , 'Binomial: Wilson score/continuity interval' ) ,
                 ( geff_11 , 'Binomial: Arcsin interval'                  ) ,
                 ( geff_12 , 'Binomial: Agresti-Coull interval'           ) ,
                 ( geff_13 , 'Binomial: Jeffreys interval'                ) , 
                 ( geff_14 , 'Binomial: Clopper-Pearson interval'         ) ] :
        
        
        if isinstance ( e , ROOT.TH1 ) :
            logger.info ( "%43s: %s" % ( n , [ e[i]*100 for i in e[:4] ] )  )
        else :
            vals = [ (e[i][3]-abs(e[i][4]),e[i][3]+e[i][5]) for i in e[:3] ]
            vals = [  "(%7.4f,%7.4f)" % ( e[0]*100 , e[1]*100 ) for e in vals  ]
            logger.info ( "%43s: %s" % ( n , vals ) ) 
            
                    
            

# =============================================================================
if '__main__' == __name__ :

    
    test_basic_1D   ()
    ## test_basic_2D   ()
    ## test_efficiency () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
