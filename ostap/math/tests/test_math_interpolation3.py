#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# =============================================================================
## @file ostap/math/tests/test_math_interpolation3.py
#  Test module for the file ostap/math/interpolation.py
# ============================================================================= 
""" Test module for ostap/math/interpolation.py
"""
# ============================================================================= 
from   ostap.utils.utils        import vrange 
import ostap.core.pyrouts 
import ostap.fitting.models     as M 
from   ostap.math.interpolation import ( Berrut1st           ,
                                         Berrut2nd           ,
                                         Barycentric         ,
                                         FloaterHormann      ,
                                         uniform_abscissas   ,
                                         chebyshev_abscissas ,
                                         lobatto_abscissas   ) 
                                        
from   ostap.core.core          import Ostap
from   ostap.math.models        import f1_draw 
from   ostap.utils.utils        import wait
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
import ROOT, random, math
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_interpolation2' ) 
else                       : logger = getLogger ( __name__                         )
# =============================================================================

# =============================================================================
## def test_fitresults ()  :
if 1 < 2 : 

    sigma = lambda x : 5 * math.exp  ( x / 100.0 )

    logger = getLogger ( 'test_fitresults' )
    logger.info ( 'Test interpolation for fit results' )
    
    g1 = M.Gauss_pdf ( 'G' ,
                       xvar  = ( 0 ,    100   ) ,
                       mean  = ( 50 , 1 , 99  ) ,
                       sigma = ( 10 ,  1 , 20  ) )
    g1.mean.fix ( 50 )

    results = {}
    
    low, high, N = 0 , 100 ,  8
    
    gr    = ROOT.TGraphErrors ( N + 1 ) 
    ## for i , m in enumerate ( vrange ( low , high , N ) ) :
    ## for i , m in enumerate ( chebyshev_abscissas ( low , high , N ) ) :
    for i , m in enumerate ( lobatto_abscissas ( low , high , N ) ) :
        
        g1.sigma.fix ( sigma ( m ) )
        g1.mean.fix ( 50 )
        ds = g1.generate ( 5000 )
        g1.mean  .release () 
        g1.sigma.release () 
        
        r , f = g1.fitTo ( ds , draw = True , nbins = 20 , silent = True ) 

        s = r.sigma_G*1.0
        results [ m ] = s

        gr [i]       = m , s
        del ds
        
    interpolants = [
        Berrut1st   ( results  ) ,
        Berrut2nd   ( results  ) ,
        Barycentric ( results  )
        ] + [
        FloaterHormann ( results , i ) for i in range ( 10 ) ] 

    graphs = []
    with wait ( 5 ) , use_canvas ( 'test_fitresults' ) : 
        f1_draw ( sigma , linecolor=1 , linewidth = 4 , minimum = 0 , xmin = 0 , xmax = 100 )        
        gr.red()
        gr.draw ('pe1')
    
        K = 255 
        
        for i,f in enumerate ( interpolants ) : 
            g = ROOT.TGraphErrors ( K + 1 ) 
            for j,x in enumerate ( vrange ( 0 , 100 , K ) ) :
                g[j] = x , f ( x )
            graphs.append ( g )
                
        for i, g in enumerate ( graphs ) :
            g.draw ('lpe1', linecolor = i+3 , markercolor = i+3 ) 

        gr.draw ('pe1')
    
    
        
    
# =============================================================================
## if '__main__' == __name__ :
    
##    test_fitresults ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
