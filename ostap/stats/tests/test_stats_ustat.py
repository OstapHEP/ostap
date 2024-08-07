#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/stats/tests/test_stats_ustat.py 
#  Test uStatistics for Goodness-Of-Fit tests
# Copyright (c) Ostap developpers.
# ============================================================================= 
""" Test uStatistics for goodness-of-fit tests
"""
# =============================================================================
import ostap.stats.ustat     as     uStat 
import ostap.fitting.models  as     Models
from   ostap.utils.timing    import timing
from   ostap.core.pyrouts    import SE
from   ostap.plotting.canvas import use_canvas 
import ROOT, random, math  
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'tests_stats_ustat' )
else                       : logger = getLogger ( __name__            )
# ==============================================================================


histos = set() 

def test_stats_ustat_G2D () :

    logger = getLogger ( "test_stats_ustat_G2D" )
    
    x   = ROOT.RooRealVar ( 'x' , 'x-variable' , 0 , 10 )
    y   = ROOT.RooRealVar ( 'y' , 'y-variable' , 0 , 10 )
    
    pdf = Models.Gauss2D_pdf ( 'G2D' , x , y ,
                               muX    = ( 5 , 4   , 6   ) ,
                               muY    = ( 5 , 4   , 6   ) ,
                               sigmaX = ( 1 , 0.5 , 1.5 ) , 
                               sigmaY = ( 2 , 1.5 , 2.5 ) , 
                               theta  = ( math.pi/4 , math.pi/8 ,  math.pi/2 ) )  
    

        
    for n in ( 100 , 200 , 500 , 1000 ) : ## , 5000 , 10000  ) :


        title = "2D: N=%4d test" % n 
        with timing ( title , logger = logger ) , use_canvas ( title , wait = 5 ) :
            
            pdf.muX    = 5
            pdf.muY    = 5
            pdf.sigmaX = 1 
            pdf.sigmaY = 2  
            pdf.theta  = math.pi/4 
            
            
            data       = pdf.generate ( n )
            
            pdf.fitTo ( data , silent = True )
            
            t , h , r  = uStat.uPlot ( pdf  , data ) 
            
            histos.add ( h )  
            h.draw()
            
            data.clear()
            del data

def test_stats_ustat_G3D () :

    logger = getLogger ( "test_stats_ustat_G3D" )
    
    x   = ROOT.RooRealVar ( 'x1' , 'x-variable' , 0 , 10 )
    y   = ROOT.RooRealVar ( 'y2' , 'y-variable' , 0 , 10 )
    z   = ROOT.RooRealVar ( 'z1' , 'z-variable' , 0 , 10 )
    
    pdf = Models.Gauss3D_pdf ( 'G3D' , x , y , z , 
                               muX    = ( 5 , 4   , 6   ) ,
                               muY    = ( 5 , 4   , 6   ) ,
                               muZ    = ( 5 , 4   , 6   ) ,
                               sigmaX = ( 1 , 0.5 , 1.5 ) , 
                               sigmaY = ( 2 , 1.5 , 2.5 ) , 
                               sigmaZ = ( 3 , 2.5 , 3.5 ) , 
                               phi    = ( math.pi/4 , math.pi/8 ,  math.pi/2 ) , 
                               theta  = ( math.pi/4 , math.pi/8 ,  math.pi/2 ) , 
                               psi    = ( math.pi/4 , math.pi/8 ,  math.pi/2 ) )
    
    for n in ( 100 , 200 , 500 , 1000 ) : ## , 5000 , 10000  ) :
        
        
        title = "3D: N=%4d test" % n 
        with timing ( title , logger = logger ) , use_canvas ( title , wait = 5 ) :
            
            pdf.muX    = 5
            pdf.muY    = 5
            pdf.muZ    = 5
            pdf.sigmaX = 1 
            pdf.sigmaY = 2  
            pdf.sigmaZ = 3   
            pdf.phi    = math.pi/4 
            pdf.theta  = math.pi/4 
            pdf.psi    = math.pi/4 

            data       = pdf.generate ( n )
            
            pdf.fitTo ( data , silent = True )
            
            t , h , r  = uStat.uPlot ( pdf  , data ) 
            
            histos.add ( h )  
            h.draw()
            
            data.clear()
            del data
            
            


# =============================================================================
if '__main__' == __name__ :

    test_stats_ustat_G2D ()
    test_stats_ustat_G3D ()
    
# =============================================================================
##                                                                      The END 
# ============================================================================= 

