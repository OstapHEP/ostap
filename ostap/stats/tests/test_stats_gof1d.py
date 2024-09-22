# -*- coding: utf-8 -*-
# ==============================================================================
# @file ostap/stats/tests/test_stats_gof1d.py 
# Test Goddness-of-fits 1D 
# Copyright (c) Ostap developpers.
# ==============================================================================
""" # Test averages for inconsistend data 
"""
# ==============================================================================
import ostap.fitting.models  as     M 
import ostap.stats.gof1d     as     G1D 
from   ostap.plotting.canvas import use_canvas
import ROOT  
# ==============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'tests_stats_gof1d' )
else                       : logger = getLogger ( __name__            )
# ==============================================================================
xvar  = ROOT.RooRealVar ( 'x', '', 0, 10)
gauss = M.Gauss_pdf     ( 'G' , xvar = xvar , mean = 5 , sigma = 1 )
model = M.Fit1D         ( signal = gauss , background = 'flat'     )


ND = 200
# ==============================================================================
def test_good_fit_1 ( ) :
    
    logger = getLogger ( 'test_good_fit_1' )
    logger.info ( 'Make a test for presumably good fit' )

    gauss.mean  = 5
    gauss.sigma = 0.75
    data = gauss.generate ( ND , sample = True )


    with use_canvas ( 'test_good_fit_1' ,      wait = 1 ) :
        gauss.fitTo ( data , draw = True , nbins = 50 , quite = True ) 

        
    with use_canvas ( 'test_good_fit_1: GoF' , wait = 1 ) :
        gof = G1D.GoF1D     ( gauss , data )
        logger.info ( 'Goodness-of-fit:\n%s' % gof )

        got = G1D.GoF1DToys ( gauss , data )
        logger.info ( 'Goodness-of-fit with toys:\n%s' % got )

# =============================================================================
def test_good_fit_2 ( ) :
    
    logger = getLogger ( 'test_good_fit_2' )
    logger.info ( 'Make a test for presumably good fit' )

    gauss.mean  = 5
    gauss.sigma = 0.75
    data = gauss.generate ( ND , sample = True )

    
    with use_canvas ( 'test_good_fit_2' ,      wait = 1 ) :
        model.fitTo ( data , draw = True , nbins = 50 , refit = 5 , quite = True ) 

        
    with use_canvas ( 'test_good_fit_2: GoF' , wait = 1 ) :
        gof = G1D.GoF1D     ( model , data )
        logger.info ( 'Goodness-of-fit:\n%s' % gof )

        got = G1D.GoF1DToys ( model , data , 500 )
        got.run ( 500 ) 
        got.run ( 500 ) 
        logger.info ( 'Goodness-of-fit with toys:\n%s' % got )


# =============================================================================
def test_bad_fit_1 ( ) :
    
    logger = getLogger ( 'test_bad_fit_1' )
    logger.info ( 'Make a test for presumably bad fit' )

    gauss.mean  = 5
    gauss.sigma = 0.75
    model.S     = 0.96 * ND
    model.B     = 0.04 * ND
    
    data = model.generate ( ND , sample = True )


    with use_canvas ( 'test_bad_fit_1' ,      wait = 1 ) :
        gauss.fitTo ( data , draw = True , nbins = 50 , quite = True ) 

    with use_canvas ( 'test_bad_fit_1: GoF' , wait = 1 ) :
        gof = G1D.GoF1D     ( gauss , data )
        logger.info ( 'Goodness-of-fit:\n%s' % gof )
        gof.draw()        
        got = G1D.GoF1DToys ( gauss , data , 5000 )
        got.run ( 10000 ) 
        logger.info ( 'Goodness-of-fit with toys:\n%s' % got )

    with use_canvas ( 'test_bad_fit_1: GoF/Kolmogorov-Smirnov' , wait = 1 ) :
        dks = got.draw('Kolmogorov-Smirnov')
    with use_canvas ( 'test_bad_fit_1: GoF/Andersen-Darling' , wait = 1 ) :
        dad = got.draw('Andersen-Darling')
    with use_canvas ( 'test_bad_fit_1: GoF/Cramer-von Mises' , wait = 1 ) :
        dcm = got.draw('Cramer-von Mises')
    with use_canvas ( 'test_bad_fit_1: GoF/Kuiper'           , wait = 1 ) :
        dcm = got.draw('Kuiper')
    with use_canvas ( 'test_bad_fit_1: GoF/ZK' , wait = 1 ) :
        dzk = got.draw('ZK')
    with use_canvas ( 'test_bad_fit_1: GoF/ZA' , wait = 1 ) :
        dza = got.draw('ZA')
    with use_canvas ( 'test_bad_fit_1: GoF/ZC' , wait = 3 ) :
        dzc = got.draw('ZC')
        
# ===============================================================================
if '__main__' == __name__ :

    test_good_fit_1 ()
    test_good_fit_2 ()
    test_bad_fit_1  ()

# ===============================================================================
##                                                                        The END 
# ===============================================================================
                   
                   
