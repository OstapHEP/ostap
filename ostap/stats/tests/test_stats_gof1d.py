#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==============================================================================
# @file ostap/stats/tests/test_stats_gof1d.py 
# Test Goddness-of-fits for 1D fits 
# Copyright (c) Ostap developpers.
# ==============================================================================
""" Test Goddness-of-fits for 1D fits 
"""
# ==============================================================================
import ostap.fitting.models  as     M 
import ostap.stats.gof1d     as     G1D 
import ostap.stats.gofnd     as     GnD
import ostap.logger.table    as     T
from   ostap.plotting.canvas import use_canvas
from   ostap.logger.pretty   import pretty_float
from   ostap.utils.utils     import vrange 
import ROOT  
# ==============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'tests_stats_gof1d' )
else                       : logger = getLogger ( __name__            )
# ==============================================================================
xvar  = ROOT.RooRealVar ( 'x', '', 0, 10)
gauss = M.Gauss_pdf     ( 'G' , xvar = xvar , mean = 5 , sigma = 1 )
model = M.Fit1D         ( signal = gauss , background = 'flat'     )

ND    = 200

# =============================================================================
## data_g: pure gaussian
gauss.mean  = 5
gauss.sigma = 0.75
data_g      = gauss.generate ( ND , sample = True )
# =============================================================================
## data_b:  95% gaussian + 4% flat background 
model.S     = 0.95 * ND
model.B     = 0.05 * ND
data_b      = model.generate ( ND , sample = True )


fitconf = { 'draw' : True , 'nbins' : 50 , 'refit' : 5 , 'quiet' : True }
# ==============================================================================
## Run Point-to_pint dissimilatory Goodness-of-Fit test
def run_PPD ( pdf , data, result , logger ) :
    """ Run Point-to_pint dissimilatory Goodness-of-Fit test"""
    
    from ostap.stats.gof_np import np,sp,s2u,cdist
    if not np or not sp or not s2u or not cdist:
        logger.warning ('No numpy/scipy/s4u/cdist: skip the PPD estimate!')
        return         

    #  - Point to Point Sissimilarity test  with Gaussian distance using different "sigma"
    rows  =  [ ( 'PPD/sigma' , 't-value'  , 'x[..]', 'p-value [%]' ) ]
    Ns    = 20  
    logger.info ( 'Run Point-to-Point Dissimilarity GoF-test for %d different values of sigma' % Ns  ) 
    for sigma in vrange ( 0.01 , 2.0 , Ns ) :
        
        ppd = GnD.PPD ( Nperm = 1000 , sigma = sigma )
        pdf.load_params ( result , silent = True ) 
        tvalue         = ppd        ( pdf , data )
        tvalue, pvalue = ppd.pvalue ( pdf , data )
        
        tv , texpo = pretty_float ( tvalue )
        pvalue *= 100
        pvalue  = '%4.1f +/- %.1f' %  ( pvalue.value() , pvalue.error() )
        row = '%.2f' % sigma , tv , '10^%+d' % texpo if texpo else '' , pvalue 
        rows.append ( row ) 
        
    title= 'Goodness-of-fit PPD-test'
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )


# ==============================================================================
def test_good_fit_1 ( ) :
    """ Make a test for presumably good fit: fit Gauss to Gauss
    """
    
    logger = getLogger ( 'test_good_fit_1' )
    logger.info ( 'Make a test for presumably good fit: fit Gauss to Gauss' )

    with use_canvas ( 'test_good_fit_1: G -> G' ,      wait = 1 ) :
        r , f = gauss.fitTo ( data_g , **fitconf ) 
        
    with use_canvas ( 'test_good_fit_1: GoF' , wait = 1 ) :
        
        gauss.load_params ( r , silent = True ) 
        gof = G1D.GoF1D     ( gauss , data_g )
        logger.info ( 'Goodness-of-fit:\n%s' % gof )

        gauss.load_params ( r , silent = True ) 
        got = G1D.GoF1DToys ( gauss , data_g )
        logger.info ( 'Goodness-of-fit with toys:\n%s' % got )

    ## Try to use multidimensional methods
    run_PPD ( gauss , data_g , r , logger )
    
# =============================================================================
def test_good_fit_2 ( ) :
    """ Make a test for presumably good fit: fit Gauss+Bkg to Gauss
    """
    
    logger = getLogger ( 'test_good_fit_2' )
    logger.info ( 'Make a test for presumably good fit: fit Gauss+Bkg to Gauss' )

    with use_canvas ( 'test_good_fit_2: G+B -> G' ,      wait = 1 ) :
        r , f = model.fitTo ( data_g , **fitconf ) 
        
    with use_canvas ( 'test_good_fit_2: GoF' , wait = 1 ) :
        model.load_params ( r , silent = True ) 
        gof = G1D.GoF1D     ( model , data_g )
        logger.info ( 'Goodness-of-fit:\n%s' % gof )

        model.load_params ( r , silent = True ) 
        got = G1D.GoF1DToys ( model , data_g , 500 )
        logger.info ( 'Goodness-of-fit with toys:\n%s' % got )

    ## Try to use multidimensional methods
    run_PPD ( model , data_g , r , logger )
    

# =============================================================================
def test_good_fit_3 ( ) :
    
    logger = getLogger ( 'test_good_fit_3' )
    logger.info ( 'Make a test for presumably good fit: fit Gauss+Bkg to Gauss+Bkg' )

    with use_canvas ( 'test_good_fit_3: G+B -> G+G' ,      wait = 1 ) :
        r , f = model.fitTo ( data_b , **fitconf ) 
        
    with use_canvas ( 'test_good_fit_2: GoF' , wait = 1 ) :
        model.load_params ( r , silent = True ) 
        gof = G1D.GoF1D     ( model , data_b )
        logger.info ( 'Goodness-of-fit:\n%s' % gof )

        model.load_params ( r , silent = True ) 
        got = G1D.GoF1DToys ( model , data_b , 500 )
        logger.info ( 'Goodness-of-fit with toys:\n%s' % got )

    ## Try to use multidimensional methods
    run_PPD ( model , data_b , r , logger )
            
# =============================================================================
def test_bad_fit_1 ( ) :
    """ Make a test for presumably bad fit: fit Gauss to Gauss+Bkg
    """
    logger = getLogger ( 'test_bad_fit_1' )
    logger.info ( 'Make a test for presumably bad fit: fit Gauss to Gauss+Bkg' )

    with use_canvas ( 'test_bad_fit_1: G -> G+B' ,      wait = 1 ) :
        r , f = gauss.fitTo ( data_b , **fitconf ) 

    with use_canvas ( 'test_bad_fit_1: GoF' , wait = 1 ) :
        
        gauss.load_params ( r , silent = True ) 
        gof = G1D.GoF1D     ( gauss , data_b )
        logger.info ( 'Goodness-of-fit:\n%s' % gof )
        gof.draw()
        
        gauss.load_params ( r , silent = True ) 
        got = G1D.GoF1DToys ( gauss , data_b , 5000 )
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

    ## Try to use multidimensional methods
    run_PPD ( gauss , data_b , r , logger )
    
        
# ===============================================================================
if '__main__' == __name__ :

    test_good_fit_1 ()  ## fit Gauss       to Gauss 
    test_good_fit_2 ()  ## fit Gauss+Bkg   to Gauss 
    test_good_fit_3 ()  ## fit Gauss+Bkg   to Gauss+Bkg
    test_bad_fit_1  ()  ## fit Gauss       to Gauss+Bkg
 
# ===============================================================================
##                                                                        The END 
# ===============================================================================
                   
                   
