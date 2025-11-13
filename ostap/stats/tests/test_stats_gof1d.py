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
from   ostap.stats.ustat      import USTAT 
from   ostap.plotting.canvas  import use_canvas
from   ostap.logger.pretty    import pretty_float
from   ostap.utils.ranges     import vrange
from   ostap.utils.root_utils import batch_env
from   ostap.utils.timing     import timing 
from   ostap.math.math_ve     import significance
from   ostap.core.core        import VE 
import ostap.fitting.models   as     M 
import ostap.stats.gof1d      as     G1D 
import ostap.stats.gofnd      as     GnD
import ostap.logger.table     as     T
import ROOT  
# ==============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'tests_stats_gof1d' )
else                       : logger = getLogger ( __name__            )
# ==============================================================================
batch_env ( logger ) 
# =============================================================================

xvar  = ROOT.RooRealVar ( 'x', '', 0, 10)
gauss = M.Gauss_pdf     ( 'G' , xvar = xvar , mean = 5 , sigma = 1 )
model = M.Fit1D         ( signal = gauss , background = 'flat' , fix_norm = True )

ND    = 100

# =============================================================================
## data_g: pure gaussian
gauss.mean  = 5
gauss.sigma = 0.75
data_g      = gauss.generate ( ND , sample = True )
# =============================================================================
## data_b:  95% gaussian + 4% flat background 
model.S     = 0.90 * ND
model.B     = 0.10 * ND
data_b      = model.generate ( ND , sample = True )

fitconf = { 'draw' : True , 'nbins' : 25 , 'refit' : 5 , 'quiet' : True }

keep = set() 
# ==============================================================================
## Run Point-to_pint dissimilatity Goodness-of-Fit test
def run_PPD ( pdf , data, result , logger ) :
    """ Run Point-to_point dissimilarity Goodness-of-Fit test"""

    from ostap.math.base    import numpy, scipy 
    if not numpy or not scipy :
        logger.warning ('No numpy/scipy: skip the PPD estimate!')
        return
    from ostap.stats.gof_np import s2u,cdist
    if not s2u or not cdist:
        logger.warning ('No s4u/cdist: skip the PPD estimate!')
        return
        
    # =========================================================================
    #  - Point to Point Dissimilarity test  with Gaussian distance using different "sigma"
    rows  = [ ( 'PPD/sigma' , 't-value'  , 'x[..]', 'p-value [%]' , '#sigma') ]
    Ns    = 3 
    logger.info ( 'Run Point-to-Point Dissimilarity GoF-test for %d different values of sigma' % Ns  ) 
    for sigma in vrange ( 0.1 , 1.0 , Ns ) :

        with timing ( 'PPD-test' , logger = logger ) :             
            ppd = GnD.PPD ( nToys = 200 , sigma = sigma )
            pdf.load_params ( result , silent = True ) 
            tvalue         = ppd          ( pdf , data )
            tvalue, pvalue = ppd.pvalue   ( pdf , data )
            
            clip = 0.1 * pvalue.error()
            pv   = pvalue 
            if   1 <= pvalue.value() : pv = VE ( 1 - clip , pv.cov2() )
            elif 0 >= pvalue.value() : pv = VE (     clip , pv.cov2() )
            nsigma = significance ( pv ) ## convert  it to significace
        
        tv , texpo = pretty_float ( tvalue )
        pvalue *= 100
        pvalue  = '%4.1f +/- %.1f' %  ( pvalue.value() , pvalue.error () )
        nsigma  = '%.1f +/- %.1f'  %  ( nsigma.value() , nsigma.error () )
        row = '%.2f' % sigma , tv , '10^%+d' % texpo if texpo else '' , pvalue , nsigma 
        rows.append ( row ) 
        
    title= 'Goodness-of-fit PPD-test'
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )

import ostap.stats.ustat as U 

# ==============================================================================
## Run Distance-to-Nearest-Neighbour Goodness-of-Fit test
def run_DNN  ( pdf , data, result , logger ) :
    """ Run Distance-to-Nearest-Neighbour Goodness-of-Fit test
    """
    
    from ostap.math.base    import numpy, scipy 
    if not numpy or not scipy :
        logger.warning ('No numpy/scipy: skip the PPD estimate!')
        return
    from ostap.stats.gof_np import s2u,cdist
    if not s2u or not cdist:
        logger.warning ('No s4u/cdist: skip the PPD estimate!')
        return
    
    rows  =  [ ( 't-value'  , 'x[..]', 'p-value [%]' , '#sigma' ) ]

    with timing ( 'DNN-test' , logger = logger ) : 
        dnn = GnD.DNN ( nToys = 200 , histo = 50 )        
        pdf.load_params ( result , silent = True )        
        tvalue         = dnn          ( pdf , data )
        tvalue, pvalue = dnn.pvalue   ( pdf , data )

        clip = 0.1 * pvalue.error()
        pv   = pvalue 
        if   1 <= pvalue.value() : pv = VE ( 1 - clip , pv.cov2() )
        elif 0 >= pvalue.value() : pv = VE (     clip , pv.cov2() )
        nsigma = significance ( pv ) ## convert  it to significace

    tv , texpo = pretty_float ( tvalue )
    pvalue *= 100
    pvalue  = '%4.1f +/- %.1f' %  ( pvalue.value() , pvalue.error() )
    nsigma  = '%.1f +/- %.1f'  %  ( nsigma.value() , nsigma.error () )
    row     = tv , '10^%+d' % texpo if texpo else '' , pvalue , nsigma 
    rows.append ( row )
    
    title= 'Goodness-of-fit DNN-test'
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )

    ## u-distribution
    return dnn.histo 

# ==============================================================================
## Run USTAT Goodness-of-Fit test
def run_USTAT  ( pdf , data, result , logger ) :
    """ Run USTAT Goodness-of-Fit test
    """

    rows  =  [ ( 't-value'  , 'x[..]', 'p-value [%]' , '#sigma' ) ]
    
    ustat = USTAT ( nToys = 200 , histo = 100 )
    
    pdf.load_params ( result , silent = True )
    
    with timing ( 'uStat-test' , logger = logger ) : 
        tvalue         = ustat        ( pdf , data )    
        tvalue, pvalue = ustat.pvalue ( pdf , data )
        
        clip = 0.1 * pvalue.error()
        pv   = pvalue 
        if   1 <= pvalue.value() : pv = VE ( 1 - clip , pv.cov2() )
        elif 0 >= pvalue.value() : pv = VE (     clip , pv.cov2() )
        nsigma = significance ( pv ) ## convert  it to significace

    tv , texpo = pretty_float ( tvalue )
    pvalue *= 100
    pvalue  = '%4.1f +/- %.1f' %  ( pvalue.value() , pvalue.error() )
    nsigma  = '%.1f +/- %.1f'  %  ( nsigma.value() , nsigma.error () )
    row     = tv , '10^%+d' % texpo if texpo else '' , pvalue, nsigma 
    rows.append ( row )
    
    title= 'Goodness-of-Fit USTAT-test'
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )

    return ustat.histo

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
        with timing ( 'GoF1D-toys' , logger = logger ) : 
            toys = G1D.GoF1DToys ( gof )
            toys = toys.run ( nToys = 500 , parallel = True )
        logger.info ( 'Goodness-of-fit with %d toys:\n%s' % ( toys.nToys , toys ) ) 

    del gof
    del toys 

    ## Try to use multidimensional methods
    run_PPD ( gauss , data_g , r , logger )
    
    udist1 = run_DNN ( gauss , data_g , r , logger )
    if udist1 :
        keep.add ( udist1 ) 
        with use_canvas ( 'test_good_fit_1: DNN' , wait = 1 ) :
            udist1.draw()

    udist2 = run_USTAT ( gauss , data_g , r , logger )
    if udist2 :
        keep.add ( udist2 ) 
        with use_canvas ( 'test_good_fit_1: USTAT' , wait = 1 ) :
            udist2.draw()
    
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
        with timing ( 'GoF1D-toys' , logger = logger ) : 
            toys = G1D.GoF1DToys ( gof )
            toys = toys.run ( nToys = 200 ) 
        logger.info ( 'Goodness-of-fit with %d toys:\n%s' % ( toys.nToys , toys ) ) 

    del gof
    del toys
        
    ## Try to use multidimensional methods
    run_PPD ( model , data_g , r , logger )
    udist1 = run_DNN ( model , data_g , r , logger )
    if udist1 :
        keep.add ( udist1 ) 
        with use_canvas ( 'test_good_fit_2: DNN' , wait = 1 ) :
            udist1.draw()
            
    udist2 = run_USTAT ( model , data_g , r , logger )
    if udist2 :
        keep.add ( udist2 ) 
        with use_canvas ( 'test_good_fit_2: USTAT' , wait = 1 ) :
            udist2.draw()
    
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
        with timing ( 'GoF1D-toys' , logger = logger ) : 
            toys = G1D.GoF1DToys ( gof  )
            toys = toys.run ( 500 )
        logger.info ( 'Goodness-of-fit with %d toys:\n%s' % ( toys.nToys , toys ) ) 

    del gof
    del toys
    
    ## Try to use multidimensional methods
    run_PPD ( model , data_b , r , logger )
    udist1 = run_DNN ( model , data_b , r , logger )
    if udist1 :
        keep.add ( udist1 ) 
        with use_canvas ( 'test_good_fit_3: DNN' , wait = 1 ) :
            udist1.draw()
    udist2 = run_USTAT ( model , data_b , r , logger )
    if udist2 :
        keep.add ( udist2 ) 
        with use_canvas ( 'test_good_fit_3: USTAT' , wait = 1 ) :
            udist2.draw()
            
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
        with timing ( 'GoF1D-toys' , logger = logger ) : 
            toys = G1D.GoF1DToys ( gof )
            toys = toys.run ( nToys = 1000 ) 
        logger.info ( 'Goodness-of-fit with %d toys:\n%s' % ( toys.nToys , toys ) ) 
    
    with use_canvas ( 'test_bad_fit_1: GoF/Kolmogorov-Smirnov' , wait = 1 ) :
        dks = toys.draw('Kolmogorov-Smirnov')
    with use_canvas ( 'test_bad_fit_1: GoF/Anderson-Darling'   , wait = 1 ) :
        dad = toys.draw('Anderson-Darling')
    with use_canvas ( 'test_bad_fit_1: GoF/Cramer-von Mises'   , wait = 1 ) :
        dcm = toys.draw('Cramer-von Mises')
    with use_canvas ( 'test_bad_fit_1: GoF/Kuiper'             , wait = 1 ) :
        dcm = toys.draw('Kuiper')
    with use_canvas ( 'test_bad_fit_1: GoF/ZK' , wait = 1 ) :
        dzk = toys.draw('ZK')
    with use_canvas ( 'test_bad_fit_1: GoF/ZA' , wait = 1 ) :
        dza = toys.draw('ZA')
    with use_canvas ( 'test_bad_fit_1: GoF/ZC' , wait = 1 ) :
        dzc = toys.draw('ZC')
    
    ## Try to use multidimensional methods
    run_PPD ( gauss , data_b , r , logger )
    udist1 = run_DNN ( gauss , data_b , r , logger )
    if udist1 :
        keep.add ( udist1 ) 
        with use_canvas ( 'test_bad_fit_1: DNN' , wait = 1 ) : udist1.draw()
        
    udist2 = run_USTAT ( gauss , data_b , r , logger )
    if udist2 :
        keep.add ( udist2 ) 
        with use_canvas ( 'test_bad_fit_3: USTAT' , wait = 1 ) : udist2.draw()
            
    del gof
    del toys
    
# ===============================================================================
if '__main__' == __name__ :

    test_good_fit_1 ()  ## fit Gauss       to Gauss 
    ## test_good_fit_2 ()  ## fit Gauss+Bkg   to Gauss 
    ## test_good_fit_3 ()  ## fit Gauss+Bkg   to Gauss+Bkg
    ## test_bad_fit_1  ()  ## fit Gauss       to Gauss+Bkg

# ===============================================================================
##                                                                        The END 
# ===============================================================================
                   
                   
