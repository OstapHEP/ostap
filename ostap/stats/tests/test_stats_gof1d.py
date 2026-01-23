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
from   ostap.plotting.canvas  import use_canvas
from   ostap.utils.ranges     import vrange
from   ostap.utils.root_utils import batch_env
from   ostap.utils.basic      import numcpu 
from   ostap.utils.timing     import timing 
from   ostap.core.core        import VE
from   ostap.stats.ustat      import USTAT 
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

xvar  = ROOT.RooRealVar ( 'x', '' ,  0 , 20 )
gauss = M.Gauss_pdf     ( 'G' , xvar = xvar , mean = 10 , sigma = 1 )
model = M.Fit1D         ( signal = gauss , background = 'flat' , fix_norm = True )

ND1   = 100
ND2   = 100

# =============================================================================
## data_g: pure gaussian
gauss.mean  = 10
gauss.sigma.fix ( 1.0 )

data_g1     = gauss.generate ( ND1 , sample = True )
data_g2     = gauss.generate ( ND2 , sample = True )

# =============================================================================
## data_b: Gaussian + flat background 
model.S     = 0.0 
model.B     = ND2
data_u      = model.generate ( ND2 , sample = True )

data_g = data_g1 + data_g2
data_b = data_g1 + data_u 

fitconf = { 'draw' : True , 'nbins' : 25 , 'refit' : 5 , 'quiet' : True }

# ==============================================================================
## Run Point-to-Point dissimilatity Goodness-of-Fit test
def run_PPD ( pdf, data, result , label , logger = logger ) :
    """ Run Point-to-Point dissimilarity Goodness-of-Fit test
    """

    nToys = 10 if numcpu () < 10 else 200 
    Ns    =  5 if numcpu () < 10 else  10
    
    logger.info ( 'Run Point-to-Point Dissimilarity GoF-test %s for %d different values of sigma' % ( label , Ns  ) )

    header = ()
    rows   = [] 
    for sigma in vrange ( 0.01 , 1.5 , Ns ) :

        with timing ( 'PPD-test %s sigma=%.3f' %  ( label  , sigma ) , logger = logger ) :
            
            ppd = GnD.PPD ( nToys = nToys , sigma = sigma )
            pdf.load_params ( result , silent = True ) 
            tvalue         = ppd          ( pdf , data )
            tvalue, pvalue = ppd.pvalue   ( pdf , data )
            
            ecdf           = ppd.ecdf
            
            header , row   = ppd.the_row ( tvalue = tvalue ,
                                           pvalue = pvalue ,
                                           ecdf   = ecdf   )
            
            row = ( '%.4f' % sigma , ) + tuple ( row ) 
            rows.append ( row ) 

    header = ( 'PPD(sigma)' , ) + header
    rows   = [ header ] + rows 
    rows   = T.remove_empty_columns ( rows ) 
    title  = 'Goodness-of-fit PPD-test %s (Gaussian with various sigmas)' % label 
    table  = T.table ( rows , title = title , prefix = '# ' , alignment = 'lccccccc' )
    logger.info ( '%s:\n%s' % ( title , table ) )

# ==============================================================================
## Run Distance-to-Nearest-Neighbour Goodness-of-Fit test
def run_DNN  ( pdf , data , result , label , logger = logger ) :
    """ Run Distance-to-Nearest-Neighbour Goodness-of-Fit test
    """

    nToys = 10 if numcpu () < 10 else 200 
    
    with timing ( 'DNN-test %s' % label , logger = logger ) :
        
        dnn = GnD.DNN   ( nToys = nToys , histo = 51 )        
        pdf.load_params ( result , silent = True )
        
        tvalue         = dnn          ( pdf , data )
        tvalue, pvalue = dnn.pvalue   ( pdf , data )

    with use_canvas ( 'DNN-toys %s'  % label ) : dnn.draw ( tvalue = tvalue ) 

    title  = 'Goodness-of-Fit DNN-test %s' % label     
    table  = dnn.table ( title = title , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )

    ## u-distribution
    return dnn.histo 

# ==============================================================================
## Run USTAT Goodness-of-Fit test
def run_USTAT  ( pdf , data, result , label , logger = logger ) :
    """ Run USTAT Goodness-of-Fit test
    """

    nToys = 10 if numcpu () < 10 else 200 

    with timing ( 'uStat-test' , logger = logger ) : 
        ustat = USTAT   ( nToys = nToys , histo = 110 , parallel = True  )        
        pdf.load_params ( result , silent = True )
        tvalue         = ustat        ( pdf , data )    
        tvalue, pvalue = ustat.pvalue ( pdf , data )
        
    with use_canvas ( 'USTAT-toys %s'  % label ) :
        ustat.draw ( tvalue = tvalue ) 
        
    title  = 'Goodness-of-Fit USTAT-test %s' % label     
    table  = ustat.table ( title = title , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )

    ## u-distribution
    return ustat.histo 

# ==============================================================================
def run_fit ( pdf , dataset , label  , logger = logger ) :
    """ Make a test for presumably GOOD/BAD fits
    """

    nToys = 50 if numcpu () < 10 else 500
    
    logger.info ( 'Make a test for presumably %s fit' % label  )

    with use_canvas ( '%s_fit_1' % label ) :
        r , f = pdf.fitTo ( dataset, **fitconf ) 

    with use_canvas ( '%s_fit_1: GoF' % label   ) :
        
        gauss.load_params ( r , silent = True ) 
        gof = G1D.GoF1D   ( pdf , data_g )
        logger.info ( 'Goodness-of-fit (%s):\n%s' %  ( label , gof ) )
        
        gauss.load_params ( r , silent = True ) 
        with timing ( 'GoF1D-toys %s' % label  , logger = logger ) : 
            toys = G1D.GoF1DToys ( gof )
            toys = toys.run ( nToys = nToys , parallel = True )
        logger.info ( 'Goodness-of-fit (%s) with %d toys:\n%s' % ( label , toys.nToys , toys ) ) 

    with use_canvas ( '%s_fit_1: GoF/Kolmogorov-Smirnov' % label ) :
        dks = toys.draw ( 'Kolmogorov-Smirnov')
    with use_canvas ( '%s_fit_1: GoF/Anderson-Darling'   % label ) :
        dad = toys.draw ( 'Anderson-Darling')
    with use_canvas ( '%s_fit_1: GoF/Cramer-von Mises'   % label ) :
        dcm = toys.draw ( 'Cramer-von Mises')
    with use_canvas ( '%s_fit_1: GoF/Kuiper'             % label ) :
        dcm = toys.draw ( 'Kuiper')
    with use_canvas ( '%s_fit_1: GoF/ZK' % label ) :
        dzk = toys.draw ( 'ZK')
    with use_canvas ( '%s_fit_1: GoF/ZA' % label ) :
        dza = toys.draw ( 'ZA')
    with use_canvas ( '%s_fit_1: GoF/ZC' % label ) :
        dzc = toys.draw ( 'ZC')
    
    ## Try to use multidimensional methods:

    run_PPD   ( pdf , dataset , r , label , logger )    
    run_DNN   ( pdf , dataset , r , label , logger )
    run_USTAT ( pdf , dataset , r , label , logger )

# =====================================================================================
def test_good_fit_1 ( ) :
    logger = getLogger ( 'test_GOOD_fit_1' )
    return run_fit ( gauss , data_g , 'GOOD1' , logger )

# =====================================================================================
def test_good_fit_2 ( ) :
    logger = getLogger ( 'test_GOOD_fit_2' )
    return run_fit ( model , data_b , 'GOOD2' , logger )

# =====================================================================================
def test_bad_fit_1  ( ) :
    logger = getLogger ( 'test_BAD_fit_1'  )
    return run_fit ( gauss , data_b , 'BAD'  , logger  )

# ===============================================================================
if '__main__' == __name__ :

    test_good_fit_1 ()  ## fit Gauss       to Gauss     data
    test_good_fit_2 ()  ## fit Gauss+Bkg   to Gauss+Bkg data 
    test_bad_fit_1  ()  ## fit Gauss       to Gauss+Bkg daat 

# ===============================================================================
##                                                                        The END 
# ===============================================================================
                   
                   
