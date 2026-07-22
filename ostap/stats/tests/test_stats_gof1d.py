#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==============================================================================
# @file ostap/stats/tests/test_stats_gof1d.py 
# Test Goodness-of-fits for 1D fits 
# Copyright (c) Ostap developpers.
# ==============================================================================
""" Test Goodness-of-fits for 1D fits 
"""
# ============================================================================== 
## ATTENTION! 
import os 
os.environ [ "OMP_NUM_THREADS"      ]  = "1"
os.environ [ "MKL_NUM_THREADS"      ]  = "1"
os.environ [ "OPENBLAS_NUM_THREADS" ]  = "1"
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

ND1   = 200
ND2   = 200

# =============================================================================
## data_g: pure gaussian
gauss.mean .fix ( 10  ) 
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

gauss.mean .release () 
gauss.sigma.fix ( 1.0 )

fitconf = { 'draw' : True , 'nbins' : 25 , 'refit' : 5 , 'quiet' : True }

small   = numcpu () <= 8

# ==============================================================================
## Run Mix-Sample Goodness-of-Fit test
def run_MIX ( pdf, data, result , label , logger = logger ) :
    """ Run Mixed-Sample Goodness-of-Fit test
    """

    nToys = 10 if small else 200 

    logger.info ( 'Run Mixed-Sample GoF-test %s' % label )

    with timing ( 'MIX-test %s' %  label  , logger = logger ) :
        
        gof = GnD.MIX ( nToys = nToys )
        pdf.load_params ( result , silent = True ) 
        tvalue         = gof          ( pdf , data )
        tvalue, pvalue = gof.pvalue   ( pdf , data , tvalue = tvalue )
    
    title  = 'Goodness-of-Fit Mixed-Samples-test %s' % label     
    table  = gof.report ( title = title , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )

# ==============================================================================
## Run Point-to-Point dissimilatity Goodness-of-Fit test
def run_PPD ( pdf, data, result , label , logger = logger ) :
    """ Run Point-to-Point dissimilarity Goodness-of-Fit test
    """

    nToys = 10 if small else 200     
    Ns    =  2 if small else   5

    Ns = 1
    logger.info ( 'Run Point-to-Point Dissimilarity GoF-test %s for %d different values of sigma' % ( label , Ns + 1   ) )

    header = ()
    rows   = [] 
    for sigma in vrange ( 0.01 , 1.5 , Ns ) :

        with timing ( 'PPD-test %s sigma=%.3f' %  ( label  , sigma ) , logger = logger ) :
            
            gof = GnD.PPD ( nToys = nToys , sigma = sigma )
            pdf.load_params ( result , silent = True ) 
            tvalue         = gof          ( pdf , data )
            tvalue, pvalue = gof.pvalue   ( pdf , data , tvalue = tvalue )            
            ecdf           = gof.ecdf            
            header , row   = gof.the_row ( tvalue = tvalue ,
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

    nToys = 10 if small else 200     
    
    with timing ( 'DNN-test %s' % label , logger = logger ) :
        
        gof            = GnD.DNN   ( nToys = nToys , histo = 51 )        
        pdf.load_params ( result , silent = True )        
        tvalue         = gof          ( pdf , data )
        tvalue, pvalue = gof.pvalue   ( pdf , data , tvalue = tvalue )

    with use_canvas ( 'DNN-toys %s'  % label ) : gof.draw ( tvalue = tvalue ) 

    title  = 'Goodness-of-Fit DNN-test %s' % label     
    table  = gof.report ( title = title , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )

# ==============================================================================
## Run USTAT Goodness-of-Fit test
def run_USTAT  ( pdf , data, result , label , logger = logger ) :
    """ Run USTAT Goodness-of-Fit test
    """

    nToys = 10 if small else 200     

    with timing ( 'uStat-test %s' % label , logger = logger ) : 
        gof = USTAT   ( nToys = nToys , histo = 110 , parallel = True  )        
        pdf.load_params ( result , silent = True )
        tvalue         = gof        ( pdf , data )    
        tvalue, pvalue = gof.pvalue ( pdf , data , tvalue = tvalue )
        
    ##with use_canvas ( 'USTAT-toys %s'  % label ) :
    ##    ustat.draw ( tvalue = tvalue ) 
        
    title  = 'Goodness-of-Fit USTAT-test %s' % label     
    table  = gof.report ( title = title , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )

plots = [] 
# ==============================================================================
def run_fit ( pdf , dataset , label  , logger = logger ) :
    """ Make a test for presumably GOOD or BAD fits
    """

    nToys = 10 if small else 200     

    nToys = 20

    
    logger.info ( 'Make a test for presumably %s fit' % label  )

    with use_canvas ( title = '%s:Fit' % label ) :
        r , f = pdf.fitTo ( dataset, **fitconf ) 
        plots.append ( f )
        
    with use_canvas ( title = '%s:GoF' % label   ) :
        
        pdf.load_params ( r , silent = True ) 
        gof = G1D.GoF1D   ( pdf , dataset , fitresult = r )
        logger.info ( 'Goodness-of-fit (%s) #%d:\n%s' %  ( label , len ( dataset ) , gof ) )
        
        gauss.load_params ( r , silent = True ) 
        with timing ( 'GoF1D-toys %s' % label  , logger = logger ) : 
            toys = G1D.GoF1DToys ( gof )
            toys = toys.run ( nToys = nToys , parallel = True )
        logger.info ( 'Goodness-of-fit (%s) with %d toys:\n%s' % ( label , toys.nToys , toys ) ) 

    """
    with use_canvas ( title = '%s:GoF/Kolmogorov-Smirnov' % label ) :
        dks = toys.draw ( 'Kolmogorov-Smirnov')
        plots.append ( dks ) 
    with use_canvas ( title = '%s:GoF/Anderson-Darling'   % label ) :
        dad = toys.draw ( 'Anderson-Darling')
        plots.append ( dad ) 
    with use_canvas ( title = '%s:GoF/Cramer-von Mises'   % label ) :
        dcm = toys.draw ( 'Cramer-von Mises')
        plots.append ( dcm ) 
    with use_canvas ( title = '%s:GoF/Kuiper'             % label ) :
        dku = toys.draw ( 'Kuiper')
        plots.append ( dku ) 
    with use_canvas ( title = '%s:GoF/ZK' % label ) :
        dzk = toys.draw ( 'ZK')
        plots.append ( dzk )         
    with use_canvas ( '%s_fit_1:GoF/ZA' % label ) :
        dza = toys.draw ( 'ZA')
        plots.append ( dza )         
    with use_canvas ( title = '%s:GoF/ZC' % label ) :
        dzc = toys.draw ( 'ZC')
        plots.append ( dzc )         
    with use_canvas ( title = '%s:GoF/Berk-Jones'         % label ) :
        dbj = toys.draw ( 'Berk-Jones')
        plots.append ( dbj )         
    with use_canvas ( title = '%s:GoF/NLL'                 % label ) :
        dnl = toys.draw ( 'BLL')
        plots.append ( dnl )         
    with use_canvas ( title = '%s:GoF/AikaikeIC'           % label ) :
        dai = toys.draw ( 'AIC')
        plots.append ( dai )         
    with use_canvas ( title = '%s:GoF/BayesianIC'          % label ) :
        dbi = toys.draw ( 'BIC')
        plots.append ( dbi )
    
    """
    
    ## Try to use multidimensional methods:
    
    ## run_MIX   ( pdf , dataset , r , label , logger )    
    ## run_PPD   ( pdf , dataset , r , label , logger )    
    ## run_DNN   ( pdf , dataset , r , label , logger )
    run_USTAT ( pdf , dataset , r , label , logger )

# =====================================================================================
def test_good_fit_1 ( ) :
    logger = getLogger ( 'test_GOOD: ( G -> G )' )
    return run_fit ( gauss , data_g , 'GOOD1' , logger )

# =====================================================================================
def test_good_fit_2 ( ) :
    logger = getLogger ( 'test_GOOD: ( G + B -> G + B )' )
    return run_fit ( model , data_b , 'GOOD2' , logger )

# =====================================================================================
def test_bad_fit_1  ( ) :
    logger = getLogger ( 'test_BAD: ( G -> G + B )'  )
    return run_fit ( gauss , data_b , 'BAD'  , logger  )

# ===============================================================================
if '__main__' == __name__ :

    test_good_fit_1 ()  ## fit Gauss       to Gauss     data
    test_good_fit_2 ()  ## fit Gauss+Bkg   to Gauss+Bkg data 
    test_bad_fit_1  ()  ## fit Gauss       to Gauss+Bkg daat 

# ===============================================================================
##                                                                        The END 
# ===============================================================================
                   
                   
