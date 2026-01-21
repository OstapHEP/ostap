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
from   ostap.stats.gof_utils  import clip_pvalue
from   ostap.logger.symbols   import plus_minus, times, greek_lower_sigma
from   ostap.math.ve          import fmt_pretty_ve
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
## data_b: 
model.S     = 0.0 
model.B     = ND2
data_u      = model.generate ( ND2 , sample = True )

data_g = data_g1 + data_g2
data_b = data_g1 + data_u 

fitconf = { 'draw' : True , 'nbins' : 25 , 'refit' : 5 , 'quiet' : True }

keep = set() 
# ==============================================================================
## Run Point-to_pint dissimilatity Goodness-of-Fit test
def run_PPD ( pdf , data, result , logger ) :
    """ Run Point-to_point dissimilarity Goodness-of-Fit test
    """

    # =========================================================================
    #  - Point to Point Dissimilarity test  with Gaussian distance using different "sigma"
    rows  = [ ( 'PPD(sigma)' ,
                't-value'    ,
                't-mean'     ,
                't-rms'      ,
                't-min/max'  ,                
                '%s[..]' % times , 'p-value [%]' , '#%s' % greek_lower_sigma ) ]
    Ns    = 10 
    logger.info ( 'Run Point-to-Point Dissimilarity GoF-test for %d different values of sigma' % Ns  ) 

    header = ()
    rows   = [] 
    for sigma in vrange ( 0.01 , 1.5 , Ns ) :

        with timing ( 'PPD-test' , logger = logger ) :
            
            ppd = GnD.PPD ( nToys = 500 , sigma = sigma )
            pdf.load_params ( result , silent = True ) 
            tvalue         = ppd          ( pdf , data )
            tvalue, pvalue = ppd.pvalue   ( pdf , data )
            
            ecdf           = ppd.ecdf
            
            header , row = ppd.the_row ( tvalue = tvalue ,
                                         pvalue = pvalue ,
                                         ecdf   = ecdf   )

            row = ( '%.4f' % sigma , ) + tuple ( row ) 
            rows.append ( row ) 

    header = ( 'PPD(sigma)' , ) + header
    rows   = [ header ] + rows 
    rows   = T.remove_empty_columns ( rows ) 
    title  = 'Goodness-of-fit PPD-test (Gaussian with various sigmas)'
    table  = T.table ( rows , title = title , prefix = '# ' , alignment = 'lccccccc' )
    logger.info ( '%s:\n%s' % ( title , table ) )

import ostap.stats.ustat as U 

# ==============================================================================
## Run Distance-to-Nearest-Neighbour Goodness-of-Fit test
def run_DNN  ( pdf , data , result , logger ) :
    """ Run Distance-to-Nearest-Neighbour Goodness-of-Fit test
    """


    with timing ( 'DNN-test' , logger = logger ) :
        
        dnn = GnD.DNN ( nToys = 200 , histo = 50 )        
        pdf.load_params ( result , silent = True )
        
        tvalue         = dnn          ( pdf , data )
        tvalue, pvalue = dnn.pvalue   ( pdf , data )

    with use_canvas ( 'DNN-toys' ) : dnn.draw ( tvalue = tvalue ) 

    title  = 'Goodness-of-fit DNN-test'    
    table  = dnn.table ( title = title , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )

    ## u-distribution
    return dnn.histo 

# ==============================================================================
## Run USTAT Goodness-of-Fit test
def run_USTAT  ( pdf , data, result , logger ) :
    """ Run USTAT Goodness-of-Fit test
    """

    rows  =  [ ( 't-value'  , '%s[..]' % times , 'p-value [%]' , '#%s' % greek_lower_sigma ) ]
    
    ustat = USTAT ( nToys = 200 , histo = 100 )
    
    pdf.load_params ( result , silent = True )
    
    with timing ( 'uStat-test' , logger = logger ) : 
        tvalue         = ustat        ( pdf , data )    
        tvalue, pvalue = ustat.pvalue ( pdf , data )
        
        pv     = clip_pvalue  ( pvalue ) 
        nsigma = significance ( pv ) ## convert  it to significace

    tv , texpo = pretty_float ( tvalue )
    pvalue *= 100
    pvalue  = '%5.2f %s %.2f' % ( pvalue.value() , plus_minus , pvalue.error() )
    nsigma  = '%.2f %s %.2f'  % ( nsigma.value() , plus_minus , nsigma.error () )    
    row     = tv , '10^%+d' % texpo if texpo else '' , pvalue, nsigma 
    rows.append ( row )
    
    title = 'Goodness-of-Fit USTAT-test'
    rows  = T.remove_empty_columns ( rows ) 
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lccccccc' )
    logger.info ( '%s:\n%s' % ( title , table ) )

    return ustat.histo


# ==============================================================================
def run_fit ( pdf , dataset , label  , logger = logger ) :
    """ Make a test for presumably GOO/BAD fit: fit Gauss to Gauss
    """
    
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
            toys = toys.run ( nToys = 500 , parallel = True )
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
    
    ## Try to use multidimensional methods
    run_PPD ( pdf , dataset , r , logger )
    
    udist1 = run_DNN    ( pdf , dataset , r , logger )
    if udist1 :
        keep.add ( udist1 ) 
        with use_canvas ( '%s_fit_1: DNN' % label ) : udist1.draw()

    udist2 = run_USTAT  ( pdf , dataset , r , logger )
    if udist2 :
        keep.add ( udist2 ) 
        with use_canvas ( '%s_fit_1: USTAT' % label ) : udist2.draw()
    
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
                   
                   
