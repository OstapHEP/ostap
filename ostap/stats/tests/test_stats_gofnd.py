#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==============================================================================
# @file ostap/stats/tests/test_stats_gof1d.py 
# Test Goddness-of-fits 1D 
# Copyright (c) Ostap developpers.
# ==============================================================================
""" # Test averages for inconsistend data 
"""
# ==============================================================================
## ATTENTION! 
import os 
os.environ["OMP_NUM_THREADS"     ]  = "1"
os.environ["MKL_NUM_THREADS"     ]  = "1"
os.environ["OPENBLAS_NUM_THREADS"]  = "1"
# ==============================================================================
from   packaging.version      import Version
from   ostap.plotting.canvas  import use_canvas
from   ostap.utils.timing     import timing
from   ostap.logger.pretty    import pretty_float
from   ostap.plotting.canvas  import use_canvas
from   ostap.math.math_ve     import significance
from   ostap.utils.root_utils import batch_env
from   ostap.utils.basic      import numcpu, typename 
from   ostap.logger.symbols   import plus_minus , greek_lower_sigma 
from   ostap.stats.gof_utils  import clip_pvalue, useLightGBM, useXGBoost, useCatBoost 
import ostap.fitting.models   as     M 
import ostap.stats.gofnd      as     GnD
import ostap.logger.table     as     T 
import ROOT, random   
# ==============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'tests_stats_gofnd' )
else                       : logger = getLogger ( __name__            )
# ==============================================================================
batch_env ( logger ) 
# =============================================================================

xvar   = ROOT.RooRealVar ( 'x', '', 0, 10)
yvar   = ROOT.RooRealVar ( 'y', '', 0, 10)
varset = ROOT.RooArgSet  ( xvar , yvar   )

xgauss = M.Gauss_pdf     ( 'GX' , xvar = xvar , mean = ( 5 , 4 , 6 ) , sigma = ( 0.5 , 0.1 , 2.5 ) )
ygauss = M.Gauss_pdf     ( 'GY' , xvar = yvar , mean = ( 5 , 4 , 6 ) , sigma = ( 0.5 , 0.1 , 2.5 ) )
gauss2 = xgauss*ygauss

NG        = 100 
NG2       = 100
data_good = gauss2.generate ( NG + NG2 , sample = False )
data_bad  = gauss2.generate ( NG       , sample = False )
for i in range ( NG2 ) :
    x = random.uniform ( 0 , 10 )
    y = random.uniform ( 0 , 10 )    
    xvar.setVal ( x )
    yvar.setVal ( y )
    data_bad.add ( varset )

xgauss.sigma.fix()
ygauss.sigma.fix()
    
pdf       = gauss2 
rgood , _ = pdf.fitTo  ( data_good , quiet = True , refit = 5 )
rbad  , _ = pdf.fitTo  ( data_bad  , quiet = True , refit = 5 )

# ==============================================================================
use_lightgbm = useLightGBM  ()
use_xgboost  = useXGBoost ()
use_catboost = useCatBoost ()
if use_lightgbm : logger.attention ( 'USE LigthGBM!'               )
else            : logger.warning   ( 'LightGBM is not available!'  )
if use_xgboost  : logger.attention ( 'USE XGBoost!'                )
else            : logger.warning   ( 'XGBoost  is not available!'  )
if use_catboost : logger.attention ( 'USE CatBoost!'               )
else            : logger.warning   ( 'CatBoost  is not available!' )

keep_it = [] 
# ===============================================================================
def probe_GOF ( gof1 , gof2, tag ) :

    logger = getLogger ("test_GOF/%s" % tag )

    if gof1 is gof2 :  logger.info ( 'Testing %s (GOOD&BAD) \n%s' % ( tag , gof1.table ( prefix = '# ' ) ) )
    
    ## presumably good fit
    with timing ( "Good fit GOF/%s" % tag , logger = logger ) :
        if not gof1 is gof2 : logger.info ( 'Testing %s (GOOD)\n%s' % ( tag , gof1.table ( prefix = '# ' ) ) )
        pdf.load_params ( rgood , silent = True ) 
        tgood, pgood = gof1.pvalue ( pdf , data_good )        
        logger.info ( 'Report %s (GOOD)\n%s' % ( tag , gof1.report ( prefix = '# ' ) ) )
        
    with use_canvas ( 'GoF-toys %s/%s (GOOD)'  % ( typename ( gof1 ) , tag  ) ) :
        keep_it.append ( gof1.draw ( tvalue = tgood ) )
        
    ## presumably bad fit 
    with timing ( "Bad  fit GOF/%s" % tag  , logger = logger ) : 
        if not gof1 is gof2 : logger.info ( 'Testing %s (BAD)\n%s' % ( tag , gof2.table ( prefix = '# ' ) ) )
        pdf.load_params ( rbad  , silent = True ) 
        tbad, pbad  = gof2.pvalue ( pdf , data_bad )
        logger.info ( 'Report %s (BAD)\n%s' % ( tag , gof2.report ( prefix = '# ' )  ) )

    with use_canvas ( 'GoF-toys %s/%s (BAD)'  % ( typename ( gof2 ) , tag  ) ) :
        keep_it.append ( gof2.draw ( tvalue = tbad ) ) 
    
    gp = pgood * 100 
    bp = pbad  * 100
    
    gt , ge = pretty_float ( tgood )
    bt , be = pretty_float ( tbad  )
    
    pvg     = clip_pvalue  ( pgood )
    pvb     = clip_pvalue  ( pbad  )        
    nsg     = significance ( pvg   )
    nsb     = significance ( pvb   )
    
    nsg    = '%.2f %s %.2f' % ( nsg.value() , plus_minus , nsg.error () )
    nsb    = '%.2f %s %.2f' % ( nsb.value() , plus_minus , nsb.error () )
    
    return tag , '%5.2f %s %.2f' % ( gp.value() , plus_minus , gp.error () ) , \
        '%5.2f %s %.2f'   % ( bp.value() , plus_minus , bp.error () ) , nsg , nsb


# ===============================================================================
def test_GOF () :
    
    to_test = []

    small = numcpu () <= 8

    ## small = True
    
    ## very small number of toys     

    nToys =  40 if small else 500

    ## nToys = 400 
    
    tconf = { 'nToys' : nToys , 'parallel' : True }
    
    # ===========================================================================
    ## MIX 
    # ===========================================================================
    mix   = GnD.MIX  ( **tconf )
    entry = mix , mix, 'MIX-Sample'
    to_test.append ( entry )

    """ 
    for conf in ( { 'psi' : 'linear'     } ,
                  { 'psi' : 'logarithm'  } ,
                  { 'psi' : 'chebyshev'  } ,   
                  { 'psi' : 'coulomb'    } ,
                  ##
                  ## { 'psi' : 'inverse2'   } ,   
                  ## { 'psi' : 'squared'    } ,   
                  ## { 'psi' : 'cosine'     } ,   
                  ## { 'psi' : 'canberra'   } ,
                  ## 
                  { 'psi' : 'braycurtis' } ,   
                  { 'psi' : 'cityblock'  } ,
                  ## 
                  { 'psi' : 'gaussian' , 'sigma' : 5.00 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 2.00 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 1.00 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.50 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.10 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.05 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.02 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.01 }
                 ) : 
        
        tag = 'PPD:%s' % conf['psi']
        if 'sigma' in conf : tag = '%s/%s=%s' % ( tag , greek_lower_sigma , conf['sigma'] ) 
        conf.update ( tconf )        
        ppd = GnD.PPD ( **conf )
        entry = ppd , ppd , tag
        to_test.append ( entry )
        
    # ===========================================================================
    ## DNN 
    # ===========================================================================
    dnn   = GnD.DNN ( nToys = nToys , histo = 50  )
    entry = dnn , dnn , 'DNN'
    to_test.append ( entry )

    # ===========================================================================
    ## U-stat 
    # ===========================================================================
    from ostap.stats.ustat import USTAT 
    ust   = USTAT ( nToys = nToys  , histo = 50 )
    entry = ust , ust , 'U-stat'
    to_test.append ( entry )
    """
    
    # ===========================================================================
    ## ADVAL-based tests:
    # ===========================================================================
    
    config  = { 'nToys' : nToys , 'parallel' : True , 'progress' : True }

    # ===========================================================================
    if use_lightgbm : # =====================================================
        # ===================================================================
        import lightgbm
        from   ostap.stats.gofnd import ADVAL_LightGBM as GOF 
        ## 
        gof    = GOF ( **config )
        test   = gof , gof , 'ADVAL:LightGBM'
        to_test.append ( test )

    # ===========================================================================
    if use_xgboost : # ======================================================
        # ===================================================================
            from   ostap.stats.gofnd import ADVAL_XGBoost as GOF 
            ## 
            gof    = GOF ( **config )
            test   = gof , gof , 'ADVAL:XBGoost'
            to_test.append ( test )             
            
    # ===========================================================================
    if use_catboost : 
        # ===================================================================
        import catboost
        from   ostap.stats.gofnd import ADVAL_CatBoost as GOF 
        ## 
        gof    = GOF ( **config )
        test   = gof , gof , 'ADVAL:CatBoost'
        to_test.append ( test ) 
        
    # ===========================================================================
    from ostap.stats.gof_utils import has_sklearn # =============================
    # ===========================================================================
    
    """ 
    # ===========================================================================
    try : # =====================================================================
        # =======================================================================
        if  has_sklearn :
            ## 
            import sklearn
            from   ostap.stats.gofnd import ADVAL_HistoGBoost as GOF 
            ##            
            gof    = GOF ( **config )
            test   = gof , gof , 'ADVAL:HistGradientBoost'
            ## 
            if not small : to_test.append ( test )
        else :
            logger.warning ( "No scikit-learn is available" )            
        ## 
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        logger.error  ( "HistoGBoost is not available, skip the test" )
        
    # ==========================================================================    
    try : # ====================================================================
        # ======================================================================
        if  has_sklearn :
            ##
            import sklearn 
            from   ostap.stats.gofnd import ADVAL_GBoost as GOF 
            ## 
            gof    = GOF ( **config )
            test   = gof , gof , 'ADVAL:GradientBoost'
            ## 
            if not small : to_test.append ( test ) 
            ## 
        else :
            logger.warning ( "No scikit-learn is available" )            
        # =======================================================================
    except ImportError : # ======================================================
        # =======================================================================
        logger.error  ( "GBoost is not available, skip the test" )
        
    # ===========================================================================
    try : # =====================================================================
        # =======================================================================
        if  has_sklearn :
            ## 
            import sklearn
            from   ostap.stats.gofnd import ADVAL_RandomForest as GOF 
            ## 
            gof    = GOF (  **config   )
            test   = gof , gof , 'ADVAL:RandomForest'
            ## 
            if not small : to_test.append ( test ) 
            ##
        else :
            logger.warning ( "No scikit-learn is available" )            
        # =======================================================================
    except ImportError : # ======================================================
        # =======================================================================
        logger.error  ( "RandomForest is not available, skip the test" )

    """

    
    # ===========================================================================
    ## -log L 
    # ===========================================================================
    nll_good = GnD.NLL ( nToys = nToys , fitresult = rgood , parallel = True )
    nll_bad  = GnD.NLL ( nToys = nToys , fitresult = rbad  , parallel = True )
    entry    = nll_good , nll_bad , '-log L'
    to_test.append ( entry )

    # ===========================================================================
    ## Aikaike IC 
    # ===========================================================================
    aic_good = GnD.AikaikeIC ( nToys = nToys , fitresult = rgood , parallel = True )
    aic_bad  = GnD.AikaikeIC ( nToys = nToys , fitresult = rbad  , parallel = True )
    entry    = aic_good , aic_bad , 'Aikaike IC'
    to_test.append ( entry )

    # ===========================================================================
    ## Bayesial IC 
    # ===========================================================================
    bic_good = GnD.BayesianIC ( nToys = nToys , fitresult = rgood , data = data_good , parallel = True )
    bic_bad  = GnD.BayesianIC ( nToys = nToys , fitresult = rbad  , data = data_bad  , parallel = True )
    entry    = bic_good , bic_bad , 'Bayesian IC'
    to_test.append ( entry )
    
    # ===========================================================================
    ## Run the test and build the table
    # ===========================================================================
    
    rows  = [ ( 'Estimator' ,
                'GOOD-FIT:p-value[%]' ,
                'BAD-FIT:p-value[%]'  ,
                'GOOD-FIT:#%s' % greek_lower_sigma ,
                'BAD-FIT:#%s'  % greek_lower_sigma , 'time [s]' ) ] 

    for gof1, gof2, tag  in to_test :
        with timing ( 'Processing GoF/%s' % tag , logger = logger ) as timer :
            row = probe_GOF ( gof1 , gof2 , tag  )
        row = list ( row ) + [ '%.1f' % timer.delta ] 
        rows.append ( row )
        
        title = 'Goodness-of-Fit tests'
        rows  = T.remove_empty_columns ( rows )     
        table = T.table ( rows , title = title , prefix = '# ')
        logger.info ( '%s:\n%s' % ( title , table ) )
        
# ===============================================================================
if '__main__' == __name__ :

    test_GOF  ()    

# ===============================================================================
##                                                                        The END 
# ===============================================================================
                   
    
