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
from   ostap.plotting.canvas  import use_canvas
from   ostap.utils.timing     import timing
from   ostap.logger.pretty    import pretty_float
from   ostap.plotting.canvas  import use_canvas
from   ostap.math.math_ve     import significance
from   ostap.utils.root_utils import batch_env
from   ostap.utils.basic      import numcpu, typename 
from   ostap.logger.symbols   import plus_minus , greek_lower_sigma 
import ostap.fitting.models   as     M 
import ostap.stats.gofnd      as     GnD
from   ostap.stats.gof_utils  import clip_pvalue 
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
## with use_canvas ( 'Good fit x-projection' ) : pdf.draw1 ( data_good , nbins = 50 ) 
## with use_canvas ( 'Good fit y-projection' ) : pdf.draw2 ( data_good , nbins = 50 ) 
rbad  , _ = pdf.fitTo  ( data_bad , quiet = True , refit = 5 )
## with use_canvas ( 'Bad  fit x-projection' ) : pdf.draw1 ( data_bad  , nbins = 50 ) 
## with use_canvas ( 'Bad  fit y-projection' ) : pdf.draw2 ( data_bad  , nbins = 50 ) 

# ===============================================================================
def test_PPD () :
    
    logger = getLogger ("test_PPD")

    
    sigma = ''

    if 10 <= numcpu () : tconf = { 'nToys' :   10 , 'parallel' : True }
    else               : tconf = { 'nToys' :   10 }

    to_test = []
    
    for conf in ( { 'psi' : 'linear'     } ,
                  ## { 'psi' : 'logarithm'  } ,
                  ## { 'psi' : 'chebyshev'  } ,   
                  ## { 'psi' : 'coulomb'    } ,
                  ##
                  ## ## { 'psi' : 'inverse2'   } ,   
                  ## ## { 'psi' : 'squared'    } ,   
                  ## ## { 'psi' : 'cosine'     } ,   
                  ## ## { 'psi' : 'canberra'   } ,
                  ## 
                  ## { 'psi' : 'braycurtis' } ,   
                  ## { 'psi' : 'cityblock'  } ,
                  ## 
                  ## { 'psi' : 'gaussian' , 'sigma' : 5.00 } ,
                  ## { 'psi' : 'gaussian' , 'sigma' : 2.00 } ,
                  ## { 'psi' : 'gaussian' , 'sigma' : 1.00 } ,
                  ## { 'psi' : 'gaussian' , 'sigma' : 0.50 } ,
                  ## { 'psi' : 'gaussian' , 'sigma' : 0.10 } ,
                  ## { 'psi' : 'gaussian' , 'sigma' : 0.05 } ,
                  ## { 'psi' : 'gaussian' , 'sigma' : 0.04 } ,
                  ## { 'psi' : 'gaussian' , 'sigma' : 0.03 } ,
                  ## { 'psi' : 'gaussian' , 'sigma' : 0.02 } ,
                  { 'psi' : 'gaussian' , 'sigma' : 0.01 } ) : 
        
        tag = 'PPD:%s' % conf['psi']
        if 'sigma' in conf : tag = '%s[%s=%s]' % ( tag , greek_lower_sigma , conf['sigma'] ) 
        conf.update ( tconf )        
        ppd = GnD.PPD ( **conf )
        entry = ppd ,tag
        to_test.append ( entry )

    rows  = [ ( 'Distance' , 'p-value/good[%]' , 'p-value/bad[%]' , '#%s/good' % greek_lower_sigma , '#%s/bad' % greek_lower_sigma , 'time' ) ]

    ## run tests: 
    for gof, tag  in to_test :
        with timing ( 'Processing GoF/%s' % tag , logger = logger ) as timer :
            row = probe_GOF ( gof , tag  )
        row = list ( row )  + [  '%.f' % timer.delta ] 
        rows.append ( row )

    print ( 'ROWS' , rows  ) 
    title = 'Goodness-of-Fit %s tests' % tag 
    rows  = T.remove_empty_columns ( rows )     
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )
    
    
# ===============================================================================
def test_DNN () :
    
    logger = getLogger ("test_DNN")

    rows  = [ ( 'p-value/good[%]' , 'p-value/bad[%]' , '#%s/good' % greek_lower_sigma , '#%s/bad' % greek_lower_sigma ) ]
    
    dnn = GnD.DNN ( nToys = 200 , histo = 50  )
    
    ## presumably good fit
    with timing ( "Good fit DNN" , logger = logger ) :
        pdf.load_params ( rgood , silent = True ) 
        tgood        = dnn        ( pdf , data_good )
        tgood, pgood = dnn.pvalue ( pdf , data_good )
        
    ## presumably bad fit 
    with timing ( "Bad  fit DNN" , logger = logger ) : 
        pdf.load_params ( rbad  , silent = True ) 
        tbad        = dnn        ( pdf , data_bad )
        tbad, pbad  = dnn.pvalue ( pdf , data_bad )
        
    gp = pgood * 100 
    bp = pbad  * 100
    
    gt , ge = pretty_float ( tgood )
    bt , be = pretty_float ( tbad  )
    
    pvg    = clip_pvalue  ( pgood )
    pvb    = clip_pvalue  ( pbad  )        
    nsg    = significance ( pvg   )
    nsb    = significance ( pvb   )
    
    nsg    = '%.2f %s %.2f' % ( nsg.value() , plus_minus , nsg.error () )
    nsb    = '%.2f %s %.2f' % ( nsb.value() , plus_minus , nsb.error () )
    
    row = '%5.2f %s %.2f' % ( gp.value() , plus_minus , gp.error () ) , \
        '%5.2f %s %.2f'   % ( bp.value() , plus_minus , bp.error () ) , nsg , nsb 
    rows.append ( row )
            
    title= 'Goodness-of-Fit DNN test'
    rows  = T.remove_empty_columns ( rows )     
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )

    
# ===============================================================================
def test_USTAT () :
    
    logger = getLogger ("test_USTAT")

    
    from ostap.stats.ustat import USTAT 
    rows  = [ ( 'p-value/good[%]' , 'p-value/bad[%]' , '#sigma/good' , '#sigma/bad') ]
    
    ust = USTAT ( nToys = 200  , histo = 50 )
    
    ## presumably good fit
    with timing ( "Good fit USTAT" , logger = logger ) :
        pdf.load_params ( rgood , silent = True ) 
        tgood        = ust        ( pdf , data_good )
        tgood, pgood = ust.pvalue ( pdf , data_good )
        
    ## presumably bad fit 
    with timing ( "Bad  fit USTAT" , logger = logger ) : 
        pdf.load_params ( rbad  , silent = True ) 
        tbad        = ust        ( pdf , data_bad )
        tbad, pbad  = ust.pvalue ( pdf , data_bad )
        
    gp = pgood * 100 
    bp = pbad  * 100
    
    gt , ge = pretty_float ( tgood )
    bt , be = pretty_float ( tbad  )

    pvg    = clip_pvalue  ( pgood )
    pvb    = clip_pvalue  ( pbad  )        
    nsg    = significance ( pvg   )
    nsb    = significance ( pvb   )

    nsg    = '%.2f %s %.2f' % ( nsg.value() , plus_minus , nsg.error () )
    nsb    = '%.2f %s %.2f' % ( nsb.value() , plus_minus , nsb.error () )
    
    row = '%5.2f %s %.2f' % ( gp.value() , plus_minus , gp.error () ) , \
        '%5.2f %s %.2f'   % ( bp.value() , plus_minus , bp.error () ) , nsg , nsb 
    rows.append ( row )
            
    title = 'Goodness-of-Fit USTAT test'
    rows  = T.remove_empty_columns ( rows )     
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )

# ===============================================================================
def test_NLL () :
    
    logger = getLogger ("test_NLL")

    rows  = [ ( 'p-value/good[%]' , 'p-value/bad[%]' , '#%s/good' % greek_lower_sigma , '#%s/bad' % greek_lower_sigma ) ]
    
    nll_good = GnD.NLL ( nToys = 50 , fitresult = rgood )
    nll_bad  = GnD.NLL ( nToys = 50 , fitresult = rbad  )
    
    ## presumably good fit
    with timing ( "Good fit NLL" , logger = logger ) :
        pdf.load_params ( rgood , silent = True ) 
        tgood, pgood = nll_good.pvalue ( pdf , data_good )
        
    ## presumably bad fit 
    with timing ( "Bad  fit NLL" , logger = logger ) : 
        pdf.load_params ( rbad  , silent = True ) 
        tbad, pbad  = nll_bad.pvalue ( pdf , data_bad )
        
    with use_canvas ( title = "Good fit -log N "  , wait = 2 ) : nll_good.draw ()     
    with use_canvas ( title = "Bad  fit -log N "  , wait = 2 ) : nll_bad .draw ()     
        
    gp = pgood * 100 
    bp = pbad  * 100
    
    gt , ge = pretty_float ( tgood )
    bt , be = pretty_float ( tbad  )
    
    pvg    = clip_pvalue  ( pgood )
    pvb    = clip_pvalue  ( pbad  )        
    nsg    = significance ( pvg   )
    nsb    = significance ( pvb   )
    
    nsg    = '%.2f %s %.2f' % ( nsg.value() , plus_minus , nsg.error () )
    nsb    = '%.2f %s %.2f' % ( nsb.value() , plus_minus , nsb.error () )
    
    row = '%5.2f %s %.2f' % ( gp.value() , plus_minus , gp.error () ) , \
        '%5.2f %s %.2f'   % ( bp.value() , plus_minus , bp.error () ) , nsg , nsb 
    rows.append ( row )
            
    title= 'Goodness-of-Fit -logL test'
    rows  = T.remove_empty_columns ( rows )     
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )

# ===============================================================================
def test_AIC () :
    
    logger = getLogger ("test_AIC")

    rows  = [ ( 'p-value/good[%]' , 'p-value/bad[%]' , '#%s/good' % greek_lower_sigma , '#%s/bad' % greek_lower_sigma ) ]
    
    aic_good = GnD.AikaikeIC ( nToys = 500 , fitresult = rgood )
    aic_bad  = GnD.AikaikeIC ( nToys = 500 , fitresult = rbad  )
    
    ## presumably good fit
    with timing ( "Good fit Aikaike IC" , logger = logger ) :
        pdf.load_params ( rgood , silent = True ) 
        tgood, pgood = aic_good.pvalue ( pdf , data_good )
                
    ## presumably bad fit 
    with timing ( "Bad  fit Aikaike IC" , logger = logger ) : 
        pdf.load_params ( rbad  , silent = True ) 
        tbad, pbad  = aic_bad.pvalue ( pdf , data_bad )
        
    with use_canvas ( title = "Good fit Aikaike IC"  , wait = 2 ) : aic_good.draw ()     
    with use_canvas ( title = "Bad  fit Aikaike IC"  , wait = 2 ) : aic_bad .draw ()     

    gp = pgood * 100 
    bp = pbad  * 100
    
    gt , ge = pretty_float ( tgood )
    bt , be = pretty_float ( tbad  )
    
    pvg    = clip_pvalue  ( pgood )
    pvb    = clip_pvalue  ( pbad  )        
    nsg    = significance ( pvg   )
    nsb    = significance ( pvb   )
    
    nsg    = '%.2f %s %.2f' % ( nsg.value() , plus_minus , nsg.error () )
    nsb    = '%.2f %s %.2f' % ( nsb.value() , plus_minus , nsb.error () )
    
    row = '%5.2f %s %.2f' % ( gp.value() , plus_minus , gp.error () ) , \
        '%5.2f %s %.2f'   % ( bp.value() , plus_minus , bp.error () ) , nsg , nsb 
    rows.append ( row )
            
    title= 'Goodness-of-Fit Aikaike IC test'
    rows  = T.remove_empty_columns ( rows )     
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )

# ===============================================================================
def test_BIC () :
    
    logger = getLogger ("test_BIC")

    rows  = [ ( 'p-value/good[%]' , 'p-value/bad[%]' , '#%s/good' % greek_lower_sigma , '#%s/bad' % greek_lower_sigma ) ]
    
    bic_good = GnD.BayesianIC ( nToys = 100 , fitresult = rgood , data = data_good , parallel = True )
    bic_bad  = GnD.BayesianIC ( nToys = 100 , fitresult = rbad  , data = data_bad  , parallel = True )
    
    ## presumably good fit
    with timing ( "Good fit Bayesian IC" , logger = logger ) :
        pdf.load_params ( rgood , silent = True ) 
        tgood, pgood = bic_good.pvalue ( pdf , data_good )
            
    ## presumably bad fit 
    with timing ( "Bad  fit Bayesian IC" , logger = logger ) : 
        pdf.load_params ( rbad  , silent = True ) 
        tbad, pbad  = bic_bad.pvalue ( pdf , data_bad )
        
    with use_canvas ( title = "Good fit Bayesian IC"  , wait = 2 ) : bic_good.draw ()     
    with use_canvas ( title = "Bad  fit Bayesian IC"  , wait = 2 ) : bic_bad .draw ()
    
    gp = pgood * 100 
    bp = pbad  * 100
    
    gt , ge = pretty_float ( tgood )
    bt , be = pretty_float ( tbad  )
    
    pvg    = clip_pvalue  ( pgood )
    pvb    = clip_pvalue  ( pbad  )        
    nsg    = significance ( pvg   )
    nsb    = significance ( pvb   )
    
    nsg    = '%.2f %s %.2f' % ( nsg.value() , plus_minus , nsg.error () )
    nsb    = '%.2f %s %.2f' % ( nsb.value() , plus_minus , nsb.error () )
    
    row = '%5.2f %s %.2f' % ( gp.value() , plus_minus , gp.error () ) , \
        '%5.2f %s %.2f'   % ( bp.value() , plus_minus , bp.error () ) , nsg , nsb 
    rows.append ( row )
            
    title= 'Goodness-of-Fit Bayesian IC test'
    rows  = T.remove_empty_columns ( rows )     
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )


# ===============================================================================
def probe_GOF ( gof , tag ) :

    logger = getLogger ("test_GOF/%s" % tag )

    ## presumably good fit
    with timing ( "Good fit GOF/%s" % tag , logger = logger ) :
        pdf.load_params ( rgood , silent = True ) 
        tgood        = gof        ( pdf , data_good )
        tgood, pgood = gof.pvalue ( pdf , data_good )
        
    ## presumably bad fit 
    with timing ( "Bad  fit GOF/%s" % tag  , logger = logger ) : 
        pdf.load_params ( rbad  , silent = True ) 
        tbad        = gof        ( pdf , data_bad )
        tbad, pbad  = gof.pvalue ( pdf , data_bad )
        
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

    # ===========================================================================
    ## ADVAL-based tests: 
    config  = { 'nToys' : 100 , 'parallel' : False , 'progress' : True }

    # ===========================================================================
    try : # =====================================================================
        # =======================================================================
        import lightgbm
        from   ostap.stats.gofnd import ADVAL_LightGBM as GOF 
        ## 
        gof    = GOF ( **config )
        test   = gof , typename ( gof ) 
        to_test.append ( test )
        ## 
        # =======================================================================
    except ImportError : # ======================================================
        # =======================================================================
        logger.error  ( "LightGBM is not available, skip the test" )

    """ 
    # ===========================================================================
    try : # =====================================================================
        # =======================================================================
        import sklearn
        from   ostap.stats.gofnd import ADVAL_HistoGBoost as GOF 
        ## 
        gof    = GOF ( **config )
        test   = gof , typename ( gof ) 
        to_test.append ( test ) 
        ## 
        ## 
        # =======================================================================
    except ImportError : # ======================================================
        # =======================================================================
        logger.error  ( "HistoGBoost is not available, skip the test" )

    # ===========================================================================
    try : # =====================================================================
        # =======================================================================
        import sklearn
        from   ostap.stats.gofnd import ADVAL_GBoost as GOF 
        ## 
        gof    = GOF (  **config   )
        test   = gof , typename ( gof ) 
        to_test.append ( test ) 
        ## 
        # =======================================================================
    except ImportError : # ======================================================
        # =======================================================================
        logger.error  ( "GBoost is not available, skip the test" )
        
    # ===========================================================================
    try : # =====================================================================
        # =======================================================================
        import xgboost
        from   ostap.stats.gofnd import ADVAL_XGBoost as GOF 
        ## 
        gof    = GOF ( **config )
        test   = gof , typename ( gof ) 
        to_test.append ( test ) 
        ## 
        # =======================================================================
    except ImportError : # ======================================================
        # =======================================================================
        logger.error  ( "XGBoost is not available, skip the test" )

    """
    
    # ===========================================================================
    from   ostap.core.cpu_info      import HAS_AVX2
    # ===========================================================================
    if HAS_AVX2 :
        # =======================================================================
        try : # =================================================================
            # ===================================================================
            import catboost
            from   ostap.stats.gofnd import ADVAL_CatBoost as GOF 
            ## 
            gof    = GOF ( **config )
            test   = gof , typename ( gof ) 
            to_test.append ( test ) 
            ## 
            # ===================================================================
        except ImportError : # ==================================================
            # ===================================================================
            logger.error  ( "CatBoost is not available, skip the test" )

    ##

    rows  = [ ( '' , 'p-value/good[%]' , 'p-value/bad[%]' , '#%s/good' % greek_lower_sigma , '#%s/bad' % greek_lower_sigma , 'time' ) ] 

    import pickle
    ## run tests: 
    for gof, tag  in to_test :
        print  ( 'CHECK PICKLE for ' , tag )
        pickle.loads ( pickle.dumps ( gof ) )
        
    ## run tests: 
    for gof, tag  in to_test :
        with timing ( 'Processing GoF/%s' % tag , logger = logger ) as timer :
            row = probe_GOF ( gof , tag  )
        row = list ( row ) + [ '%.1f' % timer.delta ] 
        rows.append ( row )
            
    title = 'Goodness-of-Fit tests'
    rows  = T.remove_empty_columns ( rows )     
    table = T.table ( rows , title = title , prefix = '# ')
    logger.info ( '%s:\n%s' % ( title , table ) )
    
    
    
# ===============================================================================
if '__main__' == __name__ :

    ## test_PPD    ()
    ## test_DNN    ()
    ## test_USTAT  ()
    test_GOF  ()    
    ## test_NLL    ()
    ## test_AIC    ()
    ## test_BIC    ()


# ===============================================================================
##                                                                        The END 
# ===============================================================================
                   
    
