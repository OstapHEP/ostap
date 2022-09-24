#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_resolutions.py
# Test module for ostap/fitting/resolution.py
# - It tests various resoltuion models
# ============================================================================= 
"""Test module for ostap/fitting/resolution.py
- It tests various resoltuion models
"""
# ============================================================================= 
from   __future__               import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ostap.fitting.roofit 
from   ostap.core.core          import VE, dsID
from   builtins                 import range
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
from   ostap.logger.colorized   import attention 
import ostap.logger.table       as     T
from   ostap.core.meta_info     import root_info, python_info
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_resolutions' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make simple test mass 
mass     = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , -3 , 3 )

from   ostap.fitting.resolution import ResoHypatia
reso = ResoHypatia       ( 'H'    , mass   ,
                           sigma  =   0.5  ,
                           zeta   =   2.0  , 
                           kappa  =   0.0  ,
                           lambd  =  -1.0  ,
                           sigma0 =   0.01 )

dataset = reso.generate ( 10000 )
logger.info ('DATASET\n%s' % dataset.table ( prefix = "# " ) ) 

stats = {} 

# =============================================================================
def make_print ( pdf , fitresult , title , logger = logger ) :
    
    title2 ='%s' % title 
    logger.info ('%-20s : fit result\n%s' % ( title ,
                                              fitresult.table ( title  = title2 ,
                                                                prefix = '# '   ) ) )
    
    rows = [ ( 'Parameter' , 'Value', '(Roo)Value' ) ]

    mean, mode, median, midpoint, rms , fwhm, skewness , kurtosis = 8 * ( attention ( '<error>' ) , )
    
    model = pdf
    
    try    :
        mean     =  "%+.3g" % model.get_mean     () 
    except :
        pass
    
    try    :
        mode     =  "%+.3g" % model .mode       () 
    except :
        pass 

    try    :
        median   =  "%+.3g" % model .median     () 
    except :
        pass
    
    try    :
        midpoint =  "%+.3g" % model .mid_point  () 
    except :
        pass

    try    :
        rms      =  "%+.3g" % model .rms        () 
    except :
        pass
    
    try    :
        fwhm     =  "%+.3g" % model .fwhm       () 
    except :
        pass
    
    try    :
        skewness =  "%+.3g" % model .skewness   () 
    except :
        pass

    try    :
        kurtosis =  "%+.3g" % model .kurtosis   () 
    except :
        pass
    
    row   = model.name  , mean , mode , median , midpoint , rms , fwhm , skewness , kurtosis

    stats [ model.name ] = row

    try : 
        row = 'mean'       , mean , '%+.6g' % pdf.roo_mean  ()
        rows.append ( row )
    except :
        pass 
    
    row = 'mode'       , mode          , ''
    rows.append ( row )
    row = 'median'     , median        , ''
    rows.append ( row )
    row = 'midpoint'   , midpoint      , ''
    rows.append ( row )

    try : 
        row = 'rms'        , rms           , '%+.6g' % pdf.roo_rms  ()
        rows.append ( row )
    except :
        pass
    
    row = 'FWHM'       , fwhm          , ''
    rows.append ( row )

    try : 
        row = 'skewness'   , skewness      , '%+.6g' % pdf.roo_skewness  ()
        rows.append ( row )
    except :
        pass

    try : 
        row = 'kurtosis'   , kurtosis      , '%+.6g' % pdf.roo_kurtosis  ()  
        rows.append ( row )
    except :
        pass 
        
    table = T.table ( rows , title = title2 ,  prefix = '# ' )
    logger.info ( 'Global features for %s\n%s' % ( title2 , table ) ) 


models = set() 
plots  = set() 

# =============================================================================
## Single gauss
# =============================================================================
def test_gauss () :
## if 1 < 2 :
    
    logger = getLogger ( 'test_gauss' )
    
    logger.info ('Test ResoGauss:  single Gaussian resolution' )
    from   ostap.fitting.resolution import ResoGauss 
    reso = ResoGauss( 'Gauss' , mass ,  sigma = ( 0.2 , 0.01 , 5 )  )
    
    result, frame = reso. fitTo ( dataset , silent = True )
    result, frame = reso. fitTo ( dataset , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_gauss' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :
        make_print ( reso , result , 'Gaussian', logger )

    models.add ( reso  )
    plots .add ( frame )
    
# =============================================================================
## Double gauss
# =============================================================================
def test_2gauss () :
    
    logger = getLogger ( 'test_2gauss' )

    logger.info ('Test ResoGauss2:  double Gaussian resolution' )
    from   ostap.fitting.resolution import ResoGauss2
    reso = ResoGauss2( 'Gauss2' , mass ,
                       sigma    = ( 0.2 , 0.01 , 5     ) ,
                       fraction = ( 0.5 , 0.01 , 0.95  ) , 
                       scale    = ( 1.2 , 1.01 , 5.0   ) )
    
    result, frame = reso. fitTo ( dataset , silent = True )
    result, frame = reso. fitTo ( dataset , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_2gauss' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Double Gaussian', logger )
        
    models.add ( reso  )
    plots .add ( frame ) 

        
# =============================================================================
## Symmetric Apollonios
# =============================================================================
def test_apo2 () :
    
    logger = getLogger ( 'test_apo2' )

    logger.info ('Test ResoApo2:  symmetric Apollonios2 resolution' )
    from   ostap.fitting.resolution import ResoApo2
    reso= ResoApo2( 'Apollonious2' , mass ,
                    sigma = ( 0.4 ,  0.01 ,  5.0 ) ,
                    beta  = ( 0.4 ,  0.01 , 20.0 ) ) 
    
    result, frame = reso. fitTo ( dataset , silent = True )
    result, frame = reso. fitTo ( dataset , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_apo2' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Symmetric Apollonious/2', logger )
    
    models.add ( reso  )
    plots .add ( frame ) 
        
# =============================================================================
## Symmetric double-sided Crystal Ball 
# =============================================================================
def test_cb2 () :
    
    logger = getLogger ( 'test_cb2' )

    logger.info ('Test ResoCB2: symmetric double-sided Crystal Ball resolution' )
    from   ostap.fitting.resolution import ResoCB2
    reso = ResoCB2( 'CrystalBall2' , mass ,                        
                    sigma = ( 0.4 ,  0.01 ,  5.0 ) ,
                    n     = ( 1.0 ,  0.0  , 20.0 ) ,
                    alpha = ( 1.0 ,  0.5  ,  5.0 ) )
    
    result, frame = reso. fitTo ( dataset , silent = True  )
    result, frame = reso. fitTo ( dataset , silent = True  )    
    with wait ( 1 ) , use_canvas ( 'test_cb2' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Symmetric Crystal Ball/2', logger )
 
    models.add ( reso  )
    plots .add ( frame ) 

# =============================================================================
## Hyperbolic secant 
# =============================================================================
def test_sech () :
    
    logger = getLogger ( 'test_sech' )

    logger.info ('Test ResoSech: hyperbolic secant resolution' )
    from   ostap.fitting.resolution import ResoSech
    reso = ResoSech ( 'Sech' , mass , sigma = ( 0.5 , 0.01 , 5.0 ) )
    
    result, frame = reso. fitTo ( dataset , silent = True  )
    result, frame = reso. fitTo ( dataset , silent = True  )    
    with wait ( 1 ) , use_canvas ( 'test_sech' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Hyperbolic secant', logger )
 
    models.add ( reso  )
    plots .add ( frame ) 


# =============================================================================
## Logistic 
# =============================================================================
def test_logistic () :
    
    logger = getLogger ( 'test_logistic' )

    logger.info ('Test ResoLogistic: logistic (sech-squared) resolution' )
    from   ostap.fitting.resolution import ResoLogistic
    reso = ResoLogistic ( 'Logistic' , mass ,  sigma = ( 0.4 , 0.01 , 5.0 ) )
            
    result, frame = reso. fitTo ( dataset , silent = True  )
    result, frame = reso. fitTo ( dataset , silent = True  )    
    with wait ( 1 ) , use_canvas ( 'test_logistic' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Logistic', logger )
 
    models.add ( reso  )
    plots .add ( frame ) 

# =============================================================================
## symmetric Bukin
# =============================================================================
def test_bukin () :
    
    logger = getLogger ( 'test_bukin' )

    logger.info ('Test ResoBukin: symmetric Bukin resolution' )
    from   ostap.fitting.resolution import ResoBukin

    reso = ResoBukin ( 'Bukin' , mass ,
                       sigma = ( 0.1 , 0.01 , 5.0 ) ,
                       rho   = (   0 , 0 , 10     ) )

    result, frame = reso. fitTo ( dataset , silent = True  )
    result, frame = reso. fitTo ( dataset , silent = True  )    
    with wait ( 1 ) , use_canvas ( 'test_bukin' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Symmetric Bukin', logger )
 
    models.add ( reso  )
    plots .add ( frame ) 

# =============================================================================
## symmetric Johnson's SU 
# =============================================================================
def test_johnsonSU () :
    
    logger.info ('Test JohnsonSU: symmetric JohnsonSU  resolution' )
    from   ostap.fitting.resolution import ResoJohnsonSU 
    reso = ResoJohnsonSU ( 'JohnsonSU' , mass ,
                           delta = ( 1.7 , 1.e-6 , 1000 ) ,
                           lambd = ( 0.2 , 1.e-6 , 1000 ) ) 
    
    result, frame = reso. fitTo ( dataset , silent = True  )
    result, frame = reso. fitTo ( dataset , silent = True  )    
    with wait ( 1 ) , use_canvas ( 'test_johnsonSU' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Symmetric Jonhson-SU', logger )
 
    models.add ( reso  )
    plots .add ( frame ) 

# =============================================================================
## Sinh-Asinh
# =============================================================================
def test_sinhasinh () :
## if 1 < 2 :
    
    logger = getLogger ( 'test_sinhasinh' )

    logger.info ('Test SinhAsinh: symmetric SinhAsinh resolution' )
    from   ostap.fitting.resolution import ResoSinhAsinh
    reso = ResoSinhAsinh ( 'SinhAsinh' , mass ,  
                           delta = ( 0.7 , 1.e-5 , 1000 ) )
    
    reso.delta.fix     ( 0.73 ) 
    reso.sigma.release () 
    reso.sigma = 0.31
    
    result, frame = reso. fitTo ( dataset , silent = True  )
    reso.delta.release () 
    result, frame = reso. fitTo ( dataset , silent = True  )    
    with wait ( 1 ) , use_canvas ( 'test_sinhasinh' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Symmetric SinhAsinh', logger )
 
    models.add ( reso  )
    plots .add ( frame ) 

# =============================================================================
## Hyperbolic 
# =============================================================================
def test_hyperbolic () :
    
    logger = getLogger ( 'test_hyperbolic' )

    logger.info ('Test Hyperbolic: symmetric Hyperbolic resolution' )
    from   ostap.fitting.resolution import ResoHyperbolic
    reso = ResoHyperbolic ( 'Hyperbolic' , mass ,
                            sigma = ( 0.1 , 0.01 , 5.0 ) ,                            
                            zeta = ( 1.0 , 1.e-5 , 1.e+5 ) )

    result, frame = reso. fitTo ( dataset , silent = True  )
    result, frame = reso. fitTo ( dataset , silent = True  )    
    with wait ( 1 ) , use_canvas ( 'test_hyperbolic' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Symmetric Hyperbolic', logger )
 
    models.add ( reso  )
    plots .add ( frame ) 

# =============================================================================
## Generalized Hyperbolic 
# =============================================================================
def test_genhyperbolic () :
    
    logger = getLogger ( 'test_genhyperbolic' )

    logger.info ('Test Hyperbolic: symmetric generalised Hyperbolic resolution' )
    from   ostap.fitting.resolution import ResoGenHyperbolic
    reso = ResoGenHyperbolic ( 'GenHyperbolic' , mass ,
                               sigma = ( 0.1 , 0.01 , 5.0 ) ,
                               zeta = ( 1.0 , 1.e-5 , 1.e+5 ) ,
                               lambd =  (-20,20) )

    result, frame = reso. fitTo ( dataset , silent = True  )
    result, frame = reso. fitTo ( dataset , silent = True  )    
    with wait ( 1 ) , use_canvas ( 'test_genhyperbolic' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Symmetric GenHyperbolic', logger )
 
    models.add ( reso  )
    plots .add ( frame ) 

# =============================================================================
## Hypatia
# =============================================================================
def test_hypatia () :
    
    logger = getLogger ( 'test_hypatia' )

    logger.info ('Test Hyperbolic: symmetric generalised Hyperbolic resolution' )
    from   ostap.fitting.resolution import ResoHypatia
    reso = ResoHypatia ( 'Hypatia' , mass ,
                         sigma  = ( 0.1 , 0.01 , 5.0 ) ,
                         zeta   = ( 1.0 , 1.e-5 , 1.e+5 ) ,
                         lambd  = ( -20,20 ) ,
                         sigma0 = 0.01 )
    
    for i in range ( 6 ) :
        result, frame = reso. fitTo ( dataset , silent = True  )
        
    with wait ( 1 ) , use_canvas ( 'test_genhyperbolic' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Symmetric Hypatia', logger )
 
    models.add ( reso  )
    plots .add ( frame ) 

# =============================================================================
## Generalized Gauss v1 
# =============================================================================
def test_gengaussv1  () : 
## if 1 <2  :
    
    logger = getLogger ( 'test_generalized_gauss_v1' )
    
    logger.info ('Test GenGaussV1 : symmetric generalised gaussian resolution' )
    from   ostap.fitting.resolution import ResoGenGaussV1
    reso = ResoGenGaussV1 ( 'GenGaussV1' , mass ,
                            alpha  = ( 0.5 , 0.001 , 10.0  ) ,
                            beta   = ( 1.4 , 0.1   , 5     ) )

    reso.alpha.release () 
    reso.beta .fix     () 
    result, frame = reso. fitTo ( dataset , silent = True  )
    result, frame = reso. fitTo ( dataset , silent = True  )
    reso.beta.release  () 
    for i in range ( 4 ) :
        result, frame = reso. fitTo ( dataset , silent = True  )
        
    with wait ( 1 ) , use_canvas ( 'test_gengaussv1' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Symmetric GenGaussv1', logger )
 
    models.add ( reso  )
    plots .add ( frame ) 

# =============================================================================
## Das 
# =============================================================================
def test_das () :
    
    logger = getLogger ( 'test_das' )

    logger.info ('Test Das: Gaussian with symmetric exponential tails ' )
    from   ostap.fitting.resolution import ResoDas
    reso = ResoDas ( 'Das' , mass ,
                     k     = ( 1.33 , 1.e-5 , 20 ) , 
                     sigma = ( 0.41 , 0.01 , 5.0 ) )
    
    result, frame = reso. fitTo ( dataset , silent = True  )
    result, frame = reso. fitTo ( dataset , silent = True  )    
    with wait ( 1 ) , use_canvas ( 'test_das' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Symmetric Das', logger )
 
    models.add ( reso  )
    plots .add ( frame ) 

# =============================================================================
## Normal Laplace 
# =============================================================================
def test_normlapl () :
    
    logger = getLogger ( 'test_normlapl' )

    logger.info ('Test Normal Laplace : Gaussian with symmetric exponential tails' )
    from   ostap.fitting.resolution import ResoNormalLaplace
    reso = ResoNormalLaplace ( 'NL' , mass ,
                               k  = ( 1.0 , 1.e-5 , 20 ) , 
                               varsigma = ( 0.1 , 0.01 , 5.0 ) )
    
    result, frame = reso. fitTo ( dataset , silent = True  )
    result, frame = reso. fitTo ( dataset , silent = True  )    
    with wait ( 1 ) , use_canvas ( 'test_normlapl' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Symmetric Normal Laplace', logger )
 
    models.add ( reso  )
    plots .add ( frame ) 

# ==============================================================================
## dump all models
# ==============================================================================
def dump_models () :
    
    header =  'Model'   , \
             'mean'     , 'mode' , 'midpoint' , 'median' , \
             'rms'      , 'fwhm' ,  \
             'skewness' , 'kurtosis'
    
    rows = [ header ] 
    for m in sorted ( stats  ) :
        rows.append ( stats [ m ]  )

    table = T.table ( rows , title = "Model's features" ,  prefix = '# ' )
    logger.info ( 'Features of models\n%s' % table )
  

# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    
    logger = getLogger ( 'test_db' )
    
    if root_info < (6,20) and (3,0)<= python_info < (3,7) :
        logger.warning ( "Test is disabled for this version of ROOT+python" )
        return 
    
    logger.info('Saving all objects into DBASE')
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db : 
        db['mass'] = mass
        db['dataset'  ] = dataset
        for m in models :
            db['model:' + m.name ] = m
            db['roo_tot:%s' % m.name ] = m.pdf
            for i,s in enumerate ( m.signals ) :
                db['roo_sig%d:%s' % ( i , m.name ) ] = s
            for i, b in enumerate ( m.backgrounds ) : 
                db['roo_bkg%d:%s' % ( i , m.name ) ] = s
            for a in m.alist1 : 
                db['cmp:%s_%s' % ( m.name , a.name ) ] = a
        db[ 'models' ] = models
        for i , p in enumerate ( plots ) :
            db [ ' plot:%s_%s' % ( i , p.name ) ] = p            
        db[ 'plots'  ] = plots 
        db.ls() 
        
# =============================================================================
if '__main__' == __name__ :
    
    with timing ("Gauss"     , logger ) :  
        test_gauss          () ## single Gaussian resolution model
        
    with timing ("2-Gauss"   , logger ) :  
        test_2gauss         () ## double Gaussian resolution model
        
    with timing ("Apo2"      , logger ) :  
        test_apo2           () ## symmetric Apollonios resoltuion model
        
    with timing ("CB2"       , logger ) :  
        test_cb2            () ## double-sided Crystal Ball resoltuion model
        
    with timing ("Sech"      , logger ) :  
        test_sech           () ## hyperbolic secant resolution model

    with timing ("Logistic"  , logger ) :  
        test_logistic       () ## logistic resolution model
        
    with timing ("Bukin"     , logger ) :  
        test_bukin          () ## Bukin resolution model
    
    with timing ("SinhAsinh" , logger ) :  
        test_sinhasinh      () ## SinhAsinh resolution model

    with timing ("JohnsonSU" , logger ) :  
        test_johnsonSU      () ## JohnsonSU resolution model

    with timing ("Hyperbolic" , logger ) :  
        test_hyperbolic     () ## Hyperbolic resolution model
    
    with timing ("GenHyperbolic" , logger ) :  
        test_genhyperbolic  () ## generalised Hyperbolic resolution model

    with timing ("Hypatia" , logger ) :  
        test_hypatia        () ## generalised Hyperbolic resolution model

    with timing ("GenGaussV1" , logger ) :  
        test_gengaussv1     ()   ## Das resolution model
    
    with timing ("Das"        , logger ) :  
        test_das           ()   ## Das resolution model

    with timing ("NormalLaplace" , logger ) :  
        test_normlapl      ()   ## Normal Laplace resolution model

    ## check finally that everything is serializeable:
    with timing ("Save to DB"    , logger ) :  
        test_db ()          

    ## print all models for comparion
    with timing ("Dump models"    , logger ) :  
        dump_models ()           

        

# =============================================================================
##                                                                      The END 
# ============================================================================= 
