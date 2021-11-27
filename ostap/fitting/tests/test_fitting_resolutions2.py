#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_resolutions.py
# Test module for ostap/fitting/resolution2.py
# - It tests various asymmetric resolution models
# ============================================================================= 
"""Test module for ostap/fitting/resolution2.py
- It tests various asymmetric resolution models
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ROOT, random
ROOT.PyConfig.IgnoreCommandLineOptions = False 
# 
import ostap.fitting.roofit 
from   ostap.core.core          import VE, dsID
from   builtins                 import range
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait
import ostap.logger.table       as     T 
from   ostap.core.meta_info     import root_info, python_info
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'ostap.test_fitting_resolutions2' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make simple test mass 
mass     = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , -5 , 5 )

logger.info ( 'BATCH? %s' % ROOT.gROOT.IsBatch() ) 


from   ostap.fitting.resolution import ResoHypatia
reso = ResoHypatia       ( 'H'    , mass  ,
                           sigma  =   0.5  ,
                           zeta   =   2.0  , 
                           kappa  =   0.1  ,
                           lambd  =  -1.0  ,
                           sigma0 =   0.01 ) 

##from   ostap.fitting.resolution import ResoGenHyperbolic
##reso = ResoGenHyperbolic ( 'GH'   , mass  ,
##                           sigma  =  0.5  ,
##                           zeta   =  2.0  , 
##                           kappa  =  0.1  ,
##                           lambd  = -1.0  )

## from   ostap.fitting.resolution import ResoHyperbolic
## reso = ResoHyperbolic    ( 'HH'   , mass  ,
##                           sigma  =  0.5  ,
##                           zeta   =  10   , 
##                           kappa  =  0.1  )

dataset = reso.generate ( 10000 )

logger.info ('DATASET\n%s' % dataset.table ( prefix = "# " ) ) 

models = set() 


# =============================================================================
def make_print ( pdf , fitresult , title , logger = logger ) :

    title2 ='%s resolution model' % title 
    logger.info ('%-20s : fit result\n%s' % ( title ,
                                              fitresult.table ( title  = title2 ,
                                                                prefix = '# '   ) ) )
    
    rows = [ ( 'Parameter' , 'Value') ]
    
    rows.append ( ( 'mean'       , '%+.6g' % pdf.get_mean  () ) )
    rows.append ( ( 'mode'       , '%+.6g' % pdf.mode      () ) )
    rows.append ( ( 'median'     , '%+.6g' % pdf.median    () ) )
    rows.append ( ( 'mindpoint'  , '%+.6g' % pdf.mid_point () ) )
    rows.append ( ( 'rms'        , '%+.6g' % pdf.rms       () ) )
    rows.append ( ( 'FWHM'       , '%+.6g' % pdf.fwhm      () ) )
    rows.append ( ( 'skewness'   , '%+.6g' % pdf.skewness  () ) )
    rows.append ( ( 'kurtosis'   , '%+.6g' % pdf.kurtosis  () ) )

    table = T.table ( rows , title = title2 ,  prefix = '# ' )
    logger.info ( 'Global features for %s\n%s' % ( title2 , table ) ) 
    

# =============================================================================
## Single gauss
# =============================================================================
def test_gauss () :

    logger = getLogger ( 'test_gauss' )
    
    logger.info ('Test ResoGauss: bifurcated Gaussian resolution model' )
    from   ostap.fitting.resolution import ResoGauss 
    reso = ResoGauss( 'Gauss' , mass ,
                      mean  = ( 0.0 , -1   , +1  ) , 
                      sigma = ( 0.5 ,  0.1 , 1.0 ) ,
                      kappa = ( 0.0 , -1.0 , 1.0 ) )
    
    result, frame = reso. fitTo ( dataset , silent = True )
    result, frame = reso. fitTo ( dataset , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_gauss' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )

    
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Asymmetric Gaussian', logger )

    models.add ( reso  )
    
# =============================================================================
## Symmetric Apollonios
# =============================================================================
def test_apo2 () :
    
    logger = getLogger ( 'test_apo2' )

    logger.info ('Test ResoApo2:  asymmetric Apollonios resolution model' )
    from   ostap.fitting.resolution import ResoApo2
    reso = ResoApo2( 'Apollonios2' , mass ,
                     mean  = ( 0.0 , -0.1  , 0.1 ) , 
                     sigma = ( 0.4 ,  0.1  , 1.0 ) ,
                     beta  = ( 1   , 1.e-5 , 100 ) ,
                     kappa = ( 0.1 , -1    , +1  ) )
    
    result, frame = reso. fitTo ( dataset , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_apo2' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Asymmetric Apollonios2', logger )
        
    models.add ( reso )

# =============================================================================
## Asymmetric double-sided Crystal Ball 
# =============================================================================
def test_cb2 () :
    
    logger = getLogger ( 'test_cb2' )

    logger.info ('Test ResoCB2: asymmetric double-sided Crystal Ball resolution model' )
    from   ostap.fitting.resolution import ResoCB2
    reso = ResoCB2( 'CrystalBall2' , mass ,
                    mean   = ( 0.0 , -0.1  ,  0.1 ) , 
                    sigma  = ( 0.2 ,  0.1  ,  1.0 ) ,
                    alpha  = ( 1.0 ,  0.5  ,  3.0 ) ,
                    n      = ( 1.0 ,  0.0  , 20.0 ) ,
                    kappaN = ( 0   ,  -1   , +1   ) ,
                    kappaA = ( 0   ,  -1   , +1   ) )
    
    result, frame = reso. fitTo ( dataset , silent = True )
    result, frame = reso. fitTo ( dataset , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_cb2' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :
        make_print ( reso , result , 'Asymmetric Crystal Ball', logger )
        
    models.add ( reso )


# =============================================================================
## Asymmetric Bukin
# =============================================================================
def test_bukin () :
    
    logger = getLogger ( 'test_bukin' )

    logger.info ('Test ResoBukin: asymmetric Bukin resolution model' )
    from   ostap.fitting.resolution import ResoBukin
    reso = ResoBukin ( 'Bukin' , mass ,
                       mean   = (  0.0  , -0.1  ,  0.1 ) , 
                       sigma  = (  0.4  ,  0.1  ,  1.0 ) ,
                       rho    = (  0.4  ,  0    ,  10  ) ,
                       xi     = ( -0.04 , -10   , +10  ) ,
                       kappa  = (  0.04 , -1    , +1   ) )
    
    result, frame = reso. fitTo ( dataset , silent = True )
    result, frame = reso. fitTo ( dataset , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_bukin' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Asymmetric Bukin', logger )
        
    models.add ( reso )
    
# =============================================================================
## Asymmetric Johnson's SU 
# =============================================================================
def test_johnsonSU () :
    
    logger.info ('Test JohnsonSU: asymmetric JohnsonSU  resolution model' )
    from   ostap.fitting.resolution import ResoJohnsonSU 
    reso = ResoJohnsonSU ( 'JohnsonSU' , mass ,
                           mean  = ( 0.0 , -0.1  ,  0.1 ) ,                            
                           delta = ( 1.7 , 1.e-6 , 1000 ) ,
                           lambd = ( 0.2 , 1.e-6 , 1000 ) ,
                           gamma = ( 0   , -100  , 100  ) ) 
                           
    result, frame = reso. fitTo ( dataset , silent = True )
    result, frame = reso. fitTo ( dataset , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_johnsonSU' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Asymmetric Johnson-SU', logger )
        
    models.add ( reso )



# =============================================================================
## Sinh-Asinh
# =============================================================================
def test_sinhasinh () :
    
    logger = getLogger ( 'test_sinhasinh' )

    logger.info ('Test SinhAsinh: symmetric SinhAsinh resolution model' )
    from   ostap.fitting.resolution import ResoSinhAsinh
    reso = ResoSinhAsinh ( 'SinhAsinh' , mass ,
                           mean    = ( 0.0 , -0.1  ,  0.1 ) ,                            
                           delta   = ( 0.7 , 1.e-5 , 1000 ) ,
                           sigma   = ( 0.2 ,  0.1  ,  1.0 ) ,
                           epsilon = ( 0   , -100  ,  100 ) ) 
                           
    
    result, frame = reso. fitTo ( dataset , silent = True )
    result, frame = reso. fitTo ( dataset , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_sinhasinh' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Asymmetric SinhAsinh', logger )
        
    models.add ( reso )


# =============================================================================
## Asymmetric Hyperbolic 
# =============================================================================
def test_hyperbolic () :
    
    logger = getLogger ( 'test_hyperbolic' )

    logger.info ('Test Hyperbolic: asymmetric Hyperbolic resolution model' )
    from   ostap.fitting.resolution import ResoHyperbolic
    reso = ResoHyperbolic ( 'Hyperbolic' , mass ,
                            mu      = ( 0.1  , -5    ,  5   ) ,                            
                            sigma   = ( 0.5  ,  0.1  ,  1.0 ) ,
                            kappa   = ( -3   , -100  , +100 ) ,                            
                            zeta    = ( +30  , 1.e-5 , 1.e+5 ) )
    
    result, frame = reso. fitTo ( dataset , silent = True )
    result, frame = reso. fitTo ( dataset , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_hyperbolic' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Asymmetric Hyperbolic', logger )
        
    models.add ( reso )


# =============================================================================
## Generalized Hyperbolic 
# =============================================================================
def test_genhyperbolic () :
    
    logger = getLogger ( 'test_genhyperbolic' )

    logger.info ('Test Hyperbolic: asymmetric generalised Hyperbolic resolution model' )
    from   ostap.fitting.resolution import ResoGenHyperbolic
    reso = ResoGenHyperbolic ( 'GenHyperbolic' , mass ,
                               mu    = ( 0.1  , -5    ,  5    ) ,                            
                               sigma = ( 0.5  ,  0.1  ,  1.0  ) ,
                               kappa = ( -0   , -100  , +100  ) ,                            
                               zeta  = ( +30  , 1.e-1 , 1.e+5 ) , 
                               lambd = ( -2  , -100  ,  100   ) )
    
    result, frame = reso. fitTo ( dataset , silent = True )
    result, frame = reso. fitTo ( dataset , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_genhyperbolic' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Asymmetric Generalised Hyperbolic', logger )
        
    models.add ( reso )

# =============================================================================
## Hypatia 
# =============================================================================
def test_hypatia () :

    logger = getLogger ( 'test_hypatia' )

    logger.info ('Test Hypatia: asymmetric Hypatia resolution model' )
    from   ostap.fitting.resolution import ResoHypatia
    reso = ResoHypatia ( 'Hypatia' , mass ,
                         mu     = ( 0.1  , -5    ,  5    ) ,                            
                         sigma  = ( 0.5  ,  0.1  ,  1.0  ) ,
                         kappa  = ( -0   , -100  , +100  ) ,                            
                         zeta   = ( +30  , 1.e-1 , 1.e+5 ) , 
                         lambd  = ( -2  , -100  ,  100   ) ,
                         sigma0 = 0.01                   )
    
    result, frame = reso. fitTo ( dataset , silent = True )
    result, frame = reso. fitTo ( dataset , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_hypatia' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Asymmetric Hypatia', logger )
        
    models.add ( reso )


    

# =============================================================================
## Das 
# =============================================================================
def test_das () :

    
    logger = getLogger ( 'test_das' )

    logger.info ('Test Das: Gaussian with symmetric exponential tails ' )
    from   ostap.fitting.resolution import ResoDas
    reso = ResoDas ( 'Das' , mass ,
                     mean  = ( 0.1 , -3 , 3 ) ,                    
                     k     = ( 1.0 , 1.e-5 , 200 ) ,
                     kappa = ( 0.1 , -1 , 1  ) )  
    
    result, frame = reso. fitTo ( dataset , silent = True )
    result, frame = reso. fitTo ( dataset , silent = True )
    with wait ( 1 ) , use_canvas ( 'test_das' ) : 
        result, frame = reso. fitTo ( dataset , silent = True , draw = True )
        
    if 0 != result.status() or 3 != result.covQual() :
        logger.warning('Fit is not perfect MIGRAD=%d QUAL=%d ' % ( result.status() , result.covQual () ) )
        print(result)
    else :     
        make_print ( reso , result , 'Asymmetric Das', logger )
        
    models.add ( reso )

# ==============================================================================
## dump all models
# ==============================================================================
def dump_models () :


    header =  'Model'   , \
             'mean'     , 'mode' , 'midpoint' , 'median' , \
             'rms'      , 'fwhm' ,  \
             'skewness' , 'kurtosis'

    mods = {}
    for m in models :
        mods[ m.name ] = m

    rows = [ header ] 
    for m in sorted ( mods ) :
        model = mods [ m ] 
        row = m , \
              '%+.4g' % model.get_mean  () , \
              '%+.4g' % model.mode      () , \
              '%+.4g' % model.median    () , \
              '%+.4g' % model.mid_point () , \
              '%+.4g' % model.rms       () , \
              '%+.4g' % model.fwhm      () , \
              '%+.4g' % model.skewness  () , \
              '%+.4g' % model.kurtosis  ()
        rows.append ( row )

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
        for m in models : db['model %s' % m.name ] = m
        db['models'   ] = models
        db.ls() 

  
    
# =============================================================================
if '__main__' == __name__ :

    
    with timing ("Gauss"     , logger ) :  
        test_gauss      () ## single Gaussian resolution model
        
    with timing ("Apo2"      , logger ) :  
        test_apo2       () ## Apollonios resoltuion model
        
    with timing ("CB2"       , logger ) :  
        test_cb2        () ## double-sided Crystal Ball resoltuion model
        
    with timing ("Bukin"     , logger ) :  
        test_bukin      () ## Bukin resolution model
    
    with timing ("JohnsonSU" , logger ) :  
        test_johnsonSU  () ## JohnsonSU resolution model

    with timing ("SinhAsinh" , logger ) :  
        test_sinhasinh  () ## SinhAsinh resolution model

    with timing ("Hyperbolic" , logger ) :  
        test_hyperbolic  () ## Hyperbolic resolution model
    
    with timing ("GenHyperbolic" , logger ) :  
        test_genhyperbolic  () ## generalised Hyperbolic resolution model

    with timing ("Hypatia" , logger ) :  
        test_hypatia        () ## Hypatia resoltuion model

    with timing ("Das"      , logger ) :  
        test_das           ()   ## Das resolution model

    ## check finally that everything is serializeable:
    with timing ("Save to DB"    , logger ) :  
        test_db ()          

    ## check finally that everything is serializeable:
    with timing ("Dump models"    , logger ) :  
        dump_models ()           

# =============================================================================
##                                                                      The END 
# ============================================================================= 
