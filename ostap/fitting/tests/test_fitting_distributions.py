#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/fitting/tests/test_fitting_distributions.py
# Test module for ostap/fitting/distributions.py
# ============================================================================= 
""" Test module for ostap/fitting/distributions.py
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.root_utils   import batch_env 
import ostap.fitting.models     as     Models 
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_distributions' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

x = ROOT.RooRealVar( 'x' , 'x-variable' , 0  , 10 ) 
y = ROOT.RooRealVar( 'x' , 'y-variable' , -5 ,  5 ) 

plots  = set ()
models = set ()

# =============================================================================
def test_GammaDist () :

    logger = getLogger("test_GammaDist")

    model  = Models.GammaDist_pdf ( 'GammaDist' , x ,
                                    k     = ( 9     , 1     , 100 ) ,
                                    theta = ( 0.333 , 0.001 , 2   ) )
    
    with use_canvas ( 'GammaDist_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_GenGammaDist () :

    logger = getLogger("test_GenGammaDist")

    model  = Models.GenGammaDist_pdf ( 'GenGammaDist' , x ,
                                       k     = ( 9     , 1     , 100 ) ,
                                       theta = ( 0.333 , 0.001 , 2   ) )
    
    with use_canvas ( 'GenGammaDist_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_Amoroso () :

    logger = getLogger("test_Amoroso")

    model  = Models.Amoroso_pdf ( 'Amoroso' , x ,
                                  alpha = ( 9     , 1     , 100 ) ,
                                  theta = ( 0.333 , 0.001 , 2   ) )
    
    with use_canvas ( 'Amoroso_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_LogGammaDist () :

    logger = getLogger("test_LogGammaDist")
    
    model  = Models.LogGammaDist_pdf ( 'LogGammaDist' , y ,
                                       k     = ( 9     , 1     , 100 ) ,
                                       theta = ( 0.333 , 0.001 , 2   ) )
    
    with use_canvas ( 'LogGammaDist_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_Log10GammaDist () :

    logger = getLogger("test_Log10GammaDist")
    
    model  = Models.Log10GammaDist_pdf ( 'Log10GammaDist' , y ,
                                       k     = ( 9     , 1     , 100 ) ,
                                       theta = ( 0.333 , 0.001 , 2   ) )
    
    with use_canvas ( 'Log10GammaDist_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_LogGamma () :

    logger = getLogger("test_LogGamma")
    
    model  = Models.LogGamma_pdf ( 'LogGamma' , x , nu = 0 , lambd = 1 , alpha = 1  )
    
    with use_canvas ( 'LogGamma_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_BetaPrime () :

    logger = getLogger("test_BetaPrime")
    
    model  = Models.BetaPrime_pdf ( 'BetaPrime' , x , alpha = 1 , beta = 1 , scale = 1 , delta = 0 )
    
    with use_canvas ( 'BetaPrime_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_Landau () :

    logger = getLogger("test_Landau")
    
    model  = Models.Landau_pdf ( 'Landau' , x , scale = 1 , delta = 0 )
    
    with use_canvas ( 'BetaPrime_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 


# =============================================================================
def test_Argus () :

    logger = getLogger("test_Argus")
    
    model  = Models.Argus_pdf ( 'Argus' , x , c = 7 , mu = 8 , chi = 0.1  )
    
    with use_canvas ( 'Argus_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_GenArgus () :

    logger = getLogger("test_GenArgus")
    
    model  = Models.GenArgus_pdf ( 'GenArgus' , x , c = 7 , mu = 8 , chi = 0.1 , dp = 2  )
    
    with use_canvas ( 'GenArgus_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_TwoExpos () :

    logger = getLogger("test_TwoExpos")
    
    model  = Models.TwoExpos_pdf ( 'TwoExpos' , x , alpha = 1 , delta = 1 )
    
    with use_canvas ( 'TwoExpos_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 


# =============================================================================
def test_Gumbel () :

    logger = getLogger("test_Gumbel")
    
    model  = Models.Gumbel_pdf ( 'Gumbel' , x , mu = 0 , beta = 1  )
    
    with use_canvas ( 'Gumbel_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_Rice () :

    logger = getLogger("test_Rice")
    
    model  = Models.Rice_pdf ( 'Rice' , x , nu = 0 , varsigma = 1  )
    
    with use_canvas ( 'Rice_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_GenInvGauss () :

    logger = getLogger("test_GenInvGauss")
    
    model  = Models.GenInvGauss_pdf ( 'GenInvGauss' , x ,
                                      theta = 1 , eta = 1 , p = 0 )
    
    with use_canvas ( 'GenInvGauss_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_Weibull () :

    logger = getLogger("test_Weibull")
    
    model  = Models.Weibull_pdf ( 'Weibull' , x ,
                                  scale = 1 , shape = 1 )
    
    with use_canvas ( 'Weibull_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_Tsallis () :

    logger = getLogger("test_Tsallis")
    
    model  = Models.Tsallis_pdf ( 'Tsallis' , x ,
                                  m0 = 1   ,
                                  n  = 10  ,
                                  T  = 0.2 )
    
    with use_canvas ( 'Tsallis_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_QGSM () :

    logger = getLogger("test_QGSM")
    
    model  = Models.QGSM_pdf ( 'QGSM' , x ,
                               m0 = 1   ,
                               b  = 10  )
    
    with use_canvas ( 'QGSM_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_Hagedorn () :

    logger = getLogger("test_Hagedorn")
    
    model  = Models.Hagedorn_pdf ( 'Hagedorn' , x ,
                                   m0   = 1   ,
                                   beta = 10  )
    
    with use_canvas ( 'QGSM_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 


# =============================================================================
def test_GenPareto () :

    logger = getLogger("test_GenPareto")
    
    model  = Models.GenPareto_pdf ( 'GenPareto' , x ,
                                    mu    = 1 ,
                                    scale = 1 ,
                                    shape = 1 )
    
    with use_canvas ( 'GenPareto_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_ExGenPareto () :

    logger = getLogger("test_ExGenPareto")
    
    model  = Models.ExGenPareto_pdf ( 'ExGenPareto' , x ,
                                      mu    = 1 ,
                                      scale = 1 ,
                                      shape = 1 )
    
    with use_canvas ( 'ExGenPareto_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_GEV () :

    logger = getLogger("test_GEV")
    
    model  = Models.GEV_pdf ( 'GEV' , x ,
                              mu    = 1 ,
                              scale = 1 ,
                                      shape = 1 )
    
    with use_canvas ( 'GEV_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 

# =============================================================================
def test_MPERT () :

    logger = getLogger("test_MPERT")
    
    model  = Models.MPERT_pdf ( 'MPERT' , x ,
                                xi    = 3  ,
                                gamma = 10 )
    
    with use_canvas ( 'MPERT_pdf' ) :
        plot = model.draw()
        
    models.add ( model )
    plots .add ( plot  ) 


# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        for m in models :
            db['model:' + m.name ] = m
            db['roo:%s' % m.name ] = m.pdf
        db['models'   ] = models
        db['plots'     ] = plots 
        db.ls() 

# =============================================================================
if '__main__' == __name__ :
    
    test_GammaDist      ()
    test_GenGammaDist   ()
    test_Amoroso        ()
    test_LogGammaDist   ()
    test_Log10GammaDist ()
    test_LogGamma       ()
    test_BetaPrime      ()
    test_Landau         ()
    test_Argus          ()
    test_GenArgus       ()
    test_TwoExpos       ()
    test_Gumbel         ()
    test_Rice           ()
    test_GenInvGauss    ()
    test_Weibull        ()
    test_Tsallis        ()
    test_QGSM           ()
    test_Hagedorn       ()
    test_GenPareto      ()
    test_ExGenPareto    ()
    test_GEV            ()
    test_MPERT          ()

    ## check that everything is serializeable 
    test_db ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
