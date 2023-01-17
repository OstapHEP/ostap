#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/roostats.py
#  Set of useful basic utilities to dela with RooStats 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2023-01-17
# =============================================================================
""" Set of useful basic utilities to dela with RooStats"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2023-01-17"
__all__     = (
    ##
    'interval_PL'   , ## Profile-likelihood confidence interval or upper/lower limit 
    'interval_FC'   , ## Feldman-Cousins    confidence interval or upper/lower limit 
    'interval_BC'   , ## Bayesian           confidence interval or upper/lowee limir
    'interval_MCMC' , ## MCMC               confidence interval or upper/lower limit 
    ##
    'create_MC_WS'  , ## helper function to create (and fill) ModelConfig & RooWorkspace
    ##
    )
# =============================================================================
import ostap.fitting.roofit
from   ostap.core.ostap_types import string_types 
import ROOT
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.roostats' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================

# =============================================================================
## Create (and fill) ModelConfig & RooWorkspace
#  @code
#  pdf     = ...
#  globars = [ ... ] 
#  mc , ws = create_MC_WS ( pdf , [ 'S' , 'B' ] , dataset = dataset , globobs = () ) 
#  @endcode
#  @see RooStats::ModelConfig 
#  @see RooWorkspace  
def create_MC_WS  ( pdf , params , dataset = None , ws = None , globobs = () ) :
    """ Create (and fill) ModelConfig & RooWorkspace
    >>> pdf     = ...
    >>> globobs = [ ... ] 
    >>> mc , ws = create_MC_WS ( pdf , [ 'S' , 'B' ] , globobs = globobs , dataset = dataset ) 
    - see `ROOT.RooStats.ModelConfig` 
    - see `ROOT.RooWorkspace`  
    """
    if   isinstance ( params , ROOT.RooAbsReal ) : params = [ params ]
    elif isinstance ( params , string_types    ) : params = [ params ]
    
    parameters = [ pdf.parameter ( p , dataset ) for p in params ]
    
    if not ws : ws = ROOT.RooWorkspace         ( 'WS%s' % pdf.name  , 'workspace' )
        
    mc = ROOT.RooStats.ModelConfig ( 'MC%s' % pdf.name  , 'model-config for %s' % pdf  , ws )

    mc.SetPdf                  ( pdf.pdf     )
    mc.SetObservables          ( pdf.vars    )

    poi = ROOT.RooArgSet       ( *parameters ) 
    mc.SetParametersOfInterest (  poi        )
    
    nuis = [ v for v in pdf.params ( dataset ) if not v in pdf.vars and not v in poi ]
    nuis = ROOT.RooArgSet      ( *nuis )    
    mc.SetNuisanceParameters   (  nuis )

    pdf.aux_keep.append ( poi  )
    pdf.aux_keep.append ( nuis )

    if globobs :
        gobs = [ pdf.parameter  ( v , dataset ) for v in globobs ]
        gobs = ROOT.RooArgSet   ( *gobs )
        mc.SetGlobalObservables (  gobs )
        pdf.aux_keep.append     (  gobs )
        
    
    return mc, ws


# =============================================================================
## helper method to get the actual innterval/limits
def get_result ( interval , limit , par = None , ws = None , mc = None ) :
    """Helper methdo to get the actual innterval/limits"""
    
    if   limit and 0 < limit :
        high   = interval.UpperLimit ( par ) if par else interval.UpperLimit ()
        result = high
    elif limit and 0 > limit :
        low    = interval.LowerLimit ( par ) if par else interval.LowerLimit ()
        result = low 
    else :        
        low    = interval.LowerLimit ( par ) if par else interval.LowerLimit ()
        high   = interval.UpperLimit ( par ) if par else interval.UpperLimit ()
        result = low , high 
        
    if ws : del ws
    if mc : del mc  

    return result


# =============================================================================
## Profile Likelihood 
# =============================================================================

# =============================================================================
## Get a profile likelihood confidence interval or upper/lower limit  
#  @code
#  pdf        =  ...
#  dataset    = ...
#  low, high  = interval_PL ( pdf , 'N' , 0.95 , dataset )
#  low_limit  = interval_PL ( pdf , 'N' , 0.95 , dataset , limit = -1 )
#  high_limit = interval_PL ( pdf , 'N' , 0.95 , dataset , limit = +1 )
#  @endcode 
#  @see RooStats::ProfileLikelihoodCalculator
def interval_PL ( pdf , param , level , dataset ,
                  globobs = () , constraints = () , ws = None , limit = None ) :
    """Get a profile likelihood confidence interval
    >>> pdf        = ...
    >>> dataset    = ...
    >>> low, high  = interval_PL ( pdf , 'N' , 0.95 , dataset )
    >>> low_limit  = interval_PL ( pdf , 'N' , 0.95 , dataset , limit = -1 )
    >>> high_limit = interval_PL ( pdf , 'N' , 0.95 , dataset , limit = +1 )
    - see `ROOT.RooStats.ProfileLikelihoodCalculator`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globobs  = globobs  ,
                              ws       = ws       )
    
    plc = ROOT.RooStats.ProfileLikelihoodCalculator ( dataset , mc )
        
    plc.SetConfidenceLevel ( level )

    interval = plc.GetInterval()

    pp       = pdf.parameter  ( param   ) 
    par      = ws.var         ( pp.name )
    
    return get_result ( interval , limit , par , mc , ws ) if ws_new else get_result (  interval , limit , par ) 
                         
# =============================================================================
## Feldman-Cousins 
# =============================================================================

# =============================================================================
## Get  Feldman-Cousins  confidence interval or upper/lower limits 
#  @code
#  pdf       = ...
#  dataset   = ...
#  low, high = conf_interval_FC ( pdf , 'N' , 0.95 , dataset )
#  lower     = conf_interval_FC ( pdf , 'N' , 0.95 , dataset , limit = -1 )
#  upper     = conf_interval_FC ( pdf , 'N' , 0.95 , dataset , limit = +1 )
#  @endcode 
#  @see RooStats::FeldmanCousins
def interval_FC ( pdf , param , level , dataset , 
                  globobs           = ()    ,
                  constraints       = ()    , 
                  ws                = None  ,
                  nbins             = 200   ,
                  adaptive_sampling = True  ,
                  fluctuate         = False , 
                  limit             = None ) :
    """Get Felddman-Cousins confidence interval
    >>> pdf       = ...
    >>> dataset   = ...
    >>> low, high = conf_interval_FC ( pdf , 'N' , 0.95 , dataset )
    - see `ROOT.RooStats.FeldmanCousins`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globobs  = globobs  ,
                              ws       = ws       )
    
    fcc = ROOT.RooStats.FeldmanCousins ( dataset , mc )
    
    fcc.SetConfidenceLevel      ( level     )
    fcc.FluctuateNumDataEntries ( fluctuate ) 
    fcc.UseAdaptiveSampling     ( True if adaptive_sampling else False  )
    fcc.SetNBins                ( nbins )
    
    interval = fcc.GetInterval()

    pp       = pdf.parameter  ( param   ) 
    par      = ws.var         ( pp.name )
    
    return get_result ( interval , limit , par , mc , ws ) if ws_new else get_result ( interval , limit , par ) 


# =============================================================================
## Bayesian calcualtor 
# =============================================================================

# =============================================================================
## Get  Bayesian confidence interval or upper/lower limit
#  @code
#  pdf       = ...
#  dataset   = ...
#  low, high = interval_BC ( pdf , 'N' , 0.95 , dataset )
#  lower     = interval_BC ( pdf , 'N' , 0.95 , dataset , limit = -1 )
#  upper     = interval_BC ( pdf , 'N' , 0.95 , dataset , limit = +1 )
#  @endcode 
#  @see RooStats::BayesianCalculator
def interval_BC ( pdf , param , level , dataset , 
                  globobs           = ()   ,
                  constraints       = ()   ,
                  ws                = None ,
                  prior             = None , 
                  limit             = None ) :
    """Get Bayesian confidence interval
    >>> pdf       = ...
    >>> dataset   = ...
    >>> low, high = conf_interval_BC ( pdf , 'N' , 0.95 , dataset )
    - see `ROOT.RooStats.BayesianCalculator`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globobs  = globobs  ,
                              ws       = ws       )
    
    pp  = pdf.parameter  ( param )
    
    from   ostap.fitting.pdfbasic import APDF1 
    if   prior and isinstance ( prior , APDF1  ) :
        mc.SetPriorPdf( prior.pdf )
    elif prior and isinstance ( prior , ROOT.RooAbsPdf ) :
        mc.SetPriorPdf( prior  )
    else :
        ws.factory("Uniform::prior(%s)" % pp.name )
        mc.SetPriorPdf(ws.pdf("prior"))
    
    bc = ROOT.RooStats.BayesianCalculator ( dataset , mc )
    bc.SetConfidenceLevel  ( level )
    
    bcInt = bc.GetInterval()

    interval = bc.GetInterval()

    par      = None 
    return get_result ( interval , limit , par , mc , ws ) if ws_new else get_result ( interval , limit , par ) 


# =============================================================================
## MCMC Calculator
# =============================================================================


# =============================================================================
## Get MCMC confidence interval
#  @code
#  pdf       = ...
#  dataset   = ...
#  low, high = interval_MCMC ( pdf , 'N' , 0.95 , dataset )
#  lower     = interval_MCMC ( pdf , 'N' , 0.95 , dataset , limit = -1 )
#  upper     = interval_MCMC ( pdf , 'N' , 0.95 , dataset , limit = +1 )
#  @endcode 
#  @see RooStats::MCMCCalculator
def interval_MCMC ( pdf , param , level , dataset , 
                    globobs      = ()     ,
                    constraints  = ()     ,
                    ws           = None   ,
                    nbins        = 200    ,
                    burnsteps    = 500    ,
                    numiters     = 100000 ,
                    leftfraction = 0.5    ,
                    limit        = None   ) :    
    """Get MCMC confidence interval
    >>> pdf       = ...
    >>> dataset   = ...
    >>> low, high = interval_FC   ( pdf , 'N' , 0.95 , dataset )
    >>> lower     = interval_MCMC ( pdf , 'N' , 0.95 , dataset , limit = -1 )
    >>> upper     = interval_MCMC ( pdf , 'N' , 0.95 , dataset , limit = +1 )
    - see `ROOT.RooStats.MCMCCalculator`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globobs  = globobs  ,
                              ws       = ws       )
    
    mcmc= ROOT.RooStats.MCMCCalculator ( dataset , mc )

    mcmc.SetConfidenceLevel       ( level        )
    
    mcmc.SetNumBins               ( nbins        )
    mcmc.SetNumBurnInSteps        ( burnsteps    )
    mcmc.SetNumIters              ( numiters     )
    mcmc.SetLeftSideTailFraction  ( leftfraction )
    
    interval = mcmc.GetInterval()

    pp       = pdf.parameter  ( param   ) 
    par      = ws.var         ( pp.name )
    
    return get_result ( interval , limit , par , mc , ws ) if ws_new else get_result ( interval , limit , par ) 


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
