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
    'create_MC_WS'       , ## helper fnuction to create (and fill) ModelConfig & RooWorkspace
    ##
    'conf_interval_PL'   , ## Profile-likelihood confidence interval
    'conf_interval_FC'   , ## Feldman-Cousins    confidence interval
    'conf_interval_BC'   , ## Bayesian           confidence interval
    'conf_interval_MCMC' , ## MCMC               confidence interval
    ##
    'upper_limit_PL'     , ## Profile-likelihood upper limit 
    'upper_limit_FC'     , ## Feldman-Cousins    upper limit 
    'upper_limit_BC'     , ## Bayesian           upper limit 
    'upper_limit_MCMC'   , ## MCMC               upper limit 
    ##
    'lower_limit_PL'     , ## Profile-likelihood lower limit 
    'lower_limit_FC'     , ## Feldman-Cousins    lower limit 
    'lower_limit_BC'     , ## Bayesian           lower limit 
    'lower_limit_MCMC'   , ## MCMC               lower limit 
    ##    
    )
# =============================================================================
import ostap.fitting.roofit 
from   ostap.core.ostap_types import stgirng_types 

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
    >>> globars = [ ... ] 
    >>> mc , ws = create_MC_WS ( pdf , [ 'S' , 'B' ] , globvars = globvars , dataset = dataset ) 
    - see `ROOT.RooStats.ModelConfig` 
    - see `ROOT.RooWorkspace`  
    """
    if   isinstance ( params , ROOT.RooAbsReal ) : params = [ params ]
    elif isinstance ( params , string_types    ) : params = [ params ]
    
    parameters = [ pdf.parameter ( p , dataset ) for p in params ]
    
    if not ws : ws = ROOT.RooWorkSpace         ( 'WS%s' % pdf.name  , 'workspace' )
        
    mc = ROOT.RooStats.ModelConfig ( 'MC' % pdf.name  , 'model-config for %s' % pdf  , ws )

    mc.SetPdf                  ( pdf.pdf     )
    mc.SetObservables          ( pdf.vars    )

    poi = ROOT.RooArgSet       ( *parameters ) 
    mc.SerParametersOfInterest (  poi        )
    
    nuis = [ v for v in pdf.params ( dataset ) if not v in pdf.vars and not v in poi ]
    nuis = ROOT.RooAbsArg      ( *nuis )    
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
## Profile Likelihood 
# =============================================================================

# =============================================================================
## Get a profile likelihood confidence interval
#  @code
#  pdf       = ...
#  dataset   = ...
#  low, high = conf_interval_PL ( pdf , 'N' , 0.95 , dataset )
#  @endcode 
#  @see RooStats::ProfileLikelihoodCalculator
def conf_interval_PL ( pdf , param , level , dataset , globvars = () , constraints = () , ws = None ) :
    """Get a profile likelihood confidence interval
    >>> pdf       = ...
    >>> dataset   = ...
    >>> low, high = conf_interval_PL ( pdf , 'N' , 0.95 , dataset )
    - see `ROOT.RooStats.ProfileLikelihoodCalculator`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globvars = globvars ,
                              ws       = ws       )
    
    plc = ROOT.RooStats.ProfileLikelihoodCalculator ( dataset , mc )
        
    plc.SetConfidenceLevel ( level )        
    plInt = plc.GetInterval()
    
    pp    = self.parameter  ( param ) 
    par   = ws.var ( pp.name )
    
    low   = plInt.LowerLimit ( par )
    high  = plInt.UpperLimit ( par )

    if ws_new : 
        del ws
        del ws 
        
    return low, high 


# =============================================================================
## Get an upper limit using profile likelihood 
#  @code
#  pdf       = ...
#  dataset   = ...
#  limit     = upper_limit_PL ( pdf , 'N' , 0.95 , dataset )
#  @endcode 
#  @see RooStats::ProfileLikelihoodCalculator
def upper_limit_PL ( pdf , param  , level , dataset , globvars = () , constraints = () , ws = None ) :
    """Get a profile likelihood confidence interval
    >>> pdf       = ...
    >>> dataset   = ...
    >>> low, high = conf_interval_PL ( pdf , 'N' , 0.95 , dataset )
    - see `ROOT.RooStats.ProfileLikelihoodCalculator`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globvars = globvars ,
                              ws       = ws       )

    plc = ROOT.RooStats.ProfileLikelihoodCalculator ( dataset , mc )
        
    plc.SetConfidenceLevel ( 1 - 2 * ( 1 - level ) )   ## ATTENTION! 
    plInt = plc.GetInterval()
    
    pp    = self.parameter  ( param ) 
    par   = ws.var ( pp.name )
    
    high  = plInt.UpperLimit ( par )

    if ws_new : 
        del ws
        del ws 

    return high 


# =============================================================================
## Get a lower limit using profile likelihood 
#  @code
#  pdf       = ...
#  dataset   = ...
#  limit     = lower_limit_PL ( pdf , 'N' , 0.95 , dataset )
#  @endcode 
#  @see RooStats::ProfileLikelihoodCalculator
def lower_limit_PL ( pdf , param  , level , dataset , globvars = () , constraints = () , ws = None ) :
    """Get a profile likelihood confidence interval
    >>> pdf       = ...
    >>> dataset   = ...
    >>> limit     = lower_limit_PL ( pdf , 'N' , 0.95 , dataset )
    - see `ROOT.RooStats.ProfileLikelihoodCalculator`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globvars = globvars ,
                              ws       = ws       )

    plc = ROOT.RooStats.ProfileLikelihoodCalculator ( dataset , mc )
        
    plc.SetConfidenceLevel ( 1 - 2 * ( 1 - level ) )   ## ATTENTION! 
    plInt = plc.GetInterval()
    
    pp    = self.parameter  ( param ) 
    par   = ws.var ( pp.name )
    
    low   = plInt.UpperLimit ( par )

    if ws_new : 
        del ws
        del ws 

    return low 


# =============================================================================
## Feldman-Cousins 
# =============================================================================


# =============================================================================
## Get  Feldman-Cousins  confidence interval
#  @code
#  pdf       = ...
#  dataset   = ...
#  low, high = conf_interval_FC ( pdf , 'N' , 0.95 , dataset )
#  @endcode 
#  @see RooStats::FeldmanCousins
def conf_interval_FC ( pdf , param , level , dataset , 
                       globvars          = ()   ,
                       constraints       = ()   ,
                       ws                = None ,
                       nbins             = 100  ,
                       adaptive_sampling = True ) :
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
                              globvars = globvars ,
                              ws       = ws       )
    
    fcc = ROOT.RooStats.FeldmanCousins ( dataset , mc )
    
    fcc.SetConfidenceLevel  ( level )
    
    fcc.UseAdaptiveSampling ( True if adaptive_sampling else False  )
    fcc.SetNBins            ( nbins )
    
    fcInt = fcc.GetInterval()
    
    pp    = self.parameter  ( param ) 
    par   = ws.var ( pp.name )
    
    low   = fcInt.LowerLimit ( par )
    high  = fcInt.UpperLimit ( par )

    if ws_new : 
        del ws
        del ws 
        
    return low, high 


# =============================================================================
## Get  Feldman-Cousins upper limit
#  @code
#  pdf       = ...
#  dataset   = ...
#  limit     = upper_limit_FC ( pdf , 'N' , 0.95 , dataset )
#  @endcode 
#  @see RooStats::FeldmanCousins
def upper_limit_FC ( pdf , param , level , dataset , globvars = () , constraints = () , ws = None ) :
    """Get Felddman-Cousins upper limit
    >>> pdf       = ...
    >>> dataset   = ...
    >>> limit     = upper_limit_FC ( pdf , 'N' , 0.95 , dataset )
    - see `ROOT.RooStats.FeldmanCousins`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globvars = globvars ,
                              ws       = ws       )
    
    fcc = ROOT.RooStats.FeldmanCousins ( dataset , mc )

    
    fcc.SetConfidenceLevel ( 1  - 2.0 * ( 1 - level )  )  ## attention!
    fcInt = fcc.GetInterval()
    
    pp    = self.parameter  ( param ) 
    par   = ws.var ( pp.name )
    
    high  = fcInt.UpperLimit ( par )

    if ws_new : 
        del ws
        del ws 
        
    return high

# =============================================================================
## Get  Feldman-Cousins lowewr limit
#  @code
#  pdf       = ...
#  dataset   = ...
#  limit     = lower_limit_FC ( pdf , 'N' , 0.95 , dataset )
#  @endcode 
#  @see RooStats::FeldmanCousins
def upper_limit_FC ( pdf , param , level , dataset , globvars = () , constraints = () , ws = None ) :
    """Get Felddman-Cousins upper limit
    >>> pdf       = ...
    >>> dataset   = ...
    >>> limit     = upper_limit_FC ( pdf , 'N' , 0.95 , dataset )
    - see `ROOT.RooStats.FeldmanCousins`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globvars = globvars ,
                              ws       = ws       )
    
    fcc = ROOT.RooStats.FeldmanCousins ( dataset , mc )

    
    fcc.SetConfidenceLevel ( 1  - 2.0 * ( 1 - level )  )  ## attention!
    fcInt = fcc.GetInterval()
    
    pp    = self.parameter  ( param ) 
    par   = ws.var ( pp.name )
    
    low   = fcInt.LowerLimit ( par )

    if ws_new : 
        del ws
        del ws 
        
    return low



# =============================================================================
## Bayesian calcualtor 
# =============================================================================

# =============================================================================
## Get  Bayesian confidence interval
#  @code
#  pdf       = ...
#  dataset   = ...
#  low, high = conf_interval_BC ( pdf , 'N' , 0.95 , dataset )
#  @endcode 
#  @see RooStats::BayesianCalculator
def conf_interval_BC ( pdf , param , level , dataset , 
                       globvars          = ()   ,
                       constraints       = ()   ,
                       ws                = None ,
                       nbins             = 200  ,
                       adaptive_sampling = True ) :
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
                              globvars = globvars ,
                              ws       = ws       )

    
    pp  = self.parameter  ( param ) 
    ws.factory("Uniform::prior(%s)" % pp.name )
    mc.SetPriorPdf(ws.pdf("prior"))
    
    bc = ROOT.RooStats.Bayesian ( dataset , mc )
    bc.SetConfidenceLevel  ( level )
    
    bcInt = bc.GetInterval()
    
    par   = ws.var ( pp.name )
    
    low   = bcInt.LowerLimit ( par )
    high  = bcInt.UpperLimit ( par )

    if ws_new : 
        del ws
        del ws 
        
    return low, high 


# =============================================================================
## Get  Bayesian upper limit 
#  @code
#  pdf       = ...
#  dataset   = ...
#  limit     = upper_limit_BC ( pdf , 'N' , 0.95 , dataset )
#  @endcode 
#  @see RooStats::BayesianCalculator
def upper_limit_BC ( pdf , param , level , dataset , 
                     globvars          = ()   ,
                     constraints       = ()   ,
                     ws                = None ) :
    """Get Bayesian upper limit
    >>> pdf       = ...
    >>> dataset   = ...
    >>> limit     = upper_limit_BC ( pdf , 'N' , 0.95 , dataset )
    - see `ROOT.RooStats.BayesianCalculator`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globvars = globvars ,
                              ws       = ws       )

    pp  = self.parameter  ( param ) 
    ws.factory("Uniform::prior(%s)" % pp.name )
    mc.SetPriorPdf(ws.pdf("prior"))
    
    bc = ROOT.RooStats.Bayesian ( dataset , mc )
    bc.SetConfidenceLevel  ( 1 - 2.0 * ( 1 - level )  ) ## ATTENTION 
    
    bcInt = bc.GetInterval()
    
    high  = bcInt.UpperLimit ()

    if ws_new : 
        del ws
        del ws 
        
    return high 


# =============================================================================
## Get  Bayesian lower limit 
#  @code
#  pdf       = ...
#  dataset   = ...
#  limit     = lower_limit_BC ( pdf , 'N' , 0.95 , dataset )
#  @endcode 
#  @see RooStats::BayesianCalculator
def lower_limit_BC ( pdf , param , level , dataset , 
                     globvars          = ()   ,
                     constraints       = ()   ,
                     ws                = None )
    """Get Bayesian lowerlimit
    >>> pdf       = ...
    >>> dataset   = ...
    >>> limit     = lower_limit_BC ( pdf , 'N' , 0.95 , dataset )
    - see `ROOT.RooStats.BayesianCalculator`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globvars = globvars ,
                              ws       = ws       )

    
    pp  = self.parameter  ( param ) 
    ws.factory("Uniform::prior(%s)" % pp.name )
    mc.SetPriorPdf(ws.pdf("prior"))
    
    bc = ROOT.RooStats.Bayesian ( dataset , mc )
    bc.SetConfidenceLevel  ( 1 - 2.0 * ( 1 - level )  ) ## ATTENTION 
    
    bcInt = bc.GetInterval()
    
    low   = bcInt.LowerLimit ()

    if ws_new : 
        del ws
        del ws 
        
    return low


# =============================================================================
## MCMC Calculator
# =============================================================================



# =============================================================================
## Get MCMC confidence interval
#  @code
#  pdf       = ...
#  dataset   = ...
#  low, high = conf_interval_MCMC ( pdf , 'N' , 0.95 , dataset )
#  @endcode 
#  @see RooStats::MCMCCalculator
def conf_interval_MCMC ( pdf , param , level , dataset , 
                         globvars     = ()     ,
                         constraints  = ()     ,
                         ws           = None   ,
                         nbins        = 200    ,
                         burnsteps    = 500    ,
                         numiters     = 100000 ,
                         leftfraction = 0.5    ) :    
    """Get MCMC confidence interval
    >>> pdf       = ...
    >>> dataset   = ...
    >>> low, high = conf_interval_FC ( pdf , 'N' , 0.95 , dataset )
    - see `ROOT.RooStats.MCMCCalculator`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globvars = globvars ,
                              ws       = ws       )
    
    mcmc= ROOT.RooStats.MCMCCalculator ( dataset , mc )

    mcmc.SetConfidenceLevel       ( level        )
    
    mcmc.SetNBins                 ( nbins        )
    mcmc.SetNumBurnInSteps        ( burnsteps    )
    mcmc.SetNumIters              ( numiters     )
    mcmc.SetLedtSideTailFractgion ( leftFraction )
    
    mcInt = mcmc.GetInterval()
    
    pp    = self.parameter  ( param ) 
    par   = ws.var ( pp.name )
    
    low   = mcmcInt.LowerLimit ( par )
    high  = mcmcInt.UpperLimit ( par )

    if ws_new : 
        del ws
        del ws 
        
    return low, high 

# =============================================================================
## Get MCMC upper limit
#  @code
#  pdf       = ...
#  dataset   = ...
#  lmit      = upper_limit_MCMC ( pdf , 'N' , 0.95 , dataset )
#  @endcode 
#  @see RooStats::MCMCCalculator
def upper_limit_MCMC ( pdf , param , level , dataset , 
                       globvars     = ()     ,
                       constraints  = ()     ,
                       ws           = None   ,
                       nbins        = 200    ,
                       burnsteps    = 500    ,
                       numiters     = 100000 ,
                       leftfraction = 0.5    ) :    
    """Get MCMC upper_limit
    >>> pdf       = ...
    >>> dataset   = ...
    >>> limit     = upper_limit_FC ( pdf , 'N' , 0.95 , dataset )
    - see `ROOT.RooStats.MCMCCalculator`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globvars = globvars ,
                              ws       = ws       )
    
    mcmc= ROOT.RooStats.MCMCCalculator ( dataset , mc )

    mcmc.SetConfidenceLevel       ( 1 - 2.0 * ( 1 - level ) ) 
    
    mcmc.SetNBins                 ( nbins        )
    mcmc.SetNumBurnInSteps        ( burnsteps    )
    mcmc.SetNumIters              ( numiters     )
    mcmc.SetLedtSideTailFractgion ( leftFraction )
    
    mcInt = mcmc.GetInterval()
    
    pp    = self.parameter  ( param ) 
    par   = ws.var ( pp.name )
    
    high  = mcmcInt.UpperLimit ( par )

    if ws_new : 
        del ws
        del ws 
        
    return high 

# =============================================================================
## Get MCMC lower limit
#  @code
#  pdf       = ...
#  dataset   = ...
#  lmit      = lower_limit_MCMC ( pdf , 'N' , 0.95 , dataset )
#  @endcode 
#  @see RooStats::MCMCCalculator
def lower_limit_MCMC ( pdf , param , level , dataset , 
                       globvars     = ()     ,
                       constraints  = ()     ,
                       ws           = None   ,
                       nbins        = 200    ,
                       burnsteps    = 500    ,
                       numiters     = 100000 ,
                       leftfraction = 0.5    ) :    
    """Get MCMC upper_limit
    >>> pdf       = ...
    >>> dataset   = ...
    >>> limit     = lower_limit_FC ( pdf , 'N' , 0.95 , dataset )
    - see `ROOT.RooStats.MCMCCalculator`
    """
    assert 0 < level < 1 , 'Invalid confidence level %s' % level

    ws_new  = True if not ws else False  
    mc , ws = create_MC_WS  ( pdf      = pdf      ,
                              params   = param    ,
                              dataset  = dataset  ,
                              globvars = globvars ,
                              ws       = ws       )
    
    mcmc= ROOT.RooStats.MCMCCalculator ( dataset , mc )

    mcmc.SetConfidenceLevel       ( 1 - 2.0 * ( 1 - level ) ) 
    
    mcmc.SetNBins                 ( nbins        )
    mcmc.SetNumBurnInSteps        ( burnsteps    )
    mcmc.SetNumIters              ( numiters     )
    mcmc.SetLedtSideTailFractgion ( leftFraction )
    
    mcInt = mcmc.GetInterval()
    
    pp    = self.parameter  ( param ) 
    par   = ws.var ( pp.name )
    
    low   = mcmcInt.LowerLimit ( par )

    if ws_new : 
        del ws
        del ws 
        
    return low


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
