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
    'ProfileLikelihoodInterval', ## Profile-likelihood confidence interval or upper/lower limit 
    'FeldmanCousinsInterval'   , ## Feldman-Cousins    confidence interval or upper/lower limit 
    'BayesianInterval'         , ## Bayesian           confidence interval or upper/lowee limir
    'MCMCInterval'             , ## MCMC               confidence interval or upper/lower limit 
    ##
    )
# =============================================================================
import ostap.fitting.roofit
from   ostap.core.ostap_types import string_types, integer_types
import ROOT, abc 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.roostats' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
## Create (and fill) ModelConfig & RooWorkspace
#  @code
#  pdf     = ...
#  globars = [ ... ] 
#  mc , ws = create_MC_WS ( pdf , [ 'S' , 'B' ] , dataset = dataset , globobs = () ) 
#  @endcode
#  @see RooStats::ModelConfig 
#  @see RooWorkspace  
def create_MC_WS  ( pdf , params , dataset = None , ws = None , globobs = () , **kwargs ) :
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

    from ostap.utils.cidict import cidict
    from ostap.core.core    import cidict_fun
    
    kw_args = cidict( transform = cidict_fun , **kwargs )
    
    if not ws :
        wsname  = kw_args.pop ( 'ws_name'  , 'WS_%s'            % pdf.name )
        wstitle = kw_args.pop ( 'ws_title' , 'workspace for %s' % pdf.name )
        ws      = ROOT.RooWorkspace         ( wsname  , wstitle )

    mcname  = kw_args.pop ( 'mc_name'  , 'MC_%s'               % pdf.name )
    mctitle = kw_args.pop ( 'mc_title' , 'model-config for %s' % pdf.name )
    mc      = ROOT.RooStats.ModelConfig ( mcname , mctitle , ws )

    mc.SetPdf                  ( pdf.pdf         )
    mc.SetObservables          ( pdf.observables )

    poi = ROOT.RooArgSet       ( *parameters     ) 
    mc.SetParametersOfInterest (  poi            )
    
    nuis = [ v for v in pdf.params ( dataset ) if not v in pdf.vars and not v in poi ]
    nuis = ROOT.RooArgSet      ( *nuis )    
    mc.SetNuisanceParameters   (  nuis )

    pdf.aux_keep.append        ( poi  )
    pdf.aux_keep.append        ( nuis )

    if globobs :
        if isinstance ( globobs , ROOT.RooAbsReal ) : globobs = [ globobs ] 
        gobs = [ pdf.parameter  ( v , dataset ) for v in globobs ]
        gobs = ROOT.RooArgSet   ( *gobs )
        mc.SetGlobalObservables (  gobs )
        pdf.aux_keep.append     (  gobs )

    if kw_args :
        logger.warning ( 'create ModelConfig: Ignore keyword arguments: %s' % [ k for k in kw_args ] )
        
    return mc, ws

# ================================================================================
## Helper (abstract) base class for confidence intervals and limits 
class CLInterval(object)  :
    """Helper (abstract) base class for confidence intervals and limits"""
    
    def __init__ ( self            ,
                   pdf             ,   ## pdf
                   params          ,   ## parameter(s) of interest
                   ws      = None  ,   ## existing workspace
                   dataset = None  ,   ## dataset
                   **kwargs        ) : ## other arguments 

        self.__new_ws = True if not ws else False
        ## use function, later move it into base class 
        self.__mc , self.__ws = create_MC_WS ( pdf     = pdf     ,
                                               params  = params  ,
                                               dataset = dataset ,
                                               ws      = ws      , **kwargs )
        
        self.__dataset  =  dataset
        pp              = pdf.parameter ( params  ) 
        self.__par      = self.__ws.var ( pp.name ) ## parameter as it is in workspace 
        self.__interval = None 

    # =========================================================================
    ## Abstract method to create the interval calculator
    @abc.abstractmethod 
    def make_interval ( self , level , dataset = None ) :
        """Abstract method to create the interval calculator"""
        pass

    # =========================================================================
    ## get the confidence interal at certain confidence level
    #  @code
    #  interval = ...
    #  lower, upper = interval.interval ( 0.90 ) 
    #  @endcode 
    def interval ( self , level , dataset = None ) :
        """get the confidence interal at certain confidence level
        >>> interval = ...
        >>> lower, upper = interval.interval ( 0.90 ) 
        """
        assert 0 < level < 1 , 'Invalid confidence level!'
        self.__interval = self.make_interval ( level , dataset )
        self.__interval.SetConfidenceLevel   ( level )

        par    = self.par
        
        low    = self.__interval.LowerLimit  ( par ) if par else self.__interval.LowerLimit ()
        high   = self.__interval.UpperLimit  ( par ) if par else self.__interval.UpperLimit ()
  
        return low, high 

    # =========================================================================
    ## get the upper/lower limit at certain confidence level
    #  @code
    #  interval =...
    #  upper    = interval.limit ( 0.90 , +1 ) 
    #  lower    = interval.limit ( 0.90 , -1 ) 
    #  @endcode
    def limit    ( self , level , limit , dataset = None ) :
        """Get the upper/lower limit at certain confidence level
        >>> interval =...
        >>> upper    = interval.limit ( 0.90 , +1 ) 
        >>> lower    = interval.limit ( 0.90 , -1 ) 
        """
        assert 0 < level < 1 , 'Invalid confidence level!'
        self.__interval = self.make_interval ( level , dataset ) 
        self.__interval.SetConfidenceLevel   ( 1.0 - 0.5 * ( 1 - level ) ) ## ATTENTION!

        par = self.par 
        if   0 < limit : return self.__interval.UpperLimit ( par ) if par else self.__interval.UpperLimit ()
        elif 0 > limit : return self.__interval.LowerLimit ( par ) if par else self.__interval.LowerLimit ()
                
    # =========================================================================
    ## Get the upper limit at certain confidence level
    #  @code
    #  interval =...
    #  upper    = interval.upper_limit ( 0.90 ) 
    #  @endcode 
    def upper_limit ( self , level , dataset = None ) :
        """Get the upper limit at certain confidence level
        >>> interval =...
        >>> upper    = interval.upper_limit ( 0.90 ) 
        """
        return self.limit ( level , limit = +1 , dataset = dataset )
    
    # =========================================================================
    ## get the lower limit at certain confidence level
    #  @code
    #  interval =...
    #  lower    = interval.lower_limit ( 0.90 ) 
    #  @endcode 
    def lower_limit ( self , level , dataset = None ) :
        """Get the lower limit at certain confidence level
        >>> interval =...
        >>> lower    = interval.lower_limit ( 0.90 ) 
        """        
        return self.limit ( level , limit = -1 , dataset = dataset )

    # =========================================================================
    @property
    def ws ( self )  :
        """'ws' : RooStats workspace"""
        return self.__ws
    @property
    def workspace ( self )  :
        """'workspace' : RooStats workspace"""
        return self.ws
    @property
    def mc ( self )  :
        """'mc' : RooStats ModelConfig object"""
        return self.__mc    
    @property
    def model_config ( self )  :
        """'model_config' : RooStats ModelConfig object"""
        return self.mc
    
    @property
    def par ( self ) :
        """'par' : parameter as it is stored in workspace"""
        return self.__par 
    @par.setter
    def par ( self , value ) :
        self.__par = value 
    
        
# ================================================================================
## Profile Likelihood confidence interval
#  @code
#  interval     = ProfileLikelihoodInterval ( .... )
#  lower, upper = interval.interval ( 0.90 ) 
#  upper        = interval.limit ( 0.90 , +1 ) 
#  lower        = interval.limit ( 0.90 , -1 ) 
#  upper        = interval.upper_limit ( 0.90 ) 
#  lower        = interval.lower_limit ( 0.90 ) 
#  @endcode
#  @see RooStats::ProfileLikelihoodCalculator
class ProfileLikelihoodInterval(CLInterval) :
    """Profile Likelihood confidence interval
    >>> interval     = ProfileLikelihoodInterval ( .... )
    >>> lower, upper = interval.interval ( 0.90 ) 
    >>> upper        = interval.limit ( 0.90 , +1 ) 
    >>> lower        = interval.limit ( 0.90 , -1 ) 
    >>> upper        = interval.upper_limit ( 0.90 ) 
    >>> lower        = interval.lower_limit ( 0.90 ) 
    - see `ROOT.RooStats.ProfileLikelihoodCalculator`
    """
    
    ## create the interval 
    def make_interval ( self , level , dataset = None ) :
    ## create the interval 

        ds = dataset if dataset else self.__dataset 
        assert ds ,           'Invalid dataset!'

        self.__plc = ROOT.RooStats.ProfileLikelihoodCalculator ( ds , self.mc )
        self.__plc.SetConfidenceLevel ( level )
        
        return self.__plc.GetInterval()

    ## ===========================================================================
    #  make a plot
    #  @see RooStats::LikelihoodIntervalPlot
    def plot ( self ) :
        """Make a plot
        - see `ROOT.RooStats.LikelihoodIntervalPlot`
        """
        if self.__interval :
            return  ROOT.RooStats.LikelihoodIntervalPlot( self.__interval )
        
# ================================================================================
## Feldman-Cousins confidence interval
#  @code
#  interval     = FeldmanCousinsInterval ( .... )
#  lower, upper = interval.interval ( 0.90 ) 
#  upper        = interval.limit ( 0.90 , +1 ) 
#  lower        = interval.limit ( 0.90 , -1 ) 
#  upper        = interval.upper_limit ( 0.90 ) 
#  lower        = interval.lower_limit ( 0.90 ) 
#  @endcode
#  @see RooStats::FeldmanCousins
class FeldmanCousinsInterval(CLInterval) :
    """Feldman-Cousins confidence interval
    >>> interval     = FeldmanCousinsInterval ( .... )
    >>> lower, upper = interval.interval ( 0.90 ) 
    >>> upper        = interval.limit ( 0.90 , +1 ) 
    >>> lower        = interval.limit ( 0.90 , -1 ) 
    >>> upper        = interval.upper_limit ( 0.90 ) 
    >>> lower        = interval.lower_limit ( 0.90 ) 
    - see `ROOT.RooStats.FeldmanCousins`
    """
    
    def __init__ ( self              ,
                   pdf               ,   ## pdf
                   params            ,   ## parameter(s) of interest
                   ws        = None  ,   ## existing workspace
                   dataset   = None  ,   ## dataset
                   fluctuate = False ,   ## for RooStats.FeldmanCousins.FluctuateNumDataEntries
                   adaptive  = True  ,   ## for RooStats.FeldmanCousins.UseAdaptiveSampling 
                   nbins     = 200   ,   ## for RooStats.FeldmanCousins.SetNbins                    
                   **kwargs          ) : ## other arguments 
        

        CLInterval.__init__ ( self,
                              pdf     = pdf     ,
                              params  = params  ,
                              ws      = ws      ,
                              dataset = dataset , **kwargs )

        assert isinstance ( nbins , integer_types ) and 10 < nbins ,'Inavlid number of bins!'
        
        self.__fluctuate = True if fluctuate else False
        self.__adaptive  = True if adaptive  else False
        self.__nbins     = nbins
        
    ## create the interval 
    def make_interval ( self , level , dataset = None ) :
        """Create the interval"""
        
        ds = dataset if dataset else self.__dataset 
        assert ds ,           'Invalid dataset!'

        self.__fc = ROOT.RooStats.FeldmanCousins( ds , self.mc )
        self.__fc.SetConfidenceLevel ( level )

        self.__fc.FluctuateNumDataEntries ( self.fluctuate ) 
        self.__fc.UseAdaptiveSampling     ( self.adaptive  )
        self.__fc.SetNBins                ( self.nbins     )
        
        return self.__fc.GetInterval()

    @property
    def fluctuate ( self ) :
        """'fluctuate': parameter for `RooStats.FeldmanCousins.FluctuateNumDataEntries`"""
        return self.__fluctuate 
    @property
    def adaptive  ( self ) :
        """'adaptive': parameter for `RooStats.FeldmanCousins.UseAdaptiveSampling`"""
        return self.__adaptive
    @property
    def nbins     ( self ) :
        """'nbins': parameter for `RooStats.FeldmanCousins.SetNBins`"""
        return self.__nbins


# ================================================================================
## Bayesian confidence interval
#  @code
#  interval     = BayesianInterval ( .... )
#  lower, upper = interval.interval ( 0.90 ) 
#  upper        = interval.limit ( 0.90 , +1 ) 
#  lower        = interval.limit ( 0.90 , -1 ) 
#  upper        = interval.upper_limit ( 0.90 ) 
#  lower        = interval.lower_limit ( 0.90 ) 
#  @endcode
#  @see RooStats::BayesianCalcualtor
class BayesianInterval(CLInterval) :
    """Bayesian confidence interval
    >>> interval     = BAyesianInterval ( .... )
    >>> lower, upper = interval.interval ( 0.90 ) 
    >>> upper        = interval.limit ( 0.90 , +1 ) 
    >>> lower        = interval.limit ( 0.90 , -1 ) 
    >>> upper        = interval.upper_limit ( 0.90 ) 
    >>> lower        = interval.lower_limit ( 0.90 ) 
    - see `ROOT.RooStats.FeldmanCousins`
    """
    
    def __init__ ( self              ,
                   pdf               ,   ## pdf
                   params            ,   ## parameter(s) of interest
                   ws        = None  ,   ## existing workspace
                   dataset   = None  ,   ## dataset
                   prior     = None  ,   ## Bayesin prior 
                   **kwargs          ) : ## other arguments 
        
        ## initialize the Base class
        CLInterval.__init__ ( self,
                              pdf     = pdf     ,
                              params  = params  ,
                              ws      = ws      ,
                              dataset = dataset , **kwargs )
        
        from   ostap.fitting.pdfbasic import APDF1
        
        assert prior is None or \
               isinstance ( prior , APDF1          ) or \
               isinstance ( prior , ROOT.RooAbsPdf ) ,  \
               "Invalid prior!"

        self.__prior = prior
        
        if   prior and isinstance ( prior , APDF1  ) :
            self.mc.SetPriorPdf ( prior.pdf )
        elif prior and isinstance ( prior , ROOT.RooAbsPdf ) :
            self.mc.SetPriorPdf ( prior  )
        else :
            self.ws.factory("Uniform::prior(%s)" % self.par.name )
            self.mc.SetPriorPdf ( self.ws.pdf("prior") )

        ## important to reset it here! 
        self.par = None

    ## create the interval 
    def make_interval ( self , level , dataset = None ) :
        """Create the interval"""
        
        ds = dataset if dataset else self.__dataset 
        assert ds ,           'Invalid dataset!'

        self.__bc = ROOT.RooStats.BayesianCalculator ( ds , self.mc )
        self.__bc.SetConfidenceLevel ( level )
        
        return self.__bc.GetInterval()

    @property
    def prior     ( self ) :
        """'prior': prior for `RooStats.BayesianCalcualtor`"""
        return self.__prior

    ## ===========================================================================
    #  make a plot
    #  @see RooStats::BayesianCalcialtor::GetPosteriorPlot
    def plot ( self ) :
        """Make a plot
        - see `ROOT.RooStats.BayesianCalculator.GetPosteriorPlot`
        """
        if self.__bc : return self.__bc.GetPosteriorPlot() 

# ================================================================================
## MCMC confidence interval
#  @code
#  interval     = MCMCInterval ( .... )
#  lower, upper = interval.interval ( 0.90 ) 
#  upper        = interval.limit ( 0.90 , +1 ) 
#  lower        = interval.limit ( 0.90 , -1 ) 
#  upper        = interval.upper_limit ( 0.90 ) 
#  lower        = interval.lower_limit ( 0.90 ) 
#  @endcode
#  @see RooStats::MCMCCalculator
class MCMCInterval(CLInterval) :
    """MCMC confidence interval
    >>> interval     = MCMCInterval  ( .... )
    >>> lower, upper = interval.interval ( 0.90 ) 
    >>> upper        = interval.limit ( 0.90 , +1 ) 
    >>> lower        = interval.limit ( 0.90 , -1 ) 
    >>> upper        = interval.upper_limit ( 0.90 ) 
    >>> lower        = interval.lower_limit ( 0.90 ) 
    - see `ROOT.RooStats.MCMCCalculator`
    """
    
    def __init__ ( self                  ,
                   pdf                   ,    ## pdf
                   params                ,   ## parameter(s) of interest
                   ws           = None   ,   ## existing workspace
                   dataset      = None   ,   ## dataset
                   burnsteps    = 500    ,   ## for RooStats.MCMCCalculator.SetNumBurnInSteps 
                   iterations   = 100000 ,   ## for RooStats.MCMCCalculator.SetNumIters  
                   leftfraction = 0.5    ,   ## for RooStats.MCMCCalculator.SetNumIters  
                   nbins        = 200    ,   ## for RooStats.MCMCCalculator.SetLeftSideTailFraction
                   **kwargs              ) : ## other arguments 
        
        CLInterval.__init__ ( self,
                              pdf     = pdf     ,
                              params  = params  ,
                              ws      = ws      ,
                              dataset = dataset , **kwargs )

        assert isinstance ( nbins      , integer_types ) and  10 < nbins      ,'Inavlid number of bins!'
        assert isinstance ( burnsteps  , integer_types ) and  10 < burnsteps  , "Invalid nurnsteps!"
        assert isinstance ( iterations , integer_types ) and 100 < iterations , "Invalid nurnsteps!"
        assert 0<= leftfraction <1                                            , 'Invalid leftfraction!'

        self.__nbins        = nbins
        self.__burnsteps    = burnsteps
        self.__iterations   = iterations
        self.__leftfraction = leftfraction
        
    ## create the interval 
    def make_interval ( self , level , dataset = None ) :
        """Create the interval"""
        
        ds = dataset if dataset else self.__dataset 
        assert ds ,           'Invalid dataset!'

        self.__mcmc = ROOT.RooStats.MCMCCalculator ( ds , self.mc )
        self.__mcmc.SetConfidenceLevel ( level )

        self.__mcmc.SetNumBins              ( self.nbins        )
        self.__mcmc.SetNumBurnInSteps       ( self.burnsteps    ) 
        self.__mcmc.SetNumIters             ( self.iterations   )
        self.__mcmc.SetLeftSideTailFraction ( self.leftfraction )

        return self.__mcmc.GetInterval()

    @property
    def nbins        ( self ) :
        """'nbins': parameter for `RooStats.MCMCCalculator.SetNumBins`"""
        return self.__nbins
    @property
    def burnsteps    ( self ) :
        """'burnsteps': parameter for `RooStats.MCMCCalculator.SetNumBurnInSteps`"""
        return self.__burnsteps
    @property
    def iterations   ( self ) :
        """'iterations': parameter for `RooStats.MCMCCalculator.SetNumIter`"""
        return self.__iterations
    @property
    def leftfraction ( self ) :
        """'leftfraction': parameter for `RooStats.MCMCCalculator.SetLeftSideTailFraction`"""
        return self.__leftfraction

    ## ===========================================================================
    #  make a plot
    #  @see RooStats::MCMCIntervalPlot
    def plot ( self ) :
        """Make a plot
        - see `ROOT.RooStats.MCMCIntervallPlot`
        """
        if self.__interval :
            return ROOT.RooStats.MCMCIntervalPlot(mcInt)( self.__interval )
    

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
