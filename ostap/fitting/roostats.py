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
    'ModelConfig'              , ## creator of RooStats ModelConfig 
    'ProfileLikelihoodInterval', ## Profile-likelihood confidence interval or upper/lower limit 
    'FeldmanCousinsInterval'   , ## Feldman-Cousins    confidence interval or upper/lower limit 
    'BayesianInterval'         , ## Bayesian           confidence interval or upper/lowee limir
    'MCMCInterval'             , ## MCMC               confidence interval or upper/lower limit 
    ##
    )
# =============================================================================
from   ostap.core.meta_info    import root_info 
from   ostap.core.ostap_types import string_types, integer_types
import ostap.fitting.roofit
import ROOT, abc 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.roostats' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
## @class ModelConfig
#  Helper class to create `RooStats::ModelConfig`
#  @see RooStats::ModelConfig 
class ModelConfig(object): 
    """Helper class to create `RooStats::ModelConfig`
    - see `ROOT.RooStats.ModelConfig`
    """
    def __init__  ( self                      ,
                    pdf                       , ## PDF 
                    poi                       , ## parameter(s) of interest
                    dataset            = None , ## dataset  (optional)
                    workspace          = None , ## worspace (optional)
                    global_observables = ()   , ## global observables
                    constraints        = ()   , ## contraints 
                    **kwargs                  ) : ## other arguments

        params = poi 
        if   isinstance ( params , ROOT.RooAbsReal ) : params = [ params ]
        elif isinstance ( params , string_types    ) : params = [ params ]

        ## allow soem freedom in arguments 
        from ostap.utils.cidict import cidict
        from ostap.core.core    import cidict_fun
        
        kw_args = cidict ( transform = cidict_fun , **kwargs )
        
        ## (1) create workspace (if needed)
        ws = workspace 
        if ws and isinstance ( ws , ROOT.RooWorkspace ) : pass
        else : 
            wsname  = kw_args.pop ( 'ws_name'  , 'WS_%s'            % pdf.name )
            wstitle = kw_args.pop ( 'ws_title' , 'workspace for %s' % pdf.name )
            ws      = ROOT.RooWorkspace         ( wsname  , wstitle )

        ## (2) create ModelConfig        
        mcname  = kw_args.pop ( 'mc_name'  , 'MC_%s'               % pdf.name )
        mctitle = kw_args.pop ( 'mc_title' , 'model-config for %s' % pdf.name )
        mc      = ROOT.RooStats.ModelConfig ( mcname , mctitle , ws )

        self.__pdf = pdf
        ## (3/4) set PDF and observables 
        mc.SetPdf         ( pdf.roo_pdf     ) ## note: we uise `roo_pdf` here
        mc.SetObservables ( pdf.observables )

        ## (5) set parameters of interest
        pars = [ pdf.parameter ( p , dataset ) for p in params ]
        pois = ROOT.RooArgSet ()
        for p in pars : pois.add   ( p    )
        mc.SetParametersOfInterest ( pois )
        pdf.aux_keep.append        ( pois )
        self.__ps = params , pars, pois
        
        ## (6) Nuisance parameters
        pars = [ v for v in pdf.params ( dataset ) if not v in pdf.vars and not v in pois ]
        nuis = ROOT.RooArgSet    ()
        for p in pars : nuis.add ( p    ) 
        mc.SetNuisanceParameters ( nuis )
        pdf.aux_keep.append      ( nuis )
        self.__ns = pars, nuis 

        ## (7) global observables
        if global_observables :
            if isinstance ( global_observables , ROOT.RooAbsReal ) :
                global_observables = [ global_observables ]                
            pars = [ pdf.parameter  ( v , dataset ) for v in global_observables ]
            gobs = ROOT.RooArgSet   () 
            for p in pars : gobs.add ( p    ) 
            mc.SetGlobalObservables  ( gobs )
            pdf.aux_keep.append      ( gobs )
            self.__go = global_observables, gobs 
            
        ## (8) constraints
        if constraints :
            if isinstance ( constrainst , ROOT.RooAbsReal ) :
                constraints = [ constraints ]
            assert all ( [ isinstance ( c , ROOT.RooAbsPdf ) for c in constraints ] ) , \
                   'Invalid consraints: %s' % constraints
            cnts = ROOT.RooArgSet()
            for c in constraints : cnts.add ( c )
            mc.SetConstraintParameters  ( cnts )
            pdf.aux_keep.append         ( cnts )
            self.__cs = constraints, cnts 
            
        ## (9) define the default dataset 
        self.__dataset = dataset 

        if kw_args :
            logger.warning ( 'create ModelConfig: Ignore keyword arguments: %s' % [ k for k in kw_args ] )

        self.__ws   = ws
        self.__mc   = mc


        
    @property
    def ws ( self ) :
        """'ws' : RooStats workspace (same as 'workspace')"""
        return self.__ws
    @property
    def workspace ( self ) :
        """'workspace' : RooStats workspace (same as 'ws')"""
        return self.ws
    @property
    def mc ( self ) :
        """'mc' : RooStats ModelConfig object (same as 'model_config')"""
        return self.__mc
    @property
    def model_config ( self ) :
        """'model_config' : RooStats ModelConfig object (same as 'mc')"""
        return self.mc
    @property
    def dataset ( self ) :
        """'dataset' : defaukt dataset"""
        return self.__dataset
    
    @property
    def poi ( self ) :
        """'poi' : parameter(s) of interest from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetParametersOfInterest`
        """
        pars = self.mc.GetParametersOfInterest()
        return pars if pars and 0 < len ( pars ) else () 
    @property
    def nuisance ( self ) :
        """'nuisance' : nuisance parameters from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetNuisanceParameters`
        """
        pars = self.mc.GetNuisanceParameters()
        return pars if pars and 0 < len ( pars ) else () 
    @property
    def global_observables ( self ) :
        """'global_observables' : Global observables from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetNuisanceParameters`
        """
        pars = self.mc.GetGobalObservables()
        return pars if pars and 0 < len ( pars ) else () 
    @property
    def constraints  ( self ) :
        """'constraints' : constrain parameters from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetConstraintParameters`
        """
        pars = self.mc.GetConstraintParameters()
        return pars if pars and 0 < len ( pars ) else () 

# ================================================================================
## Helper (abstract) base class for the confidence intervals and limits 
class CLInterval(ModelConfig)  :
    """Helper (abstract) base class for confidence intervals and limits"""
    __metaclass__ = abc.ABCMeta
    
    def __init__ ( self             ,
                   pdf              ,   ## pdf
                   poi              ,   ## parameter(s) of interest
                   workspace = None ,   ## existing workspace
                   dataset   = None ,   ## dataset
                   **kwargs         ) : ## other arguments 

        ModelConfig.__init__ ( self                ,
                               pdf       = pdf     ,
                               poi       = poi     ,
                               dataset   = dataset ,
                               workspace = workspace , **kwargs )
        
        self.__interval     = None 
        self.__calculator   = None
        
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
    def interval ( self , level , par = None , dataset = None ) :
        """get the confidence interal at certain confidence level
        >>> interval = ...
        >>> lower, upper = interval.interval ( 0.90 ) 
        """
        assert 0 < level < 1 , 'Invalid confidence level!'

        if not self.the_interval or level != self.the_interval.ConfidenceLevel or dataset : 
            self.__interval = self.make_interval ( level , dataset )
            
        ##
        par = self.par_from_poi ( par ) 
            
        low    = self.the_interval.LowerLimit  ( par ) if par else self.the_interval.LowerLimit ()
        high   = self.the_interval.UpperLimit  ( par ) if par else self.the_interval.UpperLimit ()
  
        return low, high 

    # =========================================================================
    ## get the upper/lower limit at certain confidence level
    #  @code
    #  interval =...
    #  upper    = interval.limit ( 0.90 , +1 ) 
    #  lower    = interval.limit ( 0.90 , -1 ) 
    #  @endcode
    def limit    ( self , level , limit , par = None, dataset = None ) :
        """Get the upper/lower limit at certain confidence level
        >>> interval =...
        >>> upper    = interval.limit ( 0.90 , +1 ) 
        >>> lower    = interval.limit ( 0.90 , -1 ) 
        """
        assert 0 < level < 1 , 'Invalid confidence level!'

        ## attention! 
        limit_level = 1.0 - 2 * ( 1.0 - level ) 
        if not self.the_interval or limit_level != self.the_interval.ConfidenceLevel or dataset : 
            self.__interval = self.make_interval ( limit_level , dataset ) 
            
        par = self.par_from_poi ( par )
        
        if   0 < limit : return self.the_interval.UpperLimit ( par ) if par else self.the_interval.UpperLimit ()
        elif 0 > limit : return self.the_interval.LowerLimit ( par ) if par else self.the_interval.LowerLimit ()
                
    # =========================================================================
    ## Get the upper limit at certain confidence level
    #  @code
    #  interval =...
    #  upper    = interval.upper_limit ( 0.90 ) 
    #  @endcode 
    def upper_limit ( self , level , par = None , dataset = None ) :
        """Get the upper limit at certain confidence level
        >>> interval =...
        >>> upper    = interval.upper_limit ( 0.90 ) 
        """
        return self.limit ( level , limit = +1 , par = par , dataset = dataset )
    
    # =========================================================================
    ## get the lower limit at certain confidence level
    #  @code
    #  interval =...
    #  lower    = interval.lower_limit ( 0.90 ) 
    #  @endcode 
    def lower_limit ( self , level , par = None , dataset = None ) :
        """Get the lower limit at certain confidence level
        >>> interval =...
        >>> lower    = interval.lower_limit ( 0.90 ) 
        """        
        return self.limit ( level , limit = -1 , par = par , dataset = dataset )

    # =========================================================================
    ## Helper method to get the true parameter from poi 
    def par_from_poi ( self , par ) :
        """Helper method to get the true parameter from poi
        """
        poi = self.poi             
        if   par and isinstance ( par , string_types   ) and poi : return poi [ par      ]
        elif par and isinstance ( par , ROOT.RooAbsArg ) and poi : return poi [ par.name ]
        elif not par and 1 == len ( poi )                        : return poi.front ()

        raise TypeError ( "Invalid setting for 'par' %s | %s" %  ( par , poi ) )

    # =========================================================================
    @property
    def calculator ( self ) :
        """'calculator' : get the actual calculator"""
        return self.__calculator
    @calculator.setter
    def calculator ( self , value ) :
        self.__calculator = value
    @property
    def the_interval ( self ) :
        """'the_interval' : the actual interval object"""
        return self.__interval

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
    # =========================================================================
    ## create the interval 
    def make_interval ( self , level , dataset = None ) :
        """Create the interval"""

        ds = dataset if dataset else self.dataset 
        assert ds ,           'Invalid dataset!'

        pl = ROOT.RooStats.ProfileLikelihoodCalculator ( ds , self.mc )
        pl .SetConfidenceLevel ( level )
        
        self.calculator = pl        
        return pl.GetInterval()

    ## ===========================================================================
    #  make a plot
    #  @see RooStats::LikelihoodIntervalPlot
    def plot ( self ) :
        """Make a plot
        - see `ROOT.RooStats.LikelihoodIntervalPlot`
        """
        if self.the_interval :
            return  ROOT.RooStats.LikelihoodIntervalPlot( self.the_interval )
        
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
                   poi               ,   ## parameter(s) of interest
                   dataset   = None  ,   ## dataset
                   workspace = None  ,   ## existing workspace
                   fluctuate = False ,   ## for RooStats.FeldmanCousins.FluctuateNumDataEntries
                   adaptive  = True  ,   ## for RooStats.FeldmanCousins.UseAdaptiveSampling 
                   nbins     = 200   ,   ## for RooStats.FeldmanCousins.SetNbins                    
                   **kwargs          ) : ## other arguments 
        
        CLInterval.__init__ ( self,
                              pdf       = pdf       ,
                              poi       = poi       ,
                              dataset   = dataset   ,
                              workspace = workspace , **kwargs )

        assert isinstance ( nbins , integer_types ) and 10 < nbins ,'Inavlid number of bins!'
        
        self.__fluctuate = True if fluctuate else False
        self.__adaptive  = True if adaptive  else False
        self.__nbins     = nbins
        
    ## create the interval 
    def make_interval ( self , level , dataset = None ) :
        """Create the interval"""
        
        ds = dataset if dataset else self.dataset 
        assert ds ,           'Invalid dataset!'

        fc = ROOT.RooStats.FeldmanCousins( ds , self.mc )
        fc.SetConfidenceLevel      ( level          )
        fc.FluctuateNumDataEntries ( self.fluctuate ) 
        fc.UseAdaptiveSampling     ( self.adaptive  )
        fc.SetNBins                ( self.nbins     )

        self.calculator = fc         
        return fc.GetInterval()

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

    # =========================================================================
    ## visualize the interval
    #  - inspired by rs401c_FeldmanCousins.py
    #  @code
    #  fci = ...
    #  graph = fci.plot()
    #  graph.Draw('ap')
    #  @endcode 
    def plot ( self ) :
        """Visualize the interval
        - inspired by rs401c_FeldmanCousins.py
        >>> fci = ...
        >>> graph = fci.plot()
        >>> graph.Draw('ap')
        """

        if root_info < (6,18) :
            logger.warning ( "No plots from Feldman-Cousins for ROOT<(6.18)" )
            return 
            
        if self.calculator and self.the_interval :
            
            import ostap.fitting.dataset
            import ostap.histos.graphs
            
            fc  = self.calculator 
            fci = self.the_interval

            gr1 = ROOT.TGraph ()
            gr2 = ROOT.TGraph ()
            gr1.red  ()
            gr2.blue ()
        
            ps  = fc.GetPointsToScan()
            for entry in ps :            
                point = float ( entry[0] ) 
                if fci.IsInInterval ( entry ) : gr1.append ( point , 1 )
                else                          : gr2.append ( point , 0 )

            mgr = ROOT.TMultiGraph()
            mgr.Add ( gr1 )
            mgr.Add ( gr2 )

            del ps 
            return mgr
        
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
    >>> interval     = BayesianInterval ( .... )
    >>> lower, upper = interval.interval ( 0.90 ) 
    >>> upper        = interval.limit ( 0.90 , +1 ) 
    >>> lower        = interval.limit ( 0.90 , -1 ) 
    >>> upper        = interval.upper_limit ( 0.90 ) 
    >>> lower        = interval.lower_limit ( 0.90 ) 
    - see `ROOT.RooStats.BayesianCalculator`
    """
    
    def __init__ ( self              ,
                   pdf               ,   ## pdf
                   poi               ,   ## parameter(s) of interest
                   dataset   = None  ,   ## dataset
                   workspace = None  ,   ## existing workspace
                   prior     = None  ,   ## Bayesin prior 
                   **kwargs          ) : ## other arguments 
        
        ## initialize the Base class
        CLInterval.__init__ ( self,
                              pdf       = pdf       ,
                              poi       = poi       ,
                              dataset   = dataset   ,
                              workspace = workspace , **kwargs )

        np = len ( self.poi )
        assert 1 == np , 'Bayesing interval works only for 1 poi!'

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
            par = self.poi.front()
            self.ws.factory("Uniform::prior(%s)" % par.name )
            self.mc.SetPriorPdf ( self.ws.pdf("prior") )

    # =========================================================================
    ## create the interval 
    def make_interval ( self , level , dataset = None ) :
        """Create the interval"""
        
        ds = dataset if dataset else self.dataset 
        assert ds ,           'Invalid dataset!'

        bc = ROOT.RooStats.BayesianCalculator ( ds , self.mc )
        bc.SetConfidenceLevel ( level )

        self.calculator = bc 
        return bc.GetInterval()

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
        if self.calculator :
            return self.calculator.GetPosteriorPlot() 

    # =========================================================================
    ## Helper method to get the true parameter from poi 
    def par_from_poi ( self , par ) :
        """Helper method to get the true parameter from poi
        """
        return None

# ================================================================================
## Marcov Chain MC confidence interval
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
    """Markov Chain confidence interval
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
                   poi                   ,   ## parameter(s) of interest
                   dataset      = None   ,   ## dataset
                   workspace    = None   ,   ## existing workspace
                   burnsteps    = 500    ,   ## for RooStats.MCMCCalculator.SetNumBurnInSteps 
                   iterations   = 100000 ,   ## for RooStats.MCMCCalculator.SetNumIters  
                   leftfraction = 0.5    ,   ## for RooStats.MCMCCalculator.SetNumIters  
                   nbins        = 200    ,   ## for RooStats.MCMCCalculator.SetLeftSideTailFraction
                   **kwargs              ) : ## other arguments 
        
        CLInterval.__init__ ( self,
                              pdf       = pdf       ,
                              poi       = poi       ,
                              dataset   = dataset   ,
                              workspace = workspace , **kwargs )

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
        
        ds = dataset if dataset else self.dataset 
        assert ds ,           'Invalid dataset!'

        mcmc = ROOT.RooStats.MCMCCalculator ( ds , self.mc )
        mcmc.SetConfidenceLevel ( level )

        mcmc.SetNumBins              ( self.nbins        )
        mcmc.SetNumBurnInSteps       ( self.burnsteps    ) 
        mcmc.SetNumIters             ( self.iterations   )
        mcmc.SetLeftSideTailFraction ( self.leftfraction )

        self.calculator = mcmc 
        return mcmc.GetInterval()

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
        if self.the_interval :
            return ROOT.RooStats.MCMCIntervalPlot( self.the_interval )
    

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


    assert issubclass ( ProfileLikelihoodInterval , CLInterval ) , \
           'ProfileLikelihoodInterval is not subsclas of CLInterval'
    
    assert issubclass ( FeldmanCousinsInterval    , CLInterval ) , \
           'FeldmanCousinsInterval    is not subsclas of CLInterval'

    assert issubclass ( BayesianInterval         , CLInterval ) , \
           'BayesianInterval          is not subsclas of CLInterval'

    assert issubclass ( MCMCInterval             , CLInterval ) , \
           'MCMCInterval              is not subsclas of CLInterval'

    
# =============================================================================
##                                                                      The END 
# =============================================================================
