#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/roostats.py
#  Set of useful basic utilities to dela with RooStats
#
# - confidence interval estimation
# - confidence interval and limits
# - hypothesis tests
# - Brasilian  plots 
#
#  @thanks to Dmitry Golubkov for th egreat heap and examples 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2023-01-17
# =============================================================================
""" Set of useful basic utilities to deal with RooStats

- confidence interval estimation
- confidence interval and limits
- hypothesis tests
- Brasilian  plots 

@thanks to Dmitry Golubkov for
"""
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
    'AsymptoticCalculator'     , ## AsymptoticCalculator   for limits and intervals
    'FrequentistCalculator'    , ## Frequentist calcualtor for limits and nitervals 
    'HypoTestInverter'         , ## 
    ##
    'BrasilBand'               , ## utility to porduce Brasil-band plots 
    )
# =============================================================================
from   ostap.core.meta_info   import root_info 
from   ostap.core.ostap_types import string_types, integer_types, sequence_types 
import ostap.fitting.roofit
from   ostap.fitting.pdfbasic import APDF1
import ROOT, abc, sys  
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.roostats' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
if (3,0) <= sys.version_info : from itertools import  zip_longest
else                         : from itertools import izip_longest as zip_longest 
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
                    observables        = ()   , ## observables 
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
        self.__ws   = ws

        ## (2) create ModelConfig        
        mcname    = kw_args.pop ( 'name'  , 'MC_%s'               % pdf.name )
        mctitle   = kw_args.pop ( 'title' , 'model-config for %s' % pdf.name )
        mc        = ROOT.RooStats.ModelConfig ( mcname , mctitle , ws )
        self.__mc = mc


        self.__pdf = pdf
        
        if   isinstance ( pdf , APDF1          ) :
            self.__raw_pdf = pdf.pdf 
        elif isinstance ( pdf , ROOT.RooAbsPdf ) and observables :
            self.__raw_pdf = pdf
        else :
            raise TypeError ( "Inavlid seting of pdf/observables!" ) 

        ## final pdf? (not yet...) 
        self.__final_pdf = self.__raw_pdf

        ## (0) We need to start from constraints and modify our PDF on-fly 
        if constraints :
            if isinstance ( constraints , ROOT.RooAbsReal ) :
                constraints = [ constraints ]
            assert all ( [ isinstance ( c , ROOT.RooAbsPdf ) for c in constraints ] ) , \
                   'Invalid constraints: %s' % str ( constraints ) 
            cnts = ROOT.RooArgSet()
            for c in constraints : cnts.add ( c )

            clst = ROOT.RooArgList ()
            clst.add ( self.__raw_pdf )
            for c in cnts : clst.add ( c )
            self.__final_pdf = ROOT.RooProdPdf ( '%s_constrained' % pdf.name ,
                                                 "PDF with %d constraints" % len ( cnts ) , clst )
            self.__cs   = constraints, cnts, clst

            ## attention, redefine/update 
            constraints = cnts
            
        ## propagate PDF to ModelConfig/Workspace 
        if isinstance   ( pdf , ROOT.RooAbsPdf ) and observables :
            ## (3) set PDF  
            mc.SetPdf         ( self.__final_pdf ) ## note: we use bare RooAbsPdf here 
        elif isinstance ( pdf , APDF1 ) :  
            ## (3/4) set PDF and observables 
            mc.SetPdf         ( self.__final_pdf ) ## note: we use bare RooAbsPdf here 
            mc.SetObservables ( pdf.observables  )
            
        ## (4) observables
        if observables :
            if isinstance ( observables , ROOT.RooAbsReal ) :
                observables = [ observables ]                
            pars = [ self.pdf_param ( v , dataset ) for v in observables ]
            obs  = ROOT.RooArgSet   () 
            for p in pars : obs.add ( p    ) 
            mc.SetObservables  ( obs )
            self.__ob = observables, obs 

        ## get them back from workspace 
        observables = self.observables
        
        ## (5) set parameters of interest
        pars = [ pdf.parameter ( p , dataset ) for p in params ]
        pois = ROOT.RooArgSet ()
        for p in pars : pois.add   ( p    )
        mc.SetParametersOfInterest ( pois )
        self.__ps = params , pars, pois
        
        ## (6) Nuisance parameters
        pars = [ v for v in self.pdf_params ( dataset ) if  not v in observables and not v in pois ]
        nuis = ROOT.RooArgSet    ()
        for p in pars : nuis.add ( p    ) 
        mc.SetNuisanceParameters ( nuis )
        self.__ns = pars, nuis 

        ## (7) global observables
        if global_observables :
            if isinstance ( global_observables , ROOT.RooAbsReal ) :
                global_observables = [ global_observables ]                
            pars = [ self.pdf_param  ( v , dataset ) for v in global_observables ]
            gobs = ROOT.RooArgSet   () 
            for p in pars : gobs.add ( p    ) 
            mc.SetGlobalObservables  ( gobs )
            self.__go = global_observables, gobs 
            
        ## (8) constraints
        if constraints :
            mc.SetConstraintParameters  ( constraints )
            
        ## (9) define the default dataset 
        self.__dataset = dataset 

        if kw_args :
            logger.warning ( 'create ModelConfig: Ignore keyword arguments: %s' % [ k for k in kw_args ] )

        
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
    def observables ( self ) :
        """'observables' : Global observables from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetObservables`
        """
        pars = self.mc.GetObservables()
        return pars if pars and 0 < len ( pars ) else () 
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
        - see `ROOT.RooStats.ModelConfig.GetGlobalObservables`
        """
        pars = self.mc.GetGlobalObservables()
        return pars if pars and 0 < len ( pars ) else () 
    @property
    def constraints  ( self ) :
        """'constraints' : constrain parameters from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetConstraintParameters`
        """
        pars = self.mc.GetConstraintParameters()
        return pars if pars and 0 < len ( pars ) else () 

    @property
    def snapshot ( self ) :
        """'snapshot' : get/set snapshot fomr the model/workspacee
        - see `ROOT.RooStats.ModelConfig.GetSnapshot`
        - see `ROOT.RooStats.ModelConfig.SetSnapshot`
        """
        pars = self.mc.GetSnapshot()
        return pars if pars and 0 < len ( pars ) else ()
    @snapshot.setter
    def snapshot ( self , values ) :
        
        if   isinstance ( values , ROOT.RooArgSet  ) : return self.mc.SetSnapshot (  values ) 
        elif isinstance ( values , ROOT.RooAbsReal ) : return self.mc.SetSnapshot ( ROOT.RooArgSet ( values ) )
        elif isinstance ( values , sequence_types  ) :

            vv = [ v for v in values ]
            if  all ( isinstance ( v , ROOT.RooAbsReal ) for v in vv ) :
                vs = ROOT.RooArgSet()
                for v in values : vs.add ( vs )
                return self.mc.SetSnapshot ( vs )

        raise TypeError ( 'Invalid type for snapshot %s' % type ( values )  ) 

    # =========================================================================
    ## helper function to get parameters from PDF 
    def pdf_params ( self , dataset = None ) :
        """helepr function to get the parameters from PDF"""
        pdf = self.__final_pdf
        return pdf.getParameters ( 0 ) if dataset is None else  pdf.getParameters ( dataset )
        
    # =========================================================================
    ## helper function to get parameters from PDF 
    def pdf_param  ( self , param , dataset = None ) :
        """helepr function to get the parameters from PDF"""
        params = self.pdf_params ( dataset )
        ## already parameter 
        if isinstance ( param , ROOT.RooAbsReal ) and param in params : return param
        ## loop by name 
        if isinstance ( param , string_types    ) :
            for p in params :
                if p.name == param       : return p
        elif isinstance ( param , ROOT.RooAbsReal ) :
            for p in params :
                if p.name == param.name  : return p
        ## try with workspace
        vv = self.var ( param )
        if vv : return vv
        raise KeyError ( "No parameter %s/%s defined" % ( param , type ( param ) ) )
    
    # =========================================================================
    ## get a variable from workspace
    #  @code
    #  mc = ...
    #  v1 = mc.var ( 'mass' ) ## get by name
    # 
    #  var = 
    #  v1 = mc.var ( 'var   ) ## get by var     
    #  @endcode
    def var ( self , variable ) :
        """Get a variable from workspace
        >>> mc = ...
        >>>> v1 = mc.var ( 'mass' ) ## get by name
        >>> var = 
        >>> v1 = mc.var ( 'var   ) ## get by var     
        """
        if isinstance ( variable , ROOT.RooAbsReal ) :
            return self.var ( variable.name )
        return self.ws.var ( variable ) 

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

        fc = ROOT.RooStats.FeldmanCousins ( ds , self.mc )
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

    # =========================================================================
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
## Calculators
# =============================================================================
## @class Calculator
#  base class for Calculators 
class Calculator (object) :
    """base class for Calculators 
    - see `ROOT.RooStats.AsymptoticCalculator`
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__ ( self             ,
                   H1               ,   ## H1 model, e.g. background only
                   H0               ,   ## H0 model, e.g. signal+background
                   dataset          ) : ## dataset 
        
        assert isinstance ( H1 , ( ROOT.RooStats.ModelConfig , ModelConfig ) ) , \
               'Invalid type for H1!'
        assert isinstance ( H0 , ( ROOT.RooStats.ModelConfig , ModelConfig ) ) , \
               'Invalid type for H0!'

        self.__H1         = H1
        self.__H0         = H0
        self.__dataset    = dataset 
        self.__calculator = None
        
    ## abstract method to create (and configire) calculator
    @abc.abstractmethod 
    def make_calculator ( self ) :
        """Abstract method to create the calculator"""
        pass

    @property
    def calculator ( self ) :
        """'calculator' : actual RooStats calculator"""
        if not self.__calculator : self.__calculator = self.make_calculator()
        return self.__calculator
    
    @property
    def H0 ( self ) :
        """'H0' - H0-model"""
        return self.__H0 
    @property
    def H1 ( self ) :
        """'H1' - H1-model"""
        return self.__H1 
    @property
    def h0 ( self ) :
        """'h0' : H0-model as `ROOT.RooStats.ModelConfig`"""
        return self.H0 if isinstance ( self.H0 , ROOT.RooStats.ModelConfig ) else self.H0.mc
    @property
    def h1 ( self ) :
        """'h1' : H1-model as `ROOT.RooStats.ModelConfig`"""
        return self.H1 if isinstance ( self.H1 , ROOT.RooStats.ModelConfig ) else self.H1.mc
    @property
    def dataset ( self ) :
        """'dataset' : dataset used for calculations"""
        return self.__dataset
    
# =============================================================================
## @class AsymptoticCalculator
#  @see RooStats::AsymptoticCalculator
class AsymptoticCalculator (Calculator) :
    """Asymptotoc calcualator
    - see `ROOT.RooStats.AsymptoticCalculator`
    """
    
    def __init__ ( self             ,
                   H1               ,   ## H1 model, e.g. background only
                   H0               ,   ## H0 model, e.g. signal+background
                   dataset          ,   ## dataset 
                   one_sided = True ) : 

        self.__one_sided  = one_sided 
        Calculator.__init__ ( self , H1 , H0 , dataset )

    ## create and configire the calculator 
    def make_calculator ( self ) :
        calc = ROOT.RooStats.AsymptoticCalculator ( self.dataset , self.h1 , self.h0 ) 
        calc.SetOneSided ( self.one_sided )
        return calc
        
    @property
    def one_sided ( self ) :
        """'one_sided' : actual RooStats calculator"""
        return self.__one_sided

# =============================================================================
## @class FrequentistCalculator
#  @see RooStats::FreqentistCalculator
class FrequentistCalculator (Calculator) :
    """Frequentist calcualator
    - see `ROOT.RooStats.FrequentistCalculator`
    """
    
    def __init__ ( self                   ,
                   H1                     ,   ## H1 model, e.g. background only
                   H0                     ,   ## H0 model, e.g. signal+background
                   dataset                ,   ## dataset 
                   sampler         = None ,
                   ntoys_null      = -1   ,
                   ntoys_alt       = -1   ,
                   ntoys_null_tail =  0   ,
                   ntoys_alt_tail  =  0   ,
                   ) : 

        assert isinstance ( ntoys_null , integer_types ) and -1 <= ntoys_null , \
               "Invalid ntoys_null parameter!"
        assert isinstance ( ntoys_alt  , integer_types ) and -1 <= ntoys_alt  , \
               "Invalid ntoys_alt parameter!"
        assert isinstance ( ntoys_null_tail , integer_types ) and 0 <= ntoys_null_tail , \
               "Invalid ntoys_null_tail parameter!"
        assert isinstance ( ntoys_alt_tail  , integer_types ) and 0 <= ntoys_alt_tail  , \
               "Invalid ntoys_alt_tail parameter!"
        
        self.__ntoys_null      = ntoys_null 
        self.__ntoys_alt       = ntoys_alt 
        self.__ntoys_null_tail = ntoys_null_tail
        self.__ntoys_alt_tail  = ntoys_alt_tail
        
        assert sampler is None or ( sampler and isinstance ( sampler , ROOT.RooStat.TestStatSampler ) ) , \
               'Invalid sampler!'
        
        slef.__sampler = sampler 
    
        Calculator.__init__ ( self , H1 , H0 , dataset )

    # ==============================================================================================
    ## create and configire the calculator 
    def make_calculator ( self ) :
        """Create and configure the calculator"""
        if self.sampler : calc = ROOT.RooStats.FrequestistCalculator ( self.dataset , self.h1 , self.h0 , self.sampler )
        else            : calc = ROOT.RooStats.FrequestistCalculator ( self.dataset , self.h1 , self.h0 )
        
        if -1 != self.ntoys_null      or -1 != slef.ntoys_alt :
            calc.SetToys         ( self.ntoys_null      , self.ntoys_alt      )
        if  0 != self.ntoys_null_tail or -1 != slef.ntoys_alt_tail  :
            calc.SetNToysInTails ( self.ntoys_null_tail , self.ntoys_alt_tail )
            
        return calc
        
    @property
    def sampler ( self ) :
        """'sampler' : sampler used for toys"""
        return self.__sampler

    @property
    def ntoys_null ( self ) :
        """'ntoys_null' : the first parameter of `ROOT.RooStats.FrequentistCalculator.SetToys`"""
        return self.__ntoys_null
    @property
    def ntoys_alt  ( self ) :
        """'ntoys_alt' : the seconf parameter of `ROOT.RooStats.FrequentistCalculator.SetToys`"""
        return self.__ntoys_alt 
    @property
    def ntoys_null_tail ( self ) :
        """'ntoys_null_tail' : the first parameter of `ROOT.RooStats.FrequentistCalculator.SetNToysInTails`"""
        return self.__ntoys_null_tail
    @property
    def ntoys_alt_tail  ( self ) :
        """'ntoys_alt_tail' : the second parameter of `ROOT.RooStats.FrequentistCalculator.SetNToysInTails`"""
        return self.__ntoys_alt_tail 
    

# =============================================================================
# @class HypoTestInverter
# @see RooStats::HypoTestInverter
class HypoTestInverter(object) :
    """Hypo test inverter 
    -see `ROOT.RooStats.HypoTestInverter`
    """
    
    def __init__ ( self            ,
                   calculator      ,
                   level           ,
                   use_CLs         ,
                   verbose = False ) :
        
        self.__calculator = calculator
        self.__inverter   = ROOT.RooStats.HypoTestInverter ( self.calc )
        
        self.__inverter.SetConfidenceLevel   ( level )
        
        if use_CLs :  self.__inverter.UseCLs ( True )
        self.__inverter.SetVerbose ( verbose )
        self.__interval = None 
        self.__plot     = None 
        
    @property
    def calculator ( self ) :
        """'calculator' : calcualtor"""
        return self.__calculator
    
    @property
    def calc ( self ) :
        """'calc' : calculator as RooStats object"""
        c = self.__calculator 
        return  c if not isinstance ( c , Calculator ) else c.calculator 

    # =========================================================================
    ## define (or perform the actual scan)
    #  @code
    #  hti = ...
    #  hti.scan ( 100 , 0.0 , 10.0 )   ## define scan with 100 point between 0 and 10
    #  hti.scan ()                     ## define auto scan 
    #  hti.scan (  [ 0, 1, 2, 3, 4 ] ) ## define custom scane    
    #  @endcode
    def scan ( self , *values ) :
        """Define (or perform the actual scan)
        >>> hti = ...
        >>> hti.scan ( 100 , 0.0 , 10.0 )   ## define scan with 100 point between 0 and 10
        >>> hti.scan ()                     ## define auto scan 
        >>> hti.scan (  [ 0, 1, 2, 3, 4 ] ) ## define custom scane    
        """
        
        if    not values :
            self.__inverter.SetAutoScan  ()
            if self.__interval : self.__interval = None
            if self.__plot     : self.__plot = None 
        elif  3 == len ( values ) and isinstance ( values [ 0 ] , integer_types ) and 2 <= values [ 0 ] :
            self.__inverter.SetFixedScan ( *values )
            if self.__interval : self.__interval = None 
            if self.__plot     : self.__plot = None 
        elif  1 == len ( values ) and isinstance ( values [ 0 ] , sequence_types  ) :
            for v in values [ 0 ] : self.__inverter.RunOnePoint ( v )
        else :
            for v in values : self.__inverter.RunOnePoint ( v )
            
    @property
    def inverter ( self ) :
        """'inverter' : actual `HypoTestInverter` object from `RooStats`"""
        return self.__inverter
    
    @property
    def interval ( self ) :
        """'interval' : get the confidence iterval"""
        if not self.__interval : self.__interval =  self.__inverter.GetInterval()
        return self.__interval
    
    @property 
    def limits ( self ) :
        """'limits' : get the upper and lower limit"""
        ii = self.interval
        return ii.UpperLimit() , ii.UpperLimit() 

    @property 
    def upper_limit ( self ) :
        """'upper_limit' : get the upper limit"""
        ii = self.interval
        return ii.UpperLimit() 

    @property 
    def lower_limit ( self ) :
        """'lower_limit' : get the upper limit"""
        ii = self.interval
        return ii.LowerLimit() 

    @property
    def plot ( self ) :
        """'plot' : prepare the plot
        
        >>> hti  = HypoTestInverter( ... )
        >>> plot = hti.plot
        >>> plot.Draw('CLb 2CL')
        
        Possible options:
        
        - SAME : draw in the current axis
        - OBS  : draw only the observed plot
        - EXP  : draw only the expected plot
        - CLB  : draw also the CLB
        - 2CL  : draw both clsplusb and cls
        
        default draw observed + expected with 1 and 2 sigma bands
        """
        if not self.__plot :
            self.__plot = ROOT.RooStats.HypoTestInverterPlot('','', self.interval  )
        return self.__plot 
    
    # ========================================================================= 
    ## Reset/clear current interval & plot
    def reset  ( self ) :
        """Reset/clear current interval & plot
        """
        if self.__interval : self.__interval = None
        if self.__plot     : self.__plot = None 
        
# =============================================================================
## default color for Brasil-plot bands
band_colors = ( ROOT.kGreen , ROOT.kYellow  ,
                ROOT.kCyan  , ROOT.kMagenta ,
                ROOT.kBlue  , ROOT.kOrange  )
# =============================================================================
## helper class to create and keep the 'Brasil-band' plot
#  @code
#  bp = BrasilBand( sigmas = (1,2,3) )  ## draw 1,2&3-sigma bands 
#  for value  in [ ... ] :
#       ...
#       hti      = HypoTestInverter ( ... )
#       interval = hti.interval
#       ul       = hti.upper_limit
#       bp.fill ( value , limit ,  interval )
#  p = bp.plot()
#  p.draw('a')
#  @endcode 
class BrasilBand(object) :
    """Helper class to create and keep the 'Brasil-band' plot
    >>> bp = BrasilBand( sigmas = (1,2,3) )  ## draw 1,2&3-sigma bands 
    >>> for value  in [ ... ] :
    >>>     ...
    >>>     hti      = HypoTestInverter ( ... )
    >>      interval = hti.interval
    >>>     limit    = hti.upper_limit
    >>>     bp.fill ( value , limit ,  interval )
    >>> p = bp.plot()
    >>> p.draw('a')
    """
    def __init__  ( self               ,
                    sigmas = ( 1 , 2 ) ,
                    colors = ()        ) :

        ss = set()
        for s in sigmas : ss.add (  1*s ) 
        for s in sigmas : ss.add ( -1*s ) 
        if  ( not 0 in ss ) and ( not 0.0 in ss ) : ss.add ( 0 )
        
        self.__nsigmas = tuple ( sorted ( ss ) ) 
        self.__data    = {} 

        assert 1 == len ( self.__nsigmas ) % 2 , 'Invalid number of nsigmas!'

        self.__colors = [] 
        for i , j in zip_longest ( colors , band_colors , fillvalue = None ) :
            if   not i is None : self.__colors.append ( i )
            elif not j is None : self.__colors.append ( j )

        nb = len ( self.__nsigmas ) // 2
        ic = 31 
        while len ( self.__colors ) < nb :
            self.__solors.append ( ic  ) 
            ic += 1

        self.__colors = tuple ( self.__colors )
        
        self.__plot   = None
        self.__legend = None
        
    # =========================================================================
    ## Add the point to the Brasil-plot
    def fill ( self , x , observed , hti_result ) : 
        """Add the point to the Brasil-plot
        """
        expected = tuple ( hti_result.GetExpectedUpperLimit ( s ) for s in self.__nsigmas )
        self.__data [ x ] = observed, expected
        if self.__plot   : self.__plot   = None 
        if self.__legend : self.__legend = None 

    ## ========================================================================
    #  Create the actual (multi)-graph
    #  @see TMultGraph
    #  @code 
    #  bp = BrasilBand ( ... )
    #  for ... :
    #      bp.fill ( ... ) 
    #  bp.plot.draw('a')
    #  bp.legend.draw('a')
    @property 
    def plot  ( self ) :
        """'plot' : get the actual (multi)-graph
        -see `ROOT.TMultGraph` 
        >>> bp = BrasilBand ( ... )
        >>> for .... :
        >>>    bp.fill ( ... ) 
        >>> bp.plot.draw('a')
        >>> bp.legend.draw('a')
        """

        if self.__plot and self.__legend : return self.__plot

        import ostap.histos.graphs

        np  = len ( self.__data ) 
        ## bands
        ns  = len ( self.__nsigmas )
        ngb = ns // 2

        ## get sorted data 
        data = sorted ( self.__data.items() )
        
        ## create "band-graphs"
        gr_observed = ROOT.TGraph ( np )
        gr_median   = ROOT.TGraph ( np ) 
        gr_bands    = [ ROOT.TGraphAsymmErrors ( np ) for i in range ( ngb ) ]

        for i , x in  enumerate ( sorted ( self.__data ) ) :

            observed , expected = self.__data [ x ] 

            median   = expected [ ngb ] 
            
            gr_observed [ i ] = x , observed
            gr_median   [ i ] = x , median  

            for j , g in enumerate ( gr_bands ) :

                el   = expected [          j ]
                eh   = expected [ ns - 1 - j ]
                val  = median if el < median < eh else 0.5 * ( eh + el )

                errl = el - val
                errh = eh - val 
                
                g [ i ] = ( x , 0 , 0 ) , ( val, errl , errh )


        for g,c  in zip ( reversed ( gr_bands ) , self.__colors ) : g.color ( c , marker = 1 , size = 1 )

        gr_median   .SetLineStyle ( 9 )
        gr_median   .SetLineWidth ( 2 )
        
        gr_observed .red          (   ) 
        
        ## create the final (multi) graph & populate it 
        result = ROOT.TMultiGraph()
        for g in gr_bands :
            g.SetFillStyle ( 1001 ) 
            result.Add ( g , '3' ) 
            
        result.Add ( gr_median   , 'L'  ) 
        result.Add ( gr_observed , 'LP' )

        self.__plot = result 

        ## re-create the legend 
        legend = ROOT.TLegend ( 0.2 , 0.65 , 0.4 , 0.9 )
        legend.SetFillColor  ( 0 )
        legend.SetFillStyle  ( 0 )
        legend.SetBorderSize ( 0 )
        
        legend.AddEntry ( gr_observed , "observed", "lp" )
        legend.AddEntry ( gr_median   , "expected", "l"  )

        for g , ns in zip ( reversed ( gr_bands ) , self.__nsigmas[ngb+1:]  ) :
            legend.AddEntry( g ,"expected #pm%s#sigma"% abs ( ns ) , "f")
            
        self.__legend = legend 

        return self.__plot

    # =========================================================================
    ## get standard legend for the plot
    #  @code 
    #  bp = BrasilBand ( ... ) 
    #  for .... :
    #    bp.fill ( ... ) 
    #  bp.plot.draw('a')
    #  bp.legend.draw('a')
    #  @endcode 
    @property
    def legend ( self )  :
        """'legend' : get standarde legend for the plot
        >>> bp = BrasilBand ( ... ) 
        >>> for .... :
        >>>    bp.fill ( ... ) 
        >>> bp.plot.draw('a')
        >>> bp.legend.draw('a')
        """
        if not self.__legend : plot = self.plot
        return self.__legend
    
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
