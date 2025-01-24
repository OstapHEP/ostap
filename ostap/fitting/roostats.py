#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/roostats.py
#  Set of useful basic utilities to dela with RooStats
#
# - confidence interval estimation
# - confidence intervals and limits
# - hypothesis tests
# - Brasilian  plots 
#
#  @author Vanya BELYAEV   Ivan.Belyaev@itep.ru
#  @author Dmitry GOLUBKOV Dmitry.Yu.Golubkov@cern.ch
#  @date 2023-01-17
# =============================================================================
""" Set of useful basic utilities to deal with RooStats

- confidence interval estimation
- confidence interval and limits
- hypothesis tests
- Brasilian  plots 

"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru; Dmitry GOLUBKOV Dmitry.Yu.Golubkov@cern.ch"
__date__    = "2023-01-17"
__all__     = (
    #
    ## confidence intervals and limits
    # 
    'ModelConfig'                 , ## creator of RooStats ModelConfig 
    'ProfileLikelihoodInterval'   , ## Profile-likelihood confidence interval or upper/lower limit 
    'FeldmanCousinsInterval'      , ## Feldman-Cousins    confidence interval or upper/lower limit 
    'BayesianInterval'            , ## Bayesian           confidence interval or upper/lowee limir
    'MCMCInterval'                , ## MCMC               confidence interval or upper/lower limit 
    #
    ## Hypothesis tests
    # 
    'AsymptoticCalculator'        , ## AsymptoticCalculator          for limits and intervals
    'FrequentistCalculator'       , ## Frequentist calcualtor        for limits and intervals 
    'HybridCalculator'            , ## Hybrid calculator             for limits and intervals 
    'ProfileLikelihoodCalculator' , ## Profile Likelihood calculator for limits and intervals 
    'HypoTestInverter'            , ## 
    ##
    'BrasilBand'                  , ## utility to produce Brasil-band plots 
    'P0Plot'                      , ## utility to produce p0-plots 
    )
# ==============================================================================
from   ostap.core.meta_info   import root_info 
from   ostap.core.ostap_types import ( string_types    , integer_types  ,
                                        num_types      , 
                                        sequence_types , dictlike_types )
from   ostap.core.core        import ( valid_pointer   , split_string   ,
                                       loop_items      , Ostap          ,
                                       typename  )    
import ostap.fitting.roofit
from   ostap.fitting.pdfbasic import APDF1
import ROOT, abc, sys  
# ==============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.roostats' )
else                       : logger = getLogger ( __name__                 )
# ==============================================================================
if (3,0) <= sys.version_info : from itertools import  zip_longest
else                         : from itertools import izip_longest as zip_longest 
# ==============================================================================
_new_methods_ = []
# ==============================================================================
## allowed types for observables 
obs_types = ROOT.RooAbsReal , ROOT.RooAbsCategory
## allowed types for parameters 
par_types = ROOT.RooAbsReal ,
# # ============================================================================
pdf_like  = lambda o : isinstance ( o , ROOT.RooAbsPdf ) or \
    ( hasattr ( o , 'roo_pdf' ) and isinstance ( o.roo_pdf , ROOT.RooAbsPdf ) )
get_pdf   = lambda c : c if isinstance  ( c , ROOT.RooAbsPdf ) else c.roo_pdf
# ==============================================================================
if ( 6 , 28 ) <= root_info : # =================================================
    # ==========================================================================
    def make_set ( source ) : return ROOT.RooArgSet ( source )
    # ==========================================================================
else : # =======================================================================
    # ==========================================================================
    def make_set ( source ) : 
        result = ROOT.RooArgSet()
        for item in source : result.add ( item )
        return result
# ==============================================================================
## @class ModelConfig
#  Helper class to create `RooStats::ModelConfig`
#  @see RooStats::ModelConfig 
class ModelConfig(object): 
    """ Helper class to create `RooStats::ModelConfig`
    - see `ROOT.RooStats.ModelConfig`
    """
    def __init__  ( self                        ,
                    pdf                         ,   ## PDF 
                    poi                         ,   ## parameter(s) of interest
                    dataset                     ,   ## dataset  
                    workspace            = None ,   ## worspace (optional)
                    observables          = ()   ,   ## observables 
                    global_observables   = ()   ,   ## global observables
                    constraints          = ()   ,   ## constraints
                    external_constraints = ()   ,   ## external constraints 
                    constrained          = ()   ,   ## constrained
                    conditional          = ()   ,   ## conditional observables 
                    snapshot             = None ,   ## snapshot 
                    **kwargs                    ) : ## other arguments

        ## (0) allow some freedom in arguments 
        from ostap.utils.cidict import cidict
        from ostap.core.core    import cidict_fun
        ## case-insensitive kewword argumetns 
        kw_args = cidict ( transform = cidict_fun , **kwargs )
        
        ## check pdf 
        assert isinstance ( pdf , ROOT.RooAbsPdf ) or \
            ( hasattr ( pdf , 'roo_pdf' ) and isinstance ( pdf.roo_pdf , ROOT.RooAbsPdf ) ) , \
            "Invalid `pdf' type:%s" % typename ( pdf )
        
        ## check dataset 
        assert dataset and isinstance ( dataset , ROOT.RooAbsData ), \
            "Invalid dataset type:%s" % typename ( dataset ) 

        ## check workspace 
        assert ( not workspace ) or isinstance ( workspace , ROOT.RooWorkspace ) , \
            "Invalid worskpace type:%s" % typename ( workspace )

        ## (1) quick check and very simple transform. more checks follow  

        if   isinstance ( observables , string_types    ) :
            observables = split_string ( observables  , strip = True , respect_groups = True )
        elif isinstance ( observables , obs_types       ) : 
            observables = [ observables ]
        
        if   isinstance ( poi         , string_types    ) :
            poi = split_string ( poi , strip = True , respect_groups = True )
        elif isinstance ( poi         , par_types       ) : poi         = [ poi         ]
        
        if   isinstance ( constrained , string_types    ) : 
            constrained = split_string ( constrained , strip = True , respect_groups = True )        
        elif isinstance ( constrained , ROOT.RooAbsReal ) : constrained = [ constrained ]

        if   isinstance ( conditional , string_types  ) :
            conditional = split_string ( conditional , strip = True , respect_groups = True )
        elif isinstance ( conditional , ROOT.RooAbsReal ) : conditional = [ conditional ]

        if   pdf_like ( constraints          ) : constraints          = [ constraints          ]
        if   pdf_like ( external_constraints ) : external_constraints = [ external_constraints ]
        
        if   isinstance ( global_observables , string_types    ) :
            global_observables = split_string ( global_observables  , strip = True , respect_groups = True )
        elif isinstance ( global_observables , obs_types       ) : global_observables = [ global_observables ]

        
        ## save all input args
        self.__input_pdf                  = pdf
        self.__input_poi                  = poi 
        self.__input_dataset              = dataset 
        self.__input_observables          = observables
        self.__input_global_observables   = global_observables
        self.__input_constraints          = constraints
        self.__input_external_constraints = external_constraints        
        self.__input_constrained          = constrained 
        self.__input_conditional          = conditional
        self.__input_shapshot             = snapshot
        self.__input_kwargs               = kwargs 

        self.__keep = [] 

        ## (1) create workspace (if needed)
        if not workspace :
            wsname    = kw_args.pop ( 'ws_name'  , 'WS_%s'            % pdf.name )
            wstitle   = kw_args.pop ( 'ws_title' , 'workspace for %s' % pdf.name )
            workspace = ROOT.RooWorkspace         ( wsname  , wstitle )            
        self.__ws   = workspace
    
        ## (2) create ModelConfig        
        mcname    = kw_args.pop ( 'name'  , 'MC_%s'               % pdf.name )
        mctitle   = kw_args.pop ( 'title' , 'model-config for %s' % pdf.name )
        mc        = ROOT.RooStats.ModelConfig ( mcname , mctitle , self.ws  )
        self.__mc = mc

        ## (3) raw-pwd 
        raw_pdf   = pdf if isinstance ( pdf , ROOT.RooAbsPdf ) else pdf.roo_pdf

        final_pdf = raw_pdf

        ## (4) add constraints if specified 
        self.__constraints = () 
        if constraints :
            assert all ( pdf_like ( c ) for c in constraints ) , "Invalid constraint type!"                
            constraints = tuple ( get_pdf ( c )  for c in constraints )
            from ostap.fitting.constraint import make_constrained            
            final_pdf, tail = make_constrained ( raw_pdf , *constraints )
            if tail : self.__keep .append ( tail ) 
            self.__constraints = constraints
            
        self.__raw_pdf   =   raw_pdf 
        self.__final_pdf = final_pdf 
        
        ## propagate PDF to ModelConfig/Workspace 
        ## (6) set PDF  
        mc.SetPdf ( self.__final_pdf ) ## note: we use bare RooAbsPdf here

        if not observables :
            if   hasattr ( pdf , 'observables' ) : observables = pdf.observables
            else                                 : observables = final_pdf.getObservables ( dataset )

        ## if data&observables are specified, reuqire all observables in data 
        assert all ( v in dataset for v in observables ) , \
            "Specified observables are not consistent with data"

        assert all ( v in final_pdf.getObservables ( dataset ) for v in observables ) , \
            "Specified observables are not consistent with PDF-observables" 

        ## PDF for sure depends on observables 
        assert all ( v in final_pdf for v in observables ), \
            "Specified observables are not consistent with PDF"
        
        ## (7) external constraints, just for bookkeeping   
        self.__external_constraints = () 
        if external_constraints : 
            assert all ( pdf_like ( c ) for c in external_constraints ) , "Invalid external constraint type!"
            external_constraints = tuple ( get_pdf ( c ) for c in external_constraints )
            ## Inform ModelConfig on the addtitional constraints
            if ( 6 , 28 , 10 ) <= root_info :
                self.__mc.SetExternalConstraints ( ROOT.RooArgSet ( external_constraints ) ) 
            self.__external_constraints = external_constraints 
            

        observables = make_set ( self.pdf_observable ( o , dataset ) for o in observables ) 
            
        ##  (8) set observables 
        self.__mc.SetObservables ( observables ) 
        self.__observables = observables

        ##  (9) get parameters of interess
        poi = make_set ( self.pdf_param ( p , dataset ) for p in poi ) 
        mc.SetParametersOfInterest ( poi )
        self.__poi  = poi 

        ##  (10) global observables
        assert all ( self.pdf_param ( p , dataset ) for p in global_observables )
        if global_observables :
            global_observables = make_set ( v for v in global_observables )                
            ## assert all ( any ( o in c for c in constraints ) for o in global_observables ) , \
            ##    "Global observables are inconsisent with external constraints!"
            mc.SetGlobalObservables  ( global_observables ) 
        self.__global_observables = global_observables 

        ## nuisancee = all_parameters - poi - global_observables 
        nuisance = final_pdf.getParameters ( dataset ) - poi - global_observables
        if nuisance : ROOT.RooStats.RemoveConstantParameters ( nuisance )
        
        ## (11) Nuisance parameters        
        mc.SetNuisanceParameters ( nuisance )
        self.__nuisance = nuisance 
        
        ## (12) conditional observables
        assert all ( self.pdF_param ( p , dataset ) for p in conditional ) 
        if conditional :
            conditional = ROOT.RooArgSet ( v for v in conditional )
            mc.SetConditionalObservables  ( conditional ) 
        self.__conditional = conditional

        ## (13) constrained parameters (not used by RooStats) 
        assert all ( self.pdf_param ( p , dataset ) for p in constrained  )
        if constrained :
            constrained = make_set ( v for v in constrained )
            mc.SetConstraintParameters ( constrained ) 
        self.__constrained = constrained  

        ## (14) define the default dataset 
        self.__dataset = dataset 

        ## (15) use snapshot if provided 
        if snapshot : self.snapshot = snapshot
            
        ## the rest ? 
        if kw_args :
            import ostap.logger.table as T 
            rows = [ ( 'Argument' , 'Value' ) ]
            for k , v in loop_items ( kw_args ) :
                row = k , str ( v )
                rows.append ( row )
            title = 'ModelConfig: %d ignored arguments' % len ( kw_args ) 
            table = T.table ( rows , title = title , prefix = '# ' , alignment = 'll' )    
            logger.warning ( '%s\n%s' % ( title , table ) )

        ## store in Workspace
        self.ws.Import ( self.mc )
        
        ## finally print as table 
        logger.debug ( 'Created ModelConfig: %s\n%s' % ( self.name , self.table ( prefix = '# ' ) ) ) 
                       
    @property
    def name  ( self ) :
        """'name': the name of the `RooStats.ModelConfig` objxct"""
        return self.mc.name
    @property
    def title ( self ) :
        """'title': the title of the `RooStats.ModelConfig` objxct"""
        return self.mc.title
        
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
    def pdf ( self ) :
        """`pdf` : PDF from ModelConfig"""
        return self.mc.GetPdf()

    @property
    def constraints  ( self ) :
        """'constraints' : constrain parameters from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetConstraintParameters`
        """
        return self.__constraints
    
    @property
    def observables ( self ) :
        """'observables' : Global observables from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetObservables`
        """
        pars = self.mc.GetObservables()
        if not valid_pointer ( pars ) : return () 
        return pars if pars and 0 < len ( pars ) else () 
    @property
    def poi ( self ) :
        """'poi' : parameter(s) of interest from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetParametersOfInterest`
        """
        pars = self.mc.GetParametersOfInterest()
        if not valid_pointer ( pars ) : return () 
        return pars if pars and 0 < len ( pars ) else () 
    @property
    def nuisance ( self ) :
        """'nuisance' : nuisance parameters from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetNuisanceParameters`
        """
        pars = self.mc.GetNuisanceParameters()
        if not valid_pointer ( pars ) : return () 
        return pars if pars and 0 < len ( pars ) else () 
    @property
    def global_observables ( self ) :
        """'global_observables' : Global observables from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetGlobalObservables`
        """
        pars = self.mc.GetGlobalObservables()
        if not valid_pointer ( pars ) : return () 
        return pars if pars and 0 < len ( pars ) else () 
    @property
    def conditional ( self ) :
        """'conditional' : Conditional observables from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetConditionalObservables`
        """
        pars = self.mc.GetConditionalObservables()
        if not valid_pointer ( pars ) : return () 
        return pars if pars and 0 < len ( pars ) else ()
    @property 
    def constrained ( self ) :
        """'constrained' : constrained parameters from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetConstraintParameters 
        """
        pars = self.mc.GetConstraintParameters () 
        if not valid_pointer ( pars ) : return () 
        return pars if pars and 0 < len ( pars ) else () 
    
    @property
    def external_constraints  ( self ) :
        """'external_constraints' : constrain parameters from ModelConfig
        - see `ROOT.RooStats.ModelConfig.GetConstraintParameters`
        """
        if ( 6 , 28 , 10 ) <= root_info :
            pars = self.mc.GetExternalConstraints ()
            if not valid_pointer ( pars ) : return () 
            return pars if pars and 0 < len ( pars ) else () 
        return self.__external_constraints

    # =========================================================================
    ## Get/set snapshot from the model/workspacee
    #  @code 
    #  model = ...
    #  s     = model.snapshot
    # 
    #   model.snapshot = var                         ## variable
    #   model.snapshot = var1,var2,var3              ## sequence of variables 
    #   model.snapshot = ROOT.RooArgSet (... )
    #   model.snapshot = fit_result                  ## ROOT.RooFitResult
    #   model.snapshot = { 'var1' : 1 , 'var2' : 2 } ## dictionary 
    #  @see RooStats::ModelConfig::GetSnapshot
    #  @see RooStats::ModelConfig::SetSnapshot
    @property
    def snapshot ( self ) :
        """'snapshot' : get/set snapshot from  the model/workspacee
        >>> model = ...
        >>> s     = model.snapshot

        >>> model.snapshot = var                         ## variable
        >>> model.snapshot = var1,var2,var3              ## sequence of variables 
        >>> model.snapshot = ROOT.RooArgSet (... )
        >>> model.snapshot = fit_result                  ## ROOT.RooFitResult
        >>> model.snapshot = { 'var1' : 1 , 'var2' : 2 } ## dictionary 
        
        - see `ROOT.RooStats.ModelConfig.GetSnapshot`
        - see `ROOT.RooStats.ModelConfig.SetSnapshot`        
        """
        pars = self.mc.GetSnapshot()
        if not valid_pointer ( pars ) : return () 
        return pars if pars and 0 < len ( pars ) else ()
    @snapshot.setter
    def snapshot ( self , values ) :
        
        if   isinstance ( values , ROOT.RooAbsReal   ) : return self.mc.SetSnapshot ( ROOT.RooArgSet ( values ) )
        elif isinstance ( values , ROOT.RooArgSet    ) : return self.mc.SetSnapshot (  values )
        elif isinstance ( values , ROOT.RooFitResult ) :            
            vs = ROOT.RooArgSet()
            for p in values.floatParsFinal () : vs.add ( p )
            for p in values.constPars      () : vs.add ( p )
            return self.mc.SetSnapshot ( vs  )
        
        elif isinstance ( values , dictlike_types ) :

            vlst = []
            vdct = {}            
            for key in values :
                v = self.var ( key )                
                assert v , "No valid variable for '%s'" % key                
                vv = float ( values [ key ] )
                vlst.append ( v )
                vdct [ v.name ]  = vv

            with SETVAL ( *vlst ) :
                vset = ROOT.RooArgSet() 
                for v in vlst : v.setVal ( vdct [ v.name ] )
                for v in vlst : vset.add ( v )
                return self.mc.SetSnapshot ( vset  )
                
        elif isinstance ( values , sequence_types  ) :

            vs = ROOT.RooArgSet() 
            for vv in values :
                if   isinstance ( vv , ROOT.RooAbsReal       ) : vs.add ( vv )
                elif isinstance ( vv , ROOT.RooRooFitResult  ) :
                    for p in vv.floatParsFinal () : vs.add ( p )
                    for p in vv.constPars      () : vs.add ( p )
                elif isinstance ( v , ROOT.RooAbsCollection ) :
                    for v in cv : vs.add ( v )
                else :
                    raise TypeError ( 'Invalid type for snapshot %s' % type ( vv )  )
                
            return self.mc.SetSnapshot ( vs  )
        
        raise TypeError ( 'Invalid type for snapshot %s' % type ( values )  ) 

    # =========================================================================
    ## helper function to get parameters from the final PDF 
    def pdf_params ( self , dataset = None ) :
        """ Helper function to get the parameters from the final PDF"""
        pdf = self.__final_pdf
        return pdf.getParameters ( ROOT.nullptr  ) if dataset is None else pdf.getParameters ( dataset )

    # =========================================================================
    ## helper function to get observables from the final PDF 
    def pdf_observables ( self , dataset = None ) :
        """ Helper function to get the parameters from the final PDF"""
        pdf = self.__final_pdf
        return pdf.getObservables( ROOT.nullptr  ) if dataset is None else pdf.getObservables ( dataset )
    
    # =========================================================================
    ## helper function to get parameters from PDF 
    def pdf_param  ( self , param , dataset = None ) :
        """ Helper function to get the parameters from PDF"""
        params = self.pdf_params ( dataset )
        ## already parameter ?        
        if   isinstance ( param , par_types ) and param in params : return param
        ## loop by name
        if isinstance ( param , string_types    ) :
            for p in params :
                if p.name == param       : return p                
        elif isinstance ( param , par_types ) :
            for p in params :
                if p.name == param.name  : return p
        else : raise TypeError ( "Inavlid parameter type %s/%s!" % ( param , typename ( param ) ) )
        
        ## try with workspace
        vv = self.var ( param )
        ## 
        if vv : return vv
        else  : raise KeyError ( "No parameter %s/%s defined" % ( param , typename ( param ) ) )

    # =========================================================================
    ## helper function to get observables from PDF 
    def pdf_observable  ( self , param , dataset = None ) :
        """ Helper function to get the parameters from PDF"""
        params = self.pdf_observables ( dataset )
        ## already parameter ?        
        if   isinstance ( param  , obs_types   ) and param in params : return param
        ## loop by name
        if isinstance   ( param , string_types ) :
            for p in params :
                if p.name == param       : return p                
        elif isinstance ( param , obs_types    ) :
            for p in params :
                if p.name == param.name  : return p
        else : raise TypeError ( "Inavlid observable type %s/%s!" % ( param , typename ( param ) ) )
        
        ## try with workspace
        vv = self.var ( param )
        ## 
        if vv : return vv
        else  : raise KeyError ( "No observable %s/%s defined" % ( param , typename ( param ) ) )
        
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
        """ Get a variable from workspace
        >>> mc = ...
        >>> v1 = mc.var ( 'mass' ) ## get by name
        >>> var = 
        >>> v1 = mc.var ( 'var   ) ## get by var     
        """
        if isinstance ( variable , par_types ) : return self.var ( variable.name )
        return self.ws.var ( variable ) 

    # ============================================================================
    ## print the model as a table 
    def table ( self , title = '' , prefix = '' ) : 
        """ Print the model as table
        """

        rows = [ ( '' , 'Model Config' , 'Value' ) ]

        row = '*' , 'Name' , self.mc.name
        rows.append ( row )
        
        row = '*' , 'Title' , self.mc.title 
        rows.append ( row )
        
        row  = '*' , 'Workspace' , self.ws.name
        rows.append ( row )
        row  = ''  , ''          , self.ws.title 
        rows.append ( row )

        row = '' , 'Input PDF    ' , "%s: %s" % ( typename ( self.__input_pdf ) , self.__input_pdf.name )
        rows.append ( row )
        
        if not self.__raw_pdf is self.__final_pdf :
            pdf = self.__raw_pdf 
            row = '' , 'Input raw PDF' , '%s: %s/%s' % ( typename  ( pdf ) , pdf.name , pdf.title )
            rows.append ( row )
            
        pdf = self.__final_pdf 
        row = '' , 'Final PDF' , '%s: %s/%s' % ( typename  ( pdf ) , pdf.name , pdf.title )                                       
        rows.append ( row )

        pdf = self.pdf
        row = '*' , 'PDF for RooStats' , '%s: %s/%s' % ( typename  ( pdf ) , pdf.name , pdf.title )                                       
        rows.append ( row )

        for i , c in enumerate ( self.constraints , start = 1 )  :
            row = '' , 'Added multiplicative constraint #%d:' % i , '%s: %s/%s' % ( typename ( c ) , c.name , c.title )
            rows.append ( row )

        row = '*' , 'Parameters-of-interest'  , ', '.join ( v.name for v in self.poi ) 
        rows.append ( row )

        row = '*' , 'Observables:'            , ', '.join ( v.name for v in self.observables )  
        rows.append ( row )

        row = '*' , 'Nuisance parameters'     , ', '.join ( v.name for v in self.nuisance )
        rows.append ( row )

        row = '*' , 'Global observables:'     , ', '.join ( v.name for v in self.global_observables ) 
        rows.append ( row )

        row = '*' , 'Conditional parameters:' , ', '.join ( v.name for v in self.conditional ) 
        rows.append ( row )
        
        row = '*' , 'Constrained parameters:' , ', '.join ( v.name for v in self.constrained  ) 
        rows.append ( row )
        
        q = '*' if ( 6 , 28 , 10 ) <= root_info else '' 
        for i , c in enumerate ( self.external_constraints , start = 1 )  :
            row = q , 'External constraint #%d:' % i , '%s: %s/%s' % ( typename ( c ) , c.name , c.title )
            rows.append ( row )
                    
        vmax = -1
        for v in self.snapshot : vmax = max ( vmax , len ( v.name ) ) 
        vmax += 2
        vfmt = '%%%ds : %%-+7g' % vmax
        for i , v in enumerate ( self.snapshot ) :
            if 0 == i : row = '*' , 'Snapshot:' , vfmt % ( v.name , v.getVal () )
            else      : row = '*' ''            , vfmt % ( v.name , v.getVal () )
            rows.append ( row )
            
        if not self.snapshot : rows.append (  ( 'shapshot' , ) ) 
            
        import ostap.logger.table as T
        title = title if title else 'ModelConfig %s' % self.mc.title  
        return T.table ( rows, title = title , prefix = prefix , alignment = 'lw' )
    
# ================================================================================
## Helper (abstract) base class for the confidence intervals and limits 
class CLInterval(ModelConfig)  :
    """ Helper (abstract) base class for confidence intervals and limits"""
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
        """ Abstract method to create the interval calculator"""
        pass

    # =========================================================================
    ## get the confidence interval at certain confidence level
    #  @code
    #  interval = ...
    #  lower, upper = interval.interval ( 0.90 ) 
    #  @endcode 
    def interval ( self , level , par = None , dataset = None ) :
        """ Get the confidence interal at certain confidence level
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
        """ Get the upper/lower limit at certain confidence level
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
        """ Get the upper limit at certain confidence level
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
        """ Get the lower limit at certain confidence level
        >>> interval =...
        >>> lower    = interval.lower_limit ( 0.90 ) 
        """        
        return self.limit ( level , limit = -1 , par = par , dataset = dataset )

    # =========================================================================
    ## Helper method to get the true parameter from poi 
    def par_from_poi ( self , par ) :
        """ Helper method to get the true parameter from poi
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
    """ Profile Likelihood confidence interval
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
        """ Create the interval"""

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
        """ Make a plot
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
    """ Feldman-Cousins confidence interval
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
        """ Create the interval"""
        
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
        """ Visualize the interval
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
            for entry, _  in ps :                
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
    """ Bayesian confidence interval
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
        """ Create the interval"""
        
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
        """ Make a plot
        - see `ROOT.RooStats.BayesianCalculator.GetPosteriorPlot`
        """
        if self.calculator :
            return self.calculator.GetPosteriorPlot() 

    # =========================================================================
    ## Helper method to get the true parameter from poi 
    def par_from_poi ( self , par ) :
        """ Helper method to get the true parameter from poi
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
    """ Markov Chain confidence interval
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
        """ Create the interval"""
        
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
        """ Make a plot
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
    """ Base class for Calculators 
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
        """ Abstract method to create the calculator"""
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

    @property
    def hypo_test ( self )  :
        """'hypo_test' : get `HypoTestResult` from calculator"""
        return self.calculator.GetHypoTest() 

# =============================================================================
## @class AsymptoticCalculator
#  @see RooStats::AsymptoticCalculator
class AsymptoticCalculator (Calculator) :
    """ Asymptotic calcualator
    - see `ROOT.RooStats.AsymptoticCalculator`
    """
    
    def __init__ ( self              ,
                   H1                ,   ## H1 model, e.g. background only
                   H0                ,   ## H0 model, e.g. signal+background
                   dataset           ,   ## dataset 
                   asimov    = False ,   ## nominal Asimov? (not using fitted parameter values but nominal ones)
                   silent    = True  , 
                   verbose   = False ) : 

        self.__asimov     = True if asimov else False
        ## 
        if   silent  : self.__level = -1
        elif verbose : self.__level =  2
        else         : self.__level =  1
        ##
        Calculator.__init__ ( self , H1 , H0 , dataset )

    ## create and configure the Asymptotic calculator 
    def make_calculator ( self ) :
        """Create and configure the Asymptotic calculator"""
        ROOT.RooStats.AsymptoticCalculator.SetPrintLevel ( self.__level ) 
        calc = ROOT.RooStats.AsymptoticCalculator ( self.dataset , self.h1 , self.h0 , self.asimov ) 
        return calc
    
    @property
    def asimov    ( self ) :
        """'asimov' :  parameter for `ROOT.RooStats.AsymptoticCalculator(..., nominalAsimov = ...)`
        - (not using fitted parameter values but nominal ones)
        """
        return self.__asimov 

# =============================================================================
## @class FrequentistCalculator
#  @see RooStats::FreqentistCalculator
#  @warning FRequentist calculator corrupts  input data set! 
class FrequentistCalculator (Calculator) :
    """ Frequentist calcualator
    - see `ROOT.RooStats.FrequentistCalculator`
    - warning  Frequentist calculator corrupts input data set! 
    """
    
    def __init__ ( self                   ,
                   H1                     ,   ## H1 model, e.g. background only
                   H0                     ,   ## H0 model, e.g. signal+background
                   dataset                ,   ## dataset 
                   sampler         = None ,
                   ntoys_null      = -1   ,
                   ntoys_alt       = -1   ,
                   ntoys_null_tail =  0   ,
                   ntoys_alt_tail  =  0   ) : 
        
        assert sampler is None or ( sampler and isinstance ( sampler , ROOT.RooStat.TestStatSampler ) ) , \
               'Invalid sampler!'
        
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
        
        self.__sampler = sampler 

        Calculator.__init__ ( self , H1 , H0 , dataset )
        self.__cloned_dataset = dataset.clone()   

    ## Frequentist/Hybrid calculators can destroy their input dataset... Try to fix the problem 
    def __del__ ( self ) :

        if self.cloned_dataset :
            
            ds1 = self.dataset
            ds2 = self.__cloned_dataset
            
            vars1 = tuple ( v.name for v in ds1.varset () )
            
            s1 = ds1.statVars ( vars1 )
            s2 = ds2.statVars ( vars1 )

            from ostap.math.base import isequal
            
            diffs = [] 
            for k in s1 :
                ss1     = s1 [ k ]
                ss2     = s2 [ k ]
                mean1   = float ( ss1.mean() ) 
                mean2   = float ( ss2.mean() ) 
                rms1    = ss1.rms()
                rms2    = ss2.rms()
                mn1,mx1 = ss1.minmax()
                mn2,mx2 = ss2.minmax()
                if isequal ( mean1 , mean2 ) and \
                   isequal ( rms1  , rms2  ) and \
                   isequal ( mn1   , mn2   ) and \
                   isequal ( mx1   , mx2   ) :  pass
                else :
                    diffs.append ( k ) 

            if diffs :
                title = 'ORIGINAL dataset'
                logger.warning ("Frequentist calculator %s :\n%s" % ( title , ds1.table ( diffs , prefix = '# ' , title = title ) ) ) 
                title = 'CLONED   dataset'
                logger.warning ("Frequentist calculator %s :\n%s" % ( title , ds2.table ( diffs , prefix = '# ' , title = title ) ) ) 
            else :
                logger.info    ("Frequentist calculator: datasets are fine" )
                
            
            self.__cloned_dataset = None 
            if isinstance ( ds2 , ROOT.RooDataSet ) :
                ds2 = Ostap.MoreRooFit.delete_data ( ds2 )
                del ds2 

            logger.info ( 'Frequentist calculator: CLONED dataset is deleted' )

    @property
    def cloned_dataset ( self ) :
        """'cloned_dataset' : get the CLONED dataset (the one actually used in calculator)"""
        return self.__cloned_dataset
    
    # ==============================================================================================
    ## Create and configure the calculator 
    def make_calculator ( self ) :
        """ Create and configure the calculator"""
        ##
        if self.sampler : calc = ROOT.RooStats.FrequentistCalculator ( self.cloned_dataset , self.h1 , self.h0 , self.sampler )
        else            : calc = ROOT.RooStats.FrequentistCalculator ( self.cloned_dataset , self.h1 , self.h0 )
        
        if -1 != self.ntoys_null      or -1 != self.ntoys_alt :
            calc.SetToys         ( self.ntoys_null      , self.ntoys_alt      )
        if  0 != self.ntoys_null_tail or -1 != self.ntoys_alt_tail  :
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
## @class HybridCalculator
#  @see RooStats::HybridCalculator
class HybridCalculator (FrequentistCalculator) :
    """ Hybrid calcualator
    - see `ROOT.RooStats.HybridCalculator`
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
                   prior_null      = None ,
                   prior_alt       = None ) : 
        
        assert prior_null is None or isinstance ( prior_null , ( APDF1, ROOT.RooAbsPdf ) ) or \
               isinstance ( prior_null , string_types ) , "Invalid 'prior_null'!"        
        assert prior_alt  is None or isinstance ( prior_alt  , ( APDF1, ROOT.RooAbsPdf ) ) or \
               isinstance ( prior_alt  , string_types ) , "Invalid 'prior_alt'!"
        
        self.__prior_null      = prior_null
        self.__prior_alt       = prior_alt 

        ## initialize the base 
        FrequentistCalculator.__init__ ( self                              ,
                                         H1              = H1              ,
                                         H0              = H0              ,
                                         dataset         = dataset         ,
                                         sampler         = sampler         ,
                                         ntoys_null      = ntoys_null      ,
                                         ntoys_alt       = ntoys_alt       ,
                                         ntoys_null_tail = ntoys_null_tail ,
                                         ntoys_alt_tail  = ntoys_alt_tail  )  


        # 
        self.__prior_null_raw = None
        self.__prior_alt_raw  = None
        
        if   prior_null is None :  pass 
        elif isinstance ( prior_null , APDF1          ) : self.__prior_null_raw = prior_null.pdf
        elif isinstance ( prior_null , ROOT.RooAbsPdf ) : self.__prior_null_raw = prior_null 
        elif isinstance ( prior_null , string_types   ) :
            prior =  self.ws.pdf ( prior_null )
            assert valid_pointer ( prior )  and prior and isinstance ( prior , ROOT.RooAbsPdf ) , \
                   'Cannot get prior_null from Workspace!'
            self.__prior_null_raw = prior
            
        if   prior_alt is None :  pass 
        elif isinstance ( prior_alt  , APDF1          ) : self.__prior_alt_raw = prior_alt.pdf
        elif isinstance ( prior_alt  , ROOT.RooAbsPdf ) : self.__prior_alt_raw = prior_alt 
        elif isinstance ( prior_alt , string_types   ) :
            prior =  self.ws.pdf ( prior_alt  )
            assert valid_pointer ( prior )  and prior and isinstance ( prior , ROOT.RooAbsPdf ) , \
                   'Cannot get prior_alt from Workspace!'
            self.__prior_alt_raw = prior

    # ==============================================================================================
    ## Create and configure the calculator 
    def make_calculator ( self ) :
        """ Create and configure the calculator"""
        
        if self.sampler : calc = ROOT.RooStats.HybridCalculator ( self.cloned_dataset , self.h1 , self.h0 , self.sampler )
        else            : calc = ROOT.RooStats.HybridCalculator ( self.cloned_dataset , self.h1 , self.h0 )
        
        if -1 != self.ntoys_null      or -1 != self.ntoys_alt :
            calc.SetToys         ( self.ntoys_null      , self.ntoys_alt      )
        if  0 != self.ntoys_null_tail or -1 != self.ntoys_alt_tail  :
            calc.SetNToysInTails ( self.ntoys_null_tail , self.ntoys_alt_tail )

        if self.prior_null_raw : calc.ForcePriorNuisanceNull ( self.prior_null_raw )
        if self.prior_alt_raw  : calc.ForcePriorNuisanceAlt  ( self.prior_alt_raw  )
        
        return calc

    ## # =========================================================================
    ## def __del__ ( self ) :
    ##    FrequentistCalculator.__del__ ( self )
        
    @property
    def prior_null ( self ) :
        """'prior_null' : prior for null-hypothesis"""
        return self.__prior_null
    @property
    def prior_null_raw ( self ) :
        """'prior_null_raw' : prior for null-hypothesis (as `ROOT.RooAbsPdf`)"""
        return self.__prior_null_raw    
    @property
    def prior_alt ( self ) :
        """'prior_alt' : prior for alt-hypothesis"""
        return self.__prior_alt
    @property
    def prior_alt_raw ( self ) :
        """'prior_alt_raw' : prior for alt-hypothesis (as `ROOT.RooAbsPdf`)"""
        return self.__prior_alt_raw

# =============================================================================
## @class ProfileLikelihoodCalculator
#  @see RooStats::ProfileLikelihoodCalculator
class ProfileLikelihoodCalculator (Calculator) :
    """ ProfileLikelihood calculator
    - see `ROOT.RooStats.ProfileLikelihoodCalculator`
    """
    
    def __init__ ( self        ,
                   H0          ,   ## H0 model, e.g. signal+background
                   dataset     ,   ## dataset
                   null_params ) : ## null-parameters corresponding to background-only hypothesis

        Calculator.__init__ ( self , H0 , H0 , dataset )

        assert null_params , 'Invalid null-parameters!'
        
        self.__null_params_arg = null_params

        model = self.h0
        ws    = model.GetWorkspace()
        
        if   isinstance ( null_params , ROOT.RooArgSet  ) :
            self.__null_params = null_params
            
        elif isinstance ( null_params , ROOT.RooAbsReal ) :
            params = ROOT.RooArgSet ( null_params )
            self.__null_params = params
            
        elif isinstance ( null_params , dictlike_types  ) :            
            vlst = []
            vdct = {}            
            for key in null_params :                
                if   isinstance ( key , string_types ) : v = ws.var ( key      )
                else                                   : v = ws.var ( key.name )                
                assert v , "No valid variable for '%s'" % key                
                vv = float ( null_params [ key ] )
                vlst.append ( v )
                vdct [ v.name ]  = vv
            vset = ROOT.RooArgSet() 
            for v in vlst : v.setVal ( vdct [ v.name ] )
            for v in vlst : vset.add ( v )
            self.__null_params = vset

        elif isinstance ( null_params , sequence_types  ) :
            params = ROOT.RooArgSet ()
            for p in null_params : params.add ( p ) 
            self.__null_params = params 

        else :
            raise TypeError('Invalid null-params!')
        
        sname = 'null_parameters_for_%s' % model.name
        ws    .defineSet ( sname , self.__null_params )  

    ## create and configure the Profile Likelihood calculator 
    def make_calculator ( self ) :
        """ Create and configure the PofileLikelihood calculator"""
        calc = ROOT.RooStats.ProfileLikelihoodCalculator ( self.dataset ,  self.h0 )
        calc.SetNullParameters ( self.null_params )        
        return calc

    @property
    def null_params ( self ) :
        """'null_params' : set of parameters that define null-model"""
        return self.__null_params

# =============================================================================
# @class HypoTestInverter
# @see RooStats::HypoTestInverter
class HypoTestInverter(object) :
    """ Hypo test inverter 
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
        """ Define (or perform the actual scan)
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

    # =========================================================================
    ## define (or perform the actual scan)
    #  @code
    #  hti = ...
    #  hti.scan ( 100 , 0.0 , 10.0 )   ## define scan with 100 point between 0 and 10
    #  hti.scan ()                     ## define auto scan 
    #  hti.scan (  [ 0, 1, 2, 3, 4 ] ) ## define custom scane    
    #  @endcode
    def scan_with_progress ( self , *values ) :
        """ Define (or perform the actual scan)
        >>> hti = ...
        >>> hti.scan ( 100 , 0.0 , 10.0 )   ## define scan with 100 point between 0 and 10
        >>> hti.scan ()                     ## define auto scan 
        >>> hti.scan (  [ 0, 1, 2, 3, 4 ] ) ## define custom scane    
        """

        from ostap.utils.utils        import vrange
        from ostap.utils.progress_bar import progress_bar 
        
        if    not values :
            
            self.__inverter.SetAutoScan  ()
            if self.__interval : self.__interval = None
            if self.__plot     : self.__plot = None
            return
        
        elif  3 == len ( values )  and \
             isinstance ( values [ 0 ] , integer_types ) and \
             2 <= values                                 and \
             isinstance ( values [ 1 ] , num_types     ) and \
             isinstance ( values [ 2 ] , num_types     ) :            
            return self.scan_with_progress ( vrange ( values[1] , values [2] , values[0] ) ) 

        elif  1 == len ( values ) and isinstance ( values [ 0 ] , sequence_types  ) :
            for v in progress_bar ( values [ 0 ] , description = 'Scan:' ) :
                self.__inverter.RunOnePoint ( v )                
        else :            
            for v in progress_bar ( values , description = 'Scan:' ) :
                self.__inverter.RunOnePoint ( v )

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
        """ Reset/clear current interval & plot
        """
        if self.__interval : self.__interval = None
        if self.__plot     : self.__plot = None 
        
# =============================================================================
## default color for Brasil-plot bands
band_colors = ( ROOT.kGreen , ROOT.kYellow  ,
                ROOT.kCyan  , ROOT.kMagenta ,
                ROOT.kPink  , ROOT.kOrange  )
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
    """ Helper class to create and keep the 'Brasil-band' plot
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
            self.__colors.append ( ic  ) 
            ic += 1

        self.__colors = tuple ( self.__colors )
        
        self.__plot   = None
        self.__legend = None
        
    # =========================================================================
    ## Add the point to the Brasil-plot
    def fill ( self , x , observed , hti_result ) : 
        """ Add the point to the Brasil-plot
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

        ## create "band-graphs"
        gr_observed = ROOT.TGraph ( np )
        gr_median   = ROOT.TGraph ( np ) 
        gr_bands    = [ ROOT.TGraphAsymmErrors ( np ) for i in range ( ngb ) ]

        for i , x in enumerate ( sorted ( self.__data ) ) :

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
    ## get the standard legend for the plot
    #  @code 
    #  bp = BrasilBand ( ... ) 
    #  for .... :
    #    bp.fill ( ... ) 
    #  bp.plot.draw('a')
    #  bp.legend.draw('a')
    #  @endcode 
    @property
    def legend ( self )  :
        """'legend' : get the standarde legend for the plot
        >>> bp = BrasilBand ( ... ) 
        >>> for .... :
        >>>    bp.fill ( ... ) 
        >>> bp.plot.draw('a')
        >>> bp.legend.draw('a')
        """
        if not self.__legend : plot = self.plot
        return self.__legend

    # ===========================================================================
    ## number of points in the graph 
    def __len__ ( self ) :
        """ Number of points in the graph
        """
        return len ( self.__data ) 

# =============================================================================
## @class P0Plot
#  Helper class to create graphs(s) for p0-scan plot
#  @code
#  plot = P0Plot() 
#  for ,,, :
#      value       = ...
#      p0          = ...
#      p0_expected = 
#      plot.fill ( value , p0 , p0_expected )
#  ## plot p-values 
#  plot.p0.draw          ( 'ac')
#  ## plot expected p-values 
#  plot.p0_expected .draw( 'c')
#  ## plot #sigmas  
#  plot.sigmas.draw      ( 'ac')
#  @endcode
#  @thanks Dima Golubkov 
class P0Plot(object) :
    """ Helper class to create graphs(s) for p0-scan plot
    >>> plot = P0Plot() 
    >>> for ,,, :
    ...    value       = ...
    ...    p0          = ...
    ...    p0_expected = 
    ...    plot.fill ( value , p0 , p0_expected )
    
    >>> plot.p0         .draw ( 'ac') ## plot p-values 
    >>> plot.p0_expected.draw ( 'c' ) ## plot expected p-values 
    >>> plot.sigmas     .draw ( 'ac') ## plot #sigmas  
    - thanks to Dima Golubkov! 
    """
    def __init__ ( self ) :

        self.__p0           = ROOT.TGraph() 
        self.__sigmas       = ROOT.TGraph() 
        self.__p0_expected  = None
        self.__rows         = [ ( 'value' , 'p0' , '#sigma' , 'p0(exp)' ) ] 
        
        self.__p0      .blue ()
        self.__sigmas  .red  ()
        self.__p0    .SetLineWidth(2) 
        self.__sigmas.SetLineWidth(2) 

    # ==========================================================================
    ## add the point to TGraph(s)
    #  Full version 
    #  @code
    #  plot = P0Plot() 
    #  for value in ...  :
    #     calculator = ...
    #     hypo_test  = calcualtor.hypo_test
    #     p0    = ht.     NullPValue () 
    #     p0alt = ht.AlternatePValue ()
    #     p0exp = ROOT.RooStats.AsymptoticCalculator.GetExpectedPValues (  p0 , p0alt , 0 , False ) 
    #     plot.fill ( value ,p0 , p0exp ) 
    #  @endcode
    #  Alternative version:
    #  @code
    #  plot = P0Plot() 
    #  for value in ...  :
    #     calculator = ...
    #     plot.fill ( value , calculator ) 
    #  @endcode
    def fill ( self , value , *what ) :
        """ Add the point to TGraph(s
        
        - Full version 
        >>> plot = P0Plot() 
        >>> for value in ...  :
        ...     calculator = ...
        ...     hypo_test  = calcualtor.hypo_test
        ...     p0    = ht.     NullPValue () 
        ...     p0alt = ht.AlternatePValue ()
        ...     p0exp = ROOT.RooStats.AsymptoticCalculator.GetExpectedPValues (  p0 , p0alt , 0 , False ) 
        ...     plot.fill ( value ,p0 , p0exp ) 
        - Alternative version:
        >>> plot = P0Plot() 
        >>> for value in ...  :
        ...     calculator = ...
        ...     plot.fill ( value , calculator ) 
        """
        n = len ( what )
        if n == 1 and isinstance ( what [0] , num_types ) and 0 <= what [ 0 ] <= 1 :
            p0          = what
            p0_expected = None 
        elif n == 1 and isinstance ( what [0] , Calculator ) :
            calculator  = what [ 0 ] 
            ht          = calculator.hypo_test
            p0          = ht.     NullPValue ()
            p0alt       = ht.AlternatePValue ()
            p0_expected = ROOT.RooStats.AsymptoticCalculator.GetExpectedPValues (  p0 , p0alt , 0 , False ) 
        elif n == 1 and isinstance ( what [0] , ROOT.RooStats.HypoTestCalculator ) :
            calculator  = what [ 0 ] 
            ht          = calculator.GetHypoTest() 
            p0          = ht.     NullPValue ()
            p0alt       = ht.AlternatePValue ()
            p0_expected = ROOT.RooStats.AsymptoticCalculator.GetExpectedPValues (  p0 , p0alt , 0 , False ) 
        elif n == 1 and isinstance ( what [0] , ROOT.RooStats.HypoTestResult  ) :
            ht          = what [ 0 ] 
            p0          = ht.     NullPValue ()
            p0alt       = ht.AlternatePValue ()
            p0_expected = ROOT.RooStats.AsymptoticCalculator.GetExpectedPValues (  p0 , p0alt , 0 , False ) 
        elif n == 2 and  all ( ( isinstance ( i , num_types ) and 0 <= i <= 1 )  for i in what ) :
            p0          = float ( what [ 0 ] ) 
            p0_expected = float ( what [ 1 ] ) 
        else :
            raise TypeError ( "Unknown 'what' %s" % str ( what ) ) 


        n1 = len ( self.__p0  )
        ns = ROOT.RooStats.PValueToSignificance ( p0 )
        
        self.__p0    .SetPoint ( n1 , value , p0 )
        self.__sigmas.SetPoint ( n1 , value , ns )

        row = '%+.5g' % value , '%.5g' % p0 , '%+.2f' % ns 
        
        if not p0_expected is None :
            if not self.__p0_expected :
                self.__p0_expected = ROOT.TGraph()
                self.__p0_expected.green ()
                self.__p0_expected.SetLineWidth(2) 
                
            n2 = len ( self.__p0_expected )
            self.__p0_expected.SetPoint ( n2 , value , p0_expected )
            row = row + ( '%.4g' % p0_expected , )
            
        self.__rows.append ( row )

    ## number of points in the graph 
    def __len__ ( self ) :
        """Number of points in the graph
        """
        return len ( self.__graph_p0  )
    
    @property
    def sigmas ( self ) :
        """'sigmas' : graph of significances """
        return self.__sigmas 

    @property
    def p0 ( self ) :
        """'p0' : graph of p0-values """
        return self.__p0 
    
    @property
    def p0_expected ( self ) :
        """'p0_expected' : graph of expected  p0-values """
        return self.__p0_expected
    
    # =============================================================================
    ## get the summary table 
    def table ( self , title = '' , prefix = '' ) :
        """ Get the summary table """
        import ostap.logger.table as T
        title = title if title else 'p0-scan'
        return T.table ( self.__rows , title = title , prefix = prefix , alignment = 'lccc' )

# =============================================================================
## Print ModelConfig as table 
def mc_table ( mc , title = '' , prefix = '' ) :
    """ Print ModelConfig as table 
    """
    rows = [  ( 'Model Config' , 'Value' ) ]

    row = 'Name' ,  mc.name
    rows.append ( row )

    row = 'Title' ,  mc.title 
    rows.append ( row )
    
    pdf = mc.GetPdf() 
    if valid_pointer ( pdf ) :
        row = 'PDF' , '%s: %s/%s' % ( typename  ( pdf ) , pdf.name , pdf.title )                                       
        rows.append ( row )

    pars = mc.GetObservables()
    if valid_pointer ( pars ) :
        row = 'Observables:'            , ', '.join ( v.name for v in pars )  
        rows.append ( row )

    pars = mc.GetParametersOfInterest()
    if valid_pointer ( pars ) : 
        row = 'Parameters-of-interest'  , ', '.join ( v.name for v in pars  ) 
        rows.append ( row )

    pars = mc.GetNuisanceParameters()
    if valid_pointer ( pars ) : 
        row = 'Nuisance parameters'     , ', '.join ( v.name for v in pars )
        rows.append ( row )

    pars = mc.GetGlobalObservables() 
    if valid_pointer ( pars ) : 
        row = 'Global observables:'     , ', '.join ( v.name for v in pars ) 
        rows.append ( row )

    pars = mc.GetConditionalObservables()
    if valid_pointer ( pars ) : 
        row = 'Conditional observables:'     , ', '.join ( v.name for v in pars ) 
        rows.append ( row )
        
    pars = mc.GetConstraintParameters () 
    if valid_pointer ( pars ) : 
        row = 'Constraint parameters:'     , ', '.join ( v.name for v in pars ) 
        rows.append ( row )

    if ( 6 , 29 , 10 ) <= root_info :
        pars = mc.GetConstraintParameters () 
        if valid_pointer ( pars ) :
            for i , c in enumerate ( pars , start = 1 )  :
                row = 'External constraint #%d:' % i , '%s: %s/%s' % ( typename ( c ) , c.name , c.title )
                rows.append ( row )
        
    pdf = mc.GetPriorPdf() 
    if valid_pointer ( pdf ) :
        row = 'Prior-PDF' , '%s: %s/%s' % ( typename  ( pdf ) , pdf.name , pdf.title )                                       
        rows.append ( row )

    data = mc.GetProtoData() 
    if valid_pointer ( data ) : 
        vars = data.get() 
        row  = 'ProtoData' , '%s(%s,#%d){%s}' % ( typename ( data ) ,
                                                  data.name         ,
                                                  len      ( data ) ,
                                                  ', '.join ( v.name for v in vars ) ) 
    pars = mc.GetSnapshot()
    if valid_pointer ( pars ) :
        vmax = -1
        for v in pars : vmax = max ( vmax , len ( v.name ) ) 
        vmax += 2
        vfmt = '%%%ds : %%-+7g' % vmax
        for i , v in enumerate ( pars ) :
            if 0 == i : row = 'Snapshot:' , vfmt % ( v.name , v.getVal () )
            else      : row = ''          , vfmt % ( v.name , v.getVal () )
            rows.append ( row )

    import ostap.logger.table as T
    title = title if title else '%s(%s,%s)' % ( typename ( mc ) , mc.name , mc.title )
    return T.table ( rows , title = title , prefix = prefix , alignment = 'lw' )


ROOT.RooStats.ModelConfig.table    = mc_table
ROOT.RooStats.ModelConfig.__repr__ = mc_table
ROOT.RooStats.ModelConfig.__str__  = mc_table

_new_methods_ += [  
    ROOT.RooStats.ModelConfig.table    ,
    ROOT.RooStats.ModelConfig.__repr__ ,
    ROOT.RooStats.ModelConfig.__str__  ,
]
# =============================================================================
## Print RooWorkspace as table 
def ws_table ( ws , title = '' , prefix = '' ) :
    """ Print EooWorkspace  as table 
    """
    rows = [  ( 'Component' , 'Label', 'Value' ) ]

    row = 'Name'  , '' , ws.name
    rows.append ( row )

    row = 'Title' , '' , ws.title 
    rows.append ( row )

    printed = set()

    values = ws.allPdfs()
    for i,v in enumerate ( values , start = 1 ) : 
        row = 'PDF'      ,          '%s' % v.name   , '%s(%s,%s)' % ( typename ( v ) , v.name , v.title )
        rows.append ( row ) 

    values = ws.allFunctions()
    for i,v in enumerate ( values , start = 1 ) : 
        row = 'Function' ,          '%s' % v.name   , '%s(%s,%s)' % ( typename ( v ) , v.name , v.title )
        rows.append ( row ) 

    values = ws.allCatFunctions()
    for i,v in enumerate ( values , start = 1 ) : 
        row = 'Category function' , '%s' % v.name  , '%s(%s,%s)' % ( typename ( v ) , v.name , v.title )
        rows.append ( row ) 

    values = ws.allResolutionModels()
    for i,v in enumerate ( values , start = 1 ) : 
        row = 'Resolution model' , '%s' % v.name   , '%s(%s,%s)' % ( typename ( v ) , v.name , v.title )
        rows.append ( row ) 

    values = ws.allGenericObjects ()
    for i,v in enumerate ( values , start = 1 ) : 
        row = 'Generic object' , '%s' % ( v.GetName () if hasattr ( v , 'GetName'  ) else '' ) , \
            '%s(%s,%s)' % ( typename ( v ) ,  
                            v.GetName  () if hasattr ( v , 'GetName'  ) else '' ,
                            v.GetTitle () if hasattr ( v , 'GetTitle' ) else '' ) 
        rows.append ( row ) 


    snapshots = ws.getSnapshots()
    for i, snapshot  in enumerate ( snapshots , start = 1 ) :
        fmt = '%s : %-+7g'
        row = 'Snapshot'  , '%s' % snapshot.GetName() , '{ %s }' %  ( ', '.join ( ( fmt % ( v.name , v.getVal() ) ).strip() for v in snapshot ) )
        rows.append ( row ) 
        
    sets = ws.sets()
    keys = sorted (  ( str ( k.first ) for k in sets ) ) 
    for i, key in enumerate ( keys , start = 1 ) :
        the_set = ws.set (  key  )
        row = 'Named set' , '%s' % key  , '{ %s }' % ( ', '.join ( sorted ( v.name for v in the_set ) ) ) 
        rows.append ( row )
        
    values = ws.allVars()
    row = 'Variables'  , '#%d' % len ( values ) , ', '.join ( sorted ( v.name for v in values ) ) 
    rows.append ( row ) 
                                     
    values = ws.allCats()
    for i , c in enumerate ( values , start = 1 ) : 
        row = 'Category' , '%s' % c.name , '{ %s }' % ( ', '.join (  ( '%s:%s' % ( l ,k ) ) for ( l , k ) in c.items() ) )
        rows.append ( row ) 

    values = ws.allData()
    for i, v in enumerate ( values , start = 1 ) :
        vars = v.get() 
        row = 'Data' , '%s' % v.name  , '%s(%s,%s) #%d entries { %s }' % ( typename ( v ) , v.name , v.title , len ( v ) ,
                                                                           ', '.join ( sorted ( p.name for p in vars ) ) ) 
        rows.append ( row ) 
                                                            
    values = ws.allEmbeddedData()
    for i, v in enumerate ( values , start = 1 ) :
        vars = v.get() 
        row = 'Embedded Data' , '%s' % v.name  , '%s(%s,%s) #%d entries { %s }' % ( typename ( v ) , v.name , v.title , len ( v ) ,
                                                                                    ', '.join ( sorted ( p.name for p in vars ) ) ) 
        rows.append ( row ) 

        
    import ostap.logger.table as T
    title = title if title else '%s(%s,%s)' % ( typename ( ws ) , ws.name , ws.title )
    return T.table ( rows , title = title , prefix = prefix , alignment = 'lrw' )

    
ROOT.RooWorkspace.table    = ws_table
ROOT.RooWorkspace.__repr__ = ws_table
ROOT.RooWorkspace.__str__  = ws_table

_new_methods_ += [ 
    ROOT.RooWorkspace.table    , 
    ROOT.RooWorkspace.__repr__ ,
    ROOT.RooWorkspace.__str__  , 
]

if root_info < ( 6 , 28 ) and not hasattr ( ROOT.RooWorkspace , 'Import' ) :
    ROOT.RooWorkspace.Import = getattr ( ROOT.RooWorkspace , 'import' )
    _new_methods_ . append ( ROOT.RooWorkspace.Import )

_decorated_classes_ = (
    ROOT.RooStats.ModelConfig , 
    ROOT.RooWorkspace         , 
    )

_new_methods_ = tuple ( _new_methods_ ) 
    
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
