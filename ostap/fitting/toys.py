#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file ostap/fitting/toys.py
#  Simple utilities to run fitting toys, Jackknife, bootstrap, etc...    
#  @date   2020-01-18
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
# =============================================================================
""" Simple utilities to run fitting toys, Jackknife, bootstrap, ..
"""
# =============================================================================
__author__  = 'Vanya BELYAEV  Ivan.Belyaev@itep.ru'
__date__    = "2020-01-18"
__version__ = '$Revision$'
__all__     = (
    "make_toys"        , ## run fitting toys (the same PDF to generate and fit)
    "make_toys2"       , ## run fitting toys (separate models to generate and fit)
    'make_jackknife'   , ## run Jackknife analysis 
    'make_bootstrap'   , ## run Bootstrapanalysis 
    "vars_transform"   , ## helper fnuction to transform the variables
    "print_stats"      , ## print toys      statistics 
    "print_jackknife"  , ## print jackknife statistics 
    "print_bootstrap"  , ## print bootstrap statistics 
    )
# =============================================================================
import ROOT
from   builtins          import range
from   ostap.core.core   import VE
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.toys' )
else                       : logger = getLogger( __name__             )
# =============================================================================
logger.debug ( 'Utilities to run fitting toys, Jackknife, bootstrap, ... ')
# ==============================================================================
## Technical transformation  to the dictionary :  { 'name' : float_value }
def vars_transform ( vars ) :
    """Technical transformation to the dictionary :  `{ 'name' : float_value }`
    """
    
    from   ostap.core.ostap_types     import dictlike_types
    import ostap.fitting.roofitresult
    import ostap.fitting.variables 
    
    result = {}
    if isinstance   ( vars , ROOT.RooFitResult ) :
        rs = vars.dct_parsms ()
        for p in  rs  : result [ p ] = float ( rs   [ p ] )
    elif isinstance ( vars , dictlike_types ) :
        for p in vars : result [ p ] = float ( vars [ p ] ) 
    else :
        for p in vars :
            if not isinstance ( p , ROOT.RooAbsCategory ) :
                result [ p.GetName()  ] = float ( p )
            
    return result

# =============================================================================
## print statistics of pseudoexperiments
def print_stats (  stats , ntoys = '' , logger = logger ) :
    """print statistics of pseudoexperiments
    """
    
    table = [ ( 'Parameter' , '#', 'mean' , 'rms' , '%13s / %-13s' % ( 'min' , 'max' ) ) ] 

    def make_row ( c ) :
        n      = "{:^11}".format ( c.nEntries() )
        mean   = c.mean ()
        mean   = "%+13.6g +- %-13.6g" % ( mean.value() , mean.error() )
        rms    = "%13.6g"             % c.rms ()
        minmax = "%+13.6g / %-+13.6g" % ( c.min() , c.max () ) 
        return p , n , mean , rms  , minmax 
        
    for p in sorted ( stats )  :
        if    p.startswith('pull:') : continue  
        c      = stats [ p ]
        table.append (  make_row ( c )  )
        
    for p in sorted ( stats )  :
        if not p.startswith('pull:') : continue  
        c      = stats [ p ]
        table.append (  make_row ( c )  )

    if not ntoys :
        ntoys = 0
        for s in stats :
            ntoys = max ( ntoys , stats[s].nEntries() )
            
    import ostap.logger.table as Table
    table = Table.table ( table                                    ,
                          title     = "Results of %s toys" % ntoys ,
                          alignment = 'lcccc'                      ,
                          prefix    = "# "                         )
    logger.info ( 'Results of %s toys:\n%s' % ( ntoys , table ) ) 


# =============================================================================
## Jackknife estimator from jackknife statistic
#  @code
#  statistics = ....
#  jackknife              = jackknife_estimator ( statistics         )
#  jackknife , theta_jack = jackknife_esiimator ( statistics , value )
#  @endcode 
def jackknife_statistics ( statistics , theta = None ) :
    """Jackknife estimator from jackknife statistic
    >>> statistics = ....
    >>> jacknife                      = jackknife_estimator ( statistics         )
    >>> jacknife , theta_corr , bias  = jackknife_esiimator ( statistics , value )
    """
    assert isinstance ( theta , ( VE, None) ) ,\
           "jackknife_statistics: invalid type of ``value'' %s" % type ( value ) 
    
    N         = statistics . nEntries ()             ## number of jackknife samples 
    theta_dot = statistics . mean     ()             ## mean over jackknife samples 
    jvar      = statistics . variance () * ( N - 1 ) ## variance 

    jackknife = VE ( theta_dot.value () , jvar )     ## jackknife estimate 
    
    ## get Jackknife estimator from statistics 
    if theta is None : return jackknife
    
    bias       =  ( N - 1 ) * ( theta_dot.value() - theta.value() )
    
    ## theta  corrected for the bias and with Jackknife variance  
    theta_jack = VE ( theta.value() - bias , jvar     ) ## corrected value 

    return jackknife , theta_jack

# =============================================================================
## print Jackknife statistics
def print_jackknife  ( fitresult , stats , logger = logger ) :
    """print Jackknife statistics
    """
    
    header = ( 'Parameter' , 'theta' , 'theta_(.)' ,  'theta_jack' , 'bias/sigma [%]' , 'error [%]' ) 
    table  = []
    
    for p in fitresult :
        name  = p.name
        if not name in stats : continue

        statistics  = stats [ name ]
        
        theta = p * 1.0 ## fitted value

        ## jackknife estimates 
        jackknife , theta_jack = jackknife_statistics ( statistics , theta )

        bias  = theta_jack.value () - theta     .value ()        
        scale = theta     .error () / theta_jack.error () 

        row = ( name , 
                "%+13.6g +- %-13.6g" % ( theta      . value () , theta      .error () ) , 
                "%+13.6g +- %-13.6g" % ( jackknife  . value () , jackknife  .error () ) , 
                "%+13.6g +- %-13.6g" % ( theta_jack . value () , theta_jack .error () ) ,  
                '%+6.2f'             % ( bias / theta.error() * 100 ) , 
                '%+6.2f'             % ( scale * 100 - 100 )          )
        
        table.append ( row )
        
    table.sort()
    table = [ header ] + table 
 
    import ostap.logger.table as Table
    table = Table.table ( table                           ,
                          title     = "Jackknife results" ,
                          alignment = 'lcccccc'             ,
                          prefix    = "# "                )
    logger.info ( 'Jackknife results:\n%s' % table )
    

# =============================================================================
## print Bootstrap statistics
def print_bootstrap  ( fitresult , stats , logger = logger ) :
    """print Bootstrap statistics
    """
    
    header = ( 'Parameter' , 'theta' ,  'theta_boot' , 'bias/sigma [%]' , 'error [%]' ) 
    table  = []

    n = 0 
    for p in fitresult :
        name  = p.name
        if not name in stats : continue

        statistics  = stats [ name ]

        n = max ( n , statistics.nEntries() ) 
        theta      = p * 1.0 ## fitted value

        theta_boot = VE ( statistics.mean().value() , statistics.mu2() ) 
        

        bias  = theta_boot.value () - theta     .value ()        
        scale = theta     .error () / theta_boot.error () 

        row = ( name , 
                "%+13.6g +- %-13.6g" % ( theta      . value () , theta      .error () ) , 
                "%+13.6g +- %-13.6g" % ( theta_boot . value () , theta_boot .error () ) ,  
                '%+6.2f'             % ( bias / theta.error() * 100 ) , 
                '%+6.2f'             % ( scale * 100 - 100 )          )
        
        table.append ( row )
        
    table.sort()
    table = [ header ] + table 
 
    import ostap.logger.table as Table
    table = Table.table ( table                                   ,
                          title     = "Bootstrapping with #%d samples" % n ,
                          alignment = 'lcccc'             ,
                          prefix    = "# "                )
    logger.info ( 'Bootstrapping with #%d samples:\n%s' % ( n , table ) )
    

# ==============================================================================
## Default function to generate the data
#  - simple call for <code>PDF.generate</code>
def generate_data ( pdf , varset , **config ) :
    """Default function to generate the data
    - simple call for `PDF.generate`
    """
    return pdf.generate ( varset = varset , **config )

# ==============================================================================
## Default function to perform the actual fit
#  - simple call for <code>PDF.fitTo</code>
def make_fit    ( pdf , dataset , **config ) :
    """Default function to  perform the actual fit
    - simple call for `PDF.fitTo`
    """
    result , _ = pdf.fitTo ( dataset , **config )
    return result

# ==============================================================================
## Accept fit?
#  Accept the fit result?
#   - valid fit result
#   - fit status is 0 (SUCCESS)
#   - covariance matrix quality  is either 3(full accurate matrix) or -1 (unknown/externbally provided?) 
#  @param result  fit result
#  @param pdf     pdf
#  @param dataset pdf
#
def accept_fit  ( result , pdf = None , dataset = None ) :
    """Accept the fit result?
    - valid fit result
    - fit status is 0 (SUCCESS)
    - covariance matrix quality  is either 0 or -1 
    """
    return result and ( 0 == result.status () ) and ( result.covQual () in ( -1 , 3 ) ) 


# ==============================================================================
## make <code>nToys</code> pseudoexperiments
#
#  Schematically:
#  @code
#  for toy in range ( nToys )  :
#  ...  dataset = gen_fun ( pdf , ...     , **gen_config )
#  ...  result  = fit_fun ( pdf , dataset , **fit_config )
#  ...  if not accept_fun ( result , pdf , dataset ) : continue
#  .... < collect statistics here > 
#  @endcode
#
#  For each experiment
#  - generate dataset using <code>pdf</code> with variables specified
#    in <code>data</code> and configuration specified via<code>gen_config</code>
#    for each generation the parameters of <code>pdf</code> are reset
#    for their initial values and values from <code>init_pars</code>
#  - fit generated dataset  with <code>pdf</code> using configuration
#    specified via  <code>fit_config</code>
#
#  @code
#  pdf = ...
#  results , stats = make_toys ( pdf   ,   ## PDF  to use 
#     nToys      = 1000        ,           ## Number of pseudoexperiments 
#     data       = [ 'mass' ]  ,           ## variables in dataset 
#     gen_config = { 'nEvents' : 5000 } ,  ## configuration of <code>pdf.generate</code>
#     fit_config = { 'ncpus'   : 2    } ,  ## configuration of <code>pdf.fitTo</code>
#     init_pars  = { 'mean' : 0.0 , 'sigma' : 1.0 } ) ## parameters to use for generation 
#  @endcode 
#
# Derived parameters can be also retrived via <code>more_vars</code> argument:
# @code
# ratio     = lambda res,pdf : res.ratio('x','y')
# more_vars = { 'Ratio' : ratio }
#  r,  s = make_toys ( .... , more_vars = more_vars , ... ) 
# @endcode
#
# @param pdf        PDF to be used for generation and fitting
# @param nToys      number    of pseudoexperiments to generate
# @param data       variable list of variables to be used for dataset generation
# @param gen_config configuration of <code>pdf.generate</code>
# @param fit_config configuration of <code>pdf.fitTo</code>
# @param init_pars  redefine these parameters for each pseudoexperiment
# @param more_vars  calculate more variables form fit-result
# @param get_fun    specific generate action (if needed) 
# @param fit_fun    specific fitting action (if needed) 
# @param accept_fun specific accept action (if needed) 
# @param silent     silent toys?
# @param progress   show the progress?
# @return dictionary with fit results for the toys and the dictionary of statistics
#
#  - If <code>gen_fun</code>    is not specified <code>generate_data</code> is used 
#  - If <code>fit_fun</code>    is not specified <code>make_fit</code>      is used 
#  - If <code>accept_fun</code> is not specified <code>accept_fit</code>    is used   
def make_toys ( pdf                 ,
                nToys               , 
                data                , ## template for dataset/variables 
                gen_config          , ## parameters for <code>pdf.generate</code>   
                fit_config = {}     , ## parameters for <code>pdf.fitTo</code>
                init_pars  = {}     ,
                more_vars  = {}     ,
                gen_fun    = None   , ## generator function ( pdf , varset  , **config )
                fit_fun    = None   , ## fit       function ( pdf , dataset , **config )
                accept_fun = None   , ## accept    function ( fit-result, pdf, dataset )
                silent     = True   ,                
                progress   = True   ,
                logger     = logger ) :
    """Make `nToys` pseudoexperiments

    -   Schematically:
    >>> for toy in range ( nToys )  :
    >>> ...  dataset = gen_fun ( pdf , ...     , **gen_config )
    >>> ...  result  = fit_fun ( pdf , dataset , **fit_config )
    >>> ...  if not accept_fun ( result , pdf , dataset ) : continue
    >>> .... < collect statistics here > 
    
    For each pseudoexperiment:

    1. generate dataset using `pdf` with variables specified
    in `data` and configuration specified via `gen_config`
    for each generation the parameters of `pdf` are reset
    for their initial values and valeus from `init_pars`
    
    2. fit generated dataset  with `pdf` using configuration
    specified via  `fit_config`

    - pdf        PDF to be used for generation and fitting
    - nToys      number    of pseudoexperiments to generate
    - data       variable list of variables to be used for dataset generation
    - gen_config configuration of <code>pdf.generate</code>
    - fit_config configuration of <code>pdf.fitTo</code>
    - init_pars  redefine these parameters for each pseudoexperiment
    - more_vars  dictionary of functions to define the additional results
    - gen_fun    generator function
    - fit_fun    fitting   function
    - accept_fun accept    function
    - silent     silent toys?
    - progress   show progress bar? 
    
    It returns a dictionary with fit results for the toys and a dictionary of statistics
    
    >>> pdf = ...
    ... results, stats = make_toys ( pdf     , ## PDF  to use 
    ...                 1000                 , ## number of toys 
    ...                 [ 'mass' ]           , ## variables in dataset 
    ...                 { 'nEvents' : 5000 } , ## configuration of `pdf.generate`
    ...                 { 'ncpus'   : 2    } , ## configuration of `pdf.fitTo`
    ...                 { 'mean' : 0.0 , 'sigma' : 1.0 } ## parameters to use for generation 
    ...                )

    Derived parameters can be also retrived via <code>more_vars</code> argument:
    >>> ratio    = lambda res,pdf : res.ratio('x','y') 
    >>> more_vars = { 'Ratio' : ratio }
    >>> r,  s = make_toys ( .... , more_vars = more_vars , ... ) 

    - If `gen_fun`    is not specified `generate_data` is used 
    - If `fit_fun`    is not specified `make_fit`      is used 
    - If `accept_fun` is not specified `accept_fit`    is used 
    """

    from ostap.core.ostap_types import string_types, integer_types  
    
    assert isinstance ( nToys , integer_types ) and 0 < nToys,\
           'Invalid "nToys" argument %s/%s' % ( nToys , type ( nToys ) )
    
    assert gen_config and 'nEvents' in gen_config,\
           'Number of events per toy must be specified via "gen_config" %s' % gen_config

    ## 1. generator function? 
    if gen_fun is None :
        if not silent : logger.info ( "make_toys: use default ``generate_data'' function!")
        gen_fun = generate_data 
    assert gen_fun and callable ( gen_fun ) , 'Invalid generator function!'
    
    ## 2. fitting function? 
    if fit_fun is None :
        if not silent : logger.info ( "make_toys: use default ``make_fit'' function!")
        fit_fun = make_fit 
    assert fit_fun and callable ( fit_fun ) , 'Invalid fit function!'

    ## 3. accept function? 
    if accept_fun is None :
        if not silent : logger.info ( "make_toys: use default ``accept_fit'' function!")
        accept_fun = accept_fit
    assert accept_fun and callable ( accept_fun ) , 'Invalid accept function!'
        
    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult
    import ostap.fitting.basic 

    params = pdf.params ()
    varset = ROOT.RooArgSet() 
    
    if isinstance ( data , ROOT.RooAbsData       ) : varset = data.varset() 
    else :
        for v in data :
            if   isinstance ( v , ROOT.RooAbsArg ) :
                varset.add ( v )
            elif isinstance ( v , string_types   ) and v in params :
                varset.add ( params [ v ] )
            else :
                raise TypeError('Invalid variable %s/%s' % ( v , type ( v ) ) )


    fix_pars = vars_transform ( params    ) 
    fix_init = vars_transform ( init_pars ) 

    pdf.load_params ( params = fix_pars , silent = silent )
    pdf.load_params ( params = fix_init , silent = silent )

    ## save all initial parameters (needed for the final statistics)
    params  = pdf.params      ()
    fix_all = vars_transform  ( params ) 
    
    fitcnf = {}
    fitcnf.update ( fit_config )
    if not 'silent' in fitcnf : fitcnf [ 'silent' ] = silent

    from collections import defaultdict 
    results = defaultdict(list) 

    from   ostap.core.core        import SE, VE

    fits = defaultdict ( SE )  ## fit statuses 
    covs = defaultdict ( SE )  ## covariance matrix quality
    
    ## run pseudoexperiments
    from ostap.utils.progress_bar import progress_bar 
    for i in progress_bar ( range ( nToys ) , silent = not progress ) :
                
        ## 1. reset PDF parameters 
        pdf.load_params ( params = fix_pars  , silent = silent )
        pdf.load_params ( params = init_pars , silent = silent )

        ## 2. generate dataset!  
        ## dataset = pdf.generate ( varset = varset , **gen_config )  
        dataset = gen_fun ( pdf , varset = varset , **gen_config )  
        if not silent : logger.info ( 'Generated dataset #%d\n%s' % ( i , dataset ) )
        
        ## 3. fit it!
        r = fit_fun ( pdf , dataset , **fitcnf ) 

        ## fit status 
        fits [ r.status  () ] += 1

        ## covariance matrix quality
        covs [ r.covQual () ] += 1
              
        ## ok ?
        if accept_fun ( r , pdf , dataset ) : 
            
            ## 4. save results 
            rpf = r.params ( float_only = True ) 
            for p in rpf : 
                results [ p ].append ( rpf [ p ][0] ) 
                
            for v in more_vars :
                func  = more_vars[v] 
                results [ v ] .append ( func ( r , pdf ) )

            results [ '#' ] .append ( len ( dataset ) )
            
        dataset.clear()
        del dataset
        del r

    ## make a final statistics 
    stats = defaultdict ( SE )

    for par in results :
        pars = results [ par ]
        mvar = par in more_vars 
        if not mvar : a0 = fix_all.get ( par , None  )
        for v in pars : 
            v0 = float ( v )         
            stats     [ par             ] +=   v0
            if not mvar and not a0 is None and isinstance ( v , VE ) and 0 < v.error() : 
                stats [ 'pull:%s' % par ] += ( v0 - a0 ) / v.error()

    for k in fits :
        stats ['- Status  %s' % k ] = fits [ k ]
    for k in covs :
        stats ['- CovQual %s' % k ] = covs [ k ]
            
        
    if progress or not silent : print_stats ( stats , nToys , logger = logger )
    
    return results, stats 


# =============================================================================
## make <code>nToys</code> pseudoexperiments
#
#  Schematically:
#  @code
#  for toy in range ( nToys )  :
#  ...  dataset = gen_fun ( gen_pdf , ...     , **gen_config )
#  ...  result  = fit_fun ( fit_pdf , dataset , **fit_config )
#  ...  if not accept_fun ( result  , fit_pdf , dataset ) : continue
#  .... < collect statistics here > 
#  @endcode
#
#  For each experiment
#  - generate dataset using <code>pdf</code> with variables specified
#    in <code>data</code> and configuration specified via<code>gen_config</code>
#    for each generation the parameters of <code>pdf</code> are reset
#    for their initial values and values from <code>init_pars</code>
#  - fit generated dataset  with <code>pdf</code> using configuration
#    specified via  <code>fit_config</code>
#
# @code
# gen_pdf = ... ## PDF  to use to generate pseudoexperiments 
# fit_pdf = ... ## PDF  to use to fit  pseudoexperiments 
# results , stats = make_toys2(
#     gen_pdf    = gen_pdf    , ## PDF  to use to generate pseudoexperiments 
#     fit_pdf    = fit_pdf    , ## PDF  to use to fit  pseudoexperiments 
#     nToys      = 1000       , ## number of pseudoexperiments 
#     data       = [ 'mass' ] ,           ## variables in dataset 
#     gen_config = { 'nEvents' : 5000 } , ## configuration of <code>pdf.generate</code>
#     fit_config = { 'ncpus'   : 2    } , ## configuration of <code>pdf.fitTo</code>
#     gen_pars   = { 'mean' : 0.0 , 'sigma' : 1.0 } ## parameters to use for generation 
#     )
# @endcode
#
# Derived parameters can be also   retrived via <code>more_vars</code> argument:
# @code
# ratio     = lambda res,pdf : res.ratio('x','y')
# more_vars = { 'Ratio' : ratio }
#  r,  s = make_toys2 ( .... , more_vars = more_vars , ... ) 
# @code
#
# @param gen_pdf    PDF to be used for generation 
# @param fit_pdf    PDF to be used for fitting
# @param nToys      number    of pseudoexperiments to generate
# @param data       variable list of variables to be used for dataset generation
# @param gen_config configuration of <code>pdf.generate</code>
# @param fit_config configuration of <code>pdf.fitTo</code>
# @param gen_pars   redefine these parameters for each pseudoexperiment
# @param fit_pars   redefine these parameters for each pseudoexperiment
# @param more_vars  calculate more variables form fit-result
# @param gen_fun    specific generate  action (if needed) 
# @param fit_fun    specific fitting action (if needed) 
# @param accept_fun specific accept action (if needed) 
# @param silent     silent toys?
# @param progress   show progress bar?
# @param logger     logger 
# @return dictionary with fit results for the toys and the dictionary of statistics
#
#  - If <code>gen_fun</code>    is not specified <code>generate_data</code> is used 
#  - If <code>fit_fun</code>    is not specified <code>make_fit</code>      is used 
#  - If <code>accept_fun</code> is not specified <code>accept_fit</code>    is used   
def make_toys2 ( gen_pdf             , ## pdf to generate toys 
                 fit_pdf             , ## pdf to fit  
                 nToys               , ## number of pseudoexperiments 
                 data                , ## template for dataset/variables 
                 gen_config          , ## parameters for <code>pdf.generate</code>   
                 fit_config = {}     , ## parameters for <code>pdf.fitTo</code>
                 gen_pars   = {}     , ## gen-parameters to reset/use 
                 fit_pars   = {}     , ## fit-parameters to reset/use
                 more_vars  = {}     , ## additional  results to be calculated
                 gen_fun    = None   , ## generator function ( pdf , varset  , **gen_config ) 
                 fit_fun    = None   , ## fit       function ( pdf , dataset , **fit_config ) 
                 accept_fun = None   , ## accept    function ( fit-result, pdf, dataset     )
                 silent     = True   ,
                 progress   = True   ,
                 logger     = logger ) :
    """Make `ntoys` pseudoexperiments
    
    -   Schematically:
    >>> for toy in range ( nToys )  :
    >>> ...  dataset = gen_fun ( gen_pdf , ...     , **gen_config )
    >>> ...  result  = fit_fun ( fit_pdf , dataset , **fit_config )
    >>> ...  if not accept_fun ( result  , fit_pdf , dataset ) : continue
    >>> .... < collect statistics here > 
    
    For each experiment:

    1. generate dataset using `pdf` with variables specified
    in `data` and configuration specified via `gen_config`
    for each generation the parameters of `pdf` are reset
    for their initial values and valeus from `init_pars`
    
    2. fit generated dataset  with `pdf` using configuration
    specified via  `fit_config`

    - pdf        PDF to be used for generation and fitting
    - nToys      number    of pseudoexperiments to generate
    - data       variable list of variables to be used for dataset generation
    - gen_config configuration of <code>pdf.generate</code>
    - fit_config configuration of <code>pdf.fitTo</code>
    - gen_pars   redefine these parameters for generation of each pseudoexperiment
    - fit_pars   redefine these parameters for fit of each pseudoexperiment
    - silent     silent toys?
    - progress  show progress bar? 
    
    It returns a dictionary with fit results for the toys and a dictionary of statistics
    >>> pdf = ...
    ... results, stats = make_toys ( pdf     , ## PDF  to use 
    ...                 1000                 , ## number of toys 
    ...                 [ 'mass' ]           , ## varibales in dataset 
    ...                 { 'nEvents' : 5000 } , ## configuration of `pdf.generate`
    ...                 { 'ncpus'   : 2    } , ## configuration of `pdf.fitTo`
    ...                 { 'mean' : 0.0 , 'sigma' : 1.0 } ## parameters to use for generation 
    ...                )
    """

    from ostap.core.ostap_types import string_types, integer_types  
    
    assert isinstance ( nToys , integer_types ) and 0 < nToys,\
           'Invalid "nToys" argument %s/%s' % ( nToys , type ( nToys ) )
    
    assert gen_config and 'nEvents' in gen_config,\
           'Number of events per toy must be specified via "gen_config" %s' % gen_config
    
    ## 1. generator function? 
    if gen_fun is None :
        if not silent :  logger.info ( "make_toys2: use default ``generate_data'' function!")
        gen_fun = generate_data 
    assert gen_fun and callable ( gen_fun ) , 'Invalid generator function!'
    
    ## 2. fitting function? 
    if fit_fun is None :
        if not silent :  logger.info ( "make_toys2: use default ``make_fit'' function!")
        fit_fun = make_fit 
    assert fit_fun and callable ( fit_fun ) , 'Invalid fit function!'

    ## 3. accept function? 
    if accept_fun is None :
        if not silent : logger.info ( "make_toys2: use default ``accept_fit'' function!")
        accept_fun = accept_fit
    assert accept_fun and callable ( accept_fun ) , 'Invalid accept function!'

    
    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult
    import ostap.fitting.basic 

    gparams = gen_pdf.params ()
    varset  = ROOT.RooArgSet () 
    
    if isinstance ( data , ROOT.RooAbsData ) : varset = data.varset() 
    else :
        for v in data :
            if   isinstance ( v , ROOT.RooAbsArg ) :
                varset.add ( v )
            elif isinstance ( v , string_types   ) and v in gparams :
                varset.add ( gparams [ v ] )
            else :
                raise TypeError('Invalid variable %s/%s' % ( v , type ( v ) ) )

    ## parameters for generation
            
    fix_gen_init = vars_transform ( gparams  ) 
    fix_gen_pars = vars_transform ( gen_pars )
    
    ## parameters for fitting 

    fparams = fit_pdf.params ()
    fix_fit_init = vars_transform ( fparams  )     
    fix_fit_pars = vars_transform ( fit_pars )
    
    fitcnf = {}
    fitcnf.update ( fit_config )
    if not 'silent' in fitcnf : fitcnf [ 'silent' ] = silent

    from collections import defaultdict 
    results = defaultdict(list) 

    from   ostap.core.core        import SE
    
    fits = defaultdict ( SE )  ## fit statuses 
    covs = defaultdict ( SE )  ## covarinace matrix quality

    ## run pseudoexperiments
    from ostap.utils.progress_bar import progress_bar 
    for i in progress_bar ( range ( nToys ) , silent = not progress ) :

        ## 1. reset PDF parameters 
        gen_pdf.load_params ( params = fix_gen_init , silent = silent )
        gen_pdf.load_params ( params = fix_gen_pars , silent = silent )

        ## 2. generate dataset!
        dataset =  gen_fun ( gen_pdf , varset = varset , **gen_config ) 
        if not silent : logger.info ( 'Generated dataset #%d\n%s' % ( i , dataset ) )

        ## 3. reset parameters of fit_pdf
        fit_pdf.load_params ( params = fix_fit_init , silent = silent )
        fit_pdf.load_params ( params = fix_fit_pars , silent = silent )
        
        ## 4. fit it!  
        r = fit_fun ( fit_pdf , dataset , **fitcnf ) 

        ## fit status 
        fits [ r.status  () ] += 1

        ## covariance matrix quality
        covs [ r.covQual () ] += 1

        ## ok ?
        if accept_fun ( r , fit_pdf , dataset ) : 

            ## 5. save results 
            rpf = r.params ( float_only = True ) 
            for i in rpf : 
                results [ i ].append ( rpf[i][0] ) 
                
            for v in more_vars :
                func  = more_vars[v] 
                results [ v ] .append ( func ( r , fit_pdf ) )
                
            results [ '#' ] .append ( len ( dataset ) )

        dataset.clear()
        del dataset
        
    ## make a final statistics 
    stats = defaultdict ( SE )
    
    for par in results :
        pars = results [ par ]
        for v in pars : 
            v0 = float ( v )         
            stats [ par ] += v0 
            
    for k in fits :
        stats ['- Status  %s' % k ] = fits [ k ]
    for k in covs :
        stats ['- CovQual %s' % k ] = covs [ k ]
                    
    if progress or not silent :
        print_stats ( stats , nToys , logger = logger  )

    return results, stats 

# =============================================================================
## run Jackknife analysis, useful for evaluaton of fit biased and uncertainty estimates
# 
#  For each <code>i</code> remove event with index <code>i</code> from the dataset,
#  and refit it.
#  @code
#  dataset = ...
#  model   = ...
#  r , f = model.fitTo ( dataset , .... )           ## fit the whole dataset   
#  results, stats = make_jackknife ( model , data ) ## run Jackknife 
#  print_jackknife ( r , stats )                    ## print summary table 
#  @endcode
#  @see printJackknife
#  @see https://en.wikipedia.org/wiki/Jackknife_resampling
#  @param pdf         fit model
#  @param data        original dataset
#  @param fit_config  configuration of <code>pdf.FitTo( data , ... )</code>
#  @param fit_pars    redefine these parameters before each fit 
#  @param more_vars   calculate more variables from the fit-results 
#  @param fit_fun     fitting   function
#  @param accept_fun  accept    function
#  @param event_range event range to use for jackknife   
#  @param silent      silent processing 
#  @param progress    show progress bar?
#  @param logger      use this logger
#  @return statistics of jackknife experiments 
def make_jackknife ( pdf                  ,
                     data                 ,
                     fit_config  = {}     , ## parameters for <code>pdf.fitTo</code>
                     fit_pars    = {}     , ## fit-parameters to reset/use
                     more_vars   = {}     , ## additional  results to be calculated
                     fit_fun     = None   , ## fit       function ( pdf , dataset , **fit_config ) 
                     accept_fun  = None   , ## accept    function ( fit-result, pdf, dataset     )
                     event_range = ()     , ## event range for jackknife                      
                     silent      = True   ,
                     progress    = True   ,
                     logger      = logger ) :
    """Run Jackknife analysis, useful for evaluaton of fit biased and uncertainty estimates
    For each <code>i</code> remove event with index <code>i</code> from the dataset, and refit it.
    >>> dataset = ...
    >>> model   = ...
    >>> r , f = model.fitTo ( dataset , .... )           ## fit the whole dataset   
    >>> results, stats = make_jackknife ( model , data ) ## run Jackknife 
    >>> print_jackknife ( r , stats )                    ## print summary table 
    - see https://en.wikipedia.org/wiki/Jackknife_resampling
    - see print_jackknife 
    - see jackknife_statistics


    - `pdf`         : fit model
    - `data`        : original dataset
    - `fit_config`  : configuration of `pdf.FitTo( data , ... )`
    - `fit_pars`    : redefine these parameters before each fit
    - `more_vars`   : calculate more variables from the fit-results
    - `fit_fun`     : specific fitting acion (if needed) 
    - `accept_fun`  : specific accept action (if needed)
    - `event_range` : event range to use for jackknife   
    - `silent`      : silent processing?
    - `progress`    : show progress bar?
    - `logger`      : use this logger 

    """
    
    N = len ( data )
    assert 1 < N            , 'make_jackknife: invalid dataset size %s' % N

    if not event_range : event_range = 0 , N 
    assert 2 == len ( event_range ) , 'make_jackknife: invalid event range %s ' % str ( event_range )
    
    begin , end = event_range
    assert 0 <= begin and begin < end and end <= N, 'make_jackknife: invalid event range (%s,%s)/%d' % ( begin , end , N )
    
    ## 1. fitting function? 
    if fit_fun is None :
        if not silent :  logger.info ( "make_jackknife: use default ``make_fit'' function!")
        fit_fun = make_fit 
    assert fit_fun and callable ( fit_fun ) , 'Invalid fit function!'

    ## 2. accept function? 
    if accept_fun is None :
        if not silent : logger.info ( "make_jackknife: use default ``accept_fit'' function!")
        accept_fun = accept_fit
    assert accept_fun and callable ( accept_fun ) , 'Invalid accept function!'

    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult
    import ostap.fitting.basic 
       
    ## parameters for fitting 

    fparams      = pdf.params ()
    fix_fit_init = vars_transform ( fparams  )     
    fix_fit_pars = vars_transform ( fit_pars )
    
    fitcnf = {}
    fitcnf.update ( fit_config )
    if not 'silent' in fitcnf : fitcnf [ 'silent' ] = silent
    
    from collections import defaultdict 
    results = defaultdict(list) 

    from   ostap.core.core        import SE    
    fits = defaultdict ( SE )  ## fit statuses 
    covs = defaultdict ( SE )  ## covarinace matrix quality

    ## Fit the whole sample 
    pdf.load_params ( params = fix_fit_init , silent = silent )
    pdf.load_params ( params = fix_fit_pars , silent = silent )
    r_tot = fit_fun ( pdf , data , **fitcnf )
    
    from ostap.utils.progress_bar import progress_bar
    ## run jackknife  bootstrapping
    for ds in progress_bar ( data.jackknife ( begin , end ) , max_value = end - begin , silent = not progress ) :

        ## 2. reset parameters of fit_pdf
        pdf.load_params ( params = fix_fit_init , silent = silent )
        pdf.load_params ( params = fix_fit_pars , silent = silent )
 
        ## 3. fit it!  
        r = fit_fun ( pdf , ds , **fitcnf ) 

        ## 4. fit status 
        fits [ r.status  () ] += 1

        ## 5. covariance matrix quality
        covs [ r.covQual () ] += 1

        ## ok ?
        if accept_fun ( r , pdf , ds ) : 

            ## 6. save results 
            rpf = r.params ( float_only = True ) 
            for i in rpf : 
                results [ i ].append ( rpf[i][0] ) 

            ## 7. more variables to be calculated? 
            for v in more_vars :
                func  = more_vars[v] 
                results [ v ] .append ( func ( r , pdf ) )
                
            results [ '#' ] .append ( len ( ds ) )

        ds.clear()
        
    ## 8. make a final statistics 
    stats = defaultdict ( SE )
    
    for par in results :
        pars = results [ par ]
        for v in pars : 
            v0 = float ( v )         
            stats [ par ] += v0 
            
    for k in fits :
        stats ['- Status  %s' % k ] = fits [ k ]
    for k in covs :
        stats ['- CovQual %s' % k ] = covs [ k ]
        
    if progress or not silent :

        ## 9. fit total dataset (twice) 
        r_tot = fit_fun ( pdf , data , **fitcnf )
        r_tot = fit_fun ( pdf , data , **fitcnf )
        
        ## 10. the final table  
        print_jackknife ( r_tot , stats , logger = logger )
            
    return results , stats 


# =============================================================================
## Run Bootstrap analysis, useful for evaluaton of fit biased and uncertainty estimates
# 
#  In total <code>size</code> datasets are sampled (with replacement) from the orifinal dataste
#  <code>data</code> and each sampled dataset is fit
#  @code
#  dataset = ...
#  model   = ...
#  r , f = model.fitTo ( dataset , .... )                         ## fit the whole dataset   
#  results, stats = make_bootstrap ( model , data , size = 1000 ) ## run Bootstrap 
#  print_bootstrap ( r , stats )                    ## print summary table 
#  @endcode
#  @see print_bootstrap
#  @param pdf   fit model
#  @param data  original dataset
#  @param size  number of datasets to sample
#  @param fit_config configuration of <code>pdf.FitTo( data , ... )</code>
#  @param fit_pars   redefine these parameters before each fit 
#  @param more_vars  calculate more variables from the fit-results 
#  @param fit_fun    specific fitting action (if needed) 
#  @param accept_fun specific accept action (if needed) 
#  @param silent     silent processing 
#  @param progress   show progress bar?
#  @param logger     use this logger
#  @return statistics of boostrap experiments 
def make_bootstrap ( pdf                  ,
                     data                 ,
                     size        = 100    ,   ## numbere of samples 
                     fit_config  = {}     ,   ## parameters for <code>pdf.fitTo</code>
                     fit_pars    = {}     ,   ## fit-parameters to reset/use
                     more_vars   = {}     ,   ## additional  results to be calculated
                     fit_fun     = None   ,   ## fit       function ( pdf , dataset , **fit_config ) 
                     accept_fun  = None   ,   ## accept    function ( fit-result, pdf, dataset     )
                     silent      = True   ,   ## silent processing?
                     progress    = True   ,   ## shpow progress bar? 
                     logger      = logger ) : ## use this logger 

    """Run Bootstrap analysis, useful for evaluaton of fit biased and uncertainty estimates 
    In total `size` datasets are sampled (with replacement) from the orifinal dataste
    `data` and each sampled dataset is fit
    >>> dataset = ...
    >>> model   = ...
    >>> r , f = model.fitTo ( dataset , .... )                         ## fit the whole dataset   
    >>> results, stats = make_bootstrap ( model , data , size = 1000 ) ## run Bootstrap 
    >>> print_bootstrap ( r , stats )                    ## print summary table 

    - `pdf`        : fit model
    - `data`       : original dataset
    - `size`       : number of datasets to sample
    - `fit_config` : configuration of `pdf.FitTo( data , ... )`
    - `fit_pars`   : redefine these parameters before each fit
    - `more_vars`  : calculate more variables from the fit-results
    - `fit_fun`    : specific fitting acion (if needed) 
    - `accept_fun` : specific accept action (if needed) 
    - `silent`     : silent processing?
    - `progress`   : show progress bar?
    - `logger`     : use this logger 
    """
    
    N = len ( data )
    assert 1 < N            , 'make_bootstrap: invalid dataset size %s' % N

    from ostap.core.ostap_types import integer_types  
    assert isinstance ( size , integer_types ) and 0 < size, \
           "make_bootstrap: invalid ``size'' parameter %s" % size 
    
    ## 1. fitting function? 
    if fit_fun is None :
        if not silent :  logger.info ( "make_bootstrap: use default ``make_fit'' function!")
        fit_fun = make_fit 
    assert fit_fun and callable ( fit_fun ) , 'Invalid fit function!'

    ## 2. accept function? 
    if accept_fun is None :
        if not silent : logger.info ( "make_bootstrap: use default ``accept_fit'' function!")
        accept_fun = accept_fit
    assert accept_fun and callable ( accept_fun ) , 'Invalid accept function!'

    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult
    import ostap.fitting.basic 
       
    ## parameters for fitting 

    fparams      = pdf.params ()
    fix_fit_init = vars_transform ( fparams  )     
    fix_fit_pars = vars_transform ( fit_pars )
    
    fitcnf = {}
    fitcnf.update ( fit_config )
    if not 'silent' in fitcnf : fitcnf [ 'silent' ] = silent
    
    from collections import defaultdict 
    results = defaultdict(list) 

    from   ostap.core.core        import SE    
    fits = defaultdict ( SE )  ## fit statuses 
    covs = defaultdict ( SE )  ## covarinace matrix quality

    ## fit original dataset 
    pdf.load_params ( params = fix_fit_init , silent = silent )
    pdf.load_params ( params = fix_fit_pars , silent = silent )
    r_tot = fit_fun ( pdf , data , **fitcnf )

    from ostap.utils.progress_bar import progress_bar
    ## run jackknife  bootstrapping
    for ds in progress_bar ( data.bootstrap ( size ) , max_value = size , silent = not progress ) :

        ## 2. reset parameters of fit_pdf
        pdf.load_params ( params = fix_fit_init , silent = silent )
        pdf.load_params ( params = fix_fit_pars , silent = silent )
 
        ## 3. fit it!  
        r = fit_fun ( pdf , ds , **fitcnf ) 

        ## 4. fit status 
        fits [ r.status  () ] += 1

        ## 5. covariance matrix quality
        covs [ r.covQual () ] += 1

        ## ok ?
        if accept_fun ( r , pdf , ds ) : 

            ## 6. save results 
            rpf = r.params ( float_only = True ) 
            for i in rpf : 
                results [ i ].append ( rpf[i][0] ) 

            ## 7. more variables to be calculated? 
            for v in more_vars :
                func  = more_vars[v] 
                results [ v ] .append ( func ( r , pdf ) )
                
            results [ '#' ] .append ( len ( ds ) )

        ds.clear()
        
    ## 8. make a final statistics 
    stats = defaultdict ( SE )
    
    for par in results :
        pars = results [ par ]
        for v in pars : 
            v0 = float ( v )         
            stats [ par ] += v0 
            
    for k in fits :
        stats ['- Status  %s' % k ] = fits [ k ]
    for k in covs :
        stats ['- CovQual %s' % k ] = covs [ k ]
        
    if progress or not silent :

        ## 9. fit total dataset (twice) 
        r_tot = fit_fun ( pdf , data , **fitcnf )
        r_tot = fit_fun ( pdf , data , **fitcnf )
        
        ## 10. the final table  
        print_bootstrap ( r_tot , stats , logger = logger )
            
    return results , stats 


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
