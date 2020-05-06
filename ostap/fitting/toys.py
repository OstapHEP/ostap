#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file ostap/fitting/toys.py
#  Simple utilities to run fitting toys   
#  @date   2020-01-18
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
#  - thanks to Albert PUIG
#
# =============================================================================
""" Simple utilities to run fitting toys
Python interface to two major TMVA classes
"""
# =============================================================================
__author__  = 'Vanya BELYAEV  Ivan.Belyaev@itep.ru'
__date__    = "2020-01-18"
__version__ = '$Revision$'
__all__     = (
    "make_toys"      , ## run fitting toys (the same PDF to generate and fit)
    "make_toys2"     , ## run fitting toys (separate models to generate and fit)
    "vars_transform" , ## helper fnuction to transform the variables
    "print_stats"    , ## print statistics of toys 
    )
# =============================================================================
from   builtins          import range
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.toys' )
else                       : logger = getLogger( __name__             )
# =============================================================================
logger.debug ( 'Utilities to run fitting toys')
# ==============================================================================
## Technical transformation  to the dictionaty :  { 'name' : float_value }
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
            result [ p.GetName()  ] = float ( p )
            
    return result
# =============================================================================
## print statistics of pseudoexperiments
def print_stats (  stats , ntoys = '???' ) :
    """print statistics of pseudoexperiments
    """
    
    table = [ ( '' , '#', 'mean' , 'rms' , '%11s / %-11s' % ( 'min' , 'max' ) ) ] 
    keys = stats.keys()
    keys = sorted ( keys )

    def make_row ( c ) :
        n      = "{:^11}".format ( c.nEntries() )
        mean   = c.mean ()
        mean   = "%+11.4g +- %-11.4g" % ( mean.value() , mean.error() )
        rms    = "%-11.4g"            % c.rms ()
        minmax = "%+11.4g / %-+11.4g" % ( c.min() , c.max () ) 
        return p , n , mean , rms  , minmax 
        
    
    for p in sorted ( stats )  :
        if    p.startswith('pull:') : continue  
        c      = stats [ p ]
        table.append (  make_row ( c )  )
        
    for p in sorted ( stats )  :
        if not p.startswith('pull:') : continue  
        c      = stats [ p ]
        table.append (  make_row ( c )  )
        
    import ostap.logger.table as Table
    table = Table.table ( table , title = "Results of %s toys" % ntoys , prefix = "# " )
    logger.info ( 'Results of %s toys:\n%s' % ( ntoys , table ) ) 
    
# ==============================================================================
## make <code>ntoys</code> pseudoexperiments
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
# Derived parameters can be also   retrived via <code>more_vars</code> argument:
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
# @param silent     silent toys?
# @return dictionary with fit results for the toys and the dictionary of statistics  
def make_toys ( pdf                ,
                nToys              , 
                data               , ## template for dataset/variables 
                gen_config         , ## parameters for <code>pdf.generate</code>   
                fit_config = {}    , ## parameters for <code>pdf.fitTo</code>
                init_pars  = {}    ,
                more_vars  = {}    , 
                silent     = True  ,
                progress   = True  ) :
    """Make `ntoys` pseudoexperiments
    
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
    - init_pars  redefine these parameters for each pseudoexperiment
    - more_vars  dictionary of functions to define the additional results 
    - silent     silent toys?
    - progress   show progress bar? 
    
    It returns a dictionary with fit results for the toys
    
    >>> pdf = ...
    ... results, stats = make_toys ( pdf     , ## PDF  to use 
    ...                 1000                 , ## number of toys 
    ...                 [ 'mass' ]           , ## varibales in dataset 
    ...                 { 'nEvents' : 5000 } , ## configuration of `pdf.generate`
    ...                 { 'ncpus'   : 2    } , ## configuration of `pdf.fitTo`
    ...                 { 'mean' : 0.0 , 'sigma' : 1.0 } ## parameters to use for generation 
    ...                )

    Derived parameters can be also retrived via <code>more_vars</code> argument:
    >>> ratio    = lambda res,pdf : res.ratio('x','y') 
    >>> more_vars = { 'Ratio' : ratio }
    >>> r,  s = make_toys ( .... , more_vars = more_vars , ... ) 

    """

    from ostap.core.ostap_types import string_types, integer_types  
    
    assert isinstance ( nToys , integer_types ) and 0 < nToys,\
           'Invalid "nToys" argument %s/%s' % ( nToys , type ( nToys ) )
    
    assert gen_config and 'nEvents' in gen_config,\
           'Number of events per toy must be specified via "gen_config" %s' % gen_config
    
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

    pdf.load_params ( None , fix_pars , silent = silent )
    pdf.load_params ( None , fix_init , silent = silent )

    ## save all initial parameters (needed fot the final statistics)
    params  = pdf.params      ()
    fix_all = vars_transform  ( params ) 
    
    fitcnf = {}
    fitcnf.update ( fit_config )
    if not 'silent' in fitcnf : fitcnf [ 'silent' ] = silent

    from collections import defaultdict 
    results = defaultdict(list) 

    ## run pseudoexperiments
    from ostap.utils.progress_bar import progress_bar 
    for i in progress_bar ( range ( nToys ) , silent = not progress ) :

        ## 1. reset PDF parameters 
        pdf.load_params ( None , fix_pars  , silent = silent )
        pdf.load_params ( None , init_pars , silent = silent )

        ## 2. generate dataset!  
        dataset = pdf.generate ( varset = varset , **gen_config )  
        if not silent :
            logger.info ( 'Generated dataset #%d\n%s' % ( i , dataset ) )
        
        ## 3. fit it!  
        r , _ = pdf.fitTo ( dataset , **fitcnf )
        if not silent :
            logger.info ( 'Fit result #%d\n%s' % ( i , r.table ( title = 'Fit result #%d' % i , prefix = '# ' ) ) )

        ## ok ? 
        if r and 0 == r.status () :
            
            ## 4. save results 
            rpf = r.params ( float_only = True ) 
            for i in rpf : 
                results [ i ].append ( rpf[i][0] ) 
                
            for v in more_vars :
                func  = more_vars[v] 
                results [v] .append ( func ( r , pdf ) )
                
        dataset.clear()
        del dataset
        del r

        
    ## make a final statistics 
    from   ostap.core.core        import SE 
    stats = defaultdict(SE)
    
    for par in results :
        pars = results [ par ]
        mvar = par in more_vars 
        if not mvar : a0 = fix_all [ par ]
        for v in pars : 
            v0 = float ( v )         
            stats     [ par             ] +=   v0
            if not mvar : 
                stats [ 'pull:%s' % par ] += ( v0 - a0 ) / v.error()

    if progress or not silent : print_stats ( stats , nToys )
    
    return results, stats 

# =============================================================================
## make <code>ntoys</code> pseudoexperiments
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
# @param silent     silent toys?
# @return dictionary with fit results for the toys and the dictionary of statistics  
def make_toys2 ( gen_pdf            , ## pdf to generate toys 
                 fit_pdf            , ## pdf to fit  
                 nToys              , ## number of pseudoexperiments 
                 data               , ## template for dataset/variables 
                 gen_config         , ## parameters for <code>pdf.generate</code>   
                 fit_config = {}    , ## parameters for <code>pdf.fitTo</code>
                 gen_pars   = {}    , ## gen-parameters to reset/use 
                 fit_pars   = {}    , ## fit-parameters to reset/use
                 more_vars  = {}    , ## additional  results to be calculated  
                 silent     = True  ,
                 progress   = True  ) :
    """Make `ntoys` pseudoexperiments
    
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
    
    It returns a dictionary with fit results for the toys
    
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

    ## run pseudoexperiments
    from ostap.utils.progress_bar import progress_bar 
    for i in progress_bar ( range ( nToys ) , silent = not progress ) :

        ## 1. reset PDF parameters 
        gen_pdf.load_params ( None , fix_gen_init , silent = silent )
        gen_pdf.load_params ( None , fix_gen_pars , silent = silent )

        ## 2. generate dataset!
        dataset = gen_pdf.generate ( varset = varset , **gen_config )  
        if not silent :
            logger.info ( 'Generated dataset #%d\n%s' % ( i , dataset ) )

        ## 3. reset parameters of fit_pdf
        fit_pdf.load_params ( None , fix_fit_init , silent = silent )
        fit_pdf.load_params ( None , fix_fit_pars , silent = silent )
        
        ## 4. fit it!  
        r , _ = fit_pdf.fitTo ( dataset , **fitcnf )
        if not silent :
            logger.info ( 'Fit result #%d\n%s' % ( i , r.table ( title = 'Fit result #%d' % i , prefix = '# ' ) ) )

        ## skip invalid fits 
        if r.status () : continue

        ## 5. save results 
        rpf = r.params ( float_only = True ) 
        for i in rpf : 
            results [ i ].append ( rpf[i][0] ) 

        for v in more_vars :
            func  = more_vars[v] 
            results [ v ] .append ( func ( r , fit_pdf ) )

        dataset.clear()
        del dataset
        
    ## make a final statistics 
    from   ostap.core.core        import SE 
    stats = defaultdict(SE)
    
    for par in results :
        pars = results [ par ]
        for v in pars : 
            v0 = float ( v )         
            stats [ par ] += v0 
            
    if progress or not silent : print_stats ( stats , nToys )

    return results, stats 


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
