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
    "make_toys"  , ## run fitting toys (the same PDF to generate and fit)
    "make_toys2" , ## run fitting toys (separate models to generate and fit)
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
## make <code>ntoys</code> pseudoexperiments
#
#  For each experiment
#  - generate dataset using <code>pdf</code> with variables specified
#    in <code>data</code> and configuration specified via<code>toy_config</code>
#    for each generation the parameters of <code>pdf</code> are reset
#    for their initial values and values from <code>init_pars</code>
#  - fit generated dataset  with <code>pdf</code> using configuration
#    specified via  <code>fit_config</code>
#
# @code
# pdf = ...
# results , stats = make_toys ( pdf    ,           ## PDF  to use 
#                 1000                 ,           ## number of toys 
#                 [ 'mass' ]           ,           ## varianles in dataset 
#                 { 'nEvents' : 5000 } ,           ## configuration of <code>pdf.generate</code>
#                 { 'ncpus'   : 2    } ,           ## configuration of <code>pdf.fitTo</code>
#                 { 'mean' : 0.0 , 'sigma' : 1.0 } ## parameters to use for generation 
#                )
# @endcode
#
# Derived parameters can be also   retrived via <code>more_vars</code> argument:
# @code
# more_vars = { 'Ratio' : ( lambda x,y : x/y , ('x','y') ) }
#  r,  s = make_toys ( .... , more_vars = more_vars , ... ) 
# @code
#
# @param pdf        PDF to be used for generation and fitting
# @param nToys      number    of pseudoexperiments to generate
# @param data       variable list of variables to be used for dataset generation
# @param toy_config configuration of <code>pdf.generate</code>
# @param fit_config configuration of <code>pdf.fitTo</code>
# @param init_pars  redefine these parameters for each pseudoexperiment
# @param more_vars  calculate more variables form fit-result 
# @param silent     silent toys?
# @return dictionary with fit results for the toys and the dictionary of statistics  
def make_toys ( pdf                ,
                nToys              , 
                data               , ## template for dataset/variables 
                toy_config         , ## parameters for <code>pdf.generate</code>   
                fit_config = {}    , ## parameters for <code>pdf.fitTo</code>
                init_pars  = {}    ,
                more_vars  = {}    , 
                silent     = True  ,
                progress   = False ) :
    """Make `ntoys` pseudoexperiments
    
    For each experiment:

    1. generate dataset using `pdf` with variables specified
    in `data` and configuration specified via `toy_config`
    for each generation the parameters of `pdf` are reset
    for their initial values and valeus from `init_pars`
    
    2. fit generated dataset  with `pdf` using configuration
    specified via  `fit_config`

    - pdf        PDF to be used for generation and fitting
    - nToys      number    of pseudoexperiments to generate
    - data       variable list of variables to be used for dataset generation
    - toy_config configuration of <code>pdf.generate</code>
    - fit_config configuration of <code>pdf.fitTo</code>
    - init_pars  redefine these parameters for each pseudoexperiment
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
    >>> more_vars = { 'Ratio' : ( lambda x,y : x/y , ('x','y') ) }
    >>> r,  s = make_toys ( .... , more_vars = more_vars , ... ) 

    """

    from ostap.core.ostap_types import string_types, integer_types  
    
    assert isinstance ( nToys , integer_types ) and 0 < nToys,\
           'Invalid "nToys" argument %s/%s' % ( nToys , type ( nToys ) )
    
    assert toy_config and 'nEvents' in toy_config,\
           'Number of events per toy must be specified via "toy_config" %s' % toy_config
    
    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult
    import ostap.fitting.basic 

    params = pdf.pdf.getParameters ( None )
    varset = ROOT.RooArgSet() 
    
    if   isinstance ( data , ROOT.RooArgSet  ) : pass 
    elif isinstance ( data , ROOT.RooAbsData ) : varset = data.varset() 
    else :
        for v in data :
            if   isinstance ( v , ROOT.RooAbsArg ) :
                varset.add ( v )
            elif isinstance ( v , string_types   ) and v in params :
                varset.add ( params [ v ] )
            else :
                raise TypeError('Invalid variable %s/%s' % ( v , type ( v ) ) )

    fix_pars = {}
    for v in params : fix_pars [ v.name ] = float ( v ) 

    fix_init = {}
    if   isinstance ( init_pars , ROOT.RooFitResult ) :
        ps = init_pars.dct_params()
        for p in ps        : fix_init [ p      ] = float ( ps        [ p ] )
    elif isinstance ( init_pars , dict ) :
        for p in init_pars : fix_init [ p      ] = float ( init_pars [ p ] )
    else :
        for p in init_pars : fix_init [ p.name ] = float ( p )

    pdf.load_params ( None , fix_pars , silent = silent )
    pdf.load_params ( None , fix_init , silent = silent )

    ## save all initial parameters (needed fot the final statistics)
    params  = pdf.pdf.getParameters ( None )
    fix_all = {}
    for p in params : fix_all [ p.name ] = float ( p )
    
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
        dataset = pdf.generate ( varset = varset , **toy_config )  
        if not silent :
            logger.info ( 'Generated dataset #%d\n%s' % ( i , dataset ) )
        
        ## 3. fit it!  
        r , _ = pdf.fitTo ( dataset , **fitcnf )
        if not silent :
            logger.info ( 'Fit result #%d\n%s' % ( i , r.table ( title = 'Fit result #%d' % i , prefix = '# ' ) ) )

        ## 
        if r.status () : continue

        ## 4. save results 
        rpf = r.params ( float_only = True ) 
        for i in rpf : 
            results [ i ].append ( rpf[i][0] ) 

        for v in more_vars :
            results [v] .append ( r.evaluate ( *more_vars[v] ) )
            
        dataset.clear()
        del dataset
        
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
                
    return results, stats 

# =============================================================================
## make <code>ntoys</code> pseudoexperiments
#
#  For each experiment
#  - generate dataset using <code>pdf</code> with variables specified
#    in <code>data</code> and configuration specified via<code>toy_config</code>
#    for each generation the parameters of <code>pdf</code> are reset
#    for their initial values and values from <code>init_pars</code>
#  - fit generated dataset  with <code>pdf</code> using configuration
#    specified via  <code>fit_config</code>
#
# @code
# pdf = ...
# results , stats = make_toys ( pdf    ,           ## PDF  to use 
#                 1000                 ,           ## number of toys 
#                 [ 'mass' ]           ,           ## varianles in dataset 
#                 { 'nEvents' : 5000 } ,           ## configuration of <code>pdf.generate</code>
#                 { 'ncpus'   : 2    } ,           ## configuration of <code>pdf.fitTo</code>
#                 { 'mean' : 0.0 , 'sigma' : 1.0 } ## parameters to use for generation 
#                )
# @endcode
#
# @param gen_pdf    PDF to be used for generation 
# @param fit_pdf    PDF to be used for fitting
# @param nToys      number    of pseudoexperiments to generate
# @param data       variable list of variables to be used for dataset generation
# @param toy_config configuration of <code>pdf.generate</code>
# @param fit_config configuration of <code>pdf.fitTo</code>
# @param gen_pars   redefine these parameters for each pseudoexperiment
# @param fit_pars   redefine these parameters for each pseudoexperiment
# @param silent     silent toys?
# @return dictionary with fit results for the toys and the dictionary of statistics  
def make_toys2 ( gen_pdf            , ## pdf to generate toys 
                 fit_pdf            , ## pdf to fit  
                 nToys              , ## number of pseudoexperiments 
                 data               , ## template for dataset/variables 
                 toy_config         , ## parameters for <code>pdf.generate</code>   
                 fit_config = {}    , ## parameters for <code>pdf.fitTo</code>
                 gen_pars   = {}    , ## gen-parameters to reset/use 
                 fit_pars   = {}    , ## fit-parameters to reset/use 
                 silent     = True  ,
                 progress   = False ) :
    """Make `ntoys` pseudoexperiments
    
    For each experiment:

    1. generate dataset using `pdf` with variables specified
    in `data` and configuration specified via `toy_config`
    for each generation the parameters of `pdf` are reset
    for their initial values and valeus from `init_pars`
    
    2. fit generated dataset  with `pdf` using configuration
    specified via  `fit_config`

    - pdf        PDF to be used for generation and fitting
    - nToys      number    of pseudoexperiments to generate
    - data       variable list of variables to be used for dataset generation
    - toy_config configuration of <code>pdf.generate</code>
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
    
    assert toy_config and 'nEvents' in toy_config,\
           'Number of events per toy must be specified via "toy_config" %s' % toy_config
    
    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult
    import ostap.fitting.basic 

    gparams = gen_pdf.pdf.getParameters ( None )
    varset  = ROOT.RooArgSet() 
    
    if   isinstance ( data , ROOT.RooArgSet  ) : pass 
    elif isinstance ( data , ROOT.RooAbsData ) : varset = data.varset() 
    else :
        for v in data :
            if   isinstance ( v , ROOT.RooAbsArg ) :
                varset.add ( v )
            elif isinstance ( v , string_types   ) and v in gparams :
                varset.add ( gparams [ v ] )
            else :
                raise TypeError('Invalid variable %s/%s' % ( v , type ( v ) ) )

    ## parameters for generation
            
    fix_gen_init = {}
    for v in gparams : fix_gen_init [ v.name ] = float ( v ) 

    fix_gen_pars = {}
    if   isinstance ( gen_pars , ROOT.RooFitResult ) :
        ps = gen_pars.dct_params()
        for p in ps       : fix_gen_pars [ p      ] = float ( ps       [ p ] )
    elif isinstance ( gen_pars , dict ) :
        for p in gen_pars : fix_gen_pars [ p      ] = float ( gen_pars [ p ] )
    else :
        for p in gen_pars : fix_gen_pars [ p.name ] = float ( p )

    
    ## parameters for fitting 

    fix_fit_init = {}
    fparams = fit_pdf.pdf.getParameters ( None )
    for v in fparams : fix_fit_init [ v.name ] = float ( v ) 
    
    fix_fit_pars = {}
    if   isinstance ( fit_pars , ROOT.RooFitResult ) :
        ps = fit_pars.dct_params()
        for p in ps       : fix_fit_pars [ p      ] = float ( ps       [ p ] )
    elif isinstance ( gen_pars , dict ) :
        for p in fit_pars : fix_fit_pars [ p      ] = float ( fit_pars [ p ] )
    else :
        for p in fit_pars : fix_gfit_pars [ p.name ] = float ( p )
    
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
        dataset = gen_pdf.generate ( varset = varset , **toy_config )  
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
            
    return results, stats 


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
