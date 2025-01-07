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
    "make_toys3"       , ## run fitting toys with special action (separate models to generate and fit)
    'make_jackknife'   , ## run Jackknife analysis 
    'make_bootstrap'   , ## run Bootstrapanalysis 
    "vars_transform"   , ## helper fnuction to transform the variables
    "print_stats"      , ## print toys      statistics 
    "print_jackknife"  , ## print jackknife statistics 
    "print_bootstrap"  , ## print bootstrap statistics 
    )
# =============================================================================
from   builtins               import range
from   ostap.core.ostap_types import string_types, integer_types
from   ostap.core.core        import VE, SE, Ostap 
from   ostap.logger.pretty    import pretty_ve, pretty_float, fmt_pretty_float
from   ostap.logger.colorized import attention
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.toys' )
else                       : logger = getLogger( __name__             )
# =============================================================================
logger.debug ( 'Utilities to run toys, Jackknife, bootstrap, ... ')
# ==============================================================================
## Check the applicablity of Jackknife
#  - for extended fits, parameters of yield are likely wrong!
def check_jackknife ( pdf ) :
    """ Check the applicablity of Jackknife
    - for extended fits, parameters of yield are likely wrong!
    """
    if not isinstance ( pdf , ROOT.RooAbsPdf ) : pdf = pdf.pdf 
    return False if  ( pdf.mustBeExtended () or pdf.canBeExtended  () ) else True 
# ===============================================================================
## check the applicability of bootstrap
#  - for extedned fits, non-extened bootstrap may produce wrong results!
def check_bootstrap ( pdf , extended ) : 
    """ Check the applicability of bootstrap
     - for extedned fits, non-extendd bootstrap may produce wrong results!
    """
    if not isinstance ( pdf , ROOT.RooAbsPdf ) : pdf = pdf.pdf     
    ok = True 
    if   pdf.mustBeExtended () and not extended : ok = False
    elif pdf.canBeExtended  () and not extended : ok = False
    return ok

# ==============================================================================
## Technical transformation  to the dictionary :  { 'name' : float_value }
def vars_transform ( vars ) :
    """ Technical transformation to the dictionary :  `{ 'name' : float_value }`
    """
    
    from   ostap.core.ostap_types     import dictlike_types
    import ostap.fitting.roofitresult
    import ostap.fitting.variables 
    
    result = {}
    if isinstance   ( vars , ROOT.RooFitResult ) :
        rs = vars.dct_params ()
        for p in  rs  : result [ p ] = float ( rs   [ p ] )
    elif isinstance ( vars , dictlike_types ) :
        for p in vars : result [ p ] = float ( vars [ p ] ) 
    else :
        for p in vars :
            if not isinstance ( p , ROOT.RooAbsCategory ) :
                result [ p.GetName()  ] = float ( p )
            
    return result

# =============================================================================
## helper class to get serializeable accessor to pull-variable 
class PoolVar(object) :
    def __init__ ( self , name , value  ) :
        self.__name   = name 
        self.__value  = float ( value ) * 1.0 
    def __call__ ( self , r , *args ) :
        v = getattr ( r , self.__name ) * 1 
        return ( v.value () - self.__value ) / v.error()
# =============================================================================
def pull_var ( name , value ) : return PoolVar ( name , value )
# =============================================================================
## Prepare statistics from results of toys/jackknifes/boostrap studies
#  @code
#  results    = ...
#  statistics = make_stats ( results )
#  @endcode 
def make_stats ( results , fits = {} , covs = {} , accept = SE () ) :
    """ Prepare statistics from results of toys/jackknifes/boostrap studies
    >>> results    = ...
    >>> statistics = make_stats ( results )
    """
    from collections     import defaultdict 
    from ostap.core.core import VE, SE

    stats = defaultdict ( SE )
    
    for par in results :
        if par : 
            pars = results [ par ] 
            for v in pars : stats [ par ] += float ( v )
            
    if isinstance ( accept , VE ) and 2 <= accept.nEntries() and 0 <= accept.eff() <= 1  :
        stats ['- Accept '         ] = accept 
    
    for k in fits :
        stats ['- Status  %s' % k  ] = fits [ k ]
    for k in covs :
        stats ['- CovQual %s' % k  ] = covs [ k ]
        
    return stats

# =============================================================================
## Print statistics of pseudoexperiments
def print_stats ( stats , ntoys = '' , logger = logger ) :
    """ Print statistics of pseudoexperiments
    """
    ##              0          1      2       3        4        5          6           7 
    table = [ ( 'Parameter' , '#', 'mean', 'x[..]' , 'rms' , 'x[..]', 'min / max' , 'x[..]' ) ] 

    def make_row ( c ) :    
        
        n      = '%d' % c.nEntries() 
        mean   = c.mean ()
        rms    = c.rms  ()
        minmax = c.min() , c.max() 
        
        mean , expo1 = pretty_ve    ( mean , precision = 6 , width = 8 , parentheses = False )
        rms  , expo2 = pretty_float ( rms  , precision = 4 , width = 6 )        
        fmt  , expo3 = fmt_pretty_float ( max ( abs ( c.min() ) , abs ( c.max() ) ) , precision = 3 , width = 5 )
        fmt  = '%s / %s ' % ( fmt , fmt )

        if expo3 : scale = 10**expo3
        else     : scale = 1 
        minmax = fmt % ( c.min() / scale , c.max() / scale ) 
        
        if expo1 : expo1 = '10^%+d' % expo1
        else     : expo1 = ''
        
        if expo2 : expo2 = '10^%+d' % expo2
        else     : expo2 = ''
        
        if expo3 : expo3 = '10^%+d' % expo3
        else     : expo3 = '' 

        return p , n , mean, expo1, rms, expo2, minmax, expo3 

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
    table = Table.remove_empty_columns ( table ) 
        
    table = Table.table ( table                                    ,
                          title     = "Results of %s toys" % ntoys ,
                          alignment = 'lccccccc'                   ,
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
    """ Jackknife estimator from jackknife statistic
    >>> statistics = ....
    >>> jacknife                      = jackknife_estimator ( statistics         )
    >>> jacknife , theta_corr , bias  = jackknife_esiimator ( statistics , value )
    """
    
    assert isinstance ( theta , VE ) or theta is None  ,\
           "jackknife_statistics: invalid type of 'value': %s" % type ( value ) 
    
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
def print_jackknife  ( fitresult          , ## restl of the fit of th e total datasample 
                       stats              ,
                       morevars  = {}     ,                       
                       logger    = logger ,
                       title     = ''     ) :
    """ Print Jackknife statistics
    """

    header = ( 'Parameter'  ,
               'theta'      , 'x[...]' ,
               'theta_(.)'  , 'x[...]' ,
               ## 'bias'    , 'x[...]' , 'bias [%]'   ,
               's/s<j>' )
    
    table  = [ header ]

    N = 0 
    for name in sorted ( stats ) :
        
        if   name in fitresult :
            p     = fitresult [ name ]
            theta = p * 1.0
            if not isinstance ( theta , VE ) or theta.cov2() <= 0 :
                logger.warning ("print_jackknife: parameter '%s' is invalid in 'fitresult', skip %s" % ( name , theta ) )
                continue
        elif name in morevars :
            theta = morevars [ name ]
            if not isinstance ( theta , VE ) or theta.cov2() <= 0 :
                logger.warning ("print_jackknife: parameter '%s' is invalid in 'morevars',  skip %s" % ( name , theta ) ) 
                continue
        else :
            continue 

        statistics  = stats [ name ]

        N = max ( N , statistics.nEntries() )
        
        
        ## jackknife estimates 
        jackknife , theta_jack = jackknife_statistics ( statistics , theta )

        bias =  ( N - 1 ) * ( jackknife - theta ).value()
        
        ## bias  = theta_jack.value () - theta     .value ()
        
        scale = theta     .error () / theta_jack.error () 
        
        th1 , expo1 = pretty_ve ( theta      , precision = 4 , width = 6 , parentheses = False )
        th2 , expo2 = pretty_ve ( jackknife  , precision = 4 , width = 6 , parentheses = False )
        
        ## th3 , expo3 = pretty_ve ( theta_jack , precision = 4 , width = 6 , parentheses = False )
        th3 , expo3 = pretty_float ( bias       , precision = 4 , width = 6 )

        if expo1 : expo1 = '10^%+d' % expo1
        else     : expo1 = ''
        
        if expo2 : expo2 = '10^%+d' % expo2
        else     : expo2 = ''

        if expo3 : expo3 = '10^%+d' % expo3
        else     : expo3 = ''

        bias = bias / theta.error() * 100 
        
        fbias = '%+.1f' % bias 
        errs  = '%+.2f' % scale

        if 50  < abs ( bias       ) : fbias = attention ( fbias )
        if 0.5 < abs ( scale  - 1 ) : errs  = attention ( errs  )
        
        ## row = name , th1, expo1 , th2 , expo2 , th3 , expo3 , fbias , errs
        row = name , th1, expo1 , th2 , expo2 , errs

        table.append ( row )

        
    for name in sorted ( stats ) :

        if name in fitresult : continue
        if name in morevars  : continue

        statistics = stats [ name ]
        jackknife  = jackknife_statistics ( statistics ) 

        th2, expo2 =  pretty_ve ( jackknife , precision = 4 , width = 6 , parentheses = False )

        if expo2 : expo2 = '10^%+d' % expo2
        else     : expo2 = ''
        
        ## row = name , '' , '' , th2 , expo2 , '' , '' , '' , '' 
        row = name , '' , '' , th2 , expo2 , '' 
        table.append ( row )
        
    import ostap.logger.table as Table
    table = Table.remove_empty_columns ( table )
    
    title = title if title else "Jackknife results (N=%d)" % N  
    table = Table.table ( table                           ,
                          title     = title               ,
                          alignment = 'lcccccc'           ,
                          prefix    = "# "                )
    logger.info ( '%s:\n%s' % ( title , table ) ) 
    

# =============================================================================
## print Bootstrap statistics
def print_bootstrap  ( fitresult          ,
                       stats              ,
                       morevars  = {}     ,
                       logger    = logger ,
                       title     = ''     ) :
    """ Print Bootstrap statistics
    """
    
    header = ( 'Parameter' , 'theta' , 'x[...]' ,  'theta_boot' , 'x[...]', 'bias[%]' , 's/s<b>' ) 
    table  = []

    n = 0

    for name in sorted ( stats ) :
        
        if   name in fitresult :
            p     = fitresult [ name ]
            theta = p * 1.0
            if not isinstance ( theta , VE ) or theta.cov2() <= 0 :
                logger.warning ("print_bootstrap: parameter '%s' is invalid in 'fitresult, skip %s" % ( name , theta ) )
                continue
        elif name in morevars :
            theta = morevars [ name ]
            if not isinstance ( theta , VE ) or theta.cov2() <= 0 :
                logger.warning ("print_bootstrap: parameter '%s' is invalid in 'morevars',  skip %s" % ( name , theta ) ) 
                continue
        else :
            continue 

        statistics  = stats [ name ]

        n = max ( n , statistics.nEntries() ) 

        theta_boot = VE ( statistics.mean().value() , statistics.mu2() ) 
        
        bias  = theta_boot.value () - theta     .value ()        
        scale = theta     .error () / theta_boot.error ()
        
        th1 , expo1 = pretty_ve ( theta      , precision = 4 , width = 6 , parentheses = False )
        th2 , expo2 = pretty_ve ( theta_boot , precision = 4 , width = 6 , parentheses = False )

        if expo1 : expo1 = '10^%+d' % expo1
        else     : expo1 = ''
        
        if expo2 : expo2 = '10^%+d' % expo2
        else     : expo2 = ''

        bias  = bias / theta.error() * 100
        
        fbias = '%+.1f' % bias 
        errs  = '%+.2f' % scale

        if 50  < abs ( bias       ) : fbias = attention ( fbias )
        if 0.5 < abs ( scale  - 1 ) : errs  = attention ( errs  )
        
        row = name , th1, expo1 , th2 , expo2 , fbias , errs 

        table.append ( row )

    for name in sorted ( stats ) :
        
        if name in fitresult : continue
        if name in morevars  : continue
        
        statistics = stats [ name ]
        theta_boot = VE ( statistics.mean().value() , statistics.mu2() ) 

        th2 , expo2 = pretty_ve ( theta_boot , precision = 4 , width = 6 , parentheses = False )
        
        if expo2 : expo2 = '10^%+d' % expo2
        else     : expo2 = ''

        row = name , '','' , th2 , expo2 , '' , '' 
        table.append ( row )

    table = [ header ] + table

    import ostap.logger.table as Table
    table = Table.remove_empty_columns ( table ) 
    title = title if title else "Bootstrapping with #%d samples" % n 
    table = Table.table ( table                ,
                          title     = title    ,
                          alignment = 'lccccccc'  ,
                          prefix    = "# "     )
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# ==============================================================================
## Default function to generate the data
#  - simple call for <code>PDF.generate</code>
def generate_data ( pdf , varset , **config ) :
    """ Default function to generate the data
    - simple call for `PDF.generate`
    """
    return pdf.generate ( varset = varset , **config )

# ==============================================================================
## Default function to perform the actual fit
#  - simple call for <code>PDF.fitTo</code>
def make_fit    ( pdf , dataset , **config ) :
    """ Default function to  perform the actual fit
    - simple call for `PDF.fitTo`
    """
    result = pdf.fitTo ( dataset , **config )
    if isinstance ( result , tuple ) and 2 == len ( result ) :
        result = result [ 0 ] 
    return result

# ==============================================================================
## Accept fit?
#  Accept the fit result?
#   - valid fit result
#   - fit status is 0 (SUCCESS)
#   - covariance matrix quality  is either 3(full accurate matrix) or -1 (unknown/externally provided?) 
#  @param result  fit result
#  @param pdf     pdf
#  @param dataset pdf
def accept_fit  ( result , pdf = None , dataset = None ) :
    """ Accept the fit result?
    - valid fit result
    - fit status is 0 (SUCCESS)
    - covariance matrix quality  is either 3(full accurate matrix) or -1 (unknown/externally provided?) 
    """
    return result and ( 0 == result.status () ) and ( result.covQual () in ( -1 , 3 ) ) 

# ==============================================================================
## make <code>nToys</code> pseudoexperiments
#
#  Schematically:
#  @code
#  for toy in range ( nToys )  :
#  ...  dataset    = gen_fun ( pdf , ...     , **gen_config )
#  ...  fit_result = fit_fun ( pdf , dataset , **fit_config )
#  ...  if not accept_fun ( fit_result , pdf , dataset ) : continue
#  .... < collect statistics here >
#  return results, stats 
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
# @param pdf          PDF to be used for generation and fitting
# @param nToys        number    of pseudoexperiments to generate
# @param data         variable list of variables to be used for dataset generation
# @param gen_config   configuration of <code>pdf.generate</code>
# @param fit_config   configuration of <code>pdf.fitTo</code>
# @param init_pars    redefine these parameters for each pseudoexperiment
# @param more_vars    calculate more variables form fit-result
# @param add_results  add fit resutls to the output?
# @param get_fun      specific generate action (if needed) 
# @param fit_fun      specific fitting action (if needed) 
# @param accept_fun   specific accept action (if needed) 
# @param silent       silent toys?
# @param progress     show the progress?
# @param logger       use this logger 
# @param frequency    how often to dump the intermediate results ? 
# @return dictionary with fit results for the toys and the dictionary of statistics
#
#  - If <code>gen_fun</code>    is not specified <code>generate_data</code> is used 
#  - If <code>fit_fun</code>    is not specified <code>make_fit</code>      is used 
#  - If <code>accept_fun</code> is not specified <code>accept_fit</code>    is used   
def make_toys ( pdf                   ,
                nToys                 , 
                data                  , ## template for dataset/variables 
                gen_config            , ## parameters for <code>pdf.generate</code>   
                fit_config   = {}     , ## parameters for <code>pdf.fitTo</code>
                init_pars    = {}     ,
                more_vars    = {}     ,
                add_results  = False  , ## add fit-results to the output?
                gen_fun      = None   , ## generator function ( pdf , varset  , **config )
                fit_fun      = None   , ## fit       function ( pdf , dataset , **config )
                accept_fun   = None   , ## accept    function ( fit-result, pdf, dataset )
                silent       = True   ,                
                progress     = True   ,
                logger       = logger ,
                frequency    = 500    ) : ## 
    """ Make `nToys` pseudoexperiments

    -   Schematically:
    >>> for toy in range ( nToys )  :
    >>> ...  dataset = gen_fun ( pdf , ...     , **gen_config )
    >>> ...  result  = fit_fun ( pdf , dataset , **fit_config )
    >>> ...  if not accept_fun ( result , pdf , dataset ) : continue
    >>> .... < collect statistics here > 
    >>> return results, stats 
    
    For each pseudoexperiment:

    1. generate dataset using `pdf` with variables specified
    in `data` and configuration specified via `gen_config`
    for each generation the parameters of `pdf` are reset
    for their initial values and valeus from `init_pars`
    
    2. fit generated dataset  with `pdf` using configuration
    specified via  `fit_config`

    - `pdf`          : PDF to be used for generation and fitting
    - `nToys`        : number    of pseudoexperiments to generate
    - `data`         : variable list of variables to be used for dataset generation
    - `gen_config`   : configuration of <code>pdf.generate</code>
    - `fit_config`   : configuration of <code>pdf.fitTo</code>
    - `init_pars`    : redefine these parameters for each pseudoexperiment
    - `more_vars`    : dictionary of functions to define the additional results
    - `add_results`  : add fit-resutls to the output?
    - `gen_fun`      : generator function
    - `fit_fun`      : fitting   function
    - `accept_fun`   : accept    function
    - `silent`       : silent toys?
    - `progress`     : show progress bar? 
    - `logger`       : use this logger 
    - `frequency`    : how often to dump the intermediate results ? 
    
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
        if not silent : logger.info ( "make_toys: use default 'generate_data' function!")
        gen_fun = generate_data 
    assert gen_fun and callable ( gen_fun ) , 'Invalid generator function!'
    
    ## 2. fitting function? 
    if fit_fun is None :
        if not silent : logger.info ( "make_toys: use default 'make_fit' function!")
        fit_fun = make_fit 
    assert fit_fun and callable ( fit_fun ) , 'Invalid fit function!'

    ## 3. accept function? 
    if accept_fun is None :
        if not silent : logger.info ( "make_toys: use default accept_fit' function!")
        accept_fun = accept_fit
    assert accept_fun and callable ( accept_fun ) , 'Invalid accept function!'

    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult
    import ostap.fitting.pdfbasic 
    import ostap.histos.histos   

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

    pdf.load_params ( params = fix_pars , silent = True ) ## silent = silent )
    pdf.load_params ( params = fix_init , silent = True ) ## silent = silent )

    ## save all initial parameters (needed for the final statistics)
    params  = pdf.params      ()
    fix_all = vars_transform  ( params ) 
    
    fitcnf = {}
    fitcnf.update ( fit_config )
    if not 'silent' in fitcnf : fitcnf [ 'silent' ] = silent

    from collections import defaultdict 
    results = defaultdict(list) 

    from   ostap.core.core        import SE, VE

    fits   = defaultdict ( SE )  ## fit statuses 
    covs   = defaultdict ( SE )  ## covariance matrix quality
    accept = SE ()
    
    ## run pseudoexperiments
    from ostap.utils.progress_bar import progress_bar 
    for i in progress_bar ( range ( nToys ) , silent = not progress , description = 'Toys:' ) :
                
        ## 1. reset PDF parameters 
        pdf.load_params ( params = fix_pars  , silent = True ) ## silent = silent )
        pdf.load_params ( params = init_pars , silent = True ) ## silent = silent )

        ## 2. generate dataset!  
        ## dataset = pdf.generate ( varset = varset , **gen_config )  
        dataset = gen_fun ( pdf , varset = varset , **gen_config )  
        if not silent : logger.info ( 'Generated dataset #%d\n%s' % ( i , dataset ) )

        ## histogram? 
        histo = dataset if isinstance ( dataset , ROOT.TH1 ) else None 
        
        ## 3. fit it!
        fit_result = fit_fun ( pdf , dataset , **fitcnf ) 
        
        ## fit status
        st = fit_result.status()
        if 0 != st :  fits [ st ] += 1

        ## covariance matrix quality
        cq = fit_result.covQual() 
        if not cq in ( -1 , 3 ) : covs [ cq ] += 1
              
        ## ok ?
        ok      = accept_fun ( fit_result , pdf , dataset )
        accept += 1 if ok else 0
        
        if ok :
            
            ## 4.1 save results 
            rpf = fit_result.params ( float_only = True ) 
            for p in rpf : 
                results [ p ].append ( rpf [ p ][0] ) 
                
            ## 4.2 save results 
            for v in more_vars :
                func  = more_vars [ v ] 
                results [ v ] .append ( func ( fit_result , pdf ) )

            if histo : 
                results [ '#'     ] .append ( histo.GetEntries   () )
                results [ '#sumw' ] .append ( histo.the_integral () )
            else :
                results [ '#'     ] .append ( len ( dataset ) )
                results [ '#sumw' ] .append ( dataset.sumVar ( '1' ) )
                
            ## 4.3 save results 
            if  add_results : results [ '' ].append ( fit_result )

        ## if not add_results and isinstance ( r , ROOT.RooFitResult ) :
        ##    r = Ostap.MoreRooFit.delete_result ( r )                                    
        if isinstance ( dataset , ROOT.RooAbsData ) :
            dataset = Ostap.MoreRooFit.delete_data ( dataset )
            
        del dataset
        del fit_result
        
        if progress or not silent :
            if 0 < frequency and 1 <= i and 0 ==  i % frequency : 
                stats = make_stats ( results , fits , covs , accept )
                print_stats ( stats , i + 1 , logger = logger )
           
    ## make a final statistics 
    stats = make_stats ( results , fits , covs , accept )
        
    if progress or not silent :
        print_stats ( stats , nToys , logger = logger )
    
    return results, stats 

# =============================================================================
## make <code>nToys</code> pseudoexperiments
#
#  Schematically:
#  @code
#  for toy in range ( nToys )  :
#  ...  dataset    = gen_fun ( gen_pdf , ...     , **gen_config )
#  ...  fit_result = fit_fun ( fit_pdf , dataset , **fit_config )
#  ...  if not accept_fun ( fit_result  , fit_pdf , dataset ) : continue
#  ...  result = action ( fit_result , fit_pdf , dataset )
#  ...  results [''].append ( result )
#  ...  <some statistics here> 
#  return results, stats 
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
# @param gen_pdf      PDF to be used for generation 
# @param fit_pdf      PDF to be used for fitting
# @param nToys        number    of pseudoexperiments to generate
# @param data         variable list of variables to be used for dataset generation
# @param gen_config   configuration of <code>pdf.generate</code>
# @param fit_config   configuration of <code>pdf.fitTo</code>
# @param gen_pars     redefine these parameters for each pseudoexperiment
# @param fit_pars     redefine these parameters for each pseudoexperiment
# @param more_vars    calculate more variables form fit-result
# @param add_results  add fit-resutls to the output?
# @param gen_fun      specific generate  action (if needed) 
# @param fit_fun      specific fitting action (if needed) 
# @param accept_fun   specific accept action (if needed) 
# @param silent       silent toys?
# @param progress     show progress bar?
# @param logger       logger 
# @param frequency    how often to dump the intermediate results ? 
# @return dictionary with fit results for the toys and the dictionary of statistics
#
#  - If <code>gen_fun</code>    is not specified <code>generate_data</code> is used 
#  - If <code>fit_fun</code>    is not specified <code>make_fit</code>      is used 
#  - If <code>accept_fun</code> is not specified <code>accept_fit</code>    is used   
def make_toys2 ( gen_pdf               , ## pdf to generate toys 
                 fit_pdf               , ## pdf to fit  
                 nToys                 , ## number of pseudoexperiments 
                 data                  , ## template for dataset/variables 
                 gen_config            , ## parameters for <code>pdf.generate</code>   
                 fit_config   = {}     , ## parameters for <code>pdf.fitTo</code>
                 gen_pars     = {}     , ## gen-parameters to reset/use 
                 fit_pars     = {}     , ## fit-parameters to reset/use
                 more_vars    = {}     , ## additional  results to be calculated
                 add_results  = False  , ## add fit-resutls to the output?
                 gen_fun      = None   , ## generator function ( pdf , varset  , **gen_config ) 
                 fit_fun      = None   , ## fit       function ( pdf , dataset , **fit_config ) 
                 accept_fun   = None   , ## accept    function ( fit-result, pdf, dataset     )
                 silent       = True   ,
                 progress     = True   ,
                 logger       = logger ,
                 frequency    = 500    ) :
    """ Make `nToys` pseudoexperiments
    
    -   Schematically:
    >>> for toy in range ( nToys )  :
    >>> ...  dataset    = gen_fun ( gen_pdf , ...     , **gen_config )
    >>> ...  fit_result = fit_fun ( fit_pdf , dataset , **fit_config )
    >>> ...  if not accept_fun ( fit_result  , fit_pdf , dataset ) : continue
    >>> ... < collect statistics here > 
    >>> return results, stats 
    
    For each experiment:

    1. generate dataset using `pdf` with variables specified
    in `data` and configuration specified via `gen_config`
    for each generation the parameters of `pdf` are reset
    for their initial values and valeus from `init_pars`
    
    2. fit generated dataset  with `pdf` using configuration
    specified via  `fit_config`
    
    - `gen_pdf`      : PDF to be used for generation 
    - `fit_pdf`      : PDF to be used for fitting 
    - `nToys`        : number    of pseudoexperiments to generate
    - `data`         : variable list of variables to be used for dataset generation
    - `gen_config`   : configuration of <code>pdf.generate</code>
    - `fit_config`   : configuration of <code>pdf.fitTo</code>
    - `gen_pars`     : redefine these parameters for generation of each pseudoexperiment
    - `fit_pars`     : redefine these parameters for fit of each pseudoexperiment
    - `more_vars`    : dictionary of functions to define the additional results
    - `add_results`  : add fit-results to the output?
    - `gen_fun`      : generator function ( pdf , varset  , **gen_config )
    - `fit_fun`      : fitting   function ( pdf , dataset , **fit_config ) 
    - `accept_fun`   : accept    function ( fit-result, pdf, dataset     )
    - `silent`       : silent toys?
    - `progress`     : show progress bar?
    - `logger`       : use this logger 
    - `frequency`    : how often to dump the intermediate results ? 
    
    It returns a dictionary with fit results for the toys and a dictionary of statistics
    >>> pdf = ...
    ... results, stats = make_toys2 ( 
    ...                 gen_pdf                                     , ## PDF  to generate 
    ...                 fit_pdf                                     , ## PDF  to fit 
    ...                 nToys       = 1000                          , ## number of toys 
    ...                 data       = [ 'mass' ]                     , ## varibales in dataset 
    ...                 gen_config = { 'nEvents' : 5000 }           , ## configuration of `pdf.generate`
    ...                 fit_config = { 'ncpus'   : 2    }           , ## configuration of `pdf.fitTo`
    ...                 gemn_pars  = { 'mean' : 0.0 , 'sigma' : 1.0 } ## parameters to use for generation 
    ...                )
    """

    from ostap.core.ostap_types import string_types, integer_types  
    
    assert isinstance ( nToys , integer_types ) and 0 < nToys,\
           'Invalid "nToys" argument %s/%s' % ( nToys , type ( nToys ) )
    
    assert gen_config and 'nEvents' in gen_config,\
           'Number of events per toy must be specified via "gen_config" %s' % gen_config
    
    ## 1. generator function? 
    if gen_fun is None :
        if not silent :  logger.info ( "make_toys2: use default 'generate_data' function!")
        gen_fun = generate_data 
    assert gen_fun and callable ( gen_fun ) , 'Invalid generator function!'
    
    ## 2. fitting function? 
    if fit_fun is None :
        if not silent :  logger.info ( "make_toys2: use default 'make_fit' function!")
        fit_fun = make_fit 
    assert fit_fun and callable ( fit_fun ) , 'Invalid fit function!'

    ## 3. accept function? 
    if accept_fun is None :
        if not silent : logger.info ( "make_toys2: use default 'accept_fit' function!")
        accept_fun = accept_fit
    assert accept_fun and callable ( accept_fun ) , 'Invalid accept function!'

    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult
    import ostap.fitting.pdfbasic 
    import ostap.histos.histos   

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
    
    fits   = defaultdict ( SE )  ## fit statuses 
    covs   = defaultdict ( SE )  ## covarinace matrix quality
    accept = SE()

    ## run pseudoexperiments
    from ostap.utils.progress_bar import progress_bar 
    for i in progress_bar ( range ( nToys ) , silent = not progress , description = 'Toys:' ) :

        ## 1. reset PDF parameters 
        gen_pdf.load_params ( params = fix_gen_init , silent = True ) ## silent = silent )
        gen_pdf.load_params ( params = fix_gen_pars , silent = True ) ## silent = silent )
        
        ## 2. generate dataset!
        dataset =  gen_fun ( gen_pdf , varset = varset , **gen_config ) 
        if not silent : logger.info ( 'Generated dataset #%d\n%s' % ( i , dataset ) )

        ## histogram? 
        histo = dataset if isinstance ( dataset , ROOT.TH1 ) else None 
        
        ## 3. reset parameters of fit_pdf
        fit_pdf.load_params ( params = fix_fit_init , silent = True ) ## silent = silent )
        fit_pdf.load_params ( params = fix_fit_pars , silent = True ) ## silent = silent )

        ## 4. fit it!  
        fit_result = fit_fun ( fit_pdf , dataset , **fitcnf )

        ## fit status 
        st = fit_result.status()
        if 0 != st :  fits [ st ] += 1
        
        ## covariance matrix quality
        cq = fit_result.covQual() 
        if not cq in ( -1 , 3 ) : covs [ cq ] += 1
        
        ## ok ?
        ok  = accept_fun ( fit_result , fit_pdf , dataset )
        accept += 1 if ok else 0
        
        if ok : 
            
            ## 5.1 save results 
            rpf = fit_result.params ( float_only = True ) 
            for j in rpf : 
                results [ j ].append ( rpf [ j ] [ 0 ] ) 
                    
            ## 5.2 save results 
            for v in more_vars :
                func  = more_vars[v] 
                results [ v ] .append ( func ( fit_result , fit_pdf ) )
                
            if histo :
                results [ '#'     ] .append ( histo.GetEntries   () )
                results [ '#sumw' ] .append ( histo.the_integral () )
            else :                
                results [ '#'     ] .append ( len ( dataset ) )
                results [ '#sumw' ] .append ( dataset.sumVar ( '1' ) ) 
                
            ## 5.3 save results 
            if add_results  : results [ '' ].append ( fit_result )

        ## if not add_results and isinstance ( r , ROOT.RooFitResult ) :
        ##     r = Ostap.MoreRooFit.delete_result ( r )            
        if isinstance ( dataset , ROOT.RooAbsData ) :
            dataset = Ostap.MoreRooFit.delete_data ( dataset ) 
            
        del dataset
        del fit_result
        
        if progress or not silent :
            if 0 < frequency and 1 <= i and 0 == i % frequency : 
                stats = make_stats ( results , fits , covs , accept )
                print_stats ( stats , i + 1 , logger = logger )

    ## make a final statistics 
    stats = make_stats ( results , fits , covs , accept )
                    
    if progress or not silent :
        print_stats ( stats , nToys , logger = logger  )
    
    return results, stats 


# =============================================================================
## make <code>nToys</code> pseudoexperiments
#
#  Schematically:
#  @code
#  for toy in range ( nToys )  :
#  ...  dataset    = gen_fun ( gen_pdf , ...     , **gen_config )
#  ...  fit_result = fit_fun ( fit_pdf , dataset , **fit_config )
#  ...  if not accept_fun ( fit_result  , fit_pdf , dataset ) : continue
#  ...  result = action ( fit_results , fit_pdf , dataset )
#  ...  results [''].append ( result ) 
#  ...  <collect statistics>
#  return results, stats 
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
# results , stats = make_toys3(
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
# @param gen_pdf      PDF to be used for generation 
# @param fit_pdf      PDF to be used for fitting
# @param nToys        number    of pseudoexperiments to generate
# @param data         variable list of variables to be used for dataset generation
# @param action       the action: result = action ( fit_result , fit_pdf , dataset ) 
# @param gen_config   configuration of <code>pdf.generate</code>
# @param fit_config   configuration of <code>pdf.fitTo</code>
# @param gen_pars     redefine these parameters for each pseudoexperiment
# @param gen_fun      specific generate  action (if needed) 
# @param fit_fun      specific fitting action (if needed) 
# @param accept_fun   specific accept action (if needed) 
# @param silent       silent toys?
# @param progress     show progress bar?
# @param logger       logger 
# @param frequency    how often to dump the intermediate results ? 
# @return dictionary with fit results for the toys and the dictionary of statistics
#
#  - If <code>gen_fun</code>    is not specified <code>generate_data</code> is used 
#  - If <code>fit_fun</code>    is not specified <code>make_fit</code>      is used 
#  - If <code>accept_fun</code> is not specified <code>accept_fit</code>    is used   
def make_toys3 ( gen_pdf               , ## pdf to generate toys 
                 fit_pdf               , ## pdf to fit  
                 nToys                 , ## number of pseudoexperiments 
                 data                  , ## template for dataset/variables
                 action                , ## action:  result= action ( fit_result , fit_pdf , dataset ) 
                 gen_config            , ## parameters for <code>pdf.generate</code>   
                 fit_config   = {}     , ## parameters for <code>pdf.fitTo</code>
                 gen_pars     = {}     , ## gen-parameters to reset/use 
                 fit_pars     = {}     , ## fit-parameters to reset/use
                 gen_fun      = None   , ## generator function ( pdf , varset  , **gen_config ) 
                 fit_fun      = None   , ## fit       function ( pdf , dataset , **fit_config ) 
                 accept_fun   = None   , ## accept    function ( fit-result, pdf, dataset     )
                 silent       = True   ,
                 progress     = True   ,
                 logger       = logger ,
                 frequency    = 500    ) :
    """ Make `nToys` pseudoexperiments
    
    -   Schematically:
    >>> for toy in range ( nToys )  :
    >>> ...  dataset = gen_fun ( gen_pdf , ...     , **gen_config )
    >>> ...  result  = fit_fun ( fit_pdf , dataset , **fit_config )
    >>> ...  if not accept_fun ( result  , fit_pdf , dataset ) : continue
    >>> ...  results[''].append ( action ( result , fit_pdf , dataset ) )
    >>> ...  <some statistics here>
    >>> return results, stats  
    
    For each experiment:

    1. generate dataset using `pdf` with variables specified
    in `data` and configuration specified via `gen_config`
    for each generation the parameters of `pdf` are reset
    for their initial values and valeus from `init_pars`
    
    2. fit generated dataset  with `pdf` using configuration
    specified via  `fit_config`

    - `pdf`          : PDF to be used for generation and fitting
    - `nToys`        : number    of pseudoexperiments to generate
    - `data`         : variable list of variables to be used for dataset generation
    - `gen_config`   : configuration of <code>pdf.generate</code>
    - `fit_config`   : configuration of <code>pdf.fitTo</code>
    - `gen_pars`     : redefine these parameters for generation of each pseudoexperiment
    - `fit_pars`     : redefine these parameters for fit of each pseudoexperiment
    - `more_vars`    : dictionary of functions to define the additional results
    - `add_results`  : add fit-results to the output?
    - `gen_fun`      : generator function ( pdf , varset  , **gen_config )
    - `fit_fun`      : fitting   function ( pdf , dataset , **fit_config ) 
    - `accept_fun`   : accept    function ( fit-result, pdf, dataset     )
    - `silent`       : silent toys?
    - `progress`     : show progress bar?
    - `logger`       : use this logger 
    - `frequency`    : how often to dump the intermediate results ? 
    
    It returns a dictionary with fit results for the toys and a dictionary of statistics
    >>> pdf = ...
    ... results, stats = make_toys3 ( 
    ...                 gen_pdf                                     , ## PDF  to generate 
    ...                 fit_pdf                                     , ## PDF  to fit 
    ...                 nToys      = 1000                           , ## number of toys 
    ...                 data       = [ 'mass' ]                     , ## variables in dataset 
    ...                 action     = ...                            , ## the action 
    ...                 gen_config = { 'nEvents' : 5000 }           , ## configuration of `pdf.generate`
    ...                 fit_config = { 'ncpus'   : 2    }           , ## configuration of `pdf.fitTo`
    ...                 gemn_pars  = { 'mean' : 0.0 , 'sigma' : 1.0 } ## parameters to use for generation 
    ...                )
    """

    from ostap.core.ostap_types import string_types, integer_types  
    
    assert isinstance ( nToys , integer_types ) and 0 < nToys,\
           'Invalid "nToys" argument %s/%s' % ( nToys , type ( nToys ) )
    
    assert gen_config and 'nEvents' in gen_config,\
           'Number of events per toy must be specified via "gen_config" %s' % gen_config
    
    ## 1. generator function? 
    if gen_fun is None :
        if not silent :  logger.info ( "make_toys2: use default 'generate_data' function!")
        gen_fun = generate_data 
    assert gen_fun and callable ( gen_fun ) , 'Invalid generator function!'
    
    ## 2. fitting function? 
    if fit_fun is None :
        if not silent :  logger.info ( "make_toys2: use default 'make_fit' function!")
        fit_fun = make_fit 
    assert fit_fun and callable ( fit_fun ) , 'Invalid fit function!'

    ## 3. accept function? 
    if accept_fun is None :
        if not silent : logger.info ( "make_toys2: use default 'accept_fit' function!")
        accept_fun = accept_fit
    assert accept_fun and callable ( accept_fun ) , 'Invalid accept function!'

    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult
    import ostap.fitting.pdfbasic 
    import ostap.histos.histos   

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
    
    fits   = defaultdict ( SE )  ## fit statuses 
    covs   = defaultdict ( SE )  ## covarinace matrix quality
    accept = SE()

    ## run pseudoexperiments
    from ostap.utils.progress_bar import progress_bar 
    for i in progress_bar ( range ( nToys ) , silent = not progress , description = 'Toys:' ) :

        ## 1. reset PDF parameters 
        gen_pdf.load_params ( params = fix_gen_init , silent = True ) ## silent = silent )
        gen_pdf.load_params ( params = fix_gen_pars , silent = True ) ## silent = silent )
        
        ## 2. generate dataset!
        dataset =  gen_fun ( gen_pdf , varset = varset , **gen_config ) 
        if not silent : logger.info ( 'Generated dataset #%d\n%s' % ( i , dataset ) )

        ## histogram? 
        histo = dataset if isinstance ( dataset , ROOT.TH1 ) else None 
        
        ## 3. reset parameters of fit_pdf
        fit_pdf.load_params ( params = fix_fit_init , silent = True ) ## silent = silent )
        fit_pdf.load_params ( params = fix_fit_pars , silent = True ) ## silent = silent )

        ## 4. fit it!  
        fit_result = fit_fun ( fit_pdf , dataset , **fitcnf )

        ## fit status 
        st = fit_result.status()
        if 0 != st :  fits [ st ] += 1
        
        ## covariance matrix quality
        cq = fit_result.covQual() 
        if not cq in ( -1 , 3 ) : covs [ cq ] += 1
        
        ## ok ?
        ok  = accept_fun ( fit_result , fit_pdf , dataset )
        accept += 1 if ok else 0
        
        if ok : 

            ## action!
            result = action ( fit_result , fit_pdf , dataset )
            results [ '' ].append ( result )
                        
            if histo :
                results [ '#'     ] .append ( histo.GetEntries   () )
                results [ '#sumw' ] .append ( histo.the_integral () )
            else :                
                results [ '#'     ] .append ( len ( dataset ) )
                results [ '#sumw' ] .append ( dataset.sumVar ( '1' ) ) 
                                        
        if isinstance ( dataset , ROOT.RooAbsData ) :
            dataset = Ostap.MoreRooFit.delete_data ( dataset ) 
            
        del dataset
        del fit_result
        
        if progress or not silent :
            if 0 < frequency and 1 <= i and 0 == i % frequency : 
                stats = make_stats ( results , fits , covs , accept )
                print_stats ( stats , i + 1 , logger = logger )

    ## make a final statistics 
    stats = make_stats ( results , fits , covs , accept )
                    
    if progress or not silent :
        print_stats ( stats , nToys , logger = logger  )
    
    return results, stats 

# =============================================================================
## run Jackknife analysis, useful for evaluation of fit biases and uncertainty estimates
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
#
#  Derived parameters can be also retrived via <code>more_vars</code> argument:
#  @code
#  ratio     = lambda res,pdf : res.ratio('x','y')
#  more_vars = { 'Ratio' : ratio }
#  r,  s = make_jackknife ( .... , more_vars = more_vars , ... ) 
#  @endcode
#
#  @see https://en.wikipedia.org/wiki/Jackknife_resampling
#  @param pdf         fit model
#  @param data        original dataset
#  @param fit_config  configuration of <code>pdf.FitTo( data , ... )</code>
#  @param fit_pars    redefine these parameters before each fit 
#  @param more_vars   calculate more variables from the fit-results
#  @param add_results add fit results to the output?
#  @param fit_fun     fitting   function
#  @param accept_fun  accept    function
#  @param event_range event range to use for jackknife   
#  @param silent      silent processing 
#  @param progress    show progress bar?
#  @param logger      use this logger
#  @param frequency  how often to dump the intermediate results ? 
#  @return statistics of jackknife experiments 
def make_jackknife ( pdf                  ,
                     data                 ,
                     fit_config  = {}     , ## parameters for <code>pdf.fitTo</code>
                     fit_pars    = {}     , ## fit-parameters to reset/use
                     more_vars   = {}     , ## additional  results to be calculated
                     add_results = False  , ## add fit-results to the output?                  
                     fit_fun     = None   , ## - fit       function ( pdf , dataset , **fit_config ) 
                     accept_fun  = None   , ## - accept    function ( fit-result, pdf, dataset     )
                     event_range = ()     , ## event range for jackknife                      
                     silent      = True   ,
                     progress    = True   ,
                     logger      = logger ,
                     frequency   = 500    ) :
    """ Run Jackknife analysis, useful for evaluaton of fit biased and uncertainty estimates
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
    - `add_results` : add fit resulsy to the output?
    - `fit_fun`     : specific fitting acion (if needed) 
    - `accept_fun`  : specific accept action (if needed)
    - `event_range` : event range to use for jackknife   
    - `silent`      : silent processing?
    - `progress`    : show progress bar?
    - `logger`      : use this logger 
    - `frequency`   : how often to dump the intermediate results ? 
    """
    
    N = len ( data )
    assert 1 < N            , 'make_jackknife: invalid dataset size %s' % N

    if not event_range : event_range = 0 , N 
    assert 2 == len ( event_range ) , 'make_jackknife: invalid event range %s ' % str ( event_range )
    
    begin , end = event_range

    ## check begin/end range
    assert 0 <= begin and begin < end and begin < N , 'make_jackknife: invalid event range (%s,%s)/%d' % ( begin , end , N )
    ## adjust the end 
    end   = min ( end   , N )

    ## check if pdf is extended 
    ok1 = False if  pdf.pdf.mustBeExtended () or pdf.pdf.canBeExtended  () else True 
    if not ok1 and not silent : logger.warning ( "Jackknife: estimates for `yield` parameters might be wrong!")
        
    ## 1. fitting function? 
    if fit_fun is None :
        if not silent :  logger.info ( "make_jackknife: use default 'make_fit' function!")
        fit_fun = make_fit 
    assert fit_fun and callable ( fit_fun ) , 'Invalid fit function!'

    ## 2. accept function? 
    if accept_fun is None :
        if not silent : logger.info ( "make_jackknife: use default 'accept_fit' function!")
        accept_fun = accept_fit
    assert accept_fun and callable ( accept_fun ) , 'Invalid accept function!'
    
    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult
    import ostap.fitting.pdfbasic 
       
    ## parameters for fitting 

    fparams      = pdf.params ()
    fix_fit_init = vars_transform ( fparams  )     
    fix_fit_pars = vars_transform ( fit_pars )
    
    fitcnf = {}
    fitcnf.update ( fit_config )
    if not 'silent' in fitcnf : fitcnf [ 'silent' ] = silent
    
    from collections import defaultdict 
    results = defaultdict(list) 

    from     ostap.core.core        import SE    
    fits   = defaultdict ( SE )  ## fit statuses 
    covs   = defaultdict ( SE )  ## covarinace matrix quality
    accept = SE ()
    
    ## Fit the whole sample 
    pdf.load_params ( params = fix_fit_init , silent = True ) ## silent = silent )
    pdf.load_params ( params = fix_fit_pars , silent = True ) ## silent = silent )
    r_tot = fit_fun ( pdf , data , **fitcnf )

    NN = 0 
    from ostap.utils.progress_bar import progress_bar
    ## run jackknife  bootstrapping
    for i , ds in progress_bar ( enumerate ( data.jackknife ( begin , end ) ) ,
                                 max_value   = end - begin  ,
                                 description = 'Sample:'    , 
                                 silent      = not progress ) :

        ## 2. reset parameters of fit_pdf
        pdf.load_params ( params = fix_fit_init , silent = True ) ## silent = silent )
        pdf.load_params ( params = fix_fit_pars , silent = True ) ## silent = silent )
 
        ## 3. fit it!  
        r = fit_fun ( pdf , ds , **fitcnf ) 

        ## 4. fit status 
        st = r.status()
        if 0 != st :  fits [ st ] += 1

        ## 5. covariance matrix quality
        cq = r.covQual() 
        if not cq in ( -1 , 3 ) : covs [ cq ] += 1

        ## ok ?
        ok      =  accept_fun ( r , pdf , ds )
        accept += ( 1 if ok else 0 )             
        if ok : 

            ## 6.1 save results 
            rpf = r.params ( float_only = True ) 
            for j in rpf : 
                results [ j ].append ( rpf [ j ] [ 0 ] ) 

            ## 6.2. more variables to be calculated? 
            for v in more_vars :
                func  = more_vars[v] 
                results [ v ] .append ( func ( r , pdf ) )
                
            results [ '#'     ] .append ( len ( ds ) )
            results [ '#sumw' ] .append ( ds.sumVar ( '1' ) )
            
            ## 6.3 save results 
            if   add_results                          : results [ '' ].append ( r )

            NN += 1
            
        if progress or not silent :
            if 0 < frequency and 1 <= i and 0 == i % frequency : 
                stats = make_stats ( results , fits , covs , accept )
                print_stats ( stats , i + 1 , logger = logger )
                
        ## if not add_result and isinstance ( r , ROOT.RooFitResult ) : r.Delete()                                
        ## reset/remove/delete dataset 
        if isinstance ( ds , ROOT.RooAbsData ) : ds = Ostap.MoreRooFit.delete_data ( ds ) 
        
        del ds
        del r
        
    ## 8. make a final statistics 
    stats = make_stats ( results , fits , covs , accept )
        
    if progress or not silent :

        ## 9. fit total dataset (twice) 
        r_tot = fit_fun ( pdf , data , **fitcnf )
        r_tot = fit_fun ( pdf , data , **fitcnf )
        
        ## 10. the final table
        title = '' 
        if not ok1 :
            logger.warning ( "Jackknife: estimates for `yield` parameters are likely wrong!")
            title = "Jackknife results (N=%d).[Estimates for `yield` are likely wrong!]" % NN            
        print_jackknife ( r_tot   ,
                          stats   ,
                          morevars = dict ( ( k , more_vars [ k ]( r_tot , pdf ) ) for k in more_vars ) ,
                          logger   = logger ,
                          title    = title  ) 
        
    if not ok1 and not silent : logger.warning ( "Jackknife: estimates for `yield` parameters are likely wrong!")
    return results , stats 

# =============================================================================
## Run Bootstrap analysis, useful for evaluation of fit biases and uncertainty estimates
# 
#  In total <code>size</code> datasets are sampled (with replacement) from the original dataset
#  <code>data</code> and each sampled dataset is fit
#  @code
#  dataset = ...
#  model   = ...
#  r , f = model.fitTo ( dataset , .... )                         ## fit the whole dataset   
#  results, stats = make_bootstrap ( model , data , size = 1000 ) ## run Bootstrap 
#  print_bootstrap ( r , stats )                    ## print summary table 
#  @endcode
#  @see print_bootstrap
#
#  Derived parameters can be also retrived via <code>more_vars</code> argument:
#  @code
#  ratio     = lambda res,pdf : res.ratio('x','y')
#  more_vars = { 'Ratio' : ratio }
#  r,  s = make_bootstrap ( .... , more_vars = more_vars , ... ) 
#  @endcode
#
#  @param pdf   fit model
#  @param data  original dataset
#  @param size  number of datasets to sample
#  @param fit_config configuration of <code>pdf.FitTo( data , ... )</code>
#  @param fit_pars   redefine these parameters before each fit 
#  @param more_vars  calculate more variables from the fit-results 
#  @param add_results add fit results to the output?
#  @param extended   use extended bootstrap? 
#  @param fit_fun    specific fitting action (if needed) 
#  @param accept_fun specific accept action (if needed) 
#  @param silent     silent processing 
#  @param progress   show progress bar?
#  @param logger     use this logger
#  @param frequency  how often dump the intermediate results? 
#  @return statistics of boostrap experiments 
def make_bootstrap (
        pdf                  ,
        data                 ,
        size        = 100    ,   ## numbere of samples 
        fit_config  = {}     ,   ## parameters for <code>pdf.fitTo</code>
        fit_pars    = {}     ,   ## fit-parameters to reset/use
        more_vars   = {}     ,   ## additional  results to be calculated
        add_results = False  ,   ## add fit results to the output?
        extended    = True   ,   ## use extended/non-extended bootstrtap 
        fit_fun     = None   ,   ## fit       function ( pdf , dataset , **fit_config )                     
        accept_fun  = None   ,   ## accept    function ( fit-result, pdf, dataset     )
        silent      = True   ,   ## silent processing?
        progress    = True   ,   ## show progress bar? 
        logger      = logger ,   ## use this logger 
        frequency   = 500    ) :
    
    """ Run Bootstrap analysis, useful for evaluaton of fit biased and uncertainty estimates 
    In total `size` datasets are sampled (with replacement) from the original dataste
    `data` and each sampled dataset is fit
    >>> dataset = ...
    >>> model   = ...
    >>> r , f = model.fitTo ( dataset , .... )                         ## fit the whole dataset   
    >>> results, stats = make_bootstrap ( model , data , size = 1000 ) ## run Bootstrap 
    >>> print_bootstrap ( r , stats )                    ## print summary table 

    - `pdf`         : fit model
    - `data`        : original dataset
    - `size`        : number of datasets to sample
    - `fit_config`  : configuration of `pdf.FitTo( data , ... )`
    - `fit_pars`    : redefine these parameters before each fit
    - `more_vars`   : calculate more variables from the fit-results
    - `add_results` : add fit-results to the output?
    - `extended`    : use extended bootstrap? 
    - `fit_fun`     : specific fitting acion (if needed) 
    - `accept_fun`  : specific accept action (if needed) 
    - `silent`      : silent processing?
    - `progress`    : show progress bar?
    - `logger`      : use this logger 
    - `frequency`   : how often dump the intermediate results? 
    """
    
    N = len ( data )
    assert 1 < N            , 'make_bootstrap: invalid dataset size %s' % N

    from ostap.core.ostap_types import integer_types  
    assert isinstance ( size  , integer_types ) and 1 <= size ,\
           'Invalid "size"  argument %s/%s' % ( size  , type ( size ) )
    
    ## check if pdf and `extended` flag are in agreemrnt, and print warning message otherwise
    ok1 = check_bootstrap ( pdf , extended )
    if not ok1 and not silent : logger.warning ( "Bootstrap: estimates for `yield` parameters might be wrong!")
        
    ## 1. fitting function? 
    if fit_fun is None :
        if not silent :  logger.info ( "make_bootstrap: use default 'make_fit' function!")
        fit_fun = make_fit 
    assert fit_fun and callable ( fit_fun ) , 'Invalid fit function!'

    ## 2. accept function? 
    if accept_fun is None :
        if not silent : logger.info ( "make_bootstrap: use default 'accept_fit' function!")
        accept_fun = accept_fit
    assert accept_fun and callable ( accept_fun ) , 'Invalid accept function!'

    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult
    import ostap.fitting.pdfbasic 
       
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
    fits   = defaultdict ( SE )  ## fit statuses 
    covs   = defaultdict ( SE )  ## covarinace matrix quality
    accept = SE()
    
    ## fit original dataset 
    pdf.load_params ( params = fix_fit_init , silent = True ) ## silent = silent )
    pdf.load_params ( params = fix_fit_pars , silent = True ) ## silent = silent )
    r_tot = fit_fun ( pdf , data , **fitcnf )

    
    NN = 0
    from ostap.utils.progress_bar import progress_bar
    ## run bootstrapping
    for i , ds in progress_bar ( enumerate ( data.bootstrap ( size , extended = extended ) ) ,
                                 max_value   = size         ,
                                 description = 'Sample:'    , 
                                 silent      = not progress ) :
        
        ## 2. reset parameters of fit_pdf
        pdf.load_params ( params = fix_fit_init , silent = True ) ## silent = silent )
        pdf.load_params ( params = fix_fit_pars , silent = True ) ## silent = silent )
 
        ## 3. fit it!  
        r = fit_fun ( pdf , ds , **fitcnf ) 

        ## 4. fit status 
        st = r.status()
        if 0 != st :  fits [ st ] += 1

        ## 5. covariance matrix quality
        cq = r.covQual() 
        if not cq in ( -1 , 3 ) : covs [ cq ] += 1

        ## ok ?
        ok      = accept_fun ( r , pdf , ds ) 
        accept += ( 1 if ok else 0 )
        if ok : 

            ## 6.1 save results 
            rpf = r.params ( float_only = True ) 
            for j in rpf : 
                results [ j ].append ( rpf [ j ] [ 0 ] ) 

            ## 6.2 more variables to be calculated? 
            for v in more_vars :
                func  = more_vars[v] 
                results [ v ] .append ( func ( r , pdf ) )
                
            results [ '#'     ] .append ( len ( ds ) )
            results [ '#sumw' ] .append ( ds.sumVar ( '1' ) )
            
            ## 6.3 save results 
            if   add_results                          : results [ '' ].append ( r )

            NN += 1

        if progress or not silent :
            if 0 < frequency and 1 <= i and 0 ==  i % frequency : 
                stats = make_stats ( results , fits , covs , accept )                
                title = "Bootstrapping with #%d/%d samples" % ( NN , i ) 
                if extended : title = '(Extended) %s' % title
                if not ok1  : title += '[Estimates for `yield` are likely wrong!]'
                print_bootstrap ( r_tot ,
                                  stats ,
                                  morevars = dict ( ( k , more_vars [ k ] ( r_tot , pdf ) ) for k in more_vars ),
                                  logger   = logger , 
                                  title    = title  )
                if not ok1 and not silent : logger.warning ( "Bootstrap: estimates for `yield` parameters are likely wrong!")

        ## if not add_results and isinstance ( r , ROOT.RooFitResult ) : r.Delete()            
        if isinstance ( ds , ROOT.RooAbsData ) : ds = Ostap.MoreRooFit.delete_data ( ds ) 
        
        del ds
        del r
        
    ## 8. make a final statistics 
    stats = make_stats ( results , fits , covs , accept )    
        
    if progress or not silent :

        ## 9. fit total dataset (twice) 
        r_tot = fit_fun ( pdf , data , **fitcnf )
        r_tot = fit_fun ( pdf , data , **fitcnf )
        
        ## 10. the final table  
        title = "Bootstrapping with #%d samples" % NN 
        if extended : title = '(Extended) %s'    % title
        if not ok1 :
            title += '[Estimates for `yield` are likely wrong!]'
            logger.warning ( "Bootstrap: estimates for `yield` parameters are likely wrong!")
        print_bootstrap ( r_tot ,
                          stats ,
                          morevars = dict ( ( k , more_vars [ k ]( r_tot , pdf ) ) for k in more_vars ),
                          logger   = logger ,
                          title    = title  ) 
            
        if not ok1 and not silent : logger.warning ( "Bootstrap: estimates for `yield` parameters are likely wrong!")
        
    return results , stats 


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
