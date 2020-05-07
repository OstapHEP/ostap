#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file ostap/parallel/parallel_toys.py
#  Make fitting toys  in parallel
#  @see ostap.fitting.toys 
#  @date   2020-01-18
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
# =============================================================================
""" Make fitting toys in parallel
- see ostap.fitting.toys 
"""
# =============================================================================
__author__  = 'Vanya BELYAEV  Ivan.Belyaev@itep.ru'
__date__    = "2020-01-18"
__version__ = '$Revision$'
__all__     = (
    'parallel_toys'  , ## run parallel toys (single   PDF  to generate and fit)
    'parallel_toys2' , ## run parallel toys (separate PDFs to generate and fit) 
    )
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.toys' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
import ROOT
from   ostap.parallel.parallel import Task, WorkManager
from   ostap.core.ostap_types  import string_types, integer_types  
# =============================================================================
## The simple task object for parallel fitting toys 
#  - single PDF to generate and fit
#  @see ostap.fitting.toys 
#  @see ostap.fitting.make_toys 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-01-18 
class  ToysTask(Task) :
    """The simple task object for parallel fitting toys
    - single PDF to generate and fit
    - see ostap.fitting.toys.make_toys
    """
    ## 
    def __init__ ( self               ,
                   pdf                ,
                   data               ,
                   gen_config         ,
                   fit_config = {}    ,
                   init_pars  = {}    ,
                   more_vars  = {}    ,                  
                   silent     = True  ,
                   progress   = False ) :
        
        self.pdf        = pdf
                
        self.data       = data  
        self.gen_config = gen_config 
        self.fit_config = fit_config
        self.init_pars  = init_pars 
        self.more_vars  = more_vars  

        self.silent     = silent
        self.progress   = progress 
        
        self.__the_output   = () 

    @property
    def the_output ( self ) :
        return self.__the_output
    @the_output.setter 
    def the_output ( self , value ) :
        self.__the_output = value 
    
    def initialize_local   ( self ) : self.__the_output = ()

    ## initialize the remote task, treta the random numbers  
    def initialize_remote  ( self , jobid = -1 ) :
        """Initialize the remote task, treta the random numbers  
        """

        import random, ROOT
        from ostap.parallel.utils import random_random
        random_random ( jobid )
        
        return self.initialize_local() 
    
        
    ## get the results 
    def results ( self ) :
        return self.__the_output
    
    ## the actual processing 
    def process ( self , jobid , nToys ) :

        import ROOT
        from ostap.logger.logger import logWarning
        with logWarning() :
            import ostap.core.pyrouts            
            import ostap.fitting.roofit            
            import ostap.fitting.dataset            
            import ostap.fitting.roofitresult            
            import ostap.fitting.variables
            
        from   ostap.core.ostap_types import integer_types 
        assert isinstance ( nToys , integer_types ) and 0 < nToys,\
               'Jobid %s: Invalid "nToys" argument %s/%s' % ( jobid , nToys , type ( nToys ) )
        
        import ostap.fitting.toys as Toys 
        results , stats = Toys.make_toys ( pdf        = self.pdf        ,
                                           nToys      = nToys           ,
                                           data       = self.data       ,
                                           gen_config = self.gen_config , 
                                           fit_config = self.fit_config , 
                                           init_pars  = self.init_pars  ,
                                           more_vars  = self.more_vars  ,
                                           silent     = self.silent     ,
                                           progress   = self.progress   )
        
        
        self.the_output = results , stats
        
        return self.results() 

    ## merge results/datasets 
    def merge_results ( self, result ) :
        
        if result :
            results , stat = result
            if not self.the_output :
                self.__the_output = results , stat 
            else :
                results  , stat    = result
                results_ , stat_   = self.__the_output 

                rset = set()
                for p in results  : rset.add ( p )
                for p in results_ : rset.add ( p )                
                for p in rset : results_[ p ] += results [ p ]

                sset = set()
                for p in stat  : sset.add ( p )
                for p in stat_ : sset.add ( p )
                for p in sset  : stat_ [ p  ] += stat [ p ] 
                
                self.__the_output = results_ , stat_ 
        else :
            logger.error ( "No valid results for merging" )
            

    
# =============================================================================
## The simple task object for parallel fitting toys
#  - separate PDFs to generarte and fit 
#  @see ostap.fitting.toys 
#  @see ostap.fitting.toys.make_toys2 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-01-18 
class  ToysTask2(ToysTask) :
    """The simple task object for parallel fitting toys
    - separate PDFs to generarte and fit 
    - see ostap.fitting.toys.make_toys2 
    """
    ## 
    def __init__ ( self               ,
                   gen_pdf            ,
                   fit_pdf            ,
                   data               ,
                   gen_config         ,
                   fit_config = {}    ,
                   gen_pars   = {}    ,
                   fit_pars   = {}    ,
                   more_vars  = {}    ,                  
                   silent     = True  ,
                   progress   = False ) :

        ToysTask.__init__ ( self                    ,
                            pdf        = gen_pdf    ,
                            data       = data       ,
                            gen_config = gen_config ,
                            fit_config = fit_config ,
                            init_pars  = gen_pars   ,
                            more_vars  = more_vars  ,
                            silent     = silent     ,
                            progress   = progress   )
                          
        self.gen_pdf    = self.pdf 
        self.fit_pdf    = fit_pdf
        self.gen_pars   = self.init_pars 
        self.fit_pars   = fit_pars

    ## the actual processing 
    def process ( self , jobid , nToys ) :

        import ROOT
        from ostap.logger.logger import logWarning
        with logWarning() :
            import ostap.core.pyrouts            
            import ostap.fitting.roofit            
            import ostap.fitting.dataset            
            import ostap.fitting.roofitresult            
            import ostap.fitting.variables
            
        from   ostap.core.ostap_types import integer_types 
        assert isinstance ( nToys , integer_types ) and 0 < nToys,\
               'Jobid %s: Invalid "nToys" argument %s/%s' % ( jobid , nToys , type ( nToys ) )
        
        import ostap.fitting.toys as Toys 
        results , stats = Toys.make_toys2 ( gen_pdf    = self.gen_pdf    ,
                                            fit_pdf    = self.fit_pdf    ,
                                            nToys      = nToys           ,
                                            data       = self.data       ,
                                            gen_config = self.gen_config , 
                                            fit_config = self.fit_config , 
                                            gen_pars   = self.gen_pars   ,
                                            fit_pars   = self.fit_pars   ,             
                                            more_vars  = self.more_vars  ,
                                            silent     = self.silent     ,
                                            progress   = self.progress   )
                
        self.the_output = results , stats
        
        return self.results () 

# ===================================================================================
## Run fitting toys in parallel
# @code
# pdf = ...
# results , stats = parallel_toys ( pdf , ## PDF  to use 
#     nToys      = 100000      ,          ## total number of toys
#     nSplit     = 100         ,          ## split into <code>nSplit</code> subjobs 
#     data       = [ 'mass' ]  ,          ## variables in dataset 
#     gen_config = { 'nEvents' : 5000 } , ## configuration of <code>pdf.generate</code>
#     fit_config = { 'ncpus'   : 2    } , ## configuration of <code>pdf.fitTo</code>
#     init_pars  = { 'mean' : 0.0 , 'sigma' : 1.0 } ) ## parameters to use for generation 
#  @endcode
#
# Derived parameters can be also   retrived via <code>more_vars</code> argument:
# @code
# ratio     = lambda res,pdf : res.ratio('x','y')
# more_vars = { 'Ratio' : ratio }
# r,  s = parallel_toys ( .... , more_vars = more_vars , ... ) 
# @endcode
#
# Parallelization is controlled by  two arguments
#  - <code>ncpus</code>, number of local cpus to use,
#   default is <code>'autodetect'</code>, that means - use all local processors
#  - <code>ppservers</code>,  list of serevers to be used (for parallel python)
#   
# @see ostap.fitting.toys
# @see ostap.fitting.toys.make_toys
# @param pdf        PDF to be used for generation and fitting
# @param nToys      total number    of pseudoexperiments to generate
# @param nSplit     split the total number of presudoexperiemtns into <code>nSplit</code> subjobs 
# @param data       variable list of variables to be used for dataset generation
# @param gen_config configuration of <code>pdf.generate</code>
# @param fit_config configuration of <code>pdf.fitTo</code>
# @param init_pars  redefine these parameters for each pseudoexperiment
# @param more_vars  calculate more variables form fit-result 
# @param silent     silent toys?
# @return dictionary with fit results for the toys and the dictionary of statistics  
def parallel_toys ( pdf                       ,
                    nToys                     , ## total number of toys 
                    nSplit                    , ## split into  <code>nSplit</code> subjobs 
                    data                      , ## template for dataset/variables 
                    gen_config                , ## parameters for <code>pdf.generate</code>   
                    fit_config = {}           , ## parameters for <code>pdf.fitTo</code>
                    init_pars  = {}           ,
                    more_vars  = {}           ,
                    silent     = True         ,
                    progress   = False        ,
                    ncpus      = 'autodetect' ,   
                    ppservers  = ()           ) :    
    """Make `ntoys` pseudoexperiments, splitting them into `nSplit` subjobs
    to be executed in parallel
    
    For each experiment:

    1. generate dataset using `pdf` with variables specified
    in `data` and configuration specified via `gen_config`
    for each generation the parameters of `pdf` are reset
    for their initial values and valeus from `init_pars`
    
    2. fit generated dataset  with `pdf` using configuration
    specified via  `fit_config`

    - pdf        PDF to be used for generation and fitting
    - nToys      total number    of pseudoexperiments to generate
    - nSplit     split total number of pseudoexperiments into `nSplit` subjobs  
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
    >>> r,  s = parallel_toys ( .... , more_vars = more_vars , ... ) 

    Parallelization is controlled by  two arguments
    - `ncpus` :  number of local cpus to use, default is `'autodetect'`,
    that means all local processors
    - `ppservers`:  list of serevers to be used (for parallel python)
    
    """
    from   ostap.core.ostap_types import integer_types 

    assert gen_config and 'nEvents' in gen_config,\
           'Number of events per toy must be specified via "gen_config" %s' % gen_config
    
    assert isinstance ( nToys  , integer_types ) and 0 < nToys  ,\
               'Jobid %s: Invalid "nToys"  argument %s/%s' % ( jobid , nToys  , type ( nToys  ) )

    
    assert isinstance ( nSplit , integer_types ) and 0 < nSplit ,\
               'Jobid %s: Invalid "nSplit" argument %s/%s' % ( jobid , nSplit , type ( nSplit ) )

    import ostap.fitting.toys as Toys
    if 1 == nSplit :
        return Toys.make_toys ( pdf        = pdf        ,
                                nToys      = nToys      ,
                                data       = data       ,
                                gen_config = gen_config ,
                                fit_config = fit_config ,
                                init_pars  = init_pars  ,
                                more_vars  = more_vars  ,
                                silent     = silent     ,
                                progress   = progress   )
    
        
    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult

    params = pdf.params() 
    toy_data = [] 
    if  isinstance ( data , ROOT.RooAbsData ) :
        varset =   data.varset()
        for v in varset : toy_data.append ( v.GetName() )
    else :
        for v in  data :
            if   isinstance ( v , ROOT.RooAbsArg )                  : toy_data.append ( v.GetName() )
            elif isinstance ( v , string_types   ) and v in  params : toy_data.append ( v           )
            else :
                raise TypeError ( "Invalid type of variable %s/%s" % ( v , type ( v ) ) )

    toy_init_pars = Toys.vars_transform ( init_pars )
    
    # ========================================================================

    if nToys <= nSplit :
        nToy   = 1
        nSplit = nToys
        nRest  = 0     
    else :
        nToy , nRest = divmod ( nToys , nSplit )
    
    task  = ToysTask    ( pdf        = pdf            ,
                          data       = toy_data       ,
                          gen_config = gen_config     ,
                          fit_config = fit_config     ,
                          init_pars  = toy_init_pars  ,
                          more_vars  = more_vars      ,
                          silent     = silent         ,
                          progress   = progress       )
                          
    wmgr  = WorkManager ( ncpus = ncpus , ppservers  = ppservers , silent = False )

    data  = nSplit * [ nToy ]
    if nRest : data.append ( nRest )
    
    wmgr.process( task , data )

    results , stats = task.results () 
    Toys.print_stats ( stats , nToys ) 
        
    return results, stats   

# ===================================================================================
## Run fitting toys in parallel
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
# Parallelization is controlled by  two arguments
#  - <code>ncpus</code>, number of local cpus to use,
#   default is <code>'autodetect'</code>, that means - use all local processors
#  - <code>ppservers</code>,  list of serevers to be used (for parallel python)
# 
# @see ostap.fitting.toys
# @see ostap.fitting.toys.make_toys2
#
# @param gen_pdf    PDF to be used for generation 
# @param fit_pdf    PDF to be used for fitting
# @param nToys      total number    of pseudoexperiments to generate
# @param nSplit     split the total number of presudoexperiemtns into <code>nSplit</code> subjobs 
# @param data       variable list of variables to be used for dataset generation
# @param gen_config configuration of <code>pdf.generate</code>
# @param fit_config configuration of <code>pdf.fitTo</code>
# @param gen_pars   redefine these parameters for each pseudoexperiment
# @param fit_pars   redefine these parameters for each pseudoexperiment
# @param silent     silent toys?
# @return dictionary with fit results for the toys and the dictionary of statistics  
def parallel_toys2(
    gen_pdf                   , ## PDF to generate toys 
    fit_pdf                   , ## PDF to generate toys 
    nToys                     , ## total number of toys 
    nSplit                    , ## split into  <code>nSplit</code> subjobs 
    data                      , ## template for dataset/variables 
    gen_config                , ## parameters for <code>pdf.generate</code>   
    fit_config = {}           , ## parameters for <code>pdf.fitTo</code>
    gen_pars   = {}           ,
    fit_pars   = {}           ,
    more_vars  = {}           ,
    silent     = True         ,
    progress   = False        ,
    ncpus      = 'autodetect' ,   
    ppservers  = ()           ) :    
    """Make `ntoys` pseudoexperiments, splitting them into `nSplit` subjobs
    to be executed in parallel
    
    For each experiment:

    1. generate dataset using `pdf` with variables specified
    in `data` and configuration specified via `gen_config`
    for each generation the parameters of `pdf` are reset
    for their initial values and valeus from `init_pars`
    
    2. fit generated dataset  with `pdf` using configuration
    specified via  `fit_config`

    - pdf        PDF to be used for generation and fitting
    - nToys      total number    of pseudoexperiments to generate
    - nSplit     split total number of pseudoexperiments into `nSplit` subjobs  
    - data       variable list of variables to be used for dataset generation
    - gen_config configuration of <code>pdf.generate</code>
    - fit_config configuration of <code>pdf.fitTo</code>
    - gen_pars   redefine these parameters for generation of  each pseudoexperiment
    - fit_pars   redefine these parameters for fitting of each pseudoexperiment
    - more_vars  dictionary of functions to define the additional results 
    - silent     silent toys?
    - progress   show progress bar? 
    
    It returns a dictionary with fit results for the toys
    
    >>> pdf = ...
    ... results, stats = parallel_toys2 (
    ...    gen_pdf    = gen_pdf     , ## PDF  to generate toys 
    ...    fit_pdf    = gen_pdf     , ## PDF  to fit toys  
    ...    nToys      = 100000      , ## total number of toys
    ...    nSplit     = 100         , ## split them into `nSplit` subjobs 
    ...    data       = [ 'mass' ]  , ## varibales in dataset 
    ...    gen_config = { 'nEvents' : 5000 } , ## configuration of `pdf.generate`
    ...    fit_config = { 'ncpus'   : 2    } , ## configuration of `pdf.fitTo`
    ...    gen_pars   = { 'mean'  : 0.0 , 'sigma'  : 1.0 } ## parameters to use for generation 
    ...    fit_pars   = { 'meanG' : 0.0 , 'sigmaG' : 1.0 } ## parameters to use for fitting
    ...   )

    Derived parameters can be also retrived via <code>more_vars</code> argument:
    >>> ratio    = lambda res,pdf : res.ratio('x','y') 
    >>> more_vars = { 'Ratio' : ratio }
    >>> r,  s = parallel_toys2 ( .... , more_vars = more_vars , ... ) 

    Parallelization is controlled by  two arguments
    - `ncpus` :  number of local cpus to use, default is `'autodetect'`,
    that means all local processors
    - `ppservers`:  list of serevers to be used (for parallel python)
    
    """
    from   ostap.core.ostap_types import integer_types 

    assert gen_config and 'nEvents' in gen_config,\
           'Number of events per toy must be specified via "gen_config" %s' % gen_config
    
    assert isinstance ( nToys  , integer_types ) and 0 < nToys  ,\
               'Jobid %s: Invalid "nToys"  argument %s/%s' % ( jobid , nToys  , type ( nToys  ) )

    assert isinstance ( nSplit , integer_types ) and 0 < nSplit ,\
               'Jobid %s: Invalid "nSplit" argument %s/%s' % ( jobid , nSplit , type ( nSplit ) )

    import ostap.fitting.toys as Toys
    if 1 == nSplit :
        return Toys.make_toys2 (
            gen_pdf    = gen_pdf    ,
            fit_pdf    = fit_pdf    ,
            nToys      = nToys      ,
            data       = data       ,
            gen_config = gen_config ,
            fit_config = fit_config ,
            gen_pars   = gen_pars   ,
            fit_pars   = fit_pars   ,
            more_vars  = more_vars  ,
            silent     = silent     ,
            progress   = progress   )
        
    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult

    params = gen_pdf.params () 
    toy_data = [] 
    if  isinstance ( data , ROOT.RooAbsData ) :
        varset =   data.varset()
        for v in varset : toy_data.append ( v.GetName() )
    else :
        for v in  data :
            if   isinstance ( v , ROOT.RooAbsArg )                  : toy_data.append ( v.GetName() )
            elif isinstance ( v , string_types   ) and v in  params : toy_data.append ( v           )
            else :
                raise TypeError ( "Invalid type of variable %s/%s" % ( v , type ( v ) ) )

    gen_init_pars = Toys.vars_transform ( gen_pars )
    fit_init_pars = Toys.vars_transform ( fit_pars )
    
    # ========================================================================

    if nToys <= nSplit :
        nToy   = 1
        nSplit = nToys
        nRest  = 0     
    else :
        nToy , nRest = divmod ( nToys , nSplit )
    
    task  = ToysTask2   ( gen_pdf    = gen_pdf        ,
                          fit_pdf    = fit_pdf        ,
                          data       = toy_data       ,
                          gen_config = gen_config     ,
                          fit_config = fit_config     ,
                          gen_pars   = gen_init_pars  ,
                          fit_pars   = fit_init_pars  ,
                          more_vars  = more_vars      ,
                          silent     = silent         ,
                          progress   = progress       )

    wmgr  = WorkManager ( ncpus = ncpus , ppservers  = ppservers , silent = False )

    data  = nSplit * [ nToy ]
    if nRest : data.append ( nRest )
    
    wmgr.process( task , data )

    results , stats = task.results () 
    Toys.print_stats ( stats , nToys ) 

    return results, stats   

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
