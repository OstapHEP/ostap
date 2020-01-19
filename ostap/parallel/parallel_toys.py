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
#  @see ostap.fitting.toys 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-01-18 
class  ToysTask(Task) :
    """The simple task object for parallel fitting toys
    """
    ## 
    def __init__ ( self               ,
                   pdf                ,
                   data               ,
                   toy_config         ,
                   fit_config = {}    ,
                   init_pars  = {}    ,
                   more_vars  = {}    ,                  
                   silent     = True  ,
                   progress   = False ) :
        
        self.pdf        = pdf
                
        self.data       = [] 
        self.toy_config = toy_config 
        self.fit_config = fit_config
        self.init_pars  = init_pars 
        self.more_vars  = more_vars  

        import ostap.fitting.roofit
        import ostap.fitting.dataset 
        import ostap.fitting.variables
        import ostap.fitting.roofitresult
        
        self.silent     = silent
        self.progress   = progress 
        
        self.__output   = () 

    def initialize_local   ( self ) : self.__output = () 

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
                                           toy_config = self.toy_config , 
                                           fit_config = self.fit_config , 
                                           init_pars  = self.init_pars  ,
                                           more_vars  = self.more_vars  ,
                                           silent     = self.silent     ,
                                           progress   = self.progress   )
        
        
        self.__output = results , stats
        
        return self.__output 

    ## merge results/datasets 
    def merge_results ( self, result ) :
        
        if result :
            results , stat = result
            if not self.__output :
                self.__output = results , stat 
            else :
                results  , stat    = result
                results_ , stat_   = self.__output 

                rset = set()
                for p in results  : rset.add ( p )
                for p in results_ : rset.add ( p )                
                for p in rset : results_[ p ] += results [ p ]

                sset = set()
                for p in stat  : sset.add ( p )
                for p in stat_ : sset.add ( p )
                for p in sset  : stat_ [ p  ] += stat [ p ] 
                
                self.__output = results_ , stat_ 
        else :
            logger.error ( "No valid results for merging" )
            
    ## get the results 
    def results ( self ) :
        return self.__output


# ===================================================================================
## parallel toys 
#  @code
#  chain    = ...
#  selector =  ...
#  chain.pprocess ( selector )
#  @endcode 
def parallel_toys ( pdf                       ,
                    nToys                     , ## total number of toys 
                    nSplit                    , ## split into  <code>nSplit</code> subjobs 
                    data                      , ## template for dataset/variables 
                    toy_config                , ## parameters for <code>pdf.generate</code>   
                    fit_config = {}           , ## parameters for <code>pdf.fitTo</code>
                    init_pars  = {}           ,
                    more_vars  = {}           ,
                    silent     = True         ,
                    progress   = False        ,
                    ncpus      = 'autodetect' , 
                    ppservers  = ()           ) :
    
    """
    >>>chain    = ...
    >>> selector =  ...
    >>> chain.pprocess ( selector )
    """
    from   ostap.core.ostap_types import integer_types 

    assert toy_config and 'nEvents' in toy_config,\
           'Number of events per toy must be specified via "toy_config" %s' % toy_config
    
    assert isinstance ( nToys  , integer_types ) and 0 < nToys  ,\
               'Jobid %s: Invalid "nToys"  argument %s/%s' % ( jobid , nToys  , type ( nToys  ) )

    
    assert isinstance ( nSplit , integer_types ) and 0 < nSplit ,\
               'Jobid %s: Invalid "nSplit" argument %s/%s' % ( jobid , nSplit , type ( nSplit ) )

    if 1 == nSplit :
        import ostap.fitting.toys as Toys
        return Toys.make_toys ( pdf        = pdf        ,
                                nToys      = nToys      ,
                                data       = data       ,
                                toy_config = toy_config ,
                                fit_config = fit_config ,
                                init_pars  = init_pars  ,
                                more_vars  = more_vars  ,
                                silent     = silent     ,
                                progress   = progress   )
    
        
    import ostap.fitting.roofit
    import ostap.fitting.dataset
    import ostap.fitting.variables
    import ostap.fitting.roofitresult
    
    toy_data = [] 
    if  isinstance ( data , ROOT.RooAbsData ) :
        varset =   data.varset()
        for v in varset : toy_data.append ( v.GetName() )
    else :
        for v in  data :
            if   isinstance ( v , ROOT.RooAbsArg ) : toy_data.append ( v.GetName() )
            elif isinstance ( v , string_types   ) : toy_data.append ( v           )
            else :
                raise TypeError ( "Invalid type of variable %s/%s" % ( v , type ( v ) ) )

    toy_init_pars = {} 
    if isinstance ( init_pars , ROOT.RooFitResult ) :
        ps = init_pars.dct_params()
        for p in ps        : toy_init_pars [ p           ] = float ( ps        [ p ] )
    elif isinstance ( init_pars , dict ) :
        for p in init_pars : toy_init_pars [ p           ] = float ( init_pars [ p ] )
    else :
        for p in init_pars : toy_init_pars [ p.GetName() ] = float ( p )
                
    
    # ========================================================================

    if nToys <= nSplit :
        nToy   = 1
        nSplit = nToys
        nRest  = 0     
    else :
        nToy , nRest = divmod ( nToys , nSplit )
    
    task  = ToysTask    ( pdf        = pdf            ,
                          data       = toy_data       ,
                          toy_config = toy_config     ,
                          fit_config = fit_config     ,
                          init_pars  = toy_init_pars  ,
                          more_vars  = more_vars      ,
                          silent     = silent         ,
                          progress   = progress       )
                          
    wmgr  = WorkManager ( ncpus = ncpus , ppservers  = ppservers , silent = False )

    data  = nSplit * [ nToy ]
    if nRest : data.append ( nRest )
    
    wmgr.process( task , data )
    
    return task.results()  

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
