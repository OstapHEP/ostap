#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/parallel_statvar.py
#  (parallel) Use of stat-var functions for looong TChains 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
# =============================================================================
""" (parallel) Use of stat-var functions for looong TChains
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'parallel_statistic'  ,
    'parallel_get_stat'   ,
    'parallel_project'    , 
    'parallel_size'       ,
    'parallel_covariance' ,
    'parallel_slice'      ,
    'parallel_sum'        , 
    ## 
    'CHUNK_SIZE'          , ## maximal number of events in singel chunk 
    'MAX_FILES'           , ## maximal number of files
) 
# =============================================================================
from   ostap.math.base         import ( FIRST_ENTRY ,
                                        LAST_ENTRY  , 
                                        evt_range   ,
                                        all_entries )  
from   ostap.parallel.parallel import Task, WorkManager
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.statvar' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
CHUNK_SIZE = 50000
MAX_FILES  = 1 
## The simple task object collect statistics for loooooong chains 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
class StatVarTask(Task) :
    """ The simple task object collect statistics for loooooong chains 
    """
    ## constructor: histogram 
    def __init__ ( self , what , cuts = ""  , *args , **kwargs ) :
        """ Constructor        
        >>> task  = StatVarTask ( 'mass' , 'pt>0') 
        """
        self.what   = what
        self.cuts   = cuts  
        self.kwargs = {}
        self.args   = args 
        self.kwargs.update ( kwargs )        
        self.__output = None
        
    ## local initialization (executed once in parent process)
    def initialize_local   ( self ) :
        """ Local initialization (executed once in parent process)
        """
        self.__output = None

    # =============================================================
    ## the actual processing
    #   `params' is assumed to be a tuple/list :
    #  - the file name
    #  - the tree name in the file
    #  - the variable/expression/expression list of quantities to project
    #  - the selection/weighting criteria 
    #  - the first entry in tree to process
    #  - number of entries to process
    def process ( self , jobid , item ) :
        """ The actual processing
        `params' is assumed to be a tuple-like entity:
        - the file name
        - the tree name in the file
        - the variable/expression/expression list of quantities to project
        - the selection/weighting criteria 
        - the first entry in tree to process
        - number of entries to process
        """
        from ostap.logger.utils import logWarning
        with logWarning() :
            import ROOT
            import ostap.core.pyrouts 
            import ostap.trees.trees
            
        chain = item.chain 
        first = item.first
        last  = item.last  
        
        from   ostap.stats.statvars import data_statistic
        self.__output = data_statistic ( chain         ,
                                         self.what     ,
                                         self.cuts     , first , last,
                                         *self.args    , 
                                         **self.kwargs )
        
        return self.__output 
        
    ## merge results 
    def merge_results ( self , result , jobid = -1 ) :
        
        from ostap.stats.counters   import WSE
        from ostap.core.ostap_types import dictlike_types
        
        if not self.__output : self.__output = result
        else               :
            assert type ( self.__output ) == type ( result ) , 'Invalid types for merging!'
            if isinstance ( self.__output , dictlike_types ) :
                for key in result : 
                    if    key in self.output : self.__output [ key ] += result [ key ]
                    else                     : self.__output [ key ]  = result [ key ] 
            else :
                self.__output += result

    ## get the results 
    def results ( self ) : return self.__output 

# ================================================================================
## The simple task object collect statistics for loooooong chains 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
class GetStatTask(Task) :
    """ The simple task object collect statistics for loooooong chains 
    """
    ## constructor: histogram 
    def __init__ ( self             ,
                   target           ,
                   expressions      ,
                   cuts        = '' , *args ,**kwargs ) :
        """ Constructor        
        >>> task  = StatVarTask ( 'mass' , 'pt>0') 
        """
        self.what   = expressions
        self.cuts   = "" 
        self.kwargs = {}
        self.args   = () 
        self.kwargs.update ( kwargs )
        self.__target = target
        self.__output = None 

    # =============================================================
    ## the actual processing
    #   `params' is assumed to be a tuple/list :
    #  - the file name
    #  - the tree name in the file
    #  - the variable/expression/expression list of quantities to project
    #  - the selection/weighting criteria 
    #  - the first entry in tree to process
    #  - number of entries to process
    def process ( self , jobid , item ) :
        """ The actual processing
        `params' is assumed to be a tuple-like entity:
        - the file name
        - the tree name in the file
        - the variable/expression/expression list of quantities to project
        - the selection/weighting criteria 
        - the first entry in tree to process
        - number of entries to process
        """
        from ostap.logger.utils import logWarning
        with logWarning() :
            from   ostap.stats.statvars import data_get_stat, target_copy, target_reset 
            import ostap.trees.trees
            import ROOT
            
        chain = item.chain 
        first = item.first
        last  = item.last 

        # ================================
        target = target_copy ( self.__target )
        target_reset ( target ) 
        # ================================
        
        self.__output = data_get_stat ( chain         ,
                                        target        ,
                                        self.what     ,
                                        self.cuts     , first , last ,
                                        *self.args    , 
                                        **self.kwargs )
        return self.__output 
    
    ## merge results 
    def merge_results ( self , result , jobid = -1 ) :        
        if not self.__output : self.__output  = result
        else                 : self.__output += result
        
    ## get the results 
    def results ( self ) : return self.__output 

# ================================================================================
## The simple task object collect statistics for loooooong chains 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
class ProjectTask(Task) :
    """ The simple task object collect statistics for loooooong chains 
    """
    ## constructor: histogram 
    def __init__ ( self             ,
                   target           , 
                   expressions      ,
                   cuts        = "" ,
                   *args            , 
                   **kwargs         ) :
        """ Constructor        
        >>> task  = StatVarTask ( 'mass' , 'pt>0') 
        """
        self.what   = expressions
        self.cuts   = cuts 
        self.kwargs = {}
        self.args   = args 
        self.kwargs.update ( kwargs )
        self.__target = target
        self.__output = None 

    # =============================================================
    ## the actual processing
    #   `params' is assumed to be a tuple/list :
    #  - the file name
    #  - the tree name in the file
    #  - the variable/expression/expression list of quantities to project
    #  - the selection/weighting criteria 
    #  - the first entry in tree to process
    #  - number of entries to process
    def process ( self , jobid , item ) :
        """ The actual processing
        `params' is assumed to be a tuple-like entity:
        - the file name
        - the tree name in the file
        - the variable/expression/expression list of quantities to project
        - the selection/weighting criteria 
        - the first entry in tree to process
        - number of entries to process
        """
        from ostap.logger.utils import logWarning
        with logWarning() :
            from   ostap.stats.statvars import data_project, target_copy, targer_reset
            import ostap.trees.trees
            
        chain = item.chain 
        first = item.first
        last  = item.last 

        # ================================
        target = target_copy ( self.__target )
        target_reset ( target )
        # ================================
        
        self.__output = data_project ( chain         ,
                                       target        ,
                                       self.what     ,
                                       self.cuts     , first , last ,
                                       *self.args    , 
                                       **self.kwargs )
        return self.__output 
    
    ## merge results 
    def merge_results ( self , result , jobid = -1 ) :        
        if not self.__output : self.__output  = result
        else                 : self.__output += result
        
    ## get the results 
    def results ( self ) : return self.__output 


# ================================================================================
## The simple task object collect number of good entries for loooong chains 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
class SizeTask(Task) :
    """ The simple task object collect statistics for loooooong chains 
    """
    ## constructor: histogram 
    def __init__ ( self        ,
                   cuts   = '' ,
                   *args       , 
                   **kwargs    ) :
        """ Constructor        
        >>> task  = StatVarTask ( 'mass' , 'pt>0') 
        """
        self.cuts     = cuts 
        self.kwargs   = {}
        self.args     = args 
        self.kwargs.update ( kwargs )
        self.__output = 0 

    # =============================================================
    ## the actual processing
    #   `params' is assumed to be a tuple/list :
    #  - the file name
    #  - the tree name in the file
    #  - the variable/expression/expression list of quantities to project
    #  - the selection/weighting criteria 
    #  - the first entry in tree to process
    #  - number of entries to process
    def process ( self , jobid , item ) :
        """ The actual processing
        `params' is assumed to be a tuple-like entity:
        - the file name
        - the tree name in the file
        - the variable/expression/expression list of quantities to project
        - the selection/weighting criteria 
        - the first entry in tree to process
        - number of entries to process
        """
        from ostap.logger.utils import logWarning
        with logWarning() :
            from   ostap.stats.statvars import data_size 
            import ostap.trees.trees
            
        chain = item.chain 
        first = item.first
        last  = item.last 
        
        self.__output = data_size ( chain , self.cuts , first , last , *self.args , **self.kwargs )
        return self.__output 
    
    ## merge results 
    def merge_results ( self , result , jobid = -1 ) :        
        if not self.__output : self.__output  = result
        else                 : self.__output += result
        
    ## get the results 
    def results ( self ) : return self.__output 
    
# ================================================================================
## The simple task object collect number of good entries for loooong chains 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
class CovTask(Task) :
    """ The simple task object collect statistics for loooooong chains 
    """
    ## constructor: histogram 
    def __init__ ( self        ,
                   what        ,
                   cuts   = "" ,
                   *args       , 
                   **kwargs    ) :
        """ Constructor        
        >>> task  = StatVarTask ( 'mass' , 'pt>0') 
        """
        self.what     = what
        self.cuts     = cuts 
        self.kwargs   = {}
        self.args     = args 
        self.kwargs.update ( kwargs )
        self.__output = None 
        
    # =============================================================
    ## the actual processing
    #   `params' is assumed to be a tuple/list :
    #  - the file name
    #  - the tree name in the file
    #  - the variable/expression/expression list of quantities to project
    #  - the selection/weighting criteria 
    #  - the first entry in tree to process
    #  - number of entries to process
    def process ( self , jobid , item ) :
        """ The actual processing
        `params' is assumed to be a tuple-like entity:
        - the file name
        - the tree name in the file
        - the variable/expression/expression list of quantities to project
        - the selection/weighting criteria 
        - the first entry in tree to process
        - number of entries to process
        """
        from ostap.logger.utils import logWarning
        with logWarning() :
            from   ostap.stats.statvars import data_covariance 
            import ostap.trees.trees
            
        chain = item.chain 
        first = item.first
        last  = item.last 
        
        self.__output = data_covariance  ( chain         ,
                                           self.what     ,
                                           self.cuts     , first , last ,
                                           *self.args    , 
                                           **self.kwargs )
        return self.__output 
    
    ## merge results 
    def merge_results ( self , result , jobid = -1 ) :        
        if not self.__output : self.__output  = result
        else                 : self.__output += result
        
    ## get the results 
    def results ( self ) : return self.__output 


# ================================================================================
## The simple task object to collecy slices 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
class SliceTask(Task) :
    """ The simple task object collect slices 
    """
    ## constructor: histogram 
    def __init__ ( self        ,
                   what        ,
                   cuts   = "" ,
                   *args       , 
                   **kwargs    ) :
        """ Constructor        
        >>> task  = StatVarTask ( 'mass' , 'pt>0') 
        """
        self.what     = what
        self.cuts     = cuts 
        self.kwargs   = {}
        self.args     = () 
        self.kwargs.update ( kwargs )
        self.__output = None 
        
    # =============================================================
    ## the actual processing
    #   `params' is assumed to be a tuple/list :
    #  - the file name
    #  - the tree name in the file
    #  - the variable/expression/expression list of quantities to project
    #  - the selection/weighting criteria 
    #  - the first entry in tree to process
    #  - number of entries to process
    def process ( self , jobid , item ) :
        """ The actual processing
        `params' is assumed to be a tuple-like entity:
        - the file name
        - the tree name in the file
        - the variable/expression/expression list of quantities to project
        - the selection/weighting criteria 
        - the first entry in tree to process
        - number of entries to process
        """
        from ostap.logger.utils import logWarning
        with logWarning() :
            from   ostap.stats.statvars import data_slice 
            import ostap.trees.trees
            
        chain = item.chain 
        first = item.first
        last  = item.last 
        
        self.__output = data_slice ( chain         ,
                                     self.what     ,
                                     self.cuts     , first , last ,
                                     *self.args    , 
                                     **self.kwargs ) 
        return self.__output 
    
    ## merge results 
    def merge_results ( self , result , jobid = -1 ) :        
        if not self.__output : 
            self.__output  = result
            return
        ## acgtual  merging
        v1 , w1 = self.__output
        v2 , w2 =        result

        import numpy
        
        if   w1 is None and w2 is None : ww = None 
        elif w1 is None                :
            w1 = numpy.ones ( len  ( v1 ) , dtype = numpy.float64 )
            ww = numpy.concatenate ( ( w1 , w2 ) )
        elif w2 is None                :
            w2 = numpy.ones ( len  ( v2 ) , dtype = numpy.float64 )
            ww = numpy.concatenate ( ( w1 , w2 ) )
        else :
            ww = numpy.concatenate ( ( w1 , w2 ) )
                
        vv = numpy.concatenate ( ( v1 , v2 ) )            
            
        self.__output = vv , ww 
        
    ## get the results 
    def results ( self ) : return self.__output 
    
    
# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain  = ...
#  result = parallel_statistic ( .... ) 
#  @endcode
def parallel_statistic ( chain                    ,
                         expressions              ,
                         cuts       = ''          , *args , 
                         as_weight  = True        , 
                         use_frame  = True        , 
                         progress   = False       ,
                         chunk_size = CHUNK_SIZE  ,
                         max_files  = MAX_FILES   ,
                         silent     = True        , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.pstatVar( 'mass' , 'pt>1') 
    """
    ## few special/trivial cases

    from ostap.stats.statvars import data_statistic
    
    ## adjust first/last 
    first, last = evt_range ( chain , *args[:2] )
    
    ## number of events  to proecedd 
    nevents     = last - first 
    
    if nevents <= chunk_size :
        return data_statistic ( chain       ,
                                expressions ,
                                cuts        , first     , last , *args [2:] , 
                                as_weight   = as_weight ,
                                progress    = progress  ,
                                use_frame   = use_frame ,
                                parallel    = False     )

    ## The Task
    task   = StatVarTask ( expressions            ,
                           cuts      , *args[2:]  , 
                           as_weight = True       ,
                           progress  = False      ,
                           use_frame = use_frame  , 
                           parallel   = False     ) 

    ## Manager 
    wmgr   = WorkManager ( silent = silent , progress = progress or not silent , **kwargs )

    ## split data 
    from ostap.trees.utils   import Chain
    ch     = Chain    ( chain , first = first   , last = last )
    trees  = ch.split ( chunk_size = chunk_size , max_files = max_files )

    wmgr.process ( task , trees )

    del trees
    del ch    

    results = task.results()
    
    return results 

# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain          = ...
#  chain.pSumVar ( .... ) 
#  @endcode 
def parallel_sum ( chain                    ,
                   expressions              , 
                   cuts       = ''          , *args , 
                   as_weight  = True        ,
                   progress   = False       ,
                   use_frame  = True        , 
                   chunk_size = 100000      ,
                   max_files  =  1          ,
                   silent     = True        , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.psumVar( 'mass' , 'pt>1') 
    """
    ## few special/trivial cases
    
    result = parallel_statistic ( chain                    ,
                                  expressions              ,
                                  cuts        , *args      , 
                                  as_weight   = as_weight  ,
                                  progress    = progress   ,
                                  use_frame   = use_frame  , 
                                  chunk_size  = chunk_size , 
                                  max_files   = max_files  ,
                                  silent      = silent     , **kwargs )
    
    if isinstance ( result , dict ) :
        sumres = {}
        for key in result  :
            r = result [ key ] 
            sumres [ key ] = VE ( r.sum() , r.sum2() )
        return sumres
    
    return VE ( result.sum() , result.sum2() )


# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain          = ...
#  chain.pStatVar ( .... ) 
#  @endcode
def parallel_size ( chain                    ,  
                    cuts       = ''          , *args  , 
                    as_weight  = True        , 
                    progress   = False       ,
                    use_frame  = True        , 
                    chunk_size = 100000      ,
                    max_files  = 1           ,
                    silent     = True        , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.pstatVar( 'mass' , 'pt>1') 
    """
    ## few special/trivial cases

    import ostap.trees.trees
    from   ostap.stats.statvars import data_size 
    
    ## adjust first/last 
    first, last = evt_range ( chain , *args[:2] )
    
    ## number of events to proecedd 
    nevents     = last - first 
    
    if nevents <= chunk_size :
        return data_size ( chain ,
                           cuts      , *args     , 
                           progress  = progress  ,
                           use_frame = use_frame ,
                           parallel  = False     )

    ## The Task 
    task   = SizeTask ( cuts      , *args [2:] ,
                        use_frame = use_frame  , 
                        progress  = False      )

    ## Manager 
    wmgr   = WorkManager ( silent = silent , progress = progress or not silent , **kwargs )

    ## split data 
    from ostap.trees.utils import Chain
    ch     = Chain ( chain , first = first , last = last  )
    trees  = ch.split ( chunk_size = chunk_size , max_files = max_files )

    
    wmgr.process ( task , trees )

    del trees
    del ch    

    results = task.results()
    
    return results 

# ===================================================================================
## parallel processing of loooong chain/tree 
def parallel_get_stat ( chain                    ,
                        target                   , 
                        expressions              ,
                        cuts       = ''          , *args , 
                        progress   = False       ,
                        use_frame  = True        , 
                        chunk_size = 100000      ,
                        max_files  = 1           ,
                        silent     = True        ,  **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.pstatVar( 'mass' , 'pt>1') 
    """
    ## few special/trivial cases

    import ostap.trees.trees
    from   ostap.stats.statvars import data_get_stat 

    first , last = evt_range ( chain , *args[:2]  )
    
    nevents = last - first 
    
    if nevents <= chunk_size :
        return data_get_stat  ( chain       ,
                                target      , 
                                expressions ,
                                cuts        , first     , last , *args[2:] , 
                                progress    = progress  ,
                                use_frame   = use_frame ,
                                parallel    = False     )
    
    ## The Task
    task   = GetStatTask ( target                ,
                           expressions           , 
                           cuts      , *args[2:] ,
                           progress  = False     ,
                           use_frame = use_frame ,
                           parallel  = False     ) 

    ## Manager 
    wmgr   = WorkManager ( silent = silent , progress = progress or not silent , **kwargs )

    ## split data 
    from ostap.trees.utils import Chain
    ch     = Chain    ( chain , first = first , last = last  )
    trees  = ch.split ( chunk_size = chunk_size , max_files = max_files )

    wmgr.process ( task , trees )

    del trees
    del ch    

    results = task.results()
    
    return results 

# ===================================================================================
## parallel processing of loooong chain/tree 
def parallel_project ( chain                    ,
                       target                   , 
                       expressions              ,  
                       cuts       = ''          , *args , 
                       as_weight  = True        , 
                       progress   = False       ,
                       use_frame  = False       , 
                       chunk_size = 100000      ,
                       max_files  = 1           ,
                       silent     = True        , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.pstatVar( 'mass' , 'pt>1') 
    """
    ## few special/trivial cases
    
    import ostap.trees.trees
    from   ostap.stats.statvars import data_project

    first , last = evt_range ( chain , *args[:2] ) 
    nevents = last - first
    
    if nevents < chunk_size :
        return data_project ( chain       ,
                              target      , 
                              expressions ,
                              cuts        , first     , last , *args[2:] , 
                              as_weight   = as_weight , 
                              progress    = progress  ,
                              use_frame   = use_frame ,
                              parallel    = False     )
    ## The Task 
    task   = ProjectTask ( target                 ,
                           expressions            , 
                           cuts      , *args [2:] ,
                           as_weight = as_weight  , 
                           progress  = False      ,
                           use_frame = use_frame  ,
                           parallel  = False      ) 

    ## Manager 
    wmgr   = WorkManager ( silent = silent , progress = progress or not silent , **kwargs )

    ## split data 
    from ostap.trees.utils import Chain
    ch     = Chain ( chain , first = first , last = last  )    
    trees  = ch.split ( chunk_size = chunk_size , max_files = max_files )
    
    wmgr.process ( task , trees )

    del trees
    del ch    

    results = task.results()
    
    return results 

# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain          = ...
#  chain.pStatVar ( .... ) 
#  @endcode
def parallel_covariance ( chain                    ,
                          expressions              ,
                          cuts       = ''          , *args  , 
                          as_weight  = True        , 
                          progress   = False       ,
                          use_frame  = True        , 
                          chunk_size = 100000      ,
                          max_files  = 1           ,
                          silent     = True        , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.pstatVar( 'mass' , 'pt>1') 
    """
    ## few special/trivial cases

    import ostap.trees.trees
    from   ostap.stats.statvars import data_covariance 

    ## adjust first/last 
    first, last = evt_range ( chain , *args[:2] )
    
    ## number of evetns to proecedd 
    nevents     = last - first 
    
    if nevents <= chunk_size :
        return data_covariance ( chain       ,
                                 expressions ,
                                 cuts        , first     , last , *args[2:] , 
                                 as_weight   = as_weight ,
                                 progress    = progress  ,
                                 use_frame   = use_frame ,
                                 parallel    = False     )

    ## prepare task 
    task   = CovTask ( expressions           ,
                       cuts      , *args[2:] ,
                       as_weight = as_weight ,
                       use_frame = use_frame , 
                       progress  = False     ,
                       parallel  = False     )

    ## Manager
    wmgr   = WorkManager ( silent = silent , progress = progress or not silent , **kwargs )
    
    ## split data 
    from ostap.trees.utils import Chain
    ch     = Chain    ( chain , first = first , last = last )    
    trees  = ch.split ( chunk_size = chunk_size , max_files = max_files )

    ## process 
    wmgr.process ( task , trees )

    del trees
    del ch    

    results = task.results()
    
    return results 

# ===================================================================================
## parallel processing of loooong chain/tree: get slice from TTree/TChain  
#  @code
#  chain           = ...
#  result , weight = parallel_slice ( chain , 'pt,mass', 'y<4.5' ) 
#  @endcode
def parallel_slice ( chain                    ,
                     expressions              ,
                     cuts       = ''          , *args ,
                     structured = True        ,
                     transpose  = True        , 
                     progress   = False       ,
                     use_frame  = True        , 
                     chunk_size = CHUNK_SIZE  ,
                     max_files  = MAX_FILES   ,
                     silent     = True        , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain = ...
    >>> result = parallel_slice ( chain, 'mass' , 'pt>1') 
    """
    ## few special/trivial cases

    import ostap.trees.trees
    from   ostap.stats.statvars import data_slice 
    
    ## adjust first/last 
    first, last = evt_range ( chain , *args[:2] )
    
    ## nothing to process 
    if last <= first : return () , None  
    
    ## number of events to process 
    nevents     = last - first 

    """
    if nevents <= chunk_size :
        return data_slice( chain       ,
                           expressions ,
                           cuts        , first     , last , *args[2:] , 
                           structured  = transpose ,
                           transpose   = transpose , 
                           progress    = progress  ,
                           use_frame   = use_frame ,
                           parallel    = False     )
    """
    
    ## Task 
    task   = SliceTask ( expressions             ,
                         cuts       , *args[2:]  ,
                         structured = structured ,
                         transpose  = transpose  , 
                         use_frame  = use_frame  , 
                         progress   = False      ,
                         parallel   = False      )

    ## Manager 
    wmgr   = WorkManager ( silent = silent , progress = progress or not silent , **kwargs )

    ## Split data 
    from ostap.trees.utils import Chain
    ch     = Chain    ( chain , first = first , last = last  )    
    trees  = ch.split ( chunk_size = chunk_size , max_files = max_files )
    
    wmgr.process ( task , trees )

    del trees
    del ch    

    results = task.results()
    
    return results 

# =============================================================================
_decorated_classes_ = (
    )

_new_methods_       = (
    parallel_statistic  ,
    parallel_get_stat   ,
    parallel_project    ,
    parallel_covariance ,
    parallel_slice      ,
    parallel_sum        ,
)

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
