#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/parallel_statvar.py
#  (parallel) Use of stat-var for looong TChains 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
# =============================================================================
"""(parallel) Use of stat-var for looong TChains
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'parallel_statistic' ,
    'parallel_get_stat'  ,
    'parallel_sum'       , 
) 
# =============================================================================
from   ostap.parallel.parallel import Task, WorkManager
from   ostap.core.core         import VE
from   ostap.stats.statvar     import FIRST, LAST, target_reset, target_copy  
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.statvar' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
## The simple task object collect statistics for loooooong chains 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
class StatVarTask(Task) :
    """ The simple task object collect statistics for loooooong chains 
    """
    ## constructor: histogram 
    def __init__ ( self , what , **kwargs ) :
        """ Constructor        
        >>> task  = StatVarTask ( 'mass' , 'pt>0') 
        """
        self.what     = what 
        self.kwargs = {}
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
        last  = min ( LAST , first + item.nevents if 0 < item.nevents else LAST )
        
        from   ostap.stats.statvars import data_statistic
        self.__output = data_statistic ( chain , self.what , **self.kwargs )
        
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
    def __init__ ( self        ,
                   target      ,
                   expressions ,
                   **kwargs    ) :
        """ Constructor        
        >>> task  = StatVarTask ( 'mass' , 'pt>0') 
        """
        self.wjat   = expressions 
        self.kwargs = {}
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
            from   ostap.stats.statvars import data_get_stat
            import ostap.trees.trees
            import ROOT
            
        chain = item.chain 
        first = item.first
        last  = min ( LAST , first + item.nevents if 0 < item.nevents else LAST )

        # ================================
        target = target_copy ( self.__target )
        target_reset ( target ) 
        # ================================
        
        self.__output = data_get_stat ( chain     ,
                                        target    ,
                                        self.what , **self.kwargs )        
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
    def __init__ ( self        ,
                   target      ,
                   expressions ,
                   **kwargs    ) :
        """ Constructor        
        >>> task  = StatVarTask ( 'mass' , 'pt>0') 
        """
        self.wjat   = expressions 
        self.kwargs = {}
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
            from   ostap.stats.statvars import data_project
            import ostap.trees.trees
            import ROOT
            
        chain = item.chain 
        first = item.first
        last  = min ( LAST , first + item.nevents if 0 < item.nevents else LAST )

        # ================================
        target = target_copy ( self.__target )
        target_reset ( target )
        # ================================
        
        self.__output = data_project ( chain     ,
                                       target    ,
                                       self.what , **self.kwargs )        
        return self.__output 
    
    ## merge results 
    def merge_results ( self , result , jobid = -1 ) :        
        if not self.__output : self.__output  = result
        else                 : self.__output += result
        
    ## get the results 
    def results ( self ) : return self.__output 
    
# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain          = ...
#  chain.pStatVar ( .... ) 
#  @endcode
def parallel_statistic ( chain               ,
                         expressions         ,
                         cuts       = ''     ,
                         first      = FIRTS  , 
                         last       = LAST   ,
                         as_weight  = True   ,
                         progress   = False  ,
                         use_frame  = True   , 
                         chunk_size = 100000 ,
                         max_files  = 1      ,
                         silent     = True   ,
                         progress   = False  , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.pstatVar( 'mass' , 'pt>1') 
    """
    ## few special/trivial cases

    import ostap.trees.trees
    from   ostap.stats.statvars import data_statistic
    
    assert  isinstance ( first, integer_type ) \
        and isinstance ( last , integer_type ) and 0 <= first < last , \
        "Invalid first/last setting!"
    
    nevt = len ( chain )
    assert 0 <= first <= nevt , "Invalid first setting!"
    
    last    = min ( last , first + nevt )    
    nevents = nevt 
    
    if 0 <= first and 0 < nevents < chunk_size :
        return data_statistic ( chain       ,
                                expressions , first  , last
                                cuts        = cuts      , 
                                as_weight   = as_weight ,
                                progress    = progress  ,
                                use_frame   = use_frame ,
                                parallel    = False     )
    elif isinstance ( chain , ROOT.TChain ) and 0 < nevents < chunk_size :
        return data_statistic ( chain       ,
                                expressions , first  , last
                                cuts        = cuts      , 
                                as_weight   = as_weight ,
                                progress    = progress  ,
                                use_frame   = use_frame ,
                                parallel    = False     )
    
    from ostap.trees.trees import Chain
    ch     = Chain ( chain , first = first , nevents = nevents )

    task   = StatVarTask ( expressions           ,
                           cuts      = cuts      ,
                           as_weight = True      ,
                           use_frame = use_frame , 
                           progress  = False     )
    
    wmgr   = WorkManager ( silent = silent , progress = progress or not silent , **kwargs )
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
def parallel_sum ( chain               ,
                   expressions         ,
                   cuts       = ''     ,
                   first      = FIRTS  , 
                   last       = LAST   ,                   
                   as_weight  = True   ,
                   progress   = False  ,
                   use_frame  = True   , 
                   chunk_size = 100000 ,
                   max_files  =  1     ,
                   silent     = True   , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.psumVar( 'mass' , 'pt>1') 
    """
    ## few special/trivial cases
    
    result = parallel_statistic ( chain                   ,
                                  expressions             ,
                                  cuts       = cuts       ,
                                  first      = first      ,
                                  last       = last       ,
                                  as_weight  = as_weight  ,
                                  progress   = progress   ,
                                  use_frame  = use_frame  , 
                                  chunk_size = chunk_size , 
                                  max_files  = max_files  ,
                                  silent     = silent     , **kwargs )
    
    if isinstance ( result , dict ) :
        sumres = {}
        for key in result  :
            r = result [ key ] 
            sumres [ key ] = VE ( r.sum() , r.sum2() )
        return sumres
    
    return VE ( result.sum() , result.sum2() )


# ===================================================================================
## parallel processing of loooong chain/tree 
def parallel_get_stat ( chain               ,
                        target              , 
                        expressions         , 
                        cuts       = ''     ,
                        first      = FIRTS  , 
                        last       = LAST   ,
                        progress   = False  ,
                        use_frame  = True   , 
                        chunk_size = 100000 ,
                        max_files  = 1      ,
                        silent     = True   ,
                        progress   = False  , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.pstatVar( 'mass' , 'pt>1') 
    """
    ## few special/trivial cases

    import ostap.trees.trees
    from   ostap.stats.statvars import data_get_stat 

    assert  isinstance ( first, integer_type ) \
        and isinstance ( last , integer_type ) and 0 <= first < last , \
        "Invalid first/last setting!"
    
    nevt = len ( chain )
    assert 0 <= first <= nevt , "Invalid first setting!"
    
    last    = min ( last , first + nevt )    
    nevents = nevt 
    
    if 0 <= first and 0 < nevents < chunk_size :
        return data_get_stat  ( chain       ,
                                target      , 
                                expressions , first  , last
                                cuts        = cuts      , 
                                progress    = progress  ,
                                use_frame   = use_frame ,
                                parallel    = False     )
    
    elif isinstance ( chain , ROOT.TChain ) and 0 < nevents < chunk_size :
        return data_get_stat ( chain       ,
                               target      , 
                               expressions , first  , last
                               cuts        = cuts      , 
                               progress    = progress  ,
                               use_frame   = use_frame ,
                               parallel    = False     )

    from ostap.trees.trees import Chain
    ch     = Chain ( chain , first = first , nevents = nevents )

    task   = GetStatTask ( target                ,
                           expressions           , 
                           cuts      = cuts      ,
                           use_frame = use_frame ,
                           progress  = False     ) 
    
    wmgr   = WorkManager ( silent = silent , progress = progress or not silent , **kwargs )
    trees  = ch.split ( chunk_size = chunk_size , max_files = max_files )

    wmgr.process ( task , trees )

    del trees
    del ch    

    results = task.results()
    
    return results 

# ===================================================================================
## parallel processing of loooong chain/tree 
def parallel_project ( chain               ,
                       target              , 
                       expressions         , 
                       cuts       = ''     ,
                       first      = FIRTS  , 
                       last       = LAST   ,
                       progress   = False  ,
                       use_frame  = False  , 
                       chunk_size = 100000 ,
                       max_files  = 1      ,
                       silent     = True   ,
                       progress   = False  , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.pstatVar( 'mass' , 'pt>1') 
    """
    ## few special/trivial cases
    
    import ostap.trees.trees
    from   ostap.stats.statvars import data_project
    
    assert  isinstance ( first, integer_type ) \
        and isinstance ( last , integer_type ) and 0 <= first < last , \
        "Invalid first/last setting!"
    
    nevt = len ( chain )
    assert 0 <= first <= nevt , "Invalid first setting!"
    
    last    = min ( last , first + nevt )    
    nevents = nevt 
    
    if 0 <= first and 0 < nevents < chunk_size :
        return data_project  ( chain       ,
                               target      , 
                               expressions , first  , last
                               cuts        = cuts      , 
                               progress    = progress  ,
                               use_frame   = use_frame ,
                               parallel    = False     )
    
    elif isinstance ( chain , ROOT.TChain ) and 0 < nevents < chunk_size :
        return data_project ( chain       ,
                              target      , 
                              expressions , first  , last
                              cuts        = cuts      , 
                              progress    = progress  ,
                              use_frame   = use_frame ,
                              parallel    = False     )
    
    from ostap.trees.trees import Chain
    ch     = Chain ( chain , first = first , nevents = nevents )

    task   = ProjectTask ( target               ,
                           expressions           , 
                           cuts      = cuts      ,
                           use_frame = use_frame ,
                           progress  = False     ) 
    
    wmgr   = WorkManager ( silent = silent , progress = progress or not silent , **kwargs )
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
    'parallel_statistic' ,
    'parallel_get_stat'  ,
    'parallel_project'   ,
    'parallel_sum'       ,
)

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        
# =============================================================================
#                                                                       The END 
# =============================================================================
