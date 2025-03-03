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
    'pStatVar'   , ## get the statistics from loooong TChain in paralell
    ) 
# =============================================================================
from   ostap.parallel.parallel import Task, WorkManager
from   ostap.core.core         import VE
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.statvar' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
n_large = ROOT.TVirtualTreePlayer.kMaxEntries
# =============================================================================
## The simple task object collect statistics for loooooong chains 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
class StatVarTask(Task) :
    """The simple task object collect statistics for loooooong chains 
    """
    ## constructor: histogram 
    def __init__ ( self , what , cuts = '' ) :
        """Constructor        
        >>> task  = StatVarTask ( 'mass' , 'pt>0') 
        """
        self.what     = what 
        self.cuts     = str(cuts) 
        self.__output = None
        
    ## local initialization (executed once in parent process)
    def initialize_local   ( self ) :
        """Local initialization (executed once in parent process)
        """
        from ostap.stats.counters import WSE
        self.__output = None
            
    ## the actual processing
    #   ``params'' is assumed to be a tuple/list :
    #  - the file name
    #  - the tree name in the file
    #  - the variable/expression/expression list of quantities to project
    #  - the selection/weighting criteria 
    #  - the first entry in tree to process
    #  - number of entries to process
    def process ( self , jobid , item ) :
        """The actual processing
        ``params'' is assumed to be a tuple-like entity:
        - the file name
        - the tree name in the file
        - the variable/expression/expression list of quantities to project
        - the selection/weighting criteria 
        - the first entry in tree to process
        - number of entries to process
        """

        import ROOT
        from ostap.logger.utils import logWarning
        with logWarning() :
            import ostap.core.pyrouts 
            import ostap.trees.trees 
        chain   = item.chain 
        first   = item.first
        last    = min ( n_large , first + item.nevents if 0 < item.nevents else n_large )
        
        from   ostap.stats.statvars import data_statistics 
        self.__output = data_statistics ( chain , self.what , self.cuts , first , last )

        return self.__output 
        
    ## merge results 
    def merge_results ( self , result , jobid = -1 ) :
        
        from ostap.stats.counters   import WSE
        from ostap.core.ostap_types import dictlike_types
        
        if not self.__output : self.__output = result
        else               :
            assert type( self.__output ) == type ( result ) , 'Invalid types for merging!'
            if isinstance ( self.__output , dictlike_types ) :
                for key in result : 
                    if    key in self.output : self.__output [ key ] += result [ key ]
                    else                     : self.__output [ key ]  = result [ key ] 
            else :
                self.__output += result

    ## get the results 
    def results ( self ) : return self.__output 


# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain          = ...
#  chain.pStatVar ( .... ) 
#  @endcode 
def pStatVar ( chain               ,
               what                ,
               cuts       = ''     ,
               nevents    = -1     ,
               first      =  0     ,
               chunk_size = 500000 ,
               max_files  =  1     ,
               silent     = True   , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.pstatVar( 'mass' , 'pt>1') 
    """
    ## few special/trivial cases

    import ostap.trees.trees
    last = min ( n_large , first + nevents if 0 < nevents else n_large )

    if 0 <= first and 0 < nevents < chunk_size :
        return chain.statVar ( what , cuts , first , last )
    elif isinstance ( chain , ROOT.TChain ) and 1 == chain.nFiles() and len ( chain ) < chunk_size :
        return chain.statVar ( what , cuts , first , last )                         
    elif isinstance ( chain , ROOT.TTree  ) and len ( chain ) < chunk_size :
        return chain.statVar ( what , cuts , first , last )
    
    from ostap.trees.trees import Chain
    ch     = Chain ( chain , first = first , nevents = nevents )

    task   = StatVarTask ( what , cuts )
    wmgr   = WorkManager ( silent = silent , **kwargs )

    trees  = ch.split ( chunk_size = chunk_size , max_files = max_files )

    wmgr.process ( task , trees )

    del trees
    del ch    

    results = task.results()
    
    return results 

ROOT.TChain.pstatVar  = pStatVar 
ROOT.TTree .pstatVar  = pStatVar
ROOT.TChain.pstatVars = pStatVar 
ROOT.TTree .pstatVars = pStatVar


# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain          = ...
#  chain.pSumVar ( .... ) 
#  @endcode 
def pSumVar ( chain               ,
              what                ,
              cuts       = ''     ,
              nevents    = -1     ,
              first      =  0     ,
              chunk_size = 500000 ,
              max_files  =  1     ,
              silent     = True   , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.psumVar( 'mass' , 'pt>1') 
    """
    ## few special/trivial cases

    result = pStatVar ( chain      = chain      ,
                        what       = what       ,
                        cuts       = cuts       ,
                        nevents    = nevents    ,
                        first      = first      ,
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


ROOT.TChain.psumVar  = pSumVar 
ROOT.TTree .psumVar  = pSumVar

# =============================================================================
_decorated_classes_ = (
    ROOT.TTree  ,
    ROOT.TChain ,    
    )

_new_methods_       = (
    ROOT.TTree .pstatVar  ,
    ROOT.TChain.pstatVar  ,
    ROOT.TTree .pstatVars ,
    ROOT.TChain.pstatVars ,
    ROOT.TTree .psumVar   ,
    ROOT.TChain.psumVar   ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    if not ( 2**32 - 1 ) <= n_large <= ROOT.TVirtualTreePlayer.kMaxEntries :
        logger.error ( "Invalid setting of ``n_large''(%d) parameter (>%d)" % ( n_large , ROOT.TVirtualTreePlayer.kMaxEntries ) )

# =============================================================================
#                                                                       The END 
# =============================================================================
