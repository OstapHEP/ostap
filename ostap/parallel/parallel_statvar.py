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
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.statvar' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
import ROOT
from   ostap.parallel.parallel import Task, WorkManager
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

        print ('I AM HERE START/1', jobid, item ) 
        import ROOT
        from ostap.logger.utils import logWarning
        with logWarning() :
            import ostap.core.pyrouts 
            import ostap.trees.trees 

        print ('I AM HERE START/2', jobid, item )
        
        chain   = item.chain 
        first   = item.first
        last    = min ( n_large , first + item.nevents if 0 < item.nevents else n_large )
        
        print ('I AM HERE START/3', jobid, item )

        from ostap.trees.trees  import _stat_vars_

        print ('I AM HERE START/4', jobid , self.what , self.cuts , first , last , type ( chain )  )

        self.__output = _stat_vars_ ( chain , self.what , self.cuts , first , last )

        print ('I AM HERE END', jobid, item )
        
        return self.__output 
        
    ## merge results 
    def merge_results ( self , result , jobid = -1 ) :
        
        print ('I AM HERE MERGE/1', jobid )

        from ostap.stats.counters import WSE

        
        if not self.__output : self.__output = result
        else               :
            print ('I AM HERE MERGE/2', jobid )
            
            assert type( self.__output ) == type ( result ) , 'Invalid types for merging!'

            print ('I AM HERE MERGE/3', jobid )

            if isinstance ( self.__output , dict ) :
                print ('I AM HERE MERGE/4', jobid , result.keys() )

                for key in result : 
                    if    key in self.output : self.__output [ key ] += result [ key ]
                    else                     : self.__output [ key ]  = result [ key ] 
            else :
                
                print ('I AM HERE MERGE/5', jobid )
                
                self.__output += result

            print ('I AM HERE MERGE/6', jobid )

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
               chunk_size = 250000 ,
               max_files  =  1     ,
               silent     = True   , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.pstatVar( 'mass' , 'pt>1') 
    """

    print ( 'I am pStatVar/1' ) 
    ## few special/trivial cases

    last = min ( n_large , first + nevents if 0 < nevents else n_large )
    
    if 0 <= first and 0 < nevents < chunk_size :
        print ( 'I am pStatVar/1.1' ) 
        return chain.statVar ( what , cuts , first , last )
    elif isinstance ( chain , ROOT.TChain ) :
        print ( 'I am pStatVar/1.2' ) 
        if chain.nFiles() < 3 and len ( chain ) < chunk_size :
            print ( 'I am pStatVar/1.3' ) 
            return chain.statVar ( what , cuts , first , last )                         
    elif isinstance ( chain , ROOT.TTree  ) and len ( chain ) < chunk_size :
        print ( 'I am pStatVar/1.4' ) 
        return chain.statVar ( what , cuts , first , last ) 

    print ( 'I am pStatVar/2' ) 

    from ostap.trees.trees import Chain
    ch     = Chain ( chain , first = first , nevents = nevents )

    task   = StatVarTask ( what , cuts )
    wmgr   = WorkManager ( silent = silent , **kwargs )

    trees  = ch.split ( chunk_size = chunk_size , max_files = max_files )

    print ( 'I am pStatVar/3' , len(trees)  ) 

    wmgr.process ( task , trees )

    print ( 'I am pStatVar/4' , len(trees)  ) 

    del trees
    del ch    

    results = task.results()
    
    return results 

ROOT.TChain.pstatVar = pStatVar 
ROOT.TTree .pstatVar = pStatVar

# =============================================================================
_decorated_classes_ = (
    ROOT.TTree  ,
    ROOT.TChain ,    
    )

_new_methods_       = (
    ROOT.TTree .pstatVar,
    ROOT.TChain.pstatVar,
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
