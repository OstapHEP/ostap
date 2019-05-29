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
        self.what  = what 
        self.cuts  = str(cuts) 
        
    ## local initialization (executed once in parent process)
    def initializeLocal   ( self ) :
        """Local initialization (executed once in parent process)
        """
        from ostap.stats.counters import WSE
        self.output = None 
        
    def _resetOutput(self):
        self.output = None 
 
    ## remote initialization (executed for each sub-processs)
    def initializeRemote  ( self ) : pass 
    
    ## the actual processing
    #   ``params'' is assumed to be a tuple/list :
    #  - the file name
    #  - the tree name in the file
    #  - the variable/expression/expression list of quantities to project
    #  - the selection/weighting criteria 
    #  - the first entry in tree to process
    #  - number of entries to process
    def process ( self , item ) :
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
        ## last    = min ( n_large , first + item.nevents if 0 < item.nevents else n_large )
        last    = n_large
        
        from ostap.trees.trees  import _stat_vars_
        self.output = _stat_vars_ ( chain , self.what , self.cuts , first , last )
        
    ## finalization (executed at the end at parent process)
    def finalize ( self ) : pass 

    ## merge results 
    def _mergeResults ( self , result ) :
        
        from ostap.stats.counters import WSE

        if not self.output : self.output = result
        else               :
            assert type( self.output ) == type ( result ) , 'Invalid types for merging!'
            if isinstance ( self.output , dict ) : 
                for key in result : 
                    if self.output.has_key ( key ) : self.output[key] += result[key]
                    else                           : self.output[key]  = result[key] 
            else :
                self.output += result
                    


# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain    = ...
#  chain.pStatVar ( .... ) 
#  @endcode 
def pStatVar ( chain        , what , cuts = ''    ,
               nevents = -1 ,
               first   =  0 , chunk_size = 100000 , max_files = 10 , ppservers = () , silent = True ) :
    """ Parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> chain.pstatVar( 'mass' , 'pt>1') 
    """

    ## few special/trivial cases

    last = min ( n_large , first + nevents if 0 < nevents else n_large )
    

    if 0 <= first and 0 < nevents < chunk_size :
        return chain.statVar ( what , cuts , first , last )
    elif isinstance ( chain , ROOT.TChain ) : 
        if chain.nFiles() < 5 and len ( chain ) < chunk_size :
            return chain.statVar ( what , cuts , first , last )                         
    elif isinstance ( chain , ROOT.TTree  ) and len ( chain ) < chunk_size :
        return chain.statVar ( what , cuts , first , last ) 
    
    from ostap.trees.trees import Chain
    ch     = Chain ( chain , first = first , nevents = nevents )

    task   = StatVarTask ( what , cuts )
    wmgr   = WorkManager ( ppservers = ppservers , silent = silent )

    trees  = ch.split ( chunk_size = chunk_size , max_files = max_files )

    wmgr.process ( task , trees )

    del trees
    del ch    

    return task.output 

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
# The END 
# =============================================================================
