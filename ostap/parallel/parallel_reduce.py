#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/parallel_reduce.py
#  (parallel) Reduce long TChain 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-10-11
# =============================================================================
"""(parallel) Reduce long TChain 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2019-10-11"
__all__     = (
    'reduce'  , ## paralell reduce 
    ) 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.reduce' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
import ROOT
from   ostap.parallel.parallel import Task, WorkManager
# =============================================================================
## The simple task object for more efficient reduce of long TChain 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-10-11
class  ReduceTask(Task) :
    """The single task object for more efficient fill of RooDataSet from TChain 
    - for 12-core machine, clear speed-up factor of about 8 is achieved 
    """
    ## 
    def __init__ ( self               ,
                   selection  = {}    ,
                   save_vars  = ()    ,
                   new_vars   = {}    ,
                   addselvars = False ) :
        
        self.selection   = selection
        self.save_vars   = save_vars
        self.new_vars    = new_vars 
        self.addselvars  = addselvars 

        self.__output  = ()  

    def initialize_local   ( self ) : self.__output = () 

    ## the actual processing 
    def process ( self , item ) :

        import ROOT
        import ostap.core.pyrouts            
        from   ostap.trees.trees        import Chain 
        from   ostap.frames.tree_reduce import ReduceTree
        
        ## unpack the input data 
        chain = item.chain
        
        rt = ReduceTree ( chain                        ,
                          selection  = self.selection  ,
                          save_vars  = self.save_vars  ,
                          new_vars   = self.new_vars   ,
                          addselvars = self.addselvars ,
                          tmp_keep   = True            , ## attention! True is here! 
                          silent     = True            )
        
        cname = rt.chain.GetName()
        cfile = rt.output
        
        return Chain ( name = cname ,  files = [ cfile ] ) , rt.table 
    
    ## merge results/datasets 
    def merge_results ( self, result ) :
        
        if result :
            if not self.__output : self.__output = result 
            else                 :
                a , b = self.__output
                c , d = result
                m     = []
                if len ( b ) != d : logger.warning ( 'merge: mismatch in table size!!')
                for i , j in zip ( b , d ) :
                    n1, p1, a1 = i
                    n2, p2, a2 = j
                    if n1 != n2 : logger.warning ('merge: mismatch in row names!')
                    item = n1 , p1+p2 , a1 + a1                    
                    m.append ( item )
                self.__output = a + c , m 
                                                         
    ## get the results 
    def results ( self ) :
        return self.__output
    
# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain    = ...
#  chain.
#  @endcode 
def reduce ( chain               ,
             selection  = {}     ,
             save_vars  = ()     ,
             new_vars   = {}     ,
             no_vars    = ()     , 
             output     = ''     , 
             addselvars = False  ,
             silent     = False  ,
             ## 
             ppservers  = ()     ) :
    """ Parallel processing of loooong chain/tree 
    >>>chain    = ...
    >>> selector =  ...
    >>> chain.pprocess ( selector )
    """

    from ostap.trees.trees        import Chain
    from ostap.frames.tree_reduce import ReduceTree

    
    if isinstance ( chain , ROOT.TChain ) and 1 >= len ( chain.files() ) :
        return chain.reduce ( selection  = selection  ,
                              save_vars  = save_vars  ,
                              new_vars   = new_vars   ,
                              no_vars    = no_vars    ,
                              output     = output     ,
                              addselvars = addselvars ,
                              silent     = silent     )
    
    ch   = Chain ( chain ) 

    task = ReduceTask ( selection  ,
                        save_vars  ,
                        new_vars   ,
                        addselvars )
    
    wmgr  = WorkManager  ( ppservers = ppservers  , silent = silent    )
    trees = ch.split     ( max_files = 1  )
    wmgr.process         ( task , trees   )

    result , table  = task.results ()
    for i in result.files : result.trash.add ( i )

    if not silent :
        from ostap.frames.frames import report_print_table 
        title = 'Tree -> Frame -> Tree filter/transformation'
        logger.info ( 'Reduce tree:\n%s' % report_print_table ( table , title , '# ' ) )
        
    if not output : return result
    
    reduced = ReduceTree ( result.chain        ,
                           selection  = ''     ,
                           save_vars  = ()     ,
                           addselvars = False  ,
                           silent     = silent ,
                           output     = output )

    return Chain ( reduced.chain ) 

ROOT.TTree .preduce  = reduce
ROOT.TTree .ppreduce = reduce
ROOT.TTree .p_reduce = reduce

# =============================================================================

_decorated_classes_ = (
    ROOT.TTree  ,    
    )

_new_methods_       = (
    ROOT.TTree.preduce  ,
    ROOT.TTree.ppreduce ,
    ROOT.TTree.p_reduce ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
#                                                                       The END 
# =============================================================================
