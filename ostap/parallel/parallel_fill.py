#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/parallel_fill.py
#  (parallel) Fill of RooDataSet frmo looong TChain/TTree
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
# =============================================================================
"""(parallel) Fill of RooDataSet frmo looong TChain/TTree
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'pprocess'  , ## paralell processing 
    ) 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.fill' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
import ROOT
from   ostap.parallel.parallel import Task, WorkManager
# =============================================================================
## The simple task object for more efficient fill of RooDataSet from TChain 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @see Ostap.SelectorWithVars
#  For 12-core machine, clear speed-up factor of about 8 is achieved 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23 
class  FillTask(Task) :
    """The single task object for more efficient fill of RooDataSet from TChain 
    - for 12-core machine, clear speed-up factor of about 8 is achieved 
    """
    ## 
    def __init__ ( self               ,
                   variables          ,
                   selection          ,
                   trivial   = False  , 
                   use_frame = 100000 ) :
        
        self.variables = variables 
        self.selection = selection 
        self.trivial   = trivial  
        self.use_frame = use_frame 
        self.__output  = ()  

    def initialize_local   ( self ) : self.__output = () 

    ## the actual processing 
    def process ( self , jobid , item ) :

        import ROOT
        from ostap.logger.logger import logWarning
        with logWarning() :
            import ostap.core.pyrouts            
            import ostap.trees.trees

        
        ## reconstruct chain from the item 
        chain    = item.chain
        ll       = len ( chain )  
        
        first    = item.first
        nevents  = item.nevents 

        all = 0 == first and ( nevents < 0 or ll <= nevents )
        
        if self.trivial and all : 
            import ostap.fitting.pyselectors
            self.__output = chain.make_dataset ( self.variables , self.selection , silent = True ) 
            return self.__output 

        from   ostap.fitting.pyselectors import SelectorWithVars
        
        ## use selector  
        selector = SelectorWithVars ( self.variables ,
                                      self.selection ,
                                      silence = True )
        
        args = ()  
        if not all : args  = nevents , first 
            
        num = chain.process ( selector , *args                 ,
                              shortcut  = all and self.trivial ,
                              use_frame = self.use_frame       )
        
        self.__output = selector.data, selector.stat  
        
        if  num < 0 :
            logger.warning ("Processing status %s (jobid #%s)"  %  ( num % jobid ) ) 
        
        ##del selector.data
        ##del      selector        
        logger.debug ( 'Processed %s and filled %d entries (jobid #%s)' % ( item , len( self.__output ) , jobid ) )

        return self.__output 

    ## merge results/datasets 
    def merge_results ( self, result ) :
        
        if result :
            ds , stat = result
            if not self.__output or not self.__output[0] :
                self.__output = ds , stat  
            else :
                ds_ , stat_ = self.__output
                ds_.append ( ds )
                stat_.total      += stat.total     ## total 
                stat_.processed  += stat.processed ## procesed 
                stat_.skipped    += stat.skipped   ## skipped                
                self.__output = ds_ , stat_
                ds.clear () 
            del result            
            logger.debug ( 'Merging: %d entries ' % len( self.__output[0] ) )
        else :
            logger.error ( "No valid results for merging" )
            
    ## get the results 
    def results ( self ) :
        return self.__output
    
# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain    = ...
#  selector =  ...
#  chain.pprocess ( selector )
#  @endcode 
def pprocess ( chain               ,
               selector            ,
               nevents    = -1     ,
               first      = 0      ,
               shortcut   = True   ,   ## important 
               chunk_size = 100000 ,   ## important 
               max_files  = 5      ,
               ppservers  = ()     ,
               use_frame  =  20000 ,   ## important 
               silent     = False  ) :
    """ Parallel processing of loooong chain/tree 
    >>>chain    = ...
    >>> selector =  ...
    >>> chain.pprocess ( selector )
    """
    
    from ostap.trees.trees import Chain

    ch = Chain ( chain ) 

    selection = selector.selection
    variables = selector.variables

    ## trivial   = selector.trivial_vars and not selector.morecuts
    
    trivial   = selector.really_trivial and not selector.morecuts 
    
    all = 0 == first and ( 0 > nevents or len ( chain ) <= nevents )
    
    if all and trivial and 1 < len( ch.files ) :
        logger.info ("Configuration is ``trivial'': redefine ``chunk-size'' to -1")
        chunk_size = -1
        
    task  = FillTask ( variables , selection , trivial , use_frame )
    wmgr  = WorkManager ( ppservers  = ppservers  , silent    = silent    )
    trees = ch.split    ( chunk_size = chunk_size , max_files = max_files )
    wmgr.process( task , trees )
    del trees
    
    dataset, stat = task.results()  

    selector.data = dataset
    selector.stat = stat 

    from ostap.logger.logger import attention 
    skipped = 'Skipped:%d' % stat.skipped
    skipped = '/' + attention ( skipped ) if stat.skipped else ''
    logger.info (
        'Selector(%s): Events Processed:%d/Total:%d%s CUTS: "%s"\n# %s' % (
        selector.name    ,
        stat.processed   ,
        stat.total       ,
        skipped          ,
        selector.cuts()  , dataset ) )            
    
    return 1 

ROOT.TChain.pprocess = pprocess
ROOT.TTree .pprocess = pprocess

# =============================================================================

_decorated_classes_ = (
    ROOT.TTree  ,
    ROOT.TChain ,    
    )

_new_methods_       = (
    ROOT.TTree .pprocess ,
    ROOT.TChain.pprocess ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
