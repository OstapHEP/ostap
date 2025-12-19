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
    'parallel_make_dataset' , ## efficient function for parallel fill of dataset 
    'parallel_fill_dataset' , ## flexible  function for parallel fill of dataset 
    'parallel_fill'         , ## flexible  function for parallel fill of dataset 
    ) 
# =============================================================================
from   ostap.parallel.parallel import Task, WorkManager
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.fill' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
## The simplest task object for efficient fill of RooDataSet from TChain 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  For 12-core machine, clear speed-up factor of about 8 is achieved 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23 
class  MakeDSTask(Task) :
    """ The simplest task object for efficient fill of RooDataSet from TChain 
    - for 12-core machine, clear speed-up factor of about 8 is achieved 
    """
    ## 
    def __init__ ( self           ,
                   variables      ,
                   selection      ,
                   roo_cuts  = '' ,
                   name      = '' ,
                   title     = '' ) :
        
        self.variables  = variables 
        self.selection  = selection
        self.roo_cuts   = roo_cuts
        self.name       = name
        self.title      = title         
        self.the_output = ()  

    ## local initialization 
    def initialize_local   ( self ) : self.the_output = () 

    ## the actual processing 
    def process ( self , jobid , item ) :

        import ROOT
        from ostap.logger.logger import logWarning
        with logWarning() :
            import ostap.core.pyrouts            
            import ostap.trees.trees
            import ostap.fitting.roofit 
            import ostap.fitting.pyselectors
        
        ## reconstruct chain from the item 
        chain    = item.chain        
        return chain.make_dataset ( variables = self.variables ,
                                    selection = self.selection ,
                                    roo_cuts  = self.roo_cuts  ,
                                    name      = self.name      ,
                                    title     = self.title     ,                                             
                                    silent    = True           ) 
    
    ## merge results/datasets 
    def merge_results ( self , result , jobid = -1 ) :
        """ Merge results/datasets
        """
        from ostap.fitting.dataset import Ostap
        if result :
            ds , stat = result
            if not self.the_output or not self.the_output[0] :
                self.the_output = ds , stat  
            else :
                ds_ , stat_ = self.the_output
                ds_.append ( ds )
                stat_.total      += stat.total     ## total 
                stat_.processed  += stat.processed ## procesed 
                stat_.skipped    += stat.skipped   ## skipped                
                self.the_output   = ds_ , stat_
                if isinstance ( ds , ROOT.RooDataSet ) :
                    ds = Ostap.MoreRooFit.delete_data ( ds )
                    del ds 
            del result            
            logger.debug ( 'Merging: %d entries ' % len( self.the_output[0] ) )
        else :
            logger.error ( "No valid results for merging" )
            
    ## get the results 
    def results ( self ) :
        return self.the_output
    
# =================================================================================
## The simple task object for more efficient fill of RooDataSet from TChain 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @see Ostap.SelectorWithVars
#  For 12-core machine, clear speed-up factor of about 8 is achieved 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23 
class  FillDSTask(MakeDSTask) :
    """The single task object for more efficient fill of RooDataSet from TChain 
    - for 12-core machine, clear speed-up factor of about 8 is achieved 
    """
    ## 
    def __init__ ( self               ,
                   variables          ,
                   selection          ,
                   cuts      = None   ,
                   roo_cuts  = ''     ,
                   name      = ''     ,
                   title     = ''     , 
                   shortcut  = True   , 
                   use_frame = 100000 ) :
        
        MakeDSTask.__init__ ( self                  ,
                             variables = variables , 
                             selection = selection ,
                             roo_cuts  = roo_cuts  ,
                             name      = name      ,
                             title     = title     )
        
        self.cuts      = cuts     
        self.use_frame = use_frame 
        self.shortcut  = shortcut
        
    ## the actual processing 
    def process ( self , jobid , item ) :

        import ROOT
        from ostap.logger.logger import logWarning
        with logWarning() :
            import ostap.core.pyrouts            
            import ostap.trees.trees
            import ostap.fitting.roofit 
            from   ostap.fitting.pyselectors import SelectorWithVars
            
        ## reconstruct chain from the item 
        chain    = item.chain
                
        ## use selector  
        selector = SelectorWithVars ( variables = self.variables ,
                                      selection = self.selection ,
                                      cuts      = self.cuts      ,
                                      roo_cuts  = self.roo_cuts  ,
                                      name      = self.name      ,
                                      fuillname = self.title     , 
                                      silence   = True           )
        
        return chain.fill_dataset2 ( selector  ,
                                     silent    = True            , 
                                     shortcut  = self.shortcut   ,
                                     use_frame = self.use_frame  )
        
# =================================================================================
## The simple task object for more efficient fill of RooDataSet from TChain 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @see Ostap.SelectorWithVars
#  For 12-core machine, clear speed-up factor of about 8 is achieved 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23 
class  FillTask(MakeDSTask) :
    """The single task object for more efficient fill of RooDataSet from TChain 
    - for 12-core machine, clear speed-up factor of about 8 is achieved 
    """
    ## 
    def __init__ ( self               ,
                   variables          ,
                   selection          ,
                   cuts      = None   ,
                   roo_cuts  = ''     ,
                   name      = ''     ,
                   title     = ''     , 
                   trivial   = False  , 
                   use_frame = 100000 ) :

        MakeDSTask.__init__ ( self                  ,
                              variables = variables , 
                              selection = selection ,
                              roo_cuts  = roo_cuts  ,
                              name      = name      ,
                              title     = title     )
        
        self.cuts      = cuts     
        self.trivial   = trivial  
        self.use_frame = use_frame 

    ## the actual processing 
    def process ( self , jobid , item ) :

        import ROOT
        from ostap.logger.logger import logWarning
        with logWarning() :
            import ostap.core.pyrouts            
            import ostap.trees.trees
            import ostap.fitting.roofit 
        
        ## reconstruct chain from the item 
        chain    = item.chain
        ll       = len ( chain )  
        
        first    = item.first
        nevents  = item.nevents 

        all = 0 == first and ( nevents < 0 or ll <= nevents )
        
        if self.trivial and all and not self.cuts : 
            import ostap.fitting.pyselectors
            return chain.make_dataset ( self.variables , self.selection , silent = True ) 
        
        from   ostap.fitting.pyselectors import SelectorWithVars
        
        ## use selector  
        selector = SelectorWithVars ( variables = self.variables ,
                                      selection = self.selection ,
                                      roo_cuts  = self.roo_cuts  ,
                                      cuts      = self.cuts      , 
                                      silence   = True           )
        
        args = ()  
        if not all : args  = nevents , first 
        
        return chain.fill_dataset2 ( selector  ,
                                     *args     ,
                                     silent    = True                 , 
                                     shortcut  = all and self.trivial ,
                                     use_frame = self.use_frame       )



# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain    = ...
#  selector =  ...
#  parallel_fill        ( chain , selector ) 
#  chain.parallel_fill  ( selector         ) ## ditto 
#  @endcode 
def parallel_fill ( chain                  ,
                    selector               ,
                    nevents      = -1      ,
                    first        = 0       ,
                    shortcut     = True    ,   ## important 
                    chunk_size   = 1000000 ,   ## important 
                    max_files    = 5       ,
                    use_frame    =  20000  ,   ## important 
                    silent       = False   ,
                    job_chunk    = -1      ,
                    progress     = True    , **kwargs ) :
    """ Parallel processing of loooong chain/tree 
    >>>chain    = ...
    >>> selector =  ...
    >>> chain.parallel_fill ( selector )
    """
    import ostap.fitting.roofit 
    from   ostap.fitting.pyselectors import SelectorWithVars 
    from   ostap.trees.trees         import Chain
    
    assert isinstance ( selector , SelectorWithVars ) , \
           "Invalid type of ``selector'': %s" % type ( selector ) 
    
    ch = Chain ( chain ) 

    selection = selector.selection
    variables = selector.variables
    roo_cuts  = selector.roo_cuts
    
    ## trivial   = selector.trivial_vars and not selector.morecuts
    
    trivial   = selector.really_trivial and not selector.morecuts 
    
    all = ( 0 == first ) and ( 0 > nevents or len ( chain ) <= nevents )
    
    if all and trivial and 1 < len ( ch.files ) :
        logger.info ("Configuration is `trivial': redefine `chunk-size' to -1")
        chunk_size = -1
        
    task  = FillTask ( variables = variables         ,
                       selection = selection         ,
                       cuts      = selector.morecuts , 
                       roo_cuts  = roo_cuts          ,
                       trivial   = trivial           ,
                       use_frame = use_frame         ) 
    
    wmgr  = WorkManager ( silent     = silent  , progress = True , **kwargs )
    trees = ch.split    ( chunk_size = chunk_size , max_files = max_files )
    wmgr.process ( task , trees , chunk_size = job_chunk )
    del trees
    
    dataset, stat = task.results()  

    selector.data = dataset
    selector.stat = stat 

    from ostap.logger.logger import attention 
    skipped = 'Skipped:%d' % stat.skipped
    skipped = '/' + attention ( skipped ) if stat.skipped else ''
    logger.info (
        'Selector(%s): Events Processed:%d/Total:%d%s CUTS: "%s"\n%s' % (
        selector.name    ,
        stat.processed   ,
        stat.total       ,
        skipped          ,
        selector.cuts()  , dataset.table ( prefix = '# ' ) ) )             
    
    return dataset, stat  

ROOT.TChain.parallel_fill = parallel_fill 
ROOT.TTree .parallel_fill = parallel_fill 


# =============================================================================
## Create RooDataset from the chain/tree
#  @code 
#  tree = ...
#  ds, stat = tree.pfill_dataset ( [ 'px , 'py' , 'pz' ] ) 
#  @endcode
def parallel_fill_dataset  ( chain                  ,
                             variables              ,  ## list of variables 
                             selection    = ''      ,  ## TTree-cuts 
                             roo_cuts     = ''      ,  ## RooFit cuts 
                             cuts         = None    ,  ## python callable 
                             name         = ''      ,
                             title        = ''      ,
                             shortcut     = True    ,
                             use_frame    = 20000   ,   ## important 
                             silent       = True    , 
                             max_files    = 1       ,
                             job_chunk    = -1      , **kwargs ) :
    """ Create RooDataset from the chain/tree
    >>> tree = ...
    >>> ds , stat = tree.pfill_dataset ( [ 'px , 'py' , 'pz' ] )
    """
    
    import ostap.fitting.roofit 
    from   ostap.fitting.pyselectors import SelectorWithVars 
    from   ostap.trees.trees         import Chain
    
    ch = Chain ( chain ) 
    if ch.nFiles < 2 :
        from ostap.fitting.pyselectors import fill_dataset1 as _fill_dataset_
        return _fill_dataset_ ( tree                  ,
                                variables = variables , 
                                selection = selection ,                                
                                roo_cut   = roo_cuts  ,
                                cuts      = cuts      ,
                                name      = name      ,
                                title     = title     ,
                                shortcut  = shortcut  ,
                                use_frame = use_frame ,
                                silent    = silent    )
    
    task  = FillDSTask ( variables = variables  ,
                         selection = selection  ,
                         cuts      = cuts       , 
                         roo_cuts  = roo_cuts   ,
                         name      = name       ,
                         title     = title      ,
                         shortcut  = shortcut   ,
                         use_frame = use_frame  )

    
    wmgr  = WorkManager ( silent     = silent     , **kwargs )
    trees = ch.split    ( chunk_size = -1 , max_files = max_files )
    wmgr.process ( task , trees , chunk_size = job_chunk )
    del trees
    
    dataset, stat = task.results()  

    from ostap.logger.logger import attention 
    skipped = 'Skipped:%d' % stat.skipped
    skipped = '/' + attention ( skipped ) if stat.skipped else ''
    logger.info (
        'Selector(%s): Events Processed:%d/Total:%d%s CUTS: "%s"\n%s' % (
        selector.name    ,
        stat.processed   ,
        stat.total       ,
        skipped          ,
        selector.cuts()  , dataset.table ( prefix = '# ' ) ) )             
    
    return dataset, stat  


ROOT.TChain.pfill_dataset = parallel_fill_dataset 
ROOT.TTree .pfill_dataset = parallel_fill_dataset

# =============================================================================
## Create RooDataset from the tree using parallel Tree->Frame->Dataset transformation
#  @code 
#  tree = ...
#  ds   = tree.pmake_dataset ( [ 'px , 'py' , 'pz' ] ) 
#  @endcode
def parallel_make_dataset ( chain                ,
                            variables            , ## variables 
                            selection  = ''      , ## TTree selection 
                            roo_cuts   = ''      , ## Roo-Fit selection  
                            name       = ''      , 
                            title      = ''      ,
                            silent     =  False  , 
                            max_files  =  1      ,
                            job_chunk  = -1      ,
                            progress   = True    , **kwargs ) :
    """ Create RooDataset from the tree using parallel Tree->Frame->Dataset transformation 
    >>> tree = ...
    >>> ds   = tree.pmake_dataset ( [ 'px , 'py' , 'pz' ] ) 
    """
    
    import ostap.fitting.roofit 
    from   ostap.trees.trees         import Chain
    
    ch = Chain ( chain )
    
    if ch.nFiles < 2  :
        from ostap.fitting.pyselectors import make_dataset as _make_dataset_
        return _make_dataset_ ( tree                  ,
                                variables = variables , 
                                selection = selection ,
                                roo_cut   = roo_cuts  ,
                                name      = name      ,
                                title     = title     ,
                                silent    = silent    ) 
                                
    # and here we have parallel processing 

    task = MakeDSTask ( variables = variables ,
                        selection = selection , 
                        roo_cuts  = roo_cuts  ,
                        name      = name      ,
                        title     = title     )
    
    
    wmgr  = WorkManager ( silent    = silent , progress = progress , **kwargs )
    trees = ch.split    ( chunk_size = -1 , max_files = max_files )
    wmgr.process ( task , trees , chunk_size = job_chunk )
    del trees

    dataset, stat = task.results()  

    from ostap.logger.logger import attention 
    skipped = 'Skipped:%d' % stat.skipped
    skipped = '/' + attention ( skipped ) if stat.skipped else ''
    logger.info (
        'Events Processed:%d/Total:%d%s CUTS: "%s"\n%s' % (
        stat.processed   ,
        stat.total       ,
        skipped          ,
        selection        , dataset.table ( prefix = '# ' ) ) )             
    
    return dataset, stat  

ROOT.TChain.pmake_dataset = parallel_make_dataset 
ROOT.TTree .pmake_dataset = parallel_make_dataset

# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain    = ...
#  selector =  ...
#  pprocess       ( chain , selector ) 
#  chain.pprocess ( selector         ) ## ditto 
#  @endcode 
def pprocess ( chain , selector , *args , **kwargs ) :
    """parallel processing of loooong chain/tree 
    >>> chain    = ...
    >>> selector =  ...
    >>> pprocess       ( chain , selector ) 
    >>> chain.pprocess ( selector         ) ## ditto 
    """
    
    from ostap.fitting.pyselectors import SelectorWithVars
    
    assert isinstance ( selector , SelectorWithVars ) , \
           "Invalid type of ``selector'': %s" % type ( selector ) 

    logger.warning ( "Use ``parallel_fill'' instead of ``pprocess''!" ) 
    
    return parallel_fill ( chain , selector, *args , **kwargs )  
    

ROOT.TChain.pprocess = pprocess 
ROOT.TTree .pprocess = pprocess 

# =============================================================================

_decorated_classes_ = (
    ROOT.TTree  ,
    ROOT.TChain ,    
    )

_new_methods_       = (
    ROOT.TTree .parallel_fill ,
    ROOT.TChain.parallel_fill ,
    ROOT.TTree .pfill_dataset ,
    ROOT.TChain.pfill_dataset ,
    ROOT.TTree .pmake_dataset ,
    ROOT.TChain.pmake_dataset ,
    ROOT.TTree .pprocess      ,
    ROOT.TChain.pprocess      ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
