#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some multiprocessing functionality for Ostap 
#  Currently it is not loaded on default, and requires manual activation
#  @see GaudiMP.Parallel
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
# =============================================================================
""" Multiprocessing functionality for Ostap
Currently it is not loaded on default, and requires manual activation
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'GenericTask' , ## "Generic Task" from three  basic functions 
    'ProjectTask' , ## "Project task" for very looooong chains/trees 
    'FillTask'    , ## "Fill task" for loooong chains/trees  
    'cproject'    , ##  project looong TChain into historgam   
    'tproject'    , ##  project looong TTree into histogram
    'fillDataSet' ,
    'WorkManager' 
    ) 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.kisa' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
logger.debug ( 'Multiprocessing functionality for Ostap')
# =============================================================================
import operator
import ostap.parallel.parallel as Parallel
WorkManager = Parallel.WorkManager 
# =============================================================================
n_large = 2**63 - 1  ## ROOT.TVirtualTreePlayer.kMaxEntries
## n_large = ROOT.TVirtualTreePlayer.kMaxEntries
# =============================================================================
## @class GenericTask
#  Generic task for Parallel processing  
#    One needs to  define three functions/functors:
#    - processor   :<code>        output = processor   ( item )               </code>
#    - merger      :<code>updated_output = merger ( old_output , new_output ) </code>
#    - initializer :<code>        output = initializer (      )               </code> 
class GenericTask(Parallel.Task) :
    """Generic task for parallel processing.
    One needs to  define three functions/functors:
    - processor   :         output = processor   ( item ) 
    - merger      : updated_output = merger ( old_output , new_output )
    - initializer :         output = initializer (      )  
    """
    # =========================================================================
    def __init__ ( self                       ,
                   processor                  ,
                   merger      = operator.add ,
                   initializer = tuple        ) :
        """Generic task for parallel processing. One needs to  define three functions/functors
        - processor   :         output = processor   ( item ) 
        - merger      : updated_output = merger ( old_output , new_output )
        - initializer :         output = initializer (      )  
        """        
        self.__processor   = processor
        self.__merger      = merger
        self.__initializer = initializer
        
    # =========================================================================
    ## local initialization (executed once in parent process)
    def initializeLocal   ( self ) :
        """Local initialization (executed once in parent process)"""
        self.output = self.initializer () if self.initializer else () 
        
    # =========================================================================
    ## the actual processing of the single item 
    def process  ( self , item ) :
        """The actual processing of the single item"""
        self.output = self.processor ( item )
        
    # =========================================================================
    ## merge results 
    def _mergeResults ( self , result ) :
        """Merge processing results"""
        self.output = self.merger ( self.output , result )

    # =========================================================================
    @property
    def processor  ( self ) :
        """``processor'' : the actual function for each subprocess
        - Signature: output = processor ( item ) 
        """
        return self.__processor
    @property
    def merger     ( self ) :
        """``merger'' : the actual fuction to merge results
        - Signature: updated_output = merger ( old_output , new_output )         
        """
        return self.__merger
    @property
    def initializer ( self ) :
        """``initializer'' : the actual fuction to initialize local output  
        - Signature: output = initializer() 
        """
        return self.__initializer
    
# =============================================================================
## The simple task object for more efficient projection of loooong chains/trees 
#  into histogarms
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
class ProjectTask(Parallel.Task) :
    """The simple task  object for the efficient parallel
    projection of looooooong TChains/TTrees into histograms  
    """
    ## constructor: histogram 
    def __init__ ( self , histo , what , cuts = '' ) :
        """Constructor: the histogram 
        
        >>> histo = ...
        >>> task  = ProjectTask ( histo ) 
        """        
        self.histo = histo
        self.what  = what 
        self.cuts  = str(cuts) 
        self.histo.Reset()
        
    ## local initialization (executed once in parent process)
    def initializeLocal   ( self ) :
        """Local initialization (executed once in parent process)
        """
        import ROOT,ostap.core.pyrouts
        self.output = 0, self.histo.clone()
    
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
        with logWarning() : import ostap.core.pyrouts 

        import ostap.trees.trees
        
        chain    = item.chain 
        first    = item.first
        nevents  = item.nevents

        ## Create the output histogram   NB! (why here???) 
        self.output = 0 , self.histo.Clone()
        
        ## use the regular projection  
        from ostap.trees.trees import _tt_project_ 
        self.output = _tt_project_ ( chain      , self.output[1] ,
                                     self.what  , self.cuts      ,
                                     ''         ,
                                     nevents    , first          )
        del item
        
    ## finalization (executed at the end at parent process)
    def finalize ( self ) : pass 

    ## merge results 
    def _mergeResults ( self , result ) :
        filtered    = self.output[0] + result[0] 
        self.output[1].Add ( result[1] )
        self.output = filtered, self.output[1]
        result[1].Delete () 
 
   
# =============================================================================  
## make a projection of the loooooooong chain into histogram using
#  multiprocessing functionality for per-file parallelisation
#  @code
#  >>> chain = ... ## large chain
#  >>> histo = ... ## histogram template 
#  >>> project        ( chain , histo , 'mass' , 'pt>10' )
#  >>> chain.pproject ( histo , 'mass' , 'pt>0' ) ## ditto 
#  >>> chain.cproject ( histo , 'mass' , 'pt>0' ) ## ditto 
#  @endcode
#  For 12-core machine, clear speedup factor of about 8 is achieved 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
def  cproject ( chain                ,
                histo                ,
                what                 ,
                cuts                 ,
                nentries   = -1      ,
                first      =  0      ,
                chunk_size = 1000000 ,
                silent     = False   ) :
    """Make a projection of the loooong chain into histogram
    >>> chain = ... ## large chain
    >>> histo = ... ## histogram template 
    >>> cproject        ( chain , histo , 'mass' , 'pt>10' )
    >>> chain.ppropject ( histo , 'mass' , 'pt>0' ) ## ditto 
    >>> chain.cpropject ( histo , 'mass' , 'pt>0' ) ## ditto     
    For 12-core machine, clear speedup factor of about 8 is achieved     
    """
    #
    
    from ostap.trees.trees import Chain
    ch    = Chain ( chain , first = first , nevents = nentries )
    
    task  = ProjectTask            ( histo , what , cuts )
    wmgr  = Parallel.WorkManager   ( silent = silent )
    wmgr.process ( task, ch.split  ( chunk_size  = chunk_size  ) )
    
    filtered   = task.output[0] 
    histo     += task.output[1]
    
    return filtered , histo 

import ROOT 
ROOT.TChain.cproject = cproject
ROOT.TChain.pproject = cproject

# =============================================================================  
## make a projection of the loooooooong tree into histogram using
#  multiprocessing functionality for per-file parallelisation
#  @code
#
#  >>> tree  = ... ## large tree 
#  >>> histo = ... ## histogram template 
#  >>> tproject ( tree , histo , 'mass' , 'pt>10' , maxentries = 1000000 )
#  >>> tree.pproject ( histo , 'mass' , 'pt>10' ) ## ditto 
#  @endcode
#  - significant gain can be achieved for very large ttrees with complicated expressions and cuts
#  - <code>maxentries</code> parameter should be rather large
#  @param tree       the tree
#  @param histo      the histogram
#  @param what       variable/expression/varlist to be projected
#  @param cuts       selection/weighting criteria 
#  @param nentries   number of entries to process  (>0: all entries in th tree)
#  @param first      the first entry to process
#  @param maxentries chunk size for parallel processing 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
def  tproject ( tree                 ,   ## the tree 
                histo                ,   ## histogram 
                what                 ,   ## variable/expression/list to be projected 
                cuts       = ''      ,   ## selection/weighting criteria 
                nentries   = -1      ,   ## number of entries 
                first      =  0      ,   ## the first entry 
                chunk_size = 1000000 ,   ## chunk size 
                silent     = False   ) : ## silent processing 
    """Make a projection of the loooong tree into histogram
    >>> tree  = ... ## large chain
    >>> histo = ... ## histogram template 
    >>> tproject ( tree , histo , 'mass' , 'pt>10' )    
    >>> tree.pproject ( histo , 'mass' , 'pt>10' )    ## ditto 
    - significant gain can be achieved for very large TTrees with complicated expressions and cuts
    - maxentries parameter should be rather large
    Arguments:
    - tree       the tree
    - histo      the histogram
    - what       variable/expression/varlist to be projected
    - cuts       selection/weighting criteria 
    - nentries   number of entries to process  (>0: all entries in th tree)
    - first      the first entry to process
    - maxentries chunk size for parallel processing 
    """

    from ostap.trees.trees import Tree
    ch    = Tree ( tree , first = first , nevents = nevents )
    
    task  = ProjectTask            ( histo , what , cuts )
    wmgr  = Parallel.WorkManager   ( silent     = silent       )
    wmgr.process ( task, ch.split  ( chunk_size = chunk_size ) )
    
    filtered   = task.output[0] 
    histo     += task.output[1]
    
    return filtered , histo 

import ROOT 
ROOT.TTree.tproject = tproject
ROOT.TTree.pproject = tproject

 
# =============================================================================
## The simple task object for more efficient fill of RooDataSet from TChain 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @see Ostap.SelectorWithVars
#  For 12-core machine, clear speed-up factor of about 8 is achieved 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23 
class  FillTask(Parallel.Task) :
    """The single task object for more efficient fill of RooDataSet from TChain 
    - for 12-core machine, clear speed-up factor of about 8 is achieved 
    """
    ## 
    def __init__ ( self ,  variables , selection , trivial = False ) :
        
        self.variables = variables 
        self.selection = selection 
        self.trivial   = trivial  
        self.output    = ()  

    def initializeLocal   ( self ) : pass
    def initializeRemote  ( self ) : pass

    ## the actual processing 
    def process ( self , item ) :

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
            import ostap.fitting.selectors
            self.output = chain.make_dataset ( self.variables , self.selection , silent = True ) 
            return

        from   ostap.fitting.selectors import SelectorWithVars
        
        ## use selector  
        selector = SelectorWithVars ( self.variables ,
                                      self.selection ,
                                      silence = True )
        
        args = ()  
        if not all : args  = nevents , first 
            
        num = chain.process ( selector , *args , shortcut = all and self.trivial )
        self.output = selector.data, selector.stat  
        
        if  num < 0 :
            logger.warning ("Processing status %s"  % num )
        
        ##del selector.data
        ##del      selector        
        logger.debug ( 'Processed %s and filled %d entries ' % ( item , len( self.output ) ) )

        del item

    def finalize ( self ) : pass 

    ## merge results/datasets 
    def _mergeResults(self, result) :
        
        if result :
            ds , stat = result
            if not self.output or not self.output[0] :
                self.output = ds , stat  
            else :
                ds_ , stat_ = self.output
                ds_.append ( ds )
                stat_     = list(stat_)
                stat_[0] += stat[0] ## total 
                stat_[1] += stat[1] ## procesed 
                stat_[2] += stat[2] ## skipped                
                self.output = ds_ , tuple(stat_)
                ds.clear () 
            del result            
            logger.debug ( 'Merging: %d entries ' % len( self.output[0] ) )
        else :
            logger.error ( "No valid results for merging" )

# =============================================================================
## The simple task object collect statistics for loooooong chains 
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
class StatVarTask(Parallel.Task) :
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
## @class ChopperTraining
#  parallel procession for TMVA chopping
#  @see ostap/tools/chopping.py
class ChopperTraining(Parallel.Task) :
    def __init__          ( self ) : self.output = ()
    def initializeLocal   ( self ) : self.output = () 
    def process           ( self , params ) :

        import ostap.tools.tmva
        category , chopper , log = params
        trainer  = chopper.create_trainer ( category , False )
        trainer.train ( log )
        self.output = (
            [ ( category , trainer.weights_files ) ] ,
            [ ( category , trainer.  class_files ) ] ,
            [ ( category , trainer. output_file  ) ] ,
            [ ( category , trainer.    tar_file  ) ] ,
            [ ( category , trainer.    log_file  ) ] ,
            )
        
    ## merge results/datasets 
    def _mergeResults( self , result) :
        if not  self.output : self.output =  result
        else :
            weights  = list ( self.output[0] ) + list ( result[0] ) 
            classes  = list ( self.output[1] ) + list ( result[1] ) 
            outputs  = list ( self.output[2] ) + list ( result[2] ) 
            tarfiles = list ( self.output[3] ) + list ( result[3] ) 
            logfiles = list ( self.output[4] ) + list ( result[4] ) 
            weights  . sort()
            classes  . sort()
            outputs  . sort()
            tarfiles . sort()
            logfiles . sort()
            self.output = weights , classes , outputs , tarfiles, logfiles  
                            
# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain    = ...
#  selector =  ...
#  chain.pprocess ( selector )
#  @endcode 
def _pprocess_ ( chain , selector , nevents = -1 , first = 0 , shortcut = True  , chunk_size = 100000 , ppservers = () , max_files = 10 , silent = False ) :
    """ Parallel processing of loooong chain/tree 
    >>>chain    = ...
    >>> selector =  ...
    >>> chain.pprocess ( selector )
    """
    
    from ostap.trees.trees import Chain

    ch = Chain ( chain ) 

    selection = selector.selection
    variables = selector.variables

    trivial   = selector.trivial_vars and not selector.morecuts 
    
    all = 0 == first and ( 0 > nevents or len ( chain ) <= nevents )
    
    if all and trivial and 1 < len( ch.files ) :
        logger.info ("Configuration is ``trivial'': redefine ``chunk-size'' to -1")
        chunk_size = -1
        
    task  = FillTask ( variables , selection , trivial )
    wmgr  = Parallel.WorkManager( ppservers = ppservers , silent = silent )
    trees = ch.split ( chunk_size = chunk_size , max_files = max_files )
    wmgr.process( task , trees )
    del trees
    
    dataset, stat = task.output 

    selector.data = dataset
    selector.stat = stat 

    from ostap.logger.logger import attention 
    skipped = 'Skipped:%d' % stat[2]
    skipped = '/' + attention ( skipped ) if stat[2] else ''
    logger.info (
        'Selector(%s): Events Processed:%d/Total:%d%s CUTS: "%s"\n# %s' % (
        selector.name    ,
        stat[1]          ,
        stat[0]          ,
        skipped          ,
        selector.cuts()  , dataset ) )            
    
    return 1 

    
ROOT.TChain.pprocess =  _pprocess_ 
ROOT.TTree.pprocess  =  _pprocess_ 




# ===================================================================================
## parallel processing of loooong chain/tree 
#  @code
#  chain    = ...
#  chain.pStatVar ( .... ) 
#  @endcode 
def _pStatVar_ ( chain        , what , cuts = ''    ,
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
    wmgr   = Parallel.WorkManager ( ppservers = ppservers , silent = silent )

    trees  = ch.split ( chunk_size = chunk_size , max_files = max_files )

    wmgr.process ( task , trees )

    del trees
    del ch    

    return task.output 

ROOT.TChain.pstatVar = _pStatVar_ 
ROOT.TTree .pstatVar = _pStatVar_ 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    if not ( 2**32 - 1 ) <= n_large <= ROOT.TVirtualTreePlayer.kMaxEntries :
        logger.error ( "Invalid setting of ``n_large''(%d) parameter (>%d)" % ( n_large , ROOT.TVirtualTreePlayer.kMaxEntries ) )
        
# =============================================================================
# The END 
# =============================================================================
