#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/parallel_project.py
#  Paralllel "project" from loooong chain/trees objects   
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
# =============================================================================
"""Paralllel ``project'' from loooong chain/trees objects   
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'cproject' , ## parallel project from looong TChain
    'tproject' , ## parallel project from looong TTree 
    ) 
# =============================================================================
from   ostap.parallel.parallel import Task, WorkManager
import ostap.core.pyrouts 
import ostap.trees.trees
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.project'    )
else                       : logger = getLogger ( __name__     )
# =============================================================================
## The simple task object for more efficient projection of loooong chains/trees 
#  into histogarms
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
class ProjectTask(Task) :
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
        self.cuts  = str  ( cuts ) 
        self.histo.Reset()
        
    ## local initialization (executed once in parent process)
    def initialize_local   ( self ) :
        """Local initialization (executed once in parent process)
        """
        import ROOT,ostap.core.pyrouts
        self.__output = 0, self.histo.clone()
    
    ## remote initialization (executed for each sub-processs)
    def initialize_remote  ( self , jobid = -1 ) :
        """Remote initialization (executed for each sub-processs
        """
        import ROOT,ostap.core.pyrouts
        self.__output = 0, self.histo.clone()
    
    ## finalization (executed at the end at parent process)
    def finalize ( self ) : pass 

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
            import ostap.histos.histos
            import ostap.frames.frames
            from ostap.trees.trees import Chain,   Tree
            
        input    = Chain ( name    = item.name    ,
                           files   = item.files   , 
                           first   = item.first   ,
                           nevents = item.nevents )
        
        chain    = input.chain
        first    = input.first
        nevents  = input.nevents
        
        ## use the regular projection  
        from ostap.trees.trees import _tt_project_
        
        ## Create the output histogram  NB! (why here???)
        from ostap.core.core import ROOTCWD

        with ROOTCWD() :
            
            ROOT.gROOT.cd() 
            histo    =  self.histo.Clone ()
            self.__output = 0 , histo

            from ostap.trees.trees import _tt_project_            
            self.__output = _tt_project_ ( tree     = chain      , histo = histo      ,
                                           what     = self.what  , cuts  = self.cuts  ,
                                           options  = ''         ,
                                           nentries = nevents    , firstentry = first )            
        return self.__output 
        
    ## merge results 
    def merge_results ( self , result , jobid ) :

        import ostap.histos.histos
        if not self.__output : self.__output =  result
        else : 
            filtered      = self.__output[0] + result[0]            
            self.__output[1].Add ( result[1] )
            self.__output = filtered, self.__output[1]

        ## result[1].Delete () 
            
    ## get the results 
    def results (  self ) :
        return self.__output 
    
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
                chunk_size = -1      ,
                max_files  =  5      , 
                silent     = False   , **kwargs ) :
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
    
    task  = ProjectTask ( histo , what , cuts )
    wmgr  = WorkManager ( silent = silent , **kwargs )    
    wmgr.process ( task , ch.split ( chunk_size = chunk_size , max_files = max_files ) )

    ## unpack results 
    _f , _h    = task.results ()
    filtered   = _f
    histo     += _h
    del _h 
    
    return filtered , histo 

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
                max_files  = 50      ,   ## not-used .... 
                silent     = False   , **kwargs ) : ## silent processing 
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
    ch    = Tree ( tree , first = first , nevents = nentries )
    
    task  = ProjectTask            ( histo , what , cuts )
    wmgr  = WorkManager            ( silent     = silent , **kwargs )
    wmgr.process ( task, ch.split  ( chunk_size = chunk_size ) )
    
    ## unpack results 
    _f , _h = task.results ()
    filtered   = _f
    histo     += _h
    del _h 
    
    return filtered , histo 

ROOT.TTree.tproject = tproject
ROOT.TTree.pproject = tproject

# =============================================================================
_decorated_classes_ = (
    ROOT.TTree  ,
    ROOT.TChain ,    
    )

_new_methods_       = (
    ROOT.TTree .tproject ,
    ROOT.TTree .pproject ,
    ROOT.TChain.cproject ,
    ROOT.TChain.pproject ,     
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
