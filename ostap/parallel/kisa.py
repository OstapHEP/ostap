#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some multiprocessing functionality for Ostap 
#  Currently it is not loaded on default, and requires manual activation
#
#  @see GaudiMP.Parallel
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
#
# =============================================================================
""" Multiprocessing functionality for Ostap
Currently it is not loaded on default, and requires manual activation
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'ProjectTask' , ## "Project task" for very looooong chains/trees 
    'FillTask'    , ## "Fill task" for loooong chains/trees  
    'cproject'    , ##  project looong TChain into historgam   
    'tproject'    , ##  project looong TTree into histogram
    'fillDataSet' 
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
import ostap.parallel.parallel as Parallel  
# =============================================================================
n_large = 2**63 - 1  ## ROOT.TVirtualTreePlayer.kMaxEntries
## n_large = ROOT.TVirtualTreePlayer.kMaxEntries
# =============================================================================
## The simple task object for more efficient projection of loooong chains/trees 
#  into histogarms
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
class ProjectTask(Parallel.Task) :
    """ The simple task  object for efficient parallel
    projection of looooooong TChains/TTrees into histograms  
    """
    ## constructor: histogram 
    def __init__ ( self , histo ) :
        """Constructor: the histogram 
        
        >>> histo = ...
        >>> task  = ProjectTask ( histo ) 
        """        
        self.histo = histo
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
    def process ( self , params ) :
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

        if   isinstance ( params , str ) : params = ( param , 0 , n_large  )
        elif isinstance ( params , ROOT.TChainElement ) :
            params = ( params.GetTitle()  , 0 , n_large  )

        fname    = params[0] ## file name 
        tname    = params[1] ## tree name 
        what     = params[2] ## variable/expression to project 
        cuts     = params[3] if 3 < len ( params ) else ''       ## cuts    
        first    = params[4] if 4 < len ( params ) else 0        ## the first event
        nentries = params[5] if 5 < len ( params ) else n_large  ## number of events 
        
        if isinstance ( fname , ROOT.TChainElement ) : fname = fname.GetTitle() 
        
        chain = ROOT.TChain ( tname )
        chain.Add ( fname )
        
        ## Create the output histogram   NB! (why here???) 
        self.output = 0 , self.histo.Clone()
        
        ## use the regular projection  
        from ostap.trees.trees import _tt_project_ 
        self.output = _tt_project_ ( chain      , self.output[1] ,
                                     what       , cuts           ,
                                     ''         ,
                                     nentries   , first          )
        del chain 

    ## finalization (executed at the end at parent process)
    def finalize ( self ) : pass 

    ## merge results 
    def _mergeResults(self, result) :
        #
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
def  cproject ( chain , histo , what , cuts , silent = False ) :
    """Make a projection of the loooong chain into histogram
    >>> chain = ... ## large chain
    >>> histo = ... ## histogram template 
    >>> cproject        ( chain , histo , 'mass' , 'pt>10' )
    >>> chain.ppropject ( histo , 'mass' , 'pt>0' ) ## ditto 
    >>> chain.cpropject ( histo , 'mass' , 'pt>0' ) ## ditto     
    For 12-core machine, clear speedup factor of about 8 is achieved     
    """
    #
    if not chain :
        return 0 , histo
    if not histo :
        logger.error ('cproject: invalid histogram')
        return 0 , histo
    
    import ROOT
    histo.Reset()    

    if not isinstance ( chain , ROOT.TChain ) :
        logger.warning ('cproject method is TChain-specific, skip parallelization') 
        from ostap.trees.trees import _tt_project_
        return _tt_project_ ( chain , histo , what , cuts ) 

    if isinstance ( cuts , ROOT.TCut ) : cuts = str( cuts )
    ##
    if isinstance ( what  , str ) : what = what.split(',')
    if isinstance ( what  , str ) : what = what.split(';')
    if isinstance ( what  , str ) : what = [ what ] 
    
    import ostap.trees.trees
    files = chain.files()

    cname = chain.GetName() 
    
    params = [ ( f , cname , str(w) , cuts ) for f in files for w in what ] 

    task  = ProjectTask          ( histo )
    wmgr  = Parallel.WorkManager ( silent = silent )
    wmgr.process( task, params )

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
                maxentries = 1000000 ,   ## chunk size 
                silent     = False   ) : ## silent processing 
    """Make a projection of the loooong tree into histogram
    >>> tree  = ... ## large chain
    >>> histo = ... ## histogram template 
    >>> tproject ( tree , histo , 'mass' , 'pt>10' )    
    >>> tree.pproject ( histo , 'mass' , 'pt>10' )    ## ditto 
    - significant gain can be achieved for very large ttrees with complicated expressions and cuts
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
    if not tree  :
        return 0 , histo
    if not histo :
        logger.error ('tproject: invalid histogram')
        return 0 , histo
    
    import ROOT 
    histo.Reset()
    
    num = len( tree )
    if num <= first :
        return 0, histo

    if 0 >  nentries   : nentries   = n_large 

    maxentries = long(maxentries)
    if 0 >= maxentries : maxentries = n_large 
    
    if 0 > first       : first      = 0
    
    ## use the regular projection  
    from ostap.trees.trees import _tt_project_ 

    fname = None
    tname = None
    
    if isinstance ( tree , ROOT.TChain ) :
        
        if 1 == len( tree.files() ) :
            
            fname = tree.files()[0]
            tname = tree.GetName()
            
        else :
            
            logger.warning ('``tproject'' method is TTree-specific, skip parallelization')
            return _tt_project_ ( tree , histo , what , cuts , '' , nentries , first )
        
    else :         

        tdir  = tree.GetDirectory ()
        ftree = tdir.GetFile      ()
        if not ftree :
            logger.debug ('TTree is not file resident, skip parallelization') 
            return _tt_project_ ( tree ,  histo , what , cuts , '', total , first  )
        fname         = ftree.GetName     ()
        tpath         = tdir.GetPath  ()
        pr , d , path = tpath.rpartition(':')
        tname         = path + '/' + tree.GetName()

    if not fname :
        logger.info ("Can't determine fname, skip parallelization") 
        return _tt_project_ ( tree ,  histo , what , cuts , '', total , first  )
    
    if not tname :
        logger.info ("Can't determine tname, skip parallelization") 
        return _tt_project_ ( tree ,  histo , what , cuts , '', total , first  )
        
    # 
    if isinstance ( cuts , ROOT.TCut ) : cuts = str( cuts )
    if isinstance ( what , ROOT.TCut ) : what = str( what )
    ##
    if isinstance ( what , str ) : what = what.split(',')
    if isinstance ( what , str ) : what = what.split(',')
    if isinstance ( what , str ) : what = what.split(';')
    if isinstance ( what , str ) : what = [ what ] 

    ## nothing to project 
    if not what :
        return 0 , histo
        
    ## total number of events to process :
    total = min ( num - first , nentries ) 

    ## the event range is rather short, no real need  in parallel processing
    if total * len ( what ) < maxentries and len ( what ) < 4 : 
        return _tt_project_ ( tree ,  histo , what , cuts , '', total , first  )


    ## number of chunks & reminder 
    nchunks , rest = divmod ( total , maxentries )
    csize          = int ( total / nchunks ) ## chunk size 

    ## final list of parameters [ (file_name, what , cuts , first_event , num_events ) , ... ] 
    params = []

    for i in range(nchunks) :
        for w in what : 
            params.append ( ( fname , tname , str(w) , cuts , first +       i * csize , csize ) )
            
    if rest :
        nchunks +=  1
        for w in what : 
            params.append ( ( fname , tname , str(w) , cuts , first + nchunks * csize , rest  ) )

    task  = ProjectTask          ( histo )
    wmgr  = Parallel.WorkManager ( silent = silent )
    wmgr.process( task, params )

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
    def __init__ ( self ,  variables , selection ) :
        
        self.variables = variables
        self.selection = selection 
        
        self.output   = ()  

    def initializeLocal   ( self ) : pass
    def initializeRemote  ( self ) : pass

    ## the actual processing 
    def process ( self , params ) :

        import ROOT
                    
        tree,fname = params
        tree = ROOT.TChain( tree )
        tree.Add ( fname )

        from ostap.logger.logger import logWarning
        with logWarning() : 
            
            import ostap.core.pyrouts
            from   ostap.fitting.selectors import SelectorWithVars
            
        selector = SelectorWithVars ( self.variables ,
                                      self.selection ,
                                      silence = True )
        
        tree,fname = params
        tree = ROOT.TChain( tree )
        tree.Add ( fname )

        num =  tree.process ( selector )
        self.output =  selector.data, num  

        if  num < 0 :
            logger.warning ("Processing status %s"  % num )
        
        ##del selector.data
        ##del      selector
        
        logger.debug ( 'Processed %s chain and filled %d entries ' % ( fname , len( self.output ) ) ) 

    
    def finalize ( self ) : pass 

    ## methge resulsts/datasets 
    def _mergeResults(self, result) :
        #
        if result :
            ds , st = result            
            if not self.output or not self.output[0] :
                self.output = ds , st 
            else :
                ds_ , st_ = self.output
                ds_.append ( ds )
                st_ += st
                self.output = ds_ , st_
                ds.clear () 
            del result            
            logger.debug ( 'Merging: %d entries ' % len( self.output[0] ) )
        else :
            logger.error("No valid results for merging")

# ==============================================================================
## Fill dataset from looooong TChain using per-file parallelisation
#  @code
#  >>> chain =
#  >>> vars  = ...
#  >>> dset  = fillDataSet ( chain , vars , 'pt>10' )
#  @endcode
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @see Ostap.SelectorWithVars
#  For 12-core machine, clear speed-up factor of about 8 is achieved 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23 
def  fillDataSet ( chain , variables , selection , ppservers = () , silent = False ) :
    """Fill dataset from loooong TChain using per-file parallelisation
    >>> chain =
    >>> vars  = ...
    >>> dset  = fillDataSet ( chain , vars , 'pt>10' )
    - for 12-core machine, clear speed-up factor of about 8 is achieved 
    """

    task  = FillTask ( variables , selection )
    wmgr  = Parallel.WorkManager( ppservers = ppservers , silent = silent )
    
    cname = chain.GetName() 
    files = chain.files() 
    pairs = [ ( cname,i ) for i in files ] 

    wmgr.process( task, pairs )

    return task.output


# ===================================================================================
## parallel processing of loooong chain
#  @code
#  chain    = ...
#  selector =  ...
#  chain.pprocess ( selector )
#  @endcode 
def _pprocess_ ( chain , selector , ppservers = () , silent = False ) :
    
    if not isinstance ( chain , ROOT.TChain ) :
        logger.warning ('pprocess method is TChain-specific, skip parallelization') 
        return chain.process ( selector )  
    
    selection = selector.selection
    variables = selector.variables

    dataset, status = fillDataSet ( chain , variables ,  selection , ppservers , silent )
    selector.data = dataset
    
    return status 

    
ROOT.TChain.pprocess =  _pprocess_ 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    if not ( 2**32 - 1 ) <= n_large <= ROOT.TVirtualTreePlayer.kMaxEntries :
        logger.error ( "Invalid setting of ``n_large''(%d) parameter (>%d)" % ( n_large , ROOT.TVirtualTreePlayer.kMaxEntries ) )
        
# =============================================================================
# The END 
# =============================================================================
