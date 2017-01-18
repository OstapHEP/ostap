#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
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
    ## constructor: chain name, historgam , variable , cuts 
    def __init__ ( self , tree , histo , what , cuts = '' ) :
        """Constructor: chain/tree name, historgam , variable , cuts
        
        >>> histo = ...
        >>> task  = ProjectTask ( 'MyTuple' , histo , 'mass' , 'pt>10' ) 
        """        
        self.histo = histo
        self.histo.Reset()
        
        import ROOT
        if   isinstance ( tree , ROOT.TTree ) : self.tree = tree.GetName()
        elif isinstance ( tree , str        ) : self.tree = tree 
        
        self.what = what
        self.cuts = cuts

    ## local initialization (executed once in parent process)
    def initializeLocal   ( self ) :
        """Local initialization (executed once in parent process)
        """
        import ROOT,Ostap.PyRoUts
        self.output = 0, self.histo.Clone() 
        
    ## remote initialization (executed for each sub-processs)
    def initializeRemote  ( self ) : pass 
    
    ## the actual processing
    #   ``params'' is assumed to be a tuple :
    #  - the file name
    #  - the first entry in tree to process
    #  - number of entries to process
    def process ( self , params ) :
        """The actual processing
        ``params'' is assumed to be a tuple-like entity:
        0 - the file name
        1 - the first entry in tree to process 
        2 - number of entries to process        
        """

        import ROOT
        import Ostap.PyRoUts 

        if   isinstance ( params , str ) : params = ( param , 0 , n_large  )
        elif isinstance ( params , ROOT.TChainElement ) :
            params = ( params.GetTitle()  , 0 , n_large  )

        fname    = params[0]
        first    = params[1] if 1 < len(params) else 0 
        nentries = params[2] if 2 < len(params) else n_large 
        
        if isinstance ( fname , ROOT.TChainElement ) : fname = fname.GetTitle() 
        
        chain = ROOT.TChain ( self.tree )
        chain.Add ( fname )
        
        ## Create the output histogram   NB! (why here???) 
        self.output = 0 , self.histo.Clone()
        
        ## use the regular projection  
        from Ostap.TreeDeco import _tt_project_ 
        self.output = _tt_project_ ( chain      , self.output[1] ,
                                     self.what  , self.cuts      ,
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
def  cproject ( chain , histo , what , cuts ) :
    """Make a projection of the loooong chain into histogram
    >>> chain = ... ## large chain
    >>> histo = ... ## histogram template 
    >>> cproject        ( chain , histo , 'mass' , 'pt>10' )
    >>> chain.ppropject ( histo , 'mass' , 'pt>0' ) ## ditto 
    >>> chain.cpropject ( histo , 'mass' , 'pt>0' ) ## ditto     
    For 12-core machine, clear speedup factor of about 8 is achieved     
    """
    #
    
    if not tree  :
        return 0 , histo
    if not histo :
        logger.error ('cproject: invalid histogram')
        return 0 , histo
    
    import ROOT 
    histo.Reset()    

    if not isinstance ( ROOT , TChain ) :
        logger.warning ('cproject method is TChain-specific, skip parallelization') 
        from Ostap.TreeDeco import _tt_project_
        return _tt_project_ ( chain , histo , what , cuts ) 

    import Ostap.TreeDeco 
    files = chain.files()
    
    task  = ProjectTask          ( chain , histo , what , cuts )
    wmgr  = Parallel.WorkManager ()
    wmgr.process( task, files )

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
#  >>> tree  = ... ## large tree 
#  >>> histo = ... ## histogram template 
#  >>> tproject ( tree , histo , 'mass' , 'pt>10' , maxentries = 1000000 )
#  >>> tree.pproject ( histo , 'mass' , 'pt>10' ) ## ditto 
#  >>> tree.tproject ( histo , 'mass' , 'pt>10' ) ## ditto 
#  @endcode
#  - significant gain can be achieved for very large ttrees with complicated expressions and cuts
#  - <code>maxentries</code> parameter should be rather large
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
def  tproject ( tree , histo , what , cuts , nentries = -1 , first = 0 , maxentries = 1000000 ) :
    """Make a projection of the loooong tree into histogram
    >>> tree  = ... ## large chain
    >>> histo = ... ## histogram template 
    >>> tproject ( tree , histo , 'mass' , 'pt>10' )    
    >>> tree.pproject ( histo , 'mass' , 'pt>10' )    ## ditto 
    >>> tree.tproject ( histo , 'mass' , 'pt>10' )    ## ditto 
    - significant gain can be achieved for very large ttrees with complicated expressions and cuts
    - maxentries parameter should be rather large
    """
    #
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
    
    if 0 >  nentries : nentries   = n_large 

    maxentries = long(maxentries)
    if 0 >= maxentries : maxentries = n_large 
    
    if 0 > first    : first     = 0
    
    ## use the regular projection  
    from Ostap.TreeDeco import _tt_project_ 

    if isinstance ( tree , ROOT.TChain ) :
        logger.warning ('``tproject'' method is TTree-specific, skip parallelization')
        return _tt_project_ ( tree , histo , what , cuts , '' , total , first ) 
        
    ## total number of events to process :
    total = min ( num - first , nentries ) 

    ## the event range is rather short, no need  in parallel processing
    if total < maxentries : 
        return _tt_project_ ( tree ,  histo , what , cuts , '', total , first  )
    
    ## check if tree is file-resident:
    tdir = tree.GetDirectory()
    if tdir and isinstance ( tdir , ROOT.TFile ) : pass
    else :
        return _tt_project_ ( tree ,  histo , what , cuts , '', total , first  )

    fname = tdir.GetName()
    
    ## number of jobs & reminder 
    njobs, rest = divmod ( total , maxentries )
    csize       = int ( total / njobs ) ## chunk size 

    ## final list of parameters [ (file_name, first_event , num_events ) , ... ] 
    params = []
    for i in range(njobs) : 
        params.append ( ( fname , first +     i * csize , csize ) )
        
    if rest :
        params.append ( ( fname , first + njobs * csize , rest  ) )
        njobs +=  1

    task  = ProjectTask          ( tree , histo , what , cuts )
    wmgr  = Parallel.WorkManager ()
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
    """
    The single task object for more efficient fill of RooDataSet from TChain 
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
            from   ostap.roofit.selectors import SelectorWithVars
            
        selector = SelectorWithVars ( self.variables ,
                                      self.selection ,
                                      silence = True )
        
        tree,fname = params
        tree = ROOT.TChain( tree )
        tree.Add ( fname )

        tree.process ( selector , 1000 )
        self.output =  selector.data
        
        del selector.data
        del      selector
        
        logger.debug ( 'Processed %s chain and filled %d entries ' % ( fname , len( self.output ) ) ) 

    
    def finalize ( self ) : pass 

    ## methge resulsts/datasets 
    def _mergeResults(self, result) :
        #
        if not isinstance ( self.output , ROOT.RooDataSet ) :  
            self.output = result
        else :
            self.output.append ( result )
            result.Delete () 
            if result : del result
            
        logger.debug ( 'Merging: %d entries ' % len( self.output ) )


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
def  fillDataSet ( chain , variables , selection , ppservers = () ) :
    """Fill dataset from loooong TChain using per-file parallelisation
    >>> chain =
    >>> vars  = ...
    >>> dset  = fillDataSet ( chain , vars , 'pt>10' )
    - for 12-core machine, clear speed-up factor of about 8 is achieved 
    """
        
    task  = FillTask ( variables , selection )
    wmgr  = Parallel.WorkManager( ppservers = ppservers )
    
    cname = chain.GetName() 
    files = chain.files() 
    pairs = [ ( cname,i ) for i in files ] 

    wmgr.process( task, pairs )

    return task.output

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    if not ( 2**32 - 1 ) <= n_large <= ROOT.TVirtualTreePlayer.kMaxEntries :
        logger.error ( "Invalid setting of ``n_large''(%d) parameter (>%d)" % ( n_large , ROOT.TVirtualTreePlayer.kMaxEntries ) )
        
# =============================================================================
# The END 
# =============================================================================
