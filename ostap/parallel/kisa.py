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
#  
#                    $Revision$
#  Last modification $Date$
#  by                $Author$
# =============================================================================
""" Multiprocessing functionality for Ostap
Currently it is not loaded on default, and requires manual activation
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'ProjectTask' , 
    'FillTask'    , 
    'project'     , 
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

# =============================================================================
## The simple task object for more efficient projection of loooong TChains
#  into histogarms
#  @see GaudiMP.Parallel
#  @see GaudiMP.Parallel.Task
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
class ProjectTask(Parallel.Task) :
    """ The simple task  object for efficient projection of loooong TChains
    into histogarms  
    """
    ## constructor: chain name, historgam , variable , cuts 
    def __init__ ( self , tree , histo , what , cuts = '' ) :
        """Constructor: chain name, historgam , variable , cuts
        
        >>> histo = ...
        >>> task  = ProjectTask ( 'MyTuple' , histo , 'mass' , 'pt>10' ) 
        """
        
        self.histo = histo
        self.histo.Reset()
        
        if   isinstance ( tree , ROOT.TTree ) : self.tree = tree.GetName()
        elif isinstance ( tree , str        ) : self.tree = tree 
        
        self.what = what
        self.cuts = cuts

    ## local initialization (executed once in parent process)
    def initializeLocal   ( self ) :
        """
        Local initialization (executed once in parent process)
        """
        import ROOT
        import ostap.core.pyrouts 
        self.output = self.histo.Clone() 
        
    ## remote initialization (executed for each sub-processs)
    def initializeRemote  ( self ) : pass 

    ## the actual processing 
    def process ( self , params ) :

        import ROOT
        import ostap.core.pyrouts 
        
        if   isinstance ( params , str                ) : pars = [ params            ]
        elif isinstance ( params , ROOT.TChainElement ) : pars = [ params.GetTitle() ] 
        
        chain = ROOT.TChain ( self.tree )
        for p in pars : chain.Add ( p )

        ## Create the output histogram   NB! (why here???) 
        self.output = self.histo.Clone()
        
        ## use the regular projection
        from ostap.trees.trees import _tt_project_ 
        _tt_project_ ( chain ,  self.output  , self.what , self.cuts )

        del chain 

    ## finalization (execuetd at the end at parent process)
    def finalize ( self ) : pass 

    ## merge results 
    def _mergeResults(self, result) :
        #
        self.output.Add ( result )
        result.Delete () 

# =============================================================================  
## make a projection of the loooong chain into histogram using
#  multiprocessing functionality for per-file parallelisation
#  @code
#
#  >>> chain = ... ## large chain
#  >>> histo = ... ## histogram template 
#  >>> project ( chain , histo , 'mass' , 'pt>10' )
#
#  @endcode
#
#  For 12-core machine, clear speedup factor of about 8 is achieved 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
def  project ( chain , histo , what , cuts ) :
    """Make a projection of the loooong chain into histogram
    >>> chain = ... ## large chain
    >>> histo = ... ## histogram template 
    >>> project ( chain , histo , 'mass' , 'pt>10' )
    
    For 12-core machine, clear speedup factor of about 8 is achieved     
    """
    #
    histo.Reset()    
    files = chain.files() 

    task  = ProjectTask          ( chain , histo , what , cuts )
    wmgr  = Parallel.WorkManager ()
    wmgr.process( task, files )

    histo += task.output  

    return histo


import ROOT 
ROOT.TChain._project = project

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
    
    from ostap import banner
    logger.info ( __file__ + '\n' + banner )
    logger.info ( 80*'*' )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 

# =============================================================================
# The END 
# =============================================================================
