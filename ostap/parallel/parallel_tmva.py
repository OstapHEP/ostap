#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file ostap/parallel/parallel_tmva.py
#   (some) parallelisation  for TMVA treatment
#  - adding of TMVA response to looong ROOT.TChain
#  @see ostap.tools.tmva 
#  @date   2019-05-18
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.
# =============================================================================
""" Add (some) parallelisation  for TMVA treatment
- Adding of TMVA response to looong ROOT.TChain
- see ostap.tools.tmva 
"""
# =============================================================================
__author__  = 'Vanya BELYAEV  Ivan.Belyaev@itep.ru'
__date__    = "2019-05-18"
__version__ = '$Revision$'
__all__     = (
    "addTMVAResponse" , ## add TMVA response to looong ROOT.TChain 
    )
# =============================================================================
from   ostap.parallel.parallel         import Task, WorkManager
from   ostap.utils.basic               import typename 
import ostap.parallel.parallel_statvar
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.tmva'       )
else                       : logger = getLogger ( __name__     )
# =============================================================================
## @class AddTMVA
#  Add TMVA response to looong TChain 
class AddTMVA(Task) :
    """ Add TMVA response to looong TChain 
    """
    def __init__          ( self , *args , **kwargs ) :
        self.  args   =   args
        self.kwargs   = kwargs
        self.__output = ()
        
    def initialize_local  ( self )         : self.__output = () 
    def process           ( self , jobid , trees ) :
        
        import ostap.trees.trees
        from   ostap.tools.tmva  import addTMVAResponse  as add_response       
        
        if not isinstance ( trees , ( tuple , list ) ) : trees = [ trees ]

        files = set()
        
        for tree in trees : 
            add_response ( tree.chain , *self.args , **self.kwargs )
            for f in tree.files : files.add ( f )
            
        ## list of processed  files 
        self.__output = list ( files )
        
        return self.__output 
        
    ## merge results/datasets 
    def merge_results ( self , result , jobid = -1 ) :
        if not  self.__output : self.__output = result
        else :
            s = set()
            for r in self.__output : s.add ( r )
            for r in      result   : s.add ( r )
            s = list ( s )
            s.sort()
            self.__output = tuple( s ) 

    ## get the results 
    def  results ( self ) :  return self.__output
    
# =============================================================================
## Helper function to add TMVA response to loooong TChain 
#  @code
#  tar_file = trainer.tar_file
#  chain    = = ...
#  inputs   = [ 'var1' , 'var2' , 'var2' ]
#  addTMVAResponse ( chain ,  inputs , tar_file , prefix = 'tmva_' )
#  @endcode
#  @param dataset  input dataset to be updated 
#  @param inputs   input variables
#  @param weights_files files with TMVA weigths (tar/gz or xml)
#  @param prefix   prefix for TMVA-variable
#  @param suffix   suffix for TMVA-variable
#  @param options  options to be used in TMVA Reader
#  @param verbose  verbose operation?
#  @param aux       obligatory for the cuts method, where it represents the efficiency cutoff
def addTMVAResponse ( chain                         , ## input chain 
                      inputs                        , ## input variables 
                      weights_files                 , ## files with TMVA weigths (tar/gz or xml)
                      spectators      = ()          , ## spectator variables                       
                      prefix          = 'tmva_'     , ## prefix for TMVA-variable 
                      suffix          = '_response' , ## suffix for TMVA-variable
                      options         = ''          , ## TMVA-reader options
                      verbose         = False       , ## verbosity flag 
                      aux             = 0.9         , ## ## for Cuts method : efficiency cut-off
                      progress        = True        , ## show progress ? 
                      report          = True        , **kwargs ) : ## mare a report? 
    """
    Helper function to add TMVA  response into loong TChain
    >>> tar_file = trainer.tar_file
    >>> dataset  = ...
    >>> inputs = [ 'var1' , 'var2' , 'var2' ]
    >>> dataset.addTMVAResponse (  inputs , tar_file , prefix = 'tmva_' )
    """
    from ostap.tools.tmva import addTMVAResponse as _add_response_

    assert prefix or suffix , 'addTMVAResponse: invalid prefix/suffix %s/%s' % ( prefix , suffix ) 
    
    if isinstance ( chain , ROOT.TTree ) :
        import ostap.trees.trees        
        vars    = set ( chain.branches() ) | set ( chain.leaves () ) 
        matched = sorted ( v for v in vars if v.startswith ( prefix ) and v.endswith (  suffix ) ) 
        if matched :
            matched = ','.join ( matched ) 
            logger.warning ( "addTMVAResponse:: Variables '%s' already in TTree, skip" % matched ) 
            return chain
        
    if isinstance ( chain , ROOT.RooDataSet ) or ( isinstance ( chain , ROOT.TTree ) and chain.nFiles <= 1 )  : 
        return _add_response_ ( dataset       = chain         ,
                                inputs        = inputs        ,
                                weights_files = weights_files ,
                                spectators    = spectators    , 
                                prefix        = prefix        ,
                                suffix        = suffix        ,
                                verbose       = verbose       ,
                                aux           = aux           ,
                                progress      = progress      , 
                                report        = report        )
    
    # =========================================================================
    ## TTree or TChain here
    # =========================================================================
    assert isinstance ( chain , ROOT.TTree ) , "Invalid type of `chain` %s" % typename ( chain)

    from ostap.trees.trees import Chain
    ch       = Chain ( chain )
    treepath = ch.name 
    branches = set  ( chain.branches() ) | set ( chain.leaves() ) if report else set() 
    
    ## create the task 
    task  = AddTMVA     ( inputs        = inputs        ,
                          weights_files = weights_files ,
                          spectators    = spectators    , 
                          prefix        = prefix        ,
                          suffix        = suffix        ,
                          verbose       = False         ,
                          prorgess      = False         ,
                          report        = False         , 
                          aux           = aux           )
    
    wmgr  = WorkManager ( silent = False , progress = progress , **kwargs )
    trees = ch.split    ( max_files = 1  )
    
    wmgr.process( task , trees )
    
    nc = ROOT.TChain ( treepath )
    for f in ch.files :  nc.Add ( f )
    
    if report :
    
        new_branches = sorted ( ( set ( nc.branches () ) | set ( nc.leaves () ) ) - branches )
        if new_branches :
            n = len ( new_branches )
            if 1 >= n : title = "Added %s branch to TTree(%s)"   % ( n , treepath ) 
            else      : title = "Added %s branches to TTree(%s)" % ( n , treepath )  
            table = cn.table ( new_branches , title = title , prefix = '# ' )
            logger.info ( '%s:\n%s' % ( title , table ) ) 
            nc = ROOT.TChain ( ch.name  )
            for f in ch.files :  nc.Add ( f )
    
    return nc

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


    
# =============================================================================
##                                                                      The END 
# =============================================================================
