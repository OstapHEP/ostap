#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file ostap/parallel/parallel_tmva.py
#  Add (some) parallelisation  for TMVA treatment
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
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.tmva'       )
else                       : logger = getLogger ( __name__     )
# =============================================================================
import ROOT
from   ostap.parallel.parallel import Task, WorkManager
# =============================================================================
## @class AddTMVA
#  Add TMVA response to looong TChain 
class AddTMVA(Task) :
    """Add TMVA response to looong TChain 
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
def addTMVAResponse ( chain                  ,   ## input chain 
                      inputs                 ,   ## input variables 
                      weights_files          ,   ## files with TMVA weigths (tar/gz or xml)
                      prefix   = 'tmva_'     ,   ## prefix for TMVA-variable 
                      suffix   = '_response' ,   ## suffix for TMVA-variable
                      options  = ''          ,   ## TMVA-reader options
                      verbose  = True        ,   ## verbosity flag 
                      aux      = 0.9         ) : ## for Cuts method : efficiency cut-off
    """
    Helper function to add TMVA  response into loong TChain
    >>> tar_file = trainer.tar_file
    >>> dataset  = ...
    >>> inputs = [ 'var1' , 'var2' , 'var2' ]
    >>> dataset.addTMVAResponse (  inputs , tar_file , prefix = 'tmva_' )
    """
    from ostap.tools.tmva import addTMVAResponse as _add_response_
    
    if isinstance ( chain , ROOT.TChain ) and 1 < len ( chain.files () ) : pass
    else : return _add_response_ ( dataset       = chain         ,
                                   inputs        = inputs        ,
                                   weights_files = weights_files ,
                                   prefix        = prefix        ,
                                   suffix        = suffix        ,
                                   verbose       = verbose       ,
                                   aux           = aux           )

    from ostap.trees.trees import Chain
    ch       = Chain ( chain ) 
    branches = set   ( chain.branches() )
    
    ## create the task 
    task  = AddTMVA     ( inputs        = inputs        ,
                          weights_files = weights_files ,
                          prefix        = prefix        ,
                          suffix        = suffix        ,
                          verbose       = verbose       ,
                          aux           = aux           )
    
    wmgr  = WorkManager ( silent = False )
    trees = ch.split    ( max_files = 1  )
    
    wmgr.process( task , trees )
    
    nc = ROOT.TChain ( chain.name )
    for f in ch.files :  nc.Add ( f )

    nb = list ( set ( nc.branches () ) - branches )
    if nb : logger.info ( 'Added branches:\n%s' % nc.table ( variables = nb , prefix = '# ' ) ) 
    
    return nc

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


    
# =============================================================================
##                                                                      The END 
# =============================================================================
