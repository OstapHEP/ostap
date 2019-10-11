#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file ostap/parallel/parallel_chopping.py
#  Add (some) parallelisation  for TMVA/Chopping treatment
#  - Adding of TMVA/Chopping response to looong ROOT.TChain
#  @see ostap.tools.chopping 
#  @date   2019-05-18
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.
# =============================================================================
""" Add (some) parallelisation  for TMVA/Chopping treatment
- Adding of TMVA/Chopping response to looong ROOT.TChain
- see ostap.tools.chopping
"""
# =============================================================================
__author__  = 'Vanya BELYAEV  Ivan.Belyaev@itep.ru'
__date__    = "2019-05-18"
__version__ = '$Revision$'
__all__     = (
    "addChoppingResponse" , ## add TMVA/Chopping response to looong ROOT.TChain 
    )
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.chopping' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
import ROOT
from   ostap.parallel.parallel import Task, WorkManager
# =============================================================================
## @class AddChopping
#  Add TMVA/Chopping response to looong TChain 
class AddChopping(Task) :
    """Add TMVA/Chopping response to looong TChain 
    """
    def __init__          ( self , *args , **kwargs ) :
        self.  args   =   args
        self.kwargs   = kwargs
        self.__output = ()
        
    def initialize_local  ( self ) : self.__output = () 
    def process           ( self , trees ) :

        import ostap.trees.trees
        from   ostap.tools.chopping  import addChoppingResponse as add_response 
        
        if not isinstance ( trees , ( tuple , list ) ) : trees = [ trees ]
        
        files = set()

        for tree in trees :
            add_response ( tree.chain , *self.args , **self.kwargs )
            for f in tree.files : files.add ( f )
            
        ## list of processed  files 
        self.__output = list ( files )

        return self.__output 

    ## merge results/datasets 
    def merge_results ( self , result ) :
        if not  self.__output : self.__output = result
        else :
            s = set()
            for r in self.__output : s.add ( r )
            for r in      result   : s.add ( r )
            s = list ( s )
            s.sort()
            self.__output = tuple( s )
            
    ## get the results 
    def results  ( self ) : return self.__output 

# ===================================================================================
## @class ChopperTraining
#  parallel procession for TMVA chopping
#  @see ostap/tools/chopping.py
class ChopperTraining(Task) :
    def __init__          ( self ) : self.__output = ()
    def initialize_local  ( self ) : self.__output = () 
    def process           ( self , params ) :
        
        import ROOT, ostap.tools.tmva        
        category , chopper = params
        from   ostap.utils.utils import batch
        from   sys               import version_info as python_version 
        in_batch = 2 < python_version.major or 0 != category 
        with batch ( in_batch ) : 
            trainer  = chopper.create_trainer ( category , False )
            trainer.train ()
        self.__output = (
            [ ( category , trainer.weights_files ) ] ,
            [ ( category , trainer.  class_files ) ] ,
            [ ( category , trainer. output_file  ) ] ,
            [ ( category , trainer.    tar_file  ) ] ,
            [ ( category , trainer.     dirname  ) ] ,
            [ ( category , trainer.    log_file  ) ] ,
            )

        return self.__output
    
    ## merge results/datasets 
    def merge_results ( self , result) :
        if not  self.__output : self.__output =  result
        else :
            weights  = list ( self.__output[0] ) + list ( result[0] ) 
            classes  = list ( self.__output[1] ) + list ( result[1] ) 
            outputs  = list ( self.__output[2] ) + list ( result[2] ) 
            tarfiles = list ( self.__output[3] ) + list ( result[3] ) 
            dirnames = list ( self.__output[4] ) + list ( result[4] ) 
            logfiles = list ( self.__output[5] ) + list ( result[5] ) 
            weights  . sort ()
            classes  . sort ()
            outputs  . sort ()
            tarfiles . sort ()
            logfiles . sort ()
            self.__output = weights , classes , outputs , tarfiles, dirnames , logfiles  

    ## get the results 
    def results  ( self ) : return self.__output 

# =============================================================================
## Helper function to add TMVA/chopping response into loong TChain
#  @code
#  tar_file = trainer.tar_file
#  chain    = ...
#  inputs   = [ 'var1' , 'var2' , 'var2' ] ## input variables to TMVA
#  addChoppingResponse ( chain , chopper , inputs , tar_file , prefix = 'tmva_' )
#  @endcode
#  @param dataset input dataset to be updated
#  @param chopper       chopping category/formula
#  @param N             number of categories
#  @param inputs        input variables
#  @param weights_files files with TMVA weigths (tar/gz or xml)
#  @param category_name the category
#  @param prefix        prefix for TMVA-variable
#  @param suffix        suffix for TMVA-variable
#  @param options       options to be used in TMVA Reader
#  @param verbose       verbose operation?
#  @param aux           obligatory for the cuts method, where it represents the efficiency cutoff 
def addChoppingResponse ( chain                       , ## input dataset to be updated
                          chopper                     , ## chopping category/formula 
                          N                           , ## number of categrories
                          inputs                      , ## input variables 
                          weights_files               , ## files with TMVA weigths (tar/gz or xml)
                          category_name = 'chopping'  , ## category name 
                          prefix        = 'tmva_'     , ## prefix for TMVA-variable         
                          suffix        = '_response' , ## suffix for TMVA-variable 
                          options       =  ''         , ## TMVA-reader options
                          verbose       = True        , ## verbosity flag 
                          aux           = 0.9         ) :
    """
    Helper function to add TMVA/chopping  response into dataset
    >>> tar_file = trainer.tar_file
    >>> chain    = ...
    >>> inputs   = [ 'var1' , 'var2' , 'var2' ] ## input varibales to TMVA 
    >>> addChoppingResponse ( chain , chopper ,  inputs , tar_file , prefix = 'tmva_' )
    """

    from ostap.tools.chopping import addChoppingResponse as _add_response_
    
    if isinstance ( chain , ROOT.TChain ) and 1 < len ( chain.files () ) : pass
    else : return _add_response_ ( dataset       = chain         ,
                                   chopper       = chopper       ,
                                   N             = N             ,
                                   inputs        = inputs        , 
                                   weights_files = weights_files ,
                                   prefix        = prefix        ,
                                   suffix        = suffix        ,
                                   options       = options       , 
                                   verbose       = verbose       ,
                                   aux           = aux           )
    
    from ostap.trees.trees import Chain
    ch    = Chain ( chain )
    
    task  = AddChopping ( chopper       = chopper       ,
                          N             = N             ,
                          inputs        = inputs        , 
                          weights_files = weights_files ,
                          prefix        = prefix        ,
                          suffix        = suffix        ,
                          options       = options       , 
                          verbose       = verbose       ,
                          aux           = aux           )
    
    wmgr  = WorkManager ( silent = False )
    trees = ch.split    ( max_files = 1  )
    
    wmgr.process ( task , trees )
    
    nc = ROOT.TChain ( chain.name )
    for f in ch.files :  nc.Add ( f )
    
    return nc

# =============================================================================
## Perform parallel training of TMVA/Chopping
#  - internal  function for ostap.tools.chopping.Trainer
#  @see  ostap.tools.chopping.Trainer
def chopping_training ( chopper , ncpus = 'autodetect' , ppservers = () ) :
    """Perform parallel traning of TMVA/Chopping
    - internal  function for ostap.tools.chopping.Trainer
    - see  ostap.tools.chopping.Trainer
    """

    import sys
    
    task = ChopperTraining ()
    wmgr = WorkManager ( silent = False , ncpus = ncpus , ppservers = ppservers )
    
    params = [ ( i , chopper ) for i in range ( chopper.N ) ]
    
    sys.stdout.flush()
    sys.stderr.flush()
    
    wmgr.process ( task , params )
    
    sys.stdout.flush()
    sys.stderr.flush()
    
    return task.results() 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    
# =============================================================================
##                                                                      The END 
# =============================================================================
