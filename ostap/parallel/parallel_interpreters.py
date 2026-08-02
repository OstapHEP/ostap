#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/parallel/parallel_futures.py
# Task manager for concurrent.futures processing
# =============================================================================
""" Task manager for concurrent.futures processing 
"""
# =============================================================================
__all__ = (
    'Task'        , ## Base class for Task 
    'WorkManager' , ## Task-manager
    'Checker'     , ## serialization checker     
)
# =============================================================================
from   ostap.parallel.task             import Task
import ostap.parallel.parallel_futures as     PF
import sys, concurrent.futures   
#  ============================================================================
from   ostap.logger.logger          import getLogger
logger  = getLogger('ostap.parallel.parallel_interpreters')
# =============================================================================
## serialization checker 
Checker = PF.Checker
# ==============================================================================
## @class WorkManager
#  Class to in charge of managing the tasks and distributing them to
#  the workers. They can be local (using other cores) or remote
#  using other nodes in the local cluster """
class WorkManager(PF.WorkManager) :
    """ Class to in charge of managing the tasks and distributing them to the workers.
    """
    def __init__( self ,
                  ncpus            = 'autodetect', * , 
                  silent           = True        ,
                  progress         = True        ,
                  block_size       = -1          , 
                  hyper_block_size = -1          ,                                                       
                  dump_dbase       = None        ,
                  dump_jobs        = 0           ,
                  dump_freq        = 0           ,  **kwargs ) :

        if 'ppservers' in kwargs: conf.pop ( 'ppservers' )                
        ## initialize the base class 
        PF.WorkManager.__init__  ( self ,
                                   ncpus            = ncpus      ,
                                   silent           = silent     ,
                                   progress         = progress   ,
                                   block_size       = block_size       ,
                                   hyper_block_size = hyper_block_size ,                                
                                   dump_dbase       = dump_dbase ,
                                   dump_jobs        = dump_jobs  ,
                                   dump_freq        = dump_freq  , **kwargs ) 
        
    # =========================================================================
    ## helper method to create the executor
    def make_executor ( self , *args , **kwargs ) :
        """ Helper method to create the executor
        """        
        if sys.version_info < ( 3, 14 ) :
            logger.warning ( "Fallback to ProcessPoolExecutor" )
            return PF.WorkManager.make_executor ( self , *args , **kwargs )
        ## 
        return concurrent.futures.InterpreterPoolExecutor ( *args , **kwargs ) 
    
                
# =============================================================================
if '__main__' == __name__ : # =================================================
    # =========================================================================
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )    

# =============================================================================
##                                                                      The END 
# =============================================================================
