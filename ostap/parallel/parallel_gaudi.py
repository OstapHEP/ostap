#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/parallel/parallel_gaudi.py
# This is a modified verison of the
# original <code>GaudiMP.Parallel</code> module coded by Pere MATO
# @author Pere Mato (pere.mato@cern.ch)
# 
# =============================================================================
"""
This is a modified verison of the `GaudiMP.Parallel` module by Pere MATO

GaudiMP.Parallel module:
- This module provides 'parallel' processing support for GaudiPyhton.
It is adding some sugar on top of public domain packages such as
the 'multiprocessing' or the 'pp' packages. The interface can be made
independent of the underlying implementation package.
Two main class are defined: Task and WorkManager
"""
# =============================================================================
__all__ = (
    'Task'        , ## base-class for task 
    'WorkManager' , ## Task-manager 
    'Checker'     , ## check of the object can be pickled/unpickled  
   )
# =============================================================================
from   collections.abc          import Sized
from   itertools                import repeat , count
from   ostap.core.ostap_types   import sized_types
from   ostap.utils.progress_bar import progress_bar
from   ostap.parallel.task      import Task, TaskManager
from   ostap.io.checker         import PickleChecker as Checker
from   ostap.parallel.utils     import init_worker_modules 
import multiprocessing          as     MP
import sys, os
# =============================================================================
from    ostap.logger.logger       import getLogger
logger  = getLogger('ostap.parallel.parallel_gaudi')
# =============================================================================

# =============================================================================
class pool_context :
    def __init__  ( self , pool ) :
        self.__pool = pool
    def __enter__ ( self ) :
        sys.stdout .flush ()
        sys.stderr .flush ()
        return self.__pool
    def __exit__  ( self, *_ ) :
        self.__pool.close ()
        self.__pool.join  ()
        sys.stdout .flush ()
        sys.stderr .flush ()
 
# =============================================================================
class WorkManager(TaskManager) :
    """ Class to in charge of managing the tasks and distributing them to
        the workers. They can be local (using other cores) or remote
        using other nodes in the local cluster """

    def __init__( self                            , 
                  ncpus            = 'autodetect' , * , 
                  ppservers        = None         ,
                  pp               = False        ,
                  silent           = True         ,
                  progress         = True         ,
                  block_size       = -1           , 
                  hyper_block_size = -1           ,                   
                  dump_dbase       = None         ,
                  dump_jobs        = 0            ,
                  dump_freq        = 0            ,         
                  **kwargs                        ) :
        

        if pp        : logger.warning ( "WorkManager: option ``pp'' is ignored" )
        if ppservers : logger.warning ( "WorkManager: option ``ppservers'' is ignored" )
        
        ## initialize the base class 
        TaskManager.__init__  ( self ,
                                ncpus            = ncpus            ,
                                silent           = silent           ,
                                progress         = progress         ,
                                block_size       = block_size       ,
                                hyper_block_size = hyper_block_size ,                                
                                dump_dbase       = dump_dbase       ,
                                dump_jobs        = dump_jobs        ,
                                dump_freq        = dump_freq        , **kwargs ) 
        
    # =========================================================================
    ## process the bare <code>executor</code> function
    #  @param job   function to be executed
    #  @param jobs_args the arguments, one entry per job 
    #  @return iterator to results 
    #  @code
    #  mgr  = WorkManager  ( .... )
    #  job  = ...
    #  args = ...
    #  for result in mgr.iexecute ( func , args ) :
    #  ...
    #  ... 
    #  @endcode
    #  It is a "bare minimal" interface
    #  - no statistics
    #  - no summary printout 
    #  - no merging of results  
    def iexecute ( self , job , jobs_args , * ,
                   ordered  = False ,
                   progress = False , **kwargs ) :
        """ Process the bare `executor` function
        >>> mgr  = WorkManager  ( .... )
        >>> job  = ...
        >>> args = ...
        >>> for result in mgr.iexecute ( job , args ) :
        ...
        ...
        It is a ``minimal'' interface
        - no statistics
        - no summary prin
        - no merging of results  
        """
        from ostap.utils.cidict import cidict, cidict_fun 
        myargs = cidict ( self.params , transform = cidict_fun )
        myargs.update   ( kwargs      )
        ##

        ## number of jobs 
        njobs = ( myargs.pop ( 'njobs'     , None ) or 
                  myargs.pop ( 'max_value' , None ) or
                  ( len ( jobs_args ) if isinstance ( jobs_args , sized_types ) else None ) )
        
        chunk_size = myargs.pop ( 'chunk_size' , None ) or myargs.pop ( 'batch_size' , None )
        if   ordered and chunk_size is None : pass
        elif not isinstance ( chunk_size , int ) or chunk_size < 1 :
            chunk_size = self.chunksize_guess ( jobs_args , njobs ) 
            logger.debug ( "`chunksize' is chosen to be %s'" % chunk_size )
            
        ## block-size is embedded deep into multiprocessing 
        myargs.pop ( 'block_size' , self.block_size )

        ## progress-bar description
        description = myargs.pop ( 'description' , "Jobs:" )
        
        
        if myargs : self.extra_arguments ( **myargs ) 
        
        modules_to_import = myargs.pop ( "imports" , [] )
        if isinstance ( modules_to_import , str ) : modules_to_import = [ modules_to_import ]
  
        with MP.Pool ( self.ncpus ) as pool : ##  pool_context ( self.pool ) as pool :
            
            if modules_to_import:
                ## import requested modules on worked nodes  
                n_nodes = pool.nodes if pool.nodes else self.ncpus
                pool.map ( lambda _: init_worker_modules ( modules_to_import ) , range ( n_nodes ) )
                
            ## create and submit jobs
            if ordered : jobs = pool.imap           ( job , jobs_args , chunksize = chunk_size )
            else       : jobs = pool.imap_unordered ( job , jobs_args , chunksize = chunk_size )
            
            ## retrive (asynchronous) results from the jobs
            for result in progress_bar ( jobs                       ,
                                         max_value   = njobs        ,
                                         description = description  , 
                                         silent      = not progress ) :
                yield result                
                
    # ========================================================================
    ## get PP-statistics if/when possible 
    def get_pp_stat ( self ) : 
        """ Get PP-statistics if/when possible 
        """
        return None

    # =========================================================================
    ## context protocol: ENTER 
    def __enter__  ( self      ) :
        sys.stdout .flush ()
        sys.stderr .flush ()
        return self

    # =========================================================================
    ## context protocol: close/join/clear the pool 
    def __exit__   ( self , *_ ) :        
        ##  if  self.pool :
        ##    self.pool.close()
        ##    self.pool.join  ()
        sys.stdout .flush ()
        sys.stderr .flush ()
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )    
    logger.info ("Module ``%s'' is used for multiprocessing" % MP.__name__ )
        
# =============================================================================
##                                                                      The END 
# =============================================================================
