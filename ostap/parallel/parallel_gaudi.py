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
    'WorkManager' , ## Task-manager 
    'Checker'     , ## check of the object can be pickled/unpickled  
   )
# =============================================================================
from   collections.abc          import Sized
from   itertools                import repeat , count
from   ostap.utils.progress_bar import progress_bar
from   ostap.parallel.task      import TaskManager
from   ostap.io.checker         import PickleChecker as Checker 
import multiprocessing          as     MP
import sys, os, time
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

    def __init__( self                     , 
                  ncpus      = 'autodetect' ,
                  ppservers  = None         ,
                  pp         = False        ,
                  silent     = False        ,
                  progress   = True         ,
                  dump_dbase = None         ,
                  dump_jobs  = 0            ,
                  dump_freq  = 0            ,         
                  **kwargs                  ) :
        
        if not ( isinstance ( ncpus , int ) and 0 <= ncpus ) :
            from ostap.utils.basic import numcpu 
            ncpus = numcpu () 
            
        if isinstance ( ncpus , int ) and 1 <= ncpus : pass
        else                                         : ncpus = MP.cpu_count()

        ## if pp        : logger.warning ( "WorkManager: option ``pp'' is ignored" )
        ## if ppservers : logger.warning ( "WorkManager: option ``ppservers'' is ignored" )
        
        ## initialize the base class 
        TaskManager.__init__  ( self ,
                                ncpus      = ncpus      ,
                                silent     = silent     ,
                                progress   = progress   , 
                                dump_dbase = dump_dbase ,
                                dump_jobs  = dump_jobs  ,
                                dump_freq  = dump_freq  ) 
        
        self.pool   = MP.Pool ( self.ncpus )
        
        if kwargs : self.extra_arguments ( **kwargs )

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
    def iexecute ( self , job , jobs_args , progress = False , **kwargs ) :
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
        
        njobs = kwargs.pop ( 'njobs' , kwargs.pop ( 'max_value' , len ( jobs_args ) if isinstance ( jobs_args , Sized ) else None ) ) 
        with pool_context ( self.pool ) as pool :

            ## create and submit jobs 
            jobs = pool.imap_unordered ( job , jobs_args )
            
            silent = self.silent or not progress
            
            ## retrive (asynchronous) results from the jobs
            for result in progress_bar ( jobs                               ,
                                         max_value   = njobs                ,
                                         description = kwargs.pop ( 'description' , "Jobs:" ) , 
                                         silent      = silent               ) :
                yield result                
                
        if kwargs :
            import ostap.logger.table as T 
            rows = [ ( 'Argument' , 'Value' ) ]
            for k , v in loop_items ( kw ) :
                row = k , str ( v )
                rows.append ( row )
            title = 'iexecute:: %d unused arguments' % len ( kw ) 
            table = T.table ( rows , title = 'Unused arguments' , prefix = '# ' , alignment = 'll' )    
            logger.warning ( '%s\n%s' % ( title , table ) )
            
    # ========================================================================-
    ## get PP-statistics if/when possible 
    def get_pp_stat ( self ) : 
        """ Get PP-statistics if/when possible 
        """
        return None
    
    ## context protocol: restart the pool 
    def __enter__  ( self      ) :
        sys.stdout .flush ()
        sys.stderr .flush ()
        return self
    
    ## context protocol: close/join/clear the pool 
    def __exit__   ( self , *_ ) :        
        if  self.pool :
            self.pool.close()
            self.pool.join  ()
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
