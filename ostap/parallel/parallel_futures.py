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
from   itertools                    import repeat , count, islice 
from   ostap.utils.progress_bar     import progress_bar
from   ostap.utils.utils            import chunked 
from   ostap.parallel.task          import Task, TaskManager, keyboard_interrupt 
from   ostap.io.checker             import PickleChecker as Checker
from   ostap.core.ostap_types       import sized_types
import concurrent.futures
import sys
#  =============================================================================
from   ostap.logger.logger          import getLogger
logger  = getLogger('ostap.parallel.parallel_futures')
# =============================================================================
## Top-level helper function for batch execution (pickle-safe).
def _process_chunk_ ( fn , chunk ) :
  """ Top-level helper function for batch execution (pickle-safe).
  """
  return [ fn ( item ) for item in chunk ]
# =============================================================================
## Iterative Unordered Map (iumap) with backpressure control.
#
#  Compatible with Python 3.9+.
#  
#  Yields results as soon as tasks complete, ignoring original input order.
#  Prevents memory exhaustion on large streams by limiting the internal task
#  queue.
#  
# @param executor ProcessPoolExecutor or ThreadPoolExecutor instance.
# @param fn       Worker function to apply to each element.
# @param iterable Input sequence or stream of arguments.
# @param chunksize Number of items packed per task (IPC optimization).
# @param block_size  Maximum number of active futures in executor queue
#                  (Backpressure limiat).
def future_uimap ( executor   , 
                   fn         , 
                   iterable   , 
                   chunksize  ,
                   block_size ) : 
    """ Iterative Unordered Map (iumap) with backpressure control.
    
    Compatible with Python 3.9+.
    
    Yields results as soon as tasks complete, ignoring original input order.
    Prevents memory exhaustion on large streams by limiting the internal task
    queue.
    
    :param executor: ProcessPoolExecutor or ThreadPoolExecutor instance.
    :param fn: Worker function to apply to each element.
    :param iterable: Input sequence or stream of arguments.
    :param chunksize: Number of items packed per task (IPC optimization).
    :param block_size: Maximum number of active futures in executor queue
    (Backpressure limit).
    """
    assert isinstance ( chunksize  , int ) and  1<= chunksize, \
        "future_iumap:  Invalid `chunksize': %s" % chunksize
    assert isinstance ( block_size , int ) and 1 < block_size, \
        "future_iumap:  Invalid `block_size': %s" % block_size
    
    chunks_iter = chunked ( iterable  , chunksize )

    # ==========================================================================
    # Track active task Futures currently running in the executor
    # ==========================================================================
    pending_futures = set()

    # ==========================================================================
    # Helper to submit tasks safely for both Thread and Process pools
    # ==========================================================================
    def _submit_chunk_ ( chunk ) :
        if 1 == chunksize :  return executor.submit ( fn , chunk [ 0 ] )
        return executor.submit ( _process_chunk_ , fn , chunk )

    # ==========================================================================
    # Step 1: Pre-fill the execution queue up to `block_size` limit
    # ==========================================================================
    for chunk in islice ( chunks_iter , block_size ):
        pending_futures.add ( _submit_chunk_ ( chunk ) )

    # ==========================================================================
    # Step 2: Stream completed results and dispatch new tasks dynamically
    # ==========================================================================
    while pending_futures:
        # Wait for any single task to finish (as-completed pattern)
        done_future = next ( concurrent.futures.as_completed ( pending_futures ) )
        pending_futures.remove ( done_future )

        # ======================================================================
        # Submit a replacement task from the stream to maintain queue depth
        # ======================================================================        
        try: # =================================================================
            # ==================================================================
            next_chunk = next ( chunks_iter )
            pending_futures.add ( _submit_chunk_ ( next_chunk ) )
            # ==================================================================
        except StopIteration: # ================================================
            # ==================================================================
            pass  # Input stream fully consumed

        # ======================================================================        
        # Yield unpacked results to the consumer
        # ======================================================================        
        result = done_future.result()
        if 1 == chunksize : yield result
        else:
            yield from result

# ==============================================================================
## @class WorkManager
#  Class to in charge of managing the tasks and distributing them to
#  the workers. They can be local (using other cores) or remote
#  using other nodes in the local cluster """
class WorkManager(TaskManager) :
    """ Class to in charge of managing the tasks and distributing them to the workers.
    """
    def __init__( self ,
                  ncpus            = 'autodetect', * , 
                  silent           = False       ,
                  progress         = True        ,
                  block_size       = -1          , 
                  hyper_block_size = -1          ,                                                       
                  dump_dbase       = None        ,
                  dump_jobs        = 0           ,
                  dump_freq        = 0           ,  **kwargs ) :

        ## 
        if 'ppservers' in kwargs: conf.pop ( 'ppservers' )        
        ## initialize the base class 
        TaskManager.__init__  ( self ,
                                ncpus            = ncpus      ,
                                silent           = silent     ,
                                progress         = progress   ,
                                block_size       = block_size       ,
                                hyper_block_size = hyper_block_size ,                                
                                dump_dbase       = dump_dbase ,
                                dump_jobs        = dump_jobs  ,
                                dump_freq        = dump_freq  , **kwargs ) 
        
    # =====================================================================
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
        - no summary print
        - no merging of results  
        """
        from ostap.utils.cidict import cidict, cidict_fun 
        myargs = cidict ( self.params , transform = cidict_fun )
        myargs.update   ( kwargs      )

        block_size = myargs.pop ( 'buffer_size' , None ) or myargs.pop ( 'block_size' , None  ) or self.block_size 
        if not isinstance ( block_size , int ) or block_size <= 1 : block_size = self.block_size

        chunk_size = myargs.pop ( 'batch_size'  , None ) or myargs.pop ( 'chunk_size' , None  )
        if not isinstance ( chunk_size , int ) or chunk_size < 1 :
            chunk_size = self.chunksize_guess ( jobs_args ) 
            logger.info ( "`chunksize' is chosen to be %s" % chunk_size )
            
        ## progress-bar description
        description = myargs.pop ( 'description' , "Jobs:" )
        
        ## number of jobs 
        njobs = ( myargs.pop ( 'njobs'     , None ) or 
                  myargs.pop ( 'max_value' , None ) or
                  ( len ( jobs_args ) if isinstance ( jobs_args , sized_types ) else None ) )

        ## 
        progress = progress    or self.progress        
        silent   = self.silent or not progress
        ##
        done   = 0
        
        config = dict ( chunksize = chunk_size )
        if not ordered                 : config [ 'block_size'  ] = block_size
        elif ( 3 , 14 ) <= python_info : config [ 'buffer_size' ] = block_size

        # =====================================================================
        ## creaet and use executor:
        with self.make_executor ( max_workers = self.ncpus ) as executor:

            # =================================================================
            try : # ===========================================================
                # =============================================================

                if ordered : results = executor.map (            job , jobs_args , **config )
                else       : results = future_uimap ( executor , job , jobs_args , **config ) 
                
                for result in progress_bar ( results                   ,
                                            max_value   = njobs       ,
                                             description = description ,
                                             silent      = silent      ) : 
                    yield result
                    done +=1
                # ============================================================
            except KeyboardInterrupt : # =====================================
                # ============================================================
                logger.attention ( "%s only #%d jobs are processed" % ( keyboard_interrupt , done ) )
                # ===========================================================
                return
                # ============================================================ 
            except Exception : # =============================================
                # ============================================================
                logger.error ( 'Exception caught after #%d jobs processed' % done , exc_info = True )
                raise   
            
        if kwargs : self.extra_arguments ( **kwargs ) 
        
    # =========================================================================
    ## helper method to creat ethe executor
    def make_executor ( self , *args , **kwargs ) :
        """ Helper method to creat ethe executor"""
        return concurrent.futures.ProcessPoolExecutor ( *args , **kwargs ) 
                
    # =========================================================================
    
    ## get PP-statistics if/when possible 
    def get_pp_stat ( self ) : 
        """ Get PP-statistics if/when possible 
        """
        return None

    # =========================================================================
    ## context protocol: ENTER 
    def __enter__  ( self ) :
        """ Context protocol: ENTER"""
        sys.stdout .flush ()
        sys.stderr .flush ()
        return self
    
    # =========================================================================
    ## context protocol: EXIT
    def __exit__   ( self , *_ ) :        
        """ Context protocol: EXIT"""
        sys.stdout .flush ()
        sys.stderr .flush ()        

# =============================================================================
if '__main__' == __name__ : # =================================================
    # =========================================================================
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )    

# =============================================================================
##                                                                      The END 
# =============================================================================
