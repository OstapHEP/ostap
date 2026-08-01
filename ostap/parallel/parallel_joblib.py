#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/parallel/parallel_joblib.py
# Task manager for joblb  processing
# =============================================================================
""" Task manager for joblib processing 
"""
# =============================================================================
__all__ = (
    'Task'        , ## Base class for Task 
    'WorkManager' , ## Task-manager
    'Checker'     , ## check of the object can be pickled/unpickled  
)
# =============================================================================
from   packaging.version            import Version 
from   ostap.utils.progress_bar     import progress_bar
from   ostap.parallel.task          import Task, TaskManager, keyboard_interrupt 
from   ostap.io.checker             import PickleChecker as Checker
from   ostap.core.ostap_types       import sized_types 
from   queue                        import Queue
import sys
# =============================================================================
from   ostap.logger.logger          import getLogger
logger  = getLogger('ostap.parallel.parallel_joblib')
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import joblib # ===========================================================
    joblib_version = Version ( joblib.__version__ ) 
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    joblib = None # ===========================================================
# =============================================================================
## @class WorkManager
#  Class to in charge of managing the tasks and distributing them to
#  the workers.
class WorkManager(TaskManager) :
    """ Class to in charge of managing the tasks and distributing them to the workers.
    """
    def __init__( self       ,
                  ncpus            = 'autodetect', * , 
                  silent           = False       ,
                  progress         = True        ,
                  block_size       = -1          , 
                  hyper_block_size = -1          ,                                     
                  dump_dbase       = None        ,
                  dump_jobs        = 0           ,
                  dump_freq        = 0           , **kwargs ) :

        ##
        assert joblib , "No joblib is available!"

        if 'ppservers' in kwargs : kwargs.pop ( 'ppservers' )
        ## initialize the base class 
        TaskManager.__init__  ( self       ,
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
        
        if not joblib:
            logger.error ( "No joblib module is available, processing is disabled" )
            return

        from ostap.utils.cidict import cidict, cidict_fun 
        myargs = cidict ( self.params , transform = cidict_fun )
        myargs.update   ( kwargs      )

        chunk_size = myargs.pop ( 'batch_size'   , None ) or myargs.pop ( 'chunk_size' , None  ) or 'auto'        
        block_size = myargs.pop ( 'pre_dispatch' , None ) or myargs.pop ( 'block_size' , None  ) or '2*n_jobs'

        if 'auto' != chunk_size and ( not isinstance ( chunk_size , int ) or chunk_size < 1 ) :
            chunk_size = self.chunksize_guess ( jobs_args ) 
            logger.info ( "`chunksize' is chosen to be %s'" % chunk_size )
            
        myargs.pop ( 'as_generator' , None )
        
        config = dict ( batch_size = chunk_size , pre_dispatch = block_size )
        if not ordered and Version ( '1.4.0' ) <= joblib_version : config [ 'return_as' ] = 'generator_unordered'
        elif               Version ( '1.3.0' ) <= joblib_version : config [ 'return_as' ] = 'generator'
        
        ## progress-bar description
        description = myargs.pop ( 'description' , "Jobs:" )
        
        ## number of jobs 
        njobs = ( myargs.pop ( 'njobs'     , None ) or 
                  myargs.pop ( 'max_value' , None ) or
                  ( len ( jobs_args ) if isinstance ( jobs_args , sized_types ) else None ) )

        
        if myargs : self.extra_arguments ( **myargs ) 
        
        ## 
        progress = progress    or self.progress        
        silent   = self.silent or not progress
        ##
        done = 0
        # ========================================================================
        with joblib.Parallel ( n_jobs = self.ncpus , **config )  as executor: 
            # ================================================================ 
            try: # ===========================================================
                # ============================================================
                results = executor ( joblib.delayed ( job ) ( a ) for a in jobs_args ) 
                # ================================================================ 
                for result in progress_bar ( results                   ,
                                             max_value   = njobs       ,
                                             description = description ,
                                             silent      = silent      ) : 
                    yield result
                    done += 1
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
       
    # ========================================================================-
    ## get PP-statistics if/when possible 
    def get_pp_stat ( self ) : 
        """ Get PP-statistics if/when possible 
        """
        return None

    # =========================================================================
    ## context protocol
    def __enter__  ( self ) :
        sys.stdout .flush ()
        sys.stderr .flush ()
        return self
    
    # =========================================================================
    ## context protocol
    def __exit__   ( self , *_ ) :        
        sys.stdout .flush ()
        sys.stderr .flush ()        

# =============================================================================
if '__main__' == __name__ : # =================================================
    # =========================================================================
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )    
    
    if not joblib : logger.warning ( "No joblib moduel is available!" ) 
# =============================================================================
##                                                                      The END 
# =============================================================================
