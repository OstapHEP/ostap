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
)
# =============================================================================
from   packaging.version            import Version 
from   ostap.utils.progress_bar     import progress_bar
from   ostap.parallel.task          import Task, TaskManager, keyboard_interrupt 
from   ostap.io.checker             import PickleChecker as Checker
from   ostap.core.ostap_types       import sized_types 
import sys
# =============================================================================
from   ostap.logger.logger          import getLogger
logger  = getLogger('ostap.parallel.parallel_joblib')
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import joblib # ===========================================================
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    joblib = None # ===========================================================
# =============================================================================
## @class WorkManager
#  Class to in charge of managing the tasks and distributing them to
#  the workers. They can be local (using other cores) or remote
#  using other nodes in the local cluster """
class WorkManager(TaskManager) :
    """ Class to in charge of managing the tasks and distributing them to the workers.
    """
    def __init__( self       ,
                  ncpus      = 'autodetect', * , 
                  silent     = False       ,
                  progress   = True        ,
                  chunk_size = -1          , 
                  dump_dbase = None        ,
                  dump_jobs  = 0           ,
                  dump_freq  = 0           , **kwargs ) :

        ##
        if   joblib and Version ( "1.4.0" ) <= Version ( joblib.__version__ ) : 
            kwargs [ 'return_as' ] = kwargs.pop ( 'return_as' , 'generator_unordered' )
        elif joblib and Version ( "1.3.0" ) <= Version ( joblib.__version__ ) : 
            kwargs [ 'return_as' ] = kwargs.pop ( 'return_as' , 'generator' )
        else : 
            kwargs.pop ( 'return_as' )
            if not joblib : logger.error ( "No joblib is available!" ) 

        if 'ppservers' in kwargs : kwargs.pop ( 'ppservers' )
        ## initialize the base class 
        TaskManager.__init__  ( self       ,
                                ncpus      = ncpus      ,
                                silent     = silent     ,
                                progress   = progress   ,
                                chunk_size = chunk_size , 
                                dump_dbase = dump_dbase ,
                                dump_jobs  = dump_jobs  ,
                                dump_freq  = dump_freq  , **kwargs ) 
        
        
        if not self.silent : logger.info ( 'WorkManager is joblib'  )
                
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
        - no summary print
        - no merging of results  
        """

        if not joblib:
            logger.error ( "No joblib module is available, processing is disabled" )
            return
        
        njobs    = kwargs.pop ( 'njobs' , kwargs.pop ( 'max_value' , len ( jobs_args ) if isinstance ( jobs_args , sized_types ) else None ) ) 
        ## 
        progress = progress    or self.progress        
        silent   = self.silent or not progress
        ##
        done = 0
        # ========================================================================
        with joblib.Parallel ( n_jobs = self.ncpus , **self.params )  as executor: 
            # ================================================================ 
            try: # ===========================================================
                # ============================================================
                results = executor ( joblib.delayed ( job ) ( a ) for a in jobs_args ) 
                # ================================================================ 
                for result in progress_bar ( results                            ,
                                             max_value   = njobs                ,
                                             description = kwargs.pop ( 'description' , "Jobs:" ) ,
                                             silent      = silent               ) : 
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
       
        if kwargs : self.extra_arguments ( **kwargs ) 
        
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
