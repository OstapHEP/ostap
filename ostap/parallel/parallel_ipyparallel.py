#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/parallel/parallel_ipyparallel.py
# Task manager for ipyparallel processinng 
# =============================================================================
""" Task manager for ipyparallel processinng 
"""
# =============================================================================
__all__ = () 
# =============================================================================
import sys, os, time
from   itertools                    import repeat , count
from   ostap.utils.progress_bar     import progress_bar
from   ostap.parallel.task          import Task, TaskManager, keyboard_interrupt 
from   ostap.io.checker             import PickleChecker as Checker
from   ostap.core.ostap_types       import sized_types 
# =============================================================================
from   ostap.logger.logger          import getLogger
logger  = getLogger('ostap.parallel.parallel_ipyparallel')
# =============================================================================
## Try to import ipyparallel 
ipp = None
# =============================================================================
try : # =======================================================================
    import ipyparallel as ipp
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    ipp = None
# =============================================================================
## Use only relatively fresh versions of ipyparallel 
if ipp and ( 8 , 0 ) <= ipp.version_info : # ==================================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import dill
        # =====================================================================
    except ImportError : # ====================================================
        # =====================================================================
        dill = None
    # =========================================================================
    ## @class WorkManager
    #  Class to in charge of managing the tasks and distributing them to
    #  the workers. They can be local (using other cores) or remote
    #  using other nodes in the local cluster """
    class WorkManager(TaskManager) :
        """ Class to in charge of managing the tasks and distributing them to
        the workers. They can be local (using other cores) or remote
        using other nodes in the local cluster """
        
        def __init__( self ,
                      ncpus      = 'autodetect' , * ,
                      silent     = False        ,
                      progress   = True         ,
                      balanced   = True         ,
                      use_dill   = True         ,
                      chunk_size = -1           , 
                      dump_dbase = None         ,
                      dump_jobs  = 0            ,
                      dump_freq  = 0            , **kwargs ) :


            if 'ppservers' in kwargs : kwarsg.pop ( 'ppservers' )
              
            ## initialize the base class 
            TaskManager.__init__  ( self ,
                                    ncpus      = ncpus      ,
                                    silent     = silent     ,
                                    progress   = progress   ,
                                    chunk_size = chunk_size , 
                                    dump_dbase = dump_dbase ,
                                    dump_jobs  = dump_jobs  ,
                                    dump_freq  = dump_freq  , **kwargs ) 
            
            if self.silent :
                import logging 
                kwargs [ 'log_level' ] = logging.WARNING 

            ##  ipp.Cluster arguments 
            self.__kwargs   = kwargs
            self.__balanced = True if          balanced else False 
            self.__use_dill = True if dill and use_dill else False
            
            if not self.silent :
                logger.info ( 'WorkManager is ipyparallel'  )
                if self.__kwargs :
                    rows = [ ( 'Parameter' , 'value' ) ]
                    for k, v in self.__kwargs.items()  :
                        row = k , str ( v  )
                        rows.append ( row )
                    import ostap.logger.table as T
                    title = 'ipyparallel'
                    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'll' )
                    logger.info ( '%s\n%s' % ( title , table ) )
                        
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
            
            njobs    = kwargs.pop ( 'njobs' , kwargs.pop ( 'max_value' , len ( jobs_args ) if isinstance ( jobs_args , sized_types ) else None ) )
            
            progress = progress    or self.progress        
            silent   = self.silent or not progress
            done     = 0

            with ipp.Cluster ( n = self.ncpus , **self.__kwargs ) as cluster :

                if   self.__use_dill :
                    view = cluster[:]                    
                    view.use_dill ()
                    
                ## BALANCED ? 
                view = cluster.load_balanced_view() if self.__balanced else cluster[:]

                # ================================================================
                try : # ==========================================================
                    # ============================================================                    
                    ## results = view.map_async ( job , jobs_args )                    
                    ## results = view.imap ( job , jobs_args , block = False )                    
                    results = view.imap ( job , jobs_args )                    
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
                    ## ABORT!! 
                    view.abort()                    
                    return
                    # ============================================================ 
                except Exception : # =============================================
                    # ============================================================
                    logger.error ( 'Exception caught after #%d jobs processed' % done , exc_info = True )
                    ## ABORT!! 
                    view.abort()                
                    raise   
           
            if kwargs : self.extra_arguments ( **kwargs ) 
                    
        # ========================================================================
        ## get PP-statistics if/when possible 
        def get_pp_stat ( self ) : 
            """ Get PP-statistics if/when possible 
            """
            return None

        # ========================================================================
        ## Context protocol: ENTER 
        def __enter__  ( self ) :
            """ Context protocol: ENTER: """
            sys.stdout .flush ()
            sys.stderr .flush ()
            return self
    
        # ========================================================================
        ## Context protocol: EXIT 
        def __exit__   ( self , *_ ) :        
            """ Context protocol: EXIT: """
            sys.stdout .flush ()
            sys.stderr .flush ()

    # =========================================================================
    ## define "importable" content of the module:
    __all__ = __all__ + (
        'Task'        , ## Base class for Task 
        'WorkManager' , ## Task-manager 
        )
    
# =============================================================================
if '__main__' == __name__ : # =================================================
    # =========================================================================
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )    
    if ipp : logger.info    ( 'ipyparallel version is %s' % str ( ipp.version_info ) )
    else   : logger.warning ( 'No ipyparallel is available!')
    
# =============================================================================
##                                                                      The END 
# =============================================================================
