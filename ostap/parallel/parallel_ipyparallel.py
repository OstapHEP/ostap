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
from   ostap.parallel.task          import Task, TaskManager 
from   ostap.io.checker             import PickleChecker as Checker 
# =============================================================================
from   ostap.logger.logger          import getLogger
logger  = getLogger('ostap.parallel.parallel_ipyparallel')
# =============================================================================
if ( 3 , 3 ) <= sys.version_info  : from collections.abc import Sized
else                              : from collections     import Sized 
# =============================================================================
## Try to import ipyparallel 
ipp = None
# =============================================================================
if ( 3 , 6 ) <= sys.version_info : # ==========================================
    # =========================================================================
    try : # ===================================================================
        import ipyparallel as ipp
        # ====================================================================
    except ImportError : # ===================================================
        # ====================================================================
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
                      ncpus      = 'autodetect',
                      silent     = False ,
                      progress   = True  ,
                      balanced   = False ,
                      use_dill   = True  ,
                      chunk_size = -1    , 
                      dump_dbase = None  ,
                      dump_jobs  = 0     ,
                      dump_freq  = 0     ,                             
                      **kwargs ) :

            if not ( isinstance ( ncpus , int ) and 0 <= ncpus ) :
                from ostap.utils.basic import numcpu 
                ncpus = numcpu() 
                
            n = kwargs.get ( 'n' , ncpus )
            if isinstance ( n , int ) and 1 <= n :    kwargs [ 'n' ] = n
            elif 'n' in kwargs                  : del kwargs [ 'n' ]

            if not isinstance ( chunk_size , int ) or chunk_size <= 1 :
                chunk_size = 5 * ( ncpus + 1 )
              
            ## initialize the base class 
            TaskManager.__init__  ( self ,
                                    ncpus      = ncpus      ,
                                    silent     = silent     ,
                                    progress   = progress   ,
                                    chunk_size = chunk_size , 
                                    dump_dbase = dump_dbase ,
                                    dump_jobs  = dump_jobs  ,
                                    dump_freq  = dump_freq  ) 
            
            if self.silent :
                import logging 
                kwargs [ 'log_level' ] = logging.WARNING 

            ##  ipp.Cluster arguments 
            self.__kwargs   = kwargs
            self.__balanced = True if balanced else False 
            self.__use_dill = True if use_dill else False
            if self.__use_dill and not dill :
                logger.warning ( "dill is not available, switch it off!" ) 
                self.__use_dill = False

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
            
            njobs = kwargs.pop ( 'njobs' , kwargs.pop ( 'max_value' , len ( jobs_args ) if isinstance ( jobs_args , Sized ) else None ) ) 

            silent = self.silent or not progress
            with ipp.Cluster ( **self.__kwargs ) as cluster :

                if   self.__use_dill :
                    view = cluster[:]                    
                    view.use_dill ()
                elif self.__balanced : 
                    view = cluster.load_balanced_view()
                else :
                    view = cluster[:]                    
                        
                results = view.map_async ( job , jobs_args )                    
                for result in progress_bar ( results                            ,
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

        ## context protocol
        def __enter__  ( self      ) :
            sys.stdout .flush ()
            sys.stderr .flush ()
            return self
    
        ## context protocol
        def __exit__   ( self , *_ ) :        
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
