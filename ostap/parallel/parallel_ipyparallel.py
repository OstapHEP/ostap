#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/parallel/parallel_ipyparallel.py
# Task manager for ipyparallel processinng 
# =============================================================================
""" Task manager for ipyparallel processinng 
"""
# =============================================================================
__all__ = (
    'Task'        , ## base-class for task 
    'WorkManager' , ## Task-manager 
    'Checker'     , ## check of the object can be pickled/unpickled  
) 
# =============================================================================
import sys, os, time
from   itertools                    import repeat , count
from   packaging.version            import Version
from   ostap.utils.progress_bar     import progress_bar
from   ostap.parallel.task          import Task, TaskManager, keyboard_interrupt 
from   ostap.io.checker             import PickleChecker as Checker
from   ostap.core.ostap_types       import sized_types 
# =============================================================================
from   ostap.logger.logger          import getLogger
logger  = getLogger('ostap.parallel.parallel_ipyparallel')
# =============================================================================
## Try to import ipyparallel 
# =============================================================================
try : # =======================================================================
    import ipyparallel
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    ipyparallel = None
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import dill
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    dill = None
# =============================================================================
if dill : # ===================================================================
    # =========================================================================
    from ostap.parallel.dill_checker  import DillChecker as Checker
    #==========================================================================
else : # ======================================================================
    # =========================================================================
    from ostap.parallel.cloud_checker import CloudChecker as Checker 
# =============================================================================
## @class WorkManager
#  Class in charge of managing the tasks and distributing them to
#  the workers.
#  They can be local (using other cores) or remote
#  using other nodes in the local cluster """
class WorkManager(TaskManager) :
    """ Class to in charge of managing the tasks and distributing them to
    the workers. They can be local (using other cores) or remote
    using other nodes in the local cluster
    """
    def __init__( self ,
                  ncpus            = 'autodetect' , * ,
                  silent           = True         ,
                  progress         = True         ,
                  balanced         = True         ,
                  use_dill         = True         ,
                  block_size       = -1           , 
                  hyper_block_size = -1           ,                   
                  dump_dbase       = None         ,
                  dump_jobs        = 0            ,
                  dump_freq        = 0            , **kwargs ) :

        
        assert ipyparallel and Version ( ipyparallel.__version__ ) >= Version ( "8.0.0" ) , \
                "ipyparallel is not available or too old, please install ipyparallel>=8.0.0"
            
        import logging 
        config = { 'engine_timeout' : 120   ,
                   'quiet'          : True  , 
                   'create_profile' : True  , 
                   'log_level'      : logging.WARNING }
        
        if silent : 
            import logging
            logging.getLogger ( "ipyparallel" ).setLevel ( logging.WARNING )

        config.update ( kwargs )
        if 'ppservers' in config : config.pop ( 'ppservers' )
        
        ##  ipp.Cluster arguments 
        self.__balanced = True if          balanced else False 
        self.__use_dill = True if dill and use_dill else False

        ## initialize the base class 
        TaskManager.__init__  ( self             ,
                                ncpus            = ncpus      ,
                                silent           = silent     ,
                                progress         = progress   ,
                                block_size       = block_size       ,
                                hyper_block_size = hyper_block_size ,                                
                                dump_dbase       = dump_dbase ,
                                dump_jobs        = dump_jobs  ,
                                dump_freq        = dump_freq  , **config ) 

    @property
    def balanced ( self ) :
        """`balanced` : balanced processing?"""
        return self.__balanced

    @property
    def use_dill ( self ) :
        """`use_dill` : use dill for serialization?"""
        return self.__use_dill
    
    # ==================================================================================
    @property
    def config ( self ) :
        """`config` : get all configuration parameters"""
        conf = {}
        conf.update ( super().config ) 
        conf [ 'balanced' ] = self.balanced
        conf [ 'use_dill' ] = self.use_dill 
        return conf
     
    # ==================================================================================
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
        
        block_size = ( myargs.pop ( 'max_outstanding' , None ) or 
                       myargs.pop ( 'block_size'      , None ) or self.block_size )            
        if not isinstance ( block_size , int ) or block_size <= 1 : block_size = self.block_size
        
        config = dict ( max_outstanding = block_size ,
                        ordered         = ordered    )
        
        ## number of jobs 
        njobs = ( myargs.pop ( 'njobs'     , None ) or 
                  myargs.pop ( 'max_value' , None ) or
                  ( len ( jobs_args ) if isinstance ( jobs_args , sized_types ) else None ) )
        
        if not self.__balanced :
            chunk_size = myargs.pop ( 'batch_size'  , None ) or myargs.pop ( 'chunk_size' , None  )
            if not isinstance ( chunk_size , int ) or chunk_size < 1 :
                chunk_size = self.chunksize_guess ( jobs_args , njobs ) 
                logger.debug ( "`chunksize' is chosen to be %s'" % chunk_size )
            config [ 'chunksize' ] = chunk_size
                
        ## progress-bar description
        description = myargs.pop ( 'description' , "Jobs:" )
        
        done     = 0

        myargs.pop ( 'engine_timeout' , None ) 
        myargs.pop ( 'log_level'      , None ) 
        myargs.pop ( 'quiet'          , None )
        myargs.pop ( 'create_profile' , None )  
        
        modules_to_import = myargs.pop ( "imports" , [] )
        if isinstance ( modules_to_import , str ) : modules_to_import = [ modules_to_import ]

        if myargs : self.extra_arguments ( **myargs ) 

        cluster_conf = {}
        cluster_conf.update ( self.params )
        for key in ( 'max_outstanding' ,
                     'block_size'      ,
                     'batch_size'      ,
                     'chunk_size'      ,
                     'njobs'           ,
                     'max_value'       ,
                     'progress'        ,
                     'silent'          ,                     
                     'imports'         ) :  cluster_conf.pop ( key , None )

        with ipyparallel.Cluster ( n = self.ncpus , **cluster_conf  ) as client :

            client.wait_for_engines ( self.ncpus , timeout = 360 )

            if   self.use_dill : client[:].use_dill()
            
            if modules_to_import:
                command = '; '.join( f'import {mod}' for mod in modules_to_import )
                client[:].execute ( command , block = True )
                
            ## BALANCED ? 
            view = client.load_balanced_view() if self.balanced else client[:]

            # ================================================================
            try : # ==========================================================
                # ============================================================                    
                results = view.imap ( job , jobs_args  , **config )                    
                for result in progress_bar ( results                    ,
                                             max_value   = njobs        ,
                                             description = description  ,
                                             silent      = not progress ) : 
                    yield result
                    done += 1                        
                # ============================================================
            except KeyboardInterrupt : # =====================================
                # ============================================================
                logger.attention ( "%s only #%d jobs are processed" % ( keyboard_interrupt , done ) )
                # ===========================================================
                ## ABORT!! 
                try    : view.abort()
                except : pass
                ## 
                return
                # ============================================================
            except ipyparallel.RemoteError as e : # ==========================
                # ============================================================
                traceback_text = '\n'.join ( line.rstrip() for line in e.render_traceback() )
                logger.error ( 'RemoteError exeption after #%d jobs processed \n%s' %  ( done , traceback_text ) ) 
                # ============================================================
                try    : view.abort()
                except : pass
                ## 
                raise   
                # ============================================================                               
            except Exception : # =============================================
                # ============================================================
                logger.error ( 'Exception caught after #%d jobs processed' % done , exc_info = True )
                # ============================================================
                try    : view.abort()
                except : pass
                ## 
                raise   
        
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
    
    
# =============================================================================
if '__main__' == __name__ : # =================================================
    # =========================================================================
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )    
    if ipyparallel :
        logger.info    ( 'ipyparallel version is %s' % str ( ipyparallel.version_info ) )
    else           :
        logger.warning ( 'No ipyparallel is available!')
    
# =============================================================================
##                                                                      The END 
# =============================================================================
