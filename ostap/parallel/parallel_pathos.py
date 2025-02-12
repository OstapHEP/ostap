#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/parallel_pathos.py
#
#  Useful utilities for multiprocessing and parallel processing for Ostap
#  It is largely inspired by
#  the GaudiMP.Parallel module developed by Pere MATO Pere.Mato@cern.ch
#  The original module relies on some obsolete python modules (e.g. pyssh)
#  and it has limitations related to the pickling.
#
#  The upgraded module relies on <code>pathos</code> suite that
#  has very attractive functionality and solve pickling issues
#  @see https://github.com/uqfoundation/pathos
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-02-23
# =============================================================================
"""Useful utilities for multiprocessing and parallel processing for Ostap

Actualy it is just a little bit upgraded version of original
GaudiMP.Parallel module developed by Pere MATO Pere.Mato@cern.ch

The original module relies on some obsolete python modules (e.g. pyssh)
and it has limitations related to the pickling.

The upgraded module relies on <code>pathos</code> suite that
has very attractive functionality and solve the issues with pickling
@see https://github.com/uqfoundation/pathos

"""
# =============================================================================
__version__ = '$Revision$'
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = '2016-02-23'
__all__     = (
    'WorkManager' , ## task manager
    'Checker'     , ## check of the object can be pickled/unpickled  
)
# =============================================================================
import sys, os 
from   builtins                 import range
from   itertools                import repeat , count
# =============================================================================
from   ostap.parallel.task      import ( TaskManager   ,
                                         Task          , TaskMerger    , 
                                         Statistics    , StatMerger    ,
                                         task_executor , func_executor )
from   ostap.utils.basic        import loop_items 
from   ostap.utils.progress_bar import progress_bar
from   ostap.parallel.utils     import get_local_port  , pool_context  
# =============================================================================
if ( 3 , 3 ) <= sys.version_info  : from collections.abc import Sized
else                              : from collections     import Sized 
# =============================================================================
## CORE pathos 
# =============================================================================
import pathos.core as PC
# =============================================================================
## Logging
# =============================================================================
from ostap.logger.logger        import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.parallel.parallel_pathos' )
else                      : logger = getLogger ( __name__                         ) 
# =============================================================================
## helper function to access the underlyng <code>pp.Server</code> object
#  @attention It should not be abused! 
def get_pps ( pool ) :
    """ Helper function to access the underlying pp.Server object
    - It should not be abused! 
    """
    import pathos.parallel
    return pathos.parallel.__STATE.get ( pool._id , None )

# =============================================================================
## @class WorkManager
#  Class to in charge of managing the tasks and distributing them to
#  the workers. They can be local (using other cores) or remote
#  using other nodes in the local cluster
#  @code
#  wm1 = WorkManager ()                  ## use local
#  wm2 = WorkManager ( ppservers = ... ) ## use local and remote servers
#  wm3 = WorkManager ( ncpus = 0 , ppservers = ... ) ## use only remote servers
#  @endcode 
#  @author Pere MATO Pere.Meto@cern.ch
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class WorkManager (TaskManager) :
    """ Class to in charge of managing the tasks and distributing them to
    the workers. They can be local (using other cores) or remote
    using other nodes in the local cluster
    >>> wm1 = WorkManager ()                  ## use local
    >>> wm2 = WorkManager ( ppservers = ... ) ## use local and remote servers
    >>> wm3 = WorkManager ( ncpus = 0 , ppservers = ... ) ## use only remote servers
    """
    def __init__( self                     ,
                  ncpus     = 'autodetect' ,
                  ppservers = ()           ,
                  silent    = False        ,
                  progress  = True         , **kwargs ) :

        if not ( isinstance ( ncpus , int ) and 0 <= ncpus ) :
            from pathos.helpers import cpu_count
            ncpus = cpu_count ()
            
        ## initialize the base class 
        TaskManager.__init__ ( self, ncpus =  ncpus , silent = silent , progress = progress )
        
        from ostap.utils.cidict import cidict
        kwa = cidict ( **kwargs ) 

        self.__ppservers = ()
        self.__locals    = ()
        self.__pool      = ()
        
        import socket
        local_host = socket.getfqdn ().lower()  
 
        from   ostap.core.ostap_types import string_types 
        if isinstance ( ppservers , string_types ) and \
               ppservers.lower() in ( 'config' , 'auto' , '*' ) :
            from ostap.parallel.utils import get_workers 
            ppservers = get_workers( 'Pathos' , 'OSTAP_PPSERVERS' ,  local_host )

        ## use Parallel python if ppservers are specified or explicit flag
        if ppservers or kwa.pop ( 'PP' , False ) or kwa.pop ( 'Parallel' , False ) :

            ## remove duplicates (if any) - do not sort! 
            pps = []
            for p in ppservers :
                if p not in pps : 
                    pps.append ( p )
            ppservers = tuple ( pps )
            
            ## check local ports
            local_ports = []
            remotes     = []

            for p in ppservers :
                port = get_local_port ( p )
                if port : local_ports.append ( "localhost:%d" % port )
                else    : remotes.append     ( p ) 
                
            ## check alive remote hosts 
            from ostap.parallel.utils import good_pings 
            alive     = good_pings ( *remotes )
            if len ( alive ) != len ( remotes ) :
                dead = list ( set ( remotes ) - set ( alive ) )
                logger.warning ("Dead remote hosts: %s" % dead ) 
            remotes = alive
            
            environment = kwa.pop ( 'environment' , ''   )
            script      = kwa.pop ( 'script'      , None )
            profile     = kwa.pop ( 'profile'     , None ) 
            secret      = kwa.pop ( 'secret'      , None )
            timeout     = kwa.pop ( 'timeout'     , 7200 )
            
            if script :
                assert os.path.exists ( script ) and os.path.isfile ( script ) ,\
                       'WorkManager: no script %s is found' % script
                
            if secret is None :
                from ostap.utils.utils import gen_password 
                secret = gen_password ( 16 )

            from ostap.parallel.pptunnel import ppServer, show_tunnels  
            ppsrvs = [ ppServer ( remote                    ,
                                  environment = environment ,
                                  script      = script      ,
                                  profile     = profile     ,
                                  secret      = secret      ,
                                  timeout     = timeout     ,
                                  silent      = self.silent ) for remote in remotes ]
            

            ppbad  = [ p for p in ppsrvs if not p.pid ]
            ppgood = [ p for p in ppsrvs if     p.pid ]
            
            self.__ppservers = tuple ( ppgood ) 
            
            if ppbad :
                rs = [ p.remote_host for p in ppbad ] 
                logger.warning ('Failed to start remote ppservers at %s' %  rs ) 

            ## if remote servers are available, reduce a bit the load for local server
            ## if ncpus == 'autodetect' or ncpus == 'auto' :
            ##    self.ncpus = max  ( 0 , self.ncpus - 2 )
            
            ## some trick to setup the password.
            ## unfortunately ParallelPool interface does not allow it :-(
            if secret : 
                import pathos.parallel as PP
                _ds = PP.pp.Server.default_secret 
                PP.pp.Server.default_secret = secret 

            ## add remote connections to the list of local ports 
            for p in self.ppservers : local_ports.append ( p.local )
            self.__locals = tuple ( local_ports ) 

            from pathos.pools import ParallelPool 
            self.__pool      = ParallelPool ( ncpus = self.ncpus , servers = self.locals  )

            ## end of the trick
            PP.pp.Server.default_secret = _ds
            if not self.silent and self.ppservers : show_tunnels ( self.ppservers )
            
        else :
            
            ## from pathos.multiprocessing import ProcessPool 
            from pathos.pools import ProcessPool 
            self.__pool      = ProcessPool ( self.ncpus )

        ps = '%s' % self.pool
        ps = ps.replace( '<pool ' , '' ).replace  ('>','').replace ('servers','remotes')
        for p in self.ppservers : ps = ps.replace ( p.local , p.remote )
        if not self.silent : logger.info ( 'WorkManager is %s' % ps )

        if kwa : self.extra_arguments ( **kwa )

    @property
    def pool ( self ) :
        """``pool'' : the actual processing (ParallelPool or ProcessPool) pool"""
        return self.__pool

    @property
    def ppservers ( self ) :
        """``ppservers'' : the actual list of remote pp-servers"""
        return self.__ppservers
    
    @property
    def locals  ( self ) :
        """``locals'' : list of (local) tunnel ports"""
        return self.__locals 
    
    @property
    def remotes ( self ) :
        """``remotes'' : list of (remote) tunnel ports"""
        return tuple (  ( p.remote for p in self.ppservers ) ) 

    ## context protocol: restart the pool 
    def __enter__  ( self      ) :
        sys.stdout .flush ()
        sys.stderr .flush ()
        if self.pool : self.pool.restart ( True )                           
        return self
    
    ## context protocol: close/join/clear the pool 
    def __exit__   ( self , *_ ) :        
        if  self.pool :
            self.pool.close()
            self.pool.join  ()
            self.pool.clear ()
        sys.stdout .flush ()
        sys.stderr .flush ()

    ## def __del__(self):
    ##     self.__exit__ ()
        
    ##     del self.__pool
    ##     del self.__ppservers 
    ##     del self.__locals         

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
        
        with pool_context ( self.pool ) as pool :

            ## create and submit jobs 
            jobs = pool.uimap ( job , jobs_args )
            
            njobs = kwargs.pop ( 'njobs' , kwargs.pop ( 'max_value' , len ( jobs_args ) if isinstance ( jobs_args , Sized ) else None ) ) 
            silent = self.silent or not progress
            
            ## retrive (asynchronous) results from the jobs
            for result in progress_bar ( jobs ,
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
            
    # =========================================================================
    ## get the statistics from the parallel python
    def get_pp_stat ( self ) :
        
        smpp = StatMerger ()
        if self.pool and self.locals and get_pps ( self.pool ) : 
            pps  = get_pps ( self.pool )
            statistics = pps.get_stats().items()
            for pp , stat in statistics :
                srv = pp 
                for q in self.ppservers :
                    if q.local == srv :
                        srv = q.remote
                        break 
                s       = Statistics ( srv )
                s.njobs = stat.njobs
                s.time  = stat.time
                smpp   += s
        return smpp


# =============================================================================
DILL_COMMAND = """import sys, dill
with open('%s','rb') as f : dill.load ( f )"""
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import dill 
    from ostap.io.pickling import PickleChecker 
    # =========================================================================
    ## @class DillChecker
    #  Check if the object can be properly pickled/unpickled
    class DillChecker(PickleChecker) :
        """ Check if the object can be properly pickled/unpickled
        """
        def known ( self , *objtypes) :
            return super(Checker,self).known ( *objtypes) or \
                all ( o in self.EXTRA_TYPES for o in objtypes )
        # =====================================================================
        ## check if the object can be properly pickled/unpickled 
        def pickles ( self , *objects ) :
            """ Check of the object can be properly pickled/unpickled
            """
            return self._pickles ( *objects               ,
                                   fun_dumps = dill.dumps ,
                                   fun_loads = dill.loads ) 
        # =====================================================================
        ## Check pickling of an object across another (sub) process
        def pickles_process ( self , *objects , fast = False  ) :
            """ Check pickling of an object across another (sub)process
            """
            return self._pickles_process ( *objects                ,
                                           fun_dump = dill.dump    ,
                                           command  = DILL_COMMAND ,
                                           fast     = fast         ) ;
        # =========================================================================
        ## add new type into th elist of "known-types"
        def add ( self , ntype ) :
            """ Add new type into th elist of "known-types
            """
            if ntype in self : return 
            self.EXTRA_TYPES.add ( ntype )         
    # =========================================================================
    ## usefulname
    Checker = DillChecker 
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    from ostap.io.pickling import PickleChecker as Checker 
    # =========================================================================

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
