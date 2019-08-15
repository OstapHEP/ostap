#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/mp_pathos.py
#
#  Useful utilities for multiprocessing and parallel processing for Ostap
#  Actualy it is just a little bit upgraded version of original
#  GaudiMP.Parallel module developed by Pere MATO Pere.Mato@cern.ch
#
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
# ========================================_====================================
__version__ = '$Revision$'
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = '2016-02-23'
__all__     = (
    'WorkManager' , ## task manager
    'ppServer'    , ## helper class to start remote ppserver (for parallel python)
    )
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.paralllel.mp_pathos' )
else                      : logger = getLogger ( __name__                    ) 
# =============================================================================
import sys, os 
from   builtins         import range
# =============================================================================
# PATHOS components 
# =============================================================================
import dill 
import ppft 
import pathos.core              as PC
import pathos.secure            as PS
import pathos.multiprocessing
import pathos.parallel
# =============================================================================
## helper function to access the underlyng <code>pp.Server</code> object
#  @attention It should not be abused! 
def get_pps ( pool ) :
    """Helper function to access the underlying pp.Server object
    - It should not be abused! 
    """
    return pathos.parallel.__STATE.get ( pool._id , None )

# =============================================================================
from ostap.parallel.task       import ( Task          , TaskMerger    , 
                                        Statistics    , StatMerger    ,
                                        task_executor , func_executor ) 
# =============================================================================
from ostap.utils.utils import get_open_fds

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
class WorkManager ( object ) :
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
                  silent    = False        , **kwargs ) :

        if   isinstance ( ncpus , int ) and 0 <= ncpus : self.ncpus = ncpus
        else :
            from pathos.helpers import cpu_count
            self.ncpus = cpu_count ()
            
        self.__ppservers = ()
        self.__locals    = ()

        import socket
        local_host = socket.getfqdn ().lower()  
 
        from ostap.core.ostap_types    import string_types 
        if isinstance ( ppservers , string_types ) and ppservers.lower() in ( 'config' , 'auto' ) : 
            from ostap.parallel.utils import get_ppservers
            ppservers = get_ppservers ( local_host )

        if ppservers :

            ## remove duplicates (if any) 
            pps = []
            for p in ppservers :
                if p not in pps : 
                    pps.append ( p ) 
            ppservers = tuple ( pps ) 
            
            environment = kwargs.pop ( 'environment' , ''   )
            script      = kwargs.pop ( 'script'      , None )
            profile     = kwargs.pop ( 'profile'     , None ) 
            secret      = kwargs.pop ( 'secret'      , ''   )
            timeout     = kwargs.pop ( 'timeout'     , 7200 )
            
            if script :
                assert os.path.exists ( script ) and os.path.isfile ( script ) ,\
                       'WorkManager: no script %s is found' % script
                
            if not secret :
                from ostap.utils.utils import gen_password 
                secret = gen_password ( 16 )
                
            self.__ppservers = [
                ppServer ( remote                    ,
                           environment = environment ,
                           script      = script      ,
                           profile     = profile     ,
                           secret      = secret      ,
                           timeout     = timeout     ) for remote in ppservers ]

            ## if remote servers are available, reduce a bit the load for local server
            ## if ncpus == 'autodetect' or ncpus == 'auto' :
            ##    self.ncpus = max  ( 0 , self.ncpus - 2 )

            ## some trick to setup the password.
            ## unfortnuately  ParallelPool interface does not allow it :-( 
            import pathos.parallel as PP
            _ds = PP.pp.Server.default_secret 
            PP.pp.Server.default_secret = secret 
            
            from pathos.pools import ParallelPool 
            self.__pool      = ParallelPool ( ncpus = self.ncpus , servers = self.locals  )

            PP.pp.Server.default_secret = _ds 
                                    
        else :
            
            ## from pathos.multiprocessing import ProcessPool 
            from pathos.pools import ProcessPool 
            self.__pool      = ProcessPool ( self.ncpus )

        ps = '%s' % self.pool
        ps = ps.replace( '<pool ' , '' ).replace  ('>','').replace ('servers','remotes')
        for p in self.ppservers : ps = ps.replace ( p.local , p.remote )
        logger.info ( 'WorkManager is %s' % ps )
        
        self.stats  = {}
        self.silent = True if silent else  False 

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
        return tuple (  ( p.local for p in self.ppservers ) )
    
    @property
    def remotes ( self ) :
        """``remotes'' : list of (remote) tunnel ports"""
        return tuple (  ( p.remote for p in self.ppservers ) ) 

    def __enter__  ( self      ) :
        if self.pool : self.pool.restart ( True )                           
        return self
    
    def __exit__   ( self , *_ ) :        
        if  self.pool :
            self.pool.close()
            del self.__pool
            self.__pool = None 

    def __del__(self):
        self.__exit__ ()
        
        del self.__pool
        del self.__ppservers 
        del self.__locals 
        
    # =========================================================================
    ## process callable object or Task :
    #  - process <code>Task</code>:
    #  @code
    #  class MyTask(Task) : ....
    #  my_task = MyTask ( ... ) 
    #  wm = WorkManager ( ... )
    #  result = wm.process ( my_task , items )
    #  @endcode
    #  -  process callable object 
    #  @code
    #  def my_fun ( x ) :
    #      return x**2 
    #  wm = WorkManager ( ... )
    #  items = range ( 10 )
    #  ## get list of squares as a result 
    #  result1 =  wm.process ( my_fun , items , merger = TaskMerger ( lambda  a,b : a+[b] , init = [] ) )
    #  ## get sum of them 
    #  result2 =  wm.process ( my_fun , items , merger = TaskMerger () )    
    #  @endcode
    def process ( self , task , items , **kwargs ) :
        """Process callable object or Task :
        
        - process Task
        
        >>> class MyTask(Task) : ....
        >>> my_task = MyTask ( ... ) 
        >>> wm = WorkManager ( ... )
        >>> result = wm.process ( my_task , items )
        
        -  process callable object
        
        >>> def my_fun ( x ) : return x**2 
        >>> wm = WorkManager ( ... )
        >>> items = range ( 10 )
        >>> result1 =  wm.process ( my_fun , items , merger = TaskMerger ( lambda  a,b : a+[b] , init = [] ) )
        >>> result2 =  wm.process ( my_fun , items , merger = TaskMerger () )    
        
        """
        
        from more_itertools import chunked
        
        job_chunk = kwargs.pop ( 'jobs_chunk', 10000 )
        chunks    = list ( chunked ( items , job_chunk ) )
        chunks.reverse()

        sys.stdout.flush ()
        sys.stderr.flush ()
        
        ## process chunks of data
        result = self.__process ( task , chunks , **kwargs )
        
        sys.stdout.flush ()
        sys.stderr.flush ()

        return result 
        
    # ===================================================================================
    ## helper internal method to process the task with chunks of data 
    def __process ( self , task , chunks , **kwargs ) :
        """Helper internal method to process the task with chunks of data 
        """

        if isinstance ( task , Task ) :
            kwargs.pop ( 'merger' , None ) 
            return self.__process_task ( task , chunks , **kwargs ) 

        ## mergers for statistics 
        merged_stat    = StatMerger ()
        merged_stat_pp = StatMerger ()
        merger         = kwargs.pop ( 'merger' , TaskMerger ( ) ) 

        njobs = sum  ( len(c) for c in chunks ) 
        from ostap.utils.progress_bar import ProgressBar
        with ProgressBar ( max_value = njobs , silent = self.silent ) as bar :

            while chunks :
                
                chunk = chunks.pop() 
                
                from itertools      import repeat
                jobs_args = zip ( repeat ( task ) , chunk )

                self.pool.restart ( True )
                jobs      = self.pool.uimap ( func_executor , jobs_args  )
                del jobs_args 
                
                for result , stat in jobs :
                    bar         += 1
                    merged_stat += stat
                    merger      += result 

                    del result
                    del stat 
                    
                merged_stat_pp += self.get_pp_stat()  
                self.pool.close ()
                self.pool.join  ()

        ## finalize task 
        what.finalize () 
        self.print_statistics ( merged_stat_pp , merged_stat )
        ## 
        return merger.results 

    # ===================================================================================
    ## helper internal method to process the task with chunks of data 
    def __process_task  ( self , task , chunks , **kwargs ) :
        """Helper internal method to process the task with chunks of data 
        """
        assert isinstance ( task , Task ), 'Invalid task type  %s' % type (  task ) 
        
        ## inialize the task
        task.initialize_local ()
        
        ## mergers for statistics 
        merged_stat    = StatMerger ()
        merged_stat_pp = StatMerger ()

        njobs = sum  ( len(c) for c in chunks ) 
        from ostap.utils.progress_bar import ProgressBar
        with ProgressBar ( max_value = njobs , silent = self.silent ) as bar :

            while chunks :
                
                chunk = chunks.pop() 
                
                from itertools      import repeat
                jobs_args = zip ( repeat ( task ) , chunk )

                self.pool.restart ( True )
                jobs      = self.pool.uimap ( task_executor , jobs_args  )
                del jobs_args 
                
                for result , stat in jobs :
                    bar         += 1
                    merged_stat += stat
                    task.merge_results ( result )

                    del result
                    del stat 
                    
                merged_stat_pp += self.get_pp_stat()
                self.pool.close ()
                self.pool.join  ()

        task.finalize () 
        self.print_statistics ( merged_stat_pp , merged_stat )
        ## 
        return task.results ()
    
    # =========================================================================
    ## get the statistics from the paralell python
    def get_pp_stat ( self ) :
        smpp = StatMerger ()
        if self.pool and self.ppservers and get_pps ( self.pool ) : 
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
    
    # =========================================================================
    ## print the job execution statistics 
    def print_statistics ( self , stat_pp , stat_loc ) :
        """Print the job execution statistics 
        """        
        if self.silent : return

        if stat_pp.njobs == stat_loc.njobs : 
            stat_pp .print_stats ( 'pp-' )
        else : 
            stat_loc.print_stats ( 'qq-' )
                
            
# =============================================================================
## @class ppServer
#  Helper class that starts <code>ppserver</code> on remote site
#  and makes SSH tunnel
#  @attention Password-less SSH-connection is required!
#  @code
#  pps = ppServer('lxplus009.cern.ch')
#  print pps.local, pps.remote
#  del pps
#  @endcode
#  ... or as context manager :
#  @code
#  with ppServer('lxplus009.cern.ch') as pps : 
#      print pps.local, pps.remote
#      ... do something here 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class ppServer(object) :
    """ Helper class that starts ppserver on remote site and makes SSH tunnel
    - attention:  password-less SSH-connection is required!
    
    >>> pps = ppServer('lxplus009.cern.ch')
    >>> print pps.local, pps.remote
    >>> del pps
    
    ... or as context manager :

    >>> with ppServer('lxplus009.cern.ch') as pps : 
    ...  print pps.local, pps.remote
    ...  do something here 
    """
    def __init__ ( self               ,
                   remote             ,
                   environment = ''   ,
                   script      = None ,
                   profile     = None ,
                   secret      = None ,
                   timeout     =   -1 ) :
        
        ## split user@remote
        remote_user , at , remote = remote.rpartition('@')
        
        import socket        
        remote = socket.getfqdn ( remote ).lower()  

        ## ## merge user@remote 
        ## if remote_user and at : 
        ##     import getpass
        ##     local_user = getpass.getuser()
        ##     if remote_user.lower() != local_user.lower() :
        ##         remote = ''.join ( ( remote_user, at , remote ) )
                
        remote = ''.join ( ( remote_user, at , remote ) )
                
        self.__remote_host = remote 
        
        ## SSH tunnel 
        self.__tunnel       = PS.Tunnel  ( 'SSH-tunnel %s' %  self.remote_host  )
        self.tunnel.connect ( self.remote_host )

        logger.debug ('Create SSH tunnel: %s' %  self.tunnel )

        self.__local_port   = self.tunnel._lport 
        self.__remote_port  = self.tunnel._rport 
                                  
        ## SSH session
        self.__session      = PS.Pipe ( 'SSH-session' , host = self.remote_host )  

        logger.debug ('Create SSH session: %s'%  self.session )
        
        import logging
        verbose = logging.root.manager.disable <= logging.DEBUG 
        self.session.verbose = verbose
        self.tunnel .verbose = verbose 
        del logging

        # ==================================================================
        # Get configuration information for the given remote host
        if ( not environment ) or ( not script ) or ( not profile ) :
            from ostap.parallel.utils import get_remote_conf 
            e , s , p = get_remote_conf ( remote ) 
            if e and not environment : environment = e
            if s and not script      : script      = s
            if p and not profile     : profile     = p 

        if script :
            if os.path.exists ( script ) and os.path.isfile ( script ) : pass
            else :
                logger.error ("The script %s does not exist, skip it!" % script )
                script = None 
                                  
        ## temporary directory on remote host 
        self.__tmpdir = None 
        if environment or script or profile : 
            self.session  ( command =  'mktemp -d -t pathos-XXXXXXXXXX' )
            self.session.launch()
            r = self.session.response()
            if r and 1 < len ( r ) : 
                self.__tmpdir = r [:-1] 
                logger.verbose ('Created remote temporary directory %s:%s' % ( self.remote_host , self.__tmpdir ) ) 
            else :
                logger.error   ('Cannot create remote temporary directory at %s' % self.remote_host ) 
                
        self.__environment = None
        if  environment :
            import ostap.utils.cleanup as CU
            tmpfile = CU.CleanUp.tempfile ( prefix = 'env_' , suffix = '.sh' )
            with open ( tmpfile , 'w' ) as tf :
                tf.write( environment )
                tf.write('\n')
            copier      = PS.Copier ( 'SSH-copier' )
            destination = "%s:%s" % ( self.__remote_host , self.__tmpdir )
            copier ( source = tmpfile , destination = destination )
            copier.launch()
            r = copier.response()
            if r : logger.error ('SCP: response from %s : %s' % ( copier , r ) ) 
            del copier
            self.__environment = os.path.join ( self.__tmpdir , os.path.basename ( tmpfile ) )
            logger.verbose ('SCP: environment %s:%s is copied' % ( self.remote_host , self.__environment ) )
            
        self.__script = None 
        if  script and os.path.exists ( script ) and os.path.isfile ( script ) :
            copier      = PS.Copier ( 'SSH-copier' )
            destination = "%s:%s" % ( self.__remote_host , self.__remote_tmpdir )
            copier( source = script , destination = destination )
            copier.launch()
            r =  copier.response()
            if r : logger.error ('SPC: response from %s : %s' % ( copier , r ) ) 
            del copier
            self.__script = os.path.join ( self.__tmpdir , os.path.basename ( script ) )  
            logger.verbose ('SCP: script      %s:%s is copied' %  ( self.remote_host , self.__script ) )

        self.__profile = None 
        if profile :
            self.session  ( command =  '[ -f %s ] && echo "1" || echo "0" ' %  profile  )
            self.session.launch()
            r = self.session.response()
            try :
                if int ( r ) : self.__profile = profile
            except :
                pass
            if self.__profile :
                logger.verbose ("Profile %s:%s is found"     %  ( self.remote_host , profile ) )
            else : 
                logger.warning ("Profile %s:%s is NOT found" %  ( self.remote_host , profile ) )

        commands = []
        if self.__profile     : commands.append ( ' source %s ' % profile            )
        if self.__environment : commands.append ( ' source %s ' % self.__environment )
        if self.__script      : commands.append ( ' source %s ' % self.__script      )
            
        ## pp-server itself:
        if   secret and 1 < timeout :
            commands.append ( 'ppserver -p %-7d -s %s -t %d'           % ( self.remote_port , secret , timeout ) )
            pattern = '[P,p]ython *.*ppserver *-p *%d *-s *%s *-t *%d' % ( self.remote_port , secret , timeout )
        elif secret :
            commands.append ( 'ppserver -p %-7d -s %s'                 % ( self.remote_port , secret ) ) 
            pattern = '[P,p]ython *.*ppserver *-p *%d *-s *%s'         % ( self.remote_port , secret )
        elif 1 < timeout : 
            commands.append ( 'ppserver -p %-7d -t %d'                 % ( self.remote_port , timeout ) ) 
            pattern = '[P,p]ython *.*ppserver *-p *%d *-t *%d'         % ( self.remote_port , timeout )
        else             : 
            commands.append ( 'ppserver -p %-7d'                       %   self.remote_port  ) 
            pattern = '[P,p]ython *.*ppserver *-p *%d'                 %   self.remote_port 

            
        command = ' && '.join ( commands ) 
        self.session ( command = command , background = True , options = '-f' ) 
        self.session.launch()
        r = self.session.response()
        if r : logger.error ('SERVER:response from %s : %s' % ( self.session , r ) )


        self.__pid = None
        for i in range ( 15 ) :
            try :
                import time 
                time.sleep ( 1 )
                self.__pid = PC.getpid ( pattern , self.remote_host )
                logger.verbose ('PID for remote ppserver is %s:%d' % ( self.remote_host , self.__pid ) )
            except OSError : pass
            if self.__pid  : break 
        else :
            logger.warning ('Cannot acquire PID for remote ppserver at %s:%s' % ( self.remote_host , self.remote_port ) )

        sys.stdout.flush () ## needed? 
        sys.stdin .flush () ## needed ?

        logger.debug   ( "%s" % self )
        logger.verbose ( command     )

    def start ( self ) : return self
    def end   ( self ) :

        sys.stdout.flush () ## needed ?
        sys.stdin .flush () ## needed ? 
        
        import logging
        verbose = logging.root.manager.disable <= logging.DEBUG 
        del logging
        
        if self.session : self.session.verbose = verbose 
        if self.tunnel  : self.tunnel .verbose = verbose 
        
        if self.__pid :
            try  :
                if self.session :  # and verbose :
                    command = 'kill -s SIGUSR1 %d ' % self.__pid 
                    self.session ( command = command )
                    self.session.launch   ()
                    r = self.session.response()
                    if r : logger.info ( '%s : %s : %s' %  ( self.remote_host , command , r ) )
                    
                logger.debug ( 'Kill remote process %s:%s' % ( self.remote_host , self.__pid ) )
                PC.kill ( self.__pid , self.remote_host )
            except OSError :
                pass
        self.__pid = None

        ## remove unnesessary files 
        if self.session :
            commands = []
            if self.__environment : commands.append ( 'rm -f %s' %  self.__environment )
            if self.__script      : commands.append ( 'rm -f %s' %  self.__script      )
            if self.__tmpdir      : commands.append ( 'rmdir --ignore-fail-on-non-empty %s' %  self.__tmpdir ) 
            if commands :
                command = ' ; '.join ( commands )
                self.session ( command = command ) 
                self.session.launch()
                r = self.session.response () 
                if r : logger.verbose ( "Response: %s" % r )                
            
        if self.session :
            logger.debug ( 'Kill session      %s' % self.session ) 
            self.session.kill      () 

        if self.tunnel  :
            logger.debug ( 'Disconnect tunnel %s' % self.tunnel  ) 
            self.tunnel.disconnect ()

        self.__tunnel  = None
        self.__session = None

    ## CONTEXT MANAGER: empty
    def __enter__ ( self      ) : return self.start ()

    ## CONTEX MANAGER: exit 
    def __exit__  ( self, *_  ) : return self.end   () 

    ## delete, close and cleanup 
    def __del__   ( self )      :        self.end   ()

    # =========================================================================
    ## Readable printout: representation of ppServer/Tunnel
    #  @code
    #  s = ppServer ( ... )
    #  print(s)
    #  @endcode 
    def __repr__  ( self ) :
        """ Readable printout: Representation of ppServer/Tunnel
        >>> s = ppServer ( ... )
        >>> print(s) 
        """
        return "ppServer(%d:%s:%d,pid=%s)"  % ( self.local_port , self.remote_host , self.remote_port, self.pid )
    __str__ = __repr__
    
    @property
    def pid  ( self ) :
        """``pid''  : PID for remote ppserver process"""
        return self.__pid
    
    @property
    def script ( self ) :
        """``script'' : the shell script to be sourced to get the correct  environment"""
        return self.__script

    @property
    def environment ( self ) :
        """``environment'' : the environment script to be sourced to get the correct environment"""
        return self.__environment
    
    @property
    def remote_host      ( self ):
        """``remote_host'' : remote host"""
        return self.__remote_host

    @property
    def local_port  ( self ) :
        """``local_port'' : local port for SSH-tunnel"""
        return self.__local_port
    
    @property
    def remote_port ( self ) :
        """``remote_port'' : remote port for SSH-tunnel"""
        return self.__remote_port
        
    @property
    def tunnel   ( self ) :
        """``tunnel'' : SSH tunnel   local_port <---> remote_host:remote_port"""
        return self.__tunnel 
    
    @property
    def session ( self ) :
        """``tunnel'' : SSH session"""
        return self.__session 
    
    @property
    def local ( self ) :
        """``local'' : local server/tunnel address"""
        return "localhost:%s" % self.__local_port
    
    @property
    def remote ( self ) :
        """``remote'' : remote server/tunnel address"""
        return "%s:%s" % ( self.__remote_host , self.__remote_port )


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
