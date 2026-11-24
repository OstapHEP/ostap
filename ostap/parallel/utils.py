#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/utils.py
#
#  Useful utilities for multiprocessing and parallel processing for Ostap
#  Actualy it is just a little bit upgraded version of original
#  GaudiMP.Parallel module developed by Pere MATO Pere.Mato@cern.ch
#
#  The original module relies on some obsolete python modules (e.g. pyssh)
#  and it has limitations related to the pickling.
#
#  The upgraded module relies on <code>pathos</code> suite that
#  has very attratcive functionality and solve pickling issues
#  @see https://github.com/uqfoundation/pathos
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-02-23
# =============================================================================
""" Useful utilities for multiprocessing and parallel processing for Ostap

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
    'get_remote_conf'    , ## get remote configurtaion
    'get_workers'        , ## get list of remote workers/servers 
    'ping'               , ## ping remote host
    'good_pings'         , ## get alive hosts
    'get_local_port'     , ## get local port number
    'pool_context'       , ## useful context for the pathos's Pools
    ## 
    'fix_ppsrv'          , ## small fix for pp-server
    'pp_submit'          , ## fixed "submit" method for pp
    ##
    'random_python'      , ## Re-initialize random number generation for pytohon&numpy
    'random_root'        , ## Re-initialize random number generation for ROOT&RooFit     
)
# =============================================================================
import sys, os, socket, importlib
# =============================================================================
from ostap.logger.logger    import getLogger
logger = getLogger ( 'ostap.parallel.utils' )
# =============================================================================
## Get the remove configuration for the remote host from the certain section of
#  the configuration file 
#  @code
#  env , script , profile = get_remote_config ( 'Pathos' , 'lxplus701.cern.ch' ) 
#  @endcode
#  The configuration information is looked in the following sections:
#  -  host-specific section "<SECTION >:<remote-host>"
#  -  domain specific section "<SECTION>:<remote-domain>"
#  -  global section "<SECTION>"
def get_remote_conf ( section , remote ) :
    """ Get the PP-configuration for the remote host form the configuration file 
    >>> env , script , profile = get_remote_config ( 'Pathos' , 'lxplus701.cern.ch' ) 
    The configuration information is looked in the following sections:
    -  host-specific section `<SECTION>:<remote-host>'
    -  domain specific section `<SECTION>:<remote-domain>'
    -  global section `<SECTION>'
    """
    environment = ''
    script      = None
    profile     = None 
    
    import socket
    remote = socket.getfqdn ( remote ).lower()

    nl = len ( section ) + 1 
    import ostap.core.config as CONFIG
    # =====================================================================
    # 1) Try to get specific configuration for the given remote host
    if ( not environment ) or ( not script ) or ( not profile ) :
        for k in CONFIG.config :
            if not k.startswith ( '%s:' % section ) : continue
            klow   = k.lower()
            if klow [ nl: ].strip() == remote : 
                node = CONFIG.config[ k ]
                if not environment : environment = node.get ( 'environment' , ''   )
                if not script      : script      = node.get ( 'script'      , None )
                if not profile     : profile     = node.get ( 'profile'     , None )                
                break
            
    # =====================================================================
    # 2) Try to get the domain-specific configuration
    if ( not environment ) or ( not script ) or ( not profile ) :
        for k in CONFIG.config :
            if not k.startswith ( '%s:' % section ) : continue
            klow   = k.lower()
            domain = klow [ nl: ].strip()
            if domain and remote.endswith ( domain ) and domain != remote :
                node = CONFIG.config[ k ]
                if not environment : environment = node.get ( 'environment' , ''   )
                if not script      : script      = node.get ( 'script'      , None )
                if not profile     : profile     = node.get ( 'profile'     , None )
                
    # =====================================================================
    # 3) Try to get global configuration
    if ( not environment ) or ( not script ) or ( not profile ) :
        import ostap.core.config as CONFIG
        if CONFIG.config.has_section  ( section ) :
            node = CONFIG.config [ section ]
            if not environment : environment = node.get ( 'environment' , ''   )
            if not script      : script      = node.get ( 'script'      , None )
            if not profile     : profile     = node.get ( 'profile'     , None )

    return environment , script , profile

# =============================================================================
## get the defined workers through the following scan:
#  - the domain-specific configuration section `<SECTION>:<local-domain>'
#  - the global configuration section `<SECTION>'
#  - the environment variable <code>VARIABLE</code>
#
#  Informartion is taken  form node <code>workers</code> of the section
#  and the list is entended from the envirnment variable
#  @code
#  ppservers = get_workers ( 'Pathos' , env_var = 'OSTAP_PPSERVERS' ) 
#  @endcode
def get_workers ( section , envvar = '' , local_host = '' ) :
    """ Get the defined list of workers through the following scan:
    - the domain-specific configuration section `<SECTION >:<local-domain>'
    - the global configuration section `<SECTION>'
    - the environment variable <code>envvar</code>
    
    >>> ppservers = get_workers ('Pathos' , envvar = 'OSTAP_PPSERVERS' )

    Information is taken  form node `workers` of the section
    and the list is entended from the envirnment variable

    """
    import socket, string 
    local_host = socket.getfqdn ( local_host ).lower()

    import ostap.core.config as CONFIG
    
    nl = len ( section ) + 1 
    ppsvc1 = ()
    ws     = string.whitespace
    # =================================================================
    ## 1) try domain specific configuration 
    for k in CONFIG.config :
        if not k.startswith ( '%s:' % section ) : continue
        klow   = k.lower()
        domain = klow [ nl: ].strip()
        if domain and local_host.endswith ( domain ) and domain != local_host :
            node   = CONFIG.config[ k ]
            pp     = node.get ( 'workers' , '()' )
            ppsvc1 = tuple ( i.strip ( ws ) for i in pp.split ( ',' ) if i.strip ( ws ) )  

    if not ppsvc1 and CONFIG.config.has_section  ( sections ) :
        node = CONFIG.config [ section ]
        pp     = node.get( 'workers' , '()' )
        ppsvc1 = tuple ( i.strip ( ws ) for i in pp.split ( ',' ) if i.strip ( ws ) )  
            
    ## use the environment variables
    import os 
    ppsvc2    = tuple ( os.getenv ( 'OSTAP_PPSERVERS','').split( ',' ) )
    
    if ppsvc1 : logger.debug ( 'Get workers from config          : %s' % list ( ppsvc1 ) )    
    if ppsvc2 : logger.debug ( 'Get workerd from "%s"            : %s' % ( envvar , list ( ppsvc2 ) ) ) 
    
    ppsvc = set ( ppsvc1 + ppsvc2 )
    ppsvc.discard ( '' )
    
    return tuple ( ppsvc )
                
# =============================================================================
## Get the maximum size of jobs chunk
#  for large number of parallel jobs one often gets error
#  <code>OSError 24 ("Too many open files") </code>
#  The good solution is to inclrease the limits
#  - either via <code>ulimit</code>
#  - or (better) via modification in <code>/etc/security/limits.conf</code>
#  @code
#  *               soft    nofile         65535
#  *               hard    nofile         65535 
#  @endcode
#  But this idea does not work nicely due to protection/permissions.
#  Therefore the alternative idea is to split number of parallel jobs into
#  several chunks of "reasonable" size, e.g. 20% of the soft limit
#  @code
#  maxjobs = get_max_jobs_chunk ( 1000 ) 
#  @endcode
def get_max_jobs_chunk ( jobs_chunk = None ) :
    """ Get the maximum size of jobs chunk
    for large number of parallel jobs one often gets error
    ``OSError 24 ('Too many open files')''
    The good solution is to increase the limits
    - either via ``ulimit''
    - or (better) via modification in ``/etc/security/limits.conf'', such as :
    ... *               soft    nofile         65535
    ... *               hard    nofile         65535 
    But this idea does not work nicely due to protection/permissions.
    Therefore the alternative idea is to split number of parallel jobs into
    several chunks of ``reasonable'' size, e.g. 20% of the soft limit :
    
    >>> maxjobs = get_max_jobs_chunk ( 1000 )
    
    """
    import resource
    lim = min ( *resource.getrlimit ( resource.RLIMIT_NOFILE) ) // 3 
    from ostap.core.ostap_types import integer_types 
    if isinstance ( jobs_chunk , integer_types ) and 1 <= jobs_chunk :
        return min ( lim , jobs_chunk )
    return lim

# =============================================================================
## ping the remote host 
def ping ( host ) :
    """ Ping the host
    """
    logger.debug ( "Ping for %s" % host ) 
    import subprocess , shlex 
    command = "ping -c 1 -w 1 %s" % host 
    args    = shlex.split( command )
    # =========================================================================
    try: # ====================================================================
        # =====================================================================
        subprocess.check_call ( args                     ,
                                stdout = subprocess.PIPE ,
                                stderr = subprocess.PIPE )
        return True
        # =====================================================================
    except subprocess.CalledProcessError: # ===================================
        # =====================================================================
        return False 

# =============================================================================
## get avive remote hosts (hosts with a good ping)
#  @code
#  good = good_pings ( '...' , '...' , '...' , '...' ) 
#  @endocde 
def good_pings ( *remotes ) :
    """Get alive remote hosts (hosts with a good ping)
    >>> good = good_pings ( '...' , '...' , '...' , '...' ) 
    """
    from ostap.core.ostap_types import string_types
    
    good = [] 
    for rem in remotes :

        remo = rem 
        if isinstance ( rem , string_types )  : remo = [ rem ]
        
        for remote in remo  :
            user , at , host = remote.rpartition('@')
            host , _  , port = host.partition   (':')
            if ping ( host ) : good.append ( host )
            
    return tuple ( good ) 

# =============================================================================
_patterns = []

# ============================================================================
## get the local port number from expressions:
#  valid expressions (leading and trailing spacdes are ignored)
#   - positive integer number
#   - ' positive-integer-number '
#   - ' localhost:positive-integer-number ' ##  case insensitive 
#   - ' positive-integer-number:localhost ' ##  case insensitive 
def get_local_port ( expression ) :    
    """ Get the local port number from expressions
    Valid expressions (leading and trailing spacdes are ignored)
    - positive integer-number
    - `'  positive-integer-number  '`
    - `'  localhost:positive-integer-number  '` ##  case insensitive 
    - `'  positive-integer-number:localhost  '` ##  case insensitive 

    """
    from ostap.core.ostap_types import string_types, integer_types
    
    if isinstance ( expression , integer_types ) and 0 < expression :
        return expression

    if not _patterns :
        
        import re
        _patterns.append ( re.compile ( r'\A(\d+)\Z'            , re.I ) ) 
        _patterns.append ( re.compile ( r'\Alocalhost:(\d+)\Z'  , re.I ) ) 
        _patterns.append ( re.compile ( r'\A\A(\d+):localhost\Z', re.I ) ) 
        
    if isinstance ( expression , string_types ) :

        expr = expression.strip().lower()
        for pattern in _patterns :

            match = pattern.match ( expr )
            if not match : continue
            
            port = int ( match.group ( 1 ) )
            if 0 < port : return port            ## RETURN 

    return None 
                    
# =============================================================================
## Context manager for Pathos pools
#  @code
#  with PoolContext ( pool ) :
#  ...   
#  @endcode
class PoolContext(object) :
    """ Context manager for Pathos pools
    >>> with PoolContext ( pool ) :
    >>> ...   
    """
    def __init__ ( self , pool ) :
        self.__pool = pool
        
    def __enter__ ( self )   :
        sys.stdout .flush ()
        sys.stderr .flush ()
        self.__pool.restart ( True )
        return self.__pool
    
    def __exit__  ( self , *_) :
        self.__pool.close ()
        self.__pool.join  ()
        self.__pool.clear ()   ## COMMENTED OUT 
        sys.stdout .flush ()
        sys.stderr .flush ()
        
    @property
    def pool  ( self ) :
        """``pool'' the actual Pathos pool"""
        return self.__pool 

# =============================================================================
## Context manager for Pathos pools
#  @code
#  with pool_context ( pool ) :
#  ...   
#  @endcode
def pool_context  ( pool ) :
    """ Context manager for Pathos pools
    >>> with pool_context ( pool ) :
    >>> ...   
    """    
    return  PoolContext ( pool )

# =============================================================================
## Fixed "submit" method for pp 
# =============================================================================
def pp_submit ( ppsrv , func, *args , **kwargs ) :
    """ Fixed 'submit'-method for pp 
    """
    from  ostap.utils.core import  ismethod, isfunction, typename
    ##
    if   ismethod   ( func ) : the_func = func 
    elif isfunction ( func ) : the_func = func 
    elif icallable  ( func ) : the_func = func.__call__
    else :
        raise TypeError ( 'Invalid type %s' % typename ( func ) )
    ## 
    return ppsrv.submit ( the_func , *args , **kwargs )

# =============================================================================
## some "fix" for "submit" method for pp  
def fix_ppsrv ( ppsrv ) :
    """ Some "fix" for "submit" method for pp  
    """ 
    if not isinstance ( ppsrv , type ) : ppsrv = type ( ppsrv )
    
    if not hasattr ( ppsrv , '_old_submit_' ) :
        ppsrv._old_submit_ = ppsrv.submit
        
        from  ostap.utils.core import  ismethod, isfunction, typename
        
        # ======================================================================
        def _pp_new_submit_ ( self , func , *args , **kwargs ) :
            """ Fixed "submit" method
            """            
            if   ismethod   ( func ) : the_func = func 
            elif isfunction ( func ) : the_func = func 
            elif callable   ( func ) : the_func = func.__call__ 
            else :
                raise TypeError ( 'Invalid type %s' % typename ( func ) )
            ##
            return self._old_submit_ ( the_func , *args , **kwargs )
        
        ppsrv._new_submit_ = _pp_new_submit_
        ppsrv.submit       = _pp_new_submit_



# =============================================================================
## Worker initialization function.
#
#  Imports specified modules directly into the worker's global namespace
#  dynamically without requiring the 'global' keyword.
#
#  @param modules : str or list of str
#       Module name or list of module names to import on the worker (e.g.,
def init_worker_modules ( modules ):
  """ Worker initialization function.
  
  Imports specified modules directly into the worker's global namespace
  dynamically without requiring the 'global' keyword.
  
  Parameters
  ----------
  modules : str or list of str
      Module name or list of module names to import on the worker (e.g.,
      'array' or ['array', 'numpy']).
  """
  if not modules: return
  elif  isinstance ( modules, str ) : modules = [modules]
  
  # Retrieve the module-level globals dictionary for the executing worker
  gl = globals()

  import importlib 
  for mod in modules:
    # Handle submodules (e.g., 'os.path' -> assigns variable name 'os')
    base_name = mod.split(".")[0]
    gl [ base_name ] = importlib.import_module ( mod )


# ===============================================================================
## "Unordered map iterator" over the list of (jobid,job) pairs
def uimap  ( jobs ) :
    """ Unorderd-map-iterator over the list of  (jobid, job) pairs
    """
    lst = list ( jobs ) 
    while lst :
        for i, job_pair in enumerate ( lst  ) :
            jobid , job = job_pair
            if job.finished :
                lst.pop ( i ) 
                yield  jobid, job 
                break

# =============================================================================
SEED_MAX       = 2**63
SEED_MAX_NUMPY = 2**32 
SEED_MAX_ROOT  = 2**32 - 1 
# =============================================================================
## Setup random numbers for python&numpy
#  - if seed if positive integer number, if it used for python and numpy
#  
#  @see random.seed
#  @see numpy.random.set_seed
#
#  @param seed   (INPUT) input seed
#  @param salt   (INPUT) `salt`   for generation of proper seed
#  @param pepper (INPUT) `pepper` for generation of proper seed
#  @return the actual seed used to initializer python&numpy
def random_python ( seed = None , *salt , **pepper ) :
    """ Setup random numbers for python&numpy
    - if `seed` is positive integer number, it it used for python and numpy,
    otherwise `seed` if generated via `hash` from
    
    - `os.urandom`
    - `salt`
    - `pepper`
    - `sys.environ`

    
    param `seed`   (INPUT) input seed
    param `salt`   (INPUT) `salt`   for generation of proper seed
    param `pepper` (INPUT) `pepper` for generation of proper seed
    
    - see `random.seed`
    - see `numpy.random.seed`
    
    return the actual seed used to initializer python&numpy 
    """

    if not isinstance ( seed , int ) or not 1 <= seed <= SEED_MAX :

        ## use system random 
        seed = seed , os.urandom ( 16 )
        seed = hash ( seed ) 

        from   ostap.utils.basic          import main_process 
        
        ## some os-information *AND* salt 
        seed = ( seed , os.getpid () , os.getppid () , main_process() ) + salt
        seed = hash ( seed )
        
        for k , v in pepper.items () :
            seed = seed , k , v
            seed = hash ( seed )
            
        for k , v in os.environ.items () :
            seed = seed , k , v
            seed = hash ( seed )

        seed = seed % SEED_MAX + 1
    
    import random     
    random.seed ( seed )
    
    import numpy 
    numpy.random.seed ( seed % SEED_MAX_NUMPY )

    return seed 

# =============================================================================
## Setup random numbers for ROOT&RooFit 
#  - if seed is positive integer number it is used, otherwise seed is generated 
#  
#  @see random.set_seed
#  @see numpy.random.set_seed
#
#  @param seed   (INPUT) input seed
#  @param salt   (INPUT) `salt`   for generation of proper seed
#  @param pepper (INPUT) `pepper` for generation of proper seed
#  @return the actual seed used to initializer python&numpy
def random_root ( seed = None , *salt , **pepper ) :
    """" Setup random numbers for ROOT&RooFit 
     - if seed is positive integer number it is used, otherwise seed is  generated 
    
    param `seed`   (INPUT) input seed
    param `salt`   (INPUT) `salt`   for generation of proper seed
    param `pepper` (INPUT) `pepper` for generation of proper seed
    
    return the (long) seed used to initialize ROOT & RooFit 
    """

    long_seed = random_python ( seed , *salt , **pepper )

    ## clip it 
    short_seed = long_seed % SEED_MAX_ROOT + 1 

    import ROOT 
    ROOT.gRandom                     .SetSeed ( short_seed ) 
    ROOT.RooRandom.randomGenerator() .SetSeed ( short_seed )
    
    return long_seed 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )    

# =============================================================================
#                                                                       The END 
# =============================================================================

