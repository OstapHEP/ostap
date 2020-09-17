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
    'get_ppservers'    , ## get list of PP-servers  
    'get_remote_conf'  , ## get PP-configuration for remote PP-server
    'ping'             , ## ping remote host
    'good_pings'       , ## get alive hosts
    'get_local_port'    , ## get local port number 
    )
# =============================================================================
from   builtins import range
# =============================================================================
from ostap.logger.logger    import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.parallel.utils' )
else                      : logger = getLogger ( __name__               ) 
# =============================================================================
## Get the PP-configuration for the remote host from the configuration file 
#  @code
#  env , script , profile = get_remote_config ( 'lxplu701.cern.ch' ) 
#  @endcode
#  The configuration information is looked in the following sections:
#  -  host-specific section "Parallel:<remote-host>"
#  -  domain specific section "Parallel:<remote-domain>"
#  -  global section "Parallel"
def get_remote_conf ( remote ) :
    """ Get the PP-configuration for the remote host form the configuration file 
    >>> env , script , profile = get_remote_config ( 'lxplu701.cern.ch' ) 
    The configuration information is looked in the following sections:
    -  host-specific section ``Parallel:<remote-host>''
    -  domain specific section ``Parallel:<remote-domain>''
    -  global section ``Parallel''
    """
    environment = ''
    script      = None
    profile     = None 
    
    import socket
    remote = socket.getfqdn ( remote ).lower()
    
    import ostap.core.config as CONFIG
    # =====================================================================
    # 1) Try to get specific configuration for the given remote host
    if ( not environment ) or ( not script ) or ( not profile ) :
        for k in CONFIG.config :
            if not k.startswith ( 'Parallel:' ) : continue
            klow   = k.lower()
            if klow[9:].strip() == remote : 
                node = CONFIG.config[ k ]
                if not environment : environment = node.get ( 'environment' , ''   )
                if not script      : script      = node.get ( 'script'      , None )
                if not profile     : profile     = node.get ( 'profile'     , None )                
                break
            
    # =====================================================================
    # 2) Try to get the domain-specific configuration
    if ( not environment ) or ( not script ) or ( not profile ) :
        for k in CONFIG.config :
            if not k.startswith ( 'Parallel:' ) : continue
            klow   = k.lower()
            domain = klow[9:].strip()
            if domain and remote.endswith ( domain ) and domain != remote :
                node = CONFIG.config[ k ]
                if not environment : environment = node.get ( 'environment' , ''   )
                if not script      : script      = node.get ( 'script'      , None )
                if not profile     : profile     = node.get ( 'profile'     , None )
                
    # =====================================================================
    # 3) Try to get global configuration
    if ( not environment ) or ( not script ) or ( not profile ) :
        import ostap.core.config as CONFIG
        if CONFIG.config.has_section  ( 'Parallel' ) :
            node = CONFIG.config [ 'Parallel' ]
            if not environment : environment = node.get ( 'environment' , ''   )
            if not script      : script      = node.get ( 'script'      , None )
            if not profile     : profile     = node.get ( 'profile'     , None )

    return environment , script , profile

# =============================================================================
## get the defined PP-servers through the following scan:
#  - the domain-specific configuration section ``Parallel:<local-domain>''
#  - the global configuration section ``Parallel''
#  - the environment variable <code>OSTAP_PPSERVERS</code>
#  @code
#  ppservers = get_ppservers () 
#  @endcode
def get_ppservers  ( local_host = '' ) :
    """Get the defined list of PP-servers through the following scan:
    - the domain-specific configuration section ``Parallel:<local-domain>''
    - the global configuration section ``Parallel''
    - the environment variable <code>OSTAP_PPSERVERS</code>
    
    >>> ppservers = get_ppservers ()
    
    """
    import socket, string 
    local_host = socket.getfqdn ( local_host ).lower()

    import ostap.core.config as CONFIG

    ppsvc1 = ()
    ws     = string.whitespace
    # =================================================================
    ## 1) try domain specific configuration 
    for k in CONFIG.config :
        if not k.startswith ( 'Parallel:' ) : continue
        klow   = k.lower()
        domain = klow[9:].strip()
        if domain and local_host.endswith ( domain ) and domain != local_host :
            node   = CONFIG.config[ k ]
            pp     = node.get( 'ppservers' , '()' )
            ppsvc1 = tuple ( i.strip ( ws ) for i in pp.split ( ',' ) if i.strip ( ws ) )  

    if not ppsvc1 and CONFIG.config.has_section  ( 'Parallel' ) :
        node = CONFIG.config [ 'Parallel' ]
        pp     = node.get( 'ppservers' , '()' )
        ppsvc1 = tuple ( i.strip ( ws ) for i in pp.split ( ',' ) if i.strip ( ws ) )  
            
    ## use the environment variables
    import os 
    ppsvc2    = tuple ( os.getenv ( 'OSTAP_PPSERVERS','').split( ',' ) )
    
    if ppsvc1 : logger.debug ( 'Get PP-servers from config          : %s' % list ( ppsvc1 ) )    
    if ppsvc2 : logger.debug ( 'Get PP-servers from OSTAP_PPSERVERS : %s' % list ( ppsvc2 ) )

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
    """Get the maximum size of jobs chunk
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
## Random number setting for parallel jobs
#  - python
#  - ROOT.gRandom
#  - ROOT.RooRandom 
#  @code
#  jobid = ...
#  random_random ( jobid ) 
#  @endcode
def random_random ( jobid , *args ) :
    """Random number setting for parallel jobs
    - python
    - ROOT.gRandom
    - ROOT.RooRandom
    >>> jobid = ...
    >>> random_random ( jobid ) 
    """
    
    import time, random, ROOT, sys , os, socket   
    ##
    random.seed ()
    ##
    jhid = os.getpid () , os.getppid() , socket.getfqdn () , jobid , os.uname () , time.time ()
    jhid = jhid , ( id ( ROOT ) , id ( sys ) , id ( random ) ) , hash ( args ) 
    jhid = hash ( jhid ) 
    ##
    if sys.version_info.major < 3 : random.jumpahead ( jhid )
    else :
        njumps = jhid % 9967
        for j in range ( njumps ) : random.uniform ( 0 , 1 )

    ## sleep a bit (up to one second) 
    time.sleep ( random.uniform ( 0.01 , 1.0 ) )
    
    ## now  initialize ROOT
    ROOT.gRandom.SetSeed ()
    
    ## ... and Roofit
    ROOT.RooRandom.randomGenerator().SetSeed()
    
    state = random.getstate() , ROOT.gRandom.GetSeed() , ROOT.RooRandom.randomGenerator().GetSeed() 
    

# =============================================================================
## ping the remote host 
def ping ( host ) :
    """Ping the host
    """
    logger.debug ( "Ping for %s" % host ) 
    import subprocess , shlex 
    command = "ping -c 1 -w 1 %s" % host 
    args    = shlex.split( command )
    try:
        subprocess.check_call ( args                     ,
                                stdout = subprocess.PIPE ,
                                stderr = subprocess.PIPE )
        return True 
    except subprocess.CalledProcessError:
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
    """Get the local port number from expressions
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

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )    
# =============================================================================
#                                                                       The END 
# =============================================================================

