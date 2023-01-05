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
    'ping'               , ## ping remote host
    'good_pings'         , ## get alive hosts
    'get_local_port'     , ## get local port number
    'pool_context'       , ## useful context for the pathos's Pools
    )
# =============================================================================
import sys
from   builtins import range
# =============================================================================
from ostap.logger.logger    import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.parallel.utils' )
else                      : logger = getLogger ( __name__               ) 
# =============================================================================

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
def random_random ( *jobid ) :
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
    jhid = jobid 
    jhid = jhid , os.getpid () , os.getppid() , os.uname()
    jhid = jhid , socket.getfqdn ()
    jhid = jhid , time.time ()
    jhid = jhid , id ( ROOT )  , id ( sys )   , id ( random )
    jhid = jhid , os.urandom ( 32 )
    jhid = hash ( jhid ) 
    ##
    random.seed ( jhid )
    ## 
    if sys.version_info.major < 3 :
        random.jumpahead ( jhid )
    ##
    njumps = jhid % 9967
    for j in range ( njumps ) :
        random.uniform ( 0 , 1 )
    ## 

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
## Context manager for Pathos pools
#  @code
#  with PoolContext ( pool ) :
#  ...   
#  @encode
class PoolContext(object) :
    """Context manager for Pathos pools
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
        self.__pool.clear ()
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
#  @encode
def pool_context  ( pool ) :
    """Context manager for Pathos pools
    >>> with pool_context ( pool ) :
    >>> ...   
    """    
    return  PoolContext ( pool )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )    

# =============================================================================
#                                                                       The END 
# =============================================================================

