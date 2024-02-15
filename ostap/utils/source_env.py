#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/utils/source_env.py
#  'Source' the  environmen tfile
#  @code
#  from ostap.utuls.source_env import source_env
#  ## get the dict of all modified/new variables 
#  variables =  source_env ( 'my_env_script.sh' )
#  ## inspect the list and modify <code>os.environ</code>
#  for key in  variables :
#     if k in good_keys :
#        os.environ [ key ] = variables [  key ] 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-09-15
# =============================================================================
"""'Source' the environment file

>>> from ostap.utuls.source_env import source_env
>>> # get the dict of all modified/new variables 
>>> variables =  source_env ( 'my_env_script.sh' )
>>> # inspect the list and modify os.environ
>>> for key in  variables :
... if k in good_keys :
...    os.environ [ key ] = variables [  key ] 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2020-09-15"
__all__     = (  
    'source_env' , ## 'source' environment script 
  ) 
# =============================================================================
import sys , os , time
import subprocess, shlex 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.source_env' )
else                       : logger = getLogger( __name__ )
# =============================================================================
skip = ( '_' , 'SHLVL' ) 
# ============================================================================
def clip ( text  , lmax = 55 ) :
    if lmax < 5    : return text
    ltxt = len ( text ) 
    if ltxt < lmax : return text
    ll = ( lmax - 5 ) // 2 
    return  text [:ll] + '<...>' + text [-ll:]

# ============================================================================
#  'Source' the  environment file
#  @code
#  from ostap.utuls.source_env import source_env
#  ## get the dict of all modified/new variables 
#  variables =  source_env ( 'my_env_script.sh' )
#  ## inspect the list and modify <code>os.environ</code>
#  for key in  variables :
#     if k in good_keys :
#        os.environ [ key ] = variables [  key ] 
#  @endcode 
def source_env ( scripts , silent = False , lmax1 = 65 , lmax2 = 55 ) :
    """'Source' the environment file
    >>> from ostap.utuls.source_env import source_env
    >>> get the dict of all modified/new variables 
    >>> variables =  source_env ( 'my_env_script.sh' )
    >>> # inspect the list and modify os.environ
    >>> for key in  variables :
    ... if k in good_keys :
    ...    os.environ[key] = variables [  key ] 
    """
    from ostap.core.ostap_types import string_types

    if isinstance (  scripts , string_types ) : scripts = [ scripts ] 

    if not scripts : return 

    before = {}
    before.update ( os.environ )

    variables = {} 
    for script in scripts :
        if not os.path.exists ( script ) or not os.path.isfile ( script ) :
            logger.error ( "Invalid script '%s' for sourcing! skip it!" % script )
            
        command = 'sh -c "source %s && env"' % script 
        pargs   = shlex.split ( command )
        proc    = subprocess.Popen ( pargs , stdout = subprocess.PIPE ) 
    ## stderr = subprocess.PIPE )
    ## error = False 
    ## for line in proc.stderr :
    ##     if line : line = line [:-1]
    ##     logger.error ( "stderr: %s" % line )
    ##     error = True
        
    ## assert not error,\
    ##        "source_env: cannot source script %s using '%s'" % ( script , command )
        vars = {} 
        for line in proc.stdout:
            if line : line  = line [:-1]        
            key , _ , value = line.partition("=")
            if _ and not key in skip :
                vars [ key ] = value
        proc.communicate()
        
        variables.update ( vars ) 
    
    new_keys = set ( variables.keys() ) - set ( before.keys () )        
    modified = set()
    for key in variables :
        if key in before and variables [ key ] != before [ key ] :
            modified.add ( key )

    if not silent and ( new_keys or modified ) :
        
        import ostap.logger.table as Table

        if new_keys : 
            new_keys = list ( new_keys )
            new_keys.sort ()
            rows = [ (  "Variable" , 'Value' , '#') ] 
            for key in new_keys :
                value = variables [ key ]
                
                path  = value.split(':')
                while '' in path :  path.remove ('')
                path = len ( path )
                path = str ( path ) if 2<= path else ''
                
                value = clip ( value , lmax1 )
                row   =  key , value , path 
                rows.append ( row )
            table = Table.table ( rows , title = "New variables" , prefix = '# ')
            logger.info ('New %d variables:%s\n%s' %  ( len ( new_keys ) , new_keys , table ) )
        
        if modified :
            modified = list ( modified )
            modified.sort()
            rows = [ ( "Variable" , 'New value' , '#' , 'Old value' , '#') ]
            for key in modified :
                value1 = variables [key]
                value2 = before    [key]

                path1  = value1.split(':')
                while '' in path1 :  path1.remove ('')
                path1 = len ( path1 )
                path1 = str ( path1 ) if 2<= path1 else ''

                path2  = value2.split(':')
                while '' in path2 :  path2.remove ('')
                path2 = len ( path2 )
                path2 = str ( path2 ) if 2<= path2 else ''
                
                value1 = clip ( value1 , lmax2 )
                value2 = clip ( value2 , lmax2 )
                
                row = key , value1 , path1 , value2 , path2 
                
                rows.append ( row )
            table = Table.table ( rows , title = "Modified variables" , prefix = '# ')
            logger.info ('Modified %d variables:%s\n%s' % ( len ( modified ) , modified , table ) )
        
    result = {}
    for key in new_keys : result [ key ] = variables [ key ] 
    for key in modified : result [ key ] = variables [ key ]
    
    return result 

        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    import argparse
    parser = argparse.ArgumentParser  (
        prog       = "source_env" )
    parser.add_argument    (
        '-s' , '--silent' ,
        help    = "Silent processing [default: %(Default)s]" ,
        action  = "store_true")
    parser.add_argument    (
        '-l1', '--lmax1'  ,
        type    = int ,
        help    = "Maximal width of columns in 'new-variables' table [default: %(default)s]" , 
        default = 65  )
    parser.add_argument    (
        '-l2', '--lmax2'  ,
        type    = int ,
        help    = "Maximal width of columns in 'modified-variables' table [default: %(default)s]" ,
        default = 55  )
    parser.add_argument   (
        "scripts"     ,
        nargs   = '*' ,
        help    = "Scripts to be 'sourced'"
        )
    
    import sys
    config = parser.parse_args ( sys.argv[1:] )
    logger.debug ( "Configuration: %s" % config ) 

    if config.scripts :
        source_env ( config.scripts , config.silent , config.lmax1 , config.lmax2 ) 
    
# =============================================================================
##                                                                      The END
# =============================================================================
