#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/config.py
#  The basic configuration of ostap.
#  Ostap parses the following configuration files :
#   - <code>'~/.ostaprc'</code>
#   - <code>'~/.config/ostap/.ostaprc'</code>
#   - <code>'.ostaprc'</code>
#   - <code>$OSTAP_CONFIG</code>
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-05-19
# =============================================================================
"""The basic configuration of ostap
Ostap parses the following configuration files :
- ~/.ostaprc
- ~/.config/ostap/.ostaprc
- .ostaprc
- $OSTAP_CONFIG
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2019-05-19"
__all__     = (
    'config'    , ## the parsed configuration 
    )
# =============================================================================
import configparser, os, sys  

# =============================================================================
## print for configparger 
def _cp_str_ ( cp ) :
    import io 
    with io.StringIO() as o :
        config.write( o )
        return o.getvalue()

config = configparser.ConfigParser()

type(config).__str__  = _cp_str_
type(config).__repr__ = _cp_str_

## Define the major sections
config [ 'General'  ] = {
    'Quiet'   : 'False' ,
    'Verbose' : 'False' ,
    }

config [ 'Canvas'   ] = { 'Width' :  '1000' , 'Height' :  '800' } 
config [ 'Fit Draw' ] = {}
config [ 'Parallel' ] = {}

## the list of processes config files 
files_read = config.read ( [
    u'~/.ostaprc'                       ,
    u'~/.config/ostap/.ostaprc'         ,
    u'.ostaprc'                         ,
    os.environ.get ( 'OSTAP_CONFIG', '' ) ] )

# =============================================================================
## sections
general = config [ 'General' ]

quiet   = general.getboolean ( 'Quiet'  , fallback = False )
verbose = general.getboolean ( 'Verbose', fallback = False )

# =============================================================================
## section with canvas configuration
canvas  = config [ 'Canvas'    ]

# =============================================================================
## section for fit drawing options 
fit_draw = config [ 'Fit Draw' ]

# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.core.config' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
import logging
logging.disable ( ( logging.WARNING - 1 ) if quiet   else
                  ( logging.DEBUG   - 5 ) if verbose else ( logging.INFO - 1 ) )

# =============================================================================

import atexit
@atexit.register
def config_goodby () :
    import  datetime
    now = datetime.datetime.now() 
    if files_read :
        logger.info  ( 'The configuration of Ostap was read from %s' %  files_read )        
    import io 
    with io.StringIO() as o : 
        config.write( o )
        logger.verbose ( 'Ostap configuration:\n%s' % o.getvalue() )
    try :
        dump = '.ostap_config.txt'
        if os.path.exists ( dump ) : os.remove ( dump )
        with open ( dump , 'w' ) as ff :
            ff.write('#' + 78*'*' + '\n')
            ff.write('# Ostap configuration (read from %s)\n' % files_read )
            ff.write('#' + 78*'*' + '\n')                
            config.write( ff )            
            ff.write('#' + 78*'*' + '\n')
            ff.write('# Configuration saved at %s\n' % now.strftime('%c') )
            ff.write('#' + 78*'*' + '\n')            
        if os.path.exists ( dump ) and os.path.isfile ( dump ) :
            logger.info ( 'Ostap  configuration saved: %s' %  dump )
    except :
        pass
    

# =============================================================================
if '__main__' == __name__ :

    def _cp_hash_ ( cp ) : return hash ( str ( cp ) )
    type(config).__hash__ = _cp_hash_

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    cnf = '\n' + str(config)
    cnf = cnf.replace ('\n','\n# ')
    logger.info ( 'Ostap configuration is:%s' % cnf )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
