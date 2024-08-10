#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/config.py
#  The basic configuration of ostap.
#  Ostap parses the following configuration files :
#  - <code>$OSTAPDIR/.ostaprc</code>
#  - <code>$HOME/.ostaprc</code>
#  - <code>~/.ostaprc</code>
#  - <code>- ~/.config/ostap/.ostaprc<.code>
#  - <code>- $HOME/.config/ostap/.ostaprc<.code>
#  - <code>.ostaprc</code>
#  - <code>$OSTAP_CONFIG</code>
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-05-19
# =============================================================================
"""The basic configuration of ostap
Ostap parses the following configuration files :
- $OSTAPDIR/.ostaprc
- $HOME/.ostaprc
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
import ostap.core.default_config as     default_config 
from   ostap.utils.basic         import get_env        as ostap_getenv 
# =============================================================================
## print for configparger 
def _cp_str_ ( cp ) :
    """print for configparger"""
    import io 
    with io.StringIO() as o :
        config.write( o )
        return o.getvalue()

config = configparser.ConfigParser()
type(config).__str__  = _cp_str_
type(config).__repr__ = _cp_str_

## Define the major sections
config [ 'General'  ] = {
    'Quiet'              : str ( default_config.quiet   ) ,
    'Verbose'            : str ( default_config.verbose ) ,
    'Parallel'           : 'PATHOS'                       ,
    'WebDisplay'         : default_config.web             , ## ROOT.TROOT.Set/Get WebDisplay                        
}

## generic TCanvas configuration
config [ 'Canvas'      ] = { 'Width'       :  '1000' , 'Height'       :  '800' , 
                             'MarginTop'   : '0.05'  , 'MarginBottom' : '0.12' ,
                             'MarginRight' : '0.07'  , 'MarginLeft'   : '0.12' }

config [ 'Fit Draw'    ] = {} ## RooFit plotting configuration
config [ 'Tables'      ] = {} ## configuration for Tables 
config [ 'RooFit'      ] = {} ## RooFit configuration

config [ 'Pathos'      ] = {} ## PATHOS configuration  
config [ 'IPyparallel' ] = {} ## ipyparallel configuration 

## the list of processes config files
config_files = default_config.config_files + tuple ( ostap_getenv ( 'OSTAP_CONFIG', '' , True ).split( os.pathsep ) )
the_files    = [] 
for f in config_files :
    ff = f
    for i in range ( 5 ) :
        ff = os.path.expandvars ( ff )
        ff = os.path.expanduser ( ff )
    if not os.path.exists ( ff )                 : continue
    if not os.path.isfile ( ff )                 : continue
    if [ fn for fn in the_files if fn[0] == ff ] : continue
    the_files.append ( ( ff , f ) )
the_files = tuple ( the_files ) 
config_files = tuple ( f [ 0 ] for f in the_files ) 

## read the files 
files_read = config.read ( config_files )

# =============================================================================
## sections
general = config [ 'General' ]


quiet   = general.getboolean ( 'Quiet'  , fallback = False )
verbose = general.getboolean ( 'Verbose', fallback = False )


# =============================================================================
## section with canvas configuration
canvas  = config [ 'Canvas'    ]

# =============================================================================
## section for RooFit 
roofit   = config [ 'RooFit'   ]

# =============================================================================
## section for fit drawing options 
fit_draw = config [ 'Fit Draw' ]

# =============================================================================
## section for Tables  
tables   = config [ 'Tables'   ]

# =============================================================================
## Redefine webdidpay fron environment variable if set 
general [ 'WebDisplay' ] = ostap_getenv ( 'OSTAP_DISPLAY' , general [ 'WebDisplay' ] )

# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.core.config' )
else                       : logger = getLogger ( __name__  )
# =============================================================================
import logging
logging.disable ( ( logging.WARNING - 1 ) if quiet   else
                  ( logging.DEBUG   - 5 ) if verbose else ( logging.INFO - 1 ) )

# =============================================================================
## check (and remove) obsolete sections 
remove = set() 
for k in config :
    if k.startswith ( 'Parallel:' ) :
        logger.error ( 'Section "%s" is obsolete! Switch to specific "Pathos:..." or "IPyparallel:..."' % k )
        remove.add ( k ) 
    elif k.startswith ( 'Parallel' ) :
        logger.error ( 'Section "%s" is obsolete! Switch to specific "Pathos" or "IPyparallel"' % k )
        remove.add ( k ) 
for k in remove :
    del config[k]

# =============================================================================
## the final action...
import atexit
@atexit.register
def config_goodby () :
    import  datetime
    now = datetime.datetime.now() 
    if files_read :
        n = len ( files_read )
        if 1 == n :
            f = files_read[0]
            for fn in the_files :
                if f == fn[0] :
                    f = fn[1]
                    break                
            logger.info  ( 'The configuration of Ostap was read from %s' % f )
        else :
            import ostap.logger.table as T
            rows = [ ( '', 'file' ) ]
            for i, ff in enumerate ( files_read , start = 1 ) :
                f = ff 
                for fn in the_files :
                    if f == fn[0] :
                        f = fn[1]
                        break    
                row = '%d' % i , f 
                rows.append ( row )
            title = 'Configuration'
            table = T.table ( rows , title = title , prefix = '# ' , alignment = 'rl' )
            logger.info ( 'Configuration is read from\n%s' % table ) 
            
    import io 
    with io.StringIO() as o : 
        config.write( o )
        logger.verbose ( 'Ostap configuration:\n%s' % o.getvalue() )
    try :
        dump = '.ostap_config.txt'
        if os.path.exists ( dump ) : os.remove ( dump )
        with open ( dump , 'w' ) as ff :
            ff.write('# ' + 78*'*' + '\n')
            if files_read : 
                ff.write('# Ostap configuration read from:\n' )
                for i,f in enumerate ( files_read , start = 1 ) : ff.write('#  %2d. %s\n' % ( i , f ) ) 
            else :
                ff.write('# Ostap configuration:\n' )                
            ff.write('# ' + 78*'*' + '\n')    
            config.write( ff )            
            ff.write('# ' + 78*'*' + '\n')
            ff.write('# Configuration saved at %s\n' % now.strftime('%c') )
            ff.write('# ' + 78*'*' + '\n')            
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
