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
print ( 'CONFIG/1 Loading ostap.core.config' )

import configparser, os, sys  
import ostap.core.default_config as     default_config 
from   ostap.utils.env           import ( get_env            ,
                                          has_env            ,
                                          OSTAP_CONFIG       ,
                                          OSTAP_BATCH        ,
                                          OSTAP_SILENT       ,
                                          OSTAP_NCPUS        ,
                                          OSTAP_WEB_DISPLAY  , 
                                          OSTAP_SHOW_UNICODE , 
                                          OSTAP_PARALLEL     ,
                                          OSTAP_BUILD_DIR    ,
                                          OSTAP_CACHE_DIR    ,
                                          OSTAP_TMP_DIR      )
# =============================================================================

print ( 'CONFIG/2 Define configparser' ) 

## print for configparger 
def _cp_str_ ( cp ) :
    """ Print for config-parger
    """
    import io 
    with io.StringIO() as o :
        config.write( o )
        return o.getvalue()

config = configparser.ConfigParser()
type(config).__str__  = _cp_str_
type(config).__repr__ = _cp_str_

## define the major section:
config['General'] = { 
    'Quiet'       : str ( default_config.quiet        ) ,
    'Verbose'     : str ( default_config.verbose      ) ,
    'Parallel'    : str ( default_config.parallel     ) ,
    'WebDisplay'  : str ( default_config.web          ) , 
    'ShowUnicode' : str ( default_config.show_unicode ) , 
    'Parallel'    : str ( default_config.parallel     ) ,  
    'Batch'       : str ( default_config.batch        ) ,
    'BuildDir'    : str ( default_config.build_dir    ) ,
    'CacheDir'    : str ( default_config.cache_dir    ) ,
    'TmpDir'      : str ( default_config.tmp_dir      ) ,
}
# ===========================================================================

## generic TCanvas configuration
config [ 'Canvas'      ] = { 'Width'       :  '1000' , 'Height'       :  '800' , 
                             'MarginTop'   : '0.05'  , 'MarginBottom' : '0.12' ,
                             'MarginRight' : '0.07'  , 'MarginLeft'   : '0.12' }

config [ 'Fit Draw'    ] = {} ## RooFit plotting configuration
config [ 'Tables'      ] = {} ## configuration for Tables 
config [ 'RooFit'      ] = {} ## RooFit configuration

config [ 'Pathos'      ] = {} ## PATHOS configuration  
config [ 'IPyparallel' ] = {} ## ipyparallel configuration 

## the list of config files to be proessed
config_files = default_config.config_files 
if has_env ( OSTAP_CONFIG ) : 
    config_files = get_env ( OSTAP_CONFIG , '' , silent = True )
    config_files = tuple ( config_files.split ( os.pathsep ) ) 
    
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
    
the_files    = tuple ( the_files ) 
config_files = tuple ( f [ 0 ] for f in the_files ) 

## read the files 
files_read = config.read ( config_files )

# =============================================================================
## General section:
general      = config [ 'General' ]

if get_env ( OSTAP_SILENT , '' , silent = True ) : 
    general [ 'Quiet'   ] = 'True'
    general [ 'Verbose' ] = 'False'
    
if get_env ( OSTAP_SHOW_UNICODE , '' , silent = True ) : 
    general [ 'ShowUnicode' ] = 'True'
   
# ============================================================================
## redefine ncpus from the environment variable 
if has_env ( OSTAP_NCPUS ) : # ===============================================
    # ========================================================================
    ncpus_ = get_env ( OSTAP_NCPUS , '' , silent = True  )        
    # ========================================================================
    try : # ==================================================================
        # ====================================================================
        ncpus_ = int ( ncpus_ )
        if 1 <= ncpus_ : general [ 'NCPUs' ] = str ( ncpus_ ) 
        # ====================================================================
    except: # ================================================================
        # ====================================================================
        pass 

# ============================================================================
## redefine web display from the environment variable
if has_env ( OSTAP_WEB_DISPLAY ) :
    web_ = get_env ( OSTAP_WEB_DISPLAY , '' , silent = True  )        
    if web_ : general [ 'WebDisplay' ] = str ( web_ ) 

# ============================================================================
## redefine parallel from the environment variable
if has_env ( OSTAP_PARALLEL ) :
    parallel_ = get_env ( OSTAP_PARALLEL , '' , silent = True  )        
    if parallel_ : general [ 'Parallel' ] = str ( parallel_ )
    
if not general.getboolean ( 'Batch' , fallback = False ) and has_env ( OSTAP_BATCH ) :
    batch_ = get_env ( OSTAP_BATCH , '' , silent = True  )        
    if batch_ : general [ 'Batch' ] = 'True'
    
if not general.getboolean ( 'Batch' , fallback = False ) :
    batch_ = any ( a.lower() in ( '-b' , '--batch' , '--no-gui' ) for a in sys.argv )
    if batch_ : general [ 'Batch' ] = 'True'
      
if has_env ( OSTAP_BUILD_DIR ) :
    build_dir_ = get_env ( OSTAP_BUILD_DIR , '' , silent = True  )        
    if build_dir_ : general [ 'BuildDir' ] = str ( build_dir_ )
    
if has_env ( OSTAP_CACHE_DIR ) :    
    cache_dir_ = get_env ( OSTAP_CACHE_DIR , '' , silent = True  )        
    if cache_dir_ : general [ 'CacheDir' ] = str ( cache_dir_ ) 
    
if  has_env ( OSTAP_TMP_DIR ) :
    tmp_dir_ = get_env ( OSTAP_TMP_DIR , '' , silent = True  )        
    if tmp_dir_ : general [ 'TmpDir' ] = str ( tmp_dir_ )   
    
# =============================================================================
## Some explicit & important elements from the `General` section:
quiet        = general.getboolean ( 'Quiet'       , fallback = default_config.quiet        )
verbose      = general.getboolean ( 'Verbose'     , fallback = default_config.verbose      )
show_unicode = general.getboolean ( 'ShowUnicode' , fallback = default_config.show_unicode )
ncpus        = general.getint     ( 'NCPUs'       , fallback = default_config.ncpus        )
webdisplay   = general.get        ( 'WebDisplay'  , fallback = default_config.web          )
parallel     = general.get        ( 'Parallel'    , fallback = default_config.parallel     )
batch        = general.getboolean ( 'Batch'       , fallback = False                       )
build_dir    = general.get        ( 'BuildDir'    , fallback = default_config.build_dir    )
cache_dir    = general.get        ( 'CacheDir'    , fallback = default_config.cache_dir    )
tmp_dir      = general.get        ( 'TmpDir'      , fallback = default_config.tmp_dir      )

# =============================================================================
## Section with canvas configuration
canvas  = config [ 'Canvas'    ]

# =============================================================================
## Section for RooFit 
roofit   = config [ 'RooFit'   ]

# =============================================================================
## section for fit drawing options 
fit_draw = config [ 'Fit Draw' ]

# =============================================================================
## Section for tables  
tables   = config [ 'Tables'   ]
# =============================================================================

# =============================================================================
## reconfigure logging according to configuration 
import logging
logging.disable ( ( logging.WARNING - 1 ) if quiet   else
                  ( logging.DEBUG   - 5 ) if verbose else ( logging.INFO - 1 ) )

# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.core.config' )
else                       : logger = getLogger ( __name__  )
# =============================================================================

# =============================================================================
## the final action...
import atexit
@atexit.register
def config_goodby () :
    import  datetime, sys
    if not hasattr ( sys , 'ps1' ) : return 
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
    
print ( 'CONFIG/3 Ostap configuration is loaded' )

# =============================================================================
if '__main__' == __name__ :

    print   ( 'CONFIG/4 Running ostap.core.config as main' )    
    
    def _cp_hash_ ( cp ) : return hash ( str ( cp ) )
    type(config).__hash__ = _cp_hash_

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    cnf = '\n' + str(config)
    cnf = cnf.replace ('\n','\n# ')
    logger.info ( 'Ostap configuration is:%s' % cnf )
    
    print  ( 'CONFIG/5 Done' )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
