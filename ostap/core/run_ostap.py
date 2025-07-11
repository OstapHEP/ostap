#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# =============================================================================
## @file run_ostap.py
#  
#     .oooooo.                .                        
#    d8P'  `Y8b             .o8                        
#   888      888  .oooo.o .o888oo  .oooo.   oo.ooooo.  
#   888      888 d88(  "8   888   `P  )88b   888' `88b 
#   888      888 `"Y88b.    888    .oP"888   888   888 
#   `88b    d88' o.  )88b   888 . d8(  888   888   888 
#    `Y8bood8P'  8""888P'   "888" `Y888""8o  888bod8P' 
#                                            888       
#                                           o888o      
#                                                    
#  Simple interactive PyRoot-based analysis environment to provide access
#  to zillions useful decorators for ROOT (and not only ROOT!) objects&classes  
# 
#  This file is a part of 
#
#  @date   2012-02-15
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#
# =============================================================================
""" Simple interactive PyRoot-based analysis environment
to provide access to zillions useful decorators for ROOT
(and not only ROOT) objects&classes

"""
# =============================================================================
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = "2012-09-10"
__version__ = '$Revision$'
# =============================================================================
## start form partisg the command-line arguments 
from   ostap.core.parse_args import parse_args
import os, sys 
from   io      import StringIO
## ROOT.PyConfig.IgnoreCommandLineOptions = True
# =============================================================================
## parse arguments
# =============================================================================
arguments = parse_args()

# =============================================================================
## do we run IPython?
from ostap.utils.basic import with_ipython

# =============================================================================

# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if arguments.Color and not arguments.batch :
    from ostap.logger.logger import make_colors
    make_colors()
    del make_colors
logger = getLogger( 'ostap' )
# =============================================================================

# =============================================================================
## suppress extra prints from logging  
# =============================================================================
import logging

level = logging.INFO - 2

if arguments.Quiet :    
    logger.info ( '(silent) Interactive Ostap session (steroid-enhanced PyROOT)')
    level = logging.WARNING - 1
elif arguments.Debug      : level = logging.DEBUG   - 2 
elif arguments.Verbose    : level = logging.VERBOSE - 2
elif 0 == arguments.Level : level = 1 
elif 1 == arguments.Level : level = logging.VERBOSE   - 2   
elif 2 == arguments.Level : level = logging.DEBUG     - 2  
elif 3 == arguments.Level : level = logging.INFO      - 2 
elif 4 == arguments.Level : level = logging.ATTENTION - 2 
elif 5 == arguments.Level : level = logging.WARNING   - 2 
elif 6 == arguments.Level : level = logging.ERROR     - 2
elif 7 <= arguments.Level : level = logging.CRITICAL  - 2


logging.disable ( level  )

if level <= logging.INFO : 
    
    from ostap import banner
    logger.info ( "Welcome to Ostap\n" +  banner )
    logger.info ( __doc__ )
    del banner
    
    _vars = vars ( arguments )
    _keys = _vars.keys()
    rows = [ ( 'Argument' , 'Value' ) ] 
    logger.info ( 'Arguments  : ')
    for _k in sorted ( _keys ) :
        row = _k , str ( _vars [ _k  ] ) 
        rows.append ( row )
    title = 'Ostap arguments'
    import ostap.logger.table as T
    rows  = T.table ( rows , title = title , prefix = '# ' , alignment = 'rw' )
    logger.info ( '%s:\n%s' % ( title , rows ) ) 
    del _keys,_vars,_k,level  
    
# =============================================================================
if arguments.Config :
    from ostap.utils.env import get_env, OSTAP_CONFIG 
    cc  = get_env ( OSTAP_CONFIG , '' ).split( os.pathsep )
    cc += arguments.Config
    cc  = os.pathsep.join ( cc )
    os.environ [ OSTAP_CONFIG ] = cc

import ostap.core.config     

# =============================================================================
## Web Display
# =============================================================================
if arguments.web :
    ostap.core.config.general['WebDisplay'] = arguments.web

# =============================================================================
## use profiling ?
if arguments.Profile :
    import cProfile as _profile
    _pr = _profile.Profile()
    
    def _end_profile_ ( prof ) :
        logger.debug ( 'End of profiling, generate profiling report' )
        prof.disable()
        _sio   = StringIO()
        sortby = 'cumulative'
        import pstats 
        _pstat = pstats.Stats( prof , stream=_sio).sort_stats('cumulative')
        _pstat.print_stats()
        print(_sio.getvalue())
            

    import atexit 
    atexit.register ( _end_profile_ , _pr ) 
    _pr.enable()
    del _pr 
    logger.info ( 'Profiling is activated' )


# =============================================================================
## set ROOT into batch mode 
# =============================================================================
import ostap.fixes.fixes
import ROOT 
groot = ROOT.ROOT.GetROOT() 
groot.SetBatch ( arguments.batch )
if groot.IsBatch() : logger.info ('Batch processing is activated') 


# =============================================================================
# specify the build directory for ROOT 
# =============================================================================
import ostap.core.build_dir
if arguments.build_dir :
    
    from ostap.utils.basic   import make_dir
    from ostap.utils.cleanup import writeable 
    bdir = arguments.build_dir
    
    if   writeable ( bdir )                   : pass 
    elif bdir and not os.path.exists ( bdir ) : make_dir ( bdir ) 

    if writeable ( bdir ) :
        ROOT.gSystem.SetBuildDir ( bdir )
        logger.info ( 'Build directory for ROOT: %s' % ostap.core.build_dir.build_dir )
        ostap.core.build_dir.build_dir = bdir 
        
    del bdir, writeable , make_dir

# =============================================================================
## ostap startup: history, readlines, etc... 
# =============================================================================
import ostap.core.startup

# =============================================================================
## import everything from ostap
# =============================================================================
if arguments.Quiet : 
    from ostap.logger.utils import mute
    with mute () : 
        from   ostap.core.load_ostap import *
    del mute
else :
    from ostap.core.load_ostap import *

# =============================================================================
## create default canvas
# =============================================================================
if arguments.canvas :
    
    import ostap.plotting.canvas 
    logger.debug ( "Create the default canvas" )
    canvas    = ostap.plotting.canvas.getCanvas ()

# =============================================================================
## execute startup files 
# =============================================================================
## startup files to be executed:
_startups = ( '$OSTAPSTART'     ,  
              '$HOME/.ostap.py' ,
              '~/.ostap.py'     ,        
              './.ostap.py'     )
_executed = set() 
for _s in _startups : 
    import os
    _ss = _s 
    _ss =  os.path.expandvars ( _ss )
    _ss =  os.path.expandvars ( _ss )
    _ss =  os.path.expandvars ( _ss )
    _ss =  os.path.expanduser ( _ss )
    _ss =  os.path.expandvars ( _ss )
    if not os.path.exists     ( _ss ) : continue
    if not os.path.isfile     ( _ss ) : continue
    _ss =  os.path.abspath    ( _ss )
    if _ss in _executed           : continue
    ## execute it!
    # =========================================================================
    try : # ===================================================================
        # =====================================================================

        import runpy
        globs = runpy.run_path ( _ss , globals() , run_name = '__main__')
        globs = dict ( ( (k,v) for k,v in globs.items() if  not k.startswith('__') and not k.endswith('__')) )
        logger.debug('Symbols from %s: %s' % ( _s , globs.keys() ) )
        globals().update( globs )
        del globs
        
        _executed.add ( _ss )
        logger.info   ( "Startup file '%s' is executed"      % _s )
        
        # =====================================================================
    except: # =================================================================
        # =====================================================================
        logger.error  ( "Error in execution of '%s' startup" % _s , exc_info = True )
        

# =============================================================================
## cleanup a bit the context 
# =============================================================================
del parse_args, getLogger
del _startups,_executed,_s,_ss 

# =============================================================================
## execute the files, defined as arguments
# =============================================================================
root_files     = {}
root_macros    = []
python_scripts = []
PARAMETERS     = []
## _Load roo macro
def _load_macro_ ( macro , silent = True ) :
    """Load ROOT macro"""
    logger.debug  ("Try to load macro '%s'" % macro )
    groot = ROOT.ROOT.GetROOT()
    from   ostap.core.core     import rootError
    
    if silent :
        with rootError() : sc = groot.LoadMacro ( macro )
    else :
        sc = groot.LoadMacro ( macro )
        
    if sc :
        # - Interactive mode: print traceback and continue
        if arguments.batch :
            raise RuntimeError ( 'Failure to load macro "%s" code:%d' % ( macro , sc ) )
        else :
            logger.error       ('Failure to load macro "%s" code:%d' % ( macro , sc ) )
            return False
        
    logger.debug ("Loaded macro   '%s'" % macro )
    del sc
    return True
    
## Load root macro/pattern
def load_macro ( pattern ) :
    """ Load root macro/pattern
    """
    ## expand the actual file name 
    pattern  = os.path.expandvars ( pattern )
    pattern  = os.path.expanduser ( pattern )
    pattern  = os.path.expandvars ( pattern )
    patetrn  = os.path.expandvars ( pattern )
    ##
    patterns = [ (pattern,'') ]
    if   pattern.endswith ( '++' ) : patterns += [ ( pattern [:-2] , '++' ) ]
    elif pattern.endswith ( '+'  ) : patterns += [ ( pattern [:-1] , '+'  ) ]

    import glob
    _glob   = False 
    _loaded = 0
    
    for pat,ext in patterns :
        
        for m in glob.iglob ( pat ) :

            _glob = True 
            if m in root_macros  : continue

            loaded = _load_macro_ ( m )
            if loaded : root_macros.append ( m )
            _loaded += loaded 
            
    if not _glob and not pattern in root_macros :
        m = pattern
        if _load_macro_ ( m ) : 
            root_macros.append ( m )
            _loaded += 1 
        
    return _loaded 
# =============================================================================    
## perform treatment of the file 
def treat_file ( f ) :
    """Perform treatment of the file 
    """
    fok = os.path.exists ( f ) and os.path.isfile ( f ) 
    name,dot,ext  = f.rpartition('.')
    if  name and dot and ext in ('root', 'ROOT' ) :
        
        logger.debug  ("Try  to open file '%s'" % f )
        try :
            from ostap.core.core import ROOTCWD 
            with ROOTCWD() :
                _f = ROOT.TFile.Open(  f , 'READ')
                root_files [ name ] = _f 
                logger.info ("Open ROOT file '%s'" % f )
                _f.ls()
        except :
            
            # - Batch mode:       just re-raise exception
            if arguments.batch :  raise
            
            # - Interactive mode: print traceback and continue 
            logger.error ('Failure to open ROOT file "%s"'  % f  , exc_info = True ) 
            
    ## execute python file 
    elif fok and name and dot and 'py' == ext  :  
        
        logger.debug  ("Try    to execute '%s'" % f )
        
        try :

            import runpy 
            globs = runpy.run_path (
                f , globals() if arguments.WithContext else {}  ,
                run_name = '__main__' )
            globs = dict ( ( (k,v) for k,v in globs.items() if  not k.startswith('__') and not k.endswith('__')) )
            logger.debug('Symbols from %s: %s' % ( f , globs.keys() ) )
            globals().update ( globs )
            del globs 
                
            logger.info  ("Executed       '%s'" % f )
            python_scripts.append ( f )
            
        except :
            
            # - Batch mode:       just re-raise exception
            if arguments.batch :  raise
            
            # - Interactive mode: print traceback and continue 
            logger.error ('Failure to execute "%s"'     % f  , exc_info = True ) 
            
    elif name and dot and ext.upper() in ( 'C'   ,   'C+' ,   'C++' ,
                                           'CPP' , 'CPP+' , 'CPP++' ,
                                           'CXX' , 'CXX+' , 'CXX++' ) :
        
        if not load_macro ( f )  :
            if arguments.batch :
                raise RuntimeError ( "No macros are loaded for '%s' pattern" % f )
            else :
                logger.error       ( "No macros are loaded for '%s' pattern" % f )

    ## execute ostap-script 
    elif fok and name and dot and ext in ( 'ost' , 'ostp' , 'ostap' , 'opy' ) :
        
        logger.debug  ("Try    to execute '%s'" % f )

        # =====================================================================
        try : # ===============================================================
            # =================================================================
            
            ## this one is always executed with the context!!!
            import runpy
            globs = runpy.run_path ( f , globals() , run_name = '__main__')
            globs = dict ( ( (k,v) for k,v in globs.items() if  not k.startswith('__') and not k.endswith('__')) )
            logger.debug('Symbols from %s: %s' % ( f , globs.keys() ) ) 
            globals().update ( globs )
            del globs 
            
            python_scripts.append ( f )
            logger.info  ("Executed       '%s'" % f )

            # =================================================================
        except : # ============================================================
            # =================================================================
                
            # - Batch mode:       just re-raise exception
            if arguments.batch :  raise
            
            # - Interactive mode: print traceback and continue 
            logger.error ('Failure to execute "%s"'     % f  , exc_info = True ) 

        # =====================================================================
    else : # ==================================================================
        # =====================================================================
                
        ## collect the argument as parameters
        PARAMETERS.append ( f )
        logger.debug ( 'Add argument %s to PARAMETERS' % f )

# =============================================================================
## load ROOT macros 
for pattern in arguments.Macros :
    if not load_macro ( pattern )  :
        if arguments.batch :
            raise RuntimeError ( "No macros are loaded for '%s' pattern" % pattern )
        else :
            logger.error       ( "No macros are loaded for '%s' pattern" % pattern )

# =============================================================================
## treat all input arguments/files 
for pattern in arguments.files :
    import glob
    _glob = False
    ##
    pattern  = os.path.expandvars ( pattern )
    pattern  = os.path.expanduser ( pattern )
    pattern  = os.path.expandvars ( pattern )
    patetrn  = os.path.expandvars ( pattern )
    ##
    for _f in glob.iglob ( pattern ) :
        _glob = True
        treat_file ( _f      )
    if not _glob :
        treat_file ( pattern )
del treat_file

    
# =============================================================================
## list names of opened ROOT files 
if root_files :
    logger.info("ROOT FILES    : %d (with '.root' extension)" % len ( root_files ) )
    keys = root_files.keys()
    _n = 0  
    for k in sorted ( keys ) :
        _n +=  1
        rfile = root_files [k]
        okeys = list ( rfile.keys () )
        okeys.sort()
        logger.info('%2d: %s #keys: %d: %s' % ( _n , k , len(okeys) , okeys ) )

# =============================================================================
## list all loaded root macros
if root_macros    : 
    logger.info ('ROOT MACROS   : %s' % root_macros    )
# =============================================================================
## list all executed ostap/python scripts 
if python_scripts :
    logger.info ('OSTAP SCRIPTS : %s' % python_scripts )
    
# =============================================================================
## dump all loaded 'parameters'
if PARAMETERS : 
    logger.info ('PARAMETERS    : %s' % PARAMETERS )


## make `reload` command available 
from importlib import reload 

# =============================================================================
if '__main__' == __name__ : logger.info ( 'ostap is ready'  ) 
else                      : logger.info ( 'ostap is loaded' )

# =============================================================================
##                                                                      The END 
# =============================================================================

