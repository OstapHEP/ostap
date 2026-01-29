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
__version__ = '$Revision:$'
# =============================================================================
import ostap.core.config as config 
import os, sys
# =============================================================================
## ROOT.PyConfig.IgnoreCommandLineOptions = True
# =============================================================================

# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger, enabledInfo     
logger = getLogger ( 'ostap' )
# =============================================================================
## The command line argumenst 
arguments = config.arguments
# =============================================================================
if enabledInfo () : # =========================================================
    # =========================================================================
    ## Banner ?
    # =========================================================================
    from ostap import banner
    logger.info ( "Welcome to Ostap\n" +  banner )
    logger.info ( __doc__ )
    del banner
    # ========================================================================
    ## Table of command line arguments 
    rows  = [ ( 'Argument' , 'Value' ) ]
    vars_ = vars ( arguments )
    for key in sorted (vars_.keys ()  ) :
        row = key , str ( vars_ [ key ] ) 
        rows.append ( row )
        
    title = 'Command line arguments'
    import ostap.logger.table as T
    rows  = T.table ( rows , title = title , prefix = '# ' , alignment = 'rw' )
    logger.info ( '%s:\n%s' % ( title , rows ) ) 
    del rows, vars_, title,

# =============================================================================
## Basic setup for ROOT: Batch, Implicit MT, directories, profile 
# =============================================================================
import ostap.core.ostap_setup 

# =============================================================================
## ostap startup: history, readlines, etc... 
# =============================================================================
import ostap.core.startup

# =============================================================================
## import everything from ostap
# =============================================================================
if arguments.Quiet or arguments.Silent : # ====================================
    # =========================================================================
    from ostap.logger.utils import mute
    with mute () : # ==========================================================
        from   ostap.core.load_ostap import *
    del mute
    # =========================================================================
else : # ======================================================================
    # ========================================================================
    from ostap.core.load_ostap import *

# =============================================================================
## create default canvas
# =============================================================================
if arguments.Canvas :
    import ostap.plotting.canvas 
    logger.debug ( "Create default Ostap canvas" )
    canvas    = ostap.plotting.canvas.getCanvas ()

# =============================================================================
## execute startup files 
# =============================================================================
_executed = set() 
for sf in arguments.StartUp :
    _ss =  sf
    _ss =  os.path.expandvars ( _ss )
    _ss =  os.path.expandvars ( _ss )
    _ss =  os.path.expandvars ( _ss )
    _ss =  os.path.expanduser ( _ss )
    _ss =  os.path.expandvars ( _ss )
    if not os.path.exists     ( _ss ) : continue
    if not os.path.isfile     ( _ss ) : continue
    _ss =  os.path.abspath    ( _ss )
    if _ss in _executed           : continue
    # ==========================================================================
    ## execute it!
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import runpy
        globs = runpy.run_path ( _ss , globals() , run_name = '__main__' )
        globs = dict ( ( ( k , v ) for k , v in globs.items() if  not k.startswith ( '__' ) and not k.endswith ( '__' ) ) )
        logger.debug('Symbols from %s: %s' % ( sf , globs.keys() ) )
        globals().update( globs )
        del globs
        ## 
        _executed.add ( _ss )
        logger.info   ( "Startup file '%s' is executed"      % sf )
        ## 
        # =====================================================================
    except: # =================================================================
        # =====================================================================
        logger.error  ( "Error in execution of '%s' startup" % sf , exc_info = True )


# =============================================================================
## suppress excessive (?) RooFit printout
# =============================================================================
if ( arguments.Quiet or arguments.Silent ) and not arguments.Verbose and 3 < arguments.Level :
    from ostap.fitting.utils import suppress_topics
    suppress_topics ( "Plotting"           ,
                      "Caching"            ,
                      "Eval"               , 
                      "Integration"        ,
                      "NumericIntegration" ,
                      "Fitting"            ,
                      "InputArguments"     ,
                      "ObjectHandling"     )
    
# =============================================================================
## execute the files, defined as arguments
# =============================================================================
root_files     = {}
root_macros    = []
python_scripts = []
PARAMETERS     = []
# =============================================================================
## _Load roo macros
def _load_macro_ ( macro , silent = True ) :
    """ Load ROOT macro
    """
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

# =========================================================================
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
    """ Perform treatment of the file 
    """
    fok = os.path.exists ( f ) and os.path.isfile ( f ) 
    name , dot , ext  = f.rpartition('.')
    
    if  name and dot and ext in ('root', 'ROOT' ) :
        
        logger.debug  ("Try  to open ROOT file '%s'" % f )
        try :
            from  ostap.core.core      import ROOTCWD
            import ostap.io.root_file 
            with ROOTCWD() :
                _f = ROOT.TFile.Open(  f , 'READ')
                root_files [ name ] = _f 
                logger.info ("Open ROOT file '%s'" % f )
                _f.ls()
        except :
            
            # - Batch mode:       just re-raise exception
            if arguments.batch :  raise
            
            # - Interactive mode: print traceback and continue 
            logger.error ( 'Failure to open ROOT file "%s"'  % f  , exc_info = True ) 
            
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
            logger.debug ( 'Symbols from %s: %s' % ( f , globs.keys() ) ) 
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
            logger.error ( 'Failure to execute "%s"'     % f  , exc_info = True ) 

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
for command in arguments.Commands :
        logger.debug  ("Try  to execute command '%s'" % command  )
        try :
            exec ( command )
        except :
            if arguments.batch : raise
            logger.error       ( "Failure ro execute command: '%s'" % command , exc_info = True  )

            
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

