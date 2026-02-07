#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/ostap_setup.py
#  Basic setup for ROOT/Ostap
#  - Batch
#  - Implicit MR
#  - Build/Cache/TMP  directories
#  - Profiling 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
""" Basic setup for ROOT/Ostap 
- Batch
- Implicit MR
- Build/Cache/TMP  directories 
- Profiling 
"""
# ============================================================================= 
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    'executed_scripts' , ## list of successfully executed ascripts 
    'loaded_macros'    , ## list of successfully loaded   ROOT/C++ macros 
    'executed_macros'  , ## list of successfully executed ROOT/C++ macros 
    'root_files'       , ## list of ROOT files  
    'parameters'       , ## list of extra comand-line arguments
    'imported'         , ## dictionary  of successfully imported python symbols
    'run_py'           , ## function to run/import python scripts
)
# ============================================================================= 
from   ostap.core.base        import rootException 
import ostap.core.config      as     config 
import ostap.core.build_dir
import ostap.core.cache_dir
import ostap.utils.cleanup
import ROOT, math, os, sys, ctypes, runpy     
# =============================================================================
from ostap.logger.logger import getLogger, enabledInfo, enabledDebug  
if '__main__' ==  __name__ : logger = getLogger( 'ostap.core.ostap_setup' )
else                       : logger = getLogger( __name__                )
# =============================================================================
if not '.' in sys.path :
    logger.debug('Prepend sys.path with $PWD')
    sys.path = ['.'] + sys.path
    
# =============================================================================
## setup some most imporant properties:
# =============================================================================

# =============================================================================
## (1) Implicit MT ?
# =============================================================================
implicitMT = config.general.getboolean ( 'ImplicitMT' , fallback = False )
if       implicitMT and not ROOT.ROOT.IsImplicitMTEnabled() :
    ROOT.ROOT.EnableImplicitMT  ()
    if     ROOT.ROOT.IsImplicitMTEnabled() : logger.debug ( 'ImplicitMT is enabled'  )
elif not implicitMT and     ROOT.ROOT.IsImplicitMTEnabled() :
    ROOT.ROOT.DisableImplicitMT  ()
    if not ROOT.ROOT.IsImplicitMTEnabled() : logger.debug ( 'ImplicitMT is disabled' )
    
# =============================================================================
## (2)  Batch processing ?
# =============================================================================
groot = ROOT.ROOT.GetROOT() 
if   groot and     config.batch and not groot.IsBatch () :
    groot.SetBatch ( True  )
    if     groot.IsBatch() : logger.attention ( "BATCH processing is activated!"   )
elif groot and not config.batch and     groot.IsBatch () :
    groot.SetBatch ( False )
    if not groot.IsBatch() : logger.info      ( "BATCH processing is deactivated!" )

# =============================================================================
## (3) Root & RooFit print levels 
# =============================================================================
groot = ROOT.ROOT.GetROOT()
msg   = ROOT.RooMsgService.instance()
if config.silent :    
    if groot : groot.ProcessLine ( "gErrorIgnoreLevel= %d ; " % ROOT.kError )     
    msg.setGlobalKillBelow  ( ROOT.RooFit.ERROR )
    msg.setSilentMode       ( True  )    
elif config.quiet :    
    if groot : groot.ProcessLine ( "gErrorIgnoreLevel= %d ; " % ROOT.kWarning )     
    msg.setGlobalKillBelow  ( ROOT.RooFit.WARNING )
    msg.setSilentMode       ( True  )    
elif config.debug or config.verbose :    
    if groot : groot.ProcessLine ( "gErrorIgnoreLevel= %d ; " % ROOT.kPrint )     
    msg.setGlobalKillBelow  ( ROOT.RooFit.DEBUG )
    msg.setSilentMode       ( False )    
elif 0 <= config.level <= 8 :    
    d1 = ROOT.kFatal       - ROOT.kPrint
    d2 = ROOT.RooFit.FATAL - ROOT.RooFit.DEBUG 
    l1 = ROOT.kFatal       + math.floor ( ( d1 * config.level ) / 8 )
    l2 = ROOT.RooFit.DEBUG + math.floor ( ( d2 * config.level ) / 8 )
    if groot : groot.ProcessLine ( "gErrorIgnoreLevel= %d ; " % l1 )     
    msg.setGlobalKillBelow       ( l2 )
    msg.setSilentMode            ( ROOT.RooFit.INFO < l2 )

# =============================================================================
## RooFit topics..
# =============================================================================

# =============================================================================
## (5) actiavate profiling ?
# =============================================================================
if config.general.getboolean ( 'Profile' , fallback = config.profile ) :

    ## dedicated logger 
    plogger  = getLogger ( 'ostap.core.profiler' )

    ## profiler itself
    import cProfile as profile
    profiler = profile.Profile()
    profiler.enable()
    plogger .attention ( 'Profiling is activated!' )
    
    # =========================================================================
    ## "AtExit" action:
    #  - end profiling
    #  - print statistics     
    def _profile_atexit_ ( prof , logger ) :
        """ `At-Exit` action: 
        - end profiling 
        - print statistics 
        """       
        import io, pstats 
        
        prof.disable()
        logger.info  ( 'End of profiling, generate profiling report' )
        
        sio    = io.StringIO ()
        sortby = pstats.SortKey.CUMULATIVE 
        stat   = pstats.Stats ( prof , stream= sio ).sort_stats ( sortby )
        stat.print_stats()

        report = sio.getvalue()
        logger.info ( 'Profiler report:\n\n%s' % report  )
        logger.info ( 'The end of profiler report' ) 
        
    import atexit 
    atexit.register ( _profile_atexit_ , profiler , plogger ) 


# =============================================================================
## Run/execute python script
#  The overall effect somehow similar to <code>from ... import *</code>
#  @code
#  path = ...
#  globs = run_py ( path )
#  globals().update ( globs )
#  del globs
#  @endcode
def run_py ( path , run_name = '' , with_context = True , path_name = '' ) :
    """ Run/execute python script
    The overall effect somehow similar to `from ... import *`
    
    >>> path = ...
    >>> globs = run_py ( path )
    >>> globals().update ( globs )
    """
    if not path_name : path_name = path
    if not run_name :
        fbase , dot , ext = path.rpartition ( '.' )
        run_name = os.path.basename ( fbase ) if fbase and dot and ext else path 
        
    ## get all the keys for the current globals 
    _globals  = globals().copy() if with_context else {}
    ## 
    ## ATTENTION: we need to suppress __all__ 
    _globals.pop ( '__all__' , None )
    ##
    _k_globs  = frozenset ( _globals.keys() )
    _r_globs  = runpy.run_path ( path , None , run_name = run_name )

    if '__all__' in _r_globs :
        _r_all   = _r_globs.get ( '__all__' , () )
        _r_globs = {  k : v for k , v in _r_globs.items () if k in _r_all }
    else :
        _r_globs = {  k : v for k , v in _r_globs.items () if not k.startswith ( '_' ) }

    symbols_updated = frozenset ( k for k in _r_globs if     k in _k_globs )
    symbols_new     = frozenset ( k for k in _r_globs if not k in _k_globs )

    if enabledDebug () and ( symbols_updated or symbols_new ) : 
    
        from   ostap.utils.basic      import typename 
        rows   =  [ ( 'Symbol' , 'Type' , '' ) ]
    
        title1 = 'Symbols from %s' % path
        title2 = 'Symbols from %s' % path_name 
        
        if symbols_new and symbols_updated :
            for k in symbols_updated :
                row = k , typename ( _r_globs [ k ] ) , 'update'
                rows.append ( row ) 
            for k in symbols_new     :
                row = k , typename ( _r_globs [ k ] ) , 'new'
                rows.append ( row )            
        elif symbols_new :
            for k in symbols_new     :
                row = k , typename ( _r_globs [ k ] )
                rows.append ( row )                
            title1 = 'Imported from %s' % path
            title2 = 'Imported from %s' % path_name         
        elif symbols_updated:
            for k in symbols_updated :
                row = k , typename ( _r_globs [ k ] ) 
                rows.append ( row )                
            title1 = 'Updated from %s' % path
            title2 = 'Updated from %s' % path_name

        import ostap.logger.table as T
        rows = T.remove_empty_columns ( rows ) 
        logger.info ('%s:\n%s' % ( title1 , T.table ( rows , title = title2 , prefix = '# ' ) ) )
        
    elif symbols_new or symbols_updated :
        
        logger.info ('run_py(%s): Updated/New symbols %d/%d' % ( path ,
                                                                 len ( symbols_updated ) ,
                                                                 len ( symbols_new     ) ) ) 
        
    return _r_globs

# =============================================================================
## dictionary of successfully imported python symbols
imported = {}  ## dictionary of successfully imported python symbols

# =============================================================================
## (6) execute startup (python) files/scripts 
# =============================================================================
executed_scripts  =     []
_scripts          = set () 
for _su in config.startup_files :
    ## 
    _ss =  _su
    _ss =  os.path.expandvars ( _ss )
    _ss =  os.path.expandvars ( _ss )
    _ss =  os.path.expandvars ( _ss )
    _ss =  os.path.expanduser ( _ss )
    _ss =  os.path.expandvars ( _ss )
    if not os.path.exists     ( _ss ) : continue
    if not os.path.isfile     ( _ss ) : continue
    _ss =  os.path.abspath    ( _ss )
    if _ss in _scripts                : continue 
    ## 
    fname = _ss 
    fbase , dot , ext = fname.rpartition('.')
    bname = os.path.basename ( fbase )
    if not bname : bname = os.path.basename ( fname )
    ## 
    # =========================================================================
    ## execute it!
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        logger.debug ( "Execute the script: `%s` " % _su ) 
        ## 
        globs = run_py ( fname , run_name = bname , path_name = _su , with_context = False )
        imported.update ( globs )
        ##
        _scripts.add ( fname )
        executed_scripts.append ( ( _su , fname ) )
        ## 
        logger.debug ( "Script '%s' is executed"      % _su )
        ## 
        # =====================================================================
    except: # =================================================================
        # =====================================================================
        if config.batch : raise 
        logger.error  ( "Error in execution of script: '%s'" % _su , exc_info = True )

# =============================================================================
## (7) Load ROOT/C++ macros 
# =============================================================================
loaded_macros = [] 
for m in config.load_macros :
    ## 
    groot = ROOT.ROOT.GetROOT()
    if not groot : continue
    ##
    ## (1) just check the macro exists and readable
    cerror  = ctypes.c_int( 0 ) 
    sc      = groot.LoadMacro ( m , cerror , True ) ## CHECK MACRO 
    cerror  = cerror.value
    if sc or cerror : 
        msg = 'load_macro: macro `%s` is not found!' % m
        if config.batch : raise RuntimeError ( msg )
        else            : logger.error       ( msg )
        continue
    ## 
    ## (2) Load macro 
    logger.debug ( "Load  macro: `%s`" % m  )
    ##
    # ======================================================================
    try : # ================================================================
        # ==================================================================
        with rootException() :             
            cerror = ctypes.c_int( 0 )
            sc     = groot.LoadMacro ( m , cerror , False  )
            cerror = cerror.value 
            ##
        if sc or cerror :
            msg = 'Failure to load macro "%s", code:%d/%d' % ( m , sc , cerror )
            if config.batch : raise RuntimeError ( msg )
            else            : logger.error       ( msg )
            continue 
            ##
        loaded_macros.append ( m ) 
        logger.debug ( "Macro is loaded `%s`"      % m )
        # ==================================================================
    except : # =============================================================
        # ==================================================================
        msg = 'Failure to load macro "%s"' % m
        if config.batch : raise 
        else            : logger.error ( msg , exc_info = True )
                        
# =============================================================================
## (7) Execute ROOT/C++ macros 
# =============================================================================
executed_macros = [] 
for m in config.exec_macros :
    ## 
    logger.debug ( "Execute macro: `%s` " % m ) 
    groot = ROOT.ROOT.GetROOT()
    ##
    ## (1) just check the macro exists and readable
    ## 
    cerror = ctypes.c_int( 0 )
    sc     = groot.LoadMacro ( m , cerror , True )  ## CHECK MACRO 
    cerror = cerror.value
    if sc or cerror : 
        msg = 'exec_macro: macro `%s` is not found!' % m
        if config.batch : raise RuntimeError ( msg )
        else            : logger.error       ( msg )
        continue
    ## 
    ## (2) execute macro
    ##
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        with rootException () : 
            cerror  = ctypes.c_int( 0 )     
            sc      = groot.Macro ( m , cerror , not config.batch )
            cerror  = cerror.value 
        if sc or cerror :
            msg = 'Failure to execute macro "%s", code:%d/%d' % ( m , sc , cerror )
            if config.batch : raise RuntimeError ( msg )
            else            : logger.error       ( msg )
            continue 
        ## 
        executed_macros.append ( m ) 
        logger.debug ( "Macro '%s' is executed"      % _su )
        # =====================================================================
    except : # ================================================================
        # =====================================================================
        msg = 'Failure to execute macro "%s"' % m 
        if config.batch : raise 
        else            : logger.error ( msg , exc_info = True )
                
# =============================================================================
## (8) Execute the commands 
# =============================================================================
executed_commands = []
for _su in config.commands :
    ##
    if not _su : continue
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        logger.debug ( "Execute the command: `%s` " % _su ) 
        ## 
        exec ( _su , globals = globals() , locals = locals () )
        ## 
        executed_commands.append ( _su )
        logger.debug ( "Command '%s' is executed"      % _su )
        # =====================================================================
    except : # ================================================================
        # =====================================================================
        if config.batch : raise 
        logger.error ( 'Failure to execute command: "%s"' % _su )
        ##

# =============================================================================
## (10) Treat input files 
# =============================================================================
cpp_ext    = 'C' , 'CPP' , 'CXX' , 'C+'  , 'CPP+'  , 'CXX+' , 'C++' , 'CPP++' , 'CXX++'
root_files = []
parameters = [] 
for _su in config.input_files :

    _ss = _su
    _ss = os.path.expandvars ( _ss )
    _ss = os.path.expandvars ( _ss )
    _ss = os.path.expandvars ( _ss )
    _ss = os.path.expanduser ( _ss )
    _ss = os.path.expandvars ( _ss )
    _ss = os.path.abspath    ( _ss ) 

    fname = _ss 
    fbase , dot , ext = fname.rpartition('.')
    bname = os.path.basename ( fbase )
    if not bname : bname = os.path.basename ( fname )

    if fbase and dot and ext[:4] in ( 'root' , 'ROOT' ) :
        
        import ostap.io.root_file 
        ff = ROOT.TFile.Open ( fname , 'READ' , exception = False )
        if ff :
            
            fn = fname
            if os.path.exists ( fn ) and os.path.isfile ( fn ) : fn = os.path.abspath ( fn ) 
            root_files.append ( ( _su , fn ) )
            
            logger.debug ( "ROOT file '%s' checked " % _su )
            ff.Close() 
        else :
            logger.warning ( "ROOT file '%s' cannot be opened" % _su )
            
    ## python/ostap script ?
    elif fbase and dot and ext.upper() in ( 'PY' , 'OST' , 'OSTP' , 'OSTAP' , 'OPY' ) :

        # =====================================================================
        try : # ===============================================================
            # =================================================================
            globs = run_py ( fname , run_name = bname , with_context = False , path_name = _su )
            imported.update ( globs )
            ##
            executed_scripts.append ( ( _su ,fname  ) ) 
            logger.debug  ( "Python script '%s' is executed"      % _su  )
            # ==================================================================
        except : # =============================================================
            # ==================================================================
            if config.batch : raise
            logger.error  ( "Error in execution of '%s' script" % _su , exc_info = True )

    ## unknown arguments 
    else :
        
        parameters.append ( _su ) 
        logger.info  ( "Parameter '%s' is added at index %d" % ( _su , len ( parameters ) ) ) 
        
    
executed_scripts  = tuple ( executed_scripts  )
executed_macros   = tuple ( executed_macros   )
executed_commands = tuple ( executed_commands )
root_files        = tuple ( root_files        )
parameters        = tuple ( parameters        ) 

# =============================================================================
## list all executed ostap/python scripts 
if executed_scripts :
    rows = [  ( '#' , 'Script' ,  'Expanded name' ) ] 
    for i, ff in enumerate ( executed_scripts , start = 1 ) :
        fn , fa = ff 
        row = '%-2d' % i , fn , fa if fn != fa else '' 
        rows.append ( row )
    ##
    title = "Executed scripts #%d " % len ( executed_scripts ) 
    import ostap.logger.table as T
    table = T.table ( rows , title = title , alignment = 'cw' , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# =============================================================================
## list all executed ROOT/C++ macros 
if loaded_macros :
    rows = [  ( '#' , '' , ''  'Macro' ) ] 
    for i , m in enumerate ( loaded_macros , start = 1 ) :
        row = '%-2d' % i , m 
        rows.append ( row )
    ##  
    title = "Loaded macros #%d " % len ( executed_macros ) 
    import ostap.logger.table as T
    table = T.table ( rows , title = title , alignment = 'cw' , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# =============================================================================
## list all executed ROOT/C++ macros 
if executed_macros :
    rows = [  ( '#' , '' , '' , 'Macro'  ) ] 
    for i , m in enumerate ( executed_macros , start = 1 ) :
        row = '%-2d' % i , m 
        rows.append ( row )
    ##  
    title = "Executed macros #%d " % len ( executed_macros ) 
    import ostap.logger.table as T
    table = T.table ( rows , title = title , alignment = 'cw' , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# =============================================================================
## list names of opened ROOT files 
if root_files :
    rows = [  ( '#' , 'File' , 'Actual file' ) ] 
    for i, ff in enumerate ( root_files ) :
        fn , fa = ff 
        row = '%-2d' % i , fn , fa if fn != fa else '' 
        rows.append ( row )
    ## 
    title = "ROOT files #%d " % len ( root_files ) 
    import ostap.logger.table as T
    table = T.table ( rows , title = title , alignment = 'cw' , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# =============================================================================
## dump all lo
if parameters :
    rows = [  ( '#' , 'Parameter' ) ] 
    for i, f in enumerate ( parameters ) :
        row = '%-2d' % i , f
        rows.append ( row )
    ## 
    title = "Parameters #%d " % len ( parameters ) 
    import ostap.logger.table as T
    table = T.table ( rows , title = title , alignment = 'cw' , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme   import docme
    docme ( __name__ , logger = logger )
    

# =============================================================================

# ============================================================================= 
    
