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
    'executed_macros'  , ## list of successfully executed ROOT/C++ macros 
    'root_files'       , ## list of ROOT files  
    'parameters'       , ## list of extra comand-line arguments 
)
# ============================================================================= 
import ostap.core.config      as     config 
import ostap.fixes.fixes
import ostap.core.build_dir
import ostap.core.cache_dir
import ostap.utils.cleanup 
import ROOT, math, os, sys   
# =============================================================================
from ostap.logger.logger import getLogger
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
    if     groot.IsBatch() : logger.attention ( "BATCH processig is activated!"   )
elif groot and not config.batch and     groot.IsBatch () :
    groot.SetBatch ( False )
    if not groot.IsBatch() : logger.info      ( "BATCH processig is deactivated!" )

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
    d1 = ROOT.kFatal - ROOT.Print
    d2 = ROOT.RooFit.FATAL - ROOT.RooFit.DEBUG 
    l1 = ROOT.kFatal       + math.floot ( ( d1 * level ) / 8 )
    l2 = ROOT.RooFit.DEBUG + math.floor ( ( d2 * level ) / 8 )
    if groot : groot.ProcessLine ( "gErrorIgnoreLevel= %d ; " % l1 )     
    msg.setGlobalKillBelow  ( l2 )
    msg.setSilentMode       ( ROOT.RooFit.INFO < l2 )

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
    # =========================================================================
    ## execute it!
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        logger.debug ( "Execute the script: `%s` " % _su ) 
        import runpy
        globs = runpy.run_path ( _ss , globals() ) ## , run_name = '__main__' )
        globs = dict ( ( ( k , v ) for k , v in globs.items() if  not k.startswith ( '__' ) and not k.endswith ( '__' ) ) )
        logger.debug ( 'Symbols from %s: %s' % ( _su , globs.keys() ) )
        globals().update( globs )
        del globs
        ##
        _scripts.add ( _ss )

        if os.path.exists ( _ss ) and os.path.isfile ( _ss ) : _ss = os.path.abspath ( _ss ) 
        executed_scripts.append ( ( _su , _ss ) )
        
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
executed_macros =     []
_macros         = set () 
for _su in config.macros :
    ## 
    _ss = _su
    _ss = os.path.expandvars ( _ss )
    _ss = os.path.expandvars ( _ss )
    _ss = os.path.expandvars ( _ss )
    _ss = os.path.expanduser ( _ss )
    _ss = os.path.expandvars ( _ss )
    ##
    if    _ss.endswith  ( '++' ) : _mm = _ss [ : 2 ]
    elif  _ss.endswith  ( '+'  ) : _mm = _ss [ : 1 ]
    else                         : _mm = _ss 
    ##     
    if not os.path.exists     ( _mm ) : continue
    if not os.path.isfile     ( _mm ) : continue
    _mm =  os.path.abspath    ( _mm )
    ##
    if   _ss.endswith ( '++' ) : _mm = '%s++' % _mm
    elif _ss.endswith ( '+'  ) : _mm = '%s+'  % _mm
    ##
    if _mm in _macros : continue
    ## 
    logger.debug ( "Execute the macro: `%s` " % _su ) 
    groot = ROOT.ROOT.GetROOT()
    sc = groot.LoadMacro ( _mm )
    ## 
    if sc and config.batch : 
        raise RuntimeError  ( 'Failure to execute macro "%s", code:%d' % ( _su , sc ) )
    elif sc  : logger.error ( 'Failure to execute macro "%s", code:%d' % ( _su , sc ) )
    ## 
    _macros.add ( _mm )
    ##
    executed_macros.append (  ( _su , _mm ) ) 
    logger.debug ( "Macro '%s' is executed"      % _su )
    
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

    fname = _ss 
    fbase , dot , ext = fname.partition('.')
    
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
            
    ## ROOT/cpp macro ?
    elif fbase and dot and ext.upper() in cpp_ext :
        
        if   fname.endswith ( '++' ) : sname = fname [:-2] 
        elif fname.endswith ( '+'  ) : sname = fname [:-1] 
        else                         : sname = fname

        if os.path.exists ( sname ) and os.path.isfile ( sname ) :
            
            logger.debug ( "Execute the macro: `%s` " % _su ) 
            groot = ROOT.ROOT.GetROOT()
            sc = groot.LoadMacro ( fname )
            ## 
            if sc and config.batch : 
                raise RuntimeError  ( 'Failure to execute macro "%s", code:%d' % ( _su , sc ) )
            elif sc  : logger.error ( 'Failure to execute macro "%s", code:%d' % ( _su , sc ) )

            ## 
            logger.debug ( "Macro '%s' is executed"      % _su )
            
            fn = sname
            if os.path.exists ( fn ) and os.path.isfile ( fn ) : fn = os.path.abspath ( fn )            
            _macros.add ( ( _su , fn ) )
            
            executed_macros.append ( _su ) 

    ## python/ostap script ?
    elif fbase and dot and ext.upper() in ( 'PY' , 'OST' , 'OSTP' , 'OSTAP' , 'OPY' ) :

        # =====================================================================
        try : # ===============================================================
            # =================================================================
            
            import runpy
            globs = runpy.run_path ( fname , globals() ) ## , run_name = '__main__')
            globs = dict ( ( ( k , v ) for k,v in globs.items() if  not k.startswith('__') and not k.endswith('__')) )
            logger.debug ( 'Symbols from %s: %s' % ( _su , globs.keys() ) ) 
            globals().update ( globs )
            del globs 
            
            executed_scripts.append ( _su )
            logger.debug  ( "Python script '%s' is executed"      % _su )
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
if executed_macros :
    rows = [  ( '#' , 'Macro' , 'Expanded name' ) ] 
    for i, fd in enumerate ( executed_macros , start = 1 ) :
        fn , fa = ff 
        row = '%-2d' % i , fn , fa if fn != fa else '' 
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
    
