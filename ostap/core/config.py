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
from   ostap                     import __version__ 
from   ostap.utils.env           import ( get_env            ,
                                          has_env            ,
                                          ##
                                          boolean_true       , 
                                          boolean_false      , 
                                          ##
                                          OSTAP_CONFIG       , ## very special 
                                          OSTAP_ARGPARSE     , ## very special 
                                          OSTAP_BATCH        ,
                                          ## 
                                          OSTAP_SILENT       ,
                                          OSTAP_QUIET        ,
                                          OSTAP_DEBUG        ,
                                          OSTAP_VERBOSE      ,
                                          OSTAP_LEVEL        ,
                                          OSTAP_COLOR        ,
                                          OSTAP_UNICODE      ,
                                          ## 
                                          OSTAP_DUMP_CONFIG  ,
                                          ##
                                          OSTAP_BUILD_DIR    ,
                                          OSTAP_CACHE_DIR    ,                                          
                                          OSTAP_TMP_DIR      ,
                                          ## 
                                          OSTAP_WEB_DISPLAY  ,
                                          ##
                                          OSTAP_PARALLEL      ,                                          
                                          OSTAP_NCPUS         ,
                                          OSTAP_IMPLICITMT    ,                                          
                                          OSTAP_PROFILE       ,                                          
                                          ##
                                          OSTAP_TABLE         ,
                                          ## 
                                          OSTAP_STARTUP       )

# =============================================================================
import ostap.core.default_config      as     default_config 
import configparser, os, sys, logging, argparse, ast, copy      
# =============================================================================
## config-parser itself 
config = configparser.ConfigParser()
# =============================================================================

## define the major section:
config [ 'General' ] = {
    ##
    'ArgParse'    : str ( default_config.arg_parse    ) ,
    'Batch'       : str ( default_config.batch        ) ,
    ##
    'Silent'      : str ( default_config.silent        ) ,
    'Quiet'       : str ( default_config.quiet         ) ,
    'Debug'       : str ( default_config.debug         ) ,
    'Verbose'     : str ( default_config.verbose       ) ,
    'Level'       : str ( default_config.level         ) ,
    ##
    'Color'       : str ( default_config.color         ) , 
    'Unicode'     : str ( default_config.show_unicode  ) , 
    ##
    'DumpConfig'  : str ( default_config.dump_config   ) ,
    ##
    'BuildDir'    : str ( default_config.build_dir     ) ,
    'CacheDir'    : str ( default_config.cache_dir     ) ,
    'TmpDir'      : str ( default_config.tmp_dir       ) ,
    ##
    'WebDisplay'  : str ( default_config.webdisplay    ) , 
    ## 
    'Parallel'    : str ( default_config.parallel      ) ,
    'NCPus'       : str ( default_config.ncpus         ) ,
    'ImplicitMT'  : str ( default_config.implicitMT    ) ,
    'Profile'     : str ( default_config.profile       ) ,
    ##
    'Startup'     : str ( default_config.startup_files ) ,
    'Macros'      : str ( default_config.macros        ) ,
    'Commands'    : str ( default_config.commands      ) ,
}
# ===========================================================================
## generic TCanvas configuration
config [ 'Canvas'      ] = { 'Width'       :  '1000' , 'Height'       :  '800' , 
                             'MarginTop'   : '0.05'  , 'MarginBottom' : '0.12' ,
                             'MarginRight' : '0.07'  , 'MarginLeft'   : '0.12' }

config [ 'Fit Draw'    ] = {} ## RooFit plotting configuration
config [ 'Tables'      ] = { 'Style' : default_config.table_style } ## configuration for Tables 
config [ 'RooFit'      ] = {} ## RooFit configuration

config [ 'Pathos'      ] = {} ## PATHOS configuration  
config [ 'IPyparallel' ] = {} ## ipyparallel configuration 

# ============================================================================
## the list of config files to be processed
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

# =============================================================================
## read the files 
files_read = tuple ( config.read ( config_files ) ) 

# =============================================================================
## General section:
general      = config [ 'General' ]

# ==============================================================================
## ATTENTION: Redefien configuratuin using environment variables 
# ==============================================================================

if has_env ( OSTAP_ARGPARSE ) :
    value = get_env ( OSTAP_ARGPARSE , 'True' , silent = True  )
    value = boolean_true ( value )
    general['ArgParse' ] = 'True' if value else 'False'

if has_env ( OSTAP_BATCH ) :
    value_ = get_env ( OSTAP_BATCH , '' , silent = True  )        
    if value_ : general [ 'Batch' ] = 'True'
    
elif not general.getboolean ( 'Batch' , fallback = default_config.batch ) :
    value_ = any ( a.lower() in ( '-b' , '--batch' , '--no-gui' ) for a in sys.argv )
    if value_ : general [ 'Batch' ] = 'True'
    
# ============================================================================

if get_env ( OSTAP_SILENT , '' , silent = True ) : 

    general [ 'Silent'  ] = 'True'
    general [ 'Quiet'   ] = 'False'
    general [ 'Debug'   ] = 'False'
    general [ 'Verbose' ] = 'False'
    general [ 'Level'   ] = '%d'    % 4    

elif get_env ( OSTAP_QUIET , '' , silent = True ) :

    general [ 'Silent'  ] = 'False'
    general [ 'Quiet'   ] = 'True'
    general [ 'Debug'   ] = 'False'
    general [ 'Verbose' ] = 'False'
    general [ 'Level'   ] = '%d'    % 4    

elif get_env ( OSTAP_DEBUG , '' , silent = True ) :

    general [ 'Silent'  ] = 'False'
    general [ 'Quiet'   ] = 'False'
    general [ 'Debug'   ] = 'True'
    general [ 'Verbose' ] = 'False'
    general [ 'Level'   ] = '%d'    % 2 
    
elif get_env ( OSTAP_VERBOSE , '' , silent = True ) :

    general [ 'Silent'  ] = 'False'
    general [ 'Quiet'   ] = 'False'
    general [ 'Debug'   ] = 'False'
    general [ 'Verbose' ] = 'True'
    general [ 'Level'   ] = '%d'    % 1 

elif has_env ( OSTAP_LEVEL ) :
    
    value_ = get_env ( OSTAP_LEVEL ,  '' , silent = True )
    if value_ : general [ 'Level'   ] = value_

if has_env ( OSTAP_COLOR ) :
    value_ = get_env ( OSTAP_COLOR , '' , silent = True  )        
    if value_ : general [ 'Color'   ] = value_

if has_env ( OSTAP_UNICODE ) :
    value_ = get_env ( OSTAP_UNICODE , '' , silent = True  )        
    if value_ : general [ 'Unicode' ] = value_

# ============================================================================
## redefine ncpus from the environment variable 
if has_env ( OSTAP_NCPUS ) : # ===============================================
    # ========================================================================
    value_ = get_env ( OSTAP_NCPUS , '' , silent = True  )        
    # ========================================================================
    try : # ==================================================================
        # ====================================================================
        value_ = int ( value_ )
        if 1 <= value_ : general [ 'NCPUs' ] = str ( value_ ) 
        # ====================================================================
    except: # ================================================================
        # ====================================================================
        pass
    
# ============================================================================
## redefine parallel from the environment variable
if has_env ( OSTAP_PARALLEL ) :
    value_ = get_env ( OSTAP_PARALLEL , '' , silent = True  )        
    if value_ : general [ 'Parallel' ] = value_

if has_env ( OSTAP_PROFILE ) :
    value_ = get_env ( OSTAP_PROFILE , '' , silent = True  )        
    if value_ : general [ 'Profile' ] = value_
    
# ============================================================================
## redefine web display from the environment variable
if has_env ( OSTAP_WEB_DISPLAY ) :
    value_ = get_env ( OSTAP_WEB_DISPLAY , '' , silent = True  )        
    if value_ : general [ 'WebDisplay' ] = value_
    
if has_env ( OSTAP_BUILD_DIR ) :
    build_dir_ = get_env ( OSTAP_BUILD_DIR , '' , silent = True  )        
    if build_dir_ : general [ 'BuildDir' ] = str ( build_dir_ )
    
if has_env ( OSTAP_CACHE_DIR ) :    
    cache_dir_ = get_env ( OSTAP_CACHE_DIR , '' , silent = True  )        
    if cache_dir_ : general [ 'CacheDir' ] = str ( cache_dir_ ) 
    
if  has_env ( OSTAP_TMP_DIR ) :
    tmp_dir_ = get_env ( OSTAP_TMP_DIR , '' , silent = True  )        
    if tmp_dir_ : general [ 'TmpDir' ] = str ( tmp_dir_ )   

if  has_env ( OSTAP_DUMP_CONFIG ) :
    dump_ = get_env ( OSTAP_DUMP_CONFIG , '' , silent = True  )        
    if dump_ : general [ 'DumpConfig' ] = str ( dump_ )   

if has_env ( OSTAP_STARTUP ) : 
    value_    = get_env ( OSTAP_CONFIG , '' , silent = True )
    if value_ : general [ 'StartUp' ] = str ( tuple ( value_.split ( os.pathsep ) ) ) 
    
# =============================================================================
## Some explicit & important elements from the `General` section:
arg_parse    = general.getboolean ( 'ArgParse'    , fallback = default_config.arg_parse    )
batch        = general.getboolean ( 'Batch'       , fallback = default_config.batch        )
##
silent       = general.getboolean ( 'Silent'      , fallback = default_config.silent       )
quiet        = general.getboolean ( 'Quiet'       , fallback = default_config.quiet        )
debug        = general.getboolean ( 'Debug'       , fallback = default_config.debug        )
verbose      = general.getboolean ( 'Verbose'     , fallback = default_config.verbose      )
level        = general.getint     ( 'Level'       , fallback = default_config.level        )
color        = general.getboolean ( 'Color'       , fallback = default_config.color        )
show_unicode = general.getboolean ( 'Unicode'     , fallback = default_config.show_unicode )
##
dump_config  = general.get        ( 'DumpConfig'  , fallback = default_config.dump_config  )
##
build_dir    = general.get        ( 'BuildDir'    , fallback = default_config.build_dir    )
cache_dir    = general.get        ( 'CacheDir'    , fallback = default_config.cache_dir    )
tmp_dir      = general.get        ( 'TmpDir'      , fallback = default_config.tmp_dir      )
##
webdisplay   = general.get        ( 'WebDisplay'  , fallback = default_config.webdisplay   )
##
ncpus        = general.getint     ( 'NCPUs'       , fallback = default_config.ncpus        )
parallel     = general.get        ( 'Parallel'    , fallback = default_config.parallel     )
implicitMT   = general.getboolean ( 'ImplicitMT'  , fallback = default_config.implicitMT   )
profile      = general.getboolean ( 'Profile '    , fallback = default_config.profile      )

startup_files = general.get       ( 'StartUp' , fallback = str ( default_config.startup_files ) )
if startup_files : startup_files = ast.literal_eval ( startup_files )
else             : startup_files = [] 

macros        = general.get       ( 'Macros' , fallback = str ( default_config.macros ) )
if macros        : macros = ast.literal_eval ( macros )
else             : macros = []

commands      = general.get       ( 'Commands'  , fallback = str ( default_config.commands ) )
if commands      : commands = ast.literal_eval ( commands )
else             : commands = [] 

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

if has_env ( OSTAP_TABLE ) :
    value_ = get_env ( OSTAP_TABLE , '' , silent = True  )        
    if value_ : tables [ 'Style' ] = value_
    
# =============================================================================
## The final action "at-exit"
#  - print the list/table of read files  
#  - dump configuration into the output file 
def config_atexit ( config , files ) :
    """ The final `at-exit' action:
    - print the list/table of read files  
    - dump configuration into the output file 
    """
    import  datetime, sys
    if not hasattr ( sys , 'ps1' ) : return
    ##
    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.core.config' )
    ## 
    now = datetime.datetime.now()
    ##
    # ===========================================================================
    ## (1) print lift/table of read input config files 
    # ===========================================================================
    if 1 == len ( files ) : logger.info ( 'Ostap configuration is read from: %s' % files [ 0 ] )         
    elif files : 
        import ostap.logger.table as T
        rows = [ ( '', 'file' ) ]
        for i, f in enumerate ( files , start = 1 ) :
            row = '%d' % i , f 
            rows.append ( row )            
        title = 'Ostap configuration files'
        table = T.table ( rows , title = title , prefix = '# ' , alignment = 'rl' )
        logger.info ( 'Ostap configuration is read from\n%s' % table ) 

    # ===========================================================================
    ## (2) print lift/table of read input config files 
    # ===========================================================================
    import io 
    with io.StringIO() as o : 
        config.write( o )
        logger.verbose ( 'Ostap configuration:\n%s' % o.getvalue() )

    # ===========================================================================
    ## get the dumpfile
    import ostap.core.default_config as dc 
    dump = config [ 'General' ].get ( 'DumpConfig' , fallback = dc.dump_config )
    # ===========================================================================
    try : # =====================================================================
        # =======================================================================
        if os.path.exists ( dump ) : os.remove ( dump ) # =======================
        the_line = '# ' + 78 * '*' + '\n'
        # =======================================================================
        with open ( dump , 'wt' ) as fdump : # ==================================
            # ===================================================================
            fdump.write( the_line )
            if files : 
                fdump.write('# Ostap configuration read from:\n' )
                for i,f in enumerate ( files , start = 1 ) : fdump.write( '# %-3d %s\n' % ( i , f ) ) 
            fdump .write ('# Ostap configuration:\n' )                
            fdump .write ( the_line )    
            config.write ( fdump )            
            fdump .write ( the_line )
            fdump .write ( '# Configuration saved at %s\n' % now.strftime ( '%c' ) )
            fdump .write ( the_line )
        # ========================================================================
        if os.path.exists ( dump ) and os.path.isfile ( dump ) :
            logger.info ( 'Ostap  configuration saved to %s' %  dump )
        # ========================================================================
    except :
        pass

# ================================================================================
## register the final action
import atexit 
atexit.register ( config_atexit , config = config , files = files_read ) 

# ================================================================================
## Print config parser as table
def _cp_table_ ( parser , files = files_read , title = '' , prefix = '' ) :
    """ Print config parser as table
    """
    if not files and files_read : files = files_read 

    rows = [ ( 'Section' , 'Option' , 'Value' ) ]
    if files : 
        row  = 'Config files' , '' , '\n'.join ( f for f in files )
        rows.append ( row ) 
    for section in parser :
        row = '[%s]' % section , '' , '' 
        rows.append ( row )
        the_section = config [ section ] 
        options = sorted ( the_section ) 
        for option in options : 
            row = '' , option , the_section [ option ] 
            rows.append ( row )
            
    title = title if title else 'Ostap configuration'
    import ostap.logger.table as T
    return T.table ( rows , title = title , prefix = prefix , alignment = 'llw' )

# ==============================================================================
type ( config ).table    = _cp_table_
type ( config ).__str__  = _cp_table_
type ( config ).__repr__ = _cp_table_

# ==============================================================================
## Command-Line argumens
# ==============================================================================

# ==============================================================================
## Parse comand line arguments
def __parse_args ( args  = [] ) :
    """ Parse command line arguments
    """
    # =========================================================================
    ## @class Collect
    #  simple parsing action to collect multiple arguments
    #  @code
    #  parser =...
    #  parser.add_argument('--foo',
    #  ... action  = Collect ,
    #  ... nargs   = '*'     ,
    #  ... default = []      ,
    #  ...)
    #  print(parser.parse_args('a.txt b.txt --foo 1 2 3 --foo 4 -foo 5 '.split()))
    #  @endcode
    class Collect(argparse.Action):
        """ Simple parsing action to collect multiple arguments
        >>> parser =...
        >>> parser.add_argument('--foo',
        ... action  = Collect ,
        ... nargs   = '*'     ,
        ... default = []      ,
        ...)
        >>> print(parser.parse_args('a.txt b.txt --foo 1 2 3 --foo 4 -foo 5 '.split()))
        """
        def __init__(self             ,
                     option_strings   ,
                     dest             ,
                     nargs    = None  ,
                     const    = None  ,
                     default  = None  , 
                     type     = None  ,
                     choices  = None  ,
                     required = False ,
                     help     = None  ,
                     metavar  = None  ) :
            if nargs == 0:
                raise ValueError('nargs for Collect actions must be > 0; if arg '
                                 'strings are not supplying the value to append, '
                                 'the append const action may be more appropriate')
            if const is not None and nargs != argparse.OPTIONAL:
                raise ValueError('nargs must be %r to supply const' % argparse.OPTIONAL)
            super(Collect, self).__init__(
                option_strings = option_strings ,
                dest           = dest           ,
                nargs          = nargs          ,
                const          = const          ,
                default        = default        ,
                type           = type           ,
                choices        = choices        ,
                required       = required       ,
                help           = help           ,
                metavar        = metavar        )
            
        def __call__(self, parser, namespace, values, option_string=None):
            if None == getattr   ( namespace  , self.dest , None ) : setattr ( namespace , self.dest , []  )
            items   =  getattr   ( namespace , self.dest ) 
            items   =  copy.copy ( items ) 
                
            ## items = argparse._copy.copy(argparse._ensure_value(namespace, self.dest, []))
            
            ## the only one important line: 
            for v in values : items.append ( v )
            setattr( namespace, self.dest, items )
            
    from argparse import ArgumentParser 
    parser = ArgumentParser ( prog = 'ostap' )
    ## 
    parser.add_argument ( 
        '-v'    , '--version' ,
        action  = 'version'   , 
        version = 'Ostap %s' % __version__ )  

    ## 1st exclusive group
    group1  = parser.add_argument_group ( 'Control verbosity' , 'Control the general level of verbosity') 
    egroup1 = group1.add_mutually_exclusive_group()    
    egroup1.add_argument ( 
        "-q" , "--quiet"       ,
        dest    = 'Quiet'      , 
        action  = 'store_true' ,
        help    = "Quite processing, same as --print-level=4 [default: %(default)s]" ,
        default = quiet        )
    
    egroup1.add_argument ( 
        "--silent"             ,
        dest    = 'Silent'     , 
        action  = 'store_true' ,
        help    = "Verbose processing, same as --print-level=5 [default: %(default)s]" ,
        default = silent       )
    
    egroup1.add_argument ( 
        "--verbose"            ,
        dest    = 'Verbose'    , 
        action  = 'store_true' ,
        help    = "Verbose processing, same as --print-level=1 [default: %(default)s]" ,
        default = verbose      )
    
    egroup1.add_argument ( 
        "--debug"              ,
        dest    = 'Debug'      , 
        action  = 'store_true' ,
        help    = "Debug processing, same as --print-level=2 [default: %(default)s]" ,
        default = debug        )
    
    egroup1.add_argument ( 
        "-p" , "--print-level"     ,
        metavar = "LEVEL"          ,
        dest    = 'Level'          ,
        choices = range ( -1 , 9 ) , 
        type    = int              , 
        help    =  "Printout level [default: %(default)s]" ,
        default = level if 0 <= level <= 8 else -1 )    
    ##
    group2  = parser.add_argument_group ( 'Files' , '(ROOT) files, ROOT/C++ macros, (Python) scripts and commands for processing') 
    ###
    group2.add_argument ( 
        '--startup'     ,
        metavar = "STARTUPS"      ,
        dest    = 'StartUp'       ,
        nargs   = '*'             ,
        action  = Collect         , 
        help    = "List of startup files to be (pre)executed [default: %(default)s]" , 
        default = startup_files   )
    ##
    group2.add_argument ( 
        "-m"    , "--macros"      ,
        metavar = "MACROS"        ,
        dest    = 'Macros'        ,
        nargs   = '*'             ,
        action  = Collect         , 
        help    = "ROOT/C++ macros to be processed [default: %(default)s]" ,
        default = macros          )
    ###
    group2.add_argument ( 
        '-c'    , '--command'     ,
        metavar = "COMMANDS"      ,
        dest    = 'Commands'      ,
        nargs   = '*'             ,
        action  = Collect         , 
        help    = "The commands to be executed [default: %(default)s]" , 
        default = commands         )
    ## 
    group2.add_argument (
        "files"           ,
        metavar = "FILES" ,
        nargs   = '*'     , 
        help    = "ROOT files, python macros, python files to be processed one by one [default: %(default)s]" ,
        default = []      )
    ##
    group3  = parser.add_argument_group ( 'CPU/processes/parallelism' , 'Options for parallel processing') 
    group3.add_argument (
        '-n' , '--ncpus'          , 
        metavar = "NCPUS"         ,
        dest    = 'NCPUs'         ,
        type    = int             ,
        help    = 'Maximal number of CPUs [default: %(default)s]' , 
        default = ncpus           )
    ## 
    group3.add_argument (
        '--parallel'              , 
        metavar = "PARALELL"      ,
        dest    = 'Parallel'      ,
        help    = 'Machinery for parallel processing [default: %(default)s]' , 
        default = parallel        )
    ##
    group3.add_argument ( 
        '--no-mt'                ,        
        dest    = 'ImplicitMT'   , 
        action  = 'store_false'  , 
        help    = "EnableImplicitMT? [default: %(default)s" , 
        default = implicitMT     )
    ## 
    ## 4nd exclusive group
    group4  = parser.add_argument_group ( 'Web Display' , 'Use Web/ROOT display, see ROOT.TROOT.(Set/Get)WebDisplay') 
    egroup4 = group4.add_mutually_exclusive_group()
    egroup4.add_argument ( 
        '-w' , '--web'         ,
        metavar = "WEBDISPLAY" ,        
        dest    = 'WebDisplay' , 
        help    = "Use WebDisplay, see ROOT.TROOT.(Get/Set)WebDisplay ", 
        default = webdisplay   )   
    #
    egroup4.add_argument ( 
        '--no-canvas'           ,
        dest    = 'Canvas'      , 
        action  = 'store_false' , 
        help    = "Create global canvas [default: %(default)s", 
        default = True          )

    ## 5rd exclusive group
    group5  = parser.add_argument_group ( 'Session type' , 'General session type: interactive/embed/plain/batch...') 
    egroup5 = group5.add_mutually_exclusive_group()
    
    egroup5.add_argument (
        '-i' ,  '--interactive' ,
        dest    = 'batch'       , 
        action  = 'store_false' ,
        help    = "Interactive shell/start_ipython" ,        
        default =  batch         )
    
    egroup5.add_argument (
        '-e'    , '--embed'     , 
        action  = 'store_true'  ,
        help    = "Interactive embedded shell [default: %(default)s]"  , 
        default = False         )
        
    egroup5.add_argument  (
        '-s' , '--simple'      ,
        action = 'store_true'  ,
        help   = "Simple/trivial python shell [default: %(default)s] " , 
        default = False        )
        
    egroup5.add_argument (
        '-b' , '--batch'       ,
        action  = 'store_true' ,
        help    = "Batch processing: execute files and exit [default: %(default)s]" , 
        default = batch       ) 
    ## 
    group6 = parser.add_argument_group ( 'Directories' ,
                                         'Various directories for ROOT&Ostap') 
    group6.add_argument ( 
        '--build-dir'         ,
        metavar = "BUILD_DIR" ,        
        dest    = 'BuildDir'  , 
        help    = "Build directory for ROOT&Ostap [default: %(default)s]"     , 
        default = build_dir   )
    ##
    group6.add_argument ( 
        '--cache-dir'         ,         
        metavar = "CACHE_DIR" ,        
        dest    = 'CacheDir'  , 
        help    = "Cache directory for Ostap [default: %(default)s]"     , 
        default = cache_dir   )
    ##
    group6.add_argument ( 
        '--tmp-dir'           ,       
        metavar = "TMP_DIR"   ,        
        dest    = 'TmpDir'    ,
        help    = "Top-level temporary directory for Ostap [default: %(default)s]" ,
        default = tmp_dir     )

    ## 
    group7 = parser.add_argument_group ( 'Miscellaneous' ,
                                         'Various miscelalneous options ROOT&Ostap') 
    ## 
    group7.add_argument ( 
        '--no-color'              ,
        dest    = 'Color'         , 
        action  = 'store_false'   , 
        help    = 'Use colorization? [default: %(default)s]', 
        default = color           )    
    #
    group7.add_argument (
        '--unicode'               ,
        action  = "store_true"    ,
        dest    = 'Unicode'       ,
        help    = 'Use unicode in log-files? [default: %(default)s]',
        default = show_unicode    )
    #
    group7.add_argument ( 
        '--profile'               ,
        dest    = 'Profile'       , 
        action  = 'store_true'    , 
        help    = "Invoke profiler? [default: %(default)s]" , 
        default = profile         )
    #
    group7.add_argument ( 
        '--dump'                   ,
        metavar = "DUMP_FILE"      ,        
        dest    = 'DumpConfig'     , 
        help    = "Dump-file for configuration [default: %(default)s]" , 
        default = dump_config      )
    #
    
    # ===============================================================================
    ## use the parser!!
    # ===============================================================================
    if not args :
        import sys 
        args = sys.argv[1:]

    v = [ a for a in args ]
    if '--' in v : v.remove('--')

    ## ATTENTION !! 
    if not arg_parse : v = []  ## ATENTION! 

    return parser.parse_args( v )

# ================================================================================
## command-line arguments
arguments = __parse_args ()
        
# ================================================================================
## update the global configuration
general [ 'Batch'      ] = str ( arguments.batch )
## 
general [ 'Silent'     ] = str ( arguments.Silent     ) 
general [ 'Quiet'      ] = str ( arguments.Quiet      ) 
general [ 'Verbose'    ] = str ( arguments.Verbose    )
general [ 'Debug'      ] = str ( arguments.Debug      )
general [ 'Level'      ] = str ( arguments.Level      )
general [ 'Color'      ] = str ( arguments.Color      )
general [ 'Unicode'    ] = str ( arguments.Unicode    )
##
general [ 'DumpConfig' ] = str ( arguments.DumpConfig )
##
general [ 'BuildDir'    ] = arguments.BuildDir
general [ 'CacheDir'    ] = arguments.CacheDir
general [ 'TmpDir'      ] = arguments.TmpDir 
## 
general [ 'WebDisplay'  ] = str ( arguments.WebDisplay )
##
general [ 'NCPUs'       ] = str ( arguments.NCPUs      )
general [ 'Parallel'    ] = str ( arguments.Parallel   )
general [ 'ImplicitMT'  ] = str ( arguments.ImplicitMT )
general [ 'Profile'     ] = str ( arguments.Profile    )
general [ 'StartUp'     ] = str ( arguments.StartUp    )
general [ 'Macros'      ] = str ( arguments.Macros     )
general [ 'Commands'    ] = str ( arguments.Commands   )

## 
batch         = arguments.batch
##
silent        = arguments.Silent 
quiet         = arguments.Quiet
debug         = arguments.Debug
verbose       = arguments.Verbose
level         = arguments.Level
color         = arguments.Color
show_unicode  = arguments.Unicode
##
dump_config   = arguments.DumpConfig 
##     
build_dir     = arguments.BuildDir
cache_dir     = arguments.CacheDir
tmp_dir       = arguments.TmpDir
##
webdisplay    = arguments.WebDisplay
## 
parallel      = arguments.Parallel
ncpus         = arguments.NCPUs
implicitMT    = arguments.ImplicitMT
profile       = arguments.Profile 
startup_files = arguments.StartUp
macros        = arguments.Macros 
commands      = arguments.Commands
input_files   = arguments.files 

# ==============================================================================
if '__main__' == __name__ :

    # ==========================================================================
    # logging 
    # ==========================================================================
    from ostap.logger.logger import getLogger 
    logger = getLogger ( 'ostap.core.config' )

    ## needed for docme
    def _cp_hash_ ( cp ) : return hash ( str ( cp ) )
    type ( config ) .__hash__ = _cp_hash_

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    table = config.table ( files_read , prefix = '# ') 
    logger.info ( 'Ostap configuration is:\n%s' % table ) 
    
    # ========================================================================
    rows  = [ ( 'Argument' , 'Value' ) ]
    vars_ = vars ( arguments )
    for key in sorted (vars_.keys ()  ) :
        row = key , str ( vars_ [ key ] ) 
        rows.append ( row )
        
    title = 'Ostap command line arguments'
    import ostap.logger.table as T
    rows  = T.table ( rows , title = title , prefix = '# ' , alignment = 'rw' )
    logger.info ( '%s:\n%s' % ( title , rows ) ) 

# =============================================================================
##                                                                      The END 
# =============================================================================
