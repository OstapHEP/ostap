#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/parse_args.py
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
""" parse ostap arguments   
"""
# =============================================================================
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = "2012-09-10"
__version__ = '$Revision$'
__all__     = (
    'parse_args' ,  ## parse ostap arguments 
)
# =============================================================================
from   ostap     import __version__ 
## ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse 
# =============================================================================
## parse arguments 
def parse_args ( args = [] ) :
    """ Parse arguments 
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
            items = argparse._copy.copy(argparse._ensure_value(namespace, self.dest, []))
            ## the only one important line: 
            for v in values : items.append(v)
            setattr(namespace, self.dest, items)
            
    import ostap.core.config as config 
    
    from argparse import ArgumentParser 
    parser = ArgumentParser ( prog = 'ostap' )
    #
    ## 1st exclusive group
    group1  = parser.add_argument_group ( 'Control verbosity' , 'Control the general level of verbosity') 
    egroup1 = group1.add_mutually_exclusive_group()    
    egroup1.add_argument ( 
        "-q" , "--quiet"       ,
        dest    = 'Quiet'      , 
        action  = 'store_true' ,
        help    = "Quite processing, same as --print-level=5 [default: %(default)s]" ,
        default = config.quiet  )
    egroup1.add_argument ( 
        "--silent"     ,
        dest    = 'Silent'    , 
        action  = 'store_true' ,
        help    = "Verbose processing, same as --print-level=4 [default: %(default)s]" ,
        default = config.quiet )
    egroup1.add_argument ( 
        "--verbose"     ,
        dest    = 'Verbose'    , 
        action  = 'store_true' ,
        help    = "Verbose processing, same as --print-level=1 [default: %(default)s]" ,
        default = config.verbose )
    egroup1.add_argument ( 
        "--debug"        ,
        dest    = 'Debug'    , 
        action  = 'store_true' ,
        help    = "Debug processing, same as --print-level=2 [default: %(default)s]" ,
        default = False )
    egroup1.add_argument ( 
        "-p" , "--print-level"                  ,
        dest    = 'Level'          ,
        choices = range ( -1 , 8 ) , 
        type    = int              , 
        help    =  "Printout level [default: %(default)s]" ,
        default = -1       )    
    #
    parser.add_argument (
        "files" ,
        metavar = "FILE"  ,
        nargs   = '*'     , 
        help    = "ROOT/python/macro files to be opened/processed [default: %(default)s]" ,
        default = []  )
    #
    parser.add_argument ( 
        '-v'    , '--version'     ,
        dest    = 'Commands'      ,
        action  = 'version'       , 
        version = 'Ostap %s' % __version__ )  
    #
    parser.add_argument ( 
        '-c'    , '--command'     ,
        dest    = 'Commands'      ,
        nargs   = '*'             ,
        action  = Collect         , 
        help    = "The commands for `exec' [default: %(default)s]" , 
        default = []              )
    #
    parser.add_argument ( 
        "-m"    , "--macros"      ,
        metavar = "MACROS"        ,
        dest    = 'Macros'        ,
        nargs   = '*'             ,
        action  = Collect         , 
        help    = "ROOT macros to be loaded [default: %(default)s]",
        default = []  )
    #
    parser.add_argument (
        '--no-context'           ,
        action  = "store_false"  ,
        dest    = 'WithContext'  ,
        help    = "Do not use global Ostap context for the scripts",
        default = True           )
    #
    parser.add_argument ( 
        '--no-color'     ,
        dest    = 'Color'      , 
        action  = 'store_false' , 
        help    = "Do not use colorization", 
        default = True          )    
    #
    parser.add_argument ( 
        '--profile'     ,
        dest    = 'Profile'         , 
        action  = 'store_true'      , 
        help    = "Invoke profiler [default: %(default)s]" , 
        default = False             )
    # 
    group2  = parser.add_argument_group ( 'CPU/processes/parallelism' , 'Options for parallel processing') 
    group2.add_argument (
        '-n' , '--ncpus'          , 
        dest    = 'NCPUs'         ,
        type    = int             ,
        help    = 'Maximal number of CPUs [default: %(default)s]' , 
        default = config.ncpus    )
         
    group2.add_argument (
        '--parallel'         , 
        dest    = 'Parallel' ,
        help    = 'Machinery for parallel processing [default: %(default)s]' , 
        default = config.parallel  ) 
    group2.add_argument ( 
        '--no-mt'                     ,        
        dest    = 'NoImplicitMT'      , 
        action  = 'store_true'        , 
        help    = "DisableImplicitMT [default: %(default)s]" , 
        default = False               )
    
    #
    ## 3nd exclusive group
    group3  = parser.add_argument_group ( 'Web Display' , 'Use Web/ROOT display, see ROOT.TROOT.(Set/Get)WebDisplay') 
    egroup3 = group3.add_mutually_exclusive_group()
    egroup3.add_argument ( 
        '-w' , '--web'          ,
        dest    = 'web'         , 
        help    = "Use WebDisplay, see ROOT.TROOT.(Get/Set)WebDisplay ", 
        default = config.webdisplay  )   
    #
    egroup3.add_argument ( 
        '--no-canvas'           ,
        dest    = 'canvas'      , 
        action  = 'store_false' , 
        help    = "Do not create canvas", 
        default = True          )
    #
    ## 3rd exclusive group
    group4  = parser.add_argument_group ( 'Session type' , 'General session type: interactive/embed/plain/batch...') 
    egroup4 = group4.add_mutually_exclusive_group()
    egroup4.add_argument ( '-i' ,  
                           '--interactive' , dest='batch', 
                           action = 'store_false' , default = False ,
                           help = "Interactive shell/start_ipython" )
    egroup4.add_argument ( '-e' ,
                           '--embed' , 
                           action = 'store_true' ,
                           help = "Interactive embedded shell" )
    egroup4.add_argument ( '-s' ,
                           '--simple' ,
                           action = 'store_true' ,
                           help = "Simple python shell" )
    egroup4.add_argument ( '-b' ,
                           '--batch' ,
                           action = 'store_true' , default = config.batch , 
                           help = "Batch processing: execute files and exit" )

    group5 = parser.add_argument_group ( 'Directries' , 'Various directories for Ostap') 
    group5.add_argument ( 
        '--build-dir'                             ,         
        dest    = 'BuildDir'                      , 
        help    = "Build directory for ROOT&Ostap [default: %(default)s]"     , 
        default = config.build_dir                )
    group5.add_argument ( 
        '--cache-dir'                             ,         
        dest    = 'CacheDir'                      , 
        help    = "Cache directory for Ostap [default: %(default)s]"     , 
        default = config.cache_dir                )
    group5.add_argument ( 
        '--tmp-dir'                               ,       
        dest    = 'TmpDir'                        ,
        help    = "Temporary directory for Ostap [default: %(default)s]" ,
        default = config.tmp_dir                  )
    
    
    if not args :
        import sys 
        args = sys.argv[1:]
    
    v = [ a for a in args ]
    if '--' in v : v.remove('--')
    
    arguments = parser.parse_args( v )

    # =============================================================================
    
    config.general ['Quiet'      ] = str ( arguments.Quiet    ) 
    config.general ['Verbose'    ] = str ( arguments.Verbose  )
    config.general ['WebDisplay' ] = str ( arguments.web      ) 
    config.general ['Batch'      ] = str ( arguments.batch    )
    config.general ['Parallel'   ] = str ( arguments.Parallel ) 
    config.general ['NCPUs'      ] = str ( arguents.NCPUs     )
    
    config.general ['BuildDir'   ] = arguments.BuildDir
    config.general ['CacheDir'   ] = arguments.CacheDir
    config.general ['TmpDir'     ] = arguments.TmpDir 
    
    config.quiet      = arguments.Quiet
    config.verbose    = arguments.Verbose
    config.batch      = arguments.batch    
    config.webdisplay = arguments.web
    config.parallel   = arguments.Parallel
    config.ncpus      = arguments.NCPUs
    
    config.build_dir  = arguments.BuildDir
    config.cache_dir  = arguments.CacheDir
    config.tmp_dir    = arguments.TmpDir 


    return arguments


# =============================================================================
if '__main__' == __name__ :

    from ostap.logger.logger import getLogger 
    logger = getLogger( 'ostap.core.parse_args' )

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
