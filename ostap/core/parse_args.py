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
        def __init__(self            ,
                     option_strings  ,
                     dest            ,
                     nargs    =None  ,
                     const    =None  ,
                     default  =None  , 
                     type     =None  ,
                     choices  =None  ,
                     required =False ,
                     help     =None  ,
                     metavar  =None  ) :
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
            

    import ostap.core.default_config as _cnf 
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
        default = _cnf.quiet  )
    egroup1.add_argument ( 
        "--silent"     ,
        dest    = 'Silent'    , 
        action  = 'store_true' ,
        help    = "Verbose processing, same as --print-level=4 [default: %(default)s]" ,
        default = _cnf.verbose )
    egroup1.add_argument ( 
        "--verbose"     ,
        dest    = 'Verbose'    , 
        action  = 'store_true' ,
        help    = "Verbose processing, same as --print-level=1 [default: %(default)s]" ,
        default = _cnf.verbose )
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
        help    = "The commands for ``exec'' [default: %(default)s]" , 
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
        help    = "Invoke profiler" , 
        default = False             )
    # 
    parser.add_argument ( 
        '--no-mt'                     ,        
        dest    = 'NoImplicitMT'      , 
        action  = 'store_true'        , 
        help    = "DisableImplicitMT" , 
        default = False               )
    #
    parser.add_argument (
        '--config'          , 
        dest    = 'Config'  ,
        nargs   = '*'       ,
        action  = Collect   ,
        help    = "Config files to be parsed [default:  %(default)s]" ,
        default = []        , 
        )
    parser.add_argument ( 
        '--build-dir'                        ,        
        dest    = 'build_dir'                , 
        help    = "Build directory for ROOT" , 
        default = ''                         )
    #
    ## 2nd exclusive group
    group2  = parser.add_argument_group ( 'Web Display' , 'Use Web/ROOT display, see ROOT.TROOT.(Set/Get)WebDisplay') 
    egroup2 = group2.add_mutually_exclusive_group()
    egroup2.add_argument ( 
        '-w' , '--web'          ,
        dest    = 'web'         , 
        help    = "Use WebDisplay, see ROOT.TROOT.(Get/Set)WebDisplay ", 
        default = ''            ) 
    #
    egroup2.add_argument ( 
        '--no-canvas'           ,
        dest    = 'canvas'      , 
        action  = 'store_false' , 
        help    = "Do not create canvas", 
        default = True          )
    #
    ## 3rd exclusive group
    group3  = parser.add_argument_group ( 'Sesstion type' , 'General session type: interactive/embed/plain/batch...') 
    egroup3 = group3.add_mutually_exclusive_group()
    egroup3.add_argument ( '-i' ,  
                           '--interactive' , dest='batch', 
                           action = 'store_false' , default = False ,
                           help = "Interactive shell/start_ipython" )
    egroup3.add_argument ( '-e' ,
                           '--embed' , 
                           action = 'store_true' ,
                        help = "Interactive embedded shell" )
    egroup3.add_argument ( '-s' ,
                           '--simple' ,
                           action = 'store_true' ,
                           help = "Simple python shell" )
    egroup3.add_argument ( '-b' ,
                           '--batch' ,
                           action = 'store_true' , default = False , 
                           help = "Batch processing: execute files and exit" )

    if not args :
        import sys 
        args = sys.argv[1:]
    
    v = [ a for a in args ]
    if '--' in v : v.remove('--')
    
    return parser.parse_args( v )

# =============================================================================
if '__main__' == __name__ :

    from ostap.logger.logger import getLogger 
    logger = getLogger( 'ostap.core.parse_args' )

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    
# =============================================================================
##                                                                      The END 
# =============================================================================
