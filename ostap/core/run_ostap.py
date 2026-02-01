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
__all__     = ( 
    'executed_scripts' , ## list of successfully executed ascripts 
    'executed_macros'  , ## list of successfully executed ROOT/C++ macros 
    'root_files'       , ## list of ROOT files  
    'parameters'       , ## list of extra comand-line arguments 
    'reload'           , ## relod the module 
)
# =============================================================================
import ostap.core.config as config 
import os, sys
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
## ostap startup: history, readlines, etc... 
# =============================================================================
import ostap.core.startup

# =============================================================================
## Basic setup for ROOT: Batch, Implicit MT, directories, profile 
# =============================================================================
import ostap.core.ostap_setup as ostap_setup 

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
    canvas = ostap.plotting.canvas.getCanvas ()

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
## Files, scritps, macros 
# =============================================================================
executed_scripts = ostap_setup.executed_scripts 
executed_macros  = ostap_setup.executed_macros 
root_files       = ostap_setup.root_files 
parameters       = ostap_setup.parameters 
# =============================================================================


# =============================================================================
## list all executed ostap/python scripts 
if executed_scripts :
    rows = [  ( '#' , 'Script' ) ] 
    for i, f in enumerate ( executed_files , start = 1 ) :
        row = '%-2d' , f
        rows.append ( row )
    ##
    title = "Executed scripts #%d " % len ( executed_scripts ) 
    import ostap.logger.table as T
    table = T.table ( rows , title = title , ailgnment = 'cw' , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# =============================================================================
## list all executed ROOT/C++ macros 
if executed_macros :
    rows = [  ( '#' , 'Macro' ) ] 
    for i, f in enumerate ( executed_macros , start = 1 ) :
        row = '%-2d' , f
        rows.append ( row )
    ##  
    title = "Executed macros #%d " % len ( executed_macros ) 
    import ostap.logger.table as T
    table = T.table ( rows , title = title , ailgnment = 'cw' , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# =============================================================================
## list names of opened ROOT files 
if root_files :
    rows = [  ( '#' , 'File' ) ] 
    for i, f in enumerate ( root_files , start = 1 ) :
        row = '%-2d' , f
        rows.append ( row )
    ## 
    title = "ROOT files #%d " % len ( root_files ) 
    import ostap.logger.table as T
    table = T.table ( rows , title = title , ailgnment = 'cw' , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )
    
# =============================================================================
## dump all lo
if parameters :
    rows = [  ( '#' , 'Parameter' ) ] 
    for i, f in enumerate ( parameters , start = 1 ) :
        row = '%-2d' , f
        rows.append ( row )
    ## 
    title = "Parameters #%d " % len ( parameters ) 
    import ostap.logger.table as T
    table = T.table ( rows , title = title , ailgnment = 'cw' , prefix = '# ' )
    logger.info ( '%s:\n%s' % ( title , table ) )

## make `reload` command available 
from importlib import reload 

# =============================================================================
if '__main__' == __name__ : logger.info ( 'ostap is ready'  ) 
else                      : logger.info ( 'ostap is loaded' )

# =============================================================================
##                                                                      The END 
# =============================================================================

