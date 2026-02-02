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
from ostap.core.ostap_setup import * 

# =============================================================================
## import everything from ostap
# =============================================================================
if config.quiet or config.silent : # ==========================================
    # =========================================================================
    from ostap.logger.utils import mute
    with mute () : # ==========================================================
        from ostap.core.load_ostap import * 
    del mute
    # =========================================================================
else : # ======================================================================
    # =========================================================================
    from ostap.core.load_ostap import * 

# =============================================================================
## create default canvas
# =============================================================================
if config.arguments.Canvas :
    import ostap.plotting.canvas 
    logger.debug ( "Create default Ostap canvas" )
    canvas = ostap.plotting.canvas.getCanvas ()

# =============================================================================
## suppress excessive (?) RooFit printout
# =============================================================================
if ( config.quiet or config.silent ) and not config.verbose and 3 < config.level :
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
## make `reload` command available 
from importlib import reload 

# =============================================================================
if '__main__' == __name__ : logger.info ( 'ostap is ready'  ) 
else                      : logger.info ( 'ostap is loaded' )

# =============================================================================
##                                                                      The END 
# =============================================================================

