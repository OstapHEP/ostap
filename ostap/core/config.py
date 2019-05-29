#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/config.py
#  The basic configuration of ostap 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-05-19
# =============================================================================
"""The basic configuration of ostap 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2019-05-19"
__all__     = (
    )
# =============================================================================
import configparser

config = configparser.ConfigParser()

## Define the major sections
config [ 'General'  ] = {
    'Quiet'   : 'False' ,
    'Verbose' : 'False' ,
    }

config [ 'Canvas'   ] = { 'Width' :  '1000' , 'Height' :  '800' } 
config [ 'Fit Draw' ] = {}
config [ 'Parallel' ] = {}
config [ 'Style'    ] = {}

read    = config.read ( [ '~/.ostaprc'               ,
                          '~/.config/ostap/.ostaprc' ,
                          '.ostaprc '                ] )

# =============================================================================
## sections
general = config [ 'General' ]

quiet   = general.getboolean ( 'Quiet'  , fallback = False )
verbose = general.getboolean ( 'Verbose', fallback = False )

# =============================================================================
canvas  = config [ 'Canvas' ]


# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.core.config' )
else                       : logger = getLogger( __name__     )
# =============================================================================
import logging
logging.disable ( ( logging.WARNING - 1 ) if quiet   else
                  ( logging.DEBUG   - 5 ) if verbose else ( logging.INFO - 1 ) )

# =============================================================================
logger.info  ( 'The basic configuration of Ostap: %s' %   read )

import io 
with io.StringIO() as o : 
    config.write( o )
    logger.debug ( 'The basic configuration of Ostap:\n %s' % o.getvalue() )
del o, io
    
# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
