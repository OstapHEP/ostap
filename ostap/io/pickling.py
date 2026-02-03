#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file pickling.py
# Helper module to define pickling for various databases 
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2010-04-30
# =============================================================================
""" Helper module to define pickling for various databases
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2010-04-30"
__version__ = "$Revision$" 
# =============================================================================
__all__ = (
    'DEFAULT_PROTOCOL' ,
    'HIGHEST_PROTOCOL' ,
    'PROTOCOL'         ,
    'Pickler'          , 
    'Unpickler'        ,
    'BytesIO'          ,
    'dumps'            ,
    'loads'            ,
    )
# =============================================================================
from   pickle           import ( Pickler          ,
                                 Unpickler        , 
                                 DEFAULT_PROTOCOL ,
                                 HIGHEST_PROTOCOL ,
                                 dumps , loads    )  
from   io                import BytesIO 
import ostap.core.config as config 
import sys, array 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.pickling' )
else                      : logger = getLogger ( __name__               )
# =============================================================================
PROTOCOL = config.protocol
if   HIGHEST_PROTOCOL < PROTOCOL : PROTOCOL = HIGHEST_PROTOCOL
elif 0                > PROTOCOL : PROTOCOL = DEFAULT_PROTOCOL
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    logger.info ( 'Pickling protocol: %s' % PROTOCOL )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
