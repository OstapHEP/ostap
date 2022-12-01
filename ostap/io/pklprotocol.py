#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file pklprotocol.py
# Helper module to pickup proper pickling protocol 
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2010-04-30
# =============================================================================
""" Helper module to pickup proper pickling protocol 
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
    )
# =============================================================================
import os
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.pklprotocol' )
else                      : logger = getLogger ( __name__               )
# =============================================================================
try:
    from cPickle   import HIGHEST_PROTOCOL, DEFAULT_PROTOCOL
except ImportError:
    from  pickle   import HIGHEST_PROTOCOL, DEFAULT_PROTOCOL
# =============================================================================
## helper function to get the protocol 
def get_protocol ( p ) :
    """helper function to get the protocol"""    
    if   p.lower () in ( 'default' , 'def'  ) :
        return DEFAULT_PROTOCOL            
    elif p.lower () in ( 'highest' , 'high' ) :
        return HIGHEST_PROTOCOL
    elif p.lower () in ( 'compat'  , 'compatible' , 'backward_compatible' ) :
        return 2 
    else :
        return int ( p )
# =============================================================================    
## pickle protocol to be used 
PROTOCOL = None
#  (1) pickup protocol form environment variable 
try :
    pe = os.environ.get ( 'OSTAP_PROTOCOL'  , '' )
    pp = get_protocol   ( pe )
    if  0 <= pp <= HIGHEST_PROTOCOL :
        PROTOCOL = pp
        logger.debug ( "Protocol %s is picked from 'OSTAP_PROTOCOL=%s' environment" % ( PROTOCOL , pe ) )
except :
    pass
# =============================================================================
#  (2) take protocol from the configuration files 
if PROTOCOL is None :
    try : 
        import ostap.core.config as OCC
        pe = OCC.general.get ( 'Protocol' , '' )
        pp = get_protocol ( pe )
        if  0 <= pp <= HIGHEST_PROTOCOL :
            PROTOCOL = pp 
            logger.debug ( "Protocol %s is picked from 'General: protocol=%s' section" %  ( PROTOCOL , pe ) ) 
    except :
        pass 
# =============================================================================
#  (3) use default protocol 
if PROTOCOL is None :
    PROTOCOL = DEFAULT_PROTOCOL 
    logger.debug ( "Default protocol %s is used" % PROTOCOL  )
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
