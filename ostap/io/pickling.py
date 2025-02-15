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
from   pickle         import ( Pickler          ,
                               Unpickler        , 
                               DEFAULT_PROTOCOL ,
                               HIGHEST_PROTOCOL ,
                               dumps , loads    )  
from   io              import BytesIO 
import sys, array 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.pickling' )
else                      : logger = getLogger ( __name__               )
# =============================================================================
## helper function to get the protocol 
def get_protocol ( p ) :
    """ helper function to get the protocol"""    
    if   p.lower () in ( 'default' , 'def'  ) : return DEFAULT_PROTOCOL            
    elif p.lower () in ( 'highest' , 'high' ) : return HIGHEST_PROTOCOL
    elif p.lower () in ( 'compat'  , 'compatible' , 'backward_compatible' ) :
        return min ( 2 , HIGHEST_PROTOCOL )  
    else :
        return int ( p )
# =============================================================================    
## pickle protocol to be used 
PROTOCOL = None
#  (1) pickup protocol from the environment variable
# =============================================================================
try : # =======================================================================
    # =========================================================================
    from ostap.utils.basic import get_env, OSTAP_PROTOCOL  
    pe = get_env ( OSTAP_PROTOCOL  , '' )
    pp = get_protocol   ( pe )
    if pp < 0 : pp = HIGHEST_PROTOCOL 
    if 0  <= pp <= HIGHEST_PROTOCOL :
        PROTOCOL = pp
        logger.debug ( "Protocol %s is picked from 'OSTAP_PROTOCOL=%s' environment" % ( PROTOCOL , pe ) )
    # =========================================================================
except : # ====================================================================
    # =========================================================================
    pass
# =============================================================================
#  (2) take protocol from the configuration files
# =============================================================================
if PROTOCOL is None : # =======================================================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import ostap.core.config as OCC
        pe = OCC.general.get ( 'Protocol' , '' )
        pp = get_protocol ( pe )
        if pp < 0 : pp = HIGHEST_PROTOCOL 
        if 0  <= pp <= HIGHEST_PROTOCOL :
            PROTOCOL = pp 
            logger.debug ( "Protocol %s is picked from 'General: protocol=%s' section" %  ( PROTOCOL , pe ) )
        # =====================================================================
    except : # ================================================================
        # =====================================================================
        pass 
# =============================================================================
#  (3) use default protocol
# =============================================================================
if PROTOCOL is None : # =======================================================
    # =========================================================================
    PROTOCOL = DEFAULT_PROTOCOL 
    logger.debug ( "Default protocol %s is used" % PROTOCOL  )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    logger.info ( 'Pickling protocol: %s' % PROTOCOL )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
