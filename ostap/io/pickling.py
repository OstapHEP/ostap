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
import os, sys
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.pickling' )
else                      : logger = getLogger ( __name__               )
# =============================================================================
if  (3, 0 ) <= sys.version_info :
    from pickle import ( Pickler, Unpickler, 
                         DEFAULT_PROTOCOL, HIGHEST_PROTOCOL,
                         dumps, loads,
                        PicklingError, UnpicklingError ) 
else : 
    DEFAULT_PROTOCOL = 2 
    try:
        from cPickle   import ( Pickler, Unpickler, HIGHEST_PROTOCOL,
                                dumps, loads,
                                PicklingError, UnpicklingError ) 
    except ImportError:
        from  pickle   import ( Pickler, Unpickler, HIGHEST_PROTOCOL,
                                dumps, loads,
                                PicklingError, UnpicklingError ) 
    DEFAULT_PROTOCOL = min ( DEFAULT_PROTOCOL , HIGHEST_PROTOCOL )
    
# =============================================================================    
try :
    from io            import BytesIO 
except ImportError :
    try:
        from cStringIO import StringIO as BytesIO 
    except ImportError:
        from  StringIO import StringIO as BytesIO
# =============================================================================
## Primitive check if the object casn be pickled and unpickled  
def pickles ( obj ) :
    """Primitive check if the object casn be pickled and unpickled 
    """
    try:
        pkl = loads ( dumps ( obj ) )
        return pkl == obj 
    except ( PicklingError, UnpicklingError , AttributeError ) : 
        return False

# =============================================================================
## Check pickling of an object across another process
def check ( obj ):
    """Check pickling of an object across another process
    """
    import subprocess
    fail = True
    try:
        _obj = dumps ( obj )
    except ( PicklingError , AttributeError ) :
        return None 
    ## 
    msg = "%s -c import pickle; print(pickle.loads(%s))" %  ( python , repr ( _obj ) )
    return subprocess.call ( msg.split ( None , 2 ) )

# =============================================================================
## helper function to get the protocol 
def get_protocol ( p ) :
    """helper function to get the protocol"""    
    if   p.lower () in ( 'default' , 'def'  ) :
        return DEFAULT_PROTOCOL            
    elif p.lower () in ( 'highest' , 'high' ) :
        return HIGHEST_PROTOCOL
    elif p.lower () in ( 'compat'  , 'compatible' , 'backward_compatible' ) :
        return min ( 2 , HIGHEST_PROTOCOL )  
    else :
        return int ( p )
# =============================================================================    
## pickle protocol to be used 
PROTOCOL = None
#  (1) pickup protocol from the environment variable 
try :
    from ostap.utils.basic import get_env as ostap_getenv 
    ##  pe = os.environ.get ( 'OSTAP_PROTOCOL'  , '' )
    pe = ostap_getenv ( 'OSTAP_PROTOCOL'  , '' )
    pp = get_protocol   ( pe )
    if pp < 0 : pp = HIGHEST_PROTOCOL 
    if 0  <= pp <= HIGHEST_PROTOCOL :
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
        if pp < 0 : pp = HIGHEST_PROTOCOL 
        if 0  <= pp <= HIGHEST_PROTOCOL :
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
    logger.info ( 'Pickling protocol: %s' % PROTOCOL )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
