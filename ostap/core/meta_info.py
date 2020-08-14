#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/meta_info.py
#  Meta-info for ostap 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-05-30
# =============================================================================
"""Core objects for ostap 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'meta_info'        , ## some meta-information for ostap
    )
# =============================================================================
from collections import namedtuple 
MetaInfo = namedtuple ( 'MetaInfo' , ( 'User' , 'Ostap' , 'Python' , 'ROOT'  ) ) 

_meta = [] 
# =============================================================================
## get some meta-info for ostap
#  @code
#  info = meta_info ()  
#  @endcode
def meta_info ( ) :
    """Get some meta-info for ostap
    >>> info = meta_info ()  
    """
    if _meta : return _meta[0]
    
    import sys
    import getpass            
    import socket
    import ROOT
    from   ostap import version as ostap_version
    user   = "%s@%s" % ( getpass.getuser() , socket.getfqdn () )
    python_version = '%d.%d.%d' % ( sys.version_info.major ,
                                    sys.version_info.minor ,
                                    sys.version_info.micro )                         
    
    _meta.append ( MetaInfo ( user ,
                              ostap_version ,
                              python_version ,
                              ROOT.gROOT.GetVersion() ) )
    return _meta[0] 

# =============================================================================
if '__main__' == __name__ :
    
    # =========================================================================
    from ostap.logger.logger import getLogger 
    if '__main__' ==  __name__ : logger = getLogger( 'ostap.core.meta_info' )
    else                       : logger = getLogger( __name__     )
    # =========================================================================
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    info = meta_info()
    logger.info ( ' %10s : %s' % ( 'User'   , info.User   ) )
    logger.info ( ' %10s : %s' % ( 'Ostap'  , info.Ostap  ) )  
    logger.info ( ' %10s : %s' % ( 'Python' , info.Python ) )  
    logger.info ( ' %10s : %s' % ( 'ROOT'   , info.ROOT   ) )  
                  
    
# =============================================================================
##                                                                      The END 
# =============================================================================
    
