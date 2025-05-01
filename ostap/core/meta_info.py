#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/meta_info.py
#  Meta-info for ostap 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-05-30
# =============================================================================
"""Verison metainfo for Ostap
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'meta_info'            , ## some meta-information for ostap
    'user'                 , ## user
    #
    'root_version'         , ## version of ROOT
    'ostap_version'        , ## version of Ostap
    'python_version'       , ## verison of Python
    #
    'root_version_int'     , ## version of ROOT   (as int) 
    'ostap_version_int'    , ## version of ROOT   (as int) 
    'python_version_int'   , ## version of PYTHON (as int)
    #
    'ostap_info'           , ## Ostap version info 
    'root_info'            , ## ROOT version info 
    'python_info'          , ## Python version info
    #
    'old_PyROOT'           , ## do we use "old" PyROOT ?
    )
# =============================================================================
from collections import namedtuple 
from   ostap     import version      as ostap_version
from   ostap     import version_info as ostap_info  
from   ostap     import version_int  as ostap_version_int 
import os, sys, socket, ROOT 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.core.meta_info' )
else                       : logger = getLogger( __name__     )
# =============================================================================
MetaInfo = namedtuple ( 'MetaInfo'  , ( 'User' , 'Ostap' , 'Python' , 'ROOT'  ) ) 
RootInfo = namedtuple ( 'RootInfo'  , ( 'major' , 'minor' , 'patch'           ) ) 

# =============================================================================
## Who am I ?
#  @cdoe
#  print ( 'I am :' % whoami() ) 
#  @endcode 
def whoami () :
    """ Who am I ?
    >>> print ( "I am ", whoami() ) 
    """
    # ==========================================================================
    try : # ====================================================================
        # ======================================================================
        return os.getlogin() # =================================================    
    except : # =================================================================
        # ======================================================================
        pass
    
    import getpass
    return getpass.getuser() 

# ==============================================================================
## Full user name
user           = "%s@%s"    % ( whoami () , socket.getfqdn () )
#
python_version = '%d.%d.%d' % ( sys.version_info.major ,
                                sys.version_info.minor ,
                                sys.version_info.micro )

python_version_int = sys.version_info.micro              + \
                     sys.version_info.minor       * 100  + \
                     sys.version_info.major * 100 * 100

groot = ROOT.ROOT.GetROOT()
root_version     = groot.GetVersion    ()
root_version_int = groot.GetVersionInt ()
old_PyROOT       = groot.GetVersionInt () < 62200

root_major       = divmod ( root_version_int                       , 100**2 ) [0]
root_minor       = divmod ( root_version_int - root_major * 100**2 , 100    ) [0]
root_patch       = root_version_int - root_major * 100**2 - root_minor * 100 

root_info        = RootInfo (  root_major ,  root_minor , root_patch ) 
meta_info        = MetaInfo ( user           ,
                              ostap_version  ,
                              python_version ,
                              root_version   )
python_info     = sys.version_info

# =============================================================================
assert ( 3 ,  8 ) <  python_info, "No support for PYTHON %s " %  ( '.'.join ( str ( i ) for i in python_info  ) )
assert ( 6 , 24 ) <= root_info  , "No support for ROOT %s "   %  ( '.'.join ( str ( i ) for i in   root_info  ) )
# =============================================================================

# =============================================================================
if '__main__' == __name__ :
    
    # =========================================================================
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    logger.info ( ' %12s : %s' % ( 'User'        , meta_info.User     ) )
    logger.info ( ' %12s : %s' % ( 'Ostap'       , meta_info.Ostap    ) )  
    logger.info ( ' %12s : %s' % ( 'Python'      , meta_info.Python   ) )  
    logger.info ( ' %12s : %s' % ( 'ROOT'        , meta_info.ROOT     ) )
    
    logger.info ( ' %12s : %s' % ( 'Ostap/vers'  ,  ostap_version     ) )  
    logger.info ( ' %12s : %s' % ( 'Python/vers' , python_version     ) )  
    logger.info ( ' %12s : %s' % ( 'ROOT/vers'   ,   root_version     ) )  

    logger.info ( ' %12s : %d' % ( 'Ostap/int'   ,  ostap_version_int ) )  
    logger.info ( ' %12s : %d' % ( 'Python/int'  , python_version_int ) )  
    logger.info ( ' %12s : %d' % ( 'ROOT/int'    ,   root_version_int ) )  

    logger.info ( ' %12s : %s' % ( 'Ostap/info'  ,  ostap_info        ) )
    logger.info ( ' %12s : %s' % ( 'ROOT/info'   ,   root_info        ) )
    logger.info ( ' %12s : %s' % ( 'Python/info' , python_info        ) )

# =============================================================================
##                                                                      The END 
# =============================================================================
    
