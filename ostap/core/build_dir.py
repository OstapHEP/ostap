#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/build_dir.py
#  Create build directoru fro ROOT
#  - check environment varibale   OSTAP_BUILD_DIR
#  - check the configurtaion file, section [General], field BUILD_DIR
#  - use the temporary directory
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Set  build dir for ROOT
- check environment varibale   OSTAP_BUILD_DIR
- check the configurtaion file, section [General], field BUILD_DIR
- use the temporary directory
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'build_dir' , ## Build directory for ROOT 
    )
# =============================================================================
import ROOT, os 
from   ostap.utils.basic import good_dir, make_dir 

build_dir = None

# =============================================================================
## 1) check the environment variable 
if 'OSTAP_BUILD_DIR' in os.environ :
    
    bdir = os.environ.get ( 'OSTAP_BUILD_DIR' , '' )
    
    if  good_dir ( bdir )                     : build_dir = bdir
    elif bdir and not os.path.exists ( bdir ) : make_dir ( bdir )
        
    if good_dir ( bdir ) : build_dir = bdir
        
# =============================================================================
## 2) use the configuration file 
if not build_dir :
    
    import ostap.core.config as _CONFIG
    if 'BUILD_DIR' in _CONFIG.general :
        
        bdir = _CONFIG.general.get ( 'BUILD_DIR' , fallback = '' )
        
        if   good_dir ( bdir )                    : build_dir = bdir
        elif bdir and not os.path.exists ( bdir ) : make_dir ( bdir ) 
        
        if good_dir ( bdir ) : build_dir = bdir

# ==============================================================================
# 3) use the temporary directory 
if not build_dir :

    from ostap.utils.cleanup import CleanUp as _CU
    bdir = _CU.tempdir ( prefix = 'build_' )
    if good_dir ( bdir ) : build_dir = bdir
    del _CU

# =============================================================================
if build_dir and good_dir ( build_dir ) : 
    ROOT.gSystem.SetBuildDir ( build_dir ) 

build_dir = ROOT.gSystem.GetBuildDir ()
    
# =============================================================================
if '__main__' == __name__ :

    # =========================================================================
    # logging 
    # =========================================================================
    from ostap.logger.logger import getLogger 
    logger = getLogger( 'ostap.core.build_dir' )

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    logger.info ('Build directory for ROOT is %s'  %  build_dir )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
