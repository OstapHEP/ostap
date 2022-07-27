#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/build_dir.py
#  Create build directoru for ROOT
#  - check environment variable   OSTAP_BUILD_DIR
#  - check the configurtaion file, section [General], field BUILD_DIR
#  - use the temporary directory
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Set  build dir for ROOT
- check environment varibale   OSTAP_BUILD_DIR
- check the configuration file, section [General], field BUILD_DIR
- use the temporary directory
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'build_dir'      , ## Build directory for ROOT
    'make_build_dir' , ## make (temporary) build directory for ROOT
    'UseBuildDir'    , ## use  (temporary) build directory for ROOT
    'use_build_dir'  , ## use  (temporary) build directory for ROOT    
    )
# =============================================================================
from   ostap.utils.basic import make_dir, writeable 
import ROOT, os, glob  
# =============================================================================
build_dir  = None
prefix_dir = 'ostap-build-dir-'
# =============================================================================
# 1) use environment variable 
if not build_dir :

    build_dir = os.environ.get ( 'OSTAP_BUILD_DIR' , '' )
    
    if build_dir and not os.path.exists ( build_dir ) :
        build_dir = make_dir ( build_dir ) 
    if build_dir and not writeable ( build_dir ) :
        logger.warning ('Directory %s is not writeable!' % tmp_dir )
        build_dir = None

## 2) get from configuration file 
if not build_dir :
    ## 2) check the configuration file 
    import ostap.core.config as OCC 
    build_dir = OCC.general.get ( 'BUILD_DIR' , None )
    del OCC
    
    if build_dir and not os.path.exists ( build_dir ) :
        build_dir = make_dir ( build_dir ) 
    if build_dir and not writeable ( build_dir ) :
        logger.warning ('Directory %s is not writeable!' % tmp_dir )
        build_dir = None

## 3) Construct something temporary and unique
if not build_dir :
    
    from ostap.utils.cleanup import CleanUp
    build_dir = CleanUp.tempdir ( prefix = 'ostap-build-' ) 
 
    
# ==============================================================================
## Context manager to use certain build build directory 
#  (useful for multiprocessing environment)
#  @code
#  with UseBuildDir() :
#  ... 
#  @endcode
class UseBuildDir ( object ) :
    """Context manager to use certain build build directory 
    (useful for multiprocessing environment)
    >>> with UseBuildDir() :
    >>> ... 
    """
    def __init__ ( self , build = None ) :
        
        ##
        if build and writeable ( build ) :
            self.__build   = build
            self.__created = False 
        else : 
            self.__created = True             
            self.__build   = make_build_dir ()
            
        self.__prev   = ROOT.gSystem.SetBuildDir ( self.__build )
        
    ## context manager: ENTER 
    def __enter__ ( self ) :

        self.__prev   = ROOT.gSystem.GetBuildDir ()
        ROOT.gSystem.SetBuildDir ( self.__build )        
        return ROOT.gSystem.GetBuildDir ()
    
    ## context manager: EXIT  
    def __exit__ ( self , *_ ) :
        
        ROOT.gSystem.SetBuildDir ( self.__prev )
        if self.__created :  
            from ostap.utils.cleanup import CleanUp
            CleanUp.remove_dir ( self.__build ) 
            
# ==============================================================================
## Context manager to use certain build build directory 
#  (useful for multiprocessing environment)
#  @code
#  with use_build_dir() :
#  ... 
#  @endcode
def use_build_dir ( build = None ) :
    """Context manager to use certain build build directory 
    (useful for multiprocessing environment)
    >>> with Use_build_dir() :
    >>> ... 
    """
    return UseBuildDir ( build )


# ==============================================================================
## create proper temporary directry for ROOT builds 
def make_build_dir ( build = None ) :
    """Create proper temporary directory for ROOT builds
    """

    if not build or not writeable  ( build ) : 
        from ostap.utils.cleanup import CleanUp 
        build = CleanUp.tempdir ( prefix = 'ostap-build-' )

    if not os.path.exists ( build ) :
        make_dir ( build )
        
    return build


# ==============================================================================
# 3) use the temporary directory 
if not build_dir :    
    build_dir = make_build_dir()

# =============================================================================
if build_dir and writeable ( build_dir ) : 
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
