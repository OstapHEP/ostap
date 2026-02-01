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
import ostap.core.config as     config
import ROOT, os, glob  
# =============================================================================
build_dir  = config.build_dir
prefix_dir = 'ostap-build-'
# =============================================================================
## create proper build directry for ROOT builds
#  - existing directory with correct permissions ? use it
#  - otherwise create the temporary 
def make_build_dir ( build = None ) :
    """ Create proper temporary directory for ROOT builds
    - existing directory with correct permissions ? use it
    - otherwise create the temporary 
    """
    
    if build and \
       os.path.exists ( build ) and \
       os.path.isdir  ( build ) and \
       os.access      ( build , os.W_OK ) : return build
        
    from ostap.utils.cleanup import CleanUp
    return CleanUp.tempdir ( prefix = 'ostap-build-' )

# =============================================================================
## directory is specified? 
if build_dir : # ==============================================================
    # =========================================================================
    build_dir = os.path.expandvars ( build_dir )
    build_dir = os.path.expanduser ( build_dir )
    build_dir = os.path.expandvars ( build_dir )
    build_dir = os.path.expanduser ( build_dir )
    # =========================================================================
    if not os.path.exists ( build_dir ): # ====================================
        # =====================================================================
        try : # ===============================================================
            # =================================================================
            os.makedirs ( build_dir , exist_ok = True ) 
            # =================================================================
        except OSError : # ====================================================
            # =================================================================
            pass

# =============================================================================
##  final check
#  - good directoy? use it!
#  - create a temporary oherwise 
build_dir = make_build_dir ( build_dir )

# =============================================================================
## Context manager to use certain build build directory 
#  (useful for multiprocessing environment)
#  @code
#  with UseBuildDir() :
#  ... 
#  @endcode
class UseBuildDir ( object ) :
    """ Context manager to use certain build directory 
    (useful for multiprocessing environment)
    >>> with UseBuildDir() :
    >>> ... 
    """
    def __init__ ( self , build = None ) :
        
        ##
        if build and \
           os.path.exists ( build ) and \
           os.path.isdir  ( build ) and \
           os.access      ( build , os.W_OK ) :
            
            self.__build   = build
            self.__created = False
            
        else :
            
            self.__build   = make_build_dir ()
            self.__created = True             
            
        self.__prev   = ROOT.gSystem.SetBuildDir ( self.__build )
        
    # =========================================================================
    ## context manager: ENTER 
    def __enter__ ( self ) :

        self.__prev   = ROOT.gSystem.GetBuildDir ()
        ROOT.gSystem.SetBuildDir ( self.__build )        
        return ROOT.gSystem.GetBuildDir ()

    # =========================================================================
    ## context manager: EXIT  
    def __exit__ ( self , *_ ) :
        # =====================================================================
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
    """ Context manager to use certain build build directory 
    (useful for multiprocessing environment)
    >>> with Use_build_dir() :
    >>> ... 
    """
    return UseBuildDir ( build )

# =============================================================================
## rdefien the default build directory: 
if build_dir and \
   os.path.exists ( build_dir ) and \
   os.path.isdir  ( build_dir ) and \
   os.access      ( build_dir , os.W_OK ) : 
    ROOT.gSystem.SetBuildDir ( build_dir )
    
# =============================================================================
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
