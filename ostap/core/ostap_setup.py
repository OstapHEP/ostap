#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/ostap_setup.py
#  Basic setup for ROOT/Ostap
#  - Batch
#  - Implicit MR
#  - Build/Cache/TMP  directories
#  - Profiling 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
""" Basic setup for ROOT/Ostap 
- Batch
- Implicit MR
- Build/Cache/TMP  directories 
- Profiling 
"""
# ============================================================================= 
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
)
# ============================================================================= 
import ostap.core.config      as     config 
import ostap.fixes.fixes
import ostap.core.build_dir
import ostap.core.cache_dir
import ostap.utils.cleanup 
import ROOT 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.core.ostap_setup' )
else                       : logger = getLogger( __name__                )
# =============================================================================

# =============================================================================
## Implicit MT ?
# =============================================================================
implicitMT = config.general.getboolean ( 'ImplicitMT' , fallback = False )
if       implicitMT and not ROOT.ROOT.IsImplicitMTEnabled() :
    ROOT.ROOT.EnableImplicitMT  ()
    if     ROOT.ROOT.IsImplicitMTEnabled() : logger.debug ( 'ImplicitMT is enabled'  )
elif not implicitMT and     ROOT.ROOT.IsImplicitMTEnabled() :
    ROOT.ROOT.DisableImplicitMT  ()
    if not ROOT.ROOT.IsImplicitMTEnabled() : logger.debug ( 'ImplicitMT is disabled' )
    
# =============================================================================
## Batch processing ?
# =============================================================================
groot = ROOT.ROOT.GetROOT() 
if   groot and     config.batch and not groot.IsBatch () :
    groot.SetBatch ( True  )
    if     groot.IsBatch() : logger.attention ( "BATCH processig is activated!"   )
elif groot and not config.batch and     groot.IsBatch () :
    groot.SetBatch ( False )
    if not groot.IsBatch() : logger.info      ( "BATCH processig is deactivated!" )

# =============================================================================
## Profile ?
# =============================================================================
if config.general.getboolean ( 'Profile' , fallback = config.profile ) :

    ## dedicated logger 
    plogger  = getLogger ( 'ostap.core.profiler' )

    ## profiler itself
    import cProfile as profile
    profiler = profile.Profile()
    profiler.enable()
    plogger .attention ( 'Profiling is activated!' )
    
    # =========================================================================
    ## "AtExit" action:
    #  - end profiling
    #  - print statistics     
    def _profile_atexit_ ( prof , logger ) :
        """ `At-Exit` action: 
        - end profiling 
        - print statistics 
        """
       
        import io, pstats 
        
        prof.disable()
        logger.info  ( 'End of profiling, generate profiling report' )
        
        sio    = io.StringIO ()
        sortby = pstats.SortKey.CUMULATIVE 
        stat   = pstats.Stats ( prof , stream= sio ).sort_stats ( sortby )
        stat.print_stats()

        report = sio.getvalue()
        logger.info ( 'Profiler report:\n\n%s' % report  )
        logger.info ( 'The end of profiler report' ) 
        
    import atexit 
    atexit.register ( _profile_atexit_ , profiler , plogger ) 
    

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme   import docme
    docme ( __name__ , logger = logger )
    

# =============================================================================

# ============================================================================= 
    
