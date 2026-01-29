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
    if     ROOT.ROOT.IsImplicitMTEnabled() : logger.info      ( 'ImplicitMT is enabled'  )
elif not implicitMT and     ROOT.ROOT.IsImplicitMTEnabled() :
    ROOT.ROOT.DisableImplicitMT  ()
    if not ROOT.ROOT.IsImplicitMTEnabled() : logger.attention ( 'ImplicitMT is disabled' )
    
# =============================================================================
## Batch processing ?
# =============================================================================
groot = ROOT.ROOT.GetROOT() 
if   groot and     config.batch and not groot.IsBatch () :
    groot.SetBatch ( True  )
    if     groot.IsBatch() : logger.attention ( "BATCH propcessig is activated!"   )
elif groot and not config.batch and     groot.IsBatch () :
    groot.SetBatch ( False )
    if not groot.IsBatch() : logger.info      ( "BATCH propcessig is deactivated!" )

# =============================================================================
## Profile ?
if config.general.getboolean ( 'Profile' , fallback = False ) :
    
    import cProfile as _profile
    _pr = _profile.Profile()
    
    # ========================================================================
    def _profile_atexit_ ( prof ) : 

        import io, pstats 
        
        logger.debug ( 'End of profiling, generate profiling report' )
        prof.disable()
        
        _sio   = io.StringIO()
        sortby = 'cumulative'
        _pstat = pstats.Stats( prof , stream=_sio).sort_stats('cumulative')
        _pstat.print_stats()
        print(_sio.getvalue())
            
    import atexit 
    atexit.register ( _profile_at_exit_ , _pr ) 
    _pr.enable()
    del _pr 
    logger.attention ( 'Profiling is activated!' )


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme   import docme
    docme ( __name__ , logger = logger )
    

# =============================================================================

# ============================================================================= 
    
