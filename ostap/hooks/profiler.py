#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================
## @file ostap/hooks/profile.py
#  Hook: Runs the script under cProfile to measure CPU performance and function execution counts.
#  Prints top CPU-bound callers upon script exit and optionally dumps stats to a file.
#  @code
#  OSTAP_PROFILE_OUTPUT="result.prof" python -m ostap.hooks.profile my_script.py --events 100
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-31
# ==========================================================================
"""Hook: Runs the script under cProfile to measure CPU performance and function execution counts.
Prints top CPU-bound callers upon script exit and optionally dumps stats to a file.

> OSTAP_PROFILE_OUTPUT="result.prof" python -m ostap.hooks.profile my_script.py --events 100
"""
# ==========================================================================
import atexit
import cProfile
import os
import pstats
import runpy
import sys
# ==========================================================================
# Global profiler instance shared between CLI launch and chain hooks
_profiler = None

# ==========================================================================
## Enable cProfile globally for the current Python process.
def apply_hook():
    """ Enable cProfile globally for the current Python process.
    """
    global _profiler
    if _profiler is None:
        _profiler = cProfile.Profile()
        _profiler.enable()

        # Register exit handler for chain-mode execution
        atexit.register(print_results)

# ==========================================================================
## Disable profiler, print top execution statistics, and optionally save dump file."""
def print_results():
    """ Disable profiler, print top execution statistics, and optionally save dump file.
    """
    global _profiler
    if _profiler is not None:
        _profiler.disable()
        print("\n" + "=" * 30 + " PROFILING RESULTS " + "=" * 30)

        # Print top 50 functions sorted by cumulative time
        stats = pstats.Stats(_profiler).sort_stats("cumulative")
        stats.print_stats(50)

        # ====================================================================
        # Save dump to file if OSTAP_PROFILE_OUTPUT is defined
        # ====================================================================
        dump_file = os.environ.get ( "OSTAP_PROFILE_OUTPUT", "ostap_profile.results" )
        # ====================================================================
        if dump_file : # =====================================================
            # ================================================================
            try : # ==========================================================
                # ============================================================
                from pathvalidate import validate_filepath, ValidationError
                def valid_filename ( f ) :
                    # ========================================================
                    try : # ==================================================
                        # ====================================================
                        validate_filepath ( f , platform = sys.platform )
                        # ====================================================
                    except ValidationError : # ===============================
                        # ====================================================
                        return False
                    return True 
                # ============================================================
            except ImportError : # ===========================================
                # ============================================================
                valid_filepath = lambda s : True 

            # ================================================================
            if valid_filepath ( dump_file ) :
                dump_file = os.path.expandvars ( dump_file )
                dump_file = os.path.expanduser ( dump_file )
                dump_file = os.path.abspath    ( dump_file )
                dump_file = os.path.normpath   ( dump_file )
                dump_dir  = os.path.dirname    ( dump_file ) 
                if os.path.isdir ( dump_dir ) and os.access ( dump_dir , os.W_OK ) :
                    stats.dump_stats ( dump_file )
                    if os.path.isfile ( dump_file ) :
                        print ( f"\nostap.hooks.profile Saved profiler dump to: {dump_file}" )
                        
        print("=" * 79)
        _profiler = None

# ==========================================================================
if __name__ == "__main__": # ===============================================
    # ======================================================================
    
    if len ( sys.argv ) <= 1 :
        sys.stderr.write ( "Usage: OSTAP_PROFILE_OUTPUT='result.prof' python -m ostap.hooks.profile <script.py> [args...]\n" )
        sys.exit ()

    # Remove hook module name from CLI arguments
    sys.argv.pop ( 0 )

    target_script = sys.argv [ 0 ]
    if not os.path.isfile ( target_script ) :
        sys.stderr.write ( f"ostap.hooks.profile ERROR: File '{target_script}' not found.\n" )
        sys.exit ( 1 )

    # ==========================================================================
    # Start profiling using shared hook logic
    # ==========================================================================
    apply_hook()
    
    # ==========================================================================
    try: # =====================================================================
        # ======================================================================
        runpy.run_path ( target_script , run_name = "__main__" )
        # ======================================================================
    finally: # =================================================================
        # ======================================================================
        # Guarantee results output and file saving even if script fails/raises exception
        print_results()

# ==============================================================================
##                                                                      The END
# ==============================================================================
