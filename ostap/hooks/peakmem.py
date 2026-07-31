#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================
## @file ostap/hooks/peakmem.py
#  Hook: Tracks and reports maximum Resident Set Size (Peak RSS) upon script exit.
#  @code
#  python -m ostap.hooks.peakmem my_script.py --events 100
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-31
# ==========================================================================
""" Hook: Tracks and reports maximum Resident Set Size (Peak RSS) upon script exit.

> python -m ostap.hooks.peakmem my_script.py --events 100
"""
# ==========================================================================
import atexit
import os
import resource
import runpy
import sys
# =========================================================================
## Global flag to prevent duplicate output when registered via atexit
_memcheck_active = False


# ==========================================================================
## Return peak Resident Set Size (RSS) memory used by process in Megabytes.
def get_peak_memory_mb():
    """ Return peak Resident Set Size (RSS) memory used by process in Megabytes.
    """
    usage = resource.getrusage(resource.RUSAGE_SELF)
    # Linux reports maxrss in Kilobytes; macOS (darwin) in Bytes
    if sys.platform == "darwin":
        return usage.ru_maxrss / (1024.0 * 1024.0)
    return usage.ru_maxrss / 1024.0

# ==========================================================================
## Register exit report handler for peak RSS tracking.
def apply_hook():
    """Register exit report handler for peak RSS tracking.
    """
    global _memcheck_active
    if not _memcheck_active:
        _memcheck_active = True
        atexit.register(print_results)

# ==========================================================================
## Print peak process RSS memory usage.
def print_results():
    """Print peak process RSS memory usage.
    """
    global _memcheck_active
    if _memcheck_active:
        peak_os_mb = get_peak_memory_mb()
        print ( f"\nostap.hooks.memcheck Peak RAM (RSS): {peak_os_mb:.2f} MB" )
        _memcheck_active = False

# ==========================================================================
if __name__ == "__main__": # ===============================================

    
    if len ( sys.argv ) <= 1 :
        sys.stderr.write ( "Usage: python -m ostap.hooks.peakmem <script.py> [args...]\n" )
        sys.exit ( 1 )

    # Remove launcher module name from sys.argv
    sys.argv.pop(0)

    target_script = sys.argv [ 0 ]
    if not os.path.isfile(target_script):
        sys.stderr.write ( f"[ostap.hooks.peakmem ERROR: File '{target_script}' not found.\n" )
        sys.exit ( 1 )

    # Register hook handler
    apply_hook()

    # ==========================================================================
    try: # =====================================================================
        # ======================================================================
        runpy.run_path ( target_script , run_name = "__main__" )
        # ======================================================================
    finally: # =================================================================
        # ======================================================================
        # Ensure memory report is printed even if the script crashes
        print_results()

# ==============================================================================
##                                                                      The END
# ==============================================================================
