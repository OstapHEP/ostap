#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/hooks/fork.py
#  Hook: Forces the 'forkserver' start method for multiprocessing before any C++/ROOT
#  modules are initialized. Prevents deadlocks and segmentation faults when
#  forking processes with active C++ threads or ROOT contexts.
#  @code
#  python -m ostap.hooks.fork my_script.py --events 100
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-31
# ==============================================================================
""" Hook: Forces the 'forkserver' start method for multiprocessing before any C++/ROOT
modules are initialized. Prevents deadlocks and segmentation faults when
forking processes with active C++ threads or ROOT contexts.

> python -m ostap.hooks.fork my_script.py --events 100
"""
# ==============================================================================
import os, sys, runpy
# ==============================================================================
METHOD = 'fork'
# ==============================================================================
## Set process start method to 'forks' for multiprocess(ing) modules.
def apply_hook():
    """ Set process start method to 'fork' for multiprocess(ing) modules.
    """
    # ==========================================================================
    try: # =====================================================================
        # ======================================================================
        import multiprocessing
        multiprocessing.set_start_method ( METHOD , force = True )
        # ======================================================================
    except ( ImportError , RuntimeError ) : # ==================================
        # ======================================================================
        pass
    
    # ==========================================================================
    try: # =====================================================================
        # ======================================================================
        import multiprocess
        multiprocess.set_start_method ( METHOD , force=True)
        # ======================================================================        
    except ( ImportError , RuntimeError ) : # ==================================
        # ======================================================================
        pass

# ==============================================================================
if __name__ == "__main__": # ===================================================
    # ==========================================================================
    
    if len ( sys.argv ) <= 1 :
        sys.stderr.write ( "Usage: python -m ostap.hooks.fork <script.py> [args...]\n" )
        sys.exit ()

    # Remove launcher module name from sys.argv
    sys.argv.pop ( 0 )

    target_script = sys.argv[0]
    if not os.path.isfile(target_script):
        sys.stderr.write( f"ostap.hooks.spawn ERROR: File '{target_script}' not found.\n" )
        sys.exit ( 1 )

    # Apply hook right before running target script
    apply_hook()

    runpy.run_path ( target_script , run_name = "__main__" )

# ==============================================================================
##                                                                      The END
# ==============================================================================
