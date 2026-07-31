#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================
## @file ostap/hooks/batch.py
#  Hook: Forces ROOT into headless (batch) mode. Useful for running jobs
#  on remote clusters, CI/CD pipelines, or servers without a display server (X11/Wayland).
#  @code
#  python -m ostap.hooks.batch my_script.py --events 100
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-31
# ==========================================================================
""" Hook: Forces ROOT into headless (batch) mode. Useful for running jobs
on remote clusters, CI/CD pipelines, or servers without a display server (X11/Wayland).
"""
# ==========================================================================
import os, sys, runpy 
# ==========================================================================
## Configure environment variables and ROOT for headless execution.
def apply_hook():
    """ Configure environment variables and ROOT for headless execution.
    """
    # ======================================================================
    # Disable ROOT graphics display at ROOT C++ level
    try: # =================================================================
        # ==================================================================
        import ROOT    
        ROOT.gROOT.SetBatch(True)
        # ==================================================================
    except ImportError: # ==================================================
        # ==================================================================
        pass
    
# ==========================================================================
if __name__ == "__main__":
    # ======================================================================

    # Remove hook module name from CLI arguments
    sys.argv.pop ( 0 )

    apply_hook()
    
    # ======================================================================
    if len ( sys.argv ) <= 1 :
        sys.stderr.write ( "Usage: python -m ostap.hooks.batch <script.py> [args...]\n" )
        sys.exit ()

    target_script = sys.argv [ 0 ]
    if not os.path.isfile ( target_script ):
        sys.stderr.write ( f"ostap.hooks.batch ERROR: File '{target_script}' not found.\n" )
        sys.exit ( 1 )

    runpy.run_path(target_script, run_name="__main__")

# ==============================================================================
##                                                                      The END
# ==============================================================================
