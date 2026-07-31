#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/hooks/chain.py
#  Hook: Chain launcher allowing combination of multiple hooks via environment variable or CLI flag.
#  Example usage:
#  @code
#  OSTAP_HOOKS="spawn,batch,quiet" python -m ostap.hooks.chain my_script.py --events 100
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-31
# =============================================================================
"""Hook: Chain launcher allowing combination of multiple hooks via environment variable or CLI flag.

Example usage:

    OSTAP_HOOKS="spawn,batch,quiet" python -m ostap.hooks.chain my_script.py --events 100
"""
# =============================================================================
import importlib
import os
import runpy
import sys

# =============================================================================
## Dynamically import and run apply_hook() for each specified hook in order.
def apply_chain(hook_names):
    """Dynamically import and run apply_hook() for each specified hook in order."""
    for name in hook_names:
        name = name.strip()
        if not name:
            continue
        # ======================================================================
        try:
            # ==================================================================
            module = importlib.import_module( f"ostap.hooks.{name}" )
            if module and hasattr  ( module, "apply_hook" ) :
                module.apply_hook()
                print ( f"[ostap.hooks.chain] Applied hook: {name}" )
            # ==================================================================
        except ImportError:
            # ==================================================================
            sys.stderr.write ( f"ostap.hooks.chain WARNING: Hook '{name}' not found.\n" )


# ==============================================================================
if __name__ == "__main__":
    # ==========================================================================
    # Parse hook list from environment variable OSTAP_HOOKS
    hooks_env      = os.environ.get ( "OSTAP_HOOKS" , "" )
    hooks_to_apply = [ h.strip () for h in hooks_env.split ( "," ) if h.strip() ]

    if hooks_to_apply: apply_chain ( hooks_to_apply )

    # ==========================================================================
    if len ( sys.argv ) <= 1:
        sys.stderr.write ( "Usage: OSTAP_HOOKS='spawn,batch' python -m ostap.hooks.chain <script.py> [args...]\n" )
        sys.exit ()

    sys.argv.pop ( 0 )
    target_script = sys.argv [ 0 ]
    if not os.path.isfile(target_script):
        sys.stderr.write ( f"ostap.hooks.chain   ERROR: File '{target_script}' not found.\n" )
        sys.exit ( 1 )

    runpy.run_path ( target_script , run_name = "__main__" )

# ==============================================================================
##                                                                      The END
# ==============================================================================
