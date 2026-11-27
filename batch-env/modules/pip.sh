#!/usr/bin/env bash
# ==============================================================================
# modules/pip.sh
# Python package installation helper using pip in the active environment.
#
# Note:
#   Works with any active Python interpreter — either a temporary Conda
#   environment created via `conda_install_tmp` or a pre-existing system/LCG
#   Python environment (or a virtualenv).
# ==============================================================================
set -euo pipefail

# ==============================================================================
# pip_install [pip_install_args...]
#
# Upgrades pip in the currently active Python environment and installs any
# packages or requirement files passed as arguments (e.g. `-r requirements.txt`).
# Does nothing if called with zero arguments.
#
# Arguments:
#   $@ - Arguments forwarded directly to `python -m pip install`.
# ============================================================================
pip_install() {
  if [[ "$#" -eq 0 ]]; then
    return 0
  fi

  python -m pip install --upgrade pip >/dev/null
  python -m pip install "$@"
}

# ============================================================================
##                                                                     The END
# ============================================================================