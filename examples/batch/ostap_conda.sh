#!/usr/bin/env bash
set -euo pipefail

# Prevent sourcing the script so environment variables do not leak into the caller shell
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    echo "ERROR: This script must be executed, not sourced (run './$(basename "${BASH_SOURCE[0]}")' or 'bash ...')." >&2
    return 1 2>/dev/null || exit 1
fi

# Unset inherited conda state so a pre-activated conda session does not interfere
unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_PROMPT_MODIFIER CONDA_SHLVL CONDA_EXE _CE_CONDA _CE_M || true

# =============================================================================
## 1. Create a temporary working directory and ensure it is removed on exit (even on error)
# =============================================================================
WORK_DIR="$(mktemp -d -t "ostap-$(date '+%Y-%b-%d')-OSTAP-CONDA-XXXXXX")"
MINICONDA_DIR="$WORK_DIR/miniconda3"
ENV_DIR="$WORK_DIR/env"
INSTALLER="$WORK_DIR/miniconda.sh"
ENV_ACTIVATED=0

# Isolate conda configuration and package cache inside the temporary directory
export CONDARC="$MINICONDA_DIR/.condarc"
export CONDA_PKGS_DIRS="$WORK_DIR/pkgs"
export CONDA_ALWAYS_YES=true

_remove_dir_safely() {
    local target="${1:-}"
    if [[ -z "${target}" ]] || [[ "${target}" == "/" ]] || [[ "${target}" == "." ]] || [[ "${target}" == ".." ]]; then
        return 0
    fi
    if [[ -d "${target}" ]]; then
        echo ">>> Cleaning up directory: ${target}"
        rm -rf "${target}"
    fi
}

cleanup() {
    local exit_code=$?
    trap - EXIT

    # Deactivate only if this script activated its own conda environment
    if [[ "${ENV_ACTIVATED:-0}" -eq 1 ]] && declare -F conda >/dev/null 2>&1; then
        conda deactivate >/dev/null 2>&1 || true
        ENV_ACTIVATED=0
    fi

    _remove_dir_safely "${WORK_DIR:-}"

    exit "$exit_code"
}

if declare -F trap_add >/dev/null 2>&1; then
    trap_add cleanup EXIT
else
    trap cleanup EXIT
fi

# =============================================================================
## 2. Download and install Miniconda into the temporary directory (without modifying ~/.bashrc)
# =============================================================================
echo ">>> Downloading Miniconda..."
curl -fsSL "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" -o "$INSTALLER"

echo ">>> Installing Miniconda into $MINICONDA_DIR..."
bash "$INSTALLER" -b -p "$MINICONDA_DIR"
rm -f "$INSTALLER"

# Initialize conda for the current script session only
# shellcheck source=/dev/null
source "$MINICONDA_DIR/etc/profile.d/conda.sh"

# =============================================================================
## 3. Configure conda (suppress prompts and set channels in an isolated .condarc)
# =============================================================================
conda config --file "$CONDARC" --set auto_activate_base false
conda config --file "$CONDARC" --add channels conda-forge
conda config --file "$CONDARC" --set channel_priority strict

# =============================================================================
## 4. Create an isolated environment with the desired Python version
# =============================================================================
echo ">>> Creating conda environment at $ENV_DIR..."
conda create --prefix "$ENV_DIR" --yes ostaphep root ipython libdb gdbm 

# Activate the environment within the script
conda activate "$ENV_DIR"
ENV_ACTIVATED=1

# =============================================================================
## 5. Install required packages using python -m pip
# =============================================================================
echo ">>> Installing packages via python -m pip..."
python -m pip install --upgrade pip setuptools wheel              
python -m pip install cpufeature pragmastat terminaltables3 torch

# =============================================================================
## 6. Install Berkeley DB and build Python package `berkeleydb`
# =============================================================================
export BERKELEYDB_DIR="$CONDA_PREFIX"
export YES_I_HAVE_THE_RIGHT_TO_USE_THIS_BERKELEY_DB_VERSION=1
export CPATH="$CONDA_PREFIX/include${CPATH:+:$CPATH}"
export LIBRARY_PATH="$CONDA_PREFIX/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# Quick sanity check that db.h exists inside the conda prefix
test -f "$CONDA_PREFIX/include/db.h"

python -m pip install --no-build-isolation "berkeleydb" --config-settings="--berkeley-db=$CONDA_PREFIX" \
    || python -m pip install --no-build-isolation "berkeleydb"

python -m pip list

# =============================================================================
## 7. print all dependencies   
# =============================================================================
check-dependencies

# =============================================================================
## 8. Run the target Python script
# =============================================================================
echo ">>> Running target script..."
python -c '
import sys, numpy, requests
print("Python:", sys.version)
print("NumPy:", numpy.__version__)
'

# Or run an external script file:
# python /path/to/your_script.py

# =============================================================================
## 9. Deactivate the conda environment
# =============================================================================
echo ">>> Deactivating conda environment..."
conda deactivate >/dev/null 2>&1 || true
ENV_ACTIVATED=0

# On script exit, the trap cleanup handler will automatically delete $WORK_DIR

# =============================================================================
##                                                                      The END
# =============================================================================
