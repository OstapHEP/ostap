#!/usr/bin/env bash
# ==============================================================================
# batch.sh
# Main entry point for the batch-env library.
#
# Intended usage:
#   Source this file at the top of a batch job script (e.g. lxbatch/HTCondor,
#   Ganga/Grid) to initialize an isolated temporary workspace with automatic
#   cleanup and load helper modules for Git and Conda/pip.
#
# Example:
#   source "/path/to/batch-env/batch.sh"
#   repo_dir=$(step "Clone" repo_clone "https://github.com/org/repo.git" "main")
#   step "Conda" conda_install_tmp "myenv" "3.11"
#   step "Pip" conda_pip_install -r "${repo_dir}/requirements.txt"
#   step "Run" python "${repo_dir}/worker.py" --out "${BATCH_TMP_DIR}/out.root"
# ==============================================================================
set -euo pipefail

# Resolve the absolute directory path where batch.sh resides so modules can be
# sourced reliably regardless of the caller's current working directory.
BATCH_ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load core modules (SC1091 disabled because paths are resolved dynamically).

# shellcheck disable=SC1091
source "${BATCH_ENV_DIR}/modules/tmpdir.sh"
# shellcheck disable=SC1091
source "${BATCH_ENV_DIR}/modules/git.sh"
# shellcheck disable=SC1091
source "${BATCH_ENV_DIR}/modules/conda.sh"
# shellcheck disable=SC1091
source "${BATCH_ENV_DIR}/modules/pip.sh"
# shellcheck disable=SC1091
source "${BATCH_ENV_DIR}/modules/berkeleydb.sh"

# ============================================================================
# step <title> <command> [args...]
#
# Executes a pipeline step with a formatted progress header printed to stderr.
# Printing to stderr keeps stdout clean so command substitutions like
# `dir=$(step "Clone" repo_clone ...)` capture only the function's return value.
#
# Arguments:
#   $1 - Human-readable step label.
#   $@ - Command and arguments to execute.
# ============================================================================
step() {
  local title="${1:?step: step title is required}"
  shift
  printf '[step] %s\n' "$title" >&2
  "$@"
}

# Automatically initialize the temporary workspace and register cleanup traps.
batch_tmp_init

# ============================================================================
##                                                                     The END
# ============================================================================  