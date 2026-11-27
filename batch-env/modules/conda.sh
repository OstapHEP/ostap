#!/usr/bin/env bash
# ==============================================================================
# modules/conda.sh
# Isolated Miniconda provisioning inside BATCH_TMP_DIR.
#
# Exports:
#   BATCH_CONDA_PREFIX - Installation prefix of the temporary Miniconda instance.
#   BATCH_CONDA_ENV    - Name of the activated Conda environment.
#   CONDARC            - Path to isolated .condarc inside BATCH_TMP_DIR.
#   CONDA_PKGS_DIRS    - Path to temporary Conda package cache directory.
#   CONDA_ALWAYS_YES   - Set to "true" to suppress interactive prompts.
# ==============================================================================
set -euo pipefail

# ------------------------------------------------------------------------------
# conda_arch_url
#
# Internal helper that detects the host OS and CPU architecture and prints the
# matching official Miniconda3 installer download URL to stdout.
# ------------------------------------------------------------------------------
conda_arch_url() {
  local os arch
  os=$(uname -s)
  arch=$(uname -m)

  case "${os}_${arch}" in
    Linux_x86_64)  printf '%s\n' "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" ;;
    Linux_aarch64) printf '%s\n' "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh" ;;
    Darwin_arm64)  printf '%s\n' "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh" ;;
    Darwin_x86_64) printf '%s\n' "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh" ;;
    *) echo "Unsupported platform: ${os}_${arch}" >&2; return 1 ;;
  esac
}

# ------------------------------------------------------------------------------
# conda_yml_name <yml_file>
#
# Internal helper that extracts a top-level `name:` value from a Conda YAML
# file if present, stripping surrounding quotes and whitespace.
# ------------------------------------------------------------------------------
conda_yml_name() {
  local yml_file="$1"
  awk '
    /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
    /^name:[[:space:]]*/ {
      sub(/^name:[[:space:]]*/, "", $0)
      sub(/[[:space:]]*#.*$/, "", $0)
      gsub(/^[[:space:]"'\'']+|[[:space:]"'\'']+$/, "", $0)
      if (length($0) > 0) {
        print $0
        exit
      }
    }
  ' "$yml_file"
}

# ------------------------------------------------------------------------------
# conda_bootstrap_base
#
# Internal helper that clears inherited Conda variables, installs Miniconda into
# ${BATCH_TMP_DIR}/conda with isolated .condarc and package cache, and loads
# conda shell functions in the current session.
# ------------------------------------------------------------------------------
conda_bootstrap_base() {
  local prefix pkgs_dir condarc url installer

  # Unset inherited Conda session variables so a pre-existing Conda setup on
  # the worker node does not interfere with our temporary installation.
  unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_PROMPT_MODIFIER CONDA_SHLVL CONDA_EXE _CE_CONDA _CE_M || true

  prefix="${BATCH_TMP_DIR:?batch_tmp_init has not been called}/conda"
  pkgs_dir=$(batch_make_dir "conda_pkgs")
  condarc="${BATCH_TMP_DIR}/condarc"
  url=$(conda_arch_url)
  installer="${BATCH_TMP_DIR}/miniconda.sh"

  # Keep configuration and package cache inside BATCH_TMP_DIR so nothing
  # pollutes the user's home directory or persists across batch jobs.
  export CONDARC="$condarc"
  export CONDA_PKGS_DIRS="$pkgs_dir"
  export CONDA_ALWAYS_YES=true

  # Do not pre-create `prefix`: Miniconda's installer fails if the target
  # directory already exists.
  curl -fsSL "$url" -o "$installer"
  bash "$installer" -b -p "$prefix" >/dev/null
  rm -f "$installer"


# Initialize Conda shell functions in the current non-interactive bash session.
  # shellcheck disable=SC1091
  source "$prefix/etc/profile.d/conda.sh"

  # Configure isolated .condarc to use conda-forge cleanly without defaults.
  set +u 
  conda config --file "$condarc" --set auto_activate false
  conda config --file "$condarc" --remove-key channels >/dev/null 2>&1 || true
  conda config --file "$condarc" --add channels conda-forge
  conda config --file "$condarc" --add channels conda-pypi
  conda config --file "$condarc" --set channel_priority strict
  set -u 

  # Accept Anaconda Terms of Service in the isolated prefix if conda tos is present.
  if conda tos --help >/dev/null 2>&1; then
    set +u 
    conda tos accept --override-channels --channel "https://repo.anaconda.com/pkgs/main" >/dev/null 2>&1 || true
    conda tos accept --override-channels --channel "https://repo.anaconda.com/pkgs/r" >/dev/null 2>&1 || true
    set -u 
  fi

  export BATCH_CONDA_PREFIX="$prefix"
}

# ============================================================================
# conda_install_tmp [env_name] [python_version] [extra_conda_pkgs...]
#
# Installs Miniconda into BATCH_TMP_DIR and creates an environment from an
# explicit Python version and optional package list.
#
# Arguments:
#   $1 - Name of the Conda environment to create (default: "workenv").
#   $@ - Optional Conda packages to install (e.g., `libdb gdbm`).
# ============================================================================
conda_install_tmp() {
  local env_name="${1:-batchenv}"
  shift $(( $# >= 1 ? 1 : 0 )) || true

  conda_bootstrap_base

  set +u
  conda create -y -n "$env_name" "$@" >/dev/null
  set -u 

  set +u
  conda activate "$env_name"
  set -u

  export BATCH_CONDA_ENV="$env_name"
}

# ============================================================================
# conda_install_yml <yml_file> [env_name] [extra_conda_pkgs...]
#
# Installs Miniconda into BATCH_TMP_DIR, creates an environment from a Conda
# YAML specification file, activates it, and installs any extra Conda packages
# passed after the environment name. If `env_name` ($2) is omitted or empty,
# uses the `name:` field from `yml_file` (falling back to "workenv").
#
# Arguments:
#   $1 - Path to environment.yml file.
#   $2 - Optional target environment name (default: `name:` from YAML or "workenv").
#   $@ - Optional extra Conda packages to install on top of the YAML env.
# ============================================================================
conda_install_yml() {
  local yml_file="${1:?conda_install_yml: path to environment.yml is required}"
  local env_name="${2:-}"
  shift $(( $# >= 2 ? 2 : $# )) || true

  if [[ ! -f "$yml_file" ]]; then
    echo "conda_install_yml: file not found: $yml_file" >&2
    return 1
  fi

  if [[ -z "$env_name" ]]; then
    env_name=$(conda_yml_name "$yml_file")
    env_name="${env_name:-workenv}"
  fi

  conda_bootstrap_base
 
  set +u
  conda env create -y -f "$yml_file" -n "$env_name" >/dev/null
  set -u 

  set +u
  conda activate "$env_name"
  set -u

  if [[ "$#" -gt 0 ]]; then
    set +u
    conda install -y "$@" >/dev/null
    set -u 
  fi

  export BATCH_CONDA_ENV="$env_name"
}

# ============================================================================
##                                                                     The END
# ============================================================================