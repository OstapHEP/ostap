#!/usr/bin/env bash
# ==============================================================================
# modules/berkeleydb.sh
# Helper for building and installing the Python `berkeleydb` package against
# Berkeley DB headers and libraries (such as those provided via Conda `libdb`).
#
# Exports:
#   BERKELEYDB_DIR
#   YES_I_HAVE_THE_RIGHT_TO_USE_THIS_BERKELEY_DB_VERSION
#   CPATH
#   LIBRARY_PATH
#   LD_LIBRARY_PATH
# ==============================================================================
set -euo pipefail

# ============================================================================
# berkeleydb_install [db_prefix]
#
# Configures compiler/linker paths for Berkeley DB, verifies that `db.h` is
# present, and installs the `berkeleydb` Python package via pip without build
# isolation.
#
# Arguments:
#   $1 - Prefix directory containing `include/db.h` and `lib/`
#        (default: active `${CONDA_PREFIX:-${BATCH_CONDA_PREFIX:-}}`).
# ============================================================================
berkeleydb_install() {
  local db_prefix="${1:-${CONDA_PREFIX:-${BATCH_CONDA_PREFIX:-}}}"

  if [[ -z "$db_prefix" ]]; then
    echo "berkeleydb_install: Berkeley DB prefix is not set (activate Conda or pass prefix)" >&2
    return 1
  fi

  if [[ ! -f "${db_prefix}/include/db.h" ]]; then
    echo "berkeleydb_install: missing ${db_prefix}/include/db.h (install libdb in Conda)" >&2
    return 1
  fi

  # Point build and runtime search paths to Berkeley DB (e.g. Conda `libdb`)
  export BERKELEYDB_DIR="$db_prefix"
  export BERKELEYDB_INCDIR="${db_prefix}/include"
  export BERKELEYDB_LIBDIR="${db_prefix}/lib"
  export YES_I_HAVE_THE_RIGHT_TO_USE_THIS_BERKELEY_DB_VERSION=1

  case ":${CPATH:-}:" in
    *":${db_prefix}/include:"*) ;;
    *) export CPATH="${db_prefix}/include${CPATH:+:${CPATH}}" ;;
  esac

  case ":${LIBRARY_PATH:-}:" in
    *":${db_prefix}/lib:"*) ;;
    *) export LIBRARY_PATH="${db_prefix}/lib${LIBRARY_PATH:+:${LIBRARY_PATH}}" ;;
  esac

  case ":${LD_LIBRARY_PATH:-}:" in
    *":${db_prefix}/lib:"*) ;;
    *) export LD_LIBRARY_PATH="${db_prefix}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}" ;;
  esac

  python -m pip install --no-build-isolation --no-cache-dir "berkeleydb"
}
# ============================================================================
##                                                                     The END
# ============================================================================