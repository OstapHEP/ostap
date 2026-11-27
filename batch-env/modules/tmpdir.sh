#!/usr/bin/env bash
# ==============================================================================
# modules/tmpdir.sh
# Temporary workspace management and signal-safe cleanup for batch jobs.
#
# Exports:
#   BATCH_TMP_DIR - Path to the root temporary directory created for this job.
#
# Environment flags:
#   TMPDIR        - Optional base directory for temporary storage (defaults to /tmp).
#   KEEP_TMP      - Set to "1" to skip deletion on exit for debugging purposes.
# ==============================================================================
set -euo pipefail

# Global registry of all base temporary directories created in this shell session.
BATCH_TMP_DIRS=()

# Guard flag ensuring cleanup signal traps are registered only once.
BATCH_CLEANUP_REGISTERED=0

# ============================================================================
# tmp_trap_add <command> [signal] [mode]
#
# Safely appends or prepends a command to an existing shell signal trap without
# overwriting handlers previously registered by the caller or other modules.
#
# Arguments:
#   $1 - Shell command string to add to the trap.
#   $2 - Signal name (default: "EXIT").
#   $3 - Insertion order: "prepend" (default, runs before existing handlers)
#        or "append" (runs after existing handlers).
# ============================================================================
tmp_trap_add() {
  local new_cmd="${1:?tmp_trap_add: command string is required}"
  local signal="${2:-EXIT}"
  local mode="${3:-prepend}"

  local trap_stmt old_cmd=""
  trap_stmt=$(trap -p "$signal")

  # Extract the existing command string from `trap -p` output (format: trap -- 'cmd' SIGNAL).
  if [[ -n "$trap_stmt" ]]; then
    local raw="${trap_stmt#trap -- }"
    raw="${raw% "$signal"}"
    eval "old_cmd=$raw"
  fi

  local combined
  if [[ -n "$old_cmd" ]]; then
    if [[ "$mode" == "append" ]]; then
      combined="${old_cmd}; ${new_cmd}"
    else
      combined="${new_cmd}; ${old_cmd}"
    fi
  else
    combined="$new_cmd"
  fi

  trap -- "$combined" "$signal"
}

# ============================================================================
# batch_tmp_init
#
# Creates a unique root temporary directory under ${TMPDIR:-/tmp}, exports its
# path as BATCH_TMP_DIR, and registers cleanup hooks on EXIT, INT, and TERM.
# ============================================================================
batch_tmp_init() {
  local base_dir
  base_dir=$(mktemp -d -t "ostap-$(date '+%Y-%b-%d')-OSTAP-BATCH-XXXXXX")
  BATCH_TMP_DIRS+=("$base_dir")
  export BATCH_TMP_DIR="$base_dir"

  # Register cleanup handlers once per shell process:
  #   - EXIT: triggers on normal completion or any error under `set -e`.
  #   - INT:  triggers on Ctrl+C (exits with standard 128 + 2 = 130).
  #   - TERM: triggers when batch scheduler terminates the job (exits with 128 + 15 = 143).
  if [[ "$BATCH_CLEANUP_REGISTERED" -eq 0 ]]; then
    tmp_trap_add "batch_cleanup_all" EXIT
    tmp_trap_add "batch_cleanup_all; exit 130" INT
    tmp_trap_add "batch_cleanup_all; exit 143" TERM
    BATCH_CLEANUP_REGISTERED=1
  fi
}

# ============================================================================
# batch_make_dir [subdir_name]
#
# Creates a named subdirectory inside BATCH_TMP_DIR and prints its full path to
# stdout. Fails immediately if batch_tmp_init has not been called yet.
#
# Arguments:
#   $1 - Subdirectory name inside BATCH_TMP_DIR (default: "sub").
#
# Outputs:
#   Prints the created directory path to stdout.
# ============================================================================
batch_make_dir() {
  local name="${1:-sub}"
  local target="${BATCH_TMP_DIR:?batch_tmp_init has not been called}/$name"
  mkdir -p "$target"
  printf '%s\n' "$target"
}

# ============================================================================
# batch_cleanup_all
#
# Removes all temporary directories tracked in BATCH_TMP_DIRS unless KEEP_TMP=1
# is set in the environment.
# ============================================================================
batch_cleanup_all() {
  if [[ "${KEEP_TMP:-0}" == "1" ]]; then
    echo "[cleanup] KEEP_TMP=1, preserving temporary directories: ${BATCH_TMP_DIRS[*]:-}" >&2
    return 0
  fi

  local d
  for d in "${BATCH_TMP_DIRS[@]:-}"; do
    if [[ -n "$d" && -d "$d" ]]; then
      rm -rf "$d"
    fi
  done
}

# ============================================================================
##                                                                     The END
# ============================================================================