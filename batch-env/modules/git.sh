#!/usr/bin/env bash
# ==============================================================================
# modules/git.sh
# Lightweight Git helper for shallow-cloning repositories into BATCH_TMP_DIR.
#
# Exports:
#   BATCH_REPO_DIR - Path to the most recently cloned repository directory.
# ==============================================================================
set -euo pipefail

# ==============================================================================
# repo_clone <repo_url> [branch_or_tag] [dest_subdir]
#
# Performs a shallow clone (--depth 1) of a Git repository into a subdirectory
# inside BATCH_TMP_DIR. Git progress output is redirected to stderr so stdout
# returns only the destination path.
#
# Arguments:
#   $1 - Git repository URL (HTTPS or SSH).
#   $2 - Optional branch name or tag to check out (defaults to remote HEAD).
#   $3 - Subdirectory name inside BATCH_TMP_DIR (default: "repo").
#
# Outputs:
#   Prints the full path to the cloned repository to stdout.
# ============================================================================
repo_clone() {
  local url="${1:-https://github.com/OstapHEP/ostap.git}"
  local ref="${2:-master}"
  local name="${3:-ostap}"
  local dest
  dest=$(batch_make_dir "$name")

  local args=(clone --depth 1)
  if [[ -n "$ref" ]]; then
    args+=(--branch "$ref")
  fi
  args+=("$url" "$dest")

  git "${args[@]}" >&2
  export BATCH_REPO_DIR="$dest"
  printf '%s\n' "$dest"
}

# ============================================================================
##                                                                     The END 
# ============================================================================