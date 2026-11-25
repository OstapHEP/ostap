#!/usr/bin/env bash
set -euo pipefail

# Prevent sourcing the script so environment variables do not leak into the caller shell
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    echo "ERROR: This script must be executed, not sourced (run './$(basename "${BASH_SOURCE[0]}")' or 'bash ...')." >&2
    return 1 2>/dev/null || exit 1
fi

tmp_root=""
_remove_dir_safely() {
    local target="$1"
    if [[ -z "${target}" ]] || [[ "${target}" == "/" ]] || [[ "${target}" == "." ]] || [[ "${target}" == ".." ]]; then
        return 0
    fi
    if [[ -d "${target}" ]]; then
        printf ">>> Cleaning up directory: %s\n" "${target}"
        rm -rf "${target}"
    fi
}

cleanup() {
    _remove_dir_safely "${tmp_root:-}"
    _remove_dir_safely "${OSTAP_BUILD:-}"
    _remove_dir_safely "${OSTAP_INSTALL:-}"
}

if declare -F trap_add >/dev/null 2>&1; then
    trap_add cleanup EXIT
else
    trap cleanup EXIT
fi


## OSTAP_SOURCE="/afs/cern.ch/user/i/ibelyaev/cmtuser/RELEASE/ostap"

OSTAP_SOURCE="${OSTAP_SOURCE:-}"
DEFAULT_GIT_REF="${OSTAP_BRANCH:-master}"

# ==========================================================================
## Valid source directory?
# ==========================================================================
if [[ -n "${OSTAP_SOURCE}"                ]] && \
   [[ -d "${OSTAP_SOURCE}"                ]] && \
   [[ -d "${OSTAP_SOURCE}/cmake"          ]] && \
   [[ -d "${OSTAP_SOURCE}/data"           ]] && \
   [[ -d "${OSTAP_SOURCE}/docs"           ]] && \
   [[ -d "${OSTAP_SOURCE}/examples"       ]] && \
   [[ -d "${OSTAP_SOURCE}/ostap"          ]] && \
   [[ -d "${OSTAP_SOURCE}/scripts"        ]] && \
   [[ -d "${OSTAP_SOURCE}/source"         ]] && \
   [[ -d "${OSTAP_SOURCE}/ReleaseNotes"   ]] && \
   [[ -f "${OSTAP_SOURCE}/CMakeLists.txt" ]] ; then

    printf "Use OSTAP_SOURCE      : %s\n" "${OSTAP_SOURCE}"
else

    if [[ -n "${OSTAP_SOURCE}" ]]; then
        git_ref="${OSTAP_SOURCE}"
    else
        git_ref="${DEFAULT_GIT_REF}"
    fi

    printf "No valid source directory '%s' is found! Try git-clone ref='%s'\n" "${OSTAP_SOURCE}" "${git_ref}"

    tmp_root="$(mktemp -d -t "ostap-$(date '+%Y-%b-%d')-OSTAP-GIT-XXXXXX")"
    OSTAP_SOURCE="${tmp_root}/ostap"

    git clone --branch "${git_ref}" --depth=1 git@github.com:OstapHEP/ostap.git "${OSTAP_SOURCE}"

    printf "Use OSTAP_SOURCE      : %s\n" "${OSTAP_SOURCE}"
fi


USE_LCG="dev3/latest"
USE_PLATFORM="x86_64-el9-clang22-opt"

if [[ -d "${OSTAP_SOURCE}/.aux"            ]] && \
   [[ -f "${OSTAP_SOURCE}/.aux/.build-lcg" ]] ; then

    OSTAP_BUILD=TMP \
        OSTAP_INSTALL=TMP \
        OSTAP_SOURCE="${OSTAP_SOURCE}" \
        source "${OSTAP_SOURCE}/.aux/.build-lcg" "${USE_LCG}" "${USE_PLATFORM}"
else 
    echo "ERROR: cannot find .aux/.build-lcg script" >&2
    return 2 2>/dev/null || exit 2
fi

# ==========================================================================
## Use it!!! 
# ==========================================================================


# =============================================================================
##                                                                      The END
# =============================================================================

