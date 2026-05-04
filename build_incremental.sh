#!/usr/bin/env bash
#
# Build script for libsparse_rref_inc — the incremental sparse RREF shared library.
#
# Output: libsparse_rref_inc.{dylib,so} in the current directory.
#
# Env vars:
#   CXX            C++ compiler (default: g++-15 on macOS arm64, g++ elsewhere)
#   FLINT_PREFIX   FLINT install prefix. If unset, /opt/homebrew on macOS, else
#                  auto-detected: /usr if /usr/include/flint exists, otherwise
#                  /usr/local. Pass explicitly to override.
#   GMP_PREFIX     GMP install prefix (default: same as FLINT_PREFIX)
#
# The upstream code uses C++20 parallel STL for tensor ops which Apple Clang
# does not ship; we build with GCC-15 + libtbb on macOS to stay aligned with
# the standard build instructions.

set -euo pipefail

UNAME="$(uname -s)"

if [[ -z "${CXX:-}" ]]; then
    if [[ "$UNAME" == "Darwin" ]]; then
        CXX="g++-15"
    else
        CXX="g++"
    fi
fi

# FLINT prefix detection. On Linux distros that ship libflint-dev, headers are
# in /usr/include/flint and libs in a multiarch path that g++ already searches,
# so we leave -I/-L empty and rely on the compiler's defaults. On macOS arm64
# (Homebrew) headers/libs are in /opt/homebrew. /usr/local is the historical
# fallback for Linux source builds.
if [[ -z "${FLINT_PREFIX:-}" ]]; then
    if [[ "$UNAME" == "Darwin" ]]; then
        FLINT_PREFIX="/opt/homebrew"
    elif [[ -d "/usr/include/flint" ]]; then
        FLINT_PREFIX=""  # use compiler default search paths
    else
        FLINT_PREFIX="/usr/local"
    fi
fi
GMP_PREFIX="${GMP_PREFIX:-$FLINT_PREFIX}"

INCLUDE_FLAGS=""
LIB_FLAGS=""
if [[ -n "$FLINT_PREFIX" ]]; then
    INCLUDE_FLAGS="-I$FLINT_PREFIX/include"
    LIB_FLAGS="-L$FLINT_PREFIX/lib"
fi
if [[ -n "$GMP_PREFIX" && "$GMP_PREFIX" != "$FLINT_PREFIX" ]]; then
    INCLUDE_FLAGS="$INCLUDE_FLAGS -I$GMP_PREFIX/include"
    LIB_FLAGS="$LIB_FLAGS -L$GMP_PREFIX/lib"
fi

if [[ "$UNAME" == "Darwin" ]]; then
    LIB_EXT="dylib"
else
    LIB_EXT="so"
fi

OUT="libsparse_rref_inc.${LIB_EXT}"

echo "Building $OUT with $CXX (FLINT_PREFIX='${FLINT_PREFIX}', GMP_PREFIX='${GMP_PREFIX}')"

"$CXX" sprref_incremental.cpp \
    -fPIC -shared -O3 -std=c++20 \
    $INCLUDE_FLAGS $LIB_FLAGS \
    -lflint -lgmp -ltbb \
    -o "$OUT"

echo "Built: $(pwd)/$OUT"
