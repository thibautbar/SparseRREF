#!/usr/bin/env bash
#
# Build script for libsparse_rref_inc — the incremental sparse RREF shared library.
#
# Output: libsparse_rref_inc.{dylib,so} in the current directory.
#
# Env vars:
#   CXX            C++ compiler (default: g++-15 on macOS arm64, g++ elsewhere)
#   FLINT_PREFIX   FLINT install prefix (default: /opt/homebrew on macOS, /usr/local elsewhere)
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

if [[ -z "${FLINT_PREFIX:-}" ]]; then
    if [[ "$UNAME" == "Darwin" ]]; then
        FLINT_PREFIX="/opt/homebrew"
    else
        FLINT_PREFIX="/usr/local"
    fi
fi
GMP_PREFIX="${GMP_PREFIX:-$FLINT_PREFIX}"

if [[ "$UNAME" == "Darwin" ]]; then
    LIB_EXT="dylib"
else
    LIB_EXT="so"
fi

OUT="libsparse_rref_inc.${LIB_EXT}"

echo "Building $OUT with $CXX (FLINT=$FLINT_PREFIX, GMP=$GMP_PREFIX)"

"$CXX" sprref_incremental.cpp \
    -fPIC -shared -O3 -std=c++20 \
    -I"$FLINT_PREFIX/include" -I"$GMP_PREFIX/include" \
    -L"$FLINT_PREFIX/lib"     -L"$GMP_PREFIX/lib" \
    -lflint -lgmp -ltbb \
    -o "$OUT"

echo "Built: $(pwd)/$OUT"
