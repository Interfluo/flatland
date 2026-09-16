#!/usr/bin/env bash
#
# FlatLand test suite — dependency-free (bash + awk + a stock python3 for
# generating binary STL fixtures).
#
# Usage:
#   tests/run_tests.sh [path-to-flatland-binary] [suite ...]
#
# With no binary, builds one with `make`. With no suite names, runs them all.
# Suites: geometry parsing cli batch
#
# Examples:
#   tests/run_tests.sh                       # build, run everything
#   tests/run_tests.sh ./flatland            # test an existing binary
#   tests/run_tests.sh ./flatland geometry   # one suite
#   CXX=clang++ tests/run_tests.sh           # build with a specific compiler

set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

ALL_SUITES="geometry parsing cli batch"

# First argument is the binary only if it looks like a path, not a suite name.
BIN=""
case "${1:-}" in
    "")                    ;;
    geometry|parsing|cli|batch) ;;
    *) BIN="$1"; shift     ;;
esac
SUITES="${*:-$ALL_SUITES}"

if [ -z "$BIN" ]; then
    echo "Building with make ..."
    make -C "$ROOT" >/dev/null || { echo "BUILD FAILED"; exit 1; }
    BIN="$ROOT/flatland"
fi

if [ ! -x "$BIN" ]; then
    echo "error: '$BIN' is not an executable" >&2
    exit 1
fi

# Absolute path, so tests can chdir freely.
case "$BIN" in /*) ;; *) BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")" ;; esac

export ROOT TMP BIN

# shellcheck source=tests/lib.sh
. "$ROOT/tests/lib.sh"

echo "FlatLand test suite"
echo "  binary: $BIN"
echo "  suites: $SUITES"

for suite in $SUITES; do
    f="$ROOT/tests/test_$suite.sh"
    if [ ! -f "$f" ]; then
        bad "suite '$suite' does not exist"
        continue
    fi
    # shellcheck source=/dev/null
    . "$f"
done

echo
if [ "$fail" -eq 0 ]; then
    printf "%sResults: %d passed, 0 failed%s\n" "$C_PASS" "$pass" "$C_OFF"
else
    printf "%sResults: %d passed, %d failed%s\n" "$C_FAIL" "$pass" "$fail" "$C_OFF"
fi
[ "$fail" -eq 0 ]
