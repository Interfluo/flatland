#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# The Python binding: python/flatland/__init__.py, a pure-ctypes wrapper over the
# same C ABI the C test exercises.
#
# The binding's premise is zero dependencies, with NumPy used only if it happens
# to be installed. That is two code paths, and running only the one this machine
# happens to have would leave the other untested — so the Python test program is
# run TWICE: once as the environment is, and once with an import shim that makes
# `import numpy` fail. Neither run touches the installed packages.
#
# Sourced by run_tests.sh with ROOT and TMP already set.

sect "Python binding"

PY_BIN="${PYTHON:-python3}"
P="$TMP/python"; mkdir -p "$P"

if ! command -v "$PY_BIN" >/dev/null 2>&1; then
    # Not a failure: the binding is optional and the rest of the suite stands
    # on its own. Reported as a pass so the totals stay honest either way.
    ok "python3 is unavailable, so the Python binding suite is skipped"
    return 0 2>/dev/null || true
fi

if make -C "$ROOT" lib >"$P/build.log" 2>&1; then
    ok "libflatland builds for the Python binding"
else
    bad "libflatland builds for the Python binding ($(tail -3 "$P/build.log" | tr '\n' ' '))"
    return 0 2>/dev/null || true
fi

# An inherited PYTHONPATH is kept, so a numpy living in a venv or a non-default
# location is still seen by the "as installed" run.
BASE_PYPATH="$ROOT/python${PYTHONPATH:+:$PYTHONPATH}"

SHARED=""
for cand in "$ROOT/libflatland.so" "$ROOT/libflatland.dylib" "$ROOT/flatland.dll"; do
    [ -f "$cand" ] && SHARED="$cand" && break
done
[ -n "$SHARED" ] && ok "a shared library exists for ctypes to load" \
                 || bad "a shared library exists for ctypes to load"

# The binding must be loadable with no environment help at all: it looks next to
# the package and then at the repository root, which is where `make lib` writes.
( cd "$P" && PYTHONPATH="$BASE_PYPATH" "$PY_BIN" -c \
    'import flatland, sys; sys.exit(0 if flatland.library_path() else 1)' ) \
    >"$P/import.log" 2>&1
expect_ok $? "the binding finds libflatland with no environment hints"

# Syntax must stay inside the 3.8 subset even when a newer interpreter runs it.
"$PY_BIN" -c 'import ast, sys
for p in sys.argv[1:]:
    with open(p) as f:
        ast.parse(f.read(), p, feature_version=(3, 8))' \
    "$ROOT/python/flatland/__init__.py" "$ROOT/python/example.py" \
    "$ROOT/python/test_flatland.py" >"$P/syntax.log" 2>&1
expect_ok $? "the binding parses under Python 3.8 syntax rules"

# Replay the Python program's own PASS/FAIL lines into this suite's counters, so
# a failed assertion is reported here instead of as one opaque non-zero exit.
replay() {
    while IFS= read -r line; do
        case "$line" in
            "  PASS "*) pass=$((pass+1)) ;;
            "  FAIL "*) bad "${line#  FAIL }" ;;
        esac
    done < "$1"
}

# run_python <label> <expect-numpy> [extra PYTHONPATH prefix]
run_python() {
    local label="$1" expect="$2" prefix="${3:-}"
    local out="$P/$label.out"
    local pypath="$BASE_PYPATH"
    [ -n "$prefix" ] && pypath="$prefix:$pypath"
    ( cd "$P" && PYTHONPATH="$pypath" "$PY_BIN" "$ROOT/python/test_flatland.py" \
        "$ROOT" --expect-numpy "$expect" >"$out" 2>&1 )
    local rc=$?
    replay "$out"
    if [ "$rc" -eq 0 ]; then
        ok "the Python test program exits 0 ($label)"
    else
        bad "the Python test program exits $rc ($label)"
        tail -6 "$out" | sed 's/^/      /'
    fi
}

# --- run 1: the environment as it is ---------------------------------------
if "$PY_BIN" -c 'import numpy' >/dev/null 2>&1; then
    HAVE_NUMPY=yes
else
    HAVE_NUMPY=no
fi
printf "  (numpy importable here: %s)\n" "$HAVE_NUMPY"
run_python "as-installed" "$HAVE_NUMPY"

# --- run 2: numpy forced unavailable ---------------------------------------
# A sys.path shim, not an uninstall: a numpy.py that refuses to import shadows
# the real package for this run only, and leaves the interpreter untouched.
mkdir -p "$P/no_numpy_shim"
cat > "$P/no_numpy_shim/numpy.py" <<'SHIM'
raise ImportError("numpy is blocked for this FlatLand test run")
SHIM

if [ "$HAVE_NUMPY" = "yes" ]; then
    # Prove the shim actually bites before trusting the run it guards.
    ( cd "$P" && PYTHONPATH="$P/no_numpy_shim" "$PY_BIN" -c 'import numpy' ) >/dev/null 2>&1
    expect_fail $? 'the no-numpy shim makes import numpy fail'
fi
run_python "no-numpy" "no" "$P/no_numpy_shim"

# --- library discovery ------------------------------------------------------
sect "Python library discovery"

if [ -n "$SHARED" ]; then
    ( cd "$P" && FLATLAND_LIBRARY="$SHARED" PYTHONPATH="$ROOT/python" "$PY_BIN" -c \
        'import flatland, os, sys; sys.exit(0 if os.path.realpath(flatland.library_path()) == os.path.realpath(sys.argv[1]) else 1)' \
        "$SHARED" ) >"$P/override.log" 2>&1
    expect_ok $? "FLATLAND_LIBRARY overrides discovery"
fi

# A missing library must be a clear, actionable diagnostic — not an ImportError
# about ctypes internals and not a crash.
( cd "$P" && FLATLAND_LIBRARY="$P/definitely-not-here.so" PYTHONPATH="$ROOT/python" \
    "$PY_BIN" -c 'import flatland' ) >"$P/missing.log" 2>&1
expect_fail $? "a missing library is rejected"
stderr_has "$P/missing.log" "FLATLAND_LIBRARY" "...and the message names the override"
stderr_has "$P/missing.log" "make" "...and says how to build one"

# The same package copied somewhere with no library beside it and no loader hint
# must say what to do, not raise a bare OSError from ctypes.
mkdir -p "$P/iso"
cp -r "$ROOT/python/flatland" "$P/iso/flatland"
( cd "$P" && env -u FLATLAND_LIBRARY -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH \
    PYTHONPATH="$P/iso" "$PY_BIN" -c 'import flatland' ) >"$P/notfound.log" 2>&1
expect_fail  $? "an unfindable library is a clean diagnostic"
stderr_has   "$P/notfound.log" "make" "...naming the build command"
stderr_has   "$P/notfound.log" "FLATLAND_LIBRARY" "...and the environment override"
stderr_lacks "$P/notfound.log" "Traceback.*ctypes" "...without a bare ctypes error"

# Discovery with the repository root removed from the search: the loader's own
# path must still be consulted, which is what an installed binding relies on.
( cd "$P" && LD_LIBRARY_PATH="$ROOT:${LD_LIBRARY_PATH:-}" \
              DYLD_LIBRARY_PATH="$ROOT:${DYLD_LIBRARY_PATH:-}" \
              PYTHONPATH="$ROOT/python" "$PY_BIN" -c 'import flatland; flatland.angle_to_dir(0, 0)' \
) >"$P/ldpath.log" 2>&1
expect_ok $? "the binding loads with the library on LD_LIBRARY_PATH"

# --- the documented example -------------------------------------------------
sect "Python example"

( cd "$P" && PYTHONPATH="$BASE_PYPATH" "$PY_BIN" "$ROOT/python/example.py" "$ROOT" ) \
    >"$P/example.out" 2>&1
expect_ok $? "python/example.py runs end to end"
grep -qi "area" "$P/example.out" \
    && ok "...and reports a projected area" \
    || bad "...and reports a projected area ($(tail -2 "$P/example.out" | tr '\n' ' '))"
