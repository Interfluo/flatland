#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# Command-line interface: flags, validation, and the quality of diagnostics.
#
# The standing contract is that NO input reaches the engine unvalidated. Every
# rejection must be a clean exit with a message naming the offending option —
# never a crash, and never a raw standard-library exception string such as
# "vector::_M_default_append" or "stoul", which tells a user nothing.
#
# Sourced by run_tests.sh with BIN, ROOT and TMP already set.

CUBE="$ROOT/examples/cube_area/cube.obj"
C="$TMP/cli"; mkdir -p "$C"

sect "Flags and basic interface"

"$BIN" --help >/dev/null 2>&1; expect_ok $? "--help exits 0"
"$BIN" -h     >/dev/null 2>&1; expect_ok $? "-h exits 0"
"$BIN"        >/dev/null 2>&1; expect_fail $? "no arguments is an error"

# Long forms must be honoured, not silently ignored.
P=$("$BIN" "$CUBE" --view 1 0 0 --res 0.1 -j 2>/dev/null | grep -m1 '"pixels"' | grep -oE '[0-9]+' | tail -1)
equal "$P" "100" "--res long flag honoured (100 px)"

"$BIN" "$CUBE" -v 1 0 0 -p double -j >/dev/null 2>&1; expect_ok $? "double precision runs"
"$BIN" "$CUBE" -v 1 0 0 -p quad   >/dev/null 2>&1; expect_fail $? "invalid precision rejected"
"$BIN" "$CUBE" -v 1 0 0 --bogus   >/dev/null 2>&1; expect_fail $? "unknown flag rejected"
"$BIN" "$CUBE" "$CUBE" -v 1 0 0   >/dev/null 2>&1; expect_fail $? "second positional argument rejected"
"$BIN" "$C/nope.obj" -v 1 0 0     >/dev/null 2>&1; expect_fail $? "missing mesh rejected"
"$BIN" "$CUBE"                    >/dev/null 2>&1; expect_fail $? "no views specified is an error"
"$BIN" "$CUBE" -v 1 0             >/dev/null 2>&1; expect_fail $? "incomplete -v rejected"
"$BIN" "$CUBE" -a 45              >/dev/null 2>&1; expect_fail $? "incomplete -a rejected"
"$BIN" "$CUBE" -r                 >/dev/null 2>&1; expect_fail $? "flag with a missing argument rejected"

sect "View-direction validation"

"$BIN" "$CUBE" -v 0 0 0 >/dev/null 2>&1; expect_fail $? "zero view vector rejected"
for bad_dir in "inf 0 1" "nan 0 1" "1 0 inf"; do
    # shellcheck disable=SC2086
    "$BIN" "$CUBE" -v $bad_dir >/dev/null 2>"$C/err"
    expect_fail $? "non-finite view direction rejected ($bad_dir)"
done
"$BIN" "$CUBE" -a nan 0 >/dev/null 2>&1; expect_fail $? "non-finite angle rejected"

sect "Resolution validation"

"$BIN" "$CUBE" -v 1 0 0 -r 0    >/dev/null 2>&1; expect_fail $? "zero resolution rejected"
"$BIN" "$CUBE" -v 1 0 0 -r -1   >/dev/null 2>&1; expect_fail $? "negative resolution rejected"
"$BIN" "$CUBE" -v 1 0 0 -r abc  >/dev/null 2>&1; expect_fail $? "non-numeric resolution rejected"

# NaN and infinity slip past a bare `res <= 0` test, then surface deep inside
# the rasterizer as an allocator error.
for bad_res in nan inf -inf; do
    "$BIN" "$CUBE" -v 1 0 0 -r "$bad_res" >/dev/null 2>"$C/err"
    expect_fail $? "non-finite resolution rejected ($bad_res)"
    stderr_has   "$C/err" "resolution" "...and the message names the resolution ($bad_res)"
    stderr_lacks "$C/err" "_M_default_append\|bad_alloc\|length_error" \
                 "...without leaking an allocator exception ($bad_res)"
done

# A resolution fine enough to demand an impossible raster is a user error and
# deserves a user-facing message, not a signed-overflow crash.
"$BIN" "$CUBE" -v 1 0 0 -r 1e-10 >/dev/null 2>"$C/err"
expect_fail   $? "absurdly fine resolution rejected"
stderr_has    "$C/err" "resolution\|pixel\|too\|large" "...and the message explains why"
stderr_lacks  "$C/err" "_M_default_append\|bad_alloc\|length_error" \
              "...without leaking an allocator exception"

sect "Thread-count validation"

"$BIN" "$CUBE" -v 1 0 0 -t 4 -j >/dev/null 2>&1; expect_ok $? "-t 4 accepted"
"$BIN" "$CUBE" -v 1 0 0 -t 0 -j >/dev/null 2>&1; expect_ok $? "-t 0 means auto"

# A big batch plus an unbounded -t used to try to spawn tens of thousands of
# threads; the std::thread constructor throws, unwinds through a vector holding
# joinable threads, and the process aborts. expect_fail treats SIGABRT as a
# failure precisely so this cannot pass by being non-zero.
python3 -c "open('$C/huge.txt','w').write('0 0 1 0.5\n'*4000)"
for bad_t in -1 100000 2.5 abc; do
    "$BIN" "$CUBE" -b "$C/huge.txt" -t "$bad_t" >/dev/null 2>"$C/err"
    expect_fail $? "invalid thread count rejected (-t $bad_t)"
done

# A legitimate request must still not be reported dishonestly: the banner should
# describe the workers actually started, which is capped by the view count.
"$BIN" "$CUBE" -v 1 0 0 -t 32 2>"$C/banner" >/dev/null
grep -qE 'on 1 thread' "$C/banner" \
    && ok "thread banner reports the workers actually used" \
    || bad "thread banner reports the request, not reality ($(tr -d '\n' < "$C/banner" | cut -c1-70))"

sect "Angle views"

"$BIN" "$CUBE" -a 0 0 -r 0.01 -j 2>/dev/null > "$C/a.json"
near "$(jget "$C/a.json" area)" 1.0 0.01 "angle view (-a 0 0 == +X) area == 1"
"$BIN" "$CUBE" -a 90 0 -r 0.01 -j 2>/dev/null > "$C/a90.json"
near "$(jget "$C/a90.json" area)" 1.0 0.01 "angle view (-a 90 0 == +Y) area == 1"
