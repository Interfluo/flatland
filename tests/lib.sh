#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# Shared helpers for the FlatLand test suite.
#
# Every test file sources this and then uses the assertion helpers below. The
# driver (run_tests.sh) exports BIN, ROOT and TMP before sourcing anything.

# --- counters ---------------------------------------------------------------
pass=0
fail=0

if [ -t 1 ]; then
    C_PASS=$'\033[32m'; C_FAIL=$'\033[31m'; C_SECT=$'\033[1m'; C_OFF=$'\033[0m'
else
    C_PASS=''; C_FAIL=''; C_SECT=''; C_OFF=''
fi

ok()   { printf "  %sPASS%s %s\n" "$C_PASS" "$C_OFF" "$1"; pass=$((pass+1)); }
bad()  { printf "  %sFAIL%s %s\n" "$C_FAIL" "$C_OFF" "$1"; fail=$((fail+1)); }
sect() { printf "\n%s== %s ==%s\n" "$C_SECT" "$1" "$C_OFF"; }

# --- numeric assertions -----------------------------------------------------

# near <actual> <expected> <tol> <label>   — absolute tolerance
near() {
    if [ -z "${1:-}" ]; then bad "$4 (got empty output, want $2)"; return; fi
    awk -v a="$1" -v e="$2" -v t="$3" 'BEGIN{d=a-e;if(d<0)d=-d;exit !(d<=t)}' \
        && ok "$4 ($1 ~= $2)" || bad "$4 (got $1, want $2)"
}

# rel_near <actual> <expected> <rel_tol> <label>   — relative tolerance
rel_near() {
    if [ -z "${1:-}" ]; then bad "$4 (got empty output, want $2)"; return; fi
    awk -v a="$1" -v e="$2" -v t="$3" \
        'BEGIN{d=a-e;if(d<0)d=-d;s=e;if(s<0)s=-s;if(s==0)s=1;exit !(d/s<=t)}' \
        && ok "$4 ($1 ~= $2)" || bad "$4 (got $1, want $2)"
}

# equal <actual> <expected> <label>   — exact string match
equal() {
    [ "$1" = "$2" ] && ok "$3 ($1)" || bad "$3 (got '$1', want '$2')"
}

# --- process assertions -----------------------------------------------------

# expect_fail <exit_code> <label>  — asserts a NON-ZERO exit that is not a crash.
# A segfault (139) or abort (134) is a failure even though it is non-zero: the
# contract is a clean diagnostic, never a crash.
expect_fail() {
    case "$1" in
        0)        bad "$2 (expected non-zero exit, got 0)" ;;
        134)      bad "$2 (ABORTED - expected a clean error)" ;;
        139)      bad "$2 (SEGFAULT - expected a clean error)" ;;
        13[0-9])  bad "$2 (killed by signal $(( $1 - 128 )) - expected a clean error)" ;;
        *)        ok  "$2 (exit $1)" ;;
    esac
}

# expect_ok <exit_code> <label>
expect_ok() { [ "$1" -eq 0 ] && ok "$2" || bad "$2 (exit $1, expected 0)"; }

# no_crash <exit_code> <label>  — any clean exit, but must not be killed
no_crash() {
    case "$1" in
        13[0-9]) bad "$2 (killed by signal $(( $1 - 128 )))" ;;
        *)       ok  "$2 (exit $1)" ;;
    esac
}

# stderr_has <file> <pattern> <label>  — diagnostic quality check
stderr_has() {
    if grep -qi -- "$2" "$1"; then
        ok "$3"
    else
        bad "$3 (message was: $(tr -d '\n' < "$1" | cut -c1-90))"
    fi
}

# stderr_lacks <file> <pattern> <label>
stderr_lacks() {
    if grep -qi -- "$2" "$1"; then
        bad "$3 (leaked internal detail: $(tr -d '\n' < "$1" | cut -c1-90))"
    else
        ok "$3"
    fi
}

# --- JSON extraction --------------------------------------------------------

# jget <file> <key>  — numeric value of the first occurrence of "key"
jget() {
    local raw; raw="$(jraw "$1" "$2")"
    [ "$raw" = "null" ] && return 0          # null is not a number; yield nothing
    printf '%s' "$raw" | grep -oE '[-0-9.eE+]+' | tail -1
}

# jraw <file> <key>  — raw token (so `null` is distinguishable from 0)
jraw() { grep -m1 "\"$2\"" "$1" | sed -E 's/.*"'"$2"'"[[:space:]]*:[[:space:]]*//; s/,[[:space:]]*$//'; }

# jget_n <file> <key> <n>  — numeric value of the n-th occurrence (1-based)
jget_n() {
    local raw; raw="$(jraw_n "$1" "$2" "$3")"
    [ "$raw" = "null" ] && return 0
    printf '%s' "$raw" | grep -oE '[-0-9.eE+]+' | tail -1
}

# jraw_n <file> <key> <n>  — raw token of the n-th occurrence (1-based)
jraw_n() { grep "\"$2\"" "$1" | sed -n "$3p" | sed -E 's/.*"'"$2"'"[[:space:]]*:[[:space:]]*//; s/,[[:space:]]*$//'; }

# --- fixtures ---------------------------------------------------------------

# A closed unit box spanning [0,1]^3 with correct CCW outward normals.
# Used wherever a test needs an unambiguous near/far surface.
write_box() {
    cat > "$1" <<'OBJ'
v 0 0 0
v 1 0 0
v 1 1 0
v 0 1 0
v 0 0 1
v 1 0 1
v 1 1 1
v 0 1 1
f 1 4 3
f 1 3 2
f 5 6 7
f 5 7 8
f 1 2 6
f 1 6 5
f 2 3 7
f 2 7 6
f 3 4 8
f 3 8 7
f 4 1 5
f 4 5 8
OBJ
}

# A single unit square in the XY plane, CCW seen from +Z.
write_square() {
    printf 'v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nf 1 2 3\nf 1 3 4\n' > "$1"
}
