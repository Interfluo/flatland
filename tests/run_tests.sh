#!/usr/bin/env bash
#
# FlatLand test suite — dependency-free (bash + awk only).
#
# Usage: tests/run_tests.sh [path-to-flatland-binary]
# If no binary is given, builds one from src/main.cpp with the system compiler.
#
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${1:-}"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

pass=0; fail=0

ok()   { printf "  \033[32mPASS\033[0m %s\n" "$1"; pass=$((pass+1)); }
bad()  { printf "  \033[31mFAIL\033[0m %s\n" "$1"; fail=$((fail+1)); }

# near <actual> <expected> <tol> <label>
near() {
    awk -v a="$1" -v e="$2" -v t="$3" 'BEGIN{d=a-e;if(d<0)d=-d;exit !(d<=t)}' \
        && ok "$4 ($1 ~= $2)" || bad "$4 (got $1, want $2)"
}
# expect_fail <label> ; reads $? of previous command
expect_fail() { [ "$1" -ne 0 ] && ok "$2 (exit $1)" || bad "$2 (expected non-zero exit)"; }

if [ -z "$BIN" ]; then
    BIN="$TMP/flatland"
    echo "Building $BIN ..."
    "${CXX:-c++}" -std=c++17 -O2 -pthread "$ROOT/src/main.cpp" -o "$BIN" || { echo "BUILD FAILED"; exit 1; }
fi

CUBE="$ROOT/examples/cube_area/cube.obj"
SPHERE="$ROOT/examples/sphere_areas/sphere_fine.obj"

# Helper: extract a numeric JSON field value for the first result.
jget() { grep -m1 "\"$2\"" "$1" | grep -oE '[-0-9.eE+]+' | tail -1; }

echo "== Geometry / math =="
# 1. Unit cube, axis view -> area 1
"$BIN" "$CUBE" -v 1 0 0 -r 0.005 -j 2>/dev/null > "$TMP/c.json"
near "$(jget "$TMP/c.json" area)" 1.0 0.01 "cube projected area == 1"

# 2. Unit sphere -> area pi from any angle
"$BIN" "$SPHERE" -v 1 1 1 -r 0.004 -j 2>/dev/null > "$TMP/s.json"
A=$(jget "$TMP/s.json" area)
near "$A" 3.14159 0.01 "sphere projected area == pi"

# 3. Constant field integral == const * area
seq 1 "$(grep -c '^v ' "$CUBE")" | awk '{print 5.0}' > "$TMP/const.txt"
"$BIN" "$CUBE" -v 1 0 0 -r 0.005 -d "$TMP/const.txt" -j 2>/dev/null > "$TMP/f.json"
near "$(jget "$TMP/f.json" average)"  5.0 0.001 "constant field average == 5"
near "$(jget "$TMP/f.json" integral)" 5.0 0.01 "field integral == const*area"

echo "== Interface / robustness =="
# 4. Long flag --res is honored (was silently ignored before v2)
P=$("$BIN" "$CUBE" --view 1 0 0 --res 0.1 -j 2>/dev/null | grep -m1 '"pixels"' | grep -oE '[0-9]+' | tail -1)
[ "$P" = "100" ] && ok "--res long flag honored (100 px)" || bad "--res long flag (got $P px)"

# 5. Unknown flag -> error
"$BIN" "$CUBE" -v 1 0 0 --bogus >/dev/null 2>&1; expect_fail $? "unknown flag rejected"
# 6. Missing mesh -> clean error, no crash
"$BIN" "$TMP/nope.obj" -v 1 0 0 >/dev/null 2>&1; expect_fail $? "missing mesh rejected"
# 7. Zero view vector -> error
"$BIN" "$CUBE" -v 0 0 0 >/dev/null 2>&1; expect_fail $? "zero view vector rejected"
# 8. -v with too few args -> error
"$BIN" "$CUBE" -v 1 0 >/dev/null 2>&1; expect_fail $? "incomplete -v rejected"
# 9. Field size mismatch -> error
echo "1.0" > "$TMP/bad.txt"
"$BIN" "$CUBE" -v 1 0 0 -d "$TMP/bad.txt" >/dev/null 2>&1; expect_fail $? "field size mismatch rejected"
# 9b. Out-of-bounds vertex index in OBJ -> error (not a crash)
printf "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 99\n" > "$TMP/oob.obj"
"$BIN" "$TMP/oob.obj" -v 1 0 0 >/dev/null 2>&1; expect_fail $? "out-of-bounds vertex index rejected"

echo "== Batch / time-series =="
# 10. Batch with relative data path resolves against the batch file's directory
mkdir -p "$TMP/case"
seq 1 "$(grep -c '^v ' "$CUBE")" | awk '{print 2.0}' > "$TMP/case/f.txt"
printf "1 0 0 0.01 f.txt\n0 1 0 0.01 f.txt\n" > "$TMP/case/ts.txt"
"$BIN" "$CUBE" -b "$TMP/case/ts.txt" -j 2>/dev/null > "$TMP/b.json"
N=$(grep -c '"idx"' "$TMP/b.json")
[ "$N" = "2" ] && ok "batch produced 2 results" || bad "batch result count (got $N)"
near "$(jget "$TMP/b.json" average)" 2.0 0.001 "batch relative data path resolved"

# 10b. Batch data file WITHOUT an explicit resolution must still apply (token sniff)
printf "1 0 0 f.txt\n" > "$TMP/case/nores.txt"
"$BIN" "$CUBE" -b "$TMP/case/nores.txt" -j 2>/dev/null > "$TMP/nr.json"
near "$(jget "$TMP/nr.json" average)" 2.0 0.001 "batch data file without resolution applied"

# 11. float and double precision both run; invalid precision rejected
"$BIN" "$CUBE" -v 1 0 0 -p double -j >/dev/null 2>&1 && ok "double precision runs" || bad "double precision"
"$BIN" "$CUBE" -v 1 0 0 -p quad  >/dev/null 2>&1; expect_fail $? "invalid precision rejected"

echo
echo "Results: $pass passed, $fail failed"
[ "$fail" -eq 0 ]
