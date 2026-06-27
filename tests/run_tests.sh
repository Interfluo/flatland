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

echo "== Mesh formats (STL) =="
# Build ASCII + binary STL of the cube and confirm both yield area 1.
python3 - "$CUBE" "$TMP/cube_ascii.stl" "$TMP/cube_bin.stl" <<'PY'
import sys, struct
obj, ascii_out, bin_out = sys.argv[1:4]
verts=[]; faces=[]
for line in open(obj):
    p=line.split()
    if not p: continue
    if p[0]=="v": verts.append(tuple(map(float,p[1:4])))
    elif p[0]=="f":
        idx=[int(t.split("/")[0])-1 for t in p[1:]]
        for k in range(2,len(idx)): faces.append((idx[0],idx[k-1],idx[k]))
def nrm(a,b,c):
    ux,uy,uz=(b[0]-a[0],b[1]-a[1],b[2]-a[2]); vx,vy,vz=(c[0]-a[0],c[1]-a[1],c[2]-a[2])
    return (uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx)
with open(ascii_out,"w") as f:
    f.write("solid cube\n")
    for a,b,c in faces:
        f.write("facet normal %f %f %f\nouter loop\n"%nrm(verts[a],verts[b],verts[c]))
        for vi in (a,b,c): f.write("vertex %f %f %f\n"%verts[vi])
        f.write("endloop\nendfacet\n")
    f.write("endsolid cube\n")
with open(bin_out,"wb") as f:
    f.write(b"\0"*80); f.write(struct.pack("<I",len(faces)))
    for a,b,c in faces:
        f.write(struct.pack("<3f",*nrm(verts[a],verts[b],verts[c])))
        for vi in (a,b,c): f.write(struct.pack("<3f",*verts[vi]))
        f.write(struct.pack("<H",0))
PY
"$BIN" "$TMP/cube_ascii.stl" -v 1 0 0 -r 0.01 -j 2>/dev/null > "$TMP/sa.json"
near "$(jget "$TMP/sa.json" area)" 1.0 0.01 "ASCII STL projected area == 1"
"$BIN" "$TMP/cube_bin.stl" -v 1 0 0 -r 0.01 -j 2>/dev/null > "$TMP/sb.json"
near "$(jget "$TMP/sb.json" area)" 1.0 0.01 "binary STL projected area == 1"

echo "== Field matrix (time series in one file) =="
# 26 vertex rows, 3 columns with values 1,2,3 -> selecting @c gives average c+1.
awk 'BEGIN{for(i=0;i<26;i++) print 1.0, 2.0, 3.0}' > "$TMP/matrix.txt"
"$BIN" "$CUBE" -v 1 0 0 -r 0.02 -d "$TMP/matrix.txt@0" -j 2>/dev/null > "$TMP/m0.json"
near "$(jget "$TMP/m0.json" average)" 1.0 0.001 "matrix column 0 selected"
"$BIN" "$CUBE" -v 1 0 0 -r 0.02 -d "$TMP/matrix.txt@2" -j 2>/dev/null > "$TMP/m2.json"
near "$(jget "$TMP/m2.json" average)" 3.0 0.001 "matrix column 2 selected"
"$BIN" "$CUBE" -v 1 0 0 -d "$TMP/matrix.txt@9" >/dev/null 2>&1; expect_fail $? "out-of-range column rejected"

echo "== Angle views =="
# -a 0 0 is +X; cube area along an axis is 1.
"$BIN" "$CUBE" -a 0 0 -r 0.01 -j 2>/dev/null > "$TMP/a.json"
near "$(jget "$TMP/a.json" area)" 1.0 0.01 "angle view (-a 0 0 == +X) area == 1"
# Batch angle line.
printf "a 90 0 0.01\n" > "$TMP/ang.txt"
"$BIN" "$CUBE" -b "$TMP/ang.txt" -j 2>/dev/null > "$TMP/ab.json"
near "$(jget "$TMP/ab.json" area)" 1.0 0.01 "batch angle line area == 1"

echo "== Field mode override =="
awk 'BEGIN{for(i=0;i<48;i++) print 7.0}' > "$TMP/face48.txt"
"$BIN" "$CUBE" -v 1 0 0 -r 0.02 -d "$TMP/face48.txt" --field-mode face -j 2>/dev/null > "$TMP/fm.json"
near "$(jget "$TMP/fm.json" average)" 7.0 0.001 "--field-mode face applied"
"$BIN" "$CUBE" -v 1 0 0 -d "$TMP/matrix.txt@0" --field-mode face >/dev/null 2>&1; expect_fail $? "--field-mode count mismatch rejected"

echo
echo "Results: $pass passed, $fail failed"
[ "$fail" -eq 0 ]
