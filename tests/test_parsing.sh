#!/usr/bin/env bash
#
# Mesh and field-file parsing: correct files must parse, malformed files must
# produce a clean diagnostic. The cardinal rule tested throughout is that input
# is never SILENTLY dropped — a skipped vertex row or field row shifts every
# index after it and turns a parse failure into a plausible wrong answer.
#
# Sourced by run_tests.sh with BIN, ROOT and TMP already set.

CUBE="$ROOT/examples/cube_area/cube.obj"
P="$TMP/parse"; mkdir -p "$P"

# 2-unit right triangle in the XY plane, CCW from +Z; true area 2.
# Seen from -Z (the camera looking along -Z sees the +Z face) with --no-cull.
write_tri() { printf 'v 0 0 0\nv 2 0 0\nv 0 2 0\nv 0 1 0\nf 1 2 3\n' > "$1"; }
tri_area() { "$BIN" "$1" -v 0 0 -1 -r 0.01 --no-cull -j 2>/dev/null > "$P/o.json"; jget "$P/o.json" area; }

sect "OBJ parsing"

write_tri "$P/plain.obj"
BASE=$(tri_area "$P/plain.obj")
near "$BASE" 2.01 0.02 "baseline OBJ parses"

# Free-form OBJ permits leading whitespace, and plenty of exporters indent.
# Indexing line[0]/line[1] instead of tokenizing drops these rows silently.
printf 'v 0 0 0\nv 2 0 0\n v 0 2 0\nv 0 1 0\nf 1 2 3\n' > "$P/indent.obj"
near "$(tri_area "$P/indent.obj")" "$BASE" 0.001 "OBJ vertex with a leading space parses"

printf 'v 0 0 0\nv 2 0 0\nv\t0 2 0\nv 0 1 0\nf 1 2 3\n' > "$P/tab.obj"
near "$(tri_area "$P/tab.obj")" "$BASE" 0.001 "OBJ vertex with a tab separator parses"

printf '\xef\xbb\xbfv 0 0 0\nv 2 0 0\nv 0 2 0\nv 0 1 0\nf 1 2 3\n' > "$P/bom.obj"
near "$(tri_area "$P/bom.obj")" "$BASE" 0.001 "OBJ with a UTF-8 BOM parses"

printf 'v 0 0 0\nv 2 0 0\nv 0 2 0\nv 0 1 0\n f 1 2 3\n' > "$P/findent.obj"
near "$(tri_area "$P/findent.obj")" "$BASE" 0.001 "OBJ face with a leading space parses"

# A malformed vertex must be an error, never a silent skip.
printf 'v 0 0 0\nv 2 0\nv 0 2 0\nv 0 1 0\nf 1 2 3\n' > "$P/short.obj"
"$BIN" "$P/short.obj" -v 0 0 -1 >"$P/out" 2>"$P/err"; expect_fail $? "OBJ vertex with too few coordinates rejected"
stderr_has "$P/err" "coordinate\|vertex\|line" "...and the message names the problem"

printf 'v 0 0 0\nv nan 0 0\nv 0 2 0\nf 1 2 3\n' > "$P/nanv.obj"
"$BIN" "$P/nanv.obj" -v 0 0 -1 >/dev/null 2>&1; expect_fail $? "OBJ with a non-finite coordinate rejected"

# These were already correct; they are pinned so the rewrite cannot regress them.
printf 'v 0 0 0\nv 2 0 0\nv 0 2 0\nv 0 1 0\nvn 0 0 1\nvt 0 0\nf 1//1 2//1 3//1\n' > "$P/slash.obj"
near "$(tri_area "$P/slash.obj")" "$BASE" 0.001 "OBJ f v//vn form and vn/vt lines"
printf 'v 0 0 0\nv 2 0 0\nv 0 2 0\nf -3 -2 -1\n' > "$P/neg.obj"
near "$(tri_area "$P/neg.obj")" 2.01 0.02 "OBJ negative relative indices"
printf 'v 0 0 0\r\nv 2 0 0\r\nv 0 2 0\r\nv 0 1 0\r\nf 1 2 3\r\n' > "$P/crlf.obj"
near "$(tri_area "$P/crlf.obj")" "$BASE" 0.001 "OBJ with CRLF line endings"
printf 'v 0 0 0\nv 2 0 0\nv 2 2 0\nv 0 2 0\nf 1 2 3 4\n' > "$P/quad.obj"
near "$(tri_area "$P/quad.obj")" 4.02 0.06 "OBJ quad is fan-triangulated"

sect "STL parsing"

# Both encodings of the same cube must give the analytic answer.
python3 - "$CUBE" "$P/cube_ascii.stl" "$P/cube_bin.stl" <<'PY'
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
"$BIN" "$P/cube_ascii.stl" -v 1 0 0 -r 0.01 -j 2>/dev/null > "$P/sa.json"
near "$(jget "$P/sa.json" area)" 1.0 0.01 "ASCII STL projected area == 1"
"$BIN" "$P/cube_bin.stl" -v 1 0 0 -r 0.01 -j 2>/dev/null > "$P/sb.json"
near "$(jget "$P/sb.json" area)" 1.0 0.01 "binary STL projected area == 1"

# A truncated binary STL must say so, not fall through to the ASCII reader and
# report the misleading "contains no vertices".
python3 -c "
import struct
d=b'hdr'.ljust(80,b'\0')+struct.pack('<I',1)+struct.pack('<3f',0,0,1)+struct.pack('<9f',0,0,0,2,0,0,0,2,0)+b'\0\0'
open('$P/trunc.stl','wb').write(d[:-10])"
"$BIN" "$P/trunc.stl" -v 0 0 -1 >/dev/null 2>"$P/err"; expect_fail $? "truncated binary STL rejected"
stderr_has "$P/err" "truncat\|incomplete\|size\|byte" "...and the message says the file is truncated"

sect "Field-file parsing"

# 4 vertices, 3 faces — so a single dropped row silently reinterprets a node
# field as a face field and returns a plausible but wrong answer.
printf 'v 0 0 0\nv 1 0 0\nv 0 1 0\nv -1 -1 0\nf 1 2 3\nf 1 3 4\nf 1 4 2\n' > "$P/fan.obj"
fan_avg() { "$BIN" "$P/fan.obj" -v 0 0 -1 -r 0.05 --no-cull -d "$1" -j 2>/dev/null > "$P/fa.json"; jget "$P/fa.json" average; }

printf '10\n20\n30\n40\n' > "$P/node.txt"
NODE=$(fan_avg "$P/node.txt")
# 23.3 is the node answer; the face answer for this mesh is 29.9, so the
# tolerance here is deliberately far narrower than the gap between them.
near "$NODE" 23.31 0.1 "4-row node field parses as a node field"

printf '\xef\xbb\xbf10\n20\n30\n40\n' > "$P/node_bom.txt"
near "$(fan_avg "$P/node_bom.txt")" "$NODE" 0.001 "BOM'd field file parses identically"

# A row that carries no number is a defect in the data, not a blank line.
printf '10\nN/A\n30\n40\n' > "$P/na.txt"
"$BIN" "$P/fan.obj" -v 0 0 -1 -d "$P/na.txt" >/dev/null 2>"$P/err"; expect_fail $? "field row with no numeric value rejected"
stderr_has "$P/err" "line\|row\|numeric" "...and the message locates the row"

# Trailing junk means the row was not understood.
printf '10 abc\n20 abc\n30 abc\n40 abc\n' > "$P/junk.txt"
"$BIN" "$P/fan.obj" -v 0 0 -1 -d "$P/junk.txt" >/dev/null 2>&1; expect_fail $? "field row with trailing junk rejected"

# Blank lines and comments remain legal, and the row ordinal in a ragged-row
# error must be the FILE LINE, or it points the user at the wrong line.
printf '# header\n\n1 1\n2 2\n3 3 3\n4 4\n' > "$P/ragged.txt"
"$BIN" "$P/fan.obj" -v 0 0 -1 -d "$P/ragged.txt" >/dev/null 2>"$P/err"; expect_fail $? "ragged field row rejected"
stderr_has "$P/err" "line 5" "...and the error names file line 5, not row 3"

# Column selection, including the diagnostic for an unusable suffix.
awk 'BEGIN{for(i=0;i<26;i++) print 1.0, 2.0, 3.0}' > "$P/matrix.txt"
"$BIN" "$CUBE" -v 1 0 0 -r 0.02 -d "$P/matrix.txt@0" -j 2>/dev/null > "$P/m0.json"
near "$(jget "$P/m0.json" average)" 1.0 0.001 "matrix column 0 selected"
"$BIN" "$CUBE" -v 1 0 0 -r 0.02 -d "$P/matrix.txt@2" -j 2>/dev/null > "$P/m2.json"
near "$(jget "$P/m2.json" average)" 3.0 0.001 "matrix column 2 selected"
"$BIN" "$CUBE" -v 1 0 0 -d "$P/matrix.txt@9" >/dev/null 2>&1; expect_fail $? "out-of-range column rejected"
"$BIN" "$CUBE" -v 1 0 0 -d "$P/matrix.txt@99999999999999999999" >/dev/null 2>"$P/err"
expect_fail $? "absurd column index rejected"
stderr_lacks "$P/err" "^Error: stoul$" "...without leaking a bare 'stoul' from the standard library"

# A data file whose name genuinely ends in @<digits> must still be reachable.
printf '1\n2\n3\n' > "$P/run@1"
"$BIN" "$P/fan.obj" -v 0 0 -1 -r 0.05 --no-cull -d "$P/run@1" -j 2>/dev/null > "$P/at.json"
near "$(jget "$P/at.json" average)" 2.0 0.5 "a data file named '...@1' resolves to the file"

sect "Field mode"

awk 'BEGIN{for(i=0;i<48;i++) print 7.0}' > "$P/face48.txt"
"$BIN" "$CUBE" -v 1 0 0 -r 0.02 -d "$P/face48.txt" --field-mode face -j 2>/dev/null > "$P/fm.json"
near "$(jget "$P/fm.json" average)" 7.0 0.001 "--field-mode face applied"
"$BIN" "$CUBE" -v 1 0 0 -d "$P/matrix.txt@0" --field-mode face >/dev/null 2>&1
expect_fail $? "--field-mode count mismatch rejected"
echo "1.0" > "$P/bad.txt"
"$BIN" "$CUBE" -v 1 0 0 -d "$P/bad.txt" >/dev/null 2>&1; expect_fail $? "field size mismatch rejected"
