#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# Geometry, projection and rasterization correctness.
#
# Sourced by run_tests.sh with BIN, ROOT and TMP already set.

CUBE="$ROOT/examples/cube_area/cube.obj"
SPHERE="$ROOT/examples/sphere_areas/sphere_fine.obj"
G="$TMP/geom"; mkdir -p "$G"

sect "Analytical validation"

# A unit cube's projection along an axis is exactly 1.
"$BIN" "$CUBE" -v 1 0 0 -r 0.005 -j 2>/dev/null > "$G/c.json"
near "$(jget "$G/c.json" area)" 1.0 0.01 "cube projected area == 1"

# A unit sphere projects to pi from every direction.
"$BIN" "$SPHERE" -v 1 1 1 -r 0.004 -j 2>/dev/null > "$G/s.json"
near "$(jget "$G/s.json" area)" 3.14159 0.01 "sphere projected area == pi"

# The area integral of a constant field is that constant times the area.
seq 1 "$(grep -c '^v ' "$CUBE")" | awk '{print 5.0}' > "$G/const5.txt"
"$BIN" "$CUBE" -v 1 0 0 -r 0.005 -d "$G/const5.txt" -j 2>/dev/null > "$G/f.json"
near "$(jget "$G/f.json" average)"  5.0 0.001 "constant field average == 5"
near "$(jget "$G/f.json" integral)" 5.0 0.01 "field integral == const * area"

sect "Backface culling agrees with the depth test"

# Culling is an optimisation, not a change of subject: on a closed mesh it must
# resolve to the same surface the z-buffer picks with culling disabled. A field
# equal to the vertex x-coordinate distinguishes the near surface (x = -0.5 on
# this cube) from the far one; a CONSTANT field cannot, which is why this class
# of bug survived the original suite.
awk '/^v /{print $2}' "$CUBE" > "$G/fx.txt"
"$BIN" "$CUBE" -v 1 0 0 -r 0.01 -d "$G/fx.txt" -p double -j 2>/dev/null > "$G/cull.json"
"$BIN" "$CUBE" -v 1 0 0 -r 0.01 -d "$G/fx.txt" -p double --no-cull -j 2>/dev/null > "$G/nocull.json"
near "$(jget "$G/cull.json" average)" -0.5 0.001 "culled view resolves to the NEAR surface"
near "$(jget "$G/cull.json" average)" "$(jget "$G/nocull.json" average)" 0.001 \
     "culling agrees with --no-cull on a closed mesh"

# Same check from the opposite direction, so a sign flip cannot pass both.
"$BIN" "$CUBE" -v -1 0 0 -r 0.01 -d "$G/fx.txt" -p double -j 2>/dev/null > "$G/cull2.json"
near "$(jget "$G/cull2.json" average)" 0.5 0.001 "culled view from -X resolves to its near surface"

sect "Numerical conditioning"

# A CONSTANT field must interpolate to exactly that constant at every covered
# pixel: barycentric weights sum to one. If the rasterizer accumulates its edge
# functions incrementally but divides by a statically computed area, the weights
# drift and min/max bracket the constant instead of equalling it.
printf '10000000\n10000000\n10000000\n10000000\n' > "$G/const1e7.txt"
write_square "$G/sq.obj"
"$BIN" "$G/sq.obj" -v 0 0 -1 -r 0.0002 -d "$G/const1e7.txt" -p float -j 2>/dev/null > "$G/drift.json"
rel_near "$(jget "$G/drift.json" min)"     10000000 1e-6 "float: constant field min == the constant"
rel_near "$(jget "$G/drift.json" max)"     10000000 1e-6 "float: constant field max == the constant"
rel_near "$(jget "$G/drift.json" average)" 10000000 1e-6 "float: constant field average == the constant"

# Projected area must not depend on where the mesh sits in world space. The two
# meshes below are the same unit square, one at the origin and one translated by
# 5e6 — a perfectly ordinary magnitude for CAD or geospatial coordinates.
cat > "$G/near.obj" <<'OBJ'
v 0.0 0.0 0
v 0.85081111 0.525471651 0
v 0.325339458 1.376282761 0
v -0.525471651 0.85081111 0
f 1 2 3
f 1 3 4
OBJ
cat > "$G/far.obj" <<'OBJ'
v 5000000.0 5000000.0 0
v 5000000.85081111 5000000.525471651 0
v 5000000.325339458 5000001.376282761 0
v 4999999.474528349 5000000.85081111 0
f 1 2 3
f 1 3 4
OBJ
"$BIN" "$G/near.obj" -v 0 0 -1 -r 0.001 -p float -j 2>/dev/null > "$G/near.json"
"$BIN" "$G/far.obj"  -v 0 0 -1 -r 0.001 -p float -j 2>/dev/null > "$G/far.json"
rel_near "$(jget "$G/near.json" area)" 1.0 0.01 "float: unit square at the origin has area 1"
rel_near "$(jget "$G/far.json" area)"  1.0 0.01 "float: the SAME square at offset 5e6 has area 1"

# A view direction is a direction: its magnitude must be irrelevant. Computing
# the norm in the working type lets an extreme magnitude overflow or underflow
# to zero, silently yielding an empty projection.
for mag in 1e25 1e-25; do
    "$BIN" "$G/sq.obj" -v 0 0 -$mag -r 0.01 -p float -j 2>/dev/null > "$G/mag.json"
    rel_near "$(jget "$G/mag.json" area)" 1.0 0.01 "float: view magnitude $mag is irrelevant"
done

sect "Coverage is crack-free along shared edges"

# The reference cube is subdivided (26 vertices, 48 triangles), so each face is
# four quads split by diagonals -- and at these resolutions a whole run of pixel
# centres lands EXACTLY on those shared diagonals.
#
# The fill rule admits a pixel when all three edge functions are >= 0, so a
# centre lying exactly on a shared edge is claimed by BOTH neighbouring
# triangles. That only holds while edge(a,b,p) == -edge(b,a,p) is exact. Let the
# compiler fuse a multiply and an add into an FMA and the two sides round
# differently: both come out slightly negative, the pixel is claimed by NEITHER,
# and one-pixel cracks open along every diagonal. Hence -ffp-contract=off in the
# Makefile's FPFLAGS and in CMakeLists.txt.
#
# Seen head-on, a unit face at resolution 1/n must cover exactly n^2 pixels.
# Without that flag this reports 98 of 100 at 0.1 and 398 of 400 at 0.05 on any
# target that has an FMA instruction -- arm64, or x86-64 built with -mfma --
# while staying correct at 0.25, 0.125 and 0.025. The sweep is the test; no
# single resolution would have caught it.
for r in 0.5 0.25 0.2 0.125 0.1 0.05 0.025; do
    n=$(awk "BEGIN{printf \"%d\", 1/$r + 0.5}")
    "$BIN" "$CUBE" -v 1 0 0 --res "$r" -j 2>/dev/null > "$G/crack.json"
    equal "$(jget "$G/crack.json" pixels)" "$((n*n))" \
          "unit face at resolution $r covers exactly $((n*n)) pixels"
done

sect "Degenerate and empty views"

# The degenerate-triangle cutoff must scale with the geometry. An absolute
# threshold discards a legitimately small mesh outright.
printf 'v 0 0 0\nv 1e-7 0 0\nv 1e-7 1e-7 0\nv 0 1e-7 0\nf 1 2 3\nf 1 3 4\n' > "$G/tiny.obj"
"$BIN" "$G/tiny.obj" -v 0 0 -1 -r 1e-10 -p double -j 2>/dev/null > "$G/tiny.json"
rel_near "$(jget "$G/tiny.json" area)" 1e-14 0.05 "a 1e-7-scale mesh is not discarded as degenerate"

# ...and must still reject a genuinely collinear triangle, including in float
# where roundoff at ordinary coordinate magnitudes dwarfs any fixed epsilon.
printf 'v 1000.3137 1000.7193 0\nv 1000.9317339887 1001.1012660113 0\nv 1001.5497679774 1001.4832320226 0\nf 1 2 3\n' > "$G/collinear.obj"
"$BIN" "$G/collinear.obj" -v 0 0 -1 -r 0.001 -p float --no-cull -j 2>/dev/null > "$G/col.json"
near "$(jget "$G/col.json" area)" 0.0 1e-9 "float: a collinear triangle contributes no area"

# When nothing is covered the field statistics are not measurements. Reporting a
# fabricated 0 corrupts any batch-wide min/max; the honest answer is null.
printf 'v 0 0 0\nv 1 0 0\nv 1 0.000001 0\nf 1 2 3\n' > "$G/sliver.obj"
printf -- '-5\n-3\n-1\n' > "$G/neg.txt"
"$BIN" "$G/sliver.obj" -v 0 0 -1 -r 0.1 -d "$G/neg.txt" -p double --no-cull -j 2>/dev/null > "$G/sliver.json"
equal "$(jget "$G/sliver.json" pixels)" "0" "sliver view covers no pixels"
equal "$(jraw "$G/sliver.json" average)"  "null" "zero-coverage average is null, not 0"
equal "$(jraw "$G/sliver.json" integral)" "null" "zero-coverage integral is null, not 0"
equal "$(jraw "$G/sliver.json" min)"      "null" "zero-coverage min is null, not 0"
equal "$(jraw "$G/sliver.json" max)"      "null" "zero-coverage max is null, not 0"

# A view that rasterizes nothing must not inherit the previous view's raster.
# The renderer is reused across views by each worker, so skipping the resize on
# an empty view leaves stale dimensions pointing at a stale (or undersized)
# value buffer — a wrong image at best, an out-of-bounds read at worst.
printf '1\n2\n3\n4\n' > "$G/f4.txt"
rm -f "$G/stale"_*.png
printf '0 0 -1 0.01\n1 0 0 0.01\n' > "$G/stale.txt"      # view 1 is exactly edge-on
"$BIN" "$G/sq.obj" -b "$G/stale.txt" --no-cull -t 1 -o "$G/stale" -j 2>/dev/null > "$G/stale.json"
expect_ok $? "empty view after a rendered one does not crash"
equal "$(jget_n "$G/stale.json" pixels 2)" "0" "edge-on view covers no pixels"
equal "$(jraw_n "$G/stale.json" image 2)" "null" "no image file is claimed for an empty view"
[ ! -f "$G/stale_0001.png" ] && ok "no stale image written for an empty view" \
                             || bad "stale image written for an empty view"

# The same hazard with the value buffer: an empty FIRST view, then a view with a
# field, inside one worker. This is the configuration that segfaulted.
rm -f "$G/crash"_*.png
printf '0 0 1 0.01\n0 0 -1 0.01 f4.txt\n' > "$G/crash.txt"
"$BIN" "$G/sq.obj" -b "$G/crash.txt" -t 1 -o "$G/crash" >/dev/null 2>&1
expect_ok $? "empty view followed by a field view does not crash"
