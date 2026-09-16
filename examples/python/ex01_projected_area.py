#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
# for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

"""
Example 1 — Projected area, and finding the orientation that minimises it.

The question this answers: *how much of my part does a given direction see, and
which orientation exposes the least (or most) of it?* That shows up as frontal
area for drag, as an aperture for radiative load, as a footprint for packing.

    python3 examples/python/ex01_projected_area.py

Covered:
  - loading an OBJ and projecting along a few directions
  - an azimuth/elevation sweep in a single parallel batch call
  - checking the answer against a closed form, because for a cube there is one
  - choosing a resolution by watching the answer converge
"""

import math

from _setup import load, mesh_path, rule

fl = load()


def main():
    print("FlatLand %s" % fl.__version__)

    mesh = fl.Mesh.load(mesh_path("cube_area", "cube.obj"))
    print("loaded a mesh with %d vertices and %d triangles"
          % (mesh.vertex_count, mesh.face_count))

    # ---------------------------------------------------------------------
    rule("A few directions")

    # A view direction is the direction the camera LOOKS ALONG. Magnitude is
    # irrelevant, so (2,0,0) and (1,0,0) are the same view.
    for view in [(1, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)]:
        r = mesh.project(view, resolution=1e-3)
        print("  view %-12s area %.6f   (%d pixels, %dx%d raster)"
              % (str(view), r.area, r.covered_pixels, r.width, r.height))

    # ---------------------------------------------------------------------
    rule("Is that right?")

    # This mesh is a unit cube, and a convex polyhedron's projected area has a
    # closed form:  A(n) = 1/2 sum_i A_i |n . n_i|.  For a unit cube that
    # reduces to |nx| + |ny| + |nz|. Worth checking against, because a silently
    # wrong area looks exactly like a right one.
    def cube_exact(view):
        n = math.sqrt(sum(c * c for c in view))
        return sum(abs(c) / n for c in view)

    worst = 0.0
    for view in [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, -1, 3), (0.3, 0.9, -0.31)]:
        got = mesh.project(view, resolution=5e-4, precision="double").area
        want = cube_exact(view)
        rel = abs(got - want) / want
        worst = max(worst, rel)
        print("  %-20s measured %.7f   exact %.7f   rel err %.1e"
              % (str(view), got, want, rel))
    print("  worst relative error: %.1e" % worst)

    # ---------------------------------------------------------------------
    rule("Sweeping orientation to find the extremes")

    # 60 directions around the object, evaluated in ONE call. project_batch
    # parses and centres the mesh once and runs the views across all cores;
    # looping over project() would redo that work every time.
    views, labels = [], []
    for az in range(0, 360, 30):
        for el in (-60, -30, 0, 30, 60):
            views.append(fl.angle_to_dir(az, el))
            labels.append((az, el))

    results = mesh.project_batch(views, resolution=2e-3, threads=0)
    areas = [r.area for r in results]

    lo = min(range(len(areas)), key=lambda i: areas[i])
    hi = max(range(len(areas)), key=lambda i: areas[i])
    print("  %d orientations evaluated in one batch" % len(views))
    print("  minimum  %.5f at azimuth %3d, elevation %3d"
          % (areas[lo], labels[lo][0], labels[lo][1]))
    print("  maximum  %.5f at azimuth %3d, elevation %3d"
          % (areas[hi], labels[hi][0], labels[hi][1]))
    print("  mean     %.5f over all orientations" % (sum(areas) / len(areas)))

    # For a convex body the mean projected area over ALL directions is exactly
    # a quarter of the surface area (Cauchy). A cube of side 1 has S = 6, so
    # the true mean is 1.5 — this grid is coarse and biased toward the equator,
    # so it only lands near it, but a full sweep converges.
    print("  (Cauchy's formula gives S/4 = 1.5 for the true all-direction mean)")

    # ---------------------------------------------------------------------
    rule("Choosing a resolution")

    # Coverage is decided by testing whether a pixel CENTRE falls inside a
    # triangle, so the answer converges as the pixel shrinks. Sweep it and stop
    # when the digits you care about stop moving. Convergence is first order
    # and not monotonic, so look at the trend rather than one step.
    exact = cube_exact((0.3, 0.9, -0.31))
    print("  %-10s %-12s %-10s" % ("pixel", "area", "rel err"))
    for res in (2e-2, 1e-2, 5e-3, 2e-3, 1e-3, 5e-4):
        a = mesh.project((0.3, 0.9, -0.31), resolution=res, precision="double").area
        print("  %-10.5f %-12.7f %.2e" % (res, a, abs(a - exact) / exact))

    mesh.close()


if __name__ == "__main__":
    main()
