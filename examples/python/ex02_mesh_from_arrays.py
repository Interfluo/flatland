#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
# for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

"""
Example 2 — Geometry straight from memory, no files involved.

The question this answers: *my mesh is already in Python — do I have to write it
to disk?* No. `Mesh()` takes vertex and face arrays directly, which is the whole
reason the C ABI exists. Nothing here touches the filesystem.

    python3 examples/python/ex02_mesh_from_arrays.py

Covered:
  - building a mesh from plain sequences, and from NumPy arrays
  - face winding, and why it decides which surface you measure
  - sweeping a parametric shape without ever serialising it
  - reading the geometry back out
"""

import math

from _setup import load, rule

fl = load()


def tessellated_cylinder(sides, radius=1.0, height=2.0):
    """A closed prism approximating a cylinder, built as plain lists.

    Winding matters: each triangle's vertices must run counter-clockwise seen
    from OUTSIDE, so its normal points outward. FlatLand culls faces whose
    normal points away from the camera, so a mesh wound inside-out measures the
    far surface of the object instead of the near one. The areas would look
    fine; any field statistic would be wrong.
    """
    hz = height / 2.0
    verts = []
    for z in (-hz, hz):
        for i in range(sides):
            a = 2.0 * math.pi * i / sides
            verts.append((radius * math.cos(a), radius * math.sin(a), z))
    bottom_centre = len(verts); verts.append((0.0, 0.0, -hz))
    top_centre = len(verts); verts.append((0.0, 0.0, hz))

    faces = []
    for i in range(sides):
        j = (i + 1) % sides
        lo_i, lo_j = i, j
        hi_i, hi_j = sides + i, sides + j
        faces.append((lo_i, lo_j, hi_j))           # wall
        faces.append((lo_i, hi_j, hi_i))
        faces.append((bottom_centre, lo_j, lo_i))  # bottom cap, normal -z
        faces.append((top_centre, hi_i, hi_j))     # top cap, normal +z
    return verts, faces


def main():
    print("FlatLand %s   (NumPy available: %s)" % (fl.__version__, fl.HAS_NUMPY))

    # ---------------------------------------------------------------------
    rule("A mesh from plain Python sequences")

    verts, faces = tessellated_cylinder(sides=128, radius=1.0, height=2.0)
    mesh = fl.Mesh(verts, faces)
    print("  built a %d-sided prism in memory: %d vertices, %d triangles"
          % (128, mesh.vertex_count, mesh.face_count))

    # Viewed down the axis this is a disc of area pi r^2; side-on it is a
    # 2r x h rectangle. Both have closed forms, so the numbers are checkable.
    end_on = mesh.project((0, 0, 1), resolution=1e-3, precision="double")
    side_on = mesh.project((1, 0, 0), resolution=1e-3, precision="double")
    print("  end-on  %.6f   (pi r^2      = %.6f)" % (end_on.area, math.pi))
    print("  side-on %.6f   (2 r h       = %.6f)" % (side_on.area, 4.0))

    # ---------------------------------------------------------------------
    rule("The same thing with NumPy")

    if fl.HAS_NUMPY:
        import numpy as np
        # NumPy arrays go in exactly the same places. Shape (N,3) float for
        # vertices, (M,3) integer for faces; a flat array of length 3N/3M also
        # works if that is what your pipeline produces.
        v = np.asarray(verts, dtype=np.float64)
        f = np.asarray(faces, dtype=np.int32)
        with fl.Mesh(v, f) as m:
            r = m.project(np.array([0.0, 0.0, 1.0]), resolution=1e-3,
                          precision="double")
            print("  ndarray input gives the same answer: %.6f" % r.area)
            print("  (agrees with the list version: %s)"
                  % (abs(r.area - end_on.area) < 1e-12))
    else:
        print("  NumPy is not installed here, so this section is skipped.")
        print("  The binding needs it only if you want to pass arrays; plain")
        print("  sequences work either way.")

    # ---------------------------------------------------------------------
    rule("Refining a parametric shape without writing a single file")

    # The prism converges on a true cylinder as the facet count grows. Watching
    # that convergence is the sort of loop that is painful through a CLI and
    # trivial in process: the mesh never leaves memory.
    exact = math.pi          # end-on area of the smooth cylinder, r = 1
    print("  %-8s %-12s %-12s" % ("sides", "end-on area", "rel err vs pi"))
    for sides in (8, 16, 32, 64, 128, 256):
        v, f = tessellated_cylinder(sides)
        with fl.Mesh(v, f) as m:
            a = m.project((0, 0, 1), resolution=1e-3, precision="double").area
        print("  %-8d %-12.7f %.2e" % (sides, a, abs(a - exact) / exact))

    # ---------------------------------------------------------------------
    rule("Reading the geometry back")

    # vertices/faces return what you put in — FlatLand recentres meshes
    # internally for numerical conditioning, but that is invisible here.
    back_v = mesh.vertices
    back_f = mesh.faces
    print("  vertices out: %d, faces out: %d" % (len(back_v), len(back_f)))
    print("  first vertex in : %s" % (tuple(round(c, 6) for c in verts[0]),))
    print("  first vertex out: %s" % (tuple(round(c, 6) for c in back_v[0]),))
    print("  round trip exact: %s"
          % all(abs(a - b) < 1e-12 for a, b in zip(verts[0], back_v[0])))

    mesh.close()


if __name__ == "__main__":
    main()
