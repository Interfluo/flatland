#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
# for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

"""
example.py — an end-to-end tour of the FlatLand Python binding.

Run it from anywhere:

    make lib                                  # build libflatland.so first
    PYTHONPATH=python python3 python/example.py

It works with or without NumPy; pass a repository path as the first argument if
you are running a copy from somewhere else.

What it covers, in order:
  1. a mesh built from plain Python sequences, and its projected area
  2. a varying node field, and why the average proves which surface was seen
  3. a face field, and the area integral
  4. a whole time series in one parallel batch call
  5. statistics that do not exist, and are reported as None rather than 0.0
  6. a raster, its coverage mask, a PNG, and the raw values as .npy
  7. a mesh loaded from an OBJ file
"""

import os
import sys

import flatland


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))
    )

    print("FlatLand %s  (%s)" % (flatland.__version__, flatland.library_path()))
    print("NumPy in use: %s\n" % flatland.HAS_NUMPY)

    # ---------------------------------------------------------------- 1 -----
    # A closed unit box spanning [0,1]^3, CCW outward normals. Lists of tuples
    # here; NumPy arrays are accepted in exactly the same places.
    vertices = [
        (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
        (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
    ]
    faces = [
        (0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
    ]

    with flatland.Mesh(vertices, faces) as box:
        print("1. %s" % box)
        r = box.project((1, 0, 0), resolution=0.002)
        print("   looking along +X: area = %.6f over %d pixels (%dx%d raster)\n"
              % (r.area, r.covered_pixels, r.width, r.height))

        # ------------------------------------------------------------ 2 -----
        # A node field equal to each vertex's x coordinate. The camera looks
        # ALONG +X, so the visible face is the one at x = 0 and the mean must be
        # 0. Reverse the view and it is 1. That difference is the whole reason
        # to test with a field that varies.
        field_x = [v[0] for v in vertices]
        near = box.project((1, 0, 0), field=field_x, resolution=0.002)
        far = box.project((-1, 0, 0), field=field_x, resolution=0.002)
        print("2. node field = vertex x")
        print("   view (+1, 0, 0): average = %.6f   <- the near surface, x = 0"
              % near.average)
        print("   view (-1, 0, 0): average = %.6f   <- the opposite face, x = 1\n"
              % far.average)

        # ------------------------------------------------------------ 3 -----
        # A face field is constant over each triangle, so the area integral of a
        # constant field is just that constant times the area. Supply pressures
        # and the integral is a force; supply radiance and it is an intensity.
        pressure = [2.5] * box.face_count
        r = box.project((0, 0, 1), field=pressure, resolution=0.002)
        print("3. face field = 2.5 everywhere")
        print("   average  = %.6f" % r.average)
        print("   integral = %.6f   (= 2.5 x area %.6f)\n" % (r.integral, r.area))

        # ------------------------------------------------------------ 4 -----
        # One row per vertex, one column per timestep: a whole time series in a
        # single array, projected in parallel. Each view reads its own column.
        views = [
            flatland.angle_to_dir(0, 0),      # +X
            flatland.angle_to_dir(90, 0),     # +Y
            flatland.angle_to_dir(0, 90),     # +Z
        ]
        timeseries = [[10.0, 20.0, 30.0] for _ in vertices]  # 8 rows x 3 steps
        results = box.project_batch(
            views,
            field_matrix=timeseries,
            field_columns=[0, 1, 2],
            resolution=0.004,
            threads=0,  # one worker per core
        )
        print("4. batch of %d views over a %dx%d field matrix"
              % (len(views), len(timeseries), len(timeseries[0])))
        for i, res in enumerate(results):
            print("   view %d: area = %.4f  average = %.1f  integral = %.4f"
                  % (i, res.area, res.average, res.integral))
        print()

        # ------------------------------------------------------------ 5 -----
        # A view can legitimately measure nothing. Statistics that were never
        # measured come back as None, not as a plausible-looking 0.0.
        sliver = flatland.Mesh([(0, 0, 0), (1, 0, 0), (1, 1e-6, 0)], [(0, 1, 2)])
        z = sliver.project((0, 0, -1), field=[-5.0, -3.0, -1.0],
                           resolution=0.1, cull=False)
        print("5. a sliver thinner than one pixel")
        print("   covered_pixels = %d, has_field = %s, has_stats = %s"
              % (z.covered_pixels, z.has_field, z.has_stats))
        print("   average = %r, integral = %r, min = %r, max = %r"
              % (z.average, z.integral, z.min, z.max))
        print("   (None, not 0.0: nothing was measured, so there is nothing to "
              "report)\n")
        sliver.close()

        # ------------------------------------------------------------ 6 -----
        # render() keeps the raster. mask and values are flat, width*height
        # long, row 0 at the bottom, pixel [y * width + x].
        image = box.render((1, 1, 0.4), field=field_x, resolution=0.01)
        covered = sum(int(v) for v in image.mask)
        print("6. raster: %dx%d, %d covered pixels (result says %d)"
              % (image.width, image.height, covered,
                 image.result.covered_pixels))
        out = os.path.join(os.getcwd(), "example_render.png")
        image.save_png(out)
        print("   wrote %s (%d bytes, false-colour PNG)\n"
              % (out, os.path.getsize(out)))
        image.close()

    # -------------------------------------------------------------- 7 -------
    cube = os.path.join(root, "examples", "cube_area", "cube.obj")
    if os.path.exists(cube):
        with flatland.Mesh.load(cube) as mesh:
            r = mesh.project((1, 0, 0), resolution=0.005)
            print("7. %s" % mesh)
            print("   projected area = %.6f (a unit cube, so 1)" % r.area)
    else:
        print("7. (skipped: %s not found)" % cube)

    return 0


if __name__ == "__main__":
    sys.exit(main())
