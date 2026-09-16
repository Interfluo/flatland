#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
# for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

"""
Example 4 — A whole time series in one call.

The question this answers: *I have hundreds of timesteps, each with its own view
direction and its own field. How do I not wait all afternoon?*

`project_batch` takes every view and the entire field matrix across the boundary
once, then runs the views across all cores. The mesh is validated and recentred a
single time no matter how many timesteps there are.

    python3 examples/python/ex04_timeseries.py

PASS NUMPY ARRAYS IF YOU HAVE THEM. This is the one place it matters for speed.
ctypes has to walk a Python list element by element to build the C array, and on
a 36k-vertex mesh that costs more than the projection itself — enough to hide
the parallelism completely. A NumPy array is handed over as a buffer, with no
per-element work. Measured on a 4-core machine, 69k triangles, r = 1e-3:

    plain lists   6.9 ms/step one at a time; batching is not a win at all,
                  because marshalling dominates and it happens either way
    NumPy         3.2 ms/step one at a time, 0.78 ms/step batched  ->  4.1x

Both give identical answers; only one gives you the parallelism.

Covered:
  - one field matrix, one column per timestep
  - what the batch path actually buys, measured rather than asserted
  - per-view resolutions
  - guarding on has_stats before aggregating
"""

import math
import time

from _setup import load, mesh_path, rule

fl = load()

# The list path costs ~3 ms/step in ctypes marshalling, so keep the demo
# short when NumPy is absent.
STEPS = 240 if fl.HAS_NUMPY else 48


def build_field_matrix(nx, steps):
    """A travelling wave along x: one row per vertex, one column per timestep.

    Substitute your own solver output. The only requirement is the shape —
    (n_entities, n_timesteps) — and that row r is the value at mesh entity r.
    """
    if fl.HAS_NUMPY:
        import numpy as np
        t = np.arange(steps)
        return np.sin(4.0 * np.pi * nx[:, None] - 2.0 * np.pi * t[None, :] / steps)
    return [[math.sin(4.0 * math.pi * x - 2.0 * math.pi * t / steps)
             for t in range(steps)] for x in nx]


def main():
    print("FlatLand %s   (NumPy: %s)" % (fl.__version__, fl.HAS_NUMPY))
    if not fl.HAS_NUMPY:
        print("  NumPy is not installed, so this runs the slower list path and")
        print("  uses %d steps instead of 240. See the module docstring." % STEPS)

    mesh = fl.Mesh.load(mesh_path("with_fields", "bunny.obj"))
    nv = mesh.vertex_count
    print("loaded the bunny: %d vertices, %d triangles" % (nv, mesh.face_count))

    # ---------------------------------------------------------------------
    rule("Building the time series")

    views = [fl.angle_to_dir(360.0 * t / STEPS, 20.0) for t in range(STEPS)]

    xs = [v[0] for v in mesh.vertices]
    x0, x1 = min(xs), max(xs)
    span = (x1 - x0) or 1.0
    if fl.HAS_NUMPY:
        import numpy as np
        views = np.asarray(views)
        nx = (np.asarray(xs) - x0) / span
    else:
        nx = [(x - x0) / span for x in xs]

    field_matrix = build_field_matrix(nx, STEPS)
    columns = list(range(STEPS))
    print("  %d views, field matrix %d x %d (%d values)"
          % (STEPS, nv, STEPS, nv * STEPS))

    # ---------------------------------------------------------------------
    rule("Running it")

    t0 = time.perf_counter()
    results = mesh.project_batch(
        views,
        field_matrix=field_matrix,
        field_columns=columns,      # timestep t reads column t
        resolution=1e-3,
        threads=0,                  # 0 = one worker per core
    )
    batch_time = time.perf_counter() - t0
    print("  %d timesteps in %.3f s  (%.0f steps/s, %.2f ms/step)"
          % (len(results), batch_time, len(results) / batch_time,
             batch_time * 1e3 / len(results)))

    # ---------------------------------------------------------------------
    rule("What the batch call bought")

    # Time the same views one at a time, with the field columns extracted OUTSIDE
    # the timed region so this measures FlatLand and not list slicing.
    n_cmp = min(24, STEPS)
    if fl.HAS_NUMPY:
        import numpy as np
        cols = [np.ascontiguousarray(field_matrix[:, t]) for t in range(n_cmp)]
    else:
        cols = [[field_matrix[i][t] for i in range(nv)] for t in range(n_cmp)]

    t0 = time.perf_counter()
    for t in range(n_cmp):
        mesh.project(views[t], field=cols[t], resolution=1e-3)
    loop_time = time.perf_counter() - t0

    # Single-threaded batch, to separate "crossing the boundary once" from
    # "running on every core".
    t0 = time.perf_counter()
    mesh.project_batch(views[:n_cmp], field_matrix=field_matrix,
                       field_columns=columns[:n_cmp], resolution=1e-3, threads=1)
    serial_batch = time.perf_counter() - t0

    print("  one at a time        %6.2f ms/step" % (loop_time * 1e3 / n_cmp))
    print("  batched, 1 thread    %6.2f ms/step" % (serial_batch * 1e3 / n_cmp))
    print("  batched, all cores   %6.2f ms/step" % (batch_time * 1e3 / STEPS))
    print("  speedup              %6.2fx" % (loop_time / n_cmp / (batch_time / STEPS)))
    if not fl.HAS_NUMPY:
        print("  (mostly ctypes list conversion here — install NumPy to see the")
        print("   parallelism rather than the marshalling)")

    # ---------------------------------------------------------------------
    rule("The results")

    # Always check has_stats before aggregating. A view that covered no pixels
    # has no field statistics, and Python reports them as None rather than 0.0
    # precisely so they cannot be averaged in by accident.
    good = [r for r in results if r.has_stats]
    if len(good) != len(results):
        print("  %d of %d views covered no pixels and were excluded"
              % (len(results) - len(good), len(results)))

    areas = [r.area for r in results]
    integrals = [r.integral for r in good]
    print("  area      min %.6f  max %.6f  mean %.6f"
          % (min(areas), max(areas), sum(areas) / len(areas)))
    print("  integral  min %+.6f  max %+.6f  mean %+.6f"
          % (min(integrals), max(integrals), sum(integrals) / len(integrals)))

    peak = max(range(len(results)), key=lambda i: results[i].area)
    print("  largest silhouette at timestep %d (azimuth %.0f deg): %.6f"
          % (peak, 360.0 * peak / STEPS, results[peak].area))

    print("\n  first six timesteps:")
    print("  %-5s %-10s %-12s %-12s %-8s" % ("step", "area", "mean f", "integral", "pixels"))
    for t in range(6):
        r = results[t]
        print("  %-5d %-10.6f %-+12.6f %-+12.6f %-8d"
              % (t, r.area, r.average, r.integral, r.covered_pixels))

    # ---------------------------------------------------------------------
    rule("Per-view resolution")

    # Each view can carry its own pixel size — useful for a cheap coarse pass
    # before a fine one, or when some orientations need more detail.
    res_per_view = [4e-3] * 4 + [1e-3] * 4
    sub = mesh.project_batch(views[:8], resolutions=res_per_view)
    for t, r in enumerate(sub):
        print("  step %d  res %.4f  area %.6f  raster %dx%d"
              % (t, res_per_view[t], r.area, r.width, r.height))

    mesh.close()


if __name__ == "__main__":
    main()
