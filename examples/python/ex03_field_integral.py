#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
# for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

"""
Example 3 — The area integral of a scalar field over the visible projection.

The question this answers: *I have a quantity defined on my mesh. What is its
integral over the part a given direction can actually see?*

    I = integral f dA = sum over covered pixels of (value x pixel_area)

FlatLand attaches no physical meaning to f. Supply a radiance and I is a radiant
intensity; supply a pressure and it is a force; supply an emissivity-weighted
temperature and it is a thermal signature. The arithmetic is the same and the
interpretation is yours.

    python3 examples/python/ex03_field_integral.py

Covered:
  - node fields (interpolated) and face fields (piecewise constant)
  - the Lambertian sphere, checked against its closed form
  - which surface the statistics come from, and why that matters
  - has_stats: statistics that do not exist, reported as None rather than 0.0
"""

import math

from _setup import load, rule

fl = load()


def icosphere(subdivisions=3, radius=1.0):
    """An icosahedron subdivided and pushed onto a sphere."""
    t = (1.0 + math.sqrt(5.0)) / 2.0
    raw = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0),
           (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
           (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    faces = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
             (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
             (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
             (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]

    def unit(p):
        n = math.sqrt(sum(c * c for c in p))
        return (p[0] / n, p[1] / n, p[2] / n)

    verts = [list(unit(p)) for p in raw]
    for _ in range(subdivisions):
        cache, out = {}, []

        def mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                m = unit(tuple(verts[a][i] + verts[b][i] for i in range(3)))
                verts.append(list(m))
                cache[key] = len(verts) - 1
            return cache[key]

        for (a, b, c) in faces:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            out += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        faces = out

    return [tuple(c * radius for c in p) for p in verts], faces


def main():
    print("FlatLand %s" % fl.__version__)

    # ---------------------------------------------------------------------
    rule("A node field: the Lambertian sphere")

    # f = cos(angle between the surface normal and the view direction). For a
    # sphere centred at the origin the outward normal at a vertex is just the
    # vertex direction, so f = -vhat . n on the visible side.
    #
    # This has a closed form. At projected radius rho the cosine is
    # sqrt(1 - (rho/r)^2), so
    #
    #   integral f dA = int_0^r sqrt(1-(rho/r)^2) 2 pi rho drho = (2/3) pi r^2
    #   mean          = 2/3
    #
    # If you read f as a radiance, that integral is the radiant intensity of a
    # Lambertian sphere. FlatLand does not know or care.
    view = (0.0, 0.0, 1.0)
    exact_integral = (2.0 / 3.0) * math.pi
    exact_mean = 2.0 / 3.0

    print("  %-8s %-8s %-12s %-12s %-10s" % ("subdiv", "tris", "integral", "exact", "rel err"))
    for k in (1, 2, 3, 4):
        verts, faces = icosphere(k, 1.0)
        # One field value per vertex, in mesh vertex order.
        field = [-(v[0] * view[0] + v[1] * view[1] + v[2] * view[2]) for v in verts]
        with fl.Mesh(verts, faces) as mesh:
            r = mesh.project(view, field=field, resolution=2e-3, precision="double")
        rel = abs(r.integral - exact_integral) / exact_integral
        print("  %-8d %-8d %-12.7f %-12.7f %.2e" % (k, len(faces), r.integral, exact_integral, rel))
    print("  mean over the disc: %.6f   (exact 2/3 = %.6f)" % (r.average, exact_mean))
    print("  min %.4f at the limb, max %.4f at the sub-observer point" % (r.min, r.max))

    # ---------------------------------------------------------------------
    rule("A face field: constant per triangle")

    # A face field has one value per TRIANGLE and is constant across it — right
    # for something measured per facet, like a per-panel emissivity or a
    # per-element CFD cell value. The length decides which kind it is; when a
    # mesh has as many faces as vertices, pass field_mode explicitly.
    verts, faces = icosphere(2, 1.0)
    face_field = [2.5] * len(faces)
    with fl.Mesh(verts, faces) as mesh:
        r = mesh.project(view, field=face_field, resolution=2e-3, precision="double")
        print("  %d faces, all carrying 2.5" % len(faces))
        print("  mean     %.6f  (a constant field must come back exactly)" % r.average)
        print("  integral %.6f  = 2.5 x area %.6f = %.6f"
              % (r.integral, r.area, 2.5 * r.area))

    # ---------------------------------------------------------------------
    rule("Which surface do the statistics describe?")

    # The visible one — the surface facing the camera. A field that varies
    # through the object makes that concrete. On a unit box with f = x, looking
    # along +x means the near face is x = 0, so the mean must be 0.
    box_v = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
             (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    box_f = [(0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7),
             (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
             (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7)]
    fx = [v[0] for v in box_v]

    with fl.Mesh(box_v, box_f) as box:
        near = box.project((1, 0, 0), field=fx, resolution=2e-3, precision="double")
        far = box.project((-1, 0, 0), field=fx, resolution=2e-3, precision="double")
        print("  looking along +x: mean f = %.4f   (the x = 0 face)" % near.average)
        print("  looking along -x: mean f = %.4f   (the x = 1 face)" % far.average)
        print("  the areas are identical (%.4f), only the field differs" % near.area)

    # ---------------------------------------------------------------------
    rule("When there are no statistics to report")

    # A view can legitimately cover nothing: geometry thinner than a pixel, or
    # a mesh seen exactly edge-on. The four field statistics are then ABSENT,
    # not zero — Python gives you None. Zero is an ordinary field value, so
    # passing it through would be indistinguishable from a measurement and
    # would quietly drag any mean or min you compute toward it.
    sliver_v = [(0, 0, 0), (1, 0, 0), (1, 1e-6, 0)]
    sliver_f = [(0, 1, 2)]
    with fl.Mesh(sliver_v, sliver_f) as sliver:
        r = sliver.project((0, 0, -1), field=[-5.0, -3.0, -1.0],
                           resolution=0.1, cull=False)
        print("  covered pixels : %d" % r.covered_pixels)
        print("  has_field      : %s   (a field WAS supplied)" % r.has_field)
        print("  has_stats      : %s  (but nothing was measured)" % r.has_stats)
        print("  average        : %s" % r.average)
        print("  area           : %s   (still a real measurement, so not None)" % r.area)

    print("\n  When aggregating a batch, guard on has_stats:")
    print("    good = [r for r in results if r.has_stats]")
    print("    mean = sum(r.average for r in good) / len(good)")


if __name__ == "__main__":
    main()
