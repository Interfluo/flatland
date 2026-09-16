#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
# for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

"""Check the analytic formulas in cases.py against independently-known values.

This runs FIRST, before any FlatLand comparison, and it exists because of an
obvious failure mode: if the closed forms were wrong, every study built on them
would agree with the wrong answer and report success. So each formula is pinned
against a value known from elsewhere — a textbook constant, a different formula
for the same quantity, or a limit the shape must converge to.

Exits non-zero if any identity fails. No FlatLand binary is needed.
"""

import math
import sys

import cases as C

fails = 0


def chk(got, want, tol, what):
    global fails
    ok = abs(got - want) <= tol
    if not ok:
        fails += 1
    print("  %s %-54s got %.10f  want %.10f"
          % ("ok " if ok else "BAD", what, got, want))


def cauchy_mean(verts, faces, n_dirs):
    dirs = C.fibonacci_directions(n_dirs)
    return sum(C.projected_area_convex(verts, faces, d) for d in dirs)/len(dirs)


def main():
    # Cauchy's identity is exact; only the direction SAMPLING is approximate, and
    # its error falls off quickly with a Fibonacci spiral. The sample count is
    # therefore scaled to the mesh size — a 12-triangle cube can afford 20000
    # directions, a 5120-triangle sphere cannot, and neither needs more than
    # enough to sit inside the tolerance.
    print("cube")
    v, f = C.unit_cube()
    chk(C.surface_area(v, f), 6.0, 1e-12, "surface area == 6")
    chk(C.projected_area_convex(v, f, (1, 0, 0)), 1.0, 1e-12, "A(+x) == 1")
    # Projecting a cube down its space diagonal gives a regular hexagon of area sqrt(3).
    chk(C.projected_area_convex(v, f, (1, 1, 1)), math.sqrt(3), 1e-12,
        "A(body diagonal) == sqrt(3), the hexagon")
    chk(C.projected_area_convex(v, f, (1, 1, 0)), math.sqrt(2), 1e-12,
        "A(face diagonal) == sqrt(2)")
    # Two independent formulas for the same quantity must agree everywhere.
    worst = max(abs(C.projected_area_convex(v, f, d) - C.cube_projected_area(d))
                for d in C.fibonacci_directions(500))
    chk(worst, 0.0, 1e-12, "|nx|+|ny|+|nz| == the general convex formula")
    chk(cauchy_mean(v, f, 20000), 1.5, 2e-3, "Cauchy: <A> == S/4 == 1.5")

    print("tetrahedron")
    v, f = C.regular_tetrahedron(1.0)
    chk(C.norm(C.sub(v[0], v[1])), 1.0, 1e-12, "edge length == 1 as requested")
    chk(C.surface_area(v, f), math.sqrt(3), 1e-12, "S == sqrt(3) e^2")
    chk(cauchy_mean(v, f, 20000), math.sqrt(3)/4, 2e-3, "Cauchy: <A> == sqrt(3)/4")

    print("octahedron")
    v, f = C.regular_octahedron(1.0)
    edge = C.norm(C.sub(v[0], v[2]))
    chk(edge, math.sqrt(2), 1e-12, "edge == sqrt(2) at radius 1")
    chk(C.surface_area(v, f), 2*math.sqrt(3)*edge*edge, 1e-12, "S == 2 sqrt(3) a^2")
    chk(C.projected_area_convex(v, f, (0, 0, 1)), 2.0, 1e-12,
        "A(+z) == 2, the square cross-section")

    print("icosphere -> sphere")
    prev_S = prev_A = 0.0
    for k in (0, 1, 2, 3, 4):
        v, f = C.icosphere(k, 1.0)
        S = C.surface_area(v, f)
        A = C.projected_area_convex(v, f, (1, 0.3, 0.2))
        # An inscribed polyhedron must UNDERSTATE the sphere, and must improve
        # monotonically with subdivision.
        assert S < 4*math.pi and A < math.pi, "an inscribed polyhedron cannot exceed the sphere"
        assert S > prev_S and A > prev_A, "subdivision must improve the approximation"
        prev_S, prev_A = S, A
        # Cauchy is a statement about the MEAN over directions, not about any
        # single one: a 20-face icosahedron is anisotropic enough that one
        # direction misses S/4 by 2.5%. Averaged, it holds on the polyhedron
        # itself at every level of subdivision.
        if k in (0, 2):
            chk(cauchy_mean(v, f, 4000), S/4.0, 3e-3,
                "subdiv %d: Cauchy <A> == S/4 on the polyhedron" % k)
    v, f = C.icosphere(4, 1.0)
    chk(C.surface_area(v, f), 4*math.pi, 0.02, "S -> 4 pi (0.12%% low at subdiv 4)")
    chk(C.projected_area_convex(v, f, (1, 0.3, 0.2)), math.pi, 0.005, "A -> pi")
    spread = [C.projected_area_convex(v, f, d) for d in C.fibonacci_directions(120)]
    chk(max(spread) - min(spread), 0.0, 2e-3, "a sphere projects alike from every direction")

    print("cylinder -> smooth cylinder")
    for nsides, tol in ((64, 6e-3), (256, 4e-4)):
        v, f = C.cylinder(nsides, 1.0, 2.0)
        for deg in (0, 45, 90):
            th = math.radians(deg)
            got = C.projected_area_convex(v, f, (math.sin(th), 0.0, math.cos(th)))
            chk(got, C.cylinder_projected_area(1.0, 2.0, th), tol,
                "%d-gon prism at %d deg -> 2rh sin + pi r^2 cos" % (nsides, deg))

    print("planar field integrals")
    v, f = C.unit_square()
    I, A, m = C.field_integral_planar(v, f, C.coordinate_field(v, 0), (0, 0, 1))
    chk(A, 1.0, 1e-12, "unit square area == 1")
    chk(I, 0.5, 1e-12, "integral of x over the unit square == 1/2")
    chk(m, 0.5, 1e-12, "its mean == 1/2")
    v, f = C.right_triangle()
    I, A, m = C.field_integral_planar(v, f, C.coordinate_field(v, 0), (0, 0, 1))
    chk(A, 0.5, 1e-12, "right triangle area == 1/2")
    chk(I, 1.0/6, 1e-12, "integral of x over it == 1/6")
    chk(m, 1.0/3, 1e-12, "its mean == 1/3 == mean of the vertex values")

    print("Lambertian sphere")
    chk(C.lambert_integral(1.0), 2*math.pi/3, 1e-12, "closed form == 2 pi / 3")
    # Independent numerical confirmation: integrate over the visible hemisphere
    # of a fine mesh and check it lands on the analytic value.
    v, f = C.icosphere(4, 1.0)
    n = (0.0, 0.0, 1.0)
    vals = C.lambert_field(v, n, 1.0)
    tot_a = tot_i = 0.0
    for tri in f:
        c = C.face_cross(v, tri)
        if C.dot(c, n) >= 0:
            continue                       # back-facing to a camera looking along +z
        a = 0.5*abs(C.dot(n, c))
        tot_a += a
        tot_i += a*sum(vals[i] for i in tri)/3.0
    chk(tot_a, math.pi, 5e-3, "the visible hemisphere projects to pi")
    chk(tot_i, C.lambert_integral(1.0), 6e-3, "numeric integral == (2/3) pi")
    chk(tot_i/tot_a, C.LAMBERT_MEAN, 3e-3, "mean cosine == 2/3")

    print()
    if fails:
        print("%d analytic identities FAILED - the expectations themselves are wrong" % fails)
        return 1
    print("all analytic identities hold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
