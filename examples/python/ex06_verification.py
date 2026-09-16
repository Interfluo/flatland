#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
# for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

"""
Example 6 — Checking FlatLand against closed-form answers.

The question this answers: *how do I know the number is right?*

The recipe is always the same. Find a geometry whose answer you can write down,
run FlatLand on it at the resolution you intend to use for real, and compare.
Then refine the mesh and the pixel size separately, because they are different
error sources and only one of them is usually worth spending effort on.

    python3 examples/python/ex06_verification.py

Covered:
  - projected area of a cube: A(n) = |nx| + |ny| + |nz|
  - Cauchy's identity <A> = S/4, over many directions at once
  - blackbody radiant intensity, and the Stefan-Boltzmann law recovered from it
  - the phase curve of a sphere in radiative equilibrium
  - separating rasterization error from mesh error, which is the point

The same closed forms drive the project's own validation suite. See
docs/VALIDATION.md, and validation/cases.py for the derivations.
"""

import math

from _setup import load, rule

fl = load()

# sigma is exact in the 2019 SI: h, k and c are defined constants, so
# sigma = 2 pi^5 k^4 / (15 h^3 c^2) carries no experimental uncertainty.
SIGMA = 5.670374419184431e-8                       # W m^-2 K^-4


def unit(p):
    n = math.sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
    return (p[0]/n, p[1]/n, p[2]/n)


def icosphere(subdivisions=3, radius=1.0):
    """An icosahedron subdivided and pushed onto a sphere. Inscribed, so it
    genuinely projects to slightly LESS than pi r^2 — which is a property of the
    mesh, not an error in the tool. Keeping that straight is half this example."""
    t = (1.0 + math.sqrt(5.0)) / 2.0
    raw = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0),
           (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
           (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    faces = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
             (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
             (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
             (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    verts = [list(unit(p)) for p in raw]
    for _ in range(subdivisions):
        cache, out = {}, []

        def mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                verts.append(list(unit(tuple(verts[a][i] + verts[b][i] for i in range(3)))))
                cache[key] = len(verts) - 1
            return cache[key]

        for (a, b, c) in faces:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            out += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        faces = out
    return [tuple(c*radius for c in p) for p in verts], faces


def unit_cube():
    v = [(x, y, z) for x in (-.5, .5) for y in (-.5, .5) for z in (-.5, .5)]
    quads = [(0, 1, 3, 2), (4, 6, 7, 5), (0, 4, 5, 1),
             (2, 3, 7, 6), (0, 2, 6, 4), (1, 5, 7, 3)]
    return v, [t for (a, b, c, d) in quads for t in ((a, b, c), (a, c, d))]


def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])


def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def face_cross(verts, f):
    v0, v1, v2 = verts[f[0]], verts[f[1]], verts[f[2]]
    return cross(tuple(v1[i]-v0[i] for i in range(3)),
                 tuple(v2[i]-v0[i] for i in range(3)))


def projected_area_convex(verts, faces, direction):
    """A(n) = 1/4 sum |n . c_i| — exact for any closed convex polyhedron.

    Every line through a convex body crosses the surface twice, so each direction
    sees half the total projected face area. This is THE reference to compare
    against: it is the exact answer for the mesh you actually handed over, with
    no reference to whatever smooth shape that mesh approximates.
    """
    n = unit(direction)
    return 0.25*sum(abs(dot(n, face_cross(verts, f))) for f in faces)


def field_integral_convex(verts, faces, values, direction):
    """Exact area integral of a per-vertex field over the visible side.

    Barycentric interpolation reproduces a linear function exactly and the mean
    of a linear function over a triangle is the mean of its vertex values, so a
    front-facing triangle contributes (its projected area) x (that mean).
    """
    n = unit(direction)
    total = 0.0
    for f in faces:
        c = face_cross(verts, f)
        if dot(c, n) >= 0:
            continue                               # back-facing: camera looks ALONG n
        total += 0.5*abs(dot(n, c))*(values[f[0]] + values[f[1]] + values[f[2]])/3.0
    return total


def fibonacci_directions(count):
    """Quasi-uniform directions on the sphere. Deterministic, and converges far
    faster than random sampling for the direction averages below."""
    ga = math.pi*(3.0 - math.sqrt(5.0))
    out = []
    for i in range(count):
        z = 1.0 - (2.0*i + 1.0)/count
        r = math.sqrt(max(0.0, 1.0 - z*z))
        out.append((r*math.cos(i*ga), r*math.sin(i*ga), z))
    return out


def report(label, measured, exact, tol):
    rel = abs(measured - exact)/abs(exact) if exact else abs(measured)
    print("  %-34s %14.6f %14.6f %10.2e  %s"
          % (label, measured, exact, rel, "ok" if rel <= tol else "OFF"))
    return rel <= tol


def header(a="measured", b="closed form"):
    print("  %-34s %14s %14s %10s" % ("case", a, b, "rel err"))
    print("  %-34s %14s %14s %10s" % ("-"*34, "-"*14, "-"*14, "-"*10))


def main():
    ok = True

    # ---------------------------------------------------------------------
    rule("1. Projected area of a cube")
    # For a unit cube the general convex formula collapses to something you can
    # check in your head: A(n) = |nx| + |ny| + |nz| for a UNIT direction. Down a
    # face normal that is 1; down a face diagonal sqrt(2); down the body diagonal
    # sqrt(3), the regular hexagon.
    verts, faces = unit_cube()
    mesh = fl.Mesh(verts, faces)
    header()
    for d, name in (((1, 0, 0), "face normal"),
                    ((1, 1, 0), "face diagonal"),
                    ((1, 1, 1), "body diagonal (hexagon)"),
                    ((0.3, 0.9, -0.31), "generic direction")):
        n = unit(d)
        exact = abs(n[0]) + abs(n[1]) + abs(n[2])
        ok &= report(name, mesh.project(d, resolution=5e-4, precision="double").area,
                     exact, 2e-3)

    # ---------------------------------------------------------------------
    rule("2. Cauchy's identity: <A> = S/4")
    # An integral identity, so it probes many directions at once rather than a
    # few hand-picked ones. A cube has S = 6, so the mean projected area is 1.5 —
    # and note that no single direction gives 1.5.
    dirs = fibonacci_directions(160)
    areas = [r.area for r in mesh.project_batch(dirs, resolution=1e-3, precision="double")]
    header()
    ok &= report("cube, mean over 160 directions", sum(areas)/len(areas), 1.5, 3e-3)
    print("    (the spread is %.4f .. %.4f — the identity is about the MEAN)"
          % (min(areas), max(areas)))
    mesh.close()

    # ---------------------------------------------------------------------
    rule("3. Blackbody radiant intensity")
    # FlatLand integrates over the PROJECTED area, and dA_proj = cos(theta) dA.
    # So if the field is a radiance L, the integral IS the radiant intensity:
    #
    #     I(n) = integral L cos(theta) dA   [W/sr]
    #
    # A blackbody is a Lambertian emitter with L = sigma T^4 / pi, hence for an
    # isothermal body I = L x A_proj, and for a sphere the pi cancels:
    #
    #     I = (sigma T^4 / pi)(pi r^2) = sigma T^4 r^2,  the same from every side.
    T = 1200.0
    radiance = SIGMA*T**4/math.pi
    sv, sf = icosphere(4, 1.0)
    sphere = fl.Mesh(sv, sf)
    field = [radiance]*len(sv)
    exact_I = SIGMA*T**4*1.0**2

    print("  isothermal sphere at %.0f K: L = %.2f W/m2/sr\n" % (T, radiance))
    header("FlatLand", "sigma T^4 r^2")
    for d in ((0, 0, 1), (1, 1, 1), (-0.4, 0.9, 0.2)):
        r = sphere.project(d, field=field, resolution=1e-3, precision="double")
        ok &= report("I along (%.1f, %.1f, %.1f)" % d, r.integral, exact_I, 3e-3)
    # The MEAN is a different kind of check. A constant field must interpolate to
    # exactly that constant at every covered pixel, whatever shape those pixels
    # cover, so the mean carries no mesh error at all — it comes back to ten
    # digits, while the integral above is held to about three by the inscribed
    # mesh's projected-area deficit.
    r = sphere.project((0, 0, 1), field=field, resolution=1e-3, precision="double")
    ok &= report("mean radiance (interpolation)", r.average, radiance, 1e-9)

    # ---------------------------------------------------------------------
    rule("4. Stefan-Boltzmann, recovered from projected areas")
    # Integrate the intensity over all directions and apply Cauchy:
    #
    #     integral I dOmega = (sigma T^4/pi)(4 pi)(S/4) = sigma T^4 S
    #
    # The left side never uses the surface area. Getting sigma T^4 S back out of
    # a pile of projected areas is the check — it is the Stefan-Boltzmann law,
    # reassembled from geometry.
    dirs = fibonacci_directions(300)
    rs = sphere.project_batch(dirs, field_matrix=[[x] for x in field],
                              resolution=2e-3, precision="double")
    power = 4.0*math.pi*sum(x.integral for x in rs)/len(rs)
    surface = sum(0.5*math.sqrt(dot(face_cross(sv, f), face_cross(sv, f))) for f in sf)
    header("4 pi <I> [W]", "sigma T^4 S")
    ok &= report("icosphere(4), 300 directions", power, SIGMA*T**4*surface, 5e-3)
    print("    (mesh surface area %.6f, vs 4 pi = %.6f for a true sphere)"
          % (surface, 4*math.pi))
    sphere.close()

    # ---------------------------------------------------------------------
    rule("5. A sphere in radiative equilibrium — the phase curve")
    # Absorbed flux goes as the incidence cosine, so balancing it against
    # sigma T^4 gives the subsolar law T(psi) = T_sub cos^{1/4}(psi) on the lit
    # side and nothing beyond the terminator. The disc-integrated result is the
    # Lambert-sphere phase function (Russell 1916):
    #
    #     I(alpha) = (2/3) sigma T_sub^4 r^2 Phi(alpha)
    #     Phi(alpha) = [sin(alpha) + (pi - alpha) cos(alpha)] / pi
    #
    # Two references are shown, and the gap between them is the whole lesson: the
    # middle column is exact FOR THIS MESH, so the "raster" error is the tool's
    # own; "vs sphere" adds the mesh's fidelity to a true sphere on top.
    sv, sf = icosphere(4, 1.0)
    sphere = fl.Mesh(sv, sf)
    src = (0.0, 0.0, 1.0)
    vals = [radiance*max(0.0, dot(unit(p), src)) for p in sv]
    base = (2.0/3.0)*SIGMA*T**4

    print("  icosphere(4), %d triangles, %.0f K subsolar\n" % (len(sf), T))
    print("  %6s %14s %14s %14s %10s %10s"
          % ("phase", "I FlatLand", "exact (mesh)", "exact (sphere)", "raster", "vs sphere"))
    print("  %6s %14s %14s %14s %10s %10s"
          % ("-"*6, "-"*14, "-"*14, "-"*14, "-"*10, "-"*10))
    for deg in (0, 30, 60, 90, 120):
        a = math.radians(deg)
        # The observer sits at angle alpha from the source; the camera looks the
        # other way, back at the body.
        eye = (math.sin(a), 0.0, math.cos(a))
        view = (-eye[0], -eye[1], -eye[2])
        r = sphere.project(view, field=vals, resolution=1e-3, precision="double")
        mesh_exact = field_integral_convex(sv, sf, vals, view)
        phi = (math.sin(a) + (math.pi - a)*math.cos(a))/math.pi
        smooth = base*phi
        raster = abs(r.integral - mesh_exact)/mesh_exact
        print("  %5d° %14.4f %14.4f %14.4f %10.2e %10.2e"
              % (deg, r.integral, mesh_exact, smooth, raster,
                 abs(r.integral - smooth)/smooth))
        ok &= raster <= 2e-3
    sphere.close()

    # ---------------------------------------------------------------------
    rule("6. Which error should you spend effort on?")
    # Refine the mesh and the pixel size independently. Only one of them is
    # usually worth the money, and this table tells you which.
    print("  %8s %8s %8s %14s %12s %12s"
          % ("subdiv", "tris", "res", "area", "vs mesh", "vs pi"))
    print("  %8s %8s %8s %14s %12s %12s"
          % ("-"*8, "-"*8, "-"*8, "-"*14, "-"*12, "-"*12))
    for k in (1, 2, 3, 4):
        v, f = icosphere(k, 1.0)
        m = fl.Mesh(v, f)
        for res in (4e-3, 1e-3):
            r = m.project((1, 0.3, 0.2), resolution=res, precision="double")
            exact = projected_area_convex(v, f, (1, 0.3, 0.2))
            print("  %8d %8d %8.0e %14.7f %12.2e %12.2e"
                  % (k, len(f), res, r.area, abs(r.area - exact)/exact,
                     abs(r.area - math.pi)/math.pi))
        m.close()
    print("""
  Read the last two columns against each other.

  "vs mesh" is FlatLand against the exact answer for the polyhedron it was
  actually handed, so it is the tool's own error. It stays around 1e-5 to 1e-6
  throughout and bounces around rather than falling cleanly: coverage is decided
  by point-sampling pixel centres, and the boundary band that error lives in is
  not a smooth function of the pixel size. Expect a trend, not a rate.

  "vs pi" is the same runs against a true sphere, so it adds the mesh's own
  fidelity. It falls fourfold per subdivision, 7e-2 down to 1e-3, and would keep
  going until it reached the floor the first column sets.

  On this mesh the verdict is unambiguous: triangles are the limit, not pixels.
  At subdiv 1 a 4x finer raster moves the projected area in the fourth decimal
  and leaves you 7% away from a sphere. Run this table on your own geometry
  before deciding where to spend — guessing usually picks wrong.""")

    print()
    if ok:
        print("every closed-form check agreed within tolerance.")
        return 0
    print("SOME CHECKS WERE OUTSIDE TOLERANCE — see the OFF rows above.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
