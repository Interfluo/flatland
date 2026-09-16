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


# --------------------------------------------------------------------------
# Independent quadrature for the blackbody sphere cases
#
# These do the AZIMUTHAL integral analytically and the polar one numerically, a
# different route to the same number than the closed forms in cases.py. Both
# return the dimensionless integral over the visible hemisphere of
#
#     w(phat) (phat . ehat) dOmega
#
# so multiplying by L r^2 gives a radiant intensity. ehat is taken as +z and the
# source is tilted by alpha, which costs no generality.
# --------------------------------------------------------------------------

def _quad_equilibrium(alpha, n=20000):
    """w = max(phat . shat, 0) - radiative equilibrium, with a terminator.

    At polar angle theta the source cosine is B cos(phi) + A with
    A = cos(alpha) cos(theta) and B = sin(alpha) sin(theta). It is positive on
    |phi| < phi0 = arccos(-A/B), and integrating B cos(phi) + A over that arc
    gives 2[B sin(phi0) + A phi0]. The lit cap (A >= B) and the dark cap
    (A <= -B) are the two degenerate ends.
    """
    ca, sa = math.cos(alpha), math.sin(alpha)
    total = 0.0
    dmu = 1.0/n
    for i in range(n):
        mu = (i + 0.5)*dmu                       # mu = cos(theta), visible half
        a_term = ca*mu
        b_term = sa*math.sqrt(max(0.0, 1.0 - mu*mu))
        if a_term >= b_term:
            inner = 2.0*math.pi*a_term           # fully lit at this latitude
        elif a_term <= -b_term:
            inner = 0.0                          # fully dark
        else:
            phi0 = math.acos(max(-1.0, min(1.0, -a_term/b_term)))
            inner = 2.0*(b_term*math.sin(phi0) + a_term*phi0)
        total += mu*inner
    return total*dmu


def _quad_graded(alpha, n=20000):
    """w = (1 + phat . shat)/2 - affine in position, no terminator.

    The cos(phi) term integrates to zero over the full circle, so the azimuthal
    integral is just pi(1 + cos(alpha) cos(theta)).
    """
    ca = math.cos(alpha)
    total = 0.0
    dmu = 1.0/n
    for i in range(n):
        mu = (i + 0.5)*dmu
        total += mu*math.pi*(1.0 + ca*mu)
    return total*dmu



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

    print("blackbody constants")
    # sigma is not a measured number any more: the 2019 SI fixes h, k and c
    # exactly, so sigma = 2 pi^5 k^4 / (15 h^3 c^2) is exact too. Rebuild it from
    # the defining constants rather than trusting the literal in cases.py.
    k_B, h_P, c_0 = 1.380649e-23, 6.62607015e-34, 299792458.0
    sigma = 2*math.pi**5*k_B**4/(15*h_P**3*c_0**2)
    chk(C.STEFAN_BOLTZMANN/sigma, 1.0, 1e-14, "sigma == 2 pi^5 k^4 / 15 h^3 c^2")
    chk(C.blackbody_radiance(1000.0)*math.pi/C.blackbody_exitance(1000.0), 1.0, 1e-14,
        "L == M / pi")
    # Two independent textbook magnitudes, so a slipped exponent cannot hide.
    chk(C.blackbody_exitance(5778.0)/6.3196e7, 1.0, 1e-3,
        "M(5778 K, the solar Teff) == 63.2 MW/m^2")
    chk(C.blackbody_exitance(300.0)/459.3, 1.0, 1e-3, "M(300 K) == 459 W/m^2")

    print("the Lambert-sphere phase function")
    chk(C.lambert_phase(0.0), 1.0, 1e-15, "Phi(0) == 1, normalised at full phase")
    chk(C.lambert_phase(math.pi/2), 1.0/math.pi, 1e-15, "Phi(pi/2) == 1/pi")
    chk(C.lambert_phase(math.pi), 0.0, 1e-15, "Phi(pi) == 0, nothing lit is visible")
    # Monotone decreasing: a sphere cannot brighten as it moves away from full.
    phis = [C.lambert_phase(math.pi*i/200.0) for i in range(201)]
    chk(float(all(b <= a + 1e-15 for a, b in zip(phis, phis[1:]))), 1.0, 0.0,
        "Phi decreases monotonically over [0, pi]")

    print("blackbody sphere intensities vs numerical quadrature")
    # An independent route to the same numbers: do the azimuthal integral
    # ANALYTICALLY and the polar one numerically. Nothing here reuses the closed
    # forms being checked, so agreement is evidence, not circularity.
    T, R = 1200.0, 1.0
    L = C.blackbody_radiance(T)
    for deg in (0, 45, 90, 135, 175):
        a = math.radians(deg)
        chk(L*R*R*_quad_equilibrium(a)/C.equilibrium_sphere_intensity(T, R, a), 1.0, 1e-5,
            "radiative equilibrium: I(%d deg) == (2/3) sigma T^4 r^2 Phi" % deg)
    for deg in (0, 45, 90, 135, 180):
        a = math.radians(deg)
        chk(L*R*R*_quad_graded(a)/C.graded_sphere_intensity(T, R, a), 1.0, 1e-5,
            "graded T^4: I(%d deg) == (sigma T^4 r^2/2)(1 + 2cos/3)" % deg)

    print("the three sphere cases agree where they overlap")
    iso = C.isothermal_sphere_intensity(T, R)
    chk(C.equilibrium_sphere_intensity(T, R, 0.0)/iso, 2.0/3.0, 1e-14,
        "equilibrium at full phase == (2/3) x isothermal")
    chk(C.graded_sphere_intensity(T, R, math.pi/2)/iso, 0.5, 1e-14,
        "graded seen edge-on == half the isothermal value")
    chk(C.isothermal_intensity(T, math.pi*R*R)/iso, 1.0, 1e-14,
        "I = L A_proj reduces to sigma T^4 r^2 on a sphere")

    print("Stefan-Boltzmann recovered from projected areas alone")
    # The surface area is never used in the left-hand side: only projected areas
    # and Cauchy's identity. Getting sigma T^4 S back out is the check.
    T = 420.0
    # As in cauchy_mean, the direction count is scaled to the mesh: only the
    # SAMPLING is approximate, and an isotropic sphere needs far fewer samples
    # to sit inside the tolerance than a cube with six flat faces does.
    for name, (v, f), n_dirs in (("cube", C.unit_cube(), 8000),
                                 ("octahedron", C.regular_octahedron(1.0), 8000),
                                 ("icosphere(2)", C.icosphere(2, 1.0), 600)):
        dirs = C.fibonacci_directions(n_dirs)
        mean_I = sum(C.isothermal_intensity(T, C.projected_area_convex(v, f, d))
                     for d in dirs)/len(dirs)
        chk(4*math.pi*mean_I/C.blackbody_power(v, f, T), 1.0, 1e-4,
            "%s: 4 pi <I> == sigma T^4 S" % name)

    print("polyhedral integrals converge on the smooth closed forms")
    # The exact answer for the MESH, against the exact answer for the SPHERE the
    # mesh approximates. They must agree to the mesh's fidelity - and the gap
    # must shrink as the mesh improves, which is the next loop.
    src = (0.0, 0.0, 1.0)
    v, f = C.icosphere(4, 1.0)
    T, R = 1200.0, 1.0
    vg = C.graded_radiance_field(v, src, T)
    ve = C.equilibrium_radiance_field(v, src, T)
    for deg, tol in ((0, 3e-3), (60, 3e-3), (120, 3e-3), (150, 1e-2)):
        a = math.radians(deg)
        n = C.phase_view(src, a)
        Ig = C.field_integral_convex(v, f, vg, n)[0]
        Ie = C.field_integral_convex(v, f, ve, n)[0]
        chk(Ig/C.graded_sphere_intensity(T, R, a), 1.0, tol,
            "graded: subdiv-4 mesh integral -> smooth value (%d deg)" % deg)
        chk(Ie/C.equilibrium_sphere_intensity(T, R, a), 1.0, tol,
            "equilibrium: subdiv-4 mesh integral -> smooth value (%d deg)" % deg)
    # Refinement must actually help, at the first-order rate an inscribed
    # polyhedron gives. A formula error would show up as a floor instead.
    a = math.radians(60.0)
    n = C.phase_view(src, a)
    errs = []
    for k in (2, 3, 4):
        v, f = C.icosphere(k, 1.0)
        vals = C.graded_radiance_field(v, src, T)
        I = C.field_integral_convex(v, f, vals, n)[0]
        errs.append(abs(I/C.graded_sphere_intensity(T, R, a) - 1.0))
    chk(float(errs[0] > errs[1] > errs[2]), 1.0, 0.0,
        "graded: the gap to the sphere shrinks with subdivision")
    chk(errs[1]/errs[2], 4.0, 1.0, "...and does so roughly fourfold per level")

    print()
    if fails:
        print("%d analytic identities FAILED - the expectations themselves are wrong" % fails)
        return 1
    print("all analytic identities hold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
