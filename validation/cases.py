#!/usr/bin/env python3
"""Analytically-known test cases for validating FlatLand.

Dependency-free: standard library only, matching the rest of the project.

The point of this module is that every case here has a CLOSED-FORM answer, so a
FlatLand run can be checked against mathematics rather than against a previous
run of FlatLand. Regression tests catch changes; these catch being wrong.

----------------------------------------------------------------------------
THE THREE RESULTS EVERYTHING ELSE IS BUILT ON

1. Projected area of a convex polyhedron.

   For a closed convex surface, each direction sees exactly half the total
   projected face area, because every line through the body crosses the surface
   exactly twice:

       A(n) = 1/2 * sum_i A_i |n . nhat_i|

   With A_i = 1/2 |c_i| and nhat_i = c_i/|c_i| for the triangle cross product
   c_i = (v1-v0) x (v2-v0), this collapses to

       A(n) = 1/4 * sum_i |n . c_i|                              (exact)

   No integration, no discretization: for a mesh that IS the polyhedron, this is
   the exact answer FlatLand should converge to as the pixel size shrinks. That
   separation matters — it lets us test the RASTERIZER without the mesh's
   faithfulness to a curved shape getting mixed in.

2. Cauchy's projection formula.

   Averaged over directions uniform on the sphere, the projected area of any
   convex body is a quarter of its surface area:

       <A> = S/4                                                 (exact)

   This is an integral identity, so it tests the tool across many directions at
   once rather than at hand-picked ones.

3. Mean of a linear field over a triangle.

   Barycentric interpolation reproduces a linear function exactly, and the mean
   of a linear function over a triangle is the mean of its vertex values:

       (1/A) * integral f dA = (f0 + f1 + f2)/3                  (exact)

   So the area integral over a planar region is sum_tri A_tri (f0+f1+f2)/3.

----------------------------------------------------------------------------
CONVENTIONS

A view direction is the direction the camera LOOKS ALONG, so the visible surface
of a closed mesh is the one whose outward normals oppose it. Fields are given
per vertex, in mesh vertex order, one value per line.
"""

import math

# --------------------------------------------------------------------------
# Small vector helpers
# --------------------------------------------------------------------------

def sub(a, b):   return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def add(a, b):   return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def scale(a, s): return (a[0]*s, a[1]*s, a[2]*s)
def dot(a, b):   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def norm(a):     return math.sqrt(dot(a, a))
def unit(a):
    n = norm(a)
    if n == 0:
        raise ValueError("cannot normalize the zero vector")
    return scale(a, 1.0/n)

# --------------------------------------------------------------------------
# Exact analytic quantities
# --------------------------------------------------------------------------

def face_cross(verts, face):
    v0, v1, v2 = (verts[i] for i in face)
    return cross(sub(v1, v0), sub(v2, v0))


def surface_area(verts, faces):
    """Total surface area: sum of triangle areas."""
    return 0.5 * sum(norm(face_cross(verts, f)) for f in faces)


def projected_area_convex(verts, faces, direction):
    """Exact projected area of a CLOSED CONVEX mesh along `direction`.

    A(n) = 1/4 sum_i |n . c_i|.  Only valid for convex, closed surfaces; for an
    open sheet use projected_area_open().
    """
    n = unit(direction)
    return 0.25 * sum(abs(dot(n, face_cross(verts, f))) for f in faces)


def projected_area_open(verts, faces, direction):
    """Exact projected area of an open, non-self-occluding sheet.

    Each triangle contributes A_i |n . nhat_i| = 1/2 |n . c_i|, with no factor of
    one half for front/back pairing because there is only one layer.
    """
    n = unit(direction)
    return 0.5 * sum(abs(dot(n, face_cross(verts, f))) for f in faces)


def field_integral_planar(verts, faces, values, direction):
    """Exact area integral of a per-vertex field over a flat, single-layer mesh.

    Barycentric interpolation is exact for a linear function, so each triangle
    contributes (its projected area) x (the mean of its three vertex values).
    Returns (integral, area, mean).
    """
    n = unit(direction)
    total_a = 0.0
    total_i = 0.0
    for f in faces:
        a = 0.5 * abs(dot(n, face_cross(verts, f)))
        mean = (values[f[0]] + values[f[1]] + values[f[2]]) / 3.0
        total_a += a
        total_i += a * mean
    return total_i, total_a, (total_i/total_a if total_a else float("nan"))


def fibonacci_directions(count):
    """`count` quasi-uniform directions on the sphere (Fibonacci spiral).

    Used for the Cauchy mean-projected-area study. A spiral converges far faster
    than random sampling and is deterministic, so the study is reproducible.
    """
    out = []
    ga = math.pi * (3.0 - math.sqrt(5.0))      # golden angle
    for i in range(count):
        z = 1.0 - 2.0*(i + 0.5)/count
        r = math.sqrt(max(0.0, 1.0 - z*z))
        th = ga * i
        out.append((r*math.cos(th), r*math.sin(th), z))
    return out

# --------------------------------------------------------------------------
# Meshes. Every one is exact: the triangles ARE the shape, with no curved
# surface being approximated, except icosphere() and cylinder(), which
# approximate deliberately and whose exact polyhedral values are still
# available through projected_area_convex().
# --------------------------------------------------------------------------

def unit_cube():
    """Axis-aligned unit cube spanning [0,1]^3, outward CCW normals.

    Analytic: A(n) = |nx| + |ny| + |nz| for unit n.  S = 6, so <A> = 1.5.
    """
    v = [(0,0,0), (1,0,0), (1,1,0), (0,1,0),
         (0,0,1), (1,0,1), (1,1,1), (0,1,1)]
    f = [(0,3,2), (0,2,1),      # z = 0, normal -z
         (4,5,6), (4,6,7),      # z = 1, normal +z
         (0,1,5), (0,5,4),      # y = 0, normal -y
         (1,2,6), (1,6,5),      # x = 1, normal +x
         (2,3,7), (2,7,6),      # y = 1, normal +y
         (3,0,4), (3,4,7)]      # x = 0, normal -x
    return v, f


def cube_projected_area(direction):
    """Closed form for the unit cube: A(n) = |nx| + |ny| + |nz|."""
    n = unit(direction)
    return abs(n[0]) + abs(n[1]) + abs(n[2])


def regular_tetrahedron(edge=1.0):
    """Regular tetrahedron with the given edge length, outward CCW normals.

    Analytic: S = sqrt(3) * edge^2, so <A> = sqrt(3) edge^2 / 4.
    """
    s = edge / math.sqrt(8.0)      # places vertices at (+-s,+-s,+-s) with edge `edge`
    v = [( s,  s,  s), ( s, -s, -s), (-s,  s, -s), (-s, -s,  s)]
    f = [(0,1,2), (0,3,1), (0,2,3), (1,3,2)]
    return v, _orient_outward(v, f)


def regular_octahedron(radius=1.0):
    """Regular octahedron with vertices at distance `radius` along each axis."""
    r = radius
    v = [( r,0,0), (-r,0,0), (0, r,0), (0,-r,0), (0,0, r), (0,0,-r)]
    f = [(0,2,4), (2,1,4), (1,3,4), (3,0,4),
         (2,0,5), (1,2,5), (3,1,5), (0,3,5)]
    return v, _orient_outward(v, f)


def icosphere(subdivisions=2, radius=1.0):
    """Icosahedron subdivided `subdivisions` times and projected onto a sphere.

    This is the one deliberately inexact shape: it UNDERESTIMATES the sphere,
    because an inscribed polyhedron sits inside it. Its own exact projected area
    is still available from projected_area_convex(), which lets the two error
    sources — mesh fidelity and rasterization — be measured separately.
    """
    t = (1.0 + math.sqrt(5.0)) / 2.0
    v = [(-1, t, 0), ( 1, t, 0), (-1,-t, 0), ( 1,-t, 0),
         ( 0,-1, t), ( 0, 1, t), ( 0,-1,-t), ( 0, 1,-t),
         ( t, 0,-1), ( t, 0, 1), (-t, 0,-1), (-t, 0, 1)]
    f = [(0,11,5), (0,5,1), (0,1,7), (0,7,10), (0,10,11),
         (1,5,9), (5,11,4), (11,10,2), (10,7,6), (7,1,8),
         (3,9,4), (3,4,2), (3,2,6), (3,6,8), (3,8,9),
         (4,9,5), (2,4,11), (6,2,10), (8,6,7), (9,8,1)]
    v = [list(unit(p)) for p in v]

    for _ in range(subdivisions):
        midpoint = {}
        new_f = []

        def mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in midpoint:
                m = unit(add(v[a], v[b]))       # push the new point onto the sphere
                v.append(list(m))
                midpoint[key] = len(v) - 1
            return midpoint[key]

        for (a, b, c) in f:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            new_f += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        f = new_f

    v = [tuple(scale(p, radius)) for p in v]
    return v, _orient_outward(v, f)


def cylinder(sides=64, radius=1.0, height=2.0):
    """Closed right prism approximating a cylinder, axis along +z, centred.

    For the smooth cylinder viewed at angle theta to its axis the projected area
    is  2*r*h*sin(theta) + pi*r^2*cos(theta);  the prism converges to it as
    `sides` grows. The prism's own exact value comes from
    projected_area_convex().
    """
    hz = height / 2.0
    v = []
    for i in range(sides):
        a = 2.0*math.pi*i/sides
        v.append((radius*math.cos(a), radius*math.sin(a), -hz))
    for i in range(sides):
        a = 2.0*math.pi*i/sides
        v.append((radius*math.cos(a), radius*math.sin(a), hz))
    bot_c = len(v); v.append((0.0, 0.0, -hz))
    top_c = len(v); v.append((0.0, 0.0,  hz))

    f = []
    for i in range(sides):
        j = (i + 1) % sides
        f.append((i, j, sides + j))                 # side wall
        f.append((i, sides + j, sides + i))
        f.append((bot_c, j, i))                     # bottom cap
        f.append((top_c, sides + i, sides + j))     # top cap
    return v, _orient_outward(v, f)


def cylinder_projected_area(radius, height, theta_rad):
    """Smooth-cylinder closed form: 2 r h sin(theta) + pi r^2 |cos(theta)|."""
    return (2.0*radius*height*abs(math.sin(theta_rad))
            + math.pi*radius*radius*abs(math.cos(theta_rad)))


def unit_square():
    """Unit square in the XY plane, +Z normal. Open sheet, two triangles."""
    v = [(0,0,0), (1,0,0), (1,1,0), (0,1,0)]
    f = [(0,1,2), (0,2,3)]
    return v, f


def right_triangle():
    """Right triangle (0,0), (1,0), (0,1) in the XY plane, +Z normal."""
    return [(0,0,0), (1,0,0), (0,1,0)], [(0,1,2)]


def _orient_outward(verts, faces):
    """Flip any triangle whose normal points toward the centroid.

    The analytic formulas take absolute values and so do not care, but FlatLand's
    backface culling does, and a mesh with inconsistent winding is not a fair
    test of it.
    """
    cx = sum(p[0] for p in verts)/len(verts)
    cy = sum(p[1] for p in verts)/len(verts)
    cz = sum(p[2] for p in verts)/len(verts)
    centre = (cx, cy, cz)
    out = []
    for f in faces:
        c = face_cross(verts, f)
        mid = scale(add(add(verts[f[0]], verts[f[1]]), verts[f[2]]), 1.0/3.0)
        out.append(f if dot(c, sub(mid, centre)) > 0 else (f[0], f[2], f[1]))
    return out

# --------------------------------------------------------------------------
# Fields with known integrals
# --------------------------------------------------------------------------

def lambert_field(verts, direction, radius=1.0):
    """f = cos(angle between the surface normal and the viewer) on a sphere.

    For a sphere centred at the origin the outward unit normal at a vertex is
    just the vertex direction, so f = -vhat . n is the cosine on the visible
    side. The analytic result is the classic disc-integrated Lambertian value:

        integral f dA = (2/3) pi r^2        mean = 2/3

    Derivation: at projected radius rho the cosine is sqrt(1 - (rho/r)^2), so
        integral = int_0^r sqrt(1-(rho/r)^2) 2 pi rho drho
                 = 2 pi r^2 int_0^1 u sqrt(1-u^2) du = (2/3) pi r^2.

    FlatLand knows nothing about optics here — this is just the area integral of
    a particular scalar field, which is the point of keeping it domain-agnostic.
    """
    n = unit(direction)
    return [-dot(unit(p), n) for p in verts]


LAMBERT_MEAN = 2.0/3.0


def lambert_integral(radius=1.0):
    return (2.0/3.0) * math.pi * radius * radius


def coordinate_field(verts, axis):
    """f = the given coordinate of each vertex. Linear, so exactly interpolated."""
    return [p[axis] for p in verts]

# --------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------

def write_obj(path, verts, faces):
    with open(path, "w") as fh:
        fh.write("# generated by validation/cases.py\n")
        for p in verts:
            fh.write("v %.17g %.17g %.17g\n" % p)
        for f in faces:
            fh.write("f %d %d %d\n" % (f[0]+1, f[1]+1, f[2]+1))


def write_field(path, values):
    with open(path, "w") as fh:
        for x in values:
            fh.write("%.17g\n" % x)


def write_batch(path, directions, resolution, data=None):
    with open(path, "w") as fh:
        fh.write("# <nx> <ny> <nz> <resolution> [data]\n")
        for d in directions:
            line = "%.17g %.17g %.17g %.17g" % (d[0], d[1], d[2], resolution)
            if data:
                line += " " + data
            fh.write(line + "\n")
