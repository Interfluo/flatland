#!/usr/bin/env python3
"""Run FlatLand against analytically-known answers and report the errors.

    python3 validation/run_validation.py [--flatland ./flatland] [--quick] [--json out.json]

Exits non-zero if any study exceeds its tolerance. `--quick` runs the coarse
subset the test suite uses; the default runs the full study, including the
convergence sweeps that back the numbers in docs/VALIDATION.md.

Every expected value here comes from validation/cases.py, which derives them in
closed form. Nothing is compared against a stored FlatLand output, so a wrong
answer cannot be blessed into a regression baseline.
"""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases as C


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

class Report:
    def __init__(self):
        self.rows = []
        self.failures = 0

    def check(self, study, case, measured, expected, tol, note=""):
        if expected == 0:
            rel = abs(measured)
        else:
            rel = abs(measured - expected) / abs(expected)
        ok = rel <= tol
        if not ok:
            self.failures += 1
        self.rows.append(dict(study=study, case=case, measured=measured,
                              expected=expected, rel_error=rel, tol=tol,
                              ok=ok, note=note))
        return ok

    def table(self, study):
        rows = [r for r in self.rows if r["study"] == study]
        if not rows:
            return
        print(f"\n{study}")
        print(f"  {'case':<38} {'measured':>13} {'expected':>13} {'rel err':>10}  ")
        print(f"  {'-'*38} {'-'*13} {'-'*13} {'-'*10}")
        for r in rows:
            mark = " " if r["ok"] else "  <-- EXCEEDS TOLERANCE"
            print(f"  {r['case']:<38} {r['measured']:>13.8f} {r['expected']:>13.8f} "
                  f"{r['rel_error']:>10.2e}{mark}")


# --------------------------------------------------------------------------
# Driving FlatLand
# --------------------------------------------------------------------------

class Runner:
    def __init__(self, binary, workdir):
        self.binary = binary
        self.workdir = workdir

    def run(self, mesh, args):
        cmd = [self.binary, mesh] + list(args) + ["-j"]
        p = subprocess.run(cmd, capture_output=True, text=True, cwd=self.workdir)
        if p.returncode != 0:
            raise RuntimeError(f"flatland failed ({p.returncode}): {' '.join(cmd)}\n{p.stderr}")
        try:
            return json.loads(p.stdout)["results"]
        except (ValueError, KeyError) as exc:
            raise RuntimeError(f"could not parse output of {' '.join(cmd)}: {exc}\n{p.stdout[:400]}")

    def batch(self, mesh, directions, resolution, data=None, extra=()):
        """One process, many views — the same path a real time series takes."""
        bpath = os.path.join(self.workdir, "batch.txt")
        C.write_batch(bpath, directions, resolution, data)
        return self.run(mesh, ["-b", "batch.txt"] + list(extra))


# --------------------------------------------------------------------------
# Studies
# --------------------------------------------------------------------------

def study_convex_polyhedra(run, rep, quick):
    """Projected area of exact polyhedra, where A(n) = 1/4 sum |n . c_i| is exact.

    These meshes ARE the shape, so any discrepancy is rasterization error and
    nothing else.
    """
    shapes = [
        ("unit cube",     C.unit_cube(),               0.0015),
        ("tetrahedron",   C.regular_tetrahedron(1.0),  0.004),
        ("octahedron",    C.regular_octahedron(1.0),   0.002),
    ]
    dirs = [(1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,1,1), (2,-1,3), (0.3,0.9,-0.31)]
    res = 0.002 if quick else 0.0005

    for name, (v, f), tol in shapes:
        mesh = f"{name.replace(' ','_')}.obj"
        C.write_obj(os.path.join(run.workdir, mesh), v, f)
        results = run.batch(mesh, dirs, res)
        for d, r in zip(dirs, results):
            expected = C.projected_area_convex(v, f, d)
            rep.check("Projected area of exact convex polyhedra (vs 1/4 sum |n.c|)",
                      f"{name}  n={d}", r["area"], expected, tol)


def study_cube_closed_form(run, rep, quick):
    """The unit cube has an especially simple closed form: A(n) = |nx|+|ny|+|nz|."""
    v, f = C.unit_cube()
    mesh = "cube_cf.obj"
    C.write_obj(os.path.join(run.workdir, mesh), v, f)
    dirs = C.fibonacci_directions(12 if quick else 40)
    res = 0.002 if quick else 0.001
    results = run.batch(mesh, dirs, res)
    worst = 0.0
    for d, r in zip(dirs, results):
        expected = C.cube_projected_area(d)
        worst = max(worst, abs(r["area"] - expected)/expected)
    rep.check("Unit cube against |nx|+|ny|+|nz|",
              f"worst of {len(dirs)} directions (as a ratio)", 1.0 + worst, 1.0, 0.002,
              note=f"worst relative error {worst:.2e}")


def study_cauchy(run, rep, quick):
    """Cauchy: averaged over the sphere, a convex body's projected area is S/4.

    An integral identity rather than a point check, so it probes the tool over
    many directions at once.
    """
    n_dirs = 64 if quick else 400
    dirs = C.fibonacci_directions(n_dirs)
    res = 0.004 if quick else 0.002
    shapes = [
        ("unit cube",    C.unit_cube(),          0.004),
        ("tetrahedron",  C.regular_tetrahedron(1.0), 0.01),
        ("octahedron",   C.regular_octahedron(1.0),  0.006),
        ("icosphere(3)", C.icosphere(3, 1.0),    0.004),
    ]
    for name, (v, f), tol in shapes:
        mesh = f"cauchy_{name.split('(')[0].replace(' ','_')}.obj"
        C.write_obj(os.path.join(run.workdir, mesh), v, f)
        results = run.batch(mesh, dirs, res)
        mean = sum(r["area"] for r in results)/len(results)
        rep.check("Cauchy mean projected area == S/4",
                  f"{name}  ({n_dirs} directions)", mean, C.surface_area(v, f)/4.0, tol)


def study_sphere(run, rep, quick):
    """Sphere, separating the two error sources.

    An icosphere is an INSCRIBED polyhedron, so it genuinely projects to less
    than pi r^2 — that is a property of the mesh, not a defect in the tool.
    Comparing against the polyhedron's own exact value isolates rasterization
    error; comparing against pi shows the mesh converging.
    """
    levels = [1, 2, 3] if quick else [1, 2, 3, 4, 5]
    d = (0.577350269, 0.577350269, 0.577350269)
    res = 0.004 if quick else 0.001
    print("\n  sphere convergence (radius 1, view along the body diagonal)")
    print(f"  {'subdiv':>6} {'tris':>7} {'FlatLand':>12} {'exact poly':>12} "
          f"{'raster err':>11} {'vs pi':>11}")
    for k in levels:
        v, f = C.icosphere(k, 1.0)
        mesh = f"sphere{k}.obj"
        C.write_obj(os.path.join(run.workdir, mesh), v, f)
        r = run.batch(mesh, [d], res)[0]
        exact_poly = C.projected_area_convex(v, f, d)
        raster_err = abs(r["area"] - exact_poly)/exact_poly
        vs_pi = abs(r["area"] - math.pi)/math.pi
        print(f"  {k:>6} {len(f):>7} {r['area']:>12.7f} {exact_poly:>12.7f} "
              f"{raster_err:>11.2e} {vs_pi:>11.2e}")
        rep.check("Sphere: FlatLand vs the mesh's own exact projected area",
                  f"icosphere({k}), {len(f)} triangles", r["area"], exact_poly, 0.004)
    # The finest mesh should also be close to the true sphere.
    v, f = C.icosphere(levels[-1], 1.0)
    mesh = f"sphere{levels[-1]}.obj"
    r = run.batch(mesh, [d], res)[0]
    rep.check("Sphere: finest mesh vs the true sphere (pi r^2)",
              f"icosphere({levels[-1]}) vs pi", r["area"], math.pi,
              0.01 if quick else 0.004)


def study_cylinder(run, rep, quick):
    """Cylinder: 2 r h sin(theta) + pi r^2 cos(theta).

    Exercises a shape with both curved and flat surfaces at once, at angles where
    both terms contribute.
    """
    sides = 128 if quick else 512
    v, f = C.cylinder(sides, 1.0, 2.0)
    mesh = "cylinder.obj"
    C.write_obj(os.path.join(run.workdir, mesh), v, f)
    angles = [0, 15, 30, 45, 60, 75, 90]
    dirs = [(math.sin(math.radians(a)), 0.0, math.cos(math.radians(a))) for a in angles]
    res = 0.004 if quick else 0.001
    results = run.batch(mesh, dirs, res)
    for a, r in zip(angles, results):
        expected = C.cylinder_projected_area(1.0, 2.0, math.radians(a))
        rep.check("Cylinder: 2rh sin(theta) + pi r^2 cos(theta)",
                  f"theta = {a} deg", r["area"], expected, 0.01 if quick else 0.004)


def study_linear_fields(run, rep, quick):
    """Linear fields, where barycentric interpolation is EXACT.

    The mean of a linear field over a triangle is the mean of its vertex values,
    so these integrals have no discretization error at all in the field — only
    the raster's boundary sampling.
    """
    res = 0.002 if quick else 0.0005

    # Unit square, f = x.  integral x dA = 1/2, mean 1/2, min 0, max 1.
    v, f = C.unit_square()
    vals = C.coordinate_field(v, 0)
    C.write_obj(os.path.join(run.workdir, "square.obj"), v, f)
    C.write_field(os.path.join(run.workdir, "square_fx.txt"), vals)
    r = run.run("square.obj", ["-v", "0", "0", "-1", "-r", str(res),
                               "-d", "square_fx.txt", "--no-cull", "-p", "double"])[0]
    exp_i, exp_a, exp_m = C.field_integral_planar(v, f, vals, (0, 0, 1))
    rep.check("Linear field on a flat sheet (exact interpolation)",
              "unit square, f = x: area", r["area"], exp_a, 0.004)
    rep.check("Linear field on a flat sheet (exact interpolation)",
              "unit square, f = x: integral", r["integral"], exp_i, 0.004)
    rep.check("Linear field on a flat sheet (exact interpolation)",
              "unit square, f = x: mean", r["average"], exp_m, 0.004)

    # Right triangle, f = x.  integral = 1/6, area 1/2, mean 1/3.
    v, f = C.right_triangle()
    vals = C.coordinate_field(v, 0)
    C.write_obj(os.path.join(run.workdir, "tri.obj"), v, f)
    C.write_field(os.path.join(run.workdir, "tri_fx.txt"), vals)
    r = run.run("tri.obj", ["-v", "0", "0", "-1", "-r", str(res),
                            "-d", "tri_fx.txt", "--no-cull", "-p", "double"])[0]
    rep.check("Linear field on a flat sheet (exact interpolation)",
              "right triangle, f = x: integral", r["integral"], 1.0/6.0, 0.006)
    rep.check("Linear field on a flat sheet (exact interpolation)",
              "right triangle, f = x: mean", r["average"], 1.0/3.0, 0.004)

    # Unit cube, f = y, viewed along +x: the visible face is x = 0, so the mean
    # of y over it is 1/2 and the integral is 1/2. This one also depends on the
    # depth test picking the NEAR surface, so it fails loudly if culling inverts.
    v, f = C.unit_cube()
    vals = C.coordinate_field(v, 1)
    C.write_obj(os.path.join(run.workdir, "cube_fy.obj"), v, f)
    C.write_field(os.path.join(run.workdir, "cube_fy.txt"), vals)
    r = run.run("cube_fy.obj", ["-v", "1", "0", "0", "-r", str(res),
                                "-d", "cube_fy.txt", "-p", "double"])[0]
    rep.check("Linear field on a flat sheet (exact interpolation)",
              "unit cube, f = y seen along +x: mean", r["average"], 0.5, 0.004)
    rep.check("Linear field on a flat sheet (exact interpolation)",
              "unit cube, f = y seen along +x: integral", r["integral"], 0.5, 0.006)


def study_lambert(run, rep, quick):
    """The area integral of f = cos(incidence) over a sphere.

        integral cos dA = (2/3) pi r^2       mean = 2/3

    A non-trivial, smoothly varying field with a closed-form integral — the
    strongest single check of the interpolation-and-integration path. Supply a
    radiance field and this quantity is a radiant intensity; FlatLand itself
    attaches no meaning to it.
    """
    levels = [2, 3] if quick else [2, 3, 4, 5]
    d = (0.0, 0.0, 1.0)
    res = 0.004 if quick else 0.001
    print("\n  Lambertian sphere: integral of cos(incidence) dA over the visible disc")
    print(f"  {'subdiv':>6} {'tris':>7} {'integral':>12} {'exact':>12} {'mean':>10} "
          f"{'exact':>10} {'rel err':>10}")
    for k in levels:
        v, f = C.icosphere(k, 1.0)
        vals = C.lambert_field(v, d, 1.0)
        mesh = f"lambert{k}.obj"
        C.write_obj(os.path.join(run.workdir, mesh), v, f)
        C.write_field(os.path.join(run.workdir, f"lambert{k}.txt"), vals)
        r = run.run(mesh, ["-v", "0", "0", "1", "-r", str(res),
                           "-d", f"lambert{k}.txt", "-p", "double"])[0]
        exact_i = C.lambert_integral(1.0)
        rel = abs(r["integral"] - exact_i)/exact_i
        print(f"  {k:>6} {len(f):>7} {r['integral']:>12.7f} {exact_i:>12.7f} "
              f"{r['average']:>10.6f} {C.LAMBERT_MEAN:>10.6f} {rel:>10.2e}")
        tol = 0.05 if k <= 2 else (0.02 if k == 3 else 0.01)
        rep.check("Lambertian sphere: integral cos dA == (2/3) pi r^2",
                  f"icosphere({k}) integral", r["integral"], exact_i, tol)
        rep.check("Lambertian sphere: integral cos dA == (2/3) pi r^2",
                  f"icosphere({k}) mean", r["average"], C.LAMBERT_MEAN, tol)


def study_convergence(run, rep, quick):
    """Error against pixel size, to confirm the documented convergence behaviour.

    Coverage is decided by point-sampling pixel centres, so the error is
    dominated by the boundary band and should fall roughly linearly in the pixel
    size h. Reported as a fitted slope; the tolerance is loose because the error
    is not monotonic — a silhouette edge can land favourably on the pixel grid at
    one resolution and badly at the next.
    """
    v, f = C.unit_cube()
    mesh = "conv_cube.obj"
    C.write_obj(os.path.join(run.workdir, mesh), v, f)
    d = (0.3, 0.9, -0.31)
    exact = C.projected_area_convex(v, f, d)
    hs = [0.02, 0.01, 0.005, 0.0025] if quick else [0.02, 0.01, 0.005, 0.0025, 0.00125, 0.000625]

    print("\n  convergence of projected area with pixel size (unit cube, generic view)")
    print(f"  {'h':>10} {'area':>13} {'abs err':>12} {'err/h':>10}")
    pts = []
    for h in hs:
        r = run.run(mesh, ["-v", str(d[0]), str(d[1]), str(d[2]),
                           "-r", str(h), "-p", "double"])[0]
        err = abs(r["area"] - exact)
        print(f"  {h:>10.6f} {r['area']:>13.8f} {err:>12.3e} {err/h:>10.4f}")
        if err > 0:
            pts.append((math.log(h), math.log(err)))

    if len(pts) >= 3:
        n = len(pts)
        sx = sum(p[0] for p in pts); sy = sum(p[1] for p in pts)
        sxx = sum(p[0]*p[0] for p in pts); sxy = sum(p[0]*p[1] for p in pts)
        slope = (n*sxy - sx*sy)/(n*sxx - sx*sx)
        print(f"  fitted order: error ~ h^{slope:.2f}")
        rep.check("Convergence order in pixel size",
                  "fitted exponent (expect ~1, first order)", slope, 1.0, 0.5,
                  note="point sampling gives a boundary band of width ~h")

    # Whatever the order, the finest resolution must actually be accurate.
    r = run.run(mesh, ["-v", str(d[0]), str(d[1]), str(d[2]),
                       "-r", str(hs[-1]), "-p", "double"])[0]
    rep.check("Convergence order in pixel size",
              f"area at h={hs[-1]}", r["area"], exact, 0.002)


def study_invariances(run, rep, quick):
    """Properties that must hold identically, not merely approximately.

    These are the invariants whose violation produced silently wrong answers in
    earlier versions: a projection that depended on where the mesh sat in world
    space, and a cull that disagreed with the depth test.
    """
    v, f = C.unit_cube()
    d = (0.3, 0.9, -0.31)
    res = 0.002

    C.write_obj(os.path.join(run.workdir, "inv_origin.obj"), v, f)
    far = [(p[0]+5.0e6, p[1]+5.0e6, p[2]+5.0e6) for p in v]
    C.write_obj(os.path.join(run.workdir, "inv_far.obj"), far, f)

    for prec in ("float", "double"):
        a0 = run.run("inv_origin.obj", ["-v", *map(str, d), "-r", str(res), "-p", prec])[0]["area"]
        a1 = run.run("inv_far.obj",    ["-v", *map(str, d), "-r", str(res), "-p", prec])[0]["area"]
        rep.check("Invariances",
                  f"translation by 5e6 leaves the area unchanged ({prec})", a1, a0, 0.002)

    # A convex closed mesh must project identically whether or not backfaces are
    # culled, and the field statistics must agree too.
    vals = C.coordinate_field(v, 1)
    C.write_field(os.path.join(run.workdir, "inv_fy.txt"), vals)
    base = ["-v", "1", "0", "0", "-r", str(res), "-d", "inv_fy.txt", "-p", "double"]
    culled = run.run("inv_origin.obj", base)[0]
    nocull = run.run("inv_origin.obj", base + ["--no-cull"])[0]
    rep.check("Invariances", "culled vs --no-cull: area", culled["area"], nocull["area"], 1e-9)
    rep.check("Invariances", "culled vs --no-cull: mean", culled["average"], nocull["average"], 1e-9)

    # Scaling the mesh by s must scale the area by s^2 exactly.
    s = 3.0
    big = [(p[0]*s, p[1]*s, p[2]*s) for p in v]
    C.write_obj(os.path.join(run.workdir, "inv_big.obj"), big, f)
    a0 = run.run("inv_origin.obj", ["-v", *map(str, d), "-r", str(res), "-p", "double"])[0]["area"]
    a1 = run.run("inv_big.obj",    ["-v", *map(str, d), "-r", str(res*s), "-p", "double"])[0]["area"]
    rep.check("Invariances", "scaling by 3 scales the area by 9", a1, a0*s*s, 0.003)


# --------------------------------------------------------------------------

STUDIES = [
    ("convex",      study_convex_polyhedra),
    ("cube",        study_cube_closed_form),
    ("cauchy",      study_cauchy),
    ("sphere",      study_sphere),
    ("cylinder",    study_cylinder),
    ("fields",      study_linear_fields),
    ("lambert",     study_lambert),
    ("convergence", study_convergence),
    ("invariance",  study_invariances),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flatland", default=None, help="path to the flatland binary")
    ap.add_argument("--quick", action="store_true", help="coarse subset, for CI")
    ap.add_argument("--json", default=None, help="write the full result table here")
    ap.add_argument("--keep", default=None, help="keep generated meshes in this directory")
    ap.add_argument("--only", default=None,
                    help="comma-separated study names: " + ",".join(s for s, _ in STUDIES))
    args = ap.parse_args()

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    binary = args.flatland or os.path.join(root, "flatland")
    if not (os.path.isfile(binary) and os.access(binary, os.X_OK)):
        print(f"error: '{binary}' is not an executable. Build it with `make` first.",
              file=sys.stderr)
        return 2
    binary = os.path.abspath(binary)

    workdir = args.keep or tempfile.mkdtemp(prefix="flatland_validation_")
    os.makedirs(workdir, exist_ok=True)

    selected = set(args.only.split(",")) if args.only else None
    rep = Report()
    run = Runner(binary, workdir)

    print("FlatLand validation against closed-form results")
    print(f"  binary: {binary}")
    print(f"  mode:   {'quick' if args.quick else 'full'}")

    try:
        for name, fn in STUDIES:
            if selected and name not in selected:
                continue
            fn(run, rep, args.quick)
        for study in dict.fromkeys(r["study"] for r in rep.rows):
            rep.table(study)
    finally:
        if not args.keep:
            shutil.rmtree(workdir, ignore_errors=True)

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(rep.rows, fh, indent=2)
        print(f"\nwrote {args.json}")

    total = len(rep.rows)
    print(f"\n{total - rep.failures}/{total} checks within tolerance")
    if rep.failures:
        print(f"{rep.failures} FAILED")
    return 1 if rep.failures else 0


if __name__ == "__main__":
    sys.exit(main())
