#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
# for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

"""
test_flatland.py — exercises the ctypes binding the way a real consumer would.

Driven by tests/test_python.sh, which runs it twice: once as the environment is
(with NumPy if it is installed) and once with NumPy forced unavailable, because
the optional-dependency path is only genuinely covered if both are executed.

Prints "  PASS <what>" / "  FAIL <what> (<detail>)" lines for the shell driver
to replay into the suite counters, and exits non-zero if anything failed.

Usage: test_flatland.py [repo-root] [--expect-numpy yes|no]
"""

import gc
import math
import os
import shutil
import sys
import tempfile

passed = 0
failed = 0


def ok(what):
    global passed
    print("  PASS %s" % what)
    passed += 1


def bad(what, detail="condition false"):
    global failed
    print("  FAIL %s (%s)" % (what, detail))
    failed += 1


def sect(name):
    print("\n-- %s --" % name)


def check(cond, what):
    if cond:
        ok(what)
    else:
        bad(what)


def near(actual, expected, tol, what):
    if actual is None:
        bad(what, "got None, want %.9g" % expected)
    elif abs(actual - expected) <= tol:
        ok("%s (%.9g ~= %.9g)" % (what, actual, expected))
    else:
        bad(what, "got %.9g, want %.9g" % (actual, expected))


def equal(actual, expected, what):
    if actual == expected:
        ok("%s (%r)" % (what, actual))
    else:
        bad(what, "got %r, want %r" % (actual, expected))


def raises(exc_type, fn, what, also=None):
    """Assert fn() raises exc_type. `also` is an extra type it must be too."""
    try:
        fn()
    except exc_type as exc:
        if also is not None and not isinstance(exc, also):
            bad(what, "raised %s, which is not a %s" % (type(exc).__name__, also.__name__))
            return
        detail = str(exc)
        if not detail:
            bad(what, "raised %s with an empty message" % type(exc).__name__)
            return
        ok("%s [%s: %.70s]" % (what, type(exc).__name__, detail))
    except Exception as exc:  # noqa: BLE001 - reporting the wrong type is the point
        bad(what, "raised %s: %.60s" % (type(exc).__name__, exc))
    else:
        bad(what, "did not raise")


def covered(mask):
    """Number of covered pixels, for either flavour of mask.

    NumPy's uint8 arrays wrap around under the builtin sum(), so the ndarray
    case has to accumulate in a wider type — mask.sum() already does.
    """
    try:
        return int(mask.sum(dtype="int64"))
    except (AttributeError, TypeError):
        return int(sum(mask))


def flatten(rows):
    """Flat list of floats from an (N, k) NumPy array or a list of tuples."""
    out = []
    for row in rows:
        out.extend(float(v) for v in row)
    return out


# --------------------------------------------------------------------------- #
# Fixtures                                                                     #
# --------------------------------------------------------------------------- #

# A closed unit box spanning [0,1]^3, CCW outward normals. 8 vertices and 12
# faces, so node and face fields are never ambiguous.
BOX_V = [
    (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
    (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
]
BOX_F = [
    (0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7),
    (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
    (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
]

RES = 0.002  # fine enough that the unit box lands within 1% of area 1


def main():
    root = None
    expect_numpy = None
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--expect-numpy":
            expect_numpy = args[i + 1]
            i += 2
        else:
            root = args[i]
            i += 1
    if root is None:
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    import flatland
    from flatland import (
        FlatlandError,
        FlatlandIOError,
        FlatlandValueError,
        Mesh,
    )

    np = None
    if flatland.HAS_NUMPY:
        import numpy as np  # noqa: F811

    print("  (binding: %s, numpy: %s)" % (flatland.library_path(), flatland.HAS_NUMPY))

    tmp = tempfile.mkdtemp(prefix="flatland_py_")

    # ---------------------------------------------------------------- setup --
    sect("module, version and bindings")

    if expect_numpy is not None:
        check(
            flatland.HAS_NUMPY == (expect_numpy == "yes"),
            "the run uses the intended numpy state (expected %s)" % expect_numpy,
        )

    check(bool(flatland.__version__), "__version__ is non-empty")
    check(
        isinstance(flatland.version_info, tuple)
        and len(flatland.version_info) == 3
        and all(isinstance(v, int) for v in flatland.version_info),
        "version_info is a 3-tuple of ints",
    )
    check(
        flatland.__version__.startswith(
            "%d.%d.%d" % flatland.version_info
        ),
        "__version__ agrees with version_info",
    )
    # fl_options_init's defaults double as a check that the fl_options struct
    # layout declared here matches the library's.
    near(flatland.DEFAULT_RESOLUTION, 0.001, 0.0, "the default resolution is 0.001")

    # The classic ctypes bug is a missing argtype: the pointer is then passed
    # through the default int conversion and truncated to 32 bits.
    unset = [
        name
        for name, restype, argtypes in flatland._PROTOTYPES
        if getattr(flatland._lib, name).argtypes is None
        or list(getattr(flatland._lib, name).argtypes) != list(argtypes)
        or getattr(flatland._lib, name).restype is not restype
    ]
    check(not unset, "every bound function declares argtypes and restype")
    check(len(flatland._PROTOTYPES) >= 30, "the whole C ABI is bound (%d functions)"
          % len(flatland._PROTOTYPES))

    # ----------------------------------------------------------------- mesh --
    sect("mesh from Python sequences")

    box = Mesh(BOX_V, BOX_F)
    equal(box.vertex_count, 8, "vertex_count")
    equal(box.face_count, 12, "face_count")
    equal(len(box.vertices), 8, "len(mesh.vertices)")
    equal(len(box.faces), 12, "len(mesh.faces)")
    check(flatten(box.vertices) == flatten(BOX_V), "vertices round-trip in the original frame")
    check(
        [int(v) for v in flatten(box.faces)] == [int(v) for v in flatten(BOX_F)],
        "faces round-trip",
    )

    # A flat sequence is accepted as well as an (N, 3) one.
    flat_box = Mesh(flatten(BOX_V), [int(v) for v in flatten(BOX_F)])
    equal(flat_box.vertex_count, 8, "a flat vertex sequence is accepted")
    near(flat_box.project((1, 0, 0), resolution=RES).area, 1.0, 0.01,
         "...and projects identically")
    flat_box.close()

    # ----------------------------------------------------------- projection --
    sect("projection")

    r = box.project((1, 0, 0), resolution=RES, precision="double")
    near(r.area, 1.0, 0.01, "unit box projects to area 1")
    check(r.has_field is False and r.has_stats is False, "no field means no statistics")
    check(
        r.average is None and r.integral is None and r.min is None and r.max is None,
        "the four statistics are None without a field",
    )
    check(r.covered_pixels > 0 and r.width > 0 and r.height > 0,
          "raster dimensions are reported")
    near(r.area, r.covered_pixels * RES * RES, 1e-9, "area == covered_pixels * resolution^2")

    # Field = vertex x. Looking along +X the visible surface is x = 0, so the
    # mean must be 0. This is the assertion that catches an inverted backface
    # cull, and it needs a field that VARIES across the mesh.
    fx = [v[0] for v in BOX_V]
    r = box.project((1, 0, 0), field=fx, resolution=RES, precision="double")
    check(r.has_field is True and r.has_stats is True, "a covered field view has statistics")
    near(r.average, 0.0, 1e-6, "culled view resolves to the NEAR surface (average 0)")
    near(r.min, 0.0, 1e-6, "field min over the near surface")
    near(r.max, 0.0, 1e-6, "field max over the near surface")
    near(r.integral, 0.0, 1e-6, "integral over a zero-valued surface is 0")

    r_nocull = box.project((1, 0, 0), field=fx, resolution=RES, precision="double", cull=False)
    near(r_nocull.average, 0.0, 1e-6, "culling agrees with no culling")

    # The opposite view sees x = 1: a constant-field mesh could not tell these
    # two apart, which is exactly why the field varies.
    r_back = box.project((-1, 0, 0), field=fx, resolution=RES, precision="double")
    near(r_back.average, 1.0, 1e-6, "the reversed view resolves to the far face (average 1)")

    # A constant face field: the area integral is the constant times the area.
    ff = [4.0] * 12
    r = box.project((0, 0, 1), field=ff, resolution=RES, precision="double")
    near(r.average, 4.0, 1e-9, "constant face field average")
    near(r.integral, 4.0 * r.area, 1e-9, "integral == constant * area")
    near(r.min, 4.0, 1e-9, "constant face field min")
    near(r.max, 4.0, 1e-9, "constant face field max")

    # field_mode can be forced; 8 values is a node field, 12 a face field.
    near(box.project((0, 0, 1), field=ff, field_mode="face", resolution=RES).average,
         4.0, 1e-6, "field_mode='face' is honoured")
    near(box.project((1, 0, 0), field=fx, field_mode="node", resolution=RES).average,
         0.0, 1e-6, "field_mode='node' is honoured")

    a_f = box.project((1, 1, 1), resolution=RES, precision="float").area
    a_d = box.project((1, 1, 1), resolution=RES, precision="double").area
    near(a_f, a_d, 1e-3, "float and double agree on area")

    # ---------------------------------------------------------------- batch --
    sect("batch")

    views = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    matrix = [[10.0, 20.0, 30.0] for _ in range(8)]  # 8 rows (node) x 3 columns
    out = box.project_batch(views, field_matrix=matrix, field_columns=[0, 1, 2],
                            resolutions=[0.004] * 3)
    equal(len(out), 3, "project_batch returns one result per view")
    near(out[0].average, 10.0, 1e-9, "view 0 reads column 0")
    near(out[1].average, 20.0, 1e-9, "view 1 reads column 1")
    near(out[2].average, 30.0, 1e-9, "view 2 reads column 2")
    near(out[0].area, 1.0, 0.02, "batch view 0 area")

    single = box.project_batch(views, field_matrix=matrix, field_columns=[0, 1, 2],
                               resolutions=[0.004] * 3, threads=1)
    many = box.project_batch(views, field_matrix=matrix, field_columns=[0, 1, 2],
                             resolutions=[0.004] * 3, threads=8)
    check(
        [x.average for x in single] == [x.average for x in many]
        and [x.area for x in single] == [x.area for x in many]
        and [x.covered_pixels for x in single] == [x.covered_pixels for x in many],
        "batch results are identical at threads=1 and threads=8",
    )
    check([x.average for x in many] == [x.average for x in out],
          "...and identical to the default thread count")

    no_cols = box.project_batch(views, field_matrix=matrix, resolutions=[0.004] * 3)
    near(no_cols[2].average, 10.0, 1e-9, "field_columns=None means every view reads column 0")

    geom = box.project_batch(views, resolution=0.004)
    check(
        all(x.has_field is False and x.has_stats is False for x in geom)
        and all(x.average is None for x in geom),
        "a geometry-only batch has no statistics",
    )
    near(geom[1].area, 1.0, 0.02, "geometry-only batch still measures area")

    per_view = box.project_batch([(1, 0, 0), (1, 0, 0)], resolutions=[0.01, 0.005])
    check(per_view[1].width > per_view[0].width,
          "per-view resolutions are honoured (%d px vs %d px)"
          % (per_view[0].width, per_view[1].width))

    # A single view is a legitimate batch.
    equal(len(box.project_batch([(1, 0, 0)], resolution=0.01)), 1, "a one-view batch works")

    # -------------------------------------------------- zero-coverage stats --
    sect("zero coverage")

    # A sliver thinner than one pixel: covered_pixels is 0, so the field
    # statistics are not measurements and must come back as None, never 0.0.
    sliver = Mesh([(0, 0, 0), (1, 0, 0), (1, 0.000001, 0)], [(0, 1, 2)])
    z = sliver.project((0, 0, -1), field=[-5.0, -3.0, -1.0], resolution=0.1, cull=False)
    equal(z.covered_pixels, 0, "the sliver covers no pixels")
    check(z.has_field is True, "has_field is set")
    check(z.has_stats is False, "...but has_stats is not")
    check(z.average is None, "average is None, not 0.0")
    check(z.integral is None, "integral is None, not 0.0")
    check(z.min is None, "min is None, not 0.0")
    check(z.max is None, "max is None, not 0.0")
    near(z.area, 0.0, 0.0, "area is a genuine 0")
    check(z.as_dict()["average"] is None, "as_dict() reports the missing average as None")

    zb = sliver.project_batch([(0, 0, -1)], field_matrix=[[-5.0], [-3.0], [-1.0]],
                              resolution=0.1, cull=False)
    check(zb[0].has_stats is False and zb[0].average is None,
          "a zero-coverage batch view has no statistics either")
    sliver.close()

    # --------------------------------------------------------------- raster --
    sect("rasters")

    fy = [v[1] for v in BOX_V]
    img = box.render((1, 0, 0), field=fy, resolution=0.01)
    equal((img.width, img.height), (img.result.width, img.result.height),
          "image dimensions match the result")
    equal(len(img.mask), img.width * img.height, "len(mask) == width * height")
    equal(covered(img.mask), img.result.covered_pixels, "mask agrees with covered_pixels")
    check(img.values is not None, "values are exposed for a field view")
    equal(len(img.values), img.width * img.height, "len(values) == width * height")
    check(img.has_field is True, "has_field is true for a field raster")

    ppm = os.path.join(tmp, "render.ppm")
    img.save_ppm(ppm)
    with open(ppm, "rb") as fp:
        header = fp.read(2)
    equal(header, b"P6", "save_ppm writes a P6 file")
    check(os.path.getsize(ppm) > 2, "...with a body")

    plain = box.render((0, 0, 1), resolution=0.01)
    check(plain.values is None, "a fieldless raster has no values")
    equal(len(plain.mask), plain.width * plain.height, "a fieldless mask is still sized")
    plain.save_ppm(os.path.join(tmp, "silhouette.ppm"))
    ok("a fieldless raster writes a silhouette PPM")
    plain.close()

    # A view that covers nothing must still yield a valid, empty image.
    tri = Mesh([(0, 0, 0), (1, 0, 0), (1, 1, 0)], [(0, 1, 2)])
    empty = tri.render((1, 0, 0), resolution=RES)
    equal((empty.width, empty.height), (0, 0), "an empty view yields a 0x0 image")
    equal(len(empty.mask), 0, "an empty image has an empty mask")
    check(empty.result.covered_pixels == 0 and empty.result.has_stats is False,
          "an empty view reports no coverage and no statistics")
    raises(FlatlandError, lambda: empty.save_ppm(os.path.join(tmp, "nope.ppm")),
           "writing a PPM for an empty view is refused")
    check(not os.path.exists(os.path.join(tmp, "nope.ppm")), "...and writes no file")
    empty.close()
    tri.close()

    # -------------------------------------------------------------- numpy-ish --
    sect("numpy interop" if flatland.HAS_NUMPY else "list/array fallback")

    if flatland.HAS_NUMPY:
        nv = np.array(BOX_V, dtype=float)
        nf = np.array(BOX_F, dtype=np.int32)
        nmesh = Mesh(nv, nf)
        equal(nmesh.vertex_count, 8, "a mesh builds from NumPy arrays")
        check(isinstance(nmesh.vertices, np.ndarray) and nmesh.vertices.shape == (8, 3),
              "mesh.vertices is an (N, 3) ndarray")
        check(isinstance(nmesh.faces, np.ndarray) and nmesh.faces.shape == (12, 3),
              "mesh.faces is an (M, 3) ndarray")
        check(bool(np.allclose(nmesh.vertices, nv)), "ndarray vertices round-trip")
        near(nmesh.project(np.array([1.0, 0.0, 0.0]), field=nv[:, 0],
                           resolution=RES, precision="double").average,
             0.0, 1e-6, "an ndarray field gives the same near-surface average")
        # A non-contiguous column view must be copied, not passed as-is.
        strided = nv[:, 0:3:2][:, 0]
        check(not strided.flags["C_CONTIGUOUS"], "the strided field really is non-contiguous")
        near(nmesh.project((1, 0, 0), field=strided, resolution=RES, precision="double").average,
             0.0, 1e-6, "a non-contiguous ndarray field is handled")
        nout = nmesh.project_batch(np.array(views, dtype=float),
                                   field_matrix=np.tile([10.0, 20.0, 30.0], (8, 1)),
                                   field_columns=np.array([0, 1, 2], dtype=np.int32),
                                   resolutions=np.full(3, 0.004))
        check([x.average for x in nout] == [10.0, 20.0, 30.0],
              "an ndarray batch reads its per-view columns")
        # Fortran order is not row-major; it must be repacked, not misread.
        fort = nmesh.project_batch(views, field_matrix=np.asfortranarray(
            np.tile([10.0, 20.0, 30.0], (8, 1))), field_columns=[0, 1, 2],
            resolutions=[0.004] * 3)
        check([x.average for x in fort] == [10.0, 20.0, 30.0],
              "a Fortran-ordered field matrix is repacked correctly")

        nimg = nmesh.render((1, 0, 0), field=nv[:, 1], resolution=0.01)
        check(isinstance(nimg.mask, np.ndarray) and nimg.mask.dtype == np.uint8,
              "Image.mask is a uint8 ndarray")
        check(isinstance(nimg.values, np.ndarray) and nimg.values.dtype == np.float64,
              "Image.values is a float64 ndarray")
        check(nimg.mask.flags["OWNDATA"] is False and nimg.values.flags["OWNDATA"] is False,
              "the raster arrays are zero-copy views, not copies")
        check(nimg.mask.flags["WRITEABLE"] is False,
              "the borrowed raster arrays are read-only")
        check(nimg.mask is nimg.mask, "repeated access returns the same array")
        ncov = covered(nimg.mask)
        equal(ncov, nimg.result.covered_pixels, "the ndarray mask agrees with covered_pixels")
        # Dropping the Image must not free the buffer under a live array.
        held_mask, held_values = nimg.mask, nimg.values
        del nimg
        gc.collect()
        equal(covered(held_mask), ncov,
              "a mask array keeps its Image alive after the Image is dropped")
        check(math.isfinite(float(held_values[len(held_values) // 2])),
              "...and the values array is still readable")
        del held_mask, held_values
        nmesh.close()
    else:
        import array as pyarray

        check(isinstance(box.vertices, list) and isinstance(box.vertices[0], tuple),
              "mesh.vertices is a list of tuples without NumPy")
        check(isinstance(box.faces, list) and isinstance(box.faces[0], tuple),
              "mesh.faces is a list of tuples without NumPy")
        check(all(isinstance(v, int) for v in box.faces[0]),
              "face indices come back as ints")
        m2 = box.render((1, 0, 0), field=fy, resolution=0.01)
        check(isinstance(m2.mask, pyarray.array) and m2.mask.typecode == "B",
              "Image.mask is an array('B') without NumPy")
        check(isinstance(m2.values, pyarray.array) and m2.values.typecode == "d",
              "Image.values is an array('d') without NumPy")
        equal(len(m2.mask), m2.width * m2.height, "the fallback mask is sized correctly")
        # The fallback is a copy, so it is unaffected by closing the image.
        held = m2.mask
        expected_cov = m2.result.covered_pixels
        m2.close()
        equal(covered(held), expected_cov, "the fallback mask survives closing the image")
        raises(ValueError, lambda: m2.mask, "using a closed Image raises")

    # ------------------------------------------------------------- lifetime --
    sect("handle lifetime")

    cube_path = os.path.join(root, "examples", "cube_area", "cube.obj")
    loaded = Mesh.load(cube_path)
    check(loaded.vertex_count > 0 and loaded.face_count > 0, "Mesh.load reads an OBJ")
    lr = loaded.project((1, 0, 0), resolution=0.005)
    near(lr.area, 1.0, 0.02, "the cube example projects to area 1")
    near(loaded.project(flatland.angle_to_dir(90, 0), resolution=0.005).area, 1.0, 0.02,
         "...and from +Y via angle_to_dir")

    limg = loaded.render((1, 0, 0), field=[1.0] * loaded.vertex_count, resolution=0.02)
    lcov = limg.result.covered_pixels
    loaded.close()
    gc.collect()
    # The result and the raster are independent of the mesh that made them.
    near(lr.area, 1.0, 0.02, "a result outlives the mesh it came from")
    equal(covered(limg.mask), lcov, "a raster outlives the mesh it came from")
    limg.save_ppm(os.path.join(tmp, "detached.ppm"))
    ok("...and can still be written to disk")
    limg.close()
    raises(ValueError, lambda: loaded.project((1, 0, 0)), "a closed mesh refuses to project")
    loaded.close()
    ok("close() is idempotent")

    with Mesh(BOX_V, BOX_F) as ctx:
        near(ctx.project((1, 0, 0), resolution=0.01).area, 1.0, 0.01,
             "a mesh works inside a with-block")
    check(ctx.closed, "leaving the with-block closes the mesh")
    raises(ValueError, lambda: ctx.project((1, 0, 0)), "...and using it afterwards raises")

    with Mesh(BOX_V, BOX_F).render((1, 0, 0), resolution=0.01) as wimg:
        check(wimg.width > 0, "an Image works as a context manager")
    check(wimg.closed, "leaving the with-block closes the image")

    # Dropping every Python reference must free the handle without a crash.
    for _ in range(200):
        Mesh(BOX_V, BOX_F).project((1, 0, 0), resolution=0.05)
    gc.collect()
    ok("200 unreferenced meshes are created and collected without a crash")

    # ---------------------------------------------------------- field files --
    sect("field matrices from disk")

    field_path = os.path.join(tmp, "field.txt")
    with open(field_path, "w") as fp:
        fp.write("# one row per vertex, two timesteps\n")
        for _ in range(8):
            fp.write("2 5\n")
    fm = flatland.load_field(field_path, box)
    equal(fm.shape, (8, 2), "load_field reads a rows x cols matrix")
    equal(fm.field_mode, "node", "the matrix knows it is a node field")
    equal(len(fm.column(0)), 8, "a column has one value per row")
    near(float(fm.column(1)[0]), 5.0, 1e-12, "the second column reads back")
    fmb = box.project_batch([(1, 0, 0), (0, 1, 0)], field_matrix=fm, field_columns=[0, 1],
                            resolution=0.01)
    check([round(x.average, 9) for x in fmb] == [2.0, 5.0],
          "a FieldMatrix feeds project_batch directly")
    raises(FlatlandError, lambda: fm.column(99), "an out-of-range column is refused")
    fm.close()
    raises(ValueError, lambda: fm.column(0), "a closed FieldMatrix refuses to be read")

    # --------------------------------------------------------------- errors --
    sect("errors are raised, not crashes")

    raises(FlatlandValueError, lambda: box.project((1, 0, 0), field=[0.0] * 7),
           "a wrong-length field is rejected", also=ValueError)
    raises(FlatlandValueError, lambda: box.project((0, 0, 0)),
           "a zero view vector is rejected", also=ValueError)
    raises(FlatlandValueError, lambda: box.project((float("inf"), 0, 1)),
           "an infinite view direction is rejected", also=ValueError)
    raises(FlatlandValueError, lambda: box.project((float("nan"), 0, 1)),
           "a NaN view direction is rejected", also=ValueError)
    raises(FlatlandValueError, lambda: Mesh(BOX_V, [(0, 1, 99)]),
           "an out-of-range face index is rejected", also=ValueError)
    raises(FlatlandIOError, lambda: Mesh.load(os.path.join(tmp, "no_such_file.obj")),
           "a missing mesh file is rejected", also=OSError)
    raises(FlatlandValueError, lambda: box.project((1, 0, 0), resolution=0),
           "a zero resolution is rejected", also=ValueError)
    raises(FlatlandValueError, lambda: box.project((1, 0, 0), resolution=-1),
           "a negative resolution is rejected", also=ValueError)
    raises(FlatlandValueError, lambda: box.project((1, 0, 0), resolution=1e-10),
           "an impossibly fine resolution is refused", also=ValueError)
    raises(FlatlandValueError, lambda: box.project((1, 0, 0), resolution=float("nan")),
           "a non-finite resolution is rejected", also=ValueError)
    raises(FlatlandValueError,
           lambda: box.project_batch([(1, 0, 0)], field_matrix=[[1.0]] * 8, field_columns=[9]),
           "an out-of-range field column is rejected", also=ValueError)
    raises(FlatlandValueError, lambda: box.project_batch([(1, 0, 0)], threads=-1),
           "a negative thread count is rejected", also=ValueError)
    raises(FlatlandValueError,
           lambda: box.project_batch([(1, 0, 0)], field_matrix=[[1.0]] * 7),
           "a field matrix with the wrong row count is rejected", also=ValueError)
    raises(FlatlandIOError, lambda: flatland.load_field("/no/such/field.txt", box),
           "a missing field file is rejected", also=OSError)

    # Mistakes caught in Python, before the ABI: plain built-in exceptions.
    raises(ValueError, lambda: box.project((1, 0, 0), precision="quad"),
           "an unknown precision name is rejected")
    raises(ValueError, lambda: box.project((1, 0, 0), field=[1.0] * 8, field_mode="cells"),
           "an unknown field_mode name is rejected")
    raises(ValueError, lambda: box.project((1, 0)),
           "a two-component view is rejected")
    raises(ValueError, lambda: Mesh([], []), "an empty mesh is rejected")
    raises(ValueError, lambda: Mesh([(0, 0), (1, 0, 0), (0, 1, 0)], [(0, 1, 2)]),
           "a ragged vertex list is rejected")
    raises(ValueError,
           lambda: box.project_batch([(1, 0, 0), (0, 1, 0)], field_matrix=[[1.0]] * 8,
                                     field_columns=[0]),
           "field_columns shorter than the view list is rejected")
    raises(ValueError,
           lambda: box.project_batch([(1, 0, 0), (0, 1, 0)], resolutions=[0.01]),
           "resolutions shorter than the view list is rejected")
    raises(ValueError, lambda: box.project_batch([(1, 0, 0)], field_columns=[0]),
           "field_columns without a field_matrix is rejected")
    raises(ValueError, lambda: Mesh(BOX_V, [(0, 1, 2 ** 40)]),
           "a face index too large for int32 is rejected")

    # An error must carry the library's own diagnostic, not a generic string.
    try:
        box.project((1, 0, 0), field=[0.0] * 7)
    except FlatlandError as exc:
        check(bool(exc.message) and bool(exc.status_name) and exc.status == 5,
              "FlatlandError carries fl_last_error() text and the status name")
        check("field" in exc.message, "...and the message names the offending input")

    # A failure must not poison the next call.
    near(box.project((1, 0, 0), resolution=RES).area, 1.0, 0.01,
         "a mesh still works after a failed call")

    # ---------------------------------------------------------------- misc ---
    sect("angles and conversions")

    for (az, el, want, label) in (
        (0.0, 0.0, (1, 0, 0), "+X"),
        (90.0, 0.0, (0, 1, 0), "+Y"),
        (0.0, 90.0, (0, 0, 1), "+Z"),
        (180.0, 0.0, (-1, 0, 0), "-X"),
    ):
        d = flatland.angle_to_dir(az, el)
        check(len(d) == 3 and all(abs(a - b) < 1e-9 for a, b in zip(d, want)),
              "angle_to_dir(%g, %g) == %s" % (az, el, label))
    near(box.project(flatland.angle_to_dir(0, 0), resolution=0.01).area,
         box.project((1, 0, 0), resolution=0.01).area, 1e-12,
         "an angle view matches the equivalent vector view")

    check("Result(" in repr(box.project((1, 0, 0), resolution=0.05)), "Result has a repr")
    check("Mesh" in repr(box), "Mesh has a repr")

    box.close()
    shutil.rmtree(tmp, ignore_errors=True)

    print("\nPython binding: %d passed, %d failed" % (passed, failed))
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:  # a crash in the harness itself must still be visible
        import traceback

        traceback.print_exc()
        print("  FAIL the python test program itself raised")
        sys.exit(1)
