# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
# for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

"""
flatland — Python bindings for the FlatLand projected-area engine.

FlatLand computes the projected (visible) surface area of a triangle mesh from
arbitrary view directions, and — given a scalar field defined per vertex or per
face — the mean, the extremes and the area integral over the visible projection.

This module is a pure-``ctypes`` wrapper over the C ABI declared in
``include/flatland.h``. It has no build step and no mandatory dependencies.
NumPy is used when it is importable (mesh arrays are passed without copying
where possible, and rasters come back as zero-copy views); everything works with
plain lists and tuples when it is not.

Quick start
-----------

    import flatland

    mesh = flatland.Mesh.load("part.obj")
    r = mesh.project((1, 0, 0), resolution=0.01)
    print(r.area)

    with flatland.Mesh(vertices, faces) as mesh:          # deterministic free
        r = mesh.project((1, 0, 0), field=pressure)
        if r.has_stats:                                   # see "Statistics" below
            print(r.average, r.integral)

Conventions (inherited from the C ABI, unchanged)
-------------------------------------------------
Coordinates   Vertices are (x, y, z); ``(N, 3)`` nested sequences, a flat
              sequence of ``3N`` numbers, or a NumPy array are all accepted.
Faces         Triangles are triples of ZERO-BASED vertex indices. (OBJ files on
              disk are 1-based; the loader converts.)
View          A view direction is the direction the camera LOOKS ALONG, so the
              visible surface is the one facing back toward ``-direction``.
              Magnitude is irrelevant; the vector is normalized internally.
Fields        A field has one value per vertex (node mode) or per face (face
              mode). Node fields are interpolated across each triangle; face
              fields are constant over it. ``field_mode="auto"`` infers the mode
              from the length, and is ambiguous — an error — only when the mesh
              has as many faces as vertices.

Statistics
----------
``Result.average``, ``.integral``, ``.min`` and ``.max`` are ``None`` whenever
``Result.has_stats`` is false, which happens when the view carried no field or
covered no pixels (geometry thinner than one pixel, or a mesh seen exactly
edge-on). They are not zero in that case: there is no measurement to report, and
returning 0.0 would fabricate one. ``.area`` and ``.covered_pixels`` are always
meaningful.

Errors
------
Every failure reported by the C library raises ``FlatlandError`` (or a subclass)
carrying ``.message`` — the library's own diagnostic — plus ``.status`` and
``.status_name``. The subclasses also inherit from a built-in exception so the
usual ``except`` clauses work:

    FL_ERR_INVALID_ARGUMENT, FL_ERR_DIMENSION, FL_ERR_NUMERIC
        -> FlatlandValueError(FlatlandError, ValueError)
    FL_ERR_IO, FL_ERR_PARSE
        -> FlatlandIOError(FlatlandError, OSError)
    FL_ERR_OUT_OF_MEMORY
        -> FlatlandMemoryError(FlatlandError, MemoryError)
    anything else (FL_ERR_UNSUPPORTED, FL_ERR_INTERNAL)
        -> FlatlandError

FL_ERR_PARSE maps to the I/O exception rather than the value one because the
library reports an unopenable mesh file that way; "this file could not be read"
is the useful distinction for a caller, not which half of the read failed.

Argument mistakes caught in Python before the library is called (an unknown
``precision`` name, a field_columns list of the wrong length, use of a closed
handle) raise plain ``ValueError`` or ``TypeError``: they never reached the C
ABI, so there is no status code to report.

Memory
------
``Mesh``, ``Image`` and ``FieldMatrix`` own C handles. Each frees its handle in
``close()``, on ``__del__``, and on leaving a ``with`` block. ``Result`` holds
only Python numbers, so results outlive the mesh that produced them.

``Image.mask`` and ``Image.values`` are zero-copy, read-only NumPy views of the
library's buffers when NumPy is present; the arrays keep the ``Image`` alive, so
they remain valid for as long as they are reachable. Calling ``Image.close()``
explicitly while such an array is still in use is the one way to invalidate one
— copy with ``numpy.array(img.mask)`` first if that is a possibility. Without
NumPy the same properties return ``array.array`` copies, which are never
affected.

Finding the library
-------------------
``libflatland.so`` / ``.dylib`` / ``flatland.dll`` is looked for in, in order:
``$FLATLAND_LIBRARY`` (a file or a directory), the directory holding this
package, its parents up to the repository root, a ``build/`` directory there,
then the loader's own search path (``LD_LIBRARY_PATH``, ``DYLD_LIBRARY_PATH``,
``PATH`` on Windows, the system directories).
"""

import array as _array
import ctypes as _ct
import ctypes.util as _ctutil
import os as _os
import sys as _sys

# NumPy is optional, always. Nothing below may assume it is present.
try:
    import numpy as _np

    HAS_NUMPY = True
except Exception:  # pragma: no cover - exercised by the no-numpy test run
    _np = None
    HAS_NUMPY = False

__all__ = [
    "Mesh",
    "Result",
    "Image",
    "FieldMatrix",
    "FlatlandError",
    "FlatlandValueError",
    "FlatlandIOError",
    "FlatlandMemoryError",
    "FlatlandLibraryError",
    "angle_to_dir",
    "load_field",
    "library_path",
    "HAS_NUMPY",
    "DEFAULT_RESOLUTION",
    "__version__",
    "version_info",
]


# --------------------------------------------------------------------------- #
# Status codes                                                                 #
# --------------------------------------------------------------------------- #

FL_OK = 0
FL_ERR_INVALID_ARGUMENT = 1
FL_ERR_OUT_OF_MEMORY = 2
FL_ERR_IO = 3
FL_ERR_PARSE = 4
FL_ERR_DIMENSION = 5
FL_ERR_NUMERIC = 6
FL_ERR_UNSUPPORTED = 7
FL_ERR_INTERNAL = 8


class FlatlandError(Exception):
    """A failure reported by the FlatLand C library.

    Attributes:
        message      the library's diagnostic, from fl_last_error()
        status       the numeric fl_status code, or None
        status_name  the human-readable status name, or None
    """

    message = ""
    status = None
    status_name = None


class FlatlandValueError(FlatlandError, ValueError):
    """A bad argument, a dimension mismatch or a non-finite value."""


class FlatlandIOError(FlatlandError, OSError):
    """A file could not be read or written, or its contents are malformed."""


class FlatlandMemoryError(FlatlandError, MemoryError):
    """The library could not allocate."""


class FlatlandLibraryError(FlatlandError, ImportError):
    """The shared library could not be located or loaded."""


_EXC_FOR_STATUS = {
    FL_ERR_INVALID_ARGUMENT: FlatlandValueError,
    FL_ERR_DIMENSION: FlatlandValueError,
    FL_ERR_NUMERIC: FlatlandValueError,
    FL_ERR_IO: FlatlandIOError,
    FL_ERR_PARSE: FlatlandIOError,
    FL_ERR_OUT_OF_MEMORY: FlatlandMemoryError,
}


# --------------------------------------------------------------------------- #
# Library discovery                                                            #
# --------------------------------------------------------------------------- #


def _library_names():
    if _sys.platform.startswith("win"):
        return ["flatland.dll", "libflatland.dll"]
    if _sys.platform == "darwin":
        return ["libflatland.dylib", "libflatland.so"]
    return ["libflatland.so"]


def _search_directories():
    """Directories to look in, nearest first, without duplicates."""
    here = _os.path.dirname(_os.path.abspath(__file__))
    root = _os.path.dirname(_os.path.dirname(here))  # <repo>/python/flatland -> <repo>
    dirs = [
        here,  # installed alongside the package
        _os.path.join(here, "lib"),
        _os.path.dirname(here),  # <repo>/python
        root,  # the repository root, where `make lib` writes
        _os.path.join(root, "build"),  # a CMake build tree
        _os.path.join(root, "build", "lib"),
        _os.path.join(root, "lib"),
    ]
    seen = set()
    out = []
    for d in dirs:
        d = _os.path.abspath(d)
        if d not in seen:
            seen.add(d)
            out.append(d)
    return out


def _load_library():
    """Return (CDLL, path_or_name). Raises FlatlandLibraryError if not found."""
    names = _library_names()
    tried = []

    override = _os.environ.get("FLATLAND_LIBRARY")
    if override:
        candidates = []
        if _os.path.isdir(override):
            candidates = [_os.path.join(override, n) for n in names]
        else:
            candidates = [override]
        for cand in candidates:
            try:
                return _ct.CDLL(cand), cand
            except OSError as exc:
                tried.append("%s (%s)" % (cand, exc))
        root = _os.path.dirname(
            _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
        )
        raise FlatlandLibraryError(
            "FLATLAND_LIBRARY is set to %r but no FlatLand library could be "
            "loaded from it.\n"
            "Build one with:\n"
            "    make -C %s lib\n"
            "then point FLATLAND_LIBRARY at the resulting %s, or unset it to "
            "search the usual places.\nTried:\n  %s"
            % (override, root, names[0], "\n  ".join(tried))
        )

    for d in _search_directories():
        for n in names:
            cand = _os.path.join(d, n)
            if _os.path.isfile(cand):
                try:
                    return _ct.CDLL(cand), cand
                except OSError as exc:
                    tried.append("%s (%s)" % (cand, exc))
            else:
                tried.append(cand)

    # ctypes.util.find_library consults the linker cache / DYLD / PATH.
    found = _ctutil.find_library("flatland")
    if found:
        try:
            return _ct.CDLL(found), found
        except OSError as exc:
            tried.append("%s (%s)" % (found, exc))

    # Bare names: the dynamic loader's own search path.
    for n in names:
        try:
            return _ct.CDLL(n), n
        except OSError as exc:
            tried.append("%s (%s)" % (n, exc))

    root = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
    raise FlatlandLibraryError(
        "the FlatLand shared library (%s) could not be found.\n"
        "Build it with:\n"
        "    make -C %s lib\n"
        "or point the binding at an existing one:\n"
        "    export FLATLAND_LIBRARY=/path/to/%s\n"
        "Searched:\n  %s" % (" / ".join(names), root, names[0], "\n  ".join(tried))
    )


_lib, _lib_path = _load_library()


def library_path():
    """Absolute path (or bare name) of the shared library actually loaded."""
    return _lib_path


# --------------------------------------------------------------------------- #
# C types                                                                      #
# --------------------------------------------------------------------------- #

_c_double_p = _ct.POINTER(_ct.c_double)
_c_float_p = _ct.POINTER(_ct.c_float)
_c_int32_p = _ct.POINTER(_ct.c_int32)
_c_uint8_p = _ct.POINTER(_ct.c_uint8)


class _FlMesh(_ct.Structure):
    """Opaque fl_mesh."""


class _FlImage(_ct.Structure):
    """Opaque fl_image."""


class _FlFieldMatrix(_ct.Structure):
    """Opaque fl_field_matrix."""


_mesh_p = _ct.POINTER(_FlMesh)
_image_p = _ct.POINTER(_FlImage)
_matrix_p = _ct.POINTER(_FlFieldMatrix)


class _FlOptions(_ct.Structure):
    _fields_ = [
        ("resolution", _ct.c_double),
        ("cull", _ct.c_int32),
        ("precision", _ct.c_int),  # enum fl_precision
        ("field_mode", _ct.c_int),  # enum fl_field_mode
        ("threads", _ct.c_int32),
        ("reserved", _ct.c_uint32 * 8),
    ]


class _FlResult(_ct.Structure):
    _fields_ = [
        ("area", _ct.c_double),
        ("average", _ct.c_double),
        ("integral", _ct.c_double),
        ("min", _ct.c_double),
        ("max", _ct.c_double),
        ("covered_pixels", _ct.c_int64),
        ("width", _ct.c_int32),
        ("height", _ct.c_int32),
        ("has_field", _ct.c_int32),
        ("has_stats", _ct.c_int32),
        ("reserved", _ct.c_uint32 * 8),
    ]


class _FlBatchDesc(_ct.Structure):
    _fields_ = [
        ("views", _c_double_p),
        ("view_count", _ct.c_size_t),
        ("field_matrix", _c_double_p),
        ("field_rows", _ct.c_size_t),
        ("field_cols", _ct.c_size_t),
        ("field_columns", _c_int32_p),
        ("resolutions", _c_double_p),
        ("reserved", _ct.c_uint32 * 8),
    ]


# Every binding below declares BOTH argtypes and restype. A missing argtype lets
# ctypes pass a pointer through its default int conversion, which truncates to
# 32 bits on LP64 platforms and corrupts the handle — silently, and only for
# addresses above 4 GiB.
_PROTOTYPES = (
    # name,                         restype,        argtypes
    ("fl_status_string", _ct.c_char_p, [_ct.c_int]),
    ("fl_last_error", _ct.c_char_p, []),
    ("fl_version", None, [_ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int)]),
    ("fl_version_string", _ct.c_char_p, []),
    ("fl_mesh_create", _ct.c_int, [_c_double_p, _ct.c_size_t, _c_int32_p, _ct.c_size_t, _ct.POINTER(_mesh_p)]),
    ("fl_mesh_create_f32", _ct.c_int, [_c_float_p, _ct.c_size_t, _c_int32_p, _ct.c_size_t, _ct.POINTER(_mesh_p)]),
    ("fl_mesh_load", _ct.c_int, [_ct.c_char_p, _ct.POINTER(_mesh_p)]),
    ("fl_mesh_destroy", None, [_mesh_p]),
    ("fl_mesh_vertex_count", _ct.c_size_t, [_mesh_p]),
    ("fl_mesh_face_count", _ct.c_size_t, [_mesh_p]),
    ("fl_mesh_copy_vertices", _ct.c_int, [_mesh_p, _c_double_p, _ct.c_size_t]),
    ("fl_mesh_copy_faces", _ct.c_int, [_mesh_p, _c_int32_p, _ct.c_size_t]),
    ("fl_options_init", None, [_ct.POINTER(_FlOptions)]),
    ("fl_project", _ct.c_int, [_mesh_p, _c_double_p, _c_double_p, _ct.c_size_t, _ct.POINTER(_FlOptions), _ct.POINTER(_FlResult)]),
    ("fl_batch_desc_init", None, [_ct.POINTER(_FlBatchDesc)]),
    ("fl_project_batch", _ct.c_int, [_mesh_p, _ct.POINTER(_FlBatchDesc), _ct.POINTER(_FlOptions), _ct.POINTER(_FlResult)]),
    ("fl_angle_to_dir", None, [_ct.c_double, _ct.c_double, _c_double_p]),
    ("fl_project_image", _ct.c_int, [_mesh_p, _c_double_p, _c_double_p, _ct.c_size_t, _ct.POINTER(_FlOptions), _ct.POINTER(_image_p), _ct.POINTER(_FlResult)]),
    ("fl_image_width", _ct.c_int32, [_image_p]),
    ("fl_image_height", _ct.c_int32, [_image_p]),
    ("fl_image_mask", _c_uint8_p, [_image_p]),
    ("fl_image_values", _c_double_p, [_image_p]),
    ("fl_image_write_ppm", _ct.c_int, [_image_p, _ct.c_char_p]),
    ("fl_image_destroy", None, [_image_p]),
    ("fl_field_matrix_load", _ct.c_int, [_ct.c_char_p, _mesh_p, _ct.c_int, _ct.POINTER(_matrix_p)]),
    ("fl_field_matrix_rows", _ct.c_size_t, [_matrix_p]),
    ("fl_field_matrix_cols", _ct.c_size_t, [_matrix_p]),
    ("fl_field_matrix_mode", _ct.c_int, [_matrix_p]),
    ("fl_field_matrix_data", _c_double_p, [_matrix_p]),
    ("fl_field_matrix_copy_column", _ct.c_int, [_matrix_p, _ct.c_size_t, _c_double_p, _ct.c_size_t]),
    ("fl_field_matrix_destroy", None, [_matrix_p]),
)


def _bind():
    missing = []
    for name, restype, argtypes in _PROTOTYPES:
        try:
            fn = getattr(_lib, name)
        except AttributeError:
            missing.append(name)
            continue
        fn.restype = restype
        fn.argtypes = argtypes
    if missing:
        raise FlatlandLibraryError(
            "%s does not export the expected FlatLand C ABI; missing: %s. "
            "Rebuild it from this source tree with `make lib`."
            % (_lib_path, ", ".join(missing))
        )


_bind()


# --------------------------------------------------------------------------- #
# Errors                                                                       #
# --------------------------------------------------------------------------- #


def _status_name(status):
    raw = _lib.fl_status_string(int(status))
    return raw.decode("utf-8", "replace") if raw else "status %d" % status


def _check(status):
    """Raise the mapped exception unless status is FL_OK.

    fl_last_error() is read first: its buffer is only guaranteed until the next
    FlatLand call on this thread, and ctypes' c_char_p restype copies the bytes
    out at call time, so the message is safe once it is a Python object.
    """
    status = int(status)
    if status == FL_OK:
        return
    raw = _lib.fl_last_error()
    message = raw.decode("utf-8", "replace") if raw else ""
    name = _status_name(status)
    exc_type = _EXC_FOR_STATUS.get(status, FlatlandError)
    text = "%s: %s" % (name, message) if message else name
    exc = exc_type(text)  # one argument only: OSError parses 2+ args as errno/strerror
    exc.message = message
    exc.status = status
    exc.status_name = name
    raise exc


# --------------------------------------------------------------------------- #
# Version                                                                      #
# --------------------------------------------------------------------------- #


def _version_info():
    maj, mnr, pat = _ct.c_int(0), _ct.c_int(0), _ct.c_int(0)
    _lib.fl_version(_ct.byref(maj), _ct.byref(mnr), _ct.byref(pat))
    return (maj.value, mnr.value, pat.value)


version_info = _version_info()
__version__ = (_lib.fl_version_string() or b"").decode("utf-8", "replace")


# --------------------------------------------------------------------------- #
# Enumerations                                                                 #
# --------------------------------------------------------------------------- #

PRECISION_FLOAT = 0
PRECISION_DOUBLE = 1

FIELD_AUTO = 0
FIELD_NODE = 1
FIELD_FACE = 2

_PRECISION_CODES = {"float": 0, "single": 0, "f32": 0, "double": 1, "f64": 1}
_FIELD_MODE_CODES = {"auto": 0, "node": 1, "vertex": 1, "point": 1, "face": 2, "cell": 2}
_FIELD_MODE_NAMES = {0: "auto", 1: "node", 2: "face"}


def _precision_code(precision):
    if isinstance(precision, bool):
        raise TypeError("precision must be 'float' or 'double', not a bool")
    if isinstance(precision, int):
        if precision in (PRECISION_FLOAT, PRECISION_DOUBLE):
            return precision
        raise ValueError("unknown precision %r; use 'float' or 'double'" % (precision,))
    try:
        return _PRECISION_CODES[str(precision).strip().lower()]
    except KeyError:
        raise ValueError("unknown precision %r; use 'float' or 'double'" % (precision,))


def _field_mode_code(field_mode):
    if isinstance(field_mode, bool):
        raise TypeError("field_mode must be 'auto', 'node' or 'face', not a bool")
    if isinstance(field_mode, int):
        if field_mode in (FIELD_AUTO, FIELD_NODE, FIELD_FACE):
            return field_mode
        raise ValueError(
            "unknown field_mode %r; use 'auto', 'node' or 'face'" % (field_mode,)
        )
    try:
        return _FIELD_MODE_CODES[str(field_mode).strip().lower()]
    except KeyError:
        raise ValueError(
            "unknown field_mode %r; use 'auto', 'node' or 'face'" % (field_mode,)
        )


# --------------------------------------------------------------------------- #
# Input conversion                                                             #
# --------------------------------------------------------------------------- #
#
# Each helper returns a ctypes-compatible buffer that the caller must keep
# referenced for the duration of the C call. Every array FlatLand receives is
# borrowed for the call only, so a temporary is fine as long as it outlives it.


def _is_ndarray(obj):
    return HAS_NUMPY and isinstance(obj, _np.ndarray)


def _doubles(values, name):
    """Flat C array of doubles. Returns (buffer, count)."""
    if HAS_NUMPY:
        try:
            arr = _np.ascontiguousarray(values, dtype=_np.float64)
        except ValueError as exc:
            # A ragged nested sequence lands here; it is a bad value, not a bad
            # type, and must stay a ValueError for the caller.
            raise ValueError("%s must be a rectangular block of numbers (%s)" % (name, exc))
        except TypeError as exc:
            raise TypeError("%s must be a sequence of numbers (%s)" % (name, exc))
        arr = arr.reshape(-1)
        return arr, int(arr.size)
    try:
        flat = [float(v) for v in values]
    except TypeError as exc:
        raise TypeError("%s must be a sequence of numbers (%s)" % (name, exc))
    buf = (_ct.c_double * len(flat))(*flat)
    return buf, len(flat)


def _doubles_ptr(buf):
    """POINTER(c_double) for a buffer returned by _doubles/_rows/_matrix."""
    if buf is None:
        return None
    if _is_ndarray(buf):
        if buf.size == 0:
            return None
        return buf.ctypes.data_as(_c_double_p)
    return _ct.cast(buf, _c_double_p)


def _int32_ptr(buf):
    if buf is None:
        return None
    if _is_ndarray(buf):
        if buf.size == 0:
            return None
        return buf.ctypes.data_as(_c_int32_p)
    return _ct.cast(buf, _c_int32_p)


def _rows(values, ncols, name):
    """An (N, ncols) table. Accepts nested sequences, a flat sequence of N*ncols
    numbers, or a NumPy array of either shape. Returns (buffer, nrows)."""
    if HAS_NUMPY:
        try:
            arr = _np.ascontiguousarray(values, dtype=_np.float64)
        except ValueError as exc:
            raise ValueError("%s must be a rectangular block of numbers (%s)" % (name, exc))
        except TypeError as exc:
            raise TypeError("%s must be numeric (%s)" % (name, exc))
        if arr.ndim == 2:
            if arr.shape[1] != ncols:
                raise ValueError(
                    "%s must have %d columns, got shape %r" % (name, ncols, arr.shape)
                )
            nrows = int(arr.shape[0])
        elif arr.ndim == 1:
            if arr.size % ncols:
                raise ValueError(
                    "%s has %d values, which is not a multiple of %d"
                    % (name, arr.size, ncols)
                )
            nrows = int(arr.size // ncols)
        else:
            raise ValueError("%s must be 1- or 2-dimensional, got %dD" % (name, arr.ndim))
        if nrows == 0:
            raise ValueError("%s is empty" % name)
        return arr.reshape(-1), nrows

    items = list(values)
    if not items:
        raise ValueError("%s is empty" % name)
    flat = []
    if isinstance(items[0], (list, tuple)) or (
        hasattr(items[0], "__len__") and not isinstance(items[0], (str, bytes))
    ):
        for i, row in enumerate(items):
            row = list(row)
            if len(row) != ncols:
                raise ValueError(
                    "%s[%d] has %d values, expected %d" % (name, i, len(row), ncols)
                )
            flat.extend(float(v) for v in row)
        nrows = len(items)
    else:
        if len(items) % ncols:
            raise ValueError(
                "%s has %d values, which is not a multiple of %d"
                % (name, len(items), ncols)
            )
        flat = [float(v) for v in items]
        nrows = len(items) // ncols
    buf = (_ct.c_double * len(flat))(*flat)
    return buf, nrows


def _int32s(values, name, count=None):
    """Flat C array of int32. Returns (buffer, length)."""
    if HAS_NUMPY and not isinstance(values, (list, tuple)):
        try:
            arr = _np.ascontiguousarray(values).reshape(-1)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("%s must be integers that fit in int32 (%s)" % (name, exc))
        if arr.size:
            # astype(int32) would wrap silently; check the range first. A NaN
            # fails both comparisons and is rejected here too.
            lo, hi = arr.min(), arr.max()
            if not (lo >= -2147483648 and hi <= 2147483647):
                raise ValueError(
                    "%s contains values that do not fit in int32 (%r .. %r)"
                    % (name, lo, hi)
                )
        arr = _np.ascontiguousarray(arr.astype(_np.int32, copy=False))
        n = int(arr.size)
    else:
        flat = []
        for v in values:
            try:
                iv = int(v)
            except (TypeError, ValueError) as exc:
                raise ValueError("%s must contain integers (%s)" % (name, exc))
            if not (-2147483648 <= iv <= 2147483647):
                raise ValueError("%s contains %d, which does not fit in int32" % (name, iv))
            flat.append(iv)
        arr = (_ct.c_int32 * len(flat))(*flat)
        n = len(flat)
    if count is not None and n != count:
        raise ValueError("%s has %d entries, expected %d" % (name, n, count))
    return arr, n


def _int32_rows(values, ncols, name):
    """An (N, ncols) integer table, same input shapes as _rows."""
    if HAS_NUMPY and not isinstance(values, (list, tuple)):
        try:
            arr = _np.ascontiguousarray(values)
        except ValueError as exc:
            raise ValueError("%s must be a rectangular block of integers (%s)" % (name, exc))
        except TypeError as exc:
            raise TypeError("%s must be numeric (%s)" % (name, exc))
        if arr.ndim == 2:
            if arr.shape[1] != ncols:
                raise ValueError(
                    "%s must have %d columns, got shape %r" % (name, ncols, arr.shape)
                )
            nrows = int(arr.shape[0])
        elif arr.ndim == 1:
            if arr.size % ncols:
                raise ValueError(
                    "%s has %d values, which is not a multiple of %d"
                    % (name, arr.size, ncols)
                )
            nrows = int(arr.size // ncols)
        else:
            raise ValueError("%s must be 1- or 2-dimensional, got %dD" % (name, arr.ndim))
        if nrows == 0:
            raise ValueError("%s is empty" % name)
        buf, _ = _int32s(arr.reshape(-1), name)
        return buf, nrows

    items = list(values)
    if not items:
        raise ValueError("%s is empty" % name)
    flat = []
    if isinstance(items[0], (list, tuple)) or (
        hasattr(items[0], "__len__") and not isinstance(items[0], (str, bytes))
    ):
        for i, row in enumerate(items):
            row = list(row)
            if len(row) != ncols:
                raise ValueError(
                    "%s[%d] has %d values, expected %d" % (name, i, len(row), ncols)
                )
            flat.extend(row)
        nrows = len(items)
    else:
        if len(items) % ncols:
            raise ValueError(
                "%s has %d values, which is not a multiple of %d"
                % (name, len(items), ncols)
            )
        flat = items
        nrows = len(items) // ncols
    buf, _ = _int32s(flat, name)
    return buf, nrows


def _matrix(values, name):
    """A row-major (rows, cols) matrix. A flat/1-D input is one column.
    Returns (buffer, rows, cols)."""
    if HAS_NUMPY:
        try:
            arr = _np.ascontiguousarray(values, dtype=_np.float64)
        except ValueError as exc:
            raise ValueError("%s must be a rectangular block of numbers (%s)" % (name, exc))
        except TypeError as exc:
            raise TypeError("%s must be numeric (%s)" % (name, exc))
        if arr.ndim == 1:
            rows, cols = int(arr.size), 1
        elif arr.ndim == 2:
            rows, cols = int(arr.shape[0]), int(arr.shape[1])
        else:
            raise ValueError("%s must be 1- or 2-dimensional, got %dD" % (name, arr.ndim))
        if rows == 0 or cols == 0:
            raise ValueError("%s is empty" % name)
        return arr.reshape(-1), rows, cols

    items = list(values)
    if not items:
        raise ValueError("%s is empty" % name)
    if isinstance(items[0], (list, tuple)) or (
        hasattr(items[0], "__len__") and not isinstance(items[0], (str, bytes))
    ):
        rows_list = [list(r) for r in items]
        cols = len(rows_list[0])
        if cols == 0:
            raise ValueError("%s has empty rows" % name)
        flat = []
        for i, row in enumerate(rows_list):
            if len(row) != cols:
                raise ValueError(
                    "%s is ragged: row 0 has %d values, row %d has %d"
                    % (name, cols, i, len(row))
                )
            flat.extend(float(v) for v in row)
        rows = len(rows_list)
    else:
        flat = [float(v) for v in items]
        rows, cols = len(flat), 1
    buf = (_ct.c_double * len(flat))(*flat)
    return buf, rows, cols


def _field_arg(field):
    """Prepare an optional field argument.

    Returns (buffer, pointer, length). The caller MUST hold on to the buffer
    until the C call returns — the pointer alone does not necessarily own it.
    """
    if field is None:
        return None, None, 0
    buf, n = _doubles(field, "field")
    if n == 0:
        raise ValueError("field is empty; pass field=None for a geometry-only run")
    return buf, _doubles_ptr(buf), n


def _fspath(path, name="path"):
    if isinstance(path, bytes):
        return path
    try:
        path = _os.fspath(path)
    except TypeError:
        raise TypeError("%s must be a filesystem path, got %r" % (name, type(path).__name__))
    if isinstance(path, bytes):
        return path
    return path.encode(_sys.getfilesystemencoding() or "utf-8", "surrogateescape")


def _build_options(resolution, precision, cull, field_mode, threads):
    o = _FlOptions()
    _lib.fl_options_init(_ct.byref(o))  # zeroes the reserved words; do not skip
    if resolution is not None:
        try:
            o.resolution = float(resolution)
        except (TypeError, ValueError):
            raise ValueError("resolution must be a number, got %r" % (resolution,))
    o.cull = 1 if cull else 0
    o.precision = _precision_code(precision)
    o.field_mode = _field_mode_code(field_mode)
    try:
        o.threads = int(threads)
    except (TypeError, ValueError):
        raise ValueError("threads must be an integer, got %r" % (threads,))
    return o


def _default_resolution():
    o = _FlOptions()
    _lib.fl_options_init(_ct.byref(o))
    return o.resolution


DEFAULT_RESOLUTION = _default_resolution()


# --------------------------------------------------------------------------- #
# Results                                                                      #
# --------------------------------------------------------------------------- #


class Result(object):
    """What one projected view measured.

    ``average``, ``integral``, ``min`` and ``max`` are ``None`` unless
    ``has_stats`` is true — that is, unless a field was supplied AND the view
    covered at least one pixel. They are never reported as 0.0 in that case: an
    unmeasured statistic is not a measurement of zero.
    """

    __slots__ = (
        "area",
        "average",
        "integral",
        "min",
        "max",
        "covered_pixels",
        "width",
        "height",
        "has_field",
        "has_stats",
    )

    def __init__(self, raw):
        self.area = float(raw.area)
        self.covered_pixels = int(raw.covered_pixels)
        self.width = int(raw.width)
        self.height = int(raw.height)
        self.has_field = bool(raw.has_field)
        self.has_stats = bool(raw.has_stats)
        if self.has_stats:
            self.average = float(raw.average)
            self.integral = float(raw.integral)
            self.min = float(raw.min)
            self.max = float(raw.max)
        else:
            self.average = None
            self.integral = None
            self.min = None
            self.max = None

    def as_dict(self):
        """A plain dict, suitable for json.dumps (unmeasured stats are null)."""
        return dict(
            area=self.area,
            average=self.average,
            integral=self.integral,
            min=self.min,
            max=self.max,
            covered_pixels=self.covered_pixels,
            width=self.width,
            height=self.height,
            has_field=self.has_field,
            has_stats=self.has_stats,
        )

    def __repr__(self):
        if self.has_stats:
            return (
                "Result(area=%.6g, average=%.6g, integral=%.6g, min=%.6g, "
                "max=%.6g, covered_pixels=%d, %dx%d)"
                % (
                    self.area,
                    self.average,
                    self.integral,
                    self.min,
                    self.max,
                    self.covered_pixels,
                    self.width,
                    self.height,
                )
            )
        return "Result(area=%.6g, covered_pixels=%d, %dx%d, has_stats=False)" % (
            self.area,
            self.covered_pixels,
            self.width,
            self.height,
        )


# --------------------------------------------------------------------------- #
# Handles                                                                      #
# --------------------------------------------------------------------------- #


class _Handle(object):
    """Common ownership machinery: close(), `with`, and a safe __del__.

    The destroy function is stashed on the instance, because module globals can
    already be None by the time __del__ runs during interpreter shutdown.
    """

    __slots__ = ("_ptr", "_destroy", "__weakref__")

    def _adopt(self, ptr, destroy):
        self._ptr = ptr
        self._destroy = destroy

    @property
    def closed(self):
        return not getattr(self, "_ptr", None)

    def _check_open(self):
        ptr = getattr(self, "_ptr", None)
        if not ptr:
            raise ValueError(
                "this %s has been closed; its C handle is gone"
                % type(self).__name__.lower()
            )
        return ptr

    def close(self):
        """Release the C handle now. Idempotent."""
        ptr = getattr(self, "_ptr", None)
        self._ptr = None
        if ptr:
            self._destroy(ptr)

    def __enter__(self):
        self._check_open()
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    def __del__(self):
        try:
            self.close()
        except Exception:  # pragma: no cover - shutdown paths
            pass


class Image(object):
    """A rendered view: the coverage mask and, where a field was supplied, the
    interpolated per-pixel values.

    Row 0 is the BOTTOM row in mesh space and pixels are indexed
    ``[y * width + x]``, matching the C ABI. ``mask`` and ``values`` are flat
    sequences of ``width * height`` elements; with NumPy they are read-only
    zero-copy views of the library's buffers that keep this Image alive, so they
    stay valid as long as they are reachable. Calling ``close()`` explicitly
    invalidates them — copy first if you need them afterwards.

    A view that covers no pixels is a valid 0x0 image with an empty mask.
    """

    __slots__ = ("_ptr", "_destroy", "_width", "_height", "_mask", "_values", "result", "__weakref__")

    def __init__(self, ptr, result):
        self._ptr = ptr
        self._destroy = _lib.fl_image_destroy
        self._width = int(_lib.fl_image_width(ptr))
        self._height = int(_lib.fl_image_height(ptr))
        self._mask = None
        self._values = None
        self.result = result

    # -- ownership ---------------------------------------------------------
    @property
    def closed(self):
        return not self._ptr

    def _check_open(self):
        if not self._ptr:
            raise ValueError("this Image has been closed; its C buffers are gone")
        return self._ptr

    def close(self):
        """Release the raster now. Idempotent.

        Any NumPy view previously returned by ``mask``/``values`` refers to
        freed memory afterwards; take a copy before closing if that matters.
        """
        ptr = self._ptr
        self._ptr = None
        self._mask = None
        self._values = None
        if ptr:
            self._destroy(ptr)

    def __enter__(self):
        self._check_open()
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    def __del__(self):
        try:
            self.close()
        except Exception:  # pragma: no cover
            pass

    # -- geometry ----------------------------------------------------------
    @property
    def width(self):
        return self._width

    @property
    def height(self):
        return self._height

    @property
    def size(self):
        """(width, height)."""
        return (self._width, self._height)

    @property
    def has_field(self):
        """True if this raster carries per-pixel field values."""
        self._check_open()
        return bool(_lib.fl_image_values(self._ptr))

    # -- buffers -----------------------------------------------------------
    def _wrap(self, ptr, ctype, np_dtype, array_code):
        """Zero-copy view (NumPy) or a copy (array.array) of a borrowed buffer."""
        n = self._width * self._height
        if n == 0 or not ptr:
            if HAS_NUMPY:
                empty = _np.empty(0, dtype=np_dtype)
                empty.setflags(write=False)
                return empty
            return _array.array(array_code)
        address = _ct.cast(ptr, _ct.c_void_p).value
        if not HAS_NUMPY:
            out = _array.array(array_code)
            out.frombytes(_ct.string_at(address, n * _ct.sizeof(ctype)))
            return out
        buf = (ctype * n).from_address(address)
        # The C buffer belongs to the fl_image. Anchor the Image to the ctypes
        # object so a live array can never outlive the memory it points at:
        # numpy's frombuffer keeps `buf` as the array's base, and `buf` keeps
        # `self` alive through this attribute.
        buf._flatland_image = self
        arr = _np.frombuffer(buf, dtype=np_dtype, count=n)
        arr.setflags(write=False)  # the C ABI hands these out as const
        return arr

    @property
    def mask(self):
        """Coverage, 1 where a pixel centre fell on a visible triangle.

        ``len(mask) == width * height``; ``sum(mask)`` equals the result's
        ``covered_pixels``.
        """
        self._check_open()
        if self._mask is None:
            self._mask = self._wrap(_lib.fl_image_mask(self._ptr), _ct.c_uint8, "uint8", "B")
        return self._mask

    @property
    def values(self):
        """Per-pixel field values, or None when the view carried no field (or
        covered no pixels). Values at uncovered pixels are not meaningful; use
        ``mask`` to select."""
        self._check_open()
        if self._values is None:
            ptr = _lib.fl_image_values(self._ptr)
            if not ptr:
                return None
            self._values = self._wrap(ptr, _ct.c_double, "float64", "d")
        return self._values

    def save_ppm(self, path):
        """Write a false-color binary PPM (P6).

        Covered pixels are ramped across the field range; a fieldless image is
        written as a white silhouette. An empty (0x0) view has nothing to write
        and raises FlatlandValueError.
        """
        ptr = self._check_open()
        _check(_lib.fl_image_write_ppm(ptr, _fspath(path)))
        return path

    def __repr__(self):
        if not self._ptr:
            return "<Image closed>"
        return "<Image %dx%d%s>" % (
            self._width,
            self._height,
            " with field" if self.has_field else "",
        )


class FieldMatrix(_Handle):
    """A field time series read from disk: one row per mesh entity, one column
    per timestep. Pass it straight to ``Mesh.project_batch(field_matrix=...)``.
    """

    __slots__ = ("_rows", "_cols", "_mode")

    def __init__(self, path, mesh, field_mode="auto"):
        if not isinstance(mesh, Mesh):
            raise TypeError("mesh must be a flatland.Mesh")
        out = _matrix_p()
        _check(
            _lib.fl_field_matrix_load(
                _fspath(path), mesh._check_open(), _field_mode_code(field_mode), _ct.byref(out)
            )
        )
        self._adopt(out, _lib.fl_field_matrix_destroy)
        self._rows = int(_lib.fl_field_matrix_rows(out))
        self._cols = int(_lib.fl_field_matrix_cols(out))
        self._mode = int(_lib.fl_field_matrix_mode(out))

    @property
    def rows(self):
        return self._rows

    @property
    def cols(self):
        return self._cols

    @property
    def shape(self):
        return (self._rows, self._cols)

    @property
    def field_mode(self):
        """'node' or 'face' — which entity the rows correspond to."""
        return _FIELD_MODE_NAMES.get(self._mode, "auto")

    def column(self, index):
        """One timestep as a NumPy array (or a list), length ``rows``."""
        ptr = self._check_open()
        n = self._rows
        buf = (_ct.c_double * n)()
        _check(_lib.fl_field_matrix_copy_column(ptr, int(index), buf, n))
        if HAS_NUMPY:
            return _np.frombuffer(bytes(buf), dtype=_np.float64).copy()
        return list(buf)

    def to_array(self):
        """The whole matrix as a (rows, cols) NumPy array, or a list of lists."""
        ptr = self._check_open()
        data = _lib.fl_field_matrix_data(ptr)
        n = self._rows * self._cols
        if not data or n == 0:
            return _np.empty((self._rows, self._cols)) if HAS_NUMPY else []
        raw = _ct.string_at(_ct.cast(data, _ct.c_void_p), n * _ct.sizeof(_ct.c_double))
        if HAS_NUMPY:
            return _np.frombuffer(raw, dtype=_np.float64).copy().reshape(self._rows, self._cols)
        flat = _array.array("d")
        flat.frombytes(raw)
        return [list(flat[r * self._cols : (r + 1) * self._cols]) for r in range(self._rows)]

    def __len__(self):
        return self._rows

    def __repr__(self):
        if self.closed:
            return "<FieldMatrix closed>"
        return "<FieldMatrix %dx%d (%s)>" % (self._rows, self._cols, self.field_mode)


class Mesh(_Handle):
    """An immutable triangle mesh.

    ``Mesh(vertices, faces)`` copies caller arrays into the library;
    ``Mesh.load(path)`` reads an OBJ or STL file. The handle is freed by
    ``close()``, by leaving a ``with`` block, or when the object is collected.

    A Mesh is immutable once created and safe to share across threads;
    ``project_batch`` parallelizes internally.
    """

    __slots__ = ("_vertex_count", "_face_count", "_source")

    def __init__(self, vertices, faces):
        vbuf, nv = _rows(vertices, 3, "vertices")
        fbuf, nf = _int32_rows(faces, 3, "faces")
        out = _mesh_p()
        _check(
            _lib.fl_mesh_create(
                _doubles_ptr(vbuf), nv, _int32_ptr(fbuf), nf, _ct.byref(out)
            )
        )
        self._adopt(out, _lib.fl_mesh_destroy)
        self._vertex_count = int(_lib.fl_mesh_vertex_count(out))
        self._face_count = int(_lib.fl_mesh_face_count(out))
        self._source = None

    @classmethod
    def load(cls, path):
        """Read a mesh from an OBJ or STL file (binary or ASCII, auto-detected).

        Polygonal OBJ faces are fan-triangulated. STL carries no shared
        vertices, so an STL mesh has three vertices per triangle and its natural
        field mode is 'face'.
        """
        self = cls.__new__(cls)
        out = _mesh_p()
        _check(_lib.fl_mesh_load(_fspath(path), _ct.byref(out)))
        self._adopt(out, _lib.fl_mesh_destroy)
        self._vertex_count = int(_lib.fl_mesh_vertex_count(out))
        self._face_count = int(_lib.fl_mesh_face_count(out))
        self._source = _os.fspath(path) if not isinstance(path, bytes) else path
        return self

    # -- geometry ----------------------------------------------------------
    @property
    def vertex_count(self):
        return self._vertex_count

    @property
    def face_count(self):
        return self._face_count

    @property
    def vertices(self):
        """The vertices, in the frame they were given in: an (N, 3) NumPy array,
        or a list of (x, y, z) tuples. A fresh copy on every access."""
        ptr = self._check_open()
        n = self._vertex_count * 3
        buf = (_ct.c_double * n)()
        _check(_lib.fl_mesh_copy_vertices(ptr, buf, n))
        if HAS_NUMPY:
            return _np.frombuffer(bytes(buf), dtype=_np.float64).copy().reshape(-1, 3)
        return [(buf[i], buf[i + 1], buf[i + 2]) for i in range(0, n, 3)]

    @property
    def faces(self):
        """The triangles as zero-based vertex indices: an (M, 3) NumPy array, or
        a list of (a, b, c) tuples. A fresh copy on every access."""
        ptr = self._check_open()
        n = self._face_count * 3
        buf = (_ct.c_int32 * n)()
        _check(_lib.fl_mesh_copy_faces(ptr, buf, n))
        if HAS_NUMPY:
            return _np.frombuffer(bytes(buf), dtype=_np.int32).copy().reshape(-1, 3)
        return [(buf[i], buf[i + 1], buf[i + 2]) for i in range(0, n, 3)]

    # -- projection --------------------------------------------------------
    def project(
        self,
        view,
        field=None,
        resolution=None,
        precision="float",
        cull=True,
        field_mode="auto",
    ):
        """Project one view and return a Result.

        view        the direction the camera looks along, 3 numbers; the visible
                    surface is the one facing back toward -view. Need not be a
                    unit vector, but must be non-zero and finite.
        field       one value per vertex (node) or per face (face), or None for
                    a geometry-only run
        resolution  pixel edge length in mesh units (default: DEFAULT_RESOLUTION,
                    the library's own default of 0.001)
        precision   'float' (default, faster) or 'double'
        cull        True (default) to backface-cull, False to render all faces
        field_mode  'auto' (default), 'node' or 'face'
        """
        ptr = self._check_open()
        vbuf, nv = _rows(view, 3, "view")
        if nv != 1:
            raise ValueError("view must be a single (x, y, z) direction, got %d" % nv)
        fbuf, fptr, flen = _field_arg(field)
        opts = _build_options(resolution, precision, cull, field_mode, 0)
        raw = _FlResult()
        _check(
            _lib.fl_project(
                ptr, _doubles_ptr(vbuf), fptr, flen, _ct.byref(opts), _ct.byref(raw)
            )
        )
        del fbuf, vbuf  # borrowed for the call only; explicit about the lifetime
        return Result(raw)

    def project_batch(
        self,
        views,
        field_matrix=None,
        field_columns=None,
        resolutions=None,
        threads=0,
        resolution=None,
        precision="float",
        cull=True,
        field_mode="auto",
    ):
        """Project many views over one mesh in parallel; returns a list of Result
        in view order, whatever order the workers finished in.

        views          (K, 3) directions
        field_matrix   optional (rows, cols) matrix, one ROW per mesh entity and
                       one COLUMN per timestep, so a whole time series is one
                       array. A FieldMatrix, a NumPy array, a list of rows, or a
                       flat sequence (taken as a single column) all work.
        field_columns  which column each view reads, K entries; None means every
                       view reads column 0
        resolutions    per-view pixel size, K entries; None means every view uses
                       ``resolution``
        threads        worker threads; 0 (default) means one per core. Results do
                       not depend on this.

        If any view fails, the whole call fails and reports the lowest-indexed
        failure, so the error is the same run to run.
        """
        ptr = self._check_open()
        vbuf, nviews = _rows(views, 3, "views")

        desc = _FlBatchDesc()
        _lib.fl_batch_desc_init(_ct.byref(desc))
        desc.views = _doubles_ptr(vbuf)
        desc.view_count = nviews

        # Keep every temporary referenced until after the call returns.
        mbuf = cbuf = rbuf = None
        if field_matrix is not None:
            if isinstance(field_matrix, FieldMatrix):
                mptr = _lib.fl_field_matrix_data(field_matrix._check_open())
                desc.field_matrix = mptr
                desc.field_rows = field_matrix.rows
                desc.field_cols = field_matrix.cols
                mbuf = field_matrix  # hold the handle open for the call
            else:
                mbuf, rows, cols = _matrix(field_matrix, "field_matrix")
                desc.field_matrix = _doubles_ptr(mbuf)
                desc.field_rows = rows
                desc.field_cols = cols
        elif field_columns is not None:
            raise ValueError("field_columns was given without a field_matrix")

        if field_columns is not None:
            # The library reads exactly view_count entries; a short list would be
            # an out-of-bounds read, so the length is checked here.
            cbuf, _ = _int32s(field_columns, "field_columns", count=nviews)
            desc.field_columns = _int32_ptr(cbuf)

        if resolutions is not None:
            rbuf, nres = _doubles(resolutions, "resolutions")
            if nres != nviews:
                raise ValueError(
                    "resolutions has %d entries, expected one per view (%d)"
                    % (nres, nviews)
                )
            desc.resolutions = _doubles_ptr(rbuf)

        opts = _build_options(resolution, precision, cull, field_mode, threads)
        out = (_FlResult * nviews)()
        _check(_lib.fl_project_batch(ptr, _ct.byref(desc), _ct.byref(opts), out))
        results = [Result(out[i]) for i in range(nviews)]
        del mbuf, cbuf, rbuf, vbuf
        return results

    def render(
        self,
        view,
        field=None,
        resolution=None,
        precision="float",
        cull=True,
        field_mode="auto",
    ):
        """Project one view and keep the raster. Returns an Image; its
        ``.result`` is the same Result ``project()`` would have returned."""
        ptr = self._check_open()
        vbuf, nv = _rows(view, 3, "view")
        if nv != 1:
            raise ValueError("view must be a single (x, y, z) direction, got %d" % nv)
        fbuf, fptr, flen = _field_arg(field)
        opts = _build_options(resolution, precision, cull, field_mode, 0)
        img = _image_p()
        raw = _FlResult()
        _check(
            _lib.fl_project_image(
                ptr,
                _doubles_ptr(vbuf),
                fptr,
                flen,
                _ct.byref(opts),
                _ct.byref(img),
                _ct.byref(raw),
            )
        )
        del fbuf, vbuf
        return Image(img, Result(raw))

    def load_field(self, path, field_mode="auto"):
        """Read a whitespace-separated field matrix sized against this mesh."""
        return FieldMatrix(path, self, field_mode)

    def __repr__(self):
        if self.closed:
            return "<Mesh closed>"
        src = " from %r" % self._source if self._source else ""
        return "<Mesh %d vertices, %d faces%s>" % (
            self._vertex_count,
            self._face_count,
            src,
        )


# --------------------------------------------------------------------------- #
# Free functions                                                               #
# --------------------------------------------------------------------------- #


def angle_to_dir(azimuth_deg, elevation_deg):
    """Unit view direction for an azimuth/elevation pair, in degrees.

    Azimuth sweeps around +Z measured from +X; elevation rises from the XY plane
    toward +Z. Returns a 3-tuple, which every function here accepts.
    """
    out = (_ct.c_double * 3)()
    _lib.fl_angle_to_dir(float(azimuth_deg), float(elevation_deg), out)
    return (out[0], out[1], out[2])


def load_field(path, mesh, field_mode="auto"):
    """Read a field matrix from disk, sized against ``mesh``. See FieldMatrix."""
    return FieldMatrix(path, mesh, field_mode)
