# FlatLand for Python

A pure-`ctypes` binding to the FlatLand C ABI. No build step, no compiler, no
pybind11 or Cython — just `libflatland.so` (or `.dylib` / `.dll`) and the
standard library. NumPy is used if it is installed and is never required.

## Getting started

```sh
make lib                                   # builds libflatland.so at the repo root
PYTHONPATH=python python3 python/example.py
```

To use it from elsewhere, put `python/` on `PYTHONPATH` (or copy the `flatland/`
package next to your code) and make sure the shared library is findable — see
[Finding the library](#finding-the-library).

```python
import flatland

mesh = flatland.Mesh.load("part.obj")
r = mesh.project((1, 0, 0), resolution=0.001)
print(r.area)                              # projected (visible) area
```

## API

```python
flatland.__version__                       # "3.0.0", from fl_version_string()
flatland.version_info                      # (3, 0, 0)
flatland.HAS_NUMPY                         # True if the NumPy paths are active
flatland.DEFAULT_RESOLUTION                # 0.001, the library's own default

mesh = flatland.Mesh(vertices, faces)      # (N,3) floats, (M,3) ints; lists or ndarray
mesh = flatland.Mesh.load("part.obj")      # .obj / .stl, ASCII or binary
mesh.vertex_count, mesh.face_count
mesh.vertices, mesh.faces                  # copies back out, (N,3) / (M,3)

r = mesh.project((1, 0, 0), field=values, resolution=1e-3,
                 precision="float", cull=True, field_mode="auto")
r.area, r.covered_pixels, r.width, r.height
r.average, r.integral, r.min, r.max        # None unless r.has_stats
r.has_field, r.has_stats, r.as_dict()

results = mesh.project_batch(views,                    # (K,3)
                             field_matrix=matrix,      # (rows, steps)
                             field_columns=[0, 1, 2],  # or None -> column 0
                             resolutions=None,         # or K pixel sizes
                             threads=0)                # 0 = one per core

img = mesh.render((1, 0, 0), field=values)
img.width, img.height, img.mask, img.values, img.result
img.save_png("view.png")

fm = flatland.load_field("timeseries.txt", mesh)       # or mesh.load_field(...)
fm.shape, fm.rows, fm.cols, fm.field_mode, fm.column(0), fm.to_array()

flatland.angle_to_dir(azimuth_deg, elevation_deg)      # -> (x, y, z)
flatland.library_path()                                # which .so was loaded
```

A view direction is the direction the camera **looks along**, so the visible
surface is the one facing back toward `-direction`; magnitude is irrelevant.
Face indices are zero-based (OBJ files on disk are 1-based; the loader converts).
Fields carry one value per vertex (`field_mode="node"`) or per face
(`"face"`); `"auto"` picks from the length and is ambiguous only when a mesh has
as many faces as vertices.

## Statistics that do not exist

`average`, `integral`, `min` and `max` are `None` whenever `has_stats` is false —
that is, whenever the view carried no field or covered no pixels (geometry
thinner than one pixel, or a mesh seen exactly edge-on):

```python
r = mesh.project((0, 0, 1), field=values)
if r.has_stats:
    print(r.average)
else:
    print("nothing visible from here")     # r.average is None, not 0.0
```

The C ABI leaves those four fields at zero and sets `has_stats = 0`; passing the
zeros through as if they were measurements would fabricate data, so this binding
does not. `area` and `covered_pixels` are always meaningful.

## NumPy is optional

| | with NumPy | without |
|---|---|---|
| input arrays | `ndarray` accepted directly (contiguous `float64` passed without a copy) | lists, tuples, `array.array` |
| `mesh.vertices` / `.faces` | `(N, 3)` `ndarray` | list of tuples |
| `Image.mask` / `.values` | read-only zero-copy `ndarray` view | `array.array('B')` / `('d')` copy |

Both paths are exercised by the test suite on every run. Note that `Image.mask`
is `uint8` under NumPy, so count coverage with `mask.sum()` or
`np.count_nonzero(mask)` — the builtin `sum()` accumulates in `uint8` and wraps.

The raster buffers belong to the C library. Under NumPy they are exposed as
views rather than copies, and each view keeps its `Image` alive, so they stay
valid as long as they are reachable. The one way to invalidate one is to call
`Image.close()` (or leave its `with` block) while still holding it; copy with
`numpy.array(img.mask)` first if that can happen.

## Errors

Failures reported by the library raise `FlatlandError`, carrying `.message` (the
library's own `fl_last_error()` text), `.status` and `.status_name`. Subclasses
also inherit the matching built-in, so ordinary `except` clauses work:

| status | exception | also a |
|---|---|---|
| `FL_ERR_INVALID_ARGUMENT`, `FL_ERR_DIMENSION`, `FL_ERR_NUMERIC` | `FlatlandValueError` | `ValueError` |
| `FL_ERR_IO`, `FL_ERR_PARSE` | `FlatlandIOError` | `OSError` |
| `FL_ERR_OUT_OF_MEMORY` | `FlatlandMemoryError` | `MemoryError` |
| `FL_ERR_UNSUPPORTED`, `FL_ERR_INTERNAL` | `FlatlandError` | — |
| library not found at import | `FlatlandLibraryError` | `ImportError` |

`FL_ERR_PARSE` maps to the I/O exception because the library reports an
unopenable mesh file that way; "this file could not be read" is the distinction a
caller acts on.

Mistakes caught in Python before the ABI is reached — an unknown `precision`
name, a `field_columns` list that is not one entry per view, use of a closed
handle — raise plain `ValueError` or `TypeError`: no C call was made, so there is
no status code to report.

```python
try:
    r = mesh.project((0, 0, 0))
except flatland.FlatlandError as e:
    print(e.status_name, e.message)        # 'invalid argument' 'view_dir must be ...'
```

## Memory

`Mesh`, `Image` and `FieldMatrix` own C handles and free them in `close()`, on
`__del__`, and on leaving a `with` block:

```python
with flatland.Mesh(vertices, faces) as mesh:
    r = mesh.project((1, 0, 0), field=pressure)
print(r.integral)                          # Result holds only numbers; it outlives the mesh
```

A `Mesh` is immutable once created and safe to share across threads;
`project_batch` parallelizes internally and its results do not depend on
`threads`.

## Finding the library

Searched in order: `$FLATLAND_LIBRARY` (a file or a directory), the directory
holding the `flatland` package, its parents up to the repository root, a
`build/` directory there, then the dynamic loader's own path
(`LD_LIBRARY_PATH`, `DYLD_LIBRARY_PATH`, `PATH` on Windows, the system
directories). `libflatland.so`, `libflatland.dylib`, `flatland.dll` and
`libflatland.dll` are all recognised. If nothing is found, the import error says
how to build one and which paths were tried.

```sh
export FLATLAND_LIBRARY=/opt/flatland/lib/libflatland.so
```

## Tests

```sh
./tests/run_tests.sh ./flatland python      # the binding's suite alone
./tests/run_tests.sh ./flatland             # everything
PYTHONPATH=python python3 python/test_flatland.py    # the Python program directly
```

`tests/test_python.sh` runs `python/test_flatland.py` twice — once as the
environment is, and once with an import shim that makes `import numpy` fail — so
whichever path this machine would otherwise skip is still covered. It skips
cleanly if `python3` is not installed.
