# FlatLand for MATLAB

A MATLAB binding for FlatLand, built on `loadlibrary`/`calllib` over the C ABI
in `include/flatland.h`. No MEX file, no compilation of binding code — the
shared library you already build with `make lib` is loaded directly.

```matlab
fl = Flatland(vertices, faces);
r  = fl.project([1 0 0], 'Field', values, 'Resolution', 1e-3);
fprintf('area %.4f, mean field %.4f\n', r.area, r.average);
delete(fl);
```

---

## Prerequisites

| | |
|---|---|
| **MATLAB** | R2016b or newer (the binding uses `isstring`, `strlength`, `strjoin`) |
| **A C compiler configured for MATLAB** | required — see below |
| **`libflatland.so` / `.dylib` / `.dll`** | build with `make lib` in the repository root |

### The compiler prerequisite

This is the one setup step people miss. `loadlibrary` does not read a header
directly: it **preprocesses it with a C compiler** and then parses the result
with its own restricted C parser. MATLAB must therefore have a compiler
configured, even though you are not compiling any binding code yourself.

```matlab
>> mex -setup C
```

| Platform | What to install |
|---|---|
| Linux | `gcc` (`apt install build-essential`) |
| macOS | Xcode Command Line Tools (`xcode-select --install`) |
| Windows | MinGW-w64 — free, from **Home → Add-Ons → Get Add-Ons**, search "MinGW" — or MSVC |

`flatland_load` checks for a configured compiler first and raises
`Flatland:NoCompiler` with these instructions, rather than letting you hit a
parser error forty lines deep.

---

## Setup

```bash
cd /path/to/flatland
make lib                 # produces libflatland.so at the repository root
```

```matlab
>> cd /path/to/flatland/matlab
>> addpath(pwd)          % optional, to use it from elsewhere
>> flatland_test         % self-check
>> example               % runnable end-to-end demo
```

`flatland_load` finds the library by searching, in order:

1. the path in the `FLATLAND_LIB` environment variable, if set;
2. the repository root next to `matlab/` (`../libflatland.so`), then `matlab/`
   itself, then `../build/`;
3. `libflatland` on the system loader path.

---

## Files

| File | Purpose |
|---|---|
| `Flatland.m` | the main `handle` class wrapping one mesh |
| `flatland_load.m` | loads / unloads / introspects the shared library |
| `flatland_test.m` | self-check, run inside MATLAB |
| `example.m` | runnable end-to-end example |
| `flatland_matlab.h` | the header `loadlibrary` parses (a shim, see below) |
| `stubinc/stddef.h`, `stubinc/stdint.h` | minimal system-header shadows, parse-time only |

Nothing outside `matlab/` is modified, and `include/flatland.h` is used exactly
as published.

### Why there is a shim header

`loadlibrary` is pointed at `flatland_matlab.h`, never at `include/flatland.h`.
The shim adds nothing to the ABI; it just presents the real header in a form
`loadlibrary`'s parser accepts, then includes it. Two problems are neutralised:

1. **`FL_API`.** Every declaration in `flatland.h` is prefixed with it, and it
   expands to `__attribute__((visibility("default")))` under GCC/Clang or
   `__declspec(dllimport)` on Windows. Confirmed with `gcc -E` against this
   repository, the parser would otherwise see:

   ```c
   __attribute__((visibility("default"))) fl_status fl_mesh_create(...)
   ```

   Attributes are not in the C subset `loadlibrary` documents. The shim
   defines `FLATLAND_STATIC` (which empties `FL_API` on the Windows arm) *and*
   `#undef`s `__GNUC__` (which empties it on the GCC arm), so all three
   preprocessor personalities land on an empty `FL_API`. This is a header-parse
   concern only: the symbols the `.so` exports are fixed at build time by the
   real `FL_API`, and all 31 are present.

2. **`<stddef.h>`.** GCC's real one defines `max_align_t` using `long double`,
   `__attribute__` and an unnamed struct — three separate constructs outside
   `loadlibrary`'s subset. `stubinc/` shadows `<stddef.h>` and `<stdint.h>`
   with minimal, ABI-identical typedefs for the parse only. Nothing compiled
   or linked uses them.

With both in place, the whole translation unit `loadlibrary` must digest is
140 non-blank lines containing no attributes, no `long double`, no unions, no
bitfields, no variadics and no function pointers.

---

## Indexing — read this one

Everything at the MATLAB boundary is **1-based**, and the binding converts to
the C API's 0-based indices internally. This is the single easiest thing to get
wrong when moving between the C API, the CLI and this binding.

| What | Here | In `flatland.h` |
|---|---|---|
| Face vertex indices | **1-based** | 0-based |
| `'FieldColumns'` | **1-based** | 0-based |
| Image rows | 1-based, **row 1 is the BOTTOM row** | row 0 is the bottom row |

```matlab
F = [1 4 3; 1 3 2; ...];    % correct: refers to vertices 1, 3 and 4
F = [0 3 2; 0 2 1; ...];    % WRONG: raises Flatland:InvalidFaces
```

`fl.faces()` returns 1-based indices too, so `Flatland(V, fl.faces())` round
trips. OBJ files on disk are 1-based, the loader converts to 0-based for C, and
`faces()` converts back — the round trip is consistent at every layer.

If your projected areas look plausible but your field values look scrambled,
check the indexing before anything else.

---

## NaN, not zero

**When a view covers no pixels, the field statistics are `NaN`.**

A view can legitimately cover nothing: geometry thinner than one pixel, or a
mesh seen exactly edge-on. The C struct signals this with `has_stats = 0` and
leaves `average`/`integral`/`min`/`max` at **zero**. Zero is a perfectly
ordinary field value, so passing it through would be indistinguishable from a
real measurement of zero — it would quietly drag any mean or sum you compute
toward zero.

So this binding reports those four as `NaN` whenever `hasStats` is false. A
view with no coverage has no statistics; it does not have statistics equal to
zero.

```matlab
r = sliver.project([0 0 -1], 'Field', f, 'Resolution', 0.1, 'Cull', false);
r.coveredPixels   % 0
r.hasField        % true  - a field WAS supplied
r.hasStats        % false - but nothing was measured
r.average         % NaN
r.area            % 0     - still a real measurement, so NOT NaN
```

`area`, `coveredPixels`, `width` and `height` are always real measurements and
are never replaced with `NaN`.

When aggregating a batch, guard on `hasStats`:

```matlab
R    = fl.projectBatch(views, 'FieldMatrix', M, 'FieldColumns', 1:K);
good = [R.hasStats];
m    = mean([R(good).average]);       % or use mean(..., 'omitnan')
```

---

## Row-major vs column-major

The C API reads its arrays **row-major**; MATLAB stores **column-major**. You
pass the MATLAB-natural shapes and the binding transposes:

| You pass | C receives |
|---|---|
| `vertices` N-by-3 | `vertex_count*3` doubles, xyz interleaved |
| `faces` M-by-3, 1-based | `face_count*3` int32, 0-based, interleaved |
| `views` K-by-3 | `view_count*3` doubles, interleaved |
| `'FieldMatrix'` nEntities-by-nTimesteps | row-major, `data[r*nCols + c]` |

`'FieldMatrix'` is the one that bites, because getting it wrong **does not
raise an error** — it silently returns wrong numbers. One row per mesh entity,
one column per timestep, which is what you would write anyway:

```matlab
M = repmat([10 20 30], 8, 1);    % 8 vertices x 3 timesteps
R = fl.projectBatch(views, 'FieldMatrix', M, 'FieldColumns', [1 2 3]);
[R.average]                       % 10  20  30
```

Internally the binding sends `M.'`, whose column-major flattening is exactly
C's row-major order. Sending `M` unchanged returns `20  18.33336  15` for that
same input — verified by running both linearisations through the C library.
`flatland_test` asserts against those wrong values as a regression guard.

---

## API

### Construction

```matlab
fl = Flatland(vertices, faces)   % N-by-3 double, M-by-3 1-based indices
fl = Flatland.load('part.obj')   % OBJ or STL, ASCII or binary, auto-detected
delete(fl)                       % frees the C handle; deterministic
```

`Flatland` is a `handle` subclass, so `delete` is deterministic and copies
share one mesh. Both input arrays are copied by the library; you may clear
yours immediately.

### `project`

```matlab
r = fl.project(view, 'Name', value, ...)
```

`view` is the direction the camera **looks along**; the visible surface is the
one facing back toward `-view`. Magnitude is irrelevant.

| Option | Default | Meaning |
|---|---|---|
| `'Field'` | none | one value per vertex (node) or per face (face) |
| `'Resolution'` | `1e-3` | pixel edge length in mesh units |
| `'Precision'` | `'single'` | `'single'`/`'float'` or `'double'` — engine working type only; everything crossing the boundary is double |
| `'Cull'` | `true` | backface cull |
| `'FieldMode'` | `'auto'` | `'node'` or `'face'`; needed only when vertex and face counts tie |

Returns a struct: `area`, `average`, `integral`, `min`, `max`,
`coveredPixels`, `width`, `height`, `hasField`, `hasStats`.

### `projectBatch`

```matlab
R = fl.projectBatch(views, 'Name', value, ...)
```

`views` is K-by-3. Returns a 1-by-K struct array; gather with `[R.area]`.
Adds, on top of `project`'s options minus `'Field'`:

| Option | Default | Meaning |
|---|---|---|
| `'FieldMatrix'` | none | nEntities-by-nTimesteps |
| `'FieldColumns'` | all column 1 | K **1-based** column indices, one per view |
| `'Resolutions'` | `'Resolution'` | K per-view pixel sizes |
| `'Threads'` | `0` | worker threads; 0 means one per core |

The whole view array and the whole field matrix cross the boundary in **one**
`calllib`. Looping per view from MATLAB would re-validate and re-centre the
mesh every time and throw away the parallelism this path exists for.

### `render`

```matlab
img = fl.render(view, 'Field', values, ...)
```

Returns `width`, `height`, `mask` (height-by-width logical), `values`
(height-by-width double, `NaN` off the silhouette, `[]` if no field) and
`result` (the same struct `project` returns).

**Row 1 is the bottom row** in mesh space, matching the C API:

```matlab
imagesc(flipud(img.values)); axis image off
```

A view covering nothing yields a valid 0-by-0 image, not a stale raster.

### Other

```matlab
V = fl.vertices()                    % N-by-3, original input frame
F = fl.faces()                       % M-by-3, 1-based
tf = fl.hasMesh()
fl.vertexCount, fl.faceCount

d = Flatland.angleToDir(az, el)      % degrees -> unit direction, N-by-3
v = Flatland.version()               % [major minor patch]
s = Flatland.versionString()
s = Flatland.statusString(code)

info = flatland_load()               % load (idempotent) and describe
flatland_load('unload')
flatland_load('reload')              % after rebuilding the library
tf = flatland_load('isloaded')
```

---

## Errors

Every `fl_status` is checked, and the text from `fl_last_error()` is folded
into the MATLAB message:

```
Error using Flatland/project
fl_project failed: dimension mismatch (fl_status 5): field has 7 values,
but the mesh has 8 vertices / 12 faces
```

Identifiers are stable and catchable:

| Identifier | Cause |
|---|---|
| `Flatland:InvalidField` | field length, shape, finiteness, or a bad `'FieldColumns'` |
| `Flatland:InvalidVertices` / `Flatland:InvalidFaces` | bad geometry input; the 0-based-faces message says so explicitly |
| `Flatland:InvalidView` | zero, non-finite, or wrongly shaped direction |
| `Flatland:InvalidOption` | unknown name/value option or out-of-range value |
| `Flatland:AmbiguousFieldMode` | vertex count equals face count; set `'FieldMode'` |
| `Flatland:NoMesh` | the object was deleted |
| `Flatland:FileNotFound` | `Flatland.load` on a missing path |
| `Flatland:NoCompiler` | `mex -setup C` has not been run |
| `Flatland:DimensionMismatch`, `Flatland:NumericError`, `Flatland:IOError`, `Flatland:ParseError`, `Flatland:InvalidArgument`, `Flatland:OutOfMemory`, `Flatland:Unsupported`, `Flatland:InternalError` | one per `fl_status` code, raised by the library |

---

## Measured struct layout

From `sizeof`/`offsetof` compiled against `include/flatland.h` on x86-64 Linux
(GCC 13). The binding relies on MATLAB marshalling these by field name, but
the `fl_result` size is also used by a byte-buffer fallback in `projectBatch`.

```
fl_options      56 bytes   resolution 0   cull 8   precision 12
                           field_mode 16  threads 20  reserved[8] 24

fl_result       96 bytes   area 0   average 8   integral 16  min 24  max 32
                           covered_pixels 40 (int64)
                           width 48  height 52  has_field 56  has_stats 60
                           reserved[8] 64

fl_batch_desc   88 bytes   views 0   view_count 8   field_matrix 16
                           field_rows 24  field_cols 32  field_columns 40
                           resolutions 48  reserved[8] 56

size_t 8   void* 8   int32_t 4   int64_t 8   double 8
enums (fl_status, fl_precision, fl_field_mode) 4 bytes each
```

---

## Troubleshooting

**`Flatland:NoCompiler`, or a parse error from `loadlibrary`**
Run `mex -setup C`. See the compiler prerequisite above.

**`Flatland:LibraryNotFound`, or `loadlibrary` cannot find the library**
Run `make lib` in the repository root, or set `FLATLAND_LIB` to the full path
of the built library.

**Linux: `GLIBCXX_3.4.x not found` when loading**
`libflatland.so` is C++ inside and links `libstdc++`. This build requires up to
`GLIBCXX_3.4.32`, and MATLAB ships its own older `libstdc++` in
`$MATLABROOT/sys/os/glnxa64` which it prefers. Start MATLAB with the system one
in front:

```bash
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 matlab
```

Or rebuild FlatLand against an older toolchain, or with `-static-libstdc++`.
Check what your MATLAB has with:

```bash
strings $MATLABROOT/sys/os/glnxa64/libstdc++.so.6 | grep GLIBCXX | sort -V | tail -3
```

**`Flatland:SymbolsMissing`**
The header declares functions the binary does not export — usually a stale
library beside a newer header. Rebuild with `make lib`.

**`Flatland:HandleNotReturned`**
A call succeeded but no handle came back. This is the one place the binding
depends on how your MATLAB release marshals a `T**` output argument; it tries
several forms. If you hit it, please report the output of
`libfunctions('flatland', '-full')`.

**Results changed after rebuilding the library**
MATLAB caches the loaded library for the session. Run `flatland_load('reload')`.

**Octave**
Octave does **not** provide `loadlibrary`/`calllib`, so this binding cannot run
there. Use the C API directly, the CLI, or Octave's own `mex`/`oct` interface.

---

## Verification status

The numeric expectations in `flatland_test.m` were produced by compiling and
running the equivalent C against the same `libflatland.so`, so they are
cross-checks against the C ABI rather than against this binding's own output:

| Check | Value |
|---|---|
| unit box, view `[1 0 0]`, resolution 0.002, double | `area = 1`, 250000 covered pixels, 502x502 raster |
| node field = vertex x, same view | `average = min = max = 0` |
| 3-view batch, 8x3 matrix of 10/20/30 | averages `10 / 20 / 30`, areas `1 / 1 / 1` |
| same matrix without the transpose | averages `20 / 18.33336 / 15` |
| sliver thinner than one pixel | `covered = 0`, `has_field = 1`, `has_stats = 0` |
| `cube.obj` | 26 vertices, 48 faces, area 1 |

The header's parseability, the exported symbol list (31 declared = 31
exported), and the struct layout above were all verified mechanically.

**`flatland_test.m` itself has not been executed** — no MATLAB installation was
available where this binding was written. Run it before relying on the binding;
it is written to fail loudly and specifically rather than quietly.
