# FlatLand

FlatLand is a small, dependency-free C++ utility for the geometric analysis of 3D
meshes (OBJ and STL) via orthographic projection and software rasterization.

Its core job is computing the **projected (visible) surface area** of a mesh from
arbitrary viewing angles. On top of that it can map a **scalar field** (defined per
node or per face — temperature, pressure, radiance, stress, …) onto the visible
projection and report per-view statistics: average, min, max, and the **area
integral** `∫ f dA`. The integral is deliberately generic: supply a radiance field
and it is radiant intensity; supply pressure and it is a force; FlatLand itself
stays domain-agnostic.

It is built for **time-dependent batches** — sweep thousands of timesteps with a
changing view direction *and* a changing field over a fixed mesh, in seconds.

## Design goals
- **Reliable** — clean error handling (no crashes on bad input), a test suite, and CI.
- **Fast** — multi-threaded batch processing; ~1000 changing-field/changing-view
  timesteps on a 70k-triangle mesh in a couple of seconds.
- **Dependency-free** — standard C++17 and nothing else. Builds with a single
  compiler invocation, which matters on minimal or airgapped systems.

## Build

FlatLand is a single `.cpp` file. No third-party libraries.

**Simplest recipe (airgapped-friendly):**
```shell
make                 # produces ./flatland
```

**Or a single compiler call, no build system at all:**
```shell
c++ -std=c++17 -O2 -pthread src/main.cpp -o flatland
```

**Or via CMake (3.14+):**
```shell
cmake -S src -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Verify:
```shell
./flatland --help
make test            # build + run the test suite
```

> The only system requirement is a C++17 compiler with threads (`-pthread`),
> available out of the box with stock GCC or Clang.

## Usage
```shell
./flatland <mesh.obj> [options]
```

### Meshes
OBJ and STL (binary or ASCII) are accepted, auto-detected by file extension. STL
triangles carry no shared vertices, so STL fields are face-based.

### Options
| Flag | Arguments | Description |
| :--- | :-------- | :---------- |
| `-v, --view` | x y z | Add a view direction (vector). Repeatable. |
| `-a, --angle` | az el | Add a view by azimuth/elevation in degrees. Repeatable. |
| `-b, --batch` | file | Load views/timesteps from a file (see format below). Additive with `-v`/`-a`. |
| `-r, --res` | val | Default pixel resolution (default `0.001`). |
| `-d, --data` | file[@col] | Default scalar field file (see [Field files](#field-files)). |
| `--field-mode` | node\|face\|auto | Force field interpretation (default `auto`, inferred from row count). |
| `-p, --precision` | mode | `float` (default, faster) or `double` (more accurate). |
| `-t, --threads` | n | Worker threads for batch views (default: auto / all cores). |
| `--no-cull` | — | Disable backface culling (render all faces). |
| `-o, --out` | prefix | Save PPM heatmaps as `<prefix>_<idx>.ppm`. |
| `-j, --json` | — | Emit structured JSON to stdout. |
| `-h, --help` | — | Show help. |

### Views
A view is a direction the camera looks along. Specify it as a raw vector (`-v x y z`)
or by **azimuth/elevation in degrees** (`-a az el`), where azimuth sweeps around `+Z`
from `+X` and elevation rises from the XY-plane. Both forms are repeatable and can be
mixed with a batch file.

### Field files
A field file is a matrix: **one row per mesh entity** (vertex or face), and **one or
more whitespace-separated columns**. Each column is a timestep, so an entire time
series can live in a single file. Append `@<col>` to a data path to pick a column
(default `0`), e.g. `fields.txt@7`. A classic one-value-per-line file is just the
single-column case. Node-vs-face is inferred from the row count; use `--field-mode`
to force it when the vertex and face counts coincide.

### Per-view outputs
| Field | Meaning |
| :--- | :--- |
| `area` | Visible projected area = `covered_pixels × resolution²`. |
| `average` | Mean field value over visible pixels (requires a data file). |
| `integral` | Area integral of the field: `Σ value × pixel_area` (= radiant intensity, force, …). |
| `min` / `max` | Field extremes over visible pixels. |

> **Accuracy note:** coverage is computed by point-sampling pixel centers with edge
> functions and a Z-buffer. Area therefore converges as resolution increases; choose
> a resolution fine enough to resolve thin features (use `-o` to inspect).

## Time-dependent / batch workflow

The batch file is the time-series interface: **one line per timestep**. Each line is
either a direction vector or an azimuth/elevation angle, followed by an optional
resolution and/or field selector:

```text
<nx> <ny> <nz>  [resolution] [data[@col]]      # direction vector
a <az> <el>     [resolution] [data[@col]]      # azimuth / elevation (degrees)
```
- `resolution` and `data` are optional and order-independent — a pure number is read
  as a resolution, anything else as a data-file path. Omitted values fall back to
  `-r` / `-d`.
- Pair a single field matrix with `@<col>` to advance the field per timestep
  (`field.txt@0`, `field.txt@1`, …) without one file per step.
- Relative data paths resolve against the **batch file's directory**, so a case folder
  is portable and runs from any working directory.
- Lines starting with `#` are comments.

Each distinct field file is loaded once into a shared, read-only cache; workers extract
their timestep's column concurrently, so a single shared matrix is never re-parsed
per thread.

**Example** `timeseries.txt` (orbit in azimuth, advancing one field column per step):
```text
# a <azimuth> <elevation> <resolution> field.txt@<timestep>
a   0 20 0.002 field.txt@0
a  15 20 0.002 field.txt@1
a  30 20 0.002 field.txt@2
```
```shell
./flatland part.obj -b timeseries.txt -j > results.json
```

### Generating a time series
`examples/with_fields/generate_field.py` builds an orbiting-view / evolving-field case
using **only the Python standard library** (no numpy/pyvista). It writes a single
multi-column field file plus an angle-based batch:
```shell
python3 examples/with_fields/generate_field.py bunny.obj data 1000 0.002
./flatland bunny.obj -b data/timeseries.txt -j > results.json
```
A small ready-to-run demo (two files: one matrix + one batch) is committed at
`examples/with_fields/timeseries_demo/`:
```shell
./flatland examples/sphere_areas/sphere_coarse.obj \
  -b examples/with_fields/timeseries_demo/timeseries.txt
```

## Examples
**1. Basic area** — projected area of the bunny down the X axis:
```shell
./flatland bunny.obj -v 1 0 0
```
**2. STL + angle view** — an STL part viewed from azimuth 45°, elevation 30°:
```shell
./flatland part.stl -a 45 30
```
**3. Field heatmap, double precision** — map `temps.txt`, render to `heat_0000.ppm`:
```shell
./flatland engine.obj -v 0 1 1 -d temps.txt -p double -o heat
```
**4. Time series from one matrix** — orbit views, one field column per step, to JSON:
```shell
./flatland part.obj -b timeseries.txt -j > results.json
```

## Validation

Two analytical checks (also enforced by the test suite):
- **Cube** — the projected area of a unit cube along an axis is exactly `1`.
- **Sphere** — the projected area of a unit sphere is `π` from any direction; FlatLand
  converges to `π` as mesh and pixel resolution increase.

Run the suite (dependency-free, bash + awk):
```shell
make test            # or: ./tests/run_tests.sh
```
CI builds with GCC and Clang (Linux + macOS) and via CMake on every push — see
`.github/workflows/ci.yml`.

## Tips
1. **Batch over repeated CLI calls** — one `-b` run loads and parses the mesh once;
   repeatedly invoking the binary re-parses the OBJ every time.
2. **Pick resolution deliberately** — sweep `-r` with `-o` image output when setting up
   a new case so thin geometry isn't under-resolved.
3. **Threads** — batch views run in parallel automatically; cap with `-t` if you need
   to leave cores free.
4. **Culling** — `--no-cull` is for non-manifold or open surfaces where "inside" faces
   should still contribute to the silhouette.
