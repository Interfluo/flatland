# FlatLand

FlatLand is a small, dependency-free C++ utility for the geometric analysis of 3D
meshes (OBJ) via orthographic projection and software rasterization.

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

### Options
| Flag | Arguments | Description |
| :--- | :-------- | :---------- |
| `-v, --view` | x y z | Add a single view direction. Repeatable. |
| `-b, --batch` | file | Load views/timesteps from a file (see format below). Additive with `-v`. |
| `-r, --res` | val | Default pixel resolution (default `0.001`). |
| `-d, --data` | file | Default scalar field file (one value per line; length = vertex or face count). |
| `-p, --precision` | mode | `float` (default, faster) or `double` (more accurate). |
| `-t, --threads` | n | Worker threads for batch views (default: auto / all cores). |
| `--no-cull` | — | Disable backface culling (render all faces). |
| `-o, --out` | prefix | Save PPM heatmaps as `<prefix>_<idx>.ppm`. |
| `-j, --json` | — | Emit structured JSON to stdout. |
| `-h, --help` | — | Show help. |

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

The batch file is the time-series interface: **one line per timestep**, each with its
own view direction and (optionally) its own resolution and field file.

```text
<nx> <ny> <nz> [resolution] [data_file]
```
- `resolution` and `data_file` are optional; omitted values fall back to `-r` / `-d`.
- Relative `data_file` paths are resolved against the **batch file's directory**, so a
  case folder is portable and runs from any working directory.
- Lines starting with `#` are comments.

A fixed field (one `-d` file) is loaded once and shared across all timesteps; per-line
field files are streamed on demand, so memory stays bounded by the thread count rather
than the number of timesteps.

**Example** `timeseries.txt`:
```text
# nx ny nz resolution field_file
1.0 0.0 0.3 0.002 field_0000.txt
0.0 1.0 0.3 0.002 field_0001.txt
-1.0 0.0 0.3 0.002 field_0002.txt
```
```shell
./flatland part.obj -b timeseries.txt -j > results.json
```

### Generating a time series
`examples/with_fields/generate_field.py` builds a rotating-view / evolving-field case
using **only the Python standard library** (no numpy/pyvista):
```shell
python3 examples/with_fields/generate_field.py bunny.obj data 1000 0.002
./flatland bunny.obj -b data/timeseries.txt -j > results.json
```
A small ready-to-run demo is committed at `examples/with_fields/timeseries_demo/`:
```shell
./flatland examples/sphere_areas/sphere_coarse.obj \
  -b examples/with_fields/timeseries_demo/timeseries.txt
```

## Examples
**1. Basic area** — projected area of the bunny down the X axis:
```shell
./flatland bunny.obj -v 1 0 0
```
**2. Field heatmap, double precision** — map `temps.txt`, render to `heat_0000.ppm`:
```shell
./flatland engine.obj -v 0 1 1 -d temps.txt -p double -o heat
```
**3. JSON pipeline** — batch of views to JSON:
```shell
./flatland part.obj -b views.txt -j > results.json
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
