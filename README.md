# FlatLand

FlatLand computes the **projected (visible) surface area** of a triangle mesh
from arbitrary viewing directions, by orthographic projection and software
rasterization.

Given a scalar field defined per vertex or per face it also reports, over the
visible projection: the mean, the extremes, and the **area integral**

```
I  =  ∫ f dA  =  Σ (value × pixel_area)
```

The integral is deliberately generic. Supply a radiance field and `I` is a
radiant intensity; supply a pressure field and it is a force; supply a
temperature field and it is a weighted thermal area. FlatLand attaches no
physical meaning to the numbers it is given — it computes a geometric quantity
and leaves the interpretation to you.

It is built for **time-dependent batches**: thousands of timesteps with a
changing view direction *and* a changing field over a fixed mesh, in seconds.

## Use it from

| | |
| :--- | :--- |
| **C / C++** | [`include/flatland.h`](include/flatland.h) — a flat C ABI, `libflatland.so`/`.dylib`/`.dll` |
| **Python** | [`python/`](python/README.md) — pure ctypes, numpy optional, no compiled extension |
| **MATLAB** | [`matlab/`](matlab/README.md) — `loadlibrary`/`calllib`, no MEX build |
| **Command line** | `./flatland mesh.obj -v 1 0 0` |
| **Anything else** | The C ABI is consumable from Julia, R, C#, Go, Rust — anything that speaks C |

Six worked examples, written twice so the two bindings are directly comparable:
projected area and orientation sweeps, geometry straight from memory, field
integrals, a parallel time series, getting at the raster, and checking the
answers against closed forms — see **[examples/README.md](examples/README.md)**.

The library takes **arrays**, not just files. You do not have to write a mesh to
disk and shell out to use it.

```python
import flatland
mesh = flatland.Mesh(vertices, faces)          # numpy arrays or plain lists
r = mesh.project((1, 0, 0), field=temperatures)
print(r.area, r.integral)
```

```c
fl_mesh* mesh = NULL;
fl_mesh_create(vertices, n_vertices, faces, n_faces, &mesh);

fl_options opts; fl_options_init(&opts);
opts.resolution = 1e-3;

fl_result r;
fl_project(mesh, (double[]){1, 0, 0}, field, n_vertices, &opts, &r);
printf("area %g   integral %g\n", r.area, r.integral);
fl_mesh_destroy(mesh);
```

```matlab
fl = Flatland(vertices, faces);                % 1-based faces
r  = fl.project([1 0 0], 'Field', values);
```

## Design goals

- **Correct, and shown to be.** Checked against closed-form results — convex
  projection formulas, Cauchy's projection identity, analytic field integrals,
  blackbody radiant intensities and the Stefan–Boltzmann law — not against its
  own previous output. See [validation](#validation).
- **Fast.** Multi-threaded batches: ~1000 changing-field, changing-view timesteps
  on a 70k-triangle mesh in a couple of seconds.
- **Dependency-free.** Standard C++17 and nothing else. One compiler invocation
  builds it, which matters on minimal or airgapped systems.

## Install

### Prebuilt binaries

Each release attaches archives for linux-x86_64, linux-aarch64, macos-arm64,
macos-x86_64 and windows-x86_64, containing the CLI, the shared and static
libraries, and the header. Download from
[Releases](https://github.com/Interfluo/flatland/releases) and verify against the
published `SHA256SUMS`.

### From source

```shell
make                 # the CLI               -> ./flatland
make lib             # the libraries         -> libflatland.{a,so}
make all             # both
make test            # build and run the test suite
make validate        # the full validation study
make install         # to $PREFIX, default /usr/local
```

Or with CMake:

```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

The only requirement is a C++17 compiler with threads, which stock GCC, Clang
and MSVC all provide.

> Linking the **static** library on Windows requires defining `FLATLAND_STATIC`
> so the header declares imports correctly. CMake consumers get this
> automatically from the target.

## Concepts

**View direction.** The direction the camera *looks along*, so the visible
surface of a closed mesh is the one whose outward normals oppose it. Magnitude is
irrelevant. Specify it as a vector (`-v x y z`) or by azimuth/elevation in
degrees (`-a az el`), where azimuth sweeps around `+Z` from `+X` and elevation
rises from the XY plane.

**Resolution.** The pixel edge length. Coverage is decided by testing whether a
pixel *centre* falls inside a triangle, so area converges as resolution
increases. Choose one fine enough to resolve thin features, and sweep it to check
the answer has settled.

**Fields.** One value per vertex (interpolated across each triangle) or per face
(constant over it). A field file is a matrix: one row per mesh entity, one column
per timestep, so an entire time series lives in one file. Select a column with
`@<col>`, e.g. `fields.txt@7`. STL carries no shared vertices, so STL fields are
face-based.

**`has_stats`.** A view can legitimately cover no pixels — geometry thinner than
one pixel, or a mesh seen exactly edge-on. The four field statistics are then
*absent*, not zero: the C API sets `has_stats = 0`, Python returns `None`, MATLAB
returns `NaN`, and the CLI's JSON emits `null`. Check it before consuming a
batch; a fabricated zero is indistinguishable from a measurement once it reaches
your min/max.

**Rasters.** Two output kinds, because they answer different questions. `-o`
writes a **PNG** — the field mapped through a colour ramp to 8 bits per channel,
which is right for looking at and wrong for computing with. `--npy` writes the
**numbers**: NumPy `.npy`, float64, shape `(height, width)`, `NaN` wherever
nothing was covered. Row 0 is the bottom row in mesh space, so plot it with
`origin='lower'`.

The PNG encoder is built in — FlatLand links no image or compression library,
so the dependency-free property survives. On a typical heatmap that is about
10x smaller than the raw raster; on a silhouette, 100x.

## Command line

```shell
./flatland <mesh.obj|.stl> [options]
```

| Flag | Arguments | Description |
| :--- | :-------- | :---------- |
| `-v, --view` | x y z | Add a view direction. Repeatable. |
| `-a, --angle` | az el | Add a view by azimuth/elevation in degrees. Repeatable. |
| `-b, --batch` | file | Load views/timesteps from a file. Additive with `-v`/`-a`. |
| `-r, --res` | val | Default pixel resolution (default `0.001`). |
| `-d, --data` | file[@col] | Default scalar field file. |
| `--field-mode` | node\|face\|auto | Force field interpretation (default `auto`, from row count). |
| `-p, --precision` | float\|double | Working precision (default `float`). |
| `-t, --threads` | n | Worker threads for batch views (0 = one per core). |
| `--no-cull` | — | Disable backface culling (render all faces). |
| `-o, --out` | prefix | Save PNG heatmaps as `<prefix>_<idx>.png`. |
| `--npy` | prefix | Save raw per-pixel field values as `<prefix>_<idx>.npy`. |
| `-j, --json` | — | Emit structured JSON to stdout. |
| `-h, --help` | — | Show help. |

### Per-view outputs

| Field | Meaning |
| :--- | :--- |
| `area` | Visible projected area = `covered_pixels × resolution²` |
| `average` | Mean field value over visible pixels |
| `integral` | Area integral `Σ value × pixel_area` |
| `min` / `max` | Field extremes over visible pixels |

### Batch files — the time-series interface

One line per timestep. Each line is a direction vector or an azimuth/elevation
pair, followed by an optional resolution and/or field selector:

```text
<nx> <ny> <nz>  [resolution] [data[@col]]      # direction vector
a <az> <el>     [resolution] [data[@col]]      # azimuth / elevation (degrees)
```

- `resolution` and `data` are optional and order-independent. A token naming an
  existing file is data; otherwise a number is a resolution. The filesystem gets
  the deciding vote so that a data file named for its timestep (`0100`) is not
  mistaken for a resolution.
- Pair one field matrix with `@<col>` to advance the field per timestep without
  one file per step.
- Relative data paths resolve against the **batch file's** directory, so a case
  folder is portable.
- A token beginning with `#` starts a comment; a `#` inside a path does not.

```text
# a <azimuth> <elevation> <resolution> field.txt@<timestep>
a   0 20 0.002 field.txt@0
a  15 20 0.002 field.txt@1
a  30 20 0.002 field.txt@2
```

```shell
./flatland part.obj -b timeseries.txt -j > results.json
```

Each distinct field file is parsed once into a shared, bounded cache; workers
extract their timestep's column concurrently.

## Examples

```shell
./flatland bunny.obj -v 1 0 0                        # projected area down +X
./flatland part.stl -a 45 30                         # STL from azimuth 45, elevation 30
./flatland engine.obj -v 0 1 1 -d temps.txt -o heat  # field heatmap to heat_0000.png
./flatland part.obj -b timeseries.txt -j > out.json  # time series to JSON
```

`examples/with_fields/generate_field.py` builds an orbiting-view, evolving-field
case using only the Python standard library. A ready-to-run demo is committed at
`examples/with_fields/timeseries_demo/`.

For the library rather than the CLI, [examples/README.md](examples/README.md)
indexes five worked examples in Python and MATLAB, each answering a question
rather than touring the API.

## Validation

FlatLand is checked against closed-form results, not against its own previous
output — **[docs/VALIDATION.md](docs/VALIDATION.md)** has the derivations and the
measured tables. Highlights:

| Study | Result |
| :--- | :--- |
| Exact convex polyhedra (cube, tetrahedron, octahedron), 21 shape/direction pairs | worst relative error 3.2×10⁻⁴ |
| Cauchy's identity `⟨A⟩ = S/4`, four shapes, 400 directions | worst 1.2×10⁻⁵ |
| Sphere vs. its own mesh's exact projected area, 5 subdivision levels | ~10⁻⁶, flat in subdivision |
| Cylinder `2rh sin θ + πr² cos θ`, seven angles | worst 4.6×10⁻⁵ |
| Lambertian sphere `∫cos dA = (2/3)πr²` | 5.4×10⁻⁴ at 20480 triangles |
| Blackbody sphere `I = σT⁴r²`, twelve directions | 2.8–3.2×10⁻⁴; the mean radiance to 1.4×10⁻¹¹ |
| Stefan–Boltzmann `4π⟨I⟩ = σT⁴S`, three convex shapes | worst 1.1×10⁻⁵ |
| Radiative-equilibrium sphere vs. the Lambert phase function | 10⁻⁷–10⁻⁴ against its own mesh |
| Translation by 5×10⁶, and culled vs. `--no-cull` | exact |

The sphere study is the informative one: rasterization error stays near 10⁻⁶ and
does **not** improve with subdivision, while the gap to `πr²` falls fourfold per
level. FlatLand reproduces whatever mesh it is handed to about six digits — the
remaining difference from a true sphere belongs to the mesh, and closing it needs
more triangles, not smaller pixels.

Run it yourself with `make validate`.

## Testing

```shell
make test                              # everything, 548 assertions
./tests/run_tests.sh ./flatland cli    # one suite
```

Suites: `geometry`, `parsing`, `cli`, `batch`, `capi`, `python`, `validation`.
CI builds with GCC and Clang on Linux and macOS, MSVC and MinGW on Windows, via
both build systems, and runs AddressSanitizer, UndefinedBehaviorSanitizer and
ThreadSanitizer.

## Tips

1. **Batch over repeated invocations.** One `-b` run parses the mesh once;
   invoking the binary repeatedly re-parses it every time. From the library, build
   the mesh once and call `fl_project_batch`.
2. **Pick resolution deliberately.** Sweep `-r` with `-o` when setting up a new
   case, so thin geometry is not under-resolved. Convergence is first order and
   not monotonic — see the convergence table in the validation doc.
3. **Threads.** Batch views run in parallel automatically; cap with `-t` to leave
   cores free.
4. **Culling.** `--no-cull` is for non-manifold or open surfaces where "inside"
   faces should still contribute to the silhouette. On a closed mesh it changes
   nothing but the speed.

## Licensing

FlatLand is dual-licensed.

- **[GNU AGPL v3](LICENSE)** — free for any use, including commercial use,
  provided anything you build on it is also released under the AGPL. Using
  FlatLand internally, on your own geometry, in your own pipeline, is free and
  carries no obligations: the AGPL attaches to distribution and to
  network-accessible services, not to in-house use.
- **[Commercial license](COMMERCIAL-LICENSE.md)** — required to ship FlatLand
  inside a closed-source product, or to offer it as a hosted service, without
  releasing your own source. Terms are negotiated per engagement.

Releases up to and including `v2` were published under the GNU GPL v3 and remain
available under it; the AGPL applies from the relicensing commit onward.
