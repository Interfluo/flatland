# Examples

Five worked examples, written twice — once in Python, once in MATLAB — so the
two are directly comparable. Each answers a question rather than touring an API.

| | Question | Python | MATLAB |
| :-- | :--- | :--- | :--- |
| 1 | How much of my part does a direction see, and which orientation exposes least? | [`ex01_projected_area.py`](python/ex01_projected_area.py) | [`ex01_projected_area.m`](matlab/ex01_projected_area.m) |
| 2 | My mesh is already in memory — must I write it to disk? | [`ex02_mesh_from_arrays.py`](python/ex02_mesh_from_arrays.py) | [`ex02_mesh_from_arrays.m`](matlab/ex02_mesh_from_arrays.m) |
| 3 | What is the integral of my scalar field over the visible part? | [`ex03_field_integral.py`](python/ex03_field_integral.py) | [`ex03_field_integral.m`](matlab/ex03_field_integral.m) |
| 4 | I have hundreds of timesteps. How do I not wait all afternoon? | [`ex04_timeseries.py`](python/ex04_timeseries.py) | [`ex04_timeseries.m`](matlab/ex04_timeseries.m) |
| 5 | I want to see the raster, not just the number. | [`ex05_heatmap.py`](python/ex05_heatmap.py) | [`ex05_heatmap.m`](matlab/ex05_heatmap.m) |

There is also an API tour in each binding — [`python/example.py`](../python/example.py)
and [`matlab/example.m`](../matlab/example.m) — which walk the whole surface
rather than a use case.

## Running them

Build the shared library first; everything here loads it at runtime.

```shell
cd /path/to/flatland && make lib
```

**Python** — no installation needed, the examples find the package in this
checkout:

```shell
python3 examples/python/ex01_projected_area.py
python3 examples/python/ex05_heatmap.py /tmp/out      # writes images there
```

NumPy is optional. Everything runs without it; example 4 explains the one place
it matters for speed.

**MATLAB** — needs a C compiler configured (`mex -setup C`), because
`loadlibrary` preprocesses the header with one. See
[`matlab/README.md`](../matlab/README.md) for the setup and the
troubleshooting list.

```matlab
>> cd /path/to/flatland/examples/matlab
>> ex01_projected_area
```

## Things worth knowing before you start

**A view direction is the direction the camera looks along.** The visible
surface is the one facing back toward `-direction`. Magnitude is irrelevant.

**MATLAB is 1-based, C and Python are 0-based.** Face indices and
`FieldColumns` are 1-based in the MATLAB binding and converted internally. Do
not subtract one yourself; the binding checks and tells you so.

**Winding decides which surface you measure.** Triangles must run
counter-clockwise seen from outside. A mesh wound inside-out reports the *far*
surface — the area looks right and every field statistic is wrong. Example 2
demonstrates that failure deliberately.

**When a view covers no pixels there are no statistics.** Not zero — absent.
Python returns `None`, MATLAB returns `NaN`, the C API sets `has_stats = 0`, and
the CLI's JSON emits `null`. Zero is an ordinary field value, so passing it
through would be indistinguishable from a measurement. Guard on it before
aggregating a batch.

**Row 0 (Python) / row 1 (MATLAB) of a raster is the BOTTOM row** in mesh
space, which is the opposite of most image formats. Flip before writing a PNG,
or use `origin='lower'` / `axis xy` when displaying.

## Reference data

The meshes and field files the examples and tests use:

| Path | What |
| :--- | :--- |
| `cube_area/cube.obj` | unit cube, 26 vertices, 48 triangles — has a closed-form projected area |
| `sphere_areas/sphere_coarse.obj` | icosphere, 320 triangles |
| `sphere_areas/sphere_fine.obj` | icosphere, 20480 triangles |
| `with_fields/bunny.obj` | Stanford bunny, 69451 triangles — the performance case |
| `with_fields/generate_field.py` | builds an orbiting-view, evolving-field case for the CLI |
| `with_fields/timeseries_demo/` | a committed ready-to-run time series |

For the analytic machinery behind the closed-form checks these examples make,
see [`docs/VALIDATION.md`](../docs/VALIDATION.md) and `validation/cases.py`.
